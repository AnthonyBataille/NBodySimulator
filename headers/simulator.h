#ifndef SIMULATOR_H
#define SIMULATOR_H
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <random>
#include <chrono>
#include <ctime>
#include <omp.h>


using namespace boost::numeric::ublas;

template <typename T> inline void vector_delete_elem(vector<T>& v, const std::size_t i){ // remove element i of vector
	if(i < v.size()){
		v(i) = v(v.size() - 1);
		v.resize(v.size() - 1);
	}
}

// Taken from http://www.crystalclearsoftware.com/cgi-bin/boost_wiki/wiki.pl?Examples_-_How_To_Extend_UBLAS
// (op v) [i] = op( v [i] )
template<typename OP, typename E> 
BOOST_UBLAS_INLINE
typename vector_unary_traits<E, OP>::result_type
apply_to_all(const vector_expression<E> &e, const OP& op = OP()){
	typedef typename vector_unary_traits<E, OP>::expression_type expression_type;
	return expression_type (e());
}

namespace functor{
	template <typename T> class sqrt{ // apply square root to all elements
	public:
		typedef T value_type;
		typedef T result_type;
		sqrt() {}
    
		static BOOST_UBLAS_INLINE result_type apply(const value_type& x){
		return std::sqrt(x);
		}
	};
	template <typename T> class inv{ // inverse all elements
	public:
		typedef T value_type;
		typedef T result_type;
		inv() {}
    
		static BOOST_UBLAS_INLINE result_type apply(const value_type& x){
		return 1.0 / x;
		}
	};
}

template <typename real_type> class integrator;

template <typename real_type> class simulator{
	const int nthreads = 4; // number of threads for OpenMP
	
	// Reference units in S.I
	const double M_ref = 1.99e30; // Mass of the sun
	const double L_ref = 1.5e11; // 1 astronomical unit
	const double T_ref = 3600.0*24.0*365.0*10.0; // 10 years

	public:
		typedef vector<real_type> vec;
		simulator() = default;
		simulator(int N, real_type flength):N(N), flength(flength){
			G = real_type(6.674e-11 * M_ref * T_ref*T_ref / pow(L_ref, 3.0));
			
			std::mt19937 rand_generator(std::random_device{}());
			std::normal_distribution<real_type> normal_dist(0.0, flength * 0.1);
			std::uniform_real_distribution<real_type> uni_dist(-flength / 2.0, flength / 2.0);
	
			vel_x = vec(N);
			vel_y = vec(N);
			for(int i = 0; i < N; ++i){
				vel_x(i) = normal_dist(rand_generator); // generate velocities according to a normal distribution
				vel_y(i) = normal_dist(rand_generator);
			}
			point_x = vec(N);
			point_y = vec(N);
			for(int i = 0; i < N; ++i){
				point_x(i) = uni_dist(rand_generator); // generate points according to a uniform distribution in the field
				point_y(i) = uni_dist(rand_generator);
			}
			
			mass = scalar_vector<real_type>(N, 1.0);
			accel_x = vec(N, 0.0);
			accel_y = vec(N, 0.0);
		}
		
		void init_solar_system(real_type M){
			std::mt19937 rand_generator(std::random_device{}());
			std::uniform_real_distribution<real_type> uni_dist_radius(std::sqrt(M), flength / 2.0);
			std::uniform_real_distribution<real_type> uni_dist_angle(0.0, 2 * std::acos(-1));
			point_x(0) = 0.0;
			point_y(0) = 0.0;
			
			vel_x(0) = 0.0;
			vel_y(0) = 0.0;
			
			mass(0) = M;
			
			for(int i = 1; i < N; ++i){
				real_type radius = uni_dist_radius(rand_generator);
				real_type angle = uni_dist_angle(rand_generator);
				point_x(i) = radius * std::cos(angle);
				point_y(i) = radius * std::sin(angle);
				
				real_type r = std::sqrt(point_x(i)*point_x(i) + point_y(i)*point_y(i));

				vel_x(i) = - std::sqrt(G*M/r) * point_y(i)/r;
				vel_y(i) = std::sqrt(G*M/r) * point_x(i)/r;
			}
		}
		
		void merge_point(const real_type K, bool profile){ // Merge two points if they collide. K is the "real length" per pixel constant (squared)
			double begin = 0.0;
			if(profile){
				begin = omp_get_wtime();
			}
			for(int i = 0; i < N; ++i){
				for(int j = 0; j < i ; ++j){
					real_type dist_x = point_x(j) - point_x(i);
					real_type dist_y = point_y(j) - point_y(i);
					real_type r2 = dist_x*dist_x + dist_y*dist_y;
					if(r2 < K * 2.0 * std::sqrt(mass(i)*mass(j))+mass(i)+mass(j)){
						int& k = mass(i) < mass(j) ? i : j;
						int& l = mass(i) < mass(j) ? j : i;
						real_type sum_mass = mass(k) + mass(l);
						vel_x(l) = (mass(k)*vel_x(k) + mass(l)*vel_x(l)) / sum_mass;
						vel_y(l) = (mass(k)*vel_y(k) + mass(l)*vel_y(l)) / sum_mass;

						mass(l) += mass(k);
						vector_delete_elem<real_type>(mass, k);
						vector_delete_elem<real_type>(accel_x, k);
						vector_delete_elem<real_type>(accel_y, k);
						vector_delete_elem<real_type>(vel_x, k);
						vector_delete_elem<real_type>(vel_y, k);
						vector_delete_elem<real_type>(point_x, k);
						vector_delete_elem<real_type>(point_y, k);
						--N;
						--k;
						--l;
					}
				}
			}
			if(profile){
				double wtime = (omp_get_wtime() - begin) * 1000.0;
				std::cout << "merge_point : " << wtime << " ms" << std::endl;
			}
		}
		
		void update_accel(bool profile){
			double begin = 0.0;
			if(profile){
				begin = omp_get_wtime();
			}
			#pragma omp parallel for schedule(static) num_threads(nthreads)
			for(int i = 0; i < N; ++i){ // Compute acceleration of every particle i based on Newton law of gravitation
				vec p_x{point_x};
				vec p_y{point_y};

				vector_delete_elem<real_type>(p_x, i);
				vector_delete_elem<real_type>(p_y, i);

				vec r_x{p_x - scalar_vector<real_type>(N - 1, point_x(i))};
				vec r_y{p_y - scalar_vector<real_type>(N - 1, point_y(i))};

				vec r2_inv{element_prod(r_x, r_x) + element_prod(r_y, r_y)};
				r2_inv = apply_to_all(r2_inv, functor::inv<real_type>());
				vec r_inv{apply_to_all(r2_inv, functor::sqrt<real_type>())};

				vec r32_inv{element_prod(r2_inv, r_inv)};
				
				r_x = element_prod(r_x, r32_inv);
				r_y = element_prod(r_y, r32_inv);
				
				vec m{mass};
				vector_delete_elem<real_type>(m, i);

				accel_x(i) = G * inner_prod(m, r_x);
				accel_y(i) = G * inner_prod(m, r_y);
			}
			if(profile){
				double wtime = (omp_get_wtime() - begin) * 1000.0;
				std::cout << "update_accel : " << wtime << " ms" << std::endl;
			}
		}
		
		void update_point(const integrator<real_type>& in, bool profile){
			double begin = 0.0;
			if(profile){
				begin = omp_get_wtime();
			}
			
			in.integrate_step(*this);
			if(profile){
				double wtime = (omp_get_wtime() - begin) * 1000.0;
				std::cout << "update_point : " << wtime << " ms" << std::endl;
			}
		}
		
		real_type compute_energy(bool profile) const{ // Computes the total energy (Hamiltonian) of the system
			double begin = 0.0;
			if(profile){
				begin = omp_get_wtime();
			}
			
			vec T_vec{0.5 * element_prod(mass, element_prod(vel_x, vel_x) + element_prod(vel_y, vel_y))};
			real_type T = sum(T_vec);
			
			real_type V = 0.0;
			
			for(int i = 0; i < N; ++i){
				vec p_x{point_x};
				vec p_y{point_y};

				vector_delete_elem<real_type>(p_x, i);
				vector_delete_elem<real_type>(p_y, i);

				vec r_x{p_x - scalar_vector<real_type>(N - 1, point_x(i))};
				vec r_y{p_y - scalar_vector<real_type>(N - 1, point_y(i))};

				vec r2_inv{element_prod(r_x, r_x) + element_prod(r_y, r_y)};
				r2_inv = apply_to_all(r2_inv, functor::inv<real_type>());
				vec r_inv{apply_to_all(r2_inv, functor::sqrt<real_type>())};
				
				vec m{mass};
				vector_delete_elem<real_type>(m, i);
				
				 V -= G * mass(i) * sum(element_prod(m, r_inv));
			}
			
			if(profile){
				double wtime = (omp_get_wtime() - begin) * 1000.0;
				std::cout << "compute_energy : " << wtime << " ms" << std::endl;
			}
			
			return T + V;
		}
		
		int N;
		real_type flength;
		real_type G;
		
		vec point_x;
		vec point_y;
		vec vel_x;
		vec vel_y;
		vec accel_x;
		vec accel_y;
		vec mass;
};

template <typename real_type> class integrator{
	public:
		integrator() = default;
		integrator(const real_type timestep, const bool profile): dt(timestep), profile_flag(profile){}
		virtual void integrate_step(simulator<real_type>& s) const{}
	protected:
		real_type dt;
		bool profile_flag;
};

template <typename real_type> class euler_integrator : public integrator<real_type>{
	using integrator<real_type>::profile_flag;
	using integrator<real_type>::dt;
	public:
		euler_integrator() = default;
		euler_integrator(const real_type timestep, const bool profile): integrator<real_type>(timestep, profile){}
		void integrate_step(simulator<real_type>& s) const{
			s.update_accel(profile_flag);
			
			s.point_x += dt*s.vel_x;
			s.point_y += dt*s.vel_y;
		
			s.vel_x += dt*s.accel_x;
			s.vel_y += dt*s.accel_y;
		}
};

template <typename real_type> class leapfrog_integrator : public integrator<real_type>{
	using integrator<real_type>::profile_flag;
	using integrator<real_type>::dt;
	public:
		leapfrog_integrator() = default;
		leapfrog_integrator(const real_type timestep, const bool profile, simulator<real_type>& s): integrator<real_type>(timestep, profile){
			s.vel_x += dt/2.0 * s.accel_x;
			s.vel_y += dt/2.0 * s.accel_y;
		}
		void integrate_step(simulator<real_type>& s) const{
			s.point_x += dt * s.vel_x;
			s.point_y += dt * s.vel_y;
			
			s.update_accel(profile_flag);
		
			s.vel_x += dt * s.accel_x;
			s.vel_y += dt * s.accel_y;
		}
};

#endif
