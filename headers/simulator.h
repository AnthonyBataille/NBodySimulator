#ifndef SIMULATOR_H
#define SIMULATOR_H
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <random>
#include <chrono>
#include <ctime>
#include <omp.h>


using namespace boost::numeric::ublas;

const int nthreads = 4; // number of threads for OpenMP

// Reference units in S.I
const double M_ref = 1.99e30; // Mass of the sun
const double L_ref = 1.5e11; // 1 astronomical unit
const double T_ref = 3600.0*24.0*365.0*10.0; // 10 years

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

template <typename T> inline void vector_delete_elem(vector<T>& v, const std::size_t i){ // remove element i of vector
	if(i < v.size()){
		v(i) = v(v.size() - 1);
		v.resize(v.size() - 1);
	}
}

template <typename real_type> class simulator{
	typedef vector<real_type> vec;
	public:
		simulator() = default;
		simulator(int N, real_type flength):N(N), flength(flength){
			G = real_type(6.674e-11 * M_ref * T_ref*T_ref / pow(L_ref, 3.0));
			
			mass = scalar_vector<real_type>(N, 1.0);
			accel_x = vec(N);
			accel_y = vec(N);
			
			
			std::mt19937 rand_generator(std::random_device{}());
			std::normal_distribution<real_type> normal_dist(0.0, flength * 0.01);
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
		}
		
		void merge_point(const real_type K, bool profile){ // Merge two points if they collide. K is the "real length" per pixel constant (squared)
			std::clock_t begin;
			if(profile){
				begin = std::clock();
			}
			for(int i = 0; i < N; ++i){
				for(int j = 0; j < i ; ++j){
					real_type dist_x = point_x(j) - point_x(i);
					real_type dist_y = point_y(j) - point_y(i);
					real_type r2 = dist_x*dist_x + dist_y*dist_y;
					if(r2 < K * 2.0*std::sqrt(mass(i)*mass(j))+mass(i)+mass(j)){
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
				std::clock_t end{std::clock()};
				std::cout << "merge_point : " << 1000.0 * (end - begin) / CLOCKS_PER_SEC << " ms" << std::endl;
			}
		}
		
		void update_accel(bool profile){
			double begin = 0.0;
			if(profile){
				begin = omp_get_wtime();
			}
			#pragma omp parallel for schedule(static) num_threads(nthreads)
			for(int i = 0; i < N - 1; ++i){ // Compute acceleration of every particle i based on Newton law of gravitation

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
		
		void update_point(const real_type dt, bool profile){ // Update position of points using a Euler scheme (low precision). TODO : implement a more accurate scheme
			std::clock_t begin;
			if(profile){
				begin = std::clock();
			}
			point_x += dt*vel_x;
			point_y += dt*vel_y;
			
			vel_x += dt*accel_x;
			vel_y += dt*accel_y;
			if(profile){
				std::clock_t end{std::clock()};
				std::cout << "update_point : " << 1000.0 * (end - begin) / CLOCKS_PER_SEC << " ms" << std::endl;
			}
		}
		
		int N;
		real_type flength;
		real_type G;
		
		vec point_x;
		vec point_y;
		vec mass;
		
	private:
		vec accel_x;
		vec accel_y;
		vec vel_x;
		vec vel_y;
};
#endif
