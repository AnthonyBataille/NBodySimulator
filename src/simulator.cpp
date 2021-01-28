#define NDEBUG
#include <cassert>

//#define SIMULATOR_PROFILE // Uncomment to print execution time of functions

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <random>
#include <chrono>
#include <ctime>

#include "simulator.h"

using namespace boost::numeric::ublas;
// Reference units in S.I
const double M_ref = 1.99e30; // Mass of the sun
const double L_ref = 1.5e11; // 1 astronomical unit
const double T_ref = 3600.0*24.0*365.0*10.0; // 10 years

// Taken from http://www.crystalclearsoftware.com/cgi-bin/boost_wiki/wiki.pl?Examples_-_How_To_Extend_UBLAS
// (op v) [i] = op( v [i] )
template<class OP, class E> 
BOOST_UBLAS_INLINE
typename vector_unary_traits<E, OP>::result_type
apply_to_all(const vector_expression<E> &e, const OP& op = OP()){
	typedef typename vector_unary_traits<E, OP>::expression_type expression_type;
	return expression_type (e());
}

namespace functor{
	template <class T> class sqrt{ // apply square root to all elements
	public:
		typedef T value_type;
		typedef T result_type;
		sqrt() { }
    
		static BOOST_UBLAS_INLINE result_type apply(const value_type& x){
		return std::sqrt(x);
		}
	};
	template <class T> class inv{ // inverse all elements
	public:
		typedef T value_type;
		typedef T result_type;
		inv() { }
    
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

simulator::simulator() = default;

simulator::simulator(int N, float flength):N(N), flength(flength){
	G = float(6.674e-11 * M_ref * T_ref*T_ref / pow(L_ref, 3.0));
	
	mass = scalar_vector<float>(N, 1.0);
	accel_x = vector<float>(N);
	accel_y = vector<float>(N);
	
	
	std::mt19937 rand_generator(std::random_device{}());
	std::normal_distribution<float> normal_dist(0.0, flength * 0.01);
	std::uniform_real_distribution<float> uni_dist(-flength / 2.0, flength / 2.0);
	
	vel_x = vector<float>(N);
	vel_y = vector<float>(N);
	for(int i = 0; i < N; ++i){
		vel_x(i) = normal_dist(rand_generator); // generate velocities according to a normal distribution
		vel_y(i) = normal_dist(rand_generator);
	}
	point_x = vector<float>(N);
	point_y = vector<float>(N);
	for(int i = 0; i < N; ++i){
		point_x(i) = uni_dist(rand_generator); // generate points according to a uniform distribution in the field
		point_y(i) = uni_dist(rand_generator);
	}
}

void simulator::update_accel(){
#ifdef SIMULATOR_PROFILE
	std::clock_t begin{std::clock()};
#endif
	for(int i = 0; i < N - 1; ++i){ // Compute acceleration of every particle i based on Newton law of gravitation
		vector<float> p_x{point_x};
		vector<float> p_y{point_y};

		vector_delete_elem<float>(p_x, i);
		vector_delete_elem<float>(p_y, i);

		vector<float> r_x{p_x - scalar_vector<float>(N - 1, point_x(i))};
		vector<float> r_y{p_y - scalar_vector<float>(N - 1, point_y(i))};
		
		vector<float> r2_inv{element_prod(r_x, r_x) + element_prod(r_y, r_y)};
		r2_inv = apply_to_all(r2_inv, functor::inv<float>());
		vector<float> r_inv{apply_to_all(r2_inv, functor::sqrt<float>())};
		vector<float> r32_inv{element_prod(r2_inv, r_inv)};
		
		r_x = element_prod(r_x, r32_inv);
		r_y = element_prod(r_y, r32_inv);
		

		vector<float> m{mass};
		vector_delete_elem<float>(m, i);

		accel_x(i) = G * inner_prod(m, r_x);
		accel_y(i) = G * inner_prod(m, r_y);
	}
#ifdef SIMULATOR_PROFILE
	std::clock_t end{std::clock()};
	std::cout << "update_accel : " << 1000.0*(end - begin) / CLOCKS_PER_SEC << " ms" << std::endl;
#endif
}

void simulator::merge_point(const float K){ // Merge two points if they collide. K is the "real length" per pixel constant
#ifdef SIMULATOR_PROFILE
	std::clock_t begin{std::clock()};
#endif
	for(int i = 0; i < N; ++i){
		for(int j = 0; j < i ; ++j){
			float dist_x = point_x(j) - point_x(i);
			float dist_y = point_y(j) - point_y(i);
			float r = std::sqrt(dist_x*dist_x + dist_y*dist_y);
			if(r < K * (std::sqrt(mass(i)) + std::sqrt(mass(j)))){
				int& k = mass(i) < mass(j) ? i : j;
				int& l = mass(i) < mass(j) ? j : i;
				float sum_mass = mass(k) + mass(l);
				vel_x(l) = (mass(k)*vel_x(k) + mass(l)*vel_x(l)) / sum_mass;
				vel_y(l) = (mass(k)*vel_y(k) + mass(l)*vel_y(l)) / sum_mass;

				mass(l) += mass(k);
				vector_delete_elem<float>(mass, k);
				vector_delete_elem<float>(accel_x, k);
				vector_delete_elem<float>(accel_y, k);
				vector_delete_elem<float>(vel_x, k);
				vector_delete_elem<float>(vel_y, k);
				vector_delete_elem<float>(point_x, k);
				vector_delete_elem<float>(point_y, k);
				--N;
				--k;
				--l;
			}
		}
	}
#ifdef SIMULATOR_PROFILE
	std::clock_t end{std::clock()};
	std::cout << "merge_point : " << 1000.0*(end - begin) / CLOCKS_PER_SEC << " ms" << std::endl;
#endif
}


void simulator::update_point(const float dt){ // Update position of points using a Euler scheme (low precision). TODO : implement a more accurate scheme
#ifdef SIMULATOR_PROFILE
	std::clock_t begin{std::clock()};
#endif
	point_x += dt*vel_x;
	point_y += dt*vel_y;
	
	vel_x += dt*accel_x;
	vel_y += dt*accel_y;
#ifdef SIMULATOR_PROFILE
	std::clock_t end{std::clock()};
	std::cout << "update_point : " << 1000.0*(end - begin) / CLOCKS_PER_SEC << " ms" << std::endl;
#endif
}

