#ifndef SIMULATOR_H
#define SIMULATOR_H
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>

using namespace boost::numeric::ublas;

template <typename T> inline void vector_delete_elem(vector<T>& v, const std::size_t i);

class simulator{
	public:
		simulator();
		simulator(int N, float flength);
		
		void update_accel();
		void merge_point(const float K);
		void update_point(const float dt);
		
		int N;
		float flength;
		float G;
		
		vector<float> point_x;
		vector<float> point_y;
		vector<float> mass;
		
	private:
		vector<float> accel_x;
		vector<float> accel_y;
		vector<float> vel_x;
		vector<float> vel_y;
};
#endif
