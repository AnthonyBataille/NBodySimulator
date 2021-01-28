#include "simulator.h"
#include "graphics.h"

int main(){
	int N = 0;
	float flength = 0.0;
	int size = 0;
	float dt = 0.0;
	while(N < 1){
		std::cout << "Number of particles: ";
		std::cin >> N;
	}
	while(flength <= 0.0){
		std::cout << "Field length: ";
		std::cin >> flength;
	}
	while(size < 1){
		std::cout << "Screen size (pixels): ";
		std::cin >> size;
	}
	while(dt <= 0.0){
		std::cout << "Time step (dt): ";
		std::cin >> dt;
	}
	simulator s{N, flength};
	graphics g{size};
	
	while(g.quit == false){
		g.process_input();
		s.merge_point(flength / size);
		s.update_accel();
		s.update_point(dt);
		g.render_point(s.point_x, s.point_y, s.mass, s.N, s.flength);
	}
	
	return 0;
}
