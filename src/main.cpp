#include <stdlib.h>
#include <string.h>
#include <getopt.h>

#include "simulator.h"
#include "graphics.h"


using real_type = float;

void usage (const char* path) // Options parsing taken from https://riptutorial.com/c/example/30858/using-gnu-getopt-tools
{
	const char *basename = strrchr(path, '/');
	basename = basename ? basename + 1 : path;

	std::cout << "usage: " << basename << " [OPTION]" << std::endl;
	std::cout << "  -h, --help\t\t" <<
				"Print this help and exit." << std::endl;
	std::cout << "  -n, --particles-number NUMBER\t\t" <<
				"The number of particles (N) in the simulation." << std::endl;
	std::cout << "  -f, --field-length LENGTH\t\t" <<
				"The length of the simulation in unit of length." << std::endl;
	std::cout << "  -w, --window-size SIZE\t\t" <<
				"The size of the output window in pixels." << std::endl;
	std::cout << "  -t, --timestep TIME\t\t" <<
				"The timestep for numerical integration, in unit of time." << std::endl;
	std::cout << "  -p, --profile\t\t" <<
				"Enable code profiling. Print the duration of each computing step." << std::endl;
	std::cout << "  -s, --solar-system MASS\t\t" <<
				"Initialize particles to form a solar system with a central particle (star) of mass MASS relative to the unit mass." << std::endl;
}

int main(int argc, char** argv){
	if(argc < 9){
		usage(argv[0]);
		return 1;
	}
	
	bool profile = false; // Set to true using -p to profile code.
	int N = 0;
	real_type flength = 0.0;
	int size = 0;
	real_type dt = 0.0;

	int opt = 0;
	int help_flag = 0;
	bool solar_system = false;
	real_type solar_mass = 1.0;

	struct option longopts[] = {
		{"help", no_argument, &help_flag, 0},
		{"particles-number", required_argument, nullptr, 'n'},
		{"field-length", required_argument, nullptr, 'f'},
		{"window-size", required_argument, nullptr, 'w'},
		{"timestep", required_argument, nullptr, 't'},
		{"profile", no_argument, nullptr, 'p'},
		{"solar-system", required_argument, nullptr, 's'},
		{0}
	};
	
	while(true){
		opt = getopt_long(argc, argv, "hn:f:w:t:ps:", longopts, 0);
		if (opt == -1){
			break;
		}
		switch(opt){
			case 'h':
				help_flag = 1;
				break;
			case 'n':
				N = atoi(optarg);
				break;
			case 'f':
				flength = real_type(atof(optarg));
				break;
			case 'w':
				size = atoi(optarg);
				break;
			case 't':
				dt = real_type(atof(optarg));
				break;
			case 'p':
				profile = true;
				break;
			case 's':
				solar_system = true;
				solar_mass = real_type(atof(optarg));
				break;
			case '?':
				usage(argv[0]);
				return 1;
			default:
				break;
		}
	}
	
	if(help_flag){
		usage(argv[0]);
		return 0;
	}
	
	real_type K = flength / size;
	
	simulator<real_type> s{N, flength};
	if(solar_system){
		s.init_solar_system(solar_mass);
	}
	
	graphics g{size};
	//euler_integrator<real_type> euler(dt, profile);
	leapfrog_integrator<real_type> leapfrog(dt, profile, s);
	
	while(g.quit == false){
		g.process_input();
		s.merge_point(K, profile);
		s.update_point(leapfrog, profile);
		//std::cout << s.compute_energy(profile) << std::endl;
		g.render_point<real_type>(s, profile);
	}
	
	return 0;
}
