#ifndef GRAPHICS_H
#define GRAPHICS_H

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <SDL2/SDL.h>

using namespace boost::numeric::ublas;

SDL_Window* create_window(const std::string& title, const int& pos_x, const int& pos_y, const int& width, const int& height);
SDL_Renderer* create_renderer(SDL_Window* win);

class graphics{
	public:
		
		bool quit;
		
		int size;
	
		graphics();
		graphics(const int size);
		~graphics();
		void process_input();
		void render_point(const vector<float>& point_x, const vector<float>& point_y, const vector<float>& mass, const int N, const float flength);
	private:
		SDL_Window* win;
		SDL_Renderer* ren;
		SDL_Event e;
};

#endif
