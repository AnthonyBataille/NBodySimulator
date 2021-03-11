#ifndef GRAPHICS_H
#define GRAPHICS_H

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <SDL2/SDL.h>
#include <chrono>
#include <ctime>
#include "simulator.h"

using namespace boost::numeric::ublas;


void logSDLError(std::ostream& os, const std::string& msg){
	os << msg << " error: " << SDL_GetError() << std::endl;
}

SDL_Window* create_window(const std::string& title, const int& pos_x, const int& pos_y, const int& width, const int& height){
	SDL_Window* win = SDL_CreateWindow(title.c_str(), pos_x, pos_y, width, height, SDL_WINDOW_SHOWN);
	if(win == nullptr){
		logSDLError(std::cout, "SDL_CreateWindow");
	}
	return win;
}

SDL_Renderer* create_renderer(SDL_Window* win){
	SDL_Renderer* ren = SDL_CreateRenderer(win, -1, SDL_RENDERER_ACCELERATED | SDL_RENDERER_PRESENTVSYNC);
	if(ren == nullptr){
		logSDLError(std::cout, "SDL_CreateRenderer");
	}
	return ren;
}

void draw_circle(SDL_Renderer* ren, const SDL_Point& center, const int radius){
	if(radius > 0){
		int radius2 = radius*radius;
		SDL_Point pixels[1 + 4*radius + 4*radius2];
		int count = 0;
		for(int y = center.y - radius; y <= center.y + radius; ++y){
			for(int x = center.x - radius; x <= center.x + radius; ++x){
				int dx = x - center.x;
				int dy = y - center.y;
				if(dx*dx + dy*dy <= radius2){
					pixels[count].x = x;
					pixels[count].y = y;
					++count;
				}
			}
		}
		SDL_RenderDrawPoints(ren, pixels, count);
	}
}


class graphics{
	public:
		
		bool quit;
		
		int size;
	
		graphics() = default;
		
		graphics(const int size):quit(false), size(size), win(nullptr), ren(nullptr){
			if (SDL_Init(SDL_INIT_VIDEO) != 0){
				logSDLError(std::cout, "SDL_Init");
				quit = true;
			}
			
			win = create_window("NBodySimulator", SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, size, size);
			if(win == nullptr){
				quit = true;
			}

			ren = create_renderer(win);
			if(ren == nullptr){
				quit = true;
			}
		}
		
		~graphics(){
			if(ren != nullptr) SDL_DestroyRenderer(ren);
			if(win != nullptr) SDL_DestroyWindow(win);
			SDL_Quit();
		}
		
		void process_input(){
			while(SDL_PollEvent(&e)){
				if(e.type == SDL_QUIT){
					quit = true;
				}
			}
		}
		
		template <typename real_type> void render_point(const simulator<real_type>& s, bool profile){
			std::clock_t begin;
			if(profile){
				begin = std::clock();
			}
			SDL_SetRenderDrawColor(ren, 0xFF, 0xFF, 0xFF, 0xFF); // Clear screen (white)
			SDL_RenderClear(ren);
			
			SDL_SetRenderDrawColor(ren, 0x00, 0x00, 0xFF, 0xFF);
			for(int i = 0; i < s.N; ++i){
				draw_circle(ren, {int((s.point_x(i)/s.flength + 0.5)*size), int((s.point_y(i)/s.flength + 0.5)*size)}, int(std::sqrt(s.mass(i)))); // Draw every particle.
			}
			
			SDL_RenderPresent(ren);
			if(profile){
				std::clock_t end{std::clock()};
				std::cout << "render_point : " << 1000.0 * (end - begin) / CLOCKS_PER_SEC << " ms" << std::endl;
			}
		}
	private:
		SDL_Window* win;
		SDL_Renderer* ren;
		SDL_Event e;
};
#endif
