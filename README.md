# NBodySimulator
A simple 2D N-Body Simulator. Simulates the motion of N particles under the influence of Newton's Gravitation law.

Two particles merge if they collide.

Acceleration is computed directly for each particles (no octree aproximation, possible future improvement).

Integration is made with an explicit Euler scheme or a leapfrog scheme.

## Python script (draft)
The file python/simulator.py a working but very slow version of the simulator. Usage of Numba can improve performance but doesn't seem to work on every platform. 

The python script can simulate two initial configurations modes:

- Default: points and velocity are random with several distributions possible.
- Solar: simulates N-1 particles rotating around a massive one (e.g. sun, planet, ...).

## C++ code
The C++ code is equivalent to the python script with much better performances.

It uses the Boost::uBLAS library for speedup in computations and SDL2 library for rendering graphics.
Default and solar system modes are also supported.

Example of usage : `./NBodySimulator -n 2000 -f 2000 -w 600 -t 0.01 -s`

## Required Libraries
- libsdl2-2.0-0
- libsdl2-dev
- Boost::uBLAS (included)
