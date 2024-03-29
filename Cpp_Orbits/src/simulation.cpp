#include "orbit_integration.h"
#include <omp.h>
#include <iostream>

template <typename Integrator>
void run_simulation(Integrator integrator, int iterations, int report_frequency)
{
	for (auto i = 0; i < iterations; i++)
	{
		if (i % report_frequency == 0)
			record_state(integrator.get_bodies());
		integrator.compute_gravity_step();
	}
	output_states(integrator.get_bodies());
}

int main(int argc, char *argv[])
{
	std::vector<body> bodies;
	/*
	bodies.push_back(solar_system::sun);
	bodies.push_back(solar_system::mercury);
	bodies.push_back(solar_system::venus);
	bodies.push_back(solar_system::earth);
	bodies.push_back(solar_system::mars);
	bodies.push_back(solar_system::saturn);
	bodies.push_back(solar_system::jupiter);
	bodies.push_back(solar_system::uranus);
	bodies.push_back(solar_system::neptune);
	bodies.push_back(solar_system::pluto);
	*/
    bodies.push_back(benchmarking::m1);
	bodies.push_back(benchmarking::m2);
	#pragma omp parallel
	{

	  Orbit_integration::RK4 orbit(bodies, 0.01);
	  run_simulation(orbit, (int)1e4, 1);
	}

	std::cout << "Done" << std::endl;
}
