/*
 * main.cpp
 * 
 * Copyright 2019 Miquel Bernat Laporta i Granados 
 * <mlaportaigranados@gmail.com>
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA 02110-1301, USA.
 * 
 * 
 */

/*
Matplotlib link instructions:
Compile with folowing flags:
-I/usr/include/python2.7 -lpython2.7
*/

#include <iostream>
#include <cmath>
#include <omp.h>
#include "include/structures.h"
#include "include/integration.h"
#include "include/planet_data.h"
#include "include/menu.h"
#include "include/benchmark.h"
#include "parser.h"
#include "astro_constants.h"
//#include "include/astro_epochs.h"
//#include "matplotlibcpp.h"

std::vector<body> bodies;

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

void Process_Selection_Three(double timestep){
    two_body_algorithms::f_and_g();
}

void Process_Selection_One(double timestep){
    
        Orbit_integration::Euler orbit(bodies, 10);
        run_simulation(orbit, int(1e4), 1);
}

void Process_Selection_Two(double timestep){
    
        Orbit_integration::RK4 orbit(bodies, 10);
        run_simulation(orbit, int(1e4), 1);
}

int main(int argc, char *argv[]){

    //Using solar system data in planet_data.h for benchmarking
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
        //parse_file();
        //parse_data(argv[1]);


    number_of_cores();
    //spawn_menu();
    {
        Timer timer;
    }
    std::cout << "Execution terminated!!" << std::endl;

    return 0;
}
