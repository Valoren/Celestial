/*
 * integration.h
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

#pragma once

#include "structures.h"
#include "planet_data.h"

using namespace solar_system;

/*
*CLASS: Integrator
*
*DESCRIPTION: Encapsulates all basic integration algorithms.
*
*/
class Integrator
{
public:

    static void calculate_single_body_acceleration(int);
    void compute_velocity();
    void update_location();
    void compute_gravity_step();
    std::vector<body> &get_bodies();
    std::vector<body> bodies;


};

void  Integrator::calculate_single_body_acceleration(int) {
    std::cout << "Hola" << std::endl;

}

/*
*NAMESPACE: Integration_Algorithms
* 
*DESCRIPTION: Defines basic integration algorithms for the integration.
*
*/
namespace General_Integration_Algorithms
{

}


/*
*NAMESPACE: Two_Body_Algorithms
* 
*DESCRIPTION: Defines special algorithms for use with 2 body 
*systems only.
*
*/
namespace Two_Body_Algorithms {

    void F_and_G(){


        point calculate_single_body_acceleration(int);

        void compute_velocity();

        void update_location();

        void compute_gravity_step();

        std::vector<body> &get_bodies();

        std::vector<body> bodies;
    }

    void euler_forward(){
        double r[3], v[3], a[3];
        double dt = 0.00000001;
        r[0] = test_object.location.x;
        r[1] = test_object.location.y;
        r[2] = test_object.location.z;
        v[0] = test_object.velocity.x;
        v[1] = test_object.velocity.y;
        v[2] = test_object.velocity.z;
        double dt_out = 0.01;
        double t_out = dt_out;
        for (double t = 0; t < 10; t += dt) {
            double r2 = r[0] * r[0] + r[1] * r[1] + r[2] * r[2];
            for (int k = 0; k < 3; k++)a[k] = -r[k] / (r2 * sqrt(r2));
            for (int k = 0; k < 3; k++) {
                r[k] += v[k] * dt;
                v[k] += a[k] * dt;
            }
            if (t >= t_out) {
                std::cout << r[0] << " " << r[1] << " " << r[2] << " ";
                std::cout << v[0] << " " << v[1] << " " << v[2] << std::endl;
                t_out += dt_out;
            }
        }

    }
}
