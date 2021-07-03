/*
 * integration.cpp
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

#include "include/integration.h"
#include "structures.h"

#define NUMBER_OF_STEPS 100

using namespace solar_system;

/*-----------------------------------------------------------------------------
 *PROCEDURE:  get_snapshot  
 *
 *DESCRIPTION: reads a single snapshot from the input stream cin.
 * note: in this implementation, only the particle data are read in, and it
 *        is left to the main program to first read particle number and time
 *
 *RETURNS: -
 *-----------------------------------------------------------------------------
 */

void get_snapshot(real mass[], real pos[][NDIM], real vel[][NDIM], int n)
{
    for (int i = 0; i < n ; i++){
        std::cin >> mass[i];                       // mass of particle i
        for (int k = 0; k < NDIM; k++)
            std::cin >> pos[i][k];                 // position of particle i
        for (int k = 0; k < NDIM; k++)
            std::cin >> vel[i][k];                 // velocity of particle i
    }
}

/*-----------------------------------------------------------------------------
 *PROCEDURE:  put_snapshot  
 *
 *DESCRIPTION: writes a single snapshot on the output stream cout.
 * note: unlike get_snapshot(), put_snapshot handles particle number and time
 *
 *RETURNS: -
 *-----------------------------------------------------------------------------
 */
void put_snapshot(const real mass[], const real pos[][NDIM],
                  const real vel[][NDIM], int n, real t)
{
    std::cout.precision(16);                       // full double precision

    std::cout << n << std::endl;                        // N, total particle number
    std::cout << t << std::endl;                        // current time
    for (int i = 0; i < n ; i++){
        std::cout << mass[i];                      // mass of particle i
        for (int k = 0; k < NDIM; k++)
            std::cout << ' ' << pos[i][k];         // position of particle i
        for (int k = 0; k < NDIM; k++)
            std::cout << ' ' << vel[i][k];         // velocity of particle i
        std::cout << std::endl;
    }
}


/*----------------------------------------------------------------------------
 *PROCEDURE: calculate_single_body_acceleration
 *
 *DESCRIPTION: Calculates acceleration vector for a single object for Euler algorithm
 *
 *RETURNS: Acceleration vector
 *----------------------------------------------------------------------------
 */
point Orbit_integration::Euler::calculate_single_body_acceleration(int body_index)
{
    point acceleration{ 0, 0, 0 };
    body target_body = m_bodies[body_index];

    int index = 0;

    for (auto external_body = m_bodies.begin(); external_body != m_bodies.end(); *external_body++, index++)
    {
        if (index != body_index)
        {
            double r = (pow((target_body.location.x - external_body->location.x), 2) + pow((target_body.location.y - external_body->location.y), 2) + pow((target_body.location.z - external_body->location.z), 2));
            r = sqrt(r);
            auto tmp = G_const * external_body->mass / (r*r*r);
            acceleration = acceleration + (external_body->location - target_body.location) * tmp;
        }
    }
    return acceleration;
}

/*
 *PROCEDURE: compute_velocity
 *
 *DESCRIPTION: Calculates velocity vector for a single object for Euler algorithm
 *
 *RETURNS: -
 *
 */
void Orbit_integration::Euler::compute_velocity()
{
    for (int i = 0; i < m_bodies.size(); i++)
    {
        point acceleration = Orbit_integration::Euler::calculate_single_body_acceleration(i);
        m_bodies[i].velocity += acceleration * m_time_step;
    }
}

/*
 *PROCEDURE: update_location
 *
 *DESCRIPTION: Updates position vector to actual step for Euler algorithm
 *
 *RETURNS: -
 *
 */
void Orbit_integration::Euler::update_location()
{
    for (auto target_body = m_bodies.begin(); target_body != m_bodies.end(); *target_body++)
    {
        target_body->location += target_body->velocity * m_time_step;
    }
}

/*
 *PROCEDURE: compute_gravity_step
 *
 *DESCRIPTION: Merges compute_velocity and update_location methods into one. Applied to Euler algorithm.
 *
 *RETURNS: -
 *
 */
void Orbit_integration::Euler::compute_gravity_step()
{
    compute_velocity();
    update_location();
}

/*
 *PROCEDURE: calculate_single_body_acceleration
 *
 *DESCRIPTION: Calculates acceleration vector for a single object for RK4 algorithm
 *
 *RETURNS: Acceleration vector
 *
 */
point Orbit_integration::RK4::calculate_single_body_acceleration(int body_index)
{
    point acceleration{ 0, 0, 0 };
    point velocity_update{ 0, 0, 0 };
    point location_update{ 0, 0, 0 };
    body target_body = m_bodies[body_index];

    int index = 0;
    for (auto external_body = m_bodies.begin(); external_body != m_bodies.end(); *external_body++, index++)
    {
        if (index != body_index)
        {
            point k1{ 0, 0, 0 };
            point k2{ 0, 0, 0 };
            point k3{ 0, 0, 0 };
            point k4{ 0, 0, 0 };

            double r = (pow((target_body.location.x - external_body->location.x), 2) +
                        pow((target_body.location.y - external_body->location.y), 2) +
                        pow((target_body.location.z - external_body->location.z), 2));

            r = sqrt(r);

            auto tmp = G_const * external_body->mass / (r*r*r);

            //k1 - acceleration at current location
            k1 = (external_body->location - target_body.location) * tmp;

            //k2 - acceleration 0.5 timesteps in the future based on k1 acceleration value
            velocity_update = partial_step(target_body.velocity, k1, 0.5);
            location_update = partial_step(target_body.location, velocity_update, 0.5);
            k2 = (external_body->location - location_update) * tmp;

            //k3 acceleration 0.5 timesteps in the future using k2 acceleration
            velocity_update = partial_step(target_body.velocity, k2, 0.5);
            location_update = partial_step(target_body.location, velocity_update, 0.5);
            k3 = (external_body->location - location_update) * tmp;

            //k4 - location 1 timestep in the future using k3 acceleration
            velocity_update = partial_step(target_body.velocity, k3, 1);
            location_update = partial_step(target_body.location, velocity_update, 1);
            k4 = (external_body->location - location_update) * tmp;

            acceleration += (k1 + k2 * 2 + k3 * 2 + k4) / 6;
        }
    }
    return acceleration;
}

/*
 *PROCEDURE: partial_step
 *
 *DESCRIPTION: Calculates partial step constants for RK4 algorithm
 *
 *RETURNS: f.x, f.y and f.z values
 *
 */
point Orbit_integration::RK4::partial_step(point &f, point &df, double scale)
{
    return point{
            f.x + df.x * m_time_step * scale,
            f.y + df.y * m_time_step * scale,
            f.z + df.z * m_time_step * scale
    };
}

/*
 *PROCEDURE: compute_velocity
 *
 *DESCRIPTION: Calculates velocity vector for a single object for RK4 algorithm
 *
 *RETURNS: -
 *
 */
void Orbit_integration::RK4::compute_velocity()
{
    for (int i = 0; i < m_bodies.size(); i++)
    {
        point acceleration = Orbit_integration::RK4::calculate_single_body_acceleration(i);
        m_bodies[i].velocity += acceleration * m_time_step;
    }
}

/*
 *PROCEDURE: update_location
 *
 *DESCRIPTION: Updates position vector for all objects for RK4 algorithm
 *
 *RETURNS: -
 *
 */
void Orbit_integration::RK4::update_location()
{
    for (auto target_body = m_bodies.begin(); target_body != m_bodies.end(); *target_body++)
    {
        target_body->location += target_body->velocity * m_time_step;
    }
}

/*
 *PROCEDURE: compute_gravity_step
 *
 *DESCRIPTION: Merges compute_velocity and update_location methods into one. Applied to RK4 algorithm.
 *
 *RETURNS: -
 *
 */
void Orbit_integration::RK4::compute_gravity_step()
{
    compute_velocity();
    update_location();
}

/*
 *PROCEDURE: ekin
 *
 *DESCRIPTION: Calculates kinetic energy on m=1 scale
 *
 *RETURNS: double
 *
 */
double two_body_algorithms::ekin(point p1){
    return 0.5 * norm(p1) * norm(p1);
}

/*
 *PROCEDURE:epot
 *
 *DESCRIPTION: Calculates potential energy on m=1 scale
 *
 *RETURNS: double
 *
 */
double two_body_algorithms::epot(point p2){
    return 1.0 / norm(p2);
}

/*
 *PROCEDURE: e_out
 *
 *DESCRIPTION: Calculates total energy on the system
 *
 *RETURNS: double
 *
 */
double two_body_algorithms::e_out(point a, point b){
    return ekin(a)+epot(b);
}

//NOT PROPERLY TESTED
void two_body_algorithms::f_and_g(){

        point r{};
        point v{};
        double u, z, q, f, g;
        double F2, F3, F4;
        double G3, G4;
        double f_dot, g_dot;
        double mu = test_object.mass + 1;
        r.x = f_g_test.location.x;
        r.y = f_g_test.location.y;
        r.z = f_g_test.location.z;
        v.x = f_g_test.velocity.x;
        v.y = f_g_test.velocity.y;
        v.z = f_g_test.velocity.z;
        for(int k = 0; k < NUMBER_OF_STEPS; k++){
            double t = 0.0;
            u = mu / (norm(r) * norm(r) * norm(r));
            z = scalar(r, v)/(norm(r) * norm(r));
            q = (scalar(r, r)/(norm(r) * norm(r)) - u);
            F2 = -u/2;
            F3 = (u * z)/2;
            F4 = (u/24)*(3 * q - 15 * z *z + u);
            G3 = -u/6;
            G4 = (u * z)/4;
            f = 1 + F2*t*t + F3*t*t*t + F4*t*t*t*t;
            g = t + G3*t*t*t + G4*t*t*t*t;
            f_dot = 2*F2*t + 3*F3*t*t + 4*F4*t*t*t;
            g_dot = 1 + 3*G3*t*t + 4*G4*t*t*t;
            r.x += f*r.x + g*v.x;
            r.y += f*r.y + g*v.y;
            r.z += f*r.z + g*v.z;
            v.x += f_dot*r.x + g_dot*v.x;
            v.y += f_dot*r.y + g_dot*v.y;
            v.z += f_dot*r.z + g_dot*v.z;
            std::cout << r.x << " " << r.y << " " << r.z << " ";
            std::cout << v.x << " " << v.y << " " << v.z << std::endl;
        }
}

void two_body_algorithms::euler_forward(double dt){
    double r[3], v[3], a[3];
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

/*
*PROCEDURE: leapfrog
*
*DESCRIPTION: Integrates postion using leapfrog algorithm
*
*/
void two_body_algorithms::leapfrog(double dt, double integration_time){
    point r{};
    point v{};
    point a{};
    double dt_out = 0.01;
    double t_out = dt_out;
    r.x = f_g_test.location.x;
    r.y = f_g_test.location.y;
    r.z = f_g_test.location.z;
    v.x = f_g_test.velocity.x;
    v.y = f_g_test.velocity.y;
    v.z = f_g_test.velocity.z;
    ekin(v);
    epot(v);
    double r_init=norm(r);
    double r2= r.x * r.x + r.y * r.y + r.z * r.z; //Initial r without the norm(r_init)^2 
    
    //Initial accelerations
    a.x= -r.x / (r_init* r2);
    a.y= -r.y / (r_init* r2);
    a.z= -r.x / (r_init* r2);

    //Main integration loop
    for(double t=0; t < integration_time; t += dt){
        v.x += 0.5 * a.x * dt;  
        v.y += 0.5 * a.y * dt;
        v.z += 0.5 * a.z * dt;
        r.x += v.x * dt;
        r.y += v.y * dt;
        r.z += v.z * dt; 
        r2= r.x * r.x + r.y * r.y + r.z * r.z;
        a.x= -r.x / (sqrt(r2)* r2);
        a.y= -r.y / (sqrt(r2)* r2);
        a.z= -r.x / (sqrt(r2)* r2);
        v.x += 0.5 * a.x * dt;  
        v.y += 0.5 * a.y * dt;
        v.z += 0.5 * a.z * dt;
        if(t >= dt_out){
            std::cout << r.x << " " << r.y << " " << r.z << " ";
            std::cout << v.x << " " << v.y << " " << v.z << std::endl;
            t_out += dt_out;
        }
    }
    ekin = 0.5 * (v.x * v.x + v.y * v.y + v.z * v.z);
    epot = -1.0 / norm(r);
    eout
    std::cout<< "Final total energy:" << 



}
