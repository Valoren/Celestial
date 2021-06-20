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

//Functions declarations used on Piet-Hut Hermite integrator implementation
void correct_step(real pos[][NDIM], real vel[][NDIM],
                  const real acc[][NDIM], const real jerk[][NDIM],
                  const real old_pos[][NDIM], const real old_vel[][NDIM],
                  const real old_acc[][NDIM], const real old_jerk[][NDIM],
                  int n, real dt);
void evolve(const real mass[], real pos[][NDIM], real vel[][NDIM],
            int n, real & t, real dt_param, real dt_dia, real dt_out,
            real dt_tot, bool init_out, bool x_flag);
void evolve_step(const real mass[], real pos[][NDIM], real vel[][NDIM],
                 real acc[][NDIM], real jerk[][NDIM], int n, real & t,
                 real dt, real & epot, real & coll_time);
void get_acc_jerk_pot_coll(const real mass[], const real pos[][NDIM],
                           const real vel[][NDIM], real acc[][NDIM],
                           real jerk[][NDIM], int n, real & epot,
                           real & coll_time);
void get_snapshot(real mass[], real pos[][NDIM], real vel[][NDIM], int n);
void predict_step(real pos[][NDIM], real vel[][NDIM],
                  const real acc[][NDIM], const real jerk[][NDIM],
                  int n, real dt);
void put_snapshot(const real mass[], const real pos[][NDIM],
                  const real vel[][NDIM], int n, real t);
bool read_options(int argc, char *argv[], real & dt_param, real & dt_dia,
                  real & dt_out, real & dt_tot, bool & i_flag, bool & x_flag);
void write_diagnostics(const real mass[], const real pos[][NDIM],
                       const real vel[][NDIM], const real acc[][NDIM],
                       const real jerk[][NDIM], int n, real t, real epot,
                       int nsteps, real & einit, bool init_flag,
                       bool x_flag);


static const double dt = 0.00000001;

static void record_state(std::vector<body>& bodies)
{
    for (auto body_iterator = bodies.begin(); body_iterator != bodies.end(); *body_iterator++)
    {
        body_iterator->locations.push_back(body_iterator->location);
    }
}

static void output_states(const std::vector<body>& body_locations)
{

    for (auto body_iterator = body_locations.begin(); body_iterator != body_locations.end(); *body_iterator++)
    {
        std::ofstream f;
        f.open(body_iterator->name + ".dat");
        f << body_iterator->name << std::endl;
        f.close();
    }

    for (auto body_iterator = body_locations.begin(); body_iterator != body_locations.end(); *body_iterator++)
    {
        std::ofstream f;
        f.open(body_iterator->name + ".dat", std::ofstream::out | std::ofstream::app);
        for (auto location = body_iterator->locations.begin(); location < body_iterator->locations.end(); *location++)
        {
            f << location->x << ","
              << location->y << ","
              << location->z << std::endl;
        }
        f.close();
    }
}

/*
*NAMESPACE: Orbit integration
*
*DESCRIPTION: Defines basic integration algorithms for the integration.
*
*/
namespace Orbit_integration {
/*
*CLASS: Integrator
*
*DESCRIPTION: Encapsulates all basic integration algorithms.
*
*/
    class Integrator {
    public:
        virtual void compute_gravity_step() = 0;
        virtual std::vector<body> &get_bodies() = 0;
    };

    class Euler : virtual Integrator {
    public:
        Euler(std::vector<body> bodies, double time_step = 1) :
                m_bodies(bodies),
                m_time_step(time_step) {};

        std::vector<body> &get_bodies() { return m_bodies; };

        void compute_gravity_step();

    private:
        point calculate_single_body_acceleration(int);

        void compute_velocity();

        void update_location();

    protected:
        std::vector<body> m_bodies;
        double m_time_step;
    };

    class RK4 : virtual public Integrator {
    public:
        RK4(std::vector<body> bodies, double time_step = 1) :
                m_bodies(std::move(bodies)),
                m_time_step(time_step) {};

        std::vector<body> &get_bodies() { return m_bodies; };

        void compute_gravity_step();

    private:
        point calculate_single_body_acceleration(int);

        point partial_step(point &, point &, double);

        void compute_velocity();

        void update_location();

    protected:
        std::vector<body> m_bodies;
        double m_time_step;
    };
}

/*
*NAMESPACE: Two_Body_Algorithms
* 
*DESCRIPTION: Defines special algorithms for use with 2 body 
*systems only.
*
*/
namespace two_body_algorithms{
    void f_and_g();
    void euler_forward(double dt);
    void leapfrog(double dt, double integration_time);
    double e_out(point a, point b); //Final total energy
    double ekin(point p1);
    double epot(point p2);
}

/*
*NAMESPACE: revised_algorithms
*
*DESCRIPTION: Alternative implementations
*
*/
namespace revised_algorithms{
  
   void euler_def();
   void hermite();
   void RK4();
}
