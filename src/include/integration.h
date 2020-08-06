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
#define NUMBER_OF_STEPS 100

using namespace solar_system;

static const double dt = 0.00000001;

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
                m_bodies(bodies),
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
namespace Two_Body_Algorithms
{
    class Integrator {
    public:
        virtual void compute_gravity_step() = 0;
        virtual std::vector<body> &get_bodies() = 0;
    };

    class Leapfrog: virtual Integrator{
        Leapfrog(std::vector<body> bodies, double time_step = 1) :
        m_bodies(bodies),
        m_time_step(time_step) {};

        std::vector<body> &get_bodies() { return m_bodies; };

        void compute_gravity_step();

    private:
        point calculate_single_body_acceleration(int){
            point acceleration{ 0, 0, 0 };
            body target_body = [test_object, origin];

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

        void compute_velocity();

        void update_location();

    protected:
        std::vector<body> m_bodies;
        double m_time_step;

    };

}
