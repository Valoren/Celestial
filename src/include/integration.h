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

/*
*NAMESPACE: Integration_Algorithms
* 
*DESCRIPTION: Defines basic integration algorithms for the integration.
*
*/
namespace Integration_Algorithms
{

class Integrator
{
	public:
	
	point calculate_single_body_acceleration(int);
	void compute_velocity();
	void update_location();
	void compute_gravity_step();
    std::vector<body> &get_bodies();
    std::vector<body> bodies;
	

};
}


