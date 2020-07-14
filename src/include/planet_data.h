/*
 * planet_data.h
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
*NAMESPACE: solar_system
* 
*DESCRIPTION: Simple chunk of data for solar system planets. Main purpose is 
*for benchmark tests. Units given in SI standard system(kg, m, m/sec)
*
*/
namespace solar_system
{

constexpr double G = 6.673e-11;

// location(x, y, z), mass, radius, velocity(x, y, z), full_name
static body sun{ ORIGIN, 1.989e30, 695700000, ORIGIN, "Sun"};

static body mercury{{-46000000000, 0.0, 0.0}, 3.3011e23, 2439700, {0.0,-58980,0.0}, "Mercury"};

static body venus{{-107480000000,0.0,  0.0}, 4.8675e24, 6051800, {0.0,-35260,0.0}, "Venus"};

static body earth{{-147095000000, 0.0, 0.0}, 5.98e24, 6371000, {0.0, -30300, 0.0 }, "Earth"};

static body mars{{-206620000000, 0.0, 0.0}, 6.4171e23, 3389500, {0.0, -26500, 0.0}, "Mars"};

static body jupiter{{-740520000000, 0.0, 0.0}, 1898.19e24, 71492000,{0.0, -13720, 0.0}, "Jupiter"};

static body saturn{{-1352550000000, 0.0, 0.0},  568.34e24, 54364000,{0.0, -10180, 0.0}, "Saturn"};

static body uranus = {{-2741300000000, 0.0, 0.0}, 86.813e24, 24973000,{0.0, -7110, 0.0 }, "Uranus"};

static body neptune = {{-4444450000000, 0.0, 0.0}, 102.413e24, 24341000,{0.0, -5500, 0.0}, "Neptune"};

static body pluto = {{5906376272000, 0.0, 0.0}, 1.30900e22, 1188300,{0.0, -4760, 0.0}, "Pluto"};

static body test_object = {{1, 0, 0}, 1, 1, {0, 0.5, 0}, "0B545"};
}
