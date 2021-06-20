/*
 * structures.h
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

#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <iostream>
#include <algorithm>
#include <iomanip>

const int NDIM = 3;

/*
 * 	Struct: Point
 * 
 *  Members: Coordinates, basic vector operations and ">>" operator overload
 */
struct point
{
	double x;
	double y;
	double z;

	point operator*(point a) const
	{
		return point{ x*a.x, y*a.y, z*a.z };
	}

	point operator*(double a) const
	{
		return point{ x*a, y*a, z*a };
	}

	point operator/(point a) const
	{
	    return point{ x/a.x, y/a.y, z/a.z };
	}

	point operator/(double a) const
	{
		return point{ x / a, y / a, z / a };
	}

	point operator+(point a) const
	{
		return point{ x +a.x, y + a.y, z + a.z };
	}

	point operator+=(const point& a)
	{
		x += a.x;
		y += a.y;
		z += a.z;
		return point{ x, y, z };
	}

	point operator-(point a) const
	{
		return point{ x - a.x, y - a.y, z - a.z };
	}
};

const point ORIGIN{0.0, 0.0, 0.0};

/*
 *PROCEDURE: norm
 * 
 *DESCRIPTION: Calculates norm of a given cartesian Vector
 * 
 * RETURNS: norm of Vector(double)
 */
double norm(point f);


/*
 *PROCEDURE: scalar
 * 
 *DESCRIPTION: Calculates scalar product of 2 cartesian Vectors
 * 
 *RETURNS: scalar product(double)
 */
double scalar(point a, point g);

/*
 * Struct: body
 * 
 * OBJECTS: Name, position, velocity, acceleration, mass and radius
*/
struct body{

	point location;
	double mass;
	double radius;
    point velocity;
	std::string name;
	std::vector<point> locations;
};


