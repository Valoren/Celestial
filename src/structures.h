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

//Global scope variables used on the integration
double Grav_const;
double timeStep;
double Cycles;
double integration_time;
int num_bodies;

/*
 * 	CLASS: Vector
 * 
 *  Members: Coordinates, basic vector operations and ">>" operator overload
 */
class Vector{

public:

  double x;
  double y;
  double z;

  //Product of a scalar and a vector
  Vector operator* (double k){
    return Vector
    {
    k *x, k * y, k *  z};
  }
  //Vector sum
  Vector operator+ (Vector b){
    return Vector
    {
    x + b.x, y + b.y, z + b.z};
  }
  //Vector rest
  Vector operator- (Vector c){
    return Vector
    {
    x - c.x, y - c.y, z - c.z};
  }
  //+= Vector operator
  Vector operator+= (const Vector & d){
    x += d.x;
    y += d.y;
    z += d.z;
    return Vector
    {
    x, y, z};
  }
  //Cross product of 2 Vectors 
  Vector operator^ (Vector e){
    return Vector{
    y * e.z - z * e.y, z * e.x - x * e.z, x * e.y - y * e.x };
  }
  //Operator >> overloading for input file method
  friend std::istream& operator>>(std::istream&, Vector&);
};
std::istream& operator>>(std::istream& in, Vector& v) {
    return in >> v.x >> v.y >> v.z;
}


/*
 *PROCEDURE: norm
 * 
 *DESCRIPTION: Calculates norm of a given cartesian Vector
 * 
 * RETURNS: norm of Vector(double)
 */
double norm(Vector f){
    return sqrt(f.x * f.x + f.y * f.y + f.z * f.z);
}


/*
 *PROCEDURE: scalar
 * 
 *DESCRIPTION: Calculates scalar product of 2 cartesian Vectors
 * 
 * RETURNS: scalar product(double)
 */
double scalar(Vector a, Vector g){
  return a.x * g.x + a.y * g.y + a.z * g.z;
}

const Vector ORIGIN{0.0, 0.0, 0.0};

/*
 * CLASS: body
 * 
 * OBJECTS: Name, position, velocity, acceleration, mass and radius
*/
class body{
public:

  std::vector < std::string > name;
  std::vector < Vector > vel;
  std::vector < Vector > pos;
  std::vector < Vector > accel;
  std::vector < double >radius;
  std::vector < double >mass;
  void initiate_system (std::string& filename);
  void update_pos ();
  void update_vel ();
  void compute_accel ();
};

/*
PROCEDURE: initiate_system
* 
DESCRIPTION: Initiates system with initial data input from filename.txt
Input format:

Grav_constant number_of_bodies integration_timestep
name_of_body_1
mass_1
pos.x_1 pos.y_1 pos.z_1
vel.x_1 vel.y_1 vel.z_1
name_of_body_2
.
.
.
*
*  
RETURNS: void
*/
void body::initiate_system (std::string& filename){
  ifstream input_file;
  input_file.open (filename, ios::in);
  if (!input_file.is_open ())
    {
      throw std::runtime_error ("Error in file opening...");
    }
  //Read header data and compute number of cycles to simulate
  input_file >> Grav_const >> num_bodies >> timeStep;
  Cycles = floor (integration_time / timeStep);

  //Initiate structures and set all parameters to 0
  mass.resize (num_bodies);
  pos.resize (num_bodies);
  fill (pos.begin (), pos.end (), ORIGIN);
  vel.resize (num_bodies);
  fill (vel.begin (), pos.end (), ORIGIN);
  accel.resize (num_bodies);
  fill (accel.begin (), accel.end (), ORIGIN);

  for (int i = 0; i < num_bodies; ++i)
    {
      input_file >> mass[i] >> pos[i] >> vel[i];
    }
}
                
