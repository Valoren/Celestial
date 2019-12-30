/*
Vector and body data structures definition.
Basic vector operations used in numerical integration are also defined.
Miquel Bernat Laporta i Granados
UAB 2019
*/

#pragma once

#include <fstream>
#include <string>
#include <cmath>
#include <vector>

using namespace std;

//Global scope variables used on the integration
double Grav_const;
double timeStep;
double Cycles;
double integration_time;
int num_bodies;

//Vector definition. Required vector operations are also defined.
class Vector
{
  public:

  double x;
  double y;
  double z;
  /*
  //Norm of a vector 
  double norm(Vector f)
  {
    return sqrt(f.x * f.x + f.y * f.y + f.z * f.z);
  }

  //Scalar product of 2 vectors
  double scalar_p(Vector a, Vector g)
  {
    return double(x * a.x + y * a.y + z * a.z);
  }
  */
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
};

const Vector ORIGIN{0.0, 0.0, 0.0};

/*
Body object. Atributes are name, position, velocity, acceleration, mass and radius.
Spherical objects are assumed for simplicity of the simulation.
*/
class body
{
public:

  vector < string > name;
  vector < Vector > vel;
  vector < Vector > pos;
  vector < Vector > accel;
  vector < double >radius;
  vector < double >mass;
  void initiate_system (char *filename);
  void update_pos ();
  void update_vel ();
  void compute_accel ();
};

/*
PROCEDURE: initiate_system
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
.
RETURNS: void
*/
void body::initiate_system (char *filename)
{
  ifstream input_file;
  input_file.open (filename, ios::in);
  if (!input_file.is_open ())
    {
      throw runtime_error ("Error in file opening...");
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
                