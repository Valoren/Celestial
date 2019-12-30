/*
Basic integration library including Euler first and second order and
Runge-Kutta 4th order algorithms. Barnes-Hut will be aviable in the future.
Miquel Bernat Laporta i Granados
UAB 2019
*/

#include <omp.h>
#include <time.h>
#include <math.h>
#include "data_structures.h"

#define NUMBER_OF_ALGORITHMS 4

/*
Algortihm: Euler first order:

Descrption: advances particles one time step

Returns:      0 if all goes well
              1 if error occurs
*/
void tech_euler_1(double timeStep, double total_time){

  int i,j;
  double t = 0; //initial time of simulation
  do{
  //acceleration calculation loop(force)
	for(i=0;i<NUMBER_OF_BODIES;i++){
		accel[i].x = 0;
		accel[i].y = 0;
		accel[i].z = 0;
		for(j=0;j<NUMBER_OF_BODIES;j++){
			if(i!=j){
				accel[i].x = accel[i].x +   addVECTORs(accelerations[i],scaleVECTOR(GravConstant*masses[j]/pow(norm(subtractVECTORs(positions[i],positions[j])),3),subtractVECTORs(positions[j],positions[i])));
			}
		}
	}

  //positions computations
  for(i=0;i<NUMBER_OF_BODIES;i++)
		positions[i] = addVECTORs(positions[i],addVECTORs(velocities[i],scaleVECTOR(0.5,accelerations[i])));

  //velocities computations
  for(i=0;i<NUMBER_OF_BODIES;i++)
		velocities[i] = addVECTORs(velocities[i],accelerations[i]);

  //collisions resolution
  for(i=0;i<NUMBER_OF_BODIES-1;i++)
		for(j=i+1;j<NUMBER_OF_BODIES;j++){
			if(positions[i].x==positions[j].x && positions[i].y==positions[j].y && positions[i].z==positions[j].z){
				VECTOR temp = velocities[i];
				velocities[i] = velocities[j];
				velocities[j] = temp;
			}
		}
  t += suggested_timestep
  }
  while(t < total_time);
}

/*
Algortihm: Euler second order:

Descrption: advances particles one time step

Returns:      0 if all goes well
              1 if error occurs
*/
int tech_euler_2(double suggested_timestep, double *total_time);

int tech_rk4(double suggested_timestep, double *total_time);
