#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "vectors.h"
//#include "planet_data.h"

#define NMAX  1000 //Max recommended number of bodies

int bodies;
double *masses;
double GravConstant;
VECTOR *positions,*velocities,*accelerations;

void initiateSystem(char* fileName){
	int i;
	FILE* fp = fopen(fileName,"r");

	fscanf(fp,"%lf%d",&GravConstant,&bodies);
  if(bodies <= NMAX){

	masses = (double*)malloc(bodies*sizeof(double));
	positions = (VECTOR*)malloc(bodies*sizeof(VECTOR));
	velocities = (VECTOR*)malloc(bodies*sizeof(VECTOR));
	accelerations = (VECTOR*)malloc(bodies*sizeof(VECTOR));

	for(i=0;i<bodies;i++){
		fscanf(fp,"%lf",&masses[i]);
		fscanf(fp,"%lf%lf%lf",&positions[i].x,&positions[i].y,&positions[i].z);
		fscanf(fp,"%lf%lf%lf",&velocities[i].x,&velocities[i].y,&velocities[i].z);
	}

	fclose(fp);
}
else{
  printf("Exceded recommended max number of bodies\n");
  exit(0);
}
}

void resolveCollisions(){
	int i,j;

	for(i=0;i<bodies-1;i++)
		for(j=i+1;j<bodies;j++){
			if(positions[i].x==positions[j].x && positions[i].y==positions[j].y && positions[i].z==positions[j].z){
				VECTOR temp = velocities[i];
				velocities[i] = velocities[j];
				velocities[j] = temp;
			}
		}
}

void computeAccelerations(){
	int i,j;

	for(i=0;i<bodies;i++){
		accelerations[i].x = 0;
		accelerations[i].y = 0;
		accelerations[i].z = 0;
		for(j=0;j<bodies;j++){
			if(i!=j){
				accelerations[i] = addVECTORs(accelerations[i],scaleVECTOR(GravConstant*masses[j]/pow(norm(subtractVECTORs(positions[i],positions[j])),3),subtractVECTORs(positions[j],positions[i])));
			}
		}
	}
}

void computeVelocities(){
	int i;

	for(i=0;i<bodies;i++)
		velocities[i] = addVECTORs(velocities[i],accelerations[i]);
}

void computePositions(){
	int i;

	for(i=0;i<bodies;i++)
		positions[i] = addVECTORs(positions[i],addVECTORs(velocities[i],scaleVECTOR(0.5,accelerations[i])));
}

void simulate(){
	computeAccelerations();
	computePositions();
	computeVelocities();
	resolveCollisions();
}


int main(int argc, char* argv[]) {

double t=0;
double dt=0.1;
FILE *ft;
int j;
ft = fopen(argv[2],"w");
if (ft == NULL){
	printf ("File error...");
  exit(0);
}
if(argc!=3){
  printf("Usage : %s <file name containing system configuration data> <file name for result dump>\n",argv[0]);
  return -1;
}
else{
  initiateSystem(argv[1]);
  fprintf(ft,"Body   :     x      y      z   |      vx      vy      vz   ");
  while(t<10){
    fprintf(ft,"\nCycle %.3G\n",t/dt);
	  simulate();
    for(j=0;j<bodies;j++){
      fprintf(ft,"Body %d : %lf\t%f\t%lf\t|\t%lf\t%lf\t%lf\n",j+1,positions[j].x,positions[j].y,positions[j].z,velocities[j].x,velocities[j].y,velocities[j].z);
}
    t += dt;
}
}
fclose(ft);
return 0;
}
