//Exemple de n-cossos de https://www.rosettacode.org
#include<stdlib.h>
#include<stdio.h>
#include<math.h>

typedef struct{
	double x,y,z;
}VECTOR;

int bodies,timeSteps;
double *masses,GravConstant;
VECTOR *positions,*velocities,*accelerations;

VECTOR addVECTORs(VECTOR a,VECTOR b){
	VECTOR c = {a.x+b.x,a.y+b.y,a.z+b.z};

	return c;
}

VECTOR scaleVECTOR(double b,VECTOR a){
	VECTOR c = {b*a.x,b*a.y,b*a.z};

	return c;
}

VECTOR subtractVECTORs(VECTOR a,VECTOR b){
	VECTOR c = {a.x-b.x,a.y-b.y,a.z-b.z};

	return c;
}

double mod(VECTOR a){
	return sqrt(a.x*a.x + a.y*a.y + a.z*a.z);
}

void initiateSystem(char* fileName){
	int i;
	FILE* fp = fopen(fileName,"r");

	fscanf(fp,"%lf%d%d",&GravConstant,&bodies,&timeSteps);

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
				accelerations[i] = addVECTORs(accelerations[i],scaleVECTOR(GravConstant*masses[j]/pow(mod(subtractVECTORs(positions[i],positions[j])),3),subtractVECTORs(positions[j],positions[i])));
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

int main(int argC,char* argV[])
{
	int i,j;

	if(argC!=2)
		printf("Usage : %s <file name containing system configuration data>",argV[0]);
	else{
		initiateSystem(argV[1]);
		printf("Body   :     x              y               z           |           vx              vy              vz   ");
		for(i=0;i<timeSteps;i++){
			printf("\nCycle %d\n",i+1);
			simulate();
			for(j=0;j<bodies;j++)
				printf("Body %d : %lf\t%f\t%lf\t|\t%lf\t%lf\t%lf\n",j+1,positions[j].x,positions[j].y,positions[j].z,velocities[j].x,velocities[j].y,velocities[j].z);
		}
	}
	return 0;
}
 
