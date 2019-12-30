/*
3D vectors main library. Components denoted by x,y,z
*/
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

typedef struct vector {
  long double x,y,z;
} VECTOR;

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

double norm(VECTOR a){
	return sqrt(a.x*a.x + a.y*a.y + a.z*a.z);
}

double dotproduct(VECTOR a,VECTOR b){
  return (a.x*b.x+a.y*b.y+a.z*b.z);
}

VECTOR crossp(VECTOR a,VECTOR b){
  VECTOR c = {a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x};

  return c;
}

typedef struct body{

        int index[10];
        char name[255];
        double mass;
        VECTOR positions;
        VECTOR velocities;
        VECTOR accel;
}body;
