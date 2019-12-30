/*
 * main.cpp
 * 
 * Copyright 2019 Miquel Bernat Laporta i Granados <mlaportaigranados@gmail.com>
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

/*
Matplotlib link instructions:
Compile with folowing flags:
-I/usr/include/python2.7 -lpython2.7
*/

//#include "matplotlibcpp.h"
#include <iostream>
#include <stdlib.h>
#include <math.h>
//#include "cores.h"
#include "structures.h"
#include <omp.h>
//#include "integration.h"
#include <string>
//#include "planet_data.h"
#include <fstream>
#include <ctime>
#include <stdio.h>

//using namespace solar_system;


#define INVALID_ALGORTIHM -1

/*
PROCEDURE: randomize_galaxy
DESCRIPTION: Creates random galaxy-like system with central mass(black hole)
at origin. Positions and velocities are all randomized.

RETURNS : Data file in standard format
*
void randomize_galaxy(char *data_filename, int N){

        FILE *fp;
 nap/firefox/258'
 double upper = 8.56828e36;
        double lower = 7.99176e36;
        double black_hole_mass =  (rand() % (upper - lower + 1)); 
        fprintf("%d", &N);
        fprintf("%lf", &black_hole_mass);
        fprintf(fp, "0 0 0\n");
		    fprintf(fp, "0 0 0\n");
        for(int i=0; i<N; i++){
          fprintf(fp,"%li\n",rand() % 10000000000000);
          fprintf(fp,"%d %d %d \n",rand() % 10,rand() % 10,rand() % 10);
          fprintf(fp,"%d %d %d \n",rand() % 10,rand() % 10,rand() % 10);
}
}
*/

/*
List of integration algorithms used
Function for interpreting algorithm input by string comparation.
Returns invalid algorithm if not in list
*/


/* 
void compareFunction(string s1, string s2) 
{ 
    int g = s1.compare(s2); 
  
    if (x != 0){
        cout << "Error algorithm"<< s1 <<"not defined..." << endl;
        exit(EXIT_FAILURE);
    }
    else cout << "Initializing simulation with" << s1 << "algorithm." << endl;
   } 
*/

/*
Initiating the simulation
*/
int main (int argc, char *argV[])
{
  //setup_m0();
  if ((argc) = !2)
    printf ("Usage: %s <input file name>\n", argV[0]);
  else
    {
      printf ("\n \n");
      printf
	("Welcome to Celestial gravitation simulation: N-Body simulations\nand orbits calculation\nVersion 1.1\n");
      printf
	("This program is free software: you can redistribute it and/or modify it.\n");
      printf
	("Written by Miquel Bernat Laporta, Mathematics UAB \n2017-2019");
      printf ("Initiating system with given inital data\n");
      initiate_system (argV[1]);
      //compareFunction(algorithm_name, str1);  
      


}
  return 0;
}
