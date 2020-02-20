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
#include <stdlib.h>
#include <math.h>
//#include "cores.h"
//#include "structures.h"
#include <omp.h>
//#include "integration.h"
#include <string>
//#include "planet_data.h"
#include <fstream>
#include <ctime>
#include <stdio.h>
#include "menu.h"
#include "parser.h" 
#include "exit_codes.h"
#include "benchmark.h"
//using namespace solar_system;

#define INVALID_ALGORTIHM -1



/*
Initiating the simulation
*/
int main(){
  /*
  spawn_menu();
  open_files(read,write);//calling the open_files() function
  close_files(read,write);//calling the close_files() function
  return EXIT_SUCCESS;
  */
  int value = 0;
  {
        Timer timer1;
        for (int i=0; i<1000000; i++)
            value +=2;

    }
    std::cout << value << std::endl;
}
