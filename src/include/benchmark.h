/*
 * benchmark.cpp
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

#pragma once

#ifdef _WIN32
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#else
#include <unistd.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <chrono>

class Timer
{
public:
      Timer();
      ~Timer();
     
      void Stop(){
          auto endTimepoint = std::chrono::high_resolution_clock::now();
          auto start = std::chrono::time_point_cast<std::chrono::microseconds>(m_StartTimepoint).time_since_epoch().count();
          auto stop = std::chrono::time_point_cast<std::chrono::microseconds>(endTimepoint).time_since_epoch().count();
          auto duration = stop - start;
          double ms = duration * 0.001;

          std::cout << duration << "us (" << ms << "ms)\n";
}
private:
        std::chrono::time_point<std::chrono::high_resolution_clock> m_StartTimepoint;
};

Timer::~Timer(){

          Stop();
}

Timer::Timer(){

         m_StartTimepoint = std::chrono::high_resolution_clock::now();
}

/*
*PROCEDURE: number_of_cores
*
*DESCRIPTION: Detects system OS and displays number of processing cores aviable on the system
*Main purpose is for multithreading optimization with OMP library.
*
*/
void number_of_cores()
{
  long nprocs = -1;
  long nprocs_max = -1;
#ifdef _WIN32
#ifndef _SC_NPROCESSORS_ONLN
SYSTEM_INFO info;
printf("Detected Windows system\n");
GetSystemInfo(&info);
#define sysconf(a) info.dwNumberOfProcessors
#define _SC_NPROCESSORS_ONLN
#endif
#endif
#ifdef _SC_NPROCESSORS_ONLN
printf("Detected GNU-Linux system");
  nprocs = sysconf(_SC_NPROCESSORS_ONLN);
  if (nprocs < 1)
  {
    fprintf(stderr, "Could not determine number of CPUs online:\n%s\n", 
strerror (errno));
    exit (EXIT_FAILURE);
  }
  nprocs_max = sysconf(_SC_NPROCESSORS_CONF);
  if (nprocs_max < 1)
  {
    fprintf(stderr, "Could not determine number of CPUs configured:\n%s\n", 
strerror (errno));
    exit (EXIT_FAILURE);
  }
  printf ("%ld of %ld processors online\n",nprocs, nprocs_max);
  //exit (EXIT_SUCCESS);
#else
  fprintf(stderr, "Could not determine number of CPUs");
  exit (EXIT_FAILURE);
#endif
}