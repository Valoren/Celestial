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