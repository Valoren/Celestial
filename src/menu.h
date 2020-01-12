/*
 * menu.h
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
#include <iostream>
//#include "integration.h"

typedef void (*Menu_Processing_Function_Pointer)(void);

struct Menu_Option
{
  char choice;
  char const * p_selection_text;
  Menu_Processing_Function_Pointer p_processing_function;
};

//Dummy functions for testing menu options

void Process_Selection_One(){
    std::cout << "Hello World" << std::endl;
};
void Process_Selection_Two(){
    std::cout << "Hola" << std::endl;
};

static const Menu_Option main_menu[] =
{
  {'1', "Algorithm_1", Process_Selection_One},
  {'2', "Algorithm_2", Process_Selection_Two},
  {'3', "Algorithm_3",Process_Selection_Two}
};

static const size_t quantity_selections =
    sizeof(main_menu) / sizeof(main_menu[0]);

void spawn_title(){
    std::cout << "-----------------------------------------\n" << std::endl;
    std::cout << "Celestial\n" << std::endl;
    std::cout << "A stellar dynamics toolkit\n" << std::endl;
    std::cout << "Written by Miquel Bernat Laporta i Granados\n Mathematics UAB\n" << std::endl;
    std::cout << "2017-2020\n" << std::endl;
    std::cout << "-----------------------------------------\n" << std::endl;
}

void spawn_menu(){
  spawn_title();
  for (size_t i = 0; i < quantity_selections; ++i)
  {
    std::cout << main_menu[i].p_selection_text << "\n";
  }
  std::cout << "Enter selection, 0 to quit: ";
  char choice;
  std::cin >> choice;
  for (size_t i = 0; i < quantity_selections; ++i)
  {
     if (choice == main_menu[i].choice)
     {
       Menu_Processing_Function_Pointer p_function = main_menu[i].p_processing_function;
       (p_function)();
       break;
     }
  }
  }
