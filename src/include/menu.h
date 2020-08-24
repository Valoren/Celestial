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
#include "integration.h"

typedef void (*Menu_Processing_Function_Pointer)();

struct Menu_Option
{
  char choice;
  char const * p_selection_text;
  Menu_Processing_Function_Pointer p_processing_function;
};

/*
*GENERIC PROCEDURES: Process_Selection
*
*DESCRIPTION: Initializes integration with desired algorithm when prompted
 * on menu
*
*/
void Process_Selection_One();
void Process_Selection_Two();
void Process_Selection_Three();
void Process_Selection_Four();

static const Menu_Option main_menu[] =
{
  {'1', "Euler first order for 2 body systems", Process_Selection_One},
  {'2', "F and G series for 2 body systems", Process_Selection_Two},
  {'3', "Runge Kutta-Fehlberg 5th order", Process_Selection_Two},
  {'4',"Runge-Kutta 4th order, Process_Selection_Two"}
};

/*
*PROCEDURE: spawn_title
* 
*DESCRIPTION: Writes start information on program execution
*
*/
void spawn_title();

/*
*PROCEDURE: spawn_menu
* 
*DESCRIPTION: Spawns menu and prompts available integration algorithms for the user
*to choose from
*
*/
void spawn_menu();

/*
*PROCEDURE: error_message
*
*DESCRIPTION: Prints error message and terminates program execution
*
*/
void error_message(const std::string& msg);