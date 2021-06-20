/*
 * menu.cpp
 *
 * Copyright 2019 Miquel Bernat Laporta i Granados
 * <mlaportaigranados@gmail.com>
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

#include "menu.h"
#include "integration.h"

/*
 *PROCEDURE: spawn_title
 *
 *DESCRIPTION: prints basic welcome
 *
 *RETURNS: -
 *
 */
void spawn_title(){
    std::cout << "-----------------------------------------\n" << std::endl;
    std::cout << "Celestial\n" << std::endl;
    std::cout << "A stellar dynamics toolkit\n" << std::endl;
    std::cout << "Written by Miquel Bernat Laporta i Granados\n Mathematics UAB\n" << std::endl;
    std::cout << "2017-2020\n" << std::endl;
    std::cout << "-----------------------------------------\n" << std::endl;
}

/*
 *PROCEDURE: read_options
 *
 *DESCRIPTION:  reads the command line options, and implements them.
 *
 *RETURNS: bool
 *
 */
bool read_options(int argc, char *argv[], real & dt_param, real & dt_dia,
                  real & dt_out, real & dt_tot, bool & i_flag, bool & x_flag)
{
    int c;
    while ((c = getopt(argc, argv, "hd:e:o:t:ix")) != -1)
        switch(c){
            case 'h': cerr << "usage: " << argv[0]
                           << " [-h (for help)]"
                           << " [-d step_size_control_parameter]\n"
                           << "         [-e diagnostics_interval]"
                           << " [-o output_interval]\n"
                           << "         [-t total_duration]"
                           << " [-i (start output at t = 0)]\n"
                           << "         [-x (extra debugging diagnostics)]"
                           << endl;
                      return false;         // execution should stop after help
            case 'd': dt_param = atof(optarg);
                      break;
            case 'e': dt_dia = atof(optarg);
                      break;
            case 'i': i_flag = true;
                      break;
            case 'o': dt_out = atof(optarg);
                      break;
            case 't': dt_tot = atof(optarg);
                      break;
            case 'x': x_flag = true;
                      break;
            case '?': cerr << "usage: " << argv[0]
                           << " [-h (for help)]"
                           << " [-d step_size_control_parameter]\n"
                           << "         [-e diagnostics_interval]"
                           << " [-o output_interval]\n"
                           << "         [-t total_duration]"
                           << " [-i (start output at t = 0)]\n"
                           << "         [-x (extra debugging diagnostics)]"
                           << endl;
                      return false;        // execution should stop after error
            }

    return true;                         // ready to continue program execution
}




void spawn_menu(){
    spawn_title();
    for (const auto & i : main_menu)
    {
        std::cout << i.p_selection_text << "\n";
    }
    std::cout << "Enter selection, 0 to quit: ";
    char choice;
    std::cin >> choice;
    for (const auto & i : main_menu)
    {
        if (choice == i.choice)
        {
            Menu_Processing_Function_Pointer p_function = i.p_processing_function;
            (p_function)();
            break;
        }
    }
}

void error_message(const std::string& msg)
{
std::cerr<<msg<<std::endl;
exit(EXIT_FAILURE);
}

void print_cmd_options(){
    std::cout << "Usage: Celestial [input_file] [algorithm] [steps] [output_file]\n" << std::endl;
    std::cout << "---------------------------------------------------------------\n" << std::endl;
    std::cout << "Implemented algorithms:\n" << std::endl;
    std::cout << "-RK4 - Runge-Kutta 4th order\n" << std::endl;
    std::cout << "-Euler - Standard Euler integration\n" << std::endl;
    std::cout << "-f_and_g - F and G series expansion(for 2 body systems only!)" << std::endl;
    std::cout << "-taylor - Taylor series expansion(for 2 body systems only!)" << std::endl;
    std::cout << "---------------------------------------------------------------\n" << std::endl;
    std::cout << "Extra flags:" << std::endl;
    std::cout << "-test - Uses default solar system testing system" << std::endl;
    
    exit(EXIT_FAILURE);
}

void Process_Selection_One(){
  std::cout << "Caca" << std::endl;
}
void Process_Selection_Two(){
    std::cout << "Caca" << std::endl;

}
void Process_Selection_Three(){
    std::cout << "Caca" << std::endl;

}
void Process_Selection_Four(){
    std::cout << "Caca" << std::endl;

}