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
