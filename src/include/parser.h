/*
 * parser.h
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

#include <cstdio>
#include <cstdlib>
#include <iostream>

/*
 *PROCEDURE: check_comment
 *
 *DESCRIPTION: Checks for both types of comments,
 * then passes on control to below comments
 *
 *RETURNS: -
 *
 */
void check_comment (char);

/*
 *PROCEDURE: block_comment
 *
 *DESCRIPTION: Handles block or single line comments
 *
 *RETURNS: -
 *
 */
void block_comment();

/*
 *PROCEDURE: single_comment
 *
 *DESCRIPTION: Handles single line comments
 *
 *RETURNS: -
 *
 */
void single_comment();

/*
 *PROCEDURE: parse_file
 *
 *DESCRIPTION: Parses file and removes comments
 *
 *RETURNS: -
 *
 */
void parse_file();

/*
 *PROCEDURE: initiateSystem
 *
 *DESCRIPTION: Parses file and initializes system with received bodies (in C)
 *
 *RETURNS: -
 *
 */
void initiateSystem(char * filename);

/*
 *PROCEDURE: parse_data
 *
 *DESCRIPTION: Parses file and initializes system with received bodies (in C++)
 *
 *RETURNS: -
 *
 */
void parse_data(const std::string& filename);

