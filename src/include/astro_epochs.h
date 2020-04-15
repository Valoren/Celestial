/*
 * astro_epochs.h
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

#pragma once

#include <cmath>

std::string week_day[] =
{
"monday", "tuesday", "wednesday", "thursday", "friday", "saturday", "sunday"};

/*
*PROCEDURE: calendar_to_julian
*
*DESCRIPTION: Converts given date to julian calendar date
*
*/
double
calendar_to_julian (double year, double month, double day)
{

  double
    a,
    b;
  if (month <= 2)
    {
      year--;
      month += 12;
    }

  a = floor (year / 100);
  b = year + 0.01 * month + 0.0001 * day >= 1582.1015 ? 2 - a + floor (a / 4)	/* correcció per la reforma gregoriana */
    : 0;

  return floor (365.25 * (year + 4716)) + floor (30.6001 * (month + 1)) +
    day + b - 1524.5;
}

/*
*
*PROCEDURE: julian_to_calendar
*
*DESCRIPTION: Converts given julian date to calendar
*
*/
double
julian_to_calendar (double jd, double *yy, double *mm, double *dd)
{
  double z, f, b, c, d, e;
  jd += 0.5;
  z = floor (jd);
  f = jd - z;
  b = z;
  if (z >= 2299161)
    {
      double
	alf = floor ((z - 1867126.25) / (36524.5));
      b += 1 + alf - floor (alf / 4);
    }
  b += 1524;
  c = floor ((b - 122.1) / 365.25);
  d = floor (365.25 * c);
  e = floor ((b - d) / 30.6001);
  *dd = b - d - floor (30.6001 * e) + f;
  *mm = e < 14 ? e - 1 : e - 13;
  *yy = 4716;
  if (*mm <= 2)
    (*yy)--;
  *yy = c - (*yy);
}
