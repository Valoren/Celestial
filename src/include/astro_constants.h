/*
 * integration.h
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

#define _USE_MATH_DEFINES

#include <cmath>

/* Gaussian gravitational constant */
#define GAUSSK (0.01720209895)

/* values related to PI */
#ifndef PI
#define PI (139755218526789.0 / 44485467702853.0)
#endif

#define TWOPI  (2.0 * PI)
#define PIDIV2 (0.5 * PI)

/* degrees to radians and radians to degrees */
#define D2R (PI / 180.0)
#define R2D (180.0 / PI)

/* hours to radians and radians to hours */
#define H2R (PI / 12.0)
#define R2H (12.0 / PI)

/* arcseconds to radians and radians to arcseconds */
#define A2R (PI / 648000.0)
#define R2A (648000.0 / PI)

/* AU to km and km to AU */
#define AU2KM (149597870.66)
#define KM2AU (1.0 / 149597870.66)

/* speed of light in km/s and AU/day */
#define CKMS (299792.458)
#define CAUD (CKMS * 86400.0 * KM2AU)

/* Earth's angular velocity in rad/s */
#define EarAngVelRD (6.30038748685432)

/* Earth's equatorial radius in km and AU */
#define EarthRadKM (6378.137)
#define EarthRadAU (EarthRadKM * KM2AU)

/* Earth's flattening factor */
#define EFlat (1.0 / 298.257)

/* JED of standard epoch */
#define J2000 (2451545.0)

/* days in a Julian century */
#define JulCty (36525.0)

/* Earth's gravitational constant */
#define MUC (2.0 * GAUSSK * GAUSSK / (CAUD * CAUD))

/* global state vector variable */
extern double StateVector[15][15][2][6];

extern double obsr_lon;
extern double obsr_lat;
extern double obsr_ele;

