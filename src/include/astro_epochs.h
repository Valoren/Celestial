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

#include"astro_constants.h"

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
#include <sys/types.h>
#include <netinet/in.h>
#include <stdarg.h>
#include <sys/param.h>
/* #include <386/endian.h> */

#define MAX_NAME_SIZE       255
#define MAX_EXTENSION_SIZE  5
#define BIG_ENDIAN_TEST     (htonl(1) == 1)

#ifndef TRUE
#define TRUE   (1)
#define FALSE  (0)
#endif

/* Function prototypes */
void make_little_endian(char *ptr, int len);
void convert_little_endian(char *ptr, int len);
void reverse_bytes(char *ptr, int len);
void ucase(char str[]);
void left(char str[], int n, char dest[]);
void right(char str[], int n, char dest[]);
void Trim(char str[]);
void RightTrim(char str[]);
void LeftTrim(char str[]);
int  fexist(char *filename);
int  LogOpen(char *filename);
void LogClose(void);
void LogMsg(FILE *fptr, const char *format, ...);

#endif /* _SUPPORT_ */

#ifndef TRUE
#define TRUE  (1)
#define FALSE (0)
#endif

#ifndef SEEK_SET
#define SEEK_SET  (0)
 #define SEEK_CUR  (1)
 #define SEEK_END  (2)
#endif

/* double precision matrix data type */
typedef double** DMatrix;

/* function prototypes */
void Aberrate(double *p1, double *EBdot, double *p2);
double amodulo(double a, double b);
void Cal2JED(int m, int d, int y, double utc, int s, double tai_utc,
             double ut1_utc, int w, double *jed);
void constants(char *FileName);
void Conway(double uele[], double mu, double jed, double *r_cos_nu,
            double *r_sin_nu);
DMatrix createDMatrix(int n);
double deg(double x);
double dms(double x);
double DRound(double x, int n);
FILE *ephopn(char *FileName);
void Epoch2JED(char *epoch, double *jed);
void Eq2Ecl(double *a, int s, double eps, double *b);
void Eq2Hor(double *a, int s, double *b);
void errprt(int group, char *message);
double fix(double x);
void FmtDms(double x, int n, int m, char *s);
void free_matrix(double **m, int nr, int nc);
void freeDMatrix(DMatrix m, int n);
void FunArgIAU(double jed, double *funarg);
void GeocenObs(double jed, double *obsr_geo);
void GetGST(double jed, int s, double *gst);
void GetInvQMatrix(DMatrix QMatrix, DMatrix InvQMatrix);
int  GetNumde(void);
void GetPrecessParams(double jed1, double jed2, double *zeta,
                      double *z, double *theta, double *zetadot,
                      double *zdot, double *thetadot);
void GetQMatrix(double phi, double phidot, int axis, int s,
                DMatrix QMatrix);
void GetRPNmat(double jed1, double jed2, int rpn, int d,
               DMatrix m3, DMatrix m6);
void GetStateVector(double JD, int TARG, int CENT, int recpol,
                    double StateVector[15][15][2][6]);
void HelEphemeris(double *uelement, double mu, double jed, double *posvel);
void interp(int buff, double *t, int ncf, int ncm, int na, int fl,
            double pv[3][2]);
void Interp1(double *x, double *y, int L, int m, double *arg, double *v);
double Interpol(double *x, double *y, int i, double arg);
void JED2Cal(double jed, int *yr, int *mo, int *dy, double *ti);
void JED2Epoch(double jed, char *s, char *epoch);
double **matrix(int nr, int nc);
void MatXMat(DMatrix a, DMatrix b, DMatrix c, int n);
void MatXVec(DMatrix a, double *b, double *c, int l, int m);
void MRotate(double *vin, int axis, double phi, double *vout);
void GetDpsiDeps(double jed, double *dpsi, double *deps,
                 double *dpsidot, double *depsdot);
void Obliquity(double jed1, double jed2, int m, double *obl,
               double *obldot);
void pleph(double jd, int targ, int cent, double *rrd, int *inside);
void Pol2Rec(double *a, double *b);
void PrecessElements(char *eqnx1, double *element, char *eqnx2);
void QRotate(double *vin, int axis, double phi, double phidot,
             int s, double *vout);
void RayBend(double *earth_hel, double *body_geo, double *body_hel,
             double *p1);
void Rec2Pol(double *a, double *b);
void Reduce (double jed, int body, int place, double StarData[],
             double p3[]);
void Refract(double ra1, double dec1, double lha, double temp,
             double press, double *ra2, double *dec2);
void RotMat(int axis, double phi, DMatrix r, DMatrix dr);
void RST(double jed, double *ra, double *dec, double z0, double deltat,
         char *ris, char *trn, char *set);
void split(double tt, double *ipart, double *fpart);
void SplitStateVector(double *pv, double *p, double *v);
void state(double *jed, int LList[], double pv[6][13], double *nut);
void Stat2Elem(double *posvel, double mu, double jed, double *uelement);
double StumpffN(double x, int Norder);
void Transpose(DMatrix a, DMatrix b, int n);
void Uvector(double *a, double *unita);
void Vcross(double *a, double *b, double *acrossb);
void Vdot(int n, double *a, double *b, double *adotb);
double Vecmag(double *a);


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
