/*********************************************************************
Project:  tspp
Filename: tspp.c
Author:   Joe Heafner.
Purpose:  Example driver for Stat2Elem.
*********************************************************************/

/* Header Files *****************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "astrolib.h"

void GetUserData(double *state, double *mu, double *jdate);
void UseThisData(double *state, double *mu, double *jdate);

/*********************************************************************
Name:    main
Purpose: Main routine for tspp.
Inputs:  None.
Outputs: None.
Returns: 0 if execution successful.
Status:  Finished.
Errors:  None known.
*********************************************************************/
int main(void) {

  char line[1024];
  int Which;
  double element[6], state[6];
  double mu, jdate;

  printf("EXAMPLE DRIVER FOR Stat2Elem\n\n");

  do {
    do {
      fprintf(stdout, "[1] Input data from the keyboard\n");
      fprintf(stdout, "[2] Use default data\n");
      fprintf(stdout, "Your choice: ");
      fflush(stdout);
      fgets(line, sizeof(line), stdin);
    } while (sscanf(line, "%d", &Which) != 1);
  } while ((Which < 1) || (Which > 2));

  printf("\n");

  switch (Which) {
  case 1:
    GetUserData(state, &mu, &jdate);
    break;
  case 2:
    UseThisData(state, &mu, &jdate);
    break;
  }

  Stat2Elem(state, mu, jdate, element);

  printf("INITIAL STATE VECTOR\n");
  printf("%+19.15f %+19.15f %+19.15f AU\n", state[0], state[1], state[2]);
  printf("%+19.15f %+19.15f %+19.15f AU/DAY\n\n",
    state[3], state[4], state[5]);

  printf("ELEMENT SET\n");
  printf("Perihelion distance (AU) %11.8f\n", element[0]);
  printf("Eccentricity             %12.9f\n", element[1]);
  printf("Inclination      (deg) %+13.8f\n", element[2] * R2D);
  printf("Long. of node    (deg) %+13.8f\n", element[3] * R2D);
  printf("Arg. of peri.    (deg) %+13.8f\n", element[4] * R2D);
  printf("T (TDT)             %16.8f\n", element[5]);

  printf("\nOK\n");
  return (0);
}

/*********************************************************************
Name:    GetUserData
Purpose: Routine to allow user to input desired data from the keyboard.
Inputs:  None.
Outputs: state[].
         mu.
         jdate.
Returns: Nothing
Status:  Finished.
Errors:  None known.
*********************************************************************/
void GetUserData(double *state, double *mu, double *jdate) {

  char line[1024];
  double mass;

  printf("\n");

  do {
    fprintf(stdout, "Enter position vector X component: ");
    fflush(stdout);
    fgets(line, sizeof(line), stdin);
  } while (sscanf(line, "%lf", &state[0]) != 1);

  do {
    fprintf(stdout, "Enter position vector Y component: ");
    fflush(stdout);
    fgets(line, sizeof(line), stdin);
  } while (sscanf(line, "%lf", &state[1]) != 1);

  do {
    fprintf(stdout, "Enter position vector Z component: ");
    fflush(stdout);
    fgets(line, sizeof(line), stdin);
  } while (sscanf(line, "%lf", &state[2]) != 1);

  do {
    fprintf(stdout, "Enter velocity vector X component: ");
    fflush(stdout);
    fgets(line, sizeof(line), stdin);
  } while (sscanf(line, "%lf", &state[3]) != 1);

  do {
    fprintf(stdout, "Enter velocity vector Y component: ");
    fflush(stdout);
    fgets(line, sizeof(line), stdin);
  } while (sscanf(line, "%lf", &state[4]) != 1);

  do {
    fprintf(stdout, "Enter velocity vector Z component: ");
    fflush(stdout);
    fgets(line, sizeof(line), stdin);
  } while (sscanf(line, "%lf", &state[5]) != 1);

  do {
    fprintf(stdout, "Enter Julian date: ");
    fflush(stdout);
    fgets(line, sizeof(line), stdin);
  } while (sscanf(line, "%lf", jdate) != 1);

  do {
    fprintf(stdout, "Enter mass in solar masses: ");
    fflush(stdout);
    fgets(line, sizeof(line), stdin);
  } while (sscanf(line, "%lf", &mass) != 1);

  *mu = (double) (GAUSSK * GAUSSK * (1.0 + mass));
}

/*********************************************************************
Name:    UseThisData
Purpose: Subprogram to supply default input data. Data is taken from
         the paper by Mansfield (AIAA Paper No. 86-2269-CP) and is
         for Comet Halley.
Inputs:  None.
Outputs: state[].
         mu.
         jdate.
Returns: Nothing
Status:  Finished.
Errors:  None known.
*********************************************************************/
void UseThisData(double *state, double *mu, double *jdate) {

  double mass;

  /* Components referred to equinox/ecliptic of B1950 */

  /* Position vector components */
  state[0] = 0.7661010445554045;
  state[1] = 0.1412849664351819;
  state[2] = 0.1845468752429582;

  /* Velocity vector components */
  state[3] = -0.0101488586182795;
  state[4] = -2.485310205695914e-02;
  state[5] = 1.440200276392141e-03;

  /* Mass */
  mass = 0.0;
  *mu = (double) (GAUSSK * GAUSSK * (1.0 + mass));

  /* Julian day number (TDB) */
  *jdate = (double) 2446445.5;
}

/* End Of File - tspp.c *********************************************/
