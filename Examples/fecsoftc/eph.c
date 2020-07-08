/*********************************************************************
Project:  eph
Filename: eph.c
Author:   Joe Heafner.
Purpose:  Example driver for HelEphemeris.
*********************************************************************/

/* Header Files *****************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "astrolib.h"

void GetUserData(double *element, double *mass, double *StartingJD,
  double *stepsize, int *numsteps);
void UseThisData (double *element, double *mass, double *StartingJD,
  double *stepsize, int *numsteps);

/*********************************************************************
Name:    main
Purpose: Main routine for eph.
Inputs:  None.
Outputs: None.
Returns: 0 if execution successful.
Status:  Finished.
Errors:  None known.
*********************************************************************/
int main(void) {

  char line[1024];
  int Which, numsteps, i;
  double element[6], mass, StartingJD, stepsize;
  double mu, jdate, posvel[6];

  printf("EXAMPLE DRIVER FOR HelEphemeris\n\n");

  do {
    do {
      fprintf(stdout, "[1] Input data from keyboard\n");
      fprintf(stdout, "[2] Use default Comet Halley data\n");
      fprintf(stdout, "Your choice: ");
      fflush(stdout);
      fgets(line, sizeof(line), stdin);
    } while (sscanf(line, "%d", &Which) != 1);
  } while ((Which < 1) || (Which > 2));

  printf("\n");

  switch (Which) {
  case 1:
    GetUserData(element, &mass, &StartingJD, &stepsize, &numsteps);
    break;
  case 2:
    UseThisData(element, &mass, &StartingJD, &stepsize, &numsteps);
    break;
  }

  /* Compute heliocentric gravitational constant for the body */
  mu = GAUSSK * GAUSSK * (1.0 + mass);

  i = 0;
  jdate = StartingJD;
  do {
    HelEphemeris(element, mu, jdate, posvel);
    printf("date %f\n", jdate);
    printf("position %+17.13f %+17.13f %+17.13f AU\n",
      posvel[0], posvel[1], posvel[2]);
    printf("velocity %+17.13f %+17.13f %+17.13f AU/DAY\n",
      posvel[3], posvel[4], posvel[5]);
    i++;
    jdate += stepsize;
  } while (i < numsteps);

  printf("\nOK\n");
  return (0);
}

/*********************************************************************
Name:    GetUserData
Purpose: Routine to allow user to input desired data
         from the keyboard.
Inputs:  None.
Outputs: element.
         mass.
         StartingJD.
         stepsize.
         numsteps.
Returns: Nothing
Status:  Finished.
Errors:  None known.
*********************************************************************/
void GetUserData(double *element, double *mass, double *StartingJD,
  double *stepsize, int *numsteps) {

  char line[1024];

  printf("\n");
  do {
    fprintf(stdout, "Perihelion distance in AU: ");
    fflush(stdout);
    fgets(line, sizeof(line), stdin);
  } while (sscanf(line, "%lf", &element[0]) != 1);

  do {
    fprintf(stdout, "Eccentricity: ");
    fflush(stdout);
    fgets(line, sizeof(line), stdin);
  } while (sscanf(line, "%lf", &element[1]) != 1);

  do {
    fprintf(stdout, "Inclination in degrees: ");
    fflush(stdout);
    fgets(line, sizeof(line), stdin);
  } while (sscanf(line, "%lf", &element[2]) != 1);

  do {
    fprintf(stdout, "Long. of asc. node in degrees: ");
    fflush(stdout);
    fgets(line, sizeof(line), stdin);
  } while (sscanf(line, "%lf", &element[3]) != 1);

  do {
    fprintf(stdout, "Arg. of peri. in degrees: ");
    fflush(stdout);
    fgets(line, sizeof(line), stdin);
  } while (sscanf(line, "%lf", &element[4]) != 1);

  do {
    fprintf(stdout, "Julian date of peri. passage: ");
    fflush(stdout);
    fgets(line, sizeof(line), stdin);
  } while (sscanf(line, "%lf", &element[5]) != 1);

  element[2] = element[2] * D2R;
  element[3] = element[3] * D2R;
  element[4] = element[4] * D2R;

  do {
    fprintf(stdout, "Body's mass in solar units: ");
    fflush(stdout);
    fgets(line, sizeof(line), stdin);
  } while (sscanf(line, "%lf", mass) != 1);

  do {
    fprintf(stdout, "Starting Julian date: ");
    fflush(stdout);
    fgets(line, sizeof(line), stdin);
  } while (sscanf(line, "%lf", StartingJD) != 1);

  do {
    do {
      fprintf(stdout, "Number of steps: ");
      fflush(stdout);
      fgets(line, sizeof(line), stdin);
    } while (sscanf(line, "%d", numsteps) != 1);
  } while (*numsteps < 0);

  if (*numsteps > 1) {
    do {
      fprintf(stdout, "Increment in days: ");
      fflush(stdout);
      fgets(line, sizeof(line), stdin);
    } while (sscanf(line, "%lf", stepsize) != 1);
    printf("\n");
  } else {
    *stepsize = 0.0;
  }
}

/*********************************************************************
Name:    UseThisData
Purpose: Routine to supply default input data.  Data is taken
         from the paper by Mansfield (AIAA Paper No. 86-2269-CP)
         and is for Comet Halley. Elements are referred to equinox/ecliptic
         of B1950.
Inputs:  None.
Outputs: element.
         mass.
         StartingJD.
         stepsize.
         numsteps.
Returns: Nothing
Status:  Finished.
Errors:  None known.
*********************************************************************/
void UseThisData(double *element, double *mass, double *StartingJD,
  double *stepsize, int *numsteps) {

  /* Perihelion distance (AU) */
  element[0] = 0.5871047;
  /* Eccentricity */
  element[1] = 0.967276;
  /* Inclination (radians) */
  element[2] = 162.23928 * D2R;
  /* Longitude of ascending node (radians) */
  element[3] = 58.14536 * D2R;
  /* Argument of perihelion (radians) */
  element[4] = 111.84809 * D2R;
  /* Julian date of perihelion passage (TDT) */
  element[5] = 2446470.95175;
  /* Mass (solar masses) */
  *mass = (double)0.0;

  *StartingJD = (double)2446445.5;
  *stepsize   = (double)15.0;
  *numsteps   = (int)5;
}

/* End Of File - eph.c **********************************************/
