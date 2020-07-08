/*********************************************************************
Project:  geobs
Filename: geobs.c
Author:   Joe Heafner.
Purpose:  Example driver for GeocenObs.
*********************************************************************/

/* Header Files *****************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "astrolib.h"

/*********************************************************************
Name:    main
Purpose: Main routine for geobs.
Inputs:  None.
Outputs: None.
Returns: 0 if execution successful.
Status:  Finished.
Errors:  None known.
*********************************************************************/
int main(void) {

  double jdate=0.0, longitude=0.0, latitude=0.0, elevation=0.0, obsr_geo[6];
  char line[1024];

  printf("EXAMPLE DRIVER FOR GeocenObs\n\n");

  do {
    fprintf(stdout, "Enter Julian day number on TDB scale: ");
    fflush(stdout);
    fgets(line, sizeof(line), stdin);
  } while (sscanf(line, "%lf", &jdate) != 1);

  do {
    fprintf(stdout, "Enter observer's geodetic longitude (-W, dd.mmss): ");
    fflush(stdout);
    fgets(line, sizeof(line), stdin);
  } while (sscanf(line, "%lf", &longitude) != 1);

  do {
    fprintf(stdout, "Enter observer's geodetic latitude  (+N, dd.mmss): ");
    fflush(stdout);
    fgets(line, sizeof(line), stdin);
  } while (sscanf(line, "%lf", &latitude) != 1);

  do {
    fprintf(stdout, "Enter observer's elevation in meters: ");
    fflush(stdout);
    fgets(line, sizeof(line), stdin);
  } while (sscanf(line, "%lf", &elevation) != 1);

  longitude = deg(longitude) * D2R;
  latitude  = deg(latitude)  * D2R;

  obsr_lon = longitude;
  obsr_lat = latitude;
  obsr_ele = elevation;

  GeocenObs(jdate, obsr_geo);

  printf("Observer's state vector\n");
  printf("%+17.15f %+17.15f %+17.15f AU\n",
    obsr_geo[0], obsr_geo[1], obsr_geo[2]);
  printf("%+17.15f %+17.15f %+17.15f AU/day\n",
    obsr_geo[3], obsr_geo[4], obsr_geo[5]);

  printf("\nOK\n");
  return (0);
}

/* End Of File - geobs.c ********************************************/
