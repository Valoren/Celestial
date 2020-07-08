/*********************************************************************
Project:  rise
Filename: rise.c
Author:   Joe Heafner.
Purpose:  Example driver for RST.
*********************************************************************/

/* Header Files *****************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "astrolib.h"

/*********************************************************************
Name:    main
Purpose: Main routine for rise.
Inputs:  None.
Outputs: None.
Returns: 0 if execution successful.
Status:  Finished.
Errors:  None known.
*********************************************************************/
int main(void) {

  char line[1024];
  double jed, deltat, ra[3], dec[3], hp, sd, z0;
  char ris[100], trn[100], set[100];

  printf("EXAMPLE DRIVER FOR RST\n\n");

  do {
    fprintf(stdout, "Enter deltat in seconds: ");
    fflush(stdout);
    fgets(line, sizeof(line), stdin);
  } while (sscanf(line, "%lf", &deltat) != 1);

  do {
    fprintf(stdout, "Enter latitude (+N/-S, dd.mmss): ");
    fflush(stdout);
    fgets(line, sizeof(line), stdin);
  } while (sscanf(line, "%lf", &obsr_lat) != 1);
  obsr_lat = D2R * deg(obsr_lat);

  do {
    fprintf(stdout, "Enter longitude (-W/+E, dd.mmss)\n");
    fprintf(stdout, "Enter '9999' for Ephemeris Meridian: ");
    fflush(stdout);
    fgets(line, sizeof(line), stdin);
  } while (sscanf(line, "%lf", &obsr_lon) != 1);
  if (obsr_lon == 9999.0) {
    obsr_lon = 1.002738 * (deltat / 3600.0) * H2R;
  } else {
    obsr_lon = D2R * deg(obsr_lon);
  }

  /* The following data is must be replaced for specific */
  /* applications.                                       */
  jed = 2450093.5;  /* 1/11/1996 */
  ra[0]  = H2R * deg(10.1917027);
  dec[0] = D2R * deg( 6.291308);
  ra[1]  = H2R * deg(11.0558576);
  dec[1] = D2R * deg( 2.383610);
  ra[2]  = H2R * deg(11.5311971);
  dec[2] = D2R * deg(-1.223363);

  hp = D2R * deg(0.444156);
  sd = asin(0.272493 * sin(hp));
  z0 = D2R * deg(90.34) + sd - hp;

  RST(jed, ra, dec, z0, deltat, ris, trn, set);

  printf("\n");
  printf("%f\n", jed);
  printf("Rise %s Transit %s Set %s\n", ris, trn, set);

  printf("\nOK\n");
  return (0);
}

/* End Of File - rise.c *********************************************/
