/*********************************************************************
Project:  jdtest
Filename: jdtest.c
Author:   Joe Heafner.
Purpose:  Example driver for Cal2JED.
*********************************************************************/

/* Header Files *****************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "astrolib.h"

/*********************************************************************
Name:    main
Purpose: Main routine for jdtest.
Inputs:  None.
Outputs: None.
Returns: 0 if execution successful.
Status:  Finished.
Errors:  None known.
*********************************************************************/
int main(void) {

  int day, month, year, w, s;
  double utc, tai_utc, ut1_utc, jed, jd[5];
  char line[1024], m;

  printf("EXAMPLE DRIVER FOR Cal2JED\n\n");

  do {
    fprintf(stdout, "Enter Day: ");
    fflush(stdout);
    fgets(line, sizeof(line), stdin);
  } while (sscanf(line, "%d", &day) != 1);

  do {
    fprintf(stdout, "Enter Month: ");
    fflush(stdout);
    fgets(line, sizeof(line), stdin);
  } while (sscanf(line, "%d", &month) != 1);

  do {
    fprintf(stdout, "Enter Year: ");
    fflush(stdout);
    fgets(line, sizeof(line), stdin);
  } while (sscanf(line, "%d", &year) != 1);

  do {
    fprintf(stdout, "Enter UTC as hh.mmssss: ");
    fflush(stdout);
    fgets(line, sizeof(line), stdin);
  } while (sscanf(line, "%lf", &utc) != 1);

  do {
    fprintf(stdout, "Enter TAI-UTC in seconds: ");
    fflush(stdout);
    fgets(line, sizeof(line), stdin);
  } while (sscanf(line, "%lf", &tai_utc) != 1);

  do {
    fprintf(stdout, "Enter UT1-UTC in seconds: ");
    fflush(stdout);
    fgets(line, sizeof(line), stdin);
  } while (sscanf(line, "%lf", &ut1_utc) != 1);

  do {
    do {
      fprintf(stdout, "Do you want MJD's (y/n)?: ");
      fflush(stdout);
      fgets(line, sizeof(line), stdin);
    } while (sscanf(line, "%c", &m) != 1);
  } while ((toupper(m) != 'Y') && (toupper(m) != 'N'));

  if (toupper(m) == 'Y') {
    w = 1;
  } else {
    w = 0;
  }

  for (s = 1; s <= 5; s++) {
    Cal2JED(month, day, year, utc, s, tai_utc, ut1_utc, w, &jed);
    jd[s-1] = jed;
  }

  if (toupper(m) == 'Y') {
    printf("\nMODIFIED JULIAN DATES ON ALL FIVE TIME SCALES:\n\n");
  } else {
    printf("\nJULIAN DATES FOR DIFFERENT TIME SCALES:\n\n");
  }

  printf("%16.8f UT1\n", jd[0]);
  printf("%16.8f UT2\n", jd[1]);
  printf("%16.8f TDT\n", jd[2]);
  printf("%16.8f TDB\n", jd[3]);
  printf("%16.8f UTC\n", jd[4]);

  printf("\nOK\n");
  return (0);
}

/* End Of File - jdtest.c *******************************************/
