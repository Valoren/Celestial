/*********************************************************************
Project:  nuttest
Filename: nuttest.c
Author:   Joe Heafner.
Purpose:  Example driver for GetDpsiDeps.
*********************************************************************/

/* Header Files *****************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "astrolib.h"

/*********************************************************************
Name:    main
Purpose: Main routine for nuttest.
Inputs:  None.
Outputs: None.
Returns: 0 if execution successful.
Status:  Finished.
Errors:  None known.
*********************************************************************/
int main(void) {

  double jdate, dpsi, deps, dpsidot, depsdot, inc = 1.0;
  double MeanEps, MeanEpsDot, TrueEps, TrueEpsDot, TEps, MEps;
  int NumSteps = 0, i;
  char line[1024], temp[30], right_buffer[30];
  char m[30], t[30], dp[30], de[30];

  printf("EXAMPLE DRIVER FOR GetDpsiDeps\n\n");

  do {
    fprintf(stdout, "Enter Julian date on TDB scale: ");
    fflush(stdout);
    fgets(line, sizeof(line), stdin);
  } while (sscanf(line, "%lf", &jdate) != 1);

  do {
    do {
      fprintf(stdout, "Enter number of steps (1 - 15): ");
      fflush(stdout);
      fgets(line, sizeof(line), stdin);
    } while (sscanf(line, "%d", &NumSteps) != 1);
  } while ((NumSteps < 0) || (NumSteps > 15));

  if (NumSteps >= 1) {
    do {
      do {
        fprintf(stdout, "Enter increment in days: ");
        fflush(stdout);
        fgets(line, sizeof(line), stdin);
      } while (sscanf(line, "%lf", &inc) != 1);
    } while (inc < 1.0);
  } else {
    NumSteps = 0;
  }

  printf("   JED        Mean Obl.    True     dpsi      deps\n");
  printf("-----------------------------------------------------\n");

  i = 0;

  do {
    /* Get nutation parameters */
    GetDpsiDeps(jdate, &dpsi, &deps, &dpsidot, &depsdot);
    /* Get mean obliquity */
    Obliquity(J2000, jdate, 0, &MeanEps, &MeanEpsDot);
    /* Get true obliquity */
    Obliquity(J2000, jdate, 1, &TrueEps, &TrueEpsDot);

    TEps = TrueEps * R2D;
    MEps = MeanEps * R2D;
    dpsi = dpsi * R2D;
    deps = deps * R2D;

    FmtDms(MEps, 3, 0, m);

    FmtDms(TEps, 3, 0, temp);
    right(temp, 7, t);

    FmtDms(dpsi, 4, 0, temp);
    left(temp, 1, dp);
    right(temp, 8, right_buffer);
    strcat(dp, right_buffer);

    FmtDms(deps, 4, 0, temp);
    left(temp, 1, de);
    right(temp, 8, right_buffer);
    strcat(de, right_buffer);

    printf("%9.1f  %s %s %s %s\n", jdate, m, t, dp, de);

    jdate += inc;
    i++;
  } while (i < NumSteps);

  printf("\nOK\n");
  return (0);
}

/* End Of File - nuttest.c ******************************************/
