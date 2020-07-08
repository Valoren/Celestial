/*********************************************************************
Project:  sidtest
Filename: sidtest.c
Author:   Joe Heafner.
Purpose:  Example driver for GetGST.
*********************************************************************/

/* Header Files *****************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "astrolib.h"

/*********************************************************************
Name:    main
Purpose: Main routine for sidtest.
Inputs:  None.
Outputs: None.
Returns: 0 if execution successful.
Status:  Finished.
Errors:  None known.
*********************************************************************/
int main(void) {

  double jdate;
  int NumSteps = 0, i;
  char line[1024], temp1[30], temp2[30];
  double gmst, gast, inc = 1.0;

  printf("EXAMPLE DRIVER FOR GetGST\n\n");

  do {
    fprintf(stdout, "Enter starting Julian date on TDB scale: ");
    fflush(stdout);
    fgets(line, sizeof(line), stdin);
  } while (sscanf(line, "%lf", &jdate) != 1);

  do {
    do {
      fprintf(stdout, "Enter number of steps (1 - 15): ");
      fflush(stdout);
      fgets(line, sizeof(line), stdin);
    } while (sscanf(line, "%d", &NumSteps) != 1);
  } while ((NumSteps < 1) || (NumSteps > 15));

  if (NumSteps > 1) {
    do {
      do {
        fprintf(stdout, "Enter increment in days: ");
        fflush(stdout);
        fgets(line, sizeof(line), stdin);
      } while (sscanf(line, "%lf", &inc) != 1);
    } while (inc < 0.0);
  } else {
    inc = 0.0;
  }

  printf("\n");

  i = 0;
  do {
    GetGST(jdate, 0, &gmst);
    gmst = gmst * R2H;
    GetGST(jdate, 1, &gast);
    gast = gast * R2H;
    FmtDms(gmst,4,2, temp1);
    FmtDms(gast,4,2, temp2);
    printf("JED %9.1f  GMST %s  GAST %s\n", jdate, temp1, temp2);
    jdate += inc;
    i++;
  } while (i < NumSteps);

  printf("\nOK\n");
  return (0);
}

/* End Of File - sidtest.c ******************************************/
