/*********************************************************************
Project:  rpnmat
Filename: rpnmat.c
Author:   Joe Heafner.
Purpose:  Example driver for GetRPNmat.
*********************************************************************/

/* Header Files *****************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "astrolib.h"

/*********************************************************************
Name:   main
Purpose:  Main routine for rpnmat.
Inputs:   None.
Outputs:  None.
Returns:  0 if execution successful.
Status:   Finished.
Errors:   None known.
*********************************************************************/
int main(void) {

  int i, Which;
  char line[1024];
  double jed1, jed2;
  DMatrix m6, pmat;

  printf("EXAMPLE DRIVER FOR GetRPNmat\n\n");

  do {
    fprintf(stdout, "Enter JED1 on TDB scale: ");
    fflush(stdout);
    fgets(line, sizeof(line), stdin);
  } while (sscanf(line, "%lf", &jed1) != 1);

  do {
    fprintf(stdout, "Enter JED2 on TDB scale: ");
    fflush(stdout);
    fgets(line, sizeof(line), stdin);
  } while (sscanf(line, "%lf", &jed2) != 1);

  do {
    do {
      fprintf(stdout, "[1] Precession [2] Nutation [3] Combined: ");
      fflush(stdout);
      fgets(line, sizeof(line), stdin);
    } while (sscanf(line, "%d", &Which) != 1);
  } while ((Which < 1) || (Which > 3));

  printf("\n");

  pmat = createDMatrix(6);
  m6 = createDMatrix(6);

  GetRPNmat(jed1, jed2, Which, 1, pmat, m6);

  for (i = 0; i < 6; i++) {
    printf("%+11.8f %+11.8f %+11.8f %+11.8f %+11.8f %+11.8f\n",
      pmat[i][0], pmat[i][1], pmat[i][2], pmat[i][3], pmat[i][4], pmat[i][5]);
  }

  freeDMatrix(pmat, 6);
  freeDMatrix(m6, 6);

  printf("\nOK\n");
  return (0);
}

/* End Of File - rpnmat.c *******************************************/
