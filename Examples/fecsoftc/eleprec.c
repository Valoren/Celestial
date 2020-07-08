/*********************************************************************
Project:  eleprec
Filename: eleprec.c
Author:   Joe Heafner.
Purpose:  Example driver for PrecessElements.
*********************************************************************/

/* Header Files *****************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "astrolib.h"

/*********************************************************************
Name:    main
Purpose: Main routine for eleprec.
Inputs:  None.
Outputs: None.
Returns: 0 if execution successful.
Status:  Finished.
Errors:  None known.
*********************************************************************/
int main(void) {

  char line[1024];
  char Eqnx1[50], Eqnx2[50];
  double Inclination, Node, ArgPeri, element[3];

  printf("EXAMPLE DRIVER FOR PrecessElements\n\n");

  fprintf(stdout, "First equinox: ");
  fflush(stdout);
  fgets(Eqnx1, sizeof(Eqnx1), stdin);
  Trim(Eqnx1);

  fprintf(stdout, "Final equinox: ");
  fflush(stdout);
  fgets(Eqnx2, sizeof(Eqnx2), stdin);
  Trim(Eqnx2);

  do {
    fprintf(stdout, "Inclination: ");
    fflush(stdout);
    fgets(line, sizeof(line), stdin);
  } while (sscanf(line, "%lf", &Inclination) != 1);

  do {
    fprintf(stdout, "Node: ");
    fflush(stdout);
    fgets(line, sizeof(line), stdin);
  } while (sscanf(line, "%lf", &Node) != 1);

  do {
    fprintf(stdout, "Arg. peri.: ");
    fflush(stdout);
    fgets(line, sizeof(line), stdin);
  } while (sscanf(line, "%lf", &ArgPeri) != 1);

  element[0] = Inclination * D2R;
  element[1] = Node * D2R;
  element[2] = ArgPeri * D2R;

  PrecessElements(Eqnx1, element, Eqnx2);

  printf("\nPRECESSED ELEMENTS\n");

  Inclination = element[0] * R2D;
  Node = element[1] * R2D;
  ArgPeri = element[2] * R2D;

  printf("Inclination %+10.5f\n", Inclination);
  printf("Node        %+10.5f\n", Node);
  printf("Arg. peri.  %+10.5f\n", ArgPeri);

  printf("\nOK\n");
  return (0);
}

/* End Of File - eleprec.c ******************************************/
