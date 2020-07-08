/*********************************************************************
Project:  pretran
Filename: pretran.c
Author:   Joe Heafner.
Purpose:  Example driver for GetPrecessParams.
*********************************************************************/

/* Header Files *****************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "astrolib.h"

/*********************************************************************
Name:    main
Purpose: Main routine for pretran.
Inputs:  None.
Outputs: None.
Returns: 0 if execution successful.
Status:  Finished.
Errors:  None known.
*********************************************************************/
int main(void) {

  char line[1024];
  char epoch1[30], epoch2[30], *ptr;
  double jed1, jed2, r1[3], r2[3];
  double zeta=0.0, z=0.0, theta=0.0, zetadot=0.0, zdot=0.0, thetadot=0.0;

  printf("EXAMPLE DRIVER FOR GetPrecessParams\n\n");

  do {
    fprintf(stdout, "Enter initial epoch or Julian date: ");
    fflush(stdout);
    fgets(epoch1, sizeof(epoch1), stdin);
    if ((ptr = strrchr(epoch1, '\n')) != NULL)
      *ptr = '\0';
    ucase(epoch1);
    Trim(epoch1);
  } while (strlen(epoch1) < 1);

  do {
    fprintf(stdout, "Enter final epoch or Julian date: ");
    fflush(stdout);
    fgets(epoch2, sizeof(epoch2), stdin);
    if ((ptr = strrchr(epoch2, '\n')) != NULL)
      *ptr = '\0';
    ucase(epoch2);
    Trim(epoch2);
  } while (strlen(epoch2) < 1);

  if ((epoch1[0] == 'J') || (epoch1[0] == 'B'))
    Epoch2JED(epoch1, &jed1);
  else
    jed1 = atof(epoch1);

  if ((epoch2[0] == 'J') || (epoch2[0] == 'B'))
    Epoch2JED(epoch2, &jed2);
  else
    jed2 = atof(epoch2);

  do {
    fprintf(stdout, "Enter initial position vector X component: ");
    fflush(stdout);
    fgets(line, sizeof(line), stdin);
  } while (sscanf(line, "%lf", &r1[0]) != 1);

  do {
    fprintf(stdout, "Enter initial position vector Y component: ");
    fflush(stdout);
    fgets(line, sizeof(line), stdin);
  } while (sscanf(line, "%lf", &r1[1]) != 1);

  do {
    fprintf(stdout, "Enter initial position vector Z component: ");
    fflush(stdout);
    fgets(line, sizeof(line), stdin);
  } while (sscanf(line, "%lf", &r1[2]) != 1);

  printf("\nINITIAL POSITION VECTOR COMPONENTS\n");
  printf("%+13.10f %+13.10f %+13.10f\n\n", r1[0], r1[1], r1[2]);

  /* Compute the precessional angles */
  GetPrecessParams(jed1, jed2, &zeta, &z, &theta, &zetadot, &zdot, &thetadot);

  /* Perform the matrix rotations one at a time.    */
  /* Note the order of application of the matrices. */
  MRotate(r1, 3, -zeta, r1);
  MRotate(r1, 2, theta, r1);
  MRotate(r1, 3, -z, r2);

  printf("PRECESSED POSITION VECTOR COMPONENTS\n");
  printf("%+13.10f %+13.10f %+13.10f\n", r2[0], r2[1], r2[2]);

  printf("\nOK\n");
  return (0);
}

/* End Of File - pretran.c ******************************************/
