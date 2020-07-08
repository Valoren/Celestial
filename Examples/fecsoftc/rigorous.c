/*********************************************************************
Project:  rigorous
Filename: rigorous.c
Author:   Joe Heafner.
Purpose:  Example driver for rigorous element precession.
*********************************************************************/

/* Header Files *****************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "astrolib.h"

/*********************************************************************
Name:    main
Purpose: Main routine for rigorous.
Inputs:  None.
Outputs: None.
Returns: 0 if execution successful.
Status:  Finished.
Errors:  None known.
*********************************************************************/
int main(void) {

  char line[1024], *ptr;
  char epoch1[30], epoch2[30];
  double jed1, jed2;
  double I, node, arg, p[3], q[3], w[3];
  double newi, newnode, newarg;
  double zeta, z, theta, zetadot, zdot, thetadot;
  double eps1, eps2, epsdot;

  printf("EXAMPLE DRIVER FOR RIGOROUS ELEMENT PRECESSION\n\n");

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
    fprintf(stdout, "Enter inclination: ");
    fflush(stdout);
    fgets(line, sizeof(line), stdin);
  } while (sscanf(line, "%lf", &I) != 1);
  I = I * D2R;

  do {
    fprintf(stdout, "Enter node: ");
    fflush(stdout);
    fgets(line, sizeof(line), stdin);
  } while (sscanf(line, "%lf", &node) != 1);
  node = node * D2R;

  do {
    fprintf(stdout, "Enter arg peri: ");
    fflush(stdout);
    fgets(line, sizeof(line), stdin);
  } while (sscanf(line, "%lf", &arg) != 1);
  arg = arg * D2R;

  p[0] = cos(arg) * cos(node) - sin(arg) * sin(node) * cos(I);
  p[1] = cos(arg) * sin(node) + sin(arg) * cos(node) * cos(I);
  p[2] = sin(arg) * sin(I);

  q[0] = -sin(arg) * cos(node) - cos(arg) * sin(node) * cos(I);
  q[1] = -sin(arg) * sin(node) + cos(arg) * cos(node) * cos(I);
  q[2] =  cos(arg) * sin(I);

  w[0] =  sin(node) * sin(I);
  w[1] = -cos(node) * sin(I);
  w[2] =  cos(I);

  /* Compute the equatorial precessional angles */
  GetPrecessParams(jed1, jed2, &zeta, &z, &theta, &zetadot, &zdot, &thetadot);

  /* Precess p[] */
  Obliquity(J2000, jed1, 0, &eps1, &epsdot);
  Obliquity(J2000, jed2, 0, &eps2, &epsdot);
  MRotate (p, 1, -eps1  , p);
  MRotate (p, 3, -zeta  , p);
  MRotate (p, 2,  theta , p);
  MRotate (p, 3, -z     , p);
  MRotate (p, 1,  eps2  , p);

  /* Precess q[] */
  MRotate (q, 1, -eps1  , q);
  MRotate (q, 3, -zeta  , q);
  MRotate (q, 2,  theta , q);
  MRotate (q, 3, -z     , q);
  MRotate (q, 1,  eps2  , q);

  /* Precess w[] */
  MRotate (w, 1, -eps1  , w);
  MRotate (w, 3, -zeta  , w);
  MRotate (w, 2,  theta , w);
  MRotate (w, 3, -z     , w);
  MRotate (w, 1,  eps2  , w);

  /* Get new elements */
  newnode = atan2(w[0], -w[1]);
  if (newnode < 0.0)
    newnode += TWOPI;

  newi    = atan2(sqrt(p[2]*p[2] + q[2]*q[2]), w[2]);
  if (newi < 0.0)
    newi += TWOPI;

  newarg  = atan2(p[2], q[2]);
  if (newarg < 0.0)
    newarg += TWOPI;

  printf("NEW ELEMENTS\n");
  printf(" %.13f\n", newi * R2D);
  printf(" %.13f\n", newnode * R2D);
  printf(" %.13f\n", newarg * R2D);

  printf("\nOK\n");
  return (0);
}

/* End Of File - rigorous.c *****************************************/
