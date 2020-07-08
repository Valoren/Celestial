/*********************************************************************
Project:  coords
Filename: coords.c
Author:   Joe Heafner.
Purpose:  Example driver for Eq2Ecl and Eq2Hor.
*********************************************************************/

/* Header Files *****************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "astrolib.h"

/*********************************************************************
Name:    main
Purpose: Main routine for coords.
Inputs:  None.
Outputs: None.
Returns: 0 if execution successful.
Status:  Finished.
Errors:  None known.
*********************************************************************/
int main(void) {

  char line[1024];
  int   Which;
  double  Alpha, Delta, LAST, Altitude, Azimuth;
  double  Lambda, Beta, TrueEps, HourAngle, r1[6], r2[6];
  char  RA[40], Del[40], Az[40], Alt[40], Lam[40], Bet[40];

  printf("EXAMPLE DRIVER FOR Eq2Ecl and Eq2Hor\n\n");

  do {
    do {
      fprintf(stdout, "[1] EQ TO HOR [2] HOR TO EQ\n");
      fprintf(stdout, "[3] EQ TO ECL [4] ECL TO EQ\n");
      fprintf(stdout, "Select option (1-4): ");
      fflush(stdout);
      fgets(line, sizeof(line), stdin);
    } while (sscanf(line, "%d", &Which) != 1);
  } while ((Which < 1) || (Which > 4));

  switch (Which) {
  case 1:
    do {
      fprintf(stdout,
        "Enter equatorial coordinates (ra, dec) in dd.mmss format: ");
      fflush(stdout);
      fgets(line, sizeof(line), stdin);
    } while (sscanf(line, "%lf, %lf", &Alpha, &Delta) != 2);

    do {
      fprintf(stdout, "Enter observer latitude in dd.mmss format: ");
      fflush(stdout);
      fgets(line, sizeof(line), stdin);
    } while (sscanf(line, "%lf", &obsr_lat) != 1);

    do {
      fprintf(stdout, "Enter LAST in hh.mmss format: ");
      fflush(stdout);
      fgets(line, sizeof(line), stdin);
    } while (sscanf(line, "%lf", &LAST) != 1);

    Alpha = deg(Alpha) * H2R;
    Delta = deg(Delta) * D2R;
    obsr_lat = deg(obsr_lat) * D2R;
    LAST = deg(LAST) * H2R;
    HourAngle = LAST - Alpha;
    r1[0] = HourAngle;
    r1[1] = Delta;
    r1[2] = 1.0;
    r1[3] = 0.0;
    r1[4] = 0.0;
    r1[5] = 0.0;

    /* Change to rectangular variables */
    Pol2Rec(r1, r2);
    /* Transform */
    Eq2Hor(r2, 0, r2);
    /* Change to polar variables */
    Rec2Pol(r2, r2);

    Azimuth  = r2[0] * R2D;
    Altitude = r2[1] * R2D;

    FmtDms(Azimuth, 4, 0, Az);
    FmtDms(Altitude, 4, 0, Alt);
    Alpha = Alpha * R2H;
    Delta = Delta * R2D;
    FmtDms(Alpha, 4, 1, RA);
    FmtDms(Delta, 4, 0, Del);

    printf("\n");
    printf("EQUATORIAL COORDINATES\n");
    printf("RA %s DELTA %s\n", RA, Del);
    printf("HORIZON COORDINATES\n");
    printf("AZIMUTH %s ALTITUDE %s\n", Az, Alt);
    break;
  case 2:
    do {
      fprintf(stdout,
        "Enter horizon coordinates (alt, az) in dd.mmss format: ");
      fflush(stdout);
      fgets(line, sizeof(line), stdin);
    } while (sscanf(line, "%lf, %lf", &Altitude, &Azimuth) != 2);

    do {
      fprintf(stdout, "Enter observer latitude in dd.mmss format: ");
      fflush(stdout);
      fgets(line, sizeof(line), stdin);
    } while (sscanf(line, "%lf", &obsr_lat) != 1);

    do {
      fprintf(stdout, "Enter LAST in hh.mmss format: ");
      fflush(stdout);
      fgets(line, sizeof(line), stdin);
    } while (sscanf(line, "%lf", &LAST) != 1);

    Altitude = deg(Altitude) * D2R;
    Azimuth = deg(Azimuth) * D2R;
    obsr_lat = deg(obsr_lat) * D2R;
    LAST = deg(LAST) * H2R;
    r1[0] = Azimuth;
    r1[1] = Altitude;
    r1[2] = 1.0;
    r1[3] = 0.0;
    r1[4] = 0.0;
    r1[5] = 0.0;

    /* Change to rectangular variables */
    Pol2Rec(r1, r2);
    /* Transform */
    Eq2Hor(r2, 1, r2);
    /* Change to polar variables */
    Rec2Pol(r2, r2);

    HourAngle = r2[0];
    Delta = r2[1];

    Azimuth  = Azimuth * R2D;
    Altitude = Altitude * R2D;

    FmtDms(Azimuth, 4, 0, Az);
    FmtDms(Altitude, 4, 0, Alt);

    Alpha = LAST - HourAngle;
    Alpha = amodulo(Alpha,TWOPI) * R2H;
    Delta = Delta * R2D;
    FmtDms(Alpha, 4, 1, RA);
    FmtDms(Delta, 4, 0, Del);

    printf("\n");
    printf("HORIZON COORDINATES\n");
    printf("AZIMUTH %s ALTITUDE %s\n", Az, Alt);
    printf("EQUATORIAL COORDINATES\n");
    printf("RA %s DELTA %s\n", RA, Del);
    break;
  case 3:
    do {
      fprintf(stdout,
        "Enter equatorial coordinates (ra, dec) in dd.mmss format: ");
      fflush(stdout);
      fgets(line, sizeof(line), stdin);
    } while (sscanf(line, "%lf, %lf", &Alpha, &Delta) != 2);

    do {
      fprintf(stdout, "Enter true obliquity in dd.mmss format: ");
      fflush(stdout);
      fgets(line, sizeof(line), stdin);
    } while (sscanf(line, "%lf", &TrueEps) != 1);

    Alpha = deg(Alpha) * H2R;
    Delta = deg(Delta) * D2R;
    TrueEps = deg(TrueEps) * D2R;
    r1[0] = Alpha;
    r1[1] = Delta;
    r1[2] = 1.0;
    r1[3] = 0.0;
    r1[4] = 0.0;
    r1[5] = 0.0;

    /* Change to rectangular variables */
    Pol2Rec(r1, r2);
    /* Transform */
    Eq2Ecl(r2, 0, TrueEps, r2);
    /* Change to polar variables */
    Rec2Pol(r2, r2);

    Lambda = r2[0] * R2D;
    Beta = r2[1] * R2D;

    FmtDms(Lambda, 4, 0, Lam);
    FmtDms(Beta, 4, 0, Bet);

    Alpha = Alpha * R2H;
    Delta = Delta * R2D;
    FmtDms(Alpha, 4, 1, RA);
    FmtDms(Delta, 4, 0, Del);

    printf("\n");
    printf("EQUATORIAL COORDINATES\n");
    printf("RA %s DELTA %s\n", RA, Del);
    printf("ECLIPTIC COORDINATES\n");
    printf("LAMBDA %s BETA %s\n", Lam, Bet);
    break;
  case 4:
    do {
      fprintf(stdout,
        "Enter ecliptic coordinates (lon, lat) in dd.mmss format: ");
      fflush(stdout);
      fgets(line, sizeof(line), stdin);
    } while (sscanf(line, "%lf, %lf", &Lambda, &Beta) != 2);

    do {
      fprintf(stdout, "Enter true obliquity in dd.mmss format: ");
      fflush(stdout);
      fgets(line, sizeof(line), stdin);
    } while (sscanf(line, "%lf", &TrueEps) != 1);

    Lambda = deg(Lambda) * D2R;
    Beta = deg(Beta) * D2R;
    TrueEps = deg(TrueEps) * D2R;
    r1[0] = Lambda;
    r1[1] = Beta;
    r1[2] = 1.0;
    r1[3] = 0.0;
    r1[4] = 0.0;
    r1[5] = 0.0;

    /* Change to rectangular variables */
    Pol2Rec(r1, r2);
    /* Transform */
    Eq2Ecl(r2, 1, TrueEps, r2);
    /* Change to polar variables */
    Rec2Pol(r2, r2);

    Alpha = r2[0] * R2H;
    Delta = r2[1] * R2D;

    Lambda = Lambda * R2D;
    Beta = Beta * R2D;
    FmtDms(Lambda, 4, 0, Lam);
    FmtDms(Beta, 4, 0, Bet);

    FmtDms(Alpha, 4, 1, RA);
    FmtDms(Delta, 4, 0, Del);

    printf("\n");
    printf("ECLIPTIC COORDINATES\n");
    printf("LAMBDA %s BETA %s\n", Lam, Bet);
    printf("EQUATORIAL COORDINATES\n");
    printf("RA %s DELTA %s\n", RA, Del);
    break;
  }

  printf("\nOK\n");
  return (0);
}

/* End Of File - coords.c *******************************************/
