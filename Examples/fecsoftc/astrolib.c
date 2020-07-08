/*********************************************************************
Project:  fecsoftc
Filename: astrolib.c
Author:   Joe Heafner
Purpose:  Astronomical ephemeris library.
Thanks to Charles Gamble for extensive modifications.
*********************************************************************/

/* Header Files *****************************************************/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include "astrolib.h"
#include "support.h"

#ifndef MAX_FILENAME
  #define MAX_FILENAME  255
#endif

/* Globals **********************************************************/
/* Ephemeris header information */
char ttl[3][65], cnam[400][7];
int  ncon;
double SS[3], cval[400], au, emrat;
short int NUMDE, ipt[3][12], lpt[3];

/* Other data file specs not in header */
long LengthOfFile, ncoeff, BlockLength;
int LengthOfHeader, NumBlocks, bary, bsav;
double db[1100], pvsun[6];

static FILE *fpBinaryFile = NULL;

/* Temporary data used when reading in binary file */
short tmpShort;
double tmpDouble;
int tmpInt;
long tmpLong;

/*********************************************************************
Name:    Aberrate
Purpose: Function to correct an input vector for aberration.
Inputs:  p1    - Zero-offset geocentric state vector of object.
         EBdot - Zero-offset barycentric velocity of Earth.
Outputs: p2 - Zero-offset geocentric state vector of object corrected
              for aberration.
Returns: Nothing.
Status:  Finished.
Errors:  None known.
*********************************************************************/
void Aberrate(double p1[], double EBdot[], double p2[]) {

  double v[3], up[3], p[3], pdot[3], magP, magV, p1dotv, beta;
  int i;

  /* Extract the pos. portion of the state vector */
  SplitStateVector(p1, p, pdot);
  magP = Vecmag(p);

  /* Need to make Ppos() a unit vector */
  Uvector(p, up);

  for (i=0; i<3; i++) {
    v[i] = EBdot[i] / CAUD;
  }

  Vdot(3, up, v, &p1dotv);
  magV = Vecmag(v);

  beta = 1.0 / sqrt(1.0 - magV * magV);

  for (i=0; i<3; i++) {
    p2[i] = ((up[i] / beta) +
      (1.0 + p1dotv / (1.0 + (1.0 / beta))) * v[i]) / (1.0 + p1dotv);
    /* Make p2[] a non-unit vector */
    p2[i] = magP * p2[i];
    p2[i+3] = p1[i+3];
  }

  return;
}

/*********************************************************************
Name:    amodulo
Purpose: Function to reduce a to the range 0 <= a < b.
Inputs:  a - See above.
         b - See above.
Outputs: None.
Returns: See above.
Status:  Finished.
Errors:  None known.
*********************************************************************/
double amodulo(double a, double b) {

  double x;

  x = a - b * floor(a/b);
  return (x);
}

/*********************************************************************
Name:    Cal2JED
Purpose: Function to compute the Julian date on the specified time scale.
         The Julian day algorithm of Meeus is used, and the algorithm for
         converting to TDB is due to Hirayama et al. and can be found in
         REFERENCE FRAMES IN ASTRONOMY AND GEOPHYSICS ed. by Kovalevsky
         et al., 1989, pp. 439-442.
Inputs:  m       - Month.
         d       - Day.
         y       - Year.
         utc     - UTC as hh.mmssss.
         s       - 1 for UT1, 2 for UT2, 3 for TDT,
                   4 for TDB, 5 for UTC.
         tai_utc - TAI-UTC in seconds.
         ut1_utc - UT1-UTC in seconds.
         w       - 1 for MJD, 0 otherwise.
Outputs: jed - Appropriate JED.
Returns: Nothing.
Status:  Finished.
Errors:  None known.
*********************************************************************/
void Cal2JED(int m, int d, int y, double utc, int s,
  double tai_utc, double ut1_utc, int w, double *jed) {
  double time, datnum, a, b, term1, corr, T;

  /* Convert to decimal hours */
  time = deg(utc);

  if ((m == 1) || (m == 2)) {
    y = y - 1;
    m = m + 12;
  }

  /*
      Test to see if date is in Gregorian calendar.
      Converts date to a number of the form YYYY.MMDD.
  */
  datnum = (double) y + 0.01 * (double) m + 0.0001 * (double) d;

  if (datnum >= 1582.1015) {
    /* Gregorian calendar */
    a = fix(0.01 * (double) y);
    b = 2.0 - a + fix(0.25 * (double) a);
  } else {
    /* Julian calendar */
    a = 0.0;
    b = 0.0;
  }

  if (y < 0) {
    /* Handles negative years */
    term1 = fix(365.25 * (double) y - 0.75); /* change here */
  } else {
    term1 = fix(365.25 * (double) y);
  }

  *jed = term1 + fix(30.6001 * ((double) m + 1.0)) +
    (double) d + time / 24.0 + 1720994.5 + (double) b;

  switch (s) {
    case (1):
      corr = ut1_utc / 86400.0;
      break;
    case (2):
      corr = ut1_utc / 86400.0;
      *jed = *jed + corr;
      /* Compute date in Besselian years */
      T = 2000.0 + (*jed - 2451544.5333981) / 365.242198781;
      corr = 0.022 * sin(TWOPI * T);
      corr += -0.012 * cos(TWOPI * T);
      corr += -0.006 * sin(2.0 * TWOPI * T);
      corr += 0.007 * cos(2.0 * TWOPI * T);
      corr = corr / 86400.0;
      break;
    case (3):
      corr = (tai_utc + 32.184) / 86400.0;
      break;
    case (4):
      /* First convert to TDT */
      corr = (tai_utc + 32.184) / 86400.0;
      *jed = *jed + corr;
      T = (*jed - J2000) / JulCty;
      /* Now compute the new correction in microseconds */
      corr = 1656.675     * sin((35999.3729 * T + 357.5387) * D2R);
      corr +=  22.418     * sin((32964.467  * T + 246.199)  * D2R);
      corr +=  13.84      * sin((71998.746  * T + 355.057)  * D2R);
      corr +=   4.77      * sin(( 3034.906  * T +  25.463)  * D2R);
      corr +=   4.67      * sin((34777.259  * T + 230.394)  * D2R);
      corr +=  10.216 * T * sin((35999.373  * T + 243.451)  * D2R);
      corr +=   0.171 * T * sin((71998.746  * T + 240.98 )  * D2R);
      corr +=   0.027 * T * sin(( 1222.114  * T + 194.661)  * D2R);
      corr +=   0.027 * T * sin(( 3034.906  * T + 336.061)  * D2R);
      corr +=   0.026 * T * sin((  -20.186  * T +   9.382)  * D2R);
      corr +=   0.007 * T * sin((29929.562  * T + 264.911)  * D2R);
      corr +=   0.006 * T * sin((  150.678  * T +  59.775)  * D2R);
      corr +=   0.005 * T * sin(( 9037.513  * T + 256.025)  * D2R);
      corr +=   0.043 * T * sin((35999.373  * T + 151.121)  * D2R);
      /* Convert from microseconds to seconds */
      corr = corr * 0.000001;
      /* Convert to days */
      corr = corr / 86400.0;
      break;
    case (5):
      corr = 0.0;
      break;
    default:
      corr = 0.0;
  }

  *jed += corr;

  /* Return modified JED if requested */
  if (w == 1) *jed -= 2400000.5;
}

/*********************************************************************
Name:    constants
Purpose: Function to display the constants from the ephemeris file.
         Also prints other useful information about the file.
Inputs:  FileName - Filename of file to be read in.
Outputs: None.
Returns: Nothing.
Status:  Finished.
Errors:  None known.
*********************************************************************/
void constants(char *FileName) {

  int i;
  char Y[1024] = "";

  /* Open the file and read the header data */
  if (ephopn(FileName) == NULL) {
    LogMsg(stderr, "An error occurred in ephopn().\n");
    LogClose();
    exit(1);
  }

  for (i=0; i<3; i++)
    LogMsg(stdout, "%s\n", ttl[i]);

  LogMsg(stdout, "First valid JED in file  %9.1lf\n", SS[0]);
  LogMsg(stdout, "Final valid JED in file  %9.1lf\n", SS[1]);
  LogMsg(stdout, "Block interval  %2.0lf days\n",     SS[2]);
  LogMsg(stdout,
    "Length of file                 % 8.0ld bytes\n", LengthOfFile);
  LogMsg(stdout,
    "Length of header               % 8.0d bytes\n" , LengthOfHeader);
  LogMsg(stdout,
    "Length of each block           % 8.0ld bytes\n", BlockLength);
  LogMsg(stdout,
    "Coeffs per block               % 8.0ld\n"      , ncoeff);
  LogMsg(stdout,
    "Blocks in file                 % 8.0d\n"       , NumBlocks);

  LogMsg(stdout, "Bodies present in data file:  ");
  if (ipt[0][0] != 0) LogMsg(stdout, "MER ");
  if (ipt[0][1] != 0) LogMsg(stdout, "VEN ");
  if (ipt[0][2] != 0) LogMsg(stdout, "EMB ");
  if (ipt[0][3] != 0) LogMsg(stdout, "MAR ");
  if (ipt[0][4] != 0) LogMsg(stdout, "JUP ");
  if (ipt[0][5] != 0) LogMsg(stdout, "SAT ");
  if (ipt[0][6] != 0) LogMsg(stdout, "URA ");
  if (ipt[0][7] != 0) LogMsg(stdout, "NEP ");
  if (ipt[0][8] != 0) LogMsg(stdout, "PLU ");
  if (ipt[0][9] != 0) LogMsg(stdout, "MOO ");
  if (ipt[0][10] != 0) LogMsg(stdout, "SUN\n");

  if (ipt[0][11] == 0 || ipt[1][11] == 0 || ipt[2][11] == 0)
    LogMsg(stdout, "This ephemeris does not contain nutations.\n");
  else
    LogMsg(stdout, "This ephemeris contains nutations.\n");

  if (lpt[0] == 0 || lpt[1] == 0 || lpt[2] == 0)
    LogMsg(stdout, "This ephemeris does not contain librations.\n");
  else
    LogMsg(stdout, "This ephemeris contains librations.\n");

  LogMsg(stdout, "This ephemeris has %d ephemeris constants.\n", ncon);

  do {
    fprintf(stdout, "Display constant names and values (Y/N)? ");
    fflush(stdout);
    fgets(Y, sizeof(Y), stdin);
    Trim(Y);
    ucase(Y);
  } while ((strcmp(Y, "Y") != 0) && (strcmp(Y, "N") != 0));

  if (strcmp(Y, "Y") == 0) {
    LogMsg(stdout, "\n");
    for (i=0; i<ncon; i++) {
      LogMsg(stdout, "%6s  %+16.15E      ", cnam[i], cval[i]);
      if ((i+1) % 2 == 0) {
        LogMsg(stdout, "\n");
      }
    }

    if ((i % 2) != 0) {
      LogMsg(stdout, "\n");
    }
  }
}

/*********************************************************************
Name:    Conway
Purpose: Function to solve Kepler's equation using the method of
         Conway-Laguerre for universal variables.
Inputs:  uele - Array of universal orbital elements.
         mu   - Gravitational constant for the body.
         jed  - Time for which computations are to be performed.
Outputs: r_cos_nu - r*cos(true anomaly).
         r_sin_nu - r*sin(true anomaly).
Returns: Nothing.
Status:  Finished.
Errors:  None known.
*********************************************************************/
void Conway(double uele[], double mu, double jed, double *r_cos_nu,
  double *r_sin_nu) {

  double tspp, alpha, s, x, c1, c2, c3, f, fp, fpp, term, ds;
  int niter;

  /* Compute time since perihelion passage in days */
  tspp = jed - uele[5];

  /* Compute alpha */
  alpha = mu * (1.0 - uele[1]) / uele[0];

  /* Initial guess for s */
  s = 0.0;

  /* Initialize iteration counter */
  niter = 0;

  do {
    /* Compute Stumpff functions */
    x = alpha * s * s;
    c1 = StumpffN(x, 1);
    c2 = StumpffN(x, 2);
    c3 = StumpffN(x, 3);

    f = uele[0] * s + mu * uele[1] * s * s * s * c3 - tspp;
    fp = uele[0] + mu * uele[1] * s * s * c2;
    fpp = mu * uele[1] * s * c1;

    /* Evaluate Laguerre's method */
    term = 16.0 * fp * fp - 20.0 * f * fpp;
    ds = -5.0 * f / (fp + fp/fabs(fp) * sqrt(fabs(term)));

    s = s + ds;

    niter = niter + 1;

    /* Check for convergence or more than ten iterations */
    if (niter > 10)
    {
      LogMsg(stderr, "Conway: more than ten iterations.\n");
      exit(1);
    }

    if (fabs(f) < 1e-12) break;

  } while(1);

  /* Now compute r_sin_nu and r_cos_nu */
  *r_sin_nu = sqrt(mu * uele[0] * (1.0 + uele[1])) * s * c1;
  *r_cos_nu = uele[0] - mu * s * s * c2;
}

/*********************************************************************
Name:    createDMatrix
Purpose: Creates an n x n double matrix.
Inputs:  n - Dimension of matrix to create.
Outputs: None.
Returns: DMatrix.
Status:  Finished.
Errors:  None known.
*********************************************************************/
DMatrix createDMatrix(int n) {

  DMatrix m;
  int i;

  m = calloc(n, sizeof(double*));
  if (m == NULL) {
    return NULL;
  }

  for (i=0; i<n; i++) {
    m[i] = calloc(n, sizeof(double));
    if (m[i] == NULL) {
      freeDMatrix(m, i); /* Avoids garbage */
      return NULL;
    }
  }

  return m;
}

/*********************************************************************
Name:    deg
Purpose: Converts dd.mmssss to dd.dddddd.
         Converts hh.mmssss to hh.hhhhhh.
Inputs:  x - Value to convert.
Outputs: None.
Returns: Converted value.
Status:  Finished.
Errors:  None known.
*********************************************************************/
double deg(double x) {

  double dd, d, fixdd, ddfixdd;

  if (x == 0.0) return (0.0);

  dd = fabs(x);
  fixdd = floor(dd);
  ddfixdd = dd - fixdd + 5.0e-10;  /* fudge factor */
  d = fixdd + floor(100.0 * ddfixdd) / 60.0;
  d = d + (10000.0 * ddfixdd - 100.0 * floor(100.0 * ddfixdd)) / 3600.0;

  return ((x / dd) * d);
}

/*********************************************************************
Name:    dms
Purpose: Converts dd.dddddd to dd.mmssss.
         Converts hh.hhhhhh to hh.mmssss.
Inputs:  x - Value to convert.
Outputs: None.
Returns: Converted value.
Status:  Finished.
Errors:  None known.
*********************************************************************/
double dms(double x) {

  double dd, fdd, dmfd, d;

  if (x == 0.0) return (0.0);

  dd = fabs(x);
  fdd = floor(dd);
  dmfd = dd - fdd + 5.0e-10;  /* fudge factor */
  d = fdd + floor(60.0 * dmfd) / 100.0 + 0.006 *
    (60.0 * dmfd - floor(60.0 * dmfd));

  return ((x / dd) * d);
}

/*********************************************************************
Name:    DRound
Purpose: Function to round a number to a given number of decimal places.
Inputs:  x - Value to round.
         n - Number of decimal places.
Outputs: None.
Returns: Rounded number.
Status:  Finished.
Errors:  None known.
*********************************************************************/
double DRound(double x, int n) {

  double a, y;
  int sgn;

  if (x > 0) {
    sgn = 1;
  }
  else if (x < 0) {
    sgn = -1;
  } else {
    sgn = 0;
  }

  a = pow(10.0, (double) n);
  y = floor((a * x) + ((double) sgn) * 0.5) / a;

  return y;
}

/*********************************************************************
Name:    ephopn
Purpose: Function to open a binary ephemeris data file for access,
         read the header information, and store it in memory for
         later use.
Inputs:  FileName - Filename of file to be read in.
Outputs: None.
Returns: NULL on error.
         Open file pointer otherwise.
Status:  Finished.
Errors:  None known.
*********************************************************************/
FILE *ephopn(char *FileName) {

  /* Make sure all variables read from the data file are global */
  int i, j;
  static int efirst;
  long curpos;

  if (!efirst) {
    /* Default file name is JPLEPH */
    if (strlen(FileName) == 0) {
      strcpy(FileName, "JPLEPH");
    }

    /*
        The existence of the data file MUST be confirmed
        in the calling program. It is NOT checked here!
    */
    if ((fpBinaryFile = fopen(FileName, "rb")) == NULL) {
      LogMsg(stderr,
        "EPHOPN: Can't open binary data file: %s\n", strerror(errno));
      return (NULL);
    }

    /* BEGINNING OF HEADER */
    /* 195 bytes */
    for (i = 0; i < 3; i++) {
      /* char TTL[3][65] */
      if (fread(ttl[i], sizeof(char), 65, fpBinaryFile) != 65) {
        LogMsg(stderr,
          "EPHOPN: Error reading binary data file: %s\n", strerror(errno));
        return (NULL);
      }
      ttl[i][64] = '\0';
    }

    /* 2 bytes */
    if (fread(&tmpShort, sizeof(short), 1, fpBinaryFile) != 1) {
      LogMsg(stderr,
        "EPHOPN: Error reading binary data file: %s\n", strerror(errno));
      return (NULL);
    }

    convert_little_endian((char *) &tmpShort, sizeof(short));
    ncon = (int) tmpShort;

    /* ncon*6 bytes */
    for (j = 0; j < ncon; ++j) {
      /* char CNAM[NCON][7] */
      if (fread(&cnam[j], sizeof(char), 6, fpBinaryFile) != 6) {
        LogMsg(stderr,
          "EPHOPN: Error reading binary data file: %s\n", strerror(errno));
        return (NULL);
      }
      cnam[j][6] = '\0';
    }

    /* 24 bytes, 8 each */
    for (j = 0; j < 3; j++) {
      /* double SS[3] */
      if (fread(&tmpDouble, sizeof(double), 1, fpBinaryFile) != 1) {
        LogMsg(stderr,
          "EPHOPN: Error reading binary data file: %s\n", strerror(errno));
        return (NULL);
      }
      convert_little_endian((char *) &tmpDouble, sizeof(double));
      SS[j] = (double) tmpDouble;
    }

    /* 16 bytes, 8 each */
    if (fread(&tmpDouble,    sizeof(double), 1, fpBinaryFile) != 1) {
      LogMsg(stderr,
        "EPHOPN: Error reading binary data file: %s\n", strerror(errno));
      return (NULL);
    }
    convert_little_endian((char *) &tmpDouble, sizeof(double));
    au = (double) tmpDouble;

    if (fread(&tmpDouble, sizeof(double), 1, fpBinaryFile) != 1) {
      LogMsg(stderr,
        "EPHOPN: Error reading binary data file: %s\n", strerror(errno));
      return (NULL);
    }
    convert_little_endian((char *) &tmpDouble, sizeof(double));
    emrat = (double) tmpDouble;

    /* 72 bytes */
    for (i = 0; i < 3; i++) {
      /* short IPT[3][12] */
      for (j = 0; j < 12; j++) {
        if (fread(&tmpShort, sizeof(short), 1, fpBinaryFile) != 1) {
          LogMsg(stderr,
            "EPHOPN: Error reading binary data file: %s\n", strerror(errno));
          return (NULL);
        }
        convert_little_endian((char *) &tmpShort, sizeof(short));
        ipt[i][j] = (short) tmpShort;
      }
    }

    /* 2 bytes */
    if (fread(&tmpShort, sizeof(short), 1, fpBinaryFile) != 1) {
      LogMsg(stderr,
        "EPHOPN: Error reading binary data file: %s\n", strerror(errno));
      return (NULL);
    }
    convert_little_endian((char *) &tmpShort, sizeof(short));
    NUMDE = (short) tmpShort;

    /* 6 bytes, 2 each */
    for (i = 0; i < 3; i++) {
      if (fread(&tmpShort, sizeof(short), 1, fpBinaryFile) != 1) {
        LogMsg(stderr,
          "EPHOPN: Error reading binary data file: %s\n", strerror(errno));
        return (NULL);
      }
      convert_little_endian((char *) &tmpShort, sizeof(short));
      lpt[i] = (short) tmpShort;
    }

    /* ncon*8 bytes */
    for (j = 0; j < ncon; j++) {
      if (fread(&tmpDouble, sizeof(double), 1, fpBinaryFile) != 1) {
        LogMsg(stderr,
          "EPHOPN: Error reading binary data file: %s\n", strerror(errno));
        return (NULL);
      }
      convert_little_endian((char *) &tmpDouble, sizeof(double));
      cval[j] = (double) tmpDouble;
    }
    /* END OF HEADER */

    /* This block used to be a separate function */
    /* but it didn't work correctly.             */
    if ((curpos = ftell(fpBinaryFile)) == -1L) {
      LogMsg(stderr, "EPHOPN: ftell() returned -1L.\n");
      return (NULL);
    }

    if (fseek(fpBinaryFile, 0L, SEEK_END) != 0) {
      LogMsg(stderr, "EPHOPN: fseek() failed.\n");
      return (NULL);
    }

    LengthOfFile = ftell(fpBinaryFile);

    if (fseek(fpBinaryFile, curpos, SEEK_SET) != 0) {
      LogMsg(stderr, "EPHOPN: fseek() failed.\n");
      return (NULL);
    }
    /* end block */

    LengthOfHeader = 317L + 14L * (long) ncon;

    /* Length of block of coeffs = 8*NCOEFF bytes */
    /* calculate number of coeffs per block */
    ncoeff = 0L;

    for (j = 0; j < 12; j++) {
      if (j == 11)
        ncoeff = ncoeff + (long) ipt[1][j] * (long) ipt[2][j] * 2L;
      else
        ncoeff = ncoeff + (long) ipt[1][j] * (long) ipt[2][j] * 3L;
    }

    ncoeff = ncoeff + (long) lpt[1] * (long) lpt[2] * 3L + 2L;
    BlockLength = ncoeff * 8L;
    NumBlocks = (LengthOfFile - (long) LengthOfHeader) / BlockLength;
    efirst = TRUE;
  }

  return (fpBinaryFile);
}

/*********************************************************************
Name:    Epoch2JED
Purpose: Function to convert an epoch to its corresponding JED.
Inputs:  epoch - J or B epoch (e.g. J2000).
Outputs: jed - Appropriate JED.
Returns: Nothing.
Status:  Finished.
Errors:  None known.
*********************************************************************/
void Epoch2JED(char *epoch, double *jed) {

  double date;
  char buffer[80];

  strcpy(buffer, epoch);

  ucase(buffer);

  sscanf(buffer, "%*c %lf", &date);

  if (epoch[0] == 'J') {
    /* Julian epoch */
    *jed = (date - 2000.0) * 365.25 + 2451545.0;
  } else {
    /* Besselian epoch */
    *jed = (date - 1900.0) * 365.242198781731 + 2415020.31352;
  }

  return;
}

/*********************************************************************
Name:    Eq2Ecl
Purpose: Subprogram to convert an equatorial state vector to an
         ecliptic state vector or the reverse transformation.
Inputs:  a   - Zero-offset input vector.
         s   - 0 for eq  -> ecl, 1 for ecl -> eq.
         eps - Mean or apparent obliquity.
Outputs: b - Zero-offset output vector.
Returns: Nothing.
Status:  Finished.
Errors:  None known.
*********************************************************************/
void Eq2Ecl(double a[], int s, double eps, double b[]) {

  DMatrix r, dr;

  r  = createDMatrix(3);
  dr = createDMatrix(3);

  if (s == 0) {
    /* Equatorial to ecliptic */
    RotMat(1,  eps, r, dr);
  } else {
    /* Ecliptic to equatorial */
    RotMat(1, -eps, r, dr);
  }

  MatXVec(r, a, b, 3, 3);

  freeDMatrix(r, 3);
  freeDMatrix(dr,3);

  return;
}

/*********************************************************************
Name:    Eq2Hor
Purpose: Function to convert an equatorial state vector to a
         horizon state vector or the reverse transformation.
Inputs:  a - Zero-offset input vector.
         s - 0 for eq -> hor, 1 for hor -> eq.
Outputs: b - Zero-offset output vector.
Returns: Nothing.
Status:  Finished.
Errors:  None known.
*********************************************************************/
void Eq2Hor(double a[], int s, double b[]) {

  DMatrix r2, r3, r, dr2, dr3, dr;

  r2 = createDMatrix(3);
  r3 = createDMatrix(3);
  r  = createDMatrix(3);
  dr2 = createDMatrix(3);
  dr3 = createDMatrix(3);
  dr = createDMatrix(3);

  if (s == 0) {
    /* Equatorial to horizon */
    RotMat(2, PIDIV2 - obsr_lat, r2, dr2);
    RotMat(3, -PI, r3, dr3);
    MatXMat(r3, r2, r, 3);
  } else {
    /* Horizon to equatorial */
    RotMat(3, PI, r3, dr3);
    RotMat(2, obsr_lat - PIDIV2, r2, dr2);
    MatXMat(r2, r3, r, 3);
  }

  MatXVec(r, a, b, 3, 3);

  freeDMatrix(r2, 3);
  freeDMatrix(r3, 3);
  freeDMatrix(r, 3);
  freeDMatrix(dr2, 3);
  freeDMatrix(dr3, 3);
  freeDMatrix(dr, 3);

  return;
}

/*********************************************************************
Name:    errprt
Purpose: Function to print error messages.
Inputs:  group   - Error code number.
         message - Error message.
Outputs: None.
Returns: Nothing.
Status:  Finished.
Errors:  None known.
*********************************************************************/
void errprt(int group, char *message) {

  LogMsg(stderr, "\nERROR #%d %s\n", group, message);
  exit(1);
}

/*********************************************************************
Name:    fix
Purpose: Function to return the integer part of a number. Basically,
         this function emulates the PowerBasic fix() function.
Inputs:  x - Decimal value.
Outputs: None.
Returns: Integral part of a number.
Status:  Finished.
Errors:  None known.
*********************************************************************/
double fix(double x) {

  if (x < 0) {
    return ceil(x);
  } else {
    return floor(x);
  }
}

/*********************************************************************
Name:    FmtDms
Purpose: Function to format angular and time quantities. Basically,
         this function is a pretty-print version of dms().
         Modified with permission from Rex Shudde's code by Joe Heafner.
Inputs:  x - Decimal value.
         n - Number of seconds decimal digits.
         m - 0 convert decimal degrees to:
               ddd mm ss      if n = 0
               ddd mm ss.f    if n = 1
               ddd mm ss.ff   if n = 2
               ddd mm ss.fff  if n = 3
               ddd mm ss.ffff if n = 4
         m - 1 converts decimal hours to:
               hh mm ss       if n = 0
               hh mm ss.f     if n = 1
               hh mm ss.ff    if n = 2
               hh mm ss.fff   if n = 3
               hh mm ss.ffff  if n = 4
Outputs: s - Destination for output string.
Returns: Nothing.
Status:  Finished.
Errors:  None known.
*********************************************************************/
void FmtDms(double x, int n, int m, char *s) {

  double absx, deg, min, sec;
  int nf;
  unsigned int fig;
  static int defd;
  static char dh[3];
  static char mm[3];
  static char ss[3];
  char buffer[16] = "              ";
  char right_buffer[16], left_buffer[16];

  strcpy(s, "              ");

  if (!defd) {
    defd++;
    dh[0] = '\xf8'; /* f8 is the hex code for degree symbol */
    dh[1] = 'h';
    dh[2] = '\0';
    mm[0] = '\x27'; /* 27 is the hex code for arcmin symbol */
    mm[1] = 'm';
    mm[2] = '\0';
    ss[0] = '\x22'; /* 22 is the hex code for arcsec symbol */
    ss[1] = 's';
    ss[2] = '\0';
  }

  absx = fabs(x);

  /* determine how many digits before the decimal */
  if (absx < 100) {
    fig = 2;
  } else {
    fig = 3;
  }

  deg = floor(absx);
  absx = 60.0 * (absx - deg);
  min = floor(absx);
  sec = 60.0 * (absx - min);
  sec = DRound(sec, n);

  if (sec >= 60.0) {
    sec = 0.0;
    min = min + 1.0;
  }

  if (min >= 60.0) {
    min = 0.0;
    deg = deg + 1.0;
  }

  if (((deg == 24.0) && m == 1) || (deg == 360.0)) {
    deg = 0.0;
  }

  strcpy(s, " ");
  if (x < 0) {
    strcpy(s, "-");
  }
  else if (m == 0) {
    strcpy(s, "+");
  }

  /* begin building up the return string in buffer */
  sprintf(buffer, "%.0f", (1000.0+deg));
  right(buffer, fig, right_buffer);
  strcat(s, right_buffer);
  left(dh, m+1, left_buffer);
  right(left_buffer, 1, right_buffer);
  strcat(s, right_buffer);

  sprintf(buffer, "%.0f", (100.0+min));
  right(buffer, 2, right_buffer);
  strcat(s, right_buffer);
  left(mm, m+1, left_buffer);
  right(left_buffer, 1, right_buffer);
  strcat(s, right_buffer);

  switch (n) {
    case 0:
      sprintf(buffer, "%.0f", (100.0+sec+0.0000001));
      break;
    case 1:
      sprintf(buffer, "%.1f", (100.0+sec+0.0000001));
      break;
    case 2:
      sprintf(buffer, "%.2f", (100.0+sec+0.0000001));
      break;
    case 3:
      sprintf(buffer, "%.3f", (100.0+sec+0.0000001));
      break;
    case 4:
      sprintf(buffer, "%.4f", (100.0+sec+0.0000001));
  }

  if (n == 0) {
    nf = 2;
  } else {
    nf = n + 3;
  }

  right(buffer, strlen(buffer)-1, right_buffer);
  left(right_buffer, nf, left_buffer);
  strcat(s, left_buffer);

  left(ss, m+1, left_buffer);
  right(left_buffer, 1, right_buffer);
  strcat(s, right_buffer);
}

/*********************************************************************
Name:    free_matrix
Purpose: Free a double matrix allocated by matrix().
Inputs:  m    - Pointer to matrix.
         nrow - Number of rows.
         ncol - Number of columns.
Outputs: None.
Returns: Nothing.
Status:  Finished.
Errors:  None known.
*********************************************************************/
void free_matrix(double **m, int nrow, int ncol) {

  int i;

  for (i = 0;i < nrow;i++) {
    free(m[i]);
  }

  free(m);
}

/*********************************************************************
Name:    freeDMatrix
Purpose: Free all storage associated with matrix m.
Inputs:  m - Matrix to free.
         n - Dimension of matrix.
Outputs: None.
Returns: Nothing.
Status:  Finished.
Errors:  None known.
*********************************************************************/
void freeDMatrix(DMatrix m, int n)
{
  int i;

  for (i=0; i<n; i++)
  {
    free(m[i]); /* storage for doubles in row i */
  }
  free(m); /* storage for pointers to rows */

  return;
}

/*********************************************************************
Name:    FunArgIAU
Purpose: Function to compute the fundamental arguments using the IAU
         expressions.
Inputs:  jed - JED on TDB scale.
Outputs: funarg[0]   - Mean anomaly of Moon.
         funarg[1]   - Mean anomaly of Sun.
         funarg[2]   - Argument of lat. of Moon.
         funarg[3]   - Mean elongation of Moon.
         funarg[4]   - Mean longitude of Moon's ascending node.
         funarg[5]   - Mean longitude of Moon.
         funarg[i+6] - Derivative of ith fundamental argument.
         NOTE: All derivatives are in units of rad/day.
Returns: Nothing.
Status:  Finished.
Errors:  None known.
*********************************************************************/
void FunArgIAU(double jed, double *funarg) {

  static double jedz;
  double T, T2, T3, L, Ldot, lp, lpdot, F, Fdot, D, Ddot;
  double N, Ndot, LL, LLdot;

  if (jed == jedz) return;
  jedz = jed;

  T = (jed - J2000) / JulCty;
  T2 = T * T;
  T3 = T2 * T;

  /* Compute fundamental arguments */
  L = 485866.733 + (1325. * 1296000. + 715922.633) * T
    + 31.31 * T2 + 0.064 * T3;
  Ldot = (1325. * 1296000. + 715922.633) + 2. * 31.31 * T
    + 3. * 0.064 * T2;
    L = amodulo(L * A2R, TWOPI);
    Ldot = amodulo(Ldot * A2R / 36525., TWOPI);
  lp = 1287099.804 + (99. * 1296000. + 1292581.224) * T
    - 0.577 * T2 - 0.012 * T3;
  lpdot = (99. * 1296000. + 1292581.224) - 2. * 0.577 * T
    - 3. * 0.012 * T2;
    lp = amodulo(lp * A2R, TWOPI);
    lpdot = amodulo(lpdot * A2R / 36525., TWOPI);
  F = 335778.877 + (1342. * 1296000. + 295263.137) * T
    - 13.257 * T2 + 0.011 * T3;
  Fdot = (1342. * 1296000. + 295263.137) - 2. * 13.257 * T
    + 3. * 0.011 * T2;
    F = amodulo(F * A2R, TWOPI);
    Fdot = amodulo(Fdot * A2R / 36525., TWOPI);
  D = 1072261.307 + (1236. * 1296000. + 1105601.328) * T
    - 6.891 * T2 + 0.019 * T3;
  Ddot = (1236. * 1296000. + 1105601.328) - 2. * 6.891 * T
    + 3. * 0.019 * T2;
    D = amodulo(D * A2R, TWOPI);
    Ddot = amodulo(Ddot * A2R / 36525., TWOPI);
  N = 450160.28 - (5. * 1296000. + 482890.539) * T
    + 7.455 * T2 + 0.008 * T3;
  Ndot = (5. * 1296000. + 482890.539) + 2. * 7.455 * T
    + 3. * 0.008 * T2;
    N = amodulo(N * A2R, TWOPI);
    Ndot = amodulo(Ndot * A2R / 36525., TWOPI);
  LL = 785939.157 + (1336. * 1296000. + 1108372.598) * T
    - 5.802 * T2 + 0.019 * T3;
  LLdot = (1336. * 1296000. + 1108372.598) - 2. * 5.802 * T
    + 3. * 0.019 * T2;
    LL = amodulo(LL * A2R, TWOPI);
    LLdot = amodulo(LLdot * A2R / 36525., TWOPI);

  funarg[0]  = L;
  funarg[6]  = Ldot;
  funarg[1]  = lp;
  funarg[7]  = lpdot;
  funarg[2]  = F;
  funarg[8]  = Fdot;
  funarg[3]  = D;
  funarg[9]  = Ddot;
  funarg[4]  = N;
  funarg[10] = Ndot;
  funarg[5]  = LL;
  funarg[11] = LLdot;
}

/*********************************************************************
Name:    GeocenObs
Purpose: Subprogram to compute the geocentric equatorial
         state vectors of an observer.  The vectors are
         referred to the J2000.0 frame.
Inputs:  JED - Julian Date on TDB scale.
Outputs: obsr_geo - Observer's geocentric state vector
                    in the celestial reference frame.
Returns: Nothing.
Status:  Finished.
Errors:  None known.
*********************************************************************/
void GeocenObs(double jed, double obsr_geo[]) {

  double gast,h,CosLat,SinLat,Cos2Lat,Sin2Lat,CosLon,SinLon,C,S;
  double PolMotX,PolMotY,g[6],dpsi,deps,dpsidot,depsdot;
  double zeta,z,theta,zetadot,zdot,thetadot;
  double MeanEps,TrueEps,MeanEpsDot,TrueEpsDot;

  /* Elevation is expected to be in meters. */
  h = (obsr_ele / 1000.0) * KM2AU;
  CosLat  = cos(obsr_lat);
  Cos2Lat = CosLat * CosLat;
  SinLat  = sin(obsr_lat);
  Sin2Lat = SinLat * SinLat;
  CosLon  = cos(obsr_lon);
  SinLon  = sin(obsr_lon);
  C = 1.0 / sqrt(Cos2Lat + (1.0 - EFlat) * (1.0 - EFlat) * Sin2Lat);
  S = (1.0 - EFlat) * (1.0 - EFlat) * C;

  /* Compute the observer's state vector in the Earth-fixed frame. */
  g[0] = (EarthRadAU * C + h) * CosLat * CosLon;
  g[1] = (EarthRadAU * C + h) * CosLat * SinLon;
  g[2] = (EarthRadAU * S + h) * SinLat;
  g[3] = 0.0;
  g[4] = 0.0;
  g[5] = 0.0;

  /* Compute the required angular quantities and their derivatives. */
  GetPrecessParams(J2000, jed, &zeta, &z, &theta,
    &zetadot, &zdot, &thetadot);
  GetDpsiDeps(jed, &dpsi, &deps, &dpsidot, &depsdot);
  GetGST(jed, 1, &gast);
  Obliquity(J2000, jed, 0, &MeanEps, &MeanEpsDot);
  Obliquity(J2000, jed, 1, &TrueEps, &TrueEpsDot);

  /* Somehow get the quantities PolMotX, PolMotY */
  PolMotX = 0.0;
  PolMotY = 0.0;

  /* Wobble */
  QRotate(g, 2, PolMotX, 0, 1, g);
  QRotate(g, 1, PolMotY, 0, 1, g);

  /* Spin */
  QRotate(g, 3, -gast, -EarAngVelRD, 1, g);

  /* Nutation */
  QRotate(g, 1,  TrueEps,  TrueEpsDot, 1, g);
  QRotate(g, 3,  dpsi,  dpsidot, 1, g);
  QRotate(g, 1, -MeanEps, -MeanEpsDot, 1, g);

  /* Precession */
  QRotate(g, 3,  zeta,  zetadot, 1, g);
  QRotate(g, 2, -theta, -thetadot, 1, g);
  QRotate(g, 3,  z, zdot, 1, obsr_geo);
}

/*********************************************************************
Name:    GetGST
Purpose: Function to compute the Greenwich mean or Greenwich apparent
         sidereal time, depending on which is required.  Strictly speaking,
         the computed time is the DYNAMICAL SIDEREAL TIME, but for all
         practical purposes is the same as the observed sidereal time.
Inputs:  jed - JED on TDB scale.
         s   - 1 for apparent sidereal time, otherwise mean sidereal time.
Outputs: gst - GMST or GAST in radians.
Returns: Nothing.
Status:  Finished.
Errors:  None known.
*********************************************************************/
void GetGST(double jed, int s, double *gst) {

  static double jedz, sz;
  double T, dpsi, deps, dpsidot, depsdot;
  double meaneps, trueeps, meanepsdot, trueepsdot, eqeqnx, eqeqnxdot;

  if ((jed == jedz) && (s == sz)) return;
  jedz = jed;
  sz = s;

  T = (jed - J2000) / JulCty;

  /* compute GMST in seconds */
  *gst = 67310.54841 + T * ((876600.0 * 3600.0 + 8640184.812866)
    + T * (0.093104 + T * (-0.0000062)));

  /* convert to radians */
  *gst = amodulo(*gst / 3600.0, 24.0) * H2R;

  if (s == 1) {
    /* get nutation quantities */
    GetDpsiDeps(jed, &dpsi, &deps, &dpsidot, &depsdot);

    /* get mean obliquity */
    Obliquity(J2000, jed, 0, &meaneps, &meanepsdot);

    /* get true obliquity */
    Obliquity (J2000, jed, 1, &trueeps, &trueepsdot);

    /* get equation of the equinoxes */
    eqeqnx    = dpsi * cos(trueeps);
    eqeqnxdot = dpsidot * cos(trueeps) - dpsi * trueepsdot * sin(trueeps);
    *gst     += eqeqnx;
  }

  *gst = amodulo(*gst, TWOPI);
}

/*********************************************************************
Name:    GetInvQMatrix
Purpose: Function to compute the elements of the inverse of a
         given Q-matrix.
Inputs:  QMatrix - Zero-offset 6X6 Q-matrix.
Outputs: InvQMatrix - Zero-offset inverse of QMatrix.
Returns: Nothing.
Status:  Finished.
Errors:  None known.
*********************************************************************/
void GetInvQMatrix(DMatrix QMatrix, DMatrix InvQMatrix) {

  int i, j;

  for (i= 0; i<3; i++) {
    for (j=0; j<3; j++) {
      InvQMatrix[i][j] = QMatrix[j][i];
      InvQMatrix[i][j+3] = 0.0;
      InvQMatrix[i+3][j] = QMatrix[j+3][i];
      InvQMatrix[i+3][j+3] = QMatrix[j+3][i+3];
    }
  }
}

/*********************************************************************
Name:    GetPrecessParams
Purpose: Function that computes the general precession parameters and their
         derivatives for precessing equatorial rectangular coordinates and
         velocities from jed1 to jed2. Lieske's formulae are used.
Inputs:  jed1 - Initial Julian Date.
         jed2 - Final Julian Date.
Outputs: zeta     - equatorial precession parameter.
         z        - equatorial precession parameter.
         theta    - equatorial precession parameter.
         zetadot  - derivative in rad/day.
         zdot     - derivative in rad/day.
         thetadot - derivative in rad/day.
Returns: Nothing.
Status:  Finished.
Errors:  None known.
*********************************************************************/
void GetPrecessParams(double jed1, double jed2, double *zeta, double *z,
  double *theta, double *zetadot, double *zdot, double *thetadot) {

  double T1, T2, c1, c2, c3, c4, c5, c6, p1, p2, x, xdot;

  T1 = (jed1 - J2000) / JulCty;
  T2 = (jed2 - jed1) / JulCty;

  /* compute zeta, z, theta, zetadot, zdot, thetadot */
  c1 = 2306.2181;
  c2 =    1.39656;
  c3 =   -0.000139;
  c4 =    0.30188;
  c5 =   -0.000344;
  c6 =    0.017998;
  p1 = c1 + c2 * T1 + c3 * T1 * T1;
  p2 = c4 + c5 * T1;
  x = p1 * T2 + p2 * T2 * T2 + c6 * T2 * T2 * T2;
  xdot = p1 + 2.0 * p2 * T2 + 3.0 * c6 * T2 * T2;
  *zeta = x * A2R;
  *zetadot = xdot * A2R / JulCty;

  c1 = 2306.2181;
  c2 =    1.39656;
  c3 =   -0.000139;
  c4 =    1.09468;
  c5 =    0.000066;
  c6 =    0.018203;
  p1 = c1 + c2 * T1 + c3 * T1 * T1;
  p2 = c4 + c5 * T1;
  x = p1 * T2 + p2 * T2 * T2 + c6 * T2 * T2 * T2;
  xdot = p1 + 2.0 * p2 * T2 + 3.0 * c6 * T2 * T2;
  *z = x * A2R;
  *zdot = xdot * A2R / JulCty;

  c1 = 2004.3109;
  c2 =   -0.85330;
  c3 =   -0.000217;
  c4 =   -0.42665;
  c5 =   -0.000217;
  c6 =   -0.041833;
  p1 = c1 + c2 * T1 + c3 * T1 * T1;
  p2 = c4 + c5 * T1;
  x = p1 * T2 + p2 * T2 * T2 + c6 * T2 * T2 * T2;
  xdot = p1 + 2.0 * p2 * T2 + 3.0 * c6 * T2 * T2;
  *theta = x * A2R;
  *thetadot = xdot * A2R / JulCty;
}

/*********************************************************************
Name:    GetQMatrix
Purpose: Function to compute the elements of the Q-matrix for a
         given angle and its time derivative.
Inputs:  phi    - Angle in radians
         phidot - Time der. of phi.
         axis   - 1 for x-axis,
                  2 for y-axis,
                  3 for z-axis.
         s      - 0 for inertial to inertial,
                  1 for inertial to rotating.
Outputs: QMatrix - Zero-offset 6X6 Q-matrix.
Returns: Nothing.
Status:  Finished.
Errors:  None known.
*********************************************************************/
void GetQMatrix(double phi, double phidot, int axis,
  int s, DMatrix QMatrix) {

  DMatrix r, dr;
  int i, j;

  r  = createDMatrix(3);
  dr = createDMatrix(3);

  /* form the 3X3 r[] and dr[] matrices */
  RotMat(axis, phi, r, dr);

  /* form the 6X6 Q-matrix */
  for (i=0; i<3; i++) {
    for (j=0; j<3; j++) {
      QMatrix[i][j] = r[i][j];
      QMatrix[i][j+3] = 0.0;
      QMatrix[i+3][j] = phidot * dr[i][j] * s;
      QMatrix[i+3][j+3] = r[i][j];
    }
  }

  freeDMatrix(r, 3);
  freeDMatrix(dr, 3);
}

/*********************************************************************
Name:    GetRPNmat
Purpose: Function to compute the precession matrix, the nutation
         matrix, or the combined R matrix for equatorial coordinates.
Inputs:  jed1 - Initial JED on TDB scale.
         jed2 - Final JED on TDB scale.
         rpn  - 1 for precession matrix,
                2 for nutation matrix,
                3 for R matrix.
         d    - 1 for zero-offset 3X3 matrix,
                2 for zero-offset 6X6 Qmatrix.
Outputs: m3 - Requested zero-offset 3X3 matrix.
         m6 - Requested zero-offset 6X6 matrix.
Returns: Nothing.
Status:  Finished.
Errors:  None known.
*********************************************************************/
void GetRPNmat(double jed1, double jed2, int rpn, int d,
  DMatrix m3, DMatrix m6) {

  double dpsi, deps, dpsidot, depsdot, trueeps, trueepsdot,
    meaneps, meanepsdot, zeta, zetadot, z, zdot, theta, thetadot;
  DMatrix p1mat, p2mat, p3mat, pmat, n1mat, n2mat, n3mat,
    nmat, rmat, p, n;

  p1mat = createDMatrix(6);
  p2mat = createDMatrix(6);
  p3mat = createDMatrix(6);
  pmat  = createDMatrix(6);
  n1mat = createDMatrix(6);
  n2mat = createDMatrix(6);
  n3mat = createDMatrix(6);
  nmat  = createDMatrix(6);
  rmat  = createDMatrix(6);
  p     = createDMatrix(6);
  n     = createDMatrix(6);

  switch (rpn) {
    case 1:
      /* compute precession matrix */
      GetPrecessParams(jed1,jed2,&zeta,&z,&theta,&zetadot,&zdot,&thetadot);
      GetQMatrix(-zeta, -zetadot, 3, 0, p1mat);
      GetQMatrix(theta, thetadot, 2, 0, p2mat);
      GetQMatrix(-z, -zdot, 3, 0, p3mat);
      MatXMat(p2mat, p1mat, pmat, 6);
      MatXMat(p3mat, pmat, m3, 6);
      freeDMatrix(p3mat, 6);
      freeDMatrix(p2mat, 6);
      freeDMatrix(p1mat, 6);
      freeDMatrix(pmat, 6);
      break;
    case 2:
      /* compute nutation matrix */
      GetDpsiDeps(jed2, &dpsi, &deps, &dpsidot, &depsdot);
      Obliquity(jed1, jed2, 0, &meaneps, &meanepsdot);
      Obliquity(jed1, jed2, 1, &trueeps, &trueepsdot);
      GetQMatrix(meaneps, meanepsdot, 1, 0, n1mat);
      GetQMatrix(-dpsi, -dpsidot, 3, 0, n2mat);
      GetQMatrix(-trueeps, -trueepsdot, 1, 0, n3mat);
      MatXMat(n2mat, n1mat, nmat, 6);
      MatXMat(n3mat, nmat, m3, 6);
      freeDMatrix(n3mat, 6);
      freeDMatrix(n2mat, 6);
      freeDMatrix(n1mat, 6);
      freeDMatrix(nmat, 6);
      break;
    case 3:
      /* compute R matrix */
      GetPrecessParams(jed1,jed2,&zeta,&z,&theta,&zetadot,&zdot,&thetadot);
      GetQMatrix( -zeta, -zetadot, 3, 0, p1mat);
      GetQMatrix(theta, thetadot, 2, 0, p2mat);
      GetQMatrix(-z, -zdot, 3, 0, p3mat);
      MatXMat(p2mat, p1mat, pmat, 6);
      MatXMat(p3mat, pmat, p, 6);
      GetDpsiDeps(jed2, &dpsi, &deps, &dpsidot, &depsdot);
      Obliquity(jed1, jed2, 0, &meaneps, &meanepsdot);
      Obliquity(jed1, jed2, 1, &trueeps, &trueepsdot);
      GetQMatrix(meaneps, meanepsdot, 1, 0, n1mat);
      GetQMatrix(-dpsi, -dpsidot, 3, 0, n2mat);
      GetQMatrix(-trueeps, -trueepsdot, 1, 0, n3mat);
      MatXMat(n2mat, n1mat, nmat, 6);
      MatXMat(n3mat, nmat, n, 6);
      MatXMat(n, p, m3, 6);
      freeDMatrix(p3mat, 6);
      freeDMatrix(p2mat, 6);
      freeDMatrix(p1mat, 6);
      freeDMatrix(pmat, 6);
      freeDMatrix(p, 6);
      freeDMatrix(n3mat, 6);
      freeDMatrix(n2mat, 6);
      freeDMatrix(n1mat, 6);
      freeDMatrix(nmat, 6);
      freeDMatrix(n, 6);
      break;
  }
}

/*********************************************************************
Name:    GetStateVector
Purpose: Function to retrieve the state vector of a given body
         w.r.t. a given origin at a given time. This routine is
         a wrapper for PLEPH.
Inputs:  jd     - Julian date on TDB time scale.
         targ   - Body identification number.
                  1=MER,  2=VEN,  3=EAR,  4=MAR
                  5=JUP,  6=SAT,  7=URA,  8=NEP
                  9=PLU,  10=MOO, 11=SUN, 12=SSB
                  13=EMB, 14=NUT, 15=LIB
         cent   - Origin identification number.
                  (Numbering as above).
         recpol - 1 for cartesian, 2 for orthogonal polar.
Outputs: StateVector - A 15x15x2x6 array containing the requested
                       state vector.
         NOTE: StateVector has the following form:
         StateVector[targ-1][cent-1][recpol-1][component]
Returns: Nothing.
Status:  Finished.
Errors:  None known.
*********************************************************************/
void GetStateVector(double jd, int targ, int cent, int recpol,
  double StateVector[15][15][2][6]) {

  static double b[6] = {0.};
  static double rrd[6] = {0.};
  static double dpsi, deps, dpsidot, depsdot;
  int inside, i;

  /* Get the state vector from JPL ephemeris */
  pleph(jd, targ, cent, rrd, &inside);

  if (!inside) {
    LogMsg(stderr,
      "GetStateVector: requested date not covered by ephemeris file.\n");
    for (i=0; i<6; i++) {
      StateVector[targ-1][cent-1][recpol-1][i] = 999.;
    }
    exit (1);
  }

  /* Check for nutations */
  if (targ == 14) {
    dpsi    = rrd[0];
    deps    = rrd[1];
    dpsidot = rrd[2];
    depsdot = rrd[3];
    recpol  = 1;
  }

  if (recpol == 1) {
    for (i=0; i<6; i++) {
      StateVector[targ-1][cent-1][recpol-1][i] = rrd[i];
    }
  } else {
    /* Convert rectangular vector to polar vector */
    Rec2Pol(rrd, b);
    for (i=0; i<6; i++) {
      StateVector[targ-1][cent-1][recpol-1][i] = b[i];
    }
  }

  return;
}

/*********************************************************************
Name:    HelEphemeris
Purpose: Function to compute an orbiting body's heliocentric
         ecliptic state vector for a single instant of time
         given the universal orbital elements and the time.
         Reference: Mansfield. 1986. AIAA Paper 86-2269-CP.
Inputs:  uelement[] - Array of universal orbital elements:
         uelement[0] = q
         uelement[1] = e
         uelement[2] = i
         uelement[3] = node
         uelement[4] = arg. peri.
         uelement[5] = T
         mu         - Gravitational constant for the body.
         jed        - Time.
Outputs: posvel[] - State vector:
         posvel[0..2] = position vector.
         posvel[3..5] = velocity vector.
Returns: Nothing.
Status:  Finished.
Errors:  None known.
*********************************************************************/
void HelEphemeris(double *uelement, double mu, double jed, double *posvel) {

  int i;
  double LongPeri, RetLongPeri, cosi, sini, coslp, sinlp, cosrlp, sinrlp;
  double cosw, sinw, rcosnu, rsinnu, r, cosnu, sinnu, param, p[3], q[3];

  /* Compute longitude of perihelion */
  LongPeri = uelement[4] + uelement[3];
  /* Compute retrograde longitude of perihelion */
  RetLongPeri = uelement[4] - uelement[3];

  /* Compute the P vector */
  cosi   = cos(uelement[2]);
  sini   = sin(uelement[2]);
  coslp  = cos(LongPeri);
  sinlp  = sin(LongPeri);
  cosrlp = cos(RetLongPeri);
  sinrlp = sin(RetLongPeri);
  cosw   = cos(uelement[4]);
  sinw   = sin(uelement[4]);
  p[0] = 0.5 * (1.0 + cosi) * coslp + 0.5 * (1.0 - cosi) * cosrlp;
  p[1] = 0.5 * (1.0 + cosi) * sinlp - 0.5 * (1.0 - cosi) * sinrlp;
  p[2] = sinw * sini;

  /* Compute the Q vector */
  q[0] = -0.5 * (1.0 + cosi) * sinlp - 0.5 * (1.0 - cosi) * sinrlp;
  q[1] =  0.5 * (1.0 + cosi) * coslp - 0.5 * (1.0 - cosi) * cosrlp;
  q[2] = cosw * sini;

  /* Solve Kepler's equation */
  Conway(uelement, mu, jed, &rcosnu, &rsinnu);

  /* Compute magnitude of radius vector */
  r = sqrt(rcosnu * rcosnu + rsinnu * rsinnu);
  cosnu = rcosnu / r;
  sinnu = rsinnu / r;

  /* Compute heliocentric ecliptic position vector */
  for (i = 0; i < 3; i++) {
    posvel[i] = p[i] * rcosnu + q[i] * rsinnu;
  }

  /* Compute heliocentric ecliptic velocity vector */
  param = uelement[0] * (1.0 + uelement[1]);

  for (i = 3; i < 6; i++) {
    posvel[i]   = -p[i-3] * sqrt(mu / param) * sinnu;
    posvel[i] = posvel[i] + q[i-3] * sqrt(mu / param) *
      (uelement[1] + cosnu);
  }
}

/*********************************************************************
Name:    interp
Purpose: This subroutine differentiates and interpolates a set of
         Chebyshev coefficients to give position and velocity.
Inputs:  buff - 1st location of array of Chebyshev coefficients.
         t    - t[0] is fractional time in interval covered by
                coefficients at which interpolation is wanted
                (0 <= t[0] <= 1). t[1] is length of whole
                interval in input time units.
         ncf  - Number of coefficients per component.
         ncm  - Number of components per set of coefficients.
         na   - Number of sets of coefficients in full array
                (i.e., # of sub-intervals in full interval).
         fl   - Integer flag:  =1 for positions only.
                               =2 for pos and vel.
Outputs: dumpv - Interpolated quantities requested. Dimension
         expected is pv[ncm][fl].
Returns: Nothing.
Status:  Finished.
Errors:  None known.
*********************************************************************/
void interp(int buff, double *t, int ncf, int ncm, int na,
  int fl, double dumpv[3][2]) {

  int i, j, l, n, m;
  static int bcoef, np, nv;
  static double pc[18], vc[18], cbody[1200];
  static double cbuf[15][3][8]={0.};
  static double twot, dna, dt1, temp, ll, tc, vfac;

  np = 2;
  nv = 3;
  twot = 0.0;
  pc[0] = 1.0;
  pc[1] = 0.0;
  vc[1] = 1.0;

  /*
      Entry point. Get correct sub-interval number for this set
      of coefficients and then get normalized chebyshev time
      within that subinterval.
  */

  dna = (double) na;
  dt1 = floor(t[0]);
  temp = dna*t[0];
  ll = floor((temp - dt1) + 1.0);

  /* 'tc' is the normalized chebyshev time (-1 <= tc <= 1) */
  tc = 2.0*(fmod(temp, 1.0) + dt1) - 1.0;

  /*
      Check to see whether chebyshev time has changed,
      and compute new polynomial values if it has.
      (The element pc[1] is the value of 'tc' and hence
      contains the value of 'tc' on the previous call.)
  */

  if (tc != pc[1]) {
    np = 2;
    nv = 3;
    pc[1] = tc;
    twot = tc + tc;
  }

  /*
      Be sure that at least 'ncf' polynomials have been evaluated
      and are stored in the array pc[].
  */

  if (np < ncf) {
    for (i=np; i<ncf; i++) {
      pc[i] = twot * pc[i-1] - pc[i-2];
    }
    np = ncf;
  }

  /* interpolate to get position for each component */

  /* number of coefficients for body */
  bcoef= ncf*na*ncm;

  /* stored body's coefficients in an array */
  n = buff;
  for (m=0; m<bcoef; m++) {
    cbody[m] = db[n];
    n++;
  }

  /* fill the cbuf[][][] array */
  n = 0;
  /* loop for each sub-interval */
  for (l=0; l<na; l++) {
    /* loop for each component */
    for (i=0; i<ncm; i++) {
      /* loop for each set of coeffs */
      for (j=0; j<ncf; j++) {
        cbuf[j][i][l] = cbody[n];
        n++;
      }
    }
  }

  for (i=0; i<ncm; i++) {
    dumpv[i][0] = 0.0;
    for (j=ncf-1; j>=0; j--) {
      dumpv[i][0] = dumpv[i][0] + pc[j] * cbuf[j][i][((int)ll)-1];
    }
  }
  if (fl <= 1) return;

  /*
      If velocity interpolation is wanted, be sure enough
      derivative polynomials have been generated and stored.
  */

  vfac = (dna+dna)/t[1];
  vc[2] = twot+twot;
  if (nv < ncf) {
    for (i=nv; i<ncf; i++) {
      vc[i]=twot * vc[i-1] + pc[i-1] + pc[i-1] - vc[i-2];
    }
    nv = ncf;
  }

  /* interpolate to get velocity for each component */

  for (i=0; i<ncm; i++) {
    dumpv[i][1] = 0.0;
    for (j=ncf-1; j>=1; j--) {
      dumpv[i][1] = dumpv[i][1] + vc[j] * cbuf[j][i][((int)ll)-1];
    }
    dumpv[i][1] = dumpv[i][1] * vfac;
  }
}

/*********************************************************************
Name:    Interp1
Purpose: Subprogram to perform interpolation on a list
         of tabular values to find maxima/minima or
         to find the zero of a function within the limits
         of the table.  The algorithms used are those of Meeus.
Inputs:  x[] - Array of function arguments.
         y[] - Array of function values at the
               corresponding arguments.
         L   - Equals 3 for 3-pt. interpolation.
               5 for 5-pt. interpolation.
         m   - Equals 0 for zero of function.
               1 for a extremum of the function.
Outputs: arg - Argument corresponding to the
               extremum or zero of the tabular function.
         v   - The value of the function corresponding to arg.
Returns: Nothing.
Status:  Finished.
Errors:  None known.
*********************************************************************/
void Interp1(double *x, double *y, int L, int m, double *arg, double *v) {

  double a, b, c, d, e, f, g, h, j, k, nm, n0, n1;

  if (L == 3) {
    /* Use a 3 pt. table */
    a = y[1] - y[0];
    b = y[2] - y[1];
    c = b - a;
    if (m == 1) {
      /* Find extremum of the tabular function */
      *v = y[1] - (a + b) * (a + b) / (8.0 * c);
      nm = -(a + b) / (2.0 * c);
      *arg = x[1] + nm * (x[1] - x[0]);
    } else {
      n1 = 0.0; /* Find zero of the tabular function */
      do {
        n0 = -2.0 * y[1] / (a + b + c * n1);
        if (fabs(n0 - n1) < 0.000000000000001)
          break;
        n1 = n0;
      } while (1 == 1);
      *v = y[1] + 0.5 * n0 * (x[1] - x[0]) * (a + b + n0 * c);
      *arg = x[1] + n0 * (x[1] - x[0]);
    }
  } else {
    /* Use a 5 pt. table */
    a = y[1] - y[0];
    b = y[2] - y[1];
    c = y[3] - y[2];
    d = y[4] - y[3];
    e = b - a;
    f = c - b;
    g = d - c;
    h = f - e;
    j = g - f;
    k = j - h;
    if (m == 0) {
      /* Find extremum of tabular function */
      n1 = 0.0;
      do {
        nm = (6.0 * b + 6.0 * c - h - j + 3.0 * n1 * n1 *
          (h + j) + 2.0 * n1 * n1 * n1 * k) / (k - 12.0 * f);
        if (fabs(nm - n1) < 0.000000000000001)
          break;
        n1 = nm;
      } while (1 == 1);
      nm = nm * (x[3] - x[2]);
      *arg = x[2] + nm;
      *v = y[2] + 0.5 * nm * (b + c) + 0.5 * nm * nm * f + nm
        * (nm * nm - 1.0) * (h + j) / 12.0 + nm * nm * (nm * nm - 1.0)
        * k / 24.0;
    } else {
      /* Find zero of the tabular function */
      n1 = 0.0;
      do {
        n0 = (-24.0 * y[2] + n1 * n1 * (k - 12.0 * f) - 2.0 *
          n1 * n1 * n1 * (h + j) - n1 * n1 * n1 * n1 * k) / (2.0
          * (6.0 * b + 6.0 * c - h - j));
        if (fabs(n0 - n1) < 0.000000000000001)
          break;
        n1 = n0;
      } while (1 == 1);
      n0 = n0 * (x[3] - x[2]);
      *arg = x[2] + n0;
      *v = y[2] + 0.5 * n0 * (b + c) + 0.5 * n0 * n0 * f + n0 * (n0 * n0 - 1.0) * (h + j) / 12.0 + n0 * n0 * (n0 * n0 - 1.0) * k / 24.0;
    }
  }
}

/*********************************************************************
Name:    Interpol
Purpose: Function for 3 or 5 point interpolation
         using Meeus' algorithm.
Inputs:  x[] - Aarray of function arguments.
         y[] - Array of function values at
               the three arguments.
         i   - Switch indicating whether a
               3 or 5 point interpolation
               is performed (i=3 for 3 pt.,
               i=5 for 5 pt.).
         arg - Argument for which a function value
               is needed.
Outputs: None.
Returns: Interpolated value.
Status:  Finished.
Errors:  None known.
*********************************************************************/
double Interpol(double *x, double *y, int i, double arg) {

  double d1, d2, d3, d4, d5, d6, d7, d8, d9, d10, n;

  if (i == 3) {
    d1 = y[1] - y[0];
    d2 = y[2] - y[1];
    d3 = d2 - d1;
    n = (arg - x[1]) / (x[1] - x[0]);
    return (y[1] + (n / 2.0) * (d1 + d2 + n * d3));
  } else {
    d1 = y[1] - y[0];
    d2 = y[2] - y[1];
    d3 = y[3] - y[2];
    d4 = y[4] - y[3];
    d5 = d2 - d1;
    d6 = d3 - d2;
    d7 = d4 - d3;
    d8 = d6 - d5;
    d9 = d7 - d6;
    d10 = d9 - d8;
    n = (arg - x[2]) / (x[2] - x[1]);
    return (y[2] + n * ((d2 + d3) / 2.0 - (d8 + d9) / 12.0) +
      n * n * (d6 / 2.0 - d10 / 24.0) + n * n * n *
      ((d8 + d9) / 12.0) + n * n * n * n * (d10 / 24.0));
  }
}

/*********************************************************************
Name:    JED2Cal
Purpose: Function to convert a JED to an ordinary calendar date.
         The algorithm is that of Meeus.
Inputs:  jed - Julian Ephemeris Date.
Outputs: yr - Year.
         mo - Month.
         dy - Day.
         ti - Time of day in decimal hours.
Returns: Nothing.
Status:  Finished.
Errors:  None known.
*********************************************************************/
void JED2Cal(double jed, int *yr, int *mo, int *dy, double *ti) {

  static double jedz;
  double z, f, a, b, c, d, e, alpha, daywtime;

  if (jed == jedz) return;
  jedz = jed;

  z = floor(jed + 0.5);
  f = (jed + 0.5) - z + 5.0e-10; /* fudge factor */
  if (z < 2299161.0) {
    a = z;
  } else {
    alpha = floor((z - 1867216.25) / 36524.25);
    a = z + 1.0 + alpha - floor(0.25 * alpha);
  }

  b = a + 1524.0;
  c = floor((b - 122.1) / 365.25);
  d = floor(365.25 * c);
  e = floor((b - d) / 30.6001);
  daywtime = b - d - floor(30.6001 * e) + f;
  *dy = (int) floor(daywtime);
  *ti = 24.0 * (daywtime - (double) *dy);

  if (e  < 13.5) *mo = (int) floor(e - 1.0);
  if (e  > 13.5) *mo = (int) floor(e - 13.0);
  if (*mo > 2.5) *yr = (int) floor(c - 4716.0);
  if (*mo < 2.5) *yr = (int) floor(c - 4715.0);
}

/*********************************************************************
Name:    JED2Epoch
Purpose: Function to convert a given JED to its corresponding Julian
         or Besselian epoch. For example, 2451545.0 becomes J2000.
Inputs:  jed - Input jed.
         s   - J or B, whichever is desired.
Outputs: epoch - J or B epoch as a string.
Returns: Nothing.
Status:  Finished.
Errors:  None known.
*********************************************************************/
void JED2Epoch(double jed, char *s, char *epoch) {

  char d[80];
  double date;

  if (strcmp(s,"J") == 0) {
    date = 2000.0 + (jed - 2451545.0) / 365.25;
    sprintf(d, "%f", date);
  } else {
    date = 1900.0 + (jed - (double) 2415020.31352) /
      (double) 365.242198781731;
    sprintf(d, "%f", date);
  }

  strcpy(epoch, s);
  strcat(epoch, d);
}

/*********************************************************************
Name:    LightTime
Purpose: Subprogram to compute the vectors body_geo() and body_hel()
         from a barycentric ephemeris of the object, Earth, and the Sun.
Inputs:  JED       - Desired time of reduction.
         Body      - Body number, 99 for stellar reduction.
         earth_ssb - Barycentric state vector of Earth at JED.
         earth_hel - Heliocentric state vector of Earth.
         StarData  - Stellar data.
Outputs: body_geo - Geometric geocentric state vector of object.
         body_hel - Heliocentric state vector of object.
         Ltime    - Light time in days to the object.
Returns: Nothing.
Status:  Finished.
Errors:  None known.
*********************************************************************/
void LightTime(double jed, int body, double earth_ssb[],
  double earth_hel[], double StarData[], double body_geo[],
  double body_hel[], double *Ltime) {

  int inside, i;
  double E[3],Edot[3],RQB[6],RSB[6],RQ[6],RP[6],Q[3],Qdot[3],P[3],Pdot[3];
  double unit_body_hel[6], unit_body_geo[6];
  double magE,magQ,magP,trialtau,tau,xtime,RelTerm,
         TrueHelDist=0,TrueGeoDist=0;
  double AppHelDist,AppGeoDist,RA,DE,sinRA,cosRA,sinDE,cosDE;
  double plx,r,muRA,muDE,rdot,drda[3],drdd[3],DT;

  SplitStateVector(earth_hel, E, Edot);
  magE = Vecmag(E);

  if (body != 99) {
    trialtau = 0.0;
    do {
      xtime = jed - trialtau;

      /* body's barycentric state vector */
      pleph(xtime, body, 12, RQB, &inside);

      /* Sun's barycentric state vector */
      pleph(xtime, 11, 12, RSB, &inside);

      /* body's heliocentric state vector */
      for (i=0; i<6; i++) {
        RQ[i] = RQB[i] - RSB[i];
        /* body's geocentric state vector */
        RP[i] = RQB[i] - earth_ssb[i];
      }

      /* Calculate improved value of the light time */
      SplitStateVector(RQ, Q, Qdot);
      SplitStateVector(RP, P, Pdot);
      magQ = Vecmag(Q);
      magP = Vecmag(P);

      if (trialtau == 0) {
        TrueHelDist = magQ;
        TrueGeoDist = magP;
      }

      if (magQ == 0) {
        RelTerm = 0.0;
      } else {
        RelTerm = MUC * log((magE + magP + magQ) / fabs(magE - magP +
                             magQ));
      }
      tau = (magP + RelTerm) / CAUD;
      if (fabs(trialtau - tau) < 1e-9) {
        break;
      }
      trialtau = tau;
    } while(1);
  } else {
    RA = StarData[0];
    DE = StarData[1];
    sinRA = sin(RA);
    cosRA = cos(RA);
    sinDE = sin(DE);
    cosDE = cos(DE);
    if (StarData[2] == 0) {
      plx = 1e-7;
    } else {
      plx = StarData[2];
    }

    /* Star's distance in AU */
    r = R2A / plx;

    /* Star's barycentric position vector */
    RQB[0] = r*cosRA*cosDE;
    RQB[1] = r*sinRA*cosDE;
    RQB[2] = r      *sinDE;

    /* Convert proper motions and radial velocity */
    /* to units of AU/day */
    muRA = StarData[3] * 15.0 * cosDE / (JulCty * plx);
    muDE = StarData[4] / (JulCty * plx);
    rdot = StarData[5] * 86400.0 * KM2AU;

    /* Partial of unit vector wrt RA */
    drda[0] = - sinRA * cosDE;
    drda[1] =   cosRA * cosDE;
    drda[2] =   0.0;

    /* Partial of unit vector wrt DE */
    drdd[0] = - cosRA * sinDE;
    drdd[1] = - sinRA * sinDE;
    drdd[2] =           cosDE;

    /* Star's barycentric velocity vector */
    for (i=3; i<6; i++) {
      RQB[i] = muRA*drda[i-3] + muDE*drdd[i-3] + rdot*RQB[i-3]/r;
    }

    /* Correction for space motion */
    DT = jed - J2000;
    for (i=0; i<3; i++) {
      RQB[i] = RQB[i] + RQB[i+3]*DT;
    }

    /* Sun's barycentric state vector */
    pleph(jed, 11, 12, RSB, &inside);

    /* Star's heliocentric state vector */
    for (i=0; i<6; i++) {
      RQ[i] = RQB[i] - RSB[i];
    }

    /* Star's geo/topocentric state vector */
    /* This correction is annual parallax. */
    for (i=0; i<6; i++) {
      RP[i] = RQB[i] - earth_ssb[i];
    }

    SplitStateVector(RQ, Q, Qdot);
    SplitStateVector(RP, P, Pdot);
    magQ = Vecmag(Q);
    magP = Vecmag(P);

    TrueGeoDist = magP;
    TrueHelDist = magQ;

    tau = r / CAUD;
  }

  *Ltime = tau;
  AppHelDist = magQ;
  AppGeoDist = magP;

  for (i=0; i<6; i++) {
    body_hel[i] = RQ[i];
    body_geo[i] = RP[i];
  }

  Uvector(body_hel, unit_body_hel);
  Uvector(body_geo, unit_body_geo);

  /* I don't understand why this next loop is needed, but it is. */
  for (i=3; i<6; i++) {
   unit_body_hel[i] = RQ[i];
   unit_body_geo[i] = RP[i];
  }

  for (i=0; i<6; i++) {
    body_geo[i] = unit_body_geo[i] * TrueGeoDist;
    body_hel[i] = unit_body_hel[i] * TrueHelDist;
  }
}

/*********************************************************************
Name:    matrix
Purpose: Function that allocates a double matrix with dimensions
         m[nrow][ncol].
Inputs:  nrow - Number of rows.
         ncol - Number of columns.
Outputs: None.
Returns: Pointer to 2-D matrix of doubles.
Status:  Finished.
Errors:  None known.
*********************************************************************/
double **matrix(int nrow, int ncol) {

  long i;
  double **m;

  /* Allocate pointer to rows */
  m = (double **) malloc(nrow * sizeof(double *));
  if (!m) {
    LogMsg(stderr, "matrix(): allocation failure 1\n");
    exit(1);
  }

  /* Allocate rows and set pointers to them */
  for (i = 0;i < nrow;i++) {
    m[i] = (double *) malloc(ncol * sizeof(double));
    if (!m[i]) {
      LogMsg(stderr, "matrix(): allocation failure 2\n");
      exit(1);
    }
  }
  /* return pointer to array of pointers to rows */
  return m;
}

/*********************************************************************
Name:    MatXMat
Purpose: Square matrix multiplication subroutine.
         NOTE: If a[][] and c[][] have the same name in the
         calling program, the original a[][] matrix is overwritten.
Inputs:  a[][] - Zero-offset matrix a.
         b[][] - Zero-offset matrix b.
         n     - Dimension of matrix.
Outputs: c[][] - Zero-offset matrix a x b.
Returns: Nothing.
Status:  Finished.
Errors:  None known.
*********************************************************************/
void MatXMat(DMatrix a, DMatrix b, DMatrix c, int n) {

  int i, j, k;

  for (i=0; i<n; i++) {
    for (j=0; j<n; j++) {
      c[i][j] = 0.0;
      for (k=0; k<n; k++) {
        c[i][j] += a[i][k] * b[k][j];
      }
    }
  }
}

/*********************************************************************
Name:    MatXVec
Purpose: Matrix/vector multiplication subroutine.
Inputs:  a[][] - Zero-offset matrix a (l rows by m columns).
         b[]   - Zero-offset vector b (m rows by 1 column).
         l     - Dimension of vector (i.e. lx1).
         m     - Dimension of matrix (i.e. mxl).
Outputs: c[]   - Zero-offset vector a x b (l rows by 1 column).
Returns: Nothing.
Status:  Finished.
Errors:  None known.
*********************************************************************/
void MatXVec (DMatrix a, double b[], double c[], int l, int m) {

  int i, j;
  double s, temp[6];

  for (i=0; i<l; i++) {
    s = 0.0;
    for (j=0; j<m; j++) {
      s += a[i][j] * b[j];
    }
    temp[i] = s;
  }

  for (i=0; i<l; i++) {
    c[i] = temp[i];
  }
}

/*********************************************************************
Name:    MRotate
Purpose: Function to perform a matrix rotation of an input vector
         through a given angle about a desired axis. A right-handed
         orthogonal coordinate system is assumed, and a 3X3
         rotation matrix is used, not a Q-matrix.
         NOTE: The vin[] and vout[] can have the same name in the
         calling program.
Inputs:  vin[] - Zero-offset input vector.
         axis  - 1 for rot. about x-axis,
                 2 for rot. about y-axis,
                 3 for rot. about z-axis.
         phi   - Rotation angle in radians.
Outputs: vout[] - Zero-offset transformed vector.
Returns: Nothing.
Status:  Finished.
Errors:  None known.
*********************************************************************/
void MRotate(double vin[], int axis, double phi, double vout[]) {

  double cosphi, sinphi, T;

  cosphi = cos(phi);
  sinphi = sin(phi);

  switch (axis) {
    case 1:  /* rotate about x-axis */
      T = cosphi * vin[1] + sinphi * vin[2];
      vout[2] = -sinphi * vin[1] + cosphi * vin[2];
      vout[1] = T;
      vout[0] = vin[0];
      break;
    case 2:  /* rotate about y-axis */
      T = cosphi * vin[0] - sinphi * vin[2];
      vout[2] = sinphi * vin[0] + cosphi * vin[2];
      vout[0] = T;
      vout[1] = vin[1];
      break;
    case 3:  /* rotate about z-axis */
      T = cosphi * vin[0] + sinphi * vin[1];
      vout[1] = -sinphi * vin[0] + cosphi * vin[1];
      vout[0] = T;
      vout[2] = vin[2];
      break;
    default:
      LogMsg(stderr,"MRotate: axis not valid.\n");
      exit(1);
  }
}

/*********************************************************************
Name:    GetDpsiDeps
Purpose: Function to compute dpsi, deps, dpsidot, and depsdot.
Inputs:  jed - JED on TDB scale.
Outputs: dpsi    - Nutation in longitude.
         deps    - Nutation in obliquity.
         dpsidot - Derivative of dpsi in radians/day.
         depsdot - Derivative of deps in radians/day.
Returns: Nothing.
Status:  Finished.
Errors:  None known.
*********************************************************************/
void GetDpsiDeps(double jed, double *dpsi, double *deps,
  double *dpsidot, double *depsdot) {

  double L, Ldot, LP, LPdot, F, Fdot, D, Ddot, N, Ndot, LL, LLdot;
  double T, c1, c2, c3, c4, c5, lamp, oamp, ls, os, arg, argdot;
  double funarg[12];
  int j;
  static double dpsiz, depsz, dpsidotz, depsdotz, jedz;

  /* 1980 IAU nutation coefficients */
  static double nc[106][9] = {
    {-171996., -174.2, 92025., 8.9,  0., 0., 0., 0., 1.},
    {   2062.,    0.2,  -895., 0.5,  0., 0., 0., 0., 2.},
    {     46.,     0.,   -24.,  0., -2., 0., 2., 0., 1.},
    {     11.,     0.,     0.,  0.,  2., 0.,-2., 0., 0.},
    {     -3.,     0.,     1.,  0., -2., 0., 2., 0., 2.},
    {     -3.,     0.,     0.,  0.,  1.,-1., 0.,-1., 0.},
    {     -2.,     0.,     1.,  0.,  0.,-2., 2.,-2., 1.},
    {      1.,     0.,     0.,  0.,  2., 0.,-2., 0., 1.},
    { -13187.,   -1.6,  5736., -3.1, 0., 0., 2.,-2., 2.},
    {   1426.,   -3.4,    54., -0.1, 0., 1., 0., 0., 0.},
    {   -517.,    1.2,   224., -0.6, 0., 1., 2.,-2., 2.},
    {    217.,   -0.5,   -95.,  0.3, 0.,-1., 2.,-2., 2.},
    {    129.,    0.1,   -70.,  0.,  0., 0., 2.,-2., 1.},
    {     48.,     0.,     1.,  0.,  2., 0., 0.,-2., 0.},
    {    -22.,     0.,     0.,  0.,  0., 0., 2.,-2., 0.},
    {     17.,   -0.1,     0.,  0.,  0., 2., 0., 0., 0.},
    {    -15.,     0.,     9.,  0.,  0., 1., 0., 0., 1.},
    {    -16.,    0.1,     7.,  0.,  0., 2., 2.,-2., 2.},
    {    -12.,     0.,     6.,  0.,  0.,-1., 0., 0., 1.},
    {     -6.,     0.,     3.,  0., -2., 0., 0., 2., 1.},
    {     -5.,     0.,     3.,  0.,  0.,-1., 2.,-2., 1.},
    {      4.,     0.,    -2.,  0.,  2., 0., 0.,-2., 1.},
    {      4.,     0.,    -2.,  0.,  0., 1., 2.,-2., 1.},
    {     -4.,     0.,     0.,  0.,  1., 0., 0.,-1., 0.},
    {      1.,     0.,     0.,  0.,  2., 1., 0.,-2., 0.},
    {      1.,     0.,     0.,  0.,  0., 0.,-2., 2., 1.},
    {     -1.,     0.,     0.,  0.,  0., 1.,-2., 2., 0.},
    {      1.,     0.,     0.,  0.,  0., 1., 0., 0., 2.},
    {      1.,     0.,     0.,  0., -1., 0., 0., 1., 1.},
    {     -1.,     0.,     0.,  0.,  0., 1., 2.,-2., 0.},
    {  -2274.,   -0.2,   977., -0.5, 0., 0., 2., 0., 2.},
    {    712.,    0.1,    -7.,   0., 1., 0., 0., 0., 0.},
    {   -386.,   -0.4,   200.,   0., 0., 0., 2., 0., 1.},
    {   -301.,     0.,   129., -0.1, 1., 0., 2., 0., 2.},
    {   -158.,     0.,    -1.,   0., 1., 0., 0.,-2., 0.},
    {    123.,     0.,   -53.,   0.,-1., 0., 2., 0., 2.},
    {     63.,     0.,    -2.,   0., 0., 0., 0., 2., 0.},
    {     63.,    0.1,   -33.,   0., 1., 0., 0., 0., 1.},
    {    -58.,   -0.1,    32.,   0.,-1., 0., 0., 0., 1.},
    {    -59.,     0.,    26.,   0.,-1., 0., 2., 2., 2.},
    {    -51.,     0.,    27.,   0., 1., 0., 2., 0., 1.},
    {    -38.,     0.,    16.,   0., 0., 0., 2., 2., 2.},
    {     29.,     0.,    -1.,   0., 2., 0., 0., 0., 0.},
    {     29.,     0.,   -12.,   0., 1., 0., 2.,-2., 2.},
    {    -31.,     0.,    13.,   0., 2., 0., 2., 0., 2.},
    {     26.,     0.,    -1.,   0., 0., 0., 2., 0., 0.},
    {     21.,     0.,   -10.,   0.,-1., 0., 2., 0., 1.},
    {     16.,     0.,    -8.,   0.,-1., 0., 0., 2., 1.},
    {    -13.,     0.,     7.,   0., 1., 0., 0.,-2., 1.},
    {    -10.,     0.,     5.,   0.,-1., 0., 2., 2., 1.},
    {     -7.,     0.,     0.,   0., 1., 1., 0.,-2., 0.},
    {      7.,     0.,    -3.,   0., 0., 1., 2., 0., 2.},
    {     -7.,     0.,     3.,   0., 0.,-1., 2., 0., 2.},
    {     -8.,     0.,     3.,   0., 1., 0., 2., 2., 2.},
    {      6.,     0.,     0.,   0., 1., 0., 0., 2., 0.},
    {      6.,     0.,    -3.,   0., 2., 0., 2.,-2., 2.},
    {     -6.,     0.,     3.,   0., 0., 0., 0., 2., 1.},
    {     -7.,     0.,     3.,   0., 0., 0., 2., 2., 1.},
    {      6.,     0.,    -3.,   0., 1., 0., 2.,-2., 1.},
    {     -5.,     0.,     3.,   0., 0., 0., 0.,-2., 1.},
    {      5.,     0.,     0.,   0., 1.,-1., 0., 0., 0.},
    {     -5.,     0.,     3.,   0., 2., 0., 2., 0., 1.},
    {     -4.,     0.,     0.,   0., 0., 1., 0.,-2., 0.},
    {      4.,     0.,     0.,   0., 1., 0.,-2., 0., 0.},
    {     -4.,     0.,     0.,   0., 0., 0., 0., 1., 0.},
    {     -3.,     0.,     0.,   0., 1., 1., 0., 0., 0.},
    {      3.,     0.,     0.,   0., 1., 0., 2., 0., 0.},
    {     -3.,     0.,     1.,   0., 1.,-1., 2., 0., 2.},
    {     -3.,     0.,     1.,   0.,-1.,-1., 2., 2., 2.},
    {     -2.,     0.,     1.,   0.,-2., 0., 0., 0., 1.},
    {     -3.,     0.,     1.,   0., 3., 0., 2., 0., 2.},
    {     -3.,     0.,     1.,   0., 0.,-1., 2., 2., 2.},
    {      2.,     0.,    -1.,   0., 1., 1., 2., 0., 2.},
    {     -2.,     0.,     1.,   0.,-1., 0., 2.,-2., 1.},
    {      2.,     0.,    -1.,   0., 2., 0., 0., 0., 1.},
    {     -2.,     0.,     1.,   0., 1., 0., 0., 0., 2.},
    {      2.,     0.,     0.,   0., 3., 0., 0., 0., 0.},
    {      2.,     0.,    -1.,   0., 0., 0., 2., 1., 2.},
    {      1.,     0.,    -1.,   0.,-1., 0., 0., 0., 2.},
    {     -1.,     0.,     0.,   0., 1., 0., 0.,-4., 0.},
    {      1.,     0.,    -1.,   0.,-2., 0., 2., 2., 2.},
    {     -2.,     0.,     1.,   0.,-1., 0., 2., 4., 2.},
    {     -1.,     0.,     0.,   0., 2., 0., 0.,-4., 0.},
    {      1.,     0.,    -1.,   0., 1., 1., 2.,-2., 2.},
    {     -1.,     0.,     1.,   0., 1., 0., 2., 2., 1.},
    {     -1.,     0.,     1.,   0.,-2., 0., 2., 4., 2.},
    {      1.,     0.,     0.,   0.,-1., 0., 4., 0., 2.},
    {      1.,     0.,     0.,   0., 1.,-1., 0.,-2., 0.},
    {      1.,     0.,    -1.,   0., 2., 0., 2.,-2., 1.},
    {     -1.,     0.,     0.,   0., 2., 0., 2., 2., 2.},
    {     -1.,     0.,     0.,   0., 1., 0., 0., 2., 1.},
    {      1.,     0.,     0.,   0., 0., 0., 4.,-2., 2.},
    {      1.,     0.,     0.,   0., 3., 0., 2.,-2., 2.},
    {     -1.,     0.,     0.,   0., 1., 0., 2.,-2., 0.},
    {      1.,     0.,     0.,   0., 0., 1., 2., 0., 1.},
    {      1.,     0.,     0.,   0.,-1.,-1., 0., 2., 1.},
    {     -1.,     0.,     0.,   0., 0., 0.,-2., 0., 1.},
    {     -1.,     0.,     0.,   0., 0., 0., 2.,-1., 2.},
    {     -1.,     0.,     0.,   0., 0., 1., 0., 2., 0.},
    {     -1.,     0.,     0.,   0., 1., 0.,-2.,-2., 0.},
    {     -1.,     0.,     0.,   0., 0.,-1., 2., 0., 1.},
    {     -1.,     0.,     0.,   0., 1., 1., 0.,-2., 1.},
    {     -1.,     0.,     0.,   0., 1., 0.,-2., 2., 0.},
    {      1.,     0.,     0.,   0., 2., 0., 0., 2., 0.},
    {     -1.,     0.,     0.,   0., 0., 0., 2., 4., 2.},
    {      1.,     0.,     0.,   0., 0., 1., 0., 1., 0.}
  };

  if (jedz != jed) {
    jedz = jed;

    T = (jed - J2000) / JulCty;

    FunArgIAU(jed, funarg);
    L = funarg[0];
    Ldot = funarg[6];
    LP = funarg[1];
    LPdot = funarg[7];
    F  = funarg[2];
    Fdot  = funarg[8];
    D  = funarg[3];
    Ddot  = funarg[9];
    N  = funarg[4];
    Ndot  = funarg[10];
    LL = funarg[5];
    LLdot = funarg[11];

    /* evaluate the series */
    dpsiz = 0.;  /* initialize to zero */
    depsz = 0.;
    dpsidotz = 0.;
    depsdotz = 0.;
    for (j=0; j<106; j++) {
      lamp = nc[j][0];
      oamp = nc[j][2];
      ls   = nc[j][1];
      os   = nc[j][3];
      c1   = nc[j][4];
      c2   = nc[j][5];
      c3   = nc[j][6];
      c4   = nc[j][7];
      c5   = nc[j][8];
      arg  = c1 * L + c2 * LP + c3 * F + c4 * D + c5 * N;
      arg  = amodulo(arg, TWOPI);
      dpsiz += (lamp + ls * T) * sin(arg);
      depsz += (oamp + os * T) * cos(arg);
      argdot = c1 * Ldot + c2 * LPdot + c3 * Fdot + c4 * Ddot + c5 * Ndot;
      argdot = amodulo(argdot, TWOPI);
      dpsidotz += (lamp + ls * T) * argdot * cos(arg) +
        ls * sin(arg) / JulCty;
      depsdotz -= (oamp + os * T) * argdot * sin(arg) +
        os * cos(arg) / JulCty;
    }

    /* normalize and convert units */
    dpsiz *= 0.0001 * A2R;
    depsz *= 0.0001 * A2R;
    dpsidotz *= 0.0001 * A2R;
    depsdotz *= 0.0001 * A2R;
  }

  *dpsi = dpsiz;
  *deps = depsz;
  *dpsidot = dpsidotz;
  *depsdot = depsdotz;
}

/*********************************************************************
Name:    Obliquity
Purpose: Function to compute the mean obliquity of the ecliptic
         and its derivative using Lieske's formula.
Inputs:  jed1 - Initial JED on the TDB scale.
         jed2 - Final JED on the TDB scale.
         m    - 0 for mean obliquity,
                1 for true obliquity.
Outputs: obl    - Obliquity of the ecliptic in radians.
         obldot - Derivative in radians/day.
Returns: Nothing.
Status:  Finished.
Errors:  None known.
*********************************************************************/
void Obliquity(double jed1, double jed2, int m, double *obl,
  double *obldot) {

  double t1, t2, e0, e1, e2, e3, e4, e5, e6, epsbar;
  double dpsi, deps, dpsidot, depsdot;

  t1 = (jed1 - J2000) / JulCty;
  t2 = (jed2 - jed1) / JulCty;

  e0 = 84381.448;
  e1 =   -46.815;
  e2 =    -0.00059;
  e3 =     0.001813;
  epsbar = e0 + e1 * t1 + e2 * t1 * t1 + e3 * t1 * t1 * t1;
  e1 = -46.815;
  e2 =  -0.00117;
  e3 =   0.005439;
  e4 =  -0.00059;
  e5 =   0.005439;
  e6 =   0.001813;
  *obl = epsbar + (e1 + t1 * (e2 + t1 * e3)) * t2
    + (e4 + e5 * t1) * t2 * t2
    + e6 * t2 * t2 * t2;
  *obldot = e1 + t1 * (e2 + t1 * e3)
    + 2.0 * (e4 + e5 * t1) * t2
    + 3.0 * e6 * t2 * t2;

  if (m == 1) {
    /* need true obliquity */
    GetDpsiDeps(jed2, &dpsi, &deps, &dpsidot, &depsdot);
    /* Unit conversion is needed because obl is  */
    /* in arc seconds, nutations are in radians. */
    *obl = *obl + deps * R2A;
    *obldot = *obldot + depsdot * R2A;
  }

  /* convert to radians and radians/day */
  *obl = *obl * A2R;
  *obldot = *obldot * A2R / JulCty;
}

/*********************************************************************
Name:    pleph
Purpose: Function that reads the JPL planetary ephemeris and gives the
         position and velocity of the point 'targ' with respect to 'cent'.
Inputs:  jd   - JED at which interpolation is wanted.
         targ - Number of target point.
         cent - Number of center point.
         The numbering convention for 'targ' and 'cent' is:
         1 = Mercury           8 = Neptune
         2 = Venus             9 = Pluto
         3 = Earth            10 = Moon
         4 = Mars             11 = Sun
         5 = Jupiter          12 = SSB
         6 = Saturn           13 = EMB
         7 = Uranus           14 = Nutations  (if present)
         15 = Librations (if present)
         If nutations are wanted, set targ = 14. For librations,
         set targ = 15. 'cent' will be ignored on either call.
Outputs: rrd[]  - 6 element array containing the state vector of 'targ'
                  relative to 'cent'. The units are AU and AU/DAY. For
                  librations the units are RAD and RAD/DAY. For
                  nutations the first 4 elements of RRD[] are set to
                  nutations and rates, in RAD and RAD/DAY.
         inside - TRUE if 'jd' is within the ephemeris time span.
                  If not, 'inside' is set to FALSE.
Returns: Nothing.
Status:  Finished.
Errors:  None known.
*********************************************************************/
void pleph(double jd, int targ, int cent, double *rrd, int *inside) {

  int i;
  static double fac;
  static int nemb, ipv, ncmp, lme;
  static int pfirst;
  static double jdtot;
  static double pv[6][13]={0.};
  static double embf[2], ve[2], jed[2];
  static int LList[12], LLst[13];
  static int L[2], tc[2];

  /*
      pv[] is a zero-offset 6x13 matrix. The column number (0-12) specifies
      a body number. The row number (0-5) specifies a component of the
      body's state vector x,y,z,xdot,ydot,zdot in that order.
  */

  /* necessary for zero-offset arrays */
  targ = targ - 1;
  cent = cent - 1;
  /*
      From here on, 'targ' and 'cent' are one less than their
      calling values.
  */

  /* 1st time in, be sure ephemeris is initialized */
  if (!pfirst) {
    pfirst = TRUE;
    ipv = 2; /* we want pos and vel */
    ephopn("");
    ve[0] = 1.0/(1.0+emrat);
    ve[1] = emrat*ve[0];

    jed[0] = 0.0;
    jed[1] = 0.0;

    embf[0] = -1.0;
    embf[1] =  1.0;
    for (i=0; i<12; i++) {
      LList[i] = 0;
    }
    L[0] = 0;
    L[1] = 0;
    tc[0] = 0;
    tc[1] = 0;
    for (i=0; i<13; i++) {
      LLst[i] = i;
      if (i ==  2) LLst[i] = 9;
      if (i == 11) LLst[i] = 10;
      if (i == 12) LLst[i] = 2;
    }
    fac = 0.0;
    nemb = 1;
  }

  /* Initialize jed[] for state() and set up component count */
  jed[0] = jd;
  jed[1] = 0.0;

  jdtot = jed[0] + jed[1];

  if ((jdtot >= SS[0]) && (jdtot <= SS[1])) {
    *inside = TRUE;
  } else {
    *inside = FALSE;
    return;
  }

  ncmp = 3*ipv; /* total number of components */

  /* check for nutation call */
  if (targ == 13) {
    if (ipt[1][11] > 0) {
      LList[10] = ipv;
      state(jed, LList, pv, rrd);
      LList[10] = 0;
      return;
    } else {
      LogMsg(stderr, "pleph: no nutations on the ephemeris file.");
      exit(1);
    }
  }

  /* check for librations */
  if (targ == 14) {
    if (lpt[1] > 0) {
      LList[11] = ipv;
      state(jed, LList, pv, rrd);
      LList[11] = 0;
      for (i=0; i<ncmp; i++) {
        rrd[i] = pv[i][10];
      }
      return;
    } else {
      LogMsg(stderr, "pleph: no librations on the ephemeris file.");
      exit(1);
    }
  }

  /* check for targ = cent */
  if (targ == cent) {
    for (i=0; i<ncmp; i++) {
      rrd[i] = 0.0;
    }
    return;
  }

  /* force barycentric output by state() */
  bsav = bary;
  bary = TRUE;

  /* set up proper entries in LList[] array for state() call */
  tc[0] = targ;
  tc[1] = cent;
  lme = 0;

  for (i=0; i<2; i++) {
    L[i] = LLst[tc[i]];
    if (L[i] < 10) LList[L[i]] = ipv;
    if (tc[i] == 2) {
      lme = 2;
      fac = -ve[0];
    }
    else if (tc[i] == 9) {
      lme = 9;
      fac = ve[1];
    }
    else if (tc[i] == 12) {
      nemb = i;
    }
  }

  if ((LList[9] == ipv) && (L[0] != L[1])) LList[2] = ipv - LList[2];

  /* make call to state() */
  state(jed, LList, pv, rrd);

  /* case: Earth-to-Moon */
  if ((targ == 9) && (cent == 2)) {
    for (i=0; i<ncmp; i++) {
      rrd[i] = pv[i][9];
    }

    /* case: Moon-to-Earth */
  }
  else if ((targ == 2) && (cent == 9)) {
    for (i=0; i<ncmp; i++) {
      rrd[i] = -pv[i][9];
    }

    /* case: EMB-to-Moon or -Earth */
  }
  else if ((targ == 12 || cent == 12) && LList[9] == ipv) {
    for (i=0; i<ncmp; i++) {
      rrd[i] = pv[i][9]*fac*embf[nemb];
    }

    /* otherwise, get Earth or Moon vector and then get output vector */
  } else {
    for (i=0; i<ncmp; i++) {
      pv[i][10] = pvsun[i];
      pv[i][12] = pv[i][2];
      if (lme > 0) pv[i][lme] = pv[i][2]+fac*pv[i][9];
      rrd[i] = pv[i][targ] - pv[i][cent];
    }
  }

  /* clear state() body array and restore barycenter flag */
  LList[2] = 0;
  LList[L[0]] = 0;
  LList[L[1]] = 0;
  bary = bsav;
}

/*********************************************************************
Name:    Pol2Rec
Purpose: Function to convert a polar state vector into a cartesian
         state vector.
         NOTE:  THIS ROUTINE EXPECTS THE POLAR VELOCITY VECTOR TO BE
         THE TOTAL VELOCITY CORRECTED FOR THE EFFECT OF LATITUDE.
Inputs:  a[] - Zero-offset polar state vector.
Outputs: b[] - Zero-offset cartesian state vector.
Returns: Nothing.
Status:  Finished.
Errors:  None known.
*********************************************************************/
void Pol2Rec(double a[], double b[]) {

  double lambda, beta, R, v_lambda, v_beta, v_r;
  double lambda_dot, beta_dot, r_dot, CosL, SinL, CosB, SinB;

  lambda   = a[0];
  beta     = a[1];
  R        = a[2];
  v_lambda = a[3];
  v_beta   = a[4];
  v_r      = a[5];

  /* separate the angluar derivatives from the total velocity components */
  CosL = cos(lambda);
  SinL = sin(lambda);
  CosB = cos(beta);
  SinB = sin(beta);

  lambda_dot = v_lambda / (R * CosB);
  beta_dot = v_beta / R;
  r_dot = v_r;

  /* position vector components */
  b[0] = R * CosL * CosB;
  b[1] = R * SinL * CosB;
  b[2] = R * SinB;

  /* velocity vector components */
  b[3] = r_dot * CosL * CosB - R * lambda_dot * SinL * CosB -
    R * beta_dot * CosL * SinB;
  b[4] = r_dot * SinL * CosB + R * lambda_dot * CosL * CosB -
    R * beta_dot * SinL * SinB;
  b[5] = r_dot * SinB + R * beta_dot * CosB;
}

/*********************************************************************
Name:    PrecessElements
Purpose: Subprogram to transform angular orbital elements
         from equinox eqnx1 to equinox eqnx2.
Inputs:  eqnx1   - String containing initial
                   eqninox.  e.g. B1950 or 2415020.5.
         element - Array containing elements
                   referred to eqnx1 as follows:
                   element[0] = inclination.
                   element[1] = long. of asc. node.
                   element[2] = arg. of peri.
         eqnx2   - String containing final
                   eqninox.  e.g. J2000 or 2451545.
Outputs: element[] - Array containing elements referred to eqnx2.
Returns: Nothing.
Status:  Finished.
Errors:  None known.
*********************************************************************/
void PrecessElements(char *eqnx1, double *element, char *eqnx2) {

  double  jd1, jd2, tt, t, tt2, t2, t3, tmp;
  double  SmallPI, LargePI, pa, cosi, sini, sinisinn, sinicosn;
  double  newi, newnode, sindwsini, cosdwsini, dw, newargp;

  ucase(eqnx1);
  ucase(eqnx2);
  if ((eqnx1[0] == 'B') || (eqnx1[0] == 'J')) {
    Epoch2JED(eqnx1, &jd1);
  } else {
    jd1 = atof(eqnx1);
  }

  if ((eqnx2[0] == 'B') || (eqnx2[0] == 'J')) {
    Epoch2JED(eqnx2, &jd2);
  } else {
    jd2 = atof(eqnx2);
  }

  tt = (jd1 - 2451545.0) / 365250.0;
  t = (jd2 - jd1) / 365250.0;
  tt2 = tt * tt;
  t2 = t * t;
  t3 = t2 * t;

  /* Compute precessional quantities */
  SmallPI = (470.029 - 6.603 * tt + 0.598 * tt2) * t +
    (-3.302 + 0.598 * tt) * t2 + 0.06 * t3;
  LargePI = 174.8763838888889 + (32894.789 * tt) / 3600.0 +
    (60.622 * tt2) / 3600.0 + ((-8698.089 - 50.491 * tt) * t) /
    3600.0 + (3.536 * t2) / 3600.0;
  pa = (50290.966 + 222.226 * tt - 0.042 * tt2) * t +
    (111.113 - 0.042 * tt) * t2 - 0.006 * t3;
  SmallPI = SmallPI * A2R;
  LargePI = LargePI * D2R;
  pa = pa * A2R;

  /* Compute new inclination and node */
  cosi = cos(element[0]) * cos(SmallPI) + sin(element[0]) *
    sin(SmallPI) * cos(LargePI - element[1]);

  sinisinn = sin(element[0]) * sin(LargePI - element[1]);
  sinicosn = -sin(SmallPI) * cos(element[0]) + cos(SmallPI) *
    sin(element[0]) * cos(LargePI - element[1]);

  sini = sqrt(sinisinn * sinisinn + sinicosn * sinicosn);
  newi = atan2(sini, cosi);
  if (newi < 0.0)
    newi += TWOPI;

  tmp = atan2(sinisinn, sinicosn);
  if (tmp < 0.0)
    tmp += TWOPI;
  newnode = pa + LargePI - tmp;
  newnode = amodulo(newnode, TWOPI);

  /* Compute new argument of perihelion */
  sindwsini = sin(SmallPI) * sin(LargePI - element[1]);
  cosdwsini = sin(element[0]) * cos(SmallPI) - cos(element[0]) *
    sin(SmallPI) * cos(LargePI - element[1]);
  dw = atan2(sindwsini, cosdwsini);
  if (dw < 0.0)
    dw += TWOPI;

  newargp = element[2] + dw;
  newargp = amodulo(newargp, TWOPI);

  /* Put new elements in element[] array */
  element[0] = newi;
  element[1] = newnode;
  element[2] = newargp;
}

/*********************************************************************
Name:    QRotate
Purpose: Function to perform a matrix rotation of a state vector through
         a given angle about a desired axis.  A right-handed orthogonal
         coordinate system is assumed, and a 6X6 Q-matrix is used.
         NOTE: The vin[] and vout[] can have the same name in
         the calling program.
Inputs:  vin    - Zero-offset input state vector.
         axis   - 1 for x-axis,
                  2 for y-axis,
                  3 for z-axis.
         phi    - Rotation angle in radians.
         phidot - Derivative of phi.
         s      - 0 for inertial to inertial,
                  1 for inertial to rotating.
Outputs: vout[] - Zero-offset transformed state vector.
Returns: Nothing.
Status:  Finished.
Errors:  None known.
*********************************************************************/
void QRotate(double vin[], int axis, double phi, double phidot,
  int s, double vout[]) {

  int i;
  double temp[6];
  DMatrix QMatrix;

  QMatrix = createDMatrix(6);

  GetQMatrix(phi, phidot, axis, s, QMatrix);

  for (i=0; i<6; i++) {
    temp[i] = vin[i];
  }

  MatXVec(QMatrix, temp, vout, 6, 6);

  freeDMatrix(QMatrix, 6);
}

/*********************************************************************
Name:    RayBend
Purpose: Function to correct an input vector for relativistic light
         deflection due to the Sun's gravity field.
Inputs:  earth_hel[] - Zero-offset heliocentric state vector of Earth.
         body_geo[]  - Zero-offset geocentric state vector of object.
         body_hel[]  - Zero-offset heliocentric state vector of object.
Outputs: p1[] - Zero-offset geocentric state vector of object
                corrected for light deflection.
Returns: Nothing.
Status:  Finished.
Errors:  None known.
*********************************************************************/
void RayBend(double earth_hel[], double body_geo[], double body_hel[],
  double p1[6]) {

  double ue[3], up[3], uq[3], E[3], Edot[3], P[3], Pdot[3], Q[3], Qdot[3];
  double magE, magP, pdotq, edotp, qdote;
  int i;

  /* extract the pos. portions of the state vectors */
  SplitStateVector(earth_hel, E, Edot);
  SplitStateVector( body_geo, P, Pdot);
  SplitStateVector( body_hel, Q, Qdot);

  /* form unit vectors */
  Uvector(E, ue);
  Uvector(P, up);
  Uvector(Q, uq);

  /* form dot products and other quantities */
  magE = Vecmag(E);
  magP = Vecmag(P);
  Vdot(3, up, uq, &pdotq);
  Vdot(3, ue, up, &edotp);
  Vdot(3, uq, ue, &qdote);

  for (i=0; i<3; i++) {
    p1[i] = up[i] + (MUC / magE) * (pdotq * ue[i] - edotp * uq[i])
      / (1.0 + qdote);
    /* make p1[] a non-unit vector */
    p1[i] = magP * p1[i];
    p1[i+3] = body_geo[i+3];
  }
}

/*********************************************************************
Name:    Rec2Pol
Purpose: Subprogram to convert a cartesian state vector into a
         polar state vector. NOTE: THE POLAR VELOCITY VECTOR IS
         THE TOTAL VELOCITY, CORRECTED FOR THE EFFECT OF LATITUDE.
Inputs:  a[] - Zero-offset cartesian state vector.
Outputs: b[] - Zero-offset polar state vector.
Returns: Nothing.
Status:  Finished.
Errors:  None known.
*********************************************************************/
void Rec2Pol(double a[], double b[]) {

  double x, y, z, x_dot, y_dot, z_dot,
    rho, r, lambda, beta, lambda_dot,
    beta_dot, r_dot;

  x = a[0];
  y = a[1];
  z = a[2];
  x_dot = a[3];
  y_dot = a[4];
  z_dot = a[5];

  rho = sqrt(x * x + y * y);
    r = sqrt(rho * rho + z * z);

  lambda = atan2(y, x);
  if (lambda < 0.0)
    lambda += TWOPI;

  beta   = atan2(z, rho);
  if (beta < 0.0)
    beta += TWOPI;

  if (z < 0) {
    beta = beta - TWOPI;
  }

  if (rho == 0) {
    lambda_dot = 0.0;
    beta_dot = 0.0;
  } else {
    lambda_dot = (x*y_dot-y*x_dot) / (rho*rho);
    beta_dot = (z_dot*rho*rho-z*(x*x_dot+y*y_dot))/(r*r*rho);
  }

  r_dot = (x * x_dot + y * y_dot + z * z_dot) / r;

  /* position vector components */
  b[0] = lambda;
  if (b[0] >= TWOPI)
    b[0] = b[0] - TWOPI;

  b[1] = beta;
  b[2] = r;
  /* total velocity vector components */
  b[3] = r * lambda_dot * cos(beta);
  b[4] = r * beta_dot;
  b[5] = r_dot;
}

/*********************************************************************
Name:    Reduce
Purpose: Function for reducing source planetary or solar ephemerides
         to apparent, topocentric, virtual, local, or astrometric
         place. The processes are rigorous and include all corrections.
         This function is intended for use with barycentric source
         ephemerides such as DE200.
Inputs:  jed        - Desired time of reduction.
         body       - Body number (99 for stellar reduction).
         place      - 1 for apparent place, 2 for topocentric place,
                      3 for virtual place, 4 for local place,
                      5 for astrometric place.
         StarData[] - Array containing stellar data.
Outputs: p3[] - Array containing the requested state vector.
Returns: Nothing.
Status:  Finished.
Errors:  None known.
*********************************************************************/
void Reduce(double jed, int body, int place, double StarData[],
  double p3[]) {

  int i;
  double Ltime;
  double earth_ssb[6], earth_hel[6], body_geo[6];
  double body_hel[6], obsr_geo[6];
  double eb[6], ebdot[6], p1[6], p2[6];
  double zeta,z,theta,zetadot,zdot,thetadot;
  double dpsi,deps,dpsidot,depsdot;
  double TrueEps,MeanEps,TrueEpsDot,MeanEpsDot;

  if ((body == 3) || (body == 12) || (body == 13) ||
    (body == 14) || (body == 15)) {
    LogMsg(stderr,
      "Reduce: can't perform reduction for specified body.\n");
    exit(1);
  }

  /* Earth's barycentric state vector */
  GetStateVector(jed, 3, 12, 1, StateVector);
  for (i=0; i<6; i++) {
    earth_ssb[i] = StateVector[3-1][12-1][1-1][i];
  }

  /* Earth's heliocentric state vector */
  GetStateVector(jed, 3, 11, 1, StateVector);
  for (i=0; i<6; i++) {
    earth_hel[i] = StateVector[3-1][11-1][1-1][i];
  }

  if ((place == 2) || (place == 4)) {
    /* Compute geocentric state vector of observer */
    GeocenObs(jed, obsr_geo);
    /* Translate origin from geocenter to topocenter */
    for (i=0; i<6; i++) {
      earth_ssb[i] = earth_ssb[i] + obsr_geo[i];
    }
  }

  /*
      Compute geo/topocentric state vector of object corrected
      for light time.
  */
  LightTime(jed, body, earth_ssb, earth_hel, StarData, body_geo,
    body_hel, &Ltime);

  if (place == 5) {
    /* Compute the astrometric place */
    for (i=0; i<6; i++) {
      p3[i] = body_geo[i];
    }
  }

  if (place < 5) {
    /* Perform correction for relativistic light deflection */
    RayBend(earth_hel, body_geo, body_hel, p1);

    /* Perform correction for aberration */
    SplitStateVector(earth_ssb, eb, ebdot);
    Aberrate(p1, ebdot, p2);

    if (place < 3) {
      /* Correction for precession and nutation from J2000.0 */
      GetPrecessParams(J2000, jed, &zeta, &z, &theta, &zetadot,
        &zdot, &thetadot);
      GetDpsiDeps(jed, &dpsi, &deps, &dpsidot, &depsdot);
      Obliquity(J2000, jed, 0, &MeanEps, &MeanEpsDot);
      Obliquity(J2000, jed, 1, &TrueEps, &TrueEpsDot);

      /* First correct for precession */
      QRotate(p2, 3, -zeta, -zetadot, 1, p3);
      QRotate(p3, 2, theta, thetadot, 1, p3);
      QRotate(p3, 3, -z, -zdot, 1, p3);

      /* Now correct for nutation */
      QRotate(p3, 1, MeanEps, MeanEpsDot, 1, p3);
      QRotate(p3, 3, -dpsi, -dpsidot, 1, p3);
      QRotate(p3, 1, -TrueEps, -TrueEpsDot, 1, p3);
    } else {
      for (i=0; i<6; i++) {
        p3[i] = p2[i];
      }
    }
  }

  return;
}

/*********************************************************************
Name:    Refract
Purpose: Subprogram to correct right ascension and declination
         for atmospheric refraction.
         Reference:  Taff. Computational Spherical Astronomy, pp. 78-80.
Inputs:  ra1   - Uncorrected RA and DEC in radians.
         dec1  - Uncorrected RA and DEC in radians.
         lha   - Local apparent hour angle in radians.
         temp  - Temperature in deg F.
         press - Air pressure in inches of Hg.
Outputs: ra2  - Corrected RA and DEC in radians.
         dec2 - Corrected RA and DEC in radians.
Returns: Nothing.
Status:  Finished.
Errors:  None known.
*********************************************************************/
void Refract(double ra1, double dec1, double lha, double temp,
  double press, double *ra2, double *dec2) {

  double cosz, z, r1, r2, r, tanzr, rr, diff, rlocal, denom;

  cosz = sin(obsr_lat) * sin(dec1) + cos(obsr_lat) * cos(dec1) * cos(lha);
  z = acos(cosz);

  r1 = 58.294 * A2R;
  r2 = -0.0668 * A2R;

  r = 0.0;
  do {
    tanzr = tan(z - r);
    rr = r1 * tanzr + r2 * tanzr * tanzr * tanzr;
    printf("%f\n", rr * R2A);
    diff = fabs(rr - r);
    r = rr;
  } while (diff > 1e-15);

  rlocal = (17.0 * press) * r / (460.0 + temp);

  denom = 1.0 / (cos(dec1) * sin(z));
  *dec2 = dec1 + (rlocal * (sin(obsr_lat) - sin(dec1) * cos(z))) * denom;
  denom = 1.0 / (cos(*dec2) * sin(z));
  *ra2 = ra1 + (rlocal * cos(obsr_lat) * sin(lha)) * denom;
}

/*********************************************************************
Name:    RotMat
Purpose: Function to compute the 3X3 rotation matrix and its partial
         derivative for a rotation about the Ith coordinate axis
         through an angle phi.
Inputs:  axis - 1 for x-axis,
                2 for y-axis,
                3 for z-axis.
         phi  - Rotation angle in radians.
Outputs: r[]  - Zero-offset 3X3 rotation matrix.
         dr[] - Zero-offset 3X3 partial derivative of the matrix r[].
Returns: Nothing.
Status:  Finished.
Errors:  None known.
*********************************************************************/
void RotMat(int axis, double phi, DMatrix r, DMatrix dr) {

  double cosphi, sinphi;

  cosphi = cos(phi);
  sinphi = sin(phi);

  switch (axis) {
    case 1:  /* rotate about x-axis */
      r[0][0] = 1.0;
      r[0][1] = 0.0;
      r[0][2] = 0.0;
      r[1][0] = 0.0;
      r[1][1] =  cosphi;
      r[1][2] =  sinphi;
      r[2][0] = 0.0;
      r[2][1] = -sinphi;
      r[2][2] =  cosphi;
      dr[0][0] = 0.0;
      dr[0][1] = 0.0;
      dr[0][2] = 0.0;
      dr[1][0] = 0.0;
      dr[1][1] = -sinphi;
      dr[1][2] =  cosphi;
      dr[2][0] = 0.0;
      dr[2][1] = -cosphi;
      dr[2][2] = -sinphi;
      break;
    case 2:  /* rotate about y-axis */
      r[0][0] =  cosphi;
      r[0][1] = 0.0;
      r[0][2] = -sinphi;
      r[1][0] = 0.0;
      r[1][1] = 1.0;
      r[1][2] = 0.0;
      r[2][0] =  sinphi;
      r[2][1] = 0.0;
      r[2][2] =  cosphi;
      dr[0][0] = -sinphi;
      dr[0][1] = 0.0;
      dr[0][2] = -cosphi;
      dr[1][0] = 0.0;
      dr[1][1] = 0.0;
      dr[1][2] = 0.0;
      dr[2][0] =  cosphi;
      dr[2][1] = 0.0;
      dr[2][2] = -sinphi;
      break;
    case 3:  /* rotate about z-axis */
      r[0][0] =  cosphi;
      r[0][1] =  sinphi;
      r[0][2] = 0.0;
      r[1][0] = -sinphi;
      r[1][1] =  cosphi;
      r[1][2] = 0.0;
      r[2][0] = 0.0;
      r[2][1] = 0.0;
      r[2][2] = 1.0;
      dr[0][0] = -sinphi;
      dr[0][1] =  cosphi;
      dr[0][2] = 0.0;
      dr[1][0] = -cosphi;
      dr[1][1] = -sinphi;
      dr[1][2] = 0.0;
      dr[2][0] = 0.0;
      dr[2][1] = 0.0;
      dr[2][2] = 0.0;
      break;
    default:
      LogMsg(stderr,"RotMat: axis not valid.\n");
      exit(1);
  }
}

/*********************************************************************
Name:    RST_Interpolate
Purpose: Interpolation routine used by RST.
Inputs:  c      - flag
         z0     - The "standard" zenith distance for the object
                  at rise or set.  This quantity has different
                  values for different objects according to the
                  following table:
                  z0 = 90d 50'          Sun.
                  z0 = 90d 34' + s - pi Moon.
                  z0 = 90d 34'          stars and planets.
                  z0 = 108d             astronomical twilight.
                  z0 = 102d         nautical twilight.
                  z0 =  96d             civil twilight.
                  z0 should be given by the program that calls RST.
                  z0 is expressed in radians.
         oldm   - intermediate value
         gast0  - GAST at 0h
         deltat - TDT-UT1 in seconds of time.
         ra[]   - Array containing the object's apparent right
                  ascension at times jed-1, jed, and jed+1.
         dec[]  - Array containing the object's apparent
                  declination at times jed-1, jed, and jed+1.
         rx, dx - RA and DEC of object on day x
Outputs: newm   - intermediate value
Returns: Nothing.
Status:  Finished.
Errors:  None known.
*********************************************************************/
static void RST_Interpolate(int c, double z0, double oldm, double gast0,
  double deltat, double *ra, double *dec, double r1, double r2, double r3,
  double d1, double d2, double d3, double *newm) {

  double alpha, dm, h, gast, delta, alt, n;

  *newm = oldm;
  do {
    gast = gast0 + 6.300388093 * (*newm);
    gast = amodulo(gast, TWOPI);
    n = *newm + deltat / 86400.0;
    alpha = ra[1] + 0.5 * n * (r1 + r2 + n * r3);
    alpha = amodulo(alpha, TWOPI);
    delta = dec[1] + 0.5 * n * (d1 + d2 + n * d3);
    h = gast + obsr_lon - alpha;
    alt = asin(sin(delta) * sin(obsr_lat) + cos(delta) * cos(obsr_lat) * cos(h));
    if (c == 0) {
      /* h must satisfy -PI <= h <= PI */
      h = amodulo(h, TWOPI);
      if (h > PI) {
        h = h - TWOPI;
      }
      dm = -h / TWOPI;
    } else {
      dm = (alt - PIDIV2 + z0) / (TWOPI * cos(delta) * cos(obsr_lat) * sin(h));
    }
    *newm = (*newm) + dm;
  } while (fabs(dm) >= 1e-15);
}

/*********************************************************************
Name:    RST
Purpose: Subprogram to compute the times of rise, set, and transit for
         any object given the observer's location and arrays containing
         the APPARENT right ascension and declination of the object for
         three dates centered on the input JED.  The algorithm is
         completely rigorous,  and takes into account atmospheric
         refraction and the object's sidereal motion in the intervals
         between rising, setting, and transiting.  The algorithm is
         explained in Chapter 42 of Meeus' ASTRONOMICAL FORMULAE FOR
         CALCULATORS. With the appropriate values for z0, this routine
         can also be used to compute the times of civil, nautical, or
         astronomical twilight. Note that the times are on the UT1 scale!
         Reference:  Meeus.  ASTRONOMICAL FORMULAE FOR CALCULATORS,
         4TH ED., Chapter 42.
Inputs:  jed    - Julian day number at 0h UT1
         ra[]   - Array containing the object's apparent right
                  ascension at times jed-1, jed, and jed+1.
         dec[]  - Array containing the object's apparent
                  declination at times jed-1, jed, and jed+1.
         z0     - The "standard" zenith distance for the object
                  at rise or set.  This quantity has different
                  values for different objects according to the
                  following table:
                  z0 = 90d 50'          Sun.
                  z0 = 90d 34' + s - pi Moon.
                  z0 = 90d 34'          stars and planets.
                  z0 = 108d             astronomical twilight.
                  z0 = 102d         nautical twilight.
                  z0 =  96d             civil twilight.
                  z0 should be given by the program that calls RST.
                  z0 is expressed in radians.
         deltat - TDT-UT1 in seconds of time.
Outputs: ris[] - String containing rise time in
                 form hh.mm or an appropriate symbol.
         trn[] - String containing transit time
                 in form hh.mm or an appropriate symbol.
         set[] - String containing set time in
                 form hh.mm or an appropriate symbol.
Returns: Nothing.
Status:  Finished.
Errors:  None known.
*********************************************************************/
void RST(double jed, double *ra, double *dec, double z0, double deltat,
  char *ris, char *trn, char *set) {

  int rsflag, c;
  double h0, cosh0, newm, oldm, m, m0, m1, m2;
  double ristime, settime, trntime, gast0;
  double d1, d2, d3, r1, r2, r3;

  /* Make sure the ra[]'s are in continuous order */
  if ((ra[1] < ra[0]) && (ra[2] > ra[1])) {
    ra[1] = ra[1] + TWOPI;
    ra[2] = ra[2] + TWOPI;
  }
  else if ((ra[1] > ra[0]) && (ra[2] < ra[1])) {
    ra[2] = ra[2] + TWOPI;
  }

  r1 = ra[1] - ra[0];
  r2 = ra[2] - ra[1];
  r3 = r2 - r1;
  d1 = dec[1] - dec[0];
  d2 = dec[2] - dec[1];
  d3 = d2 - d1;

  rsflag = -1;

  cosh0 = (cos(z0) - sin(obsr_lat) * sin(dec[1])) / (cos(obsr_lat) * cos(dec[1]));

  if (cosh0 < -1.0) {
    /* Object is circumpolar */
    strcpy(ris, "***********");
    if ((z0 * R2D) >= 96.0) {
      strcpy(set, "**BRIGHT***");
    } else {
      strcpy(set, "**NO SET***");
    }
    rsflag = 0;
  }
  else if (cosh0 > 1.0) {
    /* Object never rises */
    if ((z0 * R2D) >= 96.0) {
      strcpy(ris, "---DARK----");
    } else {
      strcpy(ris, "--NO RISE--");
    }
    strcpy(set, "-----------");
    rsflag = 0;
  }

  GetGST(jed, 1, &gast0);

  m0 = (ra[1] - obsr_lon - gast0) / TWOPI;
  m0 = amodulo(m0, 1.0);

  if (rsflag) {
    h0 = acos(cosh0);
    h0 = amodulo(h0, PI);
    m1 = m0 - h0 / TWOPI;
    m1 = amodulo(m1, 1.0);
    m2 = m0 + h0 / TWOPI;
    m2 = amodulo(m2, 1.0);

    /* Rising */
    oldm = m1;
    c = 1;
    RST_Interpolate(c, z0, oldm, gast0, deltat, ra, dec, r1,
      r2, r3, d1, d2, d3, &newm);
    m = newm;
    ristime = 24.0 * m;
    if (ristime > 24.0) {
      ristime = ristime - 24.0;
      /* Event occurs the following day */
      FmtDms(ristime, 0, 1, ris);
      strcat(ris, "(f)");
    }
    else if (ristime < 0.0) {
      ristime = ristime + 24.0;
      /* Event occurs the previous day */
      FmtDms(ristime, 0, 1, ris);
      strcat(ris, "(p)");
    } else {
      FmtDms(ristime, 0, 1, ris);
    }

    /* Setting */
    oldm = m2;
    c = 1;
    RST_Interpolate(c, z0, oldm, gast0, deltat, ra, dec, r1,
      r2, r3, d1, d2, d3, &newm);
    m = newm;
    settime = 24.0 * m;
    if (settime > 24.0) {
      settime = settime - 24.0;
      FmtDms(settime, 0, 1, set);
      strcat(set, "(f)");
    }
    else if (settime < 0.0) {
      settime = settime + 24.0;
      FmtDms(settime, 0, 1, set);
      strcat(set, "(p)");
    } else {
      FmtDms(settime, 0, 1, set);
    }
  }

  /* Transiting */
  oldm = m0;
  c = 0;
  RST_Interpolate(c, z0, oldm, gast0, deltat, ra, dec, r1,
    r2, r3, d1, d2, d3, &newm);
  m = newm;
  trntime = 24.0 * m;
  if (trntime > 24.0) {
    trntime = trntime - 24.0;
    FmtDms(trntime, 0, 1, trn);
    strcat(trn, "(f)");
  }
  else if (trntime < 0.0) {
    trntime = trntime + 24.0;
    FmtDms(trntime, 0, 1, trn);
    strcat(trn, "(p)");
  } else {
    FmtDms(trntime, 0, 1, trn);
  }
}

/*********************************************************************
Name:    split
Purpose: Function to break a number into an integer part and a
         fractional part. For negative input numbers, 'ipart'
         contains the next more negative number and 'fpart' contains
         a positive fraction.
Inputs:  tt - Number to be split.
Outputs: ipart - Integer part.
         fpart - Fractional part.
Returns: Nothing.
Status:  Finished.
Errors:  None known.
*********************************************************************/
void split(double tt, double *ipart, double *fpart) {

  *ipart = floor(tt);
  *fpart = tt - *ipart;

  if ((tt < 0) && (*fpart != 0)) {
    *ipart = *ipart - 1.0;
    *fpart = *fpart + 1.0;
  }
}

/*********************************************************************
Name:    SplitStateVector
Purpose: Subprogram to split a state vector into its position
         and velocity components.
Inputs:  pv[] - Zero-offset 6-d state vector.
Outputs: p[] - Zero-offset 3-d position vector.
         v[] - Zero-offset 3-d velocity vector.
Returns: Nothing.
Status:  Finished.
Errors:  None known.
*********************************************************************/
void SplitStateVector(double pv[], double p[], double v[]) {

  int i;

  for (i=0; i<3; i++) {
    p[i] = pv[i];
    v[i] = pv[i+3];
  }
}

/*********************************************************************
Name:    state
Purpose: This subroutine reads and interpolates the JPL planetary
         ephemeris file.
Inputs:  jed[]   - 2 element array containing the JED epoch at which
                   interpolation is wanted. Any combination of
                   jed[0]+jed[1] which falls within the time span on
                   the file is a permissible epoch. For ease in
                   programming, the user may put the entire epoch in
                   jed[0] and set jed[1] = 0. For maximum accuracy,
                   set jed[0] = most recent midnight at or before
                   interpolation epoch and set jed[1] = fractional
                   part of a day elapsed between jed[0] and epoch.
                   As an alternative, it may prove convenient to set
                   jed[0] = some fixed epoch, such as start of
                   integration and jed[1] = elapsed interval between
                   interval between then and epoch.
         LList[] - 12 element array specifying what interpolation
                   is wanted for each of the bodies on the file.
                   LList[i] =0, no interpolation for body i,
                            =1, position only,
                            =2, position and velocity.
                   The designation of the astronomical bodies by i is:
                   i = 0: Mercury,
                     = 1: Venus,
                     = 2: Emb,
                     = 3: Mars,
                     = 4: Jupiter,
                     = 5: Saturn,
                     = 6: Uranus,
                     = 7: Neptune,
                     = 8: Pluto,
                     = 9: Moon (geocentric),
                     =10: Nutations in longitude and obliquity (if present),
                     =11: Lunar librations (if present).
Outputs: pv[]  - 6 x 13 array that will contain requested interpolated
                 quantities. the body specified by LList[i] will have its
                 state in the array starting at pv[0][i]. (On any given
                 call, only those words in pv[] which are affected by the
                 first 10 LList[] entries (and by LList[12] if librations
                 are on the file) are set. The rest of the pv[] array
                 is untouched.) the order of components starting in
                 pv[0][i] is: x,y,z,dx,dy,dz.

                 All output vectors are referenced to the earth mean
                 equator and equinox of epoch. The Moon state is always
                 geocentric; the other nine states are either heliocentric
                 or solar-system barycentric, depending on the setting of
                 common flags (see below).

                 Lunar librations, if in the data file, are put into
                 pv[k][11] if LList[11] is 1 or 2.

         nut[] - 4-word array that will contain nutations and rates,
                 depending on the setting of LList[10]. The order of
                 quantities in nut[] is:
                 dpsi     (nutation in longitude),
                 depsilon (nutation in obliquity),
                 dpsi dot,
                 depsilon dot.
Returns: Nothing.
Status:  Finished.
Errors:  None known.
*********************************************************************/
void state(double *jed, int LList[], double pv[][13], double *nut) {
/***********************************************************************
* Other important variables:                                           *
*                                                                      *
* km   Logical flag defining physical units of the output              *
*      states. km = TRUE,  units are km and km/sec,                    *
*                 = FALSE, units are au and au/day.                    *
*      default value = FALSE (km determines time unit                  *
*      for nutations and librations. Angle unit is always radians.)    *
*                                                                      *
* bary Logical flag defining output center.                            *
*      only the 9 planets are affected.                                *
*                    bary = TRUE  = center is SSB,                     *
*                         = FALSE = center is Sun.                     *
*           default value = FALSE.                                     *
*                                                                      *
* pvsun[] 6-word array containing the barycentric position and         *
*         velocity of the Sun.                                         *
***********************************************************************/

  static double t[2], jd[4], temp[6];
  static double dumpv[3][2]={{0.,0.},{0.,0.},{0.,0.}};
  static int sfirst;
  static long int nrl;

  static int buff, ncf, na, km = 0;
  static double aufac, s, ipart, fpart;

  int i, j, k, m;
  static long int nr;

  /* 1st time in, get pointer data, etc., from ephemeris file */
  if (!sfirst) {
    sfirst = TRUE;
    aufac = 1.0;
    nrl = 0L;
    ephopn("");
    if (km) {
      t[1] = SS[2]*86400.0;
    } else {
      t[1] = SS[2];
      aufac = 1.0 / au;
    }
  }

  /* main entry point -- check epoch and read right record */
  s = jed[0] - 0.5;
  split(s, &ipart, &fpart);
  jd[0] = ipart;
  jd[1] = fpart;
  split(jed[1], &ipart, &fpart);
  jd[2] = ipart;
  jd[3] = fpart;
  jd[0] = jd[0] + jd[2] + 0.5;
  jd[1] = jd[1] + jd[3];
  split(jd[1], &ipart, &fpart);
  jd[2] = ipart;
  jd[3] = fpart;
  jd[0] = jd[0] + jd[2];

  /* error return of epoch out of range */
  if ((jd[0] < SS[0]) || (jd[0]+jd[3]) > SS[1]) {
    LogMsg(stderr, "state: epoch out of range.\n");
    exit(1);
  }

  /* 'nr' is the byte index of the first coefficient */
  nr = (long) (floor((jd[0]-SS[0])/SS[2]));
  nr = LengthOfHeader + nr * BlockLength;
  /* use previous block if necessary */
  if (jd[0] == SS[1])
    nr = nr - BlockLength;
  if (nr < 1L) {
    LogMsg(stderr, "state: block not present.\n");
    exit(1);
  }

  /* calculate relative time in interval (0 <= t[0] <= 1) */
  t[0] = ((jd[0]-(((double)nr - (double)(LengthOfHeader))/(double)BlockLength
    * SS[2] + SS[0])) + jd[3]) / SS[2];
  if (t[0] < 0.0)
    t[0] = t[0]+ 1.0;

  /* read correct record if not in core */
  if (nr != nrl) {
    nrl = nr;
    if (fseek(fpBinaryFile, (long) (nr), 0) != 0) {
      LogMsg(stderr, "state: fseek() failed.\n");
      exit(1);
    }
    k = 0;
    do {
      if (k < ncoeff) {
        fread(&tmpDouble, sizeof(double), 1, fpBinaryFile);
        convert_little_endian((char *) &tmpDouble, sizeof(double));
        db[k] = (double) tmpDouble;
      }
      k++;
    } while (!feof(fpBinaryFile) && k <= ncoeff);
  }

  /* interpolate SSBARY Sun */
  buff = ipt[0][10]-1; /* location of first coeff */
  ncf  = ipt[1][10]; /* number of coeffs per component */
  na   = ipt[2][10]; /* number of sets of coeffs per 32day interval */

  interp(buff, t, ncf, 3, na, 2, dumpv);

  k = 0;
  for (j=0; j<2; j++) {
    for (i=0; i<3; i++) {
      pvsun[k] = dumpv[i][j]*aufac;
      k++;
    }
  }

  /* check and interpolate whichever bodies are requested */
  for (i=0; i<10; i++) {
    if (LList[i] <= 0)
      continue;
    if (ipt[1][i] <= 0) {
      errprt(i ,"th body requested - not on file.\n");
    }
    buff = ipt[0][i]-1; /* location of first coeff */
    ncf  = ipt[1][i]; /* number of coeffs per component */
    na   = ipt[2][i]; /* number of sets of coeffs per 32day interval */

    interp(buff, t, ncf, 3, na, LList[i], dumpv);

    /* need to re-map dumpv[1..3][1..2] --> temp[1..6] */
    k = 0;
    for (j=0; j<2; j++) {
      for (m=0; m<3; m++) {
        temp[k] = dumpv[m][j];
        k++;
      }
    }

    for (j=0; j<(LList[i]*3); j++) {
      if ((i <= 8) && (!bary)) {
        pv[j][i] = temp[j]*aufac-pvsun[j];
      } else {
        pv[j][i] = temp[j]*aufac;
      }
    }
  }

  /* do nutations if requested and if on file */
  if ((LList[10] > 0) && (ipt[1][11] > 0)) {
    buff = ipt[0][11]-1; /* location of first coeff */
    ncf  = ipt[1][11]; /* number of coeffs per component */
    na   = ipt[2][11]; /* number of sets of coeffs per 32day interval */

    interp(buff, t, ncf, 2, na, LList[10], dumpv);

    /* need to re-map dumpv(1:3,1:2) --> temp(1:6) */
    k = 0;
    for (j=0; j<2; j++) {
      for (m=0; m<2; m++) {
        nut[k] = dumpv[m][j];
        k++;
      }
    }
    nut[4] = 0.0;
    nut[5] = 0.0;
  }

  /* get librations if requested and if on file */
  if ((lpt[1] > 0) && (LList[11] > 0)) {
    buff = lpt[0]-1; /* location of first coeff */
    ncf  = lpt[1]; /* number of coeffs per component */
    na   = lpt[2]; /* number of sets of coeffs per 32day interval */

    interp(buff, t, ncf, 3, na, LList[11], dumpv);

    pv[0][10] = dumpv[0][0];
    pv[1][10] = dumpv[1][0];
    pv[2][10] = dumpv[2][0];
    pv[3][10] = dumpv[0][1];
    pv[4][10] = dumpv[1][1];
    pv[5][10] = dumpv[2][1];
  }
}

/*********************************************************************
Name:    Stat2Elem
Purpose: Subprogram to transform the components of a state vector into
         the universal orbital element set at a given time.
         Reference: Mansfield. 1986. AIAA Paper 86-2269-CP.
Inputs:  posvel[] - State vector.
         mu       - Gravitational const.
         jed      - Time.
Outputs: uelement[] - Array of universal elements:
                      uelement[0] = q.
                      uelement[1] = e.
                      uelement[2] = i.
                      uelement[3] = node.
                      uelement[4] = arg.peri.
                      uelement[5] = T.
Returns: Nothing.
Status:  Finished.
Errors:  None known.
*********************************************************************/
void Stat2Elem(double *posvel, double mu, double jed, double *uelement) {

  int i;
  double r[3], rdot[3], L[3], wvec[3], u[3], v[3], pvec[3], qvec[3];
  double magr, magL, p, vsquared, alpha, chi, psi, sini, cosi, sinn, cosn;
  double rnudot, rrdot, magrdot, esinnu, ecosnu, sinnu, cosnu;
  double z, w, h, s, x, c3, tspp;

  /* Get pos. and vel. vectors from state vector */
  r[0] = posvel[0];
  r[1] = posvel[1];
  r[2] = posvel[2];
  rdot[0] = posvel[3];
  rdot[1] = posvel[4];
  rdot[2] = posvel[5];

  /* Compute magr */
  magr = Vecmag(r);

  /* Compute angular momentum vector */
  Vcross(r, rdot, L);

  /* Compute magL */
  magL = Vecmag(L);

  /* Compute wvec[] */
  Uvector(L, wvec);

  /* Compute u[] */
  Uvector(r, u);

  /* Compute v[] */
  Vcross(wvec, u, v);

  /* Compute semilatus rectum */
  p = magL * magL / mu;

  /* Compute square of velocity */
  Vdot(3, rdot, rdot, &vsquared);

  /* Compute alpha */
  alpha = 2.0 * mu / magr - vsquared;

  /* Compute eccentricity */
  uelement[1] = sqrt(1.0 - alpha * magL * magL / (mu * mu));

  /* Compute perihelion distance */
  uelement[0] = p / (1.0 + uelement[1]);

  /* Compute node and inclination */
  /* First compute chi and psi */
  chi = wvec[0] / (1.0 + wvec[2]);
  psi = -wvec[1] / (1.0 + wvec[2]);
  /* Now get inclination */
  sini = 2.0 * sqrt(chi * chi + psi * psi) / (1.0 + chi * chi + psi * psi);
  cosi = (1.0 - chi * chi - psi * psi) / (1.0 + chi * chi + psi * psi);
  uelement[2] = atan2(sini, cosi);
  if (uelement[2] < 0.0)
    uelement[2] += TWOPI;

  /* Now get node */
  sinn = chi / sqrt(chi * chi + psi * psi);
  cosn = psi / sqrt(chi * chi + psi * psi);
  uelement[3] = atan2(sinn, cosn);
  if (uelement[3] < 0.0)
    uelement[3] += TWOPI;

  /* Compute arg. peri. */
  /* Compute rnudot */
  Vdot(3, rdot, v, &rnudot);
  /* Compute magrdot */
  Vdot(3, r, rdot, &rrdot);
  magrdot = rrdot / magr;
  esinnu = magrdot * sqrt(p / mu);
  ecosnu = rnudot * sqrt(p / mu) - 1.0;
  /* Proceed to compute arg. peri. */
  z = esinnu / (uelement[1] + ecosnu);
  sinnu = 2.0 * z / (1.0 + z * z);
  cosnu = (1.0 - z * z) / (1.0 + z * z);
  /* Get the pvec[] and qvec[] vectors */
  for (i = 0; i < 3; i++) {
    pvec[i] = u[i] * cosnu - v[i] * sinnu;
    qvec[i] = u[i] * sinnu + v[i] * cosnu;
  }

  /* Finally compute arg. peri. */
  uelement[4] = atan2(pvec[2], qvec[2]);
  if (uelement[4] < 0.0)
    uelement[4] += TWOPI;

  /* Compute time of peri. passage */
  w = sqrt(uelement[0] / (mu * (1.0 + uelement[1]))) * z;
  h = alpha * w * w;
  s = 2.0 * ((((((h / 13.0 - 1.0 / 11.0) * h + 1.0 / 9.0) * h - 1.0 / 7.0)
    * h + 1.0 / 5.0) * h - 1.0 / 3.0) * h + 1.0) * w;
  /* Compute Stumpff functions */
  x = alpha * s * s;
  c3 = StumpffN(x, 3);
  tspp = uelement[0] * s + mu * uelement[1] * s * s * s * c3;
  uelement[5] = jed - tspp;
}
/*********************************************************************
Name:    StumpffN
Purpose: Function to compute the nth order Stumpff function for any x.
         Reference: Danby. FUN, pp. 172-174.
Inputs:  x      - Argument.
         Norder - Order desired.
Outputs: None.
Returns: nth order Stumpff function.
Status:  Finished.
Errors:  None known.
*********************************************************************/
double StumpffN(double x, int Norder) {

  int n = 0;
  double a, b, c0, c1, c2, c3;

  do {
    n++;
    x = x / 4.0;
  } while (fabs(x) > 0.1);

  a = (1.0 - x * (1.0 - x / 182.0) / 132.0);
  b = (1.0 - x * (1.0 - x * a / 90.0) / 56.0);
  c2 = (1.0 - x * (1.0 - x * b / 30.0) / 12.0) / 2.0;
  a = (1.0 - x * (1.0 - x / 210.0) / 156.0);
  b = (1.0 - x * (1.0 - x * a / 110.0) / 72.0);
  c3 = (1.0 - x * (1.0  - x * b / 42.0)/ 20.0) / 6.0;

  c1 = 1.0 - x * c3;
  c0 = 1.0 - x * c2;

  do {
    n--;
    c3 = (c2 + c0 * c3) / 4.0;
    c2 = c1 * c1 / 2.0;
    c1 = c0 * c1;
    c0 = 2.0 * c0 * c0 - 1.0;
    x = x * 4.0;
  } while (n > 0);

  switch (Norder) {
    case 0:
      return (c0);
      case 1:
      return (c1);
    case 2:
      return (c2);
      case 3:
      return (c3);
    default:
      LogMsg(stderr, "StumpffN: Norder not 1, 2, or 3\n");
      exit(1);
  }

  return(-999);
}

/*********************************************************************
Name:    Transpose
Purpose: Function to compute the transpose of a matrix.
Inputs:  a[][] - Zero-offset matrix a (n rows by n columns).
Outputs: b[][] - Zero-offset matrix transpose (n rows by n columns).
Returns: Nothing.
Status:  Finished.
Errors:  None known.
*********************************************************************/
void Transpose(DMatrix a, DMatrix b, int n) {

  int i, j;

  for (i=0; i<n; i++) {
    for (j=0; j<n; j++) {
      b[i][j] = a[j][i];
    }
  }
}

/*********************************************************************
Name:    Uvector
Purpose: Unit vector subroutine.
Inputs:  a[] - Zero-offset column vector (3 rows by 1 column).
Outputs: unita[] - Zero-offset unit vector (3 rows by 1 column).
Returns: Nothing.
Status:  Finished.
Errors:  None known.
*********************************************************************/
void Uvector(double a[], double unita[]) {

  double maga;
  int i;

  maga = sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]);
  for (i=0; i<3; i++) {
    if (maga != 0) {
      unita[i] = a[i]/maga;
    } else {
      unita[i] = 0.0;
    }
  }
}

/*********************************************************************
Name:    Vcross
Purpose: Vector cross product subroutine.
Inputs:  a - Zero-offset vector a (3 rows by 1 column).
         b - Zero-offset vector b (3 rows by 1 column).
Outputs: acrossb - a X b (3 rows by 1 column).
Returns: Nothing.
Status:  Finished.
Errors:  None known.
*********************************************************************/
void Vcross(double a[], double b[], double acrossb[]) {

  acrossb[0] = a[1] * b[2] - a[2] * b[1];
  acrossb[1] = a[2] * b[0] - a[0] * b[2];
  acrossb[2] = a[0] * b[1] - a[1] * b[0];
}

/*********************************************************************
Name:    Vdot
Purpose: Vector dot product function.
Inputs:  n   - Number of rows.
         a[] - Zero-offset column vector with N rows.
         b[] - Zero-offset column vector with N rows.
Outputs: adotb - Dot product of a and b.
Returns: Nothing.
Status:  Finished.
Errors:  None known.
*********************************************************************/
void Vdot(int n, double a[], double b[], double *adotb) {

  int i;

  *adotb = 0.0;
  for (i = 0; i<n; i++) {
    *adotb += a[i] * b[i];
  }
}

/*********************************************************************
Name:    Vecmag
Purpose: Vector magnitude function.
Inputs:  a[] - Zero-offset column vector (3 rows by 1 column).
Outputs: None.
Returns: Magnitude of vector.
Status:  Finished.
Errors:  None known.
*********************************************************************/
double Vecmag (double a[]) {

  double x;

  x = sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);

  return (x);
}

/* End Of File - astrolib.c *****************************************/
