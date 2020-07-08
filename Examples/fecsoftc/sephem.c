/*********************************************************************
Project:  sephem
Filename: sephem.c
Author:   Joe Heafner.
Purpose:  A small command line ephemeris program.
Thanks to Charles Gamble for extensive modifications.
*********************************************************************/

/* Header Files *****************************************************/
#include  <stdio.h>
#include  <stdlib.h>
#include  <string.h>
#include  <errno.h>
#include  "astrolib.h"
#include  "support.h"

/* Function Prototypes **********************************************/
void ParseCmdLine(int argc, char *argv[], char *InFile);
void PrintBanner(void);
void PrintUsage(void);

/* Globals **********************************************************/
char szVersion[] = "SEPHEM v1.060300c";
char szLogFile[] = "sephem.log";
extern short int NUMDE;
extern double SS[];
extern double obsr_lon, obsr_lat, obsr_ele;

/*********************************************************************
Name:    main
Purpose: Main routine for sephem.
Inputs:  argc - Number of command-line arguments.
         argv - Pointer to array of command-line arguments.
Outputs: None.
Returns: 0 if execution successful.
Status:  Finished.
Errors:  None known.
*********************************************************************/
int main (int argc, char *argv[]) {

  char InFile[MAX_NAME_SIZE+1] = "";
  char RA[100], DEC[100];
  char line[1024];
  char PL[20];
  int RedGeo = 0, SolSysSum = 0, Place = 0;
  int i, Targ, Cent, NumSteps, counter;
  double TDB = 0.0;
  double xx, yy, zz;
  double xxdot, yydot, zzdot;
  double r, RA_hours, DEC_deg, DateInc;
  double StarData[6], p3[6];

  /* Open the log file */
  if (LogOpen(szLogFile) == FALSE) {
    fprintf(stderr,
      "Could not open log file '%s': %s\n\n", szLogFile, strerror(errno));
    exit(1);    /* Exit with an error code */
  }

  PrintBanner();
  ParseCmdLine(argc, argv, InFile);

  if (strlen(InFile) == 0) {
    /* Prompt user for file name */
    do {
      fprintf(stdout, "File name to use (XXXX to end): ");
      fflush(stdout);
      fgets(InFile, MAX_NAME_SIZE+1, stdin);

      /* Remove whitespace from either end */
      Trim(InFile);
    } while (strlen(InFile) == 0);
  }

  /* We have a filename now */
  /* ucase(InFile); */  /* May not be used for OS/2 and UNIX compiles */
                        /* if you want mixed-case filenames.          */

  /* Test for user exit request */
  if (strcmp(InFile, "XXXX") == 0) {
    LogMsg(stdout, "\nOK\n");
    LogClose();
    remove(szLogFile);
    return (0);
  }

  /* Test to see if filename exists */
  if (!fexist(InFile)) {
    LogMsg(stdout, "Requested input file does not exist.\n");
    LogMsg(stdout, "Specify another file or move to another directory.\n");
    LogMsg(stdout, "\nOK\n");
    LogClose();
    exit(1);
  }

  /* Open the input file and read in the header info */
  if (ephopn(InFile) == NULL) {
    LogMsg(stderr, "An error occurred in ephopn().\n");
    LogClose();
    exit(1);
  }

  LogMsg(stdout,
    "Using ephemeris DE%d with dates from %lf to %lf\n",
    NUMDE, SS[0], SS[1]);

  do {
    fprintf(stdout, "Enter start JED on TDB time scale: ");
    fflush(stdout);
    fgets(line, sizeof(line), stdin);
  } while (sscanf(line, "%lf", &TDB) != 1);

  do {
    do {
      fprintf(stdout, "[1] REDUCED PLACE [2] GEOMETRIC PLACE: ");
      fflush(stdout);
      fgets(line, sizeof(line), stdin);
    } while (sscanf(line, "%d", &RedGeo) != 1);
  } while ((RedGeo != 1) && (RedGeo != 2));

  if (RedGeo == 1) {
    do {
      do {
        fprintf(stdout, "[1] SOLAR SYSTEM SUMMARY [2] SPECIFIC BODY: ");
        fflush(stdout);
        fgets(line, sizeof(line), stdin);
      } while (sscanf(line, "%d", &SolSysSum) != 1);
    } while ((SolSysSum != 1) && (SolSysSum != 2));

    do {
      do {
        fprintf(stdout, "[1] APPARENT [2] TOPOCENTRIC [3] VIRTUAL "
                        "[4] LOCAL [5] ASTROMETRIC: ");
        fflush(stdout);
        fgets(line, sizeof(line), stdin);
      } while (sscanf(line, "%d", &Place) != 1);
    } while ((Place < 1) || (Place > 5));

    if ((Place == 2) || (Place == 4)) {
      do {
        obsr_lat = 0.0;
        fprintf(stdout, "Enter latitude (dd.mmss): ");
        fflush(stdout);
        fgets(line, sizeof(line), stdin);
      } while (sscanf(line, "%lf", &obsr_lat) != 1);
      obsr_lat = deg(obsr_lat) * D2R;

      do {
        obsr_lon = 0.0;
        fprintf(stdout, "Enter longitude (dd.mmss): ");
        fflush(stdout);
        fgets(line, sizeof(line), stdin);
      } while (sscanf(line, "%lf", &obsr_lon) != 1);
      obsr_lon = deg(obsr_lon) * D2R;

      do {
        obsr_ele = 0.0;
        fprintf(stdout, "Enter elevation (meters): ");
        fflush(stdout);
        fgets(line, sizeof(line), stdin);
      } while (sscanf(line, "%lf", &obsr_ele) != 1);
    }

    switch (Place) {
      case 1:
        strcpy(PL, "APPARENT");
        break;
      case 2:
        strcpy(PL, "TOPOCENTRIC");
        break;
      case 3:
        strcpy(PL, "VIRTUAL");
        break;
      case 4:
        strcpy(PL, "LOCAL");
        break;
      case 5:
        strcpy(PL, "ASTROMETRIC");
        break;
    }

    if (SolSysSum == 1) {
      LogMsg(stdout,
        "%s COORDINATES OF SOLAR SYSTEM BODIES ON %lf\n", PL, TDB);
      for (i = 1; i <= 11; i++) {
        if (i != 3) {
          Reduce(TDB, i, Place, StarData, p3);
          r = sqrt(p3[0]*p3[0]+p3[1]*p3[1]+p3[2]*p3[2]);
          RA_hours = atan2(p3[1], p3[0]) * R2H;
          /* adjust for C's atan2() function */
          if (RA_hours < 0.0) {
            RA_hours += 24.0;
          }
          FmtDms(RA_hours, 3, 1, RA);
          DEC_deg = asin(p3[2] / r) * R2D;
          FmtDms(DEC_deg, 2, 0, DEC);

          LogMsg(stdout,
            "BODY %2d ALPHA %s DELTA %s DIST %12.9lf\n", i, RA, DEC, r);
        }
      }
    } else {
      Targ = 12;
      do {
        do {
          fprintf(stdout, "Enter body number (1-11,99): ");
          fflush(stdout);
          fgets(line, sizeof(line), stdin);
        } while (sscanf(line, "%d", &Targ) != 1);
      } while ((Targ < 1 || Targ > 11) && Targ != 99);

      if (Targ == 99) {
        /* Read in star's FK5 catalog data */
        /* GetStarData(FileName, StarData[]); */
        /*
            The previous line MUST be modified by the user to
            indicate the correct name of the data file that
            contains the stellar positional data.
        */

        /* FK5 data for Betelgeuse */
        /*
            Un-comment these lines for default data
            StarData[0] = 5.9195297222 * H2R
            StarData[1] = 7.4070416667 * D2R
            StarData[2] = 0.005  (arcsec)
            StarData[3] = 0.1730 (s/cty)
            StarData[4] = 0.8700 ("/cty)
            StarData[5] = 21     (km/s)
        */

        do {
          fprintf(stdout, "Enter FK5 right ascension (hh.mmss): ");
          fflush(stdout);
          fgets(line, sizeof(line), stdin);
        } while (sscanf(line, "%lf", &StarData[0]) != 1);
        StarData[0] = deg(StarData[0]) * H2R;

        do {
          fprintf(stdout, "Enter FK5 declination (dd.mmss): ");
          fflush(stdout);
          fgets(line, sizeof(line), stdin);
        } while (sscanf(line, "%lf", &StarData[1]) != 1);
        StarData[1] = deg(StarData[1]) * D2R;

        do {
          fprintf(stdout, "Enter parallax (arcsec): ");
          fflush(stdout);
          fgets(line, sizeof(line), stdin);
        } while (sscanf(line, "%lf", &StarData[2]) != 1);

        do {
          fprintf(stdout, "Enter mu alpha (sec/cty): ");
          fflush(stdout);
          fgets(line, sizeof(line), stdin);
        } while (sscanf(line, "%lf", &StarData[3]) != 1);

        do {
          fprintf(stdout, "Enter mu delta (arcsec/cty): ");
          fflush(stdout);
          fgets(line, sizeof(line), stdin);
        } while (sscanf(line, "%lf", &StarData[4]) != 1);

        do {
          fprintf(stdout, "Enter radial velocity (km/s): ");
          fflush(stdout);
          fgets(line, sizeof(line), stdin);
        } while (sscanf(line, "%lf", &StarData[5]) != 1);
      }

      do {
        fprintf(stdout, "Enter increment in days: ");
        fflush(stdout);
        fgets(line, sizeof(line), stdin);
      } while (sscanf(line, "%lf", &DateInc) != 1);

      do {
        do {
          fprintf(stdout, "Enter number of steps: ");
          fflush(stdout);
          fgets(line, sizeof(line), stdin);
        } while (sscanf(line, "%d", &NumSteps) != 1);
      } while (NumSteps <= 0);

      LogMsg(stdout, "%s EPHEMERIS OF BODY %d\n", PL, Targ);
      counter = 0;

      do {
        Reduce(TDB, Targ, Place, StarData, p3);
        r = sqrt(p3[0] * p3[0] + p3[1] * p3[1] + p3[2] * p3[2]);
        RA_hours = atan2(p3[1], p3[0]) * R2H;
        /* adjust for C's atan2() function */
        if (RA_hours < 0.0) {
          RA_hours += 24.0;
        }
        FmtDms(RA_hours, 3, 1, RA);
        DEC_deg = asin(p3[2] / r) * R2D;
        FmtDms(DEC_deg, 2, 0, DEC);

        LogMsg(stdout, "%11.3lf ALPHA %s", TDB, RA);

        if (Targ != 99) {
          LogMsg(stdout, " DELTA %s", DEC);
          LogMsg(stdout, " DIST %11.9lf\n", r);
        } else {
          LogMsg(stdout, " DELTA %s\n", DEC);
        }

        counter++;
        TDB += DateInc;
      } while (counter < NumSteps);
    }
  } else {
    do {
      do {
        fprintf(stdout, "Enter target number (1-15): ");
        fflush(stdout);
        fgets(line, sizeof(line), stdin);
      } while (sscanf(line, "%d", &Targ) != 1);
    } while ((Targ < 1) || (Targ > 15));

    do {
      do {
        fprintf(stdout, "Enter center number (1-15): ");
        fflush(stdout);
        fgets(line, sizeof(line), stdin);
      } while (sscanf(line, "%d", &Cent) != 1);
    } while ((Cent < 1) || (Cent > 15));

    do {
      fprintf(stdout, "Enter increment in days: ");
      fflush(stdout);
      fgets(line, sizeof(line), stdin);
    } while (sscanf(line, "%lf", &DateInc) != 1);

    do {
      do {
        fprintf(stdout, "Enter number of steps: ");
        fflush(stdout);
        fgets(line, sizeof(line), stdin);
      } while (sscanf(line, "%d", &NumSteps) != 1);
    } while (NumSteps <= 0);

    LogMsg(stdout, "GEOMETRIC EPHEMERIS OF BODY %d WRT BODY %d\n",
      Targ, Cent);
    counter = 0;
    do {
      GetStateVector(TDB, Targ, Cent, 1, StateVector);
      xx = StateVector[Targ-1][Cent-1][0][0];
      yy = StateVector[Targ-1][Cent-1][0][1];
      zz = StateVector[Targ-1][Cent-1][0][2];
      xxdot = StateVector[Targ-1][Cent-1][0][3];
      yydot = StateVector[Targ-1][Cent-1][0][4];
      zzdot = StateVector[Targ-1][Cent-1][0][5];
      LogMsg(stdout, "%13.5lf %+14.10lf %+14.10lf %+14.10lf AU\n",
        TDB, xx, yy, zz);
      LogMsg(stdout, "              %+14.10lf %+14.10lf %+14.10lf AU/DAY\n",
        xxdot, yydot, zzdot);
      counter++;
      TDB += DateInc;
    } while (counter < NumSteps);
  }

  LogMsg(stdout, "\nOK\n");
  LogClose(); /* Close the log file */
  return(0);
}

/*********************************************************************
Name:    ParseCmdLine
Purpose: Routine to parse the command line.
Inputs:  argc - Number of command-line arguments.
         argv - Pointer to array of command-line arguments.
Outputs: InFile - Filename given by user for input file.
Returns: Nothing.
Status:  Finished.
Errors:  None known.
*********************************************************************/
void ParseCmdLine(int argc, char *argv[], char *InFile) {

  /*
      Command line parser ported from Basic and optimized for C
      by Varian Swieter.
  */

  int i;

  for (i = 1; i < argc; i++) {
    /*
        Find next occurance of "/" or "-" by looking at
        the first character of each argument.
    */
    if (argv[i][0] != '-' && argv[i][0] != '/') {
      /* Found an argument that did not begin with "/" or "-" */
      LogMsg(stderr,
        "Command line arguments must begin with '-' or '/'.\n");
      PrintUsage();
      LogClose();
      remove(szLogFile);
      exit(1);
    }

    switch (argv[i][1])  {
      case '\0':
        LogMsg(stderr, "Space not allowed between - and option.\n");
        LogMsg(stderr, "Type SEPHEM -h for help.\n");
        LogClose();
        remove(szLogFile);
        exit(1);
      case 'I':
      case 'i':
        if (strlen(argv[i]) == 2) {
          /* We were given nothing after the "-I" so return empty string */
          strcpy(InFile, "");
        } else {
          if (argv[i][2] != ':') {
            /* User missed out colon so return empty string */
            strcpy(InFile, "");
          } else {
            strcpy(InFile, &(argv[i][3]));
          }
        }
        break;
      case '?': /* Using this for help will cause problems under UNIX */
                /* because it is used by the shell for substitution.  */
      case 'h':
        PrintUsage();
        LogMsg(stdout, "OK\n");
        LogClose();
        remove(szLogFile);
        exit(0);
      case 'v': /* Undocumented command line option */
                /* to print version information.    */
        LogMsg(stdout, "OK\n");
        LogClose();
        remove(szLogFile);
        exit(0);
      default:
        LogMsg(stderr, "Option not recognized.\n");
        LogMsg(stderr, "Type SEPHEM -h for help.\n");
        LogClose();
        remove(szLogFile);
        exit(1);
    }
  }
}

/*********************************************************************
Name:    PrintBanner
Purpose: Prints a banner to stdout and log file.
Inputs:  None.
Outputs: None.
Returns: Nothing.
Status:  Finished.
Errors:  None known.
*********************************************************************/
void PrintBanner() {

  LogMsg(stdout, "\n");
  LogMsg(stdout, "***Program %s  ", szVersion);
  LogMsg(stdout, "                         Written by Joe Heafner\n");
  LogMsg(stdout, "***Small ephemeris program\n");
  LogMsg(stdout, "***The author can be reached via Internet:");
  LogMsg(stdout, "             heafnerj@interpath.com\n");
  LogMsg(stdout, "\n");
}

/*********************************************************************
Name:    PrintUsage
Purpose: Prints usage information.
Inputs:  None.
Outputs: None.
Returns: Nothing.
Status:  Finished.
Errors:  None known.
*********************************************************************/
void PrintUsage() {

  printf("Usage: SEPHEM [-i:FILE][-h]\n");
  printf("Valid command line options are:\n");
  printf("\n");
  printf("    -i:FILE   Use FILE as input\n");
  printf("    -h        Display this help screen\n");
  printf("\n");
  printf("Command line options may be in any order.\n");
  printf("\n");
}

/* End Of File - sephem.c *******************************************/
