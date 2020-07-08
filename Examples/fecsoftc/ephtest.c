/*********************************************************************
Project:  ephtest
Filename: ephtest.c
Author:   Joe Heafner
Purpose:  Tests a binary ephemeris file against the correct test data
          obtained from JPL.
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
char szVersion[] = "EPHTEST v1.060300c";
char szLogFile[] = "ephtest.log";
extern short int NUMDE;

/*********************************************************************
Name:    main
Purpose: Main routine for ephtest.
Inputs:  argc - Number of command-line arguments.
         argv - Pointer to array of command-line arguments.
Outputs: None.
Returns: 0 if execution successful.
Status:  Finished.
Errors:  None known.
*********************************************************************/
int main(int argc, char *argv[]) {

  char InFile[MAX_NAME_SIZE+1] = "";
  char TestFile[MAX_NAME_SIZE+1] = "";
  char line[1024] = "";
  FILE *fpTestFile = NULL;
  int  bFirst, bInside, i;
  int  ENUM = 0;
  int  Target = 0;
  int  Center = 0;
  int  Component = 0;
  double JDate = 0.0;
  double CValue = 0.0, MyValue = 0.0, Residual = 0.0;
  double RRD[6];

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
      fgets(InFile, sizeof(InFile), stdin);

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
  if (!fexist(InFile)){
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

  /* See if the test file exists */
  sprintf(TestFile, "testpo.%d", NUMDE);
  if (!fexist(TestFile)) {
    LogMsg(stderr, "Test file '%s' not found.\n", TestFile);
    LogMsg(stderr, "\nOK\n");
    LogClose();
    exit(1);
  }

  /* Open the test file */
  if ((fpTestFile = fopen(TestFile, "r")) == NULL) {
    LogMsg(stderr,
      "Could not open test file '%s': %s.\n", TestFile, strerror(errno));
    LogClose();
    exit(1);
  }

  /* Read in the first six lines, which we don't need */
  for (i = 1; i <= 6; i++) {
    if (!fgets(line, sizeof(line), fpTestFile)) {
      LogMsg(stderr,
        "An error occurred reading line %d from '%s'\n", i, TestFile);
      LogClose();
      fclose(fpTestFile);
      exit(1);
    }
  }

  /* Setup some boolean values */
  bInside = FALSE;
  bFirst = TRUE;

  /* Read in one line at a time, and parse the data */
  do {
    /* Read a line in */
    fgets(line, sizeof(line), fpTestFile);
    /* Force termination - just in case */
    line[sizeof(line)-1] = '\0';
    /* Strip off trailing whitespace */
    RightTrim(line);

    if (strlen(line) > 0) {
      /* Extract the individual values */
      if (sscanf(&line[0], "%d", &ENUM) != 1) {
        LogMsg(stderr, "Invalid data found in '%s'\n", TestFile);
        LogClose();
        fclose(fpTestFile);
        exit(1);
      }

      if (sscanf(&line[16], "%lf", &JDate) != 1) {
        LogMsg(stderr, "Invalid data found in '%s'\n", TestFile);
        LogClose();
        fclose(fpTestFile);
        exit(1);
      }

      if (sscanf(&line[26], "%d", &Target) != 1) {
        LogMsg(stderr, "Invalid data found in '%s'\n", TestFile);
        LogClose();
        fclose(fpTestFile);
        exit(1);
      }

      if (sscanf(&line[29], "%d", &Center) != 1) {
        LogMsg(stderr, "Invalid data found in '%s'\n", TestFile);
        LogClose();
        fclose(fpTestFile);
        exit(1);
      }

      if (Center == 0)
        Center = 11;

      if (sscanf(&line[33], "%d", &Component) != 1) {
        LogMsg(stderr, "Invalid data found in '%s'\n", TestFile);
        LogClose();
        fclose(fpTestFile);
        exit(1);
      }


      if (sscanf(&line[35], "%lf", &CValue) != 1) {
        LogMsg(stderr, "Invalid data found in '%s'\n", TestFile);
        LogClose();
        fclose(fpTestFile);
        exit(1);
      }

      if (ENUM != NUMDE) {
        LogMsg(stderr, "Test file and binary file don't match.\n");
        LogClose();
        fclose(fpTestFile);
        exit(1);
      }

      pleph(JDate, Target, Center, RRD, &bInside);

      if (bInside) {
        if (bFirst) {
          LogMsg(stdout,
            "---JED-----T--C-N--------JPL VALUE--"
            "---------USER VALUE----------RESIDUAL----\n");
          bFirst = FALSE;
        }

        MyValue = RRD[Component-1];
        Residual = MyValue - CValue;
        if (fabs(Residual) >= 0.000000000001) {
          LogMsg(stdout, " *****WARNING: NEXT RESIDUAL >= 1D-12 ***** \n");
        }

        LogMsg(stdout,
          "%8.1f %2d %2d %1d %+19.12f %+19.12f %+19.12f\n",
          JDate, Target, Center, Component, CValue, MyValue, Residual);
      }
    }
  } while (!feof(fpTestFile));

  LogMsg(stdout, "\nOK\n");
  LogClose(); /* Close the log file */
  fclose(fpTestFile); /* Close the test file */
  return (0);
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
        Find next occurance of "/" or "-" by looking at the
        first character of each argument.
    */
    if (argv[i][0] != '-' && argv[i][0] != '/') {
      /* Found an argument that did not begin with "/" or "-" */
      LogMsg(stderr, "Command line arguments must begin with '-' or '/'.\n");
      PrintUsage();
      LogClose();
      remove(szLogFile);
      exit(1);
    }

    switch (argv[i][1]) {
      case '\0':
        LogMsg(stderr, "Space not allowed between - and option.\n");
        LogMsg(stderr, "Type ephtest -h for help.\n");
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
        LogMsg(stderr, "Type ephtest -h for help.\n");
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
  LogMsg(stdout, "                        Written by Joe Heafner\n");
  LogMsg(stdout, "***Tests a binary ephemeris file against JPL test data\n");
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

  printf("Usage: EPHTEST [-i:FILE][-h]\n");
  printf("Valid command line options are:\n");
  printf("\n");
  printf("    -i:FILE   Use FILE as input\n");
  printf("    -h        Display this help screen\n");
  printf("\n");
  printf("Command line options may be in any order.\n");
  printf("\n");
}

/* End Of File - ephtest.c ******************************************/
