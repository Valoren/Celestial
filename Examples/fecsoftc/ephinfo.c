/*********************************************************************
Project:  ephinfo
Filename: ephinfo.c
Author:   Ported to C by Joe Heafner and Varian Swieter.
Purpose:  Opens a binary ephemeris file and prints the header information.
Thanks to Charles Gamble for extensive modifications.
*********************************************************************/

/* Header Files *****************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include "astrolib.h"
#include "support.h"

/* Function Prototypes **********************************************/
void ParseCmdLine(int argc, char *argv[], char *InFile);
void PrintBanner(void);
void PrintUsage(void);

/* Globals **********************************************************/
char szVersion[] = "EPHINFO v1.060300c";
char szLogFile[] = "ephinfo.log";

/*********************************************************************
Name:    main
Purpose: Main routine for ephinfo.
Inputs:  argc - Number of command-line arguments.
         argv - Pointer to array of command-line arguments.
Outputs: None.
Returns: 0 if execution successful.
Status:  Finished.
Errors:  None known.
*********************************************************************/
int main (int argc, char *argv[]) {

  char InFile[MAX_NAME_SIZE+1] = "";

  /* Open the log file */
  if (LogOpen(szLogFile) == FALSE) {
    fprintf(stderr,
      "Could not open log file '%s': %s\n\n", szLogFile, strerror(errno));
    exit(1);                /* Exit with an error code */
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

  /* We have a filename now. */
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

  constants(InFile);

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

  switch (argv[i][1]) {
    case '\0':
      LogMsg(stderr, "Space not allowed between - and option.\n");
      LogMsg(stderr, "Type ephinfo -h for help.\n");
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
      LogMsg(stderr, "Type ephinfo -h for help.\n");
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
  LogMsg(stdout, "***Opens a binary ephemeris file and prints the header\n");
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

  printf("Usage: EPHINFO [-i:FILE][-h]\n");
  printf("Valid command line options are:\n");
  printf("\n");
  printf("    -i:FILE   Use FILE as input\n");
  printf("    -h        Display this help screen\n");
  printf("\n");
  printf("Command line options may be in any order.\n");
  printf("\n");
}

/* End Of File - ephinfo.c ******************************************/
