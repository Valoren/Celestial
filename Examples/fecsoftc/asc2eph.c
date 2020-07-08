/*********************************************************************
Project:  asc2eph
Filename: asc2eph.c
Author:   Ported to C by Joe Heafner and Varian Swieter.
Purpose:  This program converts a JPL ephemeris ASCII data file
          into a binary data file. Note that the binary data files
          created with this program are not identical to those
          created with the JPL ephemeris software.
Thanks to Charles Gamble for extensive modifications.
*********************************************************************/

/* Header Files *****************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <errno.h>
#include <sys/types.h>
#include <netinet/in.h>
#include <sys/param.h>
#include "astrolib.h"
#include "support.h"

/* Function Prototypes **********************************************/
void ParseCmdLine(int argc, char *argv[], char *InFile, char *OutFile);
void PrintBanner(void);
void PrintUsage(void);
void ERRPRT(int group, char *message);
void NXTGRP(FILE *fptr, char *outstr);

/* Globals **********************************************************/
char szVersion[] = "ASC2EPH v1.060300c";
char szLogFile[] = "asc2eph.log";

#define MAX_KSIZE       2048
#define MAX_TTL         66
#define MAX_CNAME       7
#define DEFAULT_OUTPUT  "JPLEPH"

int DeletefpLogFile = FALSE,
    PromptForDates  = FALSE,
    HeaderIncluded  = FALSE,
    FIRST           = TRUE,
    KSIZE, NCON, tmpInt;

FILE *fpInputFile, *fpOutputFile, *fpHFile;

short IPT[3][13], LPT[3], NUMDE, tmpShort;

char InFile[MAX_NAME_SIZE+1]  = "",
     OutFile[MAX_NAME_SIZE+1] = "",
     HFile[MAX_NAME_SIZE+1]   = "",
     HEADER[MAX_NAME_SIZE+1]  = "";

double T1, T2, tmpDouble;

char Body[13][4] = { "MER", "VEN", "EMB", "MAR",
                     "JUP", "SAT", "URA", "NEP",
                     "PLU", "MOO", "SUN", "NUT",
                     "LIB" };

/*********************************************************************
Name:    main
Purpose: Main routine for asc2eph.
Inputs:  argc - Number of command-line arguments.
         argv - Pointer to array of command-line arguments.
Outputs: None.
Returns: 0 if execution successful.
Status:  Finished.
Errors:  None known.
*********************************************************************/
int main(int argc, char *argv[]) {

  char Progress[50] = "Searching for first requested record ",
       TTL[3][MAX_TTL] = { "", "", "" },
       CNAM[400][MAX_CNAME], right_buffer[4], line[1024];

  double DB[MAX_KSIZE], SS[3], CVAL[400], AU, EMRAT, DB2Z = 0;

  int DB_size, i, k, jrow, jcol, N, NROUT, NRW, NCOEFF,
      exponent, LeftOver, NumZeros, LengthOfHeader,
      LengthOfRecord, loop;

  long int fpeof, LengthOfFile;

  /* Open the log file */
  if (LogOpen(szLogFile) == FALSE)
    {
      fprintf(stdout, "Could not open log file '%s': %s\n\n",
        szLogFile, strerror(errno));
      exit(1); /* Exit with an error code */
    }

  /* Write a fingerprint to the screen and to log file */
  PrintBanner();

  ParseCmdLine(argc, argv, InFile, OutFile);

  /*
      If you don't want all the data, set T1 and T2 to the begin
      and end times of the span you desire. Units are JED.
  */
  if (PromptForDates) {
    do {
      do {
        fprintf(stdout, "Enter start JED: ");
        fflush(stdout);
        fgets(line, sizeof(line), stdin);
      } while (sscanf(line, "%lf", &T1) != 1);

      do {
        fprintf(stdout, "Enter final JED: ");
        fflush(stdout);
        fgets(line, sizeof(line), stdin);
      } while (sscanf(line, "%lf", &T2) != 1);

      if (T2 < T1) {
        fprintf(stdout, "Final JED must be later than start JED!\n");
      }
    } while (T2 < T1);

    if (T1 == T2) {
      fprintf(stdout,
        "Start and final JED's cannot be the same. Using Defaults.\n");
      T1 = 0.0;
      T2 = 9999999.0;
    }
  } else {
    T1 = 0.0;
    T2 = 9999999.0;
  }

  if (strlen(InFile) == 0) {
    /* Prompt user for file name */
    do {
      fprintf(stdout, "File name to convert (XXXX to end): ");
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

  if (strlen(OutFile) == 0) strcpy(OutFile, DEFAULT_OUTPUT);

  LogMsg(stdout, "Input will be read from %s\n", InFile);
  LogMsg(stdout, "Output will be written to %s\n", OutFile);

  /* Form the header file name */
  if (!HeaderIncluded) {
    right(InFile, 3, right_buffer);
    strcpy(HFile, "header.");
    strcat(HFile, right_buffer);

    if (!fexist(HFile)) {
      LogMsg(stdout, "Header file not present.\n");
      LogMsg(stdout, "\nOK\n");
      LogClose();
      exit(1);
    }
  }

  if (HeaderIncluded) strcpy(HFile, InFile);

  /* Open the header file for input */
  if ((fpHFile = fopen(HFile, "r")) == NULL) {
    LogMsg(stdout, "Header file not present.\n");
    LogMsg(stdout, "\nOK\n");
    LogClose();
    exit(1);
  }

  /*
      Get KSIZE, throwing away 'KSIZE=', 'NCOEFF=', and
      NCOEFF decimal constant.
  */
  if (fscanf(fpHFile, " %*s %d %*s %*d", &KSIZE) != 1) {
    LogMsg(stderr, "Error reading KSIZE\n");
    LogClose();
    fclose(fpHFile);
    exit(1);
  }
  LogMsg(stdout, "KSIZE=  %d\n", KSIZE);

  /*
      Set max size of coeficient array for this particular
      ephemeris file. If KSIZE is an odd number, then
      DB_size = (KSIZE / 2) + (KSIZE % 2).
  */
  DB_size = KSIZE / 2;

  if (DB_size > MAX_KSIZE) {
    /* Array is not large enough */
    LogMsg(stdout,
      "Adjust #define MAX_KSIZE to larger value & recompile.\n");
    LogMsg(stdout, "\nOK\n");
    LogClose();
    exit(1);
  }

  /* Initialize DB through DB_size */
  NXTGRP(fpHFile, HEADER);
  if (strcmp(HEADER, "GROUP   1010") != 0) {
    ERRPRT(1010, "NOT HEADER");
  }

  for (i = 0; i < 3; i++) {
    fscanf(fpHFile, " %65[^\n]", TTL[i]);
    LogMsg(stdout, "%s\n",   TTL[i]);
  }

  /* Read start, end and record span (GROUP 1030) */
  NXTGRP(fpHFile, HEADER);
  if (strcmp(HEADER, "GROUP   1030") != 0) {
    ERRPRT(1030, "NOT HEADER");
  }

  /* Read in values of ss[0], ss[1], ss[2] */
  for (i = 0; i < 3; i++) {
    if (fscanf(fpHFile, " %lf", &SS[i]) != 1) {
      LogMsg(stderr, "Error reading SS[%d].\n", i);
      LogClose();
      fclose(fpHFile);
      exit(1);
    }
  }

  /* Read number of constants and names of constants (GROUP 1040/4) */
  NXTGRP(fpHFile, HEADER);
  if (strcmp(HEADER, "GROUP   1040") != 0) {
    ERRPRT(1040, "NOT HEADER");
  }

  if (fscanf(fpHFile, " %d", &N) != 1) {
    LogMsg(stderr, "Error reading N.\n");
    LogClose();
    fclose(fpHFile);
    exit(1);
  }

  /* Now parse the constant names from each line of input */
  for (i = 0; i < N; i++) {
    if (fscanf(fpHFile, " %s", CNAM[i]) != 1) {
      LogMsg(stderr, "Error reading CNAM[%d].\n", i);
      LogClose();
      fclose(fpHFile);
      exit(1);
    }
  }

  NCON = N;

  /* Read number of values and values (GROUP 1041/4) */
  NXTGRP(fpHFile, HEADER);
  if (strcmp(HEADER, "GROUP   1041") != 0) {
    ERRPRT(1041, "NOT HEADER");
  }

  if (fscanf(fpHFile, " %d", &N) != 1) {
    LogMsg(stderr, "Error reading N.\n");
    LogClose();
    fclose(fpHFile);
    exit(1);
  }

  LogMsg(stdout, "\n"
    "Ephemeris Constants\n"
    "-------------------\n");

  for (i = 0; i < N; i++) {
    /*
        Read cval mask out D and read exponent
        then convert cval to reflect exponent.
    */
    if (fscanf(fpHFile, " %lfD%d", &CVAL[i], &exponent) != 2) {
      LogMsg(stderr, "Error reading CVAL[%d] and exponent.\n", i);
      LogClose();
      fclose(fpHFile);
      exit(1);
    }

    CVAL[i] *= pow(10, exponent);
    if (strcmp(CNAM[i], "AU")    == 0) AU    = CVAL[i];
    if (strcmp(CNAM[i], "EMRAT") == 0) EMRAT = CVAL[i];
    if (strcmp(CNAM[i], "DENUM") == 0) NUMDE = (int) CVAL[i];

    LogMsg(stdout, "%-6s   %+.15E\t", CNAM[i], CVAL[i]);

    if (i % 2) {
      LogMsg(stdout, "\n");
    }
  }

  LeftOver = N % 3;
  switch (LeftOver) {
    case 0:
            NumZeros = 0;
            break;
    case 1:
    case 2:
            NumZeros = 3 - LeftOver;
            break;
  }

  for (i = 0; i < LeftOver; i++) {
    /* Throw away padded values */
    fscanf(fpHFile, " %*fD%*d");
  }

  NXTGRP(fpHFile, HEADER);
  if (strcmp(HEADER, "GROUP   1050") != 0) {
    ERRPRT (1050, "NOT HEADER");
  }

  /* Read pointers from file */
  for (jrow = 0; jrow < 3; jrow++) {
    for (jcol = 0; jcol < 13; jcol++) {
      if (fscanf(fpHFile, " %hd", &IPT[jrow][jcol]) != 1) {
        LogMsg(stderr, "Error reading IPT[%d][%d].\n", jrow, jcol);
        LogClose();
        fclose(fpHFile);
        exit(1);
      }
    }
  }

  LPT[0] = IPT[0][12];
  LPT[1] = IPT[1][12];
  LPT[2] = IPT[2][12];

  LogMsg(stdout,
    "\n\n"
    "  Body    1st coeff      coefs/component      sets of coefs\n"
    "-------------------------------------------------------------\n");

  for (jcol = 0; jcol < 13; jcol++) {
    LogMsg(stdout,
      " %2d %3s     %3d               %3d                  %3d\n",
    jcol+1, Body[jcol], IPT[0][jcol], IPT[1][jcol], IPT[2][jcol]);
  }
  LogMsg(stdout, "\n");

  /* Open direct-access output file ('JPLEPH' by default) */
  if (strlen(OutFile) == 0) {
    if (fexist(DEFAULT_OUTPUT)) remove(DEFAULT_OUTPUT);
      strcpy(OutFile, DEFAULT_OUTPUT);
    }

    if ((fpOutputFile = fopen(OutFile, "wb")) == NULL) {
      LogMsg(stdout, "Can't create binary data file.\n");
      LogMsg(stdout, "\nOK\n");
      LogClose();
      fclose(fpHFile);
      exit(1);
    }

    /******  BEGINNING OF HEADER  *******/
    for (k = 0; k < 3; k++)
      fprintf (fpOutputFile, "%-65s", TTL[k]); /* 195 bytes */

    tmpShort = (short) NCON;
    make_little_endian((char *)&tmpShort, sizeof(short));
    fwrite(&tmpShort, sizeof(short), 1, fpOutputFile); /* 2 bytes */

    for (k = 0; k < NCON; k++)
      fprintf(fpOutputFile, "%-6s", CNAM[k]); /* NCON*6 bytes */

    for (loop = 0; loop < 3; loop++) {
      tmpDouble = (double) SS[loop];
      make_little_endian((char *)&tmpDouble, sizeof(double));
      fwrite(&tmpDouble, sizeof(double), 1, fpOutputFile);
    } /* 24 bytes */

    tmpDouble = (double) AU;
    make_little_endian((char *)&tmpDouble, sizeof(double));
    fwrite(&tmpDouble, sizeof(double), 1, fpOutputFile); /* 8 bytes */

    tmpDouble = (double) EMRAT;
    make_little_endian((char *)&tmpDouble, sizeof(double));
    fwrite(&tmpDouble, sizeof(double), 1, fpOutputFile); /* 8 bytes */

    for (k = 0; k < 3; k++) {
      for (loop = 0; loop < 12; loop++) {
        tmpShort = (short)IPT[k][loop];
        make_little_endian((char *)&tmpShort, sizeof(short));
        fwrite(&tmpShort, sizeof(short), 1, fpOutputFile);
      } /* 72 bytes */
    }

    tmpShort = (short) NUMDE;
    make_little_endian((char *)&tmpShort, sizeof(short));
    fwrite(&tmpShort, sizeof(short), 1, fpOutputFile); /* 2 bytes */

    for (loop = 0; loop < 3; loop++) {
      tmpShort = (short) LPT[loop];
      make_little_endian((char *)&tmpShort, sizeof(short));
      fwrite(&tmpShort, sizeof(short), 1, fpOutputFile);
    } /* 6 bytes */

    for (loop = 0; loop < NCON; loop++) {
      tmpDouble = (double) CVAL[loop];
      make_little_endian((char *)&tmpDouble, sizeof(double));
      fwrite(&tmpDouble, sizeof(double), 1, fpOutputFile);
    } /* NCON*8 bytes */

    /* Length of header = 317+NCON*14 bytes */
    /*******  END OF HEADER  *******/

    /* Read and write the ephemeris data records (GROUP 1070) */
    NXTGRP(fpHFile, HEADER);
    if (strcmp(HEADER, "GROUP   1070") != 0) {
      ERRPRT(1070, "NOT HEADER");
    }

    /*
        Close the header file and open the actual data file
        if the header is in a seperate file.
    */
    if (!HeaderIncluded) {
      fclose(fpHFile);
      if ((fpInputFile = fopen(InFile, "rt")) == NULL) {
        LogMsg(stdout, "Can't open ascii data file '%s'.\n", InFile);
        LogMsg(stdout, "\nOK\n");
        LogClose();
        exit(1);
      }
    }

    NROUT = 0;

    if (HeaderIncluded) fpInputFile = fpHFile;

    /* Read the very first record in */
    if (!feof(fpInputFile)) {
      if (fscanf(fpInputFile, " %d %d", &NRW, &NCOEFF) != 2) {
        LogMsg(stderr, "Error reading NRW and NCOEFF.\n");
        LogClose();
        fclose(fpInputFile);
        exit(1);
       }

       for (k = 0; k < NCOEFF; k++) {
         /* Read DB mask out D and read exponent */
         /* and convert DB to reflect exponent.  */
         if (fscanf(fpInputFile, " %lfD%d", &DB[k], &exponent) != 2) {
           LogMsg(stderr, "Error reading DB[%d] and exponent.\n", k);
           LogClose();
           fclose(fpInputFile);
           exit(1);
         }
         DB[k] *= pow(10, exponent);
       }

       LeftOver = NCOEFF % 3;
       switch (LeftOver) {
         case 0:
                  NumZeros = 0;
                  break;
         case 1:
         case 2:
                  NumZeros = 3 - LeftOver;
                  break;
       }

       /* Read in padded values and discard them */
       for (k = 0; k < NumZeros; k++) {
         fscanf(fpInputFile, " %*fD%*d");
       }
    }

    while (!feof(fpInputFile) && (DB[1] < T2)) {
      if ((2 * NCOEFF) != KSIZE) {
        ERRPRT(NCOEFF, " 2*NCOEFF not equal to KSIZE");
      }

      /*
          Skip this data block if the end of the interval is less
          than the specified start time or if it does not begin
          where the previous block ended.
      */
      if ((DB[1] > T1) && (DB[0] >= DB2Z)) {
        if (FIRST) {
          /*
              Don't worry about the intervals overlapping or abutting
              if this is the first applicable interval.
          */
          DB2Z = DB[0];
          FIRST = FALSE;
        }
        if (DB[0] != DB2Z) {
          /*
              Beginning of current interval is past the end
              of the previous one.
          */
          ERRPRT(NRW, "Records do not overlap or abut.");
        }

        DB2Z = DB[1];
        NROUT++;

        /* Write the numbers to binary file */
        for (loop = 0; loop < NCOEFF; loop++) {
          tmpDouble = (double) DB[loop];
          make_little_endian((char *)&tmpDouble, sizeof(double));
          fwrite(&tmpDouble, sizeof(double), 1, fpOutputFile);
        }

        /*
            Save this block's starting date, it's interval span,
            and its end date.
        */
        if (NROUT == 1) {
          SS[0] = DB[0];
          SS[2] = DB[1] - DB[0];
        }
        SS[1] = DB[1];

        /* Update the user as to our progress every 10th block */
        if ((NROUT % 10) == 1) {
          if (DB[0] >= T1) {
            fprintf(stdout,
              "  %4d EPHEMERIS RECORD(S) WRITTEN.  LAST JED = %9.1f\r",
              NROUT, DB[1]);
          }
        }
      } else {
        fprintf(stdout, "%s\r", Progress);
        strcat(Progress, ".");
        if (strlen(Progress) == 47) {
          strcpy(Progress, "Searching for first requested record ");
          fprintf(stdout,
            "Searching for first requested record           \r");
        }
      }

      /*
          Read next block of coefficients unless EOF has
          been reached.
      */
      if (!feof(fpInputFile)) {
        if (fscanf(fpInputFile, " %d %d", &NRW, &NCOEFF) != 2) {
          LogMsg(stderr, "Error reading NRW and NCOEFF.\n");
          LogClose();
          fclose(fpInputFile);
          exit(1);
        }

        for (k = 0; k < NCOEFF; k++) {
          /*
              Read DB mask out D and read exponent
              and convert DB to reflect exponent.
          */
          if (fscanf(fpInputFile, " %lfD%d", &DB[k], &exponent) != 2) {
            LogMsg(stderr, "Error reading DB[%d] and exponent.\n", k);
            LogClose();
            fclose(fpInputFile);
            exit(1);
          }
          DB[k] *= pow(10, exponent);
        }
        LeftOver = NCOEFF % 3;
        switch (LeftOver) {
          case 0:
                   NumZeros = 0;
                   break;
          case 1:
          case 2:
                   NumZeros = 3 - LeftOver;
                   break;
        }

        /* Read in padded values and discard them */
        for (k = 0; k < NumZeros; k++) {
          fscanf(fpInputFile, " %*fD%*d "); /* Need trailing space */
        }                                   /* to test eof.        */
      }
    }

    /*
        End of file, but no records yet written OR
        just no records yet written.
    */
    if ((feof(fpInputFile) && (NROUT == 0)) || (NROUT == 0)) {
      NROUT++;
      SS[0] = DB[0];
      SS[1] = DB[1];
      for (loop = 0; loop < NCOEFF; loop++) {
        tmpDouble = (double) DB[loop];
        make_little_endian((char *)&tmpDouble, sizeof(double));
        fwrite(&tmpDouble, sizeof(double), 1, fpOutputFile);
      }
    }
    /* End of file but T2 lies within most recently read record */
    else if (feof(fpInputFile)) {
      NROUT++;
      SS[1] = DB[1];
      for (loop = 0; loop < NCOEFF; loop++) {
        tmpDouble = (double) DB[loop];
        make_little_endian((char *)&tmpDouble, sizeof(double));
        fwrite(&tmpDouble, sizeof(double), 1, fpOutputFile);
      }
    }
    else if (DB[1] > T2) {
      NROUT++;
      SS[1] = DB[1];
      for (loop = 0; loop < NCOEFF; loop++) {
        tmpDouble = (double) DB[loop];
        make_little_endian((char *)&tmpDouble, sizeof(double));
        fwrite(&tmpDouble, sizeof(double), 1, fpOutputFile);
      }
    }

    /* Reset file pointer from beginning of file */
    /* and write start and final JED to header.  */
    fpeof = ftell(fpOutputFile);
    fseek(fpOutputFile, (197 + NCON * 6), SEEK_SET);

    for (loop = 0; loop < 2; loop++) {
      tmpDouble = (double) SS[loop];
      make_little_endian((char *)&tmpDouble, sizeof(double));
      fwrite(&tmpDouble, sizeof(double), 1, fpOutputFile);
    }

    LengthOfHeader = 317 + NCON * 14;
    LengthOfRecord =   8 * NCOEFF;
    LengthOfFile = (long)LengthOfHeader + (long)LengthOfRecord * (long)NROUT;
    /* LengthOfFile should equal fpeof */

    LogMsg(stdout,
      " %4d EPHEMERIS RECORD(S) WRITTEN.  LAST JED = %9.1lf\n", NROUT, DB[1]);
    LogMsg(stdout, "\nStart JED for this file is %9.1lf\n", SS[0]);
    LogMsg(stdout, "Final JED for this file is %9.1lf\n\n", SS[1]);
    LogMsg(stdout, "Header is %d bytes\n", LengthOfHeader);
    LogMsg(stdout, "Record is %d bytes\n", LengthOfRecord);
    LogMsg(stdout, " %ld bytes written\n", LengthOfFile);
    LogMsg(stdout, " %ld bytes actual file size\n", fpeof);

    LogMsg(stdout, "\nOK\n");

    /* Close all files */
    LogClose();
    fclose(fpOutputFile);
    fclose(fpInputFile);

    if (DeletefpLogFile) {
      remove(szLogFile);
    }

    return(0);
}

/*********************************************************************
Name:    ParseCmdLine
Purpose: Routine to parse the command line.
Inputs:  argc - Number of command-line arguments.
         argv - Pointer to array of command-line arguments.
Outputs: InFile  - Filename given by user for input file.
         OutFile - Filename given by user for output file.
Returns: Nothing.
Status:  Finished.
Errors:  None known.
*********************************************************************/
void ParseCmdLine(int argc, char *argv[], char *InFile, char *OutFile) {
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
      LogMsg(stdout, "Command line arguments must begin with '-' or '/'.\n");
      PrintUsage();
      LogClose();
      remove(szLogFile);
      exit(1);
    }

    switch (argv[i][1]) {
      case '\0':
                 LogMsg(stdout, "Space not allowed between - and option.\n");
                 LogMsg(stdout, "Type asc2eph -h for help.\n");
                 LogClose();
                 remove(szLogFile);
                 exit(1);
      case 'I':
      case 'i':
                 if (strlen(argv[i]) == 2) {
                   /* We were given nothing after the "-I"
                      so return empty string */
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
      case 'O':
      case 'o':
                 if (strlen(argv[i]) == 2) {
                   /* We were given nothing after the "-O"
                      so return empty string */
                   strcpy(OutFile, "");
                 } else {
                   if (argv[i][2] != ':') {
                     /* User missed out colon so return empty string */
                     strcpy(OutFile, "");
                   } else {
                     strcpy(OutFile, &(argv[i][3]));
                   }
                 }

                 /* If OutFile is an empty string then return "JPLEPH"
                    as a default */
                 if (strlen(OutFile) == 0) {
                   strcpy(OutFile, DEFAULT_OUTPUT);
                 }
                 break;
      case '?':  /* Using this for help will cause problems under UNIX */
                 /* because it is used by the shell for substitution.  */
      case 'h':
                 PrintUsage();
                 LogClose();
                 remove(szLogFile);
                 exit(0);
      case 'H':
                 LogMsg(stdout ,   "Assuming included header.\n");
                 HeaderIncluded = TRUE;
                 break;
      case 'D':
                 PromptForDates = TRUE;
                 break;
      case 'K':
                 DeletefpLogFile = TRUE;
                 break;
      case 'v':  /* Undocumented command line option */
                 /* to print version information.    */
                 LogMsg(stdout, "OK\n");
                 LogClose();
                 remove(szLogFile);
                 exit(0);
      default:
                 LogMsg(stdout, "Option not recognized.\n");
                 LogMsg(stdout, "Type asc2eph -h for help.\n");
                 LogClose();
                 remove(szLogFile);
                 exit(1);
    }
  }
}

/*********************************************************************
Name:    NXTGRP
Purpose: Finds the next data group in an ephemeris file.
Inputs:  fptr - File pointer to read from.
Outputs: outstr - Destination for the output string.
Returns: Nothing.
Status:  Finished.
Errors:  None known.
*********************************************************************/
void NXTGRP(FILE *fptr, char *outstr) {

  char string[50];
  char str[6];
  char num[5];

  fscanf(fptr, " %s %s", str, num);
  strcpy(string, str);
  strcat(string, "   ");
  strcat(string, num);

  strcpy(outstr, string);
}

/*********************************************************************
Name:    ERRPRT
Purpose: Prints error messages to stdout and log file.
Inputs:  group   - Error number.
         message - Error message.
Outputs: None.
Returns: Nothing.
Status:  Finished.
Errors:  None known.
*********************************************************************/
void ERRPRT(int group, char *message) {

  LogMsg(stdout, "\nERROR #%d %s\n", group, message);
  LogClose();
  fclose(fpHFile);
  fclose(fpOutputFile);
  exit(1);
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
  LogMsg(stdout, "***Converts JPL ASCII ephemeris data files to binary\n");
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

  printf("Usage: ASC2EPH [-i:FILE][-o:FILE][-D][-H][-K][-h]\n");
  printf("Valid command line options are:\n");
  printf("\n");
  printf("    -i:FILE   Use FILE as input assuming separate header file\n");
  printf("              Will prompt if this option is not used\n");
  printf("    -o:FILE   Use FILE as name of binary file\n");
  printf("              Default binary file name is %s\n", DEFAULT_OUTPUT);
  printf("    -D        Program will prompt for initial and final dates\n");
  printf("              Defaults are full range of data file\n");
  printf("    -H        Assume data file has header included\n");
  printf("              Default is separate header in HEADER.XXX\n");
  printf("    -K        Delete log file ASC2EPH.LOG\n");
  printf("              Log file kept by default\n");
  printf("    -h        Display this help screen\n");
  printf("\n");
  printf("Command line options may be in any order.\n");
}
/* End Of File - asc2eph.c ******************************************/
