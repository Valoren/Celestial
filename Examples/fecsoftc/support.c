/*********************************************************************
Project:        fecsoftc
Filename:       support.c
Author:         Joe Heafner
Purpose:        General support routines.
Thanks to Charles Gamble for extensive modifications.
*********************************************************************/

/* Header Files *****************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <stdarg.h>
#include <sys/types.h>
#include <netinet/in.h>
#include <sys/param.h>
#include "support.h"

/* Globals **********************************************************/
static FILE *fpLogFile = NULL;

/*********************************************************************
Name:    make_little_endian
Purpose: Reverses the byte ordering of the given block of memory IF this
         machine is BIG ENDIAN. Used to ensure that binary data written
         out to file is little-endian in nature.
Inputs:  ptr - Pointer to block of memory to reverse.
         len - Length of block of memory to reverse.
Outputs: ptr - Byte-reversed block of memory.
Returns: Nothing.
Status:  Finished.
Errors:  None known.
*********************************************************************/
void make_little_endian(char *ptr, int len) {

  if (BIG_ENDIAN_TEST) {
    reverse_bytes(ptr, len);
  }
}

/*********************************************************************
Name:    convert_little_endian
Purpose: Reverses the byte ordering of the given block of memory IF this
         machine is BIG ENDIAN. Used to convert the little-endian binary
         data read in from a file into the native endian format for this
         machine.
Inputs:  ptr - Pointer to block of memory to reverse.
         len - Length of block of memory to reverse.
Outputs: ptr - Byte-reversed block of memory.
Returns: Nothing.
Status:  Finished.
Errors:  None known.
*********************************************************************/
void convert_little_endian(char *ptr, int len) {

  if (BIG_ENDIAN_TEST) {
    reverse_bytes(ptr, len);
  }
}

/*********************************************************************
Name:    reverse_bytes
Purpose: Reverses the byte ordering of the given block of memory.
Inputs:  ptr - Pointer to block of memory to reverse.
         len - Length of block of memory to reverse.
Outputs: ptr - Byte-reversed block of memory.
Returns: Nothing.
Status:  Finished.
Errors:  None known.
*********************************************************************/
void reverse_bytes(char *ptr, int len) {

  int i;
  char *tmp;

  if ((tmp = malloc(len)) == NULL) {
    fprintf(stderr, "Memory allocation error (reverse_bytes:%d)\n", len);
    exit(1);
  }

  for (i = 0; i < len; i++) {
      tmp[i] = ptr[(len-1)-i];
  }

  memcpy(ptr, tmp, len);
  free(tmp);
}

/*********************************************************************
Name:    ucase
Purpose: Converts string to uppercase. This emulates the BASIC function
         ucase(a$).
Inputs:  str - String to convert.
Outputs: str - Converted string.
Returns: Nothing.
Status:  Finished.
Errors:  None known.
*********************************************************************/
void ucase(char str[]) {

  int i;
  int len = strlen(str);

  for (i = 0; i < len; i++)
    str[i] = toupper(str[i]);
}

/*********************************************************************
Name:    left
Purpose: Copies the first n characters of string str into dest.
Inputs:  str  - String to copy.
         n    - Number of characters to copy.
Outputs: dest - Copied string.
Returns: Nothing.
Status:  Finished.
Errors:  None known.
*********************************************************************/
void left(char str[], int n, char dest[]) {

  /* Safety test - make sure n is >= 0 */
  if (n > 0) {
    strncpy(dest, str, n);
    dest[n] = '\0'; /* Terminate just in case strlen(str) >= n */
  } else {
    strcpy(dest, "");
  }
}

/*********************************************************************
Name:    right
Purpose: Copies the last n characters of string str into dest.
Inputs:  str  - String to copy.
         n    - Number of characters to copy.
Outputs: dest - Copied string.
Returns: Nothing.
Status:  Finished.
Errors:  None known.
*********************************************************************/
void right(char str[], int n, char dest[]) {

  int len;

  if (n > 0) {
    len = strlen(str);
    if (n >= len) {
      /* Just copy whole string */
      strcpy(dest, str);
    } else {
      /* Copy section of string */
      strcpy(dest, &str[len-n]);
    }
  } else {
    strcpy(dest, "");
  }
}

/*********************************************************************
Name:    Trim
Purpose: Removes whitespace from either end of string.
Inputs:  str - String to trim.
Outputs: str - Trimmed string.
Returns: Nothing.
Status:  Finished.
Errors:  None known.
*********************************************************************/
void Trim(char str[]) {

  RightTrim(str); /* Remove trailing whitespace */
  LeftTrim(str);  /* Remove preceeding whitespace */
}

/********************************************************************
Name:    RightTrim
Purpose: Removes trailing whitespace from string.
Inputs:  str - String to trim.
Outputs: str - Trimmed string.
Returns: Nothing.
Status:  Finished.
Errors:  None known.
*********************************************************************/
void RightTrim(char str[]) {

  int i;

  i = strlen(str) - 1;
  while (i >= 0 && isspace(str[i])) {
    str[i] = '\0';
    i--;
  }
}

/*********************************************************************
Name:    LeftTrim
Purpose: Removes preceeding whitespace from string.
Inputs:  str - String to trim.
Outputs: str - Trimmed string.
Returns: Nothing.
Status:  Finished.
Errors:  None known.
*********************************************************************/
void LeftTrim(char str[]) {

  unsigned int i;

  i = 0;
  while (i < strlen(str) && isspace(str[i]))
    i++;

  if (i)
   strcpy(str, &str[i]);
}

/*********************************************************************
Name:    fexist
Purpose: Tests to see if a file exists.
Inputs:  filename - Filename to test.
Outputs: None.
Returns: 0 if the file does not exist.
         1 if the file exists.
Status:  Finished.
Errors:  None known.
*********************************************************************/
int fexist(char *filename) {

  FILE *fp;

  if ((fp = fopen(filename,"r")) == NULL) {
    /* File could not be opened */
    return(0);
  } else {
    /* File was opened successfully */
    fclose(fp);
    return(1);
  }
}

/*********************************************************************
Name:    LogOpen
Purpose: Open a log file with the specified name.
Inputs:  filename - Name of log file to open.
Outputs: None.
Returns: TRUE if successful, FALSE otherwise.
Status:  Finished.
Errors:  None known.
*********************************************************************/
int LogOpen(char *filename) {

  /* Check to see if file is already open */
  if (fpLogFile)
    return (TRUE);

  /* Open the log file */
  if ((fpLogFile = fopen(filename, "wt")) == NULL)
    return (FALSE);
  else
    return (TRUE);
}

/*********************************************************************
Name:    LogClose
Purpose: Closes the log file.
Inputs:  None.
Outputs: None.
Returns: Nothing.
Status:  Finished.
Errors:  None known.
*********************************************************************/
void LogClose() {

  /* Check to see if file is open */
  if (fpLogFile) {
    fclose(fpLogFile);
    fpLogFile = NULL;
  }
}

/*********************************************************************
Name:    LogMsg
Purpose: Prints a message to a given file pointer and to the log file.
         Log file must have been opened previously with LogOpen().
Inputs:  fptr     - File pointer to print with.
         format   - Format string (like printf).
         va_alist - Variable length argument list.
Outputs: None.
Returns: Nothing.
Status:  Finished.
Errors:  None known.
*********************************************************************/
void LogMsg(FILE* fptr, const char *format, ...) {

  char buffer[1024];
  va_list argptr;

  va_start(argptr, format);
  vsprintf(buffer, format, argptr);
  va_end(argptr);

  /* Write to the log file if it is open */
  if (fpLogFile)
    fprintf(fpLogFile, "%s", buffer);

  /* Write to the given file pointer if it is valid */
  if (fptr)
    fprintf(fptr, "%s", buffer);
}

/* End Of File - support.c ******************************************/
