#ifndef _SUPPORT_
#define _SUPPORT_

/*********************************************************************
Project:  fecsoftc
Filename: support.h
Author:   Joe Heafner
Purpose:  General support routines header file.
Thanks to Charles Gamble for extensive modifications.
*********************************************************************/

#include <sys/types.h>
#include <netinet/in.h>
#include <stdarg.h>
#include <sys/param.h>
/* #include <386/endian.h> */

#define MAX_NAME_SIZE       255
#define MAX_EXTENSION_SIZE  5
#define BIG_ENDIAN_TEST     (htonl(1) == 1)

#ifndef TRUE
 #define TRUE   (1)
 #define FALSE  (0)
#endif

/* Function prototypes */
void make_little_endian(char *ptr, int len);
void convert_little_endian(char *ptr, int len);
void reverse_bytes(char *ptr, int len);
void ucase(char str[]);
void left(char str[], int n, char dest[]);
void right(char str[], int n, char dest[]);
void Trim(char str[]);
void RightTrim(char str[]);
void LeftTrim(char str[]);
int  fexist(char *filename);
int  LogOpen(char *filename);
void LogClose(void);
void LogMsg(FILE *fptr, const char *format, ...);

#endif /* _SUPPORT_ */
