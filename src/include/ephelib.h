//
// Created by miquel on 13/6/20.
//

#ifndef CELESTIAL_EPHELIB_H
#define CELESTIAL_EPHELIB_H

#endif //CELESTIAL_EPHELIB_H

#include"astro_constants.h"

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
#include "support.h"

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

#ifndef TRUE
#define TRUE  (1)
#define FALSE (0)
#endif

#ifndef SEEK_SET
#define SEEK_SET  (0)
 #define SEEK_CUR  (1)
 #define SEEK_END  (2)
#endif

/* double precision matrix data type */
typedef double** DMatrix;

/* function prototypes */
void Aberrate(double *p1, double *EBdot, double *p2);
double amodulo(double a, double b);
void Cal2JED(int m, int d, int y, double utc, int s, double tai_utc,
             double ut1_utc, int w, double *jed);
void constants(char *FileName);
void Conway(double uele[], double mu, double jed, double *r_cos_nu,
            double *r_sin_nu);
DMatrix createDMatrix(int n);
double deg(double x);
double dms(double x);
double DRound(double x, int n);
FILE *ephopn(char *FileName);
void Epoch2JED(char *epoch, double *jed);
void Eq2Ecl(double *a, int s, double eps, double *b);
void Eq2Hor(double *a, int s, double *b);
void errprt(int group, char *message);
double fix(double x);
void FmtDms(double x, int n, int m, char *s);
void free_matrix(double **m, int nr, int nc);
void freeDMatrix(DMatrix m, int n);
void FunArgIAU(double jed, double *funarg);
void GeocenObs(double jed, double *obsr_geo);
void GetGST(double jed, int s, double *gst);
void GetInvQMatrix(DMatrix QMatrix, DMatrix InvQMatrix);
int  GetNumde(void);
void GetPrecessParams(double jed1, double jed2, double *zeta,
                      double *z, double *theta, double *zetadot,
                      double *zdot, double *thetadot);
void GetQMatrix(double phi, double phidot, int axis, int s,
                DMatrix QMatrix);
void GetRPNmat(double jed1, double jed2, int rpn, int d,
               DMatrix m3, DMatrix m6);
void GetStateVector(double JD, int TARG, int CENT, int recpol,
                    double StateVector[15][15][2][6]);
void HelEphemeris(double *uelement, double mu, double jed, double *posvel);
void interp(int buff, double *t, int ncf, int ncm, int na, int fl,
            double pv[3][2]);
void Interp1(double *x, double *y, int L, int m, double *arg, double *v);
double Interpol(double *x, double *y, int i, double arg);
void JED2Cal(double jed, int *yr, int *mo, int *dy, double *ti);
void JED2Epoch(double jed, char *s, char *epoch);
double **matrix(int nr, int nc);
void MatXMat(DMatrix a, DMatrix b, DMatrix c, int n);
void MatXVec(DMatrix a, double *b, double *c, int l, int m);
void MRotate(double *vin, int axis, double phi, double *vout);
void GetDpsiDeps(double jed, double *dpsi, double *deps,
                 double *dpsidot, double *depsdot);
void Obliquity(double jed1, double jed2, int m, double *obl,
               double *obldot);
void pleph(double jd, int targ, int cent, double *rrd, int *inside);
void Pol2Rec(double *a, double *b);
void PrecessElements(char *eqnx1, double *element, char *eqnx2);
void QRotate(double *vin, int axis, double phi, double phidot,
             int s, double *vout);
void RayBend(double *earth_hel, double *body_geo, double *body_hel,
             double *p1);
void Rec2Pol(double *a, double *b);
void Reduce (double jed, int body, int place, double StarData[],
             double p3[]);
void Refract(double ra1, double dec1, double lha, double temp,
             double press, double *ra2, double *dec2);
void RotMat(int axis, double phi, DMatrix r, DMatrix dr);
void RST(double jed, double *ra, double *dec, double z0, double deltat,
         char *ris, char *trn, char *set);
void split(double tt, double *ipart, double *fpart);
void SplitStateVector(double *pv, double *p, double *v);
void state(double *jed, int LList[], double pv[6][13], double *nut);
void Stat2Elem(double *posvel, double mu, double jed, double *uelement);
double StumpffN(double x, int Norder);
void Transpose(DMatrix a, DMatrix b, int n);
void Uvector(double *a, double *unita);
void Vcross(double *a, double *b, double *acrossb);
void Vdot(int n, double *a, double *b, double *adotb);
double Vecmag(double *a);
