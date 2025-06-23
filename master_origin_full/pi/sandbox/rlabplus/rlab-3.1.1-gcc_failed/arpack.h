// arpack.h

//
// Defines, and declarations for a
// C => fortran, RLaB => ARPACK interface.
//

#ifndef RLAB_ARPACK_H
#define RLAB_ARPACK_H

#include "fi.h"

#ifdef HAVE_FORTRAN_UND_BACK
  #define DSDRV1   dsdrv1_
  #define DNDRV1   dndrv1_
  #define DNDRV2   dndrv2_
  #define ZNDRV1   zndrv1_
  #define ZNDRV2   zndrv2_
  #define DSEUPD   dseupd_
  #define DSAUPD   dsaupd_
  #define DNEUPD   dneupd_
  #define DNAUPD   dnaupd_
  #define ZNEUPD   zneupd_
  #define ZNAUPD   znaupd_
#endif

#ifdef HAVE_FORTRAN_UND_FRONT
  #define DSDRV1   _dsdrv1
  #define DNDRV1   _dndrv1
  #define DNDRV1   _dndrv2
  #define ZNDRV1   _zndrv1
  #define ZNDRV2   _zndrv2
  #define DSEUPD   _dseupd
  #define DSAUPD   _dsaupd
  #define DNEUPD   _dneupd
  #define DNAUPD   _dnaupd
  #define ZNEUPD   _zneupd
  #define ZNAUPD   _znaupd
#endif

#ifdef HAVE_FORTRAN_UPPERCASE
    // Do nothing, the existing code is OK
#endif

#ifdef HAVE_FORTRAN_LOWERCASE
  #define DSDRV1   dsdrv1
  #define DNDRV1   dndrv1
  #define DNDRV2   dndrv2
  #define ZNDRV1   zndrv1
  #define ZNDRV2   zndrv2
  #define DSEUPD   dseupd
  #define DSAUPD   dsaupd
  #define DNEUPD   dneupd
  #define DNAUPD   dnaupd
  #define ZNEUPD   zneupd
  #define ZNAUPD   znaupd
#endif

//
// ARPACK FORTRAN77 solvers: ar1rl.f
//
extern int DSDRV1 (double *, int *, double *, int *,
                   double *, double *, double *, double *, int *, int *, int *,
                   double *, int *, int *, int *, char *, double *, int *, int *, int *);
extern int DNDRV1 (double *, int *, double *, int *, double *,
                   double *, double *, double *, double *, int *, int *, int *,
                   double *, int *, int *, int *, char *, double *, int *, int *, int *);
extern int DNDRV2 (double *, int *, double *, int *, double *,
                   double *, double *, double *, int *, int *, int *,
                   double *, double *, int *, int *, int *, char *, double *, int *, int *,
                   int *);
extern int ZNDRV1 (Complex *, int *, Complex *, int *, Complex *,
                   Complex *, Complex *, Complex *, Complex *, int *, int *, int *,
                   Complex *, int *, int *, int *, char *, Complex *, int *, int *, int *,
                   double *, double *);
extern int ZNDRV2 (Complex *, Complex *, int *, Complex *,
                   Complex *, Complex *, Complex *, int *, int *, int *,
                   Complex *, Complex *, int *, int *, int *,
                   char *, Complex *, int *, int *, int *,
                   double *, double *);
extern int DSEUPD (int *, char *, int *, double *,
                   double *, int *, double *, char *,
                   int *, char *, int *, double *,
                   double *, int *, double *, int *,
                   int *, int *, double *, double *,
                   int *, int *);
extern int DSAUPD (int *, char *, int *, char *, int *, double *, double *,
                   int *, double *, int *, int *, int *, double *, double *,
                   int *, int *);
extern int DNEUPD (int *, char *, int *, double *, double *, double *, int *,
                   double *, double *, double *, char *, int *, char *, int *, double *,
                   double *, int *, double *, int *, int *, int *, double *, double *,
                   int *, int *);
extern int DNAUPD (int *, char *, int *, char *, int *, double *, double *,
                   int *, double *, int *, int *, int *, double *, double *,
                   int *, int *);
extern int ZNEUPD (int *, char *, int *, Complex *, Complex *, int *, Complex *,
                   Complex *, char *, int *, char *, int *, double *, Complex *, int *,
                   Complex *, int *, int *, int *, Complex *, Complex *, int *,
                   double *, int *);
extern int ZNAUPD (int *, char *, int *, char *, int *, double *, Complex *,
                   int *, Complex *, int *, int *, int *, Complex *, Complex *, int *,
                   double *, int *);
#endif // ARPACK.H
