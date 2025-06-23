// powell's minimization without derivatives solvers

//
// Defines, and declarations for a
// C => fortran, RLaB => Powell solvers
//

#ifndef RLAB_LIBMJDPOWELL_H
#define RLAB_LIBMJDPOWELL_H

#include "fi.h"

#ifdef HAVE_FORTRAN_UND_BACK
  #define NEWUOA   newuoa_
  #define BOBYQA   bobyqa_
  #define LINCOA   lincoa_
  #define COBYLA   cobyla_
#endif

#ifdef HAVE_FORTRAN_UND_FRONT
  #define NEWUOA   _newuoa
  #define BOBYQA   _bobyqa
  #define LINCOA   _lincoa
  #define COBYLA   _cobyla
#endif

#ifdef HAVE_FORTRAN_UPPERCASE
    // Do nothing, the existing code is OK
#endif

#ifdef HAVE_FORTRAN_LOWERCASE
  #define NEWUOA   newuoa
  #define BOBYQA   bobyqa
  #define LINCOA   lincoa
  #define COBYLA   cobyla
#endif

//
// FORTRAN77 solvers: newuoa.f, bobyqa.f, lincoa.f, cobyla.f and their subroutines
//
extern int NEWUOA (int *, int *, double *, double *, double *, int *, int *, double *,
                   int (*)(int * idummy, double * xval, double * f), char * outfile
                  );

extern int BOBYQA (int *, int *, double *, double *, double *, double *, double *, int *, int *, double *,
                   int (*)(int * idummy, double * xval, double * f), char * outfile
                  );

extern int LINCOA (int *, int *, int *, double *, int *, double *, double *,
                   double *, double *, int *, int *, double *,
                   int (*)(int * idummy, double * xval, double * f),
                   char * outfile);

extern int COBYLA (int *, int *, double *, double *, double *,
                   int *, int *, double *, int *, 
                   int (*)(int *, int *, double *, double *, double *),
                   char *);

#endif // NEWUOA.H
