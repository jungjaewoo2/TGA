// sp.h

//
// Defines, and declarations for a
// C => fortran, RLaB => signal processing routines and packages
//

#ifndef RLAB_SP_H
#define RLAB_SP_H

#include "fi.h"

#ifdef HAVE_FORTRAN_UND_BACK
    #define GCVSPL gcvspl_
    #define SPLDER splder_
    #define STL2   stl2_
    #endif

#ifdef HAVE_FORTRAN_UND_FRONT
    #define GCVSPL _gcvspl
    #define SPLDER _splder
    #define STL2   _stl2
    #endif

#ifdef HAVE_FORTRAN_UPPERCASE
    // Do nothing, the existing code is OK
#endif

#ifdef HAVE_FORTRAN_LOWERCASE
    #define GCVSPL gcvspl
    #define SPLDER splder
    #define STL2   stl2
    #endif

//
// FORTRAN77 solvers
//
extern int GCVSPL (double *, double *, int *, double *, double *, int *,
                   int *, int *, int *, double *, double *, int *, double *, int *);
extern double SPLDER (int *, int *, int *, double *, double *, double *, int *, double *);
extern int STL2 (double *, double *, double *, int *, double *, double *, double *, int *, int *);

//
// Global variables for all signal processing routines
//

#endif // RLAB_SP_H
