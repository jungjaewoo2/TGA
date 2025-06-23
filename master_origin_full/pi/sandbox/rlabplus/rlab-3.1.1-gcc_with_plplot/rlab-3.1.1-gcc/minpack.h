// minpack.h: header file for fortran libraries
//  libminpack.a
// (C) 2006, Marijan Kostrun, project rlabplus
//
// Defines, and declarations for a
// C => fortran, RLaB => MINPACK interface.
//

#ifndef RLAB_MINPACK_H
#define RLAB_MINPACK_H

#include "fi.h"

#ifdef HAVE_FORTRAN_UND_BACK
  #define QRSOLV   qrsolv_
#endif

#ifdef HAVE_FORTRAN_UND_FRONT
  #define QRSOLV   _qrsolv
#endif

#ifdef HAVE_FORTRAN_UPPERCASE
    // Do nothing, the existing code is OK
#endif

#ifdef HAVE_FORTRAN_LOWERCASE
  #define QRSOLV   qrsolv
#endif

//
// MINPACK FORTRAN77 dense solvers: qrsolv
//
extern F_INT QRSOLV ();

//
// rlab/c functions for MINPACK solvers read README
//
int QRSOLV (int *, double *, int *, int *, double *,  double *, double *, double *, double *);

#endif // RLAB_MINPACK_H
