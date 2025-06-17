// slatec.h: header file for fortran libraries
//  libslatec.a
// (C) 2006, Marijan Kostrun, project rlabplus
//
// Defines, and declarations for a
// C => fortran, RLaB => SLATEC interface.
//

#ifndef RLAB_SLATEC_H
#define RLAB_SLATEC_H

#include "fi.h"

#ifdef HAVE_FORTRAN_UND_BACK
  #define DPOLFT   dpolft_
  #define DPCOEF   dpcoef_
  #define DLSEI    dlsei_
  #define DAVINT   davint_
#endif

#ifdef HAVE_FORTRAN_UND_FRONT
  #define DPOLFT   _dpolft
  #define DPCOEF   _dpcoef
  #define DLSEI    _dlsei
  #define DAVINT   _davint
#endif

#ifdef HAVE_FORTRAN_UPPERCASE
    // Do nothing, the existing code is OK
#endif

#ifdef HAVE_FORTRAN_LOWERCASE
  #define DPOLFT   dpolft
  #define DPCOEF   dpcoef
  #define DLSEI    dlsei
  #define DAVINT   davint
#endif

//
// SLATEC FORTRAN77 dense solvers: dpolft, dpcoef, dlsei
//
extern int DPOLFT (int *, double *, double *, double *, int *, int *, double *,
                   double *, int *, double *);
extern int DPCOEF (int *, double *, double *, double *);
extern int DLSEI  ();
extern int DAVINT (double *, double *, int *, double *, double *, double *, int *);

#endif // RLAB_SLATEC_H
