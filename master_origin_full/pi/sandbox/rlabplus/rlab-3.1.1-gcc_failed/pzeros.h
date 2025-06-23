// pzeros.h

//
// Defines, and declarations for a
// C => fortran, RLaB => PZEROS.F interface.
//

#ifndef RLAB_PZEROS_H
#define RLAB_PZEROS_H

#include "fi.h"

#ifdef HAVE_FORTRAN_UND_BACK
  #define POLZEROS   polzeros_
  #define CPSC       cpsc_
  #define CPS1FUNC   cps1func_
#endif

#ifdef HAVE_FORTRAN_UND_FRONT
  #define POLZEROS   _polzeros
  #define CPSC       _cpsc
  #define CPS1FUNC   _cps1func
#endif

#ifdef HAVE_FORTRAN_UPPERCASE
    // Do nothing, the existing code is OK
#endif

#ifdef HAVE_FORTRAN_LOWERCASE
  #define POLZEROS   polzeros
  #define CPSC       cpsc
  #define CPS1FUNC   cps1func
#endif

//
// FORTRAN77 solvers: polzeros.f, 579.f
//
extern int POLZEROS (int *, Complex *, double *, double *, double *,
                     int *, Complex *, double *, int *, int *,
                     double *, double *);
extern int CPSC (Complex *, int *, int *, double *, Complex *, double *, double *);

int CPS1FUNC (Complex *, Complex *);

#endif // PZEROS.H
