// dierckx.h

//
// Defines, and declarations for a
// C => fortran, RLaB => libfit
//

#ifndef RLAB_DIERCKX_H
#define RLAB_DIERCKX_H

#include "fi.h"

#ifdef HAVE_FORTRAN_UND_BACK
  #define CURFIT   curfit_
  #define PERCUR   percur_
  #define CONCUR   concur_
  #define CONCON   concon_
  #define COCOSP   cocosp_
  #define SURFIT   surfit_
  #define REGRID   regrid_
  #define DXSPLDER dxsplder_
  #define DXPARDER dxparder_
#endif

#ifdef HAVE_FORTRAN_UND_FRONT
  #define CURFIT   _curfit
  #define PERCUR   _percur
  #define CONCUR   _concur
  #define CONCON   _concon
  #define COCOSP   _cocosp
  #define SURFIT   _surfit
  #define REGRID   _regrid
  #define DXSPLDER _dxsplder
  #define DXPARDER _dxparder
#endif

#ifdef HAVE_FORTRAN_UPPERCASE
  // Do nothing, the existing code is OK
#endif

#ifdef HAVE_FORTRAN_LOWERCASE
  #define CURFIT   curfit
  #define PERCUR   percur
  #define CONCUR   concur
  #define CONCON   concon
  #define COCOSP   cocosp
  #define SURFIT   surfit
  #define REGRID   regrid
  #define DXSPLDER dxsplder
  #define DXPARDER dxparder
#endif

//
// FORTRAN77 solvers
//
extern int CURFIT (int *, int *, double *, double *, double *, double *, double *,
                   int *, double *, int *, int *, double *, double *, double *,
                   double *, int *, int *, int *, int *, double *);
extern int PERCUR (int *, int *, double *, double *, double *, int *, double *, int *,
                   int *, double *, double *, double *, double *, int *, int *, int *,
                   int *, double *);
extern int CONCUR (int *, int *, int *, double *, int *, double *, double *, double *,
                   int *, double *, int *, int *, double *, int *,
                   int *, double *, int *, int *, double *, int *, double *,
                   int *, double *, double *, double *, int *, int *, int *,
                   int *, double *);
extern int CONCON (int *, int *, double *, double *, double *, double *, double *, int *,
                   int *, int *, int *, double *, double *, double *, double *,
                   int *, double *, int *, int *, int *, int *);
extern int COCOSP (int *, double *, double *, double *, int *, double *, double *,
                   int *, int *, double *, double *, double *, double *,
                   double *, int *, int *, int *, int *);
extern int SURFIT (int *, int *, double *, double *, double *, double *, double *, double *,
                   double *, double *, int *, int *, double *, int *, int *, int *, double *,
                   int *, double *, int *, double *, double *, double *,
                   double *, int *, double *, int *, int *, int *, int *,
                   int *, double *);
extern int REGRID (int *, int *, double *, int *, double *, double *, double *,
                   double *, double *, double *,
                   int *, int *, double *, int *, int *, int *, double *,
                   int *, double *, double *, double *, double *, int *, int *,
                   int *, int *, int *, double *);
extern int DXSPLDER (double *, int *, double *, int *, int *, double *, double *,
                     int *, double *, int *);
extern int DXPARDER (double *, int *, double *, int *, double *, int *, int *,
                     int *, int *, double *, int *, double *, int *, double *,
                     double *, int *, int *, int *, int *);
#endif // RLAB_DIERCKX_H
