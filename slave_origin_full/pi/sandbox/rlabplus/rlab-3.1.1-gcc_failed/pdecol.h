// pdecol.h
// (C) 2006, Marijan Kostrun, project rlabplus
//
// Defines, and declarations for a
// C => fortran, RLaB => PDECOL (688.f) interface.
//

#ifndef RLAB_PDECOL_H
#define RLAB_PDECOL_H

#include "fi.h"

#ifdef HAVE_FORTRAN_UND_BACK
  #define PDECOL   pdecol_
  #define PDE1VAL  pde1val_
  #define PDE1F    pde1f_
  #define PDE1FJAC pde1fjac_
  #define PDE1BND  pde1bnd_
  #define PDE1UINT pde1uint_
  #define BACOL    bacol_
  #define BAC1VAL  bac1val_
  #define BAC1BXA  bac1bxa_
  #define BAC1DBXA bac1dbxa_
  #define BAC1BXB  bac1bxb_
  #define BAC1DBXB bac1dbxb_
#endif

#ifdef HAVE_FORTRAN_UND_FRONT
  #define PDECOL   _pdecol
  #define PDE1VAL  _pde1val
  #define PDE1F    _pde1f
  #define PDE1FJAC _pde1fjac
  #define PDE1BND  _pde1bnd
  #define PDE1UINT _pde1uint
  #define BACOL    _bacol
  #define BAC1VAL  _bac1val
  #define BAC1BXA  _bac1bxa
  #define BAC1DBXA _bac1dbxa
  #define BAC1BXB  _bac1bxb
  #define BAC1DBXB _bac1dbxb
#endif

#ifdef HAVE_FORTRAN_UPPERCASE
    // Do nothing, the existing code is OK
#endif

#ifdef HAVE_FORTRAN_LOWERCASE
  #define PDECOL   pdecol
  #define PDE1VAL  pde1val
  #define PDE1F    pde1f
  #define PDE1FJAC pde1fjac
  #define PDE1BND  pde1bnd
  #define PDE1UINT pde1uint
  #define BACOL    bacol
  #define BAC1VAL  bac1val
  #define BAC1BXA  bac1bxa
  #define BAC1DBXA bac1dbxa
  #define BAC1BXB  bac1bxb
  #define BAC1DBXB bac1dbxb
#endif

//
// PDECOL and BACOL FORTRAN77 solvers
//
extern int PDECOL (double *, double *, double *,
                   double *, double *, int *, int *,
                   int *, int *, int *, int *,
                   double *, int *,
                   int *, char *, int *, int *);
extern int BACOL  (double *, double *, double *, double *, int *, int *, int *,
                   int *, double *, int *,
                   double *, int *, int *, int *, double *, int *,
                   int *, char *);

//
// rlab/c functions for PDECOL solver
//
int PDE1F    (double * tp, double * xp, double *up, double *uxp, double *uxxp,
              double * dfdup, int * idummy
             );
int PDE1FJAC (double * tp, double * xp, double *up, double *uxp, double *uxxp,
              double * dfdup, double * dfduxp, double * dfduxxp, int * idummy
             );
int PDE1BND  (double * tp, double *xp, double *up, double *uxp,
             double *dbdup, double * dbdux, double *dzdtp, int * idummy
             );
int PDE1UINT (double * xp, double *up, int * idummy);
int PDE1VAL  (double * xp, double * usolL, double *sctch, int * ndim1,
             int * ndim2, int * npts, double * work);

//
// rlab/c function for BACOL solver
//
int BAC1VAL  (int * kcol, double * xsol, int * nint, double * x, int * npde,
              int * npts, double * usol, double * y, double * work
             );

#endif // RLAB_PDECOL_H
