// ode.h

//
// Defines, and declarations for a
// C => fortran, RLaB => ODE,ACDC,ODEBIM interface.
//

#ifndef RLAB_ODE_H
#define RLAB_ODE_H

#include "fi.h"

#ifdef HAVE_FORTRAN_UND_BACK
  #define ODE     ode_
  #define BIM     bim_
  #define COLDAE  coldae_
  #define APPSLN  appsln_
  #define ACDC    acdc_
  #define SLEIGN  sleign_
  #define SLCOUP  slcoup_
  #define MIRKDC  mirkdc_
  #define DVODE   dvode_
  #define MEBDFI  mebdfi_
  #define DDASKR  ddaskr_
  #define BIMD    bimd_
#endif

#ifdef HAVE_FORTRAN_UND_FRONT
  #define ODE     _ode
  #define BIM     _bim
  #define COLDAE  _coldae
  #define APPSLN  _appsln
  #define ACDC    _acdc
  #define SLEIGN  _sleign
  #define SLCOUP  _slcoup
  #define MIRKDC  _mirkdc
  #define DVODE   _dvode
  #define MEBDFI  _mebdfi
  #define DDASKR  _ddaskr
  #define BIMD    _bimd
#endif

#ifdef HAVE_FORTRAN_UPPERCASE
    // Do nothing, the existing code is OK
#endif

#ifdef HAVE_FORTRAN_LOWERCASE
  #define ODE     ode
  #define BIM     bim
  #define COLDAE  coldae
  #define APPSLN  appsln
  #define ACDC    acdc
  #define SLEIGN  sleign
  #define SLCOUP  slcoup
  #define MIRKDC  mirkdc
  #define DVODE   dvode
  #define MEBDFI  mebdfi
  #define DDASKR  ddaskr
  #define BIMD    bimd
#endif

//
// FORTRAN77 solvers
//
// adams ODEIV solver
extern int ODE (int (*)(double *, double *, double *),
                int *, double *, double *, double *, double *, double *, int *,
                double *, int *);
// ODEBIM solver
extern int BIM (int *,
                int (*)(int *, double *, double *, double *, int *, double * , int *),
                double *, double *, double *, double *, double *, double *,
                int (*)(int *, double *, double *, double *, int *, int *, double *, int *),
                int *, int *, int *, double *, int *, int *, int *,
                double *, int *, int *, int *);

extern int ACDC (int *, int *, int *, double *, double *, int *, double *,
                 int *, int *, double *, int *, int *, int *, int *, double *,
                 int *, double *, int *, int *, double *, int *, int *, int *,
                 double *, double *,
                 int (*)(int *, double *, double *, double *, double *),
                 int (*)(int *, double *, double *, double *, double *),
                 int (*)(int *, int *, double *, double *, double *),
                 int (*)(int *, int *, double *, double *, double *),
                 int *, char *, int *);
extern int COLDAE (int *, int *, int *, double *, double *, double *, int *, int *,
                   double *, double *, int *, double *, int *,
                   int (*)(double *, double *, double *, double *),
                   int (*)(double *, double *, double *, double *),
                   int (*)(int *, double *, double *),
                   int (*)(int *, double *, double *),
                   int (*)(), char *, int *);
extern F_INT APPSLN (); //  solution for COLDAE/COLSYS ODEBV solver
extern int MIRKDC (int *, double *, int *, int *, int *,
                   int *, double *, int *, double *, int *, int *, int *,
                   int *, double *, double *, double *, int *,
                   int (*)(int *, double *, double *, double *),
                   int (*)(int *, double *, double *, double *),
                   int (*)(int *, double *, double *, double *),
                   int (*)(int *, double *, double *, double *),
                   char *, int *);
extern int SLEIGN (double *, double *, int *, double *, double *, double *, double *, double *, double *,
                   double *, double *, int *, double *, double *, int *, int *, double *, int *, int *,
                   char *, int *);
extern int SLCOUP (double *, double *, int *, double *, double *, double *, double *, double *, double *,
                   double *, double *, int *, double *, double *, int *, int *, double *, int *, int *,
                   double *, double *, double *, double *, double *, char *, int *);
extern int DVODE (int (*)(int *, double *, double *, double *, double *, int *),
                  int *, double *, double *, double *, int *, double *, double *, int *,
                  int *, int *, double *, int *, int *, int *,
                  int (*)(int *, double *, double *, int *, int *, double *, int *, double *, int *),
                  int *, double *, int *);
extern F_INT MEBDFI (int *, double *, double *, double *, double *, double *, double *,
                     int *, int *, int *,
                     double *, int *, int *i, int *, int *,
                     int *, double *, double *, double *, int*,
                     int (*)(double * t, double * y, double * pd, int * n, double * yprime,
                             int * mbnd4, double * con, int * ipar, double * rpar, int * ierr),
                     int (*)(int * n, double * t, double * y, double * delta, double * yprime,
                             int * ipar, double * rpar, int * ierr),
                     int *, char *, int *); // MEBDF solver
extern F_INT DDASKR (int (*)(double * t, double * y, double * yprime, double * cj,
                             double * delta, int * ires, double * rpar, int * ipar),
                     int *, double *, double *, double *, double *,
                     int *, double *, double *,
                     int *, double *, int *, int *, int *, double *, int *,
                     int (*)(double * t, double * y, double * yprime, double *pd, double * cj,
                             double * rpar, int * ipar),
                     int(*)(),
                     double *, int *, int *); // DDASKR solver
extern F_INT BIMD   (int *,
                     int (*)(int * m, double * t, double * y, double * dy, int * ierr, double * rpar, int * ipar),
                     double *, double *, double *, double *,
                     double *, double *, int *,
                     int (*)(int * m, double * t, double * y, double * jac, int * ldjac, int * ierr, double * rpar, int * ipar),
                     int *, int *, int *,
                     double *, int *, int *, int *, int (*)(), int *,
                     double *, int *, int *, int *,
                     double *, int *, int *); // BIMD solver
#endif // RLAB_ODE_H
