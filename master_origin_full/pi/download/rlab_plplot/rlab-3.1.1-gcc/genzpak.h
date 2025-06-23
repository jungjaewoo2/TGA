// genzpak.h
// (C) 2006, Marijan Kostrun, project rlabplus
//
// Defines, and declarations for a
// C => fortran, RLaB => genzpak by Alan Genz
//

#ifndef RLAB_GENZPAK_H
#define RLAB_GENZPAK_H

#include "fi.h"

#ifdef HAVE_FORTRAN_UND_BACK
# define PATSYM   patsym_
# define PATFUN   patfun_
# define SMPINT   smpint_
# define DCUHRE   dcuhre_
# define DECUHR   decuhr_
# define INVLTF   invltf_
#endif

#ifdef HAVE_FORTRAN_UND_FRONT
# define PATSYM   _patsym
# define PATFUN   _patfun
# define SMPINT   _smpint
# define DCUHRE   _dcuhre
# define DECUHR   _decuhr
# define INVLTF   _invltf
#endif

#ifdef HAVE_FORTRAN_UPPERCASE
  // Do nothing, the existing code is OK
#endif

#ifdef HAVE_FORTRAN_LOWERCASE
# define PATSYM   patsym
# define PATFUN   patfun
# define SMPINT   smpint
# define DCUHRE   dcuhre
# define DECUHR   decuhr
# define INVLTF   invltf
#endif

//
// genzpak solvers
//
extern int PATSYM (int *, double *, double *, int *, int *, int *,
                     int (*)(int *ndim, double *x, int *numfun, double *funvls),
                     double *, double *, int *, double *,
                     double *, int *, int *, double *);
extern int SMPINT (int *, int *f, int *, int *,
                     int (*)(int *ndim, double *x, int *numfun, double *funvls),
                     double *, double *, int *, int *, int *, double *,
                     int *, double *, double *, int *, int *);
extern int DCUHRE (int *, int *, double *, double *, int *, int *,
                     int (*)(int *ndim, double *x, int *numfun, double *funvls),
                     double *, double *, int *, int *, int *, double *,
                     double *, int *, int *, double *);
extern int DECUHR (int *, int *, double *, double *, int *, int *,
                     int (*)(int *ndim, double *x, int *numfun, double *funvls),
                     int *, double *, int *, double *, double *, int *,
                     int *, int *, int *, int *, double *, double *,
                     int *, int *, double *, int *);

// invltf solver
extern int INVLTF (double *, double *,
                   Complex (*)(Complex *),
                   double *, double *, int *,
                   double *, double *, int *, Complex *, int *);

//
// rlab/c functions for patpack.f and simpack.f solvers
//
extern int PATFUN (int *ndim, double *x, int *numfun, double *funvls);

#endif // RLAB_GENZPAK_H
