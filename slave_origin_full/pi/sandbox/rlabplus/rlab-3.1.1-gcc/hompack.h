// hompack.h: header file for fortran libraries
//  libhompack.a, libcontin.a
// (C) 2006, Marijan Kostrun, project rlabplus
//
// Defines, and declarations for a
// C => fortran, RLaB => HOMPACK, CONTIN interface.
//

#ifndef RLAB_HOMPACK_H
#define RLAB_HOMPACK_H

#include "fi.h"

#ifdef HAVE_FORTRAN_UND_BACK
  // hompack
  #define FIXPDF   fixpdf_
  #define FIXPNF   fixpnf_
  #define FIXPQF   fixpqf_
  #define HOM1F    hom1f_
  #define HOM1FJAC hom1fjac_
  #define HOM1RHO  hom1rho_
  #define HOM1RHOJ hom1rhoj_
  // contin
  #define PITCON   pitcon_
  #define DENSLV   denslv_
#endif

#ifdef HAVE_FORTRAN_UND_FRONT
  // hompack
  #define FIXPDF   _fixpdf
  #define FIXPNF   _fixpnf
  #define FIXPQF   _fixpqf
  #define HOM1F    _hom1f
  #define HOM1FJAC _hom1fjac
  #define HOM1RHO  _hom1rho
  #define HOM1RHOJ _hom1rhoj
  // contin
  #define PITCON   _pitcon
  #define DENSLV   _denslv
#endif

#ifdef HAVE_FORTRAN_UPPERCASE
    // Do nothing, the existing code is OK
#endif

#ifdef HAVE_FORTRAN_LOWERCASE
  // hompack
  #define FIXPDF   fixpdf
  #define FIXPNF   fixpnf
  #define FIXPQF   fixpqf
  #define HOM1F    hom1f
  #define HOM1FJAC hom1fjac
  #define HOM1RHO  hom1rho
  #define HOM1RHOJ hom1rhoj
  // contin
  #define PITCON   pitcon
  #define DENSLV   denslv
#endif

//
// HOMPACK FORTRAN77 dense solvers: fixpdf, fixpnf, fixpqf
//
extern int FIXPDF (int *, double *, int *, double *, double *,
                   int *, double *, int *, int *,
                   double *, double *, double *, double *,
                   double *, double *, int *,
                   double *, double *, double *,
                   double *, int *,
                   char *, int *, int *, double *, double *);
extern int FIXPNF (int *, double *, int *, double *, double *, double *, double *,
                   int *, double *, int *, double *,
                   double *, double *, double *, double *,
                   double *, double *, int *,
                   double *, double *, double *, double *, double *, double *, int *,
                   char *, int *, int *, double *);
extern int FIXPQF (int *, double *, int *, double *, double *, double *, double *,
                   int *, double *, int *, double *,
                   double *, double *, double *, double *,
                   double *, double *, double *, double *, double *,
                   double *, double *, double *, double *, double *, int *,
                   char *, int *, int *);
extern int PITCON (int (*)(int *, double *, int *, double *, double *, int *),
                   double *,
                   int (*)(int *, double *, int *, double *, double *, int *),
                   int *, int *, int *, int *, int *, double *,
                   int *, double *,
                   int (*)(),
                   char *, int *);
extern int DENSLV ();

//
// rlab/c functions for HOMPACK solvers read README
//
int HOM1F    (double *, double *);
int HOM1FJAC (double *, double *, int *);
int HOM1RHO  (double *, double *, double *, double *, double *, int *);
int HOM1RHOJ (double *, double *, double *, double *, int *, double *, int *);
//
// rlab/c functions for CONTIN solver read README
//
int PIT1RHO (int *, double *, int *, double *, double *, int *);
int PIT1RHOJ(int *, double *, int *, double *, double *, int *);

#endif // RLAB_HOMPACK_H
