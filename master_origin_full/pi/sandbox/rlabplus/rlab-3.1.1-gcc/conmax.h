// conmax.h: header file for fortran files
//  conmax.f, 811.f
// (C) 2006, Marijan Kostrun, project rlabplus
//
// Defines, and declarations for a
// C => fortran, RLaB => conmax,811 interface.
//

#ifndef RLAB_CONMAX_H
#define RLAB_CONMAX_H

#include "fi.h"

#ifdef HAVE_FORTRAN_UND_BACK
  #define CONMAX2  conmax2_
  #define CMXFNSET cmxfnset_
  #define PBUNU    pbunu_
  #define PBUNL    pbunl_
  #define PNEWL    pnewl_
  #define PNEWU    pnewu_
  #define PMINL    pminl_
  #define PMINU    pminu_
  #define PMFUNDER pmfunder_
  #define PMFUN    pmfun_
  #define PMDER    pmder_
  #define PMHES    pmhes_
#endif

#ifdef HAVE_FORTRAN_UND_FRONT
  #define CONMAX2   _conmax2
  #define CMXFNSET _cmxfnset
  #define PMINL    _pminl
  #define PBUNU    _pbunu
  #define PNEWU    _pnewu
  #define PBUNL    _pbunl
  #define PNEWL    _pnewl
  #define PMINU    _pminu
  #define PMFUNDER _pmfunder
  #define PMFUN    _pmfun
  #define PMDER    _pmder
  #define PMHES    _pmhes
#endif

#ifdef HAVE_FORTRAN_UPPERCASE
    // Do nothing, the existing code is OK
#endif

#ifdef HAVE_FORTRAN_LOWERCASE
  #define CONMAX2  conmax2
  #define CMXFNSET cmxfnset
  #define PMINL    pminl
  #define PMINU    pminu
  #define PBUNU    pbunu
  #define PNEWU    pnewu
  #define PBUNL    pbunl
  #define PNEWL    pnewl
  #define PMFUNDER pmfunder
  #define PMFUN    pmfun
  #define PMDER    pmder
  #define PMHES    pmhes
#endif

//
// conmax FORTRAN77 solver
//
extern int CONMAX2 (int *, int *, int *, int *, double *, int *,
                    double *, int *, int *,
                    int *, int *, double *, int *,
                    int *, double *, double *,
                    char *, int *,
                    int (*)(int *, int *, double *, int *, int *, double *, int *, int *, int *, double *)
                   );

//
// proxImal bundle FORTRAN77 solvers
//
extern int PBUNU  (int *, int *, double *, int *, double *,
                   int *, double *, double *, double *, int *,
                   char *, int *);
extern int PNEWU  (int *, int *, double *, int *, double *,
                   int *, double *, double *, double *, int *, int *,
                   char *, int *);
extern int PBUNL  (int *, int *, int *, int *, double *, int *, double *, double *,
                   double *, int *, double *, double *, double *,
                   int *, double *, int *, double *, double *, double *, int *,
                   char *, int *);
extern int PNEWL  (int *, int *, int *, int *, double *, int *, double *, double *,
                   double *, int *, double *, double *, double *,
                   int *, double *, int *, double *, double *, double *,
                   int *, int *,
                   char *, int *);
extern int PMINL  (int *, int *, int *, int *, double *, int *, double *, double *,
                   double *, int *, double *, double *, double *, double *,
                   int *, double *, int *, double *, double *, double *, int *, int *,
                   char *, int *);
extern int PMINU  (int *, int *, double *, double *,
                   int *, double *, int *, double *, double *, double *, int *, int *,
                   char *, int *);

//
// rlab/c functions for 811 (proximal bundle) solver
//
int PMFUNDER (int *, double *, double *, double *);
int PMFUN    (int *, int *, double *, double *);
int PMDER    (int *, int *, double *, double *);
int PMHES    (int *, double *, double *);
#endif // RLAB_HOMPACK_H
