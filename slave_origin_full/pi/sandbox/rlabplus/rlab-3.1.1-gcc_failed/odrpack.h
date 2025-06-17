// odrpack.h

//
// Defines, and declarations for a
// C => fortran, RLaB => ODRPACK interface.
//

#ifndef RLAB_ODRPACK_H
#define RLAB_ODRPACK_H

#include "fi.h"

#ifdef HAVE_FORTRAN_UND_BACK
  #define DODRCRL dodrcrl_
  #define DODRL   dodrl_
  #define ODRFDF  odrfdf_
#endif

#ifdef HAVE_FORTRAN_UND_FRONT
  #define DODRCRL _dodrcrl
  #define DODRL   _dodrl
  #define ODRFDF  _odrfdf
#endif

#ifdef HAVE_FORTRAN_UPPERCASE
    // Do nothing, the existing code is OK
#endif

#ifdef HAVE_FORTRAN_LOWERCASE
  #define DODRCRL dodrcrl
  #define DODRL   dodrl
  #define ODRFDF  odrfdf
#endif

//
// ODRPACK solvers
//
extern int DODRCRL (int (*)(int * N, int * M, int * NP, int * NQ, int * LDN, int * LDM,
                            int * LDNP, double * BETA, double * XPLUSD, int * IFIXB,
                            int * IFIXX, int * LDIFX, int * IDEVAL, double * F,
                            double * FJACB, double * FJACD, int * ISTOP),
                    int *, int *, int *, int *,
                    double *, double *, int *, double *, int *,
                    double *, int *, int *, double *, int *, int *,
                    int *, int *, int *, int *, int *, double *,
                    double *, double *, int *, char *, int *,
                    double *, double *, int *,
                    double *, double *, int *, double *, int *, int *, int *, int *);
extern int ODRFDF(int * N, int * M, int * NP, int * NQ, int * LDN, int * LDM,
                  int * LDNP, double * BETA, double * XPLUSD, int * IFIXB,
                  int * IFIXX, int * LDIFX, int * IDEVAL, double * F,
                  double * FJACB, double * FJACD, int * ISTOP);
extern int DODRL   ();


#endif // RLAB_ODRPACK_H
