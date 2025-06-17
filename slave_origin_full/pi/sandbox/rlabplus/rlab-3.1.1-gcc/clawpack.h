// clawpack.h

//
// Defines, and declarations for a
// C => fortran, RLaB => CLAWPACK interface.
//

#ifndef RLAB_CLAW1_H
#define RLAB_CLAW1_H

#include "fi.h"

#ifdef HAVE_FORTRAN_UND_BACK
  #define CLAW1RL   claw1rl_
  #define CLAW1S    claw1s_
  #define BC1LEFT   bc1left_
  #define BC1RIGHT  bc1right_
  #define B4STEP1   b4step1_
  #define SRC1      src1_
  #define CLAW1R    claw1r_
#endif

#ifdef HAVE_FORTRAN_UND_FRONT
  #define CLAW1RL   _claw1rl
  #define CLAW1S    _claw1s
  #define BC1LEFT   _bc1left
  #define BC1RIGHT  _bc1right
  #define B4STEP1   _b4step1
  #define SRC1      _src1
  #define CLAW1R    _claw1r
#endif

#ifdef HAVE_FORTRAN_UPPERCASE
    // Do nothing, the existing code is OK
#endif

#ifdef HAVE_FORTRAN_LOWERCASE
  #define CLAW1RL   claw1rl
  #define CLAW1S    claw1s
  #define BC1LEFT   bc1left
  #define BC1RIGHT  bc1right
  #define B4STEP1   b4step1
  #define SRC1      src1
  #define CLAW1R    claw1r
#endif

//
// FORTRAN77 solvers:
//
extern int CLAW1RL (int *, int *, int *, int *, int *, int *,
                    double *, double *, double *, double *, double *, int *,
                    double *, double *, int *, double *, int *,
                    int *, int *, double *);
extern int CLAW1();
extern int RP1();

//
// rlab/c functions called from FORTRAN
//
int CLAW1S (double *, double *, double *);

int BC1LEFT (int *, int *, int *, int *, double *, double *, double *,
             int *, double *, double *, double *, int *);
int BC1RIGHT (int *, int *, int *, int *, double *, double *, double *,
              int *, double *, double *, double *, int *);
int B4STEP1 (int *, int *, int *, int *, double *, double *, double *dx,
             double *, double *, int *, double *);
int SRC1 (int *, int *, int *, int *, double *, double *, double *,
          int *, double *, double *, double *, int *);
int CLAW1R (double *, double *, double *, double *, double *,
            double *, double *);

#endif // RLAB_CLAW_H
