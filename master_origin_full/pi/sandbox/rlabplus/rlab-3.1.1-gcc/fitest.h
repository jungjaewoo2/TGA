/* fitest.h */

/*
 * Defines, and declarations for a
 * C => fortran, RLaB => BLAS interface.
 */

#ifndef RLAB_FITEST_H
#define RLAB_FITEST_H

#include "fi.h"

#ifdef HAVE_FORTRAN_UND_BACK

#define TEST test_

#endif

#ifdef HAVE_FORTRAN_UND_FRONT

#define TEST _test

#endif

#ifdef HAVE_FORTRAN_UPPERCASE
/* Do nothing, the existing code is OK */
#endif

#ifdef HAVE_FORTRAN_LOWERCASE

#define TEST test

#endif

extern F_INT TEST ();

#endif  /* RLAB_FITEST_H */
