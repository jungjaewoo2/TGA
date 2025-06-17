/* fi.h */

/*
 * Defines for C => Fortran interface.
 */

#ifndef FORTRAN_C_INT_H
#define FORTRAN_C_INT_H

/*
 * Various defines for different
 * Fortran <=> C interfaces.
 */

#ifdef USE_F2C

#define HAVE_FORTRAN_UND_BACK 1

#else

#ifdef USE_UPPER
#define HAVE_FORTRAN_UPPERCASE 1
#endif

#ifdef USE_LOWER
#define HAVE_FORTRAN_LOWERCASE 1
#endif

#ifdef USE_BACK_UND
#define HAVE_FORTRAN_UND_BACK 1
#endif

#ifdef USE_FRONT_UND
#define HAVE_FORTRAN_UND_FRONT 1
#endif

#endif  /* USE_F2C */


/*
 * The following must be set according to your 
 * platforms Fortran convention. Some Fortran
 * compilers have long int INTEGERS, some have 
 * int INTEGERS. If float is not equivalent to 
 * a REAL, and double is not equivalent to
 * DOUBLE PRECISION on your system, then you are in
 * trouble.
 */

typedef int       F_INT;
typedef float     F_REAL;
typedef double    F_DOUBLE;

#endif /* FORTRAN_C_INT_H */
