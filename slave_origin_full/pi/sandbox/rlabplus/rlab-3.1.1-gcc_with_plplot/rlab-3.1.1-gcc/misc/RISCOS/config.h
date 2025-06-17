/* config.h.  Generated automatically by configure.  */
/* RLaB system configuration header file. */

/*
 * This configuration file contains all the "defines" that
 * are necessary to compile rlab. If this file, as generated
 * by configure, does not work, then either "define" or "undef"
 * symbols to generate a working config.h. Please send any
 * changes to ians@eskimo.com.
 */

#ifndef _CONFIG_H_
#define _CONFIG_H_

/*
 * Is this a UNIX system ?
 * As opposed to an OS/2 or Mac system.
 */

/* #undef unix */

#ifdef unix
#define HAVE_PIPE
#endif

#ifdef __riscos
#define HAVE_FORTRAN_UND_BACK
#define __GNU_LIBRARY__ /* only used by getopt */
#endif

/*
 * Does the C compiler support const ?
 * If it doesn't then:
 * #define const
 */

/* #undef const */

/* Standard  C (ANSI) header files */
#define STDC_HEADERS 1

/* Does this system have stdlib.h */
#define HAVE_STDLIB_H 1
#ifdef HAVE_STDLIB_H
#include <stdlib.h>
#else
#include <malloc.h>
#endif

/* Does this system have unistd.h (for chdir()) */
#undef HAVE_UNISTD_H
#define HAVE_TIME_H 1
#define HAVE_CLOCK 1

/* Set ANSI compiler features */
#ifdef __STDC__
#define _PROTO(proto)  proto
typedef void *VPTR;
/* #undef YY_USE_PROTOS */
#else
#define _PROTO(proto)  ()
/* #undef YY_USE_PROTOS */
typedef char *VPTR;
#endif  /* __STDC__ */

/*
 * If this system does not have a size_t type,
 * then use unsigned int. Else, include the file
 * that defines size_t.
 */

/* #undef size_t */
#include <stddef.h>

/*#ifndef size_t
#include <sys/types.h>
#else
typedef unsigned int size_t;
#endif
*/

/*
 * #define HAVE_UMFPACK if you want to use the
 * UMFPACK sparse solver (see the license in the
 * umfpack directory).
 */

#undef HAVE_UMFPACK

/*
 * #define HAVE_GC if you have libgc.a
 * available (Boehm's garabage-collector.
 */

#undef HAVE_GC

/*
 * #define HAVE_SO if your system
 * has shared objects (dynamic linking).
 */

#undef HAVE_SO
#undef HAVE_DLFCN_H
#undef HAVE_DLOPEN
#undef HAVE_DLSYM
#undef HAVE_DLCLOSE

/*
 * #define HAVE_READLINE 1
 * If you want readline command editing.
 * Also edit the Makefile to add the location
 * and name of the readline library.
 */

#undef HAVE_READLINE

/*
 * #define HAVE_RLAB_PLPLOT
 * If you have the PLPLOT library and
 * you wish to use this for plotting data.
 * Also edit the Makefile to add the location
 * and name of the PLPLOT library, and any other
 * support libraries it may need.
 */

#define HAVE_RLAB_PLPLOT 1

/* Does this system have dirent.h or what ? */
#ifdef __riscos
#define DIRENT
#endif

#undef HAVE_DIRENT_H
/* #undef HAVE_SYS_NDIR_H */
/* #undef HAVE_SYS_DIR_H */
/* #undef HAVE_NDIR_H */

/* If your system has an index(3), and rindex(3) */
#undef HAVE_RINDEX

/* If your math library has rint() */
#undef HAVE_RINT
#undef HAVE_RINT_DEC
#ifdef HAVE_RINT
#ifndef HAVE_RINT_DEC
extern double rint();
#endif
#endif

/*
 * If your system has logb.
 * logb is IEEE-754.
 */

/* #undef HAVE_LOGB */

/* If your system has a difftime() function */
#define HAVE_DIFFTIME 1

/* If your system has times(2) */
#undef HAVE_TIMES
#undef HAVE_SYSCONF
#undef HAVE_SYS_TYPES_H
#undef HAVE_SYS_TIMES_H

/* If your system has a putenv() function */
#undef HAVE_PUTENV

/* This defines the byte significance of your machine */
#define RBIG_ENDIAN 4321
#define RLITTLE_ENDIAN 1234

#undef WORDS_BIGENDIAN
#ifdef WORDS_BIGENDIAN
#define RBYTE_ORDER RBIG_ENDIAN
#else
#define RBYTE_ORDER RLITTLE_ENDIAN
#endif

/* A couple of key data object sizes */
#define SIZEOF_DOUBLE 8
#define SIZEOF_LONG_INT 4
#define SIZEOF_INT 4
#define SIZEOF_SHORT_INT 2

/*
 * If your system does not have a float.h, then try building,
 * and running misc/enquire.c. Enquire will generate a suitable
 * float.h. If you don't have, or can't build a float.h, then you
 * will need to define DBL_EPSILON.
 */

#define HAVE_FLOAT_H 1
#ifndef HAVE_FLOAT_H
#define DBL_EPSILON 2.22e-16
#else
#include <float.h>
#endif

/*
 * The width of your terminal. This will be replaced later on.
 */

#define TERM_WIDTH 72

/*
 * Here is where we determine what functions and headers
 * the system offers for dealing with floating point
 * exceptions.
 */

/*
 * We can turn ALL floating-point setup ON or OFF.
 * ON if SETUP_FPE is defined.
 * OFF if SETUP_FPE is not defined.
 */

#define SETUP_FPE 1

/* SVR3.2, SVR4.x method */
/* #undef HAVE_IEEEFP_H */
/* #undef HAVE_FPSETMASK */

/* Berkeley/SunOS-4 method */
/* #undef HAVE_FLOATINGPOINT_H */
/* #undef HAVE_IEEE_HANDLER */

/* DEC-Alpha method */
/* #undef HAVE_MACHINE_FPU_H */
/* #undef HAVE_IEEE_SET_FP_CONTROL */

/* Linux method. */
#undef HAVE_FPU_CONTROL_H
#undef HAVE___SETFPUCW

/*
 * The following tells RLaB how to handle math library errors.
 * `#undef USE_MATHERR'
 * `#define errcheck( a , b ) errno_check( a , b )'
 * will tell RLaB to check errno after every call to math library
 * functions. This method should work on any UNIX platform.
 *
 * `#define USE_MATHERR'
 * `#define errcheck( a , b )    a'
 * will tell RLaB to use matherr(). This is more efficient, since
 * matherr() is only called when an exception occurs. However
 * matherr() may not be available on BSD-style systems.
 */

#undef HAVE_MATHERR
#ifdef HAVE_MATHERR
#define USE_MATHERR 1
#define errcheck( a , b )    a
#else
#define errcheck( a , b )    errno_check( a , b )
#endif

/*
 * Misc defines for systems without declarations
 * in their header files (SUN).
 */

#define HAVE_FPRINTF_DEC 1
#define HAVE_SPRINTF_DEC 1
#define HAVE_FREAD_DEC 1

#endif /* _CONFIG_H_ */
