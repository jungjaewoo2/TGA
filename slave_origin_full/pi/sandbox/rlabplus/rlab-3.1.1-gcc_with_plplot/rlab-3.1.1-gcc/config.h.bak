/* config.h.  Generated automatically by configure.  */
/* RLaB system configuration header file. */

/*
 * This configuration file contains all the "defines" that
 * are necessary to compile rlab. If this file, as generated
 * by configure, does not work, then either "define" or "undef"
 * symbols to generate a working config.h. Please send any
 * changes to searleir@yahoo.com.
 */

#ifndef _CONFIG_H_
#define _CONFIG_H_

/*
 * Is this a UNIX system ?
 * As opposed to an OS/2 or Mac system.
 */

/* #undef unix */
/* #undef winnt */

#ifdef unix
#define HAVE_PIPE
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
#define HAVE_UNISTD_H 1
#define HAVE_TIME_H 1

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
#ifndef size_t
#include <sys/types.h>
#else
typedef unsigned int size_t;
#endif

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
//  R L A B P L U S   -   E X T E N S I O N S: BEGIN
//
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


/*
 * NGSPICE IS HERE if built with --with-ngspice
 *
 * #define HAVE_NGSPICE 1
 * 
 */


/*
 * python is available unless one insists --disable-python
 *
 */
#define HAVE_PYTHON 1


/*
 * jvm is available if built with --with-jvm
 *
 */
/* #undef HAVE_JVM */


/*
 * gphoto2 is available unless one insists --disable-gphoto2
 *
 */
#define HAVE_GPHOTO2 1


/*
 * imagemagick is available unless one insists --disable-im
 *
 */
#define HAVE_IMAGEMAGICK 1

/*
 * ans - make available results of command line calculations
 *
 */
#define HAVE_ANSWER_FOR_ALL "ans"


//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
//  R L A B P L U S   -   E X T E N S I O N S: END
//
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++




/*
 * #define HAVE_METIS if you want to use the
 * Metis Graph Partitioning Library.
 */
/* #undef HAVE_METIS3 */

/*
 * #define HAVE_GC if you have libgc.a
 * available (Boehm's garabage-collector.
 */

#define HAVE_GC 1

/*
 * #define HAVE_SO if your system
 * has shared objects (dynamic linking).
 */

#define HAVE_SO 1
/* #undef HAVE_DL_H */
#define HAVE_DLFCN_H 1
#define HAVE_DLOPEN 1
#define HAVE_DLSYM 1
#define HAVE_DLCLOSE 1

/*
 * #define HAVE_READLINE 1
 * If you want readline command editing.
 * Also edit the Makefile to add the location
 * and name of the readline library.
 */

#define HAVE_READLINE 1

/*
 * #define HAVE_RLAB_PLPLOT
 * If you have the PLPLOT library and
 * you wish to use this for plotting data.
 * Also edit the Makefile to add the location
 * and name of the PLPLOT library, and any other
 * support libraries it may need.
 */
#define HAVE_RLAB_PLPLOT 1

/*
 * #define HAVE_RLAB_PGPLOT
 * If you have the PGPLOT library and
 * you wish to use this for plotting data.
 * Also edit the Makefile to add the location
 * and name of the PLPLOT library, and any other
 * support libraries it may need.
 */
/* #undef HAVE_RLAB_PGPLOT */

/*
 * this funkiness is to allow low power systems to run rlab without
 * SUITE SPARSE (became available since opensuse 12.1)
 *
 * #define HAVE_SUPERLU 1
 * #define HAVE_SUITESPARSE 1
 *  
 *  [azg] commented out
 *
 */

/*
 * this funkiness is to allow low power systems to run rlab without
 * ARPACK (became available since opensuse 12.1) and
 * GLPK
 *
 */
#define HAVE_ARPACK 1
#define HAVE_GLPK 1


/* Does this system have dirent.h or what ? */
#define HAVE_DIRENT_H 1
/* #undef HAVE_SYS_NDIR_H */
/* #undef HAVE_SYS_DIR_H */
/* #undef HAVE_NDIR_H */

/* If your system has an index(3), and rindex(3) */
#define HAVE_RINDEX 1

/* If your math library has rint() */
#define HAVE_RINT 1
#define HAVE_RINT_DEC 1
#ifdef HAVE_RINT
#ifndef HAVE_RINT_DEC
extern double rint();
#endif
#endif

/*
 * If your system has logb.
 * logb is IEEE-754.
 */

#define HAVE_LOGB 1

/* If your system has a difftime() function */
#define HAVE_DIFFTIME 1

/* If your system has times(2) */
#define HAVE_TIMES 1
#define HAVE_SYSCONF 1
#define HAVE_SYS_TYPES_H 1
#define HAVE_SYS_TIMES_H 1

/* If your system has time.h (ANSI-C). */
#define HAVE_TIME_H 1

/* If your system has clock(3) (ANSI-C). */
#define HAVE_CLOCK 1

/* If your system has a putenv() function */
#define HAVE_PUTENV 1

/* If your system has a getenv() function */
#define HAVE_GETENV 1

/* If your system has a sleep() function */
#define HAVE_SLEEP 1

/* This defines the byte significance of your machine */
#define RBIG_ENDIAN 4321
#define RLITTLE_ENDIAN 1234

/* #undef WORDS_BIGENDIAN */
#ifdef WORDS_BIGENDIAN
#define RBYTE_ORDER RBIG_ENDIAN
#define WE_ARE_BIG_ENDIAN 1
#else
#define RBYTE_ORDER RLITTLE_ENDIAN
#define WE_ARE_LITTLE_ENDIAN 1
#endif

/* A couple of key data object sizes */
#define SIZEOF_DOUBLE 8
#define SIZEOF_LONG_INT 8
#define SIZEOF_INT 4
#define SIZEOF_SHORT_INT 2

/*
 * Figure out the maximum value of a signed int.
 */

#define HAVE_VALUES_H 1

#ifndef HAVE_VALUES_H
#if SIZEOF_INT == 2
#define MAXINT  32767
#endif

#if SIZEOF_INT == 4
#define MAXINT  2147483647
#endif
#else
#include <limits.h>
#define MAXINT INT_MAX
#endif  /* HAVE_VALUES_H */

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
#define HAVE_FPU_CONTROL_H 1
/* #undef HAVE___SETFPUCW */

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

#define HAVE_MATHERR 1
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
