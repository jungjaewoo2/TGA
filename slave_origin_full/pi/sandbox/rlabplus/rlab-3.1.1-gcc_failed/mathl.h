/* mathl.h
 * Math Library Functions */

/*  This file is a part of RLaB ("Our"-LaB)
   Copyright (C) 1995  Ian R. Searle

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.

   See the file ./COPYING
   ********************************************************************** */

#ifndef RLAB_MATH_L
#define RLAB_MATH_L

#include "rlab_macros.h"
#include "complex.h"

/*
 * Functions for creating Inf and NaNs
 */

extern void init_inf_nan (void);
extern double create_inf (void);
extern double create_nan (void);

#ifdef THINK_C
#if __option(double_8)
/* 8-byte double */
#define	r__inf_val_bytes  { 0x7f, 0xf0, 0, 0, 0, 0, 0, 0 }
#define	r__nan_bytes      { 0x7f, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff }
#else
/* 12-byte native fp double */
#define	r__inf_val_bytes  { 0x7f,0xff,0,0,0,   0,   0,   0,   0,   0,   0,   0    }
#define	r__nan_bytes      { 0x7f,0xff,0,0,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff }
#endif
#else

#ifdef WORDS_BIGENDIAN
#define	r__inf_val_bytes     { 0x7f, 0xf0, 0, 0, 0, 0, 0, 0 }
#define	r__nan_bytes         { 0x7f, 0xf8, 0, 0, 0, 0, 0, 0 }
#else
#define	r__inf_val_bytes     { 0, 0, 0, 0, 0, 0, 0xf0, 0x7f }
#define	r__nan_bytes         { 0, 0, 0, 0, 0, 0, 0xf8, 0x7f }
#endif
#endif

extern double errno_check (double d, char *s);

#ifndef HAVE_RINT
#define rint Rrint
extern double Rrint _PROTO ((double x));
#endif

extern int detect_inf_r (double *a, int n);
extern int detect_inf_c (Complex * a, int n);
extern int detect_nan_r (double *a, int n);
extern int detect_nan_c (Complex * a, int n);

extern double rfrexp (double value, int *exp);

union r_ieee754_double
{
  double d;

  /* This is the IEEE 754 double-precision format.  */
  struct
  {
#if	RBYTE_ORDER == RBIG_ENDIAN
    unsigned int negative:1;
    unsigned int exponent:11;
    /* Together these comprise the mantissa.  */
    unsigned int mantissa0:20;
    unsigned int mantissa1:32;
#endif				/* Big endian.  */
#if	RBYTE_ORDER == RLITTLE_ENDIAN
    /* Together these comprise the mantissa.  */
    unsigned int mantissa1:32;
    unsigned int mantissa0:20;
    unsigned int exponent:11;
    unsigned int negative:1;
#endif				/* Little endian.  */
  }
  ieee;
};

extern int    isnand( double );
extern double rlab_dmin_vector(int, int, double *);
extern double rlab_dmax_vector(int, int, double *);
extern double rlab_dmax_mdr(MDR *x);
extern double rlab_dmin_mdr(MDR *x);
extern double rlab_mean_vector (int, double *, MDR*, MDR*,int);
extern double rlab_mean2_vector(int, int, double *, int, int *, int);
extern double rlab_var_vector (int, double *, int, MDR*, MDR*);
extern double rlab_var2_vector (int, int, double *, int, int, int *, int);
extern void   rlab_meanstat_from_valstat_vectors(int length, double *x, double *w, MDR *ig_idx, MDR * iuse_idx,
                                               double *m1, double *v1, int do_std);

extern int rlab_vector_flip (unsigned char *, int, size_t);

#endif /* RLAB_MATH_L */
