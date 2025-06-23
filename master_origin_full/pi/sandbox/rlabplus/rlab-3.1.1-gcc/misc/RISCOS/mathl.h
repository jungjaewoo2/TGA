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

#include "complex.h"

/*
 * Define Inf (Infinity) and NaN (Not a Number)
 */
#ifdef __riscos
#define	r__inf_val_bytes     { 0, 0, 0xf0, 0x7f, 0, 0, 0, 0 }
#define	r__nan_bytes         { 0, 0, 0xf8, 0x7f, 0, 0, 0, 0 }
#else
/*#ifdef THINK_C
#if __option(double_8)*/
/* 8-byte double */
/*#define	r__inf_val_bytes  { 0x7f, 0xf0, 0, 0, 0, 0, 0, 0 }
#define	r__nan_bytes      { 0x7f, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff }
#else*/
/* 12-byte native fp double */
/*#define	r__inf_val_bytes  { 0x7f,0xff,0,0,0,   0,   0,   0,   0,   0,   0,   0    }
#define	r__nan_bytes      { 0x7f,0xff,0,0,0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff }
#endif
#else*/

#ifdef WORDS_BIGENDIAN
#define	r__inf_val_bytes     { 0x7f, 0xf0, 0, 0, 0, 0, 0, 0 }
#define	r__nan_bytes         { 0x7f, 0xf8, 0, 0, 0, 0, 0, 0 }
#else
#define	r__inf_val_bytes     { 0, 0, 0, 0, 0, 0, 0xf0, 0x7f }
#define	r__nan_bytes         { 0, 0, 0, 0, 0, 0, 0xf8, 0x7f }
#endif
/*#endif*/
#endif /*__riscos*/


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

#endif /* RLAB_MATH_L */
