/* mathl.h
 * Math Library Functions */

/*  This file is a part of RLaB ("Our"-LaB)
   Copyright (C) 1995  Ian R. Searle
   Copyright (C) 2001  Marijan Kostrun

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

/*
 * Relationships
 */
static inline int lessthen(double x1, double x2, double erel, double eabs)
{
  return (x1 < ((1-erel)*x2 - eabs));
}

static inline int approxequal(double x1, double x2, double erel, double eabs)
{
  return ((((1-erel)*x2 - eabs) <= x1) && (x1 <= ((1+erel)*x2 + eabs)));
}

static inline int greaterthen(double x1, double x2, double erel, double eabs)
{
  return (x1 > ((1+erel)*x2 + eabs));
}

/*
 * Functions for creating Inf and NaNs
 */
extern void   init_inf_nan (void);
extern double create_inf (void);
extern double create_nan (void);
extern int    isnand( double );
extern double rlab_dmin_vector(int, int, double *, int);
extern double rlab_dmax_vector(int, int, double *, int);
extern double rlab_mean2_vector(int, int, double *, int, int *, int);
extern double rlab_var2_vector (int, int, double *, int, int, int *, int);
extern int    rlab_vector_flip (unsigned char *, int, size_t);


#ifndef isinf
static inline int isinf(double x)
{
  if (x == create_inf())
    return 1;

  if (x == -create_inf())
    return 1;

  return 0;
}
#endif

#ifndef isfinite
static inline int isfinite(double x)
{
  if (x == -create_inf())
    return 0;

  if (x == create_inf())
    return 0;

  if (isnand (x))
    return 0;

  return 1;
}
#endif

static inline int isposinf(double x)
{
  if (x == create_inf())
    return 1;

  return 0;
}

static inline int isneginf(double x)
{
  if (x == -create_inf())
    return 1;

  return 0;
}

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


