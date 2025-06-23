/*
 * mathl.c
 * Math Library Functions.
 * Some pieces of code borrowed from the GNU C-library
 * (some of which, in turn, were borrowed from UCB)
 */

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

#include "rlab.h"
#include "mdr.h"
#include "mathl.h"
#include "util.h"
#include "rlab_macros.h"

#include <stdio.h>
#include <math.h>
#include <errno.h>

//
//
//  system operators
//
//
// bracketing of plus and minus operators
double rlab_op_zero_abs=0.0;      // |res| <= zabs => res = 0
double rlab_op_zero_rel=0.0;      // |res| <= zrel*max(arg1,arg2) => res=0
double rlab_op_min;  // arg1 - arg2 = max(arg1-arg2, sub_min)
double rlab_op_max;  // arg1 + arg2 = min(arg1+arg2, add_max)


/* **************************************************************
 * Check errno when calling math library functions.
 * ************************************************************** */

double errno_check (double d, char *s)
{
  if (errno == EDOM)
  {
    errno = 0;
    fprintf (stderr, "%s ", s);
    rerror ("argument out of domain");
  }
  else if (errno == ERANGE)
  {
    errno = 0;
    fprintf (stderr, "%s ", s);
    rerror ("result out of range");
  }
  return (d);
}

#ifdef USE_MATHERR

/* **************************************************************
 * RLaB's matherr(). For now, use error_1 when matherr() is called.
 * This will suffice until I figure out what to do with the
 * various "standard", but incompatible math libraries.
 * ************************************************************** */

#ifdef win32
/* Contrary to their own documentation, Microsoft defines _exception, not exception. */
#define exception _exception
#endif

int matherr (struct exception *x)
{
  switch (x->type)
  {
  case DOMAIN:
    rerror ("argument out of DOMAIN");
    break;
  case SING:
    rerror ("argument SINGULARITY");
    break;
  case OVERFLOW:
    rerror ("OVERFLOW");
    break;
  case UNDERFLOW:
    return (int) (x->retval);
    break;
  case TLOSS:
    rerror ("Total LOSS of precision");
    break;
  case PLOSS:
    rerror ("Partial LOSS of precision");
    break;
  }
  return (1);
}
#endif /* USE_MATHERR */

#ifndef HAVE_RINT
/* **************************************************************
 * A replacement rint() for deficient systems. This is not an IEEE
 * compatible rint(). This is a SIMPLE-MINDED function (it is better
 * than nothing :-).
 * ************************************************************** */

double
Rrint (double x)
{
  if (x != x)			/* NaN */
    return (x);
  return (x >= 0.0 ? floor (x + 0.5) : -floor (0.5 - x));
}

#endif /* ! HAVE_RINT */

/* **************************************************************
 * A function to check a double array for the presence of an
 * Inf or a NaN. This is to "protect" the Fortran subroutines
 * from getting "hung" on this type of input.
 * ************************************************************** */
/*
 * Initialize global variables for floating point operations.
 * It is crucial that the processor generate Infs and NaNs at
 * this point in program execution. Floating point exception
 * handling can be altered at a later time.
 */

/* Our own (private) definitions of Inf and NaN. */
static double R_INF;
static double R_NAN;

void init_inf_nan (void)
{
  double num = 1.0;
  double denom = 0.0;

  /* Generate an Inf. */
  R_INF = num / denom;

#if !defined(win32)
  /* Generate a NaN. */
  num = 0.0;
  R_NAN = num / denom;
#else
  {
    /* Generate a NaN. */
    const unsigned char __nan[8] = r__nan_bytes;
    R_NAN = (*(const double *) __nan);
  }
#endif /* win32 */

  rlab_op_min = -R_INF;
  rlab_op_max =  R_INF;
}

double create_inf (void)
{
  return R_INF;
}

double create_nan (void)
{
  return R_NAN;
}

int isnand (double x)
{
  // isnan has become a macro: this is the outcome if the argument is
  // double.
  return __isnan(x);
}

/*
 * Create Infs, and NaNs
 */
int detect_inf_r (double *a, int n)
{
  int i;

  for (i = 0; i < n; i++)
  {
    if ((R_INF==a[i]) || (R_INF==-a[i]))
      return (1);
  }
  return (0);
}

int detect_inf_c (Complex * a, int n)
{
  int j;

  for (j = 0; j < n; j++)
  {
    if ( (R_INF == RE(a[j])) || (R_INF == -RE(a[j]))
          || (R_INF == IM(a[j])) || (R_INF == -IM(a[j])) )
      return (1);
  }
  return (0);
}

int detect_nan_r (double *a, int n)
{
  int i;

  for (i = 0; i < n; i++)
  {
    if (isnand(a[i]))
      return (1);
  }
  return (0);
}

int detect_nan_c (Complex * a, int n)
{
  int j;

  for (j = 0; j < n; j++)
  {
    if (isnand(RE(a[j])) || isnand(RE(a[j])))
      return (1);
  }
  return (0);
}

/* **************************************************************
 * Our own frexp(). Borrowed from the GNU libc.
 * ************************************************************** */

double rfrexp (double value, int *exp)
{
  if (value == 0)
  {
    *exp = 0;
    return 0;
  }
  else
  {
    union r_ieee754_double u;
    u.d = value;
    *exp = u.ieee.exponent - 1022;
    u.ieee.exponent = 1022;
    return u.d;
  }
}


int rlab_vector_flip(unsigned char *x, int n, size_t size)
{
  int i,k;
  unsigned char dummy;

  for (i=0; i < (n>>1); i++)
  {
    for (k=0; k<size; k++)
    {
      dummy = x[i*size + k];
      x[i*size + k]  = x[(n-i-1)*size + k];
      x[(n-i-1)*size + k] = dummy;
    }
  }

  return 0;
}

//
// determine a minimum of a matrix as quickly as possible
//
double rlab_dmin_vector(int length, int istrid, double *x, int ignore_inf)
{
  double rval=create_nan();
  int i,j;
  if (length>0)
  {
    // find first not-nan
    j=0;
    while (j<length)
    {
      rval=x[j*istrid];
      j++;

      if (ignore_inf)
      {
        if (isfinite(rval))
          break;
      }
      else
      {
        if(!isnand (rval))
          break;
      }
    }
    for (i=j; i<length; i++)
    {
      if(isnand (rval))
        continue;
      if (ignore_inf)
      {
        if (isinf(rval))
          continue;
      }
      if(rval >= x[i*istrid])
        rval = x[i*istrid];
    }
  }

  return rval;
}

// find max of a matrix, ignore nan's
double rlab_dmax_vector(int length, int istrid, double *x, int ignore_inf)
{
  double rval=create_nan();
  int i,j;
  if (length >0)
  {
    // find first not-nan
    j=0;
    while (j<length)
    {
      rval=x[j*istrid];
      j++;

      if (ignore_inf)
      {
        if (isfinite(rval))
          break;
      }
      else
      {
        if(!isnand (rval))
          break;
      }
    }
    // now do the comparison
    for (i=j; i<length; i++)
    {
      if(isnand (rval))
        continue;
      if (ignore_inf)
      {
        if (isinf(rval))
          continue;
      }
      if(rval <= x[i*istrid])
        rval = x[i*istrid];
    }
  }

  return rval;
}

double rlab_var2_vector(int length, int istrid, double *x, int bias, int lig_idx, int *ig_idx, int isw)
{
  static double s1, s2;
  static double count=0;
  int i, count2;
  double rval, s1_2, s2_2, dummy;

  // if isw=0 then calculate
  //  s1    -> sum of all x
  //  count -> how many x's there are in total including excluded indices
  // if isw>0 then use s1 and count from the earlier calculations
  if (!isw || isnand(s1)|| isnand(s1) || count==0)
  {
    if (length)
    {
      s1 = 0;
      s2 = 0;
      count = 0;
      for (i=0; i<length; i++)
      {
        dummy = x[i*istrid];

        if(isnand (dummy))
          continue;

        s1 += dummy;
        s2 += dummy * dummy;
        count ++;
      }
    }
  }

  if(!length || !count)
  {
    s1 = create_nan();
    s2 = create_nan();
    return s1;
  }

  s1_2=0;
  s2_2=0;
  count2=0;
  if (lig_idx)
  {
    MDR * r_idx = mdi_Create(1,length);
    mdr_Ones(r_idx); // book keeping array so that any index can be removed only once

    // sum all x's with the indices to be removed
    for (i=0; i<lig_idx; i++)
    {
      // tard proofing: ignore indices that go beyond the range of x
      if ((ig_idx[i] < 0) || (ig_idx[i] > length - 1))
        continue;

      dummy = x[ig_idx[i]*istrid];

      if( isnand( dummy ) )
        continue;

      // to prevent tards from having ig_idx=[ ..j .. j ] to remove same element twice
      if ( MdiV0(r_idx, ig_idx[i]) )
      {
        MdiV0(r_idx,ig_idx[i]) = 0;
        s1_2 += dummy;
        s2_2 += dummy*dummy;
        count2 ++;
      }
    }

    mdr_Destroy(r_idx);
  }

  if (count > count2 + bias)
    rval = (s2 - s2_2) / (count - count2 - bias) - (s1 - s1_2)*(s1 - s1_2) / (count - count2 -bias) / (count - count2);
  else
    rval = create_nan(); // though 0/0 should produce nan as well

  return rval;
}

void recast_double_as_int_in_double_storage(int n, void * data)
{
  double *x = (double *) data;
  int    *i = (int *) data;
  int j;
  for (j=0; j<n; j++)
  {
    i[j] = (int) x[j];
  }
  return;
}

void recast_int_as_double_in_double_storage(int n, void * data)
{
  double *x = (double *) data;
  int    *i = (int *) data;
  int j;
  for (j=n-1; j>=0; j--)
  {
    x[j] = i[j];
  }
  return;
}









