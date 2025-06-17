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

double
errno_check (double d, char *s)
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

int
matherr (struct exception *x)
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
int
detect_inf_r (double *a, int n)
{
  int i;

  for (i = 0; i < n; i++)
  {
    if ((R_INF==a[i]) || (R_INF==-a[i]))
      return (1);
  }
  return (0);
}

int
detect_inf_c (Complex * a, int n)
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

int
detect_nan_r (double *a, int n)
{
  int i;

  for (i = 0; i < n; i++)
  {
    if (isnand(a[i]))
      return (1);
  }
  return (0);
}

int
detect_nan_c (Complex * a, int n)
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

double
rfrexp (double value, int *exp)
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
double rlab_dmin_vector(int length, int istrid, double *x)
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
      if(!isnand (rval))
        break;
    }
    for (i=j; i<length; i++)
    {
      if(rval >= x[i*istrid])
        rval = x[i*istrid];
    }
  }

  return rval;
}

// find max of a matrix, ignore nan's
double rlab_dmax_vector(int length, int istrid, double *x)
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
      if(!isnand (rval))
        break;
    }
    // now do the comparison
    for (i=j; i<length; i++)
    {
      if(rval <= x[i*istrid])
        rval = x[i*istrid];
    }
  }

  return rval;
}


// find max of a matrix, ignore nan's
double rlab_dmax_mdr(MDR *x)
{
  double rval=create_nan();
  int length = MNR(x) * MNC(x);

  int i,j;

  if (length >0)
  {
    // find first not-nan
    j=0;
    while (j<length)
    {
      rval=mdrV0(x,j);
      j++;
      if(!isnand (rval))
        break;
    }
    // now do the comparison
    for (i=j; i<length; i++)
    {
      if(isnand (rval))
        continue;
      if(rval <= mdrV0(x,i))
        rval = mdrV0(x,i);
    }
  }

  return rval;
}

// find min of a matrix, ignore nan's
double rlab_dmin_mdr(MDR *x)
{
  double rval=create_nan();
  int length = MNR(x) * MNC(x);

  int i,j;

  if (length >0)
  {
    // find first not-nan
    j=0;
    while (j<length)
    {
      rval=mdrV0(x,j);
      j++;
      if(!isnand (rval))
        break;
    }
    // now do the comparison
    for (i=j; i<length; i++)
    {
      if(isnand (rval))
        continue;
      if(rval >= mdrV0(x,i))
        rval = mdrV0(x,i);
    }
  }

  return rval;
}



// find a mean of vector - ignore nan's, inf's
double rlab_mean_vector(int length, double *x, MDR *ig_idx, MDR * iuse_idx, int ig_infs)
{
  double rval = create_nan();
  int i, k;
  double count=0;
  MDR * r_idx=0;

  if (length==0)
    return rval;

  if (!ig_idx && !iuse_idx)
  {
    for (i=0; i<length; i++)
    {
      if(isnand (x[i]))
        continue;

      count ++;

      if (isnand(rval))
        rval = x[i];
      else
        rval += x[i];
    }

    if (count)
      rval /= count;

    return rval;
  }

  // book keeping array so that any index can be removed at most once
  // if provided then use only those indices
  if (iuse_idx)
    r_idx = iuse_idx;
  else
  {
    r_idx = mdi_Create(1,length);
    for(i=0; i<length;i++)
      MdiV0(r_idx,i) = 1;
  }

  if (ig_idx)
  {
    for (k=0; k<MNR(ig_idx)*MNC(ig_idx); k++)
    {
      if ((mdiV0(ig_idx,k) < 0) || (mdiV0(ig_idx,k) > length-1 ))
        continue;
      MdiV0(r_idx,mdiV0(ig_idx,k))=0;
    }
  }

  for (i=0; i<length; i++)
  {
    if(isnand (x[i]))
      continue;

    if((x[i]==R_INF) || (x[i]==-R_INF))
    {
      if (ig_infs)
        continue;
      else
      {
        rval = R_INF;
        break;
      }
    }

    if (mdiV0(r_idx,i)==0)
      continue;

    count ++;

    if (isnand(rval))
      rval = x[i];
    else
      rval += x[i];
  }

  if (!iuse_idx)
    mdr_Destroy(r_idx);

  if (count)
    rval /= count;
  return rval;
}


// find a mean of vector - ignore nan's
void rlab_meanstat_from_valstat_vectors(int length, double *x, double *w, MDR *ig_idx, MDR * iuse_idx,
                                double *m1, double *v1, int do_std)
{
  int i, k;
  double sw=0, sw2=0, dw, m=create_nan(), v=create_nan();
  MDR * r_idx=0;

  if (length==0)
    return;

  // book keeping array so that any index can be removed at most once
  // if provided then use only those indices
  if (iuse_idx)
  {
    r_idx = mdi_Create(1,length);
    mdr_Zero(r_idx);
    for (i=0; i<SIZE(iuse_idx); i++)
    {
      if (mdiV0(iuse_idx,i)<0 || mdiV0(iuse_idx,i)>length-1)
        continue;
      MdiV0(r_idx,mdiV0(iuse_idx,i))=1;
    }
  }
  else
  {
    r_idx = mdi_Create(1,length);
    for(i=0; i<length;i++)
      MdiV0(r_idx,i) = 1;
  }

  // check 'x' for nans
  for(i=0; i<length;i++)
  {
    if(isnand(x[i]) || isnand(w[i]))
      MdiV0(r_idx,i)=0;
  }

  if (ig_idx)
  {
    for (k=0; k<MNR(ig_idx)*MNC(ig_idx); k++)
    {
      if ((MdiV0(ig_idx,k) < 0) || (MdiV0(ig_idx,k) > length-1 ))
        continue;
      MdiV0(r_idx,mdiV0(ig_idx,k))=0;
    }
  }

  for (i=0; i<length; i++)
  {
    if (!MdiV0(r_idx,i))
      continue;

    if (do_std)
      dw = 1/(w[i] * w[i]);
    else
      dw = w[i];


      // update mean
    if (isnand(m))
    {
      m   = dw * x[i];
      sw  = dw;
      sw2 = dw * dw;
    }
    else
    {
      m   += dw * x[i];
      sw  += dw;
      sw2 += dw * dw;
    }
      // update variance
    if (isnand(v))
      v  = dw * x[i] * x[i];
    else
      v += dw * x[i] * x[i];
  }

  mdr_Destroy(r_idx);

  m /= sw;
  v  = (v/sw - m*m) * (sw*sw)/((sw*sw) - sw2);

  if (do_std)
  {
    if (v > 0)
      v = sqrt(v);
    else
      v = create_nan();
  }
  else
    v = 1/v;

  *m1 = m;
  *v1 = v;

  return;
}

// find a mean of vector - ignore nan's
double rlab_var_vector(int length, double *x, int bias, MDR *ig_idx, MDR * iuse_idx)
{
  double rval = create_nan(), rval2=create_nan();
  int i, k;
  double count=0;
  MDR * r_idx=0;

  if (length==0)
    return rval;

  // book keeping array so that any index can be removed at most once
  // if provided then use only those indices
  if (iuse_idx)
    r_idx = iuse_idx;
  else
  {
    r_idx = mdi_Create(1,length);
    for(i=0; i<length;i++)
      MdiV0(r_idx,i) = 1;
  }

  if (ig_idx)
  {
    for (k=0; k<MNR(ig_idx)*MNC(ig_idx); k++)
    {
      if ((mdiV0(ig_idx,k) < 0) || (mdiV0(ig_idx,k) > length-1 ))
        continue;
      MdiV0(r_idx,mdiV0(ig_idx,k))=0;
    }
  }

  for (i=0; i<length; i++)
  {
    if(isnand (x[i]))
      continue;

    if (mdiV0(r_idx,i)==0)
      continue;

    count += 1.0;
    if (isnand(rval))
    {
      rval  = x[i];
      rval2 = x[i] * x[i];
    }
    else
    {
      rval  += x[i];
      rval2 += x[i] * x[i];
    }
  }

  if (!iuse_idx)
    mdr_Destroy(r_idx);

  if (count)
    rval = rval2 / (count-bias) - (rval*rval) / (count-bias) / count;

  return rval;
}

// // find a var of matrix - ignore nan's
// double rlab_var_vector(int length, double *x, int bias, int lig_idx, int *ig_idx)
// {
//   double rval  = create_nan();
//   double rval2 = create_nan();
//   int i, k;
//   double count=0;
//
//   if (length >0)
//   {
//     for (i=0; i<length; i++)
//     {
//       if (lig_idx)
//       {
//         for (k=0; k<lig_idx; k++)
//         {
//           if( ig_idx[k]==i)
//             goto continue_loop;
//         }
//       }
//
//       if(isnand (x[i]))
//         continue;
//
//       if (isnand(rval))
//       {
//         rval  = x[i];
//         rval2 = x[i] * x[i];
//         count += 1.0;
//         continue;
//       }
//
//       rval  += x[i];
//       rval2 += x[i] * x[i];
//       count += 1.0;
//
// continue_loop:
//       continue;
//     }
//
//     if (count)
//       rval = rval2 / (count-bias) - (rval*rval) / (count-bias) / count;
//   }
//
//   return rval;
// }

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

      if ( MdiV0(r_idx, ig_idx[i]) )
      {
        MdiV0(r_idx,ig_idx[i]) = 0; // to prevent tards from having ig_idx=[ ..j .. j ] to remove same element twice
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












