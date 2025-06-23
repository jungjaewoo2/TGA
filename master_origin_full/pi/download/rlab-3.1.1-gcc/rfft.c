/* rfft.c: FFT, and related functions. */

/*  This file is a part of RLaB ("Our"-LaB)
   Copyright (C) 1995  Ian R. Searle
   Copyright (C) 2008-2020  Marijan Kostrun

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
#include "mem.h"
#include "symbol.h"
#include "bltin.h"
#include "mdr.h"
#include "mdc.h"
#include "util.h"
#include "ent.h"
#include "class.h"
#include "rfileio.h"

#include "blas.h"
#include "mathl.h"

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>

/*
 * Header files from the available classes...
 */

#include "ent.h"
#include "btree.h"
#include "bltin.h"
#include "function.h"
#include "mdr.h"
#include "mdrf1.h"
#include "mdrf2.h"
#include "mdc.h"
#include "mdcf1.h"
#include "mdcf2.h"
#include "mdr_mdc.h"
#include "mds.h"
#include "mdsf1.h"
#include "mdr_mds.h"
#include "msr.h"
#include "msrf1.h"
#include "msc.h"
#include "mscf1.h"
#include "mathl.h"
#include "complex.h"

#include "fftp.h"

#include <gsl/gsl_version.h>
#include "gsl_rlab.h"
#if ((GSL_MAJOR_VERSION > 1) && (GSL_MINOR_VERSION > 4))
#include <gsl/gsl_filter.h>
#endif

#include "rlab_solver_parameters_names.h"

static OpDef real_method[NCL];
static OpDef imag_method[NCL];
static OpDef conj_method[NCL];
static OpDef fft_method[NCL];
static OpDef ifft_method[NCL];

MDC *mdr_FFT_BF (MDR * m, int nfft);
MDC *mdc_FFT_BF (MDC * m, int nfft);

MDC *mdr_IFFT_BF (MDR * m, int nfft);
MDC *mdc_IFFT_BF (MDC * m, int nfft);

/* **************************************************************
 * Initialize the built-ins...
 * ************************************************************** */

void
class_fft_init (void)
{
  /*
   * Real ()
   */

  real_method[MATRIX_DENSE_REAL].type = MATRIX_DENSE_REAL;
  real_method[MATRIX_DENSE_REAL].op = (void *) mdr_Real_BF;

  real_method[MATRIX_DENSE_COMPLEX].type = MATRIX_DENSE_REAL;
  real_method[MATRIX_DENSE_COMPLEX].op = (void *) mdc_Real_BF;

  real_method[MATRIX_SPARSE_REAL].type = MATRIX_SPARSE_REAL;
  real_method[MATRIX_SPARSE_REAL].op = (void *) msr_Real_BF;

  real_method[MATRIX_SPARSE_COMPLEX].type = MATRIX_SPARSE_REAL;
  real_method[MATRIX_SPARSE_COMPLEX].op = (void *) msc_Real_BF;

  /*
   * Imag ()
   */

  imag_method[MATRIX_DENSE_REAL].type = MATRIX_DENSE_REAL;
  imag_method[MATRIX_DENSE_REAL].op = (void *) mdr_Imag_BF;

  imag_method[MATRIX_DENSE_COMPLEX].type = MATRIX_DENSE_REAL;
  imag_method[MATRIX_DENSE_COMPLEX].op = (void *) mdc_Imag_BF;

  imag_method[MATRIX_SPARSE_REAL].type = MATRIX_SPARSE_REAL;
  imag_method[MATRIX_SPARSE_REAL].op = (void *) msr_Imag_BF;

  imag_method[MATRIX_SPARSE_COMPLEX].type = MATRIX_SPARSE_REAL;
  imag_method[MATRIX_SPARSE_COMPLEX].op = (void *) msc_Imag_BF;

  /*
   * Conj ()
   */

  conj_method[MATRIX_DENSE_REAL].type = MATRIX_DENSE_REAL;
  conj_method[MATRIX_DENSE_REAL].op = (void *) mdr_Conj_BF;

  conj_method[MATRIX_DENSE_COMPLEX].type = MATRIX_DENSE_COMPLEX;
  conj_method[MATRIX_DENSE_COMPLEX].op = (void *) mdc_Conj_BF;

  conj_method[MATRIX_SPARSE_REAL].type = MATRIX_SPARSE_REAL;
  conj_method[MATRIX_SPARSE_REAL].op = (void *) msr_Conj;

  conj_method[MATRIX_SPARSE_COMPLEX].type = MATRIX_SPARSE_COMPLEX;
  conj_method[MATRIX_SPARSE_COMPLEX].op = (void *) msc_Conj;

  /*
   * fft ()
   */

  fft_method[MATRIX_DENSE_REAL].type = MATRIX_DENSE_COMPLEX;
  fft_method[MATRIX_DENSE_REAL].op = (void *) mdr_FFT_BF;

  fft_method[MATRIX_DENSE_COMPLEX].type = MATRIX_DENSE_COMPLEX;
  fft_method[MATRIX_DENSE_COMPLEX].op = (void *) mdc_FFT_BF;

  /*
   * ifft ()
   */

  ifft_method[MATRIX_DENSE_REAL].type = MATRIX_DENSE_COMPLEX;
  ifft_method[MATRIX_DENSE_REAL].op = (void *) mdr_IFFT_BF;

  ifft_method[MATRIX_DENSE_COMPLEX].type = MATRIX_DENSE_COMPLEX;
  ifft_method[MATRIX_DENSE_COMPLEX].op = (void *) mdc_IFFT_BF;

}

/* **************************************************************
 * Return the REAL part of an Entity.
 * ************************************************************** */

Ent *
Real (int nargs, Datum args[])
{
  int type, rtype;
  void *rval;
  void *(*vfptr) ();
  Ent *rent, *e;

  if (nargs != 1)
  {
    rerror ("real: 1 argument allowed");
  }

  e = bltin_get_ent (args[0]);
  type = ent_type (e);

  rtype = real_method[type].type;
  vfptr = (VFPTR) real_method[type].op;

  if (vfptr == 0)
  {
    fprintf (stderr, "Entity type: %s\n", etd (e));
    rerror ("real() operation not supported");
  }

  rval = (*vfptr) (ent_data (e));

  /*
   * Special thing here. Check that the return data is
   * not the same as the entities data. If it is, then
   * return the original entitiy.
   */

  if (rval == ent_data (e))
  {
    return (e);
  }

  rent = ent_Create ();
  ent_data (rent) = rval;
  ent_type (rent) = rtype;

  ent_Clean (e);
  return (rent);
}

/* **************************************************************
 * Return the Imaginary part of an Entity.
 * ************************************************************** */

Ent *
Imag (int nargs, Datum args[])
{
  int type, rtype;
  void *rval;
  void *(*vfptr) ();
  Ent *rent, *e;

  if (nargs != 1)
  {
    rerror ("imag: 1 argument allowed");
  }

  e = bltin_get_ent (args[0]);
  type = ent_type (e);

  rtype = imag_method[type].type;
  vfptr = (VFPTR) imag_method[type].op;

  if (vfptr == 0)
  {
    fprintf (stderr, "Entity type: %s\n", etd (e));
    rerror ("imag() operation not supported");
  }

  rval = (*vfptr) (ent_data (e));

  /*
   * Special thing here. Check that the return data is
   * not the same as the entities data. If it is, then
   * return the original entitiy.
   */

  if (rval == ent_data (e))
  {
    return (e);
  }

  rent = ent_Create ();
  ent_data (rent) = rval;
  ent_type (rent) = rtype;

  ent_Clean (e);
  return (rent);
}

/* **************************************************************
 * Return the Conjugate of an Entity.
 * ************************************************************** */

Ent *
Conj (int nargs, Datum args[])
{
  int type, rtype;
  void *rval;
  void *(*vfptr) ();
  Ent *rent, *e;

  if (nargs != 1)
  {
    rerror ("conj: 1 argument allowed");
  }

  e = bltin_get_ent (args[0]);
  type = ent_type (e);

  rtype = conj_method[type].type;
  vfptr = (VFPTR) conj_method[type].op;

  if (vfptr == 0)
  {
    fprintf (stderr, "Entity type: %s\n", etd (e));
    rerror ("conj() operation not supported");
  }

  rval = (*vfptr) (ent_data (e));

  /*
   * Special thing here. Check that the return data is
   * not the same as the entities data. If it is, then
   * return the original entitiy.
   */

  if (rval == ent_data (e))
  {
    return (e);
  }

  rent = ent_Create ();
  ent_data (rent) = rval;
  ent_type (rent) = rtype;

  ent_Clean (e);
  return (rent);
}

/* **************************************************************
 * Discrete Fourier Transform.
 * ************************************************************** */


Ent *
FFT (int nargs, Datum args[])
{
  int n, type, rtype;
  void *rval;
  void *(*vfptr) ();
  Ent *rent, *e, *e2;

  if ((nargs < 1) || (nargs > 2))
  {
    rerror ("fft: 1 or 2 arguments allowed");
  }

  e = bltin_get_ent (args[0]);
  type = ent_type (e);

  rtype = fft_method[type].type;
  vfptr = (VFPTR) fft_method[type].op;

  if (vfptr == 0)
  {
    fprintf (stderr, "Entity type: %s\n", etd (e));
    rerror ("fft() operation not supported");
  }

  if (nargs == 2)
  {
    e2 = bltin_get_ent (args[1]);
    n = (int) class_double (e2);
  }
  else
  {
    n = 0;			/* Default */
  }

  rval = (*vfptr) (ent_data (e), n);

  rent = ent_Create ();
  ent_data (rent) = rval;
  ent_type (rent) = rtype;

  ent_Clean (e);
  return (rent);
}

/* **************************************************************
 * Real fft, just coerce to COMPLEX.
 * ************************************************************** */

MDC *
mdr_FFT_BF (MDR * m, int N)
{
  MDC *mfft, *tmp;
  mfft = mdr_coerce_mdc (m);
  tmp = mdc_FFT_BF (mfft, N);
  mdc_Destroy (mfft);
  return (tmp);
}

/* **************************************************************
 * Complex FFT.
 * ************************************************************** */

MDC *
mdc_FFT_BF (MDC * m, int N)
{
  F_INT n;
  int i, vtype;
  MDR *wsave;
  MDC *mfft;

  if (MNR (m) == 1 || MNC (m) == 1)
  {
    /*
     * Vector fft
     */

    if (MNR (m) > MNC (m))	/* Set n, and figure out vector type */
    {
      n = (F_INT) MNR (m);
      vtype = 1;
    }
    else
    {
      n = (F_INT) MNC (m);
      vtype = 2;
    }

    mfft = mdc_Copy (m);

    if (N != 0)
    {
      /* Set length of vector to N */
      if (N > (int) n)
      {
        if (vtype == 1)   /* column vector */
          mdc_Extend (mfft, N, 1);
        else      /* row vector */
          mdc_Extend (mfft, 1, N);
        n = (F_INT) N;
      }
      else if (N < (int) n)
      {
        if (vtype == 1)   /* column vector */
          mdc_Truncate (mfft, N, 1);
        else      /* row vector */
          mdc_Truncate (mfft, 1, N);
        n = (F_INT) N;
      }
    }

    wsave = mdr_Create (1, 4 * n + 15);

    signal (SIGINT, intcatch);
    CFFTI (&n, MDRPTR (wsave));
    CFFTF (&n, (double *) MDCPTR (mfft), MDRPTR (wsave));
    signal (SIGINT, intcatch_wait);
  }
  else
  {
    /* Matrix fft, one column at a time */
    n = (F_INT) MNR (m);

    mfft = mdc_Copy (m);

    if (N != 0)
    {
      /* Set row dimension of [mfft] to N */
      if (N > (int) n)
      {
        mdc_Extend (mfft, N, MNC (mfft));
        n = (F_INT) N;
      }
      else if (N < (int) n)
      {
        mdc_Truncate (mfft, N, MNC (mfft));
        n = (F_INT) N;
      }
    }

    wsave = mdr_Create (1, 4 * n + 15);
    for (i = 0; i < MNC (m); i++)
    {
      signal (SIGINT, intcatch);
      CFFTI (&n, MDRPTR (wsave));
      CFFTF (&n, (double *) (MDCPTR(mfft) + (i * n)), MDRPTR (wsave));
      signal (SIGINT, intcatch_wait);
    }
  }

  mdr_Destroy (wsave);
  return (mfft);
}

/* **************************************************************
 * Inverse Discrete Fourier Transform
 * ************************************************************** */

Ent *
IFFT (int nargs, Datum args[])
{
  int n, type, rtype;
  void *rval;
  void *(*vfptr) ();
  Ent *rent, *e, *e2;

  if ((nargs < 1) || (nargs > 2))
  {
    rerror ("ifft: 1 or 2 arguments allowed");
  }

  e = bltin_get_ent (args[0]);
  type = ent_type (e);

  rtype = ifft_method[type].type;
  vfptr = (VFPTR) ifft_method[type].op;

  if (vfptr == 0)
  {
    fprintf (stderr, "Entity type: %s\n", etd (e));
    rerror ("ifft() operation not supported");
  }

  if (nargs == 2)
  {
    e2 = bltin_get_ent (args[1]);
    n = (int) class_double (e2);
  }
  else
  {
    n = 0;			/* Default */
  }

  rval = (*vfptr) (ent_data (e), n);

  rent = ent_Create ();
  ent_data (rent) = rval;
  ent_type (rent) = rtype;

  ent_Clean (e);
  return (rent);
}

MDC *
mdr_IFFT_BF (MDR * m, int N)
{
  MDC *mfft, *tmp;
  mfft = mdr_coerce_mdc (m);
  tmp = mdc_IFFT_BF (mfft, N);
  mdc_Destroy (mfft);
  return (tmp);
}

MDC *
mdc_IFFT_BF (MDC * m, int N)
{
  F_INT n;
  int i, vtype;
  MDR *wsave;
  MDC *mfft;

  if (MNR (m) == 1 || MNC (m) == 1)
  {
    /*
    * Vector fft
    */

    if (MNR (m) > MNC (m))  /* Set n, and figure out vector type */
    {
      n = (F_INT) MNR (m);
      vtype = 1;
    }
    else
    {
      n = (F_INT) MNC (m);
      vtype = 2;
    }

    mfft = mdc_Copy (m);

    if (N != 0)
    {
      /* Set length of vector to N */
      if (N > (int) n)
      {
        if (vtype == 1)   /* column vector */
          mdc_Extend (mfft, N, 1);
        else      /* row vector */
          mdc_Extend (mfft, 1, N);
        n = (F_INT) N;
      }
      else if (N < (int) n)
      {
        if (vtype == 1)   /* column vector */
          mdc_Truncate (mfft, N, 1);
        else      /* row vector */
          mdc_Truncate (mfft, 1, N);
        n = (F_INT) N;
      }
    }

    wsave = mdr_Create (1, 4 * n + 15);

    signal (SIGINT, intcatch);
    CFFTI (&n, MDRPTR (wsave));
    CFFTB (&n, (double *) MDCPTR (mfft), MDRPTR (wsave));
    signal (SIGINT, intcatch_wait);
  }
  else
  {
    /* Matrix ifft, one column at a time */
    n = MNR (m);

    mfft = mdc_Copy (m);

    if (N != 0)
    {
      /* Set row dimension of [mfft] to N */
      if (N > (int) n)
      {
        mdc_Extend (mfft, N, MNC (mfft));
        n = (F_INT) N;
      }
      else if (N < (int) n)
      {
        mdc_Truncate (mfft, N, MNC (mfft));
        n = (F_INT) N;
      }
    }

    wsave = mdr_Create (1, 4 * n + 15);

    for (i = 0; i < MNC (m); i++)
    {
      signal (SIGINT, intcatch);
      CFFTI (&n, MDRPTR (wsave));
      CFFTB (&n, (double *) (MDCPTR(mfft) + (i * n)), MDRPTR (wsave));
      signal (SIGINT, intcatch_wait);
    }
  }

  /*
  * Scale the result by 1/N
  */

  for (i = 0; i < MNR (mfft) * MNC (mfft); i++)
    MdcV0(mfft, i) = MdcV0(mfft, i) / ((double) n);;

  mdr_Destroy (wsave);
  return (mfft);
}


/*
 * This is the function that acutally computes the
 * response. There is no error checking done here,
 * because Filter() does it all.
 */
static int
complex_filter (Complex *b, Complex *a, Complex *x, Complex *y,
                Complex *vi, int *ntotal, int *N)
{
  int k, n;
  Complex cdummy, cdummy2;

  for (n = 0; n < *ntotal; n++)
  {
    // y[0] = b[0] * x[n] + vi[1];
    cdummy = b[0] * x[n];
    y[n] = cdummy + vi[1];
    for (k = 1; k < *N; k++)
    {
      // vi[k] = b[k] * x[n] - a[k] * y[n] + vi[k + 1];
      cdummy  = b[k] * x[n];
      cdummy2 = a[k] * y[n];
      vi[k]  = cdummy - cdummy2 + vi[k + 1];
    }
  }

  return 0;
}

static int
real_filter (double *b, double *a, double *x, double *y,
        double *vi, int *ntotal, int *N)
{
  int k, n;
  for (n = 0; n < *ntotal; n++)
  {
    y[n] = b[0] * x[n] + vi[1];
    for (k = 1; k < *N; k++)
    {
      vi[k] = b[k] * x[n] - a[k] * y[n] + vi[k + 1];
    }
  }

  return 0;
}


/* **************************************************************
 * Filter:
 * ************************************************************** */
#define THIS_SOLVER "filter"
Ent * Filter (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *e3=0, *e4=0, *rent;

  int ntotal,NN;

  int N, M, NC, k, j, i , isper=0;
  int rtype_x=UNDEF, rtype2=UNDEF, rtype3=UNDEF, rtype4=UNDEF;
  MDR  *a=0,  *b=0,  *x=0, *y=0, *rzi=0, *ry0=0;
  MDC *ca=0, *cb=0, *cx=0;
  void *y0=0, *zi=0;
  

#if (GSL_MAJOR_VERSION > 1) && (GSL_MINOR_VERSION > 4)
  // GSL:
  MDR *xmedian=0, *xsigma=0, *xoutlier=0;
  gsl_filter_end_t fendt = GSL_FILTER_END_PADZERO;
  gsl_filter_scale_t fscalet = GSL_FILTER_SCALE_MAD;
  int use_gsl=0, smedian=0, rmedian=0, impulse;
  double t_scale=0;
  size_t noutlier=0;
#endif

  ListNode *node=0;

  // Check n_args
  if (nargs != 2 && nargs != 3 && nargs!=4)
    rerror (THIS_SOLVER ": requires two, three or four arguments");

  //
  // Get args from list and error-check.
  //
  e1 = bltin_get_ent (args[0]);
  if (e1)
    rtype_x = ent_type(e1);
  if ( rtype_x==MATRIX_DENSE_REAL)
  {
    x  = class_matrix_real (e1);
  }
  else if (rtype_x==MATRIX_DENSE_COMPLEX )
  {
    cx = ent_data(e1);
  }
  else
    rerror (THIS_SOLVER ": First argument 'x' must be dense real or complex vector or matrix!");

  e2 = bltin_get_ent (args[1]);
  rtype2 = ent_type(e2);
  if (rtype2 == BTREE)
  {
    // get b:
    node = btree_FindNode (ent_data (e2), RLAB_NAME_FILTER_B);
    if (node)
    {
      rtype2 = ent_type(var_ent(node));
      if (rtype2 == MATRIX_DENSE_REAL)
      {
        b = ent_data(var_ent(node));
      }
      else if(rtype2 == MATRIX_DENSE_COMPLEX)
      {
        b = ent_data(var_ent(node));
      }
    }

    // get a:
    node = btree_FindNode (ent_data (e2), RLAB_NAME_FILTER_A);
    if (node)
    {
      rtype3 = ent_type(var_ent(node));
      if (rtype3 == MATRIX_DENSE_REAL || rtype3 == MATRIX_DENSE_COMPLEX)
        y0 = ent_data(var_ent(node));
    }

    // get y0:
    node = btree_FindNode (ent_data (e2), RLAB_NAME_FILTER_Y0);
    if (node)
    {
      rtype4 = ent_type(var_ent(node));
      if (rtype4 == MATRIX_DENSE_REAL || rtype4 == MATRIX_DENSE_COMPLEX)
        y0 = ent_data(var_ent(node));
    }

    // get zi
    node = btree_FindNode (ent_data (e2), RLAB_NAME_FILTER_ZI);
    if (node)
    {
      rtype4 = ent_type(var_ent(node));
      if (rtype4 == MATRIX_DENSE_REAL || rtype4 == MATRIX_DENSE_COMPLEX)
        zi = ent_data(var_ent(node));
    }
    node = btree_FindNode (ent_data (e2), RLAB_NAME_FILTER_PERIODIC);
    if (node)
      isper = 1;

    // user cannot provide more then one parameter
    if ((y0 != 0) + (zi !=0) + isper > 1)
      rerror (THIS_SOLVER ": Cannot provide more then one option 'zi', or 'y0' or 'periodic'\n");

#if ((GSL_MAJOR_VERSION > 1) && (GSL_MINOR_VERSION > 4))
    //
    // use GSL robust filters
    //
    // standard median filter
    node = btree_FindNode (ent_data (e2), RLAB_NAME_FILTER_GSL_SMEDIAN);
    if (node)
    {
      if (ent_type(var_ent(node)) == MATRIX_DENSE_REAL)
      {
        use_gsl = 1;
        smedian = (int) class_double(var_ent(node));
        if (smedian < 1)
        {
          use_gsl = 0;
          smedian = 0;
        }
      }
    }

    // recursive median filter
    node = btree_FindNode (ent_data (e2), RLAB_NAME_FILTER_GSL_RMEDIAN);
    if (node)
    {
      if (ent_type(var_ent(node)) == MATRIX_DENSE_REAL)
      {
        use_gsl = 1;
        rmedian = (int) class_double(var_ent(node));
        if (rmedian < 1)
        {
          use_gsl = 0;
          rmedian = 0;
        }
      }
    }

    // end points 
    node = btree_FindNode (ent_data (e2), RLAB_NAME_FILTER_GSL_END_T);
    if (node)
    {
      if (ent_type(var_ent(node)) == MATRIX_DENSE_STRING)
      {
        char* desc_end_t = class_char_pointer(var_ent(node));
        if (strstr(desc_end_t, "val")!=0)
        {
          fendt = GSL_FILTER_END_PADVALUE;
        }
        else if (strstr(desc_end_t, "tru")!=0)
        {
          fendt = GSL_FILTER_END_TRUNCATE;
        }
        else
        {
          fendt = GSL_FILTER_END_PADZERO;
        }
      }
    }

    // impulse filter
    node = btree_FindNode (ent_data (e2), RLAB_NAME_FILTER_GSL_IMPULSE);
    if (node)
    {
      if (ent_type(var_ent(node)) == MATRIX_DENSE_REAL)
      {
        use_gsl = 1;
        impulse = (int) class_double(var_ent(node));
        if (impulse < 1)
        {
          use_gsl = 0;
          impulse = 0;
        }
      }
    }

    // impulse filter: MAD
    node = btree_FindNode (ent_data (e2), RLAB_NAME_FILTER_GSL_IMPULSE_SCALE_MAD);
    if (node)
    {
      if (ent_type(var_ent(node)) == MATRIX_DENSE_REAL)
      {
        use_gsl = 1;
        fscalet = GSL_FILTER_SCALE_MAD;
        t_scale = (double) class_double(var_ent(node));
        if (t_scale < 0)
        {
          use_gsl = 0;
          t_scale = 0;
          impulse = 0;
        }
      }
    }

    // impulse filter: IQR
    node = btree_FindNode (ent_data (e2), RLAB_NAME_FILTER_GSL_IMPULSE_SCALE_IQR);
    if (node)
    {
      if (ent_type(var_ent(node)) == MATRIX_DENSE_REAL)
      {
        use_gsl = 1;
        fscalet = GSL_FILTER_SCALE_IQR;
        t_scale = (double) class_double(var_ent(node));
        if (t_scale < 0)
        {
          use_gsl = 0;
          t_scale = 0;
          impulse = 0;
        }
      }
    }

    // impulse filter: SN
    node = btree_FindNode (ent_data (e2), RLAB_NAME_FILTER_GSL_IMPULSE_SCALE_SN);
    if (node)
    {
      if (ent_type(var_ent(node)) == MATRIX_DENSE_REAL)
      {
        use_gsl = 1;
        fscalet = GSL_FILTER_SCALE_SN;
        t_scale = (double) class_double(var_ent(node));
        if (t_scale < 0)
        {
          use_gsl = 0;
          t_scale = 0;
          impulse = 0;
        }
      }
    }

    // impulse filter: QN
    node = btree_FindNode (ent_data (e2), RLAB_NAME_FILTER_GSL_IMPULSE_SCALE_QN);
    if (node)
    {
      if (ent_type(var_ent(node)) == MATRIX_DENSE_REAL)
      {
        use_gsl = 1;
        fscalet = GSL_FILTER_SCALE_QN;
        t_scale = (double) class_double(var_ent(node));
        if (t_scale < 0)
        {
          use_gsl = 0;
          t_scale = 0;
          impulse = 0;
        }
      }
    }
#endif

  }
  else if ( rtype2==MATRIX_DENSE_REAL)
  {
    b = class_matrix_real(e2);
  }
  else if (rtype2==MATRIX_DENSE_COMPLEX)
  {
    cb = ent_data(e2);
  }
  else
    rerror (THIS_SOLVER ": Second argument 'b' must be dense real or complex vector or list <<a;b;y0!");

  if (nargs>=3)
  {
    e3 = bltin_get_ent (args[2]);
    if (e3)
      rtype3 = ent_type(e3);
    if ( rtype3!=MATRIX_DENSE_REAL &&  rtype3!=MATRIX_DENSE_COMPLEX && rtype3!=UNDEF )
      rerror (THIS_SOLVER ": Third argument 'a' must be dense real or complex vector, or omitted!");
  }

  if (nargs==4)
  {
    // fourth argument determines algorithm being used:
    // y0 ->
    // zi -> use difference scheme, may be slower, as it coerces arrays 'a' and 'b' to same length
    // periodic -> assume periodicity of signal, and use circular buffer to filter data
    e4 = bltin_get_ent (args[3]);
    if (e4)
      if (ent_type(e4)==BTREE)
      {
        node = btree_FindNode (ent_data (e4), RLAB_NAME_FILTER_Y0);
        if (node)
        {
          rtype4 = ent_type(var_ent(node));
          if (rtype4 == MATRIX_DENSE_REAL || rtype4 == MATRIX_DENSE_COMPLEX)
            y0 = ent_data(var_ent(node));
        }
        node = btree_FindNode (ent_data (e4), RLAB_NAME_FILTER_ZI);
        if (node)
        {
          rtype4 = ent_type(var_ent(node));
          if (rtype4 == MATRIX_DENSE_REAL || rtype4 == MATRIX_DENSE_COMPLEX)
            zi = ent_data(var_ent(node));
        }
        node = btree_FindNode (ent_data (e4), RLAB_NAME_FILTER_PERIODIC);
        if (node)
          isper = 1;

        // user cannot provide more then one parameter
          if ((y0 != 0) + (zi !=0) + isper > 1)
            rerror (THIS_SOLVER ": Cannot provide more then one option 'zi', or 'y0' or 'periodic'\n");
      }
  }

  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  //
  // difference approach to calculating the filter
  //
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#if (GSL_MAJOR_VERSION > 1) && (GSL_MINOR_VERSION > 4)
  if (use_gsl)
  {
    if (rtype_x!=MATRIX_DENSE_REAL)
      rerror(THIS_SOLVER ": GSL Filters operate only on 'x' real!");

    gsl_vector * gsl_x = alloc_gsl_vector_mdr(x);
    y = mdr_Create (MNR(x), MNC(x));
    gsl_vector * gsl_y = alloc_gsl_vector_mdr(y);

    if (smedian)
    {
      gsl_filter_median_workspace * w = gsl_filter_median_alloc(2 * smedian + 1);
      gsl_filter_median(fendt, gsl_x, gsl_y, w);
      gsl_filter_median_free(w);
    }
    else if (rmedian)
    {
      gsl_filter_rmedian_workspace * w = gsl_filter_rmedian_alloc(2 * rmedian + 1);
      gsl_filter_rmedian(fendt, gsl_x, gsl_y, w);
      gsl_filter_rmedian_free(w);
    }
    else if (impulse)
    {
      xsigma  = mdr_Create (MNR(x), MNC(x));
      gsl_vector * gsl_xsigma = alloc_gsl_vector_mdr(xsigma);
      xmedian = mdr_Create (MNR(x), MNC(x));
      gsl_vector * gsl_xmedian = alloc_gsl_vector_mdr(xmedian);
      MDR * xout1 = mdi_Create (MNR(x), MNC(x));
      gsl_vector_int * gsl_xoutlier = alloc_gsl_vector_int_mdr(xout1);

      gsl_filter_impulse_workspace * w = gsl_filter_impulse_alloc(2 * impulse + 1);
      gsl_filter_impulse(fendt, fscalet, t_scale, gsl_x, gsl_y,
                         gsl_xmedian, gsl_xsigma, &noutlier, gsl_xoutlier, w);
      gsl_filter_impulse_free (w);

      if (noutlier>0)
      {
        int j=0;
        xoutlier = mdr_Create(1, noutlier);
        for (i=0; i<SIZE(x); i++)
        {
          if (MdiV0(xout1,i))
          {
            MdrV0(xoutlier,j) = i+1;
            j++;
          }
        }
      }

      mdr_Destroy(xout1);
      free_gsl_vector_int_mdr(gsl_xoutlier);
      free_gsl_vector_mdr(gsl_xmedian);
      free_gsl_vector_mdr(gsl_xsigma);
    }

    free_gsl_vector_mdr(gsl_x);
    free_gsl_vector_mdr(gsl_y);

    ent_Clean(e1);
    ent_Clean(e2);
    ent_Clean(e3);
    ent_Clean(e4);

    // write up results for user
    Btree *bt = btree_Create ();

    Ent *ey = ent_Create ();
    ent_SetType (ey, MATRIX_DENSE_REAL);
    ent_data (ey) = y;
    install (bt, "y", ey);

    if (xsigma)
    {
      install (bt, "sigma", ent_Assign_Rlab_MDR(xsigma));
    }

    if (xmedian)
    {
      install (bt, "median", ent_Assign_Rlab_MDR(xmedian));
    }

    if (noutlier>0)
    {
      install (bt, "outliers", ent_Assign_Rlab_MDR(xoutlier));
    }


    rent = ent_Create ();
    ent_SetType (rent, BTREE);
    ent_data (rent) = bt;
    return (rent);
  }
  else
#endif
    if (    rtype_x==MATRIX_DENSE_REAL
      &&  rtype2==MATRIX_DENSE_REAL
      && (rtype3==MATRIX_DENSE_REAL || rtype3==UNDEF)
      && (rtype4==MATRIX_DENSE_REAL || rtype4==UNDEF)   )
  {
    MDR *aa=0, *bb=0;

    //
    // all arguments are real: use the original Ian's code
    //
    if (SIZE(x)<1)
      rerror (THIS_SOLVER ": First argument 'X' must be real vector or matrix!");
    if (EQVECT(x))
    {
      ntotal = SIZE(x);
      NC = 1;
    }
    else
    {
      ntotal = MNR(x);
      NC = MNC(x);
    }

    if (!b)
      rerror (THIS_SOLVER ": Second argument 'b' must be real or complex vector!");
    M = SIZE(b);
    
    if (rtype3==MATRIX_DENSE_REAL)
    {
      a = class_matrix_real(e3);
      if (SIZE(a)<1)
        rerror (THIS_SOLVER ": Third argument 'a' must be real or complex vector, or undefined!");
      if (MdrV0(a,0)==0)
        rerror (THIS_SOLVER ": Third argument 'a' must have non-zero first entry!");
    }

    if (a)
      N = a->nrow * a->ncol;
    else
      N = 1;

    NN = MAX(M, N);

    // check sizes of zi and y0 arrays:
    if (zi)
      rzi = (MDR *) zi;
    if (rzi)
    {
      if ( (MNR(rzi)*MNC(rzi)!=NN-1) && (MNR(rzi)!=NN-1) )
        rerror (THIS_SOLVER ": mismatch dimension in entry 'zi' in list 'options' !");
    }
    if (y0)
      ry0 = (MDR *) y0;
    if (ry0)
    {
      if ( (MNR(ry0)*MNC(ry0)!=NN-1) && (MNR(ry0)!=NN-1) )
        rerror (THIS_SOLVER ": mismatch dimension in entry 'y0' in list 'options' !");
    }

    // Fix up pole and zero vectors. Make them the same length, this makes
    // filter's job much easier.
    aa = mdr_Create (NN, 1);
    MdrV0 (aa, 0) = 1;
    bb = mdr_Create (NN, 1);
    mdr_Zero (aa);
    mdr_Zero (bb);
    if (a)
    {
      // if 'a' is provided:
      for (i = 0; i < M; i++)
        MdrV0 (bb, i) = mdrV0 (b, i) / mdrV0 (a, 0);
      for (i = 0; i < N; i++)
        MdrV0 (aa, i) = mdrV0 (a, i) / mdrV0 (a, 0);
    }
    else
    {
      // absence of 'a' means a=1
      for (i = 0; i < M; i++)
        MdrV0 (bb, i) = mdrV0 (b, i);
    }

    // Create delay vectors and load inital delays. Add an extra term to vi[] to
    // make filter's job a little easier. This extra term will always be zero.
    double *vi = (double *) GC_malloc (sizeof (double) * (NN + 1));

    // Do the work...
    y = mdr_Create (ntotal, NC);

    MDR *zf = mdr_Create (NN-1, NC);
    for (j=0; j<NC; j++)
    {
      // vi has to be initialized every time
      for (i = 0; i < NN+1; i++)
        vi[i] = 0.0;

      if (rzi)
      {
        //
        // user provided 'zi'
        //
        if (SIZE(rzi) == NN-1)
        {
          for (i=1; i < NN; i++)
            vi[i] = mdrV1 (rzi, i);
        }
        else if (MNR(rzi)==NN-1)
        {
          for (i = 1; i < NN; i++)
            vi[i] = mdr1 (rzi, i, MIN(j+1,MNC(rzi)));
        }
      }
      else if (ry0)
      {
        //
        // user provided 'y0'
        //
        vi[1] = Mdr1(ry0,1,MIN(j+1,MNC(ry0))) - MdrV1 (bb, 1) * Mdr1(x,1,j+1);
        for (i=2; i< NN; i++)
        {
          vi[i] = Mdr1(ry0,i,MIN(j+1,MNC(ry0)));
          for (k=1; k<=i; k++)
            vi[i] -= MdrV1 (bb, k) * Mdr1(x  ,i-k+1,j+1);
          if (M > 1)
          {
            for (k=2; k<=i; k++)
              vi[i] += MdrV1 (aa, k) * Mdr1(ry0,i-k+1,MIN(j+1,MNC(ry0)));
          }
        }

//         for (k=1; k<NN; k++)
//           fprintf(stderr, "vi[%i] = %g\n", k, vi[k]);
      }

      real_filter (MDRPTR(bb), MDRPTR(aa), &Mdr0(x,0,j), &Mdr0(y,0,j), vi, &ntotal, &NN);
      for (i=1; i<NN; i++)
        Mdr1 (zf,i,j+1) = vi[i];
    }

    GC_FREE (vi);

    mdr_Destroy (aa);
    mdr_Destroy (bb);

    ent_Clean(e1);
    ent_Clean(e2);
    ent_Clean(e3);
    ent_Clean(e4);

    // write up results for user
    Btree *bt = btree_Create ();

    Ent *ey = ent_Create ();
    ent_SetType (ey, MATRIX_DENSE_REAL);
    ent_data (ey) = y;
    install (bt, "y", ey);

    Ent *ezf = ent_Create ();
    ent_SetType (ezf, MATRIX_DENSE_REAL);
    ent_data (ezf) = zf;
    install (bt, "zf", ezf);

    rent = ent_Create ();
    ent_SetType (rent, BTREE);
    ent_data (rent) = bt;
    return (rent);
  }

  MDC *caa=0, *cbb=0;
  Complex cdummy;
  Complex *cvi;
  MDC *cy=0, *czf=0;
  MDC *czi=0, *cy0=0;

  //
  // some of the arguments are complex. use difference equation for
  // complex filter
  //
  // some arguments are complex: coerce every argument to
  // a complex matrix

  // x:
  if (rtype_x==MATRIX_DENSE_REAL)
  {
    cx = mdr_coerce_mdc( x );
  }
  else if (rtype_x==MATRIX_DENSE_COMPLEX)
    cx = ent_data(e1);

  if (!cx)
    rerror (THIS_SOLVER ": First argument 'X' must be real vector or matrix!");
  if (MNR(cx)==1 || MNC(cx)==1)
  {
    ntotal = MNR(cx) * MNC(cx);
    NC = 1;
  }
  else
  {
    ntotal = MNR(cx);
    NC = MNC(cx);
  }

  // b:
  if (rtype2==MATRIX_DENSE_REAL)
  {
    b  = class_matrix_real (e2);
    cb = mdr_coerce_mdc ( b );
  }
  else
    cb = ent_data(e2);
  if (!cb)
    rerror (THIS_SOLVER ": Second argument 'b' must be real or complex vector!");

  // a:
  if (nargs>=3)
  {
    if (rtype3==MATRIX_DENSE_REAL)
    {
      a  = class_matrix_real (e3);
      if (SIZE(a)<1)
        rerror (THIS_SOLVER ": Third argument 'a' must be real or complex vector, or undefined!");
      if (MdrV0(a,0)==0)
        rerror (THIS_SOLVER ": Third argument 'a' must have non-zero first entry!");
      ca = mdr_coerce_mdc ( a );
    }
    else if (rtype3==MATRIX_DENSE_COMPLEX)
    {
      ca = ent_data (e3);
      if (SIZE(ca)<1)
        ca = 0;
    }

    if (ca)
      if (cabs(MdcV0(ca,0))==0)
        rerror (THIS_SOLVER ": Third argument 'a' must have non-zero first entry!");
  }

  M  = MNR(cb) * MNC(cb);
  if (ca)
    N = MNR(ca) * MNC(ca);
  else
    N = 1;
  NN = MAX(M, N);

  // Fix up pole and zero vectors. Make them the same length, this makes
  // filter's job much easier.
  caa = mdc_Create (NN, 1);
  mdc_Zero (caa);
  MdcV0(caa,0) = 1;

  cbb = mdc_Create (NN, 1);
  mdc_Zero (cbb);

  if (ca)
  {
    for (i = 0; i < M; i++)
    {
      cdummy = MdcV0 (cb, i) / MdcV0 (ca, 0);
      MdcV0 (cbb, i) = cdummy;
    }
    for (i = 1; i < N; i++)
    {
      cdummy = MdcV0 (ca, i) / MdcV0 (ca, 0);
      MdcV0 (caa, i) = cdummy;
    }
  }
  else
  {
    for (i = 0; i < M; i++)
    {
      MdcV0 (cbb, i) = MdcV0 (cb, i);
    }
  }

  // Create delay vectors and load inital delays. Add an extra term to vi[] to
  // make filter's job a little easier. This extra term will always be zero.
  cvi = (Complex *) GC_malloc_atomic_ignore_off_page (sizeof (Complex) * (NN + 1));

  // Do the work...
  cy  = mdc_Create (ntotal, NC);
  czf = mdc_Create (NN, NC);

  if (zi)
  {
    switch (rtype4)
    {
      case MATRIX_DENSE_COMPLEX:
        czi = (MDC *) zi;
        break;

      case MATRIX_DENSE_REAL:
        rzi = (MDR *) zi;
        czi = mdr_coerce_mdc (rzi);
        break;

      default:
        break;
    }
  }
  else if (y0)
  {
    switch (rtype4)
    {
      case MATRIX_DENSE_COMPLEX:
        cy0 = (MDC *) y0;
        break;

      case MATRIX_DENSE_REAL:
        ry0 = (MDR *) y0;
        cy0 = mdr_coerce_mdc (ry0);
        break;

      default:
        break;
    }
  }

  for (j=0; j<NC; j++)
  {
    // vi has to be initialized every time
    // vi has to be initialized for every column
    for (i = 0; i < NN; i++)
    {
      cvi[i] = 0.0;
    }
    cvi[NN] = 0.0;

    if (zi)
    {
      //
      // user provided 'zi'
      //
      if (MNR(czi) * MNC(czi) == NN-1)
      {
        for (i=1; i < NN; i++)
        {
          cvi[i] = MdcV1(czi, i);
        }
      }
      else if (MNR(czi)==NN-1)
      {
        for (i = 1; i < NN; i++)
        {
          cvi[i] = Mdc1 (czi, i, MIN(j+1,MNC(czi)));
        }
      }
    }
    else if (y0)
    {
      cdummy = MdcV1 (cbb, 1) * Mdc1(cx,1,j+1);
      cvi[1] = Mdc1(cy0,1,MIN(j+1,MNC(cy0))) - cdummy;

      for (i=2; i< NN; i++)
      {
        // vi[i] = Mdr1(ry0,i,MIN(j+1,MNC(ry0)));
        cvi[i] = Mdc1(cy0,i,MIN(j+1,MNC(cy0)));

        for (k=1; k<=i; k++)
        {
          cdummy  = MdcV1 (cbb, k) * Mdc1(cx,i-k+1,j+1);
          cvi[i] -=  cdummy;
        }
        if (M > 1)
        {
          for (k=2; k<=i; k++)
          {
            // vi[i] += MdrV1 (aa, k) * Mdr1(ry0,i-k+1,MIN(j,MNC(ry0)));
            cdummy  = MdcV1 (caa, k) * Mdc1(cy0,i-k+1,MIN(j+1,MNC(cy0)));
            cvi[i] += cdummy;
          }
        }
      }
    }

    complex_filter (MDCPTR(cbb), MDCPTR(caa), &Mdc0(cx,0,j), &Mdc0(cy,0,j), cvi, &ntotal, &NN);
    for (i=1; i<NN; i++)
    {
      // Mdr1 (zf,i,j) = vi[i];
      Mdc1(czf,i,j+1) = cvi[i];
    }
  }

  GC_FREE (cvi);

  mdc_Destroy (caa);
  mdc_Destroy (cbb);
  if(b)
    mdc_Destroy (cb);
  if(a)
    mdc_Destroy (ca);
  if(x)
    mdc_Destroy (cx);
  if (rzi)
    mdc_Destroy (czi);
  if (ry0)
    mdc_Destroy (cy0);

  ent_Clean(e1);
  ent_Clean(e2);
  ent_Clean(e3);
  ent_Clean(e4);

  // write up results for user
  rent = ent_Create ();
  Btree *bt = btree_Create ();
  ent_SetType (rent, BTREE);
  ent_data (rent) = bt;

  Ent *ey = ent_Create ();
  ent_SetType (ey, MATRIX_DENSE_COMPLEX);
  ent_data (ey) = cy;
  install (bt, "y", ey);

  Ent *ezf = ent_Create ();
  ent_SetType (ezf, MATRIX_DENSE_COMPLEX);
  ent_data (ezf) = czf;
  install (bt, "zf", ezf);

  return (rent);
}


//
// snippets from the ANSI C Time-Frequency toolbox by Manuel DAVY et al.
//
#include <gsl/gsl_math.h>

//
// window function
//
enum ctftbx_window_name
{
  UNKNOWN,
  RECTANG,
  HAMMING,
  HANNING,
  NUTTALL,
  BLACKMAN,
  HARRIS,
  BARTLETT,
  BARTHANN,
  PAPOULIS,
  GAUSS,
  PARZEN,
  HANNA,
  SPLINE
};

static int
create_window (enum ctftbx_window_name wtype, int wlen,
               double * param, int nb_param, double * Window)
{
  /*====================================================================*
  * Name of the function : window.c                                    *
  * Author               : Manuel DAVY                                 *
  * Date of creation     : 03 - 06 - 1999                              *
  *--------------------------------------------------------------------*
  * THE ALGORITHM                                                      *
  *                                                                    *
  * creates a window of given shapes and length.                       *
  * possible shapes and name :                                         *
  *     Rectangular -> RECTANG                                         *
  *     Hamming     -> HAMMING                                         *
  *     Hanning     -> HANNING                                         *
  *     Kaiser      -> KAISER            (1 optional parameter)        *
  *     Nuttal      -> NUTTALL                                          *
  *     Blackman    -> BLACKMAN                                       *
  *     Harris      -> HARRIS                                          *
  *     Triangular  -> BARTLETT, TRIANG                                *
  *     Barthann    -> BARTHANN                                        *
  *     Papoulis    -> PAPOULIS                                        *
  *     Gauss       -> GAUSS             (1 optional parameter)        *
  *     Parzen      -> PARZEN                                          *
  *     Hanna       -> HANNA             (1 optional parameter)        *
  *     Dolph (Dolf)-> DOLPH, DOLF                                     *
  *     Nutbess     -> NUTBESS           (2 optional parameters)       *
  *     Spline      -> SPLINE            (1 compulsary and 1 optional  *
  *                                       parameter)                   *
  *====================================================================*
  * INPUT VARIABLES                                                    *
  * Name           |                   role                            *
  * Window_type    | Type of window to compute (one of the choices     *
  *                | above)                                            *
  * Window_length  | Length of the window                              *
  * nb_param       | Number of parameters passed to compute the window *
  *                | used only for certain distances  (0<=nb_param<=2) *
  *--------------------------------------------------------------------*
  * OUTPUT VARIABLES                                                   *
  * Name           |                   role                            *
  * i              | index in the window                               *
  *--------------------------------------------------------------------*
  * INTERNAL VARIABLES                                                 *
  * Name           |                   role                            *
  * window         | vector of length 'Window_length' containing the   *
  *                | computed window                                   *
  *====================================================================*
  * SUBROUTINES USED HERE                                              *
  *--------------------------------------------------------------------*
  * Name   |                                                           *
  * Action |                                                           *
  * Place  |                                                           *
  *====================================================================*/

  int i;
  double dummy;

  /* computation according to the window type */
  switch (wtype)
  {
    case RECTANG:
      for (i = 0; i < wlen; i++)
      { Window[i] = 1; }
      break;

    case HAMMING:
      for (i = 0; i < wlen; i++)
      {
        Window[i] = 0.54 - 0.46 * cos ((2.0 * M_PI * (i + 1.0)) /
            (wlen + 1.0));
      }
      break;

    case HANNING:
      for (i = 0; i < wlen; i++)
      {
        Window[i] = 0.50 - 0.50 * cos ((2.0 * M_PI * (i + 1.0)) /
            (wlen + 1.0));
      }
      break;

    case NUTTALL:
      for (i = 0; i < wlen; i++)
      {
        dummy = ((-(wlen - 1.0)/2.0 + i) * 2.0 * M_PI) / wlen;
        Window[i] = 0.3635819 +
            0.4891775*cos(dummy) +
            0.1363995*cos(2.0*dummy) +
            0.0106411*cos(3.0*dummy) ;
      }
      break;

    case BLACKMAN:
      for (i = 0; i < wlen; i++)
      {
        dummy = ((-(wlen - 1.0)/2.0 + i) * 2.0 * M_PI) / wlen;
        Window[i] = 0.42 + 0.5 * cos(dummy) +
            0.08 * cos (2.0 * dummy);
      }
      break;

    case HARRIS:
      for (i = 0; i < wlen; i++)
      {
        dummy = (2.0 * M_PI * (i + 1.0)) /(wlen + 1.0);
        Window[i] = 0.35875   -
            0.48829 * cos(dummy) +
            0.14128 * cos(2.0*dummy) -
            0.01168 *cos(3.0*dummy);
      }
      break;

    case BARTLETT: /* or case TRIANG */
      for (i = 0; i < wlen; i++)
      {
        dummy = MIN(i + 1.0 , wlen - i);
        Window[i] = (2.0 * dummy) / (wlen + 1.0);
      }
      break;

    case BARTHANN:
      for (i = 0; i < wlen; i++)
      {
        dummy = MIN(i + 1.0 , wlen - i);
        Window[i] = 0.38 * (1.0 - cos ((2.0 * M_PI * (i + 1.0)) /
            (wlen + 1.0)))
            +
            0.48 * dummy / (wlen + 1.0);
      }
      break;

    case PAPOULIS:
      for (i = 0; i < wlen; i++)
      {
        dummy = ((i + 1.0) * M_PI) / (wlen + 1.0);
        Window[i] = sin (dummy);
      }
      break;

    case GAUSS:
     if (nb_param != 1)
      {
        fprintf(stdout, "window : one parameter required for GAUSS window\n");
        return 1;
      }
      for (i = 0; i < wlen; i++)
      {
        dummy = -1.0 + (2.0*i) / (wlen - 1.0);
        Window[i] = exp((dummy * dummy) * log(param[0]));
      }
      break;

    case PARZEN:
    {
      double   temp;
      for (i = 0; i < wlen; i++)
      {
        dummy = ABS(((-(wlen - 1.0)/2.0 + i)) * 2.0 / wlen);
        temp = 2.0 * pow(1.0 - dummy,3);
        Window[i] = MIN(temp - (1.0 - 2.0*dummy) * (1.0 - 2.0*dummy) * (1.0 - 2.0*dummy) , temp );
      }
    }
    break;

    case HANNA:
      if (nb_param != 1)
      {
        fprintf(stdout, "window : one parameter required for HANNA window\n");
        return 1;
      }
      for (i = 0; i < wlen; i++)
      {
        Window[i] = pow(
            sin(((2.0 * i + 1.0) * M_PI) / (2.0 * wlen)),
        2.0 * param[0]);
      }
      break;

    case SPLINE:
    {
      double nfreq, p;
      double inter;

      if ((nb_param != 1) && (nb_param != 2))
      {
        fprintf(stdout, "window: One or two parameters required for SPLINE window\n");
        return 1;
      }

      nfreq = param[0];
      if (nb_param == 2)
      {
        p = param[1];
      }
      else
      {
        p = M_PI * wlen * nfreq / 10.0;
      }


      for (i = 0; i < wlen; i++)
      {
        dummy = -(wlen - 1.0)/2.0 + i;
        inter = (0.5 * nfreq / p) * dummy ;

        if (inter != 0.0)
        {
          inter = sin(inter * M_PI)/(inter * M_PI);
          Window[i] = pow(inter,p);
        }
        else
        {
          Window[i] = 1.0;
        }
      }
    }
    break;


    default :
    {
      printf("window : Unknowm window type\n");
      return 1;
    }
    break;
  }

  return 0;
}

Ent *
ent_ctftbx_window (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *e3=0, *rent;
  MDR *x3=0, *w=0;

  char *wname=0;
  enum ctftbx_window_name wtype=UNKNOWN;
  int nb_param=0, wlen=0;
  double *param=0;


  /* Check n_args */
  if (nargs < 2)
    rerror ("window: requires at least 2 arguments");

  //
  // get the window name
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type(e1)!=MATRIX_DENSE_STRING)
    rerror ("window: first argument 'name' has to be string");
  wname = class_char_pointer(e1);

  //
  // get the window length
  //
  e2 = bltin_get_ent (args[1]);
  if (ent_type(e2)!=MATRIX_DENSE_REAL)
    rerror ("window: second argument 'length' has to be positive scalar");
  wlen = (int) class_double(e2);
  if (wlen <= 0)
    rerror ("window: second argument 'length' has to be positive scalar");

  //
  // get the window parameters
  //
  if(nargs > 2)
  {
    e3 = bltin_get_ent (args[2]);
    if (ent_type(e3)!=MATRIX_DENSE_REAL)
      rerror ("window: third argument are the filter parameters");
    x3 = ent_data(e3);
    nb_param = MIN(x3->nrow * x3->ncol,2);
    param = MDRPTR(x3);
  }

  // figure out what is the window enum:
  if (!strncmp(wname, "rec", 3) || !strncmp(wname, "REC", 3))
  {
    wtype = RECTANG;
  }
  else if (!strncmp(wname, "ham", 3) || !strncmp(wname, "HAM", 3))
  {
    wtype = HAMMING;
  }
  else if (!strncmp(wname, "hanni", 5) || !strncmp(wname, "HANNI", 5))
  {
    wtype = HANNING;
  }
  else if (!strncmp(wname, "nut", 3) || !strncmp(wname, "NUT", 3))
  {
    wtype = NUTTALL;
  }
  else if (!strncmp(wname, "bla", 3) || !strncmp(wname, "BLA", 3))
  {
    wtype = BLACKMAN;
  }
  else if (!strncmp(wname, "har", 3) || !strncmp(wname, "HAR", 3))
  {
    wtype = HARRIS;
  }
  else if (!strncmp(wname, "bartl", 5) || !strncmp(wname, "BARTL", 5) ||
           !strncmp(wname, "tri", 3) || !strncmp(wname, "TRI", 3))
  {
    wtype = BARTLETT;
  }
  else if (!strncmp(wname, "barth", 5) || !strncmp(wname, "BARTH", 5))
  {
    wtype = BARTHANN;
  }
  else if (!strncmp(wname, "pap", 3) || !strncmp(wname, "PAP", 3))
  {
    wtype = PAPOULIS;
  }
  else if (!strncmp(wname, "gau", 3) || !strncmp(wname, "GAU", 3))
  {
    wtype = GAUSS;
  }
  else if (!strncmp(wname, "par", 3) || !strncmp(wname, "PAR", 3))
  {
    wtype = PARZEN;
  }
  else if (!strncmp(wname, "hanna", 5) || !strncmp(wname, "HANNA", 5))
  {
    wtype = HANNA;
  }
  else if (!strncmp(wname, "spl", 3) || !strncmp(wname, "SPL", 3))
  {
    wtype = SPLINE;
  }

  if (wtype != UNKNOWN)
  {
    w = mdr_Create(wlen,1);
    if (create_window (wtype, wlen, param, nb_param, MDRPTR(w)))
    {
      mdr_Destroy(w);
      w = mdr_Create(0,0);
    }
  }
  else
    w = mdr_Create(0,0);

  ent_Clean(e1);
  ent_Clean(e2);
  ent_Clean(e3);

  rent = ent_Create ();
  ent_data (rent) = w;
  ent_type (rent) = MATRIX_DENSE_REAL;
  return (rent);
}

static int
rlabplus_stft (Complex *signal, int slen,
      double  *win, int wlen,
      Complex *tfr, int tfr_ncol_time, int tfr_nrow_freq,
      double  *tfr_time_instants, double *tfr_freq_bins)
{
 /* EXISTS AN INTERFACE PROGRAM TO MATLAB : STFTMEX.C                 *
  *====================================================================*
  * Name of the function : stft.c (void)                             *
  * Author               : Manuel DAVY - IRCyN                         *
  * Date of creation     : 10 - 02 - 1999                              *
  *--------------------------------------------------------------------*
  * THE ALGORITHM                                                   *
  *                                                              *
  *   Given a signal to analyze in time and frequency, computes the      *
  *   Short Time Fourier Transform (STFT) :                            *
  *                                                              *
  *                /            -j 2 pi f s                         *
  *    STFT(t,f) = | x(s)h(t-s)e           ds                      *
  *                /                                             *
  *                                                              *
  *   This function is complex valued. Its computation requires a window,*
  *   its displacement positions and the number of frequency bins to be  *
  *   computed with discrete variables.                              *
  *                                                              *
  *====================================================================*
  *   INPUT VARIABLES                                              *
  *   Name              |                role                          *
  *   signal            |   The signal to analyze. No field modified     *
  *                     |                                              *
  *   win               |   Vector containing the points of the window   *
  *   wlen              |   Number of points of the window (ODD number !)*
  *                     |                                              *
  *   tfr               |  Matrix containing the resulting TFR (complex) *
  *   tfr_time_instants |  positions of the smoothing window             *
  *   tfr_ncol_time        |  length of '.time_instants' = number of cols.  *
  *                     |  in the tfr matrix                             *
  *   tfr_nrow_freq        |  number of frequency bins = number of rows     *
  *                     |  in the tfr matrix                             *
  *   tfr.is_complex    |  must be set to TRUE (a stft is complex !)     *
  *--------------------------------------------------------------------*
  *   OUTPUT VARIABLES                                              *
  *   Name              |                role                            *
  *   tfr.real_part     |  the output tfr matrix  (real_part)            *
  *   tfr.imag_part     |  the output tfr matrix  (imag part)            *
  *   norm_vector       |  Value of the normalization factor applied at  *
  *                     |  the points of computation of the stft i.e.    *
  *                     |  tfr_time_instants                          *
  *--------------------------------------------------------------------*
  *   INTERNAL VARIABLES                                            *
  *   Name              |                 role                          *
  *                     |                                            *
  *  Nfft               | Next power of two to tfr_nrow_freq                *
  *                     |                                            *
  *   column, line      |   variables of displacement in the matrices    *
  *   tau               |   local time variable (the 's' in the equation *
  *                     |   above                                      *
  *   time              |   Current instant of computation of the spectro*
  *   taumin            |   local time variable bounds. Used to take into*
  *   taumax            |   accound the beginning and the end of the     *
  *                     |   signal, where the window is cut               *
  *   normh             |   current normalization value : depends on     *
  *                     |   wether the window is cut or not (near the    *
  *                     |   edges)                                     *
  *   wind_sig_real     |   Real and imaginary parts of the windowed     *
  *   wind_sig_imag     |   signal at the current position of the window *
  *   inter             |   several intermediary variables             *
  *====================================================================*
  *   SUBROUTINES USED HERE                                         *
  *--------------------------------------------------------------------*
  *   Name   | int idx(int line, int row, int nb_row)                    *
  *   Action | computes the vector index for an element in a matrix given*
  *        | the line and column indices and the number of lines       *
  *   Place  | divers.c                                                  *
  *--------------------------------------------------------------------*
  *   Name   | int irem( double x, double y)                             *
  *   Action | computes the remainder after Euclidean division of double *
  *   Place  | divers.c                                                  *
  *--------------------------------------------------------------------*
  * Name   | void fft(int n, int m, double *x, double *y)              *
  *   Action | Computes the fft                                          *
  *   Place  | divers.c                                                  *
  *--------------------------------------------------------------------*
  * Name   | int po2(int x)                                            *
  *   Action | Computes the next power of two of x                       *
  *   Place  | divers.c                                                  *
  *====================================================================*/

  int            col, row, time;
  int            taumin, taumax, tau;
  int            half_wlen, fidx;
  double         normh;
  double         inter;

  MDC *wind_sig=0;   /* windowed signal */

  // Test the input variables
  if (tfr_nrow_freq <= 0)
  {
    fprintf (stdout, "stft : The field tfr_nrow_freq is not correctly set\n");
    return 1;
  }
  if (tfr_ncol_time <= 0)
  {
    fprintf (stdout, "stft : The field tfr_ncol_time is not correctly set\n");
    return 1;
  }

  // checks that the window length is odd
  if (wlen%2 == 0)
  {
    fprintf (stdout, "stft : The window Length must be an ODD number\n");
    return 1;
  }

  half_wlen = (wlen - 1) / 2;
  inter = 0.0;
  for (row = 0; row <wlen; row++)
  {
    inter = inter + win[row] * win[row];
  }
  normh = sqrt (inter);

  // creation of the vector of frequency bins  (output)
  //Nfft = po2 (tfr_nrow_freq);
  //Nfft = tfr_nrow_freq;

  // memory allocation for the windowed signal
  wind_sig = mdc_Create(slen,1);

  // ******************************************************
  //
  // computation of the fft for the current windowed signal
  //
  // ******************************************************

  for (col = 0; col < tfr_ncol_time; col++)
  {
    // initialization of the intermediary vectors
    mdc_Zero(wind_sig);

    // time instants of interest to compute the stft
    time = (int) tfr_time_instants[col];

    // the signal is multipied by the window between the instants
    // time-taumin and time+taumax
    // when the window is wider than the number of desired frequencies (tfr_nrow_freq),
    // the range is limited to tfr_nrow_freq
    taumin = MIN(tfr_nrow_freq / 2, half_wlen);
    taumin = MIN(taumin, time);
    taumax = MIN((tfr_nrow_freq / 2 - 1), half_wlen);
    taumax = MIN(taumax, (slen - time - 1));
    // Computation of a normalization factor,
    // equal to the quadratic norm of the window
    //norm_vector[col] = 1.0 / normh;

    /* The signal is windowed around the current time */
    for (tau = -taumin; tau <= taumax; tau++)
    {
      //used to be (?)  row = ctftbx_irem( (tfr_nrow_freq+tau), tfr_nrow_freq ) ;
      row = time + tau;
      MdcV0(wind_sig,row) = signal[time + tau] * win[half_wlen + tau] / normh;
    }

    // fft of the windowed signal: fft (tfr_nrow_freq, Nfft, wind_sig);
    wind_sig = mdc_FFT_BF (wind_sig, slen);

    // only copy the first 'tfr_nrow_freq' of the fft to tfr matrix
    for (row = 0; row < tfr_nrow_freq; row++)
    {
      fidx = tfr_freq_bins[row];
      tfr[row + col*tfr_nrow_freq] = MdcV0(wind_sig,fidx);
    }
  }

  // write down the frequency bins
  for (row = 0; row < tfr_nrow_freq; row++)
  {
//     tfr_freq_bins[row] = ((double) row) / ((double) tfr_ncol_time);
    tfr_freq_bins[row] = ((double) tfr_freq_bins[row]) / ((double) tfr_ncol_time);
  }

  mdc_Destroy (wind_sig);

  return 0;
}

Ent *
ent_ctftbx_stft (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *e3=0, *rent;
  MDR *win=0, *x3=0;
  MDC *sig=0, *tfr=0;
  MDR *tfr_time_instants=0, *tfr_freq_bins=0;

  Btree *bt=0;

  int slen=0, wlen=0, i=0, ilow=0;
  int tfr_ncol_time=0, tfr_nrow_freq=0;

  /* Check n_args */
  if (nargs<2 || nargs>4)
    rerror ("stft: requires two or three arguments");

  //
  // get the signal
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type(e1)!=MATRIX_DENSE_REAL && ent_type(e1)!=MATRIX_DENSE_COMPLEX)
    rerror ("stft: first argument 'signal' has to be real/complex vector");
  if (ent_type(e1)==MATRIX_DENSE_REAL)
  { sig = mdr_coerce_mdc (ent_data(e1)); }
  else
    sig = ent_data (e1);
  slen = sig->nrow * sig->ncol;
  tfr_ncol_time  = slen;

  //
  // get the window length
  //
  e2 = bltin_get_ent (args[1]);
  if (ent_type(e2)!=MATRIX_DENSE_REAL)
    rerror ("stft: second argument 'window' has to be real vector");
  win = class_matrix_real(e2);
  if (win->nrow !=1 && win->ncol !=1)
    rerror ("stft: second argument 'window' has to be real vector");
  wlen = win->nrow * win->ncol;
  if (wlen >= slen)
    rerror ("stft: length of window should be less than length of signal");

  tfr_nrow_freq = slen;
  if (nargs == 3)
  {
    e3 = bltin_get_ent (args[2]);
    if (ent_type(e3)!=MATRIX_DENSE_REAL)
      rerror ("stft: third argument 'range' has to be Fmax or [Fmin,Fmax]");
    x3 = class_matrix_real (e3);
    if (x3->nrow * x3->ncol == 1)
    {
      tfr_nrow_freq = (int) class_double(e3);
      // silly user
      if (tfr_nrow_freq <= wlen)
        tfr_nrow_freq = wlen + 1;
      else if (tfr_nrow_freq >= slen)
        tfr_nrow_freq = slen;
    }
    else if (x3->nrow * x3->ncol == 2)
    {
      if ((MdrV0(x3,0) != (int) MdrV0(x3,0)) ||
          (MdrV0(x3,1) != (int) MdrV0(x3,1)))
        rerror ("stft: third argument 'range' has to be Fmax or [Fmin,Fmax]");
      if ((MdrV0(x3,0) > MdrV0(x3,1)) ||
          (MdrV0(x3,0) < 0) ||
          (MdrV0(x3,1) > slen))
        rerror ("stft: third argument 'range' has to be Fmax or [Fmin,Fmax]");
      tfr_nrow_freq = MdrV0(x3,1) -  MdrV0(x3,0) + 1;
      ilow = MdrV0(x3,0);
    }
    else
      rerror ("stft: third argument 'range' has to be Fmax or [Fmin,Fmax]");
  }

  tfr = mdc_Create(tfr_nrow_freq,tfr_ncol_time);
  tfr_time_instants = mdr_Create(tfr_ncol_time,1);
  for (i=0; i<tfr_ncol_time;i++)
    MdrV0(tfr_time_instants,i) = i+1;

  // tfr_freq_bins contains a range of indices i of f_i
  // which short time fourier transform we want to obtain
  tfr_freq_bins = mdr_Create(tfr_nrow_freq,1);
  for (i=0; i<tfr_nrow_freq;i++)
  { MdrV0(tfr_freq_bins,i) = ilow + i; }

  rlabplus_stft (MDCPTR(sig), slen, MDRPTR(win), wlen,
        MDCPTR(tfr), tfr_ncol_time, tfr_nrow_freq,
        MDRPTR(tfr_time_instants), MDRPTR(tfr_freq_bins));

  if (ent_type(e1)==MATRIX_DENSE_REAL)
  { mdc_Destroy(sig); }

  ent_Clean(e1);
  ent_Clean(e2);
  ent_Clean(e3);

  rent = ent_Create ();
  bt = btree_Create ();
  ent_SetType (rent, BTREE);
  ent_data (rent) = bt;

  Ent *r1 = ent_Create ();
  ent_SetType (r1, MATRIX_DENSE_REAL);
  ent_data (r1) = tfr_time_instants;
  install (bt, ("time_instants"), r1);

  Ent *r2 = ent_Create ();
  ent_SetType (r2, MATRIX_DENSE_REAL);
  ent_data (r2) = tfr_freq_bins;
  install (bt, ("freq_bins"), r2);

  Ent *r3 = ent_Create ();
  ent_SetType (r3, MATRIX_DENSE_COMPLEX);
  ent_data (r3) = tfr;
  install (bt, ("tfr"), r3);

  return (rent);
}

Ent *
ent_ctftbx_fftshift (int nargs, Datum args[])
{
  Ent *e1=0, *rent;
  MDR *x=0, *w;
  MDC *cx=0, *cw=0;

  int i=0, xlen=0, halflen= 0;
  double d1, d2;
  Complex c1, c2;


  /* Check n_args */
  if (nargs!=1)
    rerror ("fftshift: requires single argument");

  //
  // get the signal
  //
  e1 = bltin_get_ent (args[0]);
  if ((ent_type(e1) != MATRIX_DENSE_REAL) && (ent_type(e1) != MATRIX_DENSE_COMPLEX))
    rerror ("fftshift: first argument must be real or complex vector");

  if (ent_type(e1)==MATRIX_DENSE_REAL)
  {
    x = class_matrix_real(e1);
    if (x->nrow!=1 && x->ncol!=1)
      rerror ("fftshift: first argument must be vector");
    xlen = x->nrow * x->ncol;

    w = mdr_Create(x->nrow, x->ncol);

    /* computation of the half length in case of odd or even length */
    halflen = (int) (xlen/2.0);

    /* case where the length is odd */
    if ((xlen&1)==1)
    {
      d2=MdrV0(x,halflen);
      for (i=0; i<halflen; i++)
      {
        d1 = MdrV0(x,i);
        MdrV0(w,i) = MdrV0(x, halflen+i+1);
        MdrV0(w,halflen + i) = d1;
      }
      MdrV0(w, xlen-1) = d2;
    }
    else
    {
      for (i=0; i<halflen; i++)
      {
        d1 = MdrV0(x,halflen + i);
        MdrV0(w, halflen + i) = MdrV0(x,i);
        MdrV0(w, i) = d1;
      }
    }

    ent_Clean(e1);

    rent = ent_Create ();
    ent_SetType (rent, MATRIX_DENSE_REAL);
    ent_data (rent) = w;
    return (rent);
  }

  //
  // complex
  //

  cx = ent_data(e1);
  if (cx->nrow!=1 && cx->ncol!=1)
    rerror ("fftshift: first argument must be vector");
  xlen = cx->nrow * cx->ncol;


  cw = mdc_Create(cx->nrow, cx->ncol);

  /* computation of the half length in case of odd or even length */
  halflen = (int) (xlen/2.0);

  /* case where the length is odd */
  if ((xlen&1)==1)
  {
    c2=MdcV0(cx,halflen);
    for (i=0; i<halflen; i++)
    {
      c1 = MdcV0(cx,i);
      MdcV0(cw,i) = MdcV0(cx, halflen+i+1);
      MdcV0(cw,halflen + i) = c1;
    }
    MdcV0(cw, xlen-1) = c2;
  }
  else
  {
    for (i=0; i<halflen; i++)
    {
      c1 = MdcV0(cx,halflen + i);
      MdcV0(cw, halflen + i) = MdcV0(cx,i);
      MdcV0(cw, i) = c1;
    }
  }

  ent_Clean(e1);

  rent = ent_Create ();
  ent_SetType (rent, MATRIX_DENSE_COMPLEX);
  ent_data (rent) = cw;
  return (rent);
}


#include "./clibs/ctftbx/libctftbx.h"
int ctftbx_debug=1;

//
//
// time frequency distributions of a signal
//
//

// Born-Jordan
Ent *
ent_ctftbx_tfd_bj (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *e3=0, *e4=0, *rent;
  MDR *winT=0, *winF=0;
  MDC *sig=0, *tfr=0;
  MDR *tfr_time_instants=0, *tfr_freq_bins=0;

  type_signal mysignal;
  type_TFR mytfr;

  int wTlen=0, wFlen=0, i=0;
  int tfr_ncol_time=0, tfr_nrow_freq=0;

  double *time_instants=0;

  Btree *bt=0;

  /* Check n_args */
  if (nargs<3)
    rerror ("tfrd.bj: requires three arguments");

  //
  // get the signal
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type(e1)!=MATRIX_DENSE_REAL && ent_type(e1)!=MATRIX_DENSE_COMPLEX)
    rerror ("tfrd.bj: first argument 'signal' has to be real/complex vector");
  if (ent_type(e1)==MATRIX_DENSE_REAL)
  {
    sig = mdr_coerce_mdc (ent_data(e1));
    mysignal.is_complex = 0;
  }
  else
  {
    sig = mdc_Copy(ent_data (e1));
    mysignal.is_complex = 1;
  }
  mysignal.length = sig->nrow * sig->ncol;
  mysignal.sample_freq = 1;
  time_instants = GC_malloc (mysignal.length * sizeof(double));
  for (i=0;i<mysignal.length;i++)
    time_instants[i]= i+1;
  mysignal.time_instants = time_instants;
  mysignal.signal = MDCPTR(sig);

  //
  // get the windowT length
  //
  e2 = bltin_get_ent (args[1]);
  if (ent_type(e2)!=MATRIX_DENSE_REAL)
    rerror ("tfrd.bj: second argument 'windowT' has to be real vector");
  winT = class_matrix_real(e2);
  if (winT->nrow !=1 && winT->ncol !=1)
    rerror ("tfrd.bj: second argument 'windowT' has to be real vector");
  wTlen = winT->nrow * winT->ncol;
  if (wTlen >= mysignal.length)
    rerror ("tfrd.bj: length of 'windowT' should be less than length of signal");

  //
  // get the windowT length
  //
  e3 = bltin_get_ent (args[2]);
  if (ent_type(e3)!=MATRIX_DENSE_REAL)
    rerror ("tfrd.bj: third argument 'windowF' has to be real vector");
  winF = class_matrix_real(e3);
  if (winF->nrow !=1 && winF->ncol !=1)
    rerror ("tfrd.bj: second argument 'windowF' has to be real vector");
  wFlen = winF->nrow * winF->ncol;
  if (wFlen >= mysignal.length)
    rerror ("tfrd.bj: length of 'windowF' should be less than length of signal");

  // create TFR:
  tfr_ncol_time = mysignal.length;
  tfr_nrow_freq = mysignal.length;

  tfr = mdc_Create(tfr_nrow_freq,tfr_ncol_time);
  tfr_time_instants = mdr_Create(tfr_ncol_time,1);
  for (i=0; i<tfr_ncol_time;i++)
    MdrV0(tfr_time_instants,i) = i+1;
  tfr_freq_bins = mdr_Create(tfr_nrow_freq,1);

  mytfr.N_freq = tfr_nrow_freq;
  mytfr.N_time = tfr_ncol_time;
  mytfr.freq_bins     = MDRPTR(tfr_freq_bins);
  mytfr.time_instants = MDRPTR(tfr_time_instants);
  mytfr.is_complex = 0;
  mytfr.tfr = MDCPTR(tfr);

  ctftbx_bj (mysignal, MDRPTR(winT), wTlen, MDRPTR(winF), wFlen, mytfr);

  // clean-up
  ent_Clean(e1);
  ent_Clean(e2);
  ent_Clean(e3);
  ent_Clean(e4);

  // cleanup of ctftbx variables:
  // mysignal:
  GC_free (time_instants);
  mysignal.time_instants=0;
  mdc_Destroy(sig);
  // mytfr:
  mytfr.freq_bins     = 0;
  mytfr.time_instants = 0;
  mytfr.tfr = 0;


  rent = ent_Create ();
  bt = btree_Create ();
  ent_SetType (rent, BTREE);
  ent_data (rent) = bt;

  Ent *r1 = ent_Create ();
  ent_SetType (r1, MATRIX_DENSE_REAL);
  ent_data (r1) = tfr_time_instants;
  install (bt, ("time_instants"), r1);

  Ent *r2 = ent_Create ();
  ent_SetType (r2, MATRIX_DENSE_REAL);
  ent_data (r2) = tfr_freq_bins;
  install (bt, ("freq_bins"), r2);

  Ent *r3 = ent_Create ();
  ent_SetType (r3, MATRIX_DENSE_COMPLEX);
  ent_data (r3) = tfr;
  install (bt, ("tfr"), r3);

  return (rent);
}

// Butterworth:
Ent *
ent_ctftbx_tfd_butter (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *e3=0, *e4=0, *rent;
  MDR *winT=0, *winF=0;
  MDC *sig=0, *tfr=0;
  MDR *tfr_time_instants=0, *tfr_freq_bins=0;

  type_signal mysignal;
  type_TFR mytfr;

  int wTlen=0, wFlen=0, i=0;
  int tfr_ncol_time=0, tfr_nrow_freq=0;
  double sigma=0;

  double *time_instants=0;

  Btree *bt=0;

  /* Check n_args */
  if (nargs<4)
    rerror ("tfrd.butter: requires four arguments");

  //
  // get the signal
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type(e1)!=MATRIX_DENSE_REAL && ent_type(e1)!=MATRIX_DENSE_COMPLEX)
    rerror ("tfrd.butter: first argument 'signal' has to be real/complex vector");
  if (ent_type(e1)==MATRIX_DENSE_REAL)
  {
    sig = mdr_coerce_mdc (ent_data(e1));
    mysignal.is_complex = 0;
  }
  else
  {
    sig = mdc_Copy(ent_data (e1));
    mysignal.is_complex = 1;
  }
  mysignal.length = sig->nrow * sig->ncol;
  mysignal.sample_freq = 1;
  time_instants = GC_malloc (mysignal.length * sizeof(double));
  for (i=0;i<mysignal.length;i++)
    time_instants[i]= i+1;
  mysignal.time_instants = time_instants;
  mysignal.signal = MDCPTR(sig);

  //
  // get the windowT length
  //
  e2 = bltin_get_ent (args[1]);
  if (ent_type(e2)!=MATRIX_DENSE_REAL)
    rerror ("tfrd.butter: second argument 'windowT' has to be real vector");
  winT = class_matrix_real(e2);
  if (winT->nrow !=1 && winT->ncol !=1)
    rerror ("tfrd.butter: second argument 'windowT' has to be real vector");
  wTlen = winT->nrow * winT->ncol;
  if (wTlen >= mysignal.length)
    rerror ("tfrd.butter: length of 'windowT' should be less than length of signal");

  //
  // get the windowT length
  //
  e3 = bltin_get_ent (args[2]);
  if (ent_type(e3)!=MATRIX_DENSE_REAL)
    rerror ("tfrd.butter: third argument 'windowF' has to be real vector");
  winF = class_matrix_real(e3);
  if (winF->nrow !=1 && winF->ncol !=1)
    rerror ("tfrd.butter: second argument 'windowF' has to be real vector");
  wFlen = winF->nrow * winF->ncol;
  if (wFlen >= mysignal.length)
    rerror ("tfrd.butter: length of 'windowF' should be less than length of signal");

  //
  // get sigma
  //
  e4 = bltin_get_ent (args[3]);
  if (ent_type(e4)!=MATRIX_DENSE_REAL)
    rerror ("tfrd.butter: fourth argument 'sigma' has to be positive scalar");
  sigma = class_double(e4);

  // create TFR:
  tfr_ncol_time = mysignal.length;
  tfr_nrow_freq = mysignal.length;

  tfr = mdc_Create(tfr_nrow_freq,tfr_ncol_time);
  tfr_time_instants = mdr_Create(tfr_ncol_time,1);
  for (i=0; i<tfr_ncol_time;i++)
    MdrV0(tfr_time_instants,i) = i+1;
  tfr_freq_bins = mdr_Create(tfr_nrow_freq,1);

  mytfr.N_freq = tfr_nrow_freq;
  mytfr.N_time = tfr_ncol_time;
  mytfr.freq_bins     = MDRPTR(tfr_freq_bins);
  mytfr.time_instants = MDRPTR(tfr_time_instants);
  mytfr.is_complex = 0;
  mytfr.tfr = MDCPTR(tfr);

  ctftbx_bud (mysignal, MDRPTR(winT), wTlen, MDRPTR(winF), wFlen, sigma, mytfr);

  // clean-up
  ent_Clean(e1);
  ent_Clean(e2);
  ent_Clean(e3);
  ent_Clean(e4);

  // cleanup of ctftbx variables:
  // mysignal:
  GC_free (time_instants);
  mysignal.time_instants=0;
  mdc_Destroy(sig);
  // mytfr:
  mytfr.freq_bins     = 0;
  mytfr.time_instants = 0;
  mytfr.tfr = 0;


  rent = ent_Create ();
  bt = btree_Create ();
  ent_SetType (rent, BTREE);
  ent_data (rent) = bt;

  Ent *r1 = ent_Create ();
  ent_SetType (r1, MATRIX_DENSE_REAL);
  ent_data (r1) = tfr_time_instants;
  install (bt, ("time_instants"), r1);

  Ent *r2 = ent_Create ();
  ent_SetType (r2, MATRIX_DENSE_REAL);
  ent_data (r2) = tfr_freq_bins;
  install (bt, ("freq_bins"), r2);

  Ent *r3 = ent_Create ();
  ent_SetType (r3, MATRIX_DENSE_COMPLEX);
  ent_data (r3) = tfr;
  install (bt, ("tfr"), r3);

  return (rent);
}

// Choi-Williams
Ent *
ent_ctftbx_tfd_cw (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *e3=0, *e4=0, *rent;
  MDR *winT=0, *winF=0;
  MDC *sig=0, *tfr=0;
  MDR *tfr_time_instants=0, *tfr_freq_bins=0;

  type_signal mysignal;
  type_TFR mytfr;

  int wTlen=0, wFlen=0, i=0;
  int tfr_ncol_time=0, tfr_nrow_freq=0;
  double sigma=0;

  double *time_instants=0;

  Btree *bt=0;

  /* Check n_args */
  if (nargs<4)
    rerror ("tfrd.cw: requires four arguments");

  //
  // get the signal
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type(e1)!=MATRIX_DENSE_REAL && ent_type(e1)!=MATRIX_DENSE_COMPLEX)
    rerror ("tfrd.cw: first argument 'signal' has to be real/complex vector");
  if (ent_type(e1)==MATRIX_DENSE_REAL)
  {
    sig = mdr_coerce_mdc (ent_data(e1));
    mysignal.is_complex = 0;
  }
  else
  {
    sig = mdc_Copy(ent_data (e1));
    mysignal.is_complex = 1;
  }
  mysignal.length = sig->nrow * sig->ncol;
  mysignal.sample_freq = 1;
  time_instants = GC_malloc (mysignal.length * sizeof(double));
  for (i=0;i<mysignal.length;i++)
    time_instants[i]= i+1;
  mysignal.time_instants = time_instants;
  mysignal.signal = MDCPTR(sig);

  //
  // get the windowT length
  //
  e2 = bltin_get_ent (args[1]);
  if (ent_type(e2)!=MATRIX_DENSE_REAL)
    rerror ("tfrd.cw: second argument 'windowT' has to be real vector");
  winT = class_matrix_real(e2);
  if (winT->nrow !=1 && winT->ncol !=1)
    rerror ("tfrd.cw: second argument 'windowT' has to be real vector");
  wTlen = winT->nrow * winT->ncol;
  if (wTlen >= mysignal.length)
    rerror ("tfrd.cw: length of 'windowT' should be less than length of signal");

  //
  // get the windowT length
  //
  e3 = bltin_get_ent (args[2]);
  if (ent_type(e3)!=MATRIX_DENSE_REAL)
    rerror ("tfrd.cw: third argument 'windowF' has to be real vector");
  winF = class_matrix_real(e3);
  if (winF->nrow !=1 && winF->ncol !=1)
    rerror ("tfrd.cw: second argument 'windowF' has to be real vector");
  wFlen = winF->nrow * winF->ncol;
  if (wFlen >= mysignal.length)
    rerror ("tfrd.cw: length of 'windowF' should be less than length of signal");

  //
  // get sigma
  //
  e4 = bltin_get_ent (args[3]);
  if (ent_type(e4)!=MATRIX_DENSE_REAL)
    rerror ("tfrd.cw: fourth argument 'sigma' has to be positive scalar");
  sigma = class_double(e4);

  // create TFR:
  tfr_ncol_time = mysignal.length;
  tfr_nrow_freq = mysignal.length;

  tfr = mdc_Create(tfr_nrow_freq,tfr_ncol_time);
  tfr_time_instants = mdr_Create(tfr_ncol_time,1);
  for (i=0; i<tfr_ncol_time;i++)
    MdrV0(tfr_time_instants,i) = i+1;
  tfr_freq_bins = mdr_Create(tfr_nrow_freq,1);  // tfd fills these

  mytfr.N_freq = tfr_nrow_freq;
  mytfr.N_time = tfr_ncol_time;
  mytfr.freq_bins     = MDRPTR(tfr_freq_bins);
  mytfr.time_instants = MDRPTR(tfr_time_instants);
  mytfr.is_complex = 0;
  mytfr.tfr = MDCPTR(tfr);

  ctftbx_cw (mysignal, MDRPTR(winT), wTlen, MDRPTR(winF), wFlen, sigma, mytfr);

  // clean-up
  ent_Clean(e1);
  ent_Clean(e2);
  ent_Clean(e3);
  ent_Clean(e4);

  // cleanup of ctftbx variables:
  // mysignal:
  GC_free (time_instants);
  mysignal.time_instants=0;
  mdc_Destroy(sig);
  // mytfr:
  mytfr.freq_bins     = 0;
  mytfr.time_instants = 0;
  mytfr.tfr = 0;


  rent = ent_Create ();
  bt = btree_Create ();
  ent_SetType (rent, BTREE);
  ent_data (rent) = bt;

  Ent *r1 = ent_Create ();
  ent_SetType (r1, MATRIX_DENSE_REAL);
  ent_data (r1) = tfr_time_instants;
  install (bt, ("time_instants"), r1);

  Ent *r2 = ent_Create ();
  ent_SetType (r2, MATRIX_DENSE_REAL);
  ent_data (r2) = tfr_freq_bins;
  install (bt, ("freq_bins"), r2);

  Ent *r3 = ent_Create ();
  ent_SetType (r3, MATRIX_DENSE_COMPLEX);
  ent_data (r3) = tfr;
  install (bt, ("tfr"), r3);

  return (rent);
}

// Generalized-rectangular
Ent *
ent_ctftbx_tfd_grd (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *e3=0, *e4=0, *e5=0, *rent;
  MDR *winT=0, *winF=0;
  MDC *sig=0, *tfr=0;
  MDR *tfr_time_instants=0, *tfr_freq_bins=0;

  type_signal mysignal;
  type_TFR mytfr;

  int wTlen=0, wFlen=0, i=0;
  int tfr_ncol_time=0, tfr_nrow_freq=0;
  double rs=0, dsym=0;

  double *time_instants=0;

  Btree *bt=0;

  /* Check n_args */
  if (nargs<5)
    rerror ("tfrd.gr: requires five arguments");

  //
  // get the signal
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type(e1)!=MATRIX_DENSE_REAL && ent_type(e1)!=MATRIX_DENSE_COMPLEX)
    rerror ("tfrd.gr: first argument 'signal' has to be real/complex vector");
  if (ent_type(e1)==MATRIX_DENSE_REAL)
  {
    sig = mdr_coerce_mdc (ent_data(e1));
    mysignal.is_complex = 0;
  }
  else
  {
    sig = mdc_Copy(ent_data (e1));
    mysignal.is_complex = 1;
  }
  mysignal.length = sig->nrow * sig->ncol;
  mysignal.sample_freq = 1;
  time_instants = GC_malloc (mysignal.length * sizeof(double));
  for (i=0;i<mysignal.length;i++)
    time_instants[i]= i+1;
  mysignal.time_instants = time_instants;
  mysignal.signal = MDCPTR(sig);

  //
  // get the windowT length
  //
  e2 = bltin_get_ent (args[1]);
  if (ent_type(e2)!=MATRIX_DENSE_REAL)
    rerror ("tfrd.gr: second argument 'windowT' has to be real vector");
  winT = class_matrix_real(e2);
  if (winT->nrow !=1 && winT->ncol !=1)
    rerror ("tfrd.gr: second argument 'windowT' has to be real vector");
  wTlen = winT->nrow * winT->ncol;
  if (wTlen >= mysignal.length)
    rerror ("tfrd.gr: length of 'windowT' should be less than length of signal");

  //
  // get the windowT length
  //
  e3 = bltin_get_ent (args[2]);
  if (ent_type(e3)!=MATRIX_DENSE_REAL)
    rerror ("tfrd.gr: third argument 'windowF' has to be real vector");
  winF = class_matrix_real(e3);
  if (winF->nrow !=1 && winF->ncol !=1)
    rerror ("tfrd.gr: second argument 'windowF' has to be real vector");
  wFlen = winF->nrow * winF->ncol;
  if (wFlen >= mysignal.length)
    rerror ("tfrd.gr: length of 'windowF' should be less than length of signal");

  //
  // get rs
  //
  e4 = bltin_get_ent (args[3]);
  if (ent_type(e4)!=MATRIX_DENSE_REAL)
    rerror ("tfrd.gr: fourth argument 'rs' has to be positive scalar");
  rs = class_double(e4);

  //
  // get dsym
  //
  e5 = bltin_get_ent (args[4]);
  if (ent_type(e5)!=MATRIX_DENSE_REAL)
    rerror ("tfrd.gr: fifth argument 'dsym' has to be positive scalar");
  dsym = class_double(e5);

  // create TFR:
  tfr_ncol_time = mysignal.length;
  tfr_nrow_freq = mysignal.length;

  tfr = mdc_Create(tfr_nrow_freq,tfr_ncol_time);
  tfr_time_instants = mdr_Create(tfr_ncol_time,1);
  for (i=0; i<tfr_ncol_time;i++)
    MdrV0(tfr_time_instants,i) = i+1;
  tfr_freq_bins = mdr_Create(tfr_nrow_freq,1);  // tfd fills these

  mytfr.N_freq = tfr_nrow_freq;
  mytfr.N_time = tfr_ncol_time;
  mytfr.freq_bins     = MDRPTR(tfr_freq_bins);
  mytfr.time_instants = MDRPTR(tfr_time_instants);
  mytfr.is_complex = 0;
  mytfr.tfr = MDCPTR(tfr);

  ctftbx_grd (mysignal, MDRPTR(winT), wTlen, MDRPTR(winF), wFlen, rs, dsym, mytfr);

  // clean-up
  ent_Clean(e1);
  ent_Clean(e2);
  ent_Clean(e3);
  ent_Clean(e4);
  ent_Clean(e5);

  // cleanup of ctftbx variables:
  // mysignal:
  GC_free (time_instants);
  mysignal.time_instants=0;
  mdc_Destroy(sig);
  // mytfr:
  mytfr.freq_bins     = 0;
  mytfr.time_instants = 0;
  mytfr.tfr = 0;


  rent = ent_Create ();
  bt = btree_Create ();
  ent_SetType (rent, BTREE);
  ent_data (rent) = bt;

  Ent *r1 = ent_Create ();
  ent_SetType (r1, MATRIX_DENSE_REAL);
  ent_data (r1) = tfr_time_instants;
  install (bt, ("time_instants"), r1);

  Ent *r2 = ent_Create ();
  ent_SetType (r2, MATRIX_DENSE_REAL);
  ent_data (r2) = tfr_freq_bins;
  install (bt, ("freq_bins"), r2);

  Ent *r3 = ent_Create ();
  ent_SetType (r3, MATRIX_DENSE_COMPLEX);
  ent_data (r3) = tfr;
  install (bt, ("tfr"), r3);

  return (rent);
}

// Margenau-Hill
Ent *
ent_ctftbx_tfd_mh (int nargs, Datum args[])
{
  Ent *e1=0, *rent;
  MDC *sig=0, *tfr=0;
  MDR *tfr_time_instants=0, *tfr_freq_bins=0;

  type_signal mysignal;
  type_TFR mytfr;

  int i=0;
  int tfr_ncol_time=0, tfr_nrow_freq=0;

  double *time_instants=0;

  Btree *bt=0;

  /* Check n_args */
  if (nargs<1)
    rerror ("tfrd.mh: requires one argument");

  //
  // get the signal
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type(e1)!=MATRIX_DENSE_REAL && ent_type(e1)!=MATRIX_DENSE_COMPLEX)
    rerror ("tfrd.mh: first argument 'signal' has to be real/complex vector");
  if (ent_type(e1)==MATRIX_DENSE_REAL)
  {
    sig = mdr_coerce_mdc (ent_data(e1));
    mysignal.is_complex = 0;
  }
  else
  {
    sig = mdc_Copy(ent_data (e1));
    mysignal.is_complex = 1;
  }
  mysignal.length = sig->nrow * sig->ncol;
  mysignal.sample_freq = 1;
  time_instants = GC_malloc (mysignal.length * sizeof(double));
  for (i=0;i<mysignal.length;i++)
    time_instants[i]= i+1;
  mysignal.time_instants = time_instants;
  mysignal.signal = MDCPTR(sig);

  // create TFR:
  tfr_ncol_time = mysignal.length;
  tfr_nrow_freq = mysignal.length;

  tfr = mdc_Create(tfr_nrow_freq,tfr_ncol_time);
  tfr_time_instants = mdr_Create(tfr_ncol_time,1);
  for (i=0; i<tfr_ncol_time;i++)
    MdrV0(tfr_time_instants,i) = i+1;
  tfr_freq_bins = mdr_Create(tfr_nrow_freq,1);  // tfd fills these

  mytfr.N_freq = tfr_nrow_freq;
  mytfr.N_time = tfr_ncol_time;
  mytfr.freq_bins     = MDRPTR(tfr_freq_bins);
  mytfr.time_instants = MDRPTR(tfr_time_instants);
  mytfr.is_complex = 0;
  mytfr.tfr = MDCPTR(tfr);

  ctftbx_mh (mysignal, mytfr);

  // clean-up
  ent_Clean(e1);


  // cleanup of ctftbx variables:
  // mysignal:
  GC_free (time_instants);
  mysignal.time_instants=0;
  mdc_Destroy(sig);
  // mytfr:
  mytfr.freq_bins     = 0;
  mytfr.time_instants = 0;
  mytfr.tfr = 0;


  rent = ent_Create ();
  bt = btree_Create ();
  ent_SetType (rent, BTREE);
  ent_data (rent) = bt;

  Ent *r1 = ent_Create ();
  ent_SetType (r1, MATRIX_DENSE_REAL);
  ent_data (r1) = tfr_time_instants;
  install (bt, ("time_instants"), r1);

  Ent *r2 = ent_Create ();
  ent_SetType (r2, MATRIX_DENSE_REAL);
  ent_data (r2) = tfr_freq_bins;
  install (bt, ("freq_bins"), r2);

  Ent *r3 = ent_Create ();
  ent_SetType (r3, MATRIX_DENSE_COMPLEX);
  ent_data (r3) = tfr;
  install (bt, ("tfr"), r3);

  return (rent);
}

// page
Ent *
ent_ctftbx_tfd_page (int nargs, Datum args[])
{
  Ent *e1=0, *rent;
  MDC *sig=0, *tfr=0;
  MDR *tfr_time_instants=0, *tfr_freq_bins=0;

  type_signal mysignal;
  type_TFR mytfr;

  int i=0;
  int tfr_ncol_time=0, tfr_nrow_freq=0;

  double *time_instants=0;

  Btree *bt=0;

  /* Check n_args */
  if (nargs<1)
    rerror ("tfrd.page: requires one argument");

  //
  // get the signal
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type(e1)!=MATRIX_DENSE_REAL && ent_type(e1)!=MATRIX_DENSE_COMPLEX)
    rerror ("tfrd.page: first argument 'signal' has to be real/complex vector");
  if (ent_type(e1)==MATRIX_DENSE_REAL)
  {
    sig = mdr_coerce_mdc (ent_data(e1));
    mysignal.is_complex = 0;
  }
  else
  {
    sig = mdc_Copy(ent_data (e1));
    mysignal.is_complex = 1;
  }
  mysignal.length = sig->nrow * sig->ncol;
  mysignal.sample_freq = 1;
  time_instants = GC_malloc (mysignal.length * sizeof(double));
  for (i=0;i<mysignal.length;i++)
    time_instants[i]= i+1;
  mysignal.time_instants = time_instants;
  mysignal.signal = MDCPTR(sig);

  // create TFR:
  tfr_ncol_time = mysignal.length;
  tfr_nrow_freq = mysignal.length;

  tfr = mdc_Create(tfr_nrow_freq,tfr_ncol_time);
  tfr_time_instants = mdr_Create(tfr_ncol_time,1);
  for (i=0; i<tfr_ncol_time;i++)
    MdrV0(tfr_time_instants,i) = i+1;
  tfr_freq_bins = mdr_Create(tfr_nrow_freq,1);  // tfd fills these

  mytfr.N_freq = tfr_nrow_freq;
  mytfr.N_time = tfr_ncol_time;
  mytfr.freq_bins     = MDRPTR(tfr_freq_bins);
  mytfr.time_instants = MDRPTR(tfr_time_instants);
  mytfr.is_complex = 0;
  mytfr.tfr = MDCPTR(tfr);

  ctftbx_page (mysignal, mytfr);

  // clean-up
  ent_Clean(e1);

  // cleanup of ctftbx variables:
  // mysignal:
  GC_free (time_instants);
  mysignal.time_instants=0;
  mdc_Destroy(sig);
  // mytfr:
  mytfr.freq_bins     = 0;
  mytfr.time_instants = 0;
  mytfr.tfr = 0;


  rent = ent_Create ();
  bt = btree_Create ();
  ent_SetType (rent, BTREE);
  ent_data (rent) = bt;

  Ent *r1 = ent_Create ();
  ent_SetType (r1, MATRIX_DENSE_REAL);
  ent_data (r1) = tfr_time_instants;
  install (bt, ("time_instants"), r1);

  Ent *r2 = ent_Create ();
  ent_SetType (r2, MATRIX_DENSE_REAL);
  ent_data (r2) = tfr_freq_bins;
  install (bt, ("freq_bins"), r2);

  Ent *r3 = ent_Create ();
  ent_SetType (r3, MATRIX_DENSE_COMPLEX);
  ent_data (r3) = tfr;
  install (bt, ("tfr"), r3);

  return (rent);
}

// pseudo Margenau-Hill
Ent *
ent_ctftbx_tfd_pmh (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *rent;
  MDR *winT=0;
  MDC *sig=0, *tfr=0;
  MDR *tfr_time_instants=0, *tfr_freq_bins=0;

  type_signal mysignal;
  type_TFR mytfr;

  int wTlen=0, i=0;
  int tfr_ncol_time=0, tfr_nrow_freq=0;

  double *time_instants=0;

  Btree *bt=0;

  /* Check n_args */
  if (nargs<2)
    rerror ("tfrd.pmh: requires two arguments");

  //
  // get the signal
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type(e1)!=MATRIX_DENSE_REAL && ent_type(e1)!=MATRIX_DENSE_COMPLEX)
    rerror ("tfrd.pmh: first argument 'signal' has to be real/complex vector");
  if (ent_type(e1)==MATRIX_DENSE_REAL)
  {
    sig = mdr_coerce_mdc (ent_data(e1));
    mysignal.is_complex = 0;
  }
  else
  {
    sig = mdc_Copy(ent_data (e1));
    mysignal.is_complex = 1;
  }
  mysignal.length = sig->nrow * sig->ncol;
  mysignal.sample_freq = 1;
  time_instants = GC_malloc (mysignal.length * sizeof(double));
  for (i=0;i<mysignal.length;i++)
    time_instants[i]= i+1;
  mysignal.time_instants = time_instants;
  mysignal.signal = MDCPTR(sig);

  //
  // get the windowT length
  //
  e2 = bltin_get_ent (args[1]);
  if (ent_type(e2)!=MATRIX_DENSE_REAL)
    rerror ("tfrd.pmh: second argument 'windowT' has to be real vector");
  winT = class_matrix_real(e2);
  if (winT->nrow !=1 && winT->ncol !=1)
    rerror ("tfrd.pmh: second argument 'windowT' has to be real vector");
  wTlen = winT->nrow * winT->ncol;
  if (wTlen >= mysignal.length)
    rerror ("tfrd.pmh: length of 'windowT' should be less than length of signal");

  // create TFR:
  tfr_ncol_time = mysignal.length;
  tfr_nrow_freq = mysignal.length;

  tfr = mdc_Create(tfr_nrow_freq,tfr_ncol_time);
  tfr_time_instants = mdr_Create(tfr_ncol_time,1);
  for (i=0; i<tfr_ncol_time;i++)
    MdrV0(tfr_time_instants,i) = i+1;
  tfr_freq_bins = mdr_Create(tfr_nrow_freq,1);  // tfd fills these

  mytfr.N_freq = tfr_nrow_freq;
  mytfr.N_time = tfr_ncol_time;
  mytfr.freq_bins     = MDRPTR(tfr_freq_bins);
  mytfr.time_instants = MDRPTR(tfr_time_instants);
  mytfr.is_complex = 0;
  mytfr.tfr = MDCPTR(tfr);

  ctftbx_pmh (mysignal, MDRPTR(winT), wTlen, mytfr);

  // clean-up
  ent_Clean(e1);
  ent_Clean(e2);

  // cleanup of ctftbx variables:
  // mysignal:
  GC_free (time_instants);
  mysignal.time_instants=0;
  mdc_Destroy(sig);
  // mytfr:
  mytfr.freq_bins     = 0;
  mytfr.time_instants = 0;
  mytfr.tfr = 0;


  rent = ent_Create ();
  bt = btree_Create ();
  ent_SetType (rent, BTREE);
  ent_data (rent) = bt;

  Ent *r1 = ent_Create ();
  ent_SetType (r1, MATRIX_DENSE_REAL);
  ent_data (r1) = tfr_time_instants;
  install (bt, ("time_instants"), r1);

  Ent *r2 = ent_Create ();
  ent_SetType (r2, MATRIX_DENSE_REAL);
  ent_data (r2) = tfr_freq_bins;
  install (bt, ("freq_bins"), r2);

  Ent *r3 = ent_Create ();
  ent_SetType (r3, MATRIX_DENSE_COMPLEX);
  ent_data (r3) = tfr;
  install (bt, ("tfr"), r3);

  return (rent);
}

// pseudo page
Ent *
ent_ctftbx_tfd_ppage (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *rent;
  MDR *winT=0;
  MDC *sig=0, *tfr=0;
  MDR *tfr_time_instants=0, *tfr_freq_bins=0;

  type_signal mysignal;
  type_TFR mytfr;

  int wTlen=0, i=0;
  int tfr_ncol_time=0, tfr_nrow_freq=0;

  double *time_instants=0;

  Btree *bt=0;

  /* Check n_args */
  if (nargs<2)
    rerror ("tfrd.ppage: requires two arguments");

  //
  // get the signal
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type(e1)!=MATRIX_DENSE_REAL && ent_type(e1)!=MATRIX_DENSE_COMPLEX)
    rerror ("tfrd.ppage: first argument 'signal' has to be real/complex vector");
  if (ent_type(e1)==MATRIX_DENSE_REAL)
  {
    sig = mdr_coerce_mdc (ent_data(e1));
    mysignal.is_complex = 0;
  }
  else
  {
    sig = mdc_Copy(ent_data (e1));
    mysignal.is_complex = 1;
  }
  mysignal.length = sig->nrow * sig->ncol;
  mysignal.sample_freq = 1;
  time_instants = GC_malloc (mysignal.length * sizeof(double));
  for (i=0;i<mysignal.length;i++)
    time_instants[i]= i+1;
  mysignal.time_instants = time_instants;
  mysignal.signal = MDCPTR(sig);

  //
  // get the windowT length
  //
  e2 = bltin_get_ent (args[1]);
  if (ent_type(e2)!=MATRIX_DENSE_REAL)
    rerror ("tfrd.ppage: second argument 'windowT' has to be real vector");
  winT = class_matrix_real(e2);
  if (winT->nrow !=1 && winT->ncol !=1)
    rerror ("tfrd.ppage: second argument 'windowT' has to be real vector");
  wTlen = winT->nrow * winT->ncol;
  if (wTlen >= mysignal.length)
    rerror ("tfrd.ppage: length of 'windowT' should be less than length of signal");

  // create TFR:
  tfr_ncol_time = mysignal.length;
  tfr_nrow_freq = mysignal.length;

  tfr = mdc_Create(tfr_nrow_freq,tfr_ncol_time);
  tfr_time_instants = mdr_Create(tfr_ncol_time,1);
  for (i=0; i<tfr_ncol_time;i++)
    MdrV0(tfr_time_instants,i) = i+1;
  tfr_freq_bins = mdr_Create(tfr_nrow_freq,1);  // tfd fills these

  mytfr.N_freq = tfr_nrow_freq;
  mytfr.N_time = tfr_ncol_time;
  mytfr.freq_bins     = MDRPTR(tfr_freq_bins);
  mytfr.time_instants = MDRPTR(tfr_time_instants);
  mytfr.is_complex = 0;
  mytfr.tfr = MDCPTR(tfr);

  ctftbx_ppage (mysignal, MDRPTR(winT), wTlen, mytfr);

  // clean-up
  ent_Clean(e1);
  ent_Clean(e2);

  // cleanup of ctftbx variables:
  // mysignal:
  GC_free (time_instants);
  mysignal.time_instants=0;
  mdc_Destroy(sig);
  // mytfr:
  mytfr.freq_bins     = 0;
  mytfr.time_instants = 0;
  mytfr.tfr = 0;


  rent = ent_Create ();
  bt = btree_Create ();
  ent_SetType (rent, BTREE);
  ent_data (rent) = bt;

  Ent *r1 = ent_Create ();
  ent_SetType (r1, MATRIX_DENSE_REAL);
  ent_data (r1) = tfr_time_instants;
  install (bt, ("time_instants"), r1);

  Ent *r2 = ent_Create ();
  ent_SetType (r2, MATRIX_DENSE_REAL);
  ent_data (r2) = tfr_freq_bins;
  install (bt, ("freq_bins"), r2);

  Ent *r3 = ent_Create ();
  ent_SetType (r3, MATRIX_DENSE_COMPLEX);
  ent_data (r3) = tfr;
  install (bt, ("tfr"), r3);

  return (rent);
}

// pseudo Williams-Ville
Ent *
ent_ctftbx_tfd_pwv (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *rent;
  MDR *winT=0;
  MDC *sig=0, *tfr=0;
  MDR *tfr_time_instants=0, *tfr_freq_bins=0;

  type_signal mysignal;
  type_TFR mytfr;

  int wTlen=0, i=0;
  int tfr_ncol_time=0, tfr_nrow_freq=0;

  double *time_instants=0;

  Btree *bt=0;

  /* Check n_args */
  if (nargs<2)
    rerror ("tfrd.pwv: requires two arguments");

  //
  // get the signalwhich requires the time   window;
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type(e1)!=MATRIX_DENSE_REAL && ent_type(e1)!=MATRIX_DENSE_COMPLEX)
    rerror ("tfrd.pwv: first argument 'signal' has to be real/complex vector");
  if (ent_type(e1)==MATRIX_DENSE_REAL)
  {
    sig = mdr_coerce_mdc (ent_data(e1));
    mysignal.is_complex = 0;
  }
  else
  {
    sig = mdc_Copy(ent_data (e1));
    mysignal.is_complex = 1;
  }
  mysignal.length = sig->nrow * sig->ncol;
  mysignal.sample_freq = 1;
  time_instants = GC_malloc (mysignal.length * sizeof(double));
  for (i=0;i<mysignal.length;i++)
    time_instants[i]= i+1;
  mysignal.time_instants = time_instants;
  mysignal.signal = MDCPTR(sig);

  //
  // get the windowT length
  //
  e2 = bltin_get_ent (args[1]);
  if (ent_type(e2)!=MATRIX_DENSE_REAL)
    rerror ("tfrd.pwv: second argument 'windowT' has to be real vector");
  winT = class_matrix_real(e2);
  if (winT->nrow !=1 && winT->ncol !=1)
    rerror ("tfrd.pwv: second argument 'windowT' has to be real vector");
  wTlen = winT->nrow * winT->ncol;
  if (wTlen >= mysignal.length)
    rerror ("tfrd.pwv: length of 'windowT' should be less than length of signal");

  // create TFR:
  tfr_ncol_time = mysignal.length;
  tfr_nrow_freq = mysignal.length;

  tfr = mdc_Create(tfr_nrow_freq,tfr_ncol_time);
  tfr_time_instants = mdr_Create(tfr_ncol_time,1);
  for (i=0; i<tfr_ncol_time;i++)
    MdrV0(tfr_time_instants,i) = i+1;
  tfr_freq_bins = mdr_Create(tfr_nrow_freq,1);  // tfd fills these

  mytfr.N_freq = tfr_nrow_freq;
  mytfr.N_time = tfr_ncol_time;
  mytfr.freq_bins     = MDRPTR(tfr_freq_bins);
  mytfr.time_instants = MDRPTR(tfr_time_instants);
  mytfr.is_complex = 0;
  mytfr.tfr = MDCPTR(tfr);

  ctftbx_pwv (mysignal, MDRPTR(winT), wTlen, mytfr);

  // clean-up
  ent_Clean(e1);
  ent_Clean(e2);

  // cleanup of ctftbx variables:
  // mysignal:
  GC_free (time_instants);
  mysignal.time_instants=0;
  mdc_Destroy(sig);
  // mytfr:
  mytfr.freq_bins     = 0;
  mytfr.time_instants = 0;
  mytfr.tfr = 0;


  rent = ent_Create ();
  bt = btree_Create ();
  ent_SetType (rent, BTREE);
  ent_data (rent) = bt;

  Ent *r1 = ent_Create ();
  ent_SetType (r1, MATRIX_DENSE_REAL);
  ent_data (r1) = tfr_time_instants;
  install (bt, ("time_instants"), r1);

  Ent *r2 = ent_Create ();
  ent_SetType (r2, MATRIX_DENSE_REAL);
  ent_data (r2) = tfr_freq_bins;
  install (bt, ("freq_bins"), r2);

  Ent *r3 = ent_Create ();
  ent_SetType (r3, MATRIX_DENSE_COMPLEX);
  ent_data (r3) = tfr;
  install (bt, ("tfr"), r3);

  return (rent);
}

// ri
Ent *
ent_ctftbx_tfd_ri (int nargs, Datum args[])
{
  Ent *e1=0, *rent;
  MDC *sig=0, *tfr=0;
  MDR *tfr_time_instants=0, *tfr_freq_bins=0;

  type_signal mysignal;
  type_TFR mytfr;

  int i=0;
  int tfr_ncol_time=0, tfr_nrow_freq=0;

  double *time_instants=0;

  Btree *bt=0;

  /* Check n_args */
  if (nargs<1)
    rerror ("tfrd.ri: requires one argument");

  //
  // get the signal
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type(e1)!=MATRIX_DENSE_REAL && ent_type(e1)!=MATRIX_DENSE_COMPLEX)
    rerror ("tfrd.ri: first argument 'signal' has to be real/complex vector");
  if (ent_type(e1)==MATRIX_DENSE_REAL)
  {
    sig = mdr_coerce_mdc (ent_data(e1));
    mysignal.is_complex = 0;
  }
  else
  {
    sig = mdc_Copy(ent_data (e1));
    mysignal.is_complex = 1;
  }
  mysignal.length = sig->nrow * sig->ncol;
  mysignal.sample_freq = 1;
  time_instants = GC_malloc (mysignal.length * sizeof(double));
  for (i=0;i<mysignal.length;i++)
    time_instants[i]= i+1;
  mysignal.time_instants = time_instants;
  mysignal.signal = MDCPTR(sig);

  // create TFR:
  tfr_ncol_time = mysignal.length;
  tfr_nrow_freq = mysignal.length;

  tfr = mdc_Create(tfr_nrow_freq,tfr_ncol_time);
  tfr_time_instants = mdr_Create(tfr_ncol_time,1);
  for (i=0; i<tfr_ncol_time;i++)
    MdrV0(tfr_time_instants,i) = i+1;
  tfr_freq_bins = mdr_Create(tfr_nrow_freq,1);  // tfd fills these

  mytfr.N_freq = tfr_nrow_freq;
  mytfr.N_time = tfr_ncol_time;
  mytfr.freq_bins     = MDRPTR(tfr_freq_bins);
  mytfr.time_instants = MDRPTR(tfr_time_instants);
  mytfr.is_complex = 1;
  mytfr.tfr = MDCPTR(tfr);

  ctftbx_ri (mysignal, mytfr);

  // clean-up
  ent_Clean(e1);

  // cleanup of ctftbx variables:
  // mysignal:
  GC_free (time_instants);
  mysignal.time_instants=0;
  mdc_Destroy(sig);
  // mytfr:
  mytfr.freq_bins     = 0;
  mytfr.time_instants = 0;
  mytfr.tfr = 0;


  rent = ent_Create ();
  bt = btree_Create ();
  ent_SetType (rent, BTREE);
  ent_data (rent) = bt;

  Ent *r1 = ent_Create ();
  ent_SetType (r1, MATRIX_DENSE_REAL);
  ent_data (r1) = tfr_time_instants;
  install (bt, ("time_instants"), r1);

  Ent *r2 = ent_Create ();
  ent_SetType (r2, MATRIX_DENSE_REAL);
  ent_data (r2) = tfr_freq_bins;
  install (bt, ("freq_bins"), r2);

  Ent *r3 = ent_Create ();
  ent_SetType (r3, MATRIX_DENSE_COMPLEX);
  ent_data (r3) = tfr;
  install (bt, ("tfr"), r3);

  return (rent);
}

// Smoothed Pseudo-Wigner-Ville Distribution
Ent *
ent_ctftbx_tfd_spwv (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *e3=0, *rent;
  MDR *winT=0, *winF=0;
  MDC *sig=0, *tfr=0;
  MDR *tfr_time_instants=0, *tfr_freq_bins=0;

  type_signal mysignal;
  type_TFR mytfr;

  int wTlen=0, wFlen=0, i=0;
  int tfr_ncol_time=0, tfr_nrow_freq=0;

  double *time_instants=0;

  Btree *bt=0;

  /* Check n_args */
  if (nargs<3)
    rerror ("tfrd.spwv: requires three arguments");

  //
  // get the signal
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type(e1)!=MATRIX_DENSE_REAL && ent_type(e1)!=MATRIX_DENSE_COMPLEX)
    rerror ("tfrd.spwv: first argument 'signal' has to be real/complex vector");
  if (ent_type(e1)==MATRIX_DENSE_REAL)
  {
    sig = mdr_coerce_mdc (ent_data(e1));
    mysignal.is_complex = 0;
  }
  else
  {
    sig = mdc_Copy(ent_data (e1));
    mysignal.is_complex = 1;
  }
  mysignal.length = sig->nrow * sig->ncol;
  mysignal.sample_freq = 1;
  time_instants = GC_malloc (mysignal.length * sizeof(double));
  for (i=0;i<mysignal.length;i++)
    time_instants[i]= i+1;
  mysignal.time_instants = time_instants;
  mysignal.signal = MDCPTR(sig);

  //
  // get the windowT length
  //
  e2 = bltin_get_ent (args[1]);
  if (ent_type(e2)!=MATRIX_DENSE_REAL)
    rerror ("tfrd.spwv: second argument 'windowT' has to be real vector");
  winT = class_matrix_real(e2);
  if (winT->nrow !=1 && winT->ncol !=1)
    rerror ("tfrd.spwv: second argument 'windowT' has to be real vector");
  wTlen = winT->nrow * winT->ncol;
  if (wTlen >= mysignal.length)
    rerror ("tfrd.spwv: length of 'windowT' should be less than length of signal");

  //
  // get the windowT length
  //
  e3 = bltin_get_ent (args[2]);
  if (ent_type(e3)!=MATRIX_DENSE_REAL)
    rerror ("tfrd.spwv: third argument 'windowF' has to be real vector");
  winF = class_matrix_real(e3);
  if (winF->nrow !=1 && winF->ncol !=1)
    rerror ("tfrd.spwv: second argument 'windowF' has to be real vector");
  wFlen = winF->nrow * winF->ncol;
  if (wFlen >= mysignal.length)
    rerror ("tfrd.spwv: length of 'windowF' should be less than length of signal");

  // create TFR:
  tfr_ncol_time = mysignal.length;
  tfr_nrow_freq = mysignal.length;

  tfr = mdc_Create(tfr_nrow_freq,tfr_ncol_time);
  tfr_time_instants = mdr_Create(tfr_ncol_time,1);
  for (i=0; i<tfr_ncol_time;i++)
    MdrV0(tfr_time_instants,i) = i+1;
  tfr_freq_bins = mdr_Create(tfr_nrow_freq,1);  // tfd fills these

  mytfr.N_freq = tfr_nrow_freq;
  mytfr.N_time = tfr_ncol_time;
  mytfr.freq_bins     = MDRPTR(tfr_freq_bins);
  mytfr.time_instants = MDRPTR(tfr_time_instants);
  mytfr.is_complex = 0;
  mytfr.tfr = MDCPTR(tfr);

  ctftbx_spwv (mysignal, MDRPTR(winT), wTlen, MDRPTR(winF), wFlen, mytfr);

  // clean-up
  ent_Clean(e1);
  ent_Clean(e2);
  ent_Clean(e3);

  // cleanup of ctftbx variables:
  // mysignal:
  GC_free (time_instants);
  mysignal.time_instants=0;
  mdc_Destroy(sig);
  // mytfr:
  mytfr.freq_bins     = 0;
  mytfr.time_instants = 0;
  mytfr.tfr = 0;


  rent = ent_Create ();
  bt = btree_Create ();
  ent_SetType (rent, BTREE);
  ent_data (rent) = bt;

  Ent *r1 = ent_Create ();
  ent_SetType (r1, MATRIX_DENSE_REAL);
  ent_data (r1) = tfr_time_instants;
  install (bt, ("time_instants"), r1);

  Ent *r2 = ent_Create ();
  ent_SetType (r2, MATRIX_DENSE_REAL);
  ent_data (r2) = tfr_freq_bins;
  install (bt, ("freq_bins"), r2);

  Ent *r3 = ent_Create ();
  ent_SetType (r3, MATRIX_DENSE_COMPLEX);
  ent_data (r3) = tfr;
  install (bt, ("tfr"), r3);

  return (rent);
}

// Wigner-Ville
Ent *
ent_ctftbx_tfd_wv (int nargs, Datum args[])
{
  Ent *e1=0, *rent;
  MDC *sig=0, *tfr=0;
  MDR *tfr_time_instants=0, *tfr_freq_bins=0;

  type_signal mysignal;
  type_TFR mytfr;

  int i=0;
  int tfr_ncol_time=0, tfr_nrow_freq=0;

  double *time_instants=0;

  Btree *bt=0;

  /* Check n_args */
  if (nargs<1)
    rerror ("tfrd.wv: requires one argument");

  //
  // get the signal
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type(e1)!=MATRIX_DENSE_REAL && ent_type(e1)!=MATRIX_DENSE_COMPLEX)
    rerror ("tfrd.wv: first argument 'signal' has to be real/complex vector");
  if (ent_type(e1)==MATRIX_DENSE_REAL)
  {
    sig = mdr_coerce_mdc (ent_data(e1));
    mysignal.is_complex = 0;
  }
  else
  {
    sig = mdc_Copy(ent_data (e1));
    mysignal.is_complex = 1;
  }
  mysignal.length = sig->nrow * sig->ncol;
  mysignal.sample_freq = 1;
  time_instants = GC_malloc (mysignal.length * sizeof(double));
  for (i=0;i<mysignal.length;i++)
    time_instants[i]= i+1;
  mysignal.time_instants = time_instants;
  mysignal.signal = MDCPTR(sig);

  // create TFR:
  tfr_ncol_time = mysignal.length;
  tfr_nrow_freq = mysignal.length;

  tfr = mdc_Create(tfr_nrow_freq,tfr_ncol_time);
  tfr_time_instants = mdr_Create(tfr_ncol_time,1);
  for (i=0; i<tfr_ncol_time;i++)
    MdrV0(tfr_time_instants,i) = i+1;
  tfr_freq_bins = mdr_Create(tfr_nrow_freq,1);  // tfd fills these

  mytfr.N_freq = tfr_nrow_freq;
  mytfr.N_time = tfr_ncol_time;
  mytfr.freq_bins     = MDRPTR(tfr_freq_bins);
  mytfr.time_instants = MDRPTR(tfr_time_instants);
  mytfr.is_complex = 0;
  mytfr.tfr = MDCPTR(tfr);

  ctftbx_wv (mysignal, mytfr);

  // clean-up
  ent_Clean(e1);

  // cleanup of ctftbx variables:
  // mysignal:
  GC_free (time_instants);
  mysignal.time_instants=0;
  mdc_Destroy(sig);
  // mytfr:
  mytfr.freq_bins     = 0;
  mytfr.time_instants = 0;
  mytfr.tfr = 0;


  rent = ent_Create ();
  bt = btree_Create ();
  ent_SetType (rent, BTREE);
  ent_data (rent) = bt;

  Ent *r1 = ent_Create ();
  ent_SetType (r1, MATRIX_DENSE_REAL);
  ent_data (r1) = tfr_time_instants;
  install (bt, ("time_instants"), r1);

  Ent *r2 = ent_Create ();
  ent_SetType (r2, MATRIX_DENSE_REAL);
  ent_data (r2) = tfr_freq_bins;
  install (bt, ("freq_bins"), r2);

  Ent *r3 = ent_Create ();
  ent_SetType (r3, MATRIX_DENSE_COMPLEX);
  ent_data (r3) = tfr;
  install (bt, ("tfr"), r3);

  return (rent);
}

// zam
Ent *
ent_ctftbx_tfd_zam (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *e3=0, *rent;
  MDR *winT=0, *winF=0;
  MDC *sig=0, *tfr=0;
  MDR *tfr_time_instants=0, *tfr_freq_bins=0;

  type_signal mysignal;
  type_TFR mytfr;

  int wTlen=0, wFlen=0, i=0;
  int tfr_ncol_time=0, tfr_nrow_freq=0;

  double *time_instants=0;

  Btree *bt=0;

  /* Check n_args */
  if (nargs<3)
    rerror ("tfrd.zam: requires three arguments");

  //
  // get the signal
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type(e1)!=MATRIX_DENSE_REAL && ent_type(e1)!=MATRIX_DENSE_COMPLEX)
    rerror ("tfrd.zam: first argument 'signal' has to be real/complex vector");
  if (ent_type(e1)==MATRIX_DENSE_REAL)
  {
    sig = mdr_coerce_mdc (ent_data(e1));
    mysignal.is_complex = 0;
  }
  else
  {
    sig = mdc_Copy(ent_data (e1));
    mysignal.is_complex = 1;
  }
  mysignal.length = sig->nrow * sig->ncol;
  mysignal.sample_freq = 1;
  time_instants = GC_malloc (mysignal.length * sizeof(double));
  for (i=0;i<mysignal.length;i++)
    time_instants[i]= i+1;
  mysignal.time_instants = time_instants;
  mysignal.signal = MDCPTR(sig);

  //
  // get the windowT length
  //
  e2 = bltin_get_ent (args[1]);
  if (ent_type(e2)!=MATRIX_DENSE_REAL)
    rerror ("tfrd.zam: second argument 'windowT' has to be real vector");
  winT = class_matrix_real(e2);
  if (winT->nrow !=1 && winT->ncol !=1)
    rerror ("tfrd.zam: second argument 'windowT' has to be real vector");
  wTlen = winT->nrow * winT->ncol;
  if (wTlen >= mysignal.length)
    rerror ("tfrd.zam: length of 'windowT' should be less than length of signal");

  //
  // get the windowT length
  //
  e3 = bltin_get_ent (args[2]);
  if (ent_type(e3)!=MATRIX_DENSE_REAL)
    rerror ("tfrd.zam: third argument 'windowF' has to be real vector");
  winF = class_matrix_real(e3);
  if (winF->nrow !=1 && winF->ncol !=1)
    rerror ("tfrd.zam: second argument 'windowF' has to be real vector");
  wFlen = winF->nrow * winF->ncol;
  if (wFlen >= mysignal.length)
    rerror ("tfrd.zam: length of 'windowF' should be less than length of signal");

  // create TFR:
  tfr_ncol_time = mysignal.length;
  tfr_nrow_freq = mysignal.length;

  tfr = mdc_Create(tfr_nrow_freq,tfr_ncol_time);
  tfr_time_instants = mdr_Create(tfr_ncol_time,1);
  for (i=0; i<tfr_ncol_time;i++)
    MdrV0(tfr_time_instants,i) = i+1;
  tfr_freq_bins = mdr_Create(tfr_nrow_freq,1);  // tfd fills these

  mytfr.N_freq = tfr_nrow_freq;
  mytfr.N_time = tfr_ncol_time;
  mytfr.freq_bins     = MDRPTR(tfr_freq_bins);
  mytfr.time_instants = MDRPTR(tfr_time_instants);
  mytfr.is_complex = 0;
  mytfr.tfr = MDCPTR(tfr);

  ctftbx_zam (mysignal, MDRPTR(winT), wTlen, MDRPTR(winF), wFlen, mytfr);

  // clean-up
  ent_Clean(e1);
  ent_Clean(e2);
  ent_Clean(e3);

  // cleanup of ctftbx variables:
  // mysignal:
  GC_free (time_instants);
  mysignal.time_instants=0;
  mdc_Destroy(sig);
  // mytfr:
  mytfr.freq_bins     = 0;
  mytfr.time_instants = 0;
  mytfr.tfr = 0;


  rent = ent_Create ();
  bt = btree_Create ();
  ent_SetType (rent, BTREE);
  ent_data (rent) = bt;

  Ent *r1 = ent_Create ();
  ent_SetType (r1, MATRIX_DENSE_REAL);
  ent_data (r1) = tfr_time_instants;
  install (bt, ("time_instants"), r1);

  Ent *r2 = ent_Create ();
  ent_SetType (r2, MATRIX_DENSE_REAL);
  ent_data (r2) = tfr_freq_bins;
  install (bt, ("freq_bins"), r2);

  Ent *r3 = ent_Create ();
  ent_SetType (r3, MATRIX_DENSE_COMPLEX);
  ent_data (r3) = tfr;
  install (bt, ("tfr"), r3);

  return (rent);
}

//
//
// time frequency reduced interference distributions
//
//

// bessel
Ent *
ent_ctftbx_tfrid_bessel (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *e3=0, *rent;
  MDR *winT=0, *winF=0;
  MDC *sig=0, *tfr=0;
  MDR *tfr_time_instants=0, *tfr_freq_bins=0;

  type_signal mysignal;
  type_TFR mytfr;

  int wTlen=0, wFlen=0, i=0;
  int tfr_ncol_time=0, tfr_nrow_freq=0;

  double *time_instants=0;

  Btree *bt=0;

  /* Check n_args */
  if (nargs<3)
    rerror ("tfrid.bessel: requires three arguments");

  //
  // get the signal
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type(e1)!=MATRIX_DENSE_REAL && ent_type(e1)!=MATRIX_DENSE_COMPLEX)
    rerror ("tfrid.bessel: first argument 'signal' has to be real/complex vector");
  if (ent_type(e1)==MATRIX_DENSE_REAL)
  {
    sig = mdr_coerce_mdc (ent_data(e1));
    mysignal.is_complex = 0;
  }
  else
  {
    sig = mdc_Copy(ent_data (e1));
    mysignal.is_complex = 1;
  }
  mysignal.length = sig->nrow * sig->ncol;
  mysignal.sample_freq = 1;
  time_instants = GC_malloc (mysignal.length * sizeof(double));
  for (i=0;i<mysignal.length;i++)
    time_instants[i]= i+1;
  mysignal.time_instants = time_instants;
  mysignal.signal = MDCPTR(sig);

  //
  // get the windowT length
  //
  e2 = bltin_get_ent (args[1]);
  if (ent_type(e2)!=MATRIX_DENSE_REAL)
    rerror ("tfrid.bessel: second argument 'windowT' has to be real vector");
  winT = class_matrix_real(e2);
  if (winT->nrow !=1 && winT->ncol !=1)
    rerror ("tfrid.bessel: second argument 'windowT' has to be real vector");
  wTlen = winT->nrow * winT->ncol;
  if (wTlen >= mysignal.length)
    rerror ("tfrid.bessel: length of 'windowT' should be less than length of signal");

  //
  // get the windowT length
  //
  e3 = bltin_get_ent (args[2]);
  if (ent_type(e3)!=MATRIX_DENSE_REAL)
    rerror ("tfrid.bessel: third argument 'windowF' has to be real vector");
  winF = class_matrix_real(e3);
  if (winF->nrow !=1 && winF->ncol !=1)
    rerror ("tfrid.bessel: third argument 'windowF' has to be real vector");
  wFlen = winF->nrow * winF->ncol;
  if (wFlen >= mysignal.length)
    rerror ("tfrid.bessel: length of 'windowF' should be less than length of signal");

  // create TFR:
  tfr_ncol_time = mysignal.length;
  tfr_nrow_freq = mysignal.length;

  tfr = mdc_Create(tfr_nrow_freq,tfr_ncol_time);
  tfr_time_instants = mdr_Create(tfr_ncol_time,1);
  for (i=0; i<tfr_ncol_time;i++)
    MdrV0(tfr_time_instants,i) = i+1;
  tfr_freq_bins = mdr_Create(tfr_nrow_freq,1);  // tfd fills these

  mytfr.N_freq = tfr_nrow_freq;
  mytfr.N_time = tfr_ncol_time;
  mytfr.freq_bins     = MDRPTR(tfr_freq_bins);
  mytfr.time_instants = MDRPTR(tfr_time_instants);
  mytfr.is_complex = 1;
  mytfr.tfr = MDCPTR(tfr);

  ctftbx_ridb (mysignal, MDRPTR(winT), wTlen, MDRPTR(winF), wFlen, mytfr);

  // clean-up
  ent_Clean(e1);
  ent_Clean(e2);
  ent_Clean(e3);

  // cleanup of ctftbx variables:
  // mysignal:
  GC_free (time_instants);
  mysignal.time_instants=0;
  mdc_Destroy(sig);
  // mytfr:
  mytfr.freq_bins     = 0;
  mytfr.time_instants = 0;
  mytfr.tfr = 0;


  rent = ent_Create ();
  bt = btree_Create ();
  ent_SetType (rent, BTREE);
  ent_data (rent) = bt;

  Ent *r1 = ent_Create ();
  ent_SetType (r1, MATRIX_DENSE_REAL);
  ent_data (r1) = tfr_time_instants;
  install (bt, ("time_instants"), r1);

  Ent *r2 = ent_Create ();
  ent_SetType (r2, MATRIX_DENSE_REAL);
  ent_data (r2) = tfr_freq_bins;
  install (bt, ("freq_bins"), r2);

  Ent *r3 = ent_Create ();
  ent_SetType (r3, MATRIX_DENSE_COMPLEX);
  ent_data (r3) = tfr;
  install (bt, ("tfr"), r3);

  return (rent);
}

// binomial
Ent *
ent_ctftbx_tfrid_binom (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *e3=0, *rent;
  MDR *winT=0, *winF=0;
  MDC *sig=0, *tfr=0;
  MDR *tfr_time_instants=0, *tfr_freq_bins=0;

  type_signal mysignal;
  type_TFR mytfr;

  int wTlen=0, wFlen=0, i=0;
  int tfr_ncol_time=0, tfr_nrow_freq=0;

  double *time_instants=0;

  Btree *bt=0;

  /* Check n_args */
  if (nargs<3)
    rerror ("tfrid.binom: requires three arguments");

  //
  // get the signal
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type(e1)!=MATRIX_DENSE_REAL && ent_type(e1)!=MATRIX_DENSE_COMPLEX)
    rerror ("tfrid.binom: first argument 'signal' has to be real/complex vector");
  if (ent_type(e1)==MATRIX_DENSE_REAL)
  {
    sig = mdr_coerce_mdc (ent_data(e1));
    mysignal.is_complex = 0;
  }
  else
  {
    sig = mdc_Copy(ent_data (e1));
    mysignal.is_complex = 1;
  }
  mysignal.length = sig->nrow * sig->ncol;
  mysignal.sample_freq = 1;
  time_instants = GC_malloc (mysignal.length * sizeof(double));
  for (i=0;i<mysignal.length;i++)
    time_instants[i]= i+1;
  mysignal.time_instants = time_instants;
  mysignal.signal = MDCPTR(sig);

  //
  // get the windowT length
  //
  e2 = bltin_get_ent (args[1]);
  if (ent_type(e2)!=MATRIX_DENSE_REAL)
    rerror ("tfrid.binom: second argument 'windowT' has to be real vector");
  winT = class_matrix_real(e2);
  if (winT->nrow !=1 && winT->ncol !=1)
    rerror ("tfrid.binom: second argument 'windowT' has to be real vector");
  wTlen = winT->nrow * winT->ncol;
  if (wTlen >= mysignal.length)
    rerror ("tfrid.binom: length of 'windowT' should be less than length of signal");

  //
  // get the windowT length
  //
  e3 = bltin_get_ent (args[2]);
  if (ent_type(e3)!=MATRIX_DENSE_REAL)
    rerror ("tfrid.binom: third argument 'windowF' has to be real vector");
  winF = class_matrix_real(e3);
  if (winF->nrow !=1 && winF->ncol !=1)
    rerror ("tfrid.binom: third argument 'windowF' has to be real vector");
  wFlen = winF->nrow * winF->ncol;
  if (wFlen >= mysignal.length)
    rerror ("tfrid.binom: length of 'windowF' should be less than length of signal");

  // create TFR:
  tfr_ncol_time = mysignal.length;
  tfr_nrow_freq = mysignal.length;

  tfr = mdc_Create(tfr_nrow_freq,tfr_ncol_time);
  tfr_time_instants = mdr_Create(tfr_ncol_time,1);
  for (i=0; i<tfr_ncol_time;i++)
    MdrV0(tfr_time_instants,i) = i+1;
  tfr_freq_bins = mdr_Create(tfr_nrow_freq,1);  // tfd fills these

  mytfr.N_freq = tfr_nrow_freq;
  mytfr.N_time = tfr_ncol_time;
  mytfr.freq_bins     = MDRPTR(tfr_freq_bins);
  mytfr.time_instants = MDRPTR(tfr_time_instants);
  mytfr.is_complex = 1;
  mytfr.tfr = MDCPTR(tfr);

  ctftbx_ridbn (mysignal, MDRPTR(winT), wTlen, MDRPTR(winF), wFlen, mytfr);

  // clean-up
  ent_Clean(e1);
  ent_Clean(e2);
  ent_Clean(e3);

  // cleanup of ctftbx variables:
  // mysignal:
  GC_free (time_instants);
  mysignal.time_instants=0;
  mdc_Destroy(sig);
  // mytfr:
  mytfr.freq_bins     = 0;
  mytfr.time_instants = 0;
  mytfr.tfr = 0;


  rent = ent_Create ();
  bt = btree_Create ();
  ent_SetType (rent, BTREE);
  ent_data (rent) = bt;

  Ent *r1 = ent_Create ();
  ent_SetType (r1, MATRIX_DENSE_REAL);
  ent_data (r1) = tfr_time_instants;
  install (bt, ("time_instants"), r1);

  Ent *r2 = ent_Create ();
  ent_SetType (r2, MATRIX_DENSE_REAL);
  ent_data (r2) = tfr_freq_bins;
  install (bt, ("freq_bins"), r2);

  Ent *r3 = ent_Create ();
  ent_SetType (r3, MATRIX_DENSE_COMPLEX);
  ent_data (r3) = tfr;
  install (bt, ("tfr"), r3);

  return (rent);
}

// trigonal
Ent *
ent_ctftbx_tfrid_tri (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *e3=0, *rent;
  MDR *winT=0, *winF=0;
  MDC *sig=0, *tfr=0;
  MDR *tfr_time_instants=0, *tfr_freq_bins=0;

  type_signal mysignal;
  type_TFR mytfr;

  int wTlen=0, wFlen=0, i=0;
  int tfr_ncol_time=0, tfr_nrow_freq=0;

  double *time_instants=0;

  Btree *bt=0;

  /* Check n_args */
  if (nargs<3)
    rerror ("tfrid.tri: requires three arguments");

  //
  // get the signal
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type(e1)!=MATRIX_DENSE_REAL && ent_type(e1)!=MATRIX_DENSE_COMPLEX)
    rerror ("tfrid.tri: first argument 'signal' has to be real/complex vector");
  if (ent_type(e1)==MATRIX_DENSE_REAL)
  {
    sig = mdr_coerce_mdc (ent_data(e1));
    mysignal.is_complex = 0;
  }
  else
  {
    sig = mdc_Copy(ent_data (e1));
    mysignal.is_complex = 1;
  }
  mysignal.length = sig->nrow * sig->ncol;
  mysignal.sample_freq = 1;
  time_instants = GC_malloc (mysignal.length * sizeof(double));
  for (i=0;i<mysignal.length;i++)
    time_instants[i]= i+1;
  mysignal.time_instants = time_instants;
  mysignal.signal = MDCPTR(sig);

  //
  // get the windowT length
  //
  e2 = bltin_get_ent (args[1]);
  if (ent_type(e2)!=MATRIX_DENSE_REAL)
    rerror ("tfrid.tri: second argument 'windowT' has to be real vector");
  winT = class_matrix_real(e2);
  if (winT->nrow !=1 && winT->ncol !=1)
    rerror ("tfrid.tri: second argument 'windowT' has to be real vector");
  wTlen = winT->nrow * winT->ncol;
  if (wTlen >= mysignal.length)
    rerror ("tfrid.tri: length of 'windowT' should be less than length of signal");

  //
  // get the windowT length
  //
  e3 = bltin_get_ent (args[2]);
  if (ent_type(e3)!=MATRIX_DENSE_REAL)
    rerror ("tfrid.tri: third argument 'windowF' has to be real vector");
  winF = class_matrix_real(e3);
  if (winF->nrow !=1 && winF->ncol !=1)
    rerror ("tfrid.tri: third argument 'windowF' has to be real vector");
  wFlen = winF->nrow * winF->ncol;
  if (wFlen >= mysignal.length)
    rerror ("tfrid.tri: length of 'windowF' should be less than length of signal");

  // create TFR:
  tfr_ncol_time = mysignal.length;
  tfr_nrow_freq = mysignal.length;

  tfr = mdc_Create(tfr_nrow_freq,tfr_ncol_time);
  tfr_time_instants = mdr_Create(tfr_ncol_time,1);
  for (i=0; i<tfr_ncol_time;i++)
    MdrV0(tfr_time_instants,i) = i+1;
  tfr_freq_bins = mdr_Create(tfr_nrow_freq,1);  // tfd fills these

  mytfr.N_freq = tfr_nrow_freq;
  mytfr.N_time = tfr_ncol_time;
  mytfr.freq_bins     = MDRPTR(tfr_freq_bins);
  mytfr.time_instants = MDRPTR(tfr_time_instants);
  mytfr.is_complex = 1;
  mytfr.tfr = MDCPTR(tfr);

  ctftbx_ridt (mysignal, MDRPTR(winT), wTlen, MDRPTR(winF), wFlen, mytfr);

  // clean-up
  // clean-up
  ent_Clean(e1);
  ent_Clean(e2);
  ent_Clean(e3);

  // cleanup of ctftbx variables:
  // mysignal:
  GC_free (time_instants);
  mysignal.time_instants=0;
  mdc_Destroy(sig);
  // mytfr:
  mytfr.freq_bins     = 0;
  mytfr.time_instants = 0;
  mytfr.tfr = 0;


  rent = ent_Create ();
  bt = btree_Create ();
  ent_SetType (rent, BTREE);
  ent_data (rent) = bt;

  Ent *r1 = ent_Create ();
  ent_SetType (r1, MATRIX_DENSE_REAL);
  ent_data (r1) = tfr_time_instants;
  install (bt, ("time_instants"), r1);

  Ent *r2 = ent_Create ();
  ent_SetType (r2, MATRIX_DENSE_REAL);
  ent_data (r2) = tfr_freq_bins;
  install (bt, ("freq_bins"), r2);

  Ent *r3 = ent_Create ();
  ent_SetType (r3, MATRIX_DENSE_COMPLEX);
  ent_data (r3) = tfr;
  install (bt, ("tfr"), r3);

  return (rent);
}

Ent *
ent_ctftbx_tfrid_hann (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *e3=0, *rent;
  MDR *winT=0, *winF=0;
  MDC *sig=0, *tfr=0;
  MDR *tfr_time_instants=0, *tfr_freq_bins=0;

  type_signal mysignal;
  type_TFR mytfr;

  int wTlen=0, wFlen=0, i=0;
  int tfr_ncol_time=0, tfr_nrow_freq=0;

  double *time_instants=0;

  Btree *bt=0;

  /* Check n_args */
  if (nargs<3)
    rerror ("tfrid.hann: requires three arguments");

  //
  // get the signal
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type(e1)!=MATRIX_DENSE_REAL && ent_type(e1)!=MATRIX_DENSE_COMPLEX)
    rerror ("tfrid.hann: first argument 'signal' has to be real/complex vector");
  if (ent_type(e1)==MATRIX_DENSE_REAL)
  {
    sig = mdr_coerce_mdc (ent_data(e1));
    mysignal.is_complex = 0;
  }
  else
  {
    sig = mdc_Copy(ent_data (e1));
    mysignal.is_complex = 1;
  }
  mysignal.length = sig->nrow * sig->ncol;
  mysignal.sample_freq = 1;
  time_instants = GC_malloc (mysignal.length * sizeof(double));
  for (i=0;i<mysignal.length;i++)
    time_instants[i]= i+1;
  mysignal.time_instants = time_instants;
  mysignal.signal = MDCPTR(sig);

  //
  // get the windowT length
  //
  e2 = bltin_get_ent (args[1]);
  if (ent_type(e2)!=MATRIX_DENSE_REAL)
    rerror ("tfrid.hann: second argument 'windowT' has to be real vector");
  winT = class_matrix_real(e2);
  if (winT->nrow !=1 && winT->ncol !=1)
    rerror ("tfrid.hann: second argument 'windowT' has to be real vector");
  wTlen = winT->nrow * winT->ncol;
  if (wTlen >= mysignal.length)
    rerror ("tfrid.hann: length of 'windowT' should be less than length of signal");

  //
  // get the windowT length
  //
  e3 = bltin_get_ent (args[2]);
  if (ent_type(e3)!=MATRIX_DENSE_REAL)
    rerror ("tfrid.hann: third argument 'windowF' has to be real vector");
  winF = class_matrix_real(e3);
  if (winF->nrow !=1 && winF->ncol !=1)
    rerror ("tfrid.hann: third argument 'windowF' has to be real vector");
  wFlen = winF->nrow * winF->ncol;
  if (wFlen >= mysignal.length)
    rerror ("tfrid.hann: length of 'windowF' should be less than length of signal");

  // create TFR:
  tfr_ncol_time = mysignal.length;
  tfr_nrow_freq = mysignal.length;

  tfr = mdc_Create(tfr_nrow_freq,tfr_ncol_time);
  tfr_time_instants = mdr_Create(tfr_ncol_time,1);
  for (i=0; i<tfr_ncol_time;i++)
    MdrV0(tfr_time_instants,i) = i+1;
  tfr_freq_bins = mdr_Create(tfr_nrow_freq,1);  // tfd fills these

  mytfr.N_freq = tfr_nrow_freq;
  mytfr.N_time = tfr_ncol_time;
  mytfr.freq_bins     = MDRPTR(tfr_freq_bins);
  mytfr.time_instants = MDRPTR(tfr_time_instants);
  mytfr.is_complex = 1;
  mytfr.tfr = MDCPTR(tfr);

  ctftbx_ridh (mysignal, MDRPTR(winT), wTlen, MDRPTR(winF), wFlen, mytfr);

  // clean-up
  ent_Clean(e1);
  ent_Clean(e2);
  ent_Clean(e3);

  // cleanup of ctftbx variables:
  // mysignal:
  GC_free (time_instants);
  mysignal.time_instants=0;
  mdc_Destroy(sig);
  // mytfr:
  mytfr.freq_bins     = 0;
  mytfr.time_instants = 0;
  mytfr.tfr = 0;


  rent = ent_Create ();
  bt = btree_Create ();
  ent_SetType (rent, BTREE);
  ent_data (rent) = bt;

  Ent *r1 = ent_Create ();
  ent_SetType (r1, MATRIX_DENSE_REAL);
  ent_data (r1) = tfr_time_instants;
  install (bt, ("time_instants"), r1);

  Ent *r2 = ent_Create ();
  ent_SetType (r2, MATRIX_DENSE_REAL);
  ent_data (r2) = tfr_freq_bins;
  install (bt, ("freq_bins"), r2);

  Ent *r3 = ent_Create ();
  ent_SetType (r3, MATRIX_DENSE_COMPLEX);
  ent_data (r3) = tfr;
  install (bt, ("tfr"), r3);

  return (rent);
}


//
// ambiguity function
//
Ent *
ent_ctftbx_af (int nargs, Datum args[])
{
  Ent *e1=0, *rent;
  MDC *sig=0, *af=0;

  double *time_instants=0;
  type_signal mysignal;

  int i;

  type_AF     myaf;
  MDR *af_delay_bins=0, *af_doppler_bins=0;
  int N_doppler=0, N_delay=0;

  Btree *bt=0;

  /* Check n_args */
  if (nargs<1)
    rerror ("tframbfun: requires single argument");

  //
  // get the signal
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type(e1)!=MATRIX_DENSE_REAL && ent_type(e1)!=MATRIX_DENSE_COMPLEX)
    rerror ("tframbfun: first argument 'signal' has to be real/complex vector");
  if (ent_type(e1)==MATRIX_DENSE_REAL)
  {
    sig = mdr_coerce_mdc (ent_data(e1));
    mysignal.is_complex = 0;
  }
  else
  {
    sig = mdc_Copy(ent_data (e1));
    mysignal.is_complex = 1;
  }
  mysignal.length = sig->nrow * sig->ncol;
  mysignal.sample_freq = 1;
  time_instants = GC_malloc (mysignal.length * sizeof(double));
  for (i=0;i<mysignal.length;i++)
    time_instants[i]= i+1;
  mysignal.time_instants = time_instants;
  mysignal.signal = MDCPTR(sig);

  // create AF:
  N_delay   = mysignal.length;
  N_doppler = N_delay;
  af = mdc_Create(N_doppler, N_delay);
  mdc_Zero(af);
  af_doppler_bins = mdr_Create(N_doppler,1);
  af_delay_bins   = mdr_Create(N_delay,1);
  for (i=0; i<N_delay; i++)
    MdrV0(af_delay_bins,i) = i - (N_delay-(N_delay&1))/2;

  myaf.N_doppler = N_doppler;
  myaf.N_delay   = N_delay;
  myaf.doppler_bins = MDRPTR(af_doppler_bins);
  myaf.delay_bins   = MDRPTR(af_delay_bins);
  myaf.is_complex = 1;
  myaf.af = MDCPTR(af);

  ctftbx_af (mysignal, myaf);

  // clean-up
  ent_Clean(e1);

  // cleanup of ctftbx variables:
  // mysignal:
  GC_free (time_instants);
  mysignal.time_instants=0;
  mdc_Destroy(sig);
  // myaf:
  myaf.doppler_bins = 0;
  myaf.delay_bins   = 0;
  myaf.af = 0;

  rent = ent_Create ();
  bt = btree_Create ();
  ent_SetType (rent, BTREE);
  ent_data (rent) = bt;

  Ent *r1 = ent_Create ();
  ent_SetType (r1, MATRIX_DENSE_REAL);
  ent_data (r1) = af_doppler_bins;
  install (bt, ("doppler_bins"), r1);

  Ent *r2 = ent_Create ();
  ent_SetType (r2, MATRIX_DENSE_REAL);
  ent_data (r2) = af_delay_bins;
  install (bt, ("delay_bins"), r2);

  Ent *r3 = ent_Create ();
  ent_SetType (r3, MATRIX_DENSE_COMPLEX);
  ent_data (r3) = af;
  install (bt, ("af"), r3);

  return (rent);
}

//
// kernel
//
Ent *
ent_ctftbx_kernel (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *e3=0, *e4=0, *rent;
  MDR *px=0;
  MDC *af=0;

  char *kernel_name=0;

  double *params=0;

  int i, N_param=0;

  type_AF     myaf;
  MDR *af_delay_bins=0, *af_doppler_bins=0;
  int N_doppler=0, N_delay=0;

  Btree *bt=0;

  /* Check n_args */
  if (nargs<3)
    rerror ("tfrkernel: requires at least three arguments");

  //
  // N_doppler
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type(e1)!=MATRIX_DENSE_REAL)
    rerror ("tfrkernel: first argument 'ndoppler' has to be positive vector");
  N_doppler = (int) class_double(e1);


  //
  // N_delay
  //
  e2 = bltin_get_ent (args[1]);
  if (ent_type(e2)!=MATRIX_DENSE_REAL)
    rerror ("tfrkernel: second argument 'ndelay' has to be positive vector");
  N_delay = (int) class_double(e2);

  //
  // kernel name
  //
  e3 = bltin_get_ent (args[2]);
  if (ent_type(e3)!=MATRIX_DENSE_STRING)
    rerror ("tfrkernel: third argument 'kernel_name' has to be string");
  kernel_name = class_char_pointer (e3);
  if (!kernel_name)
    rerror ("tfrkernel: third argument 'kernel_name' has to be string");

  //
  // parameters for the kernel
  //
  if (nargs > 3 && kernel_name[0]!='w' && kernel_name[0]!='W')
  {
    e4 = bltin_get_ent (args[3]);
    if (ent_type(e4)==MATRIX_DENSE_REAL)
    {
      px = class_matrix_real(e4);
      if (px->nrow !=1 && px->ncol !=1)
        rerror ("tfrkernel: fourth argument 'kernel_params' has to be real vector");
      N_param = px->nrow * px->ncol;
      params  = MDRPTR(px);
    }
    else
      rerror ("tfrkernel: fourth argument 'kernel_params' has to be real vector");
    if ((kernel_name[0]=='m' || kernel_name[0]=='M') && (N_param!=7))
    {
      printf ("tfrkernel: 'mtek' type kernel requires 7 parameters");
      rerror ("tfrkernel: improper fourth argument 'kernel_params'");
    }
    else if( (kernel_name[0]=='r' || kernel_name[0]=='R') && ((N_param&1)==0))
    {
      printf ("tfrkernel: 'rgk' type kernel requires odd-number of parameters");
      rerror ("tfrkernel: improper fourth argument 'kernel_params'");
    }
    else if((kernel_name[0]=='g' || kernel_name[0]=='G')&&(N_param < 2))
    {
      printf ("tfrkernel: 'gmcw' type kernel requires at least two parameters");
      rerror ("tfrkernel: improper fourth argument 'kernel_params'");
    }
    else if((kernel_name[0]=='s' || kernel_name[0]=='S')&&((N_param&1)==1))
    {
      printf ("tfrkernel: 'spectro' type kernel requires single parameter");
      rerror ("tfrkernel: improper fourth argument 'kernel_params'");
    }
  }
  else
  {
    N_param = 0;
    params  = 0;
  }

  // create AF:
  af = mdc_Create(N_doppler, N_delay);
  mdc_Zero(af);
  af_doppler_bins = mdr_Create(N_doppler,1);
  for (i=0; i<N_doppler; i++)
    MdrV0(af_doppler_bins,i) = i - (int)(N_doppler/2.0);
  af_delay_bins   = mdr_Create(N_delay,1);
  for (i=0; i<N_delay; i++)
    MdrV0(af_delay_bins,i) = i - (int)(N_delay/2.0);

  myaf.N_doppler = N_doppler;
  myaf.N_delay   = N_delay;
  myaf.doppler_bins = MDRPTR(af_doppler_bins);
  myaf.delay_bins   = MDRPTR(af_delay_bins);
  myaf.is_complex = 1;
  myaf.af = MDCPTR(af);

  ctftbx_kernel (kernel_name, params, N_param, myaf);

  // clean-up
  ent_Clean(e1);
  ent_Clean(e2);
  ent_Clean(e3);
  ent_Clean(e4);

  // cleanup of ctftbx variables:
  // myaf:
  myaf.doppler_bins = 0;
  myaf.delay_bins   = 0;
  myaf.af = 0;

  rent = ent_Create ();
  bt = btree_Create ();
  ent_SetType (rent, BTREE);
  ent_data (rent) = bt;

  Ent *r1 = ent_Create ();
  ent_SetType (r1, MATRIX_DENSE_REAL);
  ent_data (r1) = af_doppler_bins;
  install (bt, ("doppler_bins"), r1);

  Ent *r2 = ent_Create ();
  ent_SetType (r2, MATRIX_DENSE_REAL);
  ent_data (r2) = af_delay_bins;
  install (bt, ("delay_bins"), r2);

  Ent *r3 = ent_Create ();
  ent_SetType (r3, MATRIX_DENSE_COMPLEX);
  ent_data (r3) = af;
  install (bt, ("af"), r3);

  return (rent);
}


//
// ambiguity function
//
Ent *
ent_ctftbx_af2tfr (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *rent;

  MDR *tfr_freq_bins=0, *tfr_time_instants=0;
  MDC *tfr_tfr=0;

  ListNode *node_af_doppler_bins=0, *node_af_delay_bins=0, *node_af_af=0;
  ListNode *node_ker_doppler_bins=0, *node_ker_delay_bins=0, *node_ker_af=0;

  int i;

  type_AF  myaf, myker;
  type_TFR mytfr;

  MDR *af_delay_bins=0, *af_doppler_bins=0;
  MDC *af_af=0;
  MDR *ker_delay_bins=0, *ker_doppler_bins=0;
  MDC *ker_af=0;

  Btree *bt=0;

  if (nargs<2)
    rerror ("af2tfr: requires two arguments");

  //
  // first list: af
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type(e1)!=BTREE)
    rerror ("af2tfr: first argument 'af' has to be AF-type list");
  node_af_doppler_bins = btree_FindNode (ent_data (e1), "doppler_bins");
  node_af_delay_bins   = btree_FindNode (ent_data (e1), "delay_bins");
  node_af_af           = btree_FindNode (ent_data (e1), "af");
  if (!node_af_doppler_bins || !node_af_delay_bins || !node_af_af)
    rerror ("af2tfr: first argument 'af' has to be AF-type list");
  if (ent_type(var_ent (node_af_doppler_bins))==MATRIX_DENSE_REAL)
    af_doppler_bins = class_matrix_real (var_ent (node_af_doppler_bins));
  if (ent_type(var_ent (node_af_delay_bins))==MATRIX_DENSE_REAL)
    af_delay_bins = class_matrix_real (var_ent (node_af_delay_bins));
  if (ent_type(var_ent (node_af_af))==MATRIX_DENSE_COMPLEX)
    af_af = ent_data (var_ent (node_af_af));
  if (!af_doppler_bins || !af_delay_bins || !af_af)
    rerror ("af2tfr: first argument 'af' has to be AF-type list");

  //
  // second list: ker
  //
  e2 = bltin_get_ent (args[1]);
  if (ent_type(e2)!=BTREE)
    rerror ("af2tfr: second argument 'ker' has to be AF-type list");
  node_ker_doppler_bins = btree_FindNode (ent_data (e1), "doppler_bins");
  node_ker_delay_bins   = btree_FindNode (ent_data (e1), "delay_bins");
  node_ker_af           = btree_FindNode (ent_data (e1), "af");
  if (!node_ker_doppler_bins || !node_ker_delay_bins || !node_ker_af)
    rerror ("af2tfr: first argument 'af' has to be AF-type list");
  if (ent_type(var_ent (node_ker_doppler_bins))==MATRIX_DENSE_REAL)
    ker_doppler_bins = class_matrix_real (var_ent (node_ker_doppler_bins));
  if (ent_type(var_ent (node_ker_delay_bins))==MATRIX_DENSE_REAL)
    ker_delay_bins = class_matrix_real (var_ent (node_ker_delay_bins));
  if (ent_type(var_ent (node_ker_af))==MATRIX_DENSE_COMPLEX)
    ker_af = ent_data (var_ent (node_ker_af));
  if (!ker_doppler_bins || !ker_delay_bins || !ker_af)
    rerror ("af2tfr: first argument 'af' has to be AF-type list");

  // assemble af:
  myaf.N_doppler = af_doppler_bins->nrow * af_doppler_bins->ncol;
  myaf.N_delay   = af_delay_bins->nrow * af_delay_bins->ncol;
  myaf.doppler_bins = MDRPTR(af_doppler_bins);
  myaf.delay_bins   = MDRPTR(af_delay_bins);
  myaf.is_complex = 1;
  myaf.af = MDCPTR(af_af);

  // assemble ker: force it to be real
  myker.N_doppler = ker_doppler_bins->nrow * ker_doppler_bins->ncol;
  myker.N_delay   = ker_delay_bins->nrow * ker_delay_bins->ncol;
  myker.doppler_bins = MDRPTR(ker_doppler_bins);
  myker.delay_bins   = MDRPTR(ker_delay_bins);
  myker.is_complex = 0;
  myker.af = MDCPTR(ker_af);

  if (myaf.N_doppler != myker.N_doppler ||
      myaf.N_delay != myker.N_delay)
    rerror ("af2tfr: 'af' and 'ker' must have same size");

  // create tfr
  mytfr.N_freq = myker.N_doppler;
  mytfr.N_time = myker.N_delay;
  tfr_freq_bins     = mdr_Create(myker.N_doppler,1);
  tfr_time_instants = mdr_Create(mytfr.N_time,1);
  tfr_tfr = mdc_Create(myker.N_doppler,mytfr.N_time);
  mdc_Zero(tfr_tfr);
  for (i=0; i<mytfr.N_time;i++)
    MdrV0(tfr_time_instants,i) = i+1;
  for (i=0; i<myker.N_doppler;i++)
    MdrV0(tfr_freq_bins,i) = ((double) i)/((double) myker.N_doppler);
  mytfr.freq_bins     = MDRPTR(tfr_freq_bins);
  mytfr.time_instants = MDRPTR(tfr_time_instants);
  mytfr.is_complex = 0;
  mytfr.tfr = MDCPTR(tfr_tfr);

  ctftbx_af2tfr (myaf, myker,mytfr);

  // clean-up
  ent_Clean(e1);
  ent_Clean(e2);

  // cleanup of ctftbx variables:
  // myaf:
  myaf.doppler_bins = 0;
  myaf.delay_bins   = 0;
  myaf.af = 0;
  // myker:
  myker.doppler_bins = 0;
  myker.delay_bins   = 0;
  myker.af = 0;
  // mytfr:
  mytfr.freq_bins     = 0;
  mytfr.time_instants = 0;
  mytfr.tfr = 0;

  rent = ent_Create ();
  bt = btree_Create ();
  ent_SetType (rent, BTREE);
  ent_data (rent) = bt;

  Ent *r1 = ent_Create ();
  ent_SetType (r1, MATRIX_DENSE_REAL);
  ent_data (r1) = tfr_freq_bins;
  install (bt, ("freq_bins"), r1);

  Ent *r2 = ent_Create ();
  ent_SetType (r2, MATRIX_DENSE_REAL);
  ent_data (r2) = tfr_time_instants;
  install (bt, ("time_instants"), r2);

  Ent *r3 = ent_Create ();
  ent_SetType (r3, MATRIX_DENSE_COMPLEX);
  ent_data (r3) = tfr_tfr;
  install (bt, ("tfr"), r3);

  return (rent);
}

//
//
// Tony Fisher's mkfilter in Jim Peters interpretation
//
//
#undef  THIS_SOLVER
#define THIS_SOLVER "mkfilter"
#include "clibs/fidlib/fidlib.h"
Ent * ent_mkfilter (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *e3=0, *rent=0;
  Btree *bt=0;

  int i, find_delay=0;
  double rate=-1, delay=0;
  char *err=0, *filter_spec=0, *opts=0;

  MDR *fir_b=0, *iir_a=0;

  FILE *rlab_stderr = (!RLAB_STDERR_DS) ? stderr : RLAB_STDERR_DS;

  /* Check n_args */
  if ((nargs!=1) &&  (nargs!=2) &&  (nargs!=3))
  {
    fprintf (rlab_stderr, THIS_SOLVER ": Digital filter creation tool.\n");
    fprintf (rlab_stderr, THIS_SOLVER ": Format:\n");
    fprintf (rlab_stderr, THIS_SOLVER ":  y = "THIS_SOLVER"(f_s, desc/,opts/),\n");
    fprintf (rlab_stderr, THIS_SOLVER ": where:\n");
    fprintf (rlab_stderr, THIS_SOLVER ": 'rate' is a sampling frequency in Hz, and");
    fprintf (rlab_stderr, " 'desc' is a filter description, see below, while 'opts' is a list of options.\n");
    fid_list_filters(rlab_stderr);
    rerror (THIS_SOLVER ": " RLAB_ERROR_ONE_OR_TWO_ARG_REQUIRED);
  }

  //
  // get the sampling rate
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type(e1)!=MATRIX_DENSE_REAL)
    rerror (THIS_SOLVER ": " RLAB_ERROR_ARG1_MDR_SCALAR);
  rate = class_double(e1);
  if (rate <= 0)
  {
    fprintf (rlab_stderr, " 'rate' must be positive scalar!\n");
    rerror (THIS_SOLVER ": " RLAB_ERROR_ARG1_MDR_SCALAR);
  }

  //
  // get the filter description
  //
  e2 = bltin_get_ent (args[1]);
  if (ent_type(e2)!=MATRIX_DENSE_STRING)
    rerror (THIS_SOLVER ": " RLAB_ERROR_ARG2_MDS_SCALAR);
  filter_spec = class_char_pointer(e2);
  if (!isvalidstring(filter_spec))
    rerror (THIS_SOLVER ": " RLAB_ERROR_ARG2_MDS_SCALAR);

  //
  // options:
  //  string containing options 'delay,....'
  //
  if (nargs > 2)
  {
    e3 = bltin_get_ent (args[2]);
    if (ent_type (e3) == MATRIX_DENSE_STRING)
    {
      opts = class_char_pointer(e3);
      if (isvalidstring(opts))
      {
        if (strstr(opts, "delay")!=0)
        {
          find_delay = 1;
        }
      }
    }
  }



  // Design a filter, and optionally get its long description
  FidFilter *ffl=0, *ff=0;
  err = fid_parse(rate, (char **) &filter_spec, (FidFilter **) &ffl);
  FidFilter *rff=ffl;

  int have_gain=0, have_a=0, have_b=0;

  //
  // return a list of filter parameters
  //
  rent = ent_Create ();
  bt = btree_Create ();
  ent_SetType (rent, BTREE);
  ent_data (rent) = bt;

  if (!err)
  {
    ff = fid_flatten(ffl);
    FidFilter *dff=ff;

    if (find_delay)
    {
      delay += fid_calc_delay(ff);
    }

    for (; ff->typ ; ff=FFNEXT(ff) )
    {

      if (ff->typ == 'F')
      {
        if (ff->len == 1)
        {
          if (!have_gain)
          {
            have_gain = 1;
            install (bt, "gain", ent_Create_Rlab_Double(ff->val[0]));
          }
          continue;
        }
        if (!have_b)
        {
          fir_b = mdr_Create(1, ff->len);
          for (i=0; i<ff->len; i++)
          {
            MdrV0(fir_b,i) = ff->val[i];
          }
          have_b = 1;
          install (bt, "b", ent_Assign_Rlab_MDR(fir_b));
          continue;
        }
      }

      if (ff->typ == 'I')
      {
        if (!have_a)
        {
          iir_a = mdr_Create(1, ff->len);
          for (i=0; i<ff->len; i++)
          {
            MdrV0(iir_a,i) = ff->val[i];
          }
          have_a = 1;
          install (bt, "a", ent_Assign_Rlab_MDR(iir_a));
        }
      }
    }

    free(dff);
  }
  else
  {
    fprintf (rlab_stderr, THIS_SOLVER ": %s\n", err);
  }

  free(err);
  free(rff);

  if ((!have_gain)&&(fir_b))
  {
    install (bt, "gain", ent_Create_Rlab_Double(MdrV0(fir_b,0)));
    for (i=1; i<SIZE(fir_b); i++)
    {
      MdrV0(fir_b,i) /= MdrV0(fir_b,0);
    }
    MdrV0(fir_b,0) = 1;
  }

  if (find_delay)
  {
    install (bt, "delay", ent_Create_Rlab_Double(delay));
  }

  // clean-up
  ent_Clean(e1);
  ent_Clean(e2);
  ent_Clean(e3);

  return (rent);
}


