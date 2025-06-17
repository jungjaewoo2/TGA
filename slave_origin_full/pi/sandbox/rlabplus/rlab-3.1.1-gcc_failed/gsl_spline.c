// Copyright (C) 2000-2004 Marijan Kostrun
//   part of rlab+4linux project on rlabplus.sourceforge.net
//
// cubic spline library
//
//   This program is free software; you can redistribute it and/or modify
//   it under the terms of the GNU General Public License as published by
//   the Free Software Foundation; either version 2 of the License, or
//   (at your option) any later version.
//
//   This program is distributed in the hope that it will be useful,
//   but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//   GNU General Public License for more details.
//
//   You should have received a copy of the GNU General Public License
//   along with this program; if not, write to the Free Software
//   Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
//
//   See the file ./COPYING
//   ********************************************************************** */

// rlab headers
#include "complex.h"
#include "rlab.h"
#include "ent.h"
#include "class.h"
#include "mem.h"
#include "bltin.h"
#include "util.h"
#include "mdr.h"
#include "mdc.h"
#include "mdcf1.h"
#include "mdcf2.h"
#include "mathl.h"

// gsl headers
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

#include "rlab_solver_parameters_names.h"

#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <math.h>

// find the matrix of [x_i,y(x_i),y''(x_ii)]
// cubicspline_init(x,y,yp1,ypN)
MDR *
mdr_cubicspline_init (MDR * x, MDR * y, double yp1, double ypN)
{
  int n, i;
  MDR *w=0, *a=0, *d=0, *lm=0, *mu=0;
  double dm, h, hp, hm;
  n = MNR (x);
  a = mdr_Create (n, 1);
  d = mdr_Create (n, 1);
  lm = mdr_Create (n, 1);
  mu = mdr_Create (n, 1);

  w = mdr_Create (n, 3);

  /* initialize the arrays */
  Mdr1 (d, 1, 1) = 0;
  Mdr1 (lm, 1, 1) = 0;
  if (yp1 < 0.9e30)
  {
    h = Mdr1 (x, 2, 1) - Mdr1 (x, 1, 1);
    Mdr1 (d, 1, 1) = 6. / h * ((Mdr1 (y, 2, 1) - Mdr1 (y, 1, 1)) / h - yp1);
    Mdr1 (lm, 1, 1) = 1.;
  }
  for (i = 2; i <= n - 1; i++)
  {
    hp = Mdr1 (x, i + 1, 1) - Mdr1 (x, i, 1);
    hm = Mdr1 (x, i, 1) - Mdr1 (x, i - 1, 1);
    Mdr1 (d, i, 1) = 6. / (hp + hm) * ((Mdr1 (y, i + 1, 1) -
                                        Mdr1 (y, i, 1)) / hp - (Mdr1 (y, i,
                                                                      1) -
                                                                Mdr1 (y, i - 1,
                                                                      1)) / hm);
    Mdr1 (lm, i, 1) =
      (Mdr1 (x, i + 1, 1) - Mdr1 (x, i, 1)) / (Mdr1 (x, i + 1, 1) -
                                               Mdr1 (x, i - 1, 1));
    Mdr1 (mu, i, 1) = 1. - Mdr1 (lm, i, 1);
  }
  Mdr1 (mu, n, 1) = 0.;
  Mdr1 (d, n, 1) = 0.;
  if (ypN < 0.9e30)
  {
    h = Mdr1 (x, n, 1) - Mdr1 (x, n - 1, 1);
    Mdr1 (mu, n, 1) = 1.;
    Mdr1 (d, n, 1) = 6. / h * (ypN - (Mdr1 (y, n, 1) - Mdr1 (y, n - 1, 1)) / h);
  }
  /* first run */
  Mdr1 (a, 1, 1) = 2.;
  for (i = 2; i <= n; i++)
  {
    dm = Mdr1 (mu, i, 1) / Mdr1 (a, i - 1, 1);
    Mdr1 (a, i, 1) = 2. - dm * Mdr1 (lm, i - 1, 1);
    Mdr1 (d, i, 1) -= dm * Mdr1 (d, i - 1, 1);
  }
  Mdr1 (w, n, 3) = Mdr1 (d, n, 1) / Mdr1 (a, n, 1);
  Mdr1 (w, n, 2) = Mdr1 (y, n, 1);
  Mdr1 (w, n, 1) = Mdr1 (x, n, 1);
  for (i = n - 1; i >= 1; i--)
  {
    Mdr1 (w, i, 3) =
      (Mdr1 (d, i, 1) - Mdr1 (lm, i, 1) * Mdr1 (w, i + 1, 3)) / Mdr1 (a, i, 1);
    Mdr1 (w, i, 2) = Mdr1 (y, i, 1);
    Mdr1 (w, i, 1) = Mdr1 (x, i, 1);
  }

  mdr_Destroy(a);
  mdr_Destroy(d);
  mdr_Destroy(lm);
  mdr_Destroy(mu);

  return w;
}

Ent *
ent_cubicspline_init (int nargs, Datum args[])
{
  Ent *X=0, *Y=0, *YP1=0, *YPN=0;
  MDR *x=0, *y=0, *w;
  double yp1 = 1.e31, ypN = 2.e31;
  int xr, yr, i;

  if (nargs < 2)
  {
    fprintf (stdout,
             "csplinterp: Finds matrix of second derivatives for cubic spline fit\n");
    fprintf (stdout,
             "csplinterp: for natural (yp1,ypN not given), periodic (yp1=ypN && y1=yN)\n");
    fprintf (stdout,
             "csplinterp: and standard boundary conditions (yp1,ypN given).\n");
    fprintf (stdout, "csplinterp: Format:\n");
    fprintf (stdout, "csplinterp:   m = csplinterp(y,x /,yp1,ypN/),\n");
    fprintf (stdout,
             "csplinterp: where 'y' are the fitted values, 'x' is independent variable,\n");
    fprintf (stdout,
             "csplinterp: 'yp1' is the value of the first derivative at x1 (left end),\n");
    fprintf (stdout,
             "csplinterp: while 'ypN' is the first derivative at xN (right end).\n");
    rerror ("requires at least two arguments");
  }
  //
  // get Y
  //
  Y = bltin_get_ent (args[0]);
  if (ent_type (Y) != MATRIX_DENSE_REAL)
    rerror ("csplinterp: 'y' has to be a single-column real matrix !");
  y = ent_data (Y);
  yr = MNR (y);
  if (MNC (y) != 1)
    rerror ("csplinterp: 'x' has to be a single-column real matrix !");
  //
  // get X
  //
  X = bltin_get_ent (args[1]);
  if (ent_type (X) != MATRIX_DENSE_REAL)
    rerror ("csplinterp: 'y' has to be a single-column real matrix !");
  x = ent_data (X);
  if (MNC (x) != 1)
    rerror ("csplinterp: 'x' has to be a single-column real matrix !");
  xr = MNR (x);
  if (MNR (x) != MNR (y))
    rerror ("csplinterp: 'y' and 'x' should have the same size !");
  // is y'(x_1) given? if not use default value which sets y''(x_1)=0
  if (nargs > 2)
  {
    YP1 = bltin_get_ent (args[2]);
    if (ent_type (YP1) == MATRIX_DENSE_REAL)
      yp1 = class_double (YP1);
  }
  // is y'(x_N) given? if not use default value which sets y''(x_N)=0
  if (nargs > 3)
  {
    YPN = bltin_get_ent (args[3]);
    if (ent_type (YPN) == MATRIX_DENSE_REAL)
      ypN = class_double (YPN);
  }


  // test for periodic boundary conditions
  if (yp1 == ypN && Mdr1 (y, 1, 1) == Mdr1 (y, yr, 1))
  {
    // we have a periodic boundary conditions,
    // we do not really need yp1 and yp2
    w = mdr_Create (xr, 3);
    for (i = 1; i <= xr; i++)
    {
      Mdr1 (w, i, 1) = Mdr1 (x, i, 1);
      Mdr1 (w, i, 2) = Mdr1 (y, i, 1);
    }
    gsl_interp_accel *acc = gsl_interp_accel_alloc ();

    gsl_spline *spline = gsl_spline_alloc (gsl_interp_cspline_periodic, xr);

    gsl_spline_init (spline, MDRPTR(x), MDRPTR(y), xr);

    for (i = 1; i <= xr; i++)
      Mdr1 (w, i, 3) = gsl_spline_eval_deriv2 (spline, Mdr1 (x, i, 1), acc);

    gsl_spline_free (spline);

    gsl_interp_accel_free (acc);
  }
  else
  {
    // boundary conditions: non-periodic or natural
    w = mdr_cubicspline_init (x, y, yp1, ypN);
  }

  // clean stuff
  ent_Clean (X);
  ent_Clean (Y);
  ent_Clean (YP1);
  ent_Clean (YPN);

  return ent_Assign_Rlab_MDR(w);
}

/* cs */
MDR *
mdr_cubicspline_eval (MDR * x, MDR * m)
{
  // m = [ x0, y0, y0'' ]
  int xr, i, j, j0, nm;
  MDR *w=0;
  double a=0, b=0, g=0, d=0, h=0, dx=0;
  xr = MNR (x) * MNC(x);
  nm = MNR (m);
  w = mdr_Create (xr, 1);
  mdr_Zero (w);
  j = 2;
  j0 = 0;
  for (i = 1; i <= xr; i++)
  {
    if (mdrV1 (x, i) <= mdr1 (m, j - 1, 1))
    {
      // test if the entry is below, presumable first entry in the
      // lookup table m = [x,y(x),y"(x)]
      MdrV1 (w, i) = mdr1 (m, j - 1, 2);
      continue;
    }
    while (mdr1 (m, j, 1) < mdrV1 (x, i) && j < nm)
      j++;                      // x[j-1] < x <= x[j]
    if (j == nm && mdrV1 (x, i) >= mdr1 (m, nm, 1))
    {
      // x > x[nm] use y[nm]
      MdrV1 (w, i) = mdr1 (m, nm, 2);
      continue;
    }
    if (mdr1 (m, j, 1) == mdrV1 (x, i))
    {
      // x=x[j]
      MdrV1 (w, i) = mdr1 (m, j, 2);
      continue;
    }
    if (j0 != j)
    {
      // j0!=j calculate the factors
      h = mdr1 (m, j, 1) - mdr1 (m, j - 1, 1);
      a = mdr1 (m, j - 1, 2);
      b = (mdr1 (m, j, 2) - mdr1 (m, j - 1, 2)) / h - (2 * mdr1 (m, j - 1, 3)
                                                       + mdr1 (m, j,
                                                               3)) / 6 * h;
      g = mdr1 (m, j - 1, 3) / 2;
      d = (mdr1 (m, j, 3) - mdr1 (m, j - 1, 3)) / (6 * h);
      j0 = j;
    }
    dx = mdrV1 (x, i) - mdr1 (m, j - 1, 1);
    MdrV1 (w, i) = a + dx * (b + dx * (g + d * dx));
  }
  return w;
}

Ent *
ent_cubicspline_eval (int nargs, Datum args[])
{
  Ent *X=0, *M=0, *rent;
  MDR *x, *m;
  if (nargs < 2)
  {
    printf ("cspleval: Cubic-spline interpolation.\n");
    printf ("cspleval: Format:\n");
    printf ("cspleval:   y = cspleval(x,m),\n");
    printf ("cspleval: where 'x' are the interpolant values, while 'm' is\n");
    printf
      ("cspleval: the matrix of second derivatives found by 'csplinterp'.\n");
    rerror ("requires two arguments");
  }
  //
  // get X
  //
  X = bltin_get_ent (args[0]);
  if (ent_type (X) != MATRIX_DENSE_REAL)
    rerror ("cspleval: 'x' has to be a single-column real matrix !");
  x = ent_data (X);
  if (MNC (x)!=1 && MNR(x)!=1 )
    rerror ("cspleval: 'x' has to be a single-column matrix!\n");
  //
  // get M
  //
  M = bltin_get_ent (args[1]);
  if (ent_type (M) != MATRIX_DENSE_REAL)
    rerror("cspleval: 'm' has to be a three-column matrix [x,y(x),y''(x)]!\n");
  m = ent_data (M);
  if (MNC (m) != 3)
    rerror("cspleval: 'm' has to be a three-column matrix [x,y(x),y''(x)]!\n");

  rent = ent_Assign_Rlab_MDR (mdr_cubicspline_eval (x, m) );

  ent_Clean (X);
  ent_Clean (M);

  return rent;
}




// see gsl-n.nn/interpolation/cspline.c
typedef struct
{
  double *c;
  double *g;
  double *diag;
  double *offdiag;
} cspline_state_t;

Ent *
ent_cubicspline_init_limited (int nargs, Datum args[])
{
  Ent *X=0, *Y=0, *YP1=0, *YPN=0;
  MDR *x, *y, *w;
  double yp1 = 1.e31, ypN = 2.e31;
  int xr, yr, i;

  if (nargs < 2)
  {
    fprintf (stdout,
             "csplinterp: Finds matrix of second derivatives for cubic spline fit\n");
    fprintf (stdout,
             "csplinterp: for natural (yp1,ypN not given), periodic (yp1=ypN && y1=yN)\n");
    fprintf (stdout,
             "csplinterp: and standard boundary conditions (yp1,ypN given).\n");
    fprintf (stdout, "csplinterp: Format:\n");
    fprintf (stdout, "csplinterp:   m = csplinterp(y,x /,yp1,ypN/),\n");
    fprintf (stdout,
             "csplinterp: where 'y' are the fitted values, 'x' is independent variable,\n");
    fprintf (stdout,
             "csplinterp: 'yp1' is the value of the first derivative at x1 (left end),\n");
    fprintf (stdout,
             "csplinterp: while 'ypN' is the first derivative at xN (right end).\n");
    rerror ("requires at least two arguments");
  }
  //
  // get Y
  //
  Y = bltin_get_ent (args[0]);
  if (ent_type (Y) != MATRIX_DENSE_REAL)
    rerror ("csplinterp: 'y' has to be a single-column real matrix !");
  y = ent_data (Y);
  yr = y->nrow * y->ncol;
  if (y->nrow != 1 && y->ncol != 1)
    rerror ("csplinterp: 'y' has to be a single-column real matrix !");
  //
  // get X
  //
  X = bltin_get_ent (args[1]);
  if (ent_type (X) != MATRIX_DENSE_REAL)
    rerror ("csplinterp: 'x' has to be a single-column real matrix !");
  x = ent_data (X);
  xr = x->nrow * x->ncol;
  if (x->nrow != 1 && x->ncol != 1)
    rerror ("csplinterp: 'x' has to be a single-column real matrix !");
  if (xr != yr)
    rerror ("csplinterp: 'x' and 'y' should have the same size !");

  // is y'(x_1) given? if not use default value which sets y''(x_1)=0
  if (nargs > 2)
  {
    YP1 = bltin_get_ent (args[2]);
    if (ent_type (YP1) == MATRIX_DENSE_REAL)
      yp1 = class_double (YP1);
  }

  // is y'(x_N) given? if not use default value which sets y''(x_N)=0
  if (nargs > 3)
  {
    YPN = bltin_get_ent (args[3]);
    if (ent_type (YPN) == MATRIX_DENSE_REAL)
      ypN = class_double (YPN);
  }

  // test for periodic boundary conditions
  gsl_interp_accel *acc = gsl_interp_accel_alloc ();

  gsl_interp *spline=0;
  if (yp1 == ypN && MdrV0(y,0) == MdrV0(y, yr-1))
  {
    // we have a periodic boundary conditions,
    // we do not really need yp1 and yp2
    spline = gsl_interp_alloc (gsl_interp_cspline_periodic, xr);
    gsl_interp_init (spline, MDRPTR(x), MDRPTR(y), xr);
  }
  else
  {
    // we have a periodic boundary conditions,
    // we do not really need yp1 and yp2
    spline = gsl_interp_alloc (gsl_interp_cspline, xr);
    gsl_interp_init (spline, MDRPTR(x), MDRPTR(y), xr);
  }

  if (spline)
  {
    cspline_state_t *state = spline->state;
    w = mdr_Create (xr, 3);
    for (i = 0; i < xr; i++)
    {
      Mdr0 (w, i, 0) = MdrV0 (x, i);
      Mdr0 (w, i, 1) = MdrV0 (y, i);
      Mdr0 (w, i, 2) = state->c[i];
    }
  }
  else
    w = mdr_Create(0,0);

  gsl_interp_free (spline);
  gsl_interp_accel_free (acc);

  // clean stuff
  ent_Clean (X);
  ent_Clean (Y);
  ent_Clean (YP1);
  ent_Clean (YPN);

  return ent_Assign_Rlab_MDR(w);
}

Ent *
ent_cubicspline_eval_limited (int nargs, Datum args[])
{
  Ent *X=0, *M=0, *e3=0;
  MDR *x=0, *m=0, *w=0, *ider=0;

  int xr, xc, nx, na, i, nid=1, j;

  if (nargs < 2)
  {
    printf ("cspleval: Cubic-spline interpolation.\n");
    printf ("cspleval: Format:\n");
    printf ("cspleval:   y = cspleval(x,m,ider),\n");
    printf ("cspleval: where 'x' are the interpolant values, while 'm' is\n");
    printf
      ("cspleval: the matrix of second derivatives found by 'csplinterp'.\n");
    rerror ("requires at least two arguments");
  }
  //
  // get X
  //
  X = bltin_get_ent (args[0]);
  if (ent_type (X) != MATRIX_DENSE_REAL)
    rerror ("cspleval: 'x' has to be real vector !");
  x = ent_data (X);
  if (x->nrow != 1 && x->ncol != 1)
    rerror ("cspleval: 'x' has to be real vector !");

  xr = x->nrow;
  xc = x->ncol;
  nx = xr * xc;

  //
  // get M
  //
  M = bltin_get_ent (args[1]);
  if (ent_type (M) == MATRIX_DENSE_REAL)
    m = ent_data (M);
  else
    rerror
      ("cspleval: 'm' has to be a three-column matrix [x,y(x),y''(x)]!\n");
  if (MNC (m) != 3)
    rerror
      ("cspleval: 'm' has to be a three-column matrix [x,y(x),y''(x)]!\n");
  na = m->nrow;

  //
  // is ther ider?
  //
  if (nargs > 2)
  {
    e3 = bltin_get_ent (args[2]);
    if (ent_type (e3) != MATRIX_DENSE_REAL)
      rerror ("cspleval: 'ider' has to be a vector !");
    ider = ent_data (e3);
    nid = ider->nrow * ider->ncol;
    for(i=0; i<nid; i++)
    {
      if (MdrV0(ider,i)!=0 && MdrV0(ider,i)!=1 && MdrV0(ider,i)!=2)
        rerror("cspleval: 'ider' can have entries 0,1 or 2 !\n");
    }
  }

  // create accelerator
  gsl_interp_accel *acc = gsl_interp_accel_alloc ();

  // allocate fake spline element and pass our element as its coefficients
  gsl_interp *spline=0;
  spline = gsl_interp_alloc (gsl_interp_cspline, m->nrow);
  cspline_state_t *state = spline->state;
  double *old_state_c = state->c;
  state->c = &MdrV0(m,2*na);

  w = mdr_Create(xr*xc, nid);
  if (ider)
    for (j=0; j<nid; j++)
  {
    if (MdrV0(ider,j)==0)
    {
      for (i=0; i<nx; i++)
        Mdr0(w,i,j) = gsl_interp_eval (spline, &MdrV0(m,0), &MdrV0(m,na), MdrV0(x,i), acc);
    }
    else if (MdrV0(ider,j)==1)
    {
      for (i=0; i<nx; i++)
        Mdr0(w,i,j) = gsl_interp_eval_deriv (spline, &MdrV0(m,0), &MdrV0(m,na), MdrV0(x,i), acc);
    }
    else if (MdrV0(ider,j)==2)
    {
      for (i=0; i<nx; i++)
        Mdr0(w,i,j) = gsl_interp_eval_deriv2 (spline, &MdrV0(m,0), &MdrV0(m,na), MdrV0(x,i), acc);
    }
  }
  else
    for (i=0; i<nx; i++)
      MdrV0(w,i) = gsl_interp_eval (spline, &MdrV0(m,0), &MdrV0(m,na), MdrV0(x,i), acc);

  // clean-up GSL
  state->c = old_state_c;
  gsl_interp_free (spline);
  gsl_interp_accel_free (acc);

  // clean-up rlab
  ent_Clean (X);
  ent_Clean (M);
  ent_Clean (e3);

  return ent_Assign_Rlab_MDR(w);
}




void rlabplus_find_piecewise_constant_approx(double *x,
                                            int *size_x, double *m, int *size_m, double *y)
{
  int i, i1, i2, j;
  int nm = *size_m;
  int nx = *size_x;

  i1 = 0;
  while (x[i1] <= m[0] && i1<nx)
  {
    y[i1] = m[nm];
    i1++;
  }

  i2 = nx;
  while (x[i2-1] > m[nm-1] && i2>0)
  {
    y[i2-1] = m[2*nm-1];
    i2--;
  }

  j = 0;
  for (i=i1; i<i2; i++)
  {
    // raise m[j] if necessary
    while (x[i] > m[j])
      j++;

    // now: m[j-1] < x[i] <= m[j]
    if (x[i] < m[j])
      y[i] = m[nm + j-1];
    else
      y[i] = 0.5*(m[nm + j] + m[MAX(nm + j - 1,nm)]);
  }

  return;
}

void rlabplus_find_piecewise_linear_approx(double *x,
                                            int *size_x, double *m, int *size_m, double *y)
{
  int i, i1, i2, j;
  double b;
  int nm = *size_m;
  int nx = *size_x;

  //
  // make linear extrapolation if beyond the range of table
  //
  i1 = 0;
  b  = (m[nm+1] - m[nm]) / (m[1] - m[0]);
  while (x[i1] <= m[0] && i1<nx)
  {
    y[i1] = m[nm] + b * (x[i1] - m[0]);
    i1++;
  }

  i2 = nx;
  b = (m[2*nm-1] - m[2*nm-2]) / (m[nm-1] -  m[nm-2]);
  while (x[i2-1] > m[nm-1] && i2>0)
  {
    y[i2-1] = m[2*nm-1] + b * (x[i2-1] - m[nm-1]);
    i2--;
  }

  j = 0;
  b  = (m[nm+j+1] - m[nm+j]) / (m[j+1] - m[j]);
  for (i=i1; i<i2; i++)
  {
    // raise m[j] if necessary
    while (x[i] > m[j])
    {
      j++;
      b  = (m[nm+j] - m[nm+j-1]) / (m[j] - m[j-1]);
    }

    // now: m[j-1] < x[i] <= m[j]
    y[i] = m[nm+j-1] + b * (x[i] - m[j-1]);
  }

  return;
}

Ent *
ent_lininterp_eval (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *e3=0;

  MDR *x=0, *m=0, *w=0;

  int xr, mr, i, idummy;
  int iorder=1;   // default is the linear interpolation

  ListNode *node;

  if (nargs < 2)
  {
    fprintf (stdout, "linterp: Linear interpolation on a table of values.\n");
    fprintf (stdout, "linterp: Format:\n");
    fprintf (stdout, "linterp:   y = linterp(x,m/,opts/),\n");
    fprintf (stdout,
             "linterp: where 'x' are the values for which the linear\n");
    fprintf (stdout,
             "linterp: interpolation is sought, while 'm' is an interpolation\n");
    fprintf (stdout, "linterp: table of values in format, m=[x1,y(x1);..].\n");
    rerror ("requires two arguments!\n");
  }
  //
  // get X
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) != MATRIX_DENSE_REAL)
    rerror ("linterp: 'x' has to be real vector!\n");
  x = ent_data (e1);
  if (x->ncol != 1 && x->nrow!=1)
    rerror ("linterp: 'x' has to be real vector!\n");

  //
  // get M
  //
  e2 = bltin_get_ent (args[1]);
  if (ent_type (e2) != MATRIX_DENSE_REAL)
    rerror ("linterp: 'm' has to be two-column real matrix [x,y(x)]!\n");
  m = ent_data (e2);
  if (m->ncol != 2)
    rerror ("linterp: 'm' has to be two-column real matrix [x,y(x)]!\n");

  // get opts
  if (nargs>2)
  {
    e3 = bltin_get_ent (args[2]);
    if (ent_type (e3) == BTREE)
    {
      // order of the approximation: const, or linear
      node = btree_FindNode (ent_data (e3), RLAB_NAME_LINTERP_ORDER);
      if (node != 0)
      {
        idummy = (int) class_double (var_ent (node));
        if (idummy==0 || idummy==1)
          iorder = idummy;
      }
    }
    else if (ent_type (e3) == MATRIX_DENSE_REAL)
    {
      idummy = (int) class_double (e3);
      if (idummy==0 || idummy==1)
        iorder = idummy;
    }
  }

  xr = x->nrow * x->ncol;
  mr = m->nrow;
  w  = mdr_Create (x->nrow, x->ncol);

  // if user has provided a single point, then the function
  // is a constant for all x
  if (mr == 1)
  {
    for (i=0; i<xr; i++)
      MdrV0 (w, i) = Mdr0 (m, 0, 1);
  }
  else if (mr > 1)
  {
    if (iorder == 0)
    {
      rlabplus_find_piecewise_constant_approx(MDRPTR(x), &xr, MDRPTR(m), &mr, MDRPTR(w));
    }
    else if (iorder == 1)
    {
      rlabplus_find_piecewise_linear_approx(MDRPTR(x), &xr, MDRPTR(m), &mr, MDRPTR(w));
    }
  }

  // clean up
  ent_Clean (e1);
  ent_Clean (e2);
  ent_Clean (e3);

  return ent_Assign_Rlab_MDR(w);
}


Ent *
ent_lininterp1_eval (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0;

  MDR *x, *m, *w;
  int xr, mr, i, j;
  double xj=0, x0=0, rdx=0;

  if (nargs < 2)
  {
    fprintf (stdout, "linterp1: Linear interpolation on a table of values.\n");
    fprintf (stdout, "linterp1: Format:\n");
    fprintf (stdout, "linterp1:   y = linterp1(x,m),\n");
    fprintf (stdout,
             "linterp1: where 'x' are the values for which the linear\n");
    fprintf (stdout,
             "linterp1: interpolation is sought, while 'm' is an interpolation\n");
    fprintf (stdout,
             "linterp1: table  m=[x1,y(x1);..], where  x1,x2,.. form an uniform mesh.\n");
    rerror ("requires two arguments!\n");
  }

  //
  // get X
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) != MATRIX_DENSE_REAL)
    rerror ("linterp1: 'x' has to be a single-column matrix!\n");
  x = class_matrix_real (e1);
  if (MNC (x) != 1 && MNR(x) != 1)
    rerror ("linterp1: 'x' has to be a single-column matrix!\n");

  //
  // get M
  //
  e2 = bltin_get_ent (args[1]);
  if (ent_type (e2) != MATRIX_DENSE_REAL)
    rerror ("linterp1: 'm' has to be a two-column matrix!\n");
  m = class_matrix_real (e2);
  if (MNC (m) != 2)
    rerror ("linterp1: 'm' has to be a two-column matrix!\n");

  xr = MNR (x) * MNC (x);
  mr = MNR (m);

  w = mdr_Create (MNR (x), MNC (x));

  x0 = Mdr1 (m, 1, 1);
  rdx = 1.0 / (Mdr1 (m, 2, 1) - Mdr1 (m, 1, 1));

  for (i = 1; i <= xr; i++)
  {
    xj = (MdrV1 (x, i) - x0) * rdx + 1;
    j = xj;
    if (j < 1)
    {
      MdrV1 (w, i) = Mdr1 (m, 1, 2);
      continue;
    }
    if (j >= mr)
    {
      MdrV1 (w, i) = Mdr1 (m, mr, 2);
      continue;
    }
    if (j == xj)
    {
      MdrV1 (w, i) = Mdr1 (m, j, 2);
      continue;
    }
    MdrV1 (w, i) = Mdr1 (m, j, 2)
        + (Mdr1  (m, j + 1, 2) - Mdr1 (m, j, 2)) * rdx
        * (MdrV1 (x, i) - Mdr1 (m, j, 1));
  }

  ent_Clean (e1);
  ent_Clean (e2);

  return ent_Assign_Rlab_MDR(w);
}

Ent *
ent_harmsum_eval (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *e3=0;

  MDR *x=0, *m=0, *w=0;
  double y;

  int xr, mr, i, j;

  double (*func)();
  func = sin;

  if (nargs < 2)
  {
    fprintf (stdout, "harmsum: Harmonic interpolation over tabulated values.\n");
    fprintf (stdout, "harmsum: Format:\n");
    fprintf (stdout, "harmsum:   y = harmsum(x,m),\n");
    fprintf (stdout, "harmsum: where 'x' are the values for which the harmonic\n");
    fprintf (stdout, "harmsum: interpolation is sought, while 'm' is a table containing\n");
    fprintf (stdout, "harmsum: m=[omega/,lambda/,A,phi].\n");
    rerror ("requires two arguments!\n");
  }
  //
  // get X
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) != MATRIX_DENSE_REAL)
    rerror ("harmsum: 'x' has to be real vector!\n");
  x = ent_data (e1);
  if (x->ncol != 1 && x->nrow!=1)
    rerror ("harmsum: 'x' has to be real vector!\n");

  //
  // get M
  //
  e2 = bltin_get_ent (args[1]);
  if (ent_type (e2) != MATRIX_DENSE_REAL)
    rerror ("harmsum: 'm' has to be three- or four-column real matrix [omega/,lambda/,A,phi]!\n");
  m = ent_data (e2);
  if (m->ncol != 3 && m->ncol != 4)
    rerror ("harmsum: 'm' has to be three- or four-column real matrix [omega/,lambda/,A,phi]!\n");

  xr = x->nrow * x->ncol;
  mr = m->nrow;
  w = mdr_Create (x->nrow, x->ncol);

  if (m->ncol == 3)
  {
    for (i=0; i<xr; i++)
    {
      y = 0;
      for (j=0; j<mr; j++)
      {
        y += Mdr0(m,j,1) * func( Mdr0(m,j,0) * MdrV0(x,i) + Mdr0(m,j,2) );
      }
      MdrV0(w, i) = y;
    }
  }
  else if (m->ncol == 4)
  {
    for (i=0; i<xr; i++)
    {
      y = 0;
      for (j=0; j<mr; j++)
      {
        y += Mdr0(m,j,2) * exp(Mdr0(m,j,1) * MdrV0(x,i)) * func( Mdr0(m,j,0) * MdrV0(x,i) + Mdr0(m,j,3) );
      }
      MdrV0(w, i) = y;
    }
  }

  // clean up
  ent_Clean (e1);
  ent_Clean (e2);
  ent_Clean (e3);

  return ent_Assign_Rlab_MDR(w);
}
