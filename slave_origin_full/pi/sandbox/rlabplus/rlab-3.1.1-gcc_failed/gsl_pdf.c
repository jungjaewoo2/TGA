// rlabplus (C) 2005-2009 Marijan Kostrun
//
// GSL Science Library - Probabiliry density functions
// Copyright (C) 1996, 1997, 1998, 1999, 2000, 2001, 2002 Gerard Jungman
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or (at
// your option) any later version.
//
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
//

// rlab headers, located in variable $RLAB_SDK
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
// shared object
#include <gsl/gsl_mode.h>
#include <gsl/gsl_precision.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_machine.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

// include: standard C headers
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <math.h>

#define PI 3.1415926
#define GSL_WARNINGS_OFF

#include "gsl_cpf.c"


// ********************************************************************************
// ********************************************************************************
//
// 2-parameter distributions
//
//
// ********************************************************************************
// ********************************************************************************


// ----------------------------------------------------------
//
// GSL distributions
//
// ----------------------------------------------------------

Ent *
ent_gsl_pdf_logistic (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0;
  int nr, nc, nr1, nc1, nr2, nc2, i, j, id, jd;
  double a, x;
  MDR *x1=0, *x2=0;
  MDR *w=0;
  if (nargs!=1 && nargs!=2)
  {
    printf
      ("pdf.logistic: Probability density for a logistic "
       "distribution with scale parameter 'a'.\n");
    printf ("pdf.logistic: Format:\n");
    printf ("pdf.logistic:   pdf.logistic(x,a), or pdf.logistic( [x,a] )\n");
    rerror ("One or two arguments required");
  }

  //
  // get first argument and figure out if it is 'x' or '[x,a]'
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) != MATRIX_DENSE_REAL)
    rerror ("pdf.logistic: argument 'x' has to be MATRIX-DENSE-REAL");
  x1 = class_matrix_real (e1);
  nr1 = nr = x1->nrow;
  nc1 = nc = x1->ncol;
  if (nr1 * nc1 == 0)
    rerror ("pdf.logistic: argument 'x' has to be MATRIX-DENSE-REAL");
  if (nargs == 1)
  {
    if (nr*nc!=2 && nc!=2)
      rerror ("pdf.logistic: Single argument implies [x,a] but it is not provided");

    w = mdr_Create (nr, 1);
    for (i = 0; i < nr; i++)
    {
      x = Mdr0 (x1, i, 0);
      a = Mdr0 (x1, i, 1);
      MdrV0 (w, i) = gsl_ran_logistic_pdf (x, a);
    }
  }
  else if (nargs == 2)
  {
    // get mu
    e2 = bltin_get_ent (args[1]);
    if (ent_type (e2) != MATRIX_DENSE_REAL)
      rerror ("pdf.logistic: argument 'a' has to be MATRIX-DENSE-REAL");
    nr1 = nr;
    nc1 = nc;
    x2 = class_matrix_real (e2);
    nr2 = x2->nrow;
    nc2 = x2->ncol;
    if (nr2 * nc2 == 0)
      rerror ("pdf.logistic: argument 'a' has to be MATRIX-DENSE-REAL");

    // figure out the size of output matrix
    nc = nc > nc2 ? nc : nc2;
    nr = nr > nr2 ? nr : nr2;
    w = mdr_Create (nr, nc);
    // calculate
    for (i = 0; i < nr; i++)
    {
      for (j = 0; j < nc; j++)
      {
        id = MIN(i, nr1-1);
        jd = MIN(j, nc1-1);
        x = Mdr0 (x1, id, jd);
        id = MIN(i, nr2-1);
        jd = MIN(j, nc2-1);
        a = Mdr0 (x2, id, jd);
        Mdr0 (w, i, j) = gsl_ran_logistic_pdf (x, a);
      }
    }
  }
  ent_Clean (e1);
  ent_Clean (e2);

  return ent_Assign_Rlab_MDR(w);
}

Ent *
ent_gsl_pdf_t (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0;
  int nr, nc, nr1, nc1, nr2, nc2, i, j, id, jd;
  double a, x;
  MDR *x1=0, *x2=0;
  MDR *w=0;
  if (nargs!=1 && nargs!=2)
  {
    printf("pdf.t: Probability density for a t-distribution "
           "with 'nu' degrees of freedom.\n");
    printf ("pdf.t: Format:\n");
    printf ("pdf.t:   pdf.t(x,nu), or pdf.t( [x,nu] )\n");
    rerror ("One or two arguments required");
  }

  //
  // get first argument and figure out if it is 'x' or '[x,a]'
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) != MATRIX_DENSE_REAL)
    rerror ("pdf.t: argument 'x' has to be MATRIX-DENSE-REAL");
  x1 = class_matrix_real (e1);
  nr1 = nr = x1->nrow;
  nc1 = nc = x1->ncol;
  if (nr1 * nc1 == 0)
    rerror ("pdf.t: argument 'x' has to be MATRIX-DENSE-REAL");
  if (nargs == 1)
  {
    if (nr*nc!=2 && nc!=2)
      rerror ("pdf.t: Single argument implies [x,nu] but it is not provided");

    w = mdr_Create (nr, 1);
    for (i = 0; i < nr; i++)
    {
      x = Mdr0 (x1, i, 0);
      a = Mdr0 (x1, i, 1);
      MdrV0 (w, i) = gsl_ran_tdist_pdf (x, a);
    }
  }
  else if (nargs == 2)
  {
    // get mu
    e2 = bltin_get_ent (args[1]);
    if (ent_type (e2) != MATRIX_DENSE_REAL)
      rerror ("pdf.t: argument 'nu' has to be MATRIX-DENSE-REAL");
    nr1 = nr;
    nc1 = nc;
    x2 = class_matrix_real (e2);
    nr2 = x2->nrow;
    nc2 = x2->ncol;
    if (nr2 * nc2 == 0)
      rerror ("pdf.t: argument 'nu' has to be MATRIX-DENSE-REAL");

    // figure out the size of output matrix
    nc = nc > nc2 ? nc : nc2;
    nr = nr > nr2 ? nr : nr2;
    w = mdr_Create (nr, nc);
    // calculate
    for (i = 0; i < nr; i++)
    {
      for (j = 0; j < nc; j++)
      {
        id = MIN(i, nr1-1);
        jd = MIN(j, nc1-1);
        x = Mdr0 (x1, id, jd);
        id = MIN(i, nr2-1);
        jd = MIN(j, nc2-1);
        a = Mdr0 (x2, id, jd);
        Mdr0 (w, i, j) = gsl_ran_tdist_pdf (x, a);
      }
    }
  }

  ent_Clean (e1);
  ent_Clean (e2);

  return ent_Assign_Rlab_MDR(w);
}


Ent *
    ent_gsl_pdf_chisq (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0;
  int nr, nc, nr1, nc1, nr2, nc2, i, j, id, jd;
  double a, x;
  MDR *x1=0, *x2=0;
  MDR *w=0;
  if (nargs!=1 && nargs!=2)
  {
    printf
        ("pdf.chisq: Probability density for a Chi-squared distribution with 'nu'\n");
    printf ("pdf.chisq: degrees of freedom.\n");
    printf ("pdf.chisq: Format:\n");
    printf ("pdf.chisq:   pdf.chisq(x,nu), or pdf.chisq( [x,nu] )\n");
    rerror ("One or two arguments required");
  }

  //
  // get first argument and figure out if it is 'x' or '[x,a]'
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) != MATRIX_DENSE_REAL)
    rerror ("pdf.chisq: argument 'x' has to be MATRIX-DENSE-REAL");
  x1 = class_matrix_real (e1);
  nr1 = nr = x1->nrow;
  nc1 = nc = x1->ncol;
  if (nr1 * nc1 == 0)
    rerror ("pdf.chisq: argument 'x' has to be MATRIX-DENSE-REAL");
  if (nargs == 1)
  {
    if (nr*nc!=2 && nc!=2)
      rerror ("pdf.chisq: Single argument implies [x,nu] but it is not provided");

    w = mdr_Create (nr, 1);
    for (i = 0; i < nr; i++)
    {
      x = Mdr0 (x1, i, 0);
      a = Mdr0 (x1, i, 1);
      MdrV0 (w, i) = gsl_ran_chisq_pdf (x, a);
    }
  }
  else if (nargs == 2)
  {
    // get mu
    e2 = bltin_get_ent (args[1]);
    if (ent_type (e2) != MATRIX_DENSE_REAL)
      rerror ("pdf.chisq: argument 'nu' has to be MATRIX-DENSE-REAL");
    nr1 = nr;
    nc1 = nc;
    x2 = class_matrix_real (e2);
    nr2 = x2->nrow;
    nc2 = x2->ncol;
    if (nr2 * nc2 == 0)
      rerror ("pdf.chisq: argument 'nu' has to be MATRIX-DENSE-REAL");

    // figure out the size of output matrix
    nc = nc > nc2 ? nc : nc2;
    nr = nr > nr2 ? nr : nr2;
    w = mdr_Create (nr, nc);
    // calculate
    for (i = 0; i < nr; i++)
    {
      for (j = 0; j < nc; j++)
      {
        id = MIN(i, nr1-1);
        jd = MIN(j, nc1-1);
        x = Mdr0 (x1, id, jd);
        id = MIN(i, nr2-1);
        jd = MIN(j, nc2-1);
        a = Mdr0 (x2, id, jd);
        Mdr0 (w, i, j) = gsl_ran_chisq_pdf (x, a);
      }
    }
  }

  ent_Clean (e1);
  ent_Clean (e2);

  return ent_Assign_Rlab_MDR(w);
}

Ent *
    ent_gsl_pdf_rayleigh (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0;
  int nr, nc, nr1, nc1, nr2, nc2, i, j, id, jd;
  double a, x;
  MDR *x1=0, *x2=0;
  MDR *w=0;
  if (nargs!=1 && nargs!=2)
  {
    printf
        ("pdf.rayleigh: Probability density for a Rayleigh distribution with scale\n");
    printf ("pdf.rayleigh: parameter 'sigma'.\n");
    printf ("pdf.rayleigh: Format:\n");
    printf ("pdf.rayleigh:   pdf.rayleigh(x,sigma), or pdf.rayleigh( [x,sigma] )\n");
    rerror ("One or two arguments required");
  }

  //
  // get first argument and figure out if it is 'x' or '[x,a]'
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) != MATRIX_DENSE_REAL)
    rerror ("pdf.rayleigh: argument 'x' has to be MATRIX-DENSE-REAL");
  x1 = class_matrix_real (e1);
  nr1 = nr = x1->nrow;
  nc1 = nc = x1->ncol;
  if (nr1 * nc1 == 0)
    rerror ("pdf.rayleigh: argument 'x' has to be MATRIX-DENSE-REAL");
  if (nargs == 1)
  {
    if (nr*nc!=2 && nc!=2)
      rerror ("pdf.rayleigh: Single argument implies [x,sigma] but it is not provided");

    w = mdr_Create (nr, 1);
    for (i = 0; i < nr; i++)
    {
      x = Mdr0 (x1, i, 0);
      a = Mdr0 (x1, i, 1);
      MdrV0 (w, i) = gsl_ran_rayleigh_pdf (x, a);
    }
  }
  else if (nargs == 2)
  {
    // get mu
    e2 = bltin_get_ent (args[1]);
    if (ent_type (e2) != MATRIX_DENSE_REAL)
      rerror ("pdf.rayleigh: argument 'sigma' has to be MATRIX-DENSE-REAL");
    nr1 = nr;
    nc1 = nc;
    x2 = class_matrix_real (e2);
    nr2 = x2->nrow;
    nc2 = x2->ncol;
    if (nr2 * nc2 == 0)
      rerror ("pdf.rayleigh: argument 'sigma' has to be MATRIX-DENSE-REAL");

    // figure out the size of output matrix
    nc = nc > nc2 ? nc : nc2;
    nr = nr > nr2 ? nr : nr2;
    w = mdr_Create (nr, nc);
    // calculate
    for (i = 0; i < nr; i++)
    {
      for (j = 0; j < nc; j++)
      {
        id = MIN(i, nr1-1);
        jd = MIN(j, nc1-1);
        x = Mdr0 (x1, id, jd);
        id = MIN(i, nr2-1);
        jd = MIN(j, nc2-1);
        a = Mdr0 (x2, id, jd);
        Mdr0 (w, i, j) = gsl_ran_rayleigh_pdf (x, a);
      }
    }
  }

  ent_Clean (e1);
  ent_Clean (e2);

  return ent_Assign_Rlab_MDR(w);
}

Ent *
    ent_gsl_pdf_cauchy (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0;
  int nr, nc, nr1, nc1, nr2, nc2, i, j, id, jd;
  double a, x;
  MDR *x1=0, *x2=0;
  MDR *w=0;
  if (nargs!=1 && nargs!=2)
  {
    printf
        ("pdf.cauchy: Probability density for a Cauchy distribution with scale\n");
    printf ("pdf.cauchy: parameter 'a'.\n");
    printf ("pdf.cauchy: Format:\n");
    printf ("pdf.cauchy:   pdf.cauchy(x,a), or pdf.cauchy( [x,a] )\n");
    rerror ("One or two arguments required");
  }

  //
  // get first argument and figure out if it is 'x' or '[x,a]'
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) != MATRIX_DENSE_REAL)
    rerror ("pdf.cauchy: argument 'x' has to be MATRIX-DENSE-REAL");
  x1 = class_matrix_real (e1);
  nr1 = nr = x1->nrow;
  nc1 = nc = x1->ncol;
  if (nr1 * nc1 == 0)
    rerror ("pdf.cauchy: argument 'x' has to be MATRIX-DENSE-REAL");
  if (nargs == 1)
  {
    if (nr*nc!=2 && nc!=2)
      rerror ("pdf.cauchy: Single argument implies [x,a] but it is not provided");

    w = mdr_Create (nr, 1);
    for (i = 0; i < nr; i++)
    {
      x = Mdr0 (x1, i, 0);
      a = Mdr0 (x1, i, 1);
      MdrV0 (w, i) = gsl_ran_cauchy_pdf (x, a);
    }
  }
  else if (nargs == 2)
  {
    // get mu
    e2 = bltin_get_ent (args[1]);
    if (ent_type (e2) != MATRIX_DENSE_REAL)
      rerror ("pdf.cauchy: argument 'a' has to be MATRIX-DENSE-REAL");
    nr1 = nr;
    nc1 = nc;
    x2 = class_matrix_real (e2);
    nr2 = x2->nrow;
    nc2 = x2->ncol;
    if (nr2 * nc2 == 0)
      rerror ("pdf.cauchy: argument 'a' has to be MATRIX-DENSE-REAL");

    // figure out the size of output matrix
    nc = nc > nc2 ? nc : nc2;
    nr = nr > nr2 ? nr : nr2;
    w = mdr_Create (nr, nc);
    // calculate
    for (i = 0; i < nr; i++)
    {
      for (j = 0; j < nc; j++)
      {
        id = MIN(i, nr1-1);
        jd = MIN(j, nc1-1);
        x = Mdr0 (x1, id, jd);
        id = MIN(i, nr2-1);
        jd = MIN(j, nc2-1);
        a = Mdr0 (x2, id, jd);
        Mdr0 (w, i, j) = gsl_ran_cauchy_pdf (x, a);
      }
    }
  }

  ent_Clean (e1);
  ent_Clean (e2);

  return ent_Assign_Rlab_MDR(w);
}

Ent *
ent_gsl_pdf_laplace (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0;
  int nr, nc, nr1, nc1, nr2, nc2, i, j, id, jd;
  double a, x;
  MDR *x1=0, *x2=0;
  MDR *w=0;
  if (nargs!=1 && nargs!=2)
  {
    printf
        ("pdf.laplace: Probability density for a Laplace distribution with mean 'a'.\n");
    printf
        ("pdf.laplace: Format:\n");
    printf
        ("pdf.laplace:   pdf.laplace(x,a), or pdf.laplace( [x,a] )\n");
    rerror ("One or two arguments required");
  }

  //
  // get first argument and figure out if it is 'x' or '[x,a]'
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) != MATRIX_DENSE_REAL)
    rerror ("pdf.laplace: argument 'x' has to be MATRIX-DENSE-REAL");
  x1 = class_matrix_real (e1);
  nr1 = nr = x1->nrow;
  nc1 = nc = x1->ncol;
  if (nr1 * nc1 == 0)
    rerror ("pdf.laplace: argument 'x' has to be MATRIX-DENSE-REAL");
  if (nargs == 1)
  {
    if (nr*nc!=2 && nc!=2)
      rerror ("pdf.laplace: Single argument implies [x,a] but it is not provided");

    w = mdr_Create (nr, 1);
    for (i = 0; i < nr; i++)
    {
      x = Mdr0 (x1, i, 0);
      a = Mdr0 (x1, i, 1);
      MdrV0 (w, i) = gsl_ran_laplace_pdf (x, a);
    }
  }
  else if (nargs == 2)
  {
    // get mu
    e2 = bltin_get_ent (args[1]);
    if (ent_type (e2) != MATRIX_DENSE_REAL)
      rerror ("pdf.laplace: argument 'a' has to be MATRIX-DENSE-REAL");
    nr1 = nr;
    nc1 = nc;
    x2 = class_matrix_real (e2);
    nr2 = x2->nrow;
    nc2 = x2->ncol;
    if (nr2 * nc2 == 0)
      rerror ("pdf.laplace: argument 'a' has to be MATRIX-DENSE-REAL");

    // figure out the size of output matrix
    nc = nc > nc2 ? nc : nc2;
    nr = nr > nr2 ? nr : nr2;
    w = mdr_Create (nr, nc);
    // calculate
    for (i = 0; i < nr; i++)
    {
      for (j = 0; j < nc; j++)
      {
        id = MIN(i, nr1-1);
        jd = MIN(j, nc1-1);
        x = Mdr0 (x1, id, jd);
        id = MIN(i, nr2-1);
        jd = MIN(j, nc2-1);
        a = Mdr0 (x2, id, jd);
        Mdr0 (w, i, j) = gsl_ran_laplace_pdf (x, a);
      }
    }
  }

  ent_Clean (e1);
  ent_Clean (e2);

  return ent_Assign_Rlab_MDR(w);
}

Ent *
ent_gsl_pdf_exp (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0;
  int nr, nc, nr1, nc1, nr2, nc2, i, j, id, jd;
  double a, x;
  MDR *x1=0, *x2=0;
  MDR *w=0;
  if (nargs!=1 && nargs!=2)
  {
    printf
        ("pdf.exp: Probability density for exponential distribution with mean mu.\n");
    printf
        ("pdf.exp: Format:\n");
    printf
        ("pdf.exp:   pdf.exp(x,mu), or pdf.exp( [x,mu] )\n");
    rerror ("One or two arguments required");
  }

  //
  // get first argument and figure out if it is 'x' or '[x,a]'
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) != MATRIX_DENSE_REAL)
    rerror ("pdf.exp: argument 'x' has to be MATRIX-DENSE-REAL");
  x1 = class_matrix_real (e1);
  nr1 = nr = x1->nrow;
  nc1 = nc = x1->ncol;
  if (nr * nc == 0)
    rerror ("pdf.exp: argument 'x' has to be MATRIX-DENSE-REAL");
  if (nargs == 1)
  {
    if (nr*nc!=2 && nc!=2)
      rerror ("pdf.exp: Single argument implies [x,a] but it is not provided");

    w = mdr_Create (nr, 1);
    for (i = 0; i < nr; i++)
    {
      x = Mdr0 (x1, i, 0);
      a = Mdr0 (x1, i, 1);
      MdrV0 (w, i) = gsl_ran_exponential_pdf (x, a);
    }
  }
  else if (nargs == 2)
  {
    // get mu
    e2 = bltin_get_ent (args[1]);
    if (ent_type (e2) != MATRIX_DENSE_REAL)
      rerror ("pdf.exp: argument 'a' has to be MATRIX-DENSE-REAL");
    nr1 = nr;
    nc1 = nc;
    x2 = class_matrix_real (e2);
    nr2 = x2->nrow;
    nc2 = x2->ncol;
    if (nr2 * nc2 == 0)
      rerror ("pdf.exp: argument 'a' has to be MATRIX-DENSE-REAL");

    // figure out the size of output matrix
    nc = nc > nc2 ? nc : nc2;
    nr = nr > nr2 ? nr : nr2;
    w = mdr_Create (nr, nc);
    // calculate
    for (i = 0; i < nr; i++)
    {
      for (j = 0; j < nc; j++)
      {
        id = MIN(i, nr1-1);
        jd = MIN(j, nc1-1);
        x = Mdr0 (x1, id, jd);
        id = MIN(i, nr2-1);
        jd = MIN(j, nc2-1);
        a = Mdr0 (x2, id, jd);
        Mdr0 (w, i, j) = gsl_ran_exponential_pdf (x, a);
      }
    }
  }

  ent_Clean (e1);
  ent_Clean (e2);

  return ent_Assign_Rlab_MDR(w);
}

Ent *
ent_gsl_pdf_poisson (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0;
  int nr, nc, nr1, nc1, nr2, nc2, i, j, id, jd;
  double a, x;
  MDR *x1=0, *x2=0;
  MDR *w=0;
  if (nargs!=1 && nargs!=2)
  {
    printf
        ("pdf.poisson: Probability density for a Poisson distribution with mean 'mu'.\n");
    printf ("pdf.poisson: Format:\n");
    printf ("pdf.poisson:   pdf.poisson(k,mu), or pdf.poisson( [k,mu] )\n");
    rerror ("One or two arguments required");
  }

  //
  // get first argument and figure out if it is 'x' or '[x,a]'
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) != MATRIX_DENSE_REAL)
    rerror ("pdf.poisson: argument 'k' has to be MATRIX-DENSE-REAL");
  x1 = class_matrix_real (e1);
  nr1 = nr = x1->nrow;
  nc1 = nc = x1->ncol;
  if (nr1 * nc1 == 0)
    rerror ("pdf.poisson: argument 'k' has to be MATRIX-DENSE-REAL");
  if (nargs == 1)
  {
    if (nr*nc!=2 && nc!=2)
      rerror ("pdf.poisson: Single argument implies [k,mu] but it is not provided");

    w = mdr_Create (nr, 1);
    for (i = 0; i < nr; i++)
    {
      x = Mdr0 (x1, i, 0);
      a = Mdr0 (x1, i, 1);
      MdrV0 (w, i) = gsl_ran_poisson_pdf ((int) x, a);
    }
  }
  else if (nargs == 2)
  {
    // get mu
    e2 = bltin_get_ent (args[1]);
    if (ent_type (e2) != MATRIX_DENSE_REAL)
      rerror ("pdf.poisson: argument 'mu' has to be MATRIX-DENSE-REAL");
    nr1 = nr;
    nc1 = nc;
    x2 = class_matrix_real (e2);
    nr2 = x2->nrow;
    nc2 = x2->ncol;
    if (nr2 * nc2 == 0)
      rerror ("pdf.poisson: argument 'mu' has to be MATRIX-DENSE-REAL");

    // figure out the size of output matrix
    nc = nc > nc2 ? nc : nc2;
    nr = nr > nr2 ? nr : nr2;
    w = mdr_Create (nr, nc);
    // calculate
    for (i = 0; i < nr; i++)
    {
      for (j = 0; j < nc; j++)
      {
        id = MIN(i, nr1-1);
        jd = MIN(j, nc1-1);
        x = Mdr0 (x1, id, jd); // x is integer
        id = MIN(i, nr2-1);
        jd = MIN(j, nc2-1);
        a = Mdr0 (x2, id, jd);
        Mdr0 (w, i, j) = gsl_ran_poisson_pdf ((int) x, a);
      }
    }
  }

  ent_Clean (e1);
  ent_Clean (e2);

  return ent_Assign_Rlab_MDR(w);
}

Ent *
ent_gsl_pdf_geometric (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0;
  int nr, nc, nr1, nc1, nr2, nc2, i, j, id, jd;
  double a, x;
  MDR *x1=0, *x2=0;
  MDR *w=0;
  if (nargs!=1 && nargs!=2)
  {
    printf
        ("pdf.geometric: Probability density for geometric distributions with\n");
    printf ("pdf.geometric: probability parameter 'p'.\n");
    printf ("pdf.geometric: Format:\n");
    printf ("pdf.geometric:   pdf.geometric(k,p), or pdf.geometric( [k,p] )\n");
    rerror ("One or two arguments required");
  }

  //
  // get first argument and figure out if it is 'x' or '[x,a]'
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) != MATRIX_DENSE_REAL)
    rerror ("pdf.geometric: argument 'k' has to be MATRIX-DENSE-REAL");
  x1 = class_matrix_real (e1);
  nr1 = nr = x1->nrow;
  nc1 = nc = x1->ncol;
  if (nr1 * nc1 == 0)
    rerror ("pdf.geometric: argument 'k' has to be MATRIX-DENSE-REAL");
  if (nargs == 1)
  {
    if (nr*nc!=2 && nc!=2)
      rerror ("pdf.geometric: Single argument implies [k,p] but it is not provided");

    w = mdr_Create (nr, 1);
    for (i = 0; i < nr; i++)
    {
      x = Mdr0 (x1, i, 0);
      a = Mdr0 (x1, i, 1);
      MdrV0 (w, i) = gsl_ran_geometric_pdf ((int) x, a);
    }
  }
  else if (nargs == 2)
  {
    // get mu
    e2 = bltin_get_ent (args[1]);
    if (ent_type (e2) != MATRIX_DENSE_REAL)
      rerror ("pdf.geometric: argument 'p' has to be MATRIX-DENSE-REAL");
    nr1 = nr;
    nc1 = nc;
    x2 = class_matrix_real (e2);
    nr2 = x2->nrow;
    nc2 = x2->ncol;
    if (nr2 * nc2 == 0)
      rerror ("pdf.geometric: argument 'p' has to be MATRIX-DENSE-REAL");

    // figure out the size of output matrix
    nc = nc > nc2 ? nc : nc2;
    nr = nr > nr2 ? nr : nr2;
    w = mdr_Create (nr, nc);
    // calculate
    for (i = 0; i < nr; i++)
    {
      for (j = 0; j < nc; j++)
      {
        id = MIN(i, nr1-1);
        jd = MIN(j, nc1-1);
        x = Mdr0 (x1, id, jd);
        id = MIN(i, nr2-1);
        jd = MIN(j, nc2-1);
        a = Mdr0 (x2, id, jd);
        Mdr0 (w, i, j) = gsl_ran_geometric_pdf ((int) x, a);
      }
    }
  }

  ent_Clean (e1);
  ent_Clean (e2);

  return ent_Assign_Rlab_MDR(w);
}

Ent *
ent_gsl_pdf_log (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0;
  int nr, nc, nr1, nc1, nr2, nc2, i, j, id, jd;
  double a, x;
  MDR *x1=0, *x2=0;
  MDR *w=0;
  if (nargs!=1 && nargs!=2)
  {
    printf
        ("pdf.log: Probability density for logarithmic distribution with probability\n");
    printf ("pdf.log: parameter 'p'.\n");
    printf ("pdf.log: Format:\n");
    printf ("pdf.log:   pdf.log(x,p)\n");
    rerror ("One or two arguments required");
  }

  //
  // get first argument and figure out if it is 'x' or '[x,a]'
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) != MATRIX_DENSE_REAL)
    rerror ("pdf.log: argument 'x' has to be MATRIX-DENSE-REAL");
  x1 = class_matrix_real (e1);
  nr1 = nr = x1->nrow;
  nc1 = nc = x1->ncol;
  if (nr1 * nc1 == 0)
    rerror ("pdf.log: argument 'x' has to be MATRIX-DENSE-REAL");
  if (nargs == 1)
  {
    if (nr*nc!=2 && nc!=2)
      rerror ("pdf.log: Single argument implies [x,a] but it is not provided");

    w = mdr_Create (nr, 1);
    for (i = 0; i < nr; i++)
    {
      x = Mdr0 (x1, i, 0);
      a = Mdr0 (x1, i, 1);
      MdrV0 (w, i) = gsl_ran_logarithmic_pdf ((int) x, a);
    }
  }
  else if (nargs == 2)
  {
    // get mu
    e2 = bltin_get_ent (args[1]);
    if (ent_type (e2) != MATRIX_DENSE_REAL)
      rerror ("pdf.log: argument 'a' has to be MATRIX-DENSE-REAL");
    nr1 = nr;
    nc1 = nc;
    x2 = class_matrix_real (e2);
    nr2 = x2->nrow;
    nc2 = x2->ncol;
    if (nr2 * nc2 == 0)
      rerror ("pdf.log: argument 'a' has to be MATRIX-DENSE-REAL");

    // figure out the size of output matrix
    nc = nc > nc2 ? nc : nc2;
    nr = nr > nr2 ? nr : nr2;
    w = mdr_Create (nr, nc);
    // calculate
    for (i = 0; i < nr; i++)
    {
      for (j = 0; j < nc; j++)
      {
        id = MIN(i, nr1-1);
        jd = MIN(j, nc1-1);
        x = Mdr0 (x1, id, jd);
        id = MIN(i, nr2-1);
        jd = MIN(j, nc2-1);
        a = Mdr0 (x2, id, jd);
        Mdr0 (w, i, j) = gsl_ran_logarithmic_pdf ((int) x, a);
      }
    }
  }

  ent_Clean (e1);
  ent_Clean (e2);

  return ent_Assign_Rlab_MDR(w);
}


// ********************************************************************************
// ********************************************************************************
//
//
// 3-parameter distributions
//
//
// ********************************************************************************
// ********************************************************************************

// ----------------------------------------------------------
// non gsl distributions
// ----------------------------------------------------------

double rlabplus_pdf_driftdiffusion(double t, double tbar, double r)
{
  double rval;

  if (t == 0)
    return 0;

  rval =  ABS( r ) / sqrt(PI * t) * exp( -r * r * (tbar - t) * (tbar - t) / t );
  return rval;
}

Ent *
ent_misc_pdf_driftdiffusion (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *e3=0;
  int nr, nc, nr1, nc1, nr2, nc2, nr3, nc3, i, j, id, jd;
  double x, a, b;
  MDR *x1=0, *x2=0, *x3=0;
  MDR *w=0;

  if (nargs!=1 && nargs!=2 && nargs!=3)
  {
    printf
        ("pdf.dd: Drift-diffusion probability density\n");
    printf
        ("pdf.dd: with parameters 'T' and 'r'.\n");
    printf
        ("pdf.dd: Format:\n");
    printf
        ("pdf.dd:   pdf.dd(x,T,r), or pdf.dd(x, [T,r]), or pdf.dd( [x,T,r] )\n");
    rerror ("One, two or three arguments required");
  }

  //
  // get x
  //
  if (nargs == 1)
  {
    // get [x, param_1, param_2]
    e1 = bltin_get_ent (args[0]);
    if (ent_type (e1) != MATRIX_DENSE_REAL)
      rerror ("pdf.dd: argument 'x' has to be MATRIX-DENSE-REAL");
    x1 = ent_data (e1);
    nr1 = nr = x1->nrow;
    nc1 = nc = x1->ncol;
    if (nr * nc == 0)
      rerror ("pdf.dd: argument 'x' has to be MATRIX-DENSE-REAL");
    if (nc!=3)
      rerror ("pdf.dd: Single argument implies [x,T,r] but it is not provided");

    w = mdr_Create (nr, 1);
    for (i = 0; i < nr; i++)
    {
      if (x1->type == RLAB_TYPE_INT32)
      {
        x = (double) Mdi0 (x1, i, 0);
        a = (double) Mdi0 (x1, i, 1);
        b = (double) Mdi0 (x1, i, 2);
      }
      else
      {
        x = Mdr0 (x1, i, 0);
        a = Mdr0 (x1, i, 1);
        b = Mdr0 (x1, i, 2);
      }
      MdrV0 (w, i)  = rlabplus_pdf_driftdiffusion (x, a, b);
    }
  }
  else if (nargs == 2)
  {
    // get x
    e1 = bltin_get_ent (args[0]);
    if (ent_type (e1) != MATRIX_DENSE_REAL)
      rerror ("pdf.dd: argument 'x' has to be MATRIX-DENSE-REAL");
    x1 = ent_data (e1);
    nr1 = nr = x1->nrow;
    nc1 = nc = x1->ncol;
    if (nr * nc == 0)
      rerror ("pdf.dd: argument 'x' has to be MATRIX-DENSE-REAL");

    // get [param_1, param_2]
    e2 = bltin_get_ent (args[1]);
    if (ent_type (e2) != MATRIX_DENSE_REAL)
      rerror ("pdf.dd: argument 'param' has to be MATRIX-DENSE-REAL");
    x2 = ent_data (e2);
    nr2 = x2->nrow;
    nc2 = x2->ncol;
    if (nc2 != 2)
      rerror ("pdf.dd: argument 'param' has to be [param_1, param_2]");

    nr = nr > nr2 ? nr : nr2;
    w = mdr_Create (nr, nc);
    // calculate
    for (i = 0; i < nr; i++)
    {
      for (j = 0; j < nc; j++)
      {
        id = MIN(i, nr1-1);
        jd = MIN(j, nc1-1);
        // x1:
        if (x1->type == RLAB_TYPE_INT32)
          x = (double) Mdi0 (x1, id, jd);
        else
          x = Mdr0 (x1, id, jd);
        // x2:
        id = MIN(i, nr2-1);
        if (x2->type == RLAB_TYPE_INT32)
        {
          a = (double) Mdi0 (x2, id, 0);
          b = (double) Mdi0 (x2, id, 1);
        }
        else
        {
          a = Mdr0 (x2, id, 0);
          b = Mdr0 (x2, id, 1);
        }
        Mdr0 (w, i, j) = rlabplus_pdf_driftdiffusion (x, a, b);
      }
    }
  }
  else if (nargs == 3)
  {
    // get x
    e1 = bltin_get_ent (args[0]);
    if (ent_type (e1) != MATRIX_DENSE_REAL)
      rerror ("pdf.dd: argument 'x' has to be MATRIX-DENSE-REAL");
    x1 = ent_data (e1);
    nr1 = nr = x1->nrow;
    nc1 = nc = x1->ncol;

    //
    // get a
    //
    e2 = bltin_get_ent (args[1]);
    if (ent_type (e2) != MATRIX_DENSE_REAL)
      rerror ("pdf.dd: argument 'T' has to be MATRIX-DENSE-REAL");
    x2 = ent_data (e2);
    nr2 = x2->nrow;
    nc2 = x2->ncol;
    if (nr2 * nc2 == 0)
      rerror ("pdf.dd: argument 'T' has to be MATRIX-DENSE-REAL");

    // update the size of output matrix
    nr = nr > nr2 ? nr : nr2;
    nc = nc > nc2 ? nc : nc2;

    //
    // get b
    //
    e3 = bltin_get_ent (args[2]);
    if (ent_type (e3) != MATRIX_DENSE_REAL)
      rerror ("pdf.dd: argument 'r' has to be MATRIX-DENSE-REAL");
    x3 = ent_data (e3);
    nr3 = x3->nrow;
    nc3 = x3->ncol;
    if (nr3 * nc3 == 0)
      rerror ("pdf.dd: argument 'r' has to be MATRIX-DENSE-REAL");

    // update the size of output matrix to max of all three sizes
    nr = nr > nr3 ? nr : nr3;
    nc = nc > nc3 ? nc : nc3;
    w = mdr_Create (nr, nc);

    // calculate
    for (i = 0; i < nr; i++)
    {
      for (j = 0; j < nc; j++)
      {
        // x:
        id = MIN(i, nr1-1);
        jd = MIN(j, nc1-1);
        if (x1->type == RLAB_TYPE_INT32)
          x = (double) Mdi0 (x1, id, jd);
        else
          x = Mdr0 (x1, id, jd);
        // a:
        id = MIN(i, nr2-1);
        jd = MIN(j, nc2-1);
        if (x2->type == RLAB_TYPE_INT32)
          a  = (double) Mdi0 (x2, id, jd);
        else
          a  = Mdr0 (x2, id, jd);
        // b:
        id = MIN(i, nr3-1);
        jd = MIN(j, nc3-1);
        if (x3->type == RLAB_TYPE_INT32)
          b  = (double) Mdi0 (x3, id, jd);
        else
          b  = Mdr0 (x3, id, jd);

        // w:
        Mdr0 (w, i, j) = rlabplus_pdf_driftdiffusion (x, a, b);
      }
    }
  }

  // clean-up

  ent_Clean (e1);
  ent_Clean (e2);
  ent_Clean (e3);

  return ent_Assign_Rlab_MDR(w);
}

// ----------------------------------------------------------
// GSL distributions
// ----------------------------------------------------------

Ent *
ent_gsl_pdf_normal (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *e3=0;
  int nr, nc, nr1, nc1, nr2, nc2, nr3, nc3, i, j, id, jd;
  double x, a, b;
  MDR *x1=0, *x2=0, *x3=0;
  MDR *w=0;

  if (nargs!=1 && nargs!=2 && nargs!=3)
  {
    printf ("pdf.normal: Probability density for a Gaussian distribution with mean,\n");
    printf ("pdf.normal: and standard deviation.\n");
    printf ("pdf.normal: Format:\n");
    printf ("pdf.normal:   pdf.normal(x, mean, std), or\n");
    printf ("pdf.normal:   pdf.normal(x, [mean,std] ), or\n");
    printf ("pdf.normal:   pdf.normal( [x,mean,std] ) \n");
    rerror ("One, two or three arguments required");
  }

  //
  // get x
  //
  if (nargs == 1)
  {
    // get [x, param_1, param_2]
    e1 = bltin_get_ent (args[0]);
    if (ent_type (e1) != MATRIX_DENSE_REAL)
      rerror ("pdf.normal: argument 'x' has to be MATRIX-DENSE-REAL");
    x1 = ent_data (e1);
    nr1 = nr = x1->nrow;
    nc1 = nc = x1->ncol;
    if (nr * nc == 0)
      rerror ("pdf.normal: argument 'x' has to be MATRIX-DENSE-REAL");
    if (nc!=3)
      rerror ("pdf.normal: Single argument implies [x,mean,std] but it is not provided");

    w = mdr_Create (nr, 1);
    for (i = 0; i < nr; i++)
    {
      if (x1->type == RLAB_TYPE_INT32)
      {
        x = (double) Mdi0 (x1, i, 0);
        a = (double) Mdi0 (x1, i, 1);
        b = (double) Mdi0 (x1, i, 2);
      }
      else
      {
        x = Mdr0 (x1, i, 0);
        a = Mdr0 (x1, i, 1);
        b = Mdr0 (x1, i, 2);
      }
      MdrV0 (w, i)  = gsl_ran_gaussian_pdf (x-a, b);
    }
  }
  else if (nargs == 2)
  {
    // get x
    e1 = bltin_get_ent (args[0]);
    if (ent_type (e1) != MATRIX_DENSE_REAL)
      rerror ("pdf.normal: argument 'x' has to be MATRIX-DENSE-REAL");
    x1 = ent_data (e1);
    nr1 = nr = x1->nrow;
    nc1 = nc = x1->ncol;
    if (nr * nc == 0)
      rerror ("pdf.normal: argument 'x' has to be MATRIX-DENSE-REAL");

    // get [param_1, param_2]
    e2 = bltin_get_ent (args[1]);
    if (ent_type (e2) != MATRIX_DENSE_REAL)
      rerror ("pdf.normal: argument '[mean,std]' has to be MATRIX-DENSE-REAL");
    x2 = ent_data (e2);
    nr2 = x2->nrow;
    nc2 = x2->ncol;
    if (nc2 != 2)
      rerror ("pdf.normal: second argument has to be [mean,std]");

    nr = nr > nr2 ? nr : nr2;
    w = mdr_Create (nr, nc);
    // calculate
    for (i = 0; i < nr; i++)
    {
      for (j = 0; j < nc; j++)
      {
        id = MIN(i, nr1-1);
        jd = MIN(j, nc1-1);
        // x1:
        if (x1->type == RLAB_TYPE_INT32)
          x = (double) Mdi0 (x1, id, jd);
        else
          x = Mdr0 (x1, id, jd);
        // x2:
        id = MIN(i, nr2-1);
        if (x2->type == RLAB_TYPE_INT32)
        {
          a = (double) Mdi0 (x2, id, 0);
          b = (double) Mdi0 (x2, id, 1);
        }
        else
        {
          a = Mdr0 (x2, id, 0);
          b = Mdr0 (x2, id, 1);
        }
        Mdr0 (w, i, j) = gsl_ran_gaussian_pdf (x-a, b);
      }
    }
  }
  else if (nargs == 3)
  {
    // get x
    e1 = bltin_get_ent (args[0]);
    if (ent_type (e1) != MATRIX_DENSE_REAL)
      rerror ("pdf.normal: argument 'x' has to be MATRIX-DENSE-REAL");
    x1 = ent_data (e1);
    nr1 = nr = x1->nrow;
    nc1 = nc = x1->ncol;

    //
    // get a
    //
    e2 = bltin_get_ent (args[1]);
    if (ent_type (e2) != MATRIX_DENSE_REAL)
      rerror ("pdf.normal: argument 'mean' has to be MATRIX-DENSE-REAL");
    x2 = ent_data (e2);
    nr2 = x2->nrow;
    nc2 = x2->ncol;
    if (nr2 * nc2 == 0)
      rerror ("pdf.normal: argument 'mean' has to be MATRIX-DENSE-REAL");

    // update the size of output matrix
    nr = nr > nr2 ? nr : nr2;
    nc = nc > nc2 ? nc : nc2;

    //
    // get b
    //
    e3 = bltin_get_ent (args[2]);
    if (ent_type (e3) != MATRIX_DENSE_REAL)
      rerror ("pdf.normal: argument 'std' has to be MATRIX-DENSE-REAL");
    x3 = ent_data (e3);
    nr3 = x3->nrow;
    nc3 = x3->ncol;
    if (nr3 * nc3 == 0)
      rerror ("pdf.normal: argument 'std' has to be MATRIX-DENSE-REAL");

    // update the size of output matrix to max of all three sizes
    nr = nr > nr3 ? nr : nr3;
    nc = nc > nc3 ? nc : nc3;
    w = mdr_Create (nr, nc);

    // calculate
    for (i = 0; i < nr; i++)
    {
      for (j = 0; j < nc; j++)
      {
        // x:
        id = MIN(i, nr1-1);
        jd = MIN(j, nc1-1);
        if (x1->type == RLAB_TYPE_INT32)
          x = (double) Mdi0 (x1, id, jd);
        else
          x = Mdr0 (x1, id, jd);
        // a:
        id = MIN(i, nr2-1);
        jd = MIN(j, nc2-1);
        if (x2->type == RLAB_TYPE_INT32)
          a  = (double) Mdi0 (x2, id, jd);
        else
          a  = Mdr0 (x2, id, jd);
        // b:
        id = MIN(i, nr3-1);
        jd = MIN(j, nc3-1);
        if (x3->type == RLAB_TYPE_INT32)
          b  = (double) Mdi0 (x3, id, jd);
        else
          b  = Mdr0 (x3, id, jd);

        // w:
        Mdr0 (w, i, j) = gsl_ran_gaussian_pdf (x-a, b);
      }
    }
  }

  // clean-up

  ent_Clean (e1);
  ent_Clean (e2);
  ent_Clean (e3);

  return ent_Assign_Rlab_MDR(w);
}

Ent *
ent_gsl_pdf_uniform (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *e3=0;
  int nr, nc, nr1, nc1, nr2, nc2, nr3, nc3, i, j, id, jd;
  double x, a, b;
  MDR *x1=0, *x2=0, *x3=0;
  MDR *w=0;

  if (nargs!=1 && nargs!=2 && nargs!=3)
  {
    printf
        ("pdf.uniform: Probability density for uniform distribution\n");
    printf ("pdf.uniform: with parameters 'a' and 'b'.\n");
    printf ("pdf.uniform: Format:\n");
    printf ("pdf.uniform:   pdf.uniform(x,a,b)\n");
    printf ("pdf.uniform:   pdf.uniform(x, [a,b] ), or\n");
    printf ("pdf.uniform:   pdf.uniform( [x,a,b] ) \n");
    rerror ("One, two or three arguments required");
  }

  //
  // get x
  //
  if (nargs == 1)
  {
    // get [x, param_1, param_2]
    e1 = bltin_get_ent (args[0]);
    if (ent_type (e1) != MATRIX_DENSE_REAL)
      rerror ("pdf.uniform: argument 'x' has to be MATRIX-DENSE-REAL");
    x1 = ent_data (e1);
    nr1 = nr = x1->nrow;
    nc1 = nc = x1->ncol;
    if (nr * nc == 0)
      rerror ("pdf.uniform: argument 'x' has to be MATRIX-DENSE-REAL");
    if (nc!=3)
      rerror ("pdf.uniform: Single argument implies [x,a,b] but it is not provided");

    w = mdr_Create (nr, 1);
    for (i = 0; i < nr; i++)
    {
      if (x1->type == RLAB_TYPE_INT32)
      {
        x = (double) Mdi0 (x1, i, 0);
        a = (double) Mdi0 (x1, i, 1);
        b = (double) Mdi0 (x1, i, 2);
      }
      else
      {
        x = Mdr0 (x1, i, 0);
        a = Mdr0 (x1, i, 1);
        b = Mdr0 (x1, i, 2);
      }
      MdrV0 (w, i)  = gsl_ran_flat_pdf (x, a, b);
    }
  }
  else if (nargs == 2)
  {
    // get x
    e1 = bltin_get_ent (args[0]);
    if (ent_type (e1) != MATRIX_DENSE_REAL)
      rerror ("pdf.uniform: argument 'x' has to be MATRIX-DENSE-REAL");
    x1 = ent_data (e1);
    nr1 = nr = x1->nrow;
    nc1 = nc = x1->ncol;
    if (nr * nc == 0)
      rerror ("pdf.uniform: argument 'x' has to be MATRIX-DENSE-REAL");

    // get [param_1, param_2]
    e2 = bltin_get_ent (args[1]);
    if (ent_type (e2) != MATRIX_DENSE_REAL)
      rerror ("pdf.uniform: argument '[a,b]' has to be MATRIX-DENSE-REAL");
    x2 = ent_data (e2);
    nr2 = x2->nrow;
    nc2 = x2->ncol;
    if (nc2 != 2)
      rerror ("pdf.uniform: second argument has to be [a,b]");

    nr = nr > nr2 ? nr : nr2;
    w = mdr_Create (nr, nc);
    // calculate
    for (i = 0; i < nr; i++)
    {
      for (j = 0; j < nc; j++)
      {
        id = MIN(i, nr1-1);
        jd = MIN(j, nc1-1);
        // x1:
        if (x1->type == RLAB_TYPE_INT32)
          x = (double) Mdi0 (x1, id, jd);
        else
          x = Mdr0 (x1, id, jd);
        // x2:
        id = MIN(i, nr2-1);
        if (x2->type == RLAB_TYPE_INT32)
        {
          a = (double) Mdi0 (x2, id, 0);
          b = (double) Mdi0 (x2, id, 1);
        }
        else
        {
          a = Mdr0 (x2, id, 0);
          b = Mdr0 (x2, id, 1);
        }
        Mdr0 (w, i, j) = gsl_ran_flat_pdf (x, a, b);
      }
    }
  }
  else if (nargs == 3)
  {
    // get x
    e1 = bltin_get_ent (args[0]);
    if (ent_type (e1) != MATRIX_DENSE_REAL)
      rerror ("pdf.uniform: argument 'x' has to be MATRIX-DENSE-REAL");
    x1 = ent_data (e1);
    nr1 = nr = x1->nrow;
    nc1 = nc = x1->ncol;

    //
    // get a
    //
    e2 = bltin_get_ent (args[1]);
    if (ent_type (e2) != MATRIX_DENSE_REAL)
      rerror ("pdf.uniform: argument 'a' has to be MATRIX-DENSE-REAL");
    x2 = ent_data (e2);
    nr2 = x2->nrow;
    nc2 = x2->ncol;
    if (nr2 * nc2 == 0)
      rerror ("pdf.uniform: argument 'a' has to be MATRIX-DENSE-REAL");

    // update the size of output matrix
    nr = nr > nr2 ? nr : nr2;
    nc = nc > nc2 ? nc : nc2;

    //
    // get b
    //
    e3 = bltin_get_ent (args[2]);
    if (ent_type (e3) != MATRIX_DENSE_REAL)
      rerror ("pdf.uniform: argument 'b' has to be MATRIX-DENSE-REAL");
    x3 = ent_data (e3);
    nr3 = x3->nrow;
    nc3 = x3->ncol;
    if (nr3 * nc3 == 0)
      rerror ("pdf.uniform: argument 'b' has to be MATRIX-DENSE-REAL");

    // update the size of output matrix to max of all three sizes
    nr = nr > nr3 ? nr : nr3;
    nc = nc > nc3 ? nc : nc3;
    w = mdr_Create (nr, nc);

    // calculate
    for (i = 0; i < nr; i++)
    {
      for (j = 0; j < nc; j++)
      {
        // x:
        id = MIN(i, nr1-1);
        jd = MIN(j, nc1-1);
        if (x1->type == RLAB_TYPE_INT32)
          x = (double) Mdi0 (x1, id, jd);
        else
          x = Mdr0 (x1, id, jd);
        // a:
        id = MIN(i, nr2-1);
        jd = MIN(j, nc2-1);
        if (x2->type == RLAB_TYPE_INT32)
          a  = (double) Mdi0 (x2, id, jd);
        else
          a  = Mdr0 (x2, id, jd);
        // b:
        id = MIN(i, nr3-1);
        jd = MIN(j, nc3-1);
        if (x3->type == RLAB_TYPE_INT32)
          b  = (double) Mdi0 (x3, id, jd);
        else
          b  = Mdr0 (x3, id, jd);

        // w:
        Mdr0 (w, i, j) = gsl_ran_flat_pdf (x, a, b);
      }
    }
  }

  // clean-up

  ent_Clean (e1);
  ent_Clean (e2);
  ent_Clean (e3);

  return ent_Assign_Rlab_MDR(w);
}

Ent *
ent_gsl_pdf_gumbel2 (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *e3=0;
  int nr, nc, nr1, nc1, nr2, nc2, nr3, nc3, i, j, id, jd;
  double x, a, b;
  MDR *x1=0, *x2=0, *x3=0;
  MDR *w=0;

  if (nargs!=1 && nargs!=2 && nargs!=3)
  {
    printf
        ("pdf.gumbel2: Probability density for a Type-2 Gumbel distribution\n");
    printf ("pdf.gumbel2: with parameters 'a' and 'b'.\n");
    printf ("pdf.gumbel2: Format:\n");
    printf ("pdf.gumbel2:   pdf.gumbel2(x,a,b)\n");
    printf ("pdf.gumbel2:   pdf.gumbel2(x, [a,b] ), or\n");
    printf ("pdf.gumbel2:   pdf.gumbel2( [x,a,b] ) \n");
    rerror ("One, two or three arguments required");
  }

  //
  // get x
  //
  if (nargs == 1)
  {
    // get [x, param_1, param_2]
    e1 = bltin_get_ent (args[0]);
    if (ent_type (e1) != MATRIX_DENSE_REAL)
      rerror ("pdf.gumbel2: argument 'x' has to be MATRIX-DENSE-REAL");
    x1 = ent_data (e1);
    nr1 = nr = x1->nrow;
    nc1 = nc = x1->ncol;
    if (nr * nc == 0)
      rerror ("pdf.gumbel2: argument 'x' has to be MATRIX-DENSE-REAL");
    if (nc!=3)
      rerror ("pdf.gumbel2: Single argument implies [x,a,b] but it is not provided");

    w = mdr_Create (nr, 1);
    for (i = 0; i < nr; i++)
    {
      if (x1->type == RLAB_TYPE_INT32)
      {
        x = (double) Mdi0 (x1, i, 0);
        a = (double) Mdi0 (x1, i, 1);
        b = (double) Mdi0 (x1, i, 2);
      }
      else
      {
        x = Mdr0 (x1, i, 0);
        a = Mdr0 (x1, i, 1);
        b = Mdr0 (x1, i, 2);
      }
      MdrV0 (w, i)  = gsl_ran_gumbel2_pdf (x, a, b);
    }
  }
  else if (nargs == 2)
  {
    // get x
    e1 = bltin_get_ent (args[0]);
    if (ent_type (e1) != MATRIX_DENSE_REAL)
      rerror ("pdf.gumbel2: argument 'x' has to be MATRIX-DENSE-REAL");
    x1 = ent_data (e1);
    nr1 = nr = x1->nrow;
    nc1 = nc = x1->ncol;
    if (nr * nc == 0)
      rerror ("pdf.gumbel2: argument 'x' has to be MATRIX-DENSE-REAL");

    // get [param_1, param_2]
    e2 = bltin_get_ent (args[1]);
    if (ent_type (e2) != MATRIX_DENSE_REAL)
      rerror ("pdf.gumbel2: argument '[a,b]' has to be MATRIX-DENSE-REAL");
    x2 = ent_data (e2);
    nr2 = x2->nrow;
    nc2 = x2->ncol;
    if (nc2 != 2)
      rerror ("pdf.gumbel2: second argument has to be [a,b]");

    nr = nr > nr2 ? nr : nr2;
    w = mdr_Create (nr, nc);
    // calculate
    for (i = 0; i < nr; i++)
    {
      for (j = 0; j < nc; j++)
      {
        id = MIN(i, nr1-1);
        jd = MIN(j, nc1-1);
        // x1:
        if (x1->type == RLAB_TYPE_INT32)
          x = (double) Mdi0 (x1, id, jd);
        else
          x = Mdr0 (x1, id, jd);
        // x2:
        id = MIN(i, nr2-1);
        if (x2->type == RLAB_TYPE_INT32)
        {
          a = (double) Mdi0 (x2, id, 0);
          b = (double) Mdi0 (x2, id, 1);
        }
        else
        {
          a = Mdr0 (x2, id, 0);
          b = Mdr0 (x2, id, 1);
        }
        Mdr0 (w, i, j) = gsl_ran_gumbel2_pdf (x, a, b);
      }
    }
  }
  else if (nargs == 3)
  {
    // get x
    e1 = bltin_get_ent (args[0]);
    if (ent_type (e1) != MATRIX_DENSE_REAL)
      rerror ("pdf.gumbel2: argument 'x' has to be MATRIX-DENSE-REAL");
    x1 = ent_data (e1);
    nr1 = nr = x1->nrow;
    nc1 = nc = x1->ncol;

    //
    // get a
    //
    e2 = bltin_get_ent (args[1]);
    if (ent_type (e2) != MATRIX_DENSE_REAL)
      rerror ("pdf.gumbel2: argument 'a' has to be MATRIX-DENSE-REAL");
    x2 = ent_data (e2);
    nr2 = x2->nrow;
    nc2 = x2->ncol;
    if (nr2 * nc2 == 0)
      rerror ("pdf.gumbel2: argument 'a' has to be MATRIX-DENSE-REAL");

    // update the size of output matrix
    nr = nr > nr2 ? nr : nr2;
    nc = nc > nc2 ? nc : nc2;

    //
    // get b
    //
    e3 = bltin_get_ent (args[2]);
    if (ent_type (e3) != MATRIX_DENSE_REAL)
      rerror ("pdf.gumbel2: argument 'b' has to be MATRIX-DENSE-REAL");
    x3 = ent_data (e3);
    nr3 = x3->nrow;
    nc3 = x3->ncol;
    if (nr3 * nc3 == 0)
      rerror ("pdf.gumbel2: argument 'b' has to be MATRIX-DENSE-REAL");

    // update the size of output matrix to max of all three sizes
    nr = nr > nr3 ? nr : nr3;
    nc = nc > nc3 ? nc : nc3;
    w = mdr_Create (nr, nc);

    // calculate
    for (i = 0; i < nr; i++)
    {
      for (j = 0; j < nc; j++)
      {
        // x:
        id = MIN(i, nr1-1);
        jd = MIN(j, nc1-1);
        if (x1->type == RLAB_TYPE_INT32)
          x = (double) Mdi0 (x1, id, jd);
        else
          x = Mdr0 (x1, id, jd);
        // a:
        id = MIN(i, nr2-1);
        jd = MIN(j, nc2-1);
        if (x2->type == RLAB_TYPE_INT32)
          a  = (double) Mdi0 (x2, id, jd);
        else
          a  = Mdr0 (x2, id, jd);
        // b:
        id = MIN(i, nr3-1);
        jd = MIN(j, nc3-1);
        if (x3->type == RLAB_TYPE_INT32)
          b  = (double) Mdi0 (x3, id, jd);
        else
          b  = Mdr0 (x3, id, jd);

        // w:
        Mdr0 (w, i, j) = gsl_ran_gumbel2_pdf (x, a, b);
      }
    }
  }

  // clean-up

  ent_Clean (e1);
  ent_Clean (e2);
  ent_Clean (e3);

  return ent_Assign_Rlab_MDR(w);
}

Ent *
ent_gsl_pdf_gumbel1 (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *e3=0;
  int nr, nc, nr1, nc1, nr2, nc2, nr3, nc3, i, j, id, jd;
  double x, a, b;
  MDR *x1=0, *x2=0, *x3=0;
  MDR *w=0;

  if (nargs!=1 && nargs!=2 && nargs!=3)
  {
    printf
        ("pdf.gumbel1: Probability density for a Type-1 Gumbel distribution\n");
    printf ("pdf.gumbel1: with parameters 'a' and 'b'.\n");
    printf ("pdf.gumbel1: Format:\n");
    printf ("pdf.gumbel1:   pdf.gumbel1(x,a,b)\n");
    printf ("pdf.gumbel1:   pdf.gumbel1(x, [a,b] ), or\n");
    printf ("pdf.gumbel1:   pdf.gumbel1( [x,a,b] ) \n");
    rerror ("One, two or three arguments required");
  }

  //
  // get x
  //
  if (nargs == 1)
  {
    // get [x, param_1, param_2]
    e1 = bltin_get_ent (args[0]);
    if (ent_type (e1) != MATRIX_DENSE_REAL)
      rerror ("pdf.gumbel1: argument 'x' has to be MATRIX-DENSE-REAL");
    x1 = ent_data (e1);
    nr1 = nr = x1->nrow;
    nc1 = nc = x1->ncol;
    if (nr * nc == 0)
      rerror ("pdf.gumbel1: argument 'x' has to be MATRIX-DENSE-REAL");
    if (nc!=3)
      rerror ("pdf.gumbel1: Single argument implies [x,a,b] but it is not provided");

    w = mdr_Create (nr, 1);
    for (i = 0; i < nr; i++)
    {
      if (x1->type == RLAB_TYPE_INT32)
      {
        x = (double) Mdi0 (x1, i, 0);
        a = (double) Mdi0 (x1, i, 1);
        b = (double) Mdi0 (x1, i, 2);
      }
      else
      {
        x = Mdr0 (x1, i, 0);
        a = Mdr0 (x1, i, 1);
        b = Mdr0 (x1, i, 2);
      }
      MdrV0 (w, i)  = gsl_ran_gumbel1_pdf (x, a, b);
    }
  }
  else if (nargs == 2)
  {
    // get x
    e1 = bltin_get_ent (args[0]);
    if (ent_type (e1) != MATRIX_DENSE_REAL)
      rerror ("pdf.gumbel1: argument 'x' has to be MATRIX-DENSE-REAL");
    x1 = ent_data (e1);
    nr1 = nr = x1->nrow;
    nc1 = nc = x1->ncol;
    if (nr * nc == 0)
      rerror ("pdf.gumbel1: argument 'x' has to be MATRIX-DENSE-REAL");

    // get [param_1, param_2]
    e2 = bltin_get_ent (args[1]);
    if (ent_type (e2) != MATRIX_DENSE_REAL)
      rerror ("pdf.gumbel1: argument '[a,b]' has to be MATRIX-DENSE-REAL");
    x2 = ent_data (e2);
    nr2 = x2->nrow;
    nc2 = x2->ncol;
    if (nc2 != 2)
      rerror ("pdf.gumbel1: second argument has to be [a,b]");

    nr = nr > nr2 ? nr : nr2;
    w = mdr_Create (nr, nc);
    // calculate
    for (i = 0; i < nr; i++)
    {
      for (j = 0; j < nc; j++)
      {
        id = MIN(i, nr1-1);
        jd = MIN(j, nc1-1);
        // x1:
        if (x1->type == RLAB_TYPE_INT32)
          x = (double) Mdi0 (x1, id, jd);
        else
          x = Mdr0 (x1, id, jd);
        // x2:
        id = MIN(i, nr2-1);
        if (x2->type == RLAB_TYPE_INT32)
        {
          a = (double) Mdi0 (x2, id, 0);
          b = (double) Mdi0 (x2, id, 1);
        }
        else
        {
          a = Mdr0 (x2, id, 0);
          b = Mdr0 (x2, id, 1);
        }
        Mdr0 (w, i, j) = gsl_ran_gumbel1_pdf (x, a, b);
      }
    }
  }
  else if (nargs == 3)
  {
    // get x
    e1 = bltin_get_ent (args[0]);
    if (ent_type (e1) != MATRIX_DENSE_REAL)
      rerror ("pdf.gumbel1: argument 'x' has to be MATRIX-DENSE-REAL");
    x1 = ent_data (e1);
    nr1 = nr = x1->nrow;
    nc1 = nc = x1->ncol;

    //
    // get a
    //
    e2 = bltin_get_ent (args[1]);
    if (ent_type (e2) != MATRIX_DENSE_REAL)
      rerror ("pdf.gumbel1: argument 'a' has to be MATRIX-DENSE-REAL");
    x2 = ent_data (e2);
    nr2 = x2->nrow;
    nc2 = x2->ncol;
    if (nr2 * nc2 == 0)
      rerror ("pdf.gumbel1: argument 'a' has to be MATRIX-DENSE-REAL");

    // update the size of output matrix
    nr = nr > nr2 ? nr : nr2;
    nc = nc > nc2 ? nc : nc2;

    //
    // get b
    //
    e3 = bltin_get_ent (args[2]);
    if (ent_type (e3) != MATRIX_DENSE_REAL)
      rerror ("pdf.gumbel1: argument 'b' has to be MATRIX-DENSE-REAL");
    x3 = ent_data (e3);
    nr3 = x3->nrow;
    nc3 = x3->ncol;
    if (nr3 * nc3 == 0)
      rerror ("pdf.gumbel1: argument 'b' has to be MATRIX-DENSE-REAL");

    // update the size of output matrix to max of all three sizes
    nr = nr > nr3 ? nr : nr3;
    nc = nc > nc3 ? nc : nc3;
    w = mdr_Create (nr, nc);

    // calculate
    for (i = 0; i < nr; i++)
    {
      for (j = 0; j < nc; j++)
      {
        // x:
        id = MIN(i, nr1-1);
        jd = MIN(j, nc1-1);
        if (x1->type == RLAB_TYPE_INT32)
          x = (double) Mdi0 (x1, id, jd);
        else
          x = Mdr0 (x1, id, jd);
        // a:
        id = MIN(i, nr2-1);
        jd = MIN(j, nc2-1);
        if (x2->type == RLAB_TYPE_INT32)
          a  = (double) Mdi0 (x2, id, jd);
        else
          a  = Mdr0 (x2, id, jd);
        // b:
        id = MIN(i, nr3-1);
        jd = MIN(j, nc3-1);
        if (x3->type == RLAB_TYPE_INT32)
          b  = (double) Mdi0 (x3, id, jd);
        else
          b  = Mdr0 (x3, id, jd);

        // w:
        Mdr0 (w, i, j) = gsl_ran_gumbel1_pdf (x, a, b);
      }
    }
  }

  // clean-up

  ent_Clean (e1);
  ent_Clean (e2);
  ent_Clean (e3);

  return ent_Assign_Rlab_MDR(w);
}

Ent *
ent_gsl_pdf_weibull (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *e3=0;
  int nr, nc, nr1, nc1, nr2, nc2, nr3, nc3, i, j, id, jd;
  double x, a, b;
  MDR *x1=0, *x2=0, *x3=0;
  MDR *w=0;

  if (nargs!=1 && nargs!=2 && nargs!=3)
  {
    printf ("pdf.weibull: Probability density for a Weibull distribution\n");
    printf ("pdf.weibull: with parameters 'a' and 'b'.\n");
    printf ("pdf.weibull: Format:\n");
    printf ("pdf.weibull:   pdf.weibull(x,a,b)\n");
    printf ("pdf.weibull:   pdf.weibull(x, [a,b] ), or\n");
    printf ("pdf.weibull:   pdf.weibull( [x,a,b] ) \n");
    rerror ("One, two or three arguments required");
  }

  //
  // get x
  //
  if (nargs == 1)
  {
    // get [x, param_1, param_2]
    e1 = bltin_get_ent (args[0]);
    if (ent_type (e1) != MATRIX_DENSE_REAL)
      rerror ("pdf.weibull: argument 'x' has to be MATRIX-DENSE-REAL");
    x1 = ent_data (e1);
    nr1 = nr = x1->nrow;
    nc1 = nc = x1->ncol;
    if (nr * nc == 0)
      rerror ("pdf.weibull: argument 'x' has to be MATRIX-DENSE-REAL");
    if (nc!=3)
      rerror ("pdf.weibull: Single argument implies [x,a,b] but it is not provided");

    w = mdr_Create (nr, 1);
    for (i = 0; i < nr; i++)
    {
      if (x1->type == RLAB_TYPE_INT32)
      {
        x = (double) Mdi0 (x1, i, 0);
        a = (double) Mdi0 (x1, i, 1);
        b = (double) Mdi0 (x1, i, 2);
      }
      else
      {
        x = Mdr0 (x1, i, 0);
        a = Mdr0 (x1, i, 1);
        b = Mdr0 (x1, i, 2);
      }
      MdrV0 (w, i)  = gsl_ran_weibull_pdf (x, a, b);
    }
  }
  else if (nargs == 2)
  {
    // get x
    e1 = bltin_get_ent (args[0]);
    if (ent_type (e1) != MATRIX_DENSE_REAL)
      rerror ("pdf.weibull: argument 'x' has to be MATRIX-DENSE-REAL");
    x1 = ent_data (e1);
    nr1 = nr = x1->nrow;
    nc1 = nc = x1->ncol;
    if (nr * nc == 0)
      rerror ("pdf.weibull: argument 'x' has to be MATRIX-DENSE-REAL");

    // get [param_1, param_2]
    e2 = bltin_get_ent (args[1]);
    if (ent_type (e2) != MATRIX_DENSE_REAL)
      rerror ("pdf.weibull: argument '[a,b]' has to be MATRIX-DENSE-REAL");
    x2 = ent_data (e2);
    nr2 = x2->nrow;
    nc2 = x2->ncol;
    if (nc2 != 2)
      rerror ("pdf.weibull: second argument has to be [a,b]");

    nr = nr > nr2 ? nr : nr2;
    w = mdr_Create (nr, nc);
    // calculate
    for (i = 0; i < nr; i++)
    {
      for (j = 0; j < nc; j++)
      {
        id = MIN(i, nr1-1);
        jd = MIN(j, nc1-1);
        // x1:
        if (x1->type == RLAB_TYPE_INT32)
          x = (double) Mdi0 (x1, id, jd);
        else
          x = Mdr0 (x1, id, jd);
        // x2:
        id = MIN(i, nr2-1);
        if (x2->type == RLAB_TYPE_INT32)
        {
          a = (double) Mdi0 (x2, id, 0);
          b = (double) Mdi0 (x2, id, 1);
        }
        else
        {
          a = Mdr0 (x2, id, 0);
          b = Mdr0 (x2, id, 1);
        }
        Mdr0 (w, i, j) = gsl_ran_weibull_pdf (x, a, b);
      }
    }
  }
  else if (nargs == 3)
  {
    // get x
    e1 = bltin_get_ent (args[0]);
    if (ent_type (e1) != MATRIX_DENSE_REAL)
      rerror ("pdf.weibull: argument 'x' has to be MATRIX-DENSE-REAL");
    x1 = ent_data (e1);
    nr1 = nr = x1->nrow;
    nc1 = nc = x1->ncol;

    //
    // get a
    //
    e2 = bltin_get_ent (args[1]);
    if (ent_type (e2) != MATRIX_DENSE_REAL)
      rerror ("pdf.weibull: argument 'a' has to be MATRIX-DENSE-REAL");
    x2 = ent_data (e2);
    nr2 = x2->nrow;
    nc2 = x2->ncol;
    if (nr2 * nc2 == 0)
      rerror ("pdf.weibull: argument 'a' has to be MATRIX-DENSE-REAL");

    // update the size of output matrix
    nr = nr > nr2 ? nr : nr2;
    nc = nc > nc2 ? nc : nc2;

    //
    // get b
    //
    e3 = bltin_get_ent (args[2]);
    if (ent_type (e3) != MATRIX_DENSE_REAL)
      rerror ("pdf.weibull: argument 'b' has to be MATRIX-DENSE-REAL");
    x3 = ent_data (e3);
    nr3 = x3->nrow;
    nc3 = x3->ncol;
    if (nr3 * nc3 == 0)
      rerror ("pdf.weibull: argument 'b' has to be MATRIX-DENSE-REAL");

    // update the size of output matrix to max of all three sizes
    nr = nr > nr3 ? nr : nr3;
    nc = nc > nc3 ? nc : nc3;
    w = mdr_Create (nr, nc);

    // calculate
    for (i = 0; i < nr; i++)
    {
      for (j = 0; j < nc; j++)
      {
        // x:
        id = MIN(i, nr1-1);
        jd = MIN(j, nc1-1);
        if (x1->type == RLAB_TYPE_INT32)
          x = (double) Mdi0 (x1, id, jd);
        else
          x = Mdr0 (x1, id, jd);
        // a:
        id = MIN(i, nr2-1);
        jd = MIN(j, nc2-1);
        if (x2->type == RLAB_TYPE_INT32)
          a  = (double) Mdi0 (x2, id, jd);
        else
          a  = Mdr0 (x2, id, jd);
        // b:
        id = MIN(i, nr3-1);
        jd = MIN(j, nc3-1);
        if (x3->type == RLAB_TYPE_INT32)
          b  = (double) Mdi0 (x3, id, jd);
        else
          b  = Mdr0 (x3, id, jd);

        // w:
        Mdr0 (w, i, j) = gsl_ran_weibull_pdf (x, a, b);
      }
    }
  }

  // clean-up

  ent_Clean (e1);
  ent_Clean (e2);
  ent_Clean (e3);

  return ent_Assign_Rlab_MDR(w);
}

Ent *
ent_gsl_pdf_pareto (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *e3=0;
  int nr, nc, nr1, nc1, nr2, nc2, nr3, nc3, i, j, id, jd;
  double x, a, b;
  MDR *x1=0, *x2=0, *x3=0;
  MDR *w=0;

  if (nargs!=1 && nargs!=2 && nargs!=3)
  {
    printf ("pdf.pareto: Probability density for a Pareto distribution\n");
    printf ("pdf.pareto: with parameters 'a' and 'b'.\n");
    printf ("pdf.pareto: Format:\n");
    printf ("pdf.pareto:   pdf.pareto(x,a,b)\n");
    printf ("pdf.pareto:   pdf.pareto(x, [a,b] ), or\n");
    printf ("pdf.pareto:   pdf.pareto( [x,a,b] ) \n");
    rerror ("One, two or three arguments required");
  }

  //
  // get x
  //
  if (nargs == 1)
  {
    // get [x, param_1, param_2]
    e1 = bltin_get_ent (args[0]);
    if (ent_type (e1) != MATRIX_DENSE_REAL)
      rerror ("pdf.pareto: argument 'x' has to be MATRIX-DENSE-REAL");
    x1 = ent_data (e1);
    nr1 = nr = x1->nrow;
    nc1 = nc = x1->ncol;
    if (nr * nc == 0)
      rerror ("pdf.pareto: argument 'x' has to be MATRIX-DENSE-REAL");
    if (nc!=3)
      rerror ("pdf.pareto: Single argument implies [x,a,b] but it is not provided");

    w = mdr_Create (nr, 1);
    for (i = 0; i < nr; i++)
    {
      if (x1->type == RLAB_TYPE_INT32)
      {
        x = (double) Mdi0 (x1, i, 0);
        a = (double) Mdi0 (x1, i, 1);
        b = (double) Mdi0 (x1, i, 2);
      }
      else
      {
        x = Mdr0 (x1, i, 0);
        a = Mdr0 (x1, i, 1);
        b = Mdr0 (x1, i, 2);
      }
      MdrV0 (w, i)  = gsl_ran_pareto_pdf (x, a, b);
    }
  }
  else if (nargs == 2)
  {
    // get x
    e1 = bltin_get_ent (args[0]);
    if (ent_type (e1) != MATRIX_DENSE_REAL)
      rerror ("pdf.pareto: argument 'x' has to be MATRIX-DENSE-REAL");
    x1 = ent_data (e1);
    nr1 = nr = x1->nrow;
    nc1 = nc = x1->ncol;
    if (nr * nc == 0)
      rerror ("pdf.pareto: argument 'x' has to be MATRIX-DENSE-REAL");

    // get [param_1, param_2]
    e2 = bltin_get_ent (args[1]);
    if (ent_type (e2) != MATRIX_DENSE_REAL)
      rerror ("pdf.pareto: argument '[a,b]' has to be MATRIX-DENSE-REAL");
    x2 = ent_data (e2);
    nr2 = x2->nrow;
    nc2 = x2->ncol;
    if (nc2 != 2)
      rerror ("pdf.pareto: second argument has to be [a,b]");

    nr = nr > nr2 ? nr : nr2;
    w = mdr_Create (nr, nc);
    // calculate
    for (i = 0; i < nr; i++)
    {
      for (j = 0; j < nc; j++)
      {
        id = MIN(i, nr1-1);
        jd = MIN(j, nc1-1);
        // x1:
        if (x1->type == RLAB_TYPE_INT32)
          x = (double) Mdi0 (x1, id, jd);
        else
          x = Mdr0 (x1, id, jd);
        // x2:
        id = MIN(i, nr2-1);
        if (x2->type == RLAB_TYPE_INT32)
        {
          a = (double) Mdi0 (x2, id, 0);
          b = (double) Mdi0 (x2, id, 1);
        }
        else
        {
          a = Mdr0 (x2, id, 0);
          b = Mdr0 (x2, id, 1);
        }
        Mdr0 (w, i, j) = gsl_ran_pareto_pdf (x, a, b);
      }
    }
  }
  else if (nargs == 3)
  {
    // get x
    e1 = bltin_get_ent (args[0]);
    if (ent_type (e1) != MATRIX_DENSE_REAL)
      rerror ("pdf.pareto: argument 'x' has to be MATRIX-DENSE-REAL");
    x1 = ent_data (e1);
    nr1 = nr = x1->nrow;
    nc1 = nc = x1->ncol;

    //
    // get a
    //
    e2 = bltin_get_ent (args[1]);
    if (ent_type (e2) != MATRIX_DENSE_REAL)
      rerror ("pdf.pareto: argument 'a' has to be MATRIX-DENSE-REAL");
    x2 = ent_data (e2);
    nr2 = x2->nrow;
    nc2 = x2->ncol;
    if (nr2 * nc2 == 0)
      rerror ("pdf.pareto: argument 'a' has to be MATRIX-DENSE-REAL");

    // update the size of output matrix
    nr = nr > nr2 ? nr : nr2;
    nc = nc > nc2 ? nc : nc2;

    //
    // get b
    //
    e3 = bltin_get_ent (args[2]);
    if (ent_type (e3) != MATRIX_DENSE_REAL)
      rerror ("pdf.pareto: argument 'b' has to be MATRIX-DENSE-REAL");
    x3 = ent_data (e3);
    nr3 = x3->nrow;
    nc3 = x3->ncol;
    if (nr3 * nc3 == 0)
      rerror ("pdf.pareto: argument 'b' has to be MATRIX-DENSE-REAL");

    // update the size of output matrix to max of all three sizes
    nr = nr > nr3 ? nr : nr3;
    nc = nc > nc3 ? nc : nc3;
    w = mdr_Create (nr, nc);

    // calculate
    for (i = 0; i < nr; i++)
    {
      for (j = 0; j < nc; j++)
      {
        // x:
        id = MIN(i, nr1-1);
        jd = MIN(j, nc1-1);
        if (x1->type == RLAB_TYPE_INT32)
          x = (double) Mdi0 (x1, id, jd);
        else
          x = Mdr0 (x1, id, jd);
        // a:
        id = MIN(i, nr2-1);
        jd = MIN(j, nc2-1);
        if (x2->type == RLAB_TYPE_INT32)
          a  = (double) Mdi0 (x2, id, jd);
        else
          a  = Mdr0 (x2, id, jd);
        // b:
        id = MIN(i, nr3-1);
        jd = MIN(j, nc3-1);
        if (x3->type == RLAB_TYPE_INT32)
          b  = (double) Mdi0 (x3, id, jd);
        else
          b  = Mdr0 (x3, id, jd);

        // w:
        Mdr0 (w, i, j) = gsl_ran_pareto_pdf (x, a, b);
      }
    }
  }

  // clean-up

  ent_Clean (e1);
  ent_Clean (e2);
  ent_Clean (e3);

  return ent_Assign_Rlab_MDR(w);
}

Ent *
ent_gsl_pdf_beta (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *e3=0;
  int nr, nc, nr1, nc1, nr2, nc2, nr3, nc3, i, j, id, jd;
  double x, a, b;
  MDR *x1=0, *x2=0, *x3=0;
  MDR *w=0;

  if (nargs!=1 && nargs!=2 && nargs!=3)
  {
    printf ("pdf.beta: Probability density for a Beta distribution\n");
    printf ("pdf.beta: with parameters 'a' and 'b'.\n");
    printf ("pdf.beta: Format:\n");
    printf ("pdf.beta:   pdf.beta(x,a,b)\n");
    printf ("pdf.beta:   pdf.beta(x, [a,b] ), or\n");
    printf ("pdf.beta:   pdf.beta( [x,a,b] ) \n");
    rerror ("One, two or three arguments required");
  }

  //
  // get x
  //
  if (nargs == 1)
  {
    // get [x, param_1, param_2]
    e1 = bltin_get_ent (args[0]);
    if (ent_type (e1) != MATRIX_DENSE_REAL)
      rerror ("pdf.beta: argument 'x' has to be MATRIX-DENSE-REAL");
    x1 = ent_data (e1);
    nr1 = nr = x1->nrow;
    nc1 = nc = x1->ncol;
    if (nr * nc == 0)
      rerror ("pdf.beta: argument 'x' has to be MATRIX-DENSE-REAL");
    if (nc!=3)
      rerror ("pdf.beta: Single argument implies [x,a,b] but it is not provided");

    w = mdr_Create (nr, 1);
    for (i = 0; i < nr; i++)
    {
      if (x1->type == RLAB_TYPE_INT32)
      {
        x = (double) Mdi0 (x1, i, 0);
        a = (double) Mdi0 (x1, i, 1);
        b = (double) Mdi0 (x1, i, 2);
      }
      else
      {
        x = Mdr0 (x1, i, 0);
        a = Mdr0 (x1, i, 1);
        b = Mdr0 (x1, i, 2);
      }
      MdrV0 (w, i)  = gsl_ran_beta_pdf (x, a, b);
    }
  }
  else if (nargs == 2)
  {
    // get x
    e1 = bltin_get_ent (args[0]);
    if (ent_type (e1) != MATRIX_DENSE_REAL)
      rerror ("pdf.beta: argument 'x' has to be MATRIX-DENSE-REAL");
    x1 = ent_data (e1);
    nr1 = nr = x1->nrow;
    nc1 = nc = x1->ncol;
    if (nr * nc == 0)
      rerror ("pdf.beta: argument 'x' has to be MATRIX-DENSE-REAL");

    // get [param_1, param_2]
    e2 = bltin_get_ent (args[1]);
    if (ent_type (e2) != MATRIX_DENSE_REAL)
      rerror ("pdf.beta: argument '[a,b]' has to be MATRIX-DENSE-REAL");
    x2 = ent_data (e2);
    nr2 = x2->nrow;
    nc2 = x2->ncol;
    if (nc2 != 2)
      rerror ("pdf.beta: second argument has to be [a,b]");

    nr = nr > nr2 ? nr : nr2;
    w = mdr_Create (nr, nc);
    // calculate
    for (i = 0; i < nr; i++)
    {
      for (j = 0; j < nc; j++)
      {
        id = MIN(i, nr1-1);
        jd = MIN(j, nc1-1);
        // x1:
        if (x1->type == RLAB_TYPE_INT32)
          x = (double) Mdi0 (x1, id, jd);
        else
          x = Mdr0 (x1, id, jd);
        // x2:
        id = MIN(i, nr2-1);
        if (x2->type == RLAB_TYPE_INT32)
        {
          a = (double) Mdi0 (x2, id, 0);
          b = (double) Mdi0 (x2, id, 1);
        }
        else
        {
          a = Mdr0 (x2, id, 0);
          b = Mdr0 (x2, id, 1);
        }
        Mdr0 (w, i, j) = gsl_ran_beta_pdf (x, a, b);
      }
    }
  }
  else if (nargs == 3)
  {
    // get x
    e1 = bltin_get_ent (args[0]);
    if (ent_type (e1) != MATRIX_DENSE_REAL)
      rerror ("pdf.beta: argument 'x' has to be MATRIX-DENSE-REAL");
    x1 = ent_data (e1);
    nr1 = nr = x1->nrow;
    nc1 = nc = x1->ncol;

    //
    // get a
    //
    e2 = bltin_get_ent (args[1]);
    if (ent_type (e2) != MATRIX_DENSE_REAL)
      rerror ("pdf.beta: argument 'a' has to be MATRIX-DENSE-REAL");
    x2 = ent_data (e2);
    nr2 = x2->nrow;
    nc2 = x2->ncol;
    if (nr2 * nc2 == 0)
      rerror ("pdf.beta: argument 'a' has to be MATRIX-DENSE-REAL");

    // update the size of output matrix
    nr = nr > nr2 ? nr : nr2;
    nc = nc > nc2 ? nc : nc2;

    //
    // get b
    //
    e3 = bltin_get_ent (args[2]);
    if (ent_type (e3) != MATRIX_DENSE_REAL)
      rerror ("pdf.beta: argument 'b' has to be MATRIX-DENSE-REAL");
    x3 = ent_data (e3);
    nr3 = x3->nrow;
    nc3 = x3->ncol;
    if (nr3 * nc3 == 0)
      rerror ("pdf.beta: argument 'b' has to be MATRIX-DENSE-REAL");

    // update the size of output matrix to max of all three sizes
    nr = nr > nr3 ? nr : nr3;
    nc = nc > nc3 ? nc : nc3;
    w = mdr_Create (nr, nc);

    // calculate
    for (i = 0; i < nr; i++)
    {
      for (j = 0; j < nc; j++)
      {
        // x:
        id = MIN(i, nr1-1);
        jd = MIN(j, nc1-1);
        if (x1->type == RLAB_TYPE_INT32)
          x = (double) Mdi0 (x1, id, jd);
        else
          x = Mdr0 (x1, id, jd);
        // a:
        id = MIN(i, nr2-1);
        jd = MIN(j, nc2-1);
        if (x2->type == RLAB_TYPE_INT32)
          a  = (double) Mdi0 (x2, id, jd);
        else
          a  = Mdr0 (x2, id, jd);
        // b:
        id = MIN(i, nr3-1);
        jd = MIN(j, nc3-1);
        if (x3->type == RLAB_TYPE_INT32)
          b  = (double) Mdi0 (x3, id, jd);
        else
          b  = Mdr0 (x3, id, jd);

        // w:
        Mdr0 (w, i, j) = gsl_ran_beta_pdf (x, a, b);
      }
    }
  }

  // clean-up

  ent_Clean (e1);
  ent_Clean (e2);
  ent_Clean (e3);

  return ent_Assign_Rlab_MDR(w);
}

Ent *
ent_gsl_pdf_F (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *e3=0;
  int nr, nc, nr1, nc1, nr2, nc2, nr3, nc3, i, j, id, jd;
  double x, a, b;
  MDR *x1=0, *x2=0, *x3=0;
  MDR *w=0;

  if (nargs!=1 && nargs!=2 && nargs!=3)
  {
    printf ("pdf.F: Probability density for a F-distribution\n");
    printf ("pdf.F: with 'nu1' and 'nu2' degrees of freedom.\n");
    printf ("pdf.F: Format:\n");
    printf ("pdf.F:   pdf.F(x,nu1,nu2)\n");
    printf ("pdf.F:   pdf.F(x, [nu1,nu2] )\n");
    printf ("pdf.F:   pdf.F( [x,nu1,nu2] )\n");
    rerror ("One, two or three arguments required");
  }

  //
  // get x
  //
  if (nargs == 1)
  {
    // get [x, param_1, param_2]
    e1 = bltin_get_ent (args[0]);
    if (ent_type (e1) != MATRIX_DENSE_REAL)
      rerror ("pdf.F: argument 'x' has to be MATRIX-DENSE-REAL");
    x1 = ent_data (e1);
    nr1 = nr = x1->nrow;
    nc1 = nc = x1->ncol;
    if (nr * nc == 0)
      rerror ("pdf.F: argument 'x' has to be MATRIX-DENSE-REAL");
    if (nc!=3)
      rerror ("pdf.F: Single argument implies [x,nu1,nu2] but it is not provided");

    w = mdr_Create (nr, 1);
    for (i = 0; i < nr; i++)
    {
      if (x1->type == RLAB_TYPE_INT32)
      {
        x = (double) Mdi0 (x1, i, 0);
        a = (double) Mdi0 (x1, i, 1);
        b = (double) Mdi0 (x1, i, 2);
      }
      else
      {
        x = Mdr0 (x1, i, 0);
        a = Mdr0 (x1, i, 1);
        b = Mdr0 (x1, i, 2);
      }
      MdrV0 (w, i)  = gsl_ran_fdist_pdf (x, a, b);
    }
  }
  else if (nargs == 2)
  {
    // get x
    e1 = bltin_get_ent (args[0]);
    if (ent_type (e1) != MATRIX_DENSE_REAL)
      rerror ("pdf.F: argument 'x' has to be MATRIX-DENSE-REAL");
    x1 = ent_data (e1);
    nr1 = nr = x1->nrow;
    nc1 = nc = x1->ncol;
    if (nr * nc == 0)
      rerror ("pdf.F: argument 'x' has to be MATRIX-DENSE-REAL");

    // get [param_1, param_2]
    e2 = bltin_get_ent (args[1]);
    if (ent_type (e2) != MATRIX_DENSE_REAL)
      rerror ("pdf.F: argument '[nu1,nu2]' has to be MATRIX-DENSE-REAL");
    x2 = ent_data (e2);
    nr2 = x2->nrow;
    nc2 = x2->ncol;
    if (nc2 != 2)
      rerror ("pdf.F: second argument has to be [nu1,nu2]");

    nr = nr > nr2 ? nr : nr2;
    w = mdr_Create (nr, nc);
    // calculate
    for (i = 0; i < nr; i++)
    {
      for (j = 0; j < nc; j++)
      {
        id = MIN(i, nr1-1);
        jd = MIN(j, nc1-1);
        // x1:
        if (x1->type == RLAB_TYPE_INT32)
          x = (double) Mdi0 (x1, id, jd);
        else
          x = Mdr0 (x1, id, jd);
        // x2:
        id = MIN(i, nr2-1);
        if (x2->type == RLAB_TYPE_INT32)
        {
          a = (double) Mdi0 (x2, id, 0);
          b = (double) Mdi0 (x2, id, 1);
        }
        else
        {
          a = Mdr0 (x2, id, 0);
          b = Mdr0 (x2, id, 1);
        }
        Mdr0 (w, i, j) = gsl_ran_fdist_pdf (x, a, b);
      }
    }
  }
  else if (nargs == 3)
  {
    // get x
    e1 = bltin_get_ent (args[0]);
    if (ent_type (e1) != MATRIX_DENSE_REAL)
      rerror ("pdf.F: argument 'x' has to be MATRIX-DENSE-REAL");
    x1 = ent_data (e1);
    nr1 = nr = x1->nrow;
    nc1 = nc = x1->ncol;

    //
    // get a
    //
    e2 = bltin_get_ent (args[1]);
    if (ent_type (e2) != MATRIX_DENSE_REAL)
      rerror ("pdf.F: argument 'nu1' has to be MATRIX-DENSE-REAL");
    x2 = ent_data (e2);
    nr2 = x2->nrow;
    nc2 = x2->ncol;
    if (nr2 * nc2 == 0)
      rerror ("pdf.F: argument 'nu1' has to be MATRIX-DENSE-REAL");

    // update the size of output matrix
    nr = nr > nr2 ? nr : nr2;
    nc = nc > nc2 ? nc : nc2;

    //
    // get b
    //
    e3 = bltin_get_ent (args[2]);
    if (ent_type (e3) != MATRIX_DENSE_REAL)
      rerror ("pdf.F: argument 'nu2' has to be MATRIX-DENSE-REAL");
    x3 = ent_data (e3);
    nr3 = x3->nrow;
    nc3 = x3->ncol;
    if (nr3 * nc3 == 0)
      rerror ("pdf.F: argument 'nu2' has to be MATRIX-DENSE-REAL");

    // update the size of output matrix to max of all three sizes
    nr = nr > nr3 ? nr : nr3;
    nc = nc > nc3 ? nc : nc3;
    w = mdr_Create (nr, nc);

    // calculate
    for (i = 0; i < nr; i++)
    {
      for (j = 0; j < nc; j++)
      {
        // x:
        id = MIN(i, nr1-1);
        jd = MIN(j, nc1-1);
        if (x1->type == RLAB_TYPE_INT32)
          x = (double) Mdi0 (x1, id, jd);
        else
          x = Mdr0 (x1, id, jd);
        // a:
        id = MIN(i, nr2-1);
        jd = MIN(j, nc2-1);
        if (x2->type == RLAB_TYPE_INT32)
          a  = (double) Mdi0 (x2, id, jd);
        else
          a  = Mdr0 (x2, id, jd);
        // b:
        id = MIN(i, nr3-1);
        jd = MIN(j, nc3-1);
        if (x3->type == RLAB_TYPE_INT32)
          b  = (double) Mdi0 (x3, id, jd);
        else
          b  = Mdr0 (x3, id, jd);

        // w:
        Mdr0 (w, i, j) = gsl_ran_fdist_pdf (x, a, b);
      }
    }
  }

  // clean-up

  ent_Clean (e1);
  ent_Clean (e2);
  ent_Clean (e3);

  return ent_Assign_Rlab_MDR(w);
}

Ent *
ent_gsl_pdf_lognormal (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *e3=0;
  int nr, nc, nr1, nc1, nr2, nc2, nr3, nc3, i, j, id, jd;
  double x, a, b;
  MDR *x1=0, *x2=0, *x3=0;
  MDR *w=0;

  if (nargs!=1 && nargs!=2 && nargs!=3)
  {
    printf
        ("pdf.lognormal: Probability density for a lognormal distribution\n");
    printf ("pdf.lognormal: with parameters 'zeta' and 'sigma'.\n");
    printf ("pdf.lognormal: Format:\n");
    printf ("pdf.lognormal:   pdf.lognormal(x,zeta,sigma)\n");
    printf ("pdf.lognormal:   pdf.lognormal(x, [zeta,sigma] )\n");
    printf ("pdf.lognormal:   pdf.lognormal( [x,zeta,sigma] )\n");
    rerror ("One, two or three arguments required");
  }

  //
  // get x
  //
  if (nargs == 1)
  {
    // get [x, param_1, param_2]
    e1 = bltin_get_ent (args[0]);
    if (ent_type (e1) != MATRIX_DENSE_REAL)
      rerror ("pdf.lognormal: argument 'x' has to be MATRIX-DENSE-REAL");
    x1 = ent_data (e1);
    nr1 = nr = x1->nrow;
    nc1 = nc = x1->ncol;
    if (nr * nc == 0)
      rerror ("pdf.lognormal: argument 'x' has to be MATRIX-DENSE-REAL");
    if (nc!=3)
      rerror ("pdf.lognormal: Single argument implies [x,zeta,sigma] but it is not provided");

    w = mdr_Create (nr, 1);
    for (i = 0; i < nr; i++)
    {
      if (x1->type == RLAB_TYPE_INT32)
      {
        x = (double) Mdi0 (x1, i, 0);
        a = (double) Mdi0 (x1, i, 1);
        b = (double) Mdi0 (x1, i, 2);
      }
      else
      {
        x = Mdr0 (x1, i, 0);
        a = Mdr0 (x1, i, 1);
        b = Mdr0 (x1, i, 2);
      }
      MdrV0 (w, i)  = gsl_ran_lognormal_pdf (x, a, b);
    }
  }
  else if (nargs == 2)
  {
    // get x
    e1 = bltin_get_ent (args[0]);
    if (ent_type (e1) != MATRIX_DENSE_REAL)
      rerror ("pdf.lognormal: argument 'x' has to be MATRIX-DENSE-REAL");
    x1 = ent_data (e1);
    nr1 = nr = x1->nrow;
    nc1 = nc = x1->ncol;
    if (nr * nc == 0)
      rerror ("pdf.lognormal: argument 'x' has to be MATRIX-DENSE-REAL");

    // get [param_1, param_2]
    e2 = bltin_get_ent (args[1]);
    if (ent_type (e2) != MATRIX_DENSE_REAL)
      rerror ("pdf.lognormal: argument '[zeta,sigma]' has to be MATRIX-DENSE-REAL");
    x2 = ent_data (e2);
    nr2 = x2->nrow;
    nc2 = x2->ncol;
    if (nc2 != 2)
      rerror ("pdf.lognormal: second argument has to be [zeta,sigma]");

    nr = nr > nr2 ? nr : nr2;
    w = mdr_Create (nr, nc);
    // calculate
    for (i = 0; i < nr; i++)
    {
      for (j = 0; j < nc; j++)
      {
        id = MIN(i, nr1-1);
        jd = MIN(j, nc1-1);
        // x1:
        if (x1->type == RLAB_TYPE_INT32)
          x = (double) Mdi0 (x1, id, jd);
        else
          x = Mdr0 (x1, id, jd);
        // x2:
        id = MIN(i, nr2-1);
        if (x2->type == RLAB_TYPE_INT32)
        {
          a = (double) Mdi0 (x2, id, 0);
          b = (double) Mdi0 (x2, id, 1);
        }
        else
        {
          a = Mdr0 (x2, id, 0);
          b = Mdr0 (x2, id, 1);
        }
        Mdr0 (w, i, j) = gsl_ran_lognormal_pdf (x, a, b);
      }
    }
  }
  else if (nargs == 3)
  {
    // get x
    e1 = bltin_get_ent (args[0]);
    if (ent_type (e1) != MATRIX_DENSE_REAL)
      rerror ("pdf.lognormal: argument 'x' has to be MATRIX-DENSE-REAL");
    x1 = ent_data (e1);
    nr1 = nr = x1->nrow;
    nc1 = nc = x1->ncol;

    //
    // get a
    //
    e2 = bltin_get_ent (args[1]);
    if (ent_type (e2) != MATRIX_DENSE_REAL)
      rerror ("pdf.lognormal: argument 'zeta' has to be MATRIX-DENSE-REAL");
    x2 = ent_data (e2);
    nr2 = x2->nrow;
    nc2 = x2->ncol;
    if (nr2 * nc2 == 0)
      rerror ("pdf.lognormal: argument 'zeta' has to be MATRIX-DENSE-REAL");

    // update the size of output matrix
    nr = nr > nr2 ? nr : nr2;
    nc = nc > nc2 ? nc : nc2;

    //
    // get b
    //
    e3 = bltin_get_ent (args[2]);
    if (ent_type (e3) != MATRIX_DENSE_REAL)
      rerror ("pdf.lognormal: argument 'sigma' has to be MATRIX-DENSE-REAL");
    x3 = ent_data (e3);
    nr3 = x3->nrow;
    nc3 = x3->ncol;
    if (nr3 * nc3 == 0)
      rerror ("pdf.lognormal: argument 'sigma' has to be MATRIX-DENSE-REAL");

    // update the size of output matrix to max of all three sizes
    nr = nr > nr3 ? nr : nr3;
    nc = nc > nc3 ? nc : nc3;
    w = mdr_Create (nr, nc);

    // calculate
    for (i = 0; i < nr; i++)
    {
      for (j = 0; j < nc; j++)
      {
        // x:
        id = MIN(i, nr1-1);
        jd = MIN(j, nc1-1);
        if (x1->type == RLAB_TYPE_INT32)
          x = (double) Mdi0 (x1, id, jd);
        else
          x = Mdr0 (x1, id, jd);
        // a:
        id = MIN(i, nr2-1);
        jd = MIN(j, nc2-1);
        if (x2->type == RLAB_TYPE_INT32)
          a  = (double) Mdi0 (x2, id, jd);
        else
          a  = Mdr0 (x2, id, jd);
        // b:
        id = MIN(i, nr3-1);
        jd = MIN(j, nc3-1);
        if (x3->type == RLAB_TYPE_INT32)
          b  = (double) Mdi0 (x3, id, jd);
        else
          b  = Mdr0 (x3, id, jd);

        // w:
        Mdr0 (w, i, j) = gsl_ran_lognormal_pdf (x, a, b);
      }
    }
  }

  // clean-up

  ent_Clean (e1);
  ent_Clean (e2);
  ent_Clean (e3);

  return ent_Assign_Rlab_MDR(w);
}

Ent *
ent_gsl_pdf_gamma (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *e3=0;
  int nr, nc, nr1, nc1, nr2, nc2, nr3, nc3, i, j, id, jd;
  double x, a, b;
  MDR *x1=0, *x2=0, *x3=0;
  MDR *w=0;

  if (nargs!=1 && nargs!=2 && nargs!=3)
  {
    printf ("pdf.gamma: Probability density for a Beta distribution\n");
    printf ("pdf.gamma: with parameters 'a' and 'b'.\n");
    printf ("pdf.gamma: Format:\n");
    printf ("pdf.gamma:   pdf.gamma(x,a,b)\n");
    printf ("pdf.gamma:   pdf.gamma(x, [a,b] ), or\n");
    printf ("pdf.gamma:   pdf.gamma( [x,a,b] ) \n");
    rerror ("One, two or three arguments required");
  }

  //
  // get x
  //
  if (nargs == 1)
  {
    // get [x, param_1, param_2]
    e1 = bltin_get_ent (args[0]);
    if (ent_type (e1) != MATRIX_DENSE_REAL)
      rerror ("pdf.gamma: argument 'x' has to be MATRIX-DENSE-REAL");
    x1 = ent_data (e1);
    nr1 = nr = x1->nrow;
    nc1 = nc = x1->ncol;
    if (nr * nc == 0)
      rerror ("pdf.gamma: argument 'x' has to be MATRIX-DENSE-REAL");
    if (nc!=3)
      rerror ("pdf.gamma: Single argument implies [x,a,b] but it is not provided");

    w = mdr_Create (nr, 1);
    for (i = 0; i < nr; i++)
    {
      if (x1->type == RLAB_TYPE_INT32)
      {
        x = (double) Mdi0 (x1, i, 0);
        a = (double) Mdi0 (x1, i, 1);
        b = (double) Mdi0 (x1, i, 2);
      }
      else
      {
        x = Mdr0 (x1, i, 0);
        a = Mdr0 (x1, i, 1);
        b = Mdr0 (x1, i, 2);
      }
      MdrV0 (w, i)  = gsl_ran_gamma_pdf (x, a, b);
    }
  }
  else if (nargs == 2)
  {
    // get x
    e1 = bltin_get_ent (args[0]);
    if (ent_type (e1) != MATRIX_DENSE_REAL)
      rerror ("pdf.gamma: argument 'x' has to be MATRIX-DENSE-REAL");
    x1 = ent_data (e1);
    nr1 = nr = x1->nrow;
    nc1 = nc = x1->ncol;
    if (nr * nc == 0)
      rerror ("pdf.gamma: argument 'x' has to be MATRIX-DENSE-REAL");

    // get [param_1, param_2]
    e2 = bltin_get_ent (args[1]);
    if (ent_type (e2) != MATRIX_DENSE_REAL)
      rerror ("pdf.gamma: argument '[a,b]' has to be MATRIX-DENSE-REAL");
    x2 = ent_data (e2);
    nr2 = x2->nrow;
    nc2 = x2->ncol;
    if (nc2 != 2)
      rerror ("pdf.gamma: second argument has to be [a,b]");

    nr = nr > nr2 ? nr : nr2;
    w = mdr_Create (nr, nc);
    // calculate
    for (i = 0; i < nr; i++)
    {
      for (j = 0; j < nc; j++)
      {
        id = MIN(i, nr1-1);
        jd = MIN(j, nc1-1);
        // x1:
        if (x1->type == RLAB_TYPE_INT32)
          x = (double) Mdi0 (x1, id, jd);
        else
          x = Mdr0 (x1, id, jd);
        // x2:
        id = MIN(i, nr2-1);
        if (x2->type == RLAB_TYPE_INT32)
        {
          a = (double) Mdi0 (x2, id, 0);
          b = (double) Mdi0 (x2, id, 1);
        }
        else
        {
          a = Mdr0 (x2, id, 0);
          b = Mdr0 (x2, id, 1);
        }
        Mdr0 (w, i, j) = gsl_ran_gamma_pdf (x, a, b);
      }
    }
  }
  else if (nargs == 3)
  {
    // get x
    e1 = bltin_get_ent (args[0]);
    if (ent_type (e1) != MATRIX_DENSE_REAL)
      rerror ("pdf.gamma: argument 'x' has to be MATRIX-DENSE-REAL");
    x1 = ent_data (e1);
    nr1 = nr = x1->nrow;
    nc1 = nc = x1->ncol;

    //
    // get a
    //
    e2 = bltin_get_ent (args[1]);
    if (ent_type (e2) != MATRIX_DENSE_REAL)
      rerror ("pdf.gamma: argument 'a' has to be MATRIX-DENSE-REAL");
    x2 = ent_data (e2);
    nr2 = x2->nrow;
    nc2 = x2->ncol;
    if (nr2 * nc2 == 0)
      rerror ("pdf.gamma: argument 'a' has to be MATRIX-DENSE-REAL");

    // update the size of output matrix
    nr = nr > nr2 ? nr : nr2;
    nc = nc > nc2 ? nc : nc2;

    //
    // get b
    //
    e3 = bltin_get_ent (args[2]);
    if (ent_type (e3) != MATRIX_DENSE_REAL)
      rerror ("pdf.gamma: argument 'b' has to be MATRIX-DENSE-REAL");
    x3 = ent_data (e3);
    nr3 = x3->nrow;
    nc3 = x3->ncol;
    if (nr3 * nc3 == 0)
      rerror ("pdf.gamma: argument 'b' has to be MATRIX-DENSE-REAL");

    // update the size of output matrix to max of all three sizes
    nr = nr > nr3 ? nr : nr3;
    nc = nc > nc3 ? nc : nc3;
    w = mdr_Create (nr, nc);

    // calculate
    for (i = 0; i < nr; i++)
    {
      for (j = 0; j < nc; j++)
      {
        // x:
        id = MIN(i, nr1-1);
        jd = MIN(j, nc1-1);
        if (x1->type == RLAB_TYPE_INT32)
          x = (double) Mdi0 (x1, id, jd);
        else
          x = Mdr0 (x1, id, jd);
        // a:
        id = MIN(i, nr2-1);
        jd = MIN(j, nc2-1);
        if (x2->type == RLAB_TYPE_INT32)
          a  = (double) Mdi0 (x2, id, jd);
        else
          a  = Mdr0 (x2, id, jd);
        // b:
        id = MIN(i, nr3-1);
        jd = MIN(j, nc3-1);
        if (x3->type == RLAB_TYPE_INT32)
          b  = (double) Mdi0 (x3, id, jd);
        else
          b  = Mdr0 (x3, id, jd);

        // w:
        Mdr0 (w, i, j) = gsl_ran_gamma_pdf (x, a, b);
      }
    }
  }

  // clean-up

  ent_Clean (e1);
  ent_Clean (e2);
  ent_Clean (e3);

  return ent_Assign_Rlab_MDR(w);
}

Ent *
ent_gsl_pdf_rayleightail (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *e3=0;
  int nr, nc, nr1, nc1, nr2, nc2, nr3, nc3, i, j, id, jd;
  double x, a, b;
  MDR *x1=0, *x2=0, *x3=0;
  MDR *w=0;

  if (nargs!=1 && nargs!=2 && nargs!=3)
  {
    printf
        ("pdf.rayleightail: Probability density for an Rayleigh Tail distribution\n");
    printf
        ("pdf.rayleightail: with scale parameter 'sigma' and lower limit 'a'.\n");
    printf ("pdf.rayleightail: Format:\n");
    printf ("pdf.rayleightail:   pdf.rayleightail(x,a,sigma)\n");
    printf ("pdf.rayleightail:   pdf.rayleightail(x, [a,sigma] )\n");
    printf ("pdf.rayleightail:   pdf.rayleightail( [x,a,sigma] )\n");
    rerror ("One, two or three arguments required");
  }

  //
  // get x
  //
  if (nargs == 1)
  {
    // get [x, param_1, param_2]
    e1 = bltin_get_ent (args[0]);
    if (ent_type (e1) != MATRIX_DENSE_REAL)
      rerror ("pdf.rayleightail: argument 'x' has to be MATRIX-DENSE-REAL");
    x1 = ent_data (e1);
    nr1 = nr = x1->nrow;
    nc1 = nc = x1->ncol;
    if (nr * nc == 0)
      rerror ("pdf.rayleightail: argument 'x' has to be MATRIX-DENSE-REAL");
    if (nc!=3)
      rerror ("pdf.rayleightail: Single argument implies [x,a,sigma] but it is not provided");

    w = mdr_Create (nr, 1);
    for (i = 0; i < nr; i++)
    {
      if (x1->type == RLAB_TYPE_INT32)
      {
        x = (double) Mdi0 (x1, i, 0);
        a = (double) Mdi0 (x1, i, 1);
        b = (double) Mdi0 (x1, i, 2);
      }
      else
      {
        x = Mdr0 (x1, i, 0);
        a = Mdr0 (x1, i, 1);
        b = Mdr0 (x1, i, 2);
      }
      MdrV0 (w, i)  = gsl_ran_rayleigh_tail_pdf (x, a, b);
    }
  }
  else if (nargs == 2)
  {
    // get x
    e1 = bltin_get_ent (args[0]);
    if (ent_type (e1) != MATRIX_DENSE_REAL)
      rerror ("pdf.rayleightail: argument 'x' has to be MATRIX-DENSE-REAL");
    x1 = ent_data (e1);
    nr1 = nr = x1->nrow;
    nc1 = nc = x1->ncol;
    if (nr * nc == 0)
      rerror ("pdf.rayleightail: argument 'x' has to be MATRIX-DENSE-REAL");

    // get [param_1, param_2]
    e2 = bltin_get_ent (args[1]);
    if (ent_type (e2) != MATRIX_DENSE_REAL)
      rerror ("pdf.rayleightail: argument '[a,sigma]' has to be MATRIX-DENSE-REAL");
    x2 = ent_data (e2);
    nr2 = x2->nrow;
    nc2 = x2->ncol;
    if (nc2 != 2)
      rerror ("pdf.rayleightail: second argument has to be [a,sigma]");

    nr = nr > nr2 ? nr : nr2;
    w = mdr_Create (nr, nc);
    // calculate
    for (i = 0; i < nr; i++)
    {
      for (j = 0; j < nc; j++)
      {
        id = MIN(i, nr1-1);
        jd = MIN(j, nc1-1);
        // x1:
        if (x1->type == RLAB_TYPE_INT32)
          x = (double) Mdi0 (x1, id, jd);
        else
          x = Mdr0 (x1, id, jd);
        // x2:
        id = MIN(i, nr2-1);
        if (x2->type == RLAB_TYPE_INT32)
        {
          a = (double) Mdi0 (x2, id, 0);
          b = (double) Mdi0 (x2, id, 1);
        }
        else
        {
          a = Mdr0 (x2, id, 0);
          b = Mdr0 (x2, id, 1);
        }
        Mdr0 (w, i, j) = gsl_ran_rayleigh_tail_pdf (x, a, b);
      }
    }
  }
  else if (nargs == 3)
  {
    // get x
    e1 = bltin_get_ent (args[0]);
    if (ent_type (e1) != MATRIX_DENSE_REAL)
      rerror ("pdf.rayleightail: argument 'x' has to be MATRIX-DENSE-REAL");
    x1 = ent_data (e1);
    nr1 = nr = x1->nrow;
    nc1 = nc = x1->ncol;

    //
    // get a
    //
    e2 = bltin_get_ent (args[1]);
    if (ent_type (e2) != MATRIX_DENSE_REAL)
      rerror ("pdf.rayleightail: argument 'a' has to be MATRIX-DENSE-REAL");
    x2 = ent_data (e2);
    nr2 = x2->nrow;
    nc2 = x2->ncol;
    if (nr2 * nc2 == 0)
      rerror ("pdf.rayleightail: argument 'a' has to be MATRIX-DENSE-REAL");

    // update the size of output matrix
    nr = nr > nr2 ? nr : nr2;
    nc = nc > nc2 ? nc : nc2;

    //
    // get b
    //
    e3 = bltin_get_ent (args[2]);
    if (ent_type (e3) != MATRIX_DENSE_REAL)
      rerror ("pdf.rayleightail: argument 'sigma' has to be MATRIX-DENSE-REAL");
    x3 = ent_data (e3);
    nr3 = x3->nrow;
    nc3 = x3->ncol;
    if (nr3 * nc3 == 0)
      rerror ("pdf.rayleightail: argument 'sigma' has to be MATRIX-DENSE-REAL");

    // update the size of output matrix to max of all three sizes
    nr = nr > nr3 ? nr : nr3;
    nc = nc > nc3 ? nc : nc3;
    w = mdr_Create (nr, nc);

    // calculate
    for (i = 0; i < nr; i++)
    {
      for (j = 0; j < nc; j++)
      {
        // x:
        id = MIN(i, nr1-1);
        jd = MIN(j, nc1-1);
        if (x1->type == RLAB_TYPE_INT32)
          x = (double) Mdi0 (x1, id, jd);
        else
          x = Mdr0 (x1, id, jd);
        // a:
        id = MIN(i, nr2-1);
        jd = MIN(j, nc2-1);
        if (x2->type == RLAB_TYPE_INT32)
          a  = (double) Mdi0 (x2, id, jd);
        else
          a  = Mdr0 (x2, id, jd);
        // b:
        id = MIN(i, nr3-1);
        jd = MIN(j, nc3-1);
        if (x3->type == RLAB_TYPE_INT32)
          b  = (double) Mdi0 (x3, id, jd);
        else
          b  = Mdr0 (x3, id, jd);

        // w:
        Mdr0 (w, i, j) = gsl_ran_rayleigh_tail_pdf (x, a, b);
      }
    }
  }

  // clean-up

  ent_Clean (e1);
  ent_Clean (e2);
  ent_Clean (e3);

  return ent_Assign_Rlab_MDR(w);
}

Ent *
ent_gsl_pdf_exppow (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *e3=0;
  int nr, nc, nr1, nc1, nr2, nc2, nr3, nc3, i, j, id, jd;
  double x, a, b;
  MDR *x1=0, *x2=0, *x3=0;
  MDR *w=0;

  if (nargs!=1 && nargs!=2 && nargs!=3)
  {
    printf ("pdf.exppow: Probability density for a Beta distribution\n");
    printf ("pdf.exppow: with parameters 'a' and 'b'.\n");
    printf ("pdf.exppow: Format:\n");
    printf ("pdf.exppow:   pdf.exppow(x,a,b)\n");
    printf ("pdf.exppow:   pdf.exppow(x, [a,b] ), or\n");
    printf ("pdf.exppow:   pdf.exppow( [x,a,b] ) \n");
    rerror ("One, two or three arguments required");
  }

  //
  // get x
  //
  if (nargs == 1)
  {
    // get [x, param_1, param_2]
    e1 = bltin_get_ent (args[0]);
    if (ent_type (e1) != MATRIX_DENSE_REAL)
      rerror ("pdf.exppow: argument 'x' has to be MATRIX-DENSE-REAL");
    x1 = ent_data (e1);
    nr1 = nr = x1->nrow;
    nc1 = nc = x1->ncol;
    if (nr * nc == 0)
      rerror ("pdf.exppow: argument 'x' has to be MATRIX-DENSE-REAL");
    if (nc!=3)
      rerror ("pdf.exppow: Single argument implies [x,a,b] but it is not provided");

    w = mdr_Create (nr, 1);
    for (i = 0; i < nr; i++)
    {
      if (x1->type == RLAB_TYPE_INT32)
      {
        x = (double) Mdi0 (x1, i, 0);
        a = (double) Mdi0 (x1, i, 1);
        b = (double) Mdi0 (x1, i, 2);
      }
      else
      {
        x = Mdr0 (x1, i, 0);
        a = Mdr0 (x1, i, 1);
        b = Mdr0 (x1, i, 2);
      }
      MdrV0 (w, i)  = gsl_ran_exppow_pdf (x, a, b);
    }
  }
  else if (nargs == 2)
  {
    // get x
    e1 = bltin_get_ent (args[0]);
    if (ent_type (e1) != MATRIX_DENSE_REAL)
      rerror ("pdf.exppow: argument 'x' has to be MATRIX-DENSE-REAL");
    x1 = ent_data (e1);
    nr1 = nr = x1->nrow;
    nc1 = nc = x1->ncol;
    if (nr * nc == 0)
      rerror ("pdf.exppow: argument 'x' has to be MATRIX-DENSE-REAL");

    // get [param_1, param_2]
    e2 = bltin_get_ent (args[1]);
    if (ent_type (e2) != MATRIX_DENSE_REAL)
      rerror ("pdf.exppow: argument '[a,b]' has to be MATRIX-DENSE-REAL");
    x2 = ent_data (e2);
    nr2 = x2->nrow;
    nc2 = x2->ncol;
    if (nc2 != 2)
      rerror ("pdf.exppow: second argument has to be [a,b]");

    nr = nr > nr2 ? nr : nr2;
    w = mdr_Create (nr, nc);
    // calculate
    for (i = 0; i < nr; i++)
    {
      for (j = 0; j < nc; j++)
      {
        id = MIN(i, nr1-1);
        jd = MIN(j, nc1-1);
        // x1:
        if (x1->type == RLAB_TYPE_INT32)
          x = (double) Mdi0 (x1, id, jd);
        else
          x = Mdr0 (x1, id, jd);
        // x2:
        id = MIN(i, nr2-1);
        if (x2->type == RLAB_TYPE_INT32)
        {
          a = (double) Mdi0 (x2, id, 0);
          b = (double) Mdi0 (x2, id, 1);
        }
        else
        {
          a = Mdr0 (x2, id, 0);
          b = Mdr0 (x2, id, 1);
        }
        Mdr0 (w, i, j) = gsl_ran_exppow_pdf (x, a, b);
      }
    }
  }
  else if (nargs == 3)
  {
    // get x
    e1 = bltin_get_ent (args[0]);
    if (ent_type (e1) != MATRIX_DENSE_REAL)
      rerror ("pdf.exppow: argument 'x' has to be MATRIX-DENSE-REAL");
    x1 = ent_data (e1);
    nr1 = nr = x1->nrow;
    nc1 = nc = x1->ncol;

    //
    // get a
    //
    e2 = bltin_get_ent (args[1]);
    if (ent_type (e2) != MATRIX_DENSE_REAL)
      rerror ("pdf.exppow: argument 'a' has to be MATRIX-DENSE-REAL");
    x2 = ent_data (e2);
    nr2 = x2->nrow;
    nc2 = x2->ncol;
    if (nr2 * nc2 == 0)
      rerror ("pdf.exppow: argument 'a' has to be MATRIX-DENSE-REAL");

    // update the size of output matrix
    nr = nr > nr2 ? nr : nr2;
    nc = nc > nc2 ? nc : nc2;

    //
    // get b
    //
    e3 = bltin_get_ent (args[2]);
    if (ent_type (e3) != MATRIX_DENSE_REAL)
      rerror ("pdf.exppow: argument 'b' has to be MATRIX-DENSE-REAL");
    x3 = ent_data (e3);
    nr3 = x3->nrow;
    nc3 = x3->ncol;
    if (nr3 * nc3 == 0)
      rerror ("pdf.exppow: argument 'b' has to be MATRIX-DENSE-REAL");

    // update the size of output matrix to max of all three sizes
    nr = nr > nr3 ? nr : nr3;
    nc = nc > nc3 ? nc : nc3;
    w = mdr_Create (nr, nc);

    // calculate
    for (i = 0; i < nr; i++)
    {
      for (j = 0; j < nc; j++)
      {
        // x:
        id = MIN(i, nr1-1);
        jd = MIN(j, nc1-1);
        if (x1->type == RLAB_TYPE_INT32)
          x = (double) Mdi0 (x1, id, jd);
        else
          x = Mdr0 (x1, id, jd);
        // a:
        id = MIN(i, nr2-1);
        jd = MIN(j, nc2-1);
        if (x2->type == RLAB_TYPE_INT32)
          a  = (double) Mdi0 (x2, id, jd);
        else
          a  = Mdr0 (x2, id, jd);
        // b:
        id = MIN(i, nr3-1);
        jd = MIN(j, nc3-1);
        if (x3->type == RLAB_TYPE_INT32)
          b  = (double) Mdi0 (x3, id, jd);
        else
          b  = Mdr0 (x3, id, jd);

        // w:
        Mdr0 (w, i, j) = gsl_ran_exppow_pdf (x, a, b);
      }
    }
  }

  // clean-up

  ent_Clean (e1);
  ent_Clean (e2);
  ent_Clean (e3);

  return ent_Assign_Rlab_MDR(w);
}

Ent *
ent_gsl_pdf_normaltail (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *e3=0;
  int nr, nc, nr1, nc1, nr2, nc2, nr3, nc3, i, j, id, jd;
  double x, a, b;
  MDR *x1=0, *x2=0, *x3=0;
  MDR *w=0;

  if (nargs!=1 && nargs!=2 && nargs!=3)
  {
    printf
        ("pdf.normaltail: Probability density for a Gaussian tail distribution\n");
    printf
        ("pdf.normaltail: with standard deviation 'sigma' and lower limit 'a'.\n");
    printf ("pdf.normaltail: Format:\n");
    printf ("pdf.normaltail:   pdf.normaltail(x,a,sigma)\n");
    printf ("pdf.normaltail:   pdf.normaltail(x, [a,sigma] )\n");
    printf ("pdf.normaltail:   pdf.normaltail( [x,a,sigma] )\n");
    rerror ("One, two or three arguments required");
  }

  //
  // get x
  //
  if (nargs == 1)
  {
    // get [x, param_1, param_2]
    e1 = bltin_get_ent (args[0]);
    if (ent_type (e1) != MATRIX_DENSE_REAL)
      rerror ("pdf.normaltail: argument 'x' has to be MATRIX-DENSE-REAL");
    x1 = ent_data (e1);
    nr1 = nr = x1->nrow;
    nc1 = nc = x1->ncol;
    if (nr * nc == 0)
      rerror ("pdf.normaltail: argument 'x' has to be MATRIX-DENSE-REAL");
    if (nc!=3)
      rerror ("pdf.normaltail: Single argument implies [x,a,sigma] but it is not provided");

    w = mdr_Create (nr, 1);
    for (i = 0; i < nr; i++)
    {
      if (x1->type == RLAB_TYPE_INT32)
      {
        x = (double) Mdi0 (x1, i, 0);
        a = (double) Mdi0 (x1, i, 1);
        b = (double) Mdi0 (x1, i, 2);
      }
      else
      {
        x = Mdr0 (x1, i, 0);
        a = Mdr0 (x1, i, 1);
        b = Mdr0 (x1, i, 2);
      }
      MdrV0 (w, i)  = gsl_ran_gaussian_tail_pdf (x, a, b);
    }
  }
  else if (nargs == 2)
  {
    // get x
    e1 = bltin_get_ent (args[0]);
    if (ent_type (e1) != MATRIX_DENSE_REAL)
      rerror ("pdf.normaltail: argument 'x' has to be MATRIX-DENSE-REAL");
    x1 = ent_data (e1);
    nr1 = nr = x1->nrow;
    nc1 = nc = x1->ncol;
    if (nr * nc == 0)
      rerror ("pdf.normaltail: argument 'x' has to be MATRIX-DENSE-REAL");

    // get [param_1, param_2]
    e2 = bltin_get_ent (args[1]);
    if (ent_type (e2) != MATRIX_DENSE_REAL)
      rerror ("pdf.normaltail: argument '[a,sigma]' has to be MATRIX-DENSE-REAL");
    x2 = ent_data (e2);
    nr2 = x2->nrow;
    nc2 = x2->ncol;
    if (nc2 != 2)
      rerror ("pdf.normaltail: second argument has to be [a,sigma]");

    nr = nr > nr2 ? nr : nr2;
    w = mdr_Create (nr, nc);
    // calculate
    for (i = 0; i < nr; i++)
    {
      for (j = 0; j < nc; j++)
      {
        id = MIN(i, nr1-1);
        jd = MIN(j, nc1-1);
        // x1:
        if (x1->type == RLAB_TYPE_INT32)
          x = (double) Mdi0 (x1, id, jd);
        else
          x = Mdr0 (x1, id, jd);
        // x2:
        id = MIN(i, nr2-1);
        if (x2->type == RLAB_TYPE_INT32)
        {
          a = (double) Mdi0 (x2, id, 0);
          b = (double) Mdi0 (x2, id, 1);
        }
        else
        {
          a = Mdr0 (x2, id, 0);
          b = Mdr0 (x2, id, 1);
        }
        Mdr0 (w, i, j) = gsl_ran_gaussian_tail_pdf (x, a, b);
      }
    }
  }
  else if (nargs == 3)
  {
    // get x
    e1 = bltin_get_ent (args[0]);
    if (ent_type (e1) != MATRIX_DENSE_REAL)
      rerror ("pdf.normaltail: argument 'x' has to be MATRIX-DENSE-REAL");
    x1 = ent_data (e1);
    nr1 = nr = x1->nrow;
    nc1 = nc = x1->ncol;

    //
    // get a
    //
    e2 = bltin_get_ent (args[1]);
    if (ent_type (e2) != MATRIX_DENSE_REAL)
      rerror ("pdf.normaltail: argument 'a' has to be MATRIX-DENSE-REAL");
    x2 = ent_data (e2);
    nr2 = x2->nrow;
    nc2 = x2->ncol;
    if (nr2 * nc2 == 0)
      rerror ("pdf.normaltail: argument 'a' has to be MATRIX-DENSE-REAL");

    // update the size of output matrix
    nr = nr > nr2 ? nr : nr2;
    nc = nc > nc2 ? nc : nc2;

    //
    // get b
    //
    e3 = bltin_get_ent (args[2]);
    if (ent_type (e3) != MATRIX_DENSE_REAL)
      rerror ("pdf.normaltail: argument 'sigma' has to be MATRIX-DENSE-REAL");
    x3 = ent_data (e3);
    nr3 = x3->nrow;
    nc3 = x3->ncol;
    if (nr3 * nc3 == 0)
      rerror ("pdf.normaltail: argument 'sigma' has to be MATRIX-DENSE-REAL");

    // update the size of output matrix to max of all three sizes
    nr = nr > nr3 ? nr : nr3;
    nc = nc > nc3 ? nc : nc3;
    w = mdr_Create (nr, nc);

    // calculate
    for (i = 0; i < nr; i++)
    {
      for (j = 0; j < nc; j++)
      {
        // x:
        id = MIN(i, nr1-1);
        jd = MIN(j, nc1-1);
        if (x1->type == RLAB_TYPE_INT32)
          x = (double) Mdi0 (x1, id, jd);
        else
          x = Mdr0 (x1, id, jd);
        // a:
        id = MIN(i, nr2-1);
        jd = MIN(j, nc2-1);
        if (x2->type == RLAB_TYPE_INT32)
          a  = (double) Mdi0 (x2, id, jd);
        else
          a  = Mdr0 (x2, id, jd);
        // b:
        id = MIN(i, nr3-1);
        jd = MIN(j, nc3-1);
        if (x3->type == RLAB_TYPE_INT32)
          b  = (double) Mdi0 (x3, id, jd);
        else
          b  = Mdr0 (x3, id, jd);

        // w:
        Mdr0 (w, i, j) = gsl_ran_gaussian_tail_pdf (x, a, b);
      }
    }
  }

  // clean-up

  ent_Clean (e1);
  ent_Clean (e2);
  ent_Clean (e3);

  return ent_Assign_Rlab_MDR(w);
}

Ent *
ent_gsl_pdf_binomial (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *e3=0;
  int nr, nc, nr1, nc1, nr2, nc2, nr3, nc3, i, j, id, jd;
  double x, a, b;
  MDR *x1=0, *x2=0, *x3=0;
  MDR *w=0;

  if (nargs!=1 && nargs!=3)
  {
    printf ("pdf.binomial: Probability density for a binomial distribution\n");
    printf
        ("pdf.binomial: with parameters: 'p', probability of a single trial,\n");
    printf ("pdf.binomial: and 'n', number of independent trials.\n");
    printf ("pdf.binomial: Format:\n");
    printf ("pdf.binomial:   pdf.binomial(k,p,n)\n");
    printf ("pdf.binomial:   pdf.binomial( [k,p,n] )\n");
    rerror ("One, two or three arguments required");
  }

  //
  // get x
  //
  if (nargs == 1)
  {
    // get [x, param_1, param_2]
    e1 = bltin_get_ent (args[0]);
    if (ent_type (e1) != MATRIX_DENSE_REAL)
      rerror ("pdf.binomial: first argument has to be MATRIX-DENSE-REAL");
    x1 = ent_data (e1);
    nr1 = nr = x1->nrow;
    nc1 = nc = x1->ncol;
    if (nr * nc == 0)
      rerror ("pdf.binomial: first argument has to be MATRIX-DENSE-REAL");
    if (nc!=3)
      rerror ("pdf.binomial: Single argument implies [k,p,n] but it is not provided");

    w = mdr_Create (nr, 1);
    for (i = 0; i < nr; i++)
    {
      if (x1->type == RLAB_TYPE_INT32)
      {
        x = (double) Mdi0 (x1, i, 0);
        a = (double) Mdi0 (x1, i, 1);
        b = (double) Mdi0 (x1, i, 2);
      }
      else
      {
        x = Mdr0 (x1, i, 0);
        a = Mdr0 (x1, i, 1);
        b = Mdr0 (x1, i, 2);
      }
      MdrV0 (w, i) = gsl_ran_binomial_pdf ( (int) x, a, (int) b);
    }
  }
  else if (nargs == 3)
  {
    // get x
    e1 = bltin_get_ent (args[0]);
    if (ent_type (e1) != MATRIX_DENSE_REAL)
      rerror ("pdf.binomial: argument 'k' has to be MATRIX-DENSE-REAL");
    x1 = ent_data (e1);
    nr1 = nr = x1->nrow;
    nc1 = nc = x1->ncol;

    //
    // get a
    //
    e2 = bltin_get_ent (args[1]);
    if (ent_type (e2) != MATRIX_DENSE_REAL)
      rerror ("pdf.binomial: argument 'p' has to be MATRIX-DENSE-REAL");
    x2 = ent_data (e2);
    nr2 = x2->nrow;
    nc2 = x2->ncol;
    if (nr2 * nc2 == 0)
      rerror ("pdf.binomial: argument 'p' has to be MATRIX-DENSE-REAL");

    // update the size of output matrix
    nr = nr > nr2 ? nr : nr2;
    nc = nc > nc2 ? nc : nc2;

    //
    // get b
    //
    e3 = bltin_get_ent (args[2]);
    if (ent_type (e3) != MATRIX_DENSE_REAL)
      rerror ("pdf.binomial: argument 'n' has to be MATRIX-DENSE-REAL");
    x3 = ent_data (e3);
    nr3 = x3->nrow;
    nc3 = x3->ncol;
    if (nr3 * nc3 == 0)
      rerror ("pdf.binomial: argument 'n' has to be MATRIX-DENSE-REAL");

    // update the size of output matrix to max of all three sizes
    nr = nr > nr3 ? nr : nr3;
    nc = nc > nc3 ? nc : nc3;
    w = mdr_Create (nr, nc);

    // calculate
    for (i = 0; i < nr; i++)
    {
      for (j = 0; j < nc; j++)
      {
        // x:
        id = MIN(i, nr1-1);
        jd = MIN(j, nc1-1);
        if (x1->type == RLAB_TYPE_INT32)
          x = (double) Mdi0 (x1, id, jd);
        else
          x = Mdr0 (x1, id, jd);
        // a:
        id = MIN(i, nr2-1);
        jd = MIN(j, nc2-1);
        if (x2->type == RLAB_TYPE_INT32)
          a  = (double) Mdi0 (x2, id, jd);
        else
          a  = Mdr0 (x2, id, jd);
        // b:
        id = MIN(i, nr3-1);
        jd = MIN(j, nc3-1);
        if (x3->type == RLAB_TYPE_INT32)
          b  = (double) Mdi0 (x3, id, jd);
        else
          b  = Mdr0 (x3, id, jd);

        // w:
        Mdr0 (w, i, j) = gsl_ran_binomial_pdf ( (int) x, a, (int) b);
      }
    }
  }

  // clean-up

  ent_Clean (e1);
  ent_Clean (e2);
  ent_Clean (e3);

  return ent_Assign_Rlab_MDR(w);
}

Ent *
ent_gsl_pdf_negbinomial (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *e3=0;
  int nr, nc, nr1, nc1, nr2, nc2, nr3, nc3, i, j, id, jd;
  double x, a, b;
  MDR *x1=0, *x2=0, *x3=0;
  MDR *w=0;

  if (nargs!=1 && nargs!=3)
  {
    printf ("pdf.negbinomial: Probability density for a binomial distribution\n");
    printf
        ("pdf.negbinomial: with parameters: 'p', probability of a single trial,\n");
    printf ("pdf.negbinomial: and 'n', number of independent trials.\n");
    printf ("pdf.negbinomial: Format:\n");
    printf ("pdf.negbinomial:   pdf.negbinomial(k,p,n)\n");
    printf ("pdf.negbinomial:   pdf.negbinomial( [k,p,n] )\n");
    rerror ("One, two or three arguments required");
  }

  //
  // get x
  //
  if (nargs == 1)
  {
    // get [x, param_1, param_2]
    e1 = bltin_get_ent (args[0]);
    if (ent_type (e1) != MATRIX_DENSE_REAL)
      rerror ("pdf.negbinomial: first argument has to be MATRIX-DENSE-REAL");
    x1 = ent_data (e1);
    nr1 = nr = x1->nrow;
    nc1 = nc = x1->ncol;
    if (nr * nc == 0)
      rerror ("pdf.negbinomial: first argument has to be MATRIX-DENSE-REAL");
    if (nc!=3)
      rerror ("pdf.negbinomial: Single argument implies [k,p,n] but it is not provided");

    w = mdr_Create (nr, 1);
    for (i = 0; i < nr; i++)
    {
      if (x1->type == RLAB_TYPE_INT32)
      {
        x = (double) Mdi0 (x1, i, 0);
        a = (double) Mdi0 (x1, i, 1);
        b = (double) Mdi0 (x1, i, 2);
      }
      else
      {
        x = Mdr0 (x1, i, 0);
        a = Mdr0 (x1, i, 1);
        b = Mdr0 (x1, i, 2);
      }
      MdrV0 (w, i) = gsl_ran_negative_binomial_pdf ( (int) x, a, (int) b);
    }
  }
  else if (nargs == 3)
  {
    // get x
    e1 = bltin_get_ent (args[0]);
    if (ent_type (e1) != MATRIX_DENSE_REAL)
      rerror ("pdf.negbinomial: argument 'k' has to be MATRIX-DENSE-REAL");
    x1 = ent_data (e1);
    nr1 = nr = x1->nrow;
    nc1 = nc = x1->ncol;

    //
    // get a
    //
    e2 = bltin_get_ent (args[1]);
    if (ent_type (e2) != MATRIX_DENSE_REAL)
      rerror ("pdf.negbinomial: argument 'p' has to be MATRIX-DENSE-REAL");
    x2 = ent_data (e2);
    nr2 = x2->nrow;
    nc2 = x2->ncol;
    if (nr2 * nc2 == 0)
      rerror ("pdf.negbinomial: argument 'p' has to be MATRIX-DENSE-REAL");

    // update the size of output matrix
    nr = nr > nr2 ? nr : nr2;
    nc = nc > nc2 ? nc : nc2;

    //
    // get b
    //
    e3 = bltin_get_ent (args[2]);
    if (ent_type (e3) != MATRIX_DENSE_REAL)
      rerror ("pdf.negbinomial: argument 'n' has to be MATRIX-DENSE-REAL");
    x3 = ent_data (e3);
    nr3 = x3->nrow;
    nc3 = x3->ncol;
    if (nr3 * nc3 == 0)
      rerror ("pdf.negbinomial: argument 'n' has to be MATRIX-DENSE-REAL");

    // update the size of output matrix to max of all three sizes
    nr = nr > nr3 ? nr : nr3;
    nc = nc > nc3 ? nc : nc3;
    w = mdr_Create (nr, nc);

    // calculate
    for (i = 0; i < nr; i++)
    {
      for (j = 0; j < nc; j++)
      {
        // x:
        id = MIN(i, nr1-1);
        jd = MIN(j, nc1-1);
        if (x1->type == RLAB_TYPE_INT32)
          x = (double) Mdi0 (x1, id, jd);
        else
          x = Mdr0 (x1, id, jd);
        // a:
        id = MIN(i, nr2-1);
        jd = MIN(j, nc2-1);
        if (x2->type == RLAB_TYPE_INT32)
          a  = (double) Mdi0 (x2, id, jd);
        else
          a  = Mdr0 (x2, id, jd);
        // b:
        id = MIN(i, nr3-1);
        jd = MIN(j, nc3-1);
        if (x3->type == RLAB_TYPE_INT32)
          b  = (double) Mdi0 (x3, id, jd);
        else
          b  = Mdr0 (x3, id, jd);

        // w:
        Mdr0 (w, i, j) = gsl_ran_negative_binomial_pdf ( (int) x, a, (int) b);
      }
    }
  }

  // clean-up

  ent_Clean (e1);
  ent_Clean (e2);
  ent_Clean (e3);

  return ent_Assign_Rlab_MDR(w);
}

Ent *
ent_gsl_pdf_hypergeom (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *e3=0, *e4=0;
  int nr=0, nc=0, nr1=0, nc1=0, nr2=0, nc2=0, nr3=0, nc3=0, nr4=0, nc4=0, i=0, j=0, id=0, jd=0;
  unsigned int k, n1, n2, t;
  MDR *x1=0, *x2=0, *x3=0, *x4=0;
  MDR *w=0;
  if (nargs!=1 && nargs!=4)
  {
    printf
      ("pdf.hypergeom: Probability density for a hypergeometric distribution\n");
    printf ("pdf.hypergeom: with parameters 'n1', 'n2' and 't'.\n");
    printf ("pdf.hypergeom: Format:\n");
    printf ("pdf.hypergeom:   pdf.hypergeom(k,n1,n2,t)\n");
    printf ("pdf.hypergeom:   pdf.hypergeom( [k,n1,n2,t] )\n");
    rerror ("No parameters given!");
  }

  //
  // [k,n1,n2,t]
  //
  if (nargs == 1)
  {
    // get [x, param_1, param_2]
    e1 = bltin_get_ent (args[0]);
    if (ent_type (e1) != MATRIX_DENSE_REAL)
      rerror ("pdf.hypergeom: first argument has to be MATRIX-DENSE-REAL");
    x1 = ent_data (e1);
    nr1 = nr = x1->nrow;
    nc1 = nc = x1->ncol;
    if (nr * nc == 0)
      rerror ("pdf.hypergeom: first argument has to be MATRIX-DENSE-REAL");
    if (nc!=4)
      rerror ("pdf.hypergeom: Single argument implies [k,n1,n2,t] but it is not provided");

    w = mdr_Create (nr, 1);
    for (i = 0; i < nr; i++)
    {
      if (x1->type == RLAB_TYPE_INT32)
      {
        k  = (unsigned int) Mdi0 (x1, i, 0);
        n1 = (unsigned int) Mdi0 (x1, i, 1);
        n2 = (unsigned int) Mdi0 (x1, i, 2);
        t  = (unsigned int) Mdi0 (x1, i, 3);
      }
      else
      {
        k  = (unsigned int) Mdr0 (x1, i, 0);
        n1 = (unsigned int) Mdr0 (x1, i, 1);
        n2 = (unsigned int) Mdr0 (x1, i, 2);
        t  = (unsigned int) Mdr0 (x1, i, 3);
      }
      MdrV0 (w, i) = gsl_ran_hypergeometric_pdf (k, n1, n2, t);
    }
  }
  else if (nargs == 4)
  {
    //
    // get n1
    //
    e2 = bltin_get_ent (args[1]);
    if (ent_type (e2) != MATRIX_DENSE_REAL)
      rerror ("pdf.hypergeom: argument 'n1' has to be MATRIX-DENSE-REAL");
    x2 = ent_data (e2);
    nr2 = x2->nrow;
    nc2 = x2->ncol;
    if (nr2 * nc2 == 0)
      rerror ("pdf.hypergeom: argument 'n1' has to be MATRIX-DENSE-REAL");
    // update the size of output matrix
    nr = nr > nr2 ? nr : nr2;
    nc = nc > nc2 ? nc : nc2;

    //
    // get n2
    //
    e3 = bltin_get_ent (args[2]);
    if (ent_type (e3) != MATRIX_DENSE_REAL)
      rerror ("pdf.hypergeom: argument 'n2' has to be MATRIX-DENSE-REAL");
    x3 = ent_data (e3);
    nr3 = x3->nrow;
    nc3 = x3->ncol;
    if (nr3 * nc3 == 0)
      rerror ("pdf.hypergeom: argument 'n' has to be MATRIX-DENSE-REAL");
    // update the size of output matrix
    nr = nr > nr3 ? nr : nr3;
    nc = nc > nc3 ? nc : nc3;

    //
    // get t
    //
    e4 = bltin_get_ent (args[3]);
    if (ent_type (e4) != MATRIX_DENSE_REAL)
      rerror ("pdf.hypergeom: argument 't' has to be MATRIX-DENSE-REAL");
    x4 = ent_data (e4);
    nr4 = x4->nrow;
    nc4 = x4->ncol;
    if (nr4 * nc4 == 0)
      rerror ("pdf.hypergeom: argument 't' has to be MATRIX-DENSE-REAL");
    // update the size of output matrix
    nr = nr > nr4 ? nr : nr4;
    nc = nc > nc4 ? nc : nc4;

    w = mdr_Create (nr, nc);

    // calculate
    // calculate
    for (i = 0; i < nr; i++)
    {
      for (j = 0; j < nc; j++)
      {
        // x:
        id = MIN(i, nr1-1);
        jd = MIN(j, nc1-1);
        if (x1->type == RLAB_TYPE_INT32)
          k = (unsigned int) Mdi0 (x1, id, jd);
        else
          k = (unsigned int) Mdr0 (x1, id, jd);
        // n1:
        id = MIN(i, nr2-1);
        jd = MIN(j, nc2-1);
        if (x2->type == RLAB_TYPE_INT32)
          n1  = (unsigned int) Mdi0 (x2, id, jd);
        else
          n1  = (unsigned int) Mdr0 (x2, id, jd);
        // n2:
        id = MIN(i, nr3-1);
        jd = MIN(j, nc3-1);
        if (x3->type == RLAB_TYPE_INT32)
          n2  = (unsigned int) Mdi0 (x3, id, jd);
        else
          n2  = (unsigned int) Mdr0 (x3, id, jd);
        // t:
        id = MIN(i, nr4-1);
        jd = MIN(j, nc4-1);
        if (x4->type == RLAB_TYPE_INT32)
          t  = (unsigned int) Mdi0 (x4, id, jd);
        else
          t  = (unsigned int) Mdr0 (x4, id, jd);

        // w:
        Mdr0 (w, i, j) = gsl_ran_hypergeometric_pdf (k, n1, n2, t);
      }
    }
  }

  // clean-up

  ent_Clean (e1);
  ent_Clean (e2);
  ent_Clean (e3);
  ent_Clean (e4);

  return ent_Assign_Rlab_MDR(w);
}

Ent *
ent_gsl_pdf_multinom (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0;
  int nr, nc, nr1, nc1, nr2, nc2, i, j, id;

  MDR *x1, *x2;
  MDR *w;
  if (nargs != 2)
  {
    printf
      ("pdf.multinom: Probability density for a multinomial distribution\n");
    printf ("pdf.multinom: with parameter arrays p=[p1..pK] and n=[n1..nK].\n");
    printf ("pdf.multinom: Format:\n");
    printf ("pdf.multinom:   pdf.multinom(p,n)\n");
    rerror ("No parameters given!");
  }

  //
  // get p
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) != MATRIX_DENSE_REAL)
    rerror ("pdf.multinom: argument 'p' has to be MATRIX-DENSE-REAL");
  x1 = ent_data (e1);
  nr = nr1 = x1->nrow;
  nc = nc1 = x1->ncol;
  if (nr1 * nc1 == 0)
    rerror ("pdf.multinom: argument 'p' has to be MATRIX-DENSE-REAL");

  //
  // get n
  //
  e2 = bltin_get_ent (args[1]);
  if (ent_type (e2) != MATRIX_DENSE_REAL)
    rerror ("pdf.multinom: argument 'n' has to be MATRIX-DENSE-REAL");
  x2 = ent_data (e2);
  nr2 = x2->nrow;
  nc2 = x2->ncol;
  if (nr2 * nc2 == 0)
    rerror ("pdf.multinom: argument 'n' has to be MATRIX-DENSE-REAL");

  // figure out the size of output matrix
  if (nc1 != nc2)
    rerror ("pdf.multinom: 'p' and 'n' have to have same number of columns");

  double *p = GC_malloc(nc * sizeof(double));
  unsigned int *n = GC_malloc(nc * sizeof(unsigned int));

  if (nr2 > nr)
    nr = nr2;
  w = mdr_Create (nr, 1);

  // calculate
  for (i=0; i<nr; i++)
  {
    // p:
    id = MIN(i, nr1-1);
    for (j = 0; j < nc; j++)
    {
      if (x1->type == RLAB_TYPE_INT32)
        p[j] = (double) Mdi0 (x1, id, j);
      else
        p[j] = Mdr0 (x1, id, j);
    }
    // n:

    id = MIN(i, nr2-1);
    for (j = 0; j < nc; j++)
    {
      if (x2->type == RLAB_TYPE_INT32)
        n[j] = (unsigned int) Mdi0 (x2, id, j);
      else
        n[j] = (unsigned int) Mdr0 (x2, id, j);
    }
    MdrV0 (w, i) = gsl_ran_multinomial_pdf (nc, p, (const unsigned int *)n);
  }

  GC_free(n);
  GC_free(p);

  // clean-up

  ent_Clean (e1);
  ent_Clean (e2);

  return ent_Assign_Rlab_MDR(w);
}

Ent *
ent_gsl_pdf_dirichlet (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0;
  int nr, nc, nr1, nc1, nr2, nc2, i, j, id;

  MDR *x1, *x2;
  MDR *w;
  if (nargs != 2)
  {
    printf ("pdf.dirichlet: Probability density for Dirichlet distribution\n");
    printf
        ("pdf.dirichlet: with parameter arrays alpha=[a1..aK] and theta=[t1..tK].\n");
    printf ("pdf.dirichlet: Format:\n");
    printf ("pdf.dirichlet:   pdf.dirichlet(alpha,theta)\n");
    rerror ("Two parameters required");
  }

  //
  // get p
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) != MATRIX_DENSE_REAL)
    rerror ("pdf.dirichlet: argument 'alpha' has to be MATRIX-DENSE-REAL");
  x1 = ent_data (e1);
  nr = nr1 = x1->nrow;
  nc = nc1 = x1->ncol;
  if (nr1 * nc1 == 0)
    rerror ("pdf.dirichlet: argument 'alpha' has to be MATRIX-DENSE-REAL");

  //
  // get n
  //
  e2 = bltin_get_ent (args[1]);
  if (ent_type (e2) != MATRIX_DENSE_REAL)
    rerror ("pdf.dirichlet: argument 'theta' has to be MATRIX-DENSE-REAL");
  x2 = ent_data (e2);
  nr2 = x2->nrow;
  nc2 = x2->ncol;
  if (nr2 * nc2 == 0)
    rerror ("pdf.dirichlet: argument 'theta' has to be MATRIX-DENSE-REAL");

  // figure out the size of output matrix
  if (nc1 != nc2)
    rerror ("pdf.dirichlet: 'alpha' and 'theta' have to have same number of columns");

  double *p = GC_malloc(nc * sizeof(double));
  double *t = GC_malloc(nc * sizeof(double));

  if (nr2 > nr)
    nr = nr2;
  w = mdr_Create (nr, 1);

  // calculate
  for (i=0; i<nr; i++)
  {
    // p:
    id = MIN(i, nr1-1);
    for (j = 0; j < nc; j++)
    {
      if (x1->type == RLAB_TYPE_INT32)
        p[j] = (double) Mdi0 (x1, id, j);
      else
        p[j] = Mdr0 (x1, id, j);
    }
    // t:
    id = MIN(i, nr2-1);
    for (j = 0; j < nc; j++)
    {
      if (x2->type == RLAB_TYPE_INT32)
        t[j] = (double) Mdi0 (x2, id, j);
      else
        t[j] = Mdr0 (x2, id, j);
    }
    MdrV0 (w, i) = gsl_ran_dirichlet_pdf (nc, p, t);
  }

  GC_free(t);
  GC_free(p);

  // clean-up

  ent_Clean (e1);
  ent_Clean (e2);

  return ent_Assign_Rlab_MDR(w);
}

//
// 5-parameter distributions
//
Ent *
ent_gsl_pdf_normalbiv (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *e3=0, *e4=0, *e5=0;
  int nr=0, nc=0, nr1=0, nc1=0, nr2=0, nc2=0, nr3=0, nc3=0, nc4=0, nr4=0, nc5=0, nr5=0, i=0, j=0, id=0, jd=0;
  double x, y, sx, sy, rho;
  MDR *x1=0, *x2=0, *x3=0, *x4=0, *x5=0;
  MDR *w=0;
  if (nargs!=1 && nargs!=2 && nargs!=5)
  {
    printf
      ("pdf.normalbiv: Probability density for a bivariate Gaussian distribution\n");
    printf
      ("pdf.normalbiv: with zero mean and standard deviations sigma_x, sigma_y,\n");
    printf ("pdf.normalbiv: and correlation coefficient 'rho'.\n");
    printf ("pdf.normalbiv: Format:\n");
    printf ("pdf.normalbiv:   pdf.normalbiv(x,y,sigma_x,sigma_y,rho)\n");
    printf ("pdf.normalbiv:   pdf.normalbiv( [x,y], [sigma_x,sigma_y,rho] )\n");
    printf ("pdf.normalbiv:   pdf.normalbiv( [x,y,sigma_x,sigma_y,rho] )\n");
    rerror ("No parameters given!");
  }

  if (nargs == 1)
  {
    // get [x,y,sigma_x,sigma_y,rho]
    e1 = bltin_get_ent (args[0]);
    if (ent_type (e1) != MATRIX_DENSE_REAL)
      rerror ("pdf.normalbiv: first argument has to be MATRIX-DENSE-REAL");
    x1 = ent_data (e1);
    nr1 = nr = x1->nrow;
    nc1 = nc = x1->ncol;
    if (nr * nc == 0)
      rerror ("pdf.normalbiv: first argument has to be MATRIX-DENSE-REAL");
    if (nc!=5)
      rerror ("pdf.normalbiv: Single argument implies [x,y,sigma_x,sigma_y,rho]"
          " but it is not provided");

    w = mdr_Create (nr, 1);
    for (i=0; i<nr; i++)
    {
      if (x1->type == RLAB_TYPE_INT32)
      {
        x   = (double) Mdi0 (x1, i, 0);
        y   = (double) Mdi0 (x1, i, 1);
        sx  = (double) Mdi0 (x1, i, 2);
        sy  = (double) Mdi0 (x1, i, 3);
        rho = (double) Mdi0 (x1, i, 4);
      }
      else
      {
        x   = Mdr0 (x1, i, 0);
        y   = Mdr0 (x1, i, 1);
        sx  = Mdr0 (x1, i, 2);
        sy  = Mdr0 (x1, i, 3);
        rho = Mdr0 (x1, i, 4);

      }

      // little adjustment of rho for silly users who did not read the manual
      if (rho >= 1)
        rho =  1.0 - 1e-10;
      if (rho <= -1)
        rho = -1.0 + 1e-10;

      MdrV0 (w, i) = gsl_ran_bivariate_gaussian_pdf (x, y, sx, sy, rho);
    }
  }
  else if (nargs == 2)
  {
    // get [x,y,sigma_x,sigma_y,rho]
    e1 = bltin_get_ent (args[0]);
    if (ent_type (e1) != MATRIX_DENSE_REAL)
      rerror ("pdf.normalbiv: first argument has to be MATRIX-DENSE-REAL");
    x1 = ent_data (e1);
    nr1 = nr = x1->nrow;
    nc1 = nc = x1->ncol;
    if (nr * nc == 0)
      rerror ("pdf.normalbiv: first argument has to be MATRIX-DENSE-REAL");
    if (nc!=2)
      rerror ("pdf.normalbiv: Two arguments imply [x,y]"
          " but it is not provided");

    // get [sigma_x,sigma_y,rho]
    e2 = bltin_get_ent (args[1]);
    if (ent_type (e2) != MATRIX_DENSE_REAL)
      rerror ("pdf.normalbiv: second argument has to be MATRIX-DENSE-REAL");
    x2 = ent_data (e2);
    nr2 = x2->nrow;
    nc2 = x2->ncol;
    if (nr2 * nc2 == 0)
      rerror ("pdf.normalbiv: second argument has to be MATRIX-DENSE-REAL");
    if (nc2!=3)
      rerror ("pdf.normalbiv: Two arguments imply  [sigma_x,sigma_y,rho]"
          " but it is not provided");

    // update size
    if (nr2 > nr)
      nr = nr2;

    w = mdr_Create (nr, 1);
    for (i=0; i<nr; i++)
    {
      // [x,y]
      if (x1->type == RLAB_TYPE_INT32)
      {
        x   = (double) Mdi0 (x1, i, 0);
        y   = (double) Mdi0 (x1, i, 1);
      }
      else
      {
        x   = Mdr0 (x1, i, 0);
        y   = Mdr0 (x1, i, 1);
      }
      // [sigma_x,sigma_y,rho]
      if (x2->type == RLAB_TYPE_INT32)
      {
        sx  = (double) Mdi0 (x2, i, 0);
        sy  = (double) Mdi0 (x2, i, 1);
        rho = (double) Mdi0 (x2, i, 2);
      }
      else
      {
        sx  = Mdr0 (x2, i, 0);
        sy  = Mdr0 (x2, i, 1);
        rho = Mdr0 (x2, i, 2);
      }

      // little adjustment of rho for silly users who did not read the manual
      if (rho >= 1)
        rho =  1.0 - 1e-10;
      if (rho <= -1)
        rho = -1.0 + 1e-10;

      MdrV0 (w, i) = gsl_ran_bivariate_gaussian_pdf (x, y, sx, sy, rho);
    }
  }
  else if (nargs == 5)
  {

    // x:
    e1 = bltin_get_ent (args[0]);
    if (ent_type (e1) != MATRIX_DENSE_REAL)
      rerror ("pdf.normalbiv: argument 'x' has to be MATRIX-DENSE-REAL");
    x1 = ent_data (e1);
    nr = nr1 = x1->nrow;
    nc = nc1 = x1->ncol;
    if (nr1 * nc1 == 0)
      rerror ("pdf.normalbiv: argument 'x' has to be MATRIX-DENSE-REAL!");

    // y:
    e2 = bltin_get_ent (args[1]);
    if (ent_type (e2) != MATRIX_DENSE_REAL)
      rerror ("pdf.normalbiv: argument 'y' has to be MATRIX-DENSE-REAL!");
    x2 = ent_data (e2);
    nr2 = x2->nrow;
    nc2 = x2->ncol;
    if (nr2 * nc2 == 0)
      rerror ("pdf.normalbiv: argument 'y' has to be MATRIX-DENSE-REAL");

    // figure out the size of output matrix
    if (nr < nr2)
      nr = nr2;
    if (nc < nc2)
      nc = nr2;

    // sigma_x:
    e3 = bltin_get_ent (args[2]);
    if (ent_type (e3) != MATRIX_DENSE_REAL)
      rerror ("pdf.normalbiv: argument 'sigma_x' has to be MATRIX-DENSE-REAL");
    x3 = ent_data (e3);
    nr3 = x3->nrow;
    nc3 = x3->ncol;
    if (nr3 * nc3 == 0)
      rerror ("pdf.normal argument 'sigma_x' has to be MATRIX-DENSE-REAL");
    // figure out the size of output matrix
    if (nc < nc3)
      nc = nc3;
    if (nr < nr3)
      nr = nr3;

    // sigma_y:
    e4 = bltin_get_ent (args[3]);
    if (ent_type (e4) != MATRIX_DENSE_REAL)
      rerror ("pdf.normalbiv: argument 'sigma_y' has to be MATRIX-DENSE-REAL");
    x4 = ent_data (e4);
    nr4 = x4->nrow;
    nc4 = x4->ncol;
    if (nr4 * nc4 == 0)
      rerror ("pdf.normalbiv: argument 'sigma_y' has to be MATRIX-DENSE-REAL");
    // figure out the size of output matrix
    if (nc < nc4)
      nc = nc4;
    if (nr < nr4)
      nr = nr4;

    // get rho
    e5 = bltin_get_ent (args[4]);
    if (ent_type (e5) != MATRIX_DENSE_REAL)
      rerror ("pdf.normalbiv: argument 'rho' has to be MATRIX-DENSE-REAL");
    x5 = ent_data (e5);
    nr5 = x5->nrow;
    nc5 = x5->ncol;
    if (nr5 * nc5 == 0)
      rerror ("pdf.normalbiv: argument 'rho' has to be MATRIX-DENSE-REAL");
    // figure out the size of output matrix
    if (nc < nc5)
      nc = nc5;
    if (nr < nr5)
      nr = nr5;

    w = mdr_Create (nr, nc);

    // calculate
    for (i=0; i<nr; i++)
    {
      for (j=0; j<nc; j++)
      {
        // x:
        id = MIN(i, nr1-1);
        jd = MIN(j, nc1-1);
        if (x1->type == RLAB_TYPE_INT32)
          x = (double) Mdi0 (x1, id, jd);
        else
          x = Mdr0 (x1, id, jd);

        // y:
        id = MIN(i, nr2-1);
        jd = MIN(j, nc2-1);
        if (x2->type == RLAB_TYPE_INT32)
          y  = (double) Mdi0 (x2, id, jd);
        else
          y  = Mdr0 (x2, id, jd);

        // sx:
        id = MIN(i, nr3-1);
        jd = MIN(j, nc3-1);
        if (x3->type == RLAB_TYPE_INT32)
          sx  = (double) Mdi0 (x3, id, jd);
        else
          sx  = Mdr0 (x3, id, jd);

        // sy:
        id = MIN(i, nr4-1);
        jd = MIN(j, nc4-1);
        if (x4->type == RLAB_TYPE_INT32)
          sy  = (double) Mdi0 (x4, id, jd);
        else
          sy  = Mdr0 (x4, id, jd);

        // rho:
        id = MIN(i, nr5-1);
        jd = MIN(j, nc5-1);
        if (x5->type == RLAB_TYPE_INT32)
          rho = (double) Mdi0 (x5, id, jd);
        else
          rho = Mdr0 (x5, id, jd);

        // little adjustment of rho for silly users who did not read the manual
        if (rho >= 1)
          rho =  1.0 - 1e-10;
        if (rho <= -1)
          rho = -1.0 + 1e-10;

        // do it!
        Mdr0 (w, i, j) = gsl_ran_bivariate_gaussian_pdf (x, y, sx, sy, rho);
      }
    }
  }

  // clean-up
  ent_Clean (e1);
  ent_Clean (e2);
  ent_Clean (e3);
  ent_Clean (e4);
  ent_Clean (e5);

  return ent_Assign_Rlab_MDR(w);
}
