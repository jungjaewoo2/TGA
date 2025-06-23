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

#include <gsl/gsl_cdf.h>


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
ent_gsl_cpf_logistic (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0;
  int nr, nc, nr1, nc1, nr2, nc2, i, j, id, jd;
  double a, x;
  MDR *x1=0, *x2=0;
  MDR *w=0;
  if (nargs!=1 && nargs!=2)
  {
    printf
      ("cpf.logistic: Probability density for a logistic "
       "distribution with scale parameter 'a'.\n");
    printf ("cpf.logistic: Format:\n");
    printf ("cpf.logistic:   cpf.logistic(x,a), or cpf.logistic( [x,a] )\n");
    rerror ("One or two arguments required");
  }

  //
  // get first argument and figure out if it is 'x' or '[x,a]'
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) != MATRIX_DENSE_REAL)
    rerror ("cpf.logistic: argument 'x' has to be MATRIX-DENSE-REAL");
  x1 = class_matrix_real (e1);
  nr1 = nr = x1->nrow;
  nc1 = nc = x1->ncol;
  if (nr1 * nc1 == 0)
    rerror ("cpf.logistic: argument 'x' has to be MATRIX-DENSE-REAL");
  if (nargs == 1)
  {
    if (nr*nc!=2 && nc!=2)
      rerror ("cpf.logistic: Single argument implies [x,a] but it is not provided");

    w = mdr_Create (nr, 1);
    for (i = 0; i < nr; i++)
    {
      x = Mdr0 (x1, i, 0);
      a = Mdr0 (x1, i, 1);
      MdrV0 (w, i) = gsl_cdf_logistic_P (x, a);
    }
  }
  else if (nargs == 2)
  {
    // get mu
    e2 = bltin_get_ent (args[1]);
    if (ent_type (e2) != MATRIX_DENSE_REAL)
      rerror ("cpf.logistic: argument 'a' has to be MATRIX-DENSE-REAL");
    nr1 = nr;
    nc1 = nc;
    x2 = class_matrix_real (e2);
    nr2 = x2->nrow;
    nc2 = x2->ncol;
    if (nr2 * nc2 == 0)
      rerror ("cpf.logistic: argument 'a' has to be MATRIX-DENSE-REAL");

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
        Mdr0 (w, i, j) = gsl_cdf_logistic_P (x, a);
      }
    }
  }

  ent_Clean (e1);
  ent_Clean (e2);

  return ent_Assign_Rlab_MDR(w);
}

Ent *
ent_gsl_cpf_t (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0;
  int nr, nc, nr1, nc1, nr2, nc2, i, j, id, jd;
  double a, x;
  MDR *x1=0, *x2=0;
  MDR *w=0;
  if (nargs!=1 && nargs!=2)
  {
    printf("cpf.t: Probability density for a t-distribution "
           "with 'nu' degrees of freedom.\n");
    printf ("cpf.t: Format:\n");
    printf ("cpf.t:   cpf.t(x,nu), or cpf.t( [x,nu] )\n");
    rerror ("One or two arguments required");
  }

  //
  // get first argument and figure out if it is 'x' or '[x,a]'
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) != MATRIX_DENSE_REAL)
    rerror ("cpf.t: argument 'x' has to be MATRIX-DENSE-REAL");
  x1 = class_matrix_real (e1);
  nr1 = nr = x1->nrow;
  nc1 = nc = x1->ncol;
  if (nr1 * nc1 == 0)
    rerror ("cpf.t: argument 'x' has to be MATRIX-DENSE-REAL");
  if (nargs == 1)
  {
    if (nr*nc!=2 && nc!=2)
      rerror ("cpf.t: Single argument implies [x,nu] but it is not provided");

    w = mdr_Create (nr, 1);
    for (i = 0; i < nr; i++)
    {
      x = Mdr0 (x1, i, 0);
      a = Mdr0 (x1, i, 1);
      MdrV0 (w, i) = gsl_cdf_tdist_P (x, a);
    }
  }
  else if (nargs == 2)
  {
    // get mu
    e2 = bltin_get_ent (args[1]);
    if (ent_type (e2) != MATRIX_DENSE_REAL)
      rerror ("cpf.t: argument 'nu' has to be MATRIX-DENSE-REAL");
    nr1 = nr;
    nc1 = nc;
    x2 = class_matrix_real (e2);
    nr2 = x2->nrow;
    nc2 = x2->ncol;
    if (nr2 * nc2 == 0)
      rerror ("cpf.t: argument 'nu' has to be MATRIX-DENSE-REAL");

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
        Mdr0 (w, i, j) = gsl_cdf_tdist_P (x, a);
      }
    }
  }

  ent_Clean (e1);
  ent_Clean (e2);

  return ent_Assign_Rlab_MDR(w);
}

Ent *
ent_gsl_cpf_chisq (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0;
  int nr, nc, nr1, nc1, nr2, nc2, i, j, id, jd;
  double a, x;
  MDR *x1=0, *x2=0;
  MDR *w=0;
  if (nargs!=1 && nargs!=2)
  {
    printf
        ("cpf.chisq: Probability density for a Chi-squared distribution with 'nu'\n");
    printf ("cpf.chisq: degrees of freedom.\n");
    printf ("cpf.chisq: Format:\n");
    printf ("cpf.chisq:   cpf.chisq(x,nu), or cpf.chisq( [x,nu] )\n");
    rerror ("One or two arguments required");
  }

  //
  // get first argument and figure out if it is 'x' or '[x,a]'
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) != MATRIX_DENSE_REAL)
    rerror ("cpf.chisq: argument 'x' has to be MATRIX-DENSE-REAL");
  x1 = class_matrix_real (e1);
  nr1 = nr = x1->nrow;
  nc1 = nc = x1->ncol;
  if (nr1 * nc1 == 0)
    rerror ("cpf.chisq: argument 'x' has to be MATRIX-DENSE-REAL");
  if (nargs == 1)
  {
    if (nr*nc!=2 && nc!=2)
      rerror ("cpf.chisq: Single argument implies [x,nu] but it is not provided");

    w = mdr_Create (nr, 1);
    for (i = 0; i < nr; i++)
    {
      x = Mdr0 (x1, i, 0);
      a = Mdr0 (x1, i, 1);
      MdrV0 (w, i) = gsl_cdf_chisq_P (x, a);
    }
  }
  else if (nargs == 2)
  {
    // get mu
    e2 = bltin_get_ent (args[1]);
    if (ent_type (e2) != MATRIX_DENSE_REAL)
      rerror ("cpf.chisq: argument 'nu' has to be MATRIX-DENSE-REAL");
    nr1 = nr;
    nc1 = nc;
    x2 = class_matrix_real (e2);
    nr2 = x2->nrow;
    nc2 = x2->ncol;
    if (nr2 * nc2 == 0)
      rerror ("cpf.chisq: argument 'nu' has to be MATRIX-DENSE-REAL");

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
        Mdr0 (w, i, j) = gsl_cdf_chisq_P (x, a);
      }
    }
  }
  ent_Clean (e1);
  ent_Clean (e2);

  return ent_Assign_Rlab_MDR(w);
}

Ent *
    ent_gsl_cpf_rayleigh (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0;
  int nr, nc, nr1, nc1, nr2, nc2, i, j, id, jd;
  double a, x;
  MDR *x1=0, *x2=0;
  MDR *w=0;
  if (nargs!=1 && nargs!=2)
  {
    printf
        ("cpf.rayleigh: Probability density for a Rayleigh distribution with scale\n");
    printf ("cpf.rayleigh: parameter 'sigma'.\n");
    printf ("cpf.rayleigh: Format:\n");
    printf ("cpf.rayleigh:   cpf.rayleigh(x,sigma), or cpf.rayleigh( [x,sigma] )\n");
    rerror ("One or two arguments required");
  }

  //
  // get first argument and figure out if it is 'x' or '[x,a]'
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) != MATRIX_DENSE_REAL)
    rerror ("cpf.rayleigh: argument 'x' has to be MATRIX-DENSE-REAL");
  x1 = class_matrix_real (e1);
  nr1 = nr = x1->nrow;
  nc1 = nc = x1->ncol;
  if (nr1 * nc1 == 0)
    rerror ("cpf.rayleigh: argument 'x' has to be MATRIX-DENSE-REAL");
  if (nargs == 1)
  {
    if (nr*nc!=2 && nc!=2)
      rerror ("cpf.rayleigh: Single argument implies [x,sigma] but it is not provided");

    w = mdr_Create (nr, 1);
    for (i = 0; i < nr; i++)
    {
      x = Mdr0 (x1, i, 0);
      a = Mdr0 (x1, i, 1);
      MdrV0 (w, i) = gsl_cdf_rayleigh_P (x, a);
    }
  }
  else if (nargs == 2)
  {
    // get mu
    e2 = bltin_get_ent (args[1]);
    if (ent_type (e2) != MATRIX_DENSE_REAL)
      rerror ("cpf.rayleigh: argument 'sigma' has to be MATRIX-DENSE-REAL");
    nr1 = nr;
    nc1 = nc;
    x2 = class_matrix_real (e2);
    nr2 = x2->nrow;
    nc2 = x2->ncol;
    if (nr2 * nc2 == 0)
      rerror ("cpf.rayleigh: argument 'sigma' has to be MATRIX-DENSE-REAL");

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
        Mdr0 (w, i, j) = gsl_cdf_rayleigh_P (x, a);
      }
    }
  }
  ent_Clean (e1);
  ent_Clean (e2);

  return ent_Assign_Rlab_MDR(w);
}

Ent *
    ent_gsl_cpf_cauchy (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0;
  int nr=0, nc=0, nr1=0, nc1=0, nr2=0, nc2=0, i=0, j=0, id=0, jd=0;
  double a, x;
  MDR *x1=0, *x2=0;
  MDR *w=0;
  if (nargs!=1 && nargs!=2)
  {
    printf
        ("cpf.cauchy: Probability density for a Cauchy distribution with scale\n");
    printf ("cpf.cauchy: parameter 'a'.\n");
    printf ("cpf.cauchy: Format:\n");
    printf ("cpf.cauchy:   cpf.cauchy(x,a), or cpf.cauchy( [x,a] )\n");
    rerror ("One or two arguments required");
  }

  //
  // get first argument and figure out if it is 'x' or '[x,a]'
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) != MATRIX_DENSE_REAL)
    rerror ("cpf.cauchy: argument 'x' has to be MATRIX-DENSE-REAL");
  x1 = class_matrix_real (e1);
  nr1 = nr = x1->nrow;
  nc1 = nc = x1->ncol;
  if (nr1 * nc1 == 0)
    rerror ("cpf.cauchy: argument 'x' has to be MATRIX-DENSE-REAL");
  if (nargs == 1)
  {
    if (nr*nc!=2 && nc!=2)
      rerror ("cpf.cauchy: Single argument implies [x,a] but it is not provided");

    w = mdr_Create (nr, 1);
    for (i = 0; i < nr; i++)
    {
      x = Mdr0 (x1, i, 0);
      a = Mdr0 (x1, i, 1);
      MdrV0 (w, i) = gsl_cdf_cauchy_P (x, a);
    }
  }
  else if (nargs == 2)
  {
    // get mu
    e2 = bltin_get_ent (args[1]);
    if (ent_type (e2) != MATRIX_DENSE_REAL)
      rerror ("cpf.cauchy: argument 'a' has to be MATRIX-DENSE-REAL");
    nr1 = nr;
    nc1 = nc;
    x2 = class_matrix_real (e2);
    nr2 = x2->nrow;
    nc2 = x2->ncol;
    if (nr2 * nc2 == 0)
      rerror ("cpf.cauchy: argument 'a' has to be MATRIX-DENSE-REAL");

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
        Mdr0 (w, i, j) = gsl_cdf_cauchy_P (x, a);
      }
    }
  }
  ent_Clean (e1);
  ent_Clean (e2);

  return ent_Assign_Rlab_MDR(w);
}

Ent *
ent_gsl_cpf_laplace (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0;
  int nr, nc, nr1, nc1, nr2, nc2, i, j, id, jd;
  double a, x;
  MDR *x1=0, *x2=0;
  MDR *w=0;
  if (nargs!=1 && nargs!=2)
  {
    printf
        ("cpf.laplace: Probability density for a Laplace distribution with mean 'a'.\n");
    printf
        ("cpf.laplace: Format:\n");
    printf
        ("cpf.laplace:   cpf.laplace(x,a), or cpf.laplace( [x,a] )\n");
    rerror ("One or two arguments required");
  }

  //
  // get first argument and figure out if it is 'x' or '[x,a]'
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) != MATRIX_DENSE_REAL)
    rerror ("cpf.laplace: argument 'x' has to be MATRIX-DENSE-REAL");
  x1 = class_matrix_real (e1);
  nr1 = nr = x1->nrow;
  nc1 = nc = x1->ncol;
  if (nr1 * nc1 == 0)
    rerror ("cpf.laplace: argument 'x' has to be MATRIX-DENSE-REAL");
  if (nargs == 1)
  {
    if (nr*nc!=2 && nc!=2)
      rerror ("cpf.laplace: Single argument implies [x,a] but it is not provided");

    w = mdr_Create (nr, 1);
    for (i = 0; i < nr; i++)
    {
      x = Mdr0 (x1, i, 0);
      a = Mdr0 (x1, i, 1);
      MdrV0 (w, i) = gsl_cdf_laplace_P (x, a);
    }
  }
  else if (nargs == 2)
  {
    // get mu
    e2 = bltin_get_ent (args[1]);
    if (ent_type (e2) != MATRIX_DENSE_REAL)
      rerror ("cpf.laplace: argument 'a' has to be MATRIX-DENSE-REAL");
    nr1 = nr;
    nc1 = nc;
    x2 = class_matrix_real (e2);
    nr2 = x2->nrow;
    nc2 = x2->ncol;
    if (nr2 * nc2 == 0)
      rerror ("cpf.laplace: argument 'a' has to be MATRIX-DENSE-REAL");

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
        Mdr0 (w, i, j) = gsl_cdf_laplace_P (x, a);
      }
    }
  }
  ent_Clean (e1);
  ent_Clean (e2);

  return ent_Assign_Rlab_MDR(w);
}

Ent *
ent_gsl_cpf_exp (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0;
  int nr, nc, nr1, nc1, nr2, nc2, i, j, id, jd;
  double a, x;
  MDR *x1=0, *x2=0;
  MDR *w=0;
  if (nargs!=1 && nargs!=2)
  {
    printf
        ("cpf.exp: Probability density for exponential distribution with mean mu.\n");
    printf
        ("cpf.exp: Format:\n");
    printf
        ("cpf.exp:   cpf.exp(x,mu), or cpf.exp( [x,mu] )\n");
    rerror ("One or two arguments required");
  }

  //
  // get first argument and figure out if it is 'x' or '[x,a]'
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) != MATRIX_DENSE_REAL)
    rerror ("cpf.exp: argument 'x' has to be MATRIX-DENSE-REAL");
  x1 = class_matrix_real (e1);
  nr1 = nr = x1->nrow;
  nc1 = nc = x1->ncol;
  if (nr * nc == 0)
    rerror ("cpf.exp: argument 'x' has to be MATRIX-DENSE-REAL");
  if (nargs == 1)
  {
    if (nr*nc!=2 && nc!=2)
      rerror ("cpf.exp: Single argument implies [x,a] but it is not provided");

    w = mdr_Create (nr, 1);
    for (i = 0; i < nr; i++)
    {
      x = Mdr0 (x1, i, 0);
      a = Mdr0 (x1, i, 1);
      MdrV0 (w, i) = gsl_cdf_exponential_P (x, a);
    }
  }
  else if (nargs == 2)
  {
    // get mu
    e2 = bltin_get_ent (args[1]);
    if (ent_type (e2) != MATRIX_DENSE_REAL)
      rerror ("cpf.exp: argument 'a' has to be MATRIX-DENSE-REAL");
    nr1 = nr;
    nc1 = nc;
    x2 = class_matrix_real (e2);
    nr2 = x2->nrow;
    nc2 = x2->ncol;
    if (nr2 * nc2 == 0)
      rerror ("cpf.exp: argument 'a' has to be MATRIX-DENSE-REAL");

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
        Mdr0 (w, i, j) = gsl_cdf_exponential_P (x, a);
      }
    }
  }
  ent_Clean (e1);
  ent_Clean (e2);

  return ent_Assign_Rlab_MDR(w);
}

Ent *
ent_gsl_cpf_poisson (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0;
  int nr, nc, nr1, nc1, nr2, nc2, i, j, id, jd;
  double a, x;
  MDR *x1=0, *x2=0;
  MDR *w=0;
  if (nargs!=1 && nargs!=2)
  {
    printf
        ("cpf.poisson: Probability density for a Poisson distribution with mean 'mu'.\n");
    printf ("cpf.poisson: Format:\n");
    printf ("cpf.poisson:   cpf.poisson(k,mu), or cpf.poisson( [k,mu] )\n");
    rerror ("One or two arguments required");
  }

  //
  // get first argument and figure out if it is 'x' or '[x,a]'
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) != MATRIX_DENSE_REAL)
    rerror ("cpf.poisson: argument 'k' has to be MATRIX-DENSE-REAL");
  x1 = class_matrix_real (e1);
  nr1 = nr = x1->nrow;
  nc1 = nc = x1->ncol;
  if (nr1 * nc1 == 0)
    rerror ("cpf.poisson: argument 'k' has to be MATRIX-DENSE-REAL");
  if (nargs == 1)
  {
    if (nr*nc!=2 && nc!=2)
      rerror ("cpf.poisson: Single argument implies [k,mu] but it is not provided");

    w = mdr_Create (nr, 1);
    for (i = 0; i < nr; i++)
    {
      x = Mdr0 (x1, i, 0);
      a = Mdr0 (x1, i, 1);
      MdrV0 (w, i) = gsl_cdf_poisson_P ((int) x, a);
    }
  }
  else if (nargs == 2)
  {
    // get mu
    e2 = bltin_get_ent (args[1]);
    if (ent_type (e2) != MATRIX_DENSE_REAL)
      rerror ("cpf.poisson: argument 'mu' has to be MATRIX-DENSE-REAL");
    nr1 = nr;
    nc1 = nc;
    x2 = class_matrix_real (e2);
    nr2 = x2->nrow;
    nc2 = x2->ncol;
    if (nr2 * nc2 == 0)
      rerror ("cpf.poisson: argument 'mu' has to be MATRIX-DENSE-REAL");

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
        Mdr0 (w, i, j) = gsl_cdf_poisson_P ((int) x, a);
      }
    }
  }
  ent_Clean (e1);
  ent_Clean (e2);

  return ent_Assign_Rlab_MDR(w);
}

Ent *
ent_gsl_cpf_geometric (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0;
  int nr, nc, nr1, nc1, nr2, nc2, i, j, id, jd;
  double a, x;
  MDR *x1=0, *x2=0;
  MDR *w=0;
  if (nargs!=1 && nargs!=2)
  {
    printf
        ("cpf.geometric: Probability density for geometric distributions with\n");
    printf ("cpf.geometric: probability parameter 'p'.\n");
    printf ("cpf.geometric: Format:\n");
    printf ("cpf.geometric:   cpf.geometric(k,p), or cpf.geometric( [k,p] )\n");
    rerror ("One or two arguments required");
  }

  //
  // get first argument and figure out if it is 'x' or '[x,a]'
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) != MATRIX_DENSE_REAL)
    rerror ("cpf.geometric: argument 'k' has to be MATRIX-DENSE-REAL");
  x1 = class_matrix_real (e1);
  nr1 = nr = x1->nrow;
  nc1 = nc = x1->ncol;
  if (nr1 * nc1 == 0)
    rerror ("cpf.geometric: argument 'k' has to be MATRIX-DENSE-REAL");
  if (nargs == 1)
  {
    if (nr*nc!=2 && nc!=2)
      rerror ("cpf.geometric: Single argument implies [k,p] but it is not provided");

    w = mdr_Create (nr, 1);
    for (i = 0; i < nr; i++)
    {
      x = Mdr0 (x1, i, 0);
      a = Mdr0 (x1, i, 1);
      MdrV0 (w, i) = gsl_cdf_geometric_P ((int) x, a);
    }
  }
  else if (nargs == 2)
  {
    // get mu
    e2 = bltin_get_ent (args[1]);
    if (ent_type (e2) != MATRIX_DENSE_REAL)
      rerror ("cpf.geometric: argument 'p' has to be MATRIX-DENSE-REAL");
    nr1 = nr;
    nc1 = nc;
    x2 = class_matrix_real (e2);
    nr2 = x2->nrow;
    nc2 = x2->ncol;
    if (nr2 * nc2 == 0)
      rerror ("cpf.geometric: argument 'p' has to be MATRIX-DENSE-REAL");

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
        Mdr0 (w, i, j) = gsl_cdf_geometric_P ((int) x, a);
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
// GSL distributions
// ----------------------------------------------------------

Ent *
ent_gsl_cpf_normal (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *e3=0;
  int nr, nc, nr1, nc1, nr2, nc2, nr3, nc3, i, j, id, jd;
  double x, a, b;
  MDR *x1=0, *x2=0, *x3=0;
  MDR *w=0;

  if (nargs!=1 && nargs!=2 && nargs!=3)
  {
    printf ("cpf.normal: Probability density for a Gaussian distribution with mean,\n");
    printf ("cpf.normal: and standard deviation.\n");
    printf ("cpf.normal: Format:\n");
    printf ("cpf.normal:   cpf.normal(x, mean, std), or\n");
    printf ("cpf.normal:   cpf.normal(x, [mean,std] ), or\n");
    printf ("cpf.normal:   cpf.normal( [x,mean,std] ) \n");
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
      rerror ("cpf.normal: argument 'x' has to be MATRIX-DENSE-REAL");
    x1 = ent_data (e1);
    nr1 = nr = x1->nrow;
    nc1 = nc = x1->ncol;
    if (nr * nc == 0)
      rerror ("cpf.normal: argument 'x' has to be MATRIX-DENSE-REAL");
    if (nc!=3)
      rerror ("cpf.normal: Single argument implies [x,mean,std] but it is not provided");

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
      MdrV0 (w, i)  = gsl_cdf_gaussian_P (x-a, b);
    }
  }
  else if (nargs == 2)
  {
    // get x
    e1 = bltin_get_ent (args[0]);
    if (ent_type (e1) != MATRIX_DENSE_REAL)
      rerror ("cpf.normal: argument 'x' has to be MATRIX-DENSE-REAL");
    x1 = ent_data (e1);
    nr1 = nr = x1->nrow;
    nc1 = nc = x1->ncol;
    if (nr * nc == 0)
      rerror ("cpf.normal: argument 'x' has to be MATRIX-DENSE-REAL");

    // get [param_1, param_2]
    e2 = bltin_get_ent (args[1]);
    if (ent_type (e2) != MATRIX_DENSE_REAL)
      rerror ("cpf.normal: argument '[mean,std]' has to be MATRIX-DENSE-REAL");
    x2 = ent_data (e2);
    nr2 = x2->nrow;
    nc2 = x2->ncol;
    if (nc2 != 2)
      rerror ("cpf.normal: second argument has to be [mean,std]");

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
        Mdr0 (w, i, j) = gsl_cdf_gaussian_P (x-a, b);
      }
    }
  }
  else if (nargs == 3)
  {
    // get x
    e1 = bltin_get_ent (args[0]);
    if (ent_type (e1) != MATRIX_DENSE_REAL)
      rerror ("cpf.normal: argument 'x' has to be MATRIX-DENSE-REAL");
    x1 = ent_data (e1);
    nr1 = nr = x1->nrow;
    nc1 = nc = x1->ncol;

    //
    // get a
    //
    e2 = bltin_get_ent (args[1]);
    if (ent_type (e2) != MATRIX_DENSE_REAL)
      rerror ("cpf.normal: argument 'mean' has to be MATRIX-DENSE-REAL");
    x2 = ent_data (e2);
    nr2 = x2->nrow;
    nc2 = x2->ncol;
    if (nr2 * nc2 == 0)
      rerror ("cpf.normal: argument 'mean' has to be MATRIX-DENSE-REAL");

    // update the size of output matrix
    nr = nr > nr2 ? nr : nr2;
    nc = nc > nc2 ? nc : nc2;

    //
    // get b
    //
    e3 = bltin_get_ent (args[2]);
    if (ent_type (e3) != MATRIX_DENSE_REAL)
      rerror ("cpf.normal: argument 'std' has to be MATRIX-DENSE-REAL");
    x3 = ent_data (e3);
    nr3 = x3->nrow;
    nc3 = x3->ncol;
    if (nr3 * nc3 == 0)
      rerror ("cpf.normal: argument 'std' has to be MATRIX-DENSE-REAL");

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
        Mdr0 (w, i, j) = gsl_cdf_gaussian_P (x-a, b);
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
ent_gsl_cpf_uniform (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *e3=0;
  int nr, nc, nr1, nc1, nr2, nc2, nr3, nc3, i, j, id, jd;
  double x, a, b;
  MDR *x1=0, *x2=0, *x3=0;
  MDR *w=0;

  if (nargs!=1 && nargs!=2 && nargs!=3)
  {
    printf
        ("cpf.uniform: Probability density for uniform distribution\n");
    printf ("cpf.uniform: with parameters 'a' and 'b'.\n");
    printf ("cpf.uniform: Format:\n");
    printf ("cpf.uniform:   cpf.uniform(x,a,b)\n");
    printf ("cpf.uniform:   cpf.uniform(x, [a,b] ), or\n");
    printf ("cpf.uniform:   cpf.uniform( [x,a,b] ) \n");
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
      rerror ("cpf.uniform: argument 'x' has to be MATRIX-DENSE-REAL");
    x1 = ent_data (e1);
    nr1 = nr = x1->nrow;
    nc1 = nc = x1->ncol;
    if (nr * nc == 0)
      rerror ("cpf.uniform: argument 'x' has to be MATRIX-DENSE-REAL");
    if (nc!=3)
      rerror ("cpf.uniform: Single argument implies [x,a,b] but it is not provided");

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
      MdrV0 (w, i)  = gsl_cdf_flat_P (x, a, b);
    }
  }
  else if (nargs == 2)
  {
    // get x
    e1 = bltin_get_ent (args[0]);
    if (ent_type (e1) != MATRIX_DENSE_REAL)
      rerror ("cpf.uniform: argument 'x' has to be MATRIX-DENSE-REAL");
    x1 = ent_data (e1);
    nr1 = nr = x1->nrow;
    nc1 = nc = x1->ncol;
    if (nr * nc == 0)
      rerror ("cpf.uniform: argument 'x' has to be MATRIX-DENSE-REAL");

    // get [param_1, param_2]
    e2 = bltin_get_ent (args[1]);
    if (ent_type (e2) != MATRIX_DENSE_REAL)
      rerror ("cpf.uniform: argument '[a,b]' has to be MATRIX-DENSE-REAL");
    x2 = ent_data (e2);
    nr2 = x2->nrow;
    nc2 = x2->ncol;
    if (nc2 != 2)
      rerror ("cpf.uniform: second argument has to be [a,b]");

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
        Mdr0 (w, i, j) = gsl_cdf_flat_P (x, a, b);
      }
    }
  }
  else if (nargs == 3)
  {
    // get x
    e1 = bltin_get_ent (args[0]);
    if (ent_type (e1) != MATRIX_DENSE_REAL)
      rerror ("cpf.uniform: argument 'x' has to be MATRIX-DENSE-REAL");
    x1 = ent_data (e1);
    nr1 = nr = x1->nrow;
    nc1 = nc = x1->ncol;

    //
    // get a
    //
    e2 = bltin_get_ent (args[1]);
    if (ent_type (e2) != MATRIX_DENSE_REAL)
      rerror ("cpf.uniform: argument 'a' has to be MATRIX-DENSE-REAL");
    x2 = ent_data (e2);
    nr2 = x2->nrow;
    nc2 = x2->ncol;
    if (nr2 * nc2 == 0)
      rerror ("cpf.uniform: argument 'a' has to be MATRIX-DENSE-REAL");

    // update the size of output matrix
    nr = nr > nr2 ? nr : nr2;
    nc = nc > nc2 ? nc : nc2;

    //
    // get b
    //
    e3 = bltin_get_ent (args[2]);
    if (ent_type (e3) != MATRIX_DENSE_REAL)
      rerror ("cpf.uniform: argument 'b' has to be MATRIX-DENSE-REAL");
    x3 = ent_data (e3);
    nr3 = x3->nrow;
    nc3 = x3->ncol;
    if (nr3 * nc3 == 0)
      rerror ("cpf.uniform: argument 'b' has to be MATRIX-DENSE-REAL");

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
        Mdr0 (w, i, j) = gsl_cdf_flat_P (x, a, b);
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
ent_gsl_cpf_gumbel2 (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *e3=0;
  int nr, nc, nr1, nc1, nr2, nc2, nr3, nc3, i, j, id, jd;
  double x, a, b;
  MDR *x1=0, *x2=0, *x3=0;
  MDR *w=0;

  if (nargs!=1 && nargs!=2 && nargs!=3)
  {
    printf
        ("cpf.gumbel2: Probability density for a Type-2 Gumbel distribution\n");
    printf ("cpf.gumbel2: with parameters 'a' and 'b'.\n");
    printf ("cpf.gumbel2: Format:\n");
    printf ("cpf.gumbel2:   cpf.gumbel2(x,a,b)\n");
    printf ("cpf.gumbel2:   cpf.gumbel2(x, [a,b] ), or\n");
    printf ("cpf.gumbel2:   cpf.gumbel2( [x,a,b] ) \n");
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
      rerror ("cpf.gumbel2: argument 'x' has to be MATRIX-DENSE-REAL");
    x1 = ent_data (e1);
    nr1 = nr = x1->nrow;
    nc1 = nc = x1->ncol;
    if (nr * nc == 0)
      rerror ("cpf.gumbel2: argument 'x' has to be MATRIX-DENSE-REAL");
    if (nc!=3)
      rerror ("cpf.gumbel2: Single argument implies [x,a,b] but it is not provided");

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
      MdrV0 (w, i)  = gsl_cdf_gumbel2_P (x, a, b);
    }
  }
  else if (nargs == 2)
  {
    // get x
    e1 = bltin_get_ent (args[0]);
    if (ent_type (e1) != MATRIX_DENSE_REAL)
      rerror ("cpf.gumbel2: argument 'x' has to be MATRIX-DENSE-REAL");
    x1 = ent_data (e1);
    nr1 = nr = x1->nrow;
    nc1 = nc = x1->ncol;
    if (nr * nc == 0)
      rerror ("cpf.gumbel2: argument 'x' has to be MATRIX-DENSE-REAL");

    // get [param_1, param_2]
    e2 = bltin_get_ent (args[1]);
    if (ent_type (e2) != MATRIX_DENSE_REAL)
      rerror ("cpf.gumbel2: argument '[a,b]' has to be MATRIX-DENSE-REAL");
    x2 = ent_data (e2);
    nr2 = x2->nrow;
    nc2 = x2->ncol;
    if (nc2 != 2)
      rerror ("cpf.gumbel2: second argument has to be [a,b]");

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
        Mdr0 (w, i, j) = gsl_cdf_gumbel2_P (x, a, b);
      }
    }
  }
  else if (nargs == 3)
  {
    // get x
    e1 = bltin_get_ent (args[0]);
    if (ent_type (e1) != MATRIX_DENSE_REAL)
      rerror ("cpf.gumbel2: argument 'x' has to be MATRIX-DENSE-REAL");
    x1 = ent_data (e1);
    nr1 = nr = x1->nrow;
    nc1 = nc = x1->ncol;

    //
    // get a
    //
    e2 = bltin_get_ent (args[1]);
    if (ent_type (e2) != MATRIX_DENSE_REAL)
      rerror ("cpf.gumbel2: argument 'a' has to be MATRIX-DENSE-REAL");
    x2 = ent_data (e2);
    nr2 = x2->nrow;
    nc2 = x2->ncol;
    if (nr2 * nc2 == 0)
      rerror ("cpf.gumbel2: argument 'a' has to be MATRIX-DENSE-REAL");

    // update the size of output matrix
    nr = nr > nr2 ? nr : nr2;
    nc = nc > nc2 ? nc : nc2;

    //
    // get b
    //
    e3 = bltin_get_ent (args[2]);
    if (ent_type (e3) != MATRIX_DENSE_REAL)
      rerror ("cpf.gumbel2: argument 'b' has to be MATRIX-DENSE-REAL");
    x3 = ent_data (e3);
    nr3 = x3->nrow;
    nc3 = x3->ncol;
    if (nr3 * nc3 == 0)
      rerror ("cpf.gumbel2: argument 'b' has to be MATRIX-DENSE-REAL");

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
        Mdr0 (w, i, j) = gsl_cdf_gumbel2_P (x, a, b);
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
ent_gsl_cpf_gumbel1 (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *e3=0;
  int nr, nc, nr1, nc1, nr2, nc2, nr3, nc3, i, j, id, jd;
  double x, a, b;
  MDR *x1=0, *x2=0, *x3=0;
  MDR *w=0;

  if (nargs!=1 && nargs!=2 && nargs!=3)
  {
    printf
        ("cpf.gumbel1: Probability density for a Type-1 Gumbel distribution\n");
    printf ("cpf.gumbel1: with parameters 'a' and 'b'.\n");
    printf ("cpf.gumbel1: Format:\n");
    printf ("cpf.gumbel1:   cpf.gumbel1(x,a,b)\n");
    printf ("cpf.gumbel1:   cpf.gumbel1(x, [a,b] ), or\n");
    printf ("cpf.gumbel1:   cpf.gumbel1( [x,a,b] ) \n");
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
      rerror ("cpf.gumbel1: argument 'x' has to be MATRIX-DENSE-REAL");
    x1 = ent_data (e1);
    nr1 = nr = x1->nrow;
    nc1 = nc = x1->ncol;
    if (nr * nc == 0)
      rerror ("cpf.gumbel1: argument 'x' has to be MATRIX-DENSE-REAL");
    if (nc!=3)
      rerror ("cpf.gumbel1: Single argument implies [x,a,b] but it is not provided");

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
      MdrV0 (w, i)  = gsl_cdf_gumbel1_P (x, a, b);
    }
  }
  else if (nargs == 2)
  {
    // get x
    e1 = bltin_get_ent (args[0]);
    if (ent_type (e1) != MATRIX_DENSE_REAL)
      rerror ("cpf.gumbel1: argument 'x' has to be MATRIX-DENSE-REAL");
    x1 = ent_data (e1);
    nr1 = nr = x1->nrow;
    nc1 = nc = x1->ncol;
    if (nr * nc == 0)
      rerror ("cpf.gumbel1: argument 'x' has to be MATRIX-DENSE-REAL");

    // get [param_1, param_2]
    e2 = bltin_get_ent (args[1]);
    if (ent_type (e2) != MATRIX_DENSE_REAL)
      rerror ("cpf.gumbel1: argument '[a,b]' has to be MATRIX-DENSE-REAL");
    x2 = ent_data (e2);
    nr2 = x2->nrow;
    nc2 = x2->ncol;
    if (nc2 != 2)
      rerror ("cpf.gumbel1: second argument has to be [a,b]");

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
        Mdr0 (w, i, j) = gsl_cdf_gumbel1_P (x, a, b);
      }
    }
  }
  else if (nargs == 3)
  {
    // get x
    e1 = bltin_get_ent (args[0]);
    if (ent_type (e1) != MATRIX_DENSE_REAL)
      rerror ("cpf.gumbel1: argument 'x' has to be MATRIX-DENSE-REAL");
    x1 = ent_data (e1);
    nr1 = nr = x1->nrow;
    nc1 = nc = x1->ncol;

    //
    // get a
    //
    e2 = bltin_get_ent (args[1]);
    if (ent_type (e2) != MATRIX_DENSE_REAL)
      rerror ("cpf.gumbel1: argument 'a' has to be MATRIX-DENSE-REAL");
    x2 = ent_data (e2);
    nr2 = x2->nrow;
    nc2 = x2->ncol;
    if (nr2 * nc2 == 0)
      rerror ("cpf.gumbel1: argument 'a' has to be MATRIX-DENSE-REAL");

    // update the size of output matrix
    nr = nr > nr2 ? nr : nr2;
    nc = nc > nc2 ? nc : nc2;

    //
    // get b
    //
    e3 = bltin_get_ent (args[2]);
    if (ent_type (e3) != MATRIX_DENSE_REAL)
      rerror ("cpf.gumbel1: argument 'b' has to be MATRIX-DENSE-REAL");
    x3 = ent_data (e3);
    nr3 = x3->nrow;
    nc3 = x3->ncol;
    if (nr3 * nc3 == 0)
      rerror ("cpf.gumbel1: argument 'b' has to be MATRIX-DENSE-REAL");

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
        Mdr0 (w, i, j) = gsl_cdf_gumbel1_P (x, a, b);
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
ent_gsl_cpf_weibull (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *e3=0;
  int nr, nc, nr1, nc1, nr2, nc2, nr3, nc3, i, j, id, jd;
  double x, a, b;
  MDR *x1=0, *x2=0, *x3=0;
  MDR *w=0;

  if (nargs!=1 && nargs!=2 && nargs!=3)
  {
    printf ("cpf.weibull: Probability density for a Weibull distribution\n");
    printf ("cpf.weibull: with parameters 'a' and 'b'.\n");
    printf ("cpf.weibull: Format:\n");
    printf ("cpf.weibull:   cpf.weibull(x,a,b)\n");
    printf ("cpf.weibull:   cpf.weibull(x, [a,b] ), or\n");
    printf ("cpf.weibull:   cpf.weibull( [x,a,b] ) \n");
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
      rerror ("cpf.weibull: argument 'x' has to be MATRIX-DENSE-REAL");
    x1 = ent_data (e1);
    nr1 = nr = x1->nrow;
    nc1 = nc = x1->ncol;
    if (nr * nc == 0)
      rerror ("cpf.weibull: argument 'x' has to be MATRIX-DENSE-REAL");
    if (nc!=3)
      rerror ("cpf.weibull: Single argument implies [x,a,b] but it is not provided");

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
      MdrV0 (w, i)  = gsl_cdf_weibull_P (x, a, b);
    }
  }
  else if (nargs == 2)
  {
    // get x
    e1 = bltin_get_ent (args[0]);
    if (ent_type (e1) != MATRIX_DENSE_REAL)
      rerror ("cpf.weibull: argument 'x' has to be MATRIX-DENSE-REAL");
    x1 = ent_data (e1);
    nr1 = nr = x1->nrow;
    nc1 = nc = x1->ncol;
    if (nr * nc == 0)
      rerror ("cpf.weibull: argument 'x' has to be MATRIX-DENSE-REAL");

    // get [param_1, param_2]
    e2 = bltin_get_ent (args[1]);
    if (ent_type (e2) != MATRIX_DENSE_REAL)
      rerror ("cpf.weibull: argument '[a,b]' has to be MATRIX-DENSE-REAL");
    x2 = ent_data (e2);
    nr2 = x2->nrow;
    nc2 = x2->ncol;
    if (nc2 != 2)
      rerror ("cpf.weibull: second argument has to be [a,b]");

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
        Mdr0 (w, i, j) = gsl_cdf_weibull_P (x, a, b);
      }
    }
  }
  else if (nargs == 3)
  {
    // get x
    e1 = bltin_get_ent (args[0]);
    if (ent_type (e1) != MATRIX_DENSE_REAL)
      rerror ("cpf.weibull: argument 'x' has to be MATRIX-DENSE-REAL");
    x1 = ent_data (e1);
    nr1 = nr = x1->nrow;
    nc1 = nc = x1->ncol;

    //
    // get a
    //
    e2 = bltin_get_ent (args[1]);
    if (ent_type (e2) != MATRIX_DENSE_REAL)
      rerror ("cpf.weibull: argument 'a' has to be MATRIX-DENSE-REAL");
    x2 = ent_data (e2);
    nr2 = x2->nrow;
    nc2 = x2->ncol;
    if (nr2 * nc2 == 0)
      rerror ("cpf.weibull: argument 'a' has to be MATRIX-DENSE-REAL");

    // update the size of output matrix
    nr = nr > nr2 ? nr : nr2;
    nc = nc > nc2 ? nc : nc2;

    //
    // get b
    //
    e3 = bltin_get_ent (args[2]);
    if (ent_type (e3) != MATRIX_DENSE_REAL)
      rerror ("cpf.weibull: argument 'b' has to be MATRIX-DENSE-REAL");
    x3 = ent_data (e3);
    nr3 = x3->nrow;
    nc3 = x3->ncol;
    if (nr3 * nc3 == 0)
      rerror ("cpf.weibull: argument 'b' has to be MATRIX-DENSE-REAL");

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
        Mdr0 (w, i, j) = gsl_cdf_weibull_P (x, a, b);
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
ent_gsl_cpf_pareto (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *e3=0;
  int nr, nc, nr1, nc1, nr2, nc2, nr3, nc3, i, j, id, jd;
  double x, a, b;
  MDR *x1=0, *x2=0, *x3=0;
  MDR *w=0;

  if (nargs!=1 && nargs!=2 && nargs!=3)
  {
    printf ("cpf.pareto: Probability density for a Pareto distribution\n");
    printf ("cpf.pareto: with parameters 'a' and 'b'.\n");
    printf ("cpf.pareto: Format:\n");
    printf ("cpf.pareto:   cpf.pareto(x,a,b)\n");
    printf ("cpf.pareto:   cpf.pareto(x, [a,b] ), or\n");
    printf ("cpf.pareto:   cpf.pareto( [x,a,b] ) \n");
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
      rerror ("cpf.pareto: argument 'x' has to be MATRIX-DENSE-REAL");
    x1 = ent_data (e1);
    nr1 = nr = x1->nrow;
    nc1 = nc = x1->ncol;
    if (nr * nc == 0)
      rerror ("cpf.pareto: argument 'x' has to be MATRIX-DENSE-REAL");
    if (nc!=3)
      rerror ("cpf.pareto: Single argument implies [x,a,b] but it is not provided");

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
      MdrV0 (w, i)  = gsl_cdf_pareto_P (x, a, b);
    }
  }
  else if (nargs == 2)
  {
    // get x
    e1 = bltin_get_ent (args[0]);
    if (ent_type (e1) != MATRIX_DENSE_REAL)
      rerror ("cpf.pareto: argument 'x' has to be MATRIX-DENSE-REAL");
    x1 = ent_data (e1);
    nr1 = nr = x1->nrow;
    nc1 = nc = x1->ncol;
    if (nr * nc == 0)
      rerror ("cpf.pareto: argument 'x' has to be MATRIX-DENSE-REAL");

    // get [param_1, param_2]
    e2 = bltin_get_ent (args[1]);
    if (ent_type (e2) != MATRIX_DENSE_REAL)
      rerror ("cpf.pareto: argument '[a,b]' has to be MATRIX-DENSE-REAL");
    x2 = ent_data (e2);
    nr2 = x2->nrow;
    nc2 = x2->ncol;
    if (nc2 != 2)
      rerror ("cpf.pareto: second argument has to be [a,b]");

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
        Mdr0 (w, i, j) = gsl_cdf_pareto_P (x, a, b);
      }
    }
  }
  else if (nargs == 3)
  {
    // get x
    e1 = bltin_get_ent (args[0]);
    if (ent_type (e1) != MATRIX_DENSE_REAL)
      rerror ("cpf.pareto: argument 'x' has to be MATRIX-DENSE-REAL");
    x1 = ent_data (e1);
    nr1 = nr = x1->nrow;
    nc1 = nc = x1->ncol;

    //
    // get a
    //
    e2 = bltin_get_ent (args[1]);
    if (ent_type (e2) != MATRIX_DENSE_REAL)
      rerror ("cpf.pareto: argument 'a' has to be MATRIX-DENSE-REAL");
    x2 = ent_data (e2);
    nr2 = x2->nrow;
    nc2 = x2->ncol;
    if (nr2 * nc2 == 0)
      rerror ("cpf.pareto: argument 'a' has to be MATRIX-DENSE-REAL");

    // update the size of output matrix
    nr = nr > nr2 ? nr : nr2;
    nc = nc > nc2 ? nc : nc2;

    //
    // get b
    //
    e3 = bltin_get_ent (args[2]);
    if (ent_type (e3) != MATRIX_DENSE_REAL)
      rerror ("cpf.pareto: argument 'b' has to be MATRIX-DENSE-REAL");
    x3 = ent_data (e3);
    nr3 = x3->nrow;
    nc3 = x3->ncol;
    if (nr3 * nc3 == 0)
      rerror ("cpf.pareto: argument 'b' has to be MATRIX-DENSE-REAL");

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
        Mdr0 (w, i, j) = gsl_cdf_pareto_P (x, a, b);
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
ent_gsl_cpf_beta (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *e3=0;
  int nr, nc, nr1, nc1, nr2, nc2, nr3, nc3, i, j, id, jd;
  double x, a, b;
  MDR *x1=0, *x2=0, *x3=0;
  MDR *w=0;

  if (nargs!=1 && nargs!=2 && nargs!=3)
  {
    printf ("cpf.beta: Probability density for a Beta distribution\n");
    printf ("cpf.beta: with parameters 'a' and 'b'.\n");
    printf ("cpf.beta: Format:\n");
    printf ("cpf.beta:   cpf.beta(x,a,b)\n");
    printf ("cpf.beta:   cpf.beta(x, [a,b] ), or\n");
    printf ("cpf.beta:   cpf.beta( [x,a,b] ) \n");
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
      rerror ("cpf.beta: argument 'x' has to be MATRIX-DENSE-REAL");
    x1 = ent_data (e1);
    nr1 = nr = x1->nrow;
    nc1 = nc = x1->ncol;
    if (nr * nc == 0)
      rerror ("cpf.beta: argument 'x' has to be MATRIX-DENSE-REAL");
    if (nc!=3)
      rerror ("cpf.beta: Single argument implies [x,a,b] but it is not provided");

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
      MdrV0 (w, i)  = gsl_cdf_beta_P (x, a, b);
    }
  }
  else if (nargs == 2)
  {
    // get x
    e1 = bltin_get_ent (args[0]);
    if (ent_type (e1) != MATRIX_DENSE_REAL)
      rerror ("cpf.beta: argument 'x' has to be MATRIX-DENSE-REAL");
    x1 = ent_data (e1);
    nr1 = nr = x1->nrow;
    nc1 = nc = x1->ncol;
    if (nr * nc == 0)
      rerror ("cpf.beta: argument 'x' has to be MATRIX-DENSE-REAL");

    // get [param_1, param_2]
    e2 = bltin_get_ent (args[1]);
    if (ent_type (e2) != MATRIX_DENSE_REAL)
      rerror ("cpf.beta: argument '[a,b]' has to be MATRIX-DENSE-REAL");
    x2 = ent_data (e2);
    nr2 = x2->nrow;
    nc2 = x2->ncol;
    if (nc2 != 2)
      rerror ("cpf.beta: second argument has to be [a,b]");

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
        Mdr0 (w, i, j) = gsl_cdf_beta_P (x, a, b);
      }
    }
  }
  else if (nargs == 3)
  {
    // get x
    e1 = bltin_get_ent (args[0]);
    if (ent_type (e1) != MATRIX_DENSE_REAL)
      rerror ("cpf.beta: argument 'x' has to be MATRIX-DENSE-REAL");
    x1 = ent_data (e1);
    nr1 = nr = x1->nrow;
    nc1 = nc = x1->ncol;

    //
    // get a
    //
    e2 = bltin_get_ent (args[1]);
    if (ent_type (e2) != MATRIX_DENSE_REAL)
      rerror ("cpf.beta: argument 'a' has to be MATRIX-DENSE-REAL");
    x2 = ent_data (e2);
    nr2 = x2->nrow;
    nc2 = x2->ncol;
    if (nr2 * nc2 == 0)
      rerror ("cpf.beta: argument 'a' has to be MATRIX-DENSE-REAL");

    // update the size of output matrix
    nr = nr > nr2 ? nr : nr2;
    nc = nc > nc2 ? nc : nc2;

    //
    // get b
    //
    e3 = bltin_get_ent (args[2]);
    if (ent_type (e3) != MATRIX_DENSE_REAL)
      rerror ("cpf.beta: argument 'b' has to be MATRIX-DENSE-REAL");
    x3 = ent_data (e3);
    nr3 = x3->nrow;
    nc3 = x3->ncol;
    if (nr3 * nc3 == 0)
      rerror ("cpf.beta: argument 'b' has to be MATRIX-DENSE-REAL");

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
        Mdr0 (w, i, j) = gsl_cdf_beta_P (x, a, b);
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
ent_gsl_cpf_F (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *e3=0;
  int nr, nc, nr1, nc1, nr2, nc2, nr3, nc3, i, j, id, jd;
  double x, a, b;
  MDR *x1=0, *x2=0, *x3=0;
  MDR *w=0;

  if (nargs!=1 && nargs!=2 && nargs!=3)
  {
    printf ("cpf.F: Probability density for a F-distribution\n");
    printf ("cpf.F: with 'nu1' and 'nu2' degrees of freedom.\n");
    printf ("cpf.F: Format:\n");
    printf ("cpf.F:   cpf.F(x,nu1,nu2)\n");
    printf ("cpf.F:   cpf.F(x, [nu1,nu2] )\n");
    printf ("cpf.F:   cpf.F( [x,nu1,nu2] )\n");
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
      rerror ("cpf.F: argument 'x' has to be MATRIX-DENSE-REAL");
    x1 = ent_data (e1);
    nr1 = nr = x1->nrow;
    nc1 = nc = x1->ncol;
    if (nr * nc == 0)
      rerror ("cpf.F: argument 'x' has to be MATRIX-DENSE-REAL");
    if (nc!=3)
      rerror ("cpf.F: Single argument implies [x,nu1,nu2] but it is not provided");

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
      MdrV0 (w, i)  = gsl_cdf_fdist_P (x, a, b);
    }
  }
  else if (nargs == 2)
  {
    // get x
    e1 = bltin_get_ent (args[0]);
    if (ent_type (e1) != MATRIX_DENSE_REAL)
      rerror ("cpf.F: argument 'x' has to be MATRIX-DENSE-REAL");
    x1 = ent_data (e1);
    nr1 = nr = x1->nrow;
    nc1 = nc = x1->ncol;
    if (nr * nc == 0)
      rerror ("cpf.F: argument 'x' has to be MATRIX-DENSE-REAL");

    // get [param_1, param_2]
    e2 = bltin_get_ent (args[1]);
    if (ent_type (e2) != MATRIX_DENSE_REAL)
      rerror ("cpf.F: argument '[nu1,nu2]' has to be MATRIX-DENSE-REAL");
    x2 = ent_data (e2);
    nr2 = x2->nrow;
    nc2 = x2->ncol;
    if (nc2 != 2)
      rerror ("cpf.F: second argument has to be [nu1,nu2]");

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
        Mdr0 (w, i, j) = gsl_cdf_fdist_P (x, a, b);
      }
    }
  }
  else if (nargs == 3)
  {
    // get x
    e1 = bltin_get_ent (args[0]);
    if (ent_type (e1) != MATRIX_DENSE_REAL)
      rerror ("cpf.F: argument 'x' has to be MATRIX-DENSE-REAL");
    x1 = ent_data (e1);
    nr1 = nr = x1->nrow;
    nc1 = nc = x1->ncol;

    //
    // get a
    //
    e2 = bltin_get_ent (args[1]);
    if (ent_type (e2) != MATRIX_DENSE_REAL)
      rerror ("cpf.F: argument 'nu1' has to be MATRIX-DENSE-REAL");
    x2 = ent_data (e2);
    nr2 = x2->nrow;
    nc2 = x2->ncol;
    if (nr2 * nc2 == 0)
      rerror ("cpf.F: argument 'nu1' has to be MATRIX-DENSE-REAL");

    // update the size of output matrix
    nr = nr > nr2 ? nr : nr2;
    nc = nc > nc2 ? nc : nc2;

    //
    // get b
    //
    e3 = bltin_get_ent (args[2]);
    if (ent_type (e3) != MATRIX_DENSE_REAL)
      rerror ("cpf.F: argument 'nu2' has to be MATRIX-DENSE-REAL");
    x3 = ent_data (e3);
    nr3 = x3->nrow;
    nc3 = x3->ncol;
    if (nr3 * nc3 == 0)
      rerror ("cpf.F: argument 'nu2' has to be MATRIX-DENSE-REAL");

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
        Mdr0 (w, i, j) = gsl_cdf_fdist_P (x, a, b);
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
ent_gsl_cpf_lognormal (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *e3=0;
  int nr, nc, nr1, nc1, nr2, nc2, nr3, nc3, i, j, id, jd;
  double x, a, b;
  MDR *x1=0, *x2=0, *x3=0;
  MDR *w=0;

  if (nargs!=1 && nargs!=2 && nargs!=3)
  {
    printf
        ("cpf.lognormal: Probability density for a lognormal distribution\n");
    printf ("cpf.lognormal: with parameters 'zeta' and 'sigma'.\n");
    printf ("cpf.lognormal: Format:\n");
    printf ("cpf.lognormal:   cpf.lognormal(x,zeta,sigma)\n");
    printf ("cpf.lognormal:   cpf.lognormal(x, [zeta,sigma] )\n");
    printf ("cpf.lognormal:   cpf.lognormal( [x,zeta,sigma] )\n");
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
      rerror ("cpf.lognormal: argument 'x' has to be MATRIX-DENSE-REAL");
    x1 = ent_data (e1);
    nr1 = nr = x1->nrow;
    nc1 = nc = x1->ncol;
    if (nr * nc == 0)
      rerror ("cpf.lognormal: argument 'x' has to be MATRIX-DENSE-REAL");
    if (nc!=3)
      rerror ("cpf.lognormal: Single argument implies [x,zeta,sigma] but it is not provided");

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
      MdrV0 (w, i)  = gsl_cdf_lognormal_P (x, a, b);
    }
  }
  else if (nargs == 2)
  {
    // get x
    e1 = bltin_get_ent (args[0]);
    if (ent_type (e1) != MATRIX_DENSE_REAL)
      rerror ("cpf.lognormal: argument 'x' has to be MATRIX-DENSE-REAL");
    x1 = ent_data (e1);
    nr1 = nr = x1->nrow;
    nc1 = nc = x1->ncol;
    if (nr * nc == 0)
      rerror ("cpf.lognormal: argument 'x' has to be MATRIX-DENSE-REAL");

    // get [param_1, param_2]
    e2 = bltin_get_ent (args[1]);
    if (ent_type (e2) != MATRIX_DENSE_REAL)
      rerror ("cpf.lognormal: argument '[zeta,sigma]' has to be MATRIX-DENSE-REAL");
    x2 = ent_data (e2);
    nr2 = x2->nrow;
    nc2 = x2->ncol;
    if (nc2 != 2)
      rerror ("cpf.lognormal: second argument has to be [zeta,sigma]");

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
        Mdr0 (w, i, j) = gsl_cdf_lognormal_P (x, a, b);
      }
    }
  }
  else if (nargs == 3)
  {
    // get x
    e1 = bltin_get_ent (args[0]);
    if (ent_type (e1) != MATRIX_DENSE_REAL)
      rerror ("cpf.lognormal: argument 'x' has to be MATRIX-DENSE-REAL");
    x1 = ent_data (e1);
    nr1 = nr = x1->nrow;
    nc1 = nc = x1->ncol;

    //
    // get a
    //
    e2 = bltin_get_ent (args[1]);
    if (ent_type (e2) != MATRIX_DENSE_REAL)
      rerror ("cpf.lognormal: argument 'zeta' has to be MATRIX-DENSE-REAL");
    x2 = ent_data (e2);
    nr2 = x2->nrow;
    nc2 = x2->ncol;
    if (nr2 * nc2 == 0)
      rerror ("cpf.lognormal: argument 'zeta' has to be MATRIX-DENSE-REAL");

    // update the size of output matrix
    nr = nr > nr2 ? nr : nr2;
    nc = nc > nc2 ? nc : nc2;

    //
    // get b
    //
    e3 = bltin_get_ent (args[2]);
    if (ent_type (e3) != MATRIX_DENSE_REAL)
      rerror ("cpf.lognormal: argument 'sigma' has to be MATRIX-DENSE-REAL");
    x3 = ent_data (e3);
    nr3 = x3->nrow;
    nc3 = x3->ncol;
    if (nr3 * nc3 == 0)
      rerror ("cpf.lognormal: argument 'sigma' has to be MATRIX-DENSE-REAL");

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
        Mdr0 (w, i, j) = gsl_cdf_lognormal_P (x, a, b);
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
ent_gsl_cpf_gamma (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *e3=0;
  int nr, nc, nr1, nc1, nr2, nc2, nr3, nc3, i, j, id, jd;
  double x, a, b;
  MDR *x1=0, *x2=0, *x3=0;
  MDR *w=0;

  if (nargs!=1 && nargs!=2 && nargs!=3)
  {
    printf ("cpf.gamma: Probability density for a Beta distribution\n");
    printf ("cpf.gamma: with parameters 'a' and 'b'.\n");
    printf ("cpf.gamma: Format:\n");
    printf ("cpf.gamma:   cpf.gamma(x,a,b)\n");
    printf ("cpf.gamma:   cpf.gamma(x, [a,b] ), or\n");
    printf ("cpf.gamma:   cpf.gamma( [x,a,b] ) \n");
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
      rerror ("cpf.gamma: argument 'x' has to be MATRIX-DENSE-REAL");
    x1 = ent_data (e1);
    nr1 = nr = x1->nrow;
    nc1 = nc = x1->ncol;
    if (nr * nc == 0)
      rerror ("cpf.gamma: argument 'x' has to be MATRIX-DENSE-REAL");
    if (nc!=3)
      rerror ("cpf.gamma: Single argument implies [x,a,b] but it is not provided");

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
      MdrV0 (w, i)  = gsl_cdf_gamma_P (x, a, b);
    }
  }
  else if (nargs == 2)
  {
    // get x
    e1 = bltin_get_ent (args[0]);
    if (ent_type (e1) != MATRIX_DENSE_REAL)
      rerror ("cpf.gamma: argument 'x' has to be MATRIX-DENSE-REAL");
    x1 = ent_data (e1);
    nr1 = nr = x1->nrow;
    nc1 = nc = x1->ncol;
    if (nr * nc == 0)
      rerror ("cpf.gamma: argument 'x' has to be MATRIX-DENSE-REAL");

    // get [param_1, param_2]
    e2 = bltin_get_ent (args[1]);
    if (ent_type (e2) != MATRIX_DENSE_REAL)
      rerror ("cpf.gamma: argument '[a,b]' has to be MATRIX-DENSE-REAL");
    x2 = ent_data (e2);
    nr2 = x2->nrow;
    nc2 = x2->ncol;
    if (nc2 != 2)
      rerror ("cpf.gamma: second argument has to be [a,b]");

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
        Mdr0 (w, i, j) = gsl_cdf_gamma_P (x, a, b);
      }
    }
  }
  else if (nargs == 3)
  {
    // get x
    e1 = bltin_get_ent (args[0]);
    if (ent_type (e1) != MATRIX_DENSE_REAL)
      rerror ("cpf.gamma: argument 'x' has to be MATRIX-DENSE-REAL");
    x1 = ent_data (e1);
    nr1 = nr = x1->nrow;
    nc1 = nc = x1->ncol;

    //
    // get a
    //
    e2 = bltin_get_ent (args[1]);
    if (ent_type (e2) != MATRIX_DENSE_REAL)
      rerror ("cpf.gamma: argument 'a' has to be MATRIX-DENSE-REAL");
    x2 = ent_data (e2);
    nr2 = x2->nrow;
    nc2 = x2->ncol;
    if (nr2 * nc2 == 0)
      rerror ("cpf.gamma: argument 'a' has to be MATRIX-DENSE-REAL");

    // update the size of output matrix
    nr = nr > nr2 ? nr : nr2;
    nc = nc > nc2 ? nc : nc2;

    //
    // get b
    //
    e3 = bltin_get_ent (args[2]);
    if (ent_type (e3) != MATRIX_DENSE_REAL)
      rerror ("cpf.gamma: argument 'b' has to be MATRIX-DENSE-REAL");
    x3 = ent_data (e3);
    nr3 = x3->nrow;
    nc3 = x3->ncol;
    if (nr3 * nc3 == 0)
      rerror ("cpf.gamma: argument 'b' has to be MATRIX-DENSE-REAL");

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
        Mdr0 (w, i, j) = gsl_cdf_gamma_P (x, a, b);
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
ent_gsl_cpf_exppow (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *e3=0;
  int nr, nc, nr1, nc1, nr2, nc2, nr3, nc3, i, j, id, jd;
  double x, a, b;
  MDR *x1=0, *x2=0, *x3=0;
  MDR *w=0;

  if (nargs!=1 && nargs!=2 && nargs!=3)
  {
    printf ("cpf.exppow: Probability density for a Beta distribution\n");
    printf ("cpf.exppow: with parameters 'a' and 'b'.\n");
    printf ("cpf.exppow: Format:\n");
    printf ("cpf.exppow:   cpf.exppow(x,a,b)\n");
    printf ("cpf.exppow:   cpf.exppow(x, [a,b] ), or\n");
    printf ("cpf.exppow:   cpf.exppow( [x,a,b] ) \n");
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
      rerror ("cpf.exppow: argument 'x' has to be MATRIX-DENSE-REAL");
    x1 = ent_data (e1);
    nr1 = nr = x1->nrow;
    nc1 = nc = x1->ncol;
    if (nr * nc == 0)
      rerror ("cpf.exppow: argument 'x' has to be MATRIX-DENSE-REAL");
    if (nc!=3)
      rerror ("cpf.exppow: Single argument implies [x,a,b] but it is not provided");

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
      MdrV0 (w, i)  = gsl_cdf_exppow_P (x, a, b);
    }
  }
  else if (nargs == 2)
  {
    // get x
    e1 = bltin_get_ent (args[0]);
    if (ent_type (e1) != MATRIX_DENSE_REAL)
      rerror ("cpf.exppow: argument 'x' has to be MATRIX-DENSE-REAL");
    x1 = ent_data (e1);
    nr1 = nr = x1->nrow;
    nc1 = nc = x1->ncol;
    if (nr * nc == 0)
      rerror ("cpf.exppow: argument 'x' has to be MATRIX-DENSE-REAL");

    // get [param_1, param_2]
    e2 = bltin_get_ent (args[1]);
    if (ent_type (e2) != MATRIX_DENSE_REAL)
      rerror ("cpf.exppow: argument '[a,b]' has to be MATRIX-DENSE-REAL");
    x2 = ent_data (e2);
    nr2 = x2->nrow;
    nc2 = x2->ncol;
    if (nc2 != 2)
      rerror ("cpf.exppow: second argument has to be [a,b]");

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
        Mdr0 (w, i, j) = gsl_cdf_exppow_P (x, a, b);
      }
    }
  }
  else if (nargs == 3)
  {
    // get x
    e1 = bltin_get_ent (args[0]);
    if (ent_type (e1) != MATRIX_DENSE_REAL)
      rerror ("cpf.exppow: argument 'x' has to be MATRIX-DENSE-REAL");
    x1 = ent_data (e1);
    nr1 = nr = x1->nrow;
    nc1 = nc = x1->ncol;

    //
    // get a
    //
    e2 = bltin_get_ent (args[1]);
    if (ent_type (e2) != MATRIX_DENSE_REAL)
      rerror ("cpf.exppow: argument 'a' has to be MATRIX-DENSE-REAL");
    x2 = ent_data (e2);
    nr2 = x2->nrow;
    nc2 = x2->ncol;
    if (nr2 * nc2 == 0)
      rerror ("cpf.exppow: argument 'a' has to be MATRIX-DENSE-REAL");

    // update the size of output matrix
    nr = nr > nr2 ? nr : nr2;
    nc = nc > nc2 ? nc : nc2;

    //
    // get b
    //
    e3 = bltin_get_ent (args[2]);
    if (ent_type (e3) != MATRIX_DENSE_REAL)
      rerror ("cpf.exppow: argument 'b' has to be MATRIX-DENSE-REAL");
    x3 = ent_data (e3);
    nr3 = x3->nrow;
    nc3 = x3->ncol;
    if (nr3 * nc3 == 0)
      rerror ("cpf.exppow: argument 'b' has to be MATRIX-DENSE-REAL");

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
        Mdr0 (w, i, j) = gsl_cdf_exppow_P (x, a, b);
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
ent_gsl_cpf_binomial (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *e3=0;
  int nr, nc, nr1, nc1, nr2, nc2, nr3, nc3, i, j, id, jd;
  double x, a, b;
  MDR *x1=0, *x2=0, *x3=0;
  MDR *w=0;

  if (nargs!=1 && nargs!=3)
  {
    printf ("cpf.binomial: Probability density for a binomial distribution\n");
    printf
        ("cpf.binomial: with parameters: 'p', probability of a single trial,\n");
    printf ("cpf.binomial: and 'n', number of independent trials.\n");
    printf ("cpf.binomial: Format:\n");
    printf ("cpf.binomial:   cpf.binomial(k,p,n)\n");
    printf ("cpf.binomial:   cpf.binomial( [k,p,n] )\n");
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
      rerror ("cpf.binomial: first argument has to be MATRIX-DENSE-REAL");
    x1 = ent_data (e1);
    nr1 = nr = x1->nrow;
    nc1 = nc = x1->ncol;
    if (nr * nc == 0)
      rerror ("cpf.binomial: first argument has to be MATRIX-DENSE-REAL");
    if (nc!=3)
      rerror ("cpf.binomial: Single argument implies [k,p,n] but it is not provided");

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
      MdrV0 (w, i) = gsl_cdf_binomial_P ( (int) x, a, (int) b);
    }
  }
  else if (nargs == 3)
  {
    // get x
    e1 = bltin_get_ent (args[0]);
    if (ent_type (e1) != MATRIX_DENSE_REAL)
      rerror ("cpf.binomial: argument 'k' has to be MATRIX-DENSE-REAL");
    x1 = ent_data (e1);
    nr1 = nr = x1->nrow;
    nc1 = nc = x1->ncol;

    //
    // get a
    //
    e2 = bltin_get_ent (args[1]);
    if (ent_type (e2) != MATRIX_DENSE_REAL)
      rerror ("cpf.binomial: argument 'p' has to be MATRIX-DENSE-REAL");
    x2 = ent_data (e2);
    nr2 = x2->nrow;
    nc2 = x2->ncol;
    if (nr2 * nc2 == 0)
      rerror ("cpf.binomial: argument 'p' has to be MATRIX-DENSE-REAL");

    // update the size of output matrix
    nr = nr > nr2 ? nr : nr2;
    nc = nc > nc2 ? nc : nc2;

    //
    // get b
    //
    e3 = bltin_get_ent (args[2]);
    if (ent_type (e3) != MATRIX_DENSE_REAL)
      rerror ("cpf.binomial: argument 'n' has to be MATRIX-DENSE-REAL");
    x3 = ent_data (e3);
    nr3 = x3->nrow;
    nc3 = x3->ncol;
    if (nr3 * nc3 == 0)
      rerror ("cpf.binomial: argument 'n' has to be MATRIX-DENSE-REAL");

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
        Mdr0 (w, i, j) = gsl_cdf_binomial_P ( (int) x, a, (int) b);
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
ent_gsl_cpf_negbinomial (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *e3=0;
  int nr, nc, nr1, nc1, nr2, nc2, nr3, nc3, i, j, id, jd;
  double x, a, b;
  MDR *x1=0, *x2=0, *x3=0;
  MDR *w=0;

  if (nargs!=1 && nargs!=3)
  {
    printf ("cpf.negbinomial: Probability density for a binomial distribution\n");
    printf
        ("cpf.negbinomial: with parameters: 'p', probability of a single trial,\n");
    printf ("cpf.negbinomial: and 'n', number of independent trials.\n");
    printf ("cpf.negbinomial: Format:\n");
    printf ("cpf.negbinomial:   cpf.negbinomial(k,p,n)\n");
    printf ("cpf.negbinomial:   cpf.negbinomial( [k,p,n] )\n");
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
      rerror ("cpf.negbinomial: first argument has to be MATRIX-DENSE-REAL");
    x1 = ent_data (e1);
    nr1 = nr = x1->nrow;
    nc1 = nc = x1->ncol;
    if (nr * nc == 0)
      rerror ("cpf.negbinomial: first argument has to be MATRIX-DENSE-REAL");
    if (nc!=3)
      rerror ("cpf.negbinomial: Single argument implies [k,p,n] but it is not provided");

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
      MdrV0 (w, i) = gsl_cdf_negative_binomial_P ( (int) x, a, (int) b);
    }
  }
  else if (nargs == 3)
  {
    // get x
    e1 = bltin_get_ent (args[0]);
    if (ent_type (e1) != MATRIX_DENSE_REAL)
      rerror ("cpf.negbinomial: argument 'k' has to be MATRIX-DENSE-REAL");
    x1 = ent_data (e1);
    nr1 = nr = x1->nrow;
    nc1 = nc = x1->ncol;

    //
    // get a
    //
    e2 = bltin_get_ent (args[1]);
    if (ent_type (e2) != MATRIX_DENSE_REAL)
      rerror ("cpf.negbinomial: argument 'p' has to be MATRIX-DENSE-REAL");
    x2 = ent_data (e2);
    nr2 = x2->nrow;
    nc2 = x2->ncol;
    if (nr2 * nc2 == 0)
      rerror ("cpf.negbinomial: argument 'p' has to be MATRIX-DENSE-REAL");

    // update the size of output matrix
    nr = nr > nr2 ? nr : nr2;
    nc = nc > nc2 ? nc : nc2;

    //
    // get b
    //
    e3 = bltin_get_ent (args[2]);
    if (ent_type (e3) != MATRIX_DENSE_REAL)
      rerror ("cpf.negbinomial: argument 'n' has to be MATRIX-DENSE-REAL");
    x3 = ent_data (e3);
    nr3 = x3->nrow;
    nc3 = x3->ncol;
    if (nr3 * nc3 == 0)
      rerror ("cpf.negbinomial: argument 'n' has to be MATRIX-DENSE-REAL");

    // update the size of output matrix to max of all three sizes
    nr = MAX(nr,nr3);
    nc = MAX(nc,nc3);
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
        Mdr0 (w, i, j) = gsl_cdf_negative_binomial_P ( (int) x, a, (int) b);
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
ent_gsl_cpf_hypergeom (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *e3=0, *e4=0;
  int nr=0, nc=0, nr1=0, nc1=0, nr2=0, nc2=0, nr3=0, nc3=0, nr4=0, nc4=0, i=0, j=0, id=0, jd=0;
  unsigned int k, n1, n2, t;
  MDR *x1=0, *x2=0, *x3=0, *x4=0;
  MDR *w=0;
  if (nargs!=1 && nargs!=4)
  {
    printf
      ("cpf.hypergeom: Probability density for a hypergeometric distribution\n");
    printf ("cpf.hypergeom: with parameters 'n1', 'n2' and 't'.\n");
    printf ("cpf.hypergeom: Format:\n");
    printf ("cpf.hypergeom:   cpf.hypergeom(k,n1,n2,t)\n");
    printf ("cpf.hypergeom:   cpf.hypergeom( [k,n1,n2,t] )\n");
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
      rerror ("cpf.hypergeom: first argument has to be MATRIX-DENSE-REAL");
    x1 = ent_data (e1);
    nr1 = nr = x1->nrow;
    nc1 = nc = x1->ncol;
    if (nr * nc == 0)
      rerror ("cpf.hypergeom: first argument has to be MATRIX-DENSE-REAL");
    if (nc!=4)
      rerror ("cpf.hypergeom: Single argument implies [k,n1,n2,t] but it is not provided");

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
      MdrV0 (w, i) = gsl_cdf_hypergeometric_P (k, n1, n2, t);
    }
  }
  else if (nargs == 4)
  {
    //
    // get n1
    //
    e2 = bltin_get_ent (args[1]);
    if (ent_type (e2) != MATRIX_DENSE_REAL)
      rerror ("cpf.hypergeom: argument 'n1' has to be MATRIX-DENSE-REAL");
    x2 = ent_data (e2);
    nr2 = x2->nrow;
    nc2 = x2->ncol;
    if (nr2 * nc2 == 0)
      rerror ("cpf.hypergeom: argument 'n1' has to be MATRIX-DENSE-REAL");
    // update the size of output matrix
    nr = nr > nr2 ? nr : nr2;
    nc = nc > nc2 ? nc : nc2;

    //
    // get n2
    //
    e3 = bltin_get_ent (args[2]);
    if (ent_type (e3) != MATRIX_DENSE_REAL)
      rerror ("cpf.hypergeom: argument 'n2' has to be MATRIX-DENSE-REAL");
    x3 = ent_data (e3);
    nr3 = x3->nrow;
    nc3 = x3->ncol;
    if (nr3 * nc3 == 0)
      rerror ("cpf.hypergeom: argument 'n' has to be MATRIX-DENSE-REAL");
    // update the size of output matrix
    nr = nr > nr3 ? nr : nr3;
    nc = nc > nc3 ? nc : nc3;

    //
    // get t
    //
    e4 = bltin_get_ent (args[3]);
    if (ent_type (e4) != MATRIX_DENSE_REAL)
      rerror ("cpf.hypergeom: argument 't' has to be MATRIX-DENSE-REAL");
    x4 = ent_data (e4);
    nr4 = x4->nrow;
    nc4 = x4->ncol;
    if (nr4 * nc4 == 0)
      rerror ("cpf.hypergeom: argument 't' has to be MATRIX-DENSE-REAL");
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
        Mdr0 (w, i, j) = gsl_cdf_hypergeometric_P (k, n1, n2, t);
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

