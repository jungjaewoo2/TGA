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
#include "rlab_macros.h"

// gsl headers
// shared object
#include <gsl/gsl_mode.h>
#include <gsl/gsl_precision.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_machine.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>

// include: standard C headers
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <math.h>

#define PI 3.1415926
#define GSL_WARNINGS_OFF



// ----------------------------------------------------------
//
// GSL distributions
//
// ----------------------------------------------------------

// ********************************************************************************
// ********************************************************************************
//
// 2-parameter distributions
//
//
// ********************************************************************************
// ********************************************************************************
RLAB_MDR_CALL_2ARGS_FUNC("pdf.logistic",gsl_ran_logistic_pdf,
                         "Probability density function for logistic distribution\n");
RLAB_MDR_CALL_2ARGS_FUNC("cdf.logistic",gsl_cdf_logistic_P,
                         "Cummulative density function for logistic distribution\n");
RLAB_MDR_CALL_2ARGS_FUNC("cdfinv.logistic",gsl_cdf_logistic_Pinv,
                         "Inverse cummulative density function for logistic distribution\n");

RLAB_MDR_CALL_2ARGS_FUNC("pdf.t",gsl_ran_tdist_pdf,
                         "Probability density function for t-distribution\n");
RLAB_MDR_CALL_2ARGS_FUNC("cdf.t",gsl_cdf_tdist_P,
                         "Cummulative density function for t-distribution\n");
RLAB_MDR_CALL_2ARGS_FUNC("cdfinv.t",gsl_cdf_tdist_Pinv,
                         "Inverse cummulative density function for t-distribution\n");

RLAB_MDR_CALL_2ARGS_FUNC("pdf.chisq",gsl_ran_chisq_pdf,
                         "Probability density function for chi-square distribution\n");
RLAB_MDR_CALL_2ARGS_FUNC("cdf.chisq",gsl_cdf_chisq_P,
                         "Cummulative density function for chi-square distribution\n");
RLAB_MDR_CALL_2ARGS_FUNC("cdfinv.chisq",gsl_cdf_chisq_Pinv,
                         "Inverse cummulative density function for chi-square distribution\n");

RLAB_MDR_CALL_2ARGS_FUNC("pdf.rayleigh",gsl_ran_rayleigh_pdf,
                         "Probability density function for rayleigh distribution\n");
RLAB_MDR_CALL_2ARGS_FUNC("cdf.rayleigh",gsl_cdf_rayleigh_P,
                         "Cummulative density function for rayleigh distribution\n");
RLAB_MDR_CALL_2ARGS_FUNC("cdfinv.rayleigh",gsl_cdf_rayleigh_Pinv,
                         "Inverse cummulative density function for rayleigh distribution\n");

RLAB_MDR_CALL_2ARGS_FUNC("pdf.cauchy",gsl_ran_cauchy_pdf,
                         "Probability density function for cauchy distribution\n");
RLAB_MDR_CALL_2ARGS_FUNC("cdf.cauchy",gsl_cdf_cauchy_P,
                         "Cummulative density function for cauchy distribution\n");
RLAB_MDR_CALL_2ARGS_FUNC("cdfinv.cauchy",gsl_cdf_cauchy_Pinv,
                         "Inverse cummulative density function for cauchy distribution\n");

RLAB_MDR_CALL_2ARGS_FUNC("pdf.laplace",gsl_ran_laplace_pdf,
                         "Probability density function for laplace distribution\n");
RLAB_MDR_CALL_2ARGS_FUNC("cdf.laplace",gsl_cdf_laplace_P,
                         "Cummulative density function for laplace distribution\n");
RLAB_MDR_CALL_2ARGS_FUNC("cdfinv.laplace",gsl_cdf_laplace_Pinv,
                         "Inverse cummulative density function for laplace distribution\n");

RLAB_MDR_CALL_2ARGS_FUNC("pdf.exp",gsl_ran_exponential_pdf,
                         "Probability density function for exponential distribution\n");
RLAB_MDR_CALL_2ARGS_FUNC("cdf.exp",gsl_cdf_exponential_P,
                         "Cummulative density function for exponential distribution\n");
RLAB_MDR_CALL_2ARGS_FUNC("cdfinv.exp",gsl_cdf_exponential_Pinv,
                         "Inverse cummulative density function for exponential distribution\n");

RLAB_MDR_CALL_2ARGS_FUNC("pdf.poisson",gsl_ran_poisson_pdf,
                         "Probability density function for poisson distribution\n");
RLAB_MDR_CALL_2ARGS_FUNC("cdf.poisson",gsl_cdf_poisson_P,
                         "Cummulative density function for poisson distribution\n");

RLAB_MDR_CALL_2ARGS_FUNC("pdf.geometric",gsl_ran_geometric_pdf,
                         "Probability density function for geometric distribution\n");
RLAB_MDR_CALL_2ARGS_FUNC("cdf.geometric",gsl_cdf_geometric_P,
                         "Cummulative density function for geometric distribution\n");


RLAB_MDR_CALL_2ARGS_FUNC("pdf.log",gsl_ran_logarithmic_pdf,
                         "Probability density function for logarithmic distribution\n");

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

double rlab2gsl_pdf_normal(double x, double m, double s)
{
  return gsl_ran_gaussian_pdf (x-m, s);
}
RLAB_MDR_CALL_3ARGS_FUNC("pdf.normal", rlab2gsl_pdf_normal,
                         "Probability density function for gaussian/normal distribution\n");

double rlab2gsl_cdf_normal_P(double x, double m, double s)
{
  return gsl_cdf_gaussian_P(x-m, s);
}
RLAB_MDR_CALL_3ARGS_FUNC("cdf.normal", rlab2gsl_cdf_normal_P,
                         "Cummulative density function for gaussian/normal distribution\n");

double rlab2gsl_cdf_normal_invP(double x, double m, double s)
{
  return (m+gsl_cdf_gaussian_Pinv(x, s));
}
RLAB_MDR_CALL_3ARGS_FUNC("cdfinv.normal",rlab2gsl_cdf_normal_invP,
                         "Inverse cummulative density function for gaussian/normal distribution\n");

RLAB_MDR_CALL_3ARGS_FUNC("pdf.pascal",gsl_ran_pascal_pdf,
                         "Probability density function for pascal distribution\n");
RLAB_MDR_CALL_3ARGS_FUNC("cdf.pascal",gsl_cdf_pascal_P,
                         "Cummulative density function for pascal distribution\n");

RLAB_MDR_CALL_3ARGS_FUNC("pdf.uniform",gsl_ran_flat_pdf,
                         "Probability density function for uniform distribution\n");
RLAB_MDR_CALL_3ARGS_FUNC("cdf.uniform",gsl_cdf_flat_P,
                         "Cummulative density function for uniform distribution\n");
RLAB_MDR_CALL_3ARGS_FUNC("cdfinv.uniform",gsl_cdf_flat_Pinv,
                         "Inverse cummulative density function for uniform distribution\n");

RLAB_MDR_CALL_3ARGS_FUNC("pdf.gumbel1",gsl_ran_gumbel1_pdf,
                         "Probability density function for Gumbel type-1 distribution\n");
RLAB_MDR_CALL_3ARGS_FUNC("cdf.gumbel1",gsl_cdf_gumbel1_P,
                         "Cummulative density function for Gumbel type-1  distribution\n");
RLAB_MDR_CALL_3ARGS_FUNC("cdfinv.gumbel1",gsl_cdf_gumbel1_Pinv,
                         "Inverse cummulative density function for Gumbel type-1 distribution\n");

RLAB_MDR_CALL_3ARGS_FUNC("pdf.gumbel2",gsl_ran_gumbel2_pdf,
                         "Probability density function for Gumbel type-2 distribution\n");
RLAB_MDR_CALL_3ARGS_FUNC("cdf.gumbel2",gsl_cdf_gumbel2_P,
                         "Cummulative density function for Gumbel type-2  distribution\n");
RLAB_MDR_CALL_3ARGS_FUNC("cdfinv.gumbel2",gsl_cdf_gumbel2_Pinv,
                         "Inverse cummulative density function for Gumbel type-2 distribution\n");

RLAB_MDR_CALL_3ARGS_FUNC("pdf.weibull",gsl_ran_weibull_pdf,
                         "Probability density function for weibull distribution\n");
RLAB_MDR_CALL_3ARGS_FUNC("cdf.weibull",gsl_cdf_weibull_P,
                         "Cummulative density function for weibull distribution\n");
RLAB_MDR_CALL_3ARGS_FUNC("cdfinv.weibull",gsl_cdf_weibull_Pinv,
                         "Inverse cummulative density function for weibull distribution\n");

RLAB_MDR_CALL_3ARGS_FUNC("pdf.pareto",gsl_ran_pareto_pdf,
                         "Probability density function for pareto distribution\n");
RLAB_MDR_CALL_3ARGS_FUNC("cdf.pareto",gsl_cdf_pareto_P,
                         "Cummulative density function for pareto distribution\n");
RLAB_MDR_CALL_3ARGS_FUNC("cdfinv.pareto",gsl_cdf_pareto_Pinv,
                         "Inverse cummulative density function for pareto distribution\n");

RLAB_MDR_CALL_3ARGS_FUNC("pdf.beta",gsl_ran_beta_pdf,
                         "Probability density function for beta distribution\n");
RLAB_MDR_CALL_3ARGS_FUNC("cdf.beta",gsl_cdf_beta_P,
                         "Cummulative density function for beta distribution\n");
RLAB_MDR_CALL_3ARGS_FUNC("cdfinv.beta",gsl_cdf_beta_Pinv,
                         "Inverse cummulative density function for beta distribution\n");

RLAB_MDR_CALL_3ARGS_FUNC("pdf.F",gsl_ran_fdist_pdf,
                         "Probability density function for F-distribution\n");
RLAB_MDR_CALL_3ARGS_FUNC("cdf.F",gsl_cdf_fdist_P,
                         "Cummulative density function for F-distribution\n");
RLAB_MDR_CALL_3ARGS_FUNC("cdfinv.F",gsl_cdf_fdist_Pinv,
                         "Inverse cummulative density function for F-distribution\n");

RLAB_MDR_CALL_3ARGS_FUNC("pdf.lognormal",gsl_ran_lognormal_pdf,
                         "Probability density function for lognormal distribution\n");
RLAB_MDR_CALL_3ARGS_FUNC("cdf.lognormal",gsl_cdf_lognormal_P,
                         "Cummulative density function for lognormal distribution\n");
RLAB_MDR_CALL_3ARGS_FUNC("cdfinv.lognormal",gsl_cdf_lognormal_Pinv,
                         "Inverse cummulative density function for lognormal distribution\n");

RLAB_MDR_CALL_3ARGS_FUNC("pdf.gamma",gsl_ran_gamma_pdf,
                         "Probability density function for gamma distribution\n");
RLAB_MDR_CALL_3ARGS_FUNC("cdf.gamma",gsl_cdf_gamma_P,
                         "Cummulative density function for gamma distribution\n");
RLAB_MDR_CALL_3ARGS_FUNC("cdfinv.gamma",gsl_cdf_gamma_Pinv,
                         "Inverse cummulative density function for gamma distribution\n");

RLAB_MDR_CALL_3ARGS_FUNC("pdf.rayleightail",gsl_ran_rayleigh_tail_pdf,
                         "Probability density function for rayleigh tail distribution\n");

RLAB_MDR_CALL_3ARGS_FUNC("pdf.exppow",gsl_ran_exppow_pdf,
                         "Probability density function for exponential power distribution\n");
RLAB_MDR_CALL_3ARGS_FUNC("cdf.exppow",gsl_cdf_exppow_P,
                         "Cummulative density function for exponential power distribution\n");

RLAB_MDR_CALL_3ARGS_FUNC("pdf.normaltail",gsl_ran_gaussian_tail_pdf,
                         "Probability density function for normal tail distribution\n");

RLAB_MDR_CALL_3ARGS_FUNC("pdf.binomial",gsl_ran_binomial_pdf,
                         "Probability density function for binomial distribution\n");
RLAB_MDR_CALL_3ARGS_FUNC("cdf.binomial",gsl_cdf_binomial_P,
                         "Cummulative density function for binomial distribution\n");

RLAB_MDR_CALL_3ARGS_FUNC("pdf.negbinomial",gsl_ran_negative_binomial_pdf,
                         "Probability density function for negative binomial distribution\n");
RLAB_MDR_CALL_3ARGS_FUNC("cdf.negbinomial",gsl_cdf_negative_binomial_P,
                         "Cummulative density function for negative binomial distribution\n");

RLAB_MDR_CALL_4ARGS_FUNC("pdf.hypergeom",gsl_ran_hypergeometric_pdf,
                         "Probability density for hypergeometric distribution\n",
                         "\tpdf.hypergeom(k,n1,n2,t) or  pdf.hypergeom( [k,n1,n2,t]\n");
RLAB_MDR_CALL_4ARGS_FUNC("cdf.hypergeom",gsl_cdf_hypergeometric_P,
                         "Cummulative density for hypergeometric distribution\n",
                         "\tcdf.hypergeom(k,n1,n2,t) or  cdf.hypergeom( [k,n1,n2,t]\n");

RLABENT_CALL_D_FUNC_I_DA_UA ("pdf.multinom",gsl_ran_multinomial_pdf,
                             "Probability density for multinomial distribution\n",
                             "pdf.multinom(p,n)\n");

RLABENT_CALL_D_FUNC_I_DA_DA("pdf.dirichlet",gsl_ran_dirichlet_pdf,
                            "Probability density for Dirichlet distribution\n",
                            "pdf.dirichlet(alpha,theta)\n" );

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



// Ent * ent_gsl_cpf_normal (int nargs, Datum args[])
// {
//   Ent *e1=0, *e2=0, *e3=0;
//   int nr, nc, nr1, nc1, nr2, nc2, nr3, nc3, i, j, id, jd;
//   double x, a, b;
//   MDR *x1=0, *x2=0, *x3=0;
//   MDR *w=0;
// 
//   if (nargs!=1 && nargs!=2 && nargs!=3)
//   {
//     printf ("cpf.normal: Probability density for a Gaussian distribution with mean,\n");
//     printf ("cpf.normal: and standard deviation.\n");
//     printf ("cpf.normal: Format:\n");
//     printf ("cpf.normal:   cpf.normal(x, mean, std), or\n");
//     printf ("cpf.normal:   cpf.normal(x, [mean,std] ), or\n");
//     printf ("cpf.normal:   cpf.normal( [x,mean,std] ) \n");
//     rerror ("One, two or three arguments required");
//   }
// 
//   //
//   // get x
//   //
//   if (nargs == 1)
//   {
//     // get [x, param_1, param_2]
//     e1 = bltin_get_ent (args[0]);
//     if (ent_type (e1) != MATRIX_DENSE_REAL)
//       rerror ("cpf.normal: argument 'x' has to be MATRIX-DENSE-REAL");
//     x1 = ent_data (e1);
//     nr1 = nr = x1->nrow;
//     nc1 = nc = x1->ncol;
//     if (nr * nc == 0)
//       rerror ("cpf.normal: argument 'x' has to be MATRIX-DENSE-REAL");
//     if (nc!=3)
//       rerror ("cpf.normal: Single argument implies [x,mean,std] but it is not provided");
// 
//     w = mdr_Create (nr, 1);
//     for (i = 0; i < nr; i++)
//     {
//       if (x1->type == RLAB_TYPE_INT32)
//       {
//         x = (double) Mdi0 (x1, i, 0);
//         a = (double) Mdi0 (x1, i, 1);
//         b = (double) Mdi0 (x1, i, 2);
//       }
//       else
//       {
//         x = Mdr0 (x1, i, 0);
//         a = Mdr0 (x1, i, 1);
//         b = Mdr0 (x1, i, 2);
//       }
//       MdrV0 (w, i)  = gsl_cdf_gaussian_P (x-a, b);
//     }
//   }
//   else if (nargs == 2)
//   {
//     // get x
//     e1 = bltin_get_ent (args[0]);
//     if (ent_type (e1) != MATRIX_DENSE_REAL)
//       rerror ("cpf.normal: argument 'x' has to be MATRIX-DENSE-REAL");
//     x1 = ent_data (e1);
//     nr1 = nr = x1->nrow;
//     nc1 = nc = x1->ncol;
//     if (nr * nc == 0)
//       rerror ("cpf.normal: argument 'x' has to be MATRIX-DENSE-REAL");
// 
//     // get [param_1, param_2]
//     e2 = bltin_get_ent (args[1]);
//     if (ent_type (e2) != MATRIX_DENSE_REAL)
//       rerror ("cpf.normal: argument '[mean,std]' has to be MATRIX-DENSE-REAL");
//     x2 = ent_data (e2);
//     nr2 = x2->nrow;
//     nc2 = x2->ncol;
//     if (nc2 != 2)
//       rerror ("cpf.normal: second argument has to be [mean,std]");
// 
//     nr = nr > nr2 ? nr : nr2;
//     w = mdr_Create (nr, nc);
//     // calculate
//     for (i = 0; i < nr; i++)
//     {
//       for (j = 0; j < nc; j++)
//       {
//         id = MIN(i, nr1-1);
//         jd = MIN(j, nc1-1);
//         // x1:
//         if (x1->type == RLAB_TYPE_INT32)
//           x = (double) Mdi0 (x1, id, jd);
//         else
//           x = Mdr0 (x1, id, jd);
//         // x2:
//         id = MIN(i, nr2-1);
//         if (x2->type == RLAB_TYPE_INT32)
//         {
//           a = (double) Mdi0 (x2, id, 0);
//           b = (double) Mdi0 (x2, id, 1);
//         }
//         else
//         {
//           a = Mdr0 (x2, id, 0);
//           b = Mdr0 (x2, id, 1);
//         }
//         Mdr0 (w, i, j) = gsl_cdf_gaussian_P (x-a, b);
//       }
//     }
//   }
//   else if (nargs == 3)
//   {
//     // get x
//     e1 = bltin_get_ent (args[0]);
//     if (ent_type (e1) != MATRIX_DENSE_REAL)
//       rerror ("cpf.normal: argument 'x' has to be MATRIX-DENSE-REAL");
//     x1 = ent_data (e1);
//     nr1 = nr = x1->nrow;
//     nc1 = nc = x1->ncol;
// 
//     //
//     // get a
//     //
//     e2 = bltin_get_ent (args[1]);
//     if (ent_type (e2) != MATRIX_DENSE_REAL)
//       rerror ("cpf.normal: argument 'mean' has to be MATRIX-DENSE-REAL");
//     x2 = ent_data (e2);
//     nr2 = x2->nrow;
//     nc2 = x2->ncol;
//     if (nr2 * nc2 == 0)
//       rerror ("cpf.normal: argument 'mean' has to be MATRIX-DENSE-REAL");
// 
//     // update the size of output matrix
//     nr = nr > nr2 ? nr : nr2;
//     nc = nc > nc2 ? nc : nc2;
// 
//     //
//     // get b
//     //
//     e3 = bltin_get_ent (args[2]);
//     if (ent_type (e3) != MATRIX_DENSE_REAL)
//       rerror ("cpf.normal: argument 'std' has to be MATRIX-DENSE-REAL");
//     x3 = ent_data (e3);
//     nr3 = x3->nrow;
//     nc3 = x3->ncol;
//     if (nr3 * nc3 == 0)
//       rerror ("cpf.normal: argument 'std' has to be MATRIX-DENSE-REAL");
// 
//     // update the size of output matrix to max of all three sizes
//     nr = nr > nr3 ? nr : nr3;
//     nc = nc > nc3 ? nc : nc3;
//     w = mdr_Create (nr, nc);
// 
//     // calculate
//     for (i = 0; i < nr; i++)
//     {
//       for (j = 0; j < nc; j++)
//       {
//         // x:
//         id = MIN(i, nr1-1);
//         jd = MIN(j, nc1-1);
//         if (x1->type == RLAB_TYPE_INT32)
//           x = (double) Mdi0 (x1, id, jd);
//         else
//           x = Mdr0 (x1, id, jd);
//         // a:
//         id = MIN(i, nr2-1);
//         jd = MIN(j, nc2-1);
//         if (x2->type == RLAB_TYPE_INT32)
//           a  = (double) Mdi0 (x2, id, jd);
//         else
//           a  = Mdr0 (x2, id, jd);
//         // b:
//         id = MIN(i, nr3-1);
//         jd = MIN(j, nc3-1);
//         if (x3->type == RLAB_TYPE_INT32)
//           b  = (double) Mdi0 (x3, id, jd);
//         else
//           b  = Mdr0 (x3, id, jd);
// 
//         // w:
//         Mdr0 (w, i, j) = gsl_cdf_gaussian_P (x-a, b);
//       }
//     }
//   }
// 
//   // clean-up
//   ent_Clean (e1);
//   ent_Clean (e2);
//   ent_Clean (e3);
// 
//   return ent_Assign_Rlab_MDR(w);
// }



