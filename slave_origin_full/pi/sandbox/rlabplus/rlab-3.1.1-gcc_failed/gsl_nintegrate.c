// Copyright (C) 2003-2008 Marijan Kostrun
//   part of rlab+4linux project on rlabplus.sourceforge.net
//
// GSL Science Library - numerical integration/QUADPACK
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


// rlab headers, located in variable $RLAB_SDK
#include "complex.h"
#include "rlab.h"
#include "ent.h"
#include "class.h"
#include "symbol.h"
#include "mem.h"
#include "mdr.h"
#include "mdc.h"
#include "mdrf1.h"
#include "mds.h"
#include "list.h"
#include "btree.h"
#include "bltin.h"
#include "util.h"
#include "mathl.h"
#include "function.h"
#include "lp.h"

// gsl headers
# include <gsl/gsl_mode.h>
# include <gsl/gsl_precision.h>
# include <gsl/gsl_math.h>
# include <gsl/gsl_errno.h>
# include <gsl/gsl_machine.h>
# include <gsl/gsl_integration.h>


#include "genzpak.h"
#include "slatec.h"

// standard headers
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>

// naming convention for the solver parameters
#include "rfileio.h"
#include "rlab_solver_parameters_names.h"

static int nintegrate_used=0;
//
// inverse laplace transform
//
static Complex
invlt_f (Complex * x);

static MDC *zmdc;
static Ent *zent, *zpent;
static Ent *zfname_ent;

//
// gsl integrators
//
static double
nint_gslrlab_f (double x, void *dummy);

static MDR *xmdr;
static Ent *xent, *pent;
static Ent *fname_ent;

// **************************************************************
// * Inverse laplace transform                                  *
// **************************************************************
#undef  THIS_SOLVER
#define THIS_SOLVER "invlt"
Ent *
ent_invlt (int nargs, Datum args[])
{
  Ent *e2=0, *e3=0, *e4=0;

  MDR *x=0, *w=0;
  MDC *zwork=0;

  double erel=0.01, sigma=0, ssbar=4, erel2[3], erel2max=0;
  int nmax=1000, ifail=0, ifailmax = 0, ifeval=0;

  double ddummy, valt;
  int idummy, nx, i;

  ListNode *node;
  FILE *rlab_stderr = (!RLAB_STDERR_DS) ? stderr : RLAB_STDERR_DS;

  if (nargs != 3 && nargs != 4)
  {
    fprintf (rlab_stderr,
             THIS_SOLVER ": Inverse laplace transform of a complex function of a complex\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": variable. Format:\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ":   y = invlt(f,/p/,X/,options/),\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": where  f=function(x /,p/), while X = [x1,x2..]\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": is a vector of points at which transform is required. The list\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": options=<<erel;sigma;ssbar;maxfeval>> contains relative error,\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": abscissa of convergence, integration parameter and a maximum number\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": of function evaluations.\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": See manual for more details.\n");
    rerror (THIS_SOLVER ": requires at least 3 arguments");
  }

  //
  // f = f(x/,p/)
  //
  zfname_ent = bltin_get_ent(args[0]);
  if (!isfuncent(zfname_ent))
    rerror (THIS_SOLVER ": " RLAB_ERROR_ARG1_FUNC_VAR "\n");

  //
  // parameter entity
  //
  e2 = bltin_get_ent (args[1]);
  zpent = 0;
  if (ent_type (e2) != UNDEF)
    zpent = ent_Copy (e2);

  //
  // x
  //
  e3 = bltin_get_ent (args[2]);
  if (ent_type (e3) != MATRIX_DENSE_REAL)
    rerror (THIS_SOLVER ": " RLAB_ERROR_ARG3_MDR_VECTOR);

  x = class_matrix_real (e3);
  if (!EQVECT(x))
    rerror (THIS_SOLVER ": " RLAB_ERROR_ARG3_MDR_VECTOR);

  nx = SIZE(x);

  //
  // options for the solver
  //
  if (nargs > 3)
  {
    e4 = bltin_get_ent (args[3]);
    if (ent_type (e4) == BTREE)
    {
      // erel
      node = btree_FindNode (ent_data (e4), RLAB_NAME_GSL_INVLT_EREL);
      if (node != 0)
      {
        ddummy = class_double ( var_ent (node) );
        if (ddummy > 0.0 && ddummy < 1)
          erel = ddummy;
      }
      // sigma
      node = btree_FindNode (ent_data (e4), RLAB_NAME_GSL_INVLT_SIGMA);
      if (node != 0)
      {
        ddummy = class_double ( var_ent (node) );
        if (ddummy > 0.0)
          sigma = ddummy;
      }
      // ssbar
      node = btree_FindNode (ent_data (e4), RLAB_NAME_GSL_INVLT_SSBAR);
      if (node != 0)
      {
        ddummy = class_double ( var_ent (node) );
        if (ddummy > 0.0)
          ssbar = ddummy;

      }
      // max function iterations
      node = btree_FindNode (ent_data (e4), RLAB_NAME_GSL_INVLT_FEVAL);
      if (node != 0)
      {
        idummy = (int) class_double (var_ent (node));
        if (idummy > 0)
          nmax = idummy;
      }
    }
  }

  // Set up ENTITIES for user-function.
  // Complex x
  zmdc = mdc_CreateEmpty(1,1);
  zent = ent_Assign_Rlab_MDC(zmdc);
  ent_IncRef (zent);

  w = mdr_Create(nx, 1);
  zwork = mdc_Create(2, nmax + 1);

  //
  // do the job
  //
  for (i = 0; i < nx ; i++)
  {
    valt = MdrV0(x,i);
    if (valt < 0)
    {
      MdrV0(w,i) = create_nan ();
      continue;
    }
    INVLTF(&erel, &valt, invlt_f, &sigma, &ssbar, &nmax,
           &MdrV0(w,i), erel2, &ifeval, MDCPTR(zwork), &ifail);

    erel2max = erel2[0] > erel2max ? erel2[0] : erel2max;
    ifailmax = ifail < ifailmax ? ifail : ifailmax;

    if (ifail == -1)
      fprintf(rlab_stderr, THIS_SOLVER ": For x=%g the accuracy was not reached "
          "after %i function evaluations!\n", valt, nmax);
    else if (ifail == -2)
      fprintf(rlab_stderr, THIS_SOLVER
          ": For x=%g the algorithm did not appear to be "
          "converging: possibly unsuitable ssbar=%g!\n", valt, ssbar);
  }

  mdc_Destroy(zwork);

  ent_Clean (zfname_ent);

  ent_Clean (e2);
  ent_Clean (zpent);

  MDPTR(zmdc) = 0;
  ent_DecRef  (zent);
  ent_Destroy (zent);

  ent_Clean (e3);
  ent_Clean (e4);

  return ent_Assign_Rlab_MDR (w);
}


static
Complex invlt_f (Complex * x)
{
  Ent *rent;
  Complex zval=0;

  // 'zent'
  MDPTR(zmdc) = x;

  // is there 'pent'
  if (zpent)
    rent = ent_call_rlab_script_2args(zfname_ent, zent, zpent);
  else
    rent = ent_call_rlab_script_1arg (zfname_ent, zent);

  zval = class_complex(rent);
  ent_Clean (rent);
  return zval;
}


// **************************************************************
// * Builtin interface to gsl nintegrate                        *
// **************************************************************
static double trapezoid(MDR *xy, double xlo, double xhi)
{
  double s=0, dx, ya;
  int i, n=MNR(xy);

  if (MNC(xy)!=2)
    return 0;

  if (n <= 1)
    return 0;

  for (i=1; i<n; i++)
  {
    if (mdr0(xy,i-1,0) < xlo)
      continue;
    if (mdr0(xy,i,0) > xhi)
      break;
    dx = mdr0(xy,i,0) - mdr0(xy,i-1,0);
    ya = mdr0(xy,i,1) + mdr0(xy,i-1,1);
    s += 0.5*dx*ya;
  }

  return s;
}

#undef THIS_SOLVER
#define THIS_SOLVER "nintegrate"
Ent *
ent_nintegrate (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *X=0, *eo=0;
  MDR *w=0, *x;
  double eabs=1e-6, erel=0.01, r, *dummy, re, result, xlo=0, xhi=0, ddummy;
  int cstat=0, status=0, max_iter=1000, imethod=3, i, j, xcol, xrow, idummy;
  int check_sort=1;
  ListNode *node;

  FILE *rlab_stderr = (!RLAB_STDERR_DS) ? stderr : RLAB_STDERR_DS;

  if (nintegrate_used)
    rerror (THIS_SOLVER ": " RLAB_ERROR_SOLVER_NO_RECURSIVE_CALLS);

  if (nargs > 7)
  {
    fprintf (rlab_stderr,
             THIS_SOLVER ": Numerical integration of a real function of a real\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": variable. Format:\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ":   y = nintegrate(f,/p/,X/,options/),\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": where  f=function(x /,p/), while X = [x1,x2..]\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": is a row matrix of integration intervals, which.\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": may include +/-inf() as the end points. The list\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": options=<<erel;eabs;maxi;ikey>> contains relative and absolute\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": error, maximum number of subdivisions of an interval\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": to achive the accuracy, while imeth is the index of the\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": method (number of points) used in the integration formula.\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": See manual for more details.\n");
    rerror (THIS_SOLVER ": requires at least 3 arguments");
  }

  //
  // Tabular integration
  //
  e1 = bltin_get_ent (args[0]);
  if(ent_type(e1) == MATRIX_DENSE_REAL)
  {
    MDR *x2 = 0;
    double xmin, xmax;

    // tf = [x_i,y_i]
    MDR *tf = ent_data(e1);
    if (MNC(tf) != 2)
      rerror (THIS_SOLVER ": " RLAB_ERROR_ARG1_TABFUNC "\n");
    if (MNR(tf) < 2)
      rerror (THIS_SOLVER ": " RLAB_ERROR_ARG1_MDR_TWOROWS "\n");
    xmin = Mdr1(tf, 1, 1);
    xmax = Mdr1(tf, MNR(tf), 1);
    if (xmin > xmax)
      rerror (THIS_SOLVER ": improper integration interval\n");

    // check if matrix is sorted w/r to first column
    if (check_sort)
      if (!mdr_IsSorted(tf,0,1))
        goto _exit_nintegrate_table;

    // get the second argument, and if it is a matrix read it in
    if (nargs > 1)
    {
      e2 = bltin_get_ent (args[1]);
      if (ent_type(e2) == MATRIX_DENSE_REAL)
      {
        x2 = class_matrix_real ( e2 );
        if ( (MNR(x2)!=1 && MNC(x2)!=1) &&
             (MNC(x2)!=2) )
          rerror (THIS_SOLVER ": incorrect integration interval\n");
      }
    }

    double *x = &Mdr0(tf,0,0), *y = &Mdr0(tf,0,1);
    int ierr=0, n=MNR(tf);
    double res;

    nintegrate_used=1;
    if (!x2)
    {
      //
      // single integration on [xlo,xhi]
      //
      // slatec tabular integrator
//       DAVINT (x, y, &n, &xmin, &xmax, &res, &ierr);

      //
      res = trapezoid(tf, xmin, xmax);
      w = mdr_CreateScalar (res);
    }
    else if (MNR(x2)==1 || MNC(x2)==1)
    {
      //
      // x2 a vector: perform integrations on [xlo, x2[i]]
      //
      xlo = xmin;
      w   = mdr_Create (MNR(x2), MNC(x2));
      for (i=0; i< MNR(x2)*MNC(x2); i++)
      {
        // foolproof intgration bounds
        xhi = MdrV0(x2, i);
        if (xhi < xmin)
          xhi = xmin;
        else if (xhi > xmax)
          xhi = xmax;

        // slatec tabular integrator
        //DAVINT (x, y, &n, &xlo, &xhi, &res, &ierr);
        //MdrV0 (w,i) = res;
        MdrV0 (w,i) = trapezoid(tf, xlo, xhi);
      }
    }
    else if (MNC(x2)==2)
    {
      //
      // x2 a two-column matrix: perform integrations on [x2[i;1], x2[i;2]]
      //
      w   = mdr_Create (MNR(x2), 1);
      for (i=0; i< MNR(x2); i++)
      {
        // foolproof intgration bounds
        xlo = Mdr0(x2, i, 0);
        if (xlo < xmin)
          xlo = xmin;
        xhi = Mdr0(x2, i, 1);
        if (xhi > xmax)
          xhi = xmax;
        if (xlo == xhi)
          res = 0;
        else if (xlo < xhi)
        {
          // slatec tabular integrator
          DAVINT (x, y, &n, &xlo, &xhi, &res, &ierr);
        }
        else
        {
          // slatec tabular integrator
          DAVINT (x, y, &n, &xhi, &xlo, &res, &ierr);
          res *= (-1);
        }
        MdrV0 (w,i) = res;
      }
    }
    else
    {
      nintegrate_used=0;
      rerror(THIS_SOLVER ": terrible internal error");
    }

_exit_nintegrate_table:

    nintegrate_used=0;
    ent_Clean (e1);
    ent_Clean (e2);
    ent_Clean (eo);
    return ent_Assign_Rlab_MDR (w);
  }


  //
  // Get function ptr
  //
  if (!isfuncent(e1))
    rerror (THIS_SOLVER ": First argument must be function-variable");

  fname_ent = e1;

  gsl_function F;
  F.function = &nint_gslrlab_f;
  F.params = &dummy;

  //
  // parameter entity
  //
  pent = 0;
  e2 = bltin_get_ent (args[1]);
  if (ent_type (e2) != UNDEF)
    pent = ent_Copy (e2);

  // x(lo), lower end of the search interval
  X = bltin_get_ent (args[2]);
  if (ent_type (X) != MATRIX_DENSE_REAL)
    rerror(THIS_SOLVER ": Third argument 'X' must be vector [x0, x1; ..]!");
  x = class_matrix_real (X);
  if (SIZE(x)<2)
    rerror(THIS_SOLVER ": Third argument 'X' must be vector [x0, x1; ..]!");
  if (EQVECT(x))
  {
    xcol = SIZE (x);
    xrow = 1;
  }
  else
  {
    xcol = MNC (x);
    xrow = MNR (x);
  }

  //
  // options for the solver
  //
  if (nargs > 3)
  {
    eo = bltin_get_ent (args[3]);
    if (ent_type (eo) == BTREE)
    {
      // eabs
      node = btree_FindNode (ent_data (eo), RLAB_NAME_GEN_EABS);
      if (node != 0)
      {
        ddummy = class_double (var_ent (node));
        if (ddummy > 0.0)
          eabs = ddummy;
      }
      // erel
      node = btree_FindNode (ent_data (eo), RLAB_NAME_GEN_EREL);
      if (node != 0)
      {
        ddummy = class_double ( var_ent (node) );
        if (ddummy >= 0.0)
          erel = ddummy;
      }
      // method
      node = btree_FindNode (ent_data (eo), RLAB_NAME_GSL_NINT_IKEY);
      if (node != 0)
      {
        idummy = (int) class_double (var_ent (node));
        if (idummy == 0 || idummy == 1 || idummy == 2)
          imethod = idummy;
      }
      // max iterations
      node = btree_FindNode (ent_data (eo), RLAB_NAME_GSL_NINT_MAXITER);
      if (node != 0)
      {
        idummy = (int) class_double (var_ent (node));
        if (idummy > 0)
          max_iter = idummy;
      }
    }

    // clean-up
    ent_Clean (eo);
  }

  // Set up ENTITIES for user-function.
  // x
  xmdr = mdr_CreateEmpty (1,1);
  xent = ent_Assign_Rlab_MDR(xmdr);
  ent_IncRef (xent);

  gsl_integration_workspace *wks = gsl_integration_workspace_alloc (max_iter);
  gsl_set_error_handler_off ();

  w = mdr_Create (xrow, 1);
  mdr_Zero (w);
  nintegrate_used=1;
  for (j=1; j <= xrow; j++)
  {
    //
    // Decide on what integration we do
    //
    for (i = 2; i < xcol; i++)
    {
      // go over intervals, and verify that inf's appear only as the first and last points
      if (Mdr1 (x, j, i) == -create_inf ())
        printf
            (THIS_SOLVER ": -inf() appeared but not as the first point of integration !");
      if (Mdr1 (x, j, i) == create_inf ())
        printf
            (THIS_SOLVER ": +inf() appeared but not as the last point of integration !");
    }

    result = 0;
    for (i = 2; i <= xcol; i++)
    {
      if (Mdr1 (x, j, i - 1) == -create_inf ())
      {
        // is -inf a lower bound of the interval
        if (Mdr1 (x, j, i) == create_inf ())
        {
          // use _qagi for [-inf(),+inf()]
          status = gsl_integration_qagi (&F, eabs, erel, max_iter, wks, &r, &re);
          cstat += status;
          result = result + r;
          continue;
        }
        else
        {
          // use _qagil for [-inf,x[i]]
          xhi = Mdr1 (x, j, i);
          status =
            gsl_integration_qagil (&F, xhi, eabs, erel, max_iter, wks, &r, &re);
          cstat += status;
          result = result + r;
          continue;
        }
      }


      // first point is not -inf
      // is the second point inf
      if (Mdr1 (x, j, i) == create_inf ())
      {
        // use _qagiu for [x[i-1],inf]
        xlo = Mdr1 (x, j, i - 1);
        status =
          gsl_integration_qagiu (&F, xlo, eabs, erel, max_iter, wks, &r, &re);
        cstat += status;
        result = result + r;
        continue;
      }
      // second point must be normal. use qag on it
      xlo = Mdr1 (x, j, i - 1);
      xhi = Mdr1 (x, j, i);
      status =
          gsl_integration_qag (&F, xlo, xhi, eabs, erel, max_iter, imethod, wks, &r,
                           &re);
      if (status != GSL_SUCCESS)
      {
        // if it doesn't work try qags (singularity in the interval
        status =
            gsl_integration_qags (&F, xlo, xhi, eabs, erel, max_iter, wks, &r, &re);
      }
      cstat += status;
      result = result + r;
    }
    if (cstat == GSL_SUCCESS)
      MdrV1(w,j) = result;
    else
      MdrV1(w,j) = create_nan();
    if (status == GSL_EDIVERGE)
      printf
        (THIS_SOLVER ": Integral divergent or slowly convergent. Try nintegrates !\n");
  }

  nintegrate_used=0;

  gsl_integration_workspace_free (wks);

  ent_Clean (pent);
  ent_Clean (e2);

  MDPTR(xmdr) = 0;
  ent_DecRef  (xent);
  ent_Destroy (xent);

  ent_Clean (X);
  ent_Clean (e1);

  return ent_Assign_Rlab_MDR (w);
}

#undef  THIS_SOLVER
#define THIS_SOLVER "cauchypval"
Ent *
ent_nintegrate_qawc (int nargs, Datum args[])
{
  Ent *e2, *X, *eo;
  MDR *w, *x;
  double eabs = 1e-6, erel = 0.01, xlo, xhi, s, r, aer, *dummy, ddummy;
  int status=0, max_iter = 1000, idummy;
  gsl_function F;
  F.function = &nint_gslrlab_f;
  F.params = &dummy;
  ListNode *node;

  FILE *rlab_stderr = (!RLAB_STDERR_DS) ? stderr : RLAB_STDERR_DS;

  if (nargs < 3 && nargs > 4)
  {
    fprintf (rlab_stderr,
             THIS_SOLVER ": Cauchy principal value integral.\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": Format:\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ":   y = cauchypval(f,/p/,X/,options/),\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": where  f=function(x /,p/), while X = [a,c,b]\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": contains a,b as the endpoints and  c  as a singularity\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": at which the principal value is calculated. The list\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": options=<<erel;eabs;maxi>> contain relative and absolute error,\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": and the maximum number of subdivisions of an interval to achive\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": the accuracy. See manual for more details.\n");
    rerror (THIS_SOLVER ": requires at 3 or 4 arguments");
  }

  //
  // f = f(x/,p/)
  //
  fname_ent = bltin_get_ent(args[0]);
  if (!isfuncent(fname_ent))
    rerror (THIS_SOLVER ": " RLAB_ERROR_ARG1_FUNC_VAR "\n");

  //
  // parameter entity
  //
  e2 = bltin_get_ent (args[1]);
  pent = 0;
  if (ent_type (e2) != UNDEF)
    pent = ent_Copy (e2);

  // X = [xlo, s, xhi]
  X = bltin_get_ent (args[2]);
  if (ent_type (X) != MATRIX_DENSE_REAL)
    rerror (THIS_SOLVER ": Third argument 'x' has to be real vector [xlo,sing,xhi]!");
  x = class_matrix_real (X);
  if (SIZE(x)!=3)
    rerror (THIS_SOLVER ": Third argument 'x' has to be real vector [xlo,sing,xhi]!");

  xlo = MdrV0 (x, 0);
  s   = MdrV0 (x, 1);
  xhi = MdrV0 (x, 2);
  if (nargs > 3)
  {
    eo = bltin_get_ent (args[3]);
    if (ent_type (eo) == BTREE)
    {
      // eabs
      node = btree_FindNode (ent_data (eo), RLAB_NAME_GSL_NINT_EABS);
      if (node != 0)
      {
        ddummy = class_double (var_ent (node));
        if (ddummy > 0.0)
          eabs = ddummy;
      }
      // erel
      node = btree_FindNode (ent_data (eo), RLAB_NAME_GSL_NINT_EREL);
      if (node != 0)
      {
        ddummy = class_double ( var_ent (node) );
        if (ddummy >= 0.0)
          erel = ddummy;
      }
      // max iterations
      node = btree_FindNode (ent_data (eo), RLAB_NAME_GSL_NINT_MAXITER);
      if (node != 0)
      {
        idummy = (int) class_double (var_ent (node));
        if (idummy > 0)
          max_iter = idummy;
      }
    }
    if ( ent_type (eo) != UNDEF ) ent_Clean (eo);
  }

  // Set up ENTITIES for user-function.
  // x
  xmdr = mdr_CreateEmpty (1,1);
  xent = ent_Assign_Rlab_MDR(xmdr);
  ent_IncRef (xent);

  gsl_integration_workspace *wks = gsl_integration_workspace_alloc (max_iter);

  gsl_set_error_handler_off ();

  status =
    gsl_integration_qawc (&F, xlo, xhi, s, eabs, erel, max_iter, wks, &r, &aer);
  if (status == GSL_SUCCESS)
    w = mdr_CreateScalar (r);
  else
    w = mdr_Create (0, 0);

  if (status == GSL_EDIVERGE)
    printf (THIS_SOLVER ": Integral divergent or slowly convergent!\n");
  if (status == GSL_ESING)
    printf (THIS_SOLVER ": Non-integrable singularity found!\n");
  if (status == GSL_EROUND)
    printf (THIS_SOLVER ": Cannot reach requested tolerance!\n");
  if (status == GSL_EMAXITER)
    printf (THIS_SOLVER ": The max number of subdivision exceeded!\n");

  gsl_integration_workspace_free (wks);

  ent_Clean(fname_ent);

  ent_Clean (pent);
  ent_Clean (e2);

  MDPTR(xmdr) = 0;
  ent_DecRef (xent);
  ent_Destroy (xent);

  ent_Clean (X);

  return ent_Assign_Rlab_MDR (w);
}

#undef  THIS_SOLVER
#define THIS_SOLVER "nintegrates"

Ent *
ent_nintegrate_qagp (int nargs, Datum args[])
{
  Ent *e2, *X, *eo;
  MDR *w, *x;
  double eabs = 1e-6, erel = 0.01, r, *dummy, re, ddummy;
  int status=0, max_iter = 1000, xcol, idummy;
  gsl_function F;
    F.function = &nint_gslrlab_f;
    F.params = &dummy;
  ListNode *node;

  FILE *rlab_stderr = (!RLAB_STDERR_DS) ? stderr : RLAB_STDERR_DS;

  if (nargs < 3)
  {
    fprintf (rlab_stderr,
             THIS_SOLVER ": Numerical integration of a real function of a real variable\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": with singularities in the integration interval endpoints.\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": Format:\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ":   y = nintegrate(f,/p/,X/,options/),\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": where  f=function(x /,p/), while X = [x1,x2..]\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": is a row matrix of integration intervals, which.\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": may NOT include +/-inf() as the end points. The list\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": options=<<erel;eabs;maxi>> contains relative and absolute\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": error, and a maximum number of subdivisions of an \n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": interval to achive the accuracy.\n");
    rerror (THIS_SOLVER ": requires at least 3 arguments");
  }

  //
  // f = f(x/,p/)
  //
  fname_ent = bltin_get_ent(args[0]);
  if (!isfuncent(fname_ent))
    rerror (THIS_SOLVER ": " RLAB_ERROR_ARG1_FUNC_VAR "\n");

  //
  // parameter entity
  //
  e2 = bltin_get_ent (args[1]);
  pent = 0;
  if (ent_type (e2) != UNDEF)
    pent = ent_Copy (e2);

  // x(lo), lower end of the search interval
  X = bltin_get_ent (args[2]);
  if (ent_type (X) != MATRIX_DENSE_REAL)
    rerror(THIS_SOLVER ": Third argument 'X' must be vector [x0, x1; ..]!");
  x = class_matrix_real (X);
  if (!x)
    rerror(THIS_SOLVER ": Third argument 'X' must be vector [x0, x1; ..]!");
  xcol = MNC (x);
  if (xcol< 2)
    rerror(THIS_SOLVER ": Third argument 'X' must be vector [x0, x1; ..]!");

  if (nargs > 3)
  {
    eo = bltin_get_ent (args[3]);
    if (ent_type (eo) == BTREE)
    {
      // eabs
      node = btree_FindNode (ent_data (eo), RLAB_NAME_GSL_NINT_EABS);
      if (node != 0)
      {
        ddummy = class_double (var_ent (node));
        if (ddummy > 0.0)
          eabs = ddummy;
      }
      // erel
      node = btree_FindNode (ent_data (eo), RLAB_NAME_GSL_NINT_EREL);
      if (node != 0)
      {
        ddummy = class_double ( var_ent (node) );
        if (ddummy >= 0.0)
          erel = ddummy;
      }
      // max iterations
      node = btree_FindNode (ent_data (eo), RLAB_NAME_GSL_NINT_MAXITER);
      if (node != 0)
      {
        idummy = (int) class_double (var_ent (node));
        if (idummy > 0)
          max_iter = idummy;
      }
    }
    if ( ent_type (eo) != UNDEF ) ent_Clean (eo);
  }

  // Set up ENTITIES for user-function.
  // x
  xmdr = mdr_CreateEmpty (1,1);
  xent = ent_Assign_Rlab_MDR(xmdr);
  ent_IncRef (xent);

  gsl_integration_workspace *wks = gsl_integration_workspace_alloc (max_iter);

  gsl_set_error_handler_off ();

  status =
    gsl_integration_qagp (&F, MDRPTR(x), MNC (x), eabs, erel, max_iter, wks, &r,
                          &re);
  if (status == GSL_SUCCESS)
    w = mdr_CreateScalar (r);
  else
    w = mdr_Create (0, 0);
  if (status == GSL_EDIVERGE)
    printf (THIS_SOLVER ": Integral divergent or slowly convergent!\n");

  gsl_integration_workspace_free (wks);

  ent_Clean (fname_ent);

  ent_Clean (pent);
  ent_Clean (e2);

  MDPTR(xmdr) = 0;
  ent_DecRef (xent);
  ent_Destroy (xent);

  ent_Clean (X);

  return ent_Assign_Rlab_MDR (w);
}


#undef  THIS_SOLVER
#define THIS_SOLVER "nintqaws"
Ent *
ent_nintegrate_qaws (int nargs, Datum args[])
{
  Ent *e2=0, *XLO=0, *XHI=0, *ALPHA=0, *BETA=0, *MU=0, *NU=0, *eo;
  MDR *w;
  double eabs = 1e-6, erel = 0.01, xlo, xhi, alpha, beta, r, *dummy, aer, ddummy;
  int status=0, max_iter = 1000, mu, nu, idummy;
  gsl_function F;
  F.function = &nint_gslrlab_f;
  F.params = &dummy;
  ListNode *node;

  FILE *rlab_stderr = (!RLAB_STDERR_DS) ? stderr : RLAB_STDERR_DS;

  if (nargs < 7)
    rerror (THIS_SOLVER ": requires at least 7 arguments");

  //
  // f = f(x/,p/)
  //
  fname_ent = bltin_get_ent(args[0]);
  if (!isfuncent(fname_ent))
    rerror (THIS_SOLVER ": " RLAB_ERROR_ARG1_FUNC_VAR "\n");

  //
  // parameter entity
  //
  e2 = bltin_get_ent (args[1]);
  pent = 0;
  if (ent_type (e2) != UNDEF)
    pent = ent_Copy (e2);

  // x(lo), lower end of the search interval
  XLO = bltin_get_ent (args[2]);
  if (ent_type (XLO) != MATRIX_DENSE_REAL)
    rerror(THIS_SOLVER ": Third argument 'XLO' must be scalar!");
  xlo = class_double (XLO);

  // x(hi), high end of the search interval
  XHI = bltin_get_ent (args[3]);
  if (ent_type (XHI) != MATRIX_DENSE_REAL)
    rerror(THIS_SOLVER ": Fourth argument 'XHI' must be scalar!");
  xhi = class_double (XHI);

  // alpha
  ALPHA = bltin_get_ent (args[4]);
  if (ent_type (ALPHA) != MATRIX_DENSE_REAL)
    rerror(THIS_SOLVER ": Fifth argument 'A' must be scalar!");
  alpha = class_double (ALPHA);

  // beta
  BETA = bltin_get_ent (args[5]);
  if (ent_type (BETA) != MATRIX_DENSE_REAL)
    rerror(THIS_SOLVER ": Sixth argument 'B' must be scalar!");
  beta = class_double (BETA);

  // mu
  MU = bltin_get_ent (args[6]);
  if (ent_type (MU) == MATRIX_DENSE_REAL)
    rerror(THIS_SOLVER ": Seventh argument 'MU' must be scalar!");
  mu = (int) class_double (MU);

  // nu
  NU = bltin_get_ent (args[7]);
  if (ent_type (NU) == MATRIX_DENSE_REAL)
    rerror(THIS_SOLVER ": Eigth argument 'NU' must be scalar!");
  nu = (int) class_double (NU);

  if (nargs > 7)
  {
    eo = bltin_get_ent (args[8]);
    if (ent_type (eo) == BTREE)
    {
      // eabs
      node = btree_FindNode (ent_data (eo), RLAB_NAME_GSL_NINT_EABS);
      if (node != 0)
      {
        ddummy = class_double (var_ent (node));
        if (ddummy > 0.0)
          eabs = ddummy;
      }
      // erel
      node = btree_FindNode (ent_data (eo), RLAB_NAME_GSL_NINT_EREL);
      if (node != 0)
      {
        ddummy = class_double ( var_ent (node) );
        if (ddummy >= 0.0)
          erel = ddummy;
      }
      // max iterations
      node = btree_FindNode (ent_data (eo), RLAB_NAME_GSL_NINT_MAXITER);
      if (node != 0)
      {
        idummy = (int) class_double (var_ent (node));
        if (idummy > 0)
          max_iter = idummy;
      }
    }
    if ( ent_type (eo) != UNDEF ) ent_Clean (eo);
  }

  // Set up ENTITIES for user-function.
  // x
  xmdr = mdr_CreateEmpty (1,1);
  xent = ent_Assign_Rlab_MDR(xmdr);
  ent_IncRef (xent);

  gsl_set_error_handler_off ();

  gsl_integration_workspace *wks = gsl_integration_workspace_alloc (max_iter);
  gsl_integration_qaws_table *t =
    gsl_integration_qaws_table_alloc (alpha, beta, mu, nu);

  status =
    gsl_integration_qaws (&F, xlo, xhi, t, eabs, erel, max_iter, wks, &r, &aer);
  if (status == GSL_SUCCESS)
    w = mdr_CreateScalar (r);
  else
    w = mdr_Create (0, 0);
  if (status == GSL_EDIVERGE)
    fprintf (rlab_stderr, THIS_SOLVER ": Integral divergent or slowly convergent!\n");

  gsl_integration_workspace_free (wks);

  gsl_integration_qaws_table_free (t);

  ent_Clean (fname_ent);

  ent_Clean (pent);
  ent_Clean (e2);

  MDPTR(xmdr) = 0;
  ent_DecRef (xent);
  ent_Destroy (xent);

  ent_Clean (XLO);
  ent_Clean (XHI);
  ent_Clean (ALPHA);
  ent_Clean (BETA);
  ent_Clean (MU);
  ent_Clean (NU);

  return ent_Assign_Rlab_MDR (w);
}


#undef  THIS_SOLVER
#define THIS_SOLVER "nintsin"
Ent *
ent_nintegrate_qawo_sin (int nargs, Datum args[])
{
  Ent *e2=0, *X=0, *OMEGA=0, *eo=0;
  MDR *w, *x;
  double eabs = 1e-6, erel = 0.01, xlo, r, *dummy, re, L, omega, result, ddummy;
  int status=0, maxn = 100, maxiter = 2000, xcol, i, cstat = 0, idummy;
  gsl_function F;
  F.function = &nint_gslrlab_f;
  F.params = &dummy;
  ListNode *node;

  FILE *rlab_stderr = (!RLAB_STDERR_DS) ? stderr : RLAB_STDERR_DS;

  if (nargs < 4)
  {
    fprintf (rlab_stderr,
             THIS_SOLVER ": Sinus-OMEGA transform of a real function of a real variable\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": Format:\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ":   y = nintsin(f,/p/,X,Omega/,options/),\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": where  f=function(x /,p/), while X = [x1,x2..]\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": is a row matrix of integration intervals, which.\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": my include +inf() as the right end point. The list\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": options=<<erel;eabs;maxi;maxn>> contains relative and absolute\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": error, a maximum number of subdivisions of an interval and a \n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": maximum number of intervals achive the accuracy. Omega is \n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": the angular frequency of the  sin()  function entering the integral.\n");
    rerror (THIS_SOLVER ": requires at least 4 arguments !");
  }
  //
  // f = f(x/,p/)
  //
  fname_ent = bltin_get_ent(args[0]);
  if (!isfuncent(fname_ent))
    rerror (THIS_SOLVER ": " RLAB_ERROR_ARG1_FUNC_VAR "\n");

  //
  // parameter entity
  //
  e2 = bltin_get_ent (args[1]);
  pent = 0;
  if (ent_type (e2) != UNDEF)
    pent = ent_Copy (e2);

  // X = [xlo, s, xhi]
  X = bltin_get_ent (args[2]);
  if (ent_type (X) != MATRIX_DENSE_REAL)
    rerror (THIS_SOLVER ": Third argument 'x' has to be real matrix [x0,x1;..]!");
  x = class_matrix_real (X);
  if (!EQVECT(x))
    rerror (THIS_SOLVER ": Third argument 'x' has to be real vector [x0,x1;..]!");
  xcol = SIZE(x);
  if (xcol < 2)
    rerror (THIS_SOLVER ": Third argument 'x' has to be real matrix [x0,x1;..]!");

  // omega,
  OMEGA = bltin_get_ent (args[3]);
  if (ent_type (OMEGA) != MATRIX_DENSE_REAL)
    rerror (THIS_SOLVER ": Fourth argument 'omega' must be scalar !");
  omega = class_double (OMEGA);

  if (nargs > 4)
  {
    eo = bltin_get_ent (args[4]);
    if (ent_type (eo) == BTREE)
    {
      // eabs
      node = btree_FindNode (ent_data (eo), RLAB_NAME_GSL_NINT_EABS);
      if (node != 0)
      {
        ddummy = class_double (var_ent (node));
        if (ddummy > 0.0)
          eabs = ddummy;
      }
      // erel
      node = btree_FindNode (ent_data (eo), RLAB_NAME_GSL_NINT_EREL);
      if (node != 0)
      {
        ddummy = class_double ( var_ent (node) );
        if (ddummy >= 0.0)
          erel = ddummy;
      }
      // max iterations
      node = btree_FindNode (ent_data (eo), RLAB_NAME_GSL_NINT_MAXITER);
      if (node != 0)
      {
        idummy = (int) class_double (var_ent (node));
        if (idummy > 0)
          maxiter = idummy;
      }
      // max iterations
      node = btree_FindNode (ent_data (eo), RLAB_NAME_GSL_NINT_MAXINT);
      if (node != 0)
      {
        idummy = (int) class_double (var_ent (node));
        if (idummy > 0)
          maxn = idummy;
      }
    }
    if ( ent_type (eo) != UNDEF ) ent_Clean (eo);
  }

  // Set up ENTITIES for user-function.
  // x
  xmdr = mdr_CreateEmpty (1,1);
  xent = ent_Assign_Rlab_MDR(xmdr);
  ent_IncRef (xent);

  gsl_set_error_handler_off ();

  gsl_integration_workspace *wks = gsl_integration_workspace_alloc (maxiter);
  gsl_integration_workspace *cycwks = gsl_integration_workspace_alloc (maxiter);

  //
  // check for +-inf()'s
  //
  for (i = 2; i <= xcol; i++)
  {
    // go over intervals, and verify that inf's appears as the last point
    if (MdrV1 (x, i - 1) == -create_inf ())
      rerror (THIS_SOLVER ": cannot deal with -inf() as an integration bound !");
    if (MdrV1 (x, i - 1) == create_inf ())
      rerror
        (THIS_SOLVER ": +inf() appeared but not as the last point of integration !");
  }
  //
  // start
  //
  result = 0;
  for (i = 2; i <= xcol; i++)
  {
    // first point cannot be -inf
    xlo = MdrV1 (x, i - 1);
    // is the second point inf
    if (MdrV1 (x, i) == create_inf ())
      L = 3.0 * M_PI * 0.5 * omega;
    else
      L = MdrV1 (x, i) - xlo;
    gsl_integration_qawo_table *t
      = gsl_integration_qawo_table_alloc (omega, L, GSL_INTEG_SINE, maxn);
    if (MdrV1 (x, i) == create_inf ())
    {
      status =
        gsl_integration_qawf (&F, xlo, eabs, maxiter, wks, cycwks, t, &r, &re);
    }
    else
    {
      status =
        gsl_integration_qawo (&F, xlo, eabs, erel, maxiter, wks, t, &r, &re);
    }
    gsl_integration_qawo_table_free (t);
    if (status == GSL_ETABLE)
      rerror
        (THIS_SOLVER ": number of levels 'maxn' insufficient for requested accuracy !");
    cstat += status;
    result = result + r;
  }

  if (cstat == GSL_SUCCESS)
    w = mdr_CreateScalar (result);
  else
    w = mdr_Create (0, 0);
  if (status == GSL_EDIVERGE)
    fprintf (rlab_stderr, THIS_SOLVER ": Integral divergent or slowly convergent !\n");

  gsl_integration_workspace_free (wks);

  gsl_integration_workspace_free (cycwks);

  ent_Clean (fname_ent);

  MDPTR(xmdr) = 0;
  ent_DecRef (xent);
  ent_Destroy (xent);

  ent_Clean (pent);
  ent_Clean (e2);

  ent_Clean (X);
  ent_Clean (OMEGA);

  return ent_Assign_Rlab_MDR (w);
}


#undef  THIS_SOLVER
#define THIS_SOLVER "nintcos"
Ent *
ent_nintegrate_qawo_cos (int nargs, Datum args[])
{
  Ent *e2, *X, *eo, *OMEGA;
  MDR *w, *x;
  double eabs = 1e-6, erel = 0.01, xlo, r, *dummy, re, L, omega, result, ddummy;
  int status=0, maxn = 20, maxiter = 1000, xcol, i, cstat = 0, idummy;
  gsl_function F;
  F.function = &nint_gslrlab_f;
  F.params = &dummy;
  ListNode *node;

  FILE *rlab_stderr = (!RLAB_STDERR_DS) ? stderr : RLAB_STDERR_DS;

  if (nargs < 4)
  {
    fprintf (rlab_stderr,
             THIS_SOLVER ": Cosinus-OMEGA transform of a real function of a real variable\n");
    fprintf (rlab_stderr, THIS_SOLVER ": Format:\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ":   y = nintcos(f,/p/,X,Omega/,options/),\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": where  f=function(x /,p/), while X = [x1,x2..]\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": is a row matrix of integration intervals, which.\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": my include +inf() as the right end point. The list\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": options=<<erel;eabs;maxi;maxn>> contains relative and absolute\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": error, a maximum number of subdivisions of an interval and a \n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": maximum number of intervals achive the accuracy. Omega is \n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": the angular frequency of the  cos()  function entering the integral.\n");
    rerror (THIS_SOLVER ": requires at least 4 arguments !");
  }

  //
  // f = f(x/,p/)
  //
  fname_ent = bltin_get_ent(args[0]);
  if (!isfuncent(fname_ent))
    rerror (THIS_SOLVER ": " RLAB_ERROR_ARG1_FUNC_VAR "\n");

  //
  // parameter entity
  //
  e2 = bltin_get_ent (args[1]);
  pent = 0;
  if (ent_type (e2) != UNDEF)
    pent = ent_Copy (e2);

  // X = [xlo, s, xhi]
  X = bltin_get_ent (args[2]);
  if (ent_type (X) != MATRIX_DENSE_REAL)
    rerror (THIS_SOLVER ": Third argument 'x' has to be real matrix [x0,x1;..]!");
  x = class_matrix_real (X);
  if (!x)
    rerror (THIS_SOLVER ": Third argument 'x' has to be real matrix [x0,x1;..]!");
  xcol = SIZE (x);
  if (xcol < 2)
    rerror (THIS_SOLVER ": Third argument 'x' has to be real matrix [x0,x1;..]!");

  // omega,
  OMEGA = bltin_get_ent (args[3]);
  if (ent_type (OMEGA) != MATRIX_DENSE_REAL)
    rerror (THIS_SOLVER ": Fourth argument 'omega' must be scalar !");
  omega = class_double (OMEGA);

  if (nargs > 4)
  {
    eo = bltin_get_ent (args[4]);
    if (ent_type (eo) == BTREE)
    {
      // eabs
      node = btree_FindNode (ent_data (eo), RLAB_NAME_GSL_NINT_EABS);
      if (node != 0)
      {
        ddummy = class_double (var_ent (node));
        if (ddummy > 0.0)
          eabs = ddummy;
      }
      // erel
      node = btree_FindNode (ent_data (eo), RLAB_NAME_GSL_NINT_EREL);
      if (node != 0)
      {
        ddummy = class_double ( var_ent (node) );
        if (ddummy >= 0.0)
          erel = ddummy;
      }
      // max iterations
      node = btree_FindNode (ent_data (eo), RLAB_NAME_GSL_NINT_MAXITER);
      if (node != 0)
      {
        idummy = (int) class_double (var_ent (node));
        if (idummy > 0)
          maxiter = idummy;
      }
      // max iterations
      node = btree_FindNode (ent_data (eo), RLAB_NAME_GSL_NINT_MAXINT);
      if (node != 0)
      {
        idummy = (int) class_double (var_ent (node));
        if (idummy > 0)
          maxn = idummy;
      }
    }
    if ( ent_type (eo) != UNDEF ) ent_Clean (eo);
  }

  // Set up ENTITIES for user-function.
  // x
  xmdr = mdr_CreateEmpty (1,1);
  xent = ent_Assign_Rlab_MDR(xmdr);
  ent_IncRef (xent);

  gsl_set_error_handler_off ();

  gsl_integration_workspace *wks = gsl_integration_workspace_alloc (maxiter);
  gsl_integration_workspace *cycwks = gsl_integration_workspace_alloc (maxiter);

  //
  // check for +-inf()'s
  //
  for (i = 2; i <= xcol; i++)
  {
    // go over intervals, and verify that inf's appears as the last point
    if (MdrV1 (x, i - 1) == -create_inf ())
      rerror (THIS_SOLVER ": cannot deal with -inf() as an integration bound !");
    if (MdrV1 (x, i - 1) == create_inf ())
      rerror
        (THIS_SOLVER ": +inf() appeared but not as the last point of integration !");
  }
  //
  // start
  //
  result = 0;
  for (i = 2; i <= xcol; i++)
  {
    // first point cannot be -inf
    xlo = MdrV1 (x, i - 1);
    // is the second point inf
    if (MdrV1 (x, i) == create_inf ())
      L = 3.0 * M_PI * 0.5 * omega;
    else
      L = MdrV1 (x, i) - xlo;
    gsl_integration_qawo_table *t
      = gsl_integration_qawo_table_alloc (omega, L, GSL_INTEG_COSINE, maxn);
    if (MdrV1 (x, i) == create_inf ())
    {
      status = gsl_integration_qawf
        (&F, xlo, eabs, maxiter, wks, cycwks, t, &r, &re);
    }
    else
    {
      status =
        gsl_integration_qawo (&F, xlo, eabs, erel, maxiter, wks, t, &r, &re);
    }
    gsl_integration_qawo_table_free (t);
    if (status == GSL_ETABLE)
      rerror
        (THIS_SOLVER ": number of levels 'maxn' insufficient for requested accuracy !");
    cstat += status;
    result = result + r;
  }

  if (cstat == GSL_SUCCESS)
    w = mdr_CreateScalar (result);
  else
    w = mdr_Create (0, 0);
  if (status == GSL_EDIVERGE)
    fprintf (rlab_stderr, THIS_SOLVER ": Integral divergent or slowly convergent !\n");

  gsl_integration_workspace_free (wks);

  gsl_integration_workspace_free (cycwks);

  ent_Clean (fname_ent);

  if (pent)
  {
    ent_Clean (pent);
    ent_Clean (e2);
  }

  MDPTR(xmdr) = 0;
  ent_DecRef (xent);
  ent_Destroy (xent);

  ent_Clean (X);
  ent_Clean (OMEGA);

  return ent_Assign_Rlab_MDR (w);
}


//
// The interface to the user-specified function.
//  MDR params is passed by default, no need for it to be the argument of the function
//
static double
nint_gslrlab_f (double x, void *dummy)
{
  Ent *rent;
  double rval;

  // 'xent'
  MDPTR(xmdr) = &x;

  // is there 'pent'
  if (pent)
  {
    rent = ent_call_rlab_script_2args(fname_ent, xent, pent);
  }
  else
  {
    rent = ent_call_rlab_script_1arg (fname_ent, xent);
  }

  rval = class_double(rent);
  ent_Clean (rent);
  return rval;
}

//
// genzpak: integration in multidimensions based on work of Alan Genz
//
#undef  THIS_SOLVER
#define THIS_SOLVER "nintmd"
Ent *
ent_nintegrate_genz_hypercube (int nargs, Datum args[])
{
  Ent *e2=0, *e3=0, *eo=0;

  MDR *work, *iwork, *x3, *a, *b;

  int i, ndim, nfun = 1, minpts = 1, maxpts = 100000, irestar =0, idummy;
  int neval, ifail, numsms=10000, imethod=0, ikey = 4;
  double eabs = 1e-3, erel = 1e-2, result, abserr;
  double ddummy;

  ListNode *node;

  FILE *rlab_stderr = (!RLAB_STDERR_DS) ? stderr : RLAB_STDERR_DS;

  // parameters and initial message
  if (nargs < 3 || nargs > 4)
  {
    fprintf (rlab_stderr,
             THIS_SOLVER ": Integration of a real scalar function in multidimensions.\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": Format:\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ":   nintmd(f,/p/,X/,options/),\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": where f = function(x/,p/), 'p' is its parameter array,\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": X=[xlo,xhi;..] is a matrix containing the integration intervals.\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": Check the manual for other options.\n");
    rerror ("requires three or four arguments");
  }

  //
  // Get function ptr
  //
  fname_ent = bltin_get_ent(args[0]);
  if (!isfuncent(fname_ent))
    rerror (THIS_SOLVER ": " RLAB_ERROR_ARG1_FUNC_VAR "\n");

  //
  // Get function parameters
  //
  pent = 0;
  e2 = bltin_get_ent (args[1]);
  if (ent_type (e2) != UNDEF)
    pent = ent_Copy(e2);

  //
  // X=[xlo, xhi;..] stacked integration intervals for each dimension
  //
  e3  = bltin_get_ent (args[2]);
  if (ent_type(e3) != MATRIX_DENSE_REAL)
    rerror (THIS_SOLVER ": Integration intervals need to be a real matrix!");
  x3 = class_matrix_real (e3);
  if (MNC(x3)!=2)
    rerror (THIS_SOLVER ": Integration intervals need to be a two-column real matrix!");
  ndim = MNR (x3);
  if (ndim > 20)
    rerror (THIS_SOLVER ": Maximum dimension of integration hypercube is 20!");
  a = mdr_PartitionCol(x3,1);
  b = mdr_PartitionCol(x3,2);
  for (i=0; i<ndim; i++)
    if (MdrV0(a,i) >= MdrV0(b,i))
  {
    fprintf(rlab_stderr, THIS_SOLVER ": bound %i is improper, lo > hi!", i+1);
    rerror (THIS_SOLVER ": Improper Integration interval!");
  }

  if (nargs > 3)
  {
    eo = bltin_get_ent (args[3]);
    if (ent_type (eo) == BTREE)
    {
      // imethod
      node = btree_FindNode (ent_data (eo), RLAB_NAME_GENZ_IMETHOD);
      if (node != 0)
      {
        idummy = (int) class_double ( var_ent (node) );
        if (idummy == 0 || idummy == 1)
          imethod = idummy;
      }
      // ikey
      node = btree_FindNode (ent_data (eo), RLAB_NAME_GENZ_IKEY);
      if (node != 0)
      {
        idummy = (int) class_double ( var_ent (node) );
        if (idummy >= 0 && idummy <= 4)
          ikey = idummy;
      }
      // maxi
      node = btree_FindNode (ent_data (eo), RLAB_NAME_GENZ_MAXI);
      if (node != 0)
      {
        idummy = (int) class_double ( var_ent (node) );
        if (idummy > 0 )
          maxpts = idummy;
      }
      // mini
      node = btree_FindNode (ent_data (eo), RLAB_NAME_GENZ_MINI);
      if (node != 0)
      {
        idummy = (int) class_double ( var_ent (node) );
        if (idummy > 0 )
          minpts = idummy;
      }
      // eabs
      node = btree_FindNode (ent_data (eo), RLAB_NAME_GENZ_EABS);
      if (node != 0)
      {
        ddummy = class_double (var_ent (node));
        if (ddummy > 0.0 && ddummy < 1.0)
          eabs = ddummy;
      }
      // erel
      node = btree_FindNode (ent_data (eo), RLAB_NAME_GENZ_EREL);
      if (node != 0)
      {
        ddummy = class_double (var_ent (node));
        if (ddummy > 0.0 && ddummy < 1.0)
          erel = ddummy;
      }
    }
    if ( ent_type (eo) != UNDEF ) ent_Clean (eo);
  }

  //
  // Set up ENTITIES for rlab function being integrated
  // x
  xmdr = mdr_CreateEmpty (ndim, 1);
  xent = ent_Assign_Rlab_MDR(xmdr);
  ent_IncRef (xent);

  if (imethod == 1)
  {
    //
    // patpack
    //

    work = mdr_Create(2*nfun + (nfun+1)*numsms,1);

    // integrate
    PATSYM(&ndim, MDRPTR(a), MDRPTR(b), &nfun, &minpts, &maxpts,
            PATFUN, &eabs, &erel, &irestar, &result,
            &abserr, &neval, &ifail, MDRPTR(work) );

    if (ifail == 1)
    {
      fprintf(rlab_stderr,
              THIS_SOLVER ": Number of integrations was too small for required accuracy!");
      fprintf(rlab_stderr,
              THIS_SOLVER ": Increase 'maxi' for more reliable result!");
    }
    // clean-up
    mdr_Destroy (work);
  }
  else if (imethod == 2)
  {
    //
    // decuhr: integral with singularity at bottom corner
    //

    int nw, num=0, isingul, ilogf=0, iemax=5, iwrksub=1000;
    double alpha;

    isingul = ndim;
    alpha   = -isingul;

    if ((ikey == 0 || ikey == 1) && ndim == 2)
      num = 65;
    else if ((ikey == 0 || ikey == 2) && (ndim == 3))
      num = 127;
    else if ((ikey == 0 && ndim > 3) || (ikey == 3))
      num = 1 + 2*ndim + 6*ndim*ndim + 4*ndim*(ndim-1)*(ndim-2)/3 + pow(2,ndim);
    else if (ikey == 4)
      num = 1 + 2*ndim*(ndim+2) + pow(2,ndim);

    maxpts = MAX(maxpts, 3*num);

    nw = 20+iwrksub*(ndim+2)*2+iemax+(3*iwrksub+9+iemax)+(iemax+1)*(iemax+1)
        + 3*ndim;
    work  = mdr_Create(nw,1);
    iwork = mdi_Create(2*iwrksub+ndim,1);

    // integrate
    DECUHR(&ndim, &nfun, MDRPTR(a), MDRPTR(b), &minpts, &maxpts,
           PATFUN, &isingul, &alpha, &ilogf, &eabs, &erel, &ikey,
           &iwrksub, &nw, &irestar, &iemax, &result, &abserr,
           &neval, &ifail, MDRPTR(work), MDIPTR(iwork));

    if (ifail >= 1)
    {
      fprintf(rlab_stderr,
              THIS_SOLVER ": Number of integrations was too small for required accuracy!");
      fprintf(rlab_stderr,
              THIS_SOLVER ": Increase 'maxi' for more reliable result!");
    }

    // clean-up
    mdr_Destroy (work);
    mdr_Destroy (iwork);
  }
  else
  {
    //
    // dcuhre
    //

    int nw, num=0, maxsub;

    if ((ikey == 0 || ikey == 1) && ndim == 2)
      num = 65;
    else if ((ikey == 0 || ikey == 2) && (ndim == 3))
      num = 127;
    else if ((ikey == 0 && ndim == 3) || (ikey == 3))
      num = 1 + 4*2*ndim + 2*ndim*(ndim-1) + 4*ndim*(ndim-1) +
          4*ndim*(ndim-1)*(ndim-2)/3 + pow(2,ndim);
    else if (ikey == 4)
      num = 1 + 3*2*ndim + 2*ndim*(ndim-1) + pow(2,ndim);

    maxpts = MAX(maxpts, 3*num);

    maxsub = (maxpts-num)/(2*num) + 1;
    nw = maxsub*(2*ndim+4) + 17 + 1;
    work = mdr_Create(nw,1);

    // integrate
    DCUHRE(&ndim, &nfun, MDRPTR(a), MDRPTR(b), &minpts, &maxpts,
            PATFUN, &eabs, &erel, &ikey, &nw, &irestar, &result,
            &abserr, &neval, &ifail, MDRPTR(work));

    if (ifail == 1)
    {
      fprintf(rlab_stderr,
              THIS_SOLVER ": Number of integrations was too small for required accuracy!");
      fprintf(rlab_stderr,
              THIS_SOLVER ": Increase 'maxi' for more reliable result!");
    }
    // clean-up
    mdr_Destroy (work);
  }

  //
  // clean-up
  //
  mdr_Destroy (a);
  mdr_Destroy (b);

  ent_Clean (fname_ent);

  ent_Clean (pent);

  ent_Clean (e2);

  MDPTR(xmdr) = 0;
  ent_DecRef (xent);
  ent_Destroy (xent);

  ent_Clean (e3);

  return ent_Create_Rlab_Double (result);
}


#undef  THIS_SOLVER
#define THIS_SOLVER "nintsimplex"
Ent *
ent_nintegrate_genz_simplex (int nargs, Datum args[])
{
  Ent *e2=0, *e3=0, *eo=0;

  MDR *work, *x3;

  int i, j, nd, nf = 1, mincls = 1, maxcls = 100000, irestar =0, idummy;
  int ifuncls, inform, ikey = 3, iwrklen, iwrksbs, isbrgns, irulcls;
  double eabs = 1e-3, erel = 1e-2, result, abserr;
  double ddummy;

  ListNode *node;
  FILE *rlab_stderr = (!RLAB_STDERR_DS) ? stderr : RLAB_STDERR_DS;

  // parameters and initial message
  if (nargs < 3)
  {
    fprintf (rlab_stderr,
             THIS_SOLVER ": Integration of a real scalar function in multidimensions.\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": Format:\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ":   nintsimplex(f,/p/,X/,options/),\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": where f = function(x/,p/), 'p' is its parameter array,\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": 'X' is a matrix containing simplex vertices row-wise.\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": Check the manual for other options.\n");
    rerror ("requires at least 4 arguments");
  }

  //
  // Get function ptr
  //
  fname_ent = bltin_get_ent(args[0]);
  if (!isfuncent(fname_ent))
    rerror (THIS_SOLVER ": " RLAB_ERROR_ARG1_FUNC_VAR "\n");

  //
  // Get function parameters
  //
  pent = 0;
  e2 = bltin_get_ent (args[1]);
  if (ent_type (e2) != UNDEF)
    pent = ent_Copy(e2);

  //
  // X=[xlo, xhi;..] stacked integration intervals for each dimension
  //
  e3  = bltin_get_ent (args[2]);
  if (ent_type(e3) != MATRIX_DENSE_REAL)
    rerror (THIS_SOLVER ": Integration intervals need to be a real matrix!\n");
  x3 = class_matrix_real (e3);
  nd = MNC (x3);
  if (MNR(x3) != nd+1)
    rerror (THIS_SOLVER ": Matrix of simplex vertices has to be (n+1)-by-n, with n=dim(x)!\n");

  if (nargs > 3)
  {
    eo = bltin_get_ent (args[3]);
    if (ent_type (eo) == BTREE)
    {
      // maxi
      node = btree_FindNode (ent_data (eo), RLAB_NAME_GENZ_MAXI);
      if (node != 0)
      {
        idummy = (int) class_double ( var_ent (node) );
        if (idummy > 0 )
          maxcls = idummy;
      }
      // mini
      node = btree_FindNode (ent_data (eo), RLAB_NAME_GENZ_MINI);
      if (node != 0)
      {
        idummy = (int) class_double ( var_ent (node) );
        if (idummy > 0 )
          mincls = idummy;
      }
      // ikey
      node = btree_FindNode (ent_data (eo), RLAB_NAME_GENZ_IKEY);
      if (node != 0)
      {
        idummy = (int) class_double ( var_ent (node) );
        if (idummy >= 0 && idummy < 5)
          ikey = idummy;
      }
      // eabs
      node = btree_FindNode (ent_data (eo), RLAB_NAME_GENZ_EABS);
      if (node != 0)
      {
        ddummy = class_double (var_ent (node));
        if (ddummy > 0.0 && ddummy < 1.0)
          eabs = ddummy;
      }
      // erel
      node = btree_FindNode (ent_data (eo), RLAB_NAME_GENZ_EREL);
      if (node != 0)
      {
        ddummy = class_double (var_ent (node));
        if (ddummy > 0.0 && ddummy < 1.0)
          erel = ddummy;
      }
    }
    if ( ent_type (eo) != UNDEF ) ent_Clean (eo);
  }

  //
  // Set up ENTITIES for rlab function being integrated
  // x
  xmdr = mdr_CreateEmpty (nd, 1);
  xent = ent_Assign_Rlab_MDR(xmdr);
  ent_IncRef (xent);

  // integration over one simplex
  isbrgns = 1;
  //
  if (ikey == 0)
    irulcls = (nd+4)*(nd+3)*(nd+2)/6 + (nd+2)*(nd+1);
  else if (ikey == 1)
    irulcls = 2*nd+3;
  else if (ikey == 2)
    irulcls = (nd+3)*(nd+2)/2 + 2*(nd+1);
  else if (ikey == 3)
    irulcls = (nd+4)*(nd+3)*(nd+2)/6 + (nd+2)*(nd+1);
  else
    irulcls = (nd+5)*(nd+4)*(nd+3)*(nd+2)/24 + 5*(nd+2)*(nd+1)/2;
  //
  iwrksbs = isbrgns + 3*( maxcls/irulcls - isbrgns*(1-irestar))/4;
  iwrklen = iwrksbs*( nd*(nd+1) + 2*nf + 3 ) + (nd+1)*(nd+2) + 7*nf;

  work = mdr_Create(iwrklen,1);
  mdr_Zero(work);
  for (i=0; i<MNR(x3); i++)
    for (j=0; j<nd; j++)
      MdrV0(work,j+i*nd) = Mdr0(x3,i,j);

  // integrate
  SMPINT(&nd, &nf, &mincls, &maxcls, PATFUN,
          &eabs, &erel, &ikey, &isbrgns, &iwrklen, MDRPTR(work),
         &irestar, &result, &abserr, &ifuncls, &inform );

  if (inform == 1)
  {
    fprintf(rlab_stderr,
            THIS_SOLVER ": Number of integrations was too small for required accuracy!\n");
    fprintf(rlab_stderr,
            THIS_SOLVER ": Increase 'maxi' for more reliable result!\n");
  }

  //
  // clean-up
  //
  mdr_Destroy (work);

  ent_Clean (fname_ent);

  ent_Clean (pent);
  ent_Clean (e2);

  MDPTR(xmdr) = 0;
  ent_DecRef (xent);
  ent_Destroy (xent);

  ent_Clean (e3);

  return ent_Create_Rlab_Double (result);
}

int
PATFUN (int *ndim, double *x, int *numfun, double *fval)
{
  Ent *rent;

  // 'xent'
  MDPTR(xmdr) = (void *) x;

  // is there 'pent'
  if (pent)
    rent = ent_call_rlab_script_2args(fname_ent, xent, pent);
  else
    rent = ent_call_rlab_script_1arg (fname_ent, xent);

  *fval = class_double(rent);
  ent_Clean (rent);
  return 1;
}

