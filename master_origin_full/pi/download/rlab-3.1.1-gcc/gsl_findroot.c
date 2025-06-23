// Copyright (C) 2003-2008 Marijan Kostrun
//   part of rlab+4linux project on rlabplus.sourceforge.net
//
// GSL Science Library - root finding in one dimension
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

#include "hompack.h"

// gsl headers
#include <gsl/gsl_mode.h>
#include <gsl/gsl_precision.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_machine.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_multiroots.h>

// standard headers
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <time.h>

// naming convention for the solver parameters
#include "rfileio.h"
#include "rlab_solver_parameters_names.h"

//
// findroot: link to gsl functions
//
static double
findroot_gslrlab_f (double x, void *dummy);

// static double
// findroot_gslrlab_fLR (double x, void *dummy);

static double
findroot_gslrlab_df (double x, void *dummy);

static void
findroot_gslrlab_fdf (double x, void *dummy, double *y, double *dy);

//
// findroot: other variables
//
static MDR *xmdr=0;
static Ent *xent=0;
static Ent *pent=0;
static Ent *fname=0;
static Ent *dfname=0;

static double findroot_f_const = 0;




// **************************************************************
// * Builtin interface to gsl 1d root finding routines          *
// **************************************************************
#undef  THIS_SOLVER
#define THIS_SOLVER "findroot"
Ent * ent_findroot (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0;
  findroot_f_const = 0;

  FILE *rlab_stderr = (!RLAB_STDERR_DS) ? stderr : RLAB_STDERR_DS;

  // set these to 0 just to avoid any wierdness
  xent=0;
  pent=0;
  fname=0;
  dfname=0;

  if (nargs!=1 && nargs!=2 && nargs < 3)
  {
    fprintf (rlab_stderr,
             THIS_SOLVER ": Finds a root of a function with and without derivative,\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": where function can be given in tabulated form:\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": Format:\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": (1) r = findroot( [x,y] /,y0/  ),\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": (2) r = findroot(f,/p/,X /,options/),\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": (3) r = findroot(f,df,/p/,x0 /,options/),\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": where f = function(x /,p/), X=[x1,x2] is the bracketing\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": interval, or x0 is the initial point. The list \n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": options=<<erel;eabs;imethod;maxi;rhs>> contains relative and\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": absolute error, choice of method, where in case (1)\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": imethod=0 for bisection, 1 for regula falsi, and 2 for\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": Brent-Dekker, or in case (2), 0 for Newton's, 1 for\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": Secant, and 2 for Steffenson's method. In the latter case\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": df=function(x/,p) has to be specified where df=df/dx.\n");
    rerror ("requires at three to five arguments !");
  }

  //
  // Get function ptr
  //
  e1 = bltin_get_ent(args[0]);
  if (isfuncent(e1))
    fname = e1;
  else if (ent_type (e1) == MATRIX_DENSE_REAL)
  {
    double yoff = 0.0;

    if (nargs==2)
    {
      e2 = bltin_get_ent(args[1]);
      if (ent_type (e2) == MATRIX_DENSE_REAL)
      {
        yoff = class_double(e2);
      }
    }

    // special case:
    // first argument is a two column real matrix which zeros we are
    // supposed to find!
    MDR *x1 = ent_data(e1);
    if (MNC(x1) != 2)
      rerror ("findroot: first argument has to be a function or a two-column real matrix!");

    MDR *w= mdr_FindRootsFromInterpolationTable(x1, yoff);

    ent_Clean (e1);
    ent_Clean (e2);

    return ent_Assign_Rlab_MDR(w);
  }
  else
    rerror ("findroot: first argument has to be a function or a two-column real matrix!");

  // set the right hand side of the problem being solved
  // f(..) = RHS
  findroot_f_const = 0;

  //
  // Get second datum
  //
  e2 = bltin_get_ent(args[1]);
  if (isfuncent(e2))
  {
    //
    // algorithms with derivatives
    //
    Ent *e3=0, *X0=0, *e4=0;
    MDR *w;
    double eabs = 1e-4, erel = 0, x0, x1, x2, *dummy, ddummy;
    int status, iter = 0, max_iter = 1000, imethod = 0, idummy;
    gsl_function_fdf FDF;
    FDF.f = &findroot_gslrlab_f;
    FDF.df = &findroot_gslrlab_df;
    FDF.fdf = &findroot_gslrlab_fdf;
    FDF.params = &dummy;
    const gsl_root_fdfsolver_type *T;
    gsl_root_fdfsolver *s;
    ListNode *node;

    //
    // df
    //
    dfname = e2;

    //
    // parameter entity
    //
    e3 = bltin_get_ent (args[2]);
    pent = 0;
    if (ent_type (e3) != UNDEF)
      pent = ent_Copy (e3);

    //
    // x0, initial guess for the root
    //
    X0 = bltin_get_ent (args[3]);
    if (ent_type (X0) != MATRIX_DENSE_REAL)
      rerror ("findroot: argument 'x0' must be real vector!");
    x0 = class_double (X0);
    if(!x0)
      rerror ("findroot: argument 'x0' must be real vector!");

    //
    // options for the solver
    //
    if (nargs >= 5)
    {
      e4 = bltin_get_ent (args[4]);
      if (ent_type (e4) == BTREE)
      {
        // eabs
        node = btree_FindNode (ent_data (e4), RLAB_NAME_GSL_ROOT1_EABS);
        if (node != 0)
        {
          ddummy = class_double (var_ent (node));
          if (ddummy > 0.0)
            eabs = ddummy;
        }
        // erel
        node = btree_FindNode (ent_data (e4), RLAB_NAME_GSL_ROOT1_EREL);
        if (node != 0)
        {
          ddummy = class_double ( var_ent (node) );
          if (ddummy >= 0.0)
            erel = ddummy;
        }
        // method
        node = btree_FindNode (ent_data (e4), RLAB_NAME_GSL_ROOT1_METHOD);
        if (node != 0)
        {
          idummy = (int) class_double (var_ent (node));
          if (idummy == 0 || idummy == 1 || idummy == 2)
            imethod = idummy;
        }
        // max iterations
        node = btree_FindNode (ent_data (e4), RLAB_NAME_GSL_ROOT1_MAXI);
        if (node != 0)
        {
          idummy = (int) class_double (var_ent (node));
          if (idummy > 0)
            max_iter = idummy;
        }
        // right hand side
        node = btree_FindNode (ent_data (e4), RLAB_NAME_GSL_ROOT1_RHS);
        if (node != 0)
        {
          findroot_f_const = class_double ( var_ent (node) );
        }
      }
      ent_Clean (e4);
    }
    //
    // Set up ENTITIES for user-function.
    // x
    xmdr = mdr_CreateEmpty (1,1);
    xent = ent_Assign_Rlab_MDR (xmdr);
    ent_IncRef (xent);

    x1 = x0;
    x2 = x0;

    gsl_set_error_handler_off ();
    switch (imethod)
    {
    case 0:
      T = gsl_root_fdfsolver_newton;
      break;
    case 1:
      T = gsl_root_fdfsolver_secant;
      break;
    case 2:
      T = gsl_root_fdfsolver_steffenson;
      break;
    default:
      T = gsl_root_fdfsolver_newton;
    }
    s = gsl_root_fdfsolver_alloc (T);
    gsl_root_fdfsolver_set (s, &FDF, x1);
    do
    {
      iter++;
      status = gsl_root_fdfsolver_iterate (s);
      x1 = x2;
      x2 = gsl_root_fdfsolver_root (s);
      status = gsl_root_test_delta (x2, x1, eabs, erel);
    }
    while (status == GSL_CONTINUE && iter < max_iter);
    if (status == GSL_SUCCESS)
      w = mdr_CreateScalar (x2);
    else
      w = mdr_Create (0, 0);
    gsl_root_fdfsolver_free (s);

    if (ent_type(e3) != UNDEF)
    {
      ent_Clean (pent);
      ent_Clean (e3);
    }

    MDPTR(xmdr) = 0;
    ent_DecRef (xent);
    ent_Destroy (xent);

    ent_Clean (X0);

    ent_Clean (e1);
    ent_Clean (e2);

    return ent_Assign_Rlab_MDR(w);
  }
  //
  // second argument is not function. assume parameter array (may be missing).
  // use _b class of methods
  //
  Ent *X=0, *e4=0;
  MDR *w=0, *x=0;
  double eabs = 1e-4, erel = 0, r, *dummy, xlo, xhi, ddummy;
  int status, iter = 0, max_iter = 1000, imethod = 0, idummy;
  // declare function for a bracketing algorithm
  gsl_function F;
  F.function = &findroot_gslrlab_f;
  F.params = &dummy;
  const gsl_root_fsolver_type *T;
  gsl_root_fsolver *s;
  ListNode *node;

  // params, a MDR parameter array
  pent = 0;
  if (ent_type(e2) != UNDEF)
    pent = ent_Copy (e2);

  // x(lo), lower end of the search interval
  if (nargs < 3)
    rerror ("findroot: Third argument 'x' must be real vector [xlo,xhi] !");
  X = bltin_get_ent (args[2]);
  if (ent_type (X) != MATRIX_DENSE_REAL)
    rerror ("findroot: Third argument 'x' must be real vector [xlo,xhi] !");
  x = class_matrix_real (X);
  if (!x)
    rerror ("findroot: Third argument 'x' must be real vector [xlo,xhi] !");
  if (MNR(x)*MNC(x)!=2)
    rerror ("findroot: Third argument 'x' must be real vector [xlo,xhi] !");
  xlo = MdrV1(x,1);
  xhi = MdrV1(x,2);

  //
  // options for the solver
  //
  if (nargs >= 4)
  {
    e4 = bltin_get_ent (args[3]);
    if (ent_type (e4) == BTREE)
    {
      // eabs
      node = btree_FindNode (ent_data (e4), RLAB_NAME_GSL_ROOT1_EABS);
      if (node != 0)
      {
        ddummy = class_double (var_ent (node));
        if (ddummy > 0.0)
          eabs = ddummy;
      }
      // erel
      node = btree_FindNode (ent_data (e4), RLAB_NAME_GSL_ROOT1_EREL);
      if (node != 0)
      {
        ddummy = class_double ( var_ent (node) );
        if (ddummy >= 0.0)
          erel = ddummy;
      }
      // method
      node = btree_FindNode (ent_data (e4), RLAB_NAME_GSL_ROOT1_METHOD);
      if (node != 0)
      {
        idummy = (int) class_double (var_ent (node));
        if (idummy == 0 || idummy == 1 || idummy == 2)
          imethod = idummy;
      }
      // max iterations
      node = btree_FindNode (ent_data (e4), RLAB_NAME_GSL_ROOT1_MAXI);
      if (node != 0)
      {
        idummy = (int) class_double (var_ent (node));
        if (idummy > 0)
          max_iter = idummy;
      }
      // right hand side
      node = btree_FindNode (ent_data (e4), RLAB_NAME_GSL_ROOT1_RHS);
      if (node != 0)
      {
        findroot_f_const = class_double ( var_ent (node) );
      }
    }
    if ( ent_type (e4) != UNDEF ) ent_Clean (e4);
  }

  //
  // Set up ENTITIES for user-function.
  // x
  xmdr = mdr_CreateEmpty (1, 1);
  xent = ent_Assign_Rlab_MDR (xmdr);
  ent_IncRef (xent);

  // continue
  gsl_set_error_handler_off ();
  switch (imethod)
  {
  case 0:
    T = gsl_root_fsolver_bisection;
    break;
  case 1:
    T = gsl_root_fsolver_falsepos;
    break;
  default:
    T = gsl_root_fsolver_brent;
  }
  s = gsl_root_fsolver_alloc (T);
  gsl_root_fsolver_set (s, &F, xlo, xhi);
  do
  {
    iter++;
    status = gsl_root_fsolver_iterate (s);
    r = gsl_root_fsolver_root (s);
    xlo = gsl_root_fsolver_x_lower (s);
    xhi = gsl_root_fsolver_x_upper (s);
    status = gsl_root_test_interval (xlo, xhi, eabs, erel);
  }
  while (status == GSL_CONTINUE && iter < max_iter);
  if (status == GSL_SUCCESS)
    w = mdr_CreateScalar (r);
  else
    w = mdr_Create (0, 0);

  ent_Clean (e1);
  if (pent)
  {
    ent_Clean (pent);
    ent_Clean (e2);
  }

  MDPTR(xmdr) = 0;
  ent_DecRef (xent);
  ent_Destroy (xent);

  ent_Clean (X);
  gsl_root_fsolver_free (s);

  return ent_Assign_Rlab_MDR(w);
}


//
// The interface to the user-specified function.
//  MDR params is passed by default, no need for it to be the argument of the function
//
static double
findroot_gslrlab_f (double x, void *dummy)
{
  Ent *rent = 0;
  double retv;

  MDPTR (xmdr) = &x;

  if (pent)
    rent = ent_call_rlab_script_2args(fname, xent, pent);
  else
    rent = ent_call_rlab_script_1arg (fname, xent);

  if (ent_type(rent)!=MATRIX_DENSE_REAL)
    rerror (THIS_SOLVER ": " RLAB_ERROR_RHS_FUNC_MUST_RETURN_MDR);

  retv = class_double (rent);

  ent_Clean (rent);

  // findroot_f_const = 0    for findroot(s)
  // findroot_f_const = 0.5  for bpoint on {0,1}
  return (retv - findroot_f_const);
}

static double
findroot_gslrlab_df (double x, void *dummy)
{
  Ent *rent = 0;
  double retv;

  MDPTR (xmdr) = &x;

  if (pent)
    rent = ent_call_rlab_script_2args(dfname, xent, pent);
  else
    rent = ent_call_rlab_script_1arg (dfname, xent);

  if (ent_type(rent)!=MATRIX_DENSE_REAL)
    rerror (THIS_SOLVER ": " RLAB_ERROR_RHS_FUNC_MUST_RETURN_MDR);

  retv = class_double (rent);
  ent_Clean (rent);
  return retv;
}

static void
findroot_gslrlab_fdf (double x, void *dummy, double *y, double *dy)
{
  *y = findroot_gslrlab_f (x, &dummy);
  *dy = findroot_gslrlab_df (x, &dummy);
}


