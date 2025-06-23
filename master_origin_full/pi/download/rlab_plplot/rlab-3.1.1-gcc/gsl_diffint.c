// Copyright (C) 2003-2008 Marijan Kostrun
//   part of rlab+4linux project on rlabplus.sourceforge.net
//
// GSL Science Library
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

#include "rfileio.h"
#include "rlab_solver_parameters_names.h"
#include "rlab_macros.h"
#include "rlab_macros_code.h"

// gsl headers
// shared object
#include <gsl/gsl_mode.h>
#include <gsl/gsl_precision.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sum.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_machine.h>
#include <gsl/gsl_chebyshev.h>
#include <gsl/gsl_deriv.h>

// standard headers
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>

static double
diffint_gslrlab_func (double x, void *dummy);

static MDR *xmdr;
static Ent *xent, *pent;
static Ent *fname;


// ****************************************************************
// * Builtin interface to gsl first derivative of a real function *
// ****************************************************************
#undef  THIS_SOLVER
#define THIS_SOLVER "ndiff"
Ent *
ent_ndiff (int nargs, Datum args[])
{
  Ent *e2=0, *X0=0, *XLO=0, *XHI=0, *H=0;
  MDR *w, *x0;
  double eabs=1e-6, xlo, xhi, *dummy=0, h=1e-6;
  int x0r, i, j;

  FILE *rlab_stderr = (!RLAB_STDERR_DS) ? stderr : RLAB_STDERR_DS;

  gsl_function F;
  F.function = &diffint_gslrlab_func;
  F.params = &dummy;

  if (nargs == 0 || nargs > 6)
  {
    fprintf (rlab_stderr,
             THIS_SOLVER ": Numerical differentiation of a real function of a real\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": scalar variable, or of a function given in a table form.\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": Format:\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": (1) ndiff(f,/p/,x/,xlo,xhi,h/),\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": where  f=function(x/,p/){..}, 'p' is the parameter array,\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": and  x=[x1;..;xN]  is a column-vector of points where the numerical\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": differential is being sought. 'xlo' and 'xhi', if given determine\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": the type of difference formula used at the endpoints of array 'x':\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": x1=xlo calls for forward difference, xN=xhi calls for backward,\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": otherwise central difference formula is used.\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": (2) ndiff(F),\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": where  F = [x1,y1;..]  is a table of function values. For endpoint\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": forward/backward difference formula is used, while for midpoints\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": the central difference formula is used.\n");
    rerror (THIS_SOLVER ": requires one or three arguments");
  }

  // Get first argument and check if it is a function ptr
  if (nargs >= 3)
  {
    //
    // f = f(x/,p/)
    //
    fname = bltin_get_ent(args[0]);
    if (!isfuncent(fname))
      rerror (THIS_SOLVER ": " RLAB_ERROR_ARG1_FUNC_VAR "\n");

    //
    // parameter entity
    //
    pent = 0;
    e2 = bltin_get_ent (args[1]);
    if (ent_type (e2) != UNDEF)
      pent = ent_Copy (e2);

    // x0
    X0 = bltin_get_ent (args[2]);
    x0 = class_matrix_double (X0);
    if (!md_is_vector(x0))
      rerror (THIS_SOLVER ": " RLAB_ERROR_ARG3_MDR_VECTOR_DOUBLE "\n");
    x0r = SIZE(x0);

    // xlo
    xlo = create_nan ();
    if (nargs > 3)
    {
      XLO = bltin_get_ent (args[3]);
      if (ent_type (XLO) == MATRIX_DENSE_REAL)
        xlo = class_double (XLO);
    }

    // xhi
    xhi = create_nan ();
    if (nargs > 4)
    {
      XHI = bltin_get_ent (args[4]);
      if (ent_type (XHI) == MATRIX_DENSE_REAL)
        xhi = class_double (XHI);
    }

    // h
    if (nargs > 5)
    {
      H = bltin_get_ent (args[5]);
      if (ent_type (H) == MATRIX_DENSE_REAL)
        h = class_double (H);
      h = h  > 0 ? h : -h;
      h = h == 0 ? 1e-6 : h;
    }
    //
    // Set up ENTITIES for user-function.
    // xmdr -> zero size MDR
    xmdr = mdr_CreateEmpty (1,1);
    xent = ent_Assign_Rlab_MDR (xmdr);
    ent_IncRef (xent);

    // params

    gsl_set_error_handler_off ();

    w = mdr_Create (x0r, 1);

    for (i=0; i<x0r; i++)
    {
      if (xlo == MdrV0 (x0, i))
        gsl_deriv_forward  (&F, MdrV0 (x0, i), h, &MdrV0(w,i), &eabs);
      else if (xhi == MdrV0 (x0, i))
        gsl_deriv_backward (&F, MdrV0 (x0, i), h, &MdrV0(w,i), &eabs);
      else
        gsl_deriv_central  (&F, MdrV0 (x0, i), h, &MdrV0(w,i), &eabs);
    }

    // xmdr -> zero size MDR
    MDPTR(xmdr) = 0;
    ent_DecRef (xent);
    ent_Destroy (xent);

    ent_Clean (e2);
    ent_Clean (pent);
    ent_Clean (X0);
    ent_Clean (XLO);
    ent_Clean (XHI);
    ent_Clean (H);
    ent_Clean (fname);

    return ent_Assign_Rlab_MDR(w);
  }

  //
  // single argument: two-column matrix
  //
  X0 = bltin_get_ent (args[0]);
  if (ent_type (X0) != MATRIX_DENSE_REAL)
    rerror ("ndiff: first argument has to be a real [x,y] matrix!");
  x0 = class_matrix_real (X0);

  x0r = MNR (x0);
  if (x0->ncol != 2)
    rerror ("ndiff: first argument has to be a real [x,y] matrix!");
  if (x0r <= 2)
    rerror ("ndiff: first argument has to be a real [x,y] matrix!");
  w = mdr_Create (x0r, 1);

  // forward difference at xlo
  MdrV1 (w, 1) =
    (Mdr1 (x0, 2, 2) - Mdr1 (x0, 1, 2)) / (Mdr1 (x0, 2, 1) - Mdr1 (x0, 1, 1));

  // backward difference at xhi
  MdrV1 (w, x0r) =
    (Mdr1 (x0, x0r, 2) - Mdr1 (x0, x0r - 1, 2)) / (Mdr1 (x0, x0r, 1) -
                                                     Mdr1 (x0, x0r - 1, 1));
  // central difference for intermediate values
  for (j = 2; j < x0r; j++)
  {
    Mdr1 (w, j, 1) =
    0.5 * (Mdr1 (x0, j + 1, 2) - Mdr1 (x0, j, 2)) / (Mdr1 (x0, j + 1, 1) -
    Mdr1 (x0, j, 1));
    Mdr1 (w, j, 1) +=
    0.5 * (Mdr1 (x0, j, 2) - Mdr1 (x0, j - 1, 2)) / (Mdr1 (x0, j, 1) -
    Mdr1 (x0, j - 1, 1));
  }

  ent_Clean (X0);

  return ent_Assign_Rlab_MDR(w);
}

// ****************************************************************
// * Builtin interface to gsl tchebyshev toolkit                  *
// ****************************************************************
#undef  THIS_SOLVER
#define THIS_SOLVER "tchebyfit"
Ent *
ent_tchebyshev_init (int nargs, Datum args[])
{
  Ent *e2=0, *e3=0, *N;
  MDR *w, *x=0;
  double xlo=-1, xhi=1, *dummy;
  int i, n;

//   FILE *rlab_stderr = (!RLAB_STDERR_DS) ? stderr : RLAB_STDERR_DS;

  gsl_function F;
  F.function = &diffint_gslrlab_func;
  F.params = &dummy;
  gsl_cheb_series *cs;

  if (nargs < 3)
    rerror (THIS_SOLVER ": requires three arguments");

  //
  // f = f(x/,p/)
  //
  fname = bltin_get_ent(args[0]);
  if (!isfuncent(fname))
    rerror (THIS_SOLVER ": " RLAB_ERROR_ARG1_FUNC_VAR "\n");

  //
  // parameter entity
  //
  e2 = bltin_get_ent (args[1]);
  pent = 0;
  if (ent_type (e2) != UNDEF)
    pent = ent_Copy (e2);

  //
  // n
  //
  N = bltin_get_ent (args[2]);
  if (ent_type(N) != MATRIX_DENSE_REAL)
    rerror (THIS_SOLVER ": 'n' has to be a real number");
  n = (int) class_double (N);

  // interval over which the fit is sought
  if (nargs > 3)
  {
    e3 = bltin_get_ent (args[3]);
    if (ent_type (e3) == MATRIX_DENSE_REAL)
      x = class_matrix_real (e3);
    else
      rerror (THIS_SOLVER ": domain '[xlo,xhi]' has to be real vector of size 2");
    if (x->nrow * x->ncol != 2)
      rerror (THIS_SOLVER ": domain '[xlo,xhi]' has to be real vector of size 2");
  }

  if (x)
  {
    xlo = MdrV0(x,0);
    xhi = MdrV0(x,1);
  }


  // Set up ENTITIES for user-function.
  // x
  xmdr = mdr_CreateEmpty (1,1);
  xent = ent_Assign_Rlab_MDR (xmdr);
  ent_IncRef (xent);

  gsl_set_error_handler_off ();

  w = mdr_Create (n + 1, 1);

  cs = gsl_cheb_alloc (n);
  gsl_cheb_init (cs, &F, xlo, xhi);
  for (i = 0; i <= n; i++)
    MdrV1 (w, i + 1) = cs->c[i];

  // clean-up
  gsl_cheb_free (cs);

  ent_Clean (e2);
  ent_Clean (pent);

  MDPTR(xmdr) = 0;
  ent_DecRef (xent);
  ent_Destroy (xent);

  ent_Clean (N);

  //
  // write result to output list
  //
  Btree *bw = btree_Create ();
  // coefficients of the interpolation
  install (bw, RLAB_NAME_DIFFINT_TCH_COEF, ent_Assign_Rlab_MDR(w));

  //  interval
  if (x)
    install (bw, RLAB_NAME_DIFFINT_TCH_INT, ent_Assign_Rlab_MDR(mdr_Copy(x)));
  else
  {
    x = mdr_Create(1,2);
    MdrV0(x,0) = xlo;
    MdrV0(x,1) = xhi;
    install (bw, RLAB_NAME_DIFFINT_TCH_INT, ent_Assign_Rlab_MDR(x));
  }
  ent_Clean(e3);

  // degree of the fit
  install (bw, RLAB_NAME_DIFFINT_TCH_DEG, ent_Create_Rlab_Double(n));
  // label
  install (bw, RLAB_NAME_DIFFINT_TCH_NAME, ent_Create_Rlab_String("tchebyshev"));

  return ent_Assign_Rlab_BTREE (bw);
}

Ent *
ent_tchebyshev_eval (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0;
  int n, i, xr, cr;
  MDR *x, *xrange, *c, *w;
  gsl_cheb_series *cs;

  if (nargs < 2)
    rerror ("tchebyval: two arguments required!");

  e1 = bltin_get_ent (args[0]);
  x = class_matrix_real (e1);
  if ((x->nrow != 1 && x->ncol != 1) ||(x->nrow * x->ncol == 0))
    rerror ("tchebyval: 'x' has to be column-vector of non-zero size!");
  xr = x->nrow * x->ncol;

  e2 = bltin_get_ent (args[1]);
  if (ent_type (e2) != BTREE)
    rerror ("tchebyval: second argument has to be list <<coef;interval;max_order;name>>!");

  // <<coef;interval;max_order;name>>
  ListNode *node;
  node = btree_FindNode (ent_data (e2), RLAB_NAME_DIFFINT_TCH_COEF);
  if (!node)
    rerror ("tchebyval: entry '" RLAB_NAME_DIFFINT_TCH_COEF "' has to be vector of non-zero size!");
  c = class_matrix_real (var_ent (node));

  node = btree_FindNode (ent_data (e2), RLAB_NAME_DIFFINT_TCH_INT);
  if (!node)
    rerror ("tchebyval: entry '" RLAB_NAME_DIFFINT_TCH_INT "' has to be vector [xlo,xhi]!");
  xrange = class_matrix_real (var_ent (node));

  if ((c->nrow != 1 && c->ncol != 1) ||(c->nrow * c->ncol == 0))
    rerror ("tchebyval: 'c' has to be column-vector of non-zero size!");
  cr = c->nrow * c->ncol;
  n  = cr - 1;

  cs = gsl_cheb_alloc (n);
  cs->order = n;
  cs->a = MdrV0(xrange,0);
  cs->b = MdrV0(xrange,1);

  for (i = 0; i <= n; i++)
    cs->c[i] = Mdr1 (c, i + 1, 1);

  w = mdr_Create (xr, 1);
  for (i = 1; i <= xr; i++)
    MdrV1 (w, i) = gsl_cheb_eval (cs, MdrV1 (x, i));

  gsl_cheb_free (cs);
  ent_Clean (e1);
  ent_Clean (e2);

  return ent_Assign_Rlab_MDR(w);
}

Ent *
ent_tchebyshev_diff (int nargs, Datum args[])
{
  Ent *e1=0;
  int n, i, cc;
  MDR *c=0, *w, *xrange=0;
  gsl_cheb_series *cs, *dcs;

  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) != BTREE)
    rerror ("tchebyval: second argument has to be list <<coef;interval;max_order;name>>!");

  // <<coef;interval;max_order;name>>
  ListNode *node;
  node = btree_FindNode (ent_data (e1), RLAB_NAME_DIFFINT_TCH_COEF);
  if (!node)
    rerror ("tchebyval: entry '" RLAB_NAME_DIFFINT_TCH_COEF "' has to be vector of non-zero size!");
  c = class_matrix_real (var_ent (node));

  node = btree_FindNode (ent_data (e1), RLAB_NAME_DIFFINT_TCH_INT);
  if (!node)
    rerror ("tchebyval: entry '" RLAB_NAME_DIFFINT_TCH_INT "' has to be vector [xlo,xhi]!");
  xrange = class_matrix_real (var_ent (node));

  cc = MNC (c);
  if (cc > 1)
    rerror ("tchebyder: 'c' is a column-vector!");

  n = MNR (c) - 1;

  cs = gsl_cheb_alloc (n);
  cs->order = n;
  cs->a = MdrV0(xrange,0);
  cs->b = MdrV0(xrange,1);
  dcs = gsl_cheb_alloc (n);
  dcs->order = n;
  dcs->a = MdrV0(xrange,0);
  dcs->b = MdrV0(xrange,1);

  for (i = 0; i <= n; i++)
    cs->c[i] = Mdr1 (c, i + 1, 1);

  w = mdr_Create (n + 1, 1);

  gsl_cheb_calc_deriv (dcs, cs);

  for (i = 0; i <= n; i++)
    Mdr1 (w, i + 1, 1) = dcs->c[i];

  gsl_cheb_free (cs);
  gsl_cheb_free (dcs);
  ent_Clean (e1);

  //
  // write result to output list
  //
  Btree *bw = btree_Create ();
  // coefficients of the interpolation
  install (bw, RLAB_NAME_DIFFINT_TCH_COEF, ent_Assign_Rlab_MDR(w));
  //  interval
  install (bw, RLAB_NAME_DIFFINT_TCH_INT, ent_Assign_Rlab_MDR( mdr_Copy(xrange) ));
  // degree of the fit
  install (bw, RLAB_NAME_DIFFINT_TCH_DEG, ent_Create_Rlab_Double(n));
  // label
  install (bw, RLAB_NAME_DIFFINT_TCH_NAME, ent_Create_Rlab_String( "tchebyshev" ));

  return ent_Assign_Rlab_BTREE (bw);
}

Ent *
ent_tchebyshev_int (int nargs, Datum args[])
{
  Ent *e1=0;
  int n, i, cc;
  MDR *c=0, *w, *xrange=0;
  gsl_cheb_series *cs, *dcs;

  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) != BTREE)
    rerror ("tchebyval: first argument has to be list <<coef;interval;max_order;name>>!");

  // <<coef;interval;max_order;name>>
  ListNode *node;
  node = btree_FindNode (ent_data (e1), RLAB_NAME_DIFFINT_TCH_COEF);
  if (node != 0)
    c = class_matrix_real (var_ent (node));
  else
    rerror ("tchebyval: entry '" RLAB_NAME_DIFFINT_TCH_COEF "' has to be vector of non-zero size!");

  node = btree_FindNode (ent_data (e1), RLAB_NAME_DIFFINT_TCH_INT);
  if (node != 0)
    xrange = class_matrix_real (var_ent (node));
  else
    rerror ("tchebyval: entry '" RLAB_NAME_DIFFINT_TCH_INT "' has to be vector [xlo,xhi]!");

  cc = MNC (c);
  if (cc > 1)
    rerror ("tchebyder: 'c' is a column-vector!");

  n = MNR (c) - 1;

  cs = gsl_cheb_alloc (n);
  cs->order = n;
  cs->a = MdrV0(xrange,0);
  cs->b = MdrV0(xrange,1);
  dcs = gsl_cheb_alloc (n);
  dcs->order = n;
  dcs->a = MdrV0(xrange,0);
  dcs->b = MdrV0(xrange,1);

  for (i = 0; i <= n; i++)
    cs->c[i] = Mdr1 (c, i + 1, 1);

  w = mdr_Create (n + 1, 1);

  gsl_cheb_calc_integ (dcs, cs);

  for (i = 0; i <= n; i++)
    Mdr1 (w, i + 1, 1) = dcs->c[i];

  // cleanup
  gsl_cheb_free (cs);
  gsl_cheb_free (dcs);
  ent_Clean (e1);

  //
  // write result to output list
  //
  Btree *bw = btree_Create ();
  // coefficients of the interpolation
  install (bw, RLAB_NAME_DIFFINT_TCH_COEF, ent_Assign_Rlab_MDR(w));
  //  interval
  install (bw, RLAB_NAME_DIFFINT_TCH_INT, ent_Assign_Rlab_MDR( mdr_Copy(xrange) ));
  // degree of the fit
  install (bw, RLAB_NAME_DIFFINT_TCH_DEG, ent_Create_Rlab_Double(n));
  // label
  install (bw, RLAB_NAME_DIFFINT_TCH_NAME, ent_Create_Rlab_String( "tchebyshev" ));

  return ent_Assign_Rlab_BTREE (bw);
}

//
// The interface to the user-specified function.
//  MDR params is passed by default,
//  no need for it to be the argument of the function
//
static double
diffint_gslrlab_func (double x, void *dummy)
{
  Ent *rent = 0;
  double rval;
  // set x
  MDPTR(xmdr) = (void *) &x;

  if (pent)
    rent = ent_call_rlab_script_2args(fname, xent, pent);
  else
    rent = ent_call_rlab_script_1arg (fname, xent);

  if (ent_type(rent)!=MATRIX_DENSE_REAL)
    rerror (THIS_SOLVER ": " RLAB_ERROR_RHS_FUNC_MUST_RETURN_MDR);

  rval = class_double (rent);

  ent_Clean (rent);
  return rval;
}
