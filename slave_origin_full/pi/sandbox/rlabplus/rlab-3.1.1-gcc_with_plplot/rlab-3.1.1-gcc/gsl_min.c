// Copyright (C) 2003-2006 Marijan Kostrun
//   part of rlabplus project, see http://rlabplus.sourceforge.net
//
// GSL Science Library - minimization of scalar function in any dimension
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

//
// conmax.f
//
#include "conmax.h"

//
// gsl headers
//
#include <gsl/gsl_mode.h>
#include <gsl/gsl_precision.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_machine.h>
#include <gsl/gsl_min.h>
#include <gsl/gsl_multimin.h>

#include "gsl_rlab.h"

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
// the GSL
//
// scalar function of a scalar variable
static double min_gslrlab_f (double x, void *dummy);

// arguments for rlab-scripted functions
static MDR *xmdr;
static Ent *xent;
static Ent *pent;
static Ent *min_f_name;

// **************************************************************
// * Builtin interface to gsl 1d minimization routines          *
// **************************************************************
#undef  THIS_SOLVER
#define THIS_SOLVER "findmin"
Ent *
ent_minimize_1d_bracket (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *eo=0;
  MDR *x;

  double eabs = 0, erel = 0.01, xlo, xhi, r, *dummy, ddummy;
  int status, iter = 0, max_iter = 1000, imethod = 0, idummy;
  ListNode *node;

  gsl_function F;
    F.function = &min_gslrlab_f;
    F.params = &dummy;
  const gsl_min_fminimizer_type *T;
  gsl_min_fminimizer *s=0;

  if (nargs < 3)
    rerror (THIS_SOLVER ": requires at least 3 arguments");

  // Get function ptr
  min_f_name = bltin_get_ent(args[0]);
  if (!isfuncent(min_f_name))
    rerror (THIS_SOLVER ": " RLAB_ERROR_ARG1_FUNC_VAR "\n");

  //
  // parameter entity
  //
  if (nargs > 1)
  {
    e1 = bltin_get_ent (args[1]);
    pent = 0;
    if (ent_type(e1) != UNDEF)
      pent = ent_Copy( e1 );
  }

  //
  // x = [xlo,xhi]
  //
  e2 = bltin_get_ent (args[2]);
  if (ent_type (e2) != MATRIX_DENSE_REAL)
    rerror (THIS_SOLVER ": improper x=[xlo,xhi] !");
  x = class_matrix_real (e2);
  if (!x)
    rerror (THIS_SOLVER ": vector x=[xlo,xhi] must be real!");
  if (MNR(x) * MNC(x) != 2)
    rerror (THIS_SOLVER ": vector x must be [xlo,xhi] !");
  xlo = MdrV0(x,0);
  xhi = MdrV0(x,1);

  if (nargs > 3)
  {
    eo = bltin_get_ent (args[3]);
    if (ent_type (eo) == BTREE)
    {
      // eabs
      node = btree_FindNode (ent_data (eo), RLAB_NAME_GSL_MIN1_EABS);
      if (node != 0)
      {
        ddummy = class_double (var_ent (node));
        if (ddummy > 0.0)
          eabs = ddummy;
      }
      // erel
      node = btree_FindNode (ent_data (eo), RLAB_NAME_GSL_MIN1_EREL);
      if (node != 0)
      {
        ddummy = class_double (var_ent (node));
        if (ddummy > 0.0)
          erel = ddummy;
      }
      // max iterations
      node = btree_FindNode (ent_data (eo), RLAB_NAME_GSL_MIN1_MAXITER);
      if (node != 0)
      {
        idummy = (int) class_double (var_ent (node));
        if (idummy > 0)
          max_iter = idummy;
      }
      // imethod
      node = btree_FindNode (ent_data (eo), RLAB_NAME_GSL_MIN1_IMETHOD);
      if (node != 0)
      {
        idummy = (int) class_double (var_ent (node));
        if (idummy >= 0 && idummy <= 3)
          imethod = idummy;
      }
    }
    if ( ent_type (eo) != UNDEF ) ent_Clean (eo);
  }

  // Set up ENTITIES for user-function.
  // x
  xmdr = mdr_CreateEmpty (1,1);
  xent = ent_Assign_Rlab_MDR (xmdr);
  ent_IncRef (xent);

  gsl_set_error_handler_off ();

  switch (imethod)
  {
  case 0:
    T = gsl_min_fminimizer_brent;
    break;
  case 1:
    T = gsl_min_fminimizer_goldensection;
    break;
  default:
    T = gsl_min_fminimizer_brent;
  }

  s = gsl_min_fminimizer_alloc (T);
  gsl_min_fminimizer_set (s, &F, 0.5 * xlo + 0.5 * xhi, xlo, xhi);
  do
  {
    iter++;
    status = gsl_min_fminimizer_iterate (s);
    r = gsl_min_fminimizer_x_minimum (s);
    xlo = gsl_min_fminimizer_x_lower (s);
    xhi = gsl_min_fminimizer_x_upper (s);
    status = gsl_min_test_interval (xlo, xhi, eabs, erel);
  }
  while (status == GSL_CONTINUE && iter < max_iter);

  ent_Clean (e1);
  ent_Clean (e2);

  ent_Clean (pent);

  MDPTR(xmdr)=0;
  ent_DecRef (xent);
  ent_Destroy (xent);

  ent_Clean  (min_f_name);

  gsl_min_fminimizer_free (s);

  if (status == GSL_SUCCESS)
    return ent_Create_Rlab_Double (r);
  else
    return ent_Assign_Rlab_MDR ( NULL );
}


//
// The interface to the user-specified function.
//  MDR params is passed by default, no need for it to be the argument of the function
//
static double
min_gslrlab_f (double x, void *dummy)
{
  Ent *rent = 0;
  double rval;

  MDPTR (xmdr) = &x;

  if (pent)
    rent = ent_call_rlab_script_2args(min_f_name, xent, pent);
  else
    rent = ent_call_rlab_script_1arg (min_f_name, xent);

  if (ent_type(rent)!=MATRIX_DENSE_REAL)
    rerror (THIS_SOLVER ": " RLAB_ERROR_RHS_FUNC_MUST_RETURN_MDR);

  rval = class_double (rent);

  ent_Clean (rent);
  return rval;
}


