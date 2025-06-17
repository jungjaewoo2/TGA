// Copyright (C) 2003-2005 Marijan Kostrun
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
#include <gsl/gsl_mode.h>
#include <gsl/gsl_precision.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_machine.h>
#include <gsl/gsl_integration.h>

// standard headers
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>

// naming convention for the solver parameters
#include "rfileio.h"
#include "rlab_solver_parameters_names.h"

static double
lt_gslrlab_f (double x, void *dummy);

static double
conv_gslrlab_f (double x, void *dummy);

static MDR *xmdr;
static Ent *xent=0, *pent=0, *pent2=0;
static Ent *fname, *fname2;

// laplace transform:
static double cs;
static int    ks;
// convolution:
static double XOFF;


// ******************************************************************
// ******************************************************************
// *                                                                *
// * Convolution of two functions on [a,b] using the gsl nintegrate *
// *                                                                *
// ******************************************************************
// ******************************************************************
#undef  THIS_SOLVER
#define THIS_SOLVER "nintconv"
Ent *
ent_convolution (int nargs, Datum args[])
{
  Ent *e2=0, *e4=0, *e5=0, *e6=0, *eo=0;
  MDR *w=0, *x=0, *xo=0;

  double eabs = 1e-6, erel = 0.01, r, *dummy, re, ddummy, x1, dx;
  double dzero = 0.0;

  int status=0, max_iter=1000, imethod=3, i, j, idummy;
  int ires=0;

  ListNode *node;
  FILE *rlab_stderr = (!RLAB_STDERR_DS) ? stderr : RLAB_STDERR_DS;


  gsl_function F;
  F.function = &conv_gslrlab_f;
  F.params = &dummy;

  // initialize the pointers
  xent=0;
  pent=0;
  pent2=0;

  if (nargs < 5)
    rerror (THIS_SOLVER ": requires at least five arguments!");

  // func1 pointer
  fname = bltin_get_ent(args[0]);
  if (!isfuncent(fname))
    rerror (THIS_SOLVER ": " RLAB_ERROR_ARG1_FUNC_VAR "\n");

  // its parameter entity
  e2 = bltin_get_ent (args[1]);
  pent = 0;
  if (ent_type (e2) != UNDEF)
  {
    pent = ent_Copy (e2);
//    mdr_Print(ent_data(e2), rlab_stderr);
  }

  // func1 pointer
  fname2 = bltin_get_ent(args[2]);
  if (!isfuncent(fname2))
    rerror (THIS_SOLVER ": " RLAB_ERROR_ARG3_FUNC_VAR "\n");

  // parameter entity
  e4 = bltin_get_ent (args[3]);
  pent2 = 0;
  if (ent_type (e4) != UNDEF)
  {
    pent2 = ent_Copy (e4);
//    mdr_Print(ent_data(e4), rlab_stderr);
  }

  // X = [xlo,xhi] or X=xhi, with xlo=0 implied
  e5 = bltin_get_ent (args[4]);
  if (ent_type (e5) != MATRIX_DENSE_REAL)
    rerror (THIS_SOLVER ": incorrect integration range (fourth argument)");

  x = class_matrix_real (e5);
  if (!EQVECT(x) && MNC(x) != 2)
    rerror (THIS_SOLVER ": incorrect integration range (fourth argument)");

  // is XOFF given
  if (nargs > 5)
  {
    e6 = bltin_get_ent (args[5]);
    if (ent_type (e6) == MATRIX_DENSE_REAL)
    {
      xo = class_matrix_real (e6);
      if (!EQVECT(xo))
        rerror (THIS_SOLVER ": incorrect offset (fifth argument)");
    }
  }

  // last argument (6 or 7) gives the options for the solver
  if (nargs > 5)
  {
    eo = bltin_get_ent (args[nargs-1]);
    if (ent_type (eo) == BTREE)
    {
      // eabs
      node = btree_FindNode (ent_data (eo), RLAB_NAME_GSL_CONV_EABS);
      if (node != 0)
      {
        ddummy = class_double (var_ent (node));
        if (ddummy > 0.0)
          eabs = ddummy;
      }
      // erel
      node = btree_FindNode (ent_data (eo), RLAB_NAME_GSL_CONV_EREL);
      if (node != 0)
      {
        ddummy = class_double ( var_ent (node) );
        if (ddummy >= 0.0)
          erel = ddummy;
      }
      // method
      node = btree_FindNode (ent_data (eo), RLAB_NAME_GSL_CONV_IKEY);
      if (node != 0)
      {
        idummy = (int) class_double (var_ent (node));
        if (idummy == 0 || idummy == 1 || idummy == 2)
          imethod = idummy;
      }
      // max iterations
      node = btree_FindNode (ent_data (eo), RLAB_NAME_GSL_CONV_MAXI);
      if (node != 0)
      {
        idummy = (int) class_double (var_ent (node));
        if (idummy > 0)
          max_iter = idummy;
      }
      // resolution
      node = btree_FindNode (ent_data (eo), RLAB_NAME_GSL_CONV_RES);
      if (node != 0)
      {
        idummy = (int) class_double (var_ent (node));
        if (idummy > 0)
          ires = idummy;
      }
    }
  }

  //
  // Set up ENTITIES for user-function.
  // x
  xmdr = mdr_Create(1,1);
  xent = ent_Assign_Rlab_MDR(xmdr);
  ent_IncRef (xent);

  gsl_integration_workspace *wks = gsl_integration_workspace_alloc (max_iter);
  gsl_set_error_handler_off ();

  w = mdr_Create(MNR(x), 1);
  if (MNC(x) == 1)
  {
    // xint = [a1; a2; ...] with c=[c1,c2...]
    // is intepreted as an integration on
    //  [0, a(1)] with c1
    //  [0, a(2)] with c2
    //  ..
    //  [0, a(k)] with c(k)
    for (i=0; i < MNR(x); i++)
    {
      MdrV0(w,i) = r = 0;
      XOFF = MdrV0(x,i);
      if (xo)
      {
        j = i > (xo->nrow*xo->ncol - 1) ? (xo->nrow * xo->ncol - 1): i;
        XOFF = MdrV0(xo, j);
      }
      if (ires)
      {
        // divide integration interval to ires steps
        if (i)
        {
          x1 = MdrV0(x,i-1);
          dx = (MdrV0(x,i)-MdrV0(x,i-1));
        }
        else
        {
          x1 = dzero;
          dx = MdrV0(x,i);
        }
        dx /= ires;
        for (j=0; j<ires; j++)
        {
          status =
              gsl_integration_qag (&F, x1, x1+dx, eabs, erel, max_iter, imethod, wks,
                                    &r, &re);
          x1 += dx;
          MdrV0(w,i) = MdrV0(w,i) + r;
        }
      }
      else
      {
        // integrate in a single step. cross your fingers!
        status =
            gsl_integration_qag (&F, dzero, MdrV0(x,i), eabs, erel, max_iter, imethod, wks,
                                 &r, &re);
        MdrV0(w,i) = r;
      }
    }
  }
  else
  {
    // xint = [a1,b1; a2,b2; ...] with c=[c1,c2...]
    // is intepreted as an integration on
    //  [a(1), b(1)] with c1
    //  [a(2), b(2)] with c2
    //  ..
    //  [a(k), b(k)] with c(k)
    for (i=0; i < MNR(x); i++)
    {
      MdrV0(w,i) = r = 0;
      XOFF = MdrV0(x,i);
      if (xo)
      {
        j = i > (xo->nrow * xo->ncol - 1) ? (xo->nrow * xo->ncol - 1): i;
        XOFF = MdrV0(xo,j);
      }

      if (ires)
      {
        // divide integration interval to ires steps
        x1 = Mdr0(x,i,0);
        dx = (Mdr0(x,i,1)-Mdr0(x,i,0))/ires;
        for (j=0; j<ires; j++)
        {
          status =
              gsl_integration_qag (&F, x1, x1+dx, eabs, erel, max_iter, imethod, wks,
                                    &r, &re);
          x1 += dx;
          MdrV0(w,i) += r;
        }
      }
      else
      {
        status =
            gsl_integration_qag (&F, Mdr0(x,i,0), Mdr0(x,i,1), eabs, erel,
                                  max_iter, imethod, wks, &r, &re);
        MdrV0(w,i) = r;
      }
    }
  }

  if (status == GSL_EDIVERGE)
    fprintf (rlab_stderr, THIS_SOLVER ": Integral divergent or slowly convergent!\n");

  gsl_integration_workspace_free (wks);

  // destroy entity and the MDR associated with it
  ent_DecRef (xent);
  ent_Destroy (xent);

  ent_Clean (pent);
  ent_Clean (pent2);

  ent_Clean (fname);
  ent_Clean (e2);
  ent_Clean (fname2);
  ent_Clean (e4);
  ent_Clean (e5);
  ent_Clean (e6);
  ent_Clean (eo);

  return ent_Assign_Rlab_MDR(w);
}


// **************************************************************
// * Laplace transform using the gsl nintegrate                        *
// **************************************************************
#undef  THIS_SOLVER
#define THIS_SOLVER "nintlt"
Ent *
ent_laplacetransform (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *e3=0, *eo=0;
  MDR *w, *x=0, *s=0, *k=0;
  double eabs = 1e-6, erel = 0.01, r, *dummy, re, result, xlo = 0, xhi = 0, ddummy;
  int cstat = 0, status=0, max_iter = 1000, imethod = 3, i, xcol, is, idummy;
  int ik, J, kill_x=0;

  ListNode *node;
  FILE *rlab_stderr = (!RLAB_STDERR_DS) ? stderr : RLAB_STDERR_DS;

  pent=0;

  gsl_function F;
  F.function = &lt_gslrlab_f;
  F.params = &dummy;

  if (nargs < 2)
  {
    fprintf (rlab_stderr,
             THIS_SOLVER ": Laplace transform of a real function of a real variable.\n");
    fprintf (rlab_stderr, THIS_SOLVER ": Format:\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ":   y = nintlt(s,k,func /,p,[x1,x2..],opts/),\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": where  y=y(s) is the transform of a function\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": func(x /,p/) with respect to variable x. The row-vector\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": [x1..] specifies integration intervals. If func is\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": defined for all x>0 the last element should be inf().\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": The conjugate variable s can be a column vector.\n");
    rerror (THIS_SOLVER "requires at least three arguments!");
  }

  // s
  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) != MATRIX_DENSE_REAL)
    rerror (THIS_SOLVER ": " RLAB_ERROR_ARG1_MDR_VECTOR);
  s = class_matrix_real (e1);
  if (!EQVECT(s))
    rerror (THIS_SOLVER ": " RLAB_ERROR_ARG1_MDR_VECTOR);

  // k
  J = 1;
  e2 = bltin_get_ent (args[J]);
  if (ent_type (e2) == MATRIX_DENSE_REAL)
  {
    k = class_matrix_real (e2);
    ++J;
  }

  // func pointer at 2 or 3
  fname = bltin_get_ent(args[J++]);
  if (!isfuncent(fname))
    rerror (THIS_SOLVER ": " RLAB_ERROR_ARG2_FUNC_VAR "\n");

  //
  // parameter entity
  //
  if (nargs > J)
  {
    e2 = bltin_get_ent (args[J++]);
    pent = 0;
    if (ent_type (e2) != UNDEF)
      pent = ent_Copy (e2);
  }

  // X = [x1,x2...]
  if (nargs > J)
  {
    e3 = bltin_get_ent (args[J++]);
    if (ent_type (e3) == MATRIX_DENSE_REAL)
      x = class_matrix_real (e3);
  }

  if (!x)
  {
    // default integration interval [0,inf()]
    x = mdr_Create (1, 3);
    MdrV0 (x, 0) = 0;
    MdrV0 (x, 1) = 1;
    MdrV0 (x, 2) = create_inf ();
    kill_x = 1;
  }

  xcol = SIZE(x);
  if (xcol < 2)
    rerror(THIS_SOLVER ": incorrect row-vector of integration intervals !");

  //
  // options for the solver
  //
  if (nargs > J)
  {
    eo = bltin_get_ent (args[J++]);
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
  }

  // Set up ENTITIES for user-function.
  // x
  xmdr = mdr_Create (1,1);
  xent = ent_Assign_Rlab_MDR (xmdr);
  ent_IncRef (xent);

  gsl_integration_workspace *wks = gsl_integration_workspace_alloc (max_iter);
  gsl_set_error_handler_off ();

  if (k)
  {
    w = mdr_Create (SIZE(s), SIZE(k));
    for (is = 1; is <= MNR (s) * MNC (s); is++)
      for (ik = 1; ik <= MNR (k) * MNC (k); ik++)
    {
    // this is passed to func via static storage
      ks = MdrV1 (k, ik);
      cs = MdrV1 (s, is);
      //
      // Decide on what integration we do
      //
      for (i = 2; i <= xcol; i++)
      {
      // verify integration intervals
        if (MdrV1 (x, i) == -create_inf ())
          printf (THIS_SOLVER ": -inf() cannot be an integration limit!");
        if (MdrV1 (x, i - 1) == create_inf ())
          printf
              (THIS_SOLVER ": +inf() cannot be an intermediate integration point!");
      }
      result = 0;
      for (i = 2; i <= xcol; i++)
      {
        // is the second point inf
        if (MdrV1 (x, i) == create_inf ())
        {
          // use _qagiu for [x[i-1],inf]
          xlo = MdrV1 (x, i - 1);
          status =
              gsl_integration_qagiu (&F, xlo, eabs, erel, max_iter, wks, &r, &re);
          cstat += status;
          result = result + r;
          continue;
        }
        // second point must be normal. use qag on it
        xlo = MdrV1 (x, i - 1);
        xhi = MdrV1 (x, i);
        status =
            gsl_integration_qag (&F, xlo, xhi, eabs, erel, max_iter, imethod, wks,
                                  &r, &re);
        if (status != GSL_SUCCESS)
        {
          // if it doesn't work try qags (singularity in the interval
          status =
              gsl_integration_qags (&F, xlo, xhi, eabs, erel, max_iter, wks, &r,
                                     &re);
        }
        result = result + r;
      }
      if (status == GSL_SUCCESS || status == GSL_EDIVERGE)
        Mdr1 (w, is,ik) = result;
      else
        Mdr1 (w, is, ik) = create_nan();

      if (status == GSL_EDIVERGE)
        fprintf (rlab_stderr,
                 THIS_SOLVER ": Integral divergent or slowly convergent!\n");
    }
  }
  else
  {
    w = mdr_Create (MNR(s) * MNC(s), 1);
    ks = 0.0;
    for (is = 1; is <= MNR (s) * MNC (s); is++)
    {
      cs = MdrV1 (s, is);
      //
      // Decide on what integration we do
      //
      for (i = 2; i <= xcol; i++)
      {
        // verify integration intervals
        if (MdrV1 (x, i) == -create_inf ())
          printf (THIS_SOLVER ": -inf() cannot be an integration limit!");
        if (MdrV1 (x, i - 1) == create_inf ())
          printf
              (THIS_SOLVER ": +inf() cannot be an intermediate integration point!");
      }
      result = 0;
      for (i = 2; i <= xcol; i++)
      {
        r = 0;
        xlo = MdrV1 (x, i - 1);
        // is the second point inf
        if (MdrV1 (x, i) == create_inf ())
        {
          // use _qagiu for [x[i-1],inf]
          status =
              gsl_integration_qagiu (&F, xlo, eabs, erel, max_iter, wks, &r, &re);
        }
        else
        {
          // second point must be normal. use qag on it
          xhi = MdrV1 (x, i);
          status =
              gsl_integration_qag (&F, xlo, xhi, eabs, erel, max_iter, imethod, wks,
                                   &r, &re);
          if (status != GSL_SUCCESS && status != GSL_EDIVERGE)
          {
            // if it doesn't work try qags (singularity in the interval
            status =
                gsl_integration_qags (&F, xlo, xhi, eabs, erel, max_iter, wks, &r,
                                      &re);
          }
        }
        result = result + r;
      }

      if (status == GSL_SUCCESS || status == GSL_EDIVERGE)
        MdrV1 (w, is) = result;
      else
        MdrV1 (w, is) = create_nan();

      if (status == GSL_EDIVERGE)
        fprintf (rlab_stderr,
                 THIS_SOLVER ": Integral divergent or slowly convergent!\n");
    }
  }

  gsl_integration_workspace_free (wks);

  ent_DecRef (xent);
  ent_Destroy (xent);

  if (kill_x)
    mdr_Destroy(x);

  ent_Clean (pent);

  ent_Clean (e1);
  ent_Clean (fname);
  ent_Clean (e2);
  ent_Clean (e3);
  ent_Clean (eo);

  return ent_Assign_Rlab_MDR(w);
}

//
// The interface to the user-specified function.
//  MDR params is passed by default, no need for it to be the argument of the function
//
static double
lt_gslrlab_f (double x, void *dummy)
{
  Ent *rent = 0;
  double retv;

  MdrV0(xmdr,0) = x;
  if (pent)
    rent = ent_call_rlab_script_2args (fname, xent, pent);
  else
    rent = ent_call_rlab_script_1arg  (fname, xent);

  retv = class_double(rent);
  ent_Clean (rent);

  if (cs)
  {
    retv *= exp (-cs * x);
    if (ks)
      retv *= pow(cs,ks);
  }
  else
  {
    if (ks)
      retv = 0;
  }
  return retv;
}


static double
conv_gslrlab_f (double x, void *dummy)
{
  Ent *rent = 0;
  double retv;

  //
  // call func no. 1: f1(x)
  //
  MdrV0(xmdr,0) = x;
  if (pent)
    rent = ent_call_rlab_script_2args (fname, xent, pent);
  else
    rent = ent_call_rlab_script_1arg  (fname, xent);

  retv = class_double(rent);
  ent_Clean (rent);

  //
  // call func no. 2: f2(XOFF-x)
  //
  MdrV0(xmdr,0) = XOFF - x;
  if (pent2)
    rent = ent_call_rlab_script_2args (fname2, xent, pent2);
  else
    rent = ent_call_rlab_script_1arg  (fname2, xent);

  retv *= class_double (rent);
  ent_Clean (rent);
  return retv;
}
