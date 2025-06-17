//
// pz1rl.c
// rlab <-> PZEROS.F interface
//

/*  This file is a part of rlabplus
   Copyright (C) 2005 M. Kostrun

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
#include "ent.h"
#include "btree.h"
#include "symbol.h"
#include "mdr.h"
#include "mdrf1.h"
#include "mdrf2.h"
#include "mdc.h"
#include "mdcf1.h"
#include "mdcf2.h"
#include "msr.h"
#include "msrf1.h"
#include "msrf2.h"
#include "msc.h"
#include "mscf1.h"
#include "mscf2.h"
#include "util.h"
#include "mathl.h"
#include "class.h"
#include "mem.h"
#include "mds.h"
#include "list.h"
#include "bltin.h"
#include "function.h"
#include "lp.h"

#include "rfileio.h"
#include "pzeros.h"

// naming convention for the solver parameters
#include "rlab_solver_parameters_names.h"

static MDC *x;
static Ent *ent_x, *ent_p;
static Ent *series_fname;


//
// basic solver
//
#undef THIS_SOLVER
#define THIS_SOLVER "polyroots"
Ent *
ent_pzeros (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0;
  MDR *x, *radius, *apoly=0, *apolyr=0, *ierr=0;
  MDC *cx=0, *poly=0, *root=0;

  ListNode *node;
  Btree *bw = btree_Create ();
  FILE *rlab_stderr = (!RLAB_STDERR_DS) ? stderr : RLAB_STDERR_DS;

  int idummy=0, ldeg;
  double dummy=0;

  //
  // Initialize default parameters
  //
  int pz_maxit   = 100;                    // max number of iterations
  double pz_eps  = 1.1102230246252e-16;    // machine precision
  double pz_dmin = 4.9406564584126e-324;   // smallest positive real
  double pz_dmax = 1.7976931348622e+308;   // largest positive real
  int iter=0, nzer=0, ndeg=0, i, pzer=0;

  //
  // Load and Check arguments.
  //
  if (nargs != 1 && nargs !=2)
  {
    fprintf(rlab_stderr, THIS_SOLVER ": Finds all roots of a polynomial.\n");
    fprintf(rlab_stderr, THIS_SOLVER ": Format:\n");
    fprintf(rlab_stderr, THIS_SOLVER ":   y = polyroots(P/,options/),\n");
    fprintf(rlab_stderr, THIS_SOLVER ": where 'P' is a row-vector with the coefficients of the\n");
    fprintf(rlab_stderr, THIS_SOLVER ": polynomial in standard format. The parameter 'options'\n");
    fprintf(rlab_stderr, THIS_SOLVER ": is a list with entries: eps, dmin, dmax, maxit.\n");
    rerror (THIS_SOLVER ": " RLAB_ERROR_ONE_OR_TWO_ARG_REQUIRED);
  }

  //
  // Did user specify the options? If so, get them first.
  //
  if (nargs == 2)
  {
    e2 = bltin_get_ent (args[1]);
    if (ent_type (e2) == BTREE)
    {
      // eps
      node = btree_FindNode (ent_data (e2), RLAB_NAME_POLY_EPS);
      if (node != 0)
      {
        dummy = class_double (var_ent (node));
        if (dummy >= pz_eps)
          pz_eps = dummy;
      }
      // dmin
      node = btree_FindNode (ent_data (e2), RLAB_NAME_POLY_DMIN);
      if (node != 0)
      {
        dummy = class_double ( var_ent (node) );
        if (dummy >= pz_dmin)
          pz_dmin = dummy;
      }
      // dmax
      node = btree_FindNode (ent_data (e2), RLAB_NAME_POLY_DMAX);
      if (node != 0)
      {
        dummy = class_double (var_ent (node));
        if (dummy <= pz_dmax && dummy>0)
          pz_dmax = dummy;
      }
      // max iterations
      node = btree_FindNode (ent_data (e2), RLAB_NAME_POLY_MAXITER);
      if (node != 0)
      {
        idummy = (int) class_double (var_ent (node));
        if (idummy > 0)
          pz_maxit = idummy;
      }
    }
  }

  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) == MATRIX_DENSE_REAL)
  {
    x = ent_data (e1);
    ndeg = MNR (x) * MNC (x) - 1;
    // figure out trailing zeros in polynomial [a_(n-1),..,a_(k),0,0..0]
    nzer = 0;
    for (i = ndeg ; i >= 0; i--)
    {
      if (MdrV0(x,i) < pz_dmin && -pz_dmin < MdrV0(x,i))
        nzer++;
      else
        break;
    }
    // figure out the leading zeros in polynomial [0,..0,a_(k),..,a_(0)]
    pzer = 0;
    for (i = 0; i <= ndeg; i++)
    {
      if (MdrV0(x,i) < pz_dmin && -pz_dmin < MdrV0(x,i))
        pzer++;
      else
        break;
    }
    poly = mdc_Create(1, ndeg+1-nzer-pzer);
    if(x->type == RLAB_TYPE_INT32)
    {
      for (i=1+pzer; i<=ndeg+1-nzer; i++)
      {
        Mdc1r(poly,1,ndeg+2-nzer-i) = MdiV1(x,i);
        Mdc1i(poly,1,ndeg+2-nzer-i) = 0;
      }
    }
    else
    {
      for (i=1+pzer; i<=ndeg+1-nzer; i++)
      {
        Mdc1r(poly,1,ndeg+2-nzer-i) = MdrV1(x,i);
        Mdc1i(poly,1,ndeg+2-nzer-i) = 0;
      }
    }
  }
  else if (ent_type (e1) == MATRIX_DENSE_COMPLEX)
  {
    cx = ent_data (e1);
    ndeg = SIZE (cx) - 1;
    nzer = 0;
    for (i = ndeg ; i >= 0; i--)
    {
      if ( cabs(MdcV0(cx,i)) < pz_dmin )
        nzer++;
      else
        break;
    }
    pzer = 0;
    for (i = 0; i <= ndeg; i++)
    {
      if ( cabs(MdcV0(cx,i)) < pz_dmin )
        pzer++;
      else
        break;
    }

    poly = mdc_Create(1, ndeg+1-nzer-pzer);
    for (i=1+pzer; i<=ndeg+1-nzer; i++)
    {
      Mdc1r(poly,1,ndeg+2-nzer-i) = MdcV1r(cx,i);
      Mdc1i(poly,1,ndeg+2-nzer-i) = MdcV1i(cx,i);
    }
  }
  else
  {
    rerror (THIS_SOLVER ": Unknown type of the first input argument");
  }

  //
  // allocate the variables
  //
  ldeg = ndeg - nzer - pzer;
  if (ldeg > 0)
  {
    // adjust the degree of polynomial
    ndeg = ndeg - pzer;

    root   = mdc_Create(ndeg,1);
    radius = mdr_Create(ndeg,1);
    ierr   = mdi_Create(ndeg,1);
    apoly  = mdr_Create(1, ndeg+1);
    apolyr = mdr_Create(1, ndeg+1);
    for (i=0; i<nzer; i++)
    {
      MdcV0r (root,i) = 0;
      MdcV0i (root,i) = 0;
      MdrV0  (radius,i) = pz_dmin;
      MdiV0  (ierr, i)  = 0;
    }

    if (ldeg > 1)
    {
      // poly is at least of degree two
      POLZEROS (&ldeg, MDCPTR(poly), &pz_eps, &pz_dmax, &pz_dmin,
                &pz_maxit, &MdcV0(root,nzer), &MdrV0(radius,nzer), &MdiV0(ierr,nzer), &iter,
                &MdrV0(apoly,nzer), &MdrV0(apolyr,nzer));
    }
    else
    {
      // poly is of degree 1:
      //  [a_{k},a_{k-1}]
      MdcV0(root,nzer) =  MdcV0(poly,1) / MdcV0(poly,0);
      MdcV0(root,nzer) = -MdcV0(root,nzer);
      MdrV0(radius,nzer) = pz_dmin;
      MdiV0(ierr,nzer)   = 0;
    }
  }
  else
  {
    // zeroth degree polynomial: no roots whatsoever
    root   = mdc_Create(0,0);
    radius = mdr_Create(0,0);
    ierr   = mdi_Create(0,0);
  }

  if (poly)
    mdc_Destroy (poly);
  if (apoly)
    mdr_Destroy (apoly);
  if (apolyr)
    mdr_Destroy (apolyr);

  ent_Clean (e1);
  ent_Clean (e2);

  //
  // return the result as a list
  //
  // roots of the polynomial
  install  (bw, "roots", ent_Assign_Rlab_MDC(root));
  // radius
  install  (bw, "radius", ent_Assign_Rlab_MDR(radius));
  // err: was calculation succesful
  install  (bw, "status", ent_Assign_Rlab_MDR(ierr));

  return ent_Assign_Rlab_BTREE(bw);
}

#undef  THIS_SOLVER
#define THIS_SOLVER "nseries"
//
// 579.f: taylor/laurent series of a function at a point
//
Ent *
ent_nseries (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *e3=0, *e4=0;
  Complex z, c;
  int n, ic, i;
  double ddummy, r = 0.1, eps = 1e-14;

  FILE *rlab_stderr = (!RLAB_STDERR_DS) ? stderr : RLAB_STDERR_DS;

  ListNode *node;

  MDC *rs;
  MDR *er;

  ent_p = 0;

  if (nargs < 4 || nargs > 5)
  {
    fprintf (rlab_stderr,
             THIS_SOLVER ": Finds series expansion of a function at a point.\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": Format:\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ":   a = nseries(f,/p/,x0,n /,options/),\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": where\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ":   f = function(x /,p/), is the function being expanded,\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ":   x0 is a point at which the expansion is being sought,\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ":   n is the order of the power series of the expansion, while\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ":   options=<<eps;r>> is a list containing error and radius of convergence of\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": the resulting power series.\n");
    rerror ("requires at three to five arguments !");
  }

  //
  // Get function ptr
  //
  series_fname = bltin_get_ent(args[0]);
  if (!isfuncent(series_fname))
    rerror (THIS_SOLVER ": " RLAB_ERROR_ARG1_FUNC_VAR "\n");

  // parameter array for the function
  e1  = bltin_get_ent (args[1]);
  if (ent_type(e1) != UNDEF)
    ent_p = ent_Copy (e1);

  // x
  e2 = bltin_get_ent (args[2]);
  if (ent_type (e2) == MATRIX_DENSE_REAL)
  {
    z = class_double (e2);
  }
  else if (ent_type (e2) == MATRIX_DENSE_COMPLEX)
  {
    MDC * xc = ent_data(e2);
    z = MdcV0(xc,0);
  }
  else
    rerror (THIS_SOLVER ": Second argument 'x0' must be real or complex scalar !");

  // n
  e3 = bltin_get_ent (args[3]);
  if (ent_type (e3) != MATRIX_DENSE_REAL)
    rerror (THIS_SOLVER ": Third argument 'n' must be positive integer scalar !");
  n = (int) class_double (e3) + 1;
  if (n<1)
    rerror (THIS_SOLVER ": Third argument 'n' must be positive integer scalar !");

  // radius
  if (nargs > 4)
  {
    e4 = bltin_get_ent (args[4]);
    if (ent_type (e4) == BTREE)
    {
      // eps
      node = btree_FindNode (ent_data (e4), RLAB_NAME_NSER_EPS);
      if (node != 0)
      {
        ddummy = class_double (var_ent (node));
        if (ddummy > 0)
          eps = ddummy;
      }
      // r
      node = btree_FindNode (ent_data (e4), RLAB_NAME_NSER_R);
      if (node != 0)
      {
        ddummy = class_double ( var_ent (node) );
        if (ddummy > 0)
          r = ddummy;
      }
    }
  }

  //
  // Set up ENTITIES for user-function.
  // x
  x = mdc_CreateEmpty (1,1);
  ent_x = ent_Assign_Rlab_MDC (x);
  ent_IncRef (ent_x);

  ic = 0; // do series expansion

  rs = mdc_Create(1,n);
  er = mdr_Create(1,n);

  CPSC(&z, &n, &ic, &r, MDCPTR(rs), MDRPTR(er), &eps);

  // flip the result l->r
  for (i = 1; i <= n/2; i++)
  {
    c = MdcV1(rs,i);
    MdcV1(rs,i) = MdcV1(rs,n-i+1);
    MdcV1(rs,n-i+1) = c;
  }

  // clean-up
  MDPTR(x) = 0;
  ent_DecRef (ent_x);
  ent_Destroy (ent_x);

  ent_Clean (ent_p); ent_p = 0;

  ent_Clean (e1);
  ent_Clean (e2);
  ent_Clean (e3);
  ent_Clean (e4);

  //
  // write result to output list
  //
  Btree *bw = btree_Create ();
  install  (bw, "coef", ent_Assign_Rlab_MDC(rs));
  install  (bw, "cov", ent_Assign_Rlab_MDR(er));
  return ent_Assign_Rlab_BTREE(bw);
}


int
CPS1FUNC (Complex * z, Complex *zval)
{
  Ent *rent = 0;

  // x
  MDPTR(x) = (void *) z;

  if (ent_p)
    rent = ent_call_rlab_script_2args(series_fname, ent_x, ent_p);
  else
    rent = ent_call_rlab_script_1arg (series_fname, ent_x);

  if (ent_type (rent) != MATRIX_DENSE_REAL && ent_type (rent) != MATRIX_DENSE_COMPLEX)
    rerror (THIS_SOLVER ": " RLAB_ERROR_RHS_FUNC_MUST_RETURN_MDR_OR_MDC);

  *zval = class_complex(rent);
  ent_Clean (rent);
  return 1;
}
