// odesl.c
// This file is a part of RLaB2 Rel.2
//  Copyright (C) 2005-2008 Marijan Kostrun. project rlabplus.

/*
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
   **********************************************************************
*/

#include "rlab.h"
#include "ent.h"
#include "class.h"
#include "symbol.h"
#include "mem.h"
#include "mdr.h"
#include "mdrf1.h"
#include "mds.h"
#include "list.h"
#include "btree.h"
#include "bltin.h"
#include "util.h"
#include "mathl.h"
#include "function.h"
#include "lp.h"
#include "ode.h"

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <time.h>

#include "mathl.h"

// naming convention for the solver parameters
#include "rfileio.h"
#include "rlab_solver_parameters_names.h"
//
// functions: their RLaB names and their Fortran arguments
//
double p_ (double *);
static Ent *sleignp_fname;

double q_ (double *);
static Ent *sleignq_fname;

double w_ (double *);
static Ent *sleignw_fname;

int uv_ (double *, double *, double *, double *, double *, double *, double *);
static Ent *sleignuv_fname;

static MDR *x;
static Ent *ent_param, *ent_x;


//
// sturm-liouville separated boundary value problem: sleign
//
#undef THIS_SOLVER
#define THIS_SOLVER "odesl.eign"
Ent *
ent_odesl_eign (int nargs, Datum args[])
{
  Ent *e4=0, *e5=0, *e6=0, *e7=0, *e8=0, *e9=0, *eo;
  MDR *x5=0, *x6=0, *x8=0, *x9=0, *weig=0, *wfunc=0, *slfun=0, *rslfun=0;

  double a, b, p0ata = 0, qfata = 0, p0atb = 0, qfatb =
    0, a1, a2, b1, b2, tol=-1e-6, mp;
  double lfp, rfp, eig = 0.1, ddummy;
  int intab, numeig, iflag, islfun = 0, nca = 1, ncb = 1, i, j, sf;
  int louts=0, ncol;

  ListNode *node;

  FILE *fptr = NULL;
  FILE *rlab_stderr = (!RLAB_STDERR_DS) ? stderr : RLAB_STDERR_DS;

  time_t timer;

  char *outs=0;

  //
  // Check arguments.
  //
  if (nargs < 6)
  {
    fprintf (rlab_stderr,
             THIS_SOLVER ": Sturm-Liouville Eigenvalue Problem solver (SLEIGN2), for separated\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": boundary conditions,\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ":   -(P(x)y'(x))' + Q(x) y(x) = eig W(x) y(x).\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": Format:\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ":   Y = odesl(P,Q,W,/UV/,/p/,x,bc,eix/,nc,atab/),\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": where, P = f(x/,p/), Q = f(x/,p/), and W = f(x/,p/) are the\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": functions, 'p' is the parameter array, and x = [a,x1, ..,xN, b] is\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": the interval with boundary points 'a' and 'b',\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": bc=[A1,A2;B1,B2] for regular\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": problem (UV not needed), otherwise UV(x)=[u(x),u'(x),v(x),v'(x)]\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": has to be specified by the user. 'eix' is the index of the\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": eigenvalue, while atab=[p0ata,qfata,p0atb,qfatb], nc=[nca,ncb]\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": is the type of endpoints: 1 (regular), 2 (weakly regular),\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": 3 (limit cycle, non-oscilatory), 4 (limit cycle, oscilatory),\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": 5 (limit point, regular at finite point), 6 (limit point),\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": 7 (limit point at infinity), 8 (limit point, bad behaved),\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ":  default is nc=[1,1].\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": Y=<<eigval;eigfun;info>>, the eigenvalue, tabulated eigenfunction\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": on user's mesh and the status info from the solver.\n");
    rerror ("requires at least 7 arguments");
  }

  //
  // options: the last argument
  //
  if (nargs == 11 || nargs == 10 || nargs == 9)
  {
    eo = bltin_get_ent (args[nargs-1]);
    if (ent_type (eo) == BTREE)
    {
      // tol
      node = btree_FindNode (ent_data (eo), RLAB_NAME_ODESL_TOL);
      if (node != 0)
      {
        ddummy = class_double (var_ent (node));
        if (ddummy > 0.0)
          tol = ddummy;
      }
      // standard output
      node = btree_FindNode (ent_data (eo), RLAB_NAME_ODESL_STDOUT);
      if (node != 0)
      {
        if (ent_type(var_ent (node))== MATRIX_DENSE_STRING)
        {
          outs = class_char_pointer(var_ent (node));
          if (outs)
            louts = strlen (outs);
        }
      }
    }
    ent_Clean (eo);
  }

  //
  // P(x)
  //
  sleignp_fname = bltin_get_ent(args[0]);
  if (!isfuncent(sleignp_fname))
    rerror (THIS_SOLVER ": " RLAB_ERROR_ARG1_FUNC_VAR "\n");

  //
  // Q(x)
  //
  sleignq_fname = bltin_get_ent(args[1]);
  if (!isfuncent(sleignq_fname))
    rerror (THIS_SOLVER ": " RLAB_ERROR_ARG2_FUNC_VAR "\n");

  //
  // W(x)
  //
  sleignw_fname = bltin_get_ent(args[2]);
  if (!isfuncent(sleignw_fname))
    rerror (THIS_SOLVER ": " RLAB_ERROR_ARG3_FUNC_VAR "\n");

  //
  // UV(x) if given
  //
  sleignuv_fname = bltin_get_ent(args[3]);
  if (!isfuncent(sleignuv_fname))
  {
    ent_Clean (sleignuv_fname);
    sleignuv_fname = 0;
  }

  //
  // parameters for the functions
  //
  e4 = bltin_get_ent (args[4]);
  ent_param = 0;
  if (ent_type (e4) != UNDEF)
    ent_param = ent_Copy (e4);

  //
  // [a,x1..xN,b], the relative mesh on which the of the eigenfunction is computed
  //
  e5 = bltin_get_ent (args[5]);
  if (ent_type (e5) != MATRIX_DENSE_REAL)
  {
    fprintf (stdout,
             THIS_SOLVER ": 5th argument has to be the interval [a,x1,..,xN,b].\n");
    rerror (THIS_SOLVER ": Improper bound!");
  }
  x5 = class_matrix_real (e5);
  ncol = SIZE (x5);
  if (!EQVECT(x5) || ncol < 2)
    rerror (THIS_SOLVER ": improper x=[a..,b]!");
  a = MdrV1 (x5, 1);
  b = MdrV1 (x5, ncol);
  if (a != -create_inf () && b != create_inf ())
    intab = 1;
  else if (a != -create_inf () && b == create_inf ())
    intab = 2;
  else if (a == -create_inf () && b != create_inf ())
    intab = 3;
  else
    intab = 4;

  //
  // extract the relative mesh for the eigenfunction
  //
  islfun = ncol - 2;
  slfun = mdr_Create (9 + ncol - 2, 1);
  if (ncol > 2)
  {
    wfunc = mdr_Create (ncol - 2, 2);
    for (i = 2; i < ncol; i++)
    {
      MdrV1 (slfun, 9 + i - 1) = MdrV1 (x5, i);
      Mdr1 (wfunc, i - 1, 1) = MdrV1 (x5, i);
    }
  }
  else
    wfunc = mdr_Create (0, 0);

  //
  // bc, boundary conditions
  //
  e6 = bltin_get_ent (args[6]);
  if (ent_type (e6) != MATRIX_DENSE_REAL)
  {
    fprintf (stdout, THIS_SOLVER ": Wrong boundary conditions matrix.\n");
    rerror (THIS_SOLVER ": Improper boundary conditions!");
  }
  x6 = class_matrix_real (e6);
  if (MNR (x6) != 2 || MNC (x6) != 2)
    rerror (THIS_SOLVER ": Improper boundary conditions!");
  a1 = Mdr1 (x6, 1, 1);
  a2 = Mdr1 (x6, 1, 2);
  b1 = Mdr1 (x6, 2, 1);
  b2 = Mdr1 (x6, 2, 2);

  //
  // eix, the index of the desired eigenvalue
  //
  e7 = bltin_get_ent (args[7]);
  if (ent_type (e7) != MATRIX_DENSE_REAL)
  {
    fprintf (stdout, THIS_SOLVER ": Improper index of the desired eigenvalue.\n");
    rerror (THIS_SOLVER ": Improper eigenvalue index!");
  }
  numeig = (int) class_double (e7);

  //
  // nc, type of boundary points
  //
  if (nargs > 8)
  {
    e8 = bltin_get_ent (args[8]);
    if (ent_type (e8) == MATRIX_DENSE_REAL)
    {
      x8 = class_matrix_real (e8);
      if (SIZE(x8) == 2)
      {
        nca = MdrV1 (x8, 1);
        ncb = MdrV1 (x8, 2);
      }
    }
  }

  //
  // p,q at a,b, can be user specified
  //
  if (nargs > 9)
  {
    e9 = bltin_get_ent (args[9]);
    if (ent_type (e9) == MATRIX_DENSE_REAL)
    {
      x9 = class_matrix_real (e9);
      if (SIZE(x9) == 4)
      {
        p0ata = MdrV1 (x9, 1);
        qfata = MdrV1 (x9, 2);
        p0atb = MdrV1 (x9, 3);
        qfatb = MdrV1 (x9, 4);
      }
    }
  }

  //
  // Set up ENTITIES for user-function. Inc the reference count once, for belonging
  // argument: x, local
  x = mdr_CreateEmpty (1,1);
  ent_x = ent_Assign_Rlab_MDR (x);
  ent_IncRef (ent_x);

  //
  // Figure out p,q,w at the edges  if user has not specified that
  //
  if (p0ata * qfata * p0atb * qfatb == 0)
  {
    // @ a
    if (p_ (&a) == 0)
      p0ata = 1;
    else
      p0ata = -1;
    if (q_ (&a) != create_inf () && q_ (&a) != -create_inf () &&
        w_ (&a) != create_inf () && w_ (&a) != -create_inf () &&
        q_ (&a) != create_nan () && w_ (&a) != create_nan ())
      qfata = 1;
    else
      qfata = -1;
    // @ b
    if (p_ (&b) == 0)
      p0atb = 1;
    else
      p0atb = -1;
    if (q_ (&b) != create_inf () && q_ (&b) != -create_inf () &&
        w_ (&b) != create_inf () && w_ (&b) != -create_inf () &&
        q_ (&b) != create_nan () && w_ (&b) != create_nan ())
      qfatb = 1;
    else
      qfatb = -1;
  }

  if (outs)
    fptr = fopen (outs, "a");
  if (fptr)
  {
    // start the timer
    timer = clock ();
    // write down boring info
    fprintf (fptr,
             "RLaB: using SLEIGN2 integrator for ODE eigenvalue-boundary value problem.\n");
    fprintf (fptr, "RLaB: The original messages from the solver follow.\n");
  }

  //
  // Find eigenvalue
  //
  SLEIGN (&a, &b, &intab, &p0ata, &qfata, &p0atb, &qfatb, &a1, &a2, &b1, &b2,
           &numeig, &eig, &tol, &iflag, &islfun, MDRPTR(slfun), &nca, &ncb,
           outs, &louts);

  if (iflag != 1)
  {
    if (fptr)
    {
      fprintf (fptr, "RLaB: SLEIGN2 integrator reports an error %i.\n", iflag);
      // check the time and close the output stream
      timer -= clock ();
      fprintf (fptr,
               "RLaB: ODE integration using SLEIGN2 lasted %g sec.\n",
               -timer / 1e6);
      fclose (fptr);
    }
    fprintf (stdout,
             THIS_SOLVER ": cannot solve the problem. Try changing the classification\n");
    fprintf (stdout,
             THIS_SOLVER ": of the endpoints or the index of the eigenvalue.\n");
    // return values
    ncol = 0;
    weig = mdr_Create (0, 0);
    rslfun = mdr_Create (0, 0);
  }
  else
  {
    if (fptr)
    {
      fprintf (fptr, "RLaB: SLEIGN2 integrator reports success.\n");
      // check the time and close the output stream
      timer -= clock ();
      fprintf (fptr,
               "RLaB: ODE integration using SLEIGN2 lasted %g sec.\n",
               -timer / 1e6);
      fclose (fptr);
    }

    // return values
    weig = mdr_CreateScalar (eig);
    rslfun = mdr_Create (1, 9);
    for (i = 0; i < 9; i++)
      MdrV0 (rslfun, i) = MdrV0 (slfun, i);
  }
  if (ncol > 2)
  {
    // get the eigenfunction at user's mesh
    mp = MdrV1 (slfun, 1);	// matching point
    for (i = 1; i < ncol - 1; i++)
    {
      if (Mdr1 (wfunc, i, 1) <= mp)
        Mdr1 (wfunc, i, 2) = MdrV1 (slfun, 9 + i);
      else
        break;
    }
    // first derivative should be continuous across the matching point
    sf = 1;
    lfp =
      (MdrV1 (slfun, 9 + i) - MdrV1 (slfun, 9 + i - 1)) / (Mdr1 (wfunc, i, 1) -
        Mdr1 (wfunc, i - 1, 1));
    rfp =
      (MdrV1 (slfun, 9 + i + 1) -
       MdrV1 (slfun, 9 + i)) / (Mdr1 (wfunc, i + 1, 1) - Mdr1 (wfunc, i, 1));
    if (ABS (lfp - rfp) / (ABS(lfp) + ABS (rfp)) > 0.7)
      sf = -1;
    for (j = i; j < ncol - 1; j++)
      Mdr1 (wfunc, j, 2) = sf * MdrV1 (slfun, 9 + j);
  }
  else
    wfunc = mdr_Create (0, 0);

  //
  // Clean Up
  //
  MDPTR (x) = 0;
  ent_DecRef (ent_x);
  ent_Destroy (ent_x);

  mdr_Destroy (slfun);

  ent_Clean (e4); ent_Clean (ent_param);
  ent_Clean (e5);
  ent_Clean (e6);
  ent_Clean (e7);
  ent_Clean (e8);
  ent_Clean (e9);
  ent_Clean (sleignp_fname);
  ent_Clean (sleignq_fname);
  ent_Clean (sleignw_fname);
  ent_Clean (sleignuv_fname);

  //
  // return
  //
  Btree *bw = btree_Create ();
  // slfun
  install (bw, "slfun", ent_Assign_Rlab_MDR (rslfun));
  // eigenvalue
  install (bw, "eigval", ent_Assign_Rlab_MDR (weig));
  //eigenfunction
  install (bw, "eigfun", ent_Assign_Rlab_MDR (wfunc));
  // error message
  install (bw, "flag", ent_Create_Rlab_Double (iflag) );

  // wrap it up
  return ent_Assign_Rlab_BTREE (bw);
}


//
// sturm-liouville coupled boundary value problem: slcoup
//
#undef THIS_SOLVER
#define THIS_SOLVER "odesl.coup"
Ent *
ent_odesl_coup (int nargs, Datum args[])
{
  Ent *e4=0, *e5=0, *e6=0, *e7=0, *e8=0, *e9=0, *e10=0, *e11=0;
  MDR *x5, *x6, *x7, *x10, *x11, *weig, *wfunc, *slfun, *rslfun = 0;

  double a, b, p0ata = 0, qfata = 0, p0atb = 0, qfatb =
    0, a1, a2, b1, b2, tol=-1e-6, mp;
  double lfp, rfp, eig = 0.1;
  double k11, k12, k21, k22, alfa;
  int intab, numeig, iflag, islfun = 0, nca = 1, ncb = 1, i, j, sf;
  int louts=0, ncol;

  FILE *rlab_stderr = (!RLAB_STDERR_DS) ? stderr : RLAB_STDERR_DS;

  char *outs=0;

  //
  // Check arguments.
  //
  if (nargs < 8)
  {
    fprintf (rlab_stderr,
             THIS_SOLVER ": Sturm-Liouville Eigenvalue Boundary Problem solver SLCOUP/SLEIGN2\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": for coupled boundary conditions,\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ":   -(P(x)y'(x))' + Q(x) y(x) = eig W(x) y(x).\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": Format:\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ":   Y = odesl(P,Q,W,/UV/,/p/,x,bc1,bc2,alfa,eix/,nc,atab/),\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": where, P=f(x/,p/), Q=f(x/,p/), and W=f(x/,p/) are the functions,\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": 'p' is the parameter array, x=[a,x1,.,xN, b] gives the boundary\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": points a,b, and the relative mesh for the eigenfunction,\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": bc1=[A1,A2;B1,B2] and bc2=[k11,k12;k21,k22] are the two sets of\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": boundary conditions with 'alfa', UV(x/,p/)=[u(x),u'(x),v(x),v'(x)]\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": gives the boundary behavior of the solution at the endpoints,\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": The index of the eigenvalue is 'eix', and\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": atab=[p0ata,qfata,p0atb,qfatb], nc=[nca,ncb]\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": is the type of the endpoints: 1 (regular), 2 (weakly regular),\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": 3 (limit cycle, non-oscilatory), 4 (limit cycle, oscilatory), 5 \n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": (limit point, regular at finite point), 6 (limit point), 7 (limit\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": point at infinity), 8 (limit point, bad behaved).\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": Y=<<eigval;eigfun;info>>, the eigenvalue, tabulated eigenfunction\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": on user's mesh and the status info from the solver.\n");
    rerror ("requires at least 9 arguments");
  }

  //
  // P(x)
  //
  sleignp_fname = bltin_get_ent(args[0]);
  if (!isfuncent(sleignp_fname))
    rerror (THIS_SOLVER ": " RLAB_ERROR_ARG1_FUNC_VAR "\n");

  //
  // Q(x)
  //
  sleignq_fname = bltin_get_ent(args[1]);
  if (!isfuncent(sleignq_fname))
    rerror (THIS_SOLVER ": " RLAB_ERROR_ARG2_FUNC_VAR "\n");

  //
  // W(x)
  //
  sleignw_fname = bltin_get_ent(args[2]);
  if (!isfuncent(sleignw_fname))
    rerror (THIS_SOLVER ": " RLAB_ERROR_ARG3_FUNC_VAR "\n");

  //
  // UV(x) if given
  //
  sleignuv_fname = bltin_get_ent(args[3]);
  if (!isfuncent(sleignuv_fname))
  {
    ent_Clean (sleignuv_fname);
    sleignuv_fname = 0;
  }


  //
  // parameters for the functions
  //
  e4 = bltin_get_ent (args[4]);
  ent_param = 0;
  if (ent_type (e4) != UNDEF)
    ent_param = ent_Copy (e4);

  //
  // [a,x1..xN,b], the relative mesh on which the of the eigenfunction is computed
  //
  e5 = bltin_get_ent (args[5]);
  if (ent_type (e5) != MATRIX_DENSE_REAL)
  {
    fprintf (stdout,
             THIS_SOLVER ": 5th argument has to be the interval [a,x1,..,xN,b].\n");
    rerror ("Improper bound!");
  }
  x5 = class_matrix_real (e5);
  ncol = MNC (x5);
  if (MNR (x5) != 1 || ncol < 2)
    rerror (THIS_SOLVER ": improper x=[a..,b]!");
  a = MdrV1 (x5, 1);
  b = MdrV1 (x5, ncol);
  if (a != -create_inf () && b != create_inf ())
    intab = 1;
  else if (a != -create_inf () && b == create_inf ())
    intab = 2;
  else if (a == -create_inf () && b != create_inf ())
    intab = 3;
  else
    intab = 4;
  //
  // extract the relative mesh for the eigenfunction
  //
  islfun = ncol - 2;
  slfun = mdr_Create (9 + ncol - 2, 1);
  if (ncol > 2)
  {
    wfunc = mdr_Create (ncol - 2, 2);
    for (i = 2; i < ncol; i++)
    {
      MdrV1 (slfun, 9 + i - 1) = MdrV1 (x5, i);
      Mdr1 (wfunc, i - 1, 1) = MdrV1 (x5, i);
    }
  }
  else
    wfunc = mdr_Create (0, 0);

  //
  // bc1, boundary conditions [a1,a2;b1,b2]
  //
  e6 = bltin_get_ent (args[6]);
  if (ent_type (e6) != MATRIX_DENSE_REAL)
  {
    fprintf (stdout, THIS_SOLVER ": Wrong boundary conditions matrix.\n");
    rerror ("Improper boundary conditions!");
  }
  x6 = class_matrix_real (e6);
  if (MNR (x6) != 2 || MNC (x6) != 2)
    rerror ("Improper boundary conditions!");
  a1 = Mdr1 (x6, 1, 1);
  a2 = Mdr1 (x6, 1, 2);
  b1 = Mdr1 (x6, 2, 1);
  b2 = Mdr1 (x6, 2, 2);
  //
  // bc2, boundary conditions [k11,k12;k21,k22]
  //
  e7 = bltin_get_ent (args[7]);
  if (ent_type (e7) != MATRIX_DENSE_REAL)
  {
    fprintf (stdout, THIS_SOLVER ": Wrong boundary conditions matrix.\n");
    rerror ("Improper boundary conditions!");
  }
  x7 = class_matrix_real (e7);
  if (MNR (x7) != 2 || MNC (x7) != 2)
    rerror ("Improper boundary conditions!");
  k11 = Mdr1 (x7, 1, 1);
  k12 = Mdr1 (x7, 1, 2);
  k21 = Mdr1 (x7, 2, 1);
  k22 = Mdr1 (x7, 2, 2);
  //
  // alfa, enters the boundary conditions [k11,k12;k21,k22]
  //
  e8 = bltin_get_ent (args[8]);
  if (ent_type (e8) != MATRIX_DENSE_REAL)
  {
    fprintf (stdout, THIS_SOLVER ": Improper 'alfa'.\n");
    rerror ("Improper boundary conditions!");
  }
  alfa = class_double (e8);
  //
  // eix, the index of the desired eigenvalue
  //
  e9 = bltin_get_ent (args[9]);
  if (ent_type (e9) != MATRIX_DENSE_REAL)
  {
    fprintf (stdout, THIS_SOLVER ": Improper index of the desired eigenvalue.\n");
    rerror ("Improper eigenvalue index!");
  }
  numeig = (int) class_double (e9);
  //
  // nc, type of boundary points
  //
  if (nargs > 10)
  {
    e10 = bltin_get_ent (args[10]);
    if (ent_type (e10) == MATRIX_DENSE_REAL)
    {
      x10 = class_matrix_real (e10);
      if (MNR (x10) * MNC (x10) == 2)
      {
        nca = MdrV1 (x10, 1);
        ncb = MdrV1 (x10, 2);
      }
    }
  }
  //
  // p,q at a,b, can be user specified
  //
  if (nargs > 11)
  {
    e11 = bltin_get_ent (args[11]);
    if (ent_type (e11) == MATRIX_DENSE_REAL)
    {
      x11 = class_matrix_real (e11);
      if (MNR (x11) * MNC (x11) == 4)
      {
        p0ata = MdrV1 (x11, 1);
        qfata = MdrV1 (x11, 2);
        p0atb = MdrV1 (x11, 3);
        qfatb = MdrV1 (x11, 4);
      }
    }
  }

  //
  // Set up ENTITIES for user-function. Inc the reference count once, for belonging
  // argument: x, local
  x = mdr_CreateEmpty (1,1);
  ent_x = ent_Assign_Rlab_MDR (x);
  ent_IncRef (ent_x);

  //
  // Figure out p,q,w at the edges  if user has not specified that
  //
  if (p0ata * qfata * p0atb * qfatb == 0)
  {
    // @ a
    if (p_ (&a) == 0)
      p0ata = 1;
    else
      p0ata = -1;
    if (q_ (&a) != create_inf () && q_ (&a) != -create_inf () &&
        w_ (&a) != create_inf () && w_ (&a) != -create_inf () &&
        q_ (&a) != create_nan () && w_ (&a) != create_nan ())
      qfata = 1;
    else
      qfata = -1;
    // @ b
    if (p_ (&b) == 0)
      p0atb = 1;
    else
      p0atb = -1;
    if (q_ (&b) != create_inf () && q_ (&b) != -create_inf () &&
        w_ (&b) != create_inf () && w_ (&b) != -create_inf () &&
        q_ (&b) != create_nan () && w_ (&b) != create_nan ())
      qfatb = 1;
    else
      qfatb = -1;
  }

  //
  // Find eigenvalue
  //
  SLCOUP (&a, &b, &intab, &p0ata, &qfata, &p0atb, &qfatb, &a1, &a2, &b1, &b2,
          &numeig, &eig, &tol, &iflag, &islfun, MDRPTR(slfun), &nca, &ncb,
          &alfa, &k11, &k12, &k21, &k22,
          outs, &louts);

  if (iflag != 1 && iflag != 2)
  {
    fprintf (stdout,
             "odesl: cannot solve the problem. Try changing the classification of the\n");
    fprintf (stdout, "odesl: endpoints or the index of the eigenvalue.\n");
    // return values
    ncol = 0;
    wfunc = mdr_Create (0, 0);
    weig = mdr_Create (0, 0);
    rslfun = mdr_Create (0, 0);
  }
  else
  {
    weig = mdr_CreateScalar (eig);
  }

  if (ncol > 2)
  {
    // prepare slfun for return
    rslfun = mdr_Create (1, 9);
    for (i = 0; i < 9; i++)
      MdrV0 (rslfun, i) = MdrV0 (slfun, i);
    // get the eigenfunction at user's mesh
    // get the eigenfunction at user's mesh
    mp = MdrV0 (slfun, 0);	// matching point
    for (i = 1; i < ncol - 1; i++)
    {
      if (Mdr1 (wfunc, i, 1) <= mp)
        Mdr1 (wfunc, i, 2) = MdrV1 (slfun, 9 + i);
      else
        break;
    }
    // first derivative should be continuous across the matching point
    sf = 1;
    lfp =
      (MdrV1 (slfun, 9 + i) - MdrV1 (slfun, 9 + i - 1)) / (Mdr1 (wfunc, i, 1) -
        Mdr1 (wfunc, i - 1, 1));
    rfp =
      (MdrV1 (slfun, 9 + i + 1) -
       MdrV1 (slfun, 9 + i)) / (Mdr1 (wfunc, i + 1, 1) - Mdr1 (wfunc, i, 1));
    if (ABS (lfp - rfp) > 10)
      sf = -1;
    for (j = i; j < ncol - 1; j++)
      Mdr1 (wfunc, j, 2) = sf * MdrV1 (slfun, 9 + j);
  }

  //
  // Clean Up
  //
  //ent_DecRef (ent_x); x->d = 0; ent_Destroy (ent_x);
  MDPTR (x) = 0;
  ent_DecRef (ent_x);
  ent_Destroy (ent_x);

  mdr_Destroy (slfun);

  ent_Clean (ent_param);

  ent_Clean (e4);
  ent_Clean (e5);
  ent_Clean (e6);
  ent_Clean (e7);
  ent_Clean (e8);
  ent_Clean (e9);
  ent_Clean (e10);
  ent_Clean (e11);
  ent_Clean (sleignp_fname);
  ent_Clean (sleignq_fname);
  ent_Clean (sleignw_fname);
  ent_Clean (sleignuv_fname);

  //
  // return
  //
  Btree *bw = btree_Create ();
  // slfun
  install (bw, "slfun", ent_Assign_Rlab_MDR (rslfun));
  // eigenvalue
  install (bw, "eigval", ent_Assign_Rlab_MDR (weig));
  //eigenfunction
  install (bw, "eigfun", ent_Assign_Rlab_MDR (wfunc));
  // error message
  install (bw, "flag", ent_Create_Rlab_Double (iflag) );

  // wrap it up
  return ent_Assign_Rlab_BTREE (bw);
}


//
// The interface to the user-specified function in fortran
//
double
p_ (double *newx)
{
  Ent *rent = 0;
  double rval=0;

  MDPTR (x) = (void *) newx;

  if (ent_param)
    rent = ent_call_rlab_script_2args(sleignp_fname, ent_x, ent_param);
  else
    rent = ent_call_rlab_script_1arg (sleignp_fname, ent_x);

  if (ent_type(rent)!=MATRIX_DENSE_REAL)
    rerror (THIS_SOLVER ": [p]" RLAB_ERROR_RHS_FUNC_MUST_RETURN_MDR);

  rval = class_double (rent);
  ent_Clean (rent);
  return rval;
}

//
// The interface to the user-specified function in fortran
//
double
q_ (double *newx)
{
  Ent *rent = 0;
  double rval=0;

  MDPTR (x) = (void *) newx;

  if (ent_param)
    rent = ent_call_rlab_script_2args(sleignq_fname, ent_x, ent_param);
  else
    rent = ent_call_rlab_script_1arg (sleignq_fname, ent_x);

  if (ent_type(rent)!=MATRIX_DENSE_REAL)
    rerror (THIS_SOLVER ": [q]" RLAB_ERROR_RHS_FUNC_MUST_RETURN_MDR);

  rval = class_double (rent);
  ent_Clean (rent);
  return rval;
}


//
// The interface to the user-specified function in fortran
//
double
w_ (double *newx)
{
  Ent *rent = 0;
  double rval=0;

  MDPTR (x) = (void *) newx;

  if (ent_param)
    rent = ent_call_rlab_script_2args(sleignw_fname, ent_x, ent_param);
  else
    rent = ent_call_rlab_script_1arg (sleignw_fname, ent_x);

  if (ent_type(rent)!=MATRIX_DENSE_REAL)
    rerror (THIS_SOLVER ": [w]" RLAB_ERROR_RHS_FUNC_MUST_RETURN_MDR);

  rval = class_double (rent);
  ent_Clean (rent);
  return rval;
}


//
// The interface to the user-specified function in fortran
//
int
uv_ (double *newx, double *u, double *pup,
     double *v, double *pvp, double *hu, double *hv)
{
  Ent *rent = 0;
  MDR *retm = 0;

  //
  // u is needed only for special boundary conditions
  //
  if (!sleignuv_fname)
    return 1.0;

  MDPTR (x) = (void *) newx;

  if (ent_param)
    rent = ent_call_rlab_script_2args(sleignuv_fname, ent_x, ent_param);
  else
    rent = ent_call_rlab_script_1arg (sleignuv_fname, ent_x);

  if (ent_type(rent)!=MATRIX_DENSE_REAL)
    rerror (THIS_SOLVER ": [uv]" RLAB_ERROR_RHS_FUNC_MUST_RETURN_MDR);

  retm = ent_data (rent);

  if (SIZE(retm) != 4)
    rerror (THIS_SOLVER ": [uv]" RLAB_ERROR_RHS_FUNC_MUST_RETURN_MDR);

  *u    = MdrV0 (retm, 0);
  *pup  = MdrV0 (retm, 1);
  *v    = MdrV0 (retm, 2);
  *pvp  = MdrV1 (retm, 3);
  *hu = 0;
  *hv = 0;

  ent_Clean (rent);
  return 1;
}
