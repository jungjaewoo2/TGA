// Copyright (C) 2005-2008 Marijan Kostrun
// part of rlabplus for linux project on rlabplus.sourceforge.net
//
// rlabplus <-> EPDCOL a.k.a. TOMS 688.f, BACOL
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

#include "pdecol.h"

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
// pdecol: global entities
//
static int pdecol_npde, pdecol_mf, pdecol_nint;
static double xleft, xright;
static int pde_df = 0, pde_dbl = 0, pde_dbr = 0;

static MDR *t, *x, *u, *ux, *uxx, *xbkpt;
static Ent *ent_t=0, *ent_x=0, *ent_u=0, *ent_ux=0, *ent_uxx=0;

static Ent *pdecol_u0_name;
static Ent *pdecol_f_name;
static Ent *pdecol_fjac_name;
static Ent *bacol_bleft_name;
static Ent *bacol_bright_name;
static Ent *pdecol_dbleft_name;
static Ent *pdecol_dbright_name;

#undef  THIS_SOLVER
#define THIS_SOLVER "pdecol"
Ent *
ent_pdecol (int nargs, Datum args[])
{
  Ent *e2=0, *e7=0, *e9=0, *eo=0;
  ListNode *node;
  MDR *timeT, *w;
  char *outs=0;
  FILE *fptr = NULL;
  time_t t1=clock(), t2=clock();
  int louts = 0, ncalls, i;
  double timer;

  //
  // init default pdecol parameters
  //
  int kord = 4;    // the order of piecewise polynomial to be used 2<=kord<=21
  int ncc = 2;     // the number of continuity conditions at the breakpoint, 2<=ncc<=kord
  int meth = 2;    // 1 adams method, 2 backward differentiation formula
  int miter = 1;   // 1 chord method with analytic jacobian,
                   // 2 chord method with finite differences jacobian
  int nindex = 1;  // 1 first call,
                   // 0 continue integration
                   // 2 same as 0, but hit tout exactly
                   // 3,4 internal
  int nogauss = 0;
  int maxder = 5;
  double dt=1e-4, t0, tout, erel = 1e-4;
  int lwork, liwork;
  MDR *work=0, *iwork=0, *w1=0, *w2=0, *scratch=0;
  int ncpts, idummy, nintp1;
  double ddummy;

  pde_df = 0;
  pde_dbl = 0;
  pde_dbr = 0;

  //
  // Load and Check arguments.
  //
  if (nargs < 8 || nargs > 9)
  {
    fprintf (stdout,
             THIS_SOLVER ": Integrates a system of partial differential equations in time.\n");
    fprintf (stdout,
             THIS_SOLVER ": Format:\n");
    fprintf (stdout,
             THIS_SOLVER ":   u(tf,x) = pdecol(npde,f,/Df/,DbL,DbR,x,u0,T/,options/),\n");
    fprintf (stdout,
             THIS_SOLVER ": where  f=function(t,x,u,ux,uxx), <<dfdu;dfdux;dfduxx>>=\n");
    fprintf (stdout,
             THIS_SOLVER ": Df=function(t,x,u,ux,uxx), if given, or otherwise has to\n");
    fprintf (stdout,
             THIS_SOLVER ": be explicitely omitted, <<dbdu;dbdux;dzdt>>=DbL,R=\n");
    fprintf (stdout,
             THIS_SOLVER ": function(t,u,ux) at xL,R  defines the boundary condition\n");
    fprintf (stdout,
             THIS_SOLVER ": B(u,ux)=Z(t). 'x' is a mesh used for discretizetion of the\n");
    fprintf (stdout,
             THIS_SOLVER ": spatial coordinate, while  u0=function(x) defines the initial\n");
    fprintf (stdout,
             THIS_SOLVER ": condition on function u. T=[ti,tf] is the integration time.\n");
    fprintf (stdout,
             THIS_SOLVER ": The list  options=<<nogauss;kord;erel;dt;ncc;iode;maxder;\n");
    fprintf (stdout,
             THIS_SOLVER ": stdout>>  controls the integrator. See manual for details.\n");
    fprintf (stdout,
             THIS_SOLVER ": \n");
    rerror ("requires at least 7 arguments");
  }

  //
  // npde
  //
  e2 = bltin_get_ent (args[0]);
  if (ent_type (e2) != MATRIX_DENSE_REAL)
    rerror ("pdecol: 'npde' has to be a scalar");
  pdecol_npde = (int) class_double (e2);

  //
  // Get function ptrs
  //
  // D1
  pdecol_f_name = bltin_get_ent(args[1]);
  if (!isfuncent(pdecol_f_name))
    rerror (THIS_SOLVER ": " RLAB_ERROR_ARG2_FUNC_VAR "\n");

  // D2
  pdecol_fjac_name = bltin_get_ent(args[2]);
  if (!isfuncent(pdecol_fjac_name) && ent_type(pdecol_fjac_name)!=UNDEF)
    rerror (THIS_SOLVER ": " RLAB_ERROR_ARG3_FUNC_VAR "\n");
  if (isfuncent(pdecol_fjac_name))
    miter = 1;

  // D3
  pdecol_dbleft_name = bltin_get_ent(args[3]);
  if (!isfuncent(pdecol_dbleft_name))
    rerror (THIS_SOLVER ": " RLAB_ERROR_ARG4_FUNC_VAR "\n");

  // D4
  pdecol_dbright_name = bltin_get_ent(args[4]);
  if (!isfuncent(pdecol_dbright_name))
    rerror (THIS_SOLVER ": " RLAB_ERROR_ARG5_FUNC_VAR "\n");

  //
  // x
  //
  e7 = bltin_get_ent (args[5]);
  if (ent_type (e7) != MATRIX_DENSE_REAL)
    rerror ("pdecol: 'x' has to be a column vector");
  xbkpt = class_matrix_real (e7);
  if (!xbkpt)
    rerror ("pdecol: 'x' has to be a column vector");
  pdecol_nint   = MNR(xbkpt) * MNC(xbkpt) - 1;
  xleft  = MdrV1(xbkpt, 1);
  xright = MdrV1(xbkpt, pdecol_nint + 1);

  // D6
  pdecol_u0_name = bltin_get_ent(args[6]);
  if (!isfuncent(pdecol_u0_name))
    rerror (THIS_SOLVER ": " RLAB_ERROR_ARG7_FUNC_VAR "\n");

  //
  // Get variables that may change between calls
  //
  // T = [t1, t2, ..]
  e9 = bltin_get_ent (args[7]);
  if (ent_type (e9) != MATRIX_DENSE_REAL)
    rerror ("pdecol: Integration times row-vector must be T=[t1,t2]");
  timeT = class_matrix_real (e9);
  ncalls = SIZE(timeT) - 1;
  if (ncalls < 1)
    rerror ("pdecol: Integration times row-vector must be T=[t1,t2]");

  //
  // options for the solver
  //
  if (nargs > 8)
  {
    eo = bltin_get_ent (args[8]);
    if (ent_type (eo) == BTREE)
    {
      // kord
      node = btree_FindNode (ent_data (eo), RLAB_NAME_PDECOL_KORD);
      if (node != 0)
      {
        idummy = class_double (var_ent (node));
        if (idummy > 0.0)
          kord = idummy;
      }
      // erel
      node = btree_FindNode (ent_data (eo), RLAB_NAME_PDECOL_EREL);
      if (node != 0)
      {
        ddummy =  class_double (var_ent (node));
        if (ddummy == -1 || ddummy > 0)
          if (nindex == 1) erel = ddummy;
      }
      // dt
      node = btree_FindNode (ent_data (eo), RLAB_NAME_PDECOL_DT);
      if (node != 0)
      {
        ddummy =  class_double (var_ent (node));
        if (ddummy > 0.0)
          if (nindex == 1)  dt = ddummy;
      }

      // ncc
      node = btree_FindNode (ent_data (eo), RLAB_NAME_PDECOL_NCC);
      if (node != 0)
      {
        ddummy =  class_double (var_ent (node));
        if (ddummy == -1 || ddummy > 0)
          ncc = ddummy;
      }
      // iode: how is ode integrated
      node = btree_FindNode (ent_data (eo), RLAB_NAME_PDECOL_IODE);
      if (node != 0)
      {
        idummy =  (int) class_double (var_ent (node));
        if (idummy == 1 || idummy == 2)
          meth = idummy;
      }
      // maxder
      node = btree_FindNode (ent_data (eo), RLAB_NAME_PDECOL_MAXDER);
      if (node != 0)
      {
        idummy =  (int) class_double (var_ent (node));
        if (idummy > 0 || idummy < 10)
          maxder = idummy;
      }
      // nogauss
      node = btree_FindNode (ent_data (eo), RLAB_NAME_PDECOL_NOGAUSS);
      if (node != 0)
      {
        idummy =  (int) class_double (var_ent (node));
        if (idummy == 0 || idummy == 1)
          nogauss = idummy;
      }
      // stdout
      node = btree_FindNode (ent_data (eo), RLAB_NAME_PDECOL_STDOUT);
      if (node != 0)
      {
        if (ent_type(var_ent (node))== MATRIX_DENSE_STRING)
        {
          outs = class_char_pointer(var_ent (node));
        }
        louts = isvalidstring (outs);
        if (louts < 1)
        {
          louts = 0;
          outs = 0;
        }
      }
    }
    ent_Clean (eo);
  }

  //
  // Prepare rlab function call for PDE1F and PDE1FJAC
  // = function(t,x,u,ux,uxx)
  //

  // time t
  t = mdr_CreateEmpty (1, 1);
  ent_t = ent_Assign_Rlab_MDR (t);
  ent_IncRef (ent_t);

  // position x
  x = mdr_CreateEmpty (1, 1);
  ent_x = ent_Assign_Rlab_MDR (x);
  ent_IncRef (ent_x);

  // position u
  u = mdr_CreateEmpty (pdecol_npde,1);
  ent_u = ent_Assign_Rlab_MDR (u);
  ent_IncRef (ent_u);

  // position ux
  ux = mdr_CreateEmpty (1,pdecol_npde);
  ent_ux = ent_Assign_Rlab_MDR (ux);
  ent_IncRef (ent_ux);

  // position uxx
  uxx = mdr_CreateEmpty (1,pdecol_npde);
  ent_uxx = ent_Assign_Rlab_MDR (uxx);
  ent_IncRef (ent_uxx);

  //
  // Prepare rlab function call for PDE1BND
  // = function(t,u,ux) at xleft, xright
  //
  // time t
  // function u
  // position ux

  //
  // Prepare rlab function call for PDE1UINT
  // = function(x)
  //
  // position x

  //
  // start initialization of variables
  //
  pdecol_mf = 10 * meth + miter;
  ncpts = kord*pdecol_nint - ncc*(pdecol_nint-1);

  // lwork and liwork
  lwork = kord +
      4*pdecol_npde + ncpts*(3*kord + 2) +
      pdecol_npde*ncpts*(maxder + 6) +
      (pdecol_npde*pdecol_npde)*(13+kord*(kord-ncc)*pdecol_nint);
  liwork = ncpts * (pdecol_npde + 1);
  work   = mdr_Create(lwork, 1); mdr_Zero (work);
  iwork  = mdi_Create(liwork, 1); mdr_Zero (iwork);
  MdiV1(iwork,1) = lwork;
  MdiV1(iwork,2) = liwork;


  if (outs)
    fptr = fopen (outs, "a");
  if (fptr)
  {
    fprintf (fptr, "RLaB: Using PDECOL time integrator. Messages from the solver follow.\n");
  }

  w = mdr_Create(0,0);
  nintp1  = pdecol_nint+1;
  w1      = mdr_Create (nintp1, pdecol_npde);
  scratch = mdr_Create (nintp1, 1);

  for (i = 1; i <= ncalls; i++)
  {
    t0   = MdrV1(timeT, i);

    //
    // integrate
    //

    tout = MdrV1(timeT, i+1);

    if (fptr)
      fprintf (fptr,
               "RLaB: (PDECOL) time evolution parameters t0, tout, dt = %g, %g, %g\n",
               t0, tout, dt);

    PDECOL(&t0, &tout, &dt,
           MDRPTR(xbkpt), &erel, &pdecol_nint, &kord,
           &ncc, &pdecol_npde, &pdecol_mf, &nindex,
           MDRPTR(work), MDIPTR(iwork),
           &louts, outs, &maxder, &nogauss);

    if (nindex >= 0)
    {
      PDE1VAL (MDRPTR(xbkpt), MDRPTR(w1), MDRPTR(scratch), &pdecol_npde, &nintp1,
               &nintp1, MDRPTR(work));

      // perform
      //  w <- mdr_Append(w, w1);
      w2 = mdr_Append(w, w1); mdr_Destroy(w);

      // now move w2 to w
      w = mdr_CreateEmpty(w2->nrow, w2->ncol);
      MDPTR(w) = MDPTR(w2);
      MDPTR(w2) = 0;
      mdr_Destroy(w2);

      continue;
    }
    break;
  }

  mdr_Destroy( scratch );
  mdr_Destroy( w1 );

  //
  // write the report from the solver
  //
  if (fptr)
  {
    fprintf (fptr,
             "RLaB: PDECOL time integrator reports %i ", nindex);
    switch (nindex)
    {
      case 0:
        fprintf (fptr, "(success) !\n");
        break;
      case -1:
        fprintf (fptr, "(failed to pass error test) !\n");
        break;
      case -2:
        fprintf (fptr, "(erel - requested accuracy too high) !\n");
        break;
      case -3:
        fprintf (fptr, "(corr - requested accuracy too high) !\n");
        break;
      default:
        fprintf (fptr, "(inconsistent input parameter(s)) !\n");
    }
      // check the time and close the output stream
    t2 = clock ();
    timer = (t2 - t1) / 1e6;
    fprintf (fptr, "RLaB: Integration lasted %g sec.\n", timer);
    fclose (fptr);
  }


  //
  // Cleanup
  //
  mdr_Destroy (work);
  mdr_Destroy (iwork);

  // ent_t
  MDPTR(t) = 0;
  ent_DecRef (ent_t);
  ent_Destroy (ent_t);

  // ent_x
  MDPTR(x) = 0;
  ent_DecRef (ent_x);
  ent_Destroy (ent_x);
  // ent_u
  MDPTR(u) = 0;
  ent_DecRef (ent_u);
  ent_Destroy (ent_u);
  // ent_ux
  MDPTR(ux) = 0;
  ent_DecRef (ent_ux);
  ent_Destroy (ent_ux);
  // ent_uxx
  MDPTR(uxx) = 0;
  ent_DecRef (ent_uxx);
  ent_Destroy (ent_uxx);

  //
  // entities used to store constant df, dbl, dbr
  //
  ent_Clean (pdecol_u0_name);
  ent_Clean (pdecol_f_name);
  ent_Clean (pdecol_fjac_name);
  ent_Clean (pdecol_dbleft_name);
  ent_Clean (pdecol_dbright_name);

  ent_Clean (e2);
  ent_Clean (e7);
  ent_Clean (e9);

  return ent_Assign_Rlab_MDR(w);
}

#undef  THIS_SOLVER
#define THIS_SOLVER "bacol"
Ent *
ent_bacol (int nargs, Datum args[])
{
  Ent *e2=0, *e7=0, *e9=0, *eo=0;
  ListNode *node;
  MDR *timeT, *w;
  char *outs=0;
  FILE *fptr = NULL;
  time_t t1=clock(), t2=clock();
  int louts = 0, ncalls, i;
  double timer;

  //
  // init default pdecol parameters
  //
  int kord = 4;    // the order of piecewise polynomial to be used 2<=kord<=21
  int ncc = 2;     // the number of continuity conditions at the breakpoint, 2<=ncc<=kord
  int idummy, nintp1, idir = 0, pdecol_nintmx = 10000, isover = 1;
  int kcol, lrp, lip, idid;
  double t0, tout, eabs = 1e-3, erel = 1e-2, h0 = 0.0;
  MDR *w1=0, *w2=0, *scratch=0, *mflag=0, *atol=0, *rtol=0, *xi=0, *rpar=0, *ipar=0, *y=0;

  double ddummy;

  // prepare to memorize constant jacobians
  pde_df = 0;
  pde_dbl = 0;
  pde_dbr = 0;

  //
  // Load and Check arguments.
  //
  if (nargs < 10 || nargs > 11)
  {
    fprintf (stdout,
             "bacol: Integrates a system of partial differential equations in time.\n");
    fprintf (stdout,
             "bacol: Format:\n");
    fprintf (stdout,
             "bacol:   u(tf,x) = bacol(npde,f,Df,bL,bR,dbL,dbR,x,u0,T/,options/),\n");
    fprintf (stdout,
             "bacol: where  f=function(t,x,u,ux,uxx), <<dfdu;dfdux;dfduxx>>=\n");
    fprintf (stdout,
             "bacol: Df=function(t,x,u,ux,uxx), is its jacobian; boundary conditions\n");
    fprintf (stdout,
             "bacol: are given by bL,R=0, where bL,R=function(t,u,ux), and\n");
    fprintf (stdout,
             "bacol: dbL,R=function(t,u,ux)=<<dbdu;dbdux;dbdt>> are their jacobians;\n");
    fprintf (stdout,
             "bacol: 'x' is a mesh used for discretizetion of the spatial coordinate,\n");
    fprintf (stdout,
             "bacol: while  u0=function(x)  defines the initial condition for u=u(x,t);\n");
    fprintf (stdout,
             "bacol: T=[ti,..,tf] are the integration times, while\n");
    fprintf (stdout,
             "bacol: options=<<kord;eabs;erel;stdout>>  controls the integrator.\n");
    fprintf (stdout,
             "bacol: \n");
    rerror ("requires at least 10 arguments");
  }

  //
  // npde
  //
  e2 = bltin_get_ent (args[0]);
  if (ent_type (e2) == MATRIX_DENSE_REAL)
    pdecol_npde = (int) class_double (e2);
  else
    rerror ("pdecol: 'npde' has to be a scalar");

  //
  // Get function ptrs
  //
  //
  // Get function ptrs
  //
  // D1
  pdecol_f_name = bltin_get_ent(args[1]);
  if (!isfuncent(pdecol_f_name))
    rerror (THIS_SOLVER ": " RLAB_ERROR_ARG2_FUNC_VAR "\n");

  // D2
  pdecol_fjac_name = bltin_get_ent(args[2]);
  if (!isfuncent(pdecol_fjac_name) && ent_type(pdecol_fjac_name)!=UNDEF)
    rerror (THIS_SOLVER ": " RLAB_ERROR_ARG3_FUNC_VAR "\n");

  // D3
  bacol_bleft_name = bltin_get_ent(args[3]);
  if (!isfuncent(bacol_bleft_name))
    rerror (THIS_SOLVER ": " RLAB_ERROR_ARG4_FUNC_VAR "\n");

  // D4
  bacol_bright_name = bltin_get_ent(args[4]);
  if (!isfuncent(bacol_bright_name))
    rerror (THIS_SOLVER ": " RLAB_ERROR_ARG5_FUNC_VAR "\n");

  // D5
  pdecol_dbleft_name = bltin_get_ent(args[5]);
  if (!isfuncent(pdecol_dbleft_name))
    rerror (THIS_SOLVER ": " RLAB_ERROR_ARG6_FUNC_VAR "\n");

  // D6
  pdecol_dbright_name = bltin_get_ent(args[6]);
  if (!isfuncent(pdecol_dbright_name))
    rerror (THIS_SOLVER ": " RLAB_ERROR_ARG7_FUNC_VAR "\n");

  //
  // x
  //
  e7 = bltin_get_ent (args[7]);
  if (ent_type (e7) == MATRIX_DENSE_REAL)
    xbkpt = class_matrix_real (e7);
  else
    rerror ("bacol: 'x' has to be a column vector");
  pdecol_nint   = MNR(xbkpt) * MNC(xbkpt) - 1;
  xleft  = MdrV1(xbkpt, 1);
  xright = MdrV1(xbkpt, pdecol_nint + 1);

  // D7
  pdecol_u0_name = bltin_get_ent(args[8]);
  if (!isfuncent(pdecol_u0_name))
    rerror (THIS_SOLVER ": " RLAB_ERROR_ARG9_FUNC_VAR "\n");

  //
  // Get variables that may change between calls
  //
  // T = [t1, t2, ..]
  e9 = bltin_get_ent (args[9]);
  if (ent_type (e9) != MATRIX_DENSE_REAL)
    rerror ("bacol: Integration times row-vector must be T=[t1,t2]");
  timeT = class_matrix_real (e9);
  ncalls = SIZE(timeT) - 1;
  if (ncalls < 1)
    rerror ("bacol: Integration times row-vector must be T=[t1,t2]");

  //
  // options for the solver
  //
  if (nargs > 10)
  {
    eo = bltin_get_ent (args[10]);
    if (ent_type (eo) == BTREE)
    {
      // kord
      node = btree_FindNode (ent_data (eo), RLAB_NAME_BACOL_KORD);
      if (node != 0)
      {
        idummy = class_double (var_ent (node));
        if (idummy > 0.0)
          kord = idummy;
      }
      // dirichlet
      node = btree_FindNode (ent_data (eo), RLAB_NAME_BACOL_IDIR);
      if (node != 0)
      {
        idummy = class_double (var_ent (node));
        if (idummy == 0 || idummy == 1)
          idir = idummy;
      }
      // nintmax
      node = btree_FindNode (ent_data (eo), RLAB_NAME_BACOL_IMAX);
      if (node != 0)
      {
        idummy = class_double (var_ent (node));
        if (idummy > 0)
          pdecol_nintmx = idummy;
      }
      // isover
      node = btree_FindNode (ent_data (eo), RLAB_NAME_BACOL_ISOVER);
      if (node != 0)
      {
        idummy = class_double (var_ent (node));
        if (idummy == 0 || idummy == 1)
          isover = idummy;
      }
      // dt
      node = btree_FindNode (ent_data (eo), RLAB_NAME_BACOL_DT);
      if (node != 0)
      {
        ddummy =  class_double (var_ent (node));
        if (ddummy > 0.0 || ddummy < 1.0)
          h0 = ddummy;
      }
      // atol
      node = btree_FindNode (ent_data (eo), RLAB_NAME_BACOL_EREL);
      if (node != 0)
      {
        if (ent_type(var_ent(node))==MATRIX_DENSE_REAL)
        {
          MDR *rdummy =  class_matrix_real (var_ent (node));
          if (MNR(rdummy)*MNC(rdummy)==1 && MNR(rdummy)*MNC(rdummy)==pdecol_npde)
            rtol = mdr_Copy(rdummy);
        }
      }
      // rtol
      node = btree_FindNode (ent_data (eo), RLAB_NAME_BACOL_EABS);
      if (node != 0)
      {
        if (ent_type(var_ent(node))==MATRIX_DENSE_REAL)
        {
          MDR *rdummy =  class_matrix_real (var_ent (node));
          if (MNR(rdummy)*MNC(rdummy)==1 && MNR(rdummy)*MNC(rdummy)==pdecol_npde)
            atol = mdr_Copy(rdummy);
        }
      }
      // stdout
      node = btree_FindNode (ent_data (eo), RLAB_NAME_BACOL_STDOUT);
      if (node != 0)
      {
        if (ent_type(var_ent (node))== MATRIX_DENSE_STRING)
        {
          outs = class_char_pointer(var_ent (node));
        }
        if (outs)
          louts = strlen(outs);
      }
    }
    ent_Clean (eo);
  }

  //
  // Prepare rlab function call for PDE1F and PDE1FJAC
  // = function(t,x,u,ux,uxx)
  //
  // time t
  t = mdr_CreateEmpty (1, 1);
  ent_t = ent_Assign_Rlab_MDR (t);
  ent_IncRef (ent_t);

  // position x
  x = mdr_CreateEmpty (1, 1);
  ent_x = ent_Assign_Rlab_MDR (x);
  ent_IncRef (ent_x);

  // position u
  u = mdr_CreateEmpty (pdecol_npde,1);
  ent_u = ent_Assign_Rlab_MDR (u);
  ent_IncRef (ent_u);

  // position ux
  ux = mdr_CreateEmpty (1,pdecol_npde);
  ent_ux = ent_Assign_Rlab_MDR (ux);
  ent_IncRef (ent_ux);


  // position uxx
  uxx = mdr_CreateEmpty (1,pdecol_npde);
  ent_uxx = ent_Assign_Rlab_MDR (uxx);
  ent_IncRef (ent_uxx);

  //
  // Prepare rlab function call for BAC1BXA,BAC1BXB,BAC1DBXA,BAC1DBXB
  // = function(t,u,ux) at xleft, xright
  //
  // function u
  // position ux

  //
  // Prepare rlab function call for PDE1UINT
  // = function(x)
  //
  // position x

  //
  // bacol
  //
  kcol = kord - 2;
  idid=0;

  nintp1 = pdecol_nint+1;

  xi = mdr_Create(pdecol_nintmx+1, 1);
  for (i=0; i < pdecol_nint+1; i++)
    MdrV0 (xi,i) = MdrV0 (xbkpt,i);

  lrp = 134 + pdecol_nintmx * (35 + 35*kcol+31*pdecol_npde
      +38*pdecol_npde*kcol+8*kcol*kcol) + 14*kcol
      + 79*pdecol_npde+pdecol_npde*pdecol_npde*     (21
      +4*pdecol_nintmx*kcol*kcol+12*pdecol_nintmx*kcol
      +6*pdecol_nintmx);
  rpar = mdr_Create(lrp,1);
  y = mdr_Create (pdecol_npde*(kcol*pdecol_nintmx+ncc), 1);

  lip = 115 + pdecol_npde*(pdecol_nintmx*(2*kcol+1)+4);
  ipar = mdi_Create(lip, 1);

  w = mdr_Create(0,0);
  w1      = mdr_Create (nintp1, pdecol_npde);
  scratch = mdr_Create ((kcol+ncc)+kcol*(pdecol_nintmx+1)+2*ncc,1);


  if (outs)
    fptr = fopen (outs, "a");
  if (fptr)
  {
    // start the timer
    t1 = clock ();
    fprintf (fptr, "RLaB: Using BACOL time integrator. Messages from the solver follow.\n");
  }

  mflag = mdi_Create(6, 1);
  // perform initialization first:
  MdiV1(mflag,1) = 0;
  // 1, enforce code to stop at tout, 0 don't
  if (isover)
  {
    MdiV1(mflag,3) = 1;
    MdrV1(rpar,1) = tout;
  }
  else
    MdiV1(mflag,3) = 1;
  // stop at tout, rather than at a predefined number of steps:
  MdiV1(mflag,4) = 0;
  // is at least one of the boundary conditions a dirichlet one:
  MdiV1(mflag,5) = idir;  //
  // initial stepsize:
  if (h0 > 0.0)
  {
    MdrV1(rpar, 2) = h0;
    MdiV1(mflag,6) = 1;
  }
  else
    MdiV1(mflag,6) = 0;

  //
  // set absolute and relative errors if user did not so
  //
  if (!atol)
  {
    atol = mdr_Create (1,1);
    MdrV1(atol,1) = eabs;
  }
  if (!rtol)
  {
    rtol = mdr_Create (1,1);
    MdrV1(rtol,1) = erel;
  }
  if (MNR(atol)*MNC(atol)!=MNR(rtol)*MNC(rtol) || MNR(atol)*MNC(atol)==1)
  {
    // use scalar relative and absolute tolerances
    MdiV1(mflag,2) = 0;
  }
  else
  {
    // atol and rtol are of size pdecol_npde
    MdiV1(mflag,2) = 1;
  }

  t0   = MdrV1(timeT, 1);

  for (i = 1; i <= ncalls; i++)
  {
    //
    // Integrate
    //
    tout = MdrV1(timeT, i+1);
    if (isover) MdrV1(rpar,1) = tout;
    if (fptr)
    {
      fprintf (fptr,"RLaB: (BACOL) time evolution parameters t0, tout, dt = ");
      fprintf (fptr,"%g,",  t0);
      fprintf (fptr,"%g,",  tout);
      fprintf (fptr,"%g\n", tout-t0);
    }

    BACOL (&t0, &tout, MDRPTR(atol), MDRPTR(rtol), &pdecol_npde, &kcol, &pdecol_nintmx,
            &pdecol_nint, MDRPTR(xi), MDIPTR(mflag),
            MDRPTR(rpar), &lrp, MDIPTR(ipar), &lip, MDRPTR(y), &idid,
           &louts, outs);

    if (idid >= 0)
    {
      BAC1VAL (&kcol, MDRPTR(xbkpt), &pdecol_nint, MDRPTR(xi), &pdecol_npde,
                &nintp1, MDRPTR(w1), MDRPTR(y), MDRPTR(scratch));
      // perform
      //  w = mdr_Append(w, w1);
      w2 = mdr_Append(w, w1);
      mdr_Destroy(w);
      // now move w2 to w
      w = mdr_Create(0,0);
      w->nrow = w2->nrow;
      w->ncol = w2->ncol;
      MDPTR(w) = MDPTR(w2);
      MDPTR(w2) = 0;
      mdr_Destroy(w2);

      idid = 1;
      MdiV1(mflag,1) = 1; // continue computation
      if (i == ncalls) isover = 0; // enforce the integration to stop at last t
      continue ;
    }
    break;
  }

  //
  // write the report from the solver
  //
  if (fptr)
  {
    fprintf (fptr,
             "RLaB: BACOL time integrator reports %i ", idid);
    switch (idid)
    {
      case 1:
      case 2:
      case 3:
        fprintf (fptr, "(success) !\n");
        break;
      case -1:
        fprintf (fptr, "(Computation takes too many iterations) !\n");
        break;
      case -2:
        fprintf (fptr, "(The error tolerances are too stringent) !\n");
        break;
      case -3:
        fprintf (fptr, "(Zero component in ATOL) !\n");
        break;
      case -6:
        fprintf (fptr, "(Repeated DASSL failure on last attempted step) !\n");
        break;
      case -7:
      case -9:
      case -10:
      case -11:
        fprintf (fptr, "(Corrector cannot converge) !\n");
        break;
      case -8:
        fprintf (fptr, "(Singular matrix of partial derivatives) !\n");
        break;
      case -12:
        fprintf (fptr, "(DASSL failed to compute YPRIME) !\n");
        break;
      default:
        fprintf (fptr, "(Terrible internal error) !\n");
    }
      // check the time and close the output stream
    t2 = clock ();
    timer = (t2 - t1) / 1e6;
    fprintf (fptr, "RLaB: Integration lasted %g sec.\n", timer);
    fclose (fptr);
  }

  mdr_Destroy( scratch );
  mdr_Destroy( w1 );

  mdr_Destroy (mflag);
  mdr_Destroy (atol);
  mdr_Destroy (rtol);
  mdr_Destroy (rpar);
  mdr_Destroy (ipar);
  mdr_Destroy (y);
  mdr_Destroy (xi);

  // ent_t
  MDPTR(t) = 0;
  ent_DecRef (ent_t);
  ent_Destroy (ent_t);

  // ent_x
  MDPTR(x) = 0;
  ent_DecRef (ent_x);
  ent_Destroy (ent_x);

  // ent_u
  MDPTR(u) = 0;
  ent_DecRef (ent_u);
  ent_Destroy (ent_u);

  // ent_ux
  MDPTR(ux) = 0;
  ent_DecRef (ent_ux);
  ent_Destroy (ent_ux);

  // ent_uxx
  MDPTR(uxx) = 0;
  ent_DecRef (ent_uxx);
  ent_Destroy (ent_uxx);

  //
  // entities used to store constant df, dbl, dbr
  //

  ent_Clean (e2);
  ent_Clean (e7);
  ent_Clean (e9);

  //
  // entities used to store constant df, dbl, dbr
  //
  ent_Clean (pdecol_u0_name);
  ent_Clean (pdecol_f_name);
  ent_Clean (pdecol_fjac_name);
  ent_Clean (bacol_bleft_name);
  ent_Clean (bacol_bright_name);
  ent_Clean (pdecol_dbleft_name);
  ent_Clean (pdecol_dbright_name);

  return ent_Assign_Rlab_MDR(w);
}

//
// PDECOL functions
//
int
PDE1F (double * tp, double * xp, double *up, double *uxp, double *uxxp,
       double * fval, int * idummy
      )
{
  int i;
  Ent *rent = 0;
  MDR *retm = 0;

  //
  // pass the pointers to the entities being submitted as the function arguments
  //
  // t
  MDPTR(t) = (void *) tp;
  // x
  MDPTR(x) = (void *) xp;
  // u
  MDPTR(u) = (void *) up;
  // ux
  MDPTR(ux) = (void *) uxp;
  // uxx
  MDPTR(uxx) = (void *) uxxp;

  //
  // call rlab  function(t,x,u,ux,uxx)
  //
  rent = ent_call_rlab_script_5args(pdecol_f_name, ent_t, ent_x, ent_u, ent_ux, ent_uxx);

  if (ent_type (rent) != MATRIX_DENSE_REAL)
    rerror (THIS_SOLVER ": " RLAB_ERROR_RHS_FUNC_MUST_RETURN_MDR);

  retm = ent_data (rent);

  if (SIZE(retm) != pdecol_npde)
    rerror (THIS_SOLVER ": " RLAB_ERROR_RHS_FUNC_INCORRECT_DIM);

  for (i = 0; i < pdecol_npde; i++)
    fval[i] = MdrV0 (retm, i);

  ent_Clean (rent);
  return 1;
}


//
// jacobian function
//
int
PDE1FJAC (double * tp, double * xp, double *up, double *uxp, double *uxxp,
          double * dfdup, double * dfduxp, double * dfduxxp, int * idummy
          )
{
  int i, j;
  Ent *rent = 0;
  MDR *retm = 0;
  ListNode *node;

  //
  // pass the pointers to the entities being submitted as the function arguments
  //
  // t
  MDPTR(t) = (void *) tp;
  // x
  MDPTR(x) = (void *) xp;
  // u
  MDPTR(u) = (void *) up;
  // ux
  MDPTR(ux) = (void *) uxp;
  // uxx
  MDPTR(uxx) = (void *) uxxp;

  //
  // call rlab  function(t,x,u,ux,uxx)
  //
  rent = ent_call_rlab_script_5args(pdecol_fjac_name, ent_t, ent_x, ent_u, ent_ux, ent_uxx);

  if (ent_type (rent) != BTREE)
    rerror ("pdecol/bacol: jacobian must return a list <<dfdu;dfdux;dfduxx>>");

  //
  // Pick data from the list
  //
  // dfdu
  node = btree_FindNode (ent_data (rent), "dfdu");
  if (!node)
    rerror ("pdecol/bacol: missing 'dfdu' from jacobian list");
  retm = class_matrix_real (var_ent (node));
  if (!retm)
    rerror ("pdecol/bacol: jacobian entry 'dfdu' must be real matrix");
  if (SIZE(retm) != pdecol_npde * pdecol_npde)
    rerror ("pdecol/bacol: jacobian list element 'dfdu' is incorrectly dimensioned");
  for (i = 0; i < pdecol_npde; i++)
    for (j = 0; j < pdecol_npde; j++)
      dfdup[pdecol_npde * j + i] = Mdr0(retm, i, j);


  // dfdux
  node = btree_FindNode (ent_data (rent), "dfdux");
  if (!node)
    rerror ("pdecol/bacol: missing 'dfdux' from jacobian list");
  retm = class_matrix_real (var_ent (node));
  if (!retm)
    rerror ("pdecol/bacol: jacobian entry 'dfdux' must be real matrix");
  if (SIZE(retm) != pdecol_npde * pdecol_npde)
    rerror ("pdecol/bacol: jacobian list element 'dfdux' is incorrectly dimensioned");
  for (i = 0; i < pdecol_npde; i++)
    for (j = 0; j < pdecol_npde; j++)
      dfduxp[pdecol_npde * j + i] = Mdr0(retm, i, j);

  // dfduxx
  node = btree_FindNode (ent_data (rent), "dfduxx");
  if (!node)
    rerror ("pdecol/bacol: missing 'dfduxx' from jacobian list");
  retm = class_matrix_real (var_ent (node));
  if (SIZE(retm) != pdecol_npde * pdecol_npde)
    rerror ("pdecol/bacol: jacobian list element 'dfduxx' is incorrectly dimensioned");
  for (i = 0; i < pdecol_npde; i++)
    for (j = 0; j < pdecol_npde; j++)
      dfduxxp[pdecol_npde * j + i] = Mdr0(retm, i, j);

  ent_Clean (rent);
  return 1;
}


int
PDE1BND (double * tp, double *xp, double *up, double *uxp,
         double *dbdup, double * dbduxp, double *dzdtp, int * idummy
        )
{
  int i,j;
  Ent *rent = 0;
  MDR *retm = 0;
  ListNode *node;

  //
  // pass the pointers to the entities being submitted as the function arguments
  //
  // t
  MDPTR(t) = (void *) tp;
  // u
  MDPTR(u) = (void *) up;
  // ux
  MDPTR(ux) = (void *) uxp;

  //
  // call appropriate rlab  function(t,u,ux)
  //
  if (*xp == xleft)
    rent = ent_call_rlab_script_3args(pdecol_dbleft_name, ent_t, ent_u, ent_ux);
  else if (*xp == xright)
    rent = ent_call_rlab_script_3args(pdecol_dbright_name, ent_t, ent_u, ent_ux);
  else
    rerror (THIS_SOLVER ": horrible internal error left or right boundary expected");

  if (ent_type(rent) != BTREE)
    rerror (THIS_SOLVER ": boundary function must return a list <<dbdu;dbdux;dzdt>>");

  //
  // Pick data from the list
  //
  // dbdu
  node = btree_FindNode (ent_data (rent), "dbdu");
  if (!node != 0)
    rerror ("pdecol: missing 'dfdu' from boundary list");

  retm = class_matrix_real (var_ent (node));
  if (SIZE(retm) != pdecol_npde * pdecol_npde)
    rerror ("pdecol: boundary list element 'dbdu' is incorrectly dimensioned");

  for (i = 0; i < pdecol_npde; i++)
    for (j = 0; j < pdecol_npde; j++)
      dbdup[pdecol_npde * j + i] = Mdr0(retm, i, j);


  // dfdux
  node = btree_FindNode (ent_data (rent), "dbdux");
  if (!node)
    rerror ("pdecol: missing 'dfdux' from boundary list");

  retm = class_matrix_real (var_ent (node));
  if (SIZE(retm) != pdecol_npde * pdecol_npde)
    rerror ("pdecol: boundary list element 'dbdux' is incorrectly dimensioned");

  for (i = 0; i < pdecol_npde; i++)
    for (j = 0; j < pdecol_npde; j++)
      dbduxp[pdecol_npde * j + i] = Mdr0(retm, i, j);


  // dfduxx
  node = btree_FindNode (ent_data (rent), "dzdt");
  if (!node)
    rerror ("pdecol: missing 'dzdt' from boundary list");

  retm = class_matrix_real (var_ent (node));
  if (SIZE(retm) != pdecol_npde)
    rerror ("pdecol: boundary list element 'dzdt' is incorrectly dimensioned");

  for (i = 0; i < pdecol_npde; i++)
    dzdtp[i] = MdrV0(retm, i);


  ent_Clean (rent);
  return 1;
}

int
BAC1BXA (double * tp, double *up, double *uxp,
          double *bval, int * idummy
         )
{
  int i;
  Ent *rent = 0;
  MDR *retm = 0;

  //
  // pass the pointers to the entities being submitted as the function arguments
  //
  // t
  MDPTR(t) = (void *) tp;
  // u
  MDPTR(u) = (void *) up;
  // ux
  MDPTR(ux) = (void *) uxp;

  //
  // call appropriate rlab  function(t,u,ux)
  //
  rent = ent_call_rlab_script_3args(bacol_bleft_name, ent_t, ent_u, ent_ux);

  if (ent_type (rent) != MATRIX_DENSE_REAL)
    rerror (THIS_SOLVER ": " RLAB_ERROR_RHS_FUNC_MUST_RETURN_MDR);

  retm = ent_data (rent);

  if (SIZE(retm) != pdecol_npde)
    rerror (THIS_SOLVER ": " RLAB_ERROR_RHS_FUNC_INCORRECT_DIM);

  for (i = 0; i < pdecol_npde; i++)
    bval[i] = MdrV0 (retm, i);

  ent_Clean (rent);
  return 1;
}

int
BAC1BXB (double * tp, double *up, double *uxp,
         double *bval, int * idummy
        )
{
  int i;
  Ent *rent = 0;
  MDR *retm = 0;

  //
  // pass the pointers to the entities being submitted as the function arguments
  //
  // t
  MDPTR(t) = (void *) tp;
  // u
  MDPTR(u) = (void *) up;
  // ux
  MDPTR(ux) = (void *) uxp;

  //
  // call appropriate rlab  function(t,u,ux)
  //
  rent = ent_call_rlab_script_3args(bacol_bright_name, ent_t, ent_u, ent_ux);

  if (ent_type (rent) != MATRIX_DENSE_REAL)
    rerror (THIS_SOLVER ": " RLAB_ERROR_RHS_FUNC_MUST_RETURN_MDR);

  retm = ent_data (rent);

  if (SIZE(retm) != pdecol_npde)
    rerror (THIS_SOLVER ": " RLAB_ERROR_RHS_FUNC_INCORRECT_DIM);

  for (i = 0; i < pdecol_npde; i++)
    bval[i] = MdrV0 (retm, i);

  ent_Clean (rent);
  return 1;
}

int
BAC1DBXA (double * tp, double *up, double *uxp,
          double *dbdup, double * dbduxp, double *dzdtp, int * idummy
         )
{
  int i,j;
  Ent *rent = 0;
  MDR *retm = 0;
  ListNode *node;

  //
  // pass the pointers to the entities being submitted as the function arguments
  //
  // t
  MDPTR(t) = (void *) tp;
  // u
  MDPTR(u) = (void *) up;
  // ux
  MDPTR(ux) = (void *) uxp;

  //
  // call appropriate rlab  function(t,u,ux)
  //
  rent = ent_call_rlab_script_3args(pdecol_dbleft_name, ent_t, ent_u, ent_ux);

  //
  // Get the result and check whether it is a list
  //
  if (ent_type(rent) != BTREE)
    rerror (THIS_SOLVER ": boundary function must return a list <<dbdu;dbdux;dzdt>>");

  //
  // Pick data from the list
  //
  // dbdu
  node = btree_FindNode (ent_data (rent), "dbdu");
  if (!node)
    rerror ("bacol: missing 'dfdu' from boundary list");

  retm = class_matrix_real (var_ent (node));
  if (SIZE(retm) != pdecol_npde * pdecol_npde)
    rerror ("bacol: boundary list element 'dbdu' is incorrectly dimensioned");
  for (i = 0; i < pdecol_npde; i++)
    for (j = 0; j < pdecol_npde; j++)
      dbdup[pdecol_npde * j + i] = Mdr0(retm, i, j);


  // dfdux
  node = btree_FindNode (ent_data (rent), "dbdux");
  if (!node)
    rerror ("bacol: missing 'dfdux' from boundary list");

  retm = class_matrix_real (var_ent (node));
  if (SIZE(retm) != pdecol_npde * pdecol_npde)
    rerror ("bacol: boundary list element 'dbdux' is incorrectly dimensioned");
  for (i = 0; i < pdecol_npde; i++)
    for (j = 0; j < pdecol_npde; j++)
      dbduxp[pdecol_npde * j + i] = Mdr0(retm, i, j);


  // dfduxx
  node = btree_FindNode (ent_data (rent), "dbdt");
  if (!node)
    rerror ("bacol: missing 'dzdt' from boundary list");

  retm = class_matrix_real (var_ent (node));
  if (SIZE(retm) != pdecol_npde)
    rerror ("bacol: boundary list element 'dbdt' is incorrectly dimensioned");
  for (i = 0; i < pdecol_npde; i++)
    dzdtp[i] = MdrV0(retm, i);


  ent_Clean (rent);
  return 1;
}

int
BAC1DBXB (double * tp, double *up, double *uxp,
          double *dbdup, double * dbduxp, double *dzdtp, int * idummy
         )
{
  int i,j;
  Ent *rent = 0;
  MDR *retm = 0;
  ListNode *node;

  //
  // pass the pointers to the entities being submitted as the function arguments
  //
  // t
  MDPTR(t) = (void *) tp;
  // u
  MDPTR(u) = (void *) up;
  // ux
  MDPTR(ux) = (void *) uxp;

  //
  // call appropriate rlab  function(t,u,ux)
  //
  rent = ent_call_rlab_script_3args(pdecol_dbright_name, ent_t, ent_u, ent_ux);

  //
  // Get the result and check whether it is a list
  //
  if (ent_type(rent) != BTREE)
    rerror (THIS_SOLVER ": boundary function must return a list <<dbdu;dbdux;dzdt>>");

  //
  // Pick data from the list
  //
  // dbdu
  node = btree_FindNode (ent_data (rent), "dbdu");
  if (!node)
    rerror ("bacol: missing 'dfdu' from boundary list");

  retm = class_matrix_real (var_ent (node));
  if (SIZE(retm)*MNC(retm) != pdecol_npde * pdecol_npde)
    rerror (THIS_SOLVER ": boundary list element 'dbdu' is incorrectly dimensioned");
  for (i = 0; i < pdecol_npde; i++)
    for (j = 0; j < pdecol_npde; j++)
      dbdup[pdecol_npde * j + i] = Mdr0(retm, i, j);


  // dfdux
  node = btree_FindNode (ent_data (rent), "dbdux");
  if (!node)
    rerror ("bacol: missing 'dfdux' from boundary list");

  retm = class_matrix_real (var_ent (node));
  if (MNR(retm)*MNC(retm) != pdecol_npde * pdecol_npde)
    rerror ("bacol: boundary list element 'dbdux' is incorrectly dimensioned");
  for (i = 0; i < pdecol_npde; i++)
    for (j = 0; j < pdecol_npde; j++)
      dbduxp[pdecol_npde * j + i] = Mdr0(retm, i, j);


  // dfduxx
  node = btree_FindNode (ent_data (rent), "dbdt");
  if (!node)
    rerror ("bacol: missing 'dbdt' from boundary list");

  retm = class_matrix_real (var_ent (node));
  if (MNR(retm)*MNC(retm) != pdecol_npde)
    rerror ("bacol: boundary list element 'dbdt' is incorrectly dimensioned");
  for (i = 0; i < pdecol_npde; i++)
    dzdtp[i] = MdrV0(retm, i);


  ent_Clean (rent);
  return 1;
}


int
PDE1UINT (double * xp, double *up, int * idummy)
{
  int i;
  Ent *rent = 0;
  MDR *retm = 0;

  //
  // pass the pointers to the entities being submitted as the function arguments
  //
  // x
  MDPTR(x) = (void *) xp;

  //
  // call rlab  function(x)
  //
  rent = ent_call_rlab_script_1arg (pdecol_u0_name, ent_x);

  if (ent_type (rent) != MATRIX_DENSE_REAL)
    rerror (THIS_SOLVER ": " RLAB_ERROR_RHS_FUNC_MUST_RETURN_MDR);

  retm = ent_data (rent);

  if (SIZE(retm) != pdecol_npde)
    rerror (THIS_SOLVER ": " RLAB_ERROR_RHS_FUNC_INCORRECT_DIM);

  for (i = 0; i < pdecol_npde; i++)
    up[i] = MdrV0 (retm, i);

  ent_Clean (rent);
  return 1;
}

