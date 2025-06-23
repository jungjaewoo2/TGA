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
// newuoa.f
//
#include "libmjdpowell.h"

//
// gsl headers
//
// shared object
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
// conmax.f
//
static int cmx_dimg;   // actual number of constraints (dim g)
static int cmx_dimf;   // number of functions being optimized (dim f)
static int cmx_numgr;  // total number of constraints (dim g + dim f)
static int cmx_dimx;   // dimension of the vector x
static int cmx_indfn;  // did user specify jacobians
static int cmx_typset;
static int cmx_islin;  // is the constraint an affine,   a.x in b
static int pb_isfun;
static MDR *cmx_icntyp;
static Ent *cmx_f_name, *cmx_g_name, *cmx_df_name, *cmx_dg_name;

//
// the GSL
//
// // scalar function of a vector variable

//
// rlab/c functions for conmax solver
//
static int CMXFNSET (int *, int *, double *, int *, int *, double *, int *,
                     int *, int *, double *);


//
// global arguments for rlab-scripted functions
//
static MDR *xmdr=0;
static Ent *xent=0;
static Ent *pent=0;
static MDR *cmx_A=0, *cmx_B=0, *pb_F=0;

#undef  THIS_SOLVER
#define THIS_SOLVER "findmins"
#include "gsl_findmins_utils.c"
Ent * ent_gsl_findmins (int nargs, Datum args[])
{
  Ent *e2=0, *e3=0, *eo=0;
  ListNode *node;

  double ddummy, abserr = 1e-4, target = create_nan();
  int idummy, maxiter=1000, istand=100000;

  char *outs=0;

  MDR * w=0, *ystart=0;
  FILE *rlab_stderr = (!RLAB_STDERR_DS) ? stderr : RLAB_STDERR_DS;
  pent = 0;

  if (nargs < 3 || nargs > 5)
  {
    fprintf
        (rlab_stderr, THIS_SOLVER ": Unconstrained minimization of a function in multidimension.\n");
    fprintf
        (rlab_stderr, THIS_SOLVER ": Format:\n");
    fprintf
        (rlab_stderr, THIS_SOLVER ": (1) x = findmins(f,/p/,x0 /,options/),\n");
    fprintf
        (rlab_stderr, THIS_SOLVER ": (2) x = findmins(f,df,/p/,x0 /,options/),\n");
    fprintf
        (rlab_stderr, THIS_SOLVER ": where f=function(x/,p/) is the function being minimized,\n");
    fprintf
        (rlab_stderr, THIS_SOLVER ": 'p' is its parameter array, df=function(x,p) is its\n");
    fprintf
        (rlab_stderr, THIS_SOLVER ": jacobian df/dx, 'x0' is the starting point. In the case (1)\n");
    fprintf
        (rlab_stderr, THIS_SOLVER ": options=<<eabs;ss;maxi;rlab_stderr>>. In the case (2) when \n");
    fprintf
        (rlab_stderr, THIS_SOLVER ": imethod=0..3 (GSL ) options=<<eabs;tol;h;rlab_stderr;maxi;imethod>>,\n");
    fprintf
        (rlab_stderr, THIS_SOLVER ": while if imethod=4 (CONMAX) then options=\n");
    fprintf
        (rlab_stderr, THIS_SOLVER ": <<rlab_stderr;nrkstep;tolcon;maxi;istep;encsm;limsm;imethod>>.\n");
    fprintf
        (rlab_stderr, THIS_SOLVER ": See manual for more information.\n");
    fprintf
        (rlab_stderr, THIS_SOLVER ": The function returns the best point 'x' found using criteria.\n");
    rerror
        ("requires at least 3 arguments !");
  }

  //
  // Get 1st parameter:  a function ptr
  //
  cmx_f_name = bltin_get_ent(args[0]);
  if (!isfuncent(cmx_f_name))
    rerror (THIS_SOLVER ": " RLAB_ERROR_ARG1_FUNC_VAR "\n");

  //
  // Get 2nd parameter and check its type: MDR or function
  //
  e2 = bltin_get_ent( args[1] );
  if (!isfuncent(e2))
  {
    //
    // second argument is not a function, use non-jacobian methods
    //
    MDR *h=0;
    int imethod = 0;

    //newuoa
    MDR *x_bounds=0, *x_constr=0;
    double rhobeg = 10;
    double rhoend = 0.1;
    int    npt = 0;

    //
    // parameter entity: if not UNDEF make a copy of it for user function
    //
    if (ent_type(e2) != UNDEF)
      pent = ent_Copy(e2);

    //
    // y0
    //
    if (nargs < 3)
      rerror (THIS_SOLVER ": " RLAB_ERROR_AT_LEAST_THREE_ARG_REQUIRED);
    e3 = bltin_get_ent (args[2]);
    if (not_ent_double_vector(e3))
      rerror(THIS_SOLVER ": " RLAB_ERROR_ARG3_MDR_VECTOR);
    ystart = ent_data(e3);
    cmx_dimx = SIZE (ystart);

    if (nargs > 3)
    {
      eo = bltin_get_ent (args[3]);
      if (ent_type (eo) == BTREE)
      {
        // imethod can be given explicitely
        node = btree_FindNode (ent_data (eo), RLAB_NAME_GEN_IMETHOD);
        if (node != 0)
        {
          idummy = (int) class_double (var_ent (node));
          if (idummy>=0 && idummy<=2)
            imethod = idummy;
        }

        //
        // parameters for NEWUOA or BYBUOA solver
        //
        // bounds
        node = btree_FindNode (ent_data (eo), RLAB_NAME_POW_BOUNDS);
        if (node != 0)
        {
          if (ent_type(var_ent(node)) == MATRIX_DENSE_REAL)
          {
            x_bounds = class_matrix_real( var_ent(node) );
            imethod = 1;
          }
        }
        // constraint
        node = btree_FindNode (ent_data (eo), RLAB_NAME_POW_CONSTR);
        if (node != 0)
        {
          if (ent_type(var_ent(node)) == MATRIX_DENSE_REAL)
          {
            x_constr = class_matrix_real( var_ent(node) );
            imethod = 1;
          }
          else if (isfuncent(var_ent(node)))
          {
            cmx_g_name = var_ent(node);
            cmx_dimg = -1;  // at this time we don't know what is the dimension of constraints
            imethod = 3;
          }
        }
        // rhobeg
        node = btree_FindNode (ent_data (eo), RLAB_NAME_NUO_RHOBEG);
        if (node != 0)
        {
          rhobeg = class_double (var_ent (node));
          if (imethod < 1)
            imethod = 1;
        }
        // rhoend
        node = btree_FindNode (ent_data (eo), RLAB_NAME_NUO_RHOEND);
        if (node != 0)
        {
          rhoend = class_double (var_ent (node));
          if (imethod < 1)
            imethod = 1;
        }
        if (rhoend > rhobeg)
        {
          fprintf(rlab_stderr, THIS_SOLVER ": RHOEND cannot be greater than RHOBEG!\n");
          double dummy = rhoend;
          rhoend = rhobeg;
          rhobeg = dummy;
        }
        // ntp: initialize default value:
        npt = (double) 0.5*(cmx_dimx+2 + 0.5*(cmx_dimx+1)*(cmx_dimx+2));
        node = btree_FindNode (ent_data (eo), RLAB_NAME_NUO_NTP);
        if (node != 0)
        {
          npt = (int) class_double (var_ent (node));
          if (imethod < 1)
            imethod = 1;
        }
        if (npt < cmx_dimx+2 || npt > 0.5*(cmx_dimx+1)*(cmx_dimx+2))
        {
          npt = (double) 0.5*(cmx_dimx+2 + 0.5*(cmx_dimx+1)*(cmx_dimx+2));
          fprintf(rlab_stderr, THIS_SOLVER ": NPT out of bounds. Reverting to default value %i.\n", npt);
        }

        // target
        node = btree_FindNode (ent_data (eo), RLAB_NAME_GSL_MINS_TARGET);
        if (node != 0)
          target = class_double (var_ent (node));

        // eabs
        node = btree_FindNode (ent_data (eo), RLAB_NAME_GSL_MINS_EABS);
        if (node != 0)
        {
          ddummy = class_double (var_ent (node));
          if (ddummy > 0.0)
            abserr = ddummy;
        }
        // ss
        node = btree_FindNode (ent_data (eo), RLAB_NAME_GSL_MINS_SS);
        if (node != 0)
        {
          if (ent_type(var_ent(node)) == MATRIX_DENSE_REAL)
          {
            h = class_matrix_real ( var_ent(node) );
          }
        }
         // max iterations at the same value of the objective function
        node = btree_FindNode (ent_data (eo), RLAB_NAME_GSL_MINS_STANDSTILL);
        if (node != 0)
        {
          idummy = (int) class_double (var_ent (node));
          if (idummy > 0)
            istand = idummy;
        }

        //
        // parameters that appear for both solvers
        //
        // standard output
        node = btree_FindNode (ent_data (eo), RLAB_NAME_GSL_MINS_STDOUT);
        if (node != 0)
        {
          if (ent_type(var_ent (node))== MATRIX_DENSE_STRING)
          {
            outs  = class_char_pointer (var_ent (node));
          }
        }

        // max iterations
        node = btree_FindNode (ent_data (eo), RLAB_NAME_GSL_MINS_MAXITER);
        if (node != 0)
        {
          idummy = (int) class_double (var_ent (node));
          if (idummy > 0)
            maxiter = idummy;
        }
      }
    }

    //
    // Set up ENTITIES for user-function.
    // Inc the reference count once, for belonging to multiroot
    //
    // array y[]:
    //
    xmdr = mdr_CreateEmpty (MNR(ystart), MNC(ystart));
    xent = ent_Assign_Rlab_MDR (xmdr);
    ent_IncRef (xent);

    //
    // parameter entity
    //
    if ((imethod == 1)||(imethod == 2))
    {
      w = mdr_findmins_powell(ystart, x_constr, x_bounds, outs, npt, rhobeg, rhoend, maxiter);
    }
    else if (imethod == 3)
    {
      w = mdr_findmins_powell_cobyla(ystart, outs, rhobeg, rhoend, maxiter);
    }
    else
    {
      w = mdr_findmins_simplex(ystart, h, outs, abserr, target, maxiter, istand);
    }

    // Clean Up
    ent_Clean( e2 );
    ent_Clean( pent );

    MDPTR(xmdr) = 0;
    ent_DecRef (xent);
    ent_Destroy (xent);

    ent_Clean (e3);
    ent_Clean (eo);
    ent_Clean (cmx_f_name);
    ent_Clean (cmx_g_name);

    return ent_Assign_Rlab_MDR (w);
  }

  //
  // second parameter is a function, assume dfunc and use appropriate solvers
  //
  if (nargs < 4)
    rerror (THIS_SOLVER ": " RLAB_ERROR_AT_LEAST_FOUR_ARG_REQUIRED);

  Ent *edummy = 0;

  //
  // CONMAX parameters
  //
  int nrkstep = 1, limsm = 0, istep = 0;
  double tolcon = 0.0, encsm = 0.0;

  //
  // proximal bundle (PB) parameters
  //
  int pb_met = 1, pb_mes = 1, pb_mtesx = 20;
  int pb_mtesf = 7;
  double pb_tolx = 1e-16, pb_tolf = 1e-8, pb_tolb = -1e60;
  double pb_tolg = 1e-6, pb_told = 1e-4, pb_tols = 1e-2;
  double pb_tolp = 0.5, pb_eta = 0.1, pb_xmax = 1e3;
//   double pb_tol[6] = { pb_tolx, pb_tolf, pb_tolg, pb_told, pb_tols, pb_tolp };

  //
  // the GSL
  //
  double tol = 1e-4, h = 0.01;
  int imethod = 0;

  cmx_df_name = e2;

  //
  // parameter entity: if not UNDEF just pass it to user function
  //
  e2 = bltin_get_ent( args[2] );
  if (ent_type(e2) != UNDEF)
    pent = ent_Copy(e2);

  //
  // y0
  //
  e3 = bltin_get_ent (args[3]);
  if (not_ent_double_vector(e3))
    rerror(THIS_SOLVER ": " RLAB_ERROR_ARG3_MDR_VECTOR);
  ystart = ent_data(e3);
  cmx_dimx = SIZE(ystart);

  //
  // options argument
  //
  if (nargs > 4)
  {
    eo = bltin_get_ent (args[4]);
    if (ent_type (eo) == BTREE)
    {
      // imethod
      node = btree_FindNode (ent_data (eo), RLAB_NAME_GSL_MINS_IMETHOD);
      if (node != 0)
      {
        idummy = (int) class_double (var_ent (node));
        if (idummy >= 0 && idummy <= 6)
          imethod = idummy;
      }
      // gsl: eabs
      node = btree_FindNode (ent_data (eo), RLAB_NAME_GSL_MINS_EABS);
      if (node != 0)
      {
        ddummy = class_double (var_ent (node));
        if (ddummy > 0.0)
          abserr = ddummy;
      }
      // pb: tol
      node = btree_FindNode (ent_data (eo), RLAB_NAME_PB_MINS_TOLGDSP);
      if (node != 0)
      {
        if (imethod == 5 || imethod == 6)
        {
          if (ent_type(var_ent(node)) == MATRIX_DENSE_REAL)
          {
            MDR * rdummy = class_matrix_real(var_ent (node));
            if ( MNR (rdummy) * MNC (rdummy) >= 4)
            {
              pb_tolg = MdrV0 (rdummy,0);
              pb_told = MdrV0 (rdummy,1);
              pb_tols = MdrV0 (rdummy,2);
              pb_tolp = MdrV0 (rdummy,3);
            }
          }
        }
        else
        {
          ddummy = class_double (var_ent (node));
          if (ddummy > 0.0)
            tol = ddummy;
        }
      }
      // gsl: h
      node = btree_FindNode (ent_data (eo), RLAB_NAME_GSL_MINS_STEP);
      if (node != 0)
      {
        ddummy = class_double (var_ent (node));
        if (ddummy > 0.0)
          h = ddummy;
      }
      // gsl: standard output
      node = btree_FindNode (ent_data (eo), RLAB_NAME_GSL_MINS_STDOUT);
      if (node != 0)
      {
        if (ent_type(var_ent (node))== MATRIX_DENSE_STRING)
        {
          outs = class_char_pointer (var_ent (node));
        }
      }
      // max iterations
      node = btree_FindNode (ent_data (eo), RLAB_NAME_GSL_MINS_MAXITER);
      if (node != 0)
      {
        idummy = (int) class_double (var_ent (node));
        if (idummy > 0)
          maxiter = idummy;
      }
      // conmax: nrkstep
      node = btree_FindNode (ent_data (eo), RLAB_NAME_CNM_NRKSTEP);
      if (node != 0)
      {
        idummy = (int) class_double (var_ent (node));
        if (idummy > 0)
          nrkstep = idummy;
      }
      // conmax: tolcon
      node = btree_FindNode (ent_data (eo), RLAB_NAME_CNM_TCON);
      if (node != 0)
      {
        ddummy = class_double (var_ent (node));
        if (ddummy > 0.0)
          tolcon = ddummy;
      }
      // conmax: istep
      node = btree_FindNode (ent_data (eo), RLAB_NAME_CNM_ISTEP);
      if (node != 0)
      {
        idummy = (int) class_double (var_ent (node));
        if (idummy >= 0 && idummy <= 2)
          istep = idummy;
      }
      // pb: eta
      node = btree_FindNode (ent_data (eo), RLAB_NAME_CNM_ETA);
      if (node != 0)
      {
        ddummy = class_double (var_ent (node));
        if (ddummy > 0.0)
          pb_eta = ddummy;
      }
      // pb: met
      node = btree_FindNode (ent_data (eo), RLAB_NAME_CNM_MET);
      if (node != 0)
      {
        idummy = (int) class_double (var_ent (node));
        if (idummy >= 1 && idummy <= 3)
          pb_met = idummy;
      }
      // pb: mes
      node = btree_FindNode (ent_data (eo), RLAB_NAME_CNM_MES);
      if (node != 0)
      {
        idummy = (int) class_double (var_ent (node));
        if (idummy >= 1 && idummy <= 4)
          pb_mes = idummy;
      }
      // pb: convx
      node = btree_FindNode (ent_data (eo), RLAB_NAME_CNM_CONVX);
      if (node != 0)
      {
        if (ent_type(var_ent(node)) == MATRIX_DENSE_REAL)
        {
          MDR * rdummy = class_matrix_real(var_ent (node));
          if ( MNR (rdummy) * MNC (rdummy) >= 2)
          {
              // pb: [tolx, mtesx]
            pb_tolx  = MdrV0(rdummy,0);
            pb_mtesx = (int) MdrV0(rdummy,1);
          }
        }
      }
      // conmax,pb: convf
      node = btree_FindNode (ent_data (eo), RLAB_NAME_CNM_CONVF);
      if (node != 0)
      {
        if (ent_type(var_ent(node)) == MATRIX_DENSE_REAL)
        {
          MDR * rdummy = class_matrix_real(var_ent (node));
          if ( MNR (rdummy) * MNC (rdummy) >= 2)
          {
            if (imethod == 5 || imethod == 6)
            {
              // pb: [tolf, mtesf]
              pb_tolf  = MdrV0(rdummy,0);
              pb_mtesf = (int) MdrV0(rdummy,1);
            }
            else if (imethod == 4)
            {
              // conmax: [encsm, limsm]
              encsm = MdrV0(rdummy,0);
              limsm = (int) MdrV0(rdummy,1);
            }
          }
        }
      }
    }
  }

  //
  // arguments for the rlab functions
  //
  xmdr = mdr_CreateEmpty (MNR(ystart), MNC(ystart));
  xent = ent_Assign_Rlab_MDR (xmdr);
  ent_IncRef (xent);

  // ===============================================
  //
  // CALL TO OPTIMIZATION CODES
  //
  // ===============================================
  switch(imethod)
  {
    case 0:
    case 1:
    case 2:
    case 3:
      w = mdr_findmins_f_df_gsl(imethod, ystart, h, outs, abserr, tol, target, maxiter);
      break;

    case 4:
      // global variables
      cmx_typset = 0;
      cmx_indfn = 1;
      // find dimension
      if (pent)
        edummy = ent_call_rlab_script_2args(cmx_f_name, e3, pent);
      else
        edummy = ent_call_rlab_script_1arg (cmx_f_name, e3);
      if (ent_type (edummy) != MATRIX_DENSE_REAL)
        rerror (THIS_SOLVER ": " RLAB_ERROR_RHS_FUNC_MUST_RETURN_MDR);
      MDR * retf = ent_data (edummy);
      cmx_dimf = SIZE(retf);
      ent_Clean (edummy);
      // other dimensions
      cmx_dimg  = 0;
      cmx_numgr = cmx_dimf + cmx_dimg;

      w = mdr_findmins_f_df_conmax(cmx_dimf, ystart, outs, encsm, limsm, nrkstep, tolcon, istep, maxiter);
      break;

    case 5:
    case 6:
      w = mdr_findmins_f_df_pbun(imethod, cmx_dimf, ystart, outs,
                                 pb_met, pb_mes, pb_mtesx, pb_mtesf, maxiter,
                                 pb_tolx, pb_tolf, pb_tolb,
                                 pb_tolg, pb_told, pb_tols,
                                 pb_tolp, pb_eta,  pb_xmax
                                );
  }

  // Clean Up
  ent_Clean (pent);
  pent = 0;
  ent_Clean (e2);

  MDPTR(xmdr) = 0;
  ent_DecRef (xent);
  ent_Destroy (xent);

  ent_Clean (e3);
  ent_Clean (cmx_f_name);
  ent_Clean (cmx_df_name);

  ent_Clean (eo);

  return ent_Assign_Rlab_MDR (w);
}



#undef  THIS_SOLVER
#define THIS_SOLVER "conmins"
Ent * ent_conmax (int nargs, Datum args[])
{
  Ent *ea=0, *eb=0, *edim=0, *e5=0, *e6=0, *e7=0, *eo=0, *rent;
  ListNode *node;

  pent = 0;

  //
  // conmax variables
  //
  int ioptn = 0, iptb, indb, liwrk, lwrk, iter, nparm;
  int i, j, k, m, idummy, louts=0;
  int nrkstep = 1, maxiter=1000, limsm = 0, istep = 0;
  double tolcon = 0.0, encsm = 0.0;

  //
  // ProxyBundle variables
  //
  int pb_met = 1, pb_mes = 1, pb_mtesx = 20;
  int pb_mtesf = 7;
  double pb_tolx = 1e-16, pb_tolf = 1e-8, pb_tolb = -1e60;
  double pb_tolg = 1e-6, pb_told = 1e-4, pb_tols = 1e-2;
  double pb_tolp = 0.5, pb_eta = 0.1, pb_xmax = 1e3;
//   double pb_tol[6] = { pb_tolx, pb_tolf, pb_tolg, pb_told, pb_tols, pb_tolp };

  double ddummy, timer;
  int imethod = 0;

  MDR *x7=0, *fun=0, *iwork=0, *work=0, *w=0, *err=0, *xstart=0;
  MDR *pb_cg=0, *pb_cu=0, *pb_cl=0;

  cmx_typset = 0;
  cmx_A = 0;
  cmx_B = 0;


  //
  // output
  //
  char *outs=0;
  FILE *fptr = NULL;

  time_t t1=clock(), t2=clock();
  //
  // Load and Check arguments.
  //
  if (nargs < 6 || nargs > 8)
  {
    fprintf (stdout,
             "conmins: Solves a constrained optimization problem in multidimensions,\n");
    fprintf (stdout,
             "conmins: find  x  such that  (1) min f(x) with  g(x)<=0 , or \n");
    fprintf (stdout,
             "conmins:                     (2) min|f(x) - F|  with  g(x)<=0 , or\n");
    fprintf (stdout,
             "conmins: combination of the two, where F(i) = inf() for (1), or\n");
    fprintf (stdout,
             "conmins: -inf() < F(i) < inf() for (2), with 1 <= i <= dim f.\n");
    fprintf (stdout,
             "conmins: Format:\n");
    fprintf (stdout,
             "conmins: (1)  x = conmins(f,df,g,dg,/p/,x0,/F/,/options/),\n");
    fprintf (stdout,
             "conmins: (2)  x = conmins(f,df,a, b,/p/,x0,/F/,/options/),\n");
    fprintf (stdout,
             "conmins: where  f=function(x), df=function(x) is its jacobian, \n");
    fprintf (stdout,
             "conmins: and similarly for the constraint condition 'g' and 'dg',\n");
    fprintf (stdout,
             "conmins: for g=g(x). Affine constraints are possible as well, where\n");
    fprintf (stdout,
             "conmins: b[i;1] <= (a.x)[i] <= b[i;2], i=1.. dim g, with 'a' and 'b'\n");
    fprintf (stdout,
             "conmins: real matrices. The starting point is 'x0', and 'options' is\n");
    fprintf (stdout,
             "conmins: a list <<stdout;nrkstep;tolcon;maxi;istep;encsm;limsm>>.\n");
    fprintf (stdout,
             "conmins: See the manual for description. The result 'x' is the best solution.\n");
    rerror ("requires at least 8 arguments");
  }

  //
  // f for optimization
  //
  cmx_f_name = bltin_get_ent(args[0]);
  if (!isfuncent(cmx_f_name))
    rerror (THIS_SOLVER ": " RLAB_ERROR_ARG1_FUNC_VAR "\n");

  //
  // its jacobian df = function(x)
  //
  cmx_df_name = bltin_get_ent(args[1]);
  if (!isfuncent(cmx_df_name))
    rerror (THIS_SOLVER ": " RLAB_ERROR_ARG2_FUNC_VAR "\n");

  cmx_indfn = 1;

  //
  // parameter for functions
  //
  e5 = bltin_get_ent (args[4]);
  if (ent_type(e5) != UNDEF)
    pent = ent_Copy(e5);

  //
  // x0
  //
  e6 = bltin_get_ent (args[5]);
  if (ent_type (e6) == MATRIX_DENSE_REAL)
    xstart = class_matrix_real (e6);
  else
    rerror ("conmins: 'x0' has to be a vector");
  // Extract cmx_dimx from ystart
  if (MNC (xstart) != 1 && MNR (xstart) != 1)
    rerror ("conmins: 'x0' has to be a vector");
  cmx_dimx = MNR (xstart) * MNC (xstart);


  //
  // determine dimf by calling f(x,p) once and measure the size of result array
  //
  //
  if (pent)
    edim = ent_call_rlab_script_2args(cmx_f_name, e6, pent);
  else
    edim = ent_call_rlab_script_1arg (cmx_f_name, e6);
  if (ent_type (edim) != MATRIX_DENSE_REAL)
    rerror (THIS_SOLVER ": " RLAB_ERROR_RHS_FUNC_MUST_RETURN_MDR);
  MDR * retf = ent_data (edim);
  cmx_dimf = SIZE(retf);
  ent_Clean (edim);

  //
  // g for constraints
  //
  cmx_dimg = 0;
  cmx_islin = 0;

  //
  cmx_g_name = bltin_get_ent(args[2]);
  if (isfuncent(cmx_g_name))
  {
    cmx_islin = 0;

    //
    // dimg
    //
    if (pent)
      edim = ent_call_rlab_script_2args(cmx_g_name, e6, pent);
    else
      edim = ent_call_rlab_script_1arg (cmx_g_name, e6);
    if (ent_type (edim) != MATRIX_DENSE_REAL)
      rerror (THIS_SOLVER ": " RLAB_ERROR_RHS_FUNC_MUST_RETURN_MDR);
    MDR * retg = ent_data (edim);
    cmx_dimg = SIZE(retg);
    ent_Clean (edim);

    //
    // its jacobian df = function(x)
    //
    if (cmx_indfn)
    {
      cmx_dg_name = bltin_get_ent(args[3]);
      cmx_indfn = 0;
      if (isfuncent(cmx_dg_name))
        cmx_indfn = 1;
    }
  }
  else
  {
    //
    // is it a matrix 'a', part of an affine constraint   a.x in b
    //
    ea = cmx_g_name;
    eb = bltin_get_ent ( args[3] );
    if (ent_type(ea) == MATRIX_DENSE_REAL && ent_type(eb) == MATRIX_DENSE_REAL)
    {
      MDR * a = class_matrix_real (ea);
      MDR * b = class_matrix_real (eb);
      if (MNC(b) != 2)
        rerror ("conmins: constraint bound is a two-column real matrix!");
      if (MNR (b) != MNR (a))
        rerror ("conmins: constraint and its bound need same number of rows!");
      for (i = 0; i < 2; i++)
      {
        if (Mdr0(b,i,0) > Mdr0(b,i,1))
        {
          fprintf(stdout, "conmins: inconsistent bound %i in constraint bound!\n", i);
          rerror ("conmins: error in constraint bound!");
        }
        for (j = 0; j < MNR (b); j++)
          if (Mdr0(b,i,j) !=  create_inf () &&
              Mdr0(b,i,j) != -create_inf ()) cmx_dimg ++;
      }
      //
      // copy the matrix for conmax and for pb
      //
      if (cmx_dimg > 0)
      {
        pb_cg = mdr_Copy (a);
        pb_cl = mdr_PartitionCol (b, 1);
        pb_cu = mdr_PartitionCol (b, 2);
        cmx_islin = 1;
      }
    }
    else
      cmx_dimg = 0; // no constraints

    // clean-up if necessary
    ent_Clean (ea);
    ent_Clean (eb);
  }



  cmx_numgr = cmx_dimf + cmx_dimg;

  cmx_icntyp = mdi_Create(cmx_numgr, 1);
  for (i=0; i < cmx_dimf; i++)
    MdiV0(cmx_icntyp,i) = 1;                // default optimization goal: f(x) <= w
  if (cmx_islin)
    for (i=0; i < cmx_dimg; i++)
      MdiV0(cmx_icntyp,i + cmx_dimf) = -1;  // the constraints are linear (affine, a.x in b)
  else
    for (i=0; i < cmx_dimg; i++)
      MdiV0(cmx_icntyp,i + cmx_dimf) = -2;  // default assumptions: constraints are nonlinear

  // copy the initial point as column vector
  w = mdr_Float_BF(xstart);
  w->nrow = (xstart->nrow) * (xstart->ncol);
  w->ncol = 1;

  //
  // F: optimization goal
  //
  fun = mdr_Create(cmx_dimf, 1);
  mdr_Zero (fun);
  pb_isfun = 0;

  e7 = bltin_get_ent (args[6]);
  if (ent_type (e7) == MATRIX_DENSE_REAL)
  {
    x7 = class_matrix_real (e7);
    if (MNR(x7) * MNC (x7) != cmx_dimf)
      rerror("conmins: improperly dimensioned optimization goal 'F'");
    pb_isfun = 1;
    for (i=0; i < cmx_dimf; i++)
    {
      if (MdrV0(x7,i) != create_inf ())
      {
        // if x7 is given and is not inf() than optimization goal
        // becomes  |f(x) - fun| <= w
        MdiV0(cmx_icntyp,i) = 2;
        MdrV0(fun,i) = MdrV0(x7,i);
      }
      else
        pb_isfun = 0;
    }
  }
  ent_Clean (e7);

  //
  // check whether the last argument is the options list
  //
  if (nargs == 8 || nargs == 7)
  {
    eo = bltin_get_ent (args[nargs-1]);
    if (ent_type (eo) == BTREE)
    {
      // all: imethod
      node = btree_FindNode (ent_data (eo), RLAB_NAME_CNM_IMETHOD);
      if (node != 0)
      {
        idummy = (int) class_double (var_ent (node));
        if ( (cmx_islin || cmx_dimg == 0)
              && (idummy == 0 || idummy == 1 || idummy == 2 || idummy == 3))
          imethod = idummy;
        else
          imethod = 0;
      }
      // conmax,pb: convf
      node = btree_FindNode (ent_data (eo), RLAB_NAME_CNM_CONVF);
      if (node != 0)
      {
        if (ent_type(var_ent(node)) == MATRIX_DENSE_REAL)
        {
          MDR * rdummy = class_matrix_real(var_ent (node));
          if ( MNR (rdummy) * MNC (rdummy) >= 2)
          {
            if (imethod == 5 || imethod == 6)
            {
              // pb: [tolf, mtesf]
              pb_tolf  = MdrV0(rdummy,0);
              pb_mtesf = (int) MdrV0(rdummy,1);
            }
            else if (imethod == 4)
            {
              // conmax: [encsm, limsm]
              encsm = MdrV0(rdummy,0);
              limsm = (int) MdrV0(rdummy,1);
            }
          }
        }
      }
      // pb: convx
      node = btree_FindNode (ent_data (eo), RLAB_NAME_CNM_CONVX);
      if (node != 0)
      {
        if (ent_type(var_ent(node)) == MATRIX_DENSE_REAL)
        {
          MDR * rdummy = class_matrix_real(var_ent (node));
          if ( MNR (rdummy) * MNC (rdummy) >= 2)
          {
              // pb: [tolx, mtesx]
            pb_tolx  = MdrV0(rdummy,0);
            pb_mtesx = (int) MdrV0(rdummy,1);
          }
        }
      }
      // conmax: nrkstep
      node = btree_FindNode (ent_data (eo), RLAB_NAME_CNM_NRKSTEP);
      if (node != 0)
      {
        idummy = (int) class_double (var_ent (node));
        if (idummy > 0)
          nrkstep = idummy;
      }
      // conmax: tolcon
      node = btree_FindNode (ent_data (eo), RLAB_NAME_CNM_TCON);
      if (node != 0)
      {
        ddummy = class_double (var_ent (node));
        if (ddummy > 0.0)
          tolcon = ddummy;
      }
      // all: maxi
      node = btree_FindNode (ent_data (eo), RLAB_NAME_CNM_MAXITER);
      if (node != 0)
      {
        idummy = (int) class_double (var_ent (node));
        if (idummy > 0)
          maxiter = idummy;
      }
      // conmax: istep
      node = btree_FindNode (ent_data (eo), RLAB_NAME_CNM_ISTEP);
      if (node != 0)
      {
        idummy = (int) class_double (var_ent (node));
        if (idummy >= 0 && idummy <= 2)
          istep = idummy;
      }
      // all: standard output
      node = btree_FindNode (ent_data (eo), RLAB_NAME_CNM_STDOUT);
      if (node != 0)
      {
        if (ent_type(var_ent (node))== MATRIX_DENSE_STRING)
        {
          outs = class_char_pointer (var_ent (node));
          if (outs)
            louts = strlen (outs);
        }
      }
    }
  }

  //
  // Set up ENTITIES for user-function.
  // Inc the reference count once, for belonging to conmins
  //
  // x:
  //
  xmdr = mdr_CreateEmpty (MNR(xstart), MNC(xstart));
  xent = ent_Assign_Rlab_MDR (xmdr);
  ent_IncRef (xent);

  if (imethod == 0 || (!cmx_islin&& cmx_dimg>0))
  {
    //
    // conmax: (un)constrained linear and non-linear optimization
    //

    //
    // conmax: write  a.x in b   as   a'.x + b' <= 0
    //
    if (cmx_dimg > 0 && cmx_islin)
    {
      cmx_A = mdr_Create(cmx_dimg, MNC(pb_cg));
      cmx_B = mdr_Create(cmx_dimg, 1);
      k = 0;
      for (i = 0; i < MNR (pb_cl); i++)
      {
        if (MdrV0(pb_cl,i) != -create_inf ())
        {
          for (m = 0; m < MNC (pb_cg); m++)
            Mdr0(cmx_A,k,m) = (-1.0) * Mdr0(pb_cg,i,m);
          MdrV0(cmx_B,k) = MdrV0(pb_cl,i);
          k ++;
        }
        if (MdrV0(pb_cu,i) !=  create_inf ())
        {
          for (m = 0; m < MNC (pb_cg); m++)
            Mdr0(cmx_A,k,m) =  Mdr0(pb_cg,i,m);
          MdrV0(cmx_B,k) = (-1.0) * MdrV0(pb_cu,i);
          k++;
        }
      }
    }

    nparm = cmx_dimx;

    iptb = 0;
    indb = 0;

    liwrk = 7 * cmx_numgr + 7 * nparm + 3;
    iwork = mdi_Create(liwrk,1);

    lwrk = 2 * nparm * nparm + 4 * cmx_numgr * nparm + 11 * cmx_numgr +
        27 * nparm + 13;
    work  = mdr_Create(lwrk,1);

    err = mdr_Create(cmx_numgr + 3, 1);

    //
    // the other arrays needed by the code
    //
    // set jacobians if user specified them
    if (cmx_indfn)
      ioptn = 0;
    else
      ioptn = 1;
    ioptn = ioptn + 10000; // call  CMXFNSET  once per each 'x'
    // 10
    if (limsm > 0 && encsm >0)
    {
      ioptn = ioptn + 10;
      MdrV0 (work, 0) = encsm;
      MdiV0 (iwork,0) = limsm;
    }
    // 100
    if (nrkstep > 1)
    {
      MdiV0 (iwork,1) = nrkstep;
      ioptn = ioptn + 100;
    }
    if (tolcon > 0)
    {
      MdrV0 (work, 1) = tolcon;
      ioptn = ioptn + 200;
    }
    // 1000
    ioptn = ioptn + 1000 * istep;

    if (louts > 1)
    {
      fptr = fopen (outs, "a");
      if (fptr)
      {
        t1 = clock ();
        fprintf (fptr,
                 "RLaB: using solver CONMAX for constrained/unconstrained"
                     " optimization in multidimensions.\n"
                );
        fprintf (fptr,
                 "RLaB: Messages from the solver follow.\n"
                );
      }
    }

    //
    // call CONMAX now
    //
    CONMAX2(&ioptn, &nparm, &cmx_numgr, &maxiter, MDRPTR(fun), &cmx_dimf,
             NULL, &iptb, &indb,
             MDIPTR(iwork), &liwrk, MDRPTR(work), &lwrk,
                    &iter, MDRPTR(w), MDRPTR(err),
                                  outs, &louts,CMXFNSET);

    if (fptr)
    {
      fprintf (fptr, "RLaB: solver CONMAX finished the optimization search.\n");
      t2 = clock ();
      timer = (t2 - t1) / 1e6;
      fprintf (fptr, "RLaB: optimization lasted %g sec.\n", timer);
      fclose (fptr);
    }

    // clean-up allocated memory arrays
    mdr_Destroy (err);
    mdr_Destroy (work);
    mdr_Destroy (iwork);

    if (cmx_dimg > 0)
    {
      if (cmx_A) mdr_Destroy (cmx_A);
      if (cmx_B) mdr_Destroy (cmx_B);
    }

  }
  else if ((imethod == 1 || imethod == 2) && cmx_dimf==1 )
  {
    //
    // ProxyBundle: dimf == 1
    //
    int    nb = 0, nc, nf = cmx_dimx, na;
    int    nia, nra, ipar[7], iterm , ihes=0;
    double fp = 0, gmax;
    double rpar[9];

    MDR * ic=0, * ia=0, * ra=0, * cf=0, *cg=0;

    if (cmx_dimg)
    {
      cmx_dimg = MNR(pb_cg);
      cg = mdr_Transpose (pb_cg);
    }

    nf = cmx_dimx;
    na = nf + 3;
    nc = cmx_dimg;

    nia = nf + na + 1;
    ia = mdi_Create(nia, 1);

    if (imethod == 1)
    {
      // pbunl
      nra = nf*(nf+1)/2+nf*(na+5)+5*na+4;
      ipar[1 -1] = pb_met;
    }
    else
    {
      // pnewl
      nra = nf*(nf+1)*(na+3)/2+nf*(na+6)+5*na+4;
      ipar[1 -1] = 1;   // mos - distance measure exponent, 1 or 2
    }
    ra  = mdr_Create(nra,1);

    //
    // Form linear constraint matrices
    //
    if (nc > 0)
    {
      ic = mdi_Create (nc, 1);
      cf = mdr_Create (nc, 1);
      for (i = 0; i < nc; i++)
      {
        MdiV0(ic,i) = 0;
        if (MdrV0(pb_cl,i) != -create_inf ())
          MdiV0(ic,i) = 1;
        if (MdrV0(pb_cu,i) !=  create_inf ())
          MdiV0(ic,i) = MdiV0(ic,i) + 2;
        if (MdrV0(pb_cu,i) == MdrV0(pb_cl,i))
          MdiV0(ic,i) = MdiV0(ic,i) + 2;
      }
    }

    // integer parameters
    ipar[2 -1] = pb_mes;
    ipar[3 -1] = pb_mtesx;
    ipar[4 -1] = pb_mtesf;
    ipar[5 -1] = maxiter;
    ipar[6 -1] = maxiter;
    if (louts > 1)
      ipar[7 -1] = -2;  // maximum verbosity
    else
      ipar[7 -1] =  0;  // no messages

    // double parameter
    rpar[1 -1] = pb_tolx;
    rpar[2 -1] = pb_tolf;
    rpar[3 -1] = pb_tolb;
    rpar[4 -1] = pb_tolg;
    rpar[5 -1] = pb_told;
    rpar[6 -1] = pb_tols;
    rpar[7 -1] = pb_tolp;
    rpar[8 -1] = pb_eta;
    rpar[9 -1] = pb_xmax;

    if (louts > 1)
    {
      fptr = fopen (outs, "a");
      if (fptr)
      {
        t1 = clock ();
        fprintf (fptr,
                 "RLaB: using solver PBUNL/PNEWL for constrained"
                     " optimization in multidimensions.\n"
                );
        fprintf (fptr,
                 "RLaB: Messages from the solver follow.\n"
                );
        fflush(fptr);
      }
    }
    //
    // call PBUNL/PNEWL if conditions are present
    //
    if (nc > 0)
    {
      if (imethod == 1)
        PBUNL (&nf, &na, &nb, &nc, MDRPTR(w), NULL, NULL, NULL,
                MDRPTR(cf), MDIPTR(ic), MDRPTR(pb_cl), MDRPTR(pb_cu), MDRPTR(cg),
                       MDIPTR(ia), MDRPTR(ra), ipar, rpar, &fp, &gmax, &iterm,
                              outs, &louts);
      else
        PNEWL (&nf, &na, &nb, &nc, MDRPTR(w), NULL, NULL, NULL,
                MDRPTR(cf), MDIPTR(ic), MDRPTR(pb_cl), MDRPTR(pb_cu), MDRPTR(cg),
                       MDIPTR(ia), MDRPTR(ra), ipar, rpar, &fp, &gmax,
                              &ihes, &iterm,
                              outs, &louts);
    }
    else
    {
      if (imethod == 1)
        PBUNU (&nf, &na, MDRPTR(w), MDIPTR(ia), MDRPTR(ra),
                ipar, rpar, &fp, &gmax, &iterm,
                outs, &louts);
      else
        PNEWU (&nf, &na, MDRPTR(w), MDIPTR(ia), MDRPTR(ra),
                ipar, rpar, &fp, &gmax, &ihes, &iterm,
                outs, &louts);
    }

    if (fptr)
    {
      fflush(fptr);
      fprintf (fptr, "RLaB: solver PBUNL/PNEWL finished the"
          " optimization search with code ");
      fprintf (fptr, "%i", iterm);
      if (iterm > 0)
        fprintf (fptr, " (success),\nRLaB: ");
      else
        fprintf (fptr, " (failure).\nRLaB: ");
      if (iterm == 1)
        fprintf (fptr, "|X - XO| was less than or equal to"
            " TOLX in MTESX subsequent iterations.\n");
      else if (iterm == 2)
        fprintf (fptr, "|F - FO| was less than or equal to"
            " TOLF in MTESF subsequent iterations.\n");
      else if (iterm == 3)
        fprintf (fptr, "F is less than or equal to TOLB.\n");
      else if (iterm == 4)
        fprintf (fptr, "GMAX is less than or equal to TOLG.\n");
      else if (iterm == 11)
        fprintf (fptr, "NFV exceeded MFV.\n");
      else if (iterm == 12)
        fprintf (fptr, "NIT exceeded MIT.\n");
      else
        fprintf (fptr, "Unknown code. Check the source code of 811.f.\n");
      t2 = clock ();
      timer = (t2 - t1) / 1e6;
      fprintf (fptr, "RLaB: optimization lasted %g sec.\n", timer);
      fclose (fptr);
    }

    // clean-up
    if (ia)
      mdr_Destroy (ia);
    if (ra)
      mdr_Destroy (ra);
    if (nc > 0)
    {
      mdr_Destroy (ic);
      mdr_Destroy (cf);
      mdr_Destroy (cg);
    }

    if (iterm < 0)
    {
      mdr_Destroy (w);
      w = mdr_Create (0,0);
    }
  }
  else if (imethod == 1 || imethod == 2 || imethod == 3)
  {
    //
    // ProxyBundle: dimf >= 1
    //
    int    nb = 0, nc, nf = cmx_dimx, na;
    int    nia, nra, ipar[7], iterm , iext=-1;
    double fp = 0, gmax;
    double rpar[7];

    if (pb_isfun)
    {
      iext=0; //fun is given and is non-infinity
      pb_F = mdr_Copy (fun);
    }

    MDR * ic=0, * ia=0, * ra=0, * cf=0, *cg=0, *af=0;

    if (cmx_dimg > 0)
    {
      cg = mdr_Transpose (pb_cg);
      nc = MNR(pb_cg);
    }
    else
      nc = 0;

    nf = cmx_dimx;
    na = cmx_dimf;

    af = mdr_Create (na, 1);

    nia = nf + na + 1;
    ia = mdi_Create(nia, 1);

    nra = (nf+na+8)*nf+2*na+4;
    ra  = mdr_Create(nra,1);

    //
    // Form linear constraint matrices
    //
    if (nc)
    {
      ic = mdi_Create(nc, 1);
      cf = mdr_Create (nc, 1);
      for (i = 0; i < nc; i++)
      {
        MdiV0(ic,i) = 0;
        if (MdrV0(pb_cl,i) != -create_inf ())
          MdiV0(ic,i) = 1;
        if (MdrV0(pb_cu,i) !=  create_inf ())
          MdiV0(ic,i) = MdiV0(ic,i) + 2;
        if (MdrV0(pb_cu,i) == MdrV0(pb_cl,i))
          MdiV0(ic,i) = MdiV0(ic,i) + 2;
      }
    }

    // integer parameters
    ipar[1 -1] = pb_met;
    ipar[2 -1] = 2;       // use Powell's correction if negative curvature occurs
    ipar[3 -1] = 1;       // perform a restart after unsuccessive variable metric update
    ipar[4 -1] = pb_mes;
    ipar[5 -1] = maxiter;
    ipar[6 -1] = maxiter;
    if (louts > 1)
      ipar[7 -1] = -2;  // maximum verbosity
    else
      ipar[7 -1] =  0;  // maximum verbosity

    // double parameter
    rpar[1 -1] = pb_tolx;
    rpar[2 -1] = pb_tolf;
    rpar[3 -1] = pb_tolb;
    rpar[4 -1] = pb_tolg;
    rpar[5 -1] = pb_told;
    rpar[6 -1] = pb_tols;
    rpar[7 -1] = pb_xmax;

    if (louts > 1)
    {
      fptr = fopen (outs, "a");
      if (fptr)
      {
        t1 = clock ();
        fprintf (fptr,
                 "RLaB: using solver PMIN for constrained optimization"
                     " in multidimensions.\n"
                );
        fprintf (fptr,
                 "RLaB: Messages from the solver follow.\n"
                );
        fflush(fptr);
      }
    }
    //
    // call PMINL now
    //
    if (nc)
      PMINL (&nf, &na, &nb, &nc, MDRPTR(w), NULL, NULL, NULL,
              MDRPTR(cf), MDIPTR(ic), MDRPTR(pb_cl), MDRPTR(pb_cu), MDRPTR(cg), MDRPTR(af),
                     MDIPTR(ia), MDRPTR(ra), ipar, rpar, &fp, &gmax, &iext, &iterm,
                            outs, &louts);
    else
      PMINU (&nf, &na, MDRPTR(w), MDRPTR(af),
              MDIPTR(ia), MDRPTR(ra), ipar, rpar, &fp, &gmax, &iext, &iterm,
                     outs, &louts);

    if (fptr)
    {
      fflush(fptr);
      fprintf (fptr, "RLaB: solver PMIN finished the optimization search with code ");
      fprintf (fptr, "%i", iterm);
      if (iterm > 0)
        fprintf (fptr, " (success),\nRLaB: ");
      else if (iterm == -6)
        fprintf (fptr, " (questionable),\nRLaB: ");
      else
        fprintf (fptr, " (failure).\nRLaB: ");
      if (iterm == 1)
        fprintf (fptr, "|X - XO| was less than or equal to"
            " TOLX in MTESX subsequent iterations.\n");
      else if (iterm == 2)
        fprintf (fptr, "|F - FO| was less than or equal to"
            " TOLF in MTESF subsequent iterations.\n");
      else if (iterm == 3)
        fprintf (fptr, "F is less than or equal to TOLB.\n");
      else if (iterm == 4)
        fprintf (fptr, "GMAX is less than or equal to TOLG.\n");
      else if (iterm == 11)
        fprintf (fptr, "NFV exceeded MFV.\n");
      else if (iterm == 12)
        fprintf (fptr, "NIT exceeded MIT.\n");
      else if (iterm == -6)
      {
        fprintf (fptr, "Solver is having a bad day but the result still"
            " might be reliable.\n");
        iterm = 6;
      }
      else
        fprintf (fptr, "Unknown code. Check the source code of 811.f.\n");
      t2 = clock ();
      timer = (t2 - t1) / 1e6;
      fprintf (fptr, "RLaB: optimization lasted %g sec.\n", timer);
      fclose (fptr);
    }

    // clean-up
    mdr_Destroy (ia);
    mdr_Destroy (ra);
    mdr_Destroy (af);
    if (nc)
    {
      mdr_Destroy (ic);
      mdr_Destroy (cf);
      mdr_Destroy (cg);
    }
    if(pb_isfun) mdr_Destroy (pb_F);

    if (iterm < 0)
    {
      mdr_Destroy (w);
      w = mdr_Create (0,0);
    }
  }
  else
    rerror ("conmins: terrible internal error");

  // clean-up constraint arrays
  if (cmx_islin)
  {
    mdr_Destroy (pb_cg);
    mdr_Destroy (pb_cl);
    mdr_Destroy (pb_cu);
  }

  mdr_Destroy (cmx_icntyp);
  mdr_Destroy (fun);

  // clean-up rlab stuff
  MDPTR(xmdr) = 0;
  ent_DecRef (xent);
  ent_Destroy (xent);

  ent_Clean (e5);
  ent_Clean (pent);
  ent_Clean (e6);
  ent_Clean (eo);

  rent = ent_Create ();
  ent_data (rent) = w;
  ent_SetType (rent, MATRIX_DENSE_REAL);
  return (rent);
}


//
// conmax.f
//
static int
CMXFNSET (int * nparm, int * numgr, double * pttbl, int * iptb, int * indm,
          double * param, int * ipt, int * indfn, int * icntyp, double * confun)
{
  Ent * rent=0;
  MDR * retf, * retg, * retdf, * retdg;
  int i,j;

  //
  // tell fortran code about the types of constraints once
  //
  if (!cmx_typset)
  {
    cmx_typset = 1;
    for (i = 0; i < cmx_numgr; i++)
      icntyp[i] = MdiV0(cmx_icntyp,i);
  }

  MDPTR(xmdr) = param;

  if (*ipt == 0)
  {
    // vector function being optimized
    //
    // x for f(x)
    //

    if (pent)
      rent = ent_call_rlab_script_2args(cmx_f_name, xent, pent);
    else
      rent = ent_call_rlab_script_1arg (cmx_f_name, xent);

    if (ent_type(rent)!=MATRIX_DENSE_REAL)
      rerror (THIS_SOLVER ": " RLAB_ERROR_RHS_FUNC_MUST_RETURN_MDR);

    retf = ent_data (rent);

    if (SIZE (retf) != cmx_dimf)
      rerror (THIS_SOLVER ": " RLAB_ERROR_RHS_FUNC_INCORRECT_DIM);

    for (i = 0; i < cmx_dimf; i++)
      confun[i] = MdrV0(retf,i);

    // Try to clean up after calling the script if possible.
    ent_Clean (rent);


    if (cmx_indfn)
    {
      // jacobian df/dx

      if (pent)
        rent = ent_call_rlab_script_2args(cmx_df_name, xent, pent);
      else
        rent = ent_call_rlab_script_1arg (cmx_df_name, xent);

      if (ent_type(rent)!=MATRIX_DENSE_REAL)
        rerror (THIS_SOLVER ": " RLAB_ERROR_JACOBIAN_FUNC_MUST_RETURN_MDR);

      retdf = ent_data (rent);

      if (SIZE(retdf) !=  cmx_dimf * cmx_dimx )
        rerror (THIS_SOLVER ": " RLAB_ERROR_JACOBIAN_FUNC_INCORRECT_DIM);

      for (i = 0; i < cmx_dimf; i++)
        for (j = 0; j < cmx_dimx; j++)
          confun[(j+1)*(cmx_numgr) + i] = Mdr0(retdf,i,j);

      // Try to clean up after calling the script if possible.
      ent_Clean (rent);
    }
  }

  if (cmx_dimg > 0)
  {
    if (!cmx_islin)
    {
      if (pent)
        rent = ent_call_rlab_script_2args(cmx_g_name, xent, pent);
      else
        rent = ent_call_rlab_script_1arg (cmx_g_name, xent);

      if (ent_type(rent)!=MATRIX_DENSE_REAL)
        rerror (THIS_SOLVER ": " RLAB_ERROR_RHS_FUNC_MUST_RETURN_MDR);

      retg = ent_data (rent);

      if (SIZE(retg) != cmx_dimg)
        rerror (THIS_SOLVER ": " RLAB_ERROR_RHS_FUNC_INCORRECT_DIM);

      for (i = 0; i < cmx_dimg; i++)
        confun[i+cmx_dimf] = MdrV0(retg,i);

      // Try to clean up after calling the script if possible.
      ent_Clean (rent);

      //
      // check if jacobian comes from user function
      //
      if (cmx_indfn)
      {
        // jacobian dg/dx

        if (pent)
          rent = ent_call_rlab_script_2args(cmx_dg_name, xent, pent);
        else
          rent = ent_call_rlab_script_1arg (cmx_dg_name, xent);

        if (ent_type(rent)!=MATRIX_DENSE_REAL)
          rerror (THIS_SOLVER ": " RLAB_ERROR_JACOBIAN_FUNC_MUST_RETURN_MDR);

        retdg = ent_data (rent);

        if (SIZE(retdg) != cmx_dimg*cmx_dimx)
          rerror (THIS_SOLVER ": " RLAB_ERROR_JACOBIAN_FUNC_INCORRECT_DIM);

        for (i = 0; i < cmx_dimg; i++)
          for (j = 0; j < cmx_dimx; j++)
            confun[(j+1)*(cmx_numgr) + i + cmx_dimf] = Mdr0(retdg,i,j);

        // Try to clean up after calling the script if possible.
        ent_Clean (rent);
      }
    }
    else
    {
      //
      // linear constraints:  a.x + b <= 0
      //
      for (i = 0; i < cmx_dimg; i++)
      {
        confun[i+cmx_dimf] = 0;
        for (j = 0; j < MNC(cmx_A); j++)
          confun[i+cmx_dimf] += Mdr0(cmx_A,i,j) * param[j];
        confun[i+cmx_dimf] += MdrV0(cmx_B,i);
      }
      for (i = 0; i < cmx_dimg; i++)
        for (j = 0; j < cmx_dimx; j++)
          confun[(j+1)*(cmx_numgr) + i + cmx_dimf] = Mdr0(cmx_A,i,j);
    }
  }

  return 1;
}




//
// proximal bundle (811.f) fortran wrappers
//
int
PMFUNDER (int *nf, double * xval, double * f, double *g)
{
  int i;
  Ent *rent = 0;
  MDR *retm = 0;

  // x
  MDPTR(xmdr) = xval;

  if (pent)
    rent = ent_call_rlab_script_2args(cmx_f_name, xent, pent);
  else
    rent = ent_call_rlab_script_1arg (cmx_f_name, xent);

  if (ent_type(rent)!=MATRIX_DENSE_REAL)
    rerror (THIS_SOLVER ": " RLAB_ERROR_RHS_FUNC_MUST_RETURN_MDR);

  *f = class_double (rent);

  // Try to clean up if possible.
  ent_Clean (rent);

  // call rlab function df(x/,p/)
  if (pent)
    rent = ent_call_rlab_script_2args(cmx_df_name, xent, pent);
  else
    rent = ent_call_rlab_script_1arg (cmx_df_name, xent);

  if (ent_type(rent)!=MATRIX_DENSE_REAL)
    rerror (THIS_SOLVER ": " RLAB_ERROR_JACOBIAN_FUNC_MUST_RETURN_MDR);

  retm = ent_data (rent);

  if (SIZE(retm) != *nf)
    rerror (THIS_SOLVER ": " RLAB_ERROR_JACOBIAN_FUNC_INCORRECT_DIM);

  for (i = 0; i < *nf; i++)
    g[i] = MdrV0 (retm, i);

  // Try to clean up if possible.
  ent_Clean (rent);
  return 1;
}

int
PMFUN (int *nf, int *k, double * x, double * fa)
{
  Ent *rent = 0;
  static MDR *retm = 0;

  if (*k == 1)
  {
    //
    // call the rlab script the first time. save the resulting vector
    // and deliver its entries as requested.
    //

    // x
    MDPTR(xmdr) = x;

    // f(x/,p/)
    if (pent)
      rent = ent_call_rlab_script_2args(cmx_f_name, xent, pent);
    else
      rent = ent_call_rlab_script_1arg (cmx_f_name, xent);

    if (ent_type(rent)!=MATRIX_DENSE_REAL)
      rerror (THIS_SOLVER ": " RLAB_ERROR_RHS_FUNC_MUST_RETURN_MDR);

    retm = mdr_Copy (ent_data (rent));

    if (SIZE(retm)!= cmx_dimf)
      rerror (THIS_SOLVER ": " RLAB_ERROR_RHS_FUNC_INCORRECT_DIM);

    // Try to clean up after calling the script if possible.
    ent_Clean (rent);
  }

  //
  // return the requested row of f(x/,p/)
  //
  if (pb_isfun)
    *fa = MdrV1 (retm, *k) - MdrV1(pb_F, *k);
  else
    *fa = MdrV1 (retm, *k);

  // clean up after the last column
  if (*k == cmx_dimf)
    mdr_Destroy (retm);

  return 1;
}

int
PMDER (int *nf, int *k, double * x, double * g)
{
  int i;
  Ent *rent = 0;
  static MDR *retm = 0;

  if (*k == 1)
  {
    //
    // call the rlab script the first time. save the resulting vector
    // and deliver its entries as requested.
    //

    // x
    MDPTR(xmdr) = x;

    // call rlab function df(x/,p/)
    if (pent)
      rent = ent_call_rlab_script_2args(cmx_df_name, xent, pent);
    else
      rent = ent_call_rlab_script_1arg (cmx_df_name, xent);

    if (ent_type(rent)!=MATRIX_DENSE_REAL)
      rerror (THIS_SOLVER ": " RLAB_ERROR_JACOBIAN_FUNC_MUST_RETURN_MDR);

    retm = mdr_Copy( ent_data (rent) );

    if (SIZE(retm) != cmx_dimf*cmx_dimx)
      rerror (THIS_SOLVER ": " RLAB_ERROR_JACOBIAN_FUNC_INCORRECT_DIM);

    // Try to clean up after calling the script if possible.
    ent_Clean (rent);
  }

  //
  // return the requested row of f(x/,p/)
  //
  for (i = 0; i < *nf; i++)
    g[i] = Mdr0 (retm, *k-1, i);

  // clean up after the last row
  if (*k == cmx_dimf)
    mdr_Destroy (retm);

  return 1;
}

int
PMHES (int * nf, double *x, double *h)
{
  return 1;
}

