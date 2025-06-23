// odeb.c
// This file is a part of RLaB2 Rel.2
//  Copyright (C) 2005 Marijan Kostrun. project rlabplus.

/*
   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.

   This program is disin the hope that it will be useful,
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

#include "rfileio.h"
#include "rlab_solver_parameters_names.h"

//
// functions: their RLaB names and their Fortran arguments
//
int acfsub_ (int *, double *, double *, double *, double *);
int fsub_ (int *, double *, double *, double *);
int colfsub (double *, double *, double *, double *);
Ent *fsub_fname;

int acdfsub_ (int *, double *, double *, double *, double *);
int dfsub_ (int *, double *, double *, double *);
int coldfsub (double *, double *, double *, double *);
Ent *dfsub_fname;

int acgsub_ (int *, int *, double *, double *, double *);
int gsub_ (int *, double *, double *, double *);
int colgsub (int *, double *, double *);
Ent *gsubl_fname;		// left boundary condition, gl(u,p) = 0
Ent *gsubr_fname;		// right boundary condtion, gr(u,p) = 0

int acdgsub_ (int *, int *, double *, double *, double *);
int dgsub_ (int *, double *, double *, double *);
int coldgsub (int *, double *, double *);
Ent *dgsubl_fname;		// left boundary condition, gl(u,p) = 0
Ent *dgsubr_fname;		// right boundary condtion, gr(u,p) = 0

int colguess (double *, double *, double *, double *);

//
// Global variables
//
static int giveps = 0;
static int islin = 0;
static int giveu = 0;

static MDR *z=0, *x=0, *teps=0;
static Ent *ent_z=0, *ent_param=0, *ent_x=0, *ent_teps=0;

static MDR *u=0, *xx=0;
static int nlbc, nmax=0, ncomp=0;

//
// interface to acdc.f and mirkdc.f
//
#undef THIS_SOLVER
#define THIS_SOLVER "odebv"
Ent *
ent_odebv (int nargs, Datum args[])
{
  int i, j, nr5, nc5;
  Ent *eo=0, *e1=0, *e2=0, *e4=0, *e5=0, *e10=0;

  double aleft = 0, aright = 0, eps = 0.5, epsmin = 0.5,  ddummy, tolmin=1e-16;

  int iflag, givmsh = 0, nmsh = 0, ntol = 1;
  int nfxpnt = 0, nmax_new = 8192, ncomp_new=0;

  MDR *fixpnt=0, *wrk=0, *x5=0, *w=0, *x10=0, *iwrk=0, *ltol=0;
  int iflbvp, lwrkfl, lwrkin;
  int info, idummy;

  FILE *fptr = NULL;
  FILE *rlab_stderr = (!RLAB_STDERR_DS) ? stderr : RLAB_STDERR_DS;

  time_t timer = clock ();

  char *outs=0;
  int louts=0, imethod=0;
  MDR *tol;

  ListNode *node;

  //
  // init some values
  //
  ent_teps = 0;
  giveps = 0;
  giveu = 0;

  //
  // Check arguments.
  //
  if (nargs < 9)
  {
    fprintf (rlab_stderr,
             THIS_SOLVER ": ODE boundary value problem solvers. Format:\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ":   y = odebv(dim,nlbc,f,df,/p/,X,gL,gR,dgL,dgR /,Eps,options/),\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": where, f = f(x,u/,p/), is the vector function of first\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": derivatives,  y' = f(x,y,/p,eps/) , 'p' is the parameter \n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": passed to all functions, X = [a,b]\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": or X = [a,..b] to set the mesh, X = [x,y(0)]\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": to set the mesh and the initial approximation to the \n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": solution.  The roots of the functions gL(y/,p,eps/)=0\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": and gR(y/,p,eps/)=0 determine the boundary conditions and\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": For 'acdc' method, Eps=[eps, mineps] is the range of parameters\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": with respect to which the continuation is done \n");
    fprintf (rlab_stderr,
             "See also   odeparams  list of functions.\n");
    rerror ("requires at least 9 arguments");
  }

  //
  // dimu: dimension of the problem
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) != MATRIX_DENSE_REAL)
    rerror (THIS_SOLVER ": Improper 1st argument, ncomp!");
  ncomp_new = (int) class_double (e1);
  if (ncomp_new < 1)
    rerror (THIS_SOLVER ": Dimension of the problem cannot be less than unity!");

  //
  // nlbc: number of left boundary condtions
  //
  e2 = bltin_get_ent (args[1]);
  if (ent_type (e2) != MATRIX_DENSE_REAL)
    rerror (THIS_SOLVER ": Improper 2nd argument, nlbc!");
  nlbc = (int) class_double (e2);
  if (nlbc == 0 || nlbc >= ncomp_new)
    rerror (THIS_SOLVER ": this is not an initial value ODE solver!\n");

  tol = mdr_Create (ncomp_new,1);
  for (i=0; i < ncomp_new; i++)
    MdrV0(tol,i) = tolmin;

  //
  // options
  //
  if (nargs == 12 || nargs == 11)
  {
    eo = bltin_get_ent (args[nargs-1]);
    if (ent_type (eo) == BTREE)
    {
      // method
      node = btree_FindNode (ent_data (eo), RLAB_NAME_ODEBV_NMAX);
      if (node != 0)
      {
        idummy = (int) class_double (var_ent (node));
        if (idummy > 0)
          nmax_new = idummy;
      }
      // tol
      node = btree_FindNode (ent_data (eo), RLAB_NAME_ODEBV_TOL);
      if (node != 0)
      {
        MDR * dd = class_matrix_real ( var_ent (node) );
        if (SIZE(dd) >= 1)
          for (i = 0; i < ncomp_new; i++)
          {
            ddummy = MdrV0(dd, MIN(i, SIZE(dd)-1));
            if (ddummy > 0)
              MdrV0(tol,i) = ddummy;
          }
      }
      // method
      node = btree_FindNode (ent_data (eo), RLAB_NAME_GEN_IMETHOD);
      if (node != 0)
      {
        idummy = (int) class_double (var_ent (node));
        if (idummy >= 0 && idummy <= 4)
          imethod = idummy;
      }
      // standard output
      node = btree_FindNode (ent_data (eo), RLAB_NAME_GEN_STDOUT);
      if (node != 0)
      {
        if (ent_type(var_ent (node))== MATRIX_DENSE_STRING)
        {
          outs = class_char_pointer(var_ent (node));
          louts = isvalidstring (outs);
          if (louts<1)
          {
            outs = 0;
            louts = 0;
          }
        }
      }
    }
    ent_Clean (eo);
  }

  //
  // same definitions: avoid repeated allocation of large memory blocks if user is
  // going to use this solver often in her session
  //
  if (nmax_new != nmax)
  {
    if (xx)
      mdr_Destroy (xx);
    nmax = nmax_new;
    xx = mdr_Create (1, nmax);
  }
  if ((nmax_new != nmax) || (ncomp_new != ncomp))
  {
    if (u)
      mdr_Destroy (u);
    nmax = nmax_new;
    ncomp = ncomp_new;
    u  = mdr_Create (ncomp, nmax);
  }

  //
  // fsub: u' = fsub(x,u,p)
  //
  fsub_fname = bltin_get_ent(args[2]);
  if (!isfuncent(fsub_fname))
    rerror (THIS_SOLVER ": " RLAB_ERROR_ARG3_FUNC_VAR "\n");

  //
  // dfsub: dfsub(x,u,p) = \partial fsub(x,u,p) / \partial u, jacobian of fsub.
  //
  dfsub_fname = bltin_get_ent(args[3]);
  if (!isfuncent(dfsub_fname))
    rerror (THIS_SOLVER ": " RLAB_ERROR_ARG4_FUNC_VAR "\n");

  //
  // parameter entity
  //
  ent_param = 0;
  e4 = bltin_get_ent (args[4]);
  if (ent_type (e4) != UNDEF)
    ent_param = ent_Copy (e4);

  //
  // Get interval / mesh / mesh + initial guess for u
  //
  e5 = bltin_get_ent (args[5]);
  if (ent_type (e5) != MATRIX_DENSE_REAL)
  {
    fprintf (stdout, THIS_SOLVER ": 5th argument has to be one of the following:\n");
    fprintf (stdout,
             THIS_SOLVER ":   (1) [a,b], interval on which the solution is sought,\n");
    fprintf (stdout,
             THIS_SOLVER ":   (2) [a,x2,..,b], interval with the initial mesh,\n");
    fprintf (stdout,
             THIS_SOLVER ":   (3) [ a,x2,..,b; u1(a),u1(x2),..u1(b); .. ], mesh and the\n");
    fprintf (stdout, THIS_SOLVER ":   initial guess for the solution u=u(x,p).\n");
    fprintf (stdout,
             THIS_SOLVER ": Note: if initial mash is given than it is kept (improved\n");
    fprintf (stdout, THIS_SOLVER ": upon) in the solution.\n");
    rerror ("Improper 5th argument!");
  }
  x5 = class_matrix_real (e5);
  nr5 = MNR (x5);
  nc5 = MNC (x5);
  if (nr5 != 1 && nr5 != ncomp + 1)
    rerror (THIS_SOLVER ": Improper 6th argument. Dimension mismatch!");
  aleft  = Mdr1 (x5, 1, 1);
  aright = Mdr1 (x5, 1, nc5);
  nfxpnt = nc5 - 2;
  fixpnt = mdr_Create (1, MAX(1, nfxpnt));
  if (nfxpnt > 0)
  {
    // mesh is given: make note and prepare xx
    givmsh = 1;
    nmsh = nc5;
    for (i = 1; i <= nc5; i++)
    {
      Mdr1 (xx, 1, i) = Mdr1 (x5, 1, i);
      if (i != 1 && i != nc5)
        Mdr1 (fixpnt, 1, i - 1) = Mdr1 (x5, 1, i);
    }
    if (nr5 > 1)
    {
      // initial guess for u=u(x,p) is given: make note and prepare u
      giveu = 1;
      for (j = 2; j <= nr5; j++)
        for (i = 1; i <= nc5; i++)
          Mdr1 (u, j - 1, i) = Mdr1 (x5, j, i);
    }
  }
  if (aleft == aright)
    rerror (THIS_SOLVER ": Improper boundaries in the 5th argument!");

  //
  // Get left and right boundary condition functions, gleft(u,p)=0 , gright(u,p)=0
  //
  gsubl_fname = bltin_get_ent(args[6]);
  if (!isfuncent(gsubl_fname))
    rerror (THIS_SOLVER ": " RLAB_ERROR_ARG7_FUNC_VAR "\n");

  gsubr_fname = bltin_get_ent(args[7]);
  if (!isfuncent(gsubr_fname))
    rerror (THIS_SOLVER ": " RLAB_ERROR_ARG8_FUNC_VAR "\n");

  //
  // get left/right differentials of boundary conditions:
  //      dgleft(p,u)=\partial gleft(u,p) / \partial u.
  //      dgright(p,u)=\partial gright(u,p) / \partial u.
  //
  dgsubl_fname = bltin_get_ent(args[8]);
  if (!isfuncent(dgsubl_fname))
    rerror (THIS_SOLVER ": " RLAB_ERROR_ARG9_FUNC_VAR "\n");

  dgsubr_fname = bltin_get_ent(args[9]);
  if (!isfuncent(dgsubr_fname))
    rerror (THIS_SOLVER ": " RLAB_ERROR_ARG10_FUNC_VAR "\n");

  //
  // Check whether we have parameter epsilon given
  //
  if (nargs >= 11)
  {
    e10 = bltin_get_ent (args[10]);
    if (ent_type (e10) == MATRIX_DENSE_REAL)
    {
      x10 = class_matrix_real (e10);
      if (SIZE(x10)==2)
      {
        //
        // epsilon is given: make note and take the values
        //
        imethod = 0;
        giveps = 1;
        eps = MdrV1 (x10, 1);
        epsmin = MdrV1 (x10, 2);
      }
    }
    ent_Clean (e10);
  }

  //
  // Set up ENTITIES for user-function. Inc the reference count once, for belonging
  // to odebv().  u' = f(x,p,u)
  // argument: x, local
  x = mdr_CreateEmpty (1,1);
  ent_x = ent_Assign_Rlab_MDR (x);
  ent_IncRef (ent_x);

  // argument: z = u(x,p), locally
  z = mdr_CreateEmpty (ncomp, 1);
  ent_z = ent_Assign_Rlab_MDR (z);
  ent_IncRef (ent_z);

  if (imethod == 0)
  {
    //
    // acdc
    //

    ntol = ncomp;

    //
    // prepare the arrays for acdc
    //
    ltol = mdi_Create(ncomp, 1);
    for (i = 0; i < ncomp; i++)
      MdiV0(ltol,i) = i + 1;

    lwrkfl =
        nmax * (4 * ncomp * ncomp + 11 * ncomp + 5) + 4 * ncomp * ncomp +
        11 * ncomp + 5;
    wrk  = mdr_Create (lwrkfl, 1);

    lwrkin = nmax * (ncomp + 2) + 3 * ncomp;
    iwrk = mdi_Create (lwrkin, 1);

    if (eps > epsmin)
    {
      teps = mdr_CreateEmpty (1,1);
      ent_teps = ent_Assign_Rlab_MDR (teps);
      ent_IncRef (ent_teps);
      giveps = 1;
    }

    // dealing with FORTRAN: if file name for stdout is given then
    // make it exist on file system. otherwise fortran will exit with
    // Fortran runtime error: File '.....' does not exist
    if (louts)
    {
      fptr = fopen (outs, "a");
      fclose  (fptr);
    }

    //
    // Call integrator
    //
    ACDC (&ncomp, &nlbc, &nmax, &aleft, &aright, &nfxpnt, MDRPTR(fixpnt),
           &ntol, MDIPTR(ltol), MDRPTR(tol), &islin, &givmsh, &giveu, &nmsh, MDRPTR(xx),
                                  &ncomp, MDRPTR(u), &nmax, &lwrkfl, MDRPTR(wrk), &lwrkin, MDIPTR(iwrk), &giveps,
          &eps, &epsmin, acfsub_, acdfsub_, acgsub_, acdgsub_,
          &iflbvp, outs, &louts);

    mdr_Destroy (ltol);

    //
    // print end of integration
    //
    if (louts)
      fptr = fopen (outs, "a");
    if (!iflbvp)
    {
      if (fptr)
      {
        fprintf (fptr, "RLaB: built-in integrator 'acdc' reports success.\n");
        // check the time and close the output stream
        timer -= clock ();
        fprintf (fptr,
                 "RLaB: ODE integration using 'acdc' lasted %g sec.\n",
                 -timer / 1e6);
        fclose (fptr);
      }
      w = mdr_Create (nmsh, ncomp + 1);
      for (i = 1; i <= nmsh; i++)
      {
        Mdr1 (w, i, 1) = Mdr1 (xx, 1, i);
        for (j = 1; j <= ncomp; j++)
          Mdr1 (w, i, j + 1) = Mdr1 (u, j, i);
      }
    }
    else
    {
      if (fptr)
      {
        fprintf (fptr, "RLaB: built-in integrator 'acdc' reports error %i\n",
                 iflbvp);
        // check the time and close the output stream
        timer -= clock ();
        fprintf (fptr,
                 "RLaB: ODE integration using 'acdc' lasted %g sec.\n",
                 -timer / 1e6);
        fclose (fptr);
      }
      w = mdr_Create (0, 0);
    }
  }
  else if (imethod == 4)
  {
    //
    // coldae in colsys mode
    //

    // ode problem and not algebraic
    int ny = 0;
    // it is the first order problem
    int *m = (int *) GC_MALLOC (ncomp * sizeof (int));
    for (i = 0; i < ncomp; i++)
      m[i] = 1;

    // figure out the bc. irrelevant
    MDR *zeta = mdr_Create (ncomp, 1);
    for (i=0; i<ncomp; i++)
    {
      if (i < nlbc)
        MdrV0 (zeta, i) = aleft;
      else
        MdrV0 (zeta, i) = aright;
    }
    // figure out the size of f-space and i-space
    lwrkin = nmax * (3 + 5 * ncomp);
    lwrkfl = nmax * (4 + 3 * ncomp + (5 + 4 * ncomp) * 5 * ncomp +
        (ncomp + nlbc) * 2 * ncomp +
        ncomp * (ncomp + 2) + 4 * ncomp);

    int *ipar = (int *) GC_MALLOC (12 * sizeof (int));
    ipar[1 - 1] = !islin;	// nonlinear problem
    ipar[2 - 1] = 0;		// use 4 points per collocation interval
    ipar[3 - 1] = nfxpnt + 1;	// number of points in initial mesh
    ipar[4 - 1] = ncomp;	// number of tolerances
    ipar[5 - 1] = lwrkfl;	// dimension of f-space
    ipar[6 - 1] = lwrkin;	// dimension of i-space
    // is there output messages to be written?
    if (louts > 1)
      ipar[7 - 1] = -1;		// yes, maximal
    else
      ipar[7 - 1] = 1;		// no
    // did user specify a mesh
    if (givmsh == 1)
      ipar[8 - 1] = 1;		// yes, fspace contains initial mesh
    else
      ipar[8 - 1] = 0;		// no, generate uniform initial mesh
    // did user specify initial guess
    if (giveu == 1)
      ipar[9 - 1] = 2;		// yes, fspace contains initial guess
    else
      ipar[9 - 1] = 0;		// no
    ipar[10 - 1] = 0;		// problem is regular
    // did user specify a fixed points (assumed to be initial mesh)
    if (givmsh == 1)
      ipar[11 - 1] = nfxpnt;	// yes
    else
      ipar[11 - 1] = 0;		// no
    ipar[12 - 1] = 0;		// it is not a dae but odebv problem

    ltol = mdi_Create(ncomp,1);
    for (i = 0; i < ncomp; i++)
      MdiV0(ltol,i) = i + 1;

    //
    // prepare the arrays for coldae
    //
    // iwrk <=> ispace
    iwrk = mdi_Create (lwrkin,1);
    // work <=> fspace
    wrk = mdr_Create (lwrkfl, 1);
    if (nfxpnt > 0)
    {
      // initial mesh in 'xx'
      for (i = 0; i < nmsh; i++)
        MdrV0 (wrk, i) = MdrV0 (xx, i);
      // initial approximation in 'u'
      if (giveu)
      {
        for (i = 0; i < nmsh; i++)
          for (j = 0; j < ncomp; j++)
            MdrV0 (wrk, nmsh + i * ncomp + j) = Mdr0 (u, j, i);
        MdiV0(iwrk,0) = nmsh - 1;
      }
    }

    // set-up run-time messages
    //
    if (louts > 1)
    {
      fptr = fopen (outs, "a");
      if (fptr)
        timer = clock ();
      else
        louts = 0;
    }
    //
    // Call integrator
    //
    COLDAE (&ncomp, &ny, m, &aleft, &aright, MDRPTR(zeta), ipar, MDIPTR(ltol),
             MDRPTR(tol), MDRPTR(fixpnt), MDIPTR(iwrk), MDRPTR(wrk), &iflag,
             colfsub, coldfsub, colgsub, coldgsub, NULL, outs, &louts);
    //
    // print end of integration
    //
    if (iflag == 1)
    {
      int ns = MdiV0(iwrk,0) + 1;
      w = mdr_Create (ns, ncomp + 1);
      for (i = 0; i < ns; i++)
      {
        Mdr0 (w, i, 0) = MdrV0 (wrk, i);
        for (j = 0; j < ncomp; j++)
          Mdr0 (w, i, j + 1) = MdrV0 (wrk, i * ncomp + j + ns);
      }
      if (fptr)
      {
        fprintf (fptr, "RLaB: built-in integrator 'coldae' reports success.\n");
        // check the time and close the output stream
        timer -= clock ();
        fprintf (fptr,
                 "RLaB: ODE integration using 'coldae' lasted %g sec.\n",
                 -timer / 1e6);
        fflush (fptr);
        fclose (fptr);
      }
    }
    else
    {
      if (fptr)
      {
        fprintf (fptr, "RLaB: built-in integrator 'coldae' reports error %i\n",
                 iflag);
        // check the time and close the output stream
        timer -= clock ();
        fprintf (fptr,
                 "RLaB: ODE integration using 'coldae' lasted %g sec.\n",
                 -timer / 1e6);
        fclose (fptr);
      }
      w = mdr_Create (0, 0);
    }
    //
    // clean after coldae
    //
    GC_FREE (m);
    GC_FREE (ipar);
    mdr_Destroy (ltol);
    mdr_Destroy (zeta);
  }
  else
  {
    // put imethod to 2,4,6 instead of 1,2,3
    imethod = 2 * imethod;

    //
    // prepare the arrays for mirkdc
    //
    int output_control = 0;
    int mxnsub = nmax - 1;
    int mxs = 10;		// max number of stages of rk method inside mirkdc
    wrk = mdr_Create (2 * mxnsub * mxs * ncomp
        + 5 * ncomp * mxnsub + 2 * pow (ncomp, 2) * mxnsub
        + 6 * ncomp + 4 * pow (ncomp, 2) + 2 * pow (ncomp,
        2) * mxs, 1);
    //iwrk = (int *) GC_MALLOC (ncomp * (mxnsub + 1) * sizeof (int));
    iwrk = mdi_Create (ncomp * (mxnsub + 1),1);

    //
    // set-up run-time messages
    //
    if (louts > 1)
    {
      fptr = fopen (outs, "a");
      if (fptr)
      {
        // start the timer
        timer = clock ();
        // write down boring info
        fprintf (fptr,
                 "RLaB: using built-in integrator 'mirkdc' for boundary value problem.\n");
        fprintf (fptr, "RLaB: the messages from the solver follow\n");
        output_control = 0;
        fclose (fptr);
      }
      else
        louts = 0;
    }
    //
    // Call integrator
    //
    int nsub = MAX(nmsh - 1, 2);	// if xx->d contains the mesh this is its size
    int sol_input = 2;
    MIRKDC (&imethod, MDRPTR(tol), &ncomp, &nsub, &mxnsub,
            &givmsh, MDRPTR(xx), &sol_input, MDRPTR(u), &ncomp, &output_control, &info,
            MDIPTR(iwrk), MDRPTR(wrk), &aleft, &aright, &nlbc,
            fsub_, gsub_, dfsub_, dgsub_, outs, &louts);
    //
    // print end of integration
    //
    if (louts > 1)
      fptr = fopen (outs, "a");
    if (!info)
    {
      if (fptr)
      {
        fprintf (fptr, "RLaB: built-in integrator 'mirkdc' reports success.\n");
        // check the time and close the output stream
        timer -= clock ();
        fprintf (fptr,
                 "RLaB: ODE integration using 'mirkdc' lasted %g sec.\n",
                 -timer / 1e6);
        fclose (fptr);
      }
      w = mdr_Create (nsub + 1, ncomp + 1);
      for (i = 1; i <= nsub + 1; i++)
      {
        Mdr1 (w, i, 1) = Mdr1 (xx, 1, i);
        for (j = 1; j <= ncomp; j++)
          Mdr1 (w, i, j + 1) = Mdr1 (u, j, i);
      }
    }
    else
    {
      if (fptr)
      {
        fprintf (fptr, "RLaB: built-in integrator 'mirkdc' reports error %i\n",
                 info);
        // check the time and close the output stream
        timer -= clock ();
        fprintf (fptr,
                 "RLaB: ODE integration using 'mirkdc' lasted %g sec.\n",
                 -timer / 1e6);
        fclose (fptr);
      }
      w = mdr_Create (0, 0);
    }
  }
  //
  // Clean Up
  //
  mdr_Destroy (tol);
  mdr_Destroy (wrk);
  mdr_Destroy (iwrk);
  mdr_Destroy (fixpnt);
//   mdr_Destroy (xx);
//   mdr_Destroy (u);

  MDPTR(x) = 0;
  ent_DecRef (ent_x);
  ent_Destroy (ent_x);

  MDPTR(z) = 0;
  ent_DecRef (ent_z);
  ent_Destroy (ent_z);

  if (ent_teps)
  {
    MDPTR(teps) = 0;
    ent_DecRef  (ent_teps);
    ent_Destroy (ent_teps);
  }

  ent_Clean (e1);
  ent_Clean (e2);
  ent_Clean (e4);  ent_Clean (ent_param);
  ent_Clean (e5);
  ent_Clean (fsub_fname);
  ent_Clean (dfsub_fname);
  ent_Clean (gsubl_fname);
  ent_Clean (gsubr_fname);
  ent_Clean (dgsubl_fname);
  ent_Clean (dgsubr_fname);

  return ent_Assign_Rlab_MDR(w);
}

//
// User callable functions for ACDC
//
//    AAA    CCC   DDDD   CCC
//   A   A  C   C  D   D C   C
//   AAAAA  C      D   D C
//   A   A  C   C  D   D C   C
//   A   A   CCC   DDDD   CCC
//
// The interface to the user-specified function in fortran
//
int
acfsub_ (int *dummy, double *newx, double *newz, double *f, double *epsilon)
{
  int j=0;
  Ent **e = (Ent **) GC_MALLOC(2 * sizeof(Ent *));
  Ent *rent = 0;
  MDR *retm = 0;

  //
  // Put x, u, and u' into  fsub_args.
  //
  MDPTR (x) = (void *) newx;
  MDPTR (z) = (void *) newz;

  if (ent_param)
    e[j++] = ent_param;
  if (giveps)
  {
    MDPTR (teps) = (void *) epsilon;
    e[j++] = ent_teps;
  }

  switch (j)
  {
    case 0:
      rent = ent_call_rlab_script_2args(fsub_fname, ent_x, ent_z);
      break;

    case 1:
      rent = ent_call_rlab_script_3args(fsub_fname, ent_x, ent_z, e[0]);
      break;

    case 2:
      rent = ent_call_rlab_script_4args(fsub_fname, ent_x, ent_z, e[0], e[1]);
      break;
  }

  if (ent_type(rent)!=MATRIX_DENSE_REAL)
    rerror (THIS_SOLVER ": " RLAB_ERROR_RHS_FUNC_MUST_RETURN_MDR);

  retm = ent_data(rent);

  if (SIZE(retm) != ncomp)
    rerror (THIS_SOLVER ": " RLAB_ERROR_RHS_FUNC_INCORRECT_DIM);

  for (j = 0; j < ncomp; j++)
    f[j] = MdrV0 (retm, j);

  e[0] = 0;
  e[1] = 0;
  GC_FREE (e);
  ent_Clean(rent);
  return (1);
}

//
// jacobian function for the acdc solver
//
int
acdfsub_ (int *dummy, double *newx, double *newz, double *df, double *epsilon)
{
  int j=0;
  Ent **e = (Ent **) GC_MALLOC(2 * sizeof(Ent *));
  Ent *rent = 0;
  MDR *retm = 0;

  MDPTR(x) = (void *) newx;
  MDPTR(z) = (void *) newz;

  if (ent_param)
    e[j++] = ent_param;
  if (giveps)
  {
    MDPTR (teps) = (void *) epsilon;
    e[j++] = ent_teps;
  }

  switch (j)
  {
    case 0:
      rent = ent_call_rlab_script_2args(dfsub_fname, ent_x, ent_z);
      break;

    case 1:
      rent = ent_call_rlab_script_3args(dfsub_fname, ent_x, ent_z, e[0]);
      break;

    case 2:
      rent = ent_call_rlab_script_4args(dfsub_fname, ent_x, ent_z, e[0], e[1]);
      break;
  }
  if (ent_type(rent)!=MATRIX_DENSE_REAL)
    rerror (THIS_SOLVER ": " RLAB_ERROR_RHS_FUNC_MUST_RETURN_MDR);

  retm = ent_data(rent);

  if (SIZE(retm) != ncomp*ncomp)
    rerror (THIS_SOLVER ": " RLAB_ERROR_JACOBIAN_FUNC_INCORRECT_DIM);

  for (j=0; j<ncomp*ncomp; j++)
    df[j] = MdrV0 (retm, j);

  e[0] = 0;
  e[1] = 0;
  GC_FREE (e);
  ent_Clean (rent);
  return (1);
}

//
// boundary value function for the solver
//
int
acgsub_ (int *i, int *dummy, double *newz, double *g, double *epsilon)
{
  int j=0;
  Ent **e = (Ent **) GC_MALLOC(2 * sizeof(Ent *));
  Ent *rent = 0;
  static MDR *retm = 0;

  MDPTR(z) = (void *) newz;

  if (*i == 1)
  {
    //
    // call left boundary condition function once, store tis result in static ent/mdr
    //
    if (ent_param)
      e[j++] = ent_param;
    if (giveps)
    {
      MDPTR (teps) = (void *) epsilon;
      e[j++] = ent_teps;
    }

    switch (j)
    {
      case 0:
        rent = ent_call_rlab_script_1arg (gsubl_fname, ent_z);
        break;

      case 1:
        rent = ent_call_rlab_script_2args(gsubl_fname, ent_z, e[0]);
        break;

      case 2:
        rent = ent_call_rlab_script_3args(gsubl_fname, ent_z, e[0], e[1]);
        break;
    }
    if (ent_type(rent)!=MATRIX_DENSE_REAL)
      rerror (THIS_SOLVER ": " RLAB_ERROR_RHS_FUNC_MUST_RETURN_MDR);

    retm = mdr_Copy (ent_data(rent));

    if (SIZE(retm) != nlbc)
      rerror (THIS_SOLVER ": " RLAB_ERROR_BV_FUNC_INCORRECT_DIM);

    ent_Clean (rent);
  }
  else if (*i == nlbc+1)
  {
    //
    // call right boundary condition once, store its result in static ent/mdr
    //
    if (ent_param)
      e[j++] = ent_param;
    if (giveps)
    {
      MDPTR (teps) = (void *) epsilon;
      e[j++] = ent_teps;
    }

    switch (j)
    {
      case 0:
        rent = ent_call_rlab_script_1arg (gsubr_fname, ent_z);
        break;

      case 1:
        rent = ent_call_rlab_script_2args(gsubr_fname, ent_z, e[0]);
        break;

      case 2:
        rent = ent_call_rlab_script_3args(gsubr_fname, ent_z, e[0], e[1]);
        break;
    }
    if (ent_type(rent)!=MATRIX_DENSE_REAL)
      rerror (THIS_SOLVER ": " RLAB_ERROR_RHS_FUNC_MUST_RETURN_MDR);

    retm = mdr_Copy (ent_data(rent));

    if (SIZE(retm) != ncomp-nlbc)
      rerror (THIS_SOLVER ": " RLAB_ERROR_BV_FUNC_INCORRECT_DIM);

    ent_Clean (rent);
  }

  if (*i <= nlbc)
  {
    *g = MdrV1 (retm, *i);
    if (*i == nlbc)
      mdr_Destroy (retm);
  }
  else
  {
    *g = MdrV1 (retm, *i-nlbc);
    if (*i == ncomp)
      mdr_Destroy (retm);
  }

  e[0] = 0;
  e[1] = 0;
  GC_FREE (e);
  return (1);
}

//
// differential of the boundary value function for the solver
//
int
acdgsub_ (int *i, int *dummy, double *newz, double *dg, double *epsilon)
{
  int j=0;
  Ent **e = (Ent **) GC_MALLOC(2 * sizeof(Ent *));
  Ent *rent = 0;
  static MDR *retm = 0;

  //
  // Copy 'u' to list of arguments, 'p' is already there
  //
  MDPTR(z) = (void *) newz;

  if (*i == 1)
  {
    //
    // call left boundary condition function once, store tis result in static ent/mdr
    //
    if (ent_param)
      e[j++] = ent_param;
    if (giveps)
    {
      MDPTR (teps) = (void *) epsilon;
      e[j++] = ent_teps;
    }

    switch (j)
    {
      case 0:
        rent = ent_call_rlab_script_1arg (dgsubl_fname, ent_z);
        break;

      case 1:
        rent = ent_call_rlab_script_2args(dgsubl_fname, ent_z, e[0]);
        break;

      case 2:
        rent = ent_call_rlab_script_3args(dgsubl_fname, ent_z, e[0], e[1]);
        break;
    }
    if (ent_type(rent)!=MATRIX_DENSE_REAL)
      rerror (THIS_SOLVER ": " RLAB_ERROR_RHS_FUNC_MUST_RETURN_MDR);

    retm = mdr_Copy (ent_data(rent));

    if (SIZE(retm) != nlbc*ncomp)
      rerror (THIS_SOLVER ": " RLAB_ERROR_BV_FUNC_INCORRECT_DIM);

    ent_Clean(rent);
  }
  else if (*i == nlbc+1)
  {
    //
    // call right boundary condition once, store its result in static ent/mdr
    //
    if (ent_param)
      e[j++] = ent_param;
    if (giveps)
    {
      MDPTR (teps) = (void *) epsilon;
      e[j++] = ent_teps;
    }

    switch (j)
    {
      case 0:
        rent = ent_call_rlab_script_1arg (dgsubr_fname, ent_z);
        break;

      case 1:
        rent = ent_call_rlab_script_2args(dgsubr_fname, ent_z, e[0]);
        break;

      case 2:
        rent = ent_call_rlab_script_3args(dgsubr_fname, ent_z, e[0], e[1]);
        break;
    }
    if (ent_type(rent)!=MATRIX_DENSE_REAL)
      rerror (THIS_SOLVER ": " RLAB_ERROR_RHS_FUNC_MUST_RETURN_MDR);

    retm = mdr_Copy (ent_data(rent));

    if (SIZE(retm) != (ncomp-nlbc)*ncomp)
      rerror (THIS_SOLVER ": " RLAB_ERROR_BV_FUNC_INCORRECT_DIM);

    ent_Clean(rent);
  }

  if (*i <= nlbc)
  {
    for (j = 1; j <= ncomp; j++)
      dg[j-1] = Mdr1 (retm, *i, j);
    if (*i == nlbc)
      mdr_Destroy (retm);
  }
  else
  {
    for (j = 1; j <= ncomp; j++)
      dg[j-1] = Mdr1 (retm, *i-nlbc, j);
    if (*i == ncomp)
      mdr_Destroy (retm);
  }

  e[0] = 0;
  e[1] = 0;
  GC_FREE (e);
  return (1);
}


//
// User callable functions for mirkdc.f
//
// MM MM III RRR  K   K DDDD   CCC
// M M M  I  R  R K  K  D   D C   C
// M M M  I  RRR  KKK   D   D C
// M   M  I  R  R K  K  D   D C   C
// M   M III R  R K   K DDDD   CCC
//
// The interface to the user-specified function in fortran for mirkdc
//
int
fsub_ (int *dummy, double *newx, double *newz, double *f)
{
  int j;
  Ent *rent = 0;
  MDR *retm = 0;

  MDPTR(x) = (void *) newx;
  MDPTR(z) = (void *) newz;

  if (ent_param)
    rent = ent_call_rlab_script_3args(fsub_fname, ent_x, ent_z, ent_param);
  else
    rent = ent_call_rlab_script_2args(fsub_fname, ent_x, ent_z);

  if (ent_type(rent)!=MATRIX_DENSE_REAL)
    rerror (THIS_SOLVER ": rhs function must return MATRIX-DENSE-REAL");

  retm = ent_data(rent);
  if (SIZE(retm) != ncomp)
    rerror (THIS_SOLVER ": incorrectly dimensioned derivitive (rhs vector)");

  for (j=0; j<ncomp; j++)
    f[j] = MdrV0 (retm, j);

  ent_Clean (rent);
  return (1);
}

//
// jacobian function for the mirkdc solver
//
int
dfsub_ (int *dummy, double *newx, double *newz, double *df)
{
  int j;
  Ent *rent = 0;
  MDR *retm = 0;

  MDPTR(x) = (void *) newx;
  MDPTR(z) = (void *) newz;

  if (ent_param)
    rent = ent_call_rlab_script_3args(dfsub_fname, ent_x, ent_z, ent_param);
  else
    rent = ent_call_rlab_script_2args(dfsub_fname, ent_x, ent_z);

  if (ent_type(rent)!=MATRIX_DENSE_REAL)
    rerror (THIS_SOLVER ": rhs function must return MATRIX-DENSE-REAL");

  retm = ent_data(rent);
  if (SIZE(retm) != ncomp*ncomp)
    rerror (THIS_SOLVER ": " RLAB_ERROR_JACOBIAN_FUNC_INCORRECT_DIM);

  for (j=0; j<ncomp*ncomp; j++)
    df[j] = MdrV0 (retm, j);

  ent_Clean (rent);
  return (1);
}

//
// boundary value function for the mirkdc solver
//
int
gsub_ (int *dummy, double *newza, double *newzb, double *bc)
{
  int j;
  Ent *rent = 0;
  MDR *retm = 0;

  //
  // left boundary
  //
  MDPTR(z) = (void *) newza;

  if (ent_param)
    rent = ent_call_rlab_script_2args(gsubl_fname, ent_z, ent_param);
  else
    rent = ent_call_rlab_script_1arg (gsubl_fname, ent_z);

  if (ent_type(rent)!=MATRIX_DENSE_REAL)
    rerror (THIS_SOLVER ": rhs function must return MATRIX-DENSE-REAL");

  retm = ent_data(rent);

  if (SIZE(retm)!= nlbc)
    rerror (THIS_SOLVER ": " RLAB_ERROR_BV_FUNC_INCORRECT_DIM);
  for (j=0; j<nlbc; j++)
    bc[j] = MdrV0 (retm, j);

  ent_Clean (rent);

  //
  // right boundary
  //
  MDPTR(z) = (void *) newzb;

  if (ent_param)
    rent = ent_call_rlab_script_2args(gsubr_fname, ent_z, ent_param);
  else
    rent = ent_call_rlab_script_1arg (gsubr_fname, ent_z);

  if (ent_type(rent)!=MATRIX_DENSE_REAL)
    rerror (THIS_SOLVER ": rhs function must return MATRIX-DENSE-REAL");

  retm = ent_data(rent);

  if (SIZE(retm)!= ncomp-nlbc)
    rerror (THIS_SOLVER ": " RLAB_ERROR_BV_FUNC_INCORRECT_DIM);
  for (j=0; j<nlbc; j++)
    bc[nlbc+j] = MdrV0 (retm, j);

  ent_Clean (rent);
  return (1);
}

//
// differential of the boundary value function for the mirkdc solver
//
int
dgsub_ (int *dummy, double *newza, double *newzb, double *dbc)
{
  int i,j;
  Ent *rent = 0;
  MDR *retm = 0;

  //
  // left boundary
  //
  MDPTR(z) = (void *) newza;

  if (ent_param)
    rent = ent_call_rlab_script_2args(dgsubl_fname, ent_z, ent_param);
  else
    rent = ent_call_rlab_script_1arg (dgsubl_fname, ent_z);

  if (ent_type(rent)!=MATRIX_DENSE_REAL)
    rerror (THIS_SOLVER ": rhs function must return MATRIX-DENSE-REAL");

  retm = ent_data(rent);

  if (SIZE(retm)!= nlbc*ncomp)
    rerror (THIS_SOLVER ": " RLAB_ERROR_BV_FUNC_INCORRECT_DIM);

  for (i=0; i < nlbc; i++)
    for (j = 0; j < ncomp; j++)
      dbc[i + ncomp * j] = Mdr0 (retm, i, j);

  ent_Clean (rent);

  //
  // right boundary
  //
  MDPTR(z) = (void *) newzb;

  if (ent_param)
    rent = ent_call_rlab_script_2args(dgsubr_fname, ent_z, ent_param);
  else
    rent = ent_call_rlab_script_1arg (dgsubr_fname, ent_z);

  if (ent_type(rent)!=MATRIX_DENSE_REAL)
    rerror (THIS_SOLVER ": rhs function must return MATRIX-DENSE-REAL");

  retm = ent_data(rent);

  if (SIZE(retm)!= (ncomp-nlbc)*ncomp)
    rerror (THIS_SOLVER ": " RLAB_ERROR_BV_FUNC_INCORRECT_DIM);

  for (i = nlbc; i < ncomp; i++)
    for (j = 0; j < ncomp; j++)
      dbc[i + ncomp * j] = Mdr0 (retm, i-nlbc, j);

  ent_Clean (rent);
  return (1);
}

//
// User callable functions for coldae.f
//
//
//
//
//
//
//
// The interface to the user-specified function in fortran for coldae in colsys mode
//
int
colfsub (double *newx, double *newz, double *y, double *f)
{
  int j;
  Ent *rent = 0;
  MDR *retm = 0;

  MDPTR(x) = (void *) newx;
  MDPTR(z) = (void *) newz;

  if (ent_param)
    rent = ent_call_rlab_script_3args(fsub_fname, ent_x, ent_z, ent_param);
  else
    rent = ent_call_rlab_script_2args(fsub_fname, ent_x, ent_z);

  if (ent_type(rent)!=MATRIX_DENSE_REAL)
    rerror (THIS_SOLVER ": rhs function must return MATRIX-DENSE-REAL");

  retm = ent_data(rent);
  if (SIZE(retm) != ncomp)
    rerror (THIS_SOLVER ": incorrectly dimensioned derivitive (rhs vector)");

  for (j=0; j<ncomp; j++)
    f[j] = MdrV0 (retm, j);

  ent_Clean (rent);
  return (1);
}


//
// jacobian function for coldae/colsys
//
int
coldfsub (double *newx, double *newz, double *y, double *df)
{
  int i,j;
  Ent *rent = 0;
  MDR *retm = 0;

  MDPTR(x) = (void *) newx;
  MDPTR(z) = (void *) newz;

  if (ent_param)
    rent = ent_call_rlab_script_3args(dfsub_fname, ent_x, ent_z, ent_param);
  else
    rent = ent_call_rlab_script_2args(dfsub_fname, ent_x, ent_z);

  if (ent_type(rent)!=MATRIX_DENSE_REAL)
    rerror (THIS_SOLVER ": rhs function must return MATRIX-DENSE-REAL");

  retm = ent_data(rent);
  if (SIZE(retm) != ncomp*ncomp)
    rerror (THIS_SOLVER ": " RLAB_ERROR_JACOBIAN_FUNC_INCORRECT_DIM);

  for (i = 0; i < ncomp; i++)
    for (j = 0; j < ncomp; j++)
      df[i + ncomp * j] = Mdr0 (retm, i, j);

  ent_Clean (rent);
  return (1);
}

//
// boundary value function for the solver
//
int
colgsub (int *i, double *newz, double *g)
{
  static MDR *retm = 0;

  MDPTR(z) = (void *) newz;

  if (*i == 1)
  {
    Ent *rent = 0;

    //
    // call left boundary condition function once, store tis result in static ent/mdr
    //
    if (ent_param)
      rent = ent_call_rlab_script_2args(gsubl_fname, ent_z, ent_param);
    else
      rent = ent_call_rlab_script_1arg (gsubl_fname, ent_z);

    if (ent_type(rent)!=MATRIX_DENSE_REAL)
      rerror (THIS_SOLVER ": " RLAB_ERROR_RHS_FUNC_MUST_RETURN_MDR);

    retm = mdr_Copy (ent_data(rent));

    if (SIZE(retm) != nlbc)
      rerror (THIS_SOLVER ": " RLAB_ERROR_BV_FUNC_INCORRECT_DIM);

    ent_Clean (rent);
  }
  else if (*i == nlbc+1)
  {
    Ent *rent = 0;

    //
    // call right boundary condition once, store its result in static ent/mdr
    //
    if (ent_param)
      rent = ent_call_rlab_script_2args(gsubr_fname, ent_z, ent_param);
    else
      rent = ent_call_rlab_script_1arg (gsubr_fname, ent_z);

    if (ent_type(rent)!=MATRIX_DENSE_REAL)
      rerror (THIS_SOLVER ": " RLAB_ERROR_RHS_FUNC_MUST_RETURN_MDR);

    retm = mdr_Copy (ent_data(rent));

    if (SIZE(retm) != ncomp-nlbc)
      rerror (THIS_SOLVER ": " RLAB_ERROR_BV_FUNC_INCORRECT_DIM);

    ent_Clean (rent);
  }

  if (*i <= nlbc)
  {
    *g = MdrV1 (retm, *i);
    if (*i == nlbc)
      mdr_Destroy (retm);
  }
  else
  {
    *g = MdrV1 (retm, *i-nlbc);
    if (*i == ncomp)
      mdr_Destroy (retm);
  }

  return (1);
}

//
// differential of the boundary value function for the solver
//
int
coldgsub (int *i, double *newz, double *dg)
{
  int j;
  static MDR *retm = 0;

  MDPTR(z) = (void *) newz;

  if (*i == 1)
  {
    Ent *rent = 0;

    //
    // call left boundary condition function once, store tis result in static ent/mdr
    //
    if (ent_param)
      rent = ent_call_rlab_script_2args(dgsubl_fname, ent_z, ent_param);
    else
      rent = ent_call_rlab_script_1arg (dgsubl_fname, ent_z);

    if (ent_type(rent)!=MATRIX_DENSE_REAL)
      rerror (THIS_SOLVER ": " RLAB_ERROR_RHS_FUNC_MUST_RETURN_MDR);

    retm = mdr_Copy (ent_data(rent));

    if (SIZE(retm) != nlbc*ncomp)
      rerror (THIS_SOLVER ": " RLAB_ERROR_BV_FUNC_INCORRECT_DIM);

    ent_Clean (rent);
  }
  else if (*i == nlbc+1)
  {
    Ent *rent = 0;

    //
    // call right boundary condition once, store its result in static ent/mdr
    //
    if (ent_param)
      rent = ent_call_rlab_script_2args(dgsubr_fname, ent_z, ent_param);
    else
      rent = ent_call_rlab_script_1arg (dgsubr_fname, ent_z);

    if (ent_type(rent)!=MATRIX_DENSE_REAL)
      rerror (THIS_SOLVER ": " RLAB_ERROR_RHS_FUNC_MUST_RETURN_MDR);

    retm = mdr_Copy (ent_data(rent));

    if (SIZE(retm) != (ncomp-nlbc)*ncomp)
      rerror (THIS_SOLVER ": " RLAB_ERROR_BV_FUNC_INCORRECT_DIM);

    ent_Clean (rent);
  }

  if (*i <= nlbc)
  {
    for (j = 0; j < ncomp; j++)
      dg[j] = Mdr1 (retm, *i, j+1);
    if (*i == nlbc)
      mdr_Destroy (retm);
  }
  else
  {
    for (j = 0; j < ncomp; j++)
      dg[j] = Mdr1 (retm, *i-nlbc, j+1);
    if (*i == ncomp)
      mdr_Destroy (retm);
  }

  return (1);
}

int
colguess (double *newx, double *newz, double *newy, double *dmval)
{
  return (1);
}
