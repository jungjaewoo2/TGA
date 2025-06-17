// Copyright (C) 2003-2007 Marijan Kostrun
//   part of rlab+4linux project on rlabplus.sourceforge.net
//
// odae solver
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
//
// See the file ./COPYING
// **********************************************************************


// rlab headers, located in variable $RLAB_SDK
#include "rlab.h"
#include "ent.h"
#include "class.h"
#include "mem.h"
#include "mdr.h"
#include "mdrf1.h"
#include "mds.h"
#include "symbol.h"
#include "list.h"
#include "listnode.h"
#include "list.h"
#include "btree.h"
#include "bltin.h"
#include "util.h"
#include "mathl.h"
#include "function.h"
#include "lp.h"

#include "ode.h"

// gsl headers
#include <gsl/gsl_mode.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_precision.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_machine.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>

// standard libraries
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <time.h>

// naming convention for the solver parameters
#include "rfileio.h"
#include "rlab_solver_parameters_names.h"


//
// mebdfdae functions
//
static int
meb_resid (int * n, double * t, double * y, double * delta, double * yprime,
       int * ipar, double * rpar, int * ierr);

static int
meb_pderv (double * t, double * y, double * pd, int * n, double * yprime,
           int * mbnd4, double * con, int * ipar, double * rpar, int * ierr);
//
// ddaskr functions
//
static int
dda_res (double * t, double * y, double * yprime, double * cj,
         double * delta, int * ires, double * rpar, int * ipar);

static int
dda_jac (double * t, double * y, double * yprime, double *pd, double * cj,
         double * rpar, int * ipar);

//
// bimd functions
//
static int
bimd_func (int * m, double * t, double * y, double * dy,
          int * ierr, double * rpar, int * ipar);
static int
bimd_fjac (int * m, double * t, double * y, double * jac, int * ldjac,
          int * ierr, double * rpar, int * ipar);


// basic global variables
static int    neq;
static MDR  * tmdr;
static MDR  * ymdr;
static MDR  * ypmdr;
static Ent  * tent, * yent, * ypent, * pent;
static Ent  * fname, * dfname;
static MDR  * mass=0;

static int nstep = 0;

#undef THIS_SOLVER
#define THIS_SOLVER "odaei"
Ent *
ent_odaei (int nargs, Datum args[])
{
  int i, i1, j, j1;
  int ii;

  Ent *e3=0, *e4=0, *e5=0, *e6=0, *eo=0, *e7=0;
  double t0, t, t1, ddummy;
  ListNode *node;

  int idummy, louts=0, imethod = 0;
  MDR *ystart, *y, *yprime=0, *yprime0=0, *time=0, *out, *out2=0;
  MDR *rtol=0, *atol=0, *work=0, *iwork = 0, *info;
  char *outs = 0;
  time_t timer;

  FILE *rlab_stderr = (!RLAB_STDERR_DS) ? stderr : RLAB_STDERR_DS;

  //
  // MEBDFI
  //
  double h0 = 1e-6, h;
  int ierr = 0;
  int maxder = 4; // cannot exceed 7
  int maxi = 1000;
  int index1=0, index2=0, index3=0;
  int mf = 22;  // numeric full jacobian, 21 analytic full jacobian
  int mbnd[4]   = {0, 0, 0, 0}; // irrelevant, jacobian is full
//   int masbnd[4] = {0, 0, 0, 0};

  //
  // DDASKR
  //
  int maxord = 5;

  //
  // BIMD
  //
  int imas = 1;   // dae problem
  int mlmas = 0, mumas = 0; // mass matrix is full
  int mljac = 0, mujac = 0; // jacobian is full
  int iout = 0;   // do not call 'solout'
  int ijac = 0;   // assume no jacobian

  //
  // both
  //
  int idid = 1; // first call to the integrator
  int itol = 1;

  FILE *fptr = 0;

  //
  // Load and Check arguments.
  //
  if (nargs < 7 || nargs > 8)
  {
    fprintf (rlab_stderr, THIS_SOLVER ": Ordinary differential algabraic equation solver for initial value\n");
    fprintf (rlab_stderr, THIS_SOLVER ": problem with and without function jacobian.\n");
    fprintf (rlab_stderr, THIS_SOLVER ": Format:\n");
    fprintf (rlab_stderr, THIS_SOLVER ":  y = odaei(f,/df/,/p/,/m/,T,y0/,options/),\n");
    rerror  (THIS_SOLVER ": requires seven or eight arguments!");
  }

  //
  // f = f(t, x/,p/)
  //
  fname = bltin_get_ent(args[0]);
  if (!isfuncent(fname))
    rerror (THIS_SOLVER ": " RLAB_ERROR_ARG1_FUNC_VAR "\n");

  //
  // df = Df/Dx(t, x/,p/)
  //
  dfname = bltin_get_ent(args[1]);
  if (!isfuncent(dfname))
  {
    ent_Clean(dfname);
    dfname = 0;
  }

  if (dfname)
  {
    mf = 21;  // analytic full jacobian for MEBDFI
  }

  // parameter entity
  e3 = bltin_get_ent (args[2]);
  pent = 0;
  if (ent_type (e3) != UNDEF)
    pent = ent_Copy (e3);

  // mass matrix if square,
  mass = 0;
  e4 = bltin_get_ent (args[3]);
  if (ent_type (e4) == MATRIX_DENSE_REAL)
  {
    mass = class_matrix_real (e4);
//     masbnd[0] = 1; // 0 (construct identity matrix), 1 if provided by user
  }
  else
    rerror ("odaei: improper fifth argument 'mass' or 'time'");

  ii = 4;
  if (MNR(mass) != MNC(mass))
  {
    mass = 0;
//     masbnd[0] = 0;
    --ii;
  }

  // t
  if (mass)
  {
    e5 = bltin_get_ent (args[ii]);
    if (ent_type (e5) != MATRIX_DENSE_REAL)
      rerror ("odaei: 'time' has to be a real row-matrix!");
    time = class_matrix_real (e5);
  }
  else
  {
    time = class_matrix_real (e4);
  }

  if (MNR (time) * MNC (time) <= 2)
    rerror ("odaei: 'time' has to be a real matrix!");
  if (MdrV1 (time, 1) == MdrV1 (time, MNC (time) * MNR (time)))
    rerror ("odaei: It cannot be that tstart == tend!");
  t0 = t = MdrV0(time,0);
  ++ii;

  // y0
  e6 = bltin_get_ent (args[ii]);
  if (ent_type(e6)!=MATRIX_DENSE_REAL)
    rerror ("odaei: y0 has to be real vector");
  ystart = class_matrix_real (e6);
  if (MNC (ystart)!=1 && MNR (ystart)!=1)
    rerror ("odaei: y0 has to be real vector");
  ++ii;

  neq = MNR (ystart) * MNC (ystart);
  y = mdr_Create (neq, 1);

  mbnd[3] = neq;
  index1  = neq; // number of index - 1 points

  // yprime0
  e7 = bltin_get_ent (args[ii]);
  if (ent_type(e7) != MATRIX_DENSE_REAL)
    rerror ("odaei: yprime0 has to be real vector");
  yprime0 = class_matrix_real (e7);
  if (MNC (yprime0)*MNR (yprime0) != neq)
    rerror ("odaei: dimension mismatch between y0 and yprime0");
  ++ii;

  yprime = mdr_Float_BF (yprime0);

  // options
  if (nargs > ii)
  {
    eo = bltin_get_ent (args[ii]);
    if (ent_type (eo) == BTREE)
    {
      // imethod = 0 for mebdf, 1 for ddaskr
      node = btree_FindNode (ent_data (eo), RLAB_NAME_ODAE_IMETHOD);
      if (node != 0)
      {
        idummy = (int) class_double (var_ent (node));
        if (idummy == 0 || idummy == 1)
          imethod = idummy;
        else if (idummy == 2 && mass)
          imethod = idummy;
      }
      // stepsize h
      node = btree_FindNode (ent_data (eo), RLAB_NAME_ODAE_STEP);
      if (node != 0)
      {
        ddummy = class_double (var_ent (node));
        if (ddummy > 0.0 || ddummy == -1.0)
          h0 = ddummy;
      }
      // index: default value is 0
      node = btree_FindNode (ent_data (eo), RLAB_NAME_ODAE_INDEX1);
      if (node != 0)
      {
        idummy = (int) class_double (var_ent (node));
        if (idummy >= 1 && idummy <= neq)
          index1 = idummy;
      }
      node = btree_FindNode (ent_data (eo), RLAB_NAME_ODAE_INDEX2);
      if (node != 0)
      {
        idummy = (int) class_double (var_ent (node));
        if (idummy >= 1 && idummy <= neq-index1)
          index2 = idummy;
      }
      node = btree_FindNode (ent_data (eo), RLAB_NAME_ODAE_INDEX3);
      if (node != 0)
      {
        idummy = (int) class_double (var_ent (node));
        if (idummy >= 1 && idummy <= neq-index1-index2)
          index3 = idummy;
      }
      // atol
      node = btree_FindNode (ent_data (eo), RLAB_NAME_ODAE_EABS);
      if (node != 0)
        atol = mdr_Copy(ent_data (var_ent (node)));
        // rtol
      node = btree_FindNode (ent_data (eo), RLAB_NAME_ODAE_EREL);
      if (node != 0)
        rtol = mdr_Copy(ent_data (var_ent (node)));
      // maxi
      node = btree_FindNode (ent_data (eo), RLAB_NAME_ODAE_MAXITER);
      if (node != 0)
      {
        idummy = (int) class_double (var_ent (node));
        if (idummy > 1)
          maxi = idummy;
      }
      // maxder
      node = btree_FindNode (ent_data (eo), RLAB_NAME_ODAE_ORDER);
      if (node != 0)
      {
        idummy = (int) class_double (var_ent (node));
        if (idummy >= 1 && idummy <= 7 && imethod == 0)
          maxder = idummy;
        else if (idummy >= 1 && idummy <= 5 && imethod == 1)
          maxord = idummy;
      }
      // standard output
      node = btree_FindNode (ent_data (eo), RLAB_NAME_ODAE_STDOUT);
      if (node != 0)
      {
        if (ent_type(var_ent (node))== MATRIX_DENSE_STRING)
        {
          outs = class_char_pointer ((Ent *) var_ent (node));
        }
      }
    }
  }

  nstep = MNR (time) * MNC (time);
  out = mdr_Create (nstep, neq + 1);

  if (!atol)
    atol = mdr_CreateScalar(1.0);
  if (!rtol)
    rtol = mdr_CreateScalar(1.0);

  if (imethod == 2)
  {
    //
    // BIMD
    //
    mlmas = neq;
    mumas = neq;
    mljac = neq;
    mujac = neq;
    if (dfname)
      ijac = 1;
    if (MNR(atol)*MNC(atol)==neq)
      itol = 1;
    else
      itol = 0;
  } // end of BIMD
  else if (imethod == 1)
  {
    //
    // ddaskr
    //
    info = mdi_Create(20,1);
    mdr_Zero (info);
    MdiV0(info,  1 - 1) = 0;    // before the first call
    // tolerances
    if (MNR(rtol)*MNC(rtol)==neq && MNR(atol)*MNC(atol)==neq)
      MdiV0(info,  2 - 1) = 1;    // both tolerances are vectors of length neq
    else
      MdiV0(info,  2 - 1) = 0;    // both tolerances are scalars
    MdiV0(info,  3 - 1) = 0;    // interval output mode
    MdiV0(info,  4 - 1) = 0;    // integration may go past t1
    // jacobian
    if (dfname)
      MdiV0(info,  5 - 1) = 1;  // analytic
    else
      MdiV0(info,  5 - 1) = 0;  // via finite differences
    MdiV0(info,  6 - 1) = 0;    // dense matrix structure of jac
    MdiV0(info,  7 - 1) = 1;    // max step size is given by 1/4*DeltaT
    // initial step size
    if (h > 0)
      MdiV0(info,  8 - 1) = 1;  // provided by user
    else
      MdiV0(info,  8 - 1) = 0;  // let the code decide
    MdiV0(info,  9 - 1) = 1;    // user supplied maxord
    MdiV0(info, 10 - 1) = 0;    // no special constraints on solution
    MdiV0(info, 11 - 1) = 0;    // initial values for y,y' are consistent
    MdiV0(info, 12 - 1) = 0;    // use Newton iteration with direct method
    MdiV0(info, 16 - 1) = 0;    // error control on all variables !!!!!!
  } // end of DDASKR
  else
  {
    //
    // MEBDFI: mass matrix and tolerances
    //
    // tolerances
    if (MNR(rtol)*MNC(rtol)==1 && MNR(atol)*MNC(atol)==1)
    {
      if (MdrV0(atol,0)==MdrV0(rtol,0))
        itol = 1;
      else
        itol = 2;
    }
    else if (MNR(rtol)*MNC(rtol)==1 && MNR(atol)*MNC(atol)==neq)
      itol = 3;
    else if (MNR(rtol)*MNC(rtol)==neq && MNR(atol)*MNC(atol)==1)
      itol = 4;
    else if (MNR(rtol)*MNC(rtol)==neq && MNR(atol)*MNC(atol)==neq)
      itol = 5;
    else
    {
      itol = 1;
      MdrV0(atol,0) = MdrV0(rtol,0) = 1.0;
    }
  } // end of MEBDFI

  //
  // t
  //
  tmdr = mdr_CreateEmpty (1, 1);
  tent = ent_Assign_Rlab_MDR (tmdr);
  ent_IncRef (tent);

  //
  // y
  //
  ymdr = mdr_Create (0, 0);
  ymdr->nrow = neq;
  ymdr->ncol = 1;
  yent = ent_Create ();
  ent_data (yent) = ymdr;
  ent_SetType (yent, MATRIX_DENSE_REAL);
  ent_IncRef (yent);

  //
  // if mass is given
  //    f = f(t,y/,p/), and the problem is mass yp = f(t,y/,p/)
  // else
  //    f = f(t,y,yp/,p/), and the problem is    f(t,y,yp/,p/) = 0
  //
  if(!mass)
  {
    //
    // yprime
    //
    ypmdr = mdr_CreateEmpty (neq, 1);
    ypent = ent_Assign_Rlab_MDR (ypmdr);
    ent_IncRef (ypent);
  }

  // Save initial conditions and setup y[]
  Mdr1 (out, 1, 1) = Mdr1 (time, 1, 1);
  for (j = 2; j <= neq + 1; j++)
  {
    Mdr1 (out, 1, j) = MdrV1 (ystart, j - 1);
    MdrV1 (y, j - 1) = MdrV1 (ystart, j - 1);
  }

  if (outs)
    louts = strlen (outs);
  if (louts > 0)
    fptr = fopen (outs, "a");

  if (fptr)
  {
    timer = clock ();
    fprintf (fptr, "RLaB: using built-in integrator ");
    if (imethod == 2)
      fprintf (fptr, "'BIMD'");
    else if (imethod == 1)
      fprintf (fptr, "'DDASKR'");
    else
      fprintf (fptr, "'MEBDFI'");
    fprintf (fptr, " for inital value DAE problem.\n"
                   "RLaB: the messages from the solver follow\n"
              );
    }

  if (imethod == 2)
  {
    //
    // bimd
    //
    //work and iwork arrays
    int lwork = 14 + (maxord-2) + 9 * neq + 5 * (maxord-2) * neq + 3*neq*neq;
    work  = mdr_Create(lwork,1);
    int liwork = 40 + neq;
    iwork  = mdi_Create(liwork,1);
    double h = h0;

    // work
    MdrV1(work,1)  = 1e-16;       // machine precision
    MdrV1(work,3)  = 1e-1;        // factol, method of order 4
    MdrV1(work,4)  = 1e-1;        // factol, method of order 6
    MdrV1(work,5)  = 1e-1;        // factol, method of order 8
    MdrV1(work,6)  = 1e-1;        // factol, method of order 10
    MdrV1(work,7)  = 1e-1;        // factol, method of order 12
    MdrV1(work,8)  = 1e-2;        // factol, stopping criterion
    MdrV1(work,9)  = 5e-2;        // factol, stopping criterion
    MdrV1(work,10) = 1.2e-1;      // facl,
    MdrV1(work,11) = 1e1;         // facr,
    MdrV1(work,12) = 5e-2;        // sfty,
    MdrV1(work,13) = 5e-2/2.0;    // sftyup,
    MdrV1(work,14) = 5e-2;        // sftydn,

    // iwork
    MdiV1(iwork,1) = maxi;      // max no. of iterations
    MdiV1(iwork,2) = 4;         // minimum order
    MdiV1(iwork,3) = maxord;    // maximum order 4,6,8,10,12
    MdiV1(iwork,4) = 10;        // max no. of steps  for method of order 4
    MdiV1(iwork,5) = 10;        // max no. of steps  for method of order 6
    MdiV1(iwork,6) = 10;        // max no. of steps  for method of order 8
    MdiV1(iwork,7) = 10;        // max no. of steps  for method of order 10
    MdiV1(iwork,8) = 10;        // max no. of steps  for method of order 12
    MdiV1(iwork,9)  = index1;   // no of index 1 variables
    MdiV1(iwork,10) = index2;   // no of index 2 variables
    MdiV1(iwork,11) = index3;   // no of index 3 variables

    idid = 0;
    for (i = 2; (i <= nstep) && !idid; i++)
    {
      t0 = t;
      t  = MdrV1(time,i);

      // max step size is 1/4 of the current integration interval
      MdrV1(work,2) = 0.25 * (t - t0);

      BIMD(&neq, bimd_func, &t0, &t, MDRPTR(y), &h,
            MDRPTR(rtol), MDRPTR(atol), &itol, bimd_fjac, &ijac, &mljac, &mujac,
            MDRPTR(mass), &imas, &mlmas, &mumas, NULL, &iout,
            MDRPTR(work), &lwork, MDIPTR(iwork), &liwork,
            NULL, NULL, &idid);

      Mdr1 (out, i, 1) = t0;
      for (j = 2; j <= neq + 1; j++)
        Mdr1 (out, i, j) = MdrV1 (y, j - 1);

      //
      // print progress of integration
      //
      if (fptr)
      {
        fprintf (fptr, "t = %g ", Mdr1 (out, i, 1));
        for (j = 2; j <= neq + 1; j++)
          fprintf (fptr, "  x[%i] = %14.6g", j - 1, Mdr1 (out, i, j));
        fprintf (fptr, "\n");
      }
    }
    if (fptr)
    {
      if (idid<0)
      {
        fprintf (fptr, "RLaB: odaei: Solver 'BIMD' reports the error code %i (",
                 idid);
        if (idid == -2)
          fprintf (fptr, "EXCEEDED NUMBER OF INTEGRATION STEPS).\n");
        else if (idid == -3)
          fprintf (fptr, "STEPSIZE TOO SMALL.\n");
        else if (idid == -4)
          fprintf (fptr, "THE MATRIX OF PARTIAL DERIVATIVES IS SINGULAR).\n");
        else if (idid == -5)
          fprintf (fptr, "NONLINEAR SYSTEM SOLVER FAILED TO CONVERGE).\n");
        else if (idid == -6)
          fprintf (fptr, "ERROR IN EVALUATION OF FUNCTION OR ITS JACOBIAN).\n");
        else
          fprintf (fptr, "RLaB: odaei: UNKNOWN ERROR CODE.\n");
        fprintf (fptr, "RLaB: odaei: Interrupting the integration!\n");
      }
      else
        fprintf (fptr, "RLaB: odaei: Solver 'BIMD' reports success.\n");
    }

    //
    // was calculation successful?
    //
    if (idid<0)
    {
      out2 = mdr_Create(i-1, neq+1);
      for (i1 = 0; i1 < i-1 ; i1++)
        for (j1 = 0; j1 < neq + 1; j1++)
          Mdr0(out2,i1,j1) = Mdr0(out,i1,j1);
      mdr_Destroy (out);
      out = 0;
    }

  } // end of bimd
  else if (imethod == 1)
  {
    //
    // ddaskr
    //
    int nrt = 0;

    if (fptr)
      MdiV0(info, 18 - 1) = 2;  // be much verbose

    //work and iwork arrays
    int lwork = 60 + MAX(maxord+4,7) * neq + neq*neq;
    work  = mdr_Create(lwork,1);
    int liwork = 40 + neq;
    iwork  = mdi_Create(liwork,1);

    // user supplied maxord
    MdiV0(iwork,2) = maxord;
    // user provided initial step size
    if (h0 > 0)
      MdrV0(work, 2) = h0;

    idid = 2;
    for (i = 2; (i <= nstep) && (idid==2 || idid==3); i++)
    {
      t0 = t;
      t  = MdrV1(time,i);
      if (i == nstep)
      {
        // last integration has to end at t
        MdiV0(info, 4 - 1) = 1;
        MdrV0(work, 0) = t;
      }
      MdrV0(work, 1) = 0.25 * (t - t0); // max step size is 1/4 of the current
                                        // integration interval

      DDASKR (dda_res, &neq, &t0, MDRPTR(y), MDRPTR(yprime), &t, MDIPTR(info), MDRPTR(rtol), MDRPTR(atol),
              &idid, MDRPTR(work), &lwork, MDIPTR(iwork), &liwork, NULL, NULL, dda_jac, NULL,
              NULL, &nrt, NULL);

      Mdr1 (out, i, 1) = t0;
      for (j = 2; j <= neq + 1; j++)
        Mdr1 (out, i, j) = MdrV1 (y, j - 1);

      //
      // print progress of integration
      //
      if (fptr)
      {
        fprintf (fptr, "t = %g ", Mdr1 (out, i, 1));
        for (j = 2; j <= neq + 1; j++)
          fprintf (fptr, "  x[%i] = %14.6g", j - 1, Mdr1 (out, i, j));
        fprintf (fptr, "\n");
      }
    }
    if (fptr)
    {
      if (idid<0)
      {
        fprintf (fptr, "RLaB: odaei: Solver 'DDASKR' reports the error code %i (",  idid);
        if (idid == -1)
          fprintf (fptr, "EXCEEDED NUMBER OF INTEGRATION STEPS).\n");
        else if (idid == -2)
          fprintf (fptr, "ACCURACY TOO HIGH).\n");
        else if (idid == -3)
          fprintf (fptr, "ERROR IN TOLERANCES - FAILING TO ACHIEVE CORRECTOR CONVERGENCE).\n");
        else if (idid == -5)
          fprintf (fptr, "ILLEGAL VALUES IN JACOBIAN).\n");
        else if (idid == -6)
          fprintf (fptr, "FAILED ERROR TESTS).\n");
        else if (idid == -7)
          fprintf (fptr, "NONLINEAR SYSTEM SOLVER FAILED TO CONVERGE).\n");
        else if (idid == -8)
          fprintf (fptr, "THE MATRIX OF PARTIAL DERIVATIVES IS SINGULAR).\n");
        else
          fprintf (fptr, "RLaB: odaei: UNKNOWN ERROR CODE.\n");
        fprintf (fptr, "RLaB: odaei: Interrupting the integration!\n");
      }
      else
        fprintf (fptr, "RLaB: odaei: Solver 'DDASKR' reports success.\n");
    }

    //
    // was calculation successful?
    //
    if (idid!=2 && idid!=3)
    {
      out2 = mdr_Create(i-1, neq+1);
      for (i1 = 0; i1 < i-1 ; i1++)
        for (j1 = 0; j1 < neq + 1; j1++)
          Mdr0(out2,i1,j1) = Mdr0(out,i1,j1);
      mdr_Destroy (out);
      out = 0;
    }

    mdr_Destroy (info);
  }
  else
  {
    //
    // mebdfi
    //
    int lwork  = 4 + 32 * neq + 3 * neq * neq;
    work = mdr_Create (lwork, 1);

    int liwork = neq + 14;
    iwork  = mdi_Create (liwork, 1);
    MdiV0(iwork,0)  = index1;   // no. of index 1 variables
    MdiV0(iwork,1)  = index2;   // no. of index 2 variables
    MdiV0(iwork,2)  = index3;   // no. of index 3 variables
    MdiV0(iwork,13) = maxi;     // max. no. of steps

    idid = 1;
    t0 = MdrV1 (time, 1);
    h  = MdrV1 (work, 2) = h0;        // start with user supplied stepsize
    for (i = 2; (i <= nstep) && (!ierr && idid >=0); i++)
    {
      t1 = t = MdrV1 (time, i);      // desired step
      MEBDFI (&neq, &t0, &h, MDRPTR(y), MDRPTR(yprime), &t, &t1, &mf, &idid, &lwork,
              MDRPTR(work), &liwork, MDIPTR(iwork), mbnd, &maxder,
              &itol, MDRPTR(rtol), MDRPTR(atol), NULL, NULL,
              meb_pderv, meb_resid, &ierr,
              outs, &louts);
      if (idid < 0)
        break;
      else
        idid = 2; // finish integration at t1 always

      Mdr1 (out, i, 1) = t;
      for (j = 2; j <= neq + 1; j++)
        Mdr1 (out, i, j) = MdrV1 (y, j - 1);

      //
      // print progress of integration
      //
      if (fptr)
      {
        fprintf (fptr, "t = %g ", Mdr1 (out, i, 1));
        for (j = 2; j <= neq + 1; j++)
          fprintf (fptr, "  x[%i] = %14.6g", j - 1, Mdr1 (out, i, j));
        fprintf (fptr, "\n");
      }
    }

    if (fptr)
    {
      if (ierr)
      {
        fprintf (fptr, "RLaB: odaei: Solver 'MEBDF' reports the error code %i (",  idid);
        if (idid == -1)
          fprintf (fptr, "STEPSIZE TOO SMALL).\n");
        else if (idid == -2)
          fprintf (fptr, "ACCURACY TOO HIGH).\n");
        else if (idid == -3)
          fprintf (fptr, "FAILING TO ACHIEVE CORRECTOR CONVERGENCE).\n");
        else if (idid == -4)
          fprintf (fptr, "ILLEGAL VALUES OF INPUT PARAMETERS).\n");
        else if (idid == -5)
          fprintf (fptr, "INTERNAL ERROR).\n");
        else if (idid == -6)
          fprintf (fptr, "EXCEEDED NUMBER OF INTEGRATION STEPS).\n");
        else
          fprintf (fptr, "RLaB: odeiv: UNKNOWN ERROR CODE.\n");
      }
      else
        fprintf (fptr, "RLaB: odaei: Solver 'MEBDF' reports success!\n");
    }

    //
    // was calculation successful?
    //
    if (ierr)
    {
      out2 = mdr_Create(i-1, neq+1);
      for (i1 = 0; i1 < i-1 ; i1++)
        for (j1 = 0; j1 < neq + 1; j1++)
          Mdr0(out2,i1,j1) = Mdr0(out,i1,j1);
      mdr_Destroy (out);
      out = 0;
    }
  }

  if (fptr)
  {
     // check the time and close the output stream
    timer -= clock ();
    fprintf (fptr,
             "RLaB: odaei: Integration lasted %g sec.\n",
             -timer / 1e6);
    fclose (fptr);
  }

  // clean-up
  if (work)
    mdr_Destroy (work);
  if (iwork)
    mdr_Destroy (iwork);

  if (atol)
    mdr_Destroy(atol);
  if (rtol)
    mdr_Destroy(rtol);

  // cleanup the arguments of functions f, df
  MDPTR(tmdr) = 0;
  ent_DecRef  (tent);
  ent_Destroy (tent);

  ent_Clean (e3); ent_Clean (pent);

  MDPTR(ymdr) = 0;
  ent_DecRef (yent);
  ent_Destroy (yent);

  if(!mass)
  {
    MDPTR(ypmdr) = 0;
    ent_DecRef (ypent);
    ent_Destroy (ypent);
  }

  //
  // Clean Up
  //
  mdr_Destroy (y);
  mdr_Destroy (yprime);

  // cleanup the rest
  ent_Clean (e4);
  ent_Clean (e5);
  ent_Clean (e6);
  ent_Clean (e7);
  ent_Clean (eo);
  if (mass)
    mass = 0;
  ent_Clean (fname);
  ent_Clean (dfname);

  if (out)
    return ent_Assign_Rlab_MDR(out);
  else
    return ent_Assign_Rlab_MDR(out2);
}

//
//
//
static int
meb_resid (int * n, double * t, double * y, double * delta, double * yprime,
           int * ipar, double * rpar, int * ierr)
{
  int i, j;
  Ent *rent = 0;
  MDR *retm = 0;

  MDPTR (tmdr) = t;
  MDPTR (ymdr) = (void *) y;

  if (mass)
  {
    // check 'p' and call f(t,y/,p/)
    if (pent)
      rent = ent_call_rlab_script_3args(fname, tent, yent, pent);
    else
      rent = ent_call_rlab_script_2args(fname, tent, yent);
  }
  else
  {
    MDPTR (ypmdr) = (void *) yprime;

    // check 'p' and call f(t,y,yp/,p/)
    if (pent)
      rent = ent_call_rlab_script_4args(fname, tent, yent, ypent, pent);
    else
      rent = ent_call_rlab_script_3args(fname, tent, yent, ypent);
  }

  if (ent_type(rent)!=MATRIX_DENSE_REAL)
    rerror (THIS_SOLVER ": " RLAB_ERROR_RHS_FUNC_MUST_RETURN_MDR);

  retm = ent_data (rent);

  if (SIZE(retm) != neq)
    rerror ("odaei: incorrectly dimensioned function (rhs vector)");

  if (mass)
  {
    //
    // given f = f(t,y/,p/) and y'=yprime, find delta,
    // where
    //    delta = mas*y' - f
    //
    for (i = 0; i < neq; i++)
    {
      delta[i] = -MdrV0 (retm, i);
      for (j = 0; j < neq; j++)
        delta[i] += Mdr0(mass,i,j) * yprime[j];
    }
  }
  else
  {
    //
    // here:
    //    delta = f(t,y,yp/,p/)
    //
    for (i = 0; i < neq; i++)
      delta[i] = MdrV0 (retm, i);
  }

  // Try to clean up if possible.
  ent_Clean (rent);
  return 1;
}


static int
meb_pderv (double * t, double * y, double * pd, int * n, double * yprime,
           int * mbnd4, double * con, int * ipar, double * rpar, int * ierr)
{
  int i, j;
  Ent *rent = 0;
  MDR *retm = 0;
  double rcon = 1.0/(*con);

  MDPTR (tmdr) = t;
  MDPTR (ymdr) = (void *) y;

  if (mass)
  {
    // check 'p' and call f(t,y/,p/)
    if (pent)
      rent = ent_call_rlab_script_3args(dfname, tent, yent, pent);
    else
      rent = ent_call_rlab_script_2args(dfname, tent, yent);
  }
  else
  {
    MDPTR (ypmdr) = (void *) yprime;

    // check 'p' and call f(t,y,yp/,p/)
    if (pent)
      rent = ent_call_rlab_script_4args(dfname, tent, yent, ypent, pent);
    else
      rent = ent_call_rlab_script_3args(dfname, tent, yent, ypent);
  }

  if (ent_type(rent)!=MATRIX_DENSE_REAL)
    rerror (THIS_SOLVER ": " RLAB_ERROR_JACOBIAN_FUNC_MUST_RETURN_MDR);

  retm = ent_data (rent);

  if (SIZE (retm) != ((1 + (mass==0)) * neq * neq))
    rerror (THIS_SOLVER ": " RLAB_ERROR_JACOBIAN_FUNC_INCORRECT_DIM);

  if (mass)
  {
    // we have
    //    jac[i;j] = dG(i)/dY(j) + CJ*dG(i)/dYPRIME(j),
    // where
    //    G = mass*y' - f(t,y)
    // so
    //    jac[i;j] = -df(i)/dy(j) + cj * mass(i,j)
    //
    for (i = 0; i < neq; i++)
      for (j = 0; j < neq; j++)
        pd[j*neq + i] = -Mdr0 (retm, i, j) +  rcon * Mdr0(mass,i,j);
  }
  else
  {
    // we have
    //    pd[i;j] = dG(i)/dY(j) + (1/con) * dG(i)/dYPRIME(j),
    // where
    //    dG(i)/dY(j) = retm[i;j]
    //    dG(i)/dYPRIME(j) = retm[i;neq+j]
    // i.e., G is stacked then
    for (i = 0; i < neq; i++)
      for (j = 0; j < neq; j++)
        pd[j*neq + i] = Mdr0 (retm, i, j) +  rcon * Mdr0(retm,i,neq+j);
  }

  ent_Clean (rent);
  return 1;
}



//
// ddaskr: func
//
static int
dda_res (double * t, double * y, double * yprime, double * cj,
         double * delta, int * ires, double * rpar, int * ipar)
{
  int i, j;
  Ent *rent = 0;
  MDR *retm = 0;

  MDPTR (tmdr) = t;
  MDPTR (ymdr) = (void *) y;

  if (mass)
  {
    // check 'p' and call f(t,y/,p/)
    if (pent)
      rent = ent_call_rlab_script_3args(fname, tent, yent, pent);
    else
      rent = ent_call_rlab_script_2args(fname, tent, yent);
  }
  else
  {
    MDPTR (ypmdr) = (void *) yprime;

    // check 'p' and call f(t,y,yp/,p/)
    if (pent)
      rent = ent_call_rlab_script_4args(fname, tent, yent, ypent, pent);
    else
      rent = ent_call_rlab_script_3args(fname, tent, yent, ypent);
  }

  if (ent_type(rent)!=MATRIX_DENSE_REAL)
    rerror (THIS_SOLVER ": " RLAB_ERROR_RHS_FUNC_MUST_RETURN_MDR);

  retm = ent_data (rent);

  if (SIZE(retm) != neq)
    rerror (THIS_SOLVER ": " RLAB_ERROR_RHS_FUNC_INCORRECT_DIM);

  if (mass)
  {
    //
    // for f = f(t,y,yp/,p/),
    //    delta = mas*yp - f
    //
    for (i = 0; i < neq; i++)
    {
      delta[i] = -MdrV0 (retm, i);
      for (j = 0; j < neq; j++)
        delta[i] += Mdr0(mass,i,j) * yprime[j];
    }
  }
  else
  {
    // delta = f(t,y,yp/,p/)
    for (i = 0; i < neq; i++)
      delta[i] = MdrV0 (retm, i);
  }

  ent_Clean (rent);
  return 1;
}

static int
dda_jac (double * t, double * y, double * yprime, double *jac, double * cj,
         double * rpar, int * ipar)
{
  int i, j;
  Ent *rent = 0;
  MDR *retm = 0;

  MDPTR (tmdr) = t;
  MDPTR (ymdr) = (void *) y;

  if (mass)
  {
    // check 'p' and call f(t,y/,p/)
    if (pent)
      rent = ent_call_rlab_script_3args(dfname, tent, yent, pent);
    else
      rent = ent_call_rlab_script_2args(dfname, tent, yent);
  }
  else
  {
    MDPTR (ypmdr) = (void *) yprime;

    // check 'p' and call f(t,y,yp/,p/)
    if (pent)
      rent = ent_call_rlab_script_4args(dfname, tent, yent, ypent, pent);
    else
      rent = ent_call_rlab_script_3args(dfname, tent, yent, ypent);
  }

  if (ent_type(rent)!=MATRIX_DENSE_REAL)
    rerror (THIS_SOLVER ": " RLAB_ERROR_JACOBIAN_FUNC_MUST_RETURN_MDR);

  retm = ent_data (rent);

  if (SIZE(retm) != ((1 + (mass==0)) * neq * neq))
    rerror (THIS_SOLVER ": " RLAB_ERROR_JACOBIAN_FUNC_INCORRECT_DIM);

  if (mass)
  {
    // we have
    //    jac[i;j] = dG(i)/dY(j) + CJ*dG(i)/dYPRIME(j),
    // where
    //    G = mass*y' - f(t,y)
    // so
    //    jac[i;j] = -df(i)/dy(j) + cj * mass(i,j)
    //
    for (i = 0; i < neq; i++)
      for (j = 0; j < neq; j++)
        jac[j*neq + i] = -Mdr0 (retm, i, j) + (*cj) * Mdr0(mass,i,j);
  }
  else
  {
    // we have
    //    pd[i;j] = dG(i)/dY(j) + (cj) * dG(i)/dYPRIME(j),
    // where
    //    dG(i)/dY(j) = retm[i;j]
    //    dG(i)/dYPRIME(j) = retm[i;neq+j]
    // i.e., G is stacked then
    for (i = 0; i < neq; i++)
      for (j = 0; j < neq; j++)
        jac[j*neq + i] = Mdr0 (retm, i, j) + (*cj) * Mdr0(retm,i,j+neq);
  }

  ent_Clean (rent);
  return 1;
}

//
// bimd functions
//
static int
bimd_func (int * m, double * t, double * y, double * dy,
           int * ierr, double * rpar, int * ipar)
{
  int i;
  Ent *rent = 0;
  MDR *retm = 0;

  MDPTR (tmdr) = t;
  MDPTR (ymdr) = (void *) y;

  if (pent)
    rent = ent_call_rlab_script_3args(fname, tent, yent, pent);
  else
    rent = ent_call_rlab_script_2args(fname, tent, yent);

  if (ent_type (rent) != MATRIX_DENSE_REAL)
    rerror (THIS_SOLVER ": " RLAB_ERROR_JACOBIAN_FUNC_MUST_RETURN_MDR);

  retm = ent_data (rent);

  if (SIZE(retm) != neq)
    rerror (THIS_SOLVER ": " RLAB_ERROR_RHS_FUNC_MUST_RETURN_MDR);

  for (i = 0; i < neq; i++)
    dy[i] = MdrV0 (retm, i);

  ent_Clean (rent);
  return 1;
}

static int
bimd_fjac (int * m, double * t, double * y, double * jac, int * ldjac, int * ierr,
          double * rpar, int * ipar)
{
  int i, j;
  Ent *rent = 0;
  MDR *retm = 0;

  MDPTR (tmdr) = t;
  MDPTR (ymdr) = (void *) y;

  if (pent)
    rent = ent_call_rlab_script_3args(dfname, tent, yent, pent);
  else
    rent = ent_call_rlab_script_2args(dfname, tent, yent);

  if (ent_type (rent) != MATRIX_DENSE_REAL)
    rerror (THIS_SOLVER ": " RLAB_ERROR_JACOBIAN_FUNC_MUST_RETURN_MDR);

  retm = ent_data (rent);

  if (SIZE(retm) != neq * neq)
    rerror (THIS_SOLVER ": " RLAB_ERROR_JACOBIAN_FUNC_INCORRECT_DIM);

  for (i = 0; i < neq; i++)
    for (j = 0; j < neq; j++)
      jac[j*neq + i] = Mdr0 (retm, i, j);

  *ierr = 0;
  ent_Clean (rent);
  return 1;
}

