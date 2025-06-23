// Copyright (C) 2003-2007 Marijan Kostrun
//   part of rlab+4linux project on rlabplus.sourceforge.net
//
// GSL Science Library - ode solver
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
#include <gsl/gsl_version.h>
#include <gsl/gsl_mode.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_precision.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_machine.h>
#include <gsl/gsl_matrix.h>

//
// i like when they switch headers between releases....
//
#if (GSL_MAJOR_VERSION <= 1 && GSL_MINOR_VERSION < 16)

#include <gsl/gsl_odeiv.h>
#define GSL_ODEIV_STEP_TYPE   gsl_odeiv_step_type
#define GSL_ODEIV_STEP_RK2    gsl_odeiv_step_rk2
#define GSL_ODEIV_STEP_RK4    gsl_odeiv_step_rk4
#define GSL_ODEIV_STEP_RKF45  gsl_odeiv_step_rkf45
#define GSL_ODEIV_STEP_RKCK   gsl_odeiv_step_rkck
#define GSL_ODEIV_STEP_RK8PD  gsl_odeiv_step_rk8pd
#define GSL_ODEIV_STEP_BSIMP  gsl_odeiv_step_bsimp
#define GSL_ODEIV_STEP        gsl_odeiv_step
#define GSL_ODEIV_STEP_ALLOC  gsl_odeiv_step_alloc
#define GSL_ODEIV_CONTROL     gsl_odeiv_control
#define GSL_ODEIV_CONTROL_YP_NEW        gsl_odeiv_control_yp_new
#define GSL_ODEIV_CONTROL_Y_NEW         gsl_odeiv_control_y_new
#define GSL_ODEIV_CONTROL_STANDARD_NEW  gsl_odeiv_control_standard_new
#define GSL_ODEIV_EVOLVE        gsl_odeiv_evolve
#define GSL_ODEIV_EVOLVE_ALLOC  gsl_odeiv_evolve_alloc
#define GSL_ODEIV_SYSTEM        gsl_odeiv_system 
#define GSL_ODEIV_CONTROL_NAME  gsl_odeiv_control_name
#define GSL_ODEIV_EVOLVE_APPLY  gsl_odeiv_evolve_apply
#define GSL_ODEIV_EVOLVE_RESET  gsl_odeiv_evolve_reset
#define GSL_ODEIV_EVOLVE_FREE   gsl_odeiv_evolve_free
#define GSL_ODEIV_CONTROL_FREE  gsl_odeiv_control_free
#define GSL_ODEIV_STEP_FREE     gsl_odeiv_step_free

#else

#include <gsl/gsl_odeiv2.h>
#define GSL_ODEIV_STEP_TYPE   gsl_odeiv2_step_type
#define GSL_ODEIV_STEP_RK2    gsl_odeiv2_step_rk2
#define GSL_ODEIV_STEP_RK4    gsl_odeiv2_step_rk4
#define GSL_ODEIV_STEP_RKF45  gsl_odeiv2_step_rkf45
#define GSL_ODEIV_STEP_RKCK   gsl_odeiv2_step_rkck
#define GSL_ODEIV_STEP_RK8PD  gsl_odeiv2_step_rk8pd
#define GSL_ODEIV_STEP_BSIMP  gsl_odeiv2_step_bsimp
#define GSL_ODEIV_STEP        gsl_odeiv2_step
#define GSL_ODEIV_STEP_ALLOC  gsl_odeiv2_step_alloc
#define GSL_ODEIV_CONTROL     gsl_odeiv2_control
#define GSL_ODEIV_CONTROL_YP_NEW        gsl_odeiv2_control_yp_new
#define GSL_ODEIV_CONTROL_Y_NEW         gsl_odeiv2_control_y_new
#define GSL_ODEIV_CONTROL_STANDARD_NEW  gsl_odeiv2_control_standard_new
#define GSL_ODEIV_EVOLVE        gsl_odeiv2_evolve
#define GSL_ODEIV_EVOLVE_ALLOC  gsl_odeiv2_evolve_alloc
#define GSL_ODEIV_SYSTEM        gsl_odeiv2_system 
#define GSL_ODEIV_CONTROL_NAME  gsl_odeiv2_control_name
#define GSL_ODEIV_EVOLVE_APPLY  gsl_odeiv2_evolve_apply
#define GSL_ODEIV_EVOLVE_RESET  gsl_odeiv2_evolve_reset
#define GSL_ODEIV_EVOLVE_FREE   gsl_odeiv2_evolve_free
#define GSL_ODEIV_CONTROL_FREE  gsl_odeiv2_control_free
#define GSL_ODEIV_STEP_FREE     gsl_odeiv2_step_free

#endif


// standard libraries
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <time.h>

// rlabplus
#include "rfileio.h"
#include "rlab_solver_parameters_names.h"
#include "rlab_macros.h"
#include "rlab_macros_code.h"

//
// the gsl functions
//
static int
ode_gsl_func (double t, const double y[], double yp[], void *dummy);
static int
ode_gsl_fjac (double t, const double y[], double *dfdy,
                            double *dfdt, void *dummy);

//
// ODE.F function
//
static int
ode_func (double *t, double *y, double *yp);

//
// odebim functions
//
static int
bim_func (int * m, double * t, double * y, double * dy,
          int * ierr, double * rpar, int * ipar);
static int
bim_fjac (int * m, double * t, double * y, double * jac, int * ldjac,
          int * ierr, double * rpar, int * ipar);

//
// vode functions
//
static int
vode_func (int * neq, double * t, double * y, double * ydot,
           double * rpar, int * ipar);
static int
vode_fjac (int * neq, double * t, double * y,
           int * ml, int * mu, double * pd,
           int * nrpd, double * rpar, int * ipar);

// basic global variables
static int   neq;
static MDR  *tmdr=0;
static MDR  *ymdr=0;
static Ent  *tent=0,  *yent=0, *pent=0;
static Ent  *fname_ent=0, *jac_fname_ent=0;

static int nstep = 0;


// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                              +
// + ODE integration library, one that started the project        +
// + rlabplus                                                     +
// +                                                              +
// + by Marijan Kostrun, V-2005, 2006, 2007                       +
// +                                                              +
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#undef THIS_SOLVER
#define THIS_SOLVER "odeiv"
Ent *
ent_gsl_odeiv (int nargs, Datum args[])
{
  int i, i1, j, j1;
  Ent *e2=0, *e3=0, *e4=0, *e5=0, *eo=0;
  double t1, t2, t2n;
  ListNode *node;

  int status=0;
  MDR *ystart=0, *y=0, *yerr=0, *time=0, *ysol=0, *ysol2 = 0, *yp=0;
  char * outs = 0;
  char odename[24];
  time_t timer=0;

  int ips = 0;  // return the phase space trajectory

  //
  // the gsl
  //
  double h0=1e-6, h, ay = 1, adydt = 0;

  //
  // odebim
  //
  int bim_ordmin = 4;
  int bim_ordmax = 12;

  //
  // both
  //
  double abserr = 1e-6, relerr = 1e-6, teps = 1e-9;
  int imethod = 0, notnan = 1;

  FILE *fptr = NULL;
  FILE *rlab_stderr = (!RLAB_STDERR_DS) ? stderr : RLAB_STDERR_DS;

  //
  // Load and Check arguments.
  //
  if (nargs < 4 || nargs > 6)
  {
    fprintf (rlab_stderr, THIS_SOLVER ": Ordinary differential equation solver for an initial value problem\n");
    fprintf (rlab_stderr, THIS_SOLVER ": with and without jacobian.\n");
    fprintf (rlab_stderr, THIS_SOLVER ": Format:\n");
    fprintf (rlab_stderr, THIS_SOLVER ": (1) y = odeiv(f,/p/,T,y0/,options/),\n");
    fprintf (rlab_stderr, THIS_SOLVER ":  where f=function(t,y/,p/) is a function of first derivatives\n");
    fprintf (rlab_stderr, THIS_SOLVER ":  (a column vector), 'p' the function parameter, T=[ti:tf:dt] a time\n");
    fprintf (rlab_stderr, THIS_SOLVER ":  sampling array, and  y0 the initial condition.\n");
    fprintf (rlab_stderr, THIS_SOLVER ": (2) y = odeiv(f,df,/p/,T,y0/,options/),\n");
    fprintf (rlab_stderr, THIS_SOLVER ":  where  df=function(t,y/,p/), is a Jacobian [df/dp,..,df/dt].\n");
    fprintf (rlab_stderr, THIS_SOLVER ": Check odeparams() on how to set different parameters for the solver.\n");
    rerror  (THIS_SOLVER ": requires four to six arguments!");
  }

  //
  // f = f(x/,p/)
  //
  fname_ent = bltin_get_ent(args[0]);
  if (!isfuncent(fname_ent))
    rerror (THIS_SOLVER ": " RLAB_ERROR_ARG1_FUNC_VAR "\n");

  //
  // figure out the second argument: if function use bsimp, otherwise
  // use regular odeiv solvers
  //
  e2 = bltin_get_ent(args[1]);
  if (!isfuncent(e2))
  {
    //
    // parameter entity
    //
    pent = 0;
    if (ent_type (e2) != UNDEF)
      pent = ent_Copy (e2);

    //
    // Get time sampling points
    //
    e3 = bltin_get_ent (args[2]);
    if (ent_type (e3) != MATRIX_DENSE_REAL)
      rerror ("odeiv: 'time' has to be a real row-matrix!");
    time = class_matrix_real (e3);
    nstep = SIZE (time);
    if (nstep < 2)
      rerror ("odeiv: 'time' has to have at least two entries, i.e. [tstart,tend]!");

    //
    // Get ystart
    //
    e4 = bltin_get_ent (args[3]);
    if (ent_type (e4) != MATRIX_DENSE_REAL)
      rerror ("odeiv: 'y0' has to be a real vector!");
    ystart = class_matrix_real (e4);
    if (!EQVECT(ystart))
      rerror ("odeiv: initial condition has to be a real vector");
    neq = MNR (ystart) * MNC (ystart);

    y = mdr_Create (neq, 1);
    yerr = mdr_Create (neq, 1);

    //
    // options
    //
    if (nargs > 4)
    {
      eo = bltin_get_ent (args[4]);
      if (ent_type (eo) == BTREE)
      {
        // h
        RLABCODE_PROCESS_BTREE_ENTRY_D(eo,node,RLAB_NAME_ODEIV_STEP,h0,class_double,>,0.0);

        // ay
        RLABCODE_PROCESS_BTREE_ENTRY_D(eo,node,RLAB_NAME_ODEIV_AY,ay,class_double,>,0.0);

        // adydt
        RLABCODE_PROCESS_BTREE_ENTRY_D(eo,node,RLAB_NAME_ODEIV_ADYDT,adydt,class_double,>,0.0);

        // erel
        RLABCODE_PROCESS_BTREE_ENTRY_D(eo,node,RLAB_NAME_ODEIV_EREL,relerr,class_double,>,0.0);

        // eabs
        RLABCODE_PROCESS_BTREE_ENTRY_D(eo,node,RLAB_NAME_ODEIV_EABS,abserr,class_double,>,0.0);

        // delta t for integration that cannot be completed on [t_i,t_i+1]
        // because of the user error
        RLABCODE_PROCESS_BTREE_ENTRY_D(eo,node,RLAB_NAME_ODEIV_DELTA_T,teps,class_double,>,0.0);

        // method
        RLABCODE_PROCESS_BTREE_ENTRY_D_INRANGE(eo,node,RLAB_NAME_GEN_IMETHOD,imethod,class_double,1,9);

        // phase_space
        RLABCODE_PROCESS_BTREE_ENTRY_BOOL(eo,node,RLAB_NAME_ODEIV_PHASE,ips);

        // standard output
        RLABCODE_PROCESS_BTREE_ENTRY_S(eo,node,RLAB_NAME_GEN_STDOUT,outs,1,0);
      }
      ent_Clean (eo);
    }

    if (ips)
    {
      yp  = mdr_Create (neq, 1);
      ysol = mdr_Create (nstep, 2*neq + 1);
    }
    else
      ysol = mdr_Create (nstep, neq + 1);

    //
    // t
    //
    tmdr = mdr_CreateEmpty (1,1);
    tent = ent_Assign_Rlab_MDR(tmdr);
    ent_IncRef (tent);

    //
    // y
    //
    ymdr = mdr_CreateEmpty (neq, 1);
    yent = ent_Assign_Rlab_MDR(ymdr);
    ent_IncRef (yent);

    //
    // p: nothing to do
    //

    // y0
    Mdr1 (ysol, 1, 1) = MdrV1 (time, 1);
    for (j = 2; j <= neq + 1; j++)
    {
      Mdr1  (ysol, 1, j) = MdrV1 (ystart, j - 1);
      MdrV1 (y, j - 1)  = MdrV1 (ystart, j - 1);
    }

    if (imethod > 0)
    {
      //
      // The GSL: Decide which integrator to use.
      //
      const GSL_ODEIV_STEP_TYPE *T;
      switch (imethod)
      {
      case 1:
        T = GSL_ODEIV_STEP_RK2;
        sprintf (odename, "rk2");
        break;
      case 2:
        T = GSL_ODEIV_STEP_RK4;
        sprintf (odename, "rk4");
        break;
      case 3:
        T = GSL_ODEIV_STEP_RKF45;
        sprintf (odename, "rkf45");
        break;
      case 4:
        T = GSL_ODEIV_STEP_RKCK;
        sprintf (odename, "rkck45");
        break;
      case 5:
        T = GSL_ODEIV_STEP_RK8PD;
        sprintf (odename, "rk8pd");
        break;
      default:
        T = GSL_ODEIV_STEP_RKCK;
        sprintf (odename, "rkck45");
      }

      GSL_ODEIV_STEP *s=GSL_ODEIV_STEP_ALLOC (T, neq);
      GSL_ODEIV_CONTROL *c;
      if ((ay == 0) && (adydt == 1))
        c = GSL_ODEIV_CONTROL_YP_NEW (abserr, relerr);
      else if ((ay == 1) && (adydt == 0))
        c = GSL_ODEIV_CONTROL_Y_NEW (abserr, relerr);
      else
        c = GSL_ODEIV_CONTROL_STANDARD_NEW (abserr, relerr, ay, adydt);
      GSL_ODEIV_EVOLVE *e = GSL_ODEIV_EVOLVE_ALLOC (neq);
      GSL_ODEIV_SYSTEM sys = { ode_gsl_func, NULL, neq };

      //
      // setup output stream for run-time messages
      //
      if (outs)
        fptr = fopen (outs, "a");
      if (fptr)
      {
        // start the timer
        timer = clock ();
        // write down boring info
        fprintf (fptr,
                 "RLaB: Integration of an ODE initial value problem using ");
        fprintf (fptr, "method '%s'.\n", odename);
        fprintf (fptr, "RLaB: the GSL control method is '%s'.\n",
                 GSL_ODEIV_CONTROL_NAME (c));
      }

      // phase space calculation?
      if (ips)
      {
        ode_gsl_func(MdrV0(time,0), MDRPTR(y), MDRPTR(yp), NULL);
        for (j = 1; j <= neq; j++)
          Mdr1 (ysol, 1, neq + 1 + j) = MdrV1 (yp, j);
      }
      h = h0;
      for (i = 2; i <= nstep && notnan; i++)
      {
        t1 = MdrV1 (time, i - 1);
        t2 = t2n = MdrV1 (time, i);
        do
        {
          // propagate the solution by step at most h which cannot go over t2
          // at the end of integration t1 = t2
          status =
              GSL_ODEIV_EVOLVE_APPLY (e, c, s, &sys, &t1, t2, &h, MDRPTR (y));
          if (h < h0)
            h = h0; // GSL may go bonkers and keep decreasing h so code ends up in a loop
                    // for a very long time
        }
        while (t1<t2 && status==GSL_SUCCESS);

        // check for status and if it is not GSL_SUCCESS assume that the integration bombed
        // on whatever the integration interval [t1,t2] is
        if (status != GSL_SUCCESS)
        {
          // we cannot finish the integration so we will have to reduce the output
          notnan = 0;

          // integration [t1,t2] was not successful. try manipulating t2 to t2'
          // so that integration [t1, t2'] is successfull, where t2' is determined
          // within precision teps:
          //  t2 -> attempts to contain successful upper integration bound
          //  t2n -> unsuccessful
          t2n = t2;           // t2n - this was unsuccessful
          t2  = 0.5*(t1+t2n); // reduce t2
          do
          {
            // try to integrate to new t2
            status =
                GSL_ODEIV_EVOLVE_APPLY (e, c, s, &sys, &t1, t2, &h, MDRPTR (y));
            // if integration is successful than 't1' ends up with t2

            // if unsuccessful decrease the upper bound
            if (status != GSL_SUCCESS)
              t2n = t2;
            t2  = 0.5*(t1+t2n);
          }
          while (t2n - t2 > teps);
        }

        // either way t1 is the last point of calculation
        Mdr1 (ysol, i, 1) = t1;
        // y(1) .. y(neq)
        for (j = 2; j <= neq + 1; j++)
          Mdr1 (ysol, i, j) = MdrV1 (y, j - 1);
        // phase space calculation ?
        if (ips)
        {
          ode_gsl_func(t1, MDRPTR(y), MDRPTR(yp), NULL);
          for (j = 1; j <= neq; j++)
            Mdr1 (ysol, i, neq + 1 + j) = MdrV1 (yp, j);
        }
        // print the information
        if (fptr)
        {
          fprintf (fptr, "t = %g ", Mdr1 (ysol, i, 1));
          for (j = 2; j <= neq + 1; j++)
            fprintf (fptr, "  x[%i] = %14.6g", j - 1, Mdr1 (ysol, i, j));
          if (ips)
          {
            for (j = neq+2; j <= 2*neq + 1; j++)
              fprintf (fptr, "  xdot[%i] = %14.6g", j - 1 - neq, Mdr1 (ysol, i, j));
          }
          fprintf (fptr, "\n");
        }
      }

      //
      // did we finish calculation or did it stop sooner because of NaN's
      //
      if (!notnan)
      {
        // calculation did not go to t2:
        //   resize the output matrix
        if (ips)
          ysol2 = mdr_Create(i-1, 2 * neq+1);
        else
          ysol2 = mdr_Create(i-1, neq+1);
        for (i1 = 0; i1 < i-1 ; i1++)
        {
          for (j1 = 0; j1 < neq + 1; j1++)
            Mdr0(ysol2,i1,j1) = Mdr0(ysol,i1,j1);
          if (ips)
            for (j1 = neq+1; j1 < 2*neq + 1; j1++)
              Mdr0(ysol2,i1,j1) = Mdr0(ysol,i1,j1);
        }
        mdr_Destroy (ysol);
        ysol = 0;
      }
      if (fptr)
      {
        fprintf (fptr, "RLaB: the GSL ODE integrator '%s' reports", odename);
        fprintf (fptr, " '%s' .\n", gsl_strerror (status));
        timer = clock() - timer;
        fprintf (fptr, "RLaB: ODE integration lasted %g sec.\n", timer / 1e6);
        fclose (fptr);
      }

      GSL_ODEIV_EVOLVE_FREE (e);
      GSL_ODEIV_CONTROL_FREE (c);
      GSL_ODEIV_STEP_FREE (s);
    }
    else
    {
      if (outs)
        fptr = fopen (outs, "a");
      if (fptr)
      {
        // start the timer
        timer = clock ();
        // write down boring info
        fprintf (fptr,
                 "RLaB: Integration of an ODE initial value problem using ");
        fprintf (fptr, "adams' method.\n");
      }

      //
      // Use code ODE.F to solve IV problem
      //
      int iwork[5], iflag = 1;
      int lenwrk = 100 + 21 * neq;
      MDR *work = mdr_Create (lenwrk, 1);

      // phase space calculation?
      if (ips)
      {
        ode_func(&MdrV0(time,0), MDRPTR(y), MDRPTR(yp));
        for (j = 1; j <= neq; j++)
          Mdr1 (ysol, 1, neq + 1 + j) = MdrV1 (yp, j);
      }
      //
      // propagate
      //
      for (i = 2; (i <= nstep) && notnan; i++)
      {
        t1 = MdrV1 (time, i - 1);
        t2 = t2n = MdrV1 (time, i);
        //
        // set the flag for integrator: finish integration
        // on the last time, don't go over
        //
        if (i == nstep)
          iflag = -1;
        else
          iflag =  1;
        while (t2 - t1 > teps)
        {
          // propagate the solution from t1 to t2
          ODE (ode_func, &neq, MDRPTR (y), &t1, &t2, &relerr, &abserr, &iflag,
               MDRPTR (work), iwork);
          if (iflag > 2  || iflag < -2)
          {
            // unsuccessful integration: don't go over the last
            // time step
            t2n = t2;
            t2  = t1 + 0.5 * (t2 - t1);
            iflag = -1;
            notnan = 0;
          }
          else
          {
            // successful integration: still don't go over the
            // last time step
            notnan = 1;
            t2 = t2n;
          }
        }
        // either way t1 is the last point
        Mdr1 (ysol, i, 1) = t1;
        for (j = 2; j <= neq + 1; j++)
          Mdr1 (ysol, i, j) = MdrV1 (y, j - 1);
        // phase space calculation ?
        if (ips)
        {
          ode_func(&t1, MDRPTR(y), MDRPTR(yp));
          for (j = 1; j <= neq; j++)
            Mdr1 (ysol, i, neq + 1 + j) = MdrV1 (yp, j);
        }
        //
        // print progress of integration
        //
        if (fptr)
        {
          fprintf (fptr, "t = %14.6g ", Mdr1 (ysol, i, 1));
          for (j = 2; j <= neq + 1; j++)
            fprintf (fptr, "  x[%i] = %14.6g", j - 1, Mdr1 (ysol, i, j));
          if (ips)
          {
            for (j = neq+2; j <= 2*neq + 1; j++)
              fprintf (fptr, "  xdot[%i] = %14.6g", j - 1 - neq, Mdr1 (ysol, i, j));
          }
          fprintf (fptr, "\n");
          if (iflag == 3 || iflag == 6)
            fprintf (fptr, "RLaB: ODE: NaN's appeared. "
                           "Calculation stopped.\n");
          if (iflag == 4 || iflag == 5)
            fprintf (fptr,
                     "RLaB: ODE: Could not integrate on interval"
                     " (possibly stiff problem).\n");
        }
      }
      //
      // did we finish calculation or did it stop sooner because of NaN's
      //
      if (!notnan)
      {
        // NaN's
        int is = 1;
        while (is)
        {
          --i;
          is = 0;
          for (j1 = 1; j1 < (1+ips)*neq + 1; j1++)
            is |= isnand(Mdr0(ysol,i-1,j1));
        }
        ysol2 = mdr_Create(i, (1+ips)*neq+1);
        for (i1 = 0; i1 < i; i1++)
          for (j1 = 0; j1 < (1+ips)*neq + 1; j1++)
            Mdr0(ysol2,i1,j1) = Mdr0(ysol,i1,j1);
        mdr_Destroy (ysol);
        ysol = NULL;
      }
      if (fptr)
      {
        fprintf (fptr, "RLaB: the ODE integrator reports success.\n");
        timer -= clock ();
        fprintf (fptr, "RLaB: ODE integration lasted %g sec.\n", -timer / 1e6);
        fclose (fptr);
      }
    }

    // Clean Up
    mdr_Destroy (y);
    mdr_Destroy (yerr);

    // t
    MDPTR(tmdr) = 0;
    ent_DecRef (tent);
    ent_Destroy (tent);

    // p
    ent_Clean (pent);
    ent_Clean (e2);

    // y
    MDPTR(ymdr) = 0;
    ent_DecRef (yent);
    ent_Destroy (yent);

    if (ips)
      mdr_Destroy (yp);

    ent_Clean (e3);
    ent_Clean (e4);

    if (ysol)
      return ent_Assign_Rlab_MDR(ysol);
    else
      return ent_Assign_Rlab_MDR(ysol2);
  }
  else
  {
    //
    // bsimp (the GSL) and ODEBIM solvers
    //

    // df = df(x/,p/)
    jac_fname_ent = e2;

    // options
    if (nargs > 5)
    {
      eo = bltin_get_ent (args[5]);
      if (ent_type (eo) == BTREE)
      {
        // h
        RLABCODE_PROCESS_BTREE_ENTRY_D(eo,node,RLAB_NAME_ODEIV_STEP,h,class_double,>,0.0);

        // erel
        RLABCODE_PROCESS_BTREE_ENTRY_D(eo,node,RLAB_NAME_ODEIV_EREL,relerr,class_double,>,0.0);

        // eabs
        RLABCODE_PROCESS_BTREE_ENTRY_D(eo,node,RLAB_NAME_ODEIV_EABS,abserr,class_double,>,0.0);

        // delta t for integration that cannot be completed on [t_i,t_i+1]
        // because of the user error
        RLABCODE_PROCESS_BTREE_ENTRY_D(eo,node,RLAB_NAME_ODEIV_DELTA_T,teps,class_double,>,0.0);

        // method
        RLABCODE_PROCESS_BTREE_ENTRY_D_INRANGE(eo,node,RLAB_NAME_GEN_IMETHOD,imethod,class_double,1,9);

        // phase_space
        RLABCODE_PROCESS_BTREE_ENTRY_BOOL(eo,node,RLAB_NAME_ODEIV_PHASE,ips);

        // standard output
        RLABCODE_PROCESS_BTREE_ENTRY_S(eo,node,RLAB_NAME_GEN_STDOUT,outs,1,0);
      }
      ent_Clean (eo);
    }

    // p
    e3 = bltin_get_ent (args[2]);
    pent = 0;
    if (ent_type (e3) != UNDEF)
      pent = ent_Copy (e3);

    // t
    e4 = bltin_get_ent (args[3]);
    if (ent_type (e4) != MATRIX_DENSE_REAL)
      rerror ("odeiv: " RLAB_ERROR_ARG4_MDR_VECTOR);
    time  = class_matrix_real (e4);
    nstep = SIZE(time);
    if (nstep < 2)
      rerror ("odeiv: " RLAB_ERROR_ARG4_MDR_VECTOR);

    // y0
    e5 = bltin_get_ent (args[4]);
    if (ent_type (e5) != MATRIX_DENSE_REAL)
      rerror ("odeiv: " RLAB_ERROR_ARG5_MDR_VECTOR);
    ystart = class_matrix_real (e5);
    if (!EQVECT(ystart))
      rerror ("odeiv: " RLAB_ERROR_ARG5_MDR_VECTOR);
    neq = SIZE(ystart);
    if (neq < 1)
      rerror ("odeiv: " RLAB_ERROR_ARG5_MDR_VECTOR);
    y   = mdr_Create (neq, 1);

    if (ips)
    {
      yp   = mdr_Create (neq, 1);
      ysol = mdr_Create (nstep, 2*neq + 1);
    }
    else
      ysol = mdr_Create (nstep, neq + 1);

    //
    // t
    //
    tmdr = mdr_CreateEmpty (1,1);
    tent = ent_Assign_Rlab_MDR(tmdr);
    ent_IncRef (tent);

    //
    // y
    //
    ymdr = mdr_CreateEmpty (neq, 1);
    yent = ent_Assign_Rlab_MDR(ymdr);
    ent_IncRef (yent);

    //
    // p: nothing to do
    //

    // Save initial conditions and setup y[]
    Mdr1 (ysol, 1, 1) = MdrV1 (time, 1);
    for (j = 2; j <= neq + 1; j++)
    {
      Mdr1  (ysol, 1, j) = MdrV1 (ystart, j - 1);
      MdrV1 (y, j - 1)   = MdrV1 (ystart, j - 1);
    }

    if (outs)
      fptr = fopen (outs, "a");

    if (imethod == 1)
    {
      //
      // odebim
      //

      if (fptr)
      {
        timer = clock ();
        fprintf (fptr,
                 "RLaB: using the ODEBIM integrator "
                 "for initial value problem\n");
      }

      int ijac=1, idid=0, ierr=0;
      int lwork, liwork, iout=0;

      lwork = 14 + (bim_ordmax-2) +8*neq +4*(bim_ordmax-2)*neq +2*neq*neq;
      MDR * work = mdr_Create (lwork, 1);
      MdrV1 (work, 1) = 1e-16;
      MdrV1 (work, 3) = 1e-1;
      MdrV1 (work, 4) = 1e-1;
      MdrV1 (work, 5) = 1e-1;
      MdrV1 (work, 6) = 1e-1;
      MdrV1 (work, 7) = 1e-1;
      MdrV1 (work, 8) = 5e-3;
      MdrV1 (work, 9) = 5e-2;
      MdrV1 (work, 10) = 1.2e-1;
      MdrV1 (work, 11) = 10.0;
      MdrV1 (work, 12) = 1.0/20.0;
      MdrV1 (work, 13) = MdrV1 (work, 12) / 2.0;
      MdrV1 (work, 14) = MdrV1 (work, 12);

      liwork = neq + 37;
      MDR * iwork  = mdi_Create (liwork, 1);
      MdiV1 (iwork, 1) = 100000;
      MdiV1 (iwork, 2) = bim_ordmin;
      MdiV1 (iwork, 3) = bim_ordmax;
      MdiV1 (iwork, 4) = 10;
      MdiV1 (iwork, 5) = 12;
      MdiV1 (iwork, 6) = 14;
      MdiV1 (iwork, 7) = 16;
      MdiV1 (iwork, 8) = 18;

      // phase space calculation?
      if (ips)
      {
        bim_func(&neq, &MdrV0(time,0), MDRPTR(y), MDRPTR(yp), &ierr, NULL, NULL);
        if (!ierr)
        {
          // yp successfuly calculated
          for (j = 1; j <= neq; j++)
            Mdr1 (ysol, 1, neq + 1 + j) = MdrV1 (yp, j);
        }
        else
        {
          // calculation of yp failed. is there a point to continue?
          rerror ("odeiv: calculation of RHS failed for the initial conditions!");
        }
      }
      for (i = 2; (i <= nstep) && notnan; i++)
      {
        t1 = MdrV1 (time, i - 1);
        t2 = MdrV1 (time, i);
        MdrV1 (work, 2) = (t2-t1)/8;

        BIM (&neq, bim_func, &t1, &t2, MDRPTR(y), &h, &relerr, &abserr,
              bim_fjac, &ijac, &neq, &neq,
              MDRPTR(work), &lwork, MDIPTR(iwork), &liwork,
              NULL, NULL, &iout, &idid);

        // contains last successful evaluate
        Mdr1 (ysol, i, 1) = t1;
        for (j = 2; j <= neq + 1; j++)
          Mdr1 (ysol, i, j) = MdrV1 (y, j - 1);
        if (ips)
        {
          bim_func(&neq, &t1, MDRPTR(y), MDRPTR(yp), &ierr, NULL, NULL);
          for (j = 1; (j <= neq)&&(!ierr); j++)
            Mdr1 (ysol, i, neq + 1 + j) = MdrV1 (yp, j);
        }

        //
        // print progress of integration
        //
        if (fptr)
        {
          fprintf (fptr, "t = %g ", Mdr1 (ysol, i, 1));
          for (j = 2; j <= neq + 1; j++)
            fprintf (fptr, "  x[%i] = %14.6g", j - 1, Mdr1 (ysol, i, j));
          fprintf (fptr, "\n");
        }

        //
        //
        //
        if (idid < 0)
        {
          notnan = 0;
          fprintf (fptr, "RLaB: odeiv: Solver 'odebim' reports the error code %i (",  idid);
          if (idid == -1)
            fprintf (fptr, "WRONG INPUT PARAMETERS).\n");
          else if (idid == -2)
            fprintf (fptr, "A LARGER NUMBER OF STEPS IS NEEDED).\n");
          else if (idid == -3)
            fprintf (fptr, "STEPSIZE TOO SMALL).\n");
          else if (idid == -4)
            fprintf (fptr, "REPEATEDLY SINGULAR MATRIX).\n");
          else if (idid == -5)
            fprintf (fptr, "TOO MANY CONSECUTIVE NEWTON FAILURES).\n");
          else
            fprintf (fptr, "RLaB: odeiv: UNKNOWN ERROR CODE.\n");
          fprintf (fptr, "RLaB: odeiv: Interrupting the integration!\n");
        }
      }

      //
      // was calculation successful?
      //
      if (!notnan)
      {
        // calculation did not go to t2
        ysol2 = mdr_Create(i-1, (ips+1)*neq+1);
        for (i1 = 0; i1 < i-1 ; i1++)
        {
          // solution
          for (j1 = 0; j1 < neq + 1; j1++)
            Mdr0(ysol2,i1,j1) = Mdr0(ysol,i1,j1);
          // respective time derivatives of solution (phase space)
          if (ips)
            for (j1 = 0; j1 < neq; j1++)
              Mdr0(ysol2,i1,j1 + neq + 1) = Mdr0(ysol,i1,j1);
        }
        mdr_Destroy (ysol);
        ysol = 0;
      }

      if (fptr)
      {
        fprintf (fptr, "RLaB: odeiv: Solver 'odebim' reports success.\n");
        // check the time and close the output stream
        timer -= clock ();
        fprintf (fptr,
                 "RLaB: odeiv: Integration lasted %g sec.\n",
                 -timer / 1e6);
        fclose (fptr);
      }

      // clean-up
      mdr_Destroy (work);
      mdr_Destroy (iwork);
    }
    else if (imethod == 2)
    {
      //
      // bsimp
      //

      if (fptr)
      {
        timer = clock ();
        fprintf (fptr,
                 "RLaB: using the GSL integrator 'bsimp' for initial value problem\n");
      }

      //
      // Call integrator repeatedly.
      //
      const GSL_ODEIV_STEP_TYPE * T = GSL_ODEIV_STEP_BSIMP;

      GSL_ODEIV_STEP *s = GSL_ODEIV_STEP_ALLOC (T, neq);

      GSL_ODEIV_CONTROL *c;
      if (ay == 1 && adydt == 0)
        c = GSL_ODEIV_CONTROL_Y_NEW (abserr, relerr);
      else if (ay == 0 && adydt == 1)
        c = GSL_ODEIV_CONTROL_YP_NEW (abserr, relerr);
      else
        c = GSL_ODEIV_CONTROL_STANDARD_NEW (abserr, relerr, ay, adydt);

      GSL_ODEIV_EVOLVE *e = GSL_ODEIV_EVOLVE_ALLOC (neq);
      GSL_ODEIV_SYSTEM sys = { ode_gsl_func, ode_gsl_fjac, neq, NULL };

      //
      // propagate
      //
      GSL_ODEIV_EVOLVE_RESET (e);

      // phase space calculation?
      if (ips)
      {
        ode_gsl_func(MdrV0(time,0), MDRPTR(y), MDRPTR(yp), NULL);
        for (j = 1; j <= neq; j++)
          Mdr1 (ysol, 1, neq + 1 + j) = MdrV1 (yp, j);
      }

      if (h<h0)
        h = h0;

      for (i = 2; i <= nstep && notnan; i++)
      {
        t1 = MdrV1 (time, i - 1);
        t2 = t2n = MdrV1 (time, i);
        //while (t1 < t2)
        while (t2 - t1 > teps)
        {
          // propagate the solution by step at most h which cannot go over t2
          // at the end of integration t1 = t2
          status =
            GSL_ODEIV_EVOLVE_APPLY (e, c, s, &sys, &t1, t2, &h, MDRPTR (y));

          if (status != GSL_SUCCESS)
          {
            // unsuccesful integration:
            //  make note that the current upper bound was unsuccesful
            //  then lower the upper bound using bisection
            t2n = t2;
            t2  = t1 + 0.5 * (t2 - t1);
            notnan = 0;
          }
          else
          {
            notnan = 1;
            t2 = t2n;
          }
        }
        // either way t1 is the last point of calculation
        Mdr1 (ysol, i, 1) = t1;
        // y(1) .. y(neq)
        for (j = 2; j <= neq + 1; j++)
          Mdr1 (ysol, i, j) = MdrV1 (y, j - 1);
        // phase space calculation ?
        if (ips)
        {
          ode_gsl_func(t1, MDRPTR(y), MDRPTR(yp), NULL);
          for (j = 1; j <= neq; j++)
            Mdr1 (ysol, i, neq + 1 + j) = MdrV1 (yp, j);
        }
        // print the information
        if (fptr)
        {
          fprintf (fptr, "t = %g ", Mdr1 (ysol, i, 1));
          for (j = 2; j <= neq + 1; j++)
            fprintf (fptr, "  x[%i] = %14.6g", j - 1, Mdr1 (ysol, i, j));
          if (ips)
          {
            for (j = neq+2; j <= 2*neq + 1; j++)
              fprintf (fptr, "  xdot[%i] = %14.6g", j - 1 - neq, Mdr1 (ysol, i, j));
          }
          fprintf (fptr, "\n");
        }
      }
      if (!notnan)
      {
        // calculation did not go to t2
              // phase space calculation?
        if (ips)
          ysol2 = mdr_Create(i-1, 2*neq+1);
        else
          ysol2 = mdr_Create(i-1, neq+1);
        for (i1 = 0; i1 < i-1 ; i1++)
        {
          for (j1 = 0; j1 < neq + 1; j1++)
            Mdr0(ysol2,i1,j1) = Mdr0(ysol,i1,j1);
          if (ips)
            for (j1 = neq+1; j1 < 2*neq + 1; j1++)
              Mdr0(ysol2,i1,j1) = Mdr0(ysol,i1,j1);
        }
        mdr_Destroy (ysol);
        ysol = 0;
      }

      //
      // print end of integration
      //
      if (fptr)
      {
        fprintf (fptr, "RLaB: the GSL ODE integrator 'bsimp' reports success.\n");
        // check the time and close the output stream
        timer = clock() - timer;
        fprintf (fptr,
                 "RLaB: ODE integration using 'bsimp' lasted %g sec.\n",
                 timer / 1e6);
        fclose (fptr);
      }

      GSL_ODEIV_EVOLVE_FREE (e);
      GSL_ODEIV_CONTROL_FREE (c);
      GSL_ODEIV_STEP_FREE (s);
    }
    else
    {
      //
      // vode
      //

      if (fptr)
      {
        timer = clock ();
        fprintf (fptr,
                 "RLaB: odeiv: using VODE integrator for initial value problem\n");
      }

      int itask=1, istate=1, iopt=0, lrw, liw, mf=21;
      int itol = 2;
      MDR *rwork, *iwork;

      lrw = 22 +  9*neq + 2*neq*neq;
      rwork = mdr_Create(lrw,1);
      liw = 30 + neq;
      iwork = mdi_Create(liw,1);

      t1  = MdrV1 (time, 1);
      // phase space calculation?
      if (ips)
      {
        vode_func(&neq, &t1, MDRPTR(y), MDRPTR(yp), NULL, NULL);
        for (j = 1; j <= neq; j++)
          Mdr1 (ysol, 1, neq + 1 + j) = MdrV1 (yp, j);
      }

      for (i = 2; (i <= nstep) && notnan; i++)
      {
        t2 = MdrV1 (time, i);

        // do not integrate beyond last time point
        itask = 4;
        MdrV1(rwork,1) = t2;

        DVODE (vode_func, &neq, MDRPTR(y), &t1, &t2, &itol, &relerr, &abserr, &itask,
               &istate, &iopt, MDRPTR(rwork), &lrw, MDIPTR(iwork), &liw,
               vode_fjac, &mf, NULL, NULL);

        Mdr1 (ysol, i, 1) = t1;
        for (j = 2; j <= neq + 1; j++)
          Mdr1 (ysol, i, j) = MdrV1 (y, j - 1);
        if (ips)
        {
          vode_func(&neq, &t1, MDRPTR(y), MDRPTR(yp), NULL, NULL);
          for (j = 1; j <= neq; j++)
            Mdr1 (ysol, 1, neq + 1 + j) = MdrV1 (yp, j);
        }


        //
        // print progress of integration
        //
        if (fptr)
        {
          fprintf (fptr, "t = %g", Mdr1 (ysol, i, 1));
          for (j = 2; j <= neq + 1; j++)
            fprintf (fptr, ", x[%i] = %14.6f", j - 1, Mdr1 (ysol, i, j));
          if (ips)
          {
            for (j = neq+2; j <= 2*neq + 1; j++)
              fprintf (fptr, ", xdot[%i] = %14.6f", j - 1 - neq, Mdr1 (ysol, i, j));
          }
          fprintf (fptr, "\n");
        }

        //
        // check for error conditions
        //
        if (istate < 0)
        {
          fprintf (fptr, "RLaB: odeiv: VODE solver reports error code %i (",  istate);
          if (istate == -2)
            fprintf (fptr, "Tolerances too small).\n");
          else if (istate == -5)
            fprintf (fptr, "Bad Jacobian).\n");
          else
            fprintf (fptr, "UNKNOWN ERROR CODE).\n");
          fprintf (fptr, "RLaB: odeiv: Interrupting the calculation!\n");
        }
      }

      //
      // was calculation successful?
      //
      if (!notnan)
      {
        if (ips)
          ysol2 = mdr_Create(i-1, 2*neq+1);
        else
          ysol2 = mdr_Create(i-1, neq+1);
        for (i1 = 0; i1 < i-1 ; i1++)
        {
          for (j1 = 0; j1 < neq + 1; j1++)
            Mdr0(ysol2,i1,j1) = Mdr0(ysol,i1,j1);
          if (ips)
            for (j1 = 0; j1 < neq; j1++)
              Mdr0(ysol2,i1,j1 + neq + 1) = Mdr0(ysol,i1,j1);
        }
        mdr_Destroy (ysol);
        ysol = 0;
      }

      if (fptr)
      {
        if (notnan)
          fprintf (fptr, "RLaB: odeiv: VODE solver reports success.\n");
        // check the time and close the output stream
        timer -= clock ();
        fprintf (fptr,
                 "RLaB: odeiv: Integration lasted %g sec.\n",
                 -timer / 1e6);
        fclose (fptr);
      }

      // clean-up
      mdr_Destroy (rwork);
      mdr_Destroy (iwork);
    }

    //
    // Clean Up
    //
    mdr_Destroy (y);
    if (ips)
      mdr_Destroy (yp);

    MDPTR(tmdr) = 0;
    ent_DecRef  (tent);
    ent_Destroy (tent);

    ent_Clean (fname_ent);
    ent_Clean (e2);
    ent_Clean (e3);
    ent_Clean (pent);

    MDPTR(ymdr) = 0;
    ent_DecRef (yent);
    ent_Destroy (yent);

    ent_Clean (e4);
    ent_Clean (e5);

    if (ysol)
      return ent_Assign_Rlab_MDR(ysol);
    else
      return ent_Assign_Rlab_MDR(ysol2);
  }
}

//
// derivative function for odeiv solver, conforms to gsl
//
static int
ode_gsl_func (double t, const double y[], double yp[], void *dummy)
{
  int i, notnan = 1;
  Ent *rent = 0;
  MDR *retm = 0;

  // t
  MDPTR(tmdr) = (void *) &t;

  // y
  MDPTR(ymdr) = (void *) y;

  // p
  if (pent)
    rent = ent_call_rlab_script_3args(fname_ent, tent, yent, pent);
  else
    rent = ent_call_rlab_script_2args(fname_ent, tent, yent);

  if (ent_type(rent)!=MATRIX_DENSE_REAL)
    rerror (THIS_SOLVER ": " RLAB_ERROR_RHS_FUNC_MUST_RETURN_MDR);

  retm = ent_data(rent);

  if (SIZE(retm) != neq)
  {
    fprintf(stderr, THIS_SOLVER " [gsl]: " RLAB_ERROR_RHS_FUNC_INCORRECT_DIM "\n");
    fprintf(stderr, THIS_SOLVER " [gsl]: Expected dimension %i, RHS dimension %i\n", neq, SIZE(retm));
    fprintf(stderr, THIS_SOLVER " [gsl]: Please check the RHS function.\n");
    rerror (THIS_SOLVER);
  }

  for (i = 0; (i < neq) && notnan; i++)
  {
    yp[i] = MdrV0 (retm, i);
    notnan = notnan & !isnand( yp[i] );
  }

  ent_Clean(rent);

  // check for nans
  if (notnan)
    return GSL_SUCCESS;
  else
    return GSL_ERANGE;
}

//
// jacobian function for odeiv solver bsimp
//
static int
ode_gsl_fjac (double t, const double y[], double *dfdy,
                 double *dfdt, void *dummy)
{
  int i, j;
  Ent *rent;
  MDR *retm;

  gsl_matrix_view dfdy_mat = gsl_matrix_view_array (dfdy, neq, neq);
  gsl_matrix *m = &dfdy_mat.matrix;

  // t
  MDPTR(tmdr) = (void *) &t;

  // y
  MDPTR(ymdr) = (void *) y;

  // p
  if (pent)
    rent = ent_call_rlab_script_3args(jac_fname_ent, tent, yent, pent);
  else
    rent = ent_call_rlab_script_2args(jac_fname_ent, tent, yent);

  if (ent_type(rent)!=MATRIX_DENSE_REAL)
    rerror (THIS_SOLVER ": " RLAB_ERROR_TERRIBLE_INTERNAL_ERROR);

  retm = ent_data (rent);

  if (SIZE(retm)!= neq * (neq + 1))
  {
    fprintf(stderr, THIS_SOLVER " [gsl]: " RLAB_ERROR_JACOBIAN_FUNC_INCORRECT_DIM "\n");
    fprintf(stderr, THIS_SOLVER " [gsl]: Expected dimension %ix%i, JACOBIAN dimension %ix%i\n", neq, neq + 1,
            MNR(retm), MNC(retm));
    fprintf(stderr, THIS_SOLVER " [gsl]: Please check the JACOBIAN function.\n");
    rerror (THIS_SOLVER);
  }

  for (i = 0; i < neq; i++)
  {
    for (j = 0; j < neq; j++)
    {
      gsl_matrix_set (m, i, j, Mdr0 (retm, i, j));
    }
  }

  for (i = 0; i < neq; i++)
    dfdt[i] = Mdr0 (retm, i, neq);

  ent_Clean(rent);

  return GSL_SUCCESS;
}

//
// The interface to the user-specified function for solver ODE.F
// N.B.: the solver does not have internal error control mechanism,
// so if something bad happened we just return nan()'s.
//
static int
ode_func (double *t, double *y, double *yp)
{
  int i;
  Ent *rent = 0;
  MDR *retm = 0;

  // t
  MDPTR(tmdr) = (void *) t;

  // y
  MDPTR(ymdr) = (void *) y;

  // p
  if (pent)
    rent = ent_call_rlab_script_3args(fname_ent, tent, yent, pent);
  else
    rent = ent_call_rlab_script_2args(fname_ent, tent, yent);

  if (ent_type(rent)!=MATRIX_DENSE_REAL)
    rerror ("odeiv: terrible internal error");
  retm = ent_data(rent);

  if (SIZE(retm)!=neq)
  {
    fprintf(stderr, THIS_SOLVER " [ode.f]: " RLAB_ERROR_RHS_FUNC_INCORRECT_DIM "\n");
    fprintf(stderr, THIS_SOLVER " [ode.f]: Expected dimension %i, RHS dimension %i\n", neq, SIZE(retm));
    fprintf(stderr, THIS_SOLVER " [ode.f]: Please check the RHS function.\n");
    rerror (THIS_SOLVER);
  }

  for (i = 0; (i < neq); i++)
    yp[i] = MdrV0 (retm, i);

  ent_Clean(rent);
  return (1);
}

static int
bim_func (int * m, double * t, double * y, double * dy,
          int * ierr, double * rpar, int * ipar)
{
  int i;
  Ent *rent = 0;
  MDR *retm = 0;

  // t
  MDPTR(tmdr) = (void *) t;

  // y
  MDPTR(ymdr) = (void *) y;

  // p
  if (pent)
    rent = ent_call_rlab_script_3args(fname_ent, tent, yent, pent);
  else
    rent = ent_call_rlab_script_2args(fname_ent, tent, yent);

  if (ent_type(rent)!=MATRIX_DENSE_REAL)
    rerror ("odeiv: [bim] terrible internal error");

  retm = ent_data(rent);

  if (SIZE(retm)!=neq)
  {
    fprintf(stderr, THIS_SOLVER " [bim]: " RLAB_ERROR_RHS_FUNC_INCORRECT_DIM "\n");
    fprintf(stderr, THIS_SOLVER " [bim]: Expected dimension %i, RHS dimension %i\n", neq, SIZE(retm));
    fprintf(stderr, THIS_SOLVER " [bim]: Please check the RHS function.\n");
    rerror (THIS_SOLVER);
  }

  for (i = 0; (i < neq)&&(! *ierr); i++)
  {
    dy[i] = MdrV0 (retm, i);
    *ierr |= isnand(dy[i]);
  }

  ent_Clean(rent);
  return 1;
}

static int
bim_fjac (int * m, double * t, double * y, double * jac, int * ldjac, int * ierr,
          double * rpar, int * ipar)
{
  int i, j;
  Ent *rent = 0;
  MDR *retm = 0;

  // t
  MDPTR(tmdr) = (void *) t;

  // y
  MDPTR(ymdr) = (void *) y;

  // p
  if (pent)
    rent = ent_call_rlab_script_3args(jac_fname_ent, tent, yent, pent);
  else
    rent = ent_call_rlab_script_2args(jac_fname_ent, tent, yent);

  if (ent_type(rent)!=MATRIX_DENSE_REAL)
    rerror (THIS_SOLVER " [bim]: " RLAB_ERROR_JACOBIAN_FUNC_MUST_RETURN_MDR);

  retm = ent_data(rent);

  if ((SIZE (retm) != neq * neq) && (SIZE (retm) != neq * (neq+1) ))
  {
    fprintf(stderr, THIS_SOLVER " [bim]: " RLAB_ERROR_JACOBIAN_FUNC_INCORRECT_DIM "\n");
    fprintf(stderr, THIS_SOLVER " [bim]: Expected dimension %ix%i or %ix%i, JACOBIAN dimension %ix%i\n",
            neq, neq + 1, neq, neq + 1, MNR(retm), MNC(retm));
    fprintf(stderr, THIS_SOLVER " [bim]: Please check the JACOBIAN function.\n");
    rerror (THIS_SOLVER);
  }

  for (i = 0; i < neq; i++)
    for (j = 0; j < neq; j++)
      jac[j*neq + i] = Mdr0 (retm, i, j);

  ent_Clean(rent);
  return 1;
}

static int
vode_func (int * n, double * t, double * y, double * ydot,
           double * rpar, int * ipar)
{
  int i;
  Ent *rent = 0;
  MDR *retm = 0;

  // t
  MDPTR(tmdr) = (void *) t;

  // y
  MDPTR(ymdr) = (void *) y;

  // p
  if (pent)
    rent = ent_call_rlab_script_3args(fname_ent, tent, yent, pent);
  else
    rent = ent_call_rlab_script_2args(fname_ent, tent, yent);

  if (ent_type(rent)!=MATRIX_DENSE_REAL)
    rerror (THIS_SOLVER " [vode]: " RLAB_ERROR_RHS_FUNC_MUST_RETURN_MDR);

  retm = ent_data(rent);

  if (SIZE(retm)!=neq)
  {
    fprintf(stderr, THIS_SOLVER " [vode]: " RLAB_ERROR_RHS_FUNC_INCORRECT_DIM "\n");
    fprintf(stderr, THIS_SOLVER " [vode]: Expected dimension %i, RHS dimension %i\n", neq, SIZE(retm));
    fprintf(stderr, THIS_SOLVER " [vode]: Please check the RHS function.\n");
    rerror (THIS_SOLVER);
  }

  for (i = 0; i < neq; i++)
    ydot[i] = MdrV0 (retm, i);

  ent_Clean(rent);
  return 1;
}


static int
vode_fjac (int * n, double * t, double * y,
           int * ml, int * mu, double * pd,
           int * nrpd, double * rpar, int * ipar)
{
  int i, j;

  Ent *rent = 0;
  MDR *retm = 0;

  // t
  MDPTR(tmdr) = (void *) t;

  // y
  MDPTR(ymdr) = (void *) y;

  // p
  if (pent)
    rent = ent_call_rlab_script_3args(jac_fname_ent, tent, yent, pent);
  else
    rent = ent_call_rlab_script_2args(jac_fname_ent, tent, yent);

  if (ent_type(rent)!=MATRIX_DENSE_REAL)
    rerror (THIS_SOLVER " [vode]: " RLAB_ERROR_JACOBIAN_FUNC_MUST_RETURN_MDR);

  retm = ent_data(rent);

  if ((SIZE(retm) != (*n) * (*n)) && (SIZE(retm) != (*n) * (*n+1)))
  {
    fprintf(stderr, THIS_SOLVER " [vode]: " RLAB_ERROR_JACOBIAN_FUNC_INCORRECT_DIM "\n");
    fprintf(stderr, THIS_SOLVER " [vode]: Expected dimension %ix%i or %ix%i, JACOBIAN dimension %ix%i\n",
            neq, neq + 1, neq, neq + 1, MNR(retm), MNC(retm));
    fprintf(stderr, THIS_SOLVER " [vode]: Please check the JACOBIAN function.\n");
    rerror (THIS_SOLVER);
  }

  for (i = 0; i < *n; i++)
    for (j = 0; j < *n; j++)
      pd[j*(*n) + i] = Mdr0 (retm, i, j);

  ent_Clean(rent);
  return 1;
}
