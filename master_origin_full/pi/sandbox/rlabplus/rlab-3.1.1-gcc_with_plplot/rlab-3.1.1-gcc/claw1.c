// Copyright (C) 2003-2014 Marijan Kostrun
//   part of rlab+4linux project on rlabplus.sourceforge.net
//
// RLaB2 Rel.2 driver for CLAWPACK 1-D
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

// standard libraries
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <time.h>

#include "clawpack.h"
#include "ode.h"

#include "rfileio.h"
#include "rlab_solver_parameters_names.h"

//
// user function used by clawpack
//
// functions
static MDR *
claw1k (void);

// their rlab names and parameter lists
static Ent *rname=0, *sname=0, *bclname=0, *bcrname=0, *b4name=0, *kname=0;

// communication between rlab functions and the clawpack
static Ent *qL_ent , *qR_ent, *q_ent, *x_ent, *x1_ent, *t_ent;
static MDR *qL, *qR, *q, *x, *x1, *t, *claw1_retk;

//
// clawpack 1-D internal variables
//
// claw1_order:
//  1 for 1st order Godunove,
//  2, for high resolution correction
static int claw1_order = 1;

// claw1_b4step1
// 0, do not call b4step1 function
// 1, call it
static int claw1_b4step1 = 0;

// claw1_mbc: number of cell layers at each boundary
static int claw1_mbc = 2;

// claw1_tsrc: 0 for 'adams', 1 for 'Collitz'
static int claw1_tsrc = 0;

// claw1_relerr: relative error in time integration with source function
static double claw1_relerr = 0.0;

// claw1_abserr: absolute error
static double claw1_abserr = 1e-6;

// claw1_m5_fracsplitmethod: source term 1 for Godunov, 2 for Strang splitting
static int claw1_m5_fracsplitmethod = 1;

static int claw1_method[7] =
{
  1,				// use flexible time step. 0 is for fixed
  0,				// set by claw1_order at the beginning of basic code
  0,				// not used originally, do the  b4step1  ?
  0,				// verbosity
  0,				// source fn, 0 for no, otherwise see  claw1_fracsplitmethod  above
  0,				// capacity fn, 0 for no, otherwise  mcapa
  0				// maux, size of aux array, maux >=mcapa
};

static int claw1_meqn;
static int claw1_nrow;
static int claw1_mwav;

static double claw1_cflv[4] = {
  1.0,				// max courant number allowed
  0.9,				// desired courant number for the computation
  0.0,
  0.0
};

static int claw1_nv[2] = {
  50000,			// max number of calls to clawpack
  0
};

// claw1_mthbc
//    type of boundary condition on left and right
//    left: 0 user, 1 zero-order extr, 2 periodic, 3 solid wall
//    right: 0 user, 1 zero-order extr, 2 periodic, 3 solid wall
static int claw1_mthbc[2] = {
  1,
  1
};

// claw1_mthlim: limiter for the riemann solver
static int claw1_mthlim[10] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

static double claw1_mindt = 1e-4;
static double claw1_maxdt = 1;

//
// clawparams list: limiters for fields in a solution
//
Ent *
ent_claw1par_lim (int nargs, Datum args[])
{
  int i, j;
  Ent *e1=0;
  //
  // Check arguments.
  //
  if (nargs == 0)
  {
    fprintf (stdout,
             "lim: A member of the auxiliary function list  clawparams .  Sets the\n");
    fprintf (stdout,
             "lim: type of limiter on each component in a solution.\n");
    fprintf (stdout,
             "lim: Format:\n");
    fprintf (stdout,
             "lim:   clawparams.lim(lim1,lim2,..),\n");
    fprintf (stdout,
             "lim: where  lim1,lim2,.. = 0 for no limiter (Lax-Wendroff scheme, default),\n");
    fprintf (stdout,
             "lim: 1 for minmod,  2 for superbee, 3  for van Leer's and 4 for monotonized\n");
    fprintf (stdout,
             "lim: centered limiter.\n");
    rerror ("requires at least one argument");
  }
  for (i = 0; i < nargs; i++)
  {
    e1 = bltin_get_ent (args[i]);
    if (ent_type (e1) != MATRIX_DENSE_REAL)
      rerror ("clawparams.lim: Improper argument!");
    j = (int) class_double (e1);
    if (j >= 0 && j <= 4)
      claw1_mthlim[i] = j;

    ent_Clean (e1);
  }

  return ent_Create_Rlab_Success();
}

//
// clawparams list: method and precision
//
#undef THIS_SOLVER
#define THIS_SOLVER "clawparams.b4s1"
Ent *
ent_claw1par_b4s1 (int nargs, Datum args[])
{
  FILE *rlab_stderr = (!RLAB_STDERR_DS) ? stderr : RLAB_STDERR_DS;

  //
  // Check arguments.
  //
  if (nargs != 1)
  {
    fprintf (rlab_stderr,
             THIS_SOLVER ": A member of the auxiliary function list  clawparams .  Sets the\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": function to be called during the basic stepping algorithm of\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": claw1  solver.  Format:\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ":   clawparams.b4s1(func),\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": where   func=function(t,x,q)  provides an additional modification of\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": the solution, q<-func(t,x,q), due to causes not covered by a regular\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": evolution. clawparams.b4s1(0) resets the solver not to call the function.\n");
    rerror  (THIS_SOLVER ": " RLAB_ERROR_ONE_ARG_REQUIRED);
  }

  if (b4name)
    ent_Clean (b4name);

  b4name = bltin_get_ent(args[0]);

  if (isfuncent(b4name))
    claw1_b4step1 = 1;
  else
  {
    claw1_b4step1 = 0;
    ent_Clean (b4name);
    b4name = 0;
  }

  return ent_Create_Rlab_Success();
}

//
// clawparams list: method and precision
//
#undef THIS_SOLVER
#define THIS_SOLVER "clawparams.bc"
Ent *
ent_claw1par_bc (int nargs, Datum args[])
{
  int i;
  Ent *e1=0;
  FILE *rlab_stderr = (!RLAB_STDERR_DS) ? stderr : RLAB_STDERR_DS;

  //
  // Check arguments.
  //
  if (nargs < 1 || nargs > 5)
  {
    fprintf (rlab_stderr,
             THIS_SOLVER ": A member of the auxiliary function list  clawparams .  Sets the\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": type of boundary conditions and, if necessary, the user defined\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": functions that determine them.  Format:\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ":   clawparams.bc(mbc,typeL,typeR/,funcL,funcR/),\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": where  mbc is a number of ghost cell layers at each boundary, while\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": boundary typeL,R=0 (specified by user function 'funcL,R'), 1 (zeroth\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": order extrapolation), 2 (periodic), 3 (solid wall).  The Functions\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": funcL,R = function(t,x,q)  provide  the boundary conditions, where\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": q<-funcL,R(t,x,q).\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": Current values:\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ":   mbc   = %i\n", claw1_mbc);
    fprintf (rlab_stderr,
             THIS_SOLVER ":   typeL = %i\n", claw1_mthbc[0]);
    fprintf (rlab_stderr,
             THIS_SOLVER ":   typeR = %i\n", claw1_mthbc[1]);
    rerror  (THIS_SOLVER ": " RLAB_ERROR_ONE_TO_FIVE_ARG_REQUIRED);
  }

  if (nargs == 1)
  {
    e1 = bltin_get_ent (args[0]);
    if (ent_type (e1) != MATRIX_DENSE_REAL)
      rerror (THIS_SOLVER ": " RLAB_ERROR_ARG1_MDR_INTEGER);
    i = (int) class_double (e1);
    if (i > 1)
      claw1_mbc = i;

    ent_Clean(e1);
  }
  if (nargs > 1)
  {
    e1 = bltin_get_ent (args[1]);
    if (ent_type (e1) != MATRIX_DENSE_REAL)
      rerror (THIS_SOLVER ": " RLAB_ERROR_ARG2_MDR_INTEGER);
    i = (int) class_double (e1);
    if (i > 0 && i < 4)
      claw1_mthbc[0] = i;

    ent_Clean(e1);
  }
  if (nargs > 2)
  {
    e1 = bltin_get_ent (args[2]);
    if (ent_type (e1) != MATRIX_DENSE_REAL)
      rerror (THIS_SOLVER ": " RLAB_ERROR_ARG3_MDR_INTEGER);
    i = (int) class_double (e1);
    if (i > 0 && i < 4)
      claw1_mthbc[1] = i;

    ent_Clean(e1);
  }
  //
  // Get the name of the function for the left boundary condition function
  //
  if (nargs > 3 && claw1_mthbc[0] == 0)
  {
    bclname = bltin_get_ent(args[3]);
    if (!isfuncent(bclname))
    {
      // user messed up. Set the default value on the left boundary
      claw1_mthbc[0] = 1;
      fprintf (rlab_stderr,
               THIS_SOLVER ": Serious error.  Tried to set left boundary condition function.\n");
      fprintf (rlab_stderr,
               THIS_SOLVER ": 'typeL' set to default value of 1 (zeroth order extrapolation).\n");
      ent_Clean (bclname);
      rerror (THIS_SOLVER ": " RLAB_ERROR_ARG4_FUNC_VAR "\n");
    }
  }

  //
  // Get the name of the function for the right boundary condition function
  //
  if (nargs > 4 && claw1_mthbc[1] == 0)
  {
    bcrname = bltin_get_ent(args[4]);
    if (!isfuncent(bcrname))
    {
      // user messed up. Set the default value on the left boundary
      claw1_mthbc[1] = 3;
      fprintf (rlab_stderr,
               THIS_SOLVER ": Serious error.  Tried to set right boundary condition function.\n");
      fprintf (rlab_stderr,
               THIS_SOLVER ": 'typeR' set to default value of 3 (solid wall).\n");
      ent_Clean (bcrname);
      rerror (THIS_SOLVER ": " RLAB_ERROR_ARG4_FUNC_VAR "\n");
    }
  }

  return ent_Create_Rlab_Success();
}


//
// clawparams list: method and precision
//
#undef THIS_SOLVER
#define THIS_SOLVER "clawparams.method"
Ent *
ent_claw1par_method (int nargs, Datum args[])
{
  Ent *e1=0;
  FILE *rlab_stderr = (!RLAB_STDERR_DS) ? stderr : RLAB_STDERR_DS;

  //
  // Check arguments.
  //
  if (nargs < 1 || nargs > 2)
  {
    fprintf (rlab_stderr,
             THIS_SOLVER ": A member of the auxiliary function list  clawparams.  Sets  the\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": integration method and the time splitting scheme for the source\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": term.  Format:\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ":   clawparams.method(order,sourcesp),\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": where  order=1 (1st order Godunov method, default), 2 (high\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": resolution correction, sourcesp=1 (Godunov splitting, default),\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": 2 (Strang splitting).\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": Current values:\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ":   order    = %i,\n", claw1_order);
    fprintf (rlab_stderr,
             THIS_SOLVER ":   sourcesp = %i.\n", claw1_m5_fracsplitmethod);
    rerror ("requires at least one argument");
  }

  if (nargs == 1)
  {
    e1 = bltin_get_ent (args[0]);
    if (ent_type (e1) != MATRIX_DENSE_REAL)
      rerror (THIS_SOLVER ": " RLAB_ERROR_ARG1_MDR_INTEGER);
    claw1_order = 3 && (int) class_double (e1);
    if (!claw1_order)
      claw1_order = 1;		// default value

    ent_Clean (e1);
  }
  if (nargs == 2)
  {
    e1 = bltin_get_ent (args[1]);
    if (ent_type (e1) != MATRIX_DENSE_REAL)
      rerror (THIS_SOLVER ": " RLAB_ERROR_ARG2_MDR_INTEGER);
    claw1_m5_fracsplitmethod = 3 && (int) class_double (e1);
    if (!claw1_m5_fracsplitmethod)
      claw1_m5_fracsplitmethod = 1;	// default value

    ent_Clean (e1);
  }

  return ent_Create_Rlab_Success();
}

//
// clawparams list: source term integration
//
#undef THIS_SOLVER
#define THIS_SOLVER "clawparams.src"
Ent *
ent_claw1par_tsrc (int nargs, Datum args[])
{
  Ent *e1=0;
  char *odename=0;
  FILE *rlab_stderr = (!RLAB_STDERR_DS) ? stderr : RLAB_STDERR_DS;

  //
  // Check arguments.
  //
  if (nargs < 1 || nargs > 3)
  {
    fprintf (rlab_stderr,
             THIS_SOLVER ": A member of the auxiliary function list  clawparams.  Sets  the\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": integration method and the absolute and relative errors for the time\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": integration of the source term.\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": Format:\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ":   clawparams.src(ode/,abserr,relerr/),\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": where  ode=\"adams\" (implicite/explicite, default),  \"rk2\" (2nd\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": order, explicite runge-kutta), 'abserr' and 'relerr' are the errors in\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": integration for \"adams\" integrator.\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": Current values:\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ":   ode    = ");
    if (claw1_tsrc == 0)
    {
      fprintf (rlab_stderr, "\"adams\",\n");
      fprintf (rlab_stderr, THIS_SOLVER ":   abserr = %g,\n", claw1_abserr);
      fprintf (rlab_stderr, THIS_SOLVER ":   relerr = %g.\n", claw1_relerr);
    }
    else
      fprintf (rlab_stderr, "\"rk2\",\n");

    rerror (THIS_SOLVER ": " RLAB_ERROR_AT_LEAST_ONE_ARG_REQUIRED);
  }

  if (nargs == 1)
  {
    e1 = bltin_get_ent (args[0]);
    if (ent_type (e1) == MATRIX_DENSE_STRING)
    {
      odename = class_char_pointer (e1);
      claw1_tsrc = 0;
      if (isvalidstring(odename))
        if (strcmp ("rk2", odename) == 0)
          claw1_tsrc = 1;
    }
    ent_Clean (e1);
  }
  if (nargs >= 2)
  {
    e1 = bltin_get_ent (args[1]);
    if (ent_type (e1) != MATRIX_DENSE_REAL)
      rerror (THIS_SOLVER ": " RLAB_ERROR_ARG2_MDR_SCALAR);
    claw1_abserr = class_double (e1);
    if (claw1_abserr <= 0)
      claw1_abserr = 1e-6;
    ent_Clean (e1);
  }
  if (nargs == 3)
  {
    e1 = bltin_get_ent (args[2]);
    if (ent_type (e1) != MATRIX_DENSE_REAL)
      rerror (THIS_SOLVER ": " RLAB_ERROR_ARG3_MDR_SCALAR);
    claw1_relerr = class_double (e1);
    if (claw1_relerr < 0)
      claw1_relerr = 0.0;
    ent_Clean (e1);
  }

  return ent_Create_Rlab_Success();
}

//
// clawparams list: courant numbers and maximum number of iterations
//
#undef THIS_SOLVER
#define THIS_SOLVER "clawparams.cour"
Ent *
ent_claw1par_courmaxit (int nargs, Datum args[])
{
  Ent *e1=0;
  double w;
  int i=0;
  FILE *rlab_stderr = (!RLAB_STDERR_DS) ? stderr : RLAB_STDERR_DS;

  //
  // Check arguments.
  //
  if (nargs < 1 || nargs > 3)
  {
    fprintf (rlab_stderr,
             THIS_SOLVER ": A member of the auxiliary function list  clawparams.\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": Sets courant numbers and maximum number of iterations for\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": claw1 solver.  Format:\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ":   clawparams.cour(maxcour,cour,maxit),\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": where  'maxcour' is the maximal courant number allowed, 'cour'\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": is the desired courant number and 'maxit' is the maximum number\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": of iterations allowed to reach the desired tolerance.\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": Current values:\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ":   maxc  = %g,\n", claw1_cflv[0]);
    fprintf (rlab_stderr,
             THIS_SOLVER ":   cour  = %g,\n", claw1_cflv[1]);
    fprintf (rlab_stderr,
             THIS_SOLVER ":   maxit = %i.\n", claw1_nv[0]);

    rerror (THIS_SOLVER ": " RLAB_ERROR_AT_LEAST_ONE_ARG_REQUIRED);
  }
  if (nargs == 1)
  {
    e1 = bltin_get_ent (args[0]);
    if (ent_type (e1) != MATRIX_DENSE_REAL)
      rerror (THIS_SOLVER ": " RLAB_ERROR_ARG1_MDR_SCALAR);
    w = class_double (e1);
    if (w > 0)
      claw1_cflv[0] = w;

    ent_Clean (e1);
  }
  if (nargs > 1)
  {
    e1 = bltin_get_ent (args[1]);
    if (ent_type (e1) != MATRIX_DENSE_REAL)
      rerror (THIS_SOLVER ": " RLAB_ERROR_ARG2_MDR_SCALAR);
    w = class_double (e1);
    if (w > 0 && w < claw1_cflv[0])
      claw1_cflv[1] = w;

    ent_Clean (e1);
  }
  if (nargs == 3)
  {
    e1 = bltin_get_ent (args[2]);
    if (ent_type (e1) == MATRIX_DENSE_REAL)
      i = (int) class_double (e1);

    if (i > 0)
      claw1_nv[0] = i;

    ent_Clean (e1);
  }

  return ent_Create_Rlab_Success();
}

//
// basic 1-d solver
//
#undef THIS_SOLVER
#define THIS_SOLVER "claw1"
Ent *
ent_claw1_basic (int nargs, Datum args[])
{
  Ent *e4=0, *e5=0;

  MDR *Q0, *tm, *aux, *work, *q1;
  int i, j, mx, maux = 0, mcapa = 0,  mwork;
  int mbc = claw1_mbc, info, *mthlim, mthbc[2];
  double xlower, t1, t2, dx, dtv[5];

  FILE *rlab_stderr = (!RLAB_STDERR_DS) ? stderr : RLAB_STDERR_DS;

  //
  // Initialize global variables for user function calls
  //
  claw1_method[1] = claw1_order;
  claw1_method[5] = 0;

  //
  // Load and Check arguments.
  //
  if (nargs < 5)
  {
    fprintf (rlab_stderr,
             THIS_SOLVER ": Basic CLAWPACK solver for 1-D problems.\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": Format:\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ":   y=claw1(/K/,R,/S/,Q0,t),\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": where K=function(x), if given, is a capacity function,\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": R=function(ql,qr) is a riemann solver for the flux function, while\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": S=function(t,x,q), if given, is a source function, matrix Q0=[x,q0]\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": contains mesh vector 'x' and the initial condition matrix 'q0'.\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": The matrix t=[t1,t2] gives the time integration interval.\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": The function returns  y=[x,q(t1)].\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": See also claw1params.\n");
    rerror  (THIS_SOLVER ": " RLAB_ERROR_AT_LEAST_FIVE_ARGS_REQUIRED);
  }

  //
  // Get capacity function K=K(x) or K=1
  //
  if (kname)
    ent_Clean (kname);
  kname = bltin_get_ent(args[0]);
  if (isfuncent(kname))
    mcapa = 1;

  //
  // Get Riemann solver function R = R(x,q,t)
  //
  if (rname)
    ent_Clean (rname);
  rname = bltin_get_ent(args[1]);
  if (!isfuncent(rname))
    rerror (THIS_SOLVER ": " RLAB_ERROR_ARG2_FUNC_VAR "\n");

  //
  // Get source function S=S(x,q,t) or if omitted do not have a source
  //
  claw1_method[4]=0;
  if (sname)
    ent_Clean (sname);
  sname = bltin_get_ent(args[2]);
  if (isfuncent(sname))
  {
     // S given as a function
    claw1_method[4] = claw1_m5_fracsplitmethod;
  }

  //
  // get mesh with initial conditions
  //
  e4 = bltin_get_ent (args[3]);
  if (ent_type (e4) != MATRIX_DENSE_REAL)
    rerror ("claw1: third entry (mesh/IC) should be a real matrix!");
  Q0 = mdr_Copy ( class_matrix_real (e4) );

  //
  // find the coordinates of x(i-1/2) walls including the ghost cells ('mbc' on each side)
  //
  mx = (int) MNR (Q0);
  dx = Mdr1 (Q0, 2, 1) - Mdr1 (Q0, 1, 1);
  xlower = Mdr1 (Q0, 1, 1) - 0.5 * dx;
  x = mdr_Create (mx + 2 * mbc, 1);
  for (i = 1; i <= mx; i++)
    Mdr1 (x, i + mbc, 1) = Mdr1 (Q0, i, 1);
  for (i = 1; i <= mbc; i++)
  {
    // x does not vary in the ghost cells, only Q
    Mdr1 (x, i, 1) = Mdr1 (x, 1 + mbc, 1);
    Mdr1 (x, mx + mbc + i, 1) = Mdr1 (x, mx + mbc, 1);
  }

  //
  // initial conditions with the offset of mbc rows at the top and at the bottom
  //
  claw1_meqn = MNC (Q0) - 1;
  claw1_nrow = mx + 2 * mbc;
  q1 = mdr_Create (claw1_nrow, claw1_meqn);
  for (j = 1; j <= claw1_meqn; j++)
    for (i = 1; i <= mx; i++)
      Mdr1 (q1, i + mbc, j) = Mdr1 (Q0, i, j + 1);
  claw1_mwav = claw1_meqn;	// number of waves

  //
  // get time interval
  //
  e5 = bltin_get_ent (args[4]);
  if (ent_type (e5) != MATRIX_DENSE_REAL)
    rerror ("claw1: 4th entry (time) should be a row-vector!");
  tm = class_matrix_real (e5);
  if (MNR (tm) * MNC (tm) != 2)
    rerror ("claw1: 4th entry (time) should be [t1,t2]!");
  t1 = MdrV1 (tm, 1);
  t2 = MdrV1 (tm, 2);

  // Set up ENTITIES for user-function.
  // x, a global variable for
  x_ent = ent_Assign_Rlab_MDR (x);
  ent_IncRef (x_ent);

  // x1, a single value of 'x'
  x1     = mdr_CreateEmpty (1, 1);
  x1_ent = ent_Assign_Rlab_MDR (x1);
  ent_IncRef (x1_ent);

  // q=q(x,t), a matrix
  q = mdr_CreateEmpty (claw1_nrow, claw1_meqn);
  q_ent = ent_Assign_Rlab_MDR (q);
  ent_IncRef (q_ent);

  // ql, a vector
  qL = mdr_Create (claw1_meqn, 1);
  qL_ent = ent_Assign_Rlab_MDR (qL);
  ent_IncRef (qL_ent);

  // qr, a vector
  qR = mdr_Create (claw1_meqn, 1);
  qR_ent = ent_Assign_Rlab_MDR (qR);
  ent_IncRef (qR_ent);

  // time, a scalar
  t = mdr_CreateEmpty (1, 1);
  t_ent = ent_Assign_Rlab_MDR (t);
  ent_IncRef (t_ent);

  //
  // needed only if K(x) is given
  //
  aux = mdr_Copy (q1);
  if (mcapa == 1)
  {
    // if K given than it is in the last column of aux!
    // however it needs to be a global variable for time integrator ODE
    claw1_retk = (MDR *) claw1k ();
    aux = mdr_Append (aux, claw1_retk);
    claw1_method[5] = MNC (aux);
  }
  claw1_method[6] = MNC (aux);
  maux = MNC (aux);


  //
  // initialize limiter per wave to a maximum of 10 waves
  //
  mthlim = (int *) GC_malloc (claw1_mwav * sizeof (int));
  for (i = 0; i < claw1_mwav; i++)
    mthlim[i] = claw1_mthlim[MIN(i, 9)];
  //
  // initialize types of boundary conditions on L (0) and R (1)
  //
  mthbc[0] = claw1_mthbc[0];
  mthbc[1] = claw1_mthbc[1];

  // cannot have rigid wall bc for a single equation (need at least two)
  if (claw1_meqn==1)
  {
    if (mthbc[0]==3)
      mthbc[0] = 1;
    if (mthbc[1]==3)
      mthbc[1] = 1;
  }

  //
  // dtv
  //
  dtv[0] = claw1_mindt;
  dtv[1] = claw1_maxdt;

  // working arrays
  mwork =
      (mx + 2 * mbc) * (2 + 4 * (claw1_meqn) + (claw1_mwav) +
      (claw1_meqn) * (claw1_mwav));
  work = mdr_Create (mwork, 1);
  //
  // Finally, call the solver
  //
  CLAW1RL (&claw1_meqn, &claw1_mwav, &mbc, &maux, &mwork, mthlim,
            MDRPTR(q1), MDRPTR(work), MDRPTR(aux),
            &xlower, &dx, &mx,
            &t1, &t2,
            claw1_method, claw1_cflv, claw1_nv,
            mthbc,
            &info, dtv);

  for (j = 1; j <= claw1_meqn; j++)
    for (i = 1; i <= mx; i++)
      Mdr1 (Q0, i, j + 1) = Mdr1 (q1, i + mbc, j);

  // Clean Up
  mdr_Destroy (work);
  mdr_Destroy (q1);
  mdr_Destroy (aux);

  if (mcapa == 1)
    mdr_Destroy (claw1_retk);

  // x_ent
  MDPTR (x)=0;
  ent_DecRef (x_ent);
  ent_Destroy (x_ent);

  // q_ent
  MDPTR (q)=0;
  ent_DecRef (q_ent);
  ent_Destroy (q_ent);

  // x1_ent
  MDPTR(x1) = 0;
  ent_DecRef (x1_ent);
  ent_Destroy (x1_ent);

  // qL
  ent_DecRef (qL_ent);
  ent_Destroy (qL_ent);

  // qR
  ent_DecRef (qR_ent);
  ent_Destroy (qR_ent);

  // t_ent
  MDPTR(t) = 0;
  ent_DecRef (t_ent);
  ent_Destroy (t_ent);

  GC_FREE (mthlim);

  ent_Clean (e4);
  ent_Clean (e5);

  // cleanup function entities
  ent_Clean (rname);
  ent_Clean (sname);
  ent_Clean (bclname);
  ent_Clean (bcrname);
  ent_Clean (b4name);
  ent_Clean (kname);

  //
  // return the result as a list
  //
  return ent_Assign_Rlab_MDR(Q0);
}



//
// b4step1 function: if something else (which is not a part of the claw1
// has to be calculated as well.
//    function(t,x,q)
//
int
B4STEP1 (int * maxmx, int * mbc, int * mx, int * meqn, double *qc,
         double * xlower, double * dx, double * tc, double * dt,
         int * maux, double * aux)
{
  int i, j;
  Ent *rent = 0;
  MDR *retm = 0;
  FILE *rlab_stderr = (!RLAB_STDERR_DS) ? stderr : RLAB_STDERR_DS;

  if (claw1_b4step1 == 0)
    return 1;

  //
  // time
  //
  MDPTR(t) = (void *) tc;

  //
  // x is already here x_ent
  //

  //
  // assign to MDR *q the content of qc (points to MDR *q1)
  //
  MDPTR(q) = (void *) qc;

  rent = ent_call_rlab_script_3args(b4name, t_ent, x_ent, q_ent);

  if (ent_type(rent)!=MATRIX_DENSE_REAL)
    rerror (THIS_SOLVER ": " RLAB_ERROR_RHS_FUNC_MUST_RETURN_MDR);

  retm = ent_data (rent);

  if (MNR (retm) != claw1_nrow || MNC (retm) != claw1_meqn)
  {
    fprintf (rlab_stderr,
             THIS_SOLVER ": The intermediate step function  b4s1=b4s1(t,x,q)  has to return\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": result such that  dim(b4s1) = dim(q).\n");
    rerror  (THIS_SOLVER ": Improper function b4s1=b4s1(t,x,q)!");
  }
  //
  // qc <- b4step1(t,x,qc)
  //
  for (i = 1; i <= claw1_nrow; i++)
    for (j = 1; j <= claw1_meqn; j++)
      qc[(j - 1) * claw1_nrow + (i - 1)] = Mdr1 (retm, i, j);

  ent_Clean (rent);
  return 1;
}


//
// left boundary condition
//
int
BC1LEFT (int * maxmx, int * meqn, int * mbc, int * mx,
         double * xlower, double * dx, double *qc,
         int *maux, double * aux, double * t1, double * dt,
         int *mthbc)
{
  int i, j;
  Ent *rent = 0;
  MDR *retm = 0;
  FILE *rlab_stderr = (!RLAB_STDERR_DS) ? stderr : RLAB_STDERR_DS;

  //
  // time
  //
  MDPTR(t) = (void *) t1;

  //
  // x is already here x_ent
  //

  //
  // assign to MDR *q the content of qc (points to MDR *q1)
  //
  MDPTR(q) = (void *) qc;

  rent = ent_call_rlab_script_3args(bclname, t_ent, x_ent, q_ent);

  if (ent_type(rent)!=MATRIX_DENSE_REAL)
    rerror (THIS_SOLVER ": " RLAB_ERROR_RHS_FUNC_MUST_RETURN_MDR);

  retm = ent_data (rent);

  if (MNR (retm) != claw1_nrow || MNC (retm) != claw1_meqn)
  {
    fprintf (rlab_stderr,
             THIS_SOLVER ": The left boundary condition function  bcL=bcL(t,x,q) has to return\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": result such that  dim(bcL) = dim(q).\n");
    rerror  (THIS_SOLVER ": Improper function bcL=bcL(t,x,q)!");
  }

  //
  // qc -> bcL(t,x,qc) only top mbc layers are copied
  //
  for (i = 1; i <= claw1_mbc; i++)
    for (j = 1; j <= claw1_meqn; j++)
      qc[(j - 1) * claw1_nrow + (i - 1)] = Mdr1 (retm, i, j);

  ent_Clean (rent);
  return 1;
}

//
// right boundary condition
//
int
BC1RIGHT (int * maxmx, int * meqn, int * mbc, int * mx,
          double * xlower, double * dx, double *qc,
          int *maux, double * aux, double * t1, double * dt,
          int *mthbc)
{
  int i, j;
  Ent *rent = 0;
  MDR *retm = 0;
  FILE *rlab_stderr = (!RLAB_STDERR_DS) ? stderr : RLAB_STDERR_DS;

  //
  // time
  //
  MDPTR(t) = (void *) t1;

  //
  // x is already here x_ent
  //

  //
  // assign to MDR *q the content of qc (points to MDR *q1)
  //
  MDPTR(q) = (void *) qc;

  rent = ent_call_rlab_script_3args(bcrname, t_ent, x_ent, q_ent);

  if (ent_type(rent)!=MATRIX_DENSE_REAL)
    rerror (THIS_SOLVER ": " RLAB_ERROR_RHS_FUNC_MUST_RETURN_MDR);

  retm = ent_data (rent);

  if (MNR (retm) != claw1_nrow || MNC (retm) != claw1_meqn)
  {
    fprintf (rlab_stderr,
             THIS_SOLVER ": The right boundary condition function  bcR=bcR(t,x,q) has to return\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": result such that  dim(bcR) = dim(q).\n");
    rerror  (THIS_SOLVER ": Improper function bcL=bcL(t,x,q)!");
  }

  //
  // qc -> bcR(t,x,qc) only last mbc layers are copied
  //
  for (i = claw1_nrow; i > claw1_nrow - claw1_mbc; i--)
    for (j = 1; j <= claw1_meqn; j++)
      qc[(j - 1) * claw1_nrow + (i - 1)] = Mdr1 (retm, i, j);

  ent_Clean (rent);
  return 1;
}


//
// Source function: ac <-> q' = S(t,x,q)
// or               ac <-> q' = S(t,x,q)/K(x)
// so 'claw1srk_' does not have to worry about K(x)
//
int
CLAW1S (double *t1, double *qc, double *ac)
{
  int i, j;
  Ent *rent = 0;
  MDR *retm = 0;
  FILE *rlab_stderr = (!RLAB_STDERR_DS) ? stderr : RLAB_STDERR_DS;

  //
  // time
  //
  MDPTR(t) = (void *) t1;

  //
  // x is already here x_ent
  //

  //
  // assign to MDR *q the content of qc (points to MDR *q1)
  //
  MDPTR(q) = (void *) qc;

  rent = ent_call_rlab_script_3args(sname, t_ent, x_ent, q_ent);

  if (ent_type(rent)!=MATRIX_DENSE_REAL)
    rerror (THIS_SOLVER ": " RLAB_ERROR_RHS_FUNC_MUST_RETURN_MDR);

  retm = ent_data (rent);

  if (MNR (retm) != claw1_nrow || MNC (retm) != claw1_meqn)
  {
    fprintf (rlab_stderr,
             THIS_SOLVER ": The source function  S=S(t,x,q) has to return\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": result such that  dim(S) = dim(q).\n");
    rerror  (THIS_SOLVER ": Improper function S=S(t,x,q)!");
  }

  //
  // ac = S(t,x,qc)
  //
  if (claw1_method[5] == 0)
  {
    for (i = 1; i <= claw1_nrow; i++)
      for (j = 1; j <= claw1_meqn; j++)
        ac[(j - 1) * claw1_nrow + (i - 1)] = Mdr1 (retm, i, j);
  }
  else
  {
    //
    // ac = S(t,x,qc)./K(x)
    //
    for (i = 1; i <= claw1_nrow; i++)
      for (j = 1; j <= claw1_meqn; j++)
        ac[(j - 1) * claw1_nrow + (i - 1)] =
            Mdr1 (retm, i, j) / MdrV1 (claw1_retk, i);
  }

  ent_Clean (rent);
  return 1;
}

//
// Integrate q(x)_t = S(t,x,q) for t=t1..t2 with q(t1)=qc.
// On return  qc=q(t2)
//
int
SRC1 (int * maxmx, int * meqn, int * mbc, int * mx,
      double * xlower, double * dx, double * qc,
      int * maux, double * aux, double * t, double * dt,
      int * method)
{
  double t1 = *t, t2 = *t + *dt;
  //
  // is there a source term?
  //
  if (method[5 -1] == 0) return 1;

  //
  // integrate in time
  //
  if (claw1_tsrc == 0)
  {
    //
    // claw1_tsrc = 0,  adams' method explicite/implicite. EXCELLENT !
    //
    int iwork[5], neqn = claw1_nrow * claw1_meqn, iflag = -1;
    MDR *work = mdr_Create (100 + 21 * neqn, 1);
    ODE (CLAW1S, &neqn, qc, &t1, &t2,
         &claw1_relerr, &claw1_abserr, &iflag, MDRPTR(work), iwork);
    mdr_Destroy (work);
  }
  else
  {
    //
    // claw1_tsrc=1, rk2, explicite. POOR !
    //
    double t1half = t1 + 0.5 * (*dt);
    int i, j;
    // aux = S(t+dt/2,x,qc)
    CLAW1S (&t1, qc, aux);
    // aux = q + dt/2*aux
    for (i = 1; i <= claw1_nrow; i++)
      for (j = 1; j <= claw1_meqn; j++)
        aux[(j - 1) * claw1_nrow + (i - 1)] =
            qc[(j - 1) * claw1_nrow + (i - 1)] +
            0.5 * (*dt) * aux[(j - 1) * claw1_nrow + (i - 1)];
    // aux = S(t+dt/2,x,aux)
    CLAW1S (&t1half, aux, aux);
    // q = q + dt*aux
    for (i = 1; i <= claw1_nrow; i++)
      for (j = 1; j <= claw1_meqn; j++)
        qc[(j - 1) * claw1_nrow + (i - 1)] +=
            (*dt) * aux[(j - 1) * claw1_nrow + (i - 1)];
  }
  return 1;
}


//
// Riemann solver: called for each point of the mesh.
// has to be wraped in fortran code RP1 because it needs
// variables from common block /comxt/
//
int
CLAW1R (double *t1, double *ql, double *qr, double *waves, double *speed,
         double *amdq, double *apdq)
{
  int i, j, k;
  double c, w;
  Ent *rent = 0;
  MDR *retm = 0;

  //
  // time
  //
  MDPTR(t) = (void *) t1;

  //
  // x is already here x_ent
  //

  for (i = 2; i <= claw1_nrow; i++)
  {
    // x1
    //Mdr1 (x1, 1, 1) = Mdr1 (x, i, 1);
    MDPTR(x1) = (void *) &MdrV0(x,i-1);

    // qL(i, :)
    for (j = 1; j <= claw1_meqn; j++)
       Mdr1 (qL, j, 1) = ql[(j - 1) * claw1_nrow + (i - 1)];

    // qR(i-1; :)
    for (j = 1; j <= claw1_meqn; j++)
      Mdr1 (qR, j, 1) = qr[(j - 1) * claw1_nrow + (i - 2)];

    rent = ent_call_rlab_script_4args(rname, t_ent, x_ent, qL_ent, qR_ent);

    if (ent_type(rent)!=MATRIX_DENSE_REAL)
      rerror (THIS_SOLVER ": " RLAB_ERROR_RHS_FUNC_MUST_RETURN_MDR);

    retm = ent_data (rent);

    //
    // Return options:
    // (1) retm = [waves; speed] or
    // (2) retm = [waves,amdq,apdq; [speed,0,0]]
    //
    // reset the waves for (1)
    if (MNC (retm) == claw1_mwav)
    {
      for (k = 1; k <= claw1_meqn; k++)
      {
        apdq[i - 1 + claw1_nrow * (k - 1)] = 0;
        amdq[i - 1 + claw1_nrow * (k - 1)] = 0;
      }
    }
    for (j = 1; j <= claw1_mwav; j++)
    {
      c = Mdr1 (retm, claw1_meqn + 1, j);
      //speed(i,j) = c
      speed[(j - 1) * claw1_nrow + (i - 1)] = c;
      for (k = 1; k <= claw1_meqn; k++)
      {
        //waves(i,:,j) = Mdr1(retm,:,j);
        w = Mdr1 (retm, k, j);
        waves[i - 1 + claw1_nrow * (k - 1) + claw1_meqn * claw1_nrow * (j - 1)] = w;
        if (MNC (retm) == claw1_mwav)
        {
          // (1) the fluxes are calculated here
          if (c > 0)
            apdq[i - 1 + claw1_nrow * (k - 1)] += c * w;
          else
            amdq[i - 1 + claw1_nrow * (k - 1)] += c * w;
        }
      }
    }
    //
    // (2) user calculated the fluxes
    // retm = [waves,amdq,apdq; [speed,0,0]]
    //
    if (MNC (retm) == claw1_mwav + 2)
    {
      for (k = 1; k <= claw1_meqn; k++)
      {
        amdq[i - 1 + claw1_nrow * (k - 1)] = Mdr1 (retm, k, claw1_mwav + 1);
        apdq[i - 1 + claw1_nrow * (k - 1)] = Mdr1 (retm, k, claw1_mwav + 2);
      }
    }

    ent_Clean (rent);
  }

  return 1;
}

//
// calculate user specified capacity function once
//
static MDR * claw1k (void)
{
  Ent *rent = 0;
  MDR *retm = 0, *rval = 0;

  //
  // x is already here x_ent
  //

  rent = ent_call_rlab_script_1arg (kname, x_ent);

  if (ent_type(rent)!=MATRIX_DENSE_REAL)
    rerror (THIS_SOLVER ": " RLAB_ERROR_RHS_FUNC_MUST_RETURN_MDR);

  retm = ent_data (rent);

  rval = mdr_Float_BF(retm);
  ent_Clean (rent);
  return rval;
}
