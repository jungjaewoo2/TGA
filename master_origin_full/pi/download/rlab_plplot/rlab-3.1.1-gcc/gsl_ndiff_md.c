// Copyright (C) 2003-2004 Marijan Kostrun
//   part of rlab+4linux project on rlabplus.sourceforge.net
//
// GSL Science Library - numerical differentiation in multidimensions and gradient
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

// gsl headers
// shared object
#include <gsl/gsl_mode.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_precision.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_machine.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_diff.h>


// standard libraries
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>

// naming convention for the solver parameters
#include "rfileio.h"
#include "rlab_solver_parameters_names.h"

static double
ndiff_gslrlab_f (double xj, void *pH);
static MDR *mdr_gsl_diff_central (double x, void *pH, double *abserr, int nef);

// other variables
static int neq, nex, nef;
static MDR *my;
static Ent *my_ent, *pent;
static Ent *fname;

//
// ndiffs: numerical differentiation of a vector function
//
#undef  THIS_SOLVER
#define THIS_SOLVER "ndiffs"
Ent *
ent_ndiffs (int nargs, Datum args[])
{
  // call parameters:
  //  derivative      d1
  //  parameters      e2
  //  x0(initial)     e3
  double eabs, xj;
  int i, j;
  Ent *e2=0, *e3=0, *rent=0;
  MDR *x0=0, *w=0, *retm=0, *dummy=0;

  FILE *rlab_stderr = (!RLAB_STDERR_DS) ? stderr : RLAB_STDERR_DS;

  //
  // Load and Check arguments.
  //
  if (nargs < 3)
  {
    fprintf (rlab_stderr, THIS_SOLVER ": Numerical differentiation of a vector function of a\n");
    fprintf (rlab_stderr, THIS_SOLVER ": vector variable.\n");
    fprintf (rlab_stderr, THIS_SOLVER ": Format:\n");
    fprintf (rlab_stderr, THIS_SOLVER ":   y = ndiffs(f,/p/,x),\n");
    fprintf (rlab_stderr, THIS_SOLVER ": where  f=function(x/,p/){..}, 'p' is its parameter array,\n");
    fprintf (rlab_stderr, THIS_SOLVER ": and 'x' is column-vector at which the derivative is being\n");
    fprintf (rlab_stderr, THIS_SOLVER ": calculated.\n");
    rerror (THIS_SOLVER ": requires three arguments");
  }

  //
  // f = f(x/,p/)
  //
  fname = bltin_get_ent(args[0]);
  if (!isfuncent(fname))
    rerror (THIS_SOLVER ": " RLAB_ERROR_ARG1_FUNC_VAR "\n");

  //
  // parameter entity
  //
  e2 = bltin_get_ent (args[1]);
  pent = 0;
  if (ent_type (e2) != UNDEF)
    pent = ent_Copy (e2);

  // Get x0 where we look for derivatives
  e3 = bltin_get_ent (args[2]);
  if (ent_type(e3)!=MATRIX_DENSE_REAL)
    rerror (THIS_SOLVER ": Second argument 'X0' must be real vector");
  x0 = class_matrix_real (e3);
  if (!x0)
    rerror (THIS_SOLVER ": Second argument 'X0' must be real vector");
  if (x0->nrow!= 1 && x0->ncol!= 1)
    rerror (THIS_SOLVER ": Second argument 'X0' must be real vector");
  nex = x0->ncol * x0->nrow;

  //
  // Set up ENTITIES for user-function.
  // Inc the reference count once, for belonging
  // array x0:
  my = mdr_Float_BF (x0);
  my_ent = ent_Assign_Rlab_MDR (my);
  ent_IncRef (my_ent);

  // getting the dimension of function, call it once and figure it out
  // looks annoying but one has to keep track or errors
  if (pent)
    rent = ent_call_rlab_script_2args(fname, my_ent, pent);
  else
    rent = ent_call_rlab_script_1arg (fname, my_ent);
  if (ent_type (rent) != MATRIX_DENSE_REAL)
    rerror (THIS_SOLVER ": " RLAB_ERROR_RHS_FUNC_MUST_RETURN_MDR);
  retm = ent_data (rent);
  nef = MNR (retm);
  ent_Clean (rent);
  //
  //

  // define
  gsl_set_error_handler_off ();

  w = mdr_Create (nef, nex);

  // define a parameter array
  double pH[2] ={0};

  // run over the indices of the jacobian matrix
  for (j = 1; j <= nex; j++)
  {
    pH[1] = j;          // x_j
    xj = MdrV1 (x0, j);
    dummy = mdr_gsl_diff_central (xj, pH, &eabs, nef);
    MdrV1 (my, j) = MdrV1 (x0, j);
    for (i = 1; i <= nef; i++)
      Mdr1 (w, i, j) = Mdr1 (dummy, i, 1);
  }

  // Clean Up
  ent_DecRef (my_ent);
  ent_Destroy (my_ent);

  ent_Clean (pent);
  ent_Clean (e2);
  ent_Clean (e3);
  ent_Clean (fname);

  return ent_Assign_Rlab_MDR(w);
}


#undef  THIS_SOLVER
#define THIS_SOLVER "ndiv"
Ent *
ent_div (int nargs, Datum args[])
{
  // call parameters:
  //  derivative      d1
  //  parameters      e2
  //  x0(initial)     e3
  double r, eabs, xj, dw;
  int i;
  Ent *e2=0, *X0;
  MDR *x0;
  X0 = 0;

  FILE *rlab_stderr = (!RLAB_STDERR_DS) ? stderr : RLAB_STDERR_DS;

  //
  // Load and Check arguments.
  //
  if (nargs < 3)
  {
    fprintf (rlab_stderr,
             THIS_SOLVER ": Numerical divergence of a vector function of a vector\n");
    fprintf (rlab_stderr, THIS_SOLVER ": variable.\n");
    fprintf (rlab_stderr, THIS_SOLVER ": Format:\n");
    fprintf (rlab_stderr, THIS_SOLVER ":   y = ndiffs(f,/p/,x),\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": where  f=function(x/,p/){..}, 'p' is its parameter array,\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": and 'x' is column-vector at which the derivative is being\n");
    fprintf (rlab_stderr, THIS_SOLVER ": calculated.\n");
    rerror (THIS_SOLVER " requires at three arguments");
  }

  //
  // f = f(x/,p/)
  //
  fname = bltin_get_ent(args[0]);
  if (!isfuncent(fname))
    rerror (THIS_SOLVER ": " RLAB_ERROR_ARG1_FUNC_VAR "\n");

  //
  // parameter entity
  //
  e2 = bltin_get_ent (args[1]);
  pent = 0;
  if (ent_type (e2) != UNDEF)
    pent = ent_Copy (e2);

  // Get x0 where we look for derivatives
  X0 = bltin_get_ent (args[2]);
  if (ent_type (X0) != MATRIX_DENSE_REAL)
    rerror ("div: Argument 'x0' must be real vector!");
  x0 = class_matrix_real (X0);
  if (!x0)
    rerror ("div: Argument 'x0' must be real vector!");
  // Extract neq from ystart
  if (MNC (x0) != 1 && MNR (x0)!=1)
    rerror ("div: Argument 'x0' must be real vector!");
  neq = MNR (x0) * MNC(x0);

  //
  // Set up ENTITIES for user-function.
  // Inc the reference count once, for belonging
  // array x0:
  my = mdr_Create (neq, 1);
  my_ent = ent_Assign_Rlab_MDR (my);
  ent_IncRef (my_ent);

  gsl_set_error_handler_off ();

  // define a parameter array
  double pH[2] ={0};

  // define a parameter array
  pH[0] = 0;            // f_i
  pH[1] = 0;            // x_j

  gsl_function F;
  F.function = &ndiff_gslrlab_f;
  F.params = pH;

  // run over the indices to find df_i/dx_i
  dw = 0.0;
  for (i = 0; i < neq; i++)
  {
    pH[0] = i;          // f_i
    pH[1] = i;          // x_j
    // calls to &F change value of xj. We want to make sure that
    // once derivative is found, xj has its original value.
    xj = MdrV0 (x0, i);
    gsl_diff_central (&F, xj, &r, &eabs);
    MdrV0 (my, i) = MdrV0 (x0, i);
    dw += r;
  }

  // Clean Up
  ent_DecRef (my_ent);
  ent_Destroy (my_ent);

  ent_Clean (pent);
  ent_Clean (e2);
  ent_Clean (X0);
  ent_Clean (fname);

  return ent_Create_Rlab_Double(dw);
}


//
// jacobian function modified to get f_i(x_j);
//
static double
ndiff_gslrlab_f (double xj, void *pH)
{
  int i, j;

  Ent *rent = 0;
  MDR *retm = 0;

  i = *((double *) pH);         // ph[0]=i as in f_i
  j = *((double *) pH + 1);     // ph[1]=j as in x_j

  double rval;

  MdrV0 (my, j) = xj;

  if (pent)
    rent = ent_call_rlab_script_2args(fname, my_ent, pent);
  else
    rent = ent_call_rlab_script_1arg (fname, my_ent);

  if (ent_type(rent)!=MATRIX_DENSE_REAL)
    rerror (THIS_SOLVER ": " RLAB_ERROR_RHS_FUNC_MUST_RETURN_MDR);

  retm = ent_data (rent);

  if (SIZE(retm) != neq)
    rerror (THIS_SOLVER ": " RLAB_ERROR_RHS_FUNC_INCORRECT_DIM);

  rval = MdrV0 (retm, i);
  ent_Clean (rent);
  return rval;
}

//
// function for derivatives using user's column-vector function
//
static MDR *
mdr_rks_func (double xj, void *pH)
{
  int j;
  double xjold;
  Ent *rent=0;
  MDR *retm=0, *rval=0;

  j = *((double *) pH + 1);     // ph[1]=j as in x_j

  // 'my' already contains x0, given as 'xj' is its j-th component
  // 'p' is already there and it doesn't change
  xjold = MdrV1 (my, j);
  MdrV1 (my, j) = xj;

  if (pent)
    rent = ent_call_rlab_script_2args(fname, my_ent, pent);
  else
    rent = ent_call_rlab_script_1arg (fname, my_ent);

  MdrV1 (my, j) = xjold;

  if (ent_type (rent) != MATRIX_DENSE_REAL)
    rerror (THIS_SOLVER ": " RLAB_ERROR_JACOBIAN_FUNC_MUST_RETURN_MDR);

  retm = ent_data (rent);

  if (SIZE(retm) != nef)
    rerror (THIS_SOLVER ": " RLAB_ERROR_JACOBIAN_FUNC_INCORRECT_DIM);

  rval = mdr_Copy(retm);
  MNC(rval) = 1;
  MNR(rval) = SIZE(retm);

  ent_Clean (rent);
  return rval;
}


MDR *
mdr_gsl_diff_central (double x, void *pH, double *abserr, int nef)
{
  /* Construct a divided difference table with a fairly large step
     size to get a very rough estimate of f'''.  Use this to estimate
     the step size which will minimize the error in calculating f'. */

  int i, j, k;
  double h = GSL_SQRT_DBL_EPSILON;
  double a[4], a3;
  MDR *d = mdr_Create (nef, 4), *dummy = mdr_Create (nef, 1), *retm =
    mdr_Create (nef, 1);

  // Algorithm based on description on pg. 204 of Conte and de Boor
  // (CdB) - coefficients of Newton form of polynomial of degree 3.
  for (i = 0; i < 4; i++)
  {
    a[i] = x + (i - 2.0) * h;
    dummy = mdr_rks_func (a[i], pH);
    for (j = 1; j <= nef; j++)
    {
      Mdr1 (d, j, i + 1) = Mdr1 (dummy, j, 1);
    }
  }

  for (k = 1; k < 5; k++)
  {
    for (i = 0; i < 4 - k; i++)
    {
      for (j = 1; j <= nef; j++)
      {
        Mdr1 (d, j, i + 1) =
          (Mdr1 (d, j, i + 1 + 1) - Mdr1 (d, j, i + 1)) / (a[i + k] - a[i]);
      }
    }
  }

  // Adapt procedure described on pg. 282 of CdB to find best
  // value of step size. Use kostrun's fudge:
  //  find a maximum of each column, then sum them up, then find abs of the sum
  // a3 = fabs (Mdr1 (mdr_Sum_BF (mdr_Max1 (mdr_Abs (d))), 1, 1));
  MDR *d_abs = mdr_Abs (d);
  MDR *d_minmax_i = mdr_MinMax1_ValIdx(d_abs,2,1,0,0, NULL, NULL);
  MDR *d_sum = mdr_Sum_BF ( d_minmax_i, NULL );

  a3 = fabs (Mdr1 (d_sum, 1, 1));

  mdr_Destroy(d_abs);
  mdr_Destroy(d_minmax_i);
  mdr_Destroy(d_sum);

  if (a3 < 100.0 * GSL_SQRT_DBL_EPSILON)
    a3 = 100.0 * GSL_SQRT_DBL_EPSILON;

  h = pow (GSL_SQRT_DBL_EPSILON / (2.0 * a3), 1.0 / 3.0);

  if (h > 100.0 * GSL_SQRT_DBL_EPSILON)
    h = 100.0 * GSL_SQRT_DBL_EPSILON;

  retm = mdr_Subtract (mdr_rks_func (x + h, pH), mdr_rks_func (x - h, pH));
  for (j = 1; j <= nef; j++)
    Mdr1 (retm, j, 1) /= (2.0 * h);

  *abserr = fabs (100.0 * a3 * h * h);

  mdr_Destroy(d);
  mdr_Destroy(dummy);

  return retm;
}
