// Copyright (C) 2003-2004 Marijan Kostrun
//   part of rlab+4linux project on rlabplus.sourceforge.net
//
// GSL Science Library - leastsquares solver
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

// shared object
#include <gsl/gsl_version.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_block.h>
#include <gsl/gsl_mode.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_precision.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_machine.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_fit.h>
#include <gsl/gsl_multifit.h>

#if (GSL_MAJOR_VERSION < 2)
#include <gsl/gsl_multifit_nlin.h>
#else
#include <gsl/gsl_multifit_nlinear.h>
#endif

// standard libraries
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <time.h>

#include "minpack.h"
#include "slatec.h"

// naming convention for the solver parameters
#include "rfileio.h"
#include "rlab_solver_parameters_names.h"
#include "gsl_rlab.h"

// other variables
static Ent *data_ent, *params_ent;
static MDR *xdata, *params, *ydata, *wdata, *xentry;
static Ent *fname, *dfname;

static int nparams, xr, xc;

// functions in this library
static int
ls_gslrlab_f (const gsl_vector * x, void *dummy, gsl_vector * f);
static int
ls_gslrlab_df (const gsl_vector * x, void *dummy, gsl_matrix * J);
#if (GSL_MAJOR_VERSION < 2)
static int
ls_gslrlab_fdf (const gsl_vector * x, void *dummy, gsl_vector * f,
                gsl_matrix * J)
{
  ls_gslrlab_f (x, params, f);
  ls_gslrlab_df (x, params, J);
  return GSL_SUCCESS;
}
#endif

#undef  THIS_SOLVER
#define THIS_SOLVER "polyfit"
Ent *
ent_polyfit (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *eo=0;
  MDR *y=0, *x=0, *w=0, *poln=0, *r=0, *a=0;

  int i, np, degn, maxdeg = 16, mdegn, ierr, idummy;
  double eps = 0.0, dzero = 0.0, f = 0.0, ddummy;

  ListNode * node;

  //
  // Load and Check arguments.
  //
  if (nargs != 2 && nargs != 3)
    rerror ("polyfit: requires two or three arguments");

  fname=0;
  dfname=0;

  //
  // Get observations: y, or <<val;wgt>>
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) == MATRIX_DENSE_REAL)
    y = class_matrix_real (e1);
  else if (ent_type (e1) == BTREE)
  {
    // val from <<val;wgt>>
    ListNode *node;
    node = btree_FindNode (ent_data (e1), "val");
    if (node != 0)
      y = class_matrix_real (var_ent (node));
    // wgt from <<val;wgt>>
    node = btree_FindNode (ent_data (e1), "wgt");
    if (node != 0)
    {
      w = class_matrix_real (var_ent (node));
      if (SIZE(w) != SIZE(y))
        rerror ("'val' and 'wgt' must have the same size!");
      for (i = 0; i < SIZE(w); i++)
        if (MdrV0 (w, i) < 0)
          rerror ("'wgt' must be positive !");
    }
  }
  else
    rerror ("polyfit: First argument 'Y' must be real matrix or list <<val;wgt>>!");

  if (!EQVECT(y))
    rerror ("polyfit: dependent variable 'y' has to be a real vector!");

  //
  // Get predictor variable: x, or <<val;wgt>>
  //
  e2 = bltin_get_ent (args[1]);
  if (ent_type (e2) != MATRIX_DENSE_REAL)
    rerror ("polyfit: 'X' has to be a real vector");
  x = class_matrix_real (e2);

  // check x and y
  np = SIZE (x);
  if (!EQVECT(x) || !EQVECT(y))
    rerror ("polyfit: dependent variable 'y' has to be a real vector");
  if (np != SIZE(y))
    rerror ("polyfit: 'Y' and 'X' must have the same number of rows");

  //
  // options for the solver
  //
  if (nargs == 3)
  {
    eo = bltin_get_ent (args[2]);
    if (ent_type (eo) == BTREE)
    {
      // maxdeg
      node = btree_FindNode (ent_data (eo), RLAB_NAME_GSL_LS_MAXDEG);
      if (node != 0)
      {
        idummy = (int) class_double (var_ent (node));
        if (idummy >= 1 && idummy <= np - 2)
          maxdeg = idummy;
      }
      // f
      node = btree_FindNode (ent_data (eo), RLAB_NAME_GSL_LS_F);
      if (node != 0)
      {
        ddummy = class_double ( var_ent (node) );
        if (ddummy == 0.01 || ddummy == 0.05 || ddummy == 0.1)
          f = ddummy;
      }
      // eps
      node = btree_FindNode (ent_data (eo), RLAB_NAME_GSL_LS_EPS);
      if (node != 0)
      {
        ddummy = class_double ( var_ent (node) );
        if (ddummy > 0.0)
          eps = ddummy;
      }
    }

    ent_Clean (eo);
  }

  // make sure that the degree of the polynomial is
  // smaller by two than the number of points
  maxdeg = maxdeg > np - 2 ? np - 2 : maxdeg;
  degn   = maxdeg;

  if (f > 0 && !eps) eps = -f;

  //
  // create utility arrays for fit procedure
  //
  r = mdr_Create(np,1);
  a = mdr_Create(3*degn+3*np+3,1);

  if (!w)
  {
    // fit without weights
    w = mdr_Create (np, 1);
    MdrV0(w,0) = -1.0;
    DPOLFT (&np, MDRPTR(x), MDRPTR(y), MDRPTR(w), &maxdeg, &degn, &eps,
             MDRPTR(r), &ierr, MDRPTR(a));
    mdr_Destroy (w);
  }
  else
  {
    // fit with weights provided by the user
    DPOLFT (&np, MDRPTR(x), MDRPTR(y), MDRPTR(w), &maxdeg, &degn, &eps,
             MDRPTR(r), &ierr, MDRPTR(a));
  }

  mdegn = -degn;
  poln = mdr_Create (1, degn + 1);
  DPCOEF (&mdegn, &dzero, MDRPTR(poln), MDRPTR(a));

  mdr_Destroy(a);
  mdr_Destroy(r);

  // clean-up
  ent_Clean (e1);
  ent_Clean (e2);

  //
  // write result to output list
  //
  Btree *bw = btree_Create ();
  // coefficients of the polynomial
  install (bw, "coef", ent_Assign_Rlab_MDR (poln));
  //  epsilon
  install (bw, "error", ent_Create_Rlab_Double(eps));
  // degree of the polynomial
  install (bw, "degree", ent_Create_Rlab_Double(degn));

  return ent_Assign_Rlab_BTREE (bw);
}

#undef  THIS_SOLVER
#define THIS_SOLVER "lsfit"
Ent * ent_gsl_lsfit (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0;
  MDR *y=0, *x=0, *rc=0, *rcov=0, *params0=0;
  double chisq;
  int i;

  //
  // local lsfit parameters and general fit parameters (see odr.c)
  //
  double ls_abserr = 1e-6;
  double ls_relerr = 0.0;
  int ls_convcrit = 0; //0 for delta-test (relerr,abserr), 1 for gradient test (abserr)
  int fit_maxit = 500;

//   FILE *rlab_stderr = (!RLAB_STDERR_DS) ? stderr : RLAB_STDERR_DS;

  //
  // clean-up from previous runs
  //
  xdata = 0;
  ydata = 0;
  wdata = 0;
  xentry = 0;

  //
  // Load and Check arguments.
  //
  if (nargs != 2 && nargs != 3 && nargs != 5 && nargs != 6)
    rerror (THIS_SOLVER ": " RLAB_ERROR_TWO_3_5_OR_6_ARG_REQUIRED);

  //
  // Get observations: y, or <<y;we>>
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) == MATRIX_DENSE_REAL)
  {
    // y
    y = ent_data (e1);
  }
  else if (ent_type (e1) == BTREE)
  {
    // <<val;wgt>>
    ListNode *node;
    node = btree_FindNode (ent_data (e1), RLAB_NAME_GEN_VALUE);
    if (node != 0)
      y = ent_data (var_ent (node));
    if (SIZE(y)<1)
      rerror (THIS_SOLVER ": " RLAB_ERROR_ARG1_STAT);

    node = btree_FindNode (ent_data (e1), RLAB_NAME_GEN_WEIGHT);
    if (node != 0)
    {
      wdata = ent_data (var_ent (node));
      if (SIZE(wdata) != SIZE(y))
        rerror (THIS_SOLVER ": " RLAB_ERROR_ARG1_MISMATCH_VAL_WGT);
      for (i=0; i<SIZE(wdata); i++)
        if (MdrV0 (wdata, i) < 0)
          rerror (THIS_SOLVER ": " RLAB_ERROR_ARG1_MISMATCH_VAL_WGT);
    }
  }
  else
    rerror (THIS_SOLVER ": " RLAB_ERROR_ARG1_STAT);

  if (!EQVECT(y))
    rerror (THIS_SOLVER ": " RLAB_ERROR_ARG_ENTRY_MDR_VECTOR);

  //
  // Get predictor variable: x
  //
  e2 = bltin_get_ent (args[1]);
  if (ent_type (e2) != MATRIX_DENSE_REAL)
    rerror (THIS_SOLVER ": " RLAB_ERROR_ARG2_MDR_MATRIX);
  x = ent_data(e2);

  if (SIZE(x) < 1 || SIZE(y)< 1)
    rerror (THIS_SOLVER ": " RLAB_ERROR_ARG1_ARG2_SAME_NUMBER_ROWS);

  // check x and y
  xr = MNR(x);
  xc = MNC(x);
  if (xr != MNR(y))
    rerror (THIS_SOLVER ": " RLAB_ERROR_ARG1_ARG2_SAME_NUMBER_ROWS);
  if ((xr < xc) || (xr <=1))
    rerror (THIS_SOLVER ": " RLAB_ERROR_LSFIT_TOO_FEW_ENTRIES);

//   fprintf(stderr, "nargs = %i, xc=%i\n", nargs,xc);


  if (nargs == 2)
  {
    //
    // two arguments: use the GSL's simple linear or multi-linear fit
    //
    gsl_matrix *X=0, *COV=0;
    gsl_vector *Y=0, *W=0, *C=0;

    rc   = mdr_Create (1,  xc);
    rcov = mdr_Create (xc, xc);

    X   = gsl_matrix_copy_mdr  (x);
    Y   = alloc_gsl_vector_mdr (y);
    C   = alloc_gsl_vector_mdr (rc);
    COV = alloc_gsl_matrix_mdr (rcov);

    // we either create GSL W, or repack
    // rlab WDATA as W
    if (wdata)
    {
      W = alloc_gsl_vector_mdr (wdata);
// #if (GSL_MAJOR_VERSION <= 1 && GSL_MINOR_VERSION < 16)
    }
    else
    {
      int i;
      W=gsl_vector_alloc (xr);
      for (i = 0; i < xr; i++)
        gsl_vector_set (W, i, 1.0);
    }
// #endif

      if (xc == 1)
      {
        gsl_fit_wmul (MDRPTR(x), 1, W->data, 1, MDRPTR(y), 1, xr, MDRPTR(rc), MDRPTR(rcov), &chisq);
      }
      else
      {
        //
        // multilinear fit: for final result
        //
        gsl_multifit_linear_workspace *work = gsl_multifit_linear_alloc (xr, xc);
        gsl_multifit_wlinear (X, W, Y, C, COV, &chisq, work);
        gsl_multifit_linear_free (work);
      }

// #if ((GSL_MAJOR_VERSION == 1 && GSL_MINOR_VERSION >= 16) || GSL_MAJOR_VERSION>1)
//     }
//     else
//     {
//       gsl_multifit_robust_type const *T = gsl_multifit_robust_bisquare;
//       gsl_multifit_robust_workspace *work = gsl_multifit_robust_alloc (T, xr, xc);
//       gsl_multifit_robust (X, Y, C, COV, work);
//       gsl_multifit_robust_free (work);
//     }
// #else
    if (wdata)
      free_gsl_vector_mdr (W);
    else
      gsl_vector_free (W);
// #endif

    gsl_matrix_free (X);
    free_gsl_vector_mdr (Y);
    free_gsl_vector_mdr (C);
    free_gsl_matrix_mdr (COV);
  }
  else if (nargs == 3 && xc==1)
  {
//     fprintf(stderr, "nargs = %i, xc=%i\n", nargs,xc);

    //
    // lsfit(y,x,p):  p=[1,1,1..] nonzero coefficients represent the powers that are present
    // lsfit(y,x,n):  n           degree of polynomial
    //    polynomial fit with and without template
    //
    Ent *e3=0;
    gsl_matrix *X=0, *COV=0;
    gsl_vector *Y=0, *W=0, *C=0;

    MDR *x3=0, *ipol=0, *rc1=0, *rcov1=0;
    int i, j, zn=0, degn=0;

    e3 = bltin_get_ent (args[2]);
    if (ent_type (e3) != MATRIX_DENSE_REAL)
      rerror (THIS_SOLVER ": " RLAB_ERROR_ARG3_MDR_INTEGER);

    x3 = class_matrix_real (e3);
    if (SIZE(x3) < 1)
      rerror (THIS_SOLVER ": " RLAB_ERROR_ARG3_MDR_INTEGER);

    if (SIZE(x3) == 1)
    {
      degn = Mdr1 (x3, 1, 1);
      zn   = degn + 1;
      ipol = mdi_Create (1, zn);
      for (i = 1; i <= degn + 1; i++)
        MdiV1 (ipol, i) = degn - i + 1;
    }
    else
    {
      degn = SIZE(x3) - 1;
      zn = 0;
      for (i = 1; i <= degn + 1; i++)
      {
        if (MdrV1 (x3, i) != 0)
          zn++;
      }
      if (zn<1)
        rerror (THIS_SOLVER ": " RLAB_ERROR_ARG3_MDR_POSITIVE_SCALAR);
      ipol = mdi_Create (1, zn);
      j = 1;
      for (i = 1; i <= degn + 1; i++)
      {
        if (MdrV1 (x3, i) != 0)
        {
          MdiV1 (ipol, j) = degn - i + 1;
          j ++;
        }
      }
      degn = MdiV1 (ipol, 1);
    }

    //
    // multilinear fit
    //
    rc1   = mdr_Create ( 1, zn);
    rcov1 = mdr_Create (zn, zn);

    Y = alloc_gsl_vector_mdr (y);

    // we either create GSL W, or repack
    // rlab WDATA as W
    if (!wdata)
    {
      W = gsl_vector_alloc (xr);
      for (i = 0; i < xr; i++)
        gsl_vector_set (W, i, 1.0);
    }
    else
      W = alloc_gsl_vector_mdr (wdata);

    // create X for gsl
    X = gsl_matrix_alloc (xr, zn);
    for (i = 0; i < xr; i++)
      for (j = 0; j < zn; j++)
        gsl_matrix_set (X, i, j, pow(mdrV0 (x, i), mdiV0(ipol,j)) );

    C   = alloc_gsl_vector_mdr ( rc1 );
    COV = alloc_gsl_matrix_mdr ( rcov1 );

    gsl_multifit_linear_workspace *work = gsl_multifit_linear_alloc (xr, (ipol->ncol));
    gsl_multifit_wlinear (X, W, Y, C, COV, &chisq, work);
    gsl_multifit_linear_free (work);

    rc   = mdr_Create (1, degn+1);
    rcov = mdr_Create (degn+1, degn+1);
    mdr_Zero (rc);
    mdr_Zero (rcov);
    // the gsl result is copied in reverse order as we are using the
    // matlab notation for the polynomial
    // [a_1,a_2..,a_{N+1}] = a_1 x^N + a_2 x^{N-1} + .... + a_{N+1}
    for (i = 1; i <= zn; i++)
    {
      MdrV1(rc, degn - MdiV1(ipol,i) + 1) = MdrV1(rc1, i);
      for (j = 1; j <= i; j++)
      {
        Mdr1(rcov, degn - MdiV1(ipol,i) + 1, degn - MdiV1(ipol,j) + 1) = Mdr1(rcov1, i, j);
        Mdr1(rcov, degn - MdiV1(ipol,j) + 1, degn - MdiV1(ipol,i) + 1) = Mdr1(rcov1, i, j);
      }
    }

    mdr_Destroy (ipol);
    mdr_Destroy (rc1);
    mdr_Destroy (rcov1);

    free_gsl_vector_mdr (C);
    free_gsl_matrix_mdr (COV);

    if (!wdata)
      gsl_vector_free (W);
    else
      free_gsl_vector_mdr (W);

    free_gsl_vector_mdr (Y);

    gsl_matrix_free (X);

    ent_Clean (e3);

  }
  else if (nargs >= 5)
  {
    Ent *e3=0, *eo=0;
    fname = 0;
    dfname = 0;

    //
    // five or six arguments: nonlinear fit using gsl with an options list
    //
    int i, iter = 0, status, idummy;
    double timer, ddummy;
    time_t t1, t2;
    FILE *fptr = 0;
    ListNode *node;
    char *sout=0;

    // transpose x
    xdata = x;
    if (xc>1)
      mdr_Transpose_inplace (xdata);
    ydata = y;

    //
    // Get parameter array 'p'
    //
    e3 = bltin_get_ent (args[2]);
    if (ent_type (e3) != MATRIX_DENSE_REAL)
      rerror (THIS_SOLVER ": " RLAB_ERROR_ARG3_MDR_VECTOR);
    params0 = ent_data(e3);

    nparams = SIZE(params0);
    if (nparams < 1)
      rerror (THIS_SOLVER ": " RLAB_ERROR_ARG3_MDR_VECTOR);

    //
    // Get function ptrs for 'f' and 'dfdp'
    //
    fname = bltin_get_ent(args[3]);
    if (!isfuncent(fname))
      rerror (THIS_SOLVER ": " RLAB_ERROR_ARG4_FUNC_VAR "\n");

    dfname = bltin_get_ent(args[4]);
    if (!isfuncent(dfname))
      rerror (THIS_SOLVER ": " RLAB_ERROR_ARG5_FUNC_VAR "\n");

    //
    // options for the solver
    //
    if (nargs >= 6)
    {
      eo = bltin_get_ent (args[5]);
      if (ent_type (eo) == BTREE)
      {
        // eabs
        node = btree_FindNode (ent_data (eo), RLAB_NAME_GSL_LS_EABS);
        if (node != 0)
        {
          ddummy = class_double (var_ent (node));
          if (ddummy > 0.0)
            ls_abserr = ddummy;
        }
        // erel
        node = btree_FindNode (ent_data (eo), RLAB_NAME_GSL_LS_EREL);
        if (node != 0)
        {
          ddummy = class_double ( var_ent (node) );
          if (ddummy >= 0.0)
            ls_relerr = ddummy;
        }
        // convergence test
        node = btree_FindNode (ent_data (eo), RLAB_NAME_GSL_LS_CONVT);
        if (node != 0)
        {
          idummy = (int) class_double (var_ent (node));
          if (idummy == 0 || idummy == 1)
            ls_convcrit = idummy;
        }
        // max iterations
        node = btree_FindNode (ent_data (eo), RLAB_NAME_GSL_LS_MAXITER);
        if (node != 0)
        {
          idummy = (int) class_double (var_ent (node));
          if (idummy > 0)
            fit_maxit = idummy;
        }
        // standard output
        node = btree_FindNode (ent_data (eo), RLAB_NAME_GSL_LS_STDOUT);
        if (node != 0)
        {
          if (ent_type(var_ent (node))== MATRIX_DENSE_STRING)
            sout = class_char_pointer (var_ent (node));
          if (isvalidstring(sout)<1)
            sout = 0;
        }
      }
    }

    //
    // Set up ENTITIES for user-function: f = f(xentry,p)
    //
    xentry = mdr_CreateEmpty (1,xc);
    data_ent = ent_Assign_Rlab_MDR (xentry);
    ent_IncRef (data_ent);

    //
    // parameter array p[]: this one changes in solver iterations
    //
    params = mdr_CreateEmpty(1,nparams);
    params_ent = ent_Assign_Rlab_MDR (params);
    ent_IncRef (params_ent);
    //
    // least square function
    //

#if ( (GSL_MINOR_VERSION > 16) || (GSL_MAJOR_VERSION > 1) )
    // changed definition of struct s without advertising.
    gsl_matrix *J;
#endif

#if (GSL_MAJOR_VERSION < 2)
    gsl_multifit_function_fdf f;
    f.fdf = *ls_gslrlab_fdf;
    J = gsl_matrix_alloc(xr, nparams);
#else
    gsl_multifit_nlinear_fdf f;
    f.fvv = NULL;
    f.params = NULL;
#endif
    f.f  = *ls_gslrlab_f;
    f.df = *ls_gslrlab_df;
    f.n  = xr;
    f.p  = nparams;

    //f.params = x->d; //this one creates seg fault. who would know why?
    gsl_matrix *cov=0;
    gsl_vector *g=0;

    rcov = mdr_Create (nparams, nparams);
    rc   = mdr_Create (1, nparams);
    cov  = alloc_gsl_matrix_mdr (rcov); // this is symmetric, so we do not need to transpose

    //
    // prepare iterator
    //
    gsl_set_error_handler_off ();

#if (GSL_MAJOR_VERSION < 2)

    const gsl_multifit_fdfsolver_type *T = gsl_multifit_fdfsolver_lmsder;
    gsl_multifit_fdfsolver *s = gsl_multifit_fdfsolver_alloc (T, xr, nparams);
    gsl_vector_view xvec = gsl_vector_view_array (MDRPTR(params0), nparams);
    gsl_multifit_fdfsolver_set (s, &f, &xvec.vector);
    if (ls_convcrit == 1)
      g = gsl_vector_alloc (xc);

#else

    int info;
    const gsl_multifit_nlinear_type *T = gsl_multifit_nlinear_trust;
    gsl_multifit_nlinear_workspace *w;
    gsl_multifit_nlinear_parameters fdf_params = gsl_multifit_nlinear_default_parameters();
    gsl_vector_view xvec = gsl_vector_view_array (MDRPTR(params0), nparams);

    /* allocate workspace with default parameters */
    w = gsl_multifit_nlinear_alloc (T, &fdf_params, xr, nparams);
    /* initialize solver with starting point */
    gsl_multifit_nlinear_init (&xvec.vector, &f, w);

#endif

    if (ls_convcrit == 1)
      g = gsl_vector_alloc (xc);

    //
    // iterate
    //
    if (sout)
    {
      fptr = fopen (sout, "a");
      t1 = clock ();
      if (fptr)
        fprintf (fptr,
                 "RLaB: using GSL non-linear least-square solver %s.\n",
#if (GSL_MAJOR_VERSION < 2)
                 gsl_multifit_fdfsolver_name (s));
#else
                 gsl_multifit_nlinear_name (w));
#endif
    }

    //
    //
    //
    do
    {
      iter++;
#if (GSL_MAJOR_VERSION < 2)
      status = gsl_multifit_fdfsolver_iterate (s);
      if (ls_convcrit == 0)
      {
        status = gsl_multifit_test_delta (s->dx, s->x, ls_abserr, ls_relerr);
      }
      else
      {
  #if ((GSL_MAJOR_VERSION == 1) && (GSL_MINOR_VERSION <= 16))
        gsl_multifit_gradient (s->J, s->f, g);
  #else
        gsl_multifit_gradient (J, s->f, g);
  #endif
        status = gsl_multifit_test_gradient (g, ls_abserr);
      }
      if (fptr)
      {
        fprintf (fptr, "iter = %3i: p = ", iter);
        for (i = 0; i < nparams; i++)
          fprintf (fptr, "  %g", gsl_vector_get (s->x, i));
        fprintf (fptr, "  ||df|| = %g\n", gsl_blas_dnrm2 (s->f));
      }

#else

      status = gsl_multifit_nlinear_iterate (w);
      status = gsl_multifit_nlinear_test(ls_abserr,  pow(GSL_DBL_EPSILON,0.333), ls_relerr,
                                          &info, w);
      if (fptr)
      {
        fprintf (fptr, "iter = %3i: p = ", iter);
        /*  */
        gsl_vector * residual_f = gsl_multifit_nlinear_residual(w);
        gsl_blas_ddot(residual_f, residual_f, &chisq);

        for (i = 0; i < nparams; i++)
          fprintf (fptr, "  %g", gsl_vector_get (w->x, i));
        fprintf (fptr, "  ||df|| = %g\n", sqrt(chisq));
      }
#endif
    }
    while (status == GSL_CONTINUE && iter < fit_maxit);

    if (fptr)
    {
      fprintf (fptr, "RLaB gsl non-linear least-square solver");
      fprintf (fptr, " reports:");
      fprintf (fptr, " '%s' !\n", gsl_strerror (status));
      // check the time and close the output stream
      t2 = clock ();
      timer = (t2 - t1) / 1e6;
      fprintf (fptr, "RLaB: least-square fitting lasted %g sec.\n", timer);
      fclose (fptr);
    }
    //
    // create a list with entries 'coef', 'cov', copy to rlab variables
    //
    if (status == GSL_SUCCESS)
    {

#if (GSL_MAJOR_VERSION < 2)
      for (i = 0; i < nparams; i++)
        MdrV0 (rc, i) = gsl_vector_get (s->x, i);
      chisq = pow (gsl_blas_dnrm2 (s->f), 2);
  #if ((GSL_MAJOR_VERSION == 1) && (GSL_MINOR_VERSION <= 16))
      gsl_multifit_covar (s->J, 0, cov);
  #else
      gsl_multifit_covar (J, 0, cov);
  #endif
#else
      J = gsl_multifit_nlinear_jac(w);
      gsl_multifit_nlinear_covar (J, 0.0, cov);

      for (i = 0; i < nparams; i++)
        MdrV0 (rc, i) = gsl_vector_get (w->x, i);
#endif
    }
    else
    {
      mdr_Destroy (rc);
      mdr_Destroy (rcov);
      rc = mdr_Create (0, 0);
      rcov = mdr_Create (0, 0);
      chisq = 0;
    }

    // put x back
    if (xc>1)
      mdr_Transpose_inplace (xdata);

    // Clean Up
    MDPTR(params)=0;
    ent_DecRef (params_ent);
    ent_Destroy (params_ent);

    MDPTR(xentry)=0;
    ent_DecRef (data_ent);
    ent_Destroy (data_ent);

    //gsl_matrix_free (cov);
    if (ls_convcrit == 1)
      gsl_vector_free (g);

#if (GSL_MAJOR_VERSION < 2)

    gsl_multifit_fdfsolver_free (s);

#else

    gsl_multifit_nlinear_free (w);

#endif

    free_gsl_matrix_mdr(cov);

    ent_Clean (e3);
    ent_Clean (eo);
    ent_Clean (fname);
    ent_Clean (dfname);

#if ( (GSL_MINOR_VERSION > 16) && (GSL_MAJOR_VERSION == 1) )
    gsl_matrix_free(J);
#endif

  }
  else
  {
    //
    rerror (THIS_SOLVER ": missing arguments !");
  }

  // clean-up
  ent_Clean (e1);
  ent_Clean (e2);
  //
  // write result to output list
  //
  Btree *bw = btree_Create ();
  install (bw, "coef", ent_Assign_Rlab_MDR (rc));
  install (bw, "cov", ent_Assign_Rlab_MDR (rcov));
  install (bw, "chisq", ent_Create_Rlab_Double (chisq));
  return ent_Assign_Rlab_BTREE (bw);
}

//
// rlab least-square function. has to be optimized for matrix argument x
//
static int
ls_gslrlab_f (const gsl_vector * x, void *dummy, gsl_vector * f)
{
  //double rks_func (const gsl_vector *x, void *params){
  int i;
  double rval;

//   for (i = 0; i < nparams; i++)
//     MdrV0 (params, i) = gsl_vector_get (x, i);
  MDPTR (params) = (void *) x->data;

  for (i = 0; i < xr; i++)
  {
    Ent *rent = 0;

    // prepare 'xentry'
    if (xc>1)
      MDPTR (xentry) = (void *) &Mdr0 (xdata, 0, i);
    else
      MDPTR (xentry) = (void *) &MdrV0 (xdata, i);

    rent = ent_call_rlab_script_2args(fname, data_ent, params_ent);

    if (ent_type(rent)!=MATRIX_DENSE_REAL)
      rerror (THIS_SOLVER ": " RLAB_ERROR_RHS_FUNC_MUST_RETURN_MDR);

    rval = class_double (rent);

    if (wdata)
    {
      gsl_vector_set (f, i,
                      (rval - MdrV0 (ydata, i)) * sqrt (MdrV0 (wdata, i)));
    }
    else
      gsl_vector_set (f, i, rval - MdrV0 (ydata, i));

    ent_Clean (rent);
  }

  return GSL_SUCCESS;
}

//
// jacobian for the rlab least-square function.
//
static int
ls_gslrlab_df (const gsl_vector * x, void *dummy, gsl_matrix * J)
{
  //double rks_func (const gsl_vector *x, void *params){
  int i,j;
  MDR *retm;

//   for (i = 0; i < nparams; i++)
//     MdrV0 (params, i) = gsl_vector_get (x, i);

  MDPTR (params) = (void *) x->data;

  for (i = 0; i < xr; i++)
  {
    Ent *rent = 0;

    // prepare 'xentry'
    if (xc>1)
      MDPTR (xentry) = (void *) &Mdr0 (xdata, 0, i);
    else
      MDPTR (xentry) = (void *) &MdrV0 (xdata, i);
//     for (j=0; j<xc; j++)
//       MdrV0 (xentry, j) = Mdr0 (xdata, i, j);

    rent = ent_call_rlab_script_2args(dfname, data_ent, params_ent);

    if (ent_type(rent)!=MATRIX_DENSE_REAL)
      rerror (THIS_SOLVER ": " RLAB_ERROR_RHS_FUNC_MUST_RETURN_MDR);

    retm = ent_data (rent);

    if (SIZE(retm) != nparams)
      rerror (THIS_SOLVER ": " RLAB_ERROR_RHS_FUNC_INCORRECT_DIM);


    //
    // initialize Jacobian data pointer (J has already been declared by the GSL)
    //
    for (j = 0; j < nparams; j++)
    {
      if (wdata)
        gsl_matrix_set (J, i, j,
                        MdrV0 (retm, j) * sqrt (MdrV0 (wdata, i)));
      else
        gsl_matrix_set (J, i, j, MdrV0 (retm, j));
    }

    ent_Clean (rent);
  }

  return GSL_SUCCESS;
}

