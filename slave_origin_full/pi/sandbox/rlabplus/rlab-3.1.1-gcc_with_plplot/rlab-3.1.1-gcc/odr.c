// Copyright (C) 2003-2008 Marijan Kostrun
//   part of rlab+4linux project on rlabplus.sourceforge.net
//
// RLaB2 Rel.2 driver for ODRPACK
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

#include "odrpack.h"
#include "rfileio.h"

//
// wrappers for rlab functions:
//
int odrfunc (double *, double *, double *);
int odrdxfn (double *, double *, double *);
int odrdbfn (double *, double *, double *);


// naming convention for the solver parameters
#include "rlab_solver_parameters_names.h"

// their rlab names and parameter lists
static Ent *fname, *dxfname, *dbfname;

// communication between rlab functions and the odrpack
static Ent *b1_ent = 0, *x1_ent = 0;
static MDR *b1 = 0, *x1 = 0;

static int odr_np = 0, odr_nq = 0, odr_n = 0, odr_m = 0;

static int odr_fjacb = 0;        // do not have jacobian dfdb - find it numerically
static int odr_fjacd = 0;        // do not have jacobian dfdx - find it numerically
static int odr_fjac = 0;
static int odr_ndiff = 0;        // relevant if code finds jacobians numerically:
                                 //  0 - use forward difference formula

#undef  THIS_SOLVER
#define THIS_SOLVER "odrfit"
Ent *
ent_odr_basic (int nargs, Datum args[])
{
  // call parameters:
  //  d1 - f(x,b),
  //  d2 - dfdx(x,b), x-coordinate jacobian
  //  d3 - dfdb(x,b), b-parameter jacobian

  //
  // odr parameters
  //
  int odr_method = 0;   // 0 - explicit orthogonal distance regression,
                        // 1 - implicit orthogonal distance regression,
                        // 2 - least square fit.

  odr_np = 0;
  odr_nq = 0;
  odr_n  = 0;
  odr_m  = 0;
  odr_fjacb = 0;
  odr_fjacd = 0;
  odr_fjac  = 0;
  odr_ndiff = 0;

  char *fit_sout=0;
  double odr_taufac = 1.0;  // size of the trust region for Gauss-Newton step
  double odr_sstol = -1.0;  // stopping tolerance for sum of squares (deltas)
  double odr_partol = -1.0; // stopping tolerance for parameter convergence
  int fit_maxit = 100;      // maximum number of iterations

  Ent *e1=0, *e2=0, *e3=0;
  MDR *x=0, *b=0, *beta=0, *work=0, *d=0, *e=0, *cov=0;
  int *ifixb=0, *ifixx=0, *iwork=0, liwork, info, lwork, ldstpd, idummy;
  int ndigit, job, ldifx, ld2wd, ldwd;
  int ld2we, ldwe, ldscld, lodrout=0;

  double stpb, stpd, scld, ddummy;

  int ny, i, j;

  double minus_done=-1;
  double *sclb=0;

  Btree *bw=0;

  FILE *fptr = NULL;
  time_t timer;

  FILE *rlab_stderr = (!RLAB_STDERR_DS) ? stderr : RLAB_STDERR_DS;

  MDR *difixx  = 0;
  MDR *difixb  = 0;
  MDR *dscaleb = 0;
  MDR *wd      = 0;
  MDR *we      = 0;
  MDR *y       = 0;

  int jobI4 = 0;

  ListNode *node;
  Ent *eo=0;

  fname=0;
  dxfname=0;
  dbfname=0;

  //
  // Load and Check arguments.
  //
  if (nargs < 4)
  {
    fprintf(rlab_stderr, THIS_SOLVER ": Basic orthogonal distance regression / least-squares solver (ODRPACK).\n");
    fprintf(rlab_stderr, THIS_SOLVER ": Format:\n");
    fprintf(rlab_stderr, THIS_SOLVER ":   c=odrfit(/Y/,X,p,f/,dfdp,dfdx,options/),\n");
    fprintf(rlab_stderr, THIS_SOLVER ": where 'Y' is a dependent variable, either as a row-wise matrix of values\n");
    fprintf(rlab_stderr, THIS_SOLVER ": or a structured list Y=<<val;wgt>> of values and weights \n");
    fprintf(rlab_stderr, THIS_SOLVER ": (omitting 'Y' implies an implicite model);\n");
    fprintf(rlab_stderr, THIS_SOLVER ": 'X' is an independent row-wise variable or a list\n");
    fprintf(rlab_stderr, THIS_SOLVER ": X=<<val;wgt;fix_x;dx>> of values, weights, fixed components\n");
    fprintf(rlab_stderr, THIS_SOLVER ": and the initial guess for errors 'dx'; 'p' is an initial\n");
    fprintf(rlab_stderr, THIS_SOLVER ": guess for parameter array; f=function(x,b) a function describing the\n");
    fprintf(rlab_stderr, THIS_SOLVER ": model; 'dfdp' a jacobian of 'f' with respect to  'p',\n");
    fprintf(rlab_stderr, THIS_SOLVER ": 'dfdx' a jacobian of 'f' with respect to 'x', where the jacobians\n");
    fprintf(rlab_stderr, THIS_SOLVER ": may be omitted (in that case, the code computes them internally);\n");
    fprintf(rlab_stderr, THIS_SOLVER ": options=<<stdout;fix_p;imethod;taufac;sstol;partol;maxi>> is a parameter\n");
    fprintf(rlab_stderr, THIS_SOLVER ": list which controls computation.\n");
    fprintf(rlab_stderr, THIS_SOLVER ": The function returns  c=<<coef;cov;dy;dx>>, where 'coef' are the new\n");
    fprintf(rlab_stderr, THIS_SOLVER ": values for 'p', 'cov' is their covariance matrix, 'dx' (if using\n");
    fprintf(rlab_stderr, THIS_SOLVER ": odr and not ls) is the error in 'x', while 'dy' is the error in\n");
    fprintf(rlab_stderr, THIS_SOLVER ": function evaluation 'f'.\n");
    fprintf(rlab_stderr, THIS_SOLVER ": Check the manuals for more details.\n");
    rerror (THIS_SOLVER ": " RLAB_ERROR_AT_LEAST_FOUR_ARG_REQUIRED);
  }
  //
  // y
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) == MATRIX_DENSE_REAL)
    y = class_matrix_real (e1);
  else if (ent_type (e1) == BTREE)
  {
    // must have x
    node = btree_FindNode (ent_data (e1), "val");
    if (node != 0)
      y = class_matrix_real (var_ent (node));
    // is there wgt ?
    node = btree_FindNode (ent_data (e1), "wgt");
    if (node != 0)
      we = mdr_Float_BF(class_matrix_real (var_ent (node)));
  }

  //
  // x
  //
  e2 = bltin_get_ent (args[1]);
  if (ent_type (e2) == MATRIX_DENSE_REAL)
    x = class_matrix_real (e2);
  else if (ent_type (e2) == BTREE)
  {
    // a list with at most three entries: x, wd, ifixx
    //
    // must have x
    node = btree_FindNode (ent_data (e2), "val");
    if (node == 0)
      rerror
          ("odrfit: improper structured list for independent variable: missing 'x'.");
    x = class_matrix_real (var_ent (node));
    // is there ifixx?
    node = btree_FindNode (ent_data (e2), "fix_x");
    if (node != 0)
      difixx = class_matrix_real (var_ent (node));
    // is there wd?
    node = btree_FindNode (ent_data (e2), "wgt");
    if (node != 0)
      wd = mdr_Float_BF( class_matrix_real (var_ent (node)) );
    // is there Delta(0)?
    node = btree_FindNode (ent_data (e2), "dx");
    if (node != 0)
    {
      d = mdr_Float_BF ( class_matrix_real (var_ent (node)) );
      if (MNR (d) == MNR (x) && MNC (d) == MNC (x))
        jobI4 = 1;
      else
        mdr_Destroy (d);
    }
  }
  else
    rerror ("odrfit: 'x' is a real matrix or a list <<val;wgt;fix_x;dx>>!");

  //
  // process x
  //
  if (SIZE (x) < 1)
    rerror ("odrfit: 'x' has to be a real matrix of non-zero size!");
  odr_n = MNR (x);
  odr_m = MNC (x);
  if (jobI4 == 0)
    d = mdr_Float_BF (x);

  //
  // wd
  //
  if (!wd)
  {
    ldwd = 1;
    ld2wd = 1;
    wd = mdr_CreateScalar (-1.0);
  }
  else
  {
    if ( MNR(wd)*MNC(wd) == 1)
    {
      // wd(1,1,1) relative weight compared to we
      ld2we = 1;
      ldwe = 1;
      if (MdrV0(wd,0) > 0)
        MdrV0(wd,0) *= -1.0;
      else
        MdrV0(wd,0)  = -1.0;
    }
    else if (MNR (wd) == 1 && MNC (wd) == odr_m)
    {
      // wd(1,1, 1:m )
      ldwd = 1;
      ld2wd = 1;
    }
    else if (MNR (wd) == odr_m && MNC (wd) == odr_m)
    {
      // wd(1, 1:m , 1:m )
      ldwd = 1;
      ld2wd = odr_m;
    }
    else if (MNR (wd) == odr_n && MNC (wd) == odr_m)
    {
      // wd( 1:n , 1 , 1:m )
      ldwd = odr_n;
      ld2wd = 1;
    }
    else if ( MNR(wd) * MNC (wd) == odr_n * odr_m * odr_m )
    {
      // wd( 1:n , 1:m , 1:m )
      ldwd  = odr_n;
      ld2wd = odr_m;
    }
    else
      rerror ("odrfit: improper dy!");
  }

  //
  // ifixx
  //
  if (difixx)
  {
    // ifixx as a row vector or a full matrix
    ldifx = MNR (difixx);
    if (ldifx != 1 && ldifx != odr_n)
      rerror ("odrfit: improper 'fix_x'!");
    if (MNC (difixx) != 1 && MNC (difixx) != odr_m)
      rerror ("odrfit: improper 'fix_x'!");
    ifixx = (int *) GC_malloc (ldifx * odr_m * sizeof (int));
    for (i = 0; i < ldifx; i++)
      for (j = 0; j < odr_m; j++)
    {
      if (difixx->type == RLAB_TYPE_INT32)
        ifixx[i + ldifx * j] = ! Mdi0 (difixx, MIN(0, i), MIN(0, j));
      else
        ifixx[i + ldifx * j] = !((int) Mdr0 (difixx, MIN(0, i), MIN(0, j)));
    }
  }
  else
  {
    // no ifixx
    ldifx = 1;
    ifixx = (int *) GC_malloc (ldifx * sizeof (int));
    ifixx[0] = -1;
  }


  //
  // process y
  //
  if (!y)
  {
    // implicite fit
    odr_method = 1;
    odr_nq = 1;
    ny     = 0;
    // fake y
    y = x;
  }
  else
  {
    // explicite fit
    odr_nq = MNC (y);
    ny     = MNR (y);
  }

  //
  // process we
  //
  if (!we)
  {
    we = mdr_CreateScalar (-1.0);
    ld2we = 1;
    ldwe  = 1;
  }
  else
  {
    if ( MNR(we)*MNC(we) == 1)
    {
      ld2we = 1;
      ldwe = 1;
      if (MdrV0(we,0) > 0)
        MdrV0(we,0) *= -1.0;
      else
        MdrV0(we,0)  = -1.0;
    }
    else if (MNR (we) == odr_nq && MNC (we) == odr_nq)
    {
      // we(1, 1:q , 1:q )
      ldwe = 1;
      ld2we = odr_nq;
    }
    else if (MNR (we) == odr_n && MNC (we) == odr_nq)
    {
      // we( 1:n , 1 , 1:q )
      ldwe = odr_n;
      ld2we = 1;
    }
    else if (MNR (we) * MNC (we) == odr_n * odr_nq * odr_nq)
    {
      // we( 1:n , 1:q , 1:q )
      ldwe = odr_n;
      ld2we = odr_nq;
    }
    else
      rerror ("odrfit: improper 'dy' entry in the list!");
  }

  //
  // Get parameter array
  //
  e3 = bltin_get_ent (args[2]);
  if (ent_type (e3) != MATRIX_DENSE_REAL)
    rerror ("odrfit: 'p' has to be a real matrix!");
  b = class_matrix_real (e3);
  odr_np = SIZE (b);
  if (odr_np<1)
    rerror (THIS_SOLVER ":" RLAB_ERROR_ARG3_MDR_VECTOR);
  if (!EQVECT(b))
    rerror (THIS_SOLVER ":" RLAB_ERROR_ARG3_MDR_VECTOR);
  beta = mdr_Float_BF (b);

  //
  // Get f = f(x,b)
  //
  fname = bltin_get_ent(args[3]);
  if (!isfuncent(fname))
    rerror (THIS_SOLVER ": " RLAB_ERROR_ARG4_FUNC_VAR "\n");

  //
  // Get dfdb = dfdb(x,b) if given, otherwise find it numerically
  //
  eo = 0;
  odr_fjacb = 0;
  if(nargs > 4)
  {
    dbfname = bltin_get_ent(args[4]);
    if (isfuncent(dbfname))
      odr_fjacb = 1;
    else
    {
      if (ent_type (dbfname) == BTREE)
        eo = dbfname;
      odr_fjacb = 0;
      dbfname = 0;
    }
  }

//   fprintf(stderr, "odr_fjacb = %i\n", odr_fjacb);

  //
  // Get dfdx = dfdx(x,b) if given, otherwise find it numerically
  //
  odr_fjacd = 0;
  if(nargs > 5)
  {
    dxfname = bltin_get_ent(args[5]);
    if (isfuncent(dxfname))
      odr_fjacd = 1;
    else
    {
      if (ent_type (dxfname) == BTREE)
      {
        eo = dxfname;
      }
      dxfname = 0;
      odr_fjacd = 0;
    }
  }

  if (nargs > 6)
    eo = bltin_get_ent (args[6]);

  if (eo!=0)
    if (ent_type (eo) == BTREE)
  {
    // ifixb
    node = btree_FindNode (ent_data (eo), RLAB_NAME_ODR_FIXP);
    if (node != 0)
    {
      if (ent_type(var_ent (node)) == MATRIX_DENSE_REAL)
      {
        difixb = class_matrix_real (var_ent (node));
      }
    }
    // scale parameter array
    node = btree_FindNode (ent_data (eo), RLAB_NAME_ODR_SCALEP);
    if (node != 0)
    {
      if (ent_type(var_ent (node)) == MATRIX_DENSE_REAL)
        dscaleb = class_matrix_real (var_ent (node));
    }
    // taufac
    node = btree_FindNode (ent_data (eo), RLAB_NAME_ODR_TAUFAC);
    if (node != 0)
    {
      if (ent_type(var_ent (node)) == MATRIX_DENSE_REAL)
      {
        ddummy = class_double (var_ent (node));
        if (ddummy > 0.0)
          odr_taufac = ddummy;
      }
    }
    // sstol
    node = btree_FindNode (ent_data (eo), RLAB_NAME_ODR_SSTOL);
    if (node != 0)
    {
      if (ent_type(var_ent (node)) == MATRIX_DENSE_REAL)
      {
        ddummy = class_double ( var_ent (node) );
        if (ddummy >= 0.0)
          odr_sstol = ddummy;
      }
    }
    // partol
    node = btree_FindNode (ent_data (eo), RLAB_NAME_ODR_PARTOL);
    if (node != 0)
    {
      if (ent_type(var_ent (node)) == MATRIX_DENSE_REAL)
      {
        ddummy = class_double ( var_ent (node) );
        if (ddummy >= 0.0)
          odr_partol = ddummy;
      }
    }
    // method
    node = btree_FindNode (ent_data (eo), RLAB_NAME_ODR_IMETHOD);
    if (node != 0)
    {
      if (ent_type(var_ent (node)) == MATRIX_DENSE_REAL)
      {
        idummy = (int) class_double (var_ent (node));
        if (idummy == 0 || idummy == 1)
          odr_method = 2*idummy;
      }
    }
    // max iterations
    node = btree_FindNode (ent_data (eo), RLAB_NAME_ODR_MAXITER);
    if (node != 0)
    {
      if (ent_type(var_ent (node)) == MATRIX_DENSE_REAL)
      {
        idummy = (int) class_double (var_ent (node));
        if (idummy > 0)
          fit_maxit = idummy;
      }
    }
    // standard output
    node = btree_FindNode (ent_data (eo), RLAB_NAME_ODR_STDOUT);
    if (node != 0)
    {
      if (ent_type(var_ent (node))== MATRIX_DENSE_STRING)
      {
        fit_sout = class_char_pointer (var_ent (node));
        lodrout = isvalidstring (fit_sout);
        if (lodrout < 1)
        {
          fit_sout = 0;
          lodrout = 0;
        }
      }
    }
  }

  if (odr_method == 0 || odr_method == 1)
    odr_fjac = odr_fjacb && odr_fjacd;	// for odr need dfdp and dfdx
  else
    odr_fjac = odr_fjacb;	// for ls need only dfdp

  //
  // Set up static ENTITIES for user-functions
  // x1, a single data point
  x1 = mdr_Create (1, odr_m);
  x1_ent = ent_Assign_Rlab_MDR (x1);
  ent_IncRef (x1_ent);

  // b1, a single value of parameter array
  b1 = mdr_CreateEmpty (1, odr_np);
  b1_ent = ent_Assign_Rlab_MDR (b1);
  ent_IncRef (b1_ent);

  //
  // iwork
  //
  liwork = 20 + odr_np + odr_nq * (odr_np + odr_m);
  iwork = (int *) GC_malloc (liwork * sizeof (int));

  //
  // ifixb: do not fix it
  //
  if (SIZE(difixb) == odr_np)
  {
    ifixb = (int *) GC_malloc (odr_np * sizeof (int));
    for (i = 0; i < odr_np; i++)
    {
      ifixb[i] = mdiV0(difixb,i) == 1 ? 0 : 1;
    }
  }
  else
  {
    ifixb  = (int *) GC_malloc (sizeof (int));
   *ifixb = -1;
  }

  // what to do: job = (I5)(I4)(I3)(I2)(I1)
  // -> method in 'odr_method'                (I1=odr_method)
  // -> always calculate Vbeta and sigma_beta (I3=0)
  // -> set  Delta = 0                        (jobI4=0) or is given (jobI4=1)
  // -> not a restart                         (I5=0)
  job = jobI4 * 1000;
  if (odr_fjac)
    job += odr_method + 20;	// user supplied derivatives checked by odrpack (I2=2)
  else
    job += odr_method + odr_ndiff * 10;	// calculate jacobians numerically              (I2=0,1)


  // figure from the problem the precision of calculation
  ndigit = -1;

  //
  // do not scale: user can do it when she prepares the model and the data
  //
  stpb = -1.0;
  stpd = -1.0;
  ldstpd = 1;

  if (dscaleb)
    sclb = MDRPTR(dscaleb);
  else
    sclb = &minus_done;

  scld = -1.0;
  ldscld = 1;

  //
  // work
  //
  if (odr_method == 0 || odr_method == 1)
  {
    // orthogonal distance
    lwork = 18 + 11 * odr_np + odr_np * odr_np + odr_m + odr_m * odr_m
        + 4 * odr_n * odr_nq + 6 * odr_n * odr_m + 2 * odr_n * odr_nq * odr_np
        + 2 * odr_n * odr_nq * odr_m + odr_nq * odr_nq
        + 5 * odr_nq + odr_nq * (odr_np + odr_m) + ldwe * ld2we * odr_nq;
  }
  else
  {
    // ordinary least square
    lwork = 18 + 11 * odr_np + odr_np * odr_np
        + odr_m + odr_m * odr_m + 4 * odr_n * odr_nq
        + 2 * odr_n * odr_m + 2 * odr_n * odr_nq * odr_np
        + 5 * odr_nq + odr_nq * (odr_np + odr_m) + ldwe * ld2we * odr_nq;
  }
  work = mdr_Create (lwork, 1);
  if (jobI4)
    for (i = 0; i < odr_n; i++)
      for (j = 0; j < odr_m; j++)
        MdrV0(work,i+j*odr_n) = Mdr0 (d, i, j);

  //
  // run time messages: open file and then close it. fortran has its own fopen
  //
  if (fit_sout)
  {
    fptr = fopen (fit_sout, "a");
    if (fptr)
    {
      // start the timer
      timer = clock ();
      // write down boring info
      fprintf(fptr, THIS_SOLVER ": using built-in Orthogonal Distance Regression solver (ODRPACK).\n");
      fprintf(fptr, THIS_SOLVER ": the messages from the solver follow.\n");
    }
    else
    {
      fit_sout = 0;
      lodrout = 0;
    }
  }

  //
  // Finally, call the rlab-wrapper for the solver
  //
  DODRCRL (ODRFDF, &odr_n, &odr_m, &odr_np, &odr_nq,
            MDRPTR(beta),
            MDRPTR(y), &ny, MDRPTR(x), &odr_n,
            MDRPTR(we), &ldwe, &ld2we, MDRPTR(wd), &ldwd, &ld2wd,
            ifixb, ifixx, &ldifx,
            &job, &ndigit, &odr_taufac,
            &odr_sstol, &odr_partol, &fit_maxit,
            fit_sout, &lodrout,
            &stpb, &stpd, &ldstpd,
            sclb, &scld, &ldscld, MDRPTR(work), &lwork, iwork, &liwork, &info);


  if (info <= 4)
  {
    // converged
    if (fptr)
    {
      if (info == 4)
        fprintf(fptr, THIS_SOLVER ": built-in solver ODRPACK reports iteration limit reached!\n");
      else
        fprintf(fptr, THIS_SOLVER ": built-in solver ODRPACK reports success.\n");
      // check the time and close the output stream
      timer -= clock ();
      fprintf(fptr, THIS_SOLVER ": ODRPACK calculation lasted %g sec.\n",
               -timer / 1e6);
      fflush (fptr);
      fclose (fptr);
    }
  }
  else
  {
    // did not converge
    if (fptr)
    {
      fprintf(fptr, THIS_SOLVER ": built-in solver ODRPACK reports error %i\n", info);
      // check the time and close the output stream
      timer -= clock ();
      fprintf(fptr, THIS_SOLVER ": ODRPACK calculation lasted %g sec.\n",
               -timer / 1e6);
      fclose (fptr);
    }
  }

  // Clean Up
  GC_FREE (iwork);
  GC_FREE (ifixx);
  GC_FREE (ifixb);
  mdr_Destroy(we);
  mdr_Destroy(wd);

  ent_DecRef (x1_ent);
  ent_Destroy (x1_ent);

  MDPTR (b1)=0;
  ent_DecRef (b1_ent);
  ent_Destroy (b1_ent);

  ent_Clean (e1);
  ent_Clean (e2);
  ent_Clean (e3);

  ent_Clean (fname);
  ent_Clean (dxfname);
  ent_Clean (dbfname);

  ent_Clean (eo);

  //
  // return the result as a list
  //
  bw = btree_Create ();

  // status matrix
  install (bw, "status", ent_Create_Rlab_Double((double) info));

  // coef matrix
  install (bw, "coef", ent_Assign_Rlab_MDR(beta));

  // dx for odr, not for ls
  if (odr_method == 0 || odr_method == 1)
  {
    for (i = 0; i < odr_n; i++)
      for (j = 0; j < odr_m; j++)
        Mdr0 (d, i, j) = MdrV0(work,i+j*odr_n);
    install (bw, "dx", ent_Assign_Rlab_MDR(d));
  }

  // dy
  e = mdr_Create (odr_n, odr_nq);
  for (i = 0; i < odr_n; i++)
    for (j = 0; j < odr_nq; j++)
      Mdr0 (e, i, j) = MdrV0(work,odr_n*odr_m+i+j*odr_n);
  install (bw, "dy", ent_Assign_Rlab_MDR(e));

  // cov
  cov = mdr_Create (odr_np, odr_np);
  for (i = 0; i < odr_np; i++)
    for (j = 0; j < odr_np; j++)
  {
      // entries scaled with residual covariance
    Mdr0 (cov, i, j) =
        MdrV0(work,2*odr_n*odr_m+2*odr_n*odr_nq+odr_np+i+j*odr_np);
    Mdr0 (cov, i, j) *=
        MdrV0(work,2*odr_n*odr_m+2*odr_n*odr_nq+odr_np+odr_np*odr_np);
  }
  install (bw, "cov", ent_Assign_Rlab_MDR(cov));

  // finish cleanup
  mdr_Destroy (work);

  return ent_Assign_Rlab_BTREE(bw);
}


int
ODRFDF (int * N, int * M, int * NP, int * NQ, int * LDN, int * LDM,
        int * LDNP, double * BETA, double * XPLUSD, int * IFIXB,
        int * IFIXX, int * LDIFX, int * IDEVAL, double * F,
        double * FJACB, double * FJACD, int * ISTOP)
{
  if ( *ISTOP != 0 ) return 0;

  if ( (*IDEVAL)%10 >= 1 )
  {
    odrfunc(XPLUSD, BETA, F);
    *ISTOP = 0;
  }

  if ( ((*IDEVAL)/10) %10 >= 1 )
    odrdbfn( XPLUSD, BETA, FJACB);

  if ( ((*IDEVAL)/100)%10 >= 1 )
    odrdxfn(XPLUSD, BETA, FJACD);

  return 0;
}

//
// fortran: f = f(x,b)
//
int
odrfunc (double *x, double *b, double *f)
{
  int i, j;
  //
  // b is the same for each datum x
  //
  MDPTR (b1) = (void *) b;

  for (i=0; i<odr_n; i++)
  {
    Ent *rent = 0;
    MDR *retm = 0;

    //
    // x1 = x[i;]
    //
    for (j=0; j<odr_m; j++)
      MdrV0 (x1, j) = x[i + odr_n * j];

    //
    // do my bidding
    //
    rent = ent_call_rlab_script_2args(fname, x1_ent, b1_ent);

    if (ent_type(rent)!=MATRIX_DENSE_REAL)
      rerror (THIS_SOLVER ": " RLAB_ERROR_RHS_FUNC_MUST_RETURN_MDR);

    retm = ent_data(rent);

    if (SIZE(retm) != odr_nq)
    {
//       fprintf(stderr, THIS_SOLVER ": Call to function : %i != %i\n", SIZE(retm), odr_nq);
      rerror (THIS_SOLVER ": " RLAB_ERROR_RHS_FUNC_INCORRECT_DIM);
    }

    //
    // f -> f(x,b)
    //
    for (j=0; j<odr_nq; j++)
      f[i + odr_n * j] = mdrV0 (retm, j);

    ent_Clean (rent);
  }

  return 1;
}


//
// fortran: dfdb = dfdb(x,b)
//
int
odrdbfn (double *x, double *b, double *fjacb)
{
  int i, j, k;

  // do not have jacobian dfdb - find it numerically
  if (odr_fjacb == 0)
    return 1;

  //
  // b is the same for each datum x
  //
  MDPTR (b1) = (void *) b;

  //
  // go over each x, evaluate dfdx and store result
  //
  for (i=0; i<odr_n; i++)
  {
    Ent *rent = 0;
    MDR *retm = 0;

    //
    // x1 = x[i;]
    //
    for (j=0; j<odr_m; j++)
      MdrV0 (x1, j) = x[i + odr_n * j];

    //
    // do my bidding
    //
    rent = ent_call_rlab_script_2args(dbfname, x1_ent, b1_ent);

    if (ent_type(rent)!=MATRIX_DENSE_REAL)
      rerror (THIS_SOLVER ": " RLAB_ERROR_RHS_FUNC_MUST_RETURN_MDR);

    retm = ent_data(rent);

    if (SIZE(retm) != odr_np * odr_nq)
      rerror (THIS_SOLVER ": " RLAB_ERROR_JACOBIAN_FUNC_INCORRECT_DIM);

    //
    // fjacb[i;j;k] -> dfdb(x,b)[j;k]
    //
    for (j=0; j<odr_nq; j++)
      for (k=0; k<odr_np; k++)
        fjacb[i + odr_n * j + odr_nq * odr_n * k] = mdr0 (retm, j, k);

    ent_Clean (rent);
  }

  return 1;
}


//
// fortran: dfdx = dfdx(x,b)
//
int
odrdxfn (double *x, double *b, double *fjacd)
{
  int i, j, k;

  // do not have jacobian dfdb - find it numerically
  if (odr_fjacd == 0)
    return 1;

  //
  // b is the same for each datum x
  //
  MDPTR (b1) = (void *) b;

  //
  // go over each x, evaluate dfdx and store result
  //
  for (i=0; i<odr_n; i++)
  {
    Ent *rent = 0;
    MDR *retm = 0;

    //
    // x1 = x[i;]
    //
    for (j=0; j<odr_m; j++)
      MdrV0 (x1, j) = x[i + odr_n * j];

    //
    // do my bidding
    //
    rent = ent_call_rlab_script_2args(dxfname, x1_ent, b1_ent);

    if (ent_type(rent)!=MATRIX_DENSE_REAL)
      rerror (THIS_SOLVER ": " RLAB_ERROR_RHS_FUNC_MUST_RETURN_MDR);

    retm = ent_data(rent);

    if (SIZE(retm) != odr_nq * odr_m)
      rerror (THIS_SOLVER ": " RLAB_ERROR_JACOBIAN_FUNC_INCORRECT_DIM);

    //
    // fjacd[i;j;k] -> dfdx(x,b)[j;k]
    //
    for (j = 0; j < odr_nq; j++)
      for (k = 0; k < odr_m; k++)
        fjacd[i + odr_n * j + odr_nq * odr_n * k] = mdr0 (retm, j, k);

    ent_Clean (rent);
  }

  return 1;
}
