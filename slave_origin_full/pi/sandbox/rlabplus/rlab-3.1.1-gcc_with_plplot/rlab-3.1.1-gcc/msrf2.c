// /* msrf2.c Matrix Sparse Real Functions ... */

//  This file is a part of RLaB ("Our"-LaB)
//
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

#include "rlab.h"
#include "mem.h"
#include "msr.h"
#include "mscf2.h"
#include "mdr.h"
#include "mdrf1.h"
#include "mdc.h"
#include "mds.h"
#include "btree.h"
#include "util.h"
#include "symbol.h"
#include "sort.h"
#include "mathl.h"

#include "msrf1.h"

#include "class.h"

#include "fi.h"

#include "sparskit.h"

#include <stdio.h>
#include <math.h>

#ifdef HAVE_SUITESPARSE
// #include "clibs/umfpack/UMFPACK/Include/umfpack.h"
#include <suitesparse/umfpack.h>
//
// sparse solver UMFPACK, parameters:
//
void *umfpack_Numeric;
MSR *umfpack_a;
int umfpack_n = 0;
#endif

//
// choice of method for solve function on sparse real matrices:
//      0 - UMFPACK  (umf)
//      1 - SuperLU  (slu)
//      2 - SPARSKIT (spk)
int msr_Solve_method = 2;

//
// sparse solver with SPARSKIT2, parameters:
//
int spkit_ipre = 3;   // ILUD
int spkit_iter = 6;   //GMRES
int spkit_maxits = 500;
int spkit_precstat = 3;   // preconditioning status: 0-none,1-left,2-right,3-both
int spkit_convcrit = 3;   // iterative solver chooses best convergence criterion
int spkit_m = 10;   // m dimension of subspace used in some iterative solvers (GMRES)
double spkit_abserr = 1e-6;
double spkit_relerr = 0;

Ent *
ent_sparse_tol (int nargs, Datum args[])
{
  Ent *e1 = 0, *e2 = 0, *e3 = 0, *e4 = 0;
  double ae = -1, re = -1;
  int maxit = 0, ccrit = 0;
  if (nargs == 0 || nargs > 4)
  {
    fprintf (stdout,
       "tol: A member of the function list  sparams. Sets the absolute and\n");
    fprintf (stdout,
       "tol: relative error, convergence criterion and the maximum number\n");
    fprintf (stdout,
       "tol: of iterations to reach the tolerances for a sparse solver\n");
    fprintf (stdout, "tol: from the SPARKIT package. Format:\n");
    fprintf (stdout,
       "tol:   sparams.tol( abserr, relerr, convcrit, maxit ),\n");
    fprintf (stdout,
       "tol: The admissible values for convergence criterion are:\n");
    fprintf (stdout,
       "tol: convcrit = 1 (residual norm), 2 (approximate solution),\n");
    fprintf (stdout, "tol: 3 (solver chooses).\n");
    fprintf (stdout,
       "tol: Currently, abserr = %g, relerr = %g, convcrit = %i and maxit = %i.\n",
       spkit_abserr, spkit_relerr, spkit_convcrit, spkit_maxits);
    rerror ("No tolerances given!");
  }
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) == MATRIX_DENSE_REAL)
    ae = class_double (e1);

  if (nargs >= 2)
  {
    e2 = bltin_get_ent (args[1]);
    if (ent_type (e2) == MATRIX_DENSE_REAL)
      re = class_double (e2);
  }
  if (nargs >= 3)
  {
    e3 = bltin_get_ent (args[2]);
    if (ent_type (e3) == MATRIX_DENSE_REAL)
      ccrit = (int) class_double (e3);
  }
  if (nargs == 4)
  {
    e4 = bltin_get_ent (args[3]);
    if (ent_type (e4) == MATRIX_DENSE_REAL)
      maxit = (int) class_double (e4);
  }

  if (ae > 0)
    spkit_abserr = ae;
  if (re >= 0)
    spkit_relerr = re;
  if (ccrit > 0 && ccrit < 4)
    spkit_convcrit = ccrit;
  if (maxit > 0)
    spkit_maxits = maxit;

  ent_Clean (e1);
  ent_Clean (e2);
  ent_Clean (e3);
  ent_Clean (e4);

  return ent_Create_Rlab_Success();
}


Ent *
ent_sparse_prec (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0;
  int ipre = 0;
  if (nargs == 0 || nargs > 2)
  {
    fprintf (stdout,
	     "precond: A member of the function list  sparams. Chooses an\n");
    fprintf (stdout,
	     "precond: preconditioner for a sparse real iterator from the.\n");
    fprintf (stdout, "precond: SPARKIT package. Format:\n");
    fprintf (stdout, "precond:   sparams.precond( ipre, stat ),\n");
    fprintf (stdout, "precond: Available preconditioners are:\n");
    fprintf (stdout,
	     "precond: ipre = 1, for ILUT, 2, for ILUTP, 3, for ILUD,\n");
    fprintf (stdout,
	     "precond: 4, for ILUDP, and, 5, for ILUK, all being some modification.\n");
    fprintf (stdout,
	     "precond: of an incomplete LU-factorization class of methods.\n");
    fprintf (stdout,
	     "precond: Preconditioning status is used later in iteration scheme\n");
    fprintf (stdout,
	     "precond: and the admissible values are: stat=0 (none), 1 (left),\n");
    fprintf (stdout, "precond: 2 (right), and 3 (both).\n");
    fprintf (stdout, "precond: Current choice is the preconditioner %i.\n",
	     spkit_ipre);
    fprintf (stdout, "precond: with preconditioning status %i.\n",
	     spkit_precstat);
    rerror ("No preconditioner given!");
  }
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) == MATRIX_DENSE_REAL)
    ipre = (int) class_double (e1);
  if (ipre > 0 && ipre < 6)
    spkit_ipre = ipre;
  if (ent_type (e1) != UNDEF)
    ent_Clean (e1);
  //
  if (nargs == 2)
  {
    e2 = bltin_get_ent (args[1]);
    if (ent_type (e2) == MATRIX_DENSE_REAL)
      ipre = (int) class_double (e2);
    if (ipre > 0 && ipre < 4)
      spkit_precstat = ipre;
    if (ent_type (e2) != UNDEF)
      ent_Clean (e2);
  }

  return ent_Create_Rlab_Success();
}


Ent *
ent_sparse_iter (int nargs, Datum args[])
{
  Ent *e1 = 0, *e2 = 0;
  int iter = 0, m = 0;
  if (nargs == 0)
  {
    fprintf (stdout,
	     "iterator: A member of the function list  sparams. Chooses an\n");
    fprintf (stdout,
	     "iterator: iterator and the dimension of linear subspace used in\n");
    fprintf (stdout,
	     "iterator: the iterative scheme, for a real sparse linear problem\n");
    fprintf (stdout, "iterator: solver from the SPARKIT package. Format:\n");
    fprintf (stdout, "iterator:   sparams.iterator( iter, m ),\n");
    fprintf (stdout, "iterator: Available iterators are:\n");
    fprintf (stdout,
	     "iterator: iter = 1 for CG, 2 for CGNR, 3 for BCG, 4 for BCGSTAB,\n");
    fprintf (stdout,
	     "iterator: 5 for TFQMR, 6 for GMRES(m), 7 for FGMRES(m), 8 for\n");
    fprintf (stdout, "iterator: DQGMRES(m) and, 9 for DBCG.\n");
    fprintf (stdout,
	     "iterator: Current choice is the iterator %i with m = %i.\n",
	     spkit_iter, spkit_m);
    rerror ("No iterator given!");
  }
  // get the name of the iterator
  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) == MATRIX_DENSE_REAL)
    iter = (int) class_double (e1);
  if (iter > 0 && iter < 10)
    spkit_iter = iter;
  if (ent_type (e1) != UNDEF)
    ent_Clean (e1);
  // get m
  if (nargs >= 2)
  {
    e2 = bltin_get_ent (args[1]);
    if (ent_type (e2) == MATRIX_DENSE_REAL)
      m = (int) class_double (e2);
    if (m > 0 && m < 100)
      spkit_m = m;
    if (ent_type (e2) != UNDEF)
      ent_Clean (e2);
  }

  return ent_Create_Rlab_Success();
}


MDR *
sparskit_msr_Solve (MSR * a, MDR * b, char *type)
{
  int nmax, nzmax, lwk, ierr, lfil, nwk;
  double fpar[16];
  int ipar[16];
  int i, j, n;
  //
  msr_Detect_Inf (a);
  msr_Detect_NaN (a);
  mdr_Detect_Inf (b);
  mdr_Detect_Nan (b);
  if (MNR(a) != MNC(a))
    rerror ("solve[spk]: requires square coefficient matrix");
  if ((a->nr) != (b->nrow))
    rerror ("solve[spk]: RHS row dim. must match LHS row dim.");

  lfil = 7;
  n = MNR(a);
  nmax = n;
  nzmax = a->nnz * (2 * lfil);
  nwk = nzmax;

  double *au, *wk;
  int *jau, *ju, *iw;

  au = (double *) GC_malloc (nwk * sizeof (double));
  jau = (int *) GC_malloc (nwk * sizeof (int));
  ju = (int *) GC_malloc ((n + 1) * sizeof (int));
  iw = (int *) GC_malloc (nmax * 2 * sizeof (int));

  MDR *x=mdr_Create (n, MNC(b));

  //
  // Loop over the columns of B, getting a solution for each.
  //
  for (i=0; i < MNC(b); i++)
  {
    // rhs = mdr_PartitionCol (b, i);
    double *rhs = &MdrV0(b,i*n);
    double *sol = &MdrV0(x,i*n);
    for (j = 0; j < n; j++)
      sol[j] = rand ();

    // setup iterator
//     m = spkit_m;
    switch (spkit_iter)
    {
      case 2:
        lwk = 5 * n;    // CGNR
        break;

      case 3:
        lwk = 7 * n;
        break;

      case 4:
        lwk = 8 * n;    // BCGSTAB
        break;

      case 5:
        lwk = 11 * n;   // TFQMR
        break;

      case 6:
        lwk = (n + 3) * (spkit_m + 2);  // GMRES
        break;

      case 7:
        lwk = (n + 3) * (spkit_m + 2) + (spkit_m + 1) * spkit_m / 2;  // FGMRES
        break;

      case 8:
        lwk = n + spkit_m * (2 * n + 4);  // DQGMRES
        break;

      case 9:
        lwk = 11 * n;   // DBCG
        break;

      default:
        // case 1:
        lwk = 5 * n;    // CG
        break;
    }

    if (n+2 > lwk)
      wk = (double *) GC_malloc ((n + 2) * sizeof (double));
    else
      wk = (double *) GC_malloc (lwk * sizeof (double));

    //
    // go through preconditioner
    //
    signal (SIGINT, intcatch);
    switch (spkit_ipre)
    {
//       case 1:
//         ILUT (&n, a->d, a->ja, a->ia, &lfil, &spkit_abserr, au, jau, ju, &nwk, wk,
//               iw, &ierr);
//         break;

      case 2:
        ILUTP (&n, a->d, a->ja, a->ia, &lfil, &spkit_abserr, au, jau, ju, &nwk,
               wk, iw, &ierr);
        break;

      case 3:
        ILUD (&n, a->d, a->ja, a->ia, &lfil, &spkit_abserr, au, jau, ju, &nwk, wk,
              iw, &ierr);
        break;

      case 4:
        ILUDP (&n, a->d, a->ja, a->ia, &lfil, &spkit_abserr, au, jau, ju, &nwk,
               wk, iw, &ierr);
        break;

      case 5:
        ILUK (&n, a->d, a->ja, a->ia, &lfil, &spkit_abserr, au, jau, ju, &nwk, wk,
              iw, &ierr);
        break;

      default:
        // case 1:
        ILUT (&n, a->d, a->ja, a->ia, &lfil, &spkit_abserr, au, jau, ju, &nwk, wk,
              iw, &ierr);
    }
    signal (SIGINT, intcatch_wait);

    if (ierr != 0)
    {
      fprintf (stdout, "spkit: Preconditioner reports an error %i.\n", ierr);
      fprintf (stdout,
         "spkit: Cannot calculate incomplete LU factorization.\n");
      rerror ("spkit: Terrible internal error!");
    }

    //
    // set up iterator
    //
    ipar[0] = 0;
    ipar[1] = spkit_precstat;	// preconditioning status: 0-none,1-left,2-right,3-both
    ipar[2] = spkit_convcrit;	// convergence criterion: 1,2,3
    ipar[3] = lwk;
    ipar[4] = spkit_m;
    ipar[5] = spkit_maxits;
    fpar[0] = spkit_abserr;
    fpar[1] = spkit_relerr;

    // iterate
    RUNRC (&n, rhs, sol, ipar, fpar, wk, a->d, a->ja, a->ia, au, jau, ju,
	   &spkit_iter);

    GC_free (wk);

    if (ipar[0] == -2)
    {
      fprintf (stdout, "spkit: solver did not have enough memory!\n");
      mdr_Destroy(x);
      x = mdr_Create (0, 0);
      break;
    }
    else if (ipar[0] == -3)
    {
      fprintf (stdout, "spkit: solver decided to quit!\n");
      mdr_Destroy(x);
      x = mdr_Create (0, 0);
      break;
    }
  }

  GC_free (au);
  GC_free (jau);
  GC_free (ju);
  GC_free (iw);

  return (x);
}


#ifdef HAVE_SUPERLU

//
// sparse solver SuperLU, parameters:
//
double superlu_diagpivothresh = 1.0;
int superlu_permc_spec = 3;
int superlu_n = 0;
double superlu_droptol = 0; // useless

Ent *
ent_sparse_dpt (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0;
  double dpt = 1;
  int ips = -1;

  if (nargs != 1 || nargs != 2)
  {
    fprintf (stdout,
             "dpt: A member of the function list  sparams. Sets the diagonal pivot\n");
    fprintf (stdout,
             "dpt: threshold (dpt) and the column ordering for SuperLU sparse solver.\n");
    fprintf (stdout, "dpt: Format:\n");
    fprintf (stdout, "dpt:   sparams.dpt( dpt, ord ),\n");
    fprintf (stdout,
             "dpt: Admissible values for ordering are: ord = 0 for natural ordering,\n");
    fprintf (stdout,
             "dpt: 1 for minimum degree on a'*a, 2 for minimum degree on a'+a, and \n");
    fprintf (stdout,
             "dpt: 3 for approximate minimum degree for unsymmetric matrices.\n");
    fprintf (stdout, "dpt: Currently, dpt = %g and ", superlu_diagpivothresh);
    fprintf (stdout, "ord = %i\n", superlu_permc_spec);
    rerror ("One or two parameters required");
  }

  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) == MATRIX_DENSE_REAL)
  {
    dpt = class_double (e1);
    if (dpt > 0)
      superlu_diagpivothresh = dpt;
  }

  if (nargs == 2)
  {
    e2 = bltin_get_ent (args[1]);
    if (ent_type (e2) == MATRIX_DENSE_REAL)
      ips = (int) class_double (e2);
    if (ips > -1 && ips < 4)
      superlu_permc_spec = ips;
  }

  ent_Clean (e1);
  ent_Clean (e2);
  return ent_Create_Rlab_Success();
}

// SuperLU interface
//   Tzong-Shuoh Yang (yang@isec.com) 10/3/96 solve()
//                                    1/14/97 factor(), backsub()
//   Modified for SuperLU-1.0, Ian Searle, 2/9/97.
//   Modified for SuperLU-3.0, Marijan Kostrun, 4/9/05
#include <superlu_enum_consts.h>
#include <slu_ddefs.h>
#define  Matrix SuperMatrix

extern void Destroy_CompCol_Matrix (Matrix *);
extern void Destroy_SuperNode_Matrix (Matrix *);
extern void Destroy_CompCol_Permuted (Matrix *);
extern void StatInit (SuperLUStat_t *);
extern void StatFree (SuperLUStat_t *);
extern void sp_preorder (superlu_options_t *, Matrix *, int *, int *, Matrix *);
extern int sp_ienv (int);
extern void dgstrs (trans_t, Matrix *, Matrix *, int *, int *,
                    Matrix *, SuperLUStat_t *, int *);
extern void dgscon (char *, Matrix *, Matrix *,
                    double, double *, SuperLUStat_t *, int *);
extern void dgssv (superlu_options_t *, Matrix *, int *, int *, Matrix *,
                   Matrix *, Matrix *, SuperLUStat_t *, int *);
extern double dlangs (char *, Matrix *);
extern void get_perm_c (int, Matrix *, int *);

SuperMatrix *superlu_L;
SuperMatrix *superlu_U;
int *superlu_perm_r;
int *superlu_perm_c;
//    } factors_t;
//factors_t *LUfactors;

/* this routine comes from SuperLU/MATLAB/mexsuperlu.c */
static void
    LUextract (Matrix * L, Matrix * U, double *Lval, int *Lrow, int *Lcol,
               double *Uval, int *Urow, int *Ucol, int *snnzL, int *snnzU)
{
  int i, j, k;
  int upper;
  int fsupc, istart, nsupr;
  int lastl = 0, lastu = 0;
  SCformat *Lstore;
  NCformat *Ustore;
  double *SNptr;

  Lstore = L->Store;
  Ustore = U->Store;
  Lcol[0] = 0;
  Ucol[0] = 0;

  /* for each supernode */
  for (k = 0; k <= Lstore->nsuper; ++k)
  {

    fsupc = L_FST_SUPC (k);
    istart = L_SUB_START (fsupc);
    nsupr = L_SUB_START (fsupc + 1) - istart;
    upper = 1;

    /* for each column in the supernode */
    for (j = fsupc; j < L_FST_SUPC (k + 1); ++j)
    {
      SNptr = &((double *) Lstore->nzval)[L_NZ_START (j)];

      /* Extract U */
      for (i = U_NZ_START (j); i < U_NZ_START (j + 1); ++i)
      {
        Uval[lastu] = ((double *) Ustore->nzval)[i];
        /* Matlab doesn't like explicit zero. */
        if (Uval[lastu] != 0.0)
          Urow[lastu++] = U_SUB (i);
      }
      for (i = 0; i < upper; ++i)
      {       /* upper triangle in the supernode */
        Uval[lastu] = SNptr[i];
        /* Matlab doesn't like explicit zero. */
        if (Uval[lastu] != 0.0)
          Urow[lastu++] = L_SUB (istart + i);
      }
      Ucol[j + 1] = lastu;

      /* Extract L */
      Lval[lastl] = 1.0;  /* unit diagonal */
      Lrow[lastl++] = L_SUB (istart + upper - 1);
      for (i = upper; i < nsupr; ++i)
      {
        Lval[lastl] = SNptr[i];
        /* Matlab doesn't like explicit zero. */
        if (Lval[lastl] != 0.0)
          Lrow[lastl++] = L_SUB (istart + i);
      }
      Lcol[j + 1] = lastl;

      ++upper;

    }       /* for j ... */

  }       /* for k ... */

  *snnzL = lastl;
  *snnzU = lastu;
}

//
// LU factorization of A. It performs the following steps:
//
//   1. Permute columns of A, forming A*Pc, where Pc is a permutation matrix.
//
//   2. Factor A as Pr*A*Pc=L*U using LU decomposition with Pr determined
//      by partial pivoting.
//
//   3. Calculate rcond = 1/cond(A)
//
// This routine is very much like SuperLU's dgssv. Dgssv wasn't sufficient
// because it doesn't produce a condition number estimate.
//
// modified by kmk for superlu 3.0
static void
    msr_super_fact (Matrix * A, int *perm_c, int *perm_r, int *etree,
                    Matrix * L, Matrix * U, int *info, double *rcond)
{

  char norm[1];
  superlu_options_t slu_options;
  set_default_options (&slu_options);
//   NCformat *Astore;
  Matrix AC;      /* Matrix postmultiplied by Pc */
  int lwork = 0;
  SuperLUStat_t slu_stat;
  // Set defaults
  double anorm;
  int panel_size;   /* panel size */
  int relax;      /* no of columns in a relaxed snodes */
  //
  // Options
  //
  slu_options.Fact = DOFACT;
  panel_size = sp_ienv (1);
  relax = sp_ienv (2);
  StatInit (&slu_stat);
//   Astore = A->Store;
  // Test the input parameters
  *info = 0;
  if (A->nrow != A->ncol || A->nrow < 0)
  {
    *info = -1;
    return;
  }
  //
  // Do the pre-ordering, given a permutation vector.
  //
  sp_preorder (&slu_options, A, perm_c, etree, &AC);
  //
  // Compute the LU factorization of A.
  //
  dgstrf (&slu_options, &AC, relax, panel_size, etree,
           NULL, lwork, perm_c, perm_r, L, U, &slu_stat, info);
  if (*info == 0)
  {
    // estimate the reciprocal of the condition number of A.
    *norm = '1';
    anorm = dlangs (norm, A);
    dgscon (norm, L, U, anorm, rcond, &slu_stat, info);
  }
  // Clean Up.
  Destroy_CompCol_Permuted (&AC);
  StatFree (&slu_stat);
}

// ***************************************************************
// * solve sparse problem using SuperLU
// ***************************************************************
MDR * superlu_msr_Solve (MSR * a, MDR * b, char *type_in)
{
  int i;
  SuperMatrix A, B, L, U;
  MDR *x = mdr_Float_BF (b);
  int n, ne, info, *ia, *ja;
  int *perm_r, *perm_c;
  SuperLUStat_t slu_stat;
  superlu_options_t options;
  set_default_options (&options);
  char *type, *dummy = "no";

  if (type_in)
    type = type_in;
  else
    type = dummy;
  //
  // Check for un-manageable matrix elements.
  //
  msr_Detect_Inf (a);
  msr_Detect_NaN (a);
  mdr_Detect_Inf (b);
  mdr_Detect_Nan (b);

  //
  // Check for poorly posed problem.
  //
  if (a->nr != a->nc)
    rerror ("solve[slu]: requires square coefficient matrix");
  if (a->nr != MNR (b))
    rerror ("solve[slu]: RHS row dim. must match LHS row dim.");
  n = a->nr;
  ne = a->nnz;

  //
  // Adjust indices, rlab uses 1,1, while SuperLU is 0,0
  //
  ia = a->ia;
  ja = a->ja;
  for (i = 0; i < ne; i++)
    ja[i]--;
  for (i = 0; i <= n; i++)
    ia[i]--;

  //
  // Setup the SuperLU [A] from [a]: MSR * a is row ordered.
  //
  dCreate_CompRow_Matrix (&A, n, n, ne, a->d, a->ja, a->ia, SLU_NR, SLU_D, SLU_GE);
  //
  // Setup the SuperLU RHS from [x] (a copy of [b]).
  //
  dCreate_Dense_Matrix (&B, n, MNC (b), MDRPTR(x), n, SLU_DN, SLU_D, SLU_GE);
  //
  // Allocate permutation arrays
  //
  perm_r = (int *) GC_MAIOP (n * sizeof (int));
  perm_c = (int *) GC_MAIOP (n * sizeof (int));
  //
  // Solve
  //
  StatInit (&slu_stat);
  //
  // Options
  //
  options.DiagPivotThresh = superlu_diagpivothresh;

  options.SymmetricMode = NO;
  if (type)
  {
    if (*type == 's' || *type == 'S')
    {
      options.SymmetricMode = YES;
    }
  }
  signal (SIGINT, intcatch);
  dgssv (&options, &A, perm_c, perm_r, &L, &U, &B, &slu_stat, &info);
  signal (SIGINT, intcatch_wait);
  if (info != 0)
  {
    fprintf (stderr, "error in solve[slu], info = %d\n", info);
    rerror ("SuperLU error");
  }
  //
  // index adjusting
  //
  for (i = 0; i < ne; i++)
    ja[i]++;
  for (i = 0; i <= n; i++)
    ia[i]++;

  //
  // SuperLU created these, so SuperLU must destroy them.
  //
  StatFree (&slu_stat);
  Destroy_CompCol_Matrix (&U);
  Destroy_SuperNode_Matrix (&L);

  GC_FREE (perm_c);
  GC_FREE (perm_r);
  free (A.Store);
  free (B.Store);
  return (x);
}


//
// spsolve with SuperLU
//
MDR *
superlu_msr_spsolve (MSR * a, MDR * b)
{
  int i;
  Matrix A, B;
  int *etree, n, ne, info, *ia, *ja;
  double rcond;
  trans_t trans;
  SuperLUStat_t slu_stat;

  //
  if (a->nr * a->nc != 0)
  {
    if (a->nr != a->nc)
      rerror ("spsolve: requires square coefficient matrix");

    //
    // spsolve(a) or spsolve(a,b): factorize  a  and store it
    //
#ifdef HAVE_SUITESPARSE
    if (umfpack_n != 0)
    {
      // erase old umfpack stuff if it exists
      umfpack_di_free_numeric (&umfpack_Numeric);
      umfpack_n = 0;
      msr_Destroy (umfpack_a);
    }
#endif
    if (superlu_n != 0)
    {
      superlu_n = 0;

      // the code was called earlier. Destroy old memory
      Destroy_CompCol_Matrix (superlu_U);
      Destroy_SuperNode_Matrix (superlu_L);
      GC_FREE (superlu_L);
      GC_FREE (superlu_U);
      GC_FREE (superlu_perm_c);
      GC_FREE (superlu_perm_r);
    }
    msr_Detect_Inf (a);
    msr_Detect_NaN (a);
    //
    n = a->nr;
    ne = a->nnz;
    superlu_n = n;
    //
    dCreate_CompCol_Matrix (&A, n, n, ne, a->d, a->ja, a->ia, SLU_NC, SLU_D,
                            SLU_GE);

    //
    ia = a->ia;
    ja = a->ja;
    for (i = 0; i < ne; i++)
      ja[i]--;
    for (i = 0; i <= n; i++)
      ia[i]--;

    //
    superlu_L = (Matrix *) GC_MALLOC (sizeof (Matrix));
    superlu_U = (Matrix *) GC_MALLOC (sizeof (Matrix));

    //
    etree = (int *) GC_MAIOP (n * sizeof (int));
    superlu_perm_r = (int *) GC_MAIOP (n * sizeof (int));
    superlu_perm_c = (int *) GC_MAIOP (n * sizeof (int));
    get_perm_c (superlu_permc_spec, &A, superlu_perm_c);

    //
    signal (SIGINT, intcatch);
    msr_super_fact (&A,
                    superlu_perm_c, superlu_perm_r, etree, superlu_L, superlu_U,
                    &info, &rcond);
    signal (SIGINT, intcatch_wait);
    if (info != 0)
    {
      if (info <= n)
        fprintf (stderr, "Diagonal zero, equation %d\n", info);
      else
        fprintf (stderr, "Memory allocation failure\n");
      rerror ("SuperLU error");
    }

    //
    for (i = 0; i < ne; i++)
      ja[i]++;
    for (i = 0; i <= n; i++)
      ia[i]++;

    //
    GC_FREE (etree);
    free (A.Store);
  }

  if (!b)
    return mdr_CreateScalar (1.0);
  if (b->nrow * b->ncol == 0)
    return mdr_CreateScalar (1.0);

  //
  // spsolve(b) or spsolve(,b)
  //
  mdr_Detect_Inf (b);
  mdr_Detect_Nan (b);
  MDR *x = mdr_Float_BF (b);
  StatInit (&slu_stat);
  //
  dCreate_Dense_Matrix (&B, superlu_n, MNC (b), MDRPTR(x), superlu_n, SLU_DN, SLU_D,
			SLU_GE);
  //
  trans = TRANS;
  signal (SIGINT, intcatch);
  dgstrs (trans,
	  superlu_L, superlu_U, superlu_perm_c, superlu_perm_r,
	  &B, &slu_stat, &info);
  signal (SIGINT, intcatch_wait);
  Destroy_SuperMatrix_Store (&B);
  StatFree (&slu_stat);
  if (info != 0)
  {
    fprintf (stderr, "dgstrs reports info = %d\n", info);
    rerror ("spsolve[slu] error");
  }
  return (x);
}


/* **************************************************************
 * Factor a sparse real matrix.
 * This function is an interface to SuperLU routines.
 * The return value is a list - n, l, u, perm_r, perm_c and rcond.
 * See msr_Backsub() for the meaning of the return list.
 * ************************************************************** */

Btree *
msr_SpFactor (MSR * a)
{
  Btree *bt;
  int i;
  Matrix A, L, U;
  NCformat *Astore, *Ustore;
  SCformat *Lstore;
  MSR *acopy;
  int *etree, n, ne, *perm_ctmp, *perm_r, info, *ia, *ja;
  double rcond;
  Ent *N, *PERM_C, *PERM_R, *RCOND;
  MDR *mperm_r, *mperm_c;

  /* Create the list that will contain the results. */
  bt = btree_Create ();

  /* Check data */
  msr_Detect_Inf (a);
  msr_Detect_NaN (a);

  if (a->nr != a->nc)
    rerror ("factor: requires square coefficient matrix");

  n = a->nr;
  ne = a->nnz;

  /* superLU uses column-wise storage */
  acopy = msr_Transpose (a);

  perm_r = (int *) GC_MAIOP (n * sizeof (int));

  A.Stype = SLU_NC;
  A.Dtype = SLU_D;
  A.Mtype = SLU_GE;
  A.nrow = n;
  A.ncol = n;
  A.Store = GC_MALLOC (sizeof (NCformat));
  Astore = (NCformat *) A.Store;
  Astore->nnz = ne;
  Astore->nzval = acopy->d;
  Astore->rowind = acopy->ja;
  Astore->colptr = acopy->ia;

  etree = (int *) GC_MAIOP (A.ncol * sizeof (int));

  /* index adjustment (SuperLU index starts from 0) */
  ia = acopy->ia;
  ja = acopy->ja;
  for (i = 0; i < ne; i++)
    ja[i]--;
  for (i = 0; i <= n; i++)
    ia[i]--;

  /* copy permutation vector or use mmd. */
  perm_ctmp = (int *) GC_MAIOP (n * sizeof (int));
  get_perm_c (superlu_permc_spec, &A, perm_ctmp);

  /* LU decomposition */
  signal (SIGINT, intcatch);
  msr_super_fact (&A, perm_ctmp, perm_r, etree, &L, &U, &info, &rcond);
  signal (SIGINT, intcatch_wait);
  if (info != 0)
  {
    if (info <= n)
      fprintf (stderr, "Diagonal zero, equation %d\n", info);
    else
      fprintf (stderr, "Memory allocation failure\n");
    rerror ("SuperLU error");
  }
  msr_Destroy (acopy);
  GC_FREE (etree);
  GC_FREE (A.Store);

  /*
   * Hook the results into the list.
   * Items that need to be saved:
   * n, L, U, perm_c, perm_r, rcond
   */

  /* Load n */
  N = ent_Create ();
  ent_SetType (N, DOUBLE);
  ent_double (N) = (double) n;
  install (bt, ("n"), N);

  /* convert L, U to rlab sparse matrices */
  {
    Ent *l, *u;
    MSR *rL, *rU, *rLt, *rUt;
    int nnzL, nnzU, *ia, *ja;

    Lstore = (SCformat *) L.Store;
    rL = msr_Create (n, n);
    nnzL = Lstore->nnz;
    rL->d = (double *) GC_MAIOP (nnzL * sizeof (double));
    rL->ja = (int *) GC_MAIOP (nnzL * sizeof (int));
    rL->ia = (int *) GC_MAIOP ((n + 1) * sizeof (int));

    Ustore = (NCformat *) U.Store;
    rU = msr_Create (n, n);
    nnzU = Ustore->nnz;
    rU->d = (double *) GC_MAIOP (nnzU * sizeof (double));
    rU->ja = (int *) GC_MAIOP (nnzU * sizeof (int));
    rU->ia = (int *) GC_MAIOP ((n + 1) * sizeof (int));
    LUextract (&L, &U, rL->d, rL->ja, rL->ia,
         rU->d, rU->ja, rU->ia, &nnzL, &nnzU);
    Destroy_CompCol_Matrix (&U);
    Destroy_SuperNode_Matrix (&L);

    /* indecies adjusting */
    rL->nnz = nnzL;
    ia = rL->ia;
    ja = rL->ja;
    for (i = 0; i < nnzL; i++)
      ja[i]++;
    for (i = 0; i <= n; i++)
      ia[i]++;
    rU->nnz = nnzU;
    ia = rU->ia;
    ja = rU->ja;
    for (i = 0; i < nnzU; i++)
      ja[i]++;
    for (i = 0; i <= n; i++)
      ia[i]++;

    /* transpose and store them */
    rLt = msr_Transpose (rL);
    msr_Destroy (rL);
    l = ent_Create ();
    ent_SetType (l, MATRIX_SPARSE_REAL);
    ent_data (l) = rLt;
    install (bt, ("l"), l);
    rUt = msr_Transpose (rU);
    msr_Destroy (rU);
    u = ent_Create ();
    ent_SetType (u, MATRIX_SPARSE_REAL);
    ent_data (u) = rUt;
    install (bt, ("u"), u);
  }

  /* Load perm_r */
  mperm_r = mdr_Create (1, n);
  for (i = 0; i < n; i++)
    MdrV0 (mperm_r, i) = (double) (perm_r[i] + 1);  /* starts from 1 */
  GC_FREE (perm_r);

  PERM_R = ent_Create ();
  ent_SetType (PERM_R, MATRIX_DENSE_REAL);
  ent_data (PERM_R) = mperm_r;
  install (bt, ("perm_r"), PERM_R);

  /* Load perm_c */
  mperm_c = mdr_Create (1, n);
  for (i = 0; i < n; i++)
    MdrV0 (mperm_c, i) = (double) (perm_ctmp[i] + 1);
  GC_FREE (perm_ctmp);
  PERM_C = ent_Create ();
  ent_SetType (PERM_C, MATRIX_DENSE_REAL);
  ent_data (PERM_C) = mperm_c;
  install (bt, ("perm_c"), PERM_C);

  /* Load rcond */
  RCOND = ent_Create ();
  ent_SetType (RCOND, DOUBLE);
  ent_double (RCOND) = rcond;
  install (bt, ("rcond"), RCOND);

  return (bt);
}

/* **************************************************************
 * Built-in interface to SpFactor. This is more complex/complete
 * than the normal factor() interface. There are extra arguments that
 * allow tailoring of the solution.
 * ************************************************************** */

Ent *
SpFactor_BF (int nargs, Datum args[])
{
  void *rval = 0;
  Ent *e1=0, *rent;
  rent = 0;

  /* Check n_args */
  if (nargs != 1)
  {
    fprintf (stdout,
       "spfactor: Factorize sparse linear matrix 'A' using Super-LU\n");
    fprintf (stdout, "spfactor: Format:\n");
    fprintf (stdout, "spfactor:   spfactor(A),\n");
    fprintf (stdout, "spfactor: where 'A' is a sparse matrix.\n");
    fprintf (stdout, "spfactor: See also  sparams.dpt(..) .\n");
    rerror ("single argument required");
  }

  e1 = bltin_get_ent (args[0]);
  if (!
      (ent_type (e1) == MATRIX_SPARSE_REAL
       || ent_type (e1) == MATRIX_SPARSE_COMPLEX))
    rerror ("spfactor[slu]: A must be sparse");
  if (ent_type (e1) == MATRIX_SPARSE_REAL)
    rval = (void *) msr_SpFactor (ent_data (e1));
  else if (ent_type (e1) == MATRIX_SPARSE_COMPLEX)
    rval = (void *) msc_SpFactor (ent_data (e1));
  else
    rerror ("spfactor: solution method not supported");

  ent_Clean (e1);

  rent = ent_Create ();
  ent_data (rent) = rval;
  ent_type (rent) = BTREE;
  return (rent);
}

Btree *
msr_Factor (MSR * a, char *type)
{
  Btree *x;
  x = msr_SpFactor (a);
  return (x);
}



#endif /* ifdef HAVE_SUPERLU */


#ifdef HAVE_SUITESPARSE

//
// sparse solver with UMFPACK
//
MDR *
umfpack_msr_Solve (MSR * a, MDR * b, char *type)
{
  int i, j, n;
  MDR *btmp, *x, *xtmp;
  //
  // Initialize UMFPACK
  //
  double Control[UMFPACK_CONTROL], Info[UMFPACK_INFO];
  void *Symbolic, *Numeric;

  msr_Detect_Inf (a);
  msr_Detect_NaN (a);
  mdr_Detect_Inf (b);
  mdr_Detect_Nan (b);
  if (a->nr != a->nc)
    rerror ("solve[umf]: requires square coefficient matrix");
  if (a->nr != MNR (b))
    rerror ("solve[umf]: RHS row dim. must match LHS row dim.");
  n = a->nr;
//   ne = a->nnz;

  // this is silly but I don;t know better: indexing is offset by one
  for (i = 0; i < a->nr + 1; i++)
    --a->ia[i];
  for (j = 0; j < a->nnz; j++)
    --a->ja[j];

  //
  // Factor the coefficient matrix: first symbolic then numeric
  //
  (void) umfpack_di_symbolic
      (n, n, a->ia, a->ja, a->d, &Symbolic, Control, Info);
  (void) umfpack_di_numeric
      (a->ia, a->ja, a->d, Symbolic, &Numeric, Control, Info);

  //
  // Solve Ax = b
  //
  xtmp = mdr_Create (MNR (b), 1);
  x = mdr_Create (MNR (b), MNC (b));
  //
  // Loop over the columns of B, getting a solution for each.
  //
  for (i = 1; i <= MNC (b); i++)
  {
    // Grab the right column of B.
    btmp = mdr_PartitionCol (b, i);
    // Solve, but remember A is transposed (MSR is row compressed format)
    (void) umfpack_di_solve
      (UMFPACK_Aat, a->ia, a->ja, a->d, MDRPTR(xtmp), MDRPTR(btmp), Numeric, NULL, NULL);
    // kill bill, I mean btmp
    mdr_Destroy (btmp);
    // Load result into X
    for (j = 1; j <= MNR (b); j++)
      Mdr1 (x, j, i) = Mdr1 (xtmp, j, 1);
  }
  mdr_Destroy (xtmp);

  umfpack_di_free_symbolic (&Symbolic);
  umfpack_di_free_numeric (&Numeric);

  // this is silly but I don;t know better: indexing is offset by one
  for (i = 0; i < a->nr + 1; i++)
    ++a->ia[i];
  for (j = 0; j < a->nnz; j++)
    ++a->ja[j];

  return (x);
}

MDR *
umfpack_msr_Det (MSR * a)
{
  int i, j, n;
  MDR *x;
  //
  // Initialize UMFPACK
  //
  double Control[UMFPACK_CONTROL], Info[UMFPACK_INFO];
  void *Symbolic, *Numeric;
  double Mx, Ex;

  msr_Detect_Inf (a);
  msr_Detect_NaN (a);
  if (a->nr != a->nc)
    rerror ("det[umf]: requires square matrix");
  n = a->nr;
//   ne = a->nnz;

  // this is silly but I don;t know better: indexing is offset by one
  for (i = 0; i < a->nr + 1; i++)
    --a->ia[i];
  for (j = 0; j < a->nnz; j++)
    --a->ja[j];

  //
  // Factor the coefficient matrix: first symbolic then numeric
  //
  signal (SIGINT, intcatch);
  (void) umfpack_di_symbolic (n, n, a->ia, a->ja, a->d, &Symbolic, Control,
            Info);
  signal (SIGINT, intcatch_wait);

  signal (SIGINT, intcatch);
  (void) umfpack_di_numeric (a->ia, a->ja, a->d, Symbolic, &Numeric, Control,
           Info);
  signal (SIGINT, intcatch_wait);

  // Find the determinant
  signal (SIGINT, intcatch);
  (void) umfpack_di_get_determinant (&Mx, &Ex, Numeric, Info);
  signal (SIGINT, intcatch_wait);

  umfpack_di_free_symbolic (&Symbolic);

  // this is silly but I don;t know better: indexing is offset by one
  for (i = 0; i < a->nr + 1; i++)
    ++a->ia[i];
  for (j = 0; j < a->nnz; j++)
    ++a->ja[j];

  if (Info[UMFPACK_STATUS] == UMFPACK_OK)
    x = mdr_CreateScalar (Mx * pow (10, Ex));
  else if (Info[UMFPACK_STATUS] == UMFPACK_WARNING_singular_matrix)
    x = mdr_CreateScalar (0.0);
  else if (Info[UMFPACK_STATUS] == UMFPACK_WARNING_determinant_underflow ||
     Info[UMFPACK_STATUS] == UMFPACK_WARNING_determinant_overflow)
  {
    x = mdr_Create (1, 2);
    Mdr1 (x, 1, 1) = Mx;
    Mdr1 (x, 1, 2) = Ex;
  }
  else
    x = mdr_Create (0, 0);
  return (x);
}

//
// spsolve with UMFPACK
//
MDR *
umfpack_msr_spsolve (MSR * a, MDR * b)
{
  int i, j, n;
  MDR *x;
  //
  // Initialize UMFPACK
  //
  double Control[UMFPACK_CONTROL], Info[UMFPACK_INFO];
  void *Symbolic;

  mdr_Detect_Inf (b);
  mdr_Detect_Nan (b);
  msr_Detect_Inf (a);
  msr_Detect_NaN (a);

  if (!b)
    return mdr_CreateScalar (1.0);
  if (b->nrow * b->ncol == 0)
    return mdr_CreateScalar (1.0);

  if (a->nr * a->nc == 0)
    return mdr_CreateScalar (2.0);

  if (a->nr != a->nc)
    rerror ("spsolve[umf]: requires square coefficient matrix");
  if ((a->nr) != b->nrow)
    rerror ("spsolve[umf]: In 'Ax=B', 'A' and 'B' must have same number of rows!");
#ifdef HAVE_SUPERLU
  if (superlu_n != 0)
  {
    superlu_n = 0;
    // the code was called earlier. Destroy old memory
    Destroy_CompCol_Matrix (superlu_U);
    Destroy_SuperNode_Matrix (superlu_L);
    GC_FREE (superlu_L);
    GC_FREE (superlu_U);
    GC_FREE (superlu_perm_c);
    GC_FREE (superlu_perm_r);
  }
#endif
  if (umfpack_n != 0)
  {
    umfpack_n = 0;
    umfpack_di_free_numeric (&umfpack_Numeric);
    msr_Destroy (umfpack_a);
  }

  umfpack_a = msr_Copy (a);
  n = umfpack_a->nr;
  umfpack_n = n;
  // this is silly but I don;t know better: indexing is offset by one
  for (i = 0; i < n + 1; i++)
    --umfpack_a->ia[i];
  for (j = 0; j < a->nnz; j++)
    --umfpack_a->ja[j];

  // Factor the coefficient matrix: first symbolic then numeric
  signal (SIGINT, intcatch);
  (void) umfpack_di_symbolic (umfpack_n, umfpack_n,
                              umfpack_a->ia, umfpack_a->ja, umfpack_a->d,
                              &Symbolic, Control, Info);
  (void) umfpack_di_numeric (umfpack_a->ia, umfpack_a->ja, umfpack_a->d,
                             Symbolic, &umfpack_Numeric, Control, Info);
  signal (SIGINT, intcatch_wait);
  //
  umfpack_di_free_symbolic (&Symbolic);

  //
  // Solve Ax = b
  //
  x = mdr_Create ((b->nrow), (b->ncol));

  //
  // Loop over the columns of B, getting a solution for each.
  //
  for (i = 0; i < (b->ncol); i++)
  {
    // Grab the right column of B.
    double *btmp = &MdrV0(b,i*MNR(b));
    double *xtmp = &MdrV0(x,i*MNR(x));

    // Solve, but remember A is transposed (MSR is row compressed format)
    signal (SIGINT, intcatch);
    (void) umfpack_di_solve (UMFPACK_Aat, umfpack_a->ia, umfpack_a->ja,
                             umfpack_a->d, xtmp, btmp, umfpack_Numeric,
                             NULL, NULL);
    signal (SIGINT, intcatch_wait);
  }

  return (x);
}
#endif


//
// spsolve
//
Ent *
SpSolve_BF (int nargs, Datum args[])
{
  void *rval = 0;
  void *(*vfptr) () = NULL;
  int rtype = 0;
  Ent *e1=0, *e2=0, *rent;
  rent = 0;

  // Check n_args
  if (nargs != 1 && nargs != 2)
  {
    fprintf (stdout,
             "spsolve: Solve real sparse linear problem  A*x=B using"
#ifdef HAVE_SUPERLU
             "Super-LU"
#endif
             " or UMFPACK"
             "\n");
    fprintf (stdout,
             "spsolve: packages by storing internally the result of LU factorization and\n");
    fprintf (stdout, "spsolve: reusing it thereafter.\n");
    fprintf (stdout, "spsolve: Format:\n");
    fprintf (stdout, "spsolve:   x=spsolve(a,b), or\n");
    fprintf (stdout, "spsolve:   spsolve(a); x=spsolve(b);\n");
    fprintf (stdout,
             "spsolve: where 'a' is a real sparse matrix, and 'b' is the dense RHS matrix.\n");
    fprintf (stdout,
             "spsolve: In the first case the linear problem is solved immediately for a\n");
    fprintf (stdout,
             "spsolve: given 'a' and 'b'. In the second case the matrix 'a' is first\n");
    fprintf (stdout,
             "spsolve: LU-factorized, and the result of factorization is stored internally.\n");
    fprintf (stdout,
             "spsolve: This is then used to solve the system for any dense matrix 'b'.\n");
    rerror ("one or two arguments required");
  }

  // Get each argument.
  e1 = bltin_get_ent (args[0]);
  if (nargs == 2)
    e2 = bltin_get_ent (args[1]);

  // figure out what package to use
#ifdef  HAVE_SUITESPARSE
  if (msr_Solve_method == 0)
    vfptr = (VFPTR) umfpack_msr_spsolve;
  else
#endif
#ifdef HAVE_SUPERLU
  if (msr_Solve_method == 1)
    vfptr = (VFPTR) superlu_msr_spsolve;
#endif

  if (!vfptr)
  {
    ent_Clean (e1);
    ent_Clean (e2);
    return ent_Create_Rlab_Failure ();
  }

  // then use it
  if (nargs == 1 && ent_type (e1) == MATRIX_SPARSE_REAL)
  {
    // spsolve(a)
    rval = (*vfptr) (ent_data (e1), mdr_Create (0, 0));
    rtype = MATRIX_DENSE_REAL;
  }
  else if (nargs == 1 && ent_type (e1) == MATRIX_DENSE_REAL)
  {
    // spsolve(b)
    rval = (*vfptr) (msr_Create (0, 0), ent_data (e1));
    rtype = MATRIX_DENSE_REAL;
  }
  else if (nargs == 2 && ent_type (e1) == MATRIX_SPARSE_REAL
	   && ent_type (e2) == MATRIX_DENSE_REAL)
  {
    // spsolve(a,b)
    rval = (*vfptr) (ent_data (e1), ent_data (e2));
    rtype = MATRIX_DENSE_REAL;
  }
  else
    rerror ("spsolve: solution method not supported");

  ent_Clean (e1);
  ent_Clean (e2);

  rent = ent_Create ();
  ent_data (rent) = rval;
  ent_type (rent) = rtype;
  return (rent);
}



/* **************************************************************
 * Finish the solution, via back-substitution of a system of
 * linear equations. The result of a previous factorization of
 * the coefficient matrix is provided as input, along with the
 * right-hand side(s).
 *
 * This function is an interface to SuperLU routines.
 *    A*x = b
 *    P*A*Q' = L*U, P=perm_r, Q=perm_c
 * procedure:
 *    1. L*y = P*b, solve for y
 *    2. U*z = y,   solve for z
 *    3. x = Q'*z
 * T. S. Yang (yang@isec.com) 1/14/97
 * ************************************************************** */

MDR *
msr_Backsub (Btree * bt, MDR * b)
{
  MSR *L = 0, *U = 0;
  MDR *sol;
  ListNode *ltmp;
  Ent *etmp;
  int i, j, kk, n = 0, m, *perm_c = 0, *perm_r = 0;
  MDR *tmp;
  double *x, *y, *z, *r;

  mdr_Detect_Inf (b);
  mdr_Detect_Nan (b);

  /* Extract the factorization data from the table. */

  if ((ltmp = btree_FindNode (bt, "n")))
  {
    etmp = var_ent (ltmp);
    if (ent_type (etmp) != DOUBLE)
      rerror ("backsub: terrible input list error 1");
    n = (int) ent_double (etmp);
    /* Check dimemsions */
    if (MNR (b) != n)
      rerror ("backsub: RHS row dim. must match LHS row dim.");
  }
  else
  {
    rerror ("backsub: terrible input list error 2");
  }

  /* no. of rhs */
  m = MNC (b);

  if ((ltmp = btree_FindNode (bt, "perm_c")))
  {
    etmp = var_ent (ltmp);
    if (ent_type (etmp) != MATRIX_DENSE_REAL)
      rerror ("backsub: terrible input list error 3");
    tmp = ent_data (etmp);
    perm_c = (int *) GC_MAIOP (n * sizeof (int));
    for (i = 0; i < n; i++)
    {
      perm_c[i] = (int) MdrV0 (tmp, i) - 1;
    }
  }
  else
  {
    rerror ("backsub: terrible input list error 4");
  }

  if ((ltmp = btree_FindNode (bt, "perm_r")))
  {
    etmp = var_ent (ltmp);
    if (ent_type (etmp) != MATRIX_DENSE_REAL)
      rerror ("backsub: terrible input list error 5");
    tmp = ent_data (etmp);
    perm_r = (int *) GC_MAIOP (n * sizeof (int));
    for (i = 0; i < n; i++)
    {
      perm_r[i] = (int) MdrV0 (tmp, i) - 1;
    }
  }
  else
  {
    rerror ("backsub: terrible input list error 6");
  }

  if ((ltmp = btree_FindNode (bt, "l")))
  {
    etmp = var_ent (ltmp);
    if (ent_type (etmp) != MATRIX_SPARSE_REAL)
      rerror ("backsub: terrible input list error");
    L = (MSR *) ent_data (etmp);
  }
  else
    rerror ("backsub: terrible input list error");

  if ((ltmp = btree_FindNode (bt, "u")))
  {
    etmp = var_ent (ltmp);
    if (ent_type (etmp) != MATRIX_SPARSE_REAL)
      rerror ("backsub: terrible input list error");
    U = (MSR *) ent_data (etmp);
  }
  else
    rerror ("backsub: terrible input list error");

  sol = mdr_Float_BF (b);   /* the solution */
  r = (double *) GC_MAIOP (n * sizeof (double));  /* work array */
  z = r;

  /* Backsub */

  signal (SIGINT, intcatch);
  for (kk = 0; kk < m; kk++)  /* loop thru no. of right-hand-side vectors */
  {
    double t;
    int nnzU;

    x = &MdrV0(sol,kk*n);
    y = x;

    /* r = P*b */
    for (i = 0; i < n; i++)
      r[perm_r[i]] = x[i];

    /* solve L*y = r */
    y[0] = r[0];
    for (i = 1; i < n; i++)
    {
      t = r[i];
      for (j = L->ia[i]; j < (L->ia[i + 1] - 1); j++)
      {
        t -= L->d[j - 1] * y[L->ja[j - 1] - 1];
      }
      y[i] = t;
    }
    /* solve U*z = y */
    nnzU = U->nnz;
    z[n - 1] = y[n - 1] / U->d[nnzU - 1];
    for (i = n - 2; i >= 0; i--)
    {
      t = y[i];
      for (j = U->ia[i]; j < (U->ia[i + 1] - 1); j++)
      {
        t -= U->d[j] * z[U->ja[j] - 1];
      }
      z[i] = t / U->d[U->ia[i] - 1];
    }

    /* x = Q'*z */
    for (i = 0; i < n; i++)
      x[i] = z[perm_c[i]];
  }
  signal (SIGINT, intcatch_wait);

  /* clean up */
  GC_FREE (r);
  GC_FREE (perm_r);
  GC_FREE (perm_c);

  return (sol);
}


#ifdef HAVE_METIS2
/* **************************************************************
 * Fill reducing ordering of a sparse matrix via Metis2 lib
 * ************************************************************** */

static int adjacency (MSR * a, int *xadj, int *adjncy);
extern int OMETIS (int *, int *, int *, int *, int *, int *, int *);

MDR *
msr_SpOrder (MSR * m)
{
  int i, n, numbering;
  int options[5], *p, *invp, *xadj, *adjncy;
  int retval;
  MDR *P = 0;
  MSR *mt;

  /* First, check for empty matrix. */
  if ((m->nr == 0 && m->nc == 0) || m->nnz == 0)
  {
    rerror ("sporder: cannot order empty matrix");
  }
  else
  {
    /* First, transpose the matrix so the ia and ja are correct. */
    mt = msr_Transpose (m);
    n = mt->nr;

    xadj = (int *) GC_MAIOP ((n + 1) * sizeof (int));
    adjncy = (int *) GC_MAIOP (mt->nnz * sizeof (int));

    if ((retval = adjacency (mt, xadj, adjncy)) == 0)
    {
      fprintf (stderr, "sporder: possible an empty row in sparse matrix");
      rerror ("sporder: cannot convert sparse matrix into graph");
    }

    /* ometis options array. */
    options[0] = 0;

    /* We are using ?-style array indexing. */
    numbering = 1;

    /* Permutation arrays. */
    p = (int *) GC_MAIOP (n * sizeof (int));
    invp = (int *) GC_MAIOP (n * sizeof (int));

    OMETIS (&n, xadj, adjncy, options, &numbering, p, invp);

    /* Now, load P. */
    P = mdr_Create (1, n);
    for (i = 0; i < n; i++)
    {
      MdrV0 (P, i) = (double) invp[i];
    }

    /* Clean up. */
    msr_Destroy (mt);
    GC_FREE (xadj);
    GC_FREE (adjncy);
    GC_FREE (p);
    GC_FREE (invp);
  }
  return (P);
}
#endif /* HAVE_METIS2 */

#ifdef HAVE_METIS3
/* **************************************************************
 * Fill reducing ordering of a sparse matrix via Metis3 lib
 * ************************************************************** */

static int adjacency (MSR * a, int *xadj, int *adjncy);
typedef int idxtype;
extern void METIS_EdgeND (int *, idxtype *, idxtype *, int *, int *, idxtype *,
			  idxtype *);
extern void METIS_NodeND (int *, idxtype *, idxtype *, int *, int *, idxtype *,
			  idxtype *);

MDR *
msr_SpOrder (MSR * m)
{
  int i, n, numbering;
  int options[8], *p, *invp, *xadj, *adjncy;
  int retval;
  MDR *P = 0;
  MSR *mt;

  /* First, check for empty matrix. */
  if ((m->nr == 0 && m->nc == 0) || m->nnz == 0)
  {
    rerror ("sporder: cannot order empty matrix");
  }
  else
  {
    /* First, transpose the matrix so the ia and ja are correct. */
    mt = msr_Transpose (m);
    n = mt->nr;

    xadj = (int *) GC_MAIOP ((n + 1) * sizeof (int));
    adjncy = (int *) GC_MAIOP (mt->nnz * sizeof (int));

    if ((retval = adjacency (mt, xadj, adjncy)) == 0)
    {
      fprintf (stderr, "sporder: possible an empty row in sparse matrix");
      rerror ("sporder: cannot convert sparse matrix into graph");
    }

    /* ometis options array. */
    options[0] = 0;

    /* We are using ?-style array indexing. */
    numbering = 1;

    /* Permutation arrays. */
    p = (int *) GC_MAIOP (n * sizeof (int));
    invp = (int *) GC_MAIOP (n * sizeof (int));

    METIS_NodeND (&n, xadj, adjncy, &numbering, options, p, invp);

    /* Now, load P. */
    P = mdr_Create (1, n);
    for (i = 0; i < n; i++)
    {
      MdrV0 (P, i) = (double) invp[i];
    }

    /* Clean up. */
    msr_Destroy (mt);
    GC_FREE (xadj);
    GC_FREE (adjncy);
    GC_FREE (p);
    GC_FREE (invp);
  }
  return (P);
}
#endif /* HAVE_METIS3 */

// /* **************************************************************
//  * Compute the adjacency lists for a sparse matrix.
//  *
//  * The diagonal elements are implicitly referenced.
//  * An error will occur if the graph is not complete (zero row
//  * or column).
//  *
//  * Original written by T.S. Yang, modified by Ian Searle.
//  * ************************************************************** */
// #if 0
// static int
// adjacency (MSR * a, int *xadj, int *adjncy)
// {
//   /* create the adjacency of sparse matrix 'a' */
//   int i, j, ii, jj, kk, nz, ntz, n, nr;
//
//   nr = a->nr;
//   ntz = 0;
//
//   /*
//    * loop over each row...
//    */
//
//   kk = 0;
//   xadj[0] = 1;
//
//   for (i = 0; i < nr; i++)
//   {
//     n = a->ia[i + 1] - a->ia[i];	/* Number in row. */
//     ii = a->ia[i];		/* Index of row-start. */
//
//     /* Check for empty row. */
//     if (n == 0)
//       return (0);
//
//     /* Go through each row. */
//     for (j = 0; j < n; j++)
//     {
//       nz = 0;
//       jj = a->ja[ii + j - 1];	/* Column index */
//       if ((i + 1) != jj)	/* Check for diagonal element. */
//       {
// 	adjncy[kk++] = jj;
// 	nz++;
//       }
//     }
//     xadj[i + 1] = xadj[i] + nz;
//   }
//
//   return (1);
// }
// #else
// static int
// adjacency (MSR * a, int *xadj, int *adjncy)
// {
//   /* create the adjacency of sparse matrix 'a' */
//   int i, j, ii, jj, nz, ntz, n, nr;
//
//   nr = a->nr;
//   ntz = 0;
//   for (i = 0; i < nr; i++)
//   {
//     n = a->ia[i + 1] - a->ia[i];
//     ii = a->ia[i];
//     nz = 0;
//
//     /* Check for empty row. */
//     if (n == 0)
//       return (0);
//
//     for (j = 0; j < n; j++)
//     {
//       jj = a->ja[ii + j - 1] - 1;
//       if (i != jj)
//       {
// 	adjncy[ntz + nz] = jj + 1;
// 	nz++;
//       }
//     }
//     xadj[i] = ntz + 1;
//     ntz += nz;
//   }
//   xadj[a->nr] = ntz + 1;
//   return (1);
// }
// #endif

/* **************************************************************
 * Force a sparse real matrix into a column vector format.
 * Retain sparsity.
 * ************************************************************** */

MSR *
msr_ReshapeCol (MSR * m)
{
  int i;
  MSR *col;

  /* First, check for empty matrix. */
  if ((m->nr == 0 && m->nc == 0) || m->nnz == 0)
  {
    col = msr_Create (m->nr, m->nc);
    msr_Setup (col, m->nnz);
  }
  else
  {
    /* Transpose, so we can pull the columns out of the rows. */
    col = msr_Transpose (m);

    /* JA gets set to all ones. */
    for (i = 0; i < col->nnz; i++)
    {
      col->ja[i] = 1;
    }

    /*
     * IA gets modified, 1:NR.
     * First we must re-allocate.
     */

    GC_FREE (col->ia);
    col->ia = (int *) GC_MAIOP ((m->nr * m->nc + 1) * sizeof (int));
    col->nr = m->nr * m->nc;
    col->nc = 1;

    for (i = 0; i < col->nr + 1; i++)
    {
      col->ia[i] = i + 1;
    }
  }
  return (col);
}

// ***************************************************************
// * Solve a set of equations with a real sparse coefficient matrix
// * All we do here is solve the general sparse matrix.
// ***************************************************************

Ent *
ent_sparse_realsolv (int nargs, Datum args[])
{
  // call parameters:
  //  e1   - rng name
  //  [e2] - seed
  Ent *e1=0;
  MDS *x1;
  char sparam[32] = { '\0' };
  if (nargs == 0)
  {
    fprintf (stdout, "realsolv: A member of the function list  sparams.\n");
    fprintf (stdout,
	     "realsolv: Chooses a method for solving sparse real linear systems.\n");
    fprintf (stdout, "realsolv: Format:\n");
    fprintf (stdout, "realsolv:   sparams.realsolv(\"method\" ),\n");
    fprintf (stdout, "realsolv: Available methods/packages are:\n");
#ifdef HAVE_SUITESPARSE
    fprintf (stdout, "realsolv: (0) \"UMFPACK\", \"umfpack\", \"umf\"\n");
#endif
#ifdef HAVE_SUPERLU
    fprintf (stdout, "realsolv: (1) \"SuperLU\", \"superlu\", \"slu\",\n");
#endif
    fprintf (stdout, "realsolv: (2) \"SPARSKIT\", \"sparskit\", \"spk\".\n");
    fprintf (stdout, "realsolv: Currently the method ");
#ifdef HAVE_SUITESPARSE
    if (msr_Solve_method == 0)
      fprintf (stdout, "UMFPACK");
    else
#endif
#ifdef HAVE_SUPERLU
    if (msr_Solve_method == 1)
      fprintf (stdout, "SuperLU");
    else
#endif
      fprintf (stdout, "SPARSKIT");

    fprintf (stdout, " is being used!\n");
    rerror ("No method given!");
  }
  //
  // get name of the generator
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) == MATRIX_DENSE_STRING)
  {
    x1 = ent_data (e1);
    sprintf (sparam, "%s", Mds1 (x1, 1, 1));
  }
  else
    rerror ("sparams.realsolv: Unknown method!");

#ifdef HAVE_SUITESPARSE
  if (strcmp (sparam, "umfpack") == 0)
    msr_Solve_method = 0;
  else if (strcmp (sparam, "UMFPACK") == 0)
    msr_Solve_method = 0;
  else if (strcmp (sparam, "umf") == 0)
    msr_Solve_method = 0;
  else
#endif
#ifdef HAVE_SUPERLU
  if (strcmp (sparam, "SuperLU") == 0)
    msr_Solve_method = 1;
  else if (strcmp (sparam, "superlu") == 0)
    msr_Solve_method = 1;
  else if (strcmp (sparam, "SUPERLU") == 0)
    msr_Solve_method = 1;
  else if (strcmp (sparam, "slu") == 0)
    msr_Solve_method = 1;
  else
#endif
    msr_Solve_method = 2;

  ent_Clean (e1);

  return ent_Create_Rlab_Success ();
}


MDR *
msr_Solve (MSR * a, MDR * b, char *type)
{
  void *(*vfptr) ();

#ifdef HAVE_SUITESPARSE
  if (msr_Solve_method == 0)
    vfptr = (VFPTR) umfpack_msr_Solve;
  else
#endif
#ifdef HAVE_SUPERLU
  if (msr_Solve_method == 1)
    vfptr = (VFPTR) superlu_msr_Solve;
  else
#endif
    vfptr = (VFPTR) sparskit_msr_Solve;
  return (*vfptr) (a, b, type);
}
