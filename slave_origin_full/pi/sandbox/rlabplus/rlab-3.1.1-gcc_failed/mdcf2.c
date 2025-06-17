/* mdcf2.c Matrix Dense Complex Functions */

/*  This file is a part of RLaB ("Our"-LaB)
   Copyright (C) 1995  Ian R. Searle

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.

   See the file ./COPYING
   ********************************************************************** */

#include "rlab.h"
#include "mem.h"
#include "ent.h"
#include "btree.h"
#include "symbol.h"
#include "mdr.h"
#include "mdc.h"
#include "mdcf1.h"
#include "util.h"
#include "complex.h"
#include "mathl.h"

#include "fi.h"
#include "lp.h"
#include "blas.h"

#include <stdio.h>
#include <math.h>

/*
 * Create Infs, and NaNs
 */
#include "mathl.h"

MDC *mdc_SolveEq_GE (MDC * a, MDC * b);
MDC *mdc_SolveEq_SYM (MDC * a, MDC * b);
MDC *mdc_SolveEq_SYM2 (MDC * a, MDC * b);
void mdc_Svd (MDC * M, MDC ** rsv, MDC ** lsv, MDR ** sigma, int flag);
double *mdc_Norm (MDC * m, char *type);
double matrix_Rcond (MDC * m);


/* **************************************************************
 * Solve the matrix equation: Ax = b
 * ************************************************************** */

MDC *mdc_SolveEq (MDC * a, MDC * b);
MDC *mdc_SolveEq_GE (MDC * a, MDC * b);
MDC *mdc_SolveEq_SYM (MDC * a, MDC * b);
MDC *mdc_SolveEq_SYM2 (MDC * a, MDC * b);

/*
 * Cover package for builtin function solve().
 */

MDC *
mdc_Solve (MDC * a, MDC * b, char *type)
{
  MDC *m = 0;

  if (type != 0)
  {
    if (!strncmp ("S", type, 1))
    {
      m = mdc_SolveEq_SYM (a, b);
    }
    else if (!strncmp ("s", type, 1))
    {
      m = mdc_SolveEq_SYM (a, b);
    }
    else if (!strncmp ("G", type, 1))
    {
      m = mdc_SolveEq_GE (a, b);
    }
    else if (!strncmp ("g", type, 1))
    {
      m = mdc_SolveEq_GE (a, b);
    }
  }
  else
  {
    m = mdc_SolveEq (a, b);
  }
  return (m);
}

MDC *
mdc_SolveEq (MDC * a, MDC * b)
{
  MDC *m;
  if (mdc_IsSymmetric (a))
  {
    m = mdc_SolveEq_SYM (a, b);
  }
  else
  {
    m = mdc_SolveEq_GE (a, b);
  }

  return (m);
}

/* **************************************************************
 * The general case...
 * ************************************************************** */

MDC *
mdc_SolveEq_GE (MDC * a, MDC * b)
{
  int info, lda, ldb, mm, n, nrhs, *ipiv;
  int trans, lwork, norm;
  double *anorm, rcond, *rwork;
  MDC *A, *B;
  MDC *work;

  mdc_Detect_Inf (a);
  mdc_Detect_Nan (a);

  mdc_Detect_Inf (b);
  mdc_Detect_Nan (b);

  lda = mm = MNR (a);
  n = MNC (a);
  lwork = mm;
  trans = (int) 'N';
  norm = (int) '1';

  if (mm != n)
    rerror ("matrix must be square for solve()");

  if (mm != MNR (b))
    rerror ("RHS row dim. must match LHS row dim.");

  ipiv = (int *) GC_malloc_atomic_ignore_off_page (sizeof (int) * mm);
  A = mdc_Copy (a);
  B = mdc_Copy (b);
  work = mdc_Create (1, 4 * lwork);

  signal (SIGINT, intcatch);
  XGETRF (&mm, &n, MDCPTR (A), &lda, ipiv, &info);
  signal (SIGINT, intcatch_wait);

  if ((int) info < 0)
    rerror ("bad argument to LAPACK ZGETRF");
  if ((int) info > 0)
    rerror ("matrix is singular");

  /*
   * Check input for ill-conditioning.
   */

  rwork = (double *) GC_malloc_atomic_ignore_off_page (sizeof (double) * 2 * n);
  anorm = mdc_Norm (a, "1");

  signal (SIGINT, intcatch);
  XGECON (&norm, &n, MDCPTR (A), &lda, anorm, &rcond, MDCPTR (work),
	  rwork, &info);
  signal (SIGINT, intcatch_wait);
  if (rcond <= DBL_EPSILON)
    warning_1 ("WARNING, ill-conditioned input");

  GC_FREE (anorm);
  nrhs = MNC (b);
  ldb = MNR (b);

  signal (SIGINT, intcatch);
  XGETRS (&trans, &n, &nrhs, MDCPTR (A), &lda, ipiv, MDCPTR (B), &ldb, &info);
  signal (SIGINT, intcatch_wait);

  if ((int) info < 0)
    rerror ("Bad argument(s) to LAPACK DGETRS");

  GC_FREE (ipiv);
  mdc_Destroy (A);
  mdc_Destroy (work);
  GC_FREE (rwork);

  return (B);
}

/* **************************************************************
 * The symmetric case...
 * ************************************************************** */

MDC *
mdc_SolveEq_SYM (MDC * a, MDC * b)
{
  int info, lda, ldb, mm, n, nrhs, *ipiv;
  int lwork;
  int uplo;
  double *anorm, rcond;
  MDC *A, *B, *work;

  mdc_Detect_Inf (a);
  mdc_Detect_Nan (a);

  mdc_Detect_Inf (b);
  mdc_Detect_Nan (b);

  n = lda = mm = MNR (a);
  uplo = (int) 'L';

  /*
   * Try and pick a good NB, without ILAENV.
   */

  if (n < 100)
    lwork = n;
  else
    lwork = 64 * n;

  if (mm != n)
    rerror ("solve: matrix must be square");

  if (mm != MNR (b))
    rerror ("solve: RHS row dim. must match LHS row dim.");

  ipiv = (int *) GC_malloc_atomic_ignore_off_page (sizeof (int) * n);
  A = mdc_Copy (a);
  B = mdc_Copy (b);
  work = mdc_Create (1, lwork);

  signal (SIGINT, intcatch);
  XHETRF (&uplo, &n, MDCPTR (A), &lda, ipiv, MDCPTR (work), &lwork, &info);
  signal (SIGINT, intcatch_wait);

  if ((int) info < 0)
    rerror ("bad argument to LAPACK ZHETRF");
  if ((int) info > 0)
    rerror ("matrix is singular");

  /*
   * Check input for ill-conditioning.
   */

  mdc_Destroy (work);
  work = mdc_Create (1, 2 * n);

  anorm = mdc_Norm (a, "1");
  signal (SIGINT, intcatch);
  XHECON (&uplo, &n, MDCPTR (A), &lda, ipiv, anorm, &rcond,
	  MDCPTR (work), &info);
  signal (SIGINT, intcatch_wait);

  if (rcond <= DBL_EPSILON)
    warning_1 ("WARNING, ill-conditioned input");

  GC_FREE (anorm);
  nrhs = MNC (b);
  ldb = MNR (b);

  signal (SIGINT, intcatch);
  XHETRS (&uplo, &n, &nrhs, MDCPTR (A), &lda, ipiv, MDCPTR (B), &ldb, &info);
  signal (SIGINT, intcatch_wait);

  if ((int) info < 0)
    rerror ("Bad argument(s) to LAPACK DSYTRS or ZHETRS");

  GC_FREE (ipiv);
  mdc_Destroy (A);
  mdc_Destroy (work);

  return (B);
}

MDC *
mdc_SolveEq_SYM2 (MDC * a, MDC * b)
{
  int info, lda, ldb, mm, n, nrhs, *ipiv;
  int lwork;
  int uplo;
  MDC *A, *B, *work;

  mdc_Detect_Inf (a);
  mdc_Detect_Nan (a);

  mdc_Detect_Inf (b);
  mdc_Detect_Nan (b);

  n = lda = mm = MNR (a);
  uplo = (int) 'U';
  nrhs = MNC (b);
  ldb = MAX (1, n);

  /*
   * Try and pick a good NB, without ILAENV.
   */

  if (n < 100)
    lwork = n;
  else
    lwork = 64 * n;

  if (mm != n)
    rerror ("solve: matrix must be square");

  if (mm != MNR (b))
    rerror ("solve: RHS row dim. must match LHS row dim.");

  A = mdc_Copy (a);
  B = mdc_Copy (b);
  ipiv = (int *) GC_malloc_atomic_ignore_off_page (sizeof (int) * n);
  work = mdc_Create (1, lwork);

  signal (SIGINT, intcatch);
  XHESV (&uplo, &n, &nrhs, MDCPTR (A), &lda, ipiv, MDCPTR (B),
	 &ldb, MDCPTR (work), &lwork, &info);
  signal (SIGINT, intcatch_wait);

  if ((int) info < 0)
    rerror ("bad argument to LAPACK ZHESV");
  if ((int) info > 0)
    rerror ("matrix is singular");

  GC_FREE (ipiv);
  mdc_Destroy (work);
  mdc_Destroy (A);

  return (B);
}

/* **************************************************************
 * Compute the recipricol of the condition number.
 * ************************************************************** */

double *
mdc_Rcond (MDC * m)
{
  int lda, n, info, *ipiv, norm;
  double *anorm, *rwork;
  double *rcond;
  MDC *A, *work;

  mdc_Detect_Inf (m);
  mdc_Detect_Nan (m);

  lda = MNR (m);
  n = MIN(MNC (m), lda);

  norm = (int) '1';

  work = mdc_Create (1, 4 * n);
  rwork = (double *) GC_malloc_atomic_ignore_off_page (sizeof (double) * 2 * n);
  ipiv = (int *) GC_malloc_atomic_ignore_off_page (sizeof (int) * lda);
  A = mdc_Copy (m);

  anorm = mdc_Norm (m, "1");
  rcond = (double *) GC_MALLOC (sizeof (double));

  signal (SIGINT, intcatch);
  XGETRF (&lda, &n, MDCPTR (A), &lda, ipiv, &info);
  XGECON (&norm, &n, MDCPTR (A), &lda, anorm, rcond, MDCPTR (work),
	  rwork, &info);
  signal (SIGINT, intcatch_wait);

  if ((int) info < 0)
    rerror ("illegal argument to LAPACK ZGECON");

  mdc_Destroy (work);
  GC_FREE (rwork);
  GC_FREE (ipiv);
  mdc_Destroy (A);

  return (rcond);
}

/* **************************************************************
 * Compute the norm of rows of matrix
 * ************************************************************** */
#define THIS_SOLVER "rownorm"
MDR *
mdc_RowNorm (MDC * m, char *type)
{
  double work;
  int nrow, ncol, one=1;
  int i, j, itype;
  Complex *vec=0;

  mdc_Detect_Inf (m);
  mdc_Detect_Nan (m);

  /* Argument error checking */
  if (!strcmp (type, "1"));
  else if (!strcmp (type, "M"));
  else if (!strcmp (type, "m"));
  else if (!strcmp (type, "O"));
  else if (!strcmp (type, "o"));
  else if (!strcmp (type, "I"));
  else if (!strcmp (type, "i"));
  else if (!strcmp (type, "F"));
  else if (!strcmp (type, "f"));
  else if (!strcmp (type, "E"));
  else if (!strcmp (type, "e"));
  else
    rerror ("norm: incorrect STRING specifier");

  nrow = (int) MNR (m);
  ncol = (int) MNC (m);

  MDR * norm = mdr_Create(nrow,1);

  vec  = (Complex *) GC_malloc_atomic_ignore_off_page((size_t) (ncol * sizeof (Complex)));

  itype = (int) (type[0] - '0' + '0');

  for (i=0; i<nrow; i++)
  {
    for (j=0;j<ncol;j++)
      vec[j] = Mdc0(m,i,j);

    signal (SIGINT, intcatch);
    MdrV0(norm,i) = XLANGE (&itype, &one, &ncol, vec, &one, &work);
    signal (SIGINT, intcatch_wait);
  }

  GC_FREE (vec);

  return (norm);
}

double *
mdc_Norm (MDC * m, char *type)
{
  F_DOUBLE *norm, *work;
  int lda, nrow, ncol;
  int itype;
  MDC *lsv, *rsv;
  MDR *sigma;

  norm = (F_DOUBLE *) GC_MALLOC (sizeof (F_DOUBLE));
  *norm = 0.0;

  mdc_Detect_Inf (m);
  mdc_Detect_Nan (m);

  /* Argument error checking */
  if (!strcmp (type, "1"));
  else if (!strcmp (type, "M"));
  else if (!strcmp (type, "m"));
  else if (!strcmp (type, "O"));
  else if (!strcmp (type, "o"));
  else if (!strcmp (type, "I"));
  else if (!strcmp (type, "i"));
  else if (!strcmp (type, "F"));
  else if (!strcmp (type, "f"));
  else if (!strcmp (type, "E"));
  else if (!strcmp (type, "e"));
  else if (!strcmp (type, "2"));
  else
    rerror ("norm: incorrect STRING specifier");

  nrow = (int) MNR (m);
  ncol = (int) MNC (m);
  lda = nrow;

  if (strcmp (type, "2"))
  {
    work = (double *) GC_malloc_atomic_ignore_off_page
    ((size_t) (nrow * sizeof (double)));

    itype = (int) (type[0] - '0' + '0');
    signal (SIGINT, intcatch);
    *norm = XLANGE (&itype, &nrow, &ncol, MDCPTR (m), &lda, work);
    signal (SIGINT, intcatch_wait);
    GC_FREE (work);
  }
  else
  {
    /* Compute the matrix 2-norm */
    mdc_Svd (m, &rsv, &lsv, &sigma, 3);

    /*
     * Return the largest singular value s[1]
     */

    *norm = (F_DOUBLE) Mdr1 (sigma, 1, 1);
    mdc_Destroy (rsv);
    mdc_Destroy (lsv);
    mdr_Destroy (sigma);
  }
  return (norm);
}

/* **************************************************************
 * Compute the matrix Pnorm.
 * ************************************************************** */

double *
mdc_PNorm (MDC * m, double p)
{
  double *pnorm;
  double sum = 0.0;
  int i, size;

  if (MNR (m) != 1 && MNC (m) != 1)
    rerror ("cannot compute P-norm of a matrix");

  pnorm = (double *) GC_MALLOC (sizeof (double));
  size = MNR (m) * MNC (m);
  if (detect_inf_r (&p, 1))
  {
    double maxr = cabs (MdcV0 (m, 0));
    for (i = 1; i < size; i++)
    {
      if (cabs (MdcV0 (m, i)) > maxr)
        maxr = cabs (MdcV0 (m, i));
    }
    *pnorm = maxr;
  }
  else
  {
    for (i = 0; i < size; i++)
      sum += errcheck (pow ((cabs (MdcV0 (m, i))), p), "pow");
    *pnorm = pow (sum, 1.0 / p);
  }
  return (pnorm);
}

/* **************************************************************
 * Singular value decomposition.
 * ************************************************************** */

void
mdc_Svd (MDC * M, MDC ** rsv, MDC ** lsv, MDR ** sigma, int flag)
{
  int k, lda, ldu, ldvt, lwork, m, n, info;
  int jobu, jobvt;
  MDC *A, *u, *vt, *work;
  MDR *rwork, *s;

  u = vt = 0;			/* Initialize */
  m = MNR (M);
  n = MNC (M);
  lda = m;
  k = MIN(m, n);
  lwork = 2 * MIN(m, n) + MAX (m, n);

  if (flag == 1)
  {
    jobu = (int) 'A';
    jobvt = (int) 'A';
    ldu = m;
    ldvt = n;
    u = mdc_Create (ldu, m);
    vt = mdc_Create (ldvt, n);
  }
  else if (flag == 2)
  {
    jobu = (int) 'S';
    jobvt = (int) 'S';
    ldu = m;
    ldvt = k;
    u = mdc_Create (ldu, k);
    vt = mdc_Create (ldvt, n);
  }
  else if (flag == 3)
  {
    jobu = (int) 'N';
    jobvt = (int) 'N';
    ldu = 1;
    ldvt = 1;
    u = mdc_Create (0, 0);
    vt = mdc_Create (0, 0);
  }

  A = mdc_Copy (M);
  rwork = mdr_Create (1, 5 * MAX (m, n));
  s = mdr_Create (1, k);
  work = mdc_Create (1, lwork);

  signal (SIGINT, intcatch);
  XGESVD (&jobu, &jobvt, &m, &n, MDCPTR (A), &lda, MDRPTR (s),
	  MDCPTR (u), &ldu, MDCPTR (vt), &ldvt,
	  MDCPTR (work), &lwork, MDRPTR (rwork), &info);
  signal (SIGINT, intcatch_wait);

  if ((int) info < 0)
    rerror ("bad argument to LAPACK ZGESVD");
  if ((int) info > 0)
    rerror ("svd algorithm failed to converge");

  /* Clean Up */
  mdc_Destroy (A);
  mdr_Destroy (rwork);
  mdc_Destroy (work);

  /* Set proper addresses */
  *rsv = vt;
  *lsv = u;
  *sigma = s;
}

/* **************************************************************
 * Least-Sqaures.
 * ************************************************************** */

MDC *
mdc_LS (MDC * a, MDC * b)
{
  int info, lda, ldb, lwork, m, n, nrhs, rank;
  double rcond;
  MDR *rwork, *s;
  MDC *A, *B, *Btmp, *work;

  m = MNR (a);
  n = MNC (a);
  nrhs = MNC (b);
  lda = m;
  ldb = MAX (m, n);
  rcond = -1.0;

  if (m >= n)
    lwork = 2 * n + MAX (nrhs, m);
  else
    lwork = 2 * m + MAX (nrhs, n);

  s = mdr_Create (1, MIN(m, n));
  work = mdc_Create (1, lwork);
  A = mdc_Copy (a);
  rwork = mdr_Create (1, MAX (5 * MIN(m, n) - 4, 1));

  /* Check to make sure B is O.K. */
  if (ldb > m)
  {
    int i, j;
    B = mdc_Create (MNR (b) + (ldb - m), MNC (b));
    for (i = 0; i < MNR (b); i++)
    {
      for (j = 0; j < MNC (b); j++)
      {
	Mdc0r (B, i, j) = Mdc0r (b, i, j);
	Mdc0i (B, i, j) = Mdc0i (b, i, j);
      }
    }
  }
  else
  {
    B = mdc_Copy (b);
  }

  signal (SIGINT, intcatch);
  XGELSS (&m, &n, &nrhs, MDCPTR (A), &lda, MDCPTR (B), &ldb, MDRPTR (s),
	  &rcond, &rank, MDCPTR (work), &lwork, MDRPTR (rwork), &info);
  signal (SIGINT, intcatch_wait);

  if ((int) info < 0)
  {
    fprintf (stderr, "ERROR: %ith argument to ZGELSS is bad\n", (int) -info);
    rerror ("illegal argument to LAPACK ZGELSS()");
  }
  else if ((int) info > 0)
    rerror ("SVD algorithm failed to converge");

  mdr_Destroy (s);
  mdc_Destroy (work);
  mdc_Destroy (A);
  mdr_Destroy (rwork);

  /* re-adjust B is necessary */
  if (m > n)
  {
    int i, j;
    Btmp = mdc_Create (n, MNC (B));
    for (i = 0; i < n; i++)
    {
      for (j = 0; j < MNC (B); j++)
      {
	Mdc0r (Btmp, i, j) = Mdc0r (B, i, j);
	Mdc0i (Btmp, i, j) = Mdc0i (B, i, j);
      }
    }
    mdc_Destroy (B);
    return (Btmp);
  }
  else
  {
    return (B);
  }
}

/* **************************************************************
 * Driver for Standard eigenvalue problem for complex dense matrices.
 * The two options (at this time) are general, and symmetric.
 * ************************************************************** */

static void mdc_EigS_SYM (MDC * M, void **val, void **vec);
static void mdc_EigS_GE (MDC * M, void **val,
                         void **vec, void **lvec, int *lflag);

Btree *
mdc_EigS (MDC * a, int * issym)
{
  Btree *bt;
  Ent *eval, *evec;
  void *val, *vec, *lvec;

  /* Create the list that will contain the results. */

  bt = btree_Create ();

  /*
   * Check for symmetry of [A]
   * Then use the correct routine.
   */
  if (*issym == -1) *issym = mdc_IsSymmetric (a);
  if (*issym == 1)
  {
    mdc_EigS_SYM (a, &val, &vec);

    /* Hook the results into the list. */
    eval = ent_Create ();
    ent_SetType (eval, MATRIX_DENSE_REAL);
    ent_data (eval) = val;
    install (bt, ("val"), eval);
    evec = ent_Create ();
    ent_SetType (evec, MATRIX_DENSE_COMPLEX);
    ent_data (evec) = vec;
    install (bt, ("vec"), evec);
  }
  else
  {
    mdc_EigS_GE (a, &val, &vec, &lvec, 0);

    /* Hook the results into the list. */
    eval = ent_Create ();
    ent_SetType (eval, MATRIX_DENSE_COMPLEX);
    ent_data (eval) = val;
    install (bt, ("val"), eval);
    evec = ent_Create ();
    ent_SetType (evec, MATRIX_DENSE_COMPLEX);
    ent_data (evec) = vec;
    install (bt, ("vec"), evec);
  }
//   ent_Clean( eval );
//   ent_Clean( evec );
  return (bt);
}

/* **************************************************************
 * Compute eigenvalues, and vectors for symmetric prob
 * [A]x = Lambda*x.
 * ************************************************************** */

static void
mdc_EigS_SYM (MDC * A, void **val, void **vec)
{
  int info, lda, lwork, m, n;
  int jobz, uplo;
  MDC *work;
  MDR *w, *rwork;

  w = 0;			/* Initialize */

  mdc_Detect_Inf (A);
  mdc_Detect_Nan (A);

  /* Some rudimentary checks */
  m = MNR (A);
  n = MNC (A);
  lda = m;

  jobz = (int) 'V';
  uplo = (int) 'L';

  if (m != n)
    rerror ("eig: input must be square");

  lwork = 2 * n - 1;
  w = mdr_Create (1, n);
  work = mdc_Create (1, lwork);
  rwork = mdr_Create (1, 3 * n - 2);

  XHEEV (&jobz, &uplo, &n, MDCPTR(A), &lda, MDRPTR(w),
          MDCPTR(work), &lwork, MDRPTR(rwork), &info);

  if ((int) info < 0)
    rerror ("bad argument to LAPACK DSYEV or ZHEEV");
  if ((int) info > 0)
    rerror ("Eigensolver failed to converge");

  *val = w;
  *vec = A;

  mdc_Destroy (work);
  mdr_Destroy (rwork);
  return;
}

/* **************************************************************
 * Compute eigenvalues, and vectors for non-symmetric prob
 * [A]x = Lambda*x.
 * lflag =  0 No left-eigenvectors
 * lflag != 0 Return left-eigenvectors
 * ************************************************************** */

static void
mdc_EigS_GE (MDC * A, void **val, void **vec, void **lvec, int *lflag)
{
  int info, lda, ldvl, ldvr, lwork, m, n;
  int jobvl, jobvr;
  MDC *w=0, *work=0, *vl=0, *vr=0;
  MDR *rwork=0;

  mdc_Detect_Inf (A);
  mdc_Detect_Nan (A);

  /* Some rudimentary checks */
  m = MNR (A);
  n = MNC (A);
  lda = m;
  ldvl = n;
  ldvr = n;

  if (m != n)
    rerror ("eig: input must be square");

  lwork = 2 * n;
  w = mdc_Create (1, n);
  vr = mdc_Create (ldvr, n);

  if (lflag == 0)
  {
    vl = mdc_Create (1, 1);
  }
  else
  {
    vl = mdc_Create (ldvr, n);
  }

  work = mdc_Create (1, lwork);
  rwork = mdr_Create (1, 2 * n);

  if (lflag == 0)
    jobvl = (int) 'N';
  else
    jobvl = (int) 'V';

  jobvr = (int) 'V';
  XGEEV (&jobvl, &jobvr, &n, MDCPTR(A), &lda, MDCPTR(w), MDCPTR(vl),
          &ldvl, MDCPTR(vr), &ldvr, MDCPTR(work), &lwork,
          MDRPTR(rwork), &info);

  if ((int) info < 0)
    rerror ("bad argument to LAPACK DGEEV or ZGEEV");
  if ((int) info > 0)
    rerror ("failed to compute all of the eigenvalues");

  *val = w;
  *vec = vr;
  if (lflag != 0)
  {
    *lvec = vl;
  }
  else
  {
    mdc_Destroy (vl);
  }
  mdc_Destroy (work);
  mdr_Destroy (rwork);
  return;
}

/* **************************************************************
 * Driver for Generalized eigenvalue problem for complex dense matrices.
 * The two options (at this time) are general, and symmetric.
 * ************************************************************** */

void mdc_EigG_SYM (MDC * A, MDC * B, void **val, void **vec);
void mdc_EigG_GE (MDC * M, MDC * B, void **val, void **vec);

Btree *
mdc_EigG (MDC * a, MDC * b, int *issympos)
{
  Btree *bt;
  Ent *eval, *evec;
  void *val, *vec;

  /* Create the list that will contain the results. */

  bt = btree_Create ();

  if (*issympos==1)
  {
    mdc_EigG_SYM (a, b, &val, &vec);
    // Hook the results into the list.
    eval = ent_Create ();
    ent_SetType (eval, MATRIX_DENSE_REAL);
    ent_data (eval) = val;
    install (bt, ("val"), eval);
    evec = ent_Create ();
    ent_SetType (evec, MATRIX_DENSE_COMPLEX);
    ent_data (evec) = vec;
    install (bt, ("vec"), evec);
  }
  else
  {
    mdc_EigG_GE (a, b, &val, &vec);

    /* Hook the results into the list. */
    eval = ent_Create ();
    ent_SetType (eval, MATRIX_DENSE_COMPLEX);
    ent_data (eval) = val;
    install (bt, ("val"), eval);

    evec = ent_Create ();
    ent_SetType (evec, MATRIX_DENSE_COMPLEX);
    ent_data (evec) = vec;
    install (bt, ("vec"), evec);
  }
//   ent_Clean( eval );
//   ent_Clean( evec );
  return (bt);
}

void
mdc_EigG_SYM (MDC * A, MDC * B, void **val, void **vec)
{
  int info, itype, lda, lwork, n, jobz, uplo;
  MDC *work;
  MDR *w, *rwork;

  work = 0;
  rwork = 0;

  mdc_Detect_Inf (A);
  mdc_Detect_Nan (A);
  mdc_Detect_Inf (B);
  mdc_Detect_Nan (B);

  /* Some rudimentary checks, MA, MB must be square */
  if (MNR (A) != MNC (A))
    rerror ("eig: A must be symmetric");
  if (MNR (B) != MNC (B))
    rerror ("eig: B must be symmetric");
  if (MNR (A) != MNR (B))
    rerror ("eig: A and B must be same size");

  itype = 1;
  n = MNR (A);
  lda = n;

  jobz = (int) 'V';
  uplo = (int) 'L';

  w = mdr_Create (1, n);

  lwork = MAX (1, 2 * n - 1);
  work = mdc_Create (1, lwork);
  rwork = mdr_Create (1, 3 * n - 2);

  XHEGV (&itype, &jobz, &uplo, &n, MDCPTR(A), &lda, MDCPTR(B), &lda,
          MDRPTR(w), MDCPTR(work), &lwork, MDRPTR(rwork), &info);

  if ((int) info < 0)
    rerror ("bad argument to LAPACK ZHEGV");
  if ((int) info > 0)
    rerror ("failure in eigensolver ZHEGV");

  *val = w;
  *vec = A;

  /* Clean Up */

  mdc_Destroy (work);
  mdr_Destroy (rwork);
}

void
mdc_EigG_GE (MDC * A, MDC * B, void **val, void **vecr)
{
  int i;
  int info, lda, ldvl, ldvr, lwork, n;
  int jobvl, jobvr;
  MDC *alpha, *beta;
  MDC *vl, *vr;
  MDC *work;
  MDR *rwork;
  MDC *eval;

  alpha = 0;
  vl = vr = work = beta = 0;
  rwork = 0;

  mdc_Detect_Inf (A);
  mdc_Detect_Nan (A);
  mdc_Detect_Inf (B);
  mdc_Detect_Nan (B);

  /* Some rudimentary checks, MA, MB must be square */
  if (MNR (A) != MNC (A))
    rerror ("eig: A must be symmetric");
  if (MNR (B) != MNC (B))
    rerror ("eig: B must be symmetric");
  if (MNR (A) != MNR (B))
    rerror ("eig: A and B must be same size");

  n = MNR (A);
  lda = ldvl = ldvr = n;

  jobvl = (int) 'N';
  jobvr = (int) 'V';

  lwork = MAX (1, 2 * n);
  work = mdc_Create (1, lwork);
  rwork = mdr_Create (1, 8 * n);
  alpha = mdc_Create (1, n);
  beta = mdc_Create (1, n);
  vl = mdc_Create (1, 1);
  vr = mdc_Create (ldvr, n);

  XGEGV (&jobvl, &jobvr, &n, MDCPTR(A), &lda, MDCPTR(B), &lda,
          MDCPTR(alpha), MDCPTR(beta), MDCPTR(vl), &ldvl,
          MDCPTR(vr), &ldvr, MDCPTR(work), &lwork, MDRPTR(rwork), &info);

  if ((int) info < 0)
    rerror ("bad argument to LAPACK ZGEGV");
  if ((int) info > 0)
    rerror ("failure in eigensolver ZGEGV");

  /*
   * Now compute eigenvalues
   */

  eval = mdc_Create (1, n);

  for (i = 1; i <= n; i++)
  {
    if ((Mdc1r (beta, 1, i) == 0.0) && (Mdc1i (beta, 1, i) == 0.0))
    {
      if ((Mdc1r (alpha, 1, i) == 0.0) && Mdc1i (alpha, 1, i) == 0.0)
      {
        Mdc1r (eval, 1, i) = create_nan ();
        Mdc1i (eval, 1, i) = create_nan ();
      }
      else
      {
        Mdc1r (eval, 1, i) = create_inf ();
        Mdc1i (eval, 1, i) = create_inf ();
      }
    }
    else
    {
      Mdc1(eval, 1, i) = Mdc1 (alpha, 1, i) / Mdc1 (beta, 1, i);
    }
  }

  *vecr = vr;
  *val = eval;

  /* Clean Up. */
  mdc_Destroy (work);
  mdc_Destroy (alpha);
  mdc_Destroy (beta);
  mdc_Destroy (vl);
}

/* **************************************************************
 * Matrix Factor.
 * ************************************************************** */

static void mdc_Factor_Ge (MDC * m, MDC ** lu, MDR ** pvt, double *cond);
static void mdc_Factor_Sym (MDC * m, MDC ** ldl, MDR ** pvt, double *cond);

Btree *
mdc_Factor (MDC * a, int flag)
{
  double rcond;
  Btree *bt;
  Ent *elu, *epvt, *ercond;
  MDC *lu;
  MDR *pvt;

  /* Create the list that will contain the results. */

  bt = btree_Create ();

  /*
   * Check for symmetry of [A]
   * Then use the correct routine.
   */

  if (flag == 0)
  {
    /*
     * Its up to us to check symmetry, and do
     * the right thing.
     */

    if (mdc_IsSymmetric (a))
    {
      mdc_Factor_Sym (a, &lu, &pvt, &rcond);

      /* Hook the results into the list. */
      elu = ent_Create ();
      ent_SetType (elu, MATRIX_DENSE_COMPLEX);
      ent_data (elu) = lu;
      install (bt, ("ldl"), elu);

      epvt = ent_Create ();
      ent_SetType (epvt, MATRIX_DENSE_REAL);
      ent_data (epvt) = pvt;
      install (bt, ("pvt"), epvt);

      ercond = ent_Create ();
      ent_SetType (ercond, DOUBLE);
      ent_double (ercond) = rcond;
      install (bt, ("rcond"), ercond);
    }
    else
    {
      mdc_Factor_Ge (a, &lu, &pvt, &rcond);

      /* Hook the results into the list. */
      elu = ent_Create ();
      ent_SetType (elu, MATRIX_DENSE_COMPLEX);
      ent_data (elu) = lu;
      install (bt, ("lu"), elu);

      epvt = ent_Create ();
      ent_SetType (epvt, MATRIX_DENSE_REAL);
      ent_data (epvt) = pvt;
      install (bt, ("pvt"), epvt);

      ercond = ent_Create ();
      ent_SetType (ercond, DOUBLE);
      ent_double (ercond) = rcond;
      install (bt, ("rcond"), ercond);
    }
  }
  else if (flag == 1)
  {
    /* User is forcing the general solution. */
    mdc_Factor_Ge (a, &lu, &pvt, &rcond);

    /* Hook the results into the list. */
    elu = ent_Create ();
    ent_SetType (elu, MATRIX_DENSE_COMPLEX);
    ent_data (elu) = lu;
    install (bt, ("lu"), elu);

    epvt = ent_Create ();
    ent_SetType (epvt, MATRIX_DENSE_REAL);
    ent_data (epvt) = pvt;
    install (bt, ("pvt"), epvt);

    ercond = ent_Create ();
    ent_SetType (ercond, DOUBLE);
    ent_double (ercond) = rcond;
    install (bt, ("rcond"), ercond);
  }
  else if (flag == 2)
  {
    /* User is forcing the symmetric solution. */
    /* User is forcing the general solution. */
    mdc_Factor_Sym (a, &lu, &pvt, &rcond);

    /* Hook the results into the list. */
    elu = ent_Create ();
    ent_SetType (elu, MATRIX_DENSE_COMPLEX);
    ent_data (elu) = lu;
    install (bt, ("ldl"), elu);

    epvt = ent_Create ();
    ent_SetType (epvt, MATRIX_DENSE_REAL);
    ent_data (epvt) = pvt;
    install (bt, ("pvt"), epvt);

    ercond = ent_Create ();
    ent_SetType (ercond, DOUBLE);
    ent_double (ercond) = rcond;
    install (bt, ("rcond"), ercond);
  }

  return (bt);
}

static void
mdc_Factor_Sym (MDC * m, MDC ** ldl, MDR ** pvt, double *cond)
{
  double *anorm, rcond;
  int i;
  int info, mm, n, lda, *ipiv;
  int lwork, uplo;
  MDC *tldl, *work;
  MDR *tpvt;

  lda = mm = MNR (m);
  n = MNC (m);
  uplo = (int) 'L';

  /*
   *Try and pick a good NB, without ILAENV.
   */

  if (n < 100)
    lwork = n;
  else
    lwork = 64 * n;

  if (mm != n)
    rerror ("factor: coefficient matrix must be square");

  ipiv = (int *) GC_malloc_atomic_ignore_off_page (sizeof (int) * mm);
  tpvt = mdr_Create (1, mm);
  tldl = mdc_Copy (m);
  work = mdc_Create (1, lwork);

  signal (SIGINT, intcatch);
  XHETRF (&uplo, &n, MDCPTR (tldl), &lda, ipiv, MDCPTR (work), &lwork, &info);
  signal (SIGINT, intcatch_wait);

  if ((int) info < 0)
    rerror ("bad argument to LAPACK ZHETRF");
  if ((int) info > 0)
    rerror ("matrix is singular");

  /*
   * Check input for ill-conditioning.
   */

  anorm = mdc_Norm (m, "1");
  mdc_Destroy (work);
  work = mdc_Create (1, 2 * n);

  signal (SIGINT, intcatch);
  XHECON (&uplo, &n, MDCPTR (tldl), &lda, ipiv, anorm, &rcond,
	  MDCPTR (work), &info);
  signal (SIGINT, intcatch_wait);
  GC_FREE (anorm);
  if (rcond <= DBL_EPSILON)
    warning_1 ("WARNING, ill-conditioned input to factor()");

  /* Fill pvt */
  for (i = 0; i < mm; i++)
  {
    Mdr0 (tpvt, 0, i) = ipiv[i];
  }

  *ldl = tldl;
  *pvt = tpvt;
  *cond = rcond;

  /* Clean up */
  GC_FREE (ipiv);
  mdc_Destroy (work);
}

static void
mdc_Factor_Ge (MDC * m, MDC ** lu, MDR ** pvt, double *cond)
{
  double *anorm, rcond, *rwork;
  int i;
  int info, mm, n, lda, *ipiv, lwork, norm;
  MDC *tlu, *work;
  MDR *tpvt;

  lwork = lda = mm = MNR (m);
  n = MNC (m);
  norm = (int) '1';

  if (mm != n)
    rerror ("factor: coefficient matrix must be square");

  ipiv = (int *) GC_malloc_atomic_ignore_off_page (sizeof (int) * mm);
  tpvt = mdr_Create (1, mm);
  tlu = mdc_Copy (m);
  rwork = (double *) GC_malloc_atomic_ignore_off_page (sizeof (double) * 2 * n);
  work = mdc_Create (1, 2 * lwork);

  signal (SIGINT, intcatch);
  XGETRF (&mm, &n, MDCPTR (tlu), &lda, ipiv, &info);
  signal (SIGINT, intcatch_wait);

  if ((int) info < 0)
    rerror ("bad argument to LAPACK ZGETRF");
  if ((int) info > 0)
    rerror ("matrix is singular");

  /*
   * Check input for ill-conditioning.
   */

  anorm = mdc_Norm (m, "1");
  signal (SIGINT, intcatch);
  XGECON (&norm, &n, MDCPTR (tlu), &lda, anorm, &rcond, MDCPTR (work),
	  rwork, &info);
  signal (SIGINT, intcatch_wait);
  GC_FREE (anorm);
  if (rcond <= DBL_EPSILON)
    warning_1 ("WARNING, ill-conditioned input to factor()");

  /* Fill pvt */
  for (i = 0; i < mm; i++)
  {
    Mdr0 (tpvt, 0, i) = ipiv[i];
  }

  *lu = tlu;
  *pvt = tpvt;
  *cond = rcond;

  /* Clean up */
  GC_FREE (ipiv);
  GC_FREE (rwork);
  mdc_Destroy (work);
}

/* **************************************************************
 * Matrix Backsub function.
 * ************************************************************** */

MDC *
mdc_Backsub_Sym (Btree * bl, MDC * b)
{
  int i;
  int *ipiv, lda, ldb, mm, n, nrhs, info;
  int uplo;
  Ent *etmp;
  ListNode *tmp;
  MDC *B, *lu;
  MDR *pvt;
  lu = 0;
  pvt = 0;

  /* Get the necessary elements off of the list. */

  if ((tmp = btree_FindNode (bl, "ldl")))
  {
    etmp = var_ent (tmp);
    if (ent_type (etmp) != MATRIX_DENSE_COMPLEX)
      rerror ("backsub: terrible input list error");
    lu = ent_data (etmp);
  }
  else
  {
    rerror ("backsub: terrible input list error");
  }

  if ((tmp = btree_FindNode (bl, "pvt")))
  {
    etmp = var_ent (tmp);
    if (ent_type (etmp) != MATRIX_DENSE_REAL)
      rerror ("backsub: terrible input list error");
    pvt = ent_data (etmp);
  }
  else
  {
    rerror ("backsub: terrible input list error");
  }

  lda = mm = MNR (lu);
  n = MNC (lu);
  nrhs = MNC (b);
  ldb = MNR (b);
  uplo = (int) 'L';

  if (ldb != lda)
    rerror ("backsub: b must have same number of rows as LHS");

  ipiv = (int *) GC_malloc_atomic_ignore_off_page (sizeof (int) * mm);
  B = mdc_Copy (b);

  /* Fill ipiv */
  for (i = 0; i < mm; i++)
  {
    ipiv[i] = (int) Mdr0 (pvt, 0, i);
  }

  signal (SIGINT, intcatch);
  XHETRS (&uplo, &n, &nrhs, MDCPTR (lu), &lda, ipiv, MDCPTR (B), &ldb, &info);
  signal (SIGINT, intcatch_wait);
  if ((int) info < 0)
    rerror ("Bad argument(s) to LAPACK ZHETRF");

  GC_FREE (ipiv);
  return (B);
}

MDC *
mdc_Backsub_Ge (Btree * bl, MDC * b)
{
  int i;
  int *ipiv, lda, ldb, mm, n, nrhs, info;
  int trans;
  Ent *etmp;
  ListNode *tmp;
  MDC *B, *lu;
  MDR *pvt;
  lu = 0;
  pvt = 0;

  /* Get the necessary elements off of the list. */

  if ((tmp = btree_FindNode (bl, "lu")))
  {
    etmp = var_ent (tmp);
    if (ent_type (etmp) != MATRIX_DENSE_COMPLEX)
      rerror ("backsub: terrible input list error");
    lu = ent_data (etmp);
  }
  else
  {
    rerror ("backsub: terrible input list error");
  }

  if ((tmp = btree_FindNode (bl, "pvt")))
  {
    etmp = var_ent (tmp);
    if (ent_type (etmp) != MATRIX_DENSE_REAL)
      rerror ("backsub: terrible input list error");
    pvt = ent_data (etmp);
  }
  else
  {
    rerror ("backsub: terrible input list error");
  }

  lda = mm = MNR (lu);
  n = MNC (lu);
  nrhs = MNC (b);
  ldb = MNR (b);
  trans = (int) 'N';

  if (ldb != lda)
    rerror ("backsub: b must have same number of rows as LHS");

  ipiv = (int *) GC_malloc_atomic_ignore_off_page (sizeof (int) * mm);
  B = mdc_Copy (b);

  /* Fill ipiv */
  for (i = 0; i < mm; i++)
  {
    ipiv[i] = (int) Mdr0 (pvt, 0, i);
  }

  signal (SIGINT, intcatch);
  XGETRS (&trans, &n, &nrhs, MDCPTR (lu), &lda, ipiv, MDCPTR (B), &ldb, &info);
  signal (SIGINT, intcatch_wait);
  if ((int) info < 0)
    rerror ("Bad argument(s) to LAPACK ZGETRF");

  GC_FREE (ipiv);
  return (B);
}

/* **************************************************************
 * Matrix SVD function.
 * ************************************************************** */

Btree *
mdc_Svd_BF (MDC * M, int flag)
{
  Btree *bt;
  Ent *esigma, *eu, *evt;
  MDC *rsv, *lsv;
  MDR *sigma;

  mdc_Svd (M, &rsv, &lsv, &sigma, flag);

  /* Hook the results into the list. */

  bt = btree_Create ();

  eu = ent_Create ();
  ent_SetType (eu, MATRIX_DENSE_COMPLEX);
  ent_data (eu) = lsv;
  install (bt, ("u"), eu);

  esigma = ent_Create ();
  ent_SetType (esigma, MATRIX_DENSE_REAL);
  ent_data (esigma) = sigma;
  install (bt, ("sigma"), esigma);

  evt = ent_Create ();
  ent_SetType (evt, MATRIX_DENSE_COMPLEX);
  ent_data (evt) = rsv;
  install (bt, ("vt"), evt);

  return (bt);
}

/* **************************************************************
 * Matrix Cholesky decomposition.
 * ************************************************************** */

MDC *
mdc_Chol (MDC * m)
{
  int i, j, info, lda, n, uplo;
  MDC *A;

  mdc_Detect_Inf (m);
  mdc_Detect_Nan (m);

  lda = n = MNC (m);
  uplo = (int) 'U';

  /* [m] gets overwritten, so copy it */
  A = mdc_Copy (m);

  /* Call LAPACK routine */

  signal (SIGINT, intcatch);
  XPOTRF (&uplo, &n, MDCPTR (A), &lda, &info);
  signal (SIGINT, intcatch_wait);

  if ((int) info < 0)
    rerror ("bad argument to LAPACK [DZ]OPTRF");
  if ((int) info > 0)
    rerror ("chol: input not positive definite");

  /* We must make sure and zero out the lower left triangle of [A] */

  for (i = 1; i < n; i++)
  {
    for (j = 0; j < i; j++)
    {
      Mdc0r (A, i, j) = 0.0;
      Mdc0i (A, i, j) = 0.0;
    }
  }

  return (A);
}

/* **************************************************************
 * Matrix Hessenberg reduction.
 * ************************************************************** */

Btree *
mdc_Hess (MDC * M)
{
  int i, j;
  int ilo, ihi, lda, m, n, info, lwork;
  MDC *A, *tau, *work;

  Btree *bt;
  Ent *ep, *eh;
  MDC *h;

  m = MNR (M);
  n = MNC (M);
  lda = m;
  lwork = n;
  ilo = 1;
  ihi = n;

  /* Set up work arrays */
  A = mdc_Copy (M);
  tau = mdc_Create (1, n);
  work = mdc_Create (1, lwork);

  signal (SIGINT, intcatch);
  XGEHRD (&n, &ilo, &ihi, MDCPTR (A), &lda, MDCPTR (tau),
	  MDCPTR (work), &lwork, &info);
  signal (SIGINT, intcatch_wait);

  if ((int) info < 0)
    rerror ("bad argument to LAPACK ZGEHRD");

  /* Save [h] */
  h = mdc_Copy (A);
  for (i = 2; i < m; i++)
  {
    for (j = 0; (j < i - 1) && (j < n); j++)
    {
      Mdc0r (h, i, j) = 0.0;
      Mdc0i (h, i, j) = 0.0;
    }
  }

  signal (SIGINT, intcatch);
  XUNGHR (&n, &ilo, &ihi, MDCPTR (A), &lda, MDCPTR (tau),
	  MDCPTR (work), &lwork, &info);
  signal (SIGINT, intcatch_wait);

  if ((int) info < 0)
    rerror ("bad argument to LAPACK ZORGHR");

  /* clean up */
  mdc_Destroy (tau);
  mdc_Destroy (work);

  bt = btree_Create ();

  ep = ent_Create ();
  ent_SetType (ep, MATRIX_DENSE_COMPLEX);
  ent_data (ep) = A;
  install (bt, ("p"), ep);

  eh = ent_Create ();
  ent_SetType (eh, MATRIX_DENSE_COMPLEX);
  ent_data (eh) = h;
  install (bt, ("h"), eh);

  return (bt);
}

/* **************************************************************
 * Balance a matrix.
 * ************************************************************** */

Btree *
mdc_Balance (MDC * m)
{
  double *scale;
  int i;
  Btree *bt;
  Ent *et, *eab;
  int info, nm, n, low, igh, job;
  MDC *A;
  MDR *T;

  mdc_Detect_Inf (m);
  mdc_Detect_Nan (m);

  if (MNR (m) != MNC (m))
    rerror ("input to balance must be square matrix");

  nm = (int) MNR (m);
  n = nm;
  job = (int) 'S';

  A = mdc_Copy (m);
  scale = (double *) GC_malloc_atomic_ignore_off_page (n * sizeof (double));

  /* LAPACK balance */
  signal (SIGINT, intcatch);
  XGEBAL (&job, &n, MDCPTR (A), &n, &low, &igh, scale, &info);
  signal (SIGINT, intcatch_wait);

  if (info)
    rerror ("error in argument to balance()");

  T = mdr_Create ((int) nm, (int) nm);
  mdr_Zero (T);
  for (i = 0; i < n; i++)
  {
    Mdr0 (T, i, i) = scale[i];
  }

  GC_FREE (scale);
  bt = btree_Create ();

  et = ent_Create ();
  ent_SetType (et, MATRIX_DENSE_REAL);
  ent_data (et) = T;
  install (bt, ("t"), et);

  eab = ent_Create ();
  ent_SetType (eab, MATRIX_DENSE_COMPLEX);
  ent_data (eab) = A;
  install (bt, ("ab"), eab);

  return (bt);
}

/* **************************************************************
 * QR decomposition.
 * ************************************************************** */

Btree *
mdc_QR (MDC * M)
{
  int i, j;
  Complex *tau, *work;
  int info, k, lwork, m, n;
  MDC *A, *R, *q;

  Btree *bt;
  Ent *eq, *er;

  mdc_Detect_Inf (M);
  mdc_Detect_Nan (M);

  /* Set dimensional paramaters */
  m = MNR (M);
  n = MNC (M);
  k = MIN(m, n);
  lwork = MAX (m, n);

  /*
   *  over / over determined q: M x M
   *                         r: M x N
   *                         p: N x N
   */

  /* Create work arrays */
  A = mdc_Copy (M);
  tau = (Complex *) GC_malloc_atomic_ignore_off_page (sizeof (Complex) * k);
  work =
    (Complex *) GC_malloc_atomic_ignore_off_page (sizeof (Complex) * lwork);

  signal (SIGINT, intcatch);
  XGEQRF (&m, &n, MDCPTR (A), &m, tau, work, &lwork, &info);
  signal (SIGINT, intcatch_wait);

  if ((int) info < 0)
    rerror ("bad argument to LAPACK ZGEQRF");

  /* Extract [R] */
  R = mdc_Create (m, n);
  mdc_Zero (R);
  for (j = 1; j <= MNC (R); j++)
  {
    for (i = 1; (i <= j) && (i <= MNR (R)); i++)
    {
      Mdc1r (R, i, j) = Mdc1r (A, i, j);
      Mdc1i (R, i, j) = Mdc1i (A, i, j);
    }
  }

  /* Now re-assemble [Q] */
  signal (SIGINT, intcatch);
  if ((m - n) > 0)
  {
    /* matrix_AppendColC (A, m - n); */
    mdc_Extend (A, MNR (A), MNC (A) + (m - n));
  }

  XUNGQR (&m, &m, &k, MDCPTR (A), &m, tau, work, &lwork, &info);
  signal (SIGINT, intcatch_wait);

  if ((int) info < 0)
    rerror ("bad argument to LAPACK ZUNGQR");

  if ((m - n) < 0)
  {
    /* Get the first M columns of A. */

    q = mdc_Create (MNR (A), m);
    for (i = 0; i < MNR (A); i++)
    {
      for (j = 0; j < m; j++)
      {
	Mdc0r (q, i, j) = Mdc0r (A, i, j);
	Mdc0i (q, i, j) = Mdc0i (A, i, j);
      }
    }
    mdc_Destroy (A);
  }
  else
  {
    q = A;
  }

  /* Clean - Up */
  GC_FREE (tau);
  GC_FREE (work);

  bt = btree_Create ();

  eq = ent_Create ();
  ent_SetType (eq, MATRIX_DENSE_COMPLEX);
  ent_data (eq) = q;
  install (bt, ("q"), eq);

  er = ent_Create ();
  ent_SetType (er, MATRIX_DENSE_COMPLEX);
  ent_data (er) = R;
  install (bt, ("r"), er);

  return (bt);
}

Btree *
mdc_QRP (MDC * M)
{
  int i, j, *jpvt;
  Complex *tau, *work;
  int info, k, lwork, m, n;
  MDC *A, *R, *q;
  MDR *rwork, *P;

  Btree *bt;
  Ent *eq, *er, *ep;

  mdc_Detect_Inf (M);
  mdc_Detect_Nan (M);

  /* Set dimensional paramaters */
  m = MNR (M);
  n = MNC (M);
  k = MIN(m, n);
  lwork = MAX (m, n);

  /*
   *  over / under determined q: M x M
   *                          r: M x N
   *                          p: N x N
   */

  /* Create work arrays */
  A = mdc_Copy (M);
  tau = (Complex *) GC_malloc_atomic_ignore_off_page (sizeof (Complex) * k);
  work =
    (Complex *) GC_malloc_atomic_ignore_off_page (sizeof (Complex) * lwork);
  jpvt = (int *) GC_malloc_atomic_ignore_off_page (sizeof (int) * n);

  for (i = 0; i < n; i++)
  {
    jpvt[i] = 0;		/* It is important to zero-out jpvt */
  }

  rwork = mdr_Create (2 * n, 1);

  signal (SIGINT, intcatch);
  XGEQPF (&m, &n, MDCPTR (A), &m, jpvt, tau, work, MDRPTR (rwork), &info);
  signal (SIGINT, intcatch_wait);

  if ((int) info < 0)
    rerror ("bad argument to LAPACK ZGEQPF");

  /* Extract [R] */
  R = mdc_Create (m, n);
  mdc_Zero (R);
  for (j = 1; j <= MNC (R); j++)
  {
    for (i = 1; (i <= j) && (i <= MNR (R)); i++)
    {
      Mdc1r (R, i, j) = Mdc1r (A, i, j);
      Mdc1i (R, i, j) = Mdc1i (A, i, j);
    }
  }

  /* Form [P] */
  P = mdr_Create (n, n);
  mdr_Zero (P);
  for (i = 1; i <= n; i++)
  {
    if (jpvt[i - 1] != 0)
      Mdr1 (P, jpvt[i - 1], i) = 1.0;
    else
      Mdr1 (P, i, i) = 1.0;
  }

  /* Now re-assemble [Q] */
  signal (SIGINT, intcatch);
  if ((m - n) > 0)
  {
    /* matrix_AppendColC (A, m - n); */
    mdc_Extend (A, MNR (A), MNC (A) + (m - n));
  }

  XUNGQR (&m, &m, &k, MDCPTR (A), &m, tau, work, &lwork, &info);
  signal (SIGINT, intcatch_wait);

  if ((int) info < 0)
    rerror ("bad argument to LAPACK ZUNGQR");

  if ((m - n) < 0)
  {
    /* Get the first M columns of A. */

    q = mdc_Create (MNR (A), m);
    for (i = 0; i < MNR (A); i++)
    {
      for (j = 0; j < m; j++)
      {
	Mdc0r (q, i, j) = Mdc0r (A, i, j);
	Mdc0i (q, i, j) = Mdc0i (A, i, j);
      }
    }
    mdc_Destroy (A);
  }
  else
  {
    q = A;
  }

  /* Clean - Up */
  GC_FREE (tau);
  GC_FREE (work);
  GC_FREE (jpvt);
  mdr_Destroy (rwork);

  bt = btree_Create ();

  eq = ent_Create ();
  ent_SetType (eq, MATRIX_DENSE_COMPLEX);
  ent_data (eq) = q;
  install (bt, ("q"), eq);

  er = ent_Create ();
  ent_SetType (er, MATRIX_DENSE_COMPLEX);
  ent_data (er) = R;
  install (bt, ("r"), er);

  ep = ent_Create ();
  ent_SetType (ep, MATRIX_DENSE_REAL);
  ent_data (ep) = P;
  install (bt, ("p"), ep);

  return (bt);
}

/* **************************************************************
 * Schur Decomposition
 * ************************************************************** */

static void (*rselect) ();	/* Dummy selction pointer */

Btree *
mdc_Schur (MDC * m)
{
  int *bwork;
  Btree *bt;
  Ent *et, *ez;
  int info, lda, lwork, n, sdim, jobvs, sort;
  MDC *a, *vs, *w, *work;
  MDR *rwork;

  /* m must be n-by-n */
  if (MNR (m) != MNC (m))
    rerror ("schur: argument ust be square (N-by-N)");

  n = (int) MNR (m);
  lda = n;
  sdim = 0;
  lwork = MAX (1, 3 * n);
  bwork = (int *) GC_MALLOC (2 * sizeof (int));	/* not used */
  jobvs = (int) 'V';
  sort = (int) 'N';

  a = mdc_Copy (m);
  w = mdc_Create (1, n);
  rwork = mdr_Create (1, n);
  vs = mdc_Create (n, n);
  work = mdc_Create (1, lwork);

  signal (SIGINT, intcatch);
  XGEES (&jobvs, &sort, rselect, &n, MDCPTR (a), &lda, &sdim, MDCPTR (w),
	 MDCPTR (vs), &n, MDCPTR (work), &lwork, MDRPTR (rwork), bwork, &info);
  signal (SIGINT, intcatch_wait);

  if ((int) info < 0)
  {
    rerror ("illegal argument to ZGEES");
  }
  else if ((int) info > 0)
  {
    if ((int) info <= n)
    {
      fprintf (stderr, "schur (ZGEES): Failed to compute all eigenvalues\n");
      fprintf (stderr, "               %i:%i converged eigenvalues\n",
	       (int) info, (int) n);
    }
    else if ((int) info == n + 1)
    {
      fprintf (stderr, "schur (ZGEES): Eigenvalues could not be re-ordered\n");
      fprintf (stderr, "               This problem is very ill-conditioned\n");
      rerror ("schur: cannot continue");
    }
    else if ((int) info == n + 2)
    {
      fprintf (stderr, "schur (ZGEES): Problems re-ordering\n");
      rerror ("schur: cannot continue");
    }
  }

  mdc_Destroy (w);
  mdr_Destroy (rwork);
  mdc_Destroy (work);
  GC_FREE (bwork);

  bt = btree_Create ();

  et = ent_Create ();
  ent_SetType (et, MATRIX_DENSE_COMPLEX);
  ent_data (et) = a;
  install (bt, ("t"), et);

  ez = ent_Create ();
  ent_SetType (ez, MATRIX_DENSE_COMPLEX);
  ent_data (ez) = vs;
  install (bt, ("z"), ez);

  return (bt);
}

/* **************************************************************
 * Sylvester Matrix Equation.
 * ************************************************************** */

MDC *
mdc_Sylv (MDC * a, MDC * b, MDC * c)
{
  double scale;
  int info, isgn, m, n, lda, ldb, ldc;
  int trana, tranb;
  MDC *C = 0;

  if (b == 0)
  {
    // Lyapunov equation

    isgn = 1;			/* A*X + X*A' = -C */
    scale = 1.0;

    // Set & Check dimensions
    m = lda = MDC_nrow (a);
    n = ldb = m;
    ldc = m;
    trana = (int) 'N';
    tranb = (int) 'C';	// was 'T' (E. Plischke)

    if (m != MDC_ncol (a))
      rerror ("sylv: A.nr must equal A.nc");
    if (m != MDC_nrow (c))
      rerror ("sylv: A.nr must equal C.nr");

    // Ensure that all inputs are complex

    C = mdc_Negate (c);

    signal (SIGINT, intcatch);
    XTRSYL (&trana, &tranb, &isgn, &m, &n, MDCPTR (a), &lda,
	    MDCPTR (a), &ldb, MDCPTR (C), &ldc, &scale, &info);
    signal (SIGINT, intcatch_wait);

    if ((int) info < 0)
    {
      fprintf (stderr, "sylv: %ith argument to ZTRSYL illegal\n", (int) info);
      rerror ("sylv: cannot continue");
    }
    else if ((int) info == 1)
    {
      fprintf (stderr,
	       "sylv: A and B have very close eigenvalues; perturbed values\n");
      fprintf (stderr,
	       "      were used to solve the equation (but the matrices \n");
      fprintf (stderr, "      A and B are unchanged\n");
    }
    return (C);
  }
  else
  {

    // Sylvester equation

    isgn = 1;			/* A*X + X*B = C */
    scale = 1.0;

    // Set & Check dimensions
    m = lda = MDC_nrow (a);
    n = ldb = MDC_nrow (b);
    ldc = m;
    trana = (int) 'N';
    tranb = (int) 'N';

    if (m != MDC_ncol (a))
      rerror ("sylv: A.nr must equal A.nc");
    if (n != MDC_ncol (b))
      rerror ("sylv: B.nr must equal B.nc");
    if (m != MDC_nrow (c))
      rerror ("sylv: A.nr must equal C.nr");
    if (n != MDC_ncol (c))
      rerror ("sylv: B.nr must equal C.nc");

    // Ensure that all inputs are complex

    signal (SIGINT, intcatch);
    XTRSYL (&trana, &tranb, &isgn, &m, &n, MDCPTR (a), &lda,
	    MDCPTR (b), &ldb, MDCPTR (C), &ldc, &scale, &info);
    signal (SIGINT, intcatch_wait);

    if ((int) info < 0)
    {
      fprintf (stderr, "sylv: %ith argument to ZTRSYL illegal\n", (int) info);
      rerror ("sylv: cannot continue");
    }
    else if ((int) info == 1)
    {
      fprintf (stderr,
	       "sylv: A and B have very close eigenvalues; perturbed values\n");
      fprintf (stderr,
	       "      were used to solve the equation (but the matrices \n");
      fprintf (stderr, "      A and B are unchanged\n");
    }
    return (C);
  }
}

/* **************************************************************
 * Matrix Determinant
 * ************************************************************** */

MDC *
mdc_Det (MDC * m)
{
  int i;
  int info, mm, n, lda, *ipiv;
  double d1, d2, ten;
  Complex det[2], z1, z2, rz;
  MDC *lu, *mrz;

  lda = mm = MNR (m);
  n = MNC (m);

  if (mm != n)
    rerror ("det: matrix argument must be square");

  ipiv = (int *) GC_malloc_atomic_ignore_off_page (sizeof (int) * mm);
  lu = mdc_Copy (m);

  signal (SIGINT, intcatch);
  XGETRF (&mm, &n, MDCPTR (lu), &lda, ipiv, &info);
  signal (SIGINT, intcatch_wait);

  if ((int) info < 0)
    rerror ("det: bad argument to LAPACK ZGETRF");
  if ((int) info > 0)
  {
    warning_1 ("det: matrix is singular");
    mrz = mdc_CreateScalar (0.0, 0.0);
    return (mrz);
  }

  /*
   * Now compute the determinant using the algorithm contained
   * in the LINPACK subroutines DGEDI, ZGEDI. I'm doing it this
   * way because the LAPACK does not provide a direct method for
   * computing the determinant?
   */

  det[0] = 1.0;
  det[1] = 0.0;
  ten = 10.;
  for (i = 1; i <= n; ++i)
  {
    if (ipiv[i - 1] != i)
    {
      z1 = -det[0];
      det[0] = z1;
    }
    z1 = Mdc1(lu, i, i) * det[0];
    det[0] = z1;
    /*        ...exit */
    RE(z1) =  IM(det[0]) ;
    IM(z1) = -RE(det[0]);
    if ((d1 = RE(det[0]), ABS(d1)) + (d2 = RE(z1), ABS(d2)) == 0.)
    {
      goto L60;
    }
  L10:
    RE(z1) =  IM(det[0]);
    IM(z1) = -RE(det[0]);
    if ((d1 = RE(det[0]), ABS(d1)) + (d2 = RE(z1), ABS(d2)) >= 1.)
    {
      goto L20;
    }
    z2 = ten;;
    z1 = z2 * det[0];
    det[0] = z1;
    z1 = det[1] - 1.;
    det[1] = z1;
    goto L10;
  L20:
  L30:
    RE(z1) =  IM(det[0]);
    IM(z1) = -RE(det[0]);
    if ((d1 = RE(det[0]), ABS(d1)) + (d2 = RE(z1), ABS(d2)) < ten)
    {
      goto L40;
    }
    z2 = ten;
    z1 = det[0] / z2;
    det[0] = z1;
    z1 = det[1] + 1.;
    det[1] = z1;
    goto L30;
  L40:
    ;
  }
L60:

  GC_FREE (ipiv);
  mdc_Destroy (lu);

  z1  = cpow (10.0, det[1]);
  rz  = det[0]*z1;
  mrz = mdc_CreateScalar (RE(rz), IM(rz));

  return (mrz);
}
