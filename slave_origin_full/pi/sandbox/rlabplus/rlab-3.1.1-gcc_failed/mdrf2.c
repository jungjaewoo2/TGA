/*
 * mdrf2.c
 * Matrix-Dense-Real Fortran Interfaces
 * Contains most of the computational interfaces to LAPACK.
 */

/*  This file is a part of RLaB ("Our"-LaB)
   Copyright (C) 1995  Ian R. Searle, (C) 2005 M. Kostrun

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
#include "mdrf1.h"
#include "mdc.h"
#include "util.h"
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

double *mdr_Norm (MDR * m, char *type);
void mdr_Svd (MDR * M, MDR ** rsv, MDR ** lsv, MDR ** sigma, int flag);

/* **************************************************************
 * Solve the system of equations [a]{x} = {b}
 * ************************************************************** */

MDR *mdr_SolveEq_GE (MDR * a, MDR * b);
MDR *mdr_SolveEq_SYM (MDR * a, MDR * b);
MDR *mdr_SolveEq (MDR * a, MDR * b);

/*
 * Cover package for builtin function solve().
 */

MDR *
mdr_Solve (MDR * A, MDR * B, char *type)
{
  MDR *m = 0, *a, *b;

  // no integer optimization
  if(A->type == RLAB_TYPE_INT32) a = mdr_Float_BF (A);
  else a = A;
  if(B->type == RLAB_TYPE_INT32) b = mdr_Float_BF (B);
  else b = B;

  if (type != 0)
  {
    if (!strncmp ("S", type, 1))
    {
      m = mdr_SolveEq_SYM (a, b);
    }
    else if (!strncmp ("s", type, 1))
    {
      m = mdr_SolveEq_SYM (a, b);
    }
    else if (!strncmp ("G", type, 1))
    {
      m = mdr_SolveEq_GE (a, b);
    }
    else if (!strncmp ("g", type, 1))
    {
      m = mdr_SolveEq_GE (a, b);
    }
  }
  else
  {
    m = mdr_SolveEq (a, b);
  }

  // no integer optimization
  if(A->type == RLAB_TYPE_INT32) mdr_Destroy (a);
  if(B->type == RLAB_TYPE_INT32) mdr_Destroy (b);

  return (m);
}

MDR *
mdr_SolveEq (MDR * A, MDR * B)
{
  MDR *m, *a, *b;

  // no integer optimization
  if(A->type == RLAB_TYPE_INT32) a = mdr_Float_BF (A);
  else a = A;
  if(B->type == RLAB_TYPE_INT32) b = mdr_Float_BF (B);
  else b = B;

  if (mdr_IsSymmetric (a))
  {
    m = mdr_SolveEq_SYM (a, b);
  }
  else
  {
    m = mdr_SolveEq_GE (a, b);
  }

  // no integer optimization
  if(A->type == RLAB_TYPE_INT32) mdr_Destroy (a);
  if(B->type == RLAB_TYPE_INT32) mdr_Destroy (b);

  return (m);
}

/* **************************************************************
 * Solve the general case...
 * ************************************************************** */

MDR *
mdr_SolveEq_GE (MDR * a, MDR * b)
{
  F_INT info, lda, ldb, mm, n, nrhs, *ipiv;
  F_INT trans, lwork, norm, *iwork;
  double *anorm, rcond;
  MDR *A, *B, *work;

  mdr_Detect_Inf (a);
  mdr_Detect_Nan (a);

  mdr_Detect_Inf (b);
  mdr_Detect_Nan (b);

  lda = mm = MNR (a);
  n = MNC (a);
  lwork = mm;
  trans = (F_INT) 'N';
  norm = (F_INT) '1';

  if (mm != n)
    rerror ("solve: coefficient matrix must be square");

  if (mm != MNR (b))
    rerror ("solve: RHS row dim. must match LHS row dim.");

  ipiv = (F_INT *) GC_malloc_atomic_ignore_off_page (sizeof (F_INT) * mm);
  A = mdr_Float_BF (a);
  B = mdr_Float_BF (b);
  work = mdr_Create (1, 4 * lwork);

  signal (SIGINT, intcatch);
  RGETRF (&mm, &n, MDRPTR (A), &lda, ipiv, &info);
  signal (SIGINT, intcatch_wait);

  if ((int) info < 0)
    rerror ("solve: bad argument to LAPACK DGETRF");
  if ((int) info > 0)
    rerror ("solve: matrix is singular");

  /*
   * Check input for ill-conditioning.
   */

  iwork = (F_INT *) GC_malloc_atomic_ignore_off_page (sizeof (F_INT) * n);

  anorm = mdr_Norm (a, "1");
  signal (SIGINT, intcatch);
  RGECON (&norm, &n, MDRPTR (A), &lda, anorm, &rcond,
          MDRPTR (work), iwork, &info);
  signal (SIGINT, intcatch_wait);

  if (rcond <= DBL_EPSILON)
  {
    fprintf (stderr, "RCOND = %e\n", rcond);
    warning_1 ("WARNING, ill-conditioned input");
  }

  GC_FREE (anorm);
  nrhs = MNC (b);
  ldb = MNR (b);

  signal (SIGINT, intcatch);
  RGETRS (&trans, &n, &nrhs, MDRPTR (A), &lda, ipiv, MDRPTR (B), &ldb, &info);
  signal (SIGINT, intcatch_wait);

  if ((int) info < 0)
    rerror ("solve: Bad argument(s) to LAPACK DGETRS");

  GC_FREE (ipiv);
  GC_FREE (A);
  mdr_Destroy (work);
  GC_FREE (iwork);

  return (B);
}

/* **************************************************************
 * Solve a symmetric set of equations...
 * ************************************************************** */

MDR *
mdr_SolveEq_SYM (MDR * a, MDR * b)
{
  int info, lda, ldb, mm, n, nrhs, *ipiv;
  int lwork, *iwork;
  int uplo;
  double *anorm, rcond;
  MDR *A, *B, *work;

  mdr_Detect_Inf (a);
  mdr_Detect_Nan (a);

  mdr_Detect_Inf (b);
  mdr_Detect_Nan (b);

  n = lda = mm = MNR (a);
  uplo = (F_INT) 'L';
//   trans = (F_INT) 'N';
//   norm = (F_INT) '1';

  /*
   *Try and pick a good NB, without ILAENV.
   */

  if (n < 100)
    lwork = n;
  else
    lwork = 64 * n;

  if (mm != n)
    rerror ("solve: coefficient matrix must be square");

  if (mm != MNR (b))
    rerror ("solve: RHS row dim. must match LHS row dim.");

  ipiv = (F_INT *) GC_malloc_atomic_ignore_off_page (sizeof (F_INT) * n);
  A = mdr_Float_BF (a);
  B = mdr_Float_BF (b);
  work = mdr_Create (1, lwork);

  signal (SIGINT, intcatch);
  RSYTRF (&uplo, &n, MDRPTR (A), &lda, ipiv, MDRPTR (work), &lwork, &info);
  signal (SIGINT, intcatch_wait);

  if ((int) info < 0)
    rerror ("bad argument to LAPACK DSYTRF");
  if ((int) info > 0)
    rerror ("matrix is singular");

  /*
   * Check input for ill-conditioning.
   */

  mdr_Destroy (work);
  work = mdr_Create (1, 2 * n);
  iwork = (F_INT *) GC_malloc_atomic_ignore_off_page (sizeof (F_INT) * n);

  anorm = mdr_Norm (a, "1");
  signal (SIGINT, intcatch);
  RSYCON (&uplo, &n, MDRPTR (A), &lda, ipiv, anorm, &rcond,
	  MDRPTR (work), iwork, &info);
  signal (SIGINT, intcatch_wait);

  if (rcond <= DBL_EPSILON)
    warning_1 ("WARNING, ill-conditioned input");

  GC_FREE (anorm);
  GC_FREE (iwork);
  nrhs = MNC (b);
  ldb = MNR (b);

  signal (SIGINT, intcatch);
  RSYTRS (&uplo, &n, &nrhs, MDRPTR (A), &lda, ipiv, MDRPTR (B), &ldb, &info);
  signal (SIGINT, intcatch_wait);

  if ((int) info < 0)
    rerror ("Bad argument(s) to LAPACK DSYTRS or ZHETRS");

  GC_FREE (ipiv);
  mdr_Destroy (A);
  mdr_Destroy (work);

  return (B);
}

/* **************************************************************
 * Compute the recipricol of the condition number.
 * ************************************************************** */

double *
mdr_Rcond (MDR * m)
{
  F_INT lda, n, info, *ipiv, *iwork, norm;
  double *anorm;
  double *rcond;
  MDR *A, *work;

  mdr_Detect_Inf (m);
  mdr_Detect_Nan (m);

  lda = MNR (m);
  n = MIN(MNC (m), lda);

  norm = (F_INT) '1';

  work = mdr_Create (1, 4 * n);
  iwork = (F_INT *) GC_malloc_atomic_ignore_off_page (sizeof (F_INT) * n);
  ipiv = (F_INT *) GC_malloc_atomic_ignore_off_page (sizeof (F_INT) * lda);

  // check for integer m
  A = mdr_Float_BF (m);

  anorm = mdr_Norm (m, "1");
  rcond = (double *) GC_MALLOC (sizeof (double));

  signal (SIGINT, intcatch);
  RGETRF (&lda, &n, MDRPTR (A), &lda, ipiv, &info);
  RGECON (&norm, &n, MDRPTR (A), &lda, anorm, rcond, MDRPTR (work),
           iwork, &info);
  signal (SIGINT, intcatch_wait);

  if ((int) info < 0)
    rerror ("illegal argument to LAPACK DGECON");

  mdr_Destroy (work);
  GC_FREE (iwork);
  GC_FREE (ipiv);

  mdr_Destroy (A);

  GC_FREE (anorm);

  return (rcond);
}

/* **************************************************************
 * Compute the norm of a real-matrix.
 * ************************************************************** */

MDR *
mdr_RowNorm (MDR * M, char *type)
{
  double work, *vec;
  int nrow, ncol, itype;
  MDR *m, *norm;
  int i,j, one=1;

  // no integer optimization
  if (M->type == RLAB_TYPE_INT32)
    m = mdr_Float_BF (M);
  else
    m = M;

  mdr_Detect_Inf (m);
  mdr_Detect_Nan (m);

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
    rerror ("rownorm: incorrect STRING specifier");

  /* Need the string length for FORTRAN */
  //   l = strlen (type);
  nrow = MNR (m);
  ncol = MNC (m);

  norm = mdr_Create(nrow,1);

  itype = (int) type[0];

  vec  = (double *) GC_malloc_atomic_ignore_off_page((size_t) (ncol * sizeof (double)));

  for (i=0; i<nrow; i++)
  {
    for (j=0;j<ncol;j++)
      vec[j] = Mdr0(m,i,j);

    signal (SIGINT, intcatch);
    MdrV0(norm,i) = RLANGE (&itype, &one, &ncol, vec, &one, &work);
    signal (SIGINT, intcatch_wait);
  }

  GC_FREE(vec);

  // no integer optimization
  if (M->type == RLAB_TYPE_INT32)
    mdr_Destroy (m);

  return (norm);
}


double *
mdr_Norm (MDR * M, char *type)
{
  F_DOUBLE *norm, *work;
  F_INT lda, nrow, ncol;
  F_INT itype;
  MDR *lsv, *rsv, *sigma, *m;

  // no integer optimization
  if (M->type == RLAB_TYPE_INT32) m = mdr_Float_BF (M);
  else m = M;

  norm = (F_DOUBLE *) GC_MALLOC (sizeof (F_DOUBLE));
  *norm = 0.0;

  mdr_Detect_Inf (m);
  mdr_Detect_Nan (m);

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

  /* Need the string length for FORTRAN */
  //   l = strlen (type);
  nrow = (F_INT) MNR (m);
  ncol = (F_INT) MNC (m);
  lda = nrow;

  if (strcmp (type, "2"))
  {
    work = (double *) GC_malloc_atomic_ignore_off_page
    ((size_t) (nrow * sizeof (double)));
    itype = (F_INT) (type[0] - '0' + '0');

    signal (SIGINT, intcatch);
    *norm = RLANGE (&itype, &nrow, &ncol, MDRPTR (m), &lda, work);
    signal (SIGINT, intcatch_wait);

    GC_FREE (work);
  }
  else
  {
    /* Compute the matrix 2-norm */
    mdr_Svd (m, &rsv, &lsv, &sigma, 3);

    /*
     * Return the largest singular value s[1]
     */

    *norm = (F_DOUBLE) Mdr1 (sigma, 1, 1);

    mdr_Destroy (rsv);
    mdr_Destroy (lsv);
    mdr_Destroy (sigma);
  }

  // no integer optimization
  if (M->type == RLAB_TYPE_INT32)
    mdr_Destroy (m);

  return (norm);
}

/* **************************************************************
 * Matrix P-norm
 * ************************************************************** */

double *
mdr_PNorm (MDR * M, double p)
{
  double *pnorm;
  double sum = 0.0;
  int i, size;

  MDR *m=0;

  // no integer optimization
  if (M->type == RLAB_TYPE_INT32)
    m = mdr_Float_BF (m);
  else
    m = M;

  if (MNR (m) != 1 && MNC (m) != 1)
    rerror ("cannot compute P-norm of a matrix");

  pnorm = (double *) GC_MALLOC (sizeof (double));
  size = MNR (m) * MNC (m);

  if (detect_inf_r (&p, 1))
  {
    double maxr = MdrV0 (m, 0);
    for (i = 1; i < size; i++)
    {
      if (MdrV0 (m, i) > maxr)
        maxr = MdrV0 (m, i);
    }
    *pnorm = maxr;
  }
  else
  {
    for (i = 0; i < size; i++)
      sum += errcheck (pow ((ABS(MdrV0 (m, i))), p), "pow");
    *pnorm = pow (sum, 1.0 / p);
  }

  // no integer optimization
  if (M->type == RLAB_TYPE_INT32)
    mdr_Destroy (m);

  return (pnorm);
}

/* **************************************************************
 * Matrix singular value decompostiion...
 * ************************************************************** */

void
mdr_Svd (MDR * M, MDR ** rsv, MDR ** lsv, MDR ** sigma, int flag)
{
  int k, lda, ldu, ldvt, lwork, m, n, info;
  int jobu, jobvt;
  MDR *s, *work;
  MDR *A, *u, *vt;

  u = vt = 0;			/* Initialize */

  mdr_Detect_Inf (M);
  mdr_Detect_Nan (M);

  m = (int) MNR (M);
  n = (int) MNC (M);
  lda = m;
  k = (int) MIN(m, n);
  // BLAS 2.0, Ian Searle
  // lwork = (int) MAX(3 * MIN(m, n) + MAX(m, n), 5 * MIN(m, n) - 4);
  // BLAS 3.0, Marijan Kostrun
  lwork = (int) MAX(3 * MIN(m, n) + MAX(m, n), 5 * MIN(m, n));
  if (flag == 1)
  {
    jobu = (int) 'A';
    jobvt = (int) 'A';
    ldu = m;
    ldvt = n;
    u = mdr_Create (ldu, m);
    vt = mdr_Create (ldvt, n);
  }
  else if (flag == 2)
  {
    jobu = (int) 'S';
    jobvt = (int) 'S';
    ldu = m;
    ldvt = k;
    u = mdr_Create (ldu, k);
    vt = mdr_Create (ldvt, n);
  }
  else if (flag == 3)
  {
    jobu = (int) 'N';
    jobvt = (int) 'N';
    ldu = 1;
    ldvt = 1;
    u = mdr_Create (0, 0);
    vt = mdr_Create (0, 0);
  }

  // no integer optimization
  A = mdr_Float_BF (M);

  s = mdr_Create (1, k);
  work = mdr_Create (1, lwork);

  signal (SIGINT, intcatch);
  RGESVD (&jobu, &jobvt, &m, &n, MDRPTR (A), &lda, MDRPTR (s),
           MDRPTR (u), &ldu, MDRPTR (vt), &ldvt, MDRPTR (work), &lwork, &info);
  signal (SIGINT, intcatch_wait);

  if ((int) info < 0)
  {
    fprintf (stderr, "ERROR: %ith argument to DGESVD is bad\n", (int) -info);
    rerror ("bad argument to LAPACK DGESVD");
  }
  if ((int) info > 0)
    rerror ("svd algorithm failed to converge");

  /* Clean Up */
  mdr_Destroy (A);
  mdr_Destroy (work);

  /* Set proper addresses */
  *rsv = vt;
  *lsv = u;
  *sigma = s;
}

/* **************************************************************
 * Matrix least squares solution...
 * ************************************************************** */

MDR *
mdr_LS (MDR * a, MDR * b)
{
  F_INT info, lda, ldb, lwork, m, n, nrhs, rank;
  double rcond, mtmp;
  MDR *A, *B, *Btmp, *s, *work;

  m = MNR (a);
  n = MNC (a);
  nrhs = MNC (b);
  lda = m;
  ldb = MAX(m, n);
  rcond = -1.0;

  if (m >= n)
  {
    mtmp = MAX(2 * n, nrhs);
    lwork = (F_INT) (3 * n + MAX(mtmp, m));
  }
  else
  {
    mtmp = MAX(2 * m, nrhs);
    lwork = (F_INT) (3 * m + MAX(mtmp, n));
  }

  s = mdr_Create (1, MIN(m, n));
  work = mdr_Create (1, lwork);

  A = mdr_Float_BF (a);

  /* Check to make sure B is O.K. */
  if (ldb > m)
  {
    int i, j;
    B = mdr_Create (MNR (b) + (ldb - m), MNC (b));
    for (i = 0; i < MNR (b); i++)
    {
      for (j = 0; j < MNC (b); j++)
      {
        if (b->type == RLAB_TYPE_INT32)
          Mdr0 (B, i, j) = Mdi0 (b, i, j);
        else
          Mdr0 (B, i, j) = Mdr0 (b, i, j);
      }
    }
  }
  else
  {
    B = mdr_Float_BF (b);
  }

  signal (SIGINT, intcatch);
  RGELSS (&m, &n, &nrhs, MDRPTR (A), &lda, MDRPTR (B), &ldb, MDRPTR (s),
           &rcond, &rank, MDRPTR (work), &lwork, &info);
  signal (SIGINT, intcatch_wait);

  if ((int) info < 0)
  {
    fprintf (stderr, "ERROR: %ith argument to DGELSS is bad\n", (int) -info);
    rerror ("illegal argument to LAPACK DGELSS()");
  }
  else if ((int) info > 0)
    rerror ("SVD algorithm failed to converge");

  mdr_Destroy (s);
  mdr_Destroy (work);
  mdr_Destroy (A);

  /* re-adjust B is necessary */
  if (m > n)
  {
    int i, j;
    Btmp = mdr_Create (n, MNC (B));
    for (i = 0; i < n; i++)
    {
      for (j = 0; j < MNC (B); j++)
      {
        Mdr0 (Btmp, i, j) = Mdr0 (B, i, j);
      }
    }
    mdr_Destroy (B);
    return (Btmp);
  }
  else
  {
    return (B);
  }
}

/* **************************************************************
 * Driver for Standard eigenvalue problem for real dense matrices.
 * The two options (at this time) are general, and symmetric.
 * ************************************************************** */

void mdr_EigS_SYM (MDR * M, void **val, void **vec);
void mdr_EigS_GE (MDR * M, void **val, void **vec, void **lvec,
                         int lflag);

Btree *
mdr_EigS (MDR * a, int * issym)
{
  Btree *bt;
  Ent *eval, *evec;
  void *val, *vec, *lvec;

  // Create the list that will contain the results.
  bt = btree_Create ();

  //
  // Check for symmetry of [A]
  // Then use the correct routine.
  //
  if (*issym == -1)
    *issym = mdr_IsSymmetric (a);

  if (*issym == 1)
  {
    mdr_EigS_SYM (a, &val, &vec);
    // Hook the results into the list.
    eval = ent_Create ();
    ent_SetType (eval, MATRIX_DENSE_REAL);
    ent_data (eval) = val;
    install (bt, ("val"), eval);

    evec = ent_Create ();
    ent_SetType (evec, MATRIX_DENSE_REAL);
    ent_data (evec) = vec;
    install (bt, ("vec"), evec);
  }
  else
  {
    mdr_EigS_GE (a, &val, &vec, &lvec, 0);
    // Hook the results into the list.
    eval = ent_Create ();
    ent_SetType (eval, MATRIX_DENSE_COMPLEX);
    ent_data (eval) = val;
    install (bt, ("val"), eval);
    ent_Clean(eval);
    evec = ent_Create ();
    ent_SetType (evec, MATRIX_DENSE_COMPLEX);
    ent_data (evec) = vec;
    install (bt, ("vec"), evec);
  }

  //
  // clean-up
  //
  return (bt);
}

/* **************************************************************
 * Compute eigenvalues, and vectors for symmetric prob
 * [A]x = Lambda*x.
 * ************************************************************** */

void
mdr_EigS_SYM (MDR * A, void **val, void **vec)
{
  int info, lda, lwork, m, n, jobz, uplo;
  MDR *work = 0;
  MDR *w = 0;

  mdr_Detect_Inf (A);
  mdr_Detect_Nan (A);

  /* Some rudimentary checks */
  m = MNR (A);
  n = MNC (A);
  if (m != n)
    rerror ("eig must input square");

  lda = m;

  jobz = (int) 'V';
  uplo = (int) 'L';

  if (m != n)
    rerror ("eig: input must be square");

  lwork = 3 * n - 1;
  w = mdr_Create (1, n);
  work = mdr_Create (1, lwork);

  MDR * a = mdr_Float_BF (A);

  RSYEV (&jobz, &uplo, &n, MDRPTR(a), &lda, MDRPTR(w),
          MDRPTR(work), &lwork, &info);

  if ((int) info < 0)
    rerror ("bad argument to LAPACK DSYEV or ZHEEV");
  if ((int) info > 0)
    rerror ("Eigensolver failed to converge");

  *val = w;
  *vec = a;

  mdr_Destroy (work);
  return;
}

/* **************************************************************
 * Compute eigenvalues, and vectors for non-symmetric prob
 * [A]x = Lambda*x.
 * lflag =  0 No left-eigenvectors
 * lflag != 0 Return left-eigenvectors
 * ************************************************************** */

void
mdr_EigS_GE (MDR * A, void **val, void **vec, void **lvec, int lflag)
{
  int i, j, k;
  F_INT info, lda, ldvl, ldvr, lwork, m, n;
  F_INT jobvl, jobvr;
  MDR *wr=0, *wi=0, *vr=0, *vl=0, *work=0;
  MDC *etmp=0, *ptmp=0, *ltmp=0;

  wr = wi = 0;          /* Initialize */
  ltmp = 0;         /* Initialize */

  mdr_Detect_Inf (A);
  mdr_Detect_Nan (A);

  // Some rudimentary checks
  m = MNR (A);
  n = MNC (A);
  lda = m;
  ldvl = n;
  ldvr = n;

  if (m != n)
    rerror ("eig must input square");

  // Copy [a] cause it will get destroyed
  //A = mdr_Copy (M);

  lwork = 4 * n;
  wr = mdr_Create (1, n);
  wi = mdr_Create (1, n);
  vr = mdr_Create (ldvr, n);

  if (lflag == 0)
    vl = mdr_Create (1, 1);
  else
    vl = mdr_Create (ldvr, n);

  work = mdr_Create (1, lwork);

//  signal (SIGINT, intcatch);
  if (lflag == 0)
  {
    jobvl = (F_INT) 'N';
    jobvr = (F_INT) 'V';
    RGEEV (&jobvl, &jobvr, &n, MDRPTR (A), &lda, MDRPTR (wr), MDRPTR (wi),
	   MDRPTR (vl), &ldvl, MDRPTR (vr), &ldvr, MDRPTR (work),
	   &lwork, &info);
  }
  else
  {
    jobvl = (F_INT) 'V';
    jobvr = (F_INT) 'V';
    RGEEV (&jobvl, &jobvr, &n, MDRPTR (A), &lda, MDRPTR (wr), MDRPTR (wi),
	   MDRPTR (vl), &ldvl, MDRPTR (vr), &ldvr, MDRPTR (work),
	   &lwork, &info);
  }
//  signal (SIGINT, intcatch_wait);

  if ((int) info < 0)
    rerror ("bad argument to LAPACK DGEEV or ZGEEV");
  if ((int) info > 0)
    rerror ("failed to compute all of the eigenvalues");


  // First unpack eigenvectors
  // Loop accros rows, get all vectors

  ptmp = mdc_Create (m, n);

  for (k = 0; k < n; k++)
  {
    // Right Eigenvectors
    if (MdrV0 (wi, k) == 0.0)
    {
      // Real Eigenvector
      for (j = 0; j < n; j++)
      {
        Mdc0r (ptmp, j, k) = Mdr0 (vr, j, k);
        Mdc0i (ptmp, j, k) = 0.0;
      }
    }
    else
    {
      // Imaginary Eigenvector
      if (MdrV0 (wi, k) >= 0.0)
      {
        // Construct positive part of complex conjugate pair
        for (j = 0; j < n; j++)
        {
          Mdc0r (ptmp, j, k) = Mdr0 (vr, j, k);
          Mdc0i (ptmp, j, k) = Mdr0 (vr, j, (k + 1));
        }
      }
      else
      {
        // Construct negative part of CC pair
        for (j = 0; j < n; j++)
        {
          Mdc0r (ptmp, j, k) =  Mdc0r (ptmp, j, k - 1);
          Mdc0i (ptmp, j, k) = -Mdc0i (ptmp, j, k - 1);
        }
      }
    }
  }

  if (lflag != 0)
  {
    ltmp = mdc_Create (m, n);
    for (k = 0; k < n; k++)
    {
      // Left Eigenvectors
      if (MdrV0 (wi, k) == 0.0)
      {
    // Real Eigenvector
        for (j = 0; j < n; j++)
        {
          Mdc0r (ltmp, j, k) = Mdr0 (vl, j, k);
          Mdc0i (ltmp, j, k) = 0.0;
        }
      }
      else
      {
    //* Imaginary Eigenvector
        if (MdrV0 (wi, k) >= 0.0)
        {
      // Construct positive part of complex conjugate pair
          for (j = 0; j < n; j++)
          {
            Mdc0r (ltmp, j, k) = Mdr0 (vl, j, k);
            Mdc0i (ltmp, j, k) = Mdr0 (vl, j, (k + 1));
          }
        }
        else
        {
      // Construct negative part of CC pair
          for (j = 0; j < n; j++)
          {
            Mdc0r (ltmp, j, k) = Mdc0r (ltmp, j, k - 1);
            Mdc0i (ltmp, j, k) = -Mdc0i (ltmp, j, k - 1);
          }
        }
      }
    }
  }

  // Now load up eigenvalues
  etmp = mdc_Create (1, MNC (wr));
  for (i = 0; i < MNC (wr); i++)
  {
    MdcV0r (etmp, i) = MdrV0 (wr, i);
    MdcV0i (etmp, i) = MdrV0 (wi, i);
  }

  // Set argument pointers
  *val = etmp;
  *vec = ptmp;

  if (lflag != 0)
    *lvec = ltmp;

  /* Clean Up */
  //mdr_Destroy (A);
  mdr_Destroy (wr);
  mdr_Destroy (wi);
  mdr_Destroy (vr);
  mdr_Destroy (work);
  mdr_Destroy (vl);
}

// *********************************************************************
// * Driver for Generalized eigenvalue problem for general real dense
// * matrices.
// *********************************************************************

void mdr_EigG_GE  (MDR * M, MDR * B, void **val, void **vec);
void mdr_EigG_SYM (MDR * A, MDR * B, void **val, void **vec);

Btree *
mdr_EigG (MDR * A, MDR * B, int * issympos)
{
  Btree *bt;
  Ent *eval, *evec;
  void *val, *vec;
  MDR *a=0, *b=0;

  // no integer matrices
  if (A->type == RLAB_TYPE_INT32) a = mdr_Float_BF (A);
  else a = A;
  if (B->type == RLAB_TYPE_INT32) b = mdr_Float_BF (B);
  else b = B;

  // Create the list that will contain the results.
  bt = btree_Create ();

  if(*issympos==1)
  {
    mdr_EigG_SYM (a, b, &val, &vec);
    // Hook the results into the list.
    eval = ent_Create ();
    ent_SetType (eval, MATRIX_DENSE_REAL);
    ent_data (eval) = val;
    install (bt, ("val"), eval);
    evec = ent_Create ();
    ent_SetType (evec, MATRIX_DENSE_REAL);
    ent_data (evec) = vec;
    install (bt, ("vec"), evec);
  }
  else
  {
    mdr_EigG_GE (a, b, &val, &vec);
    // Hook the results into the list.
    eval = ent_Create ();
    ent_SetType (eval, MATRIX_DENSE_COMPLEX);
    ent_data (eval) = val;
    install (bt, ("val"), eval);
    evec = ent_Create ();
    ent_SetType (evec, MATRIX_DENSE_COMPLEX);
    ent_data (evec) = vec;
    install (bt, ("vec"), evec);
  }
//   ent_Clean( evec );
//   ent_Clean( eval );

  // no integer matrices
  if (A->type == RLAB_TYPE_INT32) mdr_Destroy (a);
  if (B->type == RLAB_TYPE_INT32) mdr_Destroy (b);

  return (bt);
}

// *********************************************************************
// * Driver for Symmetric eigenvalue problem for symmetric real dense
// * matrices 'a' and 'b'. Matrix 'b' has also to be positive definite.
// *********************************************************************

void
mdr_EigG_SYM (MDR * A, MDR * B, void **val, void **vec)
{
  F_INT info, itype, lda, lwork, n, jobz, uplo;
  //MDR *A, *B;
  MDR *w=0, *work=0;
  //A = B = rwork = work = 0;	// Initialize
  mdr_Detect_Inf (A);
  mdr_Detect_Nan (A);
  mdr_Detect_Inf (B);
  mdr_Detect_Nan (B);

  // Some rudimentary checks, MA, MB must be square
  if (MNR (A) != MNC (A))
    rerror ("eig: A must be a square matrix");
  if (MNR (B) != MNC (B))
    rerror ("eig: B must be a square matrix");
  if (MNR (A) != MNR (B))
    rerror ("eig: A and B must be same size");
  itype = 1;
  n = MNR (A);
  lda = n;
  jobz = (F_INT) 'V';
  uplo = (F_INT) 'L';
  w = mdr_Create (1, n);
  lwork = MAX(1, 3 * n - 1);
  //A = mdr_Copy (MA);
  //B = mdr_Copy (MB);
  work = mdr_Create (1, lwork);

  RSYGV (&itype, &jobz, &uplo, &n, MDRPTR(A), &lda, MDRPTR(B), &lda,
          MDRPTR(w), MDRPTR(work), &lwork, &info);

  if ((int) info < 0)
    rerror ("bad argument to LAPACK DSYGV or ZHEGV");
  if ((int) info > 0)
    rerror ("failure in eigensolver DSYGV or ZHEGV");

  *val = w;
  *vec = A;

  /* Clean Up */

  //mdr_Destroy (B);
  mdr_Destroy (work);
}

void
mdr_EigG_GE (MDR * A, MDR * B, void **val, void **vecr)
{
  int i, j;
  F_INT info, k, lda, ldvl, ldvr, lwork, n;
  F_INT jobvl, jobvr;
  MDR *alphar=0, *alphai=0, *beta=0;
  MDR *vl, *vr;
  MDR *work=0, *rwork=0;
  MDC *eval=0, *rvec=0;

  vl = vr = work = rwork = beta = 0;	/* Initialize */
  rvec = 0;

  mdr_Detect_Inf (A);
  mdr_Detect_Nan (A);
  mdr_Detect_Inf (B);
  mdr_Detect_Nan (B);

  // Some rudimentary checks, MA, MB must be square
  if (MNR (A) != MNC (A))
    rerror ("eig: A must be symmetric");
  if (MNR (B) != MNC (B))
    rerror ("eig: B must be symmetric");
  if (MNR (A) != MNR (B))
    rerror ("eig: A and B must be same size");

  n = MNR (A);
  lda = ldvl = ldvr = n;

  jobvl = (F_INT) 'N';
  jobvr = (F_INT) 'V';

  lwork = MAX(1, 8 * n);
  work = mdr_Create (1, lwork);
  alphar = mdr_Create (1, n);
  alphai = mdr_Create (1, n);
  beta = mdr_Create (1, n);
  vl = mdr_Create (1, 1);
  vr = mdr_Create (ldvr, n);

  // lp.f (RLaB to lapack interface) defines RGEGV -> DGGEV
  RGEGV (&jobvl, &jobvr, &n, MDRPTR (A), &lda, MDRPTR (B), &lda,
          MDRPTR (alphar), MDRPTR (alphai), MDRPTR (beta),
          MDRPTR (vl), &ldvl, MDRPTR (vr), &ldvr, MDRPTR (work), &lwork, &info);

  if ((int) info < 0)
    rerror ("bad argument to LAPACK DGGEV or ZGGEV");
  if ((int) info > 0)
    rerror ("failure in eigensolver DGGEV or ZGGEV");

  //
  // Now sort out results
  // First unpack eigenvectors
  //

  rvec = mdc_Create (n, n);

  for (k = 0; k < n; k++)
  {
    if (MdrV0 (alphai, k) == 0.0)
    {
      // Real Eigenvector
      for (j = 0; j < n; j++)
      {
        Mdc0r (rvec, j, k) = Mdr0 (vr, j, k);
        Mdc0i (rvec, j, k) = 0.0;
      }
    }
    else
    {
      // Imaginary Eigenvector
      if (MdrV0 (alphai, k) >= 0.0)
      {
        // Construct positive part of complex conjugate pair
        for (j = 0; j < n; j++)
        {
          Mdc0r (rvec, j, k) = Mdr0 (vr, j, k);
          Mdc0i (rvec, j, k) = Mdr0 (vr, j, (k + 1));
        }
      }
      else
      {
	// Construct negative part of CC pair
	for (j = 0; j < n; j++)
	{
	  Mdc0r (rvec, j, k) = Mdc0r (rvec, j, k - 1);
	  Mdc0i (rvec, j, k) = -Mdc0i (rvec, j, k - 1);
	}
      }
    }
  }
  mdr_Destroy (vl);
  mdr_Destroy (vr);

  //
  // Now compute eigenvalues
  //

  eval = mdc_Create (1, n);

  for (i = 1; i <= n; i++)
  {
    if (Mdr1 (beta, 1, i) == 0.0)
    {
      if (Mdr1 (alphar, 1, i) == 0.0)
	Mdc1r (eval, 1, i) = create_nan ();
      else
	Mdc1r (eval, 1, i) = create_inf ();
      if (Mdr1 (alphai, 1, i) == 0.0)
	Mdc1i (eval, 1, i) = create_nan ();
      else
	Mdc1i (eval, 1, i) = create_inf ();
    }
    else
    {
      Mdc1r (eval, 1, i) = Mdr1 (alphar, 1, i) / Mdr1 (beta, 1, i);
      Mdc1i (eval, 1, i) = Mdr1 (alphai, 1, i) / Mdr1 (beta, 1, i);
    }
  }

  *vecr = rvec;
  *val = eval;

  // Clean Up
  mdr_Destroy (work);
  mdr_Destroy (alphar);
  mdr_Destroy (alphai);
  mdr_Destroy (beta);
}

/* **************************************************************
 * Matrix Factor function:
 * ************************************************************** */

static void mdr_Factor_Ge (MDR * m, MDR ** lu, MDR ** pvt, double *cond);
static void mdr_Factor_Sym (MDR * m, MDR ** ldl, MDR ** pvt, double *cond);

Btree *
mdr_Factor (MDR * A, int flag)
{
  double rcond;
  Btree *bt;
  Ent *elu, *epvt, *ercond;
  MDR *lu, *pvt, *a;

  // not integer optimized
  if (A->type == RLAB_TYPE_INT32) a = mdr_Float_BF (A);
  else a = A;

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

    if (mdr_IsSymmetric (a))
    {
      mdr_Factor_Sym (a, &lu, &pvt, &rcond);

      /* Hook the results into the list. */
      elu = ent_Create ();
      ent_SetType (elu, MATRIX_DENSE_REAL);
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
      mdr_Factor_Ge (a, &lu, &pvt, &rcond);

      /* Hook the results into the list. */
      elu = ent_Create ();
      ent_SetType (elu, MATRIX_DENSE_REAL);
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
    mdr_Factor_Ge (a, &lu, &pvt, &rcond);

    /* Hook the results into the list. */
    elu = ent_Create ();
    ent_SetType (elu, MATRIX_DENSE_REAL);
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
    mdr_Factor_Sym (a, &lu, &pvt, &rcond);

    /* Hook the results into the list. */
    elu = ent_Create ();
    ent_SetType (elu, MATRIX_DENSE_REAL);
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

  // not integer optimized
  if (A->type == RLAB_TYPE_INT32) mdr_Destroy (a);

  return (bt);
}

static void
mdr_Factor_Sym (MDR * m, MDR ** ldl, MDR ** pvt, double *cond)
{
  double *anorm, rcond;
  int i;
  F_INT info, mm, n, lda, *ipiv;
  F_INT *iwork, lwork, uplo;
  MDR *tldl, *tpvt, *work;

  uplo = (F_INT) 'L';
  lda = n = mm = MNR (m);
//   one = 1;
//   norm = (F_INT) '1';

  /*
   *Try and pick a good NB, without ILAENV.
   */

  if (n < 100)
    lwork = n;
  else
    lwork = 64 * n;

  if (MNR (m) != MNC (m))
    rerror ("Factor: coefficient matrix must be square");

  ipiv = (F_INT *) GC_malloc_atomic_ignore_off_page (sizeof (F_INT) * mm);
  tpvt = mdr_Create (1, mm);
  tldl = mdr_Float_BF (m);
  iwork = (F_INT *) GC_malloc_atomic_ignore_off_page (sizeof (F_INT) * n);
  work = mdr_Create (1, lwork);

  signal (SIGINT, intcatch);
  RSYTRF (&uplo, &n, MDRPTR (tldl), &lda, ipiv, MDRPTR (work), &lwork, &info);
  signal (SIGINT, intcatch_wait);

  if ((int) info < 0)
    rerror ("bad argument to LAPACK DSYTRF");
  if ((int) info > 0)
    rerror ("matrix is singular");

  /*
   * Check input for ill-conditioning.
   */

  anorm = mdr_Norm (m, "1");
  mdr_Destroy (work);
  work = mdr_Create (1, 2 * n);

  signal (SIGINT, intcatch);
  RSYCON (&uplo, &n, MDRPTR (tldl), &lda, ipiv, anorm, &rcond,
	  MDRPTR (work), iwork, &info);
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
  GC_FREE (iwork);
  mdr_Destroy (work);
}

static void
mdr_Factor_Ge (MDR * m, MDR ** lu, MDR ** pvt, double *cond)
{
  double *anorm, rcond;
  int i;
  F_INT info, mm, n, lda, *ipiv;
  F_INT *iwork, lwork, norm;
  MDR *tlu, *tpvt, *work;

  lwork = lda = mm = MNR (m);
  n = MNC (m);
  norm = (F_INT) '1';

  if (mm != n)
    rerror ("factor: coefficient matrix must be square");

  ipiv = (F_INT *) GC_malloc_atomic_ignore_off_page (sizeof (F_INT) * mm);
  tpvt = mdr_Create (1, mm);
  tlu = mdr_Float_BF (m);
  iwork = (F_INT *) GC_malloc_atomic_ignore_off_page (sizeof (F_INT) * n);
  work = mdr_Create (1, 4 * lwork);

  signal (SIGINT, intcatch);
  RGETRF (&mm, &n, MDRPTR (tlu), &lda, ipiv, &info);
  signal (SIGINT, intcatch_wait);

  if ((int) info < 0)
    rerror ("bad argument to LAPACK DGETRF");
  if ((int) info > 0)
    rerror ("matrix is singular");

  /*
   * Check input for ill-conditioning.
   */

  anorm = mdr_Norm (m, "1");
  signal (SIGINT, intcatch);
  RGECON (&norm, &n, MDRPTR (tlu), &lda, anorm, &rcond, MDRPTR (work),
	  iwork, &info);
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
  GC_FREE (iwork);
  mdr_Destroy (work);
}

/* **************************************************************
 * Matrix Backsub function.
 * ************************************************************** */

MDR *
mdr_Backsub_Sym (Btree * bl, MDR * b)
{
  int i;
  int *ipiv, lda, ldb, mm, n, nrhs, info;
  int uplo;
  Ent *etmp;
  ListNode *tmp;
  MDR *B, *lu, *pvt, *LU=0;
  lu = pvt = 0;

  /* Get the necessary elements off of the list. */

  if ((tmp = btree_FindNode (bl, "ldl")))
  {
    etmp = var_ent (tmp);
    if (ent_type (etmp) != MATRIX_DENSE_REAL)
      rerror ("backsub: terrible input list error");
    LU = ent_data (etmp);
    if (LU->type == RLAB_TYPE_INT32) lu = mdr_Float_BF (LU);
    else lu = LU;
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
//   trans = (F_INT) 'N';
  uplo = (F_INT) 'L';

  if (ldb != lda)
    rerror ("backsub: b must have same number of rows as LHS");

  ipiv = (F_INT *) GC_malloc_atomic_ignore_off_page (sizeof (F_INT) * mm);
  B = mdr_Float_BF (b);

  /* Fill ipiv */
  for (i = 0; i < mm; i++)
  {
    ipiv[i] = (F_INT) Mdr0 (pvt, 0, i);
  }

  /* Call dsytrs */
  signal (SIGINT, intcatch);
  RSYTRS (&uplo, &n, &nrhs, MDRPTR (lu), &lda, ipiv, MDRPTR (B), &ldb, &info);
  signal (SIGINT, intcatch_wait);
  if ((int) info < 0)
    rerror ("Bad argument(s) to LAPACK DSYTRF");

  if (LU->type == RLAB_TYPE_INT32) mdr_Destroy (lu);
  GC_FREE (ipiv);
  return (B);
}

MDR *
mdr_Backsub_Ge (Btree * bl, MDR * b)
{
  int i;
  int *ipiv, lda, ldb, mm, n, nrhs, info;
  int trans;
  Ent *etmp;
  ListNode *tmp;
  MDR *B, *lu, *pvt, *LU=0;
  lu = pvt = 0;

  /* Get the necessary elements off of the list. */

  if ((tmp = btree_FindNode (bl, "lu")))
  {
    etmp = var_ent (tmp);
    if (ent_type (etmp) != MATRIX_DENSE_REAL)
      rerror ("backsub: terrible input list error");
    LU = ent_data (etmp);
    if (LU->type == RLAB_TYPE_INT32)
      lu = mdr_Float_BF (LU);
    else
      lu = LU;
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
  trans = (F_INT) 'N';
//   uplo = (F_INT) 'L';

  if (ldb != lda)
    rerror ("backsub: b must have same number of rows as LHS");

  ipiv = (F_INT *) GC_malloc_atomic_ignore_off_page (sizeof (F_INT) * mm);
  B = mdr_Float_BF (b);

  /* Fill ipiv */
  for (i = 0; i < mm; i++)
  {
    ipiv[i] = (F_INT) Mdr0 (pvt, 0, i);
  }

  /* Call dgetrs */
  signal (SIGINT, intcatch);
  RGETRS (&trans, &n, &nrhs, MDRPTR (lu), &lda, ipiv, MDRPTR (B), &ldb, &info);
  signal (SIGINT, intcatch_wait);
  if ((int) info < 0)
    rerror ("Bad argument(s) to LAPACK DGETRF");

  if (LU->type == RLAB_TYPE_INT32)
    mdr_Destroy (lu);

  GC_FREE (ipiv);

  return (B);
}

/* **************************************************************
 * Matrix SVD function.
 * ************************************************************** */

Btree *
mdr_Svd_BF (MDR * M, int flag)
{
  Btree *bt;
  MDR *rsv, *lsv, *sigma, *m;

  if (M->type == RLAB_TYPE_INT32)
    m = mdr_Float_BF (M);
  else
    m = M;

  mdr_Svd (m, &rsv, &lsv, &sigma, flag);

  /* Hook the results into the list. */

  bt = btree_Create ();
  install (bt, "u", ent_Assign_Rlab_MDR(lsv));
  install (bt, "sigma", ent_Assign_Rlab_MDR(sigma));
  install (bt, "vt", ent_Assign_Rlab_MDR(rsv));

  if (M->type == RLAB_TYPE_INT32)
    mdr_Destroy (m);

  return (bt);
}

/* **************************************************************
 * Matrix Cholesky decomposition function.
 * ************************************************************** */
MDR *
mdr_Chol (MDR * m)
{
  F_INT i, j, info, lda, n, uplo;
  MDR *A;

  mdr_Detect_Inf (m);
  mdr_Detect_Nan (m);

  lda = n = MNC (m);
  uplo = (F_INT) 'U';

  /* [m] gets overwritten, so copy it */
  A = mdr_Float_BF (m);

  /* Call LAPACK routine */

  signal (SIGINT, intcatch);
  RPOTRF (&uplo, &n, MDRPTR (A), &lda, &info);
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
      Mdr0 (A, i, j) = 0.0;
    }
  }

  return (A);
}

/* **************************************************************
 * Matrix Hess function.
 * ************************************************************** */

Btree *
mdr_Hess (MDR * M)
{
  int i, j;
  F_INT ilo, ihi, lda, m, n, info, lwork;
  MDR *A, *tau, *work;

  Btree *bt;
  Ent *ep, *eh;
  MDR *h;

  mdr_Detect_Inf (M);
  mdr_Detect_Nan (M);

  /* Get matrix dimensions */
  m = MNR (M);
  n = MNC (M);
  lda = m;
  lwork = n;
  ilo = 1;
  ihi = n;

  /* The input to orthes() will be destroyed, so copy it */
  A = mdr_Float_BF (M);
  tau = mdr_Create (1, n);
  work = mdr_Create (1, lwork);

  signal (SIGINT, intcatch);
  RGEHRD (&n, &ilo, &ihi, MDRPTR (A), &lda, MDRPTR (tau),
	  MDRPTR (work), &lwork, &info);
  signal (SIGINT, intcatch_wait);

  if ((int) info < 0)
    rerror ("bad argument to LAPACK DGEHRD");

  /* Save [h] */
  h = mdr_Float_BF (A);
  for (i = 2; i < m; i++)
  {
    for (j = 0; (j < i - 1) && (j < n); j++)
    {
      Mdr0 (h, i, j) = (double) 0.0;
    }
  }

  signal (SIGINT, intcatch);
  RORGHR (&n, &ilo, &ihi, MDRPTR (A), &lda, MDRPTR (tau),
	  MDRPTR (work), &lwork, &info);
  signal (SIGINT, intcatch_wait);

  if ((int) info < 0)
    rerror ("bad argument to LAPACK DORGHR");

  /* clean up */

  mdr_Destroy (tau);
  mdr_Destroy (work);

  bt = btree_Create ();

  ep = ent_Create ();
  ent_SetType (ep, MATRIX_DENSE_REAL);
  ent_data (ep) = A;
  install (bt, ("p"), ep);

  eh = ent_Create ();
  ent_SetType (eh, MATRIX_DENSE_REAL);
  ent_data (eh) = h;
  install (bt, ("h"), eh);

  return (bt);
}

/* **************************************************************
 * Balance a matrix.
 * ************************************************************** */

Btree *
mdr_Balance (MDR * m)
{
  /* Matrix *m, **Ab, **t; */

  double *scale;
  int i;
  F_INT info, nm, n, low, igh, job;
  MDR *A, *T;
  Btree *bt;
  Ent *et, *eab;

  mdr_Detect_Inf (m);
  mdr_Detect_Nan (m);

  if (MNR (m) != MNC (m))
    rerror ("balance: input to must be a square matrix");

  nm = (F_INT) MNR (m);
  n = nm;
  job = (F_INT) 'S';

  A = mdr_Float_BF (m);
  scale = (double *) GC_malloc_atomic_ignore_off_page (n * sizeof (double));

  /* LAPACK balance */

  signal (SIGINT, intcatch);
  RGEBAL (&job, &n, MDRPTR (A), &n, &low, &igh, scale, &info);
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
  ent_SetType (eab, MATRIX_DENSE_REAL);
  ent_data (eab) = A;
  install (bt, ("ab"), eab);

  return (bt);
}

Btree *
mdr_QR (MDR * M)
{
  int i, j;
  double *tau, *work;
  F_INT info, k, lwork, m, n;
  MDR *A, *q, *R;

  Btree *bt;
  Ent *eq, *er;

  mdr_Detect_Inf (M);
  mdr_Detect_Nan (M);

  /* Set dimensional paramaters */
  m = MNR (M);
  n = MNC (M);
  k = MIN(m, n);
  lwork = MAX(m, n);

  /*
   * over/under determined q: M x M
   *                       r: M x N
   *                       p: N x N
   */

  /* Create work arrays */
  A = mdr_Float_BF (M);
  tau = (double *) GC_malloc_atomic_ignore_off_page (sizeof (double) * k);
  work = (double *) GC_malloc_atomic_ignore_off_page (sizeof (double) * lwork);

  signal (SIGINT, intcatch);
  RGEQRF (&m, &n, MDRPTR (A), &m, tau, work, &lwork, &info);
  signal (SIGINT, intcatch_wait);

  if ((int) info < 0)
    rerror ("bad argument to LAPACK DGEQRF");

  /* Extract [R] */
  R = mdr_Create (m, n);
  mdr_Zero (R);

  for (j = 1; j <= MNC (R); j++)
  {
    for (i = 1; (i <= j) && (i <= MNR (R)); i++)
    {
      Mdr1 (R, i, j) = Mdr1 (A, i, j);
    }
  }

  signal (SIGINT, intcatch);
  /* Grow [A] if necessary to hold computed [Q] */
  if ((m - n) > 0)
  {
    /* matrix_AppendColR (A, m - n); */
    mdr_Extend (A, MNR (A), MNC (A) + (m - n));
  }

  RORGQR (&m, &m, &k, MDRPTR (A), &m, tau, work, &lwork, &info);
  signal (SIGINT, intcatch_wait);

  if ((int) info < 0)
    rerror ("bad argument to LAPACK DORGQR");

  /*
   * Shrink Q if necessary
   */

  if ((m - n) < 0)
  {
    /* Get the first M columns of A. */

    q = mdr_Create (MNR (A), m);
    for (i = 0; i < MNR (A); i++)
    {
      for (j = 0; j < m; j++)
      {
	Mdr0 (q, i, j) = Mdr0 (A, i, j);
      }
    }
    mdr_Destroy (A);
  }
  else
  {
    q = A;
  }

  GC_FREE (tau);
  GC_FREE (work);

  bt = btree_Create ();

  eq = ent_Create ();
  ent_SetType (eq, MATRIX_DENSE_REAL);
  ent_data (eq) = q;
  install (bt, ("q"), eq);

  er = ent_Create ();
  ent_SetType (er, MATRIX_DENSE_REAL);
  ent_data (er) = R;
  install (bt, ("r"), er);

  return (bt);
}

Btree *
mdr_QRP (MDR * M)
{
  /* Matrix *M, **q, **r, **p; */

  int i, j, *jpvt;
  double *tau, *work;
  F_INT info, k, lwork, m, n;
  MDR *A, *R, *P, *q;

  Btree *bt;
  Ent *eq, *er, *ep;

  mdr_Detect_Inf (M);
  mdr_Detect_Nan (M);

  /* Set dimensional paramaters */
  m = MNR (M);
  n = MNC (M);
  k = MIN(m, n);
  lwork = 3 * MAX(m, n);

  /*
   *  over / under determined q: M x M
   *                          r: M x N
   *                          p: N x N
   */

  /* Create work arrays */
  A = mdr_Float_BF (M);
  tau = (double *) GC_malloc_atomic_ignore_off_page (sizeof (double) * k);
  work = (double *) GC_malloc_atomic_ignore_off_page (sizeof (double) * lwork);
  jpvt = (int *) GC_malloc_atomic_ignore_off_page (sizeof (int) * n);

  for (i = 0; i < n; i++)
  {
    jpvt[i] = 0;		/* It is important to zero-out jpvt */
  }

  signal (SIGINT, intcatch);
  RGEQPF (&m, &n, MDRPTR (A), &m, jpvt, tau, work, &info);
  signal (SIGINT, intcatch_wait);

  if ((int) info < 0)
    rerror ("bad argument to LAPACK DGEQPF");

  /* Extract [R] */
  R = mdr_Create (m, n);
  mdr_Zero (R);

  for (j = 1; j <= MNC (R); j++)
  {
    for (i = 1; (i <= j) && (i <= MNR (R)); i++)
    {
      Mdr1 (R, i, j) = Mdr1 (A, i, j);
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

  signal (SIGINT, intcatch);
  if ((m - n) > 0)
  {
    /* matrix_AppendColR (A, m - n); */
    mdr_Extend (A, MNR (A), MNC (A) + (m - n));
  }

  RORGQR (&m, &m, &k, MDRPTR (A), &m, tau, work, &lwork, &info);
  signal (SIGINT, intcatch_wait);

  if ((int) info < 0)
    rerror ("bad argument to LAPACK DORGQR");

  if ((m - n) < 0)
  {
    /* Get the first M columns of A. */

    q = mdr_Create (MNR (A), m);
    for (i = 0; i < MNR (A); i++)
    {
      for (j = 0; j < m; j++)
      {
	Mdr0 (q, i, j) = Mdr0 (A, i, j);
      }
    }
    mdr_Destroy (A);
  }
  else
  {
    q = A;
  }

  /* Clean - Up */
  GC_FREE (tau);
  GC_FREE (work);
  GC_FREE (jpvt);

  bt = btree_Create ();

  eq = ent_Create ();
  ent_SetType (eq, MATRIX_DENSE_REAL);
  ent_data (eq) = q;
  install (bt, ("q"), eq);

  er = ent_Create ();
  ent_SetType (er, MATRIX_DENSE_REAL);
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
mdr_Schur (MDR * m)
{
  int *bwork;
  Btree *bt;
  Ent *et, *ez;
  F_INT info, lda, lwork, n, sdim, jobvs, sort;
  MDR *a, *vs, *wr, *wi, *work;

  /* m must be n-by-n */
  if (MNR (m) != MNC (m))
    rerror ("schur: Argument must be square (N-by-N)");

  n = (F_INT) MNR (m);
  lda = n;
  sdim = 0;
  lwork = MAX(1, 3 * n);
  bwork = (int *) GC_MALLOC (2 * sizeof (int));
  jobvs = (F_INT) 'V';
  sort = (F_INT) 'N';

  a = mdr_Float_BF (m);
  wr = mdr_Create (1, n);
  wi = mdr_Create (1, n);
  vs = mdr_Create (n, n);
  work = mdr_Create (1, lwork);

  signal (SIGINT, intcatch);
  RGEES (&jobvs, &sort, rselect, &n, MDRPTR (a), &lda, &sdim, MDRPTR (wr),
	 MDRPTR (wi), MDRPTR (vs), &n, MDRPTR (work), &lwork, bwork, &info);
  signal (SIGINT, intcatch_wait);

  if ((int) info < 0)
  {
    rerror ("illegal argument to DGEES");
  }
  else if ((int) info > 0)
  {
    if ((int) info <= n)
    {
      fprintf (stderr, "schur (DGEES): Failed to compute all eigenvalues\n");
      fprintf (stderr, "               %i:%i converged eigenvalues\n",
	       (int) info, (int) n);
    }
    else if ((int) info == n + 1)
    {
      fprintf (stderr, "schur (DGEES): Eigenvalues could not be re-ordered\n");
      fprintf (stderr, "               This problem is very ill-conditioned\n");
      rerror ("schur: cannot continue");
    }
    else if ((int) info == n + 2)
    {
      fprintf (stderr, "schur (DGEES): Problems re-ordering\n");
      rerror ("schur: cannot continue");
    }
  }

  mdr_Destroy (wr);
  mdr_Destroy (wi);
  mdr_Destroy (work);
  GC_FREE (bwork);

  bt = btree_Create ();

  et = ent_Create ();
  ent_SetType (et, MATRIX_DENSE_REAL);
  ent_data (et) = a;
  install (bt, ("t"), et);

  ez = ent_Create ();
  ent_SetType (ez, MATRIX_DENSE_REAL);
  ent_data (ez) = vs;
  install (bt, ("z"), ez);

  return (bt);
}

/* **************************************************************
 * Solve the Sylvester Matrix Equation.
 * ************************************************************** */

MDR *
mdr_Sylv (MDR * A, MDR * B, MDR * D)
{
  double scale;
//   int lflag;
  F_INT info, isgn, m, n, lda, ldb, ldc;
  F_INT trana, tranb;
  MDR *a, *b, *c, *C;

  // no integer optimization
  if (A->type == RLAB_TYPE_INT32) a = mdr_Float_BF (A);
  else a = A;
  if (B->type == RLAB_TYPE_INT32) b = mdr_Float_BF (B);
  else b = B;
  if (D->type == RLAB_TYPE_INT32) c = mdr_Float_BF (D);
  else c = D;


  if (b == 0)
  {
//     lflag = 1;			/* Lyapunov equation */
    isgn = 1;			/* A*X + X*A' = -C */
    scale = 1.0;
    trana = (F_INT) 'N';
    tranb = (F_INT) 'C';

    /* Set & Check dimensions */
    m = lda = MDR_nrow (a);
    n = ldb = m;
    ldc = m;
    if (m != MDR_ncol (a))
      rerror ("sylv: A.nr must equal A.nc");
    if (m != MDR_nrow (c))
      rerror ("sylv: A.nr must equal C.nr");

    C = mdr_Negate (c);

    signal (SIGINT, intcatch);
    RTRSYL (&trana, &tranb, &isgn, &m, &n, MDRPTR (a), &lda,
	    MDRPTR (a), &ldb, MDRPTR (C), &ldc, &scale, &info);
    signal (SIGINT, intcatch_wait);
  }
  else
  {
//     lflag = 0;			/* Sylvester equation */
    isgn = 1;			/* A*X + X*B = C */
    scale = 1.0;

    /* Set & Check dimensions */
    m = lda = MDR_nrow (a);
    n = ldb = MDR_nrow (b);
    ldc = m;
    trana = (F_INT) 'N';
    tranb = (F_INT) 'N';

    if (m != MDR_ncol (a))
      rerror ("sylv: A.nr must equal A.nc");
    if (n != MDR_ncol (b))
      rerror ("sylv: B.nr must equal B.nc");
    if (m != MDR_nrow (c))
      rerror ("sylv: A.nr must equal C.nr");
    if (n != MNC (c))
      rerror ("sylv: B.nr must equal C.nc");

    C = mdr_Negate (c);

    signal (SIGINT, intcatch);
    RTRSYL (&trana, &tranb, &isgn, &m, &n, MDRPTR (a), &lda,
	    MDRPTR (b), &ldb, MDRPTR (C), &ldc, &scale, &info);
    signal (SIGINT, intcatch_wait);
  }

  if ((int) info < 0)
  {
    fprintf (stderr, "sylv: %ith argument to DTRSYL illegal\n", (int) info);
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

MDR *
mdr_Det (MDR * m)
{
  int i;
  F_INT info, mm, n, lda, *ipiv;
  double det[2], ten, rdv;
  MDR *lu, *mrv;

  lda = mm = MNR (m);
  n = MNC (m);
//   one = 1;

  if (mm != n)
    rerror ("det: matrix argument  must be square");

  ipiv = (F_INT *) GC_malloc_atomic_ignore_off_page (sizeof (F_INT) * mm);
  lu = mdr_Float_BF (m);

  RGETRF (&mm, &n, MDRPTR (lu), &lda, ipiv, &info);

  if ((int) info < 0)
    rerror ("det: bad argument to LAPACK DGETRF");
  if ((int) info > 0)
  {
    warning_1 ("det: matrix argument  is singular");
    mrv = mdr_CreateScalar (0.0);
    return (mrv);
  }

  /*
   * Now compute the determinant using the algorithm contained
   * in the LINPACK subroutines DGEDI, ZGEDI. I'm doing it this
   * way because the LAPACK does not provide a direct method for
   * computing the determinant?
   */

  det[0] = 1.;
  det[1] = 0.;
  ten = 10.;
  for (i = 1; i <= n; ++i)
  {
    if (ipiv[i - 1] != i)
    {
      det[0] = -det[0];
    }
    det[0] = Mdr1 (lu, i, i) * det[0];
    /*        ...exit */
    if (det[0] == 0.)
    {
      goto L60;
    }
  L10:
    if (ABS(det[0]) >= 1.)
    {
      goto L20;
    }
    det[0] = ten * det[0];
    det[1] += -1.;
    goto L10;
  L20:
  L30:
    if (ABS(det[0]) < ten)
    {
      goto L40;
    }
    det[0] /= ten;
    det[1] += 1.;
    goto L30;
  L40:
    ;
  }
L60:

  rdv = det[0] * pow (10.0, det[1]);
  mrv = mdr_CreateScalar (rdv);

  GC_FREE (ipiv);
  mdr_Destroy (lu);

  return (mrv);
}
