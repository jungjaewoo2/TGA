/* mscf2.c Matrix Sparse Complex Functions ... */

/*  This file is a part of RLaB ("Our"-LaB)
   Copyright (C) 1996  Ian R. Searle

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
#include "class.h"
#include "mem.h"
#include "msr.h"
#include "msc.h"
#include "mscf1.h"
#include "mdr.h"
#include "mdr_mdc.h"
#include "mdc.h"
#include "mdcf1.h"
#include "mds.h"
#include "btree.h"
#include "util.h"
#include "symbol.h"
#include "sort.h"

#include "mdrf1.h"
#include "msrf1.h"

#include "fi.h"
#include "sparse.h"

#ifdef HAVE_SUITESPARSE
#include <suitesparse/umfpack.h>
#endif

#include <stdio.h>
#include <math.h>


//
// choice of method for solve function on sparse complex matrices:
//      0 - UMFPACK  (umf)
//      1 - SuperLU  (slu)
int msc_Solve_method = 2;


#ifdef HAVE_SUITESPARSE

// ***************************************************************
// * Solve a set of equations with a sparse coefficient matrix
// * All we do here is solve the general sparse matrix.
// * Use the UMFPACK routines.
// ***************************************************************
//
// by Marijan Kostrun, IV/2005
//
MDC *
umfpack_msc_Det (MSC * a)
{
  int i, j, n, ne;
  MDC *x;
  MDR *Ax, *Az;
  //
  // Initialize UMFPACK
  //
  double Control[UMFPACK_CONTROL], Info[UMFPACK_INFO];
  void *Symbolic, *Numeric;
  double Mx, Mz, Ex;

  msc_Detect_Inf (a);
  msc_Detect_NaN (a);
  if (a->nr != a->nc)
    rerror ("det[umf]: requires square matrix");
  n = a->nr;
  ne = a->nnz;

  // this is what we pass instead of 'a'
  Ax = mdr_Create (1, ne);
  Az = mdr_Create (1, ne);
  for (i = 0; i < ne; i++)
  {
    MdrV0(Ax,i) = RE(a->c[i]);
    MdrV0(Az,i) = IM(a->c[i]);
  }

  // this is silly but I don;t know better: indexing is offset by one
  for (i = 0; i < a->nr + 1; i++)
    --a->ia[i];
  for (j = 0; j < a->nnz; j++)
    --a->ja[j];

  //
  // Factor the coefficient matrix: first symbolic then numeric
  //
  signal (SIGINT, intcatch);
  (void) umfpack_zi_symbolic (n, n, a->ia, a->ja, MDRPTR(Ax), MDRPTR(Az), &Symbolic,
			      Control, Info);
  signal (SIGINT, intcatch_wait);

  signal (SIGINT, intcatch);
  (void) umfpack_zi_numeric (a->ia, a->ja, MDRPTR(Ax), MDRPTR(Az), Symbolic, &Numeric,
			     Control, Info);
  signal (SIGINT, intcatch_wait);

  // Find the determinant
  signal (SIGINT, intcatch);
  (void) umfpack_zi_get_determinant (&Mx, &Mz, &Ex, Numeric, Info);
  signal (SIGINT, intcatch_wait);

  umfpack_zi_free_symbolic (&Symbolic);

  // this is silly but I don;t know better: indexing is offset by one
  for (i = 0; i < a->nr + 1; i++)
    ++a->ia[i];
  for (j = 0; j < a->nnz; j++)
    ++a->ja[j];

  if (Info[UMFPACK_STATUS] == UMFPACK_OK)
    x = mdc_CreateScalar (Mx * pow (10, Ex), Mz * pow (10, Ex));
  else if (Info[UMFPACK_STATUS] == UMFPACK_WARNING_singular_matrix)
    x = mdc_CreateScalar (0.0, 0.0);
  else if (Info[UMFPACK_STATUS] == UMFPACK_WARNING_determinant_underflow ||
	   Info[UMFPACK_STATUS] == UMFPACK_WARNING_determinant_overflow)
  {
    x = mdc_Create (1, 2);
    Mdc1r (x, 1, 1) = Mx;
    Mdc1i (x, 1, 1) = Mz;
    Mdc1r (x, 1, 2) = Ex;
    Mdc1i (x, 1, 2) = 0;
  }
  else
    x = mdc_Create (0, 0);
  return (x);
}


MDC *
umfpack_msc_SolveCC (MSC * a, MDC * b, char *type)
{
  int i, j, n, ne;
  MDC *x;
  MDR *Xx, *Xz, *Bx, *Bz, *Ax, *Az;
  //
  // Initialize UMFPACK
  //
  double Control[UMFPACK_CONTROL], Info[UMFPACK_INFO];
  void *Symbolic, *Numeric;

  msc_Detect_Inf (a);
  msc_Detect_NaN (a);
  mdc_Detect_Inf (b);
  mdc_Detect_Nan (b);
  if (a->nr != a->nc)
    rerror ("solve[umf]: requires square coefficient matrix");
  if (a->nr != MNR (b))
    rerror ("solve[umf]: RHS row dim. must match LHS row dim.");
  n = a->nr;
  ne = a->nnz;

  // create matrix structures for UMFPACK
  Xx = mdr_Create (n, 1);
  Xz = mdr_Create (n, 1);
  Bx = mdr_Create (n, 1);
  Bz = mdr_Create (n, 1);

  // this is what we pass instead of 'a'
  Ax = mdr_Create (1, ne);
  Az = mdr_Create (1, ne);
  for (i = 0; i < ne; i++)
  {
    MdrV0(Ax,i) = RE(a->c[i]);
    MdrV0(Az,i) = IM(a->c[i]);
  }

  // output
  x = mdc_Create (MNR (b), MNC (b));


  // this is silly but I don;t know better: indexing is offset by one
  for (i = 0; i < a->nr + 1; i++)
    --a->ia[i];
  for (j = 0; j < a->nnz; j++)
    --a->ja[j];

  //
  // Factor the coefficient matrix: first symbolic then numeric
  //
  signal (SIGINT, intcatch);
  (void) umfpack_zi_symbolic (n, n, a->ia, a->ja, MDRPTR(Ax), MDRPTR(Az), &Symbolic,
			      Control, Info);
  signal (SIGINT, intcatch_wait);

  signal (SIGINT, intcatch);
  (void) umfpack_zi_numeric (a->ia, a->ja, MDRPTR(Ax), MDRPTR(Az), Symbolic, &Numeric,
			     Control, Info);
  signal (SIGINT, intcatch_wait);

  //
  // Solve Ax = b
  //
  //
  // Loop over the columns of B, getting a solution for each.
  //
  for (j = 0; j < MNC (b); j++)
  {
    // Grab i-th column of 'b'
    for (i = 0; i < n; i++)
    {
      Mdr0 (Bx, i, 0) = RE(Mdc0 (b, i, j));
      Mdr0 (Bz, i, 0) = IM(Mdc0 (b, i, j));
    }
    // Solve, but remember A is transposed (MSR is row compressed format)
    signal (SIGINT, intcatch);
    (void) umfpack_zi_solve (UMFPACK_Aat, a->ia, a->ja, MDRPTR(Ax), MDRPTR(Az),
     MDRPTR(Xx), MDRPTR(Xz), MDRPTR(Bx), MDRPTR(Bz), Numeric, NULL, NULL);
    signal (SIGINT, intcatch_wait);
    // Load result into X
    for (i = 0; i < MNR (b); i++)
    {
      Mdc0(x, i, j) = Mdr0 (Xx, i, 0) + Mdr0 (Xz, i, 0)*1i;
    }
  }

  // clean all
  mdr_Destroy (Xx);
  mdr_Destroy (Xz);
  mdr_Destroy (Bx);
  mdr_Destroy (Bz);
  mdr_Destroy (Ax);
  mdr_Destroy (Az);

  umfpack_zi_free_symbolic (&Symbolic);
  umfpack_zi_free_numeric (&Numeric);

  // this is silly but I don;t know better: indexing is offset by one
  for (i = 0; i < a->nr + 1; i++)
    ++a->ia[i];
  for (j = 0; j < a->nnz; j++)
    ++a->ja[j];

  return (x);
}

MDC *
umfpack_msc_SolveRC (MSR * a, MDC * b, char *type)
{
  int i, j, n, ne;
  MDC *x;
  MDR *Xx, *Xz, *Bx, *Bz, *Az;
  //
  // Initialize UMFPACK
  //
  double Control[UMFPACK_CONTROL], Info[UMFPACK_INFO];
  void *Symbolic, *Numeric;

  msr_Detect_Inf (a);
  msr_Detect_NaN (a);
  mdc_Detect_Inf (b);
  mdc_Detect_Nan (b);
  if (a->nr != a->nc)
    rerror ("solve[umf]: requires square coefficient matrix");
  if (a->nr != MNR (b))
    rerror ("solve[umf]: RHS row dim. must match LHS row dim.");
  n = a->nr;
  ne = a->nnz;

  // create matrix structures for UMFPACK
  Xx = mdr_Create (n, 1);
  Xz = mdr_Create (n, 1);
  Bx = mdr_Create (n, 1);
  Bz = mdr_Create (n, 1);

  // this is what we pass instead of 'a'
  Az = mdr_Create (1, ne);
  mdr_Zero (Az);

  // output
  x = mdc_Create (MNR (b), MNC (b));

  // this is silly but I don;t know better: indexing is offset by one
  for (i = 0; i < a->nr + 1; i++)
    --a->ia[i];
  for (j = 0; j < a->nnz; j++)
    --a->ja[j];

  //
  // Factor the coefficient matrix: first symbolic then numeric
  //
  signal (SIGINT, intcatch);
  (void) umfpack_zi_symbolic (n, n, a->ia, a->ja, a->d, MDRPTR(Az), &Symbolic,
			      Control, Info);
  signal (SIGINT, intcatch_wait);

  signal (SIGINT, intcatch);
  (void) umfpack_zi_numeric (a->ia, a->ja, a->d, MDRPTR(Az), Symbolic, &Numeric,
			     Control, Info);
  signal (SIGINT, intcatch_wait);

  //
  // Solve Ax = b
  //
  //
  // Loop over the columns of B, getting a solution for each.
  //
  for (j = 0; j < MNC (b); j++)
  {
    // Grab i-th column of 'b'
    for (i = 0; i < n; i++)
    {
      Mdr0 (Bx, i, 0) = RE(Mdc0(b, i, j));
      Mdr0 (Bz, i, 0) = IM(Mdc0(b, i, j));
    }
    // Solve, but remember A is transposed (MSR is row compressed format)
    signal (SIGINT, intcatch);
    (void) umfpack_zi_solve (UMFPACK_Aat, a->ia, a->ja, a->d, MDRPTR(Az),
     MDRPTR(Xx), MDRPTR(Xz), MDRPTR(Bx), MDRPTR(Bz), Numeric, NULL, NULL);
    signal (SIGINT, intcatch_wait);
    // Load result into X
    for (i = 0; i < MNR (b); i++)
    {
      Mdc0(x, i, j) = Mdr0 (Xx, i, 0) + Mdr0 (Xz, i, 0)*I;
    }
  }

  // clean all
  mdr_Destroy (Xx);
  mdr_Destroy (Xz);
  mdr_Destroy (Bx);
  mdr_Destroy (Bz);
  mdr_Destroy (Az);

  umfpack_zi_free_symbolic (&Symbolic);
  umfpack_zi_free_numeric (&Numeric);

  // this is silly but I don;t know better: indexing is offset by one
  for (i = 0; i < a->nr + 1; i++)
    ++a->ia[i];
  for (j = 0; j < a->nnz; j++)
    ++a->ja[j];

  return (x);
}


MDC *
umfpack_msc_SolveCR (MSC * a, MDR * b, char *type)
{
  int i, j, n, ne;
  MDC *x;
  MDR *Xx, *Xz, *Bx, *Bz, *Ax, *Az;
  //
  // Initialize UMFPACK
  //
  double Control[UMFPACK_CONTROL], Info[UMFPACK_INFO];
  void *Symbolic, *Numeric;

  msc_Detect_Inf (a);
  msc_Detect_NaN (a);
  mdr_Detect_Inf (b);
  mdr_Detect_Nan (b);
  if (a->nr != a->nc)
    rerror ("solve[umf]: requires square coefficient matrix");
  if (a->nr != MNR (b))
    rerror ("solve[umf]: RHS row dim. must match LHS row dim.");
  n = a->nr;
  ne = a->nnz;

  // create matrix structures for UMFPACK
  Xx = mdr_Create (n, 1);
  Xz = mdr_Create (n, 1);
  Bx = mdr_Create (n, 1);
  Bz = mdr_Create (n, 1);

  // this is what we pass instead of 'a'
  Ax = mdr_Create (1, ne);
  Az = mdr_Create (1, ne);
  for (i = 0; i < ne; i++)
  {
    MdrV0(Ax,i) = RE(a->c[i]);
    MdrV0(Az,i) = IM(a->c[i]);
  }

  // output
  x = mdc_Create (MNR (b), MNC (b));

  // this is silly but I don;t know better: indexing is offset by one
  for (i = 0; i < a->nr + 1; i++)
    --a->ia[i];
  for (j = 0; j < a->nnz; j++)
    --a->ja[j];

  //
  // Factor the coefficient matrix: first symbolic then numeric
  //
  signal (SIGINT, intcatch);
  (void) umfpack_zi_symbolic (n, n, a->ia, a->ja, MDRPTR(Ax), MDRPTR(Az), &Symbolic,
			      Control, Info);
  signal (SIGINT, intcatch_wait);

  signal (SIGINT, intcatch);
  (void) umfpack_zi_numeric (a->ia, a->ja, MDRPTR(Ax), MDRPTR(Az), Symbolic, &Numeric,
			     Control, Info);
  signal (SIGINT, intcatch_wait);

  //
  // Solve Ax = b
  //
  //
  // Loop over the columns of B, getting a solution for each.
  //
  for (j = 0; j < MNC (b); j++)
  {
    // Grab i-th column of 'b'
    for (i = 0; i < MNR (b); i++)
    {
      Mdr0 (Bx, i, 0) = Mdr0 (b, i, j);
      Mdr0 (Bz, i, 0) = 0;
    }
    //Bx = mdr_PartitionCol(b,j);
    // Solve, but remember A is transposed (MSR is row compressed format)
    signal (SIGINT, intcatch);
    (void) umfpack_zi_solve (UMFPACK_Aat, a->ia, a->ja, MDRPTR(Ax), MDRPTR(Az),
     MDRPTR(Xx), MDRPTR(Xz), MDRPTR(Bx), MDRPTR(Bz), Numeric, NULL, NULL);
    signal (SIGINT, intcatch_wait);
    // Load result into X
    for (i = 0; i < MNR (b); i++)
    {
      Mdc0 (x, i, j) = Mdr0 (Xx, i, 0) + Mdr0 (Xz, i, 0)*I;
    }
    //mdr_Destroy(Bx);
  }

  // clean all
  mdr_Destroy (Xx);
  mdr_Destroy (Xz);
  mdr_Destroy (Bx);
  mdr_Destroy (Bz);
  mdr_Destroy (Ax);
  mdr_Destroy (Az);

  umfpack_zi_free_symbolic (&Symbolic);
  umfpack_zi_free_numeric (&Numeric);
  // this is silly but I don;t know better: indexing is offset by one
  for (i = 0; i < a->nr + 1; i++)
    ++a->ia[i];
  for (j = 0; j < a->nnz; j++)
    ++a->ja[j];

  return (x);
}

#endif

/* **************************************************************
 * Force a sparse real matrix into a column vector format.
 * Retain sparsity.
 * ************************************************************** */

MSC *
msc_ReshapeCol (MSC * m)
{
  int i;
  MSC *col;

  /* First, check for empty matrix. */
  if ((m->nr == 0 && m->nc == 0) || m->nnz == 0)
  {
    col = msc_Create (m->nr, m->nc);
    msc_Setup (col, m->nnz);
  }
  else
  {
    /* Transpose, so we can pull the columns out of the rows. */
    col = msc_NcTranspose (m);

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

#ifdef HAVE_METIS2
/* **************************************************************
 * Fill reducing ordering of a sparse matrix via Metis2 lib
 * ************************************************************** */

extern int OMETIS (int *, int *, int *, int *, int *, int *, int *);

MDR *
msc_SpOrder (MSC * m)
{
  int i, ia, j, n, numbering;
  int options[5], *p, *invp;
  MDR *P = 0;
  MSC *mt, *mt1;

  /* First, check for empty matrix. */
  if ((m->nr == 0 && m->nc == 0) || m->nnz == 0)
  {
    rerror ("sporder: cannot order empty matrix");
  }
  else
  {
    /* First, transpose the matrix so the ia and ja are correct. */
    mt = msc_Transpose (m);

    /*
     * Now, make the sparse matrix "sensible" to Metis.
     * There should be no diagonal elements (self-loops).
     * There should be at least one adjancey per row, even
     * if it is a zero.
     */

    /*
     * Remove the diagonal elements...
     * And, check for empty rows.
     */

    for (i = 0; i < mt->nr; i++)
    {
      n = mt->ia[i + 1] - mt->ia[i];
      ia = mt->ia[i];

      if (n == 0)
      {
	msc_Destroy (mt);
	rerror ("sporder: input matrix has an empty column");
      }

      for (j = 0; j < n; j++)
      {
	if (mt->ja[ia + j - 1] == i + 1)
	{
	  /* Found a diagonal element, set to zero. */
	  mt->c[ia + j - 1] = 0.0;
	}
      }
    }

    mt1 = msc_ReSparse (mt);
    msc_Destroy (mt);
    mt = mt1;

    n = mt->nr;

    /* ometis options array. */
    options[0] = 0;

    /* We are using Fortran-style array indexing. */
    numbering = 1;

    /* Permutation arrays. */
    p = (int *) GC_MAIOP (n * sizeof (int));
    invp = (int *) GC_MAIOP (n * sizeof (int));

    OMETIS (&n, mt->ia, mt->ja, options, &numbering, p, invp);

    /* Now, load P. */
    P = mdr_Create (1, n);
    for (i = 0; i < n; i++)
    {
      MdrV0 (P, i) = (double) invp[i];
    }

    /* Clean up. */
    msc_Destroy (mt);
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

typedef int idxtype;
extern void METIS_EdgeND (int *, idxtype *, idxtype *, int *, int *, idxtype *,
			  idxtype *);
extern void METIS_NodeND (int *, idxtype *, idxtype *, int *, int *, idxtype *,
			  idxtype *);

MDR *
msc_SpOrder (MSC * m)
{
  int i, ia, j, n, numbering;
  int options[8];
  idxtype *p, *invp;
  MDR *P = 0;
  MSC *mt, *mt1;

  /* First, check for empty matrix. */
  if ((m->nr == 0 && m->nc == 0) || m->nnz == 0)
  {
    rerror ("sporder: cannot order empty matrix");
  }
  else
  {
    /* First, transpose the matrix so the ia and ja are correct. */
    mt = msc_Transpose (m);

    /*
     * Now, make the sparse matrix "sensible" to Metis.
     * There should be no diagonal elements (self-loops).
     * There should be at least one adjancey per row, even
     * if it is a zero.
     */

    /*
     * Remove the diagonal elements...
     * And, check for empty rows.
     */

    for (i = 0; i < mt->nr; i++)
    {
      n = mt->ia[i + 1] - mt->ia[i];
      ia = mt->ia[i];

      if (n == 0)
      {
	msc_Destroy (mt);
	rerror ("sporder: input matrix has an empty column");
      }

      for (j = 0; j < n; j++)
      {
	if (mt->ja[ia + j - 1] == i + 1)
	{
	  /* Found a diagonal element, set to zero. */
	  mt->c[ia + j - 1] = 0.0;
	}
      }
    }

    mt1 = msc_ReSparse (mt);
    msc_Destroy (mt);
    mt = mt1;

    n = mt->nr;

    /* ometis options array. */
    options[0] = 0;

    /* We are using Fortran-style array indexing. */
    numbering = 1;

    /* Permutation arrays. */
    p = (int *) GC_MAIOP (n * sizeof (int));
    invp = (int *) GC_MAIOP (n * sizeof (int));

    METIS_NodeND (&n, mt->ia, mt->ja, &numbering, options, p, invp);

    /* Now, load P. */
    P = mdr_Create (1, n);
    for (i = 0; i < n; i++)
    {
      MdrV0 (P, i) = (double) invp[i];
    }

    /* Clean up. */
    msc_Destroy (mt);
    GC_FREE (p);
    GC_FREE (invp);
  }
  return (P);
}
#endif /* HAVE_METIS3 */


#ifdef HAVE_SUPERLU
/* SuperLU interface
   Tzong-Shuoh Yang (yang@isec.com) 10/3/96 solve()
                                    1/14/97 factor(), backsub()
   Modified for Complex sparse matrices and SuperLU-1.0.
   Ian Searle, 2/9/97.
   Modified for Complex sparse matrices and SuperLU-3.0
   Marijan Kostrun, IV/2005
*/

//
// superlu parameters
//
extern double superlu_diagpivothresh;
extern int superlu_permc_spec;
extern int superlu_n;
extern double superlu_droptol;


/* Do this because of possible conflicts on some platforms (Solaris). */
#undef _S
#undef _D
#undef _C
#undef _Z

// #include "./clibs/superlu/src/zsp_defs.h"
// #include "./clibs/superlu/src/dcomplex.h"
#include <superlu/superlu_enum_consts.h>
#include <superlu/slu_ddefs.h>
#define  Matrix SuperMatrix

extern void Destroy_CompCol_Matrix (Matrix *);
extern void Destroy_SuperNode_Matrix (Matrix *);
extern void Destroy_CompCol_Permuted (Matrix *);
extern void StatInit (SuperLUStat_t *);	//(int,int);
extern void StatFree ();
extern void sp_preorder (superlu_options_t *, Matrix *, int *, int *, Matrix *);
extern int sp_ienv (int);
extern void zgstrf (superlu_options_t *, SuperMatrix *, double,
		    int, int, int *, void *, int, int *, int *,
		    SuperMatrix *, SuperMatrix *, SuperLUStat_t *, int *);
extern void zgstrs (trans_t, SuperMatrix *, SuperMatrix *, int *, int *,
		    SuperMatrix *, SuperLUStat_t *, int *);
extern void zgscon (char *, SuperMatrix *, SuperMatrix *,
		    double, double *, SuperLUStat_t *, int *);
extern double zlangs (char *, Matrix *);
extern void get_perm_c (int, Matrix *, int *);
extern void zCreate_CompCol_Matrix(SuperMatrix *, int, int, int, Complex *,
                                   int *, int *, Stype_t, Dtype_t, Mtype_t);
extern void zCreate_CompRow_Matrix(SuperMatrix *, int, int, int, Complex *,
                                   int *, int *, Stype_t, Dtype_t, Mtype_t);
extern void zCreate_Dense_Matrix(SuperMatrix *, int, int, Complex *, int,
                                 Stype_t, Dtype_t, Mtype_t);
extern void zgssv(superlu_options_t *, SuperMatrix *, int *, int *, SuperMatrix *,
                  SuperMatrix *, SuperMatrix *, SuperLUStat_t *, int *);

/* this routine comes from SuperLU/MATLAB/mexsuperlu.c */
static void
msc_LUextract (Matrix * L, Matrix * U, Complex * Lval, int *Lrow, int *Lcol,
	       Complex * Uval, int *Urow, int *Ucol, int *snnzL, int *snnzU)
{
  int i, j, k;
  int upper;
  int fsupc, istart, nsupr;
  int lastl = 0, lastu = 0;
  SCformat *Lstore;
  NCformat *Ustore;
  Complex *SNptr;

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
      SNptr = &((Complex *) Lstore->nzval)[L_NZ_START (j)];

      /* Extract U */
      for (i = U_NZ_START (j); i < U_NZ_START (j + 1); ++i)
      {
	Uval[lastu] = ((Complex *) Ustore->nzval)[i];
	/* Matlab doesn't like explicit zero. */
	if (cabs(Uval[lastu]) != 0.0)
	  Urow[lastu++] = U_SUB (i);
      }
      for (i = 0; i < upper; ++i)
      {				/* upper triangle in the supernode */
	Uval[lastu] = SNptr[i];
	/* Matlab doesn't like explicit zero. */
	if (cabs(Uval[lastu]) != 0.0)
	  Urow[lastu++] = L_SUB (istart + i);
      }
      Ucol[j + 1] = lastu;

      /* Extract L */
      Lval[lastl] = 1.0;	/* unit diagonal */
      Lrow[lastl++] = L_SUB (istart + upper - 1);
      for (i = upper; i < nsupr; ++i)
      {
	Lval[lastl] = SNptr[i];
	/* Matlab doesn't like explicit zero. */
	if (cabs(Lval[lastl]) != 0.0)
	  Lrow[lastl++] = L_SUB (istart + i);
      }
      Lcol[j + 1] = lastl;

      ++upper;

    }				/* for j ... */

  }				/* for k ... */

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
msc_super_fact (Matrix * A, int *perm_c, int *perm_r, int *etree,
		Matrix * L, Matrix * U, int *info, double *rcond,
		double diag_pivot_thresh)
{
  char norm[1];
  Matrix AC;			/* Matrix postmultiplied by Pc */
  int lwork = 0;
  SuperLUStat_t slu_stat;	// kmk
  superlu_options_t slu_options;
  set_default_options (&slu_options);
  // Set defaults
  double anorm;
  double drop_tol = 0;
  int panel_size;
  int relax;
  //
  // Options
  //
  slu_options.Fact = DOFACT;
  panel_size = sp_ienv (1);
  relax = sp_ienv (2);
  StatInit (&slu_stat);

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
  signal (SIGINT, intcatch);
  zgstrf (&slu_options, &AC, drop_tol, relax, panel_size, etree, NULL, lwork,
	  perm_c, perm_r, L, U, &slu_stat, info);
  signal (SIGINT, intcatch_wait);
  if (*info == 0)
  {
    // estimate the reciprocal of the condition number of A.
    *norm = '1';
    anorm = zlangs (norm, A);
    zgscon (norm, L, U, anorm, rcond, &slu_stat, info);
  }
  // Clean Up.
  Destroy_CompCol_Permuted (&AC);
  StatFree (&slu_stat);
}

// ***************************************************************
// * Solve a set of equations with a sparse coefficient matrix   *
// * as quickly as possible and with the least memory            *
// ***************************************************************

MDC *
superlu_msc_SolveCC (MSC * a, MDC * b, char *type)
{
  int i;
  SuperMatrix A, B, L, U;
  MDC *x = mdc_Copy (b);
  int n, ne, *perm_c, *perm_r, info, *ia, *ja;
  SuperLUStat_t slu_stat;
  superlu_options_t options;
  set_default_options (&options);
  //
  // Check for un-manageable matrix elements.
  //
  msc_Detect_Inf (a);
  msc_Detect_NaN (a);
  mdc_Detect_Inf (b);
  mdc_Detect_Nan (b);
  //
  // Check for poorly posed problem.
  //
  if (a->nr != a->nc)
    rerror ("solve [slu]: requires square coefficient matrix");
  if (a->nr != MNR (b))
    rerror ("solve [slu]: RHS row dim. must match LHS row dim.");
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
  // Setup the SuperLU [A] from [a]: MSC * a is row ordered.
  //
  zCreate_CompRow_Matrix (&A, n, n, ne, a->c, a->ja, a->ia, SLU_NR, SLU_Z,
			  SLU_GE);
  //
  // Options
  //
  options.DiagPivotThresh = superlu_diagpivothresh;
  if (*type == 's' || *type == 'S')
    options.SymmetricMode = YES;
  else
    options.SymmetricMode = NO;

  //
  // Setup the SuperLU RHS from [x] (a copy of [b]).
  //
  zCreate_Dense_Matrix (&B, n, (b->ncol), MDCPTR(x), n, SLU_DN, SLU_Z, SLU_GE);
  //
  // Allocate permutation arrays
  //
  perm_r = (int *) GC_MAIOP (n * sizeof (int));
  perm_c = (int *) GC_MAIOP (n * sizeof (int));
  //
  // Solve
  //
  StatInit (&slu_stat);
  signal (SIGINT, intcatch);
  zgssv (&options, &A, perm_c, perm_r, &L, &U, &B, &slu_stat, &info);
  signal (SIGINT, intcatch_wait);

  StatFree (&slu_stat);
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
  Destroy_CompCol_Matrix (&U);
  Destroy_SuperNode_Matrix (&L);
  GC_FREE (perm_c);
  GC_FREE (perm_r);
  GC_FREE (A.Store);
  GC_FREE (B.Store);
  return (x);
}


MDC *
superlu_msc_SolveCR (MSC * a, MDR * b, char *type)
{
  int i;
  SuperMatrix A, B, L, U;
  MDC *x = mdr_coerce_mdc (b);
  int n, ne, *perm_c, *perm_r, info, *ia, *ja;
  SuperLUStat_t slu_stat;
  superlu_options_t options;
  set_default_options (&options);
  //
  // Check for un-manageable matrix elements.
  //
  msc_Detect_Inf (a);
  msc_Detect_NaN (a);
  mdr_Detect_Inf (b);
  mdr_Detect_Nan (b);
  //
  // Check for poorly posed problem.
  //
  if (a->nr != a->nc)
    rerror ("solve [slu]: requires square coefficient matrix");
  if (a->nr != MNR (b))
    rerror ("solve [slu]: RHS row dim. must match LHS row dim.");
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
  // Setup the SuperLU [A] from [a]: MSC * a is row ordered.
  //
  zCreate_CompRow_Matrix (&A, n, n, ne, a->c, a->ja, a->ia, SLU_NR, SLU_Z,
			  SLU_GE);
  //
  // Setup the SuperLU RHS from [x] (a copy of [b]).
  //
  zCreate_Dense_Matrix (&B, n, MNC (b), MDCPTR(x), n, SLU_DN, SLU_Z, SLU_GE);
  //
  // Allocate permutation arrays
  //
  perm_r = (int *) GC_MAIOP (n * sizeof (int));
  perm_c = (int *) GC_MAIOP (n * sizeof (int));
  //
  // Solve
  //
  if (*type == 's' || *type == 'S')
    options.SymmetricMode = YES;
  else
    options.SymmetricMode = NO;
  StatInit (&slu_stat);
  signal (SIGINT, intcatch);
  zgssv (&options, &A, perm_c, perm_r, &L, &U, &B, &slu_stat, &info);
  signal (SIGINT, intcatch_wait);

  StatFree (&slu_stat);
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
  Destroy_CompCol_Matrix (&U);
  Destroy_SuperNode_Matrix (&L);
  GC_FREE (perm_c);
  GC_FREE (perm_r);
  GC_FREE (A.Store);
  GC_FREE (B.Store);
  return (x);
}


MDC *
msc_SpSolve (MSC * a, MDC * b, double diag_pivot_thresh, int *perm_c)
{
  int i;
  Matrix A, B;
  Matrix *L, *U;
//   NCformat *Astore;
//   DNformat *Bstore;
  MDC *x;
  int *etree, n, ne, *perm_ctmp, *perm_r, info, *ia, *ja;
  double rcond;
  trans_t trans;
  SuperLUStat_t slu_stat;
  //
  // Check for un-manageable matrix elements.
  //
  msc_Detect_Inf (a);
  msc_Detect_NaN (a);
  mdc_Detect_Inf (b);
  mdc_Detect_Nan (b);
  //
  // Check for poorly posed problem.
  //
  if (a->nr != a->nc)
    rerror ("solve: requires square coefficient matrix");
  if (a->nr != MNR (b))
    rerror ("solve: RHS row dim. must match LHS row dim.");
  n = a->nr;
  ne = a->nnz;
  x = mdc_Copy (b);
  perm_r = (int *) GC_MAIOP (n * sizeof (int));
  //
  // Setup the SuperLU [A] from [a]: MSR * a is row ordered. Matrix *A
  // is column ordered. Thus A = a' and in the end we will have to use
  // transpose flag for backsubstitution.
  // kmk, IV/2005
  //
  //A.Stype = SLU_NC;
  //A.Dtype = SLU_Z;
  //A.Mtype = SLU_GE;
  //A.nrow  = n;
  //A.ncol  = n;
  //A.Store = GC_MALLOC(sizeof(NCformat));
  //Astore  = (NCformat *) A.Store;
  //Astore->nnz    = ne;
  //Astore->nzval  = a->c;
  //Astore->rowind = a->ja;
  //Astore->colptr = a->ia;
  zCreate_CompCol_Matrix (&A, n, n, ne, a->c, a->ja, a->ia, SLU_NC, SLU_Z,
			  SLU_GE);
  etree = (int *) GC_MAIOP (A.ncol * sizeof (int));
  //
  // Setup the SuperLU RHS from [x] (a copy of [b]).
  //
  //B.Stype = SLU_DN;
  //B.Dtype = SLU_Z;
  //B.Mtype = SLU_GE;
  //B.nrow  = n;
  //B.ncol  = MNC (b);
  //B.Store = (DNformat *) GC_MALLOC (sizeof(DNformat));
  //Bstore  = (DNformat *) B.Store;
  //Bstore->lda = n;
  //Bstore->nzval = x->c;
  zCreate_Dense_Matrix (&B, n, MNC (b), MDCPTR(x), n, SLU_DN, SLU_Z, SLU_GE);
  //
  // Allocate matrices L and U. kmk
  //
  L = (Matrix *) GC_MALLOC (sizeof (Matrix));
  U = (Matrix *) GC_MALLOC (sizeof (Matrix));
  // Solve Ax = b using SuperLU
  ia = a->ia;
  ja = a->ja;
  for (i = 0; i < ne; i++)
    ja[i]--;
  for (i = 0; i <= n; i++)
    ia[i]--;
  // create a permutation vector if user does not have one
  if (perm_c == 0)
  {
    int permc_spec = 1;		/* use A'*A */
    perm_ctmp = (int *) GC_MAIOP (n * sizeof (int));
    // use Liu's mmd algorithm
    get_perm_c (permc_spec, &A, perm_ctmp);
  }
  else
  {
    // use user's permutation vector
    perm_ctmp = perm_c;
  }
  // LU decomposition
  signal (SIGINT, intcatch);
  msc_super_fact (&A, perm_ctmp, perm_r, etree, L, U,
		  &info, &rcond, diag_pivot_thresh);
  signal (SIGINT, intcatch_wait);

  if (info != 0)
  {
    if (info <= n)
      fprintf (stderr, "Diagonal zero, equation %d\n", info);
    else
      fprintf (stderr, "Memory allocation failure\n");
    rerror ("SuperLU error");
  }
  // Backsub
  {
//     int panel_size, relax;
//     panel_size = sp_ienv (1);
//     relax = sp_ienv (2);
    //StatInit (panel_size, relax);
    StatInit (&slu_stat);
  }
  //
  // From above: SuperLU uses column-wise storage so Matrix *A = (MSR *a)'
  //
  trans = TRANS;

  signal (SIGINT, intcatch);
  zgstrs (trans, L, U, perm_ctmp, perm_r, &B, &slu_stat, &info);
  signal (SIGINT, intcatch_wait);

  StatFree (&slu_stat);
  if (info != 0)
  {
    fprintf (stderr, " error in backsub, info = %d\n", info);
    rerror ("SuperLU error");
  }
  // index adjusting
  for (i = 0; i < ne; i++)
    ja[i]++;
  for (i = 0; i <= n; i++)
    ia[i]++;
  // SuperLU created these, so SuperLU must destroy them.
  Destroy_CompCol_Matrix (U);
  Destroy_SuperNode_Matrix (L);
  if (perm_c == 0)
    GC_FREE (perm_ctmp);
  GC_FREE (perm_r);
  GC_FREE (etree);
  GC_FREE (A.Store);
  GC_FREE (B.Store);
  return (x);
}

/* **************************************************************
 * Solve a set of equations with a sparse coefficient matrix
 * All we do here is solve the general sparse matrix.
 * Use the SuperLU routines.
 * ************************************************************** */

MDC *
msc_Solve (MSC * a, MDC * b, char *type)
{
  double diag_pivot_thresh = 1.0;
  int *perm_c=0;
  MDC *x;
  x = msc_SpSolve (a, b, diag_pivot_thresh, perm_c);
  return (x);
}

/* **************************************************************
 * Factor a sparse real matrix.
 * This function is an interface to SuperLU routines.
 * The return value is a list - n, l, u, perm_r, perm_c and rcond.
 * See msr_Backsub() for the meaning of the return list.
 * ************************************************************** */

Btree *
msc_SpFactor (MSC * a)
{
  Btree *bt;
  int i;
  Matrix A, L, U;
  NCformat *Astore, *Ustore;
  SCformat *Lstore;
  MSC *acopy;
  int *etree, n, ne, *perm_ctmp, *perm_r, info, *ia, *ja;
  double rcond;
  Ent *N, *PERM_C, *PERM_R, *RCOND;
  MDR *mperm_r, *mperm_c;

  /* Create the list that will contain the results. */
  bt = btree_Create ();

  /* Check data */
  msc_Detect_Inf (a);
  msc_Detect_NaN (a);

  if (a->nr != a->nc)
    rerror ("factor: requires square coefficient matrix");

  n = a->nr;
  ne = a->nnz;

  /* superLU uses column-wise storage */
  acopy = msc_NcTranspose (a);

  perm_ctmp = (int *) GC_MAIOP (n * sizeof (int));
  perm_r = (int *) GC_MAIOP (n * sizeof (int));

  A.Stype = SLU_NC;
  A.Dtype = SLU_Z;
  A.Mtype = SLU_GE;
  A.nrow = n;
  A.ncol = n;
  A.Store = GC_MALLOC (sizeof (NCformat));
  Astore = (NCformat *) A.Store;
  Astore->nnz = ne;
  Astore->nzval = acopy->c;
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
  msc_super_fact (&A, perm_ctmp, perm_r, etree, &L, &U, &info,
		  &rcond, superlu_diagpivothresh);
  signal (SIGINT, intcatch_wait);
  if (info != 0)
  {
    if (info <= n)
      fprintf (stderr, "Diagonal zero, equation %d\n", info);
    else
      fprintf (stderr, "Memory allocation failure\n");
    rerror ("SuperLU error");
  }
  msc_Destroy (acopy);
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
    MSC *rL, *rU, *rLt, *rUt;
    int nnzL, nnzU, *ia, *ja;

    Lstore = (SCformat *) L.Store;
    rL = msc_Create (n, n);
    nnzL = Lstore->nnz;
    rL->c = (Complex *) GC_MAIOP (nnzL * sizeof (Complex));
    rL->ja = (int *) GC_MAIOP (nnzL * sizeof (int));
    rL->ia = (int *) GC_MAIOP ((n + 1) * sizeof (int));

    Ustore = (NCformat *) U.Store;
    rU = msc_Create (n, n);
    nnzU = Ustore->nnz;
    rU->c = (Complex *) GC_MAIOP (nnzU * sizeof (Complex));
    rU->ja = (int *) GC_MAIOP (nnzU * sizeof (int));
    rU->ia = (int *) GC_MAIOP ((n + 1) * sizeof (int));
    msc_LUextract (&L, &U, rL->c, rL->ja, rL->ia,
		   rU->c, rU->ja, rU->ia, &nnzL, &nnzU);
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
    rLt = msc_NcTranspose (rL);
    msc_Destroy (rL);
    l = ent_Create ();
    ent_SetType (l, MATRIX_SPARSE_COMPLEX);
    ent_data (l) = rLt;
    install (bt, ("l"), l);
    rUt = msc_NcTranspose (rU);
    msc_Destroy (rU);
    u = ent_Create ();
    ent_SetType (u, MATRIX_SPARSE_COMPLEX);
    ent_data (u) = rUt;
    install (bt, ("u"), u);
  }

  /* Load perm_r */
  mperm_r = mdr_Create (1, n);
  for (i = 0; i < n; i++)
  {
    MdrV0 (mperm_r, i) = (double) (perm_r[i] + 1);	/* starts from 1 */
  }

  GC_FREE (perm_r);
  PERM_R = ent_Create ();
  ent_SetType (PERM_R, MATRIX_DENSE_REAL);
  ent_data (PERM_R) = mperm_r;
  install (bt, ("perm_r"), PERM_R);

  /* Load perm_c */
  mperm_c = mdr_Create (1, n);
  for (i = 0; i < n; i++)
  {
    MdrV0 (mperm_c, i) = (double) (perm_ctmp[i] + 1);
  }

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
 * Factor
 * ************************************************************** */

Btree *
msc_Factor (MSC * a, char *type)
{
  Btree *x;

  x = msc_SpFactor (a);

  return (x);
}
#endif


MDC *
msc_Backsub (Btree * bt, MDC * b)
{
  MSC *L = 0, *U = 0;
  MDC *sol;
  ListNode *ltmp;
  Ent *etmp;
  int i, j, kk, n = 0, m, *perm_c = 0, *perm_r = 0;
  MDR *tmp;
  Complex *x, *y, *z, *r;

  mdc_Detect_Inf (b);
  mdc_Detect_Nan (b);

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
    if (ent_type (etmp) != MATRIX_SPARSE_COMPLEX)
      rerror ("backsub: terrible input list error");
    L = (MSC *) ent_data (etmp);
  }
  else
    rerror ("backsub: terrible input list error");

  if ((ltmp = btree_FindNode (bt, "u")))
  {
    etmp = var_ent (ltmp);
    if (ent_type (etmp) != MATRIX_SPARSE_COMPLEX)
      rerror ("backsub: terrible input list error");
    U = (MSC *) ent_data (etmp);
  }
  else
    rerror ("backsub: terrible input list error");

  sol = mdc_Copy (b);		/* the solution */
  r = (Complex *) GC_MAIOP (n * sizeof (Complex));	/* work array */
  z = r;

  /* Backsub */

  signal (SIGINT, intcatch);
  for (kk = 0; kk < m; kk++)	/* loop thru no. of right-hand-side vectors */
  {
    Complex ctmp, t;
    int nnzU;

    x = &MdcV0(sol,kk * n);
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
        ctmp = L->c[j - 1] * y[L->ja[j - 1] - 1];
        t -= ctmp;
      }
      y[i] = t;
    }
    /* solve U*z = y */
    nnzU = U->nnz;
    z[n - 1] = y[n - 1] / U->c[nnzU - 1];
    for (i = n - 2; i >= 0; i--)
    {
      t = y[i];
      for (j = U->ia[i]; j < (U->ia[i + 1] - 1); j++)
      {
        ctmp = U->c[j] * z[U->ja[j] - 1];
        t -= ctmp;
      }
      z[i] = t / U->c[U->ia[i] - 1];
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


Ent *
ent_sparse_compsolv (int nargs, Datum args[])
{
  // call parameters:
  //  e1   - rng name
  //  [e2] - seed
  Ent *e1 = 0, *rent;
  MDS *x1=0;
  char sparam[32] = { '\0' };
  if (nargs == 0)
  {
    fprintf (stdout, "compsolv: A member of the function list  sparams.\n");
    fprintf (stdout,
	     "compsolv: Chooses a method for solving sparse complex linear systems.\n");
    fprintf (stdout, "compsolv: Format:\n");
    fprintf (stdout, "compsolv:   sparams.compsolv(\"method\" ),\n");
    fprintf (stdout,
	     "compsolv: Available methods/packages are:\"UMFPACK\" (\"umfpack\",\"umf\")\n");
    fprintf (stdout, "compsolv: and \"SuperLU\" (\"superlu\",\"slu\").\n");
    fprintf (stdout, "compsolv: Currently the method ");
    if (msc_Solve_method == 0)
      fprintf (stdout, "UMFPACK");
    else
      fprintf (stdout, "SuperLU");
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

  if (strcmp (sparam, "umfpack") == 0)
    msc_Solve_method = 0;
  else if (strcmp (sparam, "UMFPACK") == 0)
    msc_Solve_method = 0;
  else if (strcmp (sparam, "umf") == 0)
    msc_Solve_method = 0;
#ifdef HAVE_SUPERLU
  else if (strcmp (sparam, "SuperLU") == 0)
    msc_Solve_method = 1;
  else if (strcmp (sparam, "superlu") == 0)
    msc_Solve_method = 1;
  else if (strcmp (sparam, "SUPERLU") == 0)
    msc_Solve_method = 1;
  else if (strcmp (sparam, "slu") == 0)
    msc_Solve_method = 1;
#endif
  else
    msc_Solve_method = 0;

  rent = ent_Create ();
  ent_data (rent) = mdr_CreateScalar (1);
  ent_SetType (rent, MATRIX_DENSE_REAL);
  return rent;
}

MDC *
msc_SolveCC (MSC * a, MDC * b, char *type)
{
  void *(*vfptr) ();
#ifdef HAVE_SUPERLU
  if (msc_Solve_method == 1)
    vfptr = (VFPTR) superlu_msc_SolveCC;
  else
#endif
#ifdef HAVE_SUITESPARSE
    vfptr = (VFPTR) umfpack_msc_SolveCC;
#endif

  return (*vfptr) (a, b, type);
}

MDC *
msc_SolveCR (MSC * a, MDR * b, char *type)
{
  void *(*vfptr) ();
#ifdef HAVE_SUPERLU
  if (msc_Solve_method == 1)
    vfptr = (VFPTR) superlu_msc_SolveCR;
  else
#endif
#ifdef HAVE_SUITESPARSE
    vfptr = (VFPTR) umfpack_msc_SolveCR;
#endif
  return (*vfptr) (a, b, type);
}
