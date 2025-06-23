/* mscf1.c Matrix Sparse Complex Functions ... */

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
#include "mem.h"
#include "msr.h"
#include "msrf1.h"
#include "msc.h"
#include "mdr.h"
#include "mdr_mdc.h"
#include "mdc.h"
#include "mds.h"
#include "btree.h"
#include "util.h"
#include "symbol.h"
#include "sort.h"

#include "mathl.h"
#include "fi.h"
#include "sparse.h"

#include <stdio.h>
#include <math.h>

/* **************************************************************
 * Convert a dense complex matrix to a sparse complex matrix.
 * ************************************************************** */

MSC *
msc_Sparse (MDC * m)
{
  MSC *sm;
  int i, j, k, nnz, nc, nr, size;

  nnz = 0;
  nr = m->nrow;
  nc = m->ncol;
  sm = msc_Create (nr, nc);

  /* Count the number of non-zeros */
  size = nr * nc;
  for (i = 0; i < size; i++)
  {
    if ((MdcV0r (m, i) != 0.0) || (MdcV0i (m, i) != 0.0))
    {
      nnz++;
    }
  }
  sm->nnz = nnz;

  /* Now allocate the sparse storage arrays. */
  sm->c = (Complex *) GC_MAIOP (nnz * sizeof (Complex));
  sm->ja = (int *) GC_MAIOP (nnz * sizeof (int));
  sm->ia = (int *) GC_MAIOP ((nr + 1) * sizeof (int));
  sm->ia[0] = 1;
  k = 1;

  /* Fill the sparse storage arrays with the proper values. */
  for (i = 0; i < nr; i++)
  {
    for (j = 0; j < nc; j++)
    {
      if ((Mdc0r (m, i, j) != 0.0) || (Mdc0i (m, i, j) != 0.0))
      {
        sm->c[k - 1] = Mdc0 (m, i, j);
        sm->ja[k - 1] = j + 1;
        k++;
      }
    }
    sm->ia[i + 1] = k;
  }

  sm->order = 1;
  return (sm);
}

/* **************************************************************
 * Convert a sparse complex matrix to a dense real matrix.
 * ************************************************************** */

MDC *
msc_Dense (MSC * m)
{
  MDC *dm;
  int i, ia, j, n, nr, nc, row;

  nr = m->nr;
  nc = m->nc;
  dm = mdc_Create (nr, nc);
  mdc_Zero (dm);

  if (m->nnz)
  {
    row = 1;
    for (i = 0; i < nr; i++)
    {
      n = m->ia[i + 1] - m->ia[i];
      ia = m->ia[i];
      for (j = 0; j < n; j++)
      {
        Mdc1(dm, row, m->ja[ia + j - 1]) = m->c[ia + j - 1];
      }
      row++;
    }
  }
  return (dm);
}

/* **************************************************************
 * Convert a sparse matrix to dense-4-column representation.
 * 1st column:  Row id
 * 2nd column:  Column id
 * 3rd column:  matrix elements (non-zeros, real)
 * 4th column:  matrix elements (non-zeros, imaginary)
 * ************************************************************** */

MDR *
msc_Spconvert (MSC * m)
{
  MDR *dm;
  int i, ii, ia, j, n, row;

  if (((m->nr == 0) && (m->nc == 0)) || (m->nnz == 0))
  {
    dm = mdr_Create (0, 0);
    return (dm);
  }

  dm = mdr_Create (m->nnz, 4);

  row = 1;
  ii = 1;
  for (i = 0; i < m->nr; i++)
  {
    n = m->ia[i + 1] - m->ia[i];
    ia = m->ia[i];

    for (j = 0; j < n; j++)
    {
      Mdr1 (dm, ii, 1) = row;
      Mdr1 (dm, ii, 2) = m->ja[ia + j - 1];
      Mdr1 (dm, ii, 3) = RE(m->c[ia + j - 1]);
      Mdr1 (dm, ii, 4) = IM(m->c[ia + j - 1]);
      ii++;
    }
    row++;
  }
  return (dm);
}

/* **************************************************************
 * Append: [ Sparse , Dense ]
 * ************************************************************** */

MDC *
msc_mdc_Append (MSC * m1, MDC * m2)
{
  MDC *m, *mtmp;

  mtmp = msc_Dense (m1);
  m = mdc_Append (mtmp, m2);
  mdc_Destroy (mtmp);

  return (m);
}

/* **************************************************************
 * Append: [ Dense , Sparse ]
 * ************************************************************** */

MDC *
mdc_msc_Append (MDC * m1, MSC * m2)
{
  MDC *m, *mtmp;

  mtmp = msc_Dense (m2);
  m = mdc_Append (m1, mtmp);
  mdc_Destroy (mtmp);

  return (m);
}

/* **************************************************************
 * Stack: [ Sparse ; Dense ]
 * ************************************************************** */

MDC *
msc_mdc_Stack (MSC * m1, MDC * m2)
{
  MDC *m, *mtmp;

  mtmp = msc_Dense (m1);
  m = mdc_Stack (mtmp, m2);
  mdc_Destroy (mtmp);

  return (m);
}

/* **************************************************************
 * Stack: [ Dense ; Sparse ]
 * ************************************************************** */

MDC *
mdc_msc_Stack (MDC * m1, MSC * m2)
{
  MDC *m, *mtmp;

  mtmp = msc_Dense (m2);
  m = mdc_Stack (m1, mtmp);
  mdc_Destroy (mtmp);

  return (m);
}

/* **************************************************************
 * Add two sparse matrices.
 * ************************************************************** */

MSC *
msc_Add (MSC * m1, MSC * m2)
{
  MSC *new;
  Complex *w;
  int *cia, *cja, *iw;

  /* Some special cases. */
  if (((m1->nr == 0) && (m1->nc == 0)) || ((m2->nr == 0) && (m2->nc == 0)))
  {
    new = msc_Create (0, 0);
    return (new);
  }

  /* For now, just add two matrices of the same dimension! */
  if (m1->nr != m2->nr)
  {
    fprintf (stderr, "matrix row sizes must be equal\n");
    fprintf (stderr, "matrix row sizes: %i and %i\n", m1->nr, m2->nr);
    rerror ("sparse-add: matrices must have the same row dimension");
  }

  if (m1->nc != m2->nc)
  {
    fprintf (stderr, "matrix column sizes must be equal\n");
    fprintf (stderr, "matrix column sizes: %i and %i\n", m1->nc, m2->nc);
    rerror ("sparse-add: matrices must have the same column dimension");
  }

  /*
   * Now get on with the add...
   */

  /* Take care of special cases... */
  if (m1->nnz == 0 && m2->nnz == 0)
  {
    new = msc_Create (m1->nr, m1->nc);
    return (new);
  }

  if (m1->nnz == 0)
  {
    return (msc_Copy (m2));
  }

  if (m2->nnz == 0)
  {
    return (msc_Copy (m1));
  }

  /* Create the new matrix. */
  new = msc_Create (m1->nr, m1->nc);

  /* Create the working space for the symbolic addition. */
  cia = (int *) GC_MAIOP ((m1->nr + 1) * sizeof (int));
  cja = (int *) GC_MAIOP ((m1->nnz + m2->nnz) * sizeof (int));
  iw = (int *) GC_MAIOP ((m1->nc + 1) * sizeof (int));

  /* Do the symbolic part of the add. */
  XGSADD (m1->ia, m1->ja, m2->ia, m2->ja, &m1->nr, &m1->nc, cia, cja, iw);

  /* Start forming the new matrix. */
  new->nnz = cia[m1->nr] - 1;
  new->ia = cia;

  /*
   * Allocate new space, and copy the data. new->nnz is probably less
   * than (m1->nnz + m2->nnz). Keep cja around for the transpose/order
   * operations.
   */

  new->ja = (int *) GC_MAIOP ((new->nnz) * sizeof (int));
  memcpy (new->ja, cja, new->nnz * sizeof (int));

  /* Transpose twice to order the new structure. */
  XGSTRN (new->ia, new->ja, &new->nr, &new->nc, iw, cja);
  XGSTRN (iw, cja, &new->nc, &new->nr, new->ia, new->ja);
  GC_FREE (cja);
  GC_FREE (iw);

  /* Make room for the elements. */
  new->c = (Complex *) GC_MAIOP (new->nnz * sizeof (Complex));

  /* Now do the add. */
  w = (Complex *) GC_MAIOP (new->nc * sizeof (Complex));
  ZGSADD (m1->ia, m1->ja, m1->c,
	  m2->ia, m2->ja, m2->c, &m1->nr, &m1->nc, new->ia, new->ja, new->c, w);

  GC_FREE (w);
  return (new);
}

/* **************************************************************
 * Add a dense complex matrix, and a sparse complex matrix.
 * ************************************************************** */

MDC *
mdc_msc_Add (MDC * m1, MSC * m2)
{
  MDC *m, *mtmp;

  mtmp = msc_Dense (m2);
  m = mdc_Add (m1, mtmp);
  mdc_Destroy (mtmp);

  return (m);
}

/* **************************************************************
 * Add a sparse complex matrix, and a dense complex matrix.
 * ************************************************************** */

MDC *
msc_mdc_Add (MSC * m1, MDC * m2)
{
  MDC *m, *mtmp;

  mtmp = msc_Dense (m1);
  m = mdc_Add (mtmp, m2);
  mdc_Destroy (mtmp);

  return (m);
}

/* **************************************************************
 * Subtract two sparse matrices.
 * ************************************************************** */

MSC *
msc_Subtract (MSC * m1, MSC * m2)
{
  MSC *new;
  Complex *w;
  int *cia, *cja, *iw;

  /* For now, just add two matrices of the same dimension! */
  if (m1->nr != m2->nr)
  {
    fprintf (stderr, "matrix row sizes must be equal\n");
    fprintf (stderr, "matrix row sizes: %i and %i\n", m1->nr, m2->nr);
    rerror ("sparse-subtract: matrices must have the same row dimension");
  }

  if (m1->nc != m2->nc)
  {
    fprintf (stderr, "matrix column sizes must be equal\n");
    fprintf (stderr, "matrix column sizes: %i and %i\n", m1->nc, m2->nc);
    rerror ("sparse-subtract: matrices must have the same column dimension");
  }

  /*
   * Now get on with the subtraction...
   */

  /* Take care of special cases... */
  if (m1->nnz == 0 && m2->nnz == 0)
  {
    new = msc_Create (m1->nr, m1->nc);
    return (new);
  }

  if (m1->nnz == 0)
  {
    return (msc_Negate (m2));
  }

  if (m2->nnz == 0)
  {
    return (msc_Copy (m1));
  }

  /* Create the new matrix. */
  new = msc_Create (m1->nr, m1->nc);

  /* Create the working space for the symbolic addition. */
  cia = (int *) GC_MAIOP ((m1->nr + 1) * sizeof (int));
  cja = (int *) GC_MAIOP ((m1->nnz + m2->nnz) * sizeof (int));
  iw = (int *) GC_MAIOP ((m1->nc + 1) * sizeof (int));

  /* Do the symbolic part of the subtract. */
  XGSADD (m1->ia, m1->ja, m2->ia, m2->ja, &m1->nr, &m1->nc, cia, cja, iw);

  /* Start forming the new matrix. */
  new->nnz = cia[m1->nr] - 1;
  new->ia = cia;

  /*
   * Allocate new space, and copy the data. new->nnz is probably less
   * than (m1->nnz + m2->nnz). Keep cja around for the transpose/order
   * operations.
   */

  new->ja = (int *) GC_MAIOP ((new->nnz) * sizeof (int));
  memcpy (new->ja, cja, new->nnz * sizeof (int));

  /* Transpose twice to order the new structure. */
  XGSTRN (new->ia, new->ja, &new->nr, &new->nc, iw, cja);
  XGSTRN (iw, cja, &new->nc, &new->nr, new->ia, new->ja);
  GC_FREE (cja);
  GC_FREE (iw);

  /* Make room for the elements. */
  new->c = (Complex *) GC_MAIOP (new->nnz * sizeof (Complex));

  /* Now do the subtract. */
  w = (Complex *) GC_MAIOP (new->nc * sizeof (Complex));
  ZGSSUB (m1->ia, m1->ja, m1->c,
	  m2->ia, m2->ja, m2->c, &m1->nr, &m1->nc, new->ia, new->ja, new->c, w);

  GC_FREE (w);
  return (new);
}

/* **************************************************************
 * Subtract a dense complex matrix, and a sparse complex matrix.
 * ************************************************************** */

MDC *
mdc_msc_Subtract (MDC * m1, MSC * m2)
{
  MDC *m, *mtmp;

  mtmp = msc_Dense (m2);
  m = mdc_Subtract (m1, mtmp);
  mdc_Destroy (mtmp);

  return (m);
}

/* **************************************************************
 * Subtract a sparse complex matrix, and a dense complex matrix.
 * ************************************************************** */

MDC *
msc_mdc_Subtract (MSC * m1, MDC * m2)
{
  MDC *m, *mtmp;

  mtmp = msc_Dense (m1);
  m = mdc_Subtract (mtmp, m2);
  mdc_Destroy (mtmp);

  return (m);
}

/* **************************************************************
 * Support functions for sparse matrix multiply.
 * ************************************************************** */

MSC *
msc_scalar_mult (MSC * m, Complex s)
{
  int i;
  MSC *new;

  new = msc_Create (m->nr, m->nc);
  if (m->nnz)
  {
    msc_Setup (new, m->nnz);

    memcpy (new->ia, m->ia, (m->nr + 1) * sizeof (int));
    memcpy (new->ja, m->ja, (m->nnz) * sizeof (int));

    for (i = 0; i < m->nnz; i++)
    {
      new->c[i] = s * m->c[i];
    }
  }
  return (new);
}

/* **************************************************************
 * Multiply two sparse matrices.
 * ************************************************************** */

MSC *
msc_Multiply (MSC * m1, MSC * m2, int *rtype)
{
  MSC *new;
  int maxjc;
  int *ja_tmp, *tmp;
  Complex *a_tmp;

  *rtype = MATRIX_SPARSE_COMPLEX;

  /* For now, just multiply two matrices with consistent dimensions! */
  if (m1->nc != m2->nr)
  {
    fprintf (stderr, "matrix dimensions must be consistent\n");
    fprintf (stderr, "matrix col and row sizes: %i and %i\n", m1->nc, m2->nr);
    rerror ("sparse-multiply: matrices must have the same row dimension");
  }

  /*
   * Now get on with the multiply...
   */

  new = msc_Create (m1->nr, m2->nc);
  if (m1->nnz > 0 && m2->nnz > 0)
  {
    /*
     * Do the symbolic part of the multiply.
     * This can be tricky, cause we don't know apriori how
     * much space is needed. So, we must check with XGSMUL
     * afterwards, and be prepared to try again.
     */

    new->ia = (int *) GC_MAIOP ((new->nr + 1) * sizeof (int));
    tmp = (int *) GC_MAIOP ((new->nc + 1) * sizeof (int));
    maxjc = m1->nnz + m2->nnz;
    ja_tmp = 0;

    do
    {
      maxjc *= 2;
      if (ja_tmp)
	GC_FREE (ja_tmp);
      ja_tmp = (int *) GC_MAIOP (maxjc * sizeof (int));
      XGSMUL (m1->ia, m1->ja, m2->ia, m2->ja,
	      &new->nr, &new->nc, new->ia, ja_tmp, &maxjc, tmp);
    }
    while (!new->ia[0]);

    new->nnz = new->ia[new->nr] - 1;
    new->ja = (int *) GC_MAIOP (new->nnz * sizeof (int));
    memcpy (new->ja, ja_tmp, new->nnz * sizeof (int));

    /* Now, order the representation. */
    XGSTRN (new->ia, new->ja, &new->nr, &new->nc, tmp, ja_tmp);
    XGSTRN (tmp, ja_tmp, &new->nc, &new->nr, new->ia, new->ja);
    GC_FREE (tmp);
    GC_FREE (ja_tmp);

    /* Now, finally, get the data storage, and do the multiply. */
    new->c = (Complex *) GC_MAIOP (new->nnz * sizeof (Complex));
    a_tmp = (Complex *) GC_MAIOP ((new->nc + 1) * sizeof (Complex));
    ZGSMUL (m1->ia, m1->ja, m1->c,
	    m2->ia, m2->ja, m2->c,
	    &new->nr, &new->nc, new->ia, new->ja, new->c, a_tmp);
    GC_FREE (a_tmp);
  }

  return (new);
}

MDC *
mdc_msc_Multiply (MDC * l, MSC * r, int *rtype)
{
  int i, j, k;
  MDC *m;

  *rtype = MATRIX_DENSE_COMPLEX;

  /* For now, just multiply two matrices with consistent dimensions! */
  if (MNC (l) != r->nr)
  {
    fprintf (stderr, "matrix dimensions must be consistent\n");
    fprintf (stderr, "matrix col and row sizes: %i and %i\n", MNC (l), r->nr);
    rerror ("sparse-multiply: matrices must have the same row dimension");
  }

  /*
   * Now get on with the multiply...
   */

  m = mdc_Create (MNR (l), r->nc);
  mdc_Zero (m);
  if (r->nnz == 0)
  {
    return (m);
  }

  for (i = 0; i < MNR (m); i++)
  {
    for (j = 0; j < r->nr; j++)
    {
      for (k = r->ia[j]; k < r->ia[j + 1]; k++)
      {
        Mdc0(m,i,r->ja[k - 1] - 1) += Mdc0(l,i,j) * r->c[k - 1];
      }
    }
  }
  return (m);
}

void *
msc_mdc_Multiply (MSC * l, MDC * r, int *rtype)
{
  int i, j, k;
  void *m = 0;

  /* Check sizes ... */
  if (l->nc != MNR (r))
  {
    if ((l->nr == 1) && (l->nc == 1))
    {
      /* The sparse one is scalar, return dense. */
      m = (void *) mdc_Create (MNR (r), MNC (r));
      for (i = 0; i < MNR (r) * MNC (r); i++)
      {
        MdcV0 ((MDC *) m, i) = l->c[0] * MdcV0(r, i);
      }
      *rtype = MATRIX_DENSE_COMPLEX;
    }
    else if ((MNR (r) == 1) && (MNC (r) == 1))
    {
      /* The dense one is scalar, return sparse. */
      m = (void *) msc_scalar_mult (l, Mdc0 (r, 0, 0));
      *rtype = MATRIX_SPARSE_COMPLEX;
    }
    else if ((l->nr == 0) || (l->nc == 0))
    {
      m = (void *) msc_Create (0, 0);
      *rtype = MATRIX_SPARSE_COMPLEX;
    }
    else if ((MNR (r) == 0) || (MNC (r) == 0))
    {
      m = (void *) msc_Create (0, 0);
      *rtype = MATRIX_SPARSE_COMPLEX;
    }
    else
    {
      fprintf (stderr, "matrix dimensions must be consistent\n");
      fprintf (stderr, "matrix col and row sizes: %i and %i\n", l->nc, MNR (r));
      rerror ("sparse-multiply: matrices must have the same row dimension");
    }
  }
  else
  {
    *rtype = MATRIX_DENSE_COMPLEX;
    m = (void *) mdc_Create (l->nr, MNC (r));
    mdc_Zero ((MDC *) m);

    if (l->nnz == 0)
    {
      return (m);
    }

    for (i = 0; i < MNR (m); i++)
    {
      for (j = 0; j < MNC (m); j++)
      {
        for (k = l->ia[i] - 1; k < l->ia[i + 1] - 1; k++)
        {
          Mdc0((MDC *)m,i,j) += l->c[k] * Mdc0(r,l->ja[k] - 1, j);
        }
      }
    }
  }
  return (m);
}

/* **************************************************************
 * Element Right divide for Sparse-Complex
 * ************************************************************** */

MDC *
msc_ElRdivide (MSC * m1, MSC * m2, int *rtype)
{
  MDC *l, *r;
  MDC *new;

  /* First convert both to dense. */
  l = msc_Dense (m1);
  r = msc_Dense (m2);

  new = mdc_ElRdivide (l, r, rtype);
  mdc_Destroy (l);
  mdc_Destroy (r);

  return (new);
}

/* **************************************************************
 * Element Right divide for Sparse-Complex / Dense-Real
 * ************************************************************** */

MSC *
msc_mdr_ElRdivide (MSC * m1, MDR * m2, int *rtype)
{
  int i;
  void *new=0;
  MDC *l=0;

  /*
   * Look for some common operations so that we might
   * retain sparsity.
   */

  if (MNR (m2) == 1 && MNC (m2) == 1)
  {
    *rtype = MATRIX_SPARSE_COMPLEX;
    new = msc_Create (m1->nr, m1->nc);
    ((MSR *) new)->nnz = m1->nnz;

    /* Copy the IA array. */
    if (m1->ia)
    {
      ((MSR *) new)->ia = (int *) GC_MAIOP (sizeof (int) * (m1->nr + 1));
      memcpy (((MSR *) new)->ia, m1->ia, (m1->nr + 1) * sizeof (int));
    }

    /* Copy the JA array. */
    if (m1->ja)
    {
      ((MSR *) new)->ja = (int *) GC_MAIOP (sizeof (int) * (m1->nnz));
      memcpy (((MSR *) new)->ja, m1->ja, m1->nnz * sizeof (int));
    }

    if (m1->c)
    {
      ((MSC *) new)->c = (Complex *) GC_MAIOP (sizeof (Complex) * (m1->nnz));
    }

    for (i = 0; i < m1->nnz; i++)
    {
      ((MSC *) new)->c[i] = m1->c[i] / MdrV0 (m2, 0);
    }
  }
  else
  {
    *rtype = MATRIX_DENSE_COMPLEX;
    /* First convert to dense. */
    l = msc_Dense (m1);

    new = mdc_mdr_ElRdivide (l, m2, rtype);
    mdc_Destroy (l);
  }

  return (new);
}

MSC *
msr_msc_Add (MSR * l, MSC * r)
{
  MSC *new, *msc;

  msc = msr_coerce_msc (l);
  new = msc_Add (msc, r);
  msc_Destroy (msc);

  return (new);
}

MSC *
msc_msr_Add (MSC * l, MSR * r)
{
  MSC *new, *msc;

  msc = msr_coerce_msc (r);
  new = msc_Add (l, msc);
  msc_Destroy (msc);

  return (new);
}

size_t
msc_Sizeof (MSC * m)
{
  int size = m->nnz * sizeof (int) +
    m->nnz * sizeof (Complex) + (m->nr + 1) * sizeof (int);
  return (size_t) (size);
}

/* **************************************************************
 * Get the size of a a matrix for size()...
 * ************************************************************** */

MDR *
msc_Size_BF (MSC * m)
{
  MDR *size = mdr_Create (1, 2);
  Mdr0 (size, 0, 0) = (double) (m->nr);
  Mdr0 (size, 0, 1) = (double) (m->nc);
  return (size);
}

/* **************************************************************
 * Return the type of a matrix...
 * ************************************************************** */

MDS *
msc_Type_BF (MSC * m)
{
  MDS *type = mds_CreateScalar ( "complex" );
  return (type);
}

/* **************************************************************
 * Builtin functions for any()
 * Return TRUE if ANY of the elements of A are non-zero.
 * ************************************************************** */

MDR *
msc_Any (MSC * m)
{
  int i, j, rptr, nrow;
  MDR *new;
  MSC *mtmp = 0;

  if (m->nr == 1)
  {
    /* Vector operation */
    new = mdr_Create (1, 1);
    Mdr0 (new, 0, 0) = 0.0;
    for (i = 0; i < m->nnz; i++)
    {
      if (cabs(m->c[i])>0)
      {
        MdrV0 (new, i) = 1.0;
        break;
      }
    }
  }
  else
  {
    /* Matrix operation */
    new = mdr_Create (1, m->nc);

    if (m->nnz == 0)
    {
      mdr_Zero (new);
      return (new);
    }

    /* Transpose so we can operate on the rows. */
    mtmp = msc_Transpose (m);

    for (i = 0; i < mtmp->nr; i++)
    {
      rptr = mtmp->ia[i];
      nrow = mtmp->ia[i + 1] - rptr;
      Mdr0 (new, 0, i) = 0.0;

      for (j = 0; j < nrow; j++)
      {
        if (cabs(mtmp->c[rptr + j - 1])>0)
        {
          Mdr0 (new, 0, i) = 1.0;
          break;
        }
      }
    }
  }

  if (mtmp)
    msc_Destroy (mtmp);

  return (new);
}

/* **************************************************************
 * Builtin functions for all()
 * Return TRUE if ALL of the elements of A are non-zero.
 * ************************************************************** */

MDR *
msc_All (MSC * m)
{
  int i, j;
  Complex ctmp;
  MDR *new;
  MSC *mtmp = 0;

  if (m->nr == 1)
  {
    /* Vector operation */
    new = mdr_Create (1, 1);
    if (m->nnz == 0)
    {
      Mdr0 (new, 0, 0) = 0.0;
      return (new);
    }

    Mdr0 (new, 0, 0) = 1.0;
    for (i = 1; i <= m->nc; i++)
    {
      ctmp = msc_GetEl (m, 1, i);
      if (cabs(ctmp) == 0.0)
      {
        Mdr0 (new, 0, 0) = 0.0;
        break;
      }
    }
  }
  else
  {
    /* Matrix operation */
    new = mdr_Create (1, m->nc);

    if (m->nnz == 0)
    {
      mdr_Zero (new);
      return (new);
    }

    /* Transpose so we can operate on the rows. */
    mtmp = msc_NcTranspose (m);

    for (i = 1; i <= mtmp->nr; i++)
    {
      Mdr0 (new, 0, i - 1) = 1.0;

      for (j = 1; j <= mtmp->nc; j++)
      {
        ctmp = msc_GetEl (mtmp, i, j);
        if (cabs(ctmp) == 0.0)
        {
          Mdr0 (new, 0, i - 1) = 0.0;
          break;
        }
      }
    }
  }

  if (mtmp)
    msc_Destroy (mtmp);

  return (new);
}

/* **************************************************************
 * Return the complex-conjugate of a a sparse matrix.
 * ************************************************************** */

MSC *
    msc_Conj (MSC * m)
{
  MSC *new;
  int i;

  new = msc_Copy (m);
  for (i = 0; i < new->nnz; i++)
  {
    IM(new->c[i]) = -IM(new->c[i]);
  }

  return (new);
}

/* **************************************************************
 * Return the "length" of a sparse complex matrix
 * ************************************************************** */

MDR *
msc_Length_BF (MSC * m)
{
  MDR *size = mdr_Create (1, 1);
  Mdr0 (size, 0, 0) = MAX(m->nr, m->nc);
  return (size);
}

/* **************************************************************
 * Abs function...
 * ************************************************************** */

MSR *
msc_Abs (MSC * m)
{
  int i;
  MSR *new;

  new = msr_Create (m->nr, m->nc);
  msr_Setup (new, m->nnz);

  memcpy (new->ia, m->ia, (m->nr + 1) * sizeof (int));
  memcpy (new->ja, m->ja, (m->nnz) * sizeof (int));

  for (i = 0; i < m->nnz; i++)
  {
    new->d[i] = cabs(m->c[i]);
  }
  return (new);
}

/* **************************************************************
 * Detect Infs and NaNs in sparse matrices.
 * ************************************************************** */

void
msc_Detect_Inf (MSC * m)
{
  if (detect_inf_c (m->c, m->nnz))
    rerror ("matrix contains Inf value");
}

void
msc_Detect_NaN (MSC * m)
{
  if (detect_nan_c (m->c, m->nnz))
    rerror ("matrix contains NaN value");
}

/* **************************************************************
 * Test a sparse complex matrix for symmetry.
 * ************************************************************** */

int
msc_IsSymmetric (MSC * m)
{
  int i;
  Complex ctmp;
  MSC *mt;

  if (m->nr != m->nc)
    return (0);

  /* Check that M and Mt are the same. */
  mt = msc_Transpose (m);
  for (i = 0; i < m->nnz; i++)
  {
    if ((RE(m->c[i]) != RE(mt->c[i])) || (IM(m->c[i]) != IM(mt->c[i])) ||
         (m->ja[i] != mt->ja[i]))
    {
      msc_Destroy (mt);
      return (0);
    }
  }
  msc_Destroy (mt);

  /*
   * Now we must check the matrix diagonals
   * to make sure they are real.
   */

  for (i = 1; i <= m->nr; i++)
  {
    ctmp = msc_GetEl (m, i, i);
    if (IM(ctmp) != 0.0)
      return (0);
  }

  return (1);
}

/* **************************************************************
 * Find the indices of non-zero elements in a sparse matrix.
 * This is almost too easy...
 *
 * MSPARSE[i;j]
 * The index is:  (m->nr*(j-1)) + i
 * ************************************************************** */

MDR *
msc_Find_BF (MSC * m)
{
  int i, ia, j, k, n, nnz, row, col;
  MDR *find;

  nnz = m->nnz;
  find = mdr_Create (1, nnz);

  row = 1;
  k = 0;
  for (i = 0; i < m->nr; i++)
  {
    /* Loop over each row. */
    n = m->ia[i + 1] - m->ia[i];
    ia = m->ia[i];

    for (j = 0; j < n; j++)
    {
      /* Get the column number. */
      col = m->ja[ia + j - 1];
      /* Calculate the index value. */
      MdrV0 (find, k++) = (double) ((m->nr) * (col - 1) + row);
    }
    row++;
  }
  return (find);
}

/* **************************************************************
 * Return the real part of a sparse-complex matrix.
 * ************************************************************** */

MSR *
msc_Real_BF (MSC * m)
{
  int i;
  MSR *r;

  r = msr_Create (m->nr, m->nc);
  msr_Setup (r, m->nnz);

  memcpy (r->ia, m->ia, ((m->nr) + 1) * sizeof (int));
  memcpy (r->ja, m->ja, (m->nnz) * sizeof (int));

  for (i = 0; i < m->nnz; i++)
  {
    r->d[i] = RE(m->c[i]);
  }

  return (r);
}

/* **************************************************************
 * Return the imaginary part of a sparse-complex matrix.
 * ************************************************************** */

MSR *
msc_Imag_BF (MSC * m)
{
  int j;
  MSR *i;

  i = msr_Create (m->nr, m->nc);
  msr_Setup (i, m->nnz);

  memcpy (i->ia, m->ia, ((m->nr) + 1) * sizeof (int));
  memcpy (i->ja, m->ja, (m->nnz) * sizeof (int));

  for (j = 0; j < m->nnz; j++)
  {
    i->d[j] = IM(m->c[j]);
  }

  return (i);
}

MDC *
msc_Sum_BF (MSC * m, void *h)
{
  int i, j;
  MDC *rm;

  if (m->nr == 1)
  {
    /* Single row, return single double value. */
    Complex c;
    c = 0.0 + 0.0i;
    for (i = 0; i < m->nnz; i++)
    {
      c += m->c[i];
    }
    rm = mdc_CreateScalar (RE(c), IM(c));
  }
  else
  {
    /*
     * N-by-M matrix, return a vector.
     * A sum-value for each column in the matrix.
     * Transpose the matrix, in-efficient but easy.
     */

    int ia, n;
    MSC *stmp;

    stmp = msc_Transpose (m);
    rm = mdc_Create (1, m->nr);
    mdc_Zero (rm);

    /* Now, compute the sum for each row (was column). */
    for (i = 0; i < stmp->nr; i++)
    {
      n = m->ia[i + 1] - m->ia[i];
      ia = m->ia[i];

      for (j = 0; j < n; j++)
      {
        MdcV0(rm, i) = MdcV0(rm, i) + stmp->c[ia + j - 1];
      }
    }
    msc_Destroy (stmp);
  }
  return (rm);
}

double *
msc_Norm (MSC * m, char *type)
{
  double *norm = (double *) GC_MALLOC (sizeof (double));
  *norm = 0.0;

  msc_Detect_Inf (m);
  msc_Detect_NaN (m);

  /* Argument error checking */
  if (!strcmp (type, "1") || !strcmp (type, "O") || !strcmp (type, "o"))
  {
    /* largest column sum. */
    int i, size;
    MSR *mtmp = msc_Abs (m);
    MDR *msum = msr_Sum_BF (mtmp,NULL);
    msr_Destroy (mtmp);

    size = MNR (msum) * MNC (msum);
    for (i = 0; i < size; i++)
    {
      *norm = MAX(*norm, MdrV0 (msum, i));
    }
    mdr_Destroy (msum);
  }
  else if (!strcmp (type, "I") || !strcmp (type, "i"))
  {
    /*
     * largest row sum.
     * Same as 1-norm, except transpose the matrix.
     */
    int i, size;
    MSR *mtmp1 = msc_Abs (m);
    MSR *mtmp2 = msr_Transpose (mtmp1);
    MDR *msum = msr_Sum_BF (mtmp2,NULL);
    msr_Destroy (mtmp1);
    msr_Destroy (mtmp2);

    size = MNR (msum) * MNC (msum);
    for (i = 0; i < size; i++)
    {
      *norm = MAX(*norm, MdrV0 (msum, i));
    }
    mdr_Destroy (msum);
  }
  else if (!strcmp (type, "M"))
    rerror ("norm: type M not supported for sparse matrix");
  else if (!strcmp (type, "m"))
    rerror ("norm: type m not supported for sparse matrix");
  else if (!strcmp (type, "F"));
  else if (!strcmp (type, "f"));
  else if (!strcmp (type, "E"));
  else if (!strcmp (type, "e"));
  else if (!strcmp (type, "2"));
  else
    rerror ("norm: incorrect STRING specifier");

  return (norm);
}

/* **************************************************************
 * Max function: Maximum scalar value if the input is a vector.
 * Maximum row-vector output if input is a matrix.
 * ************************************************************** */

MDC *
msc_Max1 (MSC * m)
{
  int i, j;
  Complex maxc;
  MDC *mmax = 0;

  if (m->nr == 1)
  {
    mmax = mdc_CreateScalar(0.0, 0.0);
    if (m->nnz > 0)
    {
      /* Input is a row-vector, return scalar max. */
      maxc = m->c[0];
      for (i = 1; i < m->nnz; i++)
      {
        if (cabs (m->c[i]) > cabs(maxc))
          maxc = m->c[i];
      }
      mmax = mdc_Create (1, 1);
      Mdc0 (mmax, 0, 0) = maxc;
    }
  }
  else
  {
    if (m->nnz == 0)
    {
      mmax = mdc_Create (1, m->nc);
      mdc_Zero (mmax);
    }
    else
    {
     /*
      * General matrix case.
      * Go through the columns, one-by-one.
      * Transpose the matrix, to make life easier.
      */
      int ia, n;
      MSC *mtmp = msc_NcTranspose (m);
      mmax = mdc_Create (1, mtmp->nr);

      /* Row-by-row */
      for (i = 0; i < mtmp->nr; i++)
      {
        /* N(umber) in the row. */
        n = mtmp->ia[i + 1] - mtmp->ia[i];
        ia = mtmp->ia[i];
        maxc = mtmp->c[ia - 1];

        for (j = 0; j < n; j++)
        {
          if (cabs (mtmp->c[ia + j - 1]) > cabs (maxc))
            maxc = mtmp->c[ia + j - 1];
        }
        Mdc0 (mmax, 0, i) = maxc;
      }
      msc_Destroy (mtmp);
    }
  }
  return (mmax);
}

/* **************************************************************
 * Min function: Minimum scalar value if the input is a vector.
 * Minimum row-vector output if input is a matrix.
 * ************************************************************** */

MDC *
msc_Min1 (MSC * m)
{
  int i, j;
  Complex maxc;
  MDC *mmax = 0;

  if (m->nr == 1)
  {
    mmax = mdc_CreateScalar(0.0, 0.0);
    if (m->nnz > 0)
    {
      /* Input is a row-vector, return scalar max. */
      maxc = m->c[0];
      for (i = 1; i < m->nnz; i++)
      {
        if (cabs (m->c[i]) < cabs (maxc))
          maxc = m->c[i];
      }
      mmax = mdc_Create (1, 1);
      Mdc0 (mmax, 0, 0) = maxc;
    }
  }
  else
  {
    if (m->nnz == 0)
    {
      mmax = mdc_Create (1, m->nc);
      mdc_Zero (mmax);
    }
    else
    {
    /*
      * General matrix case.
      * Go through the columns, one-by-one.
      * Transpose the matrix, to make life easier.
    */
      int ia, n;
      MSC *mtmp = msc_NcTranspose (m);
      mmax = mdc_Create (1, mtmp->nr);

      /* Row-by-row */
      for (i = 0; i < mtmp->nr; i++)
      {
        /* N(umber) in the row. */
        n = mtmp->ia[i + 1] - mtmp->ia[i];
        ia = mtmp->ia[i];
        maxc = mtmp->c[ia - 1];

        for (j = 0; j < n; j++)
        {
          if (cabs (mtmp->c[ia + j - 1]) < cabs (maxc))
            maxc = mtmp->c[ia + j - 1];
        }
        Mdc0 (mmax, 0, i) = maxc;
      }
      msc_Destroy (mtmp);
    }
  }
  return (mmax);
}

/* **************************************************************
 * Int function...
 * ************************************************************** */

MSC *
msc_Int_BF (MSC * m)
{
  int i;
  MSC *im;

  im = msc_Create (m->nr, m->nc);
  msc_Setup (im, m->nnz);

  memcpy (im->ia, m->ia, (m->nr + 1) * sizeof (int));
  memcpy (im->ja, m->ja, (m->nnz) * sizeof (int));

  for (i = 0; i < im->nnz; i++)
  {
    RE(im->c[i]) = (double) ((int) RE(m->c[i]));
    IM(im->c[i]) = (double) ((int) IM(m->c[i]));
  }
  return (im);
}

/* **************************************************************
 * Ceil function...
 * ************************************************************** */

MSC *
msc_Ceil_BF (MSC * m, MDR *b, MDR *o)
{
  int i;
  MSC *cm;

  cm = msc_Create (m->nr, m->nc);
  msc_Setup (cm, m->nnz);

  memcpy (cm->ia, m->ia, (m->nr + 1) * sizeof (int));
  memcpy (cm->ja, m->ja, (m->nnz) * sizeof (int));

  for (i = 0; i < cm->nnz; i++)
  {
    RE(cm->c[i]) = errcheck (ceil (RE(m->c[i])), "ceil");
    IM(cm->c[i]) = errcheck (ceil (IM(m->c[i])), "ceil");
  }
  return (cm);
}

/* **************************************************************
 * Floor function...
 * ************************************************************** */

MSC *
msc_Floor_BF (MSC * m, MDR *b, MDR *o)
{
  int i;
  MSC *cm;

  cm = msc_Create (m->nr, m->nc);
  msc_Setup (cm, m->nnz);

  memcpy (cm->ia, m->ia, (m->nr + 1) * sizeof (int));
  memcpy (cm->ja, m->ja, (m->nnz) * sizeof (int));

  for (i = 0; i < m->nnz; i++)
  {
    RE(cm->c[i]) = errcheck (floor (RE(m->c[i])), "floor");
    IM(cm->c[i]) = errcheck (floor (IM(m->c[i])), "floor");
  }
  return (cm);
}

/* **************************************************************
 * Round (rint) function...
 * ************************************************************** */

MSC *
msc_Round_BF (MSC * m)
{
  int i;
  MSC *cm;

  cm = msc_Create (m->nr, m->nc);
  msc_Setup (cm, m->nnz);

  memcpy (cm->ia, m->ia, (m->nr + 1) * sizeof (int));
  memcpy (cm->ja, m->ja, (m->nnz) * sizeof (int));

  for (i = 0; i < cm->nnz; i++)
  {
    RE(cm->c[i]) = errcheck (rint (RE(m->c[i])), "round");
    IM(cm->c[i]) = errcheck (rint (IM(m->c[i])), "round");
  }
  return (cm);
}
