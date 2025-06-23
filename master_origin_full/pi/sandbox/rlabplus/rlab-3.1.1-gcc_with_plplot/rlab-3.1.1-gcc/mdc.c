/* mdc.c Matrix Dense Complex */

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
#include "mdr.h"
#include "mem.h"
#include "mdc.h"
#include "mdcf2.h"
#include "mdr_mdc.h"
#include "util.h"
#include "bltin1.h"
#include "complex.h"
#include "rfileio.h"
#include "mathl.h"

#include "blas.h"
#include "fi.h"

#include <stdio.h>
#include <math.h>
#include "rlab_solver_parameters_names.h"

#ifdef __riscos
#include <limits.h>
#define MAXINT INT_MAX
#endif

void mdc_Shift (MDC *m, int ud, int lr)
{
  if ((!ud) && (!lr))
    return;

  int nr = MNR(m);
  int nc = MNC(m);

  if ((nr<1)||(nc<1))
    return;

  int i, j;


  if (ud != 0)
  {
    // shift up: ud>0
    // or down: ud<0
    if ((ud >= nr)||(ud <= -nr))
    {
      mdc_Zero (m);
      return;
    }

    if (ud > 0)
    {
      for (j=0; j<nc; j++)
      {
        for (i=ud; i<nr; i++)
          Mdc0(m,i-ud,j) = Mdc0(m,i,j);
        for (i=0; i<ud; i++)
          Mdc0(m,nr-i-1,j) = 0;
      }
    }
    if (ud < 0)
    {
      for (j=0; j<nc; j++)
      {
        for (i=nr+ud-1; i>=0; i--)
          Mdc0(m,i-ud,j) = Mdc0(m,i,j);
        for (i=0; i<-ud; i++)
          Mdc0(m,i,j) = 0.0;
      }
    }
  }

  return;
}









/* **************************************************************
 * Create a complex matrix
 * ************************************************************** */

MDC *
mdc_Create (int nrow, int ncol)
{
  MDC *new;

  if (nrow < 0 || ncol < 0)
    rerror ("cannot specify a negative matrix dimension");

  new = (MDC *) GC_MALLOC (sizeof (MDC));
  if (new == 0)
    rerror ("out of memory");

  new->nrow = nrow;
  new->ncol = ncol;
  new->type = RLAB_TYPE_COMPLEX;

  if (nrow * ncol != 0)
  {
    size_t size = nrow * ncol;
    if (size >= MAXINT)
    {
      fprintf (stderr,
               "ERROR: requested matrix size caused integer overflow\n");
      rerror ("out of memory");
    }
    MDCPTR(new) = (void *) GC_malloc_atomic_ignore_off_page(size * sizeof (Complex));
    if (!MDCPTR(new))
      rerror ("out of memory");
  }
  else
  {
    new->nrow = 0;
    new->ncol = 0;
    MDCPTR(new) = 0;
  }
  new->list = 0;

  return (new);
}

MDC * mdc_CreateEmpty (int nrow, int ncol)
{
  MDC *new = (MDC *) GC_MALLOC (sizeof (MDC));
  if (new == 0)
    rerror ("out of memory");

  if (nrow * ncol)
  {
    MNR(new) = nrow;
    MNC(new) = ncol;
  }
  new->type = RLAB_TYPE_COMPLEX;
  new->list = 0;
  return (new);
}

/* **************************************************************
 * Free a matrix, and wipe out the structure members.
 * ************************************************************** */

void * mdc_Destroy (MDC * m)
{
  return md_Destroy(m);
}

void * mdc_DestroyEmpty (MDC * m)
{
  MDCPTR(m) = 0;
  if (m->list)
  {
    btree_Destroy (m->list);
    m->list = 0;
  }
  GC_FREE (m);
  return (0);
}

/* **************************************************************
 * Copy a matrix. Create the new matrix, and return the new,
 * copied matrix as the result.
 * ************************************************************** */

MDC * mdc_Copy (MDC * m)
{
  int size;
  MDC *new = 0;

  new = mdc_Create (MNR (m), MNC (m));
  size = MNR (m) * MNC (m);

  /* Make sure the matrix is not empty. */
  if (size)
    memcpy (MDCPTR(new), MDCPTR(m), size * sizeof (Complex));

  /* Copy the list if there is one. */
  if (m->list)
  {
    new->list = btree_Copy (m->list);
  }

  return (new);
}

/* **************************************************************
 * Reshape a matrix (change nrow, ncol).
 * ************************************************************** */

MDC * mdc_Reshape (MDC * m, int nrow, int ncol)
{
  MDC *new;
  if (nrow * ncol != MNR (m) * MNC (m))
  {
    fprintf (stderr, "incompatible dimensions for complex-matrix reshape\n");
    fprintf (stderr, "nrow*ncol must equal: %i\n", MNR (m) * MNC (m));
    rerror ("error");
  }

  new = mdc_Copy (m);

  /* Now, reset the row and column sizes. */
  if ((MNR (new) == 0) && (MNC (new) == 0))
  {
    /* Leave it alone. */
    ;
  }
  else
  {
    /* Fix the row and column sizes. */
    new->nrow = nrow;
    new->ncol = ncol;
  }

  return (new);
}

/* **************************************************************
 * Extend a matrix. Realloc the data, extending the data size.
 * Zero out the newly added space.
 * ************************************************************** */

MDC *
mdc_Extend (MDC * m, int nrow, int ncol)
{
  unsigned char *new;
  int nr_old = m->nrow;
  int nc_old = m->ncol;

  new = md_extend ((unsigned char *) MDCPTR(m), nr_old, nc_old, nrow, ncol, sizeof(Complex));

  m->nrow = nrow;
  m->ncol = ncol;

  GC_FREE (MDCPTR(m));
  MDCPTR(m) = (Complex *) new;

  return (m);
}

/* **************************************************************
 * Truncate a matrix.
 * ************************************************************** */

void
    mdc_Truncate (MDC * m, int nr, int nc)
{
  unsigned char *data;
  int size_a;

  // Check first
  if (m->nrow > nr || m->ncol <= nc)
  {
    size_a = sizeof(Complex);
    data = (unsigned char *) md_extend ((unsigned char *)MDCPTR(m), m->nrow, m->ncol, nr, nc, size_a);

    GC_FREE(MDCPTR(m));

    MDCPTR(m) = (Complex *) data;
    m->nrow = nr;
    m->ncol = nc;
  }
}

/* **************************************************************
 * Zero the elements of an existing matrix.
 * ************************************************************** */

void
mdc_Zero (MDC * m)
{
  int i, size;
  size = SIZE(m);

  for (i = 0; i < size; i++)
  {
    MdcV0r(m,i) = 0.0;
    MdcV0i(m,i) = 0.0;
  }
}

/* **************************************************************
 * Return a double (scalar) value from a matrix.
 * Only works if the matrix is 1x1 in size.
 * For complex matrix, return the real part.
 * ************************************************************** */

double
mdc_scalar (MDC * m)
{
  if (MNR (m) != 1 || MNC (m) != 1)
  {
    fprintf (stderr, "complex matrix size %i by %i\n", MNR (m), MNC (m));
    rerror ("cannot coerce matrix to scalar");
  }
  return (Mdc0r (m, 0, 0));
}

/* **************************************************************
 * Return a double value, real part of 1st element.
 * ************************************************************** */

double * mdc_Double (MDC * m)
{
  double *dval;

  if (MNR (m) == 0 && MNC (m) == 0)
    rerror ("empty-matrix");

  dval = (double *) GC_MALLOC (sizeof (double));
  *dval = (double) MdcV0r(m,0);
  return (dval);
}


MDR *
mdc_MatrixReal (MDC * m)
{
  MDR *mr;
  mr = mdc_coerce_mdr (m);
  return (mr);
}

/* **************************************************************
 * Print out a matrix.
 * ************************************************************** */

#include <sys/ioctl.h>

void
mdc_Print (MDC * matrix, FILE * stream)
{
  char tmp[200];
  int i, j, k, nrow, ncol, npri, rem, start, width;
  int n_print, fwidth, fprec, swidth;

  struct winsize sz;
  ioctl(0, TIOCGWINSZ, &sz);
  swidth = sz.ws_col;
  if (swidth <= 0)
    swidth = TERM_WIDTH;

  fwidth = get_fwidth ();
  fprec = get_fprec ();

  /* Check format */
  width = 2 * (MAX(fwidth, fprec + 3)) + 6;
  if (width > swidth)
  {
    rerror ("format too large for complex-matrix print");
  }

  n_print = swidth / (width + 1);
  nrow = MNR (matrix);
  ncol = MNC (matrix);
  npri = MNC (matrix) / n_print;
  rem = MNC (matrix) % n_print;

  /* Special case, empty matrix */
  if (nrow == 0 && ncol == 0)
  {
    fprintf (stream, "\t[]\n");
    fflush (stream);
    return;
  }

  start = 1;
  for (i = 0; i < npri; i++)
  {
    if (npri >= 1)
      fprintf (stream, " matrix columns %d thru %d\n",
	       n_print * i + 1, (n_print) + (n_print * i));

    for (k = 1; k <= nrow; k++)	/* print all rows */
    {
      for (j = start; j <= n_print + start - 1; j++)
      {
	if (Mdc1i (matrix, k, j) >= 0.0)
	{
	  sprintf (tmp, "%*.*g + %.*gi", fwidth, fprec, Mdc1r (matrix, k, j),
		   fprec, Mdc1i (matrix, k, j));
	}
	else
	{
	  sprintf (tmp, "%*.*g - %.*gi", fwidth, fprec, Mdc1r (matrix, k, j),
		   fprec, -Mdc1i (matrix, k, j));
	}
	fprintf (stream, "%*s", width, tmp);
      }
      fprintf (stream, "\n");
    }
    start += n_print;		/* inc our col position */
    fprintf (stream, "\n");
    fflush (stream);
  }

  /* Now come back and write out the last few colums */
  if (!rem)
    return;
  if (npri >= 1)
    fprintf (stream, " matrix columns %d thru %d\n",
	     MNC (matrix) - rem + 1, MNC (matrix));

  for (k = 1; k <= nrow; k++)
  {
    for (i = ncol - rem + 1; i <= ncol; i++)
    {
      if (Mdc1i (matrix, k, i) >= 0.0)
      {
	sprintf (tmp, "%*.*g + %.*gi", fwidth, fprec, Mdc1r (matrix, k, i),
		 fprec, Mdc1i (matrix, k, i));
      }
      else
      {
	sprintf (tmp, "%*.*g - %.*gi", fwidth, fprec, Mdc1r (matrix, k, i),
		 fprec, -Mdc1i (matrix, k, i));
      }
      fprintf (stream, "%*s", width, tmp);
    }
    fprintf (stream, "\n");
    fflush (stream);
  }
}

/*
 * Add two matrices together.
 */

MDC *
mdc_Add (MDC * m1, MDC * m2)
{
  MDC *new, *m;
  Complex c;
  int i, nr, nc, size;

  /*
   * Check sizes.
   * Check for special case (scalar + matrix) first.
   */

  if (((MNR (m1) == 1) && (MNC (m1) == 1)) ||
      ((MNR (m2) == 1) && (MNC (m2) == 1)))
  {
    if ((MNR (m1) == 1) && (MNC (m1) == 1))
    {
      c = Mdc0 (m1, 0, 0);
      nr = MNR (m2);
      nc = MNC (m2);
      m = m2;
    }
    else
    {
      c = Mdc0 (m2, 0, 0);
      nr = MNR (m1);
      nc = MNC (m1);
      m = m1;
    }

    new = mdc_Create (nr, nc);
    size = nr * nc;

    for (i = 0; i < size; i++)
      MdcV0(new,i) = MdcV0(m,i) + c;

    return (new);
  }

  if (MNR (m1) != MNR (m2))
  {
    fprintf (stderr, "matrix row sizes must be equal\n");
    fprintf (stderr, "matrix row sizes: %i and %i\n", MNR (m1), MNR (m2));
    rerror ("matrix-addition: row size mis-match");
  }

  if (MNC (m1) != MNC (m2))
  {
    fprintf (stderr, "matrix column sizes must be equal\n");
    fprintf (stderr, "matrix column sizes: %i and %i\n", MNC (m1), MNC (m2));
    rerror ("matrix-addition: column size mis-match");
  }

  new = mdc_Create (MNR (m1), MNC (m1));
  size = MNR (m1) * MNC (m1);

  for (i = 0; i < size; i++)
  {
    MdcV0r(new,i) = MdcV0r(m1,i) + MdcV0r(m2,i);
    MdcV0i(new,i) = MdcV0i(m1,i) + MdcV0i(m2,i);
  }

  return (new);
}

/*
 * Subtract two matrices.
 */

MDC *
mdc_Subtract (MDC * m1, MDC * m2)
{
  MDC *new;
  int i, size;

  /*
   * Check sizes.
   * Check for special case (scalar + matrix) first.
   */

  if (((MNR (m1) == 1) && (MNC (m1) == 1)) ||
      ((MNR (m2) == 1) && (MNC (m2) == 1)))
  {
    if ((MNR (m1) == 1) && (MNC (m1) == 1))
    {
      new = (MDC *) mdc_Create (MNR (m2), MNC (m2));
      size = MNR (m2) * MNC (m2);

      for (i = 0; i < size; i++)
      {
        MdcV0r(new,i) = MdcV0r(m1,0) - MdcV0r(m2,i);
        MdcV0i(new,i) = MdcV0i(m1,0) - MdcV0i(m2,i);
      }
      return (new);
    }
    else
    {
      new = mdc_Create (MNR (m1), MNC (m1));
      size = MNR (m1) * MNC (m1);

      for (i=0; i<size; i++)
      {
        MdcV0r(new,i) = MdcV0r(m1,i) - MdcV0r(m2,0);
        MdcV0i(new,i) = MdcV0i(m1,i) - MdcV0i(m2,0);
      }
      return (new);
    }
  }

  if (MNR (m1) != MNR (m2))
  {
    fprintf (stderr, "\tmatrix row sizes must be equal\n");
    fprintf (stderr, "\tmatrix row sizes: %i and %i\n", MNR (m1), MNR (m2));
    rerror ("matrix-subtraction: row size mis-match");
  }

  if (MNC (m1) != MNC (m2))
  {
    fprintf (stderr, "\tmatrix column sizes must be equal\n");
    fprintf (stderr, "\tmatrix column sizes: %i and %i\n", MNC (m1), MNC (m2));
    rerror ("matrix-subtraction: column size mis-match");
  }

  new = mdc_Create (MNR (m1), MNC (m1));
  size = MNR (m1) * MNC (m1);

  for (i = 0; i < size; i++)
  {
    MdcV0r(new,i) = MdcV0r(m1,i) - MdcV0r(m2,i);
    MdcV0i(new,i) = MdcV0i(m1,i) - MdcV0i(m2,i);
  }

  return (new);
}

static Complex zalpha = 1.0;
static Complex zbeta  = 0.0 + 0.0i;

MDC *
mdc_Multiply (MDC * ma, MDC * mb, int *rtype)
{
  int i, size;
  MDC *mc = 0;
  char tra, trb;

  *rtype = MATRIX_DENSE_COMPLEX;

  /* Check [ma], [mb] dimensions */
  if (MNC (ma) != MNR (mb))
  {
    /*
     * Handle condition where one of the operands is a
     * scalar, or an empty matrix.
     */

    if (MNR (ma) == 1 && MNC (ma) == 1)
    {
      size = MNR (mb) * MNC (mb);
      mc = mdc_Create (MNR (mb), MNC (mb));
      for (i = 0; i < size; i++)
      {
        MdcV0 (mc, i) = MdcV0 (ma, 0) * MdcV0 (mb, i);
      }
    }
    else if (MNR (mb) == 1 && MNC (mb) == 1)
    {
      size = MNR (ma) * MNC (ma);
      mc = mdc_Create (MNR (ma), MNC (ma));
      for (i = 0; i < size; i++)
      {
        MdcV0 (mc, i) = MdcV0 (ma, i) * MdcV0 (mb, 0);
      }
    }
    else if (MNR (ma) == 0 || MNR (mb) == 0)
    {
      mc = mdc_Create (0, 0);
    }
    else
    {
      fprintf (stderr, "\tmatrix dimensions must be consistent\n");
      fprintf (stderr, "\tmatrix A: row = %i, column = %i\n", MNR (ma),
               MNC (ma));
      fprintf (stderr, "\tmatrix B: row = %i, column = %i\n", MNR (mb),
               MNC (mb));
      rerror ("matrix-multiplication dimensions mis-match");
    }
  }
  else
  {
    int m, n, k;
    mc = mdc_Create (MNR (ma), MNC (mb));
    tra = 'N';
    trb = 'N';
    m = (int) MNR (ma);
    n = (int) MNC (mb);
    k = (int) MNC (ma);

    /* BLAS option: Complex alpha = (1.0,0.0), beta = (0.0,0.0) */
    ZGEMM (&tra, &trb, &m, &n, &k, &zalpha, MDCPTR (ma),
           &m, MDCPTR (mb), &k, &zbeta, MDCPTR (mc), &m);

    /*
     * cmmpy (MNR (ma), MNC (ma), MNC (mb), MDPTRc (ma),
     *        MDPTRc (mb), MDPTRc (mc));
     */
  }
  return (mc);
}

/*
 * Element-by-Element Multiply two matrices: mc = ma .* mb
 */

MDC *
mdc_ElMultiply (MDC * ma, MDC * mb)
{
  int i, j, size;
  MDC *mc = 0;

  /* Check [ma], [mb] dimensions */
  if ((MNC (ma) == MNC (mb)) && (MNR (ma) == MNR (mb)))
  {
    mc = mdc_Create (MNR (ma), MNC (ma));
    size = MNR (mc) * MNC (mc);
    for (i = 0; i < size; i++)
      MdcV0 (mc, i) = MdcV0 (ma, i) *	MdcV0 (mb, i);
  }
  else if ((MNR (ma) * MNC (ma)) == 1)
  {
    mc = mdc_Create (MNR (mb), MNC (mb));
    size = MNR (mc) * MNC (mc);
    for (i = 0; i < size; i++)
    {
      MdcV0 (mc, i) = MdcV0 (ma, 0) * MdcV0 (mb, i);
    }
  }
  else if ((MNR (mb) * MNC (mb)) == 1)
  {
    mc = mdc_Create (MNR (ma), MNC (ma));
    size = MNR (mc) * MNC (mc);
    for (i = 0; i < size; i++)
    {
      MdcV0 (mc, i) = MdcV0 (ma, i) * MdcV0 (mb, 0);
    }
  }
  else if (MNR (ma) == 0 || MNR (mb) == 0)
  {
    mc = mdc_Create (0, 0);
  }

  /*
   * Handle special row/column conditions...
   */

  else if (MNR (ma) == MNR (mb) && MNC (ma) == 1)
  {
    mc = mdc_Create (MNR (mb), MNC (mb));
    for (i = 0; i < MNC (mb); i++)
    {
      for (j = 0; j < MNR (mb); j++)
      {
        Mdc0 (mc, j, i) = Mdc0 (ma, j, 0) * Mdc0 (mb, j, i);
      }
    }
  }
  else if (MNR (mb) == MNR (ma) && MNC (mb) == 1)
  {
    mc = mdc_Create (MNR (ma), MNC (ma));
    for (i = 0; i < MNC (ma); i++)
    {
      for (j = 0; j < MNR (ma); j++)
      {
        Mdc0 (mc, j, i) = Mdc0 (ma, j, i) * Mdc0 (mb, j, 0);
      }
    }
  }
  else if (MNC (ma) == MNC (mb) && MNR (ma) == 1)
  {
    mc = mdc_Create (MNR (mb), MNC (mb));
    for (i = 0; i < MNC (mb); i++)
    {
      for (j = 0; j < MNR (mb); j++)
      {
        Mdc0 (mc, j, i) = Mdc0 (ma, 0, i) * Mdc0 (mb, j, i);
      }
    }
  }
  else if (MNC (mb) == MNC (ma) && MNR (mb) == 1)
  {
    mc = mdc_Create (MNR (ma), MNC (ma));
    for (i = 0; i < MNC (ma); i++)
    {
      for (j = 0; j < MNR (ma); j++)
      {
        Mdc0 (mc, j, i) = Mdc0 (ma, j, i) * Mdc0 (mb, 0, i);
      }
    }
  }
  else
  {
    fprintf (stderr, "\tmatrix dimensions must be consistent");
    fprintf (stderr, "\tmatrix A: row = %i, column = %i\n", MNR (ma), MNC (ma));
    fprintf (stderr, "\tmatrix B: row = %i, column = %i\n", MNR (mb), MNC (mb));
    rerror ("matrix el-multiply, dimension mis-match");
  }
  return (mc);
}

/*
 * m1 / m2 or B / A
 * same as B*inv(A)
 */

MDC *
mdc_Rdivide (MDC * m1, MDC * m2)
{
  int i, size;
  MDC *new=0, *tnew=0, *t1=0, *t2=0;

  /*
   * Check for special case where denominator is a scalar.
   * The only thing we gain if m2 is a scalar, is speed.
   */

  if (MNR (m2) == 1 && MNC (m2) == 1)
  {
    new = mdc_Create (MNR (m1), MNC (m1));
    size = MNR (m1) * MNC (m1);
    for (i = 0; i < size; i++)
    {
      MdcV0 (new, i) = MdcV0 (m1, i) / MdcV0 (m2, 0);
    }
    return (new);
  }
  else if (MNR (m1) == 0 || MNR (m2) == 0)
  {
    new = mdc_Create (0, 0);
    return (new);
  }

  t1 = mdc_Transpose (m1);
  t2 = mdc_Transpose (m2);

  /* Check row dimensions */
  if (MNR (t1) != MNR (t2))
  {
    mdc_Destroy (t1);
    mdc_Destroy (t2);
    rerror ("column dimensions must agree for right-divide");
  }

  if (MNR (t2) == MNC (t2))
    tnew = mdc_SolveEq (t2, t1);
  else
    tnew = mdc_LS (t2, t1);

  new = mdc_Transpose (tnew);
  mdc_Destroy (t1);
  mdc_Destroy (t2);
  mdc_Destroy (tnew);

  return (new);
}

/*
 * Element Right divide.
 * ma / mb
 */

MDC *
mdc_ElRdivide (MDC * ma, MDC * mb, int *rtype)
{
  int i, j, size;
  MDC *mc = 0;

  *rtype = MATRIX_DENSE_COMPLEX;
  /* Check [ma], [mb] dimensions */
  if ((MNR (ma) == MNR (mb)) && (MNC (ma) == MNC (mb)))
  {
    size = MNR (ma) * MNC (ma);
    mc = mdc_Create (MNR (ma), MNC (ma));
    for (i = 0; i < size; i++)
    {
      MdcV0 (mc, i) = MdcV0 (ma, i) / MdcV0 (mb, i);
    }
  }

  /*
   * Handle condition where one of the operands is a
   * scalar, or an empty matrix.
   */

  else if (MNR (ma) == 1 && MNC (ma) == 1)
  {
    size = MNR (mb) * MNC (mb);
    mc = mdc_Create (MNR (mb), MNC (mb));
    for (i = 0; i < size; i++)
    {
      MdcV0 (mc, i) = MdcV0 (ma, 0) / MdcV0 (mb, i);
#if 0
      MdcV0 (mc, i) = complex_div (MdcV0r (ma, 0), MdcV0i (ma, 0),
				   MdcV0r (mb, i), MdcV0i (mb, i));
#endif
    }
  }
  else if (MNR (mb) == 1 && MNC (mb) == 1)
  {
    size = MNR (ma) * MNC (ma);
    mc = mdc_Create (MNR (ma), MNC (ma));
    for (i = 0; i < size; i++)
    {
      MdcV0 (mc, i) = MdcV0 (ma, i) / MdcV0 (mb, 0);
#if 0
      MdcV0 (mc, i) = complex_div (MdcV0r (ma, i), MdcV0i (ma, i),
				   MdcV0r (mb, 0), MdcV0i (mb, 0));
#endif
    }
  }
  else if (MNR (ma) == 0 || MNR (mb) == 0)
  {
    mc = mdc_Create (0, 0);
  }

  /*
   * Handle special row/column conditions...
   */

  else if (MNR (ma) == MNR (mb) && MNC (ma) == 1)
  {
    mc = mdc_Create (MNR (mb), MNC (mb));
    for (i = 0; i < MNR (mb); i++)
    {
      for (j = 0; j < MNC (mb); j++)
      {
	       Mdc0 (mc, i, j) = Mdc0 (ma, i, 0) / Mdc0 (mb, i, j);
      }
    }
  }
  else if (MNR (mb) == MNR (ma) && MNC (mb) == 1)
  {
    mc = mdc_Create (MNR (ma), MNC (ma));
    for (i = 0; i < MNR (ma); i++)
    {
      for (j = 0; j < MNC (ma); j++)
      {
        Mdc0 (mc, i, j) = Mdc0 (ma, i, j) / Mdc0 (mb, i, 0);
      }
    }
  }
  else if (MNC (ma) == MNC (mb) && MNR (ma) == 1)
  {
    mc = mdc_Create (MNR (mb), MNC (mb));
    for (i = 0; i < MNR (mb); i++)
    {
      for (j = 0; j < MNC (mb); j++)
      {
	       Mdc0 (mc, i, j) = Mdc0 (ma, 0, j) / Mdc0 (mb, i, j);
      }
    }
  }
  else if (MNC (mb) == MNC (ma) && MNR (mb) == 1)
  {
    mc = mdc_Create (MNR (ma), MNC (ma));
    for (i = 0; i < MNR (ma); i++)
    {
      for (j = 0; j < MNC (ma); j++)
      {
	       Mdc0 (mc, i, j) = Mdc0 (ma, i, j) / Mdc0 (mb, 0, j);
      }
    }
  }
  else
  {
    fprintf (stderr, "\tmatrix dimensions must be consistent\n");
    fprintf (stderr, "\tmatrix A: row = %i, column = %i\n", MNR (ma), MNC (ma));
    fprintf (stderr, "\tmatrix B: row = %i, column = %i\n", MNR (mb), MNC (mb));
    rerror ("matrix-division dimensions mis-match");
  }
  return (mc);
}

/* **************************************************************
 * Matrix left-divide
 *
 * m1 \ m2 or A \ B
 * Ax = B
 *
 * ************************************************************** */

MDC *
mdc_Ldivide (MDC * m1, MDC * m2)
{
  int i, size;
  MDC *new;

  /*
   * Check for special case where denominator (m1) is a scalar.
   * The only thing we gain if m2 is a scalar, is speed.
   */

  if (MNR (m1) == 1 && MNC (m1) == 1)
  {
    new = mdc_Create (MNR (m2), MNC (m2));
    size = MNR (m2) * MNC (m2);
    for (i = 0; i < size; i++)
    {
      MdcV0 (new, i) = MdcV0 (m2, i) / MdcV0 (m1, 0);
    }
    return (new);
  }
  else if (MNR (m1) == 0 || MNR (m2) == 0)
  {
    new = mdc_Create (0, 0);
    return (new);
  }

  /* Check dimensions */
  if (MNR (m1) != MNR (m2))
    rerror ("RHS row dim. must match LHS row dim.");

  if (MNR (m1) == MNC (m1))
    new = mdc_SolveEq_GE (m1, m2);
  else
    new = mdc_LS (m1, m2);

  return (new);
}

/* **************************************************************
 * Matrix Element Left Divide
 * ************************************************************** */

MDC *
mdc_ElLdivide (MDC * m1, MDC * m2)
{
  MDC *new;
  int i, size;

  if (MNR (m1) != MNR (m2))
  {
    fprintf (stderr, "matrix row sizes must be equal\n");
    fprintf (stderr, "matrix row sizes: %i and %i\n", MNR (m1), MNR (m2));
    rerror ("matrix-element-left-divide: row size mis-match");
  }

  if (MNC (m1) != MNC (m2))
  {
    fprintf (stderr, "matrix column sizes must be equal\n");
    fprintf (stderr, "matrix column sizes: %i and %i\n", MNC (m1), MNC (m2));
    rerror ("matrix-element-left-divide: column size mis-match");
  }

  new = mdc_Create (MNR (m1), MNC (m1));
  size = MNR (m1) * MNC (m1);

  for (i = 0; i < size; i++)
  {
    MdcV0 (new, i) = MdcV0(m2,i) / MdcV0(m1,i);
  }

  return (new);
}

/* **************************************************************
 * Matrix Power operations...
 * ************************************************************** */

/*
 * Scalar ^ Matrix
 * Do the scalar ^ [1-by-1] case only.
 */

void *mdc_elpower1 (MDC * m1, MDC * m2, int *type);

void *
mdc_power1 (MDC * m1, MDC * m2, int *type)
{
  void *m = 0;

  /* Special case: scalar ^ scalar */
  if ((MNR (m2) == 1) && (MNC (m2) == 1))
  {
    m = mdc_elpower1 (m1, m2, type);
    return (m);
  }

  rerror ("scalar ^ matrix not supported yet");
  return (m);
}

/*
 * Matrix ^ Scalar
 *
 * If scalar is an integer, use matrix multiplies.
 * If scalar is not integer use eigenvalues/vectors.
 */

void *
mdc_power2 (MDC * m1, MDC * m2, int *type)
{
  void *m = 0;

  /* Special case: scalar ^ scalar */
  if ((MNR (m1) == 1) && (MNC (m1) == 1))
  {
    m = mdc_elpower1 (m1, m2, type);
    return (m);
  }

  rerror ("matrix ^ scalar not supported yet");
  return (m);

#if 0
  int i;
  MDR *mtmp1, *mtmp2;
  MDC *eval, *evec, *mtmp3, *dum;

  if (MNR (m) != MNC (m))
    rerror ("matrix ^ scalar: operand must be square matrix");

  /* Special case: scalar ^ scalar */
  if ((MNR (m1) == 1) && (MNC (m1) == 1))
  {
    mtmp1 = mdr_elpower1 (m1, m2, type);
    return (mtmp1);
  }

  if (floor (MdrV0 (m2, 0)) == MdrV0 (m2, 0))
  {
    /*
     * Use matrix Multiply method.
     * Note that we could use a more efficient
     * algorithm that only involves about half the
     * multiplies.
     */

    if (MdrV0 (m2, 0) == 0.0)
    {
      /* Return I */
      mtmp1 = mdr_Create (MNR (m1), MNC (m1));
      mdr_Zero (mtmp1);
      for (i = 0; i < MNR (m); i++)
      {
	Mdr0 (mtmp1, i, i) = 1.0;
      }
    }
    else if (MdrV0 (s, 0) < 0.0)
    {
      if (MdrV0 (m2, 0) == -1.0)
      {
	mtmp1 = mdr_Copy (m);
      }
      else
      {
	mtmp1 = mdr_Multiply (m1, m1);
      }

      for (i = 1; i < -MdrV0 (m2, 0) - 1; i++)
      {
	mtmp2 = mdr_Multiply (mtmp1, m);
	mdr_Destroy (mtmp1);
	mtmp1 = mtmp2;
      }

      /* Invert result */
      mtmp2 = mdr_Inverse (mtmp1);
      mdr_Destroy (mtmp1);
      mtmp1 = mtmp2;
    }
    else
    {
      /* Positive integer exponent */
      if (MdrV0 (s, 0) == 1.0)
      {
	mtmp1 = mdr_Copy (m);
      }
      else
      {
	mtmp1 = mdr_Multiply (m, m);
      }
      for (i = 1; i < MdrV0 (s, 0) - 1; i++)
      {
	mtmp2 = mdr_Multiply (mtmp1, m);
	mdr_Destroy (mtmp1);
	mtmp1 = mtmp2;
      }
    }
  }
  else
  {
    /*
     * Use matrix eigenvalues method.
     */

    signal (SIGINT, intcatch);
    mdr_EigS_GE (m, &eval, &evec, &dum, 0);
    signal (SIGINT, intcatch_wait);

    mtmp1 = mdc_Pow1 (eval, m2);

    /* Now form diagonal matrix from mtmp */
    mtmp2 = matrix_CreateC (MNC (mtmp1), MNC (mtmp1));
    matrix_Zero (mtmp2);
    for (i = 1; i <= MNC (mtmp2); i++)
    {
      MATr (mtmp2, i, i) = MATcvr1 (mtmp1, i);
      MATi (mtmp2, i, i) = MATcvi1 (mtmp1, i);
    }
    matrix_Destroy (mtmp1);

    /* Now do matrix multiplies to get results */
    mtmp3 = matrix_MultiplyCC (evec, mtmp2);
    mtmp1 = matrix_Rdivide (mtmp3, evec);

    /* Clean-Up */
    matrix_Destroy (mtmp2);
    matrix_Destroy (mtmp3);
    matrix_Destroy (evec);
    matrix_Destroy (eval);
  }
  return (mtmp1);
#endif
}

/*
 * Matrix ^ Matrix
 */

void *
mdc_power3 (MDC * m1, MDC * m2, int *type)
{
  void *m = 0;

  /* Special case: scalar ^ scalar */
  if ((MNR (m1) == 1) && (MNC (m1) == 1) && (MNR (m2) == 1) && (MNC (m2) == 1))
  {
    m = mdc_elpower1 (m1, m2, type);
    return (m);
  }

  rerror ("matrix ^ matrix not supported yet");
  return (m);
}

void *
mdc_Power (MDC * m1, MDC * m2, int *type)
{
  void *m = 0;
  /* Check sizes 1st. */

  if (MNR (m1) == 1 && MNC (m1) == 1)
  {
    /* Special case: 1-by-1 ^ Matrix */
    m = mdc_power1 (m1, m2, type);
  }
  else if (MNR (m2) == 1 && MNC (m2) == 1)
  {
    /* Special case: Matrix ^ 1-by-1 */
    m = mdc_power2 (m1, m2, type);
  }
  else if (MNR (m1) != MNR (m2) || MNC (m1) != MNC (m2))
  {
    rerror ("row and column dimensions must match for ^ operation");
  }
  else
  {
    /* Last possibility: Matrix ^ Matrix. */
    fprintf (stderr, "\tinvalid matrix dimensions\n");
    fprintf (stderr, "\tmatrix A: row = %i, column = %i\n", MNR (m1), MNC (m1));
    fprintf (stderr, "\tmatrix B: row = %i, column = %i\n", MNR (m2), MNC (m2));
    rerror ("matrix ^ matrix, invalid operation");
  }
  return (m);
}


/* **************************************************************
 * Element-by-element power operations...
 * ************************************************************** */

/*
 * Scalar .^ Matrix
 *
 * Here we do the operation one of two ways:
 * etype 0: Positive operand (s), use a real power function.
 *          OR integer M
 * etype 1: Everything else (default).
 */

void *
mdc_elpower1 (MDC * m1, MDC * m2, int *type)
{
  void *m;
  int i, size;

  size = MNR (m2) * MNC (m2);

  /* Special case: 0i^A */
  if ((MdcV0r (m1, 0) == 0.0) && (MdcV0i (m1, 0) == 0.0))
  {
    m = mdr_Create (MNR (m2), MNC (m2));
    mdr_Zero (m);
    *type = MATRIX_DENSE_REAL;

    /* Handle special case of zero exponent. */
    for (i = 0; i < size; i++)
    {
      if (MdcV0r (m2, i) == 0.0)
        MdrV0 (m, i) = 1.0;
    }
    return (m);
  }

  m = (void *) mdc_Create (MNR (m2), MNC (m2));
  *type = MATRIX_DENSE_COMPLEX;
  for (i = 0; i < size; i++)
  {
    MdcV0 ((MDC *) m, i) = cpow (MdcV0 (m1, 0), MdcV0 (m2, i));
  }
  return (m);
}

/*
 * Matrix .^ Scalar
 * m1 -- matrix dimensions
 * m2 -- scalar dimensions.
 */

void *
mdc_elpower2 (MDC * m1, MDC * m2, int *type)
{
  int i, size;
  void *m;

  /* Check for zero exponent. */
  if ((MdcV0r (m2, 0) == 0.0) && (MdcV0i (m2, 0) == 0.0))
  {
    m = mdr_Create (MNR (m1), MNC (m1));
    *type = MATRIX_DENSE_REAL;
    size = MNR (m1) * MNC (m1);

    for (i = 0; i < size; i++)
    {
      MdrV0 (m, i) = 1.0;
    }
    return (m);
  }

  /* Regular operation. */

  m = (void *) mdc_Create (MNR (m1), MNC (m1));
  *type = MATRIX_DENSE_COMPLEX;
  size = MNR (m1) * MNC (m1);

  for (i = 0; i < size; i++)
  {
    /* Check for zero base. */
    if ((MdcV0r (m1, i) == 0.0) && (MdcV0i (m1, i) == 0.0))
    {
      MdcV0r ((MDC *) m, i) = 0.0;
      MdcV0i ((MDC *) m, i) = 0.0;
    }
    else
    {
      MdcV0 ((MDC *) m, i) = cpow (MdcV0(m1, i), MdcV0(m2, 0));
    }
  }

  return (m);
}

/*
 * Matrix .^ Matrix
 * Here we do the operation one of two ways:
 * type 0: Positive operand (M1), use a real power function.
 *         OR integer M2
 * type 1: Everything else (default).
 */

void *
mdc_elpower3 (MDC * m1, MDC * m2, int *type)
{
  int i, size;
  void *m;

  if (MNR (m1) != MNR (m2) || MNC (m1) != MNC (m2))
  {
    fprintf (stderr, "\tmatrix dimensions must be consistent\n");
    fprintf (stderr, "\tmatrix A: row = %i, column = %i\n", MNR (m1), MNC (m1));
    fprintf (stderr, "\tmatrix B: row = %i, column = %i\n", MNR (m2), MNC (m2));
    rerror ("Element-power (.^) operation dimension mis-match");
  }

  size = MNR (m1) * MNC (m1);
  m = (void *) mdc_Create (MNR (m1), MNC (m1));
  *type = MATRIX_DENSE_COMPLEX;

  for (i = 0; i < size; i++)
  {
    if ((MdcV0r (m1, i) == 0.0) && (MdcV0i (m1, i) == 0.0))
    {
      if ((MdcV0r (m2, i) == 0.0) && (MdcV0i (m2, i) == 0.0))
      {
        MdcV0r ((MDC *) m, i) = 1.0;
        MdcV0i ((MDC *) m, i) = 1.0;
      }
      else
      {
        MdcV0r ((MDC *) m, i) = 0.0;
        MdcV0i ((MDC *) m, i) = 0.0;
      }
    }
    else
    {
      MdcV0 ((MDC *) m, i) = cpow (MdcV0(m1, i), MdcV0(m2, i));
    }
  }
  return (m);
}

void *
mdc_ElPower (MDC * m1, MDC * m2, int *type)
{
  void *m = 0;

  /* Check sizes 1st. */

  /* Special case: 1-by-1 .^ Matrix */
  if (MNR (m1) == 1 && MNC (m1) == 1)
  {
    m = mdc_elpower1 (m1, m2, type);
  }

  /* Special case: Matrix .^ 1-by-1 */
  else if (MNR (m2) == 1 && MNC (m2) == 1)
  {
    m = mdc_elpower2 (m1, m2, type);
  }

  else if (MNR (m1) != MNR (m2) || MNC (m1) != MNC (m2))
  {
    rerror ("row and column dimensions must match for .^ operation");
  }
  else
  {
    /* Last possibility: Matrix .^ Matrix. */
    m = mdc_elpower3 (m1, m2, type);
  }
  return (m);
}


/* **************************************************************
 * Negate a matrix.
 * ************************************************************** */

MDC *
mdc_Negate (MDC * m)
{
  int i, size;
  MDC *mnew;

  mnew = mdc_Create (MNR (m), MNC (m));

  size = MNR (m) * MNC (m);
  for (i = 0; i < size; i++)
  {
    MdcV0r (mnew, i) = -MdcV0r (m, i);
    MdcV0i (mnew, i) = -MdcV0i (m, i);
  }
  return (mnew);
}

/* **************************************************************
 * Return a logical-scalar for the input matrix.
 * ************************************************************** */

int
mdc_LogicalScalar (MDC * m)
{
  if ((MNR (m) == 1) && (MNC (m) == 1))
  {
    if ((MdcV0r (m, 0) == 0.0))
      return (0);
    else
      return (1);
  }
  else if ((MNR (m) == 0) && (MNC (m) == 0))
  {
    return (0);
  }
  else
  {
    fprintf (stderr, "cannot compute a logical scalar value\n");
    fprintf (stderr, "from a %i by %i matrix\n", MNR (m), MNC (m));
    rerror ("cannot compute a logical scalar value");
  }
  return (0);			/* Shut up compiler. */
}

/* **************************************************************
 * Return a the size of the matrix.
 * ************************************************************** */

int
mdc_Size (MDC * m)
{
  return (MNR (m) * MNC (m));
}

/* **************************************************************
 * Append two matrices together:  [ m1 , m2 ]
 * ************************************************************** */

MDC *
mdc_Append (MDC * m1, MDC * m2)
{
  int i, j, nrow, ncol;
  MDC *new=0;

  /* Check for empty matrices... */
  if (SIZE(m1)<1 && SIZE(m2)<1)
  {
    new = mdc_Create(0,0);
    return new;
  }
  else if (SIZE(m1)<1)
  {
    new = mdc_Copy (m2);
    return (new);
  }
  else if (SIZE(m2)<1)
  {
    new = mdc_Copy (m1);
    return (new);
  }

  /* Do the append. */
  /* Create the new matrix, large enough to hold both. */
  int r1 = MNR(m1);
  int c1 = MNC(m1);
  int r2 = MNR(m2);
  int c2 = MNC(m2);
  nrow = MAX(r1,r2);
  ncol = c1 + c2;
  new = mdc_Create (nrow, ncol);

  for (i=0; i<nrow; i++)
  {
    for (j=0; j<c1; j++)
      Mdc0(new,i,j)     = Mdc0(m1,MIN(i,r1-1),j);
    for (j=0; j<c2; j++)
      Mdc0(new,i,j+c1)  = Mdc0(m2,MIN(i,r2-1),j);
  }

  return (new);
}

/* **************************************************************
 * Stack two matrices together:  [ m1 ; m2 ]
 * ************************************************************** */

MDC *
mdc_Stack (MDC * m1, MDC * m2)
{
  int i, j, nrow, ncol;
  MDC *new=0;

  /* Check for empty matrices... */
  if (SIZE(m1)<1 && SIZE(m2)<1)
  {
    new = mdc_Create(0,0);
    return new;
  }
  else if (SIZE(m1)<1)
  {
    new = mdc_Copy (m2);
    return (new);
  }
  else if (SIZE(m1)<1)
  {
    new = mdc_Copy (m1);
    return (new);
  }

  /* Do the stack. */
  /* Create the new matrix, large enough to hold both. */
  int r1 = MNR(m1);
  int c1 = MNC(m1);
  int r2 = MNR(m2);
  int c2 = MNC(m2);
  nrow = r1+r2;
  ncol = MAX(c1,c2);
  new = mdc_Create (nrow, ncol);

  for (j=0; j<ncol; j++)
  {
    for (i=0; i<r1; i++)
      Mdc0(new,i,j)     = Mdc0(m1,i,MIN(j,c1-1));
    for (i=0; i<r2; i++)
      Mdc0(new,i+r1,j)  = Mdc0(m2,i,MIN(j,c2-1));
  }

  return (new);
}

/* **************************************************************
 * Sub-Matrix Expression Evaluation.
 * This function must look at the index expressions, and decide
 * how to act. At this point, MATRIX_DENSE_REAL only understands
 * integer-like indices. Any non-numeric index values must be
 * handled by a coercion function.
 * ************************************************************** */

MDC *
mdc_MatrixSub (MDC * var, int *i, int *j, int *type)
{
  int m, n;
  MDC *new;

  *type = MATRIX_DENSE_COMPLEX;

  /* Handle empty matrix indices. */

  if ((i == 0) || (j == 0) || (var->nrow * var->ncol == 0))
  {
    new = mdc_Create (0, 0);
    return (new);
  }

  new = mdc_Create (i[0], j[0]);

  for (m = 1; m <= i[0]; m++)
  {
    if (i[m] > MNR (var))
    {
      mdc_Destroy (new);
      fprintf (stderr, "index exceeds matrix limits\n");
      fprintf (stderr, "index value: %i, matrix size: %i\n", i[m], MNR (var));
      rerror ("sub-matrix evaluation");
    }
    for (n = 1; n <= j[0]; n++)
    {
      if (j[n] > MNC (var))
      {
	mdc_Destroy (new);
	fprintf (stderr, "index exceeds matrix limits\n");
	fprintf (stderr, "index value: %i, matrix size: %i\n", j[n], MNC (var));
	rerror ("sub-matrix evaluation");
      }
      Mdc1 (new, m, n) = Mdc1 (var, i[m], j[n]);
    }
  }

  return (new);
}

MDC *
mdc_MatrixSubR (MDC * var, int *i, int *type)
{
  int m, n;
  MDC *new;

  *type = MATRIX_DENSE_COMPLEX;

  /* Handle empty matrix indices. */

  if ((i == 0) || (var->nrow * var->ncol == 0))
  {
    new = mdc_Create (0, 0);
    return (new);
  }

  new = mdc_Create (i[0], MNC (var));

  for (m = 1; m <= i[0]; m++)
  {
    if (i[m] > MNR (var))
    {
      mdc_Destroy (new);
      fprintf (stderr, "index exceeds matrix limits\n");
      fprintf (stderr, "index value: %i, matrix size: %i\n", i[m], MNR (var));
      rerror ("sub-matrix evaluation");
    }
    for (n = 1; n <= MNC (var); n++)
    {
      Mdc1 (new, m, n) = Mdc1 (var, i[m], n);
    }
  }

  return (new);
}

MDC *
mdc_MatrixSubC (MDC * var, int *j, int *type)
{
  int m, n;
  MDC *new;

  *type = MATRIX_DENSE_COMPLEX;

  /* Handle empty matrix indices. */

  if ((j == 0) || (var->nrow * var->ncol == 0))
  {
    new = mdc_Create (0, 0);
    return (new);
  }

  new = mdc_Create (MNR (var), j[0]);

  for (m = 1; m <= MNR (var); m++)
  {
    for (n = 1; n <= j[0]; n++)
    {
      if (j[n] > MNC (var))
      {
	mdc_Destroy (new);
	fprintf (stderr, "index exceeds matrix limits\n");
	fprintf (stderr, "index value: %i, matrix size: %i\n", j[n], MNC (var));
	rerror ("sub-matrix evaluation");
      }
      Mdc1 (new, m, n) = Mdc1 (var, m, j[n]);
    }
  }

  return (new);
}

/* **************************************************************
 * Assign to a range of a matrix. We do not create a new matrix.
 * Automatically extend the size of a matrix if I or J indices
 * exceed current bounds.
 * ************************************************************** */

MDC *
mdc_MatrixAssign (MDC * var, int *i, int *j, MDC * rhs)
{
  int m, n;
  int dflag;
  MDC *mtmp;

  if ((i == 0) || (j == 0))
  {
    if (MdcNR (rhs) * MdcNC (rhs) != 0)
    {
      rerror ("matrix-assign: LHS and RHS must have same dimensions");
    }
    else if (var != 0)
    {
      return (var);
    }
    else
    {
      rerror
	("matrix-assign: cannot assign to an UNDEFINED object without index specification");
    }
  }

  /* Check LHS. */

  if (var == 0)
  {
    var = mdc_Create (1, 1);
    mdc_Zero (var);
  }

  /* Check RHS for empty matrix. */

  if (MNR (rhs) * MNC (rhs) == 0)
  {
    fprintf (stderr, "cannot assign empty matrix to element(s) of a matrix\n");
    rerror ("matrix-assign");
  }

  /* Check the LHS, and RHS dimensions. */

  if (i[0] != MNR (rhs))
  {
    fprintf (stderr, "LHS, and RHS row dimensions must match\n");
    fprintf (stderr, "LHS row: %i, RHS row: %i\n", i[0], MNR (rhs));
    rerror ("matrix-assign");
  }

  if (j[0] != MNC (rhs))
  {
    fprintf (stderr, "LHS, and RHS column dimensions must match\n");
    fprintf (stderr, "LHS column: %i, RHS column: %i\n", j[0], MNC (rhs));
    rerror ("matrix-assign");
  }

  /* Now do the assign. */

  if (var == rhs)
  {
    mtmp = mdc_Copy (rhs);
    dflag = 1;
  }
  else
  {
    mtmp = rhs;
    dflag = 0;
  }

  if (MNR (var) == 0)
    mdc_Extend (var, 1, 1);

  for (m = 1; m <= i[0]; m++)
  {
    if (i[m] > MNR (var))
    {
      var = mdc_Extend (var, i[m], MNC (var));
    }
    for (n = 1; n <= j[0]; n++)
    {
      if (j[n] > MNC (var))
      {
	var = mdc_Extend (var, MNR (var), j[n]);
      }
      Mdc1 (var, i[m], j[n]) = Mdc1 (mtmp, m, n);
    }
  }

  if (dflag)
    mdc_Destroy (mtmp);

  return (var);
}

MDC *
mdc_MatrixAssignR (MDC * var, int *i, MDC * rhs)
{
  int m, n;
  int dflag;
  MDC *mtmp;

  if (i == 0)
  {
    if (MdcNR (rhs) * MdcNC (rhs) != 0)
    {
      rerror ("matrix-assign: LHS and RHS must have same dimensions");
    }
    else if (var != 0)
    {
      return (var);
    }
    else
    {
      rerror
	("matrix-assign: cannot assign to an UNDEFINED object without index specification");
    }
  }

  /* Check LHS. */

  if (var == 0)
  {
    var = mdc_Create (1, 1);
    mdc_Zero (var);
  }

  /*
   * Check RHS for empty matrix.
   * If RHS is empty, then eliminate the
   * rows specified in i.
   */

  if (MNR (rhs) * MNC (rhs) == 0)
  {
    int copy, j, k, n, *rs1, rs1_size, *rs2, rs2_size;

    /* First, get rid of duplicates in the row-spec. */

    rs1 = remove_duplicates (i + 1, i[0], &rs1_size);

    /* Throw out row-spec elements greater then MNR. */

    rs2 = remove_invalid (rs1, rs1_size, MNR (var), &rs2_size);
    GC_FREE (rs1);

    /* Now, create the new matrix with the right size. */
    mtmp = mdc_Create (MNR (var) - rs2_size, MNC (var));

    /* We might as well exit now, if mtmp = [] */
    if (MNR (mtmp) == 0)
    {
      GC_FREE (rs2);
      mdc_Destroy (var);
      var = mtmp;
      return (var);
    }

    /*
     * Now, copy the correct elements into the smaller matrix.
     * This must be done element at a time, cause the matrix
     * is stored columnwise.
     */

    n = 0;
    for (j = 0; j < MNR (var); j++)
    {
      copy = 1;
      for (k = 0; k < rs2_size; k++)
      {
	if (rs2[k] == (j + 1))
	{
	  copy = 0;
	  break;
	}
      }

      if (copy)
      {
	for (k = 0; k < MNC (var); k++)
	  Mdc0 (mtmp, n, k) = Mdc0 (var, j, k);
	n++;
      }
    }

    /* Free the old, and replace it with the new. */
    GC_FREE (rs2);
    mdc_Destroy (var);
    var = mtmp;
    return (var);
  }

  /* Check the LHS, and RHS dimensions. */

  if (MNR (var) != 0)
  {
    if (i[0] != MNR (rhs))
    {
      fprintf (stderr, "LHS, and RHS row dimensions must match\n");
      fprintf (stderr, "LHS row: %i, RHS row: %i\n", i[0], MNR (rhs));
      rerror ("matrix-assign");
    }
  }
  else
  {
    mdc_Extend (var, 1, 1);
  }

  /* Now do the assign. */

  if (var == rhs)
  {
    mtmp = mdc_Copy (rhs);
    dflag = 1;
  }
  else
  {
    mtmp = rhs;
    dflag = 0;
  }

  for (m = 1; m <= i[0]; m++)
  {
    if (i[m] > MNR (var))
    {
      var = mdc_Extend (var, i[m], MNC (var));
    }
    for (n = 1; n <= MNC (rhs); n++)
    {
      if (n > MNC (var))
      {
	var = mdc_Extend (var, MNR (var), n);
      }
      Mdc1 (var, i[m], n) = Mdc1 (mtmp, m, n);
    }
  }

  if (dflag)
    mdc_Destroy (mtmp);

  return (var);
}

MDC *
mdc_MatrixAssignC (MDC * var, int *j, MDC * rhs)
{
  int m, n;
  int dflag;
  MDC *mtmp;

  if (j == 0)
  {
    if (MdcNR (rhs) * MdcNC (rhs) != 0)
    {
      rerror ("matrix-assign: LHS and RHS must have same dimensions");
    }
    else if (var != 0)
    {
      return (var);
    }
    else
    {
      rerror
	("matrix-assign: cannot assign to an UNDEFINED object without index specification");
    }
  }

  /* Check LHS. */

  if (var == 0)
  {
    var = mdc_Create (1, 1);
    mdc_Zero (var);
  }

  /*
   * Check RHS for empty matrix.
   * If RHS is empty, then eliminate the
   * rows specified in i.
   */

  if (MNR (rhs) * MNC (rhs) == 0)
  {
    int copy, jj, k, n, *rs1, rs1_size, *rs2, rs2_size;

    /* First, get rid of duplicates in the row-spec. */

    rs1 = remove_duplicates (j + 1, j[0], &rs1_size);

    /* Throw out row-spec elements greater then MNR. */

    rs2 = remove_invalid (rs1, rs1_size, MNR (var), &rs2_size);
    GC_FREE (rs1);

    /* Now, create the new matrix with the right size. */
    mtmp = mdc_Create (MNR (var), MNC (var) - rs2_size);

    /* We might as well exit now, if mtmp = [] */
    if (MNR (mtmp) == 0)
    {
      GC_FREE (rs2);
      mdc_Destroy (var);
      var = mtmp;
      return (var);
    }

    /*
     * Now, copy the correct elements into the smaller matrix.
     */

    n = 0;
    for (jj = 0; jj < MNC (var); jj++)
    {
      copy = 1;
      for (k = 0; k < rs2_size; k++)
      {
	if (rs2[k] == (jj + 1))
	{
	  copy = 0;
	  break;
	}
      }

      if (copy)
      {
	for (k = 0; k < MNR (var); k++)
	  Mdc0 (mtmp, k, n) = Mdc0 (var, k, jj);
	n++;
      }
    }

    /* Free the old, and replace it with the new. */
    GC_FREE (rs2);
    mdc_Destroy (var);
    var = mtmp;
    return (var);
  }

  /* Check the LHS, and RHS dimensions. */

  if (MNR (var) != 0)
  {
    if (j[0] != MNC (rhs))
    {
      fprintf (stderr, "LHS, and RHS column dimensions must match\n");
      fprintf (stderr, "LHS column: %i, RHS column: %i\n", j[0], MNC (rhs));
      rerror ("matrix-assign");
    }
  }
  else
  {
    mdc_Extend (var, 1, 1);
  }

  /* Now do the assign. */

  if (var == rhs)
  {
    mtmp = mdc_Copy (rhs);
    dflag = 1;
  }
  else
  {
    mtmp = rhs;
    dflag = 0;
  }

  for (m = 1; m <= MNR (rhs); m++)
  {
    if (m > MNR (var))
    {
      var = mdc_Extend (var, m, MNC (var));
    }
    for (n = 1; n <= j[0]; n++)
    {
      if (j[n] > MNC (var))
      {
	var = mdc_Extend (var, MNR (var), j[n]);
      }
      Mdc1 (var, m, j[n]) = Mdc1 (mtmp, m, n);
    }
  }

  if (dflag)
    mdc_Destroy (mtmp);

  return (var);
}

/* **************************************************************
 * Vector Sub-Expression
 * ************************************************************** */

MDC *
mdc_VectorSub (MDC * m, int *i, int *type)
{
  int j, size, msize;
  MDC *new;

  *type = MATRIX_DENSE_COMPLEX;

  /* Handle empty matrix indices. */

  if (i == 0)
  {
    new = mdc_Create (0, 0);
    return (new);
  }

  size = i[0];
  msize = MNR (m) * MNC (m);

  if (MNC (m) == 1)
    new = mdc_Create (size, 1);
  else
    new = mdc_Create (1, size);

  for (j = 1; j <= size; j++)
  {
    if (i[j] > msize)
    {
      mdc_Destroy (new);
      fprintf (stderr, "\tindex exceeds matrix limits\n");
      fprintf (stderr, "\tindex value: %i, matrix size: %i\n", i[j], msize);
      rerror ("sub-matrix evaluation");
    }
    MdcV1 (new, j) = MdcV1 (m, i[j]);
  }

  return (new);
}

/* **************************************************************
 * Assign into a matrix like a vector.
 * var[j] = rhs;
 * ************************************************************** */

MDC *
mdc_VectorAssign (MDC * var, int *j, MDC * rhs)
{
  int i, size;
  int nr, nc;
  int dflag;
  MDC *mtmp;

  if (j == 0)
  {
    if (MdcNR (rhs) * MdcNC (rhs) != 0)
    {
      rerror ("vector-assign: LHS and RHS must have same dimensions");
    }
    else if (var != 0)
    {
      return (var);
    }
    else
    {
      rerror
	("vector-assign: cannot assign to an UNDEFINED object without index specification");
    }
  }

  /* Check LHS for UNDEF. */
  if (var == 0)
  {
    var = mdc_Create (1, 1);
    mdc_Zero (var);
  }

  /* Check RHS for empty matrix. */

  if (MNR (rhs) * MNC (rhs) == 0)
  {
    fprintf (stderr,
	     "\tcannot assign empty matrix to element(s) of a matrix\n");
    rerror ("vector-assign");
  }

  /* Check the dimensions. */

  if (j[0] != (MNR (rhs) * MNC (rhs)))
  {
    fprintf (stderr, "\tLHS, and RHS column sizes must match\n");
    fprintf (stderr, "\tLHS: %i, RHS: %i\n", j[0], MNR (rhs) * MNC (rhs));
    rerror ("vector-assign");
  }

  /* Now do the assign. */

  if (var == rhs)
  {
    mtmp = mdc_Copy (rhs);
    dflag = 1;
  }
  else
  {
    mtmp = rhs;
    dflag = 0;
  }

  size = j[0];

  if (MNR (var) == 0)
    mdc_Extend (var, 1, 1);

  for (i = 1; i <= size; i++)
  {
    if (j[i] > MNR (var) * MNC (var))
    {
      if (MNR (var) == 1)
      {
	nr = 1;
	nc = (int) ceil (((double) j[i]) / ((double) MNR (var)));
      }
      else if (MNC (var) == 1)
      {
	nr = (int) ceil (((double) j[i]) / ((double) MNC (var)));
	nc = 1;
      }
      else
      {
	nr = MNR (var);
	nc = (int) ceil (((double) j[i]) / ((double) MNR (var)));
      }
      var = mdc_Extend (var, nr, nc);
    }
    MdcV1 (var, (j[i])) = MdcV1 (mtmp, i);
  }

  if (dflag)
    mdc_Destroy (mtmp);

  return (var);
}

/* **************************************************************
 * Return the i'th value needed during for-loop execution.
 * Do indexing from 1 for [m].
 * ************************************************************** */

MDC *
mdc_ForLoopValue (MDC * m, int i)
{
  MDC *new = mdc_Create (1, 1);
  MdcV0 (new, 0) = MdcV1 (m, i);
  return (new);
}

/*
 * Create a "scalar" (1-by-1) COMPLEX matrix.
 */

MDC *
mdc_CreateScalar (double rval, double ival)
{
  MDC *new;

  new = mdc_Create (1, 1);
  MdcV0r (new, 0) = rval;
  MdcV0i (new, 0) = ival;
  return (new);
}

/* **************************************************************
 * Return the class of a matrix-dense-complex.
 * ************************************************************** */

char *
mdc_Class (MDC * m)
{
  return (cpstr ("num"));
}

/* **************************************************************
 * Transpose a matrix.
 * ************************************************************** */

MDC *
    mdc_Transpose (MDC * m)
{
  int i, j;
  MDC *new = mdc_Create (MNC (m), MNR (m));

  for (i = 1; i <= MNR (m); i++)
  {
    for (j = 1; j <= MNC (m); j++)
    {
      Mdc1r (new, j, i) = Mdc1r (m, i, j);
      Mdc1i (new, j, i) = -Mdc1i (m, i, j);
    }
  }

  return (new);
}

MDC *
  mdc_NcTranspose (MDC * m)
{
  int i, j;
  MDC *new = mdc_Create (MNC (m), MNR (m));

  for (i = 1; i <= MNR (m); i++)
  {
    for (j = 1; j <= MNC (m); j++)
    {
      Mdc1r (new, j, i) = Mdc1r (m, i, j);
      Mdc1i (new, j, i) = Mdc1i (m, i, j);
    }
  }

  return (new);
}

void
mdc_NcTranspose_inplace (MDC * z)
{
  int nr = z->nrow;
  int nc = z->ncol;

  md_transpose_insitu((unsigned char *) MDCPTR(z), nr, nc, sizeof(Complex));

  z->nrow = nc;
  z->ncol = nr;
}


/* **************************************************************
 * Reshape a matrix into a colum vector.
 * ************************************************************** */

MDC *
mdc_ReshapeCol (MDC * m)
{
  MDC *new = mdc_Reshape (m, MNR (m) * MNC (m), 1);
  return (new);
}

/*
 * Write out a matrix in a generic sort of format
 */

void
mdc_WriteGeneric (MDC * m, FILE * fn)
{
  int i, j, fwidth, fprec, nfmt;

  if (m->nrow ==0 || m->ncol == 0)
    return;

  char *eol = get_eol_file_ds_name (fn);
  char *csp = get_csp_file_ds_name (fn);
  MDS  *fmt = get_fmt_file_ds_name (fn);

  nfmt = fmt->nrow * fmt->ncol - 1;

  fwidth = get_fwidth ();
  fprec  = get_fprec ();

  for (i = 0; i < m->nrow; i++)
  {
    // print first entry
    if (strlen(MdsV0(fmt,0))>1)
      fprintf (fn, MdsV0(fmt,0), Mdc0r (m, i, 0));
    else
      fprintf (fn, "%*.*f", fwidth, fprec, Mdc0r (m, i, 0));
    if (Mdc0i (m, i, 0) != 0)
    {
      if (Mdc0i (m, i, 0)>0)
        fprintf (fn, "+");
      if (strlen(MdsV0(fmt,0)) > 1)
        fprintf (fn, MdsV0(fmt,0), Mdc0i (m, i, 0));
      else
        fprintf (fn, "%*.*f", fwidth, fprec, Mdc0i (m, i, 0));
      fprintf (fn, "i");
    }

    for (j = 1; j < m->ncol; j++)
    {
      // print all the other entries in the line
      fprintf (fn, csp);
      if (strlen(MdsV0(fmt,MIN(j,nfmt))) > 1)
        fprintf (fn, MdsV0(fmt,MIN(j,nfmt)), Mdc0r (m, i, j));
      else
        fprintf (fn, "%*.*f", fwidth, fprec, Mdc0r (m, i, j));
      if (Mdc0i (m, i, j) != 0)
      {
        if (Mdc0i (m, i, j)>0)
          fprintf (fn, "+");
        if (strlen(MdsV0(fmt,MIN(j,nfmt))) > 1)
          fprintf (fn, MdsV0(fmt,MIN(j,nfmt)), Mdc0i (m, i, j));
        else
          fprintf (fn, "%*.*f", fwidth, fprec, Mdc0i (m, i, j));

        fprintf (fn, "i");
      }
    }
    fprintf (fn, eol);
  }

/*  for (i = 0; i < MdcNR (m); i++)
  {
    for (j = 0; j < MdcNC (m); j++)
    {
      fprintf (fn, " %*.*g", fwidth, fprec, Mdc0i (m, i, j));
    }
    fprintf (fn, "\n");
  }*/
}

int
mdc_Sublist (MD * m, char *name)
{
  if (!strcmp (name, RLAB_MEMBER_NROW))
  {
    return (1);
  }
  else if (!strcmp (name, RLAB_MEMBER_NCOL))
  {
    return (1);
  }
  else if (!strcmp (name, RLAB_MEMBER_SIZE))
  {
    return (1);
  }
  else if (!strcmp (name, RLAB_MEMBER_CLASS))
  {
    return (1);
  }
  else if (!strcmp (name, RLAB_MEMBER_TYPE))
  {
    return (1);
  }

  return (0);
}

//
//
//
Complex mdc0(MD *m, int k, int j)
{
  Complex cval;

  switch(m->type)
  {
    case RLAB_TYPE_DOUBLE:
      RE (cval) = Mdr0(m,k,j);
      IM (cval) = 0.0;
      break;

    case RLAB_TYPE_INT32:
      RE (cval) = (double) Mdi0(m, k, j);
      IM (cval) = 0.0;
      break;

    case RLAB_TYPE_COMPLEX:
      cval = Mdc0(m, k, j);
      break;

    default:
      RE (cval) = create_nan();
      IM (cval) = create_nan();
      break;
  }
  return cval;
}

Complex mdc1(MD *m, int k, int j)
{
  Complex cval;

  switch(m->type)
  {
    case RLAB_TYPE_DOUBLE:
      RE (cval) = Mdr1(m,k,j);
      IM (cval) = 0.0;
      break;

    case RLAB_TYPE_INT32:
      RE (cval) = (double) Mdi1(m, k, j);
      IM (cval) = 0.0;
      break;

    case RLAB_TYPE_COMPLEX:
      cval = Mdc1(m, k, j);
      break;

    default:
      RE (cval) = create_nan();
      IM (cval) = create_nan();
      break;
  }
  return cval;
}

Complex mdcV0(MD *m, int k)
{
  Complex cval;

  switch(m->type)
  {
    case RLAB_TYPE_DOUBLE:
      RE (cval) = MdrV0(m,k);
      IM (cval) = 0.0;
      break;

    case RLAB_TYPE_INT32:
      RE (cval) = (double) MdiV0(m, k);
      IM (cval) = 0.0;
      break;

    case RLAB_TYPE_COMPLEX:
      cval = MdcV0(m, k);
      break;

    default:
      RE (cval) = create_nan();
      IM (cval) = create_nan();
      break;
  }
  return cval;
}

Complex mdcV1(MD *m, int k)
{
  Complex cval;

  switch(m->type)
  {
    case RLAB_TYPE_DOUBLE:
      RE (cval) = MdrV1(m,k);
      IM (cval) = 0.0;
      break;

    case RLAB_TYPE_INT32:
      RE (cval) = (double) MdiV1(m, k);
      IM (cval) = 0.0;
      break;

    case RLAB_TYPE_COMPLEX:
      cval = MdcV1(m, k);
      break;

    default:
      RE (cval) = create_nan();
      IM (cval) = create_nan();
      break;
  }
  return cval;
}



