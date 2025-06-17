/* msc.c Matrix Sparse Complex */

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
#include "mdr.h"
#include "mds.h"
#include "mdc.h"
#include "msc.h"
#include "mdr_mdc.h"
#include "complex.h"
#include "util.h"
#include "bltin1.h"
#include "mathl.h"
#include "rfileio.h"
#include "rlab_solver_parameters_names.h"

#include "fi.h"
#include "sparse.h"

#include <stdio.h>
#include <math.h>

/* **************************************************************
 * Create a sparse matrix. Do not allocate any space for the data.
 * ************************************************************** */

MSC *
msc_Create (int nrow, int ncol)
{
  MSC *new;

  if (nrow < 0 || ncol < 0)
  {
    fprintf (stderr, "illegal dimensions for matrix: ");
    fprintf (stderr, "nrow=%i, ncol=%i\n", nrow, ncol);
    rerror ("cannot specify a negative matrix dimension");
  }

  new = (MSC *) GC_MALLOC (sizeof (MSC));
  if (new == 0)
    rerror ("out of memory");

  new->nr = nrow;
  new->nc = ncol;
  new->nnz = 0;
  new->ia = 0;
  new->ja = 0;
  new->order = 1;
  new->c = 0;
  new->list = 0;

  return (new);
}

/* **************************************************************
 * Setup a sparse complex matrix structure. NR, and NC are already
 * set. Now that we know the number of non-zeros, get the storage.
 * ************************************************************** */

MSC *
msc_Setup (MSC * m, int nnz)
{
  m->nnz = nnz;

  if (nnz)
  {
    if (!m->ia)
    {
      m->ia = (int *) GC_MAIOP ((m->nr + 1) * sizeof (int));
    }

    if (!m->ja)
    {
      m->ja = (int *) GC_MAIOP ((m->nnz) * sizeof (int));
    }

    if (!m->c)
    {
      m->c = (Complex *) GC_MAIOP ((m->nnz) * sizeof (Complex));
    }
  }

  return (m);
}

/* **************************************************************
 * Free a sparse matrix, and wipe out the structure members.
 * ************************************************************** */

void *
msc_Destroy (MSC * m)
{
  m->nr = -1;
  m->nc = -1;
  m->nnz = -1;
  if (m->ia)
    GC_FREE (m->ia);

  if (m->ja)
    GC_FREE (m->ja);

  m->order = -1;

  if (m->c)
    GC_FREE (m->c);

  if (m->list)
  {
    btree_Destroy (m->list);
    m->list = 0;
  }

  GC_FREE (m);
  return (0);
}

/* **************************************************************
 * Copy a sparse matrix. Create the new matrix, and return the new,
 * copied matrix as the result.
 * ************************************************************** */

MSC *
msc_Copy (MSC * m)
{
  MSC *new = 0;

  new = msc_Create (m->nr, m->nc);
  new->nnz = m->nnz;

  /* Copy the IA array. */
  if (m->ia)
  {
    new->ia = (int *) GC_MAIOP (sizeof (int) * (m->nr + 1));
    memcpy (new->ia, m->ia, (m->nr + 1) * sizeof (int));
  }

  /* Copy the JA array. */
  if (m->ja)
  {
    new->ja = (int *) GC_MAIOP (sizeof (int) * (m->nnz));
    memcpy (new->ja, m->ja, m->nnz * sizeof (int));
  }

  new->order = m->order;

  /* Copy the element-data array. */
  if (m->c)
  {
    new->c = (Complex *) GC_MAIOP (sizeof (Complex) * (m->nnz));
    memcpy (new->c, m->c, m->nnz * sizeof (Complex));
  }

  if (m->list)
  {
    new->list = btree_Copy (m->list);
  }
  return (new);
}

/* **************************************************************
 * Re-sparse a sparse matrix. Removes the zeros...
 * ************************************************************** */

MSC *
msc_ReSparse (MSC * m)
{
  int i, ia, inrow, j, *newia, n, nnz;
  MSC *s;

  /*
   * Count the non-zeros.
   * Figure out the IA structure.
   */

  nnz = 0;
  newia = (int *) GC_MAIOP ((m->nr + 1) * sizeof (int));
  newia[0] = 1;
  inrow = 1;

  for (i = 0; i < m->nr; i++)
  {
    n = m->ia[i + 1] - m->ia[i];
    ia = m->ia[i];
    inrow = 0;

    /* Go through each row... */
    for (j = 0; j < n; j++)
    {
      if ( cabs(m->c[ia + j - 1]) > 0.0)
      {
        inrow++;
        nnz++;      /* Total number of non-zeros. */
      }
    }
    newia[i + 1] = newia[i] + inrow;
  }

  /* Create the new structure. */
  s = msc_Create (m->nr, m->nc);
  s->nnz = nnz;
  s->ia = newia;
  s->ja = (int *) GC_MAIOP (nnz * sizeof (int));
  s->c = (Complex *) GC_MAIOP (nnz * sizeof (Complex));

  /*
   * Now assign into the new structure.
   * If we did things OK, then we can just blindly
   * assign.
   */

  j = 0;
  for (i = 0; i < m->nnz; i++)
  {
    if ( cabs(m->c[i]) > 0.0 )
    {
      s->c[j] = m->c[i];
      s->ja[j] = m->ja[i];
      j++;
    }
  }

  return (s);
}

/* **************************************************************
 * Print out a sparse-complex matrix.
 * ************************************************************** */

void
msc_Print (MSC * m, FILE * stream)
{
  int i, ia, j, n, row;
  int fwidth, fprec;

  fwidth = get_fwidth ();
  fprec = get_fprec ();

  if ((m->nr == 0) && (m->nc == 0))
  {
    fprintf (stream, "\t[]\n");
    fflush (stream);
    return;
  }

  if (m->nnz == 0)
  {
    fprintf (stream, "\t[all zeros]\n");
    fflush (stream);
    return;
  }

  row = 1;
  for (i = 0; i < m->nr; i++)
  {
    n = m->ia[i + 1] - m->ia[i];
    ia = m->ia[i];

    for (j = 0; j < n; j++)
    {
      fprintf (stream, " (%i, %i)    \t%*.*g \t%*.*gi\n",
	       row, m->ja[ia + j - 1],
	       fwidth, fprec, RE(m->c[ia + j - 1]),
	       fwidth, fprec, IM(m->c[ia + j - 1]) );
    }
    row++;
  }
  fflush (stream);
  return;
}

/* **************************************************************
 * Write a sparse-complex matrix in ASCII format to a file.
 * 1st column:    row IDs
 * 2nd column:    column IDs
 * 3rd column:    REAL element value
 * 4th column:    Imaginary element value
 * ************************************************************** */

void
msc_WriteGeneric (MSC * m, FILE * fn)
{
  int i, ia, j, n, row, nfmt;
  int fwidth, fprec;

  if (m->nnz == 0)
    return;

  char *eol = get_eol_file_ds_name (fn);
  char *csp = get_csp_file_ds_name (fn);
  MDS  *fmt = get_fmt_file_ds_name (fn);

  nfmt = fmt->nrow * fmt->ncol - 1;

  fwidth = get_fwidth ();
  fprec = get_fprec ();

  row = 1;
  for (i = 0; i < m->nr; i++)
  {
    n = m->ia[i + 1] - m->ia[i];
    ia = m->ia[i];

    for (j = 0; j < n; j++)
    {
      fprintf (fn, "%i", row);
      fprintf (fn, csp);
      fprintf (fn, "%i", m->ja[ia + j - 1]);
      fprintf (fn, csp);
      if (strlen(MdsV0(fmt,0))>1)
      {
        fprintf (fn, MdsV0(fmt,0), RE(m->c[ia + j - 1]));
        fprintf (fn, csp);
        fprintf (fn, MdsV0(fmt,MIN(nfmt,1)), IM(m->c[ia + j - 1]));
      }
      else
      {
        fprintf (fn, "%*.*f", fwidth, fprec, RE(m->c[ia + j - 1]));
        fprintf (fn, csp);
        fprintf (fn, "%*.*f", fwidth, fprec, IM(m->c[ia + j - 1]));
      }
      fprintf (fn, eol);
    }
    row++;
  }

  fflush (fn);
  return;
}

/* **************************************************************
 * Write out a sparse matrix in "graph" format. At least like
 * the graph format used in Metis and Chaco.
 *
 * Format:
 * < Number of vertices (rows) >  < Number of edges (nnz) >
 * < adjancey list for the ith vertex (row) >
 *
 * ************************************************************** */

void
msc_WriteGraph (MSC * m, FILE * fn)
{
  int i, ia, j, n;

  /* Write the header line. */
  fprintf (fn, "%i %i\n", m->nr, (m->nnz) / 2);

  /* Write the edges (adjancey list) for each row. */
  for (i = 0; i < m->nr; i++)
  {
    n = m->ia[i + 1] - m->ia[i];
    ia = m->ia[i];
    for (j = 0; j < n; j++)
    {
      fprintf (fn, " %i", m->ja[ia + j - 1]);
    }
    fprintf (fn, "\n");
  }
  fflush (fn);
  return;
}

/* **************************************************************
 * Write out a sparse matrix in compressed row-wise format.
 * ************************************************************** */

void
msc_WriteSparse (MSC * m, FILE * fn)
{
  int i;

  /* Write the header line. */
  fprintf (fn, "%i %i %i\n", m->nr, m->nc, m->nnz);

  /* Write out IA. */
  for (i = 0; i < m->nr + 1; i++)
  {
    fprintf (fn, "%i\n", m->ia[i]);
  }

  /* Write out JA, and D. */
  for (i = 0; i < m->nnz; i++)
  {
    fprintf (fn, "%i %g %g\n", m->ja[i], RE(m->c[i]), IM(m->c[i]));
  }

  fflush (fn);
  return;
}

/* **************************************************************
 * Return the class of a matrix-sparse-real.
 * ************************************************************** */

char *
msc_Class (MSC * m)
{
  return (cpstr ("num"));
}

/* **************************************************************
 * Return matrix-sparse-real member references.
 * ************************************************************** */

void *
msc_MemberRef (MSC * m, char *name, int *type)
{
  Ent *ne;
  void *rptr;

  if (!strcmp (name, RLAB_MEMBER_NROW))
  {
    ne = ent_Create ();
    ent_data (ne) = mdr_CreateScalar ((double) m->nr);
    ent_SetType (ne, MATRIX_DENSE_REAL);
    rptr = (void *) ne;
    *type = ENTITY;
  }
  else if (!strcmp (name, RLAB_MEMBER_NCOL))
  {
    ne = ent_Create ();
    ent_data (ne) = mdr_CreateScalar ((double) m->nc);
    ent_SetType (ne, MATRIX_DENSE_REAL);
    rptr = (void *) ne;
    *type = ENTITY;
  }
  else if (!strcmp (name, RLAB_MEMBER_SIZE))
  {
    ne = ent_Create ();
    ent_data (ne) = mdr_CreateScalar ((double) (m->nr * m->nc));
    ent_SetType (ne, MATRIX_DENSE_REAL);
    rptr = (void *) ne;
    *type = ENTITY;
  }
  else if (!strcmp (name, RLAB_MEMBER_SIZE_NON_ZERO))
  {
    ne = ent_Create ();
    ent_data (ne) = mdr_CreateScalar ((double) m->nnz);
    ent_SetType (ne, MATRIX_DENSE_REAL);
    rptr = (void *) ne;
    *type = ENTITY;
  }
  else if (!strcmp (name, RLAB_MEMBER_CLASS))
  {
    ne = ent_Create ();
    ent_data (ne) = mds_CreateScalar ( RLAB_MEMBER_CLASS_NUM );
    ent_SetType (ne, MATRIX_DENSE_STRING);
    rptr = (void *) ne;
    *type = ENTITY;
  }
  else if (!strcmp (name, RLAB_MEMBER_TYPE))
  {
    ne = ent_Create ();
    ent_data (ne) = mds_CreateScalar ( RLAB_MEMBER_TYPE_COMPLEX );
    ent_SetType (ne, MATRIX_DENSE_STRING);
    rptr = (void *) ne;
    *type = ENTITY;
  }
  else if (!strcmp (name, RLAB_MEMBER_STORAGE))
  {
    ne = ent_Create ();
    ent_data (ne) = mds_CreateScalar ( RLAB_STORAGE_SPARSE );
    ent_SetType (ne, MATRIX_DENSE_STRING);
    rptr = (void *) ne;
    *type = ENTITY;
  }
  else
  {
    /* This object does not contain one of the above.
     * Signal the caller by passing back 0.
     */
    rptr = 0;
    *type = UNDEF;
  }
  return (rptr);
}

char **
msc_Members (MSC * m, int *n)
{
  char **marray;

  marray = (char **) GC_MALLOC (7 * sizeof (char *));
  if (marray == 0)
    rerror ("out of memory");

  marray[0] = cpstr (RLAB_MEMBER_NROW);
  marray[1] = cpstr (RLAB_MEMBER_NCOL);
  marray[2] = cpstr (RLAB_MEMBER_SIZE);
  marray[3] = cpstr (RLAB_MEMBER_SIZE_NON_ZERO);
  marray[4] = cpstr (RLAB_MEMBER_CLASS);
  marray[5] = cpstr (RLAB_MEMBER_TYPE);
  marray[6] = cpstr (RLAB_MEMBER_STORAGE);

  *n = 7;
  return (marray);
}

int
msc_Sublist (MSC * m, char *name)
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
  else if (!strcmp (name, RLAB_MEMBER_SIZE_NON_ZERO))
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
  else if (!strcmp (name, RLAB_MEMBER_STORAGE))
  {
    return (1);
  }

  return (0);
}

/* **************************************************************
 * Coerce matrix indices into ints for use later...
 * Return an int array, the 1st element is the number of elements.
 * ************************************************************** */

int *
msc_IndCoerceInt (MSC * m, MDR * i)
{
  int *ind, k, size;

  if ((size = MNR (i) * MNC (i)) == 0)
  {
    return (0);
  }

  ind = (int *) GC_MAIOP ((size + 1) * sizeof (int));
  if (ind == 0)
    rerror ("out of memory");
  ind[0] = size;

  for (k = 1; k <= size; k++)
  {
    ind[k] = (int) (MdrV1 (i, k));
    if (ind[k] <= 0)
    {
      fprintf (stderr, "matrix index <= 0 not allowed\n");
      fprintf (stderr, "index value: %i\n", ind[k]);
      GC_FREE (ind);
      rerror ("matrix index coerce");
    }
  }

  return (ind);
}

/* **************************************************************
 * Get the specified value from a sparse matrix.
 * ************************************************************** */

Complex
msc_GetEl (MSC * m, int i, int j)
{
  Complex ctmp;
  int ii, rptr, nrow;

  if (i < 1 || j < 1)
  {
    fprintf (stderr, "indices i=%i and j=%i\n", i, j);
    rerror ("Internal runtime error (msc_GetEl): indices > 0 required !");
  }

  if ((i > m->nr) || (j > m->nc))
  {
    fprintf (stderr, "indices out of bounds, i=%i,   j=%i\n", i, j);
    fprintf (stderr, "                       nr=%i, nc=%i\n", m->nr, m->nc);
    rerror
      ("Internal runtime error: matrix indices must be within existing bounds");
  }

  if (m->nnz > 0)
  {
    /* Find the row pointer. */
    rptr = m->ia[i - 1];
    nrow = m->ia[i] - rptr;

    /* Now, walk the column looking for the element. */
    for (ii = 0; ii < nrow; ii++)
    {
      if (m->ja[rptr + ii - 1] == j)
        return (m->c[rptr + ii - 1]);
    }
    /* Didn't find it, must be zero. */
    ctmp = 0.0;
    return (ctmp);
  }
  else
  {
    ctmp = 0.0;
    return (ctmp);
  }

  ctmp = 0.0;
  return (ctmp);		/* Shut up compiler. */
}

/* **************************************************************
 * Assign to the elements of a sparse complex matrix.
 * This is VERY inefficient, and probably always will be.
 * This routine is also a nightmare, and probably always will be.
 * ************************************************************** */

MSC *
msc_MatrixAssign (MSC * left, int *row, int *col, MSC * rhs)
{
  int nrow, ncol, dflag;
  int i, j, k, lnr, m, nr, nc, n;
  int c, jc, *rmap, *cmap, *icol, *right_cmap;
  int *ja_tmp;
  Complex *d_tmp;
  MSC *result, *right;

  /* Error check first. */

  nrow = row[0];
  row = row + 1;
  ncol = col[0];
  col = col + 1;

  nr = nrow;
  nc = ncol;
  lnr = left->nr;

  /* Check RHS for empty matrix. */

  if ((rhs->nr == 0) || (rhs->nc == 0))
  {
    fprintf (stderr, "cannot assign empty matrix to element(s) of a matrix\n");
    rerror ("matrix-sparse-assign");
  }

  /* Check the LHS, and RHS dimensions. */

  if (nrow != rhs->nr)
  {
    fprintf (stderr, "LHS, and RHS row dimensions must match\n");
    fprintf (stderr, "LHS row: %i, RHS row: %i\n", nrow, rhs->nr);
    rerror ("matrix-sparse-assign");
  }

  if (ncol != rhs->nc)
  {
    fprintf (stderr, "LHS, and RHS column dimensions must match\n");
    fprintf (stderr, "LHS column: %i, RHS column: %i\n", ncol, rhs->nc);
    rerror ("matrix-sparse-assign");
  }

  /*
   * Check the lhs to see if the indices exceed its size.
   */

  for (i = 0; i < nrow; i++)
  {
    if (row[i] > left->nr)
    {
      fprintf (stderr, "LHS row indices exceeds matrix size\n");
      fprintf (stderr, "LHS row index: %i, matrix nr: %i\n", row[i], left->nr);
      rerror ("matrix-sparse-assign");
    }
  }

  for (i = 0; i < ncol; i++)
  {
    if (col[i] > left->nc)
    {
      fprintf (stderr, "LHS column indices exceed matrix size\n");
      fprintf (stderr, "LHS column index: %i, matrix nc: %i\n",
	       col[i], left->nc);
      rerror ("matrix-sparse-assign");
    }
  }

  /*
   * A sort of unusual circumstance, but we have to
   * behave properly in the unlikely event.
   */

  if (left == rhs)
  {
    right = msc_Copy (rhs);
    dflag = 1;
  }
  else
  {
    right = rhs;
    dflag = 0;
  }

  /*
   * rmap[i-1] is row number of `right' from which ith row
   * of `left' is assigned.  If it's zero, no assignment
   * is made from that row of `right'.
   */

#ifdef HAVE_GC
  rmap = (int *) GC_MALLOC ((left->nr) * sizeof (int));
#else
  rmap = (int *) calloc (left->nr, sizeof (int));
#endif

  for (i = 0; i < nr; i++)
  {
    rmap[row[i] - 1] = i + 1;
  }

  /*
   * cmap[i] is column number of `right' from which ith column
   * of `left' is assigned.  If it's zero, no assignment
   * is made from that column of `right'.
   */

#ifdef HAVE_GC
  cmap = (int *) GC_MALLOC ((left->nc + 1) * sizeof (int));
#else
  cmap = (int *) calloc (left->nc + 1, sizeof (int));
#endif

  for (j = 0; j < nc; j++)
  {
    cmap[col[j]] = j + 1;
  }

  /* icol[i-1] is the index of the ith nonzero in cmap */

#ifdef HAVE_GC
  icol = (int *) GC_MALLOC ((nc + 1) * sizeof (int));
#else
  icol = (int *) calloc (nc + 1, sizeof (int));
#endif

  for (n = 0, k = 1; k < left->nc + 1; k++)
  {
    if (cmap[k])
    {
      icol[n++] = k;
    }
  }
  jc = n;

  /*
   * right_cmap[j] will give the position in `right' at
   * which column j is located.  If column j of `right'
   * is not stored for row i, then right_cmap[j] will be
   * less than right->ia[i] or greater than right->ia[i+1]-1.
   */

#ifdef HAVE_GC
  right_cmap = (int *) GC_MALLOC ((right->nc + 1) * sizeof (int));
#else
  right_cmap = (int *) calloc (right->nc + 1, sizeof (int));
#endif

  /*
   * Start setting up the result.  We'll still need to decide its
   * type and symmetry, and allocate its data space.
   */

  result = msc_Create (left->nr, left->nc);

  /* Allocate as much space as it could possibly need */
  if ((result->nnz = (left->nnz + right->nnz)))
  {
    msc_Setup (result, result->nnz);
  }

  if (result->nnz)
  {
    k = n = 0;
    result->ia[0] = 1;
    for (i = 0; i < lnr; i++)
    {
      if (rmap[i])		/* assignments in this row? */
      {
	if (right->ia)
	{
	  for (m = right->ia[rmap[i] - 1]; m < right->ia[rmap[i]]; m++)
	    right_cmap[right->ja[m - 1]] = m;
	}

	if (left->ia && left->ia[i + 1] > left->ia[i])
	{
	  for (j = 0; j < jc; j++)
	  {
	    c = icol[j];
	    while (k < left->ia[i + 1] - 1 && left->ja[k] < c)
	    {
	      result->ja[n] = left->ja[k];
	      result->c[n] = left->c[k];
	      n++;
	      k++;
	    }
	    if (k < left->ia[i + 1] - 1 && left->ja[k] == c)
	      k++;
	    if (right->ia &&
		right_cmap[cmap[c]] >=
		right->ia[rmap[i] - 1] &&
		right_cmap[cmap[c]] < right->ia[rmap[i]])
	    {
	      result->ja[n] = c;
	      result->c[n] = right->c[right_cmap[cmap[c]] - 1];
	      n++;
	    }
	  }
	  while (k < left->ia[i + 1] - 1)
	  {
	    result->ja[n] = left->ja[k];
	    result->c[n] = left->c[k];
	    n++;
	    k++;
	  }
	}
	else
	{
	  for (j = 0; j < jc; j++)
	  {
	    c = icol[j];
	    if (right->ia &&
		right_cmap[cmap[c]] >=
		right->ia[rmap[i] - 1] &&
		right_cmap[cmap[c]] < right->ia[rmap[i]])
	    {
	      result->ja[n] = c;
	      result->c[n] = right->c[right_cmap[cmap[c]] - 1];
	      n++;
	    }
	  }
	}
      }
      else
      {
	if (left->ia && left->ia[i + 1] > left->ia[i])
	{
	  for (j = left->ia[i]; j < left->ia[i + 1]; j++)
	  {
	    result->ja[n] = left->ja[j - 1];
	    result->c[n] = left->c[j - 1];
	    n++;
	    k++;
	  }
	}
      }

      result->ia[i + 1] = n + 1;
    }

    result->nnz = n;

    ja_tmp = (int *) GC_MAIOP (n * sizeof (int));
    memcpy (ja_tmp, result->ja, n * sizeof (int));
    GC_FREE (result->ja);
    result->ja = ja_tmp;

    d_tmp = (Complex *) GC_MAIOP (n * sizeof (Complex));
    memcpy (d_tmp, result->c, n * sizeof (Complex));
    GC_FREE (result->c);
    result->c = d_tmp;
  }

  /* Clean up. */
  GC_FREE (rmap);
  GC_FREE (cmap);
  GC_FREE (icol);
  GC_FREE (right_cmap);
  if (dflag)
    msc_Destroy (right);
  msc_Destroy (left);

  return (result);
}

/* **************************************************************
 * Assign using vector notation. Simple fix up the indices, and
 * use msr_MatrixAssign...
 * ************************************************************** */

MSC *
msc_VectorAssign (MSC * var, int *j, MSC * rhs)
{
  int k, *row, *col;
  MSC *lhs;

  /*
   * Error check, so we can print a sensible (vector) message
   * from this point if needed.
   */

  /* Check for empty indices. */
  if (j == 0)
    return (var);

  /* Check LHS for UNDEF. */
  if (var == 0)
  {
    var = msc_Create (1, 1);
    msc_Setup (var, 0);
  }

  /* Check RHS for empty matrix. */
  if (rhs->nr * rhs->nc == 0)
  {
    fprintf (stderr,
	     "\tcannot assign empty matrix to element(s) of a matrix\n");
    rerror ("vector-assign");
  }

  /* Check the dimensions. */
  if (j[0] != (rhs->nr * rhs->nc))
  {
    fprintf (stderr, "\tLHS, and RHS column sizes must match\n");
    fprintf (stderr, "\tLHS: %i, RHS: %i\n", j[0], rhs->nr * rhs->nc);
    rerror ("vector-assign");
  }

  /*
   * Figure out row and column index integer arrays.
   */

  row = (int *) GC_MAIOP ((j[0] + 1) * sizeof (int));
  col = (int *) GC_MAIOP ((j[0] + 1) * sizeof (int));

  row[0] = j[0];
  col[0] = j[0];

  for (k = 1; k <= j[0]; k++)
  {
    col[k] = (j[k] - 1) / var->nr + 1;
    if (col[k] == 1)
    {
      row[k] = j[k];
    }
    else
    {
      row[k] = j[k] % var->nr;
      if (row[k] == 0)
	row[k] = var->nr;
    }
  }

  /* Finally, just fo it. */
  lhs = msc_MatrixAssign (var, row, col, rhs);

  GC_FREE (row);
  GC_FREE (col);

  return (lhs);
}

/* **************************************************************
 * Assign to the elements of a sparse complex matrix.
 * Nearly the same as previous function, but uses a dense RHS.
 * ************************************************************** */

extern MSC *msc_Sparse (MDC * m);

MSC *
msc_mdc_MatrixAssign (MSC * mlhs, int *irow, int *jcol, MDC * rhs)
{
  MSC *mtmp, *result;

  mtmp = msc_Sparse (rhs);
  result = msc_MatrixAssign (mlhs, irow, jcol, mtmp);
  msc_Destroy (mtmp);
  return (result);
}

MSC *
msc_msr_MatrixAssign (MSC * mlhs, int *irow, int *jcol, MSR * rhs)
{
  MSC *mtmp, *result;

  mtmp = msr_coerce_msc (rhs);
  result = msc_MatrixAssign (mlhs, irow, jcol, mtmp);
  msc_Destroy (mtmp);
  return (result);
}

MSC *
msc_mdr_MatrixAssign (MSC * mlhs, int *irow, int *jcol, MDR * rhs)
{
  MSC *mtmp, *result;
  MSR *msrtmp;

  msrtmp = msr_Sparse (rhs);
  mtmp = msr_coerce_msc (msrtmp);
  msr_Destroy (msrtmp);

  result = msc_MatrixAssign (mlhs, irow, jcol, mtmp);
  msc_Destroy (mtmp);
  return (result);
}

MSC *
msc_mdc_VectorAssign (MSC * mlhs, int *j, MDC * rhs)
{
  MSC *mtmp, *result;

  mtmp = msc_Sparse (rhs);
  result = msc_VectorAssign (mlhs, j, mtmp);
  msc_Destroy (mtmp);
  return (result);
}

MSC *
msc_msr_VectorAssign (MSC * mlhs, int *j, MSR * rhs)
{
  MSC *mtmp, *result;

  mtmp = msr_coerce_msc (rhs);
  result = msc_VectorAssign (mlhs, j, mtmp);
  msc_Destroy (mtmp);
  return (result);
}

MSC *
msc_mdr_VectorAssign (MSC * mlhs, int *j, MDR * rhs)
{
  MSC *mtmp, *result;
  MSR *dtmp;

  dtmp = msr_Sparse (rhs);
  mtmp = msr_coerce_msc (dtmp);
  result = msc_VectorAssign (mlhs, j, mtmp);
  msc_Destroy (mtmp);
  msr_Destroy (dtmp);
  return (result);
}

/* **************************************************************
 * Transpose a sparse complex matrix.
 * ************************************************************** */
extern MSC *msc_Conj (MSC * m);

MSC *
msc_Transpose (MSC * m)
{
  MSC *new, *tmp;

  tmp = msc_Conj (m);
  new = msc_Create (m->nc, m->nr);
  if (m->nnz > 0)
  {
    msc_Setup (new, m->nnz);
    ZGSTRN (tmp->ia, tmp->ja, tmp->c, &tmp->nr, &tmp->nc,
            new->ia, new->ja, new->c);
  }

  return (new);
}

MSC *
msc_NcTranspose (MSC * m)
{
  MSC *new;

  new = msc_Create (m->nc, m->nr);
  msc_Setup (new, m->nnz);
  if (m->nnz)
  {
    ZGSTRN (m->ia, m->ja, m->c, &m->nr, &m->nc, new->ia, new->ja, new->c);
  }
  return (new);
}

/* **************************************************************
 * This function pulls the selected rows out of a matrix.
 * ************************************************************** */

MSC *
msc_RowPartition (MSC * m, int *row)
{
  int i, k, kk, n, nr, nnz;
  MSC *new;

  /* Handle empty matrix indices. */

  if (row == 0)
  {
    new = msc_Create (0, 0);
    return (new);
  }

  /* row[0] is the size of row. */
  nr = row[0];
  new = msc_Create (nr, m->nc);

  /* First, count the non-zeros. */
  if (m->nnz > 0)
  {
    nnz = 0;
    for (i = 1; i <= nr; i++)
    {
      nnz = nnz + (m->ia[row[i]] - m->ia[row[i] - 1]);
    }
    new->nnz = nnz;
  }
  else
  {
    return (new);
  }

  /* If there are any non-zeros. */
  if (nnz)
  {
    /* Now allocate the new sparse matrix. */
    new->ia = (int *) GC_MAIOP ((row[0] + 1) * sizeof (int));
    new->ja = (int *) GC_MAIOP (nnz * sizeof (int));
    new->c = (Complex *) GC_MAIOP (nnz * sizeof (Complex));

    /* Now pull out the proper values. */
    k = 0;
    for (i = 1; i <= nr; i++)
    {
      new->ia[i - 1] = k + 1;
      n = m->ia[row[i]] - m->ia[row[i] - 1];
      if (n > 0)
      {
	kk = m->ia[row[i] - 1] - 1;
	memcpy (new->ja + k, m->ja + kk, n * sizeof (int));
	memcpy (new->c + k, m->c + kk, n * sizeof (Complex));
      }
      k += n;
    }
    new->ia[i - 1] = k + 1;
  }

  return (new);
}

/* **************************************************************
 * Select rows and columns from a sparse-real matrix.
 * ************************************************************** */

void *
msc_MatrixSub (MSC * m, int *row, int *col, int *type)
{
  MSC *new, *mtmp1, *mtmp2, *mtmp3;
  MDC *dnew;

  /* Handle empty matrix indices. */

  if ((row == 0) || (col == 0))
  {
    new = msc_Create (0, 0);
    *type = MATRIX_SPARSE_COMPLEX;
    return (new);
  }

  /* Handle scalar indices. */

  if ((row[0] == 1) && (col[0] == 1))
  {
    Complex c = msc_GetEl (m, row[1], col[1]);
    dnew = mdc_CreateScalar (RE(c), IM(c));
    *type = MATRIX_DENSE_COMPLEX;
    return (dnew);
  }

  /*
   * Now things get a little unusual. Instead of
   * just pulling out the referenced values, we
   * are going to pull out the specified rows,
   * transpose the result, then pullout the
   * specified row (the original columns), transpose
   * the result back, and we have what we were after
   * in the first place.
   */

  mtmp1 = msc_RowPartition (m, row);

  /* Now transpose. */

  mtmp2 = msc_NcTranspose (mtmp1);
  msc_Destroy (mtmp1);

  /* Now select columns (as rows). */
  mtmp3 = msc_RowPartition (mtmp2, col);
  msc_Destroy (mtmp2);

  /* Now transpose back. */
  new = msc_NcTranspose (mtmp3);
  msc_Destroy (mtmp3);

  *type = MATRIX_SPARSE_COMPLEX;
  return (new);
}

/* **************************************************************
 * Select rows from a sparse-complex matrix.
 * ************************************************************** */

void *
msc_MatrixSubR (MSC * m, int *row, int *type)
{
  MSC *new;
  MDC *dnew;

  /* Handle empty matrix indices. */

  if (row == 0)
  {
    new = msc_Create (0, 0);
    *type = MATRIX_SPARSE_COMPLEX;
    return (new);
  }

  /* Handle scalar indices. */

  if ((row[0] == 1) && (m->nc == 1))
  {
    Complex c = msc_GetEl (m, row[1], 1);
    dnew = mdc_CreateScalar (RE(c), IM(c));
    *type = MATRIX_DENSE_COMPLEX;
    return (dnew);
  }

  *type = MATRIX_SPARSE_COMPLEX;
  new = msc_RowPartition (m, row);
  return (new);
}

/* **************************************************************
 * Select columns from a sparse-complex matrix.
 * ************************************************************** */

void *
msc_MatrixSubC (MSC * m, int *col, int *type)
{
  MSC *new, *mtmp1, *mtmp2;
  MDC *dnew;

  /* Handle empty matrix indices. */

  if (col == 0)
  {
    new = msc_Create (0, 0);
    *type = MATRIX_SPARSE_COMPLEX;
    return (new);
  }

  /* Handle scalar indices. */

  if ((col[0] == 1) && (m->nr == 1))
  {
    Complex c = msc_GetEl (m, 1, col[1]);
    dnew = mdc_CreateScalar (RE(c), IM(c));
    *type = MATRIX_DENSE_COMPLEX;
    return (dnew);
  }

  /* Transpose. */
  mtmp1 = msc_NcTranspose (m);

  /* Now select columns (as rows). */
  mtmp2 = msc_RowPartition (mtmp1, col);
  msc_Destroy (mtmp1);

  /* Now transpose back. */
  new = msc_NcTranspose (mtmp2);
  msc_Destroy (mtmp2);

  *type = MATRIX_SPARSE_COMPLEX;
  return (new);
}

/* **************************************************************
 * Sub-matrix expression using vector indices.
 * ************************************************************** */

void *
msc_VectorSub (MSC * m, int *i, int *type)
{
  int coln, j, rown, size, msize;
  MSC *new = 0;
  MSC *tmp = 0;
  MDC *dnew = 0;

  /* Handle empty matrix indices. */

  if (i == 0)
  {
    new = msc_Create (0, 0);
    return (new);
  }

  size = i[0];
  msize = m->nr * m->nc;

  /* If size == 1, just do it... */
  if (size == 1)
  {
    Complex c;

    /* Convert from 1D -> 2D */
    coln = (i[1] - 1) / m->nr + 1;
    if (coln == 1)
    {
      rown = i[1];
    }
    else
    {
      rown = i[1] % m->nr;
      if (rown == 0)
	rown = m->nr;
    }

    if (rown > m->nr || coln > m->nc)
    {
      fprintf (stderr, "\tindex exceeds matrix limits\n");
      fprintf (stderr, "\tindex value: %i, matrix size: %i\n", i[1], msize);
      rerror ("sub-matrix evaluation");
    }

    c = msc_GetEl (m, rown, coln);
    dnew = mdc_CreateScalar (RE(c), IM(c));

    *type = MATRIX_DENSE_COMPLEX;
    return (dnew);
  }

  /*
   * Now get on with the general case.
   * This is a real brute-force approach.
   * If this is not efficient enough
   * (i.e. someone compplains), it can certainly
   * be made more efficient.
   */

  new = msc_Create (1, size);	/* Most common. */
  msc_Setup (new, size);

  for (j = 0; j < size; j++)
  {
    if (i[j + 1] > msize)
    {
      msc_Destroy (new);
      fprintf (stderr, "\tindex exceeds matrix limits\n");
      fprintf (stderr, "\tindex value: %i, matrix size: %i\n", i[j + 1], msize);
      rerror ("sub-matrix evaluation");
    }

    /* Convert from 1D -> 2D */
    coln = (i[j + 1] - 1) / m->nr + 1;
    if (coln == 1)
    {
      rown = i[j + 1];
    }
    else
    {
      rown = i[j + 1] % m->nr;
      if (rown == 0)
	rown = m->nr;
    }

    new->ja[j] = j + 1;
    new->c[j] = msc_GetEl (m, rown, coln);
  }

  new->ia[0] = 1;
  new->ia[1] = size + 1;;

  /*
   * Transpose the result if necessary to
   * match Matlab's convention. Also,
   * resparse to get rid of unwanted zeros.
   */

  if (m->nc == 1)
  {
    tmp = msc_NcTranspose (new);
    msc_Destroy (new);
    new = msc_ReSparse (tmp);
    msc_Destroy (tmp);
  }
  else
  {
    tmp = msc_ReSparse (new);
    msc_Destroy (new);
    new = tmp;
  }

  *type = MATRIX_SPARSE_COMPLEX;
  return (new);
}

/* **************************************************************
 * Stack two sparse matrices together:  [ m1 ; m2 ]
 * ************************************************************** */

MSC *
msc_Stack (MSC * m1, MSC * m2)
{
  int i, j, k, nrow, ncol, nnz;
  MSC *new;

  /* Check for empty matrices... */
  if (m1->nr == 0 && m1->nc == 0)
  {
    new = msc_Copy (m2);
    return (new);
  }
  else if (m2->nr == 0 && m2->nc == 0)
  {
    new = msc_Copy (m1);
    return (new);
  }

  /* Check the sizes for correctness. */

  if (m1->nc != m2->nc)
  {
    fprintf (stderr, "matrix stack requires equal column-sizes\n");
    fprintf (stderr, "existing column sizes: %i and %i\n", m1->nc, m2->nc);
    rerror ("sparse matrix stack");
  }

  /* Do the stack. */
  nrow = m1->nr + m2->nr;
  ncol = m1->nc;
  new = msc_Create (nrow, ncol);

  nnz = m1->nnz + m2->nnz;

  if (nnz > 0)
  {
    msc_Setup (new, nnz);

    memcpy (new->ia, m1->ia, (m1->nr + 1) * sizeof (int));
    memcpy (new->ja, m1->ja, (m1->nnz) * sizeof (int));
    memcpy (new->c, m1->c, (m1->nnz) * sizeof (Complex));

    new->ia[m1->nr + 1] = new->ia[m1->nr] + (m2->ia[1] - m2->ia[0]);
    for (i = m1->nr + 2; i < nrow + 1; i++)
    {
      j = i - (m1->nr + 1);
      k = m2->ia[j + 1] - m2->ia[j];
      new->ia[i] = new->ia[i - 1] + k;
    }

    memcpy (new->ja + m1->nnz, m2->ja, (m2->nnz) * sizeof (int));
    memcpy (new->c + m1->nnz, m2->c, (m2->nnz) * sizeof (Complex));
  }

  return (new);
}

/* **************************************************************
 * Append two sparse-complex matrices together:  [ m1 , m2 ]
 * ************************************************************** */

MSC *
msc_Append (MSC * m1, MSC * m2)
{
  MSC *new, *mtmp1, *mtmp2, *mtmp3;

  /* Check for empty matrices... */
  if (m1->nr == 0 && m1->nc == 0)
  {
    new = msc_Copy (m2);
    return (new);
  }
  else if (m2->nr == 0 && m2->nc == 0)
  {
    new = msc_Copy (m1);
    return (new);
  }

  /* Check the sizes for correctness. */

  if (m1->nr != m2->nr)
  {
    fprintf (stderr, "matrix append requires equal row-sizes\n");
    fprintf (stderr, "existing row sizes: %i and %i\n", m1->nr, m2->nr);
    rerror ("sparse matrix append");
  }

  /* Do the append. */

  /* This is a real hack for now. */
  mtmp1 = msc_Transpose (m1);
  mtmp2 = msc_Transpose (m2);
  mtmp3 = msc_Stack (mtmp1, mtmp2);
  msc_Destroy (mtmp1);
  msc_Destroy (mtmp2);

  new = msc_Transpose (mtmp3);
  msc_Destroy (mtmp3);

  return (new);
}

/* **************************************************************
 * Negate a sparse complex matrix.
 * ************************************************************** */

MSC *
msc_Negate (MSC * m)
{
  int i, size;
  MSC *mnew = 0;

  mnew = msc_Copy (m);
  size = m->nnz;
  for (i = 0; i < size; i++)
  {
    mnew->c[i] = -(mnew->c[i]);
  }
  return (mnew);
}

/* **************************************************************
 * Compare two matrices (m1 == m2).
 * All comparisons are performed assuming there are no explicit
 * zero elements in the structure.
 * ************************************************************** */

MSR *
msc_Eq (MSC * m1, MSC * m2)
{
  int i, ia, j, n, size;
  MDR *dtmp;
  MSR *tmp = 0;
  MSR *new = 0;

  /* Special cases... */

  if (m1->nr == 0 || m2->nr == 0)
  {
    /* Empty matrices. */
    if (m1->nr == 0 && m2->nr == 0)
    {
      new = msr_CreateScalar (1.0);
      return (new);
    }
    else
    {
      new = msr_CreateScalar (0.0);
      return (new);
    }
  }
  else if ((m1->nr == 1) && (m1->nc == 1))
  {
    /* m1 is scalar. */
    if (m1->nnz)
    {
      /*
       * Result has the same sparsity pattern.
       * Worst case... re-sparse later.
       */

      size = m2->nnz;
      tmp = msr_Create (m2->nr, m2->nc);
      msr_Setup (tmp, size);
      memcpy (tmp->ia, m2->ia, (m2->nr + 1) * sizeof (int));
      memcpy (tmp->ja, m2->ja, (m2->nnz) * sizeof (int));

      for (i = 0; i < size; i++)
      {
        tmp->d[i] = (double) ((RE(m1->c[0]) == RE(m2->c[i])) &&
            (IM(m1->c[0]) == IM(m2->c[i])));
      }

      /* Now, re-sparse. */

      new = msr_ReSparse (tmp);
      msr_Destroy (tmp);
    }
    else
    {
      /*
       * Result has the inverse sparsity pattern.
       */

      new = msc_PatternInverse (m2);
      size = new->nnz;
      for (i = 0; i < size; i++)
      {
        new->d[i] = 1.0;
      }
    }
  }
  else if ((MNR (m2) == 1) && (MNC (m2) == 1))
  {
    /* m2 is scalar. */
    if (m2->nnz)
    {
      /*
       * Result has the same sparsity pattern.
       * Worst case... re-sparse later.
       */

      size = m1->nnz;
      tmp = msr_Create (m1->nr, m1->nc);
      msr_Setup (tmp, size);
      memcpy (tmp->ia, m1->ia, (m1->nr + 1) * sizeof (int));
      memcpy (tmp->ja, m1->ja, (m1->nnz) * sizeof (int));

      for (i = 0; i < size; i++)
      {
        tmp->d[i] = (double) ((RE(m2->c[0]) == RE(m1->c[i])) &&
            (IM(m2->c[0]) == IM(m1->c[i])));
      }

      /* Now, re-sparse. */

      new = msr_ReSparse (tmp);
      msr_Destroy (tmp);
    }
    else
    {
      /*
       * Result has the inverse sparsity pattern.
       */

      new = msc_PatternInverse (m1);
      size = new->nnz;
      for (i = 0; i < size; i++)
      {
        new->d[i] = 1.0;
      }
    }
  }
  else if ((MNR (m1) != MNR (m2)) || (MNC (m1) != MNC (m2)))
  {
    fprintf (stderr, "matrices must have the same dimension for comparision\n");
    fprintf (stderr, "m1: nrow = %i ncol = %i m2: nrow = %i, ncol = %i\n",
             MNR (m1), MNC (m1), MNR (m2), MNC (m2));
    rerror ("matrix-compare");
  }
  else
  {
    /*
    * Both matrices are the same dimension.
    * This is a lazy and in-efficient approach,
    * until I write msr_PatternUnion
    */

    Complex c1, c2;

    tmp = msr_Create (m1->nr, m1->nc);
    msc_PatternUnion (m1, m2, &(tmp->ia), &(tmp->ja), &(tmp->nnz));

    /* The result... make it all ones to start with, since
    * we know most of the result will be 1. */
    dtmp = mdr_Create (m1->nr, m1->nc);
    mdr_Ones (dtmp);

    for (i = 0; i < tmp->nr; i++)
    {
      n = tmp->ia[i + 1] - tmp->ia[i];
      ia = tmp->ia[i];

      for (j = 0; j < n; j++)
      {
        c1 = msc_GetEl (m1, i + 1, tmp->ja[ia + j - 1]);
        c2 = msc_GetEl (m2, i + 1, tmp->ja[ia + j - 1]);
        Mdr1 (dtmp, i + 1, tmp->ja[ia + j - 1]) =
            (double) ((RE(c1) == RE(c2)) && (IM(c1) == IM(c2)));
      }
    }
    new = msr_Sparse (dtmp);
    mdr_Destroy (dtmp);
    msr_Destroy (tmp);
  }

  return (new);
}

MSR *
    mdc_msc_Eq (MDC * m1, MSC * m2)
{
  MSC *mtmp;
  MSR *mret;

  mtmp = msc_Sparse (m1);
  mret = msc_Eq (mtmp, m2);
  msc_Destroy (mtmp);

  return (mret);
}

MSR *
    msc_mdc_Eq (MSC * m1, MDC * m2)
{
  MSC *mtmp;
  MSR *mret;

  mtmp = msc_Sparse (m2);
  mret = msc_Eq (m1, mtmp);
  msc_Destroy (mtmp);

  return (mret);
}

MSR *
    msc_mdr_Eq (MSC * m1, MDR * m2)
{
  MDC *mdc_tmp;
  MSC *mtmp;
  MSR *mret;

  mdc_tmp = mdr_coerce_mdc (m2);
  mtmp = msc_Sparse (mdc_tmp);
  mret = msc_Eq (m1, mtmp);
  msc_Destroy (mtmp);
  mdc_Destroy (mdc_tmp);

  return (mret);
}

MSR *
    mdr_msc_Eq (MDR * m1, MSC * m2)
{
  MDC *mdc_tmp;
  MSC *mtmp;
  MSR *mret;

  mdc_tmp = mdr_coerce_mdc (m1);
  mtmp = msc_Sparse (mdc_tmp);
  mret = msc_Eq (mtmp, m2);
  msc_Destroy (mtmp);
  mdc_Destroy (mdc_tmp);

  return (mret);
}

MSR *
msc_msr_Eq (MSC * m1, MSR * m2)
{
  MSC *mtmp;
  MSR *mret;

  mtmp = msr_coerce_msc (m2);
  mret = msc_Eq (m1, mtmp);
  msc_Destroy (mtmp);

  return (mret);
}

MSR *
msr_msc_Eq (MSR * m1, MSC * m2)
{
  MSC *mtmp;
  MSR *mret;

  mtmp = msr_coerce_msc (m1);
  mret = msc_Eq (mtmp, m2);
  msc_Destroy (mtmp);

  return (mret);
}

/* **************************************************************
 * Compare two matrices (m1 != m2).
 * All comparisons are performed assuming there are no explicit
 * zero elements in the structure.
 * ************************************************************** */

MSR *
msc_Ne (MSC * m1, MSC * m2)
{
  int i, ia, j, n;
  MDR *dtmp;
  MSR *tmp = 0;
  MSR *new = 0;

  /* Special cases... */

  if (m1->nr == 0 || m2->nr == 0)
  {
    if (m1->nr == 0 && m2->nr == 0)
    {
      new = msr_CreateScalar (0.0);
      return (new);
    }
    else
    {
      new = msr_CreateScalar (1.0);
      return (new);
    }
  }
  else if ((m1->nr == 1) && (m1->nc == 1))
  {
    if (m1->nnz)		/* m1 scalar, and != 0 */
    {
      /*
       * Result has undefined sparsity pattern !
       * Start dense, and sparse later.
       */

//       size = m2->nr * m2->nc;

      /* The result... make it all ones to start with, since
       * we know most of the result will be 1. */
      dtmp = mdr_Create (m2->nr, m2->nc);
      mdr_Ones (dtmp);

      /* Now, walk through the non-zeros of the sparse matrix. */
      for (i = 0; i < m2->nr; i++)
      {
        n = m2->ia[i + 1] - m2->ia[i];
        ia = m2->ia[i];

        for (j = 0; j < n; j++)
        {
          Complex ctmp = msc_GetEl (m2, i + 1, m2->ja[ia + j - 1]);
          Mdr1 (dtmp, i + 1, m2->ja[ia + j - 1]) =
              (double) ((RE(ctmp) != RE(m1->c[0])) || (IM(ctmp) != IM(m1->c[0])));
        }
      }
      new = msr_Sparse (dtmp);
      mdr_Destroy (dtmp);
    }
    else      /* m1 scalar, zero */
    {
      /*
       * Result has same sparsity pattern.
       */

      new = msr_Create (m1->nr, m1->nc);
      new->nnz = m1->nnz;
      /* Copy the IA array. */
      if (m1->ia)
      {
	new->ia = (int *) GC_MAIOP (sizeof (int) * (m1->nr + 1));
	memcpy (new->ia, m1->ia, (m1->nr + 1) * sizeof (int));
      }
      /* Copy the JA array. */
      if (m1->ja)
      {
	new->ja = (int *) GC_MAIOP (sizeof (int) * (m1->nnz));
	memcpy (new->ja, m1->ja, m1->nnz * sizeof (int));
      }
      new->order = m1->order;

      for (i = 0; i < new->nnz; i++)
      {
	new->d[i] = 1.0;
      }
    }
  }
  else if ((m2->nr == 1) && (m2->nc == 1))
  {
    if (m2->nnz)		/* m2 scalar, and != 0 */
    {
      /*
       * Result has undefined sparsity pattern !
       * Start dense, and sparse later.
       */

//       size = m1->nr * m1->nc;

      /* The result... make it all ones to start with, since
       * we know most of the result will be 1. */
      dtmp = mdr_Create (m1->nr, m1->nc);
      mdr_Ones (dtmp);

      /* Now, walk through the non-zeros of the sparse matrix. */
      for (i = 0; i < m1->nr; i++)
      {
	n = m1->ia[i + 1] - m1->ia[i];
	ia = m1->ia[i];

	for (j = 0; j < n; j++)
	{
	  Complex ctmp = msc_GetEl (m1, i + 1, m1->ja[ia + j - 1]);
	  Mdr1 (dtmp, i + 1, m1->ja[ia + j - 1]) =
	    (double) ((RE(ctmp) != RE(m2->c[0])) || (IM(ctmp) != IM(m2->c[0])));
	}
      }
      new = msr_Sparse (dtmp);
      mdr_Destroy (dtmp);
    }
    else
    {
      /*
       * Result has inverse sparsity pattern.
       */

      new = msr_Create (m1->nr, m1->nc);
      new->nnz = m1->nnz;

      /* Copy the IA array. */
      if (m1->ia)
      {
	new->ia = (int *) GC_MAIOP (sizeof (int) * (m1->nr + 1));
	memcpy (new->ia, m1->ia, (m1->nr + 1) * sizeof (int));
      }

      /* Copy the JA array. */
      if (m1->ja)
      {
	new->ja = (int *) GC_MAIOP (sizeof (int) * (m1->nnz));
	memcpy (new->ja, m1->ja, m1->nnz * sizeof (int));
      }

      new->order = m1->order;

      for (i = 0; i < new->nnz; i++)
      {
	new->d[i] = 1.0;
      }
    }
  }
  else if ((MNR (m1) != MNR (m2)) || (MNC (m1) != MNC (m2)))
  {
    fprintf (stderr, "matrices must have the same dimension for comparision\n");
    fprintf (stderr, "m1: nrow = %i ncol = %i m2: nrow = %i, ncol = %i\n",
	     MNR (m1), MNC (m1), MNR (m2), MNC (m2));
    rerror ("matrix-compare");
  }
  else
  {
    /*
     * Both matrices are the same dimension.
     */

//     size = m1->nr * m1->nc;
    tmp = msr_Create (m1->nr, m1->nc);
    msc_PatternUnion (m1, m2, &(tmp->ia), &(tmp->ja), &(tmp->nnz));
    tmp->d = (double *) GC_MAIOP (tmp->nnz * sizeof (double));

    /* Now, walk through the new (union) sparse matrix. */
    for (i = 0; i < tmp->nr; i++)
    {
      n = tmp->ia[i + 1] - tmp->ia[i];
      ia = tmp->ia[i];

      for (j = 0; j < n; j++)
      {
	Complex ctmp1 = msc_GetEl (m1, i + 1, tmp->ja[ia + j - 1]);
	Complex ctmp2 = msc_GetEl (m2, i + 1, tmp->ja[ia + j - 1]);

	tmp->d[ia + j - 1] = (double) ((RE(ctmp1) != RE(ctmp2))
				       || (IM(ctmp1) != IM(ctmp2)));
      }
    }
    new = msr_ReSparse (tmp);
    msr_Destroy (tmp);
  }

  return (new);
}

MSR *
mdc_msc_Ne (MDC * m1, MSC * m2)
{
  MSC *mtmp;
  MSR *mret;

  mtmp = msc_Sparse (m1);
  mret = msc_Ne (mtmp, m2);
  msc_Destroy (mtmp);

  return (mret);
}

MSR *
msc_mdc_Ne (MSC * m1, MDC * m2)
{
  MSC *mtmp;
  MSR *mret;

  mtmp = msc_Sparse (m2);
  mret = msc_Ne (m1, mtmp);
  msc_Destroy (mtmp);

  return (mret);
}

/* **************************************************************
 * Create a sparse scalar...
 * ************************************************************** */

MSC *
msc_CreateScalar (double valr, double vali)
{
  MSC *new;

  new = msc_Create (1, 1);
  msc_Setup (new, 1);

  new->c[0] = valr + vali * 1i;
  new->ja[0] = 1;
  new->ia[0] = 1;
  new->ia[1] = 2;

  return (new);
}

/* **************************************************************
 * Create a sparse matrix with the inverse sparsity patter...
 * non-zeros where there used to be zeros, and zeros where there
 * used to be non-zeros.
 *
 * Here, we return a Real Sparse Matrix !
 * ************************************************************** */

MSR *
msc_PatternInverse (MSC * m)
{
  int i, iadiff, iiadiff, j, k, l, n;
  int iaptr, iiaptr;
  MSR *mi;

  n = m->nr * m->nc;
  mi = msr_Create (m->nr, m->nc);
  msr_Setup (mi, n - m->nnz);

  /* First, compute IA. */
  mi->ia[0] = 1;
  for (i = 1; i < m->nr + 1; i++)
  {
    iadiff = m->ia[i] - m->ia[i - 1];
    mi->ia[i] = mi->ia[i - 1] + (m->nc - iadiff);
  }

  /* Now compute JA. */
  /* Now go through each row, and set the inverse pattern. */
  for (i = 0; i < m->nr; i++)
  {
    iadiff = m->ia[i + 1] - m->ia[i];	/* Length of old row. */
    iiadiff = mi->ia[i + 1] - mi->ia[i];	/* Length of new row. */

    iaptr = m->ia[i] - 1;
    iiaptr = mi->ia[i] - 1;

    if (iiadiff)
    {
      k = 0;
      l = 0;
      for (j = 0; j < m->nc; j++)
      {
	if (m->ja[iaptr + k] == j + 1)
	{
	  if (k + 1 < iadiff)
	    k++;		/* Skip to the next ja. */
	}
	else
	{
	  mi->ja[iiaptr + l++] = j + 1;
	}
      }
    }
  }

  return (mi);
}

/* **************************************************************
 * Create the IA and JA vectors for a sparse matrix with a
 * "pattern" of non-zeros that is the union of two sparse matrices.
 * ************************************************************** */

void
msc_PatternUnion (MSC * m1, MSC * m2, int **IA, int **JA, int *nnz)
{
  int i, j, ia1, ia2, *ja1, *ja2, n1, n2, nset;
  int *ia, *ja;
  int *janew, *jatmp;

  /*
   * Special cases...
   */

  if (m1->nnz == 0 || m2->nnz == 0)
  {
    /* Just return the IA and JA from the non-zero matrix. */
    if (m1->nnz != 0)
    {
      ia = (int *) GC_MAIOP ((m1->nr + 1) * sizeof (int));
      memcpy (ia, m1->ia, (m1->nr + 1) * sizeof (int));

      ja = (int *) GC_MAIOP ((m1->nnz) * sizeof (int));
      memcpy (ja, m1->ja, (m1->nnz) * sizeof (int));

      *IA = ia;
      *JA = ja;
      *nnz = m1->nnz;
      return;
    }
    else if (m2->nnz != 0)
    {
      ia = (int *) GC_MAIOP ((m2->nr + 1) * sizeof (int));
      memcpy (ia, m2->ia, (m2->nr + 1) * sizeof (int));

      ja = (int *) GC_MAIOP ((m2->nnz) * sizeof (int));
      memcpy (ja, m2->ja, (m2->nnz) * sizeof (int));

      *IA = ia;
      *JA = ja;
      *nnz = m2->nnz;
      return;
    }
  }

  /*
   * Normal behavior...
   * Go through row by row, and figure out
   * the number of elements in each row.
   */

  /* Guess (conservatively) at the new size. */
  ia = (int *) GC_MAIOP ((m1->nr + 1) * sizeof (int));
  jatmp = (int *) GC_MAIOP ((m1->nnz + m2->nnz) * sizeof (int));

  /* Row-by-row */
  ia[0] = 1;
  janew = 0;
  for (i = 0; i < m1->nr; i++)
  {
    /*
    * n? is the number of elements in each row.
    * ia? is the pointer to the beginning of each row.
    */

    n1 = m1->ia[i + 1] - m1->ia[i];
    ia1 = m1->ia[i];
    ja1 = &(m1->ja[ia1 - 1]);

    n2 = m2->ia[i + 1] - m2->ia[i];
    ia2 = m2->ia[i];
    ja2 = &(m2->ja[ia2 - 1]);

    /* Get the union of the elements in this row. */
    if (janew)
      GC_FREE (janew);

    /* Check for special cases. */
    if (n1 == 0 && n2 != 0)
    {
      nset = n2;
      janew = (int *) GC_MAIOP (n2 * sizeof (int));
      memcpy (janew, ja2, n2 * sizeof (int));
    }
    else if (n2 == 0 && n1 != 0)
    {
      nset = n1;
      janew = (int *) GC_MAIOP (n1 * sizeof (int));
      memcpy (janew, ja1, n1 * sizeof (int));
    }
    else if (n1 == 0 && n2 == 0)
    {
      nset = 0;
      janew = 0;
    }
    else
    {
      janew = array_union (ja1, n1, ja2, n2, &nset);
    }

    /* Set the value of ia[]. */
    ia[i + 1] = ia[i] + nset;

    /* Copy janew into jatmp[]. */
    for (j = 0; j < nset; j++)
    {
      jatmp[ia[i] + j - 1] = janew[j];
    }
  }
  if (janew)
    GC_FREE (janew);		/* Don't forget this! */

  /* All done, re-size ia[] and ja[]. */
  *nnz = ia[m1->nr] - 1;
  ja = (int *) GC_MAIOP (*nnz * sizeof (int));
  memcpy (ja, jatmp, *nnz * sizeof (int));

  *IA = ia;
  *JA = ja;

  GC_FREE (jatmp);
  return;
}

/* **************************************************************
 * Coerce a complex sparse matrix to real sparse
 * ************************************************************** */

MSR *
msc_coerce_msr (MSC * m)
{
  int i, size;

  MSR *msr = msr_Create (m->nr, m->nc);
  msr_Setup (msr, m->nnz);

  if (m->nnz)
  {
    memcpy (msr->ia, m->ia, (m->nr + 1) * sizeof (int));
    memcpy (msr->ja, m->ja, (m->nnz) * sizeof (int));

    size = m->nnz;
    for (i = 0; i < size; i++)
    {
      msr->d[i] = (double) RE(m->c[i]);
    }
  }
  return (msr);
}

/* **************************************************************
 * Coerce a real sparse matrix to complex sparse
 * ************************************************************** */

MSC *
msr_coerce_msc (MSR * m)
{
  int i, size;

  MSC *msc = msc_Create (m->nr, m->nc);
  msc_Setup (msc, m->nnz);

  if (m->nnz)
  {
    memcpy (msc->ia, m->ia, (m->nr + 1) * sizeof (int));
    memcpy (msc->ja, m->ja, (m->nnz) * sizeof (int));

    size = m->nnz;
    for (i = 0; i < size; i++)
    {
      msc->c[i] = (double) RE(m->d[i]);
    }
  }
  return (msc);
}

/* **************************************************************
 * Check a sparse matrix for consistency. Cause a core-dump
 * if anything is wrong.
 * ************************************************************** */
#include <signal.h>

void
msc_Check (MSC * m)
{
  int i;

  if ((m->nr * m->nc) > m->nnz)
  {
    fprintf (stderr, "sparse-matrix-handling: nr * nc > nnz !\n");
    fprintf (stderr, "sparse-matrix-handling: m->nr=%i, m->nc=%i, m->nnz=%i\n",
	     m->nr, m->nc, m->nnz);
    rerror ("terrible error in sparse matrix check (msc_Check)");
  }

  if (m->nnz)
  {
    for (i = 0; i < m->nnz; i++)
    {
      if ((m->ja[i] > m->nc) || (m->ja[i] < 1))
      {
	fprintf (stderr, "sparse-matrix-handling: m->ja[i] > m->nc\n");
	fprintf (stderr, "sparse-matrix-handling: m->ja[i] < 1\n");
	fprintf (stderr,
		 "sparse-matrix-handling: m->nr=%i, m->nc=%i, m->nnz=%i, i=%i, m->ja[i]=%i\n",
		 m->nr, m->nc, m->nnz, i, m->ja[i]);
	rerror ("Internal runtime error (msc_Check)");
      }
    }

    for (i = 0; i < m->nr; i++)
    {
      if ((m->ia[i] > m->nnz + 1) || (m->ia[i] < 1))
      {
	fprintf (stderr, "sparse-matrix-handling: m->ia[i] > m->nnz+1\n");
	fprintf (stderr, "sparse-matrix-handling: m->ia[i] < 1\n");
	fprintf (stderr,
		 "sparse-matrix-handling: m->nr=%i, m->nc=%i, m->nnz=%i, i=%i, m->ia[i]=%i\n",
		 m->nr, m->nc, m->nnz, i, m->ia[i]);
	rerror ("Internal runtime error (msc_Check)");
      }
    }
  }
}
