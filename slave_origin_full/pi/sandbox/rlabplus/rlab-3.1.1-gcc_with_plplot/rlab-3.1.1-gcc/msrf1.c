/* msrf1.c Matrix Sparse Real Functions ... */

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
#include "mdr.h"
#include "mem.h"
#include "msr.h"
#include "msc.h"
#include "mscf1.h"
#include "mdrf1.h"
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
 * Convert a dense real matrix to a sparse real matrix.
 * ************************************************************** */

MSR *
msr_Sparse (MDR * m)
{
  MSR *sm;
  int i, j, k, nnz, nc, nr, size;

  nnz = 0;
  nr = MNR (m);
  nc = MNC (m);
  sm = msr_Create (nr, nc);

  /* Count the number of non-zeros */
  size = nr * nc;
  for (i = 0; i < size; i++)
  {
    if (MdrV0 (m, i) != 0.0)
    {
      nnz++;
    }
  }
  sm->nnz = nnz;

  /* Now allocate the sparse storage arrays. */
  sm->d = (double *) GC_MAIOP (nnz * sizeof (double));
  sm->ja = (int *) GC_MAIOP (nnz * sizeof (int));
  sm->ia = (int *) GC_MAIOP ((nr + 1) * sizeof (int));
  sm->ia[0] = 1;
  k = 1;

  /* Fill the sparse storage arrays with the proper values. */
  for (i = 0; i < nr; i++)
  {
    for (j = 0; j < nc; j++)
    {
      if (Mdr0 (m, i, j) != 0.0)
      {
	sm->d[k - 1] = Mdr0 (m, i, j);
	sm->ja[k - 1] = j + 1;
	k++;
      }
    }
    sm->ia[i + 1] = k;
  }

  sm->order = 1;
  /* msr_Check (sm); *//* Get rid of this soon. */
  return (sm);
}

/* **************************************************************
 * Convert a sparse real matrix to a dense real matrix.
 * ************************************************************** */

MDR *
msr_Dense (MSR * m)
{
  MDR *dm;
  int i, ia, j, n, nr, nc, row;

  nr = m->nr;
  nc = m->nc;
  dm = mdr_Create (nr, nc);
  mdr_Zero (dm);

  if (m->nnz)
  {
    row = 1;
    for (i = 0; i < nr; i++)
    {
      n = m->ia[i + 1] - m->ia[i];
      ia = m->ia[i];
      for (j = 0; j < n; j++)
      {
	Mdr1 (dm, row, m->ja[ia + j - 1]) = m->d[ia + j - 1];
      }
      row++;
    }
  }
  return (dm);
}

/* **************************************************************
 * Convert a dense, 3-column matrix into a sparse matrix.
 * 1st column:  Row id
 * 2nd column:  Column id
 * 3rd column:  matrix elements (non-zeros)
 * 4th column:  (optional) imaginary numbers
 *
 * The input matrix does not have to be row-sorted, but duplicate
 * entries are not supported. I don't even know how to check for
 * duplicate entries yet.
 * ************************************************************** */

MSR *mdr_real_spconvert (MDR * m);
MSC *mdr_complex_spconvert (MDR * m);

void *
mdr_Spconvert (MDR * m, int *rtype)
{
  void *sm;
  int do_cmplx = 0;

  if (MNC (m) == 3)
    do_cmplx = 0;
  else if (MNC (m) == 4)
    do_cmplx = 1;
  else
    rerror ("spconvert: input must be 3 or 4 columns");

  if (!do_cmplx)
  {
    sm = (void *) mdr_real_spconvert (m);
    *rtype = MATRIX_SPARSE_REAL;
  }
  else
  {
    sm = (void *) mdr_complex_spconvert (m);
    *rtype = MATRIX_SPARSE_COMPLEX;
  }

  return (sm);
}

#define tsy
#ifdef tsy

/*
 * T.S. Yang, 3/97
 * use sorting for spconvert, instead of double transpose
 * pros: low memory overhead, short code, easy to write
 * cons: slow?, maybe, will know soon
*/

static int *q_ind, *q_r, *q_c;
static int q_cutoff = 16;	/* qsort cutoff threshold */

static void
indexswap (int i, int j)
{
  int tmp;
  tmp = q_ind[i];		/* swap indices */
  q_ind[i] = q_ind[j];
  q_ind[j] = tmp;
}

static void
insertsort (int left, int right)
{
  int i, j, k, i1, j1, i2, j2, t;

  /* insertion sort */
  for (i = left; i < right; i++)
  {
    i2 = q_r[q_ind[i]];
    j2 = q_c[q_ind[i]];
    k = i;
    for (j = i + 1; j <= right; j++)
    {
      i1 = q_r[q_ind[j]];
      j1 = q_c[q_ind[j]];
      /* double keys */
      if (i1 < i2)
      {
	i2 = i1;
	j2 = j1;
	k = j;
      }
      else if (i1 == i2)
	if (j1 < j2)
	{
	  i2 = i1;
	  j2 = j1;
	  k = j;
	}
    }
    if (k != i)
    {
      t = q_ind[i];
      q_ind[i] = q_ind[k];
      q_ind[k] = t;
    }
  }
}

static void
srqsort (int left, int right)
{
  int i, last, i1, j1, i2, j2, t;

  /* quick sort */

  if (left >= right)		/* Do nothing if array contains */
    return;			/* fewer than two elements */

  if ((right - left) <= q_cutoff)
    /* use insertion sort to avoid tail recursion */
    insertsort (left, right);
  else
  {
    indexswap (left, (left + right) / 2);	/* Move partitiion element */
    last = left;		/* to m[0] */

    i2 = q_r[q_ind[left]];
    j2 = q_c[q_ind[left]];
    for (i = left + 1; i <= right; i++)	/* Partition */
    {
      i1 = q_r[q_ind[i]];
      j1 = q_c[q_ind[i]];
      /* double keys */
      if (i1 < i2)
      {
	++last;
	t = q_ind[i];
	q_ind[i] = q_ind[last];
	q_ind[last] = t;
	/* indexswap (++last, i); */
      }
      else if (i1 == i2)
	if (j1 < j2)
	{
	  ++last;
	  t = q_ind[i];
	  q_ind[i] = q_ind[last];
	  q_ind[last] = t;
	  /* indexswap (++last, i); */
	}
    }
    indexswap (left, last);	/* Restore partition element */
    srqsort (left, last - 1);
    srqsort (last + 1, right);
  }
}

MSR *
mdr_real_spconvert (MDR * m)
{
  int i, j, i_old, j_old, ii, nnz, ptr, row, nterms, cnt, tmp;
  int maxr, maxc, duplicate;
  int *ind, *row_cnt, *rowid, *colid;
  double val;
  MSR *sm;

  /* sorting */
  nterms = MNR (m);
  ind = (int *) GC_MAIOP (nterms * sizeof (int));
  rowid = (int *) GC_MAIOP (nterms * sizeof (int));
  colid = (int *) GC_MAIOP (nterms * sizeof (int));
  for (i = 0; i < nterms; i++)
  {
    ind[i] = i;
    rowid[i] = (int) Mdr0 (m, i, 0);
    colid[i] = (int) Mdr0 (m, i, 1);
  }
  q_ind = ind;
  q_r = rowid;
  q_c = colid;
  srqsort (0, nterms - 1);

  /* Figure out the row dimension */
  maxr = rowid[ind[nterms - 1]];

  /* Now, figure out how many elements in each row.
   * row_cnt has the count of the non-zeros in each
   * row of the matrix.  Figure out column dimension.
   */
  row_cnt = (int *) GC_MAIOP ((maxr + 1) * sizeof (int));	/* later sm->ia */
  for (i = 0; i < maxr; i++)
    row_cnt[i] = 0;

  maxc = duplicate = 0;
  i_old = j_old = -1;
  for (ii = 0; ii < nterms; ii++)
  {
    row = ind[ii];
    i = rowid[row];
    j = colid[row];
    /* zero or duplicate ? */
    if ((Mdr0 (m, row, 2) == 0.0) || (i == i_old && j == j_old))
      duplicate++;
    else
    {
      row_cnt[i - 1]++;
      i_old = i;
      j_old = j;
    }
    if (j > maxc)
      maxc = j;
  }

  /* now we know nnz and dimensions */
  nnz = nterms - duplicate;
  sm = (void *) msr_Create (maxr, maxc);
  sm->ia = row_cnt;
  msr_Setup ((MSR *) sm, nnz);
  cnt = sm->ia[0];
  sm->ia[0] = 1;
  for (i = 1; i < maxr + 1; i++)
  {
    tmp = sm->ia[i];
    sm->ia[i] = sm->ia[i - 1] + cnt;
    cnt = tmp;
  }

  /* copy & sum up */
  ptr = i_old = j_old = -1;
  for (ii = 0; ii < nterms; ii++)
  {
    row = ind[ii];
    i = rowid[row];
    j = colid[row];
    val = Mdr0 (m, row, 2);
    /* zero or duplicate ? */
    if ((val == 0.0) || (i == i_old && j == j_old))
    {
      /* sum up */
      if (val)
	sm->d[ptr] += val;
    }
    else
    {
      /* new entry */
      ptr++;
      sm->ja[ptr] = j;
      sm->d[ptr] = val;
      i_old = i;
      j_old = j;
    }
  }

  /* clean up */
  GC_FREE (colid);
  GC_FREE (rowid);
  GC_FREE (ind);

  return sm;
}

MSC *
mdr_complex_spconvert (MDR * m)
{
  int i, j, i_old, j_old, ii, nnz, ptr, row, nterms, cnt, tmp;
  int maxr, maxc, duplicate;
  int *ind, *row_cnt, *rowid, *colid;
  double valr, vali;
  MSC *sm;

  /* sorting */
  nterms = MNR (m);
  ind = (int *) GC_MAIOP (nterms * sizeof (int));
  rowid = (int *) GC_MAIOP (nterms * sizeof (int));
  colid = (int *) GC_MAIOP (nterms * sizeof (int));
  for (i = 0; i < nterms; i++)
  {
    ind[i] = i;
    rowid[i] = (int) Mdr0 (m, i, 0);
    colid[i] = (int) Mdr0 (m, i, 1);
  }
  q_ind = ind;
  q_r = rowid;
  q_c = colid;
  srqsort (0, nterms - 1);

  /* Figure out the row dimension */
  maxr = rowid[ind[nterms - 1]];

  /* Now, figure out how many elements in each row.
   * row_cnt has the count of the non-zeros in each
   * row of the matrix.  Figure out column dimension.
   */
  row_cnt = (int *) GC_MAIOP ((maxr + 1) * sizeof (int));
  for (i = 0; i < maxr; i++)
    row_cnt[i] = 0;

  maxc = duplicate = 0;
  i_old = j_old = -1;
  for (ii = 0; ii < nterms; ii++)
  {
    row = ind[ii];
    i = rowid[row];
    j = colid[row];
    /* zero or duplicate ? */
    if ((Mdr0 (m, row, 2) == 0.0 && Mdr0 (m, row, 3) == 0.0) ||
	(i == i_old && j == j_old))
      duplicate++;
    else
    {
      row_cnt[i - 1]++;
      i_old = i;
      j_old = j;
    }
    if (j > maxc)
      maxc = j;
  }

  /* now we know nnz and dimensions */
  nnz = nterms - duplicate;
  sm = (void *) msc_Create (maxr, maxc);
  sm->ia = row_cnt;
  msc_Setup ((MSC *) sm, nnz);
  cnt = sm->ia[0];
  sm->ia[0] = 1;
  for (i = 1; i < maxr + 1; i++)
  {
    tmp = sm->ia[i];
    sm->ia[i] = sm->ia[i - 1] + cnt;
    cnt = tmp;
  }

  /* copy & sum up */
  ptr = i_old = j_old = -1;
  for (ii = 0; ii < nterms; ii++)
  {
    row = ind[ii];
    i = rowid[row];
    j = colid[row];
    valr = Mdr0 (m, row, 2);
    vali = Mdr0 (m, row, 3);
    /* zero or duplicate ? */
    if ((valr == 0.0 && vali == 0.0) || (i == i_old && j == j_old))
    {
      /* sum up */
      if (valr)
        sm->c[ptr] += valr;
      if (vali)
        sm->c[ptr] += vali*I;
    }
    else
    {
      /* new entry */
      ptr++;
      sm->ja[ptr] = j;
      sm->c[ptr] = valr + vali*I;
      i_old = i;
      j_old = j;
    }
  }

  /* clean up */
  GC_FREE (colid);
  GC_FREE (rowid);
  GC_FREE (ind);

  return sm;
}

#else

MSR *
mdr_real_spconvert (MDR * m)
{
  int i, ia, ianew, j, k, jj, n, ndup, nnz, rptr, nrow;
  int maxr, maxc, row, col;
  int *row_cnt, *ia_tmp, *ja_tmp;
  double val, *d_tmp;
  MSR *sm, *smnew, *tmp;

  /* Figure out the row and column dimensions. */
  nnz = MNR (m);
  maxr = 0;
  maxc = 0;
  for (i = 0; i < nnz; i++)
  {
    if (Mdr0 (m, i, 1) > maxc)
      maxc = Mdr0 (m, i, 1);
    if (Mdr0 (m, i, 0) > maxr)
      maxr = Mdr0 (m, i, 0);
  }

  sm = (void *) msr_Create (maxr, maxc);
  msr_Setup ((MSR *) sm, nnz);

  /* Create some working space. */
  row_cnt = (int *) GC_MAIOP ((maxr) * sizeof (int));

  /* Zero things out to be sure. */
  for (i = 0; i < maxr; i++)
  {
    sm->ia[i] = 0;
    row_cnt[i] = 0;
  }
  sm->ia[maxr] = 0;

  for (i = 0; i < nnz; i++)
    sm->ja[i] = 0;

  /*
   * Now, figure out how many elements in each row.
   * row_cnt has the count of the non-zeros in each
   * row of the matrix. We will use it to construct
   * IA later.
   */

  for (i = 0; i < MNR (m); i++)
  {
    row_cnt[((int) Mdr0 (m, i, 0)) - 1] += 1;
  }

  /*
   * Now generate IA. IA is simply the running sum
   * of row_count.
   */

  sm->ia[0] = 1;
  for (i = 1; i < maxr + 1; i++)
  {
    sm->ia[i] = sm->ia[i - 1] + row_cnt[i - 1];
  }

  /* Now fill in JA, and AN. */
  for (i = 0; i < nnz; i++)
  {
    row = Mdr0 (m, i, 0);
    col = Mdr0 (m, i, 1);
    val = Mdr0 (m, i, 2);

    /* Get the offset into JA and AN. */
    rptr = sm->ia[row - 1];
    nrow = sm->ia[row] - rptr;

    for (j = 0; j < nrow; j++)
    {
      if (sm->ja[rptr + j - 1] == 0)
      {
	/* Fill it in. */
	sm->ja[rptr + j - 1] = col;
	sm->d[rptr + j - 1] = val;
	break;
      }
    }
  }

  /* Clean Up */
  GC_FREE (row_cnt);

  /*
   * Now, go back, and clear up duplicates in each row.
   * This is a tough job, but it must be done prior
   * to re-ordering.
   *
   * Make an empty matrix structure the same size, and
   * copy into it bit by bit until it is full. Be careful
   * not to copy any duplicates into it. We will have to
   * re-adjust IA as SMNEW is loaded, and then compact it
   * when we are done.
   */

  smnew = msr_Create (sm->nr, sm->nc);
  msr_Setup (smnew, sm->nnz);
  nnz = smnew->nnz;
  smnew->ia[0] = 1;

  /* Row-by-row */
  for (i = 0; i < sm->nr; i++)
  {
    /* N in the row, and IA marker for old structure. */
    n = sm->ia[i + 1] - sm->ia[i];
    ia = sm->ia[i];

    /* IA marker for the new structure. */
    ianew = smnew->ia[i];

    /* Re-set number of duplicates for this pass. */
    ndup = 0;
    jj = 0;

    /* Element by element in each row. */
    for (j = 0; j < n; j++)
    {
      if (sm->ja[ia + j - 1] != 0)
      {
	/* Copy this element into the new struct. */
	smnew->ja[ianew + jj - 1] = sm->ja[ia + j - 1];
	smnew->d[ianew + jj - 1] = sm->d[ia + j - 1];

	/* Now check this element against the remaining. */
	for (k = j + 1; k < n; k++)
	{
	  if (smnew->ja[ianew + jj - 1] == sm->ja[ia + k - 1])
	  {
	    /* Sum the value, and zero out the old struct's
	     * JA so we won't use this one again. */
	    smnew->d[ianew + jj - 1] += sm->d[ia + k - 1];
	    sm->ja[ia + k - 1] = 0;
	    ndup++;
	    nnz--;
	  }
	  else
	  {
	    /* Just keep going. */
	  }
	}
	jj++;
      }
      else
      {
	/* Skip over this element, we already have it (dup). */
	continue;
      }
    }
    /* Set proper value of IA. */
    smnew->ia[i + 1] = smnew->ia[i] + (n - ndup);
  }

  /*
   * Now we must "compact" the reduced matrix (smnew).
   */

  ja_tmp = (int *) GC_MAIOP (nnz * sizeof (int));
  d_tmp = (double *) GC_MAIOP (nnz * sizeof (double));

  memcpy (ja_tmp, smnew->ja, nnz * sizeof (int));
  memcpy (d_tmp, smnew->d, nnz * sizeof (double));

  GC_FREE (smnew->ja);
  GC_FREE (smnew->d);

  smnew->ja = ja_tmp;
  smnew->d = d_tmp;
  smnew->nnz = nnz;

  msr_Destroy (sm);
  sm = smnew;

  /* msr_Check (sm); */

  /*
   * Order the resulting sparse matrix.
   */

  ia_tmp = (int *) GC_MAIOP (((sm->nc) + 1) * sizeof (int));
  ja_tmp = (int *) GC_MAIOP ((sm->nnz) * sizeof (int));
  d_tmp = (double *) GC_MAIOP ((sm->nnz) * sizeof (double));

  DGSTRN (sm->ia, sm->ja, sm->d, &sm->nr, &sm->nc, ia_tmp, ja_tmp, d_tmp);

  DGSTRN (ia_tmp, ja_tmp, d_tmp, &sm->nc, &sm->nr, sm->ia, sm->ja, sm->d);

  GC_FREE (ia_tmp);
  GC_FREE (ja_tmp);
  GC_FREE (d_tmp);

  /* msr_Check (sm); */

  /* Re-Sparse the matrix to get rid of any zeros... */
  tmp = msr_ReSparse (sm);

  /* msr_Check (tmp); */

  msr_Destroy (sm);

  return (tmp);
}


MSC *
mdr_complex_spconvert (MDR * m)
{
  int i, ia, ianew, j, jj, k, n, ndup, nnz, rptr, nrow;
  int maxr, maxc, row, col;
  int *row_cnt, *ia_tmp, *ja_tmp;
  double valr, vali;
  Complex *c_tmp;
  MSC *sm, *smnew, *tmp;

  /* Figure out the row and column dimensions. */
  nnz = MNR (m);
  maxr = 0;
  maxc = 0;
  for (i = 0; i < nnz; i++)
  {
    if (Mdr0 (m, i, 1) > maxc)
      maxc = Mdr0 (m, i, 1);
    if (Mdr0 (m, i, 0) > maxr)
      maxr = Mdr0 (m, i, 0);
  }

  sm = (void *) msc_Create (maxr, maxc);
  msc_Setup (sm, nnz);

  /* Create some working space. */
  row_cnt = (int *) GC_MAIOP ((maxr) * sizeof (int));

  /* Zero things out to be sure. */
  for (i = 0; i < maxr; i++)
  {
    sm->ia[i] = 0;
    row_cnt[i] = 0;
  }
  sm->ia[maxr] = 0;

  for (i = 0; i < nnz; i++)
    sm->ja[i] = 0;

  /*
   * Now, figure out how many elements in each row.
   * row_cnt has the count of the non-zeros in each
   * row of the matrix. We will use it to construct
   * IA later.
   */

  for (i = 0; i < MNR (m); i++)
  {
    row_cnt[((int) Mdr0 (m, i, 0)) - 1] += 1;
  }

  /*
   * Now generate IA. IA is simply the running sum
   * of row_count.
   */

  sm->ia[0] = 1;
  for (i = 1; i < maxr + 1; i++)
  {
    sm->ia[i] = sm->ia[i - 1] + row_cnt[i - 1];
  }

  /* Now fill in JA, and AN. */
  for (i = 0; i < nnz; i++)
  {
    row = Mdr0 (m, i, 0);
    col = Mdr0 (m, i, 1);
    valr = Mdr0 (m, i, 2);
    vali = Mdr0 (m, i, 3);

    /* Get the offset into JA and AN. */
    rptr = sm->ia[row - 1];
    nrow = sm->ia[row] - rptr;

    for (j = 0; j < nrow; j++)
    {
      if (sm->ja[rptr + j - 1] == 0)
      {
        /* Fill it in. */
        sm->ja[rptr + j - 1] = col;
        sm->c[rptr + j - 1] = valr + vali*I;
        break;
      }
    }
  }

  /* Clean Up */
  GC_FREE (row_cnt);

  /*
   * Now, go back, and clear up duplicates in each row.
   * This is a tough job, but it must be done prior
   * to re-ordering.
   *
   * Make an empty matrix structure the same size, and
   * copy into it bit by bit until it is full. Be careful
   * not to copy any duplicates into it. We will have to
   * re-adjust IA as SMNEW is loaded, and then compact it
   * when we are done.
   */

  smnew = msc_Create (sm->nr, sm->nc);
  msc_Setup (smnew, sm->nnz);
  nnz = smnew->nnz;
  smnew->ia[0] = 1;

  /* Row-by-row */
  for (i = 0; i < sm->nr; i++)
  {
    /* N in the row, and IA marker for old structure. */
    n = sm->ia[i + 1] - sm->ia[i];
    ia = sm->ia[i];

    /* IA marker for the new structure. */
    ianew = smnew->ia[i];

    /* Re-set number of duplicates for this pass. */
    ndup = 0;
    jj = 0;

    /* Element by element in each row. */
    for (j = 0; j < n; j++)
    {
      if (sm->ja[ia + j - 1] != 0)
      {
        /* Copy this element into the new struct. */
        smnew->ja[ianew + jj - 1] = sm->ja[ia + j - 1];
        smnew->c[ianew + jj - 1]  = sm->c[ia + j - 1];

        /* Now check this element against the remaining. */
        for (k = j + 1; k < n; k++)
        {
          if (smnew->ja[ianew + jj - 1] == sm->ja[ia + k - 1])
          {
            /* Sum the value, and zero out the old struct's
            * JA so we won't use this one again. */
            smnew->c[ianew + jj - 1] += sm->c[ia + k - 1];
            sm->ja[ia + k - 1] = 0;
            ndup++;
            nnz--;
          }
        }
        jj++;
      }
      else
      {
        /* Skip over this element, we already have it (dup). */
        continue;
      }
    }
    /* Set proper value of IA. */
    smnew->ia[i + 1] = smnew->ia[i] + (n - ndup);
  }

  /*
   * Now we must "compact" the reduced matrix (smnew).
   */

  ja_tmp = (int *) GC_MAIOP (nnz * sizeof (int));
  c_tmp = (Complex *) GC_MAIOP (nnz * sizeof (Complex));

  memcpy (ja_tmp, smnew->ja, nnz * sizeof (int));
  memcpy (c_tmp, smnew->c, nnz * sizeof (Complex));

  GC_FREE (smnew->ja);
  GC_FREE (smnew->c);

  smnew->ja = ja_tmp;
  smnew->c = c_tmp;
  smnew->nnz = nnz;

  msc_Destroy (sm);
  sm = smnew;

  /* Order the resulting sparse matrix. */
  ia_tmp = (int *) GC_MAIOP (((sm->nc) + 1) * sizeof (int));
  ja_tmp = (int *) GC_MAIOP ((sm->nnz) * sizeof (int));
  c_tmp = (Complex *) GC_MAIOP ((sm->nnz) * sizeof (Complex));

  ZGSTRN (sm->ia, sm->ja, sm->c, &sm->nr, &sm->nc, ia_tmp, ja_tmp, c_tmp);

  ZGSTRN (ia_tmp, ja_tmp, c_tmp, &sm->nc, &sm->nr, sm->ia, sm->ja, sm->c);

  GC_FREE (ia_tmp);
  GC_FREE (ja_tmp);
  GC_FREE (c_tmp);

  /* Re-Sparse the matrix, to get rid of any zeros... */
  tmp = msc_ReSparse (sm);
  msc_Destroy (sm);

  return (tmp);
}
#endif /* tsy */

/* **************************************************************
 * Convert a sparse matrix to dense-3-column representation.
 * 1st column:  Row id
 * 2nd column:  Column id
 * 3rd column:  matrix elements (non-zeros)
 * ************************************************************** */

MDR *
msr_Spconvert (MSR * m, int *rtype)
{
  MDR *dm;
  int i, ii, ia, j, n, row;

  *rtype = MATRIX_DENSE_REAL;
  if (((m->nr == 0) && (m->nc == 0)) || (m->nnz == 0))
  {
    dm = mdr_Create (0, 0);
    return (dm);
  }

  dm = mdr_Create (m->nnz, 3);

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
      Mdr1 (dm, ii, 3) = m->d[ia + j - 1];
      ii++;
    }
    row++;
  }
  return (dm);
}

/* **************************************************************
 * Append: [ Sparse , Dense ]
 * ************************************************************** */

MDR *
msr_mdr_Append (MSR * m1, MDR * m2)
{
  MDR *m, *mtmp;

  mtmp = msr_Dense (m1);
  m = mdr_Append (mtmp, m2);
  mdr_Destroy (mtmp);

  return (m);
}

/* **************************************************************
 * Append: [ Dense , Sparse ]
 * ************************************************************** */

MDR *
mdr_msr_Append (MDR * m1, MSR * m2)
{
  MDR *m, *mtmp;

  mtmp = msr_Dense (m2);
  m = mdr_Append (m1, mtmp);
  mdr_Destroy (mtmp);

  return (m);
}

/* **************************************************************
 * Stack: [ Sparse ; Dense ]
 * ************************************************************** */

MDR *
msr_mdr_Stack (MSR * m1, MDR * m2)
{
  MDR *m, *mtmp;

  mtmp = msr_Dense (m1);
  m = mdr_Stack (mtmp, m2);
  mdr_Destroy (mtmp);

  return (m);
}

/* **************************************************************
 * Stack: [ Dense ; Sparse ]
 * ************************************************************** */

MDR *
mdr_msr_Stack (MDR * m1, MSR * m2)
{
  MDR *m, *mtmp;

  mtmp = msr_Dense (m2);
  m = mdr_Stack (m1, mtmp);
  mdr_Destroy (mtmp);

  return (m);
}

/* **************************************************************
 * Add two sparse matrices.
 * ************************************************************** */

MSR *
msr_Add (MSR * m1, MSR * m2)
{
  MSR *new;
  double *w;
  int *cia, *cja, *iw;

  /* Some special cases. */
  if (((m1->nr == 0) && (m1->nc == 0)) || ((m2->nr == 0) && (m2->nc == 0)))
  {
    new = msr_Create (0, 0);
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
    new = msr_Create (m1->nr, m1->nc);
    return (new);
  }

  if (m1->nnz == 0)
  {
    return (msr_Copy (m2));
  }

  if (m2->nnz == 0)
  {
    return (msr_Copy (m1));
  }

  /* Create the new matrix. */
  new = msr_Create (m1->nr, m1->nc);

  /* Create the working space for the symbolic addition. */
  cia = (int *) GC_malloc_atomic_ignore_off_page ((m1->nr + 1) * sizeof (int));
  cja =
    (int *) GC_malloc_atomic_ignore_off_page ((m1->nnz + m2->nnz) *
					      sizeof (int));
  iw = (int *) GC_malloc_atomic_ignore_off_page ((m1->nc + 1) * sizeof (int));

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

  new->ja =
    (int *) GC_malloc_atomic_ignore_off_page ((new->nnz) * sizeof (int));
  memcpy (new->ja, cja, new->nnz * sizeof (int));

  /* Transpose twice to order the new structure. */
  XGSTRN (new->ia, new->ja, &new->nr, &new->nc, iw, cja);
  XGSTRN (iw, cja, &new->nc, &new->nr, new->ia, new->ja);
  GC_FREE (cja);
  GC_FREE (iw);

  /* Make room for the elements. */
  new->d =
    (double *) GC_malloc_atomic_ignore_off_page (new->nnz * sizeof (double));

  /* Now do the add. */
  w = (double *) GC_malloc_atomic_ignore_off_page (new->nc * sizeof (double));
  DGSADD (m1->ia, m1->ja, m1->d,
	  m2->ia, m2->ja, m2->d, &m1->nr, &m1->nc, new->ia, new->ja, new->d, w);

  GC_FREE (w);
  return (new);
}

/* **************************************************************
 * Add a dense real matrix, and a sparse real matrix.
 * ************************************************************** */

MDR *
mdr_msr_Add (MDR * m1, MSR * m2)
{
  MDR *m, *mtmp;

  mtmp = msr_Dense (m2);
  m = mdr_Add (m1, mtmp);
  mdr_Destroy (mtmp);

  return (m);
}

/* **************************************************************
 * Add a sparse real matrix, and a dense real matrix.
 * ************************************************************** */

MDR *
msr_mdr_Add (MSR * m1, MDR * m2)
{
  MDR *m, *mtmp;

  mtmp = msr_Dense (m1);
  m = mdr_Add (mtmp, m2);
  mdr_Destroy (mtmp);

  return (m);
}

/* **************************************************************
 * Subtract two sparse matrices.
 * ************************************************************** */

MSR *
msr_Subtract (MSR * m1, MSR * m2)
{
  MSR *new;
  double *w;
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
    new = msr_Create (m1->nr, m1->nc);
    return (new);
  }

  if (m1->nnz == 0)
  {
    return (msr_Negate (m2));
  }

  if (m2->nnz == 0)
  {
    return (msr_Copy (m1));
  }

  /* Create the new matrix. */
  new = msr_Create (m1->nr, m1->nc);

  /* Create the working space for the symbolic addition. */
  cia = (int *) GC_malloc_atomic_ignore_off_page ((m1->nr + 1) * sizeof (int));
  cja =
    (int *) GC_malloc_atomic_ignore_off_page ((m1->nnz + m2->nnz) *
					      sizeof (int));
  iw = (int *) GC_malloc_atomic_ignore_off_page ((m1->nc + 1) * sizeof (int));

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

  new->ja =
    (int *) GC_malloc_atomic_ignore_off_page ((new->nnz) * sizeof (int));
  memcpy (new->ja, cja, new->nnz * sizeof (int));

  /* Transpose twice to order the new structure. */
  XGSTRN (new->ia, new->ja, &new->nr, &new->nc, iw, cja);
  XGSTRN (iw, cja, &new->nc, &new->nr, new->ia, new->ja);
  GC_FREE (cja);
  GC_FREE (iw);

  /* Make room for the elements. */
  new->d =
    (double *) GC_malloc_atomic_ignore_off_page (new->nnz * sizeof (double));

  /* Now do the subtract. */
  w = (double *) GC_malloc_atomic_ignore_off_page (new->nc * sizeof (double));
  DGSSUB (m1->ia, m1->ja, m1->d,
	  m2->ia, m2->ja, m2->d, &m1->nr, &m1->nc, new->ia, new->ja, new->d, w);

  GC_FREE (w);
  return (new);
}

/* **************************************************************
 * Subtract a dense real matrix, and a sparse real matrix.
 * ************************************************************** */

MDR *
mdr_msr_Subtract (MDR * m1, MSR * m2)
{
  MDR *m, *mtmp;

  mtmp = msr_Dense (m2);
  m = mdr_Subtract (m1, mtmp);
  mdr_Destroy (mtmp);

  return (m);
}

/* **************************************************************
 * Subtract a sparse real matrix, and a dense real matrix.
 * ************************************************************** */

MDR *
msr_mdr_Subtract (MSR * m1, MDR * m2)
{
  MDR *m, *mtmp;

  mtmp = msr_Dense (m1);
  m = mdr_Subtract (mtmp, m2);
  mdr_Destroy (mtmp);

  return (m);
}

/* **************************************************************
 * Support functions for sparse matrix multiply.
 * ************************************************************** */

MSR *
msr_scalar_mult (MSR * m, double s)
{
  int i;
  MSR *new;

  new = msr_Create (m->nr, m->nc);
  if (m->nnz)
  {
    msr_Setup (new, m->nnz);

    memcpy (new->ia, m->ia, (m->nr + 1) * sizeof (int));
    memcpy (new->ja, m->ja, (m->nnz) * sizeof (int));

    for (i = 0; i < m->nnz; i++)
    {
      new->d[i] = s * m->d[i];
    }
  }
  return (new);
}

/* **************************************************************
 * Multiply two sparse matrices.
 * ************************************************************** */

MSR *
msr_Multiply (MSR * m1, MSR * m2, int *rtype)
{
  MSR *new = 0;
  int maxjc;
  int *ja_tmp, *tmp;
  double *a_tmp;

  *rtype = MATRIX_SPARSE_REAL;

  /* Check dimensions... */

  if (m1->nc != m2->nr)
  {
    /*
     * Handle condition where one of the operands is a
     * scalar, or an empty matrix.
     */

    if ((m1->nr == 1) && (m1->nc == 1))
    {
      /* A is scalar */
      new = msr_scalar_mult (m2, m1->d[0]);
    }
    else if ((m2->nr == 1) && (m2->nc == 1))
    {
      new = msr_scalar_mult (m1, m2->d[0]);
    }
    else if ((m1->nr == 0) || (m1->nc == 0))
    {
      new = msr_Create (0, 0);
    }
    else if ((m2->nr == 0) || (m2->nc == 0))
    {
      new = msr_Create (0, 0);
    }
    else
    {
      fprintf (stderr, "matrix dimensions must be consistent\n");
      fprintf (stderr, "matrix col and row sizes: %i and %i\n", m1->nc, m2->nr);
      rerror ("sparse-multiply: matrices must have the same row dimension");
    }
  }
  else
  {
    /*
     * Now get on with the multiply...
     */

    new = msr_Create (m1->nr, m2->nc);
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
      new->d = (double *) GC_MAIOP (new->nnz * sizeof (double));
      a_tmp = (double *) GC_MAIOP ((new->nc + 1) * sizeof (double));
      DGSMUL (m1->ia, m1->ja, m1->d,
	      m2->ia, m2->ja, m2->d,
	      &new->nr, &new->nc, new->ia, new->ja, new->d, a_tmp);
      GC_FREE (a_tmp);
    }
  }
  return (new);
}

void *
mdr_msr_Multiply (MDR * l, MSR * r, int *rtype)
{
  int i, j, k;
  void *m = 0;

  /* Check sizes ... */
  if (MNC (l) != r->nr)
  {
    if ((MNR (l) == 1) && (MNC (l) == 1))
    {
      /* The dense one is scalar, return sparse. */
      m = (void *) msr_scalar_mult (r, Mdr0 (l, 0, 0));
      *rtype = MATRIX_SPARSE_REAL;
    }
    else if ((r->nr == 1) && (r->nc == 1))
    {
      /* The sparse one is scalar, return dense. */
      m = mdr_Create (MNR (l), MNC (l));
      for (i = 0; i < MNR (l) * MNC (l); i++)
      {
	MdrV0 (m, 0) = MdrV0 (l, 0) * r->d[0];
      }
      *rtype = MATRIX_DENSE_REAL;
    }
    else if ((MNR (l) == 0) || (MNC (l) == 0))
    {
      m = msr_Create (0, 0);
      *rtype = MATRIX_SPARSE_REAL;
    }
    else if ((r->nr == 0) || (r->nc == 0))
    {
      m = msr_Create (0, 0);
      *rtype = MATRIX_SPARSE_REAL;
    }
    else
    {
      fprintf (stderr, "matrix dimensions must be consistent\n");
      fprintf (stderr, "matrix col and row sizes: %i and %i\n", MNC (l), r->nr);
      rerror ("sparse-multiply: matrices must have the same row dimension");
    }
  }
  else
  {
    /*
     * Now get on with the multiply...
     */

    *rtype = MATRIX_DENSE_REAL;
    m = mdr_Create (MNR (l), r->nc);
    mdr_Zero (m);
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
//          Mdr0(m,i,(r->ja[k - 1] - 1)) += Mdr0(l,i,j) * r->d[k - 1];
          Mdr0((MDR *)m,i,(r->ja[k - 1] - 1)) = 0;
        }
      }
    }
  }
  return (m);
}

void *
msr_mdr_Multiply (MSR * l, MDR * r, int *rtype)
{
  int i, j, k1;
  void *m = 0;
  *rtype = MATRIX_DENSE_REAL;

  /* Check sizes ... */
  if (l->nc != MNR (r))
  {
    if ((l->nr == 1) && (l->nc == 1))
    {
      /* The sparse one is scalar, return dense. */
      m = mdr_Create (MNR (r), MNC (r));
      for (i = 0; i < MNR (r) * MNC (r); i++)
      {
        MdrV0 (m, 0) = MdrV0 (r, 0) * l->d[0];
      }
      *rtype = MATRIX_DENSE_REAL;
    }
    else if ((MNR (r) == 1) && (MNC (r) == 1))
    {
      /* The dense one is scalar, return sparse. */
      m = msr_scalar_mult (l, Mdr0 (r, 0, 0));
      *rtype = MATRIX_SPARSE_REAL;
    }
    else if ((l->nr == 0) || (l->nc == 0))
    {
      m = msr_Create (0, 0);
      *rtype = MATRIX_SPARSE_REAL;
    }
    else if ((MNR (r) == 0) || (MNC (r) == 0))
    {
      m = msr_Create (0, 0);
      *rtype = MATRIX_SPARSE_REAL;
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
    /*
     * Now get on with the multiply...
     */

    *rtype = MATRIX_DENSE_REAL;
    m = mdr_Create (l->nr, MNC (r));
    mdr_Zero (m);

    if (l->nnz == 0)
    {
      return (m);
    }

    for (i = 0; i < MNR (m); i++)
    {
      for (j = 0; j < MNC (m); j++)
      {
        for (k1 = l->ia[i] - 1; k1 < l->ia[i + 1] - 1; k1++)
        {
          Mdr0( (MDR *)m,i,j) += l->d[k1] * Mdr0(r, l->ja[k1]-1, j);
        }
      }
    }
  }
  return (m);
}

void *
msr_mdc_Multiply (MSR * m1, MDC * m2, int *rtype)
{
  MSC *msc;
  void *new;

  msc = msr_coerce_msc (m1);
  new = msc_mdc_Multiply (msc, m2, rtype);
  msc_Destroy (msc);
  return (new);
}

/* **************************************************************
 * Matrix element by element multiply.
 * ************************************************************** */

void *
msr_ElMultiply (MSR * ma, MSR * mb, int *rtype)
{
  MSR *mc = 0;
  return (mc);
}

void *
mdr_msr_ElMultiply (MDR * ma, MSR * mb, int *rtype)
{
  MSR *mc = 0;
  return (mc);
}

void *
msr_mdr_ElMultiply (MSR * ma, MDR * mb, int *rtype)
{
  MSR *mc = 0;
  return (mc);
}

/* **************************************************************
 * Sparse Matrix Right Divide
 * ************************************************************** */

MDR *
msr_Rdivide (MSR * m1, MSR * m2)
{
  MDR *l, *r;
  MDR *new;

  /* First convert both to dense. */
  l = msr_Dense (m1);
  r = msr_Dense (m2);

  new = mdr_Rdivide (l, r);
  mdr_Destroy (l);
  mdr_Destroy (r);

  return (new);
}

MDR *
mdr_msr_Rdivide (MDR * m1, MSR * m2)
{
  MDR *r;
  MDR *new;

  /* First convert to dense. */
  r = msr_Dense (m2);

  new = mdr_Rdivide (m1, r);
  mdr_Destroy (r);

  return (new);
}

MDR *
msr_mdr_Rdivide (MSR * m1, MDR * m2)
{
  MDR *l;
  MDR *new;

  /* First convert to dense. */
  l = msr_Dense (m1);

  new = mdr_Rdivide (l, m2);
  mdr_Destroy (l);

  return (new);
}

/* **************************************************************
 * Sparse Matrix Element-by-Element Right Divide
 * ************************************************************** */

MDR *
msr_ElRdivide (MSR * m1, MSR * m2, int *rtype)
{
  MDR *l, *r;
  MDR *new;

  /* First convert both to dense. */
  l = msr_Dense (m1);
  r = msr_Dense (m2);

  new = mdr_ElRdivide (l, r, rtype);
  mdr_Destroy (l);
  mdr_Destroy (r);

  return (new);
}

MDR *
mdr_msr_ElRdivide (MDR * m1, MSR * m2, int *rtype)
{
  MDR *r;
  MDR *new;

  /* First convert to dense. */
  r = msr_Dense (m2);

  new = mdr_ElRdivide (m1, r, rtype);
  mdr_Destroy (r);

  return (new);
}

void *
msr_mdr_ElRdivide (MSR * m1, MDR * m2, int *rtype)
{
  int i;
  void *new;
  MDR *l;

  /*
   * Look for some common operations so that we might
   * retain sparsity.
   */

  if (MNR (m2) == 1 && MNC (m2) == 1)
  {
    *rtype = MATRIX_SPARSE_REAL;
    new = msr_Create (m1->nr, m1->nc);
    ((MSR *) new)->nnz = m1->nnz;

    /* Copy the IA array. */
    if (m1->ia)
    {
      ((MSR *) new)->ia = (int *) GC_malloc_atomic_ignore_off_page
	(sizeof (int) * (m1->nr + 1));
      memcpy (((MSR *) new)->ia, m1->ia, (m1->nr + 1) * sizeof (int));
    }

    /* Copy the JA array. */
    if (m1->ja)
    {
      ((MSR *) new)->ja = (int *) GC_malloc_atomic_ignore_off_page
	(sizeof (int) * (m1->nnz));
      memcpy (((MSR *) new)->ja, m1->ja, m1->nnz * sizeof (int));
    }

    if (m1->d)
    {
      ((MSR *) new)->d = (double *) GC_malloc_atomic_ignore_off_page
	(sizeof (double) * (m1->nnz));
    }

    for (i = 0; i < m1->nnz; i++)
      ((MSR *) new)->d[i] = m1->d[i] / MdrV0 (m2, 0);
  }
  else
  {
    *rtype = MATRIX_DENSE_REAL;
    /* First convert to dense. */
    l = msr_Dense (m1);

    new = mdr_ElRdivide (l, m2, rtype);
    mdr_Destroy (l);
  }

  return (new);
}

/* **************************************************************
 * Sparse Matrix Left Divide ( A \ B )
 * ************************************************************** */

MDR *msr_Solve (MSR * a, MDR * b, char *type);

MDR *
msr_Ldivide (MSR * a, MDR * b)
{
  double dtmp;
  int i, size;
  MDR *new = 0;

  /*
   * Check for special case where denominator (m1) is a scalar.
   * The only thing we gain if m2 is a scalar, is speed.
   */

  if (a->nr == 1 && a->nc == 1)
  {
    new = mdr_Create (MNR (b), MNC (b));
    size = MNR (b) * MNC (b);
    dtmp = a->d[0];
    for (i = 0; i < size; i++)
    {
      MdrV0 (new, i) = MdrV0 (b, i) / dtmp;
    }
    return (new);
  }
  else if (a->nr == 0 || MNR (b) == 0)
  {
    new = mdr_Create (0, 0);
    return (new);
  }

  /* Check dimensions */
  if (a->nr != MNR (b))
    rerror ("RHS row dim. must match LHS row dim.");

  if (a->nr == a->nc)
    new = msr_Solve (a, b, (char *) 0);
  else
    rerror ("cannot left divide sparse non-square systems");

  return (new);
}

MDR *
msr_msr_Ldivide (MSR * a, MSR * b)
{
  MDR *mtmp, *new;

  mtmp = msr_Dense (b);
  new = msr_Ldivide (a, mtmp);
  mdr_Destroy (mtmp);

  return (new);
}

/* **************************************************************
 * Sparse Matrix sizeof()
 * ************************************************************** */

size_t
msr_Sizeof (MSR * m)
{
  int size = m->nnz * sizeof (int) +
    m->nnz * sizeof (double) + (m->nr + 1) * sizeof (int);
  return ((size_t) (size));
}

/* **************************************************************
 * Get the size of a a matrix for size()...
 * ************************************************************** */

MDR *
msr_Size_BF (MSR * m)
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
msr_Type_BF (MSR * m)
{
  MDS *type = mds_CreateScalar ( "real" );
  return (type);
}

/* **************************************************************
 * Builtin functions for any()
 * Return TRUE if ANY of the elements of A are non-zero.
 * ************************************************************** */

MDR *
msr_Any (MSR * m)
{
  int i, j, rptr, nrow;
  MDR *new;
  MSR *mtmp = 0;

  if (m->nr == 1)
  {
    /* Vector operation */
    new = mdr_Create (1, 1);
    Mdr0 (new, 0, 0) = 0.0;
    for (i = 0; i < m->nnz; i++)
    {
      if (m->d[i] != 0.0)
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
    mtmp = msr_Transpose (m);

    for (i = 0; i < mtmp->nr; i++)
    {
      rptr = mtmp->ia[i];
      nrow = mtmp->ia[i + 1] - rptr;
      Mdr0 (new, 0, i) = 0.0;

      for (j = 0; j < nrow; j++)
      {
	if (mtmp->d[rptr + j - 1] != 0.0)
	{
	  Mdr0 (new, 0, i) = 1.0;
	  break;
	}
      }
    }
  }
  if (mtmp)
    msr_Destroy (mtmp);
  return (new);
}

/* **************************************************************
 * Builtin functions for all()
 * Return TRUE if ALL of the elements of A are non-zero.
 * ************************************************************** */

MDR *
msr_All (MSR * m)
{
  int i, j;
  MDR *new;
  MSR *mtmp = 0;

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
      if (msr_GetEl (m, 1, i) == 0.0)
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
    mtmp = msr_Transpose (m);

    for (i = 1; i <= mtmp->nr; i++)
    {
      Mdr0 (new, 0, i - 1) = 1.0;

      for (j = 1; j <= mtmp->nc; j++)
      {
	if (msr_GetEl (mtmp, i, j) == 0.0)
	{
	  Mdr0 (new, 0, i - 1) = 0.0;
	  break;
	}
      }
    }
  }
  if (mtmp)
    msr_Destroy (mtmp);
  return (new);
}

/* **************************************************************
 * Return the complex-conjugate of a a sparse matrix.
 * ************************************************************** */

MSR *
msr_Conj (MSR * m)
{
  return (m);
}

/* **************************************************************
 * Return the "length" of a sparse real matrix
 * ************************************************************** */

MDR *
msr_Length_BF (MSR * m)
{
  MDR *size = mdr_Create (1, 1);
  Mdr0 (size, 0, 0) = MAX(m->nr, m->nc);
  return (size);
}

/* **************************************************************
 * Abs function...
 * ************************************************************** */

MSR *
msr_Abs (MSR * m)
{
  int i;
  MSR *new;

  new = msr_Create (m->nr, m->nc);
  msr_Setup (new, m->nnz);

  memcpy (new->ia, m->ia, (m->nr + 1) * sizeof (int));
  memcpy (new->ja, m->ja, (m->nnz) * sizeof (int));

  for (i = 0; i < m->nnz; i++)
  {
    new->d[i] = ABS(m->d[i]);
  }
  return (new);
}

/* **************************************************************
 * Detect Infs and NaNs in sparse matrices.
 * ************************************************************** */

void
msr_Detect_Inf (MSR * m)
{
  if (detect_inf_r (m->d, m->nnz))
    rerror ("matrix contains Inf value");
}

void
msr_Detect_NaN (MSR * m)
{
  if (detect_nan_r (m->d, m->nnz))
    rerror ("matrix contains NaN value");
}

/* **************************************************************
 * Diag() function for sparse real matrices. Always returns a
 * sparse matrix.
 * ************************************************************** */

MSR *
msr_Diag (MSR * marg, int k)
{
  int i, kk, smin, size;
  MSR *m, *mtmp;

  if (marg->nr == 1 || marg->nc == 1)
  {
    /* Create a diagonal matrix */
    size = MAX(marg->nr, marg->nc) + ABS(k);
    m = msr_Create (size, size);
    msr_Setup (m, marg->nnz);

    /* Force marg into a column vector to make life easier. */
    if (marg->nr == 1)
      mtmp = msr_Transpose (marg);
    else
      mtmp = marg;

    if (k < 0)
    {
      /* Take care of IA first. */
      for (i = 0; i <= -k; i++)
        m->ia[i] = 1;
      kk = 2;
      for (i = -k + 1; i < m->nr + 1; i++)
        m->ia[i] = kk++;

      for (i = 0; i < m->nnz; i++)
        m->d[i] = mtmp->d[i];

      /* Now do JA. */
      kk = 0;
      for (i = 0; i < mtmp->nr; i++)
      {
        if (mtmp->ia[i + 1] != mtmp->ia[i])
        {
          m->ja[kk] = i + 1;
          kk++;
        }
      }
    }
    else      /* k > 0 */
    {
      /* Take care of IA first. */
      for (i = 0; i < mtmp->nr + 1; i++)
        m->ia[i] = mtmp->ia[i];
      for (i = mtmp->nr + 1; i < m->nr + 1; i++)
        m->ia[i] = mtmp->ia[mtmp->nr];

      for (i = 0; i < m->nnz; i++)
        m->d[i] = mtmp->d[i];

      /* Now do JA. */
      kk = 0;
      for (i = 0; i < mtmp->nr; i++)
      {
        if (mtmp->ia[i + 1] != mtmp->ia[i])
        {
          m->ja[kk] = i + 1 + k;
          kk++;
        }
      }
    }
    if (mtmp != marg)
      msr_Destroy (mtmp);
  }
  else
  {
    /* Extract the diagonal elements */
    smin = MIN(marg->nr, marg->nc);
    //     smax = MAX(marg->nr, marg->nc);

    if (k >= 0)
    {
      if (marg->nr >= marg->nc)
        size = smin - k;
      else
        size = smin - MAX(0, k - (marg->nc - marg->nr));
    }
    else
    {
      if (marg->nr >= marg->nc)
        size = smin - MAX(0, -k - (marg->nr - marg->nc));
      else
        size = smin + k;
    }
    if (size <= 0)
      size = 0;

    m = msr_Create (size, 1);
    msr_Setup (m, size);

    if (k >= 0)
    {
      m->ia[0] = 1;
      for (i = 1; i <= size; i++)
      {
        m->d[i - 1] = msr_GetEl (marg, i, i + k);
        m->ia[i] = i + 1;
        m->ja[i - 1] = 1;
      }
    }
    else
    {
      m->ia[0] = 1;
      for (i = 1; i <= size; i++)
      {
        m->d[i - 1] = msr_GetEl (marg, i - k, i);
        m->ia[i] = i + 1;
        m->ja[i - 1] = 1;
      }
    }
  }
  return (m);
}

/* **************************************************************
 * Test a sparse real matrix for symmetry.
 * ************************************************************** */

int
msr_IsSymmetric (MSR * m)
{
  int i;
  MSR *mt;

  if (m->nr != m->nc)
    return (0);

  /* Check that M and Mt are the same. */
  mt = msr_Transpose (m);
  for (i = 0; i < m->nnz; i++)
  {
    if ((m->d[i] != mt->d[i]) || (m->ja[i] != mt->ja[i]))
    {
      msr_Destroy (mt);
      return (0);
    }
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
msr_Find_BF (MSR * m)
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
 * Return the real part of a sparse-real-matrix.
 * ************************************************************** */

MSR *
msr_Real_BF (MSR * m)
{
  return (m);
}

/* **************************************************************
 * Return the imaginary part of a sparse-real-matrix. This is
 * always a zero matrix...
 * ************************************************************** */

MSR *
msr_Imag_BF (MSR * m)
{
  MSR *i = msr_Create (m->nr, m->nc);
  msr_Setup (i, 0);
  return (i);
}

MDR *
msr_Sum_BF (MSR * m, void *h)
{
  int i, j;
  MDR *rm;

  if (m->nr == 1)
  {
    /* Single row, return single double value. */
    double d = 0.0;
    for (i = 0; i < m->nnz; i++)
    {
      d = d + m->d[i];
    }
    rm = mdr_CreateScalar (d);
  }
  else
  {
    /*
     * N-by-M matrix, return a vector.
     * A sum-value for each column in the matrix.
     * Transpose the matrix, in-efficient but easy.
     */

    int ia, n;
    MSR *stmp;

    stmp = msr_Transpose (m);
    rm = mdr_Create (1, m->nr);
    mdr_Zero (rm);

    /* Now, compute the sum for each row (was column). */
    for (i = 0; i < stmp->nr; i++)
    {
      n = stmp->ia[i + 1] - stmp->ia[i];
      ia = stmp->ia[i];

      for (j = 0; j < n; j++)
      {
	MdrV0 (rm, i) = MdrV0 (rm, i) + stmp->d[ia + j - 1];
      }
    }
    msr_Destroy (stmp);
  }
  return (rm);
}

double *
msr_Norm (MSR * m, char *type)
{
  double *norm = (double *) GC_MALLOC (sizeof (double));
  *norm = 0.0;

  msr_Detect_Inf (m);
  msr_Detect_NaN (m);

  /* Argument error checking */
  if (!strcmp (type, "1") || !strcmp (type, "O") || !strcmp (type, "o"))
  {
    /* largest column sum. */
    int i, size;
    MSR *mtmp = msr_Abs (m);
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
    MSR *mtmp1 = msr_Abs (m);
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

MDR *
msr_Max1 (MSR * m)
{
  int i, j;
  double maxr;
  MDR *mmax = 0;

  if (m->nr == 1)
  {
    mmax = mdr_CreateScalar (0.0);
    /* Input is a row-vector, return scalar max. */
    if (m->nnz > 0)
    {
      maxr = m->d[0];
      for (i = 1; i < m->nnz; i++)
      {
        if (m->d[i] > maxr)
          maxr = m->d[i];
      }
      Mdr0 (mmax, 0, 0) = maxr;
    }
  }
  else
  {
    if (m->nnz == 0)
    {
      mmax = mdr_Create(1, m->nc);
      mdr_Zero (mmax);
    }
    else
    {
     /*
      * General matrix case.
      * Go through the columns, one-by-one.
      * Transpose the matrix, to make life easier.
      */
      int ia, n;
      MSR *mtmp = msr_Transpose (m);
      mmax = mdr_Create (1, mtmp->nr);

      /* Row-by-row */
      for (i = 0; i < mtmp->nr; i++)
      {
        /* N(umber) in the row. */
        n = mtmp->ia[i + 1] - mtmp->ia[i];
        ia = mtmp->ia[i];
        maxr = mtmp->d[ia - 1];

        for (j = 0; j < n; j++)
        {
          if (mtmp->d[ia + j - 1] > maxr)
            maxr = mtmp->d[ia + j - 1];
        }
        Mdr0 (mmax, 0, i) = maxr;
      }
      msr_Destroy (mtmp);
    }
  }
  return (mmax);
}

/* **************************************************************
 * Min function: Minimum scalar value if the input is a vector.
 * Minimum row-vector output if input is a matrix.
 * ************************************************************** */

MDR *
msr_Min1 (MSR * m)
{
  int i, j;
  double maxr;
  MDR *mmax = 0;

  if (m->nr == 1)
  {
    mmax = mdr_CreateScalar (0.0);
    /* Input is a row-vector, return scalar min. */
    if (m->nnz > 0)
    {
      maxr = m->d[0];
      for (i = 1; i < m->nnz; i++)
      {
        if (m->d[i] < maxr)
          maxr = m->d[i];
      }
      mmax = mdr_Create (1, 1);
      Mdr0 (mmax, 0, 0) = maxr;
    }
  }
  else
  {
    if (m->nnz == 0)
    {
      mmax = mdr_Create(1, m->nc);
      mdr_Zero (mmax);
    }
    else
    {
     /*
      * General matrix case.
      * Go through the columns, one-by-one.
      * Transpose the matrix, to make life easier.
      */
      int ia, n;
      MSR *mtmp = msr_Transpose (m);
      mmax = mdr_Create (1, mtmp->nr);

      /* Row-by-row */
      for (i = 0; i < mtmp->nr; i++)
      {
        /* N(umber) in the row. */
        n = mtmp->ia[i + 1] - mtmp->ia[i];
        ia = mtmp->ia[i];
        maxr = mtmp->d[ia - 1];

        for (j = 0; j < n; j++)
        {
          if (mtmp->d[ia + j - 1] < maxr)
            maxr = mtmp->d[ia + j - 1];
        }
        Mdr0 (mmax, 0, i) = maxr;
      }
      msr_Destroy (mtmp);
    }
  }
  return (mmax);
}

/* **************************************************************
 * Int function...
 * ************************************************************** */

MSR *
msr_Int_BF (MSR * m)
{
  int i;
  MSR *im;

  im = msr_Create (m->nr, m->nc);
  msr_Setup (im, m->nnz);

  memcpy (im->ia, m->ia, (m->nr + 1) * sizeof (int));
  memcpy (im->ja, m->ja, (m->nnz) * sizeof (int));

  for (i = 0; i < im->nnz; i++)
  {
    im->d[i] = (double) ((int) m->d[i]);
  }

  return (im);
}

/* **************************************************************
 * Ceil function...
 * ************************************************************** */

MSR *
msr_Ceil_BF (MSR * m)
{
  int i;
  MSR *cm;

  cm = msr_Create (m->nr, m->nc);
  msr_Setup (cm, m->nnz);

  memcpy (cm->ia, m->ia, (m->nr + 1) * sizeof (int));
  memcpy (cm->ja, m->ja, (m->nnz) * sizeof (int));

  for (i = 0; i < cm->nnz; i++)
  {
    cm->d[i] = errcheck (ceil (m->d[i]), "ceil");
  }
  return (cm);
}

/* **************************************************************
 * Floor function...
 * ************************************************************** */

MSR *
msr_Floor_BF (MSR * m)
{
  int i;
  MSR *cm;

  cm = msr_Create (m->nr, m->nc);
  msr_Setup (cm, m->nnz);

  memcpy (cm->ia, m->ia, (m->nr + 1) * sizeof (int));
  memcpy (cm->ja, m->ja, (m->nnz) * sizeof (int));

  for (i = 0; i < cm->nnz; i++)
  {
    cm->d[i] = errcheck (floor (m->d[i]), "floor");
  }
  return (cm);
}

/* **************************************************************
 * Round function...
 * ************************************************************** */

MSR *
msr_Round_BF (MSR * m)
{
  int i;
  MSR *cm;

  cm = msr_Create (m->nr, m->nc);
  msr_Setup (cm, m->nnz);

  memcpy (cm->ia, m->ia, (m->nr + 1) * sizeof (int));
  memcpy (cm->ja, m->ja, (m->nnz) * sizeof (int));

  for (i = 0; i < cm->nnz; i++)
  {
    cm->d[i] = errcheck (rint (m->d[i]), "round");
  }
  return (cm);
}
