/* sort.c - sort related functions. */

/*  This file is a part of RLaB ("Our"-LaB)
   Copyright (C) 1995  Ian R. Searle
   Copyright (C) 2014 M.Kostrun

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

#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include "math_macros.h"
#include "rlab_solver_parameters_names.h"
#include "rlab_macros.h"

// static int _rlab_settings_sort_put_nans = RLAB_SORT_NAN_INPLACE;
static int _rlab_settings_sort_put_nans = RLAB_SORT_NAN_ONBOTTOM;
// static int _rlab_settings_sort_put_nans = RLAB_SORT_NAN_ONTOP;

static int _rlab_settings_sort_method = RLAB_SORT_QUICK;
// int _rlab_settings_sort_method = RLAB_SORT_HEAP;

extern int isnand( double );

int rlab_sort_get_method(void)
{
  return _rlab_settings_sort_method;
}

int rlab_sort_set_method(int i)
{
  if (i == RLAB_SORT_QUICK || i==RLAB_SORT_HEAP)
    _rlab_settings_sort_method = i;

  return _rlab_settings_sort_method;
}

int rlab_sort_get_nans(void )
{
  return _rlab_settings_sort_put_nans;
}

int rlab_sort_set_nans(int i)
{
  if (i == RLAB_SORT_NAN_INPLACE || i==RLAB_SORT_NAN_ONBOTTOM || i==RLAB_SORT_NAN_ONTOP)
    _rlab_settings_sort_put_nans = i;

  return _rlab_settings_sort_put_nans;
}

static void swap2 (double *v, int i, int j, double *ind)
{
  if (i==j)
    return;
  SWAP (v[i],  v[j], double);
  SWAP (ind[i],  ind[j], double);
}

//
// sort function: HEAPSORT
//
static inline void downheap (double *data, const int N, int k, double *ind)
{
  while (k <= N / 2)
  {
    int j = 2 * k;

    if ((j<N) && (data[j]<data[j + 1]))
      j++;

    if (data[k]<data[j])
      swap2 (data, k, j, ind);
    else
      break;

    k = j;
  }
}

static void r_hsort (double *v, int l, int r, double *ind)
{
  int k;

  switch (_rlab_settings_sort_put_nans)
  {
    case RLAB_SORT_NAN_ONBOTTOM:
      // do one pass through the array and presort nans to the top of the array
      while (isnand(v[l]) && l<=r)
        l++;
      if (l==r)
        goto _exit_heap;
      for (k=l;k<=r;k++)
      {
        if (isnand(v[k]))
        {
          swap2 (v, l, k, ind);
          l++;
        }
      }
      break;

    case RLAB_SORT_NAN_ONTOP:
      // do one pass through the array and presort nans to the top of the array
      while (isnand(v[r]) && l<=r)
        r--;
      if (l==r)
        goto _exit_heap;
      for (k=r;k>=l;k--)
      {
        if (isnand(v[k]))
        {
          swap2 (v, r, k, ind);
          r--;
        }
      }
      break;
  }

  int count =  r - l + 1;

  if (count < 2)
    return;

  //
  // with what is left do the heapsort
  //
  double *data  = &v[l];
  int N = count - 1;

  k = N / 2;
  k++;                          /* Compensate the first use of 'k--' */
  do
  {
    k--;
    downheap (data, N, k, ind);
  }
  while (k > 0);

  while (N > 0)
  {
    /* first swap the elements */
    swap2 (data, 0, N, ind);

    /* then process the heap */
    N--;

    downheap (data, N, 0, ind);
  }

_exit_heap:

  return;
}

//
// sort function: QUICKSORT
//
static void r_qsort_nan_bottom (double *v, int l, int r, double *ind)
{
  int i, last, left=l, right=r;

  if (left >= right)		/* Do nothing if array contains */
    return;			/* fewer than two elements */

  // put nan's at the beginning of the array
  swap2 (v, left, (left + right) / 2, ind); /* Move partitiion element */
  last = left;      /* to v[0] */
  for (i = left + 1; i <= right; i++) /* Partition */
  {
    if (isnand(v[left]))
      break;
    if (v[i] < v[left] || isnand(v[i]))
      swap2 (v, ++last, i, ind);
  }
  swap2 (v, left, last, ind); /* Restore partition element */

  // continue with sorting
  r_qsort_nan_bottom (v, left, last - 1, ind);
  r_qsort_nan_bottom (v, last + 1, right, ind);
}

static void r_qsort_nan_ontop (double *v, int l, int r, double *ind)
{
  int i, last, left=l, right=r;

  if (left >= right)    /* Do nothing if array contains */
    return;     /* fewer than two elements */

  // put nan's at the end of the array
  swap2 (v, left, (left + right) / 2, ind);
  last = left;
  for (i = left + 1; i <= right; i++)
  {
    if (isnand(v[i]))
      continue;
    if (v[i] < v[left] || isnand(v[left]))
      swap2 (v, ++last, i, ind);
  }
  swap2 (v, left, last, ind);

  // continue with sorting
  r_qsort_nan_ontop (v, left, last - 1, ind);
  r_qsort_nan_ontop (v, last + 1, right, ind);
}

static void r_bsort_nan_inplace (double *v, int l, int r, double *ind)
{
  int i, j, left=l, right=r;

  if (left >= right)    // Do nothing if array contains
    return;             // fewer than two elements

  // keep nan's in place: cannot do it with quicksort
  //  lets be bubbly
  for (i=left+1; i<=right; i++)
  {
    if (isnand(v[i]))
      continue;

    for (j=0; j<i; j++)
    {
      if (isnand(v[j]))
        continue;
      if (v[i] < v[j])
        swap2 (v, i, j, ind);
    }
  }
}

void r_sort (double *v, int l, int r, double *ind)
{
  if (_rlab_settings_sort_put_nans == RLAB_SORT_NAN_INPLACE)
  {
    r_bsort_nan_inplace (v, l, r, ind); // bubbly
  }
  else if (_rlab_settings_sort_method == RLAB_SORT_HEAP)
  {
    r_hsort (v, l, r, ind);             // heapsort
  }
  else if (_rlab_settings_sort_method == RLAB_SORT_QUICK)
  {
    // quicksort
    switch (_rlab_settings_sort_put_nans)
    {
      case RLAB_SORT_NAN_ONBOTTOM:
        r_qsort_nan_bottom (v, l, r, ind);
        break;

      case RLAB_SORT_NAN_ONTOP:
        r_qsort_nan_ontop (v, l, r, ind);
        break;
    }
  }
}

//----------------------------------------------------------
//
// S O R T   I N T E G E R S
//
//----------------------------------------------------------

static void swap2_int (int *v, int i, int j, double *ind)
{
  if (i==j)
    return;
  SWAP (v[i],  v[j], int);
  SWAP (ind[i],  ind[j], double);
}

//
// sort function: QUICKSORT
//
static void qsort_int (int *v, int l, int r, double *ind)
{
  int i, last, left=l, right=r;

  if (left >= right)    /* Do nothing if array contains */
    return;     /* fewer than two elements */

  // put nan's at the beginning of the array
  swap2_int (v, left, (left + right) / 2, ind); /* Move partitiion element */
  last = left;      /* to v[0] */
  for (i = left + 1; i <= right; i++) /* Partition */
  {
    if (v[i] < v[left])
      swap2_int (v, ++last, i, ind);
  }
  swap2_int (v, left, last, ind); /* Restore partition element */

  // continue with sorting
  qsort_int (v, left, last - 1, ind);
  qsort_int (v, last + 1, right, ind);
}


static inline void downheap_int (int *data, const int N, int k, double *ind)
{
  while (k <= N / 2)
  {
    int j = 2 * k;

    if ((j<N) && (data[j]<data[j + 1]))
      j++;

    if (data[k]<data[j])
      swap2_int (data, k, j, ind);
    else
      break;

    k = j;
  }
}

static void hsort_int (int *v, int l, int r, double *ind)
{
  int k;
  int count =  r - l + 1;

  if (count < 2)
    return;

  //
  // with what is left do the heapsort
  //
  int *data  = &v[l];
  int N = count - 1;

  k = N / 2;
  k++;                          /* Compensate the first use of 'k--' */
  do
  {
    k--;
    downheap_int (data, N, k, ind);
  }
  while (k > 0);

  while (N > 0)
  {
    /* first swap the elements */
    swap2_int (data, 0, N, ind);

    /* then process the heap */
    N--;

    downheap_int (data, N, 0, ind);
  }

  return;
}

void i_sort (int *v, int l, int r, double *ind)
{
  if (_rlab_settings_sort_method == RLAB_SORT_HEAP)
  {
    hsort_int (v, l, r, ind);   // heapsort
  }
  else if (_rlab_settings_sort_method == RLAB_SORT_QUICK)
  {
    qsort_int (v, l, r, ind);   // quicksort
  }
}


//
// Simple character qsort.
//

// Interchange v[i] and v[j]
static void swap2_char (char **v, int i, int j, double *ind)
{
  if (i==j)
    return;
  SWAP(v[i],v[j],char *);
  SWAP(ind[i],ind[j],double);
}

void
csort (char *v[], int left, int right, double *ind)
{
  int i, last;

  if (left >= right)
    return;

  swap2_char (v, left, (left + right) / 2, ind);

  last = left;
  for (i = left + 1; i <= right; i++)
    if (strcmp (v[i], v[left]) < 0)
      swap2_char (v, ++last, i, ind);

  swap2_char (v, left, last, ind);

  csort (v, left, last - 1, ind);
  csort (v, last + 1, right, ind);
}

//
// integer sort for sparse matrices
//
void
i_qsort (int *v, int left, int right)
{
  int i, last;

  if (left >= right)		  /* Do nothing if array contains */
    return;			          /* fewer than two elements */

  SWAP(v[left],v[(left + right)>>1],int);   // iswap (v, left, (left + right) / 2);	// Move partitiion element */
  last = left;			                        // to v[0] */

  for (i = left + 1; i <= right; i++)       // Partition
    if (v[i] < v[left])
    {
      last++;
      SWAP(v[last],v[i],int);               // iswap (v, ++last, i);
    }

  SWAP(v[left],v[last],int);                // iswap (v, left, last);	/* Restore partition element */
  i_qsort (v, left, last - 1);
  i_qsort (v, last + 1, right);
}




