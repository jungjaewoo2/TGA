/* This file is released as a part of rlabplus project, see rlabplus.sourceforge.net
   Copyright (C) 2006,2014 M. Kostrun

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

/* The C clustering library.
 * Copyright (C) 2002 Michiel Jan Laurens de Hoon.
 *
 * This library was written at the Laboratory of DNA Information Analysis,
 * Human Genome Center, Institute of Medical Science, University of Tokyo,
 * 4-6-1 Shirokanedai, Minato-ku, Tokyo 108-8639, Japan.
 * Contact: mdehoon 'AT' gsc.riken.jp
 *
 * Permission to use, copy, modify, and distribute this software and its
 * documentation with or without modifications and for any purpose and
 * without fee is hereby granted, provided that any copyright notices
 * appear in all copies and that both those copyright notices and this
 * permission notice appear in supporting documentation, and that the
 * names of the contributors or copyright holders not be used in
 * advertising or publicity pertaining to distribution of the software
 * without specific prior permission.
 *
 * THE CONTRIBUTORS AND COPYRIGHT HOLDERS OF THIS SOFTWARE DISCLAIM ALL
 * WARRANTIES WITH REGARD TO THIS SOFTWARE, INCLUDING ALL IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS, IN NO EVENT SHALL THE
 * CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY SPECIAL, INDIRECT
 * OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS
 * OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE
 * OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE
 * OR PERFORMANCE OF THIS SOFTWARE.
 *
 */



// ------------------------------------------------------- //
//
//
// DISTANCE MATRIX
//
//
// ------------------------------------------------------- //

/* Calculates the ranks of the elements in the array data. Two elements with
 * the same value get the same rank, equal to the average of the ranks had the
 * elements different values. The ranks are returned as a newly allocated
 * array that should be freed by the calling routine. If getrank fails due to
 * a memory allocation error, it returns NULL.
 */
static void compute_rank (int n, double *data, double *rank, double *index, int offset)
{
  int i;

  for (i=0;i<n;i++)
    index[i]=i;

  // Call sort to get an index table: rank contains sorted data, but we do not care
  r_sort (rank, 0, n-1, index);

  // Build a rank table
  for (i=0; i<n; i++)
    rank[(int) index[i]] = i + offset;

  // Fix for equal ranks
  i = 0;
  while (i < n)
  {
    int m;
    double value = data[(int) index[i]];
    int j = i + 1;
    while (j < n && data[(int) index[j]] == value)
      j++;

    m = j - i; /* number of equal ranks found */
    value = rank[(int) index[i]] + 0.5 * (m-1);

    for (j = i; j < i + m; j++)
      rank[(int) index[j]] = value;
    i += m;
  }

  return;
}

static double * getrank (int n, double *data)
{
  double *rank=GC_malloc(n*sizeof(double));
  double *index=GC_malloc(n*sizeof(double));

  memcpy(rank,data,n*sizeof(double));

  compute_rank(n, data, rank, index, 0);

  GC_FREE (index);
  return rank;
}

//
// determine datum rank in a data vector or matrix
//
#undef THIS_SOLVER
#define THIS_SOLVER "rank"
Ent *
ent_rank (int nargs, Datum args[])
{

  Ent *e1=0;
  MDR *x1=0, *w=0;

  int j,n;
  double *index=0;
  FILE *rlab_stderr = (!RLAB_STDERR_DS) ? stderr : RLAB_STDERR_DS;

  //
  // Load and Check arguments.
  //
  if (nargs != 1)
  {
    fprintf(rlab_stderr,
            THIS_SOLVER ": Determine rank of data vector or matrix.\n");
    fprintf(rlab_stderr,
            THIS_SOLVER ": Format :\n");
    fprintf(rlab_stderr,
            THIS_SOLVER ":   r = rank(x)\n");
    fprintf(rlab_stderr,
            THIS_SOLVER ": where 'x' is the data.\n");
    rerror
        (   THIS_SOLVER ": " RLAB_ERROR_ONE_ARG_REQUIRED);
  }

  //
  // x
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) != MATRIX_DENSE_REAL)
    rerror (THIS_SOLVER ": 'x' has to be a real matrix!");
  x1 = class_matrix_real (e1);

  if (SIZE(x1) < 2)
    goto _getrank_exit;

  w = mdr_Copy(x1);

  // storage for sorted indices
  if (EQVECT(x1))
    n = SIZE(x1);
  else
    n = MNR(x1);
  index = GC_malloc(n*sizeof(double));

  if (EQVECT(x1))
    compute_rank(n, MDRPTR(x1), MDRPTR(w), index, 1);
  else
  {
    for (j=0; j<MNC(x1); j++)
      compute_rank(n, &Mdr0(x1,0,j), &Mdr0(w,0,j), index, 1);
  }

  if (index)
    GC_FREE (index);

_getrank_exit:

  ent_Clean (e1);
  return ent_Assign_Rlab_MDR(w);
}


/*
 * Purpose
 * =======
 *
 * The spearman routine calculates the Spearman distance between two rows or
 * columns. The Spearman distance is defined as one minus the Spearman rank
 * correlation.
 *
 */
static double
distance_spearman (int rowdominant, MDR* data1, MDR* data2, MDR * w, int idx1, int idx2)
{
  int i, nattr;
  int m=0;
  double *rank1=0, *tdata1=0;
  double *rank2=0, *tdata2=0;
  double result = 0.;
  double denom1 = 0.;
  double denom2 = 0.;
  double avgrank;

  if (rowdominant)
  {
    nattr = MNC(data1);
    tdata1 = (double *) GC_malloc(nattr * sizeof(double));
    tdata2 = (double *) GC_malloc(nattr * sizeof(double));

    for (i=0; i<nattr; i++)
    {
      if (isnand(Mdr0(data1,idx1,i)) || isnand(Mdr0(data2,idx2,i)))
        continue;
      tdata1[m] = Mdr0(data1,idx1,i);
      tdata2[m] = Mdr0(data2,idx2,i);
      m++;
    }
  }
  else
  {
    nattr = MNR(data1);
    tdata1 = (double *) GC_malloc(nattr * sizeof(double));
    tdata2 = (double *) GC_malloc(nattr * sizeof(double));

    for (i=0; i<nattr; i++)
    {
      if (isnand(Mdr0(data1,i,idx1)) || isnand(Mdr0(data2,i,idx2)))
        continue;
      tdata1[m] = Mdr0(data1,i,idx1);
      tdata2[m] = Mdr0(data2,i,idx2);
      m++;
    }
  }

  if (!m)
  {
    GC_FREE(tdata1);
    GC_FREE(tdata2);
    return 0;
  }

  rank1 = (double *) getrank(m, tdata1);
  GC_FREE(tdata1);

  rank2 = getrank(m, tdata2);
  GC_FREE(tdata2);

  avgrank = 0.5*(m-1); /* Average rank */

  for (i = 0; i < m; i++)
  {
    const double value1 = rank1[i];
    const double value2 = rank2[i];
    result += value1 * value2;
    denom1 += value1 * value1;
    denom2 += value2 * value2;
  }
  /* Note: denom1 and denom2 cannot be calculated directly from the number
   * of elements. If two elements have the same rank, the squared sum of
   * their ranks will change.
   */
  GC_FREE(rank1);
  GC_FREE(rank2);

  result /= m;
  denom1 /= m;
  denom2 /= m;

  result -= avgrank * avgrank;
  denom1 -= avgrank * avgrank;
  denom2 -= avgrank * avgrank;

  if (denom1 <= 0)
    return 1; /* include '<' to deal with roundoff errors */

  if (denom2 <= 0)
    return 1; /* include '<' to deal with roundoff errors */

  result = result / sqrt(denom1*denom2);

  result = 1. - result;
  return result;
}

/*
 *    Purpose
 *    =======
 *
 *    The euclid routine calculates the weighted Euclidean distance between two
 *    rows or columns in a matrix.
 *
 */
static double
distance_euclid (int rowdominant, MDR* data1, MDR* data2, MDR * w, int idx1, int idx2)
{
  double result = 0.;
  double term;
  int i, nattr;

  if (rowdominant) /* Calculate the distance between two rows */
  {
    nattr = MNC(data1);

    for (i=0; i<nattr; i++)
    {
      if (isnand(mdr0(data1,idx1,i)) || isnand(mdr0(data2,idx2,i)))
        continue;

      term = mdr0(data1,idx1,i) - mdr0(data2,idx2,i);
      if (w)
        result  += MdrV0(w,i)*term*term;
      else
        result += term*term;
    }
  }
  else
  {
    nattr = MNR(data1);

    for (i=0; i<nattr; i++)
    {
      if (isnand(mdr0(data1,i,idx1)) || isnand(mdr0(data2,i,idx2)))
        continue;

      term = mdr0(data1,i,idx1) - mdr0(data2,i,idx2);
      if (w)
        result  += MdrV0(w,i)*term*term;
      else
        result += term*term;
    }
  }

  return result;
}

/*
 *  Purpose
 *  =======
 *
 *  The cityblock routine calculates the weighted "City Block" distance between
 *  two rows or columns in a matrix. City Block distance is defined as the
 *  absolute value of X1-X2 plus the absolute value of Y1-Y2 plus..., which is
 *  equivalent to taking an "up and over" path.
 */
static double
distance_cityblock (int rowdominant, MDR* data1, MDR* data2, MDR * w, int idx1, int idx2)
{
  double result = 0.;
  double term;
  int i, nattr;

  if (rowdominant) /* Calculate the distance between two rows */
  {
    nattr = MNC(data1);

    for (i=0; i<nattr; i++)
    {
      if (isnand(mdr0(data1,idx1,i)) || isnand(mdr0(data2,idx2,i)))
        continue;

      term = ABS(mdr0(data1,idx1,i) - mdr0(data2,idx2,i));
      if (w)
        result  += MdrV0(w,i)*term;
      else
        result += term;
    }
  }
  else
  {
    nattr = MNR(data1);

    for (i=0; i<nattr; i++)
    {
      if (isnand(mdr0(data1,i,idx1)) || isnand(mdr0(data2,i,idx2)))
        continue;

      term = ABS(mdr0(data1,i,idx1) - mdr0(data2,i,idx2));
      if (w)
        result  += MdrV0(w,i)*term;
      else
        result += term;
    }
  }

  return result;
}

/*
 *  Purpose
 *  =======
 *
 *  The cityblock_max routine calculates the weighted "City Block" distance between
 *  two rows or columns in a matrix. City Block distance is defined as the
 *  absolute value of X1-X2 plus the absolute value of Y1-Y2 plus..., which is
 *  equivalent to taking an "up and over" path.
 */
static double
    distance_cityblock_max (int rowdominant, MDR* data1, MDR* data2, MDR * w, int idx1, int idx2)
{
  double result=0.0;
  double term;
  int i, nattr;

  if (rowdominant) /* Calculate the distance between two rows */
  {
    nattr = MNC(data1);

    for (i=0; i<nattr; i++)
    {
      if (isnand(mdr0(data1,idx1,i)) || isnand(mdr0(data2,idx2,i)))
        continue;

      term = ABS(mdr0(data1,idx1,i) - mdr0(data2,idx2,i));
      if (w)
      {
        if (MdrV0(w,i) > 0)
        {
          term  *= MdrV0(w,i)*term;
        }
      }
      if (term > result)
        result = term;
    }
  }
  else
  {
    nattr = MNR(data1);

    for (i=0; i<nattr; i++)
    {
      if (isnand(mdr0(data1,i,idx1)) || isnand(mdr0(data2,i,idx2)))
        continue;

      term = ABS(mdr0(data1,i,idx1) - mdr0(data2,i,idx2));
      if (w)
      {
        if (MdrV0(w,i) > 0)
        {
          term  *= MdrV0(w,i)*term;
        }
      }
      if (term > result)
        result = term;
    }
  }

  return result;
}

/*
 * Purpose
 * =======
 *
 * The correlation routine calculates the weighted Pearson distance between two
 * rows or columns in a matrix. We define the Pearson distance as one minus the
 * Pearson correlation.
 * This definition yields a semi-metric: d(a,b) >= 0, and d(a,b) = 0 iff a = b.
 * but the triangular inequality d(a,b) + d(b,c) >= d(a,c) does not hold
 * (e.g., choose b = a + c).
 *
 * ============================================================================
 */
static double
distance_pearson (int rowdominant, MDR* data1, MDR* data2, MDR * w, int idx1, int idx2)
{
  double result = 0.;
  double sum1 = 0.;
  double sum2 = 0.;
  double denom1 = 0.;
  double denom2 = 0.;
  double term1, term2;
  int i, nattr;

  if (rowdominant) /* Calculate the distance between two rows */
  {
    nattr = MNC(data1);

    for (i=0; i<nattr; i++)
    {
      if (isnand(mdr0(data1,idx1,i)) || isnand(mdr0(data2,idx2,i)))
        continue;

      term1 = mdr0(data1,idx1,i);
      term2 = mdr0(data2,idx2,i);

      if (w)
      {
        sum1 += mdrV0(w,i)*term1;
        sum2 += mdrV0(w,i)*term2;
        result += mdrV0(w,i)*term1*term2;
        denom1 += mdrV0(w,i)*term1*term1;
        denom2 += mdrV0(w,i)*term2*term2;
      }
      else
      {
        sum1 += term1;
        sum2 += term2;
        result += term1*term2;
        denom1 += term1*term1;
        denom2 += term2*term2;
      }
    }
  }
  else
  {
    nattr = MNR(data1);

    for (i=0; i<nattr; i++)
    {
      if (isnand(mdr0(data1,i,idx1)) || isnand(mdr0(data2,i,idx2)))
        continue;

      term1 = mdr0(data1,i,idx1);
      term2 = mdr0(data2,i,idx2);

      if (w)
      {
        sum1 += mdrV0(w,i)*term1;
        sum2 += mdrV0(w,i)*term2;
        result += mdrV0(w,i)*term1*term2;
        denom1 += mdrV0(w,i)*term1*term1;
        denom2 += mdrV0(w,i)*term2*term2;
      }
      else
      {
        sum1 += term1;
        sum2 += term2;
        result += term1*term2;
        denom1 += term1*term1;
        denom2 += term2*term2;
      }
    }
  }

  result -= sum1 * sum2;
  denom1 -= sum1 * sum1;
  denom2 -= sum2 * sum2;
  if (denom1 <= 0)
    return 1;
  if (denom2 <= 0)
    return 1;

  result = result / sqrt(denom1*denom2);

  result = 1.0 - result;
  return result;
}

/*
 * Purpose
 * =======
 *
 * The acorrelation routine calculates the weighted Pearson distance between two
 * rows or columns, using the absolute value of the correlation.
 * This definition yields a semi-metric: d(a,b) >= 0, and d(a,b) = 0 iff a = b.
 * but the triangular inequality d(a,b) + d(b,c) >= d(a,c) does not hold
 * (e.g., choose b = a + c).
 *
 *
 * ============================================================================
 */
static double
distance_pearson_abs (int rowdominant, MDR* data1, MDR* data2, MDR * w, int idx1, int idx2)
{
  double result = 0.;
  double sum1 = 0.;
  double sum2 = 0.;
  double denom1 = 0.;
  double denom2 = 0.;
  double term1, term2;
  int i, nattr;

  if (rowdominant) /* Calculate the distance between two rows */
  {
    nattr = MNC(data1);

    for (i=0; i<nattr; i++)
    {
      if (isnand(mdr0(data1,idx1,i)) || isnand(mdr0(data2,idx2,i)))
        continue;

      term1 = mdr0(data1,idx1,i);
      term2 = mdr0(data2,idx2,i);

      if (w)
      {
        sum1 += mdrV0(w,i)*term1;
        sum2 += mdrV0(w,i)*term2;
        result += mdrV0(w,i)*term1*term2;
        denom1 += mdrV0(w,i)*term1*term1;
        denom2 += mdrV0(w,i)*term2*term2;
      }
      else
      {
        sum1 += term1;
        sum2 += term2;
        result += term1*term2;
        denom1 += term1*term1;
        denom2 += term2*term2;
      }
    }
  }
  else
  {
    nattr = MNR(data1);

    for (i=0; i<nattr; i++)
    {
      if (isnand(mdr0(data1,i,idx1)) || isnand(mdr0(data2,i,idx2)))
        continue;

      term1 = mdr0(data1,i,idx1);
      term2 = mdr0(data2,i,idx2);

      if (w)
      {
        sum1 += mdrV0(w,i)*term1;
        sum2 += mdrV0(w,i)*term2;
        result += mdrV0(w,i)*term1*term2;
        denom1 += mdrV0(w,i)*term1*term1;
        denom2 += mdrV0(w,i)*term2*term2;
      }
      else
      {
        sum1 += term1;
        sum2 += term2;
        result += term1*term2;
        denom1 += term1*term1;
        denom2 += term2*term2;
      }
    }
  }

  result -= sum1 * sum2;
  denom1 -= sum1 * sum1;
  denom2 -= sum2 * sum2;
  if (denom1 <= 0)
    return 1;
  if (denom2 <= 0)
    return 1;

  result = ABS(result) / sqrt(denom1*denom2);

  result = 1.0 - result;
  return result;
}

/*
 * Purpose
 * =======
 *
 * The ucorrelation routine calculates the weighted Pearson distance between two
 * rows or columns, using the uncentered version of the Pearson correlation. In the
 * uncentered Pearson correlation, a zero mean is used for both vectors even if
 * the actual mean is nonzero.
 * This definition yields a semi-metric: d(a,b) >= 0, and d(a,b) = 0 iff a = b.
 * but the triangular inequality d(a,b) + d(b,c) >= d(a,c) does not hold
 * (e.g., choose b = a + c).
 *
 *
 * ============================================================================
 */
static double
distance_pearson_uncent (int rowdominant, MDR* data1, MDR* data2, MDR * w, int idx1, int idx2)
{
  double result = 0.;
  double sum1 = 0.;
  double sum2 = 0.;
  double denom1 = 0.;
  double denom2 = 0.;
  double term1, term2;
  int i, nattr, flag=0;

  if (rowdominant) /* Calculate the distance between two rows */
  {
    nattr = MNC(data1);

    for (i=0; i<nattr; i++)
    {
      if (isnand(mdr0(data1,idx1,i)) || isnand(mdr0(data2,idx2,i)))
        continue;

      term1 = mdr0(data1,idx1,i);
      term2 = mdr0(data2,idx2,i);

      if (w)
      {
        sum1 += mdrV0(w,i)*term1;
        sum2 += mdrV0(w,i)*term2;
        result += mdrV0(w,i)*term1*term2;
        denom1 += mdrV0(w,i)*term1*term1;
        denom2 += mdrV0(w,i)*term2*term2;
      }
      else
      {
        sum1 += term1;
        sum2 += term2;
        result += term1*term2;
        denom1 += term1*term1;
        denom2 += term2*term2;
      }
      flag=1;
    }
  }
  else
  {
    nattr = MNR(data1);

    for (i=0; i<nattr; i++)
    {
      if (isnand(mdr0(data1,i,idx1)) || isnand(mdr0(data2,i,idx2)))
        continue;

      term1 = mdr0(data1,i,idx1);
      term2 = mdr0(data2,i,idx2);

      if (w)
      {
        sum1 += mdrV0(w,i)*term1;
        sum2 += mdrV0(w,i)*term2;
        result += mdrV0(w,i)*term1*term2;
        denom1 += mdrV0(w,i)*term1*term1;
        denom2 += mdrV0(w,i)*term2*term2;
      }
      else
      {
        sum1 += term1;
        sum2 += term2;
        result += term1*term2;
        denom1 += term1*term1;
        denom2 += term2*term2;
      }
      flag=1;
    }
  }

  if (!flag)
    return 0.;
  if (denom1 <= 0)
    return 1;
  if (denom2 <= 0)
    return 1;

  result = result / sqrt(denom1*denom2);

  result = 1.0 - result;
  return result;
}


/*
 * Purpose
 * =======
 *
 * The uacorrelation routine calculates the weighted Pearson distance between two
 * rows or columns, using the absolute value of the uncentered version of the
 * Pearson correlation. In the uncentered Pearson correlation, a zero mean is used
 * for both vectors even if the actual mean is nonzero.
 * This definition yields a semi-metric: d(a,b) >= 0, and d(a,b) = 0 iff a = b.
 * but the triangular inequality d(a,b) + d(b,c) >= d(a,c) does not hold
 * (e.g., choose b = a + c).
 */
static double
distance_pearson_uncent_abs (int rowdominant, MDR* data1, MDR* data2, MDR * w, int idx1, int idx2)
{
  double result = 0.;
  double sum1 = 0.;
  double sum2 = 0.;
  double denom1 = 0.;
  double denom2 = 0.;
  double term1, term2;
  int i, nattr, flag=0;

  if (rowdominant) /* Calculate the distance between two rows */
  {
    nattr = MNC(data1);

    for (i=0; i<nattr; i++)
    {
      if (isnand(mdr0(data1,idx1,i)) || isnand(mdr0(data2,idx2,i)))
        continue;

      term1 = mdr0(data1,idx1,i);
      term2 = mdr0(data2,idx2,i);

      if (w)
      {
        sum1 += mdrV0(w,i)*term1;
        sum2 += mdrV0(w,i)*term2;
        result += mdrV0(w,i)*term1*term2;
        denom1 += mdrV0(w,i)*term1*term1;
        denom2 += mdrV0(w,i)*term2*term2;
      }
      else
      {
        sum1 += term1;
        sum2 += term2;
        result += term1*term2;
        denom1 += term1*term1;
        denom2 += term2*term2;
      }
      flag=1;
    }
  }
  else
  {
    nattr = MNR(data1);

    for (i=0; i<nattr; i++)
    {
      if (isnand(mdr0(data1,i,idx1)) || isnand(mdr0(data2,i,idx2)))
        continue;

      term1 = mdr0(data1,i,idx1);
      term2 = mdr0(data2,i,idx2);

      if (w)
      {
        sum1 += mdrV0(w,i)*term1;
        sum2 += mdrV0(w,i)*term2;
        result += mdrV0(w,i)*term1*term2;
        denom1 += mdrV0(w,i)*term1*term1;
        denom2 += mdrV0(w,i)*term2*term2;
      }
      else
      {
        sum1 += term1;
        sum2 += term2;
        result += term1*term2;
        denom1 += term1*term1;
        denom2 += term2*term2;
      }
      flag=1;
    }
  }

  if (!flag)
    return 0.;
  if (denom1 <= 0)
    return 1;
  if (denom2 <= 0)
    return 1;

  result = ABS(result) / sqrt(denom1*denom2);

  result = 1.0 - result;
  return result;
}


/*
 * Purpose
 * =======
 *
 * The kendall routine calculates the Kendall distance between two
 * rows or columns. The Kendall distance is defined as one minus Kendall's tau.
 */
static double
distance_kendall (int rowdominant, MDR* data1, MDR* data2, MDR * w, int idx1, int idx2)
{
  int con = 0;
  int dis = 0;
  int exx = 0;
  int exy = 0;
  int flag = 0;
  double denomx, denomy, tau, x1, x2, y1, y2;
  int i, j, nattr;

  if (rowdominant) /* Calculate the distance between two rows */
  {
    nattr = MNC(data1);

    for (i=0; i<nattr; i++)
    {
      if (isnand(mdr0(data1,idx1,i)) || isnand(mdr0(data2,idx2,i)))
        continue;

      for (j=0; j<i; j++)
      {
        if (isnand(mdr0(data1,idx1,j)) || isnand(mdr0(data2,idx2,j)))
          continue;

        x1 = mdr0(data1,idx1,i);
        x2 = mdr0(data1,idx1,j);
        y1 = mdr0(data2,idx1,i);
        y2 = mdr0(data2,idx1,j);

        if (x1 < x2 && y1 < y2)
          con++;
        if (x1 > x2 && y1 > y2)
          con++;
        if (x1 < x2 && y1 > y2)
          dis++;
        if (x1 > x2 && y1 < y2)
          dis++;
        if (x1 == x2 && y1 != y2)
          exx++;
        if (x1 != x2 && y1 == y2)
          exy++;
        flag = 1;
      }
    }
  }
  else
  {
    nattr = MNR(data1);

    for (i=0; i<nattr; i++)
    {
      if (isnand(mdr0(data1,i,idx1)) || isnand(mdr0(data2,i,idx2)))
        continue;

      for (j=0; j<i; j++)
      {
        if (isnand(mdr0(data1,j,idx1)) || isnand(mdr0(data2,j,idx2)))
          continue;

        x1 = mdr0(data1,i,idx1);
        x2 = mdr0(data1,j,idx1);
        y1 = mdr0(data2,i,idx1);
        y2 = mdr0(data2,j,idx1);

        if (x1 < x2 && y1 < y2)
          con++;
        if (x1 > x2 && y1 > y2)
          con++;
        if (x1 < x2 && y1 > y2)
          dis++;
        if (x1 > x2 && y1 < y2)
          dis++;
        if (x1 == x2 && y1 != y2)
          exx++;
        if (x1 != x2 && y1 == y2)
          exy++;
        flag = 1;
      }
    }
  }

  if (!flag) return 0.;
  denomx = con + dis + exx;
  denomy = con + dis + exy;
  if (denomx==0) return 1;
  if (denomy==0) return 1;
  tau = (con-dis)/sqrt(denomx*denomy);
  return 1.-tau;
}

static double(*distance_metric(unsigned char dist))(int rowdominant, MDR* data1, MDR* data2, MDR * w, int idx1, int idx2)
{
  switch(dist)
  {
    case 's':
      return &distance_spearman;

    case 'e':
      return &distance_euclid;

    case 'b':
      return &distance_cityblock;

    case 'm':
      return &distance_cityblock_max;

    case 'p':
      return &distance_pearson;

    case 'a':
      return &distance_pearson_abs;

    case 'u':
      return &distance_pearson_uncent;

    case 'x':
      return &distance_pearson_uncent_abs;

    case 'k':
      return &distance_kendall;

    default:
      return &distance_euclid;
  }

  return NULL;
}

//
// distance matrix
//
#undef THIS_SOLVER
#define THIS_SOLVER "distance"
Ent * ent_distance_matrix (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *e3=0, *e4=0;
  MDR *x=0, *y=0, *w=0;
  char *m=0;
  int nr, nc, i, j,issym=0;

  FILE *rlab_stderr = (!RLAB_STDERR_DS) ? stderr : RLAB_STDERR_DS;

  //
  // Load and Check arguments.
  //
  if (nargs != 4 && nargs != 3 && nargs != 2 && nargs != 1)
  {
    fprintf(rlab_stderr, THIS_SOLVER ": Calculates distance matrix in or between datasets\n");
    fprintf(rlab_stderr, THIS_SOLVER ": Format :\n");
    fprintf(rlab_stderr, THIS_SOLVER ":   d = distmat(data1/,data2/,m,w//)\n");
    fprintf(rlab_stderr, THIS_SOLVER
        ": where 'data1' and 'data2' are datasets with the same number of attributes (columns),\n");
    fprintf(rlab_stderr, THIS_SOLVER ": 'm' designates the metric\n");
    fprintf(rlab_stderr, THIS_SOLVER ":    with 'e' for euclidean (default), 'b' for city-block, 'm' for city-block-max,\n");
    fprintf(rlab_stderr, THIS_SOLVER ":    'p' for pearson, 'a' for absolute pearson, 'u' for uncentered pearson,\n");
    fprintf(rlab_stderr, THIS_SOLVER ":    'x' for uncentered absolute pearson, and 'k' for kendall-tau,\n");
    fprintf(rlab_stderr, THIS_SOLVER ": and 'w' is vector of weights of the attributes in distance contribution.\n");
    rerror (THIS_SOLVER ": " RLAB_ERROR_AT_MOST_THREE_ARG_REQUIRED);
  }

  //
  // first data set
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) != MATRIX_DENSE_REAL)
    rerror (THIS_SOLVER ": " RLAB_ERROR_ARG1_MDR_MATRIX);
  x = class_matrix_real (e1);

  //
  // second dataset or 'm'
  //
  if (nargs > 1)
  {
    e2 = bltin_get_ent (args[1]);
    if (ent_type (e2) == MATRIX_DENSE_REAL)
      y = class_matrix_real (e2);
    else if (ent_type (e2) == MATRIX_DENSE_STRING)
    {
      m = class_char_pointer (e2);
      if (m)
        if (strlen(m)!=1)
          m=0;
    }
    else
    {
      printf (THIS_SOLVER ": " RLAB_ERROR_ARG2_MDR_MATRIX "\n");
      rerror (THIS_SOLVER ": " RLAB_ERROR_ARG2_MDS_SCALAR);
    }
  }
  if (!y)
  {
    y = x;
    issym=1;
  }

  if(MNC(x)!=MNC(y))
    rerror (THIS_SOLVER ": " RLAB_ERROR_ARG1_ARG2_SAME_NUMBER_COLS);

  if (nargs > 2 && !m)
  {
    e3 = bltin_get_ent (args[2]);
    if (ent_type (e3) == MATRIX_DENSE_STRING)
    {
      m = class_char_pointer (e3);
      if (m)
      {
        if (isvalidstring(m)!=1)
          m=0;
      }
    }
    else if (ent_type (e3) == MATRIX_DENSE_REAL)
    {
      w = class_matrix_real(e3);
      if (w->type != RLAB_TYPE_DOUBLE)
        w = 0;
    }
    if (m)
      if (!strlen(m))
        m=0;
  }

  if (nargs > 3 && !w)
  {
    e4 = bltin_get_ent (args[3]);
    if (ent_type (e4) == MATRIX_DENSE_REAL)
    {
      w = class_matrix_real(e4);
      if (w->type != RLAB_TYPE_DOUBLE)
        w = 0;
    }
  }

  if (!m)
    m = "e";

  // Set the metric function as indicated by dist
  double (*distance)(int, MDR*, MDR*, MDR*, int, int) =  distance_metric(m[0]);

  // rows represent the data in x and y matrices
  nr = MNR(x);
  nc = MNR(y);

  MDR *r = mdr_Create(nr,nc);

  if (issym)
  {
    Mdr0(r,0,0) = 0.0;
    for (i=1; i<nr; i++)
    {
      Mdr0(r,i,i) = 0.0;
      for (j=0; j<i; j++)
      {
        Mdr0(r,j,i) = Mdr0(r,i,j) = distance(1, x, y, w, i, j) ;
      }
    }
  }
  else
  {
    for (i=0; i<nr; i++)
      for (j=0; j<nc; j++)
        Mdr0(r,i,j) = distance(1, x, y, w, i, j) ;
  }

  // clean-up
  ent_Clean (e1);
  ent_Clean (e2);
  ent_Clean (e3);
  ent_Clean (e4);

  return ent_Assign_Rlab_MDR(r);
}


// ------------------------------------------------------- //
//
//
// PAO CLUSTER
//
//
// ------------------------------------------------------- //
static void
cluster_pao_first_node (int *nclusters, MDR* b, MDR* dataset, MDR* counts, MDR *clusterid, MDR *cluster_cnt)
{
   int i;
   int nattr=MNC(b);

   for (i=0; i<nattr; i++)
   {
     Mdi0(cluster_cnt,0,i)=0;
     Mdr0(b,0,i) = Mdr0(dataset,0,i);
     if (isnand(Mdr0(b,0,i)))
       continue;
     Mdi0(cluster_cnt,0,i)=1; // count non-nan for each coordinate of the cluster center
   }

   *nclusters = 1;
   MdrV0(counts,0) = 1;
   MdrV0(clusterid,0) = (*nclusters);
}

static void
cluster_pao_form_new_node (int q, int *nclusters, MDR* b, MDR* dataset, MDR* counts,
                           MDR *clusterid, MDR *cluster_cnt)
{
   int i;
   int nattr=MNC(b);

   for (i=0; i<nattr; i++)
   {
     Mdi0(cluster_cnt,*nclusters,i)=0;
     Mdr0(b, *nclusters,i) = Mdr0(dataset,q,i);
     if (isnand(Mdr0(b,*nclusters,i)))
       continue;
     Mdi0(cluster_cnt,*nclusters,i)=1; // count non-nan for each coordinate of the cluster center
   }

   MdrV0(counts,*nclusters) = 1;
   (*nclusters)++;
   MdrV0(clusterid,q) = *nclusters;
}


static int
cluster_pao_compare_min_ed (int *nclusters, double *threshold, MDR* ed)
{
   int i, cluster_no = 0;
   double dmin = 10000.0;

   for (i=0; i<(*nclusters); i++)
   {
     if (MdrV0(ed,i) < dmin)
     {
       dmin = MdrV0(ed,i);
       cluster_no = i;
     }
   }

   if (MdrV0(ed,cluster_no) <= *threshold)
     return cluster_no;
   else
     return (-99);
}

static void
cluster_pao_update_wts (int cluster_no, int input_no, MDR *dataset, MDR* b, MDR* counts,
                        MDR *clusterid, MDR *cluster_cnt)
{
   int i;
   double n, m;
   int nattr = MNC(b);

   // nothing to do if attribute is nan in cluster center, and in new data point
   for (i=0; i<nattr; i++)
   {
     if (!isnand(Mdr0(dataset,input_no,i)))
     {
       // attribute does not exist in the cluster: copy it from the new datum
       if (isnand(Mdr0(b, cluster_no,i)))
       {
         Mdr0(b, cluster_no,i) = Mdr0(dataset,input_no,i);
         Mdi0(cluster_cnt,cluster_no,i) = 1;
         continue;
       }

       // update non-nan number of entries for particular attribute in the cluster
       n = Mdi0(cluster_cnt,cluster_no,i);
       m = n + 1.;
       Mdr0(b, cluster_no,i) = ((n / m * Mdr0(b, cluster_no,i)) + (1 / m * Mdr0(dataset,input_no,i)));
       Mdi0(cluster_cnt,cluster_no,i) += 1;
     }
   }

   MdrV0(clusterid,input_no) = cluster_no + 1;  // update id of the point added to the cluster
   MdrV0(counts,cluster_no) += 1.0; // update total number of points in the cluster
}

/*
 * CLUSTER_PAO : Perform the clustering on the data set. The found clusters
 * will be stored into a dataset returned by this function. The supplied
 * threshold control the generation of new cluster centers and a low
 * threshold value will allow easy generation of a new cluster. A larger
 * value of this threshold will prohibit the generation of new clusters.
 */
//     kmeanmed(nclusters, data, weight, 1, npass, *m, &cdata, &clusterid, &error, &counts, imean);
static void
gnu_cluster_pao(double* threshold, MDR* dataset, MDR *w, char *dist_metric,
                MDR** cdata, MDR **clusterid, MDR **fcounts)
{
  int ndata;
  int nattr;
  int nclusters = 0;

  MDR *ed=0, *b=0;
  MDR *counts=0, *cluster_cnt;

   // use euclidean metrics between the rows of the matrix
  double (*distance)(int, MDR*, MDR*, MDR*, int, int) =  distance_metric(*dist_metric);

  int i, j, q, cluster_no;

  nattr = MNC(dataset);
  ndata = MNR(dataset);

  ed        = mdr_Create(ndata,1);
  b         = mdr_Create(ndata,nattr);
  cluster_cnt = mdi_Create(ndata,nattr);
  counts    = mdr_Create(ndata,1);
 *clusterid = mdr_Create(ndata,1);

  // copy first dataset to b
  cluster_pao_first_node (&nclusters, b, dataset, counts, *clusterid, cluster_cnt);

  for (q=1; q<ndata; q++)
  {
    for (i=0; i<nclusters; i++)
      MdrV0(ed,i) = distance(1, b, dataset, w, i, q);

    cluster_no = cluster_pao_compare_min_ed (&nclusters, threshold, ed);

    if (cluster_no >= 0)
      cluster_pao_update_wts (cluster_no, q, dataset, b, counts, *clusterid, cluster_cnt);
    else
      cluster_pao_form_new_node (q, &nclusters, b, dataset, counts, *clusterid, cluster_cnt);
  }

  *cdata   = mdr_Create(nclusters,nattr);
  *fcounts = mdr_Create(nclusters,1);
  for (i=0; i < nclusters; i++)
  {
    MdrV0(*fcounts,i) = MdrV0(counts,i);
    for (j=0; j<nattr; j++)
      Mdr0(*cdata,i,j)  = Mdr0(b,i,j);
  }

  mdr_Destroy(b);
  mdr_Destroy(counts);
  mdr_Destroy(cluster_cnt);
  mdr_Destroy(ed);

  return;
}

//
// clustering according to pao
//
#undef THIS_SOLVER
#define THIS_SOLVER "cluster.pao"
Ent *
ent_cluster_pao (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *e3=0, *e4=0;
  char *m=0;
  MDR *x=0, *w=0;
  double d;

  FILE *rlab_stderr = (!RLAB_STDERR_DS) ? stderr : RLAB_STDERR_DS;

  //
  // Load and Check arguments.
  //
  if ((nargs != 2) && (nargs != 3) && (nargs != 4))
  {
    fprintf(rlab_stderr, THIS_SOLVER ": Perform clustering of the data according to Pao.\n");
    fprintf(rlab_stderr, THIS_SOLVER ": Format :\n");
    fprintf(rlab_stderr, THIS_SOLVER ":   c = cluster.pao(data,thresh/,m,w/)\n");
    fprintf(rlab_stderr, THIS_SOLVER ": where 'dist' is a maximal distance between the points in a cluster,\n");
    fprintf(rlab_stderr, THIS_SOLVER ": and 'm' designates the distance metric\n");
    fprintf(rlab_stderr, THIS_SOLVER ":    with 'e' for euclidean (default), 'b' for city-block, 'm' for city-block-max,\n");
    fprintf(rlab_stderr, THIS_SOLVER ":    'p' for pearson, 'a' for absolute pearson, 'u' for uncentered pearson,\n");
    fprintf(rlab_stderr, THIS_SOLVER ":    'x' for uncentered absolute pearson, and 'k' for kendall-tau,\n");
    fprintf(rlab_stderr, THIS_SOLVER ": and 'w' is vector of weights of the attributes in distance contribution.\n");
    fprintf(rlab_stderr, THIS_SOLVER ": The result is a matrix of cluster centres and their weights.\n");
    rerror (THIS_SOLVER ": requires two or three arguments");
  }

  //
  // x - data to be clustered
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) != MATRIX_DENSE_REAL)
    rerror (THIS_SOLVER ": 'X' has to be a real matrix!");
  x = class_matrix_real (e1);

  //
  // d - threshold distance
  //
  e2 = bltin_get_ent (args[1]);
  if (ent_type (e2) != MATRIX_DENSE_REAL)
    rerror (THIS_SOLVER ": 'd' has to be a real scalar!");
  d = class_double (e2);

  //
  // m - metric
  //
  if (nargs > 2)
  {
    e3 = bltin_get_ent (args[2]);
    if (ent_type (e3) == MATRIX_DENSE_STRING)
    {
      m = class_char_pointer (e3);
      if (isvalidstring(m)<1)
      {
        printf(THIS_SOLVER ": Invalid string for 'metrics'. It will be ignored!\n");
        m=0;
      }
    }
  }

  //
  // w - attribute weights for calculation of distance
  //
  if (nargs > 3)
  {
    e4 = bltin_get_ent (args[3]);
    if (ent_type (e4) != MATRIX_DENSE_REAL)
      rerror ("cluster.pao: 'W' has to be a real matrix!");
    w = class_matrix_real (e4);
    if (w->type != RLAB_TYPE_DOUBLE)
      w=0;
  }


  if (!m)
    m="e";

  MDR *cdata=0, *clusterid=0, *counts=0;

  gnu_cluster_pao(&d, x, w, m, &cdata, &clusterid, &counts);

  // clean-up
  ent_Clean (e1);
  ent_Clean (e2);
  ent_Clean (e3);
  ent_Clean (e4);

  //
  // write result to output list
  //
  Btree *bw = btree_Create ();

  // value: cluster centres
  install (bw, RLAB_NAME_GEN_VALUE, ent_Assign_Rlab_MDR(cdata));

  // size: cluster sizes
  install (bw, RLAB_NAME_GEN_SIZE,  ent_Assign_Rlab_MDR(counts));

  // size: feature of data (cluster id datum belongs to)
  install (bw, RLAB_NAME_STAT_FEATURE,  ent_Assign_Rlab_MDR(clusterid));

  //
  return ent_Assign_Rlab_BTREE(bw);
}



/*
 * Purpose
 * =======
 *
 * The randomassign routine performs an initial random clustering, needed for
 * k-means or k-median clustering. Elements (genes or microarrays) are randomly
 * assigned to clusters. The number of elements in each cluster is chosen
 * randomly, making sure that each cluster will receive at least one element.
 *
 */
extern int gslrngb_ ( int n, double p );

static void randomassign (int nclusters, MDR * clusterid)
{
  int i, j;
  int k = 0;
  double p;
  int ndata = MNR(clusterid) * MNC(clusterid);
  unsigned int n = ndata-nclusters;

/* Draw the number of elements in each cluster from a multinomial
 * distribution, reserving ncluster elements to set independently
 * in order to guarantee that none of the clusters are empty.
 */
  for (i = 0; i < nclusters-1; i++)
  {
    p = 1.0/(nclusters-i);
    j =  gslrngb_ (n, p);
    n -= j;
    j += k+1; /* Assign at least one element to cluster i */
    for ( ; k < j; k++)
      MdrV0(clusterid,k) = i;
  }

  /* Assign the remaining elements to the last cluster */
  for ( ; k < ndata; k++)
    MdrV0(clusterid,k) = i;

  /* Create a random permutation of the cluster assignments */
  for (i = 0; i < ndata; i++)
  {
    j = (int) (i + (ndata - i) * gslrngf_ ());
    k = MdrV0(clusterid,j);
    MdrV0(clusterid,j) = MdrV0(clusterid,i);
    MdrV0(clusterid,i) = k;
  }
  return;
}

static double
    median (int n, double x[])
/*
      Find the median of X(1), ... , X(N), using as much of the quicksort
      algorithm as is needed to isolate it.
      N.B. On exit, the array X is partially ordered.
      Based on Alan J. Miller's median.f90 routine.
*/

{ int i, j;
  int nr = n / 2;
  int nl = nr - 1;
  int even = 0;
  /* hi & lo are position limits encompassing the median. */
  int lo = 0;
  int hi = n-1;

  if (n==2*nr) even = 1;
  if (n<3)
  { if (n<1) return 0.;
  if (n == 1) return x[0];
  return 0.5*(x[0]+x[1]);
  }

  /* Find median of 1st, middle & last values. */
  do
  { int loop;
  int mid = (lo + hi)/2;
  double result = x[mid];
  double xlo = x[lo];
  double xhi = x[hi];
  if (xhi<xlo)
  { double temp = xlo;
  xlo = xhi;
  xhi = temp;
  }
  if (result>xhi) result = xhi;
  else if (result<xlo) result = xlo;
    /* The basic quicksort algorithm to move all values <= the sort key (XMED)
  * to the left-hand end, and all higher values to the other end.
    */
  i = lo;
  j = hi;
  do
  { while (x[i]<result) i++;
  while (x[j]>result) j--;
  loop = 0;
  if (i<j)
  { double temp = x[i];
  x[i] = x[j];
  x[j] = temp;
  i++;
  j--;
  if (i<=j) loop = 1;
  }
  } while (loop); /* Decide which half the median is in. */

  if (even)
  { if (j==nl && i==nr)
        /* Special case, n even, j = n/2 & i = j + 1, so the median is
    * between the two halves of the series.   Find max. of the first
    * half & min. of the second half, then average.
        */
  { int k;
  double xmax = x[0];
  double xmin = x[n-1];
  for (k = lo; k <= j; k++) xmax = MAX(xmax,x[k]);
  for (k = i; k <= hi; k++) xmin = MIN(xmin,x[k]);
  return 0.5*(xmin + xmax);
  }
  if (j<nl) lo = i;
  if (i>nr) hi = j;
  if (i==j)
  { if (i==nl) lo = nl;
  if (j==nr) hi = nr;
  }
  }
  else
  { if (j<nr) lo = i;
  if (i>nr) hi = j;
  /* Test whether median has been isolated. */
  if (i==j && i==nr) return result;
  }
  }
  while (lo<hi-1);

  if (even)
    return (0.5*(x[nl]+x[nr]));

  if (x[lo]>x[hi])
  { double temp = x[lo];
  x[lo] = x[hi];
  x[hi] = temp;
  }
  return x[nr];
}

/*
    Purpose
    =======

    The getclustermedians routine calculates the cluster centroids, given to which
    cluster each element belongs. The centroid is defined as the median over all
    elements for each dimension.
 */
static void
getclustermedians(int rowdominant, int nclusters, MDR * data, MDR * clusterid, MDR ** cdata )
{
  int i, j, k, ndata, nattr;
  double *cache=0;

  if (rowdominant)
  {
    nattr = MNC(data);
    ndata = MNR(data);

    // check if cdata exists and is proper size:
    //  if so reuse it
    //  otherwise destroy old, and create new one
    if (!(*cdata))
      *cdata = (MDR*) mdr_Create (nclusters, nattr);
    else
    {
      if ((MNC(*cdata)!=nattr) || (MNR(*cdata)!=nclusters))
      {
        mdr_Destroy(*cdata);
        *cdata = (MDR*) mdr_Create (nclusters, nattr);
      }
    }

    cache = (double *) GC_malloc(ndata * sizeof(double));

    for (i=0; i<nclusters; i++)
    {
      for (j=0; j<nattr; j++)
      {
        Mdr0(*cdata,i,j) = create_nan();

        int count = 0;
        for (k = 0; k < ndata; k++)
        {
          if (i==MdrV0(clusterid,k) && !isnand(Mdr0(data,k,j)))
          {
            cache[count] = Mdr0(data,k,j);
            count++;
          }
        }
        if (count>0)
          Mdr0(*cdata,i,j) = median(count,cache);
      }
    }
  }
  else
  {
    nattr = MNR(data);
    ndata = MNC(data);

    // check if cdata exists and is proper size:
    //  if so reuse it
    //  otherwise destroy old, and create new one
    if (!(*cdata))
      *cdata = (MDR*) mdr_Create (nattr,nclusters);
    else
    {
      if ((MNR(*cdata)!=nattr) || (MNC(*cdata)!=nclusters))
      {
        mdr_Destroy(*cdata);
        *cdata = (MDR*) mdr_Create (nattr,nclusters);
      }
    }

    cache = (double *) GC_malloc(ndata * sizeof(double));

    for (i=0; i<nclusters; i++)
    {
      for (j=0; j<nattr; j++)
      {
        Mdr0(*cdata,j,i) = create_nan();

        int count = 0;
        for (k = 0; k < ndata; k++)
        {
          if (i==MdrV0(clusterid,k) && !isnand(Mdr0(data,j,k)))
          {
            cache[count] = Mdr0(data,j,k);
            count++;
          }
        }
        if (count>0)
          Mdr0(*cdata,j,i) = median(count,cache);
      }
    }
  }

  if (cache)
    GC_FREE(cache);

  return;
}

/*
 * Purpose
 * =======
 *
 * The getclustermeans routine calculates the cluster centroids, given to which
 * cluster each element belongs. The centroid is defined as the mean over all
 * elements for each dimension.
 *
 */
static void
getclustermeans(int rowdominant, int nclusters, MDR * dataset, MDR * clusterid, MDR ** cdata)
{
  int i, j, k;
  int ndata, nattr;

  MDR *cmask=0;

  if (rowdominant)
  {
    nattr = MNC(dataset);
    ndata = MNR(dataset);

    cmask = (MDR*) mdi_Create (nclusters, nattr);

    for (i=0; i < nclusters; i++)
      for (j = 0; j < nattr; j++)
      {
        Mdr0(*cdata,i,j) = create_nan();
        Mdi0( cmask,i,j) = 0;
      }

    for (k=0; k<ndata; k++)
    {
      i = MdrV0(clusterid,k);
      for (j=0; j<nattr; j++)
      {
        if (isnand(Mdr0(dataset,k,j)))
          continue;
        if (isnand(Mdr0(*cdata,i,j)))
        {
          Mdr0(*cdata,i,j) = Mdr0(dataset,k,j);
          Mdi0( cmask,i,j) = 1;
        }
        else
        {
          Mdr0(*cdata,i,j) += Mdr0(dataset,k,j);
          Mdi0( cmask,i,j) += 1;
        }
      }
    }

    for (i=0; i<nclusters; i++)
      for (j=0; j<nattr; j++)
      {
        if (Mdi0(cmask,i,j))
          Mdr0(*cdata,i,j) /= (double) Mdi0(cmask,i,j);
      }
  }
  else
  {
    nattr = MNR(dataset);
    ndata = MNC(dataset);

    cmask = (MDR*) mdi_Create (nattr, nclusters);

    for (i=0; i < nclusters; i++)
      for (j = 0; j < nattr; j++)
      {
        Mdr0(*cdata,j,i) = create_nan();
        Mdi0( cmask,j,i) = 0;
      }

    for (k=0; k<ndata; k++)
    {
      i = MdrV0(clusterid,k);
      for (j=0; j<nattr; j++)
      {
        if (isnand(Mdr0(dataset,k,j)))
          continue;
        if (isnand(Mdr0(*cdata,j,i)))
        {
          Mdr0(*cdata,j,i) = Mdr0(dataset,j,k);
          Mdi0( cmask,j,i) = 1;
        }
        else
        {
          Mdr0(*cdata,j,i) += Mdr0(dataset,j,k);
          Mdi0( cmask,j,i) += 1;
        }
      }
    }

    for (i=0; i<nclusters; i++)
      for (j=0; j<nattr; j++)
      {
        if (Mdi0(cmask,j,i))
          Mdr0(*cdata,j,i) /= (double) Mdi0(cmask,j,i);
      }
  }

  if (cmask)
    mdr_Destroy(cmask);
}

/* ---------------------------------------------------------------------- */
#undef  THIS_SOLVER
#define THIS_SOLVER "kmeanmed"
static int kmeanmed ( int nclusters, MDR * data, MDR * weight, int rowdominant, 
  int npass, char dist, MDR ** cdata, MDR ** clusterid, double* error, MDR ** counts, int imean)
{
  int i, j, k;
  const int ndata = (rowdominant==1) ? MNR(data) : MNC(data);
  int ifound = 1;
  int ipass = 0;

//   dprintf("in\n");

  MDR *tclusterid = mdr_Float_BF(*clusterid);

  *counts = mdr_Create(nclusters,1);

  /* Set the metric function as indicated by dist */
  double (*metric)(int, MDR*, MDR*, MDR*, int, int) =  distance_metric(dist);

  /* We save the clustering solution periodically and check if it reappears */
  int *  saved = GC_malloc(ndata*sizeof(int));
//   int *mapping = GC_malloc(nclusters*sizeof(int));

  *error = DBL_MAX;

  do
  {
    double total = DBL_MAX;
    int counter = 0;
    int period = 10;

    /* Perform the EM algorithm. First, randomly assign elements to clusters. */
    if (npass!=0)
      randomassign (nclusters, tclusterid);

    for (i = 0; i < nclusters; i++)
      MdrV0(*counts,i) = 0;
    for (i = 0; i < ndata; i++)
      MdrV0(*counts,(int) MdrV0(tclusterid,i)) += 1.0;

    /* Start the loop */
    while(1)
    {
      double previous = total;
      total = 0.0;

      if (counter % period == 0) /* Save the current cluster assignments */
      {
        for (i = 0; i < ndata; i++)
          saved[i] = MdrV0(tclusterid,i);

        if (period < INT_MAX / 2)
          period *= 2;
      }
      counter++;
      /* Find the center */
      if (imean)
      {
        getclustermeans  (rowdominant, nclusters, data, tclusterid, cdata);
      }
      else
      {
        getclustermedians(rowdominant, nclusters, data, tclusterid, cdata);
      }

      /* Calculate the distances */
      for (i=0; i<ndata; i++)
      {
        double distance;
        k = MdrV0(tclusterid,i);
        if (MdrV0(*counts,k)==1)
          continue;
        /* No reassignment if that would lead to an empty cluster */
        /* Treat the present cluster as a special case */
        distance = metric(1, data, *cdata, weight, i, k);
        for (j = 0; j < nclusters; j++)
        {
          double tdistance;
          if (j==k)
            continue;
          tdistance = metric(1, data, *cdata, weight, i, j);

          if (tdistance < distance)
          {
            distance = tdistance;
            MdrV0(*counts, (int) MdrV0(tclusterid,i)) -= 1.0;
            MdrV0(tclusterid,i) = j;
            MdrV0(*counts,j) += 1.0;
          }
        }
        total += distance;
      }

      if (total>=previous)
        break;

      for (i = 0; i < ndata; i++)
      {
        if (saved[i]!=MdrV0(tclusterid,i))
          break;
      }
      if (i==ndata)
        break; /* Identical solution found; break out of this loop */
    }

    if (total <= *error)
    {
      ifound = 1;
      *error = total;
      for (j=0; j<ndata; j++)
        MdrV0(*clusterid,j) = MdrV0(tclusterid,j);
    }

    if (npass<=1)
    {
      *error = total;
      break;
    }

    if (i==ndata)
      ifound++; /* break statement not encountered */

  }
  while (++ipass < npass);

  // clean-up
  if (saved)
    GC_FREE(saved);
  mdr_Destroy (tclusterid);

  // put cluster-id in human readable form
  for (i=0; i<ndata; i++)
    MdrV0(*clusterid,i) += 1.0;

//   dprintf("out\n");

  return ifound;
}


/*
 *    Purpose
 *    =======
 *
 *    The getclustermedoids routine calculates the cluster centroids, given to which
 *    cluster each element belongs. The centroid is defined as the element with the
 *    smallest sum of distances to the other elements.
 */
void getclustermedoids(int nclusters, MDR * distance, MDR * clusterid, int *centroids, double * errors)
{
  int i, j, k, ndata = MNR(clusterid);
  double d;

  for (j = 0; j < nclusters; j++)
    errors[j] = DBL_MAX;

  for (i = 0; i < ndata; i++)
  {
    d = 0.0;
    j = MdrV0(clusterid,i);
    for (k = 0; k < ndata; k++)
    {
      if (i==k || MdrV0(clusterid,k)!=j)
        continue;

      d += (i < k ? Mdr0(distance,k,i) : Mdr0(distance,i,k));

      if (d > errors[j])
        break;
    }

    if (d < errors[j])
    {
      errors[j] = d;
      centroids[j] = i;
    }
  }

  return;
}


/*
 * Purpose
 * =======
 *
 * The kmedoids routine performs k-medoids clustering on a given set of elements,
 * using the distance matrix and the number of clusters passed by the user.
 * Multiple passes are being made to find the optimal clustering solution, each
 * time starting from a different initial clustering.
 */
static int
kmedoids ( int nclusters, MDR * data, MDR * weight, int rowdominant, int npass, char dist,
           MDR ** cdata, MDR ** clusterid, double* error, MDR **distance, MDR **counts)
{
  int i, j, icluster;
  const int ndata = (rowdominant==1) ? MNR(data) : MNC(data);
  const int nattr = (rowdominant==1) ? MNC(data) : MNR(data);
  int ifound = 1;
  int ipass = 0;
  int find_distmat=0;

  //
  MDR *tclusterid = mdr_Float_BF(*clusterid);

  if (*distance)
  {
    if (MNR(*distance)!=ndata || MNC(*distance)!=ndata)
    {
      mdr_Destroy(*distance);
     *distance = mdr_Create(ndata,ndata);
      find_distmat = 1;
    }
  }
  else
  {
    *distance = mdr_Create(ndata,ndata);
    find_distmat = 1;
  }

  /* Set the metric function as indicated by dist */
  double (*metric)(int, MDR*, MDR*, MDR*, int, int) =  distance_metric(dist);
  if (find_distmat)
    for(i=1;i<ndata;i++)
      for(j=0;j<i;j++)
        Mdr0(*distance,i,j) = metric(1, data, data, weight, i, j);

  /* We save the clustering solution periodically and check if it reappears */
  int *    saved = GC_malloc(ndata*sizeof(int));
  int *centroids = GC_malloc(nclusters*sizeof(int));
  double *errors = GC_malloc(nclusters*sizeof(double));;

  *error = DBL_MAX;
  do
  {
    double total = DBL_MAX;
    int counter = 0;
    int period = 10;

    /* Perform the EM algorithm. First, randomly assign elements to clusters. */
    if (npass!=0)
      randomassign (nclusters, tclusterid);

    /* Start the loop */
    while(1)
    {
      double previous = total;
      total = 0.0;

      if (counter % period == 0) /* Save the current cluster assignments */
      {
        for (i = 0; i < ndata; i++)
          saved[i] = MdrV0(tclusterid,i);

        if (period < INT_MAX / 2)
          period *= 2;
      }
      counter++;

      /* Find the center */
      getclustermedoids(nclusters, *distance, tclusterid, centroids, errors);

      /* Calculate the distances */
      for (i=0; i<ndata; i++)
      {
        double dist = DBL_MAX;
        for (icluster=0; icluster<nclusters; icluster++)
        {
          double tdistance;
          j = centroids[icluster];
          if (i==j)
          {
            dist = 0.0;
            MdrV0(tclusterid,i) = icluster;
            break;
          }
          tdistance = (i > j) ? Mdr0(*distance,i,j) : Mdr0(*distance,j,i);
          if (tdistance < dist)
          {
            dist = tdistance;
            MdrV0(tclusterid,i) = icluster;
          }
        }
        total += dist;
      }

      if (total>=previous)
        break;

      /* total>=previous is FALSE on some machines even if total and previous
       * are bitwise identical. */
      for (i = 0; i < ndata; i++)
        if (saved[i]!=MdrV0(tclusterid,i))
          break;

      /* Identical solution found; break out of this loop */
      if (i==ndata)
        break;
    }

    for (i = 0; i < ndata; i++)
    {
      if (MdrV0(tclusterid,i)!=centroids[(int) MdrV0(tclusterid,i)])
      {
        if (total < *error)
        {
          ifound = 1;
          *error  = total;
          /* Replace by the centroid in each cluster. */
          for (j=0; j<ndata; j++)
          {
            MdrV0(*clusterid,j) = centroids[(int) MdrV0(tclusterid,j)];
          }
        }
        break;
      }
    }

    if (i==ndata)
      ifound++; /* break statement not encountered */

  } while (++ipass < npass);

  // count each class
  *counts = mdr_Create(nclusters,1);
  for (i=0; i<nclusters; i++)
  {
    for (j=0; j<nattr; j++)
      Mdr0(*cdata,i,j) = Mdr0(data,centroids[i],j);

    MdrV0(*counts,i) = 0;
    for (j=0; j<ndata;j++)
    {
      if (MdrV0(*clusterid,j) == centroids[i])
        MdrV0(*counts,i) += 1.0;
    }
  }

  // convert index of cluster centers to index in human readable form
  for (i=0; i<ndata; i++)
  {
    for (j=0; j<nclusters; j++)
    {
      if(MdrV0(*clusterid,i) == centroids[j])
      {
        MdrV0(*clusterid,i) = j+1;
        break;
      }
    }
  }

  // clean-up
  if (saved)
    GC_FREE(saved);
  if (errors)
    GC_FREE(errors);
  if (centroids)
    GC_FREE(centroids);
  if (find_distmat)
    mdr_Destroy(*distance);
  mdr_Destroy (tclusterid);


  return ifound;
}

//
// clustering according to nearest neighbour
//
#undef THIS_SOLVER
#define THIS_SOLVER "cluster.knn"
Ent * ent_cluster_nn (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *e3=0;

  MDR *data, *clusterid=0, *cdata=0, *counts=0, *distance=0;
  double error;
  int npass=2, nclusters, imean=1;

  FILE *rlab_stderr = (!RLAB_STDERR_DS) ? stderr : RLAB_STDERR_DS;

  char *m=0;

  MDR * weight=0;

  //
  // Load and Check arguments.
  //
  if (nargs != 2 && nargs != 3)
  {
    fprintf(rlab_stderr,
            THIS_SOLVER ": nearest neighbour clustering of data.\n");
    fprintf(rlab_stderr,
            THIS_SOLVER ": Format :\n");
    fprintf(rlab_stderr,
            THIS_SOLVER ":   c = "THIS_SOLVER"(data,n/,opts/)\n");
    fprintf(rlab_stderr,
            THIS_SOLVER ": where 'x' is the data set, 'n' is the number of centers.\n");
    fprintf(rlab_stderr,
            THIS_SOLVER ": and 'opts' is a list of options.\n");
    fprintf(rlab_stderr,
            THIS_SOLVER ": The result 'c' gives the cluster centers of 'x'.\n");
    rerror (THIS_SOLVER ": requires three arguments");
  }

  //
  // x
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) != MATRIX_DENSE_REAL)
    rerror (THIS_SOLVER ": " RLAB_ERROR_ARG1_MDR_MATRIX);
  data = class_matrix_real (e1);

  //
  // nclusters
  //
  e2 = bltin_get_ent (args[1]);
  if (ent_type (e2) != MATRIX_DENSE_REAL)
    rerror (THIS_SOLVER ": " RLAB_ERROR_ARG2_MDR_SCALAR);
  nclusters = (int) class_double (e2);

  //
  // options
  //
  if (nargs == 3)
  {
    e3 = bltin_get_ent (args[2]);
    if (ent_type (e3) != BTREE)
      rerror (THIS_SOLVER ": " RLAB_ERROR_ARG3_BTREE);

    // metric
    ListNode *node=0;
    node = btree_FindNode (ent_data (e3), RLAB_NAME_STAT_METRIC);
    if (node != 0)
    {
      m = class_char_pointer (var_ent (node));
      if (strlen(m)!=1)
        m=0;
    }
    // weight for the metric function
    node = btree_FindNode (ent_data (e3), RLAB_NAME_STAT_WEIGHT);
    if (node != 0)
    {
      weight = class_matrix_real (var_ent (node));
      if (SIZE(weight)!= MNC (data))
        weight = 0;
    }
    // number of attempts to find best clustering
    node = btree_FindNode (ent_data (e3), RLAB_NAME_GEN_MAXITER);
    if (node != 0)
    {
      npass = (int) class_double (var_ent (node));
      if (npass < 0)
        npass = 10;
    }
    // number of attempts to find best clustering
    node = btree_FindNode (ent_data (e3), RLAB_NAME_GEN_METHOD);
    if (node != 0)
    {
      char * s = class_char_pointer (var_ent (node));
      if (s[0] == 'm' || s[0]=='M')
        imean=0;  // median
      else if (s[0] == 'a' || s[0]=='A')
        imean=1;  // average
      else
        imean=-1; // medoid
    }
    // use provided distance matrix of the dataset
    node = btree_FindNode (ent_data (e3), RLAB_NAME_STAT_DISTMAT);
    if (node != 0)
    {
      distance = class_matrix_real (var_ent (node));
      imean=-1;
    }
  }

  // default metric is euclidean
  if (!m)
    m="e";

  clusterid = mdr_Create(MNR(data),1);
  mdr_Zero(clusterid);
  cdata     = mdr_Create(nclusters, MNC(data));
  mdr_Zero(cdata);

  if (imean == 0 || imean == 1)
  {
//     dprintf("before\n");
    kmeanmed(nclusters, data, weight, 1, npass, *m, &cdata, &clusterid, &error,
             &counts, imean);
//     dprintf("after\n");
  }
  else
    kmedoids(nclusters, data, weight, 1, npass, *m, &cdata, &clusterid, &error,
             &distance, &counts);

  // clean-up
  ent_Clean (e1);
  ent_Clean (e2);
  ent_Clean (e3);

  //
  // write result to output list
  //
  Btree *bw = btree_Create ();
  // value: cluster centres
  install (bw, RLAB_NAME_GEN_VALUE, ent_Assign_Rlab_MDR(cdata));
  // size: cluster sizes
  install (bw, RLAB_NAME_GEN_SIZE,  ent_Assign_Rlab_MDR(counts));
  // size: feature of data (cluster id datum belongs to)
  install (bw, RLAB_NAME_STAT_FEATURE,  ent_Assign_Rlab_MDR(clusterid));

  //
  return ent_Assign_Rlab_BTREE(bw);
}


//
// clustering according to isodata
// I went through this, and fixed their code. This is just "IN MEMORIAM"
// (C) M.K. 2014
/*
 * ISODATA.C : Implementation of the basic isodata algorithm.
 *
 * References : R. Duda and P. Hart, "Pattern Classification and Scene
 *              Analysis", John Wiley and Sons, 1973.
 *
 * Made by : W.F. Schmidt
 * SccsId  : %Z%$RCSfile: isodata.c,v $ V%Q% $Revision: 1.5 $ $Date: 1996/10/07 06:50:24 $
 */

/*
 * COMPUTE_ERROR_DATASET : Compute the stopping criterion for the basic
 * ISODATA clustering algorithm. The mean squared error between the samples
 * is calculated and returned.
 */
static double
    compute_error_dataset (MDR * dset1, MDR * dset2, char *m)
{
  // Set the metric function as indicated by dist
  double (*distance)(int, MDR*, MDR*, MDR*, int, int) =  distance_metric(*m);

  double rval = create_inf();
  int i1=0, i2=0;
  int nr1 = MNR(dset1), nr2=MNR(dset2);

  if (MNC(dset1)!=MNC(dset2))
    return rval;

  while ((i1 < nr1) && (i2 < nr2))
  {
    rval += distance(1, dset1, dset2, NULL, i1, i2);
    i1++;
    i2++;
  }

  return rval;
}

static void isodata_basic (int nclusters, MDR * dset, MDR ** clusterid, MDR ** clustersize,
                           MDR ** cdata, double *target_mse, int maxi, int * status,
                           char *m)
{
  int i, j, k, nattr, ndata;
  double mse;
  double (*distance)(int, MDR*, MDR*, MDR*, int, int) =  distance_metric(*m);

  *status = 0;

  nattr = MNC(dset);
  ndata = MNR(dset);
  *clusterid    = mdr_Create(ndata,1);
  *clustersize  = mdr_Create(nclusters, 1);

  // initial assignment: make sure all clusters have some points, rather then
  // having completely random assignments (some might not get assigned)
  mdr_Zero(*clustersize);
  for(j=0;j<ndata;j++)
  {
    MdrV0(*clusterid,j) = j%nclusters;
    MdrV0(*clustersize,j%nclusters) += 1;
  }

  // find cluster means for initial data clustering
  *cdata = mdr_Create(nclusters, nattr);
  getclustermeans(1, nclusters, dset, *clusterid, cdata);

 /*
  * Perform the clustering.
  */
  i=0;
  while (i<maxi)
  {
    i++;
    int have_updates = 0;
    for (j=0; j<ndata; j++)
    {
      int best_cluster_id = -1;
      double best_cluster_d  = 1e99;
      for (k=0; k<nclusters; k++)
      {
        double d = distance(1, dset, *cdata, NULL, j, k);
        if (d<best_cluster_d)
        {
          best_cluster_id = k;
          best_cluster_d  = d;
        }
      }
      if (MdrV0(*clusterid,j) != best_cluster_id)
      {
        MdrV0(*clustersize,(int) MdrV0(*clusterid,j)) -= 1;
        MdrV0(*clusterid,j)    = best_cluster_id;
        MdrV0(*clustersize,best_cluster_id)  += 1;
        have_updates = 1;
      }
    }

    if (have_updates)
    {
      MDR *cdata_old    = mdr_Copy(*cdata);
      getclustermeans(1, nclusters, dset, *clusterid, cdata);  // update cdata with new clusterid
      mse = compute_error_dataset(*cdata, cdata_old, m);
      mdr_Destroy(cdata_old);

      if (mse > *target_mse)
        continue;
    }

    break; // we went over entire dataset and did not update any of the clusters. time to go home

  } // while ((i<maxi) && (mse > *target_mse))

  if (mse > *target_mse)
    *status = 1; // No convergence in MAXTRIES cycles

  // update id by 1
  for (i=0; i<ndata; i++)
    MdrV0(*clusterid,i) += 1.0;

  *target_mse = mse;
  return;
}


#undef THIS_SOLVER
#define THIS_SOLVER "cluster.iso"
Ent *
    ent_cluster_isodata (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *e3=0;
  MDR *x=0;
  char *m = "e";

  double target_mse=0;
  int nclusters, status, maxi = 100;
  FILE *rlab_stderr = (!RLAB_STDERR_DS) ? stderr : RLAB_STDERR_DS;

  //
  // Load and Check arguments.
  //
  if (nargs != 2 && nargs != 3)
  {
    fprintf(rlab_stderr,
            THIS_SOLVER ": Perform clustering of the data set according to Duda and Hart.\n");
    fprintf(rlab_stderr,
            THIS_SOLVER ": Format :\n");
    fprintf(rlab_stderr,
            THIS_SOLVER ":   c = cluster." THIS_SOLVER "(data,k/,opts/)\n");
    fprintf(rlab_stderr,
            THIS_SOLVER ": where 'x' is the data set, 'k' is the number of centers,\n");
    fprintf(rlab_stderr,
            THIS_SOLVER ": 'd' is the threshold value for the mean of the cluster.\n");
    fprintf(rlab_stderr,
            THIS_SOLVER ": 'm' is the optional parameter specifying metrics 'e' for euclidean,\n");
    fprintf(rlab_stderr,
            THIS_SOLVER ": 'b' for city block, and so forth.\n");
    fprintf(rlab_stderr,
            THIS_SOLVER ": The result 'c' gives the cluster centers of 'x'.\n");
    rerror
        (THIS_SOLVER ": requires three arguments");
  }

  //
  // x
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) != MATRIX_DENSE_REAL)
    rerror (THIS_SOLVER ": " RLAB_ERROR_ARG1_MDR_MATRIX);
  x = class_matrix_real (e1);

  //
  // k : number of clusters
  //
  e2 = bltin_get_ent (args[1]);
  if (ent_type (e2) != MATRIX_DENSE_REAL)
    rerror (THIS_SOLVER ": " RLAB_ERROR_ARG2_MDR_INTEGER);
  nclusters = (int) class_double (e2);
  nclusters = MAX(nclusters,2);

  //
  // options
  //
  if (nargs == 3)
  {
    e3 = bltin_get_ent (args[2]);
    if (ent_type (e3) != BTREE)
      rerror (THIS_SOLVER ": " RLAB_ERROR_ARG3_BTREE);

    // metric
    ListNode *node=0;
    node = btree_FindNode (ent_data (e3), RLAB_NAME_STAT_METRIC);
    if (node != 0)
    {
      m = class_char_pointer (var_ent (node));
      if (strlen(m)!=1)
        m=0;
    }
    // number of attempts to find best clustering
    node = btree_FindNode (ent_data (e3), RLAB_NAME_GEN_MAXITER);
    if (node != 0)
    {
      maxi = (int) class_double (var_ent (node));
      if (maxi < 0)
        maxi = 100;
    }

    // number of attempts to find best clustering
    node = btree_FindNode (ent_data (e3), RLAB_NAME_GEN_TOL);
    if (node != 0)
    {
      target_mse = class_double (var_ent (node));
      if (target_mse < 0)
        target_mse = 0;
    }
  }

  // default metric is euclidean
  if (!m)
    m="e";

  //
  // d
  //
  e3 = bltin_get_ent (args[1]);
  if (ent_type (e3) != MATRIX_DENSE_REAL)
    rerror (THIS_SOLVER ": " RLAB_ERROR_ARG1_MDR_REALPOSITIVESCALAR);
  target_mse = class_double (e3);

  MDR *clusterid=0;
  MDR *clustersize=0;
  MDR *cdata=0;

  isodata_basic (nclusters, x, &clusterid, &clustersize, &cdata, &target_mse, maxi, &status, m);

  // clean-up
  ent_Clean (e1);
  ent_Clean (e2);
  ent_Clean (e3);

  //
  // write result to output list
  //
  Btree *bw = btree_Create ();
  // value: cluster centres
  install (bw, RLAB_NAME_GEN_VALUE, ent_Assign_Rlab_MDR(cdata));
  // size: cluster sizes
  install (bw, RLAB_NAME_GEN_SIZE,  ent_Assign_Rlab_MDR(clustersize));
  // size: feature of data (cluster id datum belongs to)
  install (bw, RLAB_NAME_STAT_FEATURE,  ent_Assign_Rlab_MDR(clusterid));

  //
  return ent_Assign_Rlab_BTREE(bw);
}


static void
extract_class_from_mdr (MDR *y, int *ncls, MDR **cls, MDR *yicls)
{
  int ndata = MNR(y)*MNC(y);
  int i, j;

  // sort indices rather than data
  size_t * idx = GC_malloc(ndata * sizeof(size_t));
  gsl_sort_index ((size_t *) idx, MDRPTR(y), (size_t) 1, (size_t) ndata);

  MDR * allcls = mdr_Create(1, ndata);

  // find the number of different entries in sorted y
  // and in the process replace y with classet
  *ncls = 1;
  MdrV0(yicls, idx[0]) = *ncls - 1;
  MdrV0(allcls,0) = MdrV0(y,idx[0]);
  for (i=1; i<ndata; i++)
  {
    if ( MdrV0(y,idx[i-1]) != MdrV0(y,idx[i]) )
    {
      MdrV0(allcls,*ncls) = MdrV0(y,idx[i]);
      (*ncls)++;
    }
    MdrV0(yicls, idx[i]) = *ncls - 1;
  }

  *cls = mdr_Create(1, *ncls);
  for (j=0; j<(*ncls); j++)
    MdrV0(*cls,j) = MdrV0(allcls,j);

  GC_FREE (idx);
  mdr_Destroy(allcls);

  return;
}

static int
string_compare(const void *s1, const void *s2)
{
  /* The actual arguments to this function are "pointers to
   *              pointers to char", but strcmp(3) arguments are "pointers
   *              to char", hence the following cast plus dereference */
  return strcmp(* (char * const *) s1, * (char * const *) s2);
}

static void
extract_class_from_mds (MDS *y, int *ncls, MDS **cls, MDR *yset)
{
  int ndata = MNR(y)*MNC(y);
  int i, j;
  MDS * allcls = mds_Create(1, ndata);

  // sort indices rather than data:
  //  we use heapsort, because this is only one that GSL provides that can be used for strings
  //  as we do not want to copy string matrices, then sort them, then toss them away.
  size_t * idx = GC_malloc(ndata * sizeof(size_t));
  gsl_heapsort_index ((size_t *) idx, MDSPTR(y), (size_t) ndata, sizeof(char*), (gsl_comparison_fn_t) string_compare);

  // find the number of different entries in sorted y
  // and in the process replace y with classet
  *ncls = 1;
  MdrV0(yset, idx[0]) = *ncls - 1;
  MdsV0(allcls,0) = cpstr(MdsV0(y,idx[0]));
  for (i=1; i<ndata; i++)
  {
    if ( strcmp(MdsV0(y,idx[i-1]), MdsV0(y,idx[i])) )
    {
      MdsV0(allcls,*ncls) = cpstr(MdsV0(y,idx[i]));
      (*ncls)++;
    }
    MdrV0(yset, idx[i]) = *ncls - 1;
  }

  // reassign pointers from one matrix to the other, and nullify them in the latter.
  *cls = mds_Create(1, *ncls);
  for (j=0; j<(*ncls); j++)
  {
    MdsV0(*cls,j) = MdsV0(allcls,j);  // repoint the pointers
    MdsV0(allcls,j) = 0;              // zero these ones so they are not destroyed with allcls
  }

  GC_FREE (idx);
  mds_Destroy(allcls);
}

//
// classify.parzen
//
#undef THIS_SOLVER
#define THIS_SOLVER "classify.parzen"
#include "gnunet_parzen.c"
Ent *
ent_classify_parzen (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *rent=0;
  MDR *x1=0, *x0=0, *icls=0;
  MDR *df0=0, *dcls=0, *df1=0;
  MDS *sf0=0, *scls=0, *sf1=0;

  ListNode * node;

  int i, j, ncls=0, maxit=1000;

  double init_s=1.1, s3=0.0;
  FILE *rlab_stderr = (!RLAB_STDERR_DS) ? stderr : RLAB_STDERR_DS;

  //
  // Load and Check arguments.
  //
  if (nargs != 2)
  {
    fprintf(rlab_stderr,
            THIS_SOLVER ": Classify the new data based on the training dataset.\n");
    fprintf(rlab_stderr,
            THIS_SOLVER ": Format :\n");
    fprintf(rlab_stderr,
            THIS_SOLVER ":   fx = classify.parzen(x,dataset)\n");
    fprintf(rlab_stderr,
            THIS_SOLVER ": where 'x' is the data to be classified based on the training set\n");
    fprintf(rlab_stderr,
            THIS_SOLVER ": dataset=<<data;feature>>.\n");
    fprintf(rlab_stderr,
            THIS_SOLVER ": The result are the suggested features of 'x'.\n");
    rerror
        (THIS_SOLVER ": requires two arguments");
  }

  //
  // x
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) != MATRIX_DENSE_REAL)
    rerror (THIS_SOLVER ": 'x' has to be a real matrix!");
  x1 = class_matrix_real (e1);

  //
  // training data
  //
  e2 = bltin_get_ent (args[1]);
  if (ent_type (e2) == BTREE)
  {
    // data
    node = btree_FindNode (ent_data (e2), RLAB_NAME_PATTERN_DATA);
    if (node != 0)
    {
      if(ent_type(var_ent (node))!=MATRIX_DENSE_REAL)
        rerror(THIS_SOLVER ": entry 'data' has to be a real matrix!");
      x0 = class_matrix_real (var_ent (node));
    }
    else
      rerror (THIS_SOLVER ": missing entry 'data' !");

    //
    // feature
    //
    node = btree_FindNode (ent_data (e2), RLAB_NAME_PATTERN_FEATURE);
    if (node != 0)
    {
      if(ent_type(var_ent (node))==MATRIX_DENSE_REAL)
      {
        df0 = class_matrix_real (var_ent (node));           // features
        icls = mdr_Create (1, MNR(df0)*MNC(df0));           // features converted to int 0:ncls
        extract_class_from_mdr (df0, &ncls, &dcls, icls);   // extract ncls and fill in icls
      }
      else if(ent_type(var_ent (node))==MATRIX_DENSE_STRING)
      {
        sf0 = ent_data (var_ent (node));
        icls = mdr_Create (1, MNR(sf0)*MNC(sf0));
        extract_class_from_mds (sf0, &ncls, &scls, icls);
      }
      else
        rerror(THIS_SOLVER ": entry 'feat' has to be a real or a string vector!");
    }
    else
      rerror (THIS_SOLVER ": missing entry 'feat' !");
  }

  if (MNC(x1) != MNC(x0))
    rerror(THIS_SOLVER ": 'data' and 'x' have to have same number of columns!");

  //
  // find the number of classes and their aposteriori probabilities
  //
  MDR * icount = mdi_Create (ncls,1);
  mdr_Zero (icount);
  for (i = 0; i < MNR(x0); i++)
    MdiV0(icount, (int) MdrV0(icls,i)) += 1;

  MDR *p = mdr_Create (ncls,1);
  j = 0;
  for (i=0; i< ncls; i++)
  {
    MdrV0(p,j) = (double) MdiV0(icount,j) / (double) MNR(x0);
    j++;
  }

  //
  // find best smoothing parameters from the LearningSet
  //
  MDR *s = mdr_Create (ncls,1);
  for (i = 0; i < ncls; i++)
  {
    if (parzen_dset_best_s (x0, icls, ncls, i, init_s, &s3, maxit))
      MdrV0(s,i) = s3;
  }

  // find closest classes

  MDR *y1 = parzen_dataset_class(x0, icls, x1, s, p, ncls);

  //
  // convert y1 to proper output
  //
  if (df0)
  {
    df1 = mdr_Create (MNR(y1) * MNC(y1),1);
    for (i=0; i < MNR(y1) * MNC(y1); i++)
      MdrV0(df1,i) = MdrV0(dcls, (int) MdrV0(y1,i));
    rent = ent_Assign_Rlab_MDR(df1);
  }
  else if (sf0)
  {
    sf1 = mds_Create (MNR(y1) * MNC(y1),1);
    for (i = 0; i < MNR(y1) * MNC(y1); i++)
      MdsV0(sf1,i) = cpstr( MdsV0(scls, mdiV0(y1,i)) );
    rent = ent_Assign_Rlab_MDS(sf1);
  }

  // clean-up
  ent_Clean (e1);
  ent_Clean (e2);

  if (dcls)
    mdr_Destroy (dcls);
  if (scls)
    mds_Destroy (scls);
  mdr_Destroy (s);
  mdr_Destroy (p);
  mdr_Destroy (icount);
  mdr_Destroy (icls);
  mdr_Destroy (y1);

  return rent;
}

//
// classify.knn
//
#undef THIS_SOLVER
#define THIS_SOLVER "classify.knn"
Ent *
ent_classify_knn (int nargs, Datum args[])
{

  Ent *e1=0, *e2=0, *e3=0, *rent;
  MDR *x1=0, *x0=0, *icls=0, *y1=0, *k=0;
  MDR *df0 = 0, *dcls = 0, *df1 = 0;
  MDS *sf0 = 0, *scls = 0, *sf1 = 0;

  ListNode * node;

  int i, j, m, ncls=0, nk, rowdominant=1, ntrain, nsamp;

  char *mx = "e";

  // Set the metric function as indicated by dist
  double (*distance)(int, MDR*, MDR*, MDR*, int, int) =  distance_metric(*mx);

  FILE *rlab_stderr = (!RLAB_STDERR_DS) ? stderr : RLAB_STDERR_DS;


  //
  // Load and Check arguments.
  //
  if (nargs != 3 && nargs != 4)
  {
    fprintf(rlab_stderr,
            THIS_SOLVER ": Classify the new data based on the training dataset.\n");
    fprintf(rlab_stderr,
            THIS_SOLVER ": Format :\n");
    fprintf(rlab_stderr,
            THIS_SOLVER ":   fx = "THIS_SOLVER"(x,dataset,k)\n");
    fprintf(rlab_stderr,
            THIS_SOLVER ": where 'x' is the data to be classified based on the training set\n");
    fprintf(rlab_stderr,
            THIS_SOLVER ": dataset=<<"RLAB_NAME_PATTERN_DATA";"RLAB_NAME_PATTERN_FEATURE">>  and the number of nearest-neighbors 'k'.\n");
    fprintf(rlab_stderr,
            THIS_SOLVER ": The result are the suggested features of 'x'.\n");
    rerror
        (THIS_SOLVER ": requires at least three arguments");
  }

  //
  // x
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) != MATRIX_DENSE_REAL)
    rerror (THIS_SOLVER ": " RLAB_ERROR_ARG1_MDR_MATRIX);
  x1 = class_matrix_real (e1);

  // training data
  //
  e2 = bltin_get_ent (args[1]);
  if (ent_type (e2) != BTREE)
    rerror (THIS_SOLVER ": " RLAB_ERROR_ARG2_BTREE);

  // data
  node = btree_FindNode (ent_data (e2), RLAB_NAME_PATTERN_DATA);
  if (!node)
    rerror (THIS_SOLVER ": missing entry '"RLAB_NAME_PATTERN_DATA"' !");
  if(ent_type(var_ent (node))!=MATRIX_DENSE_REAL)
    rerror(THIS_SOLVER ": entry '"RLAB_NAME_PATTERN_DATA"' has to be a real matrix!");
  x0 = class_matrix_real (var_ent (node));

  //
  // feature
  //
  node = btree_FindNode (ent_data (e2), RLAB_NAME_PATTERN_FEATURE);
  if (!node)
    rerror (THIS_SOLVER ": missing entry '"RLAB_NAME_PATTERN_FEATURE"' !");

  if(ent_type(var_ent (node))==MATRIX_DENSE_REAL)
  {
    df0 = class_matrix_real (var_ent (node));           // features
    icls = mdr_Create (1, MNR(df0)*MNC(df0));           // features converted to int 0:ncls
    extract_class_from_mdr (df0, &ncls, &dcls, icls);   // extract ncls and fill in icls
  }
  else if(ent_type(var_ent (node))==MATRIX_DENSE_STRING)
  {
    sf0 = ent_data (var_ent (node));
    icls = mdr_Create (1, MNR(sf0)*MNC(sf0));
    extract_class_from_mds (sf0, &ncls, &scls, icls);
  }
  else
    rerror(THIS_SOLVER ": entry '"RLAB_NAME_PATTERN_FEATURE"' has to be real or string vector!");

  if (MNR(x0)!=SIZE(icls))
    rerror(THIS_SOLVER  ": entry '"RLAB_NAME_PATTERN_FEATURE"' has to match '"RLAB_NAME_PATTERN_DATA"' !");

  //
  // k
  //
  e3 = bltin_get_ent (args[2]);
  if (ent_type (e3) != MATRIX_DENSE_REAL)
    rerror (THIS_SOLVER ": " RLAB_ERROR_ARG3_MDR_INTEGER);
  k  = class_matrix_real(e3);
  nk = MNR (k) * MNC (k);

  if (MNC(x1) != MNC(x0))
    rerror(THIS_SOLVER ": ARG1 and ARG2 list entry '"RLAB_NAME_PATTERN_DATA"' must have same dimension!");

  // rows represent the data in x and y matrices
  ntrain = MNR(x0);
  nsamp  = MNR(x1);

  MDR *ds = mdr_Create(1,ntrain);

  int maxk = (int) rlab_dmax_mdr(k, 1);
  size_t *bk_idx = (size_t *) GC_malloc(maxk * sizeof(size_t)); // assume k is sorted
  int *cl_cnt    = (int *) GC_malloc(ncls * sizeof(int));

  y1 = mdr_Create(nsamp, nk);

  // go over each sample
  for (i=0; i<nsamp; i++)
  {
    // for each sample in x1[i;] construct distance vector in ds[]
    for (j=0; j<ntrain; j++)
      MdrV0(ds,j) = distance(rowdominant, x1, x0, NULL, i, j);

    // go over each 'k'
    for (m=0; m<nk; m++)
    {
      // now pick the first k[m]-smallest values from ds (that is their positions) - no hot water here
      gsl_sort_smallest_index ( (size_t*) bk_idx, (size_t) mdiV0(k,m), MDRPTR(ds), (size_t) 1, (size_t) ntrain);

      // array bk contains the smallest indices of k[m] smallest entries in ds, which are also
      // the rows of matrix dset
      for (j=0; j<ncls; j++)
        cl_cnt[j]=0;
      // count classes of the k[m] closest points
      for (j=0; j < mdiV0(k,m); j++)
        cl_cnt[ mdiV0(icls, (int) bk_idx[j])  ]++;
      // figure the class that has majority
      int best_class = 0;
      int best_count = cl_cnt[0];
      for (j=1; j<ncls; j++)
      {
        if (best_count < cl_cnt[j])
        {
          best_class = j;
          best_count = cl_cnt[j];
        }
      }

      Mdr0(y1,i,m) = best_class;
    }
  }

  GC_FREE (bk_idx);
  GC_FREE (cl_cnt);
  mdr_Destroy (ds);

  rent = ent_Create ();
  //
  // convert y1 to proper output
  //
  if (df0)
  {
    //
    df1 = mdr_Create (MNR(y1),MNC(y1));
    for (i = 0; i < MNR(y1); i++)
    {
      for (j = 0; j < MNC(y1); j++)
      {
        //fprintf(rlab_stderr, "y1  = %i\n", (int) Mdr0(y1,i,j));
        Mdr0(df1,i,j) = MdrV0(dcls, (int) Mdr0(y1,i,j));
      }
    }
    //
    mdr_Destroy (dcls);

    //
    ent_data (rent) = df1;
    ent_type (rent) = MATRIX_DENSE_REAL;
  }
  else if (sf0)
  {
    //
    sf1 = mds_Create (MNR(y1),MNC(y1));
    for (i = 0; i < MNR(y1); i++)
      for (j = 0; j < MNC(y1); j++)
        Mds0(sf1,i,j) = cpstr( MdsV0(scls, (int) Mdr0(y1,i,j)) );
    //
    mds_Destroy (scls);

    //
    ent_data (rent) = sf1;
    ent_type (rent) = MATRIX_DENSE_STRING;
  }

  // clean-up
  ent_Clean (e1);
  ent_Clean (e2);
  ent_Clean (e3);

  mdr_Destroy (icls);
  mdr_Destroy (y1);

  return rent;
}

//
// classify.fisher
//
#undef THIS_SOLVER
#define THIS_SOLVER "classify.fisher"
#include "gnunet_fisher.c"
Ent * ent_sprannlib_classify_fisher (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *e3=0, *rent=0;
  MDR *x1=0, *x0=0, *icls=0, *y1=0;
  MDR *df0=0, *dcls=0, *df1=0;
  MDS *sf0=0, *scls=0, *sf1=0;

  ListNode * node;

  int i, m, ncls=0, ic0=0, ic1=1;
  double c;
  FILE *rlab_stderr = (!RLAB_STDERR_DS) ? stderr : RLAB_STDERR_DS;

  //
  // Load and Check arguments.
  //
  if (nargs != 2 && nargs != 3)
  {
    fprintf(rlab_stderr,
            THIS_SOLVER ": Classify the new data based on the training dataset.\n");
    fprintf(rlab_stderr,
            THIS_SOLVER ": Format :\n");
    fprintf(rlab_stderr,
            THIS_SOLVER ":   fx = classify.fisher(x,dataset/,cli/)\n");
    fprintf(rlab_stderr,
            THIS_SOLVER ": where 'x' is the data to be classified based on the training set\n");
    fprintf(rlab_stderr,
            THIS_SOLVER ": dataset=<<"RLAB_NAME_PATTERN_DATA";"RLAB_NAME_PATTERN_FEATURE">>, and cli=[classA,classB], the classes with\n");
    fprintf(rlab_stderr,
            THIS_SOLVER ": respect to which Fisher linear discriminant should be used.\n");
    fprintf(rlab_stderr,
            THIS_SOLVER ": The result are the suggested features of 'x'.\n");
    rerror
        (THIS_SOLVER ": requires at least three arguments");
  }

  //
  // x
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) != MATRIX_DENSE_REAL)
    rerror (THIS_SOLVER ": " RLAB_ERROR_ARG1_MDR_MATRIX);
  x1 = class_matrix_real (e1);

  // training data
  //
  e2 = bltin_get_ent (args[1]);
  if (ent_type (e2) != BTREE)
    rerror (THIS_SOLVER ": " RLAB_ERROR_ARG2_BTREE);

  // data
  node = btree_FindNode (ent_data (e2), RLAB_NAME_PATTERN_DATA);
  if (!node)
    rerror (THIS_SOLVER ": missing entry '"RLAB_NAME_PATTERN_DATA"' !");
  if(ent_type(var_ent (node))!=MATRIX_DENSE_REAL)
    rerror(THIS_SOLVER ": entry '"RLAB_NAME_PATTERN_DATA"' has to be a real matrix!");
  x0 = class_matrix_real (var_ent (node));

  //
  // feature
  //
  node = btree_FindNode (ent_data (e2), RLAB_NAME_PATTERN_FEATURE);
  if (!node)
    rerror (THIS_SOLVER ": missing entry '"RLAB_NAME_PATTERN_FEATURE"' !");

  if(ent_type(var_ent (node))==MATRIX_DENSE_REAL)
  {
    df0 = class_matrix_real (var_ent (node));           // features
    icls = mdr_Create (1, SIZE(df0));           // features converted to int 0:(ncls-1)
    extract_class_from_mdr (df0, &ncls, &dcls, icls);   // extract ncls and fill in icls
  }
  else if(ent_type(var_ent (node))==MATRIX_DENSE_STRING)
  {
    sf0 = ent_data (var_ent (node));
    icls = mdr_Create (1, SIZE(sf0));
    extract_class_from_mds (sf0, &ncls, &scls, icls);
  }
  else
    rerror(THIS_SOLVER ": entry '"RLAB_NAME_PATTERN_FEATURE"' has to be real or string vector!");

  if (MNC(x1) != MNC(x0))
    rerror(THIS_SOLVER ": ARG1 and ARG2 list entry '"RLAB_NAME_PATTERN_DATA"' must have same dimension!");

  //
  // if number of classes is larger than 2 need classes 1 and 2
  //
  if (ncls > 2)
  {
    if (nargs >= 3)
    {
      //
      // get [classA, classB] from the user
      //
      e3 = bltin_get_ent (args[2]);
      if (ent_type (e3) == MATRIX_DENSE_REAL && !dcls)
      {
        //
        // convert classes to integers
        //
        MDR * x3 = ent_data(e3);
        for (i = 0; i < ncls; i++)
        {
          if (MdrV0(dcls,i) == MdrV0(x3,0))
          {
            ic0 = i;
            break;
          }
        }
        for (i = 0; i < ncls; i++)
        {
          if (MdrV0(dcls,i) == MdrV0(x3,1))
          {
            ic1 = i;
            break;
          }
        }
      }
      else if (ent_type (e3) == MATRIX_DENSE_STRING && !scls)
      {
        //
        // convert classes to integers
        //
        MDS * s3 = ent_data(e3);
        for (i = 0; i < ncls; i++)
        {  
          if (!strcmp(MdsV0(scls,i),MdsV0(s3,0)))
          {
            ic0 = i;
            break;
          }
        }
        for (i = 0; i < ncls; i++)
        {
          if (!strcmp(MdsV0(scls,i),MdsV0(s3,1)))
          {
            ic1 = i;
            break;
          }
        }
      }
      else
        rerror (THIS_SOLVER ": [featA,featB] does not match entry 'feat'");
    }
    else
      rerror (THIS_SOLVER ": missing [featA,featB] - not a two-class problem!");
  }
  if (ic0 == ic1)
    rerror (THIS_SOLVER ": [featA,featB] cannot be the same!");
  if (ic0 == ncls || ic1 == ncls)
    rerror (THIS_SOLVER ": featA or featB out of range!");

  if (MNC(x1) != MNC(x0))
    rerror(THIS_SOLVER ": 'data' and 'x' have to have same number of columns!");


  //
  // find the number of classes and their a posteriori probabilities
  //
  MDR * icount = mdi_Create (ncls,1);
  mdr_Zero (icount);
  for (i = 0; i < MNR(x0); i++)
    MdiV0(icount, (int) MdrV0(icls,i)) += 1;
  c = mdrV0(icount,ic1) / mdrV0(icount,ic0);

  //
  // determine fisher discriminant: w, b
  //
  double  b;
  MDR *w = mdr_Create(MNC(x1),1);

  fisher_dataset (x0, icls, ic0, ic1, c, w, &b);

  y1 = mdr_Create(MNR(x1), 1);

  //
  // classify with fisher's  w*x-b
  //
  for (i = 0; i < MNR(x1); i++)
  {
    MdrV0(y1,i) = 0.0;

    for (m=0; m < MNC(x0); m++)
      MdrV0(y1,i) += MdrV0(w,m) * mdr0(x1,i,m);

    MdrV0(y1,i) = MdrV0(y1,i) + b < 0 ? ic1 : ic0 ;
  }

  mdr_Destroy (w);
  mdr_Destroy (icount);

  rent = ent_Create ();
  //
  // convert y1 to proper output
  //
  if (df0)
  {
    //
    df1 = mdr_Create (MNR(y1),1);
    for (i = 0; i < MNR(y1); i++)
      MdrV0(df1,i) = MdrV0(dcls, (int) MdrV0(y1,i));

    //
    mdr_Destroy (dcls);

    //
    ent_data (rent) = df1;
    ent_type (rent) = MATRIX_DENSE_REAL;
  }
  else if (sf0)
  {
    //
    sf1 = mds_Create (MNR(y1),1);
    for (i = 0; i < MNR(y1); i++)
      MdsV0(sf1,i) = cpstr( MdsV0(scls, (int) MdrV0(y1,i)) );

    //
    mds_Destroy (scls);

    //
    ent_data (rent) = sf1;
    ent_type (rent) = MATRIX_DENSE_STRING;
  }

  // clean-up
  ent_Clean (e1);
  ent_Clean (e2);
  ent_Clean (e3);

  mdr_Destroy (icls);
  mdr_Destroy (y1);

  return rent;
}

//
// classify.fisher
//
#undef THIS_SOLVER
#define THIS_SOLVER "classify.fisherq"
#include "gnunet_fisherq.c"
Ent *
    ent_sprannlib_classify_fisherq (int nargs, Datum args[])
{

  Ent *e1=0, *e2=0, *e3=0, *rent=0;
  MDR *x1=0, *x0=0, *icls=0, *y1=0;
  MDR *df0=0, *dcls=0, *df1=0;
  MDS *sf0=0, *scls=0, *sf1=0;

  ListNode * node;

  int i, j, k, ncls=0, ic0=0, ic1=1;
  double c;
  FILE *rlab_stderr = (!RLAB_STDERR_DS) ? stderr : RLAB_STDERR_DS;

  //
  // Load and Check arguments.
  //
  if (nargs != 2 && nargs != 3)
  {
    fprintf(rlab_stderr,
            THIS_SOLVER ": Classify the new data based on the training dataset.\n");
    fprintf(rlab_stderr,
            THIS_SOLVER ": Format :\n");
    fprintf(rlab_stderr,
            THIS_SOLVER ":   fx = classify.fisher(x,dataset/,cli/)\n");
    fprintf(rlab_stderr,
            THIS_SOLVER ": where 'x' is the data to be classified based on the training set\n");
    fprintf(rlab_stderr,
            THIS_SOLVER ": dataset=<<"RLAB_NAME_PATTERN_DATA";"RLAB_NAME_PATTERN_FEATURE">>, and cli=[classA,classB], the classes with\n");
    fprintf(rlab_stderr,
            THIS_SOLVER ": respect to which Fisher linear discriminant should be used.\n");
    fprintf(rlab_stderr,
            THIS_SOLVER ": The result are the suggested features of 'x'.\n");
    rerror
        (THIS_SOLVER ": requires at least three arguments");
  }

  //
  // x
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) != MATRIX_DENSE_REAL)
    rerror (THIS_SOLVER ": " RLAB_ERROR_ARG1_MDR_MATRIX);
  x1 = class_matrix_real (e1);

  // training data
  //
  e2 = bltin_get_ent (args[1]);
  if (ent_type (e2) != BTREE)
    rerror (THIS_SOLVER ": " RLAB_ERROR_ARG2_BTREE);

  // data
  node = btree_FindNode (ent_data (e2), RLAB_NAME_PATTERN_DATA);
  if (!node)
    rerror (THIS_SOLVER ": missing entry '"RLAB_NAME_PATTERN_DATA"' !");
  if(ent_type(var_ent (node))!=MATRIX_DENSE_REAL)
    rerror(THIS_SOLVER ": entry '"RLAB_NAME_PATTERN_DATA"' has to be a real matrix!");
  x0 = class_matrix_real (var_ent (node));

  //
  // feature
  //
  node = btree_FindNode (ent_data (e2), RLAB_NAME_PATTERN_FEATURE);
  if (!node)
    rerror (THIS_SOLVER ": missing entry '"RLAB_NAME_PATTERN_FEATURE"' !");

  if(ent_type(var_ent (node))==MATRIX_DENSE_REAL)
  {
    df0 = class_matrix_real (var_ent (node));           // features
    icls = mdr_Create (1, MNR(df0)*MNC(df0));           // features converted to int 0:(ncls-1)
    extract_class_from_mdr (df0, &ncls, &dcls, icls);   // extract ncls and fill in icls
  }
  else if(ent_type(var_ent (node))==MATRIX_DENSE_STRING)
  {
    sf0 = ent_data (var_ent (node));
    icls = mdr_Create (1, MNR(sf0)*MNC(sf0));
    extract_class_from_mds (sf0, &ncls, &scls, icls);
  }
  else
    rerror(THIS_SOLVER ": entry '"RLAB_NAME_PATTERN_FEATURE"' has to be real or string vector!");

  if (MNC(x1) != MNC(x0))
    rerror(THIS_SOLVER ": ARG1 and ARG2 list entry '"RLAB_NAME_PATTERN_DATA"' must have same dimension!");

  //
  // if number of classes is larger than 2 need classes 1 and 2
  //
  if (ncls > 2)
  {
    if (nargs >= 3)
    {
      //
      // get [classA, classB] from the user
      //
      e3 = bltin_get_ent (args[2]);
      if (ent_type (e3) == MATRIX_DENSE_REAL && !dcls)
      {
        //
        // convert classes to integers
        //
        MDR * x3 = ent_data(e3);
        for (i = 0; i < ncls; i++)
          if (MdrV0(dcls,i) == MdrV0(x3,0))
        {
          ic0 = i;
          break;
        }
        for (i = 0; i < ncls; i++)
          if (MdrV0(dcls,i) == MdrV0(x3,1))
        {
          ic1 = i;
          break;
        }
      }
      else if (ent_type (e3) == MATRIX_DENSE_STRING && !scls)
      {
        //
        // convert classes to integers
        //
        MDS * s3 = ent_data(e3);
        for (i = 0; i < ncls; i++)
          if (!strcmp(MdsV0(scls,i),MdsV0(s3,0)))
        {
          ic0 = i;
          break;
        }
        for (i = 0; i < ncls; i++)
          if (!strcmp(MdsV0(scls,i),MdsV0(s3,1)))
        {
          ic1 = i;
          break;
        }
      }
      else
        rerror (THIS_SOLVER ": [featA,featB] does not match entry 'feat'");
    }
    else
      rerror (THIS_SOLVER ": missing [featA,featB] - not a two-class problem!");
  }
  if (ic0 == ic1)
    rerror (THIS_SOLVER ": [featA,featB] cannot be the same!");
  if (ic0 == ncls || ic1 == ncls)
    rerror (THIS_SOLVER ": featA or featB out of range!");

  if (MNC(x1) != MNC(x0))
    rerror(THIS_SOLVER ": 'data' and 'x' have to have same number of columns!");


  //
  // find the number of classes and their a posteriori probabilities
  //
  MDR * icount = mdi_Create (ncls,1);
  mdr_Zero (icount);
  for (i = 0; i < MNR(x0); i++)
    MdiV0(icount, (int) MdrV0(icls,i)) += 1;
  c = mdrV0(icount,ic1) / mdrV0(icount,ic0);

  //
  // determine fisher discriminant: w, b
  //
  double  b;
  MDR *w1 = mdr_Create(MNC(x1),1);
  MDR *w2 = mdr_Create(ncls,MNC(x1));

  fisherq_dataset (x0, icls, ic0, ic1, c, w2, w1, &b);

  y1 = mdr_Create(MNR(x1), 1);

  //
  // classify with fisher's  w*x-b
  //
  for (i = 0; i < MNR(x1); i++)
  {
    MdrV0(y1,i) = 0.0;
    for (j = 0; j < MNC(x0); j++)
    {
      for (k = 0; k < MNC(x0); k++)
        MdrV0(y1,i) += mdr0(x1,i,j) * Mdr0(w2,j,k) * mdr0(x1,i,k);

      MdrV0(y1,i) += MdrV0(w1,j) * mdr0(x1,i,j);
    }
    MdrV0(y1,i) = MdrV0(y1,i) + b < 0 ? ic1 : ic0 ;
  }

  mdr_Destroy (w1);
  mdr_Destroy (w2);
  mdr_Destroy (icount);

  rent = ent_Create ();
  //
  // convert y1 to proper output
  //
  if (df0)
  {
    //
    df1 = mdr_Create (MNR(y1),1);
    for (i = 0; i < MNR(y1); i++)
      MdrV0(df1,i) = MdrV0(dcls, (int) MdrV0(y1,i));

    //
    mdr_Destroy (dcls);

    //
    ent_data (rent) = df1;
    ent_type (rent) = MATRIX_DENSE_REAL;
  }
  else if (sf0)
  {
    //
    sf1 = mds_Create (MNR(y1),1);
    for (i = 0; i < MNR(y1); i++)
      MdsV0(sf1,i) = cpstr( MdsV0(scls, (int) MdrV0(y1,i)) );

    //
    mds_Destroy (scls);

    //
    ent_data (rent) = sf1;
    ent_type (rent) = MATRIX_DENSE_STRING;
  }

  // clean-up
  ent_Clean (e1);
  ent_Clean (e2);
  ent_Clean (e3);

  mdr_Destroy (icls);
  mdr_Destroy (y1);

  return rent;
}

//
// sammon mapping (data reduction)
//
#undef THIS_SOLVER
#define THIS_SOLVER "sammon"
#include "gnunet_sammon.c"
Ent *
    ent_sprannlib_sammon (int nargs, Datum args[])
{

  Ent *e1=0, *e2=0, *e3=0;
  MDR *x0=0, *map1=0, *map=0, *set_distance=0;

  ListNode * node;

  int idummy, map_dim=0, num_samples, max_steps=1000, min_max_update=1e-9;
  int imethod=0, idx=0;

  double lrate = 1e-3, momentum = 0.1, ddummy, constant;

  FILE *rlab_stderr = (!RLAB_STDERR_DS) ? stderr : RLAB_STDERR_DS;

  //
  // Load and Check arguments.
  //
  if (nargs < 2 || nargs > 3)
  {
    fprintf(rlab_stderr,
            THIS_SOLVER ": Data reduction using sammon mapping.\n");
    fprintf(rlab_stderr,
            THIS_SOLVER ": Format :\n");
    fprintf(rlab_stderr,
            THIS_SOLVER ": (1) rdata = sammon(<<data;feat>>,dim,options),\n");
    fprintf(rlab_stderr,
            THIS_SOLVER ": (2) rdata = sammon(data,dim,options),\n");
    fprintf(rlab_stderr,
            THIS_SOLVER ": where 'rdata' is a reduced dataset of a dimension 'dim'.\n");
    fprintf(rlab_stderr,
            THIS_SOLVER ": The list options=<<>>.\n");
    rerror
        (THIS_SOLVER ": requires at least three arguments");
  }

  //
  // input data
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) == MATRIX_DENSE_REAL)
    x0 = class_matrix_real (e1);
  else if (ent_type (e1) == BTREE)
  {
    // data
    node = btree_FindNode (ent_data (e1), RLAB_NAME_PATTERN_DATA);
    if (!node)
      rerror (THIS_SOLVER ": missing entry '"RLAB_NAME_PATTERN_DATA"' !");
    if(ent_type(var_ent (node))!=MATRIX_DENSE_REAL)
      rerror(THIS_SOLVER ": entry '"RLAB_NAME_PATTERN_DATA"' has to be a real matrix!");
    x0 = class_matrix_real (var_ent (node));
  }
  else
    rerror (THIS_SOLVER ": Incorrect first entry - a list <<data;feat>> or a matrix");

  //
  // dim
  //
  e2 = bltin_get_ent (args[1]);
  if (ent_type (e2) == MATRIX_DENSE_REAL)
    map_dim = (int) class_double (e2);
  else
    rerror (THIS_SOLVER ": Incorrect argument 'dim'");

  if (map_dim == 0 || map_dim >= MNC(x0))
    rerror (THIS_SOLVER ": Incorrect argument 'dim' - nonexistent or greater than dimension of 'data'!");

  //
  // options
  //
  if (nargs > 2)
  {
    e3 = bltin_get_ent (args[2]);
    if (ent_type (e3) == BTREE)
    {
      // imethod
      node = btree_FindNode (ent_data (e3), RLAB_NAME_GEN_IMETHOD);
      if (node != 0)
      {
        idummy = (int) class_double ( var_ent (node) );
        if (idummy == 0 || idummy == 1)
          imethod = idummy;
      }
      // lrate
      node = btree_FindNode (ent_data (e3), RLAB_NAME_PATTERN_BACKPROP_LEARNING_RATE);
      if (node != 0)
      {
        ddummy = class_double ( var_ent (node) );
        if (ddummy > 0 && ddummy < 1)
          lrate = ddummy;
      }
      // momentum
      node = btree_FindNode (ent_data (e3), RLAB_NAME_PATTERN_BACKPROP_MOMENTUM);
      if (node != 0)
      {
        ddummy = class_double ( var_ent (node) );
        if (ddummy > 0 && ddummy < 1)
          momentum = ddummy;
      }
      // max_steps
      node = btree_FindNode (ent_data (e3), RLAB_NAME_GEN_MAXITER);
      if (node != 0)
      {
        idummy = (int) class_double ( var_ent (node) );
        if (idummy > 0)
          max_steps = idummy;
      }
      // map1
      node = btree_FindNode (ent_data (e3), RLAB_NAME_PATTERN_SAMMAN_MAP_GUESS);
      if (node != 0)
      {
        if (ent_type(var_ent(node)) == MATRIX_DENSE_REAL)
          map1 = class_matrix_real ( var_ent (node) );
      }
    }
  }

  //
  // sammon: init
  //
  num_samples = MNR(x0);
  set_distance = mdr_Create(num_samples,num_samples);
  if (map1)
  {
    if (MNR(map1) == num_samples && MNC(map1) == map_dim)
      map = mdr_Float_BF (map1);
  }
  if (!map)
  {
    map = mdr_Create (num_samples, map_dim);
    gsl_random_number_generator (&num_samples, &map_dim, MDRPTR(map), &idx );
  }

  constant = sammon_set_distance (x0, &set_distance);

  // imethod = 0: SAMMON_GD;
  // imethod = 1: SAMMON_PN;
  sammon_solve (&map, set_distance, constant, imethod, lrate, momentum, max_steps, min_max_update);

  //
  // clean-up
  //
  mdr_Destroy(set_distance);
  ent_Clean (e1);
  ent_Clean (e2);
  ent_Clean (e3);

  return ent_Assign_Rlab_MDR(map);
}

//
// setup generic neural network
//
#include <time.h>

#define RLAB_FFNN_INIT_RAND     0x0001
#define RLAB_SAMANN_RANDOM_PAIR 0x0100
#define RLAB_SAMANN_ALL_PAIRS   0x0200
#include "gnunet_ffnn.c"

//
// evaluate neural network
//
#undef THIS_SOLVER
#define THIS_SOLVER "ffnneval"
Ent * ent_ffnn_eval (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0;
  MDR *y_out=0, *x=0, *NetVect=0;
  MDR **NetPtr_wgt=NULL, **NetPtr_ifwgt=NULL, **NetPtr_theta=NULL, **NetPtr_iftheta=NULL;
  MDS **NetPtr_act=NULL, **NetPtr_transf=NULL;


  ListNode * node=0, * node_wgt=0, * node_ifwgt=0, * node_theta=0;
  ListNode * node_iftheta=0, * node_act=0, * node_trans=0;

  int i, L=0;
  FILE *rlab_stderr = (!RLAB_STDERR_DS) ? stderr : RLAB_STDERR_DS;

  //
  // Load and Check arguments.
  //
  if (nargs != 2)
  {
    fprintf(rlab_stderr,
            THIS_SOLVER ": Evaluate feed forward neural network.\n"
            THIS_SOLVER ": Format :\n"
            THIS_SOLVER ":   y = "THIS_SOLVER"(x,nn)\n"
            THIS_SOLVER ": See rlabplus manual for more information.\n"
           );
    rerror (THIS_SOLVER ": requires two arguments");
  }

  //
  // x
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) != MATRIX_DENSE_REAL)
    rerror (THIS_SOLVER ": 'X' has to be a real matrix!");
  x = class_matrix_real (e1);

  //
  // second argument is the nn structured list or the outline of the network
  //
  e2 = bltin_get_ent (args[1]);
  if (ent_type (e2) != BTREE)
    rerror(THIS_SOLVER ": " RLAB_ERROR_ARG2_BTREE);

  //
  // create a feed forward neural network using information from
  // 'nodes'.
  //
  node = btree_FindNode (ent_data (e2), RLAB_NAME_FFNN_NODES);
  if (!node)
    rerror (THIS_SOLVER ": entry '"RLAB_NAME_FFNN_NODES"' has to be real vector!");
  if(ent_type(var_ent (node))!=MATRIX_DENSE_REAL)
    rerror(THIS_SOLVER ": entry 'nodes' has to be real vector!");
  NetVect = ent_data (var_ent (node));

  //
  // fill the network using data from the lists 'wgt', 'fix_wgt', 'bias', 'fix_bias',
  // 'act' and 'trans'
  //
  node_wgt = btree_FindNode (ent_data (e2), RLAB_NAME_FFNN_WEIGHTS);
  if (!node_wgt)
    goto _exit_ent_ffnn_eval;
  if(ent_type(var_ent (node_wgt)) != BTREE)
  {
    fprintf(rlab_stderr, THIS_SOLVER ": entry '"RLAB_NAME_FFNN_WEIGHTS"' has to be list of real matrices!\n");
    goto _exit_ent_ffnn_eval;
  }

  //
  node_ifwgt = btree_FindNode (ent_data (e2), RLAB_NAME_FFNN_FIX_WEIGHTS);
  if (node_ifwgt)
    if(ent_type(var_ent (node_ifwgt)) != BTREE)
    {
      fprintf(rlab_stderr, THIS_SOLVER ": optional entry '"RLAB_NAME_FFNN_FIX_WEIGHTS"' has to be list of real matrices!");
      node_ifwgt = 0;
    }

  //
  node_theta = btree_FindNode (ent_data (e2), RLAB_NAME_FFNN_BIAS);
  if (!node_theta)
    goto _exit_ent_ffnn_eval;
  if(ent_type(var_ent (node_theta)) != BTREE)
  {
    fprintf(rlab_stderr, THIS_SOLVER ": entry '"RLAB_NAME_FFNN_BIAS"' has to be list of real matrices!");
    goto _exit_ent_ffnn_eval;
  }

  //
  node_iftheta = btree_FindNode (ent_data (e2), RLAB_NAME_FFNN_FIX_BIAS);
  if (node_iftheta)
    if(ent_type(var_ent (node_iftheta)) != BTREE)
      fprintf(rlab_stderr, ": entry '"RLAB_NAME_FFNN_FIX_BIAS"' has to be list of real matrices!");

  //
  node_act = btree_FindNode (ent_data (e2), RLAB_NAME_FFNN_ACTIVATION);
  if (!node_act)
    goto _exit_ent_ffnn_eval;
  if(ent_type(var_ent (node_act)) != BTREE)
  {
    fprintf(rlab_stderr, ": entry '"RLAB_NAME_FFNN_ACTIVATION"' has to be list of real matrices!");
    goto _exit_ent_ffnn_eval;
  }

  //
  node_trans = btree_FindNode (ent_data (e2), RLAB_NAME_FFNN_TRANSFER);
  if (!node_trans)
    goto _exit_ent_ffnn_eval;
  if(ent_type(var_ent (node_trans)) != BTREE)
  {
    fprintf(rlab_stderr,THIS_SOLVER ": entry '"RLAB_NAME_FFNN_TRANSFER"' has to be list of real matrices!");
    goto _exit_ent_ffnn_eval;
  }

  L = SIZE(NetVect);
  if (L<2 || L>999)
  {
    fprintf(rlab_stderr,THIS_SOLVER ": entry 'nodes' has to be real vector!");
    goto _exit_ent_ffnn_eval;
  }

  int out_dim = mdiV0(NetVect,0);
  int inp_dim = mdiV0(NetVect,L-1); // this should be the dimension of the dataset

  if (out_dim >= inp_dim)
  {
    fprintf(rlab_stderr,THIS_SOLVER ": input dimension has to be greater than the output dimension!");
    goto _exit_ent_ffnn_eval;
  }

  if (inp_dim != MNC(x))
  {
    fprintf(rlab_stderr,THIS_SOLVER ": input dimension has to match the input data dimension!");
    goto _exit_ent_ffnn_eval;
  }

  //
  // this is used for passing data to the solver
  //
  NetPtr_wgt     = (MDR **) GC_malloc(L * sizeof(MDR*));
  NetPtr_ifwgt   = (MDR **) GC_malloc(L * sizeof(MDR*));
  NetPtr_theta   = (MDR **) GC_malloc(L * sizeof(MDR*));
  NetPtr_iftheta = (MDR **) GC_malloc(L * sizeof(MDR*));
  NetPtr_act     = (MDS **) GC_malloc(L * sizeof(MDS*));
  NetPtr_transf  = (MDS **) GC_malloc(L * sizeof(MDS*));

  //
  // fill in
  //
  for (i=0; i<(L-1); i++)
  {
    char node_name[4];
    int nr =  mdiV0(NetVect,i);
    int nc =  mdiV0(NetVect,i+1);
    sprintf(node_name,"%i",i+1);

    //
    NetPtr_wgt[i]   = 0;
    NetPtr_ifwgt[i] = 0;
    NetPtr_theta[i] = 0;
    NetPtr_iftheta[i] = 0;
    NetPtr_act[i] = 0;
    NetPtr_transf[i] = 0;

    // wgt
    node = btree_FindNode (ent_data (var_ent(node_wgt)), node_name);
    if (node)
    {
      if (ent_type(var_ent (node)) == MATRIX_DENSE_REAL)
      {
        NetPtr_wgt[i] = ent_data( var_ent(node) );
        if (MNR(NetPtr_wgt[i])!=nr || MNC(NetPtr_wgt[i])!=nc)
          NetPtr_wgt[i]=0;
      }
    }
    if (!NetPtr_wgt[i])
      goto _exit_ent_ffnn_eval;

    // ifwgt
    if (node_ifwgt)
    {
      node = btree_FindNode (ent_data (var_ent(node_ifwgt)), node_name);
      if (node)
      {
        if (ent_type(var_ent (node)) == MATRIX_DENSE_REAL)
        {
          NetPtr_ifwgt[i] = ent_data(var_ent(node) );
          if (MNR(NetPtr_ifwgt[i])!=nr || MNC(NetPtr_ifwgt[i])!=nc)
            NetPtr_ifwgt[i]=0;
        }
      }
    }

    // theta
    node = btree_FindNode (ent_data (var_ent(node_theta)), node_name);
    if (node)
    {
      if (ent_type(var_ent (node)) == MATRIX_DENSE_REAL)
      {
        NetPtr_theta[i] = ent_data(var_ent(node) );
        if (SIZE(NetPtr_theta[i])!=nr)
          NetPtr_theta[i]=0;
      }
    }
    if (!NetPtr_theta[i])
      goto _exit_ent_ffnn_eval;

    // fix_theta
    if (node_iftheta)
    {
      node = btree_FindNode (ent_data (var_ent(node_iftheta)), node_name);
      if (node)
      {
        if (ent_type(var_ent (node)) == MATRIX_DENSE_REAL)
        {
          NetPtr_iftheta[i] = ent_data(var_ent(node) );
          if (SIZE(NetPtr_theta[i])!=nr)
            NetPtr_iftheta[i]=0;
        }
      }
    }

    // act
    node = btree_FindNode (ent_data (var_ent(node_act)), node_name);
    if (node)
    {
      if (ent_type(var_ent (node)) == MATRIX_DENSE_STRING)
      {
        NetPtr_act[i] = ent_data( var_ent(node) );
        if (SIZE(NetPtr_act[i])!=nr)
          NetPtr_act[i]=0;
      }
    }
    if (!NetPtr_act[i])
      goto _exit_ent_ffnn_eval;

    // transf
    node = btree_FindNode (ent_data (var_ent(node_trans)), node_name);
    if (node)
    {
      if (ent_type(var_ent (node)) == MATRIX_DENSE_STRING)
      {
        NetPtr_transf[i] = ent_data( var_ent(node) );
        if (SIZE(NetPtr_transf[i])!=nr)
          NetPtr_transf[i]=0;
      }
    }
    if (!NetPtr_transf[i])
      goto _exit_ent_ffnn_eval;
  }

  // transpose for easier handling
  int num_samples = MNR(x);
  y_out = mdr_Create(out_dim,MNR(x));
  mdr_Transpose_inplace (x);

  //
  // need computational space
  //
  MDR ** y1 = GC_malloc(L*sizeof(MDR*));
  for (i=0; i<L-1; i++)
  {
    y1[i]  = mdr_Create(mdiV0(NetVect,i),1);
  }
  y1[L-1] = mdr_CreateEmpty(mdiV0(NetVect,L-1),1);  // input layer is special: we will use pointer from dset'

  for (i=0;i<num_samples; i++)
  {
    eval_ffnn ( NetVect, NetPtr_wgt, NetPtr_theta, NetPtr_act, NetPtr_transf, y1,
                &Mdr0(x,0,i), &Mdr0(y_out,0,i) );
  }

  // retranspose
  mdr_Transpose_inplace (y_out);
  mdr_Transpose_inplace (x);

  // cleanup computational space
  MDPTR(y1[L-1])=0;
  for (i=0; i<L; i++)
    mdr_Destroy(y1[i]);
  GC_FREE(y1);

_exit_ent_ffnn_eval:

  // clean-up
  ent_Clean (e1);
  ent_Clean (e2);

  for (i=0; i<L-1; i++)
  {
    NetPtr_wgt[i] = 0;
    if (NetPtr_ifwgt)
      NetPtr_ifwgt[i]=0;
    NetPtr_theta[i]=0;
    if (NetPtr_iftheta)
      NetPtr_iftheta[i]=0;
    NetPtr_act[i]=0;
    NetPtr_transf[i]=0;
  }
  GC_FREE(NetPtr_wgt);
  if (NetPtr_ifwgt)
    GC_FREE(NetPtr_ifwgt);
  GC_FREE(NetPtr_theta);
  if (NetPtr_iftheta)
    GC_FREE(NetPtr_iftheta);
  GC_FREE(NetPtr_act);
  GC_FREE(NetPtr_transf);

  // install final output of the network
  return ent_Assign_Rlab_MDR(y_out);
}


//
// samann neural network trainer for dimension data reduction
//
#undef THIS_SOLVER
#define THIS_SOLVER "samann"
Ent * ent_ffnn_samann (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *e3=0, *eo=0, *rent;
  MDR *x=0, *NetVect=0, *mu=0;
  FILE *rlab_stderr = (!RLAB_STDERR_DS) ? stderr : RLAB_STDERR_DS;

  //
  // backpropagation parameters for optimization
  //
  double lrate     = 0.1;
  double lmomentum = 0.9;

  int i, j, idummy, isc=1, isolver=1;
  double ddummy, delta_stress_threshold=0.0, stress_threshold=0, wgt_rcpdist=0;
  int L, options=RLAB_SAMANN_ALL_PAIRS, maxi=1000, lout=0, ncompstress=1;

  //
  // GSL/simplex minimizer
  //
  double s0 = 5e-1;
  double abserr = 1e-12;
  int num_attempts=20;

  //
  // NEWUOA minimizer
  //
  double rhobeg = 1e-1;
  double rhoend = 1e-3;

  FILE *fptr = 0;
  ListNode * node=0, * node_wgt=0, * node_ifwgt=0, * node_theta=0;
  ListNode * node_iftheta=0, * node_act=0, * node_trans=0;
  char *outstream=0;

  double scale = 1.0;

  Btree * bw = 0;
  Btree * bt = 0;
  Btree * biftheta = 0;
  Btree * bifwgt   = 0;
  Btree * bact     = 0;
  Btree * btransf  = 0;

  //
  // Load and Check arguments.
  //
  if (nargs != 2 && nargs != 3)
  {
    fprintf(rlab_stderr,
            THIS_SOLVER ": Realization of samann mapping via neural network.\n"
            THIS_SOLVER ": Format :\n"
           );
    fprintf(rlab_stderr,
            THIS_SOLVER ":   tn = ffnn.samann(x,nn/,nnopts/)\n"
            THIS_SOLVER ": See rlabplus manual for more information.\n");
    rerror
        (THIS_SOLVER ": requires three or four arguments");
  }

  //
  // x
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) != MATRIX_DENSE_REAL)
    rerror (THIS_SOLVER ": 'X' has to be a real matrix!");
  x = class_matrix_real (e1);

  //
  // second argument is the nn structured list or the outline of the network
  //
  e2 = bltin_get_ent (args[1]);
  if (ent_type (e2) == BTREE)
  {
    //
    // create a feed forward neural network using information from
    // 'nodes'.
    //
    node = btree_FindNode (ent_data (e2), RLAB_NAME_FFNN_NODES);
    if (!node)
      rerror (THIS_SOLVER ": entry '"RLAB_NAME_FFNN_NODES"' has to be real vector!");
    if(ent_type(var_ent (node))!=MATRIX_DENSE_REAL)
      rerror(THIS_SOLVER ": entry 'nodes' has to be real vector!");

    NetVect = mdr_Int_BF( ent_data (var_ent (node)) );

    //
    // fill the network using data from the lists 'wgt', 'fix_wgt', 'bias', 'fix_bias',
    // 'act' and 'trans'
    //
    node_wgt = btree_FindNode (ent_data (e2), RLAB_NAME_FFNN_WEIGHTS);
    if(ent_type(var_ent (node_wgt)) != BTREE)
      rerror(THIS_SOLVER ": entry '"RLAB_NAME_FFNN_WEIGHTS"' has to be list of real matrices!");
    //
    node_ifwgt = btree_FindNode (ent_data (e2), RLAB_NAME_FFNN_FIX_WEIGHTS);
    if(ent_type(var_ent (node_ifwgt)) != BTREE)
      node_ifwgt=0;

    //
    node_theta = btree_FindNode (ent_data (e2), RLAB_NAME_FFNN_BIAS);
    if(ent_type(var_ent (node_theta)) != BTREE)
      rerror(THIS_SOLVER ": entry '"RLAB_NAME_FFNN_BIAS"' has to be list of real matrices!");
    //
    node_iftheta = btree_FindNode (ent_data (e2), RLAB_NAME_FFNN_FIX_BIAS);
    if(ent_type(var_ent (node_iftheta)) != BTREE)
      node_iftheta=0;

    //
    node_act = btree_FindNode (ent_data (e2), RLAB_NAME_FFNN_ACTIVATION);
    if(ent_type(var_ent (node_act)) != BTREE)
      rerror(THIS_SOLVER ": entry '"RLAB_NAME_FFNN_ACTIVATION"' has to be list of real matrices!");
    //
    node_trans = btree_FindNode (ent_data (e2), RLAB_NAME_FFNN_TRANSFER);
    if(ent_type(var_ent (node_trans)) != BTREE)
      rerror(THIS_SOLVER ": entry '"RLAB_NAME_FFNN_TRANSFER"' has to be list of real matrices!");
  }
  else if (ent_type (e2) == MATRIX_DENSE_REAL)
  {
    // if the network is initialized through node count then do random initialization
    NetVect = mdr_Int_BF( ent_data (e2) );
    options |= RLAB_FFNN_INIT_RAND;
  }
  else
    rerror(THIS_SOLVER ": second argument should be array or nn-list");

  L = SIZE(NetVect);
  if (L<2)
    rerror(THIS_SOLVER ": entry 'nodes' has to be real vector!");

  int out_dim = mdiV0(NetVect,0);
  int inp_dim = mdiV0(NetVect,L-1); // this should be the dimension of the dataset

  if (out_dim >= inp_dim)
    rerror(THIS_SOLVER ": input dimension has to be greater than the output dimension!");

  if (inp_dim != MNC(x))
    rerror(THIS_SOLVER ": input dimension has to match the input data dimension!");

  //
  // this is used for passing data to the solver
  //
  MDR **NetPtr_wgt     = (MDR **) GC_malloc(L * sizeof(MDR*));
  MDR **NetPtr_ifwgt   = (MDR **) GC_malloc(L * sizeof(MDR*));
  MDR **NetPtr_theta   = (MDR **) GC_malloc(L * sizeof(MDR*));
  MDR **NetPtr_iftheta = (MDR **) GC_malloc(L * sizeof(MDR*));
  MDS **NetPtr_act     = (MDS **) GC_malloc(L * sizeof(MDS*));
  MDS **NetPtr_transf  = (MDS **) GC_malloc(L * sizeof(MDS*));

  //
  // fill in
  //
  for (i=0; i<(L-1); i++)
  {
    char node_name[4];
    int nr =  mdiV0(NetVect,i);
    int nc =  mdiV0(NetVect,i+1);
    sprintf(node_name,"%i",i+1);

    //
    NetPtr_wgt[i]   = 0;
    NetPtr_ifwgt[i] = 0;
    NetPtr_theta[i] = 0;
    NetPtr_iftheta[i] = 0;
    NetPtr_act[i] = 0;
    NetPtr_transf[i] = 0;

    // wgt
    if (node_wgt)
    {
      node = btree_FindNode (ent_data (var_ent(node_wgt)), node_name);
      if (node)
      {
        if (ent_type(var_ent (node)) == MATRIX_DENSE_REAL)
        {
          MDR *dummy = ent_data(var_ent(node));
          if (MNR(dummy)==nr && MNC(dummy)==nc)
            NetPtr_wgt[i] = mdr_Float_BF(dummy) ;
        }
      }
    }
    if (!NetPtr_wgt[i])
      NetPtr_wgt[i] = mdr_Create(nr, nc);

    // fix_wgt
    if (node_ifwgt)
    {
      node = btree_FindNode (ent_data (var_ent(node_ifwgt)), node_name);
      if (node)
      {
        if (ent_type(var_ent (node)) == MATRIX_DENSE_REAL)
        {
          MDR *dummy = ent_data( var_ent(node) );
          if (MNR(dummy)==nr && MNC(dummy)==nc)
          {
            NetPtr_ifwgt[i] = mdr_Int_BF(dummy);
          }
        }
      }
    }

    // theta
    if (node_theta)
    {
      node = btree_FindNode (ent_data (var_ent(node_theta)), node_name);
      if (node)
      {
        if (ent_type(var_ent (node)) == MATRIX_DENSE_REAL)
        {
          MDR *dummy = ent_data(var_ent(node) );
          if (SIZE(NetPtr_theta[i])==nr)
          {
            NetPtr_theta[i] = mdr_Float_BF(dummy);
          }
        }
      }
    }
    if (!NetPtr_theta[i])
      NetPtr_theta[i] = mdr_Create(nr, 1);

    // fix_theta
    if (node_iftheta)
    {
      node = btree_FindNode (ent_data (var_ent(node_iftheta)), node_name);
      if (node)
      {
        if (ent_type(var_ent (node)) == MATRIX_DENSE_REAL)
        {
          MDR *dummy = ent_data(var_ent(node) );
          if (SIZE(NetPtr_iftheta[i])==nr)
          {
            NetPtr_iftheta[i] = mdr_Int_BF(dummy);
          }
        }
      }
    }

    // act
    if (node_act)
    {
      node = btree_FindNode (ent_data (var_ent(node_act)), node_name);
      if (node)
      {
        if (ent_type(var_ent (node)) == MATRIX_DENSE_STRING)
        {
          int icopy=1;
          MDS *dummy = ent_data( var_ent(node) );
          if (SIZE(dummy)==nr)
          {
            for (j=0; j<nr; j++)
            {
              if (MdsV0(dummy,j)[0] != 'i')
              {
                icopy=0;
                break;
              }
            }
            if (icopy)
              NetPtr_act[i] = mds_Copy(dummy);
          }
        }
      }
    }
    if (!NetPtr_act[i])
    {
      NetPtr_act[i] = mds_Create(nr, 1);
      for (j=0; j<nr; j++)
        MdsV0(NetPtr_act[i],j) = cpstr("i");
    }

    // transf
    if (node_trans)
    {
      node = btree_FindNode (ent_data (var_ent(node_trans)), node_name);
      if (node)
      {
        if (ent_type(var_ent (node)) == MATRIX_DENSE_STRING)
        {
          int icopy=1;
          MDS *dummy = ent_data( var_ent(node) );
          if (SIZE(dummy)==nr)
          {
            for (j=0; j<nr; j++)
            {
              if (MdsV0(dummy,j)[0] != 's')
              {
                icopy=0;
                break;
              }
            }
            if (icopy)
              NetPtr_transf[i] = mds_Copy(dummy);
          }
        }
      }
    }
    if (!NetPtr_transf[i])
    {
      NetPtr_transf[i] = mds_Create(nr, 1);
      for (j=0; j<nr; j++)
        MdsV0(NetPtr_transf[i],j) = cpstr("s");
    }

  }

  //
  // options
  //
  if (nargs == 3)
  {
    eo = bltin_get_ent (args[2]);
    if (ent_type (eo) == BTREE)
    {
      // solver:
      // standard:
      //  0 - backpropagation
      // proof of principle
      //  1 - simplex minimizer
      //  2 - newuoa
      node = btree_FindNode (ent_data (eo), RLAB_NAME_FFNN_SAMANN_ISOLVER);
      if (node != 0)
      {
        idummy = (int) class_double ( var_ent (node) );
        if (idummy == 1)
          isolver = 1;
        else if (idummy == 2)
          isolver = 2;
        else
          isolver = 0;
      }

      // scale the dataset?
      node = btree_FindNode (ent_data (eo), RLAB_NAME_FFNN_SAMANN_SCALE);
      if (node != 0)
      {
        ddummy = class_double ( var_ent (node) );
        isc = (ddummy > 0);
      }

      //
      // backpropagation parameters
      //
      // learning rate
      node = btree_FindNode (ent_data (eo), RLAB_NAME_PATTERN_BACKPROP_MOMENTUM);
      if (node != 0)
      {
        ddummy = class_double ( var_ent (node) );
        if (ddummy > 0)
        {
          if (isolver!=0)
            fprintf(rlab_stderr, THIS_SOLVER ": provided is parameter '"
                RLAB_NAME_PATTERN_BACKPROP_MOMENTUM
                "' but the method is not backpropagation: Parameter is ignored!\n");
          else
            lrate = ddummy;
        }
      }

      // learning momentum
      node = btree_FindNode (ent_data (eo), RLAB_NAME_PATTERN_BACKPROP_LEARNING_RATE);
      if (node != 0)
      {
        ddummy = class_double ( var_ent (node) );
        if (ddummy > 0)
        {
          if (isolver!=0)
            fprintf(rlab_stderr, THIS_SOLVER ": provided is parameter '"
                RLAB_NAME_PATTERN_BACKPROP_LEARNING_RATE
                "' but the method is not backpropagation: Parameter is ignored!\n");
          else
            lmomentum = ddummy;
        }
      }

      // sammon stress threshold: stop if we reach it
      node = btree_FindNode (ent_data (eo), RLAB_NAME_FFNN_ERR_THRESHOLD);
      if (node != 0)
      {
        ddummy = (int) class_double ( var_ent (node) );
        if (ddummy > 0)
          stress_threshold = ddummy;
      }

      // sammon delta stress threshold: stop if we reach it
      node = btree_FindNode (ent_data (eo), RLAB_NAME_FFNN_ERR_DELTA);
      if (node != 0)
      {
        ddummy = class_double ( var_ent (node) );
        if (ddummy > 0)
          delta_stress_threshold = ddummy;
      }

      // how often compute sammon stress
      node = btree_FindNode (ent_data (eo), RLAB_NAME_FFNN_NSTEP_ERR_COMP);
      if (node != 0)
      {
        idummy = (int) class_double ( var_ent (node) );
        if (idummy > 0)
          ncompstress = idummy;
      }

      // chose random pairs or go sequentially
      //  0 - iterate over entire dataset
      //  1 - iterate over random pairs
      node = btree_FindNode (ent_data (eo), RLAB_NAME_FFNN_SAMANN_RAND_PAIRS);
      if (node != 0)
      {
        ddummy = class_double ( var_ent (node) );
        if (ddummy > 0)
        {
          options |=  RLAB_SAMANN_RANDOM_PAIR;
          options &= ~RLAB_SAMANN_ALL_PAIRS;
        }
        else
        {
          options |=  RLAB_SAMANN_ALL_PAIRS;
          options &= ~RLAB_SAMANN_RANDOM_PAIR;
        }
      }

      // maximum number of iterations
      node = btree_FindNode (ent_data (eo), RLAB_NAME_GEN_MAXITER);
      if (node != 0)
      {
        idummy = (int) class_double ( var_ent (node) );
        if (idummy > 0)
          maxi = idummy;
      }

      // standard output
      node = btree_FindNode (ent_data (eo), RLAB_NAME_GEN_STDOUT);
      if (node != 0)
      {
        if (ent_type(var_ent (node))== MATRIX_DENSE_STRING)
        {
          outstream = class_char_pointer(var_ent (node));
          lout      = isvalidstring( outstream );
        }
      }

      // do random initialization?
      node = btree_FindNode (ent_data (eo), RLAB_NAME_FFNN_INIT_RAND);
      if (node != 0)
      {
        idummy = (int) class_double ( var_ent (node) );
        if (idummy>0)
          options |= RLAB_FFNN_INIT_RAND;
      }

      // reciprocal distance weight: modify weight of distance in sammon stress function
      node = btree_FindNode (ent_data (eo), RLAB_NAME_FFNN_SAMANN_WGT_RCP_DIST);
      if (node != 0)
      {
        ddummy = class_double ( var_ent (node) );
        if (ddummy >= 0)
          wgt_rcpdist = ddummy;
      }


      // GSL simplex parameters: simplex size 'ss'
      node = btree_FindNode (ent_data (eo), RLAB_NAME_FFNN_SAMANN_SIMPLEX_S0);
      if (node != 0)
      {
        ddummy = class_double ( var_ent (node) );
        if (ddummy > 0)
          s0 = ddummy;
      }

      // GSL simplex parameters: initial simplex size 'ss'
      node = btree_FindNode (ent_data (eo), RLAB_NAME_FFNN_SAMANN_SIMPLEX_S0);
      if (node != 0)
      {
        ddummy = class_double ( var_ent (node) );
        if (ddummy > 0)
          s0 = ddummy;
      }

      // GSL simplex parameters: eabs
      node = btree_FindNode (ent_data (eo), RLAB_NAME_FFNN_SAMANN_SIMPLEX_EABS);
      if (node != 0)
      {
        ddummy = class_double ( var_ent (node) );
        if (ddummy > 0)
          abserr = ddummy;
      }

      // GSL simplex parameters: batch length over which sammon's stress is minimized
      node = btree_FindNode (ent_data (eo), RLAB_NAME_FFNN_SAMANN_NTRIES);
      if (node != 0)
      {
        ddummy = class_double ( var_ent (node) );
        if (ddummy > 0)
          num_attempts = ddummy;
      }

      // NEWUOA: rhobeg
      node = btree_FindNode (ent_data (eo), RLAB_NAME_NUO_RHOBEG);
      if (node != 0)
      {
        ddummy = class_double ( var_ent (node) );
        if (ddummy > 0)
          rhobeg = ddummy;
      }

      node = btree_FindNode (ent_data (eo), RLAB_NAME_NUO_RHOEND);
      if (node != 0)
      {
        ddummy = class_double ( var_ent (node) );
        if (ddummy > 0 && rhoend < rhobeg)
          rhoend = ddummy;
        else
          rhoend = 0.1 * rhobeg;
      }
    }
    ent_Clean (eo);
  }

  if (lout > 1)
  {
    fptr = fopen (outstream, "a");
  }

  // copy x to lset and transpose it, as the latter will be heavily modified:
  MDR *lset = mdr_Transpose(x);
  MDR *lset_dist = mdr_Create(MNC(lset),MNC(lset));
  double rcp_lambda=1.0, stress=1e99;

  if (options & RLAB_FFNN_INIT_RAND)
  {
    // do random init of the network weights and biases if requested
    ffnn_randinit_wgt_theta(NetVect, NetPtr_wgt, NetPtr_ifwgt, NetPtr_theta, NetPtr_iftheta);
  }

  if (isc)
  {
    // rescale learning set so that network can work on it
    scale = samann_init (lset, lset_dist, &mu, NetVect, NetPtr_act, NetPtr_transf, &rcp_lambda);
    mdr_Transpose_inplace (mu);
  }

  MDR *y_out = mdr_Create(out_dim,MNR(x));

  //
  // main backpropagation iteration loop
  //
  switch (isolver)
  {
    case 0:
      // backpropagation
      samann_learn_backprop (lset, lset_dist, NetVect,
                             NetPtr_wgt, NetPtr_ifwgt, NetPtr_theta, NetPtr_iftheta,
                             NetPtr_act, NetPtr_transf, lrate, lmomentum, maxi, options, rcp_lambda, wgt_rcpdist,
                             &stress, fptr, y_out, ncompstress, delta_stress_threshold, stress_threshold);
      break;

    case 1:
      // gsl simplex minimizer
      samann_learn_simplex (lset, lset_dist, NetVect,
                            NetPtr_wgt, NetPtr_ifwgt, NetPtr_theta, NetPtr_iftheta,
                            NetPtr_act, NetPtr_transf,
                            s0, abserr, num_attempts, maxi, options, rcp_lambda, wgt_rcpdist,
                            &stress, fptr, y_out, ncompstress, delta_stress_threshold, stress_threshold);
      break;

    case 2:
      // newuoa minimizer
      samann_learn_newuoa (lset, lset_dist, NetVect,
                           NetPtr_wgt, NetPtr_ifwgt, NetPtr_theta, NetPtr_iftheta,
                           NetPtr_act, NetPtr_transf,
                           rhobeg, rhoend, num_attempts, maxi, options, rcp_lambda, wgt_rcpdist,
                           &stress, fptr, y_out, ncompstress, delta_stress_threshold, stress_threshold);
      break;

  }
  mdr_Transpose_inplace (y_out);

  //
  // write everything to the output list
  //
  Btree * rtree = btree_Create ();
  bw       = btree_Create ();
  if (NetPtr_ifwgt)
    bifwgt   = btree_Create ();
  bt       = btree_Create ();
  if (NetPtr_iftheta)
    biftheta = btree_Create ();
  bact     = btree_Create ();
  btransf  = btree_Create ();

  for (i=0; i<L-1; i++)
  {
    char nnode[4];

    sprintf(nnode,"%i",i+1);
    install (bw, nnode, ent_Assign_Rlab_MDR (NetPtr_wgt[i]));
    if (NetPtr_ifwgt)
      if (NetPtr_ifwgt[i])
        install (bifwgt, nnode, ent_Assign_Rlab_MDR (NetPtr_ifwgt[i]));

    install (bt, nnode, ent_Assign_Rlab_MDR (NetPtr_theta[i]));
    if (NetPtr_iftheta)
      install (biftheta, nnode, ent_Assign_Rlab_MDR (NetPtr_iftheta[i]));

    install (bact,    nnode, ent_Assign_Rlab_MDS (mds_Copy(NetPtr_act[i])));
    install (btransf, nnode, ent_Assign_Rlab_MDS (mds_Copy(NetPtr_transf[i])));
  }

  install(rtree, RLAB_NAME_FFNN_WEIGHTS, ent_Assign_Rlab_BTREE(bw));
  if (NetPtr_ifwgt)
    install(rtree, RLAB_NAME_FFNN_FIX_WEIGHTS, ent_Assign_Rlab_BTREE(bifwgt));

  install(rtree, RLAB_NAME_FFNN_BIAS, ent_Assign_Rlab_BTREE(bt));
  if (NetPtr_iftheta)
    install(rtree, RLAB_NAME_FFNN_FIX_BIAS, ent_Assign_Rlab_BTREE(biftheta));

  install(rtree, RLAB_NAME_FFNN_ACTIVATION, ent_Assign_Rlab_BTREE(bact));
  install(rtree, RLAB_NAME_FFNN_TRANSFER, ent_Assign_Rlab_BTREE(btransf));

  // install nodes array
  install(rtree, RLAB_NAME_FFNN_NODES, ent_Assign_Rlab_MDR(NetVect));

  // did we rescale the sets?
  if (isc)
  {
    // install samann scale factors
    install(rtree, RLAB_NAME_FFNN_SAMANN_SCALE, ent_Create_Rlab_Double(scale));
    // install samann scale mean
    install(rtree, RLAB_NAME_FFNN_SAMANN_MEAN, ent_Assign_Rlab_MDR(mu));
  }

  // install samann stress factors
  install(rtree, RLAB_NAME_FFNN_SAMANN_STRESS, ent_Create_Rlab_Double(stress));
  // install final output of the network
  install(rtree, RLAB_NAME_FFNN_SAMANN_Y_OUT, ent_Assign_Rlab_MDR(y_out));

  // clean-up
  ent_Clean (e1);
  ent_Clean (e2);
  ent_Clean (e3);

  mdr_Destroy(lset);
  mdr_Destroy(lset_dist);
  for (i=0; i<L-1; i++)
  {
    NetPtr_wgt[i] = 0;
    if (NetPtr_ifwgt)
      NetPtr_ifwgt[i]=0;
    NetPtr_theta[i]=0;
    if (NetPtr_iftheta)
      NetPtr_iftheta[i]=0;
    NetPtr_act[i]=0;
    NetPtr_transf[i]=0;
  }
  GC_FREE(NetPtr_wgt);
  if (NetPtr_ifwgt)
    GC_FREE(NetPtr_ifwgt);
  GC_FREE(NetPtr_theta);
  if (NetPtr_iftheta)
    GC_FREE(NetPtr_iftheta);
  GC_FREE(NetPtr_act);
  GC_FREE(NetPtr_transf);

  //
  // write result to output list
  //
  rent = ent_Create ();
  ent_data (rent) = rtree;
  ent_type (rent) = BTREE;
  return rent;
}

//
// samann neural network trainer for dimension data reduction
//
#undef THIS_SOLVER
#define THIS_SOLVER "ffnn"
Ent * ent_ffnn (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *e3=0, *eo=0, *rent;
  MDR *NetVect=0, *mu=0, *tset=0, *dset=0;
  FILE *rlab_stderr = (!RLAB_STDERR_DS) ? stderr : RLAB_STDERR_DS;

  MDR **NetPtr_wgt, **NetPtr_ifwgt, **NetPtr_theta, **NetPtr_iftheta;
  MDS **NetPtr_act, **NetPtr_transf;

  //
  // backpropagation parameters for optimization
  //
  double lrate     = 0.1;
  double lmomentum = 0.9;

  int i, j, idummy, isc=1, isolver=1;
  double ddummy, delta_mse_threshold=0.0, mse_threshold=0;
  int L, options=RLAB_SAMANN_ALL_PAIRS, maxi=1000, lout=0, ncompmse=1;

  //
  // GSL simplex minimizer
  //
//   double abserr = 1e-6;
//   int num_attempts = 20;
//   double s0=1e-3;

  //
  // NEWUOA minimizer
  //
  double rhobeg = 1e-1;
  double rhoend = 1e-3;

  FILE *fptr = 0;
  ListNode * node=0, * node_wgt=0, * node_ifwgt=0, * node_theta=0;
  ListNode * node_iftheta=0, * node_act=0, * node_trans=0;
  char *outstream=0;

  double scale = 1.0;

  Btree * bw = 0;
  Btree * bt = 0;
  Btree * biftheta = 0;
  Btree * bifwgt   = 0;
  Btree * bact     = 0;
  Btree * btransf  = 0;

  //
  // Load and Check arguments.
  //
  if (nargs != 2 && nargs != 3)
  {
    fprintf(rlab_stderr,
            THIS_SOLVER ": Realization of samann mapping via neural network.\n"
            THIS_SOLVER ": Format :\n"
           );
    fprintf(rlab_stderr,
            THIS_SOLVER ":   tn = ffnn(x,nn/,nnopts/)\n"
            THIS_SOLVER ": See rlabplus manual for more information.\n");
    rerror
        (THIS_SOLVER ": requires three or four arguments");
  }

  //
  // <<data;feat>>
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) != BTREE)
    rerror (THIS_SOLVER ": " RLAB_ERROR_ARG1_FFNN);
  // get data to 'dset'
  node = btree_FindNode (ent_data (e1), RLAB_NAME_PATTERN_DATA);
  if (!node)
    rerror (THIS_SOLVER ": entry '"RLAB_NAME_PATTERN_DATA"' has to be real vector!");
  if(ent_type(var_ent (node))!=MATRIX_DENSE_REAL)
    rerror(THIS_SOLVER ": entry '"RLAB_NAME_PATTERN_DATA"' has to be real vector!");
  dset = ent_data(var_ent (node));
  // get feature to 'tset'
  node = btree_FindNode (ent_data (e1), RLAB_NAME_PATTERN_FEATURE);
  if (!node)
    rerror (THIS_SOLVER ": entry '"RLAB_NAME_PATTERN_FEATURE"' has to be real vector!");
  if(ent_type(var_ent (node))!=MATRIX_DENSE_REAL)
    rerror(THIS_SOLVER ": entry '"RLAB_NAME_PATTERN_FEATURE"' has to be real vector!");
  tset = ent_data(var_ent (node));

  //
  // second argument is the nn structured list or the outline of the network
  //
  e2 = bltin_get_ent (args[1]);
  if (ent_type (e2) == BTREE)
  {
    //
    // create a feed forward neural network using information from
    // 'nodes'.
    //
    node = btree_FindNode (ent_data (e2), RLAB_NAME_FFNN_NODES);
    if (!node)
      rerror (THIS_SOLVER ": entry '"RLAB_NAME_FFNN_NODES"' has to be real vector!");
    if(ent_type(var_ent (node))!=MATRIX_DENSE_REAL)
      rerror(THIS_SOLVER ": entry 'nodes' has to be real vector!");

    NetVect = mdr_Int_BF( ent_data (var_ent (node)) );

    //
    // fill the network using data from the lists 'wgt', 'fix_wgt', 'bias', 'fix_bias',
    // 'act' and 'trans'
    //
    node_wgt = btree_FindNode (ent_data (e2), RLAB_NAME_FFNN_WEIGHTS);
    if(ent_type(var_ent (node_wgt)) != BTREE)
      rerror(THIS_SOLVER ": entry '"RLAB_NAME_FFNN_WEIGHTS"' has to be list of real matrices!");
    //
    node_ifwgt = btree_FindNode (ent_data (e2), RLAB_NAME_FFNN_FIX_WEIGHTS);
    if(ent_type(var_ent (node_ifwgt)) != BTREE)
      node_ifwgt=0;

    //
    node_theta = btree_FindNode (ent_data (e2), RLAB_NAME_FFNN_BIAS);
    if(ent_type(var_ent (node_theta)) != BTREE)
      rerror(THIS_SOLVER ": entry '"RLAB_NAME_FFNN_BIAS"' has to be list of real matrices!");
    //
    node_iftheta = btree_FindNode (ent_data (e2), RLAB_NAME_FFNN_FIX_BIAS);
    if(ent_type(var_ent (node_iftheta)) != BTREE)
      node_iftheta=0;

    //
    node_act = btree_FindNode (ent_data (e2), RLAB_NAME_FFNN_ACTIVATION);
    if(ent_type(var_ent (node_act)) != BTREE)
      rerror(THIS_SOLVER ": entry '"RLAB_NAME_FFNN_ACTIVATION"' has to be list of real matrices!");
    //
    node_trans = btree_FindNode (ent_data (e2), RLAB_NAME_FFNN_TRANSFER);
    if(ent_type(var_ent (node_trans)) != BTREE)
      rerror(THIS_SOLVER ": entry '"RLAB_NAME_FFNN_TRANSFER"' has to be list of real matrices!");
  }
  else if (ent_type (e2) == MATRIX_DENSE_REAL)
  {
    // if the network is initialized through node count then do random initialization
    NetVect = mdr_Int_BF( ent_data (e2) );
    options |= RLAB_FFNN_INIT_RAND;
  }
  else
    rerror(THIS_SOLVER ": second argument should be array or nn-list");

  L = SIZE(NetVect);
  if (L<2 || L>999)
    rerror(THIS_SOLVER ": entry 'nodes' has to be real vector!");

  int out_dim = mdiV0(NetVect,0);
  int inp_dim = mdiV0(NetVect,L-1); // this should be the dimension of the dataset

  if (out_dim >= inp_dim)
    rerror(THIS_SOLVER ": input dimension has to be greater than the output dimension!");

  if (inp_dim != MNC(dset))
    rerror(THIS_SOLVER ": input dimension has to match the input data dimension!");

  if (out_dim != MNC(tset))
    rerror(THIS_SOLVER ": output dimension has to match the output data dimension!");

  //
  // this is used for passing data to the solver
  //
  NetPtr_wgt     = (MDR **) GC_malloc(L * sizeof(MDR*));
  NetPtr_ifwgt   = (MDR **) GC_malloc(L * sizeof(MDR*));
  NetPtr_theta   = (MDR **) GC_malloc(L * sizeof(MDR*));
  NetPtr_iftheta = (MDR **) GC_malloc(L * sizeof(MDR*));
  NetPtr_act     = (MDS **) GC_malloc(L * sizeof(MDS*));
  NetPtr_transf  = (MDS **) GC_malloc(L * sizeof(MDS*));

  //
  // fill in
  //
  for (i=0; i<(L-1); i++)
  {
    char node_name[4];
    int nr =  mdiV0(NetVect,i);
    int nc =  mdiV0(NetVect,i+1);
    sprintf(node_name,"%i",i+1);

    //
    NetPtr_wgt[i]   = 0;
    NetPtr_ifwgt[i] = 0;
    NetPtr_theta[i] = 0;
    NetPtr_iftheta[i] = 0;
    NetPtr_act[i] = 0;
    NetPtr_transf[i] = 0;

    // wgt
    if (node_wgt)
    {
      node = btree_FindNode (ent_data (var_ent(node_wgt)), node_name);
      if (node)
      {
        if (ent_type(var_ent (node)) == MATRIX_DENSE_REAL)
        {
          MDR *dummy = ent_data(var_ent(node));
          if (MNR(dummy)==nr && MNC(dummy)==nc)
            NetPtr_wgt[i] = mdr_Float_BF(dummy) ;
        }
      }
    }
    if (!NetPtr_wgt[i])
      NetPtr_wgt[i] = mdr_Create(nr, nc);

    // fix_wgt
    if (node_ifwgt)
    {
      node = btree_FindNode (ent_data (var_ent(node_ifwgt)), node_name);
      if (node)
      {
        if (ent_type(var_ent (node)) == MATRIX_DENSE_REAL)
        {
          MDR *dummy = ent_data( var_ent(node) );
          if (MNR(dummy)==nr && MNC(dummy)==nc)
          {
            NetPtr_ifwgt[i] = mdr_Int_BF(dummy);
          }
        }
      }
    }

    // theta
    if (node_theta)
    {
      node = btree_FindNode (ent_data (var_ent(node_theta)), node_name);
      if (node)
      {
        if (ent_type(var_ent (node)) == MATRIX_DENSE_REAL)
        {
          MDR *dummy = ent_data(var_ent(node) );
          if (SIZE(NetPtr_theta[i])==nr)
          {
            NetPtr_theta[i] = mdr_Float_BF(dummy);
          }
        }
      }
    }
    if (!NetPtr_theta[i])
      NetPtr_theta[i] = mdr_Create(nr, 1);

    // fix_theta
    if (node_iftheta)
    {
      node = btree_FindNode (ent_data (var_ent(node_iftheta)), node_name);
      if (node)
      {
        if (ent_type(var_ent (node)) == MATRIX_DENSE_REAL)
        {
          MDR *dummy = ent_data(var_ent(node) );
          if (SIZE(NetPtr_iftheta[i])==nr)
          {
            NetPtr_iftheta[i] = mdr_Int_BF(dummy);
          }
        }
      }
    }

    // act
    if (node_act)
    {
      node = btree_FindNode (ent_data (var_ent(node_act)), node_name);
      if (node)
      {
        if (ent_type(var_ent (node)) == MATRIX_DENSE_STRING)
        {
          int icopy=1;
          MDS *dummy = ent_data( var_ent(node) );
          if (SIZE(dummy)==nr)
          {
            for (j=0; j<nr; j++)
            {
              if (MdsV0(dummy,j)[0] != 'i')
              {
                icopy=0;
                break;
              }
            }
            if (icopy)
              NetPtr_act[i] = mds_Copy(dummy);
          }
        }
      }
    }
    if (!NetPtr_act[i])
    {
      NetPtr_act[i] = mds_Create(nr, 1);
      for (j=0; j<nr; j++)
        MdsV0(NetPtr_act[i],j) = cpstr("i");
    }

    // transf
    if (node_trans)
    {
      node = btree_FindNode (ent_data (var_ent(node_trans)), node_name);
      if (node)
      {
        if (ent_type(var_ent (node)) == MATRIX_DENSE_STRING)
        {
          int icopy=1;
          MDS *dummy = ent_data( var_ent(node) );
          if (SIZE(dummy)==nr)
          {
            for (j=0; j<nr; j++)
            {
              if (MdsV0(dummy,j)[0] != 's')
              {
                icopy=0;
                break;
              }
            }
            if (icopy)
              NetPtr_transf[i] = mds_Copy(dummy);
          }
        }
      }
    }
    if (!NetPtr_transf[i])
    {
      NetPtr_transf[i] = mds_Create(nr, 1);
      for (j=0; j<nr; j++)
        MdsV0(NetPtr_transf[i],j) = cpstr("s");
    }

  }

  //
  // options
  //
  if (nargs == 3)
  {
    eo = bltin_get_ent (args[2]);
    if (ent_type (eo) == BTREE)
    {
      // solver:
      // standard:
      //  0 - backpropagation
      // proof of principle
      //  1 - simplex minimizer
      //  2 - newuoa
      node = btree_FindNode (ent_data (eo), RLAB_NAME_FFNN_SAMANN_ISOLVER);
      if (node != 0)
      {
        idummy = (int) class_double ( var_ent (node) );
        if (idummy == 1)
          isolver = 1;
        else if (idummy == 2)
          isolver = 2;
        else
          isolver = 0;
      }

      // scale the dataset?
      node = btree_FindNode (ent_data (eo), RLAB_NAME_FFNN_SAMANN_SCALE);
      if (node != 0)
      {
        ddummy = class_double ( var_ent (node) );
        isc = (ddummy > 0);
      }

      //
      // backpropagation parameters
      //
      // learning rate
      node = btree_FindNode (ent_data (eo), RLAB_NAME_PATTERN_BACKPROP_MOMENTUM);
      if (node != 0)
      {
        ddummy = class_double ( var_ent (node) );
        if (ddummy > 0)
        {
          if (isolver!=0)
            fprintf(rlab_stderr, THIS_SOLVER ": provided is parameter '"
                RLAB_NAME_PATTERN_BACKPROP_MOMENTUM
                "' but the method is not backpropagation: Parameter is ignored!\n");
          else
            lrate = ddummy;
        }
      }

      // learning momentum
      node = btree_FindNode (ent_data (eo), RLAB_NAME_PATTERN_BACKPROP_LEARNING_RATE);
      if (node != 0)
      {
        ddummy = class_double ( var_ent (node) );
        if (ddummy > 0)
        {
          if (isolver!=0)
            fprintf(rlab_stderr, THIS_SOLVER ": provided is parameter '"
                RLAB_NAME_PATTERN_BACKPROP_LEARNING_RATE
                "' but the method is not backpropagation: Parameter is ignored!\n");
          else
            lmomentum = ddummy;
        }
      }

      // mse threshold: stop if we reach it
      node = btree_FindNode (ent_data (eo), RLAB_NAME_FFNN_ERR_THRESHOLD);
      if (node != 0)
      {
        ddummy = (int) class_double ( var_ent (node) );
        if (ddummy > 0)
          mse_threshold = ddummy;
      }

      // delta mse threshold: stop if we reach it
      node = btree_FindNode (ent_data (eo), RLAB_NAME_FFNN_ERR_DELTA);
      if (node != 0)
      {
        ddummy = class_double ( var_ent (node) );
        if (ddummy > 0)
          delta_mse_threshold = ddummy;
      }

      // how often compute mse
      node = btree_FindNode (ent_data (eo), RLAB_NAME_FFNN_NSTEP_ERR_COMP);
      if (node != 0)
      {
        idummy = (int) class_double ( var_ent (node) );
        if (idummy > 0)
          ncompmse = idummy;
      }

      // chose random pairs or go sequentially
      //  0 - iterate over entire dataset
      //  1 - iterate over random pairs
      node = btree_FindNode (ent_data (eo), RLAB_NAME_FFNN_SAMANN_RAND_PAIRS);
      if (node != 0)
      {
        ddummy = class_double ( var_ent (node) );
        if (ddummy > 0)
        {
          options |=  RLAB_SAMANN_RANDOM_PAIR;
          options &= ~RLAB_SAMANN_ALL_PAIRS;
        }
        else
        {
          options |=  RLAB_SAMANN_ALL_PAIRS;
          options &= ~RLAB_SAMANN_RANDOM_PAIR;
        }
      }

      // maximum number of iterations
      node = btree_FindNode (ent_data (eo), RLAB_NAME_GEN_MAXITER);
      if (node != 0)
      {
        idummy = (int) class_double ( var_ent (node) );
        if (idummy > 0)
          maxi = idummy;
      }

      // standard output
      node = btree_FindNode (ent_data (eo), RLAB_NAME_GEN_STDOUT);
      if (node != 0)
      {
        if (ent_type(var_ent (node))== MATRIX_DENSE_STRING)
        {
          outstream = class_char_pointer(var_ent (node));
          lout      = isvalidstring( outstream );
        }
      }

      // do random initialization?
      node = btree_FindNode (ent_data (eo), RLAB_NAME_FFNN_INIT_RAND);
      if (node != 0)
      {
        idummy = (int) class_double ( var_ent (node) );
        if (idummy>0)
          options |= RLAB_FFNN_INIT_RAND;
      }

//       // GSL simplex parameters: simplex size 'ss'
//       node = btree_FindNode (ent_data (eo), RLAB_NAME_FFNN_SAMANN_SIMPLEX_S0);
//       if (node != 0)
//       {
//         ddummy = class_double ( var_ent (node) );
//         if (ddummy > 0)
//           s0 = ddummy;
//       }

//       // GSL simplex parameters: initial simplex size 'ss'
//       node = btree_FindNode (ent_data (eo), RLAB_NAME_FFNN_SAMANN_SIMPLEX_S0);
//       if (node != 0)
//       {
//         ddummy = class_double ( var_ent (node) );
//         if (ddummy > 0)
//           s0 = ddummy;
//       }

      // GSL simplex parameters: eabs
//       node = btree_FindNode (ent_data (eo), RLAB_NAME_FFNN_SAMANN_SIMPLEX_EABS);
//       if (node != 0)
//       {
//         ddummy = class_double ( var_ent (node) );
//         if (ddummy > 0)
//           abserr = ddummy;
//       }

      // GSL simplex parameters: batch length over which sammon's mse is minimized
//       node = btree_FindNode (ent_data (eo), RLAB_NAME_FFNN_SAMANN_NTRIES);
//       if (node != 0)
//       {
//         ddummy = class_double ( var_ent (node) );
//         if (ddummy > 0)
//           num_attempts = ddummy;
//       }

      // NEWUOA: rhobeg
      node = btree_FindNode (ent_data (eo), RLAB_NAME_NUO_RHOBEG);
      if (node != 0)
      {
        ddummy = class_double ( var_ent (node) );
        if (ddummy > 0)
          rhobeg = ddummy;
      }

      node = btree_FindNode (ent_data (eo), RLAB_NAME_NUO_RHOEND);
      if (node != 0)
      {
        ddummy = class_double ( var_ent (node) );
        if (ddummy > 0 && rhoend < rhobeg)
          rhoend = ddummy;
        else
          rhoend = 0.1 * rhobeg;
      }
    }
    ent_Clean (eo);
  }

  if (lout > 1)
  {
    fptr = fopen (outstream, "a");
  }

  // transpose data and feature
  mdr_Transpose_inplace (dset);
  mdr_Transpose_inplace (tset);

  double mse=1e99;

  if (options & RLAB_FFNN_INIT_RAND)
  {
    // do random init of the network weights and biases if requested
    ffnn_randinit_wgt_theta(NetVect, NetPtr_wgt, NetPtr_ifwgt, NetPtr_theta, NetPtr_iftheta);
  }

  MDR *y_out = mdr_Create(out_dim,MNC(dset));

  //
  // main backpropagation iteration loop
  //
  switch (isolver)
  {
    default:
      // backpropagation
      ffnn_learn_backprop (dset, tset, NetVect,
                           NetPtr_wgt, NetPtr_ifwgt, NetPtr_theta, NetPtr_iftheta,
                           NetPtr_act, NetPtr_transf, lrate, lmomentum, maxi, options,
                           &mse, fptr, y_out, ncompmse, delta_mse_threshold, mse_threshold);
  }
  mdr_Transpose_inplace (y_out);
  mdr_Transpose_inplace (dset);
  mdr_Transpose_inplace (tset);

  //
  // write everything to the output list
  //
  Btree * rtree = btree_Create ();
  bw       = btree_Create ();
  if (NetPtr_ifwgt)
    bifwgt   = btree_Create ();
  bt       = btree_Create ();
  if (NetPtr_iftheta)
    biftheta = btree_Create ();
  bact     = btree_Create ();
  btransf  = btree_Create ();

  for (i=0; i<L-1; i++)
  {
    char nnode[4];

    sprintf(nnode,"%i",i+1);
    install (bw, nnode, ent_Assign_Rlab_MDR (NetPtr_wgt[i]));
    if (NetPtr_ifwgt)
      if (NetPtr_ifwgt[i])
        install (bifwgt, nnode, ent_Assign_Rlab_MDR (NetPtr_ifwgt[i]));

    install (bt, nnode, ent_Assign_Rlab_MDR (NetPtr_theta[i]));
    if (NetPtr_iftheta)
      install (biftheta, nnode, ent_Assign_Rlab_MDR (NetPtr_iftheta[i]));

    install (bact,    nnode, ent_Assign_Rlab_MDS (mds_Copy(NetPtr_act[i])));
    install (btransf, nnode, ent_Assign_Rlab_MDS (mds_Copy(NetPtr_transf[i])));
  }

  install(rtree, RLAB_NAME_FFNN_WEIGHTS, ent_Assign_Rlab_BTREE(bw));
  if (NetPtr_ifwgt)
    install(rtree, RLAB_NAME_FFNN_FIX_WEIGHTS, ent_Assign_Rlab_BTREE(bifwgt));

  install(rtree, RLAB_NAME_FFNN_BIAS, ent_Assign_Rlab_BTREE(bt));
  if (NetPtr_iftheta)
    install(rtree, RLAB_NAME_FFNN_FIX_BIAS, ent_Assign_Rlab_BTREE(biftheta));

  install(rtree, RLAB_NAME_FFNN_ACTIVATION, ent_Assign_Rlab_BTREE(bact));
  install(rtree, RLAB_NAME_FFNN_TRANSFER, ent_Assign_Rlab_BTREE(btransf));

  // install nodes array
  install(rtree, RLAB_NAME_FFNN_NODES, ent_Assign_Rlab_MDR(NetVect));

  // did we rescale the sets?
  if (isc)
  {
    // install samann scale factors
    install(rtree, RLAB_NAME_FFNN_SAMANN_SCALE, ent_Create_Rlab_Double(scale));
    // install samann scale mean
    install(rtree, RLAB_NAME_FFNN_SAMANN_MEAN, ent_Assign_Rlab_MDR(mu));
  }

  // install mse factors
  install(rtree, RLAB_NAME_FFNN_MSE, ent_Create_Rlab_Double(mse));
  // install final output of the network
  install(rtree, RLAB_NAME_FFNN_SAMANN_Y_OUT, ent_Assign_Rlab_MDR(y_out));

  // clean-up
  ent_Clean (e1);
  ent_Clean (e2);
  ent_Clean (e3);

  for (i=0; i<L-1; i++)
  {
    NetPtr_wgt[i] = 0;
    if (NetPtr_ifwgt)
      NetPtr_ifwgt[i]=0;
    NetPtr_theta[i]=0;
    if (NetPtr_iftheta)
      NetPtr_iftheta[i]=0;
    NetPtr_act[i]=0;
    NetPtr_transf[i]=0;
  }
  GC_FREE(NetPtr_wgt);
  if (NetPtr_ifwgt)
    GC_FREE(NetPtr_ifwgt);
  GC_FREE(NetPtr_theta);
  if (NetPtr_iftheta)
    GC_FREE(NetPtr_iftheta);
  GC_FREE(NetPtr_act);
  GC_FREE(NetPtr_transf);

  //
  // write result to output list
  //
  rent = ent_Create ();
  ent_data (rent) = rtree;
  ent_type (rent) = BTREE;
  return rent;
}

