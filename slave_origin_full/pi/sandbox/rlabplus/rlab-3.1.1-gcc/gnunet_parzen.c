/* This file is released as a part of rlabplus project, see rlabplus.sourceforge.net
   Copyright (C) 2006 M. Kostrun

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

// ------------------------------------------------------- //
//
//
// CLASSIFIERS: PARZEN
//
//
// ------------------------------------------------------- //

static double fparzen (double s, MDR *pts, int dim)
{
  double sk, fc, f, ff, p, dis;
  int i, j;

  int size = MNR(pts);

  // use euclidean metrics between the rows of the matrix
  double (*distance)(int, MDR*, MDR*, MDR*, int, int) =  distance_metric('e');

  sk = -1.0 / (2.0 * s * s);
  fc = 0.0;

  for (i=0; i<size; i++)
  {
    f = 0.0;
    ff = 0.0;
    for (j=0; j<size; j++)
    {
      if (i != j)
      {
        dis = distance(1, pts, pts, NULL, i, j);
        p   = exp((double) (sk * dis));
        f  += p;
        ff += dis * p;
      }
    }
    fc += ff / (f * s * s);
  }
  fc -= (double) (size * dim);

  return fc;
}

// points contains sub-set of data all of requested class
static MDR * find_points (MDR *dset, MDR *classset, int class)
{
  int size=0, i, j, k;
  int nattr = MNC(dset);
  int ndata = MNR(dset);

  MDR *points=0;

  for (i=0; i<ndata; i++)
    if (mdiV0(classset,i)==class)
      size++;

  if (size)
  {
    points = mdr_Create(size,nattr);

    k=0;
    for (i=0; i<ndata; i++)
    {
      if (mdiV0(classset,i)==class)
      {
        for (j=0; j<nattr; j++)
          Mdr0(points,k,j) = mdr0(dset,i,j);
        k++;
      }
    }
  }

  return (points);
}

/*
 * PARZEN_BEST_S Calculate optimal value of parzen smoothing parameter s for
 * a certain class. It uses the following inputs:
 *
 * class : class of the objects in the learning set. init_s: initial value of
 * the smoothing parameter.
 *
 */

static int
parzen_dset_best_s (MDR * dset, MDR *classset, int ncls, int class, double init_s, double *s3, int maxit)
{
  double s1, s2, f1, f2, f3, s1k, s2k, rmn;
  int stopcrit=0, cycles=0;

  MDR *points=find_points(dset, classset, class);

  int size = MNR(points);

  rmn = size * ncls;
  s1 = init_s;
  f1 = fparzen(s1,points,ncls);
  s2 = s1 * sqrt(f1 / rmn + 1.0);
  f2 = fparzen(s2,points,ncls);

  do
  {
    s1k = s1 * s1;
    s2k = s2 * s2;
    *s3 = s1k * s2k * (f2 - f1) / (f2 * s2k - f1 * s1k);
    *s3 = *s3 < 0 ? sqrt(f2 * s2k / rmn + s2k) : sqrt(*s3);
    f3 = fparzen(*s3,points,ncls);
    if ((fabs(1 - f3 / f2) < 1E-4) ||
         (fabs(1 - *s3 / s2) < 1E-3) ||
         (fabs(f3) < 1e-16))
    {
      stopcrit=1;
    }
    else
    {
      f1 = f2;
      f2 = f3;
      s1 = s2;
      s2 = *s3;
    }
    if (cycles++ > maxit)
      stopcrit=2;
  }
  while (!stopcrit);

  mdr_Destroy(points);

  return stopcrit;
}

/*
 * PARZEN_DATASET_CLASS: Classify a sample using the  DATASET *dset, which
 * is a label set which is used as a learning set.
 * Other parameters are the smoothing parameter S[0..NCls-1] and the a priori
 * possibility P[0..NCls-1] for each class.
 *
 * Return value : The class label of the most probable class.
 */

static MDR * parzen_dataset_class (MDR *dset, MDR *classset, MDR *samp, MDR *S, MDR *P, int NCls)
{
  double maxdens, S2;
  MDR *dmax=0;
  int i, k;

  int nattr = MNC(dset);
  int ndata = MNR(dset);
  int nsamp = MNR(samp);

  // use euclidean metrics between the rows of the matrix
  double (*distance)(int, MDR*, MDR*, MDR*, int, int) =  distance_metric('e');

  int *npts = GC_malloc (NCls * sizeof(int));
  double *parz = GC_malloc (NCls * sizeof(double));

  dmax = mdr_Create(nsamp,1);

  for (k=0; k<nsamp; k++)
  {
    for (i = 0; i < NCls; i++)
    {
      parz[i] = 0.0;
      npts[i] = 0;
    }

    //
    // Calculate densities
    //
    for (i=0; i<ndata; i++)
    {
      npts[ mdiV0(classset,i) ] ++;
      S2 = MdrV0(S, mdiV0(classset,i) );
      S2 *= -2.0 * S2;
      parz[mdiV0(classset,i)] += exp( distance(1, samp, dset, NULL, k, i) / S2);
    }

    for (i=0; i<NCls; i++)
      parz[i] *= mdrV0(P,i) * pow(mdrV0(S,i), -(double) nattr) / ((double) npts[i]);

    // Determine class with maximum probability
    maxdens = 0.0;
    for (i = 0; i < NCls; i++)
      if (parz[i] > maxdens)
      {
        maxdens = parz[i];
        MdrV0(dmax,k) = i;
      }

  } // for (k=0; k<nsamp; k++)

  GC_free (npts);
  GC_free (parz);

  return dmax;
}


