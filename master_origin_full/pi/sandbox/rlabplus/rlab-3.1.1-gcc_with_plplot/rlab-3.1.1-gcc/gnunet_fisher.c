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

//
//  gnu_fisher.c
//    Calculate the Fisher (linear) discrimint of a two class dataset.
//
//  References:
//    1. sprannlib's stat/fisher.c by W.F. Schmidt
//    2. Duin R.P.W, "Scheiding van Puntverzamelingen", Afstudeer verslag,
//    Technische Hogeschool Delft, Juni 1970. Appendix A.
//
#include "lp.h"
static int weight_vector (MDR *Ka, MDR *Kb, MDR *Ma, MDR *Mb, int n, MDR *W)
{
  int ival=RLAB_STATUS_FAILURE;
  int i, j;

  MDR *K = mdr_Create(n,n);
  if (!K)
    goto _exit_weight_vector;

  for (i=0; i<n; i++)
  {
    for (j=0; j<n; j++)
      Mdr0(K,i,j) = 0.5 * (Mdr0(Ka,i,j) + Mdr0(Kb,i,j));
  }

  // find inverse of K and store it in K
  ival = mdr_matinvd(K, NULL);

  if (ival)
    goto _exit_weight_vector;

  for (i=0; i<n; i++)
  {
    MdrV0(W,i) = 0.0;
    for (j=0; j<n; j++)
      MdrV0(W,i) += Mdr0(K,i,j) * (MdrV0(Ma,j) - MdrV0(Mb,j));
  }

_exit_weight_vector:

  if (K)
    mdr_Destroy (K);

  return ival;
}


static int fisher (MDR *Ka, MDR *Kb, MDR *Ma, MDR *Mb, double c, int n, MDR *W, double *b)
{
  int i, j;
  double WtKaW, WtKbW, WMa, WMb;

  /* Calculate the normal vector W. */

  if (weight_vector(Ka, Kb, Ma, Mb, n, W))
    return RLAB_STATUS_FAILURE;

  /* Calculate the constant term B. */
  WtKaW = WtKbW = WMa = WMb = 0.0;
  for (i=0; i<n; i++)
  {
    WMa += MdrV0(W,i) * MdrV0(Ma,i);
    WMb += MdrV0(W,i) * MdrV0(Mb,i);
    for (j=0; j<n; j++)
    {
      WtKaW += MdrV0(W,i)* Mdr0(Ka,i,j) * MdrV0(W,j);
      WtKbW += MdrV0(W,i)* Mdr0(Kb,i,j) * MdrV0(W,j);
    }
  }

  WtKaW = sqrt(WtKaW);
  WtKbW = sqrt(WtKbW);

  *b = WMb * WtKaW + WMa * WtKbW;

  if (WtKaW + WtKbW < 1e-6)
    return RLAB_STATUS_FAILURE;

  *b /= -1.0 * (WtKaW + WtKbW);
  *b -= log(c);

  return RLAB_STATUS_SUCCESS;
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
static int getcluster_mean_var(MDR * dataset, MDR * clusterid, int cluster0, MDR ** mean, MDR ** var)
{
  int j, k, l;
  int ndata, nattr, ival=RLAB_STATUS_SUCCESS;

  MDR *cmask=0;
  MDR *vmask=0;

  nattr = MNC(dataset);
  ndata = MNR(dataset);

  cmask = (MDR*) mdi_Create (1, nattr);
  vmask = (MDR*) mdi_Create (nattr, nattr);

  for (j=0; j<nattr; j++)
  {
    MdrV0(*mean,j)  = create_nan();
    MdiV0( cmask,j) = 0;
    for (k=0; k<nattr; k++)
    {
      Mdi0(vmask,j,k) = 0;
      Mdr0(*var,j,k)  = create_nan();
    }
  }

  for (l=0; l<ndata; l++)
  {
    if (cluster0 != MdrV0(clusterid,l))
      continue;

    // find mean
    for (j=0; j<nattr; j++)
    {
      if (isnand(Mdr0(dataset,l,j)))
        continue;

      if (isnand(MdrV0(*mean,j)))
      {
        MdrV0(*mean,j)  = Mdr0(dataset,l,j);
        MdiV0( cmask,j) = 1;
      }
      else
      {
        MdrV0(*mean,j) += Mdr0(dataset,l,j);
        MdiV0( cmask,j) += 1;
      }

      // find mean
      for (k=0; k<=j; k++)
      {
        if (isnand(Mdr0(dataset,l,k)))
          continue;

        if (isnand(Mdr0(*var,j,k)))
        {
          Mdr0(*var,j,k) = Mdr0(dataset,l,j) * Mdr0(dataset,l,k);
          Mdi0( vmask,j,k) = 1;
        }
        else
        {
          Mdr0(*var,j,k) += Mdr0(dataset,l,j) * Mdr0(dataset,l,k);
          Mdi0( vmask,j,k) += 1;
        }
      } // for (k=0; k<=j; k++)
    } // for (j=0; j<nattr; j++)
  }

  // mean_j = E[x_j]
  // var_j,k  = E[x_j * x_k] - E[x_j] * E[x_k] for k<=j
  for (j=0; j<nattr; j++)
  {
    if (MdiV0(cmask,j))
      MdrV0(*mean,j) /= (double) MdiV0(cmask,j);
    else
      ival = RLAB_STATUS_FAILURE;

    for (k=0; k<=j; k++)
    {
      if (Mdi0(vmask,j,k))
      {
        Mdr0(*var,j,k) /= ((double) Mdi0(vmask,j,k));
        Mdr0(*var,j,k) -= MdrV0(*mean,j) * MdrV0(*mean,k);
        if (j!=k)
          Mdr0(*var,k,j) = Mdr0(*var,j,k);
      }
      else
        ival = RLAB_STATUS_FAILURE;
    }
  }

  if (cmask)
    mdr_Destroy(cmask);
  if (vmask)
    mdr_Destroy(vmask);

  return ival;
}

int fisher_dataset (MDR *dset, MDR *cset, int Class0, int Class1, double c, MDR *W, double *b)
{
  int ival=RLAB_STATUS_FAILURE;
  int nattr = MNC(dset);

  //
  // class 0: mean and covariance
  //
  MDR *M0 = mdr_Create(1,nattr);
  if (!M0)
    goto _exit_fisher_dataset;
  MDR *C0 = mdr_Create(nattr,nattr);
  if (!C0)
    goto _exit_fisher_dataset;
  if (getcluster_mean_var(dset, cset, Class0, &M0, &C0))
    goto _exit_fisher_dataset;

  //
  // class 1: mean and covariance
  //
  MDR *M1 = mdr_Create(1,nattr);
  if (!M1)
    goto _exit_fisher_dataset;
  MDR *C1 = mdr_Create(nattr,nattr);
  if (!C1)
    goto _exit_fisher_dataset;
  if (getcluster_mean_var(dset, cset, Class1, &M1, &C1))
    goto _exit_fisher_dataset;

  //
  // do the fisher thingy
  //
  ival=fisher(C0, C1, M0, M1, c, nattr, W, b);

_exit_fisher_dataset:

  if (M0)
    mdr_Destroy(M0);
  if (C0)
    mdr_Destroy(C0);
  if (M1)
    mdr_Destroy(M1);
  if (C1)
    mdr_Destroy(C1);

  return ival;
}

