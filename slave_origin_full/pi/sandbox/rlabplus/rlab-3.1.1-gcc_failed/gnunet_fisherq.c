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
//  gnunet_quadr.c
//    Calculate the Fisher (quadratic) discrimint of a two class dataset.
//
//  References:
//    1. sprannlib's stat/quadr.c by W.F. Schmidt
//    2. Duin R.P.W, "Scheiding van Puntverzamelingen", Afstudeer verslag,
//    Technische Hogeschool Delft, Juni 1970. Appendix A.
//
static double weight0_vector (MDR *Kainv, MDR *Kbinv, MDR *Ma, MDR *Mb, double c, int n, double deta, double detb)
{
  int i;
  double b=0.0;

  MDR *bh = mdr_Create(1,n);

  mdr_rgemv(Kbinv, 1, Mb, 0, &bh);
  for (i=0; i<n; i++)
    b += MdrV0(bh,i) * MdrV0(Mb,i);

  mdr_rgemv(Kainv, 1, Ma, 0, &bh);
  for (i=0; i<n; i++)
    b -= MdrV0(bh,i) * MdrV0(Ma,i);

  mdr_Destroy(bh);

  // we sustained
  b += 2 * log(c) + 2 * log(deta / detb);

  return b;
}

static int fisherq (MDR *Ka, MDR *Kb, MDR *Ma, MDR *Mb, double c, int n, MDR *W2, MDR *W1, double *b)
{
  double deta, detb;
  int i, j, ival=RLAB_STATUS_SUCCESS;

  MDR * Kainv = mdr_Float_BF(Ka);
  MDR * Kbinv = mdr_Float_BF(Kb);

  mdr_matinvd(Kainv, &deta);
  mdr_matinvd(Kbinv, &detb);

 /*
  * Calculate the quadratic matrix W2 and calculate the normal vector
  * W1.
  */
  for (i=0; i<n; i++)
  {
    MdrV0(W1,i) = 0.0;
    for (j=0; j<n; j++)
    {
      Mdr0(W2,i,j) = Mdr0(Kbinv,i,j) - Mdr0(Kainv,i,j);
      MdrV0(W1,i) += 2 * (Mdr0(Kainv,i,j) * MdrV0(Ma,j) - Mdr0(Kbinv,i,j) * MdrV0(Mb,j));
    }
  }

  /* Calculate the constant term b. */
  *b = weight0_vector(Kainv, Kbinv, Ma, Mb, c, n, deta, detb);
  if (isnand(*b))
    ival=RLAB_STATUS_FAILURE;

  // we sustained
  mdr_Destroy(Kainv);
  mdr_Destroy(Kbinv);

  return ival;
}

int fisherq_dataset (MDR *dset, MDR *cset, int Class0, int Class1, double c, MDR *W2, MDR *W1, double *b)
{
  int nattr = MNC(dset);
  int ival = RLAB_STATUS_FAILURE;

  //
  // class 0: mean and covariance
  //

  MDR *M0 = mdr_Create(1,nattr);
  if (!M0)
    goto _exit_fisherq_dataset;
  MDR *C0 = mdr_Create(nattr,nattr);
  if (!C0)
    goto _exit_fisherq_dataset;
  if (getcluster_mean_var(dset, cset, Class0, &M0, &C0))
    goto _exit_fisherq_dataset;

  //
  // class 1: mean and covariance
  //
  MDR *M1 = mdr_Create(1,nattr);
  if (!M1)
    goto _exit_fisherq_dataset;
  MDR *C1 = mdr_Create(nattr,nattr);
  if (!C1)
    goto _exit_fisherq_dataset;
  if (getcluster_mean_var(dset, cset, Class1, &M1, &C1))
    goto _exit_fisherq_dataset;

  //
  // do the fisher quadratic thingy
  //
  ival = fisherq(C0, C1, M0, M1, c, nattr, W2, W1, b);

_exit_fisherq_dataset:

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
