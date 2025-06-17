/*
 * mdrf4.c
 * Matrix-Dense-Real Fortran Interfaces
 * Contains most of the computational interfaces to:
 * 1. expokit.f
 * 2. mpow
 */

/*  This file is a part of RLaB ("Our"-LaB)
   Copyright (C) 1995  Ian R. Searle, (C) 2005-2008 M. Kostrun

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
#include "ent.h"
#include "btree.h"
#include "symbol.h"
#include "mdr.h"
#include "mdrf1.h"
#include "mdc.h"
#include "util.h"
#include "mathl.h"

#include "fi.h"
#include "lp.h"
#include "blas.h"

#include <stdio.h>
#include <math.h>

/*
 * Create Infs, and NaNs
 */
#include "mathl.h"
#include "expokit.h"

//
// matrix exponentiation: A^j
//

int
mdr_mpow1 (MDR * A0, int j, MDR * B)
{
  int i, n;
  MDR *A=0, *C=0;

  int  m = A0->nrow;

  int k=(A0->nrow)*(A0->ncol);
  unsigned int uj, u1=1;
  double len2;

  // blas matrix multiplication
  char   ta='N', tb='N';
  double alpha=1.0, beta=0.0;

  if (j == 1 || (A0->nrow != A0->ncol))
    return 1.0;

  uj = (unsigned long) j;
  if (!(uj & u1))
  {
    // j is even: B=1
    for (i=0; i<m; i++)
    {
      for (n=0; n<i; n++)
      {
        Mdr0(B,i,n) = 0.0;
        Mdr0(B,n,i) = 0.0;
      }
      Mdr0(B,i,i) = 1.0;
    }
    if (j==0)
      return 1.0;
  }

  //
  // j > 1: need to do some calculation
  //
  A = mdr_Float_BF(A0);
  C = mdr_Float_BF(A0);

  // max no of operations
  len2 = -1;
  uj = (unsigned long) j;
  while (uj)
  {
    uj = uj>>1;
    len2++;
  }
  uj = (unsigned long) j;

  for (i=0; i<len2; i++)
  {
    // does i-th bit appear in binary representation of j
    // but skip the zeroth bit, as we have already done these manipulations
    if ( (uj & u1)  && i)
    {

      //
      // B -> B * A
      //
      // perform multiplication in two steps
      //  C = B * A
      DGEMM ( &ta, &tb, &m, &m, &m, &alpha, MDRPTR(B), &m,
               MDRPTR(A), &m, &beta, MDRPTR(C), &m);
      //  B = C
      for (n=0; n<k; n++)
        MdrV0(B,n) = MdrV0(C,n);
    } // if ( uj & u1 )

    // now shift uj one bit to the right
    uj = uj>>1;

    // prepare A for the next round:
    // A = A * A;
    DGEMM ( &ta, &tb, &m, &m, &m, &alpha, MDRPTR(A), &m,
             MDRPTR(A), &m, &beta, MDRPTR(C), &m);
    for (n=0; n<k; n++)
      MdrV0(A,n) = MdrV0(C,n);
  }

  // last multiplication as len2 bit appears in j
  // B = B * A
  DGEMM ( &ta, &tb, &m, &m, &m, &alpha, MDRPTR(B), &m,
           MDRPTR(A), &m, &beta, MDRPTR(C), &m);
  for (n=0; n<k; n++)
    MdrV0(B,n) = MdrV0(C,n);

  // cleanup
  mdr_Destroy (A);
  mdr_Destroy (C);

  return 1.0;
}

//
// matrix exponantiation: exp(A)
//

static MDR * mexp_a;
static int matvec();

MDR *
mdr_mexp1 (MDR * A, double t, int id)
{
  int ideg = id;
  int m    = MNC (A);
  int lwsp = 4*m*m + ideg + 1;
  int ns, iflag, iexph=0;
  int i,j;

  MDR * wsp=0, *ipiv=0, *rmdr=0;

  wsp  = mdr_Create(1,lwsp);
  ipiv = mdi_Create(1,m);

  rmdr = mdr_Create (m,m);

  if (mdr_IsSymmetric(A))
  {
    DSPADM(&ideg, &m, &t, MDRPTR(A), &m, MDRPTR(wsp), &lwsp, MDIPTR(ipiv),
           &iexph, &ns, &iflag);
  }
  else
  {
    DGPADM (&ideg, &m, &t, MDRPTR(A), &m, MDRPTR(wsp), &lwsp, MDIPTR(ipiv),
            &iexph, &ns, &iflag);
  }

  for (i=0; i<m; i++)
    for (j=0; j<m; j++)
      Mdr0(rmdr, i, j) = MdrV0(wsp,iexph+m*j+i-1);

  mdr_Destroy (wsp);
  mdr_Destroy (ipiv);

  return rmdr;
}

MDR *
mdr_mdr_mexp2 (MDR * A, double t, MDR * v, int imet)
{
  int iflag;

  MDR *wsp=0, *iwsp=0, *rmdr=0;

  if (MNC(A) != MNR (v))
    rerror("mexp: dimension mismatch between matrix and vector");


  if (imet == -1)
  {
    rmdr = mdr_Copy (v);

    // chebyshev method
    int m    = MNC (A);
    int lwsp = 2*m*(m+2);

    wsp  = mdr_Create(1,lwsp);
    iwsp = mdi_Create(1,m);

    if (mdr_IsSymmetric(A))
    {
      DSCHBV (&m, &t, MDRPTR(A), &m, MDRPTR(rmdr), MDRPTR(wsp), MDIPTR(iwsp), &iflag);
    }
    else
    {
      DGCHBV (&m, &t, MDRPTR(A), &m, MDRPTR(rmdr), MDRPTR(wsp), MDIPTR(iwsp), &iflag);
    }
    // done with chebyshev
  }
  else
  {
    mexp_a = A;
    int n      = MNC (A);
    int m      = MIN(imet,n-1);
    int itrace = 0;

    rmdr = mdr_Create(n,1);
    mdr_Zero(rmdr);

    int lwsp  = n*(m+1)+n+(m+2)*(m+2)+4*(m+2)*(m+2)+7;
    int liwsp = m+2;
    wsp  = mdr_Create(1,lwsp);
    iwsp = mdi_Create(1,m+2);

    double anorm = 1.0;
    double tol   = 0.0;

    if (mdr_IsSymmetric(A))
    {
      DSEXPV (&n, &m, &t, MDRPTR(v), MDRPTR(rmdr), &tol, &anorm,
               MDRPTR(wsp), &lwsp, MDIPTR(iwsp), &liwsp,
               matvec, &itrace, &iflag);
    }
    else
    {
      DGEXPV (&n, &m, &t, MDRPTR(v), MDRPTR(rmdr), &tol, &anorm,
              MDRPTR(wsp), &lwsp, MDIPTR(iwsp), &liwsp,
              matvec, &itrace, &iflag);
    }

    // done with krylov
  }

  mdr_Destroy (wsp);
  mdr_Destroy (iwsp);

  return rmdr;
}

//
// local matrix vector multiplication
//
static int
matvec (double *x, double *y)
{
  int i, j;
  int n = mexp_a->nrow;

  for (i = 0; i < n; i++)
  {
    y[i] = 0;
    for (j = 0; j < n; j++)
    {
      y[i] += Mdr0(mexp_a, i, j) * x[j];
    }
  }

  return 0;
}


