/*
 * mdcf4.c Matrix Dense Complex Functions
 * Contains most of the computational interfaces to:
 * 1. expokit.f
 */

/*  This file is a part of RLaB ("Our"-LaB)
   Copyright (C) 1995  Ian R. Searle, (C) 2005 M. Kostrun

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
#include "mdc.h"
#include "mdcf1.h"
#include "util.h"
#include "complex.h"
#include "mathl.h"
#include "mdr_mdc.h"
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

int
mdc_mpow1 (MDC * A0, int j, MDC * B)
{
  int i, n;
  MDC *A=0, *C=0;

  int  m = A0->nrow;

  int k=(A0->nrow)*(A0->ncol);
  unsigned int uj, u1=1;
  double len2;

  // blas matrix multiplication
  char ta='N', tb= 'N';
  double alpha[2]={1.0,0.0}, beta[2]={0.0,0.0};

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
        Mdc0(B,i,n) = 0.0;
        Mdc0(B,n,i) = 0.0;
      }
      Mdc0(B,i,i) = 1.0;
    }
    if (j==0)
      return 1.0;
  }


  //
  // j>1
  //
  A = mdc_Copy(A0);
  C = mdc_Copy(A0);

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
      // C = B * A
      ZGEMM ( &ta, &tb, &m, &m, &m, (Complex *)alpha, MDCPTR(B), &m,
              MDCPTR(A), &m, (Complex *)beta, MDCPTR(C), &m);
      //  B = C
      for (n=0; n<k; n++)
      {
        MdcV0(B,n) = MdcV0(C,n);
      }
    }

    // now shift uj one bit to the right
    uj = uj>>1;

    // prepare A for the next round:
    // A = A * A;
    ZGEMM ( &ta, &tb, &m, &m, &m, (Complex *)alpha, MDCPTR(A), &m,
            MDCPTR(A), &m, (Complex *)beta, MDCPTR(C), &m);
    for (n=0; n<k; n++)
    {
      MdcV0(A,n) = MdcV0(C,n);
    }
  }

  // last multiplication as len2 bit appears in j
  // B = B * A
  ZGEMM ( &ta, &tb, &m, &m, &m, (Complex *)alpha, MDCPTR(B), &m,
          MDCPTR(A), &m, (Complex *)beta, MDCPTR(C), &m);
  for (n=0; n<k; n++)
  {
    MdcV0(B,n) = MdcV0(C,n);
  }

  // cleanup
  mdc_Destroy (A);
  mdc_Destroy (C);

  return 1.0;
}


static MDC * mexp_a;
static int matvec (Complex *x, Complex *y);

MDC *
mdc_mexp1 (MDC * A, double t, int id)
{
  int ideg = id;
  int m    = MNC (A);
  int lwsp = 4*m*m+ideg+1;
  int ns, iflag, iexph=0;
  int i,j;

  MDR * ipiv=0;
  MDC * wsp=0, *rmdr=0;

  wsp  = mdc_Create(1,lwsp);
  ipiv = mdi_Create(1,m);

  rmdr = mdc_Create (m,m);

  if (mdc_IsSymmetric(A))
  {
    ZHPADM (&ideg, &m, &t, MDCPTR(A), &m, MDCPTR(wsp), &lwsp, MDIPTR(ipiv),
            &iexph, &ns, &iflag);
  }
  else
  {
    ZGPADM (&ideg, &m, &t, MDCPTR(A), &m, MDCPTR(wsp), &lwsp, MDIPTR(ipiv),
            &iexph, &ns, &iflag);
  }

  for (i=0; i<m; i++)
    for (j=0; j<m; j++)
  {
      Mdc0(rmdr, i, j) = MdcV0(wsp,iexph+m*j+i-1);
  }

  mdc_Destroy (wsp);
  mdr_Destroy (ipiv);

  return rmdr;
}

MDC *
mdc_mdc_mexp2 (MDC * A, double t, MDC * v, int imet)
{
  int iflag;

  MDR * iwsp=0;
  MDC * wsp=0, * rmdr=0;

  if (MNC(A) != MNR (v))
    rerror("mexp: dimension mismatch between matrix and vector");


  if (imet == -1)
  {
    rmdr = mdc_Copy (v);

    // chebyshev method
    int m    = MNC (A);
    int lwsp = 2*m*(m+2);

    wsp  = mdc_Create(1,lwsp);
    iwsp = mdi_Create(1,m);

    ZGCHBV (&m, &t, MDCPTR(A), &m, MDCPTR(rmdr), MDCPTR(wsp), MDIPTR(iwsp), &iflag);
    // done with chebyshev
  }
  else
  {
    mexp_a = A;
    int n      = MNC(A);
    int m      = MIN(imet,n-1);
    int itrace = 0;

    rmdr = mdc_Create(n,1);
    mdc_Zero(rmdr);

    int lwsp  = n*(m+1)+n+(m+2)*(m+2)+4*(m+2)*(m+2)+7;
    int liwsp = m+2;
    wsp  = mdc_Create(1,lwsp);
    iwsp = mdi_Create(1,m+2);

    double anorm = 1.0;
    double tol   = 0.0;

    if (mdc_IsSymmetric(A))
    {
      ZHEXPV (&n, &m, &t, MDCPTR(v), MDCPTR(rmdr), &tol, &anorm,
               MDCPTR(wsp), &lwsp, MDIPTR(iwsp), &liwsp,
               matvec, &itrace, &iflag);
    }
    else
    {
      ZGEXPV (&n, &m, &t, MDCPTR(v), MDCPTR(rmdr), &tol, &anorm,
               MDCPTR(wsp), &lwsp, MDIPTR(iwsp), &liwsp,
               matvec, &itrace, &iflag);
    }

    // done with krylov
  }

  mdc_Destroy (wsp);
  mdr_Destroy (iwsp);

  return rmdr;
}


MDC *
mdc_mdr_mexp2 (MDC * A, double t, MDR * v, int imet)
{
  int iflag;

  MDR * iwsp=0;
  MDC * wsp=0, * rmdr=0, * zv = 0;

  if (MNC(A) != MNR (v))
    rerror("mexp: dimension mismatch between matrix and vector");


  if (imet == -1)
  {
    rmdr = (MDC *) mdr_coerce_mdc (v);

    // chebyshev method
    int m    = MNC (A);
    int lwsp = 2*m*(m+2);

    wsp  = mdc_Create(1,lwsp);
    iwsp = mdi_Create(1,m);

    ZGCHBV (&m, &t, MDCPTR(A), &m, MDCPTR(rmdr), MDCPTR(wsp), MDIPTR(iwsp), &iflag);
    // done with chebyshev
  }
  else
  {
    zv = (MDC *) mdr_coerce_mdc (v);

    mexp_a = A;
    int n      = MNC (A);
    int m      = MIN(imet,n-1);
    int itrace = 0;

    rmdr = mdc_Create(n,1);
    mdc_Zero(rmdr);

    int lwsp  = n*(m+1)+n+(m+2)*(m+2)+4*(m+2)*(m+2)+7;
    int liwsp = m+2;
    wsp  = mdc_Create(1,lwsp);
    iwsp = mdi_Create(1,m+2);

    double anorm = 1.0;
    double tol   = 0.0;

    if (mdc_IsSymmetric(A))
    {
      ZHEXPV (&n, &m, &t, MDCPTR(zv), MDCPTR(rmdr), &tol, &anorm,
               MDCPTR(wsp), &lwsp, MDIPTR(iwsp), &liwsp,
               matvec, &itrace, &iflag);
    }
    else
    {
      ZGEXPV (&n, &m, &t, MDCPTR(zv), MDCPTR(rmdr), &tol, &anorm,
               MDCPTR(wsp), &lwsp, MDIPTR(iwsp), &liwsp,
               matvec, &itrace, &iflag);
    }

    mdc_Destroy (zv);
    // done with krylov
  }

  mdc_Destroy (wsp);
  mdr_Destroy (iwsp);

  return rmdr;
}

MDC *
mdr_mdc_mexp2 (MDR * A, double t, MDC * v, int imet)
{
  int iflag;

  MDR * iwsp=0;
  MDC * wsp=0, * rmdr=0, * za = 0;

  za = (MDC *) mdr_coerce_mdc (A);

  if (MNC(A) != MNR (v))
    rerror("mexp: dimension mismatch between matrix and vector");


  if (imet == -1)
  {
    rmdr = mdc_Copy (v);

    // chebyshev method
    int m    = MNC (za);
    int lwsp = 2*m*(m+2);

    wsp  = mdc_Create(1,lwsp);
    iwsp = mdi_Create(1,m);

    ZGCHBV (&m, &t, MDCPTR(za), &m, MDCPTR(rmdr), MDCPTR(wsp), MDIPTR(iwsp), &iflag);
    // done with chebyshev
  }
  else
  {
    mexp_a = za;
    int n      = MNC (za);
    int m      = MIN(imet,n-1);
    int itrace = 0;

    rmdr = mdc_Create(n,1);
    mdc_Zero(rmdr);

    int lwsp  = n*(m+1)+n+(m+2)*(m+2)+4*(m+2)*(m+2)+7;
    int liwsp = m+2;
    wsp  = mdc_Create(1,lwsp);
    iwsp = mdi_Create(1,m+2);

    double anorm = 1.0;
    double tol   = 0.0;

    if (mdc_IsSymmetric(za))
    {
      ZHEXPV (&n, &m, &t, MDCPTR(v), MDCPTR(rmdr), &tol, &anorm,
               MDCPTR(wsp), &lwsp, MDIPTR(iwsp), &liwsp,
               matvec, &itrace, &iflag);
    }
    else
    {
      ZGEXPV (&n, &m, &t, MDCPTR(v), MDCPTR(rmdr), &tol, &anorm,
               MDCPTR(wsp), &lwsp, MDIPTR(iwsp), &liwsp,
               matvec, &itrace, &iflag);
    }

    // done with krylov
  }

  mdc_Destroy (wsp);
  mdc_Destroy (za);
  mdr_Destroy (iwsp);

  return rmdr;
}

//
// local matrix vector multiplication
//
static int
matvec (Complex *x, Complex *y)
{
  int i, j;
  int n = mexp_a->nrow;

  for (i = 0; i < n; i++)
  {
    y[i] = 0;
    for (j = 0; j < n; j++)
    {
      y[i] += Mdc0(mexp_a, i, j) * x[j];
    }
  }

  return 0;
}

