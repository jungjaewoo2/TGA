/* mdrf1.c Matrix Dense Real Functions ... */

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
#include "mem.h"
#include "mdr.h"
#include "mdc.h"
#include "mds.h"
#include "btree.h"
#include "util.h"
#include "mdrf1.h"
#include "mathl.h"
#include "symbol.h"
#include "sort.h"
#include "rlab_solver_parameters_names.h"

#include <stdio.h>
#include <math.h>

int do_cmplx (MDR * m);
int do_cmplx_1 (MDR * m);

MDR * mdr_Mean(MDR *x, MDR *factor, MDR *use, MDR *ignore, int rowdominant, int igninfs, int flat)
{
  MDR *w=0, *use_r=0, *ignore_r=0, *use_idx=0, *dummy=0, *use_factors=0;
  int size_x = SIZE(x), nr_w=0, nc_w=0, stride_i, stride_j, i, i2, j, k, npts, ndim;

  if (size_x<1)
    goto _exit_mdr_mean;

  // copy two together:
  if (SIZE(use)>0)
  {
    use_r = mdr_VectorSet(use);
  }
  if (SIZE(ignore)>0)
  {
    ignore_r = mdr_VectorSet(ignore);
  }

  if (EQVECT(x) || flat)
  {
    stride_i = 0;
    stride_j = 1;
    nr_w = 1;
    nc_w = 1;
    ndim = 1;
    npts = size_x;
  }
  else
  {
    if (rowdominant)
    {
      nr_w = MNR(x);
      nc_w = 1;
      npts = MNC(x);
      ndim = MNR(x);
      stride_i = 1;
      stride_j = nr_w;
    }
    else
    {
      nr_w = 1;
      nc_w = MNC(x);
      npts = MNR(x);
      ndim = MNC(x);
      stride_i = npts;
      stride_j = 1;
    }
  }
  use_idx = mdr_Create(1, npts);
  for (i=0; i<npts; i++)
    MdrV0(use_idx,i) = i+1;

  // find indices of the data that will be used in calculation of the mean
  if (use_r)
  {
    dummy = mdr_VectorIntersect(use_idx, use_r);
    mdr_Destroy(use_idx);
    use_idx = dummy;
    dummy = 0;
  }
  if (ignore_r)
  {
    dummy = mdr_VectorComplement(use_idx, ignore_r);
    mdr_Destroy(use_idx);
    use_idx = dummy;
    dummy = 0;
  }

  // find indices of the factors (independent variables) that will be used
  // in calculation of the mean
  use_factors = mdr_Create(1, ndim);
  for (i=0; i<ndim; i++)
    MdrV0(use_factors,i) = i+1;
  if (factor)
  {
    dummy = mdr_VectorIntersect(use_factors, factor);
    mdr_Destroy(use_factors);
    use_factors = dummy;
    dummy = 0;
    ndim = SIZE(use_factors);
  }
  if (ndim<1)
    goto _exit_mdr_mean;

  // calculate mean
  w = mdr_Create(nr_w, nc_w);
  mdr_Nan(w);
  for (i2=0; i2<ndim; i2++)
  {
    int count=0;
    i = MdrV0(use_factors,i2) - 1;
    for (k=0; k<SIZE(use_idx); k++)
    {
      j = MdrV0(use_idx,k) - 1;
      double d = mdrV0(x,stride_i * i + stride_j * j);

      if (isnand(d))
        continue;

      if (igninfs)
        if (create_inf () == d)
          continue;

      if (isnand(MdrV0(w,i)))
        MdrV0(w,i) = d;
      else
        MdrV0(w,i) += d;
      count++;
    }
    MdrV0(w,i) /= (double) count;
  }

_exit_mdr_mean:

  if (use_r)
    mdr_Destroy(use_r);
  if (ignore_r)
    mdr_Destroy(ignore_r);
  if (use_idx)
    mdr_Destroy(use_idx);
  if (use_factors)
    mdr_Destroy(use_factors);

  return w;
}

MDR * mdr_Var(MDR *x, MDR *factor, MDR *use, MDR *ignore, int rowdominant, int igninfs, int unbias, MDR *m, int flat)
{
  MDR *w=0, *use_r=0, *ignore_r=0, *use_idx=0, *dummy=0, *use_factors=0;
  int size_x = SIZE(x), nr_w=0, nc_w=0, stride_i, stride_j, i, i2, j, k, npts, ndim;

  if (size_x<1)
    goto _exit_mdr_var;

  // copy two together:
  if (SIZE(use)>0)
  {
    use_r = mdr_VectorSet(use);
  }
  if (SIZE(ignore)>0)
  {
    ignore_r = mdr_VectorSet(ignore);
  }

  if (EQVECT(x) || flat)
  {
    stride_i = 0;
    stride_j = 1;
    nr_w = 1;
    nc_w = 1;
    ndim = 1;
    npts = size_x;
  }
  else
  {
    if (rowdominant)
    {
      nr_w = MNR(x);
      nc_w = 1;
      npts = MNC(x);
      ndim = MNR(x);
      stride_i = 1;
      stride_j = nr_w;
    }
    else
    {
      nr_w = 1;
      nc_w = MNC(x);
      npts = MNR(x);
      ndim = MNC(x);
      stride_i = npts;
      stride_j = 1;

    }
  }
  use_idx = mdr_Create(1, npts);
  for (i=0; i<npts; i++)
    MdrV0(use_idx,i) = i+1;

  // find indices of the data that will be used in calculation of the mean
  if (use_r)
  {
    dummy = mdr_VectorIntersect(use_idx, use_r);
    mdr_Destroy(use_idx);
    use_idx = dummy;
    dummy = 0;
  }
  if (ignore_r)
  {
    dummy = mdr_VectorComplement(use_idx, ignore_r);
    mdr_Destroy(use_idx);
    use_idx = dummy;
    dummy = 0;
  }

  // find indices of the factors (independent variables) that will be used
  // in calculation of the mean
  use_factors = mdr_Create(1, ndim);
  for (i=0; i<ndim; i++)
    MdrV0(use_factors,i) = i+1;
  if (factor)
  {
    dummy = mdr_VectorIntersect(use_factors, factor);
    mdr_Destroy(use_factors);
    use_factors = dummy;
    dummy = 0;
    ndim = SIZE(use_factors);
  }
  if (ndim<1)
    goto _exit_mdr_var;

  // calculate mean
  w = mdr_Create(nr_w, nc_w);
  mdr_Nan(w);
  for (i2=0; i2<ndim; i2++)
  {
    int count=0;
    i = MdrV0(use_factors,i2) - 1;
    double s=create_nan(),s2= create_nan();
    for (k=0; k<SIZE(use_idx); k++)
    {
      j = MdrV0(use_idx,k) - 1;
      double d = mdrV0(x,stride_i * i + stride_j * j);

      if (isnand(d))
        continue;

      if (igninfs)
        if (create_inf () == d)
          continue;

      if (isnand(s))
      {
        s  = d;
        s2 = d*d;
      }
      else
      {
        s  += d;
        s2 += d*d;
      }
      count++;
    }
    s  /= (double) count;
    s2 /= (double) (count - unbias);

    if (m)
    {
      s = mdrV1(m,MIN(SIZE(m),i+1)) * (2.0 * s - mdrV1(m,MIN(SIZE(m),i+1)));
    }
    else
    {
      s *= s;
    }

    MdrV0(w,i) = s2 - ((double) count / (double)(count - unbias)) * s;
  }

_exit_mdr_var:

  if (use_r)
    mdr_Destroy(use_r);
  if (ignore_r)
    mdr_Destroy(ignore_r);
  if (use_idx)
    mdr_Destroy(use_idx);
  if (use_factors)
    mdr_Destroy(use_factors);

  return w;
}



MDR * mdr_Covar(MDR *x, MDR *factor, MDR *use, MDR *ignore, int rowdominant, int igninfs, int unbias, MDR *m)
{
  MDR *w=0, *use_r=0, *ignore_r=0, *use_idx=0, *dummy=0, *use_factors=0;
  int size_x = SIZE(x), nr_w=0, nc_w=0, stride_i, stride_j, i, i1, i2, _i1, _i2, j, k, npts, ndim;

  if (size_x<1)
    goto _exit_mdr_var;

  // copy two together:
  if (SIZE(use)>0)
  {
    use_r = mdr_VectorSet(use);
  }
  if (SIZE(ignore)>0)
  {
    ignore_r = mdr_VectorSet(ignore);
  }

  if (EQVECT(x))
  {
    stride_i = 0;
    stride_j = 1;
    nr_w = 1;
    nc_w = 1;
    ndim = 1;
    npts = size_x;
  }
  else
  {
    if (rowdominant)
    {
      nr_w = MNR(x);
      nc_w = MNR(x);
      npts = MNC(x);
      ndim = MNR(x);
      stride_i = 1;
      stride_j = nr_w;
    }
    else
    {
      nr_w = MNC(x);
      nc_w = MNC(x);
      npts = MNR(x);
      ndim = MNC(x);
      stride_i = npts;
      stride_j = 1;

    }
  }
  use_idx = mdr_Create(1, npts);
  for (i=0; i<npts; i++)
    MdrV0(use_idx,i) = i+1;

  // find indices of the data that will be used in calculation of the mean
  if (use_r)
  {
    dummy = mdr_VectorIntersect(use_idx, use_r);
    mdr_Destroy(use_idx);
    use_idx = dummy;
    dummy = 0;
  }
  if (ignore_r)
  {
    dummy = mdr_VectorComplement(use_idx, ignore_r);
    mdr_Destroy(use_idx);
    use_idx = dummy;
    dummy = 0;
  }

  // find indices of the factors (independent variables) that will be used
  // in calculation of the mean
  use_factors = mdr_Create(1, ndim);
  for (i=0; i<ndim; i++)
    MdrV0(use_factors,i) = i+1;
  if (factor)
  {
    dummy = mdr_VectorIntersect(use_factors, factor);
    mdr_Destroy(use_factors);
    use_factors = dummy;
    dummy = 0;
    ndim = SIZE(use_factors);
  }
  if (ndim<1)
    goto _exit_mdr_var;

  // calculate mean
  w = mdr_Create(nr_w, nc_w);
  mdr_Nan(w);
  for (_i1=0; _i1<ndim; _i1++)
  {
    i1 = MdrV0(use_factors,_i1) - 1;

    for (_i2=_i1; _i2<ndim; _i2++)
    {
      i2 = MdrV0(use_factors,_i2) - 1;
      int count=0;
      double s1=create_nan(),s1_2= create_nan();
      double s2=create_nan(),s2_2= create_nan();
      double s12=create_nan();
      for (k=0; k<SIZE(use_idx); k++)
      {
        j = MdrV0(use_idx,k) - 1;
        double d1 = mdrV0(x,stride_i * i1 + stride_j * j);
        double d2 = mdrV0(x,stride_i * i2 + stride_j * j);

        if (isnand(d1) || isnand(d2))
          continue;

        if (igninfs)
          if ((create_inf () == d1) || (create_inf () == d2))
            continue;

        if (isnand(s1))
        {
          s1    = d1;
          s1_2  = d1*d1;
        }
        else
        {
          s1    += d1;
          s1_2  += d1*d1;
        }

        if (isnand(s2))
        {
          s2    = d2;
          s2_2  = d2*d2;
        }
        else
        {
          s2    += d2;
          s2_2  += d2*d2;
        }

        if (isnand(s12))
        {
          s12 = d1 * d2;
        }
        else
        {
          s12 += d1 * d2;
        }

        count++;
      }
      s1  /= (double) count;
      s2  /= (double) count;
      s12 /= (double) (count - unbias);

      if (m)
      {
        s1  *= ((double) count / (double)(count - unbias)) * mdrV1(m,MIN(SIZE(m),i2+1));
        s2  *= ((double) count / (double)(count - unbias)) * mdrV1(m,MIN(SIZE(m),i1+1));
        s12 += mdrV1(m,MIN(SIZE(m),i1+1)) * mdrV1(m,MIN(SIZE(m),i2+1));
        Mdr0(w,_i2,_i1) = Mdr0(w,_i1,_i2) = s12 - s1 - s2 ;
      }
      else
      {
        Mdr0(w,i2,i1) = Mdr0(w,i1,i2) = s12
            - (double) (count + unbias) / (double) (count - unbias) * s1 * s2 ;
      }
    }
  }

_exit_mdr_var:

  if (use_r)
    mdr_Destroy(use_r);
  if (ignore_r)
    mdr_Destroy(ignore_r);
  if (use_idx)
    mdr_Destroy(use_idx);
  if (use_factors)
    mdr_Destroy(use_factors);

  return w;
}



//
//
//
MDR * mdr_VectorSet(MDR *t1)
{
  MDR *r=0, *r_idx=0;
  MDR *rval=0;

  int i, n1, cnt;

  if (!EQVECT(t1))
  {
    rval = mdr_Create(0,0);
    return (rval);
  }

  n1 = SIZE(t1);

  if (MD_TYPE_DOUBLE(t1))
  {
    // copy two together:
    r     = mdr_Create(1, n1);
    r_idx = mdr_Create(1, n1);
    for (i=0; i<n1; i++)
    {
      MdrV0(r,i)     = MdrV0(t1,i);
      MdrV0(r_idx,i) = i;
    }

    // sort them
    r_sort ((double *) MDRPTR(r), 0, n1-1, (double *)MDRPTR(r_idx));

    // count how many different ones
    cnt=1;
    for(i=1;i<n1;i++)
    {
      if (isnand(MdrV0(r,i))&&isnand(MdrV0(r,i-1)))
        continue;

      if (MdrV0(r,i)!=MdrV0(r,i-1))
        cnt++;
    }

    // copy different ones to the output array
    rval = mdr_Create(1,cnt);
    cnt=0;
    MdrV0(rval,cnt) = MdrV0(r,0);
    for(i=1;i<n1;i++)
    {
      if (isnand(MdrV0(r,i))&&isnand(MdrV0(r,i-1)))
        continue;

      if (MdrV0(r,i)!=MdrV0(r,i-1))
      {
        cnt++;
        MdrV0(rval,cnt) = MdrV0(r,i);
      }
    }
  }
  else if (MD_TYPE_INT32(t1))
  {
    // copy two together:
    r     = mdi_Create(1, n1);
    r_idx = mdr_Create(1, n1);
    for (i=0; i<n1; i++)
    {
      MdiV0(r,i)     = MdiV0(t1,i);
      MdrV0(r_idx,i) = i;
    }

    // sort them
    i_sort ((int *) MDRPTR(r), 0, n1-1, (double *)MDRPTR(r_idx));

    // count how many different ones
    cnt=1;
    for(i=1;i<n1;i++)
    {
      if (MdiV0(r,i)!=MdiV0(r,i-1))
        cnt++;
    }

    // copy different ones to the output array
    rval = mdi_Create(1,cnt);
    cnt=0;
    MdiV0(rval,cnt) = MdiV0(r,0);
    for(i=1;i<n1;i++)
    {
      if (MdiV0(r,i)!=MdiV0(r,i-1))
      {
        cnt++;
        MdiV0(rval,cnt) = MdiV0(r,i);
      }
    }
  }

  mdr_Destroy(r);
  mdr_Destroy(r_idx);

  return (rval);
}

//
// find union of two discrete vectors
//
MDR * mdr_VectorUnion(MDR *t1, MDR *t2)
{
  MDR *r=0;
  MDR *rval=0;

  int i, n1, n2;

  // special cases
  if (SIZE(t1)==0)
  {
    rval = mdr_VectorSet((MDR *)t2);
    return rval;
  }
  if (SIZE(t2)==0)
  {
    rval = mdr_VectorSet((MDR *)t1);
    return rval;
  }

  if ((!EQVECT(t1)) || (!EQVECT(t2)))
  {
    rval = mdr_Create(0,0);
    return (rval);
  }

  n1 = SIZE(t1);
  n2 = SIZE(t2);

  if (MD_TYPE_INT32(t1) && MD_TYPE_INT32(t2))
  {
    // copy two together:
    r = mdi_Create(1, n1+n2);
    for (i=0; i<n1; i++)
      MdiV0(r,i)     = MdiV0(t1,i);
    for (i=0; i<n2; i++)
      MdiV0(r,i+n1) = MdiV0(t2,i);
  }
  else
  {
    // copy two together:
    r = mdr_Create(1, n1+n2);
    for (i=0; i<n1; i++)
      MdrV0(r,i)     = mdrV0(t1,i);
    for (i=0; i<n2; i++)
      MdrV0(r,i+n1) = mdrV0(t2,i);
  }

  rval = mdr_VectorSet((MDR *)r);
  mdr_Destroy(r);
  return (rval);
}

//
// find intersection of two discrete vectors
//
MDR * mdr_VectorIntersect(MDR *t1, MDR *t2)
{
  MDR *rval=0, *set_c=0;

  if ((!EQVECT(t1)) || (!EQVECT(t2)))
  {
    rval = mdr_Create(0,0);
    return (rval);
  }

  // find set of each vector:
  //  set are ascending and each element appears only once
  MDR * set_a = mdr_VectorSet(t1);
  MDR * set_b = mdr_VectorSet(t2);

  if (MD_TYPE_INT32(t1) && MD_TYPE_INT32(t2))
  {
    set_c = mdi_Create(1,MIN(SIZE(set_a),SIZE(set_b)));
    int ia=0, ib=0, ic=0;
    for (ia=0; ia<SIZE(set_a); ia++)
    {
      while (ib<SIZE(set_b))
      {
        if (MdiV0(set_a,ia) <= MdiV0(set_b,ib))
          break;
        ib++;
      }

      if (ib == SIZE(set_b))
        break;

      if (MdiV0(set_a,ia) == MdiV0(set_b,ib))
      {
        MdiV0(set_c,ic) = MdiV0(set_a,ia);
        ic++;
      }
    }

    if (ic > 0)
    {
      rval = mdi_Create(1,ic);
      for (ia=0; ia<ic; ia++)
        MdiV0(rval,ia) = MdiV0(set_c,ia);
    }
    else
      rval = mdr_Create(0,0);

  }
  else
  {
    set_c = mdr_Create(1,MIN(SIZE(set_a),SIZE(set_b)));
    int ia=0, ib=0, ic=0;
    for (ia=0; ia<SIZE(set_a); ia++)
    {
      while (ib<SIZE(set_b))
      {
        if (MdrV0(set_a,ia) <= MdrV0(set_b,ib))
          break;
        ib++;
      }

      if (ib == SIZE(set_b))
        break;

      if (MdrV0(set_a,ia) == MdrV0(set_b,ib))
      {
        MdrV0(set_c,ic) = MdrV0(set_a,ia);
        ic++;
      }
    }

    if (ic > 0)
    {
      rval = mdr_Create(1,ic);
      for (ia=0; ia<ic; ia++)
        MdrV0(rval,ia) = MdrV0(set_c,ia);
    }
    else
      rval = mdr_Create(0,0);
  }

  // clean up the sets
  mdr_Destroy(set_a);
  mdr_Destroy(set_b);
  mdr_Destroy(set_c);

  return (rval);
}

//
// find complement of two discrete vectors
//
MDR * mdr_VectorComplement(MDR *t1, MDR *t2)
{
  MDR *rval=0;

  // special case: if second vector is [] return the first vector
  if (EQVECT(t1))
  {
    if (SIZE(t2)==0)
    {
      rval = mdr_VectorSet((MDR *)t1);
      return (rval);
    }
  }

  // neither is vector
  if ((!EQVECT(t1)) || (!EQVECT(t2)))
  {
    rval = mdr_Create(0,0);
    return (rval);
  }

  // special case: if first vector is [] return []
  if (SIZE(t1)==0)
  {
    rval = mdr_Create(0,0);
    return (rval);
  }

  // find set of each vector:
  //  set are ascending and each element appears only once
  MDR * set_a = mdr_VectorSet(t1);
  MDR * set_b = mdr_VectorIntersect(t1, t2);
  MDR * set_c = 0;

  if (MD_TYPE_INT32(t1) && MD_TYPE_INT32(t2))
  {
    if (SIZE(set_b) > 0)
    {
      set_c = mdi_Create(1,MAX(SIZE(set_a),SIZE(set_b)));

      int ia, ib=0, ic=0;
      for (ia=0; ia<SIZE(set_a); ia++)
      {
        while (ib<SIZE(set_b))
        {
          if (MdiV0(set_a,ia) <= MdiV0(set_b,ib))
            break;
          ib++;
        }

        if (ib == SIZE(set_b))
        {
          MdiV0(set_c,ic) = MdiV0(set_a,ia);
          ic++;
          continue;
        }

        if (MdiV0(set_a,ia) != MdiV0(set_b,ib))
        {
          MdiV0(set_c,ic) = MdiV0(set_a,ia);
          ic++;
        }
      }

      if (ic > 0)
      {
        rval = mdi_Create(1,ic);
        for (ia=0; ia<ic; ia++)
          MdiV0(rval,ia) = MdiV0(set_c,ia);
      }
      else
        rval = mdi_Create(0,0);
    }
    else
    {
      rval = mdr_Int_BF(set_a);
    }
  }
  else
  {
    if (SIZE(set_b) > 0)
    {
      set_c = mdr_Create(1,MAX(SIZE(set_a),SIZE(set_b)));

      int ia, ib=0, ic=0;
      for (ia=0; ia<SIZE(set_a); ia++)
      {
        while (ib<SIZE(set_b))
        {
          if (mdrV0(set_a,ia) <= mdrV0(set_b,ib))
            break;
          ib++;
        }


        if (ib == SIZE(set_b))
        {
          MdrV0(set_c,ic) = mdrV0(set_a,ia);
          ic++;
          continue;
        }

        if (mdrV0(set_a,ia) != mdrV0(set_b,ib))
        {
          MdrV0(set_c,ic) = mdrV0(set_a,ia);
          ic++;
        }
      }

      if (ic > 0)
      {
        rval = mdr_Create(1,ic);
        for (ia=0; ia<ic; ia++)
          MdrV0(rval,ia) = mdrV0(set_c,ia);
      }
      else
        rval = mdr_Create(0,0);
    }
    else
    {
      rval = mdr_Float_BF(set_a);
    }
  }

  mdr_Destroy(set_a);
  mdr_Destroy(set_b);
  mdr_Destroy(set_c);
  return (rval);
}

//
// sort rows of matrix with respect to column
//
MDR* mdr_SortMatrix_ColIdx(MDR *x, int col_idx)
{
  int nr, nc, i, j, k;
  MDR *r=0, *r_idx=0, *x_sorted=0;

  if (EQNULL(x))
    return mdr_Copy(x);

  nr = MNR(x);
  nc = MNC(x);

  if (MNR(x)==1)
    return mdr_Copy(x);

  if (col_idx<1 || col_idx>nc)
    return mdr_Copy(x);

  // create vector from values in 'col_idx'
  r     = mdr_Create(1, nr);
  r_idx = mdr_Create(1, nr);
  for (i=0; i<nr; i++)
  {
    MdrV0(r,i)     = mdr0(x,i,col_idx-1);
    MdrV0(r_idx,i) = i;
  }

  // sort
  r_sort ((double *) MDRPTR(r), 0, nr-1, (double *)MDRPTR(r_idx));

  //
  if (MD_TYPE_INT32(x))
  {
    x_sorted = mdi_Create(nr,nc);
    for(i=0;i<nr;i++)
    {
      k = (int) MdrV0(r_idx,i);
      for (j=0;j<nc;j++)
        Mdi0(x_sorted,i,j) = Mdi0(x,k,j);
    }
  }
  else if (MD_TYPE_DOUBLE(x))
  {
    x_sorted = mdr_Create(nr,nc);
    for(i=0;i<nr;i++)
    {
      k = (int) MdrV0(r_idx,i);
      for (j=0;j<nc;j++)
        Mdr0(x_sorted,i,j) = Mdr0(x,k,j);
    }
  }
  else
  {
    x_sorted = mdr_Copy(x);
  }

  mdr_Destroy(r);
  mdr_Destroy(r_idx);

  return x_sorted;
}

//
// combine observations 'y' with observation 'x', where first column in each
// matrix contains independent coordinate (time or some such)
//
MDR * mdr_Merge(MDR *y, MDR *x)
{
  MDR *rval;
  int i,j;
  int copy_x=0, copy_y=0;
  MDR *xc=0, *yc=0;

  // special case: if first vector is [] return second vector
  if (SIZE(y)==0 || MNC(y)==1)
  {
    rval = mdr_Float_BF(x);
    return (rval);
  }
  // special case: if second vector is [] return the first vector
  if (SIZE(x)==0 || MNC(x)==1)
  {
    rval = mdr_Float_BF(y);
    return (rval);
  }

  // check that 'x' and 'y' have the first column sorted
  // if not copy them, and sort with respect to first column
  for(i=1; i<MNR(x); i++)
  {
    if (mdr0(x,i-1,0) >= mdr0(x,i,0))
    {
      copy_x=1;
      break;
    }
  }
  for(i=1; i<MNR(y); i++)
  {
    if (mdr0(y,i-1,0) >= mdr0(y,i,0))
    {
      copy_y=1;
      break;
    }
  }

  if (copy_x)
    xc = mdr_SortMatrix_ColIdx(x,1);
  else
    xc = x;

  if (copy_y)
    yc = mdr_SortMatrix_ColIdx(y,1);
  else
    yc = y;

  // find independent column in y
  MDR *ty = mdr_Create(0,0);
  MNR(ty) = MNR(yc);
  MNC(ty) = 1;
  MDPTR(ty) = &Mdr0(yc,0,0);
  // find independent column in x
  MDR *tx = mdr_Create(0,0);
  MNR(tx) = MNR(xc);
  MNC(tx) = 1;
  MDPTR(tx) = &Mdr0(xc,0,0);
  // find their union
  MDR * t_all = mdr_VectorUnion(ty,tx);
  // cleanup
  MDPTR(tx) = 0;
  mdr_Destroy(tx);
  MDPTR(ty) = 0;
  mdr_Destroy(ty);

  // this is the result:
  rval = mdr_Create(SIZE(t_all),MNC(yc) + MNC(xc) - 1);
  mdr_Nan(rval);

  //copy independent
  for (i=0; i<SIZE(t_all); i++)
  {
    Mdr0(rval,i,0) = MdrV0(t_all,i);
  }

  if (SIZE(t_all) == MNR(y))
  {
    //copy y directly into result
    for (i=0; i<MNR(yc); i++)
      for (j=1; j<MNC(yc); j++)
        Mdr0(rval,i,j) = mdr0(yc,i,j);
  }
  else
  {
    // copy rows from y at their updated positions with respect to independent coordinate
    int iy, it=0;
    for (iy=0; iy<MNR(yc); iy++)
    {
      while (     ((Mdr0(yc,iy,0) > MdrV0(t_all,it)) && it<SIZE(t_all))
               || (isnand(MdrV0(t_all,it)) && !isnand(Mdr0(yc,iy,0)))
            )
        it++;
      if (it == SIZE(t_all))
        break;

      if ((Mdr0(yc,iy,0) == MdrV0(t_all,it)) || (isnand(Mdr0(yc,iy,0)) && isnand(MdrV0(t_all,it))))
      {
        for (j=1; j<MNC(y); j++)
        {
          Mdr0(rval,it,j) = mdr0(yc,iy,j);
        }
      }
    }
  }

  if (SIZE(t_all) == MNR(xc))
  {
    //copy x directly into result
    for (i=0; i<MNR(xc); i++)
      for (j=1; j<MNC(xc); j++)
        Mdr0(rval,i,j + MNC(yc) - 1) = mdr0(xc,i,j);
  }
  else
  {
    // copy rows from y at their updated positions with respect to independent coordinate
    int ix=0, it=0;
    for (ix=0; ix<MNR(xc); ix++)
    {
      while (   ((Mdr0(xc,ix,0) > MdrV0(t_all,it)) && it<SIZE(t_all))
             || (isnand(MdrV0(t_all,it)) && !isnand(Mdr0(xc,ix,0)))
            )
        it++;
      if (it == SIZE(t_all))
        break;

      if ((Mdr0(xc,ix,0) == MdrV0(t_all,it)) || (isnand(Mdr0(xc,ix,0)) && isnand(MdrV0(t_all,it))))
      {
        for (j=1; j<MNC(xc); j++)
          Mdr0(rval,it,j+MNC(yc)-1) = mdr0(xc,ix,j);
      }
    }
  }

  mdr_Destroy(t_all);
  if (copy_x)
    mdr_Destroy(xc);
  if (copy_y)
    mdr_Destroy(yc);

  return (rval);
}

//
// compact rows of matrix with respect to column
//  replace multiple rows with the same entry with their mean/std/weight
//
MDR * mdr_Compact(MDR *x, int col_idx, MDR *idx_xw, MDR *idx_xs)
{
  MDR *rval;
  int i,j, dj, is_std=0, cx, cw;
  int copy_x=0;
  MDR *xc=0, *idx_xx=0;

  col_idx--;

  if (idx_xw && idx_xs)
  {
    fprintf(stderr, "compact: only one array of value/deviation or value/weight pairs is allowed!\n");
    rval = mdr_Create(0,0);
    return (rval);
  }
  else if (idx_xs)
  {
    is_std=1;
    idx_xx = idx_xs;
  }
  else
    idx_xx = idx_xw;


  // special case: if first vector is [] return second vector
  if (SIZE(x)==0)
  {
    rval = mdr_Create(0,0);
    return (rval);
  }

  // nothing to do
  if (MNR(x)==1)
    return mdr_Copy(x);

  // has to be double matrix
  if (MD_TYPE_INT32(x))
  {
    fprintf(stderr, "Compact: Argument matrix has to be real!\n");
    return mdr_Create(0,0);
  }

  // check that 'x' is sorted in ascending order: if not
  // create a copy that is
  for(i=1; i<MNR(x); i++)
  {
    if (Mdr0(x,i-1,col_idx) >= Mdr0(x,i,col_idx))
    {
      copy_x=1;
      break;
    }
  }
  if (copy_x)
    xc = mdr_SortMatrix_ColIdx(x,col_idx+1);
  else
    xc = x;

  // find independent column in x
  MDR *tx = mdr_Create(0,0);
  MNR(tx) = MNR(xc);
  MNC(tx) = 1;
  MDPTR(tx) = &Mdr0(xc,0,col_idx);

  // find its set
  MDR * t_set = mdr_VectorSet(tx);
  // cleanup
  MDPTR(tx) = 0;
  mdr_Destroy(tx);

  // this is the result:
  rval = mdr_Create(SIZE(t_set),MNC(xc));
  mdr_Nan(rval);

  //copy independent column
  for (i=0; i<SIZE(t_set); i++)
    Mdr0(rval,i,col_idx) = MdrV0(t_set,i);

  if (SIZE(t_set) == MNR(xc))
  {
    //same sizes: nothing to compact, thus copy x directly into result
    for (j=0; j<MNC(xc); j++)
    {
      if (j == col_idx)
        continue;
      for (i=0; i<MNR(xc); i++)
        Mdr0(rval,i,j) = Mdr0(xc,i,j);
    }
  }
  else
  {
    // copy rows from x at their updated positions with respect to independent coordinate
    int ix1=0, ix2, it=0, iy;
    while (ix1<MNR(xc))
    {
      // initialize first entry
      for (j=0; j<MNC(xc); j++)
      {
        if (j==col_idx)
          continue;
        Mdr0(rval,it,j) = Mdr0(xc,ix1,j);
      }

      // nothing to do: get out!
      if (ix1==MNR(xc)-1)
        break;

      // are there more entries with same x[;col_idx]?
      ix2=ix1+1;
      if (isnand(Mdr0(xc,ix1,col_idx)) && ix2<MNR(xc))
      {
        while (isnand(Mdr0(xc,ix2,col_idx)) && ix2<MNR(xc))
          ix2++;
      }
      else
      {
        while ((Mdr0(xc,ix1,col_idx) == Mdr0(xc,ix2,col_idx)) && ix2<MNR(xc))
          ix2++;
      }

      // process rows which have the same col_idx: ix1..(ix2)
      dj = ix2 - ix1;
      if (dj>1)
      {
        if (idx_xx)
        {
          for (i=0; i<MNR(idx_xx); i++)
          {
            cx = (int) mdr0(idx_xx, i, 0)-1;
            cw = (int) mdr0(idx_xx, i, 1)-1;
            rlab_meanstat_from_valstat_vectors(dj, &Mdr0(xc,ix1,cx),  &Mdr0(xc,ix1,cw),
                NULL, NULL, &Mdr0(rval,it,cx), &Mdr0(rval,it,cw), is_std);
          }
        }
        else
        {
          for (iy=0; iy<MNC(xc); iy++)
          {
            if (iy==col_idx)
              continue;
            Mdr0(rval,it,iy) = rlab_mean_vector(dj, &Mdr0(xc,ix1,iy), NULL, NULL, 1);
          }
        }
      }

      ix1 = ix2;
      it++;
    } // while(ix1<MNR(xc))
  }

  mdr_Destroy(t_set);
  if (copy_x)
    mdr_Destroy(xc);

  return (rval);
}



MDR *
mdr_Any (MDR * m)
{
  int i, j;
  MDR *new;

  if(MD_TYPE_INT32(m))
  {
    if (EQVECT(m))		/* Vector operation */
    {
      new = mdi_Create (1, 1);
      MdiV0 (new,0) = 0;
      for (i=0; i<SIZE(m); i++)
      {
        if (MdiV0 (m,i) != 0)
        {
          MdiV0 (new,0) = 1;
          break;
        }
      }
    }
    else
      /* Matrix operation */
    {
      new = mdi_Create (1, MNC (m));
      for (i=0; i<MNC(m); i++)
      {
        MdiV0 (new,i) = 0;
        for (j=0; j<MNR(m); j++)
        {
          if (Mdi0 (m, j, i) != 0)
          {
            MdiV0 (new, i) = 1;
            break;
          }
        }
      }
    }
  }
  else
  {
    if (MNR (m) == 1)		/* Vector operation */
    {
      new = mdr_Create (1, 1);
      Mdr1 (new, 1, 1) = 0.0;
      for (i = 1; i <= MNC (m); i++)
      {
        if (MdrV1 (m, i) != 0.0)
        {
          MdrV1 (new, 1) = 1.0;
          break;
        }
      }
    }
    else
    /* Matrix operation */
    {
      new = mdr_Create (1, MNC (m));
      for (i = 1; i <= MNC (m); i++)
      {
        Mdr1 (new, 1, i) = 0.0;
        for (j = 1; j <= MNR (m); j++)
        {
          if (Mdr1 (m, j, i) != 0.0)
          {
            Mdr1 (new, 1, i) = 1.0;
            break;
          }
        }
      }
    }
  }

  return (new);
}

MDR *
mdr_All (MDR * m)
{
  int i, j;
  MDR *new;

  if (EQVECT(m))    /* Vector operation */
  {
    if(m->type == RLAB_TYPE_INT32)
    {
      new = mdi_Create (1, 1);
      Mdi1 (new, 1, 1) = 1;
      for (i = 1; i <= SIZE(m); i++)
      {
        if (MdiV1(m,i) == 0)
        {
          Mdi1 (new, 1, 1) = 0;
          break;
        }
      }
    }
    else
    {
      new = mdr_Create (1, 1);
      Mdr1 (new, 1, 1) = 1;
      for (i = 1; i <= SIZE(m); i++)
      {
        if (MdrV1(m,i) == 0)
        {
          Mdr1 (new, 1, 1) = 0;
          break;
        }
      }
    }
    return (new);
  }

  if(m->type == RLAB_TYPE_INT32)
  {
    new = mdi_Create (1, MNC (m));
    for (i=1; i<=MNC (m); i++)
    {
      MdiV1 (new, i) = 1;
      for (j=1; j<=MNR(m); j++)
      {
        if (Mdi1 (m, j, i) == 0)
        {
          MdiV1 (new, i) = 0;
          break;
        }
      }
    }
  }
  else
  {
    new = mdr_Create (1, MNC (m));
    for (i=1; i<=MNC (m); i++)
    {
      MdrV1 (new, i) = 1.0;
      for (j=1; j<=MNR(m); j++)
      {
        if (Mdr1 (m, j, i) == 0.0)
        {
          MdrV1 (new, i) = 0.0;
          break;
        }
      }
    }
  }

  return (new);
}

/* **************************************************************
 * Get the size of a a matrix for size()...
 * ************************************************************** */

MDR *
mdr_Size_BF (MDR * m)
{
  int nr = MNR (m);
  int nc = MNC (m);
  //if (nr == 1 || nc == 1)
  //{ MDR *size = mdr_CreateScalar ((double) nr*nc);}
  //else
  //{
    MDR *size = mdr_Create (1, 2);
    Mdr0 (size, 0, 0) = (double) nr;
    Mdr0 (size, 0, 1) = (double) nc;
  //}
  return (size);
}

MDR *
mdr_Length_BF (MDR * m)
{
  MDR *size = mdr_Create (1, 1);
  MdrV0 (size, 0) = MAX(MNR (m), MNC (m));
  return (size);
}

/* **************************************************************
 * Return the type of a matrix...
 * ************************************************************** */

MDS *
mdr_Type_BF (MDR * m)
{
  MDS *type;
  if(m->type == RLAB_TYPE_INT32)
    type = mds_CreateScalar ( "int" );
  else
    type = mds_CreateScalar ( "real" );
  return (type);
}

/* **************************************************************
 * Mod function: matrix optimized
 * ************************************************************** */

MDR * mdr_Mod_BF (MDR * m1, MDR * m2)
{
  int i, j, nr, nc;
  MDR *m=0;

  nr = MAX( MNR ( m1 ), MNR ( m2 ) );
  nc = MAX( MNC ( m1 ), MNC ( m2 ) );

  if(MD_TYPE_INT32(m1) && MD_TYPE_INT32(m2))
  {
    m = mdi_Create (nr, nc);
    for (i = 0; i < nr; i++)
    {
      for (j = 0; j < nc; j++)
      {
        if (mdi0_safe(m2, i, j))
        {
          Mdi0 (m, i, j) = mdi0_safe (m1, i, j) % mdi0_safe(m2, i, j);
          if (Mdi0 (m, i, j) < 0)
            Mdi0 (m, i, j) += mdi0_safe(m2, i, j);
        }
        else
          Mdi0 (m, i, j) = mdi0_safe (m1, i, j);
      }
    }
  }
  else
  {
    m = mdr_Create (nr, nc);
    mdr_Zero(m);
    for (i = 0; i < nr; i++)
    {
      for (j = 0; j < nc; j++)
      {
        if (mdr0_safe(m2,i,j) == 0)
          Mdr0 (m, i, j) = mdr0_safe (m1,i,j);
        else if (mdr0_safe(m2,i,j) == 1)
          Mdr0 (m, i, j) = 0.0;
        else
        {
          Mdr0 (m, i, j) = fmod( mdr0_safe (m1,i,j), mdr0_safe(m2,i,j) );
          if (Mdr0 (m, i, j) < 0)
          {
            if (mdr0_safe(m2,i,j) > 0)
              Mdr0 (m, i, j) += mdr0_safe(m2,i,j);
            else
              Mdr0 (m, i, j) -= mdr0_safe(m2,i,j);
          }
        }
      }
    }
  }
  return (m);
}

/* **************************************************************
 * Copy a matrix and change its type to real.
 * ************************************************************** */

MDR * mdr_Float_BF (MDR * m)
{
  int k;
  MDR *new = 0;

  if(m->type == RLAB_TYPE_INT32)
  {
    new = mdr_Create (MNR (m), MNC (m));
    for(k=0; k<=SIZE(m); k++)
      MdrV0(new,k) = MdiV0(m,k);
  }
  else
  {
    new = mdr_Copy (m);
  }
  return (new);
}


/* **************************************************************
 * Int function:
 *  returns matrix dense integer
 * ************************************************************** */

MDR *
mdr_Int_BF (MDR * m)
{
  int k, isize;
  MDR *im;

  //im = mdr_Create (MNR (m), MNC (m));
  if(m->type == RLAB_TYPE_INT32)
  {
    im = mdr_Copy (m);
  }
  else
  {
    im = mdi_Create (MNR (m), MNC (m));
    isize = MNR (m) * MNC (m);
    for (k = 0; k < isize; k++)
    {
      MdiV0 (im, k) = (int) MdrV0 (m, k);
    }
  }
  return (im);
}

/* **************************************************************
 * Ceil/Floor function...
 * ************************************************************** */

MDR * mdr_CeilFlooRound_BF (int is_ceil, MDR * m, MDR *b, MDR *o, MDR *c)
{
  int i, j;
  MDR *bin=0, *offs=0;
  MDR *cm=0;
  double (*ceilflooround) ();

  if (!md_Isvalid(m))
    return mdr_Create(0,0);

  if (md_Isvalid(b))
    bin = b;
  else
    bin = mdr_CreateScalar(1.0);

  if (md_Isvalid(o))
    offs = o;
  else
    offs = mdr_CreateScalar(0.0);


  int r1 = MNR(m);
  int c1 = MNC(m);
  int r2 = MNR(bin);
  int c2 = MNC(bin);
  int r3 = MNR(offs);
  int c3 = MNC(offs);
  int rr = MAX(MAX(r1,r2),r3);
  int rc = MAX(MAX(c1,c2),c3);

  cm = mdr_Create (rr, rc);

  switch(is_ceil)
  {
    case -1:
      ceilflooround = round;
      break;

    case 0:
      ceilflooround = floor;
      break;

    default:
      ceilflooround = ceil;
      break;
  }

  for (i=0; i<rr; i++)
  {
    for (j=0; j<rc; j++)
    {
      double x1 = mdr0(m,   MIN(i,r1-1),MIN(j,c1-1));
      double o1 = Mdr0(offs,MIN(i,r3-1),MIN(j,c3-1));
      double b1 = Mdr0(bin, MIN(i,r2-1),MIN(j,c2-1));
      if (b1 == 0)
      {
        // bin width is zero: do nothing
        Mdr0 (cm,i,j) = x1;
        continue;
      }
      double d1 = ceilflooround((x1-o1)/b1);
      Mdr0 (cm,i,j) = o1 + b1 * d1;
      if (c)
      {
        if (Mdr0 (cm,i,j) < mdrV0(c,0))
        {
          Mdr0 (cm,i,j) = mdrV0(c,0);
          continue;
        }
        if (Mdr0 (cm,i,j) > mdrV0(c,1))
        {
          Mdr0 (cm,i,j) = mdrV0(c,1);
        }
      }
    }
  }

  if (!b)
    mdr_Destroy(bin);
  if (!o)
    mdr_Destroy(offs);

  return (cm);
}

MDR *
mdr_LogCeilFlooRound_BF (int is_ceil, MDR * m, double base)
{
  int i, j;
  MDR *cm=0;
  double (*ceilflooround) ();

  if (!md_Isvalid(m) || base<=0)
    return mdr_Create(0,0);

  double log_base = log(base);

  int r1 = MNR(m);
  int c1 = MNC(m);
  cm = mdr_Create (r1, c1);

  switch(is_ceil)
  {
    case -1:
      ceilflooround = round;
      break;

    case 0:
      ceilflooround = floor;
      break;

    default:
      ceilflooround = ceil;
      break;
  }

  for (i=0; i<r1; i++)
  {
    for (j=0; j<c1; j++)
    {
      double x1 = mdr0(m,i,j);
      if (x1>0)
      {
        double d1 = ceilflooround(log(x1)/log_base);
        Mdr0 (cm,i,j) = pow(base, d1);
      }
      else
        Mdr0 (cm,i,j) = x1;
    }
  }

  return (cm);
}

#include <gsl/gsl_histogram.h>
#include <gsl/gsl_errno.h>
struct _rlab_gsl_histrogram
{
  size_t n;
  double *range;
  double *bin;
};


MDR *
mdr_MeshCeilFlooRound_BF (int is_ceil, MDR * m, MDR *mesh)
{
  int i, j, nm;
  MDR *cm=0;
  struct _rlab_gsl_histrogram dummy_hist;
  gsl_histogram *h;

  if (!md_Isvalid(m) || !md_Isvalid(mesh))
    return mdr_Create(0,0);

  if (!mdr_vector_issorted(mesh))
    return mdr_Create(0,0);

  nm = SIZE(mesh);

  int r1 = MNR(m);
  int c1 = MNC(m);
  cm = mdr_Create (r1, c1);

  // create dummy histogram for GSL library functions
  h = (gsl_histogram *) &dummy_hist;

  dummy_hist.n = nm - 1;
  dummy_hist.range = MDRPTR(mesh);

  for (i=0; i<r1; i++)
  {
    for (j=0; j<c1; j++)
    {
      double x1 = mdr0(m,i,j);
      size_t loc;
      if(gsl_histogram_find(h,x1,&loc) == GSL_SUCCESS)
      {
        if (x1 == MdrV0(mesh,loc))
        {
          Mdr0 (cm,i,j) = MdrV0(mesh,loc);
          continue;
        }
        if (x1 == MdrV0(mesh,loc+1))
        {
          Mdr0 (cm,i,j) = MdrV0(mesh,loc+1);
          continue;
        }

        if (is_ceil==1 && loc < nm - 1)
          Mdr0 (cm,i,j) = MdrV0(mesh,loc+1);
        else if (is_ceil==0)
          Mdr0 (cm,i,j) = MdrV0(mesh,loc);
      }
      else
        Mdr0 (cm,i,j) = create_nan();
    }
  }

  return (cm);
}

/* **************************************************************
 * Floor function...
 * ************************************************************** */

MDR * mdr_Floor_BF (MDR * m, MDR *b, MDR *o)
{
  int i, j;
  MDR *bin=0, *offs=0;
  MDR *cm=0;

  if (!m)
    return mdr_Create(0,0);

  if (b)
    bin = b;
  else
    bin = mdr_CreateScalar(1.0);

  if (o)
    offs = o;
  else
    offs = mdr_CreateScalar(0.0);


  int r1 = MNR(m);
  int c1 = MNC(m);
  int r2 = MNR(bin);
  int c2 = MNC(bin);
  int r3 = MNR(offs);
  int c3 = MNC(offs);
  int rr = MAX(MAX(r1,r2),r3);
  int rc = MAX(MAX(c1,c2),c3);

  cm = mdr_Create (rr, rc);

  for (i=0; i<MNR(cm); i++)
  {
    for (j=0; j<MNC(cm); j++)
    {
      double x1 = mdr0(m,MIN(i,r1-1),MIN(j,c1-1));
      double o1 = Mdr0(offs,MIN(i,r3-1),MIN(j,c3-1));
      double b1 = Mdr0(bin, MIN(i,r2-1),MIN(j,c2-1));
      double d1 = (x1-o1)/b1;
      Mdr0 (cm,i,j) = o1 + b1 * (int)(d1);
    }
  }

  if (!b)
    mdr_Destroy(bin);
  if (!o)
    mdr_Destroy(offs);

  return (cm);
}
// (MDR * m)
// {
//   int i, size;
//   MDR *cm;
//
//   if(m->type == RLAB_TYPE_INT32)
//     cm = mdr_Float_BF (m);
//   else
//   {
//     cm = mdr_Create (MNR (m), MNC (m));
//     size = MNR (m) * MNC (m);
//     for (i = 0; i < size; i++)
//     {
//       MdrV0 (cm, i) = errcheck (floor (MdrV0 (m, i)), "floor");
//     }
//   }
//   return (cm);
// }

/* **************************************************************
 * Round function...
 * ************************************************************** */

MDR *
mdr_Round_BF (MDR * m)
{
  int i, size;
  MDR *cm;

  if(m->type == RLAB_TYPE_INT32)
    cm = mdr_Float_BF (m);
  else
  {
    cm = mdr_Create (MNR (m), MNC (m));
    size = MNR (m) * MNC (m);
    for (i = 0; i < size; i++)
    {
      MdrV0 (cm, i) = errcheck (rint (MdrV0 (m, i)), "round");
    }
  }
  return (cm);
}

/* **************************************************************
 * Abs function...
 * ************************************************************** */

MDR *
mdr_Abs (MDR * m)
{
  int i, isize;
  MDR *new;

  isize = SIZE(m);

  new = mdr_Copy (m);
  if(m->type == RLAB_TYPE_INT32)
  {
    for (i = 0; i < isize; i++)
    {
      if ( MdiV0 (m, i) >= 0 ) continue;
      MdiV0 (new, i) = -MdiV0 (m, i);
    }
  }
  else
  {
    for (i = 0; i < isize; i++)
    {
      if ( MdrV0 (m, i) >= 0 ) continue;
      MdrV0 (new, i) = -MdrV0 (m, i);
    }
  }
  return (new);
}

//
//
// {min,max}{val,idx}
//
//
/* **************************************************************
 * what==2: Max function
 * ************************************************************** */
/* **************************************************************
 * what==1: Max-Index function...
 * ************************************************************** */
/* **************************************************************
 * what==-2: Min function...
 * ************************************************************** */
/* **************************************************************
 * what==-1: Min-Index function...
 * ************************************************************** */

MDR * mdr_MinMax1_ValIdx (MDR * m, int what)
{
  int i, j, maxi=1, mini=1;
  double maxr=0, minr=0;
  MDR *mmax=0;

  if ( EQNULL(m) )
    return mdr_Create (0, 0);

  if ( EQVECT(m) )
  {
    minr = maxr = mdrV0 (m, 0);
    mini = maxi = 1;
    for (i=1; i<SIZE(m); i++)
    {
      if (isnand (mdrV0 (m, i)))
        continue;

      if (isnand (maxr))
      {
        maxr = mdrV0 (m, i);
        maxi = i+1;
      }

      if (isnand (minr))
      {
        minr = mdrV0 (m, i);
        mini = i+1;
      }

      // max val/index
      if (mdrV0 (m, i) > maxr)
      {
        maxr = mdrV0 (m, i);
        maxi = i+1;
        continue;
      }

      // min val/index
      if (mdrV0 (m, i) < minr)
      {
        minr = mdrV0 (m, i);
        mini = i+1;
        continue;
      }
    }

    switch (what)
    {
      case 2:
        mmax = mdr_CreateScalar(maxr);
        break;

      case 1:
        mmax = mdr_CreateScalar(maxi);
        break;

      case -1:
        mmax = mdr_CreateScalar(mini);
        break;

      case -2:
        mmax = mdr_CreateScalar(minr);
        break;

      default:
        break;
    }
  }
  else
  {
    mmax = mdr_Create (1, MNC (m));
    for (i = 0; i < MNC (m); i++)
    {
      minr = maxr = mdr0 (m, 0, i);
      mini = maxi = 1;
      for (j = 1; j < MNR (m); j++)
      {
        if (isnand (mdr0 (m, j, i)))
          continue;

        if (isnand (maxr))
        {
          maxr = mdr0 (m, j, i);
          maxi = j+1;
        }

        if (isnand (minr))
        {
          minr = mdr0 (m, j, i);
          mini = j+1;
        }

        if (mdr0 (m, j, i) > maxr)
        {
          maxr = mdr0 (m, j, i);
          maxi = j+1;
        }

        if (mdr0 (m, j, i) < minr)
        {
          minr = mdr0 (m, j, i);
          mini = j+1;
        }
      }
      switch (what)
      {
        case 2:
          MdrV0 (mmax, i) = maxr;
          break;

        case 1:
          MdrV0 (mmax, i) = maxi;
          break;

        case -1:
          MdrV0 (mmax, i) = mini;
          break;

        case -2:
          MdrV0 (mmax, i) = minr;
          break;

        default:
          break;
      }
    }
  }

  return (mmax);
}



/* **************************************************************
 * Max() with 2 arguments...
 * ************************************************************** */
MDR *
mdr_Max2 (MDR * m1, MDR * m2)
{
  int i1, i2, j1, j2, i, j, ir, jr;
  MDR *r = 0;

  /* Check sizes */
  if (EQNULL(m1))
    return mdr_Copy(m2);

  if (EQNULL(m2))
    return mdr_Copy(m1);

  i1 = MNR(m1);
  j1 = MNC(m1);
  i2 = MNR(m2);
  j2 = MNC(m2);
  ir = MAX(i1,i2);
  jr = MAX(j1,j2);

  if( MD_TYPE_INT32(m1) && MD_TYPE_INT32(m2))
  {
    r = mdi_Create(ir,jr);
    for (i=0; i<ir; i++)
    {
      for (j=0; j<jr; j++)
      {
        Mdi0(r, MIN(i,ir-1), MIN(j,jr-1)) =
            MAX(Mdi0(m1, MIN(i,i1-1), MIN(j,j1-1)), Mdi0(m2, MIN(i,i2-1), MIN(j,j2-1)));
      }
    }
  }
  else
  {
    r = mdr_Create(ir,jr);
    for (i=0; i<ir; i++)
    {
      for (j=0; j<jr; j++)
      {
        if (isnand(mdr0(m1, MIN(i,i1-1), MIN(j,j1-1))))
          Mdr0(r, i, j) = mdr0(m1, MIN(i,i1-1), MIN(j,j1-1));
        else if (isnand(mdr0(m2, MIN(i,i2-1), MIN(j,j2-1))))
          Mdr0(r, i, j) = mdr0(m2, MIN(i,i2-1), MIN(j,j2-1));
        else
          Mdr0(r, i, j) =
              MAX(mdr0(m1, MIN(i,i1-1), MIN(j,j1-1)), mdr0(m2, MIN(i,i2-1), MIN(j,j2-1)));
      }
    }
  }

  return (r);
}


/* **************************************************************
 * Min-Index function in 2-D
 * ************************************************************** */
/*
    This function searches the distance matrix to find the pair with the shortest
    distance between them. The indices of the pair are returned in ip and jp; the
    distance itself is returned by the function.
*/
void
find_idxs_min_max(MDR * distmatrix, int* ip, int* jp, double * distance, int *dir)
{
  int i, j, nr=MNR(distmatrix), nc=MNC(distmatrix);
  if (!(nr*nc))
  {
    *distance = create_nan ();
    *ip = -1;
    *jp = -1;
    return;
  }

  *ip = 0;
  *jp = 0;

  *distance = mdr0(distmatrix,*ip,*jp);

  for (i=0; i<nr; i++)
  {
    for (j=0; j<nc; j++)
    {
      // looking for MIN
      if (*dir==-1)
        if (mdr0(distmatrix,i,j)<*distance)
        {
          *distance = mdr0(distmatrix,i,j);
          *ip = i;
          *jp = j;
        }
      // looking for MAX
      if (*dir==1)
        if (mdr0(distmatrix,i,j)>*distance)
        {
          *distance = mdr0(distmatrix,i,j);
          *ip = i;
          *jp = j;
        }
    }
  }
  return;
}


/* **************************************************************
 * Min() with 2 arguments...
 * ************************************************************** */

MDR *
mdr_Min2 (MDR * m1, MDR * m2)
{
  int i1, i2, j1, j2, i, j, ir, jr;
  MDR *r = 0;

  /* Check sizes */
  if (EQNULL(m1))
    return mdr_Copy(m2);

  if (EQNULL(m2))
    return mdr_Copy(m1);

  i1 = MNR(m1);
  j1 = MNC(m1);
  i2 = MNR(m2);
  j2 = MNC(m2);
  ir = MAX(i1,i2);
  jr = MAX(j1,j2);

  if( MD_TYPE_INT32(m1) && MD_TYPE_INT32(m2))
  {
    r = mdi_Create(ir,jr);
    for (i=0; i<ir; i++)
    {
      for (j=0; j<jr; j++)
      {
        Mdi0(r, i, j) =
            MIN(Mdi0(m1, MIN(i,i1-1), MIN(j,j1-1)), Mdi0(m2, MIN(i,i2-1), MIN(j,j2-1)));
      }
    }
  }
  else
  {
    r = mdr_Create(ir,jr);
    for (i=0; i<ir; i++)
    {
      for (j=0; j<jr; j++)
      {
        if (isnand(mdr0(m1, MIN(i,i1-1), MIN(j,j1-1))))
          Mdr0(r, i, j) = mdr0(m1, MIN(i,i1-1), MIN(j,j1-1));
        else if (isnand(mdr0(m2, MIN(i,i2-1), MIN(j,j2-1))))
          Mdr0(r, i, j) = mdr0(m2, MIN(i,i2-1), MIN(j,j2-1));
        else
          Mdr0(r, i, j) =
              MIN(mdr0(m1, MIN(i,i1-1), MIN(j,j1-1)), mdr0(m2, MIN(i,i2-1), MIN(j,j2-1)));
      }
    }
  }

  return (r);
}

//***************************************************************
//* Compute the "sum" of a matrix.
//***************************************************************

MDR *
mdr_Sum_BF (MDR * m, void * h)
{
  int i, j;
  MDR *rm=0, *cond=0;
  int lcr = 0, lcc = 0;

  if (h)
    cond = (MDR *) h;

  if ( EQNULL(m) )
  {
    rm = mdr_CreateScalar (0);
    return (rm);
  }

  if ( EQVECT(m) )
  {
    if (cond)
      lcc = SIZE(cond) - 1;

    if(MD_TYPE_INT32(m))
    {
      // single row or single column matrix: add them all up
      int d = 0;
      for (i=0; i<SIZE(m); i++)
      {
        if (cond)
          if (!mdrV0(cond,MIN(i,lcc)))
            continue;

        d = d + MdiV0 (m, i);
      }

      rm = mdi_CreateScalar (d);
    }
    else
    {
      // single row or single column matrix: add them all up
      double d = 0.0;
      for (i=0; i<SIZE(m); i++)
      {
        if (isnand(MdrV0 (m, i)))
          continue;

        if (cond)
          if (!mdrV0(cond,MIN(i,lcc)))
            continue;

        d = d + MdrV0 (m, i);
      }

      rm = mdr_CreateScalar (d);
    }
  }
  else
  {
    if (cond)
    {
      lcr = MNR(cond) - 1;
      lcc = MNC(cond) - 1;
    }

    // a matrix: add entries column-wise
    if(MD_TYPE_INT32(m))
    {
      rm = mdi_Create (1, MNC(m));
      for (i=0; i<MNC(m); i++)
      {
        MdiV0 (rm, i) = 0;
        for (j=0; j<MNR(m); j++)
        {
          if (cond)
            if (!mdr0(cond,MIN(j,lcr),MIN(i,lcc)))
              continue;

          MdiV0 (rm, i) = MdiV0 (rm, i) + Mdi0 (m, j, i);
        }
      }
    }
    else
    {
      rm = mdr_Create (1, MNC(m));
      for (i = 0; i < m->ncol; i++)
      {
        MdrV0 (rm, i) = 0;
        for (j = 0; j < m->nrow; j++)
        {
          if (isnand(Mdr0 (m, j, i)))
            continue;

          if (cond)
            if (!mdr0(cond,MIN(j,lcr),MIN(i,lcc)))
              continue;

          MdrV0 (rm, i) = MdrV0 (rm, i) + Mdr0 (m, j, i);
        }
      }
    }
  }

  return (rm);
}

MDR *
mdr_CumSum_BF (MDR * m)
{
  int i, j;
  MDR *new;

  if (MNR (m) == 1 || MNC (m) == 1)
  {
    if( MD_TYPE_INT32(m) )
    {
      int N = MNR (m) * MNC (m);
      new = mdi_Create (MNR (m), MNC (m));
      MdiV0 (new, 0) = MdiV0 (m, 0);
      for (i = 1; i < N; i++)
      {
        MdiV0 (new, i) = MdiV0 (new, i - 1) + MdiV0 (m, i);
      }
    }
    else
    {
      int N = MNR (m) * MNC (m);
      new = mdr_Create (MNR (m), MNC (m));

      // suming non-numbers is zero!
      MdrV0 (new, 0) = 0;
      if (!isnand(MdrV0 (m, 0)))
        MdrV0 (new, 0) = MdrV0 (m, 0);

      for (i=1; i<N; i++)
      {
        if (isnand(MdrV0 (m, i)))
          MdrV0 (new, i) = MdrV0 (new, i - 1);
        else
          MdrV0 (new, i) = MdrV0 (new, i - 1) + MdrV0 (m, i);
      }
    }
  }
  else
  {
    if(m->type == RLAB_TYPE_INT32)
    {
      new = mdi_Create (MNR (m), MNC (m));
      /* Initialize first row. */
      for (i = 0; i < MNC (m); i++)
      {
        Mdi0 (new, 0, i) = Mdi0 (m, 0, i);
      }
      /* Now compute running sum. */
      for (i = 0; i < MNC (m); i++)
      {
        for (j = 1; j < MNR (m); j++)
        {
          Mdi0 (new, j, i) = Mdi0 (new, j - 1, i) + Mdi0 (m, j, i);
        }
      }
    }
    else
    {
      new = mdr_Create (MNR (m), MNC (m));
      /* Initialize first row. */
      for (i = 0; i < MNC (m); i++)
      {
        Mdr0 (new, 0, i) = 0;
        if (!isnand(Mdr0 (m, 0, i)))
          Mdr0 (new, 0, i) = Mdr0 (m, 0, i);
      }
      /* Now compute running sum. */
      for (i = 0; i < MNC (m); i++)
      {
        for (j = 1; j < MNR (m); j++)
        {
          if (isnand(Mdr0 (m, j, i)))
            Mdr0 (new, j, i) = Mdr0 (new, j - 1, i);
          else
            Mdr0 (new, j, i) = Mdr0 (new, j - 1, i) + Mdr0 (m, j, i);
        }
      }
    }
  }
  return (new);
}

/* **************************************************************
 * Compute the "product" of a matrix.
 * ************************************************************** */

MDR *
mdr_Prod_BF (MDR * m)
{
  int i, j;
  MDR *rm;

  if (MNR (m) * MNC (m) == 1)
  {
    rm = mdr_CreateScalar (0.0);
  }
  else if (MNR (m) == 1)
  {
    if(m->type == RLAB_TYPE_INT32)
    {
      int d = MdiV0 (m, 0);
      for (i = 1; i < MNC (m); i++)
      {
        d = d * MdiV0 (m, i);
      }
      rm = mdi_CreateScalar (d);
    }
    else
    {
      double d = MdrV0 (m, 0);
      for (i = 1; i < MNC (m); i++)
      {
        if (!isnand(MdrV0 (m, i)))
          d = d * MdrV0 (m, i);
      }
      rm = mdr_CreateScalar (d);
    }
  }
  else
  {
    if(m->type == RLAB_TYPE_INT32)
    {
      rm = mdi_Create (1, MNC (m));
      for (i = 0; i < MNC (m); i++)
      {
        MdiV0 (rm, i) = Mdi0 (m, 0, i);
        for (j = 1; j < MNR (m); j++)
        {
          MdiV0 (rm, i) = MdiV0 (rm, i) * Mdi0 (m, j, i);
        }
      }
    }
    else
    {
      rm = mdr_Create (1, MNC (m));
      for (i = 0; i < MNC (m); i++)
      {
        MdrV0 (rm, i) = Mdr0 (m, 0, i);
        for (j = 1; j < MNR (m); j++)
        {
          if (!isnand(Mdr0 (m, j, i)))
            MdrV0 (rm, i) = MdrV0 (rm, i) * Mdr0 (m, j, i);
        }
      }
    }
  }
  return (rm);
}

MDR *
mdr_CumProd_BF (MDR * m)
{
  int i, j;
  MDR *new;

  if (MNR (m) == 1 || MNC (m) == 1)
  {
    if(m->type == RLAB_TYPE_INT32)
    {
      int N = MNR (m) * MNC (m);
      new = mdi_Create (MNR (m), MNC (m));
      MdiV0 (new, 0) = MdiV0 (m, 0);
      for (i = 1; i < N; i++)
      {
        MdiV0 (new, i) = MdiV0 (new, i - 1) * MdiV0 (m, i);
      }
    }
    else
    {
      int N = MNR (m) * MNC (m);
      new = mdr_Create (MNR (m), MNC (m));
      MdrV0 (new, 0) = MdrV0 (m, 0);
      for (i = 1; i < N; i++)
      {
        MdrV0 (new, i) = MdrV0 (new, i - 1) * MdrV0 (m, i);
      }
    }
  }
  else
  {
    if(m->type == RLAB_TYPE_INT32)
    {
      new = mdi_Create (MNR (m), MNC (m));
      /* Initialize first row. */
      for (i = 0; i < MNC (m); i++)
      {
        Mdi0 (new, 0, i) = Mdi0 (m, 0, i);
      }
      /* Now compute running product. */
      for (i = 0; i < MNC (m); i++)
      {
        for (j = 1; j < MNR (m); j++)
        {
          Mdi0 (new, j, i) = Mdi0 (new, j - 1, i) * Mdi0 (m, j, i);
        }
      }
    }
    else
    {
      new = mdr_Create (MNR (m), MNC (m));
      /* Initialize first row. */
      for (i = 0; i < MNC (m); i++)
      {
        Mdr0 (new, 0, i) = Mdr0 (m, 0, i);
      }
      /* Now compute running product. */
      for (i = 0; i < MNC (m); i++)
      {
        for (j = 1; j < MNR (m); j++)
        {
          Mdr0 (new, j, i) = Mdr0 (new, j - 1, i) * Mdr0 (m, j, i);
        }
      }
    }
  }
  return (new);
}
/* **************************************************************
 * Test a row or column of a matrix for ordering: exclude NaNs
 * ************************************************************** */
int mdr_IsSorted (MDR * m, int row_dominant, int ifact)
{
  int i, ndata, ndim;

  if(SIZE(m)<1)
    return 0;

  if(SIZE(m)==1)
    return 1;

  if (row_dominant)
  {
    ndata = MNC (m);
    ndim  = MNR (m);
    if (ifact<1 || ifact>ndim)
      return 0;
    for (i=1; i<ndata; i++)
    {
      if (mdr0(m,ifact-1,i-1)>mdr0(m,ifact-1,i))
        return 0;
    }
    return 1;
  }
  else
  {
    ndata = MNR (m);
    ndim  = MNC (m);
    if (ifact<1 || ifact>ndim)
      return 0;
    for (i=1; i<ndata; i++)
    {
      if (mdr0(m,i-1,ifact-1)>mdr0(m,i,ifact-1))
        return 0;
    }
  }
  return 1;
}

/* **************************************************************
 * Test a matrix for symmetry.
 * ************************************************************** */

int
    mdr_IsSymmetric (MDR * m)
{
  int i, j, nr, nc;

  if (MNR (m) != MNC (m))
    return (0);

  nr = MNC (m);
  nc = MNC (m);

  if(m->type == RLAB_TYPE_INT32)
  {
    for (j = 0; j < nc; j++)
    {
      for (i = j + 1; i < nr; i++)
      {
        if (Mdi0 (m, i, j) != Mdi0 (m, j, i)) return (0);
      }
    }
    return (1);
  }
  else
  {
    for (j = 0; j < nc; j++)
    {
      for (i = j + 1; i < nr; i++)
      {
        if (Mdr0 (m, i, j) != Mdr0 (m, j, i)) return (0);
      }
    }
    return (1);
  }
}

void
mdr_Detect_Inf (MDR * m)
{
  int n = MNR (m) * MNC (m);

  if(MD_TYPE_DOUBLE(m))
    if (detect_inf_r (MDRPTR (m), n))
      rerror ("matrix contains Inf value");
}

void
mdr_Detect_Nan (MDR * m)
{
  int n = MNR (m) * MNC (m);
  if(MD_TYPE_DOUBLE(m))
    if (detect_nan_r (MDRPTR (m), n))
      rerror ("matrix contains NaN value");
}

/* **************************************************************
 * Test a matrix for Infs..., integer cannot be Inf
 * ************************************************************** */

MDR *
mdr_IsInf (MDR * m)
{
  int i, size;
  MDR *mi;

  if(m->type == RLAB_TYPE_INT32)
  {
    mi = mdr_Create (MNR (m), MNC (m));
    mdr_Zero(mi);
  }
  else
  {
    size = MNR (m) * MNC (m);
    mi = mdr_Create (MNR (m), MNC (m));
    for (i = 0; i < size; i++)
    {
      if ((create_inf () == MdrV0 (m, i)) || (create_inf () == -MdrV0 (m, i)))
        MdrV0 (mi, i) = 1.0;
      else
        MdrV0 (mi, i) = 0.0;
    }
  }
  return (mi);
}

// ***************************************************************
// * Test a matrix for Nans: fixed
// ***************************************************************

MDR *
mdr_IsNan (MDR * m)
{
  int i, size;
  MDR *mi;

  mi = mdr_Create (MNR (m), MNC (m));

  if(m->type == RLAB_TYPE_INT32)
  {
    mdr_Zero(mi);
  }
  else
  {
    size = MNR (m) * MNC (m);
    for (i = 0; i < size; i++)
    {
      if (isnand (MdrV0 (m, i)))
        MdrV0 (mi, i) = 1.0;
      else
        MdrV0 (mi, i) = 0.0;
    }
  }
  return (mi);
}

/* **************************************************************
 * Create a matrix of Nans...
 * ************************************************************** */

MDR *
mdr_CreateNan (int nr, int nc)
{
  int i = 0;
  int size = nr * nc;
  MDR *mnan = mdr_Create (nr, nc);

  for (i = 0; i < size; i++)
  {
    MdrV0 (mnan, i) = create_nan ();
  }
  return (mnan);
}

// ***************************************************************
// * Create a matrix of Infs
// ***************************************************************

MDR *
mdr_CreateInf (int nr, int nc)
{
  int i = 0;
  int size = nr * nc;
  MDR *mnan = mdr_Create (nr, nc);
  for (i = 0; i < size; i++)
    MdrV0 (mnan, i) = create_inf ();
  return (mnan);
}


MDR *
mdr_Finite (MDR * m)
{
  int i, size;
  MDR *mi;

  mi = mdr_Create (MNR (m), MNC (m));

  if(m->type == RLAB_TYPE_INT32)
    mdr_Ones (mi);
  else
  {
    size = MNR (m) * MNC (m);
    for (i = 0; i < size; i++)
    {
      if ((create_inf () != MdrV0 (m, i)) &&	  /* Test against +inf */
           (-create_inf () != MdrV0 (m, i)) &&	/* Test against -inf */
           !isnand (MdrV0 (m, i)))              /* Test against  NaN */
        MdrV0 (mi, i) = 1.0;
      else
        MdrV0 (mi, i) = 0.0;
    }
  }
  return (mi);
}

int IsFinite (MDR * m)
{
  int i, size;
  int rval=1;
  double inf=create_inf();

  if(m->type == RLAB_TYPE_INT32)
    return rval;
  size = SIZE(m);
  for (i=0; i<size; i++)
  {
    if ( (inf == MdrV0 (m, i)) || (-inf == MdrV0 (m, i))
          || isnand (MdrV0 (m, i)) )
    {
      rval = 0;
      break;
    }
  }
  return (rval);
}


/* **************************************************************
 * Sin function...
 * ************************************************************** */

MDR *
mdr_Sin (MDR * m)
{
  int i, size;
  MDR *cm;

  cm = mdr_Create (MNR (m), MNC (m));
  size = MNR (m) * MNC (m);
  if(m->type == RLAB_TYPE_INT32)
  {
    for (i = 0; i < size; i++)
    {
      MdrV0 (cm, i) = errcheck (sin ( (double) MdiV0 (m, i) ), "sin");
    }
  }
  else
  {
    for (i = 0; i < size; i++)
    {
      MdrV0 (cm, i) = errcheck (sin (MdrV0 (m, i)), "sin");
    }
  }
  return (cm);
}

/* **************************************************************
 * Cos function...
 * ************************************************************** */

MDR *
mdr_Cos (MDR * m)
{
  int i, size;
  MDR *cm;

  cm = mdr_Create (MNR (m), MNC (m));
  size = MNR (m) * MNC (m);
  if(m->type == RLAB_TYPE_INT32)
  {
    for (i = 0; i < size; i++)
    {
      MdrV0 (cm, i) = errcheck (cos ( (double) MdiV0 (m, i)), "cos");
    }
  }
  else
  {
    for (i = 0; i < size; i++)
    {
      MdrV0 (cm, i) = errcheck (cos (MdrV0 (m, i)), "cos");
    }
  }
  return (cm);
}

/* **************************************************************
 * Tan function...
 * ************************************************************** */

MDR *
mdr_Tan (MDR * m)
{
  int i, size;
  MDR *cm;

  cm = mdr_Create (MNR (m), MNC (m));
  size = MNR (m) * MNC (m);
  if(m->type == RLAB_TYPE_INT32)
  {
    for (i = 0; i < size; i++)
    {
      MdrV0 (cm, i) = errcheck (tan ((double) MdiV0 (m, i)), "tan");
    }
  }
  else
  {
    for (i = 0; i < size; i++)
    {
      MdrV0 (cm, i) = errcheck (tan (MdrV0 (m, i)), "tan");
    }
  }
  return (cm);
}

/* **************************************************************
 * Arc-Sin function...
 * ************************************************************** */

void *
mdr_ASin (MDR * m, int *type)
{
  int i, size;
  MDR *rm;
  MDC *cm;

  if (do_cmplx_1 (m))
  {
    cm = mdc_Create (MNR (m), MNC (m));
    size = MNR (m) * MNC (m);
    for (i=0; i<size; i++)
      MdcV0 (cm, i) = casin ( mdrV0 (m, i) );
    *type = MATRIX_DENSE_COMPLEX;
    return (cm);
  }
  else
  {
    rm = mdr_Create (MNR (m), MNC (m));
    size = MNR (m) * MNC (m);
    for (i=0; i<size; i++)
      MdrV0 (rm, i) = errcheck (asin (mdrV0(m,i)), "asin");
    *type = MATRIX_DENSE_REAL;
    return (rm);
  }
}

/* **************************************************************
 * Arc-Cos function...
 * ************************************************************** */

void *
mdr_ACos (MDR * m, int *type)
{
  int i, size;
  MDR *rm;
  MDC *cm;

  if (do_cmplx_1 (m))
  {
    cm = mdc_Create (MNR (m), MNC (m));
    size = MNR (m) * MNC (m);
    for (i = 0; i < size; i++)
      MdcV0 (cm, i) = cacos (mdrV0 (m, i));
    *type = MATRIX_DENSE_COMPLEX;
    return (cm);
  }
  else
  {
    rm = mdr_Create (MNR (m), MNC (m));
    size = MNR (m) * MNC (m);
    for (i = 0; i < size; i++)
      MdrV0 (rm, i) = errcheck (acos (mdrV0 (m, i)), "acos");
    *type = MATRIX_DENSE_REAL;
    return (rm);
  }
}

/* **************************************************************
 * Arc-Tan function...
 * ************************************************************** */

void *
mdr_ATan (MDR * m, int *type)
{
  int i, size;
  MDR *rm;

  rm = mdr_Create (MNR (m), MNC (m));
  size = SIZE(m);

  for (i = 0; i < size; i++)
  {
    MdrV0 (rm, i) = errcheck (atan (mdrV0 (m, i)), "atan");
  }

  *type = MATRIX_DENSE_REAL;
  return (rm);
}

void *
mdr_ATan2 (MDR * m1, MDR * m2, int *type)
{
  MDR *rm=0;

  if ((m1 != 0)&&(m2 != 0))
  {
    MATRIXOP_MDR_DARGS2MDR(rm,m1,m2,atan2);
  }
  else if (m1 != 0)
  {
    if (MNC(m1)==2)
    {
      int r1 = MNR(m1);
      rm = mdr_Create(r1,1);
      int ir;
      for (ir=0;ir<r1;ir++)
      {
        MdrV0(rm,ir) = atan2(mdr0(m1,ir,0), mdr0(m1,ir,1));
      }
    }
  }

  *type = MATRIX_DENSE_REAL;
  return (rm);
}

/* **************************************************************
 * Sqrt function...
 * ************************************************************** */

void *
mdr_Sqrt (MDR * m, int *type)
{
  int i, size;
  MDR *rm;
  MDC *cm;

  if (do_cmplx (m))
  {
    cm = mdc_Create (MNR (m), MNC (m));
    size = MNR (m) * MNC (m);
    for (i = 0; i < size; i++)
      MdcV0 (cm, i) = csqrt ( mdrV0 (m, i) );
    *type = MATRIX_DENSE_COMPLEX;
    return (cm);
  }
  else
  {
    rm = mdr_Create (MNR (m), MNC (m));
    size = MNR (m) * MNC (m);
    for (i = 0; i < size; i++)
      MdrV0 (rm, i) = errcheck (sqrt (mdrV0 (m, i)), "sqrt");
    *type = MATRIX_DENSE_REAL;
    return (rm);
  }
}

/* **************************************************************
 * Log function...
 * ************************************************************** */

void *
mdr_Log (MDR * m, int *type)
{
  int i, size;
  MDR *rm;
  MDC *cm;

  if (do_cmplx (m))
  {
    cm = mdc_Create (MNR (m), MNC (m));
    size = MNR (m) * MNC (m);
    for (i=0; i < size; i++)
    {
      MdcV0 (cm, i) = clog(mdrV0 (m, i));
    }
    *type = MATRIX_DENSE_COMPLEX;
    return (cm);
  }
  else
  {
    rm = mdr_Create (MNR (m), MNC (m));
    size = MNR (m) * MNC (m);
    for (i = 0; i < size; i++)
    {
      if (mdrV0 (m, i) == 0)
        MdrV0 (rm, i) = -create_inf ();
      else
        MdrV0 (rm, i) = errcheck (log (mdrV0 (m, i)), "log");
    }
    *type = MATRIX_DENSE_REAL;
    return (rm);
  }
}

/* **************************************************************
 * Log10 function...
 * ************************************************************** */

#define log10e 0.43429448190325182765

void *
mdr_Log10 (MDR * m, int *type)
{
  int i, size;
  MDR *rm;
  MDC *cm;

  if (do_cmplx (m))
  {
    cm = mdc_Create (MNR (m), MNC (m));
    size = MNR (m) * MNC (m);
    for (i = 0; i < size; i++)
    {
      if (MdiV0 (m, i) == 0)
      {
        MdcV0r (cm, i) = -create_inf ();
        MdcV0i (cm, i) = 0.0;
      }
      else
      {
        MdcV0 (cm, i) = log10e * clog (mdrV0 (m, i));
      }
    }
    *type = MATRIX_DENSE_COMPLEX;
    return (cm);
  }
  else
  {
    rm = mdr_Create (MNR (m), MNC (m));
    size = MNR (m) * MNC (m);
    for (i = 0; i < size; i++)
    {
      if (mdrV0 (m, i) == 0)
        MdrV0 (rm, i) = -create_inf ();
      else
        MdrV0 (rm, i) = errcheck (log10 ((double) mdrV0 (m, i)), "log10");
    }
    *type = MATRIX_DENSE_REAL;
    return (rm);
  }
}

/* **************************************************************
 * Exp function...
 * ************************************************************** */
MDR *
mdr_Exp (MDR * m)
{
  int i, size;
  MDR *rm;

  rm = mdr_Create (MNR (m), MNC (m));
  size = MNR (m) * MNC (m);
  for (i = 0; i < size; i++)
    MdrV0 (rm, i) = errcheck (exp (mdrV0 (m, i)), "exp");

  return (rm);
}

MDR *
mdr_Exp_weight (MDR * m, MDR *w)
{
  int i, j, nr, nc, iw, jw;
  MDR *new=0;

  if (w)
  {
    iw = MNR(w) - 1;
    jw = MNC(w) - 1;
  }

  nr = MNR (m);
  nc = MNC (m);

  if (w)
  {
    new = mdr_Create (nr, nc);

    for(i=0;i<nr; i++) for(j=0;j<nc; j++)
    {
      Mdr0(new,i,j) = mdr0(w,MIN(i,iw),MIN(j,jw)) / (mdr0(m,i,j) * mdr0(m,i,j));
    }
  }

  return (new);
}

/* **************************************************************
 * Diag() function...
 * ************************************************************** */

MDR *
mdr_Diag (MDR * marg, int k)
{
  int i, smin, size;
  MDR *m;

  if (MNR (marg) == 1 || MNC (marg) == 1)
  {
    /* Create a diagonal matrix */
    size = MAX(MNR (marg), MNC (marg)) + ABS(k);
    if(marg->type == RLAB_TYPE_INT32)
    {
      m = mdi_Create (size, size);
      mdr_Zero (m);
      if (k < 0)
      {
        for (i = 1 - k; i <= size; i++)
          Mdi1 (m, i, (i + k)) = MdiV1 (marg, (i + k));
      }
      else
      {
        for (i = 1; i <= size - k; i++)
          Mdi1 (m, i, (i + k)) = MdiV1 (marg, i);
      }
    }
    else
    {
      m = mdr_Create (size, size);
      mdr_Zero (m);
      if (k < 0)
      {
        for (i = 1 - k; i <= size; i++)
          Mdr1 (m, i, (i + k)) = MdrV1 (marg, (i + k));
      }
      else
      {
        for (i = 1; i <= size - k; i++)
          Mdr1 (m, i, (i + k)) = MdrV1 (marg, i);
      }
    }
  }
  else
  {
    /* Extract the diagonal elements */
    smin = MIN(MNR (marg), MNC (marg));

    if (k >= 0)
    {
      if (MNR (marg) >= MNC (marg))
        size = smin - k;
      else
        size = smin - MAX(0, k - (MNC (marg) - MNR (marg)));
    }
    else
    {
      if (MNR (marg) >= MNC (marg))
        size = smin - MAX(0, -k - (MNR (marg) - MNC (marg)));
      else
        size = smin + k;
    }
    if (size <= 0)
      size = 0;

    if(marg->type == RLAB_TYPE_INT32)
    {
      m = mdi_Create (size, 1);
      mdr_Zero (m);
      if (k >= 0)
      {
        for (i = 1; i <= size; i++)
          Mdi1 (m, i, 1) = Mdi1 (marg, i, i + k);
      }
      else
      {
        for (i = 1; i <= size; i++)
          Mdi1 (m, i, 1) = Mdi1 (marg, i - k, i);
      }
    }
    else
    {
      m = mdr_Create (size, 1);
      mdr_Zero (m);
      if (k >= 0)
      {
        for (i = 1; i <= size; i++)
          Mdr1 (m, i, 1) = Mdr1 (marg, i, i + k);
      }
      else
      {
        for (i = 1; i <= size; i++)
          Mdr1 (m, i, 1) = Mdr1 (marg, i - k, i);
      }
    }
  }

  return (m);
}

/* **************************************************************
 * Check the values of a MATRIX. We are looking to decide whether
 * or not we need to do a complex op on the matrix, or we can get
 * away with a real operation. Return TRUE (1) at the 1st sign of
 * an element < 0.0
 * ************************************************************** */

int
do_cmplx (MDR * m)
{
  int i, size;

  size = MNR (m) * MNC (m);
  for (i=0; i<size; i++)
  {
    if (mdrV0 (m, i) < 0.0)
    {
      return (1);
    }
  }
  return (0);
}

int
do_cmplx_1 (MDR * m)
{
  int i, size;

  size = (m->nrow) * (m->ncol);

  for (i=0; i<size; i++)
  {
    if (mdrV0(m, i) > 1.0 || mdrV0(m, i) < -1.0)
    {
      return (1);
    }
  }
  return (0);
}

/* **************************************************************
 * Return the real part of the matrix.
 * ************************************************************** */

MDR *
mdr_Real_BF (MDR * m)
{
  return (m);
}

/* **************************************************************
 * Return the imaginary part of a matrix.
 * ************************************************************** */

MDR *
mdr_Imag_BF (MDR * m)
{
  MDR *mz = mdr_Create (MNR (m), MNC (m));
  mdr_Zero (mz);

  return (mz);
}

/* **************************************************************
 * Return the conjugate of a matrix.
 * ************************************************************** */

MDR *
mdr_Conj_BF (MDR * m)
{
  return (m);
}

/* **************************************************************
 * Find the indices of non-zero elements in a matrix.
 * ************************************************************** */

MDR *
mdr_Find_BF (MDR * m)
{
  int i, j, size;
  double *dtmp;
  MDR *mf;

  /* Make mfind as big as it could be, reduce later */
  mf = mdr_Create (1, SIZE(m));

  j = 0;
  size = SIZE(m);
  for (i=0; i<size; i++)
  {
    MdrV0 (mf, j) = 0;
    if (mdrV0 (m, i) != 0)
      MdrV0 (mf, j++) = i + 1;
  }

  /* Now reduce mfind to correct size */
  if (j == 0)
  {
    /* No reduction here, just make dimensions correct */
    mf->nrow = 0;
    mf->ncol = 0;
  }
  else if (j < size)
  {
    dtmp = (double *) GC_malloc_atomic_ignore_off_page (sizeof (double) * j);

    memcpy (dtmp, MDPTR(mf), sizeof (double) * j);
    GC_FREE (MDPTR(mf));
    MDPTR(mf) = (void *) dtmp;
    mf->ncol = j;
  }

  return (mf);
}

/* **************************************************************
 * Sort a vector, or matrix (builtin-function).
 * ************************************************************** */

/*
 * Create an integer MxN matrix for Sort()
 */

MDR *
mdr_CreateFillSind (int nrow, int ncol)
{
  int i, j;
  MDR *new = mdr_Create (nrow, ncol);

  for (i = 1; i <= nrow; i++)
  {
    for (j = 1; j <= ncol; j++)
    {
      Mdr1 (new, i, j) = (double) i;
    }
  }

  return (new);
}


Btree *
mdr_Sort_BF (MDR * m)
{
  int i, n;
  Btree *bt;
  MDR *mcopy=0, *sind=0;

  if (EQNULL(m))
  {
    sind  = mdr_Create (0,0);
    mcopy = mdr_Create (0,0);
  }
  else if (EQVECT(m))
  {
    /* Vector sort */
    n = MAX(MNR (m), MNC (m));
    sind = mdr_CreateFillSind (n, 1);
    sind->ncol = sind->nrow;
    sind->nrow = 1;
    mcopy = mdr_Float_BF (m);
    r_sort ((double *) MDRPTR (mcopy), 0, n - 1, (double *) MDRPTR (sind));
  }
  else
  {
    /* Matrix sort (column-wise) */
    n = MNR (m);
    sind = mdr_CreateFillSind (MNR (m), MNC (m));
    mcopy = mdr_Float_BF (m);
    for (i = 0; i < MNC (m); i++)
    {
      r_sort (&MdrV0(mcopy,i*n), 0, n - 1, &MdrV0(sind,i*n));
    }
  }

  bt = btree_Create ();
  install (bt, RLAB_NAME_GEN_INDEX, ent_Assign_Rlab_MDR(sind));
  install (bt, RLAB_NAME_GEN_VALUE, ent_Assign_Rlab_MDR(mcopy));
  return (bt);
}


/* **************************************************************
 * Sign function.
 * ************************************************************** */

MDR *
mdr_Sign_BF (MDR * m)
{
  int i, size;
  MDR *sm;

  size = MNR (m) * MNC (m);

  if(m->type == RLAB_TYPE_INT32)
  {
    sm = mdi_Create (MNR (m), MNC (m));
    for (i = 0; i < size; i++)
    {
      if (MdiV0 (m, i) > 0)
      {
        MdiV0 (sm, i) = 1;
      }
      else if (MdiV0 (m, i) < 0)
      {
        MdiV0 (sm, i) = -1;
      }
      else
      {
        MdiV0 (sm, i) = 0;
      }
    }
  }
  else
  {
    sm = mdr_Create (MNR (m), MNC (m));
    for (i = 0; i < size; i++)
    {
      if (MdrV0 (m, i) > 0)
      {
        MdrV0 (sm, i) = 1.0;
      }
      else if (MdrV0 (m, i) < 0)
      {
        MdrV0 (sm, i) = -1.0;
      }
      else
      {
        MdrV0 (sm, i) = 0.0;
      }
    }
  }

  return (sm);
}

/* **************************************************************
 * Increment (++) a matrix. Do the operation in place.
 * ************************************************************** */

MDR *
mdr_Increment (MDR * m)
{
  int i, size;

  size = SIZE(m);
  if(MD_TYPE_INT32(m))
  {
    for (i = 0; i < size; i++)
    {
      MdiV0 (m, i) += 1;
    }
  }
  else
  {
    for (i = 0; i < size; i++)
    {
      MdrV0 (m, i) += 1.0;
    }
  }

  return (m);
}

MDR *
mdr_Decrement (MDR * m)
{
  int i, size;
  size = SIZE(m);
  if(MD_TYPE_INT32(m))
  {
    for (i = 0; i < size; i++)
    {
      MdiV0 (m, i) -= 1;
    }
  }
  else
  {
    for (i = 0; i < size; i++)
    {
      MdrV0 (m, i) -= 1.0;
    }
  }

  return (m);
}

size_t
mdr_Sizeof (MDR * m)
{
  int size = SIZE(m);
  if(m->type == RLAB_TYPE_INT32)
    return (size_t) (sizeof (int) * size);
  else
    return (size_t) (sizeof (double) * size);
}

Btree *
mdr_Frexp_BF (MDR * m)
{
  Btree *bt;
  MDR *e, *f;
  int i, ie;

  /* Only allowable solution. */

  f = mdr_Create (MNR (m), MNC (m));
  e = mdr_Create (MNR (m), MNC (m));

  for (i=0; i<SIZE(m); i++)
  {
    MdrV0 (f, i) = errcheck (rfrexp (mdrV0 (m,i), &ie), "frexp");
    MdrV0 (e, i) = (double) ie;
  }

  bt = btree_Create ();
  install (bt, "e", ent_Assign_Rlab_MDR(e));
  install (bt, "f", ent_Assign_Rlab_MDR(f));
  return (bt);
}

MDR *
mdr_Ldexp_BF (MDR * F, MDR * E)
{
  MDR *ans;
  int i, size;

  /* Special case, e is a scalar. */
  if (EQSCAL(E))
  {
    ans = mdr_Create (MNR (F), MdrNC (F));
    size = SIZE(F);

    for (i=0; i<size; i++)
    {
      MdrV0 (ans,i) = errcheck (ldexp (mdrV0 (F,i), mdiV0 (E,0)),	"ldexp");
    }

    return (ans);
  }

  /* Check F and E sizes for compatibility. */
  if (!EQSIZE(E,F))
  {
    fprintf (stderr, "incompatible matrix dimensions\n");
    fprintf (stderr, "F-nrow = %i, E-nrow = %i\n", MNR (F), MNR (E));
    fprintf (stderr, "F-ncol = %i, E-ncol = %i\n", MNC (F), MNC (E));
    rerror ("ldexp: invalid argument dimensions");
  }

  ans = mdr_Create (MNR (F), MdrNC (F));
  size = SIZE(F);

  for (i = 0; i < size; i++)
  {
    MdrV0 (ans,i) = errcheck (ldexp (mdrV0 (F,i), mdiV0 (E,i)), "ldexp");
  }

  return (ans);
}

#ifdef HAVE_LOGB
MDR *
mdr_Logb_BF (MDR * m)
{
  MDR *lb;
  int i, size;

  lb = mdr_Create (MNR (m), MNC (m));
  size = SIZE(m);

  for (i=0; i<size; i++)
  {
    MdrV0 (lb, i) = errcheck (logb(mdrV0 (m, i)), "logb");
  }

  return (lb);
}
#endif /* HAVE_LOGB */

/* **************************************************************
 * Convert a dense real matrix to a dense real matrix. This is
 * a no-op, just return the original matrix.
 * ************************************************************** */

MDR *
mdr_Dense (MDR * m)
{
  return (m);
}


//
//
// B I N A R Y   O P E R A T I O N S
//
//
static unsigned char c[8];
static unsigned char d[8];
static int convert_cbytes_to_int(int l, int use_lsb)
{
  int i;
  int rval = 0;

  if (use_lsb)
  {
    for (i=0; i<l; i++)
    {
      rval <<= 8;
      rval |= c[i];
    }
  }
  else
  {
    // use msb
    for (i=l-1; i>=0; i--)
    {
      rval <<= 8;
      rval |= c[i];
    }
  }

  return rval;
}
static void convert_int_to_cbytes(int b, int l, int use_lsb)
{
  int i;

  // convert a to l-bytes
  if (use_lsb)
  {
    for (i=l-1; i>=0; i--)
    {
      c[i] = b & 0xff;
      b >>= 8;
    }
  }
  else
  {
    for (i=0; i<l; i++)
    {
      c[i] = b & 0xff;
      b >>= 8;
    }
  }
}
static void convert_byte_to_dbits(int a)
{
  int i;
  a = a & 0xff;

  for(i=0; i<8; i++)
  {
    d[i] = ((0x01<<i) & a) == 0x00 ? 0:1;
  }
}
static int convert_dbits_to_byte( void )
{
  int i;
  int rval = 0;
  for(i=7; i>=0; i++)
  {
    rval <<= 1;
    rval |= d[i];
  }
  return rval;
}


//
//
//
static int int_bit_shift_right(int a, int k, int l)
{
  int i, carrydown, m, carrydown_new;
  int rval=0;

  if (k>l*8)
    return rval;

  // convert a to l-bytes
  convert_int_to_cbytes(a,l,0);

  // shift each byte to left. start from zero, have a carry over
  // use MSB
  for (m=0; m<k; m++)
  {
    carrydown = 0;
    for (i=l-1; i>=0; i--)
    {
      carrydown_new = ( (c[i] & 0x01) != 0x00 );
      c[i] >>= 1;
      c[i] |= carrydown;
      carrydown = carrydown_new == 0x01 ? 0x80 : 0x00;
    }
  }

  // l-bytes to a
  return convert_cbytes_to_int(l,0);
}
MDR * mdi_bit_shift_right (MDR * a, MDR * k, MDR *s)
{
  MDR *ans=0;

  if(!MD_TYPE_INT32(a))
    return ans;

  MATRIXOP_MDI_IARGS_MDI_2MD(ans,a,k,s,int_bit_shift_right);

  return (ans);
}

static int int_bit_shift_left(int a, int k, int l)
{
  int i, carryover, m, carryover_new;

  if (k>l*8)
    return 0;

  // convert a to l-bytes
  convert_int_to_cbytes(a,l,0);

  // shift each byte to left. start from zero, have a carry over
  for (m=0; m<k; m++)
  {
    carryover = 0;
    for (i=0; i<l; i++)
    {
      carryover_new = ( (c[i] & 0x80) != 0x00 );
      c[i] <<= 1;
      c[i] |= carryover;
      carryover = carryover_new;
    }
  }

  // l-bytes to a
  return convert_cbytes_to_int(l,0);
}
MDR * mdi_bit_shift_left (MDR * a, MDR * k, MDR *s)
{
  MDR *ans=0;

  if(!MD_TYPE_INT32(a))
    return ans;

  MATRIXOP_MDI_IARGS_MDI_2MD(ans,a,k,s,int_bit_shift_left);

  return (ans);
}

MDR * mdi_byte_split (int b, int use_lsb)
{
  MDR *res=0;

  int i;

  // convert a to l-bytes
  convert_int_to_cbytes(b,4,0);

  res = mdi_Create(1,4);
  if (use_lsb)
  {
    for (i=0; i<4; i++)
      MdiV0(res,i) = c[i];
  }
  else
  {
    for (i=0; i<4; i++)
      MdiV0(res,4-i-1) = c[i];
  }

  return (res);
}

MDR * mdi_bit_split (MDR * a, int l, int use_lsb)
{
  MDR *res=0;

  int i, j, k, nr=MNR(a);

  res = mdi_Create(nr,l*8);

  for (k=0; k<nr; k++)
  {
    int b = mdiV0(a,k);

    // convert a to l-bytes
    convert_int_to_cbytes(b,l,use_lsb);

    for (i=0; i<l; i++)
    {
      convert_byte_to_dbits(c[i]);
      for (j=0;j<8;j++)
        Mdi0(res,k,(l-i-1)*8+j) = d[7-j];
    }
  }

  return (res);
}

MDR * mdi_byte_join(MDR * a, int use_lsb)
{
  MDR *rval=0;
  int i,res=0,n=SIZE(a);

  if (use_lsb)
  {
    for (i=n-1; i>=0; i--)
    {
      res <<= 8;
      res |= (MdiV0(a,i) & 0xff);
    }
  }
  else
  {
    for (i=0; i<n; i++)
    {
      res <<= 8;
      res |= (MdiV0(a,i) & 0xff);
    }
  }

  rval = mdi_CreateScalar(res);
  return rval;

}

/*
   Derivation from the fortran version of CONREC by Paul Bourke
   d               ! matrix of data to contour
   x               ! data matrix column coordinates
   y               ! data matrix row coordinates
   z               ! contour levels in increasing order
*/

void mdr_Contour(MDR * d, MDR *x, MDR *y, MDR *z, MDR **r)
{
#define xsect(p1,p2) (h[p2]*xh[p1]-h[p1]*xh[p2])/(h[p2]-h[p1])
#define ysect(p1,p2) (h[p2]*yh[p1]-h[p1]*yh[p2])/(h[p2]-h[p1])
   int m1,m2,m3,case_value;
   double dmin,dmax,x1=0,x2=0,y1=0,y2=0;
   int i,j,k,m;
   double h[5];
   int sh[5];
   double xh[5],yh[5];
   int im[4] = {0,1,1,0},jm[4]={0,0,1,1};
   int castab[3][3][3] =
   {
     { {0,0,8},{0,2,5},{7,6,9} },
     { {0,3,4},{1,3,1},{4,3,0} },
     { {9,6,7},{5,2,0},{8,0,0} }
   };
   double temp1,temp2;

   int nc = SIZE(z);
   int *row_counter = GC_MALLOC(nc * sizeof(int));

   for (k=0;k<nc;k++)
     row_counter[k]=0;

   int jub = MNR(d);
   int iub = MNC(d);

   for (j=(jub-2);j>=0;j--)
   {
      for (i=0;i<(iub-1);i++)
      {
         temp1 = MIN(mdr0(d,j, i  ),mdr0(d,j+1,  i));
         temp2 = MIN(mdr0(d,j, i+1),mdr0(d,j+1,i+1));
         dmin  = MIN(temp1,temp2);
         temp1 = MAX(mdr0(d,j, i  ),mdr0(d,j+1,  i));
         temp2 = MAX(mdr0(d,j, i+1),mdr0(d,j+1,i+1));
         dmax  = MAX(temp1,temp2);
         if (dmax < mdrV0(z,0) || dmin > mdrV1(z,nc))
            continue;
         for (k=0;k<nc;k++)
         {
            if ((mdrV0(z,k) < dmin) || (mdrV0(z,k) > dmax))
               continue;
            for (m=4;m>=0;m--)
            {
               if (m > 0)
               {
                 h[m]  = mdr0(d,j+jm[m-1],i+im[m-1])-mdrV0(z,k);
                 if (x)
                   xh[m] = mdrV0(x,i+im[m-1]);
                 else
                   xh[m] = (double) (i+im[m-1]) + 1.0;

                 if (y)
                   yh[m] = mdrV0(y,j+jm[m-1]);
                 else
                   yh[m] = (double) (j+jm[m-1]) + 1.0;
               }
               else
               {
                  h[0]  = 0.25 * (h[1]+h[2]+h[3]+h[4]);

                  if (x)
                    xh[0] = 0.50 * (mdrV0(x,i)+mdrV0(x,i+1));
                  else
                    xh[0] = (double) i + 1.50;

                  if (y)
                    yh[0] = 0.50 * (mdrV0(y,j)+mdrV0(y,j+1));
                  else
                    yh[0] = (double) j + 1.50;
               }

               if (h[m] > 0.0)
                  sh[m] = 1;
               else if (h[m] < 0.0)
                  sh[m] = -1;
               else
                  sh[m] = 0;
            }

            /*
               Note: at this stage the relative heights of the corners and the
               centre are in the h array, and the corresponding coordinates are
               in the xh and yh arrays. The centre of the box is indexed by 0
               and the 4 corners by 1 to 4 as shown below.
               Each triangle is then indexed by the parameter m, and the 3
               vertices of each triangle are indexed by parameters m1,m2,and m3.
               It is assumed that the centre of the box is always vertex 2
               though this isimportant only when all 3 vertices lie exactly on
               the same contour level, in which case only the side of the box
               is drawn.
                  vertex 4 +-------------------+ vertex 3
                           | \               / |
                           |   \    m-3    /   |
                           |     \       /     |
                           |       \   /       |
                           |  m=2    X   m=2   |       the centre is vertex 0
                           |       /   \       |
                           |     /       \     |
                           |   /    m=1    \   |
                           | /               \ |
                  vertex 1 +-------------------+ vertex 2
            */
            /* Scan each triangle in the box */
            for (m=1;m<=4;m++)
            {
               m1 = m;
               m2 = 0;
               if (m != 4)
                  m3 = m + 1;
               else
                  m3 = 1;
               if ((case_value = castab[sh[m1]+1][sh[m2]+1][sh[m3]+1]) == 0)
                  continue;
               switch (case_value)
               {
                  case 1: /* Line between vertices 1 and 2 */
                   x1 = xh[m1];
                   y1 = yh[m1];
                   x2 = xh[m2];
                   y2 = yh[m2];
                   break;
               case 2: /* Line between vertices 2 and 3 */
                   x1 = xh[m2];
                   y1 = yh[m2];
                   x2 = xh[m3];
                   y2 = yh[m3];
                   break;
               case 3: /* Line between vertices 3 and 1 */
                   x1 = xh[m3];
                   y1 = yh[m3];
                   x2 = xh[m1];
                   y2 = yh[m1];
                   break;
               case 4: /* Line between vertex 1 and side 2-3 */
                   x1 = xh[m1];
                   y1 = yh[m1];
                   x2 = xsect(m2,m3);
                   y2 = ysect(m2,m3);
                   break;
               case 5: /* Line between vertex 2 and side 3-1 */
                   x1 = xh[m2];
                   y1 = yh[m2];
                   x2 = xsect(m3,m1);
                   y2 = ysect(m3,m1);
                   break;
               case 6: /* Line between vertex 3 and side 1-2 */
                   x1 = xh[m3];
                   y1 = yh[m3];
                   x2 = xsect(m1,m2);
                   y2 = ysect(m1,m2);
                   break;
               case 7: /* Line between sides 1-2 and 2-3 */
                   x1 = xsect(m1,m2);
                   y1 = ysect(m1,m2);
                   x2 = xsect(m2,m3);
                   y2 = ysect(m2,m3);
                   break;
               case 8: /* Line between sides 2-3 and 3-1 */
                   x1 = xsect(m2,m3);
                   y1 = ysect(m2,m3);
                   x2 = xsect(m3,m1);
                   y2 = ysect(m3,m1);
                   break;
               case 9: /* Line between sides 3-1 and 1-2 */
                   x1 = xsect(m3,m1);
                   y1 = ysect(m3,m1);
                   x2 = xsect(m1,m2);
                   y2 = ysect(m1,m2);
                   break;
               default:
                   break;
               }

               /* Finally draw the line  -> ConrecLine(x1,y1,x2,y2,z[k]) */
//                fprintf(stderr, "x1,y1 = %g,%g\n", x1,y1);
//                fprintf(stderr, "x2,y2 = %g,%g\n", x2,y2);
               // add entry to r[k] at the position row_counter[k]
               if (row_counter[k] <= MNR(r[k]))
                 mdr_Extend(r[k],MNR(r[k])+20,2);
               Mdr0(r[k],row_counter[k],0) = x1;
               Mdr0(r[k],row_counter[k],1) = y1;
               row_counter[k]++;
               Mdr0(r[k],row_counter[k],0) = x2;
               Mdr0(r[k],row_counter[k],1) = y2;
               row_counter[k]++;
            } /* m */

         } /* k - contour */

      } /* i */

   } /* j */

  // update sizes of return matrices
  for (k=0;k<nc;k++)
  {
    if (row_counter[k] < MNR(r[k]))
      mdr_Extend(r[k], row_counter[k], MNC(r[k]));
  }
  GC_FREE(row_counter);
}




