/* mdcf1.c Matrix Dense Complex Functions */

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
#include "mdc.h"
#include "mdr.h"
#include "mds.h"
#include "btree.h"
#include "util.h"
#include "complex.h"
#include "symbol.h"
#include "mathl.h"
#include "sort.h"
#include "rlab_solver_parameters_names.h"

#include <stdio.h>
#include <math.h>


/* **************************************************************
 * Convert a dense complex matrix to a dense complex matrix. This is
 * a no-op, just return the original matrix.
 * ************************************************************** */

MDC *
mdc_Dense (MDC * m)
{
  return (m);
}


/* **************************************************************
 * Return matrix-dense-complex member references.
 * ************************************************************** */

void *
mdc_MemberRef (MDC * m, char *name, int *type)
{
  Ent *ne;
  void *rptr;

  if (!strcmp (name, RLAB_MEMBER_NROW))
  {
    ne = ent_Create ();
    ent_data (ne) = mdr_CreateScalar ((double) MNR (m));
    ent_SetType (ne, MATRIX_DENSE_REAL);
    rptr = (void *) ne;
    *type = ENTITY;
  }
  else if (!strcmp (name, RLAB_MEMBER_NCOL))
  {
    ne = ent_Create ();
    ent_data (ne) = mdr_CreateScalar ((double) MNC (m));
    ent_SetType (ne, MATRIX_DENSE_REAL);
    rptr = (void *) ne;
    *type = ENTITY;
  }
  else if (!strcmp (name, RLAB_MEMBER_SIZE))
  {
    ne = ent_Create ();
    ent_data (ne) = mdr_CreateScalar ((double) (MNR (m) * MNC (m)));
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
    ent_data (ne) = mds_CreateScalar ( RLAB_STORAGE_DENSE );
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
mdc_Members (MDC * m, int *n)
{
  char **marray;

  marray = (char **) GC_MALLOC (6 * sizeof (char *));
  if (marray == 0)
    rerror ("out of memory");

  marray[0] = cpstr (RLAB_MEMBER_NROW);
  marray[1] = cpstr (RLAB_MEMBER_NCOL);
  marray[2] = cpstr (RLAB_MEMBER_SIZE);
  marray[3] = cpstr (RLAB_MEMBER_CLASS);
  marray[4] = cpstr (RLAB_MEMBER_TYPE);
  marray[5] = cpstr (RLAB_MEMBER_STORAGE);

  *n = 6;
  return (marray);
}

/* **************************************************************
 * Return the type of a matrix...
 * ************************************************************** */

MDS *
mdc_Type_BF (MDC * m)
{
  MDS *type = mds_CreateScalar (RLAB_MEMBER_TYPE_COMPLEX);
  return (type);
}

/* **************************************************************
 * Int function...
 * ************************************************************** */

MDC *
mdc_Int_BF (MDC * m)
{
  int i, size;
  MDC *im;

  im = mdc_Create (MNR (m), MNC (m));
  size = MNR (m) * MNC (m);
  for (i = 0; i < size; i++)
  {
    MdcV0r (im, i) = (double) ((int) MdcV0r (m, i));
    MdcV0i (im, i) = (double) ((int) MdcV0i (m, i));
  }
  return (im);
}

/* **************************************************************
 * Ceil function...
 * ************************************************************** */
MDC * mdc_CeilFlooRound_BF (int is_ceil, MDC * m, MDR *b, MD *o)
{
  int i, j;
  MDR *bin=0;
  MDC *offs=0, *cm=0;
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

  cm = mdc_Create (rr, rc);

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

  for (i=0; i<MNR(cm); i++)
  {
    for (j=0; j<MNC(cm); j++)
    {
      Complex x1 = mdc0(m,MIN(i,r1-1),MIN(j,c1-1));
      Complex o1 = mdc0(offs,MIN(i,r3-1),MIN(j,c3-1));
      double  b1 = Mdr0(bin, MIN(i,r2-1),MIN(j,c2-1));
      Complex d1 = (x1-o1)/b1;
      RE(d1) = ceilflooround(RE(d1));
      IM(d1) = ceilflooround(IM(d1));
      Mdc0 (cm,i,j) = o1 + b1 * d1;
    }
  }

  if (!b)
    mdr_Destroy(bin);
  if (!o)
    mdr_Destroy(offs);

  return (cm);
}

MDC *
mdc_Floor_BF (MDC * m)
{
  int i, size;
  MDC *cm;

  cm = mdc_Create (MNR (m), MNC (m));
  size = MNR (m) * MNC (m);
  for (i = 0; i < size; i++)
  {
    MdcV0r (cm, i) = errcheck (floor (MdcV0r (m, i)), "floor");
    MdcV0i (cm, i) = errcheck (floor (MdcV0i (m, i)), "floor");
  }
  return (cm);
}

/* **************************************************************
 * Round (rint) function...
 * ************************************************************** */

MDC *
mdc_Round_BF (MDC * m)
{
  int i, size;
  MDC *cm;

  cm = mdc_Create (MNR (m), MNC (m));
  size = MNR (m) * MNC (m);
  for (i = 0; i < size; i++)
  {
    MdcV0r (cm, i) = errcheck (rint (MdcV0r (m, i)), "round");
    MdcV0i (cm, i) = errcheck (rint (MdcV0i (m, i)), "round");
  }
  return (cm);
}

/* **************************************************************
 * Abs function...
 * ************************************************************** */

MDR *
mdc_Abs (MDC * m)
{
  int i, size;
  MDR *new;

  new = mdr_Create (MNR (m), MNC (m));
  size = MNR (m) * MNC (m);

  for (i = 0; i < size; i++)
    MdrV0 (new, i) = cabs(MdcV0 (m, i));

  return (new);
}

/* **************************************************************
 * Max-Index function...
 * ************************************************************** */

MDR * mdc_MaxI1 (MDC * m)
{
  Complex maxc;
  MDR *mmax = 0;
  int i, j, ind;

  if (MNR (m) == 1)
  {
    maxc = MdcV0 (m, 0);
    ind = 1;
    for (i = 1; i < MNC (m); i++)
    {
      if (cabs (MdcV0 (m, i)) > cabs (maxc))
      {
        maxc = MdcV0 (m, i);
        ind = i + 1;
      }
    }
    mmax = mdr_Create (1, 1);
    Mdr0 (mmax, 0, 0) = (double) ind;
  }
  else
  {
    mmax = mdr_Create (1, MNC (m));
    for (i = 0; i < MNC (m); i++)
    {
      maxc = Mdc0 (m, 0, i);
      ind = 1;
      for (j = 1; j < MNR (m); j++)
        if (cabs (Mdc0 (m, j, i)) > cabs (maxc))
      {
        maxc = Mdc0 (m, j, i);
        ind = j + 1;
      }
      Mdr0 (mmax, 0, i) = (double) ind;
    }
  }
  return (mmax);
}

/* **************************************************************
 * Max function with two arguments...
 * ************************************************************** */

MDC *
mdc_Max2 (MDC * m1, MDC * m2)
{
  int i, size;
  MDC *m = 0;

  /* Check sizes */

  if (MNR (m1) == MNR (m2) && MNC (m1) == MNC (m2))
  {
    m = mdc_Copy (m1);
    size = MNR (m1) * MNC (m1);
    for (i = 0; i < size; i++)
    {
      if (cabs (MdcV0 (m1, i)) < cabs (MdcV0 (m2, i)))
	MdcV0 (m, i) = MdcV0 (m2, i);
    }
  }
  else
    rerror ("max: matrix sizes must be the same");

  return (m);
}

/* **************************************************************
 * Min function...
 * ************************************************************** */

MDC *
mdc_Min1 (MDC * m)
{
  Complex maxc;
  MDC *mmax = 0;
  int i, j;

  if (MNR (m) == 1)
  {
    maxc = MdcV0 (m, 0);
    for (i = 1; i < MNC (m); i++)
    {
      if (cabs (MdcV0 (m, i)) < cabs (maxc))
	maxc = MdcV0 (m, i);
    }
    mmax = mdc_Create (1, 1);
    Mdc0 (mmax, 0, 0) = maxc;
  }
  else
  {
    mmax = mdc_Create (1, MNC (m));
    for (i = 0; i < MNC (m); i++)
    {
      maxc = Mdc0 (m, 0, i);
      for (j = 1; j < MNR (m); j++)
        if (cabs (Mdc0 (m, j, i)) < cabs (maxc))
          maxc = Mdc0 (m, j, i);
      Mdc0 (mmax, 0, i) = maxc;
    }
  }

  return (mmax);
}

/* **************************************************************
 * Min-Index function...
 * ************************************************************** */

MDR *
mdc_MinI1 (MDC * m)
{
  Complex maxc;
  MDR *mmax = 0;
  int i, j, ind;

  if (MNR (m) == 1)
  {
    maxc = MdcV0 (m, 0);
    ind = 1;
    for (i = 1; i < MNC (m); i++)
    {
      if (cabs (MdcV0 (m, i)) < cabs (maxc))
      {
        maxc = MdcV0 (m, i);
        ind = i + 1;
      }
    }
    mmax = mdr_Create (1, 1);
    Mdr0 (mmax, 0, 0) = (double) ind;
  }
  else
  {
    mmax = mdr_Create (1, MNC (m));
    for (i = 0; i < MNC (m); i++)
    {
      maxc = Mdc0 (m, 0, i);
      ind = 1;
      for (j = 1; j < MNR (m); j++)
        if (cabs (Mdc0 (m, j, i)) < cabs (maxc))
      {
        maxc = Mdc0 (m, j, i);
        ind = j + 1;
      }
      Mdr0 (mmax, 0, i) = (double) ind;
    }
  }
  return (mmax);
}

/* **************************************************************
 * Min function with two arguments...
 * ************************************************************** */

MDC *
    mdc_Min2 (MDC * m1, MDC * m2)
{
  int i, size;
  MDC *m = 0;

  /* Check sizes */

  if (MNR (m1) == MNR (m2) && MNC (m1) == MNC (m2))
  {
    m = mdc_Copy (m1);
    size = MNR (m1) * MNC (m1);
    for (i = 0; i < size; i++)
    {
      if (cabs (MdcV0 (m1, i)) < cabs (MdcV0 (m2, i)))
        MdcV0 (m, i) = MdcV0 (m2, i);
    }
  }
  else
    rerror ("min: matrix sizes must be the same");

  return (m);
}

/* **************************************************************
 * Matrix Sum functions.
 * ************************************************************** */

MDC *
mdc_Sum_BF (MDC * m, void *h)
{
  MDC *new=0;
  int i, j;
  MDR *cond=0;
  int lcr = 0, lcc = 0;

  if (h)
    cond = (MDR *) h;

  if ( EQNULL(m) )
  {
    new = mdc_CreateScalar (0,0);
    return (new);
  }

  if ( EQVECT(m) )
  {
    if (cond)
      lcc = SIZE(cond) - 1;

    // single row or single column matrix: add them all up
    Complex d = 0.0 + 0.0I;
    for (i=0; i<SIZE(m); i++)
    {
      if (isnand(MdcV0r(m, i)) || isnand(MdcV0i(m, i)))
        continue;

      if (cond)
        if (!mdrV0(cond,MIN(i,lcc)))
          continue;

      d = d + MdcV0 (m, i);
    }
    new = mdc_CreateScalar (RE(d), IM(d));
  }
  else
  {
    if (cond)
    {
      lcr = MNR(cond) - 1;
      lcc = MNC(cond) - 1;
    }

    new = mdc_Create (1, MNC(m));
    for (i=0; i<MNC(m); i++)
    {
      MdcV0(new,i) = 0.0 + 0.0I;
      for (j=0; j<MNR(m); j++)
      {
        if (isnand(Mdc0r(m,j,i)) || isnand(Mdc0i(m,j,i)))
          continue;

        if (cond)
          if (!mdr0(cond,MIN(j,lcr),MIN(i,lcc)))
            continue;

        MdcV0(new, i) += Mdc0 (m, j, i);
      }
    }
  }

  return (new);
}

MDC *
mdc_CumSum_BF (MDC * m)
{
  int i, j;
  MDC *new=0;

  if ( m->nrow == 0 || m->ncol == 0 )
  {
    new = mdc_CreateScalar (0, 0);
  }
  else if ( m->nrow == 1 || m->ncol == 1 )
  {
    new = mdc_Create (m->nrow, m->ncol);
    MdcV0r (new, 0) = MdcV0r (m, 0);
    MdcV0i (new, 0) = MdcV0i (m, 0);

    for (i = 1; i < (m->nrow * m->ncol); i++)
    {
      MdcV0r (new, i) = MdcV0r (new, i - 1) + MdcV0r (m, i);
      MdcV0i (new, i) = MdcV0i (new, i - 1) + MdcV0i (m, i);
    }
  }
  else
  {
    new = mdc_Create (m->nrow, m->ncol);

    /* Initialize first row. */
    for (i = 0; i < m->ncol; i++)
    {
      Mdc0r (new, 0, i) = Mdc0r (m, 0, i);
      Mdc0i (new, 0, i) = Mdc0i (m, 0, i);
    }

    /* Now compute running sum. */
    for (i = 0; i < m->ncol; i++)
    {
      for (j = 1; j < m->nrow; j++)
      {
        Mdc0r (new, j, i) = Mdc0r (new, j - 1, i) + Mdc0r (m, j, i);
        Mdc0i (new, j, i) = Mdc0i (new, j - 1, i) + Mdc0i (m, j, i);
      }
    }
  }
  return (new);
}

/* **************************************************************
 * Matrix Product functions.
 * ************************************************************** */

MDC *
mdc_Prod_BF (MDC * m)
{
  int i, j;
  MDC *new;

  if (MNR (m) * MNC (m) == 0)
  {
    new = mdc_CreateScalar (0, 0);
  }
  else if (MNR (m) == 1)
  {
    Complex c;
    c = Mdc0 (m, 0, 0);
    for (i = 1; i < MNC (m); i++)
    {
      c = c * MdcV0(m, i);
    }
    new = mdc_CreateScalar (RE(c), IM(c));
  }
  else
  {
    new = mdc_Create (1, MNC (m));
    for (i = 0; i < MNC (m); i++)
    {
      MdcV0 (new, i) = MdcV0 (m, i);
      for (j = 1; j < MNR (m); j++)
      {
        MdcV0(new, i) = MdcV0(new, i) * Mdc0(m, j, i);
      }
    }
  }
  return (new);
}

MDC *
mdc_CumProd_BF (MDC * m)
{
  int i, j;
  MDC *new;

  if (MNR (m) == 1 || MNC (m) == 1)
  {
    int N = MNR (m) * MNC (m);
    new = mdc_Create (MNR (m), MNC (m));
    MdcV0r (new, 0) = MdcV0r (m, 0);
    MdcV0i (new, 0) = MdcV0i (m, 0);

    for (i = 1; i < N; i++)
    {
      MdcV0(new, i) = MdcV0(new, i - 1) * MdcV0(m, i);
    }
  }
  else
  {
    new = mdc_Create (MNR (m), MNC (m));

    /* Initialize first row. */
    for (i = 0; i < MNC (m); i++)
    {
      Mdc0r (new, 0, i) = Mdc0r (m, 0, i);
      Mdc0i (new, 0, i) = Mdc0i (m, 0, i);
    }

    /* Now compute running product. */
    for (i = 0; i < MNC (m); i++)
    {
      for (j = 1; j < MNR (m); j++)
      {
        Mdc0(new, j, i) = Mdc0(new, j - 1, i) * Mdc0 (m, j, i);
      }
    }
  }
  return (new);
}

/* **************************************************************
 * Test for symmetry (Hermitian).
 * ************************************************************** */

int
mdc_IsSymmetric (MDC * m)
{
  int diag, i, j, nr, nc;

  nr = MNR (m);
  nc = MNC (m);
  diag = MIN(nr, nc);

  /*
   * First check the diagonal elements.
   * They must be real.
   */

  for (i = 0; i < diag; i++)
  {
    if (Mdc0i (m, i, i) != 0.0)
      return (0);
  }

  /*
   * Now test the rest (if we have gotten this far).
   */

  for (j = 0; j < nc; j++)
  {
    for (i = j + 1; i < nr; i++)
    {
      if (Mdc0r (m, i, j) != Mdc0r (m, j, i))
	return (0);
      if (Mdc0i (m, i, j) != -Mdc0i (m, j, i))
	return (0);
    }
  }
  return (1);
}

void
mdc_Detect_Inf (MDC * m)
{
  int n = MNR (m) * MNC (m);

  if (detect_inf_c (MDCPTR (m), n))
    rerror ("matrix contains Inf value");
}

void
mdc_Detect_Nan (MDC * m)
{
  int n = MNR (m) * MNC (m);

  if (detect_nan_c (MDCPTR (m), n))
    rerror ("matrix contains NaN value");
}

/* **************************************************************
 * Test a matrix for Infs...
 * ************************************************************** */

MDR *
mdc_IsInf (MDC * m)
{
  int i, size;
  MDR *mi;

  size = MNR (m) * MNC (m);
  mi = mdr_Create (MNR (m), MNC (m));

  for (i = 0; i < size; i++)
  {
    if (create_inf () == MdcV0r (m, i) || create_inf () == MdcV0i (m, i))
      MdrV0 (mi, i) = 1.0;
    else
      MdrV0 (mi, i) = 0.0;
  }
  return (mi);
}

/* **************************************************************
 * Test a matrix for Nans...
 * ************************************************************** */

MDR *
mdc_IsNan (MDC * m)
{
  int i, size;
  MDR *mi;

  size = MNR (m) * MNC (m);
  mi = mdr_Create (MNR (m), MNC (m));

  for (i = 0; i < size; i++)
  {
    if (MdcV0r (m, i) != MdcV0r (m, i) || MdcV0i (m, i) != MdcV0i (m, i))
      MdrV0 (mi, i) = 1.0;
    else
      MdrV0 (mi, i) = 0.0;
  }
  return (mi);
}

MDR *
mdc_Finite (MDC * m)
{
  int i, size;
  MDR *mi;

  size = MNR (m) * MNC (m);
  mi = mdr_Create (MNR (m), MNC (m));

  for (i = 0; i < size; i++)
  {
    if ((create_inf () != MdcV0r (m, i)) && (create_inf () != MdcV0i (m, i)) &&
	(-create_inf () != MdcV0r (m, i)) && (-create_inf () != MdcV0i (m, i))
	&& (MdcV0r (m, i) == MdcV0r (m, i)) && (MdcV0i (m, i) == MdcV0i (m, i)))
      MdrV0 (mi, i) = 1.0;
    else
      MdrV0 (mi, i) = 0.0;
  }
  return (mi);
}

/* **************************************************************
 * Sin function...
 * ************************************************************** */

MDC *
mdc_Sin (MDC * m)
{
  int i, size;
  MDC *mt;

  mt = mdc_Create (MNR (m), MNC (m));
  size = MNR (m) * MNC (m);
  for (i = 0; i < size; i++)
  {
    MdcV0 (mt, i) = csin(MdcV0 (m, i));
  }
  return (mt);
}

/* **************************************************************
 * Cos function...
 * ************************************************************** */

MDC *
mdc_Cos (MDC * m)
{
  int i, size;
  MDC *mt;

  mt = mdc_Create (MNR (m), MNC (m));
  size = MNR (m) * MNC (m);
  for (i = 0; i < size; i++)
  {
    MdcV0 (mt, i) = ccos (MdcV0 (m, i));
  }
  return (mt);
}

/* **************************************************************
 * Tan function...
 * ************************************************************** */

MDC *
mdc_Tan (MDC * m)
{
  int i, size;
  MDC *mt;

  mt = mdc_Create (MNR (m), MNC (m));
  size = MNR (m) * MNC (m);
  for (i = 0; i < size; i++)
  {
    MdcV0 (mt, i) = ctan (MdcV0 (m, i));
  }
  return (mt);
}

/* **************************************************************
 * Arc-Sin function...
 * ************************************************************** */

void *
mdc_ASin (MDC * m, int *type)
{
  int i, size;
  MDC *mt;

  mt = mdc_Create (MNR (m), MNC (m));
  size = MNR (m) * MNC (m);
  for (i = 0; i < size; i++)
  {
    MdcV0 (mt, i) = casin (MdcV0 (m, i));
  }
  *type = MATRIX_DENSE_COMPLEX;
  return (mt);
}

/* **************************************************************
 * Arc-Cos function...
 * ************************************************************** */

void *
mdc_ACos (MDC * m, int *type)
{
  int i, size;
  MDC *mt;

  mt = mdc_Create (MNR (m), MNC (m));
  size = MNR (m) * MNC (m);
  for (i = 0; i < size; i++)
  {
    MdcV0 (mt, i) = cacos (MdcV0 (m, i));
  }
  *type = MATRIX_DENSE_COMPLEX;
  return (mt);
}

/* **************************************************************
 * Tan function...
 * ************************************************************** */

void *
mdc_ATan (MDC * m, int *type)
{
  int i, size;
  MDC *mt;

  mt = mdc_Create (MNR (m), MNC (m));
  size = MNR (m) * MNC (m);
  for (i = 0; i < size; i++)
  {
    MdcV0 (mt, i) = catan (MdcV0 (m, i));
  }
  *type = MATRIX_DENSE_COMPLEX;
  return (mt);
}

/* **************************************************************
 * Sqrt function...
 * ************************************************************** */

void *
mdc_Sqrt (MDC * m, int *type)
{
  int i, size;
  MDC *mt;

  mt = mdc_Create (MNR (m), MNC (m));
  size = MNR (m) * MNC (m);
  for (i = 0; i < size; i++)
  {
    MdcV0 (mt, i) = csqrt (MdcV0 (m, i));
  }
  *type = MATRIX_DENSE_COMPLEX;
  return (mt);
}

/* **************************************************************
 * Log function...
 * ************************************************************** */

void *
mdc_Log (MDC * m, int *type)
{
  int i, size;
  MDC *mt;

  mt = mdc_Create (MNR (m), MNC (m));
  size = MNR (m) * MNC (m);
  for (i = 0; i < size; i++)
  {
    MdcV0 (mt, i) = clog (MdcV0 (m, i));
  }
  *type = MATRIX_DENSE_COMPLEX;
  return (mt);
}

/* **************************************************************
 * Log10 function...
 * ************************************************************** */

#define log10e 0.43429448190325182765

void *
mdc_Log10 (MDC * m, int *type)
{
  int i, size;
  MDC *mt;

  mt = mdc_Create (MNR (m), MNC (m));
  size = MNR (m) * MNC (m);
  for (i = 0; i < size; i++)
    MdcV0 (mt, i) = log10e * clog (MdcV0 (m, i));
  *type = MATRIX_DENSE_COMPLEX;
  return (mt);
}

/* **************************************************************
 * Exp function...
 * ************************************************************** */

MDC *
mdc_Exp (MDC * m)
{
  int i, size;
  MDC *mt;

  mt = mdc_Create (MNR (m), MNC (m));
  size = MNR (m) * MNC (m);
  for (i = 0; i < size; i++)
  {
    MdcV0 (mt, i) = cexp (MdcV0 (m, i));
  }
  return (mt);
}

/* **************************************************************
 * diag() function...
 * ************************************************************** */

MDC *
mdc_Diag (MDC * marg, int k)
{
  int i, smin, size;
  MDC *m;

  if (MNR (marg) == 1 || MNC (marg) == 1)
  {
    /* Create a diagonal matrix */
    size = MAX(MNR (marg), MNC (marg)) + ABS(k);
    m = mdc_Create (size, size);
    mdc_Zero (m);
    if (k < 0)
    {
      for (i = 1 - k; i <= size; i++)
      {
        Mdc1r (m, i, (i + k)) = MdcV1r (marg, (i + k));
        Mdc1i (m, i, (i + k)) = MdcV1i (marg, (i + k));
      }
    }
    else
    {
      for (i = 1; i <= size - k; i++)
      {
        Mdc1r (m, i, (i + k)) = MdcV1r (marg, i);
        Mdc1i (m, i, (i + k)) = MdcV1i (marg, i);
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

    m = mdc_Create (size, 1);
    mdc_Zero (m);
    if (k >= 0)
    {
      for (i = 1; i <= size; i++)
      {
        Mdc1r (m, i, 1) = Mdc1r (marg, i, i + k);
        Mdc1i (m, i, 1) = Mdc1i (marg, i, i + k);
      }
    }
    else
    {
      for (i = 1; i <= size; i++)
      {
        Mdc1r (m, i, 1) = Mdc1r (marg, i - k, i);
        Mdc1i (m, i, 1) = Mdc1i (marg, i - k, i);
      }
    }
  }
  return (m);
}

/* **************************************************************
 * Return the real part of a matrix.
 * ************************************************************** */

MDR *
mdc_Real_BF (MDC * m)
{
  int i, size;
  MDR *mr;

  mr = mdr_Create (MNR (m), MNC (m));
  size = MNR (m) * MNC (m);
  for (i = 0; i < size; i++)
  {
    MdrV0 (mr, i) = MdcV0r (m, i);
  }
  return (mr);
}

/* **************************************************************
 * Return the imaginary part of a matrix.
 * ************************************************************** */

MDR *
mdc_Imag_BF (MDC * m)
{
  int i, size;
  MDR *mr;

  mr = mdr_Create (MNR (m), MNC (m));
  size = MNR (m) * MNC (m);
  for (i = 0; i < size; i++)
  {
    MdrV0 (mr, i) = MdcV0i (m, i);
  }
  return (mr);
}

/* **************************************************************
 * Return the complex-conjugate of a matrix.
 * ************************************************************** */

MDC *
mdc_Conj_BF (MDC * m)
{
  int i, size;
  MDC *mc;

  mc = mdc_Create (MNR (m), MNC (m));
  size = MNR (m) * MNC (m);
  for (i = 0; i < size; i++)
  {
    MdcV0r (mc, i) = MdcV0r (m, i);
    MdcV0i (mc, i) = -MdcV0i (m, i);
  }
  return (mc);
}

/* **************************************************************
 * Find the indices of non-zero elements in a matrix.
 * ************************************************************** */

MDR *
mdc_Find_BF (MDC * m)
{
  int i, j, size;
  double *dtmp;
  MDR *mf;

  /* Make mfind as big as it could be, reduce later */
  mf = mdr_Create (1, MNR (m) * MNC (m));

  j = 0;
  size = MNR (m) * MNC (m);
  for (i = 0; i < size; i++)
  {
    if ((MdcV0r (m, i) != 0.0) || (MdcV0i (m, i) != 0.0))
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
    MDPTR(mf) = dtmp;
    mf->ncol = j;
  }

  return (mf);
}

/* **************************************************************
 * Sort a COMPLEX matrix. Run it through ABS first.
 * ************************************************************** */

extern MDR *mdr_CreateFillSind (int nrow, int ncol);

Btree *
mdc_Sort_BF (MDC * m)
{
  int i, j, n, size;
  Btree *bt;
  Ent *eind, *eval;
  MDR *sind, *mtmp;
  MDC *mcopy;

  if (MNR (m) * MNC (m) == 0)
  {
    sind  = mdr_Create (0,0);
    mcopy = mdc_Create (0,0);
  }
  else if (MNR (m) == 1 || MNC (m) == 1)
  {
    /* Vector sort */
    size = MNR (m) * MNC (m);
    n = MAX(MNR (m), MNC (m));
    sind = mdr_CreateFillSind (n, 1);
    sind->ncol = sind->nrow;
    sind->nrow = 1;
    mtmp = mdc_Abs (m);
    r_sort ((double *) MDRPTR (mtmp), 0, n - 1, (double *) MDRPTR (sind));

    /* Now sort [m] according to [sind] */
    mdr_Destroy (mtmp);
    mcopy = mdc_Create (MNR (m), MNC (m));
    for (i = 1; i <= size; i++)
    {
      MdcV1r (mcopy, i) = MdcV1r (m, ((int) MdrV1 (sind, i)));
      MdcV1i (mcopy, i) = MdcV1i (m, ((int) MdrV1 (sind, i)));
    }
  }
  else
  {
    /* Matrix sort (column-wise) */
    n = MNR (m);
    sind = mdr_CreateFillSind (MNR (m), MNC (m));
    mtmp = mdc_Abs (m);
    for (i = 0; i < MNC (m); i++)
    {
      r_sort (&MdrV0(mtmp,i*n), 0, n - 1, &MdrV0(sind,i*n));
    }

    /* Now sort [m] according to [sind] */
    mdr_Destroy (mtmp);
    mcopy = mdc_Create (MNR (m), MNC (m));
    for (i = 1; i <= MNC (m); i++)
    {
      for (j = 1; j <= MNR (m); j++)
      {
        Mdc1r (mcopy, j, i) = Mdc1r (m, ((int) Mdr1 (sind, j, i)), i);
        Mdc1i (mcopy, j, i) = Mdc1i (m, ((int) Mdr1 (sind, j, i)), i);
      }
    }
  }

  bt = btree_Create ();

  eind = ent_Create ();
  ent_SetType (eind, MATRIX_DENSE_REAL);
  ent_data (eind) = sind;
  install (bt, "ind", eind);

  eval = ent_Create ();
  ent_SetType (eval, MATRIX_DENSE_COMPLEX);
  ent_data (eval) = mcopy;
  install (bt, "val", eval);

  return (bt);
}

/* **************************************************************
 * Sign function.
 * ************************************************************** */

MDC *
mdc_Sign_BF (MDC * m)
{
  int i, size;
  MDC *sm;

  size = MNR (m) * MNC (m);

  sm = mdc_Create (MNR (m), MNC (m));
  for (i = 0; i < size; i++)
  {
    MdcV0 (sm, i) = MdcV0 (m, i) / cabs (MdcV0 (m, i));
  }
  return (sm);
}

/* **************************************************************
 * Increment (++) a matrix. Do the operation in place.
 * ************************************************************** */

MDC *
mdc_Increment (MDC * m)
{
  int i, size;

  size = MNR (m) * MNC (m);
  for (i = 0; i < size; i++)
  {
    MdcV0r (m, i) = MdcV0r (m, i) + 1.0;
    MdcV0i (m, i) = MdcV0i (m, i) + 1.0;
  }

  return (m);
}

MDC *
mdc_Decrement (MDC * m)
{
  int i, size;

  size = MNR (m) * MNC (m);
  for (i = 0; i < size; i++)
  {
    MdcV0r (m, i) = MdcV0r (m, i) - 1.0;
    MdcV0i (m, i) = MdcV0i (m, i) - 1.0;
  }

  return (m);
}

size_t
mdc_Sizeof (MDC * m)
{
  int size = MNR (m) * MNC (m);
  return (size_t) (sizeof (Complex) * size);
}

MDR *
mdc_Any (MDC * m)
{
  int i, j;
  MDR *new;

  if (MNR (m) == 1)		/* Vector operation */
  {
    new = mdr_Create (1, 1);
    Mdr1 (new, 1, 1) = 0.0;
    for (i = 1; i <= MNC (m); i++)
    {
      if ((MdcV1r (m, i) != 0.0) || (MdcV1i (m, 1) != 0.0))
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
	if ((Mdc1r (m, j, i) != 0.0) || (Mdc1i (m, j, i) != 0.0))
	{
	  Mdr1 (new, 1, i) = 1.0;
	  break;
	}
      }
    }

  }

  return (new);
}

MDR *
mdc_All (MDC * m)
{
  int i, j;
  MDR *new;

  if (EQVECT(m))   /* Vector operation */
  {
    new = mdr_Create (1, 1);
    Mdr1 (new, 1, 1) = 1.0;
    for (i = 1; i <= SIZE(m); i++)
    {
      if ((Mdc1r (m, 1, i) == 0.0) && (Mdc1i (m, 1, i) ==0))
      {
        Mdr1 (new, 1, 1) = 0.0;
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
      Mdr1 (new, 1, i) = 1.0;
      for (j = 1; j <= MNR (m); j++)
      {
        if ((Mdc1r (m, j, i) == 0.0) && (Mdc1i (m, j, i) == 0.0))
        {
          Mdr1 (new, 1, i) = 0.0;
          break;
        }
      }
    }
  }

  return (new);
}

MDC *
mdc_Mod (MDC * m1, MDC * m2)
{
  int i, size;
  MDC *new;

  if (MNR (m1) != MNR (m2) || MNC (m1) != MNC (m2))
    rerror ("mod: matrices must have the same dimensions");

  new = mdc_Create (MNR (m1), MNC (m1));
  size = MNR (m1) * MNC (m1);
  for (i = 0; i < size; i++)
  {
    MdcV0 (new, i) = complex_Mod (MdcV0 (m1, 0), MdcV0 (m2, i));
  }
  return (new);
}
