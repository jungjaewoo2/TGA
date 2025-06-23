/* mdr_mdc.c Matrix Dense Real and Matrix Dense Complex interaction. */

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
#include "mathl.h"
#include "mem.h"
#include "mdr.h"
#include "mdrf1.h"
#include "mdrf2.h"
#include "mdc.h"
#include "mdcf1.h"
#include "mdcf2.h"
#include "util.h"

#include <stdio.h>
#include <math.h>


/* **************************************************************
 * Coerce/copy a DENSE matrix into a DENSE-COMPLEX matrix.
 * ************************************************************** */
MDC * mdr_coerce_mdc (MD * m)
{
  MDC *mc;
  int i, size=SIZE(m);
  if (size>0)
  {
    mc = mdc_Create (MNR (m), MNC (m));
    for (i = 0; i < size; i++)
    {
      MdcV0 (mc, i) = mdcV0 (m, i);
    }
  }
  else
    mc = mdc_Create (0, 0);

  return (mc);
}


//
// Add two matrices.
//
MDC * md_mdc_AddTo (MD * m1, MD * m2)
{
  if (!MD_TYPE_COMPLEX(m1) && !MD_TYPE_COMPLEX(m2))
    rerror("Horrible Internal Error: We should not be here!");

  MD *new=0;

  // Quazi-emulate Matlab empty matrix behavior if matrices are empty
  if ( !SIZE(m1) )
  {
    new = mdr_coerce_mdc (m2);
    md_Destroy (m1);
    return (new);
  }
  if ( !SIZE(m2) )
  {
    new = mdr_coerce_mdc (m1);
    md_Destroy (m1);
    return (new);
  }

  int i, j, i1, i2, j1, j2;

  i1 = MNR(m1);
  j1 = MNC(m1);
  i2 = MNR(m2);
  j2 = MNC(m2);

  // m2 is bigger in size then m1:
  //  the result will not fit into m1, so we have to create a new matrix
  // or m1 is not complex matrix, but we are here because one of them is complex
  if ( (MNR(m1)<MNR(m2)) || (MNC(m1)<MNC(m2)) || (!MD_TYPE_COMPLEX(m1)) )
  {
    int nr = MAX(i1, i2);
    int nc = MAX(j1, j2);
    i1--;
    i2--;
    j1--;
    j2--;
    new = mdc_Create(nr, nc);
    for(i=0;i<nr;i++)
      for(j=0;j<nc;j++)
        Mdc0(new,i,j) = mdc0(m1, MIN(i,i1), MIN(j,j1)) + mdc0(m2, MIN(i,i2), MIN(j,j2));
    md_Destroy (m1);
  }
  else
  {
    // m2 is smaller then or equal in size to m1
    // and m1 is complex matrix
    //  the result will fit into m1, so we don't have to create a new matrix
    new = m1;
    int nr = MAX(i1, i2);
    int nc = MAX(j1, j2);
    i1--;
    i2--;
    j1--;
    j2--;
    for(i=0;i<nr;i++)
      for(j=0;j<nc;j++)
        Mdc0(new,i,j) += mdc0(m2, MIN(i,i2), MIN(j,j2));
  }

  return (new);
}

//
// Add two matrices.
//
MDC * md_mdc_SubtractFrom (MD * m1, MD * m2)
{
  if (!MD_TYPE_COMPLEX(m1) && !MD_TYPE_COMPLEX(m2))
    rerror("Horrible Internal Error: We should not be here!");

  MD *new=0;

  // Quazi-emulate Matlab empty matrix behavior if matrices are empty
  if ( !SIZE(m1) )
  {
    new = mdr_coerce_mdc (m2);
    md_Destroy (m1);
    return (new);
  }
  if ( !SIZE(m2) )
  {
    new = mdr_coerce_mdc (m1);
    md_Destroy (m1);
    return (new);
  }

  int i, j, i1, i2, j1, j2;

  i1 = MNR(m1);
  j1 = MNC(m1);
  i2 = MNR(m2);
  j2 = MNC(m2);

  // m2 is bigger in size then m1:
  //  the result will not fit into m1, so we have to create a new matrix
  // or m1 is not complex matrix, but we are here because one of them is complex
  if ( (MNR(m1)<MNR(m2)) || (MNC(m1)<MNC(m2)) || (!MD_TYPE_COMPLEX(m1)) )
  {
    int nr = MAX(i1, i2);
    int nc = MAX(j1, j2);
    i1--;
    i2--;
    j1--;
    j2--;
    new = mdc_Create(nr, nc);
    for(i=0;i<nr;i++)
      for(j=0;j<nc;j++)
        Mdc0(new,i,j) = mdc0(m1, MIN(i,i1), MIN(j,j1)) - mdc0(m2, MIN(i,i2), MIN(j,j2));
    md_Destroy (m1);
  }
  else
  {
    // m2 is smaller then or equal in size to m1
    // and m1 is complex matrix
    //  the result will fit into m1, so we don't have to create a new matrix
    new = m1;
    int nr = MAX(i1, i2);
    int nc = MAX(j1, j2);
    i1--;
    i2--;
    j1--;
    j2--;
    for(i=0;i<nr;i++)
      for(j=0;j<nc;j++)
        Mdc0(new,i,j) -= mdc0(m2, MIN(i,i2), MIN(j,j2));
  }

  return (new);
}

MDC * md_mdc_ElMultiplyBy (MD * m1, MD * m2)
{
  if (!MD_TYPE_COMPLEX(m1) && !MD_TYPE_COMPLEX(m2))
    rerror("Horrible Internal Error: We should not be here!");

  MD *new=0;

  // Quazi-emulate Matlab empty matrix behavior if matrices are empty
  if ( !SIZE(m1) )
  {
    new = mdr_coerce_mdc (m2);
    md_Destroy (m1);
    return (new);
  }
  if ( !SIZE(m2) )
  {
    new = mdr_coerce_mdc (m1);
    md_Destroy (m1);
    return (new);
  }

  int i, j, i1, i2, j1, j2;

  i1 = MNR(m1);
  j1 = MNC(m1);
  i2 = MNR(m2);
  j2 = MNC(m2);

  // m2 is bigger in size then m1:
  //  the result will not fit into m1, so we have to create a new matrix
  // or m1 is not complex matrix, but we are here because one of them is complex
  if ( (MNR(m1)<MNR(m2)) || (MNC(m1)<MNC(m2)) || (!MD_TYPE_COMPLEX(m1)) )
  {
    int nr = MAX(i1, i2);
    int nc = MAX(j1, j2);
    i1--;
    i2--;
    j1--;
    j2--;
    new = mdc_Create(nr, nc);
    for(i=0;i<nr;i++)
      for(j=0;j<nc;j++)
        Mdc0(new,i,j) = mdc0(m1, MIN(i,i1), MIN(j,j1)) * mdc0(m2, MIN(i,i2), MIN(j,j2));
    md_Destroy (m1);
  }
  else
  {
    // m2 is smaller then or equal in size to m1
    // and m1 is complex matrix
    //  the result will fit into m1, so we don't have to create a new matrix
    new = m1;
    int nr = MAX(i1, i2);
    int nc = MAX(j1, j2);
    i1--;
    i2--;
    j1--;
    j2--;
    for(i=0;i<nr;i++)
      for(j=0;j<nc;j++)
        Mdc0(new,i,j) *= mdc0(m2, MIN(i,i2), MIN(j,j2));
  }

  return (new);
}

MDC * md_mdc_ElRdivideBy (MD * m1, MD * m2)
{
  if (!MD_TYPE_COMPLEX(m1) && !MD_TYPE_COMPLEX(m2))
    rerror("Horrible Internal Error: We should not be here!");

  MD *new=0;

  // Quazi-emulate Matlab empty matrix behavior if matrices are empty
  if ( !SIZE(m1) )
  {
    new = mdr_coerce_mdc (m2);
    md_Destroy (m1);
    return (new);
  }
  if ( !SIZE(m2) )
  {
    new = mdr_coerce_mdc (m1);
    md_Destroy (m1);
    return (new);
  }

  int i, j, i1, i2, j1, j2;

  i1 = MNR(m1);
  j1 = MNC(m1);
  i2 = MNR(m2);
  j2 = MNC(m2);

  // m2 is bigger in size then m1:
  //  the result will not fit into m1, so we have to create a new matrix
  // or m1 is not complex matrix, but we are here because one of them is complex
  if ( (MNR(m1)<MNR(m2)) || (MNC(m1)<MNC(m2)) || (!MD_TYPE_COMPLEX(m1)) )
  {
    int nr = MAX(i1, i2);
    int nc = MAX(j1, j2);
    i1--;
    i2--;
    j1--;
    j2--;
    new = mdc_Create(nr, nc);
    for(i=0;i<nr;i++)
      for(j=0;j<nc;j++)
        Mdc0(new,i,j) = mdc0(m1, MIN(i,i1), MIN(j,j1)) / mdc0(m2, MIN(i,i2), MIN(j,j2));
    md_Destroy (m1);
  }
  else
  {
    // m2 is smaller then or equal in size to m1
    // and m1 is complex matrix
    //  the result will fit into m1, so we don't have to create a new matrix
    new = m1;
    int nr = MAX(i1, i2);
    int nc = MAX(j1, j2);
    i1--;
    i2--;
    j1--;
    j2--;
    for(i=0;i<nr;i++)
      for(j=0;j<nc;j++)
        Mdc0(new,i,j) /= mdc0(m2, MIN(i,i2), MIN(j,j2));
  }

  return (new);
}

/* **************************************************************
 * Coerce a DENSE-COMPLEX matrix into a DENSE-REAL matrix.
 * ************************************************************** */

MDR *
mdc_coerce_mdr (MDC * m)
{
  MDR *mr;
  int i, size=SIZE(m);

  if (size>0)
  {
    mr = mdr_Create (MNR (m), MNC (m));
    for (i = 0; i < size; i++)
      MdrV0 (mr, i) = MdcV0r (m, i);
  }
  else
    mr = mdr_Create(0,0);

  return (mr);
}

/* **************************************************************
 * Functions to append dense-real and dense-complex together.
 * ************************************************************** */

MDC *
mdr_mdc_Append (MDR * mr, MDC * mc)
{
  int i, j, nrow, ncol;
  MDC *new=0;

  /* Check for empty matrices... */
  if (SIZE(mr)<1 && SIZE(mc)<1)
  {
    new = mdc_Create(0,0);
    return new;
  }
  else if (SIZE(mr)<1)
  {
    new = mdc_Copy (mc);
    return (new);
  }
  else if (SIZE(mc)<1)
  {
    new = mdr_coerce_mdc (mr);;
    return (new);
  }

  /* Do the stack. */
  /* Create the new matrix, large enough to hold both. */
  int r1 = MNR(mr);
  int c1 = MNC(mr);
  int r2 = MNR(mc);
  int c2 = MNC(mc);
  nrow = MAX(r1,r2);
  ncol = c1 + c2;

  new = mdc_Create (nrow, ncol);
  for (i=0; i<nrow; i++)
  {
    for (j=0; j<c1; j++)
    {
      Mdc0r(new,i,j) = mdr0(mr,MIN(i,r1-1),j);
      Mdc0i(new,i,j) = 0;
    }
    for (j=0; j<c2; j++)
      Mdc0(new,i,j+c1)  = Mdc0(mc,MIN(i,r2-1),j);
  }

  return (new);
}

MDC *
mdc_mdr_Append (MDC * mc, MDR * mr)
{
  int i, j, nrow, ncol;
  MDC *new=0;

  /* Check for empty matrices... */
  if (SIZE(mr)<1 && SIZE(mc)<1)
  {
    new = mdc_Create(0,0);
    return new;
  }
  else if (SIZE(mr)<1)
  {
    new = mdc_Copy (mc);
    return (new);
  }
  else if (SIZE(mc)<1)
  {
    new = mdr_coerce_mdc (mr);
    return (new);
  }

  /* Do the stack. */
  /* Create the new matrix, large enough to hold both. */
  int r1 = MNR(mc);
  int c1 = MNC(mc);
  int r2 = MNR(mr);
  int c2 = MNC(mr);
  nrow = MAX(r1,r2);
  ncol = c1 + c2;

  new = mdc_Create (nrow, ncol);
  for (i=0; i<nrow; i++)
  {
    for (j=0; j<c1; j++)
      Mdc0(new,i,j) = Mdc0(mc,MIN(i,r1-1),j);
    for (j=0; j<c2; j++)
    {
      Mdc0r(new,i,j+c1) = mdr0(mr,MIN(i,r2-1),j);
      Mdc0i(new,i,j+c1) = 0;
    }
  }

  return (new);
}

/* **************************************************************
 * Functions to stack dense-real and dense-complex together.
 * ************************************************************** */

MDC *
mdr_mdc_Stack (MDR * mr, MDC * mc)
{
  MDC *mrc, *new;

  mrc = mdr_coerce_mdc (mr);
  new = mdc_Stack (mrc, mc);
  mdc_Destroy (mrc);

  return (new);
}

MDC *
mdc_mdr_Stack (MDC * mc, MDR * mr)
{
  MDC *mrc, *new;

  mrc = mdr_coerce_mdc (mr);
  new = mdc_Stack (mc, mrc);
  mdc_Destroy (mrc);

  return (new);
}

/* **************************************************************
 * Functions to add dense-real and dense-complex together.
 * ************************************************************** */

MDC * mdr_mdc_Add (MDR * mr, MDC * mc)
{
  MDC *mrc, *new;

  mrc = mdr_coerce_mdc (mr);
  new = mdc_Add (mrc, mc);
  mdc_Destroy (mrc);

  return (new);
}

MDC *
mdc_mdr_Add (MDC * mc, MDR * mr)
{
  MDC *mrc, *new;

  mrc = mdr_coerce_mdc (mr);
  new = mdc_Add (mc, mrc);
  mdc_Destroy (mrc);

  return (new);
}

/* **************************************************************
 * Functions to subtract dense-real and dense-complex matrices.
 * ************************************************************** */

MDC *
mdr_mdc_Subtract (MDR * mr, MDC * mc)
{
  MDC *mrc, *new;

  mrc = mdr_coerce_mdc (mr);
  new = mdc_Subtract (mrc, mc);
  mdc_Destroy (mrc);

  return (new);
}

MDC *
mdc_mdr_Subtract (MDC * mc, MDR * mr)
{
  MDC *mrc, *new;

  mrc = mdr_coerce_mdc (mr);
  new = mdc_Subtract (mc, mrc);
  mdc_Destroy (mrc);

  return (new);
}

/* **************************************************************
 * Functions to Multiply dense-real and dense-complex matrices.
 * ************************************************************** */

MDC *
mdr_mdc_Multiply (MDR * mr, MDC * mc, int *rtype)
{
  MDC *mrc, *new;

  mrc = mdr_coerce_mdc (mr);
  new = mdc_Multiply (mrc, mc, rtype);
  mdc_Destroy (mrc);

  return (new);
}

MDC *
mdc_mdr_Multiply (MDC * mc, MDR * mr, int *rtype)
{
  MDC *mrc, *new;

  mrc = mdr_coerce_mdc (mr);
  new = mdc_Multiply (mc, mrc, rtype);
  mdc_Destroy (mrc);

  return (new);
}

/* **************************************************************
 * Functions to Element-Multiply dense-real and dense-complex matrices.
 * ************************************************************** */

MDC *
mdr_mdc_ElMultiply (MDR * mr, MDC * mc)
{
  MDC *mrc, *new;

  mrc = mdr_coerce_mdc (mr);
  new = mdc_ElMultiply (mrc, mc);
  mdc_Destroy (mrc);

  return (new);
}

MDC *
mdc_mdr_ElMultiply (MDC * mc, MDR * mr)
{
  MDC *mrc, *new;

  mrc = mdr_coerce_mdc (mr);
  new = mdc_ElMultiply (mc, mrc);
  mdc_Destroy (mrc);

  return (new);
}

/* **************************************************************
 * Functions to Right-Divide dense-real and dense-complex matrices.
 * ************************************************************** */

MDC *
mdr_mdc_Rdivide (MDR * mr, MDC * mc)
{
  MDC *mrc, *new;

  mrc = mdr_coerce_mdc (mr);
  new = mdc_Rdivide (mrc, mc);
  mdc_Destroy (mrc);

  return (new);
}

MDC *
mdc_mdr_Rdivide (MDC * mc, MDR * mr)
{
  MDC *mrc, *new;

  mrc = mdr_coerce_mdc (mr);
  new = mdc_Rdivide (mc, mrc);
  mdc_Destroy (mrc);

  return (new);
}

/* **************************************************************
 * Functions to Element-Right-Divide dense-real and
 * dense-complex matrices.
 * ************************************************************** */

MDC *
mdr_mdc_ElRdivide (MDR * mr, MDC * mc, int *rtype)
{
  MDC *mrc, *new;

  mrc = mdr_coerce_mdc (mr);
  new = mdc_ElRdivide (mrc, mc, rtype);
  mdc_Destroy (mrc);

  return (new);
}

MDC *
mdc_mdr_ElRdivide (MDC * mc, MDR * mr, int *rtype)
{
  MDC *mrc, *new;

  mrc = mdr_coerce_mdc (mr);
  new = mdc_ElRdivide (mc, mrc, rtype);
  mdc_Destroy (mrc);

  return (new);
}

/* **************************************************************
 * Functions to Ledt-Divide dense-real and dense-complex matrices.
 * ************************************************************** */

MDC *
mdr_mdc_Ldivide (MDR * mr, MDC * mc)
{
  MDC *mrc, *new;

  mrc = mdr_coerce_mdc (mr);
  new = mdc_Ldivide (mrc, mc);
  mdc_Destroy (mrc);

  return (new);
}

MDC *
mdc_mdr_Ldivide (MDC * mc, MDR * mr)
{
  MDC *mrc, *new;

  mrc = mdr_coerce_mdc (mr);
  new = mdc_Ldivide (mc, mrc);
  mdc_Destroy (mrc);

  return (new);
}

/* **************************************************************
 * Functions to Element-Left-Divide dense-real and
 * dense-complex matrices.
 * ************************************************************** */

MDC *
mdr_mdc_ElLdivide (MDR * mr, MDC * mc)
{
  MDC *mrc, *new;

  mrc = mdr_coerce_mdc (mr);
  new = mdc_ElLdivide (mrc, mc);
  mdc_Destroy (mrc);

  return (new);
}

MDC *
mdc_mdr_ElLdivide (MDC * mc, MDR * mr)
{
  MDC *mrc, *new;

  mrc = mdr_coerce_mdc (mr);
  new = mdc_ElLdivide (mc, mrc);
  mdc_Destroy (mrc);

  return (new);
}

void *
mdr_mdc_ElPower (MDR * mr, MDC * mc, int *type)
{
  MDC *mrc, *new;

  mrc = mdr_coerce_mdc (mr);
  new = mdc_ElPower (mrc, mc, type);
  mdc_Destroy (mrc);

  return (new);
}

void *
mdc_mdr_ElPower (MDC * mc, MDR * mr, int *type)
{
  MDC *mrc, *new;

  mrc = mdr_coerce_mdc (mr);
  new = mdc_ElPower (mc, mrc, type);
  mdc_Destroy (mrc);

  return (new);
}

void *
mdr_mdc_Power (MDR * mr, MDC * mc, int *type)
{
  MDC *mrc, *new;

  mrc = mdr_coerce_mdc (mr);
  new = mdc_Power (mrc, mc, type);
  mdc_Destroy (mrc);

  return (new);
}

void *
mdc_mdr_Power (MDC * mc, MDR * mr, int *type)
{
  MDC *mrc, *new;

  mrc = mdr_coerce_mdc (mr);
  new = mdc_Power (mc, mrc, type);
  mdc_Destroy (mrc);

  return (new);
}

/* **************************************************************
 * Coerce matrix indices into ints for use later...
 * Return an int array, the 1st element is the number of elements.
 * ************************************************************** */

int *
mdc_IndCoerceInt (MDC * m, MDR * i)
{
  int *ind, k, size;

  if ((size = MNR (i) * MNC (i)) == 0)
  {
    return (0);
  }

  ind = (int *) GC_MAIOP ((size + 1) * sizeof (int));
  ind[0] = size;

  for (k = 1; k <= size; k++)
  {
    ind[k] = (int) (MdrV1 (i, k));
    if (ind[k] <= 0)
    {
      fprintf (stderr, "matrix index <= 0 not allowed\n");
      fprintf (stderr, "index value: %i\n", ind[k]);
      GC_FREE (ind);
      rerror ("matrix index coerce");
    }
  }

  return (ind);
}

/* **************************************************************
 * Handle:
 * MDR[i;j] = [MDC];
 * ************************************************************** */

MDC *
mdr_mdc_MatrixAssign (MDR * var, int *i, int *j, MDC * rhs)
{
  MDC *new;

  /* Coerce LHS into MATRIX_DENSE_COMPLEX. */

  new = mdr_coerce_mdc (var);

  /* Now do the assignment. */

  new = mdc_MatrixAssign (new, i, j, rhs);
  mdr_Destroy (var);

  return (new);
}

MDC *
mdr_mdc_MatrixAssignR (MDR * var, int *i, MDC * rhs)
{
  MDC *new;

  /* Coerce LHS into MATRIX_DENSE_COMPLEX. */

  new = mdr_coerce_mdc (var);

  /* Now do the assignment. */

  new = mdc_MatrixAssignR (new, i, rhs);
  mdr_Destroy (var);

  return (new);
}

MDC *
mdr_mdc_MatrixAssignC (MDR * var, int *j, MDC * rhs)
{
  MDC *new;

  /* Coerce LHS into MATRIX_DENSE_COMPLEX. */

  new = mdr_coerce_mdc (var);

  /* Now do the assignment. */

  new = mdc_MatrixAssignC (new, j, rhs);
  mdr_Destroy (var);

  return (new);
}

/* **************************************************************
 * Handle:
 * MDC[i;j] = [MDR];
 * ************************************************************** */

MDC *
mdc_mdr_MatrixAssign (MDC * var, int *i, int *j, MDR * rhs)
{
  MDC *nrhs;

  /* Coerce RHS into MATRIX_DENSE_COMPLEX. */

  nrhs = mdr_coerce_mdc (rhs);

  /* Now do the assignment. */

  var = mdc_MatrixAssign (var, i, j, nrhs);
  mdc_Destroy (nrhs);

  return (var);
}

MDC *
mdc_mdr_MatrixAssignR (MDC * var, int *i, MDR * rhs)
{
  MDC *nrhs;

  /* Coerce RHS into MATRIX_DENSE_COMPLEX. */

  nrhs = mdr_coerce_mdc (rhs);

  /* Now do the assignment. */

  var = mdc_MatrixAssignR (var, i, nrhs);
  mdc_Destroy (nrhs);

  return (var);
}

MDC *
mdc_mdr_MatrixAssignC (MDC * var, int *j, MDR * rhs)
{
  MDC *nrhs;

  /* Coerce RHS into MATRIX_DENSE_COMPLEX. */

  nrhs = mdr_coerce_mdc (rhs);

  /* Now do the assignment. */

  var = mdc_MatrixAssignC (var, j, nrhs);
  mdc_Destroy (nrhs);

  return (var);
}

/* **************************************************************
 * Handle COMPLEX-REAL cross-assignments.
 * ************************************************************** */

MDC *
mdr_mdc_VectorAssign (MDR * mlhs, int *i, MDC * mrhs)
{
  MDC *tmp;
  tmp = mdr_coerce_mdc (mlhs);
  tmp = mdc_VectorAssign (tmp, i, mrhs);
  mdr_Destroy (mlhs);
  return (tmp);
}

MDC *
mdc_mdr_VectorAssign (MDC * mlhs, int *i, MDR * mrhs)
{
  MDC *tmp;
  tmp = mdr_coerce_mdc (mrhs);
  mlhs = mdc_VectorAssign (mlhs, i, tmp);
  mdc_Destroy (tmp);
  return (mlhs);
}

/* **************************************************************
 * Compare two matrices (m1 == m2)
 * ************************************************************** */

MDR *
mdc_Eq (MDC * m1, MDC * m2)
{
  int i, size;
  MDR *new;

  /* Special cases... */

  if (EQNULL(m1) && EQNULL(m2))
  {
    new = mdr_Create (0, 0);
    return (new);
  }
  else if (EQNULL(m1))
  {
    new = mdr_Create (MNR(m2),MNC(m2));
    mdr_Zero(new);
    return (new);
  }
  else if (EQNULL(m2))
  {
    new = mdr_Create (MNR(m1),MNC(m1));
    mdr_Zero(new);
    return (new);
  }
  else if ( EQSCAL(m1) )
  {
    size = MNR (m2) * MNC (m2);
    new  = mdr_Create (MNR (m2), MNC (m2));
    for (i=0; i<size; i++)
    {
      if ( EQCMPL(MdcV0(m1, 0),MdcV0(m2, i)) )
        MdrV0 (new, i) = 1.0;
      else
        MdrV0 (new, i) = 0.0;
    }
    return (new);
  }
  else if ( EQSCAL(m2) )
  {
    size = MNR (m1) * MNC (m1);
    new = mdr_Create (MNR (m1), MNC (m1));
    for (i=0; i<size; i++)
    {
      if ( EQCMPL(MdcV0(m1, i), MdcV0(m2, 0)) )
        MdrV0 (new, i) = 1.0;
      else
        MdrV0 (new, i) = 0.0;
    }
    return (new);
  }
  else if ( !EQSIZE(m1,m2) )
  {
    fprintf (stderr, "matrices must have the same dimension for comparision\n");
    fprintf (stderr, "m1: nrow = %i ncol = %i m2: nrow = %i, ncol = %i\n",
	     MNR (m1), MNC (m1), MNR (m2), MNC (m2));
    rerror ("matrix-compare");
  }

  size = MNR (m1) * MNC (m1);
  new = mdr_Create (MNR (m1), MNC (m1));

  for (i = 0; i < size; i++)
  {
    if ( EQCMPL(MdcV0(m1,i), MdcV0(m2,i)) )
      MdrV0 (new, i) = 1.0;
    else
      MdrV0 (new, i) = 0.0;
  }

  return (new);
}

MDR *
mdc_Ne (MDC * m1, MDC * m2)
{
  int i, size;
  MDR *new;

  /* Special cases... */

  if (EQNULL(m1) && EQNULL(m2))
  {
    new = mdr_Create (0, 0);
    return (new);
  }
  else if (EQNULL(m1))
  {
    new = mdr_Create (MNR(m2),MNC(m2));
    mdr_Ones(new);
    return (new);
  }
  else if (EQNULL(m2))
  {
    new = mdr_Create (MNR(m1),MNC(m1));
    mdr_Ones(new);
    return (new);
  }
  else if ( EQSCAL(m1) )
  {
    size = MNR (m2) * MNC (m2);
    new  = mdr_Create (MNR (m2), MNC (m2));
    for (i=0; i<size; i++)
    {
      if ( EQCMPL(MdcV0(m1, 0),MdcV0(m2, i)) )
        MdrV0 (new, i) = 0.0;
      else
        MdrV0 (new, i) = 1.0;
    }
    return (new);
  }
  else if ( EQSCAL(m2) )
  {
    size = MNR (m1) * MNC (m1);
    new = mdr_Create (MNR (m1), MNC (m1));
    for (i=0; i<size; i++)
    {
      if ( EQCMPL(MdcV0(m1, i), MdcV0(m2, 0)) )
        MdrV0 (new, i) = 0.0;
      else
        MdrV0 (new, i) = 1.0;
    }
    return (new);
  }
  else if ( !EQSIZE(m1,m2) )
  {
    fprintf (stderr, "matrices must have the same dimension for comparision\n");
    fprintf (stderr, "m1: nrow = %i ncol = %i m2: nrow = %i, ncol = %i\n",
             MNR (m1), MNC (m1), MNR (m2), MNC (m2));
    rerror ("matrix-compare");
  }

  size = MNR (m1) * MNC (m1);
  new = mdr_Create (MNR (m1), MNC (m1));

  for (i = 0; i < size; i++)
  {
    if ( EQCMPL(MdcV0(m1,i), MdcV0(m2,i)) )
      MdrV0 (new, i) = 0.0;
    else
      MdrV0 (new, i) = 1.0;
  }

  return (new);
}

MDR *
mdc_Lt (MDC * m1, MDC * m2)
{
  int i, size;
  MDR *new;

  /* Error check. */

  /* Special case, m1 || m2 is a 1-by-1. */
  if ( EQSCAL(m1) && NEQNULL(m2) )
  {
    size = MNR (m2) * MNC (m2);
    new = mdr_Create (MNR (m2), MNC (m2));

    for (i = 0; i < size; i++)
      MdrV0 (new, i) = ( cabs(MdcV0(m1,0)) < cabs(MdcV0(m2,i)) );

    return (new);
  }
  else if ( EQSCAL(m2) && NEQNULL(m1) )
  {
    size = MNR (m1) * MNC (m1);
    new = mdr_Create (MNR (m1), MNC (m1));
    for (i = 0; i < size; i++)
      MdrV0 (new,i) = ( cabs(MdcV0(m1,i)) < cabs(MdcV0(m2,0)) );

    return (new);
  }
  else if ( !EQSIZE(m1,m2) )
  {
    fprintf (stderr, "matrices must have the same dimension for comparision\n");
    fprintf (stderr, "m1: nrow = %i ncol = %i m2: nrow = %i, ncol = %i\n",
	     MNR (m1), MNC (m1), MNR (m2), MNC (m2));
    rerror ("matrix-compare");
  }

  size = MNR (m1) * MNC (m1);
  new = mdr_Create (MNR (m1), MNC (m1));

  for (i = 0; i < size; i++)
    MdrV0 (new,i) = ( cabs(MdcV0(m1,i)) < cabs(MdcV0(m2,i)) );

  return (new);
}

MDR *
mdc_Le (MDC * m1, MDC * m2)
{
  int i, size;
  MDR *new;

  /* Error check. */

  /* Special case, m1 || m2 is a 1-by-1. */
  if ( EQSCAL(m1) && NEQNULL(m2) )
  {
    size = MNR (m2) * MNC (m2);
    new = mdr_Create (MNR (m2), MNC (m2));

    for (i = 0; i < size; i++)
      MdrV0 (new, i) = ( cabs(MdcV0(m1,0)) <= cabs(MdcV0(m2,i)) );

    return (new);
  }
  else if ( EQSCAL(m2) && NEQNULL(m1) )
  {
    size = MNR (m1) * MNC (m1);
    new = mdr_Create (MNR (m1), MNC (m1));
    for (i = 0; i < size; i++)
      MdrV0 (new,i) = ( cabs(MdcV0(m1,i)) <= cabs(MdcV0(m2,0)) );

    return (new);
  }
  else if ( !EQSIZE(m1,m2) )
  {
    fprintf (stderr, "matrices must have the same dimension for comparision\n");
    fprintf (stderr, "m1: nrow = %i ncol = %i m2: nrow = %i, ncol = %i\n",
             MNR (m1), MNC (m1), MNR (m2), MNC (m2));
    rerror ("matrix-compare");
  }

  size = MNR (m1) * MNC (m1);
  new = mdr_Create (MNR (m1), MNC (m1));

  for (i = 0; i < size; i++)
    MdrV0 (new,i) = ( cabs(MdcV0(m1,i)) <= cabs(MdcV0(m2,i)) );

  return (new);
}

MDR *
mdc_Gt (MDC * m1, MDC * m2)
{
  int i, size;
  MDR *new;

  /* Error check. */

  /* Special case, m1 || m2 is a 1-by-1. */
  if ( EQSCAL(m1) && NEQNULL(m2) )
  {
    size = MNR (m2) * MNC (m2);
    new = mdr_Create (MNR (m2), MNC (m2));

    for (i = 0; i < size; i++)
      MdrV0 (new, i) = ( cabs(MdcV0(m1,0)) > cabs(MdcV0(m2,i)) );

    return (new);
  }
  else if ( EQSCAL(m2) && NEQNULL(m1) )
  {
    size = MNR (m1) * MNC (m1);
    new = mdr_Create (MNR (m1), MNC (m1));
    for (i = 0; i < size; i++)
      MdrV0 (new,i) = ( cabs(MdcV0(m1,i)) > cabs(MdcV0(m2,0)) );

    return (new);
  }
  else if ( !EQSIZE(m1,m2) )
  {
    fprintf (stderr, "matrices must have the same dimension for comparision\n");
    fprintf (stderr, "m1: nrow = %i ncol = %i m2: nrow = %i, ncol = %i\n",
             MNR (m1), MNC (m1), MNR (m2), MNC (m2));
    rerror ("matrix-compare");
  }

  size = MNR (m1) * MNC (m1);
  new = mdr_Create (MNR (m1), MNC (m1));

  for (i = 0; i < size; i++)
    MdrV0 (new,i) = ( cabs(MdcV0(m1,i)) > cabs(MdcV0(m2,i)) );

  return (new);
}

MDR *
mdc_Ge (MDC * m1, MDC * m2)
{
  int i, size;
  MDR *new;

  /* Error check. */

  /* Special case, m1 || m2 is a 1-by-1. */
  if ( EQSCAL(m1) && NEQNULL(m2) )
  {
    size = MNR (m2) * MNC (m2);
    new = mdr_Create (MNR (m2), MNC (m2));

    for (i = 0; i < size; i++)
      MdrV0 (new, i) = ( cabs(MdcV0(m1,0)) >= cabs(MdcV0(m2,i)) );

    return (new);
  }
  else if ( EQSCAL(m2) && NEQNULL(m1) )
  {
    size = MNR (m1) * MNC (m1);
    new = mdr_Create (MNR (m1), MNC (m1));
    for (i = 0; i < size; i++)
      MdrV0 (new,i) = ( cabs(MdcV0(m1,i)) >= cabs(MdcV0(m2,0)) );

    return (new);
  }
  else if ( !EQSIZE(m1,m2) )
  {
    fprintf (stderr, "matrices must have the same dimension for comparision\n");
    fprintf (stderr, "m1: nrow = %i ncol = %i m2: nrow = %i, ncol = %i\n",
             MNR (m1), MNC (m1), MNR (m2), MNC (m2));
    rerror ("matrix-compare");
  }

  size = MNR (m1) * MNC (m1);
  new = mdr_Create (MNR (m1), MNC (m1));

  for (i = 0; i < size; i++)
    MdrV0 (new,i) = ( cabs(MdcV0(m1,i)) >= cabs(MdcV0(m2,i)) );

  return (new);
}

MDR *
mdc_And (MDC * m1, MDC * m2)
{
  int i, size;
  MDR *new;

  /* Error check. */

  /* Special case, m1 || m2 is a 1-by-1. */
  if ( EQSCAL(m1) && NEQNULL(m2) )
  {
    size = MNR (m2) * MNC (m2);
    new = mdr_Create (MNR (m2), MNC (m2));
    if(CMPLZERO(m1))
      mdr_Zero(new);
    else
      for (i = 0; i < size; i++)
        MdrV0 (new, i) = ( ! CMPLZERO(MdcV0(m2,i)) );

    return (new);
  }
  else if ( EQSCAL(m2) && NEQNULL(m1) )
  {
    size = MNR (m1) * MNC (m1);
    new = mdr_Create (MNR (m1), MNC (m1));

    if(CMPLZERO(m2))
      mdr_Zero(new);
    else
      for (i = 0; i < size; i++)
        MdrV0 (new, i) = ( ! CMPLZERO(MdcV0(m1,i)) );

    return (new);
  }
  else if ( !EQSIZE(m1,m2) )
  {
    fprintf (stderr, "matrices must have the same dimension for comparision\n");
    fprintf (stderr, "m1: nrow = %i ncol = %i m2: nrow = %i, ncol = %i\n",
             MNR (m1), MNC (m1), MNR (m2), MNC (m2));
    rerror ("matrix-and");
  }

  size = MNR (m1) * MNC (m1);
  new = mdr_Create (MNR (m1), MNC (m1));

  for (i = 0; i < size; i++)
    MdrV0 (new,i) = ((!CMPLZERO(MdcV0(m1,i))) && (!CMPLZERO(MdcV0(m2,i))));

  return (new);
}

MDR *
mdc_Or (MDC * m1, MDC * m2)
{
  int i, size;
  MDR *new;

  /* Error check. */

  /* Special case, m1 || m2 is a 1-by-1. */
  if ( EQSCAL(m1) && NEQNULL(m2) )
  {
    size = MNR (m2) * MNC (m2);
    new = mdr_Create (MNR (m2), MNC (m2));
    if(!CMPLZERO(m1))
      mdr_Ones(new);
    else
      for (i = 0; i < size; i++)
        MdrV0 (new, i) = ( ! CMPLZERO(MdcV0(m2,i)) );

        return (new);
  }
  else if ( EQSCAL(m2) && NEQNULL(m1) )
  {
    size = MNR (m1) * MNC (m1);
    new = mdr_Create (MNR (m1), MNC (m1));

    if(!CMPLZERO(m2))
      mdr_Ones(new);
    else
      for (i = 0; i < size; i++)
        MdrV0 (new, i) = ( ! CMPLZERO(MdcV0(m1,i)) );

        return (new);
  }
  else if ( !EQSIZE(m1,m2) )
  {
    fprintf (stderr, "matrices must have the same dimension for comparision\n");
    fprintf (stderr, "m1: nrow = %i ncol = %i m2: nrow = %i, ncol = %i\n",
             MNR (m1), MNC (m1), MNR (m2), MNC (m2));
    rerror ("matrix-or");
  }

  size = MNR (m1) * MNC (m1);
  new = mdr_Create (MNR (m1), MNC (m1));

  for (i = 0; i < size; i++)
    MdrV0 (new,i) = ((!CMPLZERO(MdcV0(m1,i))) || (!CMPLZERO(MdcV0(m2,i))));

  return (new);
}

MDR *
mdc_Not (MDC * m1)
{
  int i, size;
  MDR *new;

  if ( EQNULL(m1) )
    rerror ("matrix-not:matrix cannot be zero-size\n");

  size = MNR (m1) * MNC (m1);
  new = mdr_Create (MNR (m1), MNC (m1));

  for (i=0; i<size; i++)
    MdrV0 (new, i) = !(CMPLZERO(MdcV0(m1,i)));

  return (new);
}

/* **************************************************************
 * Functions to handle comparison between real and complex.
 * ************************************************************** */

MDR *
mdr_mdc_Eq (MDR * m1, MDC * m2)
{
  MDC *tmp;
  MDR *new;

  tmp = mdr_coerce_mdc (m1);
  new = mdc_Eq (tmp, m2);
  mdc_Destroy (tmp);
  return (new);
}

MDR *
mdc_mdr_Eq (MDC * m1, MDR * m2)
{
  MDC *tmp;
  MDR *new;

  tmp = mdr_coerce_mdc (m2);
  new = mdc_Eq (m1, tmp);
  mdc_Destroy (tmp);
  return (new);
}

MDR *
mdr_mdc_Ne (MDR * m1, MDC * m2)
{
  MDC *tmp;
  MDR *new;

  tmp = mdr_coerce_mdc (m1);
  new = mdc_Ne (tmp, m2);
  mdc_Destroy (tmp);
  return (new);
}

MDR *
mdc_mdr_Ne (MDC * m1, MDR * m2)
{
  MDC *tmp;
  MDR *new;

  tmp = mdr_coerce_mdc (m2);
  new = mdc_Ne (m1, tmp);
  mdc_Destroy (tmp);
  return (new);
}

MDR *
mdr_mdc_Lt (MDR * m1, MDC * m2)
{
  MDC *tmp;
  MDR *new;

  tmp = mdr_coerce_mdc (m1);
  new = mdc_Lt (tmp, m2);
  mdc_Destroy (tmp);
  return (new);
}

MDR *
mdc_mdr_Lt (MDC * m1, MDR * m2)
{
  MDC *tmp;
  MDR *new;

  tmp = mdr_coerce_mdc (m2);
  new = mdc_Lt (m1, tmp);
  mdc_Destroy (tmp);
  return (new);
}

MDR *
mdr_mdc_Le (MDR * m1, MDC * m2)
{
  MDC *tmp;
  MDR *new;

  tmp = mdr_coerce_mdc (m1);
  new = mdc_Le (tmp, m2);
  mdc_Destroy (tmp);
  return (new);
}

MDR *
mdc_mdr_Le (MDC * m1, MDR * m2)
{
  MDC *tmp;
  MDR *new;

  tmp = mdr_coerce_mdc (m2);
  new = mdc_Le (m1, tmp);
  mdc_Destroy (tmp);
  return (new);
}

MDR *
mdr_mdc_Gt (MDR * m1, MDC * m2)
{
  MDC *tmp;
  MDR *new;

  tmp = mdr_coerce_mdc (m1);
  new = mdc_Gt (tmp, m2);
  mdc_Destroy (tmp);
  return (new);
}

MDR *
mdc_mdr_Gt (MDC * m1, MDR * m2)
{
  MDC *tmp;
  MDR *new;

  tmp = mdr_coerce_mdc (m2);
  new = mdc_Gt (m1, tmp);
  mdc_Destroy (tmp);
  return (new);
}

MDR *
mdr_mdc_Ge (MDR * m1, MDC * m2)
{
  MDC *tmp;
  MDR *new;

  tmp = mdr_coerce_mdc (m1);
  new = mdc_Ge (tmp, m2);
  mdc_Destroy (tmp);
  return (new);
}

MDR *
mdc_mdr_Ge (MDC * m1, MDR * m2)
{
  MDC *tmp;
  MDR *new;

  tmp = mdr_coerce_mdc (m2);
  new = mdc_Ge (m1, tmp);
  mdc_Destroy (tmp);
  return (new);
}

MDR *
mdr_mdc_And (MDR * m1, MDC * m2)
{
  MDC *tmp;
  MDR *new;

  tmp = mdr_coerce_mdc (m1);
  new = mdc_And (tmp, m2);
  mdc_Destroy (tmp);
  return (new);
}

MDR *
mdc_mdr_And (MDC * m1, MDR * m2)
{
  MDC *tmp;
  MDR *new;

  tmp = mdr_coerce_mdc (m2);
  new = mdc_And (m1, tmp);
  mdc_Destroy (tmp);
  return (new);
}

/* **************************************************************
 * Get the size of a COMPLEX matrix for size()...
 * ************************************************************** */

MDR *
mdc_Size_BF (MDC * m)
{
  MDR *size = mdr_Create (1, 2);
  Mdr0 (size, 0, 0) = (double) MNR (m);
  Mdr0 (size, 0, 1) = (double) MNC (m);
  return (size);
}

MDR *
mdc_Length_BF (MDC * m)
{
  MDR *size = mdr_Create (1, 1);
  Mdr0 (size, 0, 0) = MAX (MNR (m), MNC (m));
  return (size);
}

/* **************************************************************
 * REAL-DENSE cross-max...
 * ************************************************************** */

MDC *
mdr_mdc_Max2 (MDR * m1, MDC * m2)
{
  MDC *tmp;
  MDC *new;

  tmp = mdr_coerce_mdc (m1);
  new = mdc_Max2 (tmp, m2);
  mdc_Destroy (tmp);
  return (new);
}

MDC *
mdc_mdr_Max2 (MDC * m1, MDR * m2)
{
  MDC *tmp;
  MDC *new;

  tmp = mdr_coerce_mdc (m2);
  new = mdc_Max2 (m1, tmp);
  mdc_Destroy (tmp);
  return (new);
}

/* **************************************************************
 * REAL-DENSE cross-min...
 * ************************************************************** */

MDC *
mdr_mdc_Min2 (MDR * m1, MDC * m2)
{
  MDC *tmp;
  MDC *new;

  tmp = mdr_coerce_mdc (m1);
  new = mdc_Min2 (tmp, m2);
  mdc_Destroy (tmp);
  return (new);
}

MDC *
mdc_mdr_Min2 (MDC * m1, MDR * m2)
{
  MDC *tmp;
  MDC *new;

  tmp = mdr_coerce_mdc (m2);
  new = mdc_Min2 (m1, tmp);
  mdc_Destroy (tmp);
  return (new);
}

MDC *
mdr_mdc_Solve (MDR * a, MDC * b, char *type)
{
  MDC *m, *tmp;

  tmp = mdr_coerce_mdc (a);
  m = mdc_Solve (tmp, b, type);
  mdc_Destroy (tmp);

  return (m);
}

MDC *
mdc_mdr_Solve (MDC * a, MDR * b, char *type)
{
  MDC *m, *tmp;

  tmp = mdr_coerce_mdc (b);
  m = mdc_Solve (a, tmp, type);
  mdc_Destroy (tmp);

  return (m);
}

MDC *
mdr_mdc_Mod (MDR * a, MDC * b)
{
  MDC *m, *tmp;

  tmp = mdr_coerce_mdc (a);
  m = mdc_Mod (tmp, b);
  mdc_Destroy (tmp);

  return (m);
}

MDC *
mdc_mdr_Mod (MDC * a, MDR * b)
{
  MDC *m, *tmp;

  tmp = mdr_coerce_mdc (b);
  m = mdc_Mod (a, tmp);
  mdc_Destroy (tmp);

  return (m);
}
