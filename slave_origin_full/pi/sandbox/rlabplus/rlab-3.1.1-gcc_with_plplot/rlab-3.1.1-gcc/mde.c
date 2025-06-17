/* mde.c Matrix Dense Entity */

/*  This file is a part of rlabplus
   Copyright (C) 2015 Marijan Kostrun

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
#include "rfileio.h"
#include "mem.h"
#include "mde.h"
#include "mdr.h"
#include "mdrf2.h"
#include "mdc.h"
#include "mds.h"
#include "util.h"
#include "bltin1.h"
#include "mathl.h"
#include "class.h"

#include "fi.h"
#include "blas.h"
#include "lp.h"

#include <stdio.h>
#include <math.h>

// #include <termcap.h>
// #include <termios.h>
#include <sys/fcntl.h>
#include <unistd.h>
#include <sys/ioctl.h>

#include <gsl/gsl_math.h>

#include "rlab_solver_parameters_names.h"

char * mde_Class (MDE * m)
{
  return (cpstr (RLAB_CLASS_MDE));
}

MDR * mde_Length_BF (MDE * m)
{
  MDR *size = mdr_Create (1, 1);
  MdrV0 (size, 0) = MAX(MNR (m), MNC (m));
  return (size);
}

MDR * mde_Size_BF (MDE * m)
{
  int nr = MNR (m);
  int nc = MNC (m);

  MDR *size = mdr_Create (1, 2);
  Mdr0 (size, 0, 0) = (double) nr;
  Mdr0 (size, 0, 1) = (double) nc;

  return (size);
}


MDE * mde_Create (int nrow, int ncol)
{
  int i;
  MDE *new=0;

  if (nrow < 0 || ncol < 0)
    rerror ("cannot specify a negative matrix dimension");

  new = (MDE *) GC_MALLOC (sizeof (MDE));
  if (!new)
    rerror ("out of memory");

  new->nrow = nrow;
  new->ncol = ncol;
  new->type = RLAB_TYPE_ENTITY;

  if (nrow * ncol != 0)
  {
    MDEPTR(new) = (Ent **) GC_MALLOC (sizeof (Ent *) * (nrow * ncol));
    if (!MDEPTR(new))
      rerror ("out of memory");
    for (i = 0; i < nrow * ncol; i++)
    {
      MdeV0(new,i) = ent_Create();
    }
  }
  else
  {
    new->nrow = 0;
    new->ncol = 0;
    MDEPTR(new) = 0;
  }
  new->list = 0;

  return (new);
}

MDE * mde_CreateEmpty (int nrow, int ncol)
{
  int i;
  MDE *new=0;

  if (nrow < 0 || ncol < 0)
    rerror ("cannot specify a negative matrix dimension");

  new = (MDE *) GC_MALLOC (sizeof (MDE));
  if (!new)
    rerror ("out of memory");

  new->nrow = nrow;
  new->ncol = ncol;
  new->type = RLAB_TYPE_ENTITY;

  if (nrow * ncol != 0)
  {
    MDEPTR(new) = (Ent **) GC_MALLOC (sizeof (Ent *) * (nrow * ncol));
    if (!MDEPTR(new))
      rerror ("out of memory");
    for (i = 0; i < nrow * ncol; i++)
    {
      MdeV0(new,i) = 0;
    }
  }
  else
  {
    new->nrow = 0;
    new->ncol = 0;
    MDEPTR(new) = 0;
  }
  new->list = 0;

  return (new);
}

MDE * mde_CreateDouble (double val)
{
  MDR *new;
  new = mde_Create (1, 1);
  MdeV0 (new, 0) = (Ent *) ent_Create_Rlab_Double(val);
  return (new);
}


MDE * mde_CreateInt32 (int val)
{
  MDR *new;
  new = mde_Create (1, 1);
  MdeV0 (new, 0) = (Ent *) ent_Create_Rlab_Int(val);
  return (new);
}

void * mde_Destroy (MDE * m)
{
  int i;

  if (!m)
    return (0);

  for (i=0; i<SIZE(m); i++)
  {
    if (!MdeV0(m,i))
      continue;
    ent_Destroy (MdeV0(m,i));
    MdeV0(m,i)=0;
  }

  if (MDPTR(m))
    GC_FREE (MDPTR(m));
  MDPTR(m)=0;

  m->nrow = -1;
  m->ncol = -1;

  if (m->list)
  {
    btree_Destroy (m->list);
    m->list = 0;
  }

  GC_FREE (m);
  return (0);
}

void * mde_DestroyEmpty (MDE * m)
{
  // something to destroy?
  if (!m)
    return 0;

  m->nrow = -1;
  m->ncol = -1;

  MDPTR(m) = 0;

  if (m->list)
  {
    btree_Destroy (m->list);
    m->list = 0;
  }

  GC_FREE (m);
  m = 0;

  return 0;
}


//
// copy of entity matrix
// 
MDE * mde_Copy (MDE * m)
{
  MDE *new=0;
  int i, size;

  new = mde_Create (MNR (m), MNC (m));

  size = MNR (m) * MNC (m);
  for (i = 0; i < size; i++)
  {
    if (MdeV0(m,i))
    {
      MdeV0(new,i) = ent_Copy ( MdeV0(m,i) );
    }
    else
    {
      MdeV0(new,i) = ent_Create();
    }
  }

  /* Copy the list, if there is one. */
  if (m->list)
  {
    new->list = btree_Copy (m->list);
  }

  return (new);
}

//
// duplicate of entity matrix
// 
MDE * mde_Duplicate(MDE * m)
{
  MDE *new=0;
  int i, size;

  new = mde_Create (MNR (m), MNC (m));

  size = MNR (m) * MNC (m);
  for (i = 0; i < size; i++)
  {
    if (MdeV0(m,i))
    {
      MdeV0(new,i) = ent_Duplicate ( MdeV0(m,i) );
    }
    else
    {
      MdeV0(new,i) = ent_Create();
    }
  }

  /* Copy the list, if there is one. */
  if (m->list)
  {
    new->list = btree_Copy (m->list);
  }

  return (new);
}


/* **************************************************************
 * Coerce matrix indices into ints for use later...
 * Return an int array, the 1st element is the number of elements.
 * ************************************************************** */

int * mde_IndCoerceInt (MDE * m, MDR * i)
{
  int *ind, k, size = SIZE(i);

  if (size<1)
  {
    return (0);
  }

  ind = (int *) GC_malloc_atomic_ignore_off_page ((size + 1) * sizeof (int));
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
 * Vector Sub-Expression
 * ************************************************************** */
MDE * mde_VectorSub (MDE * m, int *i, int *type)
{
  int j, size, msize;
  MDE *new;

  *type = MATRIX_DENSE_ENTITY;

  /* Handle empty matrix indices. */

  if (i == 0)
  {
    new = mde_Create (0, 0);
    return (new);
  }

  size = i[0];
  msize = SIZE(m);
  if (!msize)
  {
    new = mde_Create (0, 0);
    return (new);
  }

  if (MNC (m) == 1)
    new = mde_Create (size, 1);
  else
    new = mde_Create (1, size);

  for (j = 1; j <= size; j++)
  {
    if (i[j] > msize)
    {
      mde_Destroy (new);
      fprintf (stderr, "\tindex exceeds matrix limits\n");
      fprintf (stderr, "\tindex value: %i, matrix size: %i\n", i[j], msize);
      rerror ("sub-matrix evaluation");
    }
    if (MdeV1(m, i[j]))
      MdeV1 (new, j) = ent_Copy(MdeV1(m, i[j]));
    else
      MdeV1 (new, j) = ent_Create();
  }

  return (new);
}


int mde_Size (MDE * m)
{
  return (SIZE(m));
}

char ** mde_Members (MDE * m, int *n)
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

void * mde_MemberRef (MDE * m, char *name, int *type)
{
  Ent *ne;
  void *rptr;

  if (!strcmp (name, RLAB_MEMBER_NROW))
  {
    ne = ent_Create_Rlab_Double ((double) MNR (m));
    rptr = (void *) ne;
    *type = ENTITY;
  }
  else if (!strcmp (name, RLAB_MEMBER_NCOL))
  {
    ne = ent_Create_Rlab_Double ((double) MNC (m));
    rptr = (void *) ne;
    *type = ENTITY;
  }
  else if (!strcmp (name, RLAB_MEMBER_SIZE))
  {
    ne = ent_Create_Rlab_Double ((double) SIZE(m));
    rptr = (void *) ne;
    *type = ENTITY;
  }
  else if (!strcmp (name, RLAB_MEMBER_CLASS))
  {
    ne = ent_Create_Rlab_String ( RLAB_CLASS_MDE );
    rptr = (void *) ne;
    *type = ENTITY;
  }
  else if (!strcmp (name, RLAB_MEMBER_TYPE))
  {
    ne = ent_Create_Rlab_String ( RLAB_TYPE_MDE );
    rptr = (void *) ne;
    *type = ENTITY;
  }
  else if (!strcmp (name, RLAB_MEMBER_STORAGE))
  {
    ne = ent_Create_Rlab_String ( RLAB_STORAGE_DENSE );
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







