/*
 * ent.c
 * Some pieces of code borrowed from the GNU C-library
 * (some of which, in turn, were borrowed from UCB)
 */

/*  This file is a part of RLaB ("Our"-LaB)
   Copyright (C) 1995 Ian R. Searle
   Copyright (C) 2018 Marijan Kostrun

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
#include "ent.h"
#include "class.h"
#include "mem.h"
#include "util.h"
#include "mds.h"
#include "mdr.h"
#include "mdc.h"
#include "mde.h"
#include "mathl.h"
#include "rlab_solver_parameters_names.h"

int not_ent_double_vector(Ent *e)
{
  MDR *y=0;

  if (!e)
    return 1; // empty entity

  if (ent_type (e) != MATRIX_DENSE_REAL)
    return 2; // not MDR

  y = ent_data(e);
  if (SIZE(y)<1)
    return 2; // null size MDR

  if (!EQVECT(y))
    return 3; // dense matrix, and not vector

  return 0;
}

int ismde(Ent * x)
{
  MD *m=0;
  int rval=0, i;
  if (x)
  {
    if (ent_type (x) == MATRIX_DENSE_ENTITY)
    {
      rval = 1;
      m = ent_data(x);
      for (i=0; i<SIZE(m); i++)
      {
        Ent *e = MdeV0(m,i);
        if (ent_type (e) == MATRIX_DENSE_ENTITY)
        {
          rval = MAX(rval, 1 + ismde(e));
          break;
        }
      }
    }
  }
  return rval;
}

int isfuncent(Ent * x)
{
  int rval= 0;

  if (x)
  {
    if (ent_type (x) == U_FUNCTION || ent_type (x) == BLTIN)
      rval = 1;
  }

  return rval;
}

int isdensematrix(Ent *x)
{
  int rval= 0;
  if (x)
  {
    if (      (ent_type(x) == MATRIX_DENSE_REAL) || (ent_type(x) == MATRIX_DENSE_STRING)
          ||  (ent_type(x) == MATRIX_DENSE_COMPLEX)    )
      rval = 1;
  }
  return rval;
}

Ent * ent_Create (void)
{
  Ent *ent = (Ent *) GC_MALLOC (sizeof (Ent));
  if (!ent)
    rerror ("out of memory");

  ent->refc = 0;
  ent->type = 0;
  ent->d = 0.0;
  ent->data = 0;

  return (ent);
}

Ent * ent_Create_Rlab_Success (void)
{
  Ent *rent = ent_Create ();
  ent_data (rent) = mdr_CreateScalar (RLAB_STATUS_SUCCESS);
  ent_type (rent) = MATRIX_DENSE_REAL;
  return (rent);
}

Ent * ent_Create_Rlab_Failure (void)
{
  Ent *rent = ent_Create ();
  ent_data (rent) = mdr_CreateScalar (RLAB_STATUS_FAILURE);
  ent_type (rent) = MATRIX_DENSE_REAL;
  return (rent);
}

Ent * ent_Create_Rlab_Error (void)
{
  Ent *rent = ent_Create ();
  ent_data (rent) = mdr_CreateScalar (RLAB_STATUS_ERROR);
  ent_type (rent) = MATRIX_DENSE_REAL;
  return (rent);
}

Ent *
ent_Create_Rlab_Complex (double re, double im)
{
  Ent *rent = ent_Create ();
  ent_data (rent) = mdc_CreateScalar (re, im);
  ent_type (rent) = MATRIX_DENSE_COMPLEX;
  return (rent);
}

Ent *
ent_Create_Rlab_Double (double x)
{
  Ent *rent = ent_Create ();
  ent_data (rent) = mdr_CreateScalar (x);
  ent_type (rent) = MATRIX_DENSE_REAL;
  return (rent);
}

Ent * ent_Create_Rlab_String (char * s)
{
  Ent *rent = ent_Create ();
  ent_data (rent) = mds_CreateScalar (s);
  ent_type (rent) = MATRIX_DENSE_STRING;
  return (rent);
}

Ent * ent_Assign_Rlab_String (char * s)
{
  Ent *rent = ent_Create ();
  ent_data (rent) = mds_AssignScalar(s);
  ent_type (rent) = MATRIX_DENSE_STRING;
  return (rent);
}

Ent *
ent_Create_Rlab_Int (int x)
{
  Ent *rent = ent_Create ();
  ent_data (rent) = mdi_CreateScalar (x);
  ent_type (rent) = MATRIX_DENSE_REAL;
  return (rent);
}

Ent *
ent_Assign_Rlab_MDR (MDR *x)
{
  Ent *rent = ent_Create ();
  if (x)
    ent_data (rent) = x;
  else
    ent_data (rent) = mdr_Create(0,0);
  ent_type (rent) = MATRIX_DENSE_REAL;
  return (rent);
}

Ent *
ent_Assign_Rlab_MDC (MDC *x)
{
  Ent *rent = ent_Create ();

  if (x)
    ent_data (rent) = x;
  else
    ent_data (rent) = mdc_Create(0,0);
  ent_type (rent) = MATRIX_DENSE_COMPLEX;
  return (rent);
}

Ent *
ent_Assign_Rlab_MDS (MDS *x)
{
  Ent *rent = ent_Create ();
  if (x)
    ent_data (rent) = x;
  else
    ent_data (rent) = mds_Create(0,0);
  ent_type (rent) = MATRIX_DENSE_STRING;
  return (rent);
}

Ent *
ent_Assign_Rlab_MDE (MDE *x)
{
  Ent *rent = ent_Create ();
  if (x)
    ent_data (rent) = x;
  else
    ent_data (rent) = mde_Create(0,0);
  ent_type (rent) = MATRIX_DENSE_ENTITY;
  return (rent);
}


Ent *
ent_Assign_Rlab_BTREE (Btree *b)
{
  Ent *rent = ent_Create ();
  if (b)
    ent_data (rent) = b;
  else
    ent_data (rent) = btree_Create();
  ent_type (rent) = BTREE;
  return (rent);
}

Ent *
ent_Assign_Rlab_Rtype (void *x, int rtype)
{
  Ent *rent=0;
  if (rtype!=UNDEF)
  {
    rent = ent_Create ();
    ent_data (rent) = x;
    ent_type (rent) = rtype;
  }
  else
  {
    rent = ent_Create ();
    ent_data (rent) = mdr_Create(0,0);
    ent_type (rent) = MATRIX_DENSE_REAL;
  }

  return (rent);
}

/*
 * Destroy an entity. Just decrement the reference count.
 * Do not actually delete the entity until the reference
 * count equals zero.
 */

int ent_Destroy (Ent * ent)
{
  if (!ent)			/* Make sure there is an entity present. */
    return (0);

  if (ent_type(ent)==UNDEF)
    return (1);

  ent->refc--;		/* Decrement the ref count. */
  if (ent->refc <= 0)
  {
    if(ent_type(ent) == BTREE)
    {
      Btree * btree = ent_data(ent);
      if(btree->isconst)
      {
        ent->refc++;
        rerror("internal error: ent_Destroy() on a protected entity-list!");
      }
    }
    class_destroy (ent);
    GC_free (ent);
    return (0);
  }
  
  return (2);
}

/*
 * Destroy an entity IFF the reference count is zero.
 */

int ent_Clean (Ent * ent)
{
  // cannot clean: nothing to clean
  if(!ent)
    return (-2);

  // cannot clean: entity is undefined
//   if (ent_type(ent)==UNDEF)
//     return (-3);

  if (ent->refc == 0)
  {
    if(ent_type(ent) == BTREE)
    {
      Btree * btree = ent_data(ent);
      if(btree->isconst)
      {
//         fprintf(stderr, "ent_Clean: cannot clean protected entity-list!");
        return (-1);
      }
    }
    class_destroy (ent);
    GC_FREE (ent);

    // entity cleaned
    return (1);
  }

  // cannot clean because its reference count is not 0
  return (0);
}

/*
 * Destroy an entity that contains no data
 * other than the double value. This is primarily
 * for use with the class functions. It is nice to
 * have a function to call in this case, rather than
 * get an obnoxious error message.
 */

void
ent_double_Destroy (Ent * e)
{
  /* Do nothing. */
}

void
ent_undef_Destroy (Ent * e)
{
  /* Do nothing. */
}

void
ent_undef_Print (Ent * e, FILE * fptr)
{
  fprintf (fptr, "\tUNDEFINED\n");
}

void
ent_SetType (Ent * ent, int type)
{
  ent->type = type;
}

Ent *
ent_Copy (Ent * ent)
{
  ent->refc++;
  return (ent);
}

/* **************************************************************
 * Duplicate an entity. Return a new entity, and a COPY of its
 * data, IFF the reference count is > 1.
 * ************************************************************** */
Ent *
ent_Duplicate_old (Ent * ent)
{
  Ent *new;

  if (ent->refc > 1)
  {
    new = class_copy (ent);
    return (new);
  }

  return (ent);
}

// **************************************************************
// Duplicate an entity. Return a new entity, and a COPY of its
// data.
// **************************************************************
Ent *
ent_Duplicate (Ent * ent)
{
  Ent *new = class_copy (ent);
  new->refc = 1;
  return (new);
  }

Ent *
ent_Inc (Ent * ent)
{
  if (ent->refc <= 1)
  {
    ent->d = ent->d + 1.0;
    return (ent);
  }
  else
  {
    Ent *nent = ent_Copy (ent);
    nent->d = nent->d + 1.0;
    return (nent);
  }
}

Ent *
ent_Dec (Ent * ent)
{
  if (ent->refc <= 1)
  {
    ent->d = ent->d - 1.0;
    return (ent);
  }
  else
  {
    Ent *nent = ent_Copy (ent);
    nent->d = nent->d - 1.0;
    return (nent);
  }
}

void
ent_Double (Ent * ent, double d)
{
  ent->d = d;
}

char *
ent_Double_CharPointer (Ent * e)
{
  char *cp, tmp[100];

  sprintf (tmp, "%.6g", e->d);
  cp = cpstr (tmp);

  return (cp);
}

char *
double_Class (void *ptr)
{
  return cpstr ("num");
}

char *
undef_Class (void *ptr)
{
  return cpstr ("UNDEF");
}

MDS *
double_Type_BF (void *ptr)
{
  MDS *type = mds_CreateScalar ( "real" );
  return (type);
}

size_t
double_Sizeof (void *ptr)
{
  /* Don't de-ref the pointer ! */
  return (sizeof (double));
}

MDS *
undef_Type_BF (void *ptr)
{
  MDS *type = mds_CreateScalar ( "UNDEF" );
  return (type);
}

void *
undef_Copy (void *ptr)
{
  return (0);
}
