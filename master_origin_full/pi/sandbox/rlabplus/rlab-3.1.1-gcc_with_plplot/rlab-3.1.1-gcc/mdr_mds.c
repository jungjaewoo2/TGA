/* mdr_mds.c Matrix Dense Real and Matrix Dense Complex interaction. */

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
#include "mds.h"
#include "util.h"
#include "mathl.h"
#include <stdio.h>
#include <math.h>

/*
 * Used with string comparisons.
 */

static char null_str[] = { '\0' };

/* **************************************************************
 * Coerce matrix indices into ints for use later...
 * Return an int array, the 1st element is the number of elements.
 * ************************************************************** */

int *
mds_IndCoerceInt (MDS * m, MDR * i)
{
  int *ind, k, size;

  if ((size = MNR (i) * MNC (i)) == 0)
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
 * Coerce a MDR into a MDS.
 * ************************************************************** */

MDS *
mdr_coerce_mds (MDR * m)
{
  MDS *new = 0;

  if ((MNR (m) == 0) && (MNC (m) == 0))
  {
    new = mds_Create (0, 0);
    return (new);
  }
  else
  {
    fprintf (stderr, "\tcannot coerce numerical matrix to string matrix\n");
    rerror ("Matrix-Dense-Real coerce Matrix-Dense-String");
  }
  return (new);			/* Shut up the compiler. */
}

/* **************************************************************
 * Assign MDR[i;j] = [MDS]
 * ************************************************************** */

MDS *
mdr_mds_MatrixAssign (MDR * var, int *i, int *j, MDS * rhs)
{
  MDS *stmp;

  stmp = mdr_coerce_mds (var);

  stmp = mds_MatrixAssign (stmp, i, j, rhs);
  mdr_Destroy (var);

  return (stmp);
}

MDS *
mdr_mds_MatrixAssignR (MDR * var, int *i, MDS * rhs)
{
  MDS *stmp;

  stmp = mdr_coerce_mds (var);

  stmp = mds_MatrixAssignR (stmp, i, rhs);
  mdr_Destroy (var);

  return (stmp);
}

MDS *
mdr_mds_MatrixAssignC (MDR * var, int *j, MDS * rhs)
{
  MDS *stmp;

  stmp = mdr_coerce_mds (var);

  stmp = mds_MatrixAssignC (stmp, j, rhs);
  mdr_Destroy (var);

  return (stmp);
}

MDS *
mdr_mds_VectorAssign (MDR * var, int *j, MDS * rhs)
{
  MDS *stmp;

  stmp = mdr_coerce_mds (var);

  stmp = mds_VectorAssign (stmp, j, rhs);
  mdr_Destroy (var);

  return (stmp);
}

/* **************************************************************
 * Append matrices [m1 , m2]
 * ************************************************************** */

MDS *
mdr_mds_Append (MDR * m1, MDS * m2)
{
  MDS *new;
  MDS *mtmp = mdr_coerce_mds (m1);
  new = mds_Append (mtmp, m2);
  mds_Destroy (mtmp);
  return (new);
}

MDS *
mds_mdr_Append (MDS * m1, MDR * m2)
{
  MDS *new;
  MDS *mtmp = mdr_coerce_mds (m2);
  new = mds_Append (m1, mtmp);
  mds_Destroy (mtmp);
  return (new);
}

/* **************************************************************
 * Stack matrices [m1 , m2]
 * ************************************************************** */

MDS *
mdr_mds_Stack (MDR * m1, MDS * m2)
{
  MDS *new;
  MDS *mtmp = mdr_coerce_mds (m1);
  new = mds_Stack (mtmp, m2);
  mds_Destroy (mtmp);
  return (new);
}

MDS *
mds_mdr_Stack (MDS * m1, MDR * m2)
{
  MDS *new;
  MDS *mtmp = mdr_coerce_mds (m2);
  new = mds_Stack (m1, mtmp);
  mds_Destroy (mtmp);
  return (new);
}

/* **************************************************************
 * Compare two string matrices (m1 == m2)
 * ************************************************************** */

MDR *
mds_Eq (MDS * m1, MDS * m2)
{
  int i, size;
  MDR *new;

  /* Special case... */

  if (MNR (m1) == 0 || MNR (m2) == 0)
  {
    if (MNR (m1) == 0 && MNR (m2) == 0)
    {
      new = mdr_Create (1, 1);
      MdrV0 (new, 0) = 1.0;
      return (new);
    }
    else
    {
      new = mdr_Create (1, 1);
      MdrV0 (new, 0) = 0.0;
      return (new);
    }
  }
  else if ((MNR (m1) == 1) && (MNC (m1) == 1))
  {
    size = MNR (m2) * MNC (m2);
    new = mdr_Create (MNR (m2), MNC (m2));
    for (i = 0; i < size; i++)
    {
      if (MdsV0 (m1, 0) && MdsV0 (m2, i))
        MdrV0 (new, i) = (double) (strcmp (MdsV0 (m1, 0), MdsV0 (m2, i)) == 0);
      else
        MdrV0 (new, i) = 0;
    }
    return (new);
  }
  else if ((MNR (m2) == 1) && (MNC (m2) == 1))
  {
    size = MNR (m1) * MNC (m1);
    new = mdr_Create (MNR (m1), MNC (m1));
    for (i = 0; i < size; i++)
    {
      if (MdsV0 (m1, i) && MdsV0 (m2, 0))
        MdrV0 (new, i) = (double) (strcmp (MdsV0 (m1, i), MdsV0 (m2, 0)) == 0);
      else
        MdrV0 (new, i) = 0;
    }
    return (new);
  }
  else if ((MNR (m1) != MNR (m2)) || (MNC (m1) != MNC (m2)))
  {
    fprintf (stderr,
	     "\tmatrices must have the same dimension for comparision\n");
    fprintf (stderr, "\tm1: nrow = %i ncol = %i m2: nrow = %i, ncol = %i\n",
	     MNR (m1), MNC (m1), MNR (m2), MNC (m2));
    rerror ("matrix-compare");
  }

  size = MNR (m1) * MNC (m1);
  new = mdr_Create (MNR (m1), MNC (m1));

  for (i = 0; i < size; i++)
  {
    if (MdsV0 (m1, i) && MdsV0 (m2, i))
      MdrV0 (new, i) = (double) (strcmp (MdsV0 (m1, i), MdsV0 (m2, i)) == 0);
    else
      MdrV0 (new, i) = 0;
  }

  return (new);
}

/* **************************************************************
 * Compare a real and a string matrix (never true). But, check
 * dimensions first ???
 * ************************************************************** */

MDR *
mdr_mds_Eq (void *m1, MDS * m2)
{
  MDR *new;

  /* Special cases... */

  if (MNR (m1) == 0 || MNR (m2) == 0)
  {
    if (MNR (m1) == 0 && MNR (m2) == 0)
    {
      new = mdr_Create (1, 1);
      MdrV0 (new, 0) = 1.0;
      return (new);
    }
    else
    {
      new = mdr_Create (1, 1);
      MdrV0 (new, 0) = 0.0;
      return (new);
    }
  }
  else if ((MNR (m1) == 1) && (MNC (m1) == 1))
  {
    new = mdr_Create (MNR (m2), MNC (m2));
    mdr_Zero (new);
    return (new);
  }
  else if ((MNR (m2) == 1) && (MNC (m2) == 1))
  {
    new = mdr_Create (MNR (m1), MNC (m1));
    mdr_Zero (new);
    return (new);
  }
  else if ((MNR (m1) != MNR (m2)) || (MNC (m1) != MNC (m2)))
  {
    fprintf (stderr,
	     "\tmatrices must have the same dimension for comparision\n");
    fprintf (stderr, "\tm1: nrow = %i ncol = %i m2: nrow = %i, ncol = %i\n",
	     MNR (m1), MNC (m1), MNR (m2), MNC (m2));
    rerror ("matrix-compare");
  }

  new = mdr_Create (MNR (m1), MNC (m1));
  mdr_Zero (new);
  return (new);
}

MDR *
mds_mdr_Eq (MDS * m1, void *m2)
{
  return (mdr_mds_Eq (m2, m1));
}

/* **************************************************************
 * Compare two string matrices (m1 != m2)
 * ************************************************************** */

MDR *
mds_Ne (MDS * m1, MDS * m2)
{
  int i, size;
  MDR *new;

  /* Special case... */

  if (MNR (m1) == 0 || MNR (m2) == 0)
  {
    if (MNR (m1) == 0 && MNR (m2) == 0)
    {
      new = mdr_Create (1, 1);
      MdrV0 (new, 0) = 0.0;
      return (new);
    }
    else
    {
      new = mdr_Create (1, 1);
      MdrV0 (new, 0) = 1.0;
      return (new);
    }
  }
  else if ((MNR (m1) == 1) && (MNC (m1) == 1))
  {
    size = MNR (m2) * MNC (m2);
    new = mdr_Create (MNR (m2), MNC (m2));
    for (i = 0; i < size; i++)
    {
      if (MdsV0 (m1, 0) && MdsV0 (m2, i))
        MdrV0 (new, i) = (double) (strcmp (MdsV0 (m1, 0), MdsV0 (m2, i)) != 0);
      else
        MdrV0 (new, i) = 1;
    }
    return (new);
  }
  else if ((MNR (m2) == 1) && (MNC (m2) == 1))
  {
    size = MNR (m1) * MNC (m1);
    new = mdr_Create (MNR (m1), MNC (m1));
    for (i = 0; i < size; i++)
    {
      if (MdsV0 (m1, i) && MdsV0 (m2, 0))
        MdrV0 (new, i) = (double) (strcmp (MdsV0 (m1, i), MdsV0 (m2, 0)) != 0);
      else
        MdrV0 (new, i) = 1;
    }
    return (new);
  }
  else if ((MNR (m1) != MNR (m2)) || (MNC (m1) != MNC (m2)))
  {
    fprintf (stderr,
	     "\tmatrices must have the same dimension for comparision\n");
    fprintf (stderr, "\tm1: nrow = %i ncol = %i m2: nrow = %i, ncol = %i\n",
	     MNR (m1), MNC (m1), MNR (m2), MNC (m2));
    rerror ("matrix-compare");
  }

  size = MNR (m1) * MNC (m1);
  new = mdr_Create (MNR (m1), MNC (m1));

  for (i = 0; i < size; i++)
  {
    if (MdsV0 (m1, i) && MdsV0 (m2, i))
      MdrV0 (new, i) = (double) (strcmp (MdsV0 (m1, i), MdsV0 (m2, i)) != 0);
    else
      MdrV0 (new, i) = 1;
  }

  return (new);
}

/* **************************************************************
 * Compare a real and a string matrix (always true). But, check
 * dimensions first ???
 * ************************************************************** */

MDR *
mdr_mds_Ne (void *m1, MDS * m2)
{
  MDR *new;

  /* Special case... */

  if (MNR (m1) == 0 || MNR (m2) == 0)
  {
    if (MNR (m1) == 0 && MNR (m2) == 0)
    {
      new = mdr_Create (1, 1);
      MdrV0 (new, 0) = 0.0;
      return (new);
    }
    else
    {
      new = mdr_Create (1, 1);
      MdrV0 (new, 0) = 1.0;
      return (new);
    }
  }
  if ((MNR (m1) == 1) && (MNC (m1) == 1))
  {
    new = mdr_Create (MNR (m2), MNC (m2));
    mdr_Ones (new);
    return (new);
  }
  else if ((MNR (m2) == 1) && (MNC (m2) == 1))
  {
    new = mdr_Create (MNR (m1), MNC (m1));
    mdr_Ones (new);
    return (new);
  }
  else if ((MNR (m1) != MNR (m2)) || (MNC (m1) != MNC (m2)))
  {
    fprintf (stderr,
	     "\tmatrices must have the same dimension for comparision\n");
    fprintf (stderr, "\tm1: nrow = %i ncol = %i m2: nrow = %i, ncol = %i\n",
	     MNR (m1), MNC (m1), MNR (m2), MNC (m2));
    rerror ("matrix-compare");
  }

  new = mdr_Create (MNR (m1), MNC (m1));
  mdr_Ones (new);
  return (new);
}

MDR *
mds_mdr_Ne (MDS * m1, void *m2)
{
  return (mdr_mds_Ne (m2, m1));
}

MDR *
mds_Lt (MDS * m1, MDS * m2)
{
  int i, size;
  MDR *new;

  /* Error check. */

  /* Special case, m1 || m2 is a 1-by-1. */

  if ((MNR (m1) == 1) && (MNC (m1) == 1))
  {
    size = MNR (m2) * MNC (m2);
    new = mdr_Create (MNR (m2), MNC (m2));
    for (i = 0; i < size; i++)
    {
      MdrV0 (new, i) = (double) (strcmp (MdsV0 (m1, 0), MdsV0 (m2, i)) < 0);
    }
    return (new);
  }
  else if ((MNR (m2) == 1) && (MNC (m2) == 1))
  {
    size = MNR (m1) * MNC (m1);
    new = mdr_Create (MNR (m1), MNC (m1));
    for (i = 0; i < size; i++)
    {
      MdrV0 (new, i) = (double) (strcmp (MdsV0 (m1, i), MdsV0 (m2, 0)) < 0);
    }
    return (new);
  }

  if ((MNR (m1) != MNR (m2)) || (MNC (m1) != MNC (m2)))
  {
    fprintf (stderr,
	     "\tmatrices must have the same dimension for comparision\n");
    fprintf (stderr, "\tm1: nrow = %i ncol = %i m2: nrow = %i, ncol = %i\n",
	     MNR (m1), MNC (m1), MNR (m2), MNC (m2));
    rerror ("matrix-compare");
  }

  size = MNR (m1) * MNC (m1);
  new = mdr_Create (MNR (m1), MNC (m1));

  for (i = 0; i < size; i++)
  {
    MdrV0 (new, i) = (double) (strcmp (MdsV0 (m1, i), MdsV0 (m2, i)) < 0);
  }

  return (new);
}

MDR *
mds_Le (MDS * m1, MDS * m2)
{
  int i, size;
  MDR *new;

  /* Error check. */

  /* Special case, m1 || m2 is a 1-by-1. */

  if ((MNR (m1) == 1) && (MNC (m1) == 1))
  {
    size = MNR (m2) * MNC (m2);
    new = mdr_Create (MNR (m2), MNC (m2));
    for (i = 0; i < size; i++)
    {
      MdrV0 (new, i) = (double) (strcmp (MdsV0 (m1, 0), MdsV0 (m2, i)) <= 0);
    }
    return (new);
  }
  else if ((MNR (m2) == 1) && (MNC (m2) == 1))
  {
    size = MNR (m1) * MNC (m1);
    new = mdr_Create (MNR (m1), MNC (m1));
    for (i = 0; i < size; i++)
    {
      MdrV0 (new, i) = (double) (strcmp (MdsV0 (m1, i), MdsV0 (m2, 0)) <= 0);
    }
    return (new);
  }

  if ((MNR (m1) != MNR (m2)) || (MNC (m1) != MNC (m2)))
  {
    fprintf (stderr,
	     "\tmatrices must have the same dimension for comparision\n");
    fprintf (stderr, "\tm1: nrow = %i ncol = %i m2: nrow = %i, ncol = %i\n",
	     MNR (m1), MNC (m1), MNR (m2), MNC (m2));
    rerror ("matrix-compare");
  }

  size = MNR (m1) * MNC (m1);
  new = mdr_Create (MNR (m1), MNC (m1));

  for (i = 0; i < size; i++)
  {
    MdrV0 (new, i) = (double) (strcmp (MdsV0 (m1, i), MdsV0 (m2, i)) <= 0);
  }

  return (new);
}

MDR *
mds_Gt (MDS * m1, MDS * m2)
{
  int i, size;
  MDR *new;

  /* Error check. */

  /* Special case, m1 || m2 is a 1-by-1. */

  if ((MNR (m1) == 1) && (MNC (m1) == 1))
  {
    size = MNR (m2) * MNC (m2);
    new = mdr_Create (MNR (m2), MNC (m2));
    for (i = 0; i < size; i++)
    {
      MdrV0 (new, i) = (double) (strcmp (MdsV0 (m1, 0), MdsV0 (m2, i)) > 0);
    }
    return (new);
  }
  else if ((MNR (m2) == 1) && (MNC (m2) == 1))
  {
    size = MNR (m1) * MNC (m1);
    new = mdr_Create (MNR (m1), MNC (m1));
    for (i = 0; i < size; i++)
    {
      MdrV0 (new, i) = (double) (strcmp (MdsV0 (m1, i), MdsV0 (m2, 0)) > 0);
    }
    return (new);
  }

  if ((MNR (m1) != MNR (m2)) || (MNC (m1) != MNC (m2)))
  {
    fprintf (stderr,
	     "\tmatrices must have the same dimension for comparision\n");
    fprintf (stderr, "\tm1: nrow = %i ncol = %i m2: nrow = %i, ncol = %i\n",
	     MNR (m1), MNC (m1), MNR (m2), MNC (m2));
    rerror ("matrix-compare");
  }

  size = MNR (m1) * MNC (m1);
  new = mdr_Create (MNR (m1), MNC (m1));

  for (i = 0; i < size; i++)
  {
    MdrV0 (new, i) = (double) (strcmp (MdsV0 (m1, i), MdsV0 (m2, i)) > 0);
  }

  return (new);
}

MDR *
mds_Ge (MDS * m1, MDS * m2)
{
  int i, size;
  MDR *new;

  /* Error check. */

  /* Special case, m1 || m2 is a 1-by-1. */

  if ((MNR (m1) == 1) && (MNC (m1) == 1))
  {
    size = MNR (m2) * MNC (m2);
    new = mdr_Create (MNR (m2), MNC (m2));
    for (i = 0; i < size; i++)
    {
      MdrV0 (new, i) = (double) (strcmp (MdsV0 (m1, 0), MdsV0 (m2, i)) >= 0);
    }
    return (new);
  }
  else if ((MNR (m2) == 1) && (MNC (m2) == 1))
  {
    size = MNR (m1) * MNC (m1);
    new = mdr_Create (MNR (m1), MNC (m1));
    for (i = 0; i < size; i++)
    {
      MdrV0 (new, i) = (double) (strcmp (MdsV0 (m1, i), MdsV0 (m2, 0)) >= 0);
    }
    return (new);
  }

  if ((MNR (m1) != MNR (m2)) || (MNC (m1) != MNC (m2)))
  {
    fprintf (stderr,
	     "\tmatrices must have the same dimension for comparision\n");
    fprintf (stderr, "\tm1: nrow = %i ncol = %i m2: nrow = %i, ncol = %i\n",
	     MNR (m1), MNC (m1), MNR (m2), MNC (m2));
    rerror ("matrix-compare");
  }

  size = MNR (m1) * MNC (m1);
  new = mdr_Create (MNR (m1), MNC (m1));

  for (i = 0; i < size; i++)
  {
    MdrV0 (new, i) = (double) (strcmp (MdsV0 (m1, i), MdsV0 (m2, i)) >= 0);
  }

  return (new);
}

MDR *
mds_And (MDS * m1, MDS * m2)
{
  int i, size;
  MDR *new;

  /* Error check. */

  /* Special case, m1 || m2 is a 1-by-1. */

  if ((MNR (m1) == 1) && (MNC (m1) == 1))
  {
    size = MNR (m2) * MNC (m2);
    new = mdr_Create (MNR (m2), MNC (m2));
    for (i = 0; i < size; i++)
    {
      if (!strcmp (MdsV0 (m1, 0), null_str) &&
	  !strcmp (MdsV0 (m2, i), null_str))
      {
	MdrV0 (new, i) = 1.0;
      }
      else
      {
	MdrV0 (new, i) = 0.0;
      }
    }
    return (new);
  }
  else if ((MNR (m2) == 1) && (MNC (m2) == 1))
  {
    size = MNR (m1) * MNC (m1);
    new = mdr_Create (MNR (m1), MNC (m1));
    for (i = 0; i < size; i++)
    {
      if (!strcmp (MdsV0 (m1, i), null_str) &&
	  !strcmp (MdsV0 (m2, 0), null_str))
      {
	MdrV0 (new, i) = 1.0;
      }
      else
      {
	MdrV0 (new, i) = 0.0;
      }
    }
    return (new);
  }

  if ((MNR (m1) != MNR (m2)) || (MNC (m1) != MNC (m2)))
  {
    fprintf (stderr,
	     "\tmatrices must have the same dimension for comparision\n");
    fprintf (stderr, "\tm1: nrow = %i ncol = %i m2: nrow = %i, ncol = %i\n",
	     MNR (m1), MNC (m1), MNR (m2), MNC (m2));
    rerror ("matrix-compare");
  }

  size = MNR (m1) * MNC (m1);
  new = mdr_Create (MNR (m1), MNC (m1));

  for (i = 0; i < size; i++)
  {
    if (!strcmp (MdsV0 (m1, i), null_str) && !strcmp (MdsV0 (m2, i), null_str))
    {
      MdrV0 (new, i) = 1.0;
    }
    else
    {
      MdrV0 (new, i) = 0.0;
    }
  }

  return (new);
}

MDR *
mds_Or (MDS * m1, MDS * m2)
{
  int i, size;
  MDR *new;

  /* Error check. */

  /* Special case, m1 || m2 is a 1-by-1. */

  if ((MNR (m1) == 1) && (MNC (m1) == 1))
  {
    size = MNR (m2) * MNC (m2);
    new = mdr_Create (MNR (m2), MNC (m2));
    for (i = 0; i < size; i++)
    {
      if (!strcmp (MdsV0 (m1, 0), null_str) ||
	  !strcmp (MdsV0 (m2, i), null_str))
      {
	MdrV0 (new, i) = 1.0;
      }
      else
      {
	MdrV0 (new, i) = 0.0;
      }
    }
    return (new);
  }
  else if ((MNR (m2) == 1) && (MNC (m2) == 1))
  {
    size = MNR (m1) * MNC (m1);
    new = mdr_Create (MNR (m1), MNC (m1));
    for (i = 0; i < size; i++)
    {
      if (!strcmp (MdsV0 (m1, i), null_str) ||
	  !strcmp (MdsV0 (m2, 0), null_str))
      {
	MdrV0 (new, i) = 1.0;
      }
      else
      {
	MdrV0 (new, i) = 0.0;
      }
    }
    return (new);
  }

  if ((MNR (m1) != MNR (m2)) || (MNC (m1) != MNC (m2)))
  {
    fprintf (stderr,
	     "\tmatrices must have the same dimension for comparision\n");
    fprintf (stderr, "\tm1: nrow = %i ncol = %i m2: nrow = %i, ncol = %i\n",
	     MNR (m1), MNC (m1), MNR (m2), MNC (m2));
    rerror ("matrix-compare");
  }

  size = MNR (m1) * MNC (m1);
  new = mdr_Create (MNR (m1), MNC (m1));

  for (i = 0; i < size; i++)
  {
    if (!strcmp (MdsV0 (m1, i), null_str) || !strcmp (MdsV0 (m2, i), null_str))
    {
      MdrV0 (new, i) = 1.0;
    }
    else
    {
      MdrV0 (new, i) = 0.0;
    }
  }

  return (new);
}

MDR *
mds_Not (MDS * m1)
{
  int i, size;
  MDR *new;

  size = MNR (m1) * MNC (m1);
  new = mdr_Create (MNR (m1), MNC (m1));

  for (i = 0; i < size; i++)
  {
    if (strcmp (MdsV0 (m1, i), null_str))
      MdrV0 (new, i) = 0.0;
    else
      MdrV0 (new, i) = 1.0;
  }

  return (new);
}

/* **************************************************************
 * Get the size of a STRING matrix for size()...
 * ************************************************************** */

MDR *
mds_Size_BF (MDS * m)
{
  MDR *size = mdr_Create (1, 2);
  Mdr0 (size, 0, 0) = (double) MNR (m);
  Mdr0 (size, 0, 1) = (double) MNC (m);
  return (size);
}

MDR *
mds_Length_BF (MDS * m)
{
  return mdr_CreateScalar ( (double) MAX(MNR(m),MNC(m)) );
}
