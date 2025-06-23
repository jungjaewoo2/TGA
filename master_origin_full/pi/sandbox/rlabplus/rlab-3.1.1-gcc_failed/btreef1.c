/* btreef1.c */

/* Binary tree implementation using ListNode structure for tree nodes.
   The ListNode structure was originally designed for doubley linked
   lists, so for binary tree purposes, left == prev, and right == next. */

/*  This file is a part of RLaB ("Our"-LaB)
   Copyright (C) 1992, 1994  Ian R. Searle

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
#include "class.h"
#include "symbol.h"
#include "btree.h"
#include "function.h"
#include "util.h"
#include "mem.h"
#include "listnode.h"
#include "mdr.h"
#include "mds.h"

#include <stdio.h>

/* **************************************************************
 * Return the size of a list for size()...
 * ************************************************************** */

MDR *
btree_Size_BF (Btree * bt)
{
  MDR *size = mdr_Create (1, 1);
  Mdr0 (size, 0, 0) = (double) btree_CountNodes (bt);
  return (size);
}

/* **************************************************************
 * Return the type of a list for type()...
 * ************************************************************** */

MDS *
btree_Type_BF (Btree * bt)
{
  Ent *e;
  ListNode *ltype;
  MDS *type = mds_Create (1, 1);

  /* Try and report the contents of the type element */
  if ((ltype = btree_FindNode (bt, "type")))
  {
    e = var_ent (ltype);
    if (ent_type (e) == MATRIX_DENSE_STRING)
    {
      char *stmp = class_char_pointer (e);
      MdsV0 (type, 0) = cpstr (stmp);
      return (type);
    }
    else
    {
      type = mds_CreateScalar ( "" );
      return (type);
    }
  }

  type = mds_CreateScalar ("list");
  return (type);
}

Ent *
btree_protect (int nargs, Datum args[])
{
  Ent *e1 = 0, *rent;
  Btree *bt;

  //
  // Check arguments.
  //
  if (nargs != 1)
  {
    fprintf (stdout,
	     "protect: Protect a list variable from being overwritten or cleared.\n");
    fprintf (stdout, "protect: Format:\n");
    fprintf (stdout, "protect:   protect( varname ),\n");
    fprintf (stdout, "protect: See also  release() .\n");
  }
  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) == BTREE)
  {
    bt = ent_data (e1);
    bt->isconst = 1;
    ent_data (e1) = bt;
  }
  else
    rerror ("protect: Only a list can be protected!");
  rent = ent_Create ();
  ent_data (rent) = mdr_CreateScalar (1.0);
  ent_SetType (rent, MATRIX_DENSE_REAL);
  return (rent);
}


Ent *
btree_release (int nargs, Datum args[])
{
  Ent *e1 = 0, *rent;
  Btree *bt;

  //
  // Check arguments.
  //
  if (nargs != 1)
  {
    fprintf (stdout,
	     "release: Unprotect a content of a list variable so it can be cleared.\n");
    fprintf (stdout, "release: Format:\n");
    fprintf (stdout, "release:   release( varname ),\n");
    fprintf (stdout, "release: See also  release() .\n");
  }

  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) == BTREE)
  {
    bt = ent_data (e1);
    bt->isconst = 0;
    ent_data (e1) = bt;
  }
  else
    rerror ("protect: Only a list can be protected!");

  rent = ent_Create ();
  ent_data (rent) = mdr_CreateScalar (1.0);
  ent_SetType (rent, MATRIX_DENSE_REAL);
  return (rent);
}

Ent *
btree_IsConst (int nargs, Datum args[])
{
  Ent *e1 = 0, *rent;
  Btree *bt;
  MDR *w;
  //
  // Check arguments.
  //
  if (nargs != 1)
  {
    fprintf (stdout,
             "isprot: Check whether a list 'varname' is protected or not.\n");
    fprintf (stdout,
             "isprot: Format:\n");
    fprintf (stdout,
             "isprot:   i = isprot( varname ),\n");
    fprintf (stdout,
             "isprot: See also  protect(), release() .\n");
  }
  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) == BTREE)
  {
    bt = ent_data (e1);
    w = mdr_CreateScalar ((double) bt->isconst);
  }
  else
    w = mdr_CreateScalar (0.0);

  ent_Clean (e1);

  rent = ent_Create ();
  ent_data (rent) = w;
  ent_SetType (rent, MATRIX_DENSE_REAL);
  return (rent);
}
