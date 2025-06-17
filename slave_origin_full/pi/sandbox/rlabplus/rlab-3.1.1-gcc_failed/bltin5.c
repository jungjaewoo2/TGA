/* bltin5.c */

/*  This file is a part of RLaB ("Our"-LaB) + rlabplus
   Copyright (C) 2007-2016  Marijan Kostrun

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
#include "bltin.h"
#include "util.h"
#include "ent.h"
#include "class.h"
#include "rfileio.h"
#include "mathl.h"
#include "list.h"

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>

/*
 * Header files from the available classes...
 */

#include "ent.h"
#include "btree.h"
#include "bltin.h"
#include "function.h"
#include "mde.h"

#include "rlab_solver_parameters_names.h"

//
// operations on MDEs
//

// cell
static OpDef cell_method[NCL];


/* **************************************************************
 * Initialize the built-ins...
 * ************************************************************** */

void class_bltin5_init (void)
{
  //
  // cell: create MDE
  //
  cell_method[MATRIX_DENSE_ENTITY].type = MATRIX_DENSE_ENTITY;
  cell_method[MATRIX_DENSE_ENTITY].op   = (void *) mde_Create;
}

#undef  THIS_SOLVER
#define THIS_SOLVER "cell"
#undef  METHOD
#define METHOD cell_method
Ent * Cell(int nargs, Datum args[])
{
  int nr=0, nc=0, rtype=UNDEF;
  void *rval=0;
  void *(*vfptr) ();
  Ent *e1=0, *e2=0;

  if (nargs!=0 && nargs!=1 && nargs!=2)
  {
    rerror (THIS_SOLVER ": 1 or 2 arguments allowed");
  }

  // first argument:
  if (nargs == 1)
  {
    e1 = bltin_get_ent (args[0]);
    nr = MNR(e1);
    nc = MNC(e2);
  }
  else if (nargs == 2)
  {
    e1 = bltin_get_ent (args[0]);
    e2 = bltin_get_ent (args[1]);
    nr = (int) class_double(e1);
    nc = (int) class_double(e2);
  }
  nr = nr * (nc>0);
  nc = nc * (nr>0);

  rtype = METHOD[MATRIX_DENSE_ENTITY].type;
  vfptr = (VFPTR) METHOD[MATRIX_DENSE_ENTITY].op;
  if (vfptr)
    rval = (*vfptr) (nr, nc);

  ent_Clean (e1);
  ent_Clean (e2);

  return ent_Assign_Rlab_Rtype(rval, rtype);
}

#undef  THIS_SOLVER
#define THIS_SOLVER "iscell"
Ent * IsCell(int nargs, Datum args[])
{
  double rval=0;
  Ent *e1=0;

  if (nargs!=1)
  {
    rerror (THIS_SOLVER ": " RLAB_ERROR_ONE_ARG_REQUIRED "\n");
  }

  // first argument:
  e1 = bltin_get_ent (args[0]);
  rval = ismde(e1);
  ent_Clean (e1);
  return ent_Create_Rlab_Double(rval);
}

#undef  THIS_SOLVER
#define THIS_SOLVER "vcell"
Ent * VCell(int nargs, Datum args[])
{
  int nr=0, nc=0, i;
  Ent *e1=0;
  MDE *rval=0;

  nr = nargs;
  if (nr)
    nc = 1;
  rval=mde_Create(nr, nc);

  for (i=0; i<nr; i++)
  {
    e1 = bltin_get_ent (args[i]);
    MdeV0(rval,i) = ent_Copy(e1);
  }

  return ent_Assign_Rlab_MDE(rval);
}

#undef  THIS_SOLVER
#define THIS_SOLVER "hcell"
Ent * HCell(int nargs, Datum args[])
{
  int nr=0, nc=0, i;
  Ent *e1=0;
  MDE *rval=0;

  nr = nargs;
  if (nr)
    nc = 1;
  rval=mde_Create(nc, nr);

  for (i=0; i<nr; i++)
  {
    e1 = bltin_get_ent (args[i]);
    MdeV0(rval,i) = ent_Copy(e1);
  }

  return ent_Assign_Rlab_MDE(rval);
}
