/* bltin4.c */

/*  This file is a part of RLaB ("Our"-LaB) + rlabplus
   Copyright (C) 2007-2015  Marijan Kostrun

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

#ifdef __riscos
#include "riscos_ut.h"
#endif

/*
 * Header files from the available classes...
 */

#include "ent.h"
#include "btree.h"
#include "bltin.h"
#include "function.h"
#include "mdr.h"
#include "mdrf4.h"
#include "mdcf4.h"

#include "rlab_solver_parameters_names.h"

//
// operations
//

// mexp
static OpDef mexp1_method[NCL];
static OpDef mexp2_method[NCL][NCL];

// mpow
static OpDef mpow1_method[NCL];

/* **************************************************************
 * Initialize the built-ins...
 * ************************************************************** */

void
class_bltin4_init (void)
{
  //
  // mexp
  //
  mexp1_method[MATRIX_DENSE_REAL].type = MATRIX_DENSE_REAL;
  mexp1_method[MATRIX_DENSE_REAL].op   = (void *) mdr_mexp1;

  mexp1_method[MATRIX_DENSE_COMPLEX].type = MATRIX_DENSE_COMPLEX;
  mexp1_method[MATRIX_DENSE_COMPLEX].op   = (void *) mdc_mexp1;

  mexp2_method[MATRIX_DENSE_REAL][MATRIX_DENSE_REAL].type = MATRIX_DENSE_REAL;
  mexp2_method[MATRIX_DENSE_REAL][MATRIX_DENSE_REAL].op = (void *) mdr_mdr_mexp2;

  mexp2_method[MATRIX_DENSE_COMPLEX][MATRIX_DENSE_COMPLEX].type = MATRIX_DENSE_COMPLEX;
  mexp2_method[MATRIX_DENSE_COMPLEX][MATRIX_DENSE_COMPLEX].op = (void *) mdc_mdc_mexp2;

  mexp2_method[MATRIX_DENSE_COMPLEX][MATRIX_DENSE_REAL].type = MATRIX_DENSE_COMPLEX;
  mexp2_method[MATRIX_DENSE_COMPLEX][MATRIX_DENSE_REAL].op = (void *) mdc_mdr_mexp2;

  mexp2_method[MATRIX_DENSE_REAL][MATRIX_DENSE_COMPLEX].type = MATRIX_DENSE_COMPLEX;
  mexp2_method[MATRIX_DENSE_REAL][MATRIX_DENSE_COMPLEX].op = (void *) mdr_mdc_mexp2;

  //
  // mpow
  //
  mpow1_method[MATRIX_DENSE_REAL].type = MATRIX_DENSE_REAL;
  mpow1_method[MATRIX_DENSE_REAL].op   = (void *) mdr_mpow1;

  mpow1_method[MATRIX_DENSE_COMPLEX].type = MATRIX_DENSE_COMPLEX;
  mpow1_method[MATRIX_DENSE_COMPLEX].op   = (void *) mdc_mpow1;

}

/* **************************************************************
 * mexp function...
 * ************************************************************** */

Ent *
Mexp (int nargs, Datum args[])
{
  int typ[4], rtype=0, i, nlastarg, imet = -1;
  double t=0;
  void *rval = 0;
  void *(*vfptr) ();
  Ent *rent=0, *ent[4] = {0, 0, 0, 0};

  int ideg = 6;   // exp(A) using Pade approximant of degree 'ideg'
  int idummy;

  ListNode *node;

  if (nargs < 1 || nargs > 4 )
    rerror ("mexp: one to four arguments required");

  for (i=0; i<nargs;i++)
  {
    ent[i] = bltin_get_ent (args[i]);
    typ[i] = ent_type (ent[i]);
  }

  if (typ[0]==UNDEF)
    rerror ("mexp: one to four arguments required");

  //
  // is the last argument list of options
  //
  nlastarg = nargs;
  if (typ[nlastarg-1] == BTREE)
  {
    //
    // read in the options
    //
    // id: degree of the diagonal Pade to be used for exp(t*A) or exp(A)
    node = btree_FindNode (ent_data (ent[nlastarg-1]), "pade_degree");
    if (node != 0)
    {
      idummy = (int) class_double (var_ent (node));
      if (idummy > 0 || idummy < 21)
        ideg = idummy;
    }
    // imet: choice of method for exp(t*A)*v
    //  imet = -1 for chebyshev
    //  imet > 0 for krylov where imet is the size of the krylov
    //  subspace (imet < n is chacked in low level routine)
    node = btree_FindNode (ent_data (ent[nlastarg-1]), "krylov_dim");
    if (node != 0)
    {
      idummy = (int) class_double (var_ent (node));
      if (idummy == -1 || (idummy > 0))
        imet = idummy;
    }

    --nlastarg ;
  }

  if (nlastarg == 3)
  {
    //
    // mexp(A,t,v)
    //
    if (typ[1]==MATRIX_DENSE_REAL)
      t = class_double(ent[1]);
    else
      rerror("mexp: incorrect second argument");

    rtype = mexp2_method[typ[0]][typ[2]].type;
    vfptr = (VFPTR) mexp2_method[typ[0]][typ[2]].op;

    if (vfptr == 0)
    {
      fprintf (stderr, "Entity type: %s and %s\n", etd (ent[0]), etd (ent[2]));
      rerror ("mexp() operation not supported");
    }
    rval = (*vfptr) (ent_data (ent[0]), t, ent_data (ent[2]), imet);
  }
  else if (nlastarg == 2)
  {
    //
    // mexp(A,t)
    //
    if (typ[1]==MATRIX_DENSE_REAL)
      t = class_double(ent[1]);
    else
      rerror("mexp: incorrect second argument");

    rtype = mexp1_method[typ[0]].type;
    vfptr = (VFPTR) mexp1_method[typ[0]].op;

    if (vfptr == 0)
    {
      fprintf (stderr, "Entity type: %s\n", etd (ent[0]));
      rerror ("mexp() operation not supported");
    }

    rval = (*vfptr) (ent_data (ent[0]), t, ideg);
  }
  else if (nlastarg == 1)
  {
    //
    // mexp(A)
    //
    t = 1.0;

    rtype = mexp1_method[typ[0]].type;
    vfptr = (VFPTR) mexp1_method[typ[0]].op;

    if (vfptr == 0)
    {
      fprintf (stderr, "Entity type: %s\n", etd (ent[0]));
      rerror ("mexp() operation not supported");
    }
    rval = (*vfptr) (ent_data (ent[0]), t, ideg);

  }
  else
    rerror("mexp: terrible internal error");

  for (i=0; i<nargs;i++)
  {
    if (ent[i])
      if(typ[i] != UNDEF)
        ent_Clean(ent[i]);
  }

  rent = ent_Create ();
  ent_data (rent) = rval;
  ent_type (rent) = rtype;
  return (rent);
}

/* **************************************************************
 * mpow function:
 *  mpow(A,j) := A^j
 * ************************************************************** */

Ent *
Mpow (int nargs, Datum args[])
{
  double dj;
  void *(*vfptr) ();
  Ent *e1=0, *e2=0, *rent=0;

  if (nargs != 2)
    rerror ("mpow: two arguments required");

  // get matrix
  e1 = bltin_get_ent (args[0]);
  if (ent_type(e1) != MATRIX_DENSE_REAL &&
      ent_type(e1) != MATRIX_DENSE_COMPLEX)
    rerror("mpow: first argument has to be a dense numerical matrix");

  // get exponent
  e2 = bltin_get_ent (args[1]);
  if (ent_type(e2) != MATRIX_DENSE_REAL)
    rerror("mpow: second argument 'power' has to be integer");

  dj = class_double(e2);

  if (dj != floor(dj))
    rerror("mpow: second argument 'power' has to be integer");

  //
  // mpow(A,j)
  //
  vfptr = (VFPTR) mpow1_method[ ent_type(e1) ].op;

  if (vfptr == 0)
  {
    fprintf (stderr, "Entity type: %s\n", etd (e1));
    rerror ("mpow() operation not supported");
  }

  rent = ent_Duplicate (e1);
  (*vfptr) (ent_data (rent), (int) dj, ent_data (rent));

  ent_Clean(e1);
  ent_Clean(e2);

  return (rent);
}

#undef THIS_SOLVER
#define THIS_SOLVER "isabsmono"
Ent *
isAbsMonotone (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0;
  MDR *x;
  int j,nx,rval;
  double e=0, de;

  if (nargs != 1 && nargs != 2 )
    rerror (THIS_SOLVER ": " RLAB_ERROR_ONE_OR_TWO_ARG_REQUIRED);

  e1 = bltin_get_ent (args[0]);
  if (ent_type(e1) != MATRIX_DENSE_REAL)
    rerror(THIS_SOLVER ": " RLAB_ERROR_ARG1_MDR_VECTOR);
  x = class_matrix_real(e1);
  if (!x)
    rerror(THIS_SOLVER ": " RLAB_ERROR_ARG1_MDR_VECTOR);
  if (!(MNR(x)*MNC(x)))
    rerror(THIS_SOLVER ": " RLAB_ERROR_ARG1_MDR_VECTOR);
  if (MNR(x)!=1 && MNC(x)!=1)
    rerror(THIS_SOLVER ": " RLAB_ERROR_ARG1_MDR_VECTOR);

  if (nargs == 2)
  {
    e2 = bltin_get_ent (args[1]);
    if (ent_type(e2) != MATRIX_DENSE_REAL)
      rerror(THIS_SOLVER ": " RLAB_ERROR_ARG2_MDR_SCALAR);
    de = class_double(e2);
    if (de > 0)
      e = de;
  }

  nx = MNR(x) * MNC(x);

  rval = 1;
  for (j=1; j<nx; j++)
  {
    if (mdrV0(x,j) < mdrV0(x,j-1)-e )
    {
      rval = 0;
      break;
    }
  }

  if (!rval)
  {
    rval=-1;
    for (j=1; j<nx; j++)
    {
      if (mdrV0(x,j) > mdrV0(x,j-1)+e )
      {
        rval = 0;
        break;
      }
    }
  }

  ent_Clean(e1);
  ent_Clean(e2);

  return ent_Create_Rlab_Double((double) rval);
}

#undef THIS_SOLVER
#define THIS_SOLVER "isrelmono"
Ent *
isRelMonotone (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0;
  MDR *x;
  int j,nx,rval;
  double e=0, de;

  if (nargs != 1 && nargs != 2 )
    rerror (THIS_SOLVER ": " RLAB_ERROR_ONE_OR_TWO_ARG_REQUIRED);

  e1 = bltin_get_ent (args[0]);
  if (ent_type(e1) != MATRIX_DENSE_REAL)
    rerror(THIS_SOLVER ": " RLAB_ERROR_ARG1_MDR_VECTOR);
  x = class_matrix_real(e1);
  if (!x)
    rerror(THIS_SOLVER ": " RLAB_ERROR_ARG1_MDR_VECTOR);
  if (!(MNR(x)*MNC(x)))
    rerror(THIS_SOLVER ": " RLAB_ERROR_ARG1_MDR_VECTOR);
  if (MNR(x)!=1 && MNC(x)!=1)
    rerror(THIS_SOLVER ": " RLAB_ERROR_ARG1_MDR_VECTOR);

  if (nargs == 2)
  {
    e2 = bltin_get_ent (args[1]);
    if (ent_type(e2) != MATRIX_DENSE_REAL)
      rerror(THIS_SOLVER ": " RLAB_ERROR_ARG2_MDR_SCALAR);
    de = class_double(e2);
    if (de > 0)
      e = de;
  }

  nx = MNR(x) * MNC(x);

  rval = 1;
  for (j=1; j<nx; j++)
  {
    if (mdrV0(x,j) / mdrV0(x,j-1) < 1-e )
    {
      rval = 0;
      break;
    }
  }

  if (!rval)
  {
    rval=-1;
    for (j=1; j<nx; j++)
    {
      if (mdrV0(x,j) / mdrV0(x,j-1) > 1+e )
      {
        rval = 0;
        break;
      }
    }
  }

  ent_Clean(e1);
  ent_Clean(e2);

  return ent_Create_Rlab_Double((double) rval);
}

#undef THIS_SOLVER
#define THIS_SOLVER "localabsmaxi"
Ent *
localAbsMaxIdx (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0;
  MDR *x=0, *w=0, *s=0;
  int i,j,nx;
  double e=0, de;

  if (nargs != 1 && nargs != 2 )
    rerror (THIS_SOLVER ": " RLAB_ERROR_ONE_OR_TWO_ARG_REQUIRED);

  e1 = bltin_get_ent (args[0]);
  if (ent_type(e1) != MATRIX_DENSE_REAL)
    rerror(THIS_SOLVER ": " RLAB_ERROR_ARG1_MDR_VECTOR);
  x = class_matrix_real(e1);
  if (!x)
    rerror(THIS_SOLVER ": " RLAB_ERROR_ARG1_MDR_VECTOR);
  if (!EQVECT(x))
    rerror(THIS_SOLVER ": " RLAB_ERROR_ARG1_MDR_VECTOR);

  if (nargs == 2)
  {
    e2 = bltin_get_ent (args[1]);
    if (ent_type(e2) != MATRIX_DENSE_REAL)
      rerror(THIS_SOLVER ": " RLAB_ERROR_ARG2_MDR_SCALAR);
    de = class_double(e2);
    if (de > 0)
      e = de;
  }

  nx = MNR(x) * MNC(x) - 1;

  s = mdr_Create(nx+1,1);

  i = 0;
  j = 1;
  while (j < nx)
  {
    if ((mdrV0(x,j) > mdrV0(x,j-1) + e) && (mdrV0(x,j) > mdrV0(x,j+1) + e))
    {
      MdrV0(s,i++) = j+1;
      j++;
    }
    j++;
  }
  if (mdrV0(x,nx) > mdrV0(x,nx-1) + e)
    MdrV0(s,i++) = j+1;

  if (i)
  {
    w = mdr_Create(1,i);
    for (j=0; j<i; j++)
      MdrV0(w,j) = MdrV0(s,j);
  }

  ent_Clean(e1);
  ent_Clean(e2);
  mdr_Destroy(s);

  return ent_Assign_Rlab_MDR(w);
}

#undef THIS_SOLVER
#define THIS_SOLVER "localabsmini"
Ent *
localAbsMinIdx (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0;
  MDR *x=0, *w=0;
  int i,j,nx;
  double e=0, de;

  if (nargs != 1 && nargs != 2 )
    rerror (THIS_SOLVER ": " RLAB_ERROR_ONE_OR_TWO_ARG_REQUIRED);

  e1 = bltin_get_ent (args[0]);
  if (ent_type(e1) != MATRIX_DENSE_REAL)
    rerror(THIS_SOLVER ": " RLAB_ERROR_ARG1_MDR_VECTOR);
  x = class_matrix_real(e1);
  if (!x)
    rerror(THIS_SOLVER ": " RLAB_ERROR_ARG1_MDR_VECTOR);
  if (!EQVECT(x))
    rerror(THIS_SOLVER ": " RLAB_ERROR_ARG1_MDR_VECTOR);

  if (nargs == 2)
  {
    e2 = bltin_get_ent (args[1]);
    if (ent_type(e2) != MATRIX_DENSE_REAL)
      rerror(THIS_SOLVER ": " RLAB_ERROR_ARG2_MDR_SCALAR);
    de = class_double(e2);
    if (de > 0)
      e = de;
  }

  nx = SIZE(x)- 1;

  i=0;
  for (j=1; j<nx; j++)
  {
    if ((mdrV0(x,j) < mdrV0(x,j-1) - e) && (mdrV0(x,j) < mdrV0(x,j+1) - e))
      i++;
  }
  if (mdrV0(x,j) < mdrV0(x,j-1) - e)
    i++;

  w = mdr_Create(1*(i>0),i);

  if (i)
  {
    i=0;
    for (j=1; j<nx; j++)
    {
      if ((mdrV0(x,j) < mdrV0(x,j-1) - e) && (mdrV0(x,j) < mdrV0(x,j+1) - e))
      {
        MdrV0(w,i) = j+1;
        i++;
      }
    }
    if (mdrV0(x,j) < mdrV0(x,j-1) - e)
      MdrV0(w,i) = j+1;
  }

  ent_Clean(e1);
  ent_Clean(e2);

  return ent_Assign_Rlab_MDR(w);
}

#undef THIS_SOLVER
#define THIS_SOLVER "localrelmaxi"
Ent *
localRelMaxIdx (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0;
  MDR *x=0, *w=0;
  int i,j,nx;
  double e=0, de;

  if (nargs != 1 && nargs != 2 )
    rerror (THIS_SOLVER ": " RLAB_ERROR_ONE_OR_TWO_ARG_REQUIRED);

  e1 = bltin_get_ent (args[0]);
  if (ent_type(e1) != MATRIX_DENSE_REAL)
    rerror(THIS_SOLVER ": " RLAB_ERROR_ARG1_MDR_VECTOR);
  x = class_matrix_real(e1);
  if (!x)
    rerror(THIS_SOLVER ": " RLAB_ERROR_ARG1_MDR_VECTOR);
  if (!(MNR(x)*MNC(x)))
    rerror(THIS_SOLVER ": " RLAB_ERROR_ARG1_MDR_VECTOR);
  if (MNR(x)!=1 && MNC(x)!=1)
    rerror(THIS_SOLVER ": " RLAB_ERROR_ARG1_MDR_VECTOR);

  if (nargs == 2)
  {
    e2 = bltin_get_ent (args[1]);
    if (ent_type(e2) != MATRIX_DENSE_REAL)
      rerror(THIS_SOLVER ": " RLAB_ERROR_ARG2_MDR_SCALAR);
    de = class_double(e2);
    if (de > 0)
      e = de;
  }

  nx = MNR(x) * MNC(x) - 1;

  i=0;
  for (j=1; j<nx; j++)
  {
    if ((mdrV0(x,j)/mdrV0(x,j-1) > 1+e) && (mdrV0(x,j+1)/mdrV0(x,j) < 1-e))
      i++;
  }
  if (mdrV0(x,j)/mdrV0(x,j-1) > 1+e)
    i++;

  w = mdr_Create(1*(i>0),i);

  if (i)
  {
    i=0;
    for (j=1; j<nx; j++)
    {
      if ((mdrV0(x,j)/mdrV0(x,j-1) > 1+e) && (mdrV0(x,j+1)/mdrV0(x,j) < 1-e))
      {
        MdrV0(w,i) = j+1;
        i++;
      }
    }
    if (mdrV0(x,j)/mdrV0(x,j-1) > 1+e)
      MdrV0(w,i) = j+1;
  }

  ent_Clean(e1);
  ent_Clean(e2);

  return ent_Assign_Rlab_MDR(w);
}

#undef THIS_SOLVER
#define THIS_SOLVER "localrelmini"
Ent *
localRelMinIdx (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0;
  MDR *x=0, *w=0;
  int i,j,nx;
  double e=0, de;

  if (nargs != 1 && nargs != 2 )
    rerror (THIS_SOLVER ": " RLAB_ERROR_ONE_OR_TWO_ARG_REQUIRED);

  e1 = bltin_get_ent (args[0]);
  if (ent_type(e1) != MATRIX_DENSE_REAL)
    rerror(THIS_SOLVER ": " RLAB_ERROR_ARG1_MDR_VECTOR);
  x = class_matrix_real(e1);
  if (!x)
    rerror(THIS_SOLVER ": " RLAB_ERROR_ARG1_MDR_VECTOR);
  if (!(MNR(x)*MNC(x)))
    rerror(THIS_SOLVER ": " RLAB_ERROR_ARG1_MDR_VECTOR);
  if (MNR(x)!=1 && MNC(x)!=1)
    rerror(THIS_SOLVER ": " RLAB_ERROR_ARG1_MDR_VECTOR);

  if (nargs == 2)
  {
    e2 = bltin_get_ent (args[1]);
    if (ent_type(e2) != MATRIX_DENSE_REAL)
      rerror(THIS_SOLVER ": " RLAB_ERROR_ARG2_MDR_SCALAR);
    de = class_double(e2);
    if (de > 0)
      e = de;
  }

  nx = MNR(x) * MNC(x) - 1;

  i=0;
  for (j=1; j<nx; j++)
  {
    if ((mdrV0(x,j)/mdrV0(x,j-1) < 1- e) && (mdrV0(x,j+1)/mdrV0(x,j) > 1 + e))
      i++;
  }
  if (mdrV0(x,j)/mdrV0(x,j-1) < 1 - e)
    i++;

  w = mdr_Create(1*(i>0),i);

  if (i)
  {
    i=0;
    for (j=1; j<nx; j++)
    {
      if ((mdrV0(x,j)/mdrV0(x,j-1) < 1 - e) && (mdrV0(x,j+1)/mdrV0(x,j) >  1 + e))
      {
        MdrV0(w,i) = j+1;
        i++;
      }
    }
    if (mdrV0(x,j)/mdrV0(x,j-1) < 1 - e)
      MdrV0(w,i) = j+1;
  }

  ent_Clean(e1);
  ent_Clean(e2);

  return ent_Assign_Rlab_MDR(w);
}

