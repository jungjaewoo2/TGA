/* mde.h: Matrix Dense Entity */

/*  This file is a part of rlabplus ("Our"-LaB)
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

#ifndef RLABPLUS_MDE_H
#define RLABPLUS_MDE_H

#include "rlab.h"
#include "ent.h"
#include "btree.h"

#include <stdio.h>
#include <math.h>

typedef struct _matrix_dense MDE;

extern MDE  * mde_Create (int nr, int nc);
extern MDE  * mde_CreateEmpty (int nr, int nc);
extern MDE  * mde_CreateDouble (double val);
extern MDE  * mde_CreateInt32 (int val);
extern void * mde_Destroy (MDE * matrix);
extern void * mde_DestroyEmpty (MDE * matrix);
extern void * mde_MemberRef (MDE * m, char *name, int *type);
extern MDE  * mde_Copy (MDE * m_orig);
extern MDE  * mde_Duplicate (MDE * m_orig);

extern int * mde_IndCoerceInt (MDE * m, MDR * i);
extern MDE * mde_VectorSub (MDE * m, int *i, int *type);
extern int   mde_Size(MDE * var);
extern char ** mde_Members (MDE * m, int *n);
extern char * mde_Class (MDE * m);
extern MDR  * mde_Length_BF (MDE * m);
extern MDR  * mde_Size_BF (MDE * m);


// ------------------------------------------------------------------- //
//
// New:
//
// ------------------------------------------------------------------- //
// Index entity matrices from 0
#define Mde0(m,k,j)    (((Ent **)((MDE *)(m)->d))[(j)*(((MDE *)(m))->nrow)+(k)])
#define MdeV0(m,k)     (((Ent **)(((MDE *)(m))->d))[k])
// Index entity matrices from 1
#define Mde1(m,k,j)    (((Ent **)((MDE *)(m)->d))[(j-1)*(((MDE *)(m))->nrow)+(k-1)])
#define MdeV1(m,k)     (((Ent **)(((MDE *)(m))->d))[(k-1)])
#define MDEPTR(m)      (((MDE *)(m))->d)

#endif /* RLABPLUS_MDE_H */
