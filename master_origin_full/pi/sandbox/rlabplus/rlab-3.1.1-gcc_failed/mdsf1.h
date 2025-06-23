/* mdsf1.h: Matrix Dense String Functions */

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

#ifndef RLAB_MDSF1_H
#define RLAB_MDSF1_H

#include "rlab.h"
#include "ent.h"
#include "mdr.h"
#include "btree.h"

#include <stdio.h>
#include <math.h>

extern void *mds_MemberRef (MDS * m, char *name, int *type);
extern char **mds_Members (MDS * m, int *n);

extern MDS *mds_Type_BF (MDS * m);
extern MDS *mds_Sum     (MDS * m, void * h);
extern MDS *mds_VectorSet (MDS *);
extern MDS *mds_VectorUnion (MDS *,MDS *);
extern MDS *mds_VectorIntersect(MDS *, MDS *);
extern MDS *mds_VectorComplement(MDS *, MDS *);

extern int mds_IsSymmetric (MDS * m);

extern MDR *mds_Find_BF (MDS * m);

extern Btree *mds_Sort_BF (MDS * m);

extern MDR *mds_Strtod_BF (MDS * m);
extern MDR *mds_Strtol_BF (MDS * m, int base);

extern size_t mds_Sizeof (MDS * m);

#endif /* RLAB_MDSF1_H */
