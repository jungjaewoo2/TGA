/* mdr_mds.h: Matrix Dense Real and Matrix Dense String Interaction. */

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

#ifndef RLAB_MDRMDS_H
#define RLAB_MDRMDS_H

#include "rlab.h"
#include "ent.h"
#include "mdr.h"
#include "mdc.h"
#include "mds.h"

extern int *mds_IndCoerceInt (MDS * m1, MDR * i);
extern MDS *mdr_coerce_mds (MDR * m1, MDS * i);

extern MDS *mdr_mds_MatrixAssign (MDR * var, int *i, int *j, MDS * rhs);
extern MDS *mdr_mds_MatrixAssignR (MDR * var, int *i, MDS * rhs);
extern MDS *mdr_mds_MatrixAssignC (MDR * var, int *j, MDS * rhs);

extern MDS *mdr_mds_VectorAssign (MDR * var, int *i, MDS * rhs);

extern MDR *mdr_mds_Append (MDR * m1, MDS * m2);
extern MDR *mds_mdr_Append (MDS * m1, MDR * m2);

extern MDR *mdr_mds_Stack (MDR * m1, MDS * m2);
extern MDR *mds_mdr_Stack (MDS * m1, MDR * m2);

extern MDR *mds_Eq (MDS * m1, MDS * m2);
extern MDR *mdr_mds_Eq (void *m1, MDS * m2);
extern MDR *mds_mdr_Eq (MDS * m1, void *m2);

extern MDR *mds_Ne (MDS * m1, MDS * m2);
extern MDR *mdr_mds_Ne (MDR * m1, MDS * m2);
extern MDR *mds_mdr_Ne (MDS * m1, MDR * m2);

extern MDR *mds_Lt (MDS * m1, MDS * m2);
extern MDR *mds_Le (MDS * m1, MDS * m2);
extern MDR *mds_Gt (MDS * m1, MDS * m2);
extern MDR *mds_Ge (MDS * m1, MDS * m2);
extern MDR *mds_And (MDS * m1, MDS * m2);
extern MDR *mds_Or (MDS * m1, MDS * m2);
extern MDR *mds_Not (MDS * m);

extern MDR *mds_Size_BF (MDS * m);
extern MDR *mds_Length_BF (MDS * m);

#endif /* RLAB_MDRMDS_H */
