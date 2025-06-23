/* mds.h: Matrix Dense String */

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

#ifndef RLAB_MDS_H
#define RLAB_MDS_H

#include "rlab.h"
#include "ent.h"
#include "mdr.h"
#include "btree.h"

#include <stdio.h>
#include <math.h>

// Char array for holding tmp strings during scanning
// and other string operations
#include "buffer.h"

extern void   swap2strings      (char **, char **);
extern char * add2strings       (char  *, char  *);
extern void   add2strings2first (char **, char  *);
extern void   add2strings2last  (char  *, char **);
extern char * add3strings       (char  *, char  *, char  *);
extern void   add3strings2first (char **, char  *, char  *);
extern void   add3strings2last  (char  *, char  *, char **);
extern int    isvalidstring     (char *);
extern char * mdsV0_safe(MDS *m, int k);
extern char * mds0_safe(MDS *m, int k, int j);
extern char * mdsV1_safe(MDS *m, int k);
extern char * mds1_safe(MDS *m, int k, int j);

extern void mds_Shift(MDS *m, int ud, int lr);
extern void mds_Flip(MDS *m, int ud, int lr);
extern MDS *mds_Create (int nr, int nc);
extern MDS *mds_CreateScalar (char *string);
extern MDS *mds_AssignScalar (char *string);

extern void *mds_Destroy (MDS * matrix);
extern MDS *mds_Copy (MDS * m_orig);
extern MDS *mds_Reshape (MDS * m, int nrow, int ncol);
extern void mds_Extend (MDS * m, int nrow, int ncol);

extern MDS *mds_MatrixString (MDS *m);

extern int mds_Size (MDS * m);
extern void mds_Print (MDS * m, FILE *fptr);

extern MDS *mds_Append (MDS * m1, MDS * m2);
extern MDS *mds_Stack (MDS * m1, MDS * m2);

extern MDS *mds_MatrixSub (MDS * var, int *i, int *j, int *type);
extern MDS *mds_MatrixSubR (MDS * var, int *i, int *type);
extern MDS *mds_MatrixSubC (MDS * var, int *j, int *type);

extern MDS *mds_MatrixAssign (MDS * var, int *i, int *j, MDS * rhs);
extern MDS *mds_MatrixAssignR (MDS * var, int *i, MDS * rhs);
extern MDS *mds_MatrixAssignC (MDS * var, int *j, MDS * rhs);

extern MDS *mds_VectorSub (MDS * m, int *i, int *type);
extern MDS *mds_VectorAssign (MDS * var, int *j, MDS * rhs);

extern MDS *mds_ForLoopValue (MDS * m, int i);

extern char *mds_CharPointer (MDS * m);

extern char *mds_GetString (MDS * m);

extern MDS *mds_MakeCharMatrix (char **array, int n);

extern char *mds_Class (MDS * m);

extern MDS *mds_Transpose (MDS * m);
extern void *mds_Transpose_inplace (MDS * m);

extern MDS *mds_ReshapeCol (MDS * m);

extern void mds_Truncate_Inplace (MDS *old, int new_nrow, int new_ncol);

extern MDS *mds_Add (MDS * m1, MDS * m2);

extern void mds_intfd_WriteGeneric (char *name, MDS * m);
extern int mds_Sublist (MDS *m, char *name);

/* **************************************************************
 * Some defines to make de-referencing matrix object easier.
 * ************************************************************** */
#define MDSPTR(m)      (((MDS *)(m))->d)
#define Mds1(m,k,j)    (((char **)((MDS *)(m)->d))[(j-1)*(((MDS *)(m))->nrow)+(k-1)])
#define Mds0(m,k,j)    (((char **)((MDS *)(m)->d))[(j)*(((MDS *)(m))->nrow)+(k)])
#define MdsV0(m,k)     (((char **)(((MDS *)(m))->d))[(k  )])
#define MdsV1(m,k)     (((char **)(((MDS *)(m))->d))[(k-1)])
/* for pesky string matrices that are in fact scalars */
#define mdstring(m)    (((char **)(((MDS *)(m))->d))[0])

/* Index matrices from 1 */
// #define Mds1(m,i,j)    (((MDS *)(m))->s[(j-1)*(((MDS *)(m))->nrow)+(i-1)])
/* Index matrices from 0 */
// #define Mds0(m,i,j)    (((MDS *)(m))->s[(j)*(((MDS *)(m))->nrow)+(i)])

#define MdsNR(m)       (((MDS *)(m))->nrow)
#define MdsNC(m)       (((MDS *)(m))->ncol)

#define MDS_nrow(m)         (((MDS *)(m))->nrow)
#define MDS_ncol(m)         (((MDS *)(m))->ncol)

#endif /* RLAB_MDR_H */
