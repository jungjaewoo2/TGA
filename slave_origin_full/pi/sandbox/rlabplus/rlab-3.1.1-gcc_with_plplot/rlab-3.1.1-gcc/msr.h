/* msr.h: Matrix Sparse Real */

/*  This file is a part of RLaB ("Our"-LaB)
   Copyright (C) 1996  Ian R. Searle

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

#ifndef RLAB_MSR_H
#define RLAB_MSR_H

#include "rlab.h"
#include "ent.h"
#include "btree.h"

#include <stdio.h>
#include <math.h>

//#include "msrf3.h"

struct _matrix_sparse_real{
    int nr;               /* # of rows. */
    int nc;               /* # of columns. */
    int nnz;              /* # of non-zeros. */
    int order;            /* Ordered (1), or un-ordered (0). */
    int *ia;              /* Locations in JA where new rows start. */
    int *ja;              /* Column indices of non-zeros. */
    double *d;	         /* REAL data. */
    Btree *list;
    };

typedef struct _matrix_sparse_real MSR;

extern MSR *msr_Create (int nr, int nc);
extern MSR *msr_Setup (MSR *m, int nnz);
extern void *msr_Destroy (MSR * matrix);
extern MSR *msr_Copy (MSR * m_orig);
extern MSR *msr_ReSparse (MSR *m);
extern void msr_Print (MSR * matrix, FILE *fptr);

extern void msr_WriteGeneric (MSR * matrix, FILE *fn);
extern void msr_WriteGraph (MSR * matrix, FILE *fn);
extern void msr_WriteSparse (MSR * matrix, FILE *fn);
    extern void *msr_sparskit_prtmt  (MSR *m, char *ime);
    extern void *msr_sparskit_smms   (MSR *m, char *ime);
    extern void *msr_sparskit_dump   (MSR *m, char *ime);
    extern void *msr_sparskit_pspltm (MSR *m, char *ime);



extern char *msr_Class (MSR * m);
extern void *msr_MemberRef (MSR * m, char *name, int *type);
extern char **msr_Members (MSR * m, int *n);
extern int msr_Sublist (MSR *m, char *name);
extern int *msr_IndCoerceInt (MSR * m, MDR * i);
extern double msr_GetEl (MSR *m, int i, int j);

extern MSR *msr_MatrixAssign (MSR *m, int *i, int *j, MSR *rhs);
extern MSR *msr_VectorAssign (MSR *m, int *j, MSR *rhs);
extern MSR *msr_mdr_MatrixAssign (MSR *m, int *i, int *j, MDR *rhs);
extern MSR *msr_mdr_VectorAssign (MSR *m, int *j, MDR *rhs);

extern MSR *msr_Transpose (MSR *m);
extern MSR *msr_RowPartition (MSR *m, int *row);
extern void *msr_MatrixSub (MSR *var, int *i, int *j, int *type);
extern void *msr_MatrixSubR (MSR *var, int *i, int *type);
extern void *msr_MatrixSubC (MSR *var, int *j, int *type);
extern void *msr_VectorSub (MSR *var, int *i, int *type);
extern MSR *msr_Stack (MSR *m1, MSR *m2);
extern MSR *msr_Append (MSR *m1, MSR *m2);

extern MSR *msr_Negate (MSR *m);
extern void msr_Check (MSR *m);

extern MSR *msr_Eq (MSR *m1, MSR *m2);
extern MSR *mdr_msr_Eq (MDR *m1, MSR *m2);
extern MSR *msr_mdr_Eq (MSR *m1, MDR *m2);

extern MSR *msr_Ne (MSR *m1, MSR *m2);
extern MSR *mdr_msr_Ne (MDR *m1, MSR *m2);
extern MSR *msr_mdr_Ne (MSR *m1, MDR *m2);

extern MSR *msr_CreateScalar (double val);
extern MSR *msr_PatternInverse (MSR *m);
extern void msr_PatternUnion (MSR *m1, MSR *m2,
			      int **ia, int **ja, int *nnz);

#define MSPTR(m)                     (((MSR *)(m))->d)
#define MSRPTR(m)         ((double *)(((MSR *)(m))->d))
#define MSRPTR_JA(m)      ((int    *)(((MSR *)(m))->ja))
#define MSRPTR_JA0(m,i)  (((int    *)(((MSR *)(m))->ja))[(i)])
#define MSRPTR_IA(m)      ((int    *)(((MSR *)(m))->ia))
#define MSRPTR_IA0(m,i)  (((int    *)(((MSR *)(m))->ia))[(i)])

#endif /* RLAB_MSR_H */
