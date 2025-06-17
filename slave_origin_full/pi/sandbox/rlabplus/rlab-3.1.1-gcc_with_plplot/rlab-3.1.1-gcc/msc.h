/* msc.h: Matrix Sparse Complex */

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

#ifndef RLAB_MSC_H
#define RLAB_MSC_H

#include "rlab.h"
#include "ent.h"
#include "btree.h"
#include "mdc.h"
#include "complex.h"

#include <stdio.h>
#include <math.h>

struct _matrix_sparse_complex
{
  int nr;               /* # of rows. */
  int nc;               /* # of columns. */
  int nnz;              /* # of non-zeros. */
  int order;            /* Ordered (1), or un-ordered (0). */
  int *ia;              /* Locations in JA where new rows start. */
  int *ja;              /* Column indices of non-zeros. */
  Complex *c;	        /* Complex data. */
  Btree *list;
};

typedef struct _matrix_sparse_complex MSC;

extern MSC *msc_Create (int nr, int nc);
extern MSC *msc_Setup (MSC *m, int nnz);
extern void *msc_Destroy (MSC * matrix);
extern MSC *msc_Copy (MSC * m_orig);
extern MSC *msc_ReSparse (MSC *m);
extern void msc_Print (MSC * matrix, FILE *fptr);
extern void msc_WriteGeneric (MSC * matrix, FILE *fn);
extern void msc_WriteGraph (MSC *m, FILE *fn);
extern void msc_WriteSparse (MSC *m, FILE *fn);

extern char *msc_Class (MSC * m);
extern void *msc_MemberRef (MSC * m, char *name, int *type);
extern char **msc_Members (MSC * m, int *n);
extern int msc_Sublist (MSC *m, char *name);
extern int *msc_IndCoerceInt (MSC * m, MDR * i);
extern Complex msc_GetEl (MSC *m, int i, int j);

extern MSC *msc_MatrixAssign (MSC *m, int *i, int *j, MSC *rhs);
extern MSC *msc_VectorAssign (MSC *m, int *j, MSC *rhs);
extern MSC *msc_mdc_MatrixAssign (MSC *m, int *i, int *j, MDC *rhs);
extern MSC *msc_msr_MatrixAssign (MSC *mlhs, int *irow, int *jcol, MSR *rhs);
extern MSC *msc_mdr_MatrixAssign (MSC *mlhs, int *irow, int *jcol, MDR *rhs);

extern MSC *msc_mdc_VectorAssign (MSC *m, int *j, MDC *rhs);
extern MSC *msc_msr_VectorAssign (MSC *mlhs, int *j, MSR *rhs);
extern MSC *msc_mdr_VectorAssign (MSC *mlhs, int *j, MDR *rhs);

extern MSC *msc_Transpose (MSC *m);
extern MSC *msc_NcTranspose (MSC *m);
extern MSC *msc_RowPartition (MSC *m, int *row);
extern void *msc_MatrixSub (MSC *var, int *i, int *j, int *type);
extern void *msc_MatrixSubR (MSC *var, int *i, int *type);
extern void *msc_MatrixSubC (MSC *var, int *j, int *type);
extern void *msc_VectorSub (MSC *var, int *i, int *type);
extern MSC *msc_Stack (MSC *m1, MSC *m2);
extern MSC *msc_Append (MSC *m1, MSC *m2);

extern MSC *msc_Negate (MSC *m);
extern MSR *msc_Eq (MSC * m1, MSC * m2);
extern MSR *mdc_msc_Eq (MDC * m1, MSC * m2);
extern MSR *msc_mdc_Eq (MSC * m1, MDC * m2);

extern MSR *msc_mdr_Eq (MSC * m1, MDR * m2);
extern MSR *mdr_msc_Eq (MDR * m1, MSC * m2);
extern MSR *msc_msr_Eq (MSC * m1, MSR * m2);
extern MSR *msr_msc_Eq (MSR * m1, MSC * m2);

extern MSR *msc_Ne (MSC * m1, MSC * m2);
extern MSR *mdc_msc_Ne (MDC * m1, MSC * m2);
extern MSR *msc_mdc_Ne (MSC * m1, MDC * m2);

extern MSC *msc_CreateScalar (double valr, double vali);
extern MSR *msc_PatternInverse (MSC *m);
extern void msc_PatternUnion (MSC *m1, MSC *m2,
			      int **IA, int **JA, int *nnz);

extern MSR *msc_coerce_msr (MSC *m);
extern MSC *msr_coerce_msc (MSR *m);

extern void msc_Check (MSC *m);
#endif /* RLAB_MSC_H */
