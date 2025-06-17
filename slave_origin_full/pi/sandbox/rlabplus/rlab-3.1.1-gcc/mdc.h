/* mdc.h: Matrix Dense Complex */

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

#ifndef RLAB_MDC_H
#define RLAB_MDC_H

#include "rlab.h"
#include "ent.h"
#include "complex.h"
#include "btree.h"

#include <stdio.h>
#include <math.h>

extern MDC *mdc_Create (int nr, int nc);
extern void mdc_Shift (MDC *m, int ud, int lr);
extern MDC *mdc_CreateEmpty (int nr, int nc);
extern MDC *mdc_Empty (void);

extern void *mdc_Destroy (MDC * matrix);
extern void *mdc_DestroyEmpty (MDC * matrix);
extern MDC *mdc_Copy (MDC * m_orig);
extern MDC *mdc_Reshape (MDC * m, int nrow, int ncol);
extern MDC *mdc_Extend (MDC * m, int nrow, int ncol);
extern void mdc_Truncate (MDC * m, int nrow, int ncol);

extern void mdc_Zero (MDC * m);
extern double mdc_scalar (MDC * m);
extern double *mdc_Double (MDC * m);
extern MDR *mdc_MatrixReal (MDC * m);
extern void mdc_Print (MDC * m, FILE *fptr);

extern MDC *mdc_Add (MDC * m1, MDC * m2);
extern MDC *mdc_Subtract (MDC * m1, MDC * m2);
extern MDC *mdc_Multiply (MDC * m1, MDC * m2, int *rtype);
extern MDC *mdc_ElMultiply (MDC * ma, MDC * mb);
extern MDC *mdc_Rdivide (MDC * m1, MDC * m2);
extern MDC *mdc_ElRdivide (MDC * m1, MDC * m2, int *rtype);
extern MDC *mdc_Ldivide (MDC * m1, MDC * m2);
extern MDC *mdc_ElLdivide (MDC * m1, MDC * m2);
extern void *mdc_Power (MDC * m1, MDC * m2, int *type);
extern void *mdc_ElPower (MDC * m1, MDC * m2, int *type);

extern MDC *mdc_Negate (MDC * m1);
extern int mdc_LogicalScalar (MDC * m);
extern int mdc_Size (MDC * m);

extern MDC *mdc_Append (MDC * m1, MDC * m2);
extern MDC *mdc_Stack (MDC * m1, MDC * m2);

extern MDC *mdc_MatrixSub (MDC * var, int *i, int *j, int *type);
extern MDC *mdc_MatrixSubR (MDC * var, int *i, int *type);
extern MDC *mdc_MatrixSubC (MDC * var, int *j, int *type);

extern MDC *mdc_MatrixAssign (MDC * var, int *i, int *j, MDC * rhs);
extern MDC *mdc_MatrixAssignR (MDC * var, int *i, MDC * rhs);
extern MDC *mdc_MatrixAssignC (MDC * var, int *j, MDC * rhs);

extern MDC *mdc_VectorSub (MDC * m, int *i, int *type);
extern MDC *mdc_VectorAssign (MDC * var, int *j, MDC * rhs);

extern MDC *mdc_ForLoopValue (MDC * m, int i);
extern MDC *mdc_CreateScalar (double rval, double ival);

extern char *mdc_Class (MDC * m);

extern MDC *mdc_Transpose (MDC * m);
extern MDC *mdc_NcTranspose (MDC * m);
extern void mdc_NcTranspose_inplace (MDC * z);
extern MDC *mdc_ReshapeCol (MDC * m);

extern void mdc_WriteGeneric (MDC * m, FILE * fn);
extern int mdc_Sublist (MDC *m, char *name);

extern double mdr0 (MDR *m, int k, int j);
extern int    mdi0 (MDR *m, int k, int j);
extern double mdr1 (MDR *m, int k, int j);
extern int    mdi1 (MDR *m, int k, int j);
extern double mdrV0(MDR *m, int k);
extern int    mdiV0(MDR *m, int k);
extern double mdrV1(MDR *m, int k);
extern int    mdiV1(MDR *m, int k);

extern Complex mdc0(MD *m, int k, int j);
extern Complex mdc1(MD *m, int k, int j);
extern Complex mdcV0(MD *m, int k);
extern Complex mdcV1(MD *m, int k);


/* **************************************************************
 * Some defines to make de-referencing matrix object easier (rlabplus)
 * ************************************************************** */

// accessing real and imaginery parts directly
// valid only for our 'Complex'
#define  RE(x)    (((double *) &(x))[0])
#define  IM(x)    (((double *) &(x))[1])

// Index matrices from 0
#define Mdc0(m,k,j)     (((Complex *)((MDC *)(m)->d))[(j)*(((MDC *)(m))->nrow)+(k)])
#define Mdc0r(m,k,j)    (RE(((Complex *)((MDC *)(m)->d))[(j)*(((MDC *)(m))->nrow)+(k)]))
#define Mdc0i(m,k,j)    (IM(((Complex *)((MDC *)(m)->d))[(j)*(((MDC *)(m))->nrow)+(k)]))

// Index matrices from 0
#define Mdc1(m,k,j)    (((Complex *)((MDC *)(m)->d))[(j-1)*(((MDC *)(m))->nrow)+(k-1)])
#define Mdc1r(m,k,j)   (RE(((Complex *)((MDC *)(m)->d))[(j-1)*(((MDC *)(m))->nrow)+(k-1)]))
#define Mdc1i(m,k,j)   (IM(((Complex *)((MDC *)(m)->d))[(j-1)*(((MDC *)(m))->nrow)+(k-1)]))

// Allow programmers to index into matrix like a vector starting from 0
#define MdcV0(m,k)     (((Complex *)((MDC *)(m)->d))[(k)])
#define MdcV0r(m,k)    (RE(((Complex *)((MDC *)(m)->d))[(k)]))
#define MdcV0i(m,k)    (IM(((Complex *)((MDC *)(m)->d))[(k)]))

// Allow programmers to index into matrix like a vector starting from 1
#define MdcV1(m,k)     (((Complex *)((MDC *)(m)->d))[(k-1)])
#define MdcV1r(m,k)    (RE(((Complex *)((MDC *)(m)->d))[(k-1)]))
#define MdcV1i(m,k)    (IM(((Complex *)((MDC *)(m)->d))[(k-1)]))


#define MDCPTR(m)      (((MDC *)(m))->d)
#define MdcNR(m)       (((MDC *)(m))->nrow)
#define MdcNC(m)       (((MDC *)(m))->ncol)

#define MDC_nrow(m)         (((MDC *)(m))->nrow)
#define MDC_ncol(m)         (((MDC *)(m))->ncol)

#endif /* RLAB_MDC_H */
