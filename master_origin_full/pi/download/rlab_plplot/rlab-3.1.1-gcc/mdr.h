/* mdr.h: Matrix Dense Real */

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

#ifndef RLAB_MDR_H
#define RLAB_MDR_H

#include "rlab.h"
#include "ent.h"
#include "btree.h"

#include <stdio.h>
#include <math.h>

int rsizeof (Rtype rt);
extern double rlab_dmax_mdr(MDR *, int);
extern double rlab_dmin_mdr(MDR *, int);
extern double rlab_mean_vector(int length, double *x, MDR *ig_idx, MDR * iuse_idx, int ig_infs);
extern double rlab_var_vector(int length, double *x, int bias, MDR *ig_idx, MDR * iuse_idx);
extern void rlab_meanstat_from_valstat_vectors(int length, double *x, double *w, MDR *ig_idx, MDR * iuse_idx,
                                        double *m1, double *v1, int do_std);

// bracketing of plus and minus operators
extern double rlab_op_zero_abs; // |res| <= zabs => res = 0
extern double rlab_op_zero_rel; // |res| <= zrel*max(arg1,arg2) => res=0
extern double rlab_op_min;  // arg1 - arg2 = max(arg1-arg2, sub_min)
extern double rlab_op_max;  // arg1 + arg2 = min(arg1+arg2, add_max)


extern void mdr_Shift (MDR * m1, int ud, int lr);
extern void md_numeric_Shift (MD * m1, int ud, int lr);
extern void md_numeric_Flip  (MD * m1, int ud, int lr);
extern MD * md_Create (int nrow, int ncol, Rtype rt);

extern int  md_Isvalid(MD *m);
extern MDR *mdr_Create (int nr, int nc);
extern MDR *mdi_Create (int nr, int nc);
extern MDR *mdr_CreateEmpty (int nr, int nc);
extern MDR *mdr_CreateScalar (double val);
extern MDR *mdi_CreateScalar (int ival);

extern void * md_Destroy(MD *m);
extern void * mdr_Destroy (MDR * matrix);
extern MDR  * mdr_Copy (MDR * m_orig);
extern void   mdr_Copy1IntoPtr2 (MDR * m, MDR **m_new);
extern MDR *  mdr_Reshape (MDR * m, int nrow, int ncol);
extern unsigned char *
    md_extend (unsigned char *data_old, int nrow_old, int ncol_old,
               int nrow_new, int ncol_new, int size_a);
extern void mdr_Extend (MDR * m, int nrow, int ncol);

extern void
    mdr_Zero (MDR * m);
extern void
    mdr_Nan (MDR * m);
extern void
    mdr_Ones (MDR * m);
extern double
    mdr_scalar (MDR * m); extern double mdi_scalar (MDR * m);
extern char *
    mdr_CharPointer (MDR * m);
extern double *
    mdr_Double (MDR * m); extern int *mdr_Integer (MDR * m);
extern MDR *
    mdr_MatrixReal (MDR * m);
extern MDR *
    mdr_MatrixDouble (MDR * m);
extern void
    mdr_Print (MDR * m, FILE *fptr);

extern MDR *mdr_vector_create (MDR * m1, MDR * m2, MDR * m3);
extern MDR *mdr_CreateVector (int d1, int d2, int d3);

extern MDR *mdr_Add (MDR * m1, MDR * m2);
extern MDR *mdr_AddTo (MDR * m1, MDR * m2);
extern MDR *mdr_SubtractFrom (MDR * m1, MDR * m2);
extern MDR *mdr_Add_weight (MDR * m1, MDR *w1, MDR * m2, MDR *w2);
extern MDR *mdr_Subtract (MDR * m1, MDR * m2);
extern MDR *mdr_Multiply (MDR * m1, MDR * m2, int *rtype);
extern MDR *mdr_Multiply_weight (MDR * m1, MDR *w1, MDR * m2, MDR *w2);

extern MDR *mdr_Copy2double_and_Sort(MDR * x, int row_col_flat);
extern MDR *mdr_FindRootsFromInterpolationTable(MDR * x1, double yoff);


// <<val;wgt>> compatible
extern MDR *mdr_ElMultiply (MDR * m1, MDR * m2);
extern MDR *mdr_ElMultiply_weight (MDR * m1, MDR *w1, MDR * m2, MDR *w2);
extern MDR *mdr_ElMultiplyBy(MDR * m1, MDR * m2);
extern MDR *mdr_ElRdivide (MDR * m1, MDR * m2, int *rtype);
extern MDR *mdr_ElRdivideBy(MDR * m1, MDR * m2);
extern MDR *mdr_ElRdivide_weight  (MDR * m1, MDR *w1, MDR * m2, MDR *w2);
extern MDR *mdr_ElLdivide (MDR * m1, MDR * m2);

extern MDR *mdr_Rdivide (MDR * m1, MDR * m2);
extern MDR *mdr_Ldivide (MDR * m1, MDR * m2);
extern void *mdr_Power (MDR * m1, MDR * m2, int *type);
extern void *mdr_ElPower (MDR * m1, MDR * m2, int *type);
extern void *mdr_ElPower_weight (MDR * m1, MDR *w1, MDR * m2, MDR *w2);

extern MDR *mdr_Negate (MDR * m1);
extern int mdr_LogicalScalar (MDR * m);
extern int mdr_Size (MDR * m);

extern MDR *mdr_Append (MDR * m1, MDR * m2);
extern MDR *mdr_Stack (MDR * m1, MDR * m2);

extern MDR *mdr_MatrixSub (MDR * var, int *i, int *j, int *type);
extern MDR *mdr_MatrixSubR (MDR * var, int *i, int *type);
extern MDR *mdr_MatrixSubC (MDR * var, int *j, int *type);
extern MDR *mdr_PartitionCol (MDR *m, int col);

extern MDR *mdr_MatrixAssign (MDR * var, int *i, int *j, MDR * rhs);
extern MDR *mdr_MatrixAssignR (MDR * var, int *i, MDR * rhs);
extern MDR *mdr_MatrixAssignC (MDR * var, int *j, MDR * rhs);

extern MDR *mdr_VectorSub (MDR * m, int *i, int *type);
extern MDR *mdr_VectorAssign (MDR * var, int *i, MDR * rhs);

extern int *mdr_IndCoerceInt (MDR * m, MDR * i);
extern MDR *mdr_ForLoopValue (MDR * m, int i);

extern char *mdr_Class (MDR * m);
extern void *mdr_MemberRef (MDR * m, char *name, int *type);
extern char **mdr_Members (MDR * m, int *n);

extern MDR *mdr_Eq (MDR * m1, MDR * m2);
extern MDR *mdr_FindVector (MDR * m1, MDR * m2);
extern MDR *mdr_Ne (MDR * m1, MDR * m2);
extern MDR *mdr_Lt (MDR * m1, MDR * m2);
extern MDR *mdr_Le (MDR * m1, MDR * m2);
extern MDR *mdr_Gt (MDR * m1, MDR * m2);
extern MDR *mdr_Ge (MDR * m1, MDR * m2);
extern MDR *mdr_And (MDR * m1, MDR * m2);
extern MDR *mdr_Or (MDR * m1, MDR * m2);
extern MDR *mdr_Not (MDR * m);

extern MDR *
    mdr_Transpose (MDR * m);
extern void
    mdr_Transpose_inplace (MDR * m);
extern int
    md_transpose_insitu(unsigned char *a, int m, int n, int size_a);
extern MDR *mdr_ReshapeCol (MDR * m);
extern int mdr_is_positive (MDR * m);
// extern MDR *mdr_Truncate (MDR * m, int nr, int nc);

extern void mdi_intfd_WriteGeneric (char *name, MDR * m);


extern int mdr_Sublist (MDR *m, char *name);

extern MDR *mdr_Inverse (MDR * A);

//
//
//
extern int mdr_vector_issorted(MDR *x);
extern int mdr_vector_isbounded_from_below(MDR *x, double lb);
extern int mdr_vector_isbounded_from_above(MDR *x, double ub);
extern int mdr_vector_isnan(MDR *x);
extern int mdr_rgemv (MDR *A, double alfa, MDR *x, double beta, MDR **y);
extern int mdr_rgemtv(MDR *A, double alfa, MDR *x, double beta, MDR **y);
extern int mdr_matinvd(MDR *A, double *d);
extern MDR *mdr_Create_SameSize (MDR *x);
extern MDR *mdi_Create_SameSize (MDR *x);



/* **************************************************************
 * Some defines to make de-referencing matrix object easier.
 * ************************************************************** */
extern double mdr0 (MDR *m, int k, int j);
extern int    mdi0 (MDR *m, int k, int j);
extern double mdr1 (MDR *m, int k, int j);
extern int    mdi1 (MDR *m, int k, int j);
extern double mdrV0(MDR *m, int k);
extern int    mdiV0(MDR *m, int k);
extern double mdrV1(MDR *m, int k);
extern int    mdiV1(MDR *m, int k);

extern double mdr0_safe (MDR *m, int k, int j);
extern int    mdi0_safe (MDR *m, int k, int j);
extern double mdr1_safe (MDR *m, int k, int j);
extern int    mdi1_safe (MDR *m, int k, int j);
extern double mdrV0_safe(MDR *m, int k);
extern int    mdiV0_safe(MDR *m, int k);
extern double mdrV1_safe(MDR *m, int k);
extern int    mdiV1_safe(MDR *m, int k);


// ------------------------------------------------------------------- //
//
// Default: column dominant matrix storage
//
// ------------------------------------------------------------------- //
// Index matrices from 0
#define Mdr0(m,k,j)       (((double *)((MDR *)(m)->d))[(j)*(((MDR *)(m))->nrow)+(k)])
#define Mdr0_safe(m,k,j)  (((double *)((MDR *)(m)->d))[(MIN(j, (((MDR *)(m))->ncol)-1))*(((MDR *)(m))->nrow)+(MIN(k, (((MDR *)(m))->ncol)-1))])
#define Mdi0(m,k,j)       (((int    *)((MDR *)(m)->d))[(j)*(((MDR *)(m))->nrow)+(k)])
#define Mdi0_safe(m,k,j)  (((int    *)((MDR *)(m)->d))[(MIN(j, (((MDR *)(m))->ncol)-1))*(((MDR *)(m))->nrow)+(MIN(k, (((MDR *)(m))->ncol)-1))])
#define MdrV0(m,k)     (((double *)(((MDR *)(m))->d))[k])
#define MdiV0(m,k)     (((int    *)(((MDR *)(m))->d))[k])
// Index matrices from 1
#define Mdr1(m,k,j)    (((double *)((MDR *)(m)->d))[(j-1)*(((MDR *)(m))->nrow)+(k-1)])
#define Mdi1(m,k,j)    (((int    *)((MDR *)(m)->d))[(j-1)*(((MDR *)(m))->nrow)+(k-1)])
#define MdrV1(m,k)     (((double *)(((MDR *)(m))->d))[(k-1)])
#define MdiV1(m,k)     (((int    *)(((MDR *)(m))->d))[(k-1)])

#define MDPTR(m)                  (((MDR *)(m))->d)
#define MDRPTR(m)      ((double *)(((MDR *)(m))->d))
#define MDIPTR(m)      ((int    *)(((MDR *)(m))->d))

// ------------------------------------------------------------------- //
//
// New: row-dominant matrix storage
//
// ------------------------------------------------------------------- //
// Index matrices from 0
#define Mdr0_rd(m,k,j)        (((double *)((MDR *)(m)->d))[(k)*(((MDR *)(m))->ncol)+(j)])
#define Mdr0_rd_safe(m,k,j)   (((double *)((MDR *)(m)->d))[(MIN(k, (((MDR *)(m))->nrow)-1))*(((MDR *)(m))->ncol)+(MIN(j, (((MDR *)(m))->nrow)-1))])
#define Mdi0_rd(m,k,j)    (((int    *)((MDR *)(m)->d))[(k)*(((MDR *)(m))->ncol)+(j)])
#define Mdi0_rd_safe(m,k,j)   (((int    *)((MDR *)(m)->d))[(MIN(k, (((MDR *)(m))->nrow)-1))*(((MDR *)(m))->ncol)+(MIN(j, (((MDR *)(m))->nrow)-1))])
#define Mdr1_rd(m,k,j)    (((double *)((MDR *)(m)->d))[(k-1)*(((MDR *)(m))->ncol)+(j-1)])
#define Mdi1_rd(m,k,j)    (((int    *)((MDR *)(m)->d))[(k-1)*(((MDR *)(m))->ncol)+(j-1)])


// ------------------------------------------------------------------- //
//
// Old
//
// ------------------------------------------------------------------- //
// Index matrices from 0
// #define Mdr0(m,k,j)    (((MDR *)(m))->d[(j)*(((MDR *)(m))->nrow)+(k)])
// #define Mdi0(m,k,j)    (((MDR *)(m))->i[(j)*(((MDR *)(m))->nrow)+(k)])
// Index matrices from 1
// #define Mdr1(m,k,j)    (((MDR *)(m))->d[(j-1)*(((MDR *)(m))->nrow)+(k-1)])
// #define Mdi1(m,k,j)    (((MDR *)(m))->i[(j-1)*(((MDR *)(m))->nrow)+(k-1)])
// Allow programmers to index into matrix like a vector
// #define MdrV0(m,k)     (((MDR *)(m))->d[k])
// #define MdrV1(m,k)     (((MDR *)(m))->d[k-1])
// #define MdiV0(m,k)     (((MDR *)(m))->i[k])
// #define MdiV1(m,k)     (((MDR *)(m))->i[k-1])
// #define MDRPTR(m)      (((MDR *)(m))->d)
// #define MDIPTR(m)      (((MDR *)(m))->i)

#define MdrNR(m)       (((MDR *)(m))->nrow)
#define MdrNC(m)       (((MDR *)(m))->ncol)

#define MDR_nrow(m)         (((MDR *)(m))->nrow)
#define MDR_ncol(m)         (((MDR *)(m))->ncol)

#endif /* RLAB_MDR_H */
