/* mdrf1.h: Matrix Dense Real Functions ... */

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

#ifndef RLAB_MDRF1_H
#define RLAB_MDRF1_H

#include "rlab.h"
#include "ent.h"
#include "mdr.h"
#include "btree.h"

#include <stdio.h>
#include <math.h>

extern MDR *mdr_Var(MDR *x, MDR *fact, MDR *use, MDR *ignore, int rowdominant, int igninfs, int unbias, MDR *m, int flat);
extern MDR *mdr_Covar(MDR *x, MDR *fact, MDR *use, MDR *ignore, int rowdominant, int igninfs, int unbias, MDR *m);
extern MDR *mdr_Mean(MDR *x, MDR *fact, MDR *use, MDR *ignore, int rowdominant, int igninfs, int flat);
extern MDR *mdr_VectorSet(MDR *);
extern MDR *mdr_VectorUnion(MDR *,MDR *);
extern MDR *mdr_VectorIntersect(MDR *,MDR *);
extern MDR *mdr_VectorComplement(MDR *,MDR *);
extern MDR *mdr_Merge(MDR *,MDR *);
extern MDR *mdr_Compact(MDR *,int,MDR *, MDR *);
extern MDR *mdr_Any (MDR * m);
extern MDR *mdr_All (MDR * m);
extern MDR *mdr_Size_BF (MDR * m);
extern MDR *mdr_Length_BF (MDR * m);
//
//extern MDR *mdr_Rand (int nrow, int ncol);
//
extern MDS *mdr_Type_BF (MDR * m);
extern MDR *mdr_Mod_BF (MDR * m1, MDR * m2);
extern MDR *mdr_Float_BF (MDR *m);
extern MDR *mdr_Int_BF (MDR * m1);
extern MDR *mdr_CeilFlooRound_BF (int is_ceil, MDR * m, MDR *b, MDR *o, MDR *c);
extern MDR *mdr_LogCeilFlooRound_BF (int is_ceil, MDR * m, double b);
extern MDR *mdr_MeshCeilFlooRound_BF (int is_ceil, MDR * m, MDR *mesh);
extern MDR *mdr_Floor_BF (MDR * m, MDR *b, MDR *o);
extern MDR *mdr_Round_BF (MDR * m);
extern MDR *mdr_Abs (MDR * m);
extern MDR *mdr_MinMax1_ValIdx (MDR *, int, int, int, int, double*, double*);
extern MDR *mdr_Max1 (MDR * m);
extern MDR *mdr_MaxI1 (MDR * m);
extern MDR *mdr_Max2 (MDR * m1, MDR * m2, int ignore_inf);
extern MDR *mdr_Min1 (MDR * m);
extern MDR *mdr_MinI1 (MDR * m);
extern MDR *mdr_Min2 (MDR * m1, MDR * m2, int ignore_inf);
extern void find_idxs_min_max(MDR * distmatrix, int* ip, int* jp, double * distance, int *dir);

extern MDR *mdr_Sum_BF (MDR * m, void *h);
extern MDR *mdr_CumSum_BF (MDR * m);

extern MDR *mdr_Prod_BF (MDR * m);
extern MDR *mdr_CumProd_BF (MDR * m);

int mdr_IsSorted (MDR * m, int row_dominant, int ifact);
extern int mdr_IsSymmetric (MDR * a);
extern void mdr_Detect_Inf (MDR * m);
extern void mdr_Detect_Nan (MDR * m);
extern MDR *mdr_Finite (MDR * m);

extern MDR *mdr_IsInf (MDR * m);
extern MDR *mdr_IsNan (MDR * m);
extern int IsFinite (MDR * m);
extern MDR *mdr_CreateNan (int nr, int nc);
extern MDR *mdr_CreateInf (int nr, int nc);

extern MDR *mdr_Sin (MDR * m);
extern MDR *mdr_Cos (MDR * m);
extern MDR *mdr_Tan (MDR * m);
extern void *mdr_ASin (MDR * m, int *type);
extern void *mdr_ACos (MDR * m, int *type);
extern void *mdr_ATan (MDR * m, int *type);
extern void *mdr_ATan2 (MDR * m1, MDR *m2, int *type);

extern void *mdr_Sqrt (MDR * m, int *type);
extern void *mdr_Log (MDR * m, int *type);
extern void *mdr_Log10 (MDR * m, int *type);
extern MDR *mdr_Exp (MDR * m);
extern MDR *mdr_Exp_weight (MDR * m, MDR *w);

extern MDR *mdr_Diag (MDR * m, int k);

extern MDR *mdr_Real_BF (MDR * m);
extern MDR *mdr_Imag_BF (MDR * m);
extern MDR *mdr_Conj_BF (MDR * m);

extern MDR *mdr_Find_BF (MDR * m);
extern Btree *mdr_Sort_BF (MDR * m);
extern MDR *mdr_Sign_BF (MDR * m);

extern MDR *mdr_Increment (MDR * m);
extern MDR *mdr_Decrement (MDR * m);

extern size_t mdr_Sizeof (MDR * m);

extern Btree *mdr_Frexp_BF (MDR * m);
extern MDR *mdr_Ldexp_BF (MDR * f, MDR *e);

extern MDR * mdr_CreateFillSind (int nrow, int ncol);

#ifdef HAVE_LOGB
extern MDR *mdr_Logb_BF (MDR * m);
#endif /* HAVE_LOGB */
extern MDR *mdr_Dense (MDR *m);

// BITWISE
extern MDR * mdi_bit_shift_left (MDR * a, MDR * k, MDR *s);
extern MDR * mdi_bit_shift_right (MDR * a, MDR * k, MDR *s);
extern MDR * mdi_byte_split (int a, int use_lsb);
extern MDR * mdi_bit_split  (MDR * b, int l, int use_lsb);
extern MDR * mdi_byte_join (MDR * a, int stride, int use_lsb);

// mis
void mdr_Contour(MDR * d, MDR *x, MDR *y, MDR *z, MDR **r);
#endif /* RLAB_MDRF1_H */


