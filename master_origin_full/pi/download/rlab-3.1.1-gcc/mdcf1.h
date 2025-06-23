/* mdcf1.h: Matrix Dense Complex Functions */

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

#ifndef RLAB_MDCF1_H
#define RLAB_MDCF1_H

#include "rlab.h"
#include "mdc.h"
#include "mdr.h"
#include "btree.h"
#include "ent.h"
#include "complex.h"

#include <stdio.h>
#include <math.h>

extern MDC mdc_Dense(MDC *m);

extern void *mdc_MemberRef (MDC * m, char *name, int *type);
extern char **mdc_Members (MDC * m, int *n);

extern MDS *mdc_Type_BF (MDC * m);
extern MDC *mdc_Int_BF (MDC * m);
extern MDC *mdc_CeilFlooRound_BF (int is_ceil, MDC * m, MDR *b, MD *o);
extern MDC *mdc_Round_BF (MDC * m);
extern MDR *mdc_Abs (MDC * m);
extern MDR *mdc_MaxI1 (MDC * m);
extern MDC *mdc_Max2 (MDC * m1, MDC * m2);
extern MDC *mdc_Min1 (MDC * m);
extern MDR *mdc_MinI1 (MDC * m);
extern MDC *mdc_Min2 (MDC * m1, MDC * m2);

extern MDC *mdc_Sum_BF (MDC * m, void *h);
extern MDC *mdc_CumSum_BF (MDC * m);

extern MDC *mdc_Prod_BF (MDC * m);
extern MDC *mdc_CumProd_BF (MDC * m);

extern int mdc_IsSymmetric (MDC * m);
extern void mdc_Detect_Inf (MDC * m);
extern void mdc_Detect_Nan (MDC * m);

extern MDR *mdc_IsInf (MDC * m);
extern MDR *mdc_IsNan (MDC * m);
extern MDR *mdc_Finite (MDC * m);

extern MDC *mdc_Sin (MDC * m);
extern MDC *mdc_Cos (MDC * m);
extern MDC *mdc_Tan (MDC * m);
extern void *mdc_ASin (MDC * m, int *type);
extern void *mdc_ACos (MDC * m, int *type);
extern void *mdc_ATan (MDC * m, int *type);

extern void *mdc_Sqrt (MDC * m, int *type);
extern void *mdc_Log (MDC * m, int *type);
extern void *mdc_Log10 (MDC * m, int *type);
extern MDC *mdc_Exp (MDC * m);

extern MDC *mdc_Diag (MDC * m, int k);

extern MDR *mdc_Real_BF (MDC * m);
extern MDR *mdc_Imag_BF (MDC * m);
extern MDC *mdc_Conj_BF (MDC * m);

extern MDC *mdc_Find_BF (MDC * m);
extern Btree *mdc_Sort_BF (MDC * m);
extern MDC *mdc_Sign_BF (MDC * m);

extern MDC *mdc_Increment (MDC * m);
extern MDC *mdc_Decrement (MDC * m);

extern size_t mdc_Sizeof (MDC * m);

extern MDR *mdc_Any (MDC * m);
extern MDR *mdc_All (MDC * m);

extern MDC *mdc_Mod (MDC * m1, MDC * m2);
#endif /* RLAB_MDCF1_H */
