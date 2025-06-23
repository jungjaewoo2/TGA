/* msrf1.h: Matrix Sparse Real Functions ... */

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

#ifndef RLAB_MSRF1_H
#define RLAB_MSRF1_H

#include "rlab.h"
#include "ent.h"
#include "msr.h"
#include "mdr.h"
#include "mdc.h"
#include "btree.h"

#include <stdio.h>
#include <math.h>

extern MSR *msr_Sparse (MDR * m);
extern MDR *msr_Dense (MSR * m);
extern MSR *mdr_Spconvert (MDR *m, int *rtype);
extern MDR *msr_Spconvert (MSR *m, int *rtype);

extern MDR *msr_mdr_Append (MSR *m1, MDR *m2);
extern MDR *mdr_msr_Append (MDR *m1, MSR *m2);

extern MDR *msr_mdr_Stack (MSR *m1, MDR *m2);
extern MDR *mdr_msr_Stack (MDR *m1, MSR *m2);

extern MSR *msr_Add (MSR *m1, MSR *m2);
extern MDR *mdr_msr_Add (MDR *m1, MSR *m2);
extern MDR *msr_mdr_Add (MSR *m1, MDR *m2);

extern MSR *msr_Subtract (MSR *m1, MSR *m2);
extern MDR *mdr_msr_Subtract (MDR *m1, MSR *m2);
extern MDR *msr_mdr_Subtract (MSR *m1, MDR *m2);

extern MSR *msr_Multiply (MSR *m1, MSR *m2, int *rtype);
extern void *mdr_msr_Multiply (MDR *m1, MSR *m2, int *rtype);
extern void *msr_mdr_Multiply (MSR *m1, MDR *m2, int *rtype);
extern void *msr_mdc_Multiply (MSR *m1, MDC *m2, int *rtype);

extern MSR *msr_ElMultiply (MSR *m1, MSR *m2, int *rtype);
extern MSR *mdr_msr_ElMultiply (MDR *m1, MSR *m2, int *rtype);
extern MSR *msr_mdr_ElMultiply (MSR *m1, MDR *m2, int *rtype);

extern MDR *msr_Rdivide (MSR *m1, MSR *m2);
extern MDR *mdr_msr_Rdivide (MDR *m1, MSR *m2);
extern MDR *msr_mdr_Rdivide (MSR *m1, MDR *m2);

extern MDR *msr_ElRdivide (MSR *m1, MSR *m2, int *rtype);
extern MDR *mdr_msr_ElRdivide (MDR *m1, MSR *m2, int *rtype);
extern void *msr_mdr_ElRdivide (MSR *m1, MDR *m2, int *rtype);

extern MDR *msr_Ldivide (MSR *a, MDR *b);
extern MDR *msr_msr_Ldivide (MSR *a, MDR *b);

extern size_t msr_Sizeof (MSR *m);
extern MDR *msr_Size_BF (MSR *m);
extern MDS *msr_Type_BF (MSR *m);

extern MDR *msr_Any (MSR *m);
extern MDR *msr_All (MSR *m);

extern MSR *msr_Conj (MSR *m);
extern MDR *msr_Length_BF (MSR *m);
extern MSR *msr_Abs (MSR *m);

extern void msr_Detect_Inf (MSR *m);
extern void msr_Detect_NaN (MSR *m);

extern MSR *msr_Diag (MSR *marg, int k);

extern int msr_IsSymmetric (MSR *m);

extern MDR * msr_Find_BF (MSR *m);
extern MSR * msr_Real_BF (MSR *m);
extern MSR * msr_Imag_BF (MSR *m);

extern MDR *msr_Sum_BF (MSR *m, void *h);
extern double msr_Norm (MSR *m, char *type);

extern MDR *msr_Max1 (MSR *m);
extern MDR *msr_Min1 (MSR *m);

extern MSR *msr_Int_BF (MSR * m);
extern MSR *msr_Ceil_BF (MSR * m);
extern MSR *msr_Floor_BF (MSR * m);
extern MSR *msr_Round_BF (MSR * m);


#endif /* RLAB_MSRF1_H */
