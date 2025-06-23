/* mscf1.h: Matrix Sparse Complex Functions ... */

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

#ifndef RLAB_MSCF1_H
#define RLAB_MSCF1_H

#include "rlab.h"
#include "ent.h"
#include "msc.h"
#include "msr.h"
#include "mdr.h"
#include "btree.h"

#include <stdio.h>
#include <math.h>

extern MSC *msc_Sparse (MDC * m);
extern MDC *msc_Dense (MSC * m);
extern MSC *mdc_Spconvert (MDC *m);
extern MDC *msc_Spconvert (MSC *m);

extern MDC *msc_mdc_Append (MSC *m1, MDC *m2);
extern MDC *mdc_msc_Append (MDC *m1, MSC *m2);

extern MDC *msc_mdc_Stack (MSC *m1, MDC *m2);
extern MDC *mdc_msc_Stack (MDC *m1, MSC *m2);

extern MSC *msc_Add (MSC *m1, MSC *m2);
extern MDC *mdc_msc_Add (MDC *m1, MSC *m2);
extern MDC *msc_mdc_Add (MSC *m1, MDC *m2);
extern MDC *msr_msc_Add (MSR *l, MSC *r);
extern MDC *msc_msr_Add (MSC *l, MSR *r);

extern MSC *msc_Subtract (MSC *m1, MSC *m2);
extern MDC *mdc_msc_Subtract (MDC *m1, MSC *m2);
extern MDC *msc_mdc_Subtract (MSC *m1, MDC *m2);

extern MSC *msc_Multiply (MSC *m1, MSC *m2, int *rtype);
extern void *mdc_msc_Multiply (MDC *m1, MSC *m2, int *rtype);
extern MSC *msc_mdc_Multiply (MSC *m1, MDC *m2, int *rtype);

extern MDC *msc_ElRdivide (MSC *l, MSC *r, int *rtype);
extern MSC *msc_mdr_ElRdivide (MSC *l, MDR *r, int *rtype);

extern size_t msc_Sizeof (MSC *m);
extern MDC *msc_Size_BF (MSC *m);
extern MDS *msc_Type_BF (MSC *m);

extern MDR *msc_Any (MSC *m);
extern MDR *msc_All (MSC *m);
extern MSC *msc_Conj (MSC *m);
extern MDR *msc_Length_BF (MSC *m);
extern MSR *msc_Abs (MSC *m);

extern void msc_Detect_Inf (MSC *m);
extern void msc_Detect_NaN (MSC *m);

extern int msc_IsSymmetric (MSC *m);
extern MDR *msc_Find_BF (MSC *m);
extern MSR *msc_Real_BF (MSC *m);
extern MSR *msc_Imag_BF (MSC *m);

extern MDC *msc_Sum_BF (MSC *m, void *h);
extern MDC *msc_Norm (MSC *m, char *type);

extern MDC *msc_Max1 (MSC *m);
extern MDC *msc_Min1 (MSC *m);

extern MSC *msc_Int_BF (MSC * m);
extern MSC *msc_Ceil_BF (MSC * m, MDR *b, MDR *o);
extern MSC *msc_Floor_BF (MSC * m, MDR *b, MDR *o);
extern MSC *msc_Round_BF (MSC * m);

#endif /* RLAB_MSCF1_H */
