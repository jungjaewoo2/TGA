/* mdr_mdc.h: Matrix Dense Real and Matrix Dense Complex Interaction. */

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

#ifndef RLAB_MDRMDC_H
#define RLAB_MDRMDC_H

#include "rlab.h"
#include "ent.h"
#include "mdr.h"
#include "mdc.h"

extern MDC *md_mdc_AddTo        (MD * m1, MD * m2);
extern MDC *md_mdc_SubtractFrom (MD * m1, MD * m2);
extern MDC *md_mdc_ElMultiplyBy (MD * m1, MD * m2);
extern MDC *md_mdc_ElRdivideBy  (MD * m1, MD * m2);


extern MDC *mdr_coerce_mdc (MDR * m);
extern MDR *mdc_coerce_mdr (MDC * m);

extern MDC *mdr_mdc_Append (MDR * mr, MDC * mc);
extern MDC *mdc_mdr_Append (MDC * mc, MDR * mr);

extern MDC *mdr_mdc_Stack (MDR * mr, MDC * mc);
extern MDC *mdc_mdr_Stack (MDC * mc, MDR * mr);

extern MDC *mdr_mdc_Add (MDR * mr, MDC * mc);
extern MDC *mdc_mdr_Add (MDC * mc, MDR * mr);

extern MDC *mdr_mdc_Subtract (MDR * mr, MDC * mc);
extern MDC *mdc_mdr_Subtract (MDC * mc, MDR * mr);

extern MDC *mdr_mdc_Multiply (MDR * mr, MDC * mc, int *rtype);
extern MDC *mdc_mdr_Multiply (MDC * mc, MDR * mr, int *rtype);

extern MDC *mdr_mdc_ElMultiply (MDR * mr, MDC * mc);
extern MDC *mdc_mdr_ElMultiply (MDC * mc, MDR * mr);

extern MDC *mdr_mdc_Rdivide (MDR * mr, MDC * mc);
extern MDC *mdc_mdr_Rdivide (MDC * mc, MDR * mr);

extern MDC *mdr_mdc_ElRdivide (MDR * mr, MDC * mc, int *rtype);
extern MDC *mdc_mdr_ElRdivide (MDC * mc, MDR * mr, int *rtype);

extern MDC *mdr_mdc_Ldivide (MDR * mr, MDC * mc);
extern MDC *mdc_mdr_Ldivide (MDC * mc, MDR * mr);

extern MDC *mdr_mdc_ElLdivide (MDR * mr, MDC * mc);
extern MDC *mdc_mdr_ElLdivide (MDC * mc, MDR * mr);

extern void *mdr_mdc_ElPower (MDR * m1, MDC * m2, int *type);
extern void *mdc_mdr_ElPower (MDC * m1, MDR * m2, int *type);

extern void *mdr_mdc_Power (MDR * m1, MDC * m2, int *type);
extern void *mdc_mdr_Power (MDC * m1, MDR * m2, int *type);

extern int *mdc_IndCoerceInt (MDC * m1, MDR * i);

extern MDC *mdr_mdc_MatrixAssign (MDR * var, int *i, int *j, MDC * rhs);
extern MDC *mdr_mdc_MatrixAssignR (MDR * var, int *i, MDC * rhs);
extern MDC *mdr_mdc_MatrixAssignC (MDR * var, int *j, MDC * rhs);

extern MDC *mdc_mdr_MatrixAssign (MDC * var, int *i, int *j, MDR * rhs);
extern MDC *mdc_mdr_MatrixAssignR (MDC * var, int *i, MDR * rhs);
extern MDC *mdc_mdr_MatrixAssignC (MDC * var, int *j, MDR * rhs);

extern MDC *mdr_mdc_VectorAssign (MDR * mlhs, int *i, MDC * mrhs);
extern MDC *mdc_mdr_VectorAssign (MDC * mlhs, int *i, MDR * mrhs);

extern MDR *mdc_Eq (MDC * m1, MDC * m2);
extern MDR *mdc_Ne (MDC * m1, MDC * m2);
extern MDR *mdc_Lt (MDC * m1, MDC * m2);
extern MDR *mdc_Le (MDC * m1, MDC * m2);
extern MDR *mdc_Gt (MDC * m1, MDC * m2);
extern MDR *mdc_Ge (MDC * m1, MDC * m2);
extern MDR *mdc_And (MDC * m1, MDC * m2);
extern MDR *mdc_Or (MDC * m1, MDC * m2);
extern MDR *mdc_Not (MDC * m);

extern MDR *mdr_mdc_Eq (MDR * m1, MDC * m2);
extern MDR *mdr_mdc_Ne (MDR * m1, MDC * m2);
extern MDR *mdr_mdc_Lt (MDR * m1, MDC * m2);
extern MDR *mdr_mdc_Le (MDR * m1, MDC * m2);
extern MDR *mdr_mdc_Gt (MDR * m1, MDC * m2);
extern MDR *mdr_mdc_Ge (MDR * m1, MDC * m2);
extern MDR *mdr_mdc_And (MDR * m1, MDC * m2);
extern MDR *mdr_mdc_Or (MDR * m1, MDC * m2);

extern MDR *mdc_mdr_Eq (MDC * m1, MDR * m2);
extern MDR *mdc_mdr_Ne (MDC * m1, MDR * m2);
extern MDR *mdc_mdr_Lt (MDC * m1, MDR * m2);
extern MDR *mdc_mdr_Le (MDC * m1, MDR * m2);
extern MDR *mdc_mdr_Gt (MDC * m1, MDR * m2);
extern MDR *mdc_mdr_Ge (MDC * m1, MDR * m2);
extern MDR *mdc_mdr_And (MDC * m1, MDR * m2);
extern MDR *mdc_mdr_Or (MDC * m1, MDR * m2);

extern MDR *mdc_Size_BF (MDC * m);
extern MDR *mdc_Length_BF (MDC * m);

extern MDC *mdr_mdc_Max2 (MDR * mr, MDC * mc);
extern MDC *mdc_mdr_Max2 (MDC * mr, MDR * mc);

extern MDC *mdr_mdc_Min2 (MDR * mr, MDC * mc);
extern MDC *mdc_mdr_Min2 (MDC * mr, MDR * mc);

extern MDC *mdr_mdc_Solve (MDR * a, MDC * b, char *type);
extern MDC *mdc_mdr_Solve (MDC * a, MDR * b, char *type);

extern MDC *mdr_mdc_Mod (MDR * a, MDC * b);
extern MDC *mdc_mdr_Mod (MDC * a, MDR * b);

#endif /* RLAB_MDRMDC_H */
