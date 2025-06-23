/* mdcf2.h: Matrix Dense Complex Functions */

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

#ifndef RLAB_MDCF2_H
#define RLAB_MDCF2_H

#include "rlab.h"
#include "ent.h"
#include "btree.h"
#include "mdc.h"
#include "mdcf1.h"
#include "complex.h"

#include <stdio.h>
#include <math.h>

extern MDC *mdc_Solve (MDC * a, MDC * b, char *type);
extern MDC *mdc_SolveEq (MDC * a, MDC * b);
extern MDC *mdc_SolveEq_GE (MDC * a, MDC * b);
extern MDC *mdc_SolveEq_SYM (MDC * a, MDC * b);
extern double *mdc_Rcond (MDC * m);
extern MDR *mdc_RowNorm (MDC * m, char *type);
extern double *mdc_Norm (MDC * m, char *type);
extern double *mdc_PNorm (MDC * m, double p);
extern double matrix_Rcond (MDC * m);
extern void mdc_Svd (MDC * M, MDC ** rsv, MDC ** lsv,
		     MDR ** sigma, int flag);
extern MDC *mdc_LS (MDC * a, MDC * b);

extern Btree *mdc_EigS (MDC * m, int * issym);
extern Btree *mdc_EigG (MDC * ma, MDC * mb, int * issympos);

extern Btree *mdc_Factor (MDC * m, int flag);

extern MDC *mdc_Backsub_Sym (Btree * bl, MDC * b);
extern MDC *mdc_Backsub_Ge (Btree * bl, MDC * b);

extern Btree *mdc_Svd_BF (MDC * M, int flag);
extern MDC *mdc_Chol (MDC *m);
extern Btree *mdc_Hess (MDC * M);
extern Btree *mdc_Balance (MDC * M);

extern Btree *mdc_QR (MDC * m);
extern Btree *mdc_QRP (MDC * m);
extern Btree *mdc_Schur (MDC * m);
extern MDC *mdc_Sylv (MDC * a, MDC *b, MDC *c);

extern MDC *mdc_Det (MDC *m);

#endif /* RLAB_MDCF2_H */
