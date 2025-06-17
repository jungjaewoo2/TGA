/* mdrf2.h: Matrix Dense Real Functions ... */

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

#ifndef RLAB_MDRF2_H
#define RLAB_MDRF2_H

#include "rlab.h"
#include "ent.h"
#include "btree.h"
#include "mdr.h"
#include "mdrf1.h"

#include <stdio.h>
#include <math.h>

extern MDR *mdr_Solve (MDR * a, MDR * b, char *type);
extern MDR *mdr_SolveEq (MDR * a, MDR * b);
extern MDR *mdr_SolveEq_GE (MDR * a, MDR * b);
extern MDR *mdr_SolveEq_SYM (MDR * a, MDR * b);
extern double *mdr_Rcond (MDR * m);
extern double *mdr_Norm (MDR * m, char *type);
extern double *mdr_RowNorm (MDR * m, char *type);
extern double *mdr_PNorm (MDR * m, double p);
extern MDR *mdr_LS (MDR * a, MDR * b);

extern Btree *mdr_EigS (MDR * m, int * issym);
extern Btree *mdr_EigG (MDR * ma, MDR * mb, int * issympos);
extern Btree *mdr_EigG_SYMPD (MDR * ma, MDR * mb);

extern Btree *mdr_Factor (MDR * m, int flag);

extern MDR *mdr_Backsub_Sym (Btree * bl, MDR * b);
extern MDR *mdr_Backsub_Ge (Btree * bl, MDR * b);

extern Btree *mdr_Svd_BF (MDR * M, int flag);
extern MDR *mdr_Chol (MDR *m);
extern Btree *mdr_Hess (MDR * m);
extern Btree *mdr_Balance (MDR * m);
extern Btree *mdr_QR (MDR * m);
extern Btree *mdr_QRP (MDR * m);
extern Btree *mdr_Schur (MDR * m);
extern MDR *mdr_Sylv (MDR * a, MDR *b, MDR *c);
extern MDR *mdr_Det (MDR *m);

#endif /* RLAB_MDRF2_H */
