/* mscf2.h: Matrix Sparse Complex Functions ... */

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

#ifndef RLAB_MSCF2_H
#define RLAB_MSCF2_H

#include "rlab.h"
#include "ent.h"
#include "msc.h"
#include "msr.h"
#include "mdr.h"
#include "btree.h"

#include <stdio.h>
#include <math.h>

// General
extern MSC *msc_ReshapeCol (MDC *m);
extern MDR *msc_SpOrder (MSC *m);

// UMFPack specific
extern MDC * umfpack_msc_SolveRC(MSR * a, MDC * b, char *type);
extern MDC * umfpack_msc_Det(MSC *a);

// SuperLU specific
extern Btree *msc_Factor (MSC *a, int flag);
extern MDC *msc_Backsub (Btree *bt, MDC *b);
extern MDC *msc_SpSolve (MSC * a, MDC * b, double diag_pivot_thresh, int *perm_c);
extern MDC *msc_Solve (MSC * a, MDC * b, char *type);
extern Btree *msc_SpFactor (MSC * a);

// SolveCC
extern Ent * ent_sparse_compsolv(int nargs, Datum args[]);
extern MDC * msc_SolveCC(MSC * a, MDC * b, char *type);
extern MDC * umfpack_msc_SolveCC(MSC * a, MDC * b, char *type);// UMFPack
extern MDC * superlu_msc_SolveCC(MSC * a, MDC * b, char *type);// SuperLU

// SolveCR
extern MDC * msc_SolveCR(MSC * a, MDC * b, char *type);
extern MDC * umfpack_msc_SolveCR(MSC * a, MDR * b, char *type);// UMFPack
extern MDC * superlu_msc_SolveCR(MSC * a, MDC * b, char *type);// SuperLU



#endif /* RLAB_MSCF2_H */
