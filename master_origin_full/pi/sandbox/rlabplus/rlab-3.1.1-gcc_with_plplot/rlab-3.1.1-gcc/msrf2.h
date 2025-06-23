/* msrf2.h: Matrix Sparse Real Functions ... */

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

#ifndef RLAB_MSRF2_H
#define RLAB_MSRF2_H

#include "rlab.h"
#include "ent.h"
#include "msr.h"
#include "mdr.h"
#include "btree.h"

#include <stdio.h>
#include <math.h>

//
// access to solver for sparse linear problems
//
extern MDR * msr_Solve (MSR *a, MDR *b, char *type);
extern Ent * ent_sparse_realsolv (int nargs, Datum args[]);
extern Ent * ent_sparse_prec     (int nargs, Datum args[]);
extern Ent * ent_sparse_iter     (int nargs, Datum args[]);
extern Ent * ent_sparse_tol      (int nargs, Datum args[]);
extern Ent * ent_sparse_dpt      (int nargs, Datum args[]);


// UMFPACK
    extern MDR *   umfpack_msr_Solve (MSR *a, MDR *b, char *type);
    extern MDR *   umfpack_msr_Det   (MSR *a);
    extern Btree * umfpack_msr_Factor(MSR *a, MDR *b, char *type);

// Super-LU
    extern MDR *   superlu_msr_Solve (MSR *a, MDR *b, char *type);
    extern MDR *   superlu_msr_spsolve(MSR *a, MDR *b);
    extern Btree * msr_Factor (MSR *A, int flag);
    extern MDR *   msr_Backsub (Btree *bt, MDR *B);

extern MSR *   msr_ReshapeCol (MSR *m);

extern MDR *   msr_SpOrder (MSR *m);

extern Ent *   SpSolve_BF (int nargs, Datum args[]);
extern Ent *   SpFactor_BF (int nargs, Datum args[]);

#endif /* RLAB_MSRF2_H */
