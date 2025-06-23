/* bltin2.h */

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

#ifndef RLAB_BLTIN2_H
#define RLAB_BLTIN2_H

#include "rlab.h"
#include "bltin.h"

#include <stdio.h>

extern void class_bltin2_init (void);

extern Ent *Mod (int nargs, Datum args[]);
extern Ent *Int (int nargs, Datum args[]);
extern Ent *Ceil (int nargs, Datum args[]);
extern Ent *Floor (int nargs, Datum args[]);
extern Ent *Round (int nargs, Datum args[]);
extern Ent *IsSymm (int nargs, Datum args[]);
extern Ent *Abs (int nargs, Datum args[]);
extern Ent *MinMax (int nargs, Datum args[]);
extern Ent *Max (int nargs, Datum args[]);
extern Ent *Min (int nargs, Datum args[]);
extern Ent *MaxI (int nargs, Datum args[]);
extern Ent *Max2I (int nargs, Datum args[]);
extern Ent *MinI (int nargs, Datum args[]);
extern Ent *Min2I (int nargs, Datum args[]);

extern Ent *Sin (int nargs, Datum args[]);
extern Ent *Cos (int nargs, Datum args[]);
extern Ent *Tan (int nargs, Datum args[]);
extern Ent *ASin (int nargs, Datum args[]);
extern Ent *ACos (int nargs, Datum args[]);
extern Ent *ATan (int nargs, Datum args[]);
extern Ent *ATan2 (int nargs, Datum args[]);

extern Ent *Sqrt (int nargs, Datum args[]);
extern Ent *Log (int nargs, Datum args[]);
extern Ent *Log10 (int nargs, Datum args[]);
extern Ent *Exp (int nargs, Datum args[]);

extern Ent *Diag (int nargs, Datum args[]);

extern Ent *RowNorm (int nargs, Datum args[]);
extern Ent *MNorm (int nargs, Datum args[]);
extern Ent *PNorm (int nargs, Datum args[]);
extern Ent *Eig (int nargs, Datum args[]);
extern Ent *Eigs (int nargs, Datum args[]);
extern Ent *EigS (int nargs, Datum args[]);
extern Ent *EigN (int nargs, Datum args[]);
extern Ent *EigSPSYM (int nargs, Datum args[]);
extern Ent *Svd (int nargs, Datum args[]);
extern Ent *Chol (int nargs, Datum args[]);

extern Ent *Factor (int nargs, Datum args[]);
extern Ent *Backsub (int nargs, Datum args[]);
extern Ent *Solve (int nargs, Datum args[]);
extern Ent *Rcond (int nargs, Datum args[]);

extern Ent *Hess (int nargs, Datum args[]);
extern Ent *Balance (int nargs, Datum args[]);
extern Ent *QR (int nargs, Datum args[]);

extern Ent *Schur (int nargs, Datum args[]);
extern Ent *Sylv (int nargs, Datum args[]);
extern Ent *Det (int nargs, Datum args[]);

// rlabplus extensions: polynomial toolkit
extern Ent *ent_poly_init       (int nargs, Datum args[]);
extern Ent *ent_poly_eval_diff  (int nargs, Datum args[]);
extern Ent *ent_poly_eval_int   (int nargs, Datum args[]);
extern Ent *ent_poly_eval       (int nargs, Datum args[]);
// rlabplus extensions: hyperbolic functions and their inverses
extern Ent *ent_Cosh            (int nargs, Datum args[]);
extern Ent *ent_Sinh            (int nargs, Datum args[]);
extern Ent *ent_Tanh            (int nargs, Datum args[]);
extern Ent *ent_Acosh           (int nargs, Datum args[]);
extern Ent *ent_Asinh           (int nargs, Datum args[]);
extern Ent *ent_Atanh           (int nargs, Datum args[]);



#endif /* RLAB_BLTIN2_H */
