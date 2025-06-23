/* bltin1.h */

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

#ifndef RLAB_BLTIN1_H
#define RLAB_BLTIN1_H

#include "rlab.h"
#include "list.h"
#include "bltin.h"

#include <stdio.h>

extern void class_bltin1_init (void);
extern Ent * ent_FlipLR (int nargs, Datum args[]);
extern Ent * ent_FlipUD (int nargs, Datum args[]);
extern Ent * ent_ShiftU (int nargs, Datum args[]);
extern Ent * ent_ShiftD (int nargs, Datum args[]);
extern Ent * ent_ShiftL (int nargs, Datum args[]);
extern Ent * ent_ShiftR (int nargs, Datum args[]);
extern Ent * ent_Rot90 (int nargs, Datum args[]);
extern Ent * ent_Chomp (int nargs, Datum args[]);
extern Ent *Any (int nargs, Datum args[]);
extern Ent *All (int nargs, Datum args[]);
extern Ent *Group (int nargs, Datum args[]);
extern Ent *Members (int nargs, Datum args[]);
extern Ent *Size (int nargs, Datum args[]);
extern Ent *Length (int nargs, Datum args[]);
extern Ent *Type (int nargs, Datum args[]);
extern Ent *Load (int nargs, Datum args[]);
extern Ent *Debug (int nargs, Datum args[]);
extern Ent *Reshape (int nargs, Datum args[]);
extern Ent *ent_reshape_resize (int nargs, Datum args[]);
extern Ent *Zeros (int nargs, Datum args[]);
extern Ent *Ones (int nargs, Datum args[]);
extern Ent *Inf (int nargs, Datum args[]);
extern Ent *Nan (int nargs, Datum args[]);
extern Ent *IsInf (int nargs, Datum args[]);
extern Ent *IsNan (int nargs, Datum args[]);
extern Ent *Finite (int nargs, Datum args[]);
extern Ent *Clear (int nargs, Datum args[]);
extern Ent *Sum (int nargs, Datum args[]);
extern Ent *CumSum (int nargs, Datum args[]);
extern Ent *Find (int nargs, Datum args[]);
extern Ent *FindVector (int nargs, Datum args[]);
extern Ent *Sort (int nargs, Datum args[]);
extern Ent *Sort (int nargs, Datum args[]);
extern Ent *VectorSet (int nargs, Datum args[]);
extern Ent *VectorUnion (int nargs, Datum args[]);
extern Ent *VectorIntersect (int nargs, Datum args[]);
extern Ent *VectorComplement (int nargs, Datum args[]);
extern Ent *Merge (int nargs, Datum args[]);
extern Ent *Compact (int nargs, Datum args[]);
extern Ent *Strtol (int nargs, Datum args[]);
extern Ent *Sign (int nargs, Datum args[]);
extern Ent *Sizeof (int nargs, Datum args[]);
extern Ent *System (int nargs, Datum args[]);
extern Ent *Fork   (int nargs, Datum args[]);
extern Ent *Prod (int nargs, Datum args[]);
extern Ent *CumProd (int nargs, Datum args[]);
extern Ent *Frexp (int nargs, Datum args[]);
extern Ent *Ldexp (int nargs, Datum args[]);
extern Ent *Logb (int nargs, Datum args[]);

extern Ent *Tic (int nargs, Datum args[]);
extern Ent *Toc (int nargs, Datum args[]);
extern Ent *Sleep (int nargs, Datum args[]);

extern int get_fwidth (void);
extern int get_fprec (void);
extern Ent *Format (int nargs, Datum args[]);
extern Ent *Tmpnam (int nargs, Datum args[]);

extern Ent *Putenv (int nargs, Datum args[]);
extern Ent *Getenv (int nargs, Datum args[]);
extern Ent *Cd (int nargs, Datum args[]);

extern Ent *ent_copyentity (int nargs, Datum args[]);

// BITWISE
extern Ent *ent_bit_shift_left (int nargs, Datum args[]);
extern Ent *ent_bit_shift_right (int nargs, Datum args[]);
extern Ent *ent_byte_split (int nargs, Datum args[]);
extern Ent *ent_bit_split (int nargs, Datum args[]);
extern Ent *ent_bit_join (int nargs, Datum args[]);
extern Ent *ent_byte_join (int nargs, Datum args[]);

// MISC
extern Ent * ent_Contour(int nargs, Datum args[]);
extern Ent * ent_IsVec(int nargs, Datum args[]);
extern Ent * ent_IsMat(int nargs, Datum args[]);
extern Ent * ent_IsScalar(int nargs, Datum args[]);
extern Ent * ent_IsNumber(int nargs, Datum args[]);
extern Ent * ent_IsString(int nargs, Datum args[]);
extern Ent * ent_Filename (int nargs, Datum args[]);
extern Ent * ent_IsEmpty (int nargs, Datum args[]);

#endif /* RLAB_BLTIN1_H */
