/* bltin3.h: More built-in functions. */

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

#ifndef RLAB_BLTIN3_H
#define RLAB_BLTIN3_H

#include "rlab.h"
#include "ent.h"
#include "msr.h"
#include "mdr.h"
#include "btree.h"

#include <stdio.h>
#include <math.h>

extern void class_bltin3_init (void);

extern Ent *Sparse (int nargs, Datum args[]);
extern Ent *Dense (int nargs, Datum args[]);
extern Ent *Spconvert (int nargs, Datum args[]);
extern Ent *SpOrder (int nargs, Datum args[]);
extern Ent *SpWrite (int nargs, Datum args[]);

extern Ent *ReadGraph (int nargs, Datum args[]);

#endif /* RLAB_BLTIN3_H */
