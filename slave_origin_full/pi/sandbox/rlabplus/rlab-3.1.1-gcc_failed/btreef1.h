/* btreef1.h */

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

#ifndef BTREEF1_H
#define BTREEF1_H

#include "rlab.h"
#include "btree.h"
#include "listnode.h"
#include <stdio.h>

extern MDR *btree_Size_BF (Btree *bt);
extern MDS *btree_Type_BF (Btree *bt);

extern Ent *btree_protect(int nargs, Datum args[]);
extern Ent *btree_release(int nargs, Datum args[]);
extern Ent *btree_IsConst(int nargs, Datum args[]);

#endif /* BTREE_H */
