/* mdrf4.h: Matrix Dense Real Functions ... */

/*  This file is a part of RLaB ("Our"-LaB) + rlabplus
   Copyright (C) 2007  Marijan Kostrun

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

#ifndef RLAB_MDRF4_H
#define RLAB_MDRF4_H

#include "rlab.h"
#include "ent.h"
#include "btree.h"
#include "mdr.h"

#include <stdio.h>
#include <math.h>

// mexp
extern MDR * mdr_mexp1 (MDR * a, double t, int ideg);
extern MDR * mdr_mdr_mexp2 (MDR * a, double t, MDR * v, int ideg);

// mpow
//extern MDR * mdr_mpow1 (MDR * a, int j);
extern int mdr_mpow1 (MDR * a, int j, MDR * b);

#endif /* RLAB_MDRF4_H */
