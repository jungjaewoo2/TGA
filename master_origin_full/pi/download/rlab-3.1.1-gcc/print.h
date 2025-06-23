/* print.h */

/*  This file is a part of RLaB ("Our"-LaB)
   Copyright () 1995  Ian R. Searle

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

#ifndef RLAB_PRINT_H
#define RLAB_PRINT_H

#include "rlab.h"

extern Ent *Printf (int nargs, Datum args[]);
extern Ent *FPrintf (int nargs, Datum args[]);
extern Ent *SPrintf (int nargs, Datum args[]);
extern Ent *Close (int nargs, Datum args[]);
extern Ent *Open (int nargs, Datum args[]);

#endif /* RLAB_PRINT_H */
