/* bltin4.h */

/*  This file is a part of RLaB ("Our"-LaB) + rlabplus
   Copyright (C) 2007-2016  Marijan Kostrun

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

#ifndef RLAB_BLTIN4_H
#define RLAB_BLTIN4_H
extern void class_bltin4_init (void);

extern Ent *Mexp (int nargs, Datum args[]);
extern Ent *Mpow (int nargs, Datum args[]);
extern Ent *isAbsMonotone (int nargs, Datum args[]);
extern Ent *isRelMonotone (int nargs, Datum args[]);
extern Ent *localExtremaeIdx (int nargs, Datum args[]);
extern Ent *ent_blas_givens (int nargs, Datum args[]);

#endif /* RLAB_BLTIN4_H */
