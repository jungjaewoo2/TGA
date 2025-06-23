/* rjava.h */

/* This file is a part of RLaB ("Our"-LaB)
   Copyright (C) 2013  M. Kostrun

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

extern Ent *ent_jvm_CreateJVM (int nargs, Datum args[]);
extern Ent *ent_jvm_DestroyJVM (int nargs, Datum args[]);
extern Ent *ent_jvm_CallMethod (int nargs, Datum args[]);
Bltin rlab_java_bltin[] = {
  {BLTIN, "create", ent_jvm_CreateJVM},
  {BLTIN, "destroy", ent_jvm_DestroyJVM},
  {BLTIN, "call_method", ent_jvm_CallMethod},
  {0, 0, 0},
};
