/* rlabplus_glpk.h: Gnu Programming Kit for rlabplus */

/* This file is a part of rlabplus
   Copyright (C) 2017 M. Kostrun

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

#ifndef RLABPLUS_GLPK_H
#define RLABPLUS_GLPK_H

// rlabplus extensions: glpk_simplex.c
extern Ent * ent_glpk_read_file (int nargs, Datum args[]);
extern Ent * ent_glpk_write_file (int nargs, Datum args[]);
extern Ent * ent_glpk_solve_lp  (int nargs, Datum args[]);

Bltin rlab_glpk_bltin[] = {
  // dloess
//   {BLTIN, "read",  ent_glpk_read_file},
//   {BLTIN, "write", ent_glpk_write_file},
  {BLTIN, "solve", ent_glpk_solve_lp},
  {0, 0, 0},
};

#endif