/* bltin.h */

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

#ifndef RLAB_BLTIN
#define RLAB_BLTIN

#include "rlab.h"
#include "list.h"
#include "mds.h"

#include <stdio.h>

#ifdef THINK_C
/* conflict with Mac library */
#define Size Sizex
#endif

struct _bltin
{
  int type;
  char *name;
  Ent *(*func) ();
};
typedef struct _bltin Bltin;

struct _doubleconst
{
    char *name;
    double val;
};
typedef struct _doubleconst Doubleconst;


extern Ent *Error (int nargs, Datum args[]);
extern Ent *Stop (int nargs, Datum args[]);
extern Ent *Exist (int nargs, Datum args[]);
extern Ent *EntInfo (int nargs, Datum args[]);
extern Ent *FixSymTable (int nargs, Datum args[]);

extern int get_fwidth (void);
extern int get_fprec (void);

extern void bltin_Print (Bltin *dont_use, FILE *fptr);

extern char *bltin_Class (Bltin * b);

extern MDS *bltin_Type_BF (Bltin * b);

extern size_t bltin_Sizeof (Bltin * b);

extern Ent *bltin_MemberRef (Bltin *bf, char *name, int *type);
extern char **bltin_Members (Bltin *bf, int *n);

extern Ent *Diary (int nargs, Datum args[]);
extern int get_write_diary (void);
extern FILE *get_diary_file_ptr (void);

#define bltin_GetFuncPtr(a)    (a)->func;

#endif /* RLAB_BLTIN */
