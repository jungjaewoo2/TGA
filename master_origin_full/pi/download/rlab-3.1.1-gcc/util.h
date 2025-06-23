/* util.h */

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

#ifndef RLAB_UTIL_H
#define RLAB_UTIL_H

#include "rlab.h"
#include "list.h"
#include "listnode.h"

#include <stdio.h>
#include <setjmp.h>

extern void set_progname (char *value);
extern char * get_progname ();

extern void set_util_line_nos  (int val);

/* Error and warning functions */
extern void rerror (char *s);
extern void nerror (char *s);
extern void stop_script (char *s);
extern void warning_1 (char *s);

extern jmp_buf *jmp_inc_buff (void);
extern jmp_buf *jmp_dec_buff (void);
extern int inc_buff (void);
extern int dec_buff (void);
extern int get_ijmp (void);

/* Functions for signal handling */
//extern void intcatch (int);
extern void intcatch (int);
extern void intcatch_wait (int);
extern void fpecatch (int);
extern void pipecatch (int);

extern void swap (void *base, size_t size, size_t i, size_t j);
extern char *cpstr (char *string);
extern char *cpnstr (char *string, int n);
extern char *cpstr_strip (char *string);
extern char *strappend (char *s1, char *s2);
extern char *string_add (char *s1, char *s2);
extern void  string_concat (char **str1, char *str2);

extern int get_write_diary (void);
extern FILE * get_diary_file_ptr (void);

extern Ent *convert_datum_to_ent (Datum d);
extern Ent *bltin_get_ent (Datum d);

// extern Datum  call_rlab_script (char *fname, Datum *args, int nargs);
extern Ent  * get_ent_from_rlab_script (char *, Btree *);
extern Datum  call_rlab_script_ent (Ent * ent, Datum *args, int nargs);
extern int    isfuncdatum(Datum x);
extern Ent  * ent_call_rlab_script_1arg  (Ent *efunc, Ent * e1);
extern Ent  * ent_call_rlab_script_2args (Ent *efunc, Ent * e1, Ent * e2);
extern Ent  * ent_call_rlab_script_3args (Ent *efunc, Ent * e1, Ent * e2, Ent * e3);
extern Ent  * ent_call_rlab_script_4args (Ent *efunc, Ent * e1, Ent * e2, Ent * e3, Ent * e4);
extern Ent  * ent_call_rlab_script_5args (Ent *efunc, Ent * e1, Ent * e2, Ent * e3, Ent * e4, Ent * e5);


extern int *remove_duplicates (int *ia, int size, int *nsize);
extern int *remove_invalid (int *ia, int size, int A, int *nsize);

extern int *array_set (int *input, int input_size, int *output_size);
extern int *array_union (int *in1, int in_size1, int *in2, int in_size2, int *output_size);

#endif /* RLAB_UTIL_H */
