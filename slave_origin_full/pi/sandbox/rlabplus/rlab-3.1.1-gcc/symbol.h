/* symbol.h */

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

#ifndef RLAB_SYMBOL
#define RLAB_SYMBOL

#include "rlab.h"
#include "list.h"
#include "btree.h"
#include <stdio.h>

extern void symbol_table_create (void);
extern Btree *get_symtab_ptr (void);
extern Ent *get_symtab_ent (void);

extern ListNode *install (void *sym_table, char *name, Ent *ent);
extern ListNode *lookup (List * sym_table, char *s);
extern ListNode *lookup_private (List * sym_table, char *s);
extern void reset_global_symbol_table_ptr (void);

extern Var *gst (void);
extern Var *null (void);

#endif /* RLAB_SYMBOL */
