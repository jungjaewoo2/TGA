/* btree.h */

/*  This file is a part of RLaB ("Our"-LaB)
   Copyright (C) 1995  Ian R. Searle
   Copyright (C) 2005  Marijan Kostrun

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

#ifndef BTREE_H
#define BTREE_H

#include "rlab.h"
#include "listnode.h"
#include "mdr.h"
#include "mds.h"

#include <stdio.h>

extern Btree *btree_Create (void);
extern void btree_Destroy (Btree *);

extern ListNode *btree_AddNode (Btree *, ListNode *);
extern ListNode *btree_FindNode (Btree *, char *);

extern void btree_Print (Btree *bt, FILE *fptr);
extern Btree *btree_Copy (Btree *bt);
extern void btree_copy (Btree * new, ListNode * node);
extern void btree_duplicate (Btree * new, ListNode * node);
extern ListNode *btree_MemberRef (Btree *b, char *name, int *type);
extern char **btree_Members (Btree *bt, int *n);
extern char *btree_Class (Btree *bt);
extern int btree_CountNodes (Btree *bt);
extern int btree_GetRealNumNodes (Btree *root);
extern char **get_btree_node_names (ListNode *node, char **names);

extern char **
    btree_get_node_names (ListNode * node, char **names);

extern size_t btree_Sizeof (Btree *bt);
extern int btree_Sublist (Btree *bt, char *name);

extern int  btree_CountAllNodes (Btree * bt);
extern void btree_PrintAll (Btree *bt, FILE *fptr);

#define btree_NumNodes(bt)   (((Btree *)(bt))->numnodes)

#endif /* BTREE_H */
