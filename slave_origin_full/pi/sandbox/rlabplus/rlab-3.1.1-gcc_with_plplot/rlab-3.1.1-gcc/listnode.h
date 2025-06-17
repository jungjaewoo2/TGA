/* listnode.h */

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

#ifndef  LIST_NODE_H
#define  LIST_NODE_H

#include "config.h"
#include <stdio.h>
#include <string.h>

typedef struct _listNode ListNode;

struct _listNode
{
  char *key;             /* Key value; Uniquely identifies list_node. */
  void *ent;             /* Data. */
  void *listent;         /* Optional data members. */
  int scope;
  ListNode *next;        /* Pointer to the next ListNode in list. */
  ListNode *prev;        /* Pointer to the previous ListNode in list. */
};

extern ListNode *listNode_Create (void);
extern void listNode_Destroy (ListNode *);
extern int listNode_DestroyNodeOnly (ListNode *listNode);
extern void listNode_DestroyLoner (ListNode *listNode);

extern void listNode_SetOwned (ListNode *var);
extern ListNode *listNode_Copy (ListNode *lnode);
extern ListNode *listNode_Duplicate (ListNode *lnode);

extern int listNode_AttachAhead (ListNode *, ListNode *);
extern int listNode_AttachBehind (ListNode *, ListNode *);

extern ListNode *listNode_AttachEnt (ListNode *lnode, void *ent);
extern ListNode *listNode_AttachListEnt (ListNode *lnode, void *lent);

extern double var_Double (ListNode *lnode);
extern void *var_Data (ListNode *lnode);

extern int listNode_Detach (ListNode *);

extern void listNode_SetKey (ListNode *, char *);

#define listNode_IsNodeAttached(lnode) ((lnode->next || lnode->prev) ? 1 : 0)

#define listNode_GetNodeAhead(lnode)     ( lnode ? (lnode->next) : (0) );
#define listNode_GetNodeBehind(lnode)    ( lnode ? (lnode->prev) : (0) );

#define listNode_GetNextNode(lnode)      (lnode->next)
#define listNode_GetPrevNode(lnode)      (lnode->prev)

#define listNode_GetKey(a)               ( a ? (((ListNode *) (a))->key) : (0) )
#define var_key(a)                       ( a ? (((ListNode *) (a))->key) : (0) )

#define listNode_GetEnt(a)               ( a ? (((ListNode *) (a))->ent) : (0) )
#define var_ent(a)                       (((ListNode *) (a))->ent)

#define listNode_GetListEnt(a)           ( a ? (((ListNode *) (a))->listent) : (0) )
#define var_listent(a)                   (((ListNode *) (a))->listent)

#define var_name(a)                      ( a ? (a->key) : (0) )

#define var_scope(a)                     (((ListNode *) (a))->scope)
#define get_var_scope(a)                 ( a ? (((ListNode *) (a))->scope) : (0) )

#endif /* LIST_NODE_H */
