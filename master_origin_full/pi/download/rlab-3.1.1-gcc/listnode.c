/* listnode.c */

/*  This file is a part of RLaB ("Our"-LaB)
   Copyright (C) 1992  Ian R. Searle

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

#include "rlab.h"
#include "ent.h"
#include "mem.h"
#include "listnode.h"
#include "util.h"

#include <stdio.h>

/* ************************************************************************
 * PURPOSE: To dynamically allocate and initialize a ListNode.
 * ************************************************************************ */

ListNode *
listNode_Create (void)
{
  ListNode *new = (ListNode *) GC_MALLOC (sizeof (ListNode));
  if (new == 0)
    rerror ("out of memory");

  new->key = 0;

  new->ent = 0;
  new->listent = 0;
  new->scope = 0;
  new->next = 0;
  new->prev = 0;

  return (new);
}

void listNode_Destroy (ListNode * listn)
{
  if (listn->key)
    GC_FREE (listn->key);
  GC_FREE (listn);
}

/* ************************************************************************
 * PURPOSE: Free all of dynamic memory ascociated with a list-node,
 *          but NOT the attached data.
 * ************************************************************************ */

int listNode_DestroyNodeOnly (ListNode * listn)
{
  if (listn->key)
    GC_FREE (listn->key);
  GC_FREE (listn);
  return (1);
}

/* **************************************************************
 * Destroy listnodes that are loners (don't belong to anybody)
 * ************************************************************** */

void
listNode_DestroyLoner (ListNode * var)
{
  if (var->scope == 0)
  {
    if (var_listent (var))
      ent_Destroy (var_listent (var));

    listNode_DestroyNodeOnly (var);
  }
}

void
listNode_SetOwned (ListNode * var)
{
  var->scope = 1;
}

/* **************************************************************
 * Copy a ListNode.
 * ************************************************************** */

ListNode * listNode_Copy (ListNode * lnode)
{
  ListNode *new = listNode_Create ();

  new->key = cpstr (lnode->key);
  new->ent = lnode->ent;
  new->listent = lnode->listent;
  new->next = 0;
  new->prev = 0;
  return (new);
}

ListNode * listNode_Duplicate (ListNode * lnode)
{
  ListNode *new = listNode_Create ();

  new->key = cpstr (lnode->key);
  new->ent = ent_Duplicate(lnode->ent);
  new->listent = lnode->listent;
  new->next = 0;
  new->prev = 0;
  return (new);
}

/* **************************************************************
 * Attach a ListNode ahead of another ListNode.
 * listNode is put ahead of nodeBehind.  If nodeBehind already
 * had a ListNode ahead of it, then we hook up pointers
 * as necessary.
 *
 * Note that a ListNode has no concept as to what a list is.
 * All it knows how to do is to attach and detach itself from
 * other ListNode's.  Lists will use ListNode's to build
 * and destroy lists.
 * ************************************************************** */

int
listNode_AttachAhead (ListNode * listNode, ListNode * nodeBehind)
{
  ListNode *nodeAhead = listNode_GetNodeAhead (nodeBehind);

  listNode->prev = nodeBehind;
  nodeBehind->next = listNode;

  if (nodeAhead)
  {
    nodeAhead->prev = listNode;
    listNode->next = nodeAhead;
  }
  return (1);
}

/* **************************************************************
 * Attach a ListNode behind another ListNode.
 * listNode is put behind of nodeAhead.  If nodeAhead already
 * had a ListNode behind it, then we hook up pointers
 * as necessary.
 *
 * Note that a ListNode has no concept as to what a list is.
 * All it knows how to do is to attach and detach itself from
 * other ListNode's.  Lists will use ListNode's to build
 * and destroy lists.
 * ************************************************************** */

int
listNode_AttachBehind (ListNode * listNode, ListNode * nodeAhead)
{
  ListNode *nodeBehind = listNode_GetNodeBehind (nodeAhead);

  listNode->next = nodeAhead;
  nodeAhead->prev = listNode;

  if (nodeBehind)
  {
    nodeBehind->next = listNode;
    listNode->prev = nodeBehind;
  }
  return (1);
}

/*
 * Attach an entity to a listnode (variable).
 * We do _not_ increment the reference count here, cause we
 * often attach an entity to a variable (listnode) temporarily.
 * Only when
 * entity may belong to someother variable.
 */

ListNode *
listNode_AttachEnt (ListNode * lnode, void *ent)
{
  lnode->ent = ent;
  return (lnode);
}

ListNode *
listNode_AttachListEnt (ListNode * lnode, void *lent)
{
  lnode->listent = lent;
  return (lnode);
}

/*
 * Get the double value from the entity attached
 * to a variable (listNode).
 * Return 0.0 if there is no entity attached.
 */

double
var_Double (ListNode * var)
{
  Ent *ent = var->ent;
  if (ent)
  {
    return (ent->d);
  }
  else
  {
    return (0.0);
  }
}

/*
 * Get the void data pointer from the entity
 * attached to a variable (listNode).
 * Return NULL if no entity or data attached.
 */

void *
var_Data (ListNode * var)
{
  Ent *ent = var->ent;
  if (ent)
  {
    return (ent->data);
  }
  else
    return (0);
}

/* ************************************************************************
 * PURPOSE: Detach a ListNode from a list of ListNodes.
 *          The list-node is detached ONLY, not destroyed. The adjoining
 *          nodes (if there are any) are updated to reflect the
 *          detachment.
 * ************************************************************************ */

int listNode_Detach (ListNode * listNode)
{
  ListNode *nodeAhead = listNode_GetNextNode (listNode);
  ListNode *nodeBehind = listNode_GetPrevNode (listNode);

  if (nodeAhead)
    nodeAhead->prev = nodeBehind ? nodeBehind : NULL;

  if (nodeBehind)
    nodeBehind->next = nodeAhead ? nodeAhead : NULL;

  listNode->prev = listNode->next = NULL;

  return (1);
}

/* ************************************************************************
 * Set the "key" identifier in a list-node to the value of string.
 * ************************************************************************ */
void listNode_SetKey (ListNode * listn, char *string_ptr)
{
  if (listn->key)
    GC_FREE (listn->key);
//   listn->key = string_ptr;
  listn->key = cpstr(string_ptr);
}
