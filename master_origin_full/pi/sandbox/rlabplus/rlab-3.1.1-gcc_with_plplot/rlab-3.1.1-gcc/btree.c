/* btree.c */

/* Binary tree implementation using ListNode structure for tree nodes.
   The ListNode structure was originally designed for doubley linked
   lists, so for binary tree purposes, left == prev, and right == next. */

/*  This file is a part of RLaB ("Our"-LaB)
   Copyright (C) 1992, 1994  Ian R. Searle

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
#include "class.h"
#include "symbol.h"
#include "btree.h"
#include "function.h"
#include "util.h"
#include "mem.h"
#include "listnode.h"
#include "mdr.h"
#include "mds.h"
#include "rlab_solver_parameters_names.h"

#undef  THIS_FILE
#define THIS_FILE "btree.c"

#include <stdio.h>

/* **************************************************************
 * Create a binary tree structure and initialize it.
 * ************************************************************** */

Btree *
btree_Create (void)
{
  Btree *new;
  new = (Btree *) GC_MALLOC (sizeof (Btree));
  if (new == 0)
    rerror ("out of memory");
  new->numNodes  = 0;
  new->root_node = 0;
  new->isconst   = 0;
  return (new);
}

/* **************************************************************
 * Destroy a binary tree, and all of the nodes in it.
 * ************************************************************** */

static void btree_destroy_nodes (ListNode * node)
{
  if (node != 0)
  {
    btree_destroy_nodes (node->prev);
    btree_destroy_nodes (node->next);
    ent_Destroy (var_ent (node));
    ent_Destroy (var_listent (node));
    listNode_Destroy (node);
  }
}

void btree_Destroy (Btree * root)
{
  ASSERT (root);
  {
    if (root->isconst)
    {
      fprintf(stderr,"btree_Destroy: protected list cannot be destroyed!");
      return;
    }
    btree_destroy_nodes (root->root_node);
    root->numNodes = 0;
    root->root_node = 0;
    GC_FREE (root);
    root = 0;
  }
}

void btree_listnode_Destroy (ListNode * rnode)
{
  btree_destroy_nodes (rnode);
}

/* **************************************************************
 * Count the number of nodes in a binary tree. We must count
 * ourselves, cause we don't want to include UNDEF nodes.
 * ************************************************************** */

static int nnodes = 0;
static void btree_count_nodes (ListNode * node)
{
  if (node != 0)
  {
    btree_count_nodes (node->prev);
    if (ent_type (var_ent (node)) != UNDEF)
    {
      nnodes++;
    }
    btree_count_nodes (node->next);
  }
}

int btree_CountNodes (Btree * bt)
{
  nnodes = 0;
  btree_count_nodes (bt->root_node);
  return (nnodes);
}

int btree_listnode_CountNodes(ListNode *rnode)
{
  nnodes = 0;
  btree_count_nodes (rnode);
  return (nnodes);
}

static int nnodes_all = 0;
static void btree_count_all_nodes (ListNode * node)
{
  if (node != 0)
  {
    btree_count_nodes (node->prev);
    nnodes_all++;
    btree_count_nodes (node->next);
  }
}

int btree_CountAllNodes (Btree * bt)
{
  nnodes_all = 0;
  btree_count_all_nodes (bt->root_node);
  return (nnodes);
}

/* **************************************************************
 * Add a new node to the given binary tree.
 * ************************************************************** */

static ListNode *btree_add_node (ListNode * rnode, ListNode * anode);

ListNode * btree_AddNode (Btree * root, ListNode * node)
{
  ListNode *add;

  if (root->root_node == 0)
  {
    /* The first node to be added, special case. */
    root->numNodes++;

    /* Mark the node as "owned". */
    node->scope = 1;

    return (root->root_node = node);
  }
  else
  {
    /* All the rest. */
    if ((add = btree_add_node (root->root_node, node)) != 0)
      root->numNodes++;

    /* Mark the node as "owned". */
    node->scope = 1;

    return (add);
  }
}

static ListNode * btree_add_node (ListNode * rnode, ListNode * anode)
{
  int cond;

  if (!anode)
    return (0);
  char *akey = var_key (anode);
  if (!akey)
    return (0);

  if (!rnode)
  {
    /* End of tree, add node */
    return (anode);
  }
  char *rkey = var_key (rnode);
  if (!rkey)
    return (0);

  if ((cond = strcmp (akey, rkey)) == 0)
  {
    /* Error, found duplicate */
    return (0);
  }
  else if (cond < 0)
  {
    /* Traverse left branch */
    rnode->prev = btree_add_node (rnode->prev, anode);
    anode->next = 0; // enough of this singly linked list crap. daddy has GB's to spare
  }
  else
  {
    /* Traverse right branch */
    rnode->next = btree_add_node (rnode->next, anode);
    anode->prev = 0; // enough of this singly linked list crap. daddy has GB's to spare
  }
  return (rnode);
}

/* **************************************************************
 * Search a tree for the node with the given key.
 * ************************************************************** */

static ListNode *btree_find_node (ListNode * node, char *key);

ListNode * btree_FindNode (Btree * root, char *key)
{
  if (!root)
    return (0);

  return (btree_find_node (root->root_node, key));
}

static ListNode * btree_find_node (ListNode * node, char *key)
{
  int cond;
  if (!key)
    return (0);

  if (!node)
    return (0);

  char *nkey = var_key (node);
  if (!nkey)
    return (0);

  if ((cond = strcmp (key, nkey)) == 0)
    /* Found it */
    return (node);
  else if (cond < 0)
    /* Traverse left branch */
    return (btree_find_node (node->prev, key));
  else
    /* Traverse right branch */
    return (btree_find_node (node->next, key));
}

/* **************************************************************
 * Print binary tree.
 * ************************************************************** */

#define N_WORD     5
static int count = 0;

char ** btree_get_node_names (ListNode * node, char **names)
{
  if (node != 0)
  {
    names = btree_get_node_names (node->prev, names);
    if (ent_type (var_ent (node)) != UNDEF)
    {
      *(names++) = cpstr (listNode_GetKey (node));
      count++;
    }
    names = btree_get_node_names (node->next, names);
  }
  return (names);
}

void btree_Print (Btree * root, FILE * fn)
{
  int i;
  char **names;

  count = 0;

  if (root->numNodes != 0)
  {
    names = (char **) GC_MALLOC ((root->numNodes) * sizeof (char *));
    if (names == 0)
      rerror ("out of memory");
    btree_get_node_names (root->root_node, names);

    if (count)
    {
      /* Now print out names array */
      fprintf (fn, "   ");
      for (i = 1; i <= count; i++)
      {
        fprintf (fn, "%-13s\t", names[i - 1]);
        GC_FREE (names[i - 1]);
        if ((i % N_WORD) == 0)
          fprintf (fn, "\n   ");
      }
      fprintf (fn, "\n");
      GC_FREE (names);
    }
    else
    {
      fprintf (fn, "\t<<>>\n");
    }
  }
  else
    fprintf (fn, "\t<<>>\n");
}

char ** btree_get_all_node_names (ListNode * node, char **names)
{
  if (node != 0)
  {
    names = btree_get_node_names (node->prev, names);
  *(names++) = cpstr (listNode_GetKey (node));
    count++;
    names = btree_get_node_names (node->next, names);
  }
  return (names);
}

void btree_PrintAll (Btree * root, FILE * fn)
{
  int i;
  char **names;

  int nnodes_all = btree_CountAllNodes(root);

  count = 0;

  if (nnodes_all != 0)
  {
    names = (char **) GC_MALLOC (nnodes_all * sizeof (char *));
    if (names == 0)
      rerror ("out of memory");
    btree_get_all_node_names (root->root_node, names);

    if (count)
    {
      /* Now print out names array */
      fprintf (fn, "   ");
      for (i = 1; i <= count; i++)
      {
        fprintf (fn, "%-13s\t", names[i - 1]);
        GC_FREE (names[i - 1]);
        if ((i % N_WORD) == 0)
          fprintf (fn, "\n   ");
      }
      fprintf (fn, "\n");
      GC_FREE (names);
    }
    else
    {
      fprintf (fn, "\t<<>>\n");
    }
  }
  else
    fprintf (fn, "\t<<>>\n");
}


/* **************************************************************
 * Copy a Btree. Just make a new tree, and copy all of the ListNodes,
 * Do not copy the entities, just point to them.
 * ************************************************************** */
void btree_copy (Btree * new, ListNode * node)
{
  Btree *nbt;
  Ent *ent, *new_ent;
  ListNode *new_node;

  if (new->isconst)
    rerror("internal error: cannot add items to a protected tree!");

  if (node == 0)
  {
    /* End of tree */
    return;
  }
  else
  {
    /*
     * Copy the ListNode...
     * But, first, check to see if it is already a list?
     */
    ent = var_ent (node);
    if (ent_type (ent) == BTREE)
    {
      nbt = btree_Copy (ent_data (ent));
      new_ent = ent_Create ();
      ent_data (new_ent) = nbt;
      ent_type (new_ent) = BTREE;
      new_node = listNode_Create ();
      listNode_SetKey (new_node, var_key (node));
      listNode_AttachEnt (new_node, new_ent);
    }
    else
    {
      new_node = listNode_Copy (node);
    }

    if (var_listent (node))
    {
      ent = var_listent (node);
      nbt = btree_Copy (ent_data (ent));
      new_ent = ent_Create ();
      ent_data (new_ent) = nbt;
      ent_type (new_ent) = BTREE;
      listNode_AttachListEnt (new_node, new_ent);
      ent_IncRef (new_ent);
    }

    /* Increment the reference count on the new node. */
    ent_IncRef (var_ent (new_node));

    /* Add newnode to the new tree. */
    btree_AddNode (new, new_node);

    /* Traverse left branch */
    btree_copy (new, node->prev);

    /* Traverse right branch */
    btree_copy (new, node->next);
  }
}

void btree_duplicate (Btree * new, ListNode * node)
{
  Btree *nbt;
  Ent *ent, *new_ent;
  ListNode *new_node;

  if (new->isconst)
    rerror("internal error: cannot add items to a protected tree!");

  if (node == 0)
  {
    /* End of tree */
    return;
  }
  else
  {
    /*
     * Copy the ListNode...
     * But, first, check to see if it is already a list?
     */
    ent = var_ent (node);
    if (ent_type (ent) == BTREE)
    {
      nbt = btree_Copy (ent_data (ent));
      new_ent = ent_Create ();
      ent_data (new_ent) = nbt;
      ent_type (new_ent) = BTREE;
      new_node = listNode_Create ();
      listNode_SetKey (new_node, var_key (node));
      listNode_AttachEnt (new_node, new_ent);
    }
    else
    {
      new_node = listNode_Duplicate (node);
    }

    if (var_listent (node))
    {
      ent = var_listent (node);
      nbt = btree_Copy (ent_data (ent));
      new_ent = ent_Create ();
      ent_data (new_ent) = nbt;
      ent_type (new_ent) = BTREE;
      listNode_AttachListEnt (new_node, new_ent);
      ent_IncRef (new_ent);
    }

    /* Increment the reference count on the new node. */
    ent_IncRef (var_ent (new_node));

    /* Add newnode to the new tree. */
    btree_AddNode (new, new_node);

    /* Traverse left branch */
    btree_duplicate (new, node->prev);

    /* Traverse right branch */
    btree_duplicate (new, node->next);
  }
}

Btree * btree_Copy (Btree * btree)
{
  //if (btree->isconst)
  //  rerror("internal error: cannot copy a protected tree!");

  Btree *new = btree_Create ();
  btree_copy (new, btree->root_node);
  new->isconst = 0;

  return (new);
}


/* **************************************************************
 * Return a list's member.
 * ************************************************************** */

ListNode *
btree_MemberRef (Btree * b, char *name, int *type)
{
  Ent *ne;
  ListNode *nl;

  if ((nl = btree_FindNode (b, name)) == 0)
  {
    /* Didn't find it. */
    ne = ent_Create ();
    ent_SetType (ne, UNDEF);
    nl = install (b, cpstr (name), ne);
  }

  *type = VAR;
  return (nl);
}

/* **************************************************************
 * Return an array of character strings, that contain the names
 * of the lists members.
 * ************************************************************** */

char **
btree_Members (Btree * bt, int *n)
{
  char **members;
  int nnode;

  nnode = btree_CountNodes (bt);
  if (nnode != 0)
  {
    /* We might be overdoing it a bit here. */
    count = 0;
    members = (char **) GC_MALLOC (nnode * sizeof (char *));
    if (members == 0)
      rerror (RLAB_ERROR_OUT_OF_MEMORY);

    btree_get_node_names (bt->root_node, members);

    *n = count;
    return (members);
  }
  *n = 0;
  return (0);
}

/* **************************************************************
 * Return the class of a List.
 * ************************************************************** */

char *
btree_Class (Btree * bt)
{
  return (cpstr (RLAB_MEMBER_CLASS_LIST));
}

/*
 * Return the "real" number of nodes. Count only the ones
 * that are not UNDEF.
 */

static int btree_get_real_num_nodes (ListNode * node, int num);

int
btree_GetRealNumNodes (Btree * root)
{
  int num = 0;
  return (btree_get_real_num_nodes (root->root_node, num));
}

static int
btree_get_real_num_nodes (ListNode * node, int num)
{
  Ent *ent;
  if (node != 0)
  {
    num = btree_get_real_num_nodes (node->prev, num);
    ent = var_ent (node);
    if (ent_type (ent) != UNDEF)
    {
      num++;
    }
    num = btree_get_real_num_nodes (node->next, num);
  }
  return (num);
}

/*
 * Get a BTREEs element or member names.
 * The strings are copied in this function, so
 * the return'ed object is the caller's.
 */

char **
get_btree_node_names (ListNode * node, char **names)
{
  if (node != 0)
  {
    names = get_btree_node_names (node->prev, names);
    if (ent_type (var_ent (node)) != UNDEF)
    {
      *(names++) = cpstr (listNode_GetKey (node));
      count++;
    }
    names = get_btree_node_names (node->next, names);
  }
  return (names);
}

static size_t btree_get_sizeof (ListNode *);

/*
 * Get the sizeof a list, and all of its elements.
 * Like the C-language function (sizeof()).
 */

size_t
btree_Sizeof (Btree * bt)
{
  return (btree_get_sizeof (bt->root_node));
}

static size_t
btree_get_sizeof (ListNode * node)
{
  size_t size = 0;
  if (node != 0)
  {
    size += btree_get_sizeof (node->prev);
    size += class_sizeof (var_ent (node));
    if (var_listent (node))
    {
      size += class_sizeof (var_listent (node));
    }
    size += btree_get_sizeof (node->next);
  }
  return size;
}

int
btree_Sublist (Btree * bt, char *name)
{
  /*
   * For the list class, don't do anyting with name.
   */

  return (0);
}
