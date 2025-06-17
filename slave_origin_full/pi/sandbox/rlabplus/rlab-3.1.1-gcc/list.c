/* list.c */

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
#include "util.h"
#include "list.h"
#include "mem.h"

#include <stdio.h>
#include <string.h>

#undef  THIS_FILE
#define THIS_FILE "list.c"

static int list_setNumNodes _PROTO ((List * list, int num));
static int list_setFirstNode _PROTO ((List * list, ListNode * first));
static int list_setLastNode _PROTO ((List * list, ListNode * last));

/* ************************************************************************
 * To dynamically allocate and initialize an instance of a List.
 * ************************************************************************ */
List * list_Create (void)
{
  List *newlist = (List *) GC_MALLOC (sizeof (List));
  if (newlist == 0)
    rerror ("out of memory");
  list_Initialize (newlist);
  return (newlist);
}

/* ************************************************************************
 * Free all of the dynamically allocated memory associated with
 * a List.  If the List has any ListNode's on it, then
 * list_DestroyAllNodes() is called.
 * ************************************************************************ */

void
list_Destroy (List * list)
{
  ASSERT (list);
  {
    if (list_GetNumNodes (list) > 0)
      list_DestroyAllNodes (list);
    list->name = 0;
    list->type = 0;
    GC_FREE (list);
  }
}

/* ************************************************************************
 * Free all memory associated with all of the ListNode's that
 * currently reside on a List.  The difference between
 * this function and list_Destroy() is that the List itself
 * is not deallocated. Since we're destroying all the list members
 * lets try and be fast.
 * ************************************************************************ */

void
list_DestroyAllNodes (List * list)
{
  ASSERT (list);
  {
    ListNode *lnode, *nextNode;
    int i, numNodes;

    numNodes = list_GetNumNodes (list);
    lnode = list_GetFirstNode (list);

    for (i = 1; i <= numNodes; i++)
    {
      nextNode = listNode_GetNextNode (lnode);
      listNode_Destroy (lnode);
      lnode = nextNode;
    }

    list_setNumNodes (list, 0);
    list_setFirstNode (list, 0);
    list_setLastNode (list, 0);
  }
}

/* ************************************************************************
 *  Destroy all the nodes on the list, but not the node's data.
 *  There are situations where we want the ListNode's destroyed,
 *  but not the data associated with the ListNode's.
 *  Note that this function is identical to list_DestroyAllNodes(),
 *  except listNode_DestroyNodeOnly() is called in liue of
 *  listNode_Destroy().
 *  Note that this function does not destroy the List itself.
 *  In order to do so, one must call list_Destroy(), after
 *  calling this function.
 * ************************************************************************ */

int
list_DestroyNodesOnly (List * list)
{
  ASSERT (list);
  {
    int numNodes;

    while ((numNodes = list_GetNumNodes (list)))
    {
      /* detach the last node on the list until there are no more nodes */
      ListNode *nextNode = list_DetachNodeByPos (list, numNodes);

      if (nextNode)
	listNode_DestroyNodeOnly (nextNode);
    }
    return (1);
  }
}

/* ************************************************************************
 * Destroy the ListNode that is at position pos in a List.
 * First the list-node is detached, then destroyed, and then
 * the list count is updated.
 * ************************************************************************ */

int
list_DestroyNodeByPos (List * list, int pos)
{
  ASSERT (list);
  {
    if (list_IsPosValid (list, pos))
    {
      ListNode *node2Destroy = list_DetachNodeByPos (list, pos);
      if (node2Destroy)
      {
	listNode_Destroy (node2Destroy);
	return (1);
      }
    }
    return (0);
  }
}

/* ************************************************************************
 * Given a pointer to a list-node, detach and destroy the
 * list-node if it currently resides on the list.
 * ************************************************************************ */

int
list_DestroyNodeByAddr (List * list, ListNode * node2Destroy)
{
  ASSERT (list);
  ASSERT (node2Destroy);
  {
    ListNode *destroyMe = list_DetachNodeByAddr (list, node2Destroy);

    if (destroyMe)
    {
      listNode_Destroy (destroyMe);
      return (1);
    }
  }
  return (0);
}

/* **************************************************************
 * Destory a node on a list, but do NOT destroy the data.
 * ************************************************************** */
int list_DestroyNodeOnlyByAddr (List * list, ListNode * lnode)
{
  ASSERT (list);
  ASSERT (lnode);
  {
    ListNode *destroyMe = list_DetachNodeByAddr (list, lnode);

    if (destroyMe)
    {
      listNode_Destroy (destroyMe);
      return (1);
    }
  }
  return (0);
}

/* ************************************************************************
 * PURPOSE: Given a List, a ListNode, and a position to insert the
 *          ListNode in the List, insert the ListNode into the List.
 *
 *          It may seem rather odd that this function returns a pointer
 *          which was passed into it (the ListNode pointer).  Just
 *          because the caller passed the pointer in, does not mean
 *          that the caller has the pointer's value stored away, so
 *          we return it just in case.
 *
 *          This function gets two pointers first thing:  A pointer to
 *          the node that will be behind the nodeToInsert, and a
 *          pointer to the node that will be ahead of the nodeToInsert.
 *
 *          If there is a nodeBehind, then we know that there are some
 *          ListNode's currently on the list, and also that we
 *          are inserting at a position other than 1.  Therefore, we
 *          insert the node ahead of the nodeBehind.
 *          If there is a nodeBehind but not a nodeAhead, then that
 *          means that nodeToInsert is the last node in the list,
 *          and therefore we call list_setLastNode() to take care of it.
 *
 *          If there is no nodeBehind, then we are either inserting
 *          at pos 1, or there are no ListNode's currently on the
 *          list.  So, if:
 *
 *          nodeAhead is not NULL:  pos must be 1.
 *          nodeAhead is NULL:      There are no nodes on the list.
 *
 *          If nodeAhead is not NULL, then we attach nodeToInsert
 *          behind nodeAhead.
 *          If nodeAhead is NULL, then we set the last node to
 *          nodeToInsert (because there were no nodes on the list,
 *          nodeToInsert becomes the first and also the last).
 *
 *          Regardless of nodeAhead, if nodeBehind is NULL, we
 *          set the first node equal to nodeToInsert.
 *
 *          Finally, we set the number of nodes in the list equal
 *          to 1 more than the current number of nodes.
 * ************************************************************************ */

ListNode * list_InsertNode (List * list, ListNode * nodeToInsert, int pos)
{
  ASSERT (list);
  ASSERT (nodeToInsert);
  ASSERT (pos > 0 && pos <= (list_GetNumNodes (list) + 1));
  ASSERT (!listNode_IsNodeAttached (nodeToInsert));
  {
    ListNode *nodeBehind = list_GetNodeByPos (list, pos - 1);
    ListNode *nodeAhead = list_GetNodeByPos (list, pos);

    if (nodeBehind)
    {
      /* There is at least 1 node on list go ahead an stuff node in place */
      listNode_AttachAhead (nodeToInsert, nodeBehind);
      if (!nodeAhead)
	list_setLastNode (list, nodeToInsert);
    }
    else
    {
      /* We are either inserting at pos = 1, or there are no nodes on list */
      if (nodeAhead)
      {
	listNode_AttachBehind (nodeToInsert, nodeAhead);
      }
      else
      {
	list_setLastNode (list, nodeToInsert);
	list_setFirstNode (list, nodeToInsert);
      }
    }
    list_setNumNodes (list, list_GetNumNodes (list) + 1);
    return (nodeToInsert);
  }
}

/* ************************************************************************
 * Search the list for the node having a key that matches the searchkey
 * Return a pointer to the ListNode. If the node with the matching key is
 * NOT found return a NULL pointer.
 * ************************************************************************ */

#undef  THIS_SOLVER
#define THIS_SOLVER "list_GetNodeByKey"
ListNode * list_GetNodeByKey (List * list, char *key)
{
  if (!list)
    return (0);

  ListNode *lnode = list_GetFirstNode (list);
  while (lnode)
  {
//     printf(THIS_FILE ": " THIS_SOLVER ": key=%s ? lnode->key=%s\n", key, lnode->key);

    if (strcmp (key, lnode->key) == 0)
      return (lnode);

    lnode = listNode_GetNextNode (lnode);
  }
  return (0);
}

#undef  THIS_SOLVER
#define THIS_SOLVER "list_GetNodeByKey"
int list_GetPosByKey (List * list, char *key)
{
  if (!list)
    return (-1);

  int pos=0;
  ListNode *lnode = list_GetFirstNode (list);
  while (lnode)
  {
    pos++;
    if (strcmp (key, lnode->key) == 0)
      return (pos);
    lnode = listNode_GetNextNode (lnode);
  }
  return (0);
}

#undef  THIS_SOLVER
#define THIS_SOLVER "list_Install"
ListNode * list_Install (List * list, char *key, Ent *ent)
{
//   printf(THIS_FILE ": " THIS_SOLVER ": IN\n");
  int comp;
  ListNode *node = list_GetFirstNode (list), *prev=0;

  ListNode *rval = listNode_Create();
  rval->key = cpstr(key);
  rval->ent = ent;

  // empty list
  if (!node)
  {
//     printf(THIS_FILE ": " THIS_SOLVER ": NULL LIST Insert %s at 1st place\n", key);
    return list_PushNode (list, rval);
  }
    

  // list where key is before the first entry
  comp = strcmp (key, node->key);
  if (comp < 0)
  {
//     printf(THIS_FILE ": " THIS_SOLVER ": Insert %s at 1st place\n", key);
    return list_PushNode (list, rval);
  }
  else if (!comp) // can't have nodes with same name
    return (0);

//   printf(THIS_FILE ": " THIS_SOLVER ": SEARCH for %s\n", key);

  prev = node;
  node = listNode_GetNextNode (node);
  while (node)
  {

    comp = strcmp (key, node->key);
    if (!comp) // can't have nodes with same name
      return (0);

//     printf(THIS_FILE ": " THIS_SOLVER ": SEARCH comparison between %s and %s = %i\n",
//            key, node->key, comp);


    if (comp < 0)
    {
//       printf(THIS_FILE ": " THIS_SOLVER ": SEARCH Insert between %s and %s\n",
//              prev->key, node->key);

      rval->prev = prev;
      rval->next = node;
      prev->next = rval;

      if (node)
        node->prev = rval;

      list->numNodes++;
      return rval;
    }

    prev = node;
    node = listNode_GetNextNode (node);
  }

//   printf(THIS_FILE ": " THIS_SOLVER ": SEARCH: Insert %s at the last place\n", key);

  // add behind last node
  rval->prev = prev;
  rval->next = 0;
  prev->next = rval;
  list->numNodes++;
  list->lastNode = rval;
  return (rval);
}

ListNode * list_InstallListNode (List * list, ListNode *rval)
{
  int comp;
  ListNode *node = list_GetFirstNode (list), *prev=0;
  char *key = rval->key;

  // empty list
  if (!node)
  {
    return list_PushNode (list, rval);
  }

  // list where key is before the first entry
  comp = strcmp (key, node->key);
  if (comp < 0)
  {
    return list_PushNode (list, rval);
  }
  else if (!comp) // can't have nodes with same name
    return (0);

  prev = node;
  node = listNode_GetNextNode (node);
  while (node)
  {

    comp = strcmp (key, node->key);
    if (!comp) // can't have nodes with same name
      return (0);

    if (comp < 0)
    {
      rval->prev = prev;
      rval->next = node;
      prev->next = rval;

      if (node)
        node->prev = rval;

      list->numNodes++;
      return rval;
    }

    prev = node;
    node = listNode_GetNextNode (node);
  }
  // add behind last node
  rval->prev = prev;
  rval->next = 0;
  prev->next = rval;
  list->numNodes++;
  list->lastNode = rval;
  return (rval);
}

/* ************************************************************************
 * Detach a ListNode from a List, given the ListNode's
 * name.  If the position is not valid, * NULL is returned.
 * ************************************************************************ */
ListNode * list_DetachNodeByKey (List * list, char * key)
{
  if (!list)
    return (0);

  if (!key)
    return (0);

  ListNode *lnode = list_GetFirstNode (list), *prev=0, *next=0;

  if (!lnode)
    return (0);

  do
  {
    char *lkey = listNode_GetKey (lnode);
    if (lkey)
    {
      if (!strcmp (key,lkey))
      {
        prev = lnode->prev;
        next = lnode->next;

        if (prev)
          prev->next = next;
        else
          list->firstNode = next;

        if (next)
          next->prev = prev;
        else
          list->lastNode = prev;

        list->numNodes--;
        return lnode;
      }
    }

    lnode = listNode_GetNextNode (lnode);
  }
  while (lnode);

  return (0);
}


/* ************************************************************************
 * PURPOSE:     This function inserts a ListNode onto a List at the end
 *              of the List.  Note that list_EnqueueNode() and
 *              list_DequeueNode() may be used to emulate a queue with a
 *              linked list.
 * ************************************************************************ */

ListNode *
list_EnqueueNode (List * list, ListNode * nodeToInsert)
{
  ASSERT (list);
  ASSERT (nodeToInsert);
  return list_InsertNode (list, nodeToInsert, list_GetNumNodes (list) + 1);
}

/* ************************************************************************
 * PURPOSE: "Dequeue" a ListNode from a List.  list_EnqueueNode() and
 *          list_DequeueNode() may be used together to make a List
 *          emulate a queue.  list_EnqueueNode() adds a ListNode
 *          to the end of the List, and list_DequeueNode() detachs
 *          the first ListNode in the List.
 *
 * ************************************************************************ */

ListNode *
list_DequeueNode (List * list)
{
  ASSERT (list);
  return list_DetachNodeByPos (list, 1);
}

/* ************************************************************************
 * Push a ListNode onto the front of a List. Do Not use list_InsertNode,
 * we want this to be faster, so we are more specific.
 * ************************************************************************ */

ListNode * list_PushNode (List * list, ListNode * pnode)
{
  ASSERT (list);
  {
    if (list->numNodes == 0)
    {
      list->firstNode = list->lastNode = pnode;
      pnode->next = 0;
    }
    else
    {
      ListNode *lnode = list->firstNode;
      pnode->next = lnode;
      lnode->prev = pnode;
      list->firstNode = pnode;
    }
    pnode->prev = 0;
    list->numNodes++;
    return (pnode);
  }
}

/* ************************************************************************
 * PURPOSE: Pop a ListNode from a List.  Popping a ListNode, means
 *          detaching the FIRST ListNode in the List. Push node puts
 *          on the front of the list, and pop removes from the front.
 * ************************************************************************ */

ListNode *
list_PopNode (List * list)
{
  ASSERT (list);
  return (list_DetachNodeByPos (list, 1));
}

/* ************************************************************************
 * Detach a ListNode from a List, given the ListNode's
 * position in the List.  If the position is not valid,
 * NULL is returned.
 * ************************************************************************ */
ListNode * list_DetachNodeByPos (List * list, int pos)
{
  ASSERT (list);
  {
    ListNode *detachedNode = NULL;
    int numNodesOnList = list_GetNumNodes (list);

    if (list_IsPosValid (list, pos))
    {
      if ((detachedNode = list_GetNodeByPos (list, pos)))
      {
        if (pos == 1)
          list_setFirstNode (list, list_GetNodeByPos (list, 2));

        if (pos == numNodesOnList)
          list_setLastNode (list, list_GetNodeByPos (list, numNodesOnList - 1));

        listNode_Detach (detachedNode);
        list_setNumNodes (list, list_GetNumNodes (list) - 1);
      }
      return (detachedNode);
    }
    return (NULL);
  }
}

/* ************************************************************************
 * Given a pointer to a ListNode, detach the ListNode from
 * the List if the ListNode currently resides in the list.
 * Return a the pointer to the detached node, NULL if failure.
 * ************************************************************************ */

ListNode *
list_DetachNodeByAddr (List * list, ListNode * lnode)
{
  ASSERT (list);
  ASSERT (lnode);
  {
    int n = list_GetNumNodes (list);

    if (lnode == list_GetFirstNode (list))
    {
      list_setFirstNode (list, listNode_GetNextNode (lnode));
      if (lnode == list_GetLastNode (list))
	list_setLastNode (list, listNode_GetPrevNode (lnode));
      listNode_Detach (lnode);
      list->numNodes = n - 1;
      return (lnode);
    }
    else if (lnode == list_GetLastNode (list))
    {
      list_setLastNode (list, listNode_GetPrevNode (lnode));
      listNode_Detach (lnode);
      list->numNodes = n - 1;
      return (lnode);
    }
    else
    {
      listNode_Detach (lnode);
      list->numNodes = n - 1;
      return (lnode);
    }
  }
}

/* ************************************************************************
 * Return a pointer to the ListNode that resides at the specified position
 * in the list. Return a NULL if the position does not exist.
 * ************************************************************************ */
#undef  THIS_SOLVER
#define THIS_SOLVER "list_GetNodeByPos"
ListNode * list_GetNodeByPos (List * list, int pos)
{
  ASSERT (list);
  ASSERT (pos >= 0);
  {
//     printf(THIS_FILE ": " THIS_SOLVER ": IN\n");
    ListNode *nextNode;
    int i;

    if (!list_IsPosValid (list, pos))
      return (0);

    nextNode = list_GetFirstNode (list);

    for (i = 1; i != pos; i++)
    {
      nextNode = listNode_GetNextNode (nextNode);
    }
//     printf(THIS_FILE ": " THIS_SOLVER ": %s at %i\n", nextNode->key, pos);

//     printf(THIS_FILE ": " THIS_SOLVER ": OUT\n");
    return (nextNode);
  }
}

/* ************************************************************************
 * PURPOSE: Return the data contained in the ListNode at position
 *          pos in a List.  Notice that, unlike list_PopNodeData()
 *          and list_DequeueNodeData(), this function does not detach
 *          the ListNode and free it.  The List is untouched
 *          by this operation.
 *
 * Get the pointer to the list-node's data. Let list_GetNodeByPos do the
 * error checking.
 * ************************************************************************ */

void *
list_GetNodeDataByPos (List * list, int pos)
{
  ListNode *listNode;
  void *data;

  if (list_IsPosValid (list, pos))
  {
    listNode = list_GetNodeByPos (list, pos);
    return (data = listNode_GetEnt (listNode));
  }
  return (0);
}

/* ***********************************************************************
 * PURPOSE:
 *          Print the contents of a List to a file.
 * *********************************************************************** */
extern void class_print (Ent * e);
void list_PrintToFile (List * list, FILE * stream)
{
  ASSERT (list);
  ASSERT (stream);
  {
    ListNode *nextNode;
    int j, pos;

    pos = 1;
    while (pos <= list_GetNumNodes (list))
    {
      fprintf (stream, "   ");
      for (j = 0; j < 5; j++)
      {
        if ((nextNode = list_GetNodeByPos (list, pos++)))
        {
          fprintf (stream, "%-15s", listNode_GetKey (nextNode));
          class_print((Ent *) var_ent(nextNode));
        }
      }
      fprintf (stream, "\n");
    }
  }
}

void list_PrintScopeToFile (List * list, FILE * stream)
{
  if (!list)
    return;

  if (!stream)
    return;

  ListNode *nextNode;
  int j, pos;

  pos = 1;
  while (pos <= list_GetNumNodes (list))
  {
    fprintf (stream, "   ");
    for (j = 0; j < 5; j++)
    {
      if ((nextNode = list_GetNodeByPos (list, pos++)))
        fprintf (stream, "%i ", nextNode->scope);
    }
    fprintf (stream, "\n");
  }
}

/* ************************************************************************
 * Return a pointer to the last node in a _list struct.
 * ************************************************************************ */

ListNode *
list_GetLastNode (List * list)
{
  if (list)
    return list->lastNode;
  return (0);
}

/* ************************************************************************
 * Return a pointer to the first node in a List struct.
 * ************************************************************************ */

ListNode *
list_GetFirstNode (List * list)
{
  if (list)
    return list->firstNode;
  return (0);
}

/* ************************************************************************
 * PURPOSE:     Using private list-class functions,
 *              set all members a List to initial values.
 * ************************************************************************ */

int
list_Initialize (List * list)
{
  ASSERT (list);
  {
    list_setNumNodes (list, 0);
    list_setFirstNode (list, NULL);
    list_setLastNode (list, NULL);
    list->type = LIST;
    list->name = 0;
  }
  return (1);
}

/* ************************************************************************
 * Set the number of nodes in a List struct.
 * ************************************************************************ */

static int
list_setNumNodes (List * list, int numNodes)
{
  ASSERT (list);
  ASSERT (numNodes >= 0);
  list->numNodes = numNodes;
  return (1);
}

/* ************************************************************************
 * Set the firstNode pointer in a List.
 * ************************************************************************ */

static int
list_setFirstNode (List * list, ListNode * firstNode)
{
  ASSERT (list);
  list->firstNode = firstNode;
  return (1);
}

/* ************************************************************************
 * Set the lastNode pointer in a List.
 * ************************************************************************ */

static int
list_setLastNode (List * list, ListNode * lastNode)
{
  ASSERT (list);
  list->lastNode = lastNode;
  return (1);
}

/* ************************************************************************
 * Check to see if a given position is valid for a given List.
 * ************************************************************************ */

int
list_IsPosValid (List * list, int pos)
{
  return (pos > 0 && pos <= list_GetNumNodes (list));
}
