/* list.h */

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

#ifndef  LIST_H
#define  LIST_H

#include "rlab.h"
#include "listnode.h"

struct _list
{
  int type;			/* identifies the object */
  char *name;
  int numNodes;			/* The number of ListNode's in the list   */
  ListNode *firstNode;		/* The first ListNode in the list         */
  ListNode *lastNode;		/* The last ListNode in the list          */
};

typedef struct _list List;

extern List *list_Create (void);
extern void list_Destroy (List *);
extern void list_DestroyAllNodes (List *);
extern int list_DestroyNodesOnly (List *);
extern int list_DestroyNodeByPos (List *, int);
extern int list_DestroyNodeByAddr (List *, ListNode *);
extern int list_DestroyNodeOnlyByAddr (List *, ListNode *);
extern int list_Initialize (List * list);

extern ListNode *list_InsertNode (List *, ListNode *, int);
extern ListNode *list_Install    (List *, char *, Ent *);
extern ListNode *list_InstallListNode (List *, ListNode *);
extern int list_GetPosByKey (List * list, char *key);

extern ListNode *list_DetachNodeByKey (List *, char *);

extern ListNode *list_EnqueueNode (List *, ListNode *);
extern ListNode *list_DequeueNode (List *);

extern ListNode *list_PushNode (List *, ListNode *);
extern ListNode *list_PopNode (List *);

extern ListNode *list_DetachNodeByPos (List *, int);
extern ListNode *list_DetachNodeByAddr (List *, ListNode *);

extern ListNode *list_GetNodeByPos (List *, int);
extern void *list_GetNodeDataByPos (List *, int);

extern void list_PrintScopeToFile (List * list, FILE * stream);
extern void list_PrintToFile (List *, FILE *);
extern int list_IsPosValid (List *, int);

extern ListNode *list_GetLastNode (List *);
extern ListNode *list_GetFirstNode (List *);
extern ListNode *list_GetNodeByKey (List *, char *);

#define list_GetNumNodes(list)  (list ? (((List *)(list))->numNodes) : 0)

#endif /* LIST_H */
