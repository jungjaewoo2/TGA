/* symbol.c */

/*  This file is a part of RLaB ("Our"-LaB)
   Copyright (C) 1992, 1994, 1995  Ian R. Searle

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
#include "symbol.h"
#include "list.h"
#include "listnode.h"
#include "mem.h"
#include "function.h"
#include "btree.h"
#include "util.h"

#ifdef __STDC__
#include <stdlib.h>
#else
#include <malloc.h>
#endif

/*
 * We'll have to have the above until I find out how to fix all
 * of the config header files.
 */
#ifdef HAVE_STDLIB_H
#include <stdlib.h>
#endif

#undef  THIS_FILE
#define THIS_FILE "symbol.c"

extern char *etd (Ent * e);

Btree *global_symbol_table;	/* holds ptr to current global sym tab */
static Btree *symlist = 0;	/* sym table */

/*
 * Place holders for the global symbol table variable.
 */

static ListNode gst_var;
static Ent gst_ent;

/* **************************************************************
 * Create RLaB global symbol-table. Also create a List that hangs
 * on the global symbol table for storing temporary entities.
 * ************************************************************** */

void symbol_table_create (void)
{
  symlist = btree_Create ();
  global_symbol_table = symlist;
  //global_symbol_table->isconst = 1;

  //
  // Fix up the `$$' variable.
  //
  listNode_AttachEnt (&gst_var, &gst_ent);
  listNode_SetOwned (&gst_var);

  ent_data (&gst_ent) = global_symbol_table;
  ent_SetType (&gst_ent, BTREE);
  ent_IncRef (&gst_ent);
  listNode_SetKey (&gst_var, "$$");

}

Btree * get_symtab_ptr (void)
{
  return (global_symbol_table);
}

Ent * get_symtab_ent (void)
{
  return (&gst_ent);
}

/* **************************************************************
 * Install a new entry in symbol table.
 * ************************************************************** */
ListNode * install (void *sym_table, char *s, Ent * ent)
{
  ListNode *lnode, *addnode;

  /* Create the ListNode and set the key */
  lnode = listNode_Create ();
  listNode_SetKey (lnode, s);
  if (ent == UNDEF)
  {
    ent = ent_Create ();
    ent_SetType (ent, UNDEF);
  }

  listNode_AttachEnt (lnode, ent);
  ent_IncRef (ent);

  /* Decide how, and where to install the new ListNode */
  if (sym_table == 0)
    addnode = btree_AddNode (global_symbol_table, lnode);
  else if (sym_table != 0)
    addnode = btree_AddNode ((Btree *) sym_table, lnode);
  else
    rerror ("unanticipated args to install");

  if (!addnode)
    fprintf(stderr, "Terrible internal error: btree_AddNode failed\n");

  return (lnode);
}

/* **************************************************************
 * Find s in symbol table. Return 0 id we DO NOT find s.
 * ************************************************************** */

ListNode * lookup (List * sym_table, char *s)
{
  ListNode *listNode;

  if (sym_table == 0)
  {
    if ((listNode = btree_FindNode (global_symbol_table, s)))
      return (listNode);
  }
  else
  {
    if ((listNode = list_GetNodeByKey (sym_table, s)))
      return (listNode);
  }
  return (0);
}

ListNode * lookup_private (List * sym_table, char *s)
{
  ListNode *listNode = NULL;
  if (sym_table)
    listNode = list_GetNodeByKey (sym_table, s);
  return (listNode);
}

/* **************************************************************
 * Reset the pointer to the global symbol table.
 * ************************************************************** */

void
reset_global_symbol_table_ptr (void)
{
  global_symbol_table = symlist;
}

/* **************************************************************
 * Return a ENTITY that contains a pointer to the global symbol
 * table.
 * ************************************************************** */

Var * gst (void)
{
  Var *retval = (Var *) GC_MALLOC (sizeof (Var));
  if (retval == 0)
    rerror ("out of memory");

  retval->type = GLOBAL;
  retval->var = &gst_var;
  retval->offset = 0;
  retval->name = 0;

  return (retval);
}

/* **************************************************************
 * Fix up the symbol table... straigten up all of the reference
 * counts for each entity. Simply go through the entire symbol table
 * counting how many variables point to each entity. Then make a pass
 * through all the entities, setting the reference counts properly.
 * ************************************************************** */

/*
 * This is a simple struct for a singly-linked list.
 * This simple structure is OK, since we are not in
 * _any_ hurry, this is an error situation, remember.
 */

typedef struct _fix_ent_rc FixEntRC;

struct _fix_ent_rc
{
  Ent *eptr;
  int nc;
  FixEntRC *next;
};

static FixEntRC *fix = 0;
static int nfix = 0;

static FixEntRC *fx_find (Ent * e);
static void count_ent (Ent * e);
static void count_entities (ListNode * ln);
static void fix_entities (ListNode * ln, int debug);
static void count_static_entities_list (List * list);
// static void count_static_entities (Btree * bt);

extern List *get_static_tree (void);

/*
 * The parent function for this operation.
 */
#undef  THIS_SOLVER
#define THIS_SOLVER "fix_symbol_table"
void fix_symbol_table (int debug)
{
  Btree *gst = global_symbol_table;

  /* Dive into the symbol table, counting the entities. */
//   fprintf(stderr, "before count_entities\n");
  count_entities (gst->root_node);
//   fprintf(stderr, "before count_static_entities_list\n");
  count_static_entities_list (get_static_tree ());
//   fprintf(stderr, "before fix_entities\n");
  fix_entities (gst->root_node, debug);
//   fprintf(stderr, "after fix_entities\n");
}

/*
 * Count the entities in the global-symbol-table.
 */

static void count_entities (ListNode * ln)
{
  if (ln != 0)
  {
    if (ent_type (var_ent (ln)) == BTREE)
    {
      /* Special consideration for lists (bin-trees). */
      Ent *etmp = var_ent (ln);
      Btree *bt = ent_data (etmp);
      if (bt == global_symbol_table)
        count_ent (etmp);
      else
        count_entities (bt->root_node);
    }
    count_ent (ln->ent);
    count_entities (ln->next);
    count_entities (ln->prev);
  }
}

/*
 * Track each entity, and how many vars point to it.
 */

static void count_ent (Ent * ent)
{
  FixEntRC *fx, *new;

  fx = fx_find (ent);
  if (fx == 0)
  {
    /* Add new node. */
    new = (FixEntRC *) GC_MALLOC (sizeof (FixEntRC));
    if (new == 0)
      rerror ("TERRIBLE out of memory error");
    new->eptr = ent;
    new->next = fix;
    new->nc = 1;
    fix = new;
    nfix++;
  }
  else
  {
    /* Bump count. */
    fx->nc++;
  }
}

/*
 * We must make special provisions for counting
 * file-static entities.
 */

// static void count_static_ent (ListNode * ln)
// {
//   Btree *btmp;
// 
// //   printf("count_static_ent: etd = %s\n", (char*) etd(ent) );
//   if (ln != 0)
//   {
// //     printf("count_static_ent 1\n");
//     btmp = var_Data (ln);
// //     printf("count_static_ent 2\n");
//     count_entities (btmp->root_node);
// //     printf("count_static_ent 3\n");
//     count_static_ent (ln->next);
// //     printf("count_static_ent 4\n");
//     count_static_ent (ln->prev);
// //     printf("count_static_ent 5\n");
//   }
// }

// static void count_static_entities (Btree * sbt)
// {
//   count_static_ent (sbt->root_node);
// }

static void count_static_entities_list (List * list)
{
  if (list)
  {
    int n = list->numNodes, i;
//     printf("count_static_entities: \tIN list=%p size %i\n", list, n);

    for (i=1; i<=n; i++)
    {
      ListNode *ln = list_GetNodeByPos (list, i);
      if (ln)
      {
//         printf("count_static_entities: \t%s\n", ln->key);
        Btree *btmp = var_Data (ln);
        count_entities (btmp->root_node);
      }
    }
  }

//   printf("count_static_entities: \tOUT\n");
}


/*
 * Find an element in the singly-linked list we are
 * using to track the entity counts.
 * Return 0 if we cannot find it.
 */

static FixEntRC *
fx_find (Ent * e)
{
  FixEntRC *fx, *next;
  fx = fix;
  if (fx == 0)
    return (0);

  if (fx->eptr == e)
    return (fx);

  next = fx->next;
  while (next)
  {
    if (fx->eptr == e)
      return (fx);
    fx = next;
    next = fx->next;
  }
  return (0);
}

/*
 * Go through the list, and fix the reference count on
 * each entity.
 *
 * For the time being, print out some diagnostic
 * information as things are "fixed".
 */


static void
fix_entities (ListNode * ln, int debug)
{
  int i = 1;
  FixEntRC *fx, *next, *ofx;
  fx = fix;
  next = fx->next;

  if (debug)
  {
    fprintf (stderr, "\tChecking/Fixing entity reference counts...\n");
    fprintf (stderr, "\tEntity\t\t\tOld-RC\tNew-RC\tType\n");
  }

  while (next)
  {
    if (fx->eptr != 0)
    {
      /* Print some information for testing/debugging. */
      if (debug)
      {
        fprintf (stderr, "%i\t%p\t\t%i\t%i\t%s\n", i++,
                 fx->eptr,
                 (unsigned int) (fx->eptr->refc), (unsigned int) fx->nc, etd (fx->eptr));
      }
      fx->eptr->refc = fx->nc;
    }
    ofx = fx;
    fx = next;
    next = fx->next;

    GC_FREE (ofx);
  }

  /* Get the last one. */
  if (fx->eptr != 0)
  {
    /* Print some information for testing/debugging. */
    if (debug)
    {
      fprintf (stderr, "%i\t%p\t\t%i\t%i\t%s\n", i++,
	       fx->eptr,
	       (unsigned int) (fx->eptr->refc), fx->nc, etd (fx->eptr));
    }
    fx->eptr->refc = fx->nc;
  }
  GC_FREE (fx);
  fix = 0;
  nfix = 0;

  if (debug)
  {
    fprintf (stderr, ".......Done\n");
  }
}
