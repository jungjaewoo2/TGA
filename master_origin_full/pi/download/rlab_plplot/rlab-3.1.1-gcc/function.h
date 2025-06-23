/* function.h */

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

#ifndef RLAB_FUNCTION
#define RLAB_FUNCTION

#include "rlab.h"
#include "list.h"
#include "mds.h"
#include "btree.h"

struct _function
{
  int type;
  char *name;
  char *classdef_filename;
  char *r_script;
  int n_args;       /* number of arguments */
  List *args;       /* ptr to List of arguments */
  int n_local;      /* number of local variables */
  List *local;      /* ptr to List of local variables */
  int n_global;     /* number of global variables */
  List *global;     /* ptr to List of global variables */
  List *stat;       /* ptr to list of static variables */
  int n_stat;       /* ptr to list of static variables */
  int ncode;        /* size of code */
  Inst *code;       /* ptr to code segment */
  Btree *list;
};

typedef struct _function Function;

/* We have two identical structs at the moment for debugging clarity */

struct _func_arg
{
  int type;
  char *name;
  int offset;			/* stack offset for variable location */
};

typedef struct _func_arg FArg;

struct _local_var
{
  int type;
  char *name;
  int offset;			/* stack offset for variable location */
};

typedef struct _local_var LVar;

// static declaration inside function
extern void  function_SetStatPtr (Function * f, List * stat_list);
extern List *function_GetStatPtr (Function * f);
extern int   function_GetNstat (Function * f);
extern void  function_SetNstat (Function * f, int n);
extern int   function_HasStatVar (Function *f);


extern Function *function_Create (int);
extern void function_Destroy (Function * f);
extern void function_Destroy_bltin (Function * f);
extern Function *function_Copy (Function * f);
extern void function_Print (Function *f, FILE *fptr);
extern void subprog_Print (Function *f, FILE *fptr);

extern void function_SetName (Function * f, char *name);
extern void function_SetCodePtr (Function * f, Inst * ptr);
extern Inst *function_GetCodePtr (Function * f);
extern void function_SetArgPtr (Function * f, List * arg_list);
extern void function_SetLocalPtr (Function * f, List * local_var_list);
extern void function_SetGlobalPtr (Function * f, List * local_var_list);
extern int function_HasLocalVar (Function * f);
extern int function_SetCodeSize (Function * f, int size);
extern int function_GetCodeSize (Function * f);

extern List *function_GetArgPtr (Function * f);
extern void function_SetNargs (Function * f, int n_args);
extern int function_GetNargs (Function * f);
extern void function_SetNlocal (Function * f, int n_local);
extern int function_GetNlocal (Function * f);
extern void function_SetNglobal (Function * f, int n_local);
extern int function_GetNglobal (Function * f);

extern List * function_GetLocalPtr (Function *f);
extern List * function_GetGlobalPtr (Function *f);

extern char *function_GetName (Function * f);
extern void function_SetName (Function * f, char *name);

extern int function_Sublist (Function *f, char *name);

extern LVar *lvar_Create (void);
extern void lvar_Destroy (LVar * l_var);
extern void lvar_SetOffset (LVar * l_var, int offset);
extern void lvar_SetName (LVar * l_var, char *name);
extern char *lvar_GetName (LVar * l_var);

void  var_push (List ** , char *);
List *local_var_push (List * list, char *s);
List *global_var_push (List * list, char *s);
List *arg_var_push (List * list, char *s);

#define lvar_GetOffset(lvar)   (((LVar *)(lvar))->offset)

extern int function_setup1 (int lsave, char *fn);
extern ListNode *function_setup2 (List *a4, List *a10, List *a13, List *, int poff);
extern void function_setup_clear_nf ();
extern ListNode *classdef_setup (List *args, char *fname, char *r_script, int lsave);

extern char *function_Class (Function *f);
extern MDS *function_Type_BF (Function *f);
extern size_t function_Sizeof (Function *f);
extern void *function_MemberRef (Function *f, char *name, int *type);
extern char **function_Members (Function *f, int *n);

#endif /* RLAB_FUNCTION */
