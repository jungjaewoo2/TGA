// /* function.c */

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
#include "code.h"
#include "function.h"
#include "util.h"
#include "mem.h"
#include "symbol.h"
#include "mds.h"
#include "list.h"
#include "listnode.h"
#include "btree.h"
#include "rlab_solver_parameters_names.h"
#include "rlab_macros.h"

#undef  THIS_FILE
#define THIS_FILE "function.c"


void fstatic_var_push (char *file_name, char *name); /* from rlab_parser.y */

/* **************************************************************
 * Create an RLaB Function entity.
 * ************************************************************** */
#undef  THIS_SOLVER
#define THIS_SOLVER "function_Destroy"
Function * function_Create (int is_classdef)
{
  Function *new = (Function *) GC_MALLOC (sizeof (Function));
  if (new == 0)
    rerror ("out of memory");

  if (is_classdef==1)
  {
    new->type = U_CLASS;
  }
  else
  {
    new->type = U_FUNCTION;
  }
  new->name = 0;
  new->classdef_filename=0;

  new->r_script = 0;

  new->n_args = 0;
  new->args = 0;

  new->n_local = 0;
  new->local = 0;

  new->n_global = 0;
  new->global = 0;

  new->stat = 0;
  new->n_stat = 0;

  new->ncode = 0;
  new->code = 0;

  new->list = 0;

  return (new);
}

/* **************************************************************
 * Destroy an RLaB Function entity.
 * ************************************************************** */
#undef  THIS_SOLVER
#define THIS_SOLVER "function_Destroy"
void function_Destroy (Function * f)
{
  int i, n;
  Ent *e;

  f->type = U_FUNCTION;
  if (f->name)
    GC_FREE (f->name);
  if (f->r_script)
    GC_FREE (f->r_script);
  if (f->classdef_filename)
    GC_FREE (f->classdef_filename);

  /*
   * Take care of the function's arguments.
   */

  if (f->args)
  {
    n = list_GetNumNodes (f->args);
    for (i = 1; i <= n; i++)
    {
      e = list_GetNodeDataByPos (f->args, i);
      lvar_Destroy (ent_data (e));
      ent_Destroy (e);
    }
    list_Destroy (f->args);
  }
  f->args = 0;
  f->n_args = 0;

  /*
   * Take care of the function's local variables.
   */

  if (f->local)
  {
    n = list_GetNumNodes (f->local);
    for (i = 1; i <= n; i++)
    {
      e = list_GetNodeDataByPos (f->local, i);
      lvar_Destroy (ent_data (e));
      ent_Destroy (e);
    }
    list_Destroy (f->local);
  }
  f->n_local = 0;
  f->local = 0;

  if (f->stat)
  {
    list_Destroy(f->stat);
  }
  f->n_stat = 0;
  f->stat = 0;

  /*
   * Take care of the function's global variables.
   */
  if (f->global)
  {
    n = list_GetNumNodes (f->global);
    for (i = 1; i <= n; i++)
    {
      e = list_GetNodeDataByPos (f->global, i);
      lvar_Destroy (ent_data (e));
      ent_Destroy (e);
    }
    list_Destroy (f->global);
  }
  f->n_global = 0;
  f->global = 0;

  /*
   * Kill static variables
   */

  f->ncode = 0;
  if (f->code)
    GC_FREE (f->code);
  f->code = 0;

  if (f->list)
  {
    btree_Destroy (f->list);
    f->list = 0;
  }

  GC_FREE (f);
}

void function_Destroy_bltin (Function * f)
{
  ;
}

#undef  THIS_SOLVER
#define THIS_SOLVER "function_Copy"
Function * function_Copy (Function * f)
{
  char *key=0;
  int i;
  Function *new = function_Create (f->type == U_CLASS);

  /* Create/Copy the argument list */
  new->n_args = f->n_args;
  new->args = list_Create ();
  if (new->n_args > 0)
  {
    for (i = 1; i <= f->n_args; i++)
    {
      key = listNode_GetKey (list_GetNodeByPos (f->args, f->n_args - i + 1));
      arg_var_push (new->args, cpstr (key));
    }
  }

  /* Create/Copy the local var list */
  new->n_local = f->n_local;
  new->local = list_Create ();
  if (new->n_local > 0)
  {
    for (i = 1; i <= f->n_local; i++)
    {
      key = listNode_GetKey (list_GetNodeByPos (f->local, f->n_local - i + 1));
      local_var_push (new->local, cpstr (key));
    }
  }

  /* Create/Copy the global var list */
  new->n_global = f->n_global;
  new->global = list_Create ();
  if (new->n_global > 0)
  {
    for (i = 1; i <= f->n_global; i++)
    {
      key = listNode_GetKey (list_GetNodeByPos (f->global, f->n_global - i + 1));
      local_var_push (new->global, cpstr (key));
    }
  }

  /* Create/Copy the function static var list */
  new->n_stat = f->n_stat;
  if (f->n_stat > 0)
  {
    List *new_stat_list=list_Create ();
    ListNode *lnode=list_GetLastNode(f->stat);

    for (i = 1; i <= f->n_stat; i++)
    {
      list_PushNode (new_stat_list, listNode_Duplicate(lnode));
      lnode=listNode_GetNodeBehind(lnode);
    }
    new->stat = new_stat_list;
    new_stat_list=0;
  }

  new->ncode = f->ncode;
  new->code = (Inst *) GC_MALLOC (sizeof (Inst) * (new->ncode + 1));
  if (new->code == 0)
    rerror (THIS_FILE ": " THIS_SOLVER ": " RLAB_ERROR_OUT_OF_MEMORY);

  /* Copy the compiled code */
  for (i = 0; i < f->ncode; i++)
    new->code[i] = f->code[i];

  if (f->list)
  {
    new->list = btree_Copy (f->list);
  }

  return (new);
}


extern void diss_assemble (Inst * p, int pstop);

void function_Print (Function * f, FILE * fptr)
{
  if (get_print_machine ())
  {
    diss_assemble (function_GetCodePtr (f), function_GetCodeSize (f));
  }
  else
  {
    fprintf (fptr, "\t<user-function>\n");
  }
}

void subprog_Print (Function * f, FILE * fptr)
{
  if (get_print_machine ())
  {
    diss_assemble (function_GetCodePtr (f), function_GetCodeSize (f));
  }
  else
  {
    fprintf (fptr, "\t<user-subprogram>\n");
  }
}

/*--------------------------------------------------------------------*/

char *
function_GetName (Function *f)
{
  ASSERT (f);
  {
    return (f->name);
  }
}

void function_SetName (Function *f, char *name)
{
  ASSERT (f);
  {
    if (f->name)
      GC_FREE (f->name);
    f->name = name;
  }
}

int
function_Sublist (Function * f, char *name)
{
  if (!strcmp (name, "file"))
  {
    return (1);
  }
  else if (!strcmp (name, "class"))
  {
    return (1);
  }
  else if (!strcmp (name, "type"))
  {
    return (1);
  }
  return (0);
}

/*--------------------------------------------------------------------*/

int function_HasLocalVar (Function *f)
{
  if (!f)
    return 0;
  if (f->n_local > 0)
    return (1);
  else
    return (0);
}

int function_HasStatVar (Function *f)
{
  if (f)
    if (f->stat)
      if (f->n_stat > 0)
        return (1);
  return (0);
}

/* **************************************************************
 * Set the function code pointer
 * ************************************************************** */
void function_SetCodePtr (Function *f, Inst *ptr)
{
  f->code = ptr;
}

Inst * function_GetCodePtr (Function *f)
{
  return (f->code);
}

int function_SetCodeSize (Function *f, int size)
{
  ASSERT (f);
  {
#ifdef HAVE_GC
    f->code = (Inst *) GC_MALLOC (size * (sizeof (Inst)));
#else
    f->code = (Inst *) CALLOC (size, (sizeof (Inst)));
#endif
    if (f->code == 0)
      rerror ("out of memory");

    f->ncode = size;
    return (1);
  }
}

int function_GetCodeSize (Function *f)
{
  ASSERT (f);
  {
    return (f->ncode);
  }
}

/*--------------------------------------------------------------------*/

void function_SetArgPtr (Function * f, List * arg_list)
{
  f->args = arg_list;
}

void
function_SetLocalPtr (Function * f, List * local_var_list)
{
  f->local = local_var_list;
}

void
function_SetGlobalPtr (Function * f, List * global_var_list)
{
  f->global = global_var_list;
}

void function_SetStatPtr (Function * f, List * stat_list)
{
  f->stat = stat_list;
}


/*--------------------------------------------------------------------*/

List *
function_GetArgPtr (Function * f)
{
  return (f->args);
}

List *
function_GetLocalPtr (Function * f)
{
  return (f->local);
}

List *
function_GetGlobalPtr (Function * f)
{
  return (f->global);
}

List * function_GetStatPtr (Function * f)
{
  return (f->stat);
}


/*--------------------------------------------------------------------*/

int function_GetNlocal (Function * f)
{
  return (f->n_local);
}

int function_GetNglobal (Function * f)
{
  return (f->n_global);
}

int function_GetNstat(Function * f)
{
  return (f->n_stat);
}

/*--------------------------------------------------------------------*/

void
function_SetNargs (Function * f, int n_args)
{
  f->n_args = n_args;
}

int
function_GetNargs (Function * f)
{
  return (f->n_args);
}

void
function_SetNlocal (Function * f, int n_local)
{
  f->n_local = n_local;
}

void
function_SetNglobal (Function * f, int n_global)
{
  f->n_global = n_global;
}

void function_SetNstat (Function * f, int n)
{
  f->n_stat = n;
}



/* **************************************************************
 * Create a local variable object.
 * ************************************************************** */

LVar * lvar_Create (void)
{
  LVar *new = (LVar *) GC_MALLOC (sizeof (LVar));
  if (new == 0)
    rerror ("out of memory");

  new->type = LOCAL_VAR;
  new->name = 0;
  new->offset = 0;

  return (new);
}

/* **************************************************************
 * Destroy an instance of a LVar.
 * ************************************************************** */

void
lvar_Destroy (LVar * l_var)
{
  if (l_var)
  {
    l_var->type = 0;
    if (l_var->name)
      GC_FREE (l_var->name);
    l_var->offset = 0;

    GC_FREE (l_var);
  }
}

/* **************************************************************
 * Set the offset for a LVar
 * ************************************************************** */

void
lvar_SetOffset (LVar * l_var, int offset)
{
  l_var->offset = offset;
}

void
lvar_SetName (LVar * l_var, char *name)
{
  if (l_var->name)
    GC_FREE (l_var->name);
  l_var->name = cpstr(name);
}

char *
lvar_GetName (LVar * l_var)
{
  return (l_var->name);
}

/* **************************************************************
 * Push a function arg onto a list. At the same time set the arg
 * offset (on the stack).
 * ************************************************************** */

#undef  THIS_SOLVER
#define THIS_SOLVER   "arg_var_push"
List * arg_var_push (List * list, char *s)
{
  Ent *arg_ent=0;
  ListNode *new_node=0;
  LVar *l_var = lvar_Create ();

  if (list == 0)
    list = list_Create ();

  arg_ent = ent_Create ();
  ent_data (arg_ent) = l_var;

  new_node = listNode_Create ();
  listNode_AttachEnt (new_node, arg_ent);
  listNode_SetKey (new_node, s);

  list_PushNode (list, new_node);

  lvar_SetName (l_var, s);
  lvar_SetOffset (l_var, list_GetNumNodes (list));

  return (list);
}

#undef  THIS_SOLVER
#define THIS_SOLVER "var_push"
void var_push (List ** list, char *key)
{
  Ent *local_ent;
  ListNode *new_node;

  if (!*list)
    *list = list_Create();

  local_ent = ent_Create ();
  ent_type(local_ent) = UNDEF;
  ent_data(local_ent) = 0;

  new_node = listNode_Create ();
  listNode_AttachEnt (new_node, local_ent);
  listNode_SetKey (new_node, key);

  list_PushNode (*list, new_node); // put node at position 1 in the list
  var_scope(new_node) = list_GetNumNodes (*list); // put its offset to the n-i
}


/* **************************************************************
 * Push a local var onto a list. At the same time set the variable
 * offset (on the stack).
 * ************************************************************** */

#undef  THIS_SOLVER
#define THIS_SOLVER "local_var_push"
List * local_var_push (List * list, char *s)
{
  Ent *local_ent;
  ListNode *new_node;
  LVar *l_var = lvar_Create ();

  if (list == 0)
    rerror ("Something terrible has happened, no function symbol table");

  local_ent = ent_Create ();
  ent_data (local_ent) = l_var;

  new_node = listNode_Create ();
  listNode_AttachEnt (new_node, local_ent);
  listNode_SetKey (new_node, s);

  list_PushNode (list, new_node);

  lvar_SetName (l_var, s);
  lvar_SetOffset (l_var, list_GetNumNodes (list));

  return (list);
}

List * global_var_push (List * list, char *s)
{
  Ent *global_ent;
  ListNode *new_node;
  LVar *l_var = lvar_Create ();

  if (list == 0)
    rerror ("Something terrible has happened, no function symbol table");

  global_ent = ent_Create ();
  ent_data (global_ent) = l_var;

  new_node = listNode_Create ();
  listNode_AttachEnt (new_node, global_ent);
  listNode_SetKey (new_node, s);

  list_PushNode (list, new_node);

  lvar_SetName (l_var, s);
  lvar_SetOffset (l_var, list_GetNumNodes (list));

  return (list);
}

/* **************************************************************
 * Functions that help create RLaB user-functions. Re-direct program
 * code generation, and help handle local functions someday.
 * ************************************************************** */

extern Program *program_Get (void);
extern Inst *get_program_counter (void);
extern void set_program_counter (Inst * prgm);

#define NFUNC 10		/* Max number of recursive function definition */
static Inst *oldpc[NFUNC];
static Program *np[NFUNC];
static Program *op[NFUNC];
static int nf = 0;


// 1. push classdef arguments on stack
// Now eval(R_SCRIPT)
// 2. push R_SCRIPT onto stack as string variable
// code(OP_PUSH_STRING);
// codep($7);
// 3. find eval function in global symbol table - it is built-in
// if (!var_eval)
//   var_eval = (Var *) lookup (0, "eval");
// code(OP_PUSH_VAR);
// codep(var_eval->var);
// 4. do it
// code(OP_FUNCTION_CALL);
// code(1);   /* only one argument: R_SCRIPT */
// code (OP_DEF_CLASS_RET);
// code (STOP);
#undef  THIS_SOLVER
#define THIS_SOLVER "classdef_setup"
ListNode * classdef_setup (List *args, char *curr_file_name, char *r_script, int lsave)
{
  char class_ptr_name[12];

  Program *npp;
  int counter = 0;
  ListNode * eval_node = (ListNode *) lookup (0, "eval");

  if (!eval_node)
    rerror("Horrible Internal Error: Help ! I need somebody! Help! Just anybody! Help!\n");

  np[nf] = program_Create (15);
  op[nf] = program_Get ();
  oldpc[nf] = get_program_counter ();
  program_Set (np[nf]);

  npp = np[nf];

  //
  // Create unique classdef program:
  //    given r_script
  //    do: eval(r_script)
  //    while collecting class-public, and class-private variables
  // 0:
  npp->prog[counter++].op_code = OP_FILE_NAME;
  // 1:
  npp->prog[counter++].ptr = curr_file_name;
  // 2:
  npp->prog[counter++].op_code = OP_LINE_NO;
  // 3:
  if (lsave == 0 || lsave == 1)
    npp->prog[counter++].op_code = 1;
  else
    npp->prog[counter++].op_code = lsave;

  // 4:
  npp->prog[counter++].op_code = OP_DEF_CLASS_START;

  // 5:
  //  code(OP_PUSH_STRING);
  //  codep(r_script);
  npp->prog[counter++].op_code = OP_PUSH_STRING;
  npp->prog[counter++].ptr = r_script;
  // 6: 'eval'
  //  code(OP_PUSH_VAR);
  //  codep(var_eval->var);
  npp->prog[counter++].op_code = OP_PUSH_VAR;
  npp->prog[counter++].ptr = eval_node;
  // 7: execute function eval(r_script)
  npp->prog[counter++].op_code = OP_FUNCTION_CALL;
  npp->prog[counter++].op_code = 1; // nargs -> 1, that is r_script
  // 8: return from class call
  npp->prog[counter++].op_code = OP_DEF_CLASS_RET;
  // stop
  //npp->prog[counter++].op_code = OP_SAVE_EVAL;
  npp->prog[counter++].op_code = STOP;


  // 10:
  npp->progp = &(npp->prog[counter]);
  set_progoff ();

  // create classdef as function: only thing we know about it
  // is its arguments.
  Function *f = function_Create (1);
  function_SetArgPtr (f, args);
  function_SetNargs  (f, list_GetNumNodes (args));

  // post-processing classdef-private btree in static_tree
  // with key being the address of the Function *f pointer
  // this is the token under which in static_tree the class-private
  // variables will be stored
  sprintf(class_ptr_name,"%p", f);
  f->name = cpstr(class_ptr_name);
  fstatic_var_push (class_ptr_name, NULL);
  f->classdef_filename = cpstr(curr_file_name);

  // Hook the program array into the function.
  f->r_script = cpstr(r_script);
  f->code = npp->prog;
  f->ncode = npp->ncode;

  void *result = (Inst *) GC_MALLOC (sizeof (Inst) * (counter + 1));
  if (!result)
    rerror ("out of memory, save and quit!");
  memcpy (result, f->code, sizeof (Inst) * (counter + 1));
  GC_FREE (f->code);
  f->code = result;

  /* Free up the original program. */
  npp->prog = 0;
  npp->ncode = 0;
  npp->off = 0;
  npp->progp = (Inst *) 0;
  GC_FREE (npp);

  /* Reset program stuff. */
  program_Set (op[nf]);
  set_progoff ();
  set_program_counter (oldpc[nf]);

  Ent *ent = ent_Create ();
  ent_SetType (ent, U_CLASS);
  ent_data (ent) = (void *) f;

  ListNode *var = listNode_Create ();
  listNode_AttachEnt (var, ent);

  /* ent_IncRef (ent); */
  return (var);
}

#undef  THIS_SOLVER
#define THIS_SOLVER "function_setup1"
int function_setup1 (int lsave, char *curr_file_name)
{
//   dprintf("In: nf=%i\n",nf);

  Program *npp;

  /*
   * Create new program, save ptr to existing program
   * and point code generation at new program
   * array.
   */

  np[nf] = program_Create (50);
  op[nf] = program_Get ();
  oldpc[nf] = get_program_counter ();
  program_Set (np[nf]);

  npp = np[nf];

  /*
   * Setup new program array.
   */

  npp->prog[0].op_code = OP_FILE_NAME;
  npp->prog[1].ptr = curr_file_name;
  npp->prog[2].op_code = OP_LINE_NO;
  if (lsave == 0 || lsave == 1)
    npp->prog[3].op_code = 1;
  else
    npp->prog[3].op_code = lsave;

  npp->progp = &(npp->prog[4]);
  set_progoff ();

  nf++;
//   dprintf("OUT\n");

  return (1);
}

// in case the function code has errors, function_setup2 is never
// reached, which may cause segfault
void function_setup_clear_nf ()
{
  Program *npp = np[--nf];

  /* Free up the original program. */
  npp->prog = 0;
  npp->ncode = 0;
  npp->off = 0;
  npp->progp = (Inst *) 0;
  GC_FREE (npp);

  /* Reset program stuff. */
  program_Set (op[nf]);
  set_progoff ();
  set_program_counter (oldpc[nf]);

  return;
}

#undef  THIS_SOLVER
#define THIS_SOLVER "function_setup2"
ListNode * function_setup2 (List * arg, List * local, List * global, List *stat, int poff)
{
  Ent *ent;
  ListNode *var;
  Function *f;
  Program *npp = np[--nf];

  f = function_Create (0);
  ent = ent_Create ();
  ent_SetType (ent, U_FUNCTION);
  ent_data (ent) = (void *) f;

  function_SetArgPtr (f, arg);
  function_SetNargs (f, list_GetNumNodes (arg));
  function_SetGlobalPtr (f, global);
  function_SetNglobal (f, list_GetNumNodes (global));
  function_SetLocalPtr (f, local);
  function_SetNlocal (f, list_GetNumNodes (local));

  // function static variables
  function_SetStatPtr (f, stat);
  function_SetNstat (f, list_GetNumNodes (stat));

  /* Hook the program array into the function. */
  f->code = npp->prog;
  f->ncode = npp->ncode;

 /*
  * Clean up the back of the program array
  * (remove the extra STOPs).
  */
  void *result = (Inst *) GC_MALLOC (sizeof (Inst) * (poff + 1));
  if (!result)
    rerror ("out of memory, save and quit!");
  memcpy (result, f->code, sizeof (Inst) * (poff + 1));
  GC_FREE (f->code);
  f->code = result;

  /* Free up the original program. */
  npp->prog = 0;
  npp->ncode = 0;
  npp->off = 0;
  npp->progp = (Inst *) 0;
  GC_FREE (npp);

  /* Reset program stuff. */
  program_Set (op[nf]);
  set_progoff ();
  set_program_counter (oldpc[nf]);

  var = listNode_Create ();
  listNode_AttachEnt (var, ent);

//   dprintf("OUT\n");

  /* ent_IncRef (ent); */
  return (var);
}

/* **************************************************************
 * Return the class of a function.
 * ************************************************************** */

char *
function_Class (Function * f)
{
  return (cpstr ("function"));
}

MDS *
function_Type_BF (Function * f)
{
  MDS *type = mds_CreateScalar ( "user" );
  return (type);
}

size_t
function_Sizeof (Function * f)
{
  return (size_t) ((double) (f->ncode * sizeof (Inst)));
}

void *
function_MemberRef (Function * f, char *name, int *type)
{
  Ent *ne;
  void *rptr;

  if (!strcmp (name, "file"))
  {
    char *str;
    str = (char *) f->code[1].ptr;
    ne = ent_Create ();
    ent_data (ne) = mds_CreateScalar ( str );
    ent_SetType (ne, MATRIX_DENSE_STRING);
    rptr = (void *) ne;
    *type = ENTITY;
  }
  else if (!strcmp (name, "class"))
  {
    ne = ent_Create ();
    ent_data (ne) = mds_CreateScalar ( "function" );
    ent_SetType (ne, MATRIX_DENSE_STRING);
    rptr = (void *) ne;
    *type = ENTITY;
  }
  else if (!strcmp (name, "type"))
  {
    ne = ent_Create ();
    ent_data (ne) = mds_CreateScalar ( "user" );
    ent_SetType (ne, MATRIX_DENSE_STRING);
    rptr = (void *) ne;
    *type = ENTITY;
  }
  else
  {
    /* This object does not contain one of the above.
     * Signal the caller by passing back 0.
     */
    rptr = 0;
    *type = UNDEF;
  }
  return (rptr);
}

char **
function_Members (Function * f, int *n)
{
  char **marray;

  marray = (char **) GC_MALLOC (3 * sizeof (char *));

  if (marray == 0)
    rerror ("out of memory");

  marray[0] = cpstr ("file");
  marray[1] = cpstr ("class");
  marray[2] = cpstr ("type");

  *n = 3;
  return (marray);
}
