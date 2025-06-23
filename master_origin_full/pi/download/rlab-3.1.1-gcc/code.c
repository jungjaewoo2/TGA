/* code.c */

/* code.c: This file contains the majority of the virtual machine.
   The main switch statement is here, as well as most of the
   op-code functions. Some of the class related op-code functions
   are off in files like op.c, op_matrix.c, etc. */

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

#include "rlab.h"
#include "code.h"
#include "list.h"
#include "listnode.h"
#include "class.h"
#include "btree.h"
#include "mem.h"
#include "y.tab.h"
#include "bltin.h"
#include "function.h"
#include "util.h"
#include "rfile.h"
#include "mdc.h"

#include "config.h"

#include <stdio.h>
#include <string.h>


#ifdef THINK_C
#define frame framex
#endif

#undef  THIS_FILE
#define THIS_FILE "code.c"
#include "rlab_solver_parameters_names.h"
#include "rlab_macros.h"
extern FILE *RLAB_STDERR_DS;


/*
 * A data structure to facilitate execution of for-loops.
 * Each for-loop is represented by a struct, and they
 * are linked together as we encounter each loop.
 */

typedef struct _forloop ForLoop;

struct _forloop
{
  int jmp_offset; /* The jump_distance offset because of 'then' keyword */
  int cnt;        /* Current location in vector. */
  int n;          /* Size of vector. */
  int jmp;        /* The jump distance. */
  ListNode *var;  /* The for-loop variable. */
  Ent *ent;       /* The for-loop vector-whatever. */
  ForLoop *next;
};

/*
 * A structure to facilitate user-function execution.
 * This structure is used to form a stack of "frame"
 * pointers. Each pointer representing a function
 * call.
 */

typedef struct _frame Frame;

struct _frame
{
  ListNode *sp;     /* symbol ptr for function */
  Ent *esp;         /* Entity ptr for function */
  Inst *retpc;      /* ptr to machine position for function return */
  Datum *argn;      /* ptr to args on stack */
  int nargs;        /* # of args function is called with */
  List *arg_list;   /* list of user-function arguments */
  Datum *localn;    /* ptr to static and local vars on stack */
  int n_local;      /* number of local variables */
  Datum *statn;     // function static variables on stack
  int n_stat;       // number of function static variables
  int self;
  ForLoop *flstack;		/* Ptr to for-loop list for a particular func. */
};

#define NFRAME 1024 /* Fixed frame stack size, old val 512 */
static Frame *frame;      /* Frame stack */
static Frame *fp;         /* Frame pointer */

static int NSTACK = 1024;   /* Fixed stack size, old val 1024 */
static Datum *stack;        /* the stack */
static Datum *stackp;       /* next free spot on stack */

static Program *program;    /* The current program for code generation */
static Inst *pc;            /* Program counter during execution */
static int progoff;         /* The current program offset, from beginning */

//static int debug_execute = 0;

extern char *rl_histfile;
extern int   rl_histsize;

/* dissassem.c */
extern void diss_assemble (Inst * p, int progoff);

/* scan.l */
extern char *curr_file_name;
extern int lineno;		/* scan.l */
extern int loff;		/* rlab.y */
extern void destroy_fstatic_tree (void);	/* rlab.y */

static void class_return (void);

Datum op_pull_datum;
extern Btree * get_class_publ_symtab();
extern void null_class_publ_symtab();
extern int inc_class_scope();
extern int dec_class_scope();
extern int get_class_scope();

/* Protos for code.c functions */
void push (Datum d);
Datum pop (void);
void extern_push (Datum d);
Datum extern_pop (void);

void load_local_var (Function * f, int nargs);
void load_stat_var (Function * f);
void print_and_assign (Datum d);
void print (Datum d);
void inc (void);
void dec (void);
Datum assign (Datum d1, Datum d2);
Datum assign_op (Datum d1, Datum d2, int op_idx);
void bltin (Ent * esp, int narg, int popf);
void function_call (void);
void function_call_1 (void);
void function_call_2 (void);
void function_call_self (void);
void userf (Ent * esp, int nargs, int self);
void user_class (Ent * esp, int nargs);

void list_create (void);
void list_member (void);
void list_assign (void);
void list_el_create (void);
Datum olist_assign (Datum d1, Datum d2);

/* print.c */
extern void destroy_file_list (void);

/* relation.c */
extern Datum eq (Datum d1, Datum d2);
extern Datum ne (Datum d1, Datum d2);
extern Datum le (Datum d1, Datum d2);
extern Datum lt (Datum d1, Datum d2);
extern Datum ge (Datum d1, Datum d2);
extern Datum gt (Datum d1, Datum d2);
extern Datum and (Datum d1, Datum d2);
extern Datum or (Datum d1, Datum d2);
extern Datum not (Datum d);

/* print_code.c */
extern void print_code (Inst * p);

/* op.c */
extern Datum addition_op (Datum d1, Datum d2, int reuse_d1_ent);
extern Datum subtraction_op (Datum d1, Datum d2);
extern Datum multiply_op (Datum d1, Datum d2);
extern Datum rdivide (Datum d1, Datum d2);
extern Datum el_rdivide (Datum d1, Datum d2);
extern Datum ldivide (Datum d1, Datum d2);
extern Datum el_ldivide (Datum d1, Datum d2);
extern Datum negate_op (Datum d);
extern Datum power_op (Datum d1, Datum d2);
extern Datum el_power (Datum d1, Datum d2);
extern Datum empty_create (void);

/* op_vector.c */
extern Datum vector_create_op (int n, Datum d1, Datum d2, Datum d3);
extern Datum vector_append (Datum d1, Datum d2);

/* opmat.c */
extern Datum matrix_create (Datum d);
extern Datum matrix_stack (Datum d1, Datum d2);
extern Datum matrix_sub_1 (Datum i1, Datum i2, Datum var);
extern Datum matrix_sub_2 (Datum i1, Datum var);
extern Datum matrix_sub_3 (Datum i1, Datum var);
extern Datum matrix_vector_sub (Datum d1, Datum d2);
extern Datum matrix_assign_1 (Datum var, Datum i1, Datum i2, Datum a);
extern Datum matrix_assign_2 (Datum var, Datum i1, Datum a);
extern Datum matrix_assign_3 (Datum var, Datum i2, Datum a);

extern Datum matrix_vector_assign (Datum m, Datum i1, Datum rhs);
extern Datum matrix_transpose (Datum d);
extern Datum matrix_el_transpose (Datum d);
extern Datum element_multiply (Datum d1, Datum d2);
extern Datum matrix_reshape_col (Datum d);

/* scan.l */
extern void scanner_cleanup (void);

/* **************************************************************
 * Miscellaneous functions to avoid using global data. Global data
 * is not that bad, but we try to avoid it here cause we might end
 * up using this file as part of a library.
 * ************************************************************** */

static int print_machine;	/* if TRUE print out machine contents */

void set_print_machine (int val)
{
  print_machine = val;
}

int
get_print_machine (void)
{
  return (print_machine);
}

static int line_nos;		/* if TRUE include line numbers in code */

void set_line_nos (int val)
{
  line_nos = val;
}

static int use_pager;		/* if TRUE use popen () */

void set_use_pager (int val)
{
  use_pager = val;
}

static char *pager;

void set_code_pager (char *value)
{
  pager = cpstr (value);
}

int get_progoff (void)
{
  return progoff;
}

int
set_progoff (void)
{
  progoff = program->progp - program->prog;
  return (progoff);
}


/* **************************************************************
 * Catch Ctrl-Cs
 * ************************************************************** */

static int interrupt = 0;

void intcatch (int tmp)
{
#ifdef __EMX__
  /* to re-enable signals under EMX */
  signal (SIGINT, SIG_ACK);
#endif

  signal (SIGINT, SIG_IGN);
  rerror ("user-generated interrupt");
}

void intcatch_wait (int tmp)
{

#ifdef __EMX__
  /* to re-enable signals under EMX */
  signal (SIGINT, SIG_ACK);
#endif

  /*
   * Ignore interrupts until we finish with this one.
   */

  signal (SIGINT, SIG_IGN);

  /*
   * Set the flag so we stop next time we enter the
   * machine switch.
   */

  interrupt = 1;

  /*
   * Print something out to re-assure the user
   * that something is indeed happening.
   */

  fprintf (stderr, "caught user-generated interrupt, just a moment...\n");
  fflush (stderr);
}

/* **************************************************************
 * Initialization... This stuff is done this way so that machines
 * with small static memory sizes (PCs and MACs) can
 * ************************************************************** */

/*
 * Initialize the Frame stack
 */

void
init_frame (void)
{
  frame = (Frame *) GC_MALLOC (sizeof (Frame) * NFRAME);
  if (frame == 0)
    rerror ("out of memory, save and quit!");
  fp = frame;
}

/*
 * Initialize the data stack.
 */

void
init_stack (void)
{
  stack = (Datum *) GC_MALLOC (sizeof (Datum) * NSTACK);
  if (stack == 0)
    rerror ("out of memory, save and quit!");
  stackp = stack;
}

void
init_flstack (void)
{
  frame->flstack = 0;
}

/* **************************************************************
 * Functions that create and manipulate program arrays.
 * ************************************************************** */
Program * program_Create (int n)
{
  Program *new;
  new = (Program *) GC_MALLOC (sizeof (Program));
  if (new == 0)
    rerror ("out of memory");

  new->ncode = n;
#ifdef HAVE_GC
  new->prog = (Inst *) GC_MALLOC (n * sizeof (Inst));
#else
  new->prog = (Inst *) CALLOC (n, sizeof (Inst));
#endif
  if (new->prog == 0)
    rerror ("out of memory");

  new->progp = new->prog;
  new->off = 0;

  return (new);
}

/*
 * Destroy a Program...
 */
void program_Destroy (Program * program)
{
  /* int i; */

  /*
   * Try and free() up the string-literals.
   */

  /* Don't use this piece of code until the
     offsets (jumps) for for, whiles, and if
     statements won't collide with the OP-CODES.
     for (i = 0; i < program->ncode; i++)
     {
     if (program->prog[i].op_code == OP_PUSH_STRING)
     {
     GC_FREE (program->prog[i + 1].ptr);
     }
     }
   */

  /*
   * Go ahead, and get rid of the rest...
   */

  program->ncode = 0;
  program->off = 0;
  program->progp = (Inst *) 0;
  GC_FREE (program->prog);
  GC_FREE (program);
}

/*
 * Grow the current program space.
 */

void
program_Grow (void)
{
  size_t new_size, old_size;
  void *result;

#ifdef HAVE_GC
  old_size = GC_size (program->prog);
  new_size = 3 * old_size;
  result = (Inst *) GC_MALLOC (new_size);
#else
  old_size = sizeof (Inst) * program->ncode;
  new_size = 3 * old_size;
  result = (Inst *) CALLOC (3 * (program->ncode), sizeof (Inst));
#endif

  if (result == 0)
    rerror ("out of memory, save and quit!");

  if (program->prog == 0)
    rerror ("terrible program internal error, save and quit");

  memcpy (result, program->prog, old_size);
  GC_FREE (program->prog);
  program->prog = result;

  /* Increase ncode size to keep up with prog[] */
  program->ncode *= 3;

  /* Reset progp using old offset */
  program->progp = program->prog + progoff;
}

/*
 * Set program to point at the supplied program.
 */

void
program_Set (Program * progptr)
{
  program = progptr;
}

Program *
program_Get (void)
{
  return (program);
}

/*
 * Get and restore, the program-counter (pc).
 * Add one so that when we restart, we are not
 * in the same place.
 */

Inst *
get_program_counter (void)
{
  return (pc);
}

void
set_program_counter (Inst * program_counter)
{
  pc = program_counter;
}

int get_function_scope (void) 
{
  if ((fp - frame) > 0)
    return (LOCAL);
  else
    return (GLOBAL);
}

List * get_function_arglist (void)
{
  if ((fp - frame) > 0)
  {
    return (function_GetArgPtr (ent_data (fp->esp)));
  }
  else
  {
    return (0);
  }
}

List * get_function_statlist (void)
{
  if ((fp - frame) > 0)
  {
    return (function_GetStatPtr (ent_data (fp->esp)));
  }
  else
  {
    return (0);
  }
}


List *
get_function_locallist (void)
{
  if ((fp - frame) > 0)
  {
    return (function_GetLocalPtr (ent_data (fp->esp)));
  }
  else
  {
    return (0);
  }
}

List *
get_function_globallist (void)
{
  if ((fp - frame) > 0)
  {
    return (function_GetGlobalPtr (ent_data (fp->esp)));
  }
  else
  {
    return (0);
  }
}

void
reset_frame_ptr (void)
{
  fp = frame;
}

Datum *
get_stackp (void)
{
  return (stackp);
}

void
reset_stack_ptr (void)
{
  stackp = stack;
}

void
set_stackp (Datum * new_stackp)
{
  stackp = new_stackp;
}

/* **************************************************************
 * Initialize the program for code generation.
 * ************************************************************** */
void initcode (void)
{
  /* Reset current program pointer */
  program->progp = program->prog;

  /* Reset parser/program offset */
  progoff = 0;

  if (line_nos)
  {
    /* Set the file name for the next section of code */
    (program->progp++)->op_code = OP_FILE_NAME;
    (program->progp++)->ptr = curr_file_name;
    (program->progp++)->op_code = OP_LINE_NO;
    (program->progp++)->op_code = lineno + loff;
    (program->progp)->op_code = STOP;
  }

  /* Compute the current code offset, for the parser */
  progoff = program->progp - program->prog;
}

/* **************************************************************
 * Push d (Datum) onto the stack. Update the stack pointer (stackp).
 * ************************************************************** */

void
push (Datum d)
{
  if (stackp >= &stack[NSTACK - 1])
  {
    rerror ("interpreter data stack overflow");
  }
  *stackp++ = d;
}

void
extern_push (Datum d)
{
  if (stackp >= &stack[NSTACK - 1])
  {
    rerror ("interpreter data stack overflow");
  }
  *stackp++ = d;
}

/* **************************************************************
 * Pop and return the top element from the stack.
 * ************************************************************** */

Datum
pop (void)
{
  Datum d;

  if (stackp <= stack)
  {
    rerror ("interpreter data stack underflow");
  }

  d = *--stackp;
  (*stackp).u.ptr = 0;

  return (d);
}

Datum
extern_pop (void)
{
  Datum d;

  if (stackp <= stack)
  {
    rerror ("interpreter data stack underflow");
  }

  d = *--stackp;
  (*stackp).u.ptr = 0;

  return (d);
}

/* **************************************************************
 * Clean up the data-stack after an error. Decrement any entities
 * that we encounter, and set all pointers to 0.
 * ************************************************************** */

void
datum_stack_clean (void)
{
  Datum d;

  while (stackp != stack)
  {
    d = *--stackp;
    if (d.type == ENTITY)
    {
      ent_Destroy (d.u.ptr);
      d.u.ptr = 0;
    }
    else if (d.type == VAR)
    {
      d.u.ptr = 0;
    }
  }
}

/* **************************************************************
 * Install one instruction or operand. Put an instruction into the
 * next free spot. Return an offset from the base.
 * Becasue the machine instructions are a union, we have a "code()"
 * function for each type of instruction we want to install.
 * ************************************************************** */
int code (int f)
{
  Inst *oprogp = program->progp;

  if (program->progp - program->prog >= program->ncode - 1)
  {
    program_Grow ();
    oprogp = program->progp;
  }

  (program->progp++)->op_code = f;
  progoff = program->progp - program->prog;

  /* Return the offset from prog, of the last instruction */
  return (oprogp - program->prog);
}

/* **************************************************************
 * Install a pointer to another op-code.
 * ************************************************************** */
int codep (void *f)
{
  Inst *oprogp = program->progp;

  if (program->progp - program->prog >= program->ncode - 1)
  {
    program_Grow ();
    oprogp = program->progp;
  }

  (program->progp++)->ptr = f;
  progoff = program->progp - program->prog;
  return (oprogp - program->prog);
}

int coded (double f)
{
  Inst *oprogp = program->progp;

  if (program->progp - program->prog >= program->ncode - 1)
  {
    program_Grow ();
    oprogp = program->progp;
  }

  (program->progp++)->d_val = f;
  progoff = program->progp - program->prog;
  return (oprogp - program->prog);
}

/* **************************************************************
 * Install op-codes (jump offsets) at special places in the
 * program array.
 * ************************************************************** */
void code_sp (int offset, int value)
{
  program->prog[offset].op_code = value;
}

int get_code_sp (int offset)
{
  return (program->prog[offset].op_code);
}

/* **************************************************************
 * function_call() serves as the front-end for both user-functions
 * and built-in functions.
 * ************************************************************** */

#undef THIS_SOLVER
#define THIS_SOLVER "function_call"
void function_call (void)
{
  int nargs;
  Datum fst;
  Ent *esp;
  ListNode *sp;
  FILE *rlab_stderr = (!RLAB_STDERR_DS) ? stderr : RLAB_STDERR_DS;

  /* Number of function arguments */
  nargs = (*pc).op_code;

//   dprintf("nargs=%i\n", nargs);

  /* Variable holding user/builtin function. */

  fst = pop ();

  sp = (ListNode *) (fst.u.ptr);

//   dprintf("sp=%p\n", sp);
//   dprintf("name=%s\n", sp->key);

  esp = var_ent (sp);

  /*
   * Verify that we have either user-function or built-in
   */

  if ((ent_type (esp) != U_FUNCTION) && (ent_type (esp) != BLTIN) && (ent_type (esp) != U_CLASS))
  {
    if (ent_type (esp) == UNDEF)
    {
      fprintf (rlab_stderr, THIS_SOLVER ": name = %s \n", sp->key);
      fprintf (rlab_stderr, THIS_SOLVER ": " RLAB_ERROR_FUNC_CANT_USE_UNDEF_VAR "\n");
      rerror (RLAB_ERROR_FUNC_UNDEFINED);
    }
    else
    {
      rerror (RLAB_ERROR_FUNC_INVALID);
    }
  }

  /*
   * Call bltin() if we have a built-in function
   */
  switch (ent_type (esp))
  {
    case BLTIN:
      // Call the built-in from here, all the
      // right info is on the stack
      bltin (esp, nargs, 0);
      break;

    case U_FUNCTION:
      // Call the user-function from here, all the
      // right info is on the stack
      userf (esp, nargs, 0);
      break;

    case U_CLASS:
      user_class (esp, nargs);
      break;
  }
  pc++;

//   dprintf("Out\n");
}

/*
 * Call a function from inside a list.
 */
#undef THIS_SOLVER
#define THIS_SOLVER "function_call_1"
void function_call_1 (void)
{
  int nargs;
  Ent *esp=0;
  ListNode *sp=0;
  FILE *rlab_stderr = (!RLAB_STDERR_DS) ? stderr : RLAB_STDERR_DS;

//   ListNode *spm1=0, *spp1=0;

  // Number of function arguments
  nargs = (*pc).op_code;

  // Entity holding user/builtin function
  sp = (ListNode *) ((stackp - nargs - 1)->u.ptr);
  // Either this is a list node or undefined (E. Plischke):
  // for function lists this does not produce correct results!
  // removed since rlab-2.2.11.6
  //   if ((ent_type (sp) == UNDEF))
  //   {
  //     rerror ("list undefined");
  //   }

  //
  // Verify that we have either user-function or built-in
  //
  esp = var_ent (sp);

  if (!esp)
  {
    fprintf (rlab_stderr, THIS_SOLVER ": " RLAB_ERROR_FUNC_CANT_USE_UNDEF_VAR "\n");
    rerror (RLAB_ERROR_FUNC_UNDEFINED);
  }
  if (ent_type (esp) == UNDEF)
  {
    fprintf (rlab_stderr, THIS_SOLVER ": " RLAB_ERROR_FUNC_CANT_USE_UNDEF_VAR "\n");
    rerror (RLAB_ERROR_FUNC_UNDEFINED);
  }
  if ((ent_type (esp) != U_FUNCTION) && (ent_type (esp) != BLTIN))
  {
    rerror (RLAB_ERROR_FUNC_INVALID);
  }

  //
  // Call bltin() if we have a built-in function
  // 
  if (ent_type (esp) == BLTIN)
  {
    // Call the built-in from here, all the
    //  right info is on the stack
    bltin (esp, nargs, 1);
    pc++;
    return;
  }
  else
  {
    // Call the user-function from here, all the
    // right info is on the stack
    userf (esp, nargs, 1);
    pc++;
    return;
  }
}

#undef THIS_SOLVER
#define THIS_SOLVER "function_call_2"
void function_call_2 (void)
{
  int nargs;
  Ent *esp;
  FILE *rlab_stderr = (!RLAB_STDERR_DS) ? stderr : RLAB_STDERR_DS;

  /* Number of function arguments */
  nargs = (*pc).op_code;

  if ((stackp - nargs - 1)->type != ENTITY)
  {
    fprintf (rlab_stderr, THIS_SOLVER ": " RLAB_ERROR_FUNC_CANT_USE_CONST "\n");
    rerror (RLAB_ERROR_FUNC_UNDEFINED);
  }

  /* Variable holding user/builtin function. */
  esp = (Ent *) ((stackp - nargs - 1)->u.ptr);

  /*
  * Verify that we have either user-function or built-in
  */
  if ((ent_type (esp) != U_FUNCTION) && (ent_type (esp) != BLTIN))
  {
    if (ent_type (esp) == UNDEF)
    {
      fprintf (rlab_stderr, THIS_SOLVER ": " RLAB_ERROR_FUNC_CANT_USE_UNDEF_VAR "\n");
      rerror (RLAB_ERROR_FUNC_UNDEFINED);
    }
    else
    {
      rerror (RLAB_ERROR_FUNC_INVALID);
    }
  }

  /*
  * Call bltin() if we have a built-in function
  */

  if (ent_type (esp) == BLTIN)
  {
    /*
    * Call the built-in from here, all the
    * right info is on the stack
    */

    bltin (esp, nargs, 0);
    pc++;
    return;
  }
  else
  {
    /*
    * Call the user-function from here, all the
    * right info is on the stack
    */

    userf (esp, nargs, 0);
    pc++;
    return;
  }
}

/*
 * Do the work of re-directing execution into the
 * function stack.
 *
 * Function arguments are kept on the stack, to support
 * recursion. Reference counts keep the arguments alive
 * as long as necessary.
 */

extern char * set_priv_class_name(char * );
extern char * set_classdef_filename(char * );

#undef THIS_SOLVER
#define THIS_SOLVER "user_class"
void user_class (Ent * esp, int nargs)
{
  Function *f = ent_data (esp);

  // prepare btree in static_tree for collecting and lookup of class-private variables
  set_priv_class_name(f->name);

  // remind parser what was the filename in which classdef was defined.
  // this file is used for file-static variables
  set_classdef_filename(f->classdef_filename);

  int i;
  Datum d;
  Ent *etmp;
  ListNode *lnode;
  List *fargs;

  //
  // Check frame stack for overflow, then inc, and load.
  // 
  if (fp >= &frame[NFRAME - 1])
  {
    rerror ("class call nested too deeply");
  }

  fp++;
  fp->nargs = nargs;
  fp->sp = 0;     /* Just for now. */
  fp->esp = esp;

  //
  // Check the # of args on the stack vs. function requ.
  // If there are fewer args present than required, push UNDEFs
  // onto the stack
  //
//   dprintf("args = %p\n", function_GetArgPtr (f));
//   dprintf("n_args = %i\n", function_GetNargs (f));
  if (fp->nargs < function_GetNargs (f))
  {
    for (i = 0; i < function_GetNargs (f) - fp->nargs; i++)
    {
      /* Create an UNDEF entity. */
      etmp = ent_Create ();
      ent_SetType (etmp, UNDEF);

      /* Push the UNDEF entity on the stack. */
      d.u.ptr = etmp;
      d.type = ENTITY;
      push (d);

      /* Bump the ref-count. */
      ent_IncRef (etmp);
    }
    fp->nargs = function_GetNargs (ent_data (esp));
  }
  else if (fp->nargs > function_GetNargs (f))
  {
    /* Error and return */
    fp->sp = 0;
    fp->esp = 0;
    fp->nargs = 0;
    --fp;
    rerror ("too many arguments to function");
  }
  fp->retpc = pc;          /* return pointer */

  /*
  * Go through the function args.
  * Make sure each argument is a variable.
  * Copy any existing variable, so they don't get dumped on.
  * Each arg must be a var, cause it might get assigned
  * to during the course of a function.
  */
  fp->argn = stackp - 1; // arguments start from this
  fargs = function_GetArgPtr (f);
  lnode = list_GetLastNode (fargs);

  for (i = 1; i <= fp->nargs; i++)
  {
    if (fp->argn[i - fp->nargs].type == VAR)
    {
      Datum dtmp = fp->argn[i - fp->nargs];
      ListNode *var = listNode_Create ();
      Ent *ent = var_ent (dtmp.u.ptr);
      Ent *lent = 0;
      listNode_SetOwned (var);

      /* Check for a variables list-entity. Copy if necessary. */
      if ((lent = var_listent (dtmp.u.ptr)))
      {
        listNode_AttachListEnt (var, lent);
        ent_IncRef (lent);
      }

      listNode_AttachEnt (var, ent);
      listNode_SetKey (var, var_key (lnode));

      ent_IncRef (ent);

      dtmp.type = VAR;
      dtmp.u.ptr = var;

      fp->argn[i - fp->nargs] = dtmp;
    }
    else if (fp->argn[i - fp->nargs].type == ENTITY)
    {
      Datum dtmp = fp->argn[i - fp->nargs];
      ListNode *var = listNode_Create ();
      Ent *ent = dtmp.u.ptr;

      listNode_SetOwned (var);
      listNode_AttachEnt (var, ent);
      listNode_SetKey (var, var_key (lnode));
      /* ent_IncRef (ent); */
      dtmp.type = VAR;
      dtmp.u.ptr = var;

      fp->argn[i - fp->nargs] = dtmp;
    }
    else if (fp->argn[i - fp->nargs].type == CONSTANT)
    {
      Datum dtmp = fp->argn[i - fp->nargs];
      Ent *ent = convert_datum_to_ent (dtmp);
      ListNode *var = listNode_Create ();
      listNode_SetOwned (var);

      listNode_AttachEnt (var, ent);
      listNode_SetKey (var, var_key (lnode));
      dtmp.type = VAR;
      dtmp.u.ptr = var;
      ent_IncRef (ent);

      fp->argn[i - fp->nargs] = dtmp;
    }
    else if (fp->argn[i - fp->nargs].type == iCONSTANT)
    {
      Datum dtmp = fp->argn[i - fp->nargs];
      Ent *ent = convert_datum_to_ent (dtmp);
      ListNode *var = listNode_Create ();
      listNode_SetOwned (var);

      listNode_AttachEnt (var, ent);
      listNode_SetKey (var, var_key (lnode));
      dtmp.type = VAR;
      dtmp.u.ptr = var;
      ent_IncRef (ent);

      fp->argn[i - fp->nargs] = dtmp;
    }
    lnode = listNode_GetPrevNode (lnode);
  }

  fp->n_stat = 0; // not used for user_class
  fp->localn = 0;
  fp->n_local = 0;

  // will be used for collecting the names of public variables
  // initialize list for collecting return from the class call
  /* Initialize the ptr for for-loop tracking. */
  fp->flstack = 0;
  fp->self = 0;

 /*
  * Now execute the classdef function
 */
  execute ((Inst *) function_GetCodePtr (ent_data (esp)));

  pc = fp->retpc;
  --fp;

  // clean-up
}

/*
 * Do the work of re-directing execution into the
 * function stack.
 *
 * Function arguments are kept on the stack, to support
 * recursion. Reference counts keep the arguments alive
 * as long as necessary.
 */
#undef THIS_SOLVER
#define THIS_SOLVER "userf"
void userf (Ent * esp, int nargs, int self)
{
  int i;
  Datum d;
  Ent *etmp;
  ListNode *lnode;
  List *fargs;

  /*
   * Check frame stack for overflow, then inc, and load.
   */

  if (fp >= &frame[NFRAME - 1])
  {
    rerror ("function call nested too deeply");
  }

  fp++;
  fp->nargs = nargs;
  fp->sp = 0;			/* Just for now. */
  fp->esp = esp;

  /*
   * Check the # of args on the stack vs. function requ.
   * If there are fewer args present than required, push UNDEFs
   * onto the stack
   */

  if (fp->nargs < function_GetNargs (ent_data (esp)))
  {
    for (i = 0; i < function_GetNargs (ent_data (esp)) - fp->nargs; i++)
    {
      /* Create an UNDEF entity. */
      etmp = ent_Create ();
      ent_SetType (etmp, UNDEF);

      /* Push the UNDEF entity on the stack. */
      d.u.ptr = etmp;
      d.type = ENTITY;
      push (d);

      /* Bump the ref-count. */
      ent_IncRef (etmp);
    }
    fp->nargs = function_GetNargs (ent_data (esp));
  }
  else if (fp->nargs > function_GetNargs (ent_data (esp)))
  {
    /* Error and return */
    fp->sp = 0;
    fp->esp = 0;
    fp->nargs = 0;
    fp--;
    rerror ("too many arguments to function");
  }

  //printf(THIS_FILE ": " THIS_SOLVER ": 1\n");


  fp->retpc = pc;		       /* return pointer */

  /*
   * Go through the function args.
   * Make sure each argument is a variable.
   * Copy any existing variable, so they don't get dumped on.
   * Each arg must be a var, cause it might get assigned
   * to during the course of a function.
   */
  fp->argn = stackp - 1; // arguments start from this 
  fargs = function_GetArgPtr (ent_data (esp));
  lnode = list_GetLastNode (fargs);
  for (i = 1; i <= fp->nargs; i++)
  {
    if (fp->argn[i - fp->nargs].type == VAR)
    {
      Datum dtmp = fp->argn[i - fp->nargs];
      ListNode *var = listNode_Create ();
      Ent *ent = var_ent (dtmp.u.ptr);
      Ent *lent = 0;
      listNode_SetOwned (var);

      /* Check for a variables list-entity. Copy if necessary. */
      if ((lent = var_listent (dtmp.u.ptr)))
      {
        listNode_AttachListEnt (var, lent);
        ent_IncRef (lent);
      }

      listNode_AttachEnt (var, ent);
      listNode_SetKey (var, var_key (lnode));

      ent_IncRef (ent);

      dtmp.type = VAR;
      dtmp.u.ptr = var;

      fp->argn[i - fp->nargs] = dtmp;
    }
    else if (fp->argn[i - fp->nargs].type == ENTITY)
    {
      Datum dtmp = fp->argn[i - fp->nargs];
      ListNode *var = listNode_Create ();
      Ent *ent = dtmp.u.ptr;

      listNode_SetOwned (var);
      listNode_AttachEnt (var, ent);
      listNode_SetKey (var, var_key (lnode));
      /* ent_IncRef (ent); */
      dtmp.type = VAR;
      dtmp.u.ptr = var;

      fp->argn[i - fp->nargs] = dtmp;
    }
    else if (fp->argn[i - fp->nargs].type == CONSTANT)
    {
      Datum dtmp = fp->argn[i - fp->nargs];
      Ent *ent = convert_datum_to_ent (dtmp);
      ListNode *var = listNode_Create ();
      listNode_SetOwned (var);

      listNode_AttachEnt (var, ent);
      listNode_SetKey (var, var_key (lnode));
      dtmp.type = VAR;
      dtmp.u.ptr = var;
      ent_IncRef (ent);

      fp->argn[i - fp->nargs] = dtmp;
    }
    else if (fp->argn[i - fp->nargs].type == iCONSTANT)
    {
      Datum dtmp = fp->argn[i - fp->nargs];
      Ent *ent = convert_datum_to_ent (dtmp);
      ListNode *var = listNode_Create ();
      listNode_SetOwned (var);

      listNode_AttachEnt (var, ent);
      listNode_SetKey (var, var_key (lnode));
      dtmp.type = VAR;
      dtmp.u.ptr = var;
      ent_IncRef (ent);

      fp->argn[i - fp->nargs] = dtmp;
    }
    lnode = listNode_GetPrevNode (lnode);
  }

  fp->n_local = 0;
  fp->n_stat = 0;
  fp->localn = 0;

 /*
  * Set local variable data.
  */
  if (function_HasLocalVar (ent_data (esp)))
  {
    fp->n_local += function_GetNlocal (ent_data (esp));
    load_local_var (ent_data (esp), nargs);
    fp->localn = stackp - 1;
  }

 /*
  * Check for func static variables, and attach them at the end of the local list
  */
  if (function_HasStatVar (ent_data (esp)))
  {
    load_stat_var (ent_data (esp));
    fp->n_stat = function_GetNstat (ent_data (esp));
    fp->statn  = stackp - 1;
  }

  /* Initialize the ptr for for-loop tracking. */
  fp->flstack = 0;

  fp->self = self;

  /*
   * Now execute the user-function.
   */
  //printf(THIS_FILE ": " THIS_SOLVER ": 3\n");

  execute ((Inst *) function_GetCodePtr (ent_data (esp)));

  //printf(THIS_FILE ": " THIS_SOLVER ": OUT\n");

  pc = fp->retpc;
  --fp;
}

/*
 * Do the work of re-directing execution into the
 * function stack.
 *
 * Function arguments are kept on the stack, to support
 * recursion. Reference counts keep the arguments alive
 * as long as necessary.
 */
/*
 * Special, a function is calling itself.
 */

void
function_call_self (void)
{
  int nargs;

  /* Number of args on the stack */
  nargs = (*pc).op_code;

  userf (fp->esp, nargs, 0);
  pc++;
}

/* **************************************************************
 * Load a functions local variables onto the stack. All local
 * variables are initialized to UNDEFINED entities, in order to be
 * consistent with global, and argument variable rules.
 * ************************************************************** */
#undef  THIS_SOLVER
#define THIS_SOLVER "load_local_var"
void load_local_var (Function * f, int nargs)
{
  int i, n_var=function_GetNlocal (f);
  char *name;
  Datum d;
  Ent *ent;
  List *local_list=function_GetLocalPtr (f);
  ListNode *lnode;

  lnode = list_GetLastNode (local_list);
  for (i = 0; i < n_var; i++)
  {
    /*
     * Create an UNDEF variable to push on stack.
     * We must push a variable, because it may be assigned
     * to. Inc the ref-count to keep it around long enough.
     */

    d.u.ptr = listNode_Create ();
    listNode_SetOwned ((ListNode *) (d.u.ptr));
    d.type = VAR;

    ent = ent_Create ();
    listNode_AttachEnt (d.u.ptr, ent);

    name = var_key (lnode);
    if (!(strcmp ("nargs", name)))
    {
      // Set local variable nargs (number of arguments)
#if 0
      ent_SetType (ent, DOUBLE);
      ent_double (ent) = (double) nargs;
      ent_IncRef (ent);		/* Owned by the function. */
#else
      // this is done because all commands are defined for
      // MDRs and similar objects, and not doubles
      ent_data (ent) = mdr_CreateScalar(nargs);
      ent_type (ent) = MATRIX_DENSE_REAL;
#endif
    }
    else
    {
      ent_SetType (ent, UNDEF);
    }
    ent_IncRef (ent);   /* Owned by the function. */

    listNode_SetKey (d.u.ptr, name);
    push (d);
    lnode = listNode_GetPrevNode (lnode);
  }
}

#undef  THIS_SOLVER
#define THIS_SOLVER "load_stat_var"
void load_stat_var (Function * f)
{
  int i, n_var=function_GetNstat (f);
  char *name;
  Datum d;
  List *stat_list=function_GetStatPtr (f);
  ListNode *lnode=list_GetLastNode (stat_list);

  for (i = 0; i < n_var; i++)
  {
    name = var_key (lnode);
    if (!(strcmp ("nargs", name)))
    {
      printf(THIS_FILE ": " THIS_SOLVER
          ":WARNING: Cannot have function-static variable named 'nargs'. The name is reserved!\n");
      continue;
    }

   /*
    * Create an UNDEF variable to push on stack.
    * We must push a variable, because it may be assigned
    * to. Inc the ref-count to keep it around long enough.
    */
    d.u.ptr = listNode_Create ();
    listNode_SetKey (d.u.ptr, name);
    listNode_SetOwned ((ListNode *) (d.u.ptr));
    d.type = VAR;

    Ent * ent = var_ent(lnode);
    if (!ent)
    {
      var_ent(lnode) = ent_Create ();
      ent_SetType (var_ent(lnode), UNDEF);
      ent = var_ent(lnode);
    }

    listNode_AttachEnt (d.u.ptr, ent);
    ent_IncRef (ent);   /* Owned by the function. */
    push (d);

    lnode = listNode_GetPrevNode (lnode);
  }
}


// **************************************************************
//
//  rlab3: Return from an RLaB User-Class
//
// **************************************************************
#undef  THIS_SOLVER
#define THIS_SOLVER "class_return"
extern Btree * just_get_class_stat_symtab();
extern void    null_class_stat_symtab();
extern Btree * get_class_publ_symtab();
extern void    null_class_publ_symtab();
void class_return (void)
{
  int i;
  Ent * etmp=0;

  pop(); // this is the same as btree below!

  //
  // class-member public variables are in a btree.
  //
  Datum d;
  ListNode *nvar=listNode_Create();
  d.type = VAR;
  d.u.ptr = nvar;
  var_ent(nvar) = ent_Assign_Rlab_BTREE(get_class_publ_symtab());
  ent_Ref(var_ent(nvar)) = 1;
  null_class_publ_symtab();

  //
  // class-member static (private) variables are in a btree, as well
  //
  if (just_get_class_stat_symtab())
  {
    Ent *etmp = ent_Assign_Rlab_BTREE(just_get_class_stat_symtab());
    ent_Ref(etmp) = 1;
    var_listent(nvar) = etmp;
    null_class_stat_symtab();
  }
  else
    var_listent(nvar) = 0;

  //
  // Now pop the args off of the stack.
  // Destroy each argument. We must check for
  // variable or entity before we destroy.
  //
  for (i = 0; i < fp->nargs; i++)
  {
    Datum tmp = pop ();
    switch(tmp.type)
    {
      case ENTITY:
        ent_Destroy(tmp.u.ptr);
        break;

      case VAR:
        etmp = var_ent (tmp.u.ptr);
        ent_Destroy (etmp);    /* For losing the func-var. */
        break;
    }
    if (var_listent (tmp.u.ptr))
      ent_Destroy (var_listent (tmp.u.ptr));
  }

  push (d); // put d back on the stack
}




/* **************************************************************
 * Return from an RLaB User-Function
 * ************************************************************** */

void ret_from_func (void);

#undef  THIS_SOLVER
#define THIS_SOLVER "function_return"
void function_return (void)
{
//   dprintf("In\n");

  Datum d;
  Ent *ent;
  ListNode *nvar;

  /*
   * Pop the function return value from the stack.
   * ret_from_func () needs to clean the stack up.
   */

  d = pop ();

  /*
   * Now, fix up the return thing, so that it is
   * an entity.
   */

  if (d.type == VAR)
  {
    ent = var_ent (d.u.ptr);

    nvar = listNode_Create ();
    listNode_AttachEnt (nvar, ent);

    if (var_listent (d.u.ptr))
    {
      Ent *lent = ent_Create ();
      listNode_AttachListEnt (nvar, lent);
      ent_data (lent) = btree_Copy (ent_data (var_listent (d.u.ptr)));
      ent_type (lent) = BTREE;
      ent_IncRef (lent);
    }

    /* Don't set this variable as being owned by anything... its not! */
    ent_IncRef (ent);

    d.type = VAR;
    d.u.ptr = nvar;
  }

  /*
   * Now call ret_from_func to clean up the stack.
   */

  ret_from_func ();

  push (d);
//   dprintf("Out\n");
}

/* **************************************************************
 * Take care of removing arguments and local variables from the
 * stack now that the function is returning.
 * ************************************************************** */
extern void class_print(Ent *);
#undef  THIS_SOLVER
#define THIS_SOLVER "ret_from_func"
void ret_from_func (void)
{
  int i;
  Datum tmp;
  Ent *ent, *etmp, *listent;

  // reatach function static variables in case they have changed
  List *stat = function_GetStatPtr (ent_data (fp->esp));
  ListNode * lnode = list_GetFirstNode(stat);
  for (i=0; i<fp->n_stat; i++)
  {
    /* Pop the variable. */
    tmp = pop ();

    // this should be reattached to the original list of static variables
    // Entity ptr for function */
    etmp = var_ent (tmp.u.ptr);
    listent = var_listent (tmp.u.ptr);
    if (ent_type(etmp)!=UNDEF)
    {
      ent_Clean(var_ent(lnode));
      var_ent(lnode) = etmp;
    }
    if (listent)
    {
      if (var_listent(lnode))
      {
        ent_Clean(var_listent(lnode));
        var_listent(lnode) = listent;
      }
    }

    var_ent (tmp.u.ptr) = 0;
    var_listent (tmp.u.ptr) = 0;
    /* Destroy the variable structure. */
    listNode_Destroy (tmp.u.ptr);

    lnode = listNode_GetNodeAhead(lnode);
  }

 /*
  * Pop all local-vars off of the stack.
  * Destroy each local variable.
  */
  for (i=0; i<fp->n_local; i++)
  {
    /* Pop the variable. */
    tmp = pop ();

    etmp = var_ent (tmp.u.ptr);
    listent = var_listent (tmp.u.ptr);

    /* Check for, and free the listent. */
    if (listent)
      ent_Destroy (listent);

    /* Destroy the entity. */
    ent_Destroy (etmp);

    /* Destroy the variable structure. */
    listNode_Destroy (tmp.u.ptr);
  }

  /*
   * Now pop the args off of the stack.
   * Destroy each argument. We must check for
   * variable or entity before we destroy.
   */

  for (i = 0; i < fp->nargs; i++)
  {
    tmp = pop ();
    ent = var_ent (tmp.u.ptr);

    /* Check for, and free the listent. */
    if (var_listent (tmp.u.ptr))
      ent_Destroy (var_listent (tmp.u.ptr));

    listNode_Destroy (tmp.u.ptr);
    ent_Destroy (ent);		/* For losing the func-var. */
  }

  if (fp->self)
  {
    tmp = pop ();
    if (tmp.type == ENTITY)
    {
      ent_DecRef (tmp.u.ptr);
    }
    else if (tmp.type != VAR)
      rerror ("terrible function error");
  }

  /* Don't clean up the for-loop stack, the OP_FOR...
   * do that themselves.
   */
}

/* **************************************************************
 * Default way to return from a function if there is no explicit
 * return statement. In this case we clean the args and local
 * variables off the stack. Push a 0 scalar onto the stack.
 * ************************************************************** */

void
function_default_return (void)
{
  int i;
  Datum tmp, d;
  Ent *etmp, *ent;

  /*
   * Pop all local-vars off of the stack.
   * Destroy each local variable.
   */

  for (i = 0; i < fp->n_local; i++)
  {
    tmp = pop ();
    etmp = var_ent (tmp.u.ptr);

    /* Check for, and free the listent. */
    if (var_listent (tmp.u.ptr))
      ent_Destroy (var_listent (tmp.u.ptr));

    listNode_Destroy (tmp.u.ptr);
    ent_Destroy (etmp);
  }

  /*
   * Now pop the args off of the stack.
   * Destroy each argument. We must check for
   * variable or entity before we destroy.
   */

  for (i = 0; i < fp->nargs; i++)
  {
    tmp = pop ();
    ent = var_ent (tmp.u.ptr);

    /* Check for, and free the listent. */
    if (var_listent (tmp.u.ptr))
      ent_Destroy (var_listent (tmp.u.ptr));

    listNode_Destroy (tmp.u.ptr);
    ent_Destroy (ent);		/* For losing the func-var. */
  }

  if (fp->self)
  {
    tmp = pop ();
    if (tmp.type == ENTITY)
    {
      ent_DecRef (tmp.u.ptr);
    }
    else if (tmp.type != VAR)
      rerror ("terrible function error");
  }

  /* Don't clean up the for-loop stack, the OP_FOR...
   * do that themselves.
   */

  /*
   * Create a 0 scalar and push it on the stack so print()
   * will have something to do.
   */

  d.type = CONSTANT;
  d.u.val = 0.0;
  push (d);

  //printf(THIS_FILE ": " THIS_SOLVER ": OUT\n");
}

/* **************************************************************
 * Run the machine. p is the pointer to the machine code where
 * we will start execution.
 * ************************************************************** */

/*
 * Global var to handle return value from eval().
 */

static Datum eval_ret;

Ent * get_eval_ret (void)
{
  Ent *ent=0;
  ent = convert_datum_to_ent (eval_ret);
  return (ent);
}


#undef  THIS_SOLVER
#define THIS_SOLVER "execute"
static int exit_file=0;
static char *file2exit=0;

void execute (Inst * p)
{
  int n, offset, op_code;
  Datum d, d1, d2, d3, new;

  pc = p;

//   printf(THIS_FILE ": " THIS_SOLVER ": Entering 'execute'\n");

  while (1)
  {
    if (interrupt)
    {
      signal (SIGINT, intcatch_wait);
      interrupt = 0;
      rerror ("user-generated interrupt processed");
    }

    if (exit_file)
    {
      if (!curr_file_name)
      {
        GC_FREE (file2exit);
        file2exit = 0;
        exit_file = 0;
      }
      else if (strcmp(file2exit,curr_file_name))
      {
        GC_FREE (file2exit);
        file2exit = 0;
        exit_file = 0;
      }
      else
        return;
    }

    op_code = (*pc++).op_code;
//     if(debug_execute)
//     {
//       dprintf("op_code = %i\n", op_code);
//     }

    //switch ((*pc++).op_code)
    switch (op_code)
    {
      case OP_JMP_ENDOFFILE:
        file2exit = cpstr(curr_file_name);
        exit_file = 1;
      case STOP:
      case OP_ENDOFFILE:
        /* pc--; */
        return;
        break;

      case OP_PUSH_VAR:
      {
        new.u.ptr = (ListNode *) ((*pc++).ptr);
        new.type = VAR;
        *stackp++ = new;
        break;
      }

      case OP_PUSH_ARG:
        offset = (*pc++).op_code;
        new = fp->argn[offset - fp->nargs];
        *stackp++ = new;
        break;

      case OP_PUSH_STATIC_VAR:
        offset = (*pc++).op_code;
        new = fp->statn[offset - fp->n_stat];
        *stackp++ = new;
        break;

      case OP_PUSH_LOCAL_VAR:
        offset = (*pc++).op_code;
        new = fp->localn[offset - fp->n_local];
        *stackp++ = new;
        break;

      case OP_ADD:
        d2 = *--stackp;
        d1 = *--stackp;;
        new = addition_op (d1, d2,0);
        *stackp++ = new;
        break;

      case OP_SUB:
        d2 = *--stackp;
        d1 = *--stackp;;
        new = subtraction_op (d1, d2);
        *stackp++ = new;
        break;

      case OP_MUL:
        d2 = *--stackp;
        d1 = *--stackp;;
        new = multiply_op (d1, d2);
        *stackp++ = new;
        break;

      case OP_DIV:
        d2 = *--stackp;
        d1 = *--stackp;;
        new = rdivide (d1, d2);
        *stackp++ = new;
        break;

      case OP_LDIV:
        d2 = *--stackp;
        d1 = *--stackp;;
        new = ldivide (d1, d2);
        *stackp++ = new;
        break;

      case OP_NEGATE:
        d = pop ();
        new = negate_op (d);
        push (new);
        break;

      case OP_POWER:
        d2 = *--stackp;
        d1 = *--stackp;;
        new = power_op (d1, d2);
        *stackp++ = new;
        break;

      case OP_ASSIGN:
        d1 = *--stackp;
        d2 = *--stackp;;
        new = assign (d2, d1);
        *stackp++ = new;
        break;

      case OP_FOR_LOOP_I:
      {
        Datum vec, var;
        Ent *evec;
        ForLoop *fl = (ForLoop *) GC_MALLOC (sizeof (ForLoop));
        if (fl == 0)
          rerror ("out of memory");

        /* Get the vector and the var from the stack */
        vec = pop ();
        var = pop ();

        /* Get the jump distance from the program. */
        fl->jmp = (int) (*pc).op_code;
        fl->jmp_offset = 0;

        /* Keep track of the for-loop variable (i). */
        fl->var = (ListNode *) var.u.ptr;

        /* Copy the vector into flstack */
        evec = convert_datum_to_ent (vec);
        fl->ent = evec;
        ent_IncRef (evec);  /* The for-loop now owns a copy. */

        /*
        * Set the cnt part of the ForLoop to 0. We will use it to
        * mark our location in the vector as we step through.
        */
        fl->cnt = 1;

        //
        // Set the size of the loop-vector.
        //
        fl->n = class_size (fl->ent);
        if (fl->n != 0)
        {
          /* Initialize the value of the loop-vector. */
          Ent *ed = class_forloop_value (fl->ent, fl->cnt++);
          if (ed)
          {
            ent_Destroy (var_ent (fl->var));
            var_ent (fl->var) = ed;
            ent_IncRef (var_ent (fl->var));
          }
          else
            rerror (RLAB_ERROR_FORLOOP_CANNOT_ITERATE);

          /* Add our new loop struct to FRONT of the list */
          fl->next = fp->flstack;
          fp->flstack = fl;

          /* Continue loop execution. */
          pc++;
        }
        else
        {
          /*
           * Loop-vector has zero length, jump to end.
           * Here we need to set the loop-variable (i) to
           * UNDEF, and move on.
           */

          ent_Destroy (var_ent (fl->var));
          listNode_AttachEnt (fl->var, ent_Create ());
          ent_SetType (var_ent (fl->var), UNDEF);
          ent_IncRef (var_ent (fl->var));

          /*
           * Jump to end and let OP_FOR_LOOP cleanup.
           * Push something on the list so that for-loop-done
           * doesn't create trouble.
           */

          fl->next = fp->flstack;
          fp->flstack = fl;
          pc = pc + fl->jmp;
        }
        break;
      }
      case OP_FOR_LOOP_DONE:
      {
        ForLoop *ftmp;

        /* Destroy the loop vector. */
        ent_Destroy (fp->flstack->ent);

        /* Pop the front of flstack, and cleanup. */
        if (fp->flstack == 0)
          rerror ("Terrible for-loop error");

        /* Pop the forloop struct off the front of the list. */
        ftmp = fp->flstack;
        fp->flstack = fp->flstack->next;
        GC_FREE (ftmp);

        break;
      }
      case OP_FOR_LOOP_DONE_RET:
      {
        ForLoop *ftmp;

        /* Destroy the loop vector. */
        ent_Destroy (fp->flstack->ent);

        /* Pop the front flstack, and cleanup. */
        if (fp->flstack == 0)
          rerror ("Terrible for-loop error");

        ftmp = fp->flstack;
        fp->flstack = fp->flstack->next;
        GC_FREE (ftmp);

        return;
        break;
      }
      case OP_EL_MUL:
        d2 = *--stackp;
        d1 = *--stackp;;
        new = element_multiply (d1, d2);
        *stackp++ = new;
        break;

      case OP_EL_DIV:
        d2 = *--stackp;
        d1 = *--stackp;;
        new = el_rdivide (d1, d2);
        *stackp++ = new;
        break;

      case OP_EL_LDIV:
        d2 = pop ();
        d1 = pop ();
        new = el_ldivide (d1, d2);
        push (new);
        break;

      case OP_EL_POWER:
        d2 = pop ();
        d1 = pop ();
        new = el_power (d1, d2);
        push (new);
        break;

      case OP_PUSH_CONSTANT:
        new.type = CONSTANT;
        new.u.val = (*pc++).d_val;
        push (new);
        break;

      case OP_PUSH_INTEGER:
        new.type = ENTITY;
        new.u.ptr = ent_Create_Rlab_Int((int)(*pc++).d_val);
        ent_IncRef(new.u.ptr);
        push (new);
        break;

      case OP_PUSH_iCONSTANT:
        new.type = iCONSTANT;
        new.u.val = (*pc++).d_val;
        push (new);
        break;

      case OP_PULL_DATUM:
        op_pull_datum = pop ();
        if (op_pull_datum.type == ENTITY)
        {
          Ent *ent = (Ent *) (op_pull_datum.u.ptr);
          ent_DecRef (ent);
        }
        else if(op_pull_datum.type == VAR)
        {
          Ent *ent = var_ent(op_pull_datum.u.ptr);
          ent_Clean (ent);
        }
        break;

      case OP_PUSH_DATUM:
        push(op_pull_datum);
        break;

      case OP_PRINT:
        d = pop ();
        print (d);
        break;

      case OP_PRINT_ASSIGN:
        d = pop ();
        print_and_assign (d);
        break;

      case OP_GT:
        d2 = pop ();
        d1 = pop ();
        new = gt (d1, d2);
        push (new);
        break;

      case OP_LT:
        d2 = pop ();
        d1 = pop ();
        new = lt (d1, d2);
        push (new);
        break;

      case OP_EQ:
        d2 = pop ();
        d1 = pop ();
        new = eq (d1, d2);
        push (new);
        break;

      case OP_GE:
        d2 = pop ();
        d1 = pop ();
        new = ge (d1, d2);
        push (new);
        break;

      case OP_LE:
        d2 = pop ();
        d1 = pop ();
        new = le (d1, d2);
        push (new);
        break;

      case OP_NE:
        d2 = pop ();
        d1 = pop ();
        new = ne (d1, d2);
        push (new);
        break;

      case OP_AND:
        d2 = pop ();
        d1 = pop ();
        new = and (d1, d2);
        push (new);
        break;

      case OP_OR:
        d2 = pop ();
        d1 = pop ();
        new = or (d1, d2);
        push (new);
        break;

      case OP_ADDTO:
        d1 = *--stackp;
        d2 = *--stackp;;
        new = assign_op (d2, d1, 1);
        *stackp++ = new;
        break;

      case OP_SUBFROM:
        d1 = *--stackp;
        d2 = *--stackp;;
        new = assign_op (d2, d1, 2);
        *stackp++ = new;
        break;

      case OP_EL_MUL_BY:
        d1 = *--stackp;
        d2 = *--stackp;;
        new = assign_op (d2, d1, 3);
        *stackp++ = new;
        break;

      case OP_EL_DIV_BY:
        d1 = *--stackp;
        d2 = *--stackp;;
        new = assign_op (d2, d1, 4);
        *stackp++ = new;
        break;

      case OP_NOT:
        d = pop ();
        new = not (d);
        push (new);
        break;

      case OP_IFSJMP:
      {
        int d_val, ijmp;

        /*
         * Jump logic:
         * Pop stack. If pop'ed value is TRUE continue execution at same
         * point. If pop'ed value is FALSE jump to the specified value
         * one place behind the instruction.
         */

        d = pop ();   /* The truth value */
        ijmp = (int) (*pc).op_code; /* The jump value */
        d_val = 0;    /* Initialize */

        /*
         * Now figure out what the truth-datum means...
         * coerce it to either a zero, or non-zero.
         */

        switch (d.type)
        {
          case CONSTANT:
          case iCONSTANT:
            d_val = (int) d.u.val;
            break;

          case VAR:
          {
            Ent *ent = var_ent (d.u.ptr);
            switch (ent_type (ent))
            {
              case DOUBLE:
                d_val = (int) ent_double (ent);
                break;

              default:
                d_val = (int) class_logical_scalar (ent);
                break;
            }
            break;
          }

          case ENTITY:
          {
            Ent *ent = (Ent *) d.u.ptr;
            switch (ent_type (ent))
            {
              case DOUBLE:
                d_val = (int) ent_double (ent);
                ent_Destroy (ent);
                break;

              default:
                d_val = (int) class_logical_scalar (ent);
                ent_Destroy (ent);
                break;
            }
            break;
          }

          default:
            rerror ("terrible IFSJMP error");
            break;
        }

        /*
         * Finally, do the jump if called for.
         */

        if (d_val)
          /* Test value is non-zero, continue execution. */
          pc++;
        else
          /* Test value is zero, jump. */
          pc = pc + (ijmp);

        break;
      }

      case OP_FOR_THEN:
      {
        /* Get the current for-loop */
        ForLoop *fl = fp->flstack;
        fl->jmp_offset = (*pc++).op_code;

//         printf(THIS_FILE ": " THIS_SOLVER ": OP_FOR_THEN fl->jmp_offset = %i\n", fl->jmp_offset);
      }

      case OP_FOR_LOOP:
      {

        /* Get the current for-loop */
        ForLoop *fl = fp->flstack;

        /* Check for end-of-loop. */
        if (fl->cnt <= fl->n)
        {
          ent_Destroy (var_ent (fl->var));
          var_ent (fl->var) = class_forloop_value (fl->ent, fl->cnt++);
          ent_IncRef (var_ent (fl->var));

//           printf(THIS_FILE ": " THIS_SOLVER ": OP_FOR_LOOP fl->jmp        = %i\n", fl->jmp);
//           printf(THIS_FILE ": " THIS_SOLVER ": OP_FOR_LOOP fl->jmp_offset = %i\n", fl->jmp_offset);

          /* Jump back to top of loop. */
          if (fl->jmp_offset)
            pc = pc - fl->jmp + fl->jmp_offset;
          else
            pc = pc - fl->jmp;
        }

        /* else continue on to OP_FOR_LOOP_DONE */
        break;
      }
      case OP_SWAP:
        d2 = pop ();
        d1 = pop ();
        push (d2);
        push (d1);
        break;

      case OP_INC:
        inc ();
        break;

      case OP_DEC:
        dec ();
        break;

      case OP_POP:
        --stackp;
        if ((*stackp).type == ENTITY)
          ent_DecRef ((*stackp).u.ptr);
        break;

      case OP_POP_CLEAN:
        --stackp;
        if ((*stackp).type == ENTITY)
          ent_Destroy ((*stackp).u.ptr);
        break;

      case OP_VECTOR_CREATE:
        n = (*pc++).op_code;  /* get create-code from machine */
        switch (n)
        {
          case 2:
            d1 = pop ();
            d2 = pop ();
            new = vector_create_op (2, d1, d2, d3);
            break;
          case 3:
            d1 = pop ();
            d2 = pop ();
            d3 = pop ();
            new = vector_create_op (3, d1, d2, d3);
            break;
          default:
            rerror ("error in vector creation");
            break;
        }
        push (new);
        break;

      case OP_VEC_APPEND:
        d2 = pop ();
        d1 = pop ();
        new = vector_append (d1, d2);
        push (new);
        break;

      case OP_MATRIX_VEC_SUB:
        d2 = pop ();
        d1 = pop ();
        new = matrix_vector_sub (d1, d2);
        push (new);
        break;

      case OP_MATRIX_VEC_ASSIGN:
        d3 = pop ();
        d2 = pop ();
        d1 = pop ();
        new = matrix_vector_assign (d1, d2, d3);
        push (new);
        break;

      case OP_MATRIX_CREATE:
        d = pop ();
        new = matrix_create (d);
        push (new);
        break;

      case OP_MATRIX_APPEND:
        d2 = pop ();
        d1 = pop ();        
        new = matrix_stack (d1, d2);
        push (new);
        break;

      case OP_MATRIX_ASSIGN:
        n = (*pc++).op_code;  /* Get key to # of indices on stack. */
        switch (n)
        {
          case 1:     /* Both row and column indices. */
            d = pop ();
            d1 = pop ();
            d2 = pop ();
            d3 = pop ();
            new = matrix_assign_1 (d3, d2, d1, d);
            break;
            case 2:     /* Row index only. */
              d = pop ();
              d1 = pop ();
              d2 = pop ();
              new = matrix_assign_2 (d2, d1, d);
              break;
              case 3:     /* Column index only. */
                d = pop ();
                d1 = pop ();
                d2 = pop ();
                new = matrix_assign_3 (d2, d1, d);
                break;
        }
        push (new);
        break;

      case OP_MATRIX_SUB:
        n = (*pc++).op_code;  /* get key to # of indices on stack */
        switch (n)
        {
          case 1:     /* both row and column indices */
            d1 = pop ();
            d2 = pop ();
            d3 = pop ();
            new = matrix_sub_1 (d2, d1, d3);
            break;

          case 2:     /* row index only */
            d1 = pop ();
            d2 = pop ();
            new = matrix_sub_2 (d1, d2);
            break;

          case 3:     /* column index only */
            d1 = pop ();
            d2 = pop ();
            new = matrix_sub_3 (d1, d2);
            break;
        }
        push (new);
        break;

      case OP_LIST_CREATE:
        list_create ();
        break;

      case OP_LIST_MEMB:
        list_member ();
        break;

      case OP_LIST_ASSIGN:
        list_assign ();
        break;

      case OP_LIST_EL_CREATE:
        list_el_create ();
        break;

      case OP_FUNCTION_CALL:
        function_call ();
        break;

      case OP_FUNCTION_CALL_SELF:
        function_call_self ();
        break;

      case OP_FUNCTION_RETURN:
        function_return ();
        break;

      case OP_DEF_FUNC_RET:
        function_default_return ();
        break;

      case OP_TRANSPOSE:
        d = pop ();
        new = matrix_transpose (d);
        push (new);
        break;

      case OP_PUSH_STRING:
      {
        Ent *ent = ent_Create ();
        new.u.ptr = ent;
        new.type = ENTITY;
        ent_IncRef (ent);

        ent_data (ent) = mds_CreateScalar ( (*pc++).ptr );
        ent_SetType (ent, MATRIX_DENSE_STRING);
        *stackp++ = new;
        break;
      }
      case OP_QUIT:
        quit_code ();
        break;

      case OP_LINE_NO:
        /* No operation */
        pc++;     /* Inc program counter past line number */
        break;

      case OP_FILE_NAME:
        /* No operation */
        pc++;     /* Inc program counter past file name */
        break;

      case OP_JMP:
        pc = pc + (*pc).op_code;
        break;

      case OP_EMPTY_MATRIX_CREATE:
        new = empty_create ();
        *stackp++ = new;
        break;

      case OP_MATRIX_COL:
        d = pop ();
        new = matrix_reshape_col (d);
        push (new);
        break;

      case OP_EL_TRANSPOSE:
        d = pop ();
        new = matrix_el_transpose (d);
        push (new);
        break;

      case OP_RFILE:
        rfile ();
        break;

      case OP_RFILE_NAME:
        rfile_load (cpstr ((*pc++).ptr));
        break;

      case OP_REQ_NAME:
        require (cpstr ((*pc++).ptr));
        break;

      case OP_OLIST_ASSIGN:
                                /*
        * Pop the RHS (the list), N, and
        * the N LHS entities.
        * Do the assignment.
                                */

        d = *--stackp;
        d1 = *--stackp;
        /* Let olist_assign() do the rest of the popping. */
        new = olist_assign (d, d1);
        *stackp++ = new;
        break;

      case OP_OLIST:
                                /*
        * All the elements are already on the stack,
        * simply push N.
                                */

        n = (*pc++).op_code;
        new.type = CONSTANT;
        new.u.val = (double) n;
        *stackp++ = new;
        break;

      case OP_HELP:
        help ();
        break;

      case OP_HELP_NAME:
        help_name (cpstr ((*pc++).ptr));
        break;

      case OP_PUSH_UNDEF:
      {
        Ent *ent;
        new.u.ptr = listNode_Create ();
        new.type = VAR;
        ent = ent_Create ();
        ent_SetType (ent, UNDEF);
        listNode_AttachEnt (new.u.ptr, ent);
        push (new);
        break;
      }

      case OP_DEF_CLASS_START:
        inc_class_scope();
        break;

      case OP_DEF_CLASS_RET:
        class_return ();
        dec_class_scope();
        break;

      case OP_SAVE_EVAL:
        //
        // pop() the return value and store it
        // in a global var, Eval() will get it later.
        //
        eval_ret = pop ();
        break;

      case OP_FUNCTION_CALL_1:
        function_call_1 ();
        break;

      case OP_FUNCTION_CALL_2:
        function_call_2 ();
        break;

      case OP_SYS_CMD:
        new.u.val = (double) system (cpstr ((*pc++).ptr));
        new.type = CONSTANT;
        push (new);
        break;

      default:
        fprintf (stderr, "Invalid op-code: %d\n", (*--pc).op_code);
        fflush (stderr);
        rerror ("invalid op-code for interpreter");
        break;
    }
  }

//   printf(THIS_FILE ": " THIS_SOLVER ": EXITING 'execute'\n");
}

void
execute_debug (Inst * p)
{
  rerror ("can't debug yet");
}

/* **************************************************************
 * Increment a variable (postifx operation).
 *
 * - Pop the stack, getting the entity to increment....
 * - Copy the entity...
 * - Increment the original...
 * - Push the original back on the stack.
 * ************************************************************** */

void
inc (void)
{
  Datum d, new;
  ListNode *var = 0;		/* The original variable. */
  Ent *ent = 0;			/* The original entity. */
  Ent *cent = 0;		/* The entity copy. */

  d = pop ();			/* Pop the original off the stack. */
  if (d.type != VAR)
    rerror ("cannot increment (++) a CONSTANT");

  var = (ListNode *) (d.u.ptr);
  ent = var_ent (var);

  /* Copy the entity, so we can push it back on the stack. */
  cent = class_copy (ent);

  /* Inc the ref-count, cause we're going to put it on the stack. */
  ent_IncRef (cent);

  /* Now do the increment.
   * Check the ref-count first though.
   */

  if (ent_Ref (ent) > 1)
  {
    /* Make a copy, so we don't affect the others. */
    ent = ent_Create ();
    ent = class_copy (var_ent (var));
    ent_DecRef (var_ent (var));
    ent_IncRef (ent);
    /* Re-attach. */
    listNode_AttachEnt (var, ent);
  }

  /* Do the increment. */
  ent = class_increment (ent);

  /* Set up the return Datum. */
  new.u.ptr = cent;
  new.type = ENTITY;

  push (new);
}

/* **************************************************************
 * Decrement a variable
 * ************************************************************** */

void
dec (void)
{
  Datum d, new;
  ListNode *var = 0;
  Ent *ent = 0;
  Ent *cent = 0;

  d = pop ();
  if (d.type != VAR)
    rerror ("cannot increment (++) a CONSTANT");

  var = (ListNode *) (d.u.ptr);
  ent = var_ent (var);

  /* Copy the entity, so we can push it back on the stack. */
  cent = class_copy (ent);

  /* Inc the ref-count, cause we're going to put it on the stack. */
  ent_IncRef (cent);

  /* Now do the increment.
   * Check the ref-count first though.
   */

  if (ent_Ref (ent) > 1)
  {
    /* Make a copy, so we don't affect the others. */
    ent = ent_Create ();
    ent = class_copy (var_ent (var));
    ent_DecRef (var_ent (var));
    ent_IncRef (ent);
    /* Re-attach. */
    listNode_AttachEnt (var, ent);
  }

  /* Now do the increment. */
  ent = class_decrement (ent);

  /* Set up the return Datum. */
  new.u.ptr = cent;
  new.type = ENTITY;

  push (new);
}

/* **************************************************************
 * Assign the value of d2 to d1,  d1 = d2
 * d1 must be in the symbol table.
 * d2 must be an acceptable type.
 * ************************************************************** */

#undef  THIS_SOLVER
#define THIS_SOLVER "assign"
Datum assign (Datum d1, Datum d2)
{
  ListNode *lhs_var;
  Ent *lhs_ent, *nent, *rhs_ent;
  int killrhsvar;

  rhs_ent = 0;

  //
  // Check for bad or protected LHS values.
  //
  lhs_var = (ListNode *) (d1.u.ptr);
  lhs_ent = var_ent (lhs_var);

  if (ent_type (lhs_ent) == BLTIN)
  {
    if (ent_Ref(lhs_ent)<2)
      rerror ("cannot assign to a built-in function");
  }
  else if (ent_type (lhs_ent) == BTREE)
  {
    if (var_name (lhs_var) != 0)
    {
      Btree *bt = ent_data (lhs_ent);
      if (bt->isconst)
        rerror ("cannot assign to a protected variable!");
      //
      // protected variables: $$, mks, _rlab_config
      //
      //            if (!strncmp (var_name (lhs_var), "$$", 2))
      //              rerror ("cannot destroy global symbol table");
      //            if (!strncmp (var_name (lhs_var), "mks", 3))
      //              rerror ("cannot destroy protected variable 'mks'");
      //            if (!strncmp (var_name (lhs_var), "_rlab_config", 11))
      //              rerror ("cannot destroy protected variable '_rlab_config'");
    }
  }



  //
  // Now switch on the RHS value.
  //
  switch (d2.type)
  {
    /* RHS */
    case CONSTANT:
      ent_Destroy (lhs_ent);  /* Destroy old lhs. */
      nent = ent_Create (); /* Create new entity. */
      ent_data (nent) = mdr_CreateScalar (d2.u.val);
      listNode_AttachEnt (lhs_var, nent); /* Attach new entity. */
      ent_SetType (nent, MATRIX_DENSE_REAL);
      nent->refc += 1;
      break;

    case iCONSTANT:
      ent_Destroy (lhs_ent);
      nent = ent_Create ();
      ent_data (nent) = mdc_CreateScalar (0.0, d2.u.val);
      listNode_AttachEnt (lhs_var, nent);
      ent_SetType (nent, MATRIX_DENSE_COMPLEX);
      nent->refc += 1;
      break;

    case ENTITY:
      rhs_ent = (Ent *) (d2.u.ptr);
      if (rhs_ent)
        ent_DecRef (rhs_ent);
    case VAR:
      if (rhs_ent == 0)
        rhs_ent = var_ent (d2.u.ptr);

      switch (ent_type (rhs_ent))
      {
        case DOUBLE:
          ent_Destroy (lhs_ent);  /* Destroy old lhs. */
          listNode_AttachEnt (lhs_var, rhs_ent);  /* Attach rhs entity. */
          ent_SetType (rhs_ent, DOUBLE);
          rhs_ent->refc += 1;
          break;

        case iDOUBLE:
          ent_Destroy (lhs_ent);  /* Destroy old lhs. */
          listNode_AttachEnt (lhs_var, rhs_ent);  /* Attach rhs entity. */
          ent_SetType (rhs_ent, iDOUBLE);
          rhs_ent->refc += 1;
          break;

        default:
          /* Check the RHS-VAR for later destroy. */
          if (d2.type == VAR && (var_scope (d2.u.ptr) == 0))
            killrhsvar = 1;
          else
            killrhsvar = 0;

          /* Do the data. */
          listNode_AttachEnt (lhs_var, rhs_ent);
          rhs_ent->refc += 1;
          ent_Destroy (lhs_ent);

          /* Do the variable's list-entity. */
          if (var_listent (lhs_var))
          {
            ent_Destroy (var_listent (lhs_var));
            lhs_var->listent = 0;
          }

          /* Copy the list-structure of the RHS...
           * IFF the rhs-listent exists. */
          if ((d2.type == VAR) && (var_listent (d2.u.ptr)))
          {
            Ent *rent=var_listent(d2.u.ptr);
            var_listent(lhs_var) = rent;
            //Ent *lent = ent_Create ();
            //listNode_AttachListEnt (lhs_var, lent);
            //ent_data (lent) = btree_Copy (ent_data (rent));
            //ent_type (lent) = BTREE;
            //ent_IncRef (lent);
            ent_IncRef (rent);
          }

          if (killrhsvar)
            listNode_DestroyLoner ((ListNode *) (d2.u.ptr));
          break;
      }
      break;
  }

  return (d1);
}

/* **************************************************************
 * Assign the value of d2 to d1,  d1 = d2
 * d1 must be in the symbol table.
 * d2 must be an acceptable type.
 * op_idx:
 *  1:  d1 += d2
 *  2:  d1 -= d2
 *  3:  d1 *= d2
 *  4:  d1 /= d2
 * ************************************************************** */

Datum assign_op (Datum d1, Datum d2, int op_idx)
{
  ListNode *lhs_var;
  Ent *lhs_ent=0, *enew=0, *rhs_ent=0;

  //
  // Check for bad or protected LHS values.
  //
  lhs_var = (ListNode *) (d1.u.ptr);
  lhs_ent = var_ent (lhs_var);
  if (ent_type (lhs_ent) == BLTIN)
    rerror ("cannot assign to a built-in function");
  else if (ent_type (lhs_ent) == BTREE)
  {
    if (var_name (lhs_var) != 0)
    {
      Btree *bt = ent_data (lhs_ent);
      if (bt->isconst)
        rerror ("cannot assign to a protected variable!");
    }
  }

  rhs_ent = convert_datum_to_ent (d2);

  if (ent_Ref(lhs_ent)<=1)
  {
    switch(op_idx)
    {
      case 1:
        class_addto (lhs_ent, rhs_ent);
        break;

      case 2:
        class_subtractfrom (lhs_ent, rhs_ent);
        break;

      case 3:
        class_elmultiplyby(lhs_ent, rhs_ent);
        break;

      case 4:
        class_elrdivideby(lhs_ent, rhs_ent);
        break;
    }
    var_ent (lhs_var) = lhs_ent;
    ent_IncRef (lhs_ent);
  }
  else
  {
    switch(op_idx)
    {
      case 1:
        enew = class_add (lhs_ent, rhs_ent);
        break;

      case 2:
        enew = class_subtract (lhs_ent, rhs_ent);
        break;

      case 3:
        enew = class_el_multiply (lhs_ent, rhs_ent);
        break;

      case 4:
        enew = class_el_rdivide (lhs_ent, rhs_ent);
        break;

      default:
        rerror ("Horrible internal error: Help! I need somebody! Just Anybody! Help!\n");
    }
    ent_Clean (lhs_ent);
    var_ent (lhs_var) = enew;
    ent_IncRef (enew);
  }

  ent_Clean (rhs_ent);
  return (d1);
}

/* **************************************************************
 * Pop top value from stack, and print it. If the object is a
 * matrix then page it to stdout.
 * ************************************************************** */

extern void rpclose (FILE * fp);

static void assign_var (ListNode *lhs_var, Datum d)
{
  // get entity associated with the left variable
  Ent *lhs_ent = var_ent (lhs_var);
  ent_DecRef(lhs_ent); // dont kill it yet!

  // get entity associated with the right datum
  Ent *rhs_ent = bltin_get_ent (d);

  // Attach rhs to lhs as new entity
  if (d.type != ENTITY)
    ent_IncRef(rhs_ent);

  listNode_AttachEnt (lhs_var, rhs_ent);
  ent_Clean (lhs_ent);
  return;
}

#undef  THIS_SOLVER
#define THIS_SOLVER "print_and_assign"
void print_and_assign (Datum d)
{
  FILE *rlab_stderr = (!RLAB_STDERR_DS) ? stderr : RLAB_STDERR_DS;
  Ent *ent=0;

  ListNode *lf=btree_FindNode (get_symtab_ptr (), HAVE_ANSWER_FOR_ALL);
  if (!lf)
  {
    Ent *e=ent_Create();
    ent_type(e) = UNDEF;
    install(0, HAVE_ANSWER_FOR_ALL, e);
    lf= btree_FindNode (get_symtab_ptr (), HAVE_ANSWER_FOR_ALL);
  }
  assign_var (lf, d);
  fprintf(rlab_stderr, HAVE_ANSWER_FOR_ALL " =\n");

  ent = (Ent *) var_ent (lf);
  class_print (ent);
  ent_Clean (ent);

  fprintf(rlab_stderr, "\n");
  fflush (rlab_stderr);
}

void print (Datum d)
{
  FILE *rlab_stderr = (!RLAB_STDERR_DS) ? stderr : RLAB_STDERR_DS;
  Ent *ent = 0;
  int fwidth, fprec;

  fwidth = get_fwidth ();
  fprec = get_fprec ();

  switch (d.type)
  {
    case CONSTANT:
      fprintf (rlab_stderr, "%*.*g\n", fwidth, fprec, d.u.val);
      break;

    case iCONSTANT:
      fprintf (rlab_stderr, "%*.*gi\n", fwidth, fprec, d.u.val);
      break;

    case ENTITY:
      ent = (Ent *) (d.u.ptr);
      ent_DecRef (ent);

    case VAR:
      if (ent == 0)
        ent = var_ent (d.u.ptr);
      if(ent)
      {
        switch (ent_type (ent))
        {
          case DOUBLE:
          {
            double dv;
            dv = ent_double (ent);
            fprintf (rlab_stderr, "%*.*g\n", fwidth, fprec, dv);
            break;
          }

          default:
            class_print (ent);
            break;
        }
      }
      ent_Clean (ent);
      break;
  }

  fflush (rlab_stderr);
}


/* **************************************************************
 * Evaluate built-in on top of stack. Builtins get their arguments
 * in an array of entities. A single entity pointer is
 * returned.
 * ************************************************************** */

#define BLTIN_MAX_ARGS  32

void bltin (Ent * esp, int nargs, int popf)
{
  int i;
  Ent *(*built_in_func) ();
  Bltin *built_in;
  Datum args[BLTIN_MAX_ARGS], d;
  Ent *rent;

  built_in = (Bltin *) ent_data (esp);
  built_in_func = ((Bltin *) built_in)->func;

  /* Pop args off stack and install in entity array. */
  if (nargs > BLTIN_MAX_ARGS)
    rerror ("too many arguments to function");

  for (i = nargs - 1; i >= 0; i--)
  {
    args[i] = pop ();

    /* Decrement the reference count. */
    if (args[i].type == ENTITY)
      ent_DecRef (args[i].u.ptr);
  }

  /* Sometimes we have to pop the function itself off the stack. */
  if (popf)
    extern_pop ();

  /* call built-in function */
  rent = (Ent *) (*built_in_func) (nargs, args);

  /* Push the return entity back on the stack. */
  ent_IncRef (rent);
  d.type = ENTITY;
  d.u.ptr = rent;
  push (d);
}

/* **************************************************************
 * Create a list. The next instruction contains the number of
 * items on the stack that we need to install on the new list.
 *
 * For each entity we install, create a ListNode, with the proper
 * key, and inc the ref-count. Datums of type: LVAR (List-Var) are
 * to be moved directly into the list, variable and all.
 * ************************************************************** */

void
list_create (void)
{
  char jstr[20];
  int i, j, n;
  Datum d, nlist;
  Ent *ent;

  /*
   * Create the new list.
   */

  nlist.u.ptr = ent_Create ();
  nlist.type = ENTITY;
  ent_SetType (nlist.u.ptr, BTREE);
  ent_data (nlist.u.ptr) = btree_Create ();

  /* Get the # of items on the stack */
  n = (*pc++).op_code;

  j = n;
  /* Now pop each entity off the stack */
  for (i = 0; i < n; i++)
  {
    d = pop ();
    sprintf (jstr, "%d", j--);

    switch (d.type)
    {
      case CONSTANT:
        ent = ent_Create ();
        ent_data (ent) = mdr_CreateScalar (d.u.val);
        ent_SetType (ent, MATRIX_DENSE_REAL);
        install (ent_data (nlist.u.ptr), jstr, ent);
        break;

      case iCONSTANT:
        ent = ent_Create ();
        ent_data (ent) = mdc_CreateScalar (0.0, d.u.val);
        ent_SetType (ent, MATRIX_DENSE_COMPLEX);
        install (ent_data (nlist.u.ptr), jstr, ent);
        break;

      case ENTITY:
        ent_DecRef (d.u.ptr);
        install (ent_data (nlist.u.ptr), jstr, d.u.ptr);
        break;

      case VAR:
        install (ent_data (nlist.u.ptr), jstr, var_ent (d.u.ptr));
        break;

      case LVAR:
        listNode_SetOwned ((ListNode *) (d.u.ptr));
        btree_AddNode (ent_data (nlist.u.ptr), d.u.ptr);
        break;
    }
  }
  ent_IncRef (nlist.u.ptr);
  push (nlist);
}

/* **************************************************************
 * Resolve a list member reference.
 * code 1: The argument is an expr
 * code 2: The argument is a string (NAME)
 * Either way the stack must be popped to get the variable
 * id pointer, and maybe the vec_expr.
 *
 * This function is responsible for breaking apart lists that
 * are referenced by more than one variable, so that if/when
 * an assign happens, or the list-member is passed to a function,
 * the member's immediate parent has a ref-count of 1.
 *
 * This function handles more than just lists. It must also handle
 * generic class-member references. Each class, can have data members
 * (read-only ?) that can be referenced with the syntax C.m
 * Additionally, each variable can have a sublist that can be refereced
 * with the same syntax.
 * ************************************************************** */
#undef  THIS_SOLVER
#define THIS_SOLVER "list_member"
void list_member (void)
{

  char name[100], *ntmp=0;
  int code, rtype;
  void *vtmp;
  Datum new, vexpr, vid;
  Ent *ent, *tent;
  ListNode *var;

  code = (*pc++).op_code;

  if (code == 1)
    vexpr = pop ();

  /*
   * Pop the list-root off of the stack.
   */

  vid = pop ();

  if (code == 1)
  {
    /*
     * An expr was used as index.
     * Evaluate the expr, and convert result to a string
     * The string value goes in `name'.
     */

    switch (vexpr.type)
    {
      case CONSTANT:
        sprintf (name, "%.6g", vexpr.u.val);
        break;

      case iCONSTANT:
        sprintf (name, "%.6gi", vexpr.u.val);
        break;

      case ENTITY:
        ntmp = class_char_pointer (vexpr.u.ptr);
        if (isvalidstring(ntmp)<1)
          rerror ("null-string cannot be used for list-member reference!");
        sprintf (name, "%s", ntmp);
        ent_DecRef (vexpr.u.ptr); /* For coming off the stack. */
        ent_Clean (vexpr.u.ptr);
        break;

      case VAR:
        ntmp = class_char_pointer (var_ent (vexpr.u.ptr));
        if (isvalidstring(ntmp)<1)
          rerror ("null-string cannot be used for list-member reference!");
        sprintf (name, "%s", ntmp);
        break;

      default:
        rerror ("invalid list-member reference ");
        break;
    }
  }
  else if (code == 2)
  {
    /* list.NAME was used to index. Use name to index into list */
    strcpy (name, (char *) (*pc++).ptr);
  }

  /*
   * Now check type of root entity and extract the proper
   * object. The root entity can be a normal object like
   * a matrix, in this case we look for certain predefined
   * object members.
   */

  switch (vid.type)
  {
    case CONSTANT:
    case iCONSTANT:
      rerror ("cannot use CONSTANTs as lists");
      break;

    case ENTITY:
      /* We should only be using lists this way if they haven't yet
       * been attached to a variable. I.E: eig(a).vec. So there should
       * only be one reference count (after coming off the stack). If
       * this is not true, then there is a serious logic/design error.
       * Thus, the following test. If an entity shows up with a refc
       * greater than one, stop execution, and issue a nasty message.
       * mk: after introducing MDE's this warning is removed as it trips
       * every time a list that is attached to a MDE is accessed.
       */
      ent = (Ent *) (vid.u.ptr);
      ent_DecRef (ent);		/* For coming off the stack. */
//       if (ent->refc > 1)
//       {
//         fprintf (stderr, "!!!! ref-count = %i, %s\n", ent->refc, etd(ent));
//         rerror ("Stop in member-list reference, tried to reference an ENTITY!");
//       }
      if ((vtmp = class_member_ref (ent, name, &rtype)))
      {
        /* Object/Entity had the member. */
        new.type = rtype;
        new.u.ptr = vtmp;
        if (rtype == ENTITY)
          ent_IncRef (vtmp);
      }
      else
      {
        /* Object did not have the member, we can't do anything,
         * Since we don't have access to the parent (if there is one).
         */
        new.type = ENTITY;
        new.u.ptr = ent_Create ();
        ent_type (new.u.ptr) = UNDEF;
        ent_IncRef (new.u.ptr);
      }
      break;

  case VAR:

    /* This is the most common case. If we run across a list,
     * with a reference  count > 1, then we copy the list into a
     * new entity (we don't copy the data, just the list structure.
     * the data get their reference counts bumped up). This way
     * operations on sub-lists and list-elements are safe (I hope).
     */

    var = (ListNode *) (vid.u.ptr);
    ent = var_ent (var);

    if ((ent_type (ent) == BTREE) && (ent_Ref (ent) > 1))
    {
      /* Copy the tree, and make a new entity for it. */
      tent = ent_Create ();
      ent_data (tent) = btree_Copy (ent_data (ent));
      ent_type (tent) = BTREE;
      listNode_AttachEnt (var, tent);
      ent_IncRef (tent);
      ent_DecRef (ent);
    }
    else
    {
      /* Do nothing. */
      tent = ent;
    }

    if ((vtmp = class_member_ref (tent, name, &rtype)))
    {
      /* Check the list-members, or if its not a list,
       * the read-only members. */
      new.type = rtype;
      new.u.ptr = vtmp;
      if (rtype == ENTITY)
        ent_IncRef (vtmp);
    }
    else
    {
      /* Check the variable's sub-list. */
      if (var->listent)
      {
        vtmp = class_member_ref (var->listent, name, &rtype);
        new.type = rtype;
        new.u.ptr = vtmp;
        if (rtype == ENTITY)
          ent_IncRef (vtmp);
      }
      else
      {
        /* Nothing there, return an UNDEF entity. */
        new.type = ENTITY;
        new.u.ptr = ent_Create ();
        ent_type (new.u.ptr) = UNDEF;
        ent_IncRef (new.u.ptr);
      }
    }
    break;

  default:
    rerror ("invalid type for list member reference");
    break;
  }

  push (new);
}

/* **************************************************************
 * Assign to a list member. This can be kind of tricky.
 * For:        A.B = C
 *
 * First we must check A. If it belongs to someone else then we
 * must make a copy (of A only), Then, before we do any assigning, we
 * must check B. If B belongs to someone else (refc > 1), then
 * we must copy B also.
 * ************************************************************** */

void
list_assign (void)
{
  char name[100], *sname, *ntmp;
  int code;
  Datum d, expr, listd;
  Btree *btree = 0;
  ListNode *list, *listm;
  Ent *listent, *lent, *ent, *nent;
  Ent *rhs_ent = 0;

  list = 0;
  lent = 0;
  code = (*pc++).op_code;	/* Flag that tells us to use an
				   expression index, or a NAME
				   index */

  d = pop ();			/* The RHS */

  /*
   * Now take care of the list index.
   */

  if (code == 1)
  {
    /*
     * A vec_expr was used as index.
     * Evaluate the vexpr, and convert result to a string
     */

    expr = pop ();		/* The index */

    switch (expr.type)
    {
    case CONSTANT:
      sprintf (name, "%.6g", expr.u.val);
      break;

    case iCONSTANT:
      sprintf (name, "%.6gi", expr.u.val);
      break;

    case ENTITY:
      ntmp = class_char_pointer (expr.u.ptr);
      sprintf (name, "%s", ntmp);
      ent_Destroy (expr.u.ptr);
      break;

    case VAR:
      ntmp = class_char_pointer (var_ent (expr.u.ptr));
      sprintf (name, "%s", ntmp);
      break;

    default:
      rerror ("invalid list-member reference ");
      break;
    }
  }
  else if (code == 2)
  {
    /* list.NAME was used to index list */
    sname = (char *) (*pc++).ptr;
    strcpy (name, sname);
  }

  /*
   * Set up the list (LHS).
   */

  listd = pop ();		/* List variable or entity. */
  if (listd.type == VAR)
  {
    list = listd.u.ptr;
    lent = var_ent (list);
  }
  else
  {
    rerror ("cannot perform list-assignment");
  }

  /*
   * Check the list-parent (A), to see if anyone
   * else owns it. If refc > 1, then copy the object
   * before making the assignment.
   * If refc <= 1, then no copy is necessary.
   */

  if (ent_Ref (lent) > 1)
  {
    /* We must make a new entity. */
    ent_DecRef (lent);
    ent = class_copy (lent);
    listNode_AttachEnt (list, ent);
    ent_IncRef (ent);
  }
  else
  {
    /* Nobody else own the entity, we may go-ahead and modify. */
    ent = lent;
  }

  /*
   * Make sure the list-parent (LHS) is not UNDEFINED.
   * If it is, make it into a list...
   */

  if (ent_type (ent) == UNDEF)
  {
    ent_data (ent) = btree_Create ();
    ent_SetType (ent, BTREE);
  }

  /*
   * Now get on with the assignment...
   * Find the root of the tree we are going to assign to.
   * Check (in order):
   *    list: list-root
   *   other: read-only elements
   *          variable-sublist
   */

  if (ent_type (ent) == BTREE)
  {
    btree = ent_data (ent);
  }
  else
  {
    if (!class_attribute (ent, name))
    {
      if (list->listent)
      {
        btree = ent_data (list->listent);
      }
      else
      {
        btree = btree_Create ();
        listent = ent_Create ();
        ent_data (listent) = btree;
        ent_type (listent) = BTREE;
        listNode_AttachListEnt (list, listent);
        ent_IncRef (listent);
      }
    }
    else
    {
      /* Error, we cannot assign to the read-only
       * elements (attributes) of an object.
       */
      fprintf (stderr, "object element: %s is a read-only attribute\n", name);
      rerror ("list-member-assign: cannot assign to read-only elements");
    }
  }


  /* See if we can find the referenced member. */

  if ((listm = btree_FindNode (btree, name)) == 0)
  {
    /*
     * We can't find it, so create a new listnode.
     */

    listm = listNode_Create ();
    listNode_SetKey (listm, name);
    btree_AddNode (btree, listm);
  }
  else
  {
    /* Check the ref-count of listm. */

    if (ent_Ref (var_ent (listm)) > 1)
      ent_DecRef (var_ent (listm));
  }

  /* Now do the assign. */

  switch (d.type)
  {
  case CONSTANT:
    nent = ent_Create ();
    ent_data (nent) = mdr_CreateScalar (d.u.val);
    listNode_AttachEnt (listm, nent);
    ent_SetType (nent, MATRIX_DENSE_REAL);
    ent_IncRef (nent);
    break;

  case iCONSTANT:
    nent = ent_Create ();
    ent_data (nent) = mdc_CreateScalar (0.0, d.u.val);
    listNode_AttachEnt (listm, nent);
    ent_SetType (nent, MATRIX_DENSE_COMPLEX);
    ent_IncRef (nent);
    break;

  case ENTITY:
    rhs_ent = (Ent *) (d.u.ptr);
    ent_DecRef (rhs_ent);
    /* Fall through */

  case VAR:
    if (rhs_ent == 0)
      rhs_ent = var_ent (d.u.ptr);

    switch (ent_type (rhs_ent))
    {
    case DOUBLE:
      listNode_AttachEnt (listm, rhs_ent);
      ent_IncRef (rhs_ent);
      break;

    case iDOUBLE:
      listNode_AttachEnt (listm, rhs_ent);
      ent_IncRef (rhs_ent);
      break;

    default:
      listNode_AttachEnt (listm, rhs_ent);
      ent_IncRef (rhs_ent);
      break;
    }
    break;
  }

  /* Push the LHS back on the stack */
  push (listd);
}

/* **************************************************************
 * Given: NAME '=' vec_expr
 *
 * Create a new variable (with key = NAME) that points to vec_expr,
 * Push this var back on the stack, as type: LVAR. Type LVAR signifies
 * that the Datum is only to be used for list construction.
 * ************************************************************** */

void
list_el_create (void)
{
  char *name;
  Datum d, rhs;
  Ent *ent;
  ListNode *lnode;

  lnode = 0;

  /* Get the list-var-name from the machine */
  name = cpstr ((*pc++).ptr);

  /* Pop the rhs from the stack */
  rhs = pop ();

  /* Create the LVAR. */
  lnode = listNode_Create ();
  listNode_SetKey (lnode, name);

  /* Now attach an entity to the List-VAR. */
  switch (rhs.type)
  {
  case CONSTANT:
    ent = ent_Create ();
    ent_data (ent) = mdr_CreateScalar (rhs.u.val);
    ent_SetType (ent, MATRIX_DENSE_REAL);
    listNode_AttachEnt (lnode, ent);
    break;

  case iCONSTANT:
    ent = ent_Create ();
    ent_data (ent) = mdc_CreateScalar (0.0, rhs.u.val);
    ent_SetType (ent, MATRIX_DENSE_COMPLEX);
    listNode_AttachEnt (lnode, ent);
    break;

  case ENTITY:
    ent_DecRef (rhs.u.ptr);	/* For coming off the stack. */
    listNode_AttachEnt (lnode, rhs.u.ptr);
    break;

  case VAR:
    listNode_AttachEnt (lnode, var_ent (rhs.u.ptr));
    /* Check for a variable's list-entity. */
    if (var_listent ((rhs.u.ptr)))
    {
      Ent *lent, *rent;
      rent = var_listent (rhs.u.ptr);
      lent = ent_Create ();
      listNode_AttachListEnt (lnode, lent);
      ent_data (lent) = btree_Copy (ent_data (rent));
      ent_type (lent) = BTREE;
      ent_IncRef (lent);
    }
    break;
  }

  d.type = LVAR;
  d.u.ptr = lnode;
  ent_IncRef (var_ent (lnode));	/* For being attached to a new var. */
  push (d);
}

/*
 * < s1; s2; s3; ... > = LIST;
 *
 * drhs:   The RHS list structure.
 * dN:     The number of LHS entities to assign.
 */

Datum
olist_assign (Datum drhs, Datum dN)
{
  char **rhs_names;
  int N, i, nrhs;
  Btree *rhs;
  Datum d1, d2, new;
  ListNode *ltmp;

  N = 0;
  rhs = 0;			/* Initialize */

  /* Check N. */
  if (dN.type == CONSTANT)
  {
    N = (int) dN.u.val;		/* The number of assignments to do. */
  }
  else
  {
    rerror ("olist_assign: terrible error");
  }

  /*
   * Now loop N times, popping the stack, and
   * doing the assignment. If there are more RHS
   * elements than LHS, that is OK. If there are
   * more LHS, than RHS, error.
   */

  /* Check rhs type, and extract rhs. */

  if ((drhs.type != ENTITY) && (drhs.type != VAR))
  {
    rerror ("invalid RHS for open-list assign");
  }
  else if (drhs.type == VAR)
  {
    if (ent_type (var_ent (drhs.u.ptr)) != BTREE)
      rerror ("invalid RHS for open-list assign, must be LIST");
    rhs = (Btree *) ent_data (var_ent (drhs.u.ptr));
  }
  else if (drhs.type == ENTITY)
  {
    if (ent_type (drhs.u.ptr) != BTREE)
      rerror ("invalid RHS for open-list assign, must be LIST");
    rhs = (Btree *) ent_data (drhs.u.ptr);
    ent_DecRef (drhs.u.ptr);
  }
  else
  {
    rerror ("invalid RHS for open-list assign");
  }

  /* Check N-rhs and N-lhs. */
  if ((nrhs = btree_GetRealNumNodes (rhs)) >= N)
  {
    /* OK */ ;
  }
  else
  {
    rerror ("RHS length must be >= LHS length");
  }

  /*
   * Get the names in the RHS tree, so that we can extract the
   * elements alphabetically.
   */

  rhs_names = (char **) GC_MALLOC (nrhs * sizeof (char *));
  get_btree_node_names (rhs->root_node, rhs_names);

  /*
   * Now loop over the list(s) and do the assigns.
   * The rhs-list will not necessarily be destroyed.
   * So, we must be careful to handle the reference
   * counts properly...
   */

  for (i = N - 1; i >= 0; i--)
  {
    d1 = pop ();
    ltmp = btree_FindNode (rhs, rhs_names[i]);

    d2.u.ptr = ltmp;
    d2.type = VAR;

    /*
     * We must use assign(), and do a copy, cause
     * we don't know if drhs is a function return or
     * a regular list.
     */

    new = assign (d1, d2);
  }

  /* Clean Up. */

  for (i = 0; i < nrhs; i++)
  {
    GC_FREE (rhs_names[i]);
  }
  GC_FREE (rhs_names);

  if (drhs.type == ENTITY)
    ent_Clean (drhs.u.ptr);
  return (new);
}

/* **************************************************************
 * Find the line number ascociated with the current error. Go
 * FORWARD from the current program instruction untill we find a
 * OP_LINE_NO.
 * ************************************************************** */

int
find_lineno_new (void)
{
  Inst *p;

  p = pc;

  if (p)
  {
    while ((*p).op_code <= OP_SYS_CMD)
    {
      if ( (*p).op_code == OP_LINE_NO )
        return ((*(p + 1)).op_code - loff);
      p++;
    }
    p = pc;
    while ((*p).op_code <= OP_SYS_CMD)
    {
      if ( (*p).op_code == OP_LINE_NO )
        return ((*(p + 1)).op_code - loff);
      p--;
    }
  }

  return (1);
}

int find_lineno (void)
{
  Inst *p;

  p = pc;
  if (!p)
    return (1);     /* This is very rare */

  while ((*p).op_code != OP_LINE_NO)
  {
    p++;
  }

  return ((*(p + 1)).op_code - loff);
}

char *
find_file_name (void)
{
  char *fptr;

  if ((fp - frame) == 0)
  {
    fptr = (char *) program->prog[1].ptr;
  }
  else
  {
    Ent *esp = fp->esp;
    Function *func = ent_data (esp);
    fptr = (char *) func->code[1].ptr;
  }
  return (fptr);
}

/* **************************************************************
 * We leave RLaB via this indirect route in order to make `quit'
 * an executable statement. Also we try and clean up memory
 * in order to make debugging memory leaks easier.
 * ************************************************************** */

#ifdef HAVE_READLINE
extern void stifle_history (int);
extern int write_history (const char *);
#endif

void
quit_code (void)
{

#ifdef HAVE_READLINE
  if (rl_histfile)
  {
    stifle_history (rl_histsize);
    write_history (rl_histfile);
  }
#endif /* HAVE_READLINE */

  exit (1);
}
