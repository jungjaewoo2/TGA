/* util.c */

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
#include "mem.h"
#include "list.h"
#include "bltin.h"
#include "btree.h"
#include "mdr.h"
#include "mdc.h"
#include "function.h"
#include "sort.h"
#include "rlab_solver_parameters_names.h"
#include "rlab_macros.h"

//#define MIN(x,y)  ((x) < (y) ? (x) :  (y))

#undef  THIS_FILE
#define THIS_FILE "util.c"

#include <string.h>
#include <setjmp.h>

#ifdef __riscos

#define ThrowbackWarning      0
#define ThrowbackError        1

void throwback (int, int, char *);
void throwback_finish (void);
void throwback_init (int, char *);
void set_throwback_enable (int);

static int use_throwback;

/* **************************************************************
 * Set/unset throwback enabling
 * ************************************************************** */
void set_throwback_enable (int value)
{
  use_throwback = value;
}
#endif


/* **************************************************************
 * Set the program name...
 * ************************************************************** */

static char *progname;
char *cpstr (char *string);

void set_progname (char *value)
{
#ifndef WIN32
  progname = cpstr (value);
#else
  if (!strcmp ("rlab", value))
    progname = "rlab.exe";
  else
    progname = cpstr (value);
#endif
}

char *get_progname ()
{
  return progname;
}

/* **************************************************************
 * Run-time error handling and signal catching...
 * ************************************************************** */
/*
 * These hold the env for longjmp()s back to the prompt,
 * and error exit
 */

#define NJMP 40			/* Depth of load() recursion */
static jmp_buf jmp[NJMP];
static int ijmp;		/* Counter for tracking jmp[] */

static char no_file[] = "no file info available";

/* **************************************************************
 * Jump buffer handling routines. These routines merely inc, or
 * an integer, whilst checking to see that the bounds of the array
 * have not been exceeded .
 * ************************************************************** */

int get_ijmp (void)
{
  return (ijmp);
}

int inc_buff (void)
{
  if (ijmp >= NJMP)
    fprintf (stderr, "NJMP too small\n");

  return (ijmp++);
}

int dec_buff (void)
{
  if (ijmp <= 0)
    fprintf (stderr, "error while decrementing ijmp\n");

  return (--ijmp);
}

jmp_buf *jmp_inc_buff (void)
{
  return (&jmp[inc_buff ()]);
}

jmp_buf *jmp_dec_buff (void)
{
  return (&jmp[dec_buff ()]);
}

/* **************************************************************
 * Recover from a run-time error. 1 name, 1 message
 * ************************************************************** */

void warning_1 (char *s);
static int line_nos;		/* if TRUE include line numbers in code */

void set_util_line_nos (int val)
{
  line_nos = val;
}

/* symbol.c */
extern void fix_symbol_table (int debug);
extern int symtab_debug;

#undef  THIS_SOLVER
#define THIS_SOLVER "rerror"
void rerror (char *s)
{
  warning_1 (s);

  /*
   * Clean up the Datum stack.
   */
  datum_stack_clean ();


  /*
   * Make necessary repairs to the symbol table.
   */
  fix_symbol_table (0);
  longjmp (jmp[dec_buff ()], 1);
}

void nerror (char *s)
{
  warning_1 (s);
  datum_stack_clean ();
  longjmp (jmp[dec_buff ()], 1);
}

void stop_script ( char *s )
{
  if(s!=0)
    fprintf(stderr,"%s\n",s);

  datum_stack_clean ();

  longjmp (jmp[dec_buff ()], 1);
}

/* **************************************************************
 * Print warning messages.
 * ************************************************************** */

void warning_1 (char *s)
{
  char *fn;
  int lineno;

  /* Print file and line info if available */
  if (line_nos)
  {
    fn = find_file_name ();
    lineno = find_lineno ();
  }
  else
  {
    fn = no_file;
    lineno = 0;
  }

  fprintf (stderr, "ERROR: %s: %s\n", progname, s);

  if (strcmp ("stdin", fn))
    fprintf (stderr, "\tnear line %d, file: %s\n", lineno, fn);

#ifdef __riscos
  if (use_throwback)
  {
    throwback_init (1, fn);
    throwback (ThrowbackError, lineno, fn);
  }
#endif

  fflush (stderr);
}

/* **************************************************************
 * Signal catching functions.
 * ************************************************************** */

#include "fpe.h"

#undef  THIS_FUNCTION
#define THIS_FUNCTION "fpecatch"
void fpecatch (int tmp)
{
#ifdef __EMX__
  /* to re-enable signals under EMX */
  signal (SIGFPE, SIG_ACK);
#endif

  setup_fpe_handling ();

  /*
   * This may be over-cautious (paranoid), but it is
   * possible that we could get another one of these
   * on our way back to the prompt, so lets ignore
   * any more that come along.
   */

  signal (SIGFPE, SIG_IGN);
  rerror (THIS_FUNCTION ": " RLAB_ERROR_FPE);
}

void pipecatch (int tmp)
{
#ifdef __EMX__
  /* to re-enable signals under EMX */
  signal (SIGPIPE, SIG_ACK);
#endif

#if !defined(__riscos) && !defined(WIN32)
  signal (SIGPIPE, SIG_IGN);
#endif
  /* Don't print anything, just jump */
  longjmp (jmp[dec_buff ()], 1);
}

/* GSL: sort.c
   Inline swap function for moving objects around.
   For number types use SWAP macro in sort.c instead! */
void swap (void *base, size_t size, size_t i, size_t j)
{
  register unsigned char *a = size * i + (unsigned char *) base;
  register unsigned char *b = size * j + (unsigned char *) base;
  register size_t s = size;
  register unsigned char tmp;

  if (i == j)
    return;

  do
  {
    tmp = *a;
    *a++ = *b;
    *b++ = tmp;
  }
  while (--s > 0);
}


/* **************************************************************
 * 'str1'-> 'str1str2'
 * ************************************************************** */
#undef  THIS_FUNCTION
#define THIS_FUNCTION "string_concat"
void string_concat (char **str1, char *str2)
{
  int len2 = isvalidstring(str2);
  if (len2<1)
    return; // nothing to do, or string 'str2' not defined

  int len1=isvalidstring(*str1);
  if (len1>0)
  {
    *str1 = (char *) GC_REALLOC (*str1, ((len1+len2) + 1) * sizeof (char));
  }
  else
  {
    *str1 = (char *) GC_MALLOC ((len2 + 1) * sizeof (char));
  }
  if (!*str1)
    rerror (THIS_FUNCTION ": " RLAB_ERROR_OUT_OF_MEMORY);

  strcpy(*str1 + len1, str2);

  return;
}

/* **************************************************************
 * Copy a string to a new char ptr, Return the ptr to the newly
 * created string.
 * ************************************************************** */
#undef  THIS_FUNCTION
#define THIS_FUNCTION "cpstr"
char *cpstr (char *string)
{
  char *new_string=0;

  if (string)
  {
    new_string = (char *) GC_MALLOC ((strlen (string) + 1) * sizeof (char));

    if (!new_string)
      rerror (THIS_FUNCTION ": " RLAB_ERROR_OUT_OF_MEMORY);

    strcpy (new_string, string);
  }

  return (new_string);
}

/* **************************************************************
 * Copy  a string to a new char ptr, Return the ptr to the newly
 * created string. Copy at most N chars. Always terminate with a
 * null character.
 * ************************************************************** */
#undef  THIS_FUNCTION
#define THIS_FUNCTION "cpnstr"

char *cpnstr (char *string, int n)
{
  char *new_string;
  if (string != 0)
  {
    int len = MIN(strlen(string),n);

    new_string = (char *) GC_MALLOC ((len + 1) * sizeof (char));

    if (!new_string)
      rerror (THIS_FUNCTION ": " RLAB_ERROR_OUT_OF_MEMORY);

    strncpy(new_string, string, len);
    new_string[len] = '\0';

    return (new_string);
  }
  return (0);
}

/* **************************************************************
 * Copy a string to a new char *, strip the surrounding "".
 * This function is primarily used by the RLaB scanner.
 * ************************************************************** */
#undef  THIS_FUNCTION
#define THIS_FUNCTION "cpstr_strip"

char *cpstr_strip (char *string)
{
  int len = isvalidstring(string);
  char *new_string;
  if (len>0)
  {
    string[len - 1] = '\0';	/* get rid of trailing " */
    string++;			/* get rid of 1st " */

    new_string = (char *) GC_MALLOC ((size_t) (len - 1));
    if (new_string == 0)
      rerror (THIS_FUNCTION ": " RLAB_ERROR_OUT_OF_MEMORY);

    strcpy (new_string, string);
    return (new_string);
  }

  return (0);
}

/* **************************************************************
 * Append a string to another, separated with whitespace.
 * Return the new result.
 * ************************************************************** */

#undef  THIS_FUNCTION
#define THIS_FUNCTION "strappend"
char *strappend (char *s1, char *s2)
// {
//   int l1=isvalidstring (s1);
//   int l2=isvalidstring (s2);
// 
//   char *s=0;
//   int l = (l1 > 0 ? l1 : 0) + (l2 > 0 ? l2 : 0);
// 
// #ifdef HAVE_GC
//   s = GC_MALLOC ((l + 2) * sizeof (char));
// #else
//   s = CALLOC ((l + 2), sizeof (char));
// #endif
// 
//   if (s == 0)
//     rerror (THIS_FUNCTION ": " RLAB_ERROR_OUT_OF_MEMORY);
// 
//   s[0] = '\0';
//   if (l1 > 0)
//   {
//     strcpy (s, s1);
//   }
//   if (l2 > 0)
//   {
//     strcat (s, " ");    /* Separate with whitespace. */
//     strcat (s, s2);
//   }
// 
//   return (s);
// }
{
  char *snew=0;
  int l1=isvalidstring (s1);
  int l2=isvalidstring (s2);

  /* Exercise a little caution. */
  if (l1 < 1 || l2 < 1)
  {
    if (l2 > 0)
    {
      snew = cpstr (s2);
    }
    else if (l1 > 0)
    {
      snew = cpstr (s1);
    }
    else
    {
#ifdef HAVE_GC
      snew = (char *) GC_MALLOC (sizeof (char));
#else
      snew = (char *) CALLOC (1, sizeof (char));
#endif
      if (snew == 0)
        rerror (THIS_FUNCTION ": " RLAB_ERROR_OUT_OF_MEMORY);
      snew[0] = '\0';
    }
    return (snew);
  }

#ifdef HAVE_GC
  snew = (char *) GC_MALLOC ((l1 + l2 + 2) * sizeof (char));
#else
  snew = (char *) CALLOC ((l1 + l2 + 2), sizeof (char));
#endif

  if (snew == 0)
    rerror ("out of memory");

  strcpy (snew, s1);
  strcat (snew, " ");    /* Separate with whitespace. */
  strcat (snew, s2);
  return (snew);
}

/* **************************************************************
 * Add two character strings together, return ptr to the newly
 * created string.
 * ************************************************************** */

#undef  THIS_FUNCTION
#define THIS_FUNCTION "string_add"
char *string_add (char *s1, char *s2)
{
  char *snew=0;
  int l1=isvalidstring (s1);
  int l2=isvalidstring (s2);

  /* Exercise a little caution. */
  if (l1 < 1 || l2 < 1)
  {
    if (l2 > 0)
    {
      snew = cpstr (s2);
    }
    else if (l1 > 0)
    {
      snew = cpstr (s1);
    }
    else
    {
#ifdef HAVE_GC
      snew = (char *) GC_MALLOC (sizeof (char));
#else
      snew = (char *) CALLOC (1, sizeof (char));
#endif
      if (snew == 0)
        rerror (THIS_FUNCTION ": " RLAB_ERROR_OUT_OF_MEMORY);
      snew[0] = '\0';
    }
    return (snew);
  }

#ifdef HAVE_GC
  snew = (char *) GC_MALLOC ((l1 + l2 + 1) * sizeof (char));
#else
  snew = (char *) CALLOC ((l1 + l2 + 1), sizeof (char));
#endif

  if (snew == 0)
    rerror ("out of memory");

  strcpy (snew, s1);
  strcpy (snew + l1, s2);
  return (snew);
}

/* **************************************************************
 * Convert a Datum to an Entity. This function is intended to be
 * used by the interpreter's operations. This function decrements
 * The entities reference count by one, since whenever an entity
 * comes off of the stack its count is decremented, and vice-versa.
 * ************************************************************** */
Ent *convert_datum_to_ent (Datum d)
{
  Ent *ent=0;

//   fprintf(stderr, "convert_datum_to_ent: d.type = %i\n", d.type);

  switch (d.type)
  {
    case GLOBAL:
      // this is here so that command line option containing
      // rlab2 -e "rfile somerfile.r"
      // doesn't cause segfault
      break;

    case CONSTANT:
      /* Refc is already zero. */
      return ent_Create_Rlab_Double(d.u.val);
      break;

    case iCONSTANT:
      /* Refc is already zero. */
      return ent_Create_Rlab_Complex(0.0, d.u.val);
      break;

    case ENTITY:
      ent = (Ent *) (d.u.ptr);
      if (ent)
        ent_DecRef (ent);
      else
      {
/*        ent = ent_Create();
        ent_type(ent) = UNDEF;
        return (ent);*/
        rerror ("Internal Error: cannot use UNDEF as ENTITY");
        break;
      }

    case VAR:
      {
        if (ent == 0)
        {
          if (d.u.ptr)
            ent = var_ent (d.u.ptr);
          else
          {
            ent = ent_Create();
            ent_type(ent) = UNDEF;
            return (ent);
//             rerror ("Internal Error: cannot use UNDEF as VAR");
            break;
          }
        }

//         fprintf(stderr, "ent_type(rent) = %i\n", ent_type (ent));

        switch (ent_type (ent))
        {
          case DOUBLE:
            ent_data (ent) = mdr_CreateScalar (ent_double (ent));
            if (d.type == VAR)
              listNode_DestroyLoner ((ListNode *) (d.u.ptr));
            return (ent);
            break;

          case iDOUBLE:
            ent_data (ent) = mdc_CreateScalar (0.0, ent_double (ent));
            if (d.type == VAR)
              listNode_DestroyLoner ((ListNode *) (d.u.ptr));
            return (ent);
            break;

          case UNDEF:
            if (d.type == VAR)
              listNode_DestroyLoner ((ListNode *) (d.u.ptr));
            ent = ent_Create();
//             rerror ("Internal Error: cannot use UNDEF as VAR");
            break;

          default:
            if (d.type == VAR)
              listNode_DestroyLoner ((ListNode *) (d.u.ptr));
            return (ent);
            break;
        }
        break;
      }

    default:
      rerror ("Terrible Error: cannot convert CONST to ENTITY");
      break;
  }

  return (ent);     /* Shut up the compiler. */
}

Ent *bltin_get_ent (Datum d)
{
  Ent *ent = 0;

  switch (d.type)
  {
    case CONSTANT:
      ent = ent_Create ();
      ent_data (ent) = mdr_CreateScalar (d.u.val);
      ent_SetType (ent, MATRIX_DENSE_REAL);

      /* Refc is already zero. */
      return (ent);
      break;

    case iCONSTANT:
      ent = ent_Create ();
      ent_data (ent) = mdc_CreateScalar (0.0, d.u.val);
      ent_SetType (ent, MATRIX_DENSE_COMPLEX);

      /* Refc is already zero. */
      return (ent);
      break;

    case ENTITY:
      ent = (Ent *) (d.u.ptr);
      /* Fall through */

    case VAR:
    {
      if (ent == 0)
        ent = var_ent (d.u.ptr);

      switch (ent_type (ent))
      {
        case DOUBLE:
          ent_data (ent) = mdr_CreateScalar (ent_double (ent));
          ent_double (ent) = 0.0;
          ent_SetType (ent, MATRIX_DENSE_REAL);
          if (d.type == VAR)
            listNode_DestroyLoner ((ListNode *) (d.u.ptr));

          return (ent);
          break;

        default:
          if (d.type == VAR)
            listNode_DestroyLoner ((ListNode *) (d.u.ptr));
          return (ent);
          break;
      }
      break;
    }

    default:
      rerror ("Terrible Error: cannot convert CONST to ENTITY");
      break;
  }

  return (ent);     /* Shut up the compiler. */
}

/* **************************************************************
 * Call a RLaB User-Function from C source code. The arguments are:
 * fname: A null-terminated character string containing the
 *        User-Function name.
 * args: A pointer to an array of the function arguments.
 * nargs: The number of arguments.
 *
 * call_rlab_script returns a void pointer that points to the
 * entity returned by the User-Function.
 * ************************************************************** */

extern void extern_push (Datum d);
extern Datum extern_pop (void);
extern void userf (Ent * esp, int nargs, int self);
extern void bltin (Ent * esp, int nargs, int popf);

Ent * get_ent_from_rlab_script (char *fname, Btree * gst_ptr)
{
  Ent *efunc=0;

  if (isvalidstring(fname)<1)
  {
    fprintf (stderr, "Invalid function name\n");
    rerror ("calling user-function");
  }

  ListNode *lfunc=0;

  if (!gst_ptr)
    lfunc = btree_FindNode (get_symtab_ptr(), fname);
  else
    lfunc = btree_FindNode (gst_ptr, fname);

  /* Check to see if fname exists */
  if (!lfunc)
  {
    fprintf (stderr, "Function %s does not exist in provided symbol-table\n", fname);
    rerror ("calling user-function");
  }

  efunc = var_ent (lfunc);
  if (!efunc)
  {
    fprintf (stderr, "%s is not a function\n", fname);
    rerror ("must be class \"function\"");
  }

  // Error-check the function
  if (ent_type (efunc) != U_FUNCTION && ent_type (efunc) != BLTIN)
  {
    fprintf (stderr, "%s is not a function\n", fname);
    rerror ("must be class \"function\"");
  }

  return (efunc);
}


Datum call_rlab_script_ent (Ent *efunc, Datum * args, int nargs)
{
  int i;

  Datum ret_datum;

  if (!efunc)
    rerror ("Horrible internal error: Null pointer for function encountered!");

  for (i = 0; i < nargs; i++)
    extern_push (args[i]);

  /*
   * Execute the User/Bltin Function code.
   * Do not increment the program counter (pc)
   * cause we are NOT calling this function
   * from an environment that is executing
   * other code.
   */

  if (ent_type (efunc) == BLTIN)
  {
    bltin (efunc, nargs, 0);
  }
  else
  {
    userf (efunc, nargs, 0);
  }

  // Pop the return entity
  ret_datum = extern_pop ();

  return (ret_datum);
}

Ent * ent_call_rlab_script_1arg (Ent *efunc, Ent * e1)
{
  Datum fun_args[1];

  Datum ret_datum;

  if (!efunc)
    rerror ("Horrible internal error: Null pointer for function encountered!");

  if (!e1)
    rerror ("Horrible internal error: Null pointer for arg1 encountered!");

  // assign input entities to datums:
  // 1:
  fun_args[0].u.ptr = e1;
  fun_args[0].type  = ENTITY;
  ent_IncRef (e1);

  // push them on stack
  extern_push (fun_args[0]);

  if (ent_type (efunc) == BLTIN)
    bltin (efunc, 1, 0);
  else
    userf (efunc, 1, 0);

  // Pop the return entity from stack
  ret_datum = extern_pop ();

  return bltin_get_ent(ret_datum);
}

Ent * ent_call_rlab_script_2args (Ent *efunc, Ent * e1, Ent * e2)
{
  int i;
  Datum fun_args[2];

  Datum ret_datum;

  if (!efunc)
    rerror ("Horrible internal error: Null pointer for function encountered!");

  if (!e1)
    rerror ("Horrible internal error: Null pointer for arg1 encountered!");

  if (!e2)
    rerror ("Horrible internal error: Null pointer for arg2 encountered!");

  // assign input entities to datums:
  // 1:
  fun_args[0].u.ptr = e1;
  fun_args[0].type  = ENTITY;
  ent_IncRef (e1);
  // 2:
  fun_args[1].u.ptr = e2;
  fun_args[1].type  = ENTITY;
  ent_IncRef (e2);

  // push them on stack
  for (i = 0; i < 2; i++)
    extern_push (fun_args[i]);

  if (ent_type (efunc) == BLTIN)
    bltin (efunc, 2, 0);
  else
    userf (efunc, 2, 0);

  // Pop the return entity from stack
  ret_datum = extern_pop ();

  return bltin_get_ent(ret_datum);
}

Ent * ent_call_rlab_script_3args (Ent *efunc, Ent * e1, Ent * e2, Ent * e3)
{
  int i;
  Datum fun_args[3];

  Datum ret_datum;

  if (!efunc)
    rerror ("Horrible internal error: Null pointer for function encountered!");

  if (!e1)
    rerror ("Horrible internal error: Null pointer for arg1 encountered!");

  if (!e2)
    rerror ("Horrible internal error: Null pointer for arg2 encountered!");

  if (!e3)
    rerror ("Horrible internal error: Null pointer for arg3 encountered!");

  // assign input entities to datums:
  // 1:
  fun_args[0].u.ptr = e1;
  fun_args[0].type  = ENTITY;
  ent_IncRef (e1);
  // 2:
  fun_args[1].u.ptr = e2;
  fun_args[1].type  = ENTITY;
  ent_IncRef (e2);
  // 3:
  fun_args[2].u.ptr = e3;
  fun_args[2].type  = ENTITY;
  ent_IncRef (e3);

  // push them on stack
  for (i = 0; i < 3; i++)
    extern_push (fun_args[i]);

  if (ent_type (efunc) == BLTIN)
    bltin (efunc, 3, 0);
  else
    userf (efunc, 3, 0);

  // Pop the return entity from stack
  ret_datum = extern_pop ();

  return bltin_get_ent(ret_datum);
}

Ent * ent_call_rlab_script_4args (Ent *efunc, Ent * e1, Ent * e2, Ent * e3, Ent * e4)
{
  int i;
  Datum fun_args[4];

  Datum ret_datum;

  if (!efunc)
    rerror ("Horrible internal error: Null pointer for function encountered!");

  if (!e1)
    rerror ("Horrible internal error: Null pointer for arg1 encountered!");

  if (!e2)
    rerror ("Horrible internal error: Null pointer for arg2 encountered!");

  if (!e3)
    rerror ("Horrible internal error: Null pointer for arg3 encountered!");

  if (!e4)
    rerror ("Horrible internal error: Null pointer for arg4 encountered!");

  // assign input entities to datums:
  // 1:
  fun_args[0].u.ptr = e1;
  fun_args[0].type  = ENTITY;
  ent_IncRef (e1);
  // 2:
  fun_args[1].u.ptr = e2;
  fun_args[1].type  = ENTITY;
  ent_IncRef (e2);
  // 3:
  fun_args[2].u.ptr = e3;
  fun_args[2].type  = ENTITY;
  ent_IncRef (e3);

  // 4:
  fun_args[3].u.ptr = e4;
  fun_args[3].type  = ENTITY;
  ent_IncRef (e4);

  // push them on stack
  for (i = 0; i < 4; i++)
    extern_push (fun_args[i]);

  if (ent_type (efunc) == BLTIN)
    bltin (efunc, 4, 0);
  else
    userf (efunc, 4, 0);

  // Pop the return entity from stack
  ret_datum = extern_pop ();

  return bltin_get_ent(ret_datum);
}

Ent * ent_call_rlab_script_5args (Ent *efunc, Ent * e1, Ent * e2, Ent * e3, Ent * e4, Ent *e5)
{
  int i;
  Datum fun_args[5];

  Datum ret_datum;

  if (!efunc)
    rerror ("Horrible internal error: Null pointer for function encountered!");

  if (!e1)
    rerror ("Horrible internal error: Null pointer for arg1 encountered!");

  if (!e2)
    rerror ("Horrible internal error: Null pointer for arg2 encountered!");

  if (!e3)
    rerror ("Horrible internal error: Null pointer for arg3 encountered!");

  if (!e4)
    rerror ("Horrible internal error: Null pointer for arg4 encountered!");

  if (!e5)
    rerror ("Horrible internal error: Null pointer for arg5 encountered!");

  // assign input entities to datums:
  // 1:
  fun_args[0].u.ptr = e1;
  fun_args[0].type  = ENTITY;
  ent_IncRef (e1);
  // 2:
  fun_args[1].u.ptr = e2;
  fun_args[1].type  = ENTITY;
  ent_IncRef (e2);
  // 3:
  fun_args[2].u.ptr = e3;
  fun_args[2].type  = ENTITY;
  ent_IncRef (e3);

  // 4:
  fun_args[3].u.ptr = e4;
  fun_args[3].type  = ENTITY;
  ent_IncRef (e4);

  // 5:
  fun_args[4].u.ptr = e5;
  fun_args[4].type  = ENTITY;
  ent_IncRef (e5);

  // push them on stack
  for (i = 0; i < 5; i++)
    extern_push (fun_args[i]);

  if (ent_type (efunc) == BLTIN)
    bltin (efunc, 5, 0);
  else
    userf (efunc, 5, 0);

  // Pop the return entity from stack
  ret_datum = extern_pop ();

  return bltin_get_ent(ret_datum);
}



Datum call_rlab_script (char *fname, Datum * args, int nargs)
{
  Ent *efunc = get_ent_from_rlab_script(fname, NULL);
  return call_rlab_script_ent(efunc, args, nargs);
}

int isfuncdatum(Datum x)
{
  int ival = 0;

  if (x.type == VAR)
  {
    Ent *e = bltin_get_ent(x);

    if (ent_type(e)==U_FUNCTION || ent_type(e)==BLTIN)
      ival = 1;

    ent_Clean(e);
  }

  return ival;
}




/* **************************************************************
 * Remove duplicate items from an integer array.
 * A pointer to the new array is returned.
 *
 * ia:     pointer to the input array of integers
 * size:   the size of the input array
 * nsize:  the size of the output array
 * ************************************************************** */

int *remove_duplicates (int *ia, int size, int *nsize)
{
  int dup, i, j, k, *nia, *tmp;

  /* Create the new array. */
  nia = (int *) GC_MAIOP (size * sizeof (int));

  k = 0;
  for (i = 0; i < size - 1; i++)
  {
    dup = 0;
    for (j = i + 1; j < size; j++)
    {
      if (ia[i] == ia[j])
      {
        dup = 1;
        break;
      }
    }
    if (!dup)
      nia[k++] = ia[i];
  }

  /* We always get the last one. */
  nia[k] = ia[size - 1];

  /* Now, fix the size of nia. */

  tmp = (int *) GC_MAIOP ((k + 1) * sizeof (int));
  memcpy (tmp, nia, (k + 1) * sizeof (int));

  GC_FREE (nia);
  *nsize = k + 1;
  return (tmp);
}

/* **************************************************************
 * Remove elements that are greater than A.
 * ************************************************************** */

int *remove_invalid (int *ia, int size, int A, int *nsize)
{
  int i, k, *nia, *tmp;

  /* Create the new array. */
    nia = (int *) GC_MAIOP (size * sizeof (int));

    k = 0;
  for (i = 0; i < size; i++)
  {
    if (ia[i] <= A)
    {
      /* Copy it. */
      nia[k++] = ia[i];
    }
  }

  /* Now, fix the size of nia. */

  tmp = (int *) GC_MAIOP (k * sizeof (int));
  memcpy (tmp, nia, k * sizeof (int));
  GC_FREE (nia);

  *nsize = k;
  return (tmp);
}

/* **************************************************************
 * Return the "set". That is, a vector with unique, ordered
 * elements.
 * ************************************************************** */

int *array_set (int *input, int input_size, int *output_size)
{
  int *tmp;

  /* First, remove the duplicate elements. */
    tmp = remove_duplicates (input, input_size, output_size);

  /* Order (sort) the result. */
    i_qsort (tmp, 0, *output_size - 1);

    return (tmp);
}

/* **************************************************************
 * Return the union of two vectors.
 * ************************************************************** */

int *array_union (int *in1, int in_size1, int *in2, int in_size2,
		  int *output_size)
{
  int *tmp, *utmp;

  /* Combine the two arrays into one. */
    tmp = (int *) GC_MAIOP ((in_size1 + in_size2) * sizeof (int));
    memcpy (tmp, in1, in_size1 * sizeof (int));
    memcpy (tmp + in_size1, in2, in_size2 * sizeof (int));

  /* Now, compute the set (the result). */
    utmp = array_set (tmp, in_size1 + in_size2, output_size);
    GC_FREE (tmp);
    return (utmp);
}
