/* bltin.c */

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
#include "symbol.h"
#include "mem.h"
#include "bltin.h"
#include "util.h"
#include "mathl.h"
#include "class.h"
#include "rlab_solver_parameters_names.h"

#include "mds.h"
#include "function.h"

#define  THIS_FILE "bltin.c"

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>

/* **************************************************************
 * User's error() function. longjmp() back to prompt.
 * ************************************************************** */

Ent *
Error (int nargs, Datum args[])
{
  char *emsg;
  Ent *e = 0;
  Ent *rent = ent_Create ();

  if (nargs == 0)
  {
    /* Default error message */
    rerror ("USER-RAISED-ERROR");
  }
  else
  {
    e = bltin_get_ent (args[0]);
    emsg = class_char_pointer (e);
    rerror (emsg);
  }

  /*
   * Shut up compiler, execution will never get this far.
   */

  ent_Clean (e);
  return (rent);
}

Ent *
Stop (int nargs, Datum args[])
{
  char *s=0;
  Ent *e1=0;

  if (nargs==1)
  {
    e1 = bltin_get_ent ( args[0] );
    s = class_char_pointer(e1);
  }
  stop_script( s );

  ent_Clean (e1);

  return ent_Create_Rlab_Success();
}

Ent *
Exist (int nargs, Datum args[])
{
  Ent *e1=0;
  double rval=0;

  if (nargs == 1)
  {
    e1 = bltin_get_ent (args[0]);
    if (ent_type (e1) != UNDEF)
      rval = 1;
  }

  ent_Clean (e1);
  return ent_Create_Rlab_Double(rval);
}

extern Btree *static_tree;
#undef   THIS_SOLVER
#define  THIS_SOLVER "entinfo"
Ent * EntInfo (int nargs, Datum args[])
{
  char stmp[32];
  Btree *bt=0, *bt2=0;
  Ent *e=0, *rent;
  int i;

  if (nargs != 1)
    rerror (THIS_FILE ": " THIS_SOLVER ": " RLAB_ERROR_ONE_ARG_REQUIRED);

  dprintf("datum type = %i (VAR=%i)\n", args[0].type, VAR);

  e = bltin_get_ent (args[0]);
  if (ent_type (e) == UNDEF)
  {
    ent_Clean (e);

    rent = ent_Create ();
    ent_data (rent) = mdr_CreateScalar (0.0);
    ent_SetType (rent, MATRIX_DENSE_REAL);
    return (rent);
  }

  /* Set up the return-list. */
  bt = btree_Create ();

  /* Set up the string holding the entity address. */
  sprintf (stmp, "%p", e);
  Ent *addr1 = ent_Create ();
  ent_data (addr1) = mds_CreateScalar(stmp);
  ent_SetType (addr1, MATRIX_DENSE_STRING);
  install (bt, RLAB_ENTINFO_ADDR, addr1);

  /* Set up the string holding the entity address. */
  sprintf (stmp, "%p", ent_data(e));
  Ent *addr2 = ent_Create ();
  ent_data (addr2) = mds_CreateScalar(stmp);
  ent_SetType (addr2, MATRIX_DENSE_STRING);
  install (bt, RLAB_ENTINFO_DATAPTR, addr2);

  // if argument is user function, do some more digging
  if (ent_type (e) == U_FUNCTION)
  {
    Function *f = ent_data(e);
    List *stat = f->stat;
    if (stat)
    {
      bt2 = btree_Create();
      ListNode *lnode = list_GetLastNode(stat);
      for (i=0; i<f->n_stat; i++)
      {
        Ent  *ent = ent_Copy(var_ent(lnode));
        char *key = var_key(lnode);
        install (bt2, key, ent);
        lnode = listNode_GetNodeBehind(lnode);
      }
      Ent *s = ent_Create();
      ent_data (s) = bt2;
      ent_SetType (s, BTREE);
      install (bt, RLAB_ENTINFO_STAT, s);
    }
  }

  // if variable check more
  if (args[0].type==VAR)
  {
    Var *var = args[0].u.ptr;
    if (var_listent(var))
    {
      Ent *etmp = var_listent(var);
      dprintf("var_listent(var) = %p\n",etmp);
      dprintf("ent_type = %i (BTREE=%i)\n", ent_type(etmp), BTREE);
      if (ent_data(etmp))
      {
        dprintf("ent_data = %p\n", ent_data(etmp));
      }
    }
    
  }

  /* Set up the double, holding the entity ref-count. */
  Ent *refc = ent_Create ();
  ent_data (refc) = mdr_CreateScalar ((double) (e->refc));
  ent_SetType (refc, MATRIX_DENSE_REAL);
  install (bt, RLAB_ENTINFO_REFC, refc);


  ent_Clean (e);

  rent = ent_Create ();
  ent_data (rent) = bt;
  ent_SetType (rent, BTREE);
  return (rent);
}

/*
 * This is an undocumented function.
 * I was using this to fix/debug the symbol
 * table once apon a time, and decided to leave
 * it in.
 */
extern void fix_symbol_table (int debug);

Ent *
FixSymTable (int nargs, Datum args[])
{
  Ent *rent = 0;

  if (nargs != 0)
  {
    /* Default error message */
    rerror ("fsymtab: no arguments allowed");
  }

  fix_symbol_table (1);

  rent = ent_Create ();

  ent_data (rent) = mdr_CreateScalar (1.0);
  ent_SetType (rent, MATRIX_DENSE_REAL);

  return (rent);
}

void
bltin_Print (Bltin * dont_use, FILE * fptr)
{
  fprintf (fptr, "\t<bltin-function>\n");
}

/* **************************************************************
 * Return the class of a bltin function.
 * ************************************************************** */

char *
bltin_Class (Bltin * b)
{
  return (cpstr ("function"));
}

MDS *
bltin_Type_BF (Bltin * b)
{
  MDS *type = mds_CreateScalar ( RLAB_MEMBER_TYPE_F_BUILTIN );
  return (type);
}

size_t
bltin_Sizeof (Bltin * b)
{
  return ((size_t) 0);
}

Ent *
bltin_MemberRef (Bltin * bf, char *name, int *type)
{
  void *rptr = 0;
  Ent *ne;
  ne = ent_Create ();

  if (!strcmp (name, RLAB_MEMBER_CLASS))
  {
    ent_data (ne) = mds_CreateScalar ( RLAB_MEMBER_CLASS_F );
    ent_SetType (ne, MATRIX_DENSE_STRING);
  }
  else if (!strcmp (name, RLAB_MEMBER_TYPE))
  {
    ent_data (ne) = mds_CreateScalar ( RLAB_MEMBER_TYPE_F_BUILTIN );
    ent_SetType (ne, MATRIX_DENSE_STRING);
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
bltin_Members (Bltin * bf, int *n)
{
  char **marray = (char **) GC_MALLOC (2 * sizeof (char *));
  if (marray == 0)
    rerror ("out of memory");
  marray[0] = cpstr (RLAB_MEMBER_CLASS);
  marray[1] = cpstr (RLAB_MEMBER_TYPE);

  *n = 2;
  return (marray);
}

/* **************************************************************
 * Support for Diary function.
 * ************************************************************** */

#include "rfileio.h"

#include <sys/types.h>

#ifdef __STDC__
#include <time.h>
#endif /* __STDC__ */

#ifdef HAVE_TIME_H
#include <time.h>
#endif

static int write_diary = 0;
static FILE *diary_file_ptr;
static char *diary_filenm;

Ent *
Diary (int nargs, Datum args[])
{
  char *string;
  Ent *e1;
  Ent *rent = 0;
  FILE *fn;
  time_t r_time;

  if (nargs > 1)
  {
    rerror ("diary: 1 argument allowed");
  }

  /* Handle diary() args */
  if (nargs == 1)		/* 1 arg */
  {
    if (write_diary)
      rerror ("diary: only one diary file may be open at a time");

    e1 = bltin_get_ent (args[0]);
    string = class_char_pointer (e1);
  }
  else
  {				/* No args */
    if (write_diary)
    {
      /*
       * nargs = 0, and a diary file is already open
       * close the existing diary file, and return
       */

      close_file_ds (diary_filenm);
      write_diary = 0;
      diary_file_ptr = 0;

      rent = ent_Create ();
      ent_data (rent) = mdr_CreateScalar (1.0);
      ent_type (rent) = MATRIX_DENSE_REAL;
      return (rent);
    }
    else
    {
      /*
       * n_args = 0, and no diary file open, open the
       * default diary file "DIARY"
       */
      string = cpstr ("DIARY");
    }
  }

  if ((fn = get_file_ds (string, "w", 0)) == 0)
  {
    warning_1 ("cannot open for write");
    rent = ent_Create ();
    ent_data (rent) = mdr_CreateScalar (0.0);
    ent_type (rent) = MATRIX_DENSE_REAL;
    return (rent);
  }

  /* Set static variables */
  write_diary = 1;
  diary_file_ptr = fn;
  diary_filenm = cpstr (string);

  /* Write out a header to diary file */
  r_time = time (0);
  fprintf (fn, "// RLaB diary file: %s. Opened %s\n", string, ctime (&r_time));

  rent = ent_Create ();
  ent_data (rent) = mdr_CreateScalar (1.0);
  ent_type (rent) = MATRIX_DENSE_REAL;
  return (rent);
}

int
get_write_diary ()
{
  return (write_diary);
}

FILE *
get_diary_file_ptr ()
{
  return (diary_file_ptr);
}
