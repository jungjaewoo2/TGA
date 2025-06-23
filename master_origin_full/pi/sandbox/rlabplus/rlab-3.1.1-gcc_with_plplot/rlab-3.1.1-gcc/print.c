/* rfileio.c */

/*  This file is a part of RLaB ("Our"-LaB)
    Copyright (C) 1995  Ian R. Searle
    Copyright (C) 2007  Marijan Kostrun

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
 ***********************************************************************/

#include "rlab.h"
#include "symbol.h"
#include "bltin.h"
#include "listnode.h"
#include "mem.h"
#include "list.h"
#include "util.h"
#include "rfileio.h"
#include "mds.h"
#include "class.h"
#include "scan.h"

#include "rlab_solver_parameters_names.h"

#include <stdio.h>

/* getline.c */
extern char scan_code[];

static char *endp;

/*---------- types and defs for doing printf ------------*/
#define  PF_C     0
#define  PF_S     1
#define  PF_D     2 /* int conversion */
#define  PF_LD    3 /* int */
#define  PF_F     4 /* float conversion */

typedef int (*PRINTER) (VPTR, char *, ...);

/* for switch on number of '*' and type */
#define  AST(num,type)  (5*(num)+(type))

static int do_printf (FILE * fp, char *format, int n_args,
		      int arg_cnt, Datum args[], char **strv);

/* **************************************************************
 * RlaB builting functions for emulating ANSI-C printf, fprintf,
 * sprintf, open, close.
 * ************************************************************** */


/*
 * Emulate C printf().
 */
#undef  THIS_SOLVER
#define THIS_SOLVER "printf"
Ent *
Printf (int nargs, Datum args[])
{
  FILE *fp = stdout;;
  Ent *e1=0;
  int retval=-1;
  char *format=0;

  if (nargs == 0)
  {
    fprintf(stderr, THIS_SOLVER ": " RLAB_ERROR_ARG1_MDS_SCALAR "\n");
    goto _exit_printf;
  }

  // 1st arg MUST be format
  e1 = bltin_get_ent (args[0]);
  format = class_char_pointer (e1);
  if (!format)
  {
    fprintf(stderr, THIS_SOLVER ": " RLAB_ERROR_ARG1_MDS_SCALAR "\n");
    goto _exit_printf;
  }

  retval = do_printf (fp, format, nargs, 0, args, (char **) 0);

_exit_printf:

  ent_Clean (e1);
  return ent_Create_Rlab_Double((double) retval);
}

/*
 * Emulate C fprintf().
 */
#undef  THIS_SOLVER
#define THIS_SOLVER "fprintf"
Ent *
FPrintf (int nargs, Datum args[])
{
  FILE *fp=0;
  Ent *e1=0, *e2=0;
  char *filenm=0, *format=0;
  int retval=-1;

  if (nargs < 2)
  {
    fprintf(stderr, THIS_SOLVER ": " RLAB_ERROR_AT_LEAST_TWO_ARG_REQUIRED "\n");
    goto _exit_fprintf;
  }

  // filename
  e1 = bltin_get_ent (args[0]);
  filenm = class_char_pointer (e1);
  if (!filenm)
  {
    fprintf(stderr, THIS_SOLVER ": " RLAB_ERROR_ARG1_MDS_SCALAR "\n");
    goto _exit_fprintf;
  }

  // format string
  e2 = bltin_get_ent (args[1]);
  format = class_char_pointer (e2);
  if (!format)
  {
    fprintf(stderr, THIS_SOLVER ": " RLAB_ERROR_ARG2_MDS_SCALAR "\n");
    goto _exit_fprintf;
  }

  if ((fp = get_file_ds (filenm, "w", 0)))
    retval = do_printf (fp, format, nargs, 1, args, (char **) 0);

_exit_fprintf:

  ent_Clean (e1);
  ent_Clean (e2);
  return ent_Create_Rlab_Double((double) retval);
}

/*
 * Emulate C sprintf()
 */
#undef  THIS_SOLVER
#define THIS_SOLVER "sprintf"
Ent *
SPrintf (int nargs, Datum args[])
{
  char *format=0, **strv=0;
  int retval=-1;
  Ent *e1=0, *e2=0;
  ListNode *var=0;
  MDS *ms=0;

  if (nargs < 2)
  {
    fprintf(stderr, THIS_SOLVER ": " RLAB_ERROR_ARG2_MDS_SCALAR "\n");
    goto _exit_sprintf;
  }

  // 1st arg MUST be variable defined or undefined, but with name!
  if (args[0].type != VAR)
  {
    fprintf(stderr, THIS_SOLVER ": " RLAB_ERROR_ARG1_VAR "\n");
    goto _exit_sprintf;
  }
  if (!var_key(args[0].u.ptr))
  {
    fprintf(stderr, THIS_SOLVER ": Cannot leave first argument empty in call to function '" THIS_SOLVER "'\n");
    goto _exit_sprintf;
  }

  var = (ListNode *) (args[0].u.ptr);
  e1 = var_ent (var);

  // `e1' must be a variable...
  if (ent_type (e1) != MATRIX_DENSE_STRING)
  {
    // Destroy the old e1, and make a string-scalar
    if (e1->refc > 1)
      ent_DecRef (e1);
    else
      ent_Destroy (e1);

    e1 = ent_Create_Rlab_String("");
    listNode_AttachEnt (var, e1);
    ent_IncRef (e1);
  }
  else
  {
    if (e1->refc > 1)
    {
      // We must create a new entity to use.
      // First, dec the ref on the supplied entity.
      ent_DecRef (e1);
      e1 = ent_Create_Rlab_String("");
      listNode_AttachEnt (var, e1);
      ent_IncRef (e1);
    }

    // Make sure e1 is a 1-by-1 string-matrix.
    ms = ent_data (e1);
    if (!EQSCAL(ms))
    {
      mds_Destroy (ms);
      e1 = ent_Create_Rlab_String("");
    }
  }

  // Pass ptr-to-ptr to do_printf().
  ms = ent_data (e1);
  strv = MDSPTR(ms);

  e2 = bltin_get_ent (args[1]);
  format = class_char_pointer (e2);

  retval = do_printf ((FILE *) 0, format, nargs, 1, args, strv);

_exit_sprintf:

  ent_Clean (e2);
  return ent_Create_Rlab_Double((double) retval);
}

#define MAX_SPRINTF_SIZE  32768
char sprintf_buff[MAX_SPRINTF_SIZE];
char *sprintf_limit = sprintf_buff + MAX_SPRINTF_SIZE;

static int
do_printf (FILE * fp, char *format, int nargs, int arg_cnt, Datum args[], char **strv)
{
  char save;
  char *p;
  register char *q;
  char *target;
  int l_flag, h_flag;		/* seen %ld or %hd  */
  int ast_cnt;
  int ast[2];
  double dval = 0.0;
  int ival = 0;
  char sval[MAX_SPRINTF_SIZE];
  char *stmp;
  int num_conversion = 0;	/* for error messages */
  int pf_type = 0;		/* conversion type */
  PRINTER printer;		/* pts at fprintf() or sprintf() */
  int argcnt;
  int retval;			/* Attempt to return C-printf() return-val */
  char *nothing_to_print="";
  Ent *earg;

  argcnt = nargs - 1;		/* To count backward */
  q = format;
  retval = 0;

  if (fp == (FILE *) 0)   /* doing sprintf */
  {
    if (*strv == 0)
      rerror ("terrible error in [fs]printf()");
    target = sprintf_buff;
    printer = (PRINTER) sprintf;
  }
  else

  {
    /* doing printf */
    target = (char *) fp; /* will never change */
    printer = (PRINTER) fprintf;
  }

  /* Traverse format string, doing printf(). */
  while (1)
  {
    if (fp)     /* printf */
    {
      while (*q != '%')
      {
        if (*q == 0)
        {
          fflush (fp);
          return (retval);
        }
        else
        {
          putc (*q, fp);
          q++;
          retval++;
        }
      }
    }
    else
    {
      /* sprintf() */
      while (*q != '%')
      {
        if (*q == 0)
        {
          if (target > sprintf_limit) /* damaged */
            rerror ("sprintf problem, buffer too small");
          else
          {
            /* really done */
            int len = target - sprintf_buff;
            GC_FREE (*strv);
            *strv = (char *) GC_MALLOC (sizeof (char) * (len + 1));
            memcpy (*strv, sprintf_buff, (size_t) len);
            (*strv)[len] = '\0';
            return (retval);
          }
        }
        else
        {
          *target++ = *q++;
          retval++;
        }
      }
    }

    num_conversion++;

    if (*++q == '%')    /* %% */
    {
      if (fp)
      {
        putc (*q, fp);
        retval++;
      }
      else
      {
        *target++ = *q;
      }
      q++;
      continue;
    }

    /*
     * We have found a conversion specifier, figure it out,
     * then print the data asociated with it.
     */

    if (argcnt <= 0)
    {
      (*printer) ((VPTR) target, "\n");
      fflush ((FILE *) target);
      warning_1 ("printf: not enough arguments");
      return (-1);
    }
    else
    {
      /* Get The data object from the arg-list */
      earg = bltin_get_ent (args[++arg_cnt]);
      argcnt--;
    }

    /* mark the '%' with p */
    p = q - 1;

    /* eat the flags */
    while (*q == '-' || *q == '+' || *q == ' ' || *q == '#' || *q == '0')
      q++;

    ast_cnt = 0;    /* asterisk count */
    if (*q == '*')
    {
      /* Use current arg as field width spec */
      ast[ast_cnt++] = (int) class_double (earg);
      q++;

      if (argcnt <= 0)
      {
        (*printer) ((VPTR) target, "\n");
        fflush ((FILE *) target);
        warning_1 ("printf: not enough arguments");
        return (-1);
      }
      else
      {
        /* Get next arg */
        earg = bltin_get_ent (args[++arg_cnt]);
        argcnt--;
      }
    }
    else
      while (scan_code[*(unsigned char *) q] == SC_DIGIT)
        q++;
    /* width is done */

    if (*q == '.')    /* have precision */
    {
      q++;
      if (*q == '*')
      {
        /* Use current arg as precision spec */
        ast[ast_cnt++] = (int) class_double (earg);
        q++;

        if (argcnt <= 0)
        {
          (*printer) ((VPTR) target, "\n");
          fflush ((FILE *) target);
          warning_1 ("printf: not enough arguments");
          return (-1);
        }
        else
        {
          /* Get next arg */
          earg = bltin_get_ent (args[++arg_cnt]);
          argcnt--;
        }
      }
      else
        while (scan_code[*(unsigned char *) q] == SC_DIGIT)
          q++;
    }

    if (argcnt < 0)
    {
      (*printer) ((VPTR) target, "\n");
      fflush ((FILE *) target);
      warning_1 ("printf: not enough arguments");
      return (-1);
    }

    l_flag = h_flag = 0;

    if (*q == 'l')
    {
      q++;
      l_flag = 1;
    }
    else if (*q == 'h')
    {
      q++;
      h_flag = 1;
    }

    /* Set pf_type and load val */
    switch (*q++)
    {
      case 's':
        if (l_flag + h_flag)
          rerror ("printf: bad conversion");

        stmp = class_char_pointer (earg);
        if (!stmp)
          stmp = nothing_to_print;
        if (strlen (stmp) > MAX_SPRINTF_SIZE)
          rerror ("printf: string too long for printf buffer!");
        strcpy (sval, stmp);

        pf_type = PF_S;
        break;

      case 'c':
        ival = (int) strtod (class_char_pointer (earg), (char **) &endp);
        pf_type = PF_C;
        break;

      case 'o':
        rerror ("printf: \"o\" format not allowed");
        break;

      case 'd':
      case 'x':
      case 'X':
        ival = class_int (earg);
        pf_type = PF_D;
        break;

      case 'i':
      case 'u':
        /* use strod() here */
        ival = class_int (earg);
        pf_type = l_flag ? PF_LD : PF_D;
        break;

      case 'e':
      case 'g':
      case 'f':
      case 'E':
      case 'G':
        if (h_flag + l_flag)
          rerror ("printf: bad conversion");
        /* use strod() here */
        dval = class_double (earg);
        pf_type = PF_F;
        break;

      default:
        rerror ("printf: bad conversion");
    }

    save = *q;
    *q = 0;

    /* ready to call printf() */
    /*
     * target:   The output file (or variable for sprintf())
     * p:        the beginning of the format
     * ast:      array with asterisk values
     */
    switch (AST (ast_cnt, pf_type))
    {
      case AST (0, PF_C):
        retval += (*printer) ((VPTR) target, p, ival);
        break;

      case AST (1, PF_C):
        retval += (*printer) ((VPTR) target, p, ast[0], ival);
        break;

      case AST (2, PF_C):
        retval += (*printer) ((VPTR) target, p, ast[0], ast[1], ival);
        break;

      case AST (0, PF_S):
        retval += (*printer) ((VPTR) target, p, sval);
        break;

      case AST (1, PF_S):
        retval += (*printer) ((VPTR) target, p, ast[0], sval);
        break;

      case AST (2, PF_S):
        retval += (*printer) ((VPTR) target, p, ast[0], ast[1], sval);
        break;

      case AST (0, PF_D):
        retval += (*printer) ((VPTR) target, p, ival);
        break;

      case AST (1, PF_D):
        retval += (*printer) ((VPTR) target, p, ast[0], ival);
        break;

      case AST (2, PF_D):
        retval += (*printer) ((VPTR) target, p, ast[0], ast[1], ival);
        break;

      case AST (0, PF_LD):
        retval += (*printer) ((VPTR) target, p, ival);
        break;

      case AST (1, PF_LD):
        retval += (*printer) ((VPTR) target, p, ast[0], ival);
        break;

      case AST (2, PF_LD):
        retval += (*printer) ((VPTR) target, p, ast[0], ast[1], ival);
        break;

      case AST (0, PF_F):
        retval += (*printer) ((VPTR) target, p, dval);
        break;

      case AST (1, PF_F):
        retval += (*printer) ((VPTR) target, p, ast[0], dval);
        break;

      case AST (2, PF_F):
        retval += (*printer) ((VPTR) target, p, ast[0], ast[1], dval);
        break;
    }
    if (fp == (FILE *) 0)
      while (*target)
        target++;
    *q = save;
  }
}

/* **************************************************************
 * Some thoughts on a scanf() for rlab. Nothing implemented yet.
 *
 * Emulate the C-scanf family of functions.
 * Args:
 *    fp:      The file pointer. If this is 0, then do scanf from
 *             a string (see strv).
 *    format:  The scanf format string.
 *    n_args:  The number of arguments to the rlab scanf.
 *    arg_cnt: The current position in the argument array.
 *    d_arg:   The argument array.
 *    strv:    The string to do sscanf() from.
 *
 * Logic: Read the format string, until we find a vaid
 *        format/conversion string (something with a `%'). Then
 *        Use that format string, to use the libc scanf to do the
 *        work. We must look at the format string in order to get
 *        the scanf data pointer of the correct type. Not to
 *        difficult, it is either numeric, or string for us.
 *
 * Returns: do_scanf returns a list of the elements read in. If
 *          EOF or and error is encountered, then the list will
 *          be empty.
 *
 * At this point I am not sure why I should go through the trouble.
 * This scanf will not offer anything getline doesn't.
 * ************************************************************** */
