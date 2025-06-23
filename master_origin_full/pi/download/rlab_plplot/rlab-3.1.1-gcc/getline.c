/* getline.c */

/* Source code for the RLaB getline function.
 *
 * Syntax: getline ( "filename" )
 *
 * Description:
 *
 * SEE THE GETLINE HELPFILE FOR UP-TO-DATE DESCRITPION
 */

/*
 * Much of the getline scanner is code that came from mawk,
 * by Mike Brennan (GNU Copylefted). I have un-mercifily chopped,
 * sliced, and diced it to meet RLaB's needs, so any bugs are
 * mine. The original logic is Mike's.
 */

/*  This file is a part of RLaB ("Our"-LaB)
    Copyright (C) 1992, 1993, 1994  Ian R. Searle
    Copyright (C) 2015, 2016 M. Kostrun

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
#include "mathl.h"
#include "symbol.h"
#include "mem.h"
#include "listnode.h"
#include "btree.h"
#include "bltin.h"
#include "bltin1.h"
#include "util.h"
#include "scan.h"
#include "print.h"
#include "class.h"
#include "mds.h"
#include "mdr.h"
#include "mdrf1.h"
#include "rfileio.h"
#include "mdc.h"
#include "mdr_mdc.h"
#include "mdrf1.h"
#include "mdcf1.h"
#include "mdcf2.h"
#include "ent.h"
#include "mde.h"

#include <stdio.h>
#include <fcntl.h>
#include <errno.h>
#include <string.h>
#include <ctype.h>
#include <errno.h>
#include <math.h>
#include <time.h>
#include <pty.h>

#define _SLPTIME_ 10

#include "rfileio.h"
#include "rlab_solver_parameters_names.h"
#include "rlab_macros.h"
#include "rlab_macros_code.h"

extern char string_buff[MAX_STRING_BUFF];

/*
 * The getline function will be a simple C-scanner, that reads a line
 * and splits into tokens. It would be nice in the future to be able to
 * change the token delimiter (say to a comma), or be able to read fixed
 * field data (ala FORTRAN).
 *
 * The token types getline will recognize are NUMBER and STRING.
 * For now getline will skip blanks and tabs.
 * Everything other than NUMBER will be a string.
 * Getline should recognize quotes (") and get everthing between them.
 */

static int next (FILE * fn, char *chr);
static int getline_scanner (FILE * fn);
static int collect_decimal (FILE * fn, int c);
static int collect_string (FILE * fn);

/*
 * Union for getline() scanner
 */

static union _uscan
{
  double dval;
  char *str;
}
uscan;

#define NUMBER 1000
#define STRING 1001
#define UNEXPECTED  -1
#define BAD_DECIMAL -2

char scan_code[256] = {
  0, 34, 34, 34, 34, 34, 34, 34, 34, 1, 2, 1, 1, 1, 34, 34,
  34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34,
  1, 27, 23, 25, 33, 15, 10, 34, 17, 18, 13, 11, 30, 12, 31, 14,
  22, 22, 22, 22, 22, 22, 22, 22, 22, 22, 8, 3, 28, 26, 29, 7,
  34, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 19, 24, 20, 16, 21,
  34, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 5, 9, 6, 32, 34,
  34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34,
  34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34,
  34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34,
  34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34,
  34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34,
  34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34,
  34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34,
  34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34
};

//
// RLaB builtin function interface for getline().
//
#undef  THIS_SOLVER
#define THIS_SOLVER "getline"
Ent *
Getline(int nargs, Datum args[])
{
  char istr[32], *retv;
  char *fname=0;
  int N = MAX_STRING_BUFF;
  int i, token_type;
  Ent *e1=0, *e2=0;
  FILE *fn=0;
  FILE *rlab_stderr = (!RLAB_STDERR_DS) ? stderr : RLAB_STDERR_DS;
  Btree *btree=0;

  /* Check n_args */
  if (nargs != 1 && nargs != 2)
    rerror (RLAB_ERROR_ONE_OR_TWO_ARG_REQUIRED);

  // get file name
  RLABCODE_PROCESS_ARG1_S(THIS_SOLVER,0,e1,fname,1);

  if ((fn = get_file_ds (fname, "r", 0)) == 0)
  {
    ent_Clean (e1);
    fprintf (rlab_stderr, "%s, cannot open for read\n", fname);
    rerror (THIS_SOLVER ": " RLAB_ERROR_CANNOT_OPEN_FILE_FOR_READ);
  }

  if (nargs == 2)
  {
    /*
     * The second arg is either 0, in which case we read
     * a NLINE length line, or N in which case we read a N
     * length line (MAX)
     */

    e2 = bltin_get_ent (args[1]);
    if (ent_type(e2)!=MATRIX_DENSE_REAL)
      rerror (THIS_SOLVER ": " RLAB_ERROR_ARG2_MDR_INTEGER);
    N = (int) class_double (e2);
    if (N > MAX_STRING_BUFF || N <=0)
      N = MAX_STRING_BUFF;

    retv = fgets (string_buff, N, fn);

    ent_Clean(e1);
    ent_Clean(e2);

    if (!retv)
      return ent_Create_Rlab_Failure();

    return ent_Create_Rlab_String(string_buff);
  }

  //
  // "Normal" getline usage...
  // Call scanner until newline is returned.
  //

  // Create the list
  btree = btree_Create ();

  i = 1;
  while ((token_type = getline_scanner (fn)) != '\n' && token_type != 0)
  {
    sprintf (istr, "%d", i++);
    if (token_type < 0)		/* Scanner error */
      warning_1 ("scanner error, index skipped");
    else if (token_type == NUMBER)
      install (btree, istr, ent_Create_Rlab_Double(uscan.dval));
    else if (token_type == STRING)
      install (btree, istr, ent_Create_Rlab_String(uscan.str));
  }

  /* Now check for blank line */

  if (i == 1 && token_type == '\n')
  {
    /* Add a NULL string to the list */
    install (btree, "1", ent_Create_Rlab_String(""));
  }

  ent_Clean(e1);
  ent_Clean(e2);

  return ent_Assign_Rlab_BTREE(btree);
}

/*
 * Scanner code.
 */

static int
getline_scanner (FILE *fn)
{
  char c;

reswitch:

  switch (next (fn, &c))
  {
    case 0:     /* EOF */
      return (0);

    case SC_SPACE:
      goto reswitch;

    case SC_NL:
      return ('\n');

    case SC_MINUS:    /* Numbers */
    case SC_PLUS:
    case SC_DIGIT:
    case SC_DOT:
      return (collect_decimal (fn, c));

    case SC_DQUOTE:   /* Double quoted strings */
      return (collect_string (fn));

    default:      /* Anything else is a whitespace delim string */
    {
      char *p = (char *) string_buff + 1;
      string_buff[0] = c;

      while (next (fn, p++) != 0 &&
        scan_code[(int) p[-1]] != SC_SPACE &&
        scan_code[(int) p[-1]] != SC_NL)
        ;

      ungetc (*--p, fn);
      *p = '\0';
      uscan.str = cpstr (string_buff);
      return (STRING);
    }

  }
}

/* **************************************************************
 * Scanner support functions.
 * ************************************************************** */

static int
    next (FILE *fn, char *chr)
{
  int c;
  c = getc (fn);
  if (c == EOF)
  {
    *chr = 0;
    return (0);
  }
  else
  {
    *chr = (char) c;
    return (scan_code[c]);
  }
}

/*
 * Collect a decimal constant in string_buff.
 * If the number turns out to be a string, then
 * return a string (not an error).
 */

static int
    collect_decimal (FILE *fn, int c)
{
  double d;
  char *p = (char *) string_buff + 1;
  char *endp;

  string_buff[0] = c;

  if (c == '+' || c == '-')
    next (fn, p++);

  if (p[-1] == '.')
  {
    if (next (fn, p++) == 0)
      return (0);
    else if (scan_code[(int) p[-1]] != SC_DIGIT)
      goto string;
  }
  else
  {
    while (next (fn, p++) == SC_DIGIT)
      ;
    if (scan_code[(int) p[-1]] == 0)
      goto finish;
    else if (p[-1] != '.')
    {
      ungetc (p[-1], fn);
      p--;
    }
  }

  /* get rest of digits after decimal point */
  while (next (fn, p++) == SC_DIGIT)
    ;
  if (scan_code[(int) p[-1]] == 0)
    goto finish;

  /* check for exponent */
  if (p[-1] != 'e' && p[-1] != 'E')
  {
    if (scan_code[(int) p[-1]] != SC_SPACE && scan_code[(int) p[-1]] != SC_NL)
      goto string;
    ungetc (p[-1], fn);
    *--p = 0;
  }
  else if (next (fn, p++) != SC_DIGIT && p[-1] != '-' && p[-1] != '+')
  {
    goto string;
  }
  else if (scan_code[(int) p[-1]] == 0)
  {
    goto finish;
  }
  else
  {
    /* get the rest of the exponent */
    while (next (fn, p++) == SC_DIGIT)
      ;
    ungetc (p[-1], fn);
    *--p = 0;
  }

finish:

  errno = 0;			/* check for overflow/underflow */
  d = strtod (string_buff, (char **) &endp);

  if (errno)
    rerror ("decimal over/under flow problem");

  if (endp < (p - 1))
    goto string;

  uscan.dval = d;
  return (NUMBER);

string:

  /* Not a number, collect a string. */
  while (next (fn, p++) != 0 &&
	 scan_code[(int) p[-1]] != SC_SPACE && scan_code[(int) p[-1]] != SC_NL)
    ;
  ungetc (*--p, fn);
  *p = '\0';
  uscan.str = cpstr (string_buff);
  return (STRING);
}

/*
 * Collect a doubly quoted string in string_buff
 * Stuff the result into the scanner union, and
 * return the token type.
 */

extern char *rm_escape (char *s);

static int
collect_string (FILE *fn)
{
  char *p = (char *) string_buff;
  char c;
  int e_flag = 0;		/* on if have an escape char */

  while (1)
  {
    switch (next (fn, p++))
    {
      case 0:     /* EOF, unterminated string */
        rerror ("runaway string constant");

        case SC_DQUOTE:   /* done */
          *--p = 0;
          goto out;

      case SC_NL:
        p[-1] = 0;
        /* fall thru */

      case SC_ESCAPE:
        if (next (fn, &c) == '\n')
        {
          p--;
        }
        else
        {
          if (c == 0)
          {
            rerror ("runaway string constant");
          }
          else
          {
            *p++ = c;
            e_flag = 1;
          }
        }
        break;

      default:
        break;
    }
  }

out:

  uscan.str = cpstr (e_flag ? rm_escape (string_buff) : string_buff);
  return STRING;
}

//
// rlabplus (C) Marijan Kostrun, 2005-2015
//
#undef THIS_SOLVER
#define THIS_SOLVER "strlen"
Ent * Strlen (int nargs, Datum args[])
{
  Ent *e1=0;
  MDS *m=0;
  MDR *ms=0;
  int i;

  // Check number of arguments

  if (nargs != 1)
  {
    printf (THIS_SOLVER ": " RLAB_ERROR_ONE_ARG_REQUIRED "\n");
    goto _exit_strlen;
  }

  /* Get the first argument */
  e1 = bltin_get_ent (args[0]);

  if (ent_type(e1) == MATRIX_DENSE_STRING)
  {
    m  = ent_data (e1);
    ms = mdr_Create (MNR (m), MNC (m));
    for (i=0; i<SIZE(m); i++)
      MdrV0 (ms, i) = (double) (isvalidstring (MdsV0 (m, i)));
  }
  else
  {
    ms = mdr_CreateScalar (-1.0);
  }

_exit_strlen:

  ent_Clean (e1);
  return ent_Assign_Rlab_MDR(ms);
}

int strindex (char *s, char *c, int lenc)
{
  int i, j, agree;
  // need the length of the string pointed by *c
  for (i=0; s[i] != 0; i++)
  {
    agree=0;
    for (j=0; j<lenc; j++)
    {
      if (s[i + j] == 0)
        return -1;
      if (s[i + j] != c[j])
        break;
      agree++;
    }
    if (agree == lenc)
      return i;
  }
  return -1;
}


int
string_substitute (char *r, int lr, char *s, int ls, char *c, int lc, char *newc)
{
  int i = strlen (c), ioff, j, loldc, k = 0;
  char *oldc = string_buff + (MAX_STRING_BUFF>>1);

  strcpy (oldc, c);
  strcpy (newc, c);
  loldc = strlen (oldc);
  ioff = 0;
  j = strindex (oldc, s, ls);
  while ((j > -1) && (loldc > 0))
  {
    k++;
    for (i = 0; i < j; i++)
      newc[i + ioff] = oldc[i];
    for (i = 0; i < lr; i++)
      newc[j + i + ioff] = r[i];
    ioff += j + lr;
    newc[ioff] = '\0';
    for (i = 0; oldc[i + j + ls] != '\0'; i++)
    {
      oldc[i] = oldc[i + j + ls];
      oldc[i + j + ls] = '\0';
    }
    oldc[i] = '\0';
    j = strindex (oldc, s, ls);
    if (j == -1)
      strcat (newc, oldc);
    loldc = strlen (oldc);
  }
  oldc[i] = '\0';

  return k;
}



#undef THIS_SOLVER
#define THIS_SOLVER "substr"
Ent *
ent_string_substr (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0;
  int i, j, k, nr, nc, ls, li;
  MDR *in=0;
  MDS *st=0, *re=0;

  if (nargs != 2)
  {
    fprintf (stderr,
             THIS_SOLVER ": Extracting the indexed characters from a string matrix.\n");
    fprintf (stderr,
             THIS_SOLVER ": Format:\n");
    fprintf (stderr,
             THIS_SOLVER ":   substr(s, idxs),\n");
    fprintf (stderr,
             THIS_SOLVER ": where 's' is a string matrix and 'idxs' is a row-matrix of indices\n");
    fprintf (stderr,
             THIS_SOLVER ": of the characters in a string entry to be extracted.\n");
    rerror (RLAB_ERROR_TWO_ARG_REQUIRED);
  }

  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) != MATRIX_DENSE_STRING)
    rerror (RLAB_ERROR_ARG1_MDS_MATRIX);
  st = ent_data (e1);
  if (SIZE(st)<1)
    goto _exit_substr;
  nr = MNR (st);
  nc = MNC (st);

  e2 = bltin_get_ent (args[1]);
  if (ent_type (e2) != MATRIX_DENSE_REAL)
    rerror (RLAB_ERROR_ARG2_MDR_INTEGER_VECTOR);
  in = class_matrix_real (e2);
  li = SIZE(in);
  if (li<1)
    goto _exit_substr;

  re = mds_Create (nr, nc);

  for (j=0; j < nr * nc; j++)
  {
    ls = isvalidstring(MdsV0 (st, j));

    for (i=0; i<li; i++)
    {
      string_buff[i] = 0;
      k = mdiV0(in, i);
      if (k < 1 || k>ls)
        break;
      if (k <= ls)
        string_buff[i] = MdsV0(st,j)[k-1];
    }
    string_buff[i] = 0;
    MdsV0(re, j) = cpstr (string_buff);
  }

_exit_substr:

  ent_Clean (e1);
  ent_Clean (e2);

  return ent_Assign_Rlab_MDS(re);
}

Ent *
ent_spinner (int nargs, Datum args[])
{
  Ent *e1=0;
  int idiv=1, j=0;
  static char s[] = {'|' , '/', '-', '\\' };
  static int is=0, j_prev=-1;

  if (nargs)
  {
    e1 = bltin_get_ent (args[0]);
    if (ent_type (e1) == MATRIX_DENSE_REAL)
    {
      idiv = class_int(e1);
    }
    if (idiv < 1)
      idiv = 1;
  }

  ++is;
  if (idiv == 1)
  {
    j = is % 4;
  }
  else
  {
    j = ((int) (is/idiv)) % 4;
  }

  if (j != j_prev)
  {
    fprintf(stdout, "%c\b", s[j]);
    fflush (stdout);

    j_prev = j;
  }

  ent_Clean(e1);

  return ent_Create_Rlab_Success();
}

Ent *
ent_smiley (int nargs, Datum args[])
{
  static char *s[] = { ":)" , ";)", ";|", "=|", "=(", ";(", ":(" };
  static int is=0;

  is++;
  is = is % 7;

  fprintf(stdout,"%2s\b\b", s[is]);
  fflush (stdout);

  return ent_Create_Rlab_Success();
}

#undef THIS_SOLVER
#define THIS_SOLVER "capitalize"
Ent *
ent_string_capitalize (int nargs, Datum args[])
{
  Ent *M=0;
  int i, k;
  MDS *s=0, *m=0;
  char *s1;
  if (nargs != 1)
  {
    fprintf (stderr,
             THIS_SOLVER ": Capitalizes the first letter of all words in the string.\n");
    fprintf (stderr,
             THIS_SOLVER ": Format:\n");
    fprintf (stderr,
             THIS_SOLVER ":   capitalize(s),\n");
    fprintf (stderr,
             THIS_SOLVER ": where 's' is a string matrix.\n");
    rerror (RLAB_ERROR_ONE_ARG_REQUIRED);
  }
  M = bltin_get_ent (args[0]);
  if (ent_type (M) != MATRIX_DENSE_STRING)
    rerror (THIS_SOLVER ": " RLAB_ERROR_ARG1_MDS_MATRIX);
  m = ent_data (M);

  s = mds_Copy (m);
  for (i = 0; i < MNR(s) * MNC(s); i++)
  {
    s1 = MdsV0 (s,i);

    for (k = 0; k < strlen (s1); k++)
    {
      // first letter goes uppercase
      if (k == 0)
      {
        if (s1[k] > 96 && s1[k] < 123)
        {
          s1[k] -= 32;
        }
        continue;
      }

      if (s1[k-1] == ' ')
      {
        if (s1[k] > 96 && s1[k] < 123)
        {
          s1[k] -= 32;
        }
        continue;
      }

      // lowercase all the others
      if (s1[k] > 64 && s1[k] < 91)
      {
        s1[k] += 32;
      }
    }
  }

  ent_Clean (M);

  return ent_Assign_Rlab_MDS(s);
}

#undef THIS_SOLVER
#define THIS_SOLVER "tolower"
Ent *
ent_string_tolower (int nargs, Datum args[])
{
  Ent *M=0;
  int i, k;
  MDS *s=0, *m=0;
  char *s1;
  if (nargs != 1)
  {
    fprintf (stderr,
             THIS_SOLVER ": Converts upper-case letters of alphabet to lower-case.\n");
    fprintf (stderr,
             THIS_SOLVER ": Format:\n");
    fprintf (stderr,
             THIS_SOLVER ":   tolower(s),\n");
    fprintf (stderr,
             THIS_SOLVER ": where 's' is a string matrix.\n");
    rerror (RLAB_ERROR_ONE_ARG_REQUIRED);
  }
  M = bltin_get_ent (args[0]);
  if (ent_type (M) != MATRIX_DENSE_STRING)
    rerror (THIS_SOLVER ": " RLAB_ERROR_ARG1_MDS_MATRIX);
  m = ent_data (M);

  s = mds_Copy (m);
  for (i = 0; i < MNR(s) * MNC(s); i++)
  {
    s1 = MdsV0 (s,i);
    for (k = 0; k < strlen (s1); k++)
      if (s1[k] > 64 && s1[k] < 91)
        s1[k] += 32;
  }

  ent_Clean (M);

  return ent_Assign_Rlab_MDS(s);
}

#undef THIS_SOLVER
#define THIS_SOLVER "toupper"
Ent *
ent_string_toupper (int nargs, Datum args[])
{
  Ent *e1=0;
  int i, k;
  MDS *s=0, *m=0;
  char *s1=0;

  if (nargs != 1)
  {
    fprintf (stderr,
             THIS_SOLVER ": Converts lower-case letters of alphabet to upper-case.\n");
    fprintf (stderr,
             THIS_SOLVER ": Format:\n");
    fprintf (stderr,
             THIS_SOLVER ":   toupper(s),\n");
    fprintf (stderr,
             THIS_SOLVER ": where 's' is a string matrix.\n");
    rerror (RLAB_ERROR_ONE_ARG_REQUIRED);
  }

  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) != MATRIX_DENSE_STRING)
    rerror (THIS_SOLVER ": " RLAB_ERROR_ARG1_MDS_MATRIX);
  m = ent_data (e1);

  s = mds_Copy (m);
  for (i = 0; i < MNR(s) * MNC(s); i++)
  {
    s1 = MdsV0 (s, i);
    for (k = 0; k < strlen (MdsV0 (s, i)); k++)
    {
      if (s1[k] > 96 && s1[k] < 123)
      {
        s1[k] -= 32;
      }
    }
  }
  ent_Clean (e1);

  return ent_Assign_Rlab_MDS(s);
}


#undef THIS_SOLVER
#define THIS_SOLVER "num2str"
Ent *
ent_string_text (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *e3=0, *e4=0, *rent;
  char fmt[30] = { "\0" }, pfmt[30] = { "\0" };
  char onenum [128] = {'\0'};
  char p1i[32] = {'+', '1', 'i', '*'};
  char n1i[32] = {'-', '1', 'i', '*'};
  char *sep=0, *ci=0;
  double ds;
  int i, j, nr=0, nc=0, f_nr = 0, f_nc = 0, is, k, k2;
  MDR *x=0;
  MDC *z=0;
  MDS *s=0, *f=0;

  if (nargs < 1 || nargs > 5)
    rerror (THIS_SOLVER ": " RLAB_ERROR_ONE_TO_FIVE_ARG_REQUIRED "\n");

  //
  // x,z, dense real/complex matrix, data
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) == MATRIX_DENSE_STRING)
  {
    //
    // nothing to do here and no need to raise an error.
    //
    rent = ent_Copy (e1);
    ent_Clean(e1);
    return rent;
  }

  //
  // f, string matrix, format of the data
  //
  if (nargs > 1)
  {
    e2 = bltin_get_ent (args[1]);
    if (ent_type (e2) == MATRIX_DENSE_STRING)
    {
      f = ent_data (e2);
      if (SIZE(f)< 1)
        rerror (THIS_SOLVER ": improper second argument");
      f_nr = MNR (f);
      f_nc = MNC (f);
    }
  }

  //
  // sep, string matrix, separator field
  //
  if (nargs > 2)
  {
    e3 = bltin_get_ent (args[2]);
    if (ent_type (e3) == MATRIX_DENSE_STRING)
      sep = class_char_pointer (e3);
    if (isvalidstring(sep)<0)
      sep = 0;
  }

  //
  // I, complex unity
  //
  if (nargs > 3)
  {
    e4 = bltin_get_ent (args[3]);
    if (ent_type (e4) == MATRIX_DENSE_STRING)
    {
      ci = class_char_pointer (e4);
      if (isvalidstring(ci)<1)
        ci = 0;

      if (ci)
      {
        for (k=0; k<=MIN(strlen(ci),62); k++)
        {
          p1i[k] = ci[k];
          n1i[k] = ci[k];
        }
          // prepend the sign to complex unity
        if (p1i[0]=='+')
          n1i[0] = '-';
        if (n1i[0]=='-')
          p1i[0] = '+';
          // append the sign to complex unity if they are not prepended
        if (n1i[0]!='-')
        {
          k = strlen(n1i);
          n1i[k]   = '-';
          n1i[k+1] = 0;
        }
      }
    }
  }

  if (ent_type (e1) == MATRIX_DENSE_REAL)
  {
    x = ent_data (e1);
    nr = MNR(x);
    nc = MNC(x);
  }
  else if (ent_type (e1) == MATRIX_DENSE_COMPLEX)
  {
    z = ent_data (e1);
    nr = MNR(z);
    nc = MNR(z);
  }
  else
    rerror ("text: improper first argument");

  // nothing really to do
  if (nr * nc == 0)
  {
    ent_Clean (e1);
    ent_Clean (e2);
    ent_Clean (e3);
    ent_Clean (e4);

    return ent_Assign_Rlab_MDS(s);
  }




  if (x)
  {
    // default format
    if (!f)
    {
      fmt[0] = '%';
      if(x->type == RLAB_TYPE_INT32)
        fmt[1] = 'i';
      else
        fmt[1] = 'g';
      fmt[2] = '\0';
      sprintf (pfmt, "%s", fmt);
    }

    // separator exists
    if (sep)
    {
      s = mds_Create (nr, 1);
      for (i = 0; i < nr; i++)
      {
        string_buff[0] = 0;

        //
        // first column
        //
        if (f)
        {
          char *cf = Mds0(f, MIN(i,f_nr-1), 0);
          if (cf[0]=='+' || cf[0]=='-')
          {
            sprintf(fmt,  "%s", &cf[1]);
            sprintf(pfmt,"+%s", &cf[1]);
          }
          else
          {
            sprintf(fmt,  "%s", &cf[0]);
            sprintf(pfmt, "%s", &cf[0]);
          }
        }

        switch (x->type)
        {
          case RLAB_TYPE_INT32:
            if (Mdi0(x,i,0) >=0)
              sprintf(string_buff, pfmt, Mdi0(x,i,0));
            else
              sprintf(string_buff,  fmt, Mdi0(x,i,0));
            break;

          case RLAB_TYPE_DOUBLE:
            if (isnand(Mdr0(x,i,0)))
            {
              sprintf(string_buff, RLAB_SPRINTF_NAN);
            }
            else
            {
              if (Mdr0(x,i,0) >=0)
                sprintf(string_buff, pfmt, Mdr0(x,i,0));
              else
                sprintf(string_buff,  fmt, Mdr0(x,i,0));
            }
            break;

          default:
            break;
        }

        //
        // all the other columns
        //
        for (j = 1;  j < nc; j++)
        {
          // format
          if (f)
          {
            char * cf = Mds0(f, MIN(i,f_nr-1), MIN(j,f_nc-1));
            if (cf[0]=='+' || cf[0]=='-')
            {
              sprintf(fmt,  "%s", &cf[1]);
              sprintf(pfmt,"+%s", &cf[1]);
            }
            else
            {
              sprintf(fmt,  "%s", &cf[0]);
              sprintf(pfmt, "%s", &cf[0]);
            }
          }

          // copy separator to 'oneline'
          k = strlen(string_buff);
          if (isvalidstring(sep)>0)
            for (k2=0; k2<=isvalidstring(sep); k2++)
              string_buff[k + k2] = sep[k2];

          // copy number to 'oneline'
          switch (x->type)
          {
            case RLAB_TYPE_INT32:
              is = Mdi0 (x, i, j);
              if (is >= 0)
                sprintf (onenum, pfmt, is);
              else
                sprintf (onenum,  fmt, is);
              break;

            case RLAB_TYPE_DOUBLE:
              ds = Mdr0 (x, i, j);
              if (isnand(ds))
              {
                sprintf(onenum, RLAB_SPRINTF_NAN);
              }
              else
              {
                if (ds >= 0)
                  sprintf (onenum, pfmt, ds);
                else
                  sprintf (onenum,  fmt, ds);
              }
              break;

            default:
              break;
          }

          k = strlen(string_buff);
          for (k2 = 0; k2 <= strlen(onenum); k2 ++)
            string_buff[k + k2] = onenum[k2];
          onenum[0] = 0;
        }
        MdsV0(s, i) = cpstr (string_buff);
      }
    }
    else
    {
      s = mds_Create (nr, nc);
      for (i = 0; i < nr; i++)
        for (j = 0; j < nc; j++)
      {
        // format
        if (f)
        {
          char *cf = Mds0(f, MIN(i,f_nr-1), MIN(j,f_nc-1));
          if (cf[0]=='+' || cf[0]=='-')
          {
            sprintf(fmt,  "%s", &cf[1]);
            sprintf(pfmt,"+%s", &cf[1]);
          }
          else
          {
            sprintf(fmt,  "%s", &cf[0]);
            sprintf(pfmt, "%s", &cf[0]);
          }
        }

        // copy number to 'oneline'
        switch (x->type)
        {
          case RLAB_TYPE_INT32:
            is = Mdi0 (x, i, j);
            if (is >= 0)
              sprintf (onenum, pfmt, is);
            else
              sprintf (onenum,  fmt, is);
            break;

          case RLAB_TYPE_DOUBLE:
            ds = Mdr0 (x, i, j);
            if (isnand(ds))
            {
              sprintf(onenum, RLAB_SPRINTF_NAN);
            }
            else
            {
              if (ds >= 0)
                sprintf (onenum, pfmt, ds);
              else
                sprintf (onenum,  fmt, ds);
            }
            break;

          default:
            break;
        }

        Mds0(s, i, j) = cpstr( onenum );
        onenum[0] = 0;
      }
    }
  }
  else if (z)
  {
    // default format
    if (!f)
    {
      fmt[0] = '%';
      fmt[1] = 'g';
      fmt[2] = '\0';
      sprintf (pfmt, "%s", fmt);
    }

    // separator exists
    if (sep)
    {
      s = mds_Create (nr, 1);
      for (i = 0; i < nr; i++)
      {
        //
        // first column
        //
        if (f)
        {
          char *cf = Mds0(f, MIN(i,f_nr-1), 0);
          if (cf[0]=='+' || cf[0]=='-')
          {
            sprintf(fmt,  "%s", &cf[1]);
            sprintf(pfmt,"+%s", &cf[1]);
          }
          else
          {
            sprintf(fmt,  "%s", &cf[0]);
            sprintf(pfmt, "%s", &cf[0]);
          }
        }

        // if NaN, just print them
        if (isnand(Mdc0r(z,i,0)) || isnand(Mdc0i(z,i,0)))
        {
          sprintf(string_buff,RLAB_SPRINTF_NAN);
        }
        else
        {
          // print Re(z)
          if (Mdc0r(z,i,0) >=0)
            sprintf(string_buff, pfmt, Mdc0r(z,i,0));
          else
            sprintf(string_buff,  fmt, Mdc0r(z,i,0));

          // print "+1i*" or "-1i*" to 'oneline'
          if (Mdc0i(z,i,0)>0)
          {
          // copy  '+1i*' to 'oneline'
            k = strlen(string_buff);
            for (k2 = 0; k2 <= strlen(p1i); k2 ++)
              string_buff[k + k2] = p1i[k2];
          }
          else
          {
          // copy  '-1i*' to 'oneline'
            k = strlen(string_buff);
            for (k2 = 0; k2 <= strlen(n1i); k2 ++)
              string_buff[k + k2] = n1i[k2];
          }

          // print Abs(Im(z)) to 'onenum' and then append it to 'oneline'
          ds = Mdc0i (z, i, 0) >= 0 ? Mdc0i (z, i, 0) : -Mdc0i (z, i, 0);
          sprintf (onenum,  fmt, ds);
          k = strlen(string_buff);
          for (k2 = 0; k2 <= strlen(onenum); k2 ++)
            string_buff[k + k2] = onenum[k2];
          onenum[0] = 0;
        }

        //
        // all the other columns
        //
        for (j = 1;  j < nc; j++)
        {
          // format
          if (f)
          {
            char * cf = Mds0(f, MIN(i,f_nr-1), MIN(j,f_nc-1));
            if (cf[0]=='+' || cf[0]=='-')
            {
              sprintf(fmt,  "%s", &cf[1]);
              sprintf(pfmt,"+%s", &cf[1]);
            }
            else
            {
              sprintf(fmt,  "%s", &cf[0]);
              sprintf(pfmt, "%s", &cf[0]);
            }
          }

          // copy separator to 'oneline'
          k = strlen(string_buff);
          for (k2 = 0; k2 <= strlen(sep); k2 ++)
            string_buff[k + k2] = sep[k2];

          if (isnand(Mdc0r(z,i,j)) || isnand(Mdc0i(z,i,j)))
          {
            sprintf(onenum, RLAB_SPRINTF_NAN);
          }
          else
          {
            // print Re(z) to 'onenum' and then to 'oneline'
            if (Mdc0r(z,i,j) >=0)
              sprintf(onenum, pfmt, Mdc0r(z,i,j));
            else
              sprintf(onenum,  fmt, Mdc0r(z,i,j));
            k = strlen(string_buff);
            for (k2 = 0; k2 <= strlen(onenum); k2 ++)
              string_buff[k + k2] = onenum[k2];
            onenum[0] = 0;

            // print "+1i*" or "-1i*" to 'oneline'
            if (Mdc0i(z,i,j)>=0)
            {
            // copy  '+1i*' to 'oneline'
              k = strlen(string_buff);
              for (k2 = 0; k2 <= strlen(p1i); k2 ++)
                string_buff[k + k2] = p1i[k2];
            }
            else
            {
            // copy  '-1i*' to 'oneline'
              k = strlen(string_buff);
              for (k2 = 0; k2 <= strlen(n1i); k2 ++)
                string_buff[k + k2] = n1i[k2];
            }

            // print Abs(Im(z)) to 'onenum' and then append it to 'oneline'
            ds = Mdc0i (z, i, j) > 0 ? Mdc0i (z, i, j) : -Mdc0i (z, i, j);
            sprintf (onenum,  fmt, ds);
          }

          k = strlen(string_buff);
          for (k2 = 0; k2 <= strlen(onenum); k2 ++)
            string_buff[k + k2] = onenum[k2];
          onenum[0] = 0;
        }
        MdsV0(s, i) = cpstr (string_buff);
      }
    }
    else
    {
      s = mds_Create (nr, nc);
      for (i = 0; i < nr; i++)
        for (j = 0; j < nc; j++)
      {
        // format
        if (f)
        {
          char * cf = Mds0(f, MIN(i,f_nr-1), MIN(j,f_nc-1));
          if (cf[0]=='+' || cf[0]=='-')
          {
            sprintf(fmt,  "%s", &cf[1]);
            sprintf(pfmt,"+%s", &cf[1]);
          }
          else
          {
            sprintf(fmt,  "%s", &cf[0]);
            sprintf(pfmt, "%s", &cf[0]);
          }
        }

        // print Re(z) to 'oneline'
        if (Mdc0r(z,i,j) >=0)
          sprintf(string_buff, pfmt, Mdc0r(z,i,j));
        else
          sprintf(string_buff,  fmt, Mdc0r(z,i,j));

        // print "+1i*" or "-1i*" to 'oneline'
        if (Mdc0i(z,i,j)>=0)
        {
          // copy  '+1i*' to 'oneline'
          k = strlen(string_buff);
          for (k2 = 0; k2 <= strlen(p1i); k2 ++)
            string_buff[k + k2] = p1i[k2];
        }
        else
        {
          // copy  '-1i*' to 'oneline'
          k = strlen(string_buff);
          for (k2 = 0; k2 <= strlen(n1i); k2 ++)
            string_buff[k + k2] = n1i[k2];
        }

        // print Abs(Im(z)) to 'onenum' and then append it to 'oneline'
        ds = Mdc0i (z, i, j) > 0 ? Mdc0i (z, i, j) : -Mdc0i (z, i, j);
        sprintf (onenum,  fmt, ds);
        k = strlen(string_buff);
        for (k2 = 0; k2 <= strlen(onenum); k2 ++)
          string_buff[k + k2] = onenum[k2];
        onenum[0] = 0;

        Mds0(s, i, j) = cpstr( string_buff );
      }
    }
  }
  else
    rerror("text: terrible internal error");

  ent_Clean (e1);
  ent_Clean (e2);
  ent_Clean (e3);
  ent_Clean (e4);
  return ent_Assign_Rlab_MDS(s);
}

#undef  THIS_SOLVER
#define THIS_SOLVER "blank"
Ent *
ent_string_blank (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *e3=0;
  int nr=0, nc=0, i;
  char *s0=0;
  MDS *s=0;
  FILE *rlab_stderr = (!RLAB_STDERR_DS) ? stderr : RLAB_STDERR_DS;

  if (nargs==0)
    goto _exit_blank;

  if (nargs < 1 || nargs > 3)
  {
    fprintf (rlab_stderr, THIS_SOLVER ": Create a string matrix and set its initial value.\n");
    fprintf (rlab_stderr, THIS_SOLVER ": Format:\n");
    fprintf (rlab_stderr, THIS_SOLVER ":   blank(nrow,ncol/,s/),\n");
    fprintf (rlab_stderr, THIS_SOLVER ": where 'nrow' and 'ncol' determine the size of the string matrix, while\n");
    fprintf (rlab_stderr, THIS_SOLVER ": 's' is its content (default \"\").\n");
    rerror ("requires at least one argument!");
  }

  // determine the first entity
  e1 = bltin_get_ent (args[0]);
  if ((ent_type(e1) == MATRIX_DENSE_REAL) || (ent_type(e1) == MATRIX_DENSE_COMPLEX) || (ent_type(e1) == MATRIX_DENSE_STRING))
  {
    MD * x1 = ent_data (e1);
    nr = MNR( x1 );
    nc = MNC( x1 );
  }
  else
  {
    fprintf (rlab_stderr, THIS_SOLVER ": blank() is not defined for %s\n", etd(e1));
    rerror  (THIS_SOLVER ": cannot create a dense string matrix!");
  }

  if (nargs >= 2)
  {
    e2 = bltin_get_ent (args[1]);
    if (ent_type (e2) == MATRIX_DENSE_STRING)
      s0 = class_char_pointer(e2);
    else if (ent_type(e2)==MATRIX_DENSE_REAL)
    {
      nr = (int) class_double (e1);
      nc = (int) class_double (e2);
    }
  }

  nr = nr * (nc>0);
  nc = nc * (nr>0);

  // third argument is the fill-up string
  if (nargs == 3)
  {
    e3 = bltin_get_ent (args[2]);
    if (ent_type (e3) == MATRIX_DENSE_STRING)
      s0 = class_char_pointer(e3);
  }

  s = mds_Create (nr, nc);
  if (s0)
  {
    for (i=0; i<nr*nc; i++)
      MdsV0 (s, i) = cpstr (s0);
  }
  else
  {
    for (i=0; i<nr*nc; i++)
      MdsV0 (s, i) = cpstr ("");
  }

_exit_blank:

  ent_Clean (e1);
  ent_Clean (e2);
  ent_Clean (e3);

  return ent_Assign_Rlab_MDS(s);
}


#undef THIS_SOLVER
#define THIS_SOLVER "char"
Ent *
ent_string_char (int nargs, Datum args[])
{
  Ent *e1=0;
  int i, j, nr, nc;
  char *s1 = string_buff;

  MDR *m=0;
  MDS *s=0;
  if (nargs != 1)
  {
    fprintf (stdout,
             THIS_SOLVER ": Create a string vector from an integer matrix.\n");
    fprintf (stdout,
             THIS_SOLVER ": Format:\n");
    fprintf (stdout,
             THIS_SOLVER ":  s = char( x )\n");
    fprintf (stdout,
             THIS_SOLVER ": 'x' is an integer matrix. The rows of 'x' are words in 's'.\n");
    rerror (RLAB_ERROR_ONE_ARG_REQUIRED);
  }

  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) != MATRIX_DENSE_REAL)
    rerror (RLAB_ERROR_ARG1_MDR_INTEGER_MATRIX);

  m = ent_data (e1);
  nr = MNR (m);
  nc = MNC (m);
  if (nc > MAX_STRING_BUFF)
    nc = MAX_STRING_BUFF;

  s = mds_Create (nr, 1);
  if (nr * nc > 0)
  {
    for (i = 0; i < nr; i++)
    {
      for (j=0; j< nc; j++)
        s1[j] = (char) mdi0 (m, i, j);
      s1[nc] = 0;
      MdsV0 (s, i) = cpstr (s1);
    }
  }

  ent_Clean (e1);
  return ent_Assign_Rlab_MDS(s);
}

#undef THIS_SOLVER
#define THIS_SOLVER "ascii"
Ent *
ent_string_ascii (int nargs, Datum args[])
{
  Ent *e1=0;
  int i, nc;
  MDR *w=0;
  char *s=0;

  if (nargs != 1)
  {
    fprintf (stdout,
             THIS_SOLVER ": Create a integer vector from a single string.\n");
    fprintf (stdout,
             THIS_SOLVER ":   x = ascii( s )\n");
    fprintf (stdout,
             THIS_SOLVER ": 's' is a string. 'x' is a row-vector containing ascii codes\n");
    fprintf (stdout,
             THIS_SOLVER ": of the characters in 's'.\n");
    rerror (RLAB_ERROR_ONE_ARG_REQUIRED);
  }

  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) != MATRIX_DENSE_STRING)
    rerror (RLAB_ERROR_ARG1_MDS_SCALAR);
  s = class_char_pointer ( e1 );

  if (s)
  {
    nc = strlen ( s );
    if (nc)
    {
      w  = mdi_Create (1, nc);
      for (i = 0; i < nc; i++)
        MdiV0(w,i) = s[i];
    }
  }

  ent_Clean (e1);

  return ent_Assign_Rlab_MDR (w);
}

#undef  THIS_SOLVER
#define THIS_SOLVER "ifelse"
Ent * ent_ifelse (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *e3=0, *e4=0;
  MDR *x1=0;
  MDR *r2=0, *r3=0, *r4=0, * w=0;
  MDC *c2=0, *c3=0, *c4=0, *cw=0;
  MDS *s2=0, *s3=0, *s4=0, *sw=0;

  int cr, cc, wr, wc, tr, tc, fr, fc;
  int i, j;

  if ((nargs != 3) && (nargs != 4))
  {
    printf (THIS_SOLVER ": Spreadsheet-If command.\n");
    printf (THIS_SOLVER ": Format:\n");
    printf (THIS_SOLVER ":   ifelse(cond,iftrue,iffalse),\n");
    printf (THIS_SOLVER ":   ifelse(cond,ifpos,ifzero,ifneg),\n");
    printf
        (THIS_SOLVER ": where 'cond' is a either a condition matrix containing 0's and 1's,\n");
    printf
        (THIS_SOLVER ": or real numbers in which sign(x) is used to pick the selection\n");
    printf
        (THIS_SOLVER ": 'iftrue' is the entry if the condition is true, and 'iffalse'\n");
    printf
        (THIS_SOLVER ": if condition is false. The function is matrix optimized.\n");
    rerror (THIS_SOLVER ": three or four arguments required!\n");
  }

  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) != MATRIX_DENSE_REAL)
    rerror ("ifelse: 'cond' has to be a matrix with entries 0 or 1!\n");

  e2 = bltin_get_ent (args[1]);
  if (    ent_type (e2) != MATRIX_DENSE_REAL
      &&  ent_type (e2) != MATRIX_DENSE_COMPLEX
      &&  ent_type (e2) != MATRIX_DENSE_STRING
     )
    rerror ("ifelse: 'iftrue/ifpos' has to be a dense matrix!\n");

  e3 = bltin_get_ent (args[2]);
  if (    ent_type (e3) != MATRIX_DENSE_REAL
      &&  ent_type (e3) != MATRIX_DENSE_COMPLEX
      &&  ent_type (e3) != MATRIX_DENSE_STRING
     )
    rerror ("ifelse: 'iffalse/ifneg' has to be a dense matrix!\n");

  if (nargs > 3)
  {
    e4 = bltin_get_ent (args[3]);
    if ( ent_type (e4) != MATRIX_DENSE_REAL
            &&  ent_type (e4) != MATRIX_DENSE_COMPLEX
            &&  ent_type (e4) != MATRIX_DENSE_STRING
       )
      rerror ("ifelse: 'ifneg' has to be a dense matrix!\n");
  }

  x1 = ent_data (e1);
  wr = cr = MNR(x1);
  wc = cc = MNC(x1);

  if (    ent_type (e2) == MATRIX_DENSE_REAL
      &&  ent_type (e3) == MATRIX_DENSE_REAL
     )
  {
    // both matrices 'iftrue' and 'iffalse' are real
    r2 = ent_data (e2);
    tr = MNR(r2);
    tc = MNC(r2);
    //
    r3 = ent_data (e3);
    fr = MNR(r3);
    fc = MNC(r3);

    if (e4)
    {
      if (ent_type(e4) == MATRIX_DENSE_REAL)
      {
        r4 = ent_data (e4);
        fr = MAX(fr, MNR(r4));
        fc = MAX(fc, MNC(r4));
      }
    }

    //
    wr = MAX( MAX(wr, tr), fr);
    wc = MAX( MAX(wc, tc), fc);

    if ( MD_TYPE_INT32(r2) &&  MD_TYPE_INT32(r3) )
    {
      if (!r4)
      {
        // ifelse on two integer matrices produces an integer matrix
        w = mdi_Create (wr, wc);

        // iterate
        for (i = 1; i <= wr; i++)
          for (j = 1; j <= wc; j++)
        {
          // reset it if necessary
          if ( mdr1_safe (x1, i, j) > 0 )
          {
            if (SIZE(r2))
              Mdi1 (w, i, j) = mdi1_safe (r2, i, j);
            else
              Mdi1 (w, i, j) = 0;
          }
          else
          {
            if (SIZE(r3))
              Mdi1 (w, i, j) = mdi1_safe (r3, i, j);
            else
              Mdi1 (w, i, j) = 0;
          }
        }
      }
      else if (MD_TYPE_INT32(r4))
      {
              // ifelse on two integer matrices produces an integer matrix
        w = mdi_Create (wr, wc);

        // iterate
        for (i = 1; i <= wr; i++)
          for (j = 1; j <= wc; j++)
        {
          // reset it if necessary
          if ( mdr1_safe (x1, i, j) > 0 )
          {
            if (SIZE(r2))
              Mdi1 (w, i, j) = mdi1_safe (r2, i, j);
            else
              Mdi1 (w, i, j) = 0;
          }
          else if ( mdr1_safe (x1, i, j) == 0 )
          {
            if (SIZE(r3))
              Mdi1 (w, i, j) = mdi1_safe (r3, i, j);
            else
              Mdi1 (w, i, j) = 0;
          }
          else
          {
            if (SIZE(r4))
              Mdi1 (w, i, j) = mdi1_safe (r4, i, j);
            else
              Mdi1 (w, i, j) = 0;
          }
        }
      }
    }
    else
    {
      // ifelse on two double matrices
      w = mdr_Create (wr, wc);

      // iterate
      if (!r4)
      {
        for (i = 1; i <= wr; i++)
          for (j = 1; j <= wc; j++)
        {
          if (mdr1_safe (x1, i, j) > 0)
          {
            if (SIZE(r2))
              Mdr1 (w, i, j) = mdr1_safe (r2, i, j);
            else
              Mdr1 (w, i, j) = create_nan();
          }
          else
          {
            if (SIZE(r3))
              Mdr1 (w, i, j) = mdr1_safe (r3, i, j);
            else
              Mdr1 (w, i, j) = create_nan();
          }
        }
      }
      else
      {
        for (i = 1; i <= wr; i++)
          for (j = 1; j <= wc; j++)
        {
          printf("x1[%i,%i] = %g\n", i, j, mdr1_safe (x1, i, j));
          if (mdr1_safe (x1, i, j) > 0)
          {
            if (SIZE(r2))
              Mdr1 (w, i, j) = mdr1_safe (r2, i, j);
            else
              Mdr1 (w, i, j) = create_nan();
          }
          else if (mdr1_safe (x1, i, j) == 0)
          {
            if (SIZE(r3))
              Mdr1 (w, i, j) = mdr1_safe (r3, i, j);
            else
              Mdr1 (w, i, j) = create_nan();
          }
          else
          {
            if (SIZE(r4))
              Mdr1 (w, i, j) = mdr1_safe (r4, i, j);
            else
              Mdr1 (w, i, j) = create_nan();
          }
        }
      }

    }

    ent_Clean (e1);
    ent_Clean (e2);
    ent_Clean (e3);
    ent_Clean (e4);
    return ent_Assign_Rlab_MDR(w);
  }
  else if (   (ent_type (e2) == MATRIX_DENSE_COMPLEX || ent_type (e2) == MATRIX_DENSE_REAL)
           || (ent_type (e3) == MATRIX_DENSE_COMPLEX || ent_type (e3) == MATRIX_DENSE_REAL)  )
  {
    //
    // one complex and one real matrix
    //
    // one of the matrices 'iftrue' or 'iffalse' is complex
    c2 = ent_data (e2);
    tr = MNR(c2);
    tc = MNC(c2);
    //
    c3 = ent_data (e3);
    fr = MNR(c3);
    fc = MNC(c3);

    if (e4)
    {
      if ((ent_type(e4) == MATRIX_DENSE_REAL) || (ent_type(e4) == MATRIX_DENSE_COMPLEX))
      {
        c4 = ent_data (e4);
        fr = MAX(fr, MNR(r4));
        fc = MAX(fc, MNC(r4));
      }
    }

    //
    wr = MAX(MAX (wr, tr), fr);
    wc = MAX(MAX (wc, tc), fc);

    cw = mdc_Create (wr, wc);

    if (!r4)
    {
      for (i = 1; i <= wr; i++)
        for (j = 1; j <= wc; j++)
      {
        if ( mdr1_safe (x1, i, j) > 0 )
        {
          if (SIZE(c2))
            Mdc1(cw, i, j) = mdc1_safe (c2, i, j);
          else
          {
            Mdc1r (cw, i, j) = create_nan();
            Mdc1i (cw, i, j) = 0.0;
          }
        }
        else
        {
          if (SIZE(c3))
          {
            Mdc1(cw, i, j) = mdc1_safe (c3, i, j);
          }
          else
          {
            Mdc1r (cw, i, j) = create_nan();
            Mdc1i (cw, i, j) = 0.0;
          }
        }
      }
    }
    else
    {
      for (i = 1; i <= wr; i++)
        for (j = 1; j <= wc; j++)
      {
        if ( mdr1_safe (x1, i, j) > 0 )
        {
          if (SIZE(c2))
            Mdc1(cw, i, j) = mdc1_safe (c2, i, j);
          else
          {
            Mdc1r (cw, i, j) = create_nan();
            Mdc1i (cw, i, j) = 0.0;
          }
        }
        else if ( mdr1_safe (x1, i, j) == 0 )
        {
          if (SIZE(c3))
          {
            Mdc1(cw, i, j) = mdc1_safe  (c3, i, j);
          }
          else
          {
            Mdc1r (cw, i, j) = create_nan();
            Mdc1i (cw, i, j) = 0.0;
          }
        }
        else
        {
          if (SIZE(c4))
          {
            Mdc1(cw, i, j) = mdc1_safe (c4, i, j);
          }
          else
          {
            Mdc1r (cw, i, j) = create_nan();
            Mdc1i (cw, i, j) = 0.0;
          }
        }
      }
    }

    ent_Clean (e1);
    ent_Clean (e2);
    ent_Clean (e3);
    ent_Clean (e4);
    return ent_Assign_Rlab_MDC(cw);
  }
  else if (    ent_type (e2) == MATRIX_DENSE_STRING
           &&  ent_type (e3) == MATRIX_DENSE_STRING
     )
  {
    //
    // two string matrices
    //
    s2 = ent_data (e2);
    tr = s2->nrow;
    tc = s2->ncol;
    //
    s3 = ent_data (e3);
    fr = s3->nrow;
    fc = s3->ncol;

    if (e4)
    {
      if (ent_type(e4) == MATRIX_DENSE_STRING)
      {
        s4 = ent_data (e4);
        fr = MAX(fr, MNR(r4));
        fc = MAX(fc, MNC(r4));
      }
    }

    //
    wr = MAX (MAX (wr, tr), fr);
    wc = MAX (MAX (wc, tc), fc);

    sw = mds_Create (wr, wc);
    if (!s4)
    {
      for (i = 1; i <= wr; i++)
        for (j = 1; j <= wc; j++)
      {
        Mds1 (sw, i, j) = cpstr (mds1_safe (s3, i, j));
        if (mdr1_safe (x1, i, j) > 0)
        {
          if (SIZE(s2))
            Mds1 (sw, i, j) = cpstr (mds1_safe (s2, i, j));
          else
            Mds1 (sw, i, j) = cpstr("UNDEF");
        }
        else
        {
          if (SIZE(s3))
            Mds1 (sw, i, j) = cpstr (mds1_safe (s3, i, j));
          else
            Mds1 (sw, i, j) = cpstr("UNDEF");
        }
      }
    }
    else
    {
      for (i = 1; i <= wr; i++)
        for (j = 1; j <= wc; j++)
      {
        Mds1 (sw, i, j) = cpstr (mds1_safe (s3, i, j));
        if (mdr1_safe (x1, i, j) > 0)
        {
          if (SIZE(s2))
            Mds1 (sw, i, j) = cpstr (mds1_safe (s2, i, j));
          else
            Mds1 (sw, i, j) = cpstr("UNDEF");
        }
        else if (mdr1_safe (x1, i, j) == 0)
        {
          if (SIZE(s3))
            Mds1 (sw, i, j) = cpstr (mds1_safe (s3, i, j));
          else
            Mds1 (sw, i, j) = cpstr("UNDEF");
        }
        else
        {
          if (SIZE(s4))
            Mds1 (sw, i, j) = cpstr (mds1_safe (s4, i, j));
          else
            Mds1 (sw, i, j) = cpstr("UNDEF");
        }
      }
    }
    ent_Clean (e1);
    ent_Clean (e2);
    ent_Clean (e3);
    ent_Clean (e4);
    return ent_Assign_Rlab_MDS(sw);
  }

  rerror("ifelse: 'iftrue' and 'iffalse' have to be the same type (string or number)!\n");
  return 0;
}


#include "getline_gnutime.c"

Ent *
ent_iseed (void)
{
  unsigned int lui;
  FILE *fh=0;

  fh = fopen ("/dev/random", "r");
  fread (&lui, sizeof (lui), 1, fh);
  fclose (fh);

  while (lui > 2082007225)
    lui >>= 1;

  return ent_Create_Rlab_Int (lui);
}

Ent *
ent_openpty_test (int nargs, Datum args[])
{
  int masterfd, slavefd;
  char *s = string_buff;

  extern int close(int);

  if (openpty (&masterfd, &slavefd, s, NULL, NULL) < 0)
  {
    perror  (NULL);
    fprintf (stderr, "Could not find the pty number!");
    s = 0;
  }

  close (masterfd);
  close (slavefd);

//   pid_t pid = vfork();
//   if(pid == 0)
//   {
//     // this is the child
//     xconsole (s);   // start process
//     exit(0);        // child exits
//   }

  return ent_Create_Rlab_String (s);
}


#include <sys/types.h>
#include <dirent.h>
#undef THIS_SOLVER
#define THIS_SOLVER "lsdir"

Ent *
ent_LsDir (int nargs, Datum args[])
{
  Ent *e1=0;

  char *s=0;
  MDS *rval=0;
  int nc, i;

  struct stat fileStat;
  DIR *dp = 0;
  struct dirent *ep=0;

  if (nargs != 1 && nargs != 0)
  {
    fprintf(stdout, THIS_SOLVER ": List content of a directory.\n");
    fprintf(stdout, THIS_SOLVER ": Format:\n");
    fprintf(stdout, THIS_SOLVER ":   t = " THIS_SOLVER "(/d/)\n");
    fprintf(stdout, THIS_SOLVER ": where 'd' is the directory name.\n");
    rerror ("one argument required");
  }

  if (nargs == 1)
  {
    e1 = bltin_get_ent (args[0]);
    if (ent_type (e1) != MATRIX_DENSE_STRING)
      rerror (THIS_SOLVER ": " RLAB_ERROR_ARG1_MDS_SCALAR);
    s = class_char_pointer (e1);
  }
  else
    s = "./";

  if(stat(s,&fileStat) >= 0)
  {
    if (S_ISDIR(fileStat.st_mode))
    {
      dp = opendir (s);
      nc = 0;
      if (dp != NULL)
      {
        while ((ep=readdir(dp)))
          nc++;
        closedir (dp);

        if (nc)
        {
          rval = mds_Create(nc,1);
          dp = opendir (s);
          i = 0;
          while ((ep=readdir (dp)) && (i<nc))
          {
            MdsV0(rval,i) = cpstr(ep->d_name);
            i++;
          }
          closedir (dp);
        }
      }
    }
  }

  ent_Clean(e1);
  return ent_Assign_Rlab_MDS(rval);
}

// *******************************************************************
//
// TERMCAP
//
// *******************************************************************
#include <curses.h>
#include <termcap.h>
static char *termcap_clr=0, *termcap_clre=0;

//
// clrscr(): clear screen
//
Ent *
ent_termcap_clrscr (int nargs, Datum args[])
{
  if (!termcap_clr)
    termcap_clr = (char *) tgetstr("cl", 0);

  printf("%s", termcap_clr);

  return ent_Create_Rlab_Success();
}

//
// mvcrsr(y,x); move cursor to line y and column x
//
Ent *
ent_termcap_mvcrsr (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0;
  int ix = 0, iy = 0;
  int my = tgetnum("li"), mx = tgetnum("co");

  if (nargs != 2)
    rerror ("mvcrsr: requires two arguments\n");

  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) != MATRIX_DENSE_REAL)
    rerror ("mvcrsr: improper first argument 'y'\n");
  iy = (int) class_double (e1);

  e2 = bltin_get_ent (args[1]);
  if (ent_type (e2) != MATRIX_DENSE_REAL)
    rerror ("mvcrsr: improper second argument 'x'\n");
  ix = (int) class_double (e2);

  ix = MIN(ix, mx);
  iy = MIN(iy, my);

  printf("\033[%d;%dH", iy, ix);
  fflush(stdout);

  ent_Clean (e1);
  ent_Clean (e2);

  return ent_Create_Rlab_Success();
}

//
// clrpos(y)  :  clear line y
// clrpos(y,x):  clear from row-y,column-x j number of characters
// y,x are either scalars or ranges
Ent *
ent_termcap_clrpos (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0;
  int i;
  int ix = 0, iy = 0, ny = 0, j, jd, nx=0;
  int my = tgetnum("li")+1, mx = tgetnum("co")+1;
  MDR *x, *y;

  if (nargs > 2)
    rerror ("clrpos: requires none, one or two arguments\n");

  // get row
  if (nargs >= 1)
  {
    e1 = bltin_get_ent (args[0]);
    if (ent_type (e1) != MATRIX_DENSE_REAL)
      rerror ("clrpos: improper first argument 'y'\n");
    y = ent_data (e1);
    ny = MNR(y) * MNC(y);
  }

  if (!termcap_clre)
    termcap_clre  = (char *) tgetstr("ce", 0); // clear from cursor to the end of line

  if (nargs == 2)
  {
    // get column
    e2 = bltin_get_ent (args[1]);
    if (ent_type (e1) != MATRIX_DENSE_REAL)
      rerror ("clrpos: improper second argument 'x'\n");
    x = ent_data (e2);
    nx = (x->nrow) * (x->ncol);
  }

  if (nargs == 0)
  {
    printf("%s", termcap_clre);
    fflush(stdout);
  }
  else if (nargs == 1)
  {
    for (j = 0; j < ny; j++)
    {
      // where to go
      if(y->type == RLAB_TYPE_INT32)
        jd = MIN(MdiV0(y,j),my);
      else
        jd = MIN((int) MdrV0(y,j),my);

      printf("\033[%d;%dH", jd, 1);
      printf("%s", termcap_clre);
      fflush(stdout);
    }
  }
  else
  {
    for (j = 0; j < ny ; j++)
    {
      for (i = 0; i < nx; i++)
      {
        // get the position
        // y:
        if(y->type == RLAB_TYPE_INT32)
          iy = MIN(MdiV0(y,j),my);
        else
          iy = MIN((int) MdrV0(y,j),my);
        // x:
        if(x->type == RLAB_TYPE_INT32)
          ix = MIN(MdiV0(x,i),mx);
        else
          ix = MIN((int) MdrV0(x,i),mx);

        // move there
        printf("\033[%d;%dH", iy, ix);
        printf(" ");
        fflush(stdout);
      }
    }
  }

  ent_Clean (e1);
  ent_Clean (e2);

  return ent_Create_Rlab_Success();
}

//
// colors(if,ib,bold):
//   set terminal to color mode given by
//   if - foreground color,
//   ib - background color
//   bold - bold
//
Ent * ent_termcap_colors (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0;
  int ic = 0, ib = 0;
  char *cmd = string_buff;
  char ESC[5]="\033";

  if (nargs !=0 && nargs != 1 && nargs !=2)
    rerror ("colors: requires two arguments\n");

  //
  // get foreground color
  //
  if (nargs >= 1)
  {
    e1 = bltin_get_ent (args[0]);
    if (ent_type (e1) == MATRIX_DENSE_REAL)
      ic = (int) class_double (e1);
    else if (ent_type (e1) == MATRIX_DENSE_STRING)
    {
      char *s1 = class_char_pointer(e1);
      if (!strcmp(s1, "default"))
        ic = 0;
      else if (!strcmp(s1, "red"))
        ic = 1;
      else if (!strcmp(s1, "green"))
        ic = 2;
      else if (!strcmp(s1, "yellow") || !strcmp(s1, "brown"))
        ic = 3;
      else if (!strcmp(s1, "blue"))
        ic = 4;
      else if (!strcmp(s1, "magenta"))
        ic = 5;
      else if (!strcmp(s1, "cyan"))
        ic = 6;
      else if (!strcmp(s1, "grey") || !strcmp(s1, "white"))
        ic = 7;
      else
        ic = 0;
    }
    else
      rerror ("colors: improper first argument 'color'\n");

    while (ic < 0)
      ic = ic + 8;
    while (ic > 7)
      ic = ic%8;
  }

  //
  // get background color
  //
  if (nargs == 2)
  {
    e2 = bltin_get_ent (args[1]);
    if (ent_type (e2) == MATRIX_DENSE_REAL)
      ib = (int) class_double (e2);
    else if (ent_type (e2) == MATRIX_DENSE_STRING)
    {
      char *s2 = class_char_pointer(e2);
      if (!strcmp(s2, "default"))
        ib = 0;
      else if (!strcmp(s2, "red"))
        ib = 1;
      else if (!strcmp(s2, "green"))
        ib = 2;
      else if (!strcmp(s2, "yellow") || !strcmp(s2, "brown"))
        ib = 3;
      else if (!strcmp(s2, "blue"))
        ib = 4;
      else if (!strcmp(s2, "magenta"))
        ib = 5;
      else if (!strcmp(s2, "cyan"))
        ib = 6;
      else if (!strcmp(s2, "grey") || !strcmp(s2, "white"))
        ib = 7;
      else
        ib = 0;
    }
    else
      rerror ("colors: improper second argument 'color'\n");

    while (ib < 0)
      ib = ib + 8;
    while (ib > 7)
      ib = ib%8;
  }

  if (nargs == 1)
    sprintf(cmd,"%s[0m%s[1;%2im", ESC, ESC, ic+30);
  else if (nargs == 2)
    sprintf(cmd,"%s[0m%s[%2im%s[1;%2im", ESC, ESC, ib+40, ESC, ic+30);
  else
    sprintf(cmd,"%s[0m", ESC);

  printf (cmd);
  fflush(stdout);

  ent_Clean (e1);
  ent_Clean (e2);

  return ent_Create_Rlab_Success();
}


/*
 * showkey.c -- display cooked key sequences
 *
Copyright (c) 2002, Eric S. Raymond

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:<P>

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.<P>

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.<P>

 *
 * Invoke this (no arguments needed) to see keycap-to-keystrokes mappings.
 *
 * by Eric S. Raymond <esr@snark.thyrsus.com>, 1 Nov 88
 * - fix for little-endian machines (version 1.1), 21 Oct 1996.
 * - cleanup and modern packaging (version 1.2), 1 Aug 2002.
 * - changed to use termios (version 1.3), 26 Aug 2002.
 * See the RPM spec file changelog for more recent stuff.
 */
#include <stdio.h>
#include <termios.h>
#include <signal.h>
#include <fcntl.h>   /* File control definitions */
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>

#define META  0x80
#define ESC 0x1b

Ent *
ent_showkey(int nargs, Datum args[])
{
  Ent *e1=0;
  struct termios  cooked, raw;
  unsigned char c;
  int i, meta;
  double dwait=0;

  MDR *ival=0;

  if (nargs >= 1)
  {
    e1 = bltin_get_ent (args[0]);
    if (ent_type (e1) == MATRIX_DENSE_REAL)
    {
      dwait = class_double (e1);
      if (dwait < 0)
        dwait = 0;
    }
  }

  // clear the buffer
  tcflush( STDIN_FILENO, TCIFLUSH );

  // Get the state of the tty
  tcgetattr(STDIN_FILENO, &cooked);

  // Make a copy we can mess with
  memcpy(&raw, &cooked, sizeof(struct termios));

  // Turn off echoing, linebuffering, and special-character processing,
  // but not the SIGINT or SIGQUIT keys.
  raw.c_lflag &= ~(ICANON|ECHO|ISIG);
  raw.c_iflag &= ~(INLCR|IGNCR|ICRNL|IUCLC);

  // Ship the raw control blts
  if (dwait)
  {
    raw.c_cc[VINTR] = 0377;
    raw.c_cc[VTIME] = (int) (10.0 * dwait);
    raw.c_cc[VMIN]  = 0;
  }
  tcsetattr(STDIN_FILENO, TCSANOW, &raw);

  // read one character from the input
  read(STDIN_FILENO, &c, 1);
  if ((meta = (c & META)))
    c &= ~META;

  //
  //  c==0: nothing was received. return empty array
  //  c!=0: check if there are more characters input as some keyboard strokes
  //        (e.g., cursor keys, function keys, and special keys)
  //        generate, so called, escape sequences.
  if (c)
  {
    int bytesWaiting;
    ioctl(STDIN_FILENO, FIONREAD, &bytesWaiting);
    ival = mdi_Create(1, bytesWaiting + 1);
    MdiV0(ival,0) = c;
    if (bytesWaiting)
    {
      for (i=1;i<=bytesWaiting; i++)
      {
        read(STDIN_FILENO, &c, 1);
        if ((meta = (c & META)))
          c &= ~META;
        MdiV0(ival,i) = c;
      }
    }
  }

  fflush(stdout);

  // Restore the cooked state
  tcsetattr(STDIN_FILENO, TCSANOW, &cooked);

  ent_Clean(e1);

  return ent_Assign_Rlab_MDR(ival);
}

//
// RLaB builtin function interface for getline().
//
#undef THIS_SOLVER
#define THIS_SOLVER "scanf"
Ent * Scanf (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0;
  char *s=0, *fmt=0;

  float rval=0;

  /* Check n_args */
  if (nargs != 2)
    rerror (RLAB_ERROR_TWO_ARG_REQUIRED);

  // get file name: handle user mistakes gracefully
  e1 = bltin_get_ent (args[0]);
  if (ent_type(e1)!=MATRIX_DENSE_STRING)
    rerror (THIS_SOLVER ": " RLAB_ERROR_ARG1_MDS_SCALAR);
  s = class_char_pointer (e1);
  if (isvalidstring(s)<1)
  {
    printf (THIS_SOLVER ": " RLAB_ERROR_ARG1_MDS_SCALAR "\n");
    goto _exit_scanf;
  }

  // get file name: handle user mistakes gracefully
  e2 = bltin_get_ent (args[1]);
  if (ent_type(e2)!=MATRIX_DENSE_STRING)
    rerror (THIS_SOLVER ": " RLAB_ERROR_ARG1_MDS_SCALAR);
  fmt = class_char_pointer (e2);
  if (isvalidstring(fmt)<1)
  {
    printf (THIS_SOLVER ": " RLAB_ERROR_ARG2_MDS_SCALAR "\n");
    goto _exit_scanf;
  }

  sscanf(s, fmt, &rval);

_exit_scanf:

  /* Now check for blank line */
  ent_Clean(e1);
  ent_Clean(e2);

  return ent_Create_Rlab_Double((double) rval);
}

//
// RLaB builtin function interface for getline().
//
#define       MAX_NUM_MATCH_VARID 32
struct _match_var
{
  char *name;
  unsigned long id;
  unsigned long offset;
  unsigned long len;
};
static struct _match_var
    match_var[MAX_NUM_MATCH_VARID];
static int
    num_match_variables = 0;

unsigned long GetVariableId(char * name, int len)
{
  static unsigned long rval;
  unsigned long i;

  // have we seen it before?
  for (i=0; i<num_match_variables; i++)
  {
    if (!strncmp (name, match_var[i].name,len))
    {
      rval = i+1;
//       fprintf(stderr, "GetVariableId:(found it) vid=%i, name = %s\n", rval, match_var[i].name);
      return rval;
    }
  }

  //  copy name
  match_var[num_match_variables].name = GC_malloc((len+1)*sizeof(char));
  for (i=0; i<len; i++)
    match_var[num_match_variables].name[i] = name[i];
  match_var[num_match_variables].name[len]='\0';

  // we don't have it
  num_match_variables++;

  //  return the latest addition
  rval = num_match_variables;

//   fprintf(stderr, "GetVariableId: vid=%i. nmv=%i, name = %s\n", rval, num_match_variables, match_var[num_match_variables-1].name);

  return rval;
}

void * AssignVariable(unsigned long VariableId, unsigned long Offset, unsigned long Length)
{
  match_var[VariableId].offset = Offset + 2;
  match_var[VariableId].len    = Length;
//   fprintf(stderr, "Assign: vid=%i, offset=%i, len=%i\n", VariableId, Offset, Length);
  return 0;
}

void * DeAssignVariable(unsigned long VariableId)
{
//   fprintf(stderr, "DeAssign: vid=%u\n", VariableId);
  match_var[VariableId].offset = 0;
  match_var[VariableId].len    = 0;
  return 0;
}


#include "sizeof.h"
#include "match.c"

#include "getline_match.c"

Ent * Match (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *e3=0;
  MDS *s=0;
  int i, ptr_off=0;
  MDR *rval=0;
  char *pattern=0;
  const char *processed_pattern;

  int istatus;
  int iptr=0;

  Btree *bw=0;

  /* Check n_args */
  if (nargs != 2 && nargs != 3)
    nerror (RLAB_ERROR_TWO_TO_THREE_ARG_REQUIRED);

  // clean-up the variables
  if (num_match_variables)
  {
    for (i=0; i<num_match_variables; i++)
    {
      GC_free(match_var[i].name);
      match_var[i].name=0;
      match_var[i].offset = 0;
      match_var[i].len    = 0;
    }
    num_match_variables=0;
  }

  //
  //  string matrix
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type(e1)!=MATRIX_DENSE_STRING)
    rerror (THIS_SOLVER ": " RLAB_ERROR_ARG1_MDS_SCALAR);
  s = ent_data(e1);
  if (EQNULL(s))
    rerror (THIS_SOLVER ": " RLAB_ERROR_ARG1_MDS_SCALAR);

  //
  // pattern for inspection of the string matrix
  //
  e2 = bltin_get_ent (args[1]);
  if (ent_type(e2)!=MATRIX_DENSE_STRING)
    rerror (THIS_SOLVER ": " RLAB_ERROR_ARG1_MDS_SCALAR);
  pattern = class_char_pointer (e2);
  if (!pattern)
    rerror (THIS_SOLVER ": " RLAB_ERROR_ARG1_MDS_SCALAR);
  if (!strlen(pattern))
    rerror (THIS_SOLVER ": " RLAB_ERROR_ARG1_MDS_SCALAR);

  //
  // offset
  //
  if (nargs > 2)
  {
    e3 = bltin_get_ent (args[2]);
    if (ent_type(e3)==MATRIX_DENSE_REAL)
    {
      ptr_off = class_int(e3) - 1;
      ptr_off = ptr_off > 0 ? ptr_off : 0;
    }
  }


  processed_pattern = patmaker(pattern);
  if ( MatchError == MATCH_SUCCESS)
  {
    rval = mdr_Create(MNR(s),MNC(s));
    for (i=0; i<SIZE(s); i++)
    {
      iptr = 0;
      if (ptr_off > strlen(MdsV0(s,i)))
      {
        MdrV0(rval,i) = 0;
        continue;
      }
      char *w = (char *) (MdsV0(s,i) + ptr_off);

      istatus = match (strlen(w), w,  &iptr, processed_pattern);

      if (istatus == MATCH_SUCCESS)
        MdrV0(rval,i) = (iptr + ptr_off);
      else
        MdrV0(rval,i) = 0;
    }
  }
  else
    fprintf(stderr, THIS_SOLVER ": Unable to process the pattern. Please read the manual while I stop for a moment!");

  // cleanup
  ent_Clean(e1);
  ent_Clean(e2);
  ent_Clean(e3);
  freedm();

  // clean-up the variables
  istatus=0;
  bw = btree_Create ();
  if (num_match_variables)
  {
    for (i=0; i<num_match_variables; i++)
    {
      if (match_var[i].name && match_var[i].offset>0 && match_var[i].len>0)
      {
        install (bw, "begin", ent_Create_Rlab_Double ((double) match_var[i].offset));
      }
      GC_free(match_var[i].name);
      match_var[i].name=0;
      match_var[i].offset = 0;
      match_var[i].len    = 0;
      istatus++;
      if (istatus > 1)
      {
        fprintf(stderr, THIS_SOLVER ": Using more then one assignment may result in unstable code. Refer to the manual while I stop for a moment!");
      }
    }
    num_match_variables=0;
  }

  install (bw, "end", ent_Assign_Rlab_MDR(rval));
  return ent_Assign_Rlab_BTREE(bw);
}

/* **************************************************************
 * The string-split function (strsplt)
 * ************************************************************** */

#undef THIS_SOLVER
#define THIS_SOLVER "strsplt"
//
// these are given further down
//
Ent * Strsplt (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0;
  MDS *mrv=0, *comment_pattern=0, *comment=0, *note=0, *note_pattern=0, *lstrip=0, *lstrip_pattern=0;
  char *str=0, *str_clean=0, *delim=0, c=0;
  FILE *rlab_stderr = (!RLAB_STDERR_DS) ? stderr : RLAB_STDERR_DS;
  int i, len, len_clean, j, min_line_len=1, ifoundcomment=-1;
  ListNode *node=0;

  // Check nargs
  if (nargs < 1 || nargs > 2)
  {
    fprintf (rlab_stderr, THIS_SOLVER ": " RLAB_ERROR_ONE_OR_TWO_ARG_REQUIRED "\n");
    goto _exit_strsplt;
  }

  //
  // get first argument: make sure it is a string of length greater than 1
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) != MATRIX_DENSE_STRING)
  {
    fprintf (rlab_stderr, THIS_SOLVER ": " RLAB_ERROR_ARG1_MDS_SCALAR "\n");
    goto _exit_strsplt;
  }
  if (SIZE (ent_data(e1)) != 1)
  {
    fprintf (rlab_stderr, THIS_SOLVER ": Warning: " RLAB_ERROR_ARG1_MDS_SCALAR "\n");
  }
  str = class_char_pointer(e1);
  len = isvalidstring(str);

  // failure: string not provided
  // or,
  // gracious exit: string of length 1 provided: nothing to split
  if (len<=1)
  {
    mrv = mds_CreateScalar(str);
    goto _exit_strsplt;
  }

  /* Split into single characters. */
  /* Create string matrix to hold split results */
  if (nargs == 1)
  {
    mrv = mds_Create (1, len);
    for (i=0; i<len; i++)
    {
      MdsV0(mrv,i) = (char *) GC_MALLOC (2 * sizeof (char));
      MdsV0(mrv,i)[0] = str[i];
      MdsV0(mrv,i)[1] = '\0';
    }
    goto _exit_strsplt;
  }

  //
  // is there a second argument
  //
  e2 = bltin_get_ent (args[1]);
  if (ent_type(e2) == MATRIX_DENSE_REAL)
  {
    // 2nd arg is a vector: assume fixed field widths
    mrv = string_split_mdr (str, len, (MDR*) ent_data(e2));
    goto _exit_strsplt;
  }
  if ( (ent_type(e2) != MATRIX_DENSE_STRING)
        && (ent_type(e2) != BTREE) )
  {
    fprintf(rlab_stderr, THIS_SOLVER ": " RLAB_ERROR_ARG2_STRING_OR_REAL_VECTOR);
    mrv = mds_CreateScalar(str);
    goto _exit_strsplt;
  }

  if (ent_type(e2) == MATRIX_DENSE_STRING)
  {
    delim = (char *) class_char_pointer (e2);
    if (isvalidstring(delim)<1)
    {
      fprintf (rlab_stderr, THIS_SOLVER ": " RLAB_ERROR_ARG2_MDS_SCALAR);
      goto _exit_strsplt;
    }
  }
  else if (ent_type(e2) == BTREE)
  {
    // second argument is a match-compatible list
    node = btree_FindNode (ent_data (e2), RLAB_NAME_READM_CSP);
    if (node)
    {
      delim = (char *) class_char_pointer (var_ent (node));
      if (isvalidstring(delim)<1)
      {
        fprintf (rlab_stderr, THIS_SOLVER ": " RLAB_ERROR_ARG2_MDS_SCALAR);
        goto _exit_strsplt;
      }
    }
    node = btree_FindNode (ent_data (e2), RLAB_NAME_READM_COMMENT);
    if (node)
    {
      if (ent_type(var_ent (node)) == MATRIX_DENSE_STRING)
      {
        comment = ent_data(var_ent (node));
        if (SIZE(comment)<1)
          comment=0;

          // process strings to create patterns
        if (comment)
        {
          comment_pattern = mds_Create(1,SIZE(comment));
          j=0;
          for (i=0; i<SIZE(comment);i++)
          {
            if (isvalidstring(MdsV0(comment,i))>0)
              MdsV0(comment_pattern,j++) = process_rlab_string_pattern(MdsV0(comment,i));
          }
          if (j)
            mds_Extend(comment_pattern,1,j);
          else
          {
            mds_Destroy(comment_pattern);
            comment_pattern = 0;
          }
        }
      }
    }
    node = btree_FindNode (ent_data (e2), RLAB_NAME_READM_NOTE);
    if (node)
    {
      if (ent_type(var_ent (node)) == MATRIX_DENSE_STRING)
      {
        note = ent_data(var_ent (node));
        if (SIZE(note)<1)
          note=0;

         // process strings to create patterns
        if (note)
        {
          note_pattern = mds_Create(1,SIZE(note));
          j=0;
          for (i=0; i<SIZE(note);i++)
          {
            if (isvalidstring(MdsV0(note,i))>0)
              MdsV0(note_pattern,j++) = process_rlab_string_pattern(MdsV0(note,i));
          }
          if (j)
            mds_Extend(note_pattern,1,j);
          else
          {
            mds_Destroy(note_pattern);
            note_pattern = 0;
          }
        }
      }
    }
    node = btree_FindNode (ent_data (e2), RLAB_NAME_READM_MIN_LINE_LEN);
    if (node)
    {
      min_line_len = class_int (var_ent (node));
      if (min_line_len < 1 )
        min_line_len = 1;
    }
    node = btree_FindNode (ent_data (e2), RLAB_NAME_READM_LSTRIP);
    if (node)
    {
      if (ent_type(var_ent (node)) == MATRIX_DENSE_STRING)
      {
        lstrip = ent_data(var_ent (node));
        if (SIZE(lstrip)<1)
          lstrip=0;

          // process strings to create patterns
        if (lstrip)
        {
          lstrip_pattern = mds_Create(1,SIZE(lstrip));
          j=0;
          for (i=0; i<SIZE(lstrip);i++)
          {
            if (isvalidstring(MdsV0(lstrip,i))>0)
              MdsV0(lstrip_pattern,j++) = process_rlab_string_pattern(MdsV0(lstrip,i));
          }
          if (j)
            mds_Extend(lstrip_pattern,1,j);
          else
          {
            mds_Destroy(lstrip_pattern);
            lstrip_pattern = 0;
          }
        }
      }
    }
  }

  // did we implement comments?
  //    if so, temporarily replace the first character of comment with '\n'
  //    to prevent match processing past it
  if (comment_pattern)
  {
    ifoundcomment = mds_find_first_left_pattern_in_string(str,comment_pattern,NULL);
    if (ifoundcomment >= 0)
    {
      // we have found comment
      if (ifoundcomment < min_line_len)
      {
        //  what is left after the comment is removed is not long enough
        //  to warrant further processing
        goto _exit_strsplt;
      }

      // we have modified string; update its length
      c = str[ifoundcomment];
      str[ifoundcomment] = '\0';
      len = isvalidstring(str);
    }
  }

  //
  // processing depends whether the notes are provided
  //
  if (note_pattern)
  {
    // notes are removed from string before the strings are processed
    str_clean = remove_pattern_from_string(note_pattern, str);
    len_clean = isvalidstring(str_clean);

    if (len_clean < min_line_len)
    {
      if (str_clean)
        GC_FREE(str_clean);

      goto _exit_strsplt;
    }

    // we implemented comments
    //   revert string to its previous value
    if (ifoundcomment > -1)
    {
      str[ifoundcomment] = c;
      ifoundcomment = -1;
    }

    str = str_clean;
    len = len_clean;
  }

  // finally: apply lstrip to remove anything superficial in the
  // string prior to processing
  int iblank = 0;
  if (lstrip_pattern)
  {
    iblank = processed_pattern_starts_string(str, lstrip_pattern);
  }

  // process string
  mrv = string_split_string (str+iblank, len-iblank, delim);

  if (str_clean)
    GC_FREE(str_clean);

  // we implemented comments
  //   revert string to its previous value
  if (ifoundcomment > -1)
  {
    str[ifoundcomment] = c;
    ifoundcomment = -1;
  }

_exit_strsplt:

  if (note_pattern)
    mds_Destroy(note_pattern);
  if (comment_pattern)
    mds_Destroy(comment_pattern);

  ent_Clean (e1);
  ent_Clean (e2);

  return ent_Assign_Rlab_MDS(mrv);
}

/* **************************************************************
 * Findstr(). This function mimics the Matlab findstr(). Matlab
 * does this with a m-file. But we have to use a builtin function
 * cause Rlab uses real string matrices.
 *
 * From the Matlab documentation
 * function k = findstr(s1,s2)
 * FINDSTR Find one string within another.
 * K = FINDSTR(S1,S2) returns the starting indices of any occurrences
 *  of the shorter of the two strings in the longer.
 *  Example:
 *      s = 'How much wood would a woodchuck chuck?';
 *      findstr(s,'a')    returns  21
 *      findstr(s,'wood') returns  [10 23]
 *      findstr(s,'Wood') returns  []
 *      findstr(s,' ')    returns  [4 9 14 20 22 32]
 * ************************************************************** */
#undef  THIS_SOLVER
#define THIS_SOLVER "findstr"
Ent * Findstr (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0;
  int len, dlen;
  MDR *re=0;
  char *str=0, *dlm=0;
  FILE *rlab_stderr = (!RLAB_STDERR_DS) ? stderr : RLAB_STDERR_DS;

  if (nargs != 2 && nargs != 3)
  {
    fprintf (rlab_stderr, THIS_SOLVER
        ": Finds first left location of pattern 'p' in entries of string matrix 's'.\n");
    fprintf (rlab_stderr, THIS_SOLVER
        ": Format:\n");
    fprintf (rlab_stderr, THIS_SOLVER
        ":   strindex(s, p),\n");
    fprintf (rlab_stderr, THIS_SOLVER
        ": where 's' and 'p' are string matrices.\n");
    rerror (RLAB_ERROR_TWO_ARG_REQUIRED);
  }

  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) != MATRIX_DENSE_STRING)
  {
    fprintf (rlab_stderr, THIS_SOLVER ": " RLAB_ERROR_ARG1_MDS_SCALAR);
    goto _findstr_exit;
  }
  str = class_char_pointer (e1);
  len = isvalidstring(str);

  e2 = bltin_get_ent (args[1]);
  if (ent_type (e2) != MATRIX_DENSE_STRING)
  {
    fprintf (rlab_stderr, THIS_SOLVER ": " RLAB_ERROR_ARG2_MDS_SCALAR);
    goto _findstr_exit;
  }
  dlm = class_char_pointer (e2);
  dlen = isvalidstring(dlm);

  if (len>0 && dlen>0)
  {
    re = mdr_Create(1,10);
    int i1=0, icol=0, ilast=0;
    char *w=str;
    int i0 = string_strindex_string(w+ilast, len, dlm, &i1);
    while (i0)
    {
      ilast += i0;
      MdrV0(re,icol++) = ilast;
      ilast += i1-1;
      if (icol==MNC(re))
        mdr_Extend(re,1,icol+5);
      i0 = string_strindex_string(w+ilast, strlen(w+ilast), dlm, &i1);
    }
    if (icol)
      mdr_Extend(re,1,icol);
    else
    {
      mdr_Destroy(re);
      re=0;
    }
  }

_findstr_exit:

  ent_Clean (e1);
  ent_Clean (e2);
  return ent_Assign_Rlab_MDR(re);
}


#undef  THIS_SOLVER
#define THIS_SOLVER "grep"
Ent * Grep (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0;
  int i, j, nr;
  MDS *re=0, *s=0, *s_pattern=0, *t=0;
  FILE *rlab_stderr = (!RLAB_STDERR_DS) ? stderr : RLAB_STDERR_DS;

  if (nargs != 2 && nargs != 3)
  {
    fprintf (rlab_stderr, THIS_SOLVER
        ": Return string entries from a string vector that contain desired pattern.\n");
    fprintf (rlab_stderr, THIS_SOLVER
        ": Format:\n");
    fprintf (rlab_stderr, THIS_SOLVER
        ":   "THIS_SOLVER"(s, p),\n");
    fprintf (rlab_stderr, THIS_SOLVER
        ": where 's' is a string vector, and 'p' the pattern.\n");
    rerror (RLAB_ERROR_TWO_ARG_REQUIRED);
  }

  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) != MATRIX_DENSE_STRING)
  {
    fprintf (rlab_stderr, THIS_SOLVER ": " RLAB_ERROR_ARG1_MDS_MATRIX);
    goto _strindex_exit;
  }
  t = ent_data (e1);
  nr = SIZE (t);
  if (nr<1)
  {
    fprintf (rlab_stderr, THIS_SOLVER ": " RLAB_ERROR_ARG1_MDS_MATRIX);
    goto _strindex_exit;
  }

  e2 = bltin_get_ent (args[1]);
  if (ent_type (e2) != MATRIX_DENSE_STRING)
  {
    fprintf (rlab_stderr, THIS_SOLVER ": " RLAB_ERROR_ARG2_MDS_MATRIX);
    goto _strindex_exit;
  }

  // process pattern in 's'
  s = ent_data (e2);
  s_pattern = mds_Create(1,SIZE(s));
  j=0;
  for (i=0; i<SIZE(s);i++)
  {
    if (isvalidstring(MdsV0(s,i))>0)
      MdsV0(s_pattern,j++) = process_rlab_string_pattern(MdsV0(s,i));
  }
  if (!j)
    goto _strindex_exit;

  int irow = 0;
  re = mds_Create (nr, 1);
  for (i=0; i<nr; i++)
  {
    char *str = MdsV0(t, i);
    if (mds_find_first_left_pattern_in_string(str, s_pattern, NULL)>-1)
    {
      MdsV0(re,irow) = cpstr(str);
      irow++;
    }
  }
  if (irow)
    mds_Extend(re, irow, 1);
  else
  {
    mds_Destroy(re);
    re = 0;
  }

_strindex_exit:

  if (s_pattern)
    mds_Destroy(s_pattern);

  ent_Clean (e1);
  ent_Clean (e2);
  return ent_Assign_Rlab_MDS(re);
}

#undef  THIS_SOLVER
#define THIS_SOLVER "strindex"
Ent * ent_string_index (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0;
  int i, j, nr, nc, nr_s, nc_s, nr_c, nc_c;
  MDR *re=0;
  MDS *s=0, *c=0;
  FILE *rlab_stderr = (!RLAB_STDERR_DS) ? stderr : RLAB_STDERR_DS;

  if (nargs != 2 && nargs != 3)
  {
    fprintf (rlab_stderr, THIS_SOLVER
        ": Finds first left location of pattern 'p' in entries of string matrix 's'.\n");
    fprintf (rlab_stderr, THIS_SOLVER
        ": Format:\n");
    fprintf (rlab_stderr, THIS_SOLVER
        ":   strindex(s, p),\n");
    fprintf (rlab_stderr, THIS_SOLVER
        ": where 's' and 'p' are string matrices.\n");
    rerror (RLAB_ERROR_TWO_ARG_REQUIRED);
  }

  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) != MATRIX_DENSE_STRING)
  {
    re = mdr_CreateScalar(-1.0);
    goto _strindex_exit;
  }

  s = ent_data (e1);
  if (SIZE(s)<1)
  {
    re = mdr_CreateScalar (0);
    goto _strindex_exit;
  }

  nr_s = MNR (s);
  nc_s = MNC (s);

  e2 = bltin_get_ent (args[1]);
  if (ent_type (e2) != MATRIX_DENSE_STRING)
  {
    re = mdr_CreateScalar(-1.0);
    goto _strindex_exit;
  }
  c = ent_data (e2);
  nr_c = MNR (c);
  nc_c = MNC (c);

  nr = MAX(nr_s, nr_c);
  nc = MAX(nc_s, nc_c);

  re = mdr_Create (nr, nc);
  for (i=0; i<nr; i++)
  {
    for (j=0; j<nc; j++)
    {
      char *str = Mds0(s, MIN(i,nr_s-1), MIN(j,nc_s-1));
      int   len = isvalidstring(str);
      char *dlm = Mds0(c, MIN(i,nr_c-1), MIN(j,nc_c-1));
      int  dlen = isvalidstring(dlm);
      if (len>0 && dlen>0)
        Mdr0 (re,i,j) = string_strindex_string(str, len, dlm, NULL);
      else
        Mdr0 (re,i,j) = 0.0;
    }
  }

_strindex_exit:

  ent_Clean (e1);
  ent_Clean (e2);
  return ent_Assign_Rlab_MDR(re);
}

#undef THIS_SOLVER
#define THIS_SOLVER "strindexr"
Ent * ent_string_index_pattern (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0;
  int i,i0,i1,len;
  MDR *re=0;
  char *str=0, *delim=0;
  FILE *rlab_stderr = (!RLAB_STDERR_DS) ? stderr : RLAB_STDERR_DS;

  if (nargs != 2 && nargs != 3)
  {
    fprintf (rlab_stderr, THIS_SOLVER
        ": Finds indices of the first left location of pattern 'p' in entries of string matrix 's'.\n");
    fprintf (rlab_stderr, THIS_SOLVER
        ": Format:\n");
    fprintf (rlab_stderr, THIS_SOLVER
        ":   "THIS_SOLVER"(s, p),\n");
    fprintf (rlab_stderr, THIS_SOLVER
        ": where 's' and 'p' are string scalars.\n");
    rerror (RLAB_ERROR_TWO_ARG_REQUIRED);
  }

  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) != MATRIX_DENSE_STRING)
  {
    fprintf (rlab_stderr, THIS_SOLVER ": " RLAB_ERROR_ARG1_MDS_MATRIX);
    goto _strindex_exit;
  }

  str = class_char_pointer (e1);
  len = isvalidstring(str);
  if (len<1)
  {
    fprintf (rlab_stderr, THIS_SOLVER ": " RLAB_ERROR_ARG1_MDS_SCALAR);
    goto _strindex_exit;
  }

  e2 = bltin_get_ent (args[1]);
  if (ent_type (e2) != MATRIX_DENSE_STRING)
  {
    fprintf (rlab_stderr, THIS_SOLVER ": " RLAB_ERROR_ARG2_MDS_MATRIX);
    goto _strindex_exit;
  }
  delim = class_char_pointer (e2);
  if (isvalidstring(delim)<1)
  {
    fprintf (rlab_stderr, THIS_SOLVER ": " RLAB_ERROR_ARG2_MDS_SCALAR);
    goto _strindex_exit;
  }

  i0 = string_strindex_string(str, len, delim, &i1);
  if (i0)
  {
    re = mdr_Create (1, i1);
    for (i=0;i<i1;i++)
      MdrV0(re,i) = i0 + i;
  }

_strindex_exit:

  ent_Clean (e1);
  ent_Clean (e2);
  return ent_Assign_Rlab_MDR(re);
}



// **************************************************************
//
// Strtod (string-to-decimal) function...
//
// **************************************************************
#undef THIS_SOLVER
#define THIS_SOLVER "strtod"
Ent * Strtod (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0;
  MDS *s=0, *comment=0, *comment_pattern=0, *note=0, *note_pattern=0, *lstrip=0, *lstrip_pattern=0;
  MDS *start=0, *start_pattern=0, *stop=0, *stop_pattern=0;
  MDR *rval=0, *userows=0, *set_userows=0, *usecols=0, *set_usecols=0;
  char *delim=0, *end, *join_csp=0;
  int i,j,nr,nc,min_line_len=1,iskiprows=0, join_rows=1, join_rows_counter=0;
  int start_pattern_notfound=1, stop_pattern_notfound=1, max_icols=0;
  FILE *rlab_stderr = (!RLAB_STDERR_DS) ? stderr : RLAB_STDERR_DS;
  ListNode *node=0;

  if (nargs != 1 && nargs != 2)
  {
    fprintf (rlab_stderr, THIS_SOLVER
        ": Converts string scalar or vector to a real number matrix using a delimiter.\n");
    fprintf (rlab_stderr, THIS_SOLVER
        ": Format:\n");
    fprintf (rlab_stderr, THIS_SOLVER
        ":   "THIS_SOLVER"(s/, d/),\n");
    fprintf (rlab_stderr, THIS_SOLVER
        ": where 's' is a string matrix or vector of data, while 'd' is an optional string delimiter.\n");
    rerror (RLAB_ERROR_TWO_ARG_REQUIRED);
  }

  //
  // first argument: string matrix
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) == MATRIX_DENSE_REAL)
  {
    rval = mdr_Copy(ent_data(e1));
    goto _exit_strtod;
  }
  if (ent_type (e1) != MATRIX_DENSE_STRING)
  {
    fprintf (rlab_stderr, THIS_SOLVER ": " RLAB_ERROR_ARG1_MDS_MATRIX);
    goto _exit_strtod;
  }
  s = ent_data (e1);
  nr = MNR(s);
  nc = MNC(s);

  //
  // second argument: delimiter (string scalar), or list
  //
  if (nargs == 2)
  {
    e2 = bltin_get_ent (args[1]);
    if (ent_type (e2) == MATRIX_DENSE_STRING)
    {
      delim = class_char_pointer (e2);
      if (isvalidstring(delim)<1)
      {
        fprintf (rlab_stderr, THIS_SOLVER ": " RLAB_ERROR_ARG2_MDS_SCALAR);
        goto _exit_strtod;
      }
    }
    else if (ent_type (e2) == BTREE)
    {
      RLABCODE_PROCESS_BTREE_ENTRY_S(e2,node,RLAB_NAME_READM_CSP,delim,1,NULL);

      RLABCODE_PROCESS_BTREE_ENTRY_D(e2,node,RLAB_NAME_READM_MIN_LINE_LEN,min_line_len,class_double,>,0);

      RLABCODE_PROCESS_BTREE_ENTRY_D(e2,node,RLAB_NAME_READM_JOINROWS,join_rows,class_double,>,0);

      RLABCODE_PROCESS_BTREE_ENTRY_S(e2,node,RLAB_NAME_READM_JOINCSP,join_csp,1,NULL);

      node = btree_FindNode (ent_data (e2), RLAB_NAME_READM_COMMENT);
      if (node)
      {
        if (ent_type(var_ent (node)) == MATRIX_DENSE_STRING)
        {
          comment = ent_data(var_ent (node));
          if (SIZE(comment)<1)
            comment=0;

          // process strings to create patterns
          if (comment)
          {
            comment_pattern = mds_Create(1,SIZE(comment));
            j=0;
            for (i=0; i<SIZE(comment);i++)
            {
              if (isvalidstring(MdsV0(comment,i))>0)
                MdsV0(comment_pattern,j++) = process_rlab_string_pattern(MdsV0(comment,i));
            }
            if (j)
              mds_Extend(comment_pattern,1,j);
            else
            {
              mds_Destroy(comment_pattern);
              comment_pattern = 0;
            }
          }
        }
      }
      node = btree_FindNode (ent_data (e2), RLAB_NAME_READM_NOTE);
      if (node)
      {
        if (ent_type(var_ent (node)) == MATRIX_DENSE_STRING)
        {
          note = ent_data(var_ent (node));
          if (SIZE(note)<1)
            note=0;

          // process strings to create patterns
          if (note)
          {
            note_pattern = mds_Create(1,SIZE(note));
            j=0;
            for (i=0; i<SIZE(note);i++)
            {
              if (isvalidstring(MdsV0(note,i))>0)
                MdsV0(note_pattern,j++) = process_rlab_string_pattern(MdsV0(note,i));
            }
            if (j)
              mds_Extend(note_pattern,1,j);
            else
            {
              mds_Destroy(note_pattern);
              note_pattern = 0;
            }
          }
        }
      }
      node = btree_FindNode (ent_data (e2), RLAB_NAME_READM_LSTRIP);
      if (node)
      {
        if (ent_type(var_ent (node)) == MATRIX_DENSE_STRING)
        {
          lstrip = ent_data(var_ent (node));
          if (SIZE(lstrip)<1)
            lstrip=0;

          // process strings to create patterns
          if (lstrip)
          {
            lstrip_pattern = mds_Create(1,SIZE(lstrip));
            j=0;
            for (i=0; i<SIZE(lstrip);i++)
            {
              if (isvalidstring(MdsV0(lstrip,i))>0)
                MdsV0(lstrip_pattern,j++) = process_rlab_string_pattern(MdsV0(lstrip,i));
            }
            if (j)
              mds_Extend(lstrip_pattern,1,j);
            else
            {
              mds_Destroy(lstrip_pattern);
              lstrip_pattern = 0;
            }
          }
        }
      }
      node = btree_FindNode (ent_data (e2), RLAB_NAME_READM_START);
      if (node)
      {
        if (ent_type(var_ent (node)) == MATRIX_DENSE_STRING)
        {
          start = ent_data(var_ent (node));
          if (SIZE(start)<1)
            start=0;

          // process strings to create patterns
          if (start)
          {
            start_pattern = mds_Create(1,SIZE(start));
            j=0;
            for (i=0; i<SIZE(start);i++)
            {
              if (isvalidstring(MdsV0(start,i))>0)
                MdsV0(start_pattern,j++) = process_rlab_string_pattern(MdsV0(start,i));
            }
            if (j)
              mds_Extend(start_pattern,1,j);
            else
            {
              mds_Destroy(start_pattern);
              start_pattern = 0;
            }
          }
        }
      }
      node = btree_FindNode (ent_data (e2), RLAB_NAME_READM_STOP);
      if (node)
      {
        if (ent_type(var_ent (node)) == MATRIX_DENSE_STRING)
        {
          stop = ent_data(var_ent (node));
          if (SIZE(stop)<1)
            stop=0;

          // process strings to create patterns
          if (stop)
          {
            stop_pattern = mds_Create(1,SIZE(stop));
            j=0;
            for (i=0; i<SIZE(stop);i++)
            {
              if (isvalidstring(MdsV0(stop,i))>0)
                MdsV0(stop_pattern,j++) = process_rlab_string_pattern(MdsV0(stop,i));
            }
            if (j)
              mds_Extend(stop_pattern,1,j);
            else
            {
              mds_Destroy(stop_pattern);
              stop_pattern = 0;
            }
          }
        }
      }
      node = btree_FindNode (ent_data (e2), RLAB_NAME_READM_SKIPROWS);
      if (node)
      {
        iskiprows = class_int (var_ent (node));
        if (iskiprows <= 0 )
          iskiprows = 0;
      }
      node = btree_FindNode (ent_data (e2), RLAB_NAME_READM_USEROWS);
      if (node)
      {
        if (ent_type(var_ent (node)) == MATRIX_DENSE_REAL)
        {
          userows = ent_data(var_ent (node));
          if (SIZE(userows)>0)
            set_userows = mdr_VectorSet(userows);
          if (SIZE(set_userows)<1)
          {
            mdr_Destroy(set_userows);
            set_userows=0;
            userows=0;
          }
        }
      }
      node = btree_FindNode (ent_data (e2), RLAB_NAME_READM_USECOLS);
      if (node)
      {
        if (ent_type(var_ent (node)) == MATRIX_DENSE_REAL)
        {
          usecols = ent_data(var_ent (node));
          if (SIZE(usecols)>0)
            set_usecols = mdr_VectorSet(usecols);
          if (SIZE(set_usecols)<1)
          {
            mdr_Destroy(set_usecols);
            set_usecols=0;
            usecols=0;
          }
        }
      }
    }
    else
    {
      fprintf (rlab_stderr, THIS_SOLVER ": " RLAB_ERROR_ARG2_MDS_MATRIX);
      goto _exit_strtod;
    }
  }

  if (join_rows>1 && !join_csp)
    join_csp = " ";

  if ( EQCOLVECT(s) )
  {
    //
    // column vector
    //
    if (isvalidstring(delim)<1)
      delim = "'(BLANK|',')";

    int nr_max = nr;
    int idx_set_userows=0;
    if (set_userows)
      nr_max = MIN( nr_max, mdiV1(set_userows,SIZE(set_userows)) );

    if (set_usecols)
    {
      rval = mdr_Create(nr_max,SIZE(set_usecols));
    }
    else
      rval = mdr_Create(nr_max,20);

    mdr_Nan(rval);
    int irow=0,icol=0, ifoundcomment=-1;
    char c=0, *joined_str=0;

    for (j=0; j<nr; j++)
    {
      // do we skip it?
      if (j<iskiprows)
        continue;

      // did user provide list of rows they wants processed?
      if (set_userows)
      {

        if (j==0)
        {
          while ((j>(mdiV0(set_userows,idx_set_userows)-1)) && (idx_set_userows < SIZE(set_userows)) )
            idx_set_userows++;
        }

        if (idx_set_userows == SIZE(set_userows))
          break;

        if (j!=(mdiV0(set_userows,idx_set_userows)-1))
          continue;

        idx_set_userows++;
      }

      char *str=MdsV0(s,j), *str_clean=0;
      int len_clean,len=isvalidstring(str);

      // do we check for minimal length of the string
      if (len < min_line_len)
        continue;

      // if user provided start pattern: then we look for it before we can
      // process the lines
      if (start_pattern && start_pattern_notfound==1)
      {
        if (processed_pattern_starts_string(str, start_pattern) < 1)
          continue;

        start_pattern_notfound=0;
        if (stop_pattern)
          stop_pattern_notfound=1;
        continue;
      }

      // if user provided stop pattern: then we stop processing the lines
      if (stop_pattern && stop_pattern_notfound==1)
      {
        if (processed_pattern_starts_string(str, stop_pattern) > 0)
        {
          if (!start_pattern)
            break;

          stop_pattern_notfound=0;
          start_pattern_notfound=1;
        }
      }

      if (!stop_pattern_notfound)
        continue;

      // did we implement comments?
      //    if so, temporarily replace the first character of comment with '\n'
      //    to prevent match processing past it
      if (comment_pattern)
      {
        ifoundcomment = mds_find_first_left_pattern_in_string(str,comment_pattern,NULL);
        if (ifoundcomment >= 0)
        {
          // we have found comment
          if (ifoundcomment < min_line_len)
          {
            //  what is left after the comment is removed is not long enough
            //  to warrant further processing
            continue;
          }

          // we have modified string; update its length
          c = str[ifoundcomment];
          str[ifoundcomment] = '\0';
          len = isvalidstring(str);
        }
      }

      // at this point we have a row
      if (join_rows > 1)
      {
        join_rows_counter++;
        if (join_rows_counter > 1)
          string_concat(&joined_str, join_csp);
        string_concat(&joined_str, str);
        if (join_rows_counter < join_rows)
          continue;

        // fix 'str' before replacing it with another pointer
        if (ifoundcomment > -1)
        {
          str[ifoundcomment] = c;
          ifoundcomment = -1;
        }
        str = joined_str;
        join_rows_counter=0;
      }

      //
      // processing depends whether the notes are provided
      //
      if (note_pattern)
      {
        // notes are removed from string before the strings are processed
        str_clean = remove_pattern_from_string(note_pattern, str);
        len_clean = isvalidstring(str_clean);
        if (len_clean < min_line_len)
        {
          GC_FREE (str_clean);
          continue;
        }

        // fix 'str' before replacing it with another pointer
        if (ifoundcomment > -1)
        {
          str[ifoundcomment] = c;
          ifoundcomment = -1;
        }

        str = str_clean;
        len = len_clean;
      }

      // finally: apply lstrip to remove anything superficial in the
      // string prior to processing
      int iblank = 0;
      if (lstrip_pattern)
      {
        iblank = processed_pattern_starts_string(str, lstrip_pattern);
      }

      // process string
      strtod_split_string (&rval, &irow, &icol, str+iblank, len-iblank, delim, set_usecols, &max_icols);

      // we implemented comments
      //   revert string to its previous value but only if notes and join_lines have not affected
      //   its position
      if (ifoundcomment > -1)
      {
        str[ifoundcomment] = c;
        ifoundcomment = -1;
      }

      if (str_clean)
      {
        GC_FREE (str_clean);
        str_clean = 0;
      }
      if (joined_str)
      {
        GC_FREE(joined_str);
        joined_str = 0;
      }

      // we check for minimal length of the string
      if (min_line_len)
      {
        // we ignore empty lines
        if (!icol)
          continue;
      }
      irow++;

      if (!rval)
        break;
    }

    if ((SIZE(rval)<1) || irow==0)
    {
      mdr_Destroy(rval);
      rval = 0;
    }
    else
      mdr_Extend(rval, irow, MNC(rval));
  }
  else
  {
    // string matrix is parsed field by field, where to each field
    //  strtod is applied without delimiters
    rval = mdr_Create(nr,nc);
    mdr_Nan(rval);

    for (i=0; i<nr; i++)
    {
      for (j=0; j<nc; j++)
      {
        char *str=Mds0(s,i,j), *str_clean=0;
        int   len_clean;

        //
        // processing depends whether the notes are provided
        //
        if (note_pattern)
        {
          // notes are removed from string before the strings are processed
          str_clean = remove_pattern_from_string(note_pattern, str);
          len_clean = isvalidstring(str_clean);
          if (len_clean < min_line_len)
          {
            GC_FREE (str_clean);
            str_clean = 0;
            continue;
          }

          str = str_clean;
        }

        // finally: apply lstrip to remove anything superficial in the
        // string prior to processing
        int iblank = 0;
        if (lstrip_pattern)
        {
          iblank = processed_pattern_starts_string(str, lstrip_pattern);
        }

        Mdr0(rval,i,j) = strtod(str + iblank, &end);
        if ( (end == str+iblank) && (Mdr0(rval,i,j)==0) )
        {
          Mdr0(rval,i,j) = create_nan();
        }

        if (str_clean)
          GC_FREE (str_clean);
      }
    }
  }

_exit_strtod:

  if (comment_pattern)
    mds_Destroy(comment_pattern);
  if (note_pattern)
    mds_Destroy(note_pattern);
  if (lstrip_pattern)
    mds_Destroy(lstrip_pattern);
  if(stop_pattern)
    mds_Destroy(stop_pattern);
  if(start_pattern)
    mds_Destroy(start_pattern);
  if (set_userows)
    mdr_Destroy(set_userows);
  if (set_usecols)
    mdr_Destroy(set_usecols);

  if (!rval)
  {
    rval = mdr_Create(1,1);
    MdrV0(rval,0) = create_nan();
  }

  ent_Clean (e1);
  ent_Clean (e2);
  return ent_Assign_Rlab_MDR(rval);
}


#undef  THIS_SOLVER
#define THIS_SOLVER "gsub"
//
// general substitution
//  of processed source 's' with 'r' in target string 't'
//
Ent *
ent_string_substitute (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *e3=0;
  int i, j, k, nr, nc;
  MDS  *t=0, *re=0, *r=0, *s_pattern=0;
  char *s=0;
  MDR *rk=0;

  if ((nargs != 2) && (nargs != 3))
  {
    fprintf (stdout,
             THIS_SOLVER ": General sub-string substition with counting.\n");
    fprintf (stdout,
             THIS_SOLVER ": Format:\n");
    fprintf (stdout,
             THIS_SOLVER ":   gsub(r,s,t),\n");
    fprintf (stdout,
             THIS_SOLVER ":   gsub(s,t),\n");
    fprintf (stdout,
             THIS_SOLVER ": where 't' is the original string matrix, in each entry of which\n");
    fprintf (stdout,
             THIS_SOLVER ": string 's' is replaced by a string 'r', or removed if 'r' is absent.\n");
    rerror (RLAB_ERROR_THREE_ARG_REQUIRED);
  }

  if (nargs == 3)
  {
    // get replacement string
    e1 = bltin_get_ent (args[0]);
    if (ent_type (e1) != MATRIX_DENSE_STRING)
    {
      rerror(THIS_SOLVER ": " RLAB_ERROR_ARG1_MDS_MATRIX "\n");
    }
    r = ent_data (e1);
    if (SIZE(r)<1)
      r=0;

    // get search string
    e2 = bltin_get_ent (args[1]);
    if (ent_type (e2) != MATRIX_DENSE_STRING)
    {
      rerror(THIS_SOLVER ": " RLAB_ERROR_ARG2_MDS_SCALAR "\n");
    }
    s = ent_data (e2);

    // get string matrix upon which substitution is going to be performed
    e3 = bltin_get_ent (args[2]);
    if (ent_type (e3) != MATRIX_DENSE_STRING)
    {
      rerror(THIS_SOLVER ": " RLAB_ERROR_ARG3_MDS_MATRIX "\n");
    }
    t = ent_data (e3);
  }
  else if (nargs ==2)
  {
    // get search string
    e1 = bltin_get_ent (args[0]);
    if (ent_type (e1) != MATRIX_DENSE_STRING)
    {
      rerror(THIS_SOLVER ": " RLAB_ERROR_ARG1_MDS_SCALAR "\n");
    }
    s = ent_data (e1);

    // get string matrix upon which substitution is going to be performed
    e2 = bltin_get_ent (args[1]);
    if (ent_type (e2) != MATRIX_DENSE_STRING)
    {
      rerror(THIS_SOLVER ": " RLAB_ERROR_ARG2_MDS_MATRIX "\n");
    }
    t = ent_data (e2);
  }

  // nothing to do:
  if (SIZE(t)<1 || SIZE(s)<1)
  {
    re = mds_Copy(t);
    rk = mdr_CreateScalar(0);
    goto _exit_gsub;
  }

  // process pattern in 's'
  s_pattern = mds_Create(1,SIZE(s));
  j=0;
  for (i=0; i<SIZE(s);i++)
  {
    if (isvalidstring(MdsV0(s,i))>0)
      MdsV0(s_pattern,j++) = process_rlab_string_pattern(MdsV0(s,i));
  }
  if (!j)
  {
    re = mds_Copy(t);
    rk = mdr_CreateScalar(0);
    goto _exit_gsub;
  }
  mds_Extend(s_pattern,1,j);

  nr = MNR (t);
  nc = MNC (t);
  rk = mdi_Create (nr, nc);
  mdr_Zero (rk);

  re = mds_Create (nr, nc);

  for (i=0; i<nr*nc; i++)
  {
    MdsV0 (re, i) = cpstr(MdsV0(t, i));
    for (j=0; j<SIZE(s_pattern); j++)
    {
      char * str = MdsV0 (re, i);
      if (isvalidstring(str)<1)
        continue;

      char * dlm = MdsV0(s_pattern,j);
      char * rep = 0;
      if (r)
        rep = MdsV0(r,MIN(j, SIZE(r)-1));
      MdsV0 (re, i)  = replace_processed_pattern_in_string(rep, dlm, str, &k);
      GC_free(str);
      MdiV0 (rk, i) += k;
    }
  }

_exit_gsub:

  ent_Clean (e1);
  ent_Clean (e2);
  ent_Clean (e3);

  if (s_pattern)
    mds_Destroy(s_pattern);

  //
  // output
  //
  Btree *bw = btree_Create ();
  install  (bw, "string", ent_Assign_Rlab_MDS(re)); // new string matrix
  install  (bw, "count", ent_Assign_Rlab_MDR(rk));  //replacement count
  return ent_Assign_Rlab_BTREE(bw);
}


#undef THIS_SOLVER
#define THIS_SOLVER "lstrip"
Ent * Lstrip (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0;
  MDS *t=0, *r=0, *s=0, *s_pattern=0;
  int i, j;

  // Check number of arguments

  if (nargs != 2)
  {
    printf (THIS_SOLVER ": " RLAB_ERROR_TWO_ARG_REQUIRED "\n");
    goto _exit_lstrip;
  }

  /* Get the first argument */
  e1 = bltin_get_ent (args[0]);
  if (ent_type(e1) != MATRIX_DENSE_STRING)
    goto _exit_lstrip;
  t  = ent_data (e1);
  if (SIZE(t)<1)
    goto _exit_lstrip;

  // get pattern in 's'
  e2 = bltin_get_ent (args[1]);
  if (ent_type (e2) != MATRIX_DENSE_STRING)
  {
    r = mds_Copy(t);
    goto _exit_lstrip;
  }
  s = ent_data (e2);
  if (SIZE(s)<1)
  {
    r = mds_Copy(t);
    goto _exit_lstrip;
  }

  s_pattern = mds_Create(1,SIZE(s));
  j=0;
  for (i=0; i<SIZE(s);i++)
  {
    if (isvalidstring(MdsV0(s,i))>0)
    {
      MdsV0(s_pattern,j) = process_rlab_string_pattern(MdsV0(s,i));
      j++;
    }
  }
  if (!j)
  {
    r = mds_Copy(t);
    goto _exit_lstrip;
  }
  if (j != SIZE(s))
    mds_Extend(s_pattern,1,j);

  r = mds_Create(MNR(t), MNC(t));

  for (i=0; i<SIZE(t); i++)
  {
    char * str = MdsV0(t, i);
    int iend = processed_pattern_starts_string(str, s_pattern);
    MdsV0 (r, i) = cpstr(str+iend);
  }

_exit_lstrip:

  if (s_pattern)
    mds_Destroy(s_pattern);

  ent_Clean (e1);
  ent_Clean (e2);
  return ent_Assign_Rlab_MDS(r);
}

#undef THIS_SOLVER
#define THIS_SOLVER "rstrip"
Ent * Rstrip (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0;
  MDS *t=0, *r=0, *s=0, *s_pattern=0;
  int i, j;

  // Check number of arguments

  if (nargs != 2)
  {
    printf (THIS_SOLVER ": " RLAB_ERROR_TWO_ARG_REQUIRED "\n");
    goto _exit_rstrip;
  }

  /* Get the first argument */
  e1 = bltin_get_ent (args[0]);
  if (ent_type(e1) != MATRIX_DENSE_STRING)
    goto _exit_rstrip;
  t  = ent_data (e1);
  if (SIZE(t)<1)
    goto _exit_rstrip;

  // get pattern in 's'
  e2 = bltin_get_ent (args[1]);
  if (ent_type (e2) != MATRIX_DENSE_STRING)
  {
    r = mds_Copy(t);
    goto _exit_rstrip;
  }
  s = ent_data (e2);
  if (SIZE(s)<1)
  {
    r = mds_Copy(t);
    goto _exit_rstrip;
  }

  s_pattern = mds_Create(1,SIZE(s));
  j=0;
  for (i=0; i<SIZE(s);i++)
  {
    if (isvalidstring(MdsV0(s,i))>0)
      MdsV0(s_pattern,j++) = process_rlab_string_pattern(MdsV0(s,i));
  }
  if (!j)
  {
    r = mds_Copy(t);
    goto _exit_rstrip;
  }
  if (j != SIZE(s))
    mds_Extend(s_pattern,1,j);

  r = mds_Create(MNR(t), MNC(t));

  for (i=0; i<SIZE(t); i++)
  {
    char * str = MdsV0(t, i);
    int iend = does_processed_pattern_end_string(str, s_pattern);
    MdsV0 (r, i) = cpnstr(str,iend);
  }

_exit_rstrip:

  if (s_pattern)
    mds_Destroy(s_pattern);

  ent_Clean (e1);
  ent_Clean (e2);
  return ent_Assign_Rlab_MDS(r);
}

#undef THIS_SOLVER
#define THIS_SOLVER "strip"
Ent * Strip (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *e3=0;
  MDS *t=0, *r=0, *l=0, *l_pattern=0,  *r_pattern=0, *s=0;
  int i, j;

  // Check number of arguments

  if (nargs != 3)
  {
    printf (THIS_SOLVER ": " RLAB_ERROR_THREE_ARG_REQUIRED "\n");
    goto _exit_strip;
  }

  /* Get the first argument */
  e1 = bltin_get_ent (args[0]);
  if (ent_type(e1) != MATRIX_DENSE_STRING)
    goto _exit_strip;
  t  = ent_data (e1);
  if (SIZE(t)<1)
    goto _exit_strip;

  // get pattern in 's'
  e2 = bltin_get_ent (args[1]);
  if (ent_type (e2) == MATRIX_DENSE_STRING)
    l = ent_data (e2);
  if (SIZE(l)>0)
  {
    l_pattern = mds_Create(1,SIZE(l));
    j=0;
    for (i=0; i<SIZE(l);i++)
    {
      if (isvalidstring(MdsV0(l,i))>0)
        MdsV0(l_pattern,j++) = process_rlab_string_pattern(MdsV0(l,i));
    }
    if (!j)
    {
      mds_Destroy(l_pattern);
      l_pattern = 0;
    }
    else if (j != SIZE(l))
      mds_Extend(l_pattern,1,j);
  }

  // get pattern in 's'
  e3 = bltin_get_ent (args[2]);
  if (ent_type (e3) == MATRIX_DENSE_STRING)
    r = ent_data (e3);
  if (SIZE(r)>0)
  {
    r_pattern = mds_Create(1,SIZE(r));
    j=0;
    for (i=0; i<SIZE(r);i++)
    {
      if (isvalidstring(MdsV0(r,i))>0)
        MdsV0(r_pattern,j++) = process_rlab_string_pattern(MdsV0(r,i));
    }
    if (!j)
    {
      mds_Destroy(r_pattern);
      r_pattern = 0;
    }
    else if (j != SIZE(l))
      mds_Extend(r_pattern,1,j);
  }

  s = mds_Create(MNR(t), MNC(t));

  int istart;
  int iend;
  for (i=0; i<SIZE(t); i++)
  {
    char * str = MdsV0(t, i);
    iend = isvalidstring(str);
    if (iend<1)
    {
      MdsV0 (s, i) = cpstr("");
      continue;
    }
    istart = 0;
    if (l_pattern)
      istart = processed_pattern_starts_string  (str, l_pattern);
    if (r_pattern)
      iend   = does_processed_pattern_end_string(str, r_pattern);

    MdsV0 (s, i) = cpnstr(str+istart,iend-istart);
  }

_exit_strip:

  if (l_pattern)
    mds_Destroy(l_pattern);
  if (r_pattern)
    mds_Destroy(r_pattern);

  ent_Clean (e1);
  ent_Clean (e2);
  return ent_Assign_Rlab_MDS(s);
}

#undef THIS_SOLVER
#define THIS_SOLVER "argv"
extern int    rlab_argc;
extern char **rlab_argv;
Ent * Argv (int nargs, Datum args[])
{
  MDS *s = 0;
  int i;

  // Get the first argument */
  if (rlab_argc>0)
  {
    s = mds_Create(1,rlab_argc);
    for (i=0; i<rlab_argc; i++)
      MdsV0(s,i) = cpstr(rlab_argv[i]);
  }

  return ent_Assign_Rlab_MDS(s);
}

#undef THIS_SOLVER
#define THIS_SOLVER "interactive"
extern int interactive;
Ent * Interactive (int nargs, Datum args[])
{
  return ent_Create_Rlab_Double(interactive);
}

#undef THIS_SOLVER
#define THIS_SOLVER "basename"
extern char * progname;
Ent * Basename (int nargs, Datum args[])
{
  return ent_Create_Rlab_String(progname);
}

#undef THIS_SOLVER
#define THIS_SOLVER "assign"
Ent * ent_Assign(int nargs, Datum args[])
{
  Ent *e1=0, *e2=0;
  if (nargs!=2)
  {
    return ent_Create_Rlab_Failure();
  }
  e2 = bltin_get_ent (args[1]);
  if ((args[0].type != VAR)||(ent_type(e2)==UNDEF))
  {
    ent_Clean (e2);
    return ent_Create_Rlab_Failure();
  }

  ListNode *var = (ListNode *) (args[0].u.ptr);
  e1 = var_ent (var);
//   if (ent_type(e1)!=UNDEF)
//   {
//     if (e1->refc > 1)
//     {
//       ent_DecRef (e1);
//     }
//   }
//   listNode_AttachEnt (var, ent_Copy(e2));

//   if (ent_type(e1)!=UNDEF)
//   {
//     if (e1->refc > 1)
//     {
//       ent_DecRef (e1);
//       e1 = ent_Duplicate (e1);
//       listNode_AttachEnt (var, e1);
//     }
//   }

//   ent_data(e1) = ent_data(e2);
//   ent_type(e1) = ent_type(e2);
//   ent_IncRef(e2);

  ent_data_class_copy(e1,e2);
  ent_Clean (e1);
  ent_Clean (e2);
  return ent_Create_Rlab_Success();
}

#include "./clibs/jsmn/jsmn.h"

#define JSON_MAX_TOKENS 128
static Btree * rlab_unpack_jsmn_object(char * entire_parse_string, int start_parse, int end_parse);
static void * rlab_unpack_jsmn_array(char * entire_parse_string, int start_parse, int end_parse, int *rtype)
{
  void * rval=0;
  *rtype = UNDEF;
  int rtype2=UNDEF;

  int i, num_top_objects, j, count=0;
  jsmn_parser p;
  jsmntok_t t[JSON_MAX_TOKENS]; /* We expect no more than 128 tokens */

  char * parse_string = cpnstr(entire_parse_string+start_parse, end_parse - start_parse);

  jsmn_init(&p);
  num_top_objects = jsmn_parse(&p, parse_string, isvalidstring(parse_string), t, JSON_MAX_TOKENS);

  if (num_top_objects<0)
    return NULL;

  char * datum_s=0;
  double datum_d=0.0;

  // skip root object
  i = 1;
  // first object determines what kind of array it is:
  //  MDS
  //  MDR
  //  MDE
  while (i<num_top_objects)
  {
    switch (t[i].type)
    {
      case JSMN_PRIMITIVE:
        if (*rtype == UNDEF)
        {
          count = 0;

          // get the data point
          *rtype = MATRIX_DENSE_REAL;
          rval = (void *) mdr_Create(1, t[0].size);
        }
        // get data point
        if (parse_string[t[i].start] == 't' || parse_string[t[i].start] == 'T')
        {
          datum_d = 1.0;
        }
        else if (parse_string[t[i].start] == 'f' || parse_string[t[i].start] == 'F')
        {
          rval = (void *) mdr_Create(1, t[0].size);
          datum_d = 0.0;
        }
        else if (parse_string[t[i].start] == '-' || (parse_string[t[i].start] >= '0' && parse_string[t[i].start] <= '9'))
        {
          char *who_cares;
          datum_d = strtod(parse_string+t[i].start, &who_cares);
        }
        MdrV0(rval,count) = datum_d;
        count++;
        break;

      case JSMN_OBJECT:
        if (*rtype == UNDEF)
        {
          count = 0;

          // get the data point
          *rtype = MATRIX_DENSE_ENTITY;
          rval = (void *) mde_Create(1,t[0].size);
        }
        // do something
        MdeV0(rval,count) = (Ent *) ent_Assign_Rlab_BTREE ((Btree*) rlab_unpack_jsmn_object(parse_string, t[i].start, t[i].end));
        count++;
        // skip all objects belonging to this one just examined
        j = i + 1;
        while (t[j].end < t[i].end) j++;
        i = j;
        continue;
        break;

      case JSMN_ARRAY:
        rval = (void *) rlab_unpack_jsmn_array(parse_string, t[i].start, t[i].end, &rtype2);
        MdeV0(rval,count) = (Ent *) ent_Assign_Rlab_Rtype(rval,rtype2);
        count++;
        j = i + 1;
        while (t[j].end < t[i].end) j++;
        i = j;
        continue;
        break;

      case JSMN_STRING:
        if (*rtype == UNDEF)
        {
          count = 0;

          // get the data point
          *rtype = MATRIX_DENSE_STRING;
          rval = (void *) mds_Create(1, t[0].size);
        }
        if (t[i].end > t[i].start)
          datum_s=cpnstr(parse_string+t[i].start, t[i].end - t[i].start);
        else
          datum_s = cpstr("(null)");
        MdsV0(rval,count) = datum_s;
        count++;
        break;

    } // switch (t[i].type)
    i++;

  } // while (i<num_top_objects)


  GC_free(parse_string);

  return (rval);
}

static Btree * rlab_unpack_jsmn_object(char * entire_parse_string, int start_parse, int end_parse)
{
  int i, num_top_objects, j;
  jsmn_parser p;
  jsmntok_t t[JSON_MAX_TOKENS]; /* We expect no more than 128 tokens */

  char * parse_string = cpnstr(entire_parse_string+start_parse, end_parse - start_parse);

  jsmn_init(&p);

  num_top_objects = jsmn_parse(&p, parse_string, isvalidstring(parse_string), t, JSON_MAX_TOKENS);

  if (num_top_objects<0)
    return NULL;

  Btree *b=btree_Create();

  char *enode=0, *datum_s=0;
  int have_node_name=0; /* 0-completed node, 1-node known collecting data */
  double datum_d=0.0;

  // skip root object
  i = 1;
  while (i<num_top_objects)
  {
    switch (t[i].type)
    {
      case JSMN_PRIMITIVE:
        if (have_node_name)
        {
          // get the data point
          if (parse_string[t[i].start] == 't' || parse_string[t[i].start] == 'T')
          {
            datum_d = 1.0;
          }
          else if (parse_string[t[i].start] == 'f' || parse_string[t[i].start] == 'F')
          {
            datum_d = 0.0;
          }
          else if (parse_string[t[i].start] == '-' || (parse_string[t[i].start] >= '0' && parse_string[t[i].start] <= '9'))
          {
            char *who_cares;
            datum_d = strtod(parse_string+t[i].start, &who_cares);
          }

          // what to do with it?
          install(b,enode,ent_Create_Rlab_Double (datum_d));
          have_node_name = 0;
          GC_free(enode);
        }
        break;

      case JSMN_OBJECT:
        if (have_node_name)
        {
          if (t[i].size==0)
          {
            install(b,enode,ent_Assign_Rlab_BTREE(NULL));
            have_node_name = 0;
            GC_free(enode);
            break;
          }
        }
        install(b,enode,ent_Assign_Rlab_BTREE((Btree*) rlab_unpack_jsmn_object(parse_string, t[i].start, t[i].end)));
        have_node_name = 0;
        GC_free(enode);
        j = i + 1;
        while (t[j].end < t[i].end) j++;
        i = j;
        continue;
        break;

      case JSMN_ARRAY:
        if (have_node_name)
        {
          if (t[i].size==0)
          {
            install(b,enode,ent_Assign_Rlab_MDS(NULL));
            have_node_name = 0;
            GC_free(enode);
            break;
          }
        }
        int rtype;
        void * rval = rlab_unpack_jsmn_array(parse_string, t[i].start, t[i].end, &rtype);
        install(b,enode,ent_Assign_Rlab_Rtype(rval,rtype));
        have_node_name = 0;
        GC_free(enode);

        j = i + 1;
        while (t[j].end < t[i].end) j++;
        i = j;
        continue;
        break;

      case JSMN_STRING:
        if (!have_node_name)
        {
          // have name of the node
          enode=cpnstr(parse_string+t[i].start, t[i].end - t[i].start);
          have_node_name = 1;
        }
        else
        {
          if (t[i].end > t[i].start)
            datum_s=cpnstr(parse_string+t[i].start, t[i].end - t[i].start);
          else
            datum_s = cpstr("(null)");

          // have string datum that belongs to the node
          install(b,enode,ent_Assign_Rlab_String (datum_s));
          have_node_name = 0;
          GC_free(enode);
        }
        break;

      default:
        break;

    } // switch

    i++;

  } // while (i<num_top_objects)

  GC_free(parse_string);

  return (b);
}

#undef  THIS_SOLVER
#define THIS_SOLVER "unpack"
Ent * ent_Unpack(int nargs, Datum args[])
{
  int i,what=1; /*1-json, 2-xml*/
  char * method=0, *parse_string=0;
  Ent *e1=0, *e2=0;

  if (nargs<1 || nargs>2)
  {
    return ent_Create_Rlab_Failure();
  }

  // process first argument
  RLABCODE_PROCESS_ARG1_S(THIS_SOLVER,0,e1,parse_string,2);

  if (nargs == 2)
  {
    RLABCODE_PROCESS_ARG_S(THIS_SOLVER,1,e2,method,1,RLAB_ERROR_ARG2_MDS_SCALAR);
    if (method[0] == 'x' || method[0] == 'X')
      what = 2;
  }

  Ent *rval=0;
  Btree * b = 0;
  void *r=0;
  int rtype=UNDEF;

  i=0;
  while (parse_string[i] == ' ') i++;

  if (what == 1)
  {
    // root object is either 
    //    json object
    // or
    //    json array
    if (parse_string[i] == '{')
    {
      b = (Btree *) rlab_unpack_jsmn_object(parse_string, i, isvalidstring(parse_string));
      rval = (Ent *) ent_Assign_Rlab_BTREE (b);
    }
    else if (parse_string[i] == '[')
    {
      r = (void *) rlab_unpack_jsmn_array(parse_string, i, isvalidstring(parse_string), &rtype);
      rval = (Ent *) ent_Assign_Rlab_Rtype (r, rtype);
    }
  }

  ent_Clean (e1);
  ent_Clean (e2);

  return rval;
}

