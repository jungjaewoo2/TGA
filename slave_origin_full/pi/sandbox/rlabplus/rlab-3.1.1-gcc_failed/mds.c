/* mds.c Matrix Dense String */

/*  This file is a part of RLaB ("Our"-LaB)
   Copyright (C) 1995  Ian R. Searle
   Copyright (C) 2017  Marijan Kostrun

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
#include "btree.h"
#include "mem.h"
#include "mds.h"
#include "mdr.h"
#include "util.h"
#include "bltin1.h"
#include "rfileio.h"
#include "mathl.h"
#include "rlab_solver_parameters_names.h"

#include <stdio.h>
#include <math.h>

#include <termios.h>
#include <sys/fcntl.h>
#include <unistd.h>
#include <sys/ioctl.h>

//
// buffer for all string operations:
//
char string_buff[MAX_STRING_BUFF];

//
// string operations that are always useful to have
//
int isvalidstring(char *s)
{
  if (!s)
    return -1;
  return strlen(s);
}

void swap2strings (char **s1, char **s2)
{
  char *tmp=*s1;
  *s1 = *s2;
  *s2 = tmp;
  return;
}

char * add2strings (char *s1, char *s2)
{
  size_t len2=0, len1=0;
  char *result = 0;
  if (s1)
    len1 = strlen(s1);
  if (s2)
    len2 = strlen(s2);
  result = GC_MALLOC(len1+len2+1);
  if (len1)
    memcpy(result, s1, len1);
  if (len2)
    memcpy(result+len1, s2, len2);
  result[len1+len2] = '\0';

  return result;
}

void add2strings2first (char **s1, char *s2)
{
  size_t len2=0, len1=0;
  char *first  = *s1;
  char *result = 0;
  if (first)
    len1 = strlen(first);
  if (s2)
    len2 = strlen(s2);
  result = GC_MALLOC(len1+len2+1);
  if (len1)
    memcpy(result, first, len1);
  if (len2)
    memcpy(result+len1, s2, len2);
  result[len1+len2] = '\0';

  if (*s1)
    GC_FREE(*s1);

  *s1 = result;
  return;
}

void add2strings2last (char *s1, char **s2)
{
  size_t len2=0, len1=0;
  char *result = 0;
  char *last   = *s2;
  if (s1)
    len1 = strlen(s1);
  if (last)
    len2 = strlen(last);
  result = GC_MALLOC(len1+len2+1);
  if (len1)
    memcpy(result, s1, len1);
  if (len2)
    memcpy(result+len1, last, len2);
  result[len1+len2] = '\0';

  if (*s2)
    GC_FREE(*s2);

  *s2 = result;
  return;
}

char * add3strings (char *s1, char *s2, char *s3)
{
  size_t len3=0, len2=0, len1=0;
  char *result = 0;
  if (s1)
    len1 = strlen(s1);
  if (s2)
    len2 = strlen(s2);
  if (s3)
    len3 = strlen(s3);
  result = GC_MALLOC(len1+len2+len3+1);
  if (len1)
    memcpy(result, s1, len1);
  if (len2)
    memcpy(result+len1, s2, len2);
  if (len3)
    memcpy(result+len1+len2, s3, len3);
  result[len1+len2+len3] = '\0';
  return result;
}

void add3strings2first (char **s, char *s2, char *s3)
{
  size_t len3=0, len2=0, len1=0;
  char *result = 0;
  char *s1  = *s;
  if (s1)
    len1 = strlen(s1);
  if (s2)
    len2 = strlen(s2);
  if (s3)
    len3 = strlen(s3);
  result = GC_MALLOC(len1+len2+len3+1);
  if (len1)
    memcpy(result, s1, len1);
  if (len2)
    memcpy(result+len1, s2, len2);
  if (len3)
    memcpy(result+len1+len2, s3, len3);
  result[len1+len2+len3] = '\0';

  if (*s)
    GC_FREE(*s);

  *s = result;
  return;
}

void add3strings2last (char *s1, char *s2, char **s)
{
  size_t len3=0, len2=0, len1=0;
  char *result = 0;
  char *s3 = *s;
  if (s1)
    len1 = strlen(s1);
  if (s2)
    len2 = strlen(s2);
  if (s3)
    len3 = strlen(s3);
  result = GC_MALLOC(len1+len2+len3+1);
  if (len1)
    memcpy(result, s1, len1);
  if (len2)
    memcpy(result+len1, s2, len2);
  if (len3)
    memcpy(result+len1+len2, s3, len3);
  result[len1+len2+len3] = '\0';

  if (*s)
    GC_FREE(*s);

  *s = result;
  return;
}

char * mdsV0_safe(MDS *m, int k)
{
  if (SIZE(m)>0)
  {
    int i = MAX(MIN(SIZE(m)-1, k),0);
    return MdsV0(m, i);
  }
  return NULL;
}

char * mdsV1_safe(MDS *m, int k)
{
  if (SIZE(m)<1)
  {
    int i = MAX(MIN(SIZE(m)-1, k-1),0);
    return MdsV0(m, i);
  }
  return NULL;
}

char * mds0_safe(MDS *m, int k, int j)
{
  if (SIZE(m)>0)
  {
    int r = MAX(MIN(MNR(m)-1, k),0);
    int c = MAX(MIN(MNC(m)-1, j),0);
    return Mds0(m, r, c);
  }
  return NULL;
}

char * mds1_safe(MDS *m, int k, int j)
{
  if (SIZE(m)>0)
  {
    int r = MAX(MIN(MNR(m)-1, k-1),0);
    int c = MAX(MIN(MNC(m)-1, j-1),0);
    return Mds0(m, r, c);
  }
  return NULL;
}


void mds_Flip (MDS *m, int ud, int lr)
{
  if ((!ud) && (!lr))
    return;

  int nr = MNR(m);
  int nc = MNC(m);

  if ((nr<1)||(nc<1))
    return;

  int i, j;
  char *tmp=0;


  if (ud)
  {
    for (j=0; j<nc; j++)
    {
      for (i=0; i<(nr/2); i++)
      {
        tmp = Mds0(m,nr-i-1,j);
        Mds0(m,nr-i-1,j) = Mds0(m,i,j);
        Mds0(m,i,j) = tmp;
      }
    }
  }

  if (lr)
  {
    for (i=0; i<nr; i++)
    {
      for (j=0; j<(nc/2); j++)
      {
        tmp = Mds0(m,i,nc-j-1);
        Mds0(m,i,nc-j-1) = Mds0(m,i,j);
        Mds0(m,i,j) = tmp;
      }
    }
  }

  return;
}


void mds_Shift (MDS *m, int ud, int lr)
{
  if ((!ud) && (!lr))
    return;

  int nr = MNR(m);
  int nc = MNC(m);

  if ((nr<1)||(nc<1))
    return;

  int i, j;


  if (ud != 0)
  {
    // shift up: ud>0
    // or down: ud<0
    if ((ud >= nr)||(ud <= -nr))
    {
      for (i=0; i<SIZE(m); i++)
      {
        if (MdsV0(m,i))
        {
          GC_FREE (MdsV0(m,i));
          MdsV0(m,i)=0;
        }
      }
      return;
    }

    if (ud > 0)
    {
      for (j=0; j<nc; j++)
      {
        for (i=ud; i<nr; i++)
        {
          if (Mds0(m,i-ud,j))
          {
            GC_FREE (Mds0(m,i-ud,j));
            Mds0(m,i-ud,j) = 0;
          }
          Mds0(m,i-ud,j) = Mds0(m,i,j);
          Mds0(m,i,j) = 0;
        }
      }
    }
    if (ud < 0)
    {
      for (j=0; j<nc; j++)
      {
        for (i=nr+ud-1; i>=0; i--)
        {
          if (Mds0(m,i-ud,j))
          {
            GC_FREE(Mds0(m,i-ud,j));
            Mds0(m,i-ud,j) = 0;
          }
          Mds0(m,i-ud,j) = Mds0(m,i,j);
          Mds0(m,i,j) = 0;
        }
        for (i=0; i<-ud; i++)
        {
          if (Mds0(m,i,j))
          {
            GC_FREE(Mds0(m,i,j));
            Mds0(m,i,j) = 0;
          }
        }
      }
    }
  }

  if (lr != 0)
  {
    // shift up: ud>0
    // or down: ud<0
    if ((lr>= nc)||(lr <= -nc))
    {
      for (i=0; i<SIZE(m); i++)
      {
        if (MdsV0(m,i))
        {
          GC_FREE (MdsV0(m,i));
          MdsV0(m,i)=0;
        }
      }
      return;
    }

    if (lr > 0)
    {
      for (i=0; i<nr; i++)
      {
        for (j=lr; j<nc; j++)
        {
          if (Mds0(m,i,j-lr))
          {
            GC_FREE (Mds0(m,i,j-lr));
            Mds0(m,i,j-lr) = 0;
          }
          Mds0(m,i,j-lr) = Mds0(m,i,j);
          Mds0(m,i,j) = 0;
        }
      }
    }
    if (lr < 0)
    {
      for (i=0; i<nr; i++)
      {
        for (j=nc+lr-1; j>=0; j--)
        {
          if (Mds0(m,i,j-lr))
          {
            GC_FREE(Mds0(m,i,j-lr));
            Mds0(m,i,j-lr) = 0;
          }
          Mds0(m,i,j-lr) = Mds0(m,i,j);
          Mds0(m,i,j) = 0;
        }
        for (j=0; j<-lr; j++)
        {
          if (Mds0(m,i,j))
          {
            GC_FREE(Mds0(m,i,j));
            Mds0(m,i,j) = 0;
          }
        }
      }
    }
  }

  return;
}




/* **************************************************************
 * Create a matrix.
 * ************************************************************** */
MDS *
mds_Create (int nrow, int ncol)
{
  int i;
  MDS *new=0;

  if (nrow < 0 || ncol < 0)
    rerror ("cannot specify a negative matrix dimension");

  new = (MDS *) GC_MALLOC (sizeof (MDS));
  if (!new)
    rerror ("out of memory");

  MNR(new) = nrow;
  MNC(new) = ncol;
  new->type = RLAB_TYPE_STRING;

  if (nrow * ncol != 0)
  {
    MDSPTR(new) = (char **) GC_MALLOC (sizeof (char *) * (nrow * ncol));
    if (!MDSPTR(new))
      rerror ("out of memory");
    for (i = 0; i < nrow * ncol; i++)
    {
      MdsV0(new,i) = 0;
    }
  }
  else
  {
    new->nrow = 0;
    new->ncol = 0;
    MDSPTR(new) = 0;
  }
  new->list = 0;

  return (new);
}

MDS *
mds_CreateScalar (char *string)
{
  MDS *new = 0;

  if (string)
  {
    new = mds_Create (1, 1);
    MdsV0(new,0) = cpstr(string);
  }
  else
  {
    new = mds_Create (0, 0);
  }

  return (new);
}

/* **************************************************************
 * Free a matrix, and wipe out the structure members.
 * ************************************************************** */

void *
mds_Destroy (MDS * m)
{
  int i;

  if (!m)
    return (0);

  for (i=0; i<SIZE(m); i++)
  {
    if (MdsV0(m,i))
    {
      GC_FREE (MdsV0(m,i));
      MdsV0(m,i)=0;
    }
  }

  if (MDPTR(m))
    GC_FREE (MDPTR(m));
  MDPTR(m)=0;

  m->nrow = -1;
  m->ncol = -1;

  if (m->list)
  {
    btree_Destroy (m->list);
    m->list = 0;
  }

  GC_FREE (m);
  return (0);
}

/* **************************************************************
 * Copy a matrix. Create the new matrix, and return the new,
 * copied matrix as the result.
 * ************************************************************** */

MDS *
mds_Copy (MDS * m)
{
  int i, size;
  MDS *new = 0;

  new = mds_Create (MNR (m), MNC (m));

  size = MNR (m) * MNC (m);
  for (i = 0; i < size; i++)
  {
    if (MdsV0(m,i))
    {
      MdsV0(new,i) = cpstr (MdsV0(m,i));
    }
    else
    {
      MdsV0(new,i) = 0;
    }
  }

  /* Copy the list, if there is one. */
  if (m->list)
  {
    new->list = btree_Copy (m->list);
  }

  return (new);
}

/* **************************************************************
 * Reshape a matrix (change nrow, ncol).
 * ************************************************************** */

MDS *
mds_Reshape (MDS * m, int nrow, int ncol)
{
  MDS *new;
  if (nrow * ncol != MNR (m) * MNC (m))
  {
    fprintf (stderr, "\tincompatible dimensions for reshape\n");
    fprintf (stderr, "\tnrow*ncol must equal: %i\n", MNR (m) * MNC (m));
    rerror ("error");
  }

  new = mds_Copy (m);

  /* Now, reset the row and column sizes. */
  if ((MNR (new) == 0) && (MNC (new) == 0))
  {
    /* Leave it alone. */
    ;
  }
  else
  {
    /* Fix the row and column sizes. */
    new->nrow = nrow;
    new->ncol = ncol;
  }

  return (new);
}

/* **************************************************************
 * Extend a matrix. Realloc the data, modify the data size.
 *    just copy the pointers,
 *    delete pointers, and what they point to if necessary
 * ************************************************************** */
void mds_Extend (MDS *old, int new_nrow, int new_ncol)
{
  int i,j;
  MDS *new;

  if (SIZE(old)==new_nrow*new_ncol)
    return;

  new = mds_Create (new_nrow, new_ncol);
  if (!new)
    rerror ("out of memory");

  // now go over old matrix:
  //   copy pointers if they belong to new matrix
  //   or free locations they point to
  for (i=0; i<MNR(old); i++)
  {
    for (j=0; j<MNC(old); j++)
    {
      if ((i<new_nrow)&& (j<new_ncol))
        Mds0(new,i,j) = Mds0(old,i,j);
      else
        GC_FREE (Mds0(old,i,j));
      Mds0(old,i,j) = 0;
    }
  }

  // kill old->d
  GC_FREE(old->d);

  // replace with new one:
  old->d = new->d;

  // kill new matrix
  new->d = 0;
  new->nrow = 0;
  new->ncol = 0;
  mds_Destroy(new);

  // update size
  old->nrow = new_nrow;
  old->ncol = new_ncol;

  return;
}

MDS * mds_Extend_old (MDS * m, int nrow, int ncol)
{
  int i, j, nr_old, nc_old;
  char **old;

  if ((MNR (m) * MNC (m)) != 0)
  {
    if (nrow > MNR (m))
    {
      /* Save the original data, for now. */

      old = MDSPTR(m);
      nr_old = MNR (m);
      nc_old = MNC (m);

      /* Re-configure the matrix. */

      MDSPTR(m) = GC_MALLOC (nrow * nc_old * sizeof (char *));
      if (!MDSPTR(m))
        rerror ("out of memory");
      m->nrow = nrow;
      m->ncol = ncol;

      /* Copy the old into the new. */

      for (i = 0; i < nr_old; i++)
      {
        for (j = 0; j < nc_old; j++)
        {
          if (old[j * nr_old + i] != 0)
          {
            Mds0 (m, i, j) = cpstr (old[j * nr_old + i]);
            GC_FREE (old[j * nr_old + i]);
          }
          else
          {
            Mds0 (m, i, j) = 0;
          }
        }
      }

      /* Zero (sort-of) out the new extension. */

      for (i = nr_old; i < MNR (m); i++)
      {
        for (j = 0; j < MNC (m); j++)
        {
          Mds0 (m, i, j) = cpstr ("");
        }
      }

      GC_FREE (old);
    }

    if (ncol > MNC (m))
    {
      /* Save the original data, for now. */

      old = MDSPTR(m);
      nr_old = MNR (m);
      nc_old = MNC (m);

      /* Re-configure the matrix. */

      MDSPTR(m) = (char **) GC_MALLOC (nr_old * ncol * sizeof (char *));
      if (!MDSPTR(m))
        rerror ("out of memory");
      m->nrow = nrow;
      m->ncol = ncol;

      /* Copy the old into the new. */

      for (i = 0; i < nr_old; i++)
      {
        for (j = 0; j < nc_old; j++)
        {
          if (old[j * nr_old + i] != 0)
          {
            Mds0 (m, i, j) = cpstr (old[j * nr_old + i]);
            GC_FREE (old[j * nr_old + i]);
          }
          else
          {
            Mds0 (m, i, j) = 0;
          }
        }
      }

      /* Zero out the new extension. */

      for (i = 0; i < MNR (m); i++)
      {
        for (j = nc_old; j < MNC (m); j++)
        {
          Mds0 (m, i, j) = cpstr ("");
        }
      }

      GC_FREE (old);
    }
  }
  else
  {
    /* Re-configure the matrix. */

    MDSPTR(m) = (char **) GC_MALLOC (nrow * ncol * sizeof (char *));
    if (!MDSPTR(m))
      rerror ("out of memory");
    m->nrow = nrow;
    m->ncol = ncol;

    /* Zero (sort-of) out the new extension. */

    for (i = 0; i < MNR (m); i++)
    {
      for (j = 0; j < MNC (m); j++)
      {
        Mds0 (m, i, j) = cpstr ("");
      }
    }
  }
  return (m);
}

/* **************************************************************
 * Return a string matrix pointer.
 * ************************************************************** */

MDS *
mds_MatrixString (MDS * m)
{
  return (m);
}

/* **************************************************************
 * Print out a matrix.
 * ************************************************************** */
#include <curses.h>
#include <term.h>
#include <sys/ioctl.h>
#undef  THIS_SOLVER
#define THIS_SOLVER "class_print:mds_Print"
void mds_Print (MDS * m, FILE * stream)
{
  int i, j, k, length, nrow, ncol, ncol_print, npri, rem;
  int start, tmp, swidth;

  /* Special case, empty matrix */
  if (SIZE(m)<1)
  {
    fprintf (stream, "\t[]\n");
    fflush (stream);
    return;
  }

  // special case
  if (SIZE(m) == 1)
  {
    fprintf (stream, "%s", MdsV0 (m, 0));
    fflush (stream);
    return;
  }

  nrow = MNR (m);
  ncol = MNC (m);

  /* Figure out the length of the longest STRING */
  length = 0;
  for (i = 0; i < nrow * ncol; i++)
  {
    if(isvalidstring(MdsV0 (m, i))<1)
      continue;

    tmp = strlen (MdsV0 (m, i));
    if (tmp > length)
      length = tmp;
  }

  /* Now figure out how many columns we can print on a line */
  struct winsize sz;
  ioctl(0, TIOCGWINSZ, &sz);
  swidth = sz.ws_col;
  if (swidth <= 0)
    swidth = TERM_WIDTH;

  ncol_print = swidth / (length + 2);

  if (ncol_print == 0)
    ncol_print = 1;

  npri = MNC (m) / ncol_print;
  rem = MNC (m) % ncol_print;

  start = 1;
  for (i = 0; i < npri; i++)
  {
    for (k = 1; k <= nrow; k++)	/* print all rows */
    {
      for (j = start; j <= ncol_print + start - 1; j++)
      {
        fprintf (stream, "%-*s", length + 2, Mds1 (m, k, j));
      }
      fprintf (stream, "\n");
    }
    start += ncol_print;	/* inc our col position */
    fprintf (stream, "\n");
    fflush (stream);
  }

  /* Now come back and write out the last few collums */
  if (!rem)
    return;
  for (k = 1; k <= nrow; k++)
  {
    for (i = ncol - rem + 1; i <= ncol; i++)
    {
      fprintf (stream, "%-*s", length + 2, Mds1 (m, k, i));
    }
    fprintf (stream, "\n");
    fflush (stream);
  }
}

/* **************************************************************
 * Return a the size of the matrix.
 * ************************************************************** */

int
mds_Size (MDS * m)
{
  return (SIZE(m));
}

/* **************************************************************
 * Append two matrices together:  [ m1 , m2 ]
 * ************************************************************** */

MDS *
mds_Append (MDS * m1, MDS * m2)
{
  int i, j, nrow, ncol;
  MDS *new=0;

  /* Check for empty matrices... */
  if (SIZE(m1)<1 && SIZE(m2)<1)
  {
    new = mds_Create(0,0);
    return new;
  }
  else if (SIZE(m1)<1)
  {
    new = mds_Copy (m2);
    return (new);
  }
  else if (SIZE(m2)<1)
  {
    new = mds_Copy (m1);
    return (new);
  }

  /* Do the append. */
  /* Create the new matrix, large enough to hold both. */
  int r1 = MNR(m1);
  int c1 = MNC(m1);
  int r2 = MNR(m2);
  int c2 = MNC(m2);
  nrow = MAX(r1,r2);
  ncol = c1 + c2;
  new = mds_Create (nrow, ncol);

  for (i=0; i<nrow; i++)
  {
    for (j=0; j<c1; j++)
      Mds0(new,i,j)     = cpstr(mds0_safe(m1,i,j));
    for (j=0; j<c2; j++)
      Mds0(new,i,j+c1)  = cpstr(mds0_safe(m2,i,j));
  }

  return (new);
}

/* **************************************************************
 * Stack two matrices together:  [ m1 ; m2 ]
 * ************************************************************** */

MDS *
mds_Stack (MDS * m1, MDS * m2)
{
  int i, j, nrow, ncol;
  MDS *new=0;

  /* Check for empty matrices... */
  if (SIZE(m1)<1 && SIZE(m2)<1)
  {
    new = mds_Create(0,0);
    return new;
  }
  else if (SIZE(m1)<1)
  {
    new = mds_Copy (m2);
    return (new);
  }
  else if (SIZE(m2)<1)
  {
    new = mds_Copy (m1);
    return (new);
  }

  /* Do the stack. */
  /* Create the new matrix, large enough to hold both. */
  int r1 = MNR(m1);
  int c1 = MNC(m1);
  int r2 = MNR(m2);
  int c2 = MNC(m2);
  nrow = r1+r2;
  ncol = MAX(c1,c2);
  new = mds_Create (nrow, ncol);

  for (j=0; j<ncol; j++)
  {
    for (i=0; i<r1; i++)
      Mds0(new,i,j)     = cpstr( mds0_safe(m1,i,j) );
    for (i=0; i<r2; i++)
      Mds0(new,i+r1,j)  = cpstr( mds0_safe(m2,i,j) );
  }

  return (new);
}

MDS *
mds_MatrixSub (MDS * var, int *i, int *j, int *type)
{
  int m, n;
  MDS *new;

  *type = MATRIX_DENSE_STRING;

  /* Handle empty matrix indices. */

  if ((i == 0) || (j == 0))
  {
    new = mds_Create (0, 0);
    return (new);
  }

  new = mds_Create (i[0], j[0]);

  for (m = 1; m <= i[0]; m++)
  {
    if (i[m] > MNR (var))
    {
      mds_Destroy (new);
      fprintf (stderr, "\tindex exceeds matrix limits\n");
      fprintf (stderr, "\tindex value: %i, matrix size: %i\n", i[m], MNR (var));
      rerror ("sub-matrix evaluation");
    }
    for (n = 1; n <= j[0]; n++)
    {
      if (j[n] > MNC (var))
      {
	mds_Destroy (new);
	fprintf (stderr, "\tindex exceeds matrix limits\n");
	fprintf (stderr, "\tindex value: %i, matrix size: %i\n", j[n],
		 MNC (var));
	rerror ("sub-matrix evaluation");
      }
      Mds1 (new, m, n) = cpstr (Mds1 (var, i[m], j[n]));
    }
  }

  return (new);
}

MDS *
mds_MatrixSubR (MDS * var, int *i, int *type)
{
  int m, n;
  MDS *new;

  *type = MATRIX_DENSE_STRING;

  /* Handle empty matrix indices. */

  if (i == 0)
  {
    new = mds_Create (0, 0);
    return (new);
  }

  new = mds_Create (i[0], MNC (var));

  for (m = 1; m <= i[0]; m++)
  {
    if (i[m] > MNR (var))
    {
      mds_Destroy (new);
      fprintf (stderr, "index exceeds matrix limits\n");
      fprintf (stderr, "index value: %i, matrix size: %i\n", i[m], MNR (var));
      rerror ("sub-matrix evaluation");
    }
    for (n = 1; n <= MNC (var); n++)
    {
      Mds1 (new, m, n) = cpstr (Mds1 (var, i[m], n));
    }
  }

  return (new);
}

MDS *
mds_MatrixSubC (MDS * var, int *j, int *type)
{
  int m, n;
  MDS *new;

  *type = MATRIX_DENSE_STRING;

  /* Handle empty matrix indices. */

  if (j == 0)
  {
    new = mds_Create (0, 0);
    return (new);
  }

  new = mds_Create (MNR (var), j[0]);

  for (m = 1; m <= MNR (var); m++)
  {
    for (n = 1; n <= j[0]; n++)
    {
      if (j[n] > MNC (var))
      {
	mds_Destroy (new);
	fprintf (stderr, "index exceeds matrix limits\n");
	fprintf (stderr, "index value: %i, matrix size: %i\n", j[n], MNC (var));
	rerror ("sub-matrix evaluation");
      }
      Mds1 (new, m, n) = cpstr (Mds1 (var, m, j[n]));
    }
  }

  return (new);
}

/* **************************************************************
 * Assign to a range of a matrix. We do not create a new matrix.
 * Automatically extend the size of a matrix if I or J indices
 * exceed current bounds.
 * ************************************************************** */

MDS *
mds_MatrixAssign (MDS * var, int *i, int *j, MDS * rhs)
{
  int m, n;
  int dflag;
  MDS *mtmp;

  if ((i == 0) || (j == 0))
  {
    if (MdsNR (rhs) * MdsNC (rhs) != 0)
    {
      rerror ("matrix-assign: LHS and RHS must have same dimensions");
    }
    else if (var != 0)
    {
      return (var);
    }
    else
    {
      rerror
	("matrix-assign: cannot assign to an UNDEFINED object without index specification");
    }
  }

  /* Check the LHS for UNDEF. */

  if (var == 0)
  {
    var = mds_Create (1, 1);
    Mds0 (var, 0, 0) = cpstr ("");
  }

  /* Check RHS for empty matrix. */

  if (MNR (rhs) * MNC (rhs) == 0)
  {
    fprintf (stderr, "cannot assign empty matrix to element(s) of a matrix\n");
    rerror ("matrix-assign");
  }

  /* Check the LHS, and RHS dimensions. */

  if (MNR (var) != 0)
  {
    if (i[0] != MNR (rhs))
    {
      fprintf (stderr, "LHS, and RHS row dimensions must match\n");
      fprintf (stderr, "LHS row: %i, RHS row: %i\n", i[0], MNR (rhs));
      rerror ("matrix-assign");
    }

    if (j[0] != MNC (rhs))
    {
      fprintf (stderr, "LHS, and RHS column dimensions must match\n");
      fprintf (stderr, "LHS column: %i, RHS column: %i\n", j[0], MNC (rhs));
      rerror ("matrix-assign");
    }
  }
  else
  {
    mds_Extend (var, 1, 1);
  }

  /* Now do the assign. */

  if (var == rhs)
  {
    mtmp = mds_Copy (rhs);
    dflag = 1;
  }
  else
  {
    mtmp = rhs;
    dflag = 0;
  }

  for (m = 1; m <= i[0]; m++)
  {
    if (i[m] > MNR (var))
    {
      mds_Extend (var, i[m], MNC (var));
    }
    for (n = 1; n <= j[0]; n++)
    {
      if (j[n] > MNC (var))
      {
	mds_Extend (var, MNR (var), j[n]);
      }
      if (Mds1 (var, i[m], j[n]) != 0)
	GC_FREE (Mds1 (var, i[m], j[n]));
      Mds1 (var, i[m], j[n]) = cpstr (Mds1 (mtmp, m, n));
    }
  }

  if (dflag)
    mds_Destroy (mtmp);

  return (var);
}

MDS *
mds_MatrixAssignR (MDS * var, int *i, MDS * rhs)
{
  int m, n;
  int dflag;
  MDS *mtmp;

  if (i == 0)
  {
    if (MdsNR (rhs) * MdsNC (rhs) != 0)
    {
      rerror ("matrix-assign: LHS and RHS must have same dimensions");
    }
    else if (var != 0)
    {
      return (var);
    }
    else
    {
      rerror
	("matrix-assign: cannot assign to an UNDEFINED object without index specification");
    }
  }

  /* Check the LHS. */

  if (var == 0)
  {
    var = mds_Create (1, 1);
    Mds0 (var, 0, 0) = cpstr ("");
  }

  /* Check RHS for empty matrix. */

  if (MNR (rhs) * MNC (rhs) == 0)
  {
    fprintf (stderr, "cannot assign empty matrix to element(s) of a matrix\n");
    rerror ("matrix-assign");
  }

  /* Check the LHS, and RHS dimensions. */

  if (MNR (var) != 0)
  {
    if (i[0] != MNR (rhs))
    {
      fprintf (stderr, "LHS, and RHS row dimensions must match\n");
      fprintf (stderr, "LHS row: %i, RHS row: %i\n", i[0], MNR (rhs));
      rerror ("matrix-assign");
    }
  }
  else
  {
    mds_Extend (var, 1, 1);
  }

  /* Now do the assign. */

  if (var == rhs)
  {
    mtmp = mds_Copy (rhs);
    dflag = 1;
  }
  else
  {
    mtmp = rhs;
    dflag = 0;
  }

  for (m = 1; m <= i[0]; m++)
  {
    if (i[m] > MNR (var))
    {
      mds_Extend (var, i[m], MNC (var));
    }
    for (n = 1; n <= MNC (rhs); n++)
    {
      if (n > MNC (var))
      {
	mds_Extend (var, MNR (var), n);
      }
      if (Mds1 (var, i[m], n) != 0)
	GC_FREE (Mds1 (var, i[m], n));
      Mds1 (var, i[m], n) = cpstr (Mds1 (mtmp, m, n));
    }
  }

  if (dflag)
    mds_Destroy (mtmp);

  return (var);
}

MDS *
mds_MatrixAssignC (MDS * var, int *j, MDS * rhs)
{
  int m, n;
  int dflag;
  MDS *mtmp;

  if (j == 0)
  {
    if (MdsNR (rhs) * MdsNC (rhs) != 0)
    {
      rerror ("matrix-assign: LHS and RHS must have same dimensions");
    }
    else if (var != 0)
    {
      return (var);
    }
    else
    {
      rerror
	("matrix-assign: cannot assign to an UNDEFINED object without index specification");
    }
  }

  /* Check the LHS. */

  if (var == 0)
  {
    var = mds_Create (1, 1);
    Mds0 (var, 0, 0) = cpstr ("");
  }

  /* Check RHS for empty matrix. */

  if (MNR (rhs) * MNC (rhs) == 0)
  {
    fprintf (stderr, "cannot assign empty matrix to element(s) of a matrix\n");
    rerror ("matrix-assign");
  }

  /* Check the LHS, and RHS dimensions. */

  if (MNR (var) != 0)
  {
    if (j[0] != MNC (rhs))
    {
      fprintf (stderr, "LHS, and RHS column dimensions must match\n");
      fprintf (stderr, "LHS column: %i, RHS column: %i\n", j[0], MNC (rhs));
      rerror ("matrix-assign");
    }
  }
  else
  {
    mds_Extend (var, 1, 1);
  }

  /* Now do the assign. */

  if (var == rhs)
  {
    mtmp = mds_Copy (rhs);
    dflag = 1;
  }
  else
  {
    mtmp = rhs;
    dflag = 0;
  }

  for (m = 1; m <= MNR (rhs); m++)
  {
    if (m > MNR (var))
    {
      mds_Extend (var, m, MNC (var));
    }
    for (n = 1; n <= j[0]; n++)
    {
      if (j[n] > MNC (var))
      {
	mds_Extend (var, MNR (var), j[n]);
      }
      if (Mds1 (var, m, j[n]) != 0)
	GC_FREE (Mds1 (var, m, j[n]));
      Mds1 (var, m, j[n]) = cpstr (Mds1 (mtmp, m, n));
    }
  }

  if (dflag)
    mds_Destroy (mtmp);

  return (var);
}

/* **************************************************************
 * Vector Sub-Expression
 * ************************************************************** */

MDS *
mds_VectorSub (MDS * m, int *i, int *type)
{
  int j, size, msize;
  MDS *new;

  *type = MATRIX_DENSE_STRING;

  /* Handle empty matrix indices. */

  if (i == 0)
  {
    new = mds_Create (0, 0);
    return (new);
  }

  size = i[0];
  msize = (m->nrow) * (m->ncol);
  if (!msize)
  {
    new = mds_Create (0, 0);
    return (new);
  }

  if (MNC (m) == 1)
    new = mds_Create (size, 1);
  else
    new = mds_Create (1, size);

  for (j = 1; j <= size; j++)
  {
    if (i[j] > msize)
    {
      mds_Destroy (new);
      fprintf (stderr, "\tindex exceeds matrix limits\n");
      fprintf (stderr, "\tindex value: %i, matrix size: %i\n", i[j], msize);
      rerror ("sub-matrix evaluation");
    }
    MdsV1 (new, j) = cpstr (MdsV1 (m, i[j]));
  }

  return (new);
}

/* **************************************************************
 * Assign into a matrix like a vector.
 * var[j] = rhs;
 * ************************************************************** */

MDS *
mds_VectorAssign (MDS * var, int *j, MDS * rhs)
{
  int i, size;
  int nr, nc;
  int dflag;
  MDS *mtmp;

  if (j == 0)
  {
    if (MdsNR (rhs) * MdsNC (rhs) != 0)
    {
      rerror ("vector-assign: LHS and RHS must have same dimensions");
    }
    else if (var != 0)
    {
      return (var);
    }
    else
    {
      rerror
	("vector-assign: cannot assign to an UNDEFINED object without index specification");
    }
  }

  /* Check the LHS. */

  if (var == 0)
  {
    var = mds_Create (1, 1);
    Mds0 (var, 0, 0) = cpstr ("");
  }

  /* Check RHS for empty matrix. */

  if (MNR (rhs) * MNC (rhs) == 0)
  {
    fprintf (stderr,
	     "\tcannot assign empty matrix to element(s) of a matrix\n");
    rerror ("vector-assign");
  }

  /* Check the dimensions. */

  if (j[0] != (MNR (rhs) * MNC (rhs)))
  {
    fprintf (stderr, "\tLHS, and RHS column sizes must match\n");
    fprintf (stderr, "\tLHS: %i, RHS: %i\n", j[0], MNR (rhs) * MNC (rhs));
    rerror ("vector-assign");
  }

  /* Now do the assign. */

  if (var == rhs)
  {
    mtmp = mds_Copy (rhs);
    dflag = 1;
  }
  else
  {
    mtmp = rhs;
    dflag = 0;
  }

  size = j[0];

  /* Check for empty matrix. */
  if (MNR (var) == 0)
    mds_Extend (var, 1, 1);

  for (i = 1; i <= size; i++)
  {
    if (j[i] > MNR (var) * MNC (var))
    {
      if (MNR (var) == 1)
      {
	nr = 1;
	nc = (int) ceil (((double) j[i]) / ((double) MNR (var)));
      }
      else if (MNC (var) == 1)
      {
	nr = (int) ceil (((double) j[i]) / ((double) MNC (var)));
	nc = 1;
      }
      else
      {
	nr = MNR (var);
	nc = (int) ceil (((double) j[i]) / ((double) MNR (var)));
      }
      mds_Extend (var, nr, nc);
    }

    MdsV1 (var, (j[i])) = cpstr (MdsV1 (mtmp, i));
  }

  if (dflag)
    mds_Destroy (mtmp);

  return (var);
}

MDS * mds_ForLoopValue (MDS * m, int i)
{
  if (SIZE(m) >= i)
  {
    if (isvalidstring(MdsV1 (m, i)))
      return mds_CreateScalar(MdsV1 (m, i));
  }
  return mds_CreateScalar("");
}

char * mds_CharPointer (MDS * m)
{
  if ( (m->nrow)*(m->ncol) > 0)
    return (Mds0 (m, 0, 0));
  else
    return NULL;
}

char * mds_GetString (MDS * m)
{
  if ( (m->nrow)*(m->ncol) > 0)
    return (Mds0 (m, 0, 0));
  else
    return NULL;
}

/* **************************************************************
 * Make an array of strings into a MATRIX_DENSE_STRING.
 * ************************************************************** */

MDS *
mds_MakeCharMatrix (char **array, int n)
{
  MDS *new;

  new = mds_Create (0, 0);
  if (n != 0)
  {
    MNR(new) = 1;
    MNC(new) = n;
    MDSPTR(new) = array;
  }
  return (new);
}

/* **************************************************************
 * Return the class of a matrix-dense-string.
 * ************************************************************** */

char *
mds_Class (MDS * m)
{
  return (cpstr (RLAB_CLASS_STRING));
}

/* **************************************************************
 * Transpose a matrix.
 * ************************************************************** */

MDS * mds_Transpose (MDS * m)
{
  int i, j;
  MDS *new = mds_Create (MNC (m), MNR (m));

  for (i = 1; i <= MNR (m); i++)
  {
    for (j = 1; j <= MNC (m); j++)
    {
      Mds1 (new, j, i) = cpstr (Mds1 (m, i, j));
    }
  }

  return (new);
}

void * mds_Transpose_inplace (MDS * m)
{
  if (SIZE(m)<1)
    return NULL;

  int nr = MNR(m), nc=MNC(m);
  if (nr>1 && nc>1)
  {
    md_transpose_insitu((unsigned char *)MDSPTR(m), nr, nc, sizeof(char*));
  }
  MNR(m) = nc;
  MNC(m) = nr;

  return (0);
}

/* **************************************************************
 * Reshape a matrix into a colum vector.
 * ************************************************************** */

MDS *
mds_ReshapeCol (MDS * m)
{
  MDS *new=NULL;
  if (SIZE(m)>0)
  {
    new = mds_Reshape (m, SIZE(m), 1);
  }
  return (new);
}

/* **************************************************************
 * Add two string-matrices
 * ************************************************************** */
MDS *
mds_Add (MDS * m1, MDS * m2)
{
  MDS *new;
  int i, j, i1, i2, j1, j2;

  i1 = m1->nrow;
  j1 = m1->ncol;
  i2 = m2->nrow ;
  j2 = m2->ncol;

  // Quazi-emulate Matlab empty matrix behavior if both matrices are empty
  if ( i1 * j1 == 0 )
  {
    new = mds_Copy (m2);
    return (new);
  }
  if ( i2 * j2 == 0 )
  {
    new = mds_Copy (m1);
    return (new);
  }

  // do the matrix optimization thingy
  int nr = MAX(i1, i2);
  int nc = MAX(j1, j2);
  i1--;
  i2--;
  j1--;
  j2--;

  new = mds_Create(nr, nc);
  for(i=0;i<nr;i++)
    for(j=0;j<nc;j++)
      Mds0(new,i,j) = string_add ( mds0_safe(m1, i, j), mds0_safe(m2, i, j) );

  return (new);
}

MDS *
mds_Add_old (MDS * m1, MDS * m2)
{
  MDS *new, *m;
  char *s;
  int i, nr, nc, size;

  /*
   * Check sizes.
   * Check for special case (scalar + matrix) first.
   */

  if (((MNR (m1) == 1) && (MNC (m1) == 1)) ||
      ((MNR (m2) == 1) && (MNC (m2) == 1)))
  {
    if ((MNR (m1) == 1) && (MNC (m1) == 1))
    {
      s = Mds0 (m1, 0, 0);
      nr = MNR (m2);
      nc = MNC (m2);
      m = m2;

      new = mds_Create (nr, nc);
      size = nr * nc;

      for (i = 0; i < size; i++)
      {
	MdsV0 (new, i) = string_add (s, MdsV0 (m, i));
      }
    }
    else
    {
      s = Mds0 (m2, 0, 0);
      nr = MNR (m1);
      nc = MNC (m1);
      m = m1;

      new = mds_Create (nr, nc);
      size = nr * nc;

      for (i = 0; i < size; i++)
      {
	MdsV0 (new, i) = string_add (MdsV0 (m, i), s);
      }
    }

    return (new);
  }

  if (MNR (m1) != MNR (m2))
  {
    fprintf (stderr, "matrix row sizes must be equal\n");
    fprintf (stderr, "matrix row sizes: %i and %i\n", MNR (m1), MNR (m2));
    rerror ("matrix-addition: row size mis-match");
  }

  if (MNC (m1) != MNC (m2))
  {
    fprintf (stderr, "matrix column sizes must be equal\n");
    fprintf (stderr, "matrix column sizes: %i and %i\n", MNC (m1), MNC (m2));
    rerror ("matrix-addition: column size mis-match");
  }

  new = mds_Create (MNR (m1), MNC (m1));
  size = MNR (m1) * MNC (m1);

  for (i = 0; i < size; i++)
  {
    MdsV0 (new, i) = string_add (MdsV0 (m1, i), MdsV0 (m2, i));
  }

  return (new);
}


//
// writing a string matrix to serial port
//
extern int RLABPLUS_SERIAL_2BYTETIME_US;

void
mds_intfd_WriteGeneric (char *name, MDS * m)
{
  int i,n;
  int fd = get_int_file_ds(name);

  if (fd < 0)
  {
    printf ("writem: attempt to write to non-existent file/port\n");
    return;
  }

  if (!m)
  {
    printf ("writem: attempt to write non-existent string matrix to file/port\n");
    return;
  }

  // flush the port before trying to write to it
  usleep(10 * RLABPLUS_SERIAL_2BYTETIME_US);  // sleep 1 msec
  tcflush (fd, TCIOFLUSH);

  n = m->nrow * m->ncol;
  if (!n)
  {
    printf ("writem: attempt to write zero-size string matrix to file/port\n");
    return;
  }

  char *eol = get_eol_file_name (name);

  for (i=0; i<n; i++)
  {
    while(write(fd, MdsV0(m,i), strlen(MdsV0(m,i))) < strlen(MdsV0(m,i)))
    {
      usleep(RLABPLUS_SERIAL_2BYTETIME_US);// sleep a little
      tcflush (fd, TCIOFLUSH);          // flush the buffer and try again
    }

    if (eol)
      write(fd, eol, strlen(eol));

    usleep(RLABPLUS_SERIAL_2BYTETIME_US);  // sleep 1 msec
    tcflush (fd, TCIOFLUSH);            // flush the buffer before trying again
  }
}


int
mds_Sublist (MDS * m, char *name)
{
  if (!strcmp (name, RLAB_MEMBER_NROW))
  {
    return (1);
  }
  else if (!strcmp (name, RLAB_MEMBER_NCOL))
  {
    return (1);
  }
  else if (!strcmp (name, RLAB_MEMBER_SIZE))
  {
    return (1);
  }
  else if (!strcmp (name, RLAB_MEMBER_CLASS))
  {
    return (1);
  }
  else if (!strcmp (name, RLAB_MEMBER_TYPE))
  {
    return (1);
  }

  return (0);
}
