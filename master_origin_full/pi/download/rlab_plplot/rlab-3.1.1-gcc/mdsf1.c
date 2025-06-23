/* mdsf1.c Matrix Dense String Functions */

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
#include "mem.h"
#include "mds.h"
#include "mdr.h"
#include "btree.h"
#include "util.h"
#include "symbol.h"
#include "mathl.h"
#include "sort.h"
#include "rfileio.h"
#include "rlab_solver_parameters_names.h"

#include <stdio.h>
#include <math.h>

//
//
//
MDS * mds_VectorSet(MDS *t1)
{
  MDR *s=0, *s_idx=0;
  MDR *sval=0;

  int i, n1, cnt;

  if (!EQVECT(t1))
  {
    sval = mds_Create(0,0);
    return (sval);
  }

  n1 = SIZE(t1);

  // copy pointers
  s     = mds_Create(1, n1);
  s_idx = mdr_Create(1, n1);
  for (i=0; i<n1; i++)
  {
    MdsV0(s,i)     = MdsV0(t1,i);
    MdrV0(s_idx,i) = i;
  }

  // sort them
  csort ((char **) MDPTR(s), 0, n1-1, (double *)MDRPTR(s_idx));

  // count how many different ones
  cnt=1;
  for(i=1;i<n1;i++)
  {
    if (strcmp(MdsV0(s,i),MdsV0(s,i-1)))
      cnt++;
  }

  sval = mds_Create(1,cnt);
  cnt=0;
  MdsV0(sval,cnt++) = cpstr(MdsV0(s,0));
  for(i=1;i<n1;i++)
  {
    if (strcmp(MdsV0(s,i),MdsV0(s,i-1)))
      MdsV0(sval,cnt++) = cpstr(MdsV0(s,i));
  }

  // cleanup only pointers
  for (i=0;i<n1;i++)
    MdsV0(s,i) = 0;

  mds_Destroy(s);
  mdr_Destroy(s_idx);

  return (sval);
}

//
// find union of two discrete vectors
//
MDS * mds_VectorUnion(MDS *t1, MDS *t2)
{
  MDS *s=0;
  MDS *sval=0;

  int i, n1, n2;

  // special cases
  if (SIZE(t1)==0)
  {
    sval = mds_VectorSet((MDR *)t2);
    return sval;
  }
  if (SIZE(t2)==0)
  {
    sval = mds_VectorSet((MDR *)t1);
    return sval;
  }

  if ((!EQVECT(t1)) || (!EQVECT(t2)))
  {
    sval = mds_Create(0,0);
    return (sval);
  }

  n1 = SIZE(t1);
  n2 = SIZE(t2);

  // copy two together:
  s = mds_Create(1, n1+n2);
  for (i=0; i<n1; i++)
    MdsV0(s,i)    = MdsV0(t1,i);
  for (i=0; i<n2; i++)
    MdsV0(s,i+n1) = MdsV0(t2,i);

  sval = mds_VectorSet((MDS *)s);

  // cleanup only pointers
  for (i=0;i<(n1+n2);i++)
    MdsV0(s,i) = 0;
  mds_Destroy(s);

  return (sval);
}

//
//
//
MDS * mds_VectorIntersect(MDS *t1, MDS *t2)
{
  MDS *sval=0;

  if ((!EQVECT(t1)) || (!EQVECT(t2)))
  {
    sval = mds_Create(0,0);
    return (sval);
  }

  // find set of each vector:
  //  set are ascending and each element appears only once
  MDS * set_a = mds_VectorSet(t1);
  MDS * set_b = mds_VectorSet(t2);
  int ia=0, ib=0, ic=0, ic_max=MIN(SIZE(set_a),SIZE(set_b));
  MDS * set_c = mds_Create(1,ic_max);

  for (ia=0; ia<SIZE(set_a); ia++)
  {
    while (ib<SIZE(set_b))
    {
      if (strcmp(MdsV0(set_a,ia),MdsV0(set_b,ib))<=0)
        break;
      ib++;
    }

    if (ib == SIZE(set_b))
      break;

    if (!strcmp(MdsV0(set_a,ia),MdsV0(set_b,ib)))
    {
      MdsV0(set_c,ic) = MdsV0(set_a,ia);
      ic++;
    }
  }

  if (ic > 0)
  {
    sval = mds_Create(1,ic);
    for (ia=0; ia<ic; ia++)
      MdsV0(sval,ia) = cpstr(MdsV0(set_c,ia));
  }
  else
    sval = mds_Create(0,0);

  // clean up the sets
  mds_Destroy(set_a);
  mds_Destroy(set_b);
  for (ic=0;ic<ic_max;ic++)
    MdsV0(set_c,ic) = 0;
  mds_Destroy(set_c);

  return (sval);
}

//
// find complement of two discrete vectors
//
MDS * mds_VectorComplement(MDS *t1, MDS *t2)
{
  MDS *sval=0;
  int ia, ib=0, ic=0;

//   fprintf(stderr, "t1 =\n");
//   mds_Print(t1, stderr);
//   fprintf(stderr, "t2 =\n");
//   mds_Print(t2, stderr);
//   fprintf(stderr, "\n");

  // special case: if second vector is [] return the first vector
  if (EQVECT(t1))
  {
    if (SIZE(t2)<1)
    {
      sval = mds_VectorSet((MDS *)t1);
      return (sval);
    }
  }

  if ((!EQVECT(t1)) || (!EQVECT(t2)))
  {
    sval = mds_Create(0,0);
    return (sval);
  }

  // special case: if first vector is [] return []
  if (SIZE(t1)<1)
  {
    sval = mds_Create(0,0);
    return (sval);
  }

  // find set of each vector:
  //  set are ascending and each element appears only once
  MDS * set_a = mds_VectorSet(t1);
  MDS * set_b = mds_VectorIntersect(t1, t2);
  MDS * set_c = 0;

  if (SIZE(set_b) > 0)
  {
    set_c = mds_Create(1,MAX(SIZE(set_a),SIZE(set_b)));

    for (ia=0; ia<SIZE(set_a); ia++)
    {
      while (ib<SIZE(set_b))
      {
        if (strcmp(MdsV0(set_a,ia),MdsV0(set_b,ib))<=0)
          break;
        ib++;
      }

      if (ib == SIZE(set_b))
      {
        MdsV0(set_c,ic) = MdsV0(set_a,ia);
        ic++;
        continue;
      }

      if (strcmp(MdsV0(set_a,ia),MdsV0(set_b,ib)))
      {
        MdsV0(set_c,ic) = MdsV0(set_a,ia);
        ic++;
      }
    }

    if (ic > 0)
    {
      sval = mds_Create(1,ic);
      for (ia=0; ia<ic; ia++)
        MdsV0(sval,ia) = cpstr(MdsV0(set_c,ia));
    }
    else
      sval = mds_Create(0,0);
  }
  else
  {
    sval = mds_Copy(set_a);
  }

  // clean up the sets
  mds_Destroy(set_a);
  mds_Destroy(set_b);
  if (set_c)
  {
    // clean up the sets
    for (ic=0; ic<SIZE(set_c); ic++)
      MdsV0(set_c,ic) = 0;
    mds_Destroy(set_c);
  }

  return (sval);
}

/* **************************************************************
 * Sum a string matrix:
 *  if vector, just add terms
 *  if matrix, add rows (unlike numbers, where one adds columns)
 * ************************************************************** */

MDS *
mds_Sum (MDS * m, void *cs)
{
  int i,j,k,slen,cslen=0;
  char *s=0, *p, *colsep;
  MDS *sval=0;

  colsep = (char *) cs;
  if (colsep)
  {
    cslen=strlen(colsep);
  }

  if ( EQNULL(m) )
  {
    sval = mds_Create(0,0);
  }
  else if ( EQVECT(m) )
  {
    slen = 0;
    for (i=0; i<SIZE(m); i++)
    {
      slen += strlen(MdsV0(m,i));
    }
    if (cslen)
      slen += (SIZE(m)-1) * cslen;
    if (slen)
    {
      s = GC_MALLOC((slen+1)*sizeof(char));
      p = s;
      for (i=0; i<SIZE(m); i++)
      {
        // copy string
        if(strlen(MdsV0(m,i)))
        {
          for (j=0; j<strlen(MdsV0(m,i)); j++)
          {
            *p = MdsV0(m,i)[j];
            p++;
          }
        }
        // copy column separator for all entries but the last
        if (i<(SIZE(m)-1) && cslen)
        {
          for (j=0; j<cslen; j++)
          {
            *p = colsep[j];
            p++;
          }
        }
        *p='\0';
      }
      sval = mds_Create(1,1);
      MdsV0(sval,0) = s;
    }
    else
      sval = mds_CreateScalar ("");
  }
  else
  {
    sval = mds_Create(MNR(m),1);
    for (i=0; i<MNR(m); i++)
    {
      slen = 0;
      for (j=0; j<MNC(m); j++)
      {
        slen += strlen(Mds0(m,i,j));
      }
      if (cslen)
        slen += (MNC(m)-1) * cslen;
      if (slen)
      {
        s = GC_MALLOC((slen+1)*sizeof(char));
        p = s;
        for (j=0; j<MNC(m); j++)
        {
          if(strlen(Mds0(m,i,j)))
          {
            for (k=0; k<strlen(Mds0(m,i,j)); k++)
            {
              *p = Mds0(m,i,j)[k];
              p++;
            }
          }
          if (j<(MNC(m)-1) && cslen)
          {
            for (k=0; k<cslen; k++)
            {
              *p = colsep[k];
              p++;
            }
          }
          *p='\0';
        }
        MdsV0(sval,i) = s;
      }
      else
        MdsV0(sval,i) = cpstr("");
    }
  }

  return sval;
}


/* **************************************************************
 * Return matrix-dense-string member references.
 * ************************************************************** */

void *
mds_MemberRef (MDS * m, char *name, int *type)
{
  Ent *ne;
  void *rptr;

  if (!strcmp (name, RLAB_MEMBER_NROW))
  {
    ne = ent_Create ();
    ent_data (ne) = mdr_CreateScalar ((double) MNR (m));
    ent_SetType (ne, MATRIX_DENSE_REAL);
    rptr = (void *) ne;
    *type = ENTITY;
  }
  else if (!strcmp (name, RLAB_MEMBER_NCOL))
  {
    ne = ent_Create ();
    ent_data (ne) = mdr_CreateScalar ((double) MNC (m));
    ent_SetType (ne, MATRIX_DENSE_REAL);
    rptr = (void *) ne;
    *type = ENTITY;
  }
  else if (!strcmp (name, RLAB_MEMBER_SIZE))
  {
    ne = ent_Create ();
    ent_data (ne) = mdr_CreateScalar ((double) (MNR (m) * MNC (m)));
    ent_SetType (ne, MATRIX_DENSE_REAL);
    rptr = (void *) ne;
    *type = ENTITY;
  }
  else if (!strcmp (name, RLAB_MEMBER_CLASS))
  {
    ne = ent_Create ();
    ent_data (ne) = mds_CreateScalar ( RLAB_CLASS_STRING );
    ent_SetType (ne, MATRIX_DENSE_STRING);
    rptr = (void *) ne;
    *type = ENTITY;
  }
  else if (!strcmp (name, RLAB_MEMBER_TYPE))
  {
    ne = ent_Create ();
    ent_data (ne) = mds_CreateScalar ( RLAB_CLASS_STRING );
    ent_SetType (ne, MATRIX_DENSE_STRING);
    rptr = (void *) ne;
    *type = ENTITY;
  }
  else if (!strcmp (name, RLAB_MEMBER_STORAGE))
  {
    ne = ent_Create ();
    ent_data (ne) = mds_CreateScalar ( RLAB_STORAGE_DENSE );
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
mds_Members (MDS * m, int *n)
{
  char **marray;

  marray = (char **) GC_MALLOC (6 * sizeof (char *));
  if (marray == 0)
    rerror ("out of memory");

  marray[0] = cpstr (RLAB_MEMBER_NROW);
  marray[1] = cpstr (RLAB_MEMBER_NCOL);
  marray[2] = cpstr (RLAB_MEMBER_SIZE);
  marray[3] = cpstr (RLAB_MEMBER_CLASS);
  marray[4] = cpstr (RLAB_MEMBER_TYPE);
  marray[5] = cpstr (RLAB_MEMBER_STORAGE);

  *n = 6;
  return (marray);
}

/* **************************************************************
 * Return the "type" of a string matrix for type() builtin.
 * ************************************************************** */

MDS *
mds_Type_BF (MDS * m)
{
  MDS *type = mds_CreateScalar ( RLAB_CLASS_STRING );
  return (type);
}

/* **************************************************************
 * Test a string matrix for symmetry.
 * ************************************************************** */

int
mds_IsSymmetric (MDS * m)
{
  int i, j, nr, nc;

  nr = MNR (m);
  nc = MNC (m);

  for (j = 0; j < nc; j++)
  {
    for (i = j + 1; i < nr; i++)
    {
      if (strcmp (Mds0 (m, i, j), Mds0 (m, j, i)))
	return (0);
    }
  }
  return (1);
}

/* **************************************************************
 * Find the indices of non-zero elements in a matrix.
 * ************************************************************** */

MDR *
mds_Find_BF (MDS * m)
{
  int i, j, size;
  double *dtmp;
  MDR *mf;

  /* Make mfind as big as it could be, reduce later */
  mf = mdr_Create (1, MNR (m) * MNC (m));

  j = 0;
  size = MNR (m) * MNC (m);
  for (i = 0; i < size; i++)
  {
    if ((MdsV0 (m, i) != 0) && strncmp (MdsV0 (m, i), "\0", 1))
      MdrV0 (mf, j++) = i + 1;
  }

  /* Now reduce mfind to correct size */
  if (j == 0)
  {
    /* No reduction here, just make dimensions correct */
    mf->nrow = 0;
    mf->ncol = 0;
  }
  else if (j < size)
  {
    dtmp = (double *) GC_malloc_atomic_ignore_off_page (sizeof (double) * j);

    memcpy (dtmp, MDPTR(mf), sizeof (double) * j);
    GC_FREE (MDPTR(mf));
    MDPTR(mf) = (void *) dtmp;
    mf->ncol = j;
  }

  return (mf);
}

/* **************************************************************
 * Sort a vector, or matrix (builtin-function).
 * ************************************************************** */

extern MDR *mdr_CreateFillSind (int nrow, int ncol);

Btree *
mds_Sort_BF (MDS * m)
{
  int i, n;
  Btree *bt=0;
  MDR *sind=0;
  MDS *mcopy=0;

  if (MNR (m) == 1 || MNC (m) == 1)
  {
    /* Vector sort */
    n = MAX(MNR(m), MNC(m));
    sind = mdr_CreateFillSind (n, 1);
    sind->ncol = sind->nrow;
    sind->nrow = 1;
    mcopy = mds_Copy (m);
    csort ((char **) MDSPTR (mcopy), 0, n - 1, (double *) MDRPTR (sind));
  }
  else
  {
    /* Matrix sort (column-wise) */
    n = MNR (m);
    sind = mdr_CreateFillSind (MNR (m), MNC (m));
    mcopy = mds_Copy (m);
    for (i = 0; i < MNC (m); i++)
    {
      csort (&MdsV0(mcopy,i*n), 0, n - 1, &MdrV0(sind,i*n));
    }
  }

  bt = btree_Create ();
  install (bt, RLAB_NAME_GEN_INDEX, ent_Assign_Rlab_MDR(sind));
  install (bt, RLAB_NAME_GEN_VALUE, ent_Assign_Rlab_MDS(mcopy));
  return (bt);
}

/* **************************************************************
 * Strtod function
 * ************************************************************** */

static char *endptr;

MDR * mds_Strtod_BF (MDS * m1)
{
  char *oneline;
  char *myendptr;
  char *refstr;
  int i, j, k, nel, ncol, nrow, lenk, imy, nnz;
  MDS *m = mds_Copy (m1);
  MDR *new;

  nrow = MNR (m);
  ncol = MNC (m);

  if (MNC (m) > 1)
  {
    //
    // if more then one column, assume each column carries one entry
    //
    new = mdr_Create (nrow, ncol);
    nel = nrow * ncol;
    for (i = 0; i < nel; i++)
    {
      MdrV0 (new, i) = strtod (MdsV0 (m, i), NULL);
      if (endptr == MdsV0 (m, i) && MdrV0 (new, i) == 0.0)
      {
        /* strtod could not recognize the string as a number */
        MdrV0 (new, i) = create_nan ();
      }
    }
  }
  else if (MNC (m) == 1)
  {
    //
    // find first non-empty string
    //
    i = 0;
    while (1)
    {
      if (i == MNR (m) - 1)
        break;
      // is it an empty string
      if (strlen (MdsV0 (m, i)) == 0)
      {
        i++;
        continue;
      }
      // eliminate special characters and spaces at the end of the string
      imy = strlen (MdsV0 (m, i)) - 1;
      if (imy > 0)
      {
        for (j = imy; j > 0; j--)
        {
          if (MdsV0 (m, i)[j] < 42 || MdsV0 (m, i)[j] > 57)
            MdsV0 (m, i)[j] = '\0';
          else
            break;
        }
      }
      if (strlen (MdsV0 (m, i)) == 0)
      {
        i++;
        continue;
      }
      break;
    }
    if (strlen (MdsV0 (m, i)) == 0)
      return mdr_Create (0, 0);

    //
    // i-th entry contains non-empty string
    //
    ncol = 0;
    refstr = cpstr (MdsV0 (m, i));
    myendptr = refstr;
    oneline = 0;
    while (myendptr != oneline)
    {
      oneline = myendptr;
      strtod (oneline, &myendptr);
      if (myendptr[0] == ',' || myendptr[0] == ':' ||myendptr[0] == '\t')
        myendptr[0] = ' ';
      ncol++;
    }
    GC_FREE(refstr);
    ncol--;

    if (ncol < 1)
    {
      new = mdr_Create(0,0);
      return (new);
    }


    //
    // find number of non-empty rows
    //
    nnz = 0;
    for (j=i; j<nrow; j++)
    {
      if (isvalidstring(MdsV0(m,j)) < 1)
        continue;
      nnz++;
    }
    new = mdr_Create (nnz, ncol);

    //
    // do the strtod for non-empty rows only
    //
    int i1 = 0;
    for ( ; i < nrow; i++)
    {
      if (isvalidstring(MdsV0(m,i)) < 1)
        continue;
      refstr = cpstr (MdsV0 (m, i));
      oneline = refstr;
      lenk = strlen (oneline);
      // non(digits and dec.point) are replaced by spaces
      for (k = 0; k < lenk; k++)
        if (oneline[k] < 42 || oneline[k] > 122)
          oneline[k] = ' ';
      // chomp spaces at the end
      chomp_string(oneline);
//       for (k = lenk - 1; k > 0; k--)
//       {
//         if (oneline[k] == ' ')
//           oneline[k] = '\0';
//         else
//           break;
//       }
      myendptr = oneline;
      j = 0;
      while (strlen (oneline) > 0 && j < ncol)
      {
        Mdr0 (new, i1, j) = strtod (oneline, &myendptr);
        if (myendptr[0] == ',')
          myendptr[0] = ' ';
        else if (myendptr[0] == '\t')
          myendptr[0] = ' ';
        oneline = myendptr;
        j++;
      }
      i1++;
      GC_FREE(refstr);
    }
  }
  else
    new = mdr_Create(0,0);

  return (new);
}

//MDR *
//mds_Strtol_BF (MDS * m, int base)
//{
//  int i, nel;
//  MDR *new;

//  new = mdr_Create (MNR (m), MNC (m));
//  nel = MNR (m) * MNC (m);

//  for (i = 0; i < nel; i++)
//  {
//    endptr = 0;
//    MdrV0 (new, i) = (double) strtol (MdsV0 (m, i), &endptr, base);
//    if (endptr == MdsV0 (m, i) && MdrV0 (new, i) == 0.0)
//    {
//      /* strtol could not recognize the string */
//      MdrV0 (new, i) = create_nan ();
//    }
//  }
//  return (new);
//}

MDR *
mds_Strtol_BF (MDS * m1, int base)
{
  char *oneline;
  char *myendptr;
  int i, j, k, nel, ncol, nrow, lenk, imy;
  MDS *m = mds_Copy (m1);
  MDR *new;

  nrow = MNR (m);
  ncol = MNC (m);

  if (MNC (m) > 1)
  {
    //
    // if more then one column, assume each column carries one entry
    //
    new = mdi_Create (nrow, ncol);
    nel = nrow * ncol;
    for (i = 0; i < nel; i++)
    {
      MdiV0 (new, i) = strtol (MdsV0 (m, i), NULL, base);
      if (endptr == MdsV0 (m, i) && MdrV0 (new, i) == 0.0)
      {
        /* strtod could not recognize the string as a number */
        MdiV0 (new, i) = create_nan ();
      }
    }
  }
  else
  {
    //
    // if one column figure out how many entries there are
    //
    i = 0;
    while (1)
    {
      if (i == MNR (m) - 1)
        break;
      // is it an empty string
      if (strlen (MdsV0 (m, i)) == 0)
      {
        i++;
        continue;
      }
      // eliminate special characters and spaces at the end of the string
      imy = strlen (MdsV0 (m, i)) - 1;
      for (j = imy; j > 0; j--)
      {
        if (MdsV0 (m, i)[j] < 42 || MdsV0 (m, i)[j] > 57)
          MdsV0 (m, i)[j] = '\0';
        else
          break;
      }
      if (strlen (MdsV0 (m, i)) == 0)
      {
        i++;
        continue;
      }
      break;
    }
    if (strlen (MdsV0 (m, i)) == 0)
      return mdr_Create (0, 0);
    ncol = 0;
    myendptr = cpstr (MdsV0 (m, i));
    oneline = 0;
    while (myendptr != oneline)
    {
      oneline = myendptr;
      strtod (oneline, &myendptr);
      if (myendptr[0] == ',' || myendptr[0] == '\t')
        myendptr[0] = ' ';
      ncol++;
    }
    ncol--;
    new = mdi_Create (nrow, ncol);
    for (i = 0; i < nrow; i++)
    {
      oneline = cpstr (MdsV0 (m, i));
      lenk = strlen (oneline);
      // non(digits and dec.point) are replaced by spaces
      for (k = 0; k < lenk; k++)
        if (oneline[k] < 42 || oneline[k] > 122)
          oneline[k] = ' ';
      // chomp spaces at the end
      chomp_string(oneline);
//       for (k = lenk - 1; k > 0; k--)
//       {
//         if (oneline[k] == ' ')
//           oneline[k] = '\0';
//         else
//           break;
//       }
      myendptr = oneline;
      j = 0;
      while (strlen (oneline) > 0 && j < ncol)
      {
        Mdi0 (new, i, j) = strtol (oneline, &myendptr,base);
        if (myendptr[0] == ',')
          myendptr[0] = ' ';
        else if (myendptr[0] == '\t')
          myendptr[0] = ' ';
        oneline = myendptr;
        j++;
      }
    }
  }

  return (new);
}

size_t
mds_Sizeof (MDS * m)
{
  unsigned int i;
  size_t msize, size;

  msize = 0;
  size = (size_t) (MNR (m) * MNC (m));
  for (i = 0; i < size; i++)
  {
    if (MdsV0 (m, i) != 0)
      msize += (size_t) strlen (MdsV0 (m, i));
  }

  return (size_t) ((sizeof (char) * msize) + (sizeof (char *) * size));
}


