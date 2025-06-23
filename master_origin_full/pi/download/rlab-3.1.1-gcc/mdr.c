/* mdr.c Matrix Dense Real */

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
#include "rfileio.h"
#include "mem.h"
#include "mdr.h"
#include "mdrf2.h"
#include "mdc.h"
#include "mds.h"
#include "util.h"
#include "bltin1.h"
#include "mathl.h"

#include "fi.h"
#include "blas.h"
#include "lp.h"

#include <stdio.h>
#include <math.h>

// #include <termcap.h>
#include <termios.h>
#include <sys/fcntl.h>
#include <unistd.h>
#include <sys/ioctl.h>

#include <gsl/gsl_math.h>

#ifdef __riscos
#include <limits.h>
#define MAXINT INT_MAX
#endif

#include "rlab_solver_parameters_names.h"
#include "rlab_macros.h"
#include "rlab_macros_code.h"

#undef  THIS_FILE
#define THIS_FILE "mdr.c"

// rfileio.c:
extern char *RLABPLUS_EOL;

int ione=1;

int mdr_is_positive (MDR * m);

int md_Isvalid(MD *m)
{
  int rval = 0;

  if (m)
    if (SIZE(m)>0)
      rval = 1;

  return rval;
}

// find max of a matrix, ignore nan's
double rlab_dmax_mdr(MDR *x, int ignore_inf)
{
  double rval=create_nan();
  int length = MNR(x) * MNC(x);

  int i,j;

  if (length >0)
  {
    // find first not-nan
    j=0;
    while (j<length)
    {
      rval=mdrV0(x,j);
      j++;

      if(!isnand (rval))
        break;

      if (ignore_inf)
        if (!isinf(rval))
          break;

    }

    // now do the comparison
    for (i=j; i<length; i++)
    {
      if(isnand (mdrV0(x,i)))
        continue;

      if (ignore_inf)
        if (isinf(mdrV0(x,i)))
          continue;

      if(rval < mdrV0(x,i))
        rval = mdrV0(x,i);
    }
  }

  return rval;
}

// find min of a matrix, ignore nan's
double rlab_dmin_mdr(MDR *x, int ignore_inf)
{
  double rval=create_nan();
  int length = SIZE(x);

  int i,j;

  if (length > 0)
  {
    // find first not-nan (and non-inf)
    j=0;
    while (j<length)
    {
      rval = mdrV0(x,j);
      j++;

      if(!isnand (rval))
        break;

      if (ignore_inf)
        if (!isinf(rval))
          break;
    }

    // now do the comparison
    for (i=j; i<length; i++)
    {
      if(isnand (mdrV0(x,i)))
        continue;

      if (ignore_inf)
        if (isinf(mdrV0(x,i)))
          continue;

      if(rval > mdrV0(x,i))
        rval = mdrV0(x,i);
    }
  }

  return rval;
}

// find a mean of vector - ignore nan's, inf's
double rlab_mean_vector(int length, double *x, MDR *ig_idx, MDR * iuse_idx, int ig_infs)
{
  double rval = create_nan();
  int i, k;
  double count=0;
  MDR * r_idx=0;

  if (length==0)
    return rval;

  if (!ig_idx && !iuse_idx)
  {
    for (i=0; i<length; i++)
    {
      if(isnand (x[i]))
        continue;

      count ++;

      if (isnand(rval))
        rval = x[i];
      else
        rval += x[i];
    }

    if (count)
      rval /= count;

    return rval;
  }

  // book keeping array so that any index can be removed at most once
  // if provided then use only those indices
  if (iuse_idx)
    r_idx = iuse_idx;
  else
  {
    r_idx = mdi_Create(1,length);
    for(i=0; i<length;i++)
      MdiV0(r_idx,i) = 1;
  }

  if (ig_idx)
  {
    for (k=0; k<MNR(ig_idx)*MNC(ig_idx); k++)
    {
      if ((mdiV0(ig_idx,k) < 0) || (mdiV0(ig_idx,k) > length-1 ))
        continue;
      MdiV0(r_idx,mdiV0(ig_idx,k))=0;
    }
  }

  for (i=0; i<length; i++)
  {
    if(isnand (x[i]))
      continue;

    if( (x[i]==create_inf()) || (-x[i]==create_inf()) )
    {
      if (ig_infs)
        continue;
      else
      {
        rval = create_inf();
        break;
      }
    }

    if (mdiV0(r_idx,i)==0)
      continue;

    count ++;

    if (isnand(rval))
      rval = x[i];
    else
      rval += x[i];
  }

  if (!iuse_idx)
    mdr_Destroy(r_idx);

  if (count)
    rval /= count;
  return rval;
}

// find a mean of vector - ignore nan's
void rlab_meanstat_from_valstat_vectors(int length, double *x, double *w, MDR *ig_idx, MDR * iuse_idx,
                                        double *m1, double *v1, int do_std)
{
  int i, k;
  double sw=0, sw2=0, dw, m=create_nan(), v=create_nan();
  MDR * r_idx=0;

  if (length==0)
    return;

  // book keeping array so that any index can be removed at most once
  // if provided then use only those indices
  if (iuse_idx)
  {
    r_idx = mdi_Create(1,length);
    mdr_Zero(r_idx);
    for (i=0; i<SIZE(iuse_idx); i++)
    {
      if (mdiV0(iuse_idx,i)<0 || mdiV0(iuse_idx,i)>length-1)
        continue;
      MdiV0(r_idx,mdiV0(iuse_idx,i))=1;
    }
  }
  else
  {
    r_idx = mdi_Create(1,length);
    for(i=0; i<length;i++)
      MdiV0(r_idx,i) = 1;
  }

  // check 'x' for nans
  for(i=0; i<length;i++)
  {
    if(isnand(x[i]) || isnand(w[i]))
      MdiV0(r_idx,i)=0;
  }

  if (ig_idx)
  {
    for (k=0; k<MNR(ig_idx)*MNC(ig_idx); k++)
    {
      if ((MdiV0(ig_idx,k) < 0) || (MdiV0(ig_idx,k) > length-1 ))
        continue;
      MdiV0(r_idx,mdiV0(ig_idx,k))=0;
    }
  }

  for (i=0; i<length; i++)
  {
    if (!MdiV0(r_idx,i))
      continue;

    if (do_std)
      dw = 1/(w[i] * w[i]);
    else
      dw = w[i];


      // update mean
    if (isnand(m))
    {
      m   = dw * x[i];
      sw  = dw;
      sw2 = dw * dw;
    }
    else
    {
      m   += dw * x[i];
      sw  += dw;
      sw2 += dw * dw;
    }
      // update variance
    if (isnand(v))
      v  = dw * x[i] * x[i];
    else
      v += dw * x[i] * x[i];
  }

  mdr_Destroy(r_idx);

  m /= sw;
  v  = (v/sw - m*m) * (sw*sw)/((sw*sw) - sw2);

  if (do_std)
  {
    if (v > 0)
      v = sqrt(v);
    else
      v = create_nan();
  }
  else
    v = 1/v;

  *m1 = m;
  *v1 = v;

  return;
}

// find a mean of vector - ignore nan's
double rlab_var_vector(int length, double *x, int bias, MDR *ig_idx, MDR * iuse_idx)
{
  double rval = create_nan(), rval2=create_nan();
  int i, k;
  double count=0;
  MDR * r_idx=0;

  if (length==0)
    return rval;

  // book keeping array so that any index can be removed at most once
  // if provided then use only those indices
  if (iuse_idx)
    r_idx = iuse_idx;
  else
  {
    r_idx = mdi_Create(1,length);
    for(i=0; i<length;i++)
      MdiV0(r_idx,i) = 1;
  }

  if (ig_idx)
  {
    for (k=0; k<MNR(ig_idx)*MNC(ig_idx); k++)
    {
      if ((mdiV0(ig_idx,k) < 0) || (mdiV0(ig_idx,k) > length-1 ))
        continue;
      MdiV0(r_idx,mdiV0(ig_idx,k))=0;
    }
  }

  for (i=0; i<length; i++)
  {
    if(isnand (x[i]))
      continue;

    if (mdiV0(r_idx,i)==0)
      continue;

    count += 1.0;
    if (isnand(rval))
    {
      rval  = x[i];
      rval2 = x[i] * x[i];
    }
    else
    {
      rval  += x[i];
      rval2 += x[i] * x[i];
    }
  }

  if (!iuse_idx)
    mdr_Destroy(r_idx);

  if (count)
    rval = rval2 / (count-bias) - (rval*rval) / (count-bias) / count;

  return rval;
}



/* **************************************************************
 * Create a matrix.
 * ************************************************************** */
int rsizeof (Rtype rt)
{
  int rval=0;

  switch(rt)
  {
    case RLAB_TYPE_INT32:
      rval = sizeof(int);
      break;

    case RLAB_TYPE_DOUBLE:
      rval = sizeof(double);
      break;

    case RLAB_TYPE_COMPLEX:
      rval = 2*sizeof(double);
      break;

    case RLAB_TYPE_STRING:
      rval = sizeof(void *);
      break;

    default:
      rval = 0;
      break;
  }

  return rval;
}

void md_numeric_Flip (MD *m, int ud, int lr)
{
  if ((!ud) && (!lr))
    return;

  int nr = MNR(m);
  int nc = MNC(m);
  int nd = rsizeof(m->type);

  if ((nr<1)||(nc<1))
    return;

  int i, j, k;
  unsigned char * data = (unsigned char *) MDPTR(m);
  unsigned char tmp;

  if (ud)
  {
    for (j=0; j<nc; j++)
    {
      for (i=0; i<(nr/2); i++)
        for (k=0; k<nd; k++)
        {
          tmp = data[(nr-i-1)*nd + nr*nd*j + k];
          data[(nr-i-1)*nd + nr*nd*j + k] = data[i*nd + nr*nd*j + k];
          data[i*nd + nr*nd*j + k] = tmp;
        }
    }
  }

  if (lr)
  {
    for (i=0; i<nr; i++)
    {
      for (j=0; j<(nc/2); j++)
        for (k=0; k<nd; k++)
        {
          tmp = data[i*nd + nr*nd*(nc-j-1) + k];
          data[i*nd + nr*nd*(nc-j-1) + k] = data[i*nd + nr*nd*j + k];
          data[i*nd + nr*nd*j + k] = tmp;
        }
    }
  }

  return;
}
void md_numeric_Shift (MD *m, int ud, int lr)
{
  if ((!ud) && (!lr))
    return;

  int nr = MNR(m);
  int nc = MNC(m);
  int nd = rsizeof(m->type);

  if ((nr<1)||(nc<1))
    return;

  int i, j, k;
  unsigned char * data = (unsigned char *) MDPTR(m);

  if (ud != 0)
  {
    // shift up: ud>0
    // or down: ud<0
    if ((ud >= nr)||(ud <= -nr))
    {
      for (i=0; i<nr; i++)
        for (j=0; j<nc; j++)
          for (k=0; k<nd; k++)
            data[i*nd + nr*nd*j + k] = 0x00;
      return;
    }

    if (ud > 0)
    {
      for (j=0; j<nc; j++)
      {
        for (i=ud; i<nr; i++)
          for (k=0; k<nd; k++)
            data[(i-ud)*nd + nr*nd*j + k] = data[i*nd + nr*nd*j + k];
        for (i=0; i<ud; i++)
          for (k=0; k<nd; k++)
            data[(nr-i-1)*nd + nr*nd*j + k] = 0x00;
      }
    }
    if (ud < 0)
    {
      for (j=0; j<nc; j++)
      {
        for (i=nr+ud-1; i>=0; i--)
          for (k=0; k<nd; k++)
            data[(i-ud)*nd + nr*nd*j + k] = data[i*nd + nr*nd*j + k];
        for (i=0; i<-ud; i++)
          for (k=0; k<nd; k++)
            data[i*nd + nr*nd*j + k] = 0x00;
      }
    }
  }

  if (lr != 0)
  {
    // shift left: lr>0
    // or right: lr<0
    if ((lr >= nc)||(lr <= -nc))
    {
      for (i=0; i<nr; i++)
        for (j=0; j<nc; j++)
          for (k=0; k<nd; k++)
            data[i*nd + nr*nd*j + k] = 0x00;
      return;
    }

    if (lr > 0)
    {
      for (i=0; i<nr; i++)
      {
        for (j=lr; j<nc; j++)
          for (k=0; k<nd; k++)
            data[i*nd + nr*nd*(j-lr) + k] = data[i*nd + nr*nd*j + k];
        for (j=0; j<lr; j++)
          for (k=0; k<nd; k++)
            data[i*nd + nr*nd*(nc-j-1)+ k] = 0x00;
      }
    }
    if (lr < 0)
    {
      for (i=0; i<nr; i++)
      {
        for (j=nc+lr-1; j>=0; j--)
          for (k=0; k<nd; k++)
            data[i*nd + nr*nd*(j-lr) + k] = data[i*nd + nr*nd*j + k];
        for (j=0; j<-lr; j++)
          for (k=0; k<nd; k++)
            data[i*nd + nr*nd*j + k] = 0x00;
      }
    }
  }

  return;
}


MD * md_Create (int nrow, int ncol, Rtype rt)
{
  if (nrow < 0 || ncol < 0)
  {
    fprintf (stderr, "illegal dimensions for matrix: ");
    fprintf (stderr, "nrow=%i, ncol=%i\n", nrow, ncol);
    rerror ("cannot specify a negative matrix dimension");
  }

  MDR *new = (MDR *) GC_MALLOC (sizeof (MDR));
  if (new == 0)
    rerror ("out of memory");

  // default initialization
  new->nrow = nrow * (ncol > 0);
  new->ncol = ncol * (nrow > 0);
  new->list = 0;
  MDPTR(new) = NULL;
  new->type = rt;

  // we allow for objects of size 0x0
  if (nrow==0 || ncol==0)
    return (new);

// fprintf(stderr, "nrow=%g,ncol=%g,MAXINT=%i\n", nrow, ncol, MAXINT);

  double dsize = (double) ((double) nrow * (double) ncol);
  if (dsize >= (double) MAXINT)
  {
    fprintf (stderr,
             "ERROR: requested matrix size caused integer overflow\n");
    rerror ("out of memory");
  }

  if (rsizeof(rt))
    MDPTR(new) = (void *) GC_malloc_atomic_ignore_off_page(dsize * rsizeof(rt));

  return (new);
}

MDR * mdr_Create (int nrow, int ncol)
{
  MDR *new=0;

  if ((nrow>0)&&(ncol>0))
    new = md_Create (nrow, ncol, RLAB_TYPE_DOUBLE);
  else
    new = md_Create (0, 0, RLAB_TYPE_DOUBLE);

  new->allocated = 1;
  new->row_dominant = 0;

  return (new);
}

MDR * mdr_CreateEmpty (int nrow, int ncol)
{
  MDR *new=0;

  new = md_Create (0, 0, RLAB_TYPE_DOUBLE);
  if ((nrow>0)&&(ncol>0))
  {
    MNR(new) = nrow;
    MNC(new) = ncol;
  }

  new->allocated = 0;
  new->row_dominant = 0;
  return (new);
}


MDR * mdr_Create_SameSize (MDR *x)
{
  int nrow = 0;
  int ncol = 0;
  MDR *new=0;

  if (SIZE(x)>0)
  {
    nrow = MNR(x);
    ncol = MNC(x);
  }

  new = md_Create (nrow, ncol, RLAB_TYPE_DOUBLE);

  return (new);
}


MDR * mdi_Create (int nrow, int ncol)
{
  MDR *new=0;

  if ((nrow>0)&&(ncol>0))
  {
    new = md_Create (nrow, ncol, RLAB_TYPE_INT32);
  }
  else
  {
    new = md_Create (0, 0, RLAB_TYPE_INT32);
  }

  return (new);
}

MDR * mdi_Create_SameSize (MDR *x)
{
  int nrow = MNR(x);
  int ncol = MNC(x);
  MDR *new=0;

  if (EQNULL(x))
    new = md_Create (0, 0, RLAB_TYPE_INT32);
  else
    new = md_Create (nrow, ncol, RLAB_TYPE_INT32);

  return (new);
}


/*
 * Create a "scalar" (1-by-1) REAL matrix.
 */

MDR * mdr_CreateScalar (double val)
{
  MDR *new;

  new = mdr_Create (1, 1);
  MdrV0 (new, 0) = val;
  return (new);
}

MDR * mdi_CreateScalar (int ival)
{
  MDR *new;

  new = mdi_Create (1, 1);
  MdiV0 (new, 0) = ival;
  return (new);
}

/* **************************************************************
 * Free a matrix, and wipe out the structure members.
 * ************************************************************** */
void * md_Destroy (MD * m)
{
  // something to destroy?
  if (!m)
    return 0;

  // proceed with destruction
  m->nrow = -1;
  m->ncol = -1;

  // does it own the data to which it points?
  if (m->allocated)
  {
    if (MDPTR(m))
    {
      GC_FREE (MDPTR(m));
      MDPTR(m) = 0;
    }
  }

  if (m->list)
  {
    btree_Destroy (m->list);
    m->list = 0;
  }

  GC_FREE (m);

  return NULL;
}

void * mdr_Destroy (MDR * m)
{
  return md_Destroy(m);
}

/* **************************************************************
 * Copy a matrix. Create the new matrix, and return the new,
 * copied matrix as the result.
 * ************************************************************** */
unsigned char * data_Copy (unsigned char *data, int nrow, int ncol, int size_a)
{
  unsigned char * new=0;

  if (nrow*ncol == 0)
    return (new);

  new = GC_MALLOC (nrow * ncol * size_a);
  memcpy (new, data, nrow * ncol * size_a);

  return (new);
}

MDR * mdr_Copy (MDR * m)
{
  MDR *new = GC_MALLOC(sizeof(MDR));

  // copy everything first
  memcpy (new, m, sizeof (MDR));

  MDPTR(new) = (void *) data_Copy ((unsigned char *)MDPTR(m),
        new->nrow, new->ncol, rsizeof(m->type));

  // Copy the list if there is one
  if (m->list)
  {
    new->list = btree_Copy (m->list);
  }

  return (new);
}

void mdr_Copy1IntoPtr2 (MDR * m, MDR **new)
{
  if (! (*new))
  {
    *new = mdr_Copy(m);
    return;
  }

  if (m->type == (*new)->type)
  {
    if ((MNR(m)==MNR(*new)) && (MNC(m)==MNC(*new)))
    {
      memcpy (MDPTR(*new), MDPTR(m), rsizeof(m->type) * MNR(m) * MNC(m));
    }
  }


  return;
}



/* **************************************************************
 * Reshape a matrix (change nrow, ncol).
 * ************************************************************** */
MDR * mdr_Reshape (MDR * m, int nrow, int ncol)
{
  MDR *new;

  if (nrow * ncol != MNR (m) * MNC (m))
  {
    fprintf (stderr, "incompatible dimensions for reshape\n");
    fprintf (stderr, "nrow*ncol must equal: %i\n", MNR (m) * MNC (m));
    rerror ("error");
  }

  new = mdr_Copy (m);

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
 * Extend a matrix. Realloc the data, extending the data size.
 * Zero out the newly added space.
 * ************************************************************** */
unsigned char *
    md_extend (unsigned char *data_old, int nrow_old, int ncol_old,
               int nrow_new, int ncol_new, int size_a)
{
  int i, j;
  unsigned char * new=0;

  // refuse to shrink to empty matrix
  if (nrow_new*ncol_new == 0)
    return (new);

  new = (unsigned char *) GC_malloc_atomic_ignore_off_page (nrow_new * ncol_new * size_a);

  if (!new)
    rerror("terrible internal error");

  // copy 'data_old' to 'new'
  if (nrow_old * ncol_old > 0)
  {
    for (i=0; i < MIN(nrow_old, nrow_new); i++)
      for (j=0; j < MIN(ncol_old, ncol_new); j++)
    {
      memcpy (&new[(j * nrow_new + i)*size_a],
              &data_old[(j * nrow_old + i)*size_a], size_a);
    }
  }

  // fill the rest with zeros
  if (nrow_old < nrow_new)
  {
    for (i=nrow_old; i < nrow_new; i++)
      for (j=0; j < MIN(ncol_old, ncol_new); j++)
    {
      memset (&new[(j*nrow_new + i)*size_a], 0, size_a);
    }
  }

  if (ncol_old < ncol_new)
  {
    for (i=0; i < nrow_new; i++)
      for (j=ncol_old; j < ncol_new; j++)
    {
      memset (&new[(j*nrow_new + i)*size_a], 0, size_a);
    }
  }

  return (new);
}

void mdr_Extend (MDR * m, int nrow, int ncol)
{
  unsigned char *new;
  int size_a;
  int nr_old = m->nrow;
  int nc_old = m->ncol;

  size_a   = rsizeof(m->type);
  new = md_extend ((unsigned char *) MDPTR(m), nr_old, nc_old, nrow, ncol, size_a);

  m->nrow = nrow;
  m->ncol = ncol;
  GC_FREE (MDPTR(m));
  MDPTR(m) = (void *) new;

  return;
}

/* **************************************************************
 * Zero the elements of an existing matrix.
 * ************************************************************** */
void mdr_Zero (MDR * m)
{
  int k, isize;
  isize = MNR (m) * MNC (m);

  switch (m->type)
  {
    case RLAB_TYPE_INT32:
      for (k = 0; k < isize; k++)
        MdiV0(m,k) = 0;
      break;

    case RLAB_TYPE_DOUBLE:
      for (k = 0; k < isize; k++)
        MdrV0(m,k) = 0.0;
      break;

    default:
      break;
  }
}

/* **************************************************************
 * nan-ify the elements of an existing matrix.
 * ************************************************************** */
void mdr_Nan (MDR * m)
{
  int k, isize;
  isize = MNR (m) * MNC (m);

  switch (m->type)
  {
    case RLAB_TYPE_INT32:
      for (k = 0; k < isize; k++)
        MdiV0(m,k) = 0;
      break;

    case RLAB_TYPE_DOUBLE:
      for (k = 0; k < isize; k++)
        MdrV0(m,k) = create_nan ();
      break;

    default:
      break;
  }
}

/* **************************************************************
 * Ones the elements of an existing matrix.
 * ************************************************************** */
void mdr_Ones (MDR * m)
{
  int k, isize;
  isize = MNR (m) * MNC (m);

  switch (m->type)
  {
    case RLAB_TYPE_INT32:
      for (k = 0; k < isize; k++)
        MdiV0(m,k) = 1;
      break;

    case RLAB_TYPE_DOUBLE:
      for (k = 0; k < isize; k++)
        MdrV0(m,k) = 1.0;
      break;

    default:
      break;
  }
}

/* **************************************************************
 * Return a double (scalar) value from a matrix.
 * Only workds if the matrix is 1x1 in size.
 * ************************************************************** */

double mdr_scalar (MDR * m)
{
  if (SIZE(m) != 1)
  {
    fprintf (stderr, "matrix size %i by %i\n", MNR (m), MNC (m));
    rerror ("cannot coerce matrix to scalar");
  }

  return mdrV0 (m, 0);
}


/* **************************************************************
 * Return a character pointer (string).
 * ************************************************************** */

char * mdr_CharPointer (MDR * m)
{
  char *cp, tmp[100];
  double dtmp;
  int itmp;

  if (!m)
    rerror ("terrible internal matrix error");
  if (!MDIPTR(m) && !MDRPTR(m))
    rerror ("cannot coerce empty matrix to character");

  if(m->type == RLAB_TYPE_INT32)
  {
    itmp = Mdi0 (m, 0, 0);
    sprintf (tmp, "%i", itmp);
  }
  else
  {
    dtmp = MdrV0 (m, 0);
    sprintf (tmp, "%.6g", dtmp);
  }
  cp = cpstr (tmp);

  return (cp);
}

/* **************************************************************
 * Return a double value. 1st element, as long as matrix is not
 * empty.
 * ************************************************************** */

double * mdr_Double (MDR * m)
{
  double *rval= (double *) GC_MALLOC (sizeof (double));
  *rval = create_nan();

  if (SIZE(m)==1)
    *rval = mdrV0(m,0);

  return (rval);
}

int * mdr_Integer (MDR * m)
{
  int *ival= (int *) GC_MALLOC (sizeof (int));
  *ival = MAXINT;

  if (SIZE(m)==1)
    *ival = mdiV0(m,0);

  return (ival);
}


MDR * mdr_MatrixDouble (MDR * m)
{
  if (m)
    if (m->type == RLAB_TYPE_DOUBLE)
      return (m);

  return (NULL);
}

MDR * mdr_MatrixReal (MDR * m)
{
  if (m)
    if ((m->type == RLAB_TYPE_DOUBLE)||(m->type == RLAB_TYPE_INT32))
      return (m);

  return (NULL);
}

MDR * mdr_MatrixInteger (MDR * m)
{
  if (m)
    if (m->type == RLAB_TYPE_INT32)
      return (m);

  return (NULL);
}

/* **************************************************************
 * Print out a matrix.
 * ************************************************************** */
#include <sys/ioctl.h>

void
mdr_Print (MDR * matrix, FILE * stream)
{
  char tmp[20];
  int i, j, k, nrow, ncol, npri, rem, start, width;
  int n_print, fwidth, fprec, swidth;

  struct winsize sz;
  ioctl(0, TIOCGWINSZ, &sz);
  swidth = sz.ws_col;
  if (swidth <= 0)
    swidth = TERM_WIDTH;

  fwidth = get_fwidth ();

  fprec = get_fprec ();

  width = MAX(fwidth, fprec + 3) + 2;
  if (width > swidth)
    rerror ("format too large for REAL print");

  n_print = swidth / (width + 1);
  nrow = MNR (matrix);
  ncol = MNC (matrix);
  npri = MNC (matrix) / n_print;
  rem = MNC (matrix) % n_print;

  /* Special case, empty matrix */
  if (nrow == 0 && ncol == 0)
  {
    fprintf (stream, "\t[]\n");
    fflush (stream);
    return;
  }

  start = 1;
  for (i = 0; i < npri; i++)
  {
    if (npri >= 1)
      fprintf (stream, " matrix columns %d thru %d\n",
               n_print * i + 1, (n_print) + (n_print * i));

    for (k = 1; k <= nrow; k++)	/* print all rows */
    {
      for (j = start; j <= n_print + start - 1; j++)
      {
        if(matrix->type == RLAB_TYPE_INT32)
        {
          // print dense integer matrix
          sprintf (tmp, "%*i  ", fwidth, Mdi1 (matrix, k, j));
        }
        else
        {
          // print dense real matrix
          if (isnand(Mdr1 (matrix, k, j)))
            sprintf (tmp, RLAB_SPRINTF_NAN);
          else
            sprintf (tmp, "%*.*g  ", fwidth, fprec, Mdr1 (matrix, k, j));
        }
        fprintf (stream, "%*s", width, tmp);
      }
      fprintf (stream, "\n");
      fflush (stream);
    }
    start += n_print;		/* inc our col position */
    fprintf (stream, "\n");
    fflush (stream);
  }

  /* Now come back and write out the last few columns */
  if (!rem)
    return;
  if (npri >= 1)
    fprintf (stream, " matrix columns %d thru %d\n",
             MNC (matrix) - rem + 1, MNC (matrix));

  for (k = 1; k <= nrow; k++)
  {
    for (i = ncol - rem + 1; i <= ncol; i++)
    {
      if(matrix->type == RLAB_TYPE_INT32)
      {
        // print dense integer matrix
        sprintf (tmp, "%*i  ", fwidth, Mdi1 (matrix, k, i));
      }
      else
      {
        // print dense real matrix
        if (isnand(Mdr1 (matrix, k, i)))
          sprintf (tmp, RLAB_SPRINTF_NAN);
        else
          sprintf (tmp, "%*.*g  ", fwidth, fprec, Mdr1 (matrix, k, i));
      }
      fprintf (stream, "%*s", width, tmp);
      fflush (stream);
    }
    fprintf (stream, "\n");
    fflush (stream);
  }
  return;
}

MDR * mdr_vector_create (MDR * m1, MDR * m2, MDR * m3)
{
  MDR *new;

  double d, d1=0.0, d2=0.0, d3=1.0, dn;
  int i, n;

  if (!m1 || !m2 )
    rerror ("vector-create: UNDEFINED specification");

  /* Extract the scalar values needed... */
  if (SIZE(m1) != 1)
  {
    fprintf (stderr, "cannot use %i by %i matrix in vector stmt\n",
             MNR (m1), MNC (m1));
    rerror ("vector creation error");
  }
  mdr_Detect_Inf (m1);
  mdr_Detect_Nan (m1);
  d1 = mdrV0 (m1, 0);


  if (SIZE(m2) != 1)
  {
    fprintf (stderr, "cannot use %i by %i matrix in vector stmt\n",
             MNR (m2), MNC (m2));
    rerror ("vector creation error");
  }
  mdr_Detect_Inf (m2);
  mdr_Detect_Nan (m2);
  d2 = mdrV0 (m2, 0);

  if (m3)
  {
    mdr_Detect_Inf (m3);
    mdr_Detect_Nan (m3);

    if (SIZE (m3) != 1)
    {
      fprintf (stderr, "cannot use %i by %i matrix in vector stmt\n",
               MNR (m3), MNC (m3));
      rerror ("vector creation error");
    }
    d3 = mdrV0 (m3, 0);

    /* Check for d3 = 0, programmer messed up */
    if (d3 == 0)
      rerror ("must use non-zero increment for vector creation");
  }


  /* Check for condition where we return an empty vector */
  if ((d1 > d2 && d3 > 0.0) || (d1 < d2 && d3 < 0.0))
  {
    new = mdr_Create (0, 0);
    return (new);
  }

  /* Calculate n */
  dn = ABS ((d2 - d1) / d3) + 1;

  /* We could get integer overflow here ! */
  if (dn >= (double) MAXINT)
    rerror ("integer overflow during vector creation");
  n = (int) dn;


  if ((d2 - d1) == 0.0)
    d = d3;
  else
    d = ((d2 - d1) / (ABS (d2 - d1))) * ABS (d3);

  if (MD_TYPE_INT32(m1) && MD_TYPE_INT32(m2) && (d==(int)d))
  {
    /* Fill up vector */
    new = mdi_Create (1, n);
    for (i = 1; i <= n; i++)
      Mdi1 (new, 1, i) = d1 + (i - 1) * d;
  }
  else
  {
    /* Fill up vector */
    new = mdr_Create (1, n);
    for (i = 1; i <= n; i++)
      Mdr1 (new, 1, i) = d1 + (i - 1) * d;
  }

  return (new);
}

/*
 * Same as above, just a more convenient interface for
 * internal usage. I should probably integrate this into
 * the function above.
 */

MDR *
mdr_CreateVector (int d1, int d2, int d3)
{
  double d;
  int i, dn, n;
  MDR *new;

  /* Calculate n */
  dn = ABS ((d2 - d1) / d3) + 1;

  /* We could get integer overflow here ! */
  if (dn >= (double) MAXINT)
    rerror ("integer overflow during vector creation");
  n = (int) dn;

  if ((d2 - d1) == 0.0)
    d = d3;
  else
    d = ((d2 - d1) / (ABS(d2 - d1))) * ABS(d3);
  new = mdr_Create (1, n);

  /* Fill up vector */
  for (i = 1; i <= n; i++)
    Mdr1 (new, 1, i) = d1 + (i - 1) * d;

  return (new);
}

void mdr_BracketInPlace(MDR * new)
{
  if ( rlab_op_max==create_inf() && rlab_op_min==-create_inf() )
    return;

  int size=SIZE(new);
  int i;

  if (size<1)
    return;

  switch(new->type)
  {
    case RLAB_TYPE_INT32:
      if (rlab_op_min==-create_inf())
      {
        for (i=0; i<size; i++)
        {
          MdiV0(new, i) = MdiV0(new, i) > rlab_op_max ? (int) rlab_op_max : MdiV0(new, i);
        }
      }
      else if (rlab_op_max==create_inf())
      {
        for (i=0; i<size; i++)
        {
          MdiV0(new, i) = MdiV0(new, i) < rlab_op_min ? (int) rlab_op_min : MdiV0(new, i);
        }
      }
      else
      {
        for (i=0; i<size; i++)
        {
          MdiV0(new, i) = MdiV0(new, i) < rlab_op_min ? (int) rlab_op_min : MdiV0(new, i);
          MdiV0(new, i) = MdiV0(new, i) > rlab_op_max ? (int) rlab_op_max : MdiV0(new, i);
        }
      }
      break;

    default:
      if (rlab_op_min==-create_inf())
      {
        for (i=0; i<size; i++)
        {
          MdrV0(new, i) = MdrV0(new, i) > rlab_op_max ? (int) rlab_op_max : MdrV0(new, i);
        }
      }
      else if (rlab_op_max==create_inf())
      {
        for (i=0; i<size; i++)
        {
          MdrV0(new, i) = MdrV0(new, i) < rlab_op_min ? (int) rlab_op_min : MdrV0(new, i);
        }
      }
      else
      {
        for (i=0; i<size; i++)
        {
          MdrV0(new, i) = MdrV0(new, i) < rlab_op_min ? (int) rlab_op_min : MdrV0(new, i);
          MdrV0(new, i) = MdrV0(new, i) > rlab_op_max ? (int) rlab_op_max : MdrV0(new, i);
        }
      }
  }

  return;
}


//
// add operation that implements zero abs and rel tolerances
//  useful for adding large matrices, so to avoid additional
//  intermediate creation/deletion of same size matrices
//
#undef  THIS_SOLVER
#define THIS_SOLVER "mdr_Add"
MDR * mdr_Add (MDR * m1, MDR * m2)
{
  MDR *new;
  int i, j, i1, i2, j1, j2;

  i1 = m1->nrow;
  j1 = m1->ncol;
  i2 = m2->nrow ;
  j2 = m2->ncol;

  // Quazi-emulate Matlab empty matrix behavior if both matrices are empty
  if ( i1 * j1 == 0 )
  {
    new = mdr_Copy (m2);
    return (new);
  }
  if ( i2 * j2 == 0 )
  {
    new = mdr_Copy (m1);
    return (new);
  }

  int nr = MAX(i1, i2);
  int nc = MAX(j1, j2);
  i1--;
  i2--;
  j1--;
  j2--;
  if(m1->type == RLAB_TYPE_INT32 && m2->type == RLAB_TYPE_INT32)
  {
    int idummy;
    new = mdi_Create(nr, nc);
    if (rlab_op_zero_abs==0 && rlab_op_zero_rel==0)
    {
      for(i=0;i<nr;i++)
        for(j=0;j<nc;j++)
          Mdi0(new,i,j) = Mdi0(m1, MIN(i,i1), MIN(j,j1)) + Mdi0(m2, MIN(i,i2), MIN(j,j2));
    }
    else if (rlab_op_zero_rel==0)
    {
      for(i=0;i<nr;i++)
        for(j=0;j<nc;j++)
      {
        idummy = Mdi0(m1, MIN(i,i1), MIN(j,j1)) + Mdi0(m2, MIN(i,i2), MIN(j,j2));
        Mdi0(new,i,j) = ABS(idummy) > rlab_op_zero_abs ? idummy : 0;
      }
    }
    else
    {
      int x1, x2;
      double edummy;
      for(i=0;i<nr;i++)
        for(j=0;j<nc;j++)
      {
        x1 = Mdi0(m1, MIN(i,i1), MIN(j,j1));
        x2 = Mdi0(m2, MIN(i,i2), MIN(j,j2));
        idummy = x1 + x2;
        edummy = rlab_op_zero_abs + rlab_op_zero_rel * MAX(ABS(x1),ABS(x2));
        Mdi0(new,i,j) = ABS(idummy) > edummy ? idummy : 0;
      }
    }
  }
  else
  {
    double ddummy;
    new = mdr_Create(nr, nc);
    if (rlab_op_zero_abs==0 && rlab_op_zero_rel==0)
    {
      for(i=0;i<nr;i++)
        for(j=0;j<nc;j++)
          Mdr0(new,i,j) = mdr0(m1, MIN(i,i1), MIN(j,j1)) + mdr0(m2, MIN(i,i2), MIN(j,j2));
    }
    else if(rlab_op_zero_rel==0)
    {
      for(i=0;i<nr;i++)
        for(j=0;j<nc;j++)
      {
        ddummy = mdr0(m1, MIN(i,i1), MIN(j,j1)) + mdr0(m2, MIN(i,i2), MIN(j,j2));
        Mdr0(new,i,j) = ABS(ddummy) > rlab_op_zero_abs ? ddummy : 0.0;
      }
    }
    else
    {
      double x1, x2;
      double edummy;
      for(i=0;i<nr;i++)
        for(j=0;j<nc;j++)
      {
        x1 = mdr0(m1, MIN(i,i1), MIN(j,j1));
        x2 = mdr0(m2, MIN(i,i2), MIN(j,j2));
        ddummy = x1 + x2;
        edummy = rlab_op_zero_abs + rlab_op_zero_rel * MAX(ABS(x1),ABS(x2));
        Mdr0(new,i,j) = ABS(ddummy) > edummy ? ddummy : 0.0;
      }
    }
  }

  if (rlab_op_max!=create_inf() || rlab_op_min!=-create_inf())
  {
    mdr_BracketInPlace(new);
  }

  return (new);
}


//
// add in place m1 -> m1 + m2
//
#undef  THIS_SOLVER
#define THIS_SOLVER "mdr_AddTo"
MDR * mdr_AddTo (MDR * m1, MDR * m2)
{
  MDR *new;
  int i, j, i1, i2, j1, j2;

  i1 = m1->nrow;
  j1 = m1->ncol;
  i2 = m2->nrow ;
  j2 = m2->ncol;

  // Quazi-emulate Matlab empty matrix behavior if both matrices are empty
  if ( i1 * j1 == 0 )
  {
    new = mdr_Copy (m2);
    md_Destroy (m1);
    return (new);
  }
  if ( i2 * j2 == 0 )
  {
    new = mdr_Copy (m1);
    md_Destroy (m1);
    return (new);
  }

  // m2 is bigger in size then m1:
  //  the result will not fit into m1, so we have to create a new matrix
  if ((MNR(m1)<MNR(m2))||(MNC(m1)<MNC(m2))||(m1->type != m2->type))
  {
    new = mdr_Add (m1, m2);
    md_Destroy(m1);
    return (new);
  }

  // m2 is smaller then or equal in size to m1 and has the same type
  //  the result fits into m1, so we reuse it
  new = m1;
  int nr = MAX(i1, i2);
  int nc = MAX(j1, j2);
  i1--;
  i2--;
  j1--;
  j2--;
  int idummy;
  if(m1->type == RLAB_TYPE_INT32)
  {
    if (rlab_op_zero_abs==0 && rlab_op_zero_rel==0)
    {
      for(i=0;i<nr;i++) for(j=0;j<nc;j++)
          Mdi0(new,i,j) = Mdi0(m1, MIN(i,i1), MIN(j,j1)) + Mdi0(m2, MIN(i,i2), MIN(j,j2));
    }
    else if (rlab_op_zero_rel==0)
    {
      for(i=0;i<nr;i++) for(j=0;j<nc;j++)
      {
        idummy = Mdi0(m1, MIN(i,i1), MIN(j,j1)) + Mdi0(m2, MIN(i,i2), MIN(j,j2));
        Mdi0(new,i,j) = ABS(idummy) > rlab_op_zero_abs ? idummy : 0;
      }
    }
    else
    {
      int x1, x2;
      double edummy;
      for(i=0;i<nr;i++) for(j=0;j<nc;j++)
      {
        x1 = Mdi0(m1, MIN(i,i1), MIN(j,j1));
        x2 = Mdi0(m2, MIN(i,i2), MIN(j,j2));
        idummy = x1 + x2;
        edummy = rlab_op_zero_abs + rlab_op_zero_rel * MAX(ABS(x1),ABS(x2));
        Mdi0(new,i,j) = ABS(idummy) > edummy ? idummy : 0;
      }
    }
  }
  else
  {
    double ddummy;
    if (rlab_op_zero_abs==0 && rlab_op_zero_rel==0)
    {
      for(i=0;i<nr;i++) for(j=0;j<nc;j++)
          Mdr0(new,i,j) = mdr0(m1, MIN(i,i1), MIN(j,j1)) + mdr0(m2, MIN(i,i2), MIN(j,j2));
    }
    else if(rlab_op_zero_rel==0)
    {
      for(i=0;i<nr;i++) for(j=0;j<nc;j++)
      {
        ddummy = mdr0(m1, MIN(i,i1), MIN(j,j1)) + mdr0(m2, MIN(i,i2), MIN(j,j2));
        Mdr0(new,i,j) = ABS(ddummy) > rlab_op_zero_abs ? ddummy : 0.0;
      }
    }
    else
    {
      double x1, x2;
      double edummy;
      for(i=0;i<nr;i++) for(j=0;j<nc;j++)
      {
        x1 = mdr0(m1, MIN(i,i1), MIN(j,j1));
        x2 = mdr0(m2, MIN(i,i2), MIN(j,j2));
        ddummy = x1 + x2;
        edummy = rlab_op_zero_abs + rlab_op_zero_rel * MAX(ABS(x1),ABS(x2));
        Mdr0(new,i,j) = ABS(ddummy) > edummy ? ddummy : 0.0;
      }
    }
  }

  if (rlab_op_max!=create_inf() || rlab_op_max!=-create_inf())
    mdr_BracketInPlace(new);

  return (new);
}

MDR * mdr_ElMultiplyBy(MDR * m1, MDR * m2)
{
  MDR *new;
  int i, j, i1, i2, j1, j2;

  i1 = m1->nrow;
  j1 = m1->ncol;
  i2 = m2->nrow ;
  j2 = m2->ncol;

  // Quazi-emulate Matlab empty matrix behavior if both matrices are empty
  if ( i1 * j1 == 0 )
  {
    new = mdr_Copy (m2);
    md_Destroy (m1);
    return (new);
  }
  if ( i2 * j2 == 0 )
  {
    new = mdr_Copy (m1);
    md_Destroy (m1);
    return (new);
  }

  // m2 is bigger in size then m1:
  //  the result will not fit into m1, so we have to create a new matrix
  if ((MNR(m1)<MNR(m2))||(MNC(m1)<MNC(m2))||(m1->type != m2->type))
  {
    int nr = MAX(i1, i2);
    int nc = MAX(j1, j2);
    i1--;
    i2--;
    j1--;
    j2--;
    if(m1->type == RLAB_TYPE_INT32 && m2->type == RLAB_TYPE_INT32)
    {
      new = mdi_Create(nr, nc);
      for(i=0;i<nr;i++)
        for(j=0;j<nc;j++)
          Mdi0(new,i,j) = Mdi0(m1, MIN(i,i1), MIN(j,j1)) * Mdi0(m2, MIN(i,i2), MIN(j,j2));
    }
    else
    {
      new = mdr_Create(nr, nc);
      for(i=0;i<nr;i++)
        for(j=0;j<nc;j++)
          Mdr0(new,i,j) = mdr0(m1, MIN(i,i1), MIN(j,j1)) * mdr0(m2, MIN(i,i2), MIN(j,j2));
    }
    mdr_Destroy(m1);
  }
  else
  {
    // m2 is smaller then or equal in size to m1 and has the same type
    //  the result will not fit into m1, so we have to create a new matrix
    new = m1;
    int nr = MAX(i1, i2);
    int nc = MAX(j1, j2);
    i1--;
    i2--;
    j1--;
    j2--;
    if(m1->type == RLAB_TYPE_INT32)
    {
      for(i=0;i<nr;i++)
        for(j=0;j<nc;j++)
          Mdi0(new,i,j) *= Mdi0(m2, MIN(i,i2), MIN(j,j2));
    }
    else
    {
      for(i=0;i<nr;i++)
        for(j=0;j<nc;j++)
          Mdr0(new,i,j) *= Mdr0(m2, MIN(i,i2), MIN(j,j2));
    }
  }

  return (new);
}

MDR * mdr_ElRdivideBy(MDR * m1, MDR * m2)
{
  MDR *new;
  int i, j, i1, i2, j1, j2;

  i1 = m1->nrow;
  j1 = m1->ncol;
  i2 = m2->nrow ;
  j2 = m2->ncol;

  // Quazi-emulate Matlab empty matrix behavior if both matrices are empty
  if ( i1 * j1 == 0 )
  {
    new = mdr_Copy (m2);
    md_Destroy (m1);
    return (new);
  }
  if ( i2 * j2 == 0 )
  {
    new = mdr_Copy (m1);
    md_Destroy (m1);
    return (new);
  }

  // m2 is bigger in size then m1:
  //  the result will not fit into m1, so we have to create a new matrix
  if ((MNR(m1)<MNR(m2))||(MNC(m1)<MNC(m2))||(m1->type != m2->type || m1->type==RLAB_TYPE_INT32))
  {
    int nr = MAX(i1, i2);
    int nc = MAX(j1, j2);
    i1--;
    i2--;
    j1--;
    j2--;
    if(m1->type == RLAB_TYPE_INT32 && m2->type == RLAB_TYPE_INT32)
    {
      new = mdi_Create(nr, nc);
      for(i=0;i<nr;i++)
        for(j=0;j<nc;j++)
          Mdi0(new,i,j) = Mdi0(m1, MIN(i,i1), MIN(j,j1)) / Mdi0(m2, MIN(i,i2), MIN(j,j2));
    }
    else
    {
      new = mdr_Create(nr, nc);
      for(i=0;i<nr;i++)
        for(j=0;j<nc;j++)
          Mdr0(new,i,j) = mdr0(m1, MIN(i,i1), MIN(j,j1)) / mdr0(m2, MIN(i,i2), MIN(j,j2));
    }
    mdr_Destroy(m1);
  }
  else
  {
    // m2 is smaller then or equal in size to m1 and has the same type
    //  the result will not fit into m1, so we have to create a new matrix
    new = m1;
    int nr = MAX(i1, i2);
    int nc = MAX(j1, j2);
    i1--;
    i2--;
    j1--;
    j2--;
    if(m1->type == RLAB_TYPE_INT32)
    {
      for(i=0;i<nr;i++)
        for(j=0;j<nc;j++)
          Mdi0(new,i,j) /= Mdi0(m2, MIN(i,i2), MIN(j,j2));
    }
    else
    {
      for(i=0;i<nr;i++)
        for(j=0;j<nc;j++)
          Mdr0(new,i,j) /= Mdr0(m2, MIN(i,i2), MIN(j,j2));
    }
  }

  return (new);
}

MDR *
mdr_Add_weight (MDR * m1, MDR *w1, MDR * m2, MDR *w2)
{
  MDR *new=0;
  double dummy;
  int i, j, i1, i2, j1, j2;
  int iw1=0, jw1=0, iw2=0, jw2=0, nr, nc;

  i1 = m1->nrow;
  j1 = m1->ncol;
  i2 = m2->nrow ;
  j2 = m2->ncol;
  nr = MAX(i1, i2);
  nc = MAX(j1, j2);

  if (w1)
  {
    iw1 = MNR(w1);
    jw1 = MNC(w1);
  }
  if (w2)
  {
    iw2 = MNR(w2);
    jw2 = MNC(w2);
  }

  // Quazi-emulate Matlab empty matrix behavior if both matrices are empty
  if ( i1 * j1 == 0 )
  {
    new = mdr_Copy (w2);
    return (new);
  }
  if ( i2 * j2 == 0 )
  {
    new = mdr_Copy (w1);
    return (new);
  }

  // do the matrix optimization thingy
  // adjust indices for MIN/max operations
  i1--;
  j1--;
  i2--;
  j2--;
  iw1--;
  jw1--;
  iw2--;
  jw2--;

  new = mdr_Create(nr, nc);

  if (w1 || w2)
  {
    if (w1 && w2)
    {
      for(i=0;i<nr;i++)
        for(j=0;j<nc;j++)
        {
          dummy =  1 / Mdr0(w1,MIN(i,iw1),MIN(j,jw1))
            + 1 / Mdr0(w2,MIN(i,iw2),MIN(j,jw2));
          Mdr0(new,i,j) = 1 / dummy;
        }
    }
    else if (w1)
    {
      for(i=0;i<nr;i++)
        for(j=0;j<nc;j++)
        {
          Mdr0(new,i,j) = Mdr0(w1,MIN(i,iw1),MIN(j,jw1));
        }
    }
    else if (w2)
    {
      for(i=0;i<nr;i++)
        for(j=0;j<nc;j++)
        {
          Mdr0(new,i,j) = Mdr0(w2,MIN(i,iw2),MIN(j,jw2));
        }
    }
  }

  return (new);
}

//
// Subtract two matrices.
//
MDR * mdr_Subtract (MDR * m1, MDR * m2)
{
  MDR *new;
  int i, j, i1, i2, j1, j2;

  i1 = m1->nrow;
  j1 = m1->ncol;
  i2 = m2->nrow ;
  j2 = m2->ncol;

  // Quazi-emulate Matlab empty matrix behavior if both matrices are empty
  if ( i1 * j1 == 0 )
  {
    new = mdr_Copy (m2);
    return (new);
  }
  if ( i2 * j2 == 0 )
  {
    new = mdr_Copy (m1);
    return (new);
  }

  int nr = MAX(i1, i2);
  int nc = MAX(j1, j2);
  i1--;
  i2--;
  j1--;
  j2--;
  if(m1->type == RLAB_TYPE_INT32 && m2->type == RLAB_TYPE_INT32)
  {
    int idummy;
    new = mdi_Create(nr, nc);
    if (rlab_op_zero_abs==0 && rlab_op_zero_rel==0)
    {
      for(i=0;i<nr;i++)
        for(j=0;j<nc;j++)
          Mdi0(new,i,j) = Mdi0(m1, MIN(i,i1), MIN(j,j1)) - Mdi0(m2, MIN(i,i2), MIN(j,j2));
    }
    else if (rlab_op_zero_rel==0)
    {
      for(i=0;i<nr;i++)
        for(j=0;j<nc;j++)
      {
        idummy = Mdi0(m1, MIN(i,i1), MIN(j,j1)) - Mdi0(m2, MIN(i,i2), MIN(j,j2));
        Mdi0(new,i,j) = ABS(idummy) > rlab_op_zero_abs ? idummy : 0;
      }
    }
    else
    {
      int x1, x2;
      double edummy;
      for(i=0;i<nr;i++)
        for(j=0;j<nc;j++)
      {
        x1 = Mdi0(m1, MIN(i,i1), MIN(j,j1));
        x2 = Mdi0(m2, MIN(i,i2), MIN(j,j2));
        idummy = x1 - x2;
        edummy = rlab_op_zero_abs + rlab_op_zero_rel * MAX(ABS(x1),ABS(x2));
        Mdi0(new,i,j) = ABS(idummy) > edummy ? idummy : 0;
      }
    }
  }
  else
  {
    double ddummy;
    new = mdr_Create(nr, nc);
    if (rlab_op_zero_abs==0 && rlab_op_zero_rel==0)
    {
      for(i=0;i<nr;i++)
        for(j=0;j<nc;j++)
          Mdr0(new,i,j) = mdr0(m1, MIN(i,i1), MIN(j,j1)) - mdr0(m2, MIN(i,i2), MIN(j,j2));
    }
    else if(rlab_op_zero_rel==0)
    {
      for(i=0;i<nr;i++)
        for(j=0;j<nc;j++)
      {
        ddummy = mdr0(m1, MIN(i,i1), MIN(j,j1)) - mdr0(m2, MIN(i,i2), MIN(j,j2));
        Mdr0(new,i,j) = ABS(ddummy) > rlab_op_zero_abs ? ddummy : 0.0;
      }
    }
    else
    {
      double x1, x2;
      double edummy;
      for(i=0;i<nr;i++)
        for(j=0;j<nc;j++)
      {
        x1 = mdr0(m1, MIN(i,i1), MIN(j,j1));
        x2 = mdr0(m2, MIN(i,i2), MIN(j,j2));
        ddummy = x1 - x2;
        edummy = rlab_op_zero_abs + rlab_op_zero_rel * MAX(ABS(x1),ABS(x2));
        Mdr0(new,i,j) = ABS(ddummy) > edummy ? ddummy : 0.0;
      }
    }
  }

  if (rlab_op_max!=create_inf() || rlab_op_min!=-create_inf())
  {
    mdr_BracketInPlace(new);
  }

  return (new);
}

//
// Add two matrices.
//
MDR * mdr_SubtractFrom(MDR * m1, MDR * m2)
{
  MDR *new;
  int i, j, i1, i2, j1, j2;

  i1 = m1->nrow;
  j1 = m1->ncol;
  i2 = m2->nrow ;
  j2 = m2->ncol;

  // Quazi-emulate Matlab empty matrix behavior if both matrices are empty
  if ( i1 * j1 == 0 )
  {
    new = mdr_Copy (m2);
    md_Destroy (m1);
    return (new);
  }
  if ( i2 * j2 == 0 )
  {
    new = mdr_Copy (m1);
    md_Destroy (m1);
    return (new);
  }

  // m2 is bigger in size then m1:
  //  the result will not fit into m1, so we have to create a new matrix
  if ((MNR(m1)<MNR(m2))||(MNC(m1)<MNC(m2))||(m1->type != m2->type))
  {
    new = mdr_Subtract (m1, m2);
    mdr_Destroy(m1);
    return (new);
  }

  // m2 is smaller then or equal in size to m1 and has the same type
  //  the result fits into m1, so we reuse it
  new = m1;
  int nr = MAX(i1, i2);
  int nc = MAX(j1, j2);
  i1--;
  i2--;
  j1--;
  j2--;
  int idummy;
  if(m1->type == RLAB_TYPE_INT32)
  {
    if (rlab_op_zero_abs==0 && rlab_op_zero_rel==0)
    {
      for(i=0;i<nr;i++) for(j=0;j<nc;j++)
        Mdi0(new,i,j) -= mdi0_safe(m2,i, j);
    }
    else if (rlab_op_zero_rel==0)
    {
      for(i=0;i<nr;i++) for(j=0;j<nc;j++)
      {
        idummy = Mdi0(m1, MIN(i,i1), MIN(j,j1)) - Mdi0(m2, MIN(i,i2), MIN(j,j2));
        Mdi0(new,i,j) = ABS(idummy) > rlab_op_zero_abs ? idummy : 0;
      }
    }
    else
    {
      int x1, x2;
      double edummy;
      for(i=0;i<nr;i++) for(j=0;j<nc;j++)
      {
        x1 = mdi0_safe(m1,i,j);
        x2 = mdi0_safe(m2,i,j);
        idummy = x1 - x2;
        edummy = rlab_op_zero_abs + rlab_op_zero_rel * MAX(ABS(x1),ABS(x2));
        Mdi0(new,i,j) = ABS(idummy) > edummy ? idummy : 0;
      }
    }
  }
  else
  {
    double ddummy;
    if (rlab_op_zero_abs==0 && rlab_op_zero_rel==0)
    {
      for(i=0;i<nr;i++) for(j=0;j<nc;j++)
        Mdr0(new,i,j) -= mdr0_safe(m2, i,j);
    }
    else if(rlab_op_zero_rel==0)
    {
      for(i=0;i<nr;i++) for(j=0;j<nc;j++)
      {
        ddummy = mdr0(m1, MIN(i,i1), MIN(j,j1)) - mdr0(m2, MIN(i,i2), MIN(j,j2));
        Mdr0(new,i,j) = ABS(ddummy) > rlab_op_zero_abs ? ddummy : 0.0;
      }
    }
    else
    {
      double x1, x2;
      double edummy;
      for(i=0;i<nr;i++) for(j=0;j<nc;j++)
      {
        x1 = mdr0(m1, MIN(i,i1), MIN(j,j1));
        x2 = mdr0(m2, MIN(i,i2), MIN(j,j2));
        ddummy = x1 - x2;
        edummy = rlab_op_zero_abs + rlab_op_zero_rel * MAX(ABS(x1),ABS(x2));
        Mdr0(new,i,j) = ABS(ddummy) > edummy ? ddummy : 0.0;
      }
    }
  }

  if (rlab_op_max!=create_inf() || rlab_op_max!=-create_inf())
    mdr_BracketInPlace(new);

  return (new);
}




//
// Multiply two matrices: mc = ma * mb
//
MDR *
mdr_Multiply (MDR * ma, MDR * mb, int *rtype)
{
  MDR *mc = 0;
  int k, isize;
  double alpha = 1.0;
  double beta = 0.0;

  *rtype = MATRIX_DENSE_REAL;

  /* Check [ma], [mb] dimensions */
  if (MNC (ma) != MNR (mb))
  {
    /*
     * Handle condition where one of the operands is a
     * scalar, or an empty matrix.
     */

    if (SIZE(ma) == 1)
    {
      isize = MNR (mb) * MNC (mb);
      if (ma->type == RLAB_TYPE_INT32 && mb->type == RLAB_TYPE_INT32)
      {
        mc = mdi_Create (MNR (mb), MNC (mb));
        for (k = 0; k < isize; k++)
          MdiV0 (mc, k) = MdiV0 (ma, 0) * MdiV0 (mb, k);
      }
      else
      {
        alpha = mdrV0(ma,0);
        mc = mdr_Create (MNR (mb), MNC (mb));
        for (k = 0; k < isize; k++)
          MdrV0 (mc, k) = alpha * mdrV0 (mb, k);
      }
    }
    else if (MNR (mb) == 1 && MNC (mb) == 1)
    {
      isize = MNR (ma) * MNC (ma);
      if (ma->type == RLAB_TYPE_INT32 && mb->type == RLAB_TYPE_INT32)
      {
        mc = mdi_Create (MNR (ma), MNC (ma));
        for (k = 0; k < isize; k++)
          MdiV0 (mc, k) = MdiV0 (ma, k) * MdiV0 (mb, 0);
      }
      else
      {
        alpha = mdrV0(mb,0);
        mc = mdr_Create (MNR (ma), MNC (ma));
        for (k = 0; k < isize; k++)
          MdrV0 (mc, k) = mdrV0 (ma, k) * alpha;
      }
    }
    else if (MNR (ma) == 0 || MNR (mb) == 0)
    {
      mc = mdr_Create (0, 0);
    }
    else if (MNR (mb) == 0 || MNR (mb) == 0)
    {
      mc = mdr_Create (0, 0);
    }
    else
    {
      fprintf (stderr, "\tmatrix dimensions must be consistent\n");
      fprintf (stderr, "\tmatrix A: row = %i, column = %i\n", MNR (ma),
	       MNC (ma));
      fprintf (stderr, "\tmatrix B: row = %i, column = %i\n", MNR (mb),
	       MNC (mb));
      rerror ("matrix-multiplication dimension mis-match");
    }
  }
  else
  {
    int ic, jc, kc, m, n, k;
    char tra, trb;
    m = (int) MNR (ma);
    n = (int) MNC (mb);
    k = (int) MNC (ma);
    if (ma->type == RLAB_TYPE_INT32 && mb->type == RLAB_TYPE_INT32)
    {
      // non-optimized integer multiplication
      mc = mdi_Create (MNR (ma), MNC (mb));
      for(ic=1; ic<=m; ic++) for(jc=1; jc<=k; jc++)
      {
        Mdi1(mc,ic,jc) = 0;
        for (kc=1; kc <= n; kc++)
        {
          Mdi1(mc,ic,jc) += Mdi1(ma,ic,kc) * Mdi1(mb,kc,jc);
        }
      }
    }
    else
    {
      MDR *da, *db;
      if (ma->type == RLAB_TYPE_INT32)
        da = mdr_Float_BF (ma);
      else
        da = ma;
      if (mb->type == RLAB_TYPE_INT32)
        db = mdr_Float_BF (mb);
      else
        db = mb;
      mc = mdr_Create (MNR (da), MNC (mb));
      mdr_Zero (mc);
      tra = 'N';
      trb = 'N';

      /* BLAS option: use `double alpha = 1.0, beta = 0.0' */
      DGEMM (&tra, &trb, &m, &n, &k, &alpha, MDRPTR (da),
              &m, MDRPTR (db), &k, &beta, MDRPTR (mc), &m);

      if (ma->type == RLAB_TYPE_INT32)
        mdr_Destroy (da);
      if (mb->type == RLAB_TYPE_INT32)
        mdr_Destroy (db);
    }
  }
  return (mc);
}

MDR *
mdr_Multiply_weight (MDR * ma, MDR *wa, MDR * mb, MDR *wb)
{
  MDR *mc = 0;

  int iwa=0, jwa=0, iwb=0, jwb=0;
  int ia, ja, ib, jb, i, j, k;

  double dummy, dummya, dummyb;

  ia = MNR(ma);
  ja = MNC(ma);
  ib = MNR(mb);
  jb = MNC(mb);

  if (wa)
  {
    iwa = MNR(wa);
    jwa = MNC(wa);
  }
  if (wb)
  {
    iwb = MNR(wb);
    jwb = MNC(wb);
  }

  // Quazi-emulate Matlab empty matrix behavior if both matrices are empty
  if ( ia * ja == 0 )
  {
    mc = mdr_Copy (wb);
    return (mc);
  }
  if ( ib * jb == 0 )
  {
    mc = mdr_Copy (wb);
    return (mc);
  }

  // if ma or mb is 1-by-1 do the element multiplication
  if (ia * ja == 1 || ib * jb == 1)
    return mdr_ElMultiply_weight (ma, wa, mb, wb);

  // we should have full two matrices here
  if (ib!=ja)
    rerror("Terrible internal error. Cannot continue!");

  mc = mdr_Create (ia, jb);

  // do the matrix optimization thingy
  // adjust indices for MIN/max operations
  iwa--;
  jwa--;
  iwb--;
  jwb--;

  if (wa && wb)
  {
    for (i=0; i<ia; i++)
      for (j=0; j<jb; j++)
      {
        dummy = 0;
        for (k=0; k<ja; k++)
        {
          dummya  = mdr0(ma,i,k);
          dummya *= dummya;
          dummyb  = mdr0(mb,k,j);
          dummyb *= dummyb;
          dummy  += dummya / Mdr0(wb,MIN(k,iwb),MIN(j,jwb))
          + dummyb / Mdr0(wa,MIN(i,iwa),MIN(k,jwa));
        }
        Mdr0(mc,i,j) = 1 / dummy;
      }
  }
  else if (wa)
  {
    for (i=0; i<ia; i++)
      for (j=0; j<jb; j++)
      {
        dummy = 0;
        for (k=0; k<ja; k++)
        {
          dummyb  = mdr0(mb,k,j);
          dummyb *= dummyb;
          dummy  += dummyb / Mdr0(wa,MIN(i,iwa),MIN(k,jwa));
        }
        Mdr0(mc,i,j) = 1 / dummy;
      }
  }
  else if (wb)
  {
    for (i=0; i<ia; i++)
      for (j=0; j<jb; j++)
      {
        dummy = 0;
        for (k=0; k<ja; k++)
        {
          dummya  = mdr0(ma,i,k);
          dummya *= dummya;
          dummy += dummya / Mdr0(wb,MIN(k,iwb),MIN(j,jwb));
        }
        Mdr0(mc,i,j) = 1 / dummy;
      }
  }

  return (mc);
}

//
// Element-by-Element Multiply two matrices: mc = ma * mb
//
MDR *
mdr_ElMultiply (MDR * m1, MDR * m2)
{
  MDR *new=0;
  int i, j, i1, i2, j1, j2;

  i1 = m1->nrow;
  j1 = m1->ncol;
  i2 = m2->nrow ;
  j2 = m2->ncol;

  // Quazi-emulate Matlab empty matrix behavior if both matrices are empty
  if ( i1 * j1 == 0 )
  {
    new = mdr_Copy (m2);
    return (new);
  }
  if ( i2 * j2 == 0 )
  {
    new = mdr_Copy (m1);
    return (new);
  }

  // do the matrix optimization thingy
  int nr = MAX(i1, i2);
  int nc = MAX(j1, j2);
  i1--;
  i2--;
  j1--;
  j2--;

  new = mdr_Create(nr, nc);
  for(i=0;i<nr;i++)
    for(j=0;j<nc;j++)
      Mdr0(new,i,j) = mdr0(m1, MIN(i,i1), MIN(j,j1)) * mdr0(m2, MIN(i,i2), MIN(j,j2));

  return (new);
}

MDR *
mdr_ElMultiply_weight (MDR * m1, MDR *w1, MDR * m2, MDR *w2)
{
  MDR *new=0;
  double x=0,vx=0,y=0,vy=0;
  int i, j, i1, i2, j1, j2, nr, nc;
  int iw1=0, jw1=0, iw2=0, jw2=0;

  i1 = m1->nrow;
  j1 = m1->ncol;
  i2 = m2->nrow ;
  j2 = m2->ncol;

  if (w1)
  {
    iw1 = MNR(w1);
    jw1 = MNC(w1);
  }
  if (w2)
  {
    iw2 = MNR(w2);
    jw2 = MNC(w2);
  }

  // Quazi-emulate Matlab empty matrix behavior if both matrices are empty
  if ( i1 * j1 == 0 )
  {
    new = mdr_Copy (w2);
    return (new);
  }
  if ( i2 * j2 == 0 )
  {
    new = mdr_Copy (w1);
    return (new);
  }

  nr = MAX(i1, i2);
  nc = MAX(j1, j2);
  new = mdr_Create(nr, nc);

  // adjust indices for MIN/max operations
  i1--;
  j1--;
  i2--;
  j2--;
  iw1--;
  jw1--;
  iw2--;
  jw2--;

  if (w1 || w2)
  {
    if (w1 && w2)
    {
      for(i=0;i<nr;i++)
        for(j=0;j<nc;j++)
        {
          // x
          x  = mdr0(m1,MIN(i,i1), MIN(j,j1));
          vx = mdr0(w1,MIN(i,iw1),MIN(j,jw1));
          // y
          y  = mdr0(m2,MIN(i,i2), MIN(j,j2));
          vy = mdr0(w2,MIN(i,iw2),MIN(j,jw2));
          //
          Mdr0(new,i,j) = 1 / (gsl_pow_2(y) * vx  + gsl_pow_2(x) * vy);
        }
    }
    else if (w1)
    {
      for(i=0;i<nr;i++)
        for(j=0;j<nc;j++)
        {
          // x
          vx = mdr0(w1,MIN(i,iw1),MIN(j,jw1));
          // y
          y  = mdr0(m2,MIN(i,i2), MIN(j,j2));
          //
          Mdr0(new,i,j) = 1 / (gsl_pow_2(y) * vx);
        }
    }
    else if (w2)
    {
      for(i=0;i<nr;i++)
        for(j=0;j<nc;j++)
        {
          // x
          x  = mdr0(m1,MIN(i,i1), MIN(j,j1));
          // y
          vy = mdr0(w2,MIN(i,iw2),MIN(j,jw2));
          //
          Mdr0(new,i,j) = 1 / (gsl_pow_2(x) * vy);
        }
    }
  }

  return (new);
}


/*
 * m1 / m2 or B / A
 * same as B*inv(A)
 */

MDR *
mdr_Rdivide (MDR * M1, MDR * M2)
{
  double dtmp;
  int i, size;
  MDR *new, *tnew, *m1, *m2, *t1, *t2;

  /*
   * Check for special case where denominator is a scalar.
   * The only thing we gain if m2 is a scalar, is speed.
   */

  // not an integer optimized
  if(M1->type == RLAB_TYPE_INT32)
    m1 = mdr_Float_BF (M1);
  else
    m1 = M1;
  if(M2->type == RLAB_TYPE_INT32)
    m2 = mdr_Float_BF (M2);
  else
    m2 = M2;

  if (MNR (m2) == 1 && MNC (m2) == 1)
  {
    new = mdr_Create (MNR (m1), MNC (m1));
    size = MNR (m1) * MNC (m1);
    dtmp = MdrV0 (m2, 0);
    for (i = 0; i < size; i++)
    {
      MdrV0 (new, i) = MdrV0 (m1, i) / dtmp;
    }
    return (new);
  }
  else if (MNR (m1) == 0 || MNR (m2) == 0)
  {
    new = mdr_Create (0, 0);
    return (new);
  }

  t1 = mdr_Transpose (m1);
  t2 = mdr_Transpose (m2);

  /* Check row dimensions */
  if (MNR (t1) != MNR (t2))
  {
    mdr_Destroy (t1);
    mdr_Destroy (t2);
    rerror ("column dimensions must agree for right-divide");
  }

  if (MNR (t2) == MNC (t2))
    tnew = mdr_SolveEq (t2, t1);
  else
    tnew = mdr_LS (t2, t1);

  new = mdr_Transpose (tnew);
  mdr_Destroy (t1);
  mdr_Destroy (t2);
  mdr_Destroy (tnew);

  // not an integer optimized
  if(M1->type == RLAB_TYPE_INT32)
    mdr_Destroy(m1);
  if(M2->type == RLAB_TYPE_INT32)
    mdr_Destroy(m2);

  return (new);
}

//
// Element right-divide for two matrices.
//  m1 ./ m2
//
#undef  THIS_SOLVER
#define THIS_SOLVER "mdr_ElRdivide"
MDR * mdr_ElRdivide (MDR * m1, MDR * m2, int *rtype)
{
  MDR *new=0;
  int i, j, i1, i2, j1, j2;

  i1 = m1->nrow;
  j1 = m1->ncol;
  i2 = m2->nrow ;
  j2 = m2->ncol;

  *rtype = MATRIX_DENSE_REAL;

  // Quazi-emulate Matlab empty matrix behavior if both matrices are empty
  if (!SIZE(m1))
  {
    new = mdr_Copy (m2);
    return (new);
  }
  if (!SIZE(m2))
  {
    new = mdr_Copy (m1);
    return (new);
  }

  // do the matrix optimization thingy
  int nr = MAX(i1, i2);
  int nc = MAX(j1, j2);
  i1--;
  i2--;
  j1--;
  j2--;
  new = mdr_Create(nr, nc);

  for(i=0;i<nr;i++)
    for(j=0;j<nc;j++)
      Mdr0(new,i,j) = mdr0(m1, MIN(i,i1), MIN(j,j1)) / mdr0(m2, MIN(i,i2), MIN(j,j2));

  return (new);
}

//
// Element right-divide for two matrices.
//  m1 ./ m2
//
MDR *
mdr_ElRdivide_weight (MDR * m1, MDR *w1, MDR * m2, MDR *w2)
{
  MDR *new=0;
  double x=0,vx=0,y=0,vy=0;
  int i, j, i1, i2, j1, j2;
  int iw1=0, jw1=0, iw2=0, jw2=0, nr, nc;

  i1 = m1->nrow;
  j1 = m1->ncol;
  i2 = m2->nrow ;
  j2 = m2->ncol;
  nr = MAX(i1, i2);
  nc = MAX(j1, j2);

  if (w1)
  {
    iw1 = MNR(w1);
    jw1 = MNC(w1);
  }
  if (w2)
  {
    iw2 = MNR(w2);
    jw2 = MNC(w2);
  }

  // Quazi-emulate Matlab empty matrix behavior if both matrices are empty
  if ( i1 * j1 == 0 )
  {
    new = mdr_Copy (w2);
    return (new);
  }
  if ( i2 * j2 == 0 )
  {
    new = mdr_Copy (w1);
    return (new);
  }

  // do the matrix optimization thingy
  // adjust indices for min/max operations
  i1--;
  j1--;
  i2--;
  j2--;
  iw1--;
  jw1--;
  iw2--;
  jw2--;

  new = mdr_Create(nr, nc);

  if (w1 || w2)
  {
    if (w1 && w2)
    {
      for(i=0;i<nr;i++)
        for(j=0;j<nc;j++)
        {
          // x
          x  = mdr0(m1,MIN(i,i1), MIN(j,j1));
          vx = mdr0(w1,MIN(i,iw1),MIN(j,jw1));
          // y
          y  = mdr0(m2,MIN(i,i2), MIN(j,j2));
          vy = mdr0(w2,MIN(i,iw2),MIN(j,jw2));
          //
          Mdr0(new,i,j) = 1 / (vx / gsl_pow_2(y) + gsl_pow_2(x) * vy / gsl_pow_4(y));
        }
    }
    else if (w1)
    {
      for(i=0;i<nr;i++)
        for(j=0;j<nc;j++)
        {
          // x
          vx = mdr0(w1,MIN(i,iw1),MIN(j,jw1));
          // y
          y  = mdr0(m2,MIN(i,i2), MIN(j,j2));
          //
          Mdr0(new,i,j) =  gsl_pow_2(y) / vx;
        }
    }
    else if (w2)
    {
      for(i=0;i<nr;i++)
        for(j=0;j<nc;j++)
        {
          // x
          x  = mdr0(m1,MIN(i,i1), MIN(j,j1));
          // y
          y  = mdr0(m2,MIN(i,i2), MIN(j,j2));
          vy = mdr0(w2,MIN(i,iw2),MIN(j,jw2));
          //
          Mdr0(new,i,j) = gsl_pow_4(y) / gsl_pow_2(x) / vy;
        }
    }
  }

  return (new);
}


/* **************************************************************
 * Matrix left-divide
 *
 * m1 \ m2 or A \ B
 * Ax = B
 * ************************************************************** */

MDR *
mdr_Ldivide (MDR * M1, MDR * M2)
{
  double dtmp;
  int i, size;
  MDR *new, *m1, *m2;

  // not an integer optimized
  if(M1->type == RLAB_TYPE_INT32)
    m1 = mdr_Float_BF (M1);
  else
    m1 = M1;
  if(M2->type == RLAB_TYPE_INT32)
    m2 = mdr_Float_BF (M2);
  else
    m2 = M2;

  /*
   * Check for special case where denominator (m1) is a scalar.
   * The only thing we gain if m2 is a scalar, is speed.
   */

  if (MNR (m1) == 1 && MNC (m1) == 1)
  {
    new = mdr_Create (MNR (m2), MNC (m2));
    size = MNR (m2) * MNC (m2);
    dtmp = MdrV0 (m1, 0);
    for (i = 0; i < size; i++)
    {
      MdrV0 (new, i) = MdrV0 (m2, i) / dtmp;
    }
    return (new);
  }
  else if (MNR (m1) == 0 || MNR (m2) == 0)
  {
    new = mdr_Create (0, 0);
    return (new);
  }

  /* Check dimensions */
  if (MNR (m1) != MNR (m2))
    rerror ("RHS row dim. must match LHS row dim.");

  if (MNR (m1) == MNC (m1))
  {
    new = mdr_SolveEq_GE (m1, m2);
  }
  else
  {
    new = mdr_LS (m1, m2);
  }

  // not an integer optimized
  if(M1->type == RLAB_TYPE_INT32)
    mdr_Destroy(m1);
  if(M2->type == RLAB_TYPE_INT32)
    mdr_Destroy(m2);

  return (new);
}

/* **************************************************************
 * Matrix Element Left Divide ( m1 .\  m2 )
 * ************************************************************** */

MDR *
mdr_ElLdivide (MDR * m1, MDR * m2)
{
  MDR *new;
  int i, size;

  if (MNR (m1) != MNR (m2))
  {
    fprintf (stderr, "matrix row sizes must be equal\n");
    fprintf (stderr, "matrix row sizes: %i and %i\n", MNR (m1), MNR (m2));
    rerror ("matrix-element-left-divide: row size mis-match");
  }

  if (MNC (m1) != MNC (m2))
  {
    fprintf (stderr, "matrix column sizes must be equal\n");
    fprintf (stderr, "matrix column sizes: %i and %i\n", MNC (m1), MNC (m2));
    rerror ("matrix-element-left-divide: column size mis-match");
  }

  new = mdr_Create (MNR (m1), MNC (m1));
  size = MNR (m1) * MNC (m1);

  for (i = 0; i < size; i++)
    MdrV0(new,i) = mdrV0(m2,i) / mdrV0(m1,i);

  return (new);
}


/* **************************************************************
 * Matrix Power operations...
 * ************************************************************** */

/*
 * Scalar ^ Matrix
 * Do the scalar ^ [1-by-1] case only.
 */

void *mdr_elpower1 (MDR * m1, MDR * m2, int *type);

void *
mdr_power1 (MDR * m1, MDR * m2, int *type)
{
  void *m = 0;

  /* Special case: scalar ^ scalar */
  if ((MNR (m2) == 1) && (MNC (m2) == 1))
  {
    m = mdr_elpower1 (m1, m2, type);
    return (m);
  }

  rerror ("scalar ^ matrix not supported yet");
  return (m);
}

/*
 * Matrix ^ Scalar
 *
 * If scalar is an integer, use matrix multiplies.
 * If scalar is not integer use eigenvalues/vectors.
 */

void *
mdr_power2 (MDR * m1, MDR * m2, int *type)
{
  void *m = 0;

  /* Special case: scalar ^ scalar */
  if ((MNR (m1) == 1) && (MNC (m1) == 1))
  {
    m = mdr_elpower1 (m1, m2, type);
    return (m);
  }

  rerror ("matrix ^ scalar not supported yet");
  return (m);
}

/*
 * Matrix ^ Matrix
 */

void *
mdr_power3 (MDR * m1, MDR * m2, int *type)
{
  void *m = 0;

  /* Special case: scalar ^ scalar */
  if ((MNR (m1) == 1) && (MNC (m1) == 1) && (MNR (m2) == 1) && (MNC (m2) == 1))
  {
    m = mdr_elpower1 (m1, m2, type);
    return (m);
  }

  rerror ("matrix ^ matrix not supported yet");
  return (m);
}

void *
mdr_Power (MDR * M1, MDR * M2, int *type)
{

  void *m = 0;
  MDR *m1, *m2;
  /* Check sizes 1st. */

  // not an integer optimized
  if(M1->type == RLAB_TYPE_INT32)
    m1 = mdr_Float_BF (M1);
  else
    m1 = M1;
  if(M2->type == RLAB_TYPE_INT32)
    m2 = mdr_Float_BF (M2);
  else
    m2 = M2;

  if (MNR (m1) == 1 && MNC (m1) == 1)
  {
    /* Special case: 1-by-1 ^ Matrix */
    m = mdr_power1 (m1, m2, type);
  }
  else if (MNR (m2) == 1 && MNC (m2) == 1)
  {
    /* Special case: Matrix ^ 1-by-1 */
    m = mdr_power2 (m1, m2, type);
  }
  else if (MNR (m1) != MNR (m2) || MNC (m1) != MNC (m2))
  {
    rerror ("row and column dimensions must match for ^ operation");
  }
  else
  {
    /* Last possibility: Matrix ^ Matrix. */
    fprintf (stderr, "\tinvalid matrix dimensions\n");
    fprintf (stderr, "\tmatrix A: row = %i, column = %i\n", MNR (m1), MNC (m1));
    fprintf (stderr, "\tmatrix B: row = %i, column = %i\n", MNR (m2), MNC (m2));
    rerror ("matrix ^ matrix, invalid operation");
  }

  // not an integer optimized
  if(M1->type == RLAB_TYPE_INT32)
    mdr_Destroy(m1);
  if(M2->type == RLAB_TYPE_INT32)
    mdr_Destroy(m2);

  return (m);
}


/* **************************************************************
 * Element-by-element power operations...
 * ************************************************************** */

/*
 * Scalar .^ Matrix
 *   d    .^  m2
 *
 * Here we do the operation one of two ways:
 * etype 0: Positive operand (s), use a real power function.
 *          OR integer M
 * etype 1: Everything else (default).
 */

void *
mdr_elpower1 (MDR * m1, MDR * m2, int *type)
{
  void *m;
  double d;
  int i, size;
  int etype = 1;

  d = (double) MdrV0 (m1, 0);
  size = MNR (m2) * MNC (m2);

  /* We must do some checking 1st */
  if (d > 0.0)
  {
    /* Positive operand, use real power function */
    etype = 0;
  }
  else if (d == 0.0)
  {
    /* We already know the answer. */
    m = mdr_Create (MNR (m2), MNC (m2));
    mdr_Zero (m);
    *type = MATRIX_DENSE_REAL;

    /* Now check for 0 exponent (answer = 1). */
    for (i = 0; i < size; i++)
    {
      if (MdrV0 (m2, i) == 0.0)
        MdrV0 (m, i) = 1.0;
    }
    return (m);
  }
  else
  {
    /* Now check for integer powers */
    etype = 0;
    for (i = 0; i < size; i++)
    {
      if (floor (MdrV0 (m2, i)) != MdrV0 (m2, i))
      {
        etype = 1;
        break;
      }
    }
  }

  if (etype == 0)
  {
    m = mdr_Create (MNR (m2), MNC (m2));
    *type = MATRIX_DENSE_REAL;
    for (i = 0; i < size; i++)
    {
      MdrV0 (m, i) = errcheck (pow (d, MdrV0 (m2, i)), ".^");
    }
  }
  else
  {
    m = (MDC *) mdc_Create (MNR (m2), MNC (m2));
    *type = MATRIX_DENSE_COMPLEX;
    for (i = 0; i < size; i++)
    {
      MdcV0 ((MDC *) m, i) = cpow (d, MdrV0 (m2, i));
    }
  }
  return (m);
}

/*
 * Matrix .^ Scalar
 * m1 -- matrix dimensions
 * m2 -- scalar dimensions.
 */

void *
mdr_elpower2 (MDR * m1, MDR * m2, int *type)
{
  int i, size;
  double d;
  void *m;

  d = (double) MdrV0 (m2, 0);

  /* Zero power */
  if (d == 0.0)
  {
    m = mdr_Create (MNR (m1), MNC (m1));
    *type = MATRIX_DENSE_REAL;
    size = MNR (m1) * MNC (m1);
    for (i = 0; i < size; i++)
    {
      MdrV0 (m, i) = 1.0;
    }
  }

  /* integer power */
  else if (floor (d) == d)
  {
    m = mdr_Create (MNR (m1), MNC (m1));
    *type = MATRIX_DENSE_REAL;
    size = MNR (m1) * MNC (m1);
    for (i = 0; i < size; i++)
    {
      if (MdrV0 (m1, i) == 0.0)
	MdrV0 (m, i) = 0.0;
      else
	MdrV0 (m, i) = errcheck (pow (MdrV0 (m1, i), d), ".^");
    }
  }

  /* non-integer power, positive operand */
  else if (mdr_is_positive (m1))
  {
    m = mdr_Create (MNR (m1), MNC (m1));
    *type = MATRIX_DENSE_REAL;
    size = MNR (m1) * MNC (m1);
    for (i = 0; i < size; i++)
    {
      if (MdrV0 (m1, i) == 0.0)
	MdrV0 (m, i) = 0.0;
      else
	MdrV0 (m, i) = errcheck (pow (MdrV0 (m1, i), d), ".^");
    }
  }

  /* non-integer power, some negative operand elements */
  else
  {
    m = (MDC *) mdc_Create (MNR (m1), MNC (m1));
    *type = MATRIX_DENSE_COMPLEX;
    size = MNR (m1) * MNC (m1);
    for (i = 0; i < size; i++)
    {
      if (MdrV0 (m1, i) == 0.0)
      {
        MdcV0r ((MDC *) m, i) = 0.0;
        MdcV0i ((MDC *) m, i) = 0.0;
      }
      else
      {
        MdcV0 ((MDC *) m, i) = cpow(MdrV0 (m1, i), d);
      }
    }
  }

  return (m);
}

/*
 * Matrix .^ Matrix
 * Here we do the operation one of two ways:
 * type 0: Positive operand (M1), use a real power function.
 *         OR integer M2
 * type 1: Everything else (default).
 */

void *
mdr_elpower3 (MDR * m1, MDR * m2, int *type)
{
  int i, size;
  void *m;
  int etype = 0;

  if (MNR (m1) != MNR (m2) || MNC (m1) != MNC (m2))
  {
    fprintf (stderr, "\tmatrix dimensions must be consistent\n");
    fprintf (stderr, "\tmatrix A: row = %i, column = %i\n", MNR (m1), MNC (m1));
    fprintf (stderr, "\tmatrix B: row = %i, column = %i\n", MNR (m2), MNC (m2));
    rerror ("Element-power (.^) operation dimension mis-match");
  }

  size = MNR (m1) * MNC (m1);

  /* Do some operand checking */
  for (i = 0; i < size; i++)
  {
    if (floor (MdrV0 (m2, i)) != MdrV0 (m2, i))
    {
      etype = 1;
      break;
    }
  }

  if (etype == 0)
  {
    if (!mdr_is_positive (m1))
      etype = 1;
  }

  if (etype == 0)
  {
    m = mdr_Create (MNR (m1), MNC (m1));
    *type = MATRIX_DENSE_REAL;
    for (i = 0; i < size; i++)
    {
      MdrV0 (m, i) = errcheck (pow (MdrV0 (m1, i), MdrV0 (m2, i)), ".^");
    }
  }
  else
  {
    m = (MDC *) mdc_Create (MNR (m1), MNC (m1));
    *type = MATRIX_DENSE_COMPLEX;
    for (i = 0; i < size; i++)
    {
      if (MdrV0 (m1, i) == 0.0)
      {
        MdcV0r ((MDC *) m, i) = 0.0;
        MdcV0i ((MDC *) m, i) = 0.0;
      }
      else
      {
        MdcV0 ((MDC *) m, i) = cpow (MdrV0 (m1, i), MdrV0 (m2, i));
      }
    }
  }
  return (m);
}

void *
mdr_ElPower (MDR * M1, MDR * M2, int *type)
{
  void *m = 0;
  MDR *m1, *m2;

  // not an integer optimized
  if(M1->type == RLAB_TYPE_INT32)
    m1 = mdr_Float_BF (M1);
  else
    m1 = M1;
  if(M2->type == RLAB_TYPE_INT32)
    m2 = mdr_Float_BF (M2);
  else
    m2 = M2;

  /* Check sizes 1st. */

  /* Special case: 1-by-1 .^ Matrix */
  if (MNR (m1) == 1 && MNC (m1) == 1)
  {
    m = mdr_elpower1 (m1, m2, type);
  }

  /* Special case: Matrix .^ 1-by-1 */
  else if (MNR (m2) == 1 && MNC (m2) == 1)
  {
    m = mdr_elpower2 (m1, m2, type);
  }

  else if (MNR (m1) != MNR (m2) || MNC (m1) != MNC (m2))
  {
    rerror ("row and column dimensions must match for .^ operation");
  }
  else
  {
    /* Last possibility: Matrix .^ Matrix. */
    m = mdr_elpower3 (m1, m2, type);
  }

  // not an integer optimized
  if(M1->type == RLAB_TYPE_INT32)
    mdr_Destroy(m1);
  if(M2->type == RLAB_TYPE_INT32)
    mdr_Destroy(m2);

  return (m);
}

#include <gsl/gsl_math.h>
void *
mdr_ElPower_weight (MDR * m1, MDR *w1, MDR * m2, MDR *w2)
{
  MDR *new=0;
  double x=0, vx=0, a=0, va=0, vy=0;
  int i, j, i1, i2, j1, j2;
  int iw1=0, jw1=0, iw2=0, jw2=0, nr, nc;

  i1 = m1->nrow;
  j1 = m1->ncol;
  i2 = m2->nrow ;
  j2 = m2->ncol;
  nr = MAX(i1, i2);
  nc = MAX(j1, j2);

  if (w1)
  {
    iw1 = MNR(w1);
    jw1 = MNC(w1);
  }
  if (w2)
  {
    iw2 = MNR(w2);
    jw2 = MNC(w2);
  }

  new = mdr_Create(nr, nc);

  if (w1 || w2)
  {
    if (w1 && w2)
    {
      for(i=0;i<nr;i++)
        for(j=0;j<nc;j++)
      {
        // x
        x  = mdr0(m1, MIN(i,i1), MIN(j,j1));
        vx = 1 / mdr0(w1, MIN(i,iw1), MIN(j,jw1));
        // a
        a  = mdr0(m2, MIN(i,i2), MIN(j,j2));
        va = 1 / mdr0(w2, MIN(i,iw2), MIN(j,jw2));

        // this is x ^ a
        if (x > 0)
        {
          if (x != 1)
            vy = a * a * pow(x * x, a - 1) * vx + log(x) * log(x) * pow(x*x, a) * va;
          else
            vy = a * a * vx;
          Mdr0(new,i,j) = 1 / vy;
        }
        else if (x < 0)
        {
          vy = a * a * pow(x * x, a - 1) * vx +
              (log(ABS(x)) * log(ABS(x)) + M_PI * M_PI)  * pow(x*x, a) * va;
          Mdr0(new,i,j) = 1 / vy;
        }
        else
        {
          Mdr0(new,i,j) = 1.0 / 0.0;
        }
      }
    }
    else if (w1)
    {
      for(i=0;i<nr;i++)
        for(j=0;j<nc;j++)
      {
        // x
        x  = mdr0(m1, MIN(i,i1), MIN(j,j1));
        vx = 1 / mdr0(w1, MIN(i,iw1), MIN(j,jw1));
        // a
        a  = mdr0(m2, MIN(i,i2), MIN(j,j2));

        // this is x ^ a
        if (x > 0)
        {
          if (x != 1)
            vy = a * a * pow(x * x, a - 1) * vx;
          else
            vy = a * a * vx;
          Mdr0(new,i,j) = 1 / vy;
        }
        else if (x < 0)
        {
          vy = a * a * pow(x * x, a - 1) * vx;
          Mdr0(new,i,j) = 1 / vy;
        }
        else
        {
          Mdr0(new,i,j) = 1.0 / 0.0;
        }
      }
    }
    else if (w2)
    {
      for(i=0;i<nr;i++)
        for(j=0;j<nc;j++)
      {
        // x
        x  = mdr0(m1, MIN(i,i1), MIN(j,j1));
        // a
        a  = mdr0(m2, MIN(i,i2), MIN(j,j2));
        va = 1 / mdr0(w2, MIN(i,iw2), MIN(j,jw2));

        // this is x ^ a
        if (x > 0)
        {
          if (x != 1)
            vy = log(x) * log(x) * pow(x*x, a) * va;
          else
            vy = 0;
          Mdr0(new,i,j) = 1 / vy;
        }
        else if (x < 0)
        {
          vy = (log(ABS(x)) * log(ABS(x)) + M_PI * M_PI)  * pow(x*x, a) * va;
          Mdr0(new,i,j) = 1 / vy;
        }
        else
        {
          Mdr0(new,i,j) = 1.0 / 0.0;
        }
      }
    }
  }

  return (new);
}


/* **************************************************************
 * Negate a matrix.
 * ************************************************************** */

MDR *
mdr_Negate (MDR * m)
{
  int k, isize;
  MDR *mnew;

  isize = MNR (m) * MNC (m);

  if(m->type == RLAB_TYPE_INT32)
  {
    mnew = mdi_Create (MNR (m), MNC (m));
    for (k = 0; k < isize; k++)
    {
      MdiV0 (mnew, k) = -MdiV0 (m, k);
    }
  }
  else
  {
    mnew = mdr_Create (MNR (m), MNC (m));
    for (k = 0; k < isize; k++)
    {
      MdrV0 (mnew, k) = -MdrV0 (m, k);
    }
  }

  return (mnew);
}

/* **************************************************************
 * Return a logical-scalar for the input matrix.
 * ************************************************************** */

int
mdr_LogicalScalar (MDR * m)
{
  if ((MNR (m) == 1) && (MNC (m) == 1))
  {
    if(m->type == RLAB_TYPE_INT32)
    {
      if (MdiV0 (m, 0) == 0)
        return (0);
      else
        return (1);
    }
    else
    {
      if (MdrV0 (m, 0) == 0.0)
        return (0);
      else
        return (1);
    }
  }
  else if ((MNR (m) == 0) && (MNC (m) == 0))
  {
    return (0);
  }
  else
  {
    fprintf (stderr, "cannot compute a logical scalar value\n");
    fprintf (stderr, "from a %i by %i matrix\n", MNR (m), MNC (m));
    rerror ("cannot compute a logical scalar value");
  }

  return (0);			/* Shut up compiler. */
}

/* **************************************************************
 * Return a the size of the matrix.
 * ************************************************************** */

int mdr_Size (MDR * m)
{
  return (MNR (m) * MNC (m));
}

/* **************************************************************
 * Append two matrices together:  [ m1 , m2 ]
 * ************************************************************** */
MDR * mdr_Append (MDR * m1, MDR * m2)
{
  int i, j, nrow, ncol;
  MDR *new=0;

  /* Check for empty matrices... */
  if (SIZE(m1)<1 && SIZE(m2)<1)
  {
    new = mdr_Create(0,0);
    return new;
  }
  else if (SIZE(m1)<1)
  {
    new = mdr_Copy (m2);
    return (new);
  }
  else if (SIZE(m2)<1)
  {
    new = mdr_Copy (m1);
    return (new);
  }

  /* Create the new matrix, large enough to hold both. */
  int r1 = MNR(m1);
  int c1 = MNC(m1);
  int r2 = MNR(m2);
  int c2 = MNC(m2);
  nrow = MAX(r1,r2);
  ncol = c1 + c2;

  if(MD_TYPE_INT32(m1) && MD_TYPE_INT32(m2))
  {
    new = mdi_Create (nrow, ncol);
    for (i=0; i<nrow; i++)
    {
      for (j=0; j<c1; j++)
        Mdi0(new,i,j)     = Mdi0(m1,MIN(i,r1-1),j);
      for (j=0; j<c2; j++)
        Mdi0(new,i,j+c1)  = Mdi0(m2,MIN(i,r2-1),j);
    }
  }
  else
  {
    new = mdr_Create (nrow, ncol);
    for (i=0; i<nrow; i++)
    {
      for (j=0; j<c1; j++)
        Mdr0(new,i,j)     = mdr0(m1,MIN(i,r1-1),j);
      for (j=0; j<c2; j++)
        Mdr0(new,i,j+c1)  = mdr0(m2,MIN(i,r2-1),j);
    }
  }
  return (new);
}


/* **************************************************************
 * Stack two matrices together:  [ m1 ; m2 ]
 * ************************************************************** */

MDR *
mdr_Stack (MDR * m1, MDR * m2)
{
  int i, j, nrow, ncol;
  MDR *new=0;

  /* Check for empty matrices... */
  if (SIZE(m1)<1 && SIZE(m2)<1)
  {
    new = mdr_Create(0,0);
    return new;
  }
  else if (SIZE(m1)<1)
  {
    new = mdr_Copy (m2);
    return (new);
  }
  else if (SIZE(m2)<1)
  {
    new = mdr_Copy (m1);
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

  if(MD_TYPE_INT32(m1) && MD_TYPE_INT32(m2))
  {
    new = mdi_Create (nrow, ncol);
    for (j=0; j<ncol; j++)
    {
      for (i=0; i<r1; i++)
        Mdi0(new,i,j)     = Mdi0(m1,i,MIN(j,c1-1));
      for (i=0; i<r2; i++)
        Mdi0(new,i+r1,j)  = Mdi0(m2,i,MIN(j,c2-1));
    }
  }
  else
  {
    new = mdr_Create (nrow, ncol);
    for (j=0; j<ncol; j++)
    {
      for (i=0; i<r1; i++)
        Mdr0(new,i,j)     = mdr0(m1,i,MIN(j,c1-1));
      for (i=0; i<r2; i++)
        Mdr0(new,i+r1,j)  = mdr0(m2,i,MIN(j,c2-1));
    }
  }
  return (new);
}

/* **************************************************************
 * Sub-Matrix Expression Evaluation.
 * This function must look at the index expressions, and decide
 * how to act. At this point, MATRIX_DENSE_REAL only understands
 * integer-like indices. Any non-numeric index values must be
 * handled by a coercion function.
 * ************************************************************** */

MDR *
mdr_MatrixSub (MDR * var, int *i, int *j, int *type)
{
  int m, n;
  MDR *new;

  *type = MATRIX_DENSE_REAL;

  /* Handle empty matrix indices. */
  if ((i == 0) || (j == 0) || (var->nrow * var->ncol == 0))
  {
    new = mdr_Create (0, 0);
    return (new);
  }

  if(var->type == RLAB_TYPE_INT32)
  {
    new = mdi_Create (i[0], j[0]);
    for (m = 1; m <= i[0]; m++)
    {
      if (i[m] > MNR (var))
      {
        mdr_Destroy (new);
        fprintf (stderr, "index exceeds matrix limits\n");
        fprintf (stderr, "index value: %i, matrix size: %i\n", i[m], MNR (var));
        rerror ("sub-matrix evaluation");
      }
      for (n = 1; n <= j[0]; n++)
      {
        if (j[n] > MNC (var))
        {
          mdr_Destroy (new);
          fprintf (stderr, "index exceeds matrix limits\n");
          fprintf (stderr, "index value: %i, matrix size: %i\n", j[n], MNC (var));
          rerror ("sub-matrix evaluation");
        }
        Mdi1 (new, m, n) = Mdi1 (var, i[m], j[n]);
      }
    }
  }
  else
  {
    new = mdr_Create (i[0], j[0]);
    for (m = 1; m <= i[0]; m++)
    {
      if (i[m] > MNR (var))
      {
        mdr_Destroy (new);
        fprintf (stderr, "index exceeds matrix limits\n");
        fprintf (stderr, "index value: %i, matrix size: %i\n", i[m], MNR (var));
        rerror ("sub-matrix evaluation");
      }
      for (n = 1; n <= j[0]; n++)
      {
        if (j[n] > MNC (var))
        {
          mdr_Destroy (new);
          fprintf (stderr, "index exceeds matrix limits\n");
          fprintf (stderr, "index value: %i, matrix size: %i\n", j[n], MNC (var));
          rerror ("sub-matrix evaluation");
        }
        Mdr1 (new, m, n) = Mdr1 (var, i[m], j[n]);
      }
    }
  }
  return (new);
}

/* **************************************************************
 * Sub-Matrix Expression Evaluation-Row.
 * Extraxt the sub-matrix specified by row indices only.
 * ************************************************************** */

MDR *
mdr_MatrixSubR (MDR * var, int *i, int *type)
{
  int m, n;
  MDR *new;

  *type = MATRIX_DENSE_REAL;

  /* Handle empty matrix indices. */

  if ((i == 0) || (var->nrow * var->ncol == 0))
  {
    new = mdr_Create (0, 0);
    return (new);
  }

  if(var->type == RLAB_TYPE_INT32)
  {
    new = mdi_Create (i[0], MNC (var));
    for (m = 1; m <= i[0]; m++)
    {
      if (i[m] > MNR (var))
      {
        mdr_Destroy (new);
        fprintf (stderr, "index exceeds matrix limits\n");
        fprintf (stderr, "index value: %i, matrix size: %i\n", i[m], MNR (var));
        rerror ("sub-matrix evaluation");
      }
      for (n = 1; n <= MNC (var); n++)
      {
        Mdi1 (new, m, n) = Mdi1 (var, i[m], n);
      }
    }
  }
  else
  {
    new = mdr_Create (i[0], MNC (var));
    for (m = 1; m <= i[0]; m++)
    {
      if (i[m] > MNR (var))
      {
        mdr_Destroy (new);
        fprintf (stderr, "index exceeds matrix limits\n");
        fprintf (stderr, "index value: %i, matrix size: %i\n", i[m], MNR (var));
        rerror ("sub-matrix evaluation");
      }
      for (n = 1; n <= MNC (var); n++)
      {
        Mdr1 (new, m, n) = Mdr1 (var, i[m], n);
      }
    }
  }

  return (new);
}

/* **************************************************************
 * Sub-Matrix Expression Evaluation-Column.
 * Extraxt the sub-matrix specified by column indices only.
 * ************************************************************** */

MDR *
mdr_MatrixSubC (MDR * var, int *j, int *type)
{
  int m, n;
  MDR *new;

  *type = MATRIX_DENSE_REAL;

  /* Handle empty matrix indices. */

  if ((j == 0) || (var->nrow * var->ncol == 0))
  {
    new = mdr_Create (0, 0);
    return (new);
  }

  if (var->type == RLAB_TYPE_INT32)
  {
    new = mdi_Create (MNR (var), j[0]);
    for (m = 1; m <= MNR (var); m++)
    {
      for (n = 1; n <= j[0]; n++)
      {
        if (j[n] > MNC (var))
        {
          mdr_Destroy (new);
          fprintf (stderr, "index exceeds matrix limits\n");
          fprintf (stderr, "index value: %i, matrix size: %i\n", j[n], MNC (var));
          rerror ("sub-matrix evaluation");
        }
        Mdi1 (new, m, n) = Mdi1 (var, m, j[n]);
      }
    }
  }
  else
  {
    new = mdr_Create (MNR (var), j[0]);
    for (m = 1; m <= MNR (var); m++)
    {
      for (n = 1; n <= j[0]; n++)
      {
        if (j[n] > MNC (var))
        {
          mdr_Destroy (new);
          fprintf (stderr, "index exceeds matrix limits\n");
          fprintf (stderr, "index value: %i, matrix size: %i\n", j[n], MNC (var));
          rerror ("sub-matrix evaluation");
        }
        Mdr1 (new, m, n) = Mdr1 (var, m, j[n]);
      }
    }
  }
  return (new);
}

/* **************************************************************
 * Simple cover package to mdr_MatrixSubC. Column index starts
 * with 1.
 * ************************************************************** */

MDR *
mdr_PartitionCol (MDR * m, int col)
{
  int index[2];
  int type;

  index[0] = 1;
  index[1] = col;

  return (mdr_MatrixSubC (m, index, &type));
}

/* **************************************************************
 * Assign to a range of a matrix. We do not create a new matrix.
 * Automatically extend the size of a matrix if I or J indices
 * exceed current bounds.
 * ************************************************************** */

MDR *
mdr_MatrixAssign (MDR * var, int *i, int *j, MDR * rhs)
{
  int m, n;
  int dflag;
  MDR *mtmp=0;

  if ((i == 0) || (j == 0))
  {
    if (SIZE(rhs) != 0)
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

  /* Check LHS for UNDEF. */
  if (!var)
  {
    if(rhs->type == RLAB_TYPE_INT32)
      var = mdi_Create (1, 1);
    else
      var = mdr_Create (1, 1);
    mdr_Zero (var);
  }

  /* Check whether RHS and LHS are same real/integer. */
  if (var->type != rhs->type)
  {
    fprintf (stderr, "cannot assign integer to real matrix and vice versa\n");
    rerror ("matrix-assign");
  }

  /* Check RHS for empty matrix. */
  if (SIZE(rhs) == 0)
  {
    fprintf (stderr, "cannot assign empty matrix to element(s) of a matrix\n");
    rerror ("matrix-assign");
  }

  /* Check the LHS, and RHS dimensions. */
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

  /* Now do the assign. */

  if (var == rhs)
  {
    mtmp = mdr_Copy (rhs);
    dflag = 1;
  }
  else
  {
    mtmp = rhs;
    dflag = 0;
  }

  /* Check for empty matrix. */
  if (SIZE(var) == 0)
    mdr_Extend (var, 1, 1);

  for (m = 1; m <= i[0]; m++)
  {
    if (i[m] > MNR (var))
    {
      mdr_Extend (var, i[m], MNC (var));
    }
    if(var->type == RLAB_TYPE_INT32)
    {
      for (n = 1; n <= j[0]; n++)
      {
        if (j[n] > MNC (var))
        {
          mdr_Extend (var, MNR (var), j[n]);
        }
        Mdi1 (var, i[m], j[n]) = mdi1 (mtmp, m, n);
      }
    }
    else
    {
      for (n = 1; n <= j[0]; n++)
      {
        if (j[n] > MNC (var))
        {
          mdr_Extend (var, MNR (var), j[n]);
        }
        Mdr1 (var, i[m], j[n]) = mdr1 (mtmp, m, n);
      }
    }
  }

  if (dflag)
    mdr_Destroy (mtmp);

  return (var);
}

/* **************************************************************
 * Assign to the rows of a matrix.
 * a[ i ; ] = b
 * ************************************************************** */

MDR *
mdr_MatrixAssignR (MDR * var, int *i, MDR * rhs)
{
  int m, n;
  int dflag;
  MDR *mtmp;

  if (i == 0)
  {
    if (MdrNR (rhs) * MdrNC (rhs) != 0)
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

  /* Check LHS for UNDEF. */
  if (var == 0)
  {
    if(rhs->type == RLAB_TYPE_INT32)
      var = mdi_Create (1, 1);
    else
      var = mdr_Create (1, 1);
    mdr_Zero (var);
  }

  /*
   * Check RHS for empty matrix.
   * If RHS is empty, then eliminate the
   * rows specified in i.
   */

  if (MNR (rhs) * MNC (rhs) == 0)
  {
    int copy, j, k, n, *rs1, rs1_size, *rs2, rs2_size;

    /* First, get rid of duplicates in the row-spec. */

    rs1 = remove_duplicates (i + 1, i[0], &rs1_size);

    /* Throw out row-spec elements greater then MNR. */

    rs2 = remove_invalid (rs1, rs1_size, MNR (var), &rs2_size);
    GC_FREE (rs1);

    /* Now, create the new matrix with the right size. */
    if(var->type == RLAB_TYPE_INT32)
      mtmp = mdi_Create (MNR (var) - rs2_size, MNC (var));
    else
      mtmp = mdr_Create (MNR (var) - rs2_size, MNC (var));

    /* We might as well exit now, if mtmp = [] */
    if (MNR (mtmp) == 0)
    {
      GC_FREE (rs2);
      mdr_Destroy (var);
      var = mtmp;
      return (var);
    }

    /*
     * Now, copy the correct elements into the smaller matrix.
     * This must be done element at a time, cause the matrix
     * is stored columnwise.
     */
    n = 0;
    if(var->type == RLAB_TYPE_INT32)
    {
      for (j = 0; j < MNR (var); j++)
      {
        copy = 1;
        for (k = 0; k < rs2_size; k++)
        {
          if (rs2[k] == (j + 1))
          {
            copy = 0;
            break;
          }
        }
        if (copy)
        {
          for (k = 0; k < MNC (var); k++)
            Mdi0 (mtmp, n, k) = Mdi0 (var, j, k);
          n++;
        }
      }
    }
    else
    {
      for (j = 0; j < MNR (var); j++)
      {
        copy = 1;
        for (k = 0; k < rs2_size; k++)
        {
          if (rs2[k] == (j + 1))
          {
            copy = 0;
            break;
          }
        }
        if (copy)
        {
          for (k = 0; k < MNC (var); k++)
            Mdr0 (mtmp, n, k) = Mdr0 (var, j, k);
          n++;
        }
      }
    }

    /* Free the old, and replace it with the new. */
    GC_FREE (rs2);
    mdr_Destroy (var);
    var = mtmp;
    return (var);
  }

//   if((rhs->type == RLAB_TYPE_INT32) != (var->type == RLAB_TYPE_INT32))
//   {
//     rerror
//         ("matrix-assign: cannot assign real to integer matrix and vice versa.");
//   }

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
    mdr_Extend (var, 1, 1);
  }

  /* Now do the assign. */

  if (var == rhs)
  {
    mtmp = mdr_Copy (rhs);
    dflag = 1;
  }
  else
  {
    mtmp = rhs;
    dflag = 0;
  }

  if(var->type == RLAB_TYPE_INT32)
  {
    for (m = 1; m <= i[0]; m++)
    {
      if (i[m] > MNR (var))
      {
        mdr_Extend (var, i[m], MNC (var));
      }
      for (n = 1; n <= MNC (rhs); n++)
      {
        if (n > MNC (var))
        {
          mdr_Extend (var, MNR (var), n);
        }
        Mdi1 (var, i[m], n) = mdi1 (mtmp, m, n);
      }
    }
  }
  else
  {
    for (m = 1; m <= i[0]; m++)
    {
      if (i[m] > MNR (var))
      {
        mdr_Extend (var, i[m], MNC (var));
      }
      for (n = 1; n <= MNC (rhs); n++)
      {
        if (n > MNC (var))
        {
          mdr_Extend (var, MNR (var), n);
        }
        Mdr1 (var, i[m], n) = mdr1 (mtmp, m, n);
      }
    }
  }

  if (dflag)
    mdr_Destroy (mtmp);

  return (var);
}

MDR *
mdr_MatrixAssignC (MDR * var, int *j, MDR * rhs)
{
  int m, n;
  int dflag;
  MDR *mtmp;

  if (j == 0)
  {
    if (MdrNR (rhs) * MdrNC (rhs) != 0)
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

  /* Check LHS for UNDEF. */
  if (var == 0)
  {
    if(rhs->type == RLAB_TYPE_INT32)
      var = mdi_Create (1, 1);
    else
      var = mdr_Create (1, 1);
    mdr_Zero (var);
  }

  /*
   * Check RHS for empty matrix.
   * If RHS is empty, then eliminate the
   * rows specified in i.
   */

  if (MNR (rhs) * MNC (rhs) == 0)
  {
    int copy, jj, k, n, *rs1, rs1_size, *rs2, rs2_size;

    /* First, get rid of duplicates in the row-spec. */

    rs1 = remove_duplicates (j + 1, j[0], &rs1_size);

    /* Throw out row-spec elements greater then MNR. */

    rs2 = remove_invalid (rs1, rs1_size, MNC (var), &rs2_size);
    GC_FREE (rs1);

    /* Now, create the new matrix with the right size. */
    if(var->type == RLAB_TYPE_INT32)
      mtmp = mdi_Create (MNR (var), MNC (var) - rs2_size);
    else
      mtmp = mdr_Create (MNR (var), MNC (var) - rs2_size);

    /* We might as well exit now, if mtmp = [] */
    if (MNR (mtmp) == 0)
    {
      GC_FREE (rs2);
      mdr_Destroy (var);
      var = mtmp;
      return (var);
    }

    /*
     * Now, copy the correct elements into the smaller matrix.
     */
    if(var->type == RLAB_TYPE_INT32)
    {
      n = 0;
      for (jj = 0; jj < MNC (var); jj++)
      {
        copy = 1;
        for (k = 0; k < rs2_size; k++)
        {
          if (rs2[k] == (jj + 1))
          {
            copy = 0;
            break;
          }
        }
        if (copy)
        {
          for (k = 0; k < MNR (var); k++)
            Mdi0 (mtmp, k, n) = Mdi0 (var, k, jj);
          n++;
        }
      }
    }
    else
    {
      n = 0;
      for (jj = 0; jj < MNC (var); jj++)
      {
        copy = 1;
        for (k = 0; k < rs2_size; k++)
        {
          if (rs2[k] == (jj + 1))
          {
            copy = 0;
            break;
          }
        }
        if (copy)
        {
          for (k = 0; k < MNR (var); k++)
            Mdr0 (mtmp, k, n) = Mdr0 (var, k, jj);
          n++;
        }
      }
    }
    /* Free the old, and replace it with the new. */
    GC_FREE (rs2);
    mdr_Destroy (var);
    var = mtmp;
    return (var);
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
    mdr_Extend (var, 1, 1);
  }

  /* Now do the assign. */

  if (var == rhs)
  {
    mtmp = mdr_Copy (rhs);
    dflag = 1;
  }
  else
  {
    mtmp = rhs;
    dflag = 0;
  }

  if(var->type == RLAB_TYPE_INT32)
  {
    for (m = 1; m <= MNR (rhs); m++)
    {
      if (m > MNR (var))
      {
        mdr_Extend (var, m, MNC (var));
      }
      for (n = 1; n <= j[0]; n++)
      {
        if (j[n] > MNC (var))
        {
          mdr_Extend (var, MNR (var), j[n]);
        }

        if (mtmp->type == RLAB_TYPE_INT32)
          Mdi1 (var, m, j[n]) = Mdi1 (mtmp, m, n);
        else
          Mdi1 (var, m, j[n]) = Mdr1 (mtmp, m, n);
      }
    }
  }
  else
  {
    for (m = 1; m <= MNR (rhs); m++)
    {
      if (m > MNR (var))
      {
        mdr_Extend (var, m, MNC (var));
      }
      for (n = 1; n <= j[0]; n++)
      {
        if (j[n] > MNC (var))
        {
          mdr_Extend (var, MNR (var), j[n]);
        }

        if (mtmp->type == RLAB_TYPE_INT32)
          Mdr1 (var, m, j[n]) = Mdi1 (mtmp, m, n);
        else
          Mdr1 (var, m, j[n]) = Mdr1 (mtmp, m, n);
      }
    }
  }

  if (dflag)
    mdr_Destroy (mtmp);

  return (var);
}

/* **************************************************************
 * Vector Sub-Expression
 * ************************************************************** */

MDR *
mdr_VectorSub (MDR * m, int *i, int *type)
{
  int j, size, msize;
  MDR *new;

  *type = MATRIX_DENSE_REAL;

  /* Handle empty matrix indices. */

  if (i == 0)
  {
    new = mdr_Create (0, 0);
    return (new);
  }

  size = i[0];
  msize = MNR (m) * MNC (m);
  if (!msize)
  {
    new = mdr_Create (0, 0);
    return (new);
  }

  if(m->type == RLAB_TYPE_INT32)
  {
    if (MNC (m) == 1)
      new = mdi_Create (size, 1);
    else
      new = mdi_Create (1, size);
    for (j = 1; j <= size; j++)
    {
      if (i[j] > msize)
      {
        mdr_Destroy (new);
        fprintf (stderr, "\tindex exceeds matrix limits\n");
        fprintf (stderr, "\tindex value: %i, matrix size: %i\n", i[j], msize);
        rerror ("sub-matrix evaluation");
      }
      MdiV1 (new, j) = MdiV1 (m, i[j]);
    }

  }
  else
  {
    if (MNC (m) == 1)
      new = mdr_Create (size, 1);
    else
      new = mdr_Create (1, size);
    for (j = 1; j <= size; j++)
    {
      if (i[j] > msize)
      {
        mdr_Destroy (new);
        fprintf (stderr, "\tindex exceeds matrix limits\n");
        fprintf (stderr, "\tindex value: %i, matrix size: %i\n", i[j], msize);
        rerror ("sub-matrix evaluation");
      }
      MdrV1 (new, j) = MdrV1 (m, i[j]);
    }
  }

  return (new);
}

/* **************************************************************
 * Assign to a matrix, like it was a vector...
 * var[j] = rhs
 * ************************************************************** */

MDR *
mdr_VectorAssign (MDR * var, int *j, MDR * rhs)
{
  int dflag, i, size;
  int nr, nc;
  MDR *mtmp;

  /* Check for empty indices. */
  if (j == 0)
  {
    if (MdrNR (rhs) * MdrNC (rhs) != 0)
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

  /* Check LHS for UNDEF. */
  if (var == 0)
  {
    if(rhs->type == RLAB_TYPE_INT32)
      var = mdi_Create (1, 1);
    else
      var = mdr_Create (1, 1);
    mdr_Zero (var);
  }

  /*
   * Check RHS for empty matrix.
   * If RHS is empty, then eliminate the
   * rows specified in i.
   */

  if (MNR (rhs) * MNC (rhs) == 0)
  {
    int copy, jj, k, n, *rs1, rs1_size, *rs2, rs2_size;

    /* First, check that the LHS is a row or column matrix. */

    if ((MNR (var) != 1) && (MNC (var) != 1))
      rerror ("vector-assign: LHS must be a row or column when RHS is []");

    /* Get rid of duplicates in the row-spec. */

    rs1 = remove_duplicates (j + 1, j[0], &rs1_size);

    /* Throw out row-spec elements greater then MNR. */

    rs2 = remove_invalid (rs1, rs1_size, (MNR (var) * MNC (var)), &rs2_size);
    GC_FREE (rs1);

    /* Now figure out which way to make the new matrix. */

    if(var->type == RLAB_TYPE_INT32)
    {
      if (MNR (var) == 1)
        mtmp = mdi_Create (MNR (var), MNC (var) - rs2_size);
      else
        mtmp = mdi_Create (MNR (var) - rs2_size, MNC (var));

      /* We might as well exit now, if mtmp = [] */
      if (MNR (mtmp) == 0)
      {
        GC_FREE (rs2);
        mdr_Destroy (var);
        var = mtmp;
        return (var);
      }
      /*
      * Now, copy the correct elements into the smaller matrix.
      * This must be done element at a time, cause the matrix
      * is stored columnwise.
      */
      n = 0;
      for (jj = 0; jj < (MNR (var) * MNC (var)); jj++)
      {
        copy = 1;
        for (k = 0; k < rs2_size; k++)
        {
          if (rs2[k] == (jj + 1))
          {
            copy = 0;
            break;
          }
        }
        if (copy)
          MdiV0 (mtmp, n++) = MdiV0 (var, jj);
      }
    }
    else
    {
      if (MNR (var) == 1)
        mtmp = mdr_Create (MNR (var), MNC (var) - rs2_size);
      else
        mtmp = mdr_Create (MNR (var) - rs2_size, MNC (var));

      /* We might as well exit now, if mtmp = [] */
      if (MNR (mtmp) == 0)
      {
        GC_FREE (rs2);
        mdr_Destroy (var);
        var = mtmp;
        return (var);
      }
      /*
       * Now, copy the correct elements into the smaller matrix.
       * This must be done element at a time, cause the matrix
       * is stored columnwise.
       */
      n = 0;
      for (jj = 0; jj < (MNR (var) * MNC (var)); jj++)
      {
        copy = 1;
        for (k = 0; k < rs2_size; k++)
        {
          if (rs2[k] == (jj + 1))
          {
            copy = 0;
            break;
          }
        }
        if (copy)
          MdrV0 (mtmp, n++) = MdrV0 (var, jj);
      }
    }

    /* Free the old, and replace it with the new. */
    GC_FREE (rs2);
    mdr_Destroy (var);
    var = mtmp;
    return (var);
  }

  /* Check the dimensions. */
  if (j[0] != (MNR (rhs) * MNC (rhs)))
  {
    fprintf (stderr, "\tLHS, and RHS column sizes must match\n");
    fprintf (stderr, "\tLHS: %i, RHS: %i\n", j[0], MNR (rhs) * MNC (rhs));
    rerror ("vector-assign");
  }

  /* Now do the assign. */

  /*
   * Check to make sure the rhs is not the same object as the lhs.
   * If it is, then make a temporary copy.
   */

  if (var == rhs)
  {
    mtmp = mdr_Copy (rhs);
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
    mdr_Extend (var, 1, 1);

  if(rhs->type == RLAB_TYPE_INT32)
  {
    for (i = size; i >= 1; i--)
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
        mdr_Extend (var, nr, nc);
      }
      if (var->type == RLAB_TYPE_INT32)
        MdiV1 (var, (j[i])) = MdiV1 (mtmp, i);
      else
        MdrV1 (var, (j[i])) = MdiV1 (mtmp, i);
    }
  }
  else
  {
    for (i = size; i >= 1; i--)
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
        mdr_Extend (var, nr, nc);
      }
      if (var->type == RLAB_TYPE_INT32)
        MdiV1 (var, (j[i])) = MdrV1 (mtmp, i);
      else
        MdrV1 (var, (j[i])) = MdrV1 (mtmp, i);
    }
  }

  if (dflag)
  {
    mdr_Destroy (mtmp);
  }

  return (var);
}

/* **************************************************************
 * Coerce matrix indices into ints for use later...
 * Return an int array, the 1st element is the number of elements.
 * ************************************************************** */

int *
mdr_IndCoerceInt (MDR * m, MDR * i)
{
  int *ind, k, size;

  if ((size = MNR (i) * MNC (i)) == 0)
  {
    return (0);
  }

  ind = (int *) GC_MAIOP ((size + 1) * sizeof (int));
  if (ind == 0)
    rerror ("out of memory");
  ind[0] = size;

  if(i->type == RLAB_TYPE_INT32)
  {
    for (k = 1; k <= size; k++)
    {
      ind[k] = MdiV1 (i, k);
      if (ind[k] <= 0)
      {
        fprintf (stderr, "matrix index <= 0 not allowed\n");
        fprintf (stderr, "index value: %i\n", ind[k]);
        GC_FREE (ind);
        rerror ("matrix index coerce");
      }
    }
  }
  else
  {
    for (k = 1; k <= size; k++)
    {
      ind[k] = (int) (MdrV1 (i, k));
      if (ind[k] <= 0)
      {
        fprintf (stderr, "matrix index <= 0 not allowed\n");
        fprintf (stderr, "index value: %i\n", ind[k]);
        GC_FREE (ind);
        rerror ("matrix index coerce");
      }
    }
  }

  return (ind);
}

/* **************************************************************
 * Return the i'th value needed during for-loop execution.
 * Do indexing from 1 for [m].
 * ************************************************************** */

MDR *
mdr_ForLoopValue (MDR * m, int i)
{
  MDR *new;
  if(m->type == RLAB_TYPE_INT32)
  {
    new = mdi_Create (1, 1);
    MdiV0 (new, 0) = MdiV1 (m, i);
  }
  else
  {
    new = mdr_Create (1, 1);
    MdrV0 (new, 0) = MdrV1 (m, i);
  }

  return (new);
}

/* **************************************************************
 * Return the class of a matrix-dense-real.
 * ************************************************************** */

char *
mdr_Class (MDR * m)
{
  return (cpstr ("num"));
}

/* **************************************************************
 * Return matrix-dense-real member references.
 * ************************************************************** */

void *
mdr_MemberRef (MDR * m, char *name, int *type)
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
    ent_data (ne) = mds_CreateScalar ( "num" );
    ent_SetType (ne, MATRIX_DENSE_STRING);
    rptr = (void *) ne;
    *type = ENTITY;
  }
  else if (!strcmp (name, RLAB_MEMBER_TYPE))
  {
    ne = ent_Create ();
    if(m->type == RLAB_TYPE_INT32)
      ent_data (ne) = mds_CreateScalar ( RLAB_MEMBER_TYPE_INT32 );
    else
      ent_data (ne) = mds_CreateScalar ( RLAB_MEMBER_TYPE_REAL );
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
mdr_Members (MDR * m, int *n)
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
 * Compare two matrices (m1 == m2)
 * ************************************************************** */

MDR *
mdr_Eq (MDR * m1, MDR * m2)
{
  int i,j,nr,nc;
  MDR *new;

  // special case: both matrices are zero-size:
  //    one is too lazy to use isempty function
  if (EQNULL(m1) && EQNULL(m2))
  {
    new = mdr_CreateScalar(1.0);
    return new;
  }
  else if ((!EQNULL(m1) && EQNULL(m2)) || (!EQNULL(m2) && EQNULL(m1)))
  {
    new = mdr_CreateScalar(0.0);
    return new;
  }

  nr = MAX(MNR (m1), MNR (m2));
  nc = MAX(MNC (m1), MNC (m2));

  if (MD_TYPE_INT32(m1) && MD_TYPE_INT32(m2))
  {
    new = mdi_Create (nr, nc);
    for (i = 0; i < nr; i++) for (j = 0; j < nc; j++)
        Mdi0 (new,i,j) = (int)
            (Mdi0 (m1, MIN(i,MNR(m1)-1), MIN(j,MNC(m1)-1)) ==  Mdi0 (m2, MIN(i,MNR(m2)-1), MIN(j,MNC(m2)-1)));
  }
  else
  {
    new = mdr_Create (nr, nc);
    for (i = 0; i < nr; i++) for (j = 0; j < nc; j++)
        Mdr0 (new,i,j) = (double)
            (mdr0 (m1, MIN(i,MNR(m1)-1), MIN(j,MNC(m1)-1)) ==  mdr0 (m2, MIN(i,MNR(m2)-1), MIN(j,MNC(m2)-1)));
  }

  return (new);
}

MDR *
mdr_FindVector (MDR * m1, MDR * m2)
{
  int i,j,n, size1, size2;
  MDR *new;

  // no empty matrices
  if ((MNR (m1)*MNC (m1)==0) || (MNR (m2)*MNC (m2)==0))
  {
    new = mdr_Create(0,0);
    return new;
  }
  // no full matrices
  if ((MNR (m1)!=1 && MNC (m1)!=1) || (MNR (m2)!=1 && MNC (m2)!=1))
  {
    new = mdr_Create(0,0);
    return new;
  }

  size1 = MNR (m1) * MNC (m1);
  size2 = MNR (m2) * MNC (m2);

  // first vector has to be greater in size than the second
  if (size1 < size2)
  {
    new = mdr_Create(0,0);
    return new;
  }

  // first time just count
  i = 0;
  j = 0;
  n = 0;
  while (i<size1)
  {
    if (mdrV0(m1,i) == mdrV0(m2,j))
    {
      j++;
      if (j==size2)
      {
        n++;
        j=0;
      }
    }
    else
      j = 0;

    i++;
  }

  // now create it
  if (n > 0)
  {
    if (MNR(m1) == 1)
      new = mdr_Create(1,n);
    else
      new = mdr_Create(n,1);

    i = 0;
    j = 0;
    n = 0;
    while (i<size1)
    {
      if (mdrV0(m1,i) == mdrV0(m2,j))
      {
        j++;
        if (j==size2)
        {
          MdrV0(new,n) = i-size2+2;
          n++;
          j=0;
        }
      }
      else
        j = 0;

      i++;
    }
  }
  else
    new = mdr_Create(0,0);

  return (new);
}


MDR *
mdr_Ne (MDR * m1, MDR * m2)
{
  int i,j,nr,nc;
  MDR *new;

  // special case: both matrices are zero-size:
  //    one is too lazy to use isempty function
  if (EQNULL(m1) && EQNULL(m2))
  {
    new = mdr_CreateScalar(0.0);
    return new;
  }
  else if ((!EQNULL(m1) && EQNULL(m2)) || (!EQNULL(m2) && EQNULL(m1)))
  {
    new = mdr_CreateScalar(1.0);
    return new;
  }

  nr = MAX(MNR (m1), MNR (m2));
  nc = MAX(MNC (m1), MNC (m2));

  new = mdr_Create (nr, nc);
  for (i = 0; i < nr; i++) for (j = 0; j < nc; j++)
      Mdr0 (new,i,j) = (double)
          (mdr0 (m1, MIN(i,MNR(m1)-1), MIN(j,MNC(m1)-1)) !=  mdr0 (m2, MIN(i,MNR(m2)-1), MIN(j,MNC(m2)-1)));
  return (new);
}

MDR *
mdr_Lt (MDR * m1, MDR * m2)
{
  int i,j,nr,nc;
  MDR *new;

  if (SIZE(m2) < 1)
  {
    new = mdr_Create(0,0);
    return new;
  }
  if (SIZE(m1) < 1)
  {
    new = mdr_Copy(m2);
    return new;
  }

  nr = MAX(MNR (m1), MNR (m2));
  nc = MAX(MNC (m1), MNC (m2));

  new = mdr_Create (nr, nc);
  for (i = 0; i < nr; i++) for (j = 0; j < nc; j++)
      Mdr0 (new,i,j) = (double)
          (mdr0 (m1, MIN(i,MNR(m1)-1), MIN(j,MNC(m1)-1)) <  mdr0 (m2, MIN(i,MNR(m2)-1), MIN(j,MNC(m2)-1)));
  return (new);
}

MDR *
mdr_Le (MDR * m1, MDR * m2)
{
  int i,j,nr,nc;
  MDR *new;

  if (SIZE(m2) < 1)
  {
    new = mdr_Create(0,0);
    return new;
  }
  if (SIZE(m1) < 1)
  {
    new = mdr_Copy(m2);
    return new;
  }

  nr = MAX(MNR (m1), MNR (m2));
  nc = MAX(MNC (m1), MNC (m2));

  new = mdr_Create (nr, nc);
  for (i = 0; i < nr; i++) for (j = 0; j < nc; j++)
      Mdr0 (new,i,j) = (double)
          (mdr0 (m1, MIN(i,MNR(m1)-1), MIN(j,MNC(m1)-1)) <=  mdr0 (m2, MIN(i,MNR(m2)-1), MIN(j,MNC(m2)-1)));
  return (new);
}

MDR *
mdr_Gt (MDR * m1, MDR * m2)
{
  int i,j,nr,nc;
  MDR *new;

  if (SIZE(m1) < 1)
  {
    new = mdr_Create(0,0);
    return new;
  }
  if (SIZE(m2) < 1)
  {
    new = mdr_Copy(m1);
    return new;
  }

  nr = MAX(MNR (m1), MNR (m2));
  nc = MAX(MNC (m1), MNC (m2));

  new = mdr_Create (nr, nc);
  for (i = 0; i < nr; i++) for (j = 0; j < nc; j++)
      Mdr0 (new,i,j) = (double)
          (mdr0 (m1, MIN(i,MNR(m1)-1), MIN(j,MNC(m1)-1)) > mdr0 (m2, MIN(i,MNR(m2)-1), MIN(j,MNC(m2)-1)));
  return (new);
}

MDR *
mdr_Ge (MDR * m1, MDR * m2)
{
  int i,j,nr,nc;
  MDR *new;

  if (SIZE(m1) < 1)
  {
    new = mdr_Create(0,0);
    return new;
  }
  if (SIZE(m2) < 1)
  {
    new = mdr_Copy(m1);
    return new;
  }

  nr = MAX(MNR (m1), MNR (m2));
  nc = MAX(MNC (m1), MNC (m2));

  new = mdr_Create (nr, nc);
  for (i = 0; i < nr; i++) for (j = 0; j < nc; j++)
      Mdr0 (new,i,j) = (double)
          (mdr0 (m1, MIN(i,MNR(m1)-1), MIN(j,MNC(m1)-1)) >=  mdr0 (m2, MIN(i,MNR(m2)-1), MIN(j,MNC(m2)-1)));
  return (new);
}

MDR *
mdr_And (MDR * m1, MDR * m2)
{
  int i,j,nr,nc;
  MDR *new;

  if (SIZE(m1) < 1)
  {
    new = mdr_Copy(m2);
    mdr_Zero(new);
    return new;
  }
  if (SIZE(m2) < 1)
  {
    new = mdr_Copy(m1);
    mdr_Zero(new);
    return new;
  }

  nr = MAX(MNR (m1), MNR (m2));
  nc = MAX(MNC (m1), MNC (m2));

  if(m1->type == RLAB_TYPE_INT32 && m2->type == RLAB_TYPE_INT32)
  {
    new = mdi_Create (nr, nc);
    for (i = 0; i < nr; i++) for (j = 0; j < nc; j++)
      Mdi0 (new,i,j) = Mdi0 (m1, MIN(i,MNR(m1)-1), MIN(j,MNC(m1)-1)) & Mdi0 (m2, MIN(i,MNR(m2)-1), MIN(j,MNC(m2)-1));
  }
  else
  {
    new = mdr_Create (nr, nc);
    for (i = 0; i < nr; i++) for (j = 0; j < nc; j++)
        Mdr0 (new,i,j) = mdr0 (m1, MIN(i,MNR(m1)-1), MIN(j,MNC(m1)-1)) && mdr0 (m2, MIN(i,MNR(m2)-1), MIN(j,MNC(m2)-1));
  }

  return (new);
}

MDR * mdr_Or (MDR * m1, MDR * m2)
{
  int i,j,nr,nc;
  MDR *new;

  if (SIZE(m1) < 1)
  {
    new = mdr_Copy(m2);
    return new;
  }
  if (SIZE(m2) < 1)
  {
    new = mdr_Copy(m1);
    return new;
  }

  nr = MAX(MNR (m1), MNR (m2));
  nc = MAX(MNC (m1), MNC (m2));

  if(MD_TYPE_INT32(m1) && MD_TYPE_INT32(m2))
  {
    new = mdi_Create (nr, nc);
    for (i = 0; i < nr; i++) for (j = 0; j < nc; j++)
        Mdi0 (new,i,j) = Mdi0 (m1, MIN(i,MNR(m1)-1), MIN(j,MNC(m1)-1)) | Mdi0 (m2, MIN(i,MNR(m2)-1), MIN(j,MNC(m2)-1));
  }
  else
  {
    new = mdr_Create (nr, nc);
    for (i = 0; i < nr; i++) for (j = 0; j < nc; j++)
        Mdr0 (new,i,j) = mdr0 (m1, MIN(i,MNR(m1)-1), MIN(j,MNC(m1)-1)) || mdr0 (m2, MIN(i,MNR(m2)-1), MIN(j,MNC(m2)-1));
  }

  return (new);
}

MDR *
mdr_Not (MDR * m1)
{
  int i, size;
  MDR *new;

  size = SIZE(m1);
  if(MD_TYPE_INT32(m1))
  {
    new = mdi_Create (MNR (m1), MNC (m1));
    for (i = 0; i < size; i++)
    {
      MdiV0 (new, i) = ~MdiV0 (m1, i);
    }
  }
  else
  {
    new = mdr_Create (MNR (m1), MNC (m1));
    for (i= 0; i < size; i++)
    {
      MdrV0 (new, i) = (double) (!MdrV0 (m1, i));
    }
  }

  return (new);
}

/* **************************************************************
 * Transpose a matrix.
 * ************************************************************** */
unsigned char *
data_transpose (unsigned char *data, int nrow, int ncol, int size_a)
{
  unsigned char *new=0;
  int i, j;

  if (nrow*ncol == 0)
    return (new);

  new = (unsigned char *) GC_MALLOC (nrow * ncol * size_a);

  for (i=0; i<nrow; i++)
    for (j=0; j<ncol; j++)
  {
    memcpy (&new[(i*ncol + j)*size_a], &data[(j*nrow +i)*size_a], size_a);
  }

  return (new);
}


MDR *
mdr_Transpose (MDR * m)
{
  MDR *new = (MDR *) GC_MALLOC(sizeof(MDR));

  // first copy the old matrix to new one to preserve type of the original matrix
  memcpy(new, m, sizeof(MDR));

  // create a copy of transposed data
  MDPTR(new) = (void *) data_transpose ((unsigned char *)MDPTR(m),
        m->nrow, m->ncol, rsizeof(m->type));

  // finally switch the rows and columns sizes
  new->nrow = m->ncol;
  new->ncol = m->nrow;

  return (new);
}

/* **************************************************************
 * Reshape a matrix into a colum vector.
 * ************************************************************** */

MDR *
mdr_ReshapeCol (MDR * m)
{
  MDR *new = mdr_Reshape (m, SIZE(m), 1);
  return (new);
}

/*
 * Matrix_is_positive:
 * Return TRUE if ALL elements are >= 0.
 * This function is not meant to be used with
 * complex, or string matrices.
 */

int
mdr_is_positive (MDR * m)
{
  int i, size;
  size = SIZE(m);
  if (size)
  {
    for (i=0; i<size; i++)
    {
      if (mdrV0 (m,i) < 0)
        return (0);
    }
    return (1);
  }
  return (0); // empty matrix is not positive
}

//
// Truncate a matrix in place: meaningful only if (nr,nc) is less
// than the actual dimension of the matrix
//

MDR *
mdr_Truncate (MDR * m, int nr, int nc)
{
  unsigned char *data;
  int size_a;

  // Check first
  if (m->nrow <= nr && m->ncol <= nc)
    return (m);

  size_a = rsizeof(m->type);
  data = (unsigned char *) md_extend ((unsigned char *)MDPTR(m),
          m->nrow, m->ncol, nr, nc, size_a);
  GC_free (MDPTR(m));
  MDPTR(m) = (void *) data;

  m->nrow = nr;
  m->ncol = nc;

  return (m);
}


//
// writing an integer matrix to a serial port
//
extern int RLABPLUS_SERIAL_2BYTETIME_US;

void mdi_intfd_WriteGeneric (char *name, MDR * m)
{
  int i,n;
  unsigned char c;

  int fd = get_int_file_ds(name);
  if (fd < 0)
  {
    printf ("writem: attempt to write to non-existent file/port\n");
    return;
  }

  // flush the port before trying to write to it
  usleep(10 * RLABPLUS_SERIAL_2BYTETIME_US);// sleep a little
  tcflush (fd, TCIOFLUSH);

  n = SIZE(m);

  for (i=0; i<n; i++)
  {
    // coerce to integers
    c = mdiV0(m,i) & 0xff;
    while( write(fd, &c, 1) < 1 )
    {
      usleep(RLABPLUS_SERIAL_2BYTETIME_US); // sleep a little
      tcflush (fd, TCIOFLUSH);              // flush the buffer and try again
    }

    usleep(RLABPLUS_SERIAL_2BYTETIME_US);  // sleep a little
    tcflush (fd, TCIOFLUSH);            // flush the buffer before trying again
  }
}


int
mdr_Sublist (MDR * m, char *name)
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

MDR *
mdr_Inverse (MDR *A)
{
  int n, i;
  MDR *unity, *invA;

  n = MNR(A);
  if (n!=MNC(A)) rerror("A is not square matrix!");

  unity = mdr_Create (n,n);
  mdr_Zero (unity);
  for (i = 1; i <= n; i++)
    Mdr1 (unity, i, i) = 1;

  invA = mdr_Ldivide (A, unity);

  mdr_Destroy(unity);
  return invA;
}

//
// dereferencing MDR - slower than respective macros but
// simplifies handling: cannot be used on the LHS !!!!!
//
double mdr0(MDR *m, int k, int j)
{
  if (k<0 || j<0 || k >= MNR(m) || j >= MNC(m))
    rerror("mdr0: Reference out of bounds! Cannot continue!\n");

  switch(m->type)
  {
    case RLAB_TYPE_DOUBLE:
      if (m->row_dominant)
        return Mdr0_rd(m, k, j);
      else
        return Mdr0(m, k, j);

    case RLAB_TYPE_INT32:
      if (m->row_dominant)
        return (double) Mdi0_rd(m, k, j);
      else
        return (double) Mdi0(m, k, j);

    case RLAB_TYPE_COMPLEX:
      if (m->row_dominant)
        return (double) Mdc0r_rd(m, k, j);
      else
        return (double) Mdc0r(m, k, j);

    default:
      return create_nan();
  }
}

double mdr0_safe(MDR *m, int k, int j)
{
  k = MAX(MIN(MNR(m)-1, k),0);
  j = MAX(MIN(MNC(m)-1, j),0);
  switch(m->type)
  {
    case RLAB_TYPE_DOUBLE:
      if (m->row_dominant)
        return Mdr0_rd(m, k, j);
      else
        return Mdr0(m, k, j);

    case RLAB_TYPE_INT32:
      if (m->row_dominant)
        return (double) Mdi0_rd(m, k, j);
      else
        return (double) Mdi0(m, k, j);

    case RLAB_TYPE_COMPLEX:
      if (m->row_dominant)
        return (double) Mdc0r_rd(m, k, j);
      else
        return (double) Mdc0r(m, k, j);

    default:
      return create_nan();
  }
}

int mdi0(MDR *m, int k, int j)
{
  if (k<0 || j<0 || k >= MNR(m) || j >= MNC(m))
    rerror("mdi0: Reference out of bounds! Cannot continue!\n");

  switch(m->type)
  {
    case RLAB_TYPE_DOUBLE:
      if (m->row_dominant)
        return (int) Mdr0_rd(m, k, j);
      else
        return (int) Mdr0(m, k, j);

    case RLAB_TYPE_INT32:
      if (m->row_dominant)
        return Mdi0_rd(m, k, j);
      else
        return Mdi0(m, k, j);

    case RLAB_TYPE_COMPLEX:
      if (m->row_dominant)
        return (int) Mdc0r_rd(m, k, j);
      else
        return (int) Mdc0r(m, k, j);

    default:
      return (int) create_nan();
  }
}

int mdi0_safe(MDR *m, int k, int j)
{
  k = MAX(MIN(MNR(m)-1, k),0);
  j = MAX(MIN(MNC(m)-1, j),0);
  switch(m->type)
  {
    case RLAB_TYPE_DOUBLE:
      if (m->row_dominant)
        return (int) Mdr0_rd(m, k, j);
      else
        return (int) Mdr0(m, k, j);

    case RLAB_TYPE_INT32:
      if (m->row_dominant)
        return Mdi0_rd(m, k, j);
      else
        return Mdi0(m, k, j);

    case RLAB_TYPE_COMPLEX:
      if (m->row_dominant)
        return (int) Mdc0r_rd(m, k, j);
      else
        return (int) Mdc0r(m, k, j);

    default:
      return (int) create_nan();
  }
}

double mdr1(MDR *m, int k, int j)
{
  if (k<1 || j<1 || k > MNR(m) || j > MNC(m))
    rerror("mdi0: Reference out of bounds! Cannot continue!\n");

  switch(m->type)
  {
    case RLAB_TYPE_DOUBLE:
      return Mdr1(m, k, j);

    case RLAB_TYPE_INT32:
      return (double) Mdi1(m, k, j);

    case RLAB_TYPE_COMPLEX:
      return (double) Mdc1r(m, k, j);

    default:
      return create_nan();
  }
}

double mdr1_safe(MDR *m, int k, int j)
{
  k = MAX(MIN(MNR(m), k),1);
  j = MAX(MIN(MNC(m), j),1);
  switch(m->type)
  {
    case RLAB_TYPE_DOUBLE:
      if (m->row_dominant)
        return Mdr1_rd(m, k, j);
      else
        return Mdr1(m, k, j);

    case RLAB_TYPE_INT32:
      if (m->row_dominant)
        return (double) Mdi1_rd(m, k, j);
      else
        return (double) Mdi1(m, k, j);

    case RLAB_TYPE_COMPLEX:
      if (m->row_dominant)
        return (double) Mdc1r_rd(m, k, j);
      else
        return (double) Mdc1r(m, k, j);

    default:
      return create_nan();
  }
}

int mdi1(MDR *m, int k, int j)
{
  switch(m->type)
  {
    case RLAB_TYPE_DOUBLE:
      if (m->row_dominant)
        return (int) Mdr1_rd(m, k, j);
      else
        return (int) Mdr1(m, k, j);

    case RLAB_TYPE_INT32:
      if (m->row_dominant)
        return Mdi1_rd(m, k, j);
      else
        return Mdi1(m, k, j);

    case RLAB_TYPE_COMPLEX:
      if (m->row_dominant)
        return (int) Mdc1r_rd(m, k, j);
      else
        return (int) Mdc1r(m, k, j);

    default:
      return (int) create_nan();
  }
}

int mdi1_safe(MDR *m, int k, int j)
{
  k = MAX(MIN(MNR(m), k),1);
  j = MAX(MIN(MNC(m), j),1);
  switch(m->type)
  {
    case RLAB_TYPE_DOUBLE:
      if (m->row_dominant)
        return (int) Mdr1_rd(m, k, j);
      else
        return (int) Mdr1(m, k, j);

    case RLAB_TYPE_INT32:
      if (m->row_dominant)
        return Mdi1_rd(m, k, j);
      else
        return Mdi1(m, k, j);

    case RLAB_TYPE_COMPLEX:
      if (m->row_dominant)
        return (int) Mdc1r_rd(m, k, j);
      else
        return (int) Mdc1r(m, k, j);

    default:
      return (int) create_nan();
  }
}

double mdrV0(MDR *m, int k)
{
  if (k<0 || k >= SIZE(m))
  {
    printf("mdrV0: index %i, size %i\n", k, SIZE(m));
    rerror("mdrV0: Reference out of bounds! Cannot continue!\n");
  }

  switch(m->type)
  {
    case RLAB_TYPE_DOUBLE:
      return MdrV0(m, k);

    case RLAB_TYPE_INT32:
      return (double) MdiV0(m, k);

    case RLAB_TYPE_COMPLEX:
      return MdcV0r(m, k);

    default:
      return create_nan();
  }
}

double mdrV0_safe(MDR *m, int k)
{
  int i = MAX(MIN(SIZE(m)-1, k),0);

  switch(m->type)
  {
    case RLAB_TYPE_DOUBLE:
      return MdrV0(m, i);

    case RLAB_TYPE_INT32:
      return (double) MdiV0(m, i);

    case RLAB_TYPE_COMPLEX:
      return MdcV0r(m, i);

    default:
      return create_nan();
  }
}

int mdiV0(MDR *m, int k)
{
  if (k<0 || k >= SIZE(m))
    rerror("mdiV0: Reference out of bounds! Cannot continue!\n");

  switch(m->type)
  {
    case RLAB_TYPE_DOUBLE:
      return (int) MdrV0(m, k);

    case RLAB_TYPE_INT32:
      return MdiV0(m, k);

    case RLAB_TYPE_COMPLEX:
      return (int) MdcV0r(m, k);

    default:
      return 0;
  }
}

int mdiV0_safe(MDR *m, int k)
{
  int i = MAX(MIN(SIZE(m)-1, k),0);

  switch(m->type)
  {
    case RLAB_TYPE_DOUBLE:
      return (int) MdrV0(m, i);

    case RLAB_TYPE_INT32:
      return MdiV0(m, i);

    case RLAB_TYPE_COMPLEX:
      return (int) MdcV0r(m, i);

    default:
      return create_nan();
  }
}


double mdrV1(MDR *m, int k)
{
  if (k<1 || k > SIZE(m))
    rerror("mdrV1: Reference out of bounds! Cannot continue!\n");

  switch(m->type)
  {
    case RLAB_TYPE_DOUBLE:
      return MdrV1(m, k);

    case RLAB_TYPE_INT32:
      return (double) MdiV1(m, k);

    case RLAB_TYPE_COMPLEX:
      return MdcV1r(m, k);

    default:
      return create_nan();
  }
}

double mdrV1_safe(MDR *m, int k)
{
  int i = MAX(MIN(SIZE(m), k),1);

  switch(m->type)
  {
    case RLAB_TYPE_DOUBLE:
      return MdrV1(m, i);

    case RLAB_TYPE_INT32:
      return (double) MdiV1(m, i);

    case RLAB_TYPE_COMPLEX:
      return MdcV1r(m, i);

    default:
      return create_nan();
  }
}



int mdiV1(MDR *m, int k)
{
  if (k<1 || k > SIZE(m))
    rerror("mdrV1: Reference out of bounds! Cannot continue!\n");

  switch(m->type)
  {
    case RLAB_TYPE_DOUBLE:
      return (int) MdrV1(m, k);

    case RLAB_TYPE_INT32:
      return MdiV1(m, k);

    case RLAB_TYPE_COMPLEX:
      return (int) MdcV1r(m, k);

    default:
      return 0;
  }
}

int mdiV1_safe(MDR *m, int k)
{
  int i = MAX(MIN(SIZE(m), k),1);

  switch(m->type)
  {
    case RLAB_TYPE_DOUBLE:
      return (int) MdrV1(m, i);

    case RLAB_TYPE_INT32:
      return MdiV1(m, i);

    case RLAB_TYPE_COMPLEX:
      return (int) MdcV1r(m, i);

    default:
      return create_nan();
  }
}




int md_transpose_insitu(unsigned char *a, int m, int n, int size_a)
{
//    ACM ALGORITHM 380 - REVISED
//    by Esko G. Cate and David W. Twigg of
//    Boeing Computer Services, Inc., P.O. Box 24346, Seattle, WA
//    Publisher ACM  New York, NY, USA
// A IS A ONE-DIMENSIONAL ARRAY OF LENGTH MN=M*N, WHICH
// CONTAINS THE MXN MATRIX TO BE TRANSPOSED (STORED
// COLUMWISE). MOVE IS A ONE-DIMENSIONAL ARRAY OF LENGTH IWRK
// USED TO STORE INFORMATION TO SPEED UP THE PROCESS.  THE
// VALUE IWRK=(M+N)/2 IS RECOMMENDED. IOK INDICATES THE
// SUCCESS OR FAILURE OF THE ROUTINE.
// NORMAL RETURN  IOK=0
// ERRORS         IOK=-1 ,MN NOT EQUAL TO M*N
//                IOK=-2 ,IWRK NEGATIVE OR ZERO
//                IOK.GT.0, (SHOULD NEVER OCCUR),IN THIS CASE
// WE SET IOK EQUAL TO THE FINAL VALUE OF I WHEN THE SEARCH
// IS COMPLETED BUT SOME LOOPS HAVE NOT BEEN MOVED
// NOTE * MOVE(I) WILL STAY ZERO FOR FIXED POINTS
  int i, i1, i2, j, ir0, ir1=0, ir2, im, mymax, kmi, i1c, i2c;

  //  NOTHING TO DO
  if (m<2 || n<2)
    return 0;

  unsigned char *bcd = GC_malloc(3 * size_a);
  unsigned char *b, *c, *d;
  b = &bcd[0];
  c = &bcd[size_a];
  d = &bcd[2*size_a];

  if (m == n)
  {
    // IF MATRIX IS SQUARE,EXCHANGE ELEMENTS A(I,J) AND A(J,I).
    for (i=0; i<(n-1); i++) // DO 150 I=1,N1
    {
      // DO 140 J=J1,N
      for (j=i+1; j<n; j++)
      {
        i1 = i + j*n;
        i2 = j + i*m;
        memcpy(b, &a[i1*size_a], size_a);                     // B = A(I1)
        memcpy(&a[i1*size_a], &a[i2*size_a], size_a);     // A(I1) = A(I2)
        memcpy(&a[i2*size_a], b, size_a);                     // A(I2) = B
      }
    }

    // release storage
    GC_FREE (bcd);
    return 0;
  }

  int iwrk = (int)(((double)m + (double)n) * 0.5 + 1.0);
  int *move = GC_malloc(iwrk * sizeof(int));

  int ncount = 2;
  int k = m*n - 1;
  for (i=0; i<iwrk; i++)
        move[i] = 0;

  if (m>2 && n>2)
  {
    // CALCULATE THE NUMBER OF FIXED POINTS, EUCLIDS ALGORITHM
    // FOR GCD(M-1,N-1).
    ir2 = m - 1;
    ir1 = n - 1;
    ir0 = 1;
    while (ir0 != 0)
    {
      ir0 = ir2 % ir1; //ir0 = mod(ir2,ir1);
      ir2 = ir1;
      ir1 = ir0;
    }
    ncount = ncount + ir2 - 1;
  }

  // SET INITIAL VALUES FOR SEARCH
  i = 1;
  im = m;

  // AT LEAST ONE LOOP MUST BE RE-ARRANGED
  goto l80;

  // SEARCH FOR LOOPS TO REARRANGE
l40:
  mymax = k - i;
  i = i + 1;
  if (i > mymax)
  {
    GC_FREE (bcd);
    GC_FREE (move);
    return i;
  }
  im = im + m;
  if (im > k)
    im = im - k;
  i2 = im;
  if (i == i2)
    goto l40;
  if (i > iwrk)
    goto l60;
  if (move[i-1] == 0)
    goto l80;

  goto l40;

l50:
  i2 = m*i1 - k*((double)i1/(double)n);

l60:
  if (i2 <= i || i2>=mymax)
    goto l70;
  i1 = i2;
  goto l50;

l70:
  if (i2 != i)
    goto l40;

l80:
  // REARRANGE THE ELEMENTS OF A LOOP AND ITS COMPANION LOOP
  i1 = i;
  kmi = k - i;
  memcpy(b, &a[(i1)*size_a], size_a); // B = A(I1+1)
  i1c = kmi;
  memcpy(c, &a[(i1c)*size_a], size_a); // C = A(I1C+1)

l90:
  i2 = m*i1 - k*((int)((double)i1/(double)n));
  i2c = k - i2;
  if (i1 <= iwrk)
    move[i1-1] = 2;
  if (i1c <= iwrk)
    move[i1c-1] = 2;
  ncount = ncount + 2;
  if (i2 == i)
    goto l110;
  if (i2 == kmi)
    goto l100;
  memcpy(&a[(i1)*size_a], &a[(i2)*size_a], size_a);  // A(I1+1) = A(I2+1)
  memcpy(&a[(i1c)*size_a], &a[(i2c)*size_a], size_a);  // A(I1C+1) = A(I2C+1)
  i1 = i2;
  i1c = i2c;
  goto l90;

l100:
  // FINAL STORE AND TEST FOR FINISHED
  memcpy(d, b, size_a);  // D = B
  memcpy(b, c, size_a);  // B = C
  memcpy(c, d, size_a);  // C = D

l110:
  memcpy(&a[(i1)*size_a], b, size_a);    // A(I1+1) = B
  memcpy(&a[(i1c)*size_a], c, size_a);   // A(I1C+1) = C
  if (ncount < m*n)
    goto l40;

  // normal return
  GC_FREE (bcd);
  GC_FREE (move);
  return i;
}

void mdr_Transpose_inplace (MDR * m)
{
  if (SIZE(m)<1)
    return;

  if (EQVECT(m))
    goto exit;

  md_transpose_insitu((unsigned char *) MDPTR(m), MNR(m), MNC(m), rsizeof(m->type));

exit:
  SWAP(MNR(m),MNC(m),int);
}

int mdr_vector_issorted(MDR *x)
{
  int i;
  for(i=1; i<SIZE(x); i++)
  {
    if (mdrV0(x,i-1) >= mdrV0(x,i))
      return 0;
  }
  return 1;
}

MDR * mdr_FindRootsFromInterpolationTable(MDR * x1, double yoff)
{
  MDR *w=NULL;

  if (SIZE(x1)<1)
    goto exit;

  int nr=MNR(x1), nnz=0, i, j;
  MDR *w1 = mdr_Create(1, nr);

  // check that the first coordinate is sorted in ascending order
  for (i=0; i<(nr-1); i++)
  {
    if (mdr0(x1,i,0) >= mdr0(x1,i+1,0))
      goto exit;
  }

  i=0;
  nnz = 0;
  while (i < nr)
  {
    // special case: long sequence of [ ... x1, yoff; x2, yoff; ... ]
    if (mdr0(x1,i,1)==yoff)
    {
      j = i;
      while (mdr0(x1,j,1)==yoff)
      {
        j++;
        if (j > nr - 1)
          break;
      }
      MdrV0(w1, nnz) = 0.5 * (mdr0(x1,i,0) + mdr0(x1,j-1,0));
      nnz++;

      i = j;

      if (i < nr)
        continue;

      break;
    }

    if (i+1 < nr)
    {
      if ((mdr0(x1,i,1)-yoff)*(mdr0(x1,i+1,1)-yoff) < 0)
      {
        // linear interpolation
        MdrV0(w1, nnz) =
            (mdr0(x1,i+1,0)*(mdr0(x1,i,1)-yoff)-mdr0(x1,i,0)*(mdr0(x1,i+1,1)-yoff))
            / (mdr0(x1,i,1)-mdr0(x1,i+1,1));
        nnz++;
      }
    }

    i++;
  }

  if (nnz)
  {
    w = mdr_Create(1, nnz);
    for (i=0; i<nnz; i++)
      MdrV0(w, i) = MdrV0(w1, i);
  }

  mdr_Destroy(w1);

exit:

  return w;
}


#include "sort.h"
MDR * mdr_Copy2double_and_Sort(MDR * x, int row_col_flat)
{
  if (SIZE(x)<1)
    return NULL;

  int n, i, j, nr, nc;

  MDR *rval=0, *x_idx=0;

  switch (row_col_flat)
  {
    case 1:
    case 0:
      rval = mdr_Float_BF(x);
      nr = MNR(x);
      nc = MNC(x);
      if (row_col_flat == 1)
      {
        // 1 - row dominant (C)
        mdr_Transpose_inplace(rval);
      }
      // 0 - col dominant (Fortran)
      nr = MNR(rval);
      nc = MNC(rval);
      x_idx = mdr_Create(1, nr);
      for (j=0; j<nc; j++)
      {
        for (i=0; i<nr; i++)
        {
          MdrV0(x_idx,i) = i;
        }
        // sort
        r_sort ((double *) &Mdr0(rval,0,j), 0, nr-1, (double *)MDRPTR(x_idx));
      }
      if (row_col_flat == 1)
      {
        // 1 - row dominant (C)
        mdr_Transpose_inplace(rval);
      }
      break;

    case -1:
      // flat
      n = SIZE(rval);
      x_idx = mdr_Create(1, n);
      rval = mdr_Float_BF(x);
      for (i=0; i<n; i++)
      {
        MdrV0(x_idx,i) = i;
      }
      // sort
      r_sort ((double *) MDRPTR(rval), 0, n-1, (double *)MDRPTR(x_idx));
  }

  if (x_idx)
    mdr_Destroy (x_idx);

  return rval;
}


int mdr_vector_isbounded_from_below(MDR *x, double lb)
{
  int i;
  for(i=0; i<SIZE(x); i++)
  {
    if (mdrV0(x,i) < lb)
      return 0;
  }
  return 1;
}

int mdr_vector_isbounded_from_above(MDR *x, double ub)
{
  int i;
  for(i=0; i<SIZE(x); i++)
  {
    if (mdrV0(x,i) > ub)
      return 0;
  }
  return 1;
}

int mdr_vector_isnan(MDR *x)
{
  if (!MD_TYPE_DOUBLE(x))
    return 0;

  int i,rval=0;
  for(i=0; i<SIZE(x); i++)
  {
    if (isnand(mdrV0(x,i)))
      rval++;
  }
  return rval;
}

int mdr_rgemtv(MDR *A, double alfa, MDR *x, double beta, MDR **y)
{
  //
  // y -> alfa * A^T * x + beta * y
  //
  int n = MNR(A);
  int m = MNC(A);
  int t = 'T';
  int incxy=1;

  if (SIZE(x) != n)
    return 3;

  if (!y)
  {
    *y = mdr_Create(n,1);
    mdr_Zero(*y);
  }

  if (SIZE(*y) != m)
    return 5;

  return RGEMV(&t, &n, &m, &alfa, MDRPTR(A), &n, MDRPTR(x), &incxy, &beta, MDRPTR(*y), &incxy);
}

int mdr_rgemv(MDR *A, double alfa, MDR *x, double beta, MDR **y)
{
  //
  // y -> alfa * A * x + beta * y
  //
  int n = MNR(A);
  int m = MNC(A);
  int t = 'N';
  int incxy=1;

  if (SIZE(x) != m)
    return 3;

  if (!y)
  {
    *y = mdr_Create(n,1);
    mdr_Zero(*y);
  }

  if (SIZE(*y) != n)
    return 5;

  return RGEMV(&t, &n, &m, &alfa, MDRPTR(A), &n, MDRPTR(x), &incxy, &beta, MDRPTR(*y), &incxy);
}

int mdr_matinvd(MDR *A, double *det)
{
  // find inverse of A and store it in A
  int i, n=MNR(A);
  int info;
  int *ipiv = (int    *) GC_malloc(n  * sizeof(int));

  RGETRF(&n, &n, MDRPTR(A), &n, ipiv, &info);         // do LU factorization first

  if (!info)
  {
    //
    // find the determinant of the inverse
    //
    if (det)
    {
      *det=Mdr0(A,0,0);
      for (i=1; i<n; i++)
        (*det) = (*det) * Mdr0(A,i,i);
      (*det) = 1.0/(*det);
    }

    // query optimal size of work array
    int lwork=-1;
    double dwork;
    RGETRI(&n, MDRPTR(A), &n, ipiv, &dwork, &lwork, &info);

    // calculate inverse from LU factors
    lwork = dwork;
    double *work = (double *) GC_malloc(lwork * sizeof(double));
    RGETRI(&n, MDRPTR(A), &n, ipiv, work, &lwork, &info);

    GC_free (work);
  }


  GC_free (ipiv);

  return info;
}









