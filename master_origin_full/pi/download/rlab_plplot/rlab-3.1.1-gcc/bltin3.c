/* bltin3.c */

/*  This file is a part of RLaB ("Our"-LaB)
   Copyright (C) 1996  Ian R. Searle

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
#include "bltin.h"
#include "util.h"
#include "ent.h"
#include "class.h"
#include "rfileio.h"
#include "symbol.h"
#include "btree.h"
#include "mathl.h"

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>

/*
 * Header files from the available classes...
 */

#include "ent.h"
#include "btree.h"
#include "bltin.h"
#include "function.h"
#include "mdr.h"
#include "mdrf1.h"
#include "mdrf2.h"
#include "mdrf3.h"
#include "mdc.h"
#include "mdcf1.h"
#include "mdcf2.h"
#include "mdr_mdc.h"
#include "mds.h"
#include "mdsf1.h"
#include "mdr_mds.h"
#include "msr.h"
#include "msrf1.h"
#include "msrf2.h"
#include "msrf3.h"
#include "msc.h"
#include "mscf1.h"
#include "mscf2.h"

#include "rlab_solver_parameters_names.h"
#include "rlab_macros.h"
#include "rlab_macros_code.h"

static OpDef sparse_method[NCL];
static OpDef dense_method[NCL];
static OpDef spconvert_method[NCL];
static OpDef sporder_method[NCL];
static OpDef spwrite_sparse_method[NCL];
static OpDef spwrite_graph_method[NCL];

/* **************************************************************
 * Initialize the built-ins...
 * ************************************************************** */

void
class_bltin3_init (void)
{
  /*
   * sparse ()
   */

  sparse_method[MATRIX_DENSE_REAL].type = MATRIX_SPARSE_REAL;
  sparse_method[MATRIX_DENSE_REAL].op = (void *) msr_Sparse;

  sparse_method[MATRIX_DENSE_COMPLEX].type = MATRIX_SPARSE_COMPLEX;
  sparse_method[MATRIX_DENSE_COMPLEX].op = (void *) msc_Sparse;

  sparse_method[MATRIX_SPARSE_REAL].type = MATRIX_SPARSE_REAL;
  sparse_method[MATRIX_SPARSE_REAL].op = (void *) msr_ReSparse;

  sparse_method[MATRIX_SPARSE_COMPLEX].type = MATRIX_SPARSE_COMPLEX;
  sparse_method[MATRIX_SPARSE_COMPLEX].op = (void *) msc_ReSparse;

  /*
   * dense ()
   */

  dense_method[MATRIX_DENSE_REAL].type = MATRIX_DENSE_REAL;
  dense_method[MATRIX_DENSE_REAL].op = (void *) mdr_Dense;

  dense_method[MATRIX_DENSE_COMPLEX].type = MATRIX_DENSE_COMPLEX;
  dense_method[MATRIX_DENSE_COMPLEX].op = (void *) mdc_Dense;

  dense_method[MATRIX_SPARSE_REAL].type = MATRIX_DENSE_REAL;
  dense_method[MATRIX_SPARSE_REAL].op = (void *) msr_Dense;

  dense_method[MATRIX_SPARSE_COMPLEX].type = MATRIX_DENSE_COMPLEX;
  dense_method[MATRIX_SPARSE_COMPLEX].op = (void *) msc_Dense;

  /*
   * spconvert ()
   */

  spconvert_method[MATRIX_DENSE_REAL].type = 0;
  spconvert_method[MATRIX_DENSE_REAL].op = (void *) mdr_Spconvert;

  spconvert_method[MATRIX_SPARSE_REAL].type = MATRIX_DENSE_REAL;
  spconvert_method[MATRIX_SPARSE_REAL].op = (void *) msr_Spconvert;

  /*
   * sporder ()
   */

#ifdef HAVE_METIS3
  sporder_method[MATRIX_SPARSE_REAL].type = MATRIX_DENSE_REAL;
  sporder_method[MATRIX_SPARSE_REAL].op = (void *) msr_SpOrder;

  sporder_method[MATRIX_SPARSE_COMPLEX].type = MATRIX_DENSE_REAL;
  sporder_method[MATRIX_SPARSE_COMPLEX].op = (void *) msc_SpOrder;
#endif /* HAVE_METIS3 */

  /*
   * spwrite ()
   */

  spwrite_sparse_method[MATRIX_SPARSE_REAL].type = MATRIX_DENSE_REAL;
  spwrite_sparse_method[MATRIX_SPARSE_REAL].op = (void *) msr_WriteSparse;

  spwrite_sparse_method[MATRIX_SPARSE_COMPLEX].type = MATRIX_DENSE_REAL;
  spwrite_sparse_method[MATRIX_SPARSE_COMPLEX].op = (void *) msc_WriteSparse;

  spwrite_graph_method[MATRIX_SPARSE_REAL].type = MATRIX_DENSE_REAL;
  spwrite_graph_method[MATRIX_SPARSE_REAL].op = (void *) msr_WriteGraph;

  spwrite_graph_method[MATRIX_SPARSE_COMPLEX].type = MATRIX_DENSE_REAL;
  spwrite_graph_method[MATRIX_SPARSE_COMPLEX].op = (void *) msc_WriteGraph;

}

/* **************************************************************
 * sparse(): Convert full (dense) storage -> sparse storage
 * ************************************************************** */

Ent *
Sparse (int nargs, Datum args[])
{
  int type, rtype;
  void *rval;
  void *(*vfptr) ();
  Ent *rent, *e;

  if (nargs != 1)
  {
    rerror ("sparse: 1 argument allowed");
  }

  e = bltin_get_ent (args[0]);
  type = ent_type (e);

  rtype = sparse_method[type].type;
  vfptr = (VFPTR) sparse_method[type].op;

  if (vfptr == 0)
  {
    fprintf (stderr, "Entity type: %s\n", etd (e));
    rerror ("sparse() operation not supported");
  }

  rval = (*vfptr) (ent_data (e));

  rent = ent_Create ();
  ent_data (rent) = rval;
  ent_type (rent) = rtype;

  ent_Clean (e);
  return (rent);
}

/* **************************************************************
 * full(): Convert sparse storage -> full (dense) storage.
 * ************************************************************** */

Ent *
Dense (int nargs, Datum args[])
{
  int type, rtype;
  void *rval;
  void *(*vfptr) ();
  Ent *rent, *e;

  if (nargs != 1)
  {
    rerror ("full: 1 argument allowed");
  }

  e = bltin_get_ent (args[0]);
  type = ent_type (e);

  rtype = dense_method[type].type;
  vfptr = (VFPTR) dense_method[type].op;

  if (vfptr == 0)
  {
    fprintf (stderr, "Entity type: %s\n", etd (e));
    rerror ("full() operation not supported");
  }

  rval = (*vfptr) (ent_data (e));

  /*
   * Special thing here. Check that the return data is
   * not the same as the entities data. If it is, then
   * return the original entitiy.
   */

  if (rval == ent_data (e))
  {
    return (e);
  }

  rent = ent_Create ();
  ent_data (rent) = rval;
  ent_type (rent) = rtype;

  ent_Clean (e);
  return (rent);
}

/* **************************************************************
 * Spconvert: Convert a 3 column representation to sparse, or
 *            convert a sparse representation to 3 column.
 * Note: 4 column will come in the future for complex sparse and
 *       dense matrices.
 * ************************************************************** */

Ent *
Spconvert (int nargs, Datum args[])
{
  int type, rtype;
  void *rval;
  void *(*vfptr) ();
  Ent *rent, *e;

  if (nargs != 1)
  {
    rerror ("spconvert: 1 argument allowed");
  }

  e = bltin_get_ent (args[0]);
  type = ent_type (e);

  vfptr = (VFPTR) spconvert_method[type].op;

  if (vfptr == 0)
  {
    fprintf (stderr, "Entity type: %s\n", etd (e));
    rerror ("spconvert() operation not supported");
  }

  rval = (*vfptr) (ent_data (e), &rtype);

  rent = ent_Create ();
  ent_data (rent) = rval;
  ent_type (rent) = rtype;

  ent_Clean (e);
  return (rent);
}

/* **************************************************************
 * SpOrder: Compute the fill-reducing permutation matrix.
 * ************************************************************** */

Ent *
SpOrder (int nargs, Datum args[])
{
  int type, rtype;
  void *rval;
  void *(*vfptr) ();
  Ent *rent, *e;

  if (nargs != 1)
  {
    rerror ("sporder: 1 argument allowed");
  }

  e = bltin_get_ent (args[0]);
  type = ent_type (e);

  rtype = sporder_method[type].type;
  vfptr = (VFPTR) sporder_method[type].op;

  if (vfptr == 0)
  {
    fprintf (stderr, "Entity type: %s\n", etd (e));
    rerror ("sporder() operation not supported");
  }

  rval = (*vfptr) (ent_data (e));

  rent = ent_Create ();
  ent_data (rent) = rval;
  ent_type (rent) = rtype;

  ent_Clean (e);
  return (rent);
}

/* **************************************************************
 * SpWrite: Write out a sparse matrix in various formats.
 * write_type: 0 (default) write out the sparse-internal rep.
 *             1 write out the graph representation (structure only)
 * ************************************************************** */

Ent *
SpWrite (int nargs, Datum args[])
{
  int type, write_type = 0;
  char *string, *stype;
  void *(*vfptr) () = 0;
  Ent *rent, *e1, *e2, *e3 = 0;
  FILE *fn;

  if (nargs > 3 || nargs < 2)
  {
    fprintf (stdout, "spwrite: Export a sparse matrix from RLaB to a file.\n");
    fprintf (stdout, "spwrite: Format:\n");
    fprintf (stdout, "spwrite:   spwrite(\"filename\", x, \"fmt\")\n");
    fprintf (stdout, "spwrite: where \"filename\" is a string, x is a\n");
    fprintf (stdout, "spwrite: sparse matrix and \"fmt\" is a format in\n");
    fprintf (stdout,
	     "spwrite: which the matrix is going to be written in the\n");
    fprintf (stdout, "spwrite: file. Acceptable values are \"hb\", for \n");
    fprintf (stdout,
	     "spwrite: Harwell-Boing, \"smms\" for the SMMS-compatible,\n");
    fprintf (stdout,
	     "spwrite: \"csr\" for compressed sparse row, \"ps\" for\n");
    fprintf (stdout, "spwrite: post-script (eps) compatible plot of general\n");
    fprintf (stdout,
	     "spwrite: matrix structure, \"graph\" for metis-compatible,\n");
    fprintf (stdout, "spwrite: and \"sparse\" for default sparse format.\n");
    rerror ("3 arguments allowed");
  }
  // The output filename
  e1 = bltin_get_ent (args[0]);
  string = class_char_pointer (e1);

  // The sparse matrix.
  e2 = bltin_get_ent (args[1]);
  type = ent_type (e2);

  // The sparse matrix output format.
  if (nargs == 3)
  {
    e3 = bltin_get_ent (args[2]);
    stype = class_char_pointer (e3);
    if (!strcmp ("graph", stype))
      write_type = 1;
    if (!strcmp (stype, "hb"))
      write_type = 2;
    if (!strcmp (stype, "smms"))
      write_type = 3;
    if (!strcmp (stype, "csr"))
      write_type = 4;
    if (!strcmp (stype, "ps"))
      write_type = 5;
    else if (!strcmp ("sparse", stype))
      write_type = 0;
  }
  else
    write_type = 0;

  /* Try and open() the file */
  if ((fn = get_file_ds (string, "w", 0)) == 0)
  {
    fprintf (stderr, "%s, cannot open for read\n", string);
    rent = ent_Create ();

    ent_data (rent) = mdr_CreateScalar (0.0);
    ent_SetType (rent, MATRIX_DENSE_REAL);

    ent_Clean (e1);
    ent_Clean (e2);
    return (rent);
  }

  if (write_type == 0)
    vfptr = (VFPTR) spwrite_sparse_method[type].op;
  else if (write_type == 2)
    vfptr = (VFPTR) msr_sparskit_prtmt;
  else if (write_type == 3)
    vfptr = (VFPTR) msr_sparskit_smms;
  else if (write_type == 4)
    vfptr = (VFPTR) msr_sparskit_dump;
  else if (write_type == 5)
    vfptr = (VFPTR) msr_sparskit_pspltm;
  else if (write_type == 1)
    vfptr = (VFPTR) spwrite_graph_method[type].op;

  if (vfptr == 0)
  {
    fprintf (stderr, "Entity type: %s\n", etd (e2));
    rerror ("spwrite() operation not supported");
  }

  if (write_type == 2)
    (*vfptr) (ent_data (e2), string);
  else if (write_type == 3)
    (*vfptr) (ent_data (e2), string);
  else if (write_type == 4)
    (*vfptr) (ent_data (e2), string);
  else if (write_type == 5)
    (*vfptr) (ent_data (e2), string);
  else
    (*vfptr) (ent_data (e2), fn);
  rent = ent_Create ();

  ent_data (rent) = mdr_CreateScalar (1.0);
  ent_SetType (rent, MATRIX_DENSE_REAL);

  ent_Clean (e1);
  ent_Clean (e2);
  if (nargs == 3)
    ent_Clean (e3);

  // Close the file
  close_file_ds (string);

  return (rent);
}

/* **************************************************************
 * readgraph.c: Read a file containing a graph description ala
 *              Metis and Chaco.
 * ************************************************************** */

#include "rlab.h"
#include "mem.h"
#include "bltin.h"
#include "util.h"
#include "ent.h"
#include "class.h"
#include "rfileio.h"
#include "msr.h"

#include <stdio.h>

MSR *readgraph (FILE * fn);

Ent *
ReadGraph (int nargs, Datum args[])
{
  char *fnstring;
  Ent *rent, *e;
  FILE *fn;
  MSR *m;

  if (nargs != 1)
  {
    rerror ("readgraph: 1 argument allowed");
  }

  e = bltin_get_ent (args[0]);
  fnstring = class_char_pointer (e);

  if ((fn = get_file_ds (fnstring, "r", 0)) == 0)
  {
    fprintf (stderr, "%s, cannot open for read\n", fnstring);
    rerror ("");
  }

  m = readgraph (fn);

  /* Return the new sparse matrix. */
  rent = ent_Create ();
  ent_data (rent) = m;
  ent_type (rent) = MATRIX_SPARSE_REAL;

  ent_Clean (e);
  return (rent);
}

#define MAXLINE  128*1024

MSR *
readgraph (FILE * fn)
{
  char line[MAXLINE + 1], *oldstr, *newstr;
  int fmt, nvtxs;
  int i, itmp, k, nedges, nnz;
  int *ia, *ja;
  double *d;
  MSR *m;

  /* Read the number of vertices and edges. */
  do
  {
    fgets (line, MAXLINE, fn);
  }
  while (line[0] == '%');

  if (sscanf (line, "%d %d %d", &nvtxs, &nedges, &fmt) == 2)
    fmt = 0;

  if (fmt != 0)
    rerror ("readgraph: cannot read weighted graphs\n");

  nnz = 2 * nedges;

  /*
   * Allocate space for the graph.
   */

  ia = (int *) GC_MAIOP ((nvtxs + 1) * sizeof (int));
  ja = (int *) GC_MAIOP (nnz * sizeof (int));
  d = (double *) GC_MAIOP (nnz * sizeof (double));

  ia[0] = 1;

  for (i = 0; i < nvtxs; i++)
  {
    /* Read in the next adjancey list. */
    do
    {
      fgets (line, MAXLINE, fn);
    }
    while (line[0] == '%');

    oldstr = line;
    newstr = NULL;

    if (strlen (line) == MAXLINE)
      rerror ("readgraph: Buffer for fgets not big enough!");

    k = 0;
    for (;;)
    {
      itmp = (int) strtol (oldstr, &newstr, 10);
      if (oldstr == newstr)
	break;

      ja[(ia[i] - 1) + k] = itmp;
      d[(ia[i] - 1) + k++] = 1.0;
      oldstr = newstr;
    }
    /* Update the row pointer. */
    ia[i + 1] = ia[i] + k;
  }

  m = msr_Create (nvtxs, nvtxs);
  m->nnz = nnz;
  m->ia = ia;
  m->ja = ja;
  m->d = d;

  return m;
}


#undef THIS_SOLVER
#define THIS_SOLVER "rot3d"
Ent * ent_3d_rotate (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *e3=0;
  MDR *x=0, *k=0;
  double theta=0.0;
  int i;

  if (nargs != 3)
    goto _exit;

  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) != MATRIX_DENSE_REAL)
    goto _exit;
  x = ent_data(e1);
  if(MNC(x)!=3)
  {
    x = 0;
    goto _exit;
  }
  x = mdr_Float_BF(x); // copy and coerce 'x' to double

  e2 = bltin_get_ent (args[1]);
  if (ent_type (e2) != MATRIX_DENSE_REAL)
  {
    mdr_Destroy(x);
    x = 0;
    goto _exit;
  }
  theta = class_double(e2);

  e3 = bltin_get_ent (args[2]);
  if (ent_type (e3) != MATRIX_DENSE_REAL)
  {
    mdr_Destroy(x);
    x = 0;
    goto _exit;
  }
  k = ent_data(e3);
  if(SIZE(k)!=3)
  {
    mdr_Destroy(x);
    x = 0;
    goto _exit;
  }

  for (i=0; i<MNR(x); i++)
  {
    if (rotate3d(&Mdr0(x,i,0), &Mdr0(x,i,1), &Mdr0(x,i,2), mdrV0(k,0), mdrV0(k,1), mdrV0(k,2), theta))
    {
      mdr_Destroy(x);
      x=0;
      goto _exit;
    }
  }


_exit:

  ent_Clean (e1);
  ent_Clean (e2);
  ent_Clean (e3);
  return ent_Assign_Rlab_MDR(x);
}

#undef THIS_SOLVER
#define THIS_SOLVER "rotx"
Ent * ent_rotate_x (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0;
  MDR *x=0;
  double theta=0.0;
  int i;

  if (nargs != 2)
    goto _exit;

  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) != MATRIX_DENSE_REAL)
    goto _exit;
  x = ent_data(e1);
  if(MNC(x)!=3)
  {
    x = 0;
    goto _exit;
  }
  x = mdr_Float_BF(x); // copy and coerce 'x' to double

  e2 = bltin_get_ent (args[1]);
  if (ent_type (e2) != MATRIX_DENSE_REAL)
  {
    mdr_Destroy(x);
    x = 0;
    goto _exit;
  }
  theta = class_double(e2);

  for (i=0; i<MNR(x); i++)
  {
    if (rotate3d(&Mdr0(x,i,0), &Mdr0(x,i,1), &Mdr0(x,i,2), 1.0, 0.0, 0.0, theta))
    {
      mdr_Destroy(x);
      x=0;
      goto _exit;
    }
  }


_exit:

  ent_Clean (e1);
  ent_Clean (e2);
  return ent_Assign_Rlab_MDR(x);
}

#undef THIS_SOLVER
#define THIS_SOLVER "roty"
Ent * ent_rotate_y (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0;
  MDR *x=0;
  double theta=0.0;
  int i;

  if (nargs != 2)
    goto _exit;

  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) != MATRIX_DENSE_REAL)
    goto _exit;
  x = ent_data(e1);
  if(MNC(x)!=3)
  {
    x = 0;
    goto _exit;
  }
  x = mdr_Float_BF(x); // copy and coerce 'x' to double

  e2 = bltin_get_ent (args[1]);
  if (ent_type (e2) != MATRIX_DENSE_REAL)
  {
    mdr_Destroy(x);
    x = 0;
    goto _exit;
  }
  theta = class_double(e2);

  for (i=0; i<MNR(x); i++)
  {
    if (rotate3d(&Mdr0(x,i,0), &Mdr0(x,i,1), &Mdr0(x,i,2), 0.0, 1.0, 0.0, theta))
    {
      mdr_Destroy(x);
      x=0;
      goto _exit;
    }
  }


_exit:

  ent_Clean (e1);
  ent_Clean (e2);
  return ent_Assign_Rlab_MDR(x);
}

#undef THIS_SOLVER
#define THIS_SOLVER "rotz"
Ent * ent_rotate_z (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0;
  MDR *x=0;
  double theta=0.0;
  int i;

  if (nargs != 2)
    goto _exit;

  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) != MATRIX_DENSE_REAL)
    goto _exit;
  x = ent_data(e1);
  if(MNC(x)!=3)
  {
    x = 0;
    goto _exit;
  }
  x = mdr_Float_BF(x); // copy and coerce 'x' to double

  e2 = bltin_get_ent (args[1]);
  if (ent_type (e2) != MATRIX_DENSE_REAL)
  {
    mdr_Destroy(x);
    x = 0;
    goto _exit;
  }
  theta = class_double(e2);

  for (i=0; i<MNR(x); i++)
  {
    if (rotate3d(&Mdr0(x,i,0), &Mdr0(x,i,1), &Mdr0(x,i,2), 0.0, 0.0, 1.0, theta))
    {
      mdr_Destroy(x);
      x=0;
      goto _exit;
    }
  }


_exit:

  ent_Clean (e1);
  ent_Clean (e2);
  return ent_Assign_Rlab_MDR(x);
}

#undef   THIS_SOLVER
#define   THIS_SOLVER "bwlabel"
Ent * ent_bwlabel (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0;
  MDR *mask=0, *emask=0;
  int conn=8, min_area=1, max_area, i, j, k;
  int nbwblobs=0, pixel_count;
  ListNode *node=0;

  MDR *npix=0, *npix_new=0, *xy_bbox=0, *xy_bbox_new=0;

  if ( (nargs != 1) && (nargs != 2) )
    goto _exit;

  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) != MATRIX_DENSE_REAL)
    goto _exit;
  emask = ent_data(e1);
  if (SIZE(emask)<1)
    goto _exit;

  mask = mdr_Create_SameSize(emask); // copy and coerce 'mask' to double
  for (i=0; i<MNR(mask); i++)
  {
    for (j=0; j<MNC(mask); j++)
    {
      Mdr0(mask,i,j) = sign(Mdr0(emask,i,j));
    }
  }
  max_area = SIZE(mask);

  if (nargs > 1)
  {
    e2 = bltin_get_ent (args[1]);
    if (ent_type (e2) == MATRIX_DENSE_REAL)
    {
      conn = class_int(e2);
      if (conn != 4 && conn != 8)
        conn = 8;
    }
    else if (ent_type (e2) == BTREE)
    {
      RLABCODE_PROCESS_BTREE_ENTRY_DOUBLE(e2,node,RLAB_NAME_IMG_CONNECTIVITY,conn);
      if (conn!=4 && conn!=8)
        conn = 8;

      RLABCODE_PROCESS_BTREE_ENTRY_DOUBLE(e2,node,RLAB_NAME_SEP_XTR_MINAREA,min_area);
      if (min_area<1)
        min_area = 1;

      RLABCODE_PROCESS_BTREE_ENTRY_DOUBLE(e2,node,RLAB_NAME_SEP_XTR_MAXDEBAREA,max_area);
      if ( (max_area<1) || (max_area>SIZE(mask)) )
        max_area = SIZE(mask);
    }
  }

  // run through entire mask and find blob parameters:
  //    count
  // bbox:
  //    min pos_x
  //    max pos_x
  //    min pos_y
  //    max pos_y
  nbwblobs = mdr_bwlabel(mask, conn, &pixel_count, (MDR**) &npix, (MDR **) &xy_bbox);

//   printf("1: nbwblobs=%i\n", nbwblobs);
//   printf("mask=\n");
//   mdr_Print(mask,stdout);
  if (nbwblobs<1)
  {
//     printf("exiting!\n");
    goto _exit;
  }

  //
//   for (i=0; i<MNR(mask); i++)
//   {
//     for (j=0; j<MNC(mask); j++)
//     {
//       if (Mdr0(mask,i,j) < 1)
//         continue;
//
//       // blob index
//       k = Mdr0(mask,i,j) - 1;
//       if (k >= nbwblobs)
//       {
//         printf("Error occurred in bwlabel: mask[%i,%i] = %i > %i\n",
//                i+1, j+1, (int) k+1, nbwblobs);
//         printf("mask=\n");
//         FILE *fptr=fopen("./debug_bwlabel.txt","w");
//         mdr_Print(mask,fptr);
//         fclose(fptr);
//         rerror("This will crash code!\n");
//       }
//
//       // count
//       MdrV0(npix, k) += 1;
//
//       // bbox: min i, max i, min j, max j:
//       if (MdrV0(npix, k) == 1)
//       {
//         Mdr0(xy_bbox, k, 0) = i+1;
//         Mdr0(xy_bbox, k, 1) = i+1;
//         Mdr0(xy_bbox, k, 2) = j+1;
//         Mdr0(xy_bbox, k, 3) = j+1;
//         continue;
//       }
//
//       // min i
//       if (Mdr0(xy_bbox, k, 0) > i+1)
//         Mdr0(xy_bbox, k, 0) = i+1;
//
//       // max i :
//       if (Mdr0(xy_bbox, k, 1) < i+1)
//         Mdr0(xy_bbox, k, 1) = i+1;
//
//       // min j:
//       if (Mdr0(xy_bbox, k, 2) > j+1)
//         Mdr0(xy_bbox, k, 2) = j+1;
//
//       // max j :
//       if (Mdr0(xy_bbox, k, 3) < j+1)
//         Mdr0(xy_bbox, k, 3) = j+1;
//     }
//   }

  //
  // process min_area and max_area:
  //
//   printf("1: nbwblobs=%i\n", nbwblobs);
  if ( (min_area > 1) || (max_area < SIZE(mask)) )
  {
    int nbwblobs_new=nbwblobs, i0=0, i1=nbwblobs-1;
    // special case of one entry that does not mean min/max criterion:
    if (i0 == i1)
    {
      if ( (mdiV0(npix, i1)  < min_area) || (mdiV0(npix, i1) >  max_area) )
      {
        mdr_Destroy(xy_bbox);
        xy_bbox = 0;
        mdr_Zero(mask);
        mdr_Destroy(npix);
        npix = 0;
        nbwblobs = 0;
      }
    }
    while (i0 < i1)
    {
      if ( (mdiV0(npix, i0) >= min_area) && (mdiV0(npix, i0) <= max_area) )
      {
        i0++;
        if (i0 < i1)
          continue;
      }
      if ( (mdiV0(npix, i1)  < min_area) || (mdiV0(npix, i1) >  max_area) )
      {
        // rewrite its label
        for (i=mdi0(xy_bbox, i1, 0); i<=mdi0(xy_bbox, i1, 1); i++)
        {
          for (j=mdi0(xy_bbox, i1, 2); j<=mdi0(xy_bbox, i1, 3); j++)
          {
            if (Mdr0(mask,i-1,j-1)==i1+1)
              Mdr0(mask,i-1,j-1) = 0;
          }
        }

        // erase old i1:
        MdrV0(npix, i1) = 0;
        // move its bounding box
        Mdr0(xy_bbox, i1, 0) = 0;
        Mdr0(xy_bbox, i1, 1) = 0;
        Mdr0(xy_bbox, i1, 2) = 0;
        Mdr0(xy_bbox, i1, 3) = 0;
        i1--;
        nbwblobs_new--;
//         printf("2: nbwblobs_new=%i\n", nbwblobs_new);
        if (i0 == i1)
        {
          // rewrite its label
          for (i=mdi0(xy_bbox, i1, 0); i<=mdi0(xy_bbox, i1, 1); i++)
          {
            for (j=mdi0(xy_bbox, i1, 2); j<=mdi0(xy_bbox, i1, 3); j++)
            {
              if (Mdr0(mask,i-1,j-1)==i1+1)
                Mdr0(mask,i-1,j-1) = 0;
            }
          }
          i1--;
          nbwblobs_new--;
//           printf("3: nbwblobs_new=%i\n", nbwblobs_new);
        }
        continue;
      }
      if (i0 < i1)
      {
        // move i1 to i0
        for (i=mdr0(xy_bbox, i0, 0); i<=mdr0(xy_bbox, i0, 1); i++)
        {
          for (j=Mdr0(xy_bbox, i0, 2); j<=Mdr0(xy_bbox, i0, 3); j++)
          {
            if (Mdr0(mask,i-1,j-1)==i0+1)
              Mdr0(mask,i-1,j-1) = 0;
          }
        }
        // move its pixel count
        MdrV0(npix, i0) = MdrV0(npix, i1);
        // move its bounding box
        Mdr0(xy_bbox, i0, 0) = Mdr0(xy_bbox, i1, 0);
        Mdr0(xy_bbox, i0, 1) = Mdr0(xy_bbox, i1, 1);
        Mdr0(xy_bbox, i0, 2) = Mdr0(xy_bbox, i1, 2);
        Mdr0(xy_bbox, i0, 3) = Mdr0(xy_bbox, i1, 3);
        // rewrite its label
        for (i=Mdr0(xy_bbox, i0, 0); i<=Mdr0(xy_bbox, i0, 1); i++)
        {
          for (j=Mdr0(xy_bbox, i0, 2); j<=Mdr0(xy_bbox, i0, 3); j++)
          {
            if (Mdr0(mask,i-1,j-1)==i1+1)
              Mdr0(mask,i-1,j-1) = i0+1;
          }
        }
        // erase old i1:
        MdrV0(npix, i1) = 0;
        // move its bounding box
        Mdr0(xy_bbox, i1, 0) = 0;
        Mdr0(xy_bbox, i1, 1) = 0;
        Mdr0(xy_bbox, i1, 2) = 0;
        Mdr0(xy_bbox, i1, 3) = 0;
        i1--;
        nbwblobs_new--;
      }
    } /* while (i0 < i1) */

    if (nbwblobs_new < nbwblobs)
    {
      if (nbwblobs_new>0)
      {
        npix_new = mdr_Create(nbwblobs_new,1);
        xy_bbox_new = mdr_Create(nbwblobs_new,4);
        for (k=0; k<nbwblobs_new; k++)
        {
          MdrV0(npix_new,k) = MdrV0(npix,k);
          Mdr0(xy_bbox_new, k, 0) = Mdr0(xy_bbox, k, 0);
          Mdr0(xy_bbox_new, k, 1) = Mdr0(xy_bbox, k, 1);
          Mdr0(xy_bbox_new, k, 2) = Mdr0(xy_bbox, k, 2);
          Mdr0(xy_bbox_new, k, 3) = Mdr0(xy_bbox, k, 3);
        }
      }
      mdr_Destroy(npix);
      mdr_Destroy(xy_bbox);
      xy_bbox = xy_bbox_new;
      npix = npix_new;
      xy_bbox_new = 0;
      npix_new = 0;
    }
    nbwblobs = nbwblobs_new;
  }

_exit:

  ent_Clean (e1);
  ent_Clean (e2);

  Btree *bw=btree_Create();
  install(bw, RLAB_NAME_IMG_BLOB_COUNT, ent_Create_Rlab_Double(nbwblobs));
  install(bw, RLAB_NAME_SEP_XTR_SEGMAP, ent_Assign_Rlab_MDR(mask));
  install(bw, RLAB_NAME_SEP_XTR_NPIX,   ent_Assign_Rlab_MDR(npix));
  install(bw, RLAB_NAME_BLOBS_BBOX, ent_Assign_Rlab_MDR(xy_bbox));
  return ent_Assign_Rlab_BTREE(bw);
}

#undef  THIS_SOLVER
#define THIS_SOLVER "bwblobs"
Ent * ent_bwblobs (int nargs, Datum args[])
{
  Ent *e0=0, *e1=0, *e2=0;
  MDR *data=0, *mask=0;
  ListNode *node=0;
  int nbwblobs=0;

  // now process return matrix:
  Btree *bw=btree_Create();
  MDR *npix=0, *val_avg=0, *val_std=0, *val_min=0, *val_max=0, *xy_pos=0, *xy_bbox=0, *xy_m2=0;

  if (nargs != 2)
    goto _exit;

  //
  // Get the frame data
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) != MATRIX_DENSE_REAL)
    goto _exit;
  data = ent_data(e1);
  if (SIZE(data)<1)
    goto _exit;

  //
  // Get the result of bwlabel on mask of the frame data:
  //    <<npix;bbox;seg_map>>
  e2 = bltin_get_ent (args[1]);
  if (ent_type (e2) != BTREE)
    goto _exit;

  // npix:
  node = btree_FindNode (ent_data (e2), RLAB_NAME_SEP_XTR_NPIX);
  if (!node)
  {
    printf(THIS_SOLVER ": options entry '"RLAB_NAME_SEP_XTR_NPIX"' is missing!\n");
    goto _exit;
  }
  e0 = var_ent (node);
  if (ent_type(e0)!= MATRIX_DENSE_REAL)
  {
    printf(THIS_SOLVER ": options entry '"RLAB_NAME_SEP_XTR_NPIX"' must be 4-column real matrix!\n");
    goto _exit;
  }
  npix = ent_data(e0);
  nbwblobs = SIZE(npix);
  if (nbwblobs<1)
    goto _exit;

  // bbox:
  node = btree_FindNode (ent_data (e2), RLAB_NAME_BLOBS_BBOX);
  if (!node)
  {
    printf(THIS_SOLVER ": options entry '"RLAB_NAME_BLOBS_BBOX"' is missing!\n");
    goto _exit;
  }
  e0 = var_ent (node);
  if (ent_type(e0)!= MATRIX_DENSE_REAL)
  {
    printf(THIS_SOLVER ": options entry '"RLAB_NAME_BLOBS_BBOX"' must be 4-column real matrix!\n");
    goto _exit;
  }
  xy_bbox = ent_data(e0);
  if (SIZE(xy_bbox)<1)
    goto _exit;
  if (MNR(xy_bbox)!=SIZE(npix))
  {
    printf(THIS_SOLVER ": options entry '"RLAB_NAME_BLOBS_BBOX"' and '"RLAB_NAME_SEP_XTR_NPIX"' must have same number of rows!\n");
    goto _exit;
  }

  // segmap
  node = btree_FindNode (ent_data (e2), RLAB_NAME_SEP_XTR_SEGMAP);
  if (!node)
  {
    printf(THIS_SOLVER ": options entry '"RLAB_NAME_SEP_XTR_SEGMAP"' is missing!\n");
    goto _exit;
  }
  e0 = var_ent (node);
  if (ent_type(e0)!= MATRIX_DENSE_REAL)
  {
    printf(THIS_SOLVER ": options entry "RLAB_NAME_SEP_XTR_SEGMAP" must be real matrix!\n");
    goto _exit;
  }
  mask = ent_data(e0);

  // mask and data have to have same size
  if ( !(EQSIZE(mask,data)) )
  {
    printf(THIS_SOLVER ": options entry "RLAB_NAME_SEP_XTR_SEGMAP" and first argument must have same size!\n");
    goto _exit;
  }

  int i,j,k;
  val_avg = mdr_Create(nbwblobs, 1);
  val_std = mdr_Create(nbwblobs, 1);
  val_min= mdr_Create(nbwblobs, 1);
  val_max = mdr_Create(nbwblobs, 1);
  xy_pos = mdr_Create(nbwblobs, 2);
  xy_m2  = mdr_Create(nbwblobs, 3);

  // go over each blob identified by the bwlabel
  for (k=0; k<nbwblobs; k++)
  {
    // zero everything
    MdrV0(val_avg,k) = 0;
    MdrV0(val_std,k) = 0;
    MdrV0(val_min,k) = create_nan();
    MdrV0(val_max,k) = create_nan();
    Mdr0 (xy_pos, k, 0) = 0;
    Mdr0 (xy_pos, k, 1) = 0;
    Mdr0 (xy_m2,  k, 0) = 0;
    Mdr0 (xy_m2,  k, 1) = 0;
    Mdr0 (xy_m2,  k, 2) = 0;
    // iterate only over the bbox of the mask
    for (i=Mdr0(xy_bbox, k, 0); i<=Mdr0(xy_bbox, k, 1); i++)
    {
      for (j=Mdr0(xy_bbox, k, 2); j<=Mdr0(xy_bbox, k, 3); j++)
      {
        if (Mdr1(mask,i,j)==k+1)
        {
          MdrV0(val_avg,k) += Mdr1(data,i,j);
          Mdr0 (xy_pos, k, 0) += i * Mdr1(data,i,j);
          Mdr0 (xy_pos, k, 1) += j * Mdr1(data,i,j);
          Mdr0 (xy_m2,  k, 0) += i * i * Mdr1(data,i,j);
          Mdr0 (xy_m2,  k, 1) += j * j * Mdr1(data,i,j);
          Mdr0 (xy_m2,  k, 2) += i * j * Mdr1(data,i,j);
          MdrV0(val_std,k) += Mdr1(data,i,j) * Mdr1(data,i,j);
          if (isnand(MdrV0(val_min,k)))
            MdrV0(val_min,k) = Mdr1(data,i,j);
          else if (MdrV0(val_min,k) > Mdr1(data,i,j))
            MdrV0(val_min,k) = Mdr1(data,i,j);
          if (isnand(MdrV0(val_max,k)))
            MdrV0(val_max,k) = Mdr1(data,i,j);
          else if (MdrV0(val_max,k) < Mdr1(data,i,j))
            MdrV0(val_max,k) = Mdr1(data,i,j);
        }

      } /* j */

    } /* i */

    // finish statistics processing
    Mdr0 (xy_pos, k, 0) /= MdrV0(val_avg,k);
    Mdr0 (xy_pos, k, 1) /= MdrV0(val_avg,k);
    Mdr0 (xy_m2,  k, 0) /= MdrV0(val_avg,k);
    Mdr0 (xy_m2,  k, 0) -= Mdr0 (xy_pos, k, 0) * Mdr0 (xy_pos, k, 0);
    Mdr0 (xy_m2,  k, 1) /= MdrV0(val_avg,k);
    Mdr0 (xy_m2,  k, 1) -= Mdr0 (xy_pos, k, 1) * Mdr0 (xy_pos, k, 1);
    Mdr0 (xy_m2,  k, 2) /= MdrV0(val_avg,k);
    Mdr0 (xy_m2,  k, 2) -= Mdr0 (xy_pos, k, 0) * Mdr0 (xy_pos, k, 1);
    MdrV0(val_avg,k) /= MdrV0(npix,k);
    MdrV0(val_std,k) /= MdrV0(npix,k);
    MdrV0(val_std,k) -= MdrV0(val_avg,k) * MdrV0(val_avg,k);

  } /* for (k=0; k<nbwblobs; k++) */

  install(bw, RLAB_NAME_BLOBS_AVG_VALUE,  ent_Assign_Rlab_MDR(val_avg));
  install(bw, RLAB_NAME_BLOBS_STD_VALUE,  ent_Assign_Rlab_MDR(val_std));
  install(bw, RLAB_NAME_BLOBS_MIN_VALUE,  ent_Assign_Rlab_MDR(val_min));
  install(bw, RLAB_NAME_BLOBS_MAX_VALUE,  ent_Assign_Rlab_MDR(val_max));
  install(bw, RLAB_NAME_SEP_XTR_CENTER,   ent_Assign_Rlab_MDR(xy_pos));
  install(bw, RLAB_NAME_SEP_XTR_SECMOMNT, ent_Assign_Rlab_MDR(xy_m2));

_exit:

  ent_Clean (e1);
  ent_Clean (e2);
  return ent_Assign_Rlab_BTREE(bw);
}


