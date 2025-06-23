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
#include "msc.h"
#include "mscf1.h"
#include "mscf2.h"

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
