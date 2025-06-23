/* op_mat.c */
/* Matrix class operations for the RLaB machine. */

/*  This file is a part of RLaB ("Our"-LaB)
   Copyright (C) 1992, 1994, 1995  Ian R. Searle

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
#include "class.h"
#include "util.h"
#include "rlab_solver_parameters_names.h"

/* **************************************************************
 * Create a matrix given a Datum as input, return a new Datum that
 * contains the newly created matrix.
 * ************************************************************** */

Datum
matrix_create (Datum d)
{
  return (d);
}

/* **************************************************************
 * Stack a new row onto the bottom of an existing matrix.
 * d1 is the element we will append to ( an existing matrix ).
 * d2 is the element we will append ( the new row ).
 * [ d1 ; d2 ]
 * ************************************************************** */

Datum
matrix_stack (Datum d1, Datum d2)
{
  Datum new;
  Ent *e1=0, *e2=0, *enew=0;

  /* Get rid of any CONSTANTs. */

  e1 = convert_datum_to_ent (d1);
  e2 = convert_datum_to_ent (d2);

  enew = class_stack (e1, e2);
  ent_IncRef (enew);
  new.type = ENTITY;
  new.u.ptr = enew;

  /* Clean up if we can. */

  ent_Clean (e1);
  ent_Clean (e2);

  return (new);
}

/* **************************************************************
 * Create a sub-matrix. Promote a part of an existing matrix.
 *   i_flag = 1: both row and column indices.
 *            2: row index only.
 *            3: column index only.
 *
 *   var '[' i1 ';' i2 ']'
 * ************************************************************** */

Datum
matrix_sub_1 (Datum i1, Datum i2, Datum var)
{
  Datum new;
  Ent *e1=0, *e2=0, *evar=0, *enew=0;

//   fprintf(stderr, "matrix_sub_1 type=%i\n", var.type);


  /* Get rid of any CONSTANTs. */
  e1 = convert_datum_to_ent (i1);
  e2 = convert_datum_to_ent (i2);

  if (var.type == VAR)
  {
    evar = var_ent (var.u.ptr);
  }
  else if (var.type == ENTITY)
  {
    evar = (Ent *) var.u.ptr;
  }
  else
  {
    rerror ("illegal operand (A) for matrix indexing A[ri;ci]");
  }

  enew = class_matrix_sub_1 (evar, e1, e2);

  new.type = ENTITY;
  new.u.ptr = enew;
  ent_IncRef (enew);

  /* Clean up if we can. */
  ent_Clean (e1);
  ent_Clean (e2);
  ent_Clean (evar);

  return (new);
}

Datum
matrix_sub_2 (Datum irow, Datum var)
{
  Datum new;
  Ent *e, *evar, *enew;
  evar = 0;

  /* Get rid of any CONSTANTs. */
//   fprintf(stderr, "matrix_sub_2 type=%i\n", var.type);

  e = convert_datum_to_ent (irow);

  if (var.type == VAR)
  {
    evar = var_ent (var.u.ptr);
  }
  else if (var.type == ENTITY)
  {
    evar = (Ent *) var.u.ptr;
//     ent_DecRef (evar);
  }
  else
  {
    rerror ("illegal operand (A) for matrix indexing A[ri;]");
  }

  enew = class_matrix_sub_2 (evar, e);
  ent_IncRef (enew);
  new.type = ENTITY;
  new.u.ptr = enew;

  /* Clean up if we can. */

  ent_Clean (e);
  ent_Clean (evar);

  return (new);
}

Datum
matrix_sub_3 (Datum icol, Datum var)
{
//   fprintf(stderr, "matrix_sub_3 type=%i\n", var.type);

  Datum new;
  Ent *e, *evar, *enew;
  evar = 0;

  /* Get rid of any CONSTANTs. */

  e = convert_datum_to_ent (icol);

  if (var.type == VAR)
  {
    evar = var_ent (var.u.ptr);
  }
  else if (var.type == ENTITY)
  {
    evar = (Ent *) var.u.ptr;
  }
  else
  {
    rerror ("illegal operand (A) for matrix indexing A[;ci]");
  }

  enew = class_matrix_sub_3 (evar, e);
  ent_IncRef (enew);
  new.type = ENTITY;
  new.u.ptr = enew;

  /* Clean up if we can. */

  ent_Clean (e);
  ent_Clean (evar);
  return (new);
}

/*
 * Allow indexing of a MATRIX or SCALAR like a VECTOR
 * var [ ind ]
 */

Datum
matrix_vector_sub (Datum var, Datum ind)
{
  Datum new;
  Ent *eind, *evar, *enew;
  evar = 0;

  /* Get rid of any CONSTANTs. */
//   fprintf(stderr, "matrix_vector_sub type=%i\n", var.type);

  eind = convert_datum_to_ent (ind);

  if (var.type == VAR)
  {
    evar = var_ent (var.u.ptr);
  }
  else if (var.type == ENTITY)
  {
    evar = (Ent *) var.u.ptr;
  }
  else
  {
    rerror ("illegal operand (A) for matrix indexing A[index]");
  }

  enew = class_vector_sub (evar, eind);

  new.type = ENTITY;
  new.u.ptr = enew;
  ent_IncRef (enew);

  /* Clean up if we can. */
  ent_Clean (eind);
  ent_Clean (evar);

  return (new);
}

/* **************************************************************
 * Assign new values to the element(s) of a matrix.
 * i1:  row index(s).
 * i2:  col index(s).
 * a:   the new value(s).
 * var: The matrix that will be modified.
 *
 *      var [ i1; i2 ] = a
 *
 * ************************************************************** */
#undef THIS_SOLVER
#define THIS_SOLVER "matrix_assign_1"
Datum matrix_assign_1 (Datum var, Datum i1, Datum i2, Datum a)
{
  Datum new;
  Ent *ea, *e1, *e2, *evar;
  ListNode *lvar;
  evar = 0;
  lvar = 0;

  /* Get rid of any CONSTANTs. */

  e1 = convert_datum_to_ent (i1);
  e2 = convert_datum_to_ent (i2);
  ea = convert_datum_to_ent (a);

//   fprintf(stderr, "matrix_assign_1 type=%i\n", var.type);

  if (var.type == VAR)
  {
    lvar = (ListNode *) (var.u.ptr);
    evar = var_ent (lvar);

    //
    // Since we are going to assign to VAR, we
    // must copy it if necessary.
    // 
    if (evar->refc > 1)
    {
      evar = ent_Create ();
      evar = class_copy (var_ent (lvar));
      ent_DecRef (var_ent (lvar));
      ent_IncRef (evar);
    }

    //
    // Now do the assign.
    // 
    evar = class_matrix_assign_1 (evar, e1, e2, ea);

    // Re-attach the entity.
    listNode_AttachEnt (lvar, evar);

    // Set up the return Datum.
    new.u.ptr = lvar;
    new.type = VAR;
  }
  else if (var.type == ENTITY)
  {
    evar = (Ent *) var.u.ptr;
    if (evar->refc > 1)
    {
      ent_DecRef (evar);
    }
    evar = class_matrix_assign_1 (evar, e1, e2, ea);
    new.u.ptr = evar;
    new.type  = ENTITY;
    ent_IncRef (evar);
  }
  else
    rerror (THIS_SOLVER ": " RLAB_ERROR_OPMAT_INVALID_OBJECT_ASSIGN "\n");

  // Clean up if we can.
  ent_Clean (e1);
  ent_Clean (e2);
  ent_Clean (ea);

  return (new);
}

#undef THIS_SOLVER
#define THIS_SOLVER "matrix_assign_2"
Datum matrix_assign_2 (Datum var, Datum i1, Datum a)
{
  Datum new;
  Ent *ea, *e1, *evar;
  ListNode *lvar;
  evar = 0;
  lvar = 0;

  /* Get rid of any CONSTANTs. */

  e1 = convert_datum_to_ent (i1);
  ea = convert_datum_to_ent (a);

//   fprintf(stderr, "matrix_assign_2 type=%i\n", var.type);

  if (var.type == VAR)
  {
    lvar = (ListNode *) (var.u.ptr);
    evar = var_ent (lvar);

    //
    // Since we are going to assign to VAR, we
    // must copy it if necessary.
    //
    if (evar->refc > 1)
    {
      evar = ent_Create ();
      evar = class_copy (var_ent (lvar));
      ent_DecRef (var_ent (lvar));
      ent_IncRef (evar);
    }

    //
    // Now do the assign.
    // 
    evar = class_matrix_assign_2 (evar, e1, ea);

    // Re-attach the entity.
    listNode_AttachEnt (lvar, evar);

    // Set up the return Datum.
    new.u.ptr = lvar;
    new.type = VAR;
  }
  else if (var.type == ENTITY)
  {
    evar = (Ent *) var.u.ptr;
    if (evar->refc > 1)
    {
      ent_DecRef (evar);
    }
    evar = class_matrix_assign_2 (evar, e1, ea);

//     ent_DecRef (evar);
    new.u.ptr = evar;
    new.type  = ENTITY;
    ent_IncRef (evar);

  }
  else
    rerror (THIS_SOLVER ": " RLAB_ERROR_OPMAT_INVALID_OBJECT_ASSIGN "\n");

  /* Clean up if we can. */
  ent_Clean (e1);
  ent_Clean (ea);

  return (new);
}

#undef THIS_SOLVER
#define THIS_SOLVER "matrix_assign_3"
Datum
matrix_assign_3 (Datum var, Datum i2, Datum a)
{
  Datum new;
  Ent *ea, *e2, *evar;
  ListNode *lvar;
  evar = 0;
  lvar = 0;

  /* Get rid of any CONSTANTs. */

  e2 = convert_datum_to_ent (i2);
  ea = convert_datum_to_ent (a);

//   fprintf(stderr, "matrix_assign_3 type=%i\n", var.type);

  if (var.type == VAR)
  {
    lvar = (ListNode *) (var.u.ptr);
    evar = var_ent (lvar);

    //
    // Since we are going to assign to VAR, we
    // must copy it if necessary.
    // 
    if (evar->refc > 1)
    {
      evar = ent_Create ();
      evar = class_copy (var_ent (lvar));
      ent_DecRef (var_ent (lvar));
      ent_IncRef (evar);
    }

    //
    // Now do the assign.
    // 
    evar = class_matrix_assign_3 (evar, e2, ea);

    // Re-attach the entity.
    listNode_AttachEnt (lvar, evar);

    // Set up the return Datum
    new.u.ptr = lvar;
    new.type = VAR;
  }
  else if (var.type == ENTITY)
  {
    evar = (Ent *) var.u.ptr;
    if (evar->refc > 1)
    {
      ent_DecRef (evar);
    }
    evar = class_matrix_assign_3 (evar, e2, ea);
//     ent_DecRef (evar);
    new.u.ptr = evar;
    new.type  = ENTITY;
    ent_IncRef (evar);
  }
  else
    rerror (THIS_SOLVER ": " RLAB_ERROR_OPMAT_INVALID_OBJECT_ASSIGN "\n");

  /* Clean up if we can. */
  ent_Clean (e2);
  ent_Clean (ea);

  return (new);
}

/* **************************************************************
 * Handle assignments where the matrix is treated like a vector.
 * ************************************************************** */
#undef THIS_SOLVER
#define THIS_SOLVER "matrix_vector_assign"
Datum matrix_vector_assign (Datum var, Datum i, Datum rhs)
{
  Datum new;
  Ent *evar, *ei, *erhs;
  ListNode *lvar;
  evar = 0;
  lvar = 0;

  /* Get rid of any CONSTANTs. */

  ei = convert_datum_to_ent (i);
  erhs = convert_datum_to_ent (rhs);

//   fprintf(stderr, "matrix_vector_assign type=%i\n", var.type);

  if (var.type == VAR)
  {
    lvar = (ListNode *) (var.u.ptr);
    evar = var_ent (lvar);

    //
    // Since we are going to assign to VAR, we
    // must copy it if necessary.
    //
    if (evar->refc > 1)
    {
      evar = ent_Create ();
      evar = class_copy (var_ent (lvar));
      ent_DecRef (var_ent (lvar));
      ent_IncRef (evar);
    }

    //
    // Now do the assign.
    //
    evar = class_matrix_vector_assign (evar, ei, erhs);

    // Re-attach the entity
    listNode_AttachEnt (lvar, evar);

    // Set up the return Datum.
    new.u.ptr = lvar;
    new.type = VAR;
  }
  else if (var.type == ENTITY)
  {
    evar = (Ent *) var.u.ptr;
    if (evar->refc > 1)
    {
      ent_DecRef (evar);
    }
    evar = class_matrix_vector_assign (evar, ei, erhs);

    new.u.ptr = evar;
    new.type  = ENTITY;
    ent_IncRef (evar);
//     ent_DecRef (evar);
  }
  else
    rerror (THIS_SOLVER ": " RLAB_ERROR_OPMAT_INVALID_OBJECT_ASSIGN "\n");

  // Clean up if we can
  ent_Clean (ei);
  ent_Clean (erhs);

  return (new);
}

/* **************************************************************
 * Create and return the transpose of a matrix.
 * ************************************************************** */

Datum matrix_transpose (Datum d)
{
  Ent *e, *enew;

  e = convert_datum_to_ent (d);

  enew = class_transpose (e);
  d.u.ptr = enew;
  d.type = ENTITY;
  ent_IncRef (enew);

  ent_Clean (e);
  return (d);
}

/*
 * Matrix non-conjugate transpose
 */

Datum matrix_el_transpose (Datum d)
{
  Ent *e, *enew;

  e = convert_datum_to_ent (d);

  enew = class_nc_transpose (e);
  d.u.ptr = enew;
  d.type = ENTITY;
  ent_IncRef (enew);

  ent_Clean (e);
  return (d);
}

/* **************************************************************
 * Reshape a arbitrarily sized matrix into a column matrix
 * ************************************************************** */

Datum
matrix_reshape_col (Datum d)
{
  Ent *e, *enew;

  e = convert_datum_to_ent (d);

  enew = class_reshape_col (e);
  d.u.ptr = enew;
  d.type = ENTITY;
  ent_IncRef (enew);

  ent_Clean (e);
  return (d);
}
