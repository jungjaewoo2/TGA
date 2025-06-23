/* opapp.c */
/* Append related operations for the RLaB machine. */

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
#include "class.h"
#include "listnode.h"
#include "util.h"
#include "mdr.h"

/* **************************************************************
 * Create a vector from the command:
 *     expr ':' expr ':' expr
 *     start    end      increment
 *      d3       d2       d1         (if n = 3)
 * or   d2       d1                  (if n = 2)
 * Where the last expr is optional (default is 1)

 * For now (maybe for-ever), this function is hardwired...
 * ************************************************************** */

Datum
vector_create_op (int n, Datum d1, Datum d2, Datum d3)
{
  Datum new;
  Ent *e1, *e2, *e3, *er;

  er = ent_Create ();

  e1 = convert_datum_to_ent (d1);
  e2 = convert_datum_to_ent (d2);

  if (n == 3)
  {
    e3 = convert_datum_to_ent (d3);
    ent_data (er) = mdr_vector_create (ent_data (e3),
				       ent_data (e2), ent_data (e1));
    ent_Clean (e1);
    ent_Clean (e2);
    ent_Clean (e3);
  }
  else
  {
    ent_data (er) = mdr_vector_create (ent_data (e2), ent_data (e1), 0);
    ent_Clean (e1);
    ent_Clean (e2);
  }

  ent_type (er) = MATRIX_DENSE_REAL;
  ent_IncRef (er);
  new.type = ENTITY;
  new.u.ptr = er;

  return (new);
}

/* **************************************************************
 * Create a MATRIX. To do this in the most general way we must
 * build the MATRIX by appending elements to it.
 * d1 , d2
 * d1 is the root element.
 * d2 is the appendee (element we are adding).
 * ************************************************************** */

Datum
vector_append (Datum d1, Datum d2)
{
  Datum new;
  Ent *e1, *e2, *enew;

  /* Get rid of any CONSTANTs. */

  e1 = convert_datum_to_ent (d1);
  e2 = convert_datum_to_ent (d2);

  enew = class_append (e1, e2);
  ent_IncRef (enew);
  new.type = ENTITY;
  new.u.ptr = enew;

  /* Clean up if we can. */

  ent_Clean (e1);
  ent_Clean (e2);

  return (new);
}
