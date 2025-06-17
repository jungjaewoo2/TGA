/* op.c */

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
#include "listnode.h"
#include "class.h"
#include "util.h"

/* bltin1.c */
extern double errno_check _PROTO ((double d, char *s));

#define STYPE(l, r)    (10*(l) + (r))

#define DOUBLE_DOUBLE     (10*DOUBLE + DOUBLE)
#define DOUBLE_iDOUBLE    (10*DOUBLE + iDOUBLE)
#define iDOUBLE_DOUBLE    (10*iDOUBLE + DOUBLE)



/* **************************************************************
 * General addition operation on two Datums, return a new Datum.
 * ************************************************************** */
#undef  THIS_FILE
#define THIS_FILE "op.c"

#undef  THIS_SOLVER
#define THIS_SOLVER "add"
Datum addition_op (Datum d1, Datum d2, int reuse_d1_ent)
{
  Datum new;
  Ent *e1=convert_datum_to_ent (d1), *e2=convert_datum_to_ent (d2);
  Ent *enew=0;

  new.type = ENTITY;
  if (reuse_d1_ent)
  {
    if (ent_Ref(e1) == 1)
    {
      // nobody uses d1 data. we can safely rewrite it
      enew = e1;
      ent_IncRef (e1);
      class_addto (enew, e2);       // enew -> enew + e2
    }
    else
    {
      // d1 data is used by other vars. we must duplicate it, but then what is the point
      enew = class_add (e1, e2);  // enew = e1 + e2
    }
  }
  else
  {
    enew = class_add (e1, e2);  // enew = e1 + e2
  }

  new.u.ptr = enew;
  ent_IncRef (enew);
  ent_Clean (e1);
  ent_Clean (e2);
  return (new);
}

/* **************************************************************
 * General subtraction operation on two Datums, return a new Datum.
 * d1 - d2
 * ************************************************************** */

Datum
subtraction_op (d1, d2)
     Datum d1, d2;
{
  Datum new;
  Ent *e1, *e2, *enew;

  /*
   * Make sure each Datum contains an entity.
   */

  e1 = convert_datum_to_ent (d1);
  e2 = convert_datum_to_ent (d2);

  /*
   * Let the class controller handle the operation.
   */

  enew = class_subtract (e1, e2);

  /*
   * Now fix up the entity to push back on the stack.
   */

  new.type = ENTITY;
  new.u.ptr = enew;
  ent_IncRef (enew);

  /* Clean up if we can. */

  ent_Clean (e1);
  ent_Clean (e2);

  return (new);
}

/* **************************************************************
 * General multiply operation on two Datums, return a new Datum.
 * ************************************************************** */

Datum
multiply_op (d1, d2)
     Datum d1, d2;
{
  Datum new;
  Ent *e1, *e2, *enew;

  /*
   * Make sure each Datum contains an entity.
   */

  e1 = convert_datum_to_ent (d1);
  e2 = convert_datum_to_ent (d2);

  /*
   * Let the class controller handle the operation.
   */

  enew = class_multiply (e1, e2);

  /*
   * Now fix up the entity to push back on the stack.
   */

  new.type = ENTITY;
  new.u.ptr = enew;
  ent_IncRef (enew);

  /* Clean up if we can. */

  ent_Clean (e1);
  ent_Clean (e2);

  return (new);
}

/* **************************************************************
 * Element-by-element multiply. 
 * ************************************************************** */

Datum
element_multiply (Datum d1, Datum d2)
{
  Datum new;
  Ent *e1, *e2, *enew;

  /*
   * Make sure each Datum contains an entity.
   */

  e1 = convert_datum_to_ent (d1);
  e2 = convert_datum_to_ent (d2);

  /*
   * Let the class controller handle the operation.
   */

  enew = class_el_multiply (e1, e2);

  /*
   * Now fix up the entity to push back on the stack.
   */

  new.type = ENTITY;
  new.u.ptr = enew;
  ent_IncRef (enew);

  /* Clean up if we can. */

  ent_Clean (e1);
  ent_Clean (e2);

  return (new);
}

/* **************************************************************
 * Right Division  A / B
 * ************************************************************** */

Datum
rdivide (d1, d2)
     Datum d1, d2;
{
  Datum new;
  Ent *e1, *e2, *enew;

  /*
   * Make sure each Datum contains an entity.
   */

  e1 = convert_datum_to_ent (d1);
  e2 = convert_datum_to_ent (d2);

  /*
   * Let the class controller handle the operation.
   */

  enew = class_rdivide (e1, e2);

  /*
   * Now fix up the entity to push back on the stack.
   */

  new.type = ENTITY;
  new.u.ptr = enew;
  ent_IncRef (enew);

  /* Clean up if we can. */

  ent_Clean (e1);
  ent_Clean (e2);

  return (new);
}

/*
 * Element - by - element right divide.
 */

Datum
el_rdivide (d1, d2)
     Datum d1, d2;
{
  Datum new;
  Ent *e1, *e2, *enew;

  /*
   * Make sure each Datum contains an entity.
   */

  e1 = convert_datum_to_ent (d1);
  e2 = convert_datum_to_ent (d2);

  /*
   * Let the class controller handle the operation.
   */

  enew = class_el_rdivide (e1, e2);

  /*
   * Now fix up the entity to push back on the stack.
   */

  new.type = ENTITY;
  new.u.ptr = enew;
  ent_IncRef (enew);

  /* Clean up if we can. */

  ent_Clean (e1);
  ent_Clean (e2);

  return (new);
}

/* **************************************************************
 * Left Division  A \ B
 * ************************************************************** */

Datum
ldivide (d1, d2)
     Datum d1, d2;
{
  Datum new;
  Ent *e1, *e2, *enew;

  /*
   * Make sure each Datum contains an entity.
   */

  e1 = convert_datum_to_ent (d1);
  e2 = convert_datum_to_ent (d2);

  /*
   * Let the class controller handle the operation.
   */

  enew = class_ldivide (e1, e2);

  /*
   * Now fix up the entity to push back on the stack.
   */

  new.type = ENTITY;
  new.u.ptr = enew;
  ent_IncRef (enew);

  /* Clean up if we can. */

  ent_Clean (e1);
  ent_Clean (e2);

  return (new);
}

/*
 * Element - by element left divide.
 */

Datum
el_ldivide (d1, d2)
     Datum d1, d2;
{
  Datum new;
  Ent *e1, *e2, *enew;

  /*
   * Make sure each Datum contains an entity.
   */

  e1 = convert_datum_to_ent (d1);
  e2 = convert_datum_to_ent (d2);

  /*
   * Let the class controller handle the operation.
   */

  enew = class_el_ldivide (e1, e2);

  /*
   * Now fix up the entity to push back on the stack.
   */

  new.type = ENTITY;
  new.u.ptr = enew;
  ent_IncRef (enew);

  /* Clean up if we can. */

  ent_Clean (e1);
  ent_Clean (e2);

  return (new);
}

/* **************************************************************
 * Negate a CONSTANT or an ENTITY.
 * ************************************************************** */

Datum
negate_op (Datum d)
{
  Ent *e, *enew;

  e = convert_datum_to_ent (d);

  enew = class_negate (e);
  d.u.ptr = enew;
  d.type = ENTITY;
  ent_IncRef (enew);
  return (d);
}

/* **************************************************************
 * Return a Datum that contains d1 ^ d2
 * ************************************************************** */

Datum
power_op (d1, d2)
     Datum d1, d2;
{
  Datum new;
  Ent *e1, *e2, *enew;

  /*
   * Make sure each Datum contains an entity.
   */

  e1 = convert_datum_to_ent (d1);
  e2 = convert_datum_to_ent (d2);

  /*
   * Let the class controller handle the operation.
   */

  enew = class_power (e1, e2);

  /*
   * Now fix up the entity to push back on the stack.
   */

  new.type = ENTITY;
  new.u.ptr = enew;
  ent_IncRef (enew);

  /* Clean up if we can. */

  ent_Clean (e1);
  ent_Clean (e2);

  return (new);
}

/*
 * Element - by - element Power Operation.
 * d1 .^ d2
 */

Datum
el_power (d1, d2)
     Datum d1, d2;
{
  Datum new;
  Ent *e1, *e2, *enew;

  /*
   * Make sure each Datum contains an entity.
   */

  e1 = convert_datum_to_ent (d1);
  e2 = convert_datum_to_ent (d2);

  /*
   * Let the class controller handle the operation.
   */

  enew = class_el_power (e1, e2);

  /*
   * Now fix up the entity to push back on the stack.
   */

  new.type = ENTITY;
  new.u.ptr = enew;
  ent_IncRef (enew);

  /* Clean up if we can. */

  ent_Clean (e1);
  ent_Clean (e2);

  return (new);
}

/* **************************************************************
 * Create an empty matrix. This is difficult, cause we don't know
 * what type of empty matrix to create. Lets just assume it is a
 * dense real matrix (MDR), and rely upon coercion to handle the
 * rest...
 * ************************************************************** */

Datum
empty_create (void)
{
  Datum new;
  Ent *enew;

  enew = class_empty ();

  /*
   * Now fix up the entity to push back on the stack.
   */

  new.type = ENTITY;
  new.u.ptr = enew;
  ent_IncRef (enew);

  return (new);
}
