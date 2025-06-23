/* relation.c
   Take care of the object relational comparisons */

/*  This file is a part of RLaB ("Our"-LaB)
   Copyright (C) 1992, 1994  Ian R. Searle

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

/* **************************************************************
 * equal-to operation.
 * d1 == d2
 * ************************************************************** */

Datum
eq (d1, d2)
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

  enew = class_eq (e1, e2);

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
 * not-equal-to operation.
 * ************************************************************** */

Datum
ne (d1, d2)
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

  enew = class_ne (e1, e2);

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

Datum
le (d1, d2)
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

  enew = class_le (e1, e2);

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

Datum
ge (d1, d2)
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

  enew = class_ge (e1, e2);

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

Datum
lt (d1, d2)
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

  enew = class_lt (e1, e2);

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

Datum
gt (d1, d2)
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

  enew = class_gt (e1, e2);

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

Datum
and (d1, d2)
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

  enew = class_and (e1, e2);

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

Datum
or (d1, d2)
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

  enew = class_or (e1, e2);

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

Datum
not (d)
     Datum d;
{
  Datum new;
  Ent *e, *enew;

  /*
   * Make sure each Datum contains an entity.
   */

  e = convert_datum_to_ent (d);

  /*
   * Let the class controller handle the operation.
   */

  enew = class_not (e);

  /*
   * Now fix up the entity to push back on the stack.
   */

  new.type = ENTITY;
  new.u.ptr = enew;
  ent_IncRef (enew);

  /* Clean up if we can. */

  ent_Clean (e);

  return (new);
}
