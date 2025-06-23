/* simplest.c: A simple builtin function example. */

/*
   Compile this file with (in this directory):
   cc -fPIC -c simplest.c -I../../ -I../../gc
*/

/* Necessary header files. */
#include "rlab.h"
#include "ent.h"
#include "class.h"
#include "mem.h"
#include "bltin.h"
#include "util.h"

#include <stdio.h>
#include <string.h>

Ent *
Simplest (int nargs, Datum args[])
{
  double dtmp;
  Ent *e, *rent;
  MDR *m;

  /* Check the number of arguments. */
  if (nargs != 1)
  {
    rerror ("simplest: only 1 argument allowed");
  }

  /* Get the first (only) argument. */
  e = bltin_get_ent (args[0]);

  /* Perform the simplest operation. */
  dtmp = 2.0 * class_double (e);

  /* Create a new matrix containing the result. */
  m = mdr_CreateScalar (dtmp);

  /* Create the return entity. */
  rent = ent_Create ();

  /* Set the entity's data. */
  ent_data (rent) = m;

  /* Set the entity's data type. */
  ent_type (rent) = MATRIX_DENSE_REAL;

  /* Clean up the argument if possible. */
  ent_Clean (e);

  return (rent);
}
