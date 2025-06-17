/* **************************************************************
 * pgm.c: Routines to load and save a matrix containing the pixels
 *        of a gray scale image to or from a Portable Gray Map
 *        file, using the libpgm library from the netpbm
 *        distribution.
 *
 * To compile on a Linux/ELF system:
 *        gcc -g -fPIC -c pgm.c -I../../ -I../../gc
 *        gcc -shared -o pgm.so pgm.o -lpgm -lpbm
 * ************************************************************** */

#include "rlab.h"
#include "ent.h"
#include "class.h"
#include "mem.h"
#include "bltin.h"
#include "rfileio.h"
#include "util.h"

#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <math.h>

#include "pgm.h"

/* **************************************************************
 * RSavepgm: Save a PGM matrix to a file.

 * savepgm ( img_matrix, file_name, maximum_gray_level )
 * ************************************************************** */

Ent *
RSavepgm (int nargs, Datum args[])
{
  int gl, max_gl = 0, i, j;
  char *string;
  double dgl;
  FILE *fn;
  Ent *FN, *GL, *IMG, *rent;
  MDR *img;
  gray **new;
  gray *row_new;

  char *kluge_argv[1];	   /* we need to provide the pgm lib a dummy argv */

  kluge_argv[0]="savepgm"; /* initialize the dummy argv */

  /* Check nargs */
  if ((nargs < 2) || (nargs > 3))
    rerror ("savepgm: requires 2 or 3 arguments");

  /* Get the image. */
  /* First the image entity. */
  IMG = bltin_get_ent (args[0]);

  /* Next, get the image matrix from within the entity. */
  img = class_matrix_real (IMG);

  /* Then the filename for output. */
  FN = bltin_get_ent (args[1]);
  string = class_char_pointer (FN);

  /* If the third argument is present, get it for use as the maximum gray level */
  gl = -1;
  if (nargs == 3)
  {
    GL = bltin_get_ent (args[2]);
    dgl = class_double (GL);
    gl = dgl;
  }

  /* Open with file for binary write. */
  if ((fn = get_file_ds (string, "wb", 0)) == 0)
  {
    fprintf (stderr, "savepgm: %s: cannot open for write\n", string);

    /* Clean up the arguments when we error out. */
    ent_Clean (IMG);
    ent_Clean (FN);
    if (nargs == 3) ent_Clean (GL);

    /* Return 0 to indicate failure. */
    rent = ent_Create ();
    ent_data (rent) = mdr_CreateScalar (0.0);
    ent_type (rent) = MATRIX_DENSE_REAL;
    return (rent);
  }

  /*
   * First we need to call pgm_init to initialize the pgm library.
   * Normally this is called with argc and argv, but here we want to
   * just dummy it up.
   */
  i=1;
  pgm_init (&i, kluge_argv);

  /* Allocate a PGM image array of the correct size */
  new = pgm_allocarray (MNC (img), MNR (img));

  /*
   * Now for each row of the image we want to store the pixel values
   * for each column.  Of course PGM differs from RLaB in the choice
   * of column-major and row-major order.
   */

  for (j = 0; j < MNR (img);j++)
  {
    row_new = *(new+j);
    for (i = 0; i < MNC (img); i++)
    {
      *(row_new+i) = (gray) MdrV0 (img, i*MNR(img)+j);

      /* Keep track of the maximum pixel value in the image */
      if(*(row_new+i) > max_gl)
      {
	max_gl=*(row_new+i);
      }
    }
  }

  /*
   * If no maximum gray level was given as an argument, use the maximum
   * pixel value detected above. If the detected maximum pixel value is
   * greater than the one specified in argument 3, give a warning, and use
   * the maximum detected value.
   */

  if(gl == -1)
  {
    gl = max_gl;
  }
  else if(max_gl > gl)
  {
    fprintf (stderr,
	     "savepgm: image contains pixel values greater than specified maximum");
    fprintf (stderr, "\nusing maximum pixel value instead\n");
    gl = max_gl;
  }

  /* Now the array new contains the PGM image, so write it out */
  pgm_writepgm (fn, new, MNC (img), MNR (img),(gray)gl, 0);
  pgm_freearray (new, MNR (img));

  /* Clean up before returning. */
  ent_Clean (FN);
  ent_Clean (IMG);
  if (nargs == 3) ent_Clean (GL);
  close_file_ds (string);

  /* Everything OK, return 1 to indicate success. */
  rent = ent_Create ();
  ent_data (rent) = mdr_CreateScalar (1.0);
  ent_type (rent) = MATRIX_DENSE_REAL;
  return (rent);
}

/* **************************************************************
 * RLavepgm: Load a PGM into a matrix.

 * loadpgm ( file_name )
 * ************************************************************** */

Ent *
RLoadpgm (int nargs, Datum args[])
{
  int i, j, rows, cols;
  char *string;
  FILE *fn;
  Ent *FN, *rent;
  MDR *img;
  gray **new;
  gray *row_new;
  gray gl;
  char *kluge_argv[1];     /* we need to provide the pgm lib a dummy argv */

  kluge_argv[0]="savepgm"; /* initialize the dummy argv */

  /* Check nargs */
  if (nargs != 1)
    rerror ("loadpgm: requires 1 argument");

  /* The the filename for input. */
  FN = bltin_get_ent (args[0]);
  string = class_char_pointer (FN);

  /* Open with file for binary read. */
  if ((fn = get_file_ds (string, "rb", 0)) == 0)
  {
    fprintf (stderr, "loadpgm: %s: cannot open for write\n", string);

    /* Clean up the arguments when we error out. */
    ent_Clean (FN);

    /* Return 0 to indicate failure. */
    rent = ent_Create ();
    ent_data (rent) = mdr_CreateScalar (0.0);
    ent_type (rent) = MATRIX_DENSE_REAL;
    return (rent);
  }

  /*
   * First we need to call pgm_init to initialize the pgm library.
   * Normally this is called with argc and argv, but here we want to
   * just dummy it up.
   */

  i = 1;
  pgm_init (&i, kluge_argv);

  /* Allocate a PGM image array of the correct size */
  new = pgm_readpgm (fn, &cols, &rows, &gl);
  img = mdr_Create (rows, cols);

  /*
   * Now for each row of the image we want to store the pixel values
   * for each column.  Of course PGM differs from RLaB in the choice
   * of column-major and row-major order.
   */
  for (j = 0; j < rows;j++)
  {
    row_new = *(new+j);
    for (i = 0; i < cols; i++)
    {
       MdrV0 (img, i*MNR(img)+j) = *(row_new+i);
    }
  }

  /* Clean up before returning. */
  pgm_freearray(new, MNR (img));
  ent_Clean (FN);
  close_file_ds (string);

  /* Everything OK, return the image. */
  rent = ent_Create ();
  ent_data (rent) = img;
  ent_type (rent) = MATRIX_DENSE_REAL;
  return (rent);
}
