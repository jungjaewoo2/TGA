/* r_pgplot.c */

#include "config.h"

#include "rlab.h"
#include "symbol.h"
#include "bltin.h"
#include "listnode.h"
#include "mem.h"
#include "mdr.h"
#include "mdc.h"
#include "mds.h"
#include "ent.h"
#include "class.h"
#include "util.h"
#include "rfileio.h"

#include <stdio.h>
#include "cpgplot.h"

float *mdr_coerce_floatptr (MDR * m);

/* **************************************************************
 * Open a graphics device.
 * INTEGER FUNCTION PGBEG (UNIT, FILE, NXSUB, NYSUB)
 *     INTEGER       UNIT
 *     CHARACTER*(*) FILE
 *     INTEGER       NXSUB, NYSUB
 * ************************************************************** */

Ent *
_pg_cpgbeg (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *e3=0, *e4=0;
  int unit, nxsub, nysub;
  char *file;

  if (nargs != 4)
  {
    rerror ("pgbeg: 4 arguments are required");
  }

  e1 = bltin_get_ent (args[0]);
  e2 = bltin_get_ent (args[1]);
  e3 = bltin_get_ent (args[2]);
  e4 = bltin_get_ent (args[3]);

  unit = (int) class_double (e1);
  file = class_char_pointer (e2);
  nxsub = (int) class_double (e3);
  nysub = (int) class_double (e4);

  cpgbeg (unit, file, nxsub, nysub);

  ent_Clean (e1);
  ent_Clean (e2);
  ent_Clean (e3);
  ent_Clean (e4);

  return ent_Create_Rlab_Success();
}

/* **************************************************************
 * Open a graphics device.
 *     INTEGER FUNCTION PGOPEN (DEVICE)
 *         CHARACTER*(*) DEVICE
 * ************************************************************** */

Ent *
_pg_cpgopen (int nargs, Datum args[])
{
  Ent *e1=0;
  char *device;
  int istat;

  if (nargs != 1)
  {
    rerror ("pgopen: 1 argument allowed");
  }

  e1 = bltin_get_ent (args[0]);

  device = class_char_pointer (e1);

  istat = cpgopen (device);
  if (istat <= 0)
  {
    rerror ("pgopen: error opening PGPLOT plot device");
  }

  ent_Clean (e1);

  return ent_Create_Rlab_Double((double) istat);
}

/* **************************************************************
 * set viewport (normalized device coordinates)
 *     SUBROUTINE PGSVP (XLEFT, XRIGHT, YBOT, YTOP)
 *         REAL XLEFT, XRIGHT, YBOT, YTOP
 * ************************************************************** */

Ent *
_pg_cpgsvp (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *e3=0, *e4=0;
  float xleft, xright, ybot, ytop;

  if (nargs != 4)
  {
    rerror ("pgsvp: 4 arguments required");
  }

  e1 = bltin_get_ent (args[0]);
  e2 = bltin_get_ent (args[1]);
  e3 = bltin_get_ent (args[2]);
  e4 = bltin_get_ent (args[3]);

  xleft = (float) class_double (e1);
  xright = (float) class_double (e2);
  ybot = (float) class_double (e3);
  ytop = (float) class_double (e4);

  cpgsvp (xleft, xright, ybot, ytop);

  ent_Clean (e1);
  ent_Clean (e2);
  ent_Clean (e3);
  ent_Clean (e4);

  return ent_Create_Rlab_Success();
}

/* **************************************************************
 * set standard (default) viewport
 *     SUBROUTINE PGVSTD
 * ************************************************************** */

Ent *
_pg_cpgvstd (int nargs, Datum args[])
{
  if (nargs != 0)
    rerror ("pgvstd: no arguments allowed");

  cpgvstd ();

  return ent_Create_Rlab_Success();
}

/* **************************************************************
 * set window
 *    SUBROUTINE PGSWIN (X1, X2, Y1, Y2)
 *       REAL X1, X2, Y1, Y2
 * ************************************************************** */

Ent *
_pg_cpgswin (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *e3=0, *e4=0;
  float x1, x2, y1, y2;

  if (nargs != 4)
  {
    rerror ("pgswin: 4 arguments required");
  }

  e1 = bltin_get_ent (args[0]);
  e2 = bltin_get_ent (args[1]);
  e3 = bltin_get_ent (args[2]);
  e4 = bltin_get_ent (args[3]);

  x1 = (float) class_double (e1);
  x2 = (float) class_double (e2);
  y1 = (float) class_double (e3);
  y2 = (float) class_double (e4);

  cpgswin (x1, x2, y1, y2);

  ent_Clean (e1);
  ent_Clean (e2);
  ent_Clean (e3);
  ent_Clean (e4);

  return ent_Create_Rlab_Success();
}

/* **************************************************************
 * set window and adjust viewport to same aspect ratio
 *    SUBROUTINE PGWNAD (X1, X2, Y1, Y2)
 *       REAL X1, X2, Y1, Y2
 * ************************************************************** */

Ent *
_pg_cpgwnad (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *e3=0, *e4=0;
  float x1, x2, y1, y2;

  if (nargs != 4)
  {
    rerror ("pgwnad: 4 arguments required");
  }

  e1 = bltin_get_ent (args[0]);
  e2 = bltin_get_ent (args[1]);
  e3 = bltin_get_ent (args[2]);
  e4 = bltin_get_ent (args[3]);

  x1 = (float) class_double (e1);
  x2 = (float) class_double (e2);
  y1 = (float) class_double (e3);
  y2 = (float) class_double (e4);

  cpgwnad (x1, x2, y1, y2);

  ent_Clean (e1);
  ent_Clean (e2);
  ent_Clean (e3);
  ent_Clean (e4);

  return ent_Create_Rlab_Success();
}

/* **************************************************************
 * Select an open graphics device
 *    SUBROUTINE PGSLCT(ID)
 *       INTEGER ID
 * ************************************************************** */

Ent *
_pg_cpgslct (int nargs, Datum args[])
{
  Ent *e1=0;
  int istat;

  if (nargs != 1)
  {
    rerror ("pgslct: 1 argument allowed");
  }

  e1 = bltin_get_ent (args[0]);
  istat = (int) class_double (e1);

  cpgslct (istat);

  ent_Clean (e1);

  return ent_Create_Rlab_Success();
}

/* **************************************************************
 * Subdivide view surface into panels
 *     SUBROUTINE PGSUBP (NXSUB, NYSUB)
 *         INTEGER NXSUB, NYSUB
 * ************************************************************** */

Ent *
_pg_cpgsubp (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0;
  int nxsub, nysub;

  if (nargs != 2)
  {
    rerror ("pgsubp: 2 arguments allowed");
  }

  e1 = bltin_get_ent (args[0]);
  e2 = bltin_get_ent (args[1]);

  nxsub = (int) class_double (e1);
  nysub = (int) class_double (e2);

  cpgsubp (nxsub, nysub);

  ent_Clean (e1);
  ent_Clean (e2);

  return ent_Create_Rlab_Success();
}

/* **************************************************************
 * switch to a different panel on the view surface
 *    SUBROUTINE PGPANL(IX, IY)
 *        INTEGER IX, IY
 * ************************************************************** */

Ent *
_pg_cpgpanl (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0;
  int ix, iy;

  if (nargs != 2)
  {
    rerror ("pgpanl: 2 arguments allowed");
  }

  e1 = bltin_get_ent (args[0]);
  e2 = bltin_get_ent (args[1]);

  ix = (int) class_double (e1);
  iy = (int) class_double (e2);

  cpgpanl (ix, iy);

  ent_Clean (e1);
  ent_Clean (e2);

  return ent_Create_Rlab_Success();
}

/* **************************************************************
 * change the size of the view surface
 *     SUBROUTINE PGPAP (WIDTH, ASPECT)
 *       REAL WIDTH, ASPECT
 * ************************************************************** */

Ent *
_pg_cpgpap (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0;
  float width, aspect;

  if (nargs != 2)
  {
    rerror ("pgpap: 2 arguments allowed");
  }

  e1 = bltin_get_ent (args[0]);
  e2 = bltin_get_ent (args[1]);

  width = (float) class_double (e1);
  aspect = (float) class_double (e2);

  cpgpap (width, aspect);

  ent_Clean (e1);
  ent_Clean (e2);

  return ent_Create_Rlab_Success();
}

/* **************************************************************
 * Set window and viewport and draw labeled frame
 *   SUBROUTINE PGENV (XMIN, XMAX, YMIN, YMAX, JUST, AXIS)
 *       REAL XMIN, XMAX, YMIN, YMAX
 *       INTEGER JUST, AXIS
 * ************************************************************** */

Ent *
_pg_cpgenv (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *e3=0, *e4=0, *e5=0, *e6;
  float xmin, xmax, ymin, ymax;
  int just, axis;

  if (nargs != 6)
  {
    rerror ("pgenv: 6 arguments are required");
  }

  e1 = bltin_get_ent (args[0]);
  e2 = bltin_get_ent (args[1]);
  e3 = bltin_get_ent (args[2]);
  e4 = bltin_get_ent (args[3]);
  e5 = bltin_get_ent (args[4]);
  e6 = bltin_get_ent (args[5]);

  xmin = (float) class_double (e1);
  xmax = (float) class_double (e2);
  ymin = (float) class_double (e3);
  ymax = (float) class_double (e4);
  just = (int) class_double (e5);
  axis = (int) class_double (e6);

  cpgenv (xmin, xmax, ymin, ymax, just, axis);

  ent_Clean (e1);
  ent_Clean (e2);
  ent_Clean (e3);
  ent_Clean (e4);
  ent_Clean (e5);
  ent_Clean (e6);

  return ent_Create_Rlab_Success();
}

/* **************************************************************
 * Draw a polyline (curve defined by line-segments)
 *     SUBROUTINE PGLINE (N, XPTS, YPTS)
 *        INTEGER  N
 *        REAL     XPTS(*), YPTS(*)
 * ************************************************************** */

Ent *
_pg_cpgline (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *e3=0;
  float *xpts, *ypts;
  int n;

  if (nargs != 3)
  {
    rerror ("pgline: 3 arguments are required");
  }

  e1 = bltin_get_ent (args[0]);
  e2 = bltin_get_ent (args[1]);
  e3 = bltin_get_ent (args[2]);

  n = (int) class_double (e1);
  xpts = mdr_coerce_floatptr (class_matrix_real (e2));
  ypts = mdr_coerce_floatptr (class_matrix_real (e3));

  cpgline (n, xpts, ypts);

  ent_Clean (e1);
  ent_Clean (e2);
  ent_Clean (e3);

  GC_FREE (xpts);
  GC_FREE (ypts);

  return ent_Create_Rlab_Success();
}

/* **************************************************************
 * draw several graph markers
 *    SUBROUTINE PGPT (N, XPTS, YPTS, SYMBOL)
 *       INTEGER N
 *       REAL XPTS(*), YPTS(*)
 *       INTEGER SYMBOL
 * ************************************************************** */

Ent *
_pg_cpgpt (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *e3=0, *e4=0;
  float *xpts, *ypts;
  int n, symbol;

  if (nargs != 4)
  {
    rerror ("pgpt: 4 arguments are required");
  }

  e1 = bltin_get_ent (args[0]);
  e2 = bltin_get_ent (args[1]);
  e3 = bltin_get_ent (args[2]);
  e4 = bltin_get_ent (args[3]);

  n = (int) class_double (e1);
  xpts = mdr_coerce_floatptr (class_matrix_real (e2));
  ypts = mdr_coerce_floatptr (class_matrix_real (e3));
  symbol = (int) class_double (e4);

  cpgpt (n, xpts, ypts, symbol);

  ent_Clean (e1);
  ent_Clean (e2);
  ent_Clean (e3);
  ent_Clean (e4);

  GC_FREE (xpts);
  GC_FREE (ypts);

  return ent_Create_Rlab_Success();
}

/* **************************************************************
 * color image from a 2D data array
 *      SUBROUTINE PGIMAG (A, IDIM, JDIM, I1, I2, J1, J2,
 *     1                   A1, A2, TR)
 *          INTEGER IDIM, JDIM, I1, I2, J1, J2
 *          REAL    A(IDIM,JDIM), A1, A2, TR(6)
 * ************************************************************** */

Ent *
_pg_cpgimag (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *e3=0, *e4=0, *e5=0, *e6, *e7, *e8, *e9, *e10;
  int idim, jdim, i1, i2, j1, j2;
  float *a, a1, a2, *tr;

  if (nargs != 10)
  {
    rerror ("pgimag: 10 arguments are required");
  }

  e1 = bltin_get_ent (args[0]);
  e2 = bltin_get_ent (args[1]);
  e3 = bltin_get_ent (args[2]);
  e4 = bltin_get_ent (args[3]);
  e5 = bltin_get_ent (args[4]);
  e6 = bltin_get_ent (args[5]);
  e7 = bltin_get_ent (args[6]);
  e8 = bltin_get_ent (args[7]);
  e9 = bltin_get_ent (args[8]);
  e10 = bltin_get_ent (args[9]);

  a = mdr_coerce_floatptr (class_matrix_real (e1));
  idim = (int) class_double (e2);
  jdim = (int) class_double (e3);
  i1 = (int) class_double (e4);
  i2 = (int) class_double (e5);
  j1 = (int) class_double (e6);
  j2 = (int) class_double (e7);
  a1 = (float) class_double (e8);
  a2 = (float) class_double (e9);
  tr = mdr_coerce_floatptr (class_matrix_real (e10));

  cpgimag (a, idim, jdim, i1, i2, j1, j2, a1, a2, tr);

  GC_FREE (a);
  GC_FREE (tr);

  ent_Clean (e1);
  ent_Clean (e2);
  ent_Clean (e3);
  ent_Clean (e4);
  ent_Clean (e5);
  ent_Clean (e6);
  ent_Clean (e7);
  ent_Clean (e8);
  ent_Clean (e9);
  ent_Clean (e10);

  return ent_Create_Rlab_Success();
}

/* **************************************************************
 * gray-scale map of a 2D data array
 *   SUBROUTINE PGGRAY (A, IDIM, JDIM, I1, I2, J1, J2,
 *   1                   FG, BG, TR)
 *    INTEGER IDIM, JDIM, I1, I2, J1, J2
 *    REAL    A(IDIM,JDIM), FG, BG, TR(6)
 * ************************************************************** */

Ent *
_pg_cpggray (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *e3=0, *e4=0, *e5=0, *e6, *e7, *e8, *e9, *e10;
  int idim, jdim, i1, i2, j1, j2;
  float *a, fg, bg, *tr;

  if (nargs != 10)
  {
    rerror ("pggray: 10 arguments are required");
  }

  e1 = bltin_get_ent (args[0]);
  e2 = bltin_get_ent (args[1]);
  e3 = bltin_get_ent (args[2]);
  e4 = bltin_get_ent (args[3]);
  e5 = bltin_get_ent (args[4]);
  e6 = bltin_get_ent (args[5]);
  e7 = bltin_get_ent (args[6]);
  e8 = bltin_get_ent (args[7]);
  e9 = bltin_get_ent (args[8]);
  e10 = bltin_get_ent (args[9]);

  a = mdr_coerce_floatptr (class_matrix_real (e1));
  idim = (int) class_double (e2);
  jdim = (int) class_double (e3);
  i1 = (int) class_double (e4);
  i2 = (int) class_double (e5);
  j1 = (int) class_double (e6);
  j2 = (int) class_double (e7);
  fg = (float) class_double (e8);
  bg = (float) class_double (e9);
  tr = mdr_coerce_floatptr (class_matrix_real (e10));

  cpggray (a, idim, jdim, i1, i2, j1, j2, fg, bg, tr);

  GC_FREE (a);
  GC_FREE (tr);

  ent_Clean (e1);
  ent_Clean (e2);
  ent_Clean (e3);
  ent_Clean (e4);
  ent_Clean (e5);
  ent_Clean (e6);
  ent_Clean (e7);
  ent_Clean (e8);
  ent_Clean (e9);
  ent_Clean (e10);

  return ent_Create_Rlab_Success();
}

/* **************************************************************
 * annotate an image plot with a wedge
 *  SUBROUTINE PGWEDG(SIDE, DISP, WIDTH, FG, BG, LABEL)
 *     CHARACTER *(*) SIDE,LABEL
 *    REAL DISP, WIDTH, FG, BG
 * ************************************************************** */

Ent *
_pg_cpgwedg (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *e3=0, *e4=0, *e5=0, *e6;
  char *side, *label;
  float disp, width, fg, bg;

  if (nargs != 6)
  {
    rerror ("pgwedg: 6 arguments are required");
  }

  e1 = bltin_get_ent (args[0]);
  e2 = bltin_get_ent (args[1]);
  e3 = bltin_get_ent (args[2]);
  e4 = bltin_get_ent (args[3]);
  e5 = bltin_get_ent (args[4]);
  e6 = bltin_get_ent (args[5]);

  side = class_char_pointer (e1);
  disp = (float) class_double (e2);
  width = (float) class_double (e3);
  fg = (float) class_double (e4);
  bg = (float) class_double (e5);
  label = class_char_pointer (e6);

  cpgwedg (side, disp, width, fg, bg, label);

  ent_Clean (e1);
  ent_Clean (e2);
  ent_Clean (e3);
  ent_Clean (e4);
  ent_Clean (e5);
  ent_Clean (e6);

  return ent_Create_Rlab_Success();
}

/* **************************************************************
 * contour map of a 2D data array, with blanking
 *   SUBROUTINE PGCONB (A, IDIM, JDIM, I1, I2, J1, J2, C, NC, TR,
 *   1                   BLANK)
 *    INTEGER IDIM, JDIM, I1, I2, J1, J2, NC
 *    REAL    A(IDIM,JDIM), C(*), TR(6), BLANK
 * ************************************************************** */

Ent *
_pg_cpgconb (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *e3=0, *e4=0, *e5=0, *e6=0, *e7=0, *e8=0, *e9=0, *e10=0, *e11=0;
  int idim, jdim, i1, i2, j1, j2, nc;
  float *a, *c, *tr, blank;

  if (nargs != 11)
  {
    rerror ("pgconb: 11 arguments are required");
  }

  e1 = bltin_get_ent (args[0]);
  e2 = bltin_get_ent (args[1]);
  e3 = bltin_get_ent (args[2]);
  e4 = bltin_get_ent (args[3]);
  e5 = bltin_get_ent (args[4]);
  e6 = bltin_get_ent (args[5]);
  e7 = bltin_get_ent (args[6]);
  e8 = bltin_get_ent (args[7]);
  e9 = bltin_get_ent (args[8]);
  e10 = bltin_get_ent (args[9]);
  e11 = bltin_get_ent (args[10]);

  a = mdr_coerce_floatptr (class_matrix_real (e1));
  idim = (int) class_double (e2);
  jdim = (int) class_double (e3);
  i1 = (int) class_double (e4);
  i2 = (int) class_double (e5);
  j1 = (int) class_double (e6);
  j2 = (int) class_double (e7);
  c = mdr_coerce_floatptr (class_matrix_real (e8));
  nc = (int) class_double (e9);
  tr = mdr_coerce_floatptr (class_matrix_real (e10));
  blank = (float) class_double (e11);

  cpgconb (a, idim, jdim, i1, i2, j1, j2, c, nc, tr, blank);

  GC_FREE (a);
  GC_FREE (c);
  GC_FREE (tr);

  ent_Clean (e1);
  ent_Clean (e2);
  ent_Clean (e3);
  ent_Clean (e4);
  ent_Clean (e5);
  ent_Clean (e6);
  ent_Clean (e7);
  ent_Clean (e8);
  ent_Clean (e9);
  ent_Clean (e10);
  ent_Clean (e11);

  return ent_Create_Rlab_Success();
}

/* **************************************************************
 * label contour map of a 2D data array
 *    SUBROUTINE PGCONL (A, IDIM, JDIM, I1, I2, J1, J2, C, TR,
 *    1                   LABEL, INTVAL, MININT)
 *     INTEGER IDIM, JDIM, I1, J1, I2, J2, INTVAL, MININT
 *     REAL A(IDIM,JDIM), C, TR(6)
 *     CHARACTER*(*) LABEL
 * ************************************************************** */

Ent *
_pg_cpgconl (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *e3=0, *e4=0, *e5=0, *e6=0, *e7=0, *e8=0, *e9=0, *e10=0, *e11=0, *e12=0;
  int idim, jdim, i1, i2, j1, j2, intval, minint;
  float *a, c, *tr;
  char *label;

  if (nargs != 12)
  {
    rerror ("pgconl: 12 arguments are required");
  }

  e1 = bltin_get_ent (args[0]);
  e2 = bltin_get_ent (args[1]);
  e3 = bltin_get_ent (args[2]);
  e4 = bltin_get_ent (args[3]);
  e5 = bltin_get_ent (args[4]);
  e6 = bltin_get_ent (args[5]);
  e7 = bltin_get_ent (args[6]);
  e8 = bltin_get_ent (args[7]);
  e9 = bltin_get_ent (args[8]);
  e10 = bltin_get_ent (args[9]);
  e11 = bltin_get_ent (args[10]);
  e12 = bltin_get_ent (args[11]);

  a = mdr_coerce_floatptr (class_matrix_real (e1));
  idim = (int) class_double (e2);
  jdim = (int) class_double (e3);
  i1 = (int) class_double (e4);
  i2 = (int) class_double (e5);
  j1 = (int) class_double (e6);
  j2 = (int) class_double (e7);
  c = (float) class_double (e8);
  tr = mdr_coerce_floatptr (class_matrix_real (e9));
  label = class_char_pointer (e10);
  intval = (int) class_double (e11);
  minint = (int) class_double (e12);

  cpgconl (a, idim, jdim, i1, i2, j1, j2, c, tr, label, intval, minint);

  GC_FREE (a);
  GC_FREE (tr);

  ent_Clean (e1);
  ent_Clean (e2);
  ent_Clean (e3);
  ent_Clean (e4);
  ent_Clean (e5);
  ent_Clean (e6);
  ent_Clean (e7);
  ent_Clean (e8);
  ent_Clean (e9);
  ent_Clean (e10);
  ent_Clean (e11);
  ent_Clean (e12);

  return ent_Create_Rlab_Success();
}

/* **************************************************************
 * contour map of a 2D data array (contour-following)
 *   SUBROUTINE PGCONT (A, IDIM, JDIM, I1, I2, J1, J2, C, NC, TR)
 *     INTEGER IDIM, JDIM, I1, J1, I2, J2, NC
 *     REAL A(IDIM,JDIM), C(*), TR(6)
 * ************************************************************** */

Ent *
_pg_cpgcont (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *e3=0, *e4=0, *e5=0, *e6=0, *e7=0, *e8=0, *e9=0, *e10=0;
  int idim, jdim, i1, i2, j1, j2, nc;
  float *a, *c, *tr;

  if (nargs != 10)
  {
    rerror ("pgcont: 10 arguments are required");
  }

  e1 = bltin_get_ent (args[0]);
  e2 = bltin_get_ent (args[1]);
  e3 = bltin_get_ent (args[2]);
  e4 = bltin_get_ent (args[3]);
  e5 = bltin_get_ent (args[4]);
  e6 = bltin_get_ent (args[5]);
  e7 = bltin_get_ent (args[6]);
  e8 = bltin_get_ent (args[7]);
  e9 = bltin_get_ent (args[8]);
  e10 = bltin_get_ent (args[9]);

  a = mdr_coerce_floatptr (class_matrix_real (e1));
  idim = (int) class_double (e2);
  jdim = (int) class_double (e3);
  i1 = (int) class_double (e4);
  i2 = (int) class_double (e5);
  j1 = (int) class_double (e6);
  j2 = (int) class_double (e7);
  c = mdr_coerce_floatptr (class_matrix_real (e8));
  nc = (int) class_double (e9);
  tr = mdr_coerce_floatptr (class_matrix_real (e10));

  cpgcont (a, idim, jdim, i1, i2, j1, j2, c, nc, tr);

  GC_FREE (a);
  GC_FREE (c);
  GC_FREE (tr);

  ent_Clean (e1);
  ent_Clean (e2);
  ent_Clean (e3);
  ent_Clean (e4);
  ent_Clean (e5);
  ent_Clean (e6);
  ent_Clean (e7);
  ent_Clean (e8);
  ent_Clean (e9);
  ent_Clean (e10);

  return ent_Create_Rlab_Success();
}

/* **************************************************************
 * fill between two contours
 *   SUBROUTINE PGCONF (A, IDIM, JDIM, I1, I2, J1, J2, C1, C2, TR)
 *     INTEGER IDIM, JDIM, I1, I2, J1, J2
 *     REAL    A(IDIM,JDIM), C1, C2, TR(6)
 * ************************************************************** */

Ent *
_pg_cpgconf (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *e3=0, *e4=0, *e5=0, *e6, *e7, *e8, *e9, *e10;
  int idim, jdim, i1, i2, j1, j2;
  float *a, c1, c2, *tr;

  if (nargs != 10)
  {
    rerror ("pgconf: 10 arguments are required");
  }

  e1 = bltin_get_ent (args[0]);
  e2 = bltin_get_ent (args[1]);
  e3 = bltin_get_ent (args[2]);
  e4 = bltin_get_ent (args[3]);
  e5 = bltin_get_ent (args[4]);
  e6 = bltin_get_ent (args[5]);
  e7 = bltin_get_ent (args[6]);
  e8 = bltin_get_ent (args[7]);
  e9 = bltin_get_ent (args[8]);
  e10 = bltin_get_ent (args[9]);

  a = mdr_coerce_floatptr (class_matrix_real (e1));
  idim = (int) class_double (e2);
  jdim = (int) class_double (e3);
  i1 = (int) class_double (e4);
  i2 = (int) class_double (e5);
  j1 = (int) class_double (e6);
  j2 = (int) class_double (e7);
  c1 = (float) class_double (e8);
  c2 = (float) class_double (e8);
  tr = mdr_coerce_floatptr (class_matrix_real (e10));

  cpgconf (a, idim, jdim, i1, i2, j1, j2, c1, c2, tr);

  GC_FREE (a);
  GC_FREE (tr);

  ent_Clean (e1);
  ent_Clean (e2);
  ent_Clean (e3);
  ent_Clean (e4);
  ent_Clean (e5);
  ent_Clean (e6);
  ent_Clean (e7);
  ent_Clean (e8);
  ent_Clean (e9);
  ent_Clean (e10);

  return ent_Create_Rlab_Success();
}

/* **************************************************************
 * draw labeled frame around viewport
 *    SUBROUTINE PGBOX (XOPT, XTICK, NXSUB, YOPT, YTICK, NYSUB)
 *        CHARACTER*(*) XOPT, YOPT
 *        REAL XTICK, YTICK
 *        INTEGER NXSUB, NYSUB
 * ************************************************************** */

Ent *
_pg_cpgbox (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *e3=0, *e4=0, *e5=0, *e6;
  char *xopt, *yopt;
  float xtick, ytick;
  int nxsub, nysub;

  if (nargs != 6)
  {
    rerror ("pgbox: 6 arguments are required");
  }

  e1 = bltin_get_ent (args[0]);
  e2 = bltin_get_ent (args[1]);
  e3 = bltin_get_ent (args[2]);
  e4 = bltin_get_ent (args[3]);
  e5 = bltin_get_ent (args[4]);
  e6 = bltin_get_ent (args[5]);

  xopt = class_char_pointer (e1);
  xtick = (float) class_double (e2);
  nxsub = (int) class_double (e3);

  yopt = class_char_pointer (e4);
  ytick = (float) class_double (e5);
  nysub = (int) class_double (e6);

  cpgbox (xopt, xtick, nxsub, yopt, ytick, nysub);

  ent_Clean (e1);
  ent_Clean (e2);
  ent_Clean (e3);
  ent_Clean (e4);
  ent_Clean (e5);
  ent_Clean (e6);

  return ent_Create_Rlab_Success();
}

/* **************************************************************
 * draw frame and write (DD) HH MM SS.S labelling
 *     SUBROUTINE PGTBOX (XOPT, XTICK, NXSUB, YOPT, YTICK, NYSUB)
 *       REAL XTICK, YTICK
 *       INTEGER NXSUB, NYSUB
 *       CHARACTER XOPT*(*), YOPT*(*)
 * ************************************************************** */

Ent *
_pg_cpgtbox (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *e3=0, *e4=0, *e5=0, *e6;
  char *xopt, *yopt;
  float xtick, ytick;
  int nxsub, nysub;

  if (nargs != 6)
  {
    rerror ("pgtbox: 6 arguments are required");
  }

  e1 = bltin_get_ent (args[0]);
  e2 = bltin_get_ent (args[1]);
  e3 = bltin_get_ent (args[2]);
  e4 = bltin_get_ent (args[3]);
  e5 = bltin_get_ent (args[4]);
  e6 = bltin_get_ent (args[5]);

  xopt = class_char_pointer (e1);
  xtick = (float) class_double (e2);
  nxsub = (int) class_double (e3);

  yopt = class_char_pointer (e4);
  ytick = (float) class_double (e5);
  nysub = (int) class_double (e6);

  cpgtbox (xopt, xtick, nxsub, yopt, ytick, nysub);

  ent_Clean (e1);
  ent_Clean (e2);
  ent_Clean (e3);
  ent_Clean (e4);
  ent_Clean (e5);
  ent_Clean (e6);

  return ent_Create_Rlab_Success();
}

/* **************************************************************
 * set color representation
 *    SUBROUTINE PGSCR (CI, CR, CG, CB)
 *        INTEGER CI
 *        REAL    CR, CG, CB
 * ************************************************************** */

Ent *
_pg_cpgscr (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *e3=0, *e4=0;
  int ci;
  float cr, cg, cb;

  if (nargs != 4)
  {
    rerror ("pgscr: 4 arguments are required");
  }

  e1 = bltin_get_ent (args[0]);
  e2 = bltin_get_ent (args[1]);
  e3 = bltin_get_ent (args[2]);
  e4 = bltin_get_ent (args[3]);

  ci = (int) class_double (e1);
  cr = (float) class_double (e2);
  cg = (float) class_double (e3);
  cb = (float) class_double (e4);

  cpgscr (ci, cr, cg, cb);

  ent_Clean (e1);
  ent_Clean (e2);
  ent_Clean (e3);
  ent_Clean (e4);

  return ent_Create_Rlab_Success();
}

/* **************************************************************
 * set text background color index
 *   SUBROUTINE PGSTBG (TBCI)
 *       INTEGER  TBCI
 * ************************************************************** */

Ent *
_pg_cpgstbg (int nargs, Datum args[])
{
  Ent *e1=0;
  int tcbi;

  if (nargs != 1)
  {
    rerror ("pgstbg: 1 argument allowed");
  }

  e1 = bltin_get_ent (args[0]);
  tcbi = (int) class_double (e1);

  cpgstbg (tcbi);

  ent_Clean (e1);

  return ent_Create_Rlab_Success();
}

/* **************************************************************
 * set line style
 *    SUBROUTINE PGSLS (LS)
 *         INTEGER  LS
 * ************************************************************** */

Ent *
_pg_cpgsls (int nargs, Datum args[])
{
  Ent *e1=0;
  int ls;

  if (nargs != 1)
  {
    rerror ("pgsls: 1 argument allowed");
  }

  e1 = bltin_get_ent (args[0]);
  ls = (int) class_double (e1);

  cpgsls (ls);

  ent_Clean (e1);

  return ent_Create_Rlab_Success();
}

/* **************************************************************
 * set line width
 *    SUBROUTINE PGSLW (LW)
 *         INTEGER  LW
 * ************************************************************** */

Ent *
_pg_cpgslw (int nargs, Datum args[])
{
  Ent *e1=0;
  int lw;

  if (nargs != 1)
  {
    rerror ("pgslw: 1 argument allowed");
  }

  e1 = bltin_get_ent (args[0]);
  lw = (int) class_double (e1);

  cpgslw (lw);

  ent_Clean (e1);

  return ent_Create_Rlab_Success();
}

/* **************************************************************
 * close all open graphics devices
 *    SUBROUTINE PGEND
 * ************************************************************** */

Ent *
_pg_cpgend (int nargs, Datum args[])
{
  if (nargs != 0)
  {
    rerror ("pgend: 0 arguments allowed");
  }

  cpgend ();

  return ent_Create_Rlab_Success();
}

/* **************************************************************
 * close the selected graphics device
 *    SUBROUTINE PGCLOS
 * ************************************************************** */

Ent *
_pg_cpgclos (int nargs, Datum args[])
{
  if (nargs != 0)
  {
    rerror ("pgclos: 0 arguments allowed");
  }

  cpgclos ();

  return ent_Create_Rlab_Success();
}

/* **************************************************************
 * Advance to a new page.
 *    SUBROUTINE PGPAGE
 * ************************************************************** */

Ent *
_pg_cpgpage (int nargs, Datum args[])
{
  if (nargs != 0)
  {
    rerror ("pgpage: 0 arguments allowed");
  }

  cpgpage ();

  return ent_Create_Rlab_Success();
}

/* **************************************************************
 * erase all graphics from current page
 *    SUBROUTINE PGERAS
 * ************************************************************** */

Ent *
_pg_cpgeras (int nargs, Datum args[])
{
  if (nargs != 0)
  {
    rerror ("pgeras: 0 arguments allowed");
  }

  cpgeras ();

  return ent_Create_Rlab_Success();
}

/* **************************************************************
 * erase text from graphics display
 *    SUBROUTINE PGETXT
 * ************************************************************** */

Ent *
_pg_cpgetxt (int nargs, Datum args[])
{
  if (nargs != 0)
  {
    rerror ("pgetxt: 0 arguments allowed");
  }

  cpgetxt ();

  return ent_Create_Rlab_Success();
}

/* **************************************************************
 * write text (horizontal, left-justified)
 *    SUBROUTINE PGTEXT (X, Y, TEXT)
 *       REAL X, Y
 *       CHARACTER*(*) TEXT
 * ************************************************************** */

Ent *
_pg_cpgtext (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *e3=0;
  float x, y;
  char *text;

  if (nargs != 3)
  {
    rerror ("pgtext: 3 arguments are required");
  }

  e1 = bltin_get_ent (args[0]);
  e2 = bltin_get_ent (args[1]);
  e3 = bltin_get_ent (args[2]);

  x = (float) class_double (e1);
  y = (float) class_double (e2);
  text = class_char_pointer (e3);

  cpgtext (x, y, text);

  ent_Clean (e1);
  ent_Clean (e2);
  ent_Clean (e3);

  return ent_Create_Rlab_Success();
}

/* **************************************************************
 * write text at arbitrary position and angle
 *    SUBROUTINE PGPTXT (X, Y, ANGLE, FJUST, TEXT)
 *      REAL X, Y, ANGLE, FJUST
 *      CHARACTER*(*) TEXT
 * ************************************************************** */

Ent *
_pg_cpgptxt (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *e3=0, *e4=0, *e5=0;
  float x, y, angle, fjust;
  char *text;

  if (nargs != 5)
  {
    rerror ("pgptxt: 5 arguments are required");
  }

  e1 = bltin_get_ent (args[0]);
  e2 = bltin_get_ent (args[1]);
  e3 = bltin_get_ent (args[2]);
  e4 = bltin_get_ent (args[3]);
  e5 = bltin_get_ent (args[4]);

  x = (float) class_double (e1);
  y = (float) class_double (e2);
  angle = (float) class_double (e3);
  fjust = (float) class_double (e4);
  text = class_char_pointer (e5);

  cpgptxt (x, y, angle, fjust, text);

  ent_Clean (e1);
  ent_Clean (e2);
  ent_Clean (e3);
  ent_Clean (e4);
  ent_Clean (e5);

  return ent_Create_Rlab_Success();
}

/* **************************************************************
 *  write text at position relative to viewport
 *  SUBROUTINE PGMTXT (SIDE, DISP, COORD, FJUST, TEXT)
 *    CHARACTER*(*) SIDE, TEXT
 *    REAL DISP, COORD, FJUST
 * Arguments:
   SIDE   (input)  : must include one of the characters 'B', 'L', 'T',
                     or 'R' signifying the Bottom, Left, Top, or Right
                     margin of the viewport. If it includes 'LV' or
                     'RV', the string is written perpendicular to the
                     frame rather than parallel to it.
   DISP   (input)  : the displacement of the character string from the
                     specified edge of the viewport, measured outwards
                     from the viewport in units of the character
                     height. Use a negative value to write inside the
                     viewport, a positive value to write outside.
   COORD  (input)  : the location of the character string along the
                     specified edge of the viewport, as a fraction of
                     the length of the edge.
   FJUST  (input)  : controls justification of the string parallel to
                     the specified edge of the viewport. If
                     FJUST = 0.0, the left-hand end of the string will
                     be placed at COORD; if JUST = 0.5, the center of
                     the string will be placed at COORD; if JUST = 1.0,
                     the right-hand end of the string will be placed at
                     at COORD. Other values between 0 and 1 give inter-
                     mediate placing, but they are not very useful.
   TEXT   (input) :  the text string to be plotted. Trailing spaces are
                     ignored when justifying the string, but leading
                     spaces are significant.
 * ************************************************************** */

Ent *
_pg_cpgmtxt (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *e3=0, *e4=0, *e5=0;
  char *side, *text;
  float disp, coord, fjust;

  if (nargs != 5)
  {
    rerror ("pgmtxt: 5 arguments are required");
  }

  e1 = bltin_get_ent (args[0]);
  e2 = bltin_get_ent (args[1]);
  e3 = bltin_get_ent (args[2]);
  e4 = bltin_get_ent (args[3]);
  e5 = bltin_get_ent (args[4]);

  side = class_char_pointer (e1);
  disp = (float) class_double (e2);
  coord = (float) class_double (e3);
  fjust = (float) class_double (e4);
  text = class_char_pointer (e5);

  cpgmtxt (side, disp, coord, fjust, text);

  ent_Clean (e1);
  ent_Clean (e2);
  ent_Clean (e3);
  ent_Clean (e4);
  ent_Clean (e5);

  return ent_Create_Rlab_Success();
}

/* **************************************************************
 * set character font
 *    SUBROUTINE PGSCF (FONT)
 *      INTEGER  FONT
 * ************************************************************** */

Ent *
_pg_cpgscf (int nargs, Datum args[])
{
  Ent *e1=0;
  int font;

  if (nargs != 1)
  {
    rerror ("pgscf: 1 argument allowed");
  }

  e1 = bltin_get_ent (args[0]);

  font = (int) class_double (e1);

  cpgscf (font);

  ent_Clean (e1);

  return ent_Create_Rlab_Success();
}

/* **************************************************************
 * set character height
 *    SUBROUTINE PGSCH (SIZE)
 *       REAL SIZE
 * ************************************************************** */

Ent *
_pg_cpgsch (int nargs, Datum args[])
{
  Ent *e1=0;
  float size;

  if (nargs != 1)
  {
    rerror ("pgsch: 1 argument allowed");
  }

  e1 = bltin_get_ent (args[0]);

  size = (float) class_double (e1);

  cpgsch (size);

  ent_Clean (e1);

  return ent_Create_Rlab_Success();
}

/* **************************************************************
 * control new page prompting
 *      SUBROUTINE PGASK (FLAG)
 *         LOGICAL FLAG
 *
 * Set to false for NO prompting.
 * ************************************************************** */

Ent *
_pg_cpgask (int nargs, Datum args[])
{
  Ent *e1=0;
  int flag;

  if (nargs != 1)
  {
    rerror ("pgask: 1 argument allowed");
  }

  e1 = bltin_get_ent (args[0]);

  flag = (int) class_double (e1);

  cpgask (flag);

  ent_Clean (e1);

  return ent_Create_Rlab_Success();
}

/* **************************************************************
 * Write labels for x-axis, y-axis, and top of plot
 *     SUBROUTINE PGLAB (XLBL, YLBL, TOPLBL)
 *         CHARACTER*(*) XLBL, YLBL, TOPLBL
 * ************************************************************** */

Ent *
_pg_cpglab (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *e3=0;
  char *xlbl, *ylbl, *toplbl;

  if (nargs != 3)
  {
    rerror ("pglab: 3 arguments are required");
  }

  e1 = bltin_get_ent (args[0]);
  e2 = bltin_get_ent (args[1]);
  e3 = bltin_get_ent (args[2]);

  xlbl = class_char_pointer (e1);
  ylbl = class_char_pointer (e2);
  toplbl = class_char_pointer (e3);

  cpglab (xlbl, ylbl, toplbl);

  ent_Clean (e1);
  ent_Clean (e2);
  ent_Clean (e3);

  return ent_Create_Rlab_Success();
}

/* **************************************************************
 * update display
 *     SUBROUTINE PGUPDT
 * ************************************************************** */

Ent *
_pg_cpgupdt (int nargs, Datum args[])
{
  if (nargs != 0)
  {
    rerror ("pgupdt: no arguments allowed");
  }

  cpgupdt ();

  return ent_Create_Rlab_Success();
}

/* **************************************************************
 * set color index
 *     SUBROUTINE PGSCI (CI)
 *        INTEGER  CI
 * ************************************************************** */

Ent *
_pg_cpgsci (int nargs, Datum args[])
{
  Ent *e1=0;
  int ci;

  if (nargs != 1)
  {
    rerror ("psci: 1 argument required");
  }

  e1 = bltin_get_ent (args[0]);

  ci = (int) class_double (e1);

  cpgsci (ci);

  ent_Clean (e1);

  return ent_Create_Rlab_Success();
}

/* **************************************************************
 * set color index range
 *     SUBROUTINE PGSCIR(ICILO, ICIHI)
 *        INTEGER   ICILO, ICIHI
 * ************************************************************** */

Ent *
_pg_cpgscir (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0;
  int icilo, icihi;

  if (nargs != 2)
  {
    rerror ("pgscir: 2 arguments required");
  }

  e1 = bltin_get_ent (args[0]);
  e2 = bltin_get_ent (args[1]);

  icilo = (int) class_double (e1);
  icihi = (int) class_double (e2);

  cpgscir (icilo, icihi);

  ent_Clean (e1);
  ent_Clean (e2);

  return ent_Create_Rlab_Success();
}

/* **************************************************************
 * inquire color index
 *     SUBROUTINE PGQCI(CI)
 *        INTEGER   CI
 * ************************************************************** */

Ent *
_pg_cpgqci (int nargs, Datum args[])
{
  int ci;
  MDR *color;

  if (nargs != 0)
  {
    rerror ("pqsci: 0 arguments allowed");
  }

  cpgqci (&ci);

  color = mdr_Create (1, 1);
  MdrV0 (color, 0) = (double) ci;

  return ent_Create_Rlab_Success();
}

/* **************************************************************
 * inquire color index range
 *     SUBROUTINE PGQCIR(ICILO, ICIHI)
 *        INTEGER   ICILO, ICIHI
 * ************************************************************** */

Ent *
_pg_cpgqcir (int nargs, Datum args[])
{
  int ci1, ci2;
  MDR *color;

  if (nargs != 0)
  {
    rerror ("pgqcir: 0 arguments allowed");
  }

  cpgqcir (&ci1, &ci2);

  color = mdr_Create (1, 2);
  MdrV0 (color, 0) = (double) ci1;
  MdrV0 (color, 1) = (double) ci2;

  return ent_Create_Rlab_Success();
}

/* **************************************************************
 * begin batch of output (buffer)
 *   SUBROUTINE PGBBUF
 * ************************************************************** */

Ent *
_pg_cpgbbuf (int nargs, Datum args[])
{
  if (nargs != 0)
  {
    rerror ("pgbbuf: 0 arguments allowed");
  }

  cpgbbuf ();

  return ent_Create_Rlab_Success();
}

/* **************************************************************
 * end batch of output (buffer)
 *   SUBROUTINE PGEBUF
 * ************************************************************** */

Ent *
_pg_cpgebuf (int nargs, Datum args[])
{
  if (nargs != 0)
  {
    rerror ("pgebuf: 0 arguments allowed");
  }

  cpgebuf ();

  return ent_Create_Rlab_Success();
}

/* **************************************************************
 * inquire color capability
 *   SUBROUTINE PGQCOL (CI1, CI2)
 *      INTEGER  CI1, CI2
 *
 * Returns a 1x2 matrix: [ CI1, CI2 ]
 * ************************************************************** */

Ent *
_pg_cpgqcol (int nargs, Datum args[])
{
  Ent *rent;
  int ci1, ci2;
  MDR *color;

  if (nargs != 0)
  {
    rerror ("pgqcol: 0 arguments allowed");
  }

  cpgqcol (&ci1, &ci2);

  color = mdr_Create (1, 2);
  MdrV0 (color, 0) = (double) ci1;
  MdrV0 (color, 1) = (double) ci2;

  rent = ent_Create ();
  ent_data (rent) = color;
  ent_type (rent) = MATRIX_DENSE_REAL;
  return (rent);
}

/* **************************************************************
 * move pen (change current pen position)
 *    SUBROUTINE PGMOVE (X, Y)
 *        REAL X, Y
 * ************************************************************** */

Ent *
_pg_cpgmove (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0;
  float x, y;

  if (nargs != 2)
  {
    rerror ("pgmove: 2 arguments allowed");
  }

  e1 = bltin_get_ent (args[0]);
  e2 = bltin_get_ent (args[1]);

  x = (float) class_double (e1);
  y = (float) class_double (e2);

  cpgmove (x, y);

  ent_Clean (e1);
  ent_Clean (e2);

  return ent_Create_Rlab_Success();
}

/* **************************************************************
 *  draw a line from the current pen position to a point
 *  SUBROUTINE PGDRAW (X, Y)
 *     REAL X, Y
 * ************************************************************** */

Ent *
_pg_cpgdraw (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0;
  float x, y;

  if (nargs != 2)
  {
    rerror ("pgdraw: 2 arguments allowed");
  }

  e1 = bltin_get_ent (args[0]);
  e2 = bltin_get_ent (args[1]);

  x = (float) class_double (e1);
  y = (float) class_double (e2);

  cpgdraw (x, y);

  ent_Clean (e1);
  ent_Clean (e2);

  return ent_Create_Rlab_Success();
}

/* **************************************************************
 *  set fill-area style
 *    SUBROUTINE PGSFS (FS)
 *        INTEGER  FS

  Argument:
 FS     (input)  : the fill-area style to be used for subsequent
                   plotting:
                     FS = 1 => solid (default)
                     FS = 2 => outline
                     FS = 3 => hatched
                     FS = 4 => cross-hatched
                   Other values give an error message and are
                   treated as 2.
 * ************************************************************** */

Ent *
_pg_cpgsfs (int nargs, Datum args[])
{
  Ent *e1=0;
  int fs;

  if (nargs != 1)
  {
    rerror ("pgsfs: 1 argument allowed");
  }

  e1 = bltin_get_ent (args[0]);
  fs = (int) class_double (e1);

  if (fs < 1 || fs > 4)
  {
    rerror ("pgsfs: argument must be either 1, 2, 3, or 4");
  }

  cpgsfs (fs);

  ent_Clean (e1);

  return ent_Create_Rlab_Success();
}

/* **************************************************************
 *  set hatching style
 *    SUBROUTINE PGSHS (ANGLE, SEPN, PHASE)
 *       REAL ANGLE, SEPN, PHASE
 Arguments:
 ANGLE  (input)  : the angle the hatch lines make with the
                   horizontal, in degrees, increasing
                   counterclockwise (this is an angle on the
                   view surface, not in world-coordinate space).
 SEPN   (input)  : the spacing of the hatch lines. The unit spacing
                   is 1 percent of the smaller of the height or
                   width of the view surface. This should not be
                   zero.
 PHASE  (input)  : a real number between 0 and 1; the hatch lines
                   are displaced by this fraction of SEPN from a
                   fixed reference.  Adjacent regions hatched with the
                   same PHASE have contiguous hatch lines. To hatch
                   a region with alternating lines of two colors,
                   fill the area twice, with PHASE=0.0 for one color
                   and PHASE=0.5 for the other color.
 * ************************************************************** */

Ent *
_pg_cpgshs (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *e3=0;
  float angle, sepn, phase;

  if (nargs != 3)
  {
    rerror ("pgshs: 3 arguments are required");
  }

  e1 = bltin_get_ent (args[0]);
  e2 = bltin_get_ent (args[1]);
  e3 = bltin_get_ent (args[2]);

  angle = (float) class_double (e1);
  sepn = (float) class_double (e2);
  phase = (float) class_double (e3);

  cpgshs (angle, sepn, phase);

  ent_Clean (e1);
  ent_Clean (e2);
  ent_Clean (e3);

  return ent_Create_Rlab_Success();
}

/* **************************************************************
 * draw a polygon, using fill-area attributes
 *   SUBROUTINE PGPOLY (N, XPTS, YPTS)
 *     INTEGER N
 *     REAL XPTS(*), YPTS(*)
 Arguments:
 N      (input)  : number of points defining the polygon; the
                   line consists of N straight-line segments,
                   joining points 1 to 2, 2 to 3,... N-1 to N, N to 1.
                   N should be greater than 2 (if it is 2 or less,
                   nothing will be drawn).
 XPTS   (input)  : world x-coordinates of the vertices.
 YPTS   (input)  : world y-coordinates of the vertices.
                   Note: the dimension of arrays XPTS and YPTS must be
                   greater than or equal to N.
 * ************************************************************** */

Ent *
_pg_cpgpoly (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *e3=0;
  int n;
  float *xpts, *ypts;

  if (nargs != 3)
  {
    rerror ("pgpoly: 3 arguments are required");
  }

  e1 = bltin_get_ent (args[0]);
  e2 = bltin_get_ent (args[1]);
  e3 = bltin_get_ent (args[2]);

  n = (int) class_double (e1);
  xpts = mdr_coerce_floatptr (class_matrix_real (e2));
  ypts = mdr_coerce_floatptr (class_matrix_real (e3));

  cpgpoly (n, xpts, ypts);

  ent_Clean (e1);
  ent_Clean (e2);
  ent_Clean (e3);

  GC_FREE (xpts);
  GC_FREE (ypts);

  return ent_Create_Rlab_Success();
}

/* **************************************************************
 * histogram of binned data
 *  SUBROUTINE PGBIN (NBIN, X, DATA, CENTER)
 *    INTEGER NBIN
 *    REAL X(*), DATA(*)
 *    LOGICAL CENTER
 Arguments:
 NBIN   (input)  : number of values.
 X      (input)  : abscissae of bins.
 DATA   (input)  : data values of bins.
 CENTER (input)  : if .TRUE., the X values denote the center of the
                   bin; if .FALSE., the X values denote the lower
                   edge (in X) of the bin.
 * ************************************************************** */

Ent *
_pg_cpgbin (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *e3=0, *e4=0;
  int nbin, center;
  float *x, *data;

  if (nargs != 4)
  {
    rerror ("pgbin: 4 arguments are required");
  }

  e1 = bltin_get_ent (args[0]);
  e2 = bltin_get_ent (args[1]);
  e3 = bltin_get_ent (args[2]);
  e4 = bltin_get_ent (args[3]);

  nbin = (int) class_double (e1);
  x = mdr_coerce_floatptr (class_matrix_real (e2));
  data = mdr_coerce_floatptr (class_matrix_real (e3));
  center = (int) class_double (e4);

  cpgbin (nbin, x, data, center);

  ent_Clean (e1);
  ent_Clean (e2);
  ent_Clean (e3);
  ent_Clean (e4);

  GC_FREE (x);
  GC_FREE (data);

  return ent_Create_Rlab_Success();
}

/* **************************************************************
 * histogram of unbinned data
 * SUBROUTINE PGHIST(N, DATA, DATMIN, DATMAX, NBIN, PGFLAG)
 *    INTEGER N
 *    REAL    DATA(*)
 *    REAL    DATMIN, DATMAX
 *    INTEGER NBIN, PGFLAG
 * ************************************************************** */

Ent *
_pg_cpghist (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *e3=0, *e4=0, *e5=0, *e6;
  int n, nbin, pgflag;
  float *data, datmin, datmax;

  if (nargs != 6)
  {
    rerror ("pghist: 6 arguments are required");
  }

  e1 = bltin_get_ent (args[0]);
  e2 = bltin_get_ent (args[1]);
  e3 = bltin_get_ent (args[2]);
  e4 = bltin_get_ent (args[3]);
  e5 = bltin_get_ent (args[4]);
  e6 = bltin_get_ent (args[5]);

  n = (int) class_double (e1);
  data = mdr_coerce_floatptr (class_matrix_real (e2));
  datmin = (float) class_double (e3);
  datmax = (float) class_double (e4);
  nbin = (int) class_double (e5);
  pgflag = (int) class_double (e6);

  cpghist (n, data, datmin, datmax, nbin, pgflag);

  ent_Clean (e1);
  ent_Clean (e2);
  ent_Clean (e3);
  ent_Clean (e4);
  ent_Clean (e5);
  ent_Clean (e6);

  GC_FREE (data);

  return ent_Create_Rlab_Success();
}

/* **************************************************************
 * Coerce a real-dense-matrix into an array of floats.
 * ************************************************************** */

float *
mdr_coerce_floatptr (MDR * m)
{
  float *ret = 0;
  int i = 0;
  int size = 0;

  /* Compute the size of the new float array. */
  size = MNR (m) * MNC (m);

  /* Check for the possibility of a zero (empty) input. */
  if (size == 0)
  {
    return (0);
  }

  /* Get memory for the new array. */
  ret = (float *) GC_MAIOP (size * sizeof (float));

  /* Convert the double precision array into float. */
  if (m->type == RLAB_TYPE_INT32)
  {
    for (i = 0; i < size; i++)
    {
      ret[i] = (float) MdiV0 (m, i);
    }
  }
  else
  {
    for (i = 0; i < size; i++)
    {
      ret[i] = (float) MdrV0 (m, i);
    }
  }

  return (ret);
}

/* **************************************************************
 *  choose axis limits
 *   SUBROUTINE PGRNGE (X1, X2, XLO, XHI)
 *      REAL X1, X2, XLO, XHI
 * Arguments:
   X1, X2 (input)  : the data range (X1<X2), ie, the min and max values
                     to be plotted.
   XLO, XHI (output) : suitable values to use as the extremes of a graph
                     axis (XLO <= X1, XHI >= X2).*
 * ************************************************************** */

Ent *
_pg_cpgrnge (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0;
  float x1, x2, xlo, xhi;
  MDR *mrange;

  if (nargs != 2)
  {
    rerror ("pgrnge: 2 arguments are required");
  }

  e1 = bltin_get_ent (args[0]);
  e2 = bltin_get_ent (args[1]);

  x1  = (float) class_double (e1);
  x2  = (float) class_double (e2);

  cpgrnge (x1, x2, &xlo, &xhi);

  mrange = mdr_Create (1, 2);
  MdrV0 (mrange, 0) = xlo;
  MdrV0 (mrange, 1) = xhi;

  ent_Clean (e1);
  ent_Clean (e2);

  return ent_Create_Rlab_Success();
}

/* **************************************************************
 * save PGPLOT attributes
 *   SUBROUTINE PGSAVE
 * ************************************************************** */

Ent *
_pg_cpgsave (int nargs, Datum args[])
{
  if (nargs != 0)
  {
    rerror ("pgsave: 0 arguments allowed");
  }

  cpgsave ();

  return ent_Create_Rlab_Success();
}

/* **************************************************************
 * restore PGPLOT attributes
 *   ENTRY PGUNSA
 * ************************************************************** */

Ent *
_pg_cpgunsa (int nargs, Datum args[])
{
  if (nargs != 0)
  {
    rerror ("pgunsa: 0 arguments allowed");
  }

  cpgunsa ();

  return ent_Create_Rlab_Success();
}

/* **************************************************************
 * draw a rectangle, using fill-area attributes
 *   SUBROUTINE PGRECT (X1, X2, Y1, Y2)
 *     REAL X1, X2, Y1, Y2
 * Arguments:
   X1, X2 (input) : the horizontal range of the rectangle.
   Y1, Y2 (input) : the vertical range of the rectangle.
 * ************************************************************** */

Ent *
_pg_cpgrect (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *e3=0, *e4=0;
  float x1, x2, y1, y2;

  if (nargs != 4)
  {
    rerror ("pgrect: 4 arguments required");
  }

  e1 = bltin_get_ent (args[0]);
  e2 = bltin_get_ent (args[1]);
  e3 = bltin_get_ent (args[2]);
  e4 = bltin_get_ent (args[3]);

  x1 = (float) class_double (e1);
  x2 = (float) class_double (e2);
  y1 = (float) class_double (e3);
  y2 = (float) class_double (e4);

  cpgrect (x1, x2, y1, y2);

  ent_Clean (e1);
  ent_Clean (e2);
  ent_Clean (e3);
  ent_Clean (e4);

  return ent_Create_Rlab_Success();
}
