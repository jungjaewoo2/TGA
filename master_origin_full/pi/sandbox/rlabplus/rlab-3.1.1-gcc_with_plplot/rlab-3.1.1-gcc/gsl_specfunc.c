// Copyright (C) 2003-2004 Marijan Kostrun
//   part of rlab+4linux project on rlabplus.sourceforge.net
//
// GSL Science Library - Special Functions
// Copyright (C) 1996, 1997, 1998, 1999, 2000, 2001, 2002 Gerard Jungman
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or (at
// your option) any later version.
//
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
//

// rlab headers, located in variable $RLAB_SDK
#include "complex.h"
#include "rlab.h"
#include "ent.h"
#include "class.h"
#include "mem.h"
#include "bltin.h"
#include "util.h"
#include "mdr.h"
#include "mdc.h"
#include "mdcf1.h"
#include "mdcf2.h"
#include "mathl.h"

// gsl headers
// shared object
#include <gsl/gsl_version.h>
#include <gsl/gsl_mode.h>
#include <gsl/gsl_precision.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_machine.h>
#include <gsl/gsl_sf_result.h>
#include <gsl/gsl_const_num.h>


// include: standard C headers
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <math.h>

#define PI M_PI
#define GSL_WARNINGS_OFF

#include "rfileio.h"
#include "rlab_solver_parameters_names.h"

//
// AIRY FUNCTIONS
//

//
// AiryAi
//
#undef THIS_SOLVER
#define THIS_SOLVER "AiryAi"
Ent *
ent_gsl_sf_airy_ai (int nargs, Datum args[])
{
  Ent *e1=0;
  int nr, nc, i;
  MDR *x1, *w;

  gsl_set_error_handler_off ();

  if (nargs == 0)
  {
    fprintf (stdout, THIS_SOLVER ": Airy function Ai. Format:\n");
    fprintf (stdout, THIS_SOLVER ":   AiryAi(x)\n");
    rerror ("No parameters given!");
  }

  //
  // get x
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) != MATRIX_DENSE_REAL)
    rerror (THIS_SOLVER ": argument 'x' has to be MATRIX-DENSE-REAL");

  x1 = class_matrix_real (e1);
  if (!x1)
    rerror (THIS_SOLVER ": argument 'x' has to be MATRIX-DENSE-REAL");

  nr = MNR (x1);
  nc = MNC (x1);
  if (nr * nc == 0)
    rerror (THIS_SOLVER ": argument 'x' has to be MATRIX-DENSE-REAL");

  w = mdr_Create (nr, nc);

  //
  // calculate
  //
  for (i=0; i<nr*nc; i++)
    MdrV0 (w, i) = gsl_sf_airy_Ai (mdrV0 (x1, i), GSL_PREC_DOUBLE);

  ent_Clean (e1);

  return ent_Assign_Rlab_MDR(w);
}

//
// AiryBi
//
#undef THIS_SOLVER
#define THIS_SOLVER "AiryBi"
Ent *
ent_gsl_sf_airy_bi (int nargs, Datum args[])
{
  Ent *e1 = 0;
  int nr, nc, i;
  MDR *x1, *w;

  gsl_set_error_handler_off ();

  if (nargs == 0)
  {
    fprintf (stdout, THIS_SOLVER ": Airy function Bi. Format:\n");
    fprintf (stdout, THIS_SOLVER ":   AiryBi(x)\n");
    rerror ("No parameters given!");
  }
  //
  // get x
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) != MATRIX_DENSE_REAL)
    rerror (THIS_SOLVER ": argument 'x' has to be MATRIX-DENSE-REAL");

  x1 = class_matrix_real (e1);
  if (!x1)
    rerror (THIS_SOLVER ": argument 'x' has to be MATRIX-DENSE-REAL");

  nr = MNR (x1);
  nc = MNC (x1);
  if (nr * nc == 0)
    rerror (THIS_SOLVER ": argument 'x' has to be MATRIX-DENSE-REAL");

  w = mdr_Create (nr, nc);

  //
  // calculate
  //
  for (i=0; i<nr*nc; i++)
    MdrV0 (w, i) = gsl_sf_airy_Bi (mdrV0 (x1, i), GSL_PREC_DOUBLE);

  ent_Clean (e1);

  return ent_Assign_Rlab_MDR(w);
}

//
// AirydAi
//
#undef THIS_SOLVER
#define THIS_SOLVER "AirydAi"
Ent *
ent_gsl_sf_airy_dai (int nargs, Datum args[])
{
  Ent *e1 = 0;
  int nr, nc, i;
  MDR *x1, *w;

  gsl_set_error_handler_off ();

  if (nargs == 0)
  {
    fprintf (stdout, THIS_SOLVER ": Derivative of Airy function Ai. Format:\n");
    fprintf (stdout, THIS_SOLVER ":   AirydAi(x)\n");
    rerror ("No parameters given!");
  }

  //
  // get x
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) != MATRIX_DENSE_REAL)
    rerror (THIS_SOLVER ": argument 'x' has to be MATRIX-DENSE-REAL");

  x1 = class_matrix_real (e1);
  if (!x1)
    rerror (THIS_SOLVER ": argument 'x' has to be MATRIX-DENSE-REAL");

  nr = MNR (x1);
  nc = MNC (x1);
  if (nr * nc == 0)
    rerror (THIS_SOLVER ": argument 'x' has to be MATRIX-DENSE-REAL");

  w = mdr_Create (nr, nc);

  //
  // calculate
  //
  for (i=0; i<nr*nc; i++)
    MdrV0 (w, i) = gsl_sf_airy_Ai_deriv (mdrV0 (x1, i), GSL_PREC_DOUBLE);

  ent_Clean (e1);

  return ent_Assign_Rlab_MDR(w);
}

//
// AirydBi
//
#undef THIS_SOLVER
#define THIS_SOLVER "AirydBi"
Ent *
ent_gsl_sf_airy_dbi (int nargs, Datum args[])
{
  Ent *e1 = 0;
  int nr, nc, i;
  MDR *x1, *w;

  gsl_set_error_handler_off ();

  if (nargs == 0)
  {
    fprintf (stdout, THIS_SOLVER ": Derivative of Airy function Bi. Format:\n");
    fprintf (stdout, THIS_SOLVER ":   AirydBi(x)\n");
    rerror ("No parameters given!");
  }

  //
  // get x
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) != MATRIX_DENSE_REAL)
    rerror (THIS_SOLVER ": argument 'x' has to be MATRIX-DENSE-REAL");

  x1 = class_matrix_real (e1);
  if (!x1)
    rerror (THIS_SOLVER ": argument 'x' has to be MATRIX-DENSE-REAL");

  nr = MNR (x1);
  nc = MNC (x1);
  if (nr * nc == 0)
    rerror (THIS_SOLVER ": argument 'x' has to be MATRIX-DENSE-REAL");

  w = mdr_Create (nr, nc);

  //
  // calculate
  //
  for (i=0; i<nr*nc; i++)
    MdrV0 (w, i) = gsl_sf_airy_Bi_deriv (mdrV0 (x1, i), GSL_PREC_DOUBLE);

  ent_Clean (e1);

  return ent_Assign_Rlab_MDR(w);
}

//
// AiryZeroAi
//
#undef THIS_SOLVER
#define THIS_SOLVER "AiryZeroAi"
Ent *
ent_gsl_sf_airy_zeroai (int nargs, Datum args[])
{
  Ent *e1 = 0;
  int nr, nc, i;
  MDR *x1, *w;

  gsl_set_error_handler_off ();

  if (nargs == 0)
  {
    fprintf (stdout, THIS_SOLVER ": Zero of Airy function Ai. Format:\n");
    fprintf (stdout, THIS_SOLVER ":   AiryZeroAi(i)\n");
    rerror ("No parameters given!");
  }

  // get x
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) != MATRIX_DENSE_REAL)
    rerror (THIS_SOLVER ": argument 'x' has to be MATRIX-DENSE-REAL");

  x1 = class_matrix_real (e1);
  if (!x1)
    rerror (THIS_SOLVER ": argument 'x' has to be MATRIX-DENSE-REAL");

  nr = MNR (x1);
  nc = MNC (x1);
  if (nr * nc == 0)
    rerror (THIS_SOLVER ": argument 'x' has to be MATRIX-DENSE-REAL");

  w = mdr_Create (nr, nc);

  //
  // calculate
  //
  for (i=0; i<nr*nc; i++)
    MdrV0 (w, i) = gsl_sf_airy_zero_Ai (mdrV0 (x1, i));

  ent_Clean (e1);

  return ent_Assign_Rlab_MDR(w);
}

//
// AiryZerodAi
//
Ent *
ent_gsl_sf_airy_zerodai (int nargs, Datum args[])
{
  Ent *e1 = 0;
  int nr, nc, i;
  MDR *x1, *w;
  gsl_set_error_handler_off ();
  if (nargs == 0)
  {
    printf
      ("AiryZerodAi: Zero of first derivative of Airy function Ai. Format:\n");
    printf ("AiryZerodAi:   AiryZerodAi(i)\n");
    rerror ("No parameters given!");
  }
  //
  // get x
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) != MATRIX_DENSE_REAL)
    rerror ("AiryZerodAi: argument 'i' has to be MATRIX-DENSE-REAL");
  x1 = class_matrix_real (e1);
  nr = MNR (x1);
  nc = MNC (x1);
  if (nr * nc == 0)
    rerror ("AiryZerodAi: argument 'i' has to be MATRIX-DENSE-REAL");
  w = mdr_Create (nr, nc);
  //
  // calculate
  //
  for (i = 0; i < nr*nc; i++)
    MdrV0 (w, i) = gsl_sf_airy_zero_Ai_deriv (mdrV0 (x1, i));

  ent_Clean (e1);

  return ent_Assign_Rlab_MDR(w);
}

//
// Airyzerobi
//
Ent *
ent_gsl_sf_airy_zerobi (int nargs, Datum args[])
{
  Ent *e1 = 0;
  int nr, nc, i;
  MDR *x1, *w;
  gsl_set_error_handler_off ();
  if (nargs == 0)
  {
    printf ("Airyzerobi: Zero of Airy function Bi. Format:\n");
    printf ("Airyzerobi:   Airyzerobi(i)\n");
    rerror ("No parameters given!");
  }
  //
  // get x
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) != MATRIX_DENSE_REAL)
    rerror ("Airyzerobi: argument 'x' has to be MATRIX-DENSE-REAL");
  x1 = class_matrix_real (e1);
  nr = MNR (x1);
  nc = MNC (x1);
  if (nr * nc == 0)
    rerror ("Airyzerobi: argument 'x' has to be MATRIX-DENSE-REAL");
  w = mdr_Create (nr, nc);
  //
  // calculate
  //
for (i = 0; i < nr*nc; i++)      MdrV0 (w, i) = gsl_sf_airy_zero_Bi (mdrV0 (x1, i));

    ent_Clean (e1);

  return ent_Assign_Rlab_MDR(w);
}

//
// AiryZerodBi
//
Ent *
ent_gsl_sf_airy_zerodbi (int nargs, Datum args[])
{
  Ent *e1 = 0;
  int nr, nc, i;
  MDR *x1, *w;
  gsl_set_error_handler_off ();
  if (nargs == 0)
  {
    printf
      ("AiryZerodBi: Zero of first derivative of Airy function Bi. Format:\n");
    printf ("AiryZerodBi:   AiryZerodBi(i)\n");
    rerror ("No parameters given!");
  }
  //
  // get x
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) != MATRIX_DENSE_REAL)
    rerror ("AiryZerodBi: argument 'i' has to be MATRIX-DENSE-REAL");
  x1 = class_matrix_real (e1);
  nr = MNR (x1);
  nc = MNC (x1);
  if (nr * nc == 0)
    rerror ("AiryZerodBi: argument 'i' has to be MATRIX-DENSE-REAL");
  w = mdr_Create (nr, nc);
  //
  // calculate
  //
for (i = 0; i < nr*nc; i++)      MdrV0 (w, i) = gsl_sf_airy_zero_Bi_deriv (mdrV0 (x1, i));

    ent_Clean (e1);

  return ent_Assign_Rlab_MDR(w);
}

//
// BESSEL FUNCTIONS
//


//
// BesselI
//
Ent *
ent_gsl_sf_bessel_I (int nargs, Datum args[])
{
  Ent *e1 = 0, *e2 = 0;
  int nr, nc, nr1, nc1, nr2, nc2, i, j, id, jd;
  double x, n;
  MDR *x1, *x2;
  MDR *w;
  gsl_set_error_handler_off ();
  if (nargs == 0)
  {
    printf
      ("BesselI: regular modified cylindrical Bessel function of order 'n'.\n");
    printf ("BesselI: Format:\n");
    printf ("BesselI:   BesselI(x,n)\n");
    rerror ("No parameters given!");
  }
  //
  // get x
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) != MATRIX_DENSE_REAL)
    rerror ("BesselI: argument 'x' has to be MATRIX-DENSE-REAL");
  x1 = class_matrix_real (e1);
  nr1 = MNR (x1);
  nc1 = MNC (x1);
  nr = nr1;
  nc = nc1;
  if (nr1 * nc1 == 0)
    rerror ("BesselI: argument 'x' has to be MATRIX-DENSE-REAL");
  // get n
  e2 = bltin_get_ent (args[1]);
  if (ent_type (e2) != MATRIX_DENSE_REAL)
    rerror ("BesselI: argument 'n' has to be MATRIX-DENSE-REAL");
  x2 = class_matrix_real (e2);
  nr2 = MNR (x2);
  nc2 = MNC (x2);
  if (nr2 * nc2 == 0)
    rerror ("BesselI: argument 'n' has to be MATRIX-DENSE-REAL");
  // figure out the size of output matrix
  if (nc < nc2)
    nc = nc2;
  if (nr < nr2)
    nr = nr2;
  w = mdr_Create (nr, nc);

  // calculate
  for (i = 1; i <= nr; i++)
    for (j = 1; j <= nc; j++)
    {
      id = i;
      jd = j;
      if (id > nr1)
        id = nr1;
      if (jd > nc1)
        jd = nc1;
      x = Mdr1 (x1, id, jd);
      id = i;
      jd = j;
      if (id > nr2)
        id = nr2;
      if (jd > nc2)
        jd = nc2;
      n = Mdr1 (x2, id, jd);
      if (n == 0)
      {
        Mdr1 (w, i, j) = gsl_sf_bessel_I0 (x);
        continue;
      }
      if (n == 1)
      {
        Mdr1 (w, i, j) = gsl_sf_bessel_I1 (x);
        continue;
      }
      if ((n > 1) && (n == (int) n))
      {
        Mdr1 (w, i, j) = gsl_sf_bessel_In (n, x);
        continue;
      }
      Mdr1 (w, i, j) = gsl_sf_bessel_Inu (n, x);
    }

  ent_Clean (e1);
  ent_Clean (e2);

  return ent_Assign_Rlab_MDR(w);
}

//
// BesselJ
//
Ent *
ent_gsl_sf_bessel_J (int nargs, Datum args[])
{
  Ent *e1 = 0, *e2 = 0;
  int nr, nc, nr1, nc1, nr2, nc2, i, j, id, jd;
  double x, n;
  MDR *x1, *x2;
  MDR *w;
  gsl_set_error_handler_off ();
  if (nargs == 0)
  {
    printf ("BesselJ: regular cylindrical Bessel function of order 'n'.\n");
    printf ("BesselJ: Format:\n");
    printf ("BesselJ:   BesselJ(x,n)\n");
    rerror ("No parameters given!");
  }
  //
  // get x
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) != MATRIX_DENSE_REAL)
    rerror ("BesselJ: argument 'x' has to be MATRIX-DENSE-REAL");
  x1 = class_matrix_real (e1);
  nr1 = MNR (x1);
  nc1 = MNC (x1);
  nr = nr1;
  nc = nc1;
  if (nr1 * nc1 == 0)
    rerror ("BesselJ: argument 'x' has to be MATRIX-DENSE-REAL");
  // get n
  e2 = bltin_get_ent (args[1]);
  if (ent_type (e2) != MATRIX_DENSE_REAL)
    rerror ("BesselJ: argument 'n' has to be MATRIX-DENSE-REAL");
  x2 = class_matrix_real (e2);
  nr2 = MNR (x2);
  nc2 = MNC (x2);
  if (nr2 * nc2 == 0)
    rerror ("BesselJ: argument 'n' has to be MATRIX-DENSE-REAL");
  // figure out the size of output matrix
  if (nc < nc2)
    nc = nc2;
  if (nr < nr2)
    nr = nr2;
  w = mdr_Create (nr, nc);
  // calculate
  for (i = 1; i <= nr; i++)
    for (j = 1; j <= nc; j++)
    {
      id = i;
      jd = j;
      if (id > nr1)
        id = nr1;
      if (jd > nc1)
        jd = nc1;
      x = Mdr1 (x1, id, jd);
      id = i;
      jd = j;
      if (id > nr2)
        id = nr2;
      if (jd > nc2)
        jd = nc2;
      n = Mdr1 (x2, id, jd);
      if (n == 0)
      {
        Mdr1 (w, i, j) = gsl_sf_bessel_J0 (x);
        continue;
      }
      if (n == 1)
      {
        Mdr1 (w, i, j) = gsl_sf_bessel_J1 (x);
        continue;
      }
      if ((n > 1) && (n == (int) n))
      {
        Mdr1 (w, i, j) = gsl_sf_bessel_Jn (n, x);
        continue;
      }
      Mdr1 (w, i, j) = gsl_sf_bessel_Jnu (n, x);
    }
    ent_Clean (e1);
    ent_Clean (e2);

    return ent_Assign_Rlab_MDR(w);
}

//
// BesselK
//
Ent *
ent_gsl_sf_bessel_K (int nargs, Datum args[])
{
  Ent *e1 = 0, *e2 = 0;
  int nr, nc, nr1, nc1, nr2, nc2, i, j, id, jd;
  double x, n;
  MDR *x1, *x2;
  MDR *w;
  gsl_set_error_handler_off ();
  if (nargs == 0)
  {
    printf
      ("BesselK: irregular modified cylindrical Bessel function of order 'n'.\n");
    printf ("BesselK: Format:\n");
    printf ("BesselK:   BesselK(x,n)\n");
    rerror ("No parameters given!");
  }
  //
  // get x
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) != MATRIX_DENSE_REAL)
    rerror ("BesselK: argument 'x' has to be MATRIX-DENSE-REAL");
  x1 = class_matrix_real (e1);
  nr1 = MNR (x1);
  nc1 = MNC (x1);
  nr = nr1;
  nc = nc1;
  if (nr1 * nc1 == 0)
    rerror ("BesselK: argument 'x' has to be MATRIX-DENSE-REAL");
  // get n
  e2 = bltin_get_ent (args[1]);
  if (ent_type (e2) != MATRIX_DENSE_REAL)
    rerror ("BesselK: argument 'n' has to be MATRIX-DENSE-REAL");
  x2 = class_matrix_real (e2);
  nr2 = MNR (x2);
  nc2 = MNC (x2);
  if (nr2 * nc2 == 0)
    rerror ("BesselK: argument 'n' has to be MATRIX-DENSE-REAL");
  // figure out the size of output matrix
  if (nc < nc2)
    nc = nc2;
  if (nr < nr2)
    nr = nr2;
  w = mdr_Create (nr, nc);
  // calculate
  for (i = 1; i <= nr; i++)
    for (j = 1; j <= nc; j++)
    {
      id = i;
      jd = j;
      if (id > nr1)
        id = nr1;
      if (jd > nc1)
        jd = nc1;
      x = Mdr1 (x1, id, jd);
      id = i;
      jd = j;
      if (id > nr2)
        id = nr2;
      if (jd > nc2)
        jd = nc2;
      n = Mdr1 (x2, id, jd);
      if (n == 0)
      {
        Mdr1 (w, i, j) = gsl_sf_bessel_K0 (x);
        continue;
      }
      if (n == 1)
      {
        Mdr1 (w, i, j) = gsl_sf_bessel_K1 (x);
        continue;
      }
      if ((n > 1) && (n == (int) n))
      {
        Mdr1 (w, i, j) = gsl_sf_bessel_Kn (n, x);
        continue;
      }
      Mdr1 (w, i, j) = gsl_sf_bessel_Knu (n, x);
    }
    ent_Clean (e1);
    ent_Clean (e2);

    return ent_Assign_Rlab_MDR(w);

}

//
// BesselY
//
Ent *
ent_gsl_sf_bessel_Y (int nargs, Datum args[])
{
  Ent *e1 = 0, *e2 = 0;
  int nr, nc, nr1, nc1, nr2, nc2, i, j, id, jd;
  double x, n;
  MDR *x1, *x2;
  MDR *w;
  gsl_set_error_handler_off ();
  if (nargs == 0)
  {
    printf ("BesselY: irregular cylindrical Bessel function of order 'n'.\n");
    printf ("BesselY: Format:\n");
    printf ("BesselY:   BesselY(x,n)\n");
    rerror ("No parameters given!");
  }
  //
  // get x
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) != MATRIX_DENSE_REAL)
    rerror ("BesselY: argument 'x' has to be MATRIX-DENSE-REAL");
  x1 = class_matrix_real (e1);
  nr1 = MNR (x1);
  nc1 = MNC (x1);
  nr = nr1;
  nc = nc1;
  if (nr1 * nc1 == 0)
    rerror ("BesselY: argument 'x' has to be MATRIX-DENSE-REAL");
  // get n
  e2 = bltin_get_ent (args[1]);
  if (ent_type (e2) != MATRIX_DENSE_REAL)
    rerror ("BesselY: argument 'n' has to be MATRIX-DENSE-REAL");
  x2 = class_matrix_real (e2);
  nr2 = MNR (x2);
  nc2 = MNC (x2);
  if (nr2 * nc2 == 0)
    rerror ("BesselY: argument 'n' has to be MATRIX-DENSE-REAL");
  // figure out the size of output matrix
  if (nc < nc2)
    nc = nc2;
  if (nr < nr2)
    nr = nr2;
  w = mdr_Create (nr, nc);
  // calculate
  for (i = 1; i <= nr; i++)
    for (j = 1; j <= nc; j++)
    {
      id = i;
      jd = j;
      if (id > nr1)
        id = nr1;
      if (jd > nc1)
        jd = nc1;
      x = Mdr1 (x1, id, jd);
      id = i;
      jd = j;
      if (id > nr2)
        id = nr2;
      if (jd > nc2)
        jd = nc2;
      n = Mdr1 (x2, id, jd);
      if (n == 0)
      {
        Mdr1 (w, i, j) = gsl_sf_bessel_Y0 (x);
        continue;
      }
      if (n == 1)
      {
        Mdr1 (w, i, j) = gsl_sf_bessel_Y1 (x);
        continue;
      }
      if ((n > 1) && (n == (int) n))
      {
        Mdr1 (w, i, j) = gsl_sf_bessel_Yn (n, x);
        continue;
      }
      Mdr1 (w, i, j) = gsl_sf_bessel_Ynu (n, x);
    }
    ent_Clean (e1);
    ent_Clean (e2);

    return ent_Assign_Rlab_MDR(w);

}

//
// Besseli
//
Ent *
ent_gsl_sf_bessel_i (int nargs, Datum args[])
{
  Ent *e1 = 0, *e2 = 0;
  int nr, nc, nr1, nc1, nr2, nc2, i, j, id, jd;
  double x;
  int n;
  MDR *x1, *x2;
  MDR *w;
  gsl_set_error_handler_off ();
  if (nargs == 0)
  {
    printf
      ("Besseli: regular modified spherical Bessel function of integer order 'n'.\n");
    printf ("Besseli: Format:\n");
    printf ("Besseli:   Besseli(x,n)\n");
    rerror ("No parameters given!");
  }
  //
  // get x
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) != MATRIX_DENSE_REAL)
    rerror ("Besseli: argument 'x' has to be MATRIX-DENSE-REAL");
  x1 = class_matrix_real (e1);
  nr1 = MNR (x1);
  nc1 = MNC (x1);
  nr = nr1;
  nc = nc1;
  if (nr1 * nc1 == 0)
    rerror ("Besseli: argument 'x' has to be MATRIX-DENSE-REAL");
  // get n
  e2 = bltin_get_ent (args[1]);
  if (ent_type (e2) != MATRIX_DENSE_REAL)
    rerror ("Besseli: argument 'n' has to be MATRIX-DENSE-REAL");
  x2 = class_matrix_real (e2);
  nr2 = MNR (x2);
  nc2 = MNC (x2);
  if (nr2 * nc2 == 0)
    rerror ("Besseli: argument 'n' has to be MATRIX-DENSE-REAL");
  // figure out the size of output matrix
  if (nc < nc2)
    nc = nc2;
  if (nr < nr2)
    nr = nr2;
  w = mdr_Create (nr, nc);
  // calculate
  for (i = 1; i <= nr; i++)
    for (j = 1; j <= nc; j++)
    {
      id = i;
      jd = j;
      if (id > nr1)
        id = nr1;
      if (jd > nc1)
        jd = nc1;
      x = Mdr1 (x1, id, jd);
      id = i;
      jd = j;
      if (id > nr2)
        id = nr2;
      if (jd > nc2)
        jd = nc2;
      n = Mdr1 (x2, id, jd);
      if (n == 0)
      {
        Mdr1 (w, i, j) = gsl_sf_bessel_i0_scaled (x);
        continue;
      }
      if (n == 1)
      {
        Mdr1 (w, i, j) = gsl_sf_bessel_i1_scaled (x);
        continue;
      }
      if (n == 2)
      {
        Mdr1 (w, i, j) = gsl_sf_bessel_i2_scaled (x);
        continue;
      }
      Mdr1 (w, i, j) = gsl_sf_bessel_il_scaled (n, x);
    }
    ent_Clean (e1);
    ent_Clean (e2);

    return ent_Assign_Rlab_MDR(w);

}

//
// Besselj
//
Ent *
ent_gsl_sf_bessel_j (int nargs, Datum args[])
{
  Ent *e1 = 0, *e2 = 0;
  int nr, nc, nr1, nc1, nr2, nc2, i, j, id, jd;
  double x;
  int n;
  MDR *x1, *x2;
  MDR *w;
  gsl_set_error_handler_off ();
  if (nargs == 0)
  {
    printf
      ("Besselj: regular spherical Bessel function of integer order 'n'.\n");
    printf ("Besselj: Format:\n");
    printf ("Besselj:   Besselj(x,n)\n");
    rerror ("No parameters given!");
  }
  //
  // get x
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) != MATRIX_DENSE_REAL)
    rerror ("Besselj: argument 'x' has to be MATRIX-DENSE-REAL");
  x1 = class_matrix_real (e1);
  nr1 = MNR (x1);
  nc1 = MNC (x1);
  nr = nr1;
  nc = nc1;
  if (nr1 * nc1 == 0)
    rerror ("Besselj: argument 'x' has to be MATRIX-DENSE-REAL");
  // get n
  e2 = bltin_get_ent (args[1]);
  if (ent_type (e2) != MATRIX_DENSE_REAL)
    rerror ("Besselj: argument 'n' has to be MATRIX-DENSE-REAL");
  x2 = class_matrix_real (e2);
  nr2 = MNR (x2);
  nc2 = MNC (x2);
  if (nr2 * nc2 == 0)
    rerror ("Besselj: argument 'n' has to be MATRIX-DENSE-REAL");
  // figure out the size of output matrix
  if (nc < nc2)
    nc = nc2;
  if (nr < nr2)
    nr = nr2;
  w = mdr_Create (nr, nc);
  // calculate
  // calculate
  for (i = 1; i <= nr; i++)
    for (j = 1; j <= nc; j++)
    {
      id = i;
      jd = j;
      if (id > nr1)
        id = nr1;
      if (jd > nc1)
        jd = nc1;
      x = Mdr1 (x1, id, jd);
      id = i;
      jd = j;
      if (id > nr2)
        id = nr2;
      if (jd > nc2)
        jd = nc2;
      n = Mdr1 (x2, id, jd);
      if (n == 0)
      {
        Mdr1 (w, i, j) = gsl_sf_bessel_j0 (x);
        continue;
      }
      if (n == 1)
      {
        Mdr1 (w, i, j) = gsl_sf_bessel_j1 (x);
        continue;
      }
      if (n == 2)
      {
        Mdr1 (w, i, j) = gsl_sf_bessel_j2 (x);
        continue;
      }
      Mdr1 (w, i, j) = gsl_sf_bessel_jl (n, x);
    }
    ent_Clean (e1);
    ent_Clean (e2);

    return ent_Assign_Rlab_MDR(w);

}


//
// Besselk
//
Ent *
ent_gsl_sf_bessel_k (int nargs, Datum args[])
{
  Ent *e1 = 0, *e2 = 0;
  int nr, nc, nr1, nc1, nr2, nc2, i, j, id, jd;
  double x;
  int n;
  MDR *x1, *x2;
  MDR *w;
  gsl_set_error_handler_off ();
  if (nargs == 0)
  {
    printf
      ("Besselk: irregular modified spherical Bessel function of integer order 'n'.\n");
    printf ("Besselk: Format:\n");
    printf ("Besselk:   Besselk(x,n)\n");
    rerror ("No parameters given!");
  }
  //
  // get x
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) != MATRIX_DENSE_REAL)
    rerror ("Besselk: argument 'x' has to be MATRIX-DENSE-REAL");
  x1 = class_matrix_real (e1);
  nr1 = MNR (x1);
  nc1 = MNC (x1);
  nr = nr1;
  nc = nc1;
  if (nr1 * nc1 == 0)
    rerror ("Besselk: argument 'x' has to be MATRIX-DENSE-REAL");
  // get n
  e2 = bltin_get_ent (args[1]);
  if (ent_type (e2) != MATRIX_DENSE_REAL)
    rerror ("Besselk: argument 'n' has to be MATRIX-DENSE-REAL");
  x2 = class_matrix_real (e2);
  nr2 = MNR (x2);
  nc2 = MNC (x2);
  if (nr2 * nc2 == 0)
    rerror ("Besselk: argument 'n' has to be MATRIX-DENSE-REAL");
  // figure out the size of output matrix
  if (nc < nc2)
    nc = nc2;
  if (nr < nr2)
    nr = nr2;
  w = mdr_Create (nr, nc);
  // calculate
  for (i = 1; i <= nr; i++)
    for (j = 1; j <= nc; j++)
    {
      id = i;
      jd = j;
      if (id > nr1)
        id = nr1;
      if (jd > nc1)
        jd = nc1;
      x = Mdr1 (x1, id, jd);
      id = i;
      jd = j;
      if (id > nr2)
        id = nr2;
      if (jd > nc2)
        jd = nc2;
      n = Mdr1 (x2, id, jd);
      if (n == 0)
      {
        Mdr1 (w, i, j) = gsl_sf_bessel_k0_scaled (x);
        continue;
      }
      if (n == 1)
      {
        Mdr1 (w, i, j) = gsl_sf_bessel_k1_scaled (x);
        continue;
      }
      if (n == 2)
      {
        Mdr1 (w, i, j) = gsl_sf_bessel_k2_scaled (x);
        continue;
      }
      Mdr1 (w, i, j) = gsl_sf_bessel_kl_scaled (n, x);
    }
    ent_Clean (e1);
    ent_Clean (e2);

    return ent_Assign_Rlab_MDR(w);

}


//
// Bessely
//
Ent *
ent_gsl_sf_bessel_y (int nargs, Datum args[])
{
  Ent *e1 = 0, *e2 = 0;
  int nr, nc, nr1, nc1, nr2, nc2, i, j, id, jd;
  double x;
  int n;
  MDR *x1, *x2;
  MDR *w;
  gsl_set_error_handler_off ();
  if (nargs == 0)
  {
    printf
      ("Bessely: irregular spherical Bessel function of integer order 'n'.\n");
    printf ("Bessely: Format:\n");
    printf ("Bessely:   Bessely(x,n)\n");
    rerror ("No parameters given!");
  }
  //
  // get x
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) != MATRIX_DENSE_REAL)
    rerror ("Bessely: argument 'x' has to be MATRIX-DENSE-REAL");
  x1 = class_matrix_real (e1);
  nr1 = MNR (x1);
  nc1 = MNC (x1);
  nr = nr1;
  nc = nc1;
  if (nr1 * nc1 == 0)
    rerror ("Bessely: argument 'x' has to be MATRIX-DENSE-REAL");
  // get n
  e2 = bltin_get_ent (args[1]);
  if (ent_type (e2) != MATRIX_DENSE_REAL)
    rerror ("Bessely: argument 'n' has to be MATRIX-DENSE-REAL");
  x2 = class_matrix_real (e2);
  nr2 = MNR (x2);
  nc2 = MNC (x2);
  if (nr2 * nc2 == 0)
    rerror ("Bessely: argument 'n' has to be MATRIX-DENSE-REAL");
  // figure out the size of output matrix
  if (nc < nc2)
    nc = nc2;
  if (nr < nr2)
    nr = nr2;
  w = mdr_Create (nr, nc);
  // calculate
  for (i = 1; i <= nr; i++)
    for (j = 1; j <= nc; j++)
    {
      id = i;
      jd = j;
      if (id > nr1)
        id = nr1;
      if (jd > nc1)
        jd = nc1;
      x = Mdr1 (x1, id, jd);
      id = i;
      jd = j;
      if (id > nr2)
        id = nr2;
      if (jd > nc2)
        jd = nc2;
      n = Mdr1 (x2, id, jd);
      if (n == 0)
      {
        Mdr1 (w, i, j) = gsl_sf_bessel_y0 (x);
        continue;
      }
      if (n == 1)
      {
        Mdr1 (w, i, j) = gsl_sf_bessel_y1 (x);
        continue;
      }
      if (n == 2)
      {
        Mdr1 (w, i, j) = gsl_sf_bessel_y2 (x);
        continue;
      }
      Mdr1 (w, i, j) = gsl_sf_bessel_yl (n, x);
    }
    ent_Clean (e1);
    ent_Clean (e2);

    return ent_Assign_Rlab_MDR(w);

}



//
// BesselZeroJ
//
Ent *
ent_gsl_sf_bessel_zeroJ (int nargs, Datum args[])
{
  Ent *e1 = 0, *e2 = 0;
  int nr, nc, nr1, nc1, nr2, nc2, i, j, id, jd;
  int x;
  double n;
  MDR *x1, *x2;
  MDR *w;
  gsl_set_error_handler_off ();
  if (nargs == 0)
  {
    printf
      ("BesselZeroJ: zeros of regular cylindrical Bessel function of order 'n'.\n");
    printf ("BesselZeroJ: Format:\n");
    printf ("BesselZeroJ:   BesselZeroJ(i,n)\n");
    rerror ("No parameters given!");
  }
  //
  // get x
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) != MATRIX_DENSE_REAL)
    rerror ("BesselZeroJ: argument 'i' has to be MATRIX-DENSE-REAL");
  x1 = class_matrix_real (e1);
  nr1 = MNR (x1);
  nc1 = MNC (x1);
  nr = nr1;
  nc = nc1;
  if (nr1 * nc1 == 0)
    rerror ("BesselZeroJ: argument 'i' has to be MATRIX-DENSE-REAL");
  // get n
  e2 = bltin_get_ent (args[1]);
  if (ent_type (e2) != MATRIX_DENSE_REAL)
    rerror ("BesselZeroJ: argument 'n' has to be MATRIX-DENSE-REAL");
  x2 = class_matrix_real (e2);
  nr2 = MNR (x2);
  nc2 = MNC (x2);
  if (nr2 * nc2 == 0)
    rerror ("BesselZeroJ: argument 'n' has to be MATRIX-DENSE-REAL");
  // figure out the size of output matrix
  if (nc < nc2)
    nc = nc2;
  if (nr < nr2)
    nr = nr2;
  w = mdr_Create (nr, nc);
  // calculate
  for (i = 1; i <= nr; i++)
    for (j = 1; j <= nc; j++)
    {
      id = i;
      jd = j;
      if (id > nr1)
        id = nr1;
      if (jd > nc1)
        jd = nc1;
      x = Mdr1 (x1, id, jd);
      id = i;
      jd = j;
      if (id > nr2)
        id = nr2;
      if (jd > nc2)
        jd = nc2;
      n = Mdr1 (x2, id, jd);
      if (n == 0)
      {
        Mdr1 (w, i, j) = gsl_sf_bessel_zero_J0 (x);
        continue;
      }
      if (n == 1)
      {
        Mdr1 (w, i, j) = gsl_sf_bessel_zero_J1 (x);
        continue;
      }
      Mdr1 (w, i, j) = gsl_sf_bessel_zero_Jnu (n, x);
    }
    ent_Clean (e1);
    ent_Clean (e2);

    return ent_Assign_Rlab_MDR(w);

}


//
// Clausen function
//
Ent *
ent_gsl_sf_clausen (int nargs, Datum args[])
{
  Ent *e1 = 0;
  int nr, nc, i;
  MDR *x1, *w;
  gsl_set_error_handler_off ();
  if (nargs == 0)
  {
    printf ("clausen: The Clausen function. Format:\n");
    printf ("clausen:   claussen(x)\n");
    rerror ("No parameters given!");
  }
  //
  // get x
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) != MATRIX_DENSE_REAL)
    rerror ("clausen: argument 'x' has to be MATRIX-DENSE-REAL");
  x1 = class_matrix_real (e1);
  nr = MNR (x1);
  nc = MNC (x1);
  if (nr * nc == 0)
    rerror ("clausen: argument 'x' has to be MATRIX-DENSE-REAL");
  w = mdr_Create (nr, nc);
  //
  // calculate
  //
for (i = 0; i < nr*nc; i++)      MdrV0 (w, i) = gsl_sf_clausen (mdrV0 (x1, i));

    ent_Clean (e1);

  return ent_Assign_Rlab_MDR(w);
}

//
// Hydrogen
//
Ent *
ent_gsl_sf_hydrogen (int nargs, Datum args[])
{
  Ent *e1 = 0, *e2 = 0, *e3 = 0, *e4 = 0;
  int nr, nc, nr1, nc1, nr2, nc2, nr3, nc3, nr4, nc4, i, j, id, jd;
  int n, l;
  double Z, r;
  MDR *x1, *x2, *x3, *x4;
  MDR *w;
  gsl_set_error_handler_off ();
  if (nargs == 0)
  {
    printf
      ("Hydrogen: Normalized radial wave function for the Hydrogen bound state,\n");
    printf ("Hydrogen: with parameters 'n', principal quantum number,\n");
    printf ("Hydrogen: 'l', angular momentum quantum number, 'Z' atomic\n");
    printf ("Hydrogen: number and 'r', radial distance.\n");
    printf ("Hydrogen: Format:\n");
    printf ("Hydrogen:   Hydrogen(n,l,Z,r)\n");
    rerror ("No parameters given!");
  }
  //
  // get n
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) != MATRIX_DENSE_REAL)
    rerror ("Hydrogen: argument 'n' has to be MATRIX-DENSE-REAL");
  x1 = class_matrix_real (e1);
  nr1 = MNR (x1);
  nc1 = MNC (x1);
  nr = nr1;
  nc = nc1;
  if (nr1 * nc1 == 0)
    rerror ("Hydrogen: argument 'k' has to be MATRIX-DENSE-REAL");
  //
  // get l
  //
  e2 = bltin_get_ent (args[1]);
  if (ent_type (e2) != MATRIX_DENSE_REAL)
    rerror ("Hydrogen: argument 'l' has to be MATRIX-DENSE-REAL");
  x2 = class_matrix_real (e2);
  nr2 = MNR (x2);
  nc2 = MNC (x2);
  if (nr2 * nc2 == 0)
    rerror ("Hydrogen: argument 'l' has to be MATRIX-DENSE-REAL");
  // figure out the size of output matrix
  if (nc < nc2)
    nc = nc2;
  if (nr < nr2)
    nr = nr2;
  //
  // get n2
  //
  e3 = bltin_get_ent (args[2]);
  if (ent_type (e3) != MATRIX_DENSE_REAL)
    rerror ("Hydrogen: argument 'Z' has to be MATRIX-DENSE-REAL");
  x3 = class_matrix_real (e3);
  nr3 = MNR (x3);
  nc3 = MNC (x3);
  if (nr3 * nc3 == 0)
    rerror ("Hydrogen: argument 'Z' has to be MATRIX-DENSE-REAL");
  // figure out the size of output matrix
  if (nc < nc3)
    nc = nc3;
  if (nr < nr3)
    nr = nr3;
  //
  // get t
  //
  e4 = bltin_get_ent (args[3]);
  if (ent_type (e4) != MATRIX_DENSE_REAL)
    rerror ("Hydrogen: argument 'r' has to be MATRIX-DENSE-REAL");
  x4 = class_matrix_real (e4);
  nr4 = MNR (x4);
  nc4 = MNC (x4);
  if (nr4 * nc4 == 0)
    rerror ("Hydrogen: argument 'r' has to be MATRIX-DENSE-REAL");
  // figure out the size of output matrix
  if (nc < nc4)
    nc = nc4;
  if (nr < nr4)
    nr = nr4;
  w = mdr_Create (nr, nc);
  // calculate
  for (i = 1; i <= nr; i++)
    for (j = 1; j <= nc; j++)
    {
      id = i;
      jd = j;
      if (id > nr1)
        id = nr1;
      if (jd > nc1)
        jd = nc1;
      n = Mdr1 (x1, id, jd);
      id = i;
      jd = j;
      if (id > nr2)
        id = nr2;
      if (jd > nc2)
        jd = nc2;
      l = Mdr1 (x2, id, jd);
      id = i;
      jd = j;
      if (id > nr3)
        id = nr3;
      if (jd > nc3)
        jd = nc3;
      Z = Mdr1 (x3, id, jd);
      id = i;
      jd = j;
      if (id > nr4)
        id = nr4;
      if (jd > nc4)
        jd = nc4;
      r = Mdr1 (x4, id, jd);
      if (n == 1)
      {
        Mdr1 (w, i, j) = gsl_sf_hydrogenicR_1 (Z, r);
        continue;
      }
      Mdr1 (w, i, j) = gsl_sf_hydrogenicR (n, l, Z, r);
    }
  ent_Clean (e1);
  ent_Clean (e2);
  ent_Clean (e3);
  ent_Clean (e4);

  return ent_Assign_Rlab_MDR(w);

}

//
// Dawson
//
Ent *
ent_gsl_sf_dawson (int nargs, Datum args[])
{
  Ent *e1 = 0;
  int nr, nc, i;
  MDR *x1, *w;
  gsl_set_error_handler_off ();
  if (nargs == 0)
  {
    printf ("Dawson: The Dawson function. Format:\n");
    printf ("Dawson:   Dawson(x)\n");
    rerror ("No parameters given!");
  }
  //
  // get x
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) != MATRIX_DENSE_REAL)
    rerror ("Dawson: argument 'x' has to be MATRIX-DENSE-REAL");
  x1 = class_matrix_real (e1);
  nr = MNR (x1);
  nc = MNC (x1);
  if (nr * nc == 0)
    rerror ("Dawson: argument 'x' has to be MATRIX-DENSE-REAL");
  w = mdr_Create (nr, nc);
  //
  // calculate
  //
for (i = 0; i < nr*nc; i++)      MdrV0 (w, i) = gsl_sf_dawson (mdrV0 (x1, i));

    ent_Clean (e1);


  return ent_Assign_Rlab_MDR(w);

}

//
// Debye functions of the order 1 through 4
//
Ent *
ent_gsl_sf_debye (int nargs, Datum args[])
{
  Ent *e1 = 0, *e2 = 0;
  int nr, nc, i;
  int n;
  MDR *x1=0, *w=0;

  gsl_set_error_handler_off ();

  if (nargs == 0)
  {
    printf ("Debye: The Debye function of integer order n=1..4. Format:\n");
    printf ("Debye:   Debye(x,n)\n");
    rerror ("No parameters given!");
  }
  //
  // get x
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) != MATRIX_DENSE_REAL)
    rerror ("Debye: argument 'x' has to be MATRIX-DENSE-REAL");
  x1 = class_matrix_real (e1);
  nr = MNR (x1);
  nc = MNC (x1);
  if (nr * nc == 0)
    rerror ("Debye: argument 'x' has to be MATRIX-DENSE-REAL");
  //
  // get n
  //
  e2 = bltin_get_ent (args[1]);
  if (ent_type (e2) != MATRIX_DENSE_REAL)
    rerror ("Debye: argument 'n' has to an integer 1..4");
  n = class_double (e2);
  w = mdr_Create (nr, nc);
  //
  // calculate
  //
  for (i = 0; i < nr*nc; i++)
  {
    MdrV0 (w, i) = GSL_NAN;
    if (n == 1)
      MdrV0 (w, i) = gsl_sf_debye_1 (mdrV0 (x1, i));
    if (n == 2)
      MdrV0 (w, i) = gsl_sf_debye_2 (mdrV0 (x1, i));
    if (n == 3)
      MdrV0 (w, i) = gsl_sf_debye_3 (mdrV0 (x1, i));
    if (n == 4)
      MdrV0 (w, i) = gsl_sf_debye_4 (mdrV0 (x1, i));
  }

  ent_Clean (e1);
  ent_Clean (e2);

  return ent_Assign_Rlab_MDR(w);

}

//
// EllipticE
//
Ent *
ent_gsl_sf_elliptic_E (int nargs, Datum args[])
{
  Ent *e1 = 0, *e2 = 0;
  int nr, nc, nr1, nc1, nr2, nc2, i, j, id, jd;
  double k, fi;
  MDR *x1, *x2;
  MDR *w;
  gsl_set_error_handler_off ();
  if (nargs == 0)
  {
    printf ("EllipticE: Elliptic function E of order 'k'.\n");
    printf ("EllipticE: (1) Format:\n");
    printf ("EllipticE:   EllipticE(k)\n");
    printf ("EllipticE: (2) Format:\n");
    printf ("EllipticE:   EllipticE(k,fi)\n");
    rerror ("No parameters given!");
  }
  //
  // get e1
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) != MATRIX_DENSE_REAL)
    rerror ("EllipticE: argument 'k' has to be MATRIX-DENSE-REAL");
  x1 = class_matrix_real (e1);
  nr1 = MNR (x1);
  nc1 = MNC (x1);
  nr = nr1;
  nc = nc1;
  if (nr1 * nc1 == 0)
    rerror ("EllipticE: argument 'k' has to be MATRIX-DENSE-REAL");
  if (nargs == 2)
  {
    // get e2
    e2 = bltin_get_ent (args[1]);
    if (ent_type (e2) != MATRIX_DENSE_REAL)
      rerror ("EllipticE: argument 'fi' has to be MATRIX-DENSE-REAL");
    x2 = class_matrix_real (e2);
    nr2 = MNR (x2);
    nc2 = MNC (x2);
    if (nr2 * nc2 == 0)
      rerror ("EllipticE: argument 'fi' has to be MATRIX-DENSE-REAL");
    // figure out the size of output matrix
    if (nc < nc2)
      nc = nc2;
    if (nr < nr2)
      nr = nr2;
  }
  w = mdr_Create (nr, nc);
  // calculate
  if (nargs == 1)
    for (i = 1; i <= nr; i++)
      for (j = 1; j <= nc; j++)
      {
        id = i;
        jd = j;
        if (id > nr1)
          id = nr1;
        if (jd > nc1)
          jd = nc1;
        k = Mdr1 (x1, id, jd);
        MdrV0 (w, i) = gsl_sf_ellint_Ecomp (k, GSL_PREC_DOUBLE);
      }
  if (nargs == 2)
    for (i = 1; i <= nr; i++)
      for (j = 1; j <= nc; j++)
      {
        id = i;
        jd = j;
        if (id > nr1)
          id = nr1;
        if (jd > nc1)
          jd = nc1;
        k = Mdr1 (x1, id, jd);
        id = i;
        jd = j;
        if (id > nr2)
          id = nr2;
        if (jd > nc2)
          jd = nc2;
        fi = Mdr1 (x2, id, jd);
        MdrV0 (w, i) = gsl_sf_ellint_E (fi, k, GSL_PREC_DOUBLE);
      }
      ent_Clean (e1);
      ent_Clean (e2);

      return ent_Assign_Rlab_MDR(w);

}

//
// EllipticF
//
Ent *
ent_gsl_sf_elliptic_F (int nargs, Datum args[])
{
  Ent *e1 = 0, *e2 = 0;
  int nr, nc, nr1, nc1, nr2, nc2, i, j, id, jd;
  double k, fi;
  MDR *x1, *x2;
  MDR *w;
  gsl_set_error_handler_off ();
  if (nargs == 0)
  {
    printf ("EllipticF: Elliptic function F of order 'k'. Format:\n");
    printf ("EllipticF:   EllipticF(k,fi)\n");
    rerror ("No parameters given!");
  }
  //
  // get e1
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) != MATRIX_DENSE_REAL)
    rerror ("EllipticF: argument 'k' has to be MATRIX-DENSE-REAL");
  x1 = class_matrix_real (e1);
  nr1 = MNR (x1);
  nc1 = MNC (x1);
  nr = nr1;
  nc = nc1;
  if (nr1 * nc1 == 0)
    rerror ("EllipticF: argument 'k' has to be MATRIX-DENSE-REAL");
  // get e2
  e2 = bltin_get_ent (args[1]);
  if (ent_type (e2) != MATRIX_DENSE_REAL)
    rerror ("EllipticE: argument 'fi' has to be MATRIX-DENSE-REAL");
  x2 = class_matrix_real (e2);
  nr2 = MNR (x2);
  nc2 = MNC (x2);
  if (nr2 * nc2 == 0)
    rerror ("EllipticE: argument 'fi' has to be MATRIX-DENSE-REAL");
  // figure out the size of output matrix
  if (nc < nc2)
    nc = nc2;
  if (nr < nr2)
    nr = nr2;
  w = mdr_Create (nr, nc);
  // calculate
  for (i = 1; i <= nr; i++)
    for (j = 1; j <= nc; j++)
    {
      id = i;
      jd = j;
      if (id > nr1)
        id = nr1;
      if (jd > nc1)
        jd = nc1;
      k = Mdr1 (x1, id, jd);
      id = i;
      jd = j;
      if (id > nr2)
        id = nr2;
      if (jd > nc2)
        jd = nc2;
      fi = Mdr1 (x2, id, jd);
      Mdr1 (w, i, j) = gsl_sf_ellint_F (fi, k, GSL_PREC_DOUBLE);
    }
    ent_Clean (e1);
    ent_Clean (e2);

    return ent_Assign_Rlab_MDR(w);

}

//
// EllipticK
//
Ent *
ent_gsl_sf_elliptic_K (int nargs, Datum args[])
{
  Ent *e1 = 0;
  int nr, nc, nr1, nc1, i, j, id, jd;
  double k;
  MDR *x1;
  MDR *w;
  gsl_set_error_handler_off ();
  if (nargs == 0)
  {
    printf ("EllipticK: Elliptic function K of order 'k'. Format:\n");
    printf ("EllipticK:   EllipticK(k)\n");
    rerror ("No parameters given!");
  }
  //
  // get e1
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) != MATRIX_DENSE_REAL)
    rerror ("EllipticK: argument 'k' has to be MATRIX-DENSE-REAL");
  x1 = class_matrix_real (e1);
  nr1 = MNR (x1);
  nc1 = MNC (x1);
  nr = nr1;
  nc = nc1;
  if (nr1 * nc1 == 0)
    rerror ("EllipticK: argument 'k' has to be MATRIX-DENSE-REAL");
  w = mdr_Create (nr, nc);
  // calculate
  for (i = 1; i <= nr; i++)
    for (j = 1; j <= nc; j++)
    {
      id = i;
      jd = j;
      if (id > nr1)
        id = nr1;
      if (jd > nc1)
        jd = nc1;
      k = Mdr1 (x1, id, jd);
      Mdr1 (w, i, j) = gsl_sf_ellint_Kcomp (k, GSL_PREC_DOUBLE);
    }
    ent_Clean (e1);

    return ent_Assign_Rlab_MDR(w);

}

//
// EllipticP
//
Ent *
ent_gsl_sf_elliptic_P (int nargs, Datum args[])
{
  Ent *e1 = 0, *e2 = 0, *e3 = 0;
  int nr, nc, nr1, nc1, nr2, nc2, nr3, nc3, i, j, id, jd;
  double k, fi, n;
  MDR *x1, *x2, *x3;
  MDR *w;
  gsl_set_error_handler_off ();
  if (nargs == 0)
  {
    printf ("EllipticP: Incomplete elliptic integral. Format:\n");
    printf ("EllipticP:   EllipticP(k,fi,n)\n");
    rerror ("No parameters given!");
  }
  //
  // get e1
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) != MATRIX_DENSE_REAL)
    rerror ("EllipticP: argument 'k' has to be MATRIX-DENSE-REAL");
  x1 = class_matrix_real (e1);
  nr1 = MNR (x1);
  nc1 = MNC (x1);
  nr = nr1;
  nc = nc1;
  if (nr1 * nc1 == 0)
    rerror ("EllipticP: argument 'k' has to be MATRIX-DENSE-REAL");
  // get e2
  e2 = bltin_get_ent (args[1]);
  if (ent_type (e2) != MATRIX_DENSE_REAL)
    rerror ("EllipticP: argument 'fi' has to be MATRIX-DENSE-REAL");
  x2 = class_matrix_real (e2);
  nr2 = MNR (x2);
  nc2 = MNC (x2);
  if (nr2 * nc2 == 0)
    rerror ("EllipticP: argument 'fi' has to be MATRIX-DENSE-REAL");
  // figure out the size of output matrix
  if (nc < nc2)
    nc = nc2;
  if (nr < nr2)
    nr = nr2;
  // get e3
  e3 = bltin_get_ent (args[2]);
  if (ent_type (e3) != MATRIX_DENSE_REAL)
    rerror ("EllipticP: argument 'n' has to be MATRIX-DENSE-REAL");
  x3 = class_matrix_real (e3);
  nr3 = MNR (x3);
  nc3 = MNC (x3);
  if (nr3 * nc3 == 0)
    rerror ("EllipticP: argument 'fi' has to be MATRIX-DENSE-REAL");
  // figure out the size of output matrix
  if (nc < nc3)
    nc = nc3;
  if (nr < nr3)
    nr = nr3;
  w = mdr_Create (nr, nc);
  // calculate
  for (i = 1; i <= nr; i++)
    for (j = 1; j <= nc; j++)
    {
      id = i;
      jd = j;
      if (id > nr1)
        id = nr1;
      if (jd > nc1)
        jd = nc1;
      k = Mdr1 (x1, id, jd);
      id = i;
      jd = j;
      if (id > nr2)
        id = nr2;
      if (jd > nc2)
        jd = nc2;
      fi = Mdr1 (x2, id, jd);
      id = i;
      jd = j;
      if (id > nr3)
        id = nr3;
      if (jd > nc3)
        jd = nc3;
      n = Mdr1 (x3, id, jd);
      Mdr1 (w, i, j) = gsl_sf_ellint_P (fi, k, n, GSL_PREC_DOUBLE);
    }
    ent_Clean (e1);
  ent_Clean (e2);
  ent_Clean (e3);

  return ent_Assign_Rlab_MDR(w);

}

//
// EllipticD
//
#undef  THIS_SOLVER
#define THIS_SOLVER "EllipticD"
Ent * ent_gsl_sf_elliptic_D (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0;
  int nr, nc, i, j;
  MDR *x1=0, *x2=0;
  MDR *w=0;

  FILE *rlab_stderr = (!RLAB_STDERR_DS) ? stderr : RLAB_STDERR_DS;

  gsl_set_error_handler_off ();
  if (nargs == 0)
  {
    fprintf (rlab_stderr, THIS_SOLVER ": Incomplete elliptic integral. Format:\n");
    fprintf (rlab_stderr, THIS_SOLVER ":   EllipticD(k,fi)\n");
    rerror ("No parameters given!");
  }

  //
  // get e1
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) != MATRIX_DENSE_REAL)
    rerror (THIS_SOLVER ": " RLAB_ERROR_ARG1_MDR_MATRIX);
  x1 = ent_data (e1);
  if (SIZE(x1)<1)
    rerror (THIS_SOLVER ": " RLAB_ERROR_ARG1_MDR_MATRIX);

  // get e2
  e2 = bltin_get_ent (args[1]);
  if (ent_type (e2) != MATRIX_DENSE_REAL)
    rerror (THIS_SOLVER ": " RLAB_ERROR_ARG2_MDR_MATRIX);
  x2 = ent_data (e2);
  if (SIZE(x2)<1)
    rerror (THIS_SOLVER ": " RLAB_ERROR_ARG2_MDR_MATRIX);

  nr = MAX(MNR (x1),MNR (x2));
  nc = MAX(MNC (x1),MNC (x2));

  w = mdr_Create (nr, nc);
  // calculate
  for (i=0; i<nr; i++)
  {
    for (j=0; j<nc; j++)
    {
#if ((GSL_MAJOR_VERSION == 1) && (GSL_MINOR_VERSION <= 16))
      Mdr0 (w, i, j) = gsl_sf_ellint_D (mdr0_safe (x2, i, j), mdr0_safe (x1, i, j), 0, GSL_PREC_DOUBLE);
#else
      Mdr0 (w, i, j) = gsl_sf_ellint_D (mdr0_safe (x2, i, j), mdr0_safe (x1, i, j), GSL_PREC_DOUBLE);
#endif
    }
  }
  ent_Clean (e1);
  ent_Clean (e2);

  return ent_Assign_Rlab_MDR(w);
}

//
// EllipticRC
//
Ent *
ent_gsl_sf_elliptic_RC (int nargs, Datum args[])
{
  Ent *e1 = 0, *e2 = 0;
  int nr, nc, nr1, nc1, nr2, nc2, i, j, id, jd;
  double x, y;
  MDR *x1, *x2;
  MDR *w;
  gsl_set_error_handler_off ();
  if (nargs == 0)
  {
    printf
      ("EllipticRC: Carlson form of incomplete elliptic integral RC. Format:\n");
    printf ("EllipticRC:   EllipticRC(x,y)\n");
    rerror ("No parameters given!");
  }
  //
  // get e1
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) != MATRIX_DENSE_REAL)
    rerror ("EllipticRC: argument 'x' has to be MATRIX-DENSE-REAL");
  x1 = class_matrix_real (e1);
  nr1 = MNR (x1);
  nc1 = MNC (x1);
  nr = nr1;
  nc = nc1;
  if (nr1 * nc1 == 0)
    rerror ("EllipticRC: argument 'x' has to be MATRIX-DENSE-REAL");
  // get e2
  e2 = bltin_get_ent (args[1]);
  if (ent_type (e2) != MATRIX_DENSE_REAL)
    rerror ("EllipticRC: argument 'y' has to be MATRIX-DENSE-REAL");
  x2 = class_matrix_real (e2);
  nr2 = MNR (x2);
  nc2 = MNC (x2);
  if (nr2 * nc2 == 0)
    rerror ("EllipticRC: argument 'y' has to be MATRIX-DENSE-REAL");
  // figure out the size of output matrix
  if (nc < nc2)
    nc = nc2;
  if (nr < nr2)
    nr = nr2;
  w = mdr_Create (nr, nc);
  // calculate
  for (i = 1; i <= nr; i++)
    for (j = 1; j <= nc; j++)
    {
      id = i;
      jd = j;
      if (id > nr1)
        id = nr1;
      if (jd > nc1)
        jd = nc1;
      x = Mdr1 (x1, id, jd);
      id = i;
      jd = j;
      if (id > nr2)
        id = nr2;
      if (jd > nc2)
        jd = nc2;
      y = Mdr1 (x2, id, jd);
      Mdr1 (w, i, j) = gsl_sf_ellint_RC (x, y, GSL_PREC_DOUBLE);
    }
    ent_Clean (e1);
    ent_Clean (e2);

    return ent_Assign_Rlab_MDR(w);

}

//
// EllipticRD
//
Ent *
ent_gsl_sf_elliptic_RD (int nargs, Datum args[])
{
  Ent *e1 = 0, *e2 = 0, *e3 = 0;
  int nr, nc, nr1, nc1, nr2, nc2, nr3, nc3, i, j, id, jd;
  double x, y, z;
  MDR *x1, *x2, *x3;
  MDR *w;
  gsl_set_error_handler_off ();
  if (nargs == 0)
  {
    printf
      ("EllipticRD: Carlson form of incomplete elliptic integral RD. Format:\n");
    printf ("EllipticRD:   EllipticRD(x,y,z)\n");
    rerror ("No parameters given!");
  }
  //
  // get e1
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) != MATRIX_DENSE_REAL)
    rerror ("EllipticRD: argument 'x' has to be MATRIX-DENSE-REAL");
  x1 = class_matrix_real (e1);
  nr1 = MNR (x1);
  nc1 = MNC (x1);
  nr = nr1;
  nc = nc1;
  if (nr1 * nc1 == 0)
    rerror ("EllipticRD: argument 'x' has to be MATRIX-DENSE-REAL");
  // get e2
  e2 = bltin_get_ent (args[1]);
  if (ent_type (e2) != MATRIX_DENSE_REAL)
    rerror ("EllipticRD: argument 'y' has to be MATRIX-DENSE-REAL");
  x2 = class_matrix_real (e2);
  nr2 = MNR (x2);
  nc2 = MNC (x2);
  if (nr2 * nc2 == 0)
    rerror ("EllipticRD: argument 'y' has to be MATRIX-DENSE-REAL");
  // figure out the size of output matrix
  if (nc < nc2)
    nc = nc2;
  if (nr < nr2)
    nr = nr2;
  // get e3
  e3 = bltin_get_ent (args[2]);
  if (ent_type (e3) != MATRIX_DENSE_REAL)
    rerror ("EllipticRD: argument 'z' has to be MATRIX-DENSE-REAL");
  x3 = class_matrix_real (e3);
  nr3 = MNR (x3);
  nc3 = MNC (x3);
  if (nr3 * nc3 == 0)
    rerror ("EllipticRD: argument 'z' has to be MATRIX-DENSE-REAL");
  // figure out the size of output matrix
  if (nc < nc3)
    nc = nc3;
  if (nr < nr3)
    nr = nr3;
  w = mdr_Create (nr, nc);
  // calculate
  // calculate
  for (i = 1; i <= nr; i++)
    for (j = 1; j <= nc; j++)
    {
      id = i;
      jd = j;
      if (id > nr1)
        id = nr1;
      if (jd > nc1)
        jd = nc1;
      x = Mdr1 (x1, id, jd);
      id = i;
      jd = j;
      if (id > nr2)
        id = nr2;
      if (jd > nc2)
        jd = nc2;
      y = Mdr1 (x2, id, jd);
      id = i;
      jd = j;
      if (id > nr3)
        id = nr3;
      if (jd > nc3)
        jd = nc3;
      z = Mdr1 (x3, id, jd);
      Mdr1 (w, i, j) = gsl_sf_ellint_RD (x, y, z, GSL_PREC_DOUBLE);
    }

  ent_Clean (e1);
  ent_Clean (e2);
  ent_Clean (e3);

  return ent_Assign_Rlab_MDR(w);

}

//
// EllipticRF
//
Ent *
ent_gsl_sf_elliptic_RF (int nargs, Datum args[])
{
  Ent *e1 = 0, *e2 = 0, *e3 = 0;
  int nr, nc, nr1, nc1, nr2, nc2, nr3, nc3, i, j, id, jd;
  double x, y, z;
  MDR *x1, *x2, *x3;
  MDR *w;
  gsl_set_error_handler_off ();
  if (nargs == 0)
  {
    printf
      ("EllipticRF: Carlson form of incomplete elliptic integral RF. Format:\n");
    printf ("EllipticRF:   EllipticRF(x,y,z)\n");
    rerror ("No parameters given!");
  }
  //
  // get e1
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) != MATRIX_DENSE_REAL)
    rerror ("EllipticRF: argument 'x' has to be MATRIX-DENSE-REAL");
  x1 = class_matrix_real (e1);
  nr1 = MNR (x1);
  nc1 = MNC (x1);
  nr = nr1;
  nc = nc1;
  if (nr1 * nc1 == 0)
    rerror ("EllipticRF: argument 'x' has to be MATRIX-DENSE-REAL");
  // get e2
  e2 = bltin_get_ent (args[1]);
  if (ent_type (e2) != MATRIX_DENSE_REAL)
    rerror ("EllipticRF: argument 'y' has to be MATRIX-DENSE-REAL");
  x2 = class_matrix_real (e2);
  nr2 = MNR (x2);
  nc2 = MNC (x2);
  if (nr2 * nc2 == 0)
    rerror ("EllipticRF: argument 'y' has to be MATRIX-DENSE-REAL");
  // figure out the size of output matrix
  if (nc < nc2)
    nc = nc2;
  if (nr < nr2)
    nr = nr2;
  // get e3
  e3 = bltin_get_ent (args[2]);
  if (ent_type (e3) != MATRIX_DENSE_REAL)
    rerror ("EllipticRF: argument 'z' has to be MATRIX-DENSE-REAL");
  x3 = class_matrix_real (e3);
  nr3 = MNR (x3);
  nc3 = MNC (x3);
  if (nr3 * nc3 == 0)
    rerror ("EllipticRF: argument 'z' has to be MATRIX-DENSE-REAL");
  // figure out the size of output matrix
  if (nc < nc3)
    nc = nc3;
  if (nr < nr3)
    nr = nr3;
  w = mdr_Create (nr, nc);
  // calculate
  for (i = 1; i <= nr; i++)
    for (j = 1; j <= nc; j++)
    {
      id = i;
      jd = j;
      if (id > nr1)
        id = nr1;
      if (jd > nc1)
        jd = nc1;
      x = Mdr1 (x1, id, jd);
      id = i;
      jd = j;
      if (id > nr2)
        id = nr2;
      if (jd > nc2)
        jd = nc2;
      y = Mdr1 (x2, id, jd);
      id = i;
      jd = j;
      if (id > nr3)
        id = nr3;
      if (jd > nc3)
        jd = nc3;
      z = Mdr1 (x3, id, jd);
      Mdr1 (w, i, j) = gsl_sf_ellint_RF (x, y, z, GSL_PREC_DOUBLE);
    }
    ent_Clean (e1);
    ent_Clean (e2);
    ent_Clean (e3);

    return ent_Assign_Rlab_MDR(w);

}

//
// EllipticRJ
//
Ent *
ent_gsl_sf_elliptic_RJ (int nargs, Datum args[])
{
  Ent *e1 = 0, *e2 = 0, *e3 = 0, *e4 = 0;
  int nr, nc, nr1, nc1, nr2, nc2, nr3, nc3, nr4, nc4, i, j, id, jd;
  double x, y, z, p;
  MDR *x1, *x2, *x3, *x4;
  MDR *w;
  gsl_set_error_handler_off ();
  if (nargs == 0)
  {
    printf
      ("EllipticRJ: Carlson form of incomplete elliptic integral RJ. Format:\n");
    printf ("EllipticRJ:   EllipticRJ(x,y,z,p)\n");
    rerror ("No parameters given!");
  }
  //
  // get e1
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) != MATRIX_DENSE_REAL)
    rerror ("EllipticRJ: argument 'x' has to be MATRIX-DENSE-REAL");
  x1 = class_matrix_real (e1);
  nr1 = MNR (x1);
  nc1 = MNC (x1);
  nr = nr1;
  nc = nc1;
  if (nr1 * nc1 == 0)
    rerror ("EllipticRJ: argument 'x' has to be MATRIX-DENSE-REAL");
  // get e2
  e2 = bltin_get_ent (args[1]);
  if (ent_type (e2) != MATRIX_DENSE_REAL)
    rerror ("EllipticRJ: argument 'y' has to be MATRIX-DENSE-REAL");
  x2 = class_matrix_real (e2);
  nr2 = MNR (x2);
  nc2 = MNC (x2);
  if (nr2 * nc2 == 0)
    rerror ("EllipticRJ: argument 'y' has to be MATRIX-DENSE-REAL");
  // figure out the size of output matrix
  if (nc < nc2)
    nc = nc2;
  if (nr < nr2)
    nr = nr2;
  // get e3
  e3 = bltin_get_ent (args[2]);
  if (ent_type (e3) != MATRIX_DENSE_REAL)
    rerror ("EllipticRJ: argument 'z' has to be MATRIX-DENSE-REAL");
  x3 = class_matrix_real (e3);
  nr3 = MNR (x3);
  nc3 = MNC (x3);
  if (nr3 * nc3 == 0)
    rerror ("EllipticRJ: argument 'z' has to be MATRIX-DENSE-REAL");
  // figure out the size of output matrix
  if (nc < nc3)
    nc = nc3;
  if (nr < nr3)
    nr = nr3;
  // get e4
  e4 = bltin_get_ent (args[3]);
  if (ent_type (e4) != MATRIX_DENSE_REAL)
    rerror ("EllipticRJ: argument 'p' has to be MATRIX-DENSE-REAL");
  x4 = class_matrix_real (e4);
  nr4 = MNR (x4);
  nc4 = MNC (x4);
  if (nr4 * nc4 == 0)
    rerror ("EllipticRJ: argument 'p' has to be MATRIX-DENSE-REAL");
  // figure out the size of output matrix
  if (nc < nc4)
    nc = nc4;
  if (nr < nr4)
    nr = nr4;
  w = mdr_Create (nr, nc);
  // calculate
  for (i = 1; i <= nr; i++)
    for (j = 1; j <= nc; j++)
    {
      id = i;
      jd = j;
      if (id > nr1)
        id = nr1;
      if (jd > nc1)
        jd = nc1;
      x = Mdr1 (x1, id, jd);
      id = i;
      jd = j;
      if (id > nr2)
        id = nr2;
      if (jd > nc2)
        jd = nc2;
      y = Mdr1 (x2, id, jd);
      id = i;
      jd = j;
      if (id > nr3)
        id = nr3;
      if (jd > nc3)
        jd = nc3;
      z = Mdr1 (x3, id, jd);
      id = i;
      jd = j;
      if (id > nr4)
        id = nr4;
      if (jd > nc4)
        jd = nc4;
      p = Mdr1 (x4, id, jd);
      Mdr1 (w, i, j) = gsl_sf_ellint_RJ (x, y, z, p, GSL_PREC_DOUBLE);
    }

  ent_Clean (e1);
  ent_Clean (e2);
  ent_Clean (e3);
  ent_Clean (e4);

  return ent_Assign_Rlab_MDR(w);
}

//
// EllipticJacobi
//
Ent *
ent_gsl_sf_elliptic_jacobi (int nargs, Datum args[])
{
  Ent *e1 = 0, *e2 = 0;
  int nr, nr1, nc1, nr2, nc2, i, id;

  double u, m, sn, cn, dn;

  MDR *x1, *x2;
  MDR *w;
  gsl_set_error_handler_off ();
  if (nargs == 0)
  {
    printf
      ("EllipticJacobi: Jacobian elliptic functions [sn, cn, dn](u|m). Format:\n");
    printf ("EllipticJacobi:   EllipticJacobi(u,m)\n");
    rerror ("No parameters given!");
  }
  //
  // get e1
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) != MATRIX_DENSE_REAL)
    rerror ("EllipticJacobi: argument 'u' has to be MATRIX-DENSE-REAL");
  x1 = class_matrix_real (e1);
  nr1 = MNR (x1);
  nc1 = MNC (x1);
  nr = nr1;
  if (nr1 * nc1 == 0)
    rerror ("EllipticJacobi: argument 'u' has to be MATRIX-DENSE-REAL");
  if (nc1 != 1)
    rerror ("EllipticJacobi: argument 'u' has to be a single column");
  //
  // get e2
  //
  e2 = bltin_get_ent (args[1]);
  if (ent_type (e2) != MATRIX_DENSE_REAL)
    rerror ("EllipticJacobi: argument 'm' has to be MATRIX-DENSE-REAL");
  x2 = class_matrix_real (e2);
  nr2 = MNR (x2);
  nc2 = MNC (x2);
  if (nr2 * nc2 == 0)
    rerror ("EllipticJacobi: argument 'm' has to be MATRIX-DENSE-REAL");
  if (nc2 != 1)
    rerror ("EllipticJacobi: argument 'm' has to be a single column");
  if (nr1 < nr2)
    nr = nr2;
  w = mdr_Create (nr, 3);
  // calculate
  for (i = 1; i <= nr; i++)
  {
    id = i;
    if (id > nr1)
      id = nr1;
    u = Mdr1 (x1, id, 1);
    id = i;
    if (id > nr2)
      id = nr2;
    m = Mdr1 (x2, id, 1);
    gsl_sf_elljac_e (u, m, &sn, &cn, &dn);
    Mdr1 (w, i, 1) = sn;
    Mdr1 (w, i, 2) = cn;
    Mdr1 (w, i, 3) = dn;
  }
  ent_Clean (e1);
  ent_Clean (e2);

  return ent_Assign_Rlab_MDR(w);
}


//
// Erf
//
Ent *
ent_gsl_sf_erf (int nargs, Datum args[])
{
  Ent *e1 = 0;
  int nr, nc, i;
  MDR *x1, *w;
  gsl_set_error_handler_off ();
  if (nargs == 0)
  {
    printf ("Erf: Error function. Format:\n");
    printf ("Erf:   Erf(x)\n");
    rerror ("No parameters given!");
  }
  //
  // get x
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) != MATRIX_DENSE_REAL)
    rerror ("Erf: argument 'x' has to be MATRIX-DENSE-REAL");
  x1 = class_matrix_real (e1);
  nr = MNR (x1);
  nc = MNC (x1);
  if (nr * nc == 0)
    rerror ("Erf: argument 'x' has to be MATRIX-DENSE-REAL");
  w = mdr_Create (nr, nc);
  //
  // calculate
  //
for (i = 0; i < nr*nc; i++)      MdrV0 (w, i) = gsl_sf_erf (mdrV0 (x1, i));
  ent_Clean (e1);

  return ent_Assign_Rlab_MDR(w);
}

//
// Erfc
//
Ent *
ent_gsl_sf_erfc (int nargs, Datum args[])
{
  Ent *e1 = 0;
  int nr, nc, i;
  MDR *x1, *w;
  gsl_set_error_handler_off ();
  if (nargs == 0)
  {
    printf ("Erfc: Complementary error function. Format:\n");
    printf ("Erfc:   Erfc(x)\n");
    rerror ("No parameters given!");
  }
  //
  // get x
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) != MATRIX_DENSE_REAL)
    rerror ("Erfc: argument 'x' has to be MATRIX-DENSE-REAL");
  x1 = class_matrix_real (e1);
  nr = MNR (x1);
  nc = MNC (x1);
  if (nr * nc == 0)
    rerror ("Erfc: argument 'x' has to be MATRIX-DENSE-REAL");
  w = mdr_Create (nr, nc);
  //
  // calculate
  //
for (i = 0; i < nr*nc; i++)      MdrV0 (w, i) = gsl_sf_erfc (mdrV0 (x1, i));
  ent_Clean (e1);

  return ent_Assign_Rlab_MDR(w);
}

//
// ErfZ
//
Ent *
ent_gsl_sf_erfz (int nargs, Datum args[])
{
  Ent *e1 = 0;
  int nr, nc, i;
  MDR *x1, *w;
  gsl_set_error_handler_off ();
  if (nargs == 0)
  {
    printf ("ErfZ: Gaussian probability. Format:\n");
    printf ("ErfZ:   ErfZ(x)\n");
    rerror ("No parameters given!");
  }
  //
  // get x
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) != MATRIX_DENSE_REAL)
    rerror ("ErfZ: argument 'x' has to be MATRIX-DENSE-REAL");
  x1 = class_matrix_real (e1);
  nr = MNR (x1);
  nc = MNC (x1);
  if (nr * nc == 0)
    rerror ("ErfZ: argument 'x' has to be MATRIX-DENSE-REAL");
  w = mdr_Create (nr, nc);
  //
  // calculate
  //
for (i = 0; i < nr*nc; i++)      MdrV0 (w, i) = gsl_sf_erf_Z (mdrV0 (x1, i));
  ent_Clean (e1);

  return ent_Assign_Rlab_MDR(w);
}

//
// ErfQ
//
Ent *
ent_gsl_sf_erfq (int nargs, Datum args[])
{
  Ent *e1 = 0;
  int nr, nc, i;
  MDR *x1, *w;
  gsl_set_error_handler_off ();
  if (nargs == 0)
  {
    printf ("ErfQ: Upper tail of Gaussian probability function. Format:\n");
    printf ("ErfQ:   ErfQ(x)\n");
    rerror ("No parameters given!");
  }
  //
  // get x
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) != MATRIX_DENSE_REAL)
    rerror ("ErfQ: argument 'x' has to be MATRIX-DENSE-REAL");
  x1 = class_matrix_real (e1);
  nr = MNR (x1);
  nc = MNC (x1);
  if (nr * nc == 0)
    rerror ("ErfQ: argument 'x' has to be MATRIX-DENSE-REAL");
  w = mdr_Create (nr, nc);
  //
  // calculate
  //
for (i = 0; i < nr*nc; i++)      MdrV0 (w, i) = gsl_sf_erf_Q (mdrV0 (x1, i));
 ent_Clean (e1);

 return ent_Assign_Rlab_MDR(w);
}

//
// Exponential Integral
//
Ent *
ent_gsl_sf_expint (int nargs, Datum args[])
{
  Ent *e1 = 0, *e2 = 0;
  int nr, nc, i;

  int n;
  MDR *x1=0, *w=0;

  gsl_set_error_handler_off ();

  if (nargs == 0)
  {
    printf ("ExpIntegralE: The Exponential Integral of order n=1,2. Format:\n");
    printf ("ExpIntegralE:   ExpIntegralE(x,n)\n");
    rerror ("No parameters given!");
  }
  //
  // get x
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) != MATRIX_DENSE_REAL)
    rerror ("ExpIntegralE: argument 'x' has to be MATRIX-DENSE-REAL");
  x1 = class_matrix_real (e1);
  nr = MNR (x1);
  nc = MNC (x1);
  if (nr * nc == 0)
    rerror ("ExpIntegralE: argument 'x' has to be MATRIX-DENSE-REAL");
  //
  // get n
  //
  e2 = bltin_get_ent (args[1]);
  if (ent_type (e2) != MATRIX_DENSE_REAL)
    rerror ("ExpIntegralE: argument 'n' has to an integer 1,2");
  n = class_double (e2);
  w = mdr_Create (nr, nc);
  //
  // calculate
  //
for (i = 0; i < nr*nc; i++)    {
      MdrV0 (w, i) = GSL_NAN;
      if (n == 1)
        MdrV0 (w, i) = gsl_sf_expint_E1 (mdrV0 (x1, i));
      if (n == 2)
        MdrV0 (w, i) = gsl_sf_expint_E2 (mdrV0 (x1, i));
    }
  ent_Clean (e1);
  ent_Clean (e2);

  return ent_Assign_Rlab_MDR(w);
}

//
// ExpIntegralEi
//
Ent *
ent_gsl_sf_expint_ei (int nargs, Datum args[])
{
  Ent *e1 = 0;
  int nr, nc, i;
  MDR *x1, *w;
  gsl_set_error_handler_off ();
  if (nargs == 0)
  {
    printf ("ExpIntegralEi: Exponential Integral Ei. Format:\n");
    printf ("ExpIntegralEi:   ExpIntegralEi(x)\n");
    rerror ("No parameters given!");
  }
  //
  // get x
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) != MATRIX_DENSE_REAL)
    rerror ("ExpIntegralEi: argument 'x' has to be MATRIX-DENSE-REAL");
  x1 = class_matrix_real (e1);
  nr = MNR (x1);
  nc = MNC (x1);
  if (nr * nc == 0)
    rerror ("ExpIntegralEi: argument 'x' has to be MATRIX-DENSE-REAL");
  w = mdr_Create (nr, nc);
  //
  // calculate
  //
for (i = 0; i < nr*nc; i++)      MdrV0 (w, i) = gsl_sf_expint_Ei (mdrV0 (x1, i));
  ent_Clean (e1);

  return ent_Assign_Rlab_MDR(w);
}

//
// SinhIntegral
//
Ent *
ent_gsl_sf_shi (int nargs, Datum args[])
{
  Ent *e1 = 0;
  int nr, nc, i;
  MDR *x1, *w;

  gsl_set_error_handler_off ();

  if (nargs == 0)
  {
    printf ("SinhIntegral: Sinh Integral. Format:\n");
    printf ("SinhIntegral:   SinhIntegral(x)\n");
    rerror ("No parameters given!");
  }

  //
  // get x
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) != MATRIX_DENSE_REAL)
    rerror ("SinhIntegral: argument 'x' has to be MATRIX-DENSE-REAL");
  x1 = class_matrix_real (e1);
  nr = MNR (x1);
  nc = MNC (x1);
  if (nr * nc == 0)
    rerror ("SinhIntegral: argument 'x' has to be MATRIX-DENSE-REAL");
  w = mdr_Create (nr, nc);

  //
  // calculate
  //
  for (i = 0; i < nr*nc; i++)
    MdrV0 (w, i) = gsl_sf_Shi (mdrV0 (x1, i));

  ent_Clean (e1);

  return ent_Assign_Rlab_MDR(w);
}

//
// CoshIntegral
//
Ent *
ent_gsl_sf_chi (int nargs, Datum args[])
{
  Ent *e1 = 0;
  int nr, nc, i;
  MDR *x1, *w;
  gsl_set_error_handler_off ();
  if (nargs == 0)
  {
    printf ("CoshIntegral: Cosh Integral. Format:\n");
    printf ("CoshIntegral:   CoshIntegral(x)\n");
    rerror ("No parameters given!");
  }
  //
  // get x
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) != MATRIX_DENSE_REAL)
    rerror ("CoshIntegral: argument 'x' has to be MATRIX-DENSE-REAL");
  x1 = class_matrix_real (e1);
  nr = MNR (x1);
  nc = MNC (x1);
  if (nr * nc == 0)
    rerror ("CoshIntegral: argument 'x' has to be MATRIX-DENSE-REAL");
  w = mdr_Create (nr, nc);
  //
  // calculate
  //
  for (i = 0; i < nr*nc; i++)
    MdrV0 (w, i) = gsl_sf_Chi (mdrV0 (x1, i));

  ent_Clean (e1);

  return ent_Assign_Rlab_MDR(w);
}

//
// SinIntegral
//
Ent *
ent_gsl_sf_si (int nargs, Datum args[])
{
  Ent *e1 = 0;
  int nr, nc, i;
  MDR *x1, *w;

  gsl_set_error_handler_off ();

  if (nargs == 0)
  {
    printf ("SinIntegral: Sin Integral. Format:\n");
    printf ("SinIntegral:   SinIntegral(x)\n");
    rerror ("No parameters given!");
  }

  //
  // get x
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) != MATRIX_DENSE_REAL)
    rerror ("SinIntegral: argument 'x' has to be MATRIX-DENSE-REAL");
  x1 = class_matrix_real (e1);
  nr = MNR (x1);
  nc = MNC (x1);
  if (nr * nc == 0)
    rerror ("SinIntegral: argument 'x' has to be MATRIX-DENSE-REAL");
  w = mdr_Create (nr, nc);

  //
  // calculate
  //
  for (i = 0; i < nr*nc; i++)
    MdrV0 (w, i) = gsl_sf_Si (mdrV0 (x1, i));

  ent_Clean (e1);

  return ent_Assign_Rlab_MDR(w);
}

//
// CosIntegral
//
Ent *
ent_gsl_sf_ci (int nargs, Datum args[])
{
  Ent *e1 = 0;
  int nr, nc, i;
  MDR *x1, *w;
  gsl_set_error_handler_off ();
  if (nargs == 0)
  {
    printf ("CosIntegral: Cos Integral. Format:\n");
    printf ("CosIntegral:   CosIntegral(x)\n");
    rerror ("No parameters given!");
  }
  //
  // get x
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) != MATRIX_DENSE_REAL)
    rerror ("CosIntegral: argument 'x' has to be MATRIX-DENSE-REAL");
  x1 = class_matrix_real (e1);
  nr = MNR (x1);
  nc = MNC (x1);
  if (nr * nc == 0)
    rerror ("CosIntegral: argument 'x' has to be MATRIX-DENSE-REAL");
  w = mdr_Create (nr, nc);
  //
  // calculate
  //
  for (i = 0; i < nr*nc; i++)
    MdrV0 (w, i) = gsl_sf_Ci (mdrV0 (x1, i));

  ent_Clean (e1);

  return ent_Assign_Rlab_MDR(w);
}

//
// AtanIntegral
//
Ent *
ent_gsl_sf_atanint (int nargs, Datum args[])
{
  Ent *e1 = 0;
  int nr, nc, i;
  MDR *x1, *w;
  gsl_set_error_handler_off ();
  if (nargs == 0)
  {
    printf ("AtanIntegral: Arctan Integral. Format:\n");
    printf ("AtanIntegral:   AtanIntegral(x)\n");
    rerror ("No parameters given!");
  }
  //
  // get x
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) != MATRIX_DENSE_REAL)
    rerror ("AtanIntegral: argument 'x' has to be MATRIX-DENSE-REAL");
  x1 = class_matrix_real (e1);
  nr = MNR (x1);
  nc = MNC (x1);
  if (nr * nc == 0)
    rerror ("AtanIntegral: argument 'x' has to be MATRIX-DENSE-REAL");
  w = mdr_Create (nr, nc);
  //
  // calculate
  //
  for (i = 0; i < nr*nc; i++)
    MdrV0 (w, i) = gsl_sf_atanint (mdrV0 (x1, i));

  ent_Clean (e1);

  return ent_Assign_Rlab_MDR(w);
}

//
// FermiDiracIntegral
//
Ent *
ent_gsl_sf_fermidiracint (int nargs, Datum args[])
{
  Ent *e1 = 0, *e2 = 0;
  int nr, nc, i;

  int n;
  MDR *x1=0, *w;

  gsl_set_error_handler_off ();

  if (nargs == 0)
  {
    printf ("FermiDiracIntegral: The Complete Fermi-Dirac integral of order\n");
    printf ("FermiDiracIntegral: n= -1/2, 0, 1/2, 1, 3/2. Format:\n");
    printf ("FermiDiracIntegral:   FermiDiracIntegral(x,n)\n");
    rerror ("No parameters given!");
  }
  //
  // get x
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) != MATRIX_DENSE_REAL)
    rerror ("FermiDiracIntegral: argument 'x' has to be MATRIX-DENSE-REAL");
  x1 = class_matrix_real (e1);
  nr = MNR (x1);
  nc = MNC (x1);
  if (nr * nc == 0)
    rerror ("FermiDiracIntegral: argument 'x' has to be MATRIX-DENSE-REAL");
  //
  // get n
  //
  e2 = bltin_get_ent (args[1]);
  if (ent_type (e2) != MATRIX_DENSE_REAL)
    rerror
      ("FermiDiracIntegral: argument 'n' has to be -1/2, 0, 1/2, 1 or 3/2");
  n = class_double (e2);
  w = mdr_Create (nr, nc);

  //
  // calculate
  //
  switch (2*n)
  {
    case -2:
      for (i=0; i<nr*nc; i++)
        MdrV0 (w, i) = gsl_sf_fermi_dirac_m1 (mdrV0 (x1, i));
      break;

    case -1:
      for (i=0; i<nr*nc; i++)
        MdrV0 (w, i) = gsl_sf_fermi_dirac_mhalf (mdrV0 (x1, i));
      break;

    case 0:
      for (i=0; i<nr*nc; i++)
        MdrV0 (w, i) = gsl_sf_fermi_dirac_0 (mdrV0 (x1, i));
      break;

    case 1:
      for (i=0; i<nr*nc; i++)
        MdrV0 (w, i) = gsl_sf_fermi_dirac_half (mdrV0 (x1, i));
      break;

    case 2:
      for (i=0; i<nr*nc; i++)
        MdrV0 (w, i) = gsl_sf_fermi_dirac_1 (mdrV0 (x1, i));
      break;

    case 3:
      for (i=0; i<nr*nc; i++)
        MdrV0 (w, i) = gsl_sf_fermi_dirac_3half (mdrV0 (x1, i));

    default:
      break;
  }

  ent_Clean (e1);
  ent_Clean (e2);

  return ent_Assign_Rlab_MDR(w);
}


//
// LogGamma
//
Ent *
ent_gsl_sf_loggamma (int nargs, Datum args[])
{
  Ent *e1 = 0;
  int nr, nc, i;
  gsl_set_error_handler_off ();
  if (nargs == 0)
  {
    printf ("LogGamma: Analytic logarithm of Gamma function. Format:\n");
    printf ("LogGamma:   LogGamma(x)\n");
    rerror ("No parameters given!");
  }
  //
  // get x
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) != MATRIX_DENSE_REAL
      && ent_type (e1) != MATRIX_DENSE_COMPLEX)
    rerror ("LogGamma: argument 'x' has to be real or complex matrix");

  if (ent_type (e1) == MATRIX_DENSE_REAL)
  {
    //
    // real argument
    //
    MDR *x1 = class_matrix_real (e1);
    nr = MNR (x1);
    nc = MNC (x1);
    if (nr * nc == 0)
      rerror ("LogGamma: 'x' not defined.");
    MDR *w = mdr_Create (nr, nc);

    //
    // calculate
    //
    for (i = 0; i < nr*nc; i++)
      MdrV0 (w, i) = gsl_sf_lngamma (mdrV0 (x1, i));

    ent_Clean (e1);

    return ent_Assign_Rlab_MDR(w);
  }

  //
  // complex argument
  //
  MDR *x1 = ent_data (e1);
  nr = MNR (x1);
  nc = MNC (x1);
  if (nr * nc == 0)
    rerror ("LogGamma: 'x' not defined.");
  MDC *w = mdc_Create (nr, nc);

  //
  // calculate
  //
  gsl_sf_result zlnr, zarg;
  for (i = 0; i < nr*nc; i++)
  {
    gsl_sf_lngamma_complex_e (MdcV0r (x1, i), MdcV0i (x1, i), &zlnr,
                              &zarg);
    MdcV0r (w, i) = zlnr.val * cos (zarg.val);
    MdcV0i (w, i) = zlnr.val * sin (zarg.val);
  }

  ent_Clean (e1);

  return ent_Assign_Rlab_MDC(w);
}

//
// Gamma
//
Ent *
ent_gsl_sf_gamma (int nargs, Datum args[])
{
  Ent *e1 = 0;
  int nr, nc, i;
  MDR *x1, *w;
  gsl_set_error_handler_off ();
  if (nargs == 0)
  {
    printf ("Gamma: Gamma function. Format:\n");
    printf ("Gamma:   Gamma(x)\n");
    rerror ("No parameters given!");
  }
  //
  // get x
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) != MATRIX_DENSE_REAL)
    rerror ("Gamma: argument 'x' has to be MATRIX-DENSE-REAL");
  x1 = class_matrix_real (e1);
  nr = MNR (x1);
  nc = MNC (x1);
  if (nr * nc == 0)
    rerror ("Gamma: argument 'x' has to be MATRIX-DENSE-REAL");
  w = mdr_Create (nr, nc);
  //
  // calculate
  //
  for (i = 0; i < nr*nc; i++)
    MdrV0 (w, i) = gsl_sf_gamma (mdrV0 (x1, i));

  ent_Clean (e1);

  return ent_Assign_Rlab_MDR(w);
}

//
// GammaRegularized
//
Ent *
ent_gsl_sf_gammareg (int nargs, Datum args[])
{
  Ent *e1 = 0, *e2 = 0;
  int nr, nc, nr1, nc1, nr2, nc2, i, j, id, jd;
  double a, x;
  MDR *x1, *x2;
  MDR *w;
  gsl_set_error_handler_off ();
  if (nargs == 0)
  {
    printf
      ("GammaRegularized: Regularized (normalized) incomplete Gamma function.\n");
    printf ("GammaRegularized: Format:\n");
    printf ("GammaRegularized:   GammaRegularized(a,x)\n");
    rerror ("No parameters given!");
  }
  //
  // get e1
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) != MATRIX_DENSE_REAL)
    rerror ("GammaRegularized: argument 'a' has to be MATRIX-DENSE-REAL");
  x1 = class_matrix_real (e1);
  nr1 = MNR (x1);
  nc1 = MNC (x1);
  nr = nr1;
  nc = nc1;
  if (nr1 * nc1 == 0)
    rerror ("GammaRegularized: argument 'a' has to be MATRIX-DENSE-REAL");
  // get e2
  e2 = bltin_get_ent (args[1]);
  if (ent_type (e2) != MATRIX_DENSE_REAL)
    rerror ("GammaRegularized: argument 'x' has to be MATRIX-DENSE-REAL");
  x2 = class_matrix_real (e2);
  nr2 = MNR (x2);
  nc2 = MNC (x2);
  if (nr2 * nc2 == 0)
    rerror ("GammaRegularized: argument 'x' has to be MATRIX-DENSE-REAL");
  // figure out the size of output matrix
  if (nc < nc2)
    nc = nc2;
  if (nr < nr2)
    nr = nr2;
  w = mdr_Create (nr, nc);
  // calculate
  // calculate
  for (i = 1; i <= nr; i++)
    for (j = 1; j <= nc; j++)
    {
      id = i;
      jd = j;
      if (id > nr1)
        id = nr1;
      if (jd > nc1)
        jd = nc1;
      a = Mdr1 (x1, id, jd);
      id = i;
      jd = j;
      if (id > nr2)
        id = nr2;
      if (jd > nc2)
        jd = nc2;
      x = Mdr1 (x2, id, jd);
      Mdr1 (w, i, j) = gsl_sf_gamma_inc_Q (a, x);
    }
    ent_Clean (e1);
  ent_Clean (e2);

  return ent_Assign_Rlab_MDR(w);
}

//
// GammaRegularizedC
//
Ent *
ent_gsl_sf_gammaregc (int nargs, Datum args[])
{
  Ent *e1 = 0, *e2 = 0;
  int nr, nc, nr1, nc1, nr2, nc2, i, j, id, jd;
  double a, x;
  MDR *x1, *x2;
  MDR *w;
  gsl_set_error_handler_off ();
  if (nargs == 0)
  {
    printf
      ("GammaRegularizedC: Complement of regularized (normalized) incomplete\n");
    printf ("GammaRegularizedC: Gamma function.\n");
    printf ("GammaRegularizedC: Format:\n");
    printf ("GammaRegularizedC:   GammaRegularizedC(a,x)\n");
    rerror ("No parameters given!");
  }
  //
  // get e1
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) != MATRIX_DENSE_REAL)
    rerror ("GammaRegularizedC: argument 'a' has to be MATRIX-DENSE-REAL");
  x1 = class_matrix_real (e1);
  nr1 = MNR (x1);
  nc1 = MNC (x1);
  nr = nr1;
  nc = nc1;
  if (nr1 * nc1 == 0)
    rerror ("GammaRegularizedC: argument 'a' has to be MATRIX-DENSE-REAL");
  // get e2
  e2 = bltin_get_ent (args[1]);
  if (ent_type (e2) != MATRIX_DENSE_REAL)
    rerror ("GammaRegularizedC: argument 'x' has to be MATRIX-DENSE-REAL");
  x2 = class_matrix_real (e2);
  nr2 = MNR (x2);
  nc2 = MNC (x2);
  if (nr2 * nc2 == 0)
    rerror ("GammaRegularizedC: argument 'x' has to be MATRIX-DENSE-REAL");
  // figure out the size of output matrix
  if (nc < nc2)
    nc = nc2;
  if (nr < nr2)
    nr = nr2;
  w = mdr_Create (nr, nc);

  // calculate
  for (i = 1; i <= nr; i++)
    for (j = 1; j <= nc; j++)
    {
      id = i;
      jd = j;
      if (id > nr1)
        id = nr1;
      if (jd > nc1)
        jd = nc1;
      a = Mdr1 (x1, id, jd);
      id = i;
      jd = j;
      if (id > nr2)
        id = nr2;
      if (jd > nc2)
        jd = nc2;
      x = Mdr1 (x2, id, jd);
      Mdr1 (w, i, j) = gsl_sf_gamma_inc_P (a, x);
    }

  ent_Clean (e1);
  ent_Clean (e2);

  return ent_Assign_Rlab_MDR(w);
}


//
// Beta
//
Ent *
ent_gsl_sf_beta (int nargs, Datum args[])
{
  Ent *e1 = 0, *e2 = 0;
  int nr, nc, nr1, nc1, nr2, nc2, i, j, id, jd;
  double a, b;
  MDR *x1, *x2;
  MDR *w;
  gsl_set_error_handler_off ();
  if (nargs == 0)
  {
    printf ("Beta: Beta function.\n");
    printf ("Beta: Format:\n");
    printf ("Beta:   Beta(a,b)\n");
    rerror ("No parameters given!");
  }
  //
  // get e1
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) != MATRIX_DENSE_REAL)
    rerror ("Beta: argument 'a' has to be MATRIX-DENSE-REAL");
  x1 = class_matrix_real (e1);
  nr1 = MNR (x1);
  nc1 = MNC (x1);
  nr = nr1;
  nc = nc1;
  if (nr1 * nc1 == 0)
    rerror ("Beta: argument 'a' has to be MATRIX-DENSE-REAL");
  // get e2
  e2 = bltin_get_ent (args[1]);
  if (ent_type (e2) != MATRIX_DENSE_REAL)
    rerror ("Beta: argument 'b' has to be MATRIX-DENSE-REAL");
  x2 = class_matrix_real (e2);
  nr2 = MNR (x2);
  nc2 = MNC (x2);
  if (nr2 * nc2 == 0)
    rerror ("Beta: argument 'b' has to be MATRIX-DENSE-REAL");
  // figure out the size of output matrix
  if (nc < nc2)
    nc = nc2;
  if (nr < nr2)
    nr = nr2;
  w = mdr_Create (nr, nc);
  // calculate
  for (i = 1; i <= nr; i++)
    for (j = 1; j <= nc; j++)
    {
      id = i;
      jd = j;
      if (id > nr1)
        id = nr1;
      if (jd > nc1)
        jd = nc1;
      a = Mdr1 (x1, id, jd);
      id = i;
      jd = j;
      if (id > nr2)
        id = nr2;
      if (jd > nc2)
        jd = nc2;
      b = Mdr1 (x2, id, jd);
      Mdr1 (w, i, j) = gsl_sf_beta (a, b);
    }

    ent_Clean (e1);
    ent_Clean (e2);

    return ent_Assign_Rlab_MDR(w);

}

//
// LogBeta
//
Ent *
ent_gsl_sf_logbeta (int nargs, Datum args[])
{
  Ent *e1 = 0, *e2 = 0;
  int nr, nc, nr1, nc1, nr2, nc2, i, j, id, jd;
  double a, b;
  MDR *x1, *x2;
  MDR *w;
  gsl_set_error_handler_off ();
  if (nargs == 0)
  {
    printf ("LogBeta: Logarithm of Beta function.\n");
    printf ("LogBeta: Format:\n");
    printf ("LogBeta:   LogBeta(a,b)\n");
    rerror ("No parameters given!");
  }
  //
  // get e1
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) != MATRIX_DENSE_REAL)
    rerror ("LogBeta: argument 'a' has to be MATRIX-DENSE-REAL");
  x1 = class_matrix_real (e1);
  nr1 = MNR (x1);
  nc1 = MNC (x1);
  nr = nr1;
  nc = nc1;
  if (nr1 * nc1 == 0)
    rerror ("LogBeta: argument 'a' has to be MATRIX-DENSE-REAL");
  // get e2
  e2 = bltin_get_ent (args[1]);
  if (ent_type (e2) != MATRIX_DENSE_REAL)
    rerror ("LogBeta: argument 'b' has to be MATRIX-DENSE-REAL");
  x2 = class_matrix_real (e2);
  nr2 = MNR (x2);
  nc2 = MNC (x2);
  if (nr2 * nc2 == 0)
    rerror ("LogBeta: argument 'b' has to be MATRIX-DENSE-REAL");
  // figure out the size of output matrix
  if (nc < nc2)
    nc = nc2;
  if (nr < nr2)
    nr = nr2;
  w = mdr_Create (nr, nc);
  // calculate
  for (i = 1; i <= nr; i++)
    for (j = 1; j <= nc; j++)
    {
      id = i;
      jd = j;
      if (id > nr1)
        id = nr1;
      if (jd > nc1)
        jd = nc1;
      a = Mdr1 (x1, id, jd);
      id = i;
      jd = j;
      if (id > nr2)
        id = nr2;
      if (jd > nc2)
        jd = nc2;
      b = Mdr1 (x2, id, jd);
      Mdr1 (w, i, j) = gsl_sf_lnbeta (a, b);
    }
    ent_Clean (e1);
    ent_Clean (e2);

    return ent_Assign_Rlab_MDR(w);

}

//
// BetaRegularized
//
Ent *
ent_gsl_sf_betareg (int nargs, Datum args[])
{
  Ent *e1 = 0, *e2 = 0, *e3 = 0;
  int nr, nc, nr1, nc1, nr2, nc2, nr3, nc3, i, j, id, jd;
  double a, b, x;
  MDR *x1, *x2, *x3;
  MDR *w;
  gsl_set_error_handler_off ();
  if (nargs == 0)
  {
    printf
      ("BetaRegularized: Regularized (normalized) incomplete Beta function. Format:\n");
    printf ("BetaRegularized:   BetaRegularized(a,b,x)\n");
    rerror ("No parameters given!");
  }
  //
  // get e1
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) != MATRIX_DENSE_REAL)
    rerror ("BetaRegularized: argument 'a' has to be MATRIX-DENSE-REAL");
  x1 = class_matrix_real (e1);
  nr1 = MNR (x1);
  nc1 = MNC (x1);
  nr = nr1;
  nc = nc1;
  if (nr1 * nc1 == 0)
    rerror ("BetaRegularized: argument 'a' has to be MATRIX-DENSE-REAL");
  // get e2
  e2 = bltin_get_ent (args[1]);
  if (ent_type (e2) != MATRIX_DENSE_REAL)
    rerror ("BetaRegularized: argument 'b' has to be MATRIX-DENSE-REAL");
  x2 = class_matrix_real (e2);
  nr2 = MNR (x2);
  nc2 = MNC (x2);
  if (nr2 * nc2 == 0)
    rerror ("BetaRegularized: argument 'b' has to be MATRIX-DENSE-REAL");
  // figure out the size of output matrix
  if (nc < nc2)
    nc = nc2;
  if (nr < nr2)
    nr = nr2;
  // get e3
  e3 = bltin_get_ent (args[2]);
  if (ent_type (e3) != MATRIX_DENSE_REAL)
    rerror ("BetaRegularized: argument 'x' has to be MATRIX-DENSE-REAL");
  x3 = class_matrix_real (e3);
  nr3 = MNR (x3);
  nc3 = MNC (x3);
  if (nr3 * nc3 == 0)
    rerror ("BetaRegularized: argument 'x' has to be MATRIX-DENSE-REAL");
  // figure out the size of output matrix
  if (nc < nc3)
    nc = nc3;
  if (nr < nr3)
    nr = nr3;
  w = mdr_Create (nr, nc);
  // calculate
  // calculate
  for (i = 1; i <= nr; i++)
    for (j = 1; j <= nc; j++)
    {
      id = i;
      jd = j;
      if (id > nr1)
        id = nr1;
      if (jd > nc1)
        jd = nc1;
      a = Mdr1 (x1, id, jd);
      id = i;
      jd = j;
      if (id > nr2)
        id = nr2;
      if (jd > nc2)
        jd = nc2;
      b = Mdr1 (x2, id, jd);
      id = i;
      jd = j;
      if (id > nr3)
        id = nr3;
      if (jd > nc3)
        jd = nc3;
      x = Mdr1 (x3, id, jd);
      Mdr1 (w, i, j) = gsl_sf_beta_inc (a, b, x);
    }
    ent_Clean (e1);
  ent_Clean (e2);
  ent_Clean (e3);

  return ent_Assign_Rlab_MDR(w);
}

//
// RecGamma
//
Ent *
ent_gsl_sf_recgamma (int nargs, Datum args[])
{
  Ent *e1 = 0;
  int nr, nc, i;
  MDR *x1, *w;
  gsl_set_error_handler_off ();
  if (nargs == 0)
  {
    printf ("RecGamma: Reciprocal of Gamma function. Format:\n");
    printf ("RecGamma:   RecGamma(x)\n");
    rerror ("No parameters given!");
  }
  //
  // get x
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) != MATRIX_DENSE_REAL)
    rerror ("RecGamma: argument 'x' has to be MATRIX-DENSE-REAL");
  x1 = class_matrix_real (e1);
  nr = MNR (x1);
  nc = MNC (x1);
  if (nr * nc == 0)
    rerror ("RecGamma: argument 'x' has to be MATRIX-DENSE-REAL");
  w = mdr_Create (nr, nc);
  //
  // calculate
  //
  for (i = 0; i < nr*nc; i++)
    MdrV0 (w, i) = gsl_sf_gammainv (mdrV0 (x1, i));

  ent_Clean (e1);

  return ent_Assign_Rlab_MDR(w);
}

//
// Pochhammer
//
Ent *
ent_gsl_sf_poch (int nargs, Datum args[])
{
  Ent *e1 = 0, *e2 = 0;
  int nr, nc, nr1, nc1, nr2, nc2, i, j, id, jd;
  double a, x;
  MDR *x1, *x2;
  MDR *w;
  gsl_set_error_handler_off ();
  if (nargs == 0)
  {
    printf ("Pochhammer: Pochhammer (Apell) symbox (a)_x. Format:\n");
    printf ("Pochhammer:   Pochhammer(a,x)\n");
    rerror ("No parameters given!");
  }
  //
  // get e1
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) != MATRIX_DENSE_REAL)
    rerror ("Pochhammer: argument 'a' has to be MATRIX-DENSE-REAL");
  x1 = class_matrix_real (e1);
  nr1 = MNR (x1);
  nc1 = MNC (x1);
  nr = nr1;
  nc = nc1;
  if (nr1 * nc1 == 0)
    rerror ("Pochhammer: argument 'a' has to be MATRIX-DENSE-REAL");
  // get e2
  e2 = bltin_get_ent (args[1]);
  if (ent_type (e2) != MATRIX_DENSE_REAL)
    rerror ("Pochhammer: argument 'x' has to be MATRIX-DENSE-REAL");
  x2 = class_matrix_real (e2);
  nr2 = MNR (x2);
  nc2 = MNC (x2);
  if (nr2 * nc2 == 0)
    rerror ("Pochhammer: argument 'x' has to be MATRIX-DENSE-REAL");
  // figure out the size of output matrix
  if (nc < nc2)
    nc = nc2;
  if (nr < nr2)
    nr = nr2;
  w = mdr_Create (nr, nc);
  // calculate
  for (i = 1; i <= nr; i++)
    for (j = 1; j <= nc; j++)
    {
      id = i;
      jd = j;
      if (id > nr1)
        id = nr1;
      if (jd > nc1)
        jd = nc1;
      a = Mdr1 (x1, id, jd);
      id = i;
      jd = j;
      if (id > nr2)
        id = nr2;
      if (jd > nc2)
        jd = nc2;
      x = Mdr1 (x2, id, jd);
      Mdr1 (w, i, j) = gsl_sf_poch (a, x);
    }
    ent_Clean (e1);
  ent_Clean (e2);

  return ent_Assign_Rlab_MDR(w);
}


//
// GegenbauerC
//
Ent *
ent_gsl_sf_gegen (int nargs, Datum args[])
{
  Ent *e1 = 0, *e2 = 0, *e3 = 0;
  int nr, nc, nr1, nc1, nr2, nc2, nr3, nc3, i, j, id, jd;
  double lambda, x;
  int n;
  MDR *x1, *x2, *x3, *w;
  gsl_set_error_handler_off ();
  if (nargs == 0)
  {
    printf
      ("GegenbauerC: Gegenbauer (Ultraspherical) polynomials of order 'n'.\n");
    printf ("GegenbauerC: Format:\n");
    printf ("GegenbauerC:   GegenbauerC(n,lambda,x)\n");
    rerror ("No parameters given!");
  }
  //
  // get n
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) != MATRIX_DENSE_REAL)
    rerror ("GegenbauerC: argument 'n' has to be MATRIX-DENSE-REAL");
  x1 = class_matrix_real (e1);
  nr1 = MNR (x1);
  nc1 = MNC (x1);
  nr = nr1;
  nc = nc1;
  if (nr1 * nc1 == 0)
    rerror ("GegenbauerC: argument 'n' has to be MATRIX-DENSE-REAL");
  //
  // get lambda
  //
  e2 = bltin_get_ent (args[1]);
  if (ent_type (e2) != MATRIX_DENSE_REAL)
    rerror ("GegenbauerC: argument 'lambda' has to be MATRIX-DENSE-REAL");
  x2 = class_matrix_real (e2);
  nr2 = MNR (x2);
  nc2 = MNC (x2);
  if (nr2 * nc2 == 0)
    rerror ("GegenbauerC: argument 'lambda' has to be MATRIX-DENSE-REAL");
  if (nr < nr2)
    nr = nr2;
  if (nc < nc2)
    nc = nc2;
  //
  // get x
  //
  e3 = bltin_get_ent (args[2]);
  if (ent_type (e3) != MATRIX_DENSE_REAL)
    rerror ("GegenbauerC: argument 'x' has to be MATRIX-DENSE-REAL");
  x3 = class_matrix_real (e3);
  nr3 = MNR (x3);
  nc3 = MNC (x3);
  if (nr3 * nc3 == 0)
    rerror ("GegenbauerC: argument 'x' has to be MATRIX-DENSE-REAL");
  if (nr < nr3)
    nr = nr3;
  if (nc < nc3)
    nc = nc3;
  w = mdr_Create (nr, nc);
  //
  //
  // calculate
  //
  for (i = 1; i <= nr; i++)
    for (j = 1; j <= nc; j++)
    {
      id = i;
      jd = j;
      if (id > nr1)
        id = nr1;
      if (jd > nc1)
        jd = nc1;
      n = Mdr1 (x1, id, jd);
      id = i;
      jd = j;
      if (id > nr2)
        id = nr2;
      if (jd > nc2)
        jd = nc2;
      lambda = Mdr1 (x2, id, jd);
      id = i;
      jd = j;
      if (id > nr3)
        id = nr3;
      if (jd > nc3)
        jd = nc3;
      x = Mdr1 (x3, id, jd);
      if (n == 1)
      {
        Mdr1 (w, i, j) = gsl_sf_gegenpoly_1 (lambda, x);
        continue;
      }
      if (n == 2)
      {
        Mdr1 (w, i, j) = gsl_sf_gegenpoly_2 (lambda, x);
        continue;
      }
      if (n == 3)
      {
        Mdr1 (w, i, j) = gsl_sf_gegenpoly_3 (lambda, x);
        continue;
      }
      Mdr1 (w, i, j) = gsl_sf_gegenpoly_n (n, lambda, x);
    }

  ent_Clean (e1);
  ent_Clean (e2);
  ent_Clean (e3);

  return ent_Assign_Rlab_MDR(w);
}


//
// Hypergeometric0F1
//
Ent *
ent_gsl_sf_hyperg_0F1 (int nargs, Datum args[])
{
  Ent *e1 = 0, *e2 = 0;
  int nr, nc, nr1, nc1, nr2, nc2, i, j, id, jd;
  double c, x;
  MDR *x1, *x2;
  MDR *w;
  gsl_set_error_handler_off ();
  if (nargs == 0)
  {
    printf ("Hypergeometric0F1: Hypergeometric funtion 0F1. Format:\n");
    printf ("Hypergeometric0F1:   Hypergeometric0F1(c,x)\n");
    rerror ("No parameters given!");
  }
  //
  // get e1
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) != MATRIX_DENSE_REAL)
    rerror ("Hypergeometric0F1: argument 'c' has to be MATRIX-DENSE-REAL");
  x1 = class_matrix_real (e1);
  nr1 = MNR (x1);
  nc1 = MNC (x1);
  nr = nr1;
  nc = nc1;
  if (nr1 * nc1 == 0)
    rerror ("Hypergeometric0F1: argument 'c' has to be MATRIX-DENSE-REAL");
  // get e2
  e2 = bltin_get_ent (args[1]);
  if (ent_type (e2) != MATRIX_DENSE_REAL)
    rerror ("Hypergeometric0F1: argument 'x' has to be MATRIX-DENSE-REAL");
  x2 = class_matrix_real (e2);
  nr2 = MNR (x2);
  nc2 = MNC (x2);
  if (nr2 * nc2 == 0)
    rerror ("Hypergeometric0F1: argument 'x' has to be MATRIX-DENSE-REAL");
  // figure out the size of output matrix
  if (nc < nc2)
    nc = nc2;
  if (nr < nr2)
    nr = nr2;
  w = mdr_Create (nr, nc);
  // calculate
  for (i = 1; i <= nr; i++)
    for (j = 1; j <= nc; j++)
    {
      id = i;
      jd = j;
      if (id > nr1)
        id = nr1;
      if (jd > nc1)
        jd = nc1;
      c = Mdr1 (x1, id, jd);
      id = i;
      jd = j;
      if (id > nr2)
        id = nr2;
      if (jd > nc2)
        jd = nc2;
      x = Mdr1 (x2, id, jd);
      Mdr1 (w, i, j) = gsl_sf_hyperg_0F1 (c, x);
    }
    ent_Clean (e1);
  ent_Clean (e2);

  return ent_Assign_Rlab_MDR(w);
}


//
// Hypergeometric1F1
//
Ent *
ent_gsl_sf_hyperg_1F1 (int nargs, Datum args[])
{
  Ent *e1 = 0, *e2 = 0, *e3 = 0;
  int nr, nc, nr1, nc1, nr2, nc2, nr3, nc3, i, j, id, jd;
  double a, b, x;
  MDR *x1, *x2, *x3;
  MDR *w;
  gsl_set_error_handler_off ();
  if (nargs == 0)
  {
    printf ("Hypergeometric1F1: Hypergeometric function 1F1. Format:\n");
    printf ("Hypergeometric1F1:   Hypergeometric1F1(a,b,x)\n");
    rerror ("No parameters given!");
  }
  //
  // get e1
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) != MATRIX_DENSE_REAL)
    rerror ("Hypergeometric1F1: argument 'a' has to be MATRIX-DENSE-REAL");
  x1 = class_matrix_real (e1);
  nr1 = MNR (x1);
  nc1 = MNC (x1);
  nr = nr1;
  nc = nc1;
  if (nr1 * nc1 == 0)
    rerror ("Hypergeometric1F1: argument 'a' has to be MATRIX-DENSE-REAL");
  // get e2
  e2 = bltin_get_ent (args[1]);
  if (ent_type (e2) != MATRIX_DENSE_REAL)
    rerror ("Hypergeometric1F1: argument 'b' has to be MATRIX-DENSE-REAL");
  x2 = class_matrix_real (e2);
  nr2 = MNR (x2);
  nc2 = MNC (x2);
  if (nr2 * nc2 == 0)
    rerror ("Hypergeometric1F1: argument 'b' has to be MATRIX-DENSE-REAL");
  // figure out the size of output matrix
  if (nc < nc2)
    nc = nc2;
  if (nr < nr2)
    nr = nr2;
  // get e3
  e3 = bltin_get_ent (args[2]);
  if (ent_type (e3) != MATRIX_DENSE_REAL)
    rerror ("Hypergeometric1F1: argument 'x' has to be MATRIX-DENSE-REAL");
  x3 = class_matrix_real (e3);
  nr3 = MNR (x3);
  nc3 = MNC (x3);
  if (nr3 * nc3 == 0)
    rerror ("Hypergeometric1F1: argument 'x' has to be MATRIX-DENSE-REAL");
  // figure out the size of output matrix
  if (nc < nc3)
    nc = nc3;
  if (nr < nr3)
    nr = nr3;
  w = mdr_Create (nr, nc);
  // calculate
  for (i = 1; i <= nr; i++)
    for (j = 1; j <= nc; j++)
    {
      id = i;
      jd = j;
      if (id > nr1)
        id = nr1;
      if (jd > nc1)
        jd = nc1;
      a = Mdr1 (x1, id, jd);
      id = i;
      jd = j;
      if (id > nr2)
        id = nr2;
      if (jd > nc2)
        jd = nc2;
      b = Mdr1 (x2, id, jd);
      id = i;
      jd = j;
      if (id > nr3)
        id = nr3;
      if (jd > nc3)
        jd = nc3;
      x = Mdr1 (x3, id, jd);
      if (a == (int) a && b == (int) b)
        Mdr1 (w, i, j) = gsl_sf_hyperg_1F1_int ((int) a, (int) b, x);
      else
        Mdr1 (w, i, j) = gsl_sf_hyperg_1F1 (a, b, x);
    }
    ent_Clean (e1);
  ent_Clean (e2);
  ent_Clean (e3);

  return ent_Assign_Rlab_MDR(w);
}

//
// HypergeometricU
//
Ent *
ent_gsl_sf_hyperg_U (int nargs, Datum args[])
{
  Ent *e1 = 0, *e2 = 0, *e3 = 0;
  int nr, nc, nr1, nc1, nr2, nc2, nr3, nc3, i, j, id, jd;
  double a, b, x;
  MDR *x1, *x2, *x3;
  MDR *w;
  gsl_set_error_handler_off ();
  if (nargs == 0)
  {
    printf ("HypergeometricU: Confluent hypergeometric function U. Format:\n");
    printf ("HypergeometricU:   HypergeometricU(a,b,x)\n");
    rerror ("No parameters given!");
  }
  //
  // get e1
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) != MATRIX_DENSE_REAL)
    rerror ("HypergeometricU: argument 'a' has to be MATRIX-DENSE-REAL");
  x1 = class_matrix_real (e1);
  nr1 = MNR (x1);
  nc1 = MNC (x1);
  nr = nr1;
  nc = nc1;
  if (nr1 * nc1 == 0)
    rerror ("HypergeometricU: argument 'a' has to be MATRIX-DENSE-REAL");
  // get e2
  e2 = bltin_get_ent (args[1]);
  if (ent_type (e2) != MATRIX_DENSE_REAL)
    rerror ("HypergeometricU: argument 'b' has to be MATRIX-DENSE-REAL");
  x2 = class_matrix_real (e2);
  nr2 = MNR (x2);
  nc2 = MNC (x2);
  if (nr2 * nc2 == 0)
    rerror ("HypergeometricU: argument 'b' has to be MATRIX-DENSE-REAL");
  // figure out the size of output matrix
  if (nc < nc2)
    nc = nc2;
  if (nr < nr2)
    nr = nr2;
  // get e3
  e3 = bltin_get_ent (args[2]);
  if (ent_type (e3) != MATRIX_DENSE_REAL)
    rerror ("HypergeometricU: argument 'x' has to be MATRIX-DENSE-REAL");
  x3 = class_matrix_real (e3);
  nr3 = MNR (x3);
  nc3 = MNC (x3);
  if (nr3 * nc3 == 0)
    rerror ("HypergeometricU: argument 'x' has to be MATRIX-DENSE-REAL");
  // figure out the size of output matrix
  if (nc < nc3)
    nc = nc3;
  if (nr < nr3)
    nr = nr3;
  w = mdr_Create (nr, nc);
  // calculate
  for (i = 1; i <= nr; i++)
    for (j = 1; j <= nc; j++)
    {
      id = i;
      jd = j;
      if (id > nr1)
        id = nr1;
      if (jd > nc1)
        jd = nc1;
      a = Mdr1 (x1, id, jd);
      id = i;
      jd = j;
      if (id > nr2)
        id = nr2;
      if (jd > nc2)
        jd = nc2;
      b = Mdr1 (x2, id, jd);
      id = i;
      jd = j;
      if (id > nr3)
        id = nr3;
      if (jd > nc3)
        jd = nc3;
      x = Mdr1 (x3, id, jd);
      if (a == (int) a && b == (int) b)
        Mdr1 (w, i, j) = gsl_sf_hyperg_U_int ((int) a, (int) b, x);
      else
        Mdr1 (w, i, j) = gsl_sf_hyperg_U (a, b, x);
    }

  ent_Clean (e1);
  ent_Clean (e2);
  ent_Clean (e3);

  return ent_Assign_Rlab_MDR(w);
}

//
// Hypergeometric2F0
//
Ent *
ent_gsl_sf_hyperg_2F0 (int nargs, Datum args[])
{
  Ent *e1 = 0, *e2 = 0, *e3 = 0;
  int nr, nc, nr1, nc1, nr2, nc2, nr3, nc3, i, j, id, jd;
  double a, b, x;
  MDR *x1, *x2, *x3;
  MDR *w;
  gsl_set_error_handler_off ();
  if (nargs == 0)
  {
    printf ("Hypergeometric2F0: Hypergeometric function 2F0. Format:\n");
    printf ("Hypergeometric2F0:   Hypergeometric2F0(a,b,x)\n");
    rerror ("No parameters given!");
  }
  //
  // get e1
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) != MATRIX_DENSE_REAL)
    rerror ("Hypergeometric2F0: argument 'a' has to be MATRIX-DENSE-REAL");
  x1 = class_matrix_real (e1);
  nr1 = MNR (x1);
  nc1 = MNC (x1);
  nr = nr1;
  nc = nc1;
  if (nr1 * nc1 == 0)
    rerror ("Hypergeometric2F0: argument 'a' has to be MATRIX-DENSE-REAL");
  // get e2
  e2 = bltin_get_ent (args[1]);
  if (ent_type (e2) != MATRIX_DENSE_REAL)
    rerror ("Hypergeometric2F0: argument 'b' has to be MATRIX-DENSE-REAL");
  x2 = class_matrix_real (e2);
  nr2 = MNR (x2);
  nc2 = MNC (x2);
  if (nr2 * nc2 == 0)
    rerror ("Hypergeometric2F0: argument 'b' has to be MATRIX-DENSE-REAL");
  // figure out the size of output matrix
  if (nc < nc2)
    nc = nc2;
  if (nr < nr2)
    nr = nr2;
  // get e3
  e3 = bltin_get_ent (args[2]);
  if (ent_type (e3) != MATRIX_DENSE_REAL)
    rerror ("Hypergeometric2F0: argument 'x' has to be MATRIX-DENSE-REAL");
  x3 = class_matrix_real (e3);
  nr3 = MNR (x3);
  nc3 = MNC (x3);
  if (nr3 * nc3 == 0)
    rerror ("Hypergeometric2F0: argument 'x' has to be MATRIX-DENSE-REAL");
  // figure out the size of output matrix
  if (nc < nc3)
    nc = nc3;
  if (nr < nr3)
    nr = nr3;
  w = mdr_Create (nr, nc);
  // calculate
  // calculate
  for (i = 1; i <= nr; i++)
    for (j = 1; j <= nc; j++)
    {
      id = i;
      jd = j;
      if (id > nr1)
        id = nr1;
      if (jd > nc1)
        jd = nc1;
      a = Mdr1 (x1, id, jd);
      id = i;
      jd = j;
      if (id > nr2)
        id = nr2;
      if (jd > nc2)
        jd = nc2;
      b = Mdr1 (x2, id, jd);
      id = i;
      jd = j;
      if (id > nr3)
        id = nr3;
      if (jd > nc3)
        jd = nc3;
      x = Mdr1 (x3, id, jd);
      Mdr1 (w, i, j) = gsl_sf_hyperg_2F0 (a, b, x);
    }
    ent_Clean (e1);
  ent_Clean (e2);
  ent_Clean (e3);

  return ent_Assign_Rlab_MDR(w);
}

//
// Hypergeometric2F1
//
Ent *
ent_gsl_sf_hyperg_2F1 (int nargs, Datum args[])
{
  Ent *e1 = 0, *e2 = 0, *e3 = 0, *e4 = 0;
  int nr, nc, nr1, nc1, nr2, nc2, nr3, nc3, nr4, nc4, i, j, id, jd;
  double a, b, c, x;
  MDR *x1, *x2, *x3, *x4;
  MDR *w;
  gsl_set_error_handler_off ();
  if (nargs == 0)
  {
    printf ("Hypergeometric2F1: Gauss hypergeometric function 2F1. Format:\n");
    printf ("Hypergeometric2F1:   Hypergeometric2F1(a,b,c,x)\n");
    rerror ("No parameters given!");
  }
  //
  // get a
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) != MATRIX_DENSE_REAL)
    rerror ("Hypergeometric2F1: argument 'a' has to be MATRIX-DENSE-REAL");
  x1 = class_matrix_real (e1);
  nr1 = MNR (x1);
  nc1 = MNC (x1);
  nr = nr1;
  nc = nc1;
  if (nr1 * nc1 == 0)
    rerror ("Hypergeometric2F1: argument 'a' has to be MATRIX-DENSE-REAL");
  //
  // get l
  //
  e2 = bltin_get_ent (args[1]);
  if (ent_type (e2) != MATRIX_DENSE_REAL)
    rerror ("Hypergeometric2F1: argument 'b' has to be MATRIX-DENSE-REAL");
  x2 = class_matrix_real (e2);
  nr2 = MNR (x2);
  nc2 = MNC (x2);
  if (nr2 * nc2 == 0)
    rerror ("Hypergeometric2F1: argument 'b' has to be MATRIX-DENSE-REAL");
  // figure out the size of output matrix
  if (nc < nc2)
    nc = nc2;
  if (nr < nr2)
    nr = nr2;
  //
  // get n2
  //
  e3 = bltin_get_ent (args[2]);
  if (ent_type (e3) != MATRIX_DENSE_REAL)
    rerror ("Hypergeometric2F1: argument 'c' has to be MATRIX-DENSE-REAL");
  x3 = class_matrix_real (e3);
  nr3 = MNR (x3);
  nc3 = MNC (x3);
  if (nr3 * nc3 == 0)
    rerror ("Hypergeometric2F1: argument 'c' has to be MATRIX-DENSE-REAL");
  // figure out the size of output matrix
  if (nc < nc3)
    nc = nc3;
  if (nr < nr3)
    nr = nr3;
  //
  // get t
  //
  e4 = bltin_get_ent (args[3]);
  if (ent_type (e4) != MATRIX_DENSE_REAL)
    rerror ("Hypergeometric2F1: argument 'x' has to be MATRIX-DENSE-REAL");
  x4 = class_matrix_real (e4);
  nr4 = MNR (x4);
  nc4 = MNC (x4);
  if (nr4 * nc4 == 0)
    rerror ("Hypergeometric2F1: argument 'x' has to be MATRIX-DENSE-REAL");
  // figure out the size of output matrix
  if (nc < nc4)
    nc = nc4;
  if (nr < nr4)
    nr = nr4;
  w = mdr_Create (nr, nc);
  // calculate
  for (i = 1; i <= nr; i++)
    for (j = 1; j <= nc; j++)
    {
      id = i;
      jd = j;
      if (id > nr1)
        id = nr1;
      if (jd > nc1)
        jd = nc1;
      a = Mdr1 (x1, id, jd);
      id = i;
      jd = j;
      if (id > nr2)
        id = nr2;
      if (jd > nc2)
        jd = nc2;
      b = Mdr1 (x2, id, jd);
      id = i;
      jd = j;
      if (id > nr3)
        id = nr3;
      if (jd > nc3)
        jd = nc3;
      c = Mdr1 (x3, id, jd);
      id = i;
      jd = j;
      if (id > nr4)
        id = nr4;
      if (jd > nc4)
        jd = nc4;
      x = Mdr1 (x4, id, jd);
      Mdr1 (w, i, j) = gsl_sf_hyperg_2F1 (a, b, c, x);
      //Mdr1 (w, i, j) = hyp2f1 (a, b, c, x);
    }

  ent_Clean (e1);
  ent_Clean (e2);
  ent_Clean (e3);
  ent_Clean (e4);

  return ent_Assign_Rlab_MDR(w);
}

//
// LaguerreL
//
Ent *
ent_gsl_sf_laguerre (int nargs, Datum args[])
{
  Ent *e1 = 0, *e2 = 0, *e3 = 0;
  int nr, nc, nr1, nc1, nr2, nc2, nr3, nc3, i, j, id, jd;
  double a, x;
  int n;
  MDR *x1, *x2, *x3, *w;
  gsl_set_error_handler_off ();
  if (nargs == 0)
  {
    printf ("LaguerreL: Generalized Laguerre polynomials of order 'n'.\n");
    printf ("LaguerreL: Format:\n");
    printf ("LaguerreL:   LaguerreL(n,a,x)\n");
    rerror ("No parameters given!");
  }
  //
  // get n
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) != MATRIX_DENSE_REAL)
    rerror ("LaguerreL: argument 'n' has to be MATRIX-DENSE-REAL");
  x1 = class_matrix_real (e1);
  nr1 = MNR (x1);
  nc1 = MNC (x1);
  nr = nr1;
  nc = nc1;
  if (nr1 * nc1 == 0)
    rerror ("LaguerreL: argument 'n' has to be MATRIX-DENSE-REAL");
  //
  // get lambda
  //
  e2 = bltin_get_ent (args[1]);
  if (ent_type (e2) != MATRIX_DENSE_REAL)
    rerror ("LaguerreL: argument 'a' has to be MATRIX-DENSE-REAL");
  x2 = class_matrix_real (e2);
  nr2 = MNR (x2);
  nc2 = MNC (x2);
  if (nr2 * nc2 == 0)
    rerror ("LaguerreL: argument 'a' has to be MATRIX-DENSE-REAL");
  if (nr < nr2)
    nr = nr2;
  if (nc < nc2)
    nc = nc2;
  //
  // get x
  //
  e3 = bltin_get_ent (args[2]);
  if (ent_type (e3) != MATRIX_DENSE_REAL)
    rerror ("LaguerreL: argument 'x' has to be MATRIX-DENSE-REAL");
  x3 = class_matrix_real (e3);
  nr3 = MNR (x3);
  nc3 = MNC (x3);
  if (nr3 * nc3 == 0)
    rerror ("LaguerreL: argument 'x' has to be MATRIX-DENSE-REAL");
  if (nr < nr3)
    nr = nr3;
  if (nc < nc3)
    nc = nc3;
  w = mdr_Create (nr, nc);

  //
  // calculate
  //
  for (i = 1; i <= nr; i++)
    for (j = 1; j <= nc; j++)
    {
      id = i;
      jd = j;
      if (id > nr1)
        id = nr1;
      if (jd > nc1)
        jd = nc1;
      n = Mdr1 (x1, id, jd);
      id = i;
      jd = j;
      if (id > nr2)
        id = nr2;
      if (jd > nc2)
        jd = nc2;
      a = Mdr1 (x2, id, jd);
      id = i;
      jd = j;
      if (id > nr3)
        id = nr3;
      if (jd > nc3)
        jd = nc3;
      x = Mdr1 (x3, id, jd);
      if (n == 1)
      {
        Mdr1 (w, i, j) = gsl_sf_laguerre_1 (a, x);
        continue;
      }
      if (n == 2)
      {
        Mdr1 (w, i, j) = gsl_sf_laguerre_2 (a, x);
        continue;
      }
      if (n == 3)
      {
        Mdr1 (w, i, j) = gsl_sf_laguerre_3 (a, x);
        continue;
      }
      Mdr1 (w, i, j) = gsl_sf_laguerre_n (n, a, x);
    }

  ent_Clean (e1);
  ent_Clean (e2);
  ent_Clean (e3);

  return ent_Assign_Rlab_MDR(w);
}

//
// ProductLog
//
Ent *
ent_gsl_sf_lambert_W0 (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0;
  int nr, nc, i, j;
  MDR *x1=0, *w, *x2=0;
  gsl_set_error_handler_off ();
  if (nargs == 0)
  {
    printf ("ProductLog: Lambert's W function. Format:\n");
    printf ("ProductLog:   ProductLog(x/,m/)\n");
    printf
      ("ProductLog: where m=0 (default),-1 is the branch of the solution.\n");
    rerror ("No parameters given!");
  }
  //
  // get x
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) != MATRIX_DENSE_REAL)
    rerror ("ProductLog: argument 'x' has to be MATRIX-DENSE-REAL");
  x1 = class_matrix_real (e1);
  nr = MNR (x1);
  nc = MNC (x1);
  if (nr * nc == 0)
    rerror ("ProductLog: argument 'x' has to be MATRIX-DENSE-REAL");
  w = mdr_Create (nr, nc);
  if (nargs >= 2)
  {
    e2 = bltin_get_ent (args[1]);
    if (ent_type (e2) == MATRIX_DENSE_REAL)
      x2 = class_matrix_real (e2);
  }
  else
    x2 = mdr_CreateScalar (0.0);

  //
    // calculate
    //
    for (i = 1; i <= nr; i++)
      for (j = 1; j <= nc; j++)
      {
        if (Mdr1 (x2, MIN(i, MNR (x2)), MIN(j, MNC (x2))) == -1)
          Mdr1 (w, i, j) = gsl_sf_lambert_Wm1 (Mdr1 (x1, i, j));
        else
          Mdr1 (w, i, j) = gsl_sf_lambert_W0 (Mdr1 (x1, i, j));
      }

  //
  // clean
  //
  ent_Clean (e1);

  if (!e2)
    mdr_Destroy(x2);

  ent_Clean (e2);

  return ent_Assign_Rlab_MDR(w);
}


//
// LegendreP
//
Ent *
ent_gsl_sf_legendreP (int nargs, Datum args[])
{
  Ent *e1 = 0, *e2 = 0, *e3 = 0;
  int nr, nc, nr1, nc1, nr2, nc2, nr3, nc3, i, j, id, jd;
  double x;
  int l, m;
  MDR *x1, *x2, *x3, *w=0;

  gsl_set_error_handler_off ();
  if (nargs == 0)
  {
    printf ("LegendreP: (Generalized) Legendre polynomials P of order 'l'.\n");
    printf ("LegendreP: (1) Format:\n");
    printf ("LegendreP:   LegendreP(l,x)\n");
    printf ("LegendreP: (2) Format:\n");
    printf ("LegendreP:   LegendreP(l,m,x)\n");
    rerror ("No parameters given!");
  }
  //
  // get l
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) != MATRIX_DENSE_REAL)
    rerror ("LegendreP: argument 'l' has to be MATRIX-DENSE-REAL");
  x1 = class_matrix_real (e1);
  nr1 = MNR (x1);
  nc1 = MNC (x1);
  nr = nr1;
  nc = nc1;
  if (nr1 * nc1 == 0)
    rerror ("LegendreP: argument 'l' has to be MATRIX-DENSE-REAL");
  //
  // get m or x
  //
  e2 = bltin_get_ent (args[1]);
  if (ent_type (e2) != MATRIX_DENSE_REAL)
    rerror ("LegendreP: argument 'x' ('m') has to be MATRIX-DENSE-REAL");
  x2 = class_matrix_real (e2);
  nr2 = MNR (x2);
  nc2 = MNC (x2);
  if (nr2 * nc2 == 0)
    rerror ("LegendreP: argument 'x' ('m') has to be MATRIX-DENSE-REAL");
  if (nr < nr2)
    nr = nr2;
  if (nc < nc2)
    nc = nc2;
  if (nargs == 2)
  {
    w = mdr_Create (nr, nc);
    //
    // calculate
    //
    for (i = 1; i <= nr; i++)
      for (j = 1; j <= nc; j++)
      {
        id = i;
        jd = j;
        if (id > nr1)
          id = nr1;
        if (jd > nc1)
          jd = nc1;
        l = Mdr1 (x1, id, jd);
        id = i;
        jd = j;
        if (id > nr2)
          id = nr2;
        if (jd > nc2)
          jd = nc2;
        x = Mdr1 (x2, id, jd);
        if (l == 1)
          MdrV0 (w, i) = gsl_sf_legendre_P1 (x);
        else if (l == 2)
          MdrV0 (w, i) = gsl_sf_legendre_P2 (x);
        else if (l == 3)
          MdrV0 (w, i) = gsl_sf_legendre_P3 (x);
        else
          MdrV0 (w, i) = gsl_sf_legendre_Pl (l, x);
      }
  }
  if (nargs == 3)
  {
    //
    // get x
    //
    e3 = bltin_get_ent (args[2]);
    if (ent_type (e3) != MATRIX_DENSE_REAL)
      rerror ("LegendreP: argument 'x' has to be MATRIX-DENSE-REAL");
    x3 = class_matrix_real (e3);
    nr3 = MNR (x3);
    nc3 = MNC (x3);
    if (nr3 * nc3 == 0)
      rerror ("LegendreP: argument 'x' has to be MATRIX-DENSE-REAL");
    if (nr < nr3)
      nr = nr3;
    if (nc < nc3)
      nc = nc3;
    w = mdr_Create (nr, nc);
    //
    // calculate
    //
    for (i = 1; i <= nr; i++)
      for (j = 1; j <= nc; j++)
      {
        id = i;
        jd = j;
        if (id > nr1)
          id = nr1;
        if (jd > nc1)
          jd = nc1;
        l = Mdr1 (x1, id, jd);
        id = i;
        jd = j;
        if (id > nr2)
          id = nr2;
        if (jd > nc2)
          jd = nc2;
        m = Mdr1 (x2, id, jd);
        id = i;
        jd = j;
        if (id > nr3)
          id = nr3;
        if (jd > nc3)
          jd = nc3;
        x = Mdr1 (x3, id, jd);
        MdrV0 (w, i) = gsl_sf_legendre_Plm (l, m, x);
      }
  }
  ent_Clean (e1);
  ent_Clean (e2);
  ent_Clean (e3);

  return ent_Assign_Rlab_MDR(w);
}

//
// LegendreQ
//
Ent *
ent_gsl_sf_legendreQ (int nargs, Datum args[])
{
  Ent *e1 = 0, *e2 = 0;
  int nr, nc, nr1, nc1, nr2, nc2, i, j, id, jd;
  double x;
  int l;

  MDR *x1=0, *x2=0, *w;

  gsl_set_error_handler_off ();

  if (nargs == 0)
  {
    printf ("LegendreQ: Legendre polynomials Q of order 'l'.\n");
    printf ("LegendreQ: Format:\n");
    printf ("LegendreQ:   LegendreQ(l,x)\n");
    rerror ("No parameters given!");
  }
  //
  // get l
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) != MATRIX_DENSE_REAL)
    rerror ("LegendreQ: argument 'l' has to be MATRIX-DENSE-REAL");
  x1 = class_matrix_real (e1);
  nr1 = MNR (x1);
  nc1 = MNC (x1);
  nr = nr1;
  nc = nc1;
  if (nr1 * nc1 == 0)
    rerror ("LegendreQ: argument 'l' has to be MATRIX-DENSE-REAL");
  //
  // get x
  //
  e2 = bltin_get_ent (args[1]);
  if (ent_type (e2) != MATRIX_DENSE_REAL)
    rerror ("LegendreQ: argument 'x' has to be MATRIX-DENSE-REAL");
  x2 = class_matrix_real (e2);
  nr2 = MNR (x2);
  nc2 = MNC (x2);
  if (nr2 * nc2 == 0)
    rerror ("LegendreQ: argument 'x' has to be MATRIX-DENSE-REAL");
  if (nr < nr2)
    nr = nr2;
  if (nc < nc2)
    nc = nc2;
  w = mdr_Create (nr, nc);
  //
  // calculate
  //
  for (i = 1; i <= nr; i++)
    for (j = 1; j <= nc; j++)
    {
      id = i;
      jd = j;
      if (id > nr1)
        id = nr1;
      if (jd > nc1)
        jd = nc1;
      l = Mdr1 (x1, id, jd);
      id = i;
      jd = j;
      if (id > nr2)
        id = nr2;
      if (jd > nc2)
        jd = nc2;
      x = Mdr1 (x2, id, jd);
      if (l == 0)
      {
        Mdr1 (w, i, j) = gsl_sf_legendre_Q0 (x);
        continue;
      }
      if (l == 1)
      {
        Mdr1 (w, i, j) = gsl_sf_legendre_Q1 (x);
        continue;
      }
      Mdr1 (w, i, j) = gsl_sf_legendre_Ql (l, x);
    }

  ent_Clean (e1);
  ent_Clean (e2);

  return ent_Assign_Rlab_MDR(w);
}

//
// LegendreSphericalP
//
Ent *
ent_gsl_sf_legendre_sph_Plm (int nargs, Datum args[])
{
  Ent *e1 = 0, *e2 = 0, *e3 = 0;
  int nr, nc, nr1, nc1, nr2, nc2, nr3, nc3, i, j, id, jd;
  double x;
  int l, m;
  MDR *x1, *x2, *x3, *w;
  gsl_set_error_handler_off ();
  if (nargs == 0)
  {
    printf
      ("LegendreSphericalP: Normalized Associated Legendre polynomials P.\n");
    printf ("LegendreSphericalP: Format:\n");
    printf ("LegendreSphericalP:   LegendreSphericalP(l,m,x)\n");
    rerror ("No parameters given!");
  }
  //
  // get l
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) != MATRIX_DENSE_REAL)
    rerror ("LegendreSphericalP: argument 'l' has to be MATRIX-DENSE-REAL");
  x1 = class_matrix_real (e1);
  nr1 = MNR (x1);
  nc1 = MNC (x1);
  nr = nr1;
  nc = nc1;
  if (nr1 * nc1 == 0)
    rerror ("LegendreSphericalP: argument 'l' has to be MATRIX-DENSE-REAL");
  //
  // get m
  //
  e2 = bltin_get_ent (args[1]);
  if (ent_type (e2) != MATRIX_DENSE_REAL)
    rerror ("LegendreSphericalP: argument 'm' has to be MATRIX-DENSE-REAL");
  x2 = class_matrix_real (e2);
  nr2 = MNR (x2);
  nc2 = MNC (x2);
  if (nr2 * nc2 == 0)
    rerror ("LegendreSphericalP: argument 'm' has to be MATRIX-DENSE-REAL");
  if (nr < nr2)
    nr = nr2;
  if (nc < nc2)
    nc = nc2;
  //
  // get x
  //
  e3 = bltin_get_ent (args[2]);
  if (ent_type (e3) != MATRIX_DENSE_REAL)
    rerror ("LegendreSphericalP: argument 'x' has to be MATRIX-DENSE-REAL");
  x3 = class_matrix_real (e3);
  nr3 = MNR (x3);
  nc3 = MNC (x3);
  if (nr3 * nc3 == 0)
    rerror ("LegendreSphericalP: argument 'x' has to be MATRIX-DENSE-REAL");
  if (nr < nr3)
    nr = nr3;
  if (nc < nc3)
    nc = nc3;
  w = mdr_Create (nr, nc);
  //
  // calculate
  //
  for (i = 1; i <= nr; i++)
    for (j = 1; j <= nc; j++)
    {
      id = i;
      jd = j;
      if (id > nr1)
        id = nr1;
      if (jd > nc1)
        jd = nc1;
      l = Mdr1 (x1, id, jd);
      id = i;
      jd = j;
      if (id > nr2)
        id = nr2;
      if (jd > nc2)
        jd = nc2;
      m = Mdr1 (x2, id, jd);
      id = i;
      jd = j;
      if (id > nr3)
        id = nr3;
      if (jd > nc3)
        jd = nc3;
      x = Mdr1 (x3, id, jd);
      Mdr1 (w, i, j) = gsl_sf_legendre_sphPlm (l, m, x);
    }


  ent_Clean (e1);
  ent_Clean (e2);
  ent_Clean (e3);

  return ent_Assign_Rlab_MDR(w);}

//
// LegendreConicalP
//
Ent *
ent_gsl_sf_legendre_con_Plm (int nargs, Datum args[])
{
  Ent *e1 = 0, *e2 = 0, *e3 = 0;
  int nr, nc, nr1, nc1, nr2, nc2, nr3, nc3, i, j, id, jd;
  double x, lambda, m;
  MDR *x1, *x2, *x3, *w;
  gsl_set_error_handler_off ();
  if (nargs == 0)
  {
    printf
      ("LegendreConicalP: Regular and irregular spherical conical Legendre\n");
    printf ("LegendreConicalP: function P, with 2*m = -1,0,1,2.. . Format:\n");
    printf ("LegendreConicalP:   LegendreConicalP(m,lambda,x)\n");
    rerror ("No parameters given!");
  }
  //
  // get m
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) != MATRIX_DENSE_REAL)
    rerror ("LegendreConicalP: argument 'm' has to be MATRIX-DENSE-REAL");
  x1 = class_matrix_real (e1);
  nr1 = MNR (x1);
  nc1 = MNC (x1);
  nr = nr1;
  nc = nc1;
  if (nr1 * nc1 == 0)
    rerror ("LegendreConicalP: argument 'm' has to be MATRIX-DENSE-REAL");
  //
  // get lambda
  //
  e2 = bltin_get_ent (args[1]);
  if (ent_type (e2) != MATRIX_DENSE_REAL)
    rerror ("LegendreConicalP: argument 'lambda' has to be MATRIX-DENSE-REAL");
  x2 = class_matrix_real (e2);
  nr2 = MNR (x2);
  nc2 = MNC (x2);
  if (nr2 * nc2 == 0)
    rerror ("LegendreConicalP: argument 'lambda' has to be MATRIX-DENSE-REAL");
  if (nr < nr2)
    nr = nr2;
  if (nc < nc2)
    nc = nc2;
  //
  // get x
  //
  e3 = bltin_get_ent (args[2]);
  if (ent_type (e3) != MATRIX_DENSE_REAL)
    rerror ("LegendreConicalP: argument 'x' has to be MATRIX-DENSE-REAL");
  x3 = class_matrix_real (e3);
  nr3 = MNR (x3);
  nc3 = MNC (x3);
  if (nr3 * nc3 == 0)
    rerror ("LegendreConicalP: argument 'x' has to be MATRIX-DENSE-REAL");
  if (nr < nr3)
    nr = nr3;
  if (nc < nc3)
    nc = nc3;
  w = mdr_Create (nr, nc);

  //
  // calculate
  //
  for (i = 1; i <= nr; i++)
    for (j = 1; j <= nc; j++)
    {
      id = i;
      jd = j;
      if (id > nr1)
        id = nr1;
      if (jd > nc1)
        jd = nc1;
      m = Mdr1 (x1, id, jd);
      id = i;
      jd = j;
      if (id > nr2)
        id = nr2;
      if (jd > nc2)
        jd = nc2;
      lambda = Mdr1 (x2, id, jd);
      id = i;
      jd = j;
      if (id > nr3)
        id = nr3;
      if (jd > nc3)
        jd = nc3;
      x = Mdr1 (x3, id, jd);
      if (m == 1)
      {
        Mdr1 (w, i, j) = gsl_sf_conicalP_1 (lambda, x);
        continue;
      }
      if (m == 0.5)
      {
        Mdr1 (w, i, j) = gsl_sf_conicalP_half (lambda, x);
        continue;
      }
      if (m == 0)
      {
        Mdr1 (w, i, j) = gsl_sf_conicalP_0 (lambda, x);
        continue;
      }
      if (m == -0.5)
      {
        Mdr1 (w, i, j) = gsl_sf_conicalP_mhalf (lambda, x);
        continue;
      }
      if (m == (int) m)
      {
        Mdr1 (w, i, j) = gsl_sf_conicalP_cyl_reg (m, lambda, x);
        continue;
      }
      if (2 * m == (int) (2 * m))
      {
        Mdr1 (w, i, j) = gsl_sf_conicalP_sph_reg (m, lambda, x);
        continue;
      }
      Mdr1 (w, i, j) = GSL_NAN;
    }

  ent_Clean (e1);
  ent_Clean (e2);
  ent_Clean (e3);

  return ent_Assign_Rlab_MDR(w);
}

//
// LegendreH3d
//
Ent *
ent_gsl_sf_legendre_h3d (int nargs, Datum args[])
{
  Ent *e1 = 0, *e2 = 0, *e3 = 0;
  int nr, nc, nr1, nc1, nr2, nc2, nr3, nc3, i, j, id, jd;
  double lambda, eta;
  int l;
  MDR *x1, *x2, *x3, *w;
  gsl_set_error_handler_off ();
  if (nargs == 0)
  {
    printf
      ("LegendreH3d: Regular eigenfunctions of a 3-D Laplacian in hyperbolic.\n");
    printf ("LegendreH3d: coordinate system. Format:\n");
    printf ("LegendreH3d:   LegendreH3d(l,lambda,eta)\n");
    rerror ("No parameters given!");
  }
  //
  // get m
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) != MATRIX_DENSE_REAL)
    rerror ("LegendreH3d: argument 'l' has to be MATRIX-DENSE-REAL");
  x1 = class_matrix_real (e1);
  nr1 = MNR (x1);
  nc1 = MNC (x1);
  nr = nr1;
  nc = nc1;
  if (nr1 * nc1 == 0)
    rerror ("LegendreH3d: argument 'l' has to be MATRIX-DENSE-REAL");
  //
  // get lambda
  //
  e2 = bltin_get_ent (args[1]);
  if (ent_type (e2) != MATRIX_DENSE_REAL)
    rerror ("LegendreH3d: argument 'lambda' has to be MATRIX-DENSE-REAL");
  x2 = class_matrix_real (e2);
  nr2 = MNR (x2);
  nc2 = MNC (x2);
  if (nr2 * nc2 == 0)
    rerror ("LegendreH3d: argument 'lambda' has to be MATRIX-DENSE-REAL");
  if (nr < nr2)
    nr = nr2;
  if (nc < nc2)
    nc = nc2;
  //
  // get eta
  //
  e3 = bltin_get_ent (args[2]);
  if (ent_type (e3) != MATRIX_DENSE_REAL)
    rerror ("LegendreH3d: argument 'eta' has to be MATRIX-DENSE-REAL");
  x3 = class_matrix_real (e3);
  nr3 = MNR (x3);
  nc3 = MNC (x3);
  if (nr3 * nc3 == 0)
    rerror ("LegendreH3d: argument 'eta' has to be MATRIX-DENSE-REAL");
  if (nr < nr3)
    nr = nr3;
  if (nc < nc3)
    nc = nc3;
  w = mdr_Create (nr, nc);
  //
  // calculate
  //
  for (i = 1; i <= nr; i++)
    for (j = 1; j <= nc; j++)
    {
      id = i;
      jd = j;
      if (id > nr1)
        id = nr1;
      if (jd > nc1)
        jd = nc1;
      l = Mdr1 (x1, id, jd);
      id = i;
      jd = j;
      if (id > nr2)
        id = nr2;
      if (jd > nc2)
        jd = nc2;
      lambda = Mdr1 (x2, id, jd);
      id = i;
      jd = j;
      if (id > nr3)
        id = nr3;
      if (jd > nc3)
        jd = nc3;
      eta = Mdr1 (x3, id, jd);
      if (l == 0)
      {
        Mdr1 (w, i, j) = gsl_sf_legendre_H3d_0 (lambda, eta);
        continue;
      }
      if (l == 1)
      {
        Mdr1 (w, i, j) = gsl_sf_legendre_H3d_1 (lambda, eta);
        continue;
      }
      Mdr1 (w, i, j) = gsl_sf_legendre_H3d (l, lambda, eta);
    }
    ent_Clean (e1);
    ent_Clean (e2);
    ent_Clean (e3);

    return ent_Assign_Rlab_MDR(w);
}

//
// Digamma function
//
Ent *
ent_gsl_sf_psi (int nargs, Datum args[])
{
  Ent *e1 = 0;
  int nr, nc, i;
  double x;
  MDR *x1, *w;
  gsl_set_error_handler_off ();
  if (nargs == 0)
  {
    printf ("Digamma: Digamma or Psi function. Format:\n");
    printf ("Digamma:   Digamma(x)\n");
    rerror ("No parameters given!");
  }
  //
  // get x
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) != MATRIX_DENSE_REAL)
    rerror ("Digamma: argument 'x' has to be MATRIX-DENSE-REAL");
  x1 = class_matrix_real (e1);
  nr = MNR (x1);
  nc = MNC (x1);
  if (nr * nc == 0)
    rerror ("Digamma: argument 'x' has to be MATRIX-DENSE-REAL");
  w = mdr_Create (nr, nc);
  //
  // calculate
  //
for (i = 0; i < nr*nc; i++)    {
      x = mdrV0 (x1, i);
      if (x == (int) x)
        MdrV0 (w, i) = gsl_sf_psi_int (x);
      else
        MdrV0 (w, i) = gsl_sf_psi (x);
    }

    ent_Clean (e1);

    return ent_Assign_Rlab_MDR(w);
}

//
// Polygamma
//
Ent *
ent_gsl_sf_polygamma (int nargs, Datum args[])
{
  Ent *e1 = 0, *e2 = 0;
  int nr, nc, nr1, nc1, nr2, nc2, i, j, id, jd;
  double x;
  int m;
  MDR *x1, *x2;
  MDR *w;
  gsl_set_error_handler_off ();
  if (nargs == 0)
  {
    printf ("Polygamma: Polygamma function. Format:\n");
    printf ("Polygamma:   Polygamma(m,x)\n");
    rerror ("No parameters given!");
  }
  //
  // get e1
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) != MATRIX_DENSE_REAL)
    rerror ("Polygamma: argument 'm' has to be MATRIX-DENSE-REAL");
  x1 = class_matrix_real (e1);
  nr1 = MNR (x1);
  nc1 = MNC (x1);
  nr = nr1;
  nc = nc1;
  if (nr1 * nc1 == 0)
    rerror ("Polygamma: argument 'm' has to be MATRIX-DENSE-REAL");
  // get e2
  e2 = bltin_get_ent (args[1]);
  if (ent_type (e2) != MATRIX_DENSE_REAL)
    rerror ("Polygamma: argument 'x' has to be MATRIX-DENSE-REAL");
  x2 = class_matrix_real (e2);
  nr2 = MNR (x2);
  nc2 = MNC (x2);
  if (nr2 * nc2 == 0)
    rerror ("Polygamma: argument 'x' has to be MATRIX-DENSE-REAL");
  // figure out the size of output matrix
  if (nc < nc2)
    nc = nc2;
  if (nr < nr2)
    nr = nr2;
  w = mdr_Create (nr, nc);
  // calculate
  // calculate
  for (i = 1; i <= nr; i++)
    for (j = 1; j <= nc; j++)
    {
      id = i;
      jd = j;
      if (id > nr1)
        id = nr1;
      if (jd > nc1)
        jd = nc1;
      m = Mdr1 (x1, id, jd);
      id = i;
      jd = j;
      if (id > nr2)
        id = nr2;
      if (jd > nc2)
        jd = nc2;
      x = Mdr1 (x2, id, jd);
      Mdr1 (w, i, j) = gsl_sf_psi_n (m, x);
    }

    ent_Clean (e1);
    ent_Clean (e2);

    return ent_Assign_Rlab_MDR(w);
}


//
// TransportF
//
Ent *
ent_gsl_sf_transport (int nargs, Datum args[])
{
  Ent *e1 = 0, *e2 = 0;
  int nr, nc, i;

  int n;
  MDR *x1=0, *w;

  gsl_set_error_handler_off ();

  if (nargs == 0)
  {
    printf
      ("TransportF: The Transport function of integer order n=2..5. Format:\n");
    printf ("TransportF:   TransportF(x,n)\n");
    rerror ("No parameters given!");
  }
  //
  // get x
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) != MATRIX_DENSE_REAL)
    rerror ("TransportF: argument 'x' has to be MATRIX-DENSE-REAL");
  x1 = class_matrix_real (e1);
  nr = MNR (x1);
  nc = MNC (x1);
  if (nr * nc == 0)
    rerror ("TransportF: argument 'x' has to be MATRIX-DENSE-REAL");
  //
  // get n
  //
  e2 = bltin_get_ent (args[1]);
  if (ent_type (e2) != MATRIX_DENSE_REAL)
    rerror ("TransportF: argument 'n' has to an integer 2..5");
  n = class_double (e2);
  w = mdr_Create (nr, nc);
  //
  // calculate
  //
for (i = 0; i < nr*nc; i++)    {
      MdrV0 (w, i) = GSL_NAN;
      if (n == 2)
        MdrV0 (w, i) = gsl_sf_transport_2 (mdrV0 (x1, i));
      if (n == 3)
        MdrV0 (w, i) = gsl_sf_transport_3 (mdrV0 (x1, i));
      if (n == 4)
        MdrV0 (w, i) = gsl_sf_transport_4 (mdrV0 (x1, i));
      if (n == 5)
        MdrV0 (w, i) = gsl_sf_transport_5 (mdrV0 (x1, i));
    }
    ent_Clean (e1);
    ent_Clean (e2);

    return ent_Assign_Rlab_MDR(w);
}


//
// SynchrotronF
//
Ent *
ent_gsl_sf_synchrotron (int nargs, Datum args[])
{
  Ent *e1 = 0, *e2 = 0;
  int nr, nc, i;

  int n;

  MDR *x1=0, *w;

  gsl_set_error_handler_off ();

  if (nargs == 0)
  {
    printf
      ("SynchrotronF: The Synchrotron function of integer order n=1,2. Format:\n");
    printf ("SynchrotronF:   SynchrotronF(x,n)\n");
    rerror ("No parameters given!");
  }
  //
  // get x
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) != MATRIX_DENSE_REAL)
    rerror ("SynchrotronF: argument 'x' has to be MATRIX-DENSE-REAL");
  x1 = class_matrix_real (e1);
  nr = MNR (x1);
  nc = MNC (x1);
  if (nr * nc == 0)
    rerror ("SynchrotronF: argument 'x' has to be MATRIX-DENSE-REAL");
  //
  // get n
  //
  e2 = bltin_get_ent (args[1]);
  if (ent_type (e2) != MATRIX_DENSE_REAL)
    rerror ("SynchrotronF: argument 'n' has to an integer 2..5");
  n = class_double (e2);
  w = mdr_Create (nr, nc);
  //
  // calculate
  //
  switch (n)
  {
    case 1:
      for (i = 0; i < nr*nc; i++)
        MdrV0 (w, i) = gsl_sf_synchrotron_1 (mdrV0 (x1, i));
      break;

    case 2:
      for (i = 0; i < nr*nc; i++)
        MdrV0 (w, i) = gsl_sf_synchrotron_2 (mdrV0 (x1, i));
      break;

  }

  ent_Clean (e1);
  ent_Clean (e2);

  return ent_Assign_Rlab_MDR(w);
}

//
// Zeta and Hurwitz Zeta functions
//
Ent *
ent_gsl_sf_zeta (int nargs, Datum args[])
{
  Ent *e1 = 0, *e2 = 0;
  int nr, nc, nr1, nc1, nr2, nc2, i, j, id, jd;
  double s, q;
  MDR *x1, *x2;
  MDR *w;
  gsl_set_error_handler_off ();
  if (nargs == 0)
  {
    printf ("Zeta: Riemann and Hurwitz Zeta function for real arguments.\n");
    printf ("Zeta: (1) Format:\n");
    printf ("Zeta:   Zeta(s)\n");
    printf ("Zeta: (2) Format:\n");
    printf ("Zeta:   Zeta(s,q)\n");
    rerror ("No parameters given!");
  }
  //
  // get e1
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) != MATRIX_DENSE_REAL)
    rerror ("Zeta: argument 's' has to be MATRIX-DENSE-REAL");
  x1 = class_matrix_real (e1);
  nr1 = MNR (x1);
  nc1 = MNC (x1);
  nr = nr1;
  nc = nc1;
  if (nr1 * nc1 == 0)
    rerror ("Zeta: argument 's' has to be MATRIX-DENSE-REAL");
  if (nargs == 2)
  {
    // get e2
    e2 = bltin_get_ent (args[1]);
    if (ent_type (e2) != MATRIX_DENSE_REAL)
      rerror ("Zeta: argument 'q' has to be MATRIX-DENSE-REAL");
    x2 = class_matrix_real (e2);
    nr2 = MNR (x2);
    nc2 = MNC (x2);
    if (nr2 * nc2 == 0)
      rerror ("Zeta: argument 'q' has to be MATRIX-DENSE-REAL");
    // figure out the size of output matrix
    if (nc < nc2)
      nc = nc2;
    if (nr < nr2)
      nr = nr2;
  }
  w = mdr_Create (nr, nc);
  // calculate
  if (nargs == 1)
    for (i = 1; i <= nr; i++)
      for (j = 1; j <= nc; j++)
      {
        id = i;
        jd = j;
        if (id > nr1)
          id = nr1;
        if (jd > nc1)
          jd = nc1;
        s = Mdr1 (x1, id, jd);
        if (s == (int) s)
          MdrV0 (w, i) = gsl_sf_zeta_int ((int) s);
        else
          MdrV0 (w, i) = gsl_sf_zeta (s);
      }
  if (nargs == 2)
    for (i = 1; i <= nr; i++)
      for (j = 1; j <= nc; j++)
      {
        id = i;
        jd = j;
        if (id > nr1)
          id = nr1;
        if (jd > nc1)
          jd = nc1;
        s = Mdr1 (x1, id, jd);
        id = i;
        jd = j;
        if (id > nr2)
          id = nr2;
        if (jd > nc2)
          jd = nc2;
        q = Mdr1 (x2, id, jd);
        MdrV0 (w, i) = gsl_sf_hzeta (s, q);
      }
      ent_Clean (e1);
      ent_Clean (e2);

      return ent_Assign_Rlab_MDR(w);
}

//
// erling class of functions
//
static double double_sf_erling(double *k, double *r, double *t)
{
  double a = (*r) * (*t);
  int i;
  double rval=1;

  // a^k/(k!)
  for (i=1; i<=(*k); i++)
    rval = rval * a / (i);

  //
  rval = rval * exp(-a) * (*r);
  return rval;
}

static double double_sf_erlingn(double *k, double *r, double *s, double *t)
{
  double rval;
  double a = (*r) / (*s);
  double b = 1 + 2 * (*s) * (*t);

  switch ((int) (*k))
  {
    case 0:
      rval = pow(b,-2.5)*(a+b);
      break;

    case 1:
      rval = 0.5*(b-1)*pow(b,-4.5)*(a*a + 6 * a * b + 3 * b*b);
      break;

    case 2:
      rval = 0.125*pow(b-1,2)*pow(b,-6.5)
          * (a*a*a + 15*a*a * b + 45*a*b*b + 15 * b*b*b);
      break;

    case 3:
      rval = 1./48.*pow(b-1,3)*pow(b,-8.5)
          * (pow(a,4) + 28*pow(a,3)*b + 210*a*a*b*b + 420*a*pow(b,3) + 105*pow(b,4));
      break;

    case 4:
      rval = 1./384. * pow(b-1,4)*pow(b,-10.5)
          *(pow(a,5) + 45*pow(a,4)*b + 630*pow(a,3)*b*b + 3150*a*a*pow(b,3)
          + 4725*a*pow(b,4) + 945*pow(b,5));
      break;

    case 5:
      rval = 1./3840.*pow(b-1,5)*pow(b,-12.5)
          *(pow(a,6) + 66*pow(a,5)*b + 1485*pow(a,4)*b*b + 13860*pow(a,3)*pow(b,3)
          + 51975*a*a*pow(b,4) + 62370*a*pow(b,5) + 10395*pow(b,6));
      break;

    case 6:
      rval = 1./46080.*pow(b-1,6)*pow(b,-14.5)
          * (pow(a,7) + 91.*pow(a,6)*b + 3003.*pow(a,5)*b*b + 45045.*pow(a,4)*pow(b,3)
          + 315315.*pow(a,3)*pow(b,4) + 945945.*a*a*pow(b,5) + 945945.*a*pow(b,6)
          + 135135.*pow(b,7));

    case 7:
      rval = 1./645120.*pow(b-1,7)*pow(b,-16.5)
          *(pow(a,8) + 120.0*pow(a,7)*b + 5640.0*pow(a,6)*b*b + 120120.0*pow(a,5)*pow(b,3)
              + 1353150.0*pow(a,4)*pow(b,4) + 7567560.0*pow(a,3)*pow(b,5)
              + 18918900.0*pow(a,2)*pow(b,6) + 16216200.0*a*pow(b,7) + 2027025.0*pow(b,8));
      break;

    case 8:
      rval = 1./10321920.*pow(b-1,8)*pow(b,-18.5)
          * (pow(a,9) + 153.*pow(a,8)*b + 9180.*pow(a,7)* pow(b,2) + 278460.*pow(a,6)*pow(b,3)
              + 4594590.*pow(a,5)*pow(b,4) + 41351310.*pow(a,4)*pow(b,5)
              + 192972780.*pow(a,3)*pow(b,6) + 413513100.*pow(a,2)*pow(b,7)
              + 310134825.*a*pow(b,8) + 34459425.*pow(b,9));
      break;

    case 9:
      rval = 1./185794560.*pow(b-1,9)*pow(b,-20.5)
          * (pow(a,10) + 190.*pow(a,9)*b + 14535.*pow(a,8)*pow(b,2) + 581400.*pow(a,7)*pow(b,3)
              + 13226850.*pow(a,6)*pow(b,4) + 174594420.*pow(a,5)*pow(b,5)
              + 1309458150.*pow(a,4)*pow(b,6) + 5237832600.*pow(a,3)*pow(b,7)
              + 9820936125.*pow(a,2)*pow(b,8) + 6547290750.*a*pow(b,9) + 654729075.*pow(b,10));
          break;

    case 10:
      rval = 1./3715891200.*pow(b-1,10)*pow(b,-22.5)
          *(pow(a,11) + 231.*pow(a,10)*b + 21945.*pow(a,9)*pow(b,2) + 1119195.*pow(a,8)*pow(b,3)
              + 33575850.*pow(a,7)*pow(b,4) + 611080470.*pow(a,6)*pow(b,5)
              + 6721885170.*pow(a,5)*pow(b,6) + 43212118950.*pow(a,4)*pow(b,7)
              + 151242416325.*pow(a,3)*pow(b,8) + 252070693875.0*pow(a,2)*pow(b,9)
              + 151242416325.0*a*pow(b,10) + 13749310575.0*pow(b,11));
      break;

    default:
      rval = 2.0*pow(b,-(1.5+(*k)))*exp(-0.5*a)*(*s)*pow(b-1,(*k))*gsl_sf_gamma(1.5+(*k))
        / (M_SQRTPI*gsl_sf_gamma(1+(*k))) * gsl_sf_hyperg_1F1(1.5+(*k), 0.5, 0.5*a/b);

  } // switch ((int) (*k))

  // this is valid only for explicit formulae for ErlingN
  if((*k) < 11)
    rval *= (*s) * exp(-0.5 * a + a / (2 * b));

  //
  return rval;
}

Ent *
ent_sf_erling (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *e3=0, *e4=0;
  int nr, nc, nr1, nc1, nr2, nc2, nr3=0, nc3=0, nr4, nc4, i, j;
  double k1, r1, s1, t1;
  MDR *k=0, *r=0, *s=0, *t=0;
  MDR *w=0;
  gsl_set_error_handler_off ();
  if (nargs == 0)
  {
    printf ("ErlingN: Erling+Normal function for real arguments.\n");
    printf ("ErlingN: Format:\n");
    printf ("ErlingN: (1)  ErlingN(k,r,s,t)\n");
    printf ("ErlingN: (2)  ErlingN(k,r,t)\n");
    rerror ("No parameters given!");
  }
  //
  // get k
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) != MATRIX_DENSE_REAL)
    rerror ("ErlingN: argument 'k' has to be MATRIX-DENSE-REAL");
  k = class_matrix_real (e1);
  nr1 = MNR (k);
  nc1 = MNC (k);
  if (nr1 * nc1 == 0)
    rerror ("ErlingN: argument 'k' has to be MATRIX-DENSE-REAL");
  // figure out the size of output matrix
  nc = nc1;
  nr = nr1;

  // get r
  e2 = bltin_get_ent (args[1]);
  if (ent_type (e2) != MATRIX_DENSE_REAL)
    rerror ("ErlingN: argument 'r' has to be MATRIX-DENSE-REAL");
  r = class_matrix_real (e2);
  nr2 = MNR (r);
  nc2 = MNC (r);
  if (nr2 * nc2 == 0)
    rerror ("ErlingN: argument 'r' has to be MATRIX-DENSE-REAL");
  // figure out the size of output matrix
  nc = nc > nc2 ? nc : nc2;
  nr = nr > nr2 ? nr : nr2;

  j = 2;
  if (nargs == 4)
  {
    // get s
    e3 = bltin_get_ent (args[j]);
    if (ent_type (e2) != MATRIX_DENSE_REAL)
      rerror ("ErlingN: argument 's' has to be MATRIX-DENSE-REAL");
    s = class_matrix_real (e3);
    nr3 = MNR (s);
    nc3 = MNC (s);
    if (nr3 * nc3 == 0)
      rerror ("ErlingN: argument 's' has to be MATRIX-DENSE-REAL");
    // figure out the size of output matrix
    nc = nc > nc3 ? nc : nc3;
    nr = nr > nr3 ? nr : nr3;

    j++;
  }

  // get t
  e4 = bltin_get_ent (args[j]);
  if (ent_type (e4) != MATRIX_DENSE_REAL)
    rerror ("ErlingN: argument 't' has to be MATRIX-DENSE-REAL");
  t = class_matrix_real (e4);
  nr4 = MNR (t);
  nc4 = MNC (t);
  if (nr4 * nc4 == 0)
    rerror ("ErlingN: argument 't' has to be MATRIX-DENSE-REAL");
  // figure out the size of output matrix
  nc = nc > nc4 ? nc : nc4;
  nr = nr > nr4 ? nr : nr4;

  w = mdr_Create (nr, nc);
  // calculate
  for (i = 1; i <= nr; i++)
      for (j = 1; j <= nc; j++)
  {
    k1 = Mdr1(k, MIN(nr1,i), MIN(nc1,j) );
    r1 = Mdr1(r, MIN(nr2,i), MIN(nc2,j) );
    s1 = 0;
    t1 = Mdr1(t, MIN(nr4,i), MIN(nc4,j) );
    if (s)
    {
      // ErlingN(k,r,s,t)
      s1 = Mdr1(s, MIN(nr3,i), MIN(nc3,j) );
      if (s1>0)
        Mdr1(w, i, j) = double_sf_erlingn(&k1, &r1, &s1, &t1);
      else
        Mdr1(w, i, j) = double_sf_erling (&k1, &r1, &t1);
    }
    else
    {
      // ErlingN(k,r,t)
      Mdr1(w, i, j) = double_sf_erling (&k1, &r1, &t1);
    }
  }

  ent_Clean (e1);
  ent_Clean (e2);
  ent_Clean (e3);
  ent_Clean (e4);

  return ent_Assign_Rlab_MDR(w);
}
