// Copyright (C) 2003-2005 Marijan Kostrun
//   part of rlab+4linux project on rlabplus.sourceforge.net
//
// GSL Science Library - Levin-U transform
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
#include "symbol.h"

// gsl headers
#include <gsl/gsl_mode.h>
#include <gsl/gsl_precision.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sum.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_machine.h>

// include: standard C headers
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <math.h>

#define PI 3.1415926
#define GSL_WARNINGS_OFF


Ent *
ent_gsl_levinu (int nargs, Datum args[])
{
  Ent *X=0, *MINT=0, *MAXT=0;
  MDR *x, *sum = mdr_Create (1, 1), *err = mdr_Create (1, 1), *termsused =
    mdr_Create (1, 1);
  int min_terms = 0, max_terms = 0, nc, dummy;
  if (nargs == 0)
  {
    printf ("levinu: Levin-U transform for series summation acceleration.\n");
    printf ("levinu: Format:\n");
    printf ("levinu:   y = levinu(S/,mint,maxt/),\n");
    printf
      ("levinu: where 'S' is a row-vector with the terms of the series being\n");
    printf
      ("levinu: summed, while 'mint' and 'maxt' are the minimum and maximum\n");
    printf ("levinu: number of terms to be summed, if given.\n");
    printf
      ("levinu: Output: y = <<sum;terms;error>>, where 'sum' is the estimated\n");
    printf
      ("levinu: sum of the series 'S', 'terms' gives a number of terms needed\n");
    printf
      ("levinu: to obtain 'sum', and 'error' gives the error estimate of 'sum'.\n");
    rerror ("one or three arguments required!");
  }
  X = bltin_get_ent (args[0]);
  if (ent_type (X) != MATRIX_DENSE_REAL)
    rerror ("levinu: argument 'x' has to be a real row-vector");
  x = class_matrix_real (X);
  nc = MNC (x) * MNR (x);
  if (nargs > 1)
  {
    MINT = bltin_get_ent (args[1]);
    if (ent_type (MINT) == MATRIX_DENSE_REAL)
      min_terms = class_double (MINT);
    if (min_terms < 0)
      min_terms = 0;
  }
  if (nargs > 2)
  {
    MAXT = bltin_get_ent (args[2]);
    if (ent_type (MAXT) == MATRIX_DENSE_REAL)
      max_terms = class_double (MAXT);
    if (max_terms < 0)
      max_terms = 0;
  }
  if (max_terms < min_terms)
  {
    // users are getting creative here, assume that she made a
    // mistake and exchenge the terms to be summed
    dummy = max_terms;
    max_terms = min_terms;
    min_terms = dummy;
  }
  if (max_terms > nc)
    max_terms = nc;
  if (max_terms == min_terms)
    min_terms = 1;
  gsl_sum_levin_u_workspace *w = gsl_sum_levin_u_alloc (nc);
  if (max_terms == 0)
    gsl_sum_levin_u_accel (MDRPTR(x), nc, w, MDRPTR(sum), MDRPTR(err));
  else
  {
    //printf("min = %i, max = %i\n",min_terms, max_terms);
    gsl_sum_levin_u_minmax (MDRPTR(x), nc, min_terms, max_terms, w, MDRPTR(sum), MDRPTR(err));
  }
  Mdr1 (termsused, 1, 1) = w->terms_used;
  //
  // clean up
  //
  ent_Clean (X);
  ent_Clean (MINT);
  ent_Clean (MAXT);

  gsl_sum_levin_u_free (w);

  Btree *bw = btree_Create ();
  install (bw, "sum", ent_Assign_Rlab_MDR (sum));
  install (bw, "error", ent_Assign_Rlab_MDR (err));
  install (bw, "terms", ent_Assign_Rlab_MDR (termsused));
  return ent_Assign_Rlab_BTREE (bw);
}
