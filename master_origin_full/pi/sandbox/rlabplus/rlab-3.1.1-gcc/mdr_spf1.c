// mdr_spf1.c
// This file is a part of RLaB2 Rel.2
//  Copyright (C) 2005-2008 Marijan Kostrun. project rlabplus.

/*
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
   **********************************************************************
*/

#include "rlab.h"
#include "ent.h"
#include "class.h"
#include "symbol.h"
#include "mem.h"
#include "mdr.h"
#include "mdrf1.h"
#include "complex.h"
#include "mdc.h"
#include "mdcf1.h"
#include "mds.h"
#include "list.h"
#include "btree.h"
#include "bltin.h"
#include "util.h"
#include "mathl.h"
#include "function.h"
#include "lp.h"
#include "sp.h"

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <time.h>

// naming convention for the solver parameters
#include "rlab_solver_parameters_names.h"

//
// clibs/harminv
//
#include "blas.h"
#include "harminv.c"

#undef THIS_SOLVER
#define THIS_SOLVER "harminv"
#define NFMIN 100
#define NFMAX 300
Ent *
ent_harminv (int nargs, Datum args[])
{
  Ent *e2=0, *e1=0, *eo=0;
  MDR *y=0, *nofs=0, *F=0;
  MDC *cy = 0;
  MDR *freq=0, *decay=0, *err=0;
  MDC *amp;
  int i, j, n, nrf, idummy, nf, jmax;
  double ddummy, fmin, fmax;
  ListNode *node=0;

  //
  // harminv parameters
  //
  int    solve_mode = 0;
  double dt = 1.0;
  double density = 1;
  int    nfmin   = NFMIN;
  int    nfmax   = NFMAX;
  double err_thresh     =  0.1;
  double rel_err_thresh =  1e30;
  double amp_thresh     =  0.0;
  double rel_amp_thresh = -1.0;
  double Q_thresh       =  10.0;

  mode_ok_data ok_d;
  harminv_data hd;
  Complex *data = NULL;

  if (nargs < 2 || nargs > 3)
    rerror (THIS_SOLVER ": requires two or three arguments");

  //
  // Y(t)
  //
  e1 = bltin_get_ent (args[0]);
  if ((ent_type (e1) != MATRIX_DENSE_REAL)&& (ent_type (e1)!=MATRIX_DENSE_COMPLEX))
    rerror(THIS_SOLVER ": improper first argument (signal)");

  if (ent_type (e1) == MATRIX_DENSE_REAL)
  {
    y  = ent_data (e1);
    if ( !EQVECT(y) )
      rerror(THIS_SOLVER ": improper first argument (signal)");
    n  = SIZE(y);
  }
  else
  {
    cy = ent_data (e1);
    if ( !EQVECT(cy) )
      rerror(THIS_SOLVER ": improper first argument (signal)");
    n  = SIZE(cy);
  }

  data = (Complex*) GC_malloc(sizeof(Complex) * n);
  if (y)
  {
    for (i = 0; i < n; i++)
      data[i] = MdrV0(y,i);
  }
  else if (cy)
  {
    for (i = 0; i < n; i++)
      data[i] = MdcV0(cy,i);
  }

  //
  // F = [fmin,fmax]
  //
  e2 = bltin_get_ent (args[1]);
  if (ent_type (e2) != MATRIX_DENSE_REAL)
    rerror(THIS_SOLVER ": improper second argument (F)");
  F = ent_data(e2);
  nrf = SIZE(F);
  if ( nrf!=2 )
    rerror(THIS_SOLVER ": improper second argument (F)");

  //
  // options
  //
  if (nargs == 3)
  {
    eo = bltin_get_ent (args[2]);
    if (ent_type (eo) == BTREE)
    {
      // threshold: amplitude absolute
      node = btree_FindNode (ent_data (eo), RLAB_NAME_HARM_AMPABS);
      if (node != 0)
      {
        ddummy = class_double (var_ent (node));
        if (ddummy > 0.0)
          amp_thresh = ddummy;
      }
      // threshold: amplitude relative
      node = btree_FindNode (ent_data (eo), RLAB_NAME_HARM_AMPREL);
      if (node != 0)
      {
        ddummy = class_double (var_ent (node));
        if (ddummy > 0.0)
          rel_amp_thresh = ddummy;
      }
      // threshold: error absolute
      node = btree_FindNode (ent_data (eo), RLAB_NAME_HARM_EABS);
      if (node != 0)
      {
        ddummy = class_double (var_ent (node));
        if (ddummy > 0.0)
          err_thresh = ddummy;
      }
      // threshold: error relative
      node = btree_FindNode (ent_data (eo), RLAB_NAME_HARM_EREL);
      if (node != 0)
      {
        ddummy = class_double (var_ent (node));
        if (ddummy > 0.0)
          rel_err_thresh = ddummy;
      }
      // threshold: Q
      node = btree_FindNode (ent_data (eo), RLAB_NAME_HARM_Q);
      if (node != 0)
      {
        ddummy = class_double (var_ent (node));
        if (ddummy > 0.0)
          Q_thresh = ddummy;
      }
      // dt
      node = btree_FindNode (ent_data (eo), RLAB_NAME_HARM_DT);
      if (node != 0)
      {
        ddummy = class_double (var_ent (node));
        if (ddummy > 0.0)
          dt = ddummy;
      }
      // density
      node = btree_FindNode (ent_data (eo), RLAB_NAME_HARM_DENS);
      if (node != 0)
      {
        ddummy = class_double (var_ent (node));
        if (ddummy > 0.0)
          density = ddummy;
      }
      // nf
      node = btree_FindNode (ent_data (eo), RLAB_NAME_HARM_NF);
      if (node != 0)
      {
        if (ent_type( var_ent (node)) == MATRIX_DENSE_REAL)
        {
          nofs = ent_data( var_ent (node));
        }
        else
          nofs = 0;
      }
      // solve_mode:
      //  0 - just solve
      //  1 - solve once only
      // `2 - solve ok modes only
      node = btree_FindNode (ent_data (eo), RLAB_NAME_HARM_SOLVEM);
      if (node != 0)
      {
        idummy = (int) class_double (var_ent (node));
        if (idummy >= 0 && idummy <= 2)
          solve_mode = idummy;
      }
    }

  }

  //
  // set control structure
  //
  ok_d.only_f_inrange = 1;
  ok_d.err_thresh     = err_thresh;
  ok_d.rel_err_thresh = rel_err_thresh;
  ok_d.amp_thresh     = amp_thresh;
  ok_d.rel_amp_thresh = rel_amp_thresh;
  ok_d.Q_thresh       = Q_thresh;
  ok_d.verbose        = 0;

  //
  // process: fmin, fmax
  //
  fmin = mdrV0(F,0) * dt;
  fmax = mdrV0(F,1) * dt;
  if (fmin > fmax)
  {
    double dummy = fmin;
    fmin = fmax;
    fmax = dummy;
  }
  ok_d.fmin = fmin;
  ok_d.fmax = fmax;
  ok_d.verbose = 0;

  // process: nf
  if (nofs)
  {
    nfmin = mdrV0(nofs,0);
    if (nfmin < 2)
      nfmin = 2;
    nfmax = mdrV0(nofs,1);
  }

  nf = (fmax - fmin) * n * density;
  if (nf > nfmax)
    nf = nfmax;
  if (nf < nfmin)
    nf = nfmin;


  hd = harminv_data_create(n, data, fmin, fmax, nf);
  switch (solve_mode)
  {
    case 0:
      // just solve
      harminv_solve(hd);
      break;

    case 1:
      // solve once only
      harminv_solve_once(hd);
      break;

    default:
      // solve once
      harminv_solve_ok_modes(hd,mode_ok,&ok_d);
  }

  // initialize ok_d
  mode_ok(hd, -1, &ok_d);

  // pick up the solutions and write them down:
  jmax = harminv_get_num_freqs(hd);
  if (jmax)
  {
    freq  = mdr_Create(jmax,1);
    decay = mdr_Create(jmax,1);
    err   = mdr_Create(jmax,1);
    amp   = mdc_Create(jmax,1);

    for (j=0; j<jmax; j++)
    {
      MdrV0(freq, j) = harminv_get_freq(hd, j) / dt;
      MdrV0(decay,j) = harminv_get_decay(hd, j) / ABS(dt);
      MdcV0(amp,  j) = harminv_get_amplitude (hd, j);
      MdrV0(err,  j) = harminv_get_freq_error(hd, j);
    }
  }
  if (solve_mode != 2)
    mode_ok(hd, -2, &ok_d);

  harminv_data_destroy(hd);

  // clean stuff
  GC_free (data);

  ent_Clean (e1);
  ent_Clean (e2);
  ent_Clean (eo);

  Btree *bw = btree_Create ();
  install  (bw, "freq",   ent_Assign_Rlab_MDR(freq));
  install  (bw, "decay",  ent_Assign_Rlab_MDR(decay));
  install  (bw, "amp",    ent_Assign_Rlab_MDC(amp));
  install  (bw, "err",    ent_Assign_Rlab_MDR(err));
  return ent_Assign_Rlab_BTREE(bw);
}


Ent *
ent_gcv_init (int nargs, Datum args[])
{
  Ent *e2=0, *e1=0, *eo=0, *rent;
  MDR *x=0, *y=0, *wx=0, *wy=0, *w;
  int i, nx, ny, ione = 1, ier, idummy;
  double ddummy;
  Btree    *bw=0, *sw=0;
  ListNode *node;
  int gcv_M = 2;  // spline degree is 2*M-1: default is cubic spline

  //
  // gcv variables
  //
  int gcv_MD = 2; // MD = 1 (smoothing parameter p in gcv_VAL); 2 generalized cross validation;
                  // 3 true predicted meand mean-squared error (variance in gcv_VAL);
                  // 4 number of degrees of freedom in gcv_VAL (0 < gcv_VAL <= N-M)
  double gcv_VAL = 1; // Mode value, as described above under MD

  if (nargs != 2 && nargs != 3)
    rerror ("gcvsplinefit: requires two or three arguments");

  //
  // get Y
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) == MATRIX_DENSE_REAL)
    y = class_matrix_real (e1);
  else if (ent_type (e1) == BTREE)
  {
    // y.val
    node = btree_FindNode (ent_data (e1), "val");
    if (node != 0)
    {
      if (ent_type (var_ent (node)) == MATRIX_DENSE_REAL)
        y = class_matrix_real(var_ent (node));
      else
        rerror ("gcvsplinefit: list entry 'val' is not real");

      if(MNR(y)!=1 && MNC(y)!=1)
        rerror ("gcvsplinefit: list entry 'val' is not vector");
    }
    else
      rerror ("gcvsplinefit: list entry 'val' is missing");

    // y.wgt
    node = btree_FindNode (ent_data (e1), "wgt");
    if (node != 0)
    {
      if (ent_type (var_ent (node)) == MATRIX_DENSE_REAL)
        wy = mdr_Float_BF (ent_data(var_ent (node)));
      else
        rerror ("gcvsplinefit: list entry 'wgt' is not real");

      if (MNR(wy)*MNC(wy) != MNR(y)*MNC(y))
        rerror ("gcvsplinefit: list entries 'val' and 'wgt' mismatch");
    }
  }
  else
    rerror ("gcvsplinefit: missing  first argument");

  if(!wy)
  {
    wy = mdr_Create (MNR (y) * MNC (y), 1);
    for (i = 0; i < MNR (y) * MNC (y); i++)
      MdrV0 (wy, i) = 1.0;
  }

  //
  // get X
  //
  e2 = bltin_get_ent (args[1]);
  if (ent_type (e2) == MATRIX_DENSE_REAL)
    x = ent_data(e2);
  else if (ent_type (e2) == BTREE)
  {
    // x.val
    node = btree_FindNode (ent_data (e2), "val");
    if (node != 0)
    {
      if (ent_type (var_ent (node)) == MATRIX_DENSE_REAL)
        x = class_matrix_real(var_ent (node));
      else
        rerror ("gcvsplinefit: list entry 'val' is not real");
      if(MNR(x)!=1 && MNC(x)!=1)
        rerror ("gcvsplinefit: list entry 'val' is not vector");
    }
    else
      rerror ("gcvsplinefit: list entry 'val' is missing");
    // x.wgt
    node = btree_FindNode (ent_data (e2), "wgt");
    if (node != 0)
    {
      if (ent_type (var_ent (node)) == MATRIX_DENSE_REAL)
        wx = mdr_Float_BF (ent_data(var_ent (node)));
      else
        rerror ("gcvsplinefit: list entry 'wgt' is not real");

      if (MNR(wx)*MNC(wx) != MNR(x)*MNC(x))
        rerror ("gcvsplinefit: list entries 'val' and 'wgt' mismatch");
    }
  }
  else
    rerror ("gcvsplinefit: missing  second argument");

  if (!wx)
  {
    wx = mdr_Create (MNR (x) * MNC (x), 1);
    for (i = 0; i < MNR (x) * MNC (x); i++)
      MdrV0 (wx, i) = 1.0;
  }
  //
  // options
  //
  if (nargs > 2)
  {
    eo = bltin_get_ent (args[2]);
    if (ent_type (eo) == BTREE)
    {
      // spldeg
      node = btree_FindNode (ent_data (eo), RLAB_NAME_GCV_DEGREE);
      if (node != 0)
      {
        idummy = (int) class_double (var_ent (node));
        switch (idummy)
        {
          case 2:
            fprintf(stdout,"gcvsplinefit: requested degree=%i out of range\n", idummy);
            fprintf(stdout,"gcvsplinefit: using degree=1 instead!\n");
          case 1:
            // linear
            gcv_M = 1;
            break;

          case 4:
            fprintf(stdout,"gcvsplinefit: requested degree=%i out of range\n", idummy);
            fprintf(stdout,"gcvsplinefit: using degree=3 instead!\n");
          case 3:
            //cubic
            gcv_M = 2;
            break;

          case 6:
            fprintf(stdout,"gcvsplinefit: requested degree=%i out of range\n", idummy);
            fprintf(stdout,"gcvsplinefit: using degree=5 instead!\n");
          case 5:
            //quintic
            gcv_M = 3;
            break;

          case 7:
            //heptic
            gcv_M = 4;
            break;

          default:
            fprintf(stdout,"gcvsplinefit: requested degree=%i out of range\n", idummy);
            fprintf(stdout,"gcvsplinefit: using degree=3 instead!\n");
            //cubic
            gcv_M = 2;
        }
      }
      // smoothing parameter p
      node = btree_FindNode (ent_data (eo), RLAB_NAME_GCV_SMOOTH);
      if (node != 0)
      {
        ddummy = class_double (var_ent (node));
        if (ddummy > 0.0)
        {
          gcv_VAL = ddummy;
          gcv_MD = 1;
        }
      }
      // variance
      node = btree_FindNode (ent_data (eo), RLAB_NAME_GCV_VAR);
      if (node != 0)
      {
        ddummy = class_double (var_ent (node));
        if (ddummy > 0.0)
        {
          gcv_VAL = ddummy;
          gcv_MD = 3;
        }
      }
      // degrees of freedom
      node = btree_FindNode (ent_data (eo), RLAB_NAME_GCV_DF);
      if (node != 0)
      {
        ddummy = class_double (var_ent (node));
        if (ddummy > 0.0)
        {
          gcv_VAL = ddummy;
          gcv_MD = 4;
        }
      }
    }
    if (ent_type (eo) != UNDEF) ent_Clean(eo);
  }

  nx = MNR (x) * MNC (x);
  ny = MNR (y) * MNC (y);
  if (nx != ny)
    rerror ("gcvsplinefit: 'x' and 'y' do not match");

  MDR *c  = mdr_Create (nx, ione);
  MDR *wk = mdr_Create (nx + 6 * (nx * gcv_M + 1), 1);

  GCVSPL (MDRPTR(x), MDRPTR(y), &ny, MDRPTR(wx), MDRPTR(wy), &gcv_M,
          &nx, &ione, &gcv_MD, &gcv_VAL, MDRPTR(c), &nx, MDRPTR(wk), &ier);

  if (ier != 0)
  {
    fprintf (stdout, "gcvsplinefit: error number %i.\n", ier);
    if (ier == 1)
    {
      fprintf (stdout,
               "gcvsplinefit: spline degree < 0 or insufficient number\n");
      fprintf (stdout, "gcvsplinefit: of data points.\n");
    }
    if (ier == 2)
      fprintf (stdout,
               "gcvsplinefit: 'x' data non-increasing or negative weights.\n");
    if (ier == 3)
      fprintf (stdout, "gcvsplinefit: Wrong mode, parameter or value.\n");
    rerror ("error in parameters/datasets!");
  }

  w = mdr_Create (nx, 2);
  for (i = 0; i < nx; i++)
  {
    Mdr0 (w, i, 0) = MdrV0 (x, i);
    Mdr0 (w, i, 1) = MdrV0 (c, i);
  }
  // clean stuff
  ent_Clean (e1);
  ent_Clean (e2);

  //
  // put order of spline in 'w'
  //
  bw = btree_Create ();

  // degree
  Ent * R = ent_Create ();
  ent_data (R) = mdr_CreateScalar ( (double) (2*gcv_M - 1) );
  ent_type (R) = MATRIX_DENSE_REAL;
  install (bw, "degree", R);
  ent_Clean (R);

  // m
  Ent * S  = ent_Create ();
  ent_data (S) = w;
  ent_type (S) = MATRIX_DENSE_REAL;
  install (bw, "m", S);

  //
  // statistics information
  //
  sw = btree_Create ();
  // Generalized Cross Validation value
  Ent * T  = ent_Create ();
  ent_data (T) = mdr_CreateScalar ( MdrV0(wk,0) );
  ent_type (T) = MATRIX_DENSE_REAL;
  install (sw, "gcv_value", T);

  // Mean Squared Residual.
  Ent * O  = ent_Create ();
  ent_data (O) = mdr_CreateScalar ( MdrV0(wk,1) );
  ent_type (O) = MATRIX_DENSE_REAL;
  install (sw, "residual", O);

  // number of degrees of freedom of the residual sum of squares
  // per dataset
  Ent * P  = ent_Create ();
  ent_data (P) = mdr_CreateScalar ( MdrV0(wk,2) );
  ent_type (P) = MATRIX_DENSE_REAL;
  install (sw, "df", P);

  // Smoothing parameter p
  Ent * Q  = ent_Create ();
  ent_data (Q) = mdr_CreateScalar ( MdrV0(wk,3) );
  ent_type (Q) = MATRIX_DENSE_REAL;
  install (sw, "smooth", Q);

  // Estimate of the true mean squared error
  Ent * X  = ent_Create ();
  ent_data (X) = mdr_CreateScalar ( MdrV0(wk,4) );
  ent_type (X) = MATRIX_DENSE_REAL;
  install (sw, "mse", X);

  // Gauss-Markov error variance
  Ent * Y  = ent_Create ();
  ent_data (Y) = mdr_CreateScalar ( MdrV0(wk,5) );
  ent_type (Y) = MATRIX_DENSE_REAL;
  install (sw, "var", Y);

  //
  Ent * Z  = ent_Create ();
  ent_data (Z) = sw;
  ent_type (Z) = BTREE;
  install (bw, "info", Z);

  // finish clean-up
  mdr_Destroy (c);
  mdr_Destroy (wx);
  mdr_Destroy (wy);
  mdr_Destroy (wk);


  // return result
  rent = ent_Create ();
  ent_data (rent) = bw;
  ent_type (rent) = BTREE;
  return rent;
}

Ent *
ent_gcv_val (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *e3=0, *rent;
  MDR *x=0, *w=0, *ider=0, *q=0, *c=0;
  int i, j, nx, n, idummy, idf;
  int gcv_M = 3;

  ListNode *node;

  if (nargs < 2 || nargs > 3)
    rerror ("gcvsplineval: requires two or three arguments");

  //
  // get X
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) == MATRIX_DENSE_REAL)
    x = class_matrix_real(e1);
  else if (ent_type (e1) == BTREE)
  {
    // x.val
    node = btree_FindNode (ent_data (e1), "val");
    if (node != 0)
    {
      if (ent_type (var_ent (node)) == MATRIX_DENSE_REAL)
        x = class_matrix_real(var_ent (node));
      else
        rerror ("gcvsplineval: list entry 'val' is not real");
    }
    else
      rerror ("gcvsplineval: list entry 'val' is missing");
  }
  else
    rerror ("gcvsplineval: missing first argument");

  if(MNR(x)!=1 && MNC(x)!=1)
    rerror ("gcvsplineval: list entry 'val' is not vector");
  nx = MNR (x) * MNC (x);

  //
  // get 'm' and 'spldeg' from the second argument
  //
  e2 = bltin_get_ent (args[1]);
  if (ent_type (e2) == BTREE)
  {
    // spldeg
    node = btree_FindNode (ent_data(e2), "degree");
    if (node != 0)
    {
      idummy = (int) class_double (var_ent (node));
      switch (idummy)
      {
        case 2:
          fprintf(stdout,"gcvsplineval: degree=%i out of range\n", idummy);
          fprintf(stdout,"gcvsplineval: using degree=1 instead!\n");
        case 1:
          // linear
          gcv_M = 1;
          break;

        case 4:
          fprintf(stdout,"gcvsplineval: degree=%i out of range\n", idummy);
          fprintf(stdout,"gcvsplineval: using degree=3 instead!\n");
        case 3:
          //cubic
          gcv_M = 2;
          break;

        case 6:
          fprintf(stdout,"gcvsplineval: degree=%i out of range\n", idummy);
          fprintf(stdout,"gcvsplineval: using degree=5 instead!\n");
        case 5:
          //quintic
          gcv_M = 3;
          break;

        case 7:
          //heptic
          gcv_M = 4;
          break;

        default:
          fprintf(stdout,"gcvsplineval: degree=%i out of range\n", idummy);
          fprintf(stdout,"gcvsplineval: using degree=3 instead!\n");
          //cubic
          gcv_M = 2;
      }
    }

    else
      rerror ("gcvsplineval: missing list entry 'degree'");
    // m
    node = btree_FindNode (ent_data (e2), "m");
    if (node != 0)
      c = class_matrix_real (var_ent (node));
    else
      rerror ("gcvsplineval: missing list entry 'm'");
  }
  else
    rerror ("gcvsplineset: 'c' is not gcv-list");

  if (MNC (c) != 2)
    rerror ("gcvsplineset: list entry 'm' is improper");

  n = MNR (c);

  //
  // get IDER
  //
  if (nargs == 3)
  {
    e3 = bltin_get_ent (args[2]);
    if (ent_type (e3) == MATRIX_DENSE_REAL)
      ider = mdr_Float_BF ( ent_data(e3) );
  }
  if (!ider) ider = mdr_CreateScalar (0.0);

  w = mdr_Create (nx, MNR (ider) * MNC (ider));
  q = mdr_Create (2 * gcv_M, 1);
  idf = 0;
  for (i = 0; i < nx; i++)
  {
    for (j = 0; j < MNR (ider) * MNC (ider); j++)
    {
      idf = MdrV0 (ider, j);
      Mdr0 (w, i, j) =
          SPLDER (&idf, &gcv_M, &n, &MdrV0(x,i), MDRPTR(c), &MdrV0(c,n), &idummy, MDRPTR(q));
    }
  }

  // clean stuff
  ent_Clean (e1);
  ent_Clean (e2);
  ent_Clean (e3);
  mdr_Destroy (q);
  mdr_Destroy (ider);

  rent = ent_Create ();
  ent_data (rent) = w;
  ent_type (rent) = MATRIX_DENSE_REAL;
  return rent;
}


//
//
//
#undef  THIS_SOLVER
#define THIS_SOLVER "stl"
Ent *
ent_stl2 (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *e3=0, *e4=0;
  MDR *x1=0, *x2=0, *Y=0, *T=0, *u=0, *w=0, *v=0, *wr=0;
  double *x=0, *y=0, *e=0;

  double err_default = 1e-2;
  //
  // stl2
  //
  int stl_I1 = -1;
  int stl_I2 = 1;   // 1 without, 2 for band restriction
  int stl_I3 = 3;   // 3 continuity of line segments not required, 1 is required
  int i, K, m, IP;
  if (nargs < 2 || nargs > 4)
  {
    fprintf (stdout,
             THIS_SOLVER ": Continuous/discontinuous line segment approximation.\n");
    fprintf (stdout,
             THIS_SOLVER ": Format:\n");
    fprintf (stdout,
             THIS_SOLVER ":   m = stl(x,Y /,br , cc/ ),\n");
    fprintf (stdout,
             THIS_SOLVER ": where x is independent, and Y=[y /,tol/] is the \n");
    fprintf (stdout,
             THIS_SOLVER ": dependent coordinate and possibly its tolerance; br=0\n");
    fprintf (stdout,
             THIS_SOLVER ": enables band restriction while br=0 (default) disables it;\n");
    fprintf (stdout,
             THIS_SOLVER ": cc=1 for continuous line segment or cc=0 (default) for.\n");
    fprintf (stdout,
             THIS_SOLVER ": discontinuous.\n");
    rerror (RLAB_ERROR_TWO_ARG_REQUIRED);
  }

  //
  // get X
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) != MATRIX_DENSE_REAL)
    rerror (THIS_SOLVER ": " RLAB_ERROR_ARG1_MDR_VECTOR);
  x1 = ent_data (e1);
  if (MNC (x1) != 1 && MNR(x1)!=1)
    rerror (THIS_SOLVER ": " RLAB_ERROR_ARG1_MDR_VECTOR);
  if (x1->type == RLAB_TYPE_INT32)
    rerror (THIS_SOLVER ": " RLAB_ERROR_ARG1_MDR_VECTOR);
  m = MNC (x1) * MNR(x1);
  x = &MdrV0(x1,0);

  //
  // get Y
  //
  e2 = bltin_get_ent (args[1]);
  if (ent_type (e2) == MATRIX_DENSE_REAL)
  {
    x2 = ent_data (e2);
    if (x2->type == RLAB_TYPE_INT32)
      rerror (THIS_SOLVER ": " RLAB_ERROR_ARG2_MDR_VECTOR);
    if (MNC (x2) == 2)
    {
      if (MNR(x2)!=m)
        rerror (THIS_SOLVER ": "RLAB_ERROR_ARG1_ARG2_SAME_NUMBER_ROWS);
      y = &Mdr0(x2,0,0);
      e = &Mdr0(x2,m,0);
      stl_I1 = -1;
    }
    else if (MNC (x2)==1 || MNR (x2)==1 )
    {
      if (MNC(x2) * MNR(x2)!=m)
        rerror (THIS_SOLVER ": "RLAB_ERROR_ARG1_ARG2_SAME_LENGTH);
      y = &Mdr0(x2,0,0);
      e = &err_default;
      stl_I1 = 1;
    }
  }
  else if (ent_type (e2) == BTREE)
  {
    // mean of a weighted variable <<val;wgt>>
    ListNode *node1 = btree_FindNode (ent_data (e2), RLAB_NAME_STAT_VALUE);
    ListNode *node2 = btree_FindNode (ent_data (e2), RLAB_NAME_STAT_WEIGHT);
    if (node1 && node2)
    {
      Y = class_matrix_real (var_ent (node1));
      if (!Y)
        rerror (THIS_SOLVER ": ARG2 entries in <<"
                RLAB_NAME_STAT_VALUE ";" RLAB_NAME_STAT_WEIGHT
                ">> must have same size !");
      T = class_matrix_real (var_ent (node2));
      if (!T)
        rerror (THIS_SOLVER ": ARG2 entries in <<"
                RLAB_NAME_STAT_VALUE ";" RLAB_NAME_STAT_WEIGHT
                ">> must have same size !");
      if (MNR (Y) * MNC(Y) != MNR (T) * MNC(T))
        rerror (THIS_SOLVER ": ARG2 entries in <<"
            RLAB_NAME_STAT_VALUE ";" RLAB_NAME_STAT_WEIGHT
            ">> must have same size !");
      if (MNC(Y) * MNR(Y) != m)
        rerror (THIS_SOLVER ": "RLAB_ERROR_ARG1_ARG2_SAME_LENGTH);
    }
  }

  // convert weight to std using w = 1/s^2
  if (T)
  {
    e = (double *) GC_MALLOC (m * sizeof (double));
    for (i=0; i<m; i++)
    {
      if (MdrV0(T,i) < 1e-3)
        e[i]= 1e6;
      else if (MdrV0(T,i) > 1e3)
        e[i] = 1e-6;
      else
        e[i] = 1 / (MdrV0(T,i) * MdrV0(T,i));
    }
    stl_I1 = -1;
  }

  //
  // get br
  //
  if (nargs >= 3)
  {
    e3 = bltin_get_ent (args[2]);
    if (ent_type (e3) == MATRIX_DENSE_REAL)
      stl_I2 = (int) class_double (e3);
  }

  if (nargs == 4)
  {
    e4 = bltin_get_ent (args[3]);
    if (ent_type (e4) == MATRIX_DENSE_REAL)
    {
      if (ent_type (e4) == MATRIX_DENSE_REAL)
        stl_I3 = (int) class_double (e4);
      if (stl_I3 == 0)
        stl_I3 = 3;
      else
        stl_I3 = 1;
    }
  }

  u = mdr_Create (m, 1);
  mdr_Zero (u);
  v = mdr_Create (m, 1);
  mdr_Zero (v);
  w = mdr_Create (m, 1);
  mdr_Zero (w);
  IP = stl_I1 * stl_I2 * stl_I3;

  STL2 (x, y, e, &m, MDRPTR(u), MDRPTR(v), MDRPTR(w), &K, &IP);

  if (K)
  {
    if (stl_I3 == 3)
    {
      // continuous line segment
      wr = mdr_Create (K + 1, 2);
      for (i=0; i<K+1; i++)
      {
        Mdr0 (wr, i, 0) = Mdr0 (u, i, 0);
        Mdr0 (wr, i, 1) = Mdr0 (v, i, 0);
      }
    }
    else
    {
      // discontinuous line segment
      wr = mdr_Create (2 * K, 2);
      for (i=0; i<K; i++)
      {
        Mdr0 (wr, 2 * i    , 0) = MdrV0 (u, i);
        Mdr0 (wr, 2 * i    , 1) = MdrV0 (w, i);
        Mdr0 (wr, 2 * i + 1, 0) = MdrV0 (u, i + 1);
        Mdr0 (wr, 2 * i + 1, 1) = MdrV0 (v, i + 1);
      }
    }
  }

  if (T)
    GC_FREE (e);

  // clean stuff
  ent_Clean (e1);
  ent_Clean (e2);
  ent_Clean (e3);
  ent_Clean (e4);


  mdr_Destroy(u);
  mdr_Destroy(v);
  mdr_Destroy (w);

  // write up
  return ent_Assign_Rlab_MDR(wr);
}
