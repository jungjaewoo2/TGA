// Copyright (C) 2003-2004 Marijan Kostrun
//   part of rlab+4linux project on rlabplus.sourceforge.net
//
// GSL Science Library - random numbers and statistical functions
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
#include "mdrf1.h"
#include "mdc.h"
#include "mdcf1.h"
#include "mdcf2.h"
#include "mathl.h"
#include "rlab_solver_parameters_names.h"
#include "mathl.h"
#include "symbol.h"
#include "gsl_rng.h"
#include "rfileio.h"

#include "rlab_macros.h"

// gsl headers
#include <gsl/gsl_mode.h>
#include <gsl/gsl_precision.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_machine.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_combination.h>
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_histogram2d.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_heapsort.h>
#include "sort.h"

// include: standard C headers
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <math.h>

#define PI 3.1415926
#define GSL_WARNINGS_OFF

#undef  THIS_FILE
#define THIS_FILE "gsl_stat.c"
//
// type of statistics
//
int RLAB_GSL_STAT_BIAS = 0;          // use N-bias for statistical parameters

//
// calculate mean of a matrix when no nan,inf checking is necessary
//
void mdr_Avg(int rowdominant, MDR *x, MDR **w)
{
  int i,n,nw;

  if (rowdominant)
  {
    n  = MNR (x);
    nw = MNC (x);

    // find column-wise means
    if (! *w)
      *w = mdr_Create (1, nw);

    for (i=0; i<nw; i++)
      MdrV0 (*w, i) = gsl_stats_mean (&Mdr0(x,0,i), 1, n);
  }
  else
  {
    n  = MNC (x);
    nw = MNR (x);

    // find row-wise means
    if (! *w)
      *w = (MDR *) mdr_Create (nw,1);

    for (i=0; i<nw; i++)
      MdrV0 (*w, i) = gsl_stats_mean (&Mdr0(x,i,0), nw, n);
  }

  return;
}


//
// mean
//
#undef  THIS_SOLVER
#define THIS_SOLVER "mean"
Ent * ent_mean (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *e3=0;
  MDR *x=0, *w=0, *t=0;
  MDR *ig_idx=0, *iuse_idx=0, *factor_idx=0;

  //  int xr, xc, tr, tc;
  ListNode *node=0, *node1=0, *node2=0;
  int i, j, rowdominant=0, ignore_infs=1, flat=0;

  //
  // get X
  //
  if (nargs != 1 && nargs != 2)
    rerror (THIS_SOLVER ": " RLAB_ERROR_ONE_OR_TWO_ARG_REQUIRED);

  if (nargs==2)
  {
    // copy supplied indices to integers, and adjust to c-language numbering
    e2 = bltin_get_ent (args[1]);
    if (ent_type (e2) == BTREE)
    {
      node = btree_FindNode (ent_data (e2), RLAB_NAME_STAT_FLAT);
      if (node)
      {
        if (ent_type (var_ent (node)) == MATRIX_DENSE_REAL)
        {
          MDR *dummy = ent_data(var_ent (node));
          if (SIZE(dummy)>0)
            flat = mdiV0(dummy,0) == 0 ? 0 : 1;
        }
      }
      node = btree_FindNode (ent_data (e2), RLAB_NAME_STAT_FACTOR);
      if (node)
      {
        if (ent_type (var_ent (node)) == MATRIX_DENSE_REAL)
        {
          MDR *dummy = ent_data(var_ent (node));
          if (SIZE(dummy)>0)
            if (IsFinite(dummy))
              factor_idx = dummy;
        }
      }
      node = btree_FindNode (ent_data (e2), RLAB_NAME_STAT_IGNOREIDX);
      if (node)
      {
        if (ent_type (var_ent (node)) == MATRIX_DENSE_REAL)
        {
          MDR *dummy = ent_data(var_ent (node));
          if (SIZE(dummy)>0)
            if (IsFinite(dummy))
              ig_idx = dummy;
        }
      }
      else
      {
        node = btree_FindNode (ent_data (e2), RLAB_NAME_STAT_USEIDX);
        if (node)
        {
          if (ent_type (var_ent (node)) == MATRIX_DENSE_REAL)
          {
            MDR *dummy = ent_data(var_ent (node));
            if (SIZE(dummy)>0)
              if (IsFinite(dummy))
                iuse_idx = dummy;
          }
        }
      }
      node = btree_FindNode (ent_data (e2), RLAB_NAME_STAT_IGNOREINF);
      if (node)
      {
        if (ent_type (var_ent (node)) == MATRIX_DENSE_REAL)
        {
          MDR *dummy = ent_data(var_ent (node));
          if (SIZE(dummy)>0)
            ignore_infs = mdiV0(dummy,0) == 0 ? 0 : 1;
        }
      }
      node = btree_FindNode (ent_data (e2), RLAB_NAME_STAT_ROWDOM);
      if (node)
      {
        if (ent_type (var_ent (node)) == MATRIX_DENSE_REAL)
        {
          MDR *dummy = ent_data(var_ent (node));
          if (SIZE(dummy)>0)
            rowdominant = mdiV0(dummy,0) != 0 ? 1 : 0;
        }
      }
    }
  }

  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) == MATRIX_DENSE_REAL)
  {
    x = ent_data (e1);
    w = mdr_Mean(x, factor_idx, iuse_idx, ig_idx, rowdominant, ignore_infs, flat);
  }
  else if (ent_type (e1) == BTREE)
  {
    // mean of a weighted variable <<val;wgt>>
    node1 = btree_FindNode (ent_data (e1), RLAB_NAME_STAT_VALUE);
    node2 = btree_FindNode (ent_data (e1), RLAB_NAME_STAT_WEIGHT);
    if (node1 && node2)
    {
      x = class_matrix_real (var_ent (node1));
      t = class_matrix_real (var_ent (node2));
      if (MNR (x) != MNR (t) || MNC (x) != MNC (t))
        rerror ("mean: entries in <<"
            RLAB_NAME_STAT_VALUE ";" RLAB_NAME_STAT_WEIGHT
            ">> must have same size !");
      if (MNC(x)==1 || MNR(x)==1)
        w = mdr_CreateScalar ( gsl_stats_wmean (MDRPTR(t), 1, MDRPTR(x), 1, MNC (x) * MNR (x)) );
      else
      {
        w = mdr_Create(1,MNC(x));
        for (i=0;i<MNC(x);i++)
          MdrV0(w,i) = gsl_stats_wmean ( &MdrV0(t, i*MNR (x)), 1, &MdrV0(x, MNR(x)*i), 1, MNR (x));
      }
    }
    else if ((node1 && !node2) || (!node1 && node2))
      rerror ("mean: invalid first argument <<" RLAB_NAME_STAT_VALUE ";" RLAB_NAME_STAT_WEIGHT ">>");
    else
    {
      // mean of a gsl-type histogrammed variable <<bin;range>>
      node1 = btree_FindNode (ent_data (e1), RLAB_NAME_HIST1D_BIN);
      node2 = btree_FindNode (ent_data (e1), RLAB_NAME_HIST1D_RANGE);

      if (node1 != 0 && node2 !=0)
      {
        int nhist, nbins;
        // get the values
        MDR * bin   = class_matrix_real (var_ent (node1));
        MDR * range = class_matrix_real (var_ent (node2));
        // get the size
        nbins = MNR(range)*MNC(range) - 1;

        if (MNR(bin)==1 || MNC(bin)==1)
          nhist = 1;
        else
          nhist = MNC(bin);

        w = mdr_Create(1,nhist);

        // create a histogram
        gsl_histogram * hist = gsl_histogram_alloc(nbins);
        hist->range = MDRPTR(range);
        for (j=0; j<nhist; j++)
        {
          // create pointer to the first row in a column
          hist->bin   = &Mdr0(bin,0,j);

          MdrV0(w,j) = gsl_histogram_mean(hist);
        }
        // clean up
        hist->range = 0;
        hist->bin   = 0;
        gsl_histogram_free(hist);
      }
      else
        rerror ("mean: invalid first argument <<bin;range>>");
    }
  }
  else
    rerror ("mean: argument is a real matrix or a list <<" RLAB_NAME_STAT_VALUE ";" RLAB_NAME_STAT_WEIGHT ">>!");

  ent_Clean (e1);
  ent_Clean (e2);
  ent_Clean (e3);

  return ent_Assign_Rlab_MDR(w);
}

//
// variance
//
#undef  THIS_SOLVER
#define THIS_SOLVER "var"
Ent * ent_variance (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0;
  MDR *x, *w=0, *m=0, *fm=0, *t=0, *ig_idx=0, *iuse_idx=0, *factor_idx=0;
  ListNode *node=0, *node1=0, *node2=0;
  int nr, nc, i, j, rowdominant=0, ignore_infs=1, unbias=0, flat=0;

  //
  // get X
  //
  if (nargs != 1 && nargs != 2)
    rerror (THIS_SOLVER ": " RLAB_ERROR_ONE_OR_TWO_ARG_REQUIRED);

  if (nargs==2)
  {
    // copy supplied indices to integers, and adjust to c-language numbering
    e2 = bltin_get_ent (args[1]);
    if (ent_type (e2) == BTREE)
    {
      node = btree_FindNode (ent_data (e2), RLAB_NAME_STAT_FLAT);
      if (node)
      {
        if (ent_type (var_ent (node)) == MATRIX_DENSE_REAL)
        {
          MDR *dummy = ent_data(var_ent (node));
          if (SIZE(dummy)>0)
            flat = mdiV0(dummy,0) == 0 ? 0 : 1;
        }
      }
      node = btree_FindNode (ent_data (e2), RLAB_NAME_STAT_FACTOR);
      if (node)
      {
        if (ent_type (var_ent (node)) == MATRIX_DENSE_REAL)
        {
          MDR *dummy = ent_data(var_ent (node));
          if (SIZE(dummy)>0)
            if (IsFinite(dummy))
              factor_idx = dummy;
        }
      }
      node = btree_FindNode (ent_data (e2), RLAB_NAME_STAT_BIAS);
      if (node)
      {
        MDR *dummy = ent_data(var_ent (node));
        if (SIZE(dummy)>0)
          unbias = mdiV0(dummy,0) == 0 ? 1 : 0;
      }
      node = btree_FindNode (ent_data (e2), RLAB_NAME_STAT_MEAN);
      if (node)
      {
        MDR *dummy = ent_data(var_ent (node));
        if (SIZE(dummy)>0)
          if (IsFinite(dummy))
            m = dummy;
      }
      node = btree_FindNode (ent_data (e2), RLAB_NAME_STAT_IGNOREIDX);
      if (node)
      {
        if (ent_type (var_ent (node)) == MATRIX_DENSE_REAL)
        {
          MDR *dummy = ent_data(var_ent (node));
          if (SIZE(dummy)>0)
            if (IsFinite(dummy))
              ig_idx = dummy;
        }
      }
      else
      {
        node = btree_FindNode (ent_data (e2), RLAB_NAME_STAT_USEIDX);
        if (node)
        {
          if (ent_type (var_ent (node)) == MATRIX_DENSE_REAL)
          {
            MDR *dummy = ent_data(var_ent (node));
            if (SIZE(dummy)>0)
              if (IsFinite(dummy))
                iuse_idx = dummy;
          }
        }
      }
      node = btree_FindNode (ent_data (e2), RLAB_NAME_STAT_IGNOREINF);
      if (node)
      {
        if (ent_type (var_ent (node)) == MATRIX_DENSE_REAL)
        {
          MDR *dummy = ent_data(var_ent (node));
          if (SIZE(dummy)>0)
            ignore_infs = mdiV0(dummy,0) == 0 ? 0 : 1;
        }
      }
      node = btree_FindNode (ent_data (e2), RLAB_NAME_STAT_ROWDOM);
      if (node)
      {
        if (ent_type (var_ent (node)) == MATRIX_DENSE_REAL)
        {
          MDR *dummy = ent_data(var_ent (node));
          if (SIZE(dummy)>0)
            rowdominant = mdiV0(dummy,0) != 0 ? 1 : 0;
        }
      }
    }
  }

  //
  // process the first argument
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) == MATRIX_DENSE_REAL)
  {
    x = ent_data (e1);
    w = mdr_Var(x, factor_idx, iuse_idx, ig_idx, rowdominant, ignore_infs, unbias, m, flat);
  }
  else if (ent_type (e1) == BTREE)
  {
    // mean of a weighted variable <<val;wgt>>
    node1 = btree_FindNode (ent_data (e1), RLAB_NAME_STAT_VALUE);
    node2 = btree_FindNode (ent_data (e1), RLAB_NAME_STAT_WEIGHT);
    if (node1 && node2)
    {
      x = class_matrix_real (var_ent (node1));
      nr = MNR (x);
      nc = MNC (x);

      t = class_matrix_real (var_ent (node2));
      if (nr != MNR (t) || nc != MNC (t))
        rerror ("variance: size mismatch in list <<" RLAB_NAME_STAT_VALUE ";" RLAB_NAME_STAT_WEIGHT ">>");

      if (nr==1 || nc==1)
      {
        // mu
        if (m)
        {
          w = mdr_Create (1, (m->nrow) * (m->ncol));
          for (i=0; i < (m->nrow) * (m->ncol); i++)
          {
            MdrV0 (w, i) = gsl_stats_wvariance_m (MDRPTR(t), 1, MDRPTR(x), 1, nr*nc, MdrV0(m,i));
          }
        }
        else if (fm)
        {
          w = mdr_Create (1, (fm->nrow) * (fm->ncol));
          for (i = 0; i < (fm->nrow) * (fm->ncol); i++)
          {
            MdrV0 (w, i) = gsl_stats_wvariance_with_fixed_mean (MDRPTR(t), 1, MDRPTR(x), 1, nr*nc, MdrV0(fm,i));
          }
        }
        else
          w = mdr_CreateScalar (gsl_stats_wvariance (MDRPTR(t), 1, MDRPTR(x), 1, MNR(x) * MNC(x)));
      }
      else
      {
        w = mdr_Create (1, nc);

        // mu
        if (m)
        {
          for (i=0; i < nc; i++)
          {
            MdrV0 (w, i) = gsl_stats_wvariance_m (&MdrV0(t,i*nr), 1, &MdrV0(x,i*nr), 1, nr,
                                                  MdrV0(m,MIN(i,MNR(m)*MNC(m))));
          }
        }
        else if (fm)
        {
          for (i = 0; i < nc; i++)
          {
            MdrV0 (w, i) = gsl_stats_wvariance_with_fixed_mean (&MdrV0(t,i*nr), 1, &MdrV0(x,i*nr), 1, nr,
                                                                MdrV0(fm,MIN(i,MNR(fm)*MNC(fm))));
          }
        }
        else
          for (i=0; i < nc; i++)
          {
            MdrV0(w,i) = gsl_stats_wvariance (&MdrV0(t,i*nr), 1, &MdrV0(x,i*nr), 1, MNR(x));
          }
      }

    }

    // mean of a gsl-type histogrammed variable <<bin;range>>
    node1 = btree_FindNode (ent_data (e1), RLAB_NAME_HIST1D_BIN);
    node2 = btree_FindNode (ent_data (e1), RLAB_NAME_HIST1D_RANGE);
    if (node1 && node2)
    {
      int nhist, nbins;
      // get the values
      MDR * bin   = class_matrix_real (var_ent (node1));
      MDR * range = class_matrix_real (var_ent (node2));
      // get the size
      nbins = MNR(range)*MNC(range) - 1;

      if (MNR(bin)==1 || MNC(bin)==1)
        nhist = 1;
      else
        nhist = MNC(bin);

      w = mdr_Create(1,nhist);

      // create a histogram
      gsl_histogram * hist = gsl_histogram_alloc(nbins);
      hist->range = MDRPTR(range);
      for (j=0; j<nhist; j++)
      {
        // create pointer to the first row in a column
        hist->bin   = &Mdr0(bin,0,j);
        MdrV0(w,j)  = gsl_histogram_sigma(hist);
        MdrV0(w,j) *= MdrV0(w,j);
      }

      hist->range = 0;
      hist->bin   = 0;
      gsl_histogram_free(hist);
    }

    // if we are here and w is zero, then none of the above worked
    if (!w)
      w = mdr_Create(0,0);
  }
  else
    rerror ("variance: Need an argument real matrix!");

  ent_Clean (e1);
  ent_Clean (e2);

  return ent_Assign_Rlab_MDR(w);
}

//
// covariance
//
#undef  THIS_SOLVER
#define THIS_SOLVER "covar"
Ent * ent_covariance (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0;
  MDR *x, *w=0, *m=0, *ig_idx=0, *iuse_idx=0, *factor_idx=0;
  ListNode *node=0;
  int rowdominant=0, ignore_infs=1, unbias=0;

  //
  // get X
  //
  if (nargs != 1 && nargs != 2)
    rerror (THIS_SOLVER ": " RLAB_ERROR_ONE_OR_TWO_ARG_REQUIRED);

  if (nargs==2)
  {
    // copy supplied indices to integers, and adjust to c-language numbering
    e2 = bltin_get_ent (args[1]);
    if (ent_type (e2) == BTREE)
    {
      node = btree_FindNode (ent_data (e2), RLAB_NAME_STAT_FACTOR);
      if (node)
      {
        if (ent_type (var_ent (node)) == MATRIX_DENSE_REAL)
        {
          MDR *dummy = ent_data(var_ent (node));
          if (SIZE(dummy)>0)
            if (IsFinite(dummy))
              factor_idx = dummy;
        }
      }
      node = btree_FindNode (ent_data (e2), RLAB_NAME_STAT_BIAS);
      if (node)
      {
        MDR *dummy = ent_data(var_ent (node));
        if (SIZE(dummy)>0)
          unbias = mdiV0(dummy,0) == 0 ? 1 : 0;
      }
      node = btree_FindNode (ent_data (e2), RLAB_NAME_STAT_MEAN);
      if (node)
      {
        MDR *dummy = ent_data(var_ent (node));
        if (SIZE(dummy)>0)
          if (IsFinite(dummy))
            m = dummy;
      }
      node = btree_FindNode (ent_data (e2), RLAB_NAME_STAT_IGNOREIDX);
      if (node)
      {
        if (ent_type (var_ent (node)) == MATRIX_DENSE_REAL)
        {
          MDR *dummy = ent_data(var_ent (node));
          if (SIZE(dummy)>0)
            if (IsFinite(dummy))
              ig_idx = dummy;
        }
      }
      else
      {
        node = btree_FindNode (ent_data (e2), RLAB_NAME_STAT_USEIDX);
        if (node)
        {
          if (ent_type (var_ent (node)) == MATRIX_DENSE_REAL)
          {
            MDR *dummy = ent_data(var_ent (node));
            if (SIZE(dummy)>0)
              if (IsFinite(dummy))
                iuse_idx = dummy;
          }
        }
      }
      node = btree_FindNode (ent_data (e2), RLAB_NAME_STAT_IGNOREINF);
      if (node)
      {
        if (ent_type (var_ent (node)) == MATRIX_DENSE_REAL)
        {
          MDR *dummy = ent_data(var_ent (node));
          if (SIZE(dummy)>0)
            ignore_infs = mdiV0(dummy,0) == 0 ? 0 : 1;
        }
      }
      node = btree_FindNode (ent_data (e2), RLAB_NAME_STAT_ROWDOM);
      if (node)
      {
        if (ent_type (var_ent (node)) == MATRIX_DENSE_REAL)
        {
          MDR *dummy = ent_data(var_ent (node));
          if (SIZE(dummy)>0)
            rowdominant = mdiV0(dummy,0) != 0 ? 1 : 0;
        }
      }
    }
  }

  //
  // process the first argument
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) == MATRIX_DENSE_REAL)
  {
    x = ent_data (e1);
    w = mdr_Covar(x, factor_idx, iuse_idx, ig_idx, rowdominant, ignore_infs, unbias, m);
  }

  ent_Clean (e1);
  ent_Clean (e2);

  return ent_Assign_Rlab_MDR(w);
}


//
// (weighted) mean absolute difference with given mean
//
Ent *
ent_absdev (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0;
  MDR *x=0, *t=0, *w=0, *m=0;

  int nr, nc, tr, tc, i;

  ListNode *node;

  if (nargs >= 1)
    e1 = bltin_get_ent (args[0]);
  else
    rerror ("absdev: first argument has to be a list <<" RLAB_NAME_STAT_VALUE ";" RLAB_NAME_STAT_WEIGHT ">> or a vector!");

  // start from the end: is there 'm' mean
  if (nargs > 1)
  {
    e2 = bltin_get_ent (args[1]);
    if (ent_type (e2) == MATRIX_DENSE_REAL)
    {
      m = class_matrix_real (e2);

      // ignore empty matrix
      if (m->nrow * m->ncol == 0)
        m = 0;
    }
  }

  if (ent_type (e1) == BTREE)
  {
    // val
    node = btree_FindNode (ent_data (e1), RLAB_NAME_STAT_VALUE);
    if (node != 0)
      x = class_matrix_real (var_ent (node));
    else
      rerror ("absdev: first argument has to be a list <<" RLAB_NAME_STAT_VALUE ";" RLAB_NAME_STAT_WEIGHT ">> or a vector!");

    // wgt
    node = btree_FindNode (ent_data (e1), RLAB_NAME_STAT_WEIGHT);
    if (node != 0)
      t = class_matrix_real (var_ent (node));
    else
      rerror ("absdev: first argument has to be a list <<" RLAB_NAME_STAT_VALUE ";" RLAB_NAME_STAT_WEIGHT ">> or a vector!");

    nr = x->nrow;
    nc = x->ncol;
    tr = t->nrow;
    tc = t->ncol;
    if ((nc != 1 && nr!=1) || (tc != 1 && tc != 1))
      rerror ("absdev: " RLAB_NAME_STAT_VALUE " and " RLAB_NAME_STAT_WEIGHT " have to be real vectors!");
    if (nc*nr != tc*tr)
      rerror ("absdev: " RLAB_NAME_STAT_VALUE " and " RLAB_NAME_STAT_WEIGHT " must have same size!");

    // for weighted variable only scalar 'm' is considered
    if (m)
      w = mdr_CreateScalar ( gsl_stats_wabsdev_m ( MDRPTR(t), 1, MDRPTR(x), 1, nr, MdrV0(m,0) ) );
    else
      w = mdr_CreateScalar ( gsl_stats_wabsdev   ( MDRPTR(t), 1, MDRPTR(x), 1, nr ) );
  }
  else if (ent_type (e1) == MATRIX_DENSE_REAL)
  {
    x = class_matrix_real (e1);
    nr = x->nrow;
    nc = x->ncol;

    if (nr == 1 || nc == 1)
    {
      // column vectors
      if (m)
      {
        w = mdr_Create (1, m->nrow * m->ncol);
        for (i=0; i < m->nrow * m->ncol; i++)
          MdrV0 (w, i) = gsl_stats_absdev_m ( MDRPTR(x), 1, nr * nc, MdrV0(m, i) );
      }
      else
        w = mdr_CreateScalar (gsl_stats_absdev (MDRPTR(x), 1, nr * nc));
    }
    else
    {
      // matrix of data: assume columns of data
      w = mdr_Create (1, nc);
      for (i=1; i <= nc; i++)
      {
        // no mean given
        if (m)
          MdrV1 (w, i) =
            gsl_stats_absdev_m ( &MdrV0(x,(i - 1)*nr), 1, nr,
                                MdrV1 (m, MIN(i, MNC (m) * MNR (m))));
        else
          MdrV1 (w, i) = gsl_stats_absdev ( &MdrV0(x,(i-1)*nr), 1, nr);
      }
    }
  }
  else
    rerror ("absdev: first argument has to be a list <<" RLAB_NAME_STAT_VALUE ";" RLAB_NAME_STAT_WEIGHT ">> or a vector!");

  ent_Clean (e1);
  ent_Clean (e2);

  return ent_Assign_Rlab_MDR(w);
}

//
// (weighted) skew,  with optional fixed mean and std
//
Ent *
ent_skew (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *e3=0;
  MDR *x=0, *t=0, *w=0;
  MDR *m=0, *s=0;
  int nr, nc, tr, tc, i;

  ListNode *node;

  if (nargs >= 1)
    e1 = bltin_get_ent (args[0]);
  else
    rerror ("skew: first argument has to be a list <<" RLAB_NAME_STAT_VALUE ";" RLAB_NAME_STAT_WEIGHT ">> or a vector!");

  // did user provide 'mean'
  if (nargs > 1)
  {
    e2 = bltin_get_ent (args[1]);
    if (ent_type (e2) == MATRIX_DENSE_REAL)
    {
      m = class_matrix_real (e2);

      // ignore empty matrix
      if (m->nrow * m->ncol == 0)
        m = 0;
    }
  }

  // did user provide 'std'
  if (nargs > 2)
  {
    e3 = bltin_get_ent (args[2]);
    if (ent_type (e3) == MATRIX_DENSE_REAL)
    {
      s = class_matrix_real (e3);

      // ignore empty matrix
      if (s->nrow * s->ncol == 0)
        s = 0;
    }
  }

  // user has to specify either both m and s, or none
  if (!m || !s)
  {
    m = 0;
    s = 0;
  }

  if (ent_type (e1) == BTREE)
  {
    // val
    node = btree_FindNode (ent_data (e1), RLAB_NAME_STAT_VALUE);
    if (node != 0)
      x = class_matrix_real (var_ent (node));
    else
      rerror ("skew: first argument has to be a list <<" RLAB_NAME_STAT_VALUE ";" RLAB_NAME_STAT_WEIGHT ">> or a vector!");

    // wgt
    node = btree_FindNode (ent_data (e1), RLAB_NAME_STAT_WEIGHT);
    if (node != 0)
      t = class_matrix_real (var_ent (node));
    else
      rerror ("skew: first argument has to be a list <<" RLAB_NAME_STAT_VALUE ";" RLAB_NAME_STAT_WEIGHT ">> or a vector!");

    nr = x->nrow;
    nc = x->ncol;
    tr = t->nrow;
    tc = t->ncol;
    if ((nc != 1 && nr!=1) || (tc != 1 && tc != 1))
      rerror ("skew: " RLAB_NAME_STAT_VALUE " and " RLAB_NAME_STAT_WEIGHT " have to be real vectors!");
    if (nc*nr != tc*tr)
      rerror ("skew: " RLAB_NAME_STAT_VALUE " and " RLAB_NAME_STAT_WEIGHT " must have same size!");

    // m and s provided
    if (m && s)
      w = mdr_CreateScalar ( gsl_stats_wskew_m_sd ( MDRPTR(t), 1, MDRPTR(x), 1, nr, MdrV0(m,0), MdrV0(s,0) ) );
    else
      w = mdr_CreateScalar ( gsl_stats_wskew      ( MDRPTR(t), 1, MDRPTR(x), 1, nr ) );
  }
  else if (ent_type (e1) == MATRIX_DENSE_REAL)
  {
    x = class_matrix_real (e1);
    nr = x->nrow;
    nc = x->ncol;

    if (m && s)
    {
      if (nr == 1 || nc == 1)
        w = mdr_CreateScalar (gsl_stats_skew_m_sd (MDRPTR(x), 1, nr*nc, MdrV0(m,0), MdrV0(s,0)));
      else
      {
        w = mdr_Create (1, nc);
        for (i=0; i < nc; i++)
          MdrV0 (w, i) = gsl_stats_skew_m_sd ( &MdrV0(x,i*nr), 1, nr, MdrV0(m,0), MdrV0(s,0));
      }
    }
    else
    {
      if (nr == 1 || nc == 1)
        w = mdr_CreateScalar (gsl_stats_skew (MDRPTR(x), 1, nr*nc));
      else
      {
        w = mdr_Create (1, nc);
        for (i=0; i < nc; i++)
          MdrV0 (w, i) = gsl_stats_skew (MDRPTR(x) + i*nr, 1, nr);
      }
    }
  }
  else
    rerror ("skew: first argument has to be a list <<" RLAB_NAME_STAT_VALUE ";" RLAB_NAME_STAT_WEIGHT ">> or a vector!");

  ent_Clean (e1);
  ent_Clean (e2);
  ent_Clean (e3);

  return ent_Assign_Rlab_MDR(w);
}

//
// (weighted) kurtosis with optional fixed mean and std
//
Ent *
ent_kurtosis (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *e3=0;
  MDR *x=0, *t=0, *w=0;
  MDR *m=0, *s=0;
  int nr, nc, tr, tc, i;

  ListNode *node;

  if (nargs >= 1)
    e1 = bltin_get_ent (args[0]);
  else
    rerror ("kurtosis: first argument has to be a list <<" RLAB_NAME_STAT_VALUE ";" RLAB_NAME_STAT_WEIGHT ">> or a vector!");

  // did user provide 'mean'
  if (nargs > 1)
  {
    e2 = bltin_get_ent (args[1]);
    if (ent_type (e2) == MATRIX_DENSE_REAL)
    {
      m = class_matrix_real (e2);

      // ignore empty matrix
      if (m->nrow * m->ncol == 0)
        m = 0;
    }
  }

  // did user provide 'std'
  if (nargs > 2)
  {
    e3 = bltin_get_ent (args[2]);
    if (ent_type (e3) == MATRIX_DENSE_REAL)
    {
      s = class_matrix_real (e3);

      // ignore empty matrix
      if (s->nrow * s->ncol == 0)
        s = 0;
    }
  }

  // user has to specify either both m and s, or none
  if (!m || !s)
  {
    m = 0;
    s = 0;
  }

  if (ent_type (e1) == BTREE)
  {
    // val
    node = btree_FindNode (ent_data (e1), RLAB_NAME_STAT_VALUE);
    if (node != 0)
      x = class_matrix_real (var_ent (node));
    else
      rerror ("kurtosis: first argument has to be a list <<" RLAB_NAME_STAT_VALUE ";" RLAB_NAME_STAT_WEIGHT ">> or a vector!");

    // wgt
    node = btree_FindNode (ent_data (e1), RLAB_NAME_STAT_WEIGHT);
    if (node != 0)
      t = class_matrix_real (var_ent (node));
    else
      rerror ("kurtosis: first argument has to be a list <<" RLAB_NAME_STAT_VALUE ";" RLAB_NAME_STAT_WEIGHT ">> or a vector!");

    nr = x->nrow;
    nc = x->ncol;
    tr = t->nrow;
    tc = t->ncol;
    if ((nc != 1 && nr!=1) || (tc != 1 && tc != 1))
      rerror ("kurtosis: " RLAB_NAME_STAT_VALUE " and " RLAB_NAME_STAT_WEIGHT " have to be real vectors!");
    if (nc*nr != tc*tr)
      rerror ("kurtosis: " RLAB_NAME_STAT_VALUE " and " RLAB_NAME_STAT_WEIGHT " must have same size!");

    // m and s provided
    if (m && s)
      w = mdr_CreateScalar ( gsl_stats_wkurtosis_m_sd ( MDRPTR(t), 1, MDRPTR(x), 1, nr, MdrV0(m,0), MdrV0(s,0) ) );
    else
      w = mdr_CreateScalar ( gsl_stats_wkurtosis      ( MDRPTR(t), 1, MDRPTR(x), 1, nr ) );
  }
  else if (ent_type (e1) == MATRIX_DENSE_REAL)
  {
    x = class_matrix_real (e1);
    nr = x->nrow;
    nc = x->ncol;

    if (m && s)
    {
      if (nr == 1 || nc == 1)
        w = mdr_CreateScalar (gsl_stats_kurtosis_m_sd (MDRPTR(x), 1, nr*nc, MdrV0(m,0), MdrV0(s,0)));
      else
      {
        w = mdr_Create (1, nc);
        for (i=0; i < nc; i++)
          MdrV0 (w, i) = gsl_stats_kurtosis_m_sd ( &MdrV0(x,i*nr), 1, nr, MdrV0(m,0), MdrV0(s,0));
      }
    }
    else
    {
      if (nr == 1 || nc == 1)
        w = mdr_CreateScalar (gsl_stats_kurtosis (MDRPTR(x), 1, nr*nc));
      else
      {
        w = mdr_Create (1, nc);
        for (i=0; i < nc; i++)
          MdrV0 (w, i) = gsl_stats_kurtosis (&MdrV0(x,i*nr), 1, nr);
      }
    }
  }
  else
    rerror ("kurtosis: first argument has to be a list <<" RLAB_NAME_STAT_VALUE ";" RLAB_NAME_STAT_WEIGHT ">> or a vector!");

  ent_Clean (e1);
  ent_Clean (e2);
  ent_Clean (e3);

  return ent_Assign_Rlab_MDR(w);
}



//
// median
//
Ent *
ent_median (int nargs, Datum args[])
{
  Ent *X=0;
  MDR *x=0, *w;
  int nr, nc, i;

  if(nargs != 1)
  {
    fprintf(stdout,
            "median: Finds a column-wise median of a real data matrix or\n");
    fprintf(stdout,
            "median: a median of a data vector.\n");
    fprintf(stdout,
            "median: Format:\n");
    fprintf(stdout,
            "median:   y = median(x),\n");
    fprintf(stdout,
            "median: where 'x' is a data vector or a real data matrix.\n");
    rerror ("one argument required!\n");
  }

  X = bltin_get_ent (args[0]);
  if (ent_type (X) != MATRIX_DENSE_REAL)
    rerror ("median: argument 'x' has to be MATRIX-DENSE-REAL");

  x = class_matrix_real (X);

  nr = MNR (x);
  nc = MNC (x);

  if (nr*nc == 0)
    rerror ("median: argument 'x' has to be MATRIX-DENSE-REAL");

  if (nr == 1 || nc == 1)
  {
    // a row- or column-vector
    w = mdr_CreateScalar (gsl_stats_median_from_sorted_data (MDRPTR(x), 1, nr * nc));
  }
  else
  {
    // find column-wise means
    w = mdr_Create (1, nc);
    for (i = 1; i <= nc; i++)
    {
      MdrV1 (w, i) = gsl_stats_median_from_sorted_data ( &MdrV0(x,(i-1)*nr), 1, nr);
    }
  }

  return ent_Assign_Rlab_MDR(w);
}

//
// quantiles
//
#undef THIS_SOLVER
#define THIS_SOLVER "quantile"
Ent *
ent_quantile (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0;
  MDR *x1=0, *x=0, *w=0, *f=0;
  int nf, n, i, ikill_x=0;

  if(nargs != 2)
  {
    fprintf(stderr, THIS_SOLVER ":" RLAB_ERROR_TWO_ARG_REQUIRED);
    goto _quantile_end;
  }

  //
  // x -> data vector sorted or not
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) != MATRIX_DENSE_REAL)
  {
    fprintf(stderr, THIS_SOLVER ":" RLAB_ERROR_ARG1_MDR_VECTOR);
    goto _quantile_end;
  }
  x1 = class_matrix_real (e1);
  if (!EQVECT(x1))
  {
    fprintf(stderr, THIS_SOLVER ":" RLAB_ERROR_ARG1_MDR_VECTOR);
    goto _quantile_end;
  }
  n = SIZE(x1);
  if (!mdr_vector_issorted(x1))
  {
    // sort 'x' if not sorted
    x = mdr_Float_BF (x1);
    MDR * x_idx = mdr_Create(1, n);
    for (i=0; i<n; i++)
      MdrV0(x_idx,i) = i;
    // sort
    r_sort ((double *) MDRPTR(x), 0, n-1, (double *)MDRPTR(x_idx));
    mdr_Destroy (x_idx);
    ikill_x = 1;
  }
  else
    x = x1;

  //
  // q -> quantile vector sorted or not
  //
  e2 = bltin_get_ent (args[1]);
  if (ent_type (e2) != MATRIX_DENSE_REAL)
  {
    fprintf(stderr, THIS_SOLVER ":" RLAB_ERROR_ARG2_MDR_VECTOR);
    goto _quantile_end;
  }
  f = class_matrix_real (e2);
  if (!EQVECT(f))
  {
    fprintf(stderr, THIS_SOLVER ":" RLAB_ERROR_ARG2_MDR_VECTOR);
    goto _quantile_end;
  }
  if (!mdr_vector_isbounded_from_below(f, 0))
  {
    fprintf(stderr, THIS_SOLVER ":" RLAB_ERROR_ARG2_MDR_VECTOR);
    fprintf(stderr, THIS_SOLVER ": 2nd argument must be bounded from below by 0!\n" );
    goto _quantile_end;
  }
  if (!mdr_vector_isbounded_from_above(f, 1))
  {
    fprintf(stderr, THIS_SOLVER ":" RLAB_ERROR_ARG2_MDR_VECTOR);
    fprintf(stderr, THIS_SOLVER ": 2nd argument must be bounded from above by 1!\n" );
    goto _quantile_end;
  }
  nf = SIZE(f);

  int nn = mdr_vector_isnan(x);
  if (nn>0)
  {
    if (rlab_sort_nans()==RLAB_SORT_NAN_ONBOTTOM)
    {
      w = mdr_Create_SameSize (f);
      for (i=0; i < nf; i++)
      {
        // a row- or column-vector
        MdrV0(w,i) = gsl_stats_quantile_from_sorted_data (&MdrV0(x,nn), 1, n-nn, mdrV0(f,i));
      }
    }
    else if (rlab_sort_nans()==RLAB_SORT_NAN_ONTOP)
    {
      w = mdr_Create_SameSize (f);
      for (i=0; i < nf; i++)
      {
        // a row- or column-vector
        MdrV0(w,i) = gsl_stats_quantile_from_sorted_data (&MdrV0(x,0), 1, n-nn, mdrV0(f,i));
      }
    }
  }
  else
  {
    w = mdr_Create_SameSize (f);
    for (i = 0; i < nf; i++)
    {
      // a row- or column-vector
      MdrV0(w,i) = gsl_stats_quantile_from_sorted_data (MDRPTR(x), 1, n, mdrV0(f,i));
    }
  }

  if (ikill_x)
    mdr_Destroy(x);

_quantile_end:

  ent_Clean (e1);
  ent_Clean (e2);
  return ent_Assign_Rlab_MDR(w);
}


//
// combinatorics toolkit: permutations and combinations
//
#undef THIS_SOLVER
#define THIS_SOLVER "invperm"
Ent *
ent_invperm (int nargs, Datum args[])
{
  Ent *e1=0;
  MDR *x, *w;
  size_t i;

  gsl_set_error_handler_off ();

  gsl_permutation ix;
  ix.size = 0;
  ix.data = 0;
  gsl_permutation jx;
  jx.size = 0;
  jx.data = 0;

  //
  // get X
  //
  if (nargs != 1)
  {
    fprintf (stdout,
             THIS_SOLVER ": An inverse of a permutation.\n");
    fprintf (stdout,
             THIS_SOLVER ":   i = invperm( j ) .\n");
    rerror ("One argument required!");
  }

  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) != MATRIX_DENSE_REAL)
    rerror(THIS_SOLVER ": Argument is not valid permutation!");

  x = ent_data (e1);
  if ((x->nrow) * (x->ncol) < 2)
    rerror(THIS_SOLVER ": Invalid length of the permutation!");
  ix.size = (x->nrow) * (x->ncol);
  ix.data = (size_t *) GC_malloc((x->nrow) * (x->ncol) * sizeof(size_t));
  if (x->type == RLAB_TYPE_INT32)
  {
    for (i=0; i < ix.size; i++)
      ix.data[i] = (size_t) MdiV0(x,i) - 1L;
  }
  else
  {
    for (i=0; i < ix.size; i++)
      ix.data[i] = (size_t) MdrV0(x,i) - 1L;
  }
  if (gsl_permutation_valid(&ix))
    rerror(THIS_SOLVER ": Argument is not a valid permutation!");

  jx.size = ix.size;
  jx.data = (size_t *) GC_malloc((x->nrow) * (x->ncol) * sizeof(size_t));

  gsl_permutation_inverse(&jx, &ix);

  w = mdi_Create(1, ix.size);

  for (i=0; i < ix.size; i++)
  {
    MdiV0(w,i) = (int) jx.data[i] + 1;
  }

  GC_free(ix.data);
  GC_free(jx.data);
  ent_Clean (e1);

  return ent_Assign_Rlab_MDR(w);
}

#undef THIS_SOLVER
#define THIS_SOLVER "revperm"
Ent *
ent_revperm (int nargs, Datum args[])
{
  Ent *e1=0;
  MDR *x, *w;
  size_t i;

  gsl_set_error_handler_off ();

  gsl_permutation ix;
  ix.size = 0;
  ix.data = 0;

  //
  // get X
  //
  if (nargs != 1)
  {
    fprintf (stdout,
             THIS_SOLVER ": An reverse of a permutation.\n");
    fprintf (stdout,
             THIS_SOLVER ":   i = revperm( j ) .\n");
    rerror ("One argument required!");
  }
  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) == MATRIX_DENSE_REAL)
  {
    x = ent_data (e1);
    if ((x->nrow) * (x->ncol) < 2)
      rerror(THIS_SOLVER ": Invalid length of the permutation!");
    ix.size = (x->nrow) * (x->ncol);
    ix.data = (size_t *) GC_malloc((x->nrow) * (x->ncol) * sizeof(size_t));
    if (x->type == RLAB_TYPE_INT32)
    {
      for (i=0; i < ix.size; i++)
        ix.data[i] = (size_t) MdiV0(x,i) - 1L;
    }
    else
    {
      for (i=0; i < ix.size; i++)
        ix.data[i] = (size_t) MdrV0(x,i) - 1L;
    }
    if (gsl_permutation_valid(&ix))
      rerror(THIS_SOLVER ": Argument is not a valid permutation!");
  }
  else
    rerror(THIS_SOLVER ": Argument is not valid permutation!");

  w = mdi_Create (1, (int) ix.size);

  gsl_permutation_reverse(&ix);
  for (i=0; i < ix.size; i++)
    MdiV0(w,i) = (int) ix.data[i] + 1;

  GC_free (ix.data);
  ent_Clean (e1);

  return ent_Assign_Rlab_MDR(w);
}


#undef THIS_SOLVER
#define THIS_SOLVER "nextperm"
Ent *
ent_nextperm (int nargs, Datum args[])
{
  Ent *e1=0;
  MDR *x=0, *w;
  int i;

  gsl_set_error_handler_off ();

  gsl_permutation ix;
  ix.size = 0;
  ix.data = 0;

  //
  // get X
  //
  if (nargs != 1)
  {
    fprintf (stdout,
             THIS_SOLVER ": Next permutation in a lexixographic ordering of\n");
    fprintf (stdout,
             THIS_SOLVER ": permutations. Format:\n");
    fprintf (stdout,
             THIS_SOLVER ":   i = nextperm( j ) .\n");
    rerror ("One argument required!");
  }
  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) == MATRIX_DENSE_REAL)
  {
    x = ent_data (e1);
    if ((x->nrow) * (x->ncol) < 2)
      rerror(THIS_SOLVER ": Invalid length of the permutation!");
    ix.size = (x->nrow) * (x->ncol);
    ix.data = (size_t *) GC_malloc((x->nrow) * (x->ncol) * sizeof(size_t));
    if (x->type == RLAB_TYPE_INT32)
    {
      for (i=0; i < ix.size; i++)
        ix.data[i] = (size_t) MdiV0(x,i) - 1L;
    }
    else
    {
      for (i=0; i < ix.size; i++)
        ix.data[i] = (size_t) MdrV0(x,i) - 1L;
    }
    if (gsl_permutation_valid(&ix))
      rerror(THIS_SOLVER ": Argument is not a valid permutation!");
  }
  else
    rerror(THIS_SOLVER ": Argument is not valid permutation!");

  w = mdi_Create(1, ix.size);
  if (gsl_permutation_next(&ix)==GSL_SUCCESS)
  {
    for (i=0; i < ix.size; i++)
      MdiV0(w,i) = (int) ix.data[i] + 1;
  }

  GC_free (ix.data);
  ent_Clean (e1);

  return ent_Assign_Rlab_MDR(w);
}

#undef THIS_SOLVER
#define THIS_SOLVER "prevperm"
Ent *
ent_prevperm (int nargs, Datum args[])
{
  Ent *e1=0;
  MDR *x, *w;
  size_t i;

  gsl_set_error_handler_off ();

  gsl_permutation ix;
  ix.size = 0;
  ix.data = 0;

  //
  // get X
  //
  if (nargs != 1)
  {
    fprintf (stdout,
             THIS_SOLVER ": A previous permutation in lexixographic ordering of\n");
    fprintf (stdout,
             THIS_SOLVER ": permutations. Format:\n");
    fprintf (stdout,
             THIS_SOLVER ":   i = prevperm( j ) .\n");
    rerror ("One argument required!");
  }
  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) == MATRIX_DENSE_REAL)
  {
    x = ent_data (e1);
    if ((x->nrow) * (x->ncol) < 2)
      rerror(THIS_SOLVER ": Invalid length of the permutation!");
    ix.size = (x->nrow) * (x->ncol);
    ix.data = (size_t *) GC_malloc((x->nrow) * (x->ncol) * sizeof(size_t));
    if (x->type == RLAB_TYPE_INT32)
    {
      for (i=0; i < ix.size; i++)
        ix.data[i] = (size_t) MdiV0(x,i) - 1L;
    }
    else
    {
      for (i=0; i < ix.size; i++)
        ix.data[i] = (size_t) MdrV0(x,i) - 1L;
    }
    if (gsl_permutation_valid(&ix))
      rerror(THIS_SOLVER ": Argument is not a valid permutation!");
  }
  else
    rerror(THIS_SOLVER ": Argument is not valid permutation!");

  w = mdi_Create(1, ix.size);
  if (gsl_permutation_prev(&ix)==GSL_SUCCESS)
  {
    for (i=0; i < ix.size; i++)
      MdiV0(w,i) = (int) ix.data[i] + 1;
  }

  GC_free (ix.data);

  ent_Clean (e1);

  return ent_Assign_Rlab_MDR(w);
}

#undef THIS_SOLVER
#define THIS_SOLVER "validperm"
Ent *
ent_isaperm (int nargs, Datum args[])
{
  Ent *e1=0;
  MDR *x=0, *w;
  int i;

  gsl_permutation ix;
  ix.size = 0;
  ix.data = 0;

  gsl_set_error_handler_off ();

  //
  // get X
  //
  if (nargs != 1)
  {
    fprintf (stdout,
             THIS_SOLVER ": Test whether an integer vector is a valid permutation.\n");
    fprintf (stdout,
             THIS_SOLVER ":   i = validperm( j ) .\n");
    rerror ("One argument required!");
  }
  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) == MATRIX_DENSE_REAL)
  {
    x = ent_data (e1);
    if ((x->nrow) * (x->ncol) < 2)
      rerror(THIS_SOLVER ": Invalid length of the permutation!");
    ix.size = (x->nrow) * (x->ncol);
    ix.data = (size_t *) GC_malloc((x->nrow) * (x->ncol) * sizeof(size_t));
    if (x->type == RLAB_TYPE_INT32)
    {
      for (i=0; i < ix.size; i++)
        ix.data[i] = (size_t) MdiV0(x,i) - 1L;
    }
    else
    {
      for (i=0; i < ix.size; i++)
        ix.data[i] = (size_t) MdrV0(x,i) - 1L;
    }
  }

  w = mdr_CreateScalar( !gsl_permutation_valid(&ix) );

  GC_free (ix.data);
  ent_Clean (e1);

  return ent_Assign_Rlab_MDR(w);
}

#undef THIS_SOLVER
#define THIS_SOLVER "mulperm"
Ent *
ent_mulperm (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0;
  MDR *x=0, *y=0, *w=0;
  size_t i;

  gsl_set_error_handler_off ();

  gsl_permutation ix;
  ix.size = 0;
  ix.data = 0;
  gsl_permutation iy;
  iy.size = 0;
  iy.data = 0;
  gsl_permutation iz;
  iz.size = 0;
  iz.data = 0;

  //
  // get X
  //
  if (nargs != 2)
  {
    fprintf (stdout,
             THIS_SOLVER ": Multiply two permutation.\n");
    fprintf (stdout,
             THIS_SOLVER ": Format:\n");
    fprintf (stdout,
             THIS_SOLVER ":   k = mulperm( i, j ) .\n");
    rerror ("Two argument required!");
  }

  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) == MATRIX_DENSE_REAL)
  {
    x = ent_data (e1);
    if ((x->nrow) * (x->ncol) < 2)
      rerror(THIS_SOLVER ": Invalid length of the permutation!");
    ix.size = (x->nrow) * (x->ncol);
    ix.data = (size_t *) GC_malloc((x->nrow) * (x->ncol) * sizeof(size_t));
    if (x->type == RLAB_TYPE_INT32)
    {
      for (i=0; i < ix.size; i++)
        ix.data[i] = (size_t) MdiV0(x,i) - 1L;
    }
    else
    {
      for (i=0; i < ix.size; i++)
        ix.data[i] = (size_t) MdrV0(x,i) - 1L;
    }
    if (gsl_permutation_valid(&ix))
      rerror(THIS_SOLVER ": Argument is not a valid permutation!");
  }
  else
    rerror(THIS_SOLVER ": Argument is not valid permutation!");

  e2 = bltin_get_ent (args[1]);
  if (ent_type (e2) == MATRIX_DENSE_REAL)
  {
    y = ent_data (e2);
    if ((y->nrow) * (y->ncol) < 2)
      rerror(THIS_SOLVER ": Invalid length of the permutation!");
    iy.size = (y->nrow) * (y->ncol);
    iy.data = (size_t *) GC_malloc((y->nrow) * (y->ncol) * sizeof(size_t));
    if (y->type == RLAB_TYPE_INT32)
    {
      for (i=0; i < iy.size; i++)
        iy.data[i] = (size_t) MdiV0(y,i) - 1L;
    }
    else
    {
      for (i=0; i < iy.size; i++)
        iy.data[i] = (size_t) MdrV0(y,i) - 1L;
    }
    if (gsl_permutation_valid(&iy))
      rerror(THIS_SOLVER ": Argument is not a valid permutation!");
  }
  else
    rerror(THIS_SOLVER ": Argument is not valid permutation!");

  if ( ix.size != iy.size )
    rerror("mulperm: Permutations have to be of the same size!");

  w = mdi_Create(1, ix.size);
  iz.size = ix.size;
  iz.data = (size_t *) GC_malloc(x->nrow * x->ncol * sizeof(size_t));

  gsl_permutation_mul(&iz, &ix, &iy);
  for (i=0; i < ix.size; i++)
  {
    MdiV0(w,(int) i) = (int) iz.data[i] + 1;
  }

  GC_free(iz.data);
  GC_free(iy.data);
  GC_free(ix.data);
  ent_Clean (e1);
  ent_Clean (e2);

  return ent_Assign_Rlab_MDR(w);
}


#undef THIS_SOLVER
#define THIS_SOLVER "prevcomb"
Ent *
ent_prevcomb (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0;
  MDR *x, *w;
  size_t i;

  gsl_set_error_handler_off ();

  gsl_combination ix;
  ix.n = 0;
  ix.k = 0;
  ix.data = 0;

  //
  // get X
  //
  if (nargs != 2)
  {
    fprintf (stdout,
             THIS_SOLVER ": A previous combination in lexixographic ordering of\n");
    fprintf (stdout,
             THIS_SOLVER ": combinations from a set of size n. Format:\n");
    fprintf (stdout,
             THIS_SOLVER ":   i = prevcomb( n, k ) .\n");
    rerror ("Two argument required!");
  }

  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) == MATRIX_DENSE_REAL)
    ix.n = (size_t) class_double (e1);
  else
    rerror(THIS_SOLVER ": Invalid size of the set!");
  if (ix.n < 2)
    rerror(THIS_SOLVER ": Invalid size of the set!");

  e2 = bltin_get_ent (args[1]);
  if (ent_type (e2) != MATRIX_DENSE_REAL)
    rerror(THIS_SOLVER ": Argument is not valid combination!");

  x = ent_data (e2);
  ix.k = (size_t) (x->nrow) * (x->ncol);
  ix.data = (size_t *) GC_malloc(MNR(x) * MNC(x) * sizeof(size_t));
  if (x->type == RLAB_TYPE_INT32)
  {
    for (i=0; i < ix.k; i++)
      ix.data[i] = (size_t) MdiV0(x,(int) i) - (size_t) 1;
  }
  else
  {
    for (i=0; i < ix.k; i++)
      ix.data[i] = (size_t) MdrV0(x,(int) i) - (size_t) 1;
  }
  if (gsl_combination_valid(&ix))
    rerror("prevcomb: Argument is not a valid combination!");

  w = mdi_Create(1, ix.k);
  if (gsl_combination_prev(&ix)==GSL_SUCCESS)
  {
    for (i=0; i < ix.k; i++)
      MdiV0(w, (int) i) = (int) ix.data[i] + 1;
  }

  GC_free(ix.data);

  ent_Clean (e1);
  ent_Clean (e2);

  return ent_Assign_Rlab_MDR(w);
}

#undef THIS_SOLVER
#define THIS_SOLVER "nextcomb"
Ent *
ent_nextcomb (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0;
  MDR *x, *w;
  size_t i;

  gsl_set_error_handler_off ();

  gsl_combination ix;
  ix.n = 0;
  ix.k = 0;
  ix.data = 0;

  //
  // get X
  //
  if (nargs != 2)
  {
    fprintf (stdout,
             THIS_SOLVER ": Next combination in lexixographic ordering of\n");
    fprintf (stdout,
             THIS_SOLVER ": combinations from a set of size n. Format:\n");
    fprintf (stdout,
             THIS_SOLVER ":   i = nextcomb( n, j ) .\n");
    rerror ("Two argument required!");
  }

  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) != MATRIX_DENSE_REAL)
    rerror(THIS_SOLVER ": Invalid size of the set!");
  ix.n = (size_t) class_double (e1);
  if (ix.n < 2)
    rerror(THIS_SOLVER ": Invalid size of the set!");

  e2 = bltin_get_ent (args[1]);
  if (ent_type (e2) != MATRIX_DENSE_REAL)
    rerror(THIS_SOLVER ": Argument is not valid combination!");

  x = ent_data (e2);
  ix.k = (size_t) (x->nrow) * (x->ncol);
  ix.data = (size_t *) GC_malloc(MNR(x) * MNC(x) * sizeof(size_t));
  if (x->type == RLAB_TYPE_INT32)
  {
    for (i=0; i < ix.k; i++)
      ix.data[i] = (size_t) MdiV0(x,(int) i) - (size_t) 1;
  }
  else
  {
    for (i=0; i < ix.k; i++)
      ix.data[i] = (size_t) MdrV0(x,(int) i) - (size_t) 1;
  }
  if (gsl_combination_valid(&ix))
    rerror("prevcomb: Argument is not a valid combination!");

  w = mdi_Create(1, ix.k);
  if (gsl_combination_next(&ix)==GSL_SUCCESS)
  {
    for (i=0; i < ix.k; i++)
      MdiV0(w, (int) i) = (int) ix.data[i] + 1;
  }

  GC_free(ix.data);
  ent_Clean (e1);
  ent_Clean (e2);

  return ent_Assign_Rlab_MDR(w);
}

#undef THIS_SOLVER
#define THIS_SOLVER "validcomb"
Ent *
ent_iscomb (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0;
  MDR *x, *w;
  size_t i;

  gsl_set_error_handler_off ();

  gsl_combination ix;
  ix.n = 0;
  ix.k = 0;
  ix.data = 0;

  //
  // get X
  //
  if (nargs != 2)
  {
    fprintf (stdout,
             THIS_SOLVER ": Test whether the combination is valid, given the set size n.\n");
    fprintf (stdout,
             THIS_SOLVER ": Format:\n");
    fprintf (stdout,
             THIS_SOLVER ":   i = validcomb( n, j ) .\n");
    rerror ("Two argument required!");
  }

  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) != MATRIX_DENSE_REAL)
    rerror(THIS_SOLVER ": Invalid size of the set!");
  ix.n = (size_t) class_double (e1);
  if (ix.n < 2)
    rerror(THIS_SOLVER ": Invalid size of the set!");

  e2 = bltin_get_ent (args[1]);
  if (ent_type (e2) != MATRIX_DENSE_REAL)
    rerror(THIS_SOLVER ": Argument is not valid combination!");

  x = ent_data (e2);
  ix.k = (size_t) (x->nrow) * (x->ncol);
  ix.data = (size_t *) GC_malloc(MNR(x) * MNC(x) * sizeof(size_t));
  if (x->type == RLAB_TYPE_INT32)
  {
    for (i=0; i < ix.k; i++)
      ix.data[i] = (size_t) MdiV0(x,(int) i) - (size_t) 1;
  }
  else
  {
    for (i=0; i < ix.k; i++)
      ix.data[i] = (size_t) MdrV0(x,(int) i) - (size_t) 1;
  }
  w = mdr_CreateScalar( !gsl_combination_valid(&ix) );

  GC_free(ix.data);

  ent_Clean (e1);
  ent_Clean (e2);

  return ent_Assign_Rlab_MDR(w);
}

Ent *
ent_ncomb (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0;
  MDR *w;
  unsigned int i, j;

  gsl_set_error_handler_off ();

  //
  // get X
  //
  if (nargs != 2)
  {
    fprintf (stdout,
             "ncomb: Combinatorial factor  n choose m.\n");
    fprintf (stdout,
             "ncomb: Format:\n");
    fprintf (stdout,
             "ncomb:   i = ncomb( n, m ) .\n");
    rerror ("Two argument required!");
  }

  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) != MATRIX_DENSE_REAL)
    rerror("ncomb: Invalid first argument!");
  i = (int) class_double (e1);

  e2 = bltin_get_ent (args[1]);
  if (ent_type (e2) != MATRIX_DENSE_REAL)
    rerror("ncomb: Invalid second argument!");
  j = (int) class_double (e2);

  w = mdr_CreateScalar( gsl_sf_choose(i, j) );

  ent_Clean (e1);
  ent_Clean (e2);

  return ent_Assign_Rlab_MDR(w);
}

#undef  THIS_SOLVER
#define THIS_SOLVER "hist"
Ent *
ent_histogram_create (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *e3=0, *eo=0, *ent_info=0;
  MDR *x=0, *r=0, *range=0, *bin=0;
  int i, nbins=10, size, count_pos_inf=0, count_neg_inf=0, count_nan=0, count_below_min=0;
  int count_above_max=0, offs=0;
  int create_new_hist=0, duplicate_hist=0;
  double xmin=create_nan(), xmax=create_nan();

  FILE *rlab_stderr = (!RLAB_STDERR_DS) ? stderr : RLAB_STDERR_DS;
  Btree *bw=0, *b_info=0;

  // create a histogram structure, but do not fill it
  gsl_histogram *hist=0;
  ListNode *node_info=0;

  gsl_set_error_handler_off ();

  if (nargs != 1 && nargs != 2 && nargs != 3)
  {
    fprintf (stdout, THIS_SOLVER ": 1. Construct a gsl-type histogram of an input real vector data.\n");
    fprintf (stdout, THIS_SOLVER ": Format:\n");
    fprintf (stdout, THIS_SOLVER ":   h = hist(data, range, n) .\n");
    fprintf (stdout, THIS_SOLVER ": 2. Update an existing gsl-type histogram with new real vector data.\n");
    fprintf (stdout, THIS_SOLVER ": Format:\n");
    fprintf (stdout, THIS_SOLVER ":   hist(h, data) .\n");
    rerror (THIS_SOLVER ": " RLAB_ERROR_ONE_TO_THREE_ARG_REQUIRED);
  }

  //
  // data: vector, or matrix all in one histogram
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) == MATRIX_DENSE_REAL)
  {
    create_new_hist = 1;
    x = ent_data(e1);
    size = SIZE(x);
    if (size > 0)
    {
      xmin = rlab_dmin_mdr(x);
      xmax = rlab_dmax_mdr(x);
    }

    if (nargs > 1)
    {
      // second argument:
      // scalar: nbins
      // range: two element vector, then [min, max]
      //        longer vector, then [r[0], r[1], r[2], ... r[n]]
      e2 = bltin_get_ent (args[1]);
      if (ent_type (e2) == MATRIX_DENSE_REAL)
      {
        r = ent_data(e2);
        if (SIZE(r)<1)
          r=0;
        else
        {
          if (!IsFinite(r))
            r = 0;
        }
        if (SIZE(r)==1)
        {
          // nbins, we do not expect any more arguments
          if (mdiV0(r,0) > 1)
            nbins = mdiV0(r,0);
          offs = 1;
          r = 0;
        }
        else if (SIZE(r)==2)
        {
          // [xmin, xmax], we possibly may get 'nbins' as the third argument
          xmin = MIN(mdrV0(r,0), mdrV0(r,1));
          xmax = MAX(mdrV0(r,0), mdrV0(r,1));
          offs = 0;
          r = 0;
        }
      }
    }

    if (nargs > 2)
    {
      e3 = bltin_get_ent (args[2]);
      if (ent_type (e3) == MATRIX_DENSE_REAL)
      {
        int d = class_double(e3);
        if (d >= 2)
          nbins = d;
      }
    }

    if (!r )
    {
      // construct uniform range
      if ((!isnand(xmin)) && (!isnand(xmax)))
      {
        range = mdr_Create(1, nbins+1);
        for (i=0; i<=nbins; i++)
          MdrV0(range,i) = xmin + (xmax-xmin)/((double) (nbins-offs))*((double) i - 0.5*offs);
      }
      else
      {
        fprintf (stdout, THIS_SOLVER ": Cannot compute histogram range without data!\n");
        fprintf (stdout, THIS_SOLVER ": To create range for empty histogram, please provide array "
            " range=[r1, r2, r3, ...] as second argument !\n");
        goto _exit_hist;
      }
    }
    else
    {
      range = (MDR *) mdr_Float_BF(r);
      nbins = SIZE(range) - 1;
    }

    // construct bin matrix, a column for each one histogram
    bin = mdr_Create(1,nbins);
    mdr_Zero(bin);
  }
  else if (ent_type (e1) == BTREE)
  {
    ListNode *node=0;

    // locate second argument. this is data added to histogram
    if (nargs != 2)
      goto _exit_hist;  // nothing to do!
    e2 = bltin_get_ent (args[1]);
    if (ent_type(e2)!= MATRIX_DENSE_REAL)
      goto _exit_hist;  // nothing to do!
    x = ent_data(e2);
    if (!MD_TYPE_DOUBLE(x))
      goto _exit_hist;  // nothing to do!
    size = SIZE(x);
    if (size < 1)
      goto _exit_hist;  // nothing to do!

    bw = ent_data(e1);

    // must have h = <<range;bin>>
    //    bin:
    node = btree_FindNode (bw, RLAB_NAME_HIST1D_BIN);
    if (node == 0)
    {
      fprintf (rlab_stderr, THIS_SOLVER ": " RLAB_ERROR_ARG1_HIST "\n");
      goto _exit_hist;
    }
    bin = class_matrix_real (var_ent (node));

    //
    // duplicate entity if it is refered to by more than once. we are changing
    // its content which would affect all variables pointing to it.
    //
    ListNode *var = (ListNode *) (args[0].u.ptr);
    if (e1->refc > 1)
    {
      ent_DecRef (e1);
      e1 = ent_Duplicate (e1);
      listNode_AttachEnt (var, e1);
      duplicate_hist = 1; // we have to duplicate the data separately
    }

    // hist(h, x):
    //  update histogram h with data from x
    bw = ent_data(e1);

    // must have h = <<range;bin>>
    //    bin:
    node = btree_FindNode (bw, RLAB_NAME_HIST1D_BIN);
    if (node == 0)
    {
      fprintf (rlab_stderr, THIS_SOLVER ": " RLAB_ERROR_ARG1_HIST "\n");
      goto _exit_hist;
    }
    if (duplicate_hist)
    {
      var_ent (node) = ent_Duplicate (var_ent (node));
    }
    bin = class_matrix_real (var_ent (node));

    if (!MD_TYPE_DOUBLE(bin))
    {
      fprintf (rlab_stderr, THIS_SOLVER ": " RLAB_ERROR_ARG1_HIST "\n");
      fprintf (rlab_stderr, THIS_SOLVER ": Entry 'bin' must be real matrix!\n");
      goto _exit_hist;
    }
    nbins = SIZE(bin) - 1;
    //    range:
    node = btree_FindNode (bw, RLAB_NAME_HIST1D_RANGE);
    if (node == 0)
    {
      fprintf (rlab_stderr, THIS_SOLVER ": " RLAB_ERROR_ARG1_HIST "\n");
      goto _exit_hist;
    }
    range = class_matrix_real (var_ent (node));
    if (!MD_TYPE_DOUBLE(range))
    {
      fprintf (rlab_stderr, THIS_SOLVER ": " RLAB_ERROR_ARG1_HIST "\n");
      fprintf (rlab_stderr, THIS_SOLVER ": Entry 'range' must be real matrix!\n");
      goto _exit_hist;
    }
    if (duplicate_hist)
    {
      var_ent (node) = ent_Duplicate (var_ent (node));
    }
    range = class_matrix_real (var_ent (node));


    // is there info? if absent create it. otherwise update it
    node_info = btree_FindNode (bw, RLAB_NAME_HIST1D_INFO);
    if (node_info != 0)
    {
      ent_info = var_ent (node_info);
      b_info = ent_data (ent_info);
      node = btree_FindNode (b_info, RLAB_NAME_HIST1D_PINFS);
      if (node != 0)
      {
        if (duplicate_hist)
          var_ent (node) = ent_Duplicate (var_ent (node));
        count_pos_inf = class_int(var_ent (node));
      }

      node = btree_FindNode (b_info, RLAB_NAME_HIST1D_NINFS);
      if (node != 0)
      {
        if (duplicate_hist)
          var_ent (node) = ent_Duplicate (var_ent (node));
        count_neg_inf = class_int(var_ent (node));
      }

      node = btree_FindNode (b_info, RLAB_NAME_HIST1D_NANS);
      if (node != 0)
      {
        if (duplicate_hist)
          var_ent (node) = ent_Duplicate (var_ent (node));
        count_nan = class_int(var_ent (node));
      }

      node = btree_FindNode (b_info, RLAB_NAME_HIST1D_TRASH_MAX);
      if (node != 0)
      {
        if (duplicate_hist)
          var_ent (node) = ent_Duplicate (var_ent (node));
        count_above_max = class_int(var_ent (node));
      }

      node = btree_FindNode (b_info, RLAB_NAME_HIST1D_TRASH_MIN);
      if (node != 0)
      {
        if (duplicate_hist)
          var_ent (node) = ent_Duplicate (var_ent (node));
        count_below_min = class_int(var_ent (node));
      }
    }
  }
  else
    rerror (THIS_SOLVER ": " RLAB_ERROR_ARG1_HIST_OR_MDR_VECTOR);

  // create a fake histogram
  hist = gsl_histogram_alloc(nbins);
  hist->range = MDRPTR(range);
  hist->bin   = MDRPTR(bin);

  // now go over each entry of 'x' and update the histogram bin count
  for (i=0; i<size; i++)
  {
    // ignore nan's
    double d = mdrV0(x,i);
    if(isnand (d))
    {
      count_nan++;
      continue;
    }
    if(d == create_inf())
    {
      count_pos_inf++;
      continue;
    }
    if(d == -create_inf())
    {
      count_neg_inf++;
      continue;
    }
    if (gsl_histogram_increment(hist, d))
    {
      if (d < MdrV0(range,0))
      {
        count_below_min++;
        continue;
      }
      count_above_max++;
    }
  }

  // clean-up histogram
  if (hist)
  {
    hist->range = 0;
    hist->bin   = 0;
    gsl_histogram_free(hist);
  }

  if (create_new_hist)
  {
    //
    // return the result as a list
    //
    bw = btree_Create ();
    // bins
    install (bw, RLAB_NAME_HIST1D_BIN, ent_Assign_Rlab_MDR(bin));
    // range
    install (bw, RLAB_NAME_HIST1D_RANGE, ent_Assign_Rlab_MDR(range));

    b_info = btree_Create ();
    // any nans?
    install (b_info, RLAB_NAME_HIST1D_NANS, ent_Create_Rlab_Int(count_nan));
    // any +infs?
    install (b_info, RLAB_NAME_HIST1D_PINFS, ent_Create_Rlab_Int(count_pos_inf));
    // any -infs?
    install (b_info, RLAB_NAME_HIST1D_NINFS, ent_Create_Rlab_Int(count_neg_inf));
    // any  below min?
    install (b_info, RLAB_NAME_HIST1D_TRASH_MIN, ent_Create_Rlab_Int(count_below_min));
    // any  above max?
    install (b_info, RLAB_NAME_HIST1D_TRASH_MAX, ent_Create_Rlab_Int(count_above_max));
    // install list
    install (bw, RLAB_NAME_HIST1D_INFO, ent_Assign_Rlab_BTREE(b_info));
  }
  else
  {
    ListNode *node=0;
    // update existing histogram _info entry: if necessary create it
    if (!b_info)
    {
      b_info = btree_Create ();
      install (bw, RLAB_NAME_HIST1D_INFO, ent_Assign_Rlab_BTREE(b_info));
    }

    node = btree_FindNode (b_info, RLAB_NAME_HIST1D_PINFS);
    if (node != 0)
    {
      Ent *dummy_ent = var_ent (node);
      MDR *dummy_x   = ent_data(dummy_ent);
      if (MD_TYPE_DOUBLE(dummy_x))
        MdrV0(dummy_x,0) = count_pos_inf;
      else
        MdiV0(dummy_x,0) = count_pos_inf;
    }
    else
      install (b_info, RLAB_NAME_HIST1D_PINFS, ent_Create_Rlab_Int(count_pos_inf));

    node = btree_FindNode (b_info, RLAB_NAME_HIST1D_NINFS);
    if (node != 0)
    {
      Ent *dummy_ent = var_ent (node);
      MDR *dummy_x   = ent_data(dummy_ent);
      if (MD_TYPE_DOUBLE(dummy_x))
        MdrV0(dummy_x,0) = count_pos_inf;
      else
        MdiV0(dummy_x,0) = count_neg_inf;
    }
    else
      install (b_info, RLAB_NAME_HIST1D_NINFS, ent_Create_Rlab_Int(count_neg_inf));

    node = btree_FindNode (b_info, RLAB_NAME_HIST1D_NANS);
    if (node != 0)
    {
      Ent *dummy_ent = var_ent (node);
      MDR *dummy_x   = ent_data(dummy_ent);
      if (MD_TYPE_DOUBLE(dummy_x))
        MdrV0(dummy_x,0) = count_nan;
      else
        MdiV0(dummy_x,0) = count_nan;
    }
    else
      install (b_info, RLAB_NAME_HIST1D_NINFS, ent_Create_Rlab_Int(count_nan));

    node = btree_FindNode (b_info, RLAB_NAME_HIST1D_TRASH_MAX);
    if (node != 0)
    {
      Ent *dummy_ent = var_ent (node);
      MDR *dummy_x   = ent_data(dummy_ent);
      if (MD_TYPE_DOUBLE(dummy_x))
        MdrV0(dummy_x,0) = count_above_max;
      else
        MdiV0(dummy_x,0) = count_above_max;
    }
    else
      install (b_info, RLAB_NAME_HIST1D_TRASH_MIN, ent_Create_Rlab_Int(count_above_max));

    node = btree_FindNode (b_info, RLAB_NAME_HIST1D_TRASH_MIN);
    if (node != 0)
    {
      Ent *dummy_ent = var_ent (node);
      MDR *dummy_x   = ent_data(dummy_ent);
      if (MD_TYPE_DOUBLE(dummy_x))
        MdrV0(dummy_x,0) = count_below_min;
      else
        MdiV0(dummy_x,0) = count_below_min;
    }
    else
      install (b_info, RLAB_NAME_HIST1D_TRASH_MIN, ent_Create_Rlab_Int(count_below_min));
  }

_exit_hist:

  // clean up
  ent_Clean(eo);
  ent_Clean(e1);
  ent_Clean(e2);
  ent_Clean(e3);

  if (create_new_hist)
  {
    return ent_Assign_Rlab_BTREE(bw);
  }
  else
  {
    return ent_Create_Rlab_Success();
  }
}


#undef  THIS_SOLVER
#define THIS_SOLVER "hist2"
Ent * ent_histogram2d_create (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *e3=0, *e4=0, *e5=0;
  MDR *x=0, *r=0, *r2=0, *range1=0, *range2=0, *bin=0;
  int i, nrows, nbins1=10, nbins2=10;
  double x1min=0, x2min=0, x1max=0, x2max=0;

  int ecount=0;

  // create a histogram structure, but do not fill it
  gsl_histogram2d hist2;

  gsl_set_error_handler_off ();

  if (nargs < 1 || nargs > 5)
  {
    fprintf (stdout,
             THIS_SOLVER ": Construct a gsl-type histogram in 2D of an input real vector data.\n");
    fprintf (stdout,
             THIS_SOLVER ": Format:\n");
    fprintf (stdout,
             THIS_SOLVER ":   h = hist2(data, range, n) .\n");
    rerror ("One through five arguments required");
  }

  //
  // data: two-column matrix
  //
  e1 = bltin_get_ent (args[ecount++]);
  if (ent_type (e1) != MATRIX_DENSE_REAL)
    rerror(THIS_SOLVER "hist2: First argument has to be two-column real matrix");

  x = ent_data(e1);
  if (MNC(x) != 2)
    rerror(THIS_SOLVER ": First argument has to be two-column real matrix");
  nrows = MNR(x);
  if (nrows == 0)
    rerror(THIS_SOLVER ": First argument of zero size");

  if (nargs == 1)
  {
    // figure out x{i}min, x{i}max for i=1,2
    x1min = rlab_dmin_vector(nrows, 1, &MdrV0(x,0));
    x1max = rlab_dmax_vector(nrows, 1, &MdrV0(x,0));
    x2max = rlab_dmin_vector(nrows, 1, &MdrV0(x,nrows));
    x2max = rlab_dmax_vector(nrows, 1, &MdrV0(x,nrows));

    // construct uniform range
    // x1:
    range1 = mdr_Create(1, nbins1+1);
    for (i=0; i<=nbins1; i++)
      MdrV0(range1,i) = x1min + (x1max-x1min)/((double) nbins1)*((double) i);
    // x2:
    range2 = mdr_Create(1, nbins2+1);
    for (i=0; i<=nbins2; i++)
      MdrV0(range2,i) = x2min + (x2max-x2min)/((double) nbins2)*((double) i);
    goto _do_hist2;
  }

  // ----------------------------------
  // determining histogram for x1:
  // ----------------------------------
  // scalar: nbins
  // range: two element vector, then [min, max]
  //        longer vector, then [r[0], r[1], r[2], ... r[n]]
  e2 = bltin_get_ent (args[1]);
  if (ent_type (e2) != MATRIX_DENSE_REAL)
    rerror("hist2: Second argument has to be a real/integer vector!");
  r = ent_data(e2);
  if (SIZE(r)==1)
  {
      // nbins1, we do not expect any more arguments for x1
    nbins1 = mdrV0(r,0);
    if (nbins1 < 2)
      rerror("hist2: Second argument as 'nbins1' has to be greater than 2");
      // figure out xmin, xmax
    x1min = rlab_dmin_vector(nrows, 1, &MdrV0(x,0));
    x1max = rlab_dmax_vector(nrows, 1, &MdrV0(x,0));
      // construct uniform range
    range1 = mdr_Create(1, nbins1+1);
    for (i=0; i<=nbins1; i++)
      MdrV0(range1,i) = x1min + (x1max-x1min)/((double) nbins1)*((double) i);

    ecount = 1; // 'nbins1' obtained, 'range1' uniform between x1min and x1max
  }
  else if (SIZE(r)==2)
  {
      // [xmin, xmax], we possibly may get 'nbins' as the third argument
    x1min = MdrV0(r,0);
    x1max = MdrV0(r,1);
    if (x1min > x1max)
      rerror("hist2: Second argument as '[x1min,x1max]' is invalid as x1min>x1max");
    if (nargs > 2)
    {
      e3 = bltin_get_ent (args[2]);
      if (ent_type (e3) != MATRIX_DENSE_REAL)
        rerror("hist2: Third argument 'nbins1' has to be scalar greater than 2");
      nbins1 = (int) class_double(e3);
      if (nbins1 < 2)
        rerror("hist2: Third argument 'nbins1' has to be scalar greater than 2");
      ecount = 2; // second and third argument used for 'x1'
    }
      // construct uniform range
    range1 = mdr_Create(1, nbins1+1);
    for (i=0; i<=nbins1; i++)
      MdrV0(range1,i) = x1min + (x1max-x1min)/((double) nbins1)*((double) i);
  }
  else
  {
    range1 = (MDR *) mdr_Float_BF(r);
    nbins1 = SIZE(range1) - 1;
    ecount = 1; // 'range1' and 'nbins1' obtained from second argument
  }

  // ----------------------------------
  // determining histogram for x2:
  // ----------------------------------
  // scalar: nbins2
  // range: two element vector, then [min, max]
  //        longer vector, then [r[0], r[1], r[2], ... r[n]]
  ecount ++;
  e4 = bltin_get_ent (args[ecount]);
  if (ent_type (e4) != MATRIX_DENSE_REAL)
    rerror("hist2: improper arguments for 'x2'");
  r2 = ent_data(e4);
  if (SIZE(r2)==1)
  {
      // nbins1, we do not expect any more arguments for x1
    nbins2 = mdrV0(r2,0);
    if (nbins2 < 2)
      rerror("hist2: argument as 'nbins2' has to be greater than 2");
      // figure out xmin, xmax
    x2min = rlab_dmin_vector(nrows, 1, &MdrV0(x,nrows));
    x2max = rlab_dmax_vector(nrows, 1, &MdrV0(x,nrows));
      // construct uniform range
    range2 = mdr_Create(1, nbins2+1);
    for (i=0; i<=nbins2; i++)
      MdrV0(range2,i) = x2min + (x2max-x2min)/((double) nbins2)*((double) i);
  }
  else if (SIZE(r2)==2)
  {
      // [xmin, xmax], we possibly may get 'nbins' as the third argument
    x2min = MdrV0(r2,0);
    x2max = MdrV0(r2,1);
    if (x2min > x2max)
      rerror("hist2: argument as '[x2min,x2max]' is invalid as x2min>x2max");
    if (nargs > ecount+1)
    {
      e5 = bltin_get_ent (args[ecount+1]);
      if (ent_type (e5) != MATRIX_DENSE_REAL)
        rerror("hist2: argument 'nbins2' has to be scalar greater than 2");
      nbins2 = (int) class_double(e5);
      if (nbins2 < 2)
        rerror("hist2: argument 'nbins2' has to be scalar greater than 2");
    }
      // construct uniform range
    range2 = mdr_Create(1, nbins2+1);
    for (i=0; i<=nbins2; i++)
      MdrV0(range2,i) = x2min + (x2max-x2min)/((double) nbins2)*((double) i);
  }
  else
  {
    range2 = (MDR *) mdr_Float_BF(r2);
    nbins2 = SIZE(range2) - 1;
  }

_do_hist2:

  // construct bin matrix
  bin = mdr_Create(nbins1, nbins2);
  mdr_Zero(bin);

  // create histogram
  hist2.nx = nbins2;
  hist2.ny = nbins1;
  hist2.xrange = MDRPTR(range2);
  hist2.yrange = MDRPTR(range1);
  hist2.bin    = MDRPTR(bin);
// 
  // go over each data point and update the bin count
  for (i=0; i<nrows; i++)
  {
    // ignore nan's
    if(isnand (Mdr0(x,i,0)) || isnand (Mdr0(x,i,1)))
      continue;

    gsl_histogram2d_increment( &hist2, Mdr0(x,i,1), Mdr0(x,i,0));
  }

  // cleanup
  hist2.nx = 0;
  hist2.ny = 0;
  hist2.xrange = 0;
  hist2.yrange = 0;
  hist2.bin    = 0;

  // clean up
  ent_Clean(e1);
  ent_Clean(e2);
  ent_Clean(e3);
  ent_Clean(e4);
  ent_Clean(e5);

  //
  // return the result as a list
  //
  Btree *bw = btree_Create ();
  // bins
  install (bw, RLAB_NAME_HIST2D_BIN, ent_Assign_Rlab_MDR(bin));
  // xrange
  install (bw, RLAB_NAME_HIST2D_XRANGE, ent_Assign_Rlab_MDR(range1));
  // yrange
  install (bw, RLAB_NAME_HIST2D_YRANGE, ent_Assign_Rlab_MDR(range2));

  return ent_Assign_Rlab_BTREE(bw);
}

#include "gnunet.c"

