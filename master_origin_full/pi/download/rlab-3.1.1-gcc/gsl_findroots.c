// Copyright (C) 2003-2008 Marijan Kostrun
//   part of rlab+4linux project on rlabplus.sourceforge.net
//
// GSL Science Library - root finding in multi-dimensions
//
//   This program is free software; you can redistribute it and/or modify
//   it under the terms of the GNU General Public License as published by
//   the Free Software Foundation; either version 2 of the License, or
//   (at your option) any later version.
//
//   This program is distributed in the hope that it will be useful,
//   but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//   GNU General Public License for more details.
//
//   You should have received a copy of the GNU General Public License
//   along with this program; if not, write to the Free Software
//   Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
//
//   See the file ./COPYING
//   ********************************************************************** */


// rlab headers, located in variable $RLAB_SDK
#include "rlab.h"
#include "ent.h"
#include "class.h"
#include "symbol.h"
#include "mem.h"
#include "mdr.h"
#include "mdc.h"
#include "mdrf1.h"
#include "mds.h"
#include "list.h"
#include "btree.h"
#include "bltin.h"
#include "util.h"
#include "mathl.h"
#include "function.h"
#include "lp.h"

#include "hompack.h"

// gsl headers
#include <gsl/gsl_mode.h>
#include <gsl/gsl_precision.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_machine.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_multiroots.h>

// standard headers
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <time.h>

// naming convention for the solver parameters
#include "rfileio.h"
#include "rlab_solver_parameters_names.h"

//
// findroots: link to gsl functions
//
static int
multiroot_gslrlab_f   (const gsl_vector *x, void *params, gsl_vector *f);
static int
multiroot_gslrlab_df  (const gsl_vector *x, void *params, gsl_matrix *J);
static int
multiroot_gslrlab_fdf (const gsl_vector *x, void *params, gsl_vector *f,
                       gsl_matrix * J)
{
  multiroot_gslrlab_f  (x, params, f);
  multiroot_gslrlab_df (x, params, J);
  return GSL_SUCCESS;
}

//
// findroots: other variables
//
static MDR *params=0;
static Ent *pent=0;

static int neq=0;

static MDR *my=0;
static Ent *my_ent=0;

static Ent *fname=0, *df_fname=0;

#undef  THIS_SOLVER
#define THIS_SOLVER "findroots"
Ent *
ent_gsl_findroots (int nargs, Datum args[])
{
  Ent *e2=0;
  ListNode *node;
  //
  // Load and Check arguments.
  //
  if (nargs < 3)
  {
    fprintf (stdout,
             "findroots: Finding the root of a function with and without\n");
    fprintf (stdout,
             "findroots: function's derivative in multidimensions.\n");
    fprintf (stdout,
             "findroots: Format:\n");
    fprintf (stdout,
             "findroots: (1) r = findroots(f,/p/, x0 /, options/),\n");
    fprintf (stdout,
             "findroots: (2) r = findroots(f,df,/p/, x0 /,options/),\n");
    fprintf (stdout,
             "findroots: where f = function(x /,p/), x0 is the initial point.\n");
    fprintf (stdout,
             "findroots: options=<<eabs;imethod;maxi;stdout>>, where in case (1)\n");
    fprintf (stdout,
             "findroots: imethod=0 for 'hybrids', 1 for 'hybrid', 2 for 'newton', or\n");
    fprintf (stdout,
             "findroots: 3 for 'broyden', and in case (2), imethod=0 for 'hybridsj',\n");
    fprintf (stdout,
             "findroots: 1 for 'hybridj', 2 for 'newton', and 3 for 'gnewton'.\n");
    fprintf (stdout,
             "findroots: In case (2) df=function(x/,p) must be specified, where \n");
    fprintf (stdout,
             "findroots: df=df/dx. 'stdout' is the name of the file/device where the run\n");
    fprintf (stdout,
             "findroots: time messages are posted, and 'maxi' is the maximum number of\n");
    fprintf (stdout,
             "findroots: iterations allowed to find the root with desired accuracy.\n");
    rerror ("requires at least 3 arguments");
  }

  //
  // f = f(x/,p/)
  //
  fname = bltin_get_ent(args[0]);
  if (!isfuncent(fname))
    rerror (THIS_SOLVER ": " RLAB_ERROR_ARG1_FUNC_VAR "\n");

  //
  // df = df(x/,p/) or an entity
  //
  e2 = bltin_get_ent( args[1] );
  if (!isfuncent(e2))
  {
    //
    // second argument is an entity, use non-jacobian methods
    //
    double abserr = 1e-6, timer, ddummy;
    int i, j, status, maxiter = 1000, imethod = 0, idummy;
    Ent *e3=0, *e4=0;
    MDR *ystart, *w;
    MDS *outstream;
    char outs[128] = { '\0' };
    FILE *fptr = NULL;
    time_t t1, t2;
    //
    // set up gsl solver
    //
    const gsl_multiroot_fsolver_type *T;
    gsl_multiroot_fsolver *s;

    //
    // parameter
    //
    pent = 0;
    if (ent_type (e2) != UNDEF)
      pent = ent_Copy (e2);

    //
    // Get ystart
    //
    e3 = bltin_get_ent (args[2]);
    if (ent_type (e3) != MATRIX_DENSE_REAL)
      rerror ("findroots: 'y0' must be a vector !");
    ystart = class_matrix_real (e3);
    // Extract neq from ystart
    if (MNC (ystart) != 1 && MNR (ystart) != 1)
      rerror ("findroots: 'y0' must be a vector !");
    neq = MNR (ystart) * MNC (ystart);

    //
    // options for the solver
    //
    if (nargs >= 4)
    {
      e4 = bltin_get_ent (args[3]);
      if (ent_type (e4) == BTREE)
      {
        // eabs
        node = btree_FindNode (ent_data (e4), RLAB_NAME_GSL_ROOTS_EABS);
        if (node != 0)
        {
          ddummy = class_double (var_ent (node));
          if (ddummy > 0.0)
            abserr = ddummy;
        }
        // method
        node = btree_FindNode (ent_data (e4), RLAB_NAME_GSL_ROOTS_METHOD);
        if (node != 0)
        {
          idummy = (int) class_double (var_ent (node));
          if (idummy == 0 || idummy == 1 || idummy == 2 || idummy == 3)
            imethod = idummy;
        }
        // standard output
        node = btree_FindNode (ent_data (e4), RLAB_NAME_GSL_ROOTS_STDOUT);
        if (node != 0)
        {
          if (ent_type(var_ent (node))== MATRIX_DENSE_STRING)
          {
            outstream = ent_data (var_ent (node));
            sprintf (outs, "%s", Mds1 (outstream, 1, 1));
          }
        }
        // max iterations
        node = btree_FindNode (ent_data (e4), RLAB_NAME_GSL_ROOTS_MAXI);
        if (node != 0)
        {
          idummy = (int) class_double (var_ent (node));
          if (idummy > 0)
            maxiter = idummy;
        }
      }
      if ( ent_type (e4) != UNDEF ) ent_Clean (e4);
    }
    //
    // Set up ENTITIES for user-function.
    // Inc the reference count once, for belonging to multiroot
    // array y[]:
    my = mdr_CreateEmpty (MNR(ystart), MNC(ystart));
    my_ent = ent_Assign_Rlab_MDR (my);
    ent_IncRef (my_ent);

    // parameter array p[]:

    //
    // findroots
    //
    gsl_vector *x = gsl_vector_alloc (neq);
    for (i = 0; i < neq; i++)
      gsl_vector_set (x, i, MdrV0 (ystart, i));
    // chose method
    gsl_multiroot_function f = { &multiroot_gslrlab_f, neq, NULL };
    switch (imethod)
    {
      case 0:
        T = gsl_multiroot_fsolver_hybrids;
        break;

      case 1:
        T = gsl_multiroot_fsolver_hybrid;
        break;

      case 2:
        T = gsl_multiroot_fsolver_dnewton;
        break;

      case 3:
        T = gsl_multiroot_fsolver_broyden;
        break;

      default:
        T = gsl_multiroot_fsolver_hybrids;
    }
    s = gsl_multiroot_fsolver_alloc (T, neq);
    gsl_multiroot_fsolver_set (s, &f, x);

    if(strlen (outs) > 1)
      fptr = fopen (outs, "a");

    if (!fptr)
    {
      // do not write any run-time messages
      i = 0;
      do
      {
        i++;
        status = gsl_multiroot_fsolver_iterate (s);
        if (status)
          break;
        status = gsl_multiroot_test_residual (s->f, abserr);
      }
      while ((status == GSL_CONTINUE) && (i < maxiter));
    }
    else
    {
      // use name from 'outs' as a device to which the run-time messages are written
      // wastes time if 'outs' points to non-existing stream
      // start the timer
      t1 = clock ();
      if(fptr)
        fprintf (fptr,
                 "RLaB: using gsl '%s' solver for multidimensional root finding.\n",
                 gsl_multiroot_fsolver_name (s));
      i = 0;
      do
      {
        i++;
        status = gsl_multiroot_fsolver_iterate (s);
        if (status)
          break;
        status = gsl_multiroot_test_residual (s->f, abserr);
        if(fptr)
        {
          fprintf (fptr, "iteration = %i", i);
          for (j = 0; j < neq; j++)
            fprintf (fptr, " : x[%i] = %g", j + 1, gsl_vector_get (s->x, j));
          for (j = 0; j < neq; j++)
            fprintf (fptr, " : f[%i] = %g", j + 1, gsl_vector_get (s->f, j));
          fprintf (fptr, "\n");
        }
      }
      while ((status == GSL_CONTINUE) && (i < maxiter));
    }
    w = mdr_Float_BF (ystart);
    for (i = 0; i < neq; i++)
    {
      MdrV0 (w, i) = gsl_vector_get (s->x, i);
    }
    if (fptr)
    {
      fprintf (fptr, "RLaB gsl '%s' solver", gsl_multiroot_fsolver_name (s));
      fprintf (fptr, " for multidimensional root finding reports: ");
      fprintf (fptr, " '%s' !\n", gsl_strerror (status));
      // check the time and close the output stream
      t2 = clock ();
      timer = (t2 - t1) / 1e6;
      fprintf (fptr, "RLaB: root finding lasted %g sec.\n", timer);
      fclose (fptr);
    }
    gsl_multiroot_fsolver_free (s);
    gsl_vector_free (x);

    // Clean Up
    MDPTR(my) = 0;
    ent_DecRef (my_ent);
    ent_Destroy (my_ent);

    if (pent)
      ent_Clean (pent);
    ent_Clean (e2);
    ent_Clean (e3);

    return ent_Assign_Rlab_MDR(w);
  }
  //
  // second parameter is a function, assume dfunc and use appropriate solvers
  //
  // HOMPACK
//   double arctol = -1;
//   double eps = 1e-16;
  double facsqnp1 = 1.0;
  double cursw = 10.0;
  int maxtry = 10;
  double arcre=1e-5, arcae=1e-5, ansre=1e-10, ansae=1e-10;

  // THE GSL
  double abserr = 1e-4, dummy, timer, ddummy;
  int i, j, status, idummy, imethod=0, louts=0, maxiter=1000;
  Ent *e3=0, *e4=0, *e5=0;
  MDR *ystart, *w;
  MDS *outstream;
  char outs[128] = { '\0' };
  FILE *fptr = NULL;

  time_t t1= clock (), t2;

  //
  // df = df(x/,p/)
  //
  df_fname = e2;

  //
  // parameter entity
  //
  e3 = bltin_get_ent (args[2]);
  pent = 0;
  if (ent_type (e3) != UNDEF)
    pent = ent_Copy (e3);

  //
  // x0
  //
  e4 = bltin_get_ent (args[3]);
  if (ent_type (e4) != MATRIX_DENSE_REAL)
    rerror ("findroots: Argument 'y0' must be real vector");
  ystart = class_matrix_real (e4);
  if (!ystart)
    rerror ("findroots: Argument 'y0' must be real vector");
  // Extract neq from ystart
  if (MNC (ystart) != 1 && MNR (ystart) != 1)
    rerror ("findroots: Argument 'y0' must be real vector");
  neq = MNR (ystart)*MNC(ystart);

  //
  // options for the solver
  //
  if (nargs >= 4)
  {
    e5 = bltin_get_ent (args[4]);
    if (ent_type (e5) == BTREE)
    {
      // gsl: eabs
      node = btree_FindNode (ent_data (e5), RLAB_NAME_GSL_ROOTS_EABS);
      if (node != 0)
      {
        ddummy = class_double (var_ent (node));
        if (ddummy > 0.0)
          abserr = ddummy;
      }
      // gsl: method
      node = btree_FindNode (ent_data (e5), RLAB_NAME_GSL_ROOTS_METHOD);
      if (node != 0)
      {
        idummy = (int) class_double (var_ent (node));
        if (idummy >= 0 && idummy <= 6)
          imethod = idummy;
      }
      // hompack: arcre
      node = btree_FindNode (ent_data (e5), RLAB_NAME_GSL_ROOTS_ARCRE);
      if (node != 0)
      {
        dummy =  class_double (var_ent (node));
        if (dummy == -1 || dummy > 0)
          arcre = dummy;
      }
      // hompack: arcae
      node = btree_FindNode (ent_data (e5), RLAB_NAME_GSL_ROOTS_ARCAE);
      if (node != 0)
      {
        dummy =  class_double (var_ent (node));
        if (dummy == -1 || dummy > 0)
          arcae = dummy;
      }
      // hompack: ansre
      node = btree_FindNode (ent_data (e5), RLAB_NAME_GSL_ROOTS_ANSRE);
      if (node != 0)
      {
        dummy =  class_double (var_ent (node));
        if (dummy == -1 || dummy > 0)
          ansre = dummy;
      }
      node = btree_FindNode (ent_data (e5), RLAB_NAME_GSL_ROOTS_ANSAE);
      if (node != 0)
      {
        dummy =  class_double (var_ent (node));
        if (dummy == -1 || dummy > 0)
          ansae = dummy;
      }
      // hompack: facsqnp1
      node = btree_FindNode (ent_data (e5), RLAB_NAME_GSL_ROOTS_FACA);
      if (node != 0)
      {
        dummy =  class_double (var_ent (node));
        if (dummy > 0)
          facsqnp1 = dummy;
      }
      // hompack: cursw
      node = btree_FindNode (ent_data (e5), RLAB_NAME_GSL_ROOTS_CURS);
      if (node != 0)
      {
        dummy =  class_double (var_ent (node));
        if (dummy > 0)
          cursw = dummy;
      }
      // hompack: max number of repetitions
      node = btree_FindNode (ent_data (e5), RLAB_NAME_GSL_ROOTS_IREP);
      if (node != 0)
      {
        idummy = (int) class_double (var_ent (node));
        if (idummy > 0)
          maxtry = idummy;
      }
      // gsl/hompack: standard output
      node = btree_FindNode (ent_data (e5), RLAB_NAME_GSL_ROOTS_STDOUT);
      if (node != 0)
      {
        if (ent_type(var_ent (node))== MATRIX_DENSE_STRING)
        {
          outstream = ent_data (var_ent (node));
          sprintf (outs, "%s", Mds1 (outstream, 1, 1));
          louts = strlen (outs);
        }
      }
      // gsl/hompack: max iterations
      node = btree_FindNode (ent_data (e5), RLAB_NAME_GSL_ROOTS_MAXI);
      if (node != 0)
      {
        idummy = (int) class_double (var_ent (node));
        if (idummy > 0)
          maxiter = idummy;
      }
    }
  }
  //
  // Set up ENTITIES for user-function.
  // Inc the reference count once, for belonging to multiroot
  // array y[]:
  my = mdr_CreateEmpty (MNR(ystart), MNC(ystart));
  my_ent = ent_Assign_Rlab_MDR (my);
  ent_IncRef (my_ent);

  // parameter array p[]:

  if (imethod <= 3)
  {
    //
    // The GSL root finder in multidimensions
    //
    const gsl_multiroot_fdfsolver_type *T;
    gsl_multiroot_fdfsolver *s;

    gsl_vector *x = gsl_vector_alloc (neq);
    for (i = 0; i < neq; i++)
      gsl_vector_set (x, i, MdrV0 (ystart, i));
    gsl_multiroot_function_fdf f = { &multiroot_gslrlab_f,
      &multiroot_gslrlab_df, &multiroot_gslrlab_fdf, neq, NULL
    };

    switch (imethod)
    {
      case 0:
        T = gsl_multiroot_fdfsolver_hybridsj;
        break;

      case 1:
        T = gsl_multiroot_fdfsolver_hybridj;
        break;

      case 2:
        T = gsl_multiroot_fdfsolver_newton;
        break;

      case 3:
        T = gsl_multiroot_fdfsolver_gnewton;
        break;

      default:
        T = gsl_multiroot_fdfsolver_hybridsj;
    }
    s = gsl_multiroot_fdfsolver_alloc (T, neq);
    gsl_multiroot_fdfsolver_set (s, &f, x);

    if (louts > 1)
      fptr = fopen (outs, "a");

    if (!fptr)
    {
      // do not write any run-time messages
      i = 0;
      do
      {
        i++;
        status = gsl_multiroot_fdfsolver_iterate (s);
        if (status)
          break;
        status = gsl_multiroot_test_residual (s->f, abserr);
      }
      while ((status == GSL_CONTINUE) && (i < maxiter));
    }
    else
    {
      // use name from 'outs' as a device to which the run-time messages are written
      // wastes time if 'outs' points to non-existing stream
      // start the timer
      t1 = clock ();
      if (fptr)
        fprintf (fptr,
                 "RLaB: using gsl '%s' solver for multidimensional root finding.\n",
                 gsl_multiroot_fdfsolver_name (s));
      i = 0;
      do
      {
        i++;
        status = gsl_multiroot_fdfsolver_iterate (s);
        if (status)
          break;
        status = gsl_multiroot_test_residual (s->f, abserr);
        if(fptr)
        {
          fprintf (fptr, "iteration = %i", i);
          for (j = 0; j < neq; j++)
            fprintf (fptr, " : x[%i] = %g", j + 1, gsl_vector_get (s->x, j));
          for (j = 0; j < neq; j++)
            fprintf (fptr, " : f[%i] = %g", j + 1, gsl_vector_get (s->f, j));
          fprintf (fptr, "\n");
        }
      }
      while ((status == GSL_CONTINUE) && (i < maxiter));
    }
    if (status == GSL_SUCCESS)
    {
      w = mdr_Float_BF (ystart);
      for (i = 0; i < neq; i++)
        MdrV0 (w, i) = gsl_vector_get (s->x, i);
    }
    else
      w = mdr_Create (0,0);

    if (fptr)
    {
      fprintf (fptr, "RLaB gsl '%s' solver", gsl_multiroot_fdfsolver_name (s));
      fprintf (fptr, " for multidimensional root finding reports: ");
      fprintf (fptr, " '%s' !\n", gsl_strerror (status));
      // check the time and close the output stream
      t2 = clock ();
      timer = (t2 - t1) / 1e6;
      fprintf (fptr, "RLaB: root finding lasted %g sec.\n", timer);
      fclose (fptr);
    }
    gsl_multiroot_fdfsolver_free (s);
    gsl_vector_free (x);
  }
  else
  {
    //
    // Hompack root finder in multidimensions
    //
    MDR * y = mdr_Create(neq+1,1);
    for (i=1;i<=neq;i++)
      MdrV0(y,i) = MdrV1(ystart,i);
    int iflag = -1;
    int itrace = 0;
    int nfe = 0;
    double arclen = 0.0;

    MDR * yp = mdr_Create(neq+1,1);
    MDR * ypold = mdr_Create(neq+1,1);
    MDR * a = mdr_Create(neq+1,1);

    // PDF
    MDR *qr=0, *wt=0, *phi=0, *p=0, *tz=0, *ipivot=0, *alpha=0;
    // PNF
    MDR * w1, * wp, * yold, * z0, * z1;
    // PQF
    MDR *qt, *r, *f0, *f1, *dz, *t, *ysav;

    // use name from 'outs' as a device to which the run-time messages are written
    // wastes time if 'outs' points to non-existing stream
    // start the timer
    if (louts > 1)
    {
      t1 = clock ();
      fptr = fopen (outs, "a");
      if (!fptr)
        louts = 0;
    }
    else
      louts = 0;

    switch (imethod)
    {
      case 4:
        //
        // PDF solver
        //
        // init locals
        qr = mdr_Create(neq, neq+1);
        wt = mdr_Create(neq+1,1);
        phi = mdr_Create(neq+1,16);
        p = mdr_Create(neq+1,1);
        tz = mdr_Create(neq+1,1);
        ipivot = mdi_Create(neq+1,1);
        alpha = mdr_Create(3*neq+3,1);

        // start
        j = 0;
        while (iflag == -1 || iflag == 5 || iflag == 2 || iflag == 3)
        {
          if (iflag == 5)
          {
            arcre = arcre / 8;
            ansre = ansre / 8;
            iflag = -1;
          }
          FIXPDF(&neq, MDRPTR(y), &iflag, &arcre, &ansre,
                  &itrace, MDRPTR(a), &neq, &nfe,
                                  &arclen, MDRPTR(yp), MDRPTR(ypold), MDRPTR(qr),
                 MDRPTR(alpha), MDRPTR(tz), MDIPTR(ipivot),
                 MDRPTR(wt), MDRPTR(phi), MDRPTR(p),
                 NULL, NULL,
                 outs, &louts, &maxiter, &facsqnp1, &cursw
                );
          j++;
          if (j>maxtry) break;
        }
        // clean
        mdr_Destroy (qr);
        mdr_Destroy (wt);
        mdr_Destroy (phi);
        mdr_Destroy (p);
        mdr_Destroy (tz);
        mdr_Destroy (ipivot);
        mdr_Destroy (alpha);

        // write the report from the solver
        if (fptr)
        {
          fprintf (fptr, "RLaB HOMPACK 'PDF' solver for multidimensional root finding");
          fprintf (fptr, " reports IFLAG = %i ", iflag);
          switch (iflag)
          {
            case 1:
              fprintf (fptr, "(success) !\n");
              break;
            case 2:
              fprintf (fptr, "(cannot meet the tolerances) !\n");
              break;
            case 3:
              fprintf (fptr, "(max no. of iterations reached) !\n");
              break;
            case 4:
              fprintf (fptr, "(jacobian matrix does not have full rank) !\n");
              break;
            case 5:
              fprintf (fptr, "(tolerances too great) !\n");
              break;
            case 6:
              fprintf (fptr, "(I - DF(X0) is nearly singular) !\n");
              break;
            default:
              fprintf (fptr, "(unknown error code) !\n");
          }
        }
        break; // done with PDF
      case 5:
        //
        // PNF solver
        //
        // init locals
        w1 = mdr_Create (neq+1,1);
        wp = mdr_Create (neq+1,1);
        yold = mdr_Create (neq+1,1);
        qr = mdr_Create (neq,neq+2);
        z0 = mdr_Create (neq+1,1);
        z1 = mdr_Create (neq+1,1);
        tz = mdr_Create(neq+1,1); mdr_Zero (tz);
        ipivot = mdi_Create(neq+1,1);
        alpha = mdr_Create(3*neq+3,1);
        double sspar_pnf[8]={0,0,0,0,0,0,0,0};
        // calculate
        j = 0;

        while (iflag == -1 || iflag == 2 || iflag == 3 || iflag == 6)
        {
          if (iflag == 6) iflag = -1;
          FIXPNF(&neq, MDRPTR(y), &iflag, &arcre, &arcae, &ansre, &ansae,
                 &itrace, MDRPTR(a), &nfe, &arclen,
                 MDRPTR(yp), MDRPTR(yold), MDRPTR(ypold), MDRPTR(qr),
                 MDRPTR(alpha), MDRPTR(tz), MDIPTR(ipivot),
                 MDRPTR(w1), MDRPTR(wp), MDRPTR(z0), MDRPTR(z1), sspar_pnf, NULL, NULL,
                 outs, &louts, &maxiter, &cursw );

          j++;
          if (j>maxtry) break;
        }

        // clean-up
        mdr_Destroy (qr);
        mdr_Destroy (w1);
        mdr_Destroy (wp);
        mdr_Destroy (yold);
        mdr_Destroy (z0);
        mdr_Destroy (z1);
        mdr_Destroy (tz);
        mdr_Destroy (ipivot);
        mdr_Destroy (alpha);

        // write the report from the solver
        if (fptr)
        {
          fprintf (fptr, "RLaB HOMPACK 'PNF' solver for multidimensional root finding");
          fprintf (fptr, " reports IFLAG = %i ", iflag);
          switch (iflag)
          {
            case 1:
              fprintf (fptr, "(success) !\n");
              break;
            case 2:
              fprintf (fptr, "(cannot meet the tolerances) !\n");
              break;
            case 3:
              fprintf (fptr, "(max no. of iterations reached) !\n");
              break;
            case 4:
              fprintf (fptr, "(jacobian matrix does not have full rank) !\n");
              break;
            case 5:
              fprintf (fptr, "(tolerances too great) !\n");
              break;
            case 6:
              fprintf (fptr, "(netwon iteration failed to converge) !\n");
              break;
            default:
              fprintf (fptr, "(unknown error code) !\n");
          }
        }
        break;
      case 6:
        //
        // PQF solver
        //
        // init locals
        yold = mdr_Create (neq+1,1);
        qt = mdr_Create (neq+1,neq+1);
        r = mdr_Create ((neq+1)*(neq+2)/2+1,1);
        f0 = mdr_Create (neq+1,1);
        f1 = mdr_Create (neq+1,1);
        z0 = mdr_Create (neq+1,1);
        dz = mdr_Create (neq+1,1);
        w1 = mdr_Create (neq+1,1);
        t = mdr_Create (neq+1,1);
        ysav = mdr_Create (neq+1,1);

        double sspar_pqf[4]={0,0,0,0};
        // calculate
        j = 0;
        while (iflag == -1 || iflag == 2 || iflag == 3 || iflag == 6)
        {
          if (iflag == 6) iflag = -1;
          FIXPQF(&neq, MDRPTR(y), &iflag, &arcre, &arcae, &ansre, &ansae,
                  &itrace, MDRPTR(a), &nfe, &arclen,
                  MDRPTR(yp), MDRPTR(yold), MDRPTR(ypold), MDRPTR(qt),
                  MDRPTR(r), MDRPTR(f0), MDRPTR(f1), MDRPTR(z0), MDRPTR(dz),
                  MDRPTR(w1), MDRPTR(t), MDRPTR(ysav), sspar_pqf, NULL, NULL,
                  outs, &louts, &maxiter);
          j++;
          if (j>maxtry) break;
        }
        // clean-up
        mdr_Destroy (yold);
        mdr_Destroy (qt);
        mdr_Destroy (r);
        mdr_Destroy (f0);
        mdr_Destroy (f1);
        mdr_Destroy (z0);
        mdr_Destroy (dz);
        mdr_Destroy (w1);
        mdr_Destroy (t);
        mdr_Destroy (ysav);

        // write the report from the solver
        if (fptr)
        {
          fprintf (fptr, "RLaB HOMPACK 'PQF' solver for multidimensional root finding");
          fprintf (fptr, " reports IFLAG = %i ", iflag);
          switch (iflag)
          {
            case 1:
              fprintf (fptr, "(success) !\n");
              break;
            case 2:
              fprintf (fptr, "(cannot meet the tolerances) !\n");
              break;
            case 3:
              fprintf (fptr, "(max no. of iterations reached) !\n");
              break;
            case 4:
              fprintf (fptr, "(jacobian matrix does not have full rank) !\n");
              break;
            case 5:
              fprintf (fptr, "(tolerances too great) !\n");
              break;
            case 6:
              fprintf (fptr, "(quasi-netwon iteration failed to converge) !\n");
              break;
            default:
              fprintf (fptr, "(unknown error code) !\n");
          }
        }
        break;
      default:
        rerror ("findroots: unknown solver!");
    }

    // check the time and close the output stream
    t2 = clock ();
    timer = (t2 - t1) / 1e6;
    if (fptr)
    {
      fprintf (fptr, "RLaB: root finding lasted %g sec.\n", timer);
      fclose (fptr);
    }

    mdr_Destroy (yp);
    mdr_Destroy (a);
    mdr_Destroy (ypold);

    if (iflag == 1)
    {
      w = mdr_Float_BF(ystart);
      for (i=1;i<=neq;i++)
        MdrV1(w,i) = MdrV0(y,i);
    }
    else
      w = mdr_Create(0,0);
    // clean up y
    mdr_Destroy (y);
  }

  // Clean Up for all
  if (pent)
    ent_Clean (pent);

  MDPTR(my) = 0;
  ent_DecRef (my_ent);
  ent_Destroy (my_ent);

  ent_Clean (e2);
  ent_Clean (e3);
  ent_Clean (e4);
  ent_Clean (e5);

  return ent_Assign_Rlab_MDR(w);
}

//
// curve tracking with hompack
//
#undef  THIS_SOLVER
#define THIS_SOLVER "crvtrack"
Ent *
ent_ctrack (int nargs, Datum args[])
{
  Ent *e3=0, *e4=0;
  MDR *y0=0, *y=0, *w=0;

  // HOMPACK
  double facsqnp1 = 1.0;
  double cursw = 20.0;
  int    maxtry = 10;
  double arcre=1e-5, arcae=1e-5, ansre=1e-10, ansae=1e-10;

  // CONTIN
  double stepsize = 1e-2, minstepsize = 1e-10;

  // All
  double timer, ddummy;
  int i, j, nvar, idummy, imethod=1, louts=0, maxiter=1000;

  char *outs = 0;
  FILE *fptr = NULL;
  FILE *rlab_stderr = (!RLAB_STDERR_DS) ? stderr : RLAB_STDERR_DS;

  time_t t1=clock (), t2;

  ListNode *node=0;

  //
  // Load and Check arguments.
  //
  if (nargs < 3 || nargs > 4)
  {
    fprintf (rlab_stderr,
             THIS_SOLVER ": Curve tracking in multidimensions using homothopy and contination.\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": Format:\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ":   x1 = crvtrack(rho,Drho,x0/,options/),\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": where rho = function(x,c), such that  rho(x0,0) = 0, the\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": solution 'x1' of which satisfy  rho(x1,1) = 0, and where\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": Drho = function(x,c) = [Drho/dc,Drho/dx] is its jacobian.\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": The list 'options' contains general parameters valid for\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": all methods, options=<<imethod;faca;curs;sold;irep;rlab_stderr>>,\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": and detail tied to specific solver: imethod = 4 for 'PDF',\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": 5 for 'PNF', 6 for 'PQF' solver; 'rlab_stderr' is the output stream\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": where time messages are posted, 'maxi' is a maximum number\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": of iterations, in a single repetition of 'irep' of them, allowed\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": to perform the curve tracking with desired accuracy. For other\n");
    fprintf (rlab_stderr,
             THIS_SOLVER ": parameters see the manual.\n");
    rerror ("requires at least 3 arguments");
  }
  //
  // Get rho
  //
  fname = bltin_get_ent(args[0]);
  if (!isfuncent(fname))
    rerror (THIS_SOLVER ": " RLAB_ERROR_ARG1_FUNC_VAR "\n");

  //
  // Get Drho
  //
  df_fname = bltin_get_ent(args[1]);
  if (!isfuncent(df_fname))
    rerror (THIS_SOLVER ": " RLAB_ERROR_ARG2_FUNC_VAR "\n");

  //
  // Get y0
  //
  e3 = bltin_get_ent (args[2]);
  if (ent_type (e3) != MATRIX_DENSE_REAL)
    rerror (THIS_SOLVER ": Argument 'x0 must be real vector");
  y0 = ent_data(e3);
  if (!EQVECT(y0))
    rerror (THIS_SOLVER ": Argument 'x0 must be real vector");

  //
  // Get options for the solver
  //
  if (nargs > 3)
  {
    e4 = bltin_get_ent (args[3]);
    if (ent_type (e4) == BTREE)
    {
      // all: method
      node = btree_FindNode (ent_data (e4), RLAB_NAME_GSL_ROOTS_METHOD);
      if (node != 0)
      {
        idummy = (int) class_double (var_ent (node));
        if (idummy == 1 || idummy == 4 || idummy == 5 ||idummy == 6)
          imethod = idummy;
      }
      // contin: initial stepsize
      node = btree_FindNode (ent_data (e4), RLAB_NAME_GSL_ROOTS_STEP);
      if (node != 0)
      {
        ddummy =  class_double (var_ent (node));
        if (ddummy > minstepsize)
          stepsize = ddummy;
        else
        {
          stepsize    = minstepsize;
          minstepsize = ddummy;
        }
      }
      // hompack: arcre
      node = btree_FindNode (ent_data (e4), RLAB_NAME_GSL_ROOTS_ARCRE);
      if (node != 0)
      {
        ddummy =  class_double (var_ent (node));
        if (ddummy == -1 || ddummy > 0)
          arcre = ddummy;
      }
      // hompack: arcae
      node = btree_FindNode (ent_data (e4), RLAB_NAME_GSL_ROOTS_ARCAE);
      if (node != 0)
      {
        ddummy =  class_double (var_ent (node));
        if (ddummy == -1 || ddummy > 0)
          arcae = ddummy;
      }
      // hompack: ansre; contin: rtol
      node = btree_FindNode (ent_data (e4), RLAB_NAME_GSL_ROOTS_ANSRE);
      if (node != 0)
      {
        ddummy =  class_double (var_ent (node));
        if (ddummy == -1 || ddummy > 0)
          ansre = ddummy;
      }
      // hompack: ansae; contin: atol
      node = btree_FindNode (ent_data (e4), RLAB_NAME_GSL_ROOTS_ANSAE);
      if (node != 0)
      {
        ddummy =  class_double (var_ent (node));
        if (ddummy == -1 || ddummy > 0)
          ansae = ddummy;
      }
      // hompack, contin: facsqnp1
      node = btree_FindNode (ent_data (e4), RLAB_NAME_GSL_ROOTS_FACA);
      if (node != 0)
      {
        ddummy =  class_double (var_ent (node));
        if (ddummy > 0)
          facsqnp1 = ddummy;
      }
      // hompack: cursw
      node = btree_FindNode (ent_data (e4), RLAB_NAME_GSL_ROOTS_CURS);
      if (node != 0)
      {
        ddummy =  class_double (var_ent (node));
        if (ddummy > 0)
          cursw = ddummy;
      }
      // hompack: max number of repetitions
      node = btree_FindNode (ent_data (e4), RLAB_NAME_GSL_ROOTS_IREP);
      if (node != 0)
      {
        idummy = (int) class_double (var_ent (node));
        if (idummy > 0)
          maxtry = idummy;
      }
      // gsl/hompack: standard output
      node = btree_FindNode (ent_data (e4), RLAB_NAME_GSL_ROOTS_STDOUT);
      if (node != 0)
      {
        if (ent_type(var_ent (node))== MATRIX_DENSE_STRING)
        {
          outs  = class_char_pointer(var_ent (node));
          louts = strlen (outs);
        }
      }
      // gsl/hompack: max iterations
      node = btree_FindNode (ent_data (e4), RLAB_NAME_GSL_ROOTS_MAXI);
      if (node != 0)
      {
        idummy = (int) class_double (var_ent (node));
        if (idummy > 0)
          maxiter = idummy;
      }
    }
  }

  // figure out the sizes
  neq  = SIZE(y0);
  nvar = neq + 1;

  // starting point of continuation:
  //  y = [0,y0] or [y0,0]
  if (MNR(y0)==1)
    y = mdr_Create(1,nvar);
  else
    y = mdr_Create(nvar,1);

  if (imethod==4 || imethod==5 || imethod==6)
  {
    //
    // hompack: y = [0,y0]
    //
    MdrV0(y,0) = 0;
    for (i=0;i<neq; i++)
      MdrV0(y,i+1) = mdrV0(y0,i);
  }
  else
  {
    //
    // contin: y = [y0,0]
    //
    MdrV0(y,neq) = 0;
    for (i=0; i<neq; i++)
      MdrV0(y,i) = mdrV0(y0,i);
  }


  //
  // Set up ENTITIES for user-function.
  // Inc the reference count once, for belonging to multiroot
  //  f = f(x,c)
  // passing x to function call
  my = mdr_CreateEmpty (MNR(y0), MNC(y0));
  my_ent = ent_Assign_Rlab_MDR (my);
  ent_IncRef (my_ent);

  // passing c to function call
  params = mdr_CreateEmpty(1,1);
  pent = ent_Assign_Rlab_MDR (params);
  ent_IncRef (pent);

  int iflag = -2;
  int itrace = 0;
  int nfe = 0;
  double arclen = 0.0;
  int liw, lrw;

  // CONTIN
  MDR *iwork, *rwork;
  // HOMPACK
  MDR * yp;
  MDR * ypold;
  MDR * a;
  // PDF
  MDR * qr, * wt, * phi, * p, * tz, * ipivot, * alpha;
  // PNF
  MDR * w1, * wp, * yold, * z0, * z1;
  // PQF
  MDR *qt, *r, *f0, *f1, *dz, *t, *ysav;

  // use name from 'outs' as a device to which the run-time messages are written
  // wastes time if 'outs' points to non-existing stream
  // start the timer
  if (louts)
  {
    fptr = fopen (outs, "a");
    t1 = clock ();
  }

  switch (imethod)
  {
    case 1:
      //
      // CONTIN solver
      //
      //
      // init locals
      //
      liw = 29 + neq + 1;
      iwork = mdi_Create (liw,1); mdr_Zero (iwork);
      lrw = 29 + 6 * nvar + nvar * nvar;
      rwork = mdr_Create (lrw, 1); mdr_Zero (rwork);
      // init iwork
      MdiV1(iwork,1) = 0;
      MdiV1(iwork,2) = nvar;  // start continuation in c (its initial value, c = 0)
      MdiV1(iwork,3) = 1;     // do continuation only in c
      MdiV1(iwork,5) = nvar;  // end continuation in c
      if (fptr)
      {
        MdiV1(iwork,7) = 3; // print most intermediate messages
        MdiV1(iwork,8) = 12; // stdout device
        fprintf (fptr, "RLaB: CONTIN solver for curve tracking using homothopy.\n");
        fprintf (fptr, "RLaB: CONTIN solver messages follow:\n");
        fflush (fptr);
      }
      else
        MdiV1(iwork,8) = 0;   // no stdout
      //
      MdiV1 (iwork,9)  = 0;   // pit1rhoj supplies the jacobian
      MdiV1 (iwork,10) = 0;   // assume starting point is correct
      MdiV1 (iwork,14) = liw;
      MdiV1 (iwork,16) = lrw;
      MdiV1 (iwork,17) = maxiter;
      // init rwork
      if (ansae == -1)
        MdrV1 (rwork,1)  = 1e-5;  // atol
      else
        MdrV1 (rwork,1)  = ansae;
      MdrV1 (rwork,2)  = ansre;         // rtol
      MdrV1 (rwork,3)  = minstepsize;   // minimum stepsize
      MdrV1 (rwork,4)  = sqrt(nvar) * facsqnp1; // max allowed step in iteration
      MdrV1 (rwork,5)  = stepsize;    // starting stepsize
      MdrV1 (rwork,6)  = 1.0;         // move in positive direction
      MdrV1 (rwork,7)  = 1.0;         // do continuation until c = 1

      while ( (MdiV1(iwork,1) != 3) )
      {
        // iterate
        PITCON(PIT1RHOJ, NULL, PIT1RHO,
               &iflag, NULL, MDIPTR(iwork), &liw, &nvar, MDRPTR(rwork),
              &lrw, MDRPTR(y), DENSLV,
              outs, &louts);
        if (iflag != 0)
          break;
      }

      // clean-up
      mdr_Destroy (iwork);
      mdr_Destroy (rwork);

      // write the report from the solver
      if (fptr)
      {
        fprintf (fptr, "RLaB CONTIN solver for curve tracking using homothopy");
        fprintf (fptr, " reports IFLAG = %i ", iflag);
        switch (iflag)
        {
          case 0:
            fprintf (fptr, "(success) !\n");
            break;
          case 3:
            fprintf (fptr, "(singular matrix encountered) !\n");
            break;
          case 4:
            fprintf (fptr, "(unsuccessful corrector iteration) !\n");
            break;
          case 5:
            fprintf (fptr, "(too many corrector steps) !\n");
            break;
          case 6:
            fprintf (fptr, "(null tangent vector) !\n");
            break;
          case 7:
            fprintf (fptr, "(root finder failed to converge) !\n");
            break;
          case 8:
            fprintf (fptr, "(limit point iteration took too many steps) !\n");
            break;
          default:
            fprintf (fptr, "(unknown error code) !\n");
        }
      }
      break;
    case 4:
      if (fptr)
      {
        fprintf (fptr, "RLaB: HOMPACK 'PDF' solver for curve tracking using homothopy.\n");
        fprintf (fptr, "RLaB: HOMPACK 'PDF' solver messages follow:\n");
        fflush (fptr);
      }
      //
      // PDF solver
      //
      // init locals
      yp = mdr_Create(neq+1,1);
      ypold = mdr_Create(neq+1,1);
      a = mdr_Copy( y0 );
      qr = mdr_Create(neq, neq+1);
      wt = mdr_Create(neq+1,1);
      phi = mdr_Create(neq+1,16);
      p = mdr_Create(neq+1,1);
      tz = mdr_Create(neq+1,neq);
      ipivot = mdi_Create(neq+1,1);
      alpha = mdr_Create(3*neq+3,1);
      // start
      j = 0;
      while (iflag == -2 || iflag == 5 || iflag == 2 || iflag == 3)
      {
        if (iflag == 5)
        {
          iflag = -2;
          ansre = ansre / 8;
          arcre = arcre / 8;
        }
        FIXPDF(&neq, MDRPTR(y), &iflag, &arcre, &ansre,
               &itrace, MDRPTR(a), &neq, &nfe,
               &arclen, MDRPTR(yp), MDRPTR(ypold), MDRPTR(qr),
               MDRPTR(alpha), MDRPTR(tz), MDIPTR(ipivot),
               MDRPTR(wt), MDRPTR(phi), MDRPTR(p),
               NULL, NULL,
               outs, &louts, &maxiter, &facsqnp1, &cursw
               );
        j++;
        if (j>maxtry) break;
      }
      // clean
      mdr_Destroy (yp);
      mdr_Destroy (a);
      mdr_Destroy (ypold);
      mdr_Destroy (qr);
      mdr_Destroy (wt);
      mdr_Destroy (phi);
      mdr_Destroy (p);
      mdr_Destroy (tz);
      mdr_Destroy (ipivot);
      mdr_Destroy (alpha);

      // write the report from the solver
      if (fptr)
      {
        fprintf (fptr, "RLaB HOMPACK 'PDF' solver for curve tracking using homothopy");
        fprintf (fptr, " reports IFLAG = %i ", iflag);
        switch (iflag)
        {
          case 1:
            fprintf (fptr, "(success) !\n");
            break;
          case 2:
            fprintf (fptr, "(cannot meet the tolerances) !\n");
            break;
          case 3:
            fprintf (fptr, "(max no. of iterations reached) !\n");
            break;
          case 4:
            fprintf (fptr, "(jacobian matrix does not have full rank) !\n");
            break;
          case 5:
            fprintf (fptr, "(tolerances too great) !\n");
            break;
          case 6:
            fprintf (fptr, "(I - DF(X0) is nearly singular) !\n");
            break;
          default:
            fprintf (fptr, "(unknown error code) !\n");
        }
      }
      break; // done with PDF
    case 5:
      if (fptr)
      {
        fprintf (fptr, "RLaB: HOMPACK 'PNF' solver for curve tracking using homothopy.\n");
        fprintf (fptr, "RLaB: HOMPACK 'PNF' solver messages follow:\n");
        fflush (fptr);
      }
      //
      // PNF solver
      //
      // init locals
      yp = mdr_Create(neq+1,1);
      ypold = mdr_Create(neq+1,1);
      a = mdr_Create(neq+1,1);
      w1 = mdr_Create (neq+1,1);
      wp = mdr_Create (neq+1,1);
      yold = mdr_Create (neq+1,1);
      qr = mdr_Create (neq,neq+2);
      z0 = mdr_Create (neq+1,1);
      z1 = mdr_Create (neq+1,1);
      tz = mdr_Create(neq+1,neq); mdr_Zero (tz);
      ipivot = mdi_Create(neq+1,1);
      alpha = mdr_Create(3*neq+3,1);
      double sspar_pnf[8]={0,0,0,0,0,0,0,0};
      // calculate
      j = 0;
      while (iflag == -2 || iflag == 2 || iflag == 3 || iflag == 5 || iflag == 6)
      {
        if (iflag == 5)
        {
          // the error tolerances too lenient
          arcre = arcre / 8;
          arcae = arcae / 8;
          ansre = ansre / 8;
          ansae = ansae / 8;
          iflag = -2;
        }
        if (iflag == 6)
        {
          // the error tolerances too stringent
          arcre = arcre * 4;
          arcae = arcae * 4;
          ansre = ansre * 4;
          ansae = ansae * 4;
          iflag = -2;
        }
        FIXPNF(&neq, MDRPTR(y), &iflag, &arcre, &arcae, &ansre, &ansae,
               &itrace, MDRPTR(a), &nfe, &arclen,
               MDRPTR(yp), MDRPTR(yold), MDRPTR(ypold), MDRPTR(qr),
               MDRPTR(alpha), MDRPTR(tz), MDIPTR(ipivot),
               MDRPTR(w1), MDRPTR(wp), MDRPTR(z0), MDRPTR(z1), sspar_pnf, NULL, NULL,
               outs, &louts, &maxiter, &cursw );
        j++;
        if (j>maxtry) break;
      }
      // clean-up
      mdr_Destroy (yp);
      mdr_Destroy (a);
      mdr_Destroy (ypold);
      mdr_Destroy (qr);
      mdr_Destroy (w1);
      mdr_Destroy (wp);
      mdr_Destroy (yold);
      mdr_Destroy (z0);
      mdr_Destroy (z1);
      mdr_Destroy (tz);
      mdr_Destroy (ipivot);
      mdr_Destroy (alpha);
      // write the report from the solver
      if (fptr)
      {
        fprintf (fptr, "RLaB HOMPACK 'PNF' solver for curve tracking using homothopy");
        fprintf (fptr, " reports IFLAG = %i ", iflag);
        switch (iflag)
        {
          case 1:
            fprintf (fptr, "(success) !\n");
            break;
          case 2:
            fprintf (fptr, "(cannot meet the tolerances) !\n");
            break;
          case 3:
            fprintf (fptr, "(max no. of iterations reached) !\n");
            break;
          case 4:
            fprintf (fptr, "(jacobian matrix does not have full rank) !\n");
            break;
          case 5:
            fprintf (fptr, "(tolerances too great) !\n");
            break;
          case 6:
            fprintf (fptr, "(netwon iteration failed to converge) !\n");
            break;
          default:
            fprintf (fptr, "(unknown error code) !\n");
        }
      }
      break;
    case 6:
      if (fptr)
      {
        fprintf (fptr, "RLaB: HOMPACK 'PQF' solver for curve tracking using homothopy.\n");
        fprintf (fptr, "RLaB: HOMPACK 'PQF' solver messages follow:\n");
        fflush (fptr);
      }

      //
      // PQF solver
      //
      // init locals
      yp = mdr_Create(neq+1,1);
      ypold = mdr_Create(neq+1,1);
      a = mdr_Create(neq+1,1);
      yold = mdr_Create (neq+1,1);
      qt = mdr_Create (neq+1,neq+1);
      r = mdr_Create ((neq+1)*(neq+2)/2+1,1);
      f0 = mdr_Create (neq+1,1);
      f1 = mdr_Create (neq+1,1);
      z0 = mdr_Create (neq+1,1);
      dz = mdr_Create (neq+1,1);
      w1 = mdr_Create (neq+1,1);
      t = mdr_Create (neq+1,1);
      ysav = mdr_Create (neq+1,1);
      double sspar_pqf[4]={0,0,0,0};
      // calculate
      j = 0;
      while (iflag == -2 || iflag == 2 || iflag == 3 || iflag == 5 || iflag == 6)
      {
        if (iflag == 5)
        {
          // the error tolerances too lenient
          arcre = arcre / 8;
          arcae = arcae / 8;
          ansre = ansre / 8;
          ansae = ansae / 8;
          iflag = -2;
        }
        if (iflag == 6)
        {
          // the error tolerances too stringent
          arcre = arcre * 4;
          arcae = arcae * 4;
          ansre = ansre * 4;
          ansae = ansae * 4;
          iflag = -2;
        }
        FIXPQF(&neq, MDRPTR(y), &iflag, &arcre, &arcae, &ansre, &ansae,
               &itrace, MDRPTR(a), &nfe, &arclen,
               MDRPTR(yp), MDRPTR(yold), MDRPTR(ypold), MDRPTR(qt),
               MDRPTR(r), MDRPTR(f0), MDRPTR(f1), MDRPTR(z0), MDRPTR(dz),
               MDRPTR(w1), MDRPTR(t), MDRPTR(ysav), sspar_pqf, NULL, NULL,
               outs, &louts, &maxiter);
        j++;
        if (j>maxtry) break;
      }
      // clean-up
      mdr_Destroy (yp);
      mdr_Destroy (a);
      mdr_Destroy (ypold);
      mdr_Destroy (yold);
      mdr_Destroy (qt);
      mdr_Destroy (r);
      mdr_Destroy (f0);
      mdr_Destroy (f1);
      mdr_Destroy (z0);
      mdr_Destroy (dz);
      mdr_Destroy (w1);
      mdr_Destroy (t);
      mdr_Destroy (ysav);
      // write the report from the solver
      if (fptr)
      {
        fprintf (fptr, "RLaB HOMPACK 'PQF' solver for curve tracking using homothopy");
        fprintf (fptr, " reports IFLAG = %i ", iflag);
        switch (iflag)
        {
          case 1:
            fprintf (fptr, "(success) !\n");
            break;
          case 2:
            fprintf (fptr, "(cannot meet the tolerances) !\n");
            break;
          case 3:
            fprintf (fptr, "(max no. of iterations reached) !\n");
            break;
          case 4:
            fprintf (fptr, "(jacobian matrix does not have full rank) !\n");
            break;
          case 5:
            fprintf (fptr, "(tolerances too great) !\n");
            break;
          case 6:
            fprintf (fptr, "(quasi-netwon iteration failed to converge) !\n");
            break;
          default:
            fprintf (fptr, "(unknown error code) !\n");
        }
      }
      break;
    default:
      rerror (THIS_SOLVER ": unknown solver!");
  }

  if(fptr)
  {
    // check the time and close the output stream
    t2 = clock ();
    timer = (t2 - t1) / 1e6;
    fprintf (fptr, "RLaB: root finding lasted %g sec.\n", timer);
    fclose (fptr);
  }

  if (iflag == 1 && imethod >= 4 && imethod <= 6)
  {
    //
    // hompack
    //
    w = mdr_Float_BF (y0);
    for (i=1;i<=neq;i++)
      MdrV1(w,i) = MdrV0(y,i);
  }
  else if (iflag == 0 && imethod == 1)
  {
    //
    // contin
    //
    w = mdr_Float_BF (y0);
    for (i=0; i < neq; i++)
      MdrV0(w,i) = MdrV0(y,i);
  }
  else
    w = mdr_Create(0,0);

  // clean up y
  mdr_Destroy (y);

  // Clean Up for all
  // c
  MDPTR(params) = 0;
  ent_DecRef (pent);
  ent_Destroy (pent);

  // x
  MDPTR(my) = 0;
  ent_DecRef (my_ent);
  ent_Destroy (my_ent);

  ent_Clean (e3);
  ent_Clean (e4);
  ent_Clean (fname);
  ent_Clean (df_fname);
  
  return ent_Assign_Rlab_MDR(w);
}

//
// The GSL functions
//

static int
multiroot_gslrlab_f (const gsl_vector * x, void *params, gsl_vector * f)
{
  int i;
  Ent *rent = 0;
  MDR *retm = 0;

  // Put 'y' into rks_args, 'p' is already there
  //for (i = 0; i < neq; i++)
  //  MdrV0 (my, i) = gsl_vector_get (x, i);
  MDPTR(my) = (void *) x->data;

  if (pent)
    rent = ent_call_rlab_script_2args(fname, my_ent, pent);
  else
    rent = ent_call_rlab_script_1arg (fname, my_ent);

  if (ent_type(rent)!=MATRIX_DENSE_REAL)
    rerror (THIS_SOLVER ": " RLAB_ERROR_RHS_FUNC_MUST_RETURN_MDR);

  retm = ent_data (rent);

  if (SIZE(retm) != neq)
    rerror (THIS_SOLVER ": " RLAB_ERROR_RHS_FUNC_INCORRECT_DIM);

  for (i = 0; i < neq; i++)
    gsl_vector_set (f, i, MdrV0 (retm, i));

  ent_Clean (rent);
  return GSL_SUCCESS;
}

//
// jacobian function
//
static int
multiroot_gslrlab_df (const gsl_vector * x, void *params, gsl_matrix * J)
{
  int i, j;
  Ent *rent = 0;
  MDR *retm = 0;

  // Put 'y' into rks_args, 'p' is already there
  //for (i = 0; i < neq; i++)
  //  MdrV0 (my, i) = gsl_vector_get (x, i);
  MDPTR(my) = (void *) x->data;

  if (pent)
    rent = ent_call_rlab_script_2args(df_fname, my_ent, pent);
  else
    rent = ent_call_rlab_script_1arg (df_fname, my_ent);

  if (ent_type(rent)!=MATRIX_DENSE_REAL)
    rerror (THIS_SOLVER ": " RLAB_ERROR_RHS_FUNC_MUST_RETURN_MDR);

  retm = ent_data (rent);

  if (SIZE(retm) != neq * neq)
    rerror (THIS_SOLVER ": " RLAB_ERROR_JACOBIAN_FUNC_INCORRECT_DIM);

  // initialize Jacobian from the first neq x neq part of the
  // matrix retm
  for (i = 0; i < neq; i++)
    for (j = 0; j < neq; j++)
      gsl_matrix_set (J, i, j, Mdr0 (retm, i, j));

  // Try to clean up if possible.
  ent_Clean (rent);
  return GSL_SUCCESS;
}

//
// HOMPACK functions
//

int
HOM1F (double * x, double * v)
{
  int i;
  Ent *rent = 0;
  MDR *retm = 0;

  // Put 'y' into rks_args, 'p' is already there
  //for (i = 0; i < neq; i++)
  //  MdrV0 (my, i) = gsl_vector_get (x, i);
  MDPTR(my) = (void *) x;

  if (pent)
    rent = ent_call_rlab_script_2args(fname, my_ent, pent);
  else
    rent = ent_call_rlab_script_1arg (fname, my_ent);

  if (ent_type(rent)!=MATRIX_DENSE_REAL)
    rerror (THIS_SOLVER ": " RLAB_ERROR_RHS_FUNC_MUST_RETURN_MDR);

  retm = ent_data (rent);

  if (SIZE(retm) != neq)
    rerror (THIS_SOLVER ": " RLAB_ERROR_RHS_FUNC_INCORRECT_DIM);

  for (i = 0; i < neq; i++)
    v[i] = MdrV0 (retm, i);

  ent_Clean (rent);
  return 1;
}

//
// jacobian function
//
int
HOM1FJAC (double * x, double * v, int *k)
{
  int i;
  static MDR *retm = 0;

  if (*k == 1)
  {
    Ent *rent = 0;

    //
    // call the rlab script the first time. save the jacobian matrix
    // and deliver its columns as requested.
    //

    // Put 'y' into rks_args, 'p' is already there
    MDPTR(my) = (void *) x;

    if (pent)
      rent = ent_call_rlab_script_2args(df_fname, my_ent, pent);
    else
      rent = ent_call_rlab_script_1arg (df_fname, my_ent);

    if (ent_type(rent)!=MATRIX_DENSE_REAL)
      rerror (THIS_SOLVER ": " RLAB_ERROR_RHS_FUNC_MUST_RETURN_MDR);

    retm = mdr_Copy(ent_data (rent));

    if (SIZE(retm) != neq * neq)
      rerror (THIS_SOLVER ": " RLAB_ERROR_JACOBIAN_FUNC_INCORRECT_DIM);

    ent_Clean (rent);
  }

  //
  // return the requested column of the Jacobian
  //
  for (i = 0; i < neq; i++)
    v[i] = Mdr0 (retm, i, *k - 1);

  // clean up after the last column
  if (*k == neq)
    mdr_Destroy (retm);

  return 1;
}

//
// hompack: curve tracking routine  rho = function(x,lambda)
//
int
HOM1RHO (double * a, double *lambda, double *x, double * v,
        double *par, int * ipar)
{
  int i;
  Ent *rent = 0;
  MDR *retm = 0;

  MDPTR(my)     = (void *) x;
  MDPTR(params) = (void *) lambda;

  rent = ent_call_rlab_script_2args(fname, my_ent, pent);

  if (ent_type(rent)!=MATRIX_DENSE_REAL)
    rerror (THIS_SOLVER ": " RLAB_ERROR_RHS_FUNC_MUST_RETURN_MDR);

  retm = ent_data (rent);

  if (SIZE(retm) != neq)
    rerror (THIS_SOLVER ": " RLAB_ERROR_RHS_FUNC_INCORRECT_DIM);

  for (i = 0; i < neq; i++)
    v[i] = MdrV0 (retm, i);

  ent_Clean (rent);
  return 1;
}

int
HOM1RHOJ (double * a, double *lambda, double *x, double * v, int * k,
             double *par, int * ipar)
{
  int i;
  static MDR *retm = 0;

  if (*k == 1)
  {
    Ent *rent = 0;

    //
    // call the rlab script the first time. save the jacobian matrix
    // and deliver its columns as requested.
    //

    MDPTR(my) = (void *) x;
    MDPTR(params) = (void *) lambda;

    rent = ent_call_rlab_script_2args(df_fname, my_ent, pent);

    if (ent_type(rent)!=MATRIX_DENSE_REAL)
      rerror (THIS_SOLVER ": " RLAB_ERROR_RHS_FUNC_MUST_RETURN_MDR);

    retm = mdr_Copy(ent_data (rent));

    if (SIZE(retm) != neq * (neq+1))
      rerror (THIS_SOLVER ": " RLAB_ERROR_JACOBIAN_FUNC_INCORRECT_DIM);

    ent_Clean (rent);
  }

  //
  // deliver the requested column of the jacobian
  //
  if(*k == 1)
    for (i = 0; i < neq; i++)
    {
      // retm = [df/dx,df/dc]
      // rhojac requires [df/dc,df/dx]
      v[i] = Mdr0 (retm, i, neq);
    }
  else
    for (i = 0; i < neq; i++)
      v[i] = Mdr0 (retm, i, *k-2);

  // Try to clean up if possible.
  if (*k == neq+1)
    mdr_Destroy (retm);

  return 1;
}

int
PIT1RHO (int *nvar, double * dummy, int * idummy, double *x, double *v,
         int * ipar)
{
  int i;
  Ent *rent = 0;
  MDR *retm = 0;

  MDPTR(my) = (void *) x;
  MDPTR(params) = (void *) &x[neq];

  rent = ent_call_rlab_script_2args(fname, my_ent, pent);

  if (ent_type(rent)!=MATRIX_DENSE_REAL)
    rerror (THIS_SOLVER ": " RLAB_ERROR_RHS_FUNC_MUST_RETURN_MDR);

  retm = ent_data (rent);

  if (SIZE(retm) != neq)
    rerror (THIS_SOLVER ": " RLAB_ERROR_RHS_FUNC_INCORRECT_DIM);

  for (i=0; i < neq; i++)
    v[i] = MdrV0 (retm, i);

  ent_Clean (rent);
  return 1;
}

int
PIT1RHOJ (int *nvar, double *dummy, int *idummy, double *x, double * v,
          int * info)
{
  int i, j;
  Ent *rent = 0;
  MDR *retm = 0;

  MDPTR(my)     = (void *) x;
  MDPTR(params) = (void *) &x[neq];

  rent = ent_call_rlab_script_2args(df_fname, my_ent, pent);

  if (ent_type(rent)!=MATRIX_DENSE_REAL)
    rerror (THIS_SOLVER ": " RLAB_ERROR_JACOBIAN_FUNC_MUST_RETURN_MDR);

  retm = ent_data (rent);

  if (SIZE(retm) != neq*(neq+1) )
    rerror (THIS_SOLVER ": " RLAB_ERROR_JACOBIAN_FUNC_INCORRECT_DIM);

  for (i = 0; i < neq; i++)
    for (j = 0; j <= neq; j++)
      v[j * (neq+1) + i] = Mdr0 (retm, i, j);

  *info = 0;
  ent_Clean (rent);
  return 1;
}
