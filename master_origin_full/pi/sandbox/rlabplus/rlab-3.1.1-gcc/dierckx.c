// Copyright (C) 2000-2013 Marijan Kostrun
//   part of rlab+4linux project on rlabplus.sourceforge.net
//
// b-spline smoothing using dierckx's fitpack
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

// rlab headers
#include "complex.h"
#include "rlab.h"
#include "ent.h"
#include "class.h"
#include "mem.h"
#include "bltin.h"
#include "util.h"
#include "mdr.h"
#include "mdrf1.h"
#include "mdrf2.h"
#include "btree.h"
#include "symbol.h"
#include "mathl.h"

#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <math.h>

// gsl headers
// shared object
#include <gsl/gsl_mode.h>
#include <gsl/gsl_precision.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_machine.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_combination.h>

#include "dierckx.h"

// naming convention for the solver parameters
#include "rlab_solver_parameters_names.h"

#undef THIS_SOLVER
#define THIS_SOLVER "bsplinefit2"
Ent *
ent_bspline_fit2 (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *e3=0, *e4=0, *e5=0;
  MDR *x=0, *y=0, *z=0, *w=0, *rettx=0, *retty=0, *retc=0;
  MDR *tx=0, *ty=0, *c=0, *wrk1=0, *wrk2=0, *iwrk=0, *ckx=0, *cky=0;

  ListNode *node, *node2;

  int i, ny, nx, nz=0,nrz=0,ncz=0, idegx=3,idegy=3,iopt=0,nknotx=0,nknoty=0;
  int nxest=0,nyest=0,nmax;
  int lwrk1, lwrk2, liwrk, maxit=40, idummy, ier=0;
  int mynxest=0, mynyest=0;

  double s=0,xe=0,xb=0,ye=0,yb=0, fp, tol=0.001, ddummy, deps=1e-6;

  if (nargs != 4 && nargs != 5)
  {
    rerror (THIS_SOLVER ": Requires four or five arguments !");
  }

  // Z
  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) != MATRIX_DENSE_REAL && ent_type (e1) != BTREE)
    rerror (THIS_SOLVER ": " RLAB_ERROR_ARG1_STAT);
  if (ent_type (e1) == MATRIX_DENSE_REAL)
    z = class_matrix_real (e1);
  else
  {
    // val
    node = btree_FindNode (ent_data (e1), "val");
    if (!node)
      rerror (THIS_SOLVER ": " RLAB_ERROR_ARG1_STAT);
    z = class_matrix_real (var_ent (node));
    // wgt: may be absent
    node = btree_FindNode (ent_data (e1), "wgt");
    if (node)
      w = mdr_Float_BF (var_ent (node));
  }
  if (!z)
    rerror (THIS_SOLVER ": " RLAB_ERROR_ARG1_STAT);

  nz  = SIZE(z);
  if (nz<1)
    rerror (THIS_SOLVER ": " RLAB_ERROR_ARG1_STAT);
  nrz = MNR(z);
  ncz = MNC(z);

  if (!w)
  {
    // create weights
    w = mdr_Create(nrz,ncz);
    mdr_Ones (w);
  }
  if (SIZE(w) != nz)
    rerror (THIS_SOLVER ": " RLAB_ERROR_ARG1_MISMATCH_VAL_WGT);

  // X
  e2 = bltin_get_ent (args[1]);
  if (ent_type (e2) != MATRIX_DENSE_REAL && ent_type (e2) != BTREE)
    rerror (THIS_SOLVER ": " RLAB_ERROR_ARG2_STAT);
  if (ent_type (e2) == MATRIX_DENSE_REAL)
    x = class_matrix_real (e2);
  else
  {
    // val
    node = btree_FindNode (ent_data (e2), "val");
    if (!node)
      rerror (THIS_SOLVER ": " RLAB_ERROR_ARG2_STAT);
    x = class_matrix_real (var_ent (node));
  }
  if (!x)
    rerror (THIS_SOLVER ": " RLAB_ERROR_ARG2_STAT);
  nx = SIZE (x);
  if (nx<1)
    rerror (THIS_SOLVER ": " RLAB_ERROR_ARG2_STAT);

  // Y
  e3 = bltin_get_ent (args[2]);
  if (ent_type (e3) != MATRIX_DENSE_REAL && ent_type (e3) != BTREE)
    rerror (THIS_SOLVER ": " RLAB_ERROR_ARG3_STAT);
  if (ent_type (e3) == MATRIX_DENSE_REAL)
    y = class_matrix_real (e3);
  else
  {
    // val
    node = btree_FindNode (ent_data (e3), "val");
    if (!node)
      rerror (THIS_SOLVER ": " RLAB_ERROR_ARG3_STAT);
    y = class_matrix_real (var_ent (node));
  }
  if (!y)
    rerror (THIS_SOLVER ": " RLAB_ERROR_ARG3_STAT);
  ny = SIZE (y);
  if (ny<1)
    rerror (THIS_SOLVER ": " RLAB_ERROR_ARG3_STAT);

  //
  // smoothing parameter or knots for ls b-spline
  //
  e4 = bltin_get_ent (args[3]);
  if (ent_type (e4) == MATRIX_DENSE_REAL)
  {
    s = class_double (e4);
  }
  else if (ent_type (e4) == BTREE)
  {
    // x
    node = btree_FindNode (ent_data (e4), "1");
    if (!node)
      rerror (THIS_SOLVER ": ARG4 missing entry '1'");
    ckx = class_matrix_real (var_ent (node));

    // y
    node = btree_FindNode (ent_data (e4), "2");
    if (!node)
      rerror (THIS_SOLVER ": ARG4 missing entry '2'");
    cky = class_matrix_real (var_ent (node));


    // ckx and cky contain knots of the spline
    iopt = -1;    // least-squares fit spline

    nknotx = MNR(ckx) * MNC (ckx) + 2 * idegx;
    nxest  = nknotx;
    if (nknotx < 2*idegx+2)
      rerror (THIS_SOLVER ": ARG4 entry '1' (knot array) is improper size");

    nknoty = MNR(cky) * MNC (cky) + 2 * idegy;
    nyest  = nknoty;
    if (nknoty < 2*idegy+2)
      rerror (THIS_SOLVER ": ARG4 entry '2' (knot array) is improper size");

  }
  else
    rerror (THIS_SOLVER ": ARG4 is not list <<1;2>>");

  //
  // options
  //
  if (nargs > 4)
  {
    e5 = bltin_get_ent (args[4]);
    if (ent_type (e5) == BTREE)
    {
      // max_knots = <<1;2>>
      node = btree_FindNode (ent_data(e5), RLAB_NAME_DIERCKX_MAXKNOTS);
      if (node != 0)
      {
        if (ent_type(var_ent (node)) == BTREE)
        {
          // max_knots.1
          node2 = btree_FindNode (ent_data(var_ent(node)), "1");
          if (node2 != 0)
            mynxest = (int) class_double (var_ent (node2));
          if (mynxest<0)
            mynxest = 0;
          // max_knots.2
          node2 = btree_FindNode (ent_data(var_ent(node)), "2");
          if (node2 != 0)
            mynyest = (int) class_double (var_ent (node2));
          if (mynyest<0)
            mynyest = 0;
        }
      }
      // degree = <<x;y>>
      node = btree_FindNode (ent_data(e5), RLAB_NAME_DIERCKX_DEGREE);
      if (node != 0)
      {
        if (ent_type(var_ent (node)) == BTREE)
        {
          // degree.1
          node2 = btree_FindNode (ent_data(var_ent(node)), "1");
          if (node2 != 0)
          {
            idegx = (int) class_double (var_ent (node2));
            if (idegx != 1 && idegx != 3 && idegx != 5)
              idegx = 3;
          }
          // degree.2
          node2 = btree_FindNode (ent_data(var_ent(node)), "2");
          if (node2 != 0)
          {
            idegy = (int) class_double (var_ent (node2));
            if (idegy != 1 && idegy != 3 && idegy != 5)
              idegy = 3;
          }
        }
      }
      // range = <<1;2>>
      node = btree_FindNode (ent_data (e5), RLAB_NAME_DIERCKX_RANGE);
      if (node != 0)
      {
        if (ent_type(var_ent (node)) == BTREE)
        {
          // range.1
          node2 = btree_FindNode (ent_data(var_ent(node)), "1");
          if (node2 != 0)
          {
            if (ent_type(var_ent(node2)) == MATRIX_DENSE_REAL)
            {
              MDR *range1 = class_matrix_real (var_ent (node2));
              if (MNR(range1)*MNC(range1)==2)
              {
                xb = MdrV0(range1,0);
                xe = MdrV0(range1,1);
              }
            }
          }
          // degree.y
          node2 = btree_FindNode (ent_data(var_ent(node)), "2");
          if (node2 != 0)
          {
            if (ent_type(var_ent(node2)) == MATRIX_DENSE_REAL)
            {
              MDR *range2 = class_matrix_real (var_ent (node2));
              if (MNR(range2)*MNC(range2)==2)
              {
                yb = MdrV0(range2,0);
                ye = MdrV0(range2,1);
              }
            }
          }
        }
      }

      // tol
      node = btree_FindNode (ent_data(e5), RLAB_NAME_DIERCKX_TOL);
      if (node != 0)
      {
        ddummy = class_double (var_ent (node));
        if (ddummy > 0)
          tol = ddummy;
      }
      // maxi
      node = btree_FindNode (ent_data(e5), RLAB_NAME_DIERCKX_MAXITER);
      if (node != 0)
      {
        idummy = (int) class_double (var_ent (node));
        if (idummy > 0)
          maxit = idummy;
      }
      // rank_threshold
      node = btree_FindNode (ent_data(e5), RLAB_NAME_DIERCKX_RANKTHR);
      if (node != 0)
      {
        ddummy = class_double (var_ent (node));
        if (ddummy > 0)
          deps = ddummy;
      }
    }
  }


  // user can specify endpoints of spline fit
  if (xe==0 && xb==0 && xe==xb)
  {
    xb=MdrV0(x,0);
    xe=MdrV0(x,0);
    for (i=1;i<nx;i++)
    {
      if (xe < MdrV0(x,i))
      {
        xe = MdrV0(x,i);
        continue;
      }
      if (xb > MdrV0(x,i))
      {
        xb = MdrV0(x,i);
        continue;
      }
    }
  }
  if (ye==0 && yb==0 && ye==yb)
  {
    yb=MdrV0(y,0);
    ye=MdrV0(y,0);
    for (i=1;i<ny;i++)
    {
      if (ye < MdrV0(y,i))
      {
        ye = MdrV0(y,i);
        continue;
      }
      if (yb > MdrV0(y,i))
      {
        yb = MdrV0(y,i);
        continue;
      }
    }
  }

  //
  // call the appropriate solver
  //
  if (nx == nz && ny == nz)
  {
    //
    // scatter fit
    //
    if (s)
    {
      // smothing: find the size and position of knots
      if (mynxest)
        nxest = mynxest;
      else
        nxest = idegx + 1 + (int) (nx/2);
      if (mynyest)
        nyest = mynyest;
      else
        nyest = idegy + 1 + (int) (nx/2);

      nmax = MAX(nxest,nyest);

      tx    = mdr_Create(nmax,1); mdr_Zero(tx);
      ty    = mdr_Create(nmax,1); mdr_Zero(ty);
    }
    else
    {
      // least-squares: size and position of knots are given
      nmax = MAX(nxest, nyest);

      tx    = mdr_Create(nmax,1);
      for (i=0; i<nknotx-2*idegx; i++)
        MdrV0(tx,i+idegx) = MdrV0(ckx, i);

      ty    = mdr_Create(nmax,1);
      for (i=0; i<nknoty-2*idegy; i++)
        MdrV0(ty,i+idegy) = MdrV0(cky, i);
    }

    for (i=0; i<idegx+1; i++)
    {
      MdrV0(tx,i) = xb;
      MdrV0(tx,nxest-i-1) = xe;
    }
    for (i=0; i<idegy+1; i++)
    {
      MdrV0(ty,i) = yb;
      MdrV0(ty,nyest-i-1) = ye;
    }

    c    = mdr_Create((nxest-idegx-1)*(nyest-idegy-1),1);
    // dimension of the work array no.1
    //  u = nxest-kx-1, v = nyest-ky-1, km = MAX(kx,ky)+1,
    //  ne = MAX(nxest,nyest), bx = kx*v+ky+1, by = ky*u+kx+1,
    //  if(bx.le.by) b1 = bx, b2 = b1+v-ky
    //  if(bx.gt.by) b1 = by, b2 = b1+u-kx  then
    //  lwrk1 >= u*v*(2+b1+b2)+2*(u+v+km*(m+ne)+ne-kx-ky)+b2+1
    int iu=nxest-idegx-1, iv=nyest-idegy-1, ikm=MAX(idegx,idegy)+1;
    int ine=MAX(nxest,nyest), ibx=idegx*iv+idegy+1, iby=idegy*iu+idegx+1;
    int ib1, ib2;
    if (ibx <= iby)
    {
      ib1 = ibx;
      ib2 = ib1+iv-idegy;
    }
    else
    {
      ib1 = iby;
      ib2 = ib1+iu-idegx;
    }
    lwrk1 = iu*iv*(2+ib1+ib2)+2*(iu+iv+ikm*(nx+ine)+ine-idegx-idegy)+ib2+1;
    wrk1  = mdr_Create(lwrk1,1);
    lwrk2 = iu*iv*(ib2+1)+ib2;
    wrk2  = mdr_Create(lwrk2,1);
    liwrk = nx + (nxest-2*idegx-1)*(nyest-2*idegy-1);
    iwrk  = mdi_Create(liwrk,1);

    SURFIT ( &iopt, &nx, MDRPTR(x), MDRPTR(y), MDRPTR(z), MDRPTR(w), &xb, &xe, &yb, &ye,
             &idegx, &idegy, &s, &nxest, &nyest, &nmax, &deps, &nknotx, MDRPTR(tx),
             &nknoty, MDRPTR(ty), MDRPTR(c), &fp,
             MDRPTR(wrk1), &lwrk1, MDRPTR(wrk2), &lwrk2, MDIPTR(iwrk), &liwrk, &ier,
             &maxit, &tol
           );

    if (ier <= 3 && nknotx > 0 && nknoty > 0)
    {
      // tell user if she messed up by expecting too much:
      if (ier==2)
        fprintf(stdout,"bsplinefit: theoretical limit encountered. "
            "smoothing parameter too small\n");
      else if (ier==3)
        fprintf(stdout,"bsplinefit: maximal number of iterations reached\n");

      rettx = mdr_Create (1,nknotx-2*idegx);
      retty = mdr_Create (1,nknoty-2*idegy);
      retc  = mdr_Create (1,(nknotx-idegx-1)*(nknoty-idegy-1));

      for (i=0; i<nknotx-2*idegx; i++)
        MdrV0(rettx,i) = MdrV0(tx,i+idegx);
      for (i=0; i<nknoty-2*idegy; i++)
        MdrV0(retty,i) = MdrV0(ty,i+idegy);

      for (i=0; i<(nknotx-idegx-1)*(nknoty-idegy-1); i++)
        MdrV0(retc,i) = MdrV0(c,i);
    }
    else
    {
      // ier 10 appeared, that is, the solver crapped out.
      fprintf(stdout,"bsplinefit: solver failed with error no. %i\n", ier);
      retc = rettx = retty = mdr_Create (0,0);
      idegx = 0;
      idegy = 0;
      nknotx = 0;
      nknoty = 0;
    }
  }
  else if (nx == ncz && ny == nrz)
  {
    //
    // grid fit
    //
    if (s)
    {
      // smothing: find the size and position of knots
      if (mynxest)
        nxest = mynxest;
      else
        nxest = idegx + 1 + (int) (nx/2);
      if (mynyest)
        nyest = mynyest;
      else
        nyest = idegy + 1 + (int) (nx/2);

      tx    = mdr_Create(nxest,1);
      ty    = mdr_Create(nyest,1);
    }
    else
    {
      // least-squares: size and position of knots are given
      tx    = mdr_Create(nxest,1);
      for (i=0; i<nknotx-2*idegx; i++)
        MdrV0(tx,i+idegx) = MdrV0(ckx, i);

      ty    = mdr_Create(nyest,1);
      for (i=0; i<nknoty-2*idegy; i++)
        MdrV0(ty,i+idegy) = MdrV0(cky, i);
    }

    for (i=0; i<idegx+1; i++)
    {
      MdrV0(tx,i) = xb;
      MdrV0(tx,nxest-i-1) = xe;
    }
    for (i=0; i<idegy+1; i++)
    {
      MdrV0(ty,i) = yb;
      MdrV0(ty,nyest-i-1) = ye;
    }

    c    = mdr_Create((nxest-idegx-1)*(nyest-idegy-1),1);
    int iu=MAX(ny,nxest);

    lwrk1 = 4+nxest*(ny+2*idegx+5)+nyest*(2*idegy+5)+nx*(idegx+1)+
        ny*(idegy+1) + iu;
    wrk1  = mdr_Create(lwrk1,1);
    liwrk = 3+nx+ny+nxest+nyest;
    iwrk  = mdi_Create(liwrk,1);

    REGRID (&iopt, &nx, MDRPTR(x), &ny, MDRPTR(y), MDRPTR(z), &xb, &xe, &yb, &ye,
            &idegx, &idegy, &s, &nxest, &nyest, &nknotx, MDRPTR(tx),
            &nknoty, MDRPTR(ty), MDRPTR(c), &fp, MDRPTR(wrk1), &lwrk1, MDIPTR(iwrk),
            &liwrk, &ier, &maxit, &tol
           );

    if (ier <= 3 && nknotx > 0 && nknoty > 0)
    {
      // tell user if she messed up by expecting too much:
      if (ier==2)
        fprintf(stdout,"bsplinefit2: theoretical limit encountered. "
            "smoothing parameter too small\n");
      else if (ier==3)
        fprintf(stdout,"bsplinefit2: maximal number of iterations reached\n");

      rettx = mdr_Create (1,nknotx-2*idegx);
      retty = mdr_Create (1,nknoty-2*idegy);
      retc  = mdr_Create (1,(nknotx-idegx-1)*(nknoty-idegy-1));

      for (i=0; i<nknotx-2*idegx; i++)
        MdrV0(rettx,i) = MdrV0(tx,i+idegx);
      for (i=0; i<nknoty-2*idegy; i++)
        MdrV0(retty,i) = MdrV0(ty,i+idegy);

      for (i=0; i<(nknotx-idegx-1)*(nknoty-idegy-1); i++)
        MdrV0(retc,i) = MdrV0(c,i);
    }
    else
    {
      // ier 10 appeared, that is, the solver crapped out.
      fprintf(stdout,"bsplinefit2: solver failed with error no. %i\n", ier);
      retc = rettx = retty = mdr_Create (0,0);
      idegx = 0;
      idegy = 0;
      nknotx = 0;
      nknoty = 0;
    }
  }
  else
    rerror ("bsplinefit2: dimension mismatch");
  //
  // clean stuff
  //
  mdr_Destroy (w);
  mdr_Destroy (wrk1);
  if (wrk2)
    mdr_Destroy (wrk2);
  mdr_Destroy (iwrk);
  mdr_Destroy (tx); mdr_Destroy (ty);
  mdr_Destroy (c);

  ent_Clean (e1);
  ent_Clean (e2);
  ent_Clean (e3);
  ent_Clean (e4);
  ent_Clean (e5);

  //
  // write up
  //
  Btree *bw = btree_Create ();

  // knot = <<1;2>>
  Btree *bw2 = btree_Create ();

  install  (bw2, "1", ent_Assign_Rlab_MDR(rettx));
  install  (bw2, "2", ent_Assign_Rlab_MDR(retty));
  install  (bw, "knot", ent_Assign_Rlab_BTREE(bw2));
  // coef
  install  (bw, "coef", ent_Assign_Rlab_MDR(retc));
  // degree= <<1;2>>
  Btree *bw3 = btree_Create ();
  install  (bw3, "1", ent_Create_Rlab_Double(idegx));
  install  (bw3, "2", ent_Create_Rlab_Double(idegy));
  install  (bw, "degree", ent_Assign_Rlab_BTREE(bw3));
  // error
  install  (bw, "residual", ent_Create_Rlab_Double(fp));
  // status
  install (bw, "status", ent_Create_Rlab_Double(ier));

  return ent_Assign_Rlab_BTREE(bw);
}

#undef THIS_SOLVER
#define THIS_SOLVER "bsplineval2"

Ent *
ent_bspline_val2 (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *e3=0, *e4=0, *e5=0;

  MDR *x=0, *y=0, *w, *tx=0, *ty=0, *c=0, *wrk=0;
  MDR *ttx=0, *tty=0, *iwrk=0;

  int i, j, ntx, nty, nc;
  int idegx=0, idegy=0, nxknot=0, nyknot=0;
  int iderx=0,idery=0;

  ListNode *node, *node2;

  if (nargs != 3 && nargs != 4)
    rerror ("bsplineval2: requires three or four arguments");

  //
  // get X
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) == MATRIX_DENSE_REAL)
    x = class_matrix_real(e1);
  else if (ent_type (e1) == BTREE)
  {
    // x.val
    node = btree_FindNode (ent_data (e2), "val");
    if (node != 0)
    {
      if (ent_type (var_ent (node)) == MATRIX_DENSE_REAL)
        x = class_matrix_real(var_ent (node));
      else
        rerror ("bsplineval2: list entry 'val' is not real");
    }
    else
      rerror ("bsplineval2: list entry 'val' is missing");
  }
  else
    rerror ("bsplineval2: missing first argument");

  if(MNR(x)!=1 && MNC(x)!=1)
    rerror ("bsplineval2: list entry 'val' is not vector");

  //
  // get Y
  //
  e2 = bltin_get_ent (args[1]);
  if (ent_type (e2) == MATRIX_DENSE_REAL)
    y = class_matrix_real(e2);
  else if (ent_type (e2) == BTREE)
  {
    // y.val
    node = btree_FindNode (ent_data (e2), "val");
    if (node != 0)
    {
      if (ent_type (var_ent (node)) == MATRIX_DENSE_REAL)
        y = class_matrix_real(var_ent (node));
      else
        rerror ("bsplineval2: list entry 'val' is not real");
    }
    else
      rerror ("bsplineval2: list entry 'val' is missing");
  }
  else
    rerror ("bsplineval2: missing first argument");

  if(MNR(y)!=1 && MNC(y)!=1)
    rerror ("bsplineval2: list entry 'val' is not vector");

  //
  // get 'm' and 'spldeg' from the second argument
  //
  e3 = bltin_get_ent (args[2]);
  if (ent_type (e3) == BTREE)
  {
    // degree = <<1;2>>
    node = btree_FindNode (ent_data(e3), "degree");
    if (node != 0)
    {
      if (ent_type(var_ent (node)) == BTREE)
      {
        // degree.1
        node2 = btree_FindNode (ent_data(var_ent(node)), "1");
        if (node2 != 0)
          idegx = (int) class_double (var_ent (node2));
        else
          rerror("bsplineval2: incorrect entry 'degree.1'");
        // degree.2
        node2 = btree_FindNode (ent_data(var_ent(node)), "2");
        if (node2 != 0)
          idegy = (int) class_double (var_ent (node2));
        else
          rerror("bsplineval2: incorrect entry 'degree.2'");
      }
    }
    else
      rerror("bsplineval2: incorrect entry 'degree'");

    // coef
    node = btree_FindNode (ent_data (e3), "coef");
    if (node != 0)
      c = class_matrix_real (var_ent (node));
    else
      rerror ("bsplineval2: missing list entry 'coef'");

    // knot = <<1;2>>
    node = btree_FindNode (ent_data(e3), "knot");
    if (node != 0)
    {
      if (ent_type(var_ent (node)) == BTREE)
      {
        // degree.1
        node2 = btree_FindNode (ent_data(var_ent(node)), "1");
        if (node2 != 0)
          tx = class_matrix_real (var_ent (node2));
        else
          rerror("bsplineval2: incorrect entry 'knot.1'");
        // degree.2
        node2 = btree_FindNode (ent_data(var_ent(node)), "2");
        if (node2 != 0)
          ty = class_matrix_real (var_ent (node2));
        else
          rerror("bsplineval2: incorrect entry 'knot.2'");
      }
    }
  }
  else
    rerror ("bsplineval2: second argument is not b-spline list");

  if (!idegx || !idegy)
    rerror ("bsplineval2: missing list entry 'degree'");

  nc  = MNR(c) * MNC(c);    // nc  = (nxknot - idegx - 1)*(nyknot - idegy - 1)
  ntx = MNR(tx) * MNC(tx);  // ntx = nxknot - 2 * idegx
  nty = MNR(ty) * MNC(ty);  // nty = nyknot - 2 * idegy
  if (nc != (ntx+idegx-1)*(nty+idegy-1))
    rerror ("bsplineval2: dimension mismatch between knots and spline coefficients");

  //
  // get IDER
  //
  if (nargs == 4)
  {
    e4 = bltin_get_ent (args[3]);
    if (ent_type (e4) == MATRIX_DENSE_REAL)
      iderx = (int) class_double (e4);
  }
  if (nargs == 5)
  {
    e5 = bltin_get_ent (args[4]);
    if (ent_type (e5) == MATRIX_DENSE_REAL)
      idery = (int) class_double (e5);
  }

  w = mdr_Create ( MAX(MNR(x),MNR(y)), MAX(MNC(x),MNC(y)));

  nxknot = ntx + 2 * idegx;
  ttx = mdr_Create (1, nxknot);
  for (i=0; i<ntx; i++)
    MdrV0(ttx,i+idegx) = MdrV0(tx,i);
  for (i=0; i<idegx+1; i++)
  {
    MdrV0(ttx,i)=MdrV0(ttx,idegx);
    MdrV0(ttx,nxknot-i)=MdrV0(ttx,nxknot-idegx-1);
  }
  nyknot = nty + 2 * idegy;
  tty = mdr_Create (1, nyknot);
  for (i=0; i<nty; i++)
    MdrV0(tty,i+idegy) = MdrV0(ty,i);
  for (i=0; i<idegy+1; i++)
  {
    MdrV0(tty,i)=MdrV0(tty,idegy);
    MdrV0(tty,nyknot-i-1)=MdrV0(tty,nyknot-idegy-1);
  }

//   printf("ttx = \n");
//   mdr_Print(ttx,stdout);
//   printf("tty = \n");
//   mdr_Print(tty,stdout);

  int ier =  0;
  int lwrk = (idegx+1-iderx)+(idegy+1-idery)+(nxknot-idegx-1)*(nyknot-idegy-1);
  int liwrk = 12;
  int ione = 1;
  int ix, jx, iy, jy;
  double xd, yd, zd;

  wrk  = mdr_Create(1,lwrk);
  iwrk = mdi_Create(1,liwrk);

  for (i=0; i<MAX(MNR(x),MNR(y)); i++)
    for (j=0; j<MAX(MNC(x),MNC(y)); j++)
  {
    ix = MIN(i,MNR(x)-1);
    jx = MIN(j,MNC(x)-1);
    iy = MIN(i,MNR(y)-1);
    jy = MIN(j,MNC(y)-1);
    xd = Mdr0(x,ix,jx);
    yd = Mdr0(y,iy,jy);

    DXPARDER (MDRPTR(ttx), &nxknot, MDRPTR(tty), &nyknot, MDRPTR(c), &idegx, &idegy,
              &iderx, &idery, &xd, &ione, &yd, &ione, &zd, MDRPTR(wrk), &lwrk,
              MDIPTR(iwrk), &liwrk, &ier);

    Mdr0(w,i,j) = zd;
  }

  mdr_Destroy (wrk);
  mdr_Destroy (ttx);
  mdr_Destroy (tty);

  // clean stuff
  ent_Clean (e1);
  ent_Clean (e2);
  ent_Clean (e3);
  ent_Clean (e4);
  ent_Clean (e5);


  return ent_Assign_Rlab_MDR(w);
}

#undef THIS_SOLVER
#define THIS_SOLVER "bsplinefit"
Ent *
ent_bspline_fit (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *e3=0, *e4=0;
  MDR *x=0, *y=0, *w=0, *rett=0, *retc=0, *x3=0, *con=0;
  MDR *t=0, *c=0, *wrk=0, *iwrk=0;

  ListNode *node;

  int i, ny=0, nx=0, ideg=3, iopt=0, nknot=0, iper=0, iprimeb=0, iprimee=0, nest=0;
  int icon=0, lwrk, maxit=40, idummy, mynest=0, nx3=0;
  double s=0, xe=0, xb=0, fp, ypi=0, ypf=0, tol=0.001, ddummy;

  if (nargs != 4 && nargs != 3)
  {
    rerror (THIS_SOLVER ": Requires three or four arguments !");
  }

  // Y
  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) != MATRIX_DENSE_REAL && ent_type (e1) != BTREE)
    rerror (THIS_SOLVER ": First argument 'y' must be real vector or list <<val;wgt>> !");

  if (ent_type (e1) == MATRIX_DENSE_REAL)
  {
    y = class_matrix_real (e1);
    if (!y)
      rerror (THIS_SOLVER ": First argument 'y' must be real vector or list <<val;wgt>> !");
    ny = MNR (y) * MNC (y);
    if (!ny)
      rerror (THIS_SOLVER ": First argument 'y' must be real vector or list <<val;wgt>> !");
  }
  else
  {
    // val
    node = btree_FindNode (ent_data (e1), "val");
    if (!node)
      rerror (THIS_SOLVER ": First argument 'y' must be real vector or list <<val;wgt>> !");
    y = class_matrix_real (var_ent (node));
    if (!y)
      rerror (THIS_SOLVER ": First argument 'y' must be real vector or list <<val;wgt>> !");
    ny = MNR (y) * MNC (y);
    if (!ny)
      rerror (THIS_SOLVER ": First argument 'y' must be real vector or list <<val;wgt>> !");

    // wgt
    node = btree_FindNode (ent_data (e1), "wgt");
    if (!node)
      rerror (THIS_SOLVER ": First argument 'y' must be real vector or list <<val;wgt>> !");
    w = mdr_Float_BF (var_ent (node));
    if (ny != MNR (w)* MNC (w))
      rerror (THIS_SOLVER ": First argument 'y' must be real vector or list <<val;wgt>> !");
  }

  // X
  e2 = bltin_get_ent (args[1]);
  if (ent_type (e2) != MATRIX_DENSE_REAL && ent_type (e2) != BTREE)
    rerror (THIS_SOLVER ": Second argument 'x' must be real vector or list <<val;wgt>> !");

  if (ent_type (e2) == MATRIX_DENSE_REAL)
  {
    //
    // mean of a vector or of a matrix (columnwise)
    //
    x = class_matrix_real (e2);
    if (!x)
      rerror (THIS_SOLVER ": Second argument 'x' must be real vector or list <<val;wgt>> !");
    nx = MNR (x) * MNC (x);
    if (!nx)
      rerror (THIS_SOLVER ": Second argument 'x' must be real vector or list <<val;wgt>> !");
  }
  else
  {
    // val
    node = btree_FindNode (ent_data (e2), "val");
    if (!node)
      rerror (THIS_SOLVER ": Second argument 'x' must be real vector or list <<val;wgt>> !");
    x = class_matrix_real (var_ent (node));
    if (!x)
      rerror (THIS_SOLVER ": Second argument 'x' must be real vector or list <<val;wgt>> !");
    nx = MNR (x) * MNC (x);
    if (!nx)
      rerror (THIS_SOLVER ": Second argument 'x' must be real vector or list <<val;wgt>> !");
  }

  //
  // smoothing parameter or knots for ls b-spline
  //
  e3 = bltin_get_ent (args[2]);
  if (ent_type (e3) != MATRIX_DENSE_REAL)
    rerror (THIS_SOLVER ": Third argument 's' must be positive scalar !");
  x3 = ent_data (e3);
  if (MNR(x3) * MNC(x3) == 1)
    s = MdrV0(x3,0);

  //
  // options
  //
  if (nargs > 3)
  {
    e4 = bltin_get_ent (args[3]);
    if (ent_type (e4) == BTREE)
    {
      // degree
      node = btree_FindNode (ent_data(e4), RLAB_NAME_DIERCKX_DEGREE);
      if (node != 0)
      {
        ideg = (int) class_double (var_ent (node));
        if (ideg != 1 && ideg != 3 && ideg != 5)
          ideg = 3;
      }
      // range
      node = btree_FindNode (ent_data (e4), RLAB_NAME_DIERCKX_RANGE);
      if (node != 0)
      {
        if (ent_type(var_ent(node)) == MATRIX_DENSE_REAL)
        {
          MDR *xrange = class_matrix_real (var_ent (node));
          if (MNR(xrange)*MNC(xrange)==2)
          {
            xb = MdrV0(xrange,0);
            xe = MdrV0(xrange,1);
          }
        }
      }
      // periodic
      node = btree_FindNode (ent_data(e4), RLAB_NAME_DIERCKX_PER);
      if (node != 0)
      {
        iper = (int) class_double (var_ent (node));
        if (iper != 0 && iper != 1)
          iper = 0;
      }
      // tol
      node = btree_FindNode (ent_data(e4), RLAB_NAME_DIERCKX_TOL);
      if (node != 0)
      {
        ddummy = class_double (var_ent (node));
        if (ddummy > 0)
          tol = ddummy;
      }
      // maxi
      node = btree_FindNode (ent_data(e4), RLAB_NAME_DIERCKX_MAXITER);
      if (node != 0)
      {
        idummy = (int) class_double (var_ent (node));
        if (idummy > 0)
          maxit = idummy;
      }
      // max no. of knots: only for smooth fit (s!=0)
      node = btree_FindNode (ent_data(e4), RLAB_NAME_DIERCKX_MAXKNOTS);
      if (node != 0)
      {
        idummy = (int) class_double (var_ent (node));
        if (idummy > 0 && s)
          mynest = idummy;
      }
      // derivatives at endpoints
      node = btree_FindNode (ent_data(e4), RLAB_NAME_DIERCKX_YP1);
      if (node != 0)
      {
        iprimeb=1;
        ypi = class_double (var_ent (node));
      }
      node = btree_FindNode (ent_data(e4), RLAB_NAME_DIERCKX_YP2);
      if (node != 0)
      {
        iprimee=1;
        ypf = class_double (var_ent (node));
      }
      // constraint
      node = btree_FindNode (ent_data (e4), RLAB_NAME_DIERCKX_CONVEX);
      if (node != 0)
      {
        if (ent_type(var_ent(node)) == MATRIX_DENSE_REAL)
        {
          con = class_matrix_real (var_ent (node));
          if (MNR(con) * MNC(con) == nx)
          {
            // constraint on x
            icon = 1;
          }
          else if (MNR(con) * MNC(con) == MNR(x3) * MNC(x3))
          {
            // constraint on knots, if given
            icon = 2;
          }
        }
      }
    }
  }

  // check weight matrix
  if (!w)
  {
    w = mdr_Create(ny,1);
    mdr_Ones (w);
  }


  int ier=0;

  if (s==0)
  {
    // least-squares fit through the knots
    // x3 is a vector containing the position of the knots: two endpoints
    // represent xb and xe, while the points in between represent the
    // actual knots
    iopt = -1;    // least-squares fit spline
    nx3  = MNR(x3) * MNC (x3);
    xb = MdrV1(x3,1);
    xe = MdrV1(x3,nx3);
    if (nx3 < 2*ideg)
      rerror ("bsplinefit: size of the knot array too small");
  }
  else
  {
    // smooth fit
    // if user did not specify the endpoints, these are chosen
    // from the data set
    if (xe==0 && xb==0 && xe==xb)
    {
      xb=MdrV0(x,0);
      xe=MdrV0(x,nx-1);
    }
  }

  //
  // call the appropriate solver
  //
  if (iprimee||iprimeb)
  {
    //
    // boundary conditions at the endpoint(s)
    //
    double de[2]={0.0};
    double db[2]={0.0};

    // is there a boundary condition at xlo?
    db[0] = MdrV0(y,0);
    if (iprimeb)
    {
      // yes
      iprimeb = 2;
      db[1] = ypi;
    }
    else
      iprimeb = 1; // no

    // is there a boundary condition at xhi?
    de[0] = MdrV0(y,nx-1);
    if (iprimee)
    {
      // yes
      iprimee = 2;
      de[1] = ypf;
    }
    else
      iprimee = 1; // no
    //
    // one or both endpoint derivatives are used
    //
    if(s)
    {
        // s!=0: find knots via smoothing parameter
      if (mynest)
        nest = mynest;
      else
        nest = nx + ideg + 1;
    }
    else
    {
      // s==0: user specified the knots, do the least-squares fit
      nest  = nx3 + 2*ideg;
      nknot = nest;
    }

    t = mdr_Create(nest,1);
    c = mdr_Create(nest,1);
    if (s==0)
    {
      for (i=0; i<nx3; i++)
        MdrV0(t,i+ideg) = MdrV0(x3, i);
    }

    int ione=1;
    int ncp = 2*(ideg+1);
    MDR *xx = mdr_Create (1, nx);
    MDR *cp = mdr_Create (1, ncp); mdr_Zero(cp);

    lwrk = nx * (ideg + 1) + nest * (7+3*ideg);
    wrk  = mdr_Create(lwrk,1);
    iwrk = mdi_Create(nest,1);

    CONCUR ( &iopt, &ione, &nx, MDRPTR(x), &nx, MDRPTR(y), MDRPTR(xx), MDRPTR(w),
             &iprimeb, db, &iprimeb, &iprimee, de, &iprimee,
             &ideg, &s, &nest, &nknot, MDRPTR(t), &nest, MDRPTR(c),
             &ncp, MDRPTR(cp), &fp, MDRPTR(wrk), &lwrk, MDIPTR(iwrk), &ier,
             &maxit, &tol);

    mdr_Destroy(xx);
    mdr_Destroy(cp);
  }
  else
  {
    //
    // no endpoint derivatives are given
    //
    if (iper)
    {
      //
      // periodic
      //
      if(s)
      {
        // s!=0: find knots via smoothing parameter
        if (mynest)
          nest = mynest;
        else
          nest = nx + 2*ideg;
      }
      else
      {
        // s==0: user specified the knots, do the least-squares fit
        nest  = nx3 + 2*ideg;
      }
      nknot = nest;

      lwrk = nx*(ideg+1)+nest*(8+5*ideg);

      wrk  = mdr_Create(lwrk,1);
      iwrk = mdi_Create(nest,1);
      t    = mdr_Create(nest,1);
      c    = mdr_Create(nest,1);
      if (s==0)
      {
        for (i=0; i<nx3; i++)
          MdrV0(t,i+ideg) = MdrV0(x3, i);
      }
      else
      {
        for (i=0; i<ideg+1; i++)
        {
          MdrV0(t,i) = MdrV0(x, 0);
          MdrV0(t,nx-i-1) = MdrV0(x, nx-1);
        }
      }

      PERCUR ( &iopt, &nx, MDRPTR(x), MDRPTR(y), MDRPTR(w), &ideg, &s, &nest, &nknot,
                MDRPTR(t), MDRPTR(c), &fp, MDRPTR(wrk), &lwrk, MDIPTR(iwrk), &ier,
                &maxit, &tol);
    }
    else
    {
      // are the convexity constraints given
      if (icon == 0)
      {
        //
        // non-constrained:
        //
        if(s)
        {
          // s!=0: did user specify maximal number of knots in the search
          if (mynest)
            nest = mynest;
          else
            nest = nx + ideg + 1;
        }
        else
        {
          // s==0: user specified the knots, do the least-squares fit
          nest  = nx3 + 2*ideg;
        }
        nknot = nest;

        lwrk = nx * (ideg + 1) + nest * (7+3*ideg);

        t    = mdr_Create(nest,1); mdr_Zero(t);
        c    = mdr_Create(nest,1); mdr_Zero(c);
        wrk  = mdr_Create(lwrk,1);
        iwrk = mdi_Create(nest,1);
        if (s==0)
        {
          for (i=0; i<nx3; i++)
            MdrV0(t,i+ideg) = MdrV0(x3, i);
        }

        CURFIT ( &iopt, &nx, MDRPTR(x), MDRPTR(y), MDRPTR(w), &xb, &xe, &ideg, &s, &nest,
                  &nknot, MDRPTR(t), MDRPTR(c), &fp, MDRPTR(wrk), &lwrk, MDIPTR(iwrk), &ier,
                  &maxit, &tol);
      }
      else if (icon == 1 && con)
      {
        //
        // local convexity constraints for x
        //
        // did user specify maximal number of knots in the search
        if(s)
        {
          // s!=0: did user specify maximal number of knots in the search
          if (mynest)
            nest = mynest;
          else
            nest = nx + 4;
        }
        else
        {
          // s==0: user specified the knots, do the least-squares fit
          nest  = nx3 + 2*ideg;
          nknot = nest;
        }

        t = mdr_Create(nest,1);
        c = mdr_Create(nest,1);
        if (s==0)
        {
          for (i=0; i<nx3; i++)
            MdrV0(t,i+ideg) = MdrV0(x3, i);
        }

        int maxbin = nest - 6;
        int maxtr  = gsl_sf_choose(nest-5, (int)((nest-6.0)/2.0+0.5)) +
            gsl_sf_choose(nest-6, (int)((nest-6.0)/2.0+0.5)+1);
        MDR *sx   = mdr_Create (1, nx);
        MDR *bind = mdi_Create (1, nest);

        lwrk = nx*4 + nest*8 + maxbin*(maxbin+nest+1);
        wrk  = mdr_Create(lwrk,1);
        int liwrk = maxtr*4 + 2*(maxbin+1);
        iwrk = mdi_Create(liwrk,1);

        CONCON ( &iopt, &nx, MDRPTR(x), MDRPTR(y), MDRPTR(w), MDRPTR(con), &s, &nest,
                  &maxtr, &maxbin, &nknot, MDRPTR(t), MDRPTR(c), &fp, MDRPTR(sx),
                  MDIPTR(bind), MDRPTR(wrk), &lwrk, MDIPTR(iwrk), &liwrk, &ier);

        mdr_Destroy (sx);
        mdr_Destroy (bind);

      }
      else if (icon == 2 && con)
      {
        //
        // local convexity constraints for x
        //
        // did user specify maximal number of knots in the search
        if(s)
        {
          // s!=0: did user specify maximal number of knots in the search
          if (mynest)
            nest = mynest;
          else
            nest = nx + ideg + 1;
        }
        else
        {
          // s==0: user specified the knots, do the least-squares fit
          nest  = nx3 + 2*ideg;
          nknot = nest;
        }

        t = mdr_Create(nest,1);
        c = mdr_Create(nest,1);
        if (s==0)
        {
          for (i=0; i<nx3; i++)
            MdrV0(t,i+ideg) = MdrV0(x3, i);
        }

        // local convexity constraints for knots
        int maxbin = nknot - 6;
        int maxtr  = gsl_sf_choose(nest-5, (int)((nknot-6.0)/2.0+0.5)) +
            gsl_sf_choose(nknot-6, (int)((nknot-6.0)/2.0+0.5)+1);
        MDR *sx   = mdr_Create (1, nknot);
        MDR *bind = mdi_Create (1, nknot);

        lwrk = nx*4 + nknot*8 + maxbin*(maxbin+nknot+1);
        wrk  = mdr_Create(lwrk,1);
        int liwrk = maxtr*4 + 2*(maxbin+1);
        iwrk = mdi_Create(liwrk,1);

        COCOSP ( &nx, MDRPTR(x), MDRPTR(y), MDRPTR(w), &nknot, MDRPTR(t), MDRPTR(con),
                  &maxtr, &maxbin, MDRPTR(c), &fp, MDRPTR(sx), MDRPTR(bind),
                  MDRPTR(wrk), &lwrk, MDIPTR(iwrk), &liwrk, &ier);

        mdr_Destroy (sx);
        mdr_Destroy (bind);
      }
      else
      {
        rerror("bsplinefit: terrible internal error");
      }
    }
  }

  if (ier <= 3 && nknot > 0)
  {
    // tell user if she messed up by expecting too much:
    if (ier==2)
      fprintf(stdout,"bsplinefit: theoretical limit encountered."
          " smoothing parameter too small\n");
    else if (ier==3)
      fprintf(stdout,"bsplinefit: maximal number of iterations reached\n");

    rett = mdr_Create (1,nknot);
    for (i=0; i<nknot; i++)
      MdrV0(rett,i) = MdrV0(t,i);
    retc = mdr_Create (1,nknot-ideg-1);
    for (i=0; i<nknot-ideg-1; i++)
      MdrV0(retc,i) = MdrV0(c,i);
  }
  else
  {
    // ier 10 appeared, that is, the solver crapped out.
    fprintf(stdout,"bsplinefit: solver failed with error no. %i\n", ier);
    retc = rett = mdr_Create (0,0);
    ideg = 0;
    nknot = 0;
  }

  //
  // clean stuff
  //
  mdr_Destroy (w);
  mdr_Destroy (wrk);
  mdr_Destroy (iwrk);
  mdr_Destroy (t);
  mdr_Destroy (c);

  ent_Clean (e1);
  ent_Clean (e2);
  ent_Clean (e3);
  ent_Clean (e4);

  //
  // write up
  //
  Btree *bw = btree_Create ();
  // knot
  install (bw, "knot",      ent_Assign_Rlab_MDR(rett));
  // coef
  install (bw, "coef",      ent_Assign_Rlab_MDR(retc));
  // degree
  install  (bw, "degree",   ent_Create_Rlab_Double(ideg));
  // error
  install  (bw, "residual", ent_Create_Rlab_Double(fp));
  // status
  install  (bw, "status",   ent_Create_Rlab_Double(ier));

  return ent_Assign_Rlab_BTREE(bw);
}

#undef THIS_SOLVER
#define THIS_SOLVER "bsplineval"
Ent *
ent_bspline_val (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *e3=0;
  MDR *x=0, *w, *ider=0, *t=0, *c=0, *wrk=0;
  int j, nx, idf, nc;
  int ideg = 0, nknot = 3;

  ListNode *node;

  if (nargs != 2 && nargs != 3)
    rerror (THIS_SOLVER ": Requires two or three arguments");

  //
  // get X
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) != MATRIX_DENSE_REAL && ent_type (e1) != BTREE)
    rerror (THIS_SOLVER ": First argument must be real vector or list !");

  if (ent_type (e1) == MATRIX_DENSE_REAL)
    x = class_matrix_real(e1);
  else
  {
    // x.val
    node = btree_FindNode (ent_data (e1), "val");
    if (!node)
      rerror (THIS_SOLVER ": First argument must be list <<val;wgt>> !");
    if (ent_type (var_ent (node)) != MATRIX_DENSE_REAL)
      rerror (THIS_SOLVER ": First argument must be list <<val;wgt>> !");
    x = class_matrix_real(var_ent (node));
  }
  if(MNR(x)!=1 && MNC(x)!=1)
    rerror (THIS_SOLVER ": list entry 'val' is not vector");
  nx = MNR (x) * MNC (x);

  //
  // get 'm' and 'spldeg' from the second argument
  //
  e2 = bltin_get_ent (args[1]);
  if (ent_type (e2) != BTREE)
    rerror (THIS_SOLVER ": Second argument must be list <<degree;coef;knot>> !");

  // degree
  node = btree_FindNode (ent_data(e2), "degree");
  if (!node)
    rerror (THIS_SOLVER ": Second argument must be list <<degree;coef;knot>> !");
  ideg = (int) class_double (var_ent (node));
  if (!ideg)
    rerror (THIS_SOLVER ": Second argument must be list <<degree;coef;knot>> !");

  // coef
  node = btree_FindNode (ent_data (e2), "coef");
  if (!node )
    rerror (THIS_SOLVER ": Second argument must be list <<degree;coef;knot>> !");
  c = class_matrix_real (var_ent (node));
  if (!c)
    rerror (THIS_SOLVER ": Second argument must be list <<degree;coef;knot>> !");
  nc    = MNR(c) * MNC(c); // nc = nknot - ideg - 1
  if (!nc)
    rerror (THIS_SOLVER ": Entry  'coef' in second argument must be real vector of non-zero size !");

  // knot
  node = btree_FindNode (ent_data (e2), "knot");
  if (!node )
    rerror (THIS_SOLVER ": Second argument must be list <<degree;coef;knot>> !");
  t = class_matrix_real (var_ent (node));
  if (!t)
    rerror (THIS_SOLVER ": Second argument must be list <<degree;coef;knot>> !");
  nknot = MNR(t) * MNC(t); // nknot
  if (!nknot)
    rerror (THIS_SOLVER ": Entry  'coef' in second argument must be real vector of non-zero size");

  if (nknot != nc+ideg+1)
    rerror (THIS_SOLVER ": Dimension mismatch in second argument between entries 'coef' and 'knot'");

  //
  // get IDER
  //
  if (nargs == 3)
  {
    e3 = bltin_get_ent (args[2]);
    if (ent_type (e3) == MATRIX_DENSE_REAL)
      ider = (MDR *) mdr_Float_BF ( ent_data(e3) );
  }
  if (!ider) ider = mdr_CreateScalar (0.0);

  w = mdr_Create (nx, MNR (ider) * MNC (ider));

  idf = 0;
  int ier =0;
  wrk = mdr_Create(nknot,1);

  for (j = 0; j < MNR (ider) * MNC (ider); j++)
  {
    idf = MdrV0 (ider, j);
    DXSPLDER (MDRPTR(t), &nknot, MDRPTR(c), &ideg, &idf, MDRPTR(x), &MdrV0(w,j*nx), &nx, MDRPTR(wrk), &ier);
  }

  mdr_Destroy (wrk);
  mdr_Destroy (ider);

  // clean stuff
  ent_Clean (e1);
  ent_Clean (e2);
  ent_Clean (e3);

  return ent_Assign_Rlab_MDR(w);
}

