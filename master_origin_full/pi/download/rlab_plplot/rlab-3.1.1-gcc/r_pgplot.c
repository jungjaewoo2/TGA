/* r_pgplot.c */

/*
 * This file contains the PGPLOT interface routines.
 */

/*  This file is a part of RLaB ("Our"-LaB)
  Copyright (C) 1995  Ian R. Searle
  Copyright (C) 2014-2017  Marijan Kostrun

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

#include "rlab_solver_parameters_names.h"
#include "rlab_macros.h"
#include "rlab_macros_code.h"

#include <stdio.h>
#include "pgplot/cpgplot.h"

/* **************************************************************
* PLPLOT related functions and interface(s)
* ************************************************************** */
#define RLABPLUS_PGPLOT_IDX_LINEPOINT   0
#define RLABPLUS_PGPLOT_IDX_LINETYPE    1
#define RLABPLUS_PGPLOT_IDX_LINECOLOR   2
#define RLABPLUS_PGPLOT_IDX_LINEWIDTH   3
#define RLABPLUS_PGPLOT_IDX_POINTTYPE   4
#define RLABPLUS_PGPLOT_IDX_POINTCOLOR  5
#define RLABPLUS_PGPLOT_IDX_POINTSIZE   6

extern int isnand( double );
static float *mdr_coerce_floatptr (MDR * m);


#undef  THIS_SOLVER
#define THIS_SOLVER "_pgdrawgrad"
Ent * ent_pgplot_draw_grad_plane (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *e3=0, *e4=0, *e5=0;
  MDR *fmt=0;
  char *desc=0;
  FILE *rlab_stderr = (!RLAB_STDERR_DS) ? stderr : RLAB_STDERR_DS;
  float scale=1.0;

  if ((nargs != 3)&&(nargs != 4)&&(nargs != 5))
  {
    fprintf (rlab_stderr, THIS_SOLVER ": " RLAB_ERROR_THREE_ARG_REQUIRED "\n");
    goto _exit_here;
  }
  e1 = bltin_get_ent (args[0]);
  if (ent_type(e1)!=MATRIX_DENSE_REAL)
    goto _exit_here;
  fmt = ent_data(e1);
  if (SIZE(fmt)!=10)
    goto _exit_here;

  e2 = bltin_get_ent (args[1]);
  if (ent_type(e2)!=MATRIX_DENSE_STRING)
    goto _exit_here;
  desc = class_char_pointer(e2);

  e3 = bltin_get_ent (args[2]);
  if (ent_type(e3)!=MATRIX_DENSE_REAL)
    goto _exit_here;
  scale = class_double(e3);

  // store worldview
  float x1_w,x2_w,y1_w,y2_w,ch_old;
  cpgqch  (&ch_old);
  cpgqwin(&x1_w,&x2_w,&y1_w,&y2_w);

#define GRAD_DIMY (49)
#define GRAD_DIMX (2)
  float a[GRAD_DIMY*GRAD_DIMX]={0.0}, grad_slope, zmin, zmax;
  int i,j,nc,nr,palette_param[7];
  float tr[6];

  zmin = mdrV0(fmt,1);
  zmax = mdrV0(fmt,2);

//   printf(THIS_SOLVER ": zmin,zmax= %f,%f\n", zmin, zmax);
//   printf(THIS_SOLVER ": fmt = %i\n", mdiV0(fmt,0));

  if ((mdiV0(fmt,0) == 3) && !isnand(zmin) && !isnand(zmax))
  {
    //
    // pm3d legend is vertical with ticks on the righthand side
    //
    nc=GRAD_DIMY;
    nr=GRAD_DIMX;
    grad_slope=(zmax-zmin)/(GRAD_DIMY-1.0);

    for (i=0;i<7;i++)
      palette_param[i] = mdiV0(fmt,3+i);

    // set worldview for this plot
    cpgswin(0.0, 1.0, zmin, zmax);

    // x_i = [0,1/2,1], i=1,2,3
    tr[0] = -1/(nr-1);
    tr[1] =  1/(nr-1);
    tr[2] =  0.0;
    // y_i = [0,1/N,..1], i=1,2,3,GRAD_DIMY
    tr[3] =  zmin - grad_slope;
    tr[4] =  0.0;
    tr[5] =  grad_slope;

    for (i=0;i<nc;i++)
    {
      a[nr*i]=mdrV0(fmt,1) + grad_slope * i;
      for (j=1;j<nr;j++)
        a[j+nr*i]=a[nr*i];
    }
    cpgrgbfor (a, nr, nc, 1, nr, 1, nc, mdrV0(fmt,2), mdrV0(fmt,1), tr, palette_param);

    if (isvalidstring(desc))
    {
      cpgsch (scale);
      cpgmtxt("T", 1.125, 0.5, 0.5, desc);
    }

    int idummy[]={0,-1,-1,-1};
    float fdummy[]={0.0,-1.0,-1.0,-1.0};

    cpgbox ("bc", fdummy, idummy, "bcmv", fdummy, idummy);
  }
  else if (mdiV0(fmt,0) == 4)
  {
    MDR *surf=0;
    if (nargs>3)
    {
      e4 = bltin_get_ent(args[3]);
      if (ent_type(e4) == MATRIX_DENSE_REAL)
      {
        surf = ent_data(e4); // surf=[x1,x2,y1,y2] region where the contour lines will be printed
      }
    }
    if (SIZE(surf)!=4)
      goto _exit_here;

    MDR *zlev=0;
    if (nargs>4)
    {
      e5 = bltin_get_ent(args[4]);
      if (ent_type(e5) == MATRIX_DENSE_REAL)
      {
        zlev = ent_data(e5);
      }
    }
    if (SIZE(zlev)<1)
      zlev = 0;

    //
    // contour legend comprise many short vertical lines
    //
    nc=GRAD_DIMX;
    nr=GRAD_DIMY;
    grad_slope=(zmax-zmin)/(GRAD_DIMY-1.0);
    float x_slope = (mdrV0(surf,1) - mdrV0(surf,0))/(nr-1.0);
    float y_slope = (mdrV0(surf,3) - mdrV0(surf,2))/(nc-1.0);

    for (i=0;i<6;i++)
      palette_param[i] = mdiV0(fmt,4+i);

    // z-range provided?
    int zstp = ABS(mdrV0(fmt,3));
    float *c = GC_MALLOC((zstp+2) * sizeof(float));
    if (zlev)
    {
      for (i=0; i<zstp; i++)
        c[i] = mdrV0(zlev,i);

      if (c[0]==zmin)
        c[0] += 0.1 * (c[1]-c[0]);

      if (c[zstp-1]==zmax)
        c[zstp-1] -= 0.1 * (c[zstp-1]-c[zstp-2]);
    }
    else
    {
      float z_slope = (zmax - zmin)/(zstp - 1.0);
      c[0] = zmin + 0.1 * z_slope;
      for (i=1; i<(zstp-1); i++)
      {
        c[i] = c[i-1] + z_slope;
      }
      c[zstp-1] = zmax - 0.1 * z_slope;
    }
    c[zstp] = zmin;
    c[zstp+1] = zmax;

    // x_i = [0,1/N,..1], i=1,2,3,GRAD_DIMX
    tr[0] =  mdrV0(surf,0) - x_slope;
    tr[1] =  x_slope;
    tr[2] =  0.0;
    // y_i = [0,1/2,1], i=1,2,3
    tr[3] =  mdrV0(surf,2) - y_slope;
    tr[4] =  0.0;
    tr[5] =  y_slope;

//     printf(THIS_SOLVER ": zmin,zmax,zstp = %f,%f,%i\n", zmin, zmax,zstp);
//     printf(THIS_SOLVER ": tr =");
//     for (i=0;i<6;i++)
//       printf(" %f", tr[i]);
//     printf("\n");

    for (i=0;i<nr;i++)
    {
      a[i]= zmin + grad_slope * (nr-1-i);
      for (j=1;j<nc;j++)
        a[i+nr*j]=a[i];
    }

    cpgcont(a, nr, nc, 1, nr, 1, nc, c, zstp, tr, palette_param);

    if (c)
      GC_FREE(c);
  }

  // restore:
  cpgswin(x1_w,x2_w,y1_w,y2_w);
  cpgsch (ch_old);

_exit_here:

 ent_Clean (e1);
 ent_Clean (e2);
 ent_Clean (e3);
 ent_Clean (e4);
 ent_Clean (e5);

 return ent_Create_Rlab_Success();
}



#undef  THIS_SOLVER
#define THIS_SOLVER "_pgqcr"
Ent * ent_pgplot_pgqcr (int nargs, Datum args[])
{
  Ent *e1=0;
  int i, nr;
  MDR *color=0, *rgb=0;
  FILE *rlab_stderr = (!RLAB_STDERR_DS) ? stderr : RLAB_STDERR_DS;

  if (nargs != 0 && nargs!=1)
  {
    fprintf (rlab_stderr, THIS_SOLVER ": " RLAB_ERROR_NONE_OR_ONE_ARG_REQUIRED "\n");
    goto _exit_pgqcr;
  }

  if (nargs == 1)
  {
    e1 = bltin_get_ent (args[0]);
    if (ent_type(e1)==MATRIX_DENSE_REAL)
      color = ent_data(e1);
    nr = SIZE(color);
    if (nr>0)
    {
      float r,g,b;
      rgb = mdi_Create(nr,3);
      for (i=0; i<nr; i++)
      {
        cpgqcr(mdiV0(color,i), &r, &g, &b);
        Mdi0(rgb,i,0) = (int) r * 255.0;
        Mdi0(rgb,i,1) = (int) g * 255.0;
        Mdi0(rgb,i,2) = (int) b * 255.0;
      }

      goto _exit_pgqcr;
    }
  }

  nr = 16;
  rgb = mdi_Create(nr,3);
  for (i=0; i<nr; i++)
  {
    float r,g,b;
    rgb = mdi_Create(nr,3);
    for (i=0; i<nr; i++)
    {
      cpgqcr(mdiV0(color,i), &r, &g, &b);
      Mdi0(rgb,i,0) = (int) r * 255.0;
      Mdi0(rgb,i,1) = (int) g * 255.0;
      Mdi0(rgb,i,2) = (int) b * 255.0;
    }
  }

_exit_pgqcr:

  ent_Clean (e1);

  return ent_Assign_Rlab_MDR(rgb);
}


// 2D plot:
//  _pgprintf (data, format), where data is multicolumn real matrix
//  _pgprintf (data, format), where data is a list <<x;y;z>> for intensity plot
#undef THIS_SOLVER
#define THIS_SOLVER "_pgprintf"
Ent * ent_pgplot_printf (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *e3=0, *e4=0, *e5=0, *e6=0, *e7=0, *e8=0;
  int i, ndata, logx=0, logy=0, clean_col_idx=0;
  FILE *rlab_stderr = (!RLAB_STDERR_DS) ? stderr : RLAB_STDERR_DS;
  float *x=0, *y=0;

  MDR *plfmt=0, *pldata=0, *col_idx=0;
  MDS *labels=0;

  // store values before the call
  int old_lt;
  cpgqls (&old_lt);
  int old_lc;
  cpgqci(&old_lc);
  int old_lw;
  cpgqlw(&old_lw);
  float old_ps;
  cpgqch(&old_ps);

  if ((nargs < 2) || (nargs > 8))
  {
    fprintf(rlab_stderr, THIS_SOLVER ": " RLAB_ERROR_TWO_TO_EIGHT_ARG_REQUIRED "\n");
    goto _exit_pgprintf;
  }

  //
  // data set to be plotted
  //
  e1 = bltin_get_ent(args[0]);

  //
  // how to plot it
  //
  e2 = bltin_get_ent(args[1]);
  if (ent_type(e2) != MATRIX_DENSE_REAL)
  {
    fprintf(rlab_stderr, THIS_SOLVER ": " RLAB_ERROR_ARG2_MDR_VECTOR "\n");
    fprintf(rlab_stderr, THIS_SOLVER ": data type '%s' is not supported !\n", etd(e2));
    goto _exit_pgprintf;
  }
  plfmt = ent_data(e2);
  if (SIZE(plfmt)!=10)
  {
    fprintf(rlab_stderr, THIS_SOLVER ": " RLAB_ERROR_ARG2_MDR_VECTOR "\n");
    fprintf(rlab_stderr, THIS_SOLVER ": 10-component vector expected but %i-component was provided.\n"
        THIS_SOLVER ": Cannot continue!\n", MNC(plfmt));
    goto _exit_pgprintf;
  }

  if (nargs > 2)
  {
    e3 = bltin_get_ent(args[2]);
    if (ent_type(e3) == MATRIX_DENSE_REAL)
    {
      logx = (class_double(e3)>0);
    }
  }
  if (nargs > 3)
  {
    e4 = bltin_get_ent(args[3]);
    if (ent_type(e4) == MATRIX_DENSE_REAL)
    {
      logy = (class_double(e4)>0);
    }
  }

  if (nargs > 4)
  {
    e5 = bltin_get_ent(args[4]);
    if (ent_type(e5) == MATRIX_DENSE_REAL)
    {
      col_idx = ent_data(e5);
    }
    if (SIZE(col_idx)<1)
      col_idx = 0;
  }

  //
  // plot data list <<x;y,z>>
  // plfmt = [lpb, zmin, zmax, palette_mode, 0, 0, 0, 0, 0];
  if (ent_type(e1) == BTREE)
  {
    ListNode * node=0;
    int nr, nc;
    MDR *x=0, *y=0, *z=0;
    float tr[6], *a=0, zmin, zmax;

    RLABCODE_PROCESS_BTREE_ENTRY_MD(e1,node,"x",x,MDR,SIZE,1,NULL);
    if (!x)
      goto _exit_pgprintf;
    nr = SIZE(x);

    RLABCODE_PROCESS_BTREE_ENTRY_MD(e1,node,"y",y,MDR,SIZE,1,NULL);
    if (!y)
      goto _exit_pgprintf;
    nc = SIZE(y);

    RLABCODE_PROCESS_BTREE_ENTRY_MD(e1,node,"z",z,MDR,SIZE,1,NULL);
    if (!z)
      goto _exit_pgprintf;
    if ((MNR(z) != nr) || (MNC(z)!=nc))
      goto _exit_pgprintf;

    // z-range provided?
    zmin = mdrV0(plfmt,1);
    zmax = mdrV0(plfmt,2);

    if (isnand(zmin) || isnand(zmax))
    {
      zmax = zmin = mdrV0(z,0);
      for (i=1;i<nr*nc;i++)
      {
        if (zmax < mdrV0(z,i))
        {
          zmax = mdrV0(z,i);
        }
        else if (zmin > mdrV0(z,i))
        {
          zmin = mdrV0(z,i);
        }
      }
    }
    if (zmin == zmax)
    {
      zmin -= 0.5;
      zmax += 0.5;
    }

    // equidistant mesh for x:
    tr[2] = 0;
    tr[1] = ABS(mdrV0(x,0) - mdrV0(x,1));
    if (mdrV0(x,0) < mdrV0(x,nr-1))
    {
      tr[0] = mdrV0(x,0) - 0.5 * tr[1];
    }
    else
    {
      tr[0] = mdrV0(x,nr-1) - 0.5 * tr[1];
    }

    // equidistant mesh for y:
    tr[4] = 0;
    tr[5] = ABS(mdrV0(y,0) - mdrV0(y,1));
    if (mdrV0(y,0) < mdrV0(y,nc-1))
    {
      tr[3] = mdrV0(y,0) - 0.5 * tr[5];
    }
    else
    {
      tr[3] = mdrV0(y,nc-1) - 0.5 * tr[5];
    }
    a = mdr_coerce_floatptr (z);

    if (mdiV0(plfmt,0)==3)
    {
      // pm3d grey,map,grad
      // check mode and determine parameters
      int palette_param[7];
      for (i=0;i<7;i++)
        palette_param[i] = mdiV0(plfmt,3+i);

      cpgrgbfor (a, nr, nc, 1, nr, 1, nc, zmax, zmin, tr, palette_param);
    }
    else if (mdiV0(plfmt,0)==4)
    {
      int intval=150, minint = 80;
      float ts=-1.0;
      if (nargs > 5)
      {
        e6 = bltin_get_ent (args[5]);
        if (ent_type(e6)==MATRIX_DENSE_STRING)
        {
          labels = ent_data(e6);
        }
      }
      if (nargs > 6)
      {
        e7 = bltin_get_ent (args[6]);
        if (ent_type(e7)==MATRIX_DENSE_REAL)
        {
          ts = class_double(e7);
        }
      }
      if (nargs > 7)
      {
        e8 = bltin_get_ent (args[7]);
        if (ent_type(e8)==MATRIX_DENSE_REAL)
        {
          MDR * vals = ent_data(e8);
          if (SIZE(vals)==2)
          {
            intval = mdiV0(vals,0);
            minint = mdiV0(vals,1);
          }
        }
      }

      // col_idx is used for passing the zlevels
      // contour
      int contour_param[6];
      for (i=0;i<6;i++)
        contour_param[i] = mdiV0(plfmt,4+i);
      int zstp = mdrV0(plfmt,3);
      float *c = GC_MALLOC((zstp+2) * sizeof(float));
      c[zstp] = zmin;
      c[zstp+1] = zmax;
      if (col_idx)
      {
        for (i=0; i<zstp; i++)
        {
          c[i] = mdrV0(col_idx,i);
        }
      }
      else
      {
        // define our own contour
        float z_slope = (zmax - zmin)/(zstp - 1.0);
        c[0] = zmin;
        for (i=1; i<zstp; i++)
        {
          c[i] = zmin + z_slope*i;
        }
      }

      cpgcont(a, nr, nc, 1, nr, 1, nc, c, -zstp, tr, contour_param);

      // do we label?
      if (labels)
      {
        char *lx=0;
        char lab[33]={'\0'};
        for (i=0; i<=zstp; i++)
        {
          lx = MdsV0(labels, MIN(i,SIZE(labels)-1));
          if ((isvalidstring(lx)>0)&&(isvalidstring(lx)<33))
          {
            sprintf(lab, lx, c[i]);
            cpgconl(a, nr, nc, 1, nr, 1, nc, c, -zstp, tr, contour_param, i+1, lab, ts, intval, minint);
          }
        }
      }

      if (c)
        GC_FREE(c);
    }

    if (a)
      GC_FREE (a);

    goto _exit_pgprintf;
  }

  //
  // plot data matrix
  //
  if (ent_type(e1) != MATRIX_DENSE_REAL)
  {
    fprintf(rlab_stderr, THIS_SOLVER ": " RLAB_ERROR_ARG1_MDR_MATRIX "\n");
    goto _exit_pgprintf;
  }
  pldata = ent_data(e1);
  ndata = MNR(pldata);
  if (ndata < 0)
    goto _exit_pgprintf;

  int lpb = mdiV0(plfmt, RLABPLUS_PGPLOT_IDX_LINEPOINT);
  int lt  = mdiV0(plfmt, RLABPLUS_PGPLOT_IDX_LINETYPE);
  int lc  = mdiV0(plfmt, RLABPLUS_PGPLOT_IDX_LINECOLOR);
  double lw  = mdrV0(plfmt, RLABPLUS_PGPLOT_IDX_LINEWIDTH);
  int pt  = mdiV0(plfmt, RLABPLUS_PGPLOT_IDX_POINTTYPE);
  int pc  = mdiV0(plfmt, RLABPLUS_PGPLOT_IDX_POINTCOLOR);
  double ps  = mdrV0(plfmt, RLABPLUS_PGPLOT_IDX_POINTSIZE);
  int err = mdiV0(plfmt, 7);
  int err_xy = mdiV0(plfmt, 8);
  int err_lp = mdiV0(plfmt, 9);

  if (!col_idx)
  {
    clean_col_idx = 1;
    col_idx = mdr_Create(1,MNC(pldata));
    for (i=0; i<MNC(pldata); i++)
      MdrV0(col_idx,i) = i+1;
  }
  else
  {
    for (i=0; i<SIZE(col_idx); i++)
    {
      if (MdrV0(col_idx,i) > MNC(pldata))
      {
        fprintf(rlab_stderr, THIS_SOLVER ": " RLAB_ERROR_ARG5_MDR_VECTOR "\n");
        fprintf(rlab_stderr, THIS_SOLVER ": column indices cannot exceed number of columns of the data matrix!\n");
        goto _exit_pgprintf;
      }
    }
  }

  x = GC_MALLOC(ndata * sizeof(float));
  y = GC_MALLOC(ndata * sizeof(float));
  if (SIZE(col_idx)==1)
  {
    if (logx)
    {
      for (i=0; i<ndata; i++)
        x[i] = log10(i + 1);
    }
    else
    {
      for (i=0; i<ndata; i++)
        x[i] = i + 1;
    }
    if (logy)
    {
      for (i=0; i<ndata; i++)
        y[i] = log10(mdr0(pldata,i,mdiV0(col_idx,0)-1));
    }
    else
      for (i=0; i<ndata; i++)
        y[i] = mdr0(pldata,i,mdiV0(col_idx,0)-1);
  }
  else
  {
    if (logx)
    {
      for (i=0; i<ndata; i++)
        x[i] = log10(mdr0(pldata,i,mdiV0(col_idx,0)-1));
    }
    else
      for (i=0; i<ndata; i++)
        x[i] = mdr0(pldata,i,mdiV0(col_idx,0)-1);

    if (logy)
    {
      for (i=0; i<ndata; i++)
        y[i] = log10(mdr0(pldata,i,mdiV0(col_idx,1)-1));
    }
    else
      for (i=0; i<ndata; i++)
        y[i] = mdr0(pldata,i,mdiV0(col_idx,1)-1);
  }

  //
  // plot line first
  //
  if ((lpb==0) || (lpb==2) || ((err==1)&&(err_lp==0)))
  {
    if ((err==0) && (SIZE(col_idx)!=2))
    {
      fprintf(rlab_stderr, THIS_SOLVER ": " RLAB_ERROR_ARG1_MDR_MATRIX "\n");
      fprintf(rlab_stderr, THIS_SOLVER ": 1- or 2-column real matrix expected but not provided. Cannot continue!\n");
      goto _exit_pgprintf;
    }
    else if ((err==1) && (SIZE(col_idx)!=4) && (SIZE(col_idx)!=6))
    {
      fprintf(rlab_stderr, THIS_SOLVER ": " RLAB_ERROR_ARG1_MDR_MATRIX "\n");
      fprintf(rlab_stderr, THIS_SOLVER ": 4- or 6-column real matrix expected but not provided. Cannot continue!\n");
      goto _exit_pgprintf;
    }

    if (lt >=0 )
      cpgsls (((lt-1) % 5) + 1);  // set line style

    if (lc >=0 )
      cpgsci(lc); // set line color

    if (lw>0)
    {
      lw = lw > 201 ? 201 : lw;
      cpgslw (lw);// set line width
    }

    cpgline(ndata, (float *) x, (float *) y);

  }

  //
  // plot points next
  //
  if ((lpb==1) || (lpb==2)  || ((err==1)&&(err_lp==1)))
  {
    if ((err==0) && (SIZE(col_idx)!=2))
    {
      fprintf(rlab_stderr, THIS_SOLVER ": " RLAB_ERROR_ARG1_MDR_MATRIX "\n");
      fprintf(rlab_stderr, THIS_SOLVER ": 1- or 2-column real matrix expected but not provided. Cannot continue!\n");
      goto _exit_pgprintf;
    }
    else if ((err==1) && (SIZE(col_idx)!=4) && (SIZE(col_idx)!=6))
    {
      fprintf(rlab_stderr, THIS_SOLVER ": " RLAB_ERROR_ARG1_MDR_MATRIX "\n");
      fprintf(rlab_stderr, THIS_SOLVER ": 4- or 6-column real matrix expected but not provided. Cannot continue!\n");
      goto _exit_pgprintf;
    }

    if (pc >=0 )
      cpgsci(pc); // set symbol color

    if (ps >=0 )
      cpgsch(ps); // set symbol size

    cpgpt(ndata, (float *) x, (float *) y, pt);
  }

  //
  // plot error bars/lines if directed
  //
  if (err)
  {
    float x1, y1, xmin, xmax, ymin, ymax, t=1.0;
    for (i=0; i<ndata; i++)
    {
      if ((err_xy==1)||(err_xy==2))
      {
        xmin = mdr0(pldata,i,mdiV0(col_idx,2)-1);
        xmax = mdr0(pldata,i,mdiV0(col_idx,3)-1);
        y1 = mdr0(pldata,i,mdiV0(col_idx,1)-1);
        if (logx)
        {
          xmin = log10(xmin);
          xmax = log10(xmax);
        }
        if (logy)
        {
          y1 = log10(y1);
        }

        cpgerrx (1, (float *) &xmin, (float *) &xmax, (float *) &y1, t);
      }
      if (err_xy==3)
      {
        x1   = mdr0(pldata,i,mdiV0(col_idx,0)-1);
        ymin = mdr0(pldata,i,mdiV0(col_idx,2)-1);
        ymax = mdr0(pldata,i,mdiV0(col_idx,3)-1);
        if (logx)
        {
          x1 = log10(x1);
        }
        if (logy)
        {
          ymin = log10(ymin);
          ymax = log10(ymax);
        }
        cpgerry(1, (float *) &x1, (float *) &ymin, (float *) &ymax, t);
      }
      else if (err_xy==2)
      {
        x1 = mdr0(pldata,i,mdiV0(col_idx,0)-1);
        ymin = mdr0(pldata,i,mdiV0(col_idx,4)-1);
        ymax = mdr0(pldata,i,mdiV0(col_idx,5)-1);
        if (logx)
        {
          x1 = log10(x1);
        }
        if (logy)
        {
          ymin = log10(ymin);
          ymax = log10(ymax);
        }
        cpgerry(1, (float *) &x1, (float *) &ymin, (float *) &ymax, t);
      }
    }
  }

_exit_pgprintf:

  // reset values
  cpgslw (old_lw);  // line width
  cpgsci (old_lc);  // line color
  cpgsls (old_lt);  // line style
  cpgsch (old_ps);  // symbol size

  if (clean_col_idx)
    mdr_Destroy(col_idx);

  ent_Clean (e1);
  ent_Clean (e2);
  ent_Clean (e3);
  ent_Clean (e4);
  ent_Clean (e5);
  ent_Clean (e6);
  ent_Clean (e7);
  ent_Clean (e8);

  if (x)
    GC_FREE(x);
  if (y)
    GC_FREE(y);

  return ent_Create_Rlab_Success();
}


/* **************************************************************
 * Open a graphics device.
 * INTEGER FUNCTION PGBEG (UNIT, FILE, NXSUB, NYSUB)
 *     INTEGER       UNIT=0
 *     CHARACTER*(*) FILE
 *     INTEGER       NXSUB, NYSUB
 * ************************************************************** */
#undef  THIS_SOLVER
#define THIS_SOLVER "pgbeg"
Ent * _pg_cpgbeg (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *e3=0;
  int nxsub=1, nysub=1;
  char *file=0;

  if ((nargs != 1)&&(nargs != 3))
  {
    printf ("pgbeg: 1 or 3 arguments are required\n");
    goto _exit_pgbeg;
  }

  e1 = bltin_get_ent (args[0]);
  file = class_char_pointer (e1);
  if (isvalidstring(file)<1)
    goto _exit_pgbeg;

  if (nargs == 3)
  {
    e2 = bltin_get_ent (args[1]);
    e3 = bltin_get_ent (args[2]);
    if (ent_type(e2)==MATRIX_DENSE_REAL && ent_type(e3)==MATRIX_DENSE_REAL)
    {
      nxsub = class_int (e2);
      nysub = class_int (e3);
    }
  }

  cpgbeg (0, file, nxsub, nysub);

_exit_pgbeg:

  ent_Clean (e1);
  ent_Clean (e2);
  ent_Clean (e3);
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
  char *device=0;
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
#undef  THIS_SOLVER
#define THIS_SOLVER "pgsvp"
Ent * _pg_cpgsvp (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *e3=0, *e4=0;
  float xleft=-1.0, xright=-1.0, ybot=-1.0, ytop=-1.0;

  if ((nargs != 1)&&(nargs != 4))
  {
    rerror (THIS_SOLVER ": " RLAB_ERROR_ONE_OR_FOUR_ARG_REQUIRED "\n");
  }

  if (nargs == 1)
  {
    e1 = bltin_get_ent (args[0]);
    if (ent_type(e1)==MATRIX_DENSE_REAL)
    {
      MDR *x = ent_data(e1);
      if (SIZE(x)==4)
      {
        xleft  = (float) mdrV0(x,0);
        xright = (float) mdrV0(x,1);
        ybot   = (float) mdrV0(x,2);
        ytop   = (float) mdrV0(x,3);
      }
    }
  }
  else
  {
    e1 = bltin_get_ent (args[0]);
    if (ent_type(e1)==MATRIX_DENSE_REAL)
      xleft = (float) class_double (e1);

    e2 = bltin_get_ent (args[1]);
    if (ent_type(e2)==MATRIX_DENSE_REAL)
      xright = (float) class_double (e2);

    e3 = bltin_get_ent (args[2]);
    if (ent_type(e2)==MATRIX_DENSE_REAL)
      ybot = (float) class_double (e3);

    e4 = bltin_get_ent (args[3]);
    if (ent_type(e2)==MATRIX_DENSE_REAL)
      ytop = (float) class_double (e4);
  }

  if ((xleft == -1.0)||(xright == -1.0)||(ybot == -1.0)||(ytop == -1.0))
    cpgvstd ();
  else
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

Ent * _pg_cpgvstd (int nargs, Datum args[])
{
  if (nargs != 0)
    rerror ("pgvstd: no arguments allowed");

  cpgvstd ();

  return ent_Create_Rlab_Success();
}

/* **************************************************************
 * set window
 *    SUBROUTINE PGVSIZ (XL, XR, YB, YT)
 *       REAL XL, XR, YB, YT
 * ************************************************************** */

Ent * _pg_cpgvsiz (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *e3=0, *e4=0;
  float x1, x2, y1, y2;

  if (nargs != 4)
  {
    rerror ("_pgvsiz: 4 arguments required");
  }

  e1 = bltin_get_ent (args[0]);
  e2 = bltin_get_ent (args[1]);
  e3 = bltin_get_ent (args[2]);
  e4 = bltin_get_ent (args[3]);

  x1 = (float) class_double (e1);
  x2 = (float) class_double (e2);
  y1 = (float) class_double (e3);
  y2 = (float) class_double (e4);

  cpgvsiz (x1, x2, y1, y2);

  ent_Clean (e1);
  ent_Clean (e2);
  ent_Clean (e3);
  ent_Clean (e4);

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
Ent * _pg_cpgslct (int nargs, Datum args[])
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
#undef  THIS_SOLVER
#define THIS_SOLVER "pgline"
Ent * _pg_cpgline (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *e3=0;
  MDR *x1=0, *x2=0, *x3=0;
  float x[2], y[2];
  int i, n;

  if ((nargs != 1)&&(nargs != 3))
  {
    rerror (THIS_SOLVER ": 1 or 3 arguments are required");
  }

  e1 = bltin_get_ent (args[0]);
  if (ent_type(e1)!=MATRIX_DENSE_REAL)
    goto _exit_pgline;
  x1 = ent_data (e1);
  if (nargs == 1)
  {
    if (MNC(x1)!=2)
      goto _exit_pgline;
    if (MNR(x1)<2)
      goto _exit_pgline;

    x[1] = mdr0(x1,0,0);
    y[1] = mdr0(x1,0,1);
    for (i=1; i<MNR(x1); i++)
    {
      x[0] = x[1];
      x[1] = mdr0(x1,i,0);
      y[0] = y[1];
      y[1] = mdr0(x1,i,1);
      cpgline (2, x, y);
    }
  }
  else
  {
    if (SIZE(x1)!=1)
      goto _exit_pgline;
    n = mdrV0(x1,0);
    if (n<2)
      goto _exit_pgline;

    e2 = bltin_get_ent (args[1]);
    if (ent_type(e2)!=MATRIX_DENSE_REAL)
      goto _exit_pgline;
    x2 = ent_data(e2);

    e3 = bltin_get_ent (args[2]);
    if (ent_type(e3)!=MATRIX_DENSE_REAL)
      goto _exit_pgline;
    x3 = ent_data(e3);

    if (SIZE(x2)!=SIZE(x3))
      goto _exit_pgline;
    if (SIZE(x2)!=n)
      goto _exit_pgline;

    x[1] = mdrV0(x2,0);
    y[1] = mdrV0(x3,0);
    for (i=1; i<n; i++)
    {
      x[0] = x[1];
      x[1] = mdrV0(x2,i);
      y[0] = y[1];
      y[1] = mdrV0(x3,i);
      cpgline (2, x, y);
    }
  }

_exit_pgline:

  ent_Clean (e1);
  ent_Clean (e2);
  ent_Clean (e3);
  return ent_Create_Rlab_Success();
}

/* **************************************************************
 * draw several graph markers
 *    SUBROUTINE PGPT (N, XPTS, YPTS, SYMBOL)
 *       INTEGER N
 *       REAL XPTS(*), YPTS(*)
 *       INTEGER SYMBOL
 * ************************************************************** */
#undef  THIS_SOLVER
#define THIS_SOLVER "pgpt"
Ent * _pg_cpgpt (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *e3=0, *e4=0;
  MDR *x1=0, *x2=0, *x3=0;
  float x, y;
  int i, n, symbol;

  if ((nargs != 2)&&(nargs != 4))
  {
    rerror (THIS_SOLVER ": 1 or 3 arguments are required");
  }

  e1 = bltin_get_ent (args[0]);
  if (ent_type(e1)!=MATRIX_DENSE_REAL)
    goto _exit_pgpt;
  x1 = ent_data (e1);
  if (nargs == 2)
  {
    if (MNC(x1)!=2)
      goto _exit_pgpt;

    e2 = bltin_get_ent (args[1]);
    if (ent_type(e2)!=MATRIX_DENSE_REAL)
      goto _exit_pgpt;
    symbol = (int) class_double (e2);

    for (i=0; i<MNR(x1); i++)
    {
      x = mdr0(x1,i,0);
      y = mdr0(x1,i,1);
      cpgpt (1, &x, &y, symbol);
    }
  }
  else
  {
    if (SIZE(x1)!=1)
      goto _exit_pgpt;
    n = mdrV0(x1,0);

    e2 = bltin_get_ent (args[1]);
    if (ent_type(e2)!=MATRIX_DENSE_REAL)
      goto _exit_pgpt;
    x2 = ent_data(e2);

    e3 = bltin_get_ent (args[2]);
    if (ent_type(e3)!=MATRIX_DENSE_REAL)
      goto _exit_pgpt;
    x3 = ent_data(e3);

    if (SIZE(x2)!=SIZE(x3))
      goto _exit_pgpt;
    if (SIZE(x2)!=n)
      goto _exit_pgpt;

    e4 = bltin_get_ent (args[3]);
    if (ent_type(e4)!=MATRIX_DENSE_REAL)
      goto _exit_pgpt;
    symbol = (int) class_double (e4);

    for (i=0; i<n; i++)
    {
      x = mdrV0(x2,i);
      y = mdrV0(x3,i);
      cpgpt (1, &x, &y, symbol);
    }
  }

_exit_pgpt:

  ent_Clean (e1);
  ent_Clean (e2);
  ent_Clean (e3);
  ent_Clean (e4);
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
#undef THIS_SOLVER
#define THIS_SOLVER "pggrey"
Ent * _pg_cpggray (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *e3=0, *e4=0;
  MDR *x=0, *y=0, *z=0, *zmm=0;
  int nr, nc, i;
  float *a=0, zmin, zmax, tr[6]={0};

  if ((nargs != 3)&& (nargs != 4))
  {
    rerror (THIS_SOLVER ": 3 or 4 arguments required!\n");
  }

  // x
  e1 = bltin_get_ent (args[0]);
  if (ent_type(e1)!=MATRIX_DENSE_REAL)
    goto _exit_cpggray;
  x = ent_data(e1);
  nr = SIZE(x);

  e2 = bltin_get_ent (args[1]);
  if (ent_type(e2)!=MATRIX_DENSE_REAL)
    goto _exit_cpggray;
  y = ent_data(e2);
  nc = SIZE(y);

  e3 = bltin_get_ent (args[2]);
  if (ent_type(e3)!=MATRIX_DENSE_REAL)
    goto _exit_cpggray;
  z = ent_data(e3);
  if (nr != MNR(z))
    goto _exit_cpggray;
  nc = MNC(z);
  if (nc != MNC(z))
    goto _exit_cpggray;

  if (nargs>3)
  {
    e4 = bltin_get_ent (args[3]);
    if (ent_type(e4)==MATRIX_DENSE_REAL)
    {
      zmm = ent_data(e4);
      if (SIZE(zmm)==2)
      {
        zmin = mdrV0(zmm,0);
        zmax = mdrV0(zmm,1);
      }
      else
        zmm=0;
    }
  }

  if (!zmm)
  {
    zmax = zmin = mdrV0(z,0);
    for (i=1;i<nr*nc;i++)
    {
      if (zmax < mdrV0(z,i))
      {
        zmax = mdrV0(z,i);
      }
      else if (zmin > mdrV0(z,i))
      {
        zmin = mdrV0(z,i);
      }
    }
  }

  // equidistant mesh for x:
  tr[2] = 0;
  tr[1] = ABS(mdrV0(x,0) - mdrV0(x,1));
  if (mdrV0(x,0) < mdrV0(x,nr))
  {
    tr[0] = mdrV0(x,0) - tr[1];
  }
  else
  {
    tr[0] = mdrV0(x,nr) - tr[1];
  }

  // equidistant mesh for y:
  tr[4] = 0;
  tr[5] = ABS(mdrV0(y,0) - mdrV0(y,1));
  if (mdrV0(y,0) < mdrV0(y,nc))
  {
    tr[3] = mdrV0(y,0) - tr[5];
  }
  else
  {
    tr[3] = mdrV0(y,nc) - tr[5];
  }
  a = mdr_coerce_floatptr (z);

  cpggray (a, nr, nc, 1, nr, 1, nc, zmax, zmin, tr);

_exit_cpggray:

  if (a)
    GC_FREE (a);

  ent_Clean (e1);
  ent_Clean (e2);
  ent_Clean (e3);
  ent_Clean (e4);

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

  cpgconl (a, idim, jdim, i1, i2, j1, j2, &c, 1, tr, (int *)NULL, 1, label, 1.0, intval, minint);

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

  cpgcont (a, idim, jdim, i1, i2, j1, j2, c, nc, tr, NULL);

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
#undef THIS_SOLVER
#define THIS_SOLVER "pgbox"
Ent * _pg_cpgbox (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *e3=0, *e4=0, *e5=0, *e6;
  char *xopt, *yopt;
  float xtick[4]={0,-1,-1,-1}, ytick[4]={0,-1,-1,-1};
  int i, nxsub[4]={0,-1,-1,-1}, nysub[4]={0,-1,-1,-1};

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
  if (ent_type(e2)==MATRIX_DENSE_REAL)
  {
    MDR *d = ent_data(e2);
    if (SIZE(d)>0)
    {
      for (i=0;(i<SIZE(d))&&(i<4);i++)
      {
        xtick[i] = mdrV0(d,i);
      }
    }
  }
  if (ent_type(e3)==MATRIX_DENSE_REAL)
  {
    MDR *d = ent_data(e3);
    if (SIZE(d)>0)
    {
      for (i=0;(i<SIZE(d))&&(i<4);i++)
      {
        nxsub[i] = mdiV0(d,i);
      }
    }
  }

  yopt = class_char_pointer (e4);
  if (ent_type(e5)==MATRIX_DENSE_REAL)
  {
    MDR *d = ent_data(e5);
    if (SIZE(d)>0)
    {
      for (i=0;(i<SIZE(d))&&(i<4);i++)
      {
        ytick[i] = mdrV0(d,i);
      }
    }
  }
  if (ent_type(e6)==MATRIX_DENSE_REAL)
  {
    MDR *d = ent_data(e6);
    if (SIZE(d)>0)
    {
      for (i=0;(i<SIZE(d))&&(i<4);i++)
      {
        nysub[i] = mdiV0(d,i);
      }
    }
  }

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

Ent * _pg_color (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *e3=0, *e4=0;
  int ci;
  float cr, cg, cb;
  MDR *rval=0;

  if (nargs != 1 && nargs != 2 && nargs != 4)
  {
    printf ("_pgcolor: 1,2 or 4 arguments are required");
    goto _exit_pg_color;
  }

  e1 = bltin_get_ent (args[0]);
  ci = class_int (e1);
  rval = mdi_Create(1,3);

  if (nargs == 1)
  {
    cpgqcr (ci, &cr, &cg, &cb);
  }
  else if (nargs == 2)
  {
    e2 = bltin_get_ent (args[1]);
    if (ent_type(e2)==MATRIX_DENSE_REAL)
    {
      MDR *rgb = ent_data(e2);
      if (SIZE(rgb)==3)
      {
        cr = mdrV0(rgb,0) / 255.0;
        cg = mdrV0(rgb,1) / 255.0;
        cb = mdrV0(rgb,2) / 255.0;
      }
    }
  }
  else if (nargs == 4)
  {
    e2 = bltin_get_ent (args[1]);
    e3 = bltin_get_ent (args[2]);
    e4 = bltin_get_ent (args[3]);
    cr = (float) class_double (e2) / 255.0;
    cg = (float) class_double (e3) / 255.0;
    cb = (float) class_double (e4) / 255.0;
  }

  cpgscr (ci, cr, cg, cb);

  MdiV0(rval,0) = cr * 255.0;
  MdiV0(rval,1) = cg * 255.0;
  MdiV0(rval,2) = cb * 255.0;

_exit_pg_color:

  ent_Clean (e1);
  ent_Clean (e2);
  ent_Clean (e3);
  ent_Clean (e4);
  return ent_Assign_Rlab_MDR(rval);
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

  cpgsls (((ls-1) % 5) + 1);

  ent_Clean (e1);

  return ent_Create_Rlab_Success();
}

/* **************************************************************
 * set line width
 *    SUBROUTINE PGSLW (LW)
 *         INTEGER  LW
 * ************************************************************** */

Ent * _pg_cpgslw (int nargs, Datum args[])
{
  Ent *e1=0;
  int lw;

  if (nargs != 1)
  {
    rerror ("pgslw: 1 argument allowed");
  }

  e1 = bltin_get_ent (args[0]);
  lw = (int) class_double (e1);
  lw = lw > 201 ? 201 : lw;

  cpgslw (lw);

  ent_Clean (e1);

  return ent_Create_Rlab_Success();
}

/* **************************************************************
 * close all open graphics devices
 *    SUBROUTINE PGEND
 * ************************************************************** */
Ent * _pg_cpgend (int nargs, Datum args[])
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
Ent * _pg_cpgclos (int nargs, Datum args[])
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

Ent * _pg_cpgpage (int nargs, Datum args[])
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
Ent * _pg_cpgeras (int nargs, Datum args[])
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
Ent * _pg_cpgetxt (int nargs, Datum args[])
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
Ent * _pg_cpgtext (int nargs, Datum args[])
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
#undef THIS_SOLVER
#define THIS_SOLVER "pgptxt"
Ent * _pg_cpgptxt (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *e3=0, *e4=0, *e5=0;
  float x, y, angle, fjust;
  char *text;

  if (nargs != 5)
  {
    rerror (THIS_SOLVER ": 5 arguments are required\n");
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

//   fprintf(stderr,THIS_SOLVER ": [x,y]=%f,%f,%s\n",x,y,text);

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
Ent * _pg_cpgmtxt (int nargs, Datum args[])
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
Ent * _pg_cpgscf (int nargs, Datum args[])
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
Ent * _pg_cpgsch (int nargs, Datum args[])
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
Ent * _pg_cpgask (int nargs, Datum args[])
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
Ent * _pg_cpglab (int nargs, Datum args[])
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
Ent * _pg_cpgupdt (int nargs, Datum args[])
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
Ent * _pg_cpgsci (int nargs, Datum args[])
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
Ent * _pg_cpgscir (int nargs, Datum args[])
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
Ent * _pg_cpgqci (int nargs, Datum args[])
{
  int ci;
  MDR *color;

  if (nargs != 0)
  {
    rerror ("pqsci: 0 arguments allowed");
  }

  cpgqci (&ci);

  color = mdr_CreateScalar( ci );

  return ent_Assign_Rlab_MDR(color);
}

/* **************************************************************
 * inquire color index range
 *     SUBROUTINE PGQCIR(ICILO, ICIHI)
 *        INTEGER   ICILO, ICIHI
 * ************************************************************** */
Ent * _pg_cpgqcir (int nargs, Datum args[])
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
Ent * _pg_cpgbbuf (int nargs, Datum args[])
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
Ent * _pg_cpgebuf (int nargs, Datum args[])
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
Ent * _pg_cpgqcol (int nargs, Datum args[])
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
Ent * _pg_cpgmove (int nargs, Datum args[])
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
Ent * _pg_cpgdraw (int nargs, Datum args[])
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
Ent * _pg_cpgsfs (int nargs, Datum args[])
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

float * mdr_coerce_floatptr (MDR * m)
{
  float *ret = 0;
  int i = 0;
  int size = 0;

  /* Compute the size of the new float array. */
  size = SIZE(m);

  /* Check for the possibility of a zero (empty) input. */
  if (size<1)
  {
    return (NULL);
  }

  /* Get memory for the new array. */
  ret = (float *) GC_malloc (size * sizeof (float));

  /* Convert the double precision array into float. */
  for (i = 0; i < size; i++)
  {
    ret[i] = (float) mdrV0_safe(m, i);
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


/* **************************************************************
 * inquire color index range
 *    SUBROUTINE PGQVP (UNITS, X1, X2, Y1, Y2)
 *    INTEGER UNITS
 *    REAL    X1, X2, Y1, Y2
 * ************************************************************** */
#undef  THIS_SOLVER
#define THIS_SOLVER "pgqvp"
Ent * _pg_cpgqvp (int nargs, Datum args[])
{
  Ent *e1=0;
  int unit;
  float x1,x2,y1,y2;

  MDR *rval=0;

  if (nargs != 1)
  {
    rerror (THIS_SOLVER ": " RLAB_ERROR_ARG1_MDR_INTEGER);
  }

  e1 = bltin_get_ent (args[0]);
  unit = class_int (e1);
  if (unit < 0)
    unit = 0;
  else if (unit > 3)
    unit = 3;

  cpgqvp (unit, &x1, &x2, &y1, &y2);

  rval = mdr_Create (1, 4);
  MdrV0 (rval, 0) = (double) x1;
  MdrV0 (rval, 1) = (double) x2;
  MdrV0 (rval, 2) = (double) y1;
  MdrV0 (rval, 3) = (double) y2;

  ent_Clean (e1);

  return ent_Assign_Rlab_MDR( rval );
}

/* **************************************************************
 * inquire color index range
 *    SUBROUTINE PGQWIN (X1, X2, Y1, Y2)
 *    REAL X1, X2, Y1, Y2
 * ************************************************************** */
#undef  THIS_SOLVER
#define THIS_SOLVER "pgqwin"
Ent * _pg_pgqwin(int nargs, Datum args[])
{
  float x1,x2,y1,y2;

  MDR *rval=0;

  cpgqwin (&x1, &x2, &y1, &y2);

  rval = mdr_Create (1, 4);
  MdrV0 (rval, 0) = (double) x1;
  MdrV0 (rval, 1) = (double) x2;
  MdrV0 (rval, 2) = (double) y1;
  MdrV0 (rval, 3) = (double) y2;

  return ent_Assign_Rlab_MDR( rval );
}

/* **************************************************************
 * inquire color index range
 *    SUBROUTINE PGLEN (UNITS, STRING, XL, YL)
 *    REAL XL, YL
 *    INTEGER UNITS
 *    CHARACTER*(*) STRING
 * ************************************************************** */
#undef  THIS_SOLVER
#define THIS_SOLVER "pglen"
Ent * _pg_cpglen (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0;
  int unit,i;
  float x,y;
  MDS *s=0;
  MDR *rval=0;

  if (nargs != 2)
  {
    rerror (THIS_SOLVER ": " RLAB_ERROR_ARG1_MDR_INTEGER);
  }

  e1 = bltin_get_ent (args[0]);
  unit = class_int (e1);
  if (unit < 0)
    unit = 0;
  else if (unit > 3)
    unit = 3;

  e2 = bltin_get_ent (args[1]);
  if (ent_type(e2)!=MATRIX_DENSE_STRING)
    goto _exit_pglen;
  s = ent_data(e2);
  if (SIZE(s)<1)
    goto _exit_pglen;

  rval = mdr_Create(SIZE(s),2);
  for (i=0; i<SIZE(s); i++)
  {
    Mdr0(rval,i,0) = 0;
    Mdr0(rval,i,1) = 0;
    if (isvalidstring(MdsV0(s,i)))
    {
      cpglen (unit, MdsV0(s,i), &x, &y);
      Mdr0 (rval, i, 0) = (double) x;
      Mdr0 (rval, i, 1) = (double) y;
    }
  }
_exit_pglen:

  ent_Clean (e1);
  ent_Clean (e2);
  return ent_Assign_Rlab_MDR( rval );
}

