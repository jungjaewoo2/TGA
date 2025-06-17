/* r_plot.c */

/*
* This file contains the PLPLOT interface routines.
*/

/*  This file is a part of RLaB ("Our"-LaB)
  Copyright (C) 1995  Ian R. Searle
  Copyright (C) 2014-2016  Marijan Kostrunz

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

#ifdef HAVE_RLAB_PLPLOT

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
#include "mathl.h"
extern FILE *RLAB_STDERR_DS;


#include "rlab_solver_parameters_names.h"

#include <stdio.h>

#undef HAVE_CONFIG_H
#include <plplot.h>
#include <pdf.h>
#include <plstrm.h>

static int init_plot_device = 0;	/* Has the plot device been initialized */

// static PLFLT **convert_matrix_to_array (MDR * m);

/* **************************************************************
* PLPLOT related functions and interface(s)
* ************************************************************** */
#define RLABPLUS_PLPLOT_IDX_LINEPOINT 0
#define RLABPLUS_PLPLOT_IDX_LINETYPE  1
#define RLABPLUS_PLPLOT_IDX_LINECOLOR 2
#define RLABPLUS_PLPLOT_IDX_LINEWIDTH 3
#define RLABPLUS_PLPLOT_IDX_POINTTYPE  4
#define RLABPLUS_PLPLOT_IDX_POINTCOLOR 5
#define RLABPLUS_PLPLOT_IDX_POINTSIZE  6

//
// plprint (data, format)
//
#undef THIS_SOLVER
#define THIS_SOLVER "_plprintf"
Ent * ent_plplot_printf (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *e3=0, *e4=0, *e5=0;
  int i, ndata, logx=0, logy=0, clean_col_idx=0;
  FILE *rlab_stderr = (!RLAB_STDERR_DS) ? stderr : RLAB_STDERR_DS;
  double *x=0, *y=0;

  MDR *plfmt=0, *pldata=0, *col_idx=0;

  if ((nargs < 2) || (nargs > 5))
  {
    fprintf(rlab_stderr, THIS_SOLVER ": " RLAB_ERROR_TWO_TO_FIVE_ARG_REQUIRED "\n");
    goto _exit_plprint;
  }

  //
  // data set to be plotted
  //
  e1 = bltin_get_ent(args[0]);
  if (ent_type(e1) != MATRIX_DENSE_REAL)
  {
    fprintf(rlab_stderr, THIS_SOLVER ": " RLAB_ERROR_ARG1_MDR_MATRIX "\n");
    goto _exit_plprint;
  }
  pldata = ent_data(e1);
  ndata = MNR(pldata);
  if (ndata < 0)
    goto _exit_plprint;

  //
  // how to plot it
  //
  e2 = bltin_get_ent(args[1]);
  if (ent_type(e2) != MATRIX_DENSE_REAL)
  {
    fprintf(rlab_stderr, THIS_SOLVER ": " RLAB_ERROR_ARG2_MDR_VECTOR "\n");
    goto _exit_plprint;
  }
  plfmt = ent_data(e2);
  if (SIZE(plfmt)!=10)
  {
    fprintf(rlab_stderr, THIS_SOLVER ": " RLAB_ERROR_ARG2_MDR_VECTOR "\n");
    fprintf(rlab_stderr, THIS_SOLVER ": 10-component vector expected but %i-component was provided.\n"
        THIS_SOLVER ": Cannot continue!\n", MNC(plfmt));
    goto _exit_plprint;
  }
  int lpb = mdiV0(plfmt, RLABPLUS_PLPLOT_IDX_LINEPOINT);
  int lt  = mdiV0(plfmt, RLABPLUS_PLPLOT_IDX_LINETYPE);
  int lc  = mdiV0(plfmt, RLABPLUS_PLPLOT_IDX_LINECOLOR);
  double lw  = mdrV0(plfmt, RLABPLUS_PLPLOT_IDX_LINEWIDTH);
  int pt  = mdiV0(plfmt, RLABPLUS_PLPLOT_IDX_POINTTYPE);
  int pc  = mdiV0(plfmt, RLABPLUS_PLPLOT_IDX_POINTCOLOR);
  double ps  = mdrV0(plfmt, RLABPLUS_PLPLOT_IDX_POINTSIZE);
  int err = mdiV0(plfmt, 7);
  int err_xy = mdiV0(plfmt, 8);
  int err_lp = mdiV0(plfmt, 9);

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
      if (MdrV0(col_idx,i) > MNC(pldata))
    {
      fprintf(rlab_stderr, THIS_SOLVER ": " RLAB_ERROR_ARG5_MDR_VECTOR "\n");
      fprintf(rlab_stderr, THIS_SOLVER ": column indices cannot exceed number of columns of the data matrix!\n");
      goto _exit_plprint;
    }
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
      goto _exit_plprint;
    }
    else if ((err==1) && (SIZE(col_idx)!=4) && (SIZE(col_idx)!=6))
    {
      fprintf(rlab_stderr, THIS_SOLVER ": " RLAB_ERROR_ARG1_MDR_MATRIX "\n");
      fprintf(rlab_stderr, THIS_SOLVER ": 4- or 6-column real matrix expected but not provided. Cannot continue!\n");
      goto _exit_plprint;
    }

    // with lines, with both
    if (lt >=0 )
      pllsty(lt);   // set line style
    if (lc >=0 )
      plcol0(lc);   // set line color
    if (lw>0)
#ifdef plwidth
      plwidth (lw);// set line width
#else
      plwid (lw);
#endif

    if (MD_TYPE_DOUBLE(pldata))
    {
      if (SIZE(col_idx)==1)
      {
        x = GC_MALLOC(ndata * sizeof(double));
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
          y = GC_MALLOC(ndata * sizeof(double));
          for (i=0; i<ndata; i++)
            y[i] = log10(Mdr0(pldata,i,mdiV0(col_idx,0)-1));
        }
        else
          y = &MdrV0(pldata,mdiV0(col_idx,0)-1);
      }
      else
      {
        if (logx)
        {
          x = GC_MALLOC(ndata * sizeof(double));
          for (i=0; i<ndata; i++)
            x[i] = log10(Mdr0(pldata,i,mdiV0(col_idx,0)-1));
        }
        else
          x = &Mdr0(pldata,0,mdiV0(col_idx,0)-1);

        if (logy)
        {
          y = GC_MALLOC(ndata * sizeof(double));
          for (i=0; i<ndata; i++)
            y[i] = log10(Mdr0(pldata,i,mdiV0(col_idx,1)-1));
        }
        else
          y = &Mdr0(pldata,0,mdiV0(col_idx,1)-1);
      }

      plline(ndata, (PLFLT *) x, (PLFLT *) y);

      if (logx || MNC(pldata)==1)
        GC_FREE (x);
      if (logy)
        GC_FREE (y);
    }
    else if (MD_TYPE_INT32(pldata))
    {
      double x1, y1, x2, y2;
      if (SIZE(col_idx)==1)
      {
        if (ndata > 1)
        {
          for (i=0; i<ndata-1; i++)
          {
            x1 = i+1;
            y1 = Mdi0(pldata,i,mdiV0(col_idx,0)-1);
            x2 = x1 + 1;
            y2 = Mdi0(pldata,i+1,mdiV0(col_idx,0)-1);
            if (logx)
            {
              x1 = log10(x1);
              x2 = log10(x2);
            }
            if (logy)
            {
              y1 = log10(y1);
              y2 = log10(y2);
            }
            pljoin(x1,y1,x2,y2);
          }
        }
        else
        {
          x1 = 1;
          y1 = Mdi0(pldata,0,mdiV0(col_idx,0)-1);
          x2 = 1;
          y2 = Mdi0(pldata,0,mdiV0(col_idx,0)-1);
          if (logx)
          {
            x1 = log10(x1);
            x2 = log10(x2);
          }
          if (logy)
          {
            y1 = log10(y1);
            y2 = log10(y2);
          }
          pljoin(x1,y1,x2,y2);
        }
      }
      else
      {
        if (ndata > 1)
        {
          for (i=0; i<ndata-1; i++)
          {
            x1 = Mdi0(pldata,i,mdiV0(col_idx,0)-1);
            y1 = Mdi0(pldata,i,mdiV0(col_idx,1)-1);
            x2 = Mdi0(pldata,i+1,mdiV0(col_idx,0)-1);
            y2 = Mdi0(pldata,i+1,mdiV0(col_idx,1)-1);
            if (logx)
            {
              x1 = log10(x1);
              x2 = log10(x2);
            }
            if (logy)
            {
              y1 = log10(y1);
              y2 = log10(y2);
            }
            pljoin(x1,y1,x2,y2);
          }
        }
        else
        {
          x1 = Mdi0(pldata,0,mdiV0(col_idx,0)-1);
          y1 = Mdi0(pldata,0,mdiV0(col_idx,1)-1);
          x2 = Mdi0(pldata,0,mdiV0(col_idx,0)-1);
          y2 = Mdi0(pldata,0,mdiV0(col_idx,1)-1);
          if (logx)
          {
            x1 = log10(x1);
            x2 = log10(x2);
          }
          if (logy)
          {
            y1 = log10(y1);
            y2 = log10(y2);
          }
          pljoin(x1,y1,x2,y2);
        }
      }
    }
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
      goto _exit_plprint;
    }
    else if ((err==1) && (SIZE(col_idx)!=4) && (SIZE(col_idx)!=6))
    {
      fprintf(rlab_stderr, THIS_SOLVER ": " RLAB_ERROR_ARG1_MDR_MATRIX "\n");
      fprintf(rlab_stderr, THIS_SOLVER ": 4- or 6-column real matrix expected but not provided. Cannot continue!\n");
      goto _exit_plprint;
    }

    // with points, with both
    plssym(0,ps);  // set symbol size
    if (lc >=0 )
      plcol0(pc);     // set symbol color
    if (MD_TYPE_DOUBLE(pldata))
    {
      if (SIZE(col_idx)==1)
      {
        double *x = GC_MALLOC(ndata * sizeof(double));
        for (i=0; i<ndata; i++)
          x[i] = i + 1;
        plpoin(ndata, (PLFLT *) x, (PLFLT *) &Mdr0(pldata,0,mdiV0(col_idx,0)-1), pt);
        GC_FREE (x);
      }
      else
      {
        plpoin(ndata, (PLFLT *) &Mdr0(pldata,0,mdiV0(col_idx,0)-1), (PLFLT *) &Mdr0(pldata,0,mdiV0(col_idx,1)-1), pt);
      }
    }
    else if (MD_TYPE_INT32(pldata))
    {
      double x1, y1;
      if (SIZE(col_idx)==1)
      {
        for (i=0; i<ndata; i++)
        {
          x1 = i+1;
          y1 = Mdi0(pldata,i,mdiV0(col_idx,0)-1);
          if (logx)
          {
            x1 = log10(x1);
          }
          if (logy)
          {
            y1 = log10(y1);
          }
          plpoin(1, &x1, &y1, pt);
        }
      }
      else
      {
        for (i=0; i<ndata-1; i++)
        {
          x1 = Mdi0(pldata,i,mdiV0(col_idx,0)-1);
          y1 = Mdi0(pldata,i,mdiV0(col_idx,1)-1);
          if (logx)
          {
            x1 = log10(x1);
          }
          if (logy)
          {
            y1 = log10(y1);
          }
          plpoin(1, &x1, &y1, pt);
        }
      }
    }
  }

  //
  // plot error bars/lines if directed
  //
  if (err)
  {
    double x, y, xmin, xmax, ymin, ymax;
    for (i=0; i<ndata; i++)
    {
      if ((err_xy==1)||(err_xy==2))
      {
        xmin = mdr0(pldata,i,mdiV0(col_idx,2)-1);
        xmax = mdr0(pldata,i,mdiV0(col_idx,3)-1);
        y = mdr0(pldata,i,mdiV0(col_idx,1)-1);
        if (logx)
        {
          xmin = log10(xmin);
          xmax = log10(xmax);
        }
        if (logy)
        {
          y = log10(y);
        }

        plerrx(1, (PLFLT *) &xmin, (PLFLT *) &xmax, (PLFLT *) &y);
      }
      if (err_xy==3)
      {
        x    = mdr0(pldata,i,mdiV0(col_idx,0)-1);
        ymin = mdr0(pldata,i,mdiV0(col_idx,2)-1);
        ymax = mdr0(pldata,i,mdiV0(col_idx,3)-1);
        if (logx)
        {
          x = log10(x);
        }
        if (logy)
        {
          ymin = log10(ymin);
          ymax = log10(ymax);
        }
        plerry(1, (PLFLT *) &x, (PLFLT *) &ymin, (PLFLT *) &ymax);
      }
      else if (err_xy==2)
      {
        x = mdr0(pldata,i,mdiV0(col_idx,0)-1);
        ymin = mdr0(pldata,i,mdiV0(col_idx,4)-1);
        ymax = mdr0(pldata,i,mdiV0(col_idx,5)-1);
        if (logx)
        {
          x = log10(x);
        }
        if (logy)
        {
          ymin = log10(ymin);
          ymax = log10(ymax);
        }
        plerry(1, (PLFLT *) &x, (PLFLT *) &ymin, (PLFLT *) &ymax);
      }
    }
  }

_exit_plprint:

  if (clean_col_idx)
    mdr_Destroy(col_idx);
  ent_Clean (e1);
  ent_Clean (e2);
  ent_Clean (e3);
  ent_Clean (e4);
  ent_Clean (e5);

  return ent_Create_Rlab_Success();
}



/*
* This function handles erros that occur in the plplot library.
* We longjmp from here back to the prompt because the library
* will exit otherwise.
*/

/* in main.c */
extern int dec_buff (void);
extern jmp_buf jmp[];

static int
rlab_plplot_exit (message)
    char *message;
{
  rerror (message);
  /* longjmp( jmp[dec_buff()], 1 ); */

  return (0);			/* shut up compiler */
}

/*
* Make a hardcopy of the current plot-window.
*/

static char default_device[] =
#ifndef __riscos
    "ps"; // Default to BW Postscript
#else
    "arcdraw";
#endif

extern void plsfile (FILE * file);

#undef THIS_SOLVER
#define THIS_SOLVER "_plprint"

Ent *
_plot_plprint (int nargs, Datum args[])
{
  char *fname, *device=default_device;
  FILE *sfile;
  PLINT cur_pls, new_pls;
  Ent *e1=0, *e2=0;

  if (!(nargs == 1 || nargs == 2))
    rerror (THIS_SOLVER ": " RLAB_ERROR_ONE_OR_TWO_ARG_REQUIRED);

  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) != MATRIX_DENSE_STRING)
    rerror (THIS_SOLVER ": " RLAB_ERROR_ARG1_MDS_SCALAR);
  fname = class_char_pointer (e1);

  if (nargs == 2)
  {
    e2 = bltin_get_ent (args[1]);
    if (ent_type (e2) != MATRIX_DENSE_STRING)
      rerror (THIS_SOLVER ": " RLAB_ERROR_ARG2_MDS_SCALAR);
    device = class_char_pointer (e2);
  }

  if (!device)
    device=default_device;
  else if (!(!strcmp (device, "ps") ||
    !strcmp (device, "psc") ||
    !strcmp (device, "xfig") || !strcmp (device, "plmeta") ||
#ifdef __riscos
    !strcmp (device, "arcdraw") ||
#endif
    !strcmp (device, "ljii") ||
    !strcmp (device, "lj_hpgl") || !strcmp (device, "pbm")))
  {
    ent_Clean (e1);
    ent_Clean (e2);
    rerror ("plprint: invalid plot device");
  }

  /*
  * Get the current stream number.
  */
  plgstrm (&cur_pls);

  /*
  * Create a new stream and make it the default.
  * The new stream number is returned (ipls).
  */
  plmkstrm (&new_pls);

  /* Open file for writes */
  sfile = fopen (fname, "w");

  if (sfile)
  {
    /* Initialize stream */
    plsdev (device);
    plsfile (sfile);

    /* Copy the current stream to the newly-made stream */
    plcpstrm (cur_pls, 0);
    pladv (0);

    /* Remake current plot and then switch back to original stream */
    plreplot ();
    plflush ();
    plend1 ();  // this will close 'sfile' too, mk
    plsstrm (cur_pls);

    ent_Clean (e1);
    ent_Clean (e2);
    return ent_Create_Rlab_Success();
  }
  else
  {
    ent_Clean (e1);
    ent_Clean (e2);
    return ent_Create_Rlab_Failure();
  }
}

/*
* Re-plot the current stream
*/

Ent *
_plot_plreplot (int nargs, Datum args[])
{
  plreplot ();

  return ent_Create_Rlab_Success();
}

/*
* Set the number of the current output stream
*/
#undef THIS_SOLVER
#define THIS_SOLVER "_plsstrm "
Ent * _plot_plsstrm (int nargs, Datum args[])
{
  int n;

  if (nargs == 1)
  {
    Ent *e1=bltin_get_ent (args[0]);

    if (ent_type (e1) != MATRIX_DENSE_REAL)
      rerror (THIS_SOLVER ": " RLAB_ERROR_ARG1_MDR_SCALAR);

    n = (int) class_double (e1);
    plsstrm ((PLINT) n);

    ent_Clean (e1);

    return ent_Create_Rlab_Int(n);
  }

  plgstrm (&n);
  return ent_Create_Rlab_Int(n);
}

#undef THIS_SOLVER
#define THIS_SOLVER "_plmkstrm "
Ent * ent_plplot_plmkstrm (int nargs, Datum args[])
{
  PLINT n;
  plmkstrm (&n);
  return ent_Create_Rlab_Double(n);
}

#undef THIS_SOLVER
#define THIS_SOLVER "_plssub"
Ent *
_plot_plssub (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0;
  double n1, n2;

  if (nargs != 2)
    rerror (THIS_SOLVER ": " RLAB_ERROR_TWO_ARG_REQUIRED);

  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) != MATRIX_DENSE_REAL)
    rerror (THIS_SOLVER ": " RLAB_ERROR_ARG1_MDR_SCALAR);
  e2 = bltin_get_ent (args[1]);
  if (ent_type (e2) != MATRIX_DENSE_REAL)
    rerror (THIS_SOLVER ": " RLAB_ERROR_ARG2_MDR_SCALAR);

  n1 = (double) class_double (e1);
  n2 = (double) class_double (e2);

  plssub (n1, n2);

  ent_Clean (e1);
  ent_Clean (e2);

  return ent_Create_Rlab_Success();
}

Ent *
_plot_plinit (int nargs, Datum args[])
{
  plinit ();

  return ent_Create_Rlab_Success();
}

#undef THIS_SOLVER
#define THIS_SOLVER "_plstart"
Ent *
_plot_plstart (int nargs, Datum args[])
{
  char *str;
  PLFLT n1, n2;
  PLStream *strptr;
  Ent *e1=0, *e2=0, *e3=0;

  if (nargs != 3)
    rerror (THIS_SOLVER ": " RLAB_ERROR_THREE_ARG_REQUIRED);

  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) != MATRIX_DENSE_STRING)
    rerror (THIS_SOLVER ": " RLAB_ERROR_ARG1_MDS_SCALAR);
  str = class_char_pointer (e1);

  e2 = bltin_get_ent (args[1]);
  if (ent_type (e2) != MATRIX_DENSE_REAL)
    rerror (THIS_SOLVER ": " RLAB_ERROR_ARG2_MDR_SCALAR);
  n1 = (double) class_double (e2);

  e3 = bltin_get_ent (args[2]);
  if (ent_type (e3) != MATRIX_DENSE_REAL)
    rerror (THIS_SOLVER ": " RLAB_ERROR_ARG3_MDR_SCALAR);
  n2 = (double) class_double (e3);

  if (!init_plot_device)
    plsexit (rlab_plplot_exit);

  plgpls (&strptr);
  strptr->server_nokill = 1;

  plstart (str, n1, n2);

  if (!init_plot_device)
    init_plot_device = 1;

  ent_Clean (e1);
  ent_Clean (e2);
  ent_Clean (e3);

  return ent_Create_Rlab_Success();
}

#undef THIS_SOLVER
#define THIS_SOLVER "_plsfile"
Ent *
_plot_plsfile (int nargs, Datum args[])
{
  char *str;
  Ent *e1=0;
  FILE *fn;

  if (nargs != 1)
    rerror ("_plsfile: 1 argument allowed");

  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) != MATRIX_DENSE_STRING)
    rerror (THIS_SOLVER ": " RLAB_ERROR_ARG1_MDS_SCALAR);
  str = class_char_pointer (e1);

  if (!(fn = get_file_ds (str, "w", 0)))
  {
    fprintf (stderr, "%s, cannot open for read\n", str);

    ent_Clean (e1);

    ent_Create_Rlab_Failure ();
  }

  plsfile (fn);

  return ent_Create_Rlab_Success();
}

#undef THIS_SOLVER
#define THIS_SOLVER "_plenv"
Ent * _plot_plenv (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *e3=0, *e4=0, *e5=0, *e6=0;
  PLFLT xmin, xmax, ymin, ymax;
  PLINT just, axis;

  if (nargs != 6)
  {
    rerror ("_plenv: 6 arguments are required");
  }

  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) != MATRIX_DENSE_REAL)
    rerror (THIS_SOLVER ": " RLAB_ERROR_ARG1_MDR_SCALAR);
  xmin = (PLFLT) class_double (e1);

  e2 = bltin_get_ent (args[1]);
  if (ent_type (e2) != MATRIX_DENSE_REAL)
    rerror (THIS_SOLVER ": " RLAB_ERROR_ARG2_MDR_SCALAR);
  xmax = (PLFLT) class_double (e2);

  e3 = bltin_get_ent (args[2]);
  if (ent_type (e3) != MATRIX_DENSE_REAL)
    rerror (THIS_SOLVER ": " RLAB_ERROR_ARG3_MDR_SCALAR);
  ymin = (PLFLT) class_double (e3);

  e4 = bltin_get_ent (args[3]);
  if (ent_type (e4) != MATRIX_DENSE_REAL)
    rerror (THIS_SOLVER ": " RLAB_ERROR_ARG4_MDR_SCALAR);
  ymax = (PLFLT) class_double (e4);

  e5 = bltin_get_ent (args[4]);
  if (ent_type (e5) != MATRIX_DENSE_REAL)
    rerror (THIS_SOLVER ": " RLAB_ERROR_ARG5_MDR_SCALAR);
  just = (PLINT) class_double (e5);

  e6 = bltin_get_ent (args[5]);
  if (ent_type (e6) != MATRIX_DENSE_REAL)
    rerror (THIS_SOLVER ": " RLAB_ERROR_ARG6_MDR_SCALAR);
  axis = (PLINT) class_double (e6);

  plenv (xmin, xmax, ymin, ymax, just, axis);

  ent_Clean (e1);
  ent_Clean (e2);
  ent_Clean (e3);
  ent_Clean (e4);
  ent_Clean (e5);
  ent_Clean (e6);

  return ent_Create_Rlab_Success();
}

Ent *
_plot_plline (int nargs, Datum args[])
{
  int npoints;
  Ent *e1=0, *e2=0, *e3=0;
  MDR *x, *y;

  if (nargs != 3)
  {
    rerror ("_plline: 3 arguments required");
  }

  e1 = bltin_get_ent (args[0]);
  e2 = bltin_get_ent (args[1]);
  e3 = bltin_get_ent (args[2]);

  npoints = (int) class_double (e1);
  x = class_matrix_real (e2);
  y = class_matrix_real (e3);

  plline (npoints, (PLFLT *) MDRPTR (x), (PLFLT *) MDRPTR (y));

  ent_Clean (e1);
  ent_Clean (e2);
  ent_Clean (e3);

  return ent_Create_Rlab_Success();
}

Ent *
_plot_plline3 (int nargs, Datum args[])
{
  int npoints;
  Ent *e1=0, *e2=0, *e3=0, *e4=0;
  MDR *x=0, *y=0, *z=0;

  if (nargs != 4)
  {
    rerror ("_plline3: 4 arguments required");
  }

  e1 = bltin_get_ent (args[0]);
  e2 = bltin_get_ent (args[1]);
  e3 = bltin_get_ent (args[2]);
  e4 = bltin_get_ent (args[3]);

  npoints = (int) class_double (e1);

  x = class_matrix_real (e2);
  y = class_matrix_real (e3);
  z = class_matrix_real (e4);

  plline3 (npoints, (PLFLT *) MDRPTR (x),
    (PLFLT *) MDRPTR (y), (PLFLT *) MDRPTR (z));

  ent_Clean (e1);
  ent_Clean (e2);
  ent_Clean (e3);
  ent_Clean (e4);

  return ent_Create_Rlab_Success();
}

Ent *
_plot_plpoly3 (int nargs, Datum args[])
{
  int i, npoints;
  Ent *e1=0, *e2=0, *e3, *e4, *e5;
  MDR *x, *y, *z, *k;
  int *ik, size, ifcc=1;

  if (nargs != 5)
  {
    rerror ("_plpoly3: 5 arguments required");
  }

  e1 = bltin_get_ent (args[0]);
  e2 = bltin_get_ent (args[1]);
  e3 = bltin_get_ent (args[2]);
  e4 = bltin_get_ent (args[3]);
  e5 = bltin_get_ent (args[4]);

  npoints = (int) class_double (e1);
  x = class_matrix_real (e2);
  y = class_matrix_real (e3);
  z = class_matrix_real (e4);
  k = class_matrix_real (e5);

  /* Copy the draw/no-draw vector, cause it needs to be int. */
  ik = (PLINT *) GC_malloc_atomic_ignore_off_page
    (sizeof (PLINT) * MNR (k) * MNC (k));
  size = MNR (k) * MNC (k);

  for (i = 0; i < size; i++)
  {
    ik[i] = MdrV0 (k, i);
  }

  c_plpoly3 (npoints, (const PLFLT *) MDRPTR (x), (const PLFLT *) MDRPTR (y), (const PLFLT *) MDRPTR (z),
            (PLINT *) ik, ifcc);

  GC_FREE (ik);

  ent_Clean (e1);
  ent_Clean (e2);
  ent_Clean (e3);
  ent_Clean (e4);
  ent_Clean (e5);

  return ent_Create_Rlab_Success();
}

Ent *
_plot_plend (int nargs, Datum args[])
{
  if (nargs != 0)
  {
    rerror ("_plend: no arguments allowed");
  }

  plend ();

  return ent_Create_Rlab_Success();
}

Ent *
_plot_plend1 (int nargs, Datum args[])
{
  if (nargs != 0)
  {
    rerror ("_plend1: no arguments allowed");
  }

  plend1 ();

  return ent_Create_Rlab_Success();
}

Ent *
_plot_pllab (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *e3;
  char *sx, *sy, *st;

  if (nargs != 3)
  {
    rerror ("_pllab: 3 arguments required");
  }

  e1 = bltin_get_ent (args[0]);
  e2 = bltin_get_ent (args[1]);
  e3 = bltin_get_ent (args[2]);

  sx = class_char_pointer (e1);
  sy = class_char_pointer (e2);
  st = class_char_pointer (e3);

  pllab (sx, sy, st);

  ent_Clean (e1);
  ent_Clean (e2);
  ent_Clean (e3);

  return ent_Create_Rlab_Success();
}

/*
* Set the color for subsequent lines. You may specify which
* colors are associated with which number by the plancol
* routine if supported by the output driver.
*/

Ent *
_plot_plcol (int nargs, Datum args[])
{
  Ent *e1=0;
  PLINT color;

  if (nargs != 1)
  {
    rerror ("_plcol: 1 argument required");
  }

  e1 = bltin_get_ent (args[0]);
  color = (PLINT) class_double (e1);

  if (color >=0 )
    plcol0 (color);

  return ent_Create_Rlab_Success();
}

/*
* Set the background color by 8 bit RGB value
* Note: for some drivers this corresponds to a cmap 0 color
*/
#undef  THIS_SOLVER
#define THIS_SOLVER "_plscolbg"
Ent * _plot_plscolbg (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *e3=0;
  PLINT r, g, b;
  FILE *rlab_stderr = (!RLAB_STDERR_DS) ? stderr : RLAB_STDERR_DS;

  if (nargs != 3 && nargs!=1)
  {
    fprintf (rlab_stderr, THIS_SOLVER ": " RLAB_ERROR_ONE_OR_THREE_ARG_REQUIRED "\n");
    goto _exit_plscolbg;
  }

  e1 = bltin_get_ent (args[0]);
  if (ent_type(e1)!=MATRIX_DENSE_REAL)
    goto _exit_plscolbg;

  if (nargs == 3)
  {
    e2 = bltin_get_ent (args[1]);
    if (ent_type(e2)!=MATRIX_DENSE_REAL)
      goto _exit_plscolbg;

    e3 = bltin_get_ent (args[2]);
    if (ent_type(e3)!=MATRIX_DENSE_REAL)
      goto _exit_plscolbg;

    r = (PLINT) class_double (e1);
    g = (PLINT) class_double (e2);
    b = (PLINT) class_double (e3);
  }
  else
  {
    MDR *x = ent_data(e1);
    if (SIZE(x)!=3)
      goto _exit_plscolbg;

    r = (PLINT) mdiV0(x,0);
    g = (PLINT) mdiV0(x,1);
    b = (PLINT) mdiV0(x,2);
  }

  plscolbg (r, g, b);

_exit_plscolbg:

  ent_Clean (e1);
  ent_Clean (e2);
  ent_Clean (e3);
  return ent_Create_Rlab_Success();
}

/*
* Set a given color from color map 0 by 8 bit RGB value
* Does not result in any additional cells to be allocated.
*/
#undef  THIS_SOLVER
#define THIS_SOLVER "_plscol0"
Ent * _plot_plscol0 (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *e3=0, *e4=0;
  int i, nr;
  MDR *color=0, *r=0, *g=0, *b=0;
  FILE *rlab_stderr = (!RLAB_STDERR_DS) ? stderr : RLAB_STDERR_DS;

  if (nargs != 4 && nargs!=2)
  {
    fprintf (rlab_stderr, THIS_SOLVER ": " RLAB_ERROR_TWO_OR_FOUR_ARG_REQUIRED "\n");
    goto _exit_plscol0;
  }

  e1 = bltin_get_ent (args[0]);
  if (ent_type(e1)!=MATRIX_DENSE_REAL)
    goto _exit_plscol0;
  color = ent_data(e1);
  nr = SIZE(color);
  if (nr<1)
    goto _exit_plscol0;

  if (nargs == 4)
  {
    e2 = bltin_get_ent (args[1]);
    if (ent_type(e2)!=MATRIX_DENSE_REAL)
      goto _exit_plscol0;

    e3 = bltin_get_ent (args[2]);
    if (ent_type(e3)!=MATRIX_DENSE_REAL)
      goto _exit_plscol0;

    e4 = bltin_get_ent (args[3]);
    if (ent_type(e4)!=MATRIX_DENSE_REAL)
      goto _exit_plscol0;

    r = ent_data (e2);
    g = ent_data (e3);
    b = ent_data (e4);

    if ((SIZE(r)!=SIZE(g)) && (SIZE(b)!=SIZE(g)))
      goto _exit_plscol0;

    for (i=0; i<nr; i++)
    {
      plscol0 ( mdiV0(color,i), mdiV0(r,MIN(i,SIZE(r)-1)),  mdiV0(g,MIN(i,SIZE(g)-1)),  mdiV0(b,MIN(i,SIZE(b)-1)));
    }
  }
  else
  {
    e2 = bltin_get_ent (args[1]);
    if (ent_type(e2)!=MATRIX_DENSE_REAL)
      goto _exit_plscol0;

    MDR *x = ent_data(e2);
    if (MNC(x)!=3)
      goto _exit_plscol0;

    for (i=0; i<nr; i++)
    {
      plscol0 ( mdiV0(color,i), mdi0(x,MIN(i,MNR(x)-1),0),  mdi0(x,MIN(i,MNR(x)-1),1),  mdi0(x,MIN(i,MNR(x)-1),2));
    }
  }

_exit_plscol0:

  ent_Clean (e1);
  ent_Clean (e2);
  ent_Clean (e3);
  ent_Clean (e4);

  return ent_Create_Rlab_Success();
}


#undef  THIS_SOLVER
#define THIS_SOLVER "_plgcol0"
Ent * ent_plplot_plgcol0 (int nargs, Datum args[])
{
  Ent *e1=0;
  int i, nr;
  MDR *color=0, *rgb=0;
  FILE *rlab_stderr = (!RLAB_STDERR_DS) ? stderr : RLAB_STDERR_DS;

  if (nargs != 0 && nargs!=1)
  {
    fprintf (rlab_stderr, THIS_SOLVER ": " RLAB_ERROR_NONE_OR_ONE_ARG_REQUIRED "\n");
    goto _exit_plgcol0;
  }

  if (nargs == 1)
  {
    e1 = bltin_get_ent (args[0]);
    if (ent_type(e1)==MATRIX_DENSE_REAL)
      color = ent_data(e1);
    nr = SIZE(color); 
    if (nr>0)
    {
      rgb = mdi_Create(nr,3);
      for (i=0; i<nr; i++)
      {
        plgcol0(mdiV0(color,i),&Mdi0(rgb,i,0),&Mdi0(rgb,i,1),&Mdi0(rgb,i,2));
      }
      goto _exit_plgcol0;
    }
  }

  nr = 16;
  rgb = mdi_Create(nr,3);
  for (i=0; i<nr; i++)
  {
    plgcol0(i,&Mdi0(rgb,i,0),&Mdi0(rgb,i,1),&Mdi0(rgb,i,2));
  }

_exit_plgcol0:

  ent_Clean (e1);
  return ent_Assign_Rlab_MDR(rgb);
}



/*
* Sets the line style according to one of eight predefined patters.
*/

Ent *
_plot_pllsty (int nargs, Datum args[])
{
  Ent *e1=0;
  PLINT style;

  if (nargs != 1)
  {
    rerror ("_pllsty: 1 argument required");
  }

  e1 = bltin_get_ent (args[0]);

  style = (PLINT) class_double (e1);
  pllsty (style);

  ent_Clean (e1);
  return ent_Create_Rlab_Success();
}

Ent *
_plot_plclr (int nargs, Datum args[])
{
  if (nargs != 0)
  {
    rerror ("_plclr: no argument allowed");
  }

  pleop ();

  return ent_Create_Rlab_Success();
}

Ent *
_plot_plpoin (int nargs, Datum args[])
{
  PLINT n, code;
  Ent *e1=0, *e2=0, *e3, *e4;
  MDR *X, *Y;

  if (nargs != 4)
  {
    rerror ("_plpoin: 4 arguments required");
  }

  e1 = bltin_get_ent (args[0]);
  e2 = bltin_get_ent (args[1]);
  e3 = bltin_get_ent (args[2]);
  e4 = bltin_get_ent (args[3]);

  n = (PLINT) class_double (e1);
  X = class_matrix_real (e2);
  Y = class_matrix_real (e3);
  code = (PLINT) class_double (e4);

  plpoin (n, MDRPTR (X), MDRPTR (Y), code);

  ent_Clean (e1);
  ent_Clean (e2);
  ent_Clean (e3);
  ent_Clean (e4);

  return ent_Create_Rlab_Success();
}

Ent *
_plot_plpoin3 (int nargs, Datum args[])
{
  PLINT n, code;
  Ent *e1=0, *e2=0, *e3, *e4, *e5;
  MDR *X, *Y, *Z;

  if (nargs != 5)
  {
    rerror ("_plpoin3: 5 arguments required");
  }

  e1 = bltin_get_ent (args[0]);
  e2 = bltin_get_ent (args[1]);
  e3 = bltin_get_ent (args[2]);
  e4 = bltin_get_ent (args[3]);
  e5 = bltin_get_ent (args[4]);

  n = (PLINT) class_double (e1);
  X = class_matrix_real (e2);
  Y = class_matrix_real (e3);
  Z = class_matrix_real (e4);
  code = (PLINT) class_double (e5);

  plpoin3 (n, MDRPTR (X), MDRPTR (Y), MDRPTR (Z), code);

  ent_Clean (e1);
  ent_Clean (e2);
  ent_Clean (e3);
  ent_Clean (e4);
  ent_Clean (e5);

  return ent_Create_Rlab_Success();
}

Ent *
_plot_plhist (int nargs, Datum args[])
{
  PLINT n, nbin, oldwin;
  PLFLT datmin, datmax;
  PLFLT *data;
  Ent *e1=0, *e2=0, *e3, *e4, *e5, *e6;
  MDR *X;

  if (nargs != 6)
  {
    rerror ("_plhist: 6 arguments required");
  }

  e1 = bltin_get_ent (args[0]);
  e2 = bltin_get_ent (args[1]);
  e3 = bltin_get_ent (args[2]);
  e4 = bltin_get_ent (args[3]);
  e5 = bltin_get_ent (args[4]);
  e6 = bltin_get_ent (args[5]);

  n = (PLINT) class_double (e1);
  X = class_matrix_real (e2);
  data = (PLFLT *) MDRPTR (X);
  datmin = (PLFLT) class_double (e3);
  datmax = (PLFLT) class_double (e4);
  nbin = (PLINT) class_double (e5);
  oldwin = (PLINT) class_double (e6);

  plhist (n, data, datmin, datmax, nbin, oldwin);

  ent_Clean (e1);
  ent_Clean (e2);
  ent_Clean (e3);
  ent_Clean (e4);
  ent_Clean (e5);
  ent_Clean (e6);

  return ent_Create_Rlab_Success();
}

Ent *
_plot_plspage (int nargs, Datum args[])
{
  PLINT xp, yp, xleng, yleng, xoff, yoff;
  Ent *e1=0, *e2=0, *e3, *e4, *e5, *e6;

  if (nargs != 6)
  {
    rerror ("_plspage: 6 arguments required");
  }

  e1 = bltin_get_ent (args[0]);
  e2 = bltin_get_ent (args[1]);
  e3 = bltin_get_ent (args[2]);
  e4 = bltin_get_ent (args[3]);
  e5 = bltin_get_ent (args[4]);
  e6 = bltin_get_ent (args[5]);

  xp = (PLINT) class_double (e1);
  yp = (PLINT) class_double (e2);
  xleng = (PLINT) class_double (e3);
  yleng = (PLINT) class_double (e4);
  xoff = (PLINT) class_double (e5);
  yoff = (PLINT) class_double (e6);

  plspage (xp, yp, xleng, yleng, xoff, yoff);

  ent_Clean (e1);
  ent_Clean (e2);
  ent_Clean (e3);
  ent_Clean (e4);
  ent_Clean (e5);
  ent_Clean (e6);

  return ent_Create_Rlab_Success();
}

/*
* Advance to the next subpage if sub=0, performing a page
* advance if there are no remaining subpages on the current page.
* If subwindowing is not being used, pladv(0) will always advance
* the page. If sub>0, PLPLOT switches to the specified subpage.
* Note that this alows you to overwrite a plot on the specified
* subpage; if this is not what you intended, use plclr followed by
* plpage to first advance the page. This routine is called automatically
* (with sub=0) by plenv, but if plenv is not used, pladv must be called
* after initializing PLPLOT but before defining the viewport.
*/

Ent *
_plot_pladv (int nargs, Datum args[])
{
  Ent *e1;
  PLINT sub;

  if (nargs != 1)
  {
    rerror ("_pladv: 1 argument allowed");
  }

  e1 = bltin_get_ent (args[0]);
  sub = (PLINT) class_double (e1);
  pladv (sub);

  ent_Clean (e1);

  return ent_Create_Rlab_Success();
}

Ent *
_plot_plgra (int nargs, Datum args[])
{
  if (nargs != 0)
  {
    rerror ("_plgra: no arguments allowed");
  }

  plgra ();

  return ent_Create_Rlab_Success();
}

Ent * _plot_pltext (int nargs, Datum args[])
{
  if (nargs != 0)
  {
    rerror ("_pltext: no arguments allowed");
  }

  pltext ();

  return ent_Create_Rlab_Success();
}

Ent *
_plot_plflush (int nargs, Datum args[])
{
  if (nargs != 0)
  {
    rerror ("_plflush: no arguments allowed");
  }

  plflush ();

  return ent_Create_Rlab_Success();
}

Ent *
_plot_plbox (int nargs, Datum args[])
{
  char *xopt, *yopt;
  PLFLT xtick, ytick;
  PLINT nxsub, nysub;
  Ent *e1=0, *e2=0, *e3, *e4, *e5, *e6;

  if (nargs != 6)
  {
    rerror ("_plbox: 6 arguments required");
  }

  e1 = bltin_get_ent (args[0]);
  e2 = bltin_get_ent (args[1]);
  e3 = bltin_get_ent (args[2]);
  e4 = bltin_get_ent (args[3]);
  e5 = bltin_get_ent (args[4]);
  e6 = bltin_get_ent (args[5]);

  xopt = class_char_pointer (e1);
  yopt = class_char_pointer (e4);

  xtick = (PLFLT) class_double (e2);
  ytick = (PLFLT) class_double (e5);

  nxsub = (PLINT) class_double (e3);
  nysub = (PLINT) class_double (e6);

  plbox (xopt, xtick, nxsub, yopt, ytick, nysub);

  ent_Clean (e1);
  ent_Clean (e2);
  ent_Clean (e3);
  ent_Clean (e4);
  ent_Clean (e5);
  ent_Clean (e6);

  return ent_Create_Rlab_Success();
}

/*
* Draw axes, numeric and text labels for a
* three-dimensional surface plot.
*/

Ent *
_plot_plbox3 (int nargs, Datum args[])
{
  char *xopt, *yopt, *zopt;
  char *xlabel, *ylabel, *zlabel;
  PLFLT xtick, ytick, ztick;
  PLINT nxsub, nysub, nzsub;
  Ent *e1=0, *e2=0, *e3, *e4, *e5, *e6;
  Ent *e7, *e8, *e9, *e10, *e11, *e12;

  if (nargs != 12)
  {
    rerror ("_plbox3: 12 arguments required");
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

  xopt = class_char_pointer (e1);
  yopt = class_char_pointer (e5);
  zopt = class_char_pointer (e9);

  xlabel = class_char_pointer (e2);
  ylabel = class_char_pointer (e6);
  zlabel = class_char_pointer (e10);

  xtick = (PLFLT) class_double (e3);
  ytick = (PLFLT) class_double (e7);
  ztick = (PLFLT) class_double (e11);

  nxsub = (PLINT) class_double (e4);
  nysub = (PLINT) class_double (e8);
  nzsub = (PLINT) class_double (e12);

  plbox3 (xopt, xlabel, xtick, nxsub,
    yopt, ylabel, ytick, nysub, zopt, zlabel, ztick, nzsub);

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

Ent *
_plot_plvsta (int nargs, Datum args[])
{
  if (nargs != 0)
  {
    rerror ("_plvsta: no arguments allowed");
  }

  plvsta ();

  return ent_Create_Rlab_Success();
}

Ent *
_plot_plvasp (int nargs, Datum args[])
{
  Ent *e1;
  PLFLT aspect;

  if (nargs != 1)
  {
    rerror ("_plvasp: 1 argument allowed");
  }

  e1 = bltin_get_ent (args[0]);
  aspect = (PLFLT) class_double (e1);

  plvasp (aspect);

  ent_Clean (e1);

  return ent_Create_Rlab_Success();
}

Ent *
_plot_plwind (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *e3, *e4;
  PLFLT xmin, xmax, ymin, ymax;

  if (nargs != 4)
  {
    rerror ("_plwind: 4 arguments required");
  }

  e1 = bltin_get_ent (args[0]);
  e2 = bltin_get_ent (args[1]);
  e3 = bltin_get_ent (args[2]);
  e4 = bltin_get_ent (args[3]);

  xmin = (PLFLT) class_double (e1);
  xmax = (PLFLT) class_double (e2);
  ymin = (PLFLT) class_double (e3);
  ymax = (PLFLT) class_double (e4);

  plwind (xmin, xmax, ymin, ymax);

  ent_Clean (e1);
  ent_Clean (e2);
  ent_Clean (e3);
  ent_Clean (e4);

  return ent_Create_Rlab_Success();
}

/*
* Create a three-dimensional surface plot...
*/

Ent *
_plot_plot3d (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *e3=0, *e4=0, *e5=0, *e6=0, *e7=0;
  PLFLT *x, *y;
  PLINT nx, ny, opt, side;
  MDR *z;

  if (nargs != 7)
  {
    rerror ("_plot3d: 7 arguments required");
  }

  e1 = bltin_get_ent (args[0]);
  e2 = bltin_get_ent (args[1]);
  e3 = bltin_get_ent (args[2]);
  e4 = bltin_get_ent (args[3]);
  e5 = bltin_get_ent (args[4]);
  e6 = bltin_get_ent (args[5]);
  e7 = bltin_get_ent (args[6]);

  x   = (PLFLT *) MDRPTR (class_matrix_real (e1));
  if (!x)
    rerror (THIS_SOLVER ": " RLAB_ERROR_ARG1_MDR_VECTOR);
  y   = (PLFLT *) MDRPTR (class_matrix_real (e2));
  if (!y)
    rerror (THIS_SOLVER ": " RLAB_ERROR_ARG2_MDR_VECTOR);
  z   = class_matrix_real (e3);
  if (!z)
    rerror (THIS_SOLVER ": " RLAB_ERROR_ARG3_MDR_MATRIX);

  nx  = (PLINT) class_double (e4);
  ny  = (PLINT) class_double (e5);
  opt = (PLINT) class_double (e6);

  PLfGrid zp;
  zp.nx = nx;
  zp.ny = ny;
  zp.f = (PLFLT *) MDPTR(z);

  side = (PLINT) class_double (e7);

  plfplot3d(x, y, (PLF2OPS) plf2ops_grid_col_major(), &zp, nx, ny, opt, side );

  zp.f = 0;

  ent_Clean (e1);
  ent_Clean (e2);
  ent_Clean (e3);
  ent_Clean (e4);
  ent_Clean (e5);
  ent_Clean (e6);
  ent_Clean (e7);

  return ent_Create_Rlab_Success();
}

#undef THIS_SOLVER
#define THIS_SOLVER "_plmesh"
Ent *
_plot_plmesh (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *e3=0, *e4=0, *e5=0, *e6=0;
  MDR *z=0;

  PLFLT *x, *y;
  PLINT nx, ny, opt;


  if (nargs != 6)
    rerror (THIS_SOLVER ": 6 arguments required");

  e1 = bltin_get_ent (args[0]);
  e2 = bltin_get_ent (args[1]);
  e3 = bltin_get_ent (args[2]);
  e4 = bltin_get_ent (args[3]);
  e5 = bltin_get_ent (args[4]);
  e6 = bltin_get_ent (args[5]);

  x   = (PLFLT *) MDRPTR (class_matrix_real (e1));
  if (!x)
    rerror (THIS_SOLVER ": " RLAB_ERROR_ARG1_MDR_VECTOR);
  y   = (PLFLT *) MDRPTR (class_matrix_real (e2));
  if (!y)
    rerror (THIS_SOLVER ": " RLAB_ERROR_ARG2_MDR_VECTOR);
  z   = class_matrix_real (e3);
  if (!z)
    rerror (THIS_SOLVER ": " RLAB_ERROR_ARG3_MDR_MATRIX);

  nx  = (PLINT) class_double (e4);
  ny  = (PLINT) class_double (e5);
  opt = (PLINT) class_double (e6);

  PLfGrid zp;
  zp.nx = nx;
  zp.ny = ny;
  zp.f = (PLFLT *) MDPTR(z);

  plfmesh (x, y, (PLF2OPS) plf2ops_grid_col_major(), &zp, nx, ny, opt);
  zp.f = 0;

  ent_Clean (e1);
  ent_Clean (e2);
  ent_Clean (e3);
  ent_Clean (e4);
  ent_Clean (e5);
  ent_Clean (e6);

  return ent_Create_Rlab_Success();
}

Ent *
_plot_plw3d (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *e3, *e4, *e5, *e6, *e7;
  Ent *e8, *e9, *e10, *e11;
  PLFLT basex, basey, height, xmin, xmax;
  PLFLT ymin, ymax, zmin, zmax, alt, az;

  if (nargs != 11)
  {
    rerror ("_plw3d: 11 arguments required");
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

  basex = (PLFLT) class_double (e1);
  basey = (PLFLT) class_double (e2);
  height = (PLFLT) class_double (e3);
  xmin = (PLFLT) class_double (e4);
  xmax = (PLFLT) class_double (e5);
  ymin = (PLFLT) class_double (e6);
  ymax = (PLFLT) class_double (e7);
  zmin = (PLFLT) class_double (e8);
  zmax = (PLFLT) class_double (e9);
  alt = (PLFLT) class_double (e10);
  az = (PLFLT) class_double (e11);

  plw3d (basex, basey, height, xmin, xmax, ymin, ymax, zmin, zmax, alt, az);

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

//
// plprint3d (data, format)
//  data:
//    (z,format)          z is full matrix, vector for description
//    ([x,y,z],format)    3-column matrix for line, vector for description
//    (x,y,z,format)      3 x column vector, vector for description
// additional options:
//    1 - DRAW_LINEX   draw lines parallel to the X axis
//    2 - DRAW_LINEY   draw lines parallel to the Y axis
//    4 - DRAW_LINEXY  draw lines parallel to both the X and Y axis
//    8 - MAG_COLOR    draw the mesh with a color dependent of the magnitude
//    16 - BASE_CONT    draw contour plot at bottom xy plane
//    32 - TOP_CONT     draw contour plot at top xy plane (not yet)
//    64 - DRAW_SIDES   draw sides
#undef THIS_SOLVER
#define THIS_SOLVER "ent_plplot3_plprintf"
Ent * ent_plplot3_printf (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *e3=0, *e4=0, *e5=0, *e6=0,*e7=0;
  int i, j, ndata_x, ndata_y, ndata_z=0, logx=0, logy=0, logz=0, idx_ent=0;
  FILE *rlab_stderr = (!RLAB_STDERR_DS) ? stderr : RLAB_STDERR_DS;
  double *x=0, *x1=0, *y=0, *y1=0, *z=0, *z1=0;
  MDR *plfmt=0, *plxdata=0, *plydata=0, *plzdata=0;

  int clear_x=0, clear_y=0, clear_z=0, opt=DRAW_LINEXY;

  if ((nargs != 5) && (nargs != 7))
  {
    fprintf(rlab_stderr, THIS_SOLVER ": " RLAB_ERROR_TWO_OR_FOUR_ARG_REQUIRED "\n");
    goto _exit_plprint3d;
  }

  //
  // (x,y,z,format), where x,y,z have to be double matrices! no ints!
  //
  e1 = bltin_get_ent(args[idx_ent++]);
  if (ent_type(e1) != MATRIX_DENSE_REAL)
  {
    fprintf(rlab_stderr, THIS_SOLVER ": " RLAB_ERROR_ARG1_MDR_VECTOR"\n");
    goto _exit_plprint3d;
  }
  plxdata = ent_data(e1);
  if (!MD_TYPE_DOUBLE(plxdata))
    goto _exit_plprint3d;
  ndata_x = SIZE(plxdata);
  if (ndata_x < 1)
    goto _exit_plprint3d;

  if (nargs == 7)
  {
    // x:
    x = MDPTR(plxdata);

    // y:
    e2 = bltin_get_ent(args[idx_ent++]);
    if (ent_type(e2) != MATRIX_DENSE_REAL)
    {
      fprintf(rlab_stderr, THIS_SOLVER ": " RLAB_ERROR_ARG2_MDR_VECTOR "\n");
      goto _exit_plprint3d;
    }
    plydata = ent_data(e2);
    ndata_y = SIZE(plydata);
    if (!MD_TYPE_DOUBLE(plydata))
      goto _exit_plprint3d;
    if (ndata_y < 1)
      goto _exit_plprint3d;
    y = MDPTR(plydata);

    // z:
    e3 = bltin_get_ent(args[idx_ent++]);
    if (ent_type(e3) != MATRIX_DENSE_REAL)
    {
      fprintf(rlab_stderr, THIS_SOLVER ": " RLAB_ERROR_ARG3_MDR_VECTOR "\n");
      goto _exit_plprint3d;
    }
    plzdata = ent_data(e3);
    if (ndata_x != MNR(plzdata))
      goto _exit_plprint3d;
    if (ndata_y != MNC(plzdata))
      goto _exit_plprint3d;
    if (!MD_TYPE_DOUBLE(plzdata))
      goto _exit_plprint3d;
    z = MDPTR(plzdata);
    ndata_z = ndata_x * ndata_y;
  }
  else
  {
    //
    // (z, format), or
    // ([xyz], format)
    //
    if (MNC(plxdata)==3)
    {
      x = &Mdr0(plxdata,0,0);
      y = &Mdr0(plxdata,0,1);
      z = &Mdr0(plxdata,0,2);
      ndata_y = ndata_z = ndata_x = MNR(plxdata);
    }
    else
    {
      ndata_x = MNR(plxdata);
      ndata_y = MNC(plxdata);
      ndata_z = ndata_x * ndata_y;
      clear_x = 1;
      x = GC_MALLOC(ndata_x * sizeof(double));
      for (i=0;i<ndata_x; i++)
        x[i] = i+1;
      clear_y = 1;
      y = GC_MALLOC(ndata_y * sizeof(double));
      for (i=0;i<ndata_y; i++)
        y[i] = i+1;
      z = MDPTR(plxdata);
    }
  }

  //
  // next argument: how to plot it
  //
  e4 = bltin_get_ent(args[idx_ent++]);
  if (ent_type(e4) != MATRIX_DENSE_REAL)
  {
    fprintf(rlab_stderr, THIS_SOLVER ": " RLAB_ERROR_ARG4_MDR_VECTOR "\n");
    goto _exit_plprint3d;
  }
  plfmt = ent_data(e4);
  if (SIZE(plfmt)!=10)
  {
    fprintf(rlab_stderr, THIS_SOLVER ": " RLAB_ERROR_ARG4_MDR_VECTOR "\n");
    fprintf(rlab_stderr, THIS_SOLVER ": 10-component vector expected but %i-component was provided.\n"
        THIS_SOLVER ": Cannot continue!\n", SIZE(plfmt));
    goto _exit_plprint3d;
  }

  //
  // data set to be plotted
  //
  int lpb = mdiV0(plfmt, RLABPLUS_PLPLOT_IDX_LINEPOINT);
  int lt  = mdiV0(plfmt, RLABPLUS_PLPLOT_IDX_LINETYPE);
  int lc  = mdiV0(plfmt, RLABPLUS_PLPLOT_IDX_LINECOLOR);
  double lw  = mdrV0(plfmt, RLABPLUS_PLPLOT_IDX_LINEWIDTH);
  int pt  = mdiV0(plfmt, RLABPLUS_PLPLOT_IDX_POINTTYPE);
  int pc  = mdiV0(plfmt, RLABPLUS_PLPLOT_IDX_POINTCOLOR);
  double ps  = mdrV0(plfmt, RLABPLUS_PLPLOT_IDX_POINTSIZE);
  int dr = mdrV0(plfmt, 8);

  //    1 - DRAW_LINEX   draw lines parallel to the X axis
//    2 - DRAW_LINEY   draw lines parallel to the Y axis
//    4 - DRAW_LINEXY  draw lines parallel to both the X and Y axis
//    8 - MAG_COLOR    draw the mesh with a color dependent of the magnitude
//    16 - BASE_CONT    draw contour plot at bottom xy plane
//    32 - TOP_CONT     draw contour plot at top xy plane (not yet)
//    64 - DRAW_SIDES   draw sides

  if (dr < 0)
  {
    opt = 0;
    dr = -dr;
    if (dr & 1)
      opt |= DRAW_LINEX;
    if (dr & 2)
      opt |= DRAW_LINEY;
    if (dr & 4)
      opt |= DRAW_LINEXY;
    if (dr & 8)
      opt |= MAG_COLOR;
    if (dr & 16)
      opt |= BASE_CONT;
    if (dr & 64)
      opt |= DRAW_SIDES;
  }

  //
  // any log axes?
  //
  if (nargs >= idx_ent)
  {
    e5 = bltin_get_ent(args[idx_ent++]);
    if (ent_type(e5) == MATRIX_DENSE_REAL)
    {
      logx = (class_double(e5)>0);
    }
  }
  if (nargs >= idx_ent)
  {
    e6 = bltin_get_ent(args[idx_ent++]);
    if (ent_type(e6) == MATRIX_DENSE_REAL)
    {
      logy = (class_double(e6)>0);
    }
  }
  if (nargs >= idx_ent)
  {
    e7 = bltin_get_ent(args[idx_ent]);
    if (ent_type(e7) == MATRIX_DENSE_REAL)
    {
      logz = (class_double(e7)>0);
    }
  }

  //
  // process logarithms next
  //
  if (logx)
  {
    if (clear_x)
    {
      for (i=0; i<ndata_x; i++)
        x[i] = log10(x[i]);
    }
    else
    {
      clear_x = 1;
      x1 = GC_MALLOC(ndata_x * sizeof(double));
      for (i=0;i<ndata_x; i++)
        x1[i] = log10(x[i]);
      x  = x1;
      x1 = 0;
    }
  }
  if (logy)
  {
    if (clear_y)
    {
      for (i=0; i<ndata_y; i++)
        y[i] = log10(y[i]);
    }
    else
    {
      clear_y = 1;
      y1 = GC_MALLOC(ndata_y * sizeof(double));
      for (i=0;i<ndata_y; i++)
        y1[i] = log10(y[i]);
      y  = y1;
      y1 = 0;
    }
  }
  if (logz)
  {
    if (clear_z)
    {
      for (i=0; i<ndata_z; i++)
        z[i] = log10(z[i]);
    }
    else
    {
      clear_z = 1;
      z1 = GC_MALLOC(ndata_z * sizeof(double));
      for (i=0;i<ndata_z; i++)
        z1[i] = log10(z[i]);
      z  = z1;
      z1 = 0;
    }
  }

  //
  // plot line first
  //
  if ((lpb==0) || (lpb>=2))
  {
    // with lines, with both
    pllsty(lt);   // set line style
    plcol0(lc);   // set line color
#ifdef plwidth
    plwidth (lw);// set line width
#else
    plwid (lw);
#endif
    if (ndata_x * ndata_y == ndata_z)
    {
      PLfGrid zp;
      zp.nx = ndata_x;
      zp.ny = ndata_y;
      zp.f  = z;
      plfmesh (x, y, (PLF2OPS) plf2ops_grid_col_major(), &zp, ndata_x, ndata_y, opt);
      zp.f = 0;
    }
    else if ((ndata_x == ndata_z) && (ndata_y == ndata_z))
    {
      plline3(ndata_x, x, y, z);
    }
  }

  //
  // plot points next
  //
  if ((lpb==1) || (lpb==2))
  {
    // with points, with both
    plssym(0,ps); // set symbol size
    plcol0(pc);   // set symbol color

    if (ndata_x * ndata_y == ndata_z)
    {
      for (i=0; i<ndata_x; i++)
        for (j=0; j<ndata_y; j++)
          plpoin3(1, &x[i], &y[j], &z[ndata_y*i + j], pt);
    }
    else if ((ndata_x == ndata_z) && (ndata_y == ndata_z))
    {
      plpoin3(ndata_x, x, y, z, pt);
    }
  }


_exit_plprint3d:

  if (clear_x)
    GC_FREE(x);
  if (clear_y)
    GC_FREE(y);
  if (clear_z)
    GC_FREE(z);

  ent_Clean (e1);
  ent_Clean (e2);
  ent_Clean (e3);
  ent_Clean (e4);
  ent_Clean (e5);
  ent_Clean (e6);
  ent_Clean (e7);

  return ent_Create_Rlab_Success();
}





/*
* Writes text at a specified position relative to the
* viewport boundaries. Text may be written inside or
* outside the viewport, but is clipped at the subpage
* boundaries.
*/

Ent *
_plot_plmtex (int nargs, Datum args[])
{
  char *side=0, *text=0;
  Ent *e1=0, *e2=0, *e3=0, *e4=0, *e5=0;
  PLFLT disp, pos, just;

  if (nargs != 5)
  {
    rerror ("_plmtex: 5 arguments required");
  }

  e1 = bltin_get_ent (args[0]);
  e2 = bltin_get_ent (args[1]);
  e3 = bltin_get_ent (args[2]);
  e4 = bltin_get_ent (args[3]);
  e5 = bltin_get_ent (args[4]);

  side = class_char_pointer (e1);
  text = class_char_pointer (e5);

  disp  = (PLFLT) class_double (e2);
  pos   = (PLFLT) class_double (e3);
  just  = (PLFLT) class_double (e4);

  plmtex (side, disp, pos, just, text);

  ent_Clean (e1);
  ent_Clean (e2);
  ent_Clean (e3);
  ent_Clean (e4);
  ent_Clean (e5);
  return ent_Create_Rlab_Success();
}

/*
* Turn PLPLOT pause off (arg=1)
* Turn PLPLOT pause on  (arg=0)
*/

Ent *
_plot_plspause (int nargs, Datum args[])
{
  Ent *e1;
  PLINT pause;

  if (nargs != 1)
  {
    rerror ("_plspause: 1 argument required");
  }

  e1 = bltin_get_ent (args[0]);
  pause = (PLINT) class_double (e1);

  c_plspause (pause);

  return ent_Create_Rlab_Success();
}

/*
* Convert a RLaB matrix to a two-dimensional C array.
*/

// static PLFLT **
// convert_matrix_to_array (MDR * m)
// {
//   PLFLT **a;
//   int i, j, nrow, ncol;
//
//   nrow = MNR (m);
//   ncol = MNC (m);
//
//   /* Create the new array */
//   a = (PLFLT **) GC_MALLOC (nrow * sizeof (PLFLT *));
//   for (i = 0; i < nrow; i++)
//   {
//     a[i] = (PLFLT *) GC_malloc_atomic_ignore_off_page (ncol * sizeof (PLFLT));
//   }
//
//   /* Now load it up */
//   for (i = 0; i < nrow; i++)
//   {
//     for (j = 0; j < ncol; j++)
//     {
//       a[i][j] = Mdr0 (m, i, j);
//     }
//   }
//
//   return (a);
// }

/*
* Specify the desired pen width. The width must be
* between 1 and a device dependent maximum value.
*/

Ent *
_plot_plwid (int nargs, Datum args[])
{
  Ent *e1;
  PLINT width;

  if (nargs != 1)
  {
    rerror ("_plwid: 1 argument required");
  }

  e1 = bltin_get_ent (args[0]);
  width = (PLINT) class_double (e1);

#ifdef plwidth
  plwidth (width);
#else
  plwid (width);
#endif

  return ent_Create_Rlab_Success();
}

/*
* Writes text at a specified position and inclination within the viewport.
* Text is clipped at the viewport boundaries. The reference point of a string
* lies along a line passing through the string at half the height of a capital
* letter. The position of the reference point along this line is determined by
* just.
*/

Ent *
_plot_plptex (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *e3=0, *e4=0, *e5=0, *e6=0;
  PLFLT x, y, dx, dy, just;
  char *text = 0;

  if (nargs != 6)
    goto _exit_plptex;

  e1 = bltin_get_ent (args[0]);
  if (ent_type(e1)!=MATRIX_DENSE_REAL)
    goto _exit_plptex;

  e2 = bltin_get_ent (args[1]);
  if (ent_type(e2)!=MATRIX_DENSE_REAL)
    goto _exit_plptex;

  e3 = bltin_get_ent (args[2]);
  if (ent_type(e3)!=MATRIX_DENSE_REAL)
    goto _exit_plptex;

  e4 = bltin_get_ent (args[3]);
  if (ent_type(e4)!=MATRIX_DENSE_REAL)
    goto _exit_plptex;

  e5 = bltin_get_ent (args[4]);
  if (ent_type(e5)!=MATRIX_DENSE_REAL)
    goto _exit_plptex;

  e6 = bltin_get_ent (args[5]);
  if (ent_type(e6)!=MATRIX_DENSE_STRING)
    goto _exit_plptex;

  x = (PLFLT) class_double (e1);
  y = (PLFLT) class_double (e2);
  dx = (PLFLT) class_double (e3);
  dy = (PLFLT) class_double (e4);
  just = (PLFLT) class_double (e5);
  text = class_char_pointer (e6);

  if (isvalidstring(text))
    plptex (x, y, dx, dy, just, text);

_exit_plptex:

  ent_Clean (e1);
  ent_Clean (e2);
  ent_Clean (e3);
  ent_Clean (e4);
  ent_Clean (e5);
  ent_Clean (e6);

  return ent_Create_Rlab_Success();
}

/*
* Set the character set to use for subsequent characer drawing.
* May be called before initializing PLPLOT.
* 0 Standard character set.
* 1 Extended character set.
*/

Ent *
_plot_plfontld (int nargs, Datum args[])
{
  Ent *e1;
  PLINT set;

  if (nargs != 1)
  {
    rerror ("_plfontld: 1 argument required");
  }

  e1 = bltin_get_ent (args[0]);
  set = (PLINT) class_double (e1);

  plfontld (set);

  ent_Clean (e1);

  return ent_Create_Rlab_Success();
}

/*
* Sets the default characer font for subsequent character drawing.
* Also affects symbols produced by plpoin.
* 1 Normal font (simples and fastest)
* 2 Roman font
* 3 Italic font
* 4 Script font
*/

Ent *
_plot_plfont (int nargs, Datum args[])
{
  Ent *e1;
  PLINT font;

  if (nargs != 1)
  {
    rerror ("_plfont: 1 argument required");
  }

  e1 = bltin_get_ent (args[0]);
  font = (PLINT) class_double (e1);

  plfont (font);

  ent_Clean (e1);

  return ent_Create_Rlab_Success();
}

/*
* Set the current orientation.
* ori:  0 landscape
* ori:  1 portrait
* Not supported by all device drivers.
*/

Ent *
_plot_plsori (int nargs, Datum args[])
{
  Ent *e1;
  PLINT ori;

  if (nargs != 1)
  {
    rerror ("_plsori: 1 argument required");
  }

  e1 = bltin_get_ent (args[0]);
  ori = (PLINT) class_double (e1);

  plsori (ori);

  ent_Clean (e1);

  return ent_Create_Rlab_Success();
}

/*
* Globally turn color output on/off
*/

Ent *
_plot_plscolor (int nargs, Datum args[])
{
  Ent *e1;
  PLINT color;

  if (nargs != 1)
  {
    rerror ("_plscolor: 1 argument required");
  }

  e1 = bltin_get_ent (args[0]);
  color = (PLINT) class_double (e1);

  plscolor (color);

  ent_Clean (e1);

  return ent_Create_Rlab_Success();
}

/*
* Contour plotting.
*/

/*
* linear interpolation from singly dimensioned coord arrays.
*/

static int HasContour;

static void
plot_mapgrid (x, y, tx, ty, pltr_data)
    PLFLT x, y, *tx, *ty;
    PLPointer pltr_data;
{
  PLINT ul, ur, vl, vr;
  PLFLT du, dv;
  PLFLT xl, xr, yl, yr;

  PLcGrid *grid = (PLcGrid *) pltr_data;
  PLFLT *xg = grid->xg;
  PLFLT *yg = grid->yg;
  PLINT nx = grid->nx;
  PLINT ny = grid->ny;

  ul = (PLINT) x;
  ur = ul + 1;
  du = x - ul;

  vl = (PLINT) y;
  vr = vl + 1;
  dv = y - vl;

  if (x < 0 || x > nx - 1 || y < 0 || y > ny - 1)
  {
    rerror ("_mapgrid: Invalid coordinates");
  }

  xl = xg[ul];
  yl = yg[vl];

  if (ur == nx)
  {
    *tx = xl;
  }
  else
  {
    xr = xg[ur];
    *tx = xl * (1 - du) + xr * du;
  }
  if (vr == ny)
  {
    *ty = yl;
  }
  else
  {
    yr = yg[vr];
    *ty = yl * (1 - dv) + yr * dv;
  }
  HasContour = 1;
}

Ent *
_plot_plcont (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *e3, *e4, *e5, *e6, *e7, *e8;
  PLFLT *x, *y, *clevel;
  PLINT nx, ny, nzx, nzy, kx, lx, ky, ly, nlevel;
  PLcGrid *pltr_data;
  MDR *X, *Y, *clevelm, *zm;

  PLDLLIMPEXP PLFLT plf2eval( PLINT ix, PLINT iy, PLPointer plf2eval_data );

  if (nargs != 8)
  {
    rerror ("_plcont: 8 arguments required");
  }

  e1 = bltin_get_ent (args[0]);
  e2 = bltin_get_ent (args[1]);
  e3 = bltin_get_ent (args[2]);
  e4 = bltin_get_ent (args[3]);
  e5 = bltin_get_ent (args[4]);
  e6 = bltin_get_ent (args[5]);
  e7 = bltin_get_ent (args[6]);
  e8 = bltin_get_ent (args[7]);

  X = class_matrix_real (e1);
  Y = class_matrix_real (e2);
  x = (PLFLT *) MDRPTR (X);
  y = (PLFLT *) MDRPTR (Y);

  zm = class_matrix_real (e3);
  PLfGrid zp;
  zp.nx = nzx = (PLINT) MNR (zm);
  zp.ny = nzy = (PLINT) MNC (zm);
  zp.f = (PLFLT *) MDRPTR(zm);

  kx = (PLINT) class_double (e4);
  lx = (PLINT) class_double (e5);
  ky = (PLINT) class_double (e6);
  ly = (PLINT) class_double (e7);
  clevelm = class_matrix_real (e8);

  clevel = (PLFLT *) MDRPTR (clevelm);
  nlevel = MNR (clevelm) * MNC (clevelm);

  pltr_data = (PLcGrid *) GC_malloc (sizeof (PLcGrid));
  nx = (PLINT) MNR (X) * MNC (X);
  ny = (PLINT) MNR (Y) * MNC (Y);

  pltr_data->xg = x;
  pltr_data->nx = nx;
  pltr_data->yg = y;
  pltr_data->ny = ny;
  pltr_data->zg = 0;
  pltr_data->nz = 0;

  HasContour = 0;
  plfcont( plf2eval, (PLPointer) &zp, nzx, nzy, kx, lx, ky, ly, clevel, nlevel, plot_mapgrid, pltr_data);

  zp.f = 0;
  pltr_data->xg = 0;
  pltr_data->yg = 0;

  GC_FREE (pltr_data);

  ent_Clean (e1);
  ent_Clean (e2);
  ent_Clean (e3);
  ent_Clean (e4);
  ent_Clean (e5);
  ent_Clean (e6);
  ent_Clean (e7);
  ent_Clean (e8);

  return ent_Create_Rlab_Double(HasContour);
}

Ent *
_plot_plvpor (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *e3=0, *e4=0;
  PLFLT xmin, xmax, ymin, ymax;

  if ((nargs != 4) && (nargs != 1))
    goto _exit_plvpor;

  e1 = bltin_get_ent (args[0]);
  if (ent_type(e1)!=MATRIX_DENSE_REAL)
    goto _exit_plvpor;

  if (nargs == 1)
  {
    MDR *x = ent_data(e1);
    if (SIZE(x)!=4)
      goto _exit_plvpor;

    xmin = (PLFLT) mdrV0(x,0);
    xmax = (PLFLT) mdrV0(x,1);
    ymin = (PLFLT) mdrV0(x,2);
    ymax = (PLFLT) mdrV0(x,3);
  }
  else
  {
    e2 = bltin_get_ent (args[1]);
    if (ent_type(e2)!=MATRIX_DENSE_REAL)
      goto _exit_plvpor;

    e3 = bltin_get_ent (args[2]);
    if (ent_type(e3)!=MATRIX_DENSE_REAL)
      goto _exit_plvpor;

    e4 = bltin_get_ent (args[3]);
    if (ent_type(e4)!=MATRIX_DENSE_REAL)
      goto _exit_plvpor;

    xmin = (PLFLT) class_double (e1);
    xmax = (PLFLT) class_double (e2);
    ymin = (PLFLT) class_double (e3);
    ymax = (PLFLT) class_double (e4);
  }

  plvpor (xmin, xmax, ymin, ymax);

_exit_plvpor:

  ent_Clean (e1);
  ent_Clean (e2);
  ent_Clean (e3);
  ent_Clean (e4);

  return ent_Create_Rlab_Success();
}

Ent *
_plot_plvpas (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *e3, *e4, *e5;
  PLFLT xmin, xmax, ymin, ymax, aspect;

  if (nargs != 5)
  {
    rerror ("_plvpas: 5 arguments required");
  }

  e1 = bltin_get_ent (args[0]);
  e2 = bltin_get_ent (args[1]);
  e3 = bltin_get_ent (args[2]);
  e4 = bltin_get_ent (args[3]);
  e5 = bltin_get_ent (args[4]);

  xmin = (PLFLT) class_double (e1);
  xmax = (PLFLT) class_double (e2);
  ymin = (PLFLT) class_double (e3);
  ymax = (PLFLT) class_double (e4);
  aspect = (PLFLT) class_double (e5);

  plvpas (xmin, xmax, ymin, ymax, aspect);

  ent_Clean (e1);
  ent_Clean (e2);
  ent_Clean (e3);
  ent_Clean (e4);
  ent_Clean (e5);

  return ent_Create_Rlab_Success();
}

//
// _plerry() error bar plot,
//  first version by t.s.yang  6/14/94  uc berkeley
//  modified by mk to accept input
//    [x,y,y_min,y_max]
//  or
//    [x,y,y_sigma]
#undef THIS_SOLVER
#define THIS_SOLVER "_plerry"
Ent * _plot_plerry (int nargs, Datum args[])
{
  Ent *e1=0;
  int npoints, i;
  double rval = RLAB_STATUS_FAILURE;
  MDR *data=0;

  double *x=0, *y=0, *ymin=0, *ymax=0;

  if (nargs != 1)
  { rerror (THIS_SOLVER ": " RLAB_ERROR_ONE_ARG_REQUIRED "\n"); }

  e1 = bltin_get_ent (args[0]);
  data = ent_data(e1);
  if (!MD_TYPE_DOUBLE(data))
  { rerror (THIS_SOLVER ": " RLAB_ERROR_ARG1_MDR_MATRIX " with 3 or 4 columns\n"); }

  if (MNC(data)!=3 && MNC(data)!=4)
    goto _exit_plerry;

  npoints = MNR(data);
  if (npoints==0)
    goto _exit_plerry;

  // first column is 'x'
  x = &Mdr0(data,0,0);
  // second column is 'y'
  y = &Mdr0(data,0,1);

  if (MNC(data)==4)
  {
    // third column is 'y_min'
    ymin = &Mdr0(data,0,2);
    // fournt column is 'y_max'
    ymax = &Mdr0(data,0,3);
  }
  else
  {
    // compute third column 'y_min' from y and y_std
    ymin = (double *) GC_MALLOC(npoints * sizeof(double));
    // compute fourth column 'y_max' from y and y_std
    ymax = (double *) GC_MALLOC(npoints * sizeof(double));
    for (i=0; i<npoints; i++)
    {
      ymin[i] = y[i] - ABS(Mdr0(data,i,2));
      ymax[i] = y[i] + ABS(Mdr0(data,i,2));
    }
  }

  plerry (npoints, (PLFLT *) x, (PLFLT *) ymin, (PLFLT *) ymax);
  rval = RLAB_STATUS_SUCCESS;

  if (MNC(data)==3)
  {
    GC_FREE(ymin);
    GC_FREE(ymax);
  }

_exit_plerry:

  ent_Clean (e1);
  return ent_Create_Rlab_Double(rval);
}

//
// _plerrx() error bar plot,
// accepts input
//    [x,x_min,x_max,y]
//  or
//    [x,x_sigma,y]
//
#undef THIS_SOLVER
#define THIS_SOLVER "_plerrx"
Ent * _plot_plerrx (int nargs, Datum args[])
{
  Ent *e1=0;
  int npoints, i;
  double rval = RLAB_STATUS_FAILURE;
  MDR *data=0;

  double *x=0, *y=0, *xmin=0, *xmax=0;

  if (nargs != 1)
  { rerror (THIS_SOLVER ": " RLAB_ERROR_ONE_ARG_REQUIRED "\n"); }

  e1 = bltin_get_ent (args[0]);
  data = ent_data(e1);
  if (!MD_TYPE_DOUBLE(data))
  { rerror (THIS_SOLVER ": " RLAB_ERROR_ARG1_MDR_MATRIX " with 3 or 4 columns\n"); }

  if (MNC(data)!=3 && MNC(data)!=4)
    goto _exit_plerrx;

  npoints = MNR(data);
  if (npoints==0)
    goto _exit_plerrx;

  // first column is 'x'
  x = &Mdr0(data,0,0);
  // last column is 'y'
  y = &Mdr0(data,0,MNC(data)-1);

  if (MNC(data)==4)
  {
    // third column is 'y_min'
    xmin = &Mdr0(data,0,1);
    // fournt column is 'y_max'
    xmax = &Mdr0(data,0,2);
  }
  else
  {
    // compute third column 'y_min' from y and y_std
    xmin = (double *) GC_MALLOC(npoints * sizeof(double));
    // compute fourth column 'y_max' from y and y_std
    xmax = (double *) GC_MALLOC(npoints * sizeof(double));
    for (i=0; i<npoints; i++)
    {
      xmin[i] = x[i] - ABS(Mdr0(data,i,1));
      xmax[i] = x[i] + ABS(Mdr0(data,i,1));
    }
  }

  plerrx (npoints, (PLFLT *) xmin, (PLFLT *) xmax, (PLFLT *) y);
  rval = RLAB_STATUS_SUCCESS;

  if (MNC(data)==3)
  {
    GC_FREE(xmin);
    GC_FREE(xmax);
  }

_exit_plerrx:

  ent_Clean (e1);
  return ent_Create_Rlab_Double(rval);
}

//
//  _plerr(data, format, lt, pt) error bar plot for a single data set
//  accepts inputs:
//    [x, y]                              format ignored
//    [x, x_min, x_max, y]                "errorx"
//    [x, y, y_min, y_max]                "errory"
//    [x, x_min, x_max, y, y_min, y_max]  "errorxy"
#undef THIS_SOLVER
#define THIS_SOLVER "_plerr"
Ent * _plot_plerr (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *e3=0, *e4=0;
  int npoints, lt=-1, pt=-1, idummy;
  double rval = RLAB_STATUS_FAILURE;
  MDR *data=0;
  char *plfmt=0;

  double *x=0, *y=0, *ymin=0, *ymax=0, *xmin=0, *xmax=0;

  if (nargs != 2 && nargs != 3 && nargs != 4)
  { rerror (THIS_SOLVER ": " RLAB_ERROR_TWO_TO_FOUR_ARG_REQUIRED "\n"); }

  // get data
  e1 = bltin_get_ent (args[0]);
  data = ent_data(e1);
  if (!MD_TYPE_DOUBLE(data))
  { rerror (THIS_SOLVER ": " RLAB_ERROR_ARG1_MDR_MATRIX " with 2, 4 or 6 columns\n"); }
  if (MNC(data)!=2 && MNC(data)!=4 && MNC(data)!=6)
    goto _exit_plerr;

  // get description of data
  e2 = bltin_get_ent (args[1]);
  plfmt = class_char_pointer (e2);
  if (isvalidstring(plfmt)<1)
    goto _exit_plerr;

  if (nargs>2)
  {
    e3 = bltin_get_ent (args[2]);
    idummy = class_int(e3);
    if (idummy>=0)
      lt = idummy;
  }

  if (nargs>3)
  {
    e4 = bltin_get_ent (args[3]);
    idummy = class_int(e4);
    if (idummy>=0)
      pt = idummy;
  }

  // we need line or point or both
  if (lt==-1 && pt==-1)
    goto _exit_plerr;

  npoints = MNR(data);
  if (npoints==0)
    goto _exit_plerr;

  // first column is always 'x'
  x = &Mdr0(data,0,0);

  if (MNC(data) == 2)
  {
    y = &Mdr0(data,0,1);
  }
  else if (!strcmp(plfmt, "errorxy"))
  {
    if (MNC(data)!= 6)
      goto _exit_plerr;

    // second column is 'x_min'
    xmin = &Mdr0(data,0,1);
    // third column is 'x_max'
    xmax = &Mdr0(data,0,2);
    // third column is 'y'
    y = &Mdr0(data,0,3);
    // fourth column is 'y_min'
    ymin = &Mdr0(data,0,4);
    // sixth column is 'y_max'
    ymax = &Mdr0(data,0,5);
  }
  else if (!strcmp(plfmt, "errorx"))
  {
    if (MNC(data)!= 4)
      goto _exit_plerr;

    // second column is 'x_min'
    xmin = &Mdr0(data,0,1);
    // third column is 'x_max'
    xmax = &Mdr0(data,0,2);
    // third column is 'y'
    y = &Mdr0(data,0,3);
  }
  else if (!strcmp(plfmt, "errory"))
  {
    if (MNC(data)!= 4)
      goto _exit_plerr;

    // third column is 'y'
    y = &Mdr0(data,0,1);
    // fourth column is 'y_min'
    ymin = &Mdr0(data,0,2);
    // sixth column is 'y_max'
    ymax = &Mdr0(data,0,3);
  }

  if (xmin && xmax)
  {
    plerrx (npoints, (PLFLT *) xmin, (PLFLT *) xmax, (PLFLT *) y);
  }
  if (ymin && ymax)
  {
    plerry (npoints, (PLFLT *) x, (PLFLT *) ymin, (PLFLT *) ymax);
  }
  if (pt)
  {
    plpoin(npoints, (PLFLT *) x, (PLFLT *) y, pt);
  }
  if (lt)
  {
    plline(npoints, (PLFLT *) x, (PLFLT *) y);
  }
  rval = RLAB_STATUS_SUCCESS;

_exit_plerr:

  ent_Clean (e1);
  ent_Clean (e2);
  ent_Clean (e3);
  ent_Clean (e4);
  return ent_Create_Rlab_Double(rval);
}


/*
* _plschr(): Set the character height.
* _plschr (double def, double scale);
*/
#undef  THIS_SOLVER
#define THIS_SOLVER "_plchr"
Ent *
ent_plplot_plchr (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0;
  MDR *d=0, *rval=0;
  double d1=-1, d2=-1;

  if (nargs==1)
  {
    e1 = bltin_get_ent (args[0]);
    if (ent_type(e1)!=MATRIX_DENSE_REAL)
      goto _exit_plchr;

    d = ent_data(e1);
    if (SIZE(d)==2)
    {
      d1 = mdrV0(d,0);
      d2 = mdrV0(d,1);
      if (d1>=0 && d2>=0)
        plschr (d1, d2);
    }
    else
      goto _exit_plchr;
  }
  else if (nargs==2)
  {
    e1 = bltin_get_ent (args[0]);
    if (ent_type(e1)!=MATRIX_DENSE_REAL)
      goto _exit_plchr;
    e2 = bltin_get_ent (args[1]);
    if (ent_type(e2)!=MATRIX_DENSE_REAL)
      goto _exit_plchr;
    d1 = (double) class_double (e1);
    d2 = (double) class_double (e2);
    if (d1>=0 && d2>=0)
      plschr (d1, d2);
  }
  else if (nargs==0)
  {
    plgchr (&d1, &d2);
  }

  rval = mdr_Create(1,2);
  MdrV0(rval,0) = d1;
  if (d1 > 0)
  {
    d2 = d2 / d1;
  }
  if (d2==0)
  {
    d2 = 1.0;
  }
  MdrV0(rval,1) = d2;

_exit_plchr:

  ent_Clean (e1);
  ent_Clean (e2);

  return ent_Assign_Rlab_MDR(rval);
}

/*
* _plssym(): Set the symbol height.
* _plssym (double def, double scale);
*/
#undef  THIS_SOLVER
#define THIS_SOLVER "_plssym"
Ent * _plot_plssym (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0;
  MDR *d=0, *rval=0;
  double d1=-1, d2=-1;

  if (nargs==1)
  {
    e1 = bltin_get_ent (args[0]);
    if (ent_type(e1)!=MATRIX_DENSE_REAL)
      goto _exit_plssym;

    d = ent_data(e1);
    if (SIZE(d)==2)
    {
      d1 = mdrV0(d,0);
      d2 = mdrV0(d,1);
    }
    else
      goto _exit_plssym;
  }
  else if (nargs==2)
  {
    e1 = bltin_get_ent (args[0]);
    if (ent_type(e1)!=MATRIX_DENSE_REAL)
      goto _exit_plssym;
    e2 = bltin_get_ent (args[1]);
    if (ent_type(e2)!=MATRIX_DENSE_REAL)
      goto _exit_plssym;
    d1 = (double) class_double (e1);
    d2 = (double) class_double (e2);
  }

  if (d1>=0 && d2>=0)
    plssym (d1, d2);

  rval = mdr_Create(1,2);
  MdrV0(rval,0) = d1;
  MdrV0(rval,1) = d2;

_exit_plssym:

  ent_Clean (e1);
  ent_Clean (e2);

  return ent_Assign_Rlab_MDR(rval);
}

/*
* _plpsty(): Select one of the eight predefined area fill patterns to use.
* _plpsty (double n);
*/

Ent *
_plot_plpsty (int nargs, Datum args[])
{
  Ent *e1;
  PLINT n1;

  if (nargs != 1)
  {
    rerror ("_plpsty: 1 arguments required");
  }

  e1 = bltin_get_ent (args[0]);

  n1 = (PLINT) class_double (e1);

  plpsty (n1);

  ent_Clean (e1);

  return ent_Create_Rlab_Success();
}

//
// _plfill() fills the polygon with a pattarn in the color
//
#undef THIS_SOLVER
#define THIS_SOLVER "_plfill"
Ent * ent_plplot_plfill (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *e3=0;
  PLINT npoints=0, i, free_xy=0;
  MDR *xy=0;
  PLFLT *x=0, *y=0;
  FILE *rlab_stderr = (!RLAB_STDERR_DS) ? stderr : RLAB_STDERR_DS;

  if (nargs != 1 && nargs != 2 && nargs !=3)
  {
    fprintf(rlab_stderr, THIS_SOLVER ": " RLAB_ERROR_ONE_ARG_REQUIRED "\n");
    goto _exit_plfill;
  }

  e1 = bltin_get_ent (args[0]);
  if (ent_type(e1)!=MATRIX_DENSE_REAL)
    goto _exit_plfill;

  if (nargs>1)
  {
    e2 = bltin_get_ent (args[1]);
    if (ent_type(e2)==MATRIX_DENSE_REAL)
    {
      int ipat = class_int(e2);
      plpsty(ipat);
    }
  }

  if (nargs>2)
  {
    e3 = bltin_get_ent (args[2]);
    if (ent_type(e3)==MATRIX_DENSE_REAL)
    {
      int icol = class_int(e3);
      plcol0 (icol);
    }
  }

  xy = ent_data(e1);
  if ((MNC(xy)!=2) && (MNC(xy)!=4 || MNR(xy)!=1))
  {
    fprintf(rlab_stderr, THIS_SOLVER ": " RLAB_ERROR_ONE_ARG_REQUIRED "\n");
    fprintf(rlab_stderr, THIS_SOLVER ": " RLAB_ERROR_ARG1_MDR_MATRIX_2COL "\n");
    fprintf(rlab_stderr, THIS_SOLVER ": " RLAB_ERROR_ARG1_MDR_VECTOR_4COL "\n");
    goto  _exit_plfill;
  }

  if (MNR(xy) < 1)
  {
    goto _exit_plfill;
  }
  else if ((MNR(xy)==1) && (MNC(xy)==4))
  {
    // assume two opposite corners of the rectangle are given
    //  [x1,y1,x2,y2]
    // e.g., from multiplot command
    // create 4 points then:
    //  [x1,y1]
    //  [x2,y1]
    //  [x2,y2]
    //  [x1,y2]
    npoints = 4;
    x = GC_MALLOC(npoints * sizeof(PLFLT));
    y = GC_MALLOC(npoints * sizeof(PLFLT));
    x[0] = mdrV0(xy,0);
    x[1] = mdrV0(xy,2);
    x[2] = mdrV0(xy,2);
    x[3] = mdrV0(xy,0);
    y[0] = mdrV0(xy,1);
    y[1] = mdrV0(xy,1);
    y[2] = mdrV0(xy,3);
    y[3] = mdrV0(xy,3);
    free_xy = 1;
  }
  else
  {
    npoints = MNR(xy);

    if (MD_TYPE_INT32(xy))
    {
      x = GC_MALLOC(npoints * sizeof(PLFLT));
      y = GC_MALLOC(npoints * sizeof(PLFLT));
      for (i=0; i<npoints; i++)
      {
        x[i] = Mdi0(xy,i,0);
        y[i] = Mdi0(xy,i,1);
      }
      free_xy = 1;
    }
    else
    {
      x = &Mdr0(xy,0,0);
      y = &Mdr0(xy,0,1);
    }
  }

  plfill (npoints, (PLFLT *) x, (PLFLT *) y);

_exit_plfill:

  ent_Clean (e1);
  if (free_xy)
  {
    GC_FREE(x);
    GC_FREE(y);
  }
  ent_Clean (e2);
  ent_Clean (e3);
  return ent_Create_Rlab_Success();
}

//
//"#(728)"
//
#define RLAB_PLPLOT_NUM_HERSEYS 16
static char *hershey_string[RLAB_PLPLOT_NUM_HERSEYS] =
{ "#(841)",
  "#(210)", "#(225)", "#(228)", "#(840)",
  "#(227)", "#(841)", "#(842)", "#(904)",
  "#(840)", "#(735)", "#(743)", "#(844)",
  "#(841)", "#(866)", "#(868)"
};

#undef THIS_SOLVER
#define THIS_SOLVER "_pllegend"
Ent * ent_plplot_pllegend (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *e3=0, *e4=0, *e5=0, *e6=0, *e7=0, *e8=0, *e9=0;
  PLFLT legend_width, legend_height;

  char *copt=0;
  PLINT opt=PL_LEGEND_BACKGROUND | PL_LEGEND_BOUNDING_BOX | PL_LEGEND_TEXT_LEFT;
  char *cpos=0;
  PLINT pos=0, bg_color=0, bb_color=1, bb_style=1, i, nlegend;
  FILE *rlab_stderr = (!RLAB_STDERR_DS) ? stderr : RLAB_STDERR_DS;
  MDR *pos_xy=0;
  PLFLT pos_x=0, pos_y=0, plot_width=0.1;

  MDS *legend=0;
  MDR *legend_opts=0;

  PLINT *opt_array=0, *text_colors=0, *line_colors=0, *line_styles=0;
  PLFLT *line_widths=0, *symbol_scales=0;
  PLINT *symbol_colors=0, *symbol_numbers=0, nsymbols=0;

  PLFLT text_offset=1.0;
  PLFLT text_scale=1.0;
  PLFLT text_spacing=2.0;
  PLFLT text_justification=1.0;

  MDS *symbols=0;

  if (nargs < 2)
  {
    fprintf(rlab_stderr, THIS_SOLVER ": " RLAB_ERROR_AT_LEAST_TWO_ARG_REQUIRED "\n");
    fprintf(rlab_stderr, THIS_SOLVER ": " RLAB_WARNING_IGNORE "\n");
    goto _exit_pllegend;
  }

  //
  // legend is provided as a vector of length nlegend and as a matrix 10 cols by nlegend rows
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type(e1)!=MATRIX_DENSE_STRING)
  {
    fprintf(rlab_stderr, THIS_SOLVER ": " RLAB_ERROR_ARG1_MDS_VECTOR "\n");
    fprintf(rlab_stderr, THIS_SOLVER ": " RLAB_WARNING_IGNORE "\n");
    goto _exit_pllegend;
  }
  legend = ent_data(e1);
  nlegend = SIZE(legend);
  if (nlegend < 1)
  {
    fprintf(rlab_stderr, THIS_SOLVER ": " RLAB_ERROR_ARG1_MDS_VECTOR "\n");
    fprintf(rlab_stderr, THIS_SOLVER ": " RLAB_WARNING_IGNORE "\n");
    goto _exit_pllegend;
  }

  //
  // how to plot it
  //
  e2 = bltin_get_ent (args[1]);
  if (ent_type(e2)!=MATRIX_DENSE_REAL)
  {
    fprintf(rlab_stderr, THIS_SOLVER ": " RLAB_ERROR_ARG2_MDR_MATRIX "\n");
    fprintf(rlab_stderr, THIS_SOLVER ": " RLAB_WARNING_IGNORE "\n");
    goto _exit_pllegend;
  }
  legend_opts = ent_data(e2);
  if ((MNC(legend_opts)!=10) || (MNR(legend_opts)!=nlegend))
  {
    fprintf(rlab_stderr, THIS_SOLVER ": " RLAB_ERROR_ARG2_MDR_MATRIX "\n");
    fprintf(rlab_stderr, THIS_SOLVER ": " RLAB_WARNING_IGNORE "\n");
    goto _exit_pllegend;
  }

  //
  // text legend scale
  //
  if (nargs > 2)
  {
    e3 = bltin_get_ent (args[2]);
    if (ent_type(e3)==MATRIX_DENSE_REAL)
    {
      MDR * text_data = ent_data(e3);
      if (SIZE(text_data)>0)
        text_scale = mdrV0(text_data,0);
      if (SIZE(text_data)>1)
        text_spacing = mdrV0(text_data,1);
      if (SIZE(text_data)>2)
        text_justification = mdrV0(text_data,2);
      if (SIZE(text_data)>3)
        text_offset = mdrV0(text_data,3);
    }
  }

  //
  // position
  //
  if (nargs > 3)
  {
    e4 = bltin_get_ent (args[3]);
    cpos = class_char_pointer(e4);
    if (isvalidstring(cpos)<1)
      cpos=0;
    if (cpos)
    {
      pos = 0;
      if (strstr(cpos, "l")!=0)
        pos |= PL_POSITION_LEFT;
      if (strstr(cpos, "r")!=0)
        pos |= PL_POSITION_RIGHT;
      if (strstr(cpos, "t")!=0)
        pos |= PL_POSITION_TOP;
      if (strstr(cpos, "b")!=0)
        pos |= PL_POSITION_BOTTOM;
      if (strstr(cpos, "i")!=0)
        pos |= PL_POSITION_INSIDE;
      if (strstr(cpos, "o")!=0)
        pos |= PL_POSITION_OUTSIDE;
      if (strstr(cpos, "s")!=0)
        pos |= PL_POSITION_SUBPAGE;
      if (strstr(cpos, "v")!=0)
        pos |= PL_POSITION_VIEWPORT;
    }
  }

  //
  // position x,y
  //
  if (nargs > 4)
  {
    e5 = bltin_get_ent (args[4]);
    if (ent_type(e5)==MATRIX_DENSE_REAL)
    {
      pos_xy = ent_data(e5);
      if (SIZE(pos_xy)==2)
      {
        pos_x = mdrV0(pos_xy,0);
        pos_y = mdrV0(pos_xy,1);
      }
    }
  }

  //
  // options that control the overall legend:
  //  b - for bounding box
  //  s - for semitransparent background
  //  l - legend on left
  if (nargs > 5)
  {
    e6 = bltin_get_ent (args[5]);
    copt = class_char_pointer(e6);
    if (isvalidstring(copt)<1)
      copt=0;
    if (copt)
    {
      opt = 0;
      if (strstr(copt, "b")!=0)
        opt |= PL_LEGEND_BOUNDING_BOX;
      if (strstr(copt, "s")!=0)
        opt |= PL_LEGEND_BACKGROUND;
      if (strstr(copt, "l")!=0)
        opt |= PL_LEGEND_TEXT_LEFT;
    }
  }

  //
  // width of the plot legend area
  //
  if (nargs > 6)
  {
    e7 = bltin_get_ent (args[6]);
    if (ent_type(e7)==MATRIX_DENSE_REAL)
    {
      plot_width = class_double(e7);
    }
  }

  //
  // background color
  //
  if (nargs > 6)
  {
    e7 = bltin_get_ent (args[6]);
    if (ent_type(e7)==MATRIX_DENSE_REAL)
    {
      bg_color = class_int(e7);
    }
  }

  //
  // bounding box color
  //
  if (nargs > 7)
  {
    e8 = bltin_get_ent (args[7]);
    if (ent_type(e8)==MATRIX_DENSE_REAL)
    {
      bb_color = class_int(e8);
    }
  }

  //
  // bounding box style
  //
  if (nargs > 8)
  {
    e9 = bltin_get_ent (args[8]);
    if (ent_type(e9)==MATRIX_DENSE_REAL)
    {
      bb_style = class_int(e9);
    }
  }

  opt_array = GC_MALLOC(nlegend * sizeof(PLINT));
  text_colors = GC_MALLOC(nlegend * sizeof(PLINT));
  line_colors = GC_MALLOC(nlegend * sizeof(PLINT));
  line_styles = GC_MALLOC(nlegend * sizeof(PLINT));
  line_widths = GC_MALLOC(nlegend * sizeof(PLFLT));
  symbol_colors = GC_MALLOC(nlegend * sizeof(PLINT));
  symbol_scales = GC_MALLOC(nlegend * sizeof(PLFLT));
  symbol_numbers = GC_MALLOC(nlegend * sizeof(PLINT));
  symbols = mds_Create(1,nlegend);

  for (i=0; i<nlegend; i++)
  {
    // line or point or both?
    if (mdi0(legend_opts,i,RLABPLUS_PLPLOT_IDX_LINEPOINT)==0)
    {
      opt_array[i] = PL_LEGEND_LINE;
      nsymbols = 0;
    }
    else if (mdi0(legend_opts,i,RLABPLUS_PLPLOT_IDX_LINEPOINT)==1)
    {
      nsymbols = 3;
      opt_array[i] = PL_LEGEND_SYMBOL;
      int hershey_idx = mdi0(legend_opts,i,RLABPLUS_PLPLOT_IDX_POINTTYPE) % RLAB_PLPLOT_NUM_HERSEYS;
      MdsV0(symbols,i) = cpstr( hershey_string[hershey_idx] );
    }
    else
    {
      nsymbols = 2;
      opt_array[i] = PL_LEGEND_LINE | PL_LEGEND_SYMBOL;
      int hershey_idx = mdi0(legend_opts,i,RLABPLUS_PLPLOT_IDX_POINTTYPE) % RLAB_PLPLOT_NUM_HERSEYS;
      MdsV0(symbols,i) = cpstr( hershey_string[hershey_idx] );
    }

    // color of the text
    text_colors[i] = 1; /* black for now */

    // line colors
    line_colors[i] = mdi0(legend_opts,i,RLABPLUS_PLPLOT_IDX_LINECOLOR);
    // line colors
    line_styles[i] = mdi0(legend_opts,i,RLABPLUS_PLPLOT_IDX_LINETYPE);
    // line colors
    line_widths[i] = mdr0(legend_opts,i,RLABPLUS_PLPLOT_IDX_LINEWIDTH);

    // symbol colors
    symbol_colors[i] = mdi0(legend_opts,i,RLABPLUS_PLPLOT_IDX_POINTCOLOR);
    // symbol size
    symbol_scales[i] = mdr0(legend_opts,i,RLABPLUS_PLPLOT_IDX_POINTSIZE);
    // symbol type
    symbol_numbers[i] = nsymbols;
  }

  pllegend( &legend_width, &legend_height, opt, pos, pos_x, pos_y, plot_width,
             bg_color, bb_color, bb_style, 0, 0,
             nlegend, opt_array,
             text_offset, text_scale, text_spacing, text_justification,
             text_colors, (const char **) MDPTR(legend),
             NULL, NULL, NULL, NULL,
             line_colors, line_styles, line_widths,
             symbol_colors, symbol_scales, symbol_numbers, (const char **) MDPTR(symbols));

_exit_pllegend:

  ent_Clean (e1);
  ent_Clean (e2);
  ent_Clean (e3);
  ent_Clean (e4);
  ent_Clean (e5);
  ent_Clean (e6);
  ent_Clean (e7);
  ent_Clean (e8);
  ent_Clean (e9);

  if (opt_array)
    GC_FREE (opt_array);
  if (text_colors)
    GC_FREE (text_colors);
  if (line_colors)
    GC_FREE (line_colors);
  if (line_styles)
    GC_FREE (line_styles);
  if (line_widths)
    GC_FREE (line_widths);
  if (symbol_colors)
    GC_FREE (symbol_colors);
  if (symbol_scales)
    GC_FREE (symbol_scales);
  if (symbol_numbers)
    GC_FREE (symbol_numbers);

  mds_Destroy(symbols);

  return ent_Create_Rlab_Success();
}

#undef THIS_SOLVER
#define THIS_SOLVER "_plsetopts"
Ent * ent_plplot_plsetopt (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0;
  if (nargs != 2)
    goto _exit_plsetopts;

  e1 = bltin_get_ent (args[0]);
  if (ent_type(e1)!=MATRIX_DENSE_STRING)
    goto _exit_plsetopts;

  e2 = bltin_get_ent (args[1]);
  if (ent_type(e2)!=MATRIX_DENSE_STRING)
    goto _exit_plsetopts;

  char *opt = class_char_pointer(e1);
  if (isvalidstring(opt)<1)
    goto _exit_plsetopts;

  char *optarg = class_char_pointer(e2);
  if (isvalidstring(optarg)<1)
    goto _exit_plsetopts;

  plsetopt (opt, optarg);

_exit_plsetopts:

  ent_Clean (e1);
  ent_Clean (e2);
  return ent_Create_Rlab_Success();
}


#endif /* HAVE_RLAB_PLPLOT */
