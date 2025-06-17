// This file is a part of RLaB ("Our"-LaB)
// Copyright (C) 2003 Marijan Kostrun
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

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <time.h>

// naming convention for the solver parameters
#include "rlab_solver_parameters_names.h"

static char
chaos_find_nearest (MDR * s, MDR * box, MDR * list, long n, unsigned int dim,
                    double eps,
                    int ibox,
                    unsigned long theiler,
                    double *aveps,
                    double *vareps,
                    unsigned long *toolarge,
                    unsigned int *delay, double varianz, double rt)
{
  // Author: Rainer Hegger Last modified: March 1st, 1998
  // rlab modification by Marijan Kostrun, 2004
  int x, y, x1, x2, y1, i, i1;
  long element, which = -1;
  double dx, maxdx, mindx = 1.1, factor;

  x = ((int) (Mdr0 (s, n - (dim - 1) * (*delay), 0) / eps)) & ibox;
  y = ((int) (Mdr0 (s, n, 0) / eps)) & ibox;

  for (x1 = x - 1; x1 <= x + 1; x1++)
  {
    x2 = x1 & ibox;
    for (y1 = y - 1; y1 <= y + 1; y1++)
    {
      element = Mdr0 (box, x2, y1 & ibox);
      while (element != -1)
      {
        if (labs (element - n) > theiler)
        {
          maxdx = fabs (Mdr0 (s, n, 0) - Mdr0 (s, element, 0));
          for (i = 1; i < dim; i++)
          {
            i1 = i * (*delay);
            dx = fabs (Mdr0 (s, n - i1, 0) - Mdr0 (s, element - i1, 0));
            if (dx > maxdx)
              maxdx = dx;
          }
          if ((maxdx < mindx) && (maxdx > 0.0))
          {
            which = element;
            mindx = maxdx;
          }
        }
        element = Mdr0 (list, element, 0);
      }
    }
  }
  if ((which != -1) && (mindx <= eps) && (mindx <= varianz / rt))
  {
    *aveps += mindx;
    *vareps += mindx * mindx;
    factor = fabs (Mdr0 (s, n + 1, 0) - Mdr0 (s, which + 1, 0)) / mindx;
    if (factor > rt)
      (*toolarge)++;
    return 1;
  }
  return 0;
}


static void
chaos_make_box (
		 // Author: Rainer Hegger Last modified: March 1st, 1998
		 // rlab modification by Marijan Kostrun, 2004
		 MDR * s, MDR * box, MDR * list, unsigned int dim,
		 unsigned int del, double eps)
{
  int i, x, y, bs, ib, n;

  bs = MNR (box);		// = MNC(box);
  ib = bs - 1;
  n = MNR (s);

  for (x = 0; x < bs; x++)
    for (y = 0; y < bs; y++)
      Mdr0 (box, x, y) = -1;
  for (i = (dim - 1) * del; i < n - 1; i++)
  {
    x = (int) (Mdr0 (s, i - (dim - 1) * del, 0) / eps) & ib;
    y = (int) (Mdr0 (s, i, 0) / eps) & ib;
    Mdr0 (list, i, 0) = (int) Mdr0 (box, x, y);
    Mdr0 (box, x, y) = i;
  }
}


static void
chaos_variance (MDR * x, double *av, double *sd)
{
  unsigned long i, n;
  double h;

  n = MNR (x);
  *av = *sd = 0.0;
  for (i = 0; i < n; i++)
  {
    h = Mdr0 (x, i, 0);
    *av += h;
    *sd += h * h;
  }
  *av /= (double) n;
  *sd = sqrt (fabs ((*sd) / (double) n - (*av) * (*av)));
  if (*sd == 0.0)
    rerror
      ("Terrible internal error: variance of the data is zero. Exiting!\n\n");
}

MDR *
chaos_rescale_data (MDR * x, double *smin, double *sint)
{
  // Author: Rainer Hegger Last modified: March 1st, 1998
  // rlab modification by Marijan Kostrun, 2004
  int i, xr;
  double smax;
  MDR *w;
  *smin = smax = MdrV1 (x, 1);
  xr = MNR (x) * MNC (x);
  w = mdr_Create (xr, 1);
  for (i = 2; i <= xr; i++)
  {
    if (MdrV1 (x, i) < *smin)
      *smin = MdrV1 (x, i);
    if (Mdr1 (x, i, 1) > smax)
      smax = MdrV1 (x, i);
  }
  *sint = smax - *smin;
  if (*sint != 0.0)
    for (i = 1; i <= xr; i++)
      MdrV1 (w, i) = (MdrV1 (x, i) - *smin) / *sint;
  else
    rerror ("ami: Terible internal error: zero entropy time series!\n");
  return w;
}


static double
chaos_make_cond_entropy (long t, MDR * rx, MDR * h1, MDR * h2, MDR * h11, int p)
{
  // Author: Rainer Hegger Last modified: March 1st, 1998
  // rlab modification by Marijan Kostrun, 2004
  long i, j, hi, hii, count = 0, n;
  double hpi, hpj, pij, cond_ent = 0.0, norm;

  n = MNR (rx);

  mdr_Zero (h1);
  mdr_Zero (h11);
  mdr_Zero (h2);

  for (i = t + 1; i <= n; i++)
  {
    hii = Mdr1 (rx, i, 1);
    hi = Mdr1 (rx, i - t, 1);
    Mdr1 (h1, hi, 1) += 1;
    Mdr1 (h11, hii, 1) += 1;
    Mdr1 (h2, hi, hii) += 1;
    count++;
  }

  norm = 1.0 / (double) count;
  cond_ent = 0.0;

  for (i = 1; i <= p; i++)
  {
    hpi = Mdr1 (h1, i, 1) * norm;
    if (hpi > 0.0)
    {
      for (j = 1; j <= p; j++)
      {
	hpj = Mdr1 (h11, j, 1) * norm;
	if (hpj > 0.0)
	{
	  pij = Mdr1 (h2, i, j) * norm;
	  if (pij > 0.0)
	    cond_ent += pij * log (pij / hpj / hpi);
	}
      }
    }
  }
  return cond_ent;
}


//
// false nearest neighbours
// Author: Rainer Hegger Last modified: March 1st, 1998
// rlab modification by Marijan Kostrun, 2004, 2006
//
Ent *
ent_falsenn (int nargs, Datum args[])
{
  Ent *X=0, *ED=0, *D=0, *TW=0, *RT=0, *FOUT=0, *rent;
  MDR *x, *ed, *rx, *w, *list, *box;

  int BOX = 1024, ibox = BOX - 1;
  int n, neds;
  double smin, sint, epsilon, av, eps0 = 1.0e-5, rt = 5.0;
  double aveps, vareps, varianz, rtime;
  char *nearest, alldone, fnnstderr[256] = { '\0' };
  //long i, *list;
  long i, j;
  unsigned int dim, delay;
  unsigned long theiler = 0.0, donesofar, toolarge;
  FILE *fptr = NULL;
  time_t t1, t2;

  if (nargs < 2)
  {
    printf
      ("falsenn: False nearest neighbor statistics of an uniformly sampled\n");
    printf ("falsenn: time series given as a column-vector. Format:\n");
    printf ("falsenn:   falsenn(x, [e1..eN], d /, fout, escfac, theiler/ )\n");
    printf ("falsenn:     x        - 1d time series,\n");
    printf ("falsenn:     [e1..eN] - embedding dimensions of interest,\n");
    printf ("falsenn:     d        - delay,\n");
    printf ("falsenn:     fout     - name of the output stream,\n");
    printf ("falsenn:     escfac   - escape factor,\n");
    printf ("falsenn:     theiler  - size of the theiler window,\n");
    rerror ("falsenn: wrong arguments !");
  }

  //
  // get the time series in 'x'
  //
  X = bltin_get_ent (args[0]);
  if (ent_type (X) != MATRIX_DENSE_REAL)
    rerror ("falsenn: 'x' is a real vector, the time series!");
  x = class_matrix_real (X);
  if (MNC (x) != 1 && MNR (x) != 1)
    rerror ("falsenn: 'x' is a real vector, the time series!");

  //
  // get the range of embedded dimensions
  //
  ED = bltin_get_ent (args[1]);
  if (ent_type (ED) != MATRIX_DENSE_REAL)
    rerror ("falsenn: Argument 'embdim' must be vector!");
  ed = class_matrix_real (ED);
  if (!ed)
    rerror ("falsenn: Argument 'embdim' must be vector!");
  if (MNR (ed) != 1 && MNC (ed) != 1)
    rerror ("falsenn: 'embdim' is a row-vector with range of interest!");

  // get the delay array 'd'
  D = bltin_get_ent (args[2]);
  if (ent_type (D) != MATRIX_DENSE_REAL)
    rerror ("falsenn: 'delay' is a scalar!");
  delay = class_double (D);

  if (nargs > 3)
  {
    FOUT = bltin_get_ent (args[3]);
    if (ent_type (FOUT) == MATRIX_DENSE_STRING)
      sprintf (fnnstderr, "%s", MdsV0((MDS*)ent_data(FOUT),0));
    else
      fnnstderr[0] = '\0';
  }
  if (strlen (fnnstderr) > 1)
    fptr = fopen (fnnstderr, "a");
  if (nargs > 4)
  {
    TW = bltin_get_ent (args[4]);
    if (ent_type (TW) == MATRIX_DENSE_REAL)
      theiler = (long) class_double (TW);
  }
  if (nargs > 5)
  {
    RT = bltin_get_ent (args[5]);
    if (ent_type (RT) == MATRIX_DENSE_REAL)
      rt = class_double (RT);
  }

  neds = MNC (ed) * MNR (ed);
  w = mdr_Create (1, neds);
  mdr_Zero (w);
  n = MNR (x) * MNC (x);
  rx = mdr_Create (n, 1);
  list = mdr_Create (n, 1);
  box = mdr_Create (BOX, BOX);
  nearest = (char *) GC_malloc (n);

  t1 = clock ();

  rx = chaos_rescale_data (x, &smin, &sint);
  chaos_variance (rx, &av, &varianz);
  //list    = (long*)  malloc(sizeof(long)*length);
  //nearest = (char*)  malloc(n);
  //box     = (long**) malloc(sizeof(long*)*BOX);
  //for (i=0;i<BOX;i++)
  //    box[i]=(long*)malloc(sizeof(long)*BOX);


  for (j = 0; j < neds; j++)
  {
    dim = (int) MdrV0 (ed, j);
    epsilon = eps0;
    toolarge = 0;
    alldone = 0;
    donesofar = 0;
    aveps = 0.0;
    vareps = 0.0;
    for (i = 0; i < n; i++)
      nearest[i] = 0;
    if (fptr)
      fprintf (fptr, "falsenn: embedded dimension = %u\n", dim);
    while (!alldone && (epsilon < 2. * varianz / rt))
    {
      //fprintf(stdout,"."); fflush(stdout);
      alldone = 1;
      chaos_make_box (rx, box, list, dim, delay, epsilon);
      for (i = (dim - 1) * delay; i < n - 1; i++)
	if (!nearest[i])
	{
	  nearest[i] =
	    chaos_find_nearest (rx, box, list, i, dim, epsilon, ibox, theiler,
				&aveps, &vareps, &toolarge, &delay, varianz,
				rt);
	  alldone &= nearest[i];
	  donesofar += (unsigned long) nearest[i];
	  //fprintf(stdout, "nearest   = %u\n", (unsigned long)nearest[i]);
	  //fprintf(stdout, "donesofar = %u\n", (unsigned long)donesofar);
	}
      if (fptr)
        fprintf (fptr, "falsenn: found %lu points for epsilon = %16.12f\n",
                 donesofar, epsilon * sint);
      epsilon *= sqrt (2.0);
      if (!donesofar)
        eps0 = epsilon;
    }
    if (donesofar == 0)
    {
      if (fptr)
      {
      	fprintf (fptr, "falsenn: Not enough points found in phase space!\n");
	       fflush (fptr);
      }
      aveps = 0;
      vareps = 0;
      Mdr0 (w, 0, j) = 1;
    }
    else
    {
      aveps *= (1. / (double) donesofar);
      vareps *= (1. / (double) donesofar);
      Mdr0 (w, 0, j) = (double) toolarge / (double) donesofar;
      if (fptr)
      {
	fprintf (fptr, "falsenn: search results for dim = %u", dim);
	fprintf (fptr, " is fnn = %12.6f", Mdr0 (w, 0, j));
	fprintf (fptr, " with epsilon = %f", aveps);
	fprintf (fptr, " +/- %f\n", vareps);
	fflush (fptr);
      }
    }
  }

  t2 = clock ();
  rtime = t2 - t1;
  rtime /= 1e6;

  if (fptr)
  {
    fprintf (fptr, "falsenn: Calculation lasted %lg sec.\n", rtime);
    fclose (fptr);
  }

  // clean up a little
  GC_free (nearest);
  mdr_Destroy (rx);
  mdr_Destroy (list);
  mdr_Destroy (box);


  ent_Clean (X);
  ent_Clean (ED);
  ent_Clean (D);
  ent_Clean (FOUT);
  ent_Clean (TW);
  ent_Clean (RT);

  rent = ent_Create ();
  ent_data (rent) = w;
  ent_SetType (rent, MATRIX_DENSE_REAL);
  return (rent);
}


//
// average mutual information
// Author: Rainer Hegger Last modified: March 1st, 1998
// rlab modification by Marijan Kostrun, 2004
//
Ent *
ent_avgmutinfo (int nargs, Datum args[])
{
  Ent *X=0, *D=0, *P=0, *rent;
  MDR *w=0, *x=0, *d=0, *h1=0, *h11=0, *h2=0, *rx=0;
  int i, k, n, dc;
  double smin, sint, p = 16.0;
  if (nargs < 2)
  {
    fprintf (stdout,
             "ami: Average mutual information of an uniformly sampled time series.\n");
    fprintf (stdout,
             "ami: Format:\n");
    fprintf (stdout,
             "ami:   y = ami(x, d /, p/)\n");
    fprintf (stdout,
             "ami: where 'x' is the time series, 'd' is a vector of embedding\n");
    fprintf (stdout,
             "ami: dimensions, and 'p' is the number of partitions of a range\n");
    fprintf (stdout,
             "ami: of 'x'.\n");
    rerror ("three arguments required!");
  }

  // get the time series in 'x'
  X = bltin_get_ent (args[0]);
  if (ent_type (X) == MATRIX_DENSE_REAL)
    x = class_matrix_real (X);
  else
    rerror ("ami: 'x' is a real column-vector, the time series!");
  if (MNC (x) != 1 && MNR (x) != 1)
    rerror ("ami: 'x' is a real column-vector, the time series!");

  // get the delay array 'd'
  D = bltin_get_ent (args[1]);
  if (ent_type (D) == MATRIX_DENSE_REAL)
    d = class_matrix_real (D);
  else
    rerror ("ami: 'range' is a vector of dimensions of interest!");
  if (MNR (d) != 1 && MNC (d) != 1)
    rerror ("ami: 'range' is a vector of dimensions of interest!");

  if (nargs > 2)
  {
    P = bltin_get_ent (args[2]);
    if (ent_type (P) == MATRIX_DENSE_REAL)
      p = class_double (P);
    else
      rerror ("ami: 'p' is a number of partitions!");
  }

  dc = MNC (d) * MNR (d);
  n  = MNR (x) * MNC (x);
  rx = mdr_Create (n, 1);
  w  = mdr_Create (1, dc);

  //
  // preparing the bins for counting
  //
  h1  = mdr_Create (p, 1);
  h11 = mdr_Create (p, 1);
  h2  = mdr_Create (p, p);

  //
  // rescale the data and put it in the bins 1..p
  //
  rx = chaos_rescale_data (x, &smin, &sint);

  for (i = 1; i <= n; i++)
    Mdr1 (rx, i, 1) = (int) (Mdr1 (rx, i, 1) * p + 1);

  //
  // doing the runs over delays in 'd'
  //
  for (k = 1; k <= dc; k++)
    MdrV1 (w, k) =
      chaos_make_cond_entropy (MdrV1 (d, k), rx, h1, h2, h11, p);

  // clean-up a little bit
  mdr_Destroy (rx);
  mdr_Destroy (h1);
  mdr_Destroy (h11);
  mdr_Destroy (h2);

  ent_Clean (X);
  ent_Clean (D);
  ent_Clean (P);

  rent = ent_Create ();
  ent_data (rent) = w;
  ent_SetType (rent, MATRIX_DENSE_REAL);
  return (rent);
}



//
// running variance
//
Ent *
ent_runvar (int nargs, Datum args[])
{
  Ent *e1, *e2, *rent;
  MDR *w, *x, *d;
  int i, k, n, dc, del;
  double sy, sy2, rdel, mu;
  if (nargs < 2)
  {
    fprintf (stdout,
             "runvar Running variance of an uniformly sampled time series.\n");
    fprintf (stdout,
             "runvar: Format:\n");
    fprintf (stdout,
             "runvar:   v = runvar( x, T )\n");
    fprintf (stdout,
             "runvar: where 'x' is the time series and 'T' is the averaging\n");
    fprintf (stdout,
             "runvar: period.\n");
    rerror ("runvar: two arguments required!");
  }

  // get the time series in 'x'
  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) != MATRIX_DENSE_REAL)
    rerror ("runvar: First argument 'x' must be real vector");
  x = ent_data (e1);
  if (!x)
    rerror ("runvar: First argument 'x' must be real vector");
  if (MNC (x) != 1 && MNR (x) != 1)
    rerror ("runvar: First argument 'x' must be real vector");

  // get the delay array 'd'
  e2 = bltin_get_ent (args[1]);
  if (ent_type (e2) != MATRIX_DENSE_REAL)
    rerror ("runvar: Second argument 'd' must be real vector");
  d = ent_data (e2);
  if (!d)
    rerror ("runvar: Second argument 'd' must be real vector");
  if (MNC (d) != 1 && MNR (d) != 1)
    rerror ("runvar: Second argument 'd' must be real vector");

  dc = MNC (d) * MNR (d);
  n  = MNR (x) * MNC (x);
  w = mdr_Create (n, dc);

  for (k = 0; k < dc; k++)
  {
    sy  = 0.0;
    sy2 = 0.0;
    if (d->type == RLAB_TYPE_INT32)
    {
      del  = (double) MdiV0 (d, k);
      rdel = 1.0 / ((double) MdiV0 (d, k));
    }
    else
    {
      del  = MdrV0 (d, k);
      rdel = 1.0 / MdrV0 (d, k);
    }

    for (i = 0; i < del; i++)
    {
      Mdr0 (w, i, k) = 0;
      if (x->type == RLAB_TYPE_INT32)
      {
        sy  += (double) MdiV0 (x, i);
        sy2 += (double)(MdiV0 (x, i) * MdiV0 (x, i));
      }
      else
      {
        sy  += MdrV0 (x, i);
        sy2 += MdrV0 (x, i) * MdrV0 (x, i);
      }
    }
    mu = sy  * rdel;
    Mdr0 (w, i-1, k) = sy2 * rdel - mu * mu;
    for (i = del ; i < n; i++)
    {
      if (x->type == RLAB_TYPE_INT32)
      {
        sy  += (double)(MdiV0 (x, i) - MdiV0 (x, i - del));
        sy2 += (double)(MdiV0 (x, i) * MdiV0 (x, i));
        sy2 -= (double)(MdiV0 (x, i - del) * MdiV0 (x, i - del));
      }
      else
      {
        sy  += MdrV0 (x, i) - MdrV0 (x, i - del);
        sy2 += MdrV0 (x, i) * MdrV0 (x, i);
        sy2 -= MdrV0 (x, i - del) * MdrV0 (x, i - del);
      }
      mu = sy  * rdel;
      Mdr0 (w, i, k) = sy2 * rdel - mu * mu;
    }
  }

  if (e1)
    if (ent_type(e1) != UNDEF)
      ent_Clean (e1);
  if (e2)
    if (ent_type(e2) != UNDEF)
      ent_Clean (e2);

  rent = ent_Create ();
  ent_data (rent) = w;
  ent_SetType (rent, MATRIX_DENSE_REAL);
  return (rent);
}

//
// cross-correlation function
//
#define THIS_SOLVER "xcorr"
Ent *
ent_xcorr (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *e3=0, *e4=0, *rent;
  MDR  *w=0,  *rx=0,  *ry=0,  *d=0;
  MDC *cw=0, *cx=0, *cy=0;
  int i, k, n, dc, del, have_nans=0;
  void *x=0, *y=0;
  double exy, vx, vy, xi, yj;
  Complex cexy, cxi, cyj;

  char *choi=0;

  int inorm = 1;  // no normalization

  if (nargs!=3 && nargs!=4)
  {
    fprintf(stderr, THIS_SOLVER ": cross-correlation function of two uniformly sampled time series\n");
    fprintf(stderr, THIS_SOLVER ": given as a column-vector. Format:\n");
    fprintf(stderr, THIS_SOLVER ":   xcorr( x, y, d, opt)\n");
    fprintf(stderr, THIS_SOLVER ": where\n");
    fprintf(stderr, THIS_SOLVER ":  x, y  are same-length data vectors, while\n");
    fprintf(stderr, THIS_SOLVER ":  d = [d1..dN]  are delays of interest,\n");
    fprintf(stderr, THIS_SOLVER ":  o = 'c', 'n', 'b', 'u' designates crosscorrelation being sought\n");
    rerror (THIS_SOLVER ": requires three or four arguments !");
  }

  // get the time series in 'x'
  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) == MATRIX_DENSE_REAL || ent_type (e1) == MATRIX_DENSE_COMPLEX)
    x = (void *) ent_data(e1);
  if (!x)
    rerror (THIS_SOLVER ": first argument 'x' must be real or complex vector");

  e2 = bltin_get_ent (args[1]);
  if (ent_type (e2) == MATRIX_DENSE_REAL || ent_type (e2) == MATRIX_DENSE_COMPLEX)
    y = (void *) ent_data(e2);
  if (!y)
    rerror (THIS_SOLVER ": second argument 'y' must be real or complex vector");

  if (ent_type (e1) != ent_type (e2))
    rerror (THIS_SOLVER ": first and second argument must be same type, real or complex");

  e3 = bltin_get_ent (args[2]);
  if (ent_type (e3) == MATRIX_DENSE_REAL)
    d = ent_data (e3);
  if (!d)
    rerror (THIS_SOLVER ": third argument must be real");

  // get normalization information
  if (nargs > 3)
  {
    e4 = bltin_get_ent (args[3]);
    choi = class_char_pointer(e4);

    if(choi)
    {
      if(choi[0] == 'c' || choi[0] == 'C')
        inorm = 0;  // coeff
      if(choi[0] == 'n' || choi[0] == 'N')
        inorm = 1;  // none
      else if(choi[0] == 'b' || choi[0] == 'B')
        inorm = 2;  // biased
      else if(choi[0] == 'u' || choi[0] == 'U')
        inorm = 3;  // unbiased
    }
  }

  if (ent_type (e1) == MATRIX_DENSE_REAL)
  {
    rx = (MDR *) x;
    ry = (MDR *) y;

    dc = d->nrow * d->ncol;
    w  = mdr_Create (1, dc);

    n  = rx->nrow * rx->ncol;
    if (n != (ry->nrow * ry->ncol))
      rerror (THIS_SOLVER ": 'x' and 'y' have to be same size");

    for (k = 0; k < dc; k++)
    {
      exy = 0.;
      vx  = 0.;
      vy  = 0.;

      del = mdiV0 (d, k);

      if (n < del)
      {
        MdrV0 (w, k) = create_nan();
        continue;
      }

      for (i = 0; i < n - del; i++)
      {
        // x[i]
        xi = mdrV0 (rx, i);
        if (isnand(xi))
        {
          have_nans++;
          continue;
        }

        // y[i+del]
        yj = mdrV0 (ry, i+del);
        if (isnand(yj))
        {
          have_nans++;
          continue;
        }

        // now do it
        exy += xi * yj;
        vx  += xi * xi;
        vy  += yj * yj;
      }

      if (inorm == 0)
        MdrV0 (w, k) = exy / sqrt (vx * vy);  // acf at 0 lag is unity
      else
      {
        MdrV0 (w, k) = exy;
        if (inorm == 2)
          MdrV0 (w, k) /= (n - have_nans);          // biased estimate
        else if (inorm == 3)
          MdrV0 (w, k) /= (n - have_nans - del);    // unbiased estimate
      }
    }
  }
  else if (ent_type (e1) == MATRIX_DENSE_COMPLEX)
  {
    cx = (MDC *) x;
    cy = (MDC *) y;

    dc = d->nrow * d->ncol;
    cw = mdc_Create (1, dc);

    n = cx->nrow * cx->ncol;
    if (n != (cy->nrow * cy->ncol))
      rerror (THIS_SOLVER ": 'x' and 'y' have to be same size");

    for (k = 0; k < dc; k++)
    {
      cexy = 0+0i;
      vx  = 0.;
      vy  = 0.;

      del = mdiV0 (d, k);

      if (n < del)
      {
        MdcV0(cw, k) = create_nan();
        continue;
      }

      for (i = 0; i < n - del; i++)
      {
        // x[i]
        cxi = MdcV0 (cx, i);
        if (isnand(RE(cxi)) || isnand(IM(cxi)))
        {
          have_nans++;
          continue;
        }

        // y[i+del]
        cyj = MdcV0 (cy, i+del);
        if (isnand(RE(cyj)) || isnand(RE(cyj)))
        {
          have_nans++;
          continue;
        }

        // now do it
        cexy = cexy + cxi * cyj;
        vx  += pow(cabs(cxi),2);
      }

      if (inorm == 0)
      {
        MdcV0(cw, k) = cexy / sqrt (vx*vy);  // acf at 0 lag is unity
      }
      else
      {
        MdcV0(cw, k) = cexy;
        if (inorm == 2)
        {
          MdcV0 (w, k) /= (double) (n - have_nans);         // biased estimate
        }
        else if (inorm == 3)
        {
          MdcV0(w, k) /= (double) (n - have_nans - del); // unbiased estimate
        }
      }
    }
  }

  ent_Clean (e1);
  ent_Clean (e2);
  ent_Clean (e3);
  ent_Clean (e4);

  rent = ent_Create ();
  if (w)
  {
    ent_data (rent) = w;
    ent_SetType (rent, MATRIX_DENSE_REAL);
  }
  else if (cw)
  {
    ent_data (rent) = cw;
    ent_SetType (rent, MATRIX_DENSE_COMPLEX);
  }
  else
  {
    ent_data (rent) = mdr_Create(0,0);
    ent_SetType (rent, MATRIX_DENSE_REAL);
  }
  return (rent);
}


static void rec_map_func (MDR * x);

static int nx;
static MDR *xmdr=0;
static Ent *xent=0, *pent=0;
static Ent *fname=0;

//
// recurrence_map( mapfunc, funcparams, x0, n0 /,nth/ )
//
#undef THIS_SOLVER
#define THIS_SOLVER "rmap"
Ent *
ent_recurrence_map (int nargs, Datum args[])
{
  Ent *e2=0, *e3=0, *e4=0, *e5=0, *rent;
  MDR *w, *x, *x0;
  int n0, nt=0, i, j;

  if (nargs < 4)
  {
    printf ("rmap: Recurent map. Format:\n");
    printf ("rmap:   rmap(func,param,x0,N /,Nth/)\n");
    printf ("rmap:     func - function relating x_{i+1} and x_i,\n");
    printf ("rmap:     pars - paramater array for func,\n");
    printf ("rmap:     x0   - initial value,\n");
    printf ("rmap:     N    - number of iterations,\n");
    printf
      ("rmap:     Nth  - number of thermalization (disregarded) iterations,\n");
    rerror ("wrong arguments");
  }

  //
  // Get function ptr
  //
  fname = bltin_get_ent(args[0]);
  if (!isfuncent(fname))
    rerror (THIS_SOLVER ": " RLAB_ERROR_ARG1_FUNC_VAR "\n");

  //
  // parameter entity
  //
  e2 = bltin_get_ent (args[1]);
  if (ent_type (e2) != UNDEF)
    pent = ent_Copy (e2);

  // x0, an initial condition
  e3 = bltin_get_ent (args[2]);
  if (ent_type (e3) != MATRIX_DENSE_REAL)
    rerror ("rmap: Third argument 'x0' must be vector!");
  x0 = class_matrix_real (e3);
  if (!x0)
    rerror ("rmap: Third argument 'x0' must be vector!");
  nx = MNR (x0) * MNC (x0);

  // length of iteration
  e4 = bltin_get_ent (args[3]);
  if (ent_type (e4) != MATRIX_DENSE_REAL)
    rerror ("rmap: Fourth argument 'N0' must be scalar!");
  n0 = (int) class_double (e4);

  // length of thermalization
  if (nargs > 4)
  {
    e5 = bltin_get_ent (args[4]);
    if (ent_type (e5) == MATRIX_DENSE_REAL)
      nt = (int) class_double (e5);
  }
  w = mdr_Create (n0, nx);

  // Set up ENTITIES for user-function.
  // x
  xmdr = mdr_CreateEmpty (1, nx);
  xent = ent_Assign_Rlab_MDR (xmdr);
  ent_IncRef (xent);

  // params

  // do 'nt' thermalization (discarded) steps
  x = mdr_Float_BF (x0);

  for (i = 0; i < nt; i++)
    rec_map_func (x);

  // do 'n0' steps and record them
  for (i = 0; i < n0; i++)
  {
    rec_map_func (x);
    for (j = 0; j < nx; j++)
      Mdr0 (w, i, j) = MdrV0 (x, j);
  }

  // user function: clean parameter entity
  if (pent)
    ent_Clean (pent);
  pent = 0;

  // user function: clean x
  MDPTR(xmdr) = 0;
  ent_DecRef (xent);
  ent_Destroy (xent);

  // clean input
  ent_Clean(fname);
  ent_Clean (e2);
  ent_Clean (e3);
  ent_Clean (e4);
  ent_Clean (e5);

  rent = ent_Create ();
  ent_data (rent) = w;
  ent_SetType (rent, MATRIX_DENSE_REAL);
  return (rent);
}

//
// The interface to the user-specified function.
//
static void
rec_map_func (MDR * x)
{
  Ent *rent = 0;
  MDR *retm = 0;
  int i;

  // get x from the list of arguments
  MDPTR(xmdr) = MDPTR(x);

  if (pent)
    rent = ent_call_rlab_script_2args(fname, xent, pent);
  else
    rent = ent_call_rlab_script_1arg (fname, xent);

  if (ent_type(rent)!=MATRIX_DENSE_REAL)
    rerror (THIS_SOLVER ": " RLAB_ERROR_RHS_FUNC_MUST_RETURN_MDR);

  retm = ent_data(rent);

  if (SIZE(retm)!=nx)
    rerror (THIS_SOLVER ": " RLAB_ERROR_RHS_FUNC_INCORRECT_DIM);

  for (i=0; i<nx;i++)
    MdrV0(x,i) = MdrV0(retm,i);

  ent_Clean (rent);
  return;
}

//
// TISEAN: poincare
//
static MDR *
mdr_poincare (MDR * series, double delay, int dim, int comp, double where,
              int dir)
{
  MDR *w, *wret;
  unsigned long i, count = 1;
  int length, icomp;
  long j, jd;
  double delta, xcut;
  double time = 0.0, lasttime = 0.0;

  length = MNR (series);
  w = mdr_Create (length, dim);

  if (dir == 0)
  {
    for (i = (comp - 1) * delay; i < length - (dim - comp) * delay - 1; i++)
    {
      if ((Mdr0 (series, i, 0) < where) && (Mdr0 (series, i + 1, 0) >= where))
      {
        delta = Mdr0 (series, i, 0) - where;
        delta = delta / (Mdr0 (series, i, 0) - Mdr0 (series, i + 1, 0));
        time = (double) i + delta;
        if (lasttime > 0.0)
        {
          icomp = 0;
          for (j = -(comp - 1); j <= dim - comp; j++)
          {
            if (j != 0)
            {
              jd = i + j * delay;
              xcut = Mdr0 (series, jd, 0);
              xcut += delta * (Mdr0 (series, jd + 1, 0) - Mdr0 (series, jd, 0));
              icomp++;
              Mdr1 (w, count, icomp) = xcut;
            }
          }
          Mdr1 (w, count, dim) = time - lasttime;
          count++;
        }
        lasttime = time;
      }
    }
  }
  else
  {
    for (i = (comp - 1) * delay; i < length - (dim - comp) * delay - 1; i++)
    {
      if ((Mdr0 (series, i, 0) > where) && (Mdr0 (series, i + 1, 0) <= where))
      {
        delta = Mdr0 (series, i, 0) - where;
        delta = delta / (Mdr0 (series, i, 0) - Mdr0 (series, i + 1, 0));
        time = (double) i + delta;
        if (lasttime > 0.0)
        {
          icomp = 0;
          for (j = -(comp - 1); j <= dim - comp; j++)
          {
            if (j != 0)
            {
              jd = i + j * delay;
              xcut = Mdr0 (series, jd, 0);
              xcut += delta * (Mdr0 (series, jd + 1, 0) - Mdr0 (series, jd, 0));
              icomp++;
              Mdr1 (w, count, icomp) = xcut;
            }
          }
          Mdr1 (w, count, dim) = time - lasttime;
          count++;
        }
        lasttime = time;
      }
    }
  }

  // at this point w contains relevant poincare section up to the row
  // count-1. This is what is returned
  wret = mdr_Create (count - 1, dim);
  for (i = 1; i < count; i++)
    for (j = 1; j <= dim; j++)
      Mdr1 (wret, i, j) = Mdr1 (w, i, j);
    return wret;
}

Ent *
ent_poincare (int nargs, Datum args[])
{
  Ent *X=0, *ED=0, *D=0, *Q=0, *C=0, *A=0, *rent;
  MDR *x;
  int dim, comp, dir;
  double where, delay;
  if (nargs != 6)
  {
    printf ("poincare: Poincare section of an uniformly sampled time series\n");
    printf ("poincare: given as a column-vector. Format:\n");
    printf ("poincare:   poincare(x,m,d,q,c,a)\n");
    printf ("poincare:     x - 1d time series,\n");
    printf ("poincare:     m - embedding dimension,\n");
    printf ("poincare:     d - delay,\n");
    printf ("poincare:     q - component for the crossing,\n");
    printf
      ("poincare:     c - crossing direction (0: from below, 1: from above),\n");
    printf ("poincare:     a - position of the crossing,\n");
    rerror ("wrong arguments!");
  }
  // get the time series in 'x'
  X = bltin_get_ent (args[0]);
  if (ent_type (X) != MATRIX_DENSE_REAL)
    rerror ("poincare: 'x' is a real column-vector, the time series!");
  x = class_matrix_real (X);
  if (!x)
    rerror ("poincare: First argument 'x' must be a real vector!");
  if (x->ncol != 1 && x->nrow != 1)
    rerror ("poincare: First argument 'x' must be a real vector!");

  // get the embedded dimensions
  ED = bltin_get_ent (args[1]);
  if (ent_type (ED) != MATRIX_DENSE_REAL)
    rerror ("poincare: Second argument 'q' must be scalar!");
  dim = (int) class_double (ED);

  // get the delay 'd'
  D = bltin_get_ent (args[2]);
  if (ent_type (D) != MATRIX_DENSE_REAL)
    rerror ("poincare: Third argument 'd' must be scalar!");
  delay = class_double (D);

  // get the component of the crossing 'q'
  Q = bltin_get_ent (args[3]);
  if (ent_type (Q) != MATRIX_DENSE_REAL)
    rerror ("poincare: Fourth argument 'q' must be scalar!");
  comp = (int) class_double (Q);

  // get the direction of the crossing
  C = bltin_get_ent (args[4]);
  if (ent_type (C) != MATRIX_DENSE_REAL)
    rerror ("poincare: Fifth argument  'c' must be scalar!");
  dir = class_double (C);

  // get the crossing position 'a'
  A = bltin_get_ent (args[5]);
  if (ent_type (A) != MATRIX_DENSE_REAL)
    rerror ("poincare: Sixth argument 'a' must be scalar!");
  where = class_double (A);

  if (comp > dim)
  {
    printf ("poincare: 'q' must be smaller than 'm' !\n");
    rerror ("poincare: dimension-component mismatch !");
  }

  ent_Clean (X);
  ent_Clean (ED);
  ent_Clean (D);
  ent_Clean (Q);
  ent_Clean (C);
  ent_Clean (A);

  rent = ent_Create ();
  ent_data (rent) = mdr_poincare (x, delay, dim, comp, where, dir);
  ent_SetType (rent, MATRIX_DENSE_REAL);
  return (rent);
}
