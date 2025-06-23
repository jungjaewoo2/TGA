// deoptim.c
// This file is a part of RLaB2 Rel.3
//  Copyright (C) 2013-2014 Marijan Kostrun. project rlabplus.
//  Based on
/***************************************************************
**                                                            **
**        D I F F E R E N T I A L     E V O L U T I O N       **
**                                                            **
** Program: de.c                                              **
** Version: 3.6                                               **
**                                                            **
** Authors: Dr. Rainer Storn                                  **
**          c/o ICSI, 1947 Center Street, Suite 600           **
**          Berkeley, CA 94707                                **
**          Tel.:   510-642-4274 (extension 192)              **
**          Fax.:   510-643-7684                              **
**          E-mail: storn@icsi.berkeley.edu                   **
**          WWW: http://http.icsi.berkeley.edu/~storn/        **
**          on leave from                                     **
**          Siemens AG, ZFE T SN 2, Otto-Hahn Ring 6          **
**          D-81739 Muenchen, Germany                         **
**          Tel:    636-40502                                 **
**          Fax:    636-44577                                 **
**          E-mail: rainer.storn@zfe.siemens.de               **
**                                                            **
**          Kenneth Price                                     **
**          836 Owl Circle                                    **
**          Vacaville, CA 95687                               **
**          E-mail: kprice@solano.community.net               **
**                                                            **
** This program implements some variants of Differential      **
** Evolution (DE) as described in part in the techreport      **
** tr-95-012.ps of ICSI. You can get this report either via   **
** ftp.icsi.berkeley.edu/pub/techreports/1995/tr-95-012.ps.Z  **
** or via WWW: http://http.icsi.berkeley.edu/~storn/litera.html*
** A more extended version of tr-95-012.ps is submitted for   **
** publication in the Journal Evolutionary Computation.       **
**                                                            **
** You may use this program for any purpose, give it to any   **
** person or change it according to your needs as long as you **
** are referring to Rainer Storn and Ken Price as the origi-  **
** nators of the the DE idea.                                 **
** If you have questions concerning DE feel free to contact   **
** us. We also will be happy to know about your experiences   **
** with DE and your suggestions of improvement.               **
**                                                            **
***************************************************************/

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
#include "mds.h"
#include "list.h"
#include "btree.h"
#include "bltin.h"
#include "util.h"
#include "mathl.h"
#include "function.h"
#include "lp.h"
#include "ode.h"

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <time.h>

#include "rfileio.h"
#include "rlab_solver_parameters_names.h"

static MDR *z;
static Ent *zent, *pent;
static Ent *de_eval_fname;
static int de_solver(
  int *istrategy,
  int *D, double *tmp, double *best, double *bestit,
  int *NP,              // scalar, population size
  double *cost,         // dim NP, costs in population
  double *cmin,         // scalar, minimum
  int *genmax,          // scalar
  double *b,            // dim Dx2,  bounds or NULL
  double *c,            // dim DxNP, initial guess,   modified
  double *d,            // dim DxNP, next iteration,  modified
  double *F,            // scalar
  double *FDITHER,      // scalar, should be positive and much smaller than F
  double *CR,           // scalar
  double *eabs,         // scalar, convergence criterion, when cvar below eabs stop, or NULL
  double *xwidth,       // scalar, convergence criterion, spread of population, or NULL
  double *xwwgt,        // dim D,  relative weight for different parameter entries in xwidth, or NULL
  int    *mfc,          // scalar, convergence criterion, max number of failed attempts in a row, or NULL
  double *target,       // scalar, known target of optimization otherwise NULL
  double *terr,         // scalar, how close to target we can be to stop iterating
  int *nfeval,          // scalar, return value, who cares?
  FILE  *fpout_ptr,     // file pointer,
  int   *refresh,       // scalar, how often to print out the values to file pointer
  double urng(void),    // uniform[0,1> random number generator
  double dcost(int *vlen, double *v)  // function for evaluating cost
);

//
// call to the rlab function for evaluating the cost
//
static double
de_eval (int *idummy, double *vals)
{
  Ent *rent=0;
  double rval;

  // prepare for call to rlab function
  MDPTR(z) = (void *) vals;

  if (pent)
    rent = ent_call_rlab_script_2args(de_eval_fname, zent, pent);
  else
    rent = ent_call_rlab_script_1arg (de_eval_fname, zent);

  if (ent_type (rent) != MATRIX_DENSE_REAL)
    rerror (THIS_SOLVER ": " RLAB_ERROR_RHS_FUNC_MUST_RETURN_MDR);

  rval = class_double (rent);
  ent_Clean (rent);
  return rval;
}

//
// interface to differential evolution solver
//
#undef THIS_SOLVER
#define THIS_SOLVER "diffevol"
Ent * ent_gsl_diffevol(int nargs, Datum args[])
{
  Ent *e2=0, *e3=0, *e4=0;

  MDR *x=0, *c=0, *d=0, *weights=0;
  int D,i;
  int istr=1;       // strategy
  FILE *fout_ptr=0; // output for messages
  int irefresh=10;
  double F = 0.7;
  double FDITHER=0;
  double CR = 0.5;
  int    NP = 0;
  int    genmax = 1000;
  int    nfeval;
  MDR *dbest=0, *bounds=0, *dcost=0;

  char *fname = 0;;

  double ddummy;
  double *eabs=0;
  double  known_eabs;
  int    *mfc=0;
  int     known_mfc=0;
  double *target=0;
  double  known_target;
  double  terr=1e-2;
  double *xwidth=0;
  double  known_xwidth=0;
  double *b=0;
  double *wgt=0;
  double  cmin;

  ListNode *node;

  if (nargs != 3 && nargs != 4)
  {
    fprintf (stdout,
             THIS_SOLVER ": Optimization/minimization using differential evolution\n");
    fprintf (stdout,
             THIS_SOLVER ": Format:\n");
    fprintf (stdout,
             THIS_SOLVER ":   y = " THIS_SOLVER "(cost,/p/,X/,options/),\n");
    fprintf (stdout,
             THIS_SOLVER ": where\n");
    fprintf (stdout,
             THIS_SOLVER ":   cost=function(x0 /,p/), and X=[x0_1,x0_2..]\n");
    fprintf (stdout,
             THIS_SOLVER ": is matrix of stacked column vectors containing initial guesses for the optimal solution,\n");
    fprintf (stdout,
             THIS_SOLVER ": while\n");
    fprintf (stdout,
             THIS_SOLVER ":   options=<<"
             RLAB_NAME_DIFFEVOL_STRATEGY ";"
             RLAB_NAME_DIFFEVOL_F        ";"
             RLAB_NAME_DIFFEVOL_FDITHER  ";"
             RLAB_NAME_DIFFEVOL_CR       ";"
             RLAB_NAME_DIFFEVOL_GENMAX   ";"
             RLAB_NAME_DIFFEVOL_BOUNDS   ";"
             RLAB_NAME_DIFFEVOL_STDOUT   ";"
             RLAB_NAME_DIFFEVOL_IPRINT   ";"
             RLAB_NAME_DIFFEVOL_TARGET   ";"
             RLAB_NAME_DIFFEVOL_TERR     ";"
             RLAB_NAME_DIFFEVOL_MFC      ";"
             RLAB_NAME_DIFFEVOL_XWIDTH   ";"
             RLAB_NAME_DIFFEVOL_WIDTHWGT ";"
             RLAB_NAME_DIFFEVOL_EABS
             ">>.\n");
    fprintf (stdout,
             THIS_SOLVER ": See manual for details.\n");
    rerror ( THIS_SOLVER ": requires at least 3 arguments");
  }

  //
  // Get function ptr
  //
  de_eval_fname = bltin_get_ent(args[0]);
  if (!isfuncent(de_eval_fname))
    rerror (THIS_SOLVER ": " RLAB_ERROR_ARG1_FUNC_VAR "\n");

  //
  // parameter entity: if provided
  //
  pent = 0;
  e2 = bltin_get_ent (args[1]);
  if (ent_type (e2) != UNDEF)
    pent = ent_Copy (e2);

  //
  // x0: initial guess
  //
  e3 = bltin_get_ent (args[2]);
  if (ent_type (e3) == MATRIX_DENSE_REAL)
    x = class_matrix_real (e3);
  if (!x)
    rerror(THIS_SOLVER ": 3rd argument must be a real matrix of dim DxNP!");
  if ((x->nrow) == 1 && (x->ncol) == 1)
    rerror(THIS_SOLVER ": 3rd argument must be a real matrix of dim DxNP!");
  if ((x->ncol) < 7)
    rerror(THIS_SOLVER ": 3rd argument must be a real matrix of dim DxNP, with NP>6!");

  D  = x->nrow;
  NP = x->ncol;
  c = mdr_Float_BF(x);
  d = mdr_Float_BF(x);

  //
  // options for the solver
  //
  if (nargs > 3)
  {
    e4 = bltin_get_ent (args[3]);
    if (ent_type (e4) == BTREE)
    {
      // strategy
      node = btree_FindNode (ent_data (e4), RLAB_NAME_DIFFEVOL_STRATEGY);
      if (node != 0)
      {
        istr = (int) class_double (var_ent (node));
        if (istr < 1 && istr > 10)
          istr = 1;
      }
      // F
      node = btree_FindNode (ent_data (e4), RLAB_NAME_DIFFEVOL_F);
      if (node != 0)
      {
        ddummy = class_double ( var_ent (node) );
        if (ddummy > 0.0)
          F = ddummy;
      }
      // FDITHER
      node = btree_FindNode (ent_data (e4), RLAB_NAME_DIFFEVOL_FDITHER);
      if (node != 0)
      {
        ddummy = class_double ( var_ent (node) );
        if (ddummy > 0.0)
          FDITHER = ddummy;
      }
      // CR
      node = btree_FindNode (ent_data (e4), RLAB_NAME_DIFFEVOL_CR);
      if (node != 0)
      {
        ddummy = class_double ( var_ent (node) );
        if (ddummy >= 0.0 && ddummy <= 1.0)
          CR = ddummy;
      }
      // eabs
      node = btree_FindNode (ent_data (e4), RLAB_NAME_DIFFEVOL_EABS);
      if (node != 0)
      {
        ddummy = class_double ( var_ent (node) );
        if (ddummy >= 0.0)
        {
          known_eabs = ddummy;
          eabs = &known_eabs;
        }
      }
      // mfc
      node = btree_FindNode (ent_data (e4), RLAB_NAME_DIFFEVOL_MFC);
      if (node != 0)
      {
        known_mfc = (int) class_double (var_ent (node));
        if (known_mfc > 0)
          mfc = &known_mfc;
      }
      // target
      node = btree_FindNode (ent_data (e4), RLAB_NAME_DIFFEVOL_TARGET);
      if (node != 0)
      {
        known_target = class_double ( var_ent (node) );
        target = &known_target;
      }
      // width
      node = btree_FindNode (ent_data (e4), RLAB_NAME_DIFFEVOL_XWIDTH);
      if (node != 0)
      {
        known_xwidth = class_double ( var_ent (node) );
        if (known_xwidth > 0)
          xwidth = &known_xwidth;
      }
      // weights of different entries in calculation of width
      if (xwidth)
      {
        node = btree_FindNode (ent_data (e4), RLAB_NAME_DIFFEVOL_WIDTHWGT);
        if (node != 0)
        {
          if (ent_type(var_ent (node)) == MATRIX_DENSE_REAL)
          {
            weights  = ent_data (var_ent (node));
            if ((weights->nrow) * (weights->ncol) != D)
              rerror (THIS_SOLVER ": '" RLAB_NAME_DIFFEVOL_WIDTHWGT "' entry should be vector of dim D");
            if (weights->type == RLAB_TYPE_DOUBLE)
              wgt = MDRPTR(weights);
          }
        }
      }
      // near target
      if (target)
      {
        node = btree_FindNode (ent_data (e4), RLAB_NAME_DIFFEVOL_TERR);
        if (node != 0)
        {
          terr = class_double ( var_ent (node) );
        }
      }
      // genmax
      node = btree_FindNode (ent_data (e4), RLAB_NAME_DIFFEVOL_GENMAX);
      if (node != 0)
      {
        genmax = (int) class_double (var_ent (node));
        if (genmax < D)
          genmax = 1000;
      }
      //
      // bounds
      //
      node = btree_FindNode (ent_data (e4), RLAB_NAME_DIFFEVOL_BOUNDS);
      if (node != 0)
      {
        if (ent_type(var_ent (node)) == MATRIX_DENSE_REAL)
        {
          bounds  = ent_data (var_ent (node));
          if ((bounds->nrow) != D || (bounds->ncol) != 2)
            rerror (THIS_SOLVER ": '"RLAB_NAME_DIFFEVOL_BOUNDS
                "' entry should be matrix [min_i, MAX_i], for i=1,..D\n");
        }
      }
      // standard output
      node = btree_FindNode (ent_data (e4), RLAB_NAME_DIFFEVOL_STDOUT);
      if (node != 0)
      {
        if (ent_type(var_ent (node))== MATRIX_DENSE_STRING)
        {
          fname = class_char_pointer(var_ent (node));
        }
      }
      // refresh
      if (fname)
      {
        node = btree_FindNode (ent_data (e4), RLAB_NAME_DIFFEVOL_IPRINT);
        if (node != 0)
        {
          irefresh = (int) class_double (var_ent (node));
          if (irefresh < 2)
            irefresh = 10;
        }
      }
    }
  }

  if (bounds)
  {
    b = GC_malloc(2 * D * sizeof(double));
    if (bounds->type == RLAB_TYPE_INT32)
      for (i=0; i<D; i++)
      {
        b[i]   = Mdi0(bounds,i,0);
        b[i+D] = Mdi0(bounds,i,1);
      }
    else
      for (i=0; i<D; i++)
      {
        b[i]   = Mdr0(bounds,i,0);
        b[i+D] = Mdr0(bounds,i,1);
      }
  }
  else
    b = 0;

  // Set up ENTITIES for user-function.
  //
  z = mdr_CreateEmpty (1,D);
  zent = ent_Assign_Rlab_MDR(z);
  ent_IncRef (zent);

  //
  // create all the memory
  //
  double *dtmp    = GC_malloc(D  * sizeof(double));
  dbest   = mdr_Create(D,1);
  double *dbestit = GC_malloc(D  * sizeof(double));
  dcost   = mdr_Create(1, NP);

  if (fname)
    fout_ptr = fopen (fname, "a");

  //
  // call the solver now
  //
  de_solver(&istr, &D, dtmp, MDRPTR(dbest), dbestit, &NP, MDRPTR(dcost), &cmin, &genmax,
             b, MDRPTR(c), MDRPTR(d), &F, &FDITHER, &CR, eabs, xwidth, wgt, mfc, target,
            &terr, &nfeval, fout_ptr, &irefresh,
             gslrnguf_, de_eval);

  if (fout_ptr)
    fclose (fout_ptr);

  //
  // clean-up the mess
  //
  mdr_Destroy(d);
  if (b)
    GC_free(b);
  GC_free (dtmp);
  GC_free (dbestit);

  // clean zent, but first detach MDPTR(z)
  MDPTR(z) = 0;
  ent_DecRef (zent);
  ent_Destroy (zent);

  ent_Clean (pent);
  ent_Clean (e2);
  ent_Clean (e3);
  ent_Clean (e4);
  ent_Clean (de_eval_fname);

  // create output list
  Btree * bw = btree_Create ();
  install (bw, "coef", ent_Assign_Rlab_MDR(dbest));
  install (bw, "coef_cost", ent_Create_Rlab_Double(cmin));
  install (bw, "pop", ent_Assign_Rlab_MDR(c));
  install (bw, "pop_cost", ent_Assign_Rlab_MDR(dcost));

  return ent_Assign_Rlab_BTREE(bw);
}


//
// Copyright (C) 2013 Marijan Kostrun. project rlabplus:
// sampling without replacement, where the drawing order is important
// input:
//    N - the length of array from which we sample
//    k - length of the sample
//    random() - function that returns uniform random variate on [0,1>
// input/output:
//    r[k+1] - input r[0] first index that is known, see de_solver
//           - output r[1:k] indices from {0 .. (N-1)} \ {r[0]} without replacement
//    s[k+1] - input s[0]=r[0]
//           - output s[1:k] sorted indices from r[1:k]
// all variables are passed through their pointers, so fortran can use this as well.
static int de_sample(int * N, int * k, int * r, int * s, double random(void) )
{
  int l,j,u,c;

  for (j=1; j <= (*k); j++)
  {
    u = ((*N) - j) * random();
    r[j] = u;
    for (l=0; l<j; l++)
      r[j] += (r[j] >= s[l]);
    s[j] = r[j];

    // now sort s
    for (l=j-1; l>=0; l--)
    {
      if (s[l+1] < s[l])
      {
        c = s[l+1];
        s[l+1] =s[l];
        s[l] = c;
      }
      else
        break;
    }
  }

  return 0;
}


static char  *strat[] =                           // strategy names
{
  "",                                             // 0
  "DE/best/1/exp",                                // 1
  "DE/rand/1/exp",                                // 2
  "DE/rand-to-best/1/exp",                        // 3
  "DE/best/2/exp",                                // 4
  "DE/rand/2/exp",                                // 5
  "DE/best/1/bin",                                // 6
  "DE/rand/1/bin",                                // 7
  "DE/rand/1/bin/either-or-algorithm",            // 8
  "DE/rand-to-best/1/bin",                        // 9
  "DE/best/2/bin",                                // 10
  "DE/rand/2/bin"                                 // 11
};


static int de_solver(
  int *istrategy,
  int *D, double *tmp, double *best, double *bestit,
  int *NP,
  double *cost,       // dim NP, costs in population
  double *cmin,       // scalar, minimum
  int *genmax,        // scalar
  double *b,          // dim Dx2,  bounds or NULL
  double *c,          // dim DxNP, initial guess,   modified
  double *d,          // dim DxNP, next iteration,  modified
  double *F,          // scalar
  double *FDITHER,    // scalar, should be positive and much smaller than F
  double *CR,         // scalar
  double *eabs,       // scalar, convergence criterion, spread in cost of population
  double *xwidth,     // scalar, convergence criterion, spread of population
  double *xwwgt,      // dim D, relative weights for different components in xwidth
  int    *mfc,        // scalar, convergence criterion, max number of failed attempts in a row
  double *target,     // scalar, known target of optimization otherwise NULL
  double *terr,       // scalar, how close to target to call it quits
  int *nfeval,        // scalar, return value
  FILE  *fpout_ptr,   // file pointer,
  int   *refresh,     // scalar, how often to print out the values to file pointer
  double urng(void),  // uniform[0,1> random number generator
  double cost_fun(int *vlen, double *v)  // function for evaluating cost
)
{
  int     i, j, L, n=0;      /* counting variables                 */
  int     num_rs=5;           // How many random columns we generate depends on strategy. Max is five.
  int     r[12];              /* placeholders for random indices    */
  int     imin;               /* index to member with lowest energy */
  int     gen;
  double  trial_cost;        /* buffer variable                    */
  double  cvar=0;              /* computes the cost variance         */
  double  cmean;             /* mean cost                          */
//   double  cost_delta;
  double  s, s2, v2=0;
  double *pold, *pnew, *pswap;
  int     fail_count=0;

  pold = c; /* old population (generation G)   */
  pnew = d; /* new population (generation G+1) */

  // evaluate cost of all samples in original population
  (*nfeval) = 0;
  for (i=0; i<(*NP); i++)
  {
    (*nfeval)++;
    cost[i] = cost_fun(D, &pold[i*(*D)]); /* obj. funct. value */
  }
  // find the best cost value in the initial population
  (*cmin)=cost[0];
  imin=0;
  for (i=1; i<(*NP); i++)
  {
    if (cost[i]<(*cmin))
    {
      (*cmin) = cost[i];
      imin = i;
    }
  }
  // find the best vector in the initial population
  memcpy(best,  &pold[imin*(*D)], (*D)*sizeof(double));
  memcpy(bestit,&pold[imin*(*D)], (*D)*sizeof(double));

  // based on the strategy adjust how many indices are needed
  switch ((*istrategy))
  {

    case 1:
      //  DE/best/1/exp
      //  Our oldest strategy but still not bad. However, we have found several
      //  optimization problems where misconvergence occurs
      num_rs = 2; // adjust how many r's are generated
      break;

    case 2:
      //  DE/rand/1/exp
      //  This is one of my favourite strategies. It works especially well when the
      //  "bestit[]"-schemes experience misconvergence. Try e.g. F=0.7 and CR=0.5
      //  as a first guess
      num_rs = 3; // adjust how many r's are generated
      break;

    case 3:
      //  DE/rand-to-best/1/exp
      //  This strategy seems to be one of the best strategies. Try F=0.85 and CR=1
      //  If you get misconvergence try to increase NP. If this doesn't help you
      //  should play around with all three control variables
      num_rs = 2; // adjust how many r's are generated
      break;

    case 4:
      //  DE/rand/2/exp seems to be a robust optimizer for many functions
      num_rs = 4; // adjust how many r's are generated
      break;

    case 5:
      num_rs = 5; // adjust how many r's are generated
      break;

    case 6:
      //  DE/best/1/bin
      num_rs = 2; // adjust how many r's are generated
      break;

    case 8:
    case 7:
      //  DE/rand/1/bin
      num_rs = 3; // adjust how many r's are generated
      break;

    case 9:
      //  DE/rand-to-best/1/bin
      num_rs = 2; // adjust how many r's are generated
      break;

    case 10:
      //  DE/best/2/bin
      num_rs = 4; // adjust how many r's are generated
      break;

    default:
      //  DE/rand/2/bin
      num_rs = 5; // adjust how many r's are generated
      break;
  }


/*=======================================================================*/
/*=========Iteration loop================================================*/
/*=======================================================================*/
  gen = 0;
  while (gen < (*genmax))
  {
    gen++;
    imin = 0;

    for (i=0; i<(*NP); i++)         /* Start of loop through population  */
    {
      // draw a random sample from population without replacement, order is important
      //  r[0] = i
      //  r[1], r[2], r[3], r[4], r[5], r[5] \in {0,..(*NP)-1} \ {i}
      //  r[6:11] contains 'r[0:5]' sorted in ascending order
      r[6] = r[0] = i;
      de_sample(NP, &num_rs, &r[0], &r[6], urng);

/*    for (n=0;n<6;n++)
        fprintf(stderr, " %i", r[n]);
      fprintf(stderr,"\n");*/

/*=======Choice of strategy===============================================================*/
/*=======We have tried to come up with a sensible naming-convention: DE/x/y/z=============*/
/*=======DE :  stands for Differential Evolution==========================================*/
/*=======x  :  a string which denotes the vector to be perturbed==========================*/
/*=======y  :  number of difference vectors taken for perturbation of x===================*/
/*=======z  :  crossover method (exp = exponential, bin = binomial)=======================*/
/*                                                                                        */
/*=======There are some simple rules which are worth following:===========================*/
/*=======1)  F is usually between 0.5 and 1 (in rare cases > 1)===========================*/
/*=======2)  CR is between 0 and 1 with 0., 0.3, 0.7 and 1. being worth to be tried first=*/
/*=======3)  To start off NP = 10*D is a reasonable choice. Increase NP if misconvergence=*/
/*           happens.                                                                     */
/*=======4)  If you increase NP, F usually has to be decreased============================*/
/*=======5)  When the DE/best... schemes fail DE/rand... usually works and vice versa=====*/

      memcpy(tmp,  &pold[i*(*D)], (*D)*sizeof(double));
      n = (int)( urng() * (*D) );
      L = 0;
      switch ((*istrategy))
      {

        case 1:
          //  DE/best/1/exp
          //  Our oldest strategy but still not bad. However, we have found several
          //  optimization problems where misconvergence occurs
          do
          {
            tmp[n] = bestit[n] + ((*F) + (*FDITHER) * urng()) * ( pold[ r[1] * (*D) + n ] - pold[ r[2] * (*D) + n ] );
            n = (n+1) % (*D);
            L++;
          }
          while( (urng() < (*CR) ) && (L < (*D)) );
          break;

        case 2:
          //  DE/rand/1/exp
          //  This is one of my favourite strategies. It works especially well when the
          //  "bestit[]"-schemes experience misconvergence. Try e.g. F=0.7 and CR=0.5
          //  as a first guess
          do
          {
            tmp[n] = pold[ r[1] * (*D) + n ]
            + ((*F) + (*FDITHER) * urng()) * ( pold[ r[2] * (*D) + n ]- pold[ r[3] * (*D) + n ] );
            n = (n+1) % (*D);
            L++;
          }
          while( (urng() < (*CR) ) && (L < (*D)) );
          break;

        case 3:
          //  DE/rand-to-best/1/exp
          //  This strategy seems to be one of the best strategies. Try F=0.85 and CR=1
          //  If you get misconvergence try to increase NP. If this doesn't help you
          //  should play around with all three control variables
          do
          {
            tmp[n] = tmp[n] + ((*F) + (*FDITHER) * urng()) * (bestit[n] - tmp[n])
                + (*F) * ( pold[ r[1] * (*D) + n ] - pold[ r[2] * (*D) + n ] );
            n = (n+1) % (*D);
            L++;
          }
          while( (urng() < (*CR) ) && (L < (*D)) );
          break;

        case 4:
          //  DE/rand/2/exp seems to be a robust optimizer for many functions
          do
          {
            tmp[n] = bestit[n]
                + ((*F) + (*FDITHER) * urng()) * ( pold[r[1] * (*D) + n] + pold[r[2] * (*D) + n]
                    - pold[r[3] * (*D) + n] - pold[r[4] * (*D) + n] );
            n = (n+1) % (*D);
            L++;
          }
          while( (urng() < (*CR) ) && (L < (*D)) );
          break;

        case 5:
          do
          {
            tmp[n] = pold[r[5] * (*D) + n]
                + ((*F) + (*FDITHER) * urng()) * ( pold[r[1] * (*D) + n] + pold[r[2] * (*D) + n]
                  - pold[r[3] * (*D) + n] - pold[r[4] * (*D) + n] );
            n = (n+1) % (*D);
            L++;
          }
          while( (urng() < (*CR)) && (L < (*D)) );
          break;

        case 6:
          //  DE/best/1/bin
          for (L=0; L<(*D); L++)
          {
            // perform D binomial trials
            if ( (urng() < (*CR)) || (L == ((*D)-1)) )
            {
              /* change at least one parameter */
              tmp[n] = bestit[n] + ((*F) + (*FDITHER) * urng()) * ( pold[r[1] * (*D) + n] - pold[r[2] * (*D) +n] );
            }
            n = (n+1) % (*D);
          }
          break;

        case 8:
          //  DE/rand/1/bin : either-or-algorithm
          if (urng() < 0.5)
          {
            for (L=0; L<(*D); L++)
            {
              // perform D binomial trials
              if ( (urng() < (*CR)) || (L == ((*D)-1)) )
              {
                //  change at least one parameter
                tmp[n] = pold[r[1] * (*D) + n] + 0.5 * (((*F) + (*FDITHER) * urng()) +  1.0) *
                  ( pold[r[2] * (*D) + n] + pold[r[3] * (*D) + n]
                      - 2 * pold[r[1] * (*D) + n] );
              }
              n = (n+1) % (*D);
            }
            break;
          }
        case 7:
          //  DE/rand/1/bin
          for (L=0; L<(*D); L++)
          {
            // perform D binomial trials
            if ( (urng() < (*CR)) || (L == ((*D)-1)) )
            {
              //  change at least one parameter
              tmp[n] = pold[r[1] * (*D) + n] + ((*F) + (*FDITHER) * urng()) * ( pold[r[2] * (*D) + n]
                  - pold[r[3] * (*D) + n] );
            }
            n = (n+1) % (*D);
          }
          break;

        case 9:
          //  DE/rand-to-best/1/bin
          for (L=0; L<(*D); L++)
          {
            // perform D binomial trials
            if ( (urng() < (*CR)) || (L == ((*D)-1)) )
            {
              tmp[n] = tmp[n] + ((*F) + (*FDITHER) * urng()) * ( bestit[n] - tmp[n] )
                  + ((*F) + (*FDITHER) * urng()) * ( pold[r[1] * (*D) + n] - pold[r[2] * (*D) + n] );
            }
            n = (n+1) % (*D);
          }
          break;

        case 10:
          //  DE/best/2/bin
          for (L=0; L<(*D); L++)
          {
            // perform D binomial trials
            if ( (urng() < (*CR)) || (L == ((*D)-1)) )
            {
              tmp[n] = bestit[n] +
                  ((*F) + (*FDITHER) * urng()) * ( pold[r[1] * (*D) + n] + pold[r[2] * (*D) + n]
                      - pold[r[3] * (*D) + n] - pold[r[4] * (*D) + n] );
            }
            n = (n+1) % (*D);
          }
          break;

        default:
          //  DE/rand/2/bin
          for (L=0; L<(*D); L++)
          {
            // perform D binomial trials
            if ( (urng() < (*CR)) || (L == ((*D)-1)) )
            {
              tmp[n] = pold[r[5] * (*D) + n]  +
                ((*F) + (*FDITHER) * urng()) * ( pold[r[1] * (*D) + n] + pold[r[2] * (*D) + n]
                    - pold[r[3] * (*D) + n] - pold[r[4] * (*D) + n] );
            }
            n = (n+1) % (*D);
          }
          break;
      }
      // adjust for bounds
      if (b)
      {
        for (L=0; L<(*D); L++)
        {
          if (tmp[L] < b[L])
            tmp[L] = b[L] + urng() * (pold[i * (*D) + L] - b[L]);
          else if (tmp[L] > b[(*D) + L])
            tmp[L] = b[(*D) + L] + urng() * (pold[i * (*D) + L] - b[(*D) + L]);
        }
      }

      //  Trial mutation now in tmp[]. Test how good this choice really was
      (*nfeval)++;
      trial_cost = cost_fun(D, tmp);
      if (trial_cost <= cost[i])   /* improved objective function value ? */
      {
//         cost_delta = cost[i] - trial_cost;  // new improvement
        cost[i] = trial_cost;
        memcpy(&pnew[i*(*D)], tmp, (*D)*sizeof(double));
        if (trial_cost<(*cmin))
        {
          (*cmin) = trial_cost;
          imin = i;
          memcpy(best, tmp, (*D)*sizeof(double));
        }
        fail_count = 0;
      }
      else
      {
        memcpy(&pnew[i*(*D)], &pold[i*(*D)], (*D)*sizeof(double));
        fail_count++;
      }

      // test whether we have reached the target within terr
      if (target)
      {
        if ( ABS(trial_cost-(*target)) < (*terr) )
          return 0;
      }
      // did we specify how many failures we would tolerate
      if (mfc)
      {
        if (fail_count > (*mfc))
          return 0;
      }

    } //for (i=0; i<(*NP); i++)

    // store the current best guess per one iteration over generation
    memcpy(bestit, best, (*D)*sizeof(double));

    // swap population arrays. New generation becomes old one
    pswap = pold;
    pold  = pnew;
    pnew  = pswap;

    // Compute the energy variance: its width can be used as a stopping criteron
    if (eabs || fpout_ptr)
    {
      cmean = 0.0;
      cvar  = 0.0;
      for (j=0; j<(*NP); j++)
      {
        cmean += cost[j];
        cvar  += cost[j] * cost[j];
      }
      cmean = cmean / (*NP);
      cvar  = sqrt( ABS(cvar / (*NP) - cmean * cmean) * (*NP) / ((*NP) - 1) );
    }

    // compute the spread of population in parameter space
    if (xwidth)
    {
      v2 = 0;
      for (L=0; L<(*D); L++)
      {
        s = 0;
        s2 = 0;
        for (i=0; i<(*NP); i++)
        {
          s  += pold[i*(*D) + L];
          s2 += pold[i*(*D) + L] * pold[i*(*D) + L];
        }
        if (xwwgt)
          v2 += xwwgt[L] * ABS(s2 - s * s / (*NP)) / (*NP);
        else
          v2 += ABS(s2 - s * s / (*NP)) / (*NP);
      }
      v2 = sqrt(v2);
    }

    if (fpout_ptr)
    {
      //  messaging: print progress report
      if ( gen % (*refresh) ==1 )
      {
        fprintf(fpout_ptr, THIS_SOLVER ": min cost = %-15.10g",(*cmin));
        for (j=0;j<(*D);j++)
        {
          fprintf(fpout_ptr, "best[%d]=%-15.10g", j+1, best[j]);
        }
        if (xwidth)
          fprintf(fpout_ptr, "width=%g", v2);
        fprintf(fpout_ptr, "\n" THIS_SOLVER ": Generation=%d, NFEs=%i, Strategy=%s, ",
                gen, *nfeval, strat[(*istrategy)] );
        fprintf(fpout_ptr, "NP=%d, F=%-4.2g", (*NP), (*F));
        if ((*FDITHER))
          fprintf(fpout_ptr, "(FDITHER=%g)", (*FDITHER) );
        fprintf(fpout_ptr, ", CR=%-4.2g, cost-variance=%-10.5g\n", (*CR), cvar);
      }
    } // if (fpout_ptr)

    // what is the spread of the population
    if (xwidth)
      if (v2 < (*xwidth))
        return 0;

    // check whether the new improvement over cost is worth while
    if (eabs)
    {
      if (cvar < (*eabs))
        return 0;
    }

  } // while (gen < (*genmax)) // next generation

  return 0;

} // while ((gen < (*genmax))


