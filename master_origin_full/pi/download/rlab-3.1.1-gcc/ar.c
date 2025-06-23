//
// ar.c
// rlab <-> arpack interface
//

/*  This file is a part of rlabplus
   Copyright (C) 2005 M. Kostrun

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

#include "rlab.h"
#include "mem.h"
#include "ent.h"
#include "btree.h"
#include "symbol.h"
#include "mdr.h"
#include "mdrf1.h"
#include "mdrf2.h"
#include "mdc.h"
#include "mdcf1.h"
#include "mdcf2.h"
#include "msr.h"
#include "msrf1.h"
#include "msrf2.h"
#include "msc.h"
#include "mscf1.h"
#include "mscf2.h"
#include "util.h"
#include "mathl.h"

#include "ar.h"
#include "fi.h"
#include "lp.h"
// #include "blas.h"

#ifdef HAVE_ARPACK

#include "arpack.h"

#include <stdio.h>
#include <math.h>

#include "./clibs/superlu/src/dcomplex.h"
#include "./clibs/superlu/src/dsp_defs.h"

//
// from msrf2.c: need for spsolve
//
extern SuperMatrix *superlu_L;
extern SuperMatrix *superlu_U;
extern int *superlu_perm_r;
extern int *superlu_perm_c;

extern int umfpack_n;
extern void *umfpack_Numeric;
extern MSR *umfpack_a;

// arpack parameters
int ar_maxiter = 3000;


void
mdc_ar_eigs_ge2 (MDC * A, MDC * B, void **val, void **vec, int nev, char * which, Complex sigma)
{
  int info=0, ierr=0, lworkl, n, i, j, ncv;
  int *iselect, iparam[11], ipntr[14], *iwrkd=0;
  MDC *workd, *d, *workl, *v, *resid, *workev;
  MDC *w, *v2, *A2;
  MDR *rwork, *rd;

  mdc_Detect_Inf (A);
  mdc_Detect_Nan (A);
  mdc_Detect_Inf (B);
  mdc_Detect_Nan (B);
  //
  // Some rudimentary checks
  //
  n = MNC (A);

  if (MNR (A) != MNC (A)) rerror ("eigs: A must be a square matrix!");
  if (MNR (B) != MNC (B)) rerror ("eigs: B must be a square matrix!");
  if (MNR (A) != MNC (B)) rerror ("eigs: A and B must be the matrices of the same size!");

  ncv = MIN(6*nev,n);

  iselect = (int *)GC_MALLOC(ncv*sizeof(int));

  lworkl = 3*ncv*ncv+5*ncv;
  workl  = mdc_Create(1, lworkl);
  d      = mdc_Create(1,ncv);
  resid  = mdc_Create(1,n);
  v      = mdc_Create(n,ncv);
  workd  = mdc_Create(1,3*n);
  workev = mdc_Create(1,2*ncv);

  rwork  = mdr_Create(1,n);
  rd     = mdr_Create(ncv,3);

  iparam[1 -1] = 1;             // shift
  iparam[3 -1] = ar_maxiter;    // maximum number of iterations
  iparam[7 -1] = 3;             // use mode 3

  //   isigma==0,rsigma!=0; use (A-rsigma*I)^{-1}
  A2 = mdc_Copy( A );
  iwrkd  = (int *)GC_MALLOC(n*sizeof(int));
  ZNDRV2(MDCPTR(v), MDCPTR(workl), &lworkl, MDCPTR(workev),
         MDCPTR(workd), MDCPTR(d), MDCPTR(resid), iselect, iparam, ipntr,
         MDCPTR(A2), MDCPTR(B), &n, &nev, &ncv,
         which, &sigma, &info, &ierr, iwrkd,
         MDRPTR(rwork), MDRPTR(rd));
  if ( info < 0 )
  {
    fprintf(stdout, "Internal error in ARPACK 'znaupd', info = %i\n", info);
    iparam[4] = 0;
  }
  if ( ierr < 0 )
  {
    fprintf(stdout, "Internal error in ARPACK 'zneupd', ierr = %i\n", ierr);
    iparam[4] = 0;
  }

  if ( info == 1 )
    fprintf(stdout, "'zsaupd' reached the maximum number of iterations.");

  if ( info == 3 )
  {
    fprintf(stdout, "No shifts could be applied during implicit Arnoldi update,\n");
    fprintf(stdout, "try increasing 'NCV'.\n");
  }

  ncv = iparam[4];    // number of converged eigenvalues

  if(ncv>0)
  {
    w = mdc_Create(1, ncv);
    for(i=0;i<ncv;i++)
    {
      MdcV0(w,i) = MdcV0(d,i);
    }
    v2 = mdc_Create(n,ncv);
    for(j=0;j<ncv;j++)
    {
      Complex dummy;
      dummy =  conj(Mdc0(v,0,j));
      if( cabs(dummy) > 0)
      {
        //
        // normalize eigenvectors so that first coefficient is real (LAPACK)
        //
        Mdc0(v2,0,j) = cabs(dummy);
        for(i=1;i<n;i++)
        {
          Mdc0(v2,i,j) = Mdc0(v,i,j)  * dummy / Mdc0(v2,0,j);
        }
      }
      else
      {
        //
        // don't do anything. First coefficient of eigenvector is zero
        //
        for(i=0;i<n;i++)
        {
          Mdc0(v2,i,j)  = Mdc0(v,i,j);
        }
      }
    }
  }
  else
  {
    w  = mdc_Create(0,0);
    v2 = mdc_Create(0,0);
  }

  *val = w;
  *vec = v2;

  mdc_Destroy( A2 );
  GC_FREE(iwrkd);

  mdc_Destroy (workl); mdc_Destroy (d); mdc_Destroy (resid);
  mdc_Destroy (workd); mdc_Destroy (workev); mdc_Destroy (v);
  mdr_Destroy(rd); mdr_Destroy(rwork);

  GC_FREE(iselect);

  return;
}


Btree * mdc_ar_EigsG (MDC * a, MDC * b, int nev, char *which, Complex sigma)
{
  Btree *bt;
  Ent *eval, *evec;
  void *val, *vec;

  //
  // Create the list that will contain the results.
  //
  bt = btree_Create ();

  //
  // Check for symmetry of [A]
  // Then use the correct routine.
  //
  mdc_ar_eigs_ge2 (a, b, &val, &vec, nev, which, sigma);

  //
  // Hook the results into the list.
  //
  eval = ent_Create ();
  ent_SetType (eval, MATRIX_DENSE_COMPLEX);
  ent_data (eval) = val;
  install (bt,  "val", eval);
  evec = ent_Create ();
  ent_SetType (evec, MATRIX_DENSE_COMPLEX);
  ent_data (evec) = vec;
  install (bt,  "vec", evec);

//   ent_Clean( eval );
//   ent_Clean( evec );

  return (bt);
}

void
mdr_ar_eigs_ge2 (MDR * A, MDR *B, void **val, void **vec, int nev, char * which,
                    double sigma)
{
  int info=0, ierr=0, ldv, lworkl, n, i, j, ncv;
  int *iselect, iparam[11], ipntr[11], *iwrkd=0;
  MDR *workd, *d, *workl, *v, *resid, *workev;
  MDR *w, *v2;

  mdr_Detect_Inf (A);
  mdr_Detect_Nan (A);

  mdr_Detect_Inf (B);
  mdr_Detect_Nan (B);

  if (B->type == RLAB_TYPE_INT32)
    rerror ("eigs: B must be a real/complex matrix!");

  //
  // Some rudimentary checks
  //

  n = MNC (A);

  if (MNR (A) != MNC (A)) rerror ("eigs: A must be a square matrix!");
  if (MNR (B) != MNC (B)) rerror ("eigs: B must be a square matrix!");
  if (MNR (A) != MNC (B)) rerror ("eigs: A and B must be the matrices of the same size!");

  ldv = n;
  ncv = MIN(6*nev,n);

  iselect = (int *)GC_MALLOC(ncv*sizeof(int));

  lworkl = 3*ncv*ncv+6*ncv;
  workl  = mdr_Create(1, lworkl);
  d      = mdr_Create(ncv,3);
  resid  = mdr_Create(1,n);
  v      = mdr_Create(ldv,ncv);
  workd  = mdr_Create(1,3*n);
  workev = mdr_Create(1,3*ncv);

  iparam[1 -1] = 1;            // shift
  iparam[3 -1] = ar_maxiter;   // maximum number of iterations
  iparam[7 -1] = 3;

  //
  MDR * A2 = mdr_Float_BF (A);
  iwrkd  = (int *)GC_MALLOC(n*sizeof(int));
  DNDRV2(MDRPTR(v), &ldv, MDRPTR(workl), &lworkl, MDRPTR(workev),
         MDRPTR(workd), MDRPTR(d), MDRPTR(resid), iselect, iparam, ipntr,
         MDRPTR(A2), MDRPTR(B), &n, &nev, &ncv, which, &sigma, &info, &ierr,
         iwrkd);
  mdr_Destroy( A2 );
  GC_FREE(iwrkd);

  if ( info < 0 )
  {
    fprintf(stdout, "Internal error in ARPACK 'dnaupd', info = %i\n", info);
    iparam[4] = 0;
  }
  if ( ierr < 0 )
  {
    fprintf(stdout, "Internal error in ARPACK 'dneupd', ierr = %i\n", ierr);
    iparam[4] = 0;
  }
  if ( info == 1 )
    fprintf(stdout, "'dsaupd' reached the maximum number of iterations.");
  if ( info == 3 )
  {
    fprintf(stdout, "No shifts could be applied during implicit Arnoldi update,\n");
    fprintf(stdout, "try increasing 'NCV'.\n");
  }

  ncv = iparam[4];    // number of converged eigenvalues

  if(ncv>0)
  {
    w = mdr_Create(1, ncv);
    for(i=0;i<ncv;i++) MdrV0(w,i) = MdrV0(d,i);
    {
      // fix the eigenvalues
      //MdrV0(w,i) = MdrV0(d,i);
      //else MdrV0(w,i) = rsigma + 1.0/MdrV0(d,i);
    }
    v2 = mdr_Create(n,ncv);
    for(j=0;j<ncv;j++) for(i=0;i<n;i++)
    {
      Mdr0(v2,i,j) = Mdr0(v,i,j);
    }
  }
  else
  {
    w  = mdr_Create(0,0);
    v2 = mdr_Create(0,0);
  }

  *val = w;
  *vec = v2;

  mdr_Destroy (workl); mdr_Destroy (d); mdr_Destroy (resid);
  mdr_Destroy (workd); mdr_Destroy (workev); mdr_Destroy (v);
  GC_FREE(iselect);
  return;
}

Btree * mdr_ar_EigsG (MDR * a, MDR * b, int nev, char *which, double sigma)
{
  Btree *bt;
  Ent *eval, *evec;
  void *val, *vec;

  //
  // Create the list that will contain the results.
  //
  bt = btree_Create ();

  //
  mdr_ar_eigs_ge2 (a, b, &val, &vec, nev, which, sigma);

  //
  // Hook the results into the list.
  //
  eval = ent_Create ();
  ent_SetType (eval, MATRIX_DENSE_REAL);
  ent_data (eval) = val;
  install (bt,  "val", eval);
  evec = ent_Create ();
  ent_SetType (evec, MATRIX_DENSE_REAL);
  ent_data (evec) = vec;
  install (bt,  "vec", evec);

//   ent_Clean( eval );
//   ent_Clean( evec );
  return (bt);
}

void mdr_ar_eigs_sym
    (
    MDR * M, void **val, void **vec, int nev, char *which,
    double rsigma
    );
void mdr_ar_eigs_ge
    (
    MDR * M, void **val, void **vec, int nev, char * which,
    double rsigma
    );

Btree * mdr_ar_EigsS (MDR * a, int nev, char *which, double rsigma)
{
  Btree *bt;
  Ent *eval, *evec;
  void *val, *vec;
  //
  // Create the list that will contain the results.
  //
  bt = btree_Create ();
  //
  // Check for symmetry of [A]
  // Then use the correct routine.
  //
  if (mdr_IsSymmetric (a))
  {
    // ignore complex part of sigma
    mdr_ar_eigs_sym (a, &val, &vec, nev, which, rsigma);
    //
    // Hook the results into the list.
    //
    fprintf(stderr, "here\n");
    eval = ent_Create ();
    ent_SetType (eval, MATRIX_DENSE_REAL);
    ent_data (eval) = val;
    install (bt,  "val", eval);
    evec = ent_Create ();
    ent_SetType (evec, MATRIX_DENSE_REAL);
    ent_data (evec) = vec;
    install (bt,  "vec", evec);
  }
  else
  {
    mdr_ar_eigs_ge (a, &val, &vec, nev, which, rsigma);
    //
    // Hook the results into the list.
    //
    eval = ent_Create ();
    ent_SetType (eval, MATRIX_DENSE_REAL);
    ent_data (eval) = val;
    install (bt,  "val", eval);
    evec = ent_Create ();
    ent_SetType (evec, MATRIX_DENSE_REAL);
    ent_data (evec) = vec;
    install (bt,  "vec", evec);
  }

//   ent_Clean( eval );
//   ent_Clean( evec );
  return (bt);
}

/***************************************************************
 * Compute eigenvalues, and vectors for symmetric prob
 * [A]x = Lambda*x, or
 * [A]x = Lambda*x, closest to sigma
 ***************************************************************/
void
mdr_ar_eigs_sym (MDR * M, void **val, void **vec, int nev, char * which, double rsigma)
{
  int info=0, ierr=0, ldv, lworkl, m, n, i, j, ncv, maxncv;
  int *iselect=0, iparam[11], ipntr[11], *iwrkd=0;
  MDR *w=0, *workd=0, *ax=0, *d=0, *workl=0, *v=0, *v2=0, *resid=0;

  mdr_Detect_Inf (M);
  mdr_Detect_Nan (M);

  //
  // Some rudimentary checks
  //
  m = MNR (M);
  n = MNC (M);

  if (m != n)
    rerror ("eigs: input must be square");

  ldv = n;
  maxncv = MIN(4*nev,n);
  ncv    = maxncv;

  iselect = (int *) GC_MALLOC(maxncv*sizeof(int));

  lworkl = maxncv*(maxncv+8);
  workl  = mdr_Create(1, lworkl);
  d      = mdr_Create(1,2*maxncv);
  resid  = mdr_Create(1,n);
  v      = mdr_Create(ldv,maxncv);
  workd  = mdr_Create(1,3*n);
  ax     = mdr_Create(1,n);

  iparam[0] = 1;            // shift
  iparam[2] = ar_maxiter;   // maximum number of iterations

  if(rsigma != 0)
  {
    // regular mode
    iparam[6] = 3;
    //   rsigma != 0; use (A-rsigma*I)^{-1}
    MDR * A = mdr_Float_BF (M);
    iwrkd  = (int *)GC_MALLOC(n*sizeof(int));
    DSDRV1(MDRPTR(v), &ldv, MDRPTR(workl), &lworkl,
           MDRPTR(workd), MDRPTR(d), MDRPTR(resid), MDRPTR(ax), iselect, iparam, ipntr,
           MDRPTR(A), &n, &nev, &ncv, which, &rsigma, &info, &ierr, iwrkd);
    mdr_Destroy( A );
    GC_FREE(iwrkd);
  }
  else
  {
    // regular mode
    iparam[6] = 1;
    //   rsigma == 0; use A
    DSDRV1(MDRPTR(v), &ldv, MDRPTR(workl), &lworkl,
           MDRPTR(workd), MDRPTR(d), MDRPTR(resid), MDRPTR(ax), iselect, iparam, ipntr,
           MDRPTR(M), &n, &nev, &ncv, which, &rsigma, &info, &ierr, iwrkd);
  }

  if ( info < 0 )
  {
    fprintf(stdout, "Internal error in ARPACK 'dsaupd', info = %i\n", info);
    iparam[4] = 0;
  }
  if ( ierr < 0 )
  {
    fprintf(stdout, "Internal error in ARPACK 'dseupd', ierr = %i\n", ierr);
    iparam[4] = 0;
  }
  if ( info == 1 )
    fprintf(stdout, "'dsaupd' reached the maximum number of iterations.");
  if ( info == 3 )
  {
    fprintf(stdout, "No shifts could be applied during implicit Arnoldi update,\n");
    fprintf(stdout, "try increasing 'NCV'.\n");
  }

  ncv = iparam[4];    // number of converged eigenvalues

  if(ncv>0)
  {
    w = mdr_Create(1, ncv);
    for(i=0;i<ncv;i++) MdrV0(w,i) = MdrV0(d,i);

    v2 = mdr_Create(n,ncv);
    for(j=0;j<ncv;j++) for(i=0;i<n;i++) Mdr0(v2,i,j) = Mdr0(v,i,j);
  }
  else
  {
    w  = mdr_Create(0,0);
    v2 = mdr_Create(0,0);
  }

  *val = w;
  *vec = v2;

  mdr_Destroy (workl); mdr_Destroy (d); mdr_Destroy (resid);
  mdr_Destroy (workd); mdr_Destroy (v);
  mdr_Destroy (ax);
  GC_FREE(iselect);
  return;
}

void
mdr_ar_eigs_ge (MDR * M, void **val, void **vec, int nev, char * which,
                double rsigma)
{
  int info=0, ierr=0, ldv, lworkl, m, n, i, j, ncv;
  int *iselect, iparam[11], ipntr[11], *iwrkd=0;
  MDR *workd, *ax, *d, *workl, *v, *resid, *workev;
  MDR *w, *v2;

  mdr_Detect_Inf (M);
  mdr_Detect_Nan (M);

  //
  // Some rudimentary checks
  //
  m = MNR (M);
  n = MNC (M);

  if (m != n) rerror ("eigs: input must be square");

  ldv = n;
  ncv = MIN(4*nev,n);

  iselect = (int *)GC_MALLOC(ncv*sizeof(int));

  lworkl = 3*ncv*ncv+6*ncv;
  workl  = mdr_Create(1, lworkl);
  d      = mdr_Create(ncv,3);
  resid  = mdr_Create(1,n);
  v      = mdr_Create(ldv,ncv);
  workd  = mdr_Create(1,3*n);
  workev = mdr_Create(1,3*ncv);
  ax     = mdr_Create(1,n);

  iparam[1 -1] = 1;            // shift
  iparam[3 -1] = ar_maxiter;   // maximum number of iterations

  if(rsigma!=0)
  {
    // shift-invert mode
    iparam[7 -1] = 3;
    //   isigma==0,rsigma!=0; use (A-rsigma*I)^{-1}
    MDR * A = mdr_Float_BF ( M );
    iwrkd  = (int *)GC_MALLOC(n*sizeof(int));
    DNDRV1(MDRPTR(v), &ldv, MDRPTR(workl), &lworkl, MDRPTR(workev),
           MDRPTR(workd), MDRPTR(d), MDRPTR(resid), MDRPTR(ax), iselect, iparam, ipntr,
           MDRPTR(A), &n, &nev, &ncv, which, &rsigma, &info, &ierr, iwrkd);
    mdr_Destroy( A );
    GC_FREE(iwrkd);
  }
  else
  {
    // regular mode
    iparam[7 -1] = 1;
    //   isigma==0,rsigma==0; use A
    DNDRV1(MDRPTR(v), &ldv, MDRPTR(workl), &lworkl, MDRPTR(workev),
           MDRPTR(workd), MDRPTR(d), MDRPTR(resid), MDRPTR(ax), iselect, iparam, ipntr,
           MDRPTR(M), &n, &nev, &ncv, which, &rsigma, &info, &ierr, iwrkd);
  }

  if ( info < 0 )
  {
    fprintf(stdout, "Internal error in ARPACK 'dnaupd', info = %i\n", info);
    iparam[4] = 0;
  }
  if ( ierr < 0 )
  {
    fprintf(stdout, "Internal error in ARPACK 'dneupd', ierr = %i\n", ierr);
    iparam[4] = 0;
  }
  if ( info == 1 )
    fprintf(stdout, "'dsaupd' reached the maximum number of iterations.");
  if ( info == 3 )
  {
    fprintf(stdout, "No shifts could be applied during implicit Arnoldi update,\n");
    fprintf(stdout, "try increasing 'NCV'.\n");
  }

  ncv = iparam[4];    // number of converged eigenvalues

  if(ncv>0)
  {
    w = mdr_Create(1, ncv);
    for(i=0;i<ncv;i++) MdrV0(w,i) = MdrV0(d,i);
    {
      // fix the eigenvalues
      //MdrV0(w,i) = MdrV0(d,i);
      //else MdrV0(w,i) = rsigma + 1.0/MdrV0(d,i);
    }
    v2 = mdr_Create(n,ncv);
    for(j=0;j<ncv;j++) for(i=0;i<n;i++)
    {
      Mdr0(v2,i,j) = Mdr0(v,i,j);
    }
  }
  else
  {
    w  = mdr_Create(0,0);
    v2 = mdr_Create(0,0);
  }

  *val = w;
  *vec = v2;

  mdr_Destroy (workl); mdr_Destroy (d); mdr_Destroy (resid);
  mdr_Destroy (workd); mdr_Destroy (workev); mdr_Destroy (v);
  mdr_Destroy (ax); GC_FREE(iselect);
  return;
}


//
// eigs(a)
//
void
mdc_ar_eigs_ge (MDC * M, void **val, void **vec, int nev, char * which, Complex sigma);

Btree * mdc_ar_EigsS (MDC * a, int nev, char *which, Complex sigma)
{
  Btree *bt=0;
  Ent *eval, *evec;
  void *val=0, *vec=0;

  //
  // Create the list that will contain the results.
  //
  bt = btree_Create ();

  //
  // Check for symmetry of [A]
  // Then use the correct routine.
  //
  mdc_ar_eigs_ge (a, &val, &vec, nev, which, sigma);

  //
  // Hook the results into the list.
  //
  eval = ent_Create ();
  ent_SetType (eval, MATRIX_DENSE_COMPLEX);
  ent_data (eval) = val;
  install (bt,  "val", eval);
  evec = ent_Create ();
  ent_SetType (evec, MATRIX_DENSE_COMPLEX);
  ent_data (evec) = vec;
  install (bt,  "vec", evec);

//   ent_Clean( eval );
//   ent_Clean( evec );

  return (bt);
}

void
mdc_ar_eigs_ge (MDC * M, void **val, void **vec, int nev, char * which, Complex sigma)
{
  int info=0, ierr=0, ldv, lworkl, m, n, i, j, ncv;
  int *iselect, iparam[11], ipntr[14], *iwrkd=0;
  MDC *workd, *ax, *d, *workl, *v, *resid, *workev;
  MDC *w, *v2;
  MDR *rwork, *rd;

  mdc_Detect_Inf (M);
  mdc_Detect_Nan (M);

  //
  // Some rudimentary checks
  //
  m = MNR (M);
  n = MNC (M);

  if (m != n) rerror ("eigs: input must be square");

  ldv = n;
  ncv = MIN(4*nev,n);

  iselect = (int *) GC_MALLOC(ncv*sizeof(int));

  lworkl = 3*ncv*ncv+5*ncv;
  workl  = mdc_Create(1, lworkl);
  d      = mdc_Create(ncv,3);
  resid  = mdc_Create(1,n);
  v      = mdc_Create(ldv,ncv);
  workd  = mdc_Create(1,3*n);
  workev = mdc_Create(1,3*ncv);
  ax     = mdc_Create(1,n);

  rwork  = mdr_Create(1,ncv);
  rd     = mdr_Create(ncv,3);

  iparam[1 -1] = 1;            // shift
  iparam[3 -1] = ar_maxiter;   // maximum number of iterations

  if(cabs(sigma)==0)
  {
    // regular mode
    iparam[7 -1] = 1;
    //   isigma==0,rsigma==0; use A
    ZNDRV1( MDCPTR(v), &ldv, MDCPTR(workl), &lworkl, MDCPTR(workev),
            MDCPTR(workd), MDCPTR(d), MDCPTR(resid), MDCPTR(ax), iselect, iparam, ipntr,
            MDCPTR(M), &n, &nev, &ncv, which, &sigma, &info, &ierr, iwrkd,
            MDRPTR(rwork), MDRPTR(rd)
          );
  }
  else
  {
    // shift-invert mode
    iparam[7 -1] = 3;
    //   isigma==0,rsigma!=0; use (A-rsigma*I)^{-1}
    MDC * A = mdc_Copy( M );
    iwrkd  = (int *)GC_MALLOC(n*sizeof(int));
    ZNDRV1(MDCPTR(v), &ldv, MDCPTR(workl), &lworkl, MDCPTR(workev),
            MDCPTR(workd), MDCPTR(d), MDCPTR(resid), MDCPTR(ax), iselect, iparam, ipntr,
            MDCPTR(A), &n, &nev, &ncv, which, &sigma, &info, &ierr, iwrkd,
            MDRPTR(rwork), MDRPTR(rd)
          );
    mdc_Destroy( A );
    GC_FREE(iwrkd);
  }


  if ( info < 0 )
  {
    fprintf(stdout, "Internal error in ARPACK 'znaupd', info = %i\n", info);
    iparam[4] = 0;
  }
  if ( ierr < 0 )
  {
    fprintf(stdout, "Internal error in ARPACK 'zneupd', ierr = %i\n", ierr);
    iparam[4] = 0;
  }
  if ( info == 1 )
    fprintf(stdout, "'zsaupd' reached the maximum number of iterations.");
  if ( info == 3 )
  {
    fprintf(stdout, "No shifts could be applied during implicit Arnoldi update,\n");
    fprintf(stdout, "try increasing 'NCV'.\n");
  }

  ncv = iparam[4];    // number of converged eigenvalues

  if(ncv>0)
  {
    w = mdc_Create(1, ncv);
    for(i=0;i<ncv;i++)
      MdcV0(w,i) = MdcV0(d,i);
    v2 = mdc_Create(n,ncv);
    for(j=0;j<ncv;j++)
    {
      Complex dummy;
      dummy =  conj(Mdc0(v,0,j));
      if( cabs( dummy ) > 0)
      {
        //
        // normalize eigenvectors so that first coefficient is real (LAPACK)
        //
        Mdc0(v2,0,j) = cabs(dummy);
        for(i=1;i<n;i++)
        {
          Mdc0(v2,i,j)  = Mdc0(v,i,j) * dummy / Mdc0(v2,0,j);
        }
      }
      else
      {
        //
        // don't do anything. First coefficient of eigenvector is zero
        //
        for(i=0;i<n;i++)
        {
          Mdc0(v2,i,j)  = Mdc0(v,i,j);
        }
      }
    }
  }
  else
  {
    w  = mdc_Create(0,0);
    v2 = mdc_Create(0,0);
  }

  *val = w;
  *vec = v2;

  mdc_Destroy (workl); mdc_Destroy (d); mdc_Destroy (resid);
  mdc_Destroy (workd); mdc_Destroy (workev); mdc_Destroy (v);
  mdc_Destroy (ax); mdr_Destroy(rd); mdr_Destroy(rwork);
  GC_FREE(iselect);

  return;
}



void
msr_ar_eigs_sym (MSR * M, void **val, void **vec, int nev, char * which, double sigma)
{
  int info=0, ierr=0, ldv, lworkl, m, n, i, j, ncv;
  int *iselect, iparam[11], ipntr[11];
  MDR *w, *workd, *d, *workl, *v, *v2, *resid;

  char bmat='I';
  double tol=0;
  int rvec=1, ido=0;

  msr_Detect_Inf (M);

  //
  // Some rudimentary checks
  //
  m = MNR (M);
  n = MNC (M);

  if (m != n) rerror ("eigs: input must be square");

  ldv = n;
  ncv = MIN(4*nev,n);

  iselect = (int *)GC_MALLOC(ncv*sizeof(int));

  lworkl = ncv*(ncv+8);
  workl  = mdr_Create(1, lworkl);
  d      = mdr_Create(1,2*ncv);
  resid  = mdr_Create(1,n);
  v      = mdr_Create(ldv,ncv);
  workd  = mdr_Create(1,3*n);

  iparam[0] = 1;            // shift
  iparam[2] = ar_maxiter;   // maximum number of iterations

  MDR *ax=0;
  MDR *x  = mdr_Create(0,0);

  if(sigma != 0)
  {
    //
    // shift-invert mode
    //
    iparam[7 -1] = 3;

    //
    // first call to spsolve, just remember the LU factorization of 'a - sigma*I'
    //
    MSR * msr_sigma;
    MSR * aminuss;
    msr_sigma = msr_Create(n,n);
    msr_Setup(msr_sigma, n);
    msr_sigma->ia[0] = 1;
    for (i = 1; i <= n; i++)
    {
      msr_sigma->ia[i]   = i+1;
      msr_sigma->ja[i-1] = i;
      msr_sigma->d [i-1] = sigma;
    }

    aminuss = msr_Subtract(M, msr_sigma);
    superlu_msr_spsolve(aminuss, mdr_Create(0,0) );

tenone:

    DSAUPD (&ido, &bmat, &n, which, &nev, &tol, MDRPTR(resid),
             &ncv, MDRPTR(v), &ldv, iparam, ipntr, MDRPTR(workd), MDRPTR(workl),
             &lworkl, &info );
    if (ido == -1 || ido == 1)
    {
      // x is a dummy variable
      x->nrow = n;
      x->ncol = 1;
      MDPTR(x) = (void *) &MdrV0(workd, ipntr[1 -1] -1);
      ax = (MDR *) superlu_msr_spsolve( msr_Create(0,0), x);
      //
      // ax = A.x presumably worked. Copy result back.
      //
      for(i = 0; i<n; i++)
        MdrV0(workd, ipntr[2 -1] -1 + i) = MdrV0(ax,i);
      goto tenone;
    }
    // clean-up the local sparse matrices
    msr_Destroy( aminuss );
    msr_Destroy( msr_sigma );
  }
  else
  {
    //
    // regular mode
    //
    int rtype;
    iparam[7 -1] = 1;

tentwo:
    DSAUPD (&ido, &bmat, &n, which, &nev, &tol, MDRPTR(resid),
            &ncv, MDRPTR(v), &ldv, iparam, ipntr, MDRPTR(workd), MDRPTR(workl),
            &lworkl, &info );
    if (ido == -1 || ido == 1)
    {
      // x is a dummy variable
      x->nrow = n;
      x->ncol = 1;
      MDPTR(x) = (void *) &MdrV0(workd, ipntr[1 -1] -1);
      ax = msr_mdr_Multiply(M, x, &rtype);
      if(rtype==MATRIX_DENSE_REAL)
      {
        //
        // ax = A.x worked. Copy result back.
        //
        for(i = 0; i<n; i++)
          MdrV0(workd, ipntr[2 -1] -1 + i) = MdrV0(ax,i);
      }
      else
        rerror("eigs: step  ax = A . x  did not work!");
      goto tentwo;
    }
  }
  mdr_Destroy(ax);

  if (info >= 0 )
  {
    DSEUPD ( &rvec, "All", iselect, MDRPTR(d),
              MDRPTR(v), &ldv, &sigma, &bmat,
              &n, which, &nev, &tol,
              MDRPTR(resid), &ncv, MDRPTR(v), &ldv,
              iparam, ipntr, MDRPTR(workd), MDRPTR(workl),
              &lworkl, &ierr );
    if ( ierr < 0 )
    {
      fprintf(stdout, "Internal error in ARPACK 'dseupd', ierr = %i\n", ierr);
      iparam[5 -1] = 0;
    }
  }

  if ( info == 1 )
    fprintf(stdout, "'dsaupd' reached the maximum number of iterations.");
  if ( info == 3 )
  {
    fprintf(stdout, "No shifts could be applied during implicit Arnoldi update,\n");
    fprintf(stdout, "try increasing 'NCV'.\n");
  }

  ncv = iparam[5 -1];    // number of converged eigenvalues
  if(ncv>0)
  {
    w = mdr_Create(1, ncv);
    for(i=0;i<ncv;i++) MdrV0(w,i) = MdrV0(d,i);
    v2 = mdr_Create(n,ncv);
    for(j=0;j<ncv;j++) for(i=0;i<n;i++) Mdr0(v2,i,j) = Mdr0(v,i,j);
  }
  else
  {
    w  = mdr_Create(0,0);
    v2 = mdr_Create(0,0);
  }

  *val = w;
  *vec = v2;

  mdr_Destroy (workl);
  mdr_Destroy (d);
  mdr_Destroy (resid);
  mdr_Destroy (workd);
  mdr_Destroy (v);
  GC_FREE(iselect);
  return;
}

void
msr_ar_eigs_ge (MSR * M, void **val, void **vec, int nev, char * which, double sigma)
{
  int info=0, ierr=0, ldv, lworkl, m, n, i, j, ncv;
  int *iselect, iparam[11], ipntr[11];
  MDR *w, *workd, *d, *workl, *v, *v2, *resid, *workev;

  char bmat='I', what='A';
  double tol=0, zero=0;
  int rvec=1, ido=0;

  msr_Detect_Inf (M);

  //
  // Some rudimentary checks
  //
  m = MNR (M);
  n = MNC (M);

  if (m != n) rerror ("eigs: input must be square");

  ldv = n;
  ncv = MIN(4*nev,n);

  iselect = (int *)GC_MALLOC(ncv*sizeof(int));

  lworkl = 3*ncv*ncv+6*ncv;
  workl  = mdr_Create(1, lworkl);
  d      = mdr_Create(ncv,3);
  resid  = mdr_Create(1,n);
  v      = mdr_Create(ldv,ncv);
  workev = mdr_Create(1,3*ncv);
  workd  = mdr_Create(1,3*n);

  iparam[0] = 1;            // shift
  iparam[2] = ar_maxiter;   // maximum number of iterations

  MDR *ax=0;
  MDR *x  = mdr_Create(0,0);

  if(sigma != 0)
  {
    //
    // shift-invert mode
    //
    iparam[7 -1] = 3;

    //
    // first call to spsolve, just remember the LU factorization of 'a - sigma*I'
    //
    MSR * msr_sigma;
    MSR * aminuss;
    msr_sigma = msr_Create(n,n);
    msr_Setup(msr_sigma, n);
    msr_sigma->ia[0] = 1;
    for (i = 1; i <= n; i++)
    {
      msr_sigma->ia[i]   = i+1;
      msr_sigma->ja[i-1] = i;
      msr_sigma->d [i-1] = sigma;
    }

    aminuss = msr_Subtract(M, msr_sigma);
    superlu_msr_spsolve(aminuss, mdr_Create(0,0) );

tenone:

    DNAUPD (&ido, &bmat, &n, which, &nev, &tol, MDRPTR(resid),
             &ncv, MDRPTR(v), &ldv, iparam, ipntr, MDRPTR(workd), MDRPTR(workl),
             &lworkl, &info );

    if (ido == -1 || ido == 1)
    {
      // x is a dummy variable
      x->nrow = n;
      x->ncol = 1;
      MDPTR(x) = (void *) &MdrV0(workd, ipntr[1 -1] -1);
      ax = (MDR *) superlu_msr_spsolve( msr_Create(0,0), x);
      //
      // ax = A.x presumably worked. Copy result back.
      //
      for(i = 0; i<n; i++)
        MdrV0(workd, ipntr[2 -1] -1 + i) = MdrV0(ax,i);
      goto tenone;
    }
    // clean-up the local sparse matrices
    msr_Destroy( aminuss );
    msr_Destroy( msr_sigma );
  }
  else
  {
    //
    // regular mode
    //
    int rtype;
    iparam[7 -1] = 1;

tentwo:
    DNAUPD (&ido, &bmat, &n, which, &nev, &tol, MDRPTR(resid),
             &ncv, MDRPTR(v), &ldv, iparam, ipntr, MDRPTR(workd), MDRPTR(workl),
             &lworkl, &info );
    if (ido == -1 || ido == 1)
    {
      // x is a dummy variable
      x->nrow = n;
      x->ncol = 1;
      MDPTR(x) = (void *) &MdrV0(workd, ipntr[1 -1] -1);
      ax = msr_mdr_Multiply(M, x, &rtype);
      if(rtype==MATRIX_DENSE_REAL)
      {
        //
        // ax = A.x worked. Copy result back.
        //
        for(i = 0; i<n; i++)
          MdrV0(workd, ipntr[2 -1] -1 + i) = MdrV0(ax,i);
      }
      else
        rerror("eigs: step  ax = A . x  did not work!");
      goto tentwo;
    }
  }
  mdr_Destroy(ax);

  if (info >= 0 )
  {
    DNEUPD (&rvec, &what, iselect, MDRPTR(d), &Mdr1(d,1,2), MDRPTR(v), &ldv,
            &sigma, &zero, MDRPTR(workev), &bmat, &n, which, &nev, &tol,
             MDRPTR(resid), &ncv, MDRPTR(v), &ldv, iparam, ipntr, MDRPTR(workd), MDRPTR(workl),
            &lworkl, &ierr );
    if ( ierr < 0 )
    {
      fprintf(stdout, "Internal error in ARPACK 'dseupd', ierr = %i\n", ierr);
      iparam[5 -1] = 0;
    }
  }

  if ( info == 1 )
    fprintf(stdout, "'dsaupd' reached the maximum number of iterations.");
  if ( info == 3 )
  {
    fprintf(stdout, "No shifts could be applied during implicit Arnoldi update,\n");
    fprintf(stdout, "try increasing 'NCV'.\n");
  }

  ncv = iparam[5 -1];    // number of converged eigenvalues
  if(ncv>0)
  {
    w = mdr_Create(1, ncv);
    for(i=0;i<ncv;i++) MdrV0(w,i) = MdrV0(d,i);
    v2 = mdr_Create(n,ncv);
    for(j=0;j<ncv;j++) for(i=0;i<n;i++) Mdr0(v2,i,j) = Mdr0(v,i,j);
  }
  else
  {
    w  = mdr_Create(0,0);
    v2 = mdr_Create(0,0);
  }

  *val = w;
  *vec = v2;

  mdr_Destroy (workl);
  mdr_Destroy (d);
  mdr_Destroy (resid);
  mdr_Destroy (workd);
  mdr_Destroy (workev);
  mdr_Destroy (v);
  GC_FREE(iselect);
  return;
}


Btree * msr_ar_EigsS (MSR * a, int nev, char *which, double sigma)
{
  Btree *bt;
  Ent *eval, *evec;
  void *val, *vec;

  //
  // Create the list that will contain the results.
  //
  bt = btree_Create ();
  //
  // Check for symmetry of [A]
  // Then use the correct routine.
  //
  if (msr_IsSymmetric (a))
  {
    // ignore complex part of sigma
    msr_ar_eigs_sym (a, &val, &vec, nev, which, sigma);
    //
    // Hook the results into the list.
    //
    eval = ent_Create ();
    ent_SetType (eval, MATRIX_DENSE_REAL);
    ent_data (eval) = val;
    install (bt,  "val", eval);
    evec = ent_Create ();
    ent_SetType (evec, MATRIX_DENSE_REAL);
    ent_data (evec) = vec;
    install (bt,  "vec", evec);
  }
  else
  {
    msr_ar_eigs_ge (a, &val, &vec, nev, which, sigma);
    //
    // Hook the results into the list.
    //
    eval = ent_Create ();
    ent_SetType (eval, MATRIX_DENSE_REAL);
    ent_data (eval) = val;
    install (bt,  "val", eval);
    evec = ent_Create ();
    ent_SetType (evec, MATRIX_DENSE_REAL);
    ent_data (evec) = vec;
    install (bt,  "vec", evec);
  }

//   ent_Clean( eval );
//   ent_Clean( evec );
  return (bt);
}

void msc_ar_eigs_ge (MSC * M, void **val, void **vec, int nev, char * which, Complex sigma)
{
  int info=0, ierr=0, ldv, lworkl, m, n, i, j, ncv;
  int *iselect, iparam[11], ipntr[11];
  MDC *w, *workd, *d, *workl, *v, *v2, *resid, *workev;
  MDR *rwork;

  char bmat='I', what='A';
  double tol=0;
  int rvec=1, ido=0;

  msc_Detect_Inf (M);
  msc_Detect_NaN (M);

  //
  // Some rudimentary checks
  //
  m = MNR (M);
  n = MNC (M);

  if (m != n) rerror ("eigs: input must be square");

  ldv = n;
  ncv = MIN(4*nev,n);

  iselect = (int *)GC_MALLOC(ncv*sizeof(int));

  lworkl = 3*ncv*ncv+6*ncv;
  workl  = mdc_Create(1, lworkl);
  d      = mdc_Create(ncv,3);
  resid  = mdc_Create(1,n);
  v      = mdc_Create(ldv,ncv);
  workev = mdc_Create(1,3*ncv);
  workd  = mdc_Create(1,3*n);

  rwork  = mdr_Create(1,ncv);

  iparam[0] = 1;            // shift
  iparam[2] = ar_maxiter;   // maximum number of iterations

  MDC *ax=0;
  MDC *x  = mdc_Create(0,0);

  if( cabs(sigma) != 0)
  {
    //
    // shift-invert mode
    //
    iparam[7 -1] = 3;

    //
    // first call to spsolve, just remember the LU factorization of 'a - sigma*I'
    //
    MSC * msc_sigma;
    MSC * aminuss;
    msc_sigma = msc_Create(n,n);
    msc_Setup(msc_sigma, n);
    msc_sigma->ia[0] = 1;
    for (i = 1; i <= n; i++)
    {
      msc_sigma->ia[i]   = i+1;
      msc_sigma->ja[i-1] = i;
      msc_sigma->c [i-1] = sigma;
    }
    aminuss = msc_Subtract(M, msc_sigma);
    //superlu_msc_spsolve(aminuss, mdc_Create(0,0) );

tenone:
    ZNAUPD (&ido, &bmat, &n, which, &nev, &tol, MDCPTR(resid),
             &ncv, MDCPTR(v), &ldv, iparam, ipntr, MDCPTR(workd), MDCPTR(workl), &lworkl,
             MDRPTR(rwork), &info );

    if (ido == -1 || ido == 1)
    {
      // x is a dummy variable
      x->nrow = n;
      x->ncol = 1;
      MDCPTR(x) = &MdcV0(workd, ipntr[1 -1] -1);
      ax = (MDC *) umfpack_msc_SolveCC( aminuss, x, 0);
      //
      // ax = A.x presumably worked. Copy result back.
      //
      for(i = 0; i<n; i++)
        MdcV0(workd, ipntr[2 -1] -1 + i) = MdcV0(ax,i);
      goto tenone;
    }
    // clean-up the local sparse matrices
    msc_Destroy( aminuss );
    msc_Destroy( msc_sigma );
  }
  else
  {
    //
    // regular mode
    //
    int rtype;
    iparam[7 -1] = 1;

tentwo:
    ZNAUPD (&ido, &bmat, &n, which, &nev, &tol, MDCPTR(resid),
             &ncv, MDCPTR(v), &ldv, iparam, ipntr, MDCPTR(workd), MDCPTR(workl), &lworkl,
             MDRPTR(rwork), &info );

    if (ido == -1 || ido == 1)
    {
      // x is a dummy variable
      x->nrow = n;
      x->ncol = 1;
      MDCPTR(x) = MDCPTR(workd) + (ipntr[1 -1] -1);
      ax = (MDC *) msc_mdc_Multiply(M, x, &rtype);
      if(rtype==MATRIX_DENSE_COMPLEX)
      {
        //
        // ax = A.x worked. Copy result back.
        //
        for(i = 0; i<n; i++)
          MdcV0(workd, ipntr[2 -1] -1 + i) = MdcV0(ax,i);
      }
      else
        rerror("eigs: step  ax = A . x  did not work!");
      goto tentwo;
    }
  }
  mdc_Destroy(ax);

  if (info >= 0 )
  {
    ZNEUPD ( &rvec, &what, iselect, MDCPTR(d), MDCPTR(v), &ldv, &sigma,
              MDCPTR(workev), &bmat, &n, which, &nev, &tol, MDCPTR(resid), &ncv,
              MDCPTR(v), &ldv, iparam, ipntr, MDCPTR(workd), MDCPTR(workl), &lworkl,
              MDRPTR(rwork), &ierr );
    if ( ierr < 0 )
    {
      fprintf(stdout, "Internal error in ARPACK 'zneupd', ierr = %i\n", ierr);
      iparam[5 -1] = 0;
    }
  }

  if ( info == 1 )
    fprintf(stdout, "'znaupd' reached the maximum number of iterations.");
  if ( info == 3 )
  {
    fprintf(stdout, "No shifts could be applied during implicit Arnoldi update,\n");
    fprintf(stdout, "try increasing 'NCV'.\n");
  }

  ncv = iparam[5 -1];    // number of converged eigenvalues
  if(ncv>0)
  {
    w = mdc_Create(1, ncv);
    for(i=0;i<ncv;i++) MdcV0(w,i) = MdcV0(d,i);
    v2 = mdc_Create(n,ncv);
    for(j=0;j<ncv;j++) for(i=0;i<n;i++) Mdc0(v2,i,j) = Mdc0(v,i,j);
  }
  else
  {
    w  = mdc_Create(0,0);
    v2 = mdc_Create(0,0);
  }

  *val = w;
  *vec = v2;

  mdc_Destroy (workl);
  mdc_Destroy (d);
  mdc_Destroy (resid);
  mdc_Destroy (workd);
  mdc_Destroy (workev);
  mdc_Destroy (v);
  mdr_Destroy (rwork);
  GC_FREE(iselect);
  return;
}

Btree * msc_ar_EigsS (MSC * a, int nev, char *which, Complex sigma)
{
  Btree *bt = btree_Create ();
  Ent *eval=0, *evec=0;
  void *val=0, *vec=0;

  //
  // Create the list that will contain the results.
  //
  msc_ar_eigs_ge (a, &val, &vec, nev, which, sigma);

  //
  // Hook the results into the list.
  //
  eval = ent_Create ();
  ent_SetType (eval, MATRIX_DENSE_COMPLEX);
  ent_data (eval) = val;
  install (bt,  "val", eval);
  evec = ent_Create ();
  ent_SetType (evec, MATRIX_DENSE_COMPLEX);
  ent_data (evec) = vec;
  install (bt,  "vec", evec);
  return (bt);
}

void
msr_ar_eigs_ge2 (MSR * A, MSR *B, void **val, void **vec, int nev, char * which, double sigma)
{
  int info=0, ierr=0, ldv, lworkl, n, i, j, ncv;
  int *iselect, iparam[11], ipntr[11];
  MDR *w, *workd, *d, *workl, *v, *v2, *resid, *workev;

  char bmat='I', what='A';
  double tol=0, zero=0;
  int rvec=1, ido=0;

  n = MNC (A);

  //
  // Some rudimentary checks
  //
  msr_Detect_Inf (A);
  msr_Detect_Inf (B);
  msr_Detect_NaN (A);
  msr_Detect_NaN (B);

  if (MNR (A) != MNC (A)) rerror ("eigs: A must be a square matrix!");
  if (MNR (B) != MNC (B)) rerror ("eigs: B must be a square matrix!");
  if (MNR (A) != MNC (B)) rerror ("eigs: A and B must be the matrices of the same size!");

  ldv = n;
  ncv = MIN(4*nev,n);

  iselect = (int *)GC_MALLOC(ncv*sizeof(int));

  lworkl = 3*ncv*ncv+6*ncv;
  workl  = mdr_Create(1, lworkl);
  d      = mdr_Create(ncv,3);
  resid  = mdr_Create(1,n);
  v      = mdr_Create(ldv,ncv);
  workev = mdr_Create(1,3*ncv);
  workd  = mdr_Create(1,3*n);

  iparam[0] = 1;            // shift
  iparam[2] = ar_maxiter;   // maximum number of iterations

  MDR *ax=0;
  MDR *x  = mdr_Create(0,0);

  //
  // shift-invert mode
  //
  iparam[7 -1] = 3;

  //
  // first call to spsolve, just remember the LU factorization of 'a - sigma*b'
  //
  if(sigma != 0)
  {
    //
    // sigma !=0
    //
    MSR * msr_sigmab = msr_Copy( B );
    for (i = 1; i <= B->nnz; i++)
      msr_sigmab->d [i-1] *= sigma;
    MSR *aminussb = msr_Subtract(A, msr_sigmab);
    superlu_msr_spsolve(aminussb, mdr_Create(0,0) );
    msr_Destroy( aminussb );
    msr_Destroy( msr_sigmab );
  }
  else
    superlu_msr_spsolve( A, mdr_Create(0,0) );

tenone:

    DNAUPD (&ido, &bmat, &n, which, &nev, &tol, MDRPTR(resid),
            &ncv, MDRPTR(v), &ldv, iparam, ipntr, MDRPTR(workd), MDRPTR(workl),
          &lworkl, &info );

  if (ido == 2)
  {
    int rtype;
    //
    // (1) <- B . (1)
    //
    x->nrow = n;
    x->ncol = 1;
    MDPTR(x) = (void *) &MdrV0(workd, ipntr[1 -1] -1);
    ax = msr_mdr_Multiply(B, x, &rtype);
    if(rtype==MATRIX_DENSE_REAL)
    {
      //
      // ax -> (1)
      //
      for(i = 0; i<n; i++)
        MdrV0(workd, ipntr[1 -1] -1 + i) = MdrV0(ax,i);
    }
    else rerror("eigs: terrible internal error!");
    //
    // (2) <- (a-SIGMA*b)^(-1) . (1)
    //
    ax = (MDR *) superlu_msr_spsolve( msr_Create(0,0), ax);
    //
    // ax = A.x presumably worked. Copy result back.
    //
    for(i = 0; i<n; i++)
      MdrV0(workd, ipntr[2 -1] -1 + i) = MdrV0(ax,i);
    goto tenone;
  }
  else if (ido == 1)
  {
    //
    // (2) <- (A-sigma*B)^(-1) . (3)
    //
    // x is a dummy variable
    x->nrow = n;
    x->ncol = 1;
    MDPTR(x) = (void *) &MdrV0(workd, ipntr[3 -1] -1);
    ax = (MDR *) superlu_msr_spsolve( msr_Create(0,0), x);
    //
    // ax = A.x presumably worked. Copy result back.
    //
    for(i = 0; i<n; i++)
      MdrV0(workd, ipntr[2 -1] -1 + i) = MdrV0(ax,i);
    goto tenone;
  }
  else if (ido == 2)
  {
    //
    // (2) <- B . (1)
    //
    int rtype;

    // x is a dummy variable
    x->nrow = n;
    x->ncol = 1;
    MDPTR(x) = (void *) &MdrV0(workd, ipntr[1 -1] -1);
    ax = msr_mdr_Multiply(B, x, &rtype);

    if(rtype==MATRIX_DENSE_REAL)
    {
      //
      // ax = A.x worked. Copy result back.
      //
      for(i = 0; i<n; i++)
        MdrV0(workd,ipntr[2 -1] -1 + i) = MdrV0(ax,i);
    }
    else rerror("eigs: terrible internal error!");
  }

  if (info >= 0 )
  {
    DNEUPD ( &rvec, &what, iselect, MDRPTR(d), &Mdr1(d,1,2), MDRPTR(v), &ldv,
              &sigma, &zero, MDRPTR(workev), &bmat, &n, which, &nev, &tol,
              MDRPTR(resid), &ncv, MDRPTR(v), &ldv, iparam, ipntr, MDRPTR(workd), MDRPTR(workl),
              &lworkl, &ierr );
    if ( ierr < 0 )
    {
      fprintf(stdout, "Internal error in ARPACK 'dseupd', ierr = %i\n", ierr);
      iparam[5 -1] = 0;
    }
  }

  if ( info == 1 )
    fprintf(stdout, "'dsaupd' reached the maximum number of iterations.");
  if ( info == 3 )
  {
    fprintf(stdout, "No shifts could be applied during implicit Arnoldi update,\n");
    fprintf(stdout, "try increasing 'NCV'.\n");
  }

  ncv = iparam[5 -1];    // number of converged eigenvalues
  if(ncv>0)
  {
    w = mdr_Create(1, ncv);
    for(i=0;i<ncv;i++) MdrV0(w,i) = MdrV0(d,i);
    v2 = mdr_Create(n,ncv);
    for(j=0;j<ncv;j++) for(i=0;i<n;i++) Mdr0(v2,i,j) = Mdr0(v,i,j);
  }
  else
  {
    w  = mdr_Create(0,0);
    v2 = mdr_Create(0,0);
  }

  *val = w;
  *vec = v2;

  mdr_Destroy (workl);
  mdr_Destroy (d);
  mdr_Destroy (resid);
  mdr_Destroy (workd);
  mdr_Destroy (workev);
  mdr_Destroy (v);
  mdr_Destroy(ax);
  GC_FREE(iselect);
  return;
}


Btree * msr_ar_EigsG (MSR * a, MSR *b, int nev, char *which, double sigma)
{
  Btree *bt;
  Ent *eval, *evec;
  void *val, *vec;

  //
  // Create the list that will contain the results.
  //
  bt = btree_Create ();
  //
  // Check for symmetry of [A]
  // Then use the correct routine.
  //
  msr_ar_eigs_ge2 (a, b, &val, &vec, nev, which, sigma);
  //
  // Hook the results into the list.
  //
  eval = ent_Create ();
  ent_SetType (eval, MATRIX_DENSE_REAL);
  ent_data (eval) = val;
  install (bt,  "val", eval);
  evec = ent_Create ();
  ent_SetType (evec, MATRIX_DENSE_REAL);
  ent_data (evec) = vec;
  install (bt,  "vec", evec);

//   ent_Clean( eval );
//   ent_Clean( evec );
  return (bt);
}

void msc_ar_eigs_ge2 (MSC * A, MSC * B, void **val, void **vec, int nev, char * which, Complex sigma)
{
  int info=0, ierr=0, ldv, lworkl, n, i, j, ncv;
  int *iselect, iparam[11], ipntr[11];
  MDC *w, *workd, *d, *workl, *v, *v2, *resid, *workev;
  MDR *rwork;

  char bmat='I', what='A';
  double tol=0;
  int rvec=1, ido=0;

  //
  // Some rudimentary checks
  //
  msc_Detect_Inf (A);
  msc_Detect_Inf (B);
  msc_Detect_NaN (A);
  msc_Detect_NaN (B);

  if (MNR (A) != MNC (A)) rerror ("eigs: A must be a square matrix!");
  if (MNR (B) != MNC (B)) rerror ("eigs: B must be a square matrix!");
  if (MNR (A) != MNC (B)) rerror ("eigs: A and B must be the matrices of the same size!");

  n = MNC (A);

  ldv = n;
  ncv = MIN(4*nev,n);

  iselect = (int *)GC_MALLOC(ncv*sizeof(int));

  lworkl = 3*ncv*ncv+6*ncv;
  workl  = mdc_Create(1, lworkl);
  d      = mdc_Create(ncv,3);
  resid  = mdc_Create(1,n);
  v      = mdc_Create(ldv,ncv);
  workev = mdc_Create(1,3*ncv);
  workd  = mdc_Create(1,3*n);

  rwork  = mdr_Create(1,ncv);

  iparam[0] = 1;            // shift
  iparam[2] = ar_maxiter;   // maximum number of iterations

  MDC *ax=0;
  MDC *x  = mdc_Create(0,0);

  //
  // shift-invert mode for any sigma
  //
  iparam[7 -1] = 3;

  MSC * aminussb;

  if(cabs(sigma) != 0)
  {
    MSC * msc_sigmab = msc_Copy( B );
    for (i = 1; i <= B->nnz; i++)
    {
      msc_sigmab->c[i-1] = sigma * msc_sigmab->c[i-1];
    }
    aminussb = msc_Subtract(A, msc_sigmab);
    msc_Destroy( msc_sigmab );
  }
  else aminussb = A;

tenone:

  ZNAUPD (&ido, &bmat, &n, which, &nev, &tol, MDCPTR(resid),
          &ncv, MDCPTR(v), &ldv, iparam, ipntr, MDCPTR(workd), MDCPTR(workl), &lworkl,
          MDRPTR(rwork), &info );

    if (ido == -1)
    {
      //
      // (1) <- B . (1)
      //
      int rtype;
      x->nrow = n;
      x->ncol = 1;
      MDCPTR(x) = MDCPTR(workd) + (ipntr[1 -1] -1);
      ax = (MDC *) msc_mdc_Multiply(B, x, &rtype);
      //
      // ax = A.x presumably worked. Copy result back to (1)
      //
      for(i = 0; i<n; i++)
        MdcV0(workd, ipntr[1 -1] -1 + i) = MdcV0(ax,i);
      //
      // (2) <- (A-sigma*B)^(-1) . (1)
      //
      ax = (MDC *) umfpack_msc_SolveCC( aminussb, ax, 0);
      if(rtype==MATRIX_DENSE_COMPLEX)
      {
        //
        // ax = A.x worked. Copy result back to (2)
        //
        for(i = 0; i<n; i++)
          MdcV0(workd, ipntr[2 -1] -1 + i) = MdcV0(ax,i);
      }
      else
        rerror("eigs: step  ax = A . x  did not work!");
      goto tenone;
    }
    else if (ido == 1)
    {
      //
      // (2) <- (A-sigma*B)^(-1) . (3)
      //
      x->nrow = n;
      x->ncol = 1;
      MDCPTR(x) = &MdcV0(workd, ipntr[3 -1] -1);
      ax = (MDC *) umfpack_msc_SolveCC( aminussb, x, 0);
      //
      // ax = A.x worked. Copy result back to (2)
      //
      for(i = 0; i<n; i++)
        MdcV0(workd, ipntr[2 -1] -1 + i) = MdcV0(ax,i);
      goto tenone;
    }
    else if (ido == 2)
    {
      //
      // (2) <- B . (1)
      //
      int rtype;
      x->nrow = n;
      x->ncol = 1;
      MDCPTR(x) = &MdcV0(workd, ipntr[1 -1] -1);
      ax = (MDC *) msc_mdc_Multiply(B, x, &rtype);
      //
      // ax = A.x presumably worked. Copy result back to (2)
      //
      for(i = 0; i<n; i++)
        MdcV0(workd, ipntr[2 -1] -1 + i) = MdcV0(ax,i);
      goto tenone;
    }


  if (info >= 0 )
  {
    ZNEUPD ( &rvec, &what, iselect, MDCPTR(d), MDCPTR(v), &ldv, &sigma,
              MDCPTR(workev), &bmat, &n, which, &nev, &tol, MDCPTR(resid), &ncv,
              MDCPTR(v), &ldv, iparam, ipntr, MDCPTR(workd), MDCPTR(workl), &lworkl,
              MDRPTR(rwork), &ierr );
    if ( ierr < 0 )
    {
      fprintf(stdout, "Internal error in ARPACK 'zneupd', ierr = %i\n", ierr);
      iparam[5 -1] = 0;
    }
  }

  if ( info == 1 )
    fprintf(stdout, "'znaupd' reached the maximum number of iterations.");
  if ( info == 3 )
  {
    fprintf(stdout, "No shifts could be applied during implicit Arnoldi update,\n");
    fprintf(stdout, "try increasing 'NCV'.\n");
  }

  ncv = iparam[5 -1];    // number of converged eigenvalues
  if(ncv>0)
  {
    w = mdc_Create(1, ncv);
    for(i=0;i<ncv;i++) MdcV0(w,i) = MdcV0(d,i);
    v2 = mdc_Create(n,ncv);
    for(j=0;j<ncv;j++) for(i=0;i<n;i++) Mdc0(v2,i,j) = Mdc0(v,i,j);
  }
  else
  {
    w  = mdc_Create(0,0);
    v2 = mdc_Create(0,0);
  }

  *val = w;
  *vec = v2;

  mdc_Destroy (workl);
  mdc_Destroy (d);
  mdc_Destroy (resid);
  mdc_Destroy (workd);
  mdc_Destroy (workev);
  mdc_Destroy (v);
  mdc_Destroy(ax);
  mdr_Destroy (rwork);
  GC_FREE(iselect);
  return;
}

Btree * msc_ar_EigsG (MSC * a, MSC * b, int nev, char *which, Complex sigma)
{
  Btree *bt;
  Ent *eval, *evec;
  void *val, *vec;

  //
  // Create the list that will contain the results.
  //
  bt = btree_Create ();
  msc_ar_eigs_ge2 (a, b, &val, &vec, nev, which, sigma);

  //
  // Hook the results into the list.
  //
  eval = ent_Create ();
  ent_SetType (eval, MATRIX_DENSE_COMPLEX);
  ent_data (eval) = val;
  install (bt,  "val", eval);
  evec = ent_Create ();
  ent_SetType (evec, MATRIX_DENSE_COMPLEX);
  ent_data (evec) = vec;
  install (bt,  "vec", evec);

  return (bt);
}

#endif
