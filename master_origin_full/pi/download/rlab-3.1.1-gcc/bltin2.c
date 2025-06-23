/* bltin2.c */

/*  This file is a part of RLaB ("Our"-LaB)
   Copyright (C) 1995  Ian R. Searle

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
#include "bltin.h"
#include "util.h"
#include "ent.h"
#include "class.h"
#include "rfileio.h"
#include "mathl.h"
#include "list.h"
#include "symbol.h"

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>

#ifdef __riscos
#include "riscos_ut.h"
#endif

/*
 * Header files from the available classes...
 */

#include "ent.h"
#include "btree.h"
#include "bltin.h"
#include "function.h"
#include "mdr.h"
#include "mdrf1.h"
#include "mdrf2.h"
#include "ar.h"
#include "mdc.h"
#include "mdcf1.h"
#include "mdcf2.h"
#include "mdr_mdc.h"
#include "mds.h"
#include "mdsf1.h"
#include "mdr_mds.h"
#include "msr.h"
#include "msrf1.h"
#include "msrf2.h"
#include "msc.h"
#include "mscf1.h"
#include "mscf2.h"

#include "rlab_solver_parameters_names.h"
#include "rlab_macros_code.h"

static OpDef mod_method[NCL][NCL];
static OpDef int_method[NCL];
static OpDef round_method[NCL];
static OpDef symm_method[NCL];
static OpDef abs_method[NCL];
static OpDef maxi_1_method[NCL];
static OpDef mini_1_method[NCL];
static OpDef sin_method[NCL];
static OpDef cos_method[NCL];
static OpDef tan_method[NCL];
static OpDef asin_method[NCL];
static OpDef acos_method[NCL];
static OpDef atan_method[NCL];
static OpDef atan2_method[NCL][NCL];
static OpDef sqrt_method[NCL];
static OpDef log_method[NCL];
static OpDef log10_method[NCL];
static OpDef exp_method[NCL];
static OpDef exp_method_weight[NCL];
static OpDef diag_method[NCL];
static OpDef mnorm_method[NCL];
static OpDef rownorm_method[NCL];
static OpDef pnorm_method[NCL];
// static OpDef eig_s_method[NCL];
// static OpDef eig_g_method[NCL][NCL];

// static OpDef ar_s_method[NCL][NCL];
//static OpDef ar_g_method[NCL][NCL];


static OpDef factor_method[NCL];
static OpDef backsub_ge_method[NCL][NCL];
static OpDef backsub_sy_method[NCL][NCL];
static OpDef backsub_sparse_method[NCL][NCL];
static OpDef svd_method[NCL];
static OpDef chol_method[NCL];
static OpDef solve_method[NCL][NCL];
static OpDef rcond_method[NCL];
static OpDef hess_method[NCL];
static OpDef balance_method[NCL];
static OpDef qr_method[NCL];
static OpDef qrp_method[NCL];
static OpDef schur_method[NCL];
static OpDef sylv_method[NCL][NCL];
static OpDef det_method[NCL];

/* **************************************************************
 * Initialize the built-ins...
 * ************************************************************** */

void
class_bltin2_init (void)
{
  /*
   * mod ()
   */

  mod_method[MATRIX_DENSE_REAL][MATRIX_DENSE_REAL].type = MATRIX_DENSE_REAL;
  mod_method[MATRIX_DENSE_REAL][MATRIX_DENSE_REAL].op = (void *) mdr_Mod_BF;

  mod_method[MATRIX_DENSE_COMPLEX][MATRIX_DENSE_COMPLEX].type =
    MATRIX_DENSE_COMPLEX;
  mod_method[MATRIX_DENSE_COMPLEX][MATRIX_DENSE_COMPLEX].op = (void *) mdc_Mod;

  mod_method[MATRIX_DENSE_REAL][MATRIX_DENSE_COMPLEX].type =
    MATRIX_DENSE_COMPLEX;
  mod_method[MATRIX_DENSE_REAL][MATRIX_DENSE_COMPLEX].op = (void *) mdr_mdc_Mod;

  mod_method[MATRIX_DENSE_COMPLEX][MATRIX_DENSE_REAL].type =
    MATRIX_DENSE_COMPLEX;
  mod_method[MATRIX_DENSE_COMPLEX][MATRIX_DENSE_REAL].op = (void *) mdc_mdr_Mod;

  /*
   * int ()
   */

  int_method[MATRIX_DENSE_REAL].type = MATRIX_DENSE_REAL;
  int_method[MATRIX_DENSE_REAL].op = (void *) mdr_Int_BF;

  int_method[MATRIX_DENSE_COMPLEX].type = MATRIX_DENSE_COMPLEX;
  int_method[MATRIX_DENSE_COMPLEX].op = (void *) mdc_Int_BF;

  int_method[MATRIX_SPARSE_REAL].type = MATRIX_SPARSE_REAL;
  int_method[MATRIX_SPARSE_REAL].op = (void *) msr_Int_BF;

  int_method[MATRIX_SPARSE_COMPLEX].type = MATRIX_SPARSE_COMPLEX;
  int_method[MATRIX_SPARSE_COMPLEX].op = (void *) msc_Int_BF;

  /*
   * round ()
   */

  round_method[MATRIX_DENSE_REAL].type = MATRIX_DENSE_REAL;
  round_method[MATRIX_DENSE_REAL].op = (void *) mdr_Round_BF;

  round_method[MATRIX_DENSE_COMPLEX].type = MATRIX_DENSE_COMPLEX;
  round_method[MATRIX_DENSE_COMPLEX].op = (void *) mdc_Round_BF;

  round_method[MATRIX_SPARSE_REAL].type = MATRIX_SPARSE_REAL;
  round_method[MATRIX_SPARSE_REAL].op = (void *) msr_Round_BF;

  round_method[MATRIX_SPARSE_COMPLEX].type = MATRIX_SPARSE_COMPLEX;
  round_method[MATRIX_SPARSE_COMPLEX].op = (void *) msc_Round_BF;

  /*
   * symm (), always returns a MATRIX_DENSE_REAL..
   */

  symm_method[MATRIX_DENSE_REAL].type = 0;
  symm_method[MATRIX_DENSE_REAL].op = (void *) mdr_IsSymmetric;

  symm_method[MATRIX_DENSE_COMPLEX].type = 0;
  symm_method[MATRIX_DENSE_COMPLEX].op = (void *) mdc_IsSymmetric;

  symm_method[MATRIX_DENSE_STRING].type = 0;
  symm_method[MATRIX_DENSE_STRING].op = (void *) mds_IsSymmetric;

  symm_method[MATRIX_SPARSE_REAL].type = 0;
  symm_method[MATRIX_SPARSE_REAL].op = (void *) msr_IsSymmetric;

  symm_method[MATRIX_SPARSE_COMPLEX].type = 0;
  symm_method[MATRIX_SPARSE_COMPLEX].op = (void *) msc_IsSymmetric;

  /*
   * abs ()
   */

  abs_method[MATRIX_DENSE_REAL].type = MATRIX_DENSE_REAL;
  abs_method[MATRIX_DENSE_REAL].op = (void *) mdr_Abs;

  abs_method[MATRIX_DENSE_COMPLEX].type = MATRIX_DENSE_REAL;
  abs_method[MATRIX_DENSE_COMPLEX].op = (void *) mdc_Abs;

  abs_method[MATRIX_SPARSE_REAL].type = MATRIX_SPARSE_REAL;
  abs_method[MATRIX_SPARSE_REAL].op = (void *) msr_Abs;

  abs_method[MATRIX_SPARSE_COMPLEX].type = MATRIX_SPARSE_REAL;
  abs_method[MATRIX_SPARSE_COMPLEX].op = (void *) msc_Abs;

  /*
   * max ()
   */

  /*
   * min ()
   */

  /*
   * maxi ()
   */

  maxi_1_method[MATRIX_DENSE_REAL].type = MATRIX_DENSE_REAL;
  maxi_1_method[MATRIX_DENSE_REAL].op = (void *) mdr_MinMax1_ValIdx;

  maxi_1_method[MATRIX_DENSE_COMPLEX].type = MATRIX_DENSE_COMPLEX;
  maxi_1_method[MATRIX_DENSE_COMPLEX].op = (void *) mdc_MaxI1;

  /*
   * mini ()
   */

  mini_1_method[MATRIX_DENSE_REAL].type = MATRIX_DENSE_REAL;
  mini_1_method[MATRIX_DENSE_REAL].op = (void *) mdr_MinMax1_ValIdx;

  mini_1_method[MATRIX_DENSE_COMPLEX].type = MATRIX_DENSE_COMPLEX;
  mini_1_method[MATRIX_DENSE_COMPLEX].op = (void *) mdc_MinI1;

  /*
   * sin ()
   */

  sin_method[MATRIX_DENSE_REAL].type = MATRIX_DENSE_REAL;
  sin_method[MATRIX_DENSE_REAL].op = (void *) mdr_Sin;

  sin_method[MATRIX_DENSE_COMPLEX].type = MATRIX_DENSE_COMPLEX;
  sin_method[MATRIX_DENSE_COMPLEX].op = (void *) mdc_Sin;

  /*
   * cos ()
   */

  cos_method[MATRIX_DENSE_REAL].type = MATRIX_DENSE_REAL;
  cos_method[MATRIX_DENSE_REAL].op = (void *) mdr_Cos;

  cos_method[MATRIX_DENSE_COMPLEX].type = MATRIX_DENSE_COMPLEX;
  cos_method[MATRIX_DENSE_COMPLEX].op = (void *) mdc_Cos;

  /*
   * tan ()
   */

  tan_method[MATRIX_DENSE_REAL].type = MATRIX_DENSE_REAL;
  tan_method[MATRIX_DENSE_REAL].op = (void *) mdr_Tan;

  tan_method[MATRIX_DENSE_COMPLEX].type = MATRIX_DENSE_COMPLEX;
  tan_method[MATRIX_DENSE_COMPLEX].op = (void *) mdc_Tan;

  /*
   * asin ()
   */

  asin_method[MATRIX_DENSE_REAL].type = 0;
  asin_method[MATRIX_DENSE_REAL].op = (void *) mdr_ASin;

  asin_method[MATRIX_DENSE_COMPLEX].type = 0;
  asin_method[MATRIX_DENSE_COMPLEX].op = (void *) mdc_ASin;

  /*
   * acos ()
   */

  acos_method[MATRIX_DENSE_REAL].type = 0;
  acos_method[MATRIX_DENSE_REAL].op = (void *) mdr_ACos;

  acos_method[MATRIX_DENSE_COMPLEX].type = 0;
  acos_method[MATRIX_DENSE_COMPLEX].op = (void *) mdc_ACos;

  /*
   * atan ()
   */

  atan_method[MATRIX_DENSE_REAL].type = 0;
  atan_method[MATRIX_DENSE_REAL].op = (void *) mdr_ATan;

  atan_method[MATRIX_DENSE_COMPLEX].type = 0;
  atan_method[MATRIX_DENSE_COMPLEX].op = (void *) mdc_ATan;

  /*
   * atan2 ()
   */

  atan2_method[MATRIX_DENSE_REAL][MATRIX_DENSE_REAL].type = 0;
  atan2_method[MATRIX_DENSE_REAL][MATRIX_DENSE_REAL].op = (void *) mdr_ATan2;
  atan2_method[MATRIX_DENSE_REAL][UNDEF].type = 0;
  atan2_method[MATRIX_DENSE_REAL][UNDEF].op = (void *) mdr_ATan2;

  /*
   * sqrt ()
   */

  sqrt_method[MATRIX_DENSE_REAL].type = 0;
  sqrt_method[MATRIX_DENSE_REAL].op = (void *) mdr_Sqrt;

  sqrt_method[MATRIX_DENSE_COMPLEX].type = 0;
  sqrt_method[MATRIX_DENSE_COMPLEX].op = (void *) mdc_Sqrt;

  /*
   * log ()
   */

  log_method[MATRIX_DENSE_REAL].type = 0;
  log_method[MATRIX_DENSE_REAL].op = (void *) mdr_Log;

  log_method[MATRIX_DENSE_COMPLEX].type = 0;
  log_method[MATRIX_DENSE_COMPLEX].op = (void *) mdc_Log;

  /*
   * log10 ()
   */

  log10_method[MATRIX_DENSE_REAL].type = 0;
  log10_method[MATRIX_DENSE_REAL].op = (void *) mdr_Log10;

  log10_method[MATRIX_DENSE_COMPLEX].type = 0;
  log10_method[MATRIX_DENSE_COMPLEX].op = (void *) mdc_Log10;

  /*
   * exp ()
   */

  exp_method[MATRIX_DENSE_REAL].type = MATRIX_DENSE_REAL;
  exp_method[MATRIX_DENSE_REAL].op = (void *) mdr_Exp;

  exp_method[MATRIX_DENSE_COMPLEX].type = MATRIX_DENSE_COMPLEX;
  exp_method[MATRIX_DENSE_COMPLEX].op = (void *) mdc_Exp;

  exp_method_weight[MATRIX_DENSE_REAL].type = MATRIX_DENSE_REAL;
  exp_method_weight[MATRIX_DENSE_REAL].op = (void *) mdr_Exp_weight;

  /*
   * diag ()
   */

  diag_method[MATRIX_DENSE_REAL].type = MATRIX_DENSE_REAL;
  diag_method[MATRIX_DENSE_REAL].op = (void *) mdr_Diag;

  diag_method[MATRIX_DENSE_COMPLEX].type = MATRIX_DENSE_COMPLEX;
  diag_method[MATRIX_DENSE_COMPLEX].op = (void *) mdc_Diag;

  diag_method[MATRIX_SPARSE_REAL].type = MATRIX_SPARSE_REAL;
  diag_method[MATRIX_SPARSE_REAL].op = (void *) msr_Diag;

  /*
   * rownorm ()
   */
  rownorm_method[MATRIX_DENSE_REAL].type = 0;
  rownorm_method[MATRIX_DENSE_REAL].op = (void *) mdr_RowNorm;

  rownorm_method[MATRIX_DENSE_COMPLEX].type = 0;
  rownorm_method[MATRIX_DENSE_COMPLEX].op = (void *) mdc_RowNorm;

  /*
   * mnorm ()
   */
  mnorm_method[MATRIX_DENSE_REAL].type = 0;
  mnorm_method[MATRIX_DENSE_REAL].op = (void *) mdr_Norm;

  mnorm_method[MATRIX_DENSE_COMPLEX].type = 0;
  mnorm_method[MATRIX_DENSE_COMPLEX].op = (void *) mdc_Norm;

  mnorm_method[MATRIX_SPARSE_REAL].type = 0;
  mnorm_method[MATRIX_SPARSE_REAL].op = (void *) msr_Norm;

  mnorm_method[MATRIX_SPARSE_COMPLEX].type = 0;
  mnorm_method[MATRIX_SPARSE_COMPLEX].op = (void *) msc_Norm;

  /*
   * pnorm ()
   */

  pnorm_method[MATRIX_DENSE_REAL].type = 0;
  pnorm_method[MATRIX_DENSE_REAL].op = (void *) mdr_PNorm;

  pnorm_method[MATRIX_DENSE_COMPLEX].type = 0;
  pnorm_method[MATRIX_DENSE_COMPLEX].op = (void *) mdc_PNorm;

  /*
   * eig ()

  eig_s_method[MATRIX_DENSE_REAL].type = BTREE;
  eig_s_method[MATRIX_DENSE_REAL].op = (void *) mdr_EigS;

  eig_s_method[MATRIX_DENSE_COMPLEX].type = BTREE;
  eig_s_method[MATRIX_DENSE_COMPLEX].op = (void *) mdc_EigS;

  eig_g_method[MATRIX_DENSE_REAL][MATRIX_DENSE_REAL].type = BTREE;
  eig_g_method[MATRIX_DENSE_REAL][MATRIX_DENSE_REAL].op = (void *) mdr_EigG;

  eig_g_method[MATRIX_DENSE_COMPLEX][MATRIX_DENSE_COMPLEX].type = BTREE;
  eig_g_method[MATRIX_DENSE_COMPLEX][MATRIX_DENSE_COMPLEX].op = mdc_EigG;
  */


  /*
  * eigs ()
  */

  //ar_s_method[MATRIX_DENSE_REAL].type = BTREE;
  //ar_s_method[MATRIX_DENSE_REAL].op = (void *) mdr_ar_EigsS;
  //ar_s_method[MATRIX_DENSE_REAL][UNDEF].type = BTREE;
  //ar_s_method[MATRIX_DENSE_REAL][UNDEF].op = (void *) mdr_ar_EigsS;

  //ar_s_method[MATRIX_DENSE_COMPLEX].type = BTREE;
  //ar_s_method[MATRIX_DENSE_COMPLEX].op = (void *) mdc_ar_EigS;

  //ar_g_method[MATRIX_DENSE_REAL][MATRIX_DENSE_REAL].type = BTREE;
  //ar_g_method[MATRIX_DENSE_REAL][MATRIX_DENSE_REAL].op = (void *) mdr_ar_EigG;

  //ar_g_method[MATRIX_DENSE_COMPLEX][MATRIX_DENSE_COMPLEX].type = BTREE;
  //ar_g_method[MATRIX_DENSE_COMPLEX][MATRIX_DENSE_COMPLEX].op = mdc_ar_EigG;

  //
  // factor ()
  //

  factor_method[MATRIX_DENSE_REAL].type = BTREE;
  factor_method[MATRIX_DENSE_REAL].op = (void *) mdr_Factor;

  factor_method[MATRIX_DENSE_COMPLEX].type = BTREE;
  factor_method[MATRIX_DENSE_COMPLEX].op = (void *) mdc_Factor;
  factor_method[MATRIX_SPARSE_REAL].type = BTREE;
  factor_method[MATRIX_SPARSE_REAL].op = (void *) msr_Factor;
  //factor_method[MATRIX_SPARSE_REAL].type = BTREE;
  //factor_method[MATRIX_SPARSE_REAL].op = (void *) umfpack_msr_Factor;

  factor_method[MATRIX_SPARSE_COMPLEX].type = BTREE;
  factor_method[MATRIX_SPARSE_COMPLEX].op = (void *) msc_Factor;

  /*
   * backsub ()
   */

  backsub_ge_method[MATRIX_DENSE_REAL][MATRIX_DENSE_REAL].type =
    MATRIX_DENSE_REAL;
  backsub_ge_method[MATRIX_DENSE_REAL][MATRIX_DENSE_REAL].op =
    (void *) mdr_Backsub_Ge;

  backsub_ge_method[MATRIX_DENSE_COMPLEX][MATRIX_DENSE_COMPLEX].type =
    MATRIX_DENSE_COMPLEX;
  backsub_ge_method[MATRIX_DENSE_COMPLEX][MATRIX_DENSE_COMPLEX].op =
    (void *) mdc_Backsub_Ge;

  backsub_sy_method[MATRIX_DENSE_REAL][MATRIX_DENSE_REAL].type =
    MATRIX_DENSE_REAL;
  backsub_sy_method[MATRIX_DENSE_REAL][MATRIX_DENSE_REAL].op =
    (void *) mdr_Backsub_Sym;

  backsub_sy_method[MATRIX_DENSE_COMPLEX][MATRIX_DENSE_COMPLEX].type =
    MATRIX_DENSE_COMPLEX;
  backsub_sy_method[MATRIX_DENSE_COMPLEX][MATRIX_DENSE_COMPLEX].op =
    (void *) mdc_Backsub_Sym;

  backsub_sparse_method[MATRIX_DENSE_REAL][MATRIX_DENSE_REAL].type =
    MATRIX_DENSE_REAL;
  backsub_sparse_method[MATRIX_DENSE_REAL][MATRIX_DENSE_REAL].op =
    (void *) msr_Backsub;

  backsub_sparse_method[MATRIX_DENSE_REAL][MATRIX_DENSE_COMPLEX].type =
    MATRIX_DENSE_COMPLEX;
  backsub_sparse_method[MATRIX_DENSE_REAL][MATRIX_DENSE_COMPLEX].op =
    (void *) msc_Backsub;

  /*
   * svd ()
   */

  svd_method[MATRIX_DENSE_REAL].type = BTREE;
  svd_method[MATRIX_DENSE_REAL].op = (void *) mdr_Svd_BF;

  svd_method[MATRIX_DENSE_COMPLEX].type = BTREE;
  svd_method[MATRIX_DENSE_COMPLEX].op = (void *) mdc_Svd_BF;

  /*
   * chol ()
   */

  chol_method[MATRIX_DENSE_REAL].type = MATRIX_DENSE_REAL;
  chol_method[MATRIX_DENSE_REAL].op = (void *) mdr_Chol;

  chol_method[MATRIX_DENSE_COMPLEX].type = MATRIX_DENSE_COMPLEX;
  chol_method[MATRIX_DENSE_COMPLEX].op = (void *) mdc_Chol;

  /*
   * solve ()
   */

  solve_method[MATRIX_DENSE_REAL][MATRIX_DENSE_REAL].type = MATRIX_DENSE_REAL;
  solve_method[MATRIX_DENSE_REAL][MATRIX_DENSE_REAL].op = (void *) mdr_Solve;

  solve_method[MATRIX_DENSE_REAL][MATRIX_DENSE_COMPLEX].type =
    MATRIX_DENSE_COMPLEX;
  solve_method[MATRIX_DENSE_REAL][MATRIX_DENSE_COMPLEX].op =
    (void *) mdr_mdc_Solve;

  solve_method[MATRIX_DENSE_COMPLEX][MATRIX_DENSE_REAL].type =
    MATRIX_DENSE_COMPLEX;
  solve_method[MATRIX_DENSE_COMPLEX][MATRIX_DENSE_REAL].op =
    (void *) mdc_mdr_Solve;

  solve_method[MATRIX_DENSE_COMPLEX][MATRIX_DENSE_COMPLEX].type =
    MATRIX_DENSE_COMPLEX;
  solve_method[MATRIX_DENSE_COMPLEX][MATRIX_DENSE_COMPLEX].op =
    (void *) mdc_Solve;

  //
  // Solve for sparse matrices
  //
  // real
  solve_method[MATRIX_SPARSE_REAL][MATRIX_DENSE_REAL].type = MATRIX_DENSE_REAL;
  solve_method[MATRIX_SPARSE_REAL][MATRIX_DENSE_REAL].op = (void *) msr_Solve;
  // complex CC
  solve_method[MATRIX_SPARSE_COMPLEX][MATRIX_DENSE_COMPLEX].type =
    MATRIX_DENSE_COMPLEX;
  solve_method[MATRIX_SPARSE_COMPLEX][MATRIX_DENSE_COMPLEX].op =
    (void *) msc_SolveCC;
  // complex CR
  solve_method[MATRIX_SPARSE_COMPLEX][MATRIX_DENSE_REAL].type =
    MATRIX_DENSE_COMPLEX;
  solve_method[MATRIX_SPARSE_COMPLEX][MATRIX_DENSE_REAL].op =
    (void *) msc_SolveCR;
  // complex RC
  solve_method[MATRIX_SPARSE_REAL][MATRIX_DENSE_COMPLEX].type =
    MATRIX_DENSE_COMPLEX;
  solve_method[MATRIX_SPARSE_REAL][MATRIX_DENSE_COMPLEX].op =
    (void *) umfpack_msc_SolveRC;

  /*
   * rcond ()
   */

  rcond_method[MATRIX_DENSE_REAL].type = 0;
  rcond_method[MATRIX_DENSE_REAL].op = (void *) mdr_Rcond;

  rcond_method[MATRIX_DENSE_COMPLEX].type = 0;
  rcond_method[MATRIX_DENSE_COMPLEX].op = (void *) mdc_Rcond;

  /*
   * hess ()
   */

  hess_method[MATRIX_DENSE_REAL].type = BTREE;
  hess_method[MATRIX_DENSE_REAL].op = (void *) mdr_Hess;

  hess_method[MATRIX_DENSE_COMPLEX].type = BTREE;
  hess_method[MATRIX_DENSE_COMPLEX].op = (void *) mdc_Hess;

  /*
   * balance ()
   */

  balance_method[MATRIX_DENSE_REAL].type = BTREE;
  balance_method[MATRIX_DENSE_REAL].op = (void *) mdr_Balance;

  balance_method[MATRIX_DENSE_COMPLEX].type = BTREE;
  balance_method[MATRIX_DENSE_COMPLEX].op = (void *) mdc_Balance;

  /*
   * qr ()
   */

  qr_method[MATRIX_DENSE_REAL].type = BTREE;
  qr_method[MATRIX_DENSE_REAL].op = (void *) mdr_QR;

  qr_method[MATRIX_DENSE_COMPLEX].type = BTREE;
  qr_method[MATRIX_DENSE_COMPLEX].op = (void *) mdc_QR;

  qrp_method[MATRIX_DENSE_REAL].type = BTREE;
  qrp_method[MATRIX_DENSE_REAL].op = (void *) mdr_QRP;

  qrp_method[MATRIX_DENSE_COMPLEX].type = BTREE;
  qrp_method[MATRIX_DENSE_COMPLEX].op = (void *) mdc_QRP;

  /*
   * schur ()
   */

  schur_method[MATRIX_DENSE_REAL].type = BTREE;
  schur_method[MATRIX_DENSE_REAL].op = (void *) mdr_Schur;

  schur_method[MATRIX_DENSE_COMPLEX].type = BTREE;
  schur_method[MATRIX_DENSE_COMPLEX].op = (void *) mdc_Schur;

  /*
   * sylv ()
   */

  sylv_method[MATRIX_DENSE_REAL][MATRIX_DENSE_REAL].type = MATRIX_DENSE_REAL;
  sylv_method[MATRIX_DENSE_REAL][MATRIX_DENSE_REAL].op = (void *) mdr_Sylv;

  sylv_method[MATRIX_DENSE_COMPLEX][MATRIX_DENSE_COMPLEX].type =
    MATRIX_DENSE_COMPLEX;
  sylv_method[MATRIX_DENSE_COMPLEX][MATRIX_DENSE_COMPLEX].op =
    (void *) mdc_Sylv;

  //
  // det ()
  //

  det_method[MATRIX_DENSE_REAL].type = MATRIX_DENSE_REAL;
  det_method[MATRIX_DENSE_REAL].op = (void *) mdr_Det;
  det_method[MATRIX_SPARSE_REAL].type = MATRIX_DENSE_REAL;
  det_method[MATRIX_SPARSE_REAL].op = (void *) umfpack_msr_Det;
  det_method[MATRIX_DENSE_COMPLEX].type = MATRIX_DENSE_COMPLEX;
  det_method[MATRIX_DENSE_COMPLEX].op = (void *) mdc_Det;
  det_method[MATRIX_SPARSE_COMPLEX].type = MATRIX_DENSE_COMPLEX;
  det_method[MATRIX_SPARSE_COMPLEX].op = (void *) umfpack_msc_Det;
}

/* **************************************************************
 * Mod function...
 * ************************************************************** */

Ent *
Mod (int nargs, Datum args[])
{
  int type1, type2, rtype;
  void *rval;
  void *(*vfptr) ();
  Ent *rent, *e1, *e2;

  if (nargs != 2)
  {
    rerror ("mod: 2 arguments required");
  }

  e1 = bltin_get_ent (args[0]);
  e2 = bltin_get_ent (args[1]);
  type1 = ent_type (e1);
  type2 = ent_type (e2);

  rtype = mod_method[type1][type2].type;
  vfptr = (VFPTR) mod_method[type1][type2].op;

  if (vfptr == 0)
  {
    fprintf (stderr, "Entity types: %s and %s\n", etd (e1), etd (e2));
    rerror ("mod() operation not supported");
  }

  rval = (*vfptr) (ent_data (e1), ent_data (e2));

  rent = ent_Create ();
  ent_data (rent) = rval;
  ent_type (rent) = rtype;

  ent_Clean (e1);
  ent_Clean (e2);
  return (rent);
}

/* **************************************************************
 * Int function...
 * ************************************************************** */

Ent *
Int (int nargs, Datum args[])
{
  int type=UNDEF, rtype;
  void *rval=0;
  void *(*vfptr) ();
  Ent *rent, *e;

  if (nargs != 1)
  {
    rerror ("int: 1 argument allowed");
  }

  e = bltin_get_ent (args[0]);
  type = ent_type (e);

  rtype = int_method[type].type;
  vfptr = (VFPTR) int_method[type].op;

  if (vfptr == 0)
  {
    fprintf (stderr, "Entity type: %s\n", etd (e));
    rerror ("int() operation not supported");
  }

  rval = (*vfptr) (ent_data (e));

  rent = ent_Create ();
  ent_data (rent) = rval;
  ent_type (rent) = rtype;

  ent_Clean (e);
  return (rent);
}

/* **************************************************************
 * Ceil function...
 * ************************************************************** */
#undef THIS_SOLVER
#define THIS_SOLVER "ceil"
Ent * Ceil (int nargs, Datum args[])
{
  int rtype=UNDEF;
  MDR *mesh=0, *bin=0, *clamp=0;
  MD  *off=0;
  void *rval=0;
  Ent *rent, *e1=0, *e2=0;
  FILE *rlab_stderr = (!RLAB_STDERR_DS) ? stderr : RLAB_STDERR_DS;
  double base=-1.0;
  ListNode *node=0;

  /* Check n_args */
  if (nargs != 1 && nargs != 2)
    rerror (THIS_SOLVER ": " RLAB_ERROR_ONE_OR_TWO_ARG_REQUIRED "\n");

  e1 = bltin_get_ent (args[0]);
  if (   ent_type (e1) != MATRIX_DENSE_REAL && ent_type (e1) != MATRIX_DENSE_COMPLEX
      && ent_type (e1) != MATRIX_SPARSE_REAL && ent_type (e1) != MATRIX_SPARSE_COMPLEX  )
  {
    fprintf (rlab_stderr, "%s : %s\n", THIS_SOLVER, RLAB_ERROR_ARG1_NUM);
    goto _exit_ceil;
  }
  if (nargs == 2)
  {
    e2 = bltin_get_ent (args[1]);
    if (ent_type (e2) == BTREE)
    {
      node = btree_FindNode (ent_data (e2), RLAB_NAME_MATH_OFFSET);
      if (node)
      {
        if (ent_type (var_ent (node)) != MATRIX_DENSE_REAL)
        {
          fprintf (rlab_stderr, THIS_SOLVER ": list entry " RLAB_NAME_MATH_OFFSET " of wrong type\n");
          goto _exit_ceil;
        }
        off = ent_data (var_ent (node));
        if (SIZE(off)<1)
          off = 0;
      }

      node = btree_FindNode (ent_data (e2), RLAB_NAME_MATH_BIN);
      if (node)
      {
        if (ent_type (var_ent (node)) != MATRIX_DENSE_REAL)
        {
          fprintf (rlab_stderr, THIS_SOLVER ": list entry " RLAB_NAME_MATH_BIN " of wrong type\n");
          goto _exit_ceil;
        }
        bin = ent_data (var_ent (node));
        if (SIZE(bin)<1)
          bin = 0;
      }

      node = btree_FindNode (ent_data (e2), RLAB_NAME_MATH_BASE);
      if (node)
      {
        if (ent_type (var_ent (node)) != MATRIX_DENSE_REAL)
        {
          fprintf (rlab_stderr, THIS_SOLVER ": list entry " RLAB_NAME_MATH_BASE " of wrong type\n");
          goto _exit_ceil;
        }
        base = class_double(var_ent (node));
        if (base<=0)
        {
          fprintf (rlab_stderr, THIS_SOLVER ": list entry " RLAB_NAME_MATH_BASE " of wrong type\n");
          goto _exit_ceil;
        }
      }

      node = btree_FindNode (ent_data (e2), RLAB_NAME_MATH_MESH);
      if (node)
      {
        if (ent_type (var_ent (node)) != MATRIX_DENSE_REAL)
        {
          fprintf (rlab_stderr, THIS_SOLVER ": list entry " RLAB_NAME_MATH_MESH " of wrong type\n");
          goto _exit_ceil;
        }
        mesh = ent_data (var_ent (node));
        if (SIZE(mesh)<1)
          mesh = 0;
      }

      node = btree_FindNode (ent_data (e2), RLAB_NAME_MATH_CLAMP);
      if (node)
      {
        if (ent_type (var_ent (node)) != MATRIX_DENSE_REAL)
        {
          fprintf (rlab_stderr, THIS_SOLVER ": list entry " RLAB_NAME_MATH_CLAMP " of wrong type\n");
          goto _exit_ceil;
        }
        clamp = ent_data (var_ent (node));
        if (SIZE(clamp)!=2)
          clamp = 0;
      }
    }
  }

  switch(ent_type(e1))
  {
    case MATRIX_DENSE_REAL:
      rtype = MATRIX_DENSE_REAL;
      if (mesh)
        rval = mdr_MeshCeilFlooRound_BF (1, ent_data(e1), mesh);
      else if (base>0)
        rval = mdr_LogCeilFlooRound_BF (1, ent_data(e1), base);
      else
        rval = mdr_CeilFlooRound_BF (1, ent_data(e1), bin, (MD *) off, clamp);
      break;

    case MATRIX_DENSE_COMPLEX:
      rtype = MATRIX_DENSE_COMPLEX;
      rval = mdc_CeilFlooRound_BF (1, ent_data(e1), bin, (MD *) off);
      break;

    default:
      break;
  }

_exit_ceil:

  ent_Clean (e1);
  ent_Clean (e2);

  rent = ent_Create ();
  ent_data (rent) = rval;
  ent_type (rent) = rtype;
  return (rent);
}

/* **************************************************************
 * Floor function...
 * ************************************************************** */
#undef THIS_SOLVER
#define THIS_SOLVER "floor"
Ent * Floor (int nargs, Datum args[])
{
  int rtype=UNDEF;
  MDR *mesh=0, *bin=0, *clamp=0;
  MD  *off=0;
  void *rval=0;
  Ent *rent, *e1=0, *e2=0;
  FILE *rlab_stderr = (!RLAB_STDERR_DS) ? stderr : RLAB_STDERR_DS;
  double base=-1.0;
  ListNode *node=0;

  /* Check n_args */
  if (nargs != 1 && nargs != 2)
    rerror (THIS_SOLVER ": " RLAB_ERROR_ONE_OR_TWO_ARG_REQUIRED "\n");

  e1 = bltin_get_ent (args[0]);
  if (   ent_type (e1) != MATRIX_DENSE_REAL && ent_type (e1) != MATRIX_DENSE_COMPLEX
         && ent_type (e1) != MATRIX_SPARSE_REAL && ent_type (e1) != MATRIX_SPARSE_COMPLEX  )
  {
    fprintf (rlab_stderr, THIS_SOLVER ": " RLAB_ERROR_ARG1_MDS_MATRIX);
    goto _exit_ceil;
  }
  if (nargs == 2)
  {
    e2 = bltin_get_ent (args[1]);
    if (ent_type (e2) == BTREE)
    {
      node = btree_FindNode (ent_data (e2), RLAB_NAME_MATH_OFFSET);
      if (node)
      {
        if (ent_type (var_ent (node)) != MATRIX_DENSE_REAL)
        {
          fprintf (rlab_stderr, THIS_SOLVER ": list entry " RLAB_NAME_MATH_OFFSET " of wrong type\n");
          goto _exit_ceil;
        }
        off = ent_data (var_ent (node));
        if (SIZE(off)<1)
          off = 0;
      }

      node = btree_FindNode (ent_data (e2), RLAB_NAME_MATH_BIN);
      if (node)
      {
        if (ent_type (var_ent (node)) != MATRIX_DENSE_REAL)
        {
          fprintf (rlab_stderr, THIS_SOLVER ": list entry " RLAB_NAME_MATH_BIN " of wrong type\n");
          goto _exit_ceil;
        }
        bin = ent_data (var_ent (node));
        if (SIZE(bin)<1)
          bin = 0;
      }

      node = btree_FindNode (ent_data (e2), RLAB_NAME_MATH_BASE);
      if (node)
      {
        if (ent_type (var_ent (node)) != MATRIX_DENSE_REAL)
        {
          fprintf (rlab_stderr, THIS_SOLVER ": list entry " RLAB_NAME_MATH_BASE " of wrong type\n");
          goto _exit_ceil;
        }
        base = class_double(var_ent (node));
        if (base<=0)
        {
          fprintf (rlab_stderr, THIS_SOLVER ": list entry " RLAB_NAME_MATH_BASE " of wrong type\n");
          goto _exit_ceil;
        }
      }

      node = btree_FindNode (ent_data (e2), RLAB_NAME_MATH_MESH);
      if (node)
      {
        if (ent_type (var_ent (node)) != MATRIX_DENSE_REAL)
        {
          fprintf (rlab_stderr, THIS_SOLVER ": list entry " RLAB_NAME_MATH_MESH " of wrong type\n");
          goto _exit_ceil;
        }
        mesh = ent_data (var_ent (node));
        if (SIZE(mesh)<1)
          mesh = 0;
        if (SIZE(mesh)<3)
        {
          fprintf (rlab_stderr, THIS_SOLVER ": list entry " RLAB_NAME_MATH_MESH " of wrong size\n");
          goto _exit_ceil;
        }
      }

      node = btree_FindNode (ent_data (e2), RLAB_NAME_MATH_CLAMP);
      if (node)
      {
        if (ent_type (var_ent (node)) != MATRIX_DENSE_REAL)
        {
          fprintf (rlab_stderr, THIS_SOLVER ": list entry " RLAB_NAME_MATH_CLAMP " of wrong type\n");
          goto _exit_ceil;
        }
        clamp = ent_data (var_ent (node));
        if (SIZE(clamp)!=2)
          clamp = 0;
      }
    }
  }

  switch(ent_type(e1))
  {
    case MATRIX_DENSE_REAL:
      rtype = MATRIX_DENSE_REAL;
      if (mesh)
      {
        rval = mdr_MeshCeilFlooRound_BF (0, ent_data(e1), mesh);
      }
      else if (base>0)
      {
        rval = mdr_LogCeilFlooRound_BF (0, ent_data(e1), base);
      }
      else
      {
        rval = mdr_CeilFlooRound_BF (0, ent_data(e1), bin, (MD *) off, clamp);
      }
      break;

    case MATRIX_DENSE_COMPLEX:
      rtype = MATRIX_DENSE_COMPLEX;
      rval = mdc_CeilFlooRound_BF (0, ent_data(e1), bin, (MD *) off);
      break;

    default:
      break;
  }

_exit_ceil:

  ent_Clean (e1);
  ent_Clean (e2);

  rent = ent_Create ();
  ent_data (rent) = rval;
  ent_type (rent) = rtype;
  return (rent);
}

/* **************************************************************
 * Round (rint) function...
 * ************************************************************** */
#undef THIS_SOLVER
#define THIS_SOLVER "round"
Ent * Round (int nargs, Datum args[])
{
  int rtype=UNDEF;
  MDR *mesh=0, *bin=0, *clamp=0;
  MD  *off=0;
  void *rval=0;
  Ent *rent, *e1=0, *e2=0;
  FILE *rlab_stderr = (!RLAB_STDERR_DS) ? stderr : RLAB_STDERR_DS;
  double base=-1.0;
  ListNode *node=0;

  /* Check n_args */
  if (nargs != 1 && nargs != 2)
    rerror (THIS_SOLVER ": " RLAB_ERROR_ONE_OR_TWO_ARG_REQUIRED "\n");

  e1 = bltin_get_ent (args[0]);
  if (   ent_type (e1) != MATRIX_DENSE_REAL && ent_type (e1) != MATRIX_DENSE_COMPLEX
         && ent_type (e1) != MATRIX_SPARSE_REAL && ent_type (e1) != MATRIX_SPARSE_COMPLEX  )
  {
    fprintf (rlab_stderr, THIS_SOLVER ": " RLAB_ERROR_ARG1_MDS_MATRIX);
    goto _exit_ceil;
  }
  if (nargs == 2)
  {
    e2 = bltin_get_ent (args[1]);
    if (ent_type (e2) == BTREE)
    {
      node = btree_FindNode (ent_data (e2), RLAB_NAME_MATH_OFFSET);
      if (node)
      {
        if (ent_type (var_ent (node)) != MATRIX_DENSE_REAL)
        {
          fprintf (rlab_stderr, THIS_SOLVER ": list entry " RLAB_NAME_MATH_OFFSET " of wrong type\n");
          goto _exit_ceil;
        }
        off = ent_data (var_ent (node));
        if (SIZE(off)<1)
          off = 0;
      }

      node = btree_FindNode (ent_data (e2), RLAB_NAME_MATH_BIN);
      if (node)
      {
        if (ent_type (var_ent (node)) != MATRIX_DENSE_REAL)
        {
          fprintf (rlab_stderr, THIS_SOLVER ": list entry " RLAB_NAME_MATH_BIN " of wrong type\n");
          goto _exit_ceil;
        }
        bin = ent_data (var_ent (node));
        if (SIZE(bin)<1)
          bin = 0;
      }

      node = btree_FindNode (ent_data (e2), RLAB_NAME_MATH_BASE);
      if (node)
      {
        if (ent_type (var_ent (node)) != MATRIX_DENSE_REAL)
        {
          fprintf (rlab_stderr, THIS_SOLVER ": list entry " RLAB_NAME_MATH_BASE " of wrong type\n");
          goto _exit_ceil;
        }
        base = class_double(var_ent (node));
        if (base<=0)
        {
          fprintf (rlab_stderr, THIS_SOLVER ": list entry " RLAB_NAME_MATH_BASE " of wrong type\n");
          goto _exit_ceil;
        }
      }

      node = btree_FindNode (ent_data (e2), RLAB_NAME_MATH_MESH);
      if (node)
      {
        if (ent_type (var_ent (node)) != MATRIX_DENSE_REAL)
        {
          fprintf (rlab_stderr, THIS_SOLVER ": list entry " RLAB_NAME_MATH_MESH " of wrong type\n");
          goto _exit_ceil;
        }
        mesh = ent_data (var_ent (node));
        if (SIZE(mesh)<1)
          mesh = 0;
      }

      node = btree_FindNode (ent_data (e2), RLAB_NAME_MATH_CLAMP);
      if (node)
      {
        if (ent_type (var_ent (node)) != MATRIX_DENSE_REAL)
        {
          fprintf (rlab_stderr, THIS_SOLVER ": list entry " RLAB_NAME_MATH_CLAMP " of wrong type\n");
          goto _exit_ceil;
        }
        clamp = ent_data (var_ent (node));
        if (SIZE(clamp)!=2)
          clamp = 0;
      }
    }
  }

  switch(ent_type(e1))
  {
    case MATRIX_DENSE_REAL:
      rtype = MATRIX_DENSE_REAL;
      if (mesh)
      {
        rval = mdr_MeshCeilFlooRound_BF (-1, ent_data(e1), mesh);
      }
      else if (base>0)
      {
        rval = mdr_LogCeilFlooRound_BF (-1, ent_data(e1), base);
      }
      else
      {
        rval = mdr_CeilFlooRound_BF (-1, ent_data(e1), bin, (MD *) off, clamp);
      }
      break;

    case MATRIX_DENSE_COMPLEX:
      rtype = MATRIX_DENSE_COMPLEX;
      rval = mdc_CeilFlooRound_BF (0, ent_data(e1), bin, (MD *) off);
      break;

    default:
      break;
  }

_exit_ceil:

  ent_Clean (e1);
  ent_Clean (e2);

  rent = ent_Create ();
  ent_data (rent) = rval;
  ent_type (rent) = rtype;
  return (rent);
}

// {
//   int type=UNDEF, rtype;
//   void *rval;
//   void *(*vfptr) ();
//   Ent *rent, *e;
//
//   if (nargs != 1)
//   {
//     rerror ("round: 1 argument allowed");
//   }
//
//   e = bltin_get_ent (args[0]);
//   type = ent_type (e);
//
//   rtype = round_method[type].type;
//   vfptr = (VFPTR) round_method[type].op;
//
//   if (vfptr == 0)
//   {
//     fprintf (stderr, "Entity type: %s\n", etd (e));
//     rerror ("round() operation not supported");
//   }
//
//   rval = (*vfptr) (ent_data (e));
//
//   rent = ent_Create ();
//   ent_data (rent) = rval;
//   ent_type (rent) = rtype;
//
//   ent_Clean (e);
//   return (rent);
// }

/* **************************************************************
 * Symm function...
 * ************************************************************** */

Ent *
IsSymm (int nargs, Datum args[])
{
  int type;
//   int rtype;
  int rval;
  void *(*vfptr) ();
  Ent *rent, *e;

  if (nargs != 1)
  {
    rerror ("symm: 1 argument allowed");
  }

  e = bltin_get_ent (args[0]);
  type = ent_type (e);

//   rtype = symm_method[type].type;
  vfptr = (VFPTR) symm_method[type].op;

  if (vfptr == 0)
  {
    fprintf (stderr, "Entity type: %s\n", etd (e));
    rerror ("symm() operation not supported");
  }

  rval = (long int) (*vfptr) (ent_data (e));

  rent = ent_Create ();
  ent_data (rent) = mdr_CreateScalar ((double) rval);
  ent_type (rent) = MATRIX_DENSE_REAL;

  ent_Clean (e);
  return (rent);
}

/* **************************************************************
 * Abs function...
 * ************************************************************** */

Ent *
Abs (int nargs, Datum args[])
{
  int type;
  int rtype=UNDEF;
  void *rval;
  void *(*vfptr) ();
  Ent *e=0;

  if (nargs != 1)
  {
    rerror ("abs: 1 argument allowed");
  }

  e = bltin_get_ent (args[0]);
  type = ent_type (e);

  rtype = abs_method[type].type;
  vfptr = (VFPTR) abs_method[type].op;

  if (vfptr == 0)
  {
    fprintf (stderr, "Entity type: %s\n", etd (e));
    rerror ("abs() operation not supported");
  }

  rval = (*vfptr) (ent_data (e));

  ent_Clean (e);

  return ent_Assign_Rlab_Rtype(rval, rtype);
}

/* **************************************************************
 * Min function...
 * ************************************************************** */
#undef THIS_SOLVER
#undef THIS_SOLVER_ENT
#define THIS_SOLVER "min"
#define THIS_SOLVER_ENT Min
Ent * THIS_SOLVER_ENT (int nargs, Datum args[])
{
  MDR *rval=0;
  Ent **e=0;
  int loc_nargs=nargs, ignore_inf=0, flat=0, row_dom=0, i;
  double *reject_below=NULL, *reject_above=NULL;
  double rej_bel,rej_ab;
  ListNode *node=0;

  if (nargs < 1)
    goto exit;

  // assign all arguments:
  e = (Ent **) GC_MALLOC(nargs * sizeof(Ent *));
  for (i=0; i<nargs; i++)
  {
    e[i] = bltin_get_ent (args[i]);
  }

  // is the last argument list of options, if so load it
  if (ent_type(e[loc_nargs-1]) == BTREE)
  {
    // reject above?
    node = btree_FindNode (ent_data (e[loc_nargs-1]), RLAB_NAME_STAT_REJECT_ABOVE);
    if (node != 0)
    {
      if (ent_type(var_ent (node))==MATRIX_DENSE_REAL)
      {
        rej_ab = class_double(var_ent (node));
        reject_above = (double *) &rej_ab;
      }
    }
    // reject below?
    node = btree_FindNode (ent_data (e[loc_nargs-1]), RLAB_NAME_STAT_REJECT_BELOW);
    if (node != 0)
    {
      if (ent_type(var_ent (node))==MATRIX_DENSE_REAL)
      {
        rej_bel = class_double(var_ent (node));
        reject_below = (double *) &rej_bel;
      }
    }
    // ignore infs
    node = btree_FindNode (ent_data (e[loc_nargs-1]), RLAB_NAME_STAT_IGNOREINF);
    if (node != 0)
    {
      if (ent_type(var_ent (node))==MATRIX_DENSE_REAL)
      {
        ignore_inf = class_int(var_ent (node));
        if (ignore_inf)
          ignore_inf = 1;
      }
    }
    // flatten a matrix?
    node = btree_FindNode (ent_data (e[loc_nargs-1]), RLAB_NAME_STAT_FLAT);
    if (node != 0)
    {
      if (ent_type(var_ent (node))==MATRIX_DENSE_REAL)
      {
        flat = class_int(var_ent (node));
        if (flat)
          flat = 1;
      }
    }
    else
    {
      // row dominant ? mutually exclusive from 'flat'
      node = btree_FindNode (ent_data (e[loc_nargs-1]), RLAB_NAME_STAT_ROWDOM);
      if (node != 0)
      {
        if (ent_type(var_ent (node))==MATRIX_DENSE_REAL)
        {
          row_dom = class_int(var_ent (node));
          if (row_dom)
            row_dom = 1;
        }
      }
    }

    // argument less one
    loc_nargs--;
  }

  if (loc_nargs == 1)
  {
    if (ent_type(e[0]) == MATRIX_DENSE_REAL)
    {
      rval = mdr_MinMax1_ValIdx(ent_data (e[0]), -2, ignore_inf, flat, row_dom,
                                reject_above, reject_below);
    }
    else
    {
      fprintf (stderr, THIS_SOLVER  ": Entity type: %s\n", etd (e[0]));
      rerror (THIS_SOLVER ": operation not supported");
    }
  }
  else if (loc_nargs == 2)
  {
    if ((ent_type(e[0]) == MATRIX_DENSE_REAL) && (ent_type(e[1]) == MATRIX_DENSE_REAL))
    {
      rval = mdr_Min2 (ent_data(e[0]), ent_data(e[1]), ignore_inf);
    }
    else
    {
      fprintf (stderr, THIS_SOLVER  ": Entity types: %s and %s\n", etd (e[0]), etd (e[1]));
      rerror (THIS_SOLVER ": operation not supported");
    }
  }
  else
  {
    rerror (THIS_SOLVER ": operation not supported");
  }

exit:

    if (e)
    {
      for (i=0; i<nargs; i++)
        ent_Clean(e[i]);
      GC_FREE(e);
    }

    return ent_Assign_Rlab_MDR(rval);
}

/* **************************************************************
 * Max function...
 * ************************************************************** */
#undef THIS_SOLVER
#undef THIS_SOLVER_ENT
#define THIS_SOLVER "max"
#define THIS_SOLVER_ENT Max
Ent * THIS_SOLVER_ENT (int nargs, Datum args[])
{
  MDR *rval=0;
  Ent **e=0;
  int loc_nargs=nargs, ignore_inf=0, flat=0, row_dom=0, i;
  double *reject_below=NULL, *reject_above=NULL;
  double rej_bel,rej_ab;
  ListNode *node=0;

  if (nargs < 1)
    goto exit;

  // assign all arguments:
  e = (Ent **) GC_MALLOC(nargs * sizeof(Ent *));
  for (i=0; i<nargs; i++)
  {
    e[i] = bltin_get_ent (args[i]);
  }

  // is the last argument list of options, if so load it
  if (ent_type(e[loc_nargs-1]) == BTREE)
  {
    // reject above?
    node = btree_FindNode (ent_data (e[loc_nargs-1]), RLAB_NAME_STAT_REJECT_ABOVE);
    if (node != 0)
    {
      if (ent_type(var_ent (node))==MATRIX_DENSE_REAL)
      {
        rej_ab = class_double(var_ent (node));
        reject_above = (double *) &rej_ab;
      }
    }
    // reject below?
    node = btree_FindNode (ent_data (e[loc_nargs-1]), RLAB_NAME_STAT_REJECT_BELOW);
    if (node != 0)
    {
      if (ent_type(var_ent (node))==MATRIX_DENSE_REAL)
      {
        rej_bel = class_double(var_ent (node));
        reject_below = (double *) &rej_bel;
      }
    }
    // ignore infs
    node = btree_FindNode (ent_data (e[loc_nargs-1]), RLAB_NAME_STAT_IGNOREINF);
    if (node != 0)
    {
      if (ent_type(var_ent (node))==MATRIX_DENSE_REAL)
      {
        ignore_inf = class_int(var_ent (node));
        if (ignore_inf)
          ignore_inf = 1;
      }
    }
    // flatten a matrix?
    node = btree_FindNode (ent_data (e[loc_nargs-1]), RLAB_NAME_STAT_FLAT);
    if (node != 0)
    {
      if (ent_type(var_ent (node))==MATRIX_DENSE_REAL)
      {
        flat = class_int(var_ent (node));
        if (flat)
          flat = 1;
      }
    }
    else
    {
      // row dominant ? mutually exclusive from 'flat'
      node = btree_FindNode (ent_data (e[loc_nargs-1]), RLAB_NAME_STAT_ROWDOM);
      if (node != 0)
      {
        if (ent_type(var_ent (node))==MATRIX_DENSE_REAL)
        {
          row_dom = class_int(var_ent (node));
          if (row_dom)
            row_dom = 1;
        }
      }
    }

    // argument less one
    loc_nargs--;
  }

  if (loc_nargs == 1)
  {
    if (ent_type(e[0]) == MATRIX_DENSE_REAL)
      rval = mdr_MinMax1_ValIdx(ent_data (e[0]), 2, ignore_inf, flat, row_dom,
                                reject_above, reject_below);
    else
    {
      fprintf (stderr, THIS_SOLVER  ": Entity type: %s\n", etd (e[0]));
      rerror (THIS_SOLVER ": operation not supported");
    }
  }
  else if (loc_nargs == 2)
  {
    if ((ent_type(e[0]) == MATRIX_DENSE_REAL) && (ent_type(e[1]) == MATRIX_DENSE_REAL))
    {
      rval = mdr_Max2 (ent_data(e[0]), ent_data(e[1]), ignore_inf);
    }
    else
    {
      fprintf (stderr, THIS_SOLVER  ": Entity types: %s and %s\n", etd (e[0]), etd (e[1]));
      rerror (THIS_SOLVER ": operation not supported");
    }
  }
  else
  {
    rerror (THIS_SOLVER ": operation not supported");
  }

exit:

  if (e)
  {
    for (i=0; i<nargs; i++)
      ent_Clean(e[i]);
    GC_FREE(e);
  }

  return ent_Assign_Rlab_MDR(rval);
}

/* **************************************************************
 * MinMax function...
 * ************************************************************** */
#undef THIS_SOLVER
#undef THIS_SOLVER_ENT
#define THIS_SOLVER "minmax"
#define THIS_SOLVER_ENT MinMax
Ent * THIS_SOLVER_ENT (int nargs, Datum args[])
{
  MDR *rval=0;
  Ent **e=0;
  int loc_nargs=nargs, ignore_inf=0, flat=0, row_dom=0, i;
  double *reject_below=NULL, *reject_above=NULL;
  double rej_bel,rej_ab;
  ListNode *node=0;

  if (nargs < 1)
    goto exit;

  // assign all arguments:
  e = (Ent **) GC_MALLOC(nargs * sizeof(Ent *));
  for (i=0; i<nargs; i++)
  {
    e[i] = bltin_get_ent (args[i]);
  }

  // is the last argument list of options, if so load it
  if (ent_type(e[loc_nargs-1]) == BTREE)
  {
    // reject above?
    node = btree_FindNode (ent_data (e[loc_nargs-1]), RLAB_NAME_STAT_REJECT_ABOVE);
    if (node != 0)
    {
      if (ent_type(var_ent (node))==MATRIX_DENSE_REAL)
      {
        rej_ab = class_double(var_ent (node));
        reject_above = (double *) &rej_ab;
      }
    }
    // reject below?
    node = btree_FindNode (ent_data (e[loc_nargs-1]), RLAB_NAME_STAT_REJECT_BELOW);
    if (node != 0)
    {
      if (ent_type(var_ent (node))==MATRIX_DENSE_REAL)
      {
        rej_bel = class_double(var_ent (node));
        reject_below = (double *) &rej_bel;
      }
    }
    // ignore infs
    node = btree_FindNode (ent_data (e[loc_nargs-1]), RLAB_NAME_STAT_IGNOREINF);
    if (node != 0)
    {
      if (ent_type(var_ent (node))==MATRIX_DENSE_REAL)
      {
        ignore_inf = class_int(var_ent (node));
        if (ignore_inf)
          ignore_inf = 1;
      }
    }
    // flatten a matrix?
    node = btree_FindNode (ent_data (e[loc_nargs-1]), RLAB_NAME_STAT_FLAT);
    if (node != 0)
    {
      if (ent_type(var_ent (node))==MATRIX_DENSE_REAL)
      {
        flat = class_int(var_ent (node));
        if (flat)
          flat = 1;
      }
    }
    else
    {
      // row dominant ? mutually exclusive from 'flat'
      node = btree_FindNode (ent_data (e[loc_nargs-1]), RLAB_NAME_STAT_ROWDOM);
      if (node != 0)
      {
        if (ent_type(var_ent (node))==MATRIX_DENSE_REAL)
        {
          row_dom = class_int(var_ent (node));
          if (row_dom)
            row_dom = 1;
        }
      }
    }

    // argument less one
    loc_nargs--;
  }

  if (loc_nargs == 1)
  {
    if (ent_type(e[0]) == MATRIX_DENSE_REAL)
      rval = mdr_MinMax1_ValIdx(ent_data (e[0]), 3, ignore_inf, flat, row_dom,
                                reject_above, reject_below);
    else
    {
      fprintf (stderr, THIS_SOLVER  ": Entity type: %s\n", etd (e[0]));
      rerror (THIS_SOLVER ": operation not supported");
    }
  }
  else
  {
    fprintf (stderr, THIS_SOLVER  ": Entity types: %s and %s\n", etd (e[0]), etd (e[1]));
    rerror (THIS_SOLVER ": operation not supported");
  }

exit:

    if (e)
    {
      for (i=0; i<nargs; i++)
        ent_Clean(e[i]);
      GC_FREE(e);
    }

    return ent_Assign_Rlab_MDR(rval);
}

/* **************************************************************
 * Max-Index function...
 * ************************************************************** */
#undef THIS_SOLVER
#undef THIS_SOLVER_ENT
#define THIS_SOLVER "maxi"
#define THIS_SOLVER_ENT MaxI
Ent * THIS_SOLVER_ENT (int nargs, Datum args[])
{
  MDR *rval=0;
  Ent **e=0;
  int loc_nargs=nargs, ignore_inf=0, flat=0, row_dom=0, i;
  double *reject_below=NULL, *reject_above=NULL;
  double rej_bel,rej_ab;
  ListNode *node=0;

  if (nargs < 1)
    goto exit;

  // assign all arguments:
  e = (Ent **) GC_MALLOC(nargs * sizeof(Ent *));
  for (i=0; i<nargs; i++)
  {
    e[i] = bltin_get_ent (args[i]);
  }

  // is the last argument list of options, if so load it
  if (ent_type(e[loc_nargs-1]) == BTREE)
  {
    // reject above?
    node = btree_FindNode (ent_data (e[loc_nargs-1]), RLAB_NAME_STAT_REJECT_ABOVE);
    if (node != 0)
    {
      if (ent_type(var_ent (node))==MATRIX_DENSE_REAL)
      {
        rej_ab = class_double(var_ent (node));
        reject_above = (double *) &rej_ab;
      }
    }
    // reject below?
    node = btree_FindNode (ent_data (e[loc_nargs-1]), RLAB_NAME_STAT_REJECT_BELOW);
    if (node != 0)
    {
      if (ent_type(var_ent (node))==MATRIX_DENSE_REAL)
      {
        rej_bel = class_double(var_ent (node));
        reject_below = (double *) &rej_bel;
      }
    }
    // ignore infs
    node = btree_FindNode (ent_data (e[loc_nargs-1]), RLAB_NAME_STAT_IGNOREINF);
    if (node != 0)
    {
      if (ent_type(var_ent (node))==MATRIX_DENSE_REAL)
      {
        ignore_inf = class_int(var_ent (node));
        if (ignore_inf)
          ignore_inf = 1;
      }
    }
    // flatten a matrix?
    node = btree_FindNode (ent_data (e[loc_nargs-1]), RLAB_NAME_STAT_FLAT);
    if (node != 0)
    {
      if (ent_type(var_ent (node))==MATRIX_DENSE_REAL)
      {
        flat = class_int(var_ent (node));
        if (flat)
          flat = 1;
      }
    }
    else
    {
      // row dominant ? mutually exclusive from 'flat'
      node = btree_FindNode (ent_data (e[loc_nargs-1]), RLAB_NAME_STAT_ROWDOM);
      if (node != 0)
      {
        if (ent_type(var_ent (node))==MATRIX_DENSE_REAL)
        {
          row_dom = class_int(var_ent (node));
          if (row_dom)
            row_dom = 1;
        }
      }
    }

    // argument less one
    loc_nargs--;
  }

  if (loc_nargs == 1)
  {
    if (ent_type(e[0]) == MATRIX_DENSE_REAL)
      rval = mdr_MinMax1_ValIdx(ent_data (e[0]), 1, ignore_inf, flat, row_dom,
                                reject_above, reject_below);
    else
    {
      fprintf (stderr, THIS_SOLVER  ": Entity type: %s\n", etd (e[0]));
      rerror (THIS_SOLVER ": operation not supported");
    }
  }
  else
  {
    rerror (THIS_SOLVER ": operation not supported");
  }

exit:

    if (e)
    {
      for (i=0; i<nargs; i++)
        ent_Clean(e[i]);
      GC_FREE(e);
    }

    return ent_Assign_Rlab_MDR(rval);
}

// #undef THIS_SOLVER
// #undef THIS_SOLVER_ENT
// #undef METHOD
// #define THIS_SOLVER "maxi()"
// #define THIS_SOLVER_ENT MaxI
// #define METHOD  maxi_1_method
// Ent * THIS_SOLVER_ENT (int nargs, Datum args[])
// {
//   int type1=UNDEF, rtype=UNDEF;
//   void *rval;
//   void *(*vfptr) ();
//   Ent *e1=0;
// 
//   if (nargs != 1)
//     rerror (THIS_SOLVER ": 1 argument allowed");
// 
//   e1 = bltin_get_ent (args[0]);
//   type1 = ent_type (e1);
// 
//   rtype = METHOD[type1].type;
//   vfptr = (VFPTR) METHOD[type1].op;
// 
//   if (vfptr == 0)
//   {
//     fprintf (stderr, "Entity type: %s\n", etd (e1));
//     rerror (THIS_SOLVER " operation not supported!\n");
//   }
// 
//   rval = (*vfptr) (ent_data (e1), 1);
// 
//   ent_Clean (e1);
// 
//   return ent_Assign_Rlab_Rtype(rval,rtype);
// }

/* **************************************************************
 * Min-Index function...
 * ************************************************************** */
#undef THIS_SOLVER
#undef THIS_SOLVER_ENT
#define THIS_SOLVER "mini"
#define THIS_SOLVER_ENT MinI
Ent * THIS_SOLVER_ENT (int nargs, Datum args[])
{
  MDR *rval=0;
  Ent **e=0;
  int loc_nargs=nargs, ignore_inf=0, flat=0, row_dom=0, i;
  double *reject_below=NULL, *reject_above=NULL;
  double rej_bel,rej_ab;
  ListNode *node=0;

  if (nargs < 1)
    goto exit;

  // assign all arguments:
  e = (Ent **) GC_MALLOC(nargs * sizeof(Ent *));
  for (i=0; i<nargs; i++)
  {
    e[i] = bltin_get_ent (args[i]);
  }

  // is the last argument list of options, if so load it
  if (ent_type(e[loc_nargs-1]) == BTREE)
  {
    // reject above?
    node = btree_FindNode (ent_data (e[loc_nargs-1]), RLAB_NAME_STAT_REJECT_ABOVE);
    if (node != 0)
    {
      if (ent_type(var_ent (node))==MATRIX_DENSE_REAL)
      {
        rej_ab = class_double(var_ent (node));
        reject_above = (double *) &rej_ab;
      }
    }
    // reject below?
    node = btree_FindNode (ent_data (e[loc_nargs-1]), RLAB_NAME_STAT_REJECT_BELOW);
    if (node != 0)
    {
      if (ent_type(var_ent (node))==MATRIX_DENSE_REAL)
      {
        rej_bel = class_double(var_ent (node));
        reject_below = (double *) &rej_bel;
      }
    }
    // ignore infs
    node = btree_FindNode (ent_data (e[loc_nargs-1]), RLAB_NAME_STAT_IGNOREINF);
    if (node != 0)
    {
      if (ent_type(var_ent (node))==MATRIX_DENSE_REAL)
      {
        ignore_inf = class_int(var_ent (node));
        if (ignore_inf)
          ignore_inf = 1;
      }
    }
    // flatten a matrix?
    node = btree_FindNode (ent_data (e[loc_nargs-1]), RLAB_NAME_STAT_FLAT);
    if (node != 0)
    {
      if (ent_type(var_ent (node))==MATRIX_DENSE_REAL)
      {
        flat = class_int(var_ent (node));
        if (flat)
          flat = 1;
      }
    }
    else
    {
      // row dominant ? mutually exclusive from 'flat'
      node = btree_FindNode (ent_data (e[loc_nargs-1]), RLAB_NAME_STAT_ROWDOM);
      if (node != 0)
      {
        if (ent_type(var_ent (node))==MATRIX_DENSE_REAL)
        {
          row_dom = class_int(var_ent (node));
          if (row_dom)
            row_dom = 1;
        }
      }
    }

    // argument less one
    loc_nargs--;
  }

  if (loc_nargs == 1)
  {
    if (ent_type(e[0]) == MATRIX_DENSE_REAL)
      rval = mdr_MinMax1_ValIdx(ent_data (e[0]), -1, ignore_inf, flat, row_dom,
                                reject_above, reject_below);
    else
    {
      fprintf (stderr, THIS_SOLVER  ": Entity type: %s\n", etd (e[0]));
      rerror (THIS_SOLVER ": operation not supported");
    }
  }
  else
  {
    rerror (THIS_SOLVER ": operation not supported");
  }

exit:

    if (e)
    {
      for (i=0; i<nargs; i++)
        ent_Clean(e[i]);
      GC_FREE(e);
    }

    return ent_Assign_Rlab_MDR(rval);
}
// #undef THIS_SOLVER
// #undef THIS_SOLVER_ENT
// #undef METHOD
// #define THIS_SOLVER "mini()"
// #define THIS_SOLVER_ENT MinI
// #define METHOD  mini_1_method
// Ent * THIS_SOLVER_ENT (int nargs, Datum args[])
// {
//   int type1=UNDEF, rtype=UNDEF;
//   void *rval;
//   void *(*vfptr) ();
//   Ent *e1=0;
// 
//   if (nargs != 1)
//     rerror (THIS_SOLVER ": 1 argument allowed");
// 
//   e1 = bltin_get_ent (args[0]);
//   type1 = ent_type (e1);
// 
//   rtype = METHOD[type1].type;
//   vfptr = (VFPTR) METHOD[type1].op;
// 
//   if (vfptr == 0)
//   {
//     fprintf (stderr, "Entity type: %s\n", etd (e1));
//     rerror (THIS_SOLVER " operation not supported!\n");
//   }
// 
//   rval = (*vfptr) (ent_data (e1), -1);
// 
//   ent_Clean (e1);
// 
//   return ent_Assign_Rlab_Rtype(rval,rtype);
// }


#undef THIS_SOLVER
#define THIS_SOLVER "mini2"
Ent *
Min2I (int nargs, Datum args[])
{
  MDR *x;
  Ent *e1=0;
  int idir=-1, i, j;
  double d;

  if (nargs != 1)
    rerror (THIS_SOLVER ": 1 argument allowed");

  //
  // x
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) != MATRIX_DENSE_REAL)
    rerror (THIS_SOLVER ": 'X' has to be a real matrix!");

  x = class_matrix_real (e1);

  find_idxs_min_max(x, &i, &j, &d, &idir);

  //
  // write result to output list
  //
  Btree *bw = btree_Create ();
  install (bw, "x", ent_Create_Rlab_Double(i+1.0));
  install (bw, "y", ent_Create_Rlab_Double(j+1.0));
  install (bw, "val", ent_Create_Rlab_Double(d));
  return ent_Assign_Rlab_BTREE(bw);
}

#undef THIS_SOLVER
#define THIS_SOLVER "maxi2"
Ent *
Max2I (int nargs, Datum args[])
{
  MDR *x;
  Ent *e1=0;
  int idir=1, i, j;
  double d;

  if (nargs != 1)
    rerror (THIS_SOLVER ": 1 argument allowed");

  //
  // x
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) != MATRIX_DENSE_REAL)
    rerror ("cluster.iso: 'X' has to be a real matrix!");
  x = class_matrix_real (e1);

  find_idxs_min_max(x, &i, &j, &d, &idir);

  //
  // write result to output list
  //
  Btree *bw = btree_Create ();
  install (bw, "x",   ent_Create_Rlab_Double(i+1.0));
  install (bw, "y",   ent_Create_Rlab_Double(j+1.0));
  install (bw, "val", ent_Create_Rlab_Double(d));
  return ent_Assign_Rlab_BTREE(bw);
}
/* **************************************************************
 * Trigonometric functions...
 * ************************************************************** */

Ent * Sin (int nargs, Datum args[])
{
  int type=UNDEF, rtype=UNDEF;
  void *rval;
  void *(*vfptr) ();
  Ent *rent, *e;

  if (nargs != 1)
  {
    rerror ("sin: 1 argument allowed");
  }

  e = bltin_get_ent (args[0]);
  type = ent_type (e);

  rtype = sin_method[type].type;
  vfptr = (VFPTR) sin_method[type].op;

  if (vfptr == 0)
  {
    fprintf (stderr, "Entity type: %s\n", etd (e));
    rerror ("sin() operation not supported");
  }

  rval = (*vfptr) (ent_data (e));

  rent = ent_Create ();
  ent_data (rent) = rval;
  ent_type (rent) = rtype;

  ent_Clean (e);
  return (rent);
}

Ent * Cos (int nargs, Datum args[])
{
  int type=UNDEF, rtype=UNDEF;
  void *rval;
  void *(*vfptr) ();
  Ent *rent, *e;

  if (nargs != 1)
  {
    rerror ("cos: 1 argument allowed");
  }

  e = bltin_get_ent (args[0]);
  type = ent_type (e);

  rtype = cos_method[type].type;
  vfptr = (VFPTR) cos_method[type].op;

  if (vfptr == 0)
  {
    fprintf (stderr, "Entity type: %s\n", etd (e));
    rerror ("cos() operation not supported");
  }

  rval = (*vfptr) (ent_data (e));

  rent = ent_Create ();
  ent_data (rent) = rval;
  ent_type (rent) = rtype;

  ent_Clean (e);
  return (rent);
}

Ent *
Tan (int nargs, Datum args[])
{
  int type=UNDEF, rtype=UNDEF;
  void *rval;
  void *(*vfptr) ();
  Ent *rent, *e;

  if (nargs != 1)
  {
    rerror ("tan: 1 argument allowed");
  }

  e = bltin_get_ent (args[0]);
  type = ent_type (e);

  rtype = tan_method[type].type;
  vfptr = (VFPTR) tan_method[type].op;

  if (vfptr == 0)
  {
    fprintf (stderr, "Entity type: %s\n", etd (e));
    rerror ("tan() operation not supported");
  }

  rval = (*vfptr) (ent_data (e));

  rent = ent_Create ();
  ent_data (rent) = rval;
  ent_type (rent) = rtype;

  ent_Clean (e);
  return (rent);
}


Ent *
ASin (int nargs, Datum args[])
{
  int type=UNDEF, rtype=UNDEF;
  void *rval;
  void *(*vfptr) ();
  Ent *rent, *e;

  if (nargs != 1)
  {
    rerror ("asin: 1 argument allowed");
  }

  e = bltin_get_ent (args[0]);
  type = ent_type (e);

  vfptr = (VFPTR) asin_method[type].op;

  if (vfptr == 0)
  {
    fprintf (stderr, "Entity type: %s\n", etd (e));
    rerror ("asin() operation not supported");
  }

  rval = (*vfptr) (ent_data (e), &rtype);

  rent = ent_Create ();
  ent_data (rent) = rval;
  ent_type (rent) = rtype;

  ent_Clean (e);
  return (rent);
}

Ent *
ACos (int nargs, Datum args[])
{
  int type=UNDEF, rtype=UNDEF;
  void *rval;
  void *(*vfptr) ();
  Ent *rent, *e;

  if (nargs != 1)
  {
    rerror ("acos: 1 argument allowed");
  }

  e = bltin_get_ent (args[0]);
  type = ent_type (e);

  vfptr = (VFPTR) acos_method[type].op;

  if (vfptr == 0)
  {
    fprintf (stderr, "Entity type: %s\n", etd (e));
    rerror ("acos() operation not supported");
  }

  rval = (*vfptr) (ent_data (e), &rtype);

  rent = ent_Create ();
  ent_data (rent) = rval;
  ent_type (rent) = rtype;

  ent_Clean (e);
  return (rent);
}

Ent *
ATan (int nargs, Datum args[])
{
  int type=UNDEF, rtype=UNDEF;
  void *rval=0;
  void *(*vfptr) ();
  Ent *rent, *e;

  if (nargs != 1)
  {
    rerror ("atan: 1 argument allowed");
  }

  e = bltin_get_ent (args[0]);
  type = ent_type (e);

  vfptr = (VFPTR) atan_method[type].op;

  if (vfptr == 0)
  {
    fprintf (stderr, "Entity type: %s\n", etd (e));
    rerror ("atan() operation not supported");
  }

  rval = (*vfptr) (ent_data (e), &rtype);

  rent = ent_Create ();
  ent_data (rent) = rval;
  ent_type (rent) = rtype;

  ent_Clean (e);
  return (rent);
}

#undef THIS_SOLVER
#define THIS_SOLVER "atan2"
Ent *
ATan2 (int nargs, Datum args[])
{
  int type1=UNDEF, type2=UNDEF, rtype=UNDEF;
  void *rval=0;
  void *(*vfptr) ();
  Ent *rent, *e1=0, *e2=0;
  MDR *m1=0, *m2=0;

  if ((nargs != 2) &&(nargs != 1))
  {
    rerror (THIS_SOLVER ": " RLAB_ERROR_ONE_OR_TWO_ARG_REQUIRED "\n");
  }

  e1 = bltin_get_ent (args[0]);
  type1 = ent_type (e1);
  if (type1 == MATRIX_DENSE_REAL)
  {
    m1 = ent_data(e1);
  }

  if (nargs > 1)
  {
    e2 = bltin_get_ent (args[1]);
    type2 = ent_type (e2);
    if (type2 == MATRIX_DENSE_REAL)
    {
      m2 = ent_data(e2);
    }
  }

  vfptr = (VFPTR) atan2_method[type1][type2].op;

  if (vfptr == 0)
  {
    fprintf (stderr, THIS_SOLVER ": Entity types: %s and %s\n", etd(e1), etd(e2));
    rerror (THIS_SOLVER ": Operation not supported");
  }

  rval = (*vfptr) (m1, m2, &rtype);

  rent = ent_Create ();
  ent_data (rent) = rval;
  ent_type (rent) = rtype;

  ent_Clean (e1);
  ent_Clean (e2);
  return (rent);
}

/* **************************************************************
 * Sqrt function...
 * ************************************************************** */

Ent *
Sqrt (int nargs, Datum args[])
{
  int type=UNDEF, rtype=UNDEF;
  void *rval;
  void *(*vfptr) ();
  Ent *rent, *e;

  if (nargs != 1)
  {
    rerror ("sqrt: 1 argument allowed");
  }

  e = bltin_get_ent (args[0]);
  type = ent_type (e);

  vfptr = (VFPTR) sqrt_method[type].op;

  if (vfptr == 0)
  {
    fprintf (stderr, "Entity type: %s\n", etd (e));
    rerror ("sqrt() operation not supported");
  }

  rval = (*vfptr) (ent_data (e), &rtype);

  rent = ent_Create ();
  ent_data (rent) = rval;
  ent_type (rent) = rtype;

  ent_Clean (e);
  return (rent);
}

/* **************************************************************
 * Log function...
 * ************************************************************** */

Ent *
Log (int nargs, Datum args[])
{
  int type=UNDEF, rtype=UNDEF;
  void *rval;
  void *(*vfptr) ();
  Ent *rent=0, *e=0;

  if (nargs != 1)
  {
    rerror ("log: 1 argument allowed");
  }

  e = bltin_get_ent (args[0]);
  type = ent_type (e);

  vfptr = (VFPTR) log_method[type].op;

  if (vfptr == 0)
  {
    fprintf (stderr, "Entity type: %s\n", etd (e));
    rerror ("log() operation not supported");
  }

  rval = (*vfptr) (ent_data (e), &rtype);

  rent = ent_Create ();
  ent_data (rent) = rval;
  ent_type (rent) = rtype;

  ent_Clean (e);
  return (rent);
}

/* **************************************************************
 * Log function...
 * ************************************************************** */
Ent * Log10 (int nargs, Datum args[])
{
  int type=UNDEF, rtype=UNDEF;
  void *rval;
  void *(*vfptr) ();
  Ent *rent, *e;

  if (nargs != 1)
  {
    rerror ("log10: 1 argument allowed");
  }

  e = bltin_get_ent (args[0]);
  type = ent_type (e);

  vfptr = (VFPTR) log10_method[type].op;

  if (vfptr == 0)
  {
    fprintf (stderr, "Entity type: %s\n", etd (e));
    rerror ("log10() operation not supported");
  }

  rval = (*vfptr) (ent_data (e), &rtype);

  rent = ent_Create ();
  ent_data (rent) = rval;
  ent_type (rent) = rtype;

  ent_Clean (e);
  return (rent);
}

/* **************************************************************
 * Exp function...
 * ************************************************************** */

#undef THIS_SOLVER
#undef METHOD
#undef METHOD_WEIGHT
#define THIS_SOLVER "operation: exp"
#define METHOD exp_method
#define METHOD_WEIGHT exp_method_weight
Ent *
Exp (int nargs, Datum args[])
{
  int type1=UNDEF, rtype=UNDEF;
  void *rval=0;
  void *(*vfptr) ();
  void *(*vfptr_wgt) ();
  Ent *rent, *e1=0, *ex1=0;
  MDR *w=0, *w1=0;
  ListNode *node1=0, *node2=0;
  Btree *bw=0;

  if (nargs != 1)
  {
    rerror ("exp: 1 argument allowed");
  }

  e1 = bltin_get_ent (args[0]);
  type1 = ent_type (e1);

  if (type1 == BTREE)
  {
    // fetch x1
    node1 = btree_FindNode (ent_data (e1), RLAB_NAME_STAT_VALUE);
    if (!node1)
      rerror (THIS_SOLVER ": invalid first argument <<" RLAB_NAME_STAT_VALUE ";" RLAB_NAME_STAT_WEIGHT ">>");
    ex1   = var_ent (node1);
    type1 = ent_type(ex1);

    // fetch its weights
    node2 = btree_FindNode (ent_data (e1), RLAB_NAME_STAT_WEIGHT);
    if (node2)
    {
      if (ent_type(var_ent (node2)) != MATRIX_DENSE_REAL)
        rerror (THIS_SOLVER ": First argument list entry " RLAB_NAME_STAT_WEIGHT " must be real matrix !");
      w1 = class_matrix_real (var_ent (node2));
      if (SIZE(w1)<1)
        w1=0;
    }
  }
  else if (   type1 == MATRIX_DENSE_REAL
              || type1 == MATRIX_DENSE_COMPLEX
              || type1 == MATRIX_SPARSE_COMPLEX
              || type1 == MATRIX_SPARSE_REAL  )
  {
    ex1 = e1;
  }

  rtype = METHOD[type1].type;
  vfptr = (VFPTR) METHOD[type1].op;
  if (vfptr)
  {
    rval = (*vfptr) (ent_data (ex1));
  }
  else
  {
    fprintf (stderr, "Entity type: %s\n", etd (ex1));
    rerror ("exp() operation not supported");
  }

  if (w1)
  {
    vfptr_wgt = (VFPTR) METHOD_WEIGHT[type1].op;
    if (vfptr_wgt)
      w = (*vfptr_wgt) (ent_data (ex1), w1);
  }

  rent = ent_Create ();

  if (w)
  {
          // calculation of weights was successful
    bw = btree_Create();

      // copy result to list entry 'val'
    install  (bw, RLAB_NAME_STAT_VALUE, ent_Assign_Rlab_MDR(rval));

      // copy result to list entry 'val'
    install  (bw, RLAB_NAME_STAT_WEIGHT, ent_Assign_Rlab_MDR(w));

    ent_data (rent) = bw;
    ent_type (rent) = BTREE;
  }
  else
  {
    ent_data (rent) = rval;
    ent_type (rent) = rtype;
  }

  ent_Clean (e1);
  return (rent);
}

/* **************************************************************
 * diag() function...
 * ************************************************************** */

Ent *
Diag (int nargs, Datum args[])
{
  int k, type=UNDEF, rtype;
  Ent *e1=0, *e2=0, *rent;
  void *rval;
  void *(*vfptr) ();
  e2 = 0;

  if (nargs != 1 && nargs != 2)
  {
    rerror ("diag: 1 or 2 arguments allowed");
  }

  e1 = bltin_get_ent (args[0]);
  if (nargs == 1)
  {
    k = 0;
  }
  else
  {
    e2 = bltin_get_ent (args[1]);
    k = (int) class_double (e2);
  }

  type = ent_type (e1);

  rtype = diag_method[type].type;
  vfptr = (VFPTR) diag_method[type].op;

  if (vfptr == 0)
  {
    fprintf (stderr, "Entity type: %s\n", etd (e1));
    rerror ("diag() operation not supported");
  }

  rval = (*vfptr) (ent_data (e1), k);

  rent = ent_Create ();
  ent_data (rent) = rval;
  ent_type (rent) = rtype;

  ent_Clean (e1);
  ent_Clean (e2);

  return (rent);
}

//***************************************************************
//* Row Norm function...
//* accepts numerical second argument. recommended by E. Plischke
//**************************************************************
#undef THIS_SOLVER
#define THIS_SOLVER "rownorm"
Ent *
RowNorm (int nargs, Datum args[])
{
  char *ntype=0;
  double dtmp;
  int type;
  void *(*vfptr) ();
  Ent *e1=0, *e2=0;

  MDR *w=0;

  // Check n_args
  if ((nargs < 1) || (nargs > 2))
  {
    fprintf (stdout, THIS_SOLVER ": norm of rows of data matrix. Format:\n");
    fprintf (stdout, THIS_SOLVER ":  rownorm(x /,t/),\n");
    fprintf (stdout,
             THIS_SOLVER ": where 'x' is a vector or matrix, and t=inf(),1,2 is a\n");
    fprintf (stdout,
             THIS_SOLVER ": type of norm (number) or t=\"I\",\"1\",\"2\" (character).\n");
             rerror ("1 or 2 args allowed");
  }

  //
  // get the matrix
  //
  e1 = bltin_get_ent (args[0]);
  type = ent_type (e1);
  vfptr = (VFPTR) rownorm_method[type].op;
  if (vfptr == 0)
  {
    fprintf (stderr, "Entity type: %s\n", etd (e1));
    rerror (THIS_SOLVER ": Operation not supported");
  }

  if (nargs == 2)
  {
    e2 = bltin_get_ent (args[1]);
    if (ent_type (e2) == MATRIX_DENSE_REAL)
    {
      dtmp = class_double (e2);
      if (dtmp == create_inf ())
        ntype = cpstr ("I");
      else if (dtmp == 1)
        ntype = cpstr ("1");
      else if (dtmp == 2)
        ntype = cpstr ("f");
      else
        rerror (THIS_SOLVER ": second argument invalid");
    }
    else if (ent_type (e2) == MATRIX_DENSE_STRING)
    {
      ntype = cpstr(class_char_pointer(e2));
    }
    else
      rerror (THIS_SOLVER ": second argument invalid");
  }

  if (!ntype)
    w = (MDR *) (*vfptr) (ent_data (e1), "f");
  else
    w = (MDR *) (*vfptr) (ent_data (e1), ntype);

  ent_Clean (e1);
  ent_Clean (e2);

  return ent_Assign_Rlab_MDR(w);
}

//***************************************************************
//* Matrix Norm function...
//* accepts numerical second argument. recommended by E. Plischke
//**************************************************************
Ent *
MNorm (int nargs, Datum args[])
{
  char *ntype = 0;
  double dtmp, *rval=0;
  int type;
  void *(*vfptr) ();
  Ent *e1=0, *e2=0, *rent;
  rent = 0;

  // Check n_args
  if ((nargs < 1) || (nargs > 2))
  {
    fprintf (stdout, "mnorm: Matrix/vector norm. Format:\n");
    fprintf (stdout, "mnorm:  mnorm(x /,t/),\n");
    fprintf (stdout,
             "mnorm: where 'x' is a vector or matrix, and t=inf(),1,2 is a\n");
    fprintf (stdout,
             "mnorm: type of norm (number) or t=\"I\",\"1\",\"2\" (character).\n");
    rerror ("1 or 2 args allowed");
  }

  rent = ent_Create ();

  if (nargs == 1)
  {
    // Compute the matrix 1-norm.
    e1 = bltin_get_ent (args[0]);
    type = ent_type (e1);
    vfptr = (VFPTR) mnorm_method[type].op;
    if (vfptr == 0)
    {
      fprintf (stderr, "Entity type: %s\n", etd (e1));
      rerror ("mnorm() operation not supported");
    }

    rval = (double *) (*vfptr) (ent_data (e1), "1");
  }
  else if (nargs == 2)
  {
    e1 = bltin_get_ent (args[0]);
    e2 = bltin_get_ent (args[1]);
    if (ent_type (e2) == MATRIX_DENSE_REAL)
    {
      dtmp = class_double (e2);
      if (dtmp == create_inf ())
        ntype = cpstr ("I");
      else if (dtmp == 1)
        ntype = cpstr ("1");
      else if (dtmp == 2)
        ntype = cpstr ("2");
      else
        rerror ("mnorm(): second argument invalid");
    }
    else if (ent_type (e2) == MATRIX_DENSE_STRING)
    {
      ntype = cpstr(class_char_pointer(e2));
    }
    else
      rerror ("mnorm(): second argument invalid");

    type = ent_type (e1);
    vfptr = (VFPTR) mnorm_method[type].op;
    if (vfptr == 0)
    {
      fprintf (stderr, "Entity type: %s\n", etd (e1));
      rerror ("mnorm() operation not supported");
    }

    rval = (double *) (*vfptr) (ent_data (e1), ntype);
  }

  if (ntype)
    GC_FREE (ntype);

  ent_Clean (e1);
  ent_Clean (e2);

  ent_data (rent) = mdr_CreateScalar (*rval);
  if (rval)
    GC_FREE (rval);

  ent_type (rent) = MATRIX_DENSE_REAL;
  return (rent);
}

/* **************************************************************
 * P-Norm function...
 * ************************************************************** */

Ent *
PNorm (int nargs, Datum args[])
{
  double *rval=0, pval;
  int type;
  void *(*vfptr) ();
  Ent *e1, *e2, *rent;

  /* Check n_args */
  if ((nargs != 2))
    rerror ("pnorm: 2 arguments required");

  e1 = bltin_get_ent (args[0]);
  e2 = bltin_get_ent (args[1]);

  /* Get the P-value. */
  pval = class_double (e2);

  type = ent_type (e1);
  vfptr = (VFPTR) pnorm_method[type].op;

  if (vfptr == 0)
  {
    fprintf (stderr, "Entity type: %s\n", etd (e1));
    rerror ("pnorm() operation not supported");
  }

  rval = (double *) (*vfptr) (ent_data (e1), pval);

  ent_Clean (e1);
  ent_Clean (e2);

  rent = ent_Create ();
  ent_data (rent) = mdr_CreateScalar (*rval);
  if (rval)
    GC_FREE (rval);

  ent_type (rent) = MATRIX_DENSE_REAL;
  return (rent);
}

// ***************************************************************
// * Eigenvalues: General problem
// ***************************************************************

Ent *
Eig (int nargs, Datum args[])
{
  int t1=0, t2=0, t3, issym=-1, issympos=0, nm=0;
  Ent *e1=0, *e2=0, *e3=0, *rent;
  rent = 0;
  MDR *m1, *m2;
  MDC *c1, *c2;

  /* Check n_args */
  if ((nargs < 1) || (nargs > 3))
  {
    fprintf (stdout,
             "eig: Eigenvalues of a real or complex dense matrix.\n");
    fprintf (stdout,
             "eig: Format:\n");
    fprintf (stdout,
             "eig:   (1) y=eig(a/,flag/),\n");
    fprintf (stdout,
             "eig:   (2) y=eig(a,b/,flag/),\n");
    fprintf (stdout,
             "eig: where 'a' and 'b' are the square matrices of same\n");
    fprintf (stdout,
             "eig: size. The result is a list  y=<<val;vec>> , where  'val' is a\n");
    fprintf (stdout,
             "eig: row-vector, while  'vec'  is a matrix of column eigenvectors.\n");
    fprintf (stdout,
             "eig: For a single matrix argument the parameter flags=\"S\", \"G\"\n");
    fprintf (stdout,
             "eig: if  'a'  is symmetric or non-symmetric, and for two matrix\n");
    fprintf (stdout,
             "eig: flags=\"SP\"  if  'a'  is symmetric and  'b'  is positive-definite.\n");
    rerror ("one, two or three arguments allowed");
  }

  if (nargs >= 1)
  {
    //
    // Get the matrix 'a'
    //
    e1 = bltin_get_ent (args[0]);
    t1 = ent_type (e1);
    if( t1 == MATRIX_DENSE_REAL ||
        t1 == MATRIX_DENSE_COMPLEX) nm = 1;
    else
      rerror("eig: first argument has to be a dense real or complex matrix!");
  }
  if (nargs >= 2)
  {
    //
    // Get the second argument: a matrix of size 'a' or a string
    //
    e2 = bltin_get_ent (args[1]);
    t2 = ent_type (e2);
    if (t2 == MATRIX_DENSE_STRING)
    {
      if (strcmp (class_char_pointer(e2), "s") == 0||
          strcmp (class_char_pointer(e2), "S") == 0)
        issym = 1; // 'a' symmetric
      if (strcmp (class_char_pointer(e2), "g") == 0||
          strcmp (class_char_pointer(e2), "G") == 0)
        issym = 0; // 'a' non-symmetric, don't do test
    }
    else if (t2 == MATRIX_DENSE_REAL || t2 == MATRIX_DENSE_COMPLEX)
      nm = 2;
    else
      rerror("eig: second argument has to be a dense real"
             "or complex matrix or a string!");
  }
  if (nargs == 3)
  {
    //
    // Get the third argument: a string
    //
    e3 = bltin_get_ent (args[2]);
    t3 = ent_type (e3);
    if (t3 == MATRIX_DENSE_STRING)
    {
      if (strcmp (class_char_pointer(e3), "sp") == 0||
          strcmp (class_char_pointer(e3), "SP") == 0)
        issympos = 1; // 'a' symmetric, 'b' positive
    }
    else
      rerror("eig: third argument, if given, has to be a string!");
  }
  //
  // Now find the eigenvalues
  //
  if (nm == 1)
  {
    //
    // eig(a/,flag/)
    //
    if(t1 == MATRIX_DENSE_REAL)
    {
      // copy 'a' to m1
      m1 = mdr_Float_BF ( class_matrix_real (e1) );
      rent = ent_Create ();
      ent_data (rent) = ((void *) mdr_EigS(m1, &issym));
      ent_type (rent) = BTREE;
      mdr_Destroy(m1);
    }
    else if(t1 == MATRIX_DENSE_COMPLEX)
    {
      // complex 'a'
      c1 = mdc_Copy( ent_data(e1) );
      rent = ent_Create ();
      ent_data (rent) = (void *) mdc_EigS(c1, &issym);
      ent_type (rent) = BTREE;
      if (issym == 0) mdc_Destroy(c1);
    }
  }
  else if (nm == 2)
  {
    //
    // eig(a,b)
    //
    if(t1==MATRIX_DENSE_REAL && t2==MATRIX_DENSE_REAL)
    {
      m1 = mdr_Float_BF ( class_matrix_real (e1) );
      m2 = mdr_Float_BF ( class_matrix_real (e2) );
      rent = ent_Create ();
      ent_data (rent) = (void *) mdr_EigG(m1, m2, &issympos);
      ent_type (rent) = BTREE;
      if (issympos == 0) mdr_Destroy(m1);
      mdr_Destroy( m2 );
    }
    else
    {
      if (t1==MATRIX_DENSE_REAL && t2==MATRIX_DENSE_COMPLEX)
      {
        c1 = mdr_coerce_mdc(ent_data (e1));
        c2 = mdc_Copy( ent_data(e2) );
      }
      else if (t2==MATRIX_DENSE_REAL && t1==MATRIX_DENSE_COMPLEX)
      {
        c1 = mdc_Copy( ent_data (e1) );
        c2 = mdr_coerce_mdc( ent_data (e2) );
      }
      else
      {
        c1 = mdc_Copy( ent_data (e1) );
        c2 = mdc_Copy( ent_data (e2) );
      }
      rent = ent_Create ();
      ent_data (rent) = (void *) mdc_EigG(c1, c2, &issympos);
      ent_type (rent) = BTREE;
      // clean-up
      if (issympos == 0) mdc_Destroy( c1 );
      mdc_Destroy( c2 );
    }
  }

  ent_Clean(e1);
  ent_Clean(e2);
  ent_Clean(e3);
  return (rent);
}

/*
 *  eigs()
 */
#ifdef HAVE_ARPACK

Ent * Eigs (int nargs, Datum args[])
{
  void *rval=0;
  int t1, t2=UNDEF, t3, t4;
  char which[3]="SM";
  int nev=1, havem=0;
  Ent *e1=0, *e2=0, *e3=0, *e4=0, *rent;
  rent = 0;
  Complex sigma = 0.0 + 0.0i;

  /* Check n_args */
  if ((nargs < 1) || (nargs > 4))
  {
    fprintf(stdout,
            "eigs: Selected eigenvalues, standard or general, of a dense or sparse\n");
    fprintf(stdout,
            "eigs: matrix. Format:\n");
    fprintf(stdout,
            "eigs:   y=eigs(a,/b/,sigma,k/),\n");
    fprintf(stdout,
            "eigs: where 'a' and 'b' are the square matrices of the same\n");
    fprintf(stdout,
            "eigs: size. The result  y=<<val;vec>>  is a list, where 'val' is the\n");
    fprintf(stdout,
            "eigs: row-vector of eigenvalues, and  'val'  is the matrix of column\n");
    fprintf(stdout,
            "eigs: eigenvector.\n");
    rerror ("at least one argument required");
  }

  //
  // eigs(A, )
  //
  e1 = bltin_get_ent (args[0]);
  t1 = ent_type(e1);
  if(
     t1 != MATRIX_DENSE_REAL &&
     t1 != MATRIX_DENSE_COMPLEX &&
     t1 != MATRIX_SPARSE_REAL &&
     t1 != MATRIX_SPARSE_COMPLEX
    ) rerror("eigs: improper first argument!");
  havem = 1;

  //
  // check second arguments
  //
  if(nargs >= 2)
  {
    e2 = bltin_get_ent( args[1] );
    t2 = ent_type(e2);
    if( MNR(ent_data(e1))==MNR(ent_data(e2)) && MNC(ent_data(e1))==MNC(ent_data(e2)) )
      havem = 2;
    else
    {
      if(t2==MATRIX_DENSE_REAL)
      {
        // eigs(A,nev, )  or  eigs(A,1, )
        nev = (int) class_double(e2);
      }
    }
  }

  //
  // check third argument
  //
  if(nargs >=3 )
  {
    e3 = bltin_get_ent( args[2] );
    t3 = ent_type(e3);
    if(havem==2)
    {
      if(t3 == MATRIX_DENSE_REAL)
      {
        // eigs(A,B,nev, )  or  eigs(A,B,1, )
        nev = (int) class_double( e3 );
      }
    }
    else
    {
      // eigs(A,nev,sigma)
      if ((t3 == MATRIX_DENSE_REAL) || (t3 == MATRIX_DENSE_COMPLEX))
      {
        // sigma is a single real number
        sigma = class_complex(e3);
        sprintf(which,"LM");
      }
      else if(t3 == MATRIX_DENSE_STRING)
      {
        // sigma is a two-character string
        sprintf(which,"%s", class_char_pointer(e3) );
      }
      else if(t3 == BTREE)
      {
        // sigma is a list: <<sigma;which>>
        ListNode *snode, *wnode;
        snode = btree_FindNode (ent_data (e3), "sigma");
        wnode = btree_FindNode (ent_data (e3), "which");
        if (snode == 0 && wnode==0)
          rerror("Missing 'sigma' and 'which' entries in the list!");
        if (wnode != 0)
        {
          sprintf(which,"%s", class_char_pointer(var_ent(wnode)) );
        }
        if (snode !=  0)
        {
          if ((    ent_type(var_ent(snode))==MATRIX_DENSE_REAL)
               || (ent_type(var_ent(snode))==MATRIX_DENSE_COMPLEX))
          {
            sigma = class_complex(var_ent(snode));
          }
          else rerror("Improper type of entry 'sigma'!");
        }
      }
    }
  }

  //
  // check fourth argument
  //
  if(nargs ==4 )
  {
    e4 = bltin_get_ent( args[3] );
    t4 = ent_type(e4);
    // eigs(A,B,nev,sigma)
    if ((t4 == MATRIX_DENSE_REAL) || (t4 == MATRIX_DENSE_COMPLEX))
    {
      // sigma is real number
      sigma = class_complex( e4 );
    }
    else if(t4 == MATRIX_DENSE_STRING)
    {
      // sigma is a two-character string
      sprintf(which,"%s", class_char_pointer(e4) );
    }
    else if(t4 == BTREE)
    {
      // sigma is a list: <<sigma;which>>
      ListNode *snode, *wnode;
      snode = btree_FindNode (ent_data (e4), "sigma");
      wnode = btree_FindNode (ent_data (e4), "which");
      if (snode == 0 && wnode==0)
        rerror("Missing 'sigma' and 'which' entries in the list!");
      if (wnode != 0)
        sprintf(which,"%s", class_char_pointer(var_ent(wnode)) );
      if (snode !=  0)
      {
        if ( (ent_type(var_ent(snode))==MATRIX_DENSE_REAL) || (ent_type(var_ent(snode))==MATRIX_DENSE_COMPLEX))
          sigma = class_complex(var_ent(snode));
        else rerror("Improper type of entry 'sigma'!");
      }
    }
  }

  //
  // check which: only the following choices are allowed
  //
  which[0] = toupper( which[0] );
  which[1] = toupper( which[1] );
  if ( ! (
      (strcmp (which, "LM") == 0) ||
      (strcmp (which, "SM") == 0) ||
      (strcmp (which, "LA") == 0) ||
      (strcmp (which, "SA") == 0) ||
      (strcmp (which, "BE") == 0) ||
      (strcmp (which, "LR") == 0) ||
      (strcmp (which, "SR") == 0) ||
      (strcmp (which, "LI") == 0) ||
      (strcmp (which, "SI") == 0)
         ) ) sprintf(which,"SM");

  if(havem==1)
  {
    //
    // eigs(A,nev,sigma)
    //
    if (t1 == MATRIX_DENSE_REAL && IM(sigma)==0)
    {
      //
      // real dense problem
      //
      rval = (void *) mdr_ar_EigsS(ent_data (e1), nev, which, RE(sigma));
    }
    else if ((t1 == MATRIX_DENSE_REAL && IM(sigma)!=0)|| (t1 == MATRIX_DENSE_COMPLEX ))
    {
      //
      // complex dense problem
      //
      MDC * c1;
      if (t1 == MATRIX_DENSE_REAL)
        c1 = mdr_coerce_mdc(ent_data (e1));
      else
        c1 = ent_data (e1);

      rval = (void *) mdc_ar_EigsS(c1, nev, which, sigma);

      if (t1 == MATRIX_DENSE_REAL)
        mdc_Destroy( c1 );
    }
    else if(t1 == MATRIX_SPARSE_REAL && IM(sigma)==0)
    {
      //
      // real sparse problem
      //
      rval = (void *) msr_ar_EigsS(ent_data (e1), nev, which, RE(sigma));
    }
    else if ((t1 == MATRIX_SPARSE_REAL && IM(sigma)!=0)||
              (t1 == MATRIX_SPARSE_COMPLEX ))
    {
      //
      // complex sparse problem
      //
      MSC * c1;
      if (t1 == MATRIX_SPARSE_REAL)
        c1 = msr_coerce_msc(ent_data (e1));
      else
        c1 = ent_data (e1);
      //
      rval = (void *) msc_ar_EigsS(c1, nev, which, sigma);
      //
      if (t1 == MATRIX_SPARSE_REAL) msc_Destroy( c1 );
    }
    else rerror("eigs: Unsupported data types!");
  }
  else if (havem == 2)
  {
    // eigs(A,B,nev,sigma)
    if (t1 == MATRIX_DENSE_REAL &&
        t2 == MATRIX_DENSE_REAL &&
        IM(sigma) == 0)
    {
      //
      // real dense problem
      //
      rval = (void *) mdr_ar_EigsG(ent_data (e1), ent_data (e2), nev, which, RE(sigma));
    }
    else if ((t1 == MATRIX_DENSE_COMPLEX && t2 == MATRIX_DENSE_REAL) ||
             (t2 == MATRIX_DENSE_COMPLEX && t1 == MATRIX_DENSE_REAL) ||
             (t1 == MATRIX_DENSE_COMPLEX && t2 == MATRIX_DENSE_COMPLEX) ||
             (t1 == MATRIX_DENSE_REAL && t2 == MATRIX_DENSE_REAL && IM(sigma)!= 0))
    {
      //
      // complex dense problem
      //
      MDC *c1=0, *c2=0;
      //
      if (t1 == MATRIX_DENSE_REAL)
        c1 = mdr_coerce_mdc(ent_data (e1));
      else
        c1 = ent_data( e1 );
      //
      if (t2 == MATRIX_DENSE_REAL)
        c2 = mdr_coerce_mdc(ent_data (e2));
      else
        c2 = ent_data( e2 );
      rval = (void *) mdc_ar_EigsG(c1, c2, nev, which, sigma);

      //
      if (t1 == MATRIX_DENSE_REAL)
        mdc_Destroy( c1 );
      if (t2 == MATRIX_DENSE_REAL)
        mdc_Destroy( c2 );
    }
    else if (t1 == MATRIX_SPARSE_REAL &&
             t2 == MATRIX_SPARSE_REAL &&
             IM(sigma) == 0)
    {
      //
      // real sparse problem
      //
      rval = (void *) msr_ar_EigsG(ent_data (e1), ent_data (e2), nev, which, RE(sigma));
    }
    else if ((t1 == MATRIX_SPARSE_COMPLEX && t2 == MATRIX_SPARSE_REAL) ||
             (t2 == MATRIX_SPARSE_COMPLEX && t1 == MATRIX_SPARSE_REAL) ||
             (t1 == MATRIX_SPARSE_COMPLEX && t2 == MATRIX_DENSE_COMPLEX) ||
             (t1 == MATRIX_SPARSE_REAL    && t2 == MATRIX_SPARSE_REAL && IM(sigma)!= 0))
    {
      //
      // complex sparse problem
      //
      MSC *c1=0, *c2=0;
      //
      if (t1 == MATRIX_SPARSE_REAL)
        c1 = msr_coerce_msc(ent_data (e1));
      else
        c1 = ent_data( e1 );
      //
      if (t2 == MATRIX_SPARSE_REAL)
        c2 = msr_coerce_msc(ent_data (e2));
      else
        c2 = ent_data( e2 );
      rval = (void *) msc_ar_EigsG(c1, c2, nev, which, sigma);
      //
      if (t1 == MATRIX_SPARSE_REAL)
        msc_Destroy( c1 );
      if (t2 == MATRIX_SPARSE_REAL)
        msc_Destroy( c2 );
    }
    else rerror("eigs: A,B have to be both dense or sparse!");
  }

  ent_Clean (e1);
  ent_Clean (e2);
  ent_Clean (e3);
  ent_Clean (e4);

  rent = ent_Create ();
  ent_data (rent) = rval;
  ent_type (rent) = BTREE;
  return (rent);
}

#endif

/* **************************************************************
 * Factor
 * ************************************************************** */

Ent *
Factor (int nargs, Datum args[])
{
  void *rval;
  int type1, rtype;
  void *(*vfptr) ();
  Ent *e1=0, *e2=0, *rent;
  int flag;
  char *str;

  /* Check n_args */
  if ((nargs < 1) || (nargs > 2))
    rerror ("factor: 1 or 2 args allowed");

  /* Get the coefficient matrix. */
  e1 = bltin_get_ent (args[0]);
  type1 = ent_type (e1);

  flag = 0;
  if (nargs == 2)
  {
    e2 = bltin_get_ent (args[1]);
    if (ent_type (e2) != MATRIX_DENSE_STRING)
      rerror ("factor: Invalid 2nd argument");

    str = class_char_pointer (e2);
    switch (*str)
    {
      case 's':
      case 'S':
        flag = 2;
        break;

      case 'g':
      case 'G':
        flag = 1;
        break;

      default:
        break;
    }
  }

  rtype = factor_method[type1].type;
  vfptr = (VFPTR) factor_method[type1].op;

  if (vfptr == 0)
  {
    fprintf (stderr, "Entity type: %s\n", etd (e1));
    rerror ("factor() operation not supported");
  }

  rval = (*vfptr) (ent_data (e1), flag);

  rent = ent_Create ();
  ent_data (rent) = rval;
  ent_type (rent) = rtype;

  ent_Clean (e1);
  ent_Clean (e2);

  return (rent);
}

/* **************************************************************
 * Backsub: This function takes two arguments:
 *    A list that is output from the function factor()
 *    A matrix of vectors that is the new right-hand-side of
 *                    A * x = B
 *
 * To figure out who to call to do the work, Backsub must look at
 * the supplied list, and check out the names and the types of
 * the elements. If the list contains:
 *
 *   lu     General Solution (LU Decomposition was performed).
 *   ldl    Symmetric Solution (Cholesky decomposition was performed).
 *
 * Once we have make the General versus Symmetric determination, then
 * we can look at the type (REAL or COMPLEX) of the lu, or ldl matrix
 * to make the final determination. Note that we must pay attention
 * to the type of the right-hand-side (B).
 * ************************************************************** */

Ent *
Backsub (int nargs, Datum args[])
{
  void *rval;
  int rtype=UNDEF, type2, ltype;
  void *(*vfptr) ();
  Btree *bl;
  Ent *e1, *e2, *rent;
  ListNode *tmp;
  rtype = 0;
  vfptr = 0;

  /* Check n_args */
  if (nargs != 2)
    rerror ("backsub: requires 2 arguments");

  /* Get the arguments, and their types. */
  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) != BTREE)
    rerror ("backsub: 1st argument must be a LIST");
  bl = ent_data (e1);

  // get the second argument
  e2 = bltin_get_ent (args[1]);
  type2 = ent_type (e2);

  /* Make the Ge or Sym determination. */
  if ((tmp = btree_FindNode (bl, "lu")))
  {
    /* Now decide which function to call. */
    ltype = ent_type (var_ent (tmp));
    rtype = backsub_ge_method[ltype][type2].type;
    vfptr = (VFPTR) backsub_ge_method[ltype][type2].op;
  }
  else if ((tmp = btree_FindNode (bl, "ldl")))
  {
    /* Now decide which function to call. */
    ltype = ent_type (var_ent (tmp));
    rtype = backsub_sy_method[ltype][type2].type;
    vfptr = (VFPTR) backsub_sy_method[ltype][type2].op;
  }
  else if ((tmp = btree_FindNode (bl, "value")))
  {
    /* UMFPACK structure for sparse factorization. */
    ltype = ent_type (var_ent (tmp));
    rtype = backsub_sparse_method[ltype][type2].type;
    vfptr = (VFPTR) backsub_sparse_method[ltype][type2].op;
  }
  else if ((tmp = btree_FindNode (bl, "perm_c")))
  {
    /* SuperLU structure for sparse factorization. */
    ltype = ent_type (var_ent (tmp));
    rtype = backsub_sparse_method[ltype][type2].type;
    vfptr = (VFPTR) backsub_sparse_method[ltype][type2].op;
  }
  else
  {
    rerror ("backsub: list must contain either \"lu\" or \"ldl\" elements");
  }

  /* Now call the correct backsub function. */
  if (vfptr == 0)
  {
    fprintf (stderr, "Entity type: %s\n", etd (var_ent (tmp)));
    rerror ("backsub() operation not supported");
  }
  rval = (*vfptr) (ent_data (e1), ent_data (e2));

  rent = ent_Create ();
  ent_data (rent) = rval;
  ent_type (rent) = rtype;
  ent_Clean (e1);
  ent_Clean (e2);

  return (rent);
}

/* **************************************************************
 * SVD
 * ************************************************************** */

Ent *
Svd (int nargs, Datum args[])
{
  void *rval=0;
  int rtype=UNDEF, type1;
  void *(*vfptr) ();
  Ent *e1=0, *e2=0, *rent;
  int flag=1;
  char *str=0;

  /* Check n_args */
  if ((nargs < 1) || (nargs > 2))
    rerror ("svd: 1 or 2 args allowed");

  /* Get the coefficient matrix. */
  e1 = bltin_get_ent (args[0]);
  type1 = ent_type (e1);

  if (nargs == 2)
  {
    e2 = bltin_get_ent (args[1]);
    if (ent_type(e2) == MATRIX_DENSE_STRING)
    {
      str = class_char_pointer (e2);

      if (str)
      {
        switch (*str)
        {
          case 's':
          case 'S':
            flag = 2;
            break;

          case 'a':
          case 'A':
            flag = 1;
            break;

          case 'n':
          case 'N':
            flag = 3;
            break;

          default:
            break;
        }
      }
    }
  }

  rtype = svd_method[type1].type;
  vfptr = (VFPTR) svd_method[type1].op;

  if (vfptr == 0)
  {
    fprintf (stderr, "Entity type: %s\n", etd (e1));
    rerror ("svd() operation not supported");
  }

  rval = (*vfptr) (ent_data (e1), flag);

  rent = ent_Create ();
  ent_data (rent) = rval;
  ent_type (rent) = rtype;

  ent_Clean (e1);
  ent_Clean (e2);

  return (rent);
}

/* **************************************************************
 * Cholesky Decomposition
 * ************************************************************** */

Ent *
Chol (int nargs, Datum args[])
{
  void *rval;
  int rtype=UNDEF, type1;
  void *(*vfptr) ();
  Ent *e1, *rent;

  /* Check n_args */
  if (nargs != 1)
    rerror ("chol: 1 argument allowed");

  /* Get the coefficient matrix. */
  e1 = bltin_get_ent (args[0]);
  type1 = ent_type (e1);

  rtype = chol_method[type1].type;
  vfptr = (VFPTR) chol_method[type1].op;

  if (vfptr == 0)
  {
    fprintf (stderr, "Entity type: %s\n", etd (e1));
    rerror ("chol() operation not supported");
  }

  rval = (*vfptr) (ent_data (e1));

  rent = ent_Create ();
  ent_data (rent) = rval;
  ent_type (rent) = rtype;

  ent_Clean (e1);
  return (rent);
}

/* **************************************************************
 * Solve
 * ************************************************************** */
#undef THIS_SOLVER
#define THIS_SOLVER "solve"
Ent * Solve (int nargs, Datum args[])
{
  void *rval=0;
  int rtype=UNDEF, type1=UNDEF, type2=UNDEF;
  char *ctype=0;
  void *(*vfptr) ();
  Ent *e1=0, *e2=0, *e3=0, *rent=0;

  /* Check n_args */
  if ((nargs < 2) || (nargs > 3))
    rerror ("solve: 2 or 3 args allowed");

  e1 = bltin_get_ent (args[0]);
  if (!e1)
    rerror (THIS_SOLVER ": " RLAB_ERROR_ARG1_NUM "\n");
  type1 = ent_type (e1);

  e2 = bltin_get_ent (args[1]);
  if (!e2)
    rerror (THIS_SOLVER ": " RLAB_ERROR_ARG2_NUM "\n");
  type2 = ent_type (e2);

  rtype = solve_method[type1][type2].type;
  vfptr = (VFPTR) solve_method[type1][type2].op;

  if (vfptr == 0)
  {
    fprintf (stderr, "Entity type: %s and %s\n", etd (e1), etd (e2));
    rerror ("solve() operation not supported");
  }

  if (nargs == 3)
  {
    e3 = bltin_get_ent (args[2]);
    ctype = class_char_pointer (e3);
    if (isvalidstring(ctype)!=1)
      ctype = 0;
  }
  else
    ctype = 0;

  rval = (*vfptr) (ent_data (e1), ent_data (e2), ctype);

  ent_Clean (e1);
  ent_Clean (e2);
  ent_Clean (e3);

  rent = ent_Create ();
  ent_data (rent) = rval;
  ent_type (rent) = rtype;
  return (rent);
}

/* **************************************************************
 * Rcond function...
 * ************************************************************** */

Ent *
Rcond (int nargs, Datum args[])
{
  int type;
  double *rval=0;
  void *(*vfptr) ();
  Ent *e1=0, *rent;

  /* Check n_args */
  if ((nargs != 1))
    rerror ("rcond: 1 argument required");

  e1 = bltin_get_ent (args[0]);

  type = ent_type (e1);
  vfptr = (VFPTR) rcond_method[type].op;

  if (vfptr == 0)
  {
    fprintf (stderr, "Entity type: %s\n", etd (e1));
    rerror ("rcond() operation not supported");
  }

  rval = (double *) (*vfptr) (ent_data (e1));

  ent_Clean (e1);

  rent = ent_Create ();

  ent_data (rent) = mdr_CreateScalar (*rval);
  if (rval)
    GC_FREE (rval);


  ent_type (rent) = MATRIX_DENSE_REAL;
  return (rent);
}

/* **************************************************************
 * Compute the Hessenberg form of a matrix. Also compute an
 * orthogonal matrix [p], such that, [m] = [p][h][p]'.
 * ************************************************************** */

Ent *
Hess (int nargs, Datum args[])
{
  void *rval;
  int rtype=UNDEF, type;
  void *(*vfptr) ();
  Ent *e, *rent;
  rent = 0;

  /* Check n_args */
  if (nargs != 1)
    rerror ("hess: 1 argument allowed");

  e = bltin_get_ent (args[0]);
  type = ent_type (e);

  rtype = hess_method[type].type;
  vfptr = (VFPTR) hess_method[type].op;

  if (vfptr == 0)
  {
    fprintf (stderr, "Entity type: %s\n", etd (e));
    rerror ("hess() operation not supported");
  }

  rval = (*vfptr) (ent_data (e));

  rent = ent_Create ();
  ent_data (rent) = rval;
  ent_type (rent) = rtype;

  ent_Clean (e);
  return (rent);
}

/* **************************************************************
 * ************************************************************** */

Ent *
Balance (int nargs, Datum args[])
{
  void *rval;
  int rtype=UNDEF, type;
  void *(*vfptr) ();
  Ent *e, *rent;
  rent = 0;

  /* Check n_args */
  if (nargs != 1)
    rerror ("balance: 1 argument allowed");

  e = bltin_get_ent (args[0]);
  type = ent_type (e);

  rtype = balance_method[type].type;
  vfptr = (VFPTR) balance_method[type].op;

  if (vfptr == 0)
  {
    fprintf (stderr, "Entity type: %s\n", etd (e));
    rerror ("balance() operation not supported");
  }

  rval = (*vfptr) (ent_data (e));

  rent = ent_Create ();
  ent_data (rent) = rval;
  ent_type (rent) = rtype;

  ent_Clean (e);
  return (rent);
}

/* **************************************************************
 * QR
 * ************************************************************** */

Ent *
QR (int nargs, Datum args[])
{
  void *rval=0;
  int rtype=0, type1;
  void *(*vfptr) ();
  Ent *e1=0, *e2=0, *rent;
  char *str=0;

  /* Check n_args */
  if ((nargs < 1) || (nargs > 2))
    rerror ("qr: 1 or 2 args allowed");

  /* Get the coefficient matrix. */
  e1 = bltin_get_ent (args[0]);
  type1 = ent_type (e1);

  rent = ent_Create ();

  if (nargs == 1)
  {
    rtype = qr_method[type1].type;
    vfptr = (VFPTR) qr_method[type1].op;

    if (vfptr == 0)
    {
      fprintf (stderr, "Entity type: %s\n", etd (e1));
      rerror ("qr() operation not supported");
    }

    rval = (*vfptr) (ent_data (e1));
  }
  else if (nargs == 2)
  {
    e2 = bltin_get_ent (args[1]);
    if (ent_type (e2) == MATRIX_DENSE_STRING)
    {
      str = class_char_pointer (e2);
      if (str)
        switch(*str)
        {
          case 'p':
          case 'P':
            break;

          default:
            break;
        }
    }

    rtype = qrp_method[type1].type;
    vfptr = (VFPTR) qrp_method[type1].op;
    if (vfptr == 0)
    {
      fprintf (stderr, "Entity type: %s\n", etd (e1));
      rerror ("qr() operation not supported");
    }

    rval = (*vfptr) (ent_data (e1));

  }

  ent_Clean (e1);
  ent_Clean (e2);

  ent_data (rent) = rval;
  ent_type (rent) = rtype;
  return (rent);
}

/* **************************************************************
 * Schur decomposition
 * ************************************************************** */

Ent *
Schur (int nargs, Datum args[])
{
  void *rval;
  int rtype=UNDEF, type;
  void *(*vfptr) ();
  Ent *e, *rent;
  rent = 0;

  /* Check n_args */
  if (nargs != 1)
    rerror ("schur: 1 argument allowed");

  e = bltin_get_ent (args[0]);
  type = ent_type (e);

  rtype = schur_method[type].type;
  vfptr = (VFPTR) schur_method[type].op;

  if (vfptr == 0)
  {
    fprintf (stderr, "Entity type: %s\n", etd (e));
    rerror ("schur() operation not supported");
  }

  rval = (*vfptr) (ent_data (e));

  rent = ent_Create ();
  ent_data (rent) = rval;
  ent_type (rent) = rtype;

  ent_Clean (e);
  return (rent);
}

/* **************************************************************
 * Sylvester Matrix Equations
 * ************************************************************** */

Ent *
Sylv (int nargs, Datum args[])
{
  void *rval;
  int rtype=UNDEF, type1, type2;
  void *(*vfptr) ();
  Ent *e1, *e2, *e3, *rent;
  rent = 0;

  /* Check n_args */
  if (!(nargs == 2 || nargs == 3))
    rerror ("sylv: 2 or 3 arguments are required");

  /* get arg from list */
  if (nargs == 2)
  {
    /* Lyapunov */
    e1 = bltin_get_ent (args[0]);
    e2 = bltin_get_ent (args[1]);

    type1 = ent_type (e1);
    type2 = ent_type (e2);

    rtype = sylv_method[type1][type2].type;
    vfptr = (VFPTR) sylv_method[type1][type2].op;

    if (vfptr == 0)
    {
      fprintf (stderr, "Entity type: %s and %s\n", etd (e1), etd (e2));
      rerror ("sylv operation not supported");
    }

    rval = (*vfptr) (ent_data (e1), 0, ent_data (e2));

    rent = ent_Create ();
    ent_data (rent) = rval;
    ent_type (rent) = rtype;

    ent_Clean (e1);
    ent_Clean (e2);
  }
  else if (nargs == 3)
  {
    /* Sylvester */
    e1 = bltin_get_ent (args[0]);
    e2 = bltin_get_ent (args[1]);
    e3 = bltin_get_ent (args[2]);

    type1 = ent_type (e1);
    type2 = ent_type (e2);

    rtype = sylv_method[type1][type2].type;
    vfptr = (VFPTR) sylv_method[type1][type2].op;

    if (vfptr == 0)
    {
      fprintf (stderr, "Entity type: %s, %s and %s\n",
	       etd (e1), etd (e2), etd (e3));
      rerror ("sylv operation not supported");
    }

    rval = (*vfptr) (ent_data (e1), ent_data (e2), ent_data (e3));

    rent = ent_Create ();
    ent_data (rent) = rval;
    ent_type (rent) = rtype;

    ent_Clean (e1);
    ent_Clean (e2);
    ent_Clean (e3);
  }
  return (rent);
}

/* **************************************************************
 * Matrix determinant
 * ************************************************************** */

Ent *
Det (int nargs, Datum args[])
{
  void *rval;
  int rtype=UNDEF, type;
  void *(*vfptr) ();
  Ent *e, *rent;
  rent = 0;

  // Check n_args
  if (nargs != 1)
    rerror ("det: 1 argument allowed");

  e = bltin_get_ent (args[0]);
  type = ent_type (e);

  rtype = det_method[type].type;
  vfptr = (VFPTR) det_method[type].op;

  if (vfptr == 0)
  {
    fprintf (stderr, "Entity type: %s\n", etd (e));
    rerror ("det() operation not supported");
  }

  rval = (*vfptr) (ent_data (e));

  rent = ent_Create ();
  ent_data (rent) = rval;
  ent_type (rent) = rtype;

  ent_Clean (e);
  return (rent);
}

//
// rlabplus (C) Marijan Kostrun, 2005
// polynomial toolkit and hyperbolic functions

Ent *
ent_poly_eval_diff (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0;
  int zp, pc, i, k, n=1;
  MD *p=0;
  MDR *w=0;
  MDC *wc=0;

  if (nargs != 1 && nargs != 2)
  {
    fprintf (stdout, "Finds n-th derivative of polynomial.\n");
    fprintf (stdout, "Format:\n");
    fprintf (stdout, "  c = polyder(a/,n/),\n");
    fprintf (stdout, "where  'a'  is polynomial, and n>=0 the order of derivative.\n");
    rerror ("requires one or two arguments!\n");
  }

  e1 = bltin_get_ent (args[0]);
  if ((ent_type (e1) != MATRIX_DENSE_REAL) && (ent_type (e1) != MATRIX_DENSE_COMPLEX))
    rerror ("polyder: 'p' has to be a single-row matrix [a(N-1), .., a(1)]!\n");
  p = ent_data (e1);
  if (!EQVECT(p))
    rerror("polyder: 'p' has to be a single-row matrix [a(N-1), .., a(1)]!\n");
  pc = SIZE(p);

  if (nargs > 1)
  {
    e2 = bltin_get_ent (args[1]);
    if (ent_type (e2) == MATRIX_DENSE_REAL)
    {
      n = class_int(e2);
      if (n<0)
        n=0;
    }
  }

  // determine leading non-zero coefficient of the polynomial in p
  zp = 0;
  for (i=0; i<pc; i++)
  {
    if (cabs(mdcV0(p,i))!=0)
      break;
    zp++;
  }

  if (ent_type (e1) == MATRIX_DENSE_REAL)
  {
    if (!n)
    {
      w = mdr_Copy(p);
      goto _exit_polyder_real;
    }
    if (MD_TYPE_INT32(p))
    {
      //
      // we maintain integerness of polynomial
      //
      if (pc - zp < n + 1)
      {
        w = mdi_CreateScalar (0);
        goto _exit_polyder_real;
      }

      w = mdi_Create (1, pc - zp - n);
      mdr_Zero (w);
      for (i=1+zp; i<=pc - n; i++)
      {
        if (MdiV1 (p, i))
        {
          MdiV1 (w, i - zp) = MdiV1 (p, i);
          for (k=0; k<n;k++)
          {
            MdiV1 (w, i - zp) *= (pc - i - k);
          }
        }
      }
    }
    else
    {
      if (pc - zp < n + 1)
      {
        w = mdr_CreateScalar (0.0);
        goto _exit_polyder_real;
      }

      w = mdr_Create (1, pc - zp - n);
      mdr_Zero (w);
      for (i=1+zp; i<=pc - n; i++)
      {
        if (MdrV1 (p, i))
        {
          MdrV1 (w, i - zp) = MdrV1 (p, i);
          for (k=0; k<n;k++)
          {
            MdrV1 (w, i - zp) *= (pc - i - k);
          }
        }
      }
    }

_exit_polyder_real:

    ent_Clean (e1);
    ent_Clean (e2);
    return ent_Assign_Rlab_MDR(w);
  }

  if (!n)
  {
    wc = mdc_Copy(p);
    goto _exit_polyder_cmpl;
  }

  if (pc - zp < n + 1)
  {
    wc = mdc_CreateScalar (0.0,0.0);
    goto _exit_polyder_cmpl;
  }

  wc = mdc_Create (1, pc - zp - n);
  mdc_Zero (wc);
  for (i=1+zp; i<=pc - n; i++)
  {
    if (cabs(mdcV1 (p, i))>0)
    {
      MdcV1 (wc, i - zp) = mdcV1 (p, i);
      for (k=0; k<n;k++)
      {
        MdcV1 (wc, i - zp) *= (pc - i - k);
      }
    }
  }

_exit_polyder_cmpl:

  ent_Clean (e1);
  ent_Clean (e2);
  return ent_Assign_Rlab_MDC(wc);
}


Ent *
ent_poly_eval_int (int nargs, Datum args[])
{
  Ent *P=0;
  int pc, zp, i;
  MDR *p=0, *w=0;
  MDC *p2=0, *w2=0;

  if (nargs != 1)
  {
    fprintf (stdout, "Integration of a polynomial.\n");
    fprintf (stdout, "Format:\n");
    fprintf (stdout, "  c = polyint(a),\n");
    fprintf (stdout,
             "where  'c'  is a definite integral of 'a' with respect\n");
    fprintf (stdout,
             "to the argument. All polynomials are row-vectors where the\n");
    fprintf (stdout,
             "first non-zero entry is the higest non-vanishing power of\n");
    fprintf (stdout,
             "argument, down to the last entry, which is a constant term.\n");
    fprintf (stdout, "Constant of integration is zero.\n");
    rerror ("requires one argument!\n");
  }

  P = bltin_get_ent (args[0]);
  if (ent_type (P) != MATRIX_DENSE_REAL && ent_type (P) != MATRIX_DENSE_COMPLEX)
    rerror ("polyint: 'p' has to be a single-row matrix [a(N-1), .., a(1)]!\n");

  if (ent_type (P) == MATRIX_DENSE_REAL)
  {
    p = class_matrix_real (P);
    if (MNR (p) != 1)
      rerror
      ("polyint: 'p' has to be a single-row matrix [a(N-1), .., a(1)]!\n");
    pc = MNC (p);
    zp = 0;
    for (i = 1; i <= pc; i++)
    {
      if (Mdr1 (p, 1, i) != 0)
        break;
      zp++;
    }
    w = mdr_Create (1, pc - zp + 1);
    mdr_Zero (w);
    for (i = 1 + zp; i <= pc; i++)
      Mdr1 (w, 1, i - zp) = Mdr1 (p, i, 1) / (pc + 1 - i);
    Mdr1 (w, 1, pc + 1) = 0;

    ent_Clean (P);

    return ent_Assign_Rlab_MDR(w);
  }

  p2 = ent_data (P);
  if (MNR (p2) != 1)
    rerror
    ("polyint: 'p' has to be a single-row matrix [a(N-1), .., a(1)]!\n");
  pc = MNC (p2);
  zp = 0;
  for (i = 1; i <= pc; i++)
  {
    if (Mdc1r (p, 1, i) != 0 || Mdc1i (p, 1, i) != 0)
      break;
    zp++;
  }
  w2 = mdc_Create (1, pc - zp + 1);
  mdc_Zero (w2);
  for (i = 1 + zp; i <= pc; i++)
  {
    Mdc1r (w2, 1, i - zp) = Mdc1r (p, i, 1) / (pc + 1 - i);
    Mdc1i (w2, 1, i - zp) = Mdc1i (p, i, 1) / (pc + 1 - i);
  }
  Mdc1r (w2, 1, pc + 1) = 0;
  Mdc1i (w2, 1, pc + 1) = 0;

  ent_Clean (P);

  return ent_Assign_Rlab_MDC(w2);
}


Ent *
ent_poly_eval (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0;
  MD *x=0, *p=0;
  int xr, xc, zp, pc, i, j, k;
  double dummy;
  int idummy;
  Complex dummy2;

  if (nargs != 2)
  {
    fprintf (stdout, "polyval: Evaluation of polynomial.\n");
    fprintf (stdout, "polyval: Format:\n");
    fprintf (stdout, "polyval:   y = polyval(x,a),\n");
    fprintf (stdout,
             "polyval: where 'y' are the values of a polynomial 'a' for entries\n");
    fprintf (stdout,
             "polyval: of matrix 'x'.\n");
    rerror ("polyval: function requires two arguments!\n");
  }

  e1 = bltin_get_ent (args[0]);
  if ((ent_type (e1) != MATRIX_DENSE_REAL) && (ent_type (e1) != MATRIX_DENSE_COMPLEX))
    rerror ("polyval: 'x' has to be real or complex matrix!\n");
  x = ent_data (e1);
  if (SIZE(x)<1)
    rerror ("polyval: 'x' has to be real or complex matrix!\n");
  xr = MNR (x);
  xc = MNC (x);

  e2 = bltin_get_ent (args[1]);
  if ((ent_type (e2) != MATRIX_DENSE_REAL) && (ent_type (e2) != MATRIX_DENSE_COMPLEX))
    rerror ("polyval: " RLAB_ERROR_ARG2_VECTOR "\n");
  p = ent_data (e2);
  if (!EQVECT(p))
    rerror ("polyval: " RLAB_ERROR_ARG2_VECTOR "\n");
  pc = SIZE(p);

  // determine the highest non-zero coefficient in p
  zp = 0;
  for (i = 1; i <= pc; i++)
  {
    if (mdcV1 (p, i) != 0)
      break;

    zp++;
  }

  MDR *w=0;
  if (ent_type (e2) == MATRIX_DENSE_REAL && ent_type (e1) == MATRIX_DENSE_REAL)
  {
    if (x->type == RLAB_TYPE_INT32 && p->type == RLAB_TYPE_INT32)
    {
      w = mdi_Create (xr, xc);
      idummy = 0;
      for (i = 1; i <= xr; i++)
        for (j = 1; j <= xc; j++)
      {
        if (zp < pc)
        {
          idummy = MdiV1 (p, 1 + zp);
          for (k = 2 + zp; k <= pc; k++)
            idummy = Mdi1 (x, i, j) * idummy + MdiV1 (p, k);
        }
        Mdi1 (w, i, j) = idummy;
      }
    }
    else
    {
      w = mdr_Create (xr, xc);
      dummy = 0;
      for (i = 1; i <= xr; i++)
        for (j = 1; j <= xc; j++)
      {
        if (zp < pc)
        {
          dummy = mdrV1 (p, 1 + zp);
          for (k = 2 + zp; k <= pc; k++)
            dummy = mdr1 (x, i, j) * dummy + mdrV1 (p, k);
        }
        Mdr1 (w, i, j) = dummy;
      }
    }

    ent_Clean (e1);
    ent_Clean (e2);
    return ent_Assign_Rlab_MDR(w);
  }

  MDC * wc = mdc_Create (xr, xc);
  mdc_Zero (wc);
  for (i = 1; i <= xr; i++)
  {
    for (j = 1; j <= xc; j++)
    {
      dummy2 = mdcV1 (p, 1 + zp);
      for (k = 2 + zp; k <= pc; k++)
      {
        dummy2  = mdc1(x, i, j) * dummy2;
        dummy2 += mdcV1(p, k);
      }
      Mdc1(wc, i, j) = dummy2;
    }
  }

  ent_Clean (e1);
  ent_Clean (e2);
  return ent_Assign_Rlab_MDC(wc);
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
//
// H Y P E R B O L I C   F U N C T I O N S
//
//
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

// if all elements are less then 1 then say 'yes'
int do_mdr_1 (MDR * m)
{
  int i;
  for (i=0; i<SIZE(m); i++)
  {
    if ((mdrV0(m, i) > 1.0) || (mdrV0(m, i) < -1.0))
    {
      return (0);
    }
  }

  return (1);
}

// if all elements are non-negative then say 'yes'
int do_mdr_pos (MDR * m)
{
  int i;
  for (i=0; i<SIZE(m); i++)
  {
    if ((mdrV0(m, i) < 0.0))
    {
      return (0);
    }
  }

  return (1);
}


// Cosinus hyperbolicus */
MDR * mdr_mdr_cosh (MDR * x)
{
  int nr, nc, i;
  MDR *w;

  nr = MNR (x);
  nc = MNC (x);
  w = mdr_Create (nr, nc);
  mdr_Zero (w);

  for (i=0; i<nr*nc; i++)
    MdrV0(w, i) = cosh(mdrV0(x,i));

  return w;
}

MDC * mdc_mdc_cosh (MDC * x)
{
  int nr, nc, i;
  MDC *w;

  nr = MNR (x);
  nc = MNC (x);
  w = mdc_Create (nr, nc);
  mdc_Zero (w);

  for (i=0; i<nr*nc; i++)
    MdcV0(w,i) = ccosh(MdcV0(x,i));

  return w;
}

Ent *
ent_Cosh (int nargs, Datum args[])
{
  Ent *e1=0, *rent=0;

  e1 = bltin_get_ent (args[0]);

  if (ent_type (e1) == MATRIX_DENSE_REAL)
  {
    MDR *xr = ent_data (e1);
    rent = ent_Assign_Rlab_MDR( mdr_mdr_cosh(xr) );
  }
  else if (ent_type (e1) == MATRIX_DENSE_COMPLEX)
  {
    MDC *xc = ent_data (e1);
    rent = ent_Assign_Rlab_MDC( mdc_mdc_cosh(xc) );
  }
  else
    rerror("cosh: unknown argument type");

  ent_Clean(e1);

  return rent;
}

// sinus hyperbolicus */
MDR * mdr_mdr_sinh (MDR * x)
{
  int nr, nc, i;
  MDR *w;

  nr = MNR (x);
  nc = MNC (x);
  w = mdr_Create (nr, nc);
  mdr_Zero (w);

  for (i=0; i<nr*nc; i++)
    MdrV0 (w,i) = sinh(mdrV0(x,i));

  return w;
}

MDC * mdc_mdc_sinh (MDC * x)
{
  int nr, nc, i;
  MDC *w;

  nr = MNR (x);
  nc = MNC (x);
  w = mdc_Create (nr, nc);
  mdc_Zero (w);

  for (i=0; i<nr*nc; i++)
    MdcV0(w,i) = csinh(MdcV0(x,i));

  return w;
}

Ent *
ent_Sinh (int nargs, Datum args[])
{
  Ent *e1=0, *rent=0;

  e1 = bltin_get_ent (args[0]);

  if (ent_type (e1) == MATRIX_DENSE_REAL)
  {
    MDR *xr = ent_data (e1);
    rent = ent_Assign_Rlab_MDR( mdr_mdr_sinh(xr) );
  }
  else if (ent_type (e1) == MATRIX_DENSE_COMPLEX)
  {
    MDC *xc = ent_data (e1);
    rent = ent_Assign_Rlab_MDC( mdc_mdc_sinh(xc) );
  }
  else
    rerror("sinh: unknown argument type");

  ent_Clean(e1);

  return rent;
}

// tangens hyperbolicus */
MDR * mdr_mdr_tanh (MDR * x)
{
  int nr, nc, i;
  MDR *w;

  nr = MNR (x);
  nc = MNC (x);
  w = mdr_Create (nr, nc);
  mdr_Zero (w);

  for (i=0; i<nr*nc; i++)
    MdrV0(w,i) = tanh(mdrV0(x,i));

  return w;
}

MDC * mdc_mdc_tanh (MDC * x)
{
  int nr, nc, i;
  MDC *w;

  nr = MNR (x);
  nc = MNC (x);
  w = mdc_Create (nr, nc);
  mdc_Zero (w);

  for (i=0; i<nr*nc; i++)
    MdcV0(w,i) = ctanh(MdcV0(x, i));

  return w;
}

Ent *
ent_Tanh (int nargs, Datum args[])
{
  Ent *e1=0, *rent=0;

  e1 = bltin_get_ent (args[0]);

  if (ent_type (e1) == MATRIX_DENSE_REAL)
  {
    MDR *xr = ent_data (e1);
    rent = ent_Assign_Rlab_MDR( mdr_mdr_tanh(xr) );
  }
  else if (ent_type (e1) == MATRIX_DENSE_COMPLEX)
  {
    MDC *xc = ent_data (e1);
    rent = ent_Assign_Rlab_MDC( mdc_mdc_tanh(xc) );
  }
  else
    rerror("tanh: unknown argument type");

  ent_Clean(e1);

  return rent;
}

// arcus cosinus hyperbolicus */
MDR * mdr_mdr_acosh (MDR * x)
{
  int nr, nc, i;
  MDR *w;

  nr = MNR (x);
  nc = MNC (x);
  w = mdr_Create (nr, nc);
  mdr_Zero (w);

  for (i=0; i<nr*nc; i++)
    MdrV0 (w, i) = acosh(mdrV0(x,i));

  return w;
}

MDC * mdc_mdr_acosh (MDR * x)
{
  int nr, nc, i;
  MDC *w;

  nr = MNR (x);
  nc = MNC (x);
  w = mdc_Create (nr, nc);
  mdc_Zero (w);

  for (i=0; i<nr*nc; i++)
    MdcV0 (w, i) = cacosh(mdrV0(x, i));

  return w;
}

MDC * mdc_mdc_acosh (MDC * x)
{
  int nr, nc, i;
  MDC *w;

  nr = MNR (x);
  nc = MNC (x);
  w = mdc_Create (nr, nc);
  mdc_Zero (w);

  for (i=0; i<nr*nc; i++)
    MdcV0 (w, i) = cacosh(MdcV0(x, i));

  return w;
}

Ent *
ent_Acosh (int nargs, Datum args[])
{
  Ent *e1=0, *rent=0;

  e1 = bltin_get_ent (args[0]);

  if (ent_type (e1) == MATRIX_DENSE_REAL)
  {
    MDR *xr = ent_data (e1);
    if (do_mdr_pos(xr))
      rent = ent_Assign_Rlab_MDR( mdr_mdr_acosh(xr) );
    else
      rent = ent_Assign_Rlab_MDC( mdc_mdr_acosh(xr) );
  }
  else if (ent_type (e1) == MATRIX_DENSE_COMPLEX)
  {
    MDC *xc = ent_data (e1);
    rent = ent_Assign_Rlab_MDC( mdc_mdc_acosh(xc) );
  }
  else
    rerror("acosh: unknown argument type");

  ent_Clean(e1);

  return rent;
}

// arcus sinus hyperbolicus */
MDR * mdr_mdr_asinh(MDR * x)
{
  int nr, nc, i;
  MDR *w;

  nr = MNR (x);
  nc = MNC (x);
  w = mdr_Create (nr, nc);
  mdr_Zero (w);

  for (i=0; i<nr*nc; i++)
    MdrV0(w,i) = asinh(mdrV0(x,i));

  return w;
}

MDC * mdc_mdc_asinh(MDC * x)
{
  int nr, nc, i;
  MDC *w;

  nr = MNR (x);
  nc = MNC (x);
  w = mdc_Create (nr, nc);
  mdc_Zero (w);

  for (i=0; i<nr*nc; i++)
    MdcV0 (w, i) = casinh(MdcV0(x,i));

  return w;
}

Ent *
ent_Asinh (int nargs, Datum args[])
{
  Ent *e1=0, *rent=0;

  e1 = bltin_get_ent (args[0]);

  if (ent_type (e1) == MATRIX_DENSE_REAL)
  {
    MDR *xr = ent_data (e1);
    rent = ent_Assign_Rlab_MDR( mdr_mdr_asinh(xr) );
  }
  else if (ent_type (e1) == MATRIX_DENSE_COMPLEX)
  {
    MDC *xc = ent_data (e1);
    rent = ent_Assign_Rlab_MDC( mdc_mdc_asinh(xc) );
  }
  else
    rerror("asinh: unknown argument type");

  ent_Clean(e1);

  return rent;
}

// arcus tangens hyperbolicus */
MDR * mdr_mdr_atanh(MDR * x)
{
  int nr, nc, i;

  MDR *w;
  nr = MNR (x);
  nc = MNC (x);
  w = mdr_Create (nr, nc);
  mdr_Zero (w);

  for (i=0; i<nr*nc; i++)
    MdrV0(w,i) = atanh(mdrV0(x,i));

  return w;
}

MDC * mdc_mdr_atanh(MDR * x)
{
  int nr, nc, i;
  MDC *w;

  nr = MNR (x);
  nc = MNC (x);
  w = mdc_Create (nr, nc);
  mdc_Zero (w);

  for (i=0; i<nr*nc; i++)
    MdcV0(w,i) = catanh(mdrV0(x,i));

  return w;
}

MDC * mdc_mdc_atanh(MDC * x)
{
  int nr, nc, i;
  MDC *w;

  nr = MNR (x);
  nc = MNC (x);
  w = mdc_Create (nr, nc);
  mdc_Zero (w);

  for (i=0; i<nr*nc; i++)
    MdcV0(w,i) = catanh(MdcV0(x,i));

  return w;
}

#undef  THIS_SOLVER
#define THIS_SOLVER "atanh"
Ent * ent_Atanh (int nargs, Datum args[])
{
  Ent *e1=0, *rent=0;

  if (nargs != 1)
  {
    fprintf (stdout, THIS_SOLVER ": Inverse hyperbolic tangens.\n");
    fprintf (stdout, THIS_SOLVER ": Format:\n" RLAB_ERROR_ARG_1);
    fprintf (stdout, THIS_SOLVER ":   y = " THIS_SOLVER "(x),\n");
    rerror (THIS_SOLVER ": " );
  }

  e1 = bltin_get_ent (args[0]);

  if (ent_type (e1) == MATRIX_DENSE_REAL)
  {
    MDR *xr = ent_data (e1);
    if (do_mdr_1(xr))
      rent = ent_Assign_Rlab_MDR( mdr_mdr_atanh(xr) );
    else
      rent = ent_Assign_Rlab_MDC( mdc_mdr_atanh(xr) );
  }
  else if (ent_type (e1) == MATRIX_DENSE_COMPLEX)
  {
    MDC *xc = ent_data (e1);
    rent = ent_Assign_Rlab_MDC( mdc_mdc_atanh (xc) );
  }

  ent_Clean(e1);
  return rent;
}
