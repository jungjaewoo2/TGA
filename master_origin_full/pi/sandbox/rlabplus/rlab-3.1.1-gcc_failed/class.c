/* class.c */

/*
 * Provide the framework for class and class-method functionality.
 */

/* This file is a part of RLaB ("Our"-LaB)
   Copyright (C) 1995  Ian R. Searle
   Copyright (C) 2015  Marijan Kostrun

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
#include "util.h"
#include "class.h"
#include "mem.h"
#include "symbol.h"

/*
 * Header files from the available classes...
 */

#include "ent.h"
#include "btree.h"
#include "bltin.h"
#include "function.h"
#include "mdr.h"
#include "mdrf1.h"
#include "mdc.h"
#include "mdcf1.h"
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
#include "mde.h"
#include "rfileio.h"

#include "rlab_solver_parameters_names.h"

/* **************************************************************
 * Setup the class operation arrays...
 * ************************************************************** */

static OpDef copy_method[NCL];	/* How to Copy. */
static OpDef destroy_method[NCL];	/* How to destroy an object. */
static OpDef print_method[NCL];	/* How to print an object. */

static OpDef negate_method[NCL];	/* How to negate an object. */
static OpDef add_method[NCL][NCL];
static OpDef addto_method[NCL][NCL];
static OpDef subtractfrom_method[NCL][NCL];
static OpDef add_method_weight[NCL][NCL];
static OpDef subtract_method[NCL][NCL];
static OpDef subtract_method_weight[NCL][NCL];
static OpDef multiply_method[NCL][NCL];
static OpDef multiply_method_weight[NCL][NCL];


// <<val;wgt>> compatible methods
static OpDef el_multiply_method[NCL][NCL];
static OpDef elmultiplyby_method[NCL][NCL];
static OpDef elrdivideby_method[NCL][NCL];
static OpDef el_multiply_method_weight[NCL][NCL];
static OpDef el_rdivide_method[NCL][NCL];
static OpDef el_rdivide_method_weight[NCL][NCL];

static OpDef rdivide_method[NCL][NCL];
static OpDef ldivide_method[NCL][NCL];
static OpDef el_ldivide_method[NCL][NCL];
static OpDef power_method[NCL][NCL];
static OpDef el_power_method[NCL][NCL];
static OpDef el_power_method_weight[NCL][NCL];
static OpDef logical_scalar_method[NCL];
static OpDef sizeof_method[NCL];
static OpDef char_pointer_method[NCL];
static OpDef double_method[NCL];
static OpDef int_method[NCL];
static OpDef matrix_real_method[NCL];
static OpDef matrix_string_method[NCL];

static OpDef append_method[NCL][NCL];	/* Append ([ , ]) 2 objects. */
static OpDef stack_method[NCL][NCL];	/* Stack ([ ; ]) 2 objects. */

static OpDef matrix_sub_method[NCL];
static OpDef matrix_sub_r_method[NCL];
static OpDef matrix_sub_c_method[NCL];

static OpDef matrix_assign_method[NCL][NCL];
static OpDef matrix_assign_r_method[NCL][NCL];
static OpDef matrix_assign_c_method[NCL][NCL];

static OpDef matrix_vector_sub_method[NCL];
static OpDef matrix_assign_vector_method[NCL][NCL];

static OpDef ind_coerce_int_method[NCL][NCL];

static OpDef size_method[NCL];
static OpDef forloop_value_method[NCL];
static OpDef empty_method[NCL];
static OpDef member_ref_method[NCL];

static OpDef eq_method[NCL][NCL];
static OpDef ne_method[NCL][NCL];
static OpDef lt_method[NCL][NCL];
static OpDef le_method[NCL][NCL];
static OpDef gt_method[NCL][NCL];
static OpDef ge_method[NCL][NCL];
static OpDef and_method[NCL][NCL];
static OpDef or_method[NCL][NCL];
static OpDef not_method[NCL];

static OpDef transpose_method[NCL];
static OpDef nc_transpose_method[NCL];
static OpDef reshape_col_method[NCL];

static OpDef increment_method[NCL];
static OpDef decrement_method[NCL];

static OpDef attribute_method[NCL];

/* **************************************************************
 * Initialize the class-method lookup data.
 * ************************************************************** */

void
class_init (void)
{
  /*
   * Object copy.
   */

  copy_method[UNDEF].type = UNDEF;
  copy_method[UNDEF].op = (void *) undef_Copy;

  copy_method[BTREE].type = BTREE;
  copy_method[BTREE].op = (void *) btree_Copy;

  copy_method[MATRIX_DENSE_REAL].type = MATRIX_DENSE_REAL;
  copy_method[MATRIX_DENSE_REAL].op = (void *) mdr_Copy;

  copy_method[MATRIX_DENSE_COMPLEX].type = MATRIX_DENSE_COMPLEX;
  copy_method[MATRIX_DENSE_COMPLEX].op = (void *) mdc_Copy;

  copy_method[MATRIX_DENSE_STRING].type = MATRIX_DENSE_STRING;
  copy_method[MATRIX_DENSE_STRING].op = (void *) mds_Copy;

  copy_method[MATRIX_DENSE_ENTITY].type = MATRIX_DENSE_ENTITY;
  copy_method[MATRIX_DENSE_ENTITY].op = (void *) mde_Copy;

  copy_method[MATRIX_SPARSE_REAL].type = MATRIX_SPARSE_REAL;
  copy_method[MATRIX_SPARSE_REAL].op = (void *) msr_Copy;

  copy_method[MATRIX_SPARSE_COMPLEX].type = MATRIX_SPARSE_COMPLEX;
  copy_method[MATRIX_SPARSE_COMPLEX].op = (void *) msc_Copy;

  copy_method[U_FUNCTION].type = U_FUNCTION;
  copy_method[U_FUNCTION].op = (void *) function_Copy;
  

  /*
   * Object destruction.
   */

  destroy_method[UNDEF].type = -1;
  destroy_method[UNDEF].op = (void *) ent_undef_Destroy;

  destroy_method[DOUBLE].type = -1;
  destroy_method[DOUBLE].op = (void *) ent_double_Destroy;

  destroy_method[BTREE].type = -1;
  destroy_method[BTREE].op = (void *) btree_Destroy;

  destroy_method[U_FUNCTION].type = -1;
  destroy_method[U_FUNCTION].op = (void *) function_Destroy;

  destroy_method[U_CLASS].type = -1;
  destroy_method[U_CLASS].op = (void *) function_Destroy;

  destroy_method[BLTIN].type = -1;
  destroy_method[BLTIN].op = (void *) function_Destroy_bltin;

  destroy_method[MATRIX_DENSE_REAL].type = -1;
  destroy_method[MATRIX_DENSE_REAL].op = (void *) mdr_Destroy;

  destroy_method[MATRIX_SPARSE_REAL].type = -1;
  destroy_method[MATRIX_SPARSE_REAL].op = (void *) msr_Destroy;

  destroy_method[MATRIX_DENSE_COMPLEX].type = -1;
  destroy_method[MATRIX_DENSE_COMPLEX].op = (void *) mdc_Destroy;

  destroy_method[MATRIX_SPARSE_COMPLEX].type = -1;
  destroy_method[MATRIX_SPARSE_COMPLEX].op = (void *) msc_Destroy;

  destroy_method[MATRIX_DENSE_STRING].type = -1;
  destroy_method[MATRIX_DENSE_STRING].op = (void *) mds_Destroy;

  destroy_method[MATRIX_DENSE_ENTITY].type = -1;
  destroy_method[MATRIX_DENSE_ENTITY].op = (void *) mde_Destroy;

  /*
   * Object printing.
   */

  print_method[UNDEF].type = -1;
  print_method[UNDEF].op = (void *) ent_undef_Print;

  print_method[BTREE].type = -1;
  print_method[BTREE].op = (void *) btree_Print;

  print_method[U_FUNCTION].type = -1;
  print_method[U_FUNCTION].op = (void *) function_Print;

  print_method[U_CLASS].type = -1;
  print_method[U_CLASS].op = (void *) subprog_Print;

  print_method[BLTIN].type = -1;
  print_method[BLTIN].op = (void *) bltin_Print;

  print_method[MATRIX_DENSE_REAL].type = -1;
  print_method[MATRIX_DENSE_REAL].op = (void *) mdr_Print;

  print_method[MATRIX_DENSE_COMPLEX].type = -1;
  print_method[MATRIX_DENSE_COMPLEX].op = (void *) mdc_Print;

  print_method[MATRIX_DENSE_STRING].type = -1;
  print_method[MATRIX_DENSE_STRING].op = (void *) mds_Print;

  print_method[MATRIX_SPARSE_REAL].type = -1;
  print_method[MATRIX_SPARSE_REAL].op = (void *) msr_Print;

  print_method[MATRIX_SPARSE_COMPLEX].type = -1;
  print_method[MATRIX_SPARSE_COMPLEX].op = (void *) msc_Print;

  /*
   * Object negation.
   */

  negate_method[MATRIX_DENSE_REAL].type = MATRIX_DENSE_REAL;
  negate_method[MATRIX_DENSE_REAL].op = (void *) mdr_Negate;

  negate_method[MATRIX_DENSE_COMPLEX].type = MATRIX_DENSE_COMPLEX;
  negate_method[MATRIX_DENSE_COMPLEX].op = (void *) mdc_Negate;

  negate_method[MATRIX_SPARSE_REAL].type = MATRIX_SPARSE_REAL;
  negate_method[MATRIX_SPARSE_REAL].op = (void *) msr_Negate;

  negate_method[MATRIX_SPARSE_COMPLEX].type = MATRIX_SPARSE_COMPLEX;
  negate_method[MATRIX_SPARSE_COMPLEX].op = (void *) msc_Negate;

  /*
   * Object right-division.
   */

  rdivide_method[MATRIX_DENSE_REAL][MATRIX_DENSE_REAL].type = MATRIX_DENSE_REAL;
  rdivide_method[MATRIX_DENSE_REAL][MATRIX_DENSE_REAL].op =
    (void *) mdr_Rdivide;

  rdivide_method[MATRIX_DENSE_REAL][MATRIX_DENSE_COMPLEX].type =
    MATRIX_DENSE_COMPLEX;
  rdivide_method[MATRIX_DENSE_REAL][MATRIX_DENSE_COMPLEX].op =
    (void *) mdr_mdc_Rdivide;

  rdivide_method[MATRIX_DENSE_COMPLEX][MATRIX_DENSE_COMPLEX].type =
    MATRIX_DENSE_COMPLEX;
  rdivide_method[MATRIX_DENSE_COMPLEX][MATRIX_DENSE_COMPLEX].op =
    (void *) mdc_Rdivide;

  rdivide_method[MATRIX_DENSE_COMPLEX][MATRIX_DENSE_REAL].type =
    MATRIX_DENSE_COMPLEX;
  rdivide_method[MATRIX_DENSE_COMPLEX][MATRIX_DENSE_REAL].op =
    (void *) mdc_mdr_Rdivide;

  rdivide_method[MATRIX_SPARSE_REAL][MATRIX_SPARSE_REAL].type =
    MATRIX_DENSE_REAL;
  rdivide_method[MATRIX_SPARSE_REAL][MATRIX_SPARSE_REAL].op =
    (void *) msr_Rdivide;

  rdivide_method[MATRIX_DENSE_REAL][MATRIX_SPARSE_REAL].type =
    MATRIX_DENSE_REAL;
  rdivide_method[MATRIX_DENSE_REAL][MATRIX_SPARSE_REAL].op =
    (void *) mdr_msr_Rdivide;

  rdivide_method[MATRIX_SPARSE_REAL][MATRIX_DENSE_REAL].type =
    MATRIX_DENSE_REAL;
  rdivide_method[MATRIX_SPARSE_REAL][MATRIX_DENSE_REAL].op =
    (void *) msr_mdr_Rdivide;

  /*
   * Object element right-division.
   */

  el_rdivide_method[MATRIX_DENSE_REAL][MATRIX_DENSE_REAL].type = -1;
  el_rdivide_method[MATRIX_DENSE_REAL][MATRIX_DENSE_REAL].op =
    (void *) mdr_ElRdivide;

  el_rdivide_method[MATRIX_DENSE_REAL][MATRIX_DENSE_COMPLEX].type = -1;
  el_rdivide_method[MATRIX_DENSE_REAL][MATRIX_DENSE_COMPLEX].op =
    (void *) mdr_mdc_ElRdivide;

  el_rdivide_method[MATRIX_DENSE_COMPLEX][MATRIX_DENSE_COMPLEX].type = -1;
  el_rdivide_method[MATRIX_DENSE_COMPLEX][MATRIX_DENSE_COMPLEX].op =
    (void *) mdc_ElRdivide;

  el_rdivide_method[MATRIX_DENSE_COMPLEX][MATRIX_DENSE_REAL].type = -1;
  el_rdivide_method[MATRIX_DENSE_COMPLEX][MATRIX_DENSE_REAL].op =
    (void *) mdc_mdr_ElRdivide;

  el_rdivide_method[MATRIX_SPARSE_REAL][MATRIX_SPARSE_REAL].type = -1;
  el_rdivide_method[MATRIX_SPARSE_REAL][MATRIX_SPARSE_REAL].op =
    (void *) msr_ElRdivide;

  el_rdivide_method[MATRIX_SPARSE_COMPLEX][MATRIX_SPARSE_COMPLEX].type = -1;
  el_rdivide_method[MATRIX_SPARSE_COMPLEX][MATRIX_SPARSE_COMPLEX].op =
    (void *) msc_ElRdivide;

  el_rdivide_method[MATRIX_DENSE_REAL][MATRIX_SPARSE_REAL].type = -1;
  el_rdivide_method[MATRIX_DENSE_REAL][MATRIX_SPARSE_REAL].op =
    (void *) mdr_msr_ElRdivide;

  el_rdivide_method[MATRIX_SPARSE_REAL][MATRIX_DENSE_REAL].type = -1;
  el_rdivide_method[MATRIX_SPARSE_REAL][MATRIX_DENSE_REAL].op =
    (void *) msr_mdr_ElRdivide;

  el_rdivide_method[MATRIX_SPARSE_COMPLEX][MATRIX_DENSE_REAL].type = -1;
  el_rdivide_method[MATRIX_SPARSE_COMPLEX][MATRIX_DENSE_REAL].op =
    (void *) msc_mdr_ElRdivide;

  el_rdivide_method_weight[MATRIX_DENSE_REAL][MATRIX_DENSE_REAL].type = MATRIX_DENSE_REAL;
  el_rdivide_method_weight[MATRIX_DENSE_REAL][MATRIX_DENSE_REAL].op = (void *) mdr_ElRdivide_weight;


  /*
   * Object left-division.
   */

  ldivide_method[MATRIX_DENSE_REAL][MATRIX_DENSE_REAL].type = MATRIX_DENSE_REAL;
  ldivide_method[MATRIX_DENSE_REAL][MATRIX_DENSE_REAL].op =
    (void *) mdr_Ldivide;

  ldivide_method[MATRIX_DENSE_REAL][MATRIX_DENSE_COMPLEX].type =
    MATRIX_DENSE_COMPLEX;
  ldivide_method[MATRIX_DENSE_REAL][MATRIX_DENSE_COMPLEX].op =
    (void *) mdr_mdc_Ldivide;

  ldivide_method[MATRIX_DENSE_COMPLEX][MATRIX_DENSE_COMPLEX].type =
    MATRIX_DENSE_COMPLEX;
  ldivide_method[MATRIX_DENSE_COMPLEX][MATRIX_DENSE_COMPLEX].op =
    (void *) mdc_Ldivide;

  ldivide_method[MATRIX_DENSE_COMPLEX][MATRIX_DENSE_REAL].type =
    MATRIX_DENSE_COMPLEX;
  ldivide_method[MATRIX_DENSE_COMPLEX][MATRIX_DENSE_REAL].op =
    (void *) mdc_mdr_Ldivide;

  ldivide_method[MATRIX_SPARSE_REAL][MATRIX_DENSE_REAL].type =
    MATRIX_DENSE_REAL;
  ldivide_method[MATRIX_SPARSE_REAL][MATRIX_DENSE_REAL].op =
    (void *) msr_Ldivide;

  ldivide_method[MATRIX_SPARSE_REAL][MATRIX_SPARSE_REAL].type =
    MATRIX_DENSE_REAL;
  ldivide_method[MATRIX_SPARSE_REAL][MATRIX_SPARSE_REAL].op =
    (void *) msr_msr_Ldivide;

  /*
   * Object element left-division.
   */

  el_ldivide_method[MATRIX_DENSE_REAL][MATRIX_DENSE_REAL].type =
    MATRIX_DENSE_REAL;
  el_ldivide_method[MATRIX_DENSE_REAL][MATRIX_DENSE_REAL].op =
    (void *) mdr_ElLdivide;

  el_ldivide_method[MATRIX_DENSE_REAL][MATRIX_DENSE_COMPLEX].type =
    MATRIX_DENSE_COMPLEX;
  el_ldivide_method[MATRIX_DENSE_REAL][MATRIX_DENSE_COMPLEX].op =
    (void *) mdr_mdc_ElLdivide;

  el_ldivide_method[MATRIX_DENSE_COMPLEX][MATRIX_DENSE_COMPLEX].type =
    MATRIX_DENSE_COMPLEX;
  el_ldivide_method[MATRIX_DENSE_COMPLEX][MATRIX_DENSE_COMPLEX].op =
    (void *) mdc_ElLdivide;

  el_ldivide_method[MATRIX_DENSE_COMPLEX][MATRIX_DENSE_REAL].type =
    MATRIX_DENSE_COMPLEX;
  el_ldivide_method[MATRIX_DENSE_COMPLEX][MATRIX_DENSE_REAL].op =
    (void *) mdc_mdr_ElLdivide;

  /*
   * Object power operation.
   */

  power_method[MATRIX_DENSE_REAL][MATRIX_DENSE_REAL].type = MATRIX_DENSE_REAL;
  power_method[MATRIX_DENSE_REAL][MATRIX_DENSE_REAL].op = (void *) mdr_Power;

  power_method[MATRIX_DENSE_REAL][MATRIX_DENSE_COMPLEX].type =
    MATRIX_DENSE_COMPLEX;
  power_method[MATRIX_DENSE_REAL][MATRIX_DENSE_COMPLEX].op =
    (void *) mdr_mdc_Power;

  power_method[MATRIX_DENSE_COMPLEX][MATRIX_DENSE_COMPLEX].type =
    MATRIX_DENSE_COMPLEX;
  power_method[MATRIX_DENSE_COMPLEX][MATRIX_DENSE_COMPLEX].op =
    (void *) mdc_Power;

  power_method[MATRIX_DENSE_COMPLEX][MATRIX_DENSE_REAL].type =
    MATRIX_DENSE_COMPLEX;
  power_method[MATRIX_DENSE_COMPLEX][MATRIX_DENSE_REAL].op =
    (void *) mdc_mdr_Power;

  /*
   * Object element-power operation.
   */

  el_power_method[MATRIX_DENSE_REAL][MATRIX_DENSE_REAL].type =
    MATRIX_DENSE_REAL;
  el_power_method[MATRIX_DENSE_REAL][MATRIX_DENSE_REAL].op =
    (void *) mdr_ElPower;

  el_power_method[MATRIX_DENSE_REAL][MATRIX_DENSE_COMPLEX].type =
    MATRIX_DENSE_COMPLEX;
  el_power_method[MATRIX_DENSE_REAL][MATRIX_DENSE_COMPLEX].op =
    (void *) mdr_mdc_ElPower;

  el_power_method[MATRIX_DENSE_COMPLEX][MATRIX_DENSE_COMPLEX].type =
    MATRIX_DENSE_COMPLEX;
  el_power_method[MATRIX_DENSE_COMPLEX][MATRIX_DENSE_COMPLEX].op =
    (void *) mdc_ElPower;

  el_power_method[MATRIX_DENSE_COMPLEX][MATRIX_DENSE_REAL].type =
    MATRIX_DENSE_COMPLEX;
  el_power_method[MATRIX_DENSE_COMPLEX][MATRIX_DENSE_REAL].op =
    (void *) mdc_mdr_ElPower;

  // result can be real or imaginery
  el_power_method_weight[MATRIX_DENSE_REAL][MATRIX_DENSE_REAL].op =
      (void *) mdr_ElPower_weight;


  /*
   * Object multiply
   */
  // val:
  multiply_method[MATRIX_DENSE_REAL][MATRIX_DENSE_REAL].type = -1;
  multiply_method[MATRIX_DENSE_REAL][MATRIX_DENSE_REAL].op = (void *) mdr_Multiply;
  // wgt:
  multiply_method_weight[MATRIX_DENSE_REAL][MATRIX_DENSE_REAL].type = MATRIX_DENSE_REAL;
  multiply_method_weight[MATRIX_DENSE_REAL][MATRIX_DENSE_REAL].op = (void *) mdr_Multiply_weight;


  multiply_method[MATRIX_DENSE_REAL][MATRIX_DENSE_COMPLEX].type = -1;
  multiply_method[MATRIX_DENSE_REAL][MATRIX_DENSE_COMPLEX].op =
    (void *) mdr_mdc_Multiply;

  multiply_method[MATRIX_DENSE_COMPLEX][MATRIX_DENSE_COMPLEX].type = -1;
  multiply_method[MATRIX_DENSE_COMPLEX][MATRIX_DENSE_COMPLEX].op =
    (void *) mdc_Multiply;

  multiply_method[MATRIX_DENSE_COMPLEX][MATRIX_DENSE_REAL].type = -1;
  multiply_method[MATRIX_DENSE_COMPLEX][MATRIX_DENSE_REAL].op =
    (void *) mdc_mdr_Multiply;


  multiply_method[MATRIX_SPARSE_REAL][MATRIX_SPARSE_REAL].type = -1;
  multiply_method[MATRIX_SPARSE_REAL][MATRIX_SPARSE_REAL].op =
    (void *) msr_Multiply;

  multiply_method[MATRIX_DENSE_REAL][MATRIX_SPARSE_REAL].type = -1;
  multiply_method[MATRIX_DENSE_REAL][MATRIX_SPARSE_REAL].op =
    (void *) mdr_msr_Multiply;

  multiply_method[MATRIX_SPARSE_REAL][MATRIX_DENSE_REAL].type = -1;
  multiply_method[MATRIX_SPARSE_REAL][MATRIX_DENSE_REAL].op =
    (void *) msr_mdr_Multiply;

  multiply_method[MATRIX_SPARSE_COMPLEX][MATRIX_SPARSE_COMPLEX].type = -1;
  multiply_method[MATRIX_SPARSE_COMPLEX][MATRIX_SPARSE_COMPLEX].op =
    (void *) msc_Multiply;

  multiply_method[MATRIX_DENSE_COMPLEX][MATRIX_SPARSE_COMPLEX].type = -1;
  multiply_method[MATRIX_DENSE_COMPLEX][MATRIX_SPARSE_COMPLEX].op =
    (void *) mdc_msc_Multiply;

  multiply_method[MATRIX_SPARSE_COMPLEX][MATRIX_DENSE_COMPLEX].type = -1;
  multiply_method[MATRIX_SPARSE_COMPLEX][MATRIX_DENSE_COMPLEX].op =
    (void *) msc_mdc_Multiply;

  multiply_method[MATRIX_SPARSE_REAL][MATRIX_DENSE_COMPLEX].type = -1;
  multiply_method[MATRIX_SPARSE_REAL][MATRIX_DENSE_COMPLEX].op =
    (void *) msr_mdc_Multiply;

  /*
   * Object element-by-element multiply
   */
  el_multiply_method[MATRIX_DENSE_REAL][MATRIX_DENSE_REAL].type =
    MATRIX_DENSE_REAL;
  el_multiply_method[MATRIX_DENSE_REAL][MATRIX_DENSE_REAL].op =
    (void *) mdr_ElMultiply;

  el_multiply_method[MATRIX_DENSE_COMPLEX][MATRIX_DENSE_COMPLEX].type =
    MATRIX_DENSE_COMPLEX;
  el_multiply_method[MATRIX_DENSE_COMPLEX][MATRIX_DENSE_COMPLEX].op =
    (void *) mdc_ElMultiply;

  el_multiply_method[MATRIX_DENSE_REAL][MATRIX_DENSE_COMPLEX].type =
    MATRIX_DENSE_COMPLEX;
  el_multiply_method[MATRIX_DENSE_REAL][MATRIX_DENSE_COMPLEX].op =
    (void *) mdr_mdc_ElMultiply;

  el_multiply_method[MATRIX_DENSE_COMPLEX][MATRIX_DENSE_REAL].type =
    MATRIX_DENSE_COMPLEX;
  el_multiply_method[MATRIX_DENSE_COMPLEX][MATRIX_DENSE_REAL].op =
    (void *) mdc_mdr_ElMultiply;

  //
  // Object element-by-element multiply
  //
  el_multiply_method_weight[MATRIX_DENSE_REAL][MATRIX_DENSE_REAL].type = MATRIX_DENSE_REAL;
  el_multiply_method_weight[MATRIX_DENSE_REAL][MATRIX_DENSE_REAL].op = (void *) mdr_ElMultiply_weight;

  /*
   * Object subtraction.
   */
  subtract_method[MATRIX_DENSE_REAL][MATRIX_DENSE_REAL].type = MATRIX_DENSE_REAL;
  subtract_method[MATRIX_DENSE_REAL][MATRIX_DENSE_REAL].op = (void *) mdr_Subtract;
  subtract_method_weight[MATRIX_DENSE_REAL][MATRIX_DENSE_REAL].type = MATRIX_DENSE_REAL;
  subtract_method_weight[MATRIX_DENSE_REAL][MATRIX_DENSE_REAL].op = (void *) mdr_Add_weight;


  subtract_method[MATRIX_DENSE_REAL][MATRIX_DENSE_COMPLEX].type =
    MATRIX_DENSE_COMPLEX;
  subtract_method[MATRIX_DENSE_REAL][MATRIX_DENSE_COMPLEX].op =
    (void *) mdr_mdc_Subtract;

  subtract_method[MATRIX_DENSE_COMPLEX][MATRIX_DENSE_COMPLEX].type =
    MATRIX_DENSE_COMPLEX;
  subtract_method[MATRIX_DENSE_COMPLEX][MATRIX_DENSE_COMPLEX].op =
    (void *) mdc_Subtract;

  subtract_method[MATRIX_DENSE_COMPLEX][MATRIX_DENSE_REAL].type =
    MATRIX_DENSE_COMPLEX;
  subtract_method[MATRIX_DENSE_COMPLEX][MATRIX_DENSE_REAL].op =
    (void *) mdc_mdr_Subtract;

  subtract_method[MATRIX_SPARSE_REAL][MATRIX_SPARSE_REAL].type =
    MATRIX_SPARSE_REAL;
  subtract_method[MATRIX_SPARSE_REAL][MATRIX_SPARSE_REAL].op =
    (void *) msr_Subtract;

  subtract_method[MATRIX_SPARSE_REAL][MATRIX_DENSE_REAL].type =
    MATRIX_DENSE_REAL;
  subtract_method[MATRIX_SPARSE_REAL][MATRIX_DENSE_REAL].op =
    (void *) msr_mdr_Subtract;

  subtract_method[MATRIX_DENSE_REAL][MATRIX_SPARSE_REAL].type =
    MATRIX_DENSE_REAL;
  subtract_method[MATRIX_DENSE_REAL][MATRIX_SPARSE_REAL].op =
    (void *) mdr_msr_Subtract;

  subtract_method[MATRIX_SPARSE_COMPLEX][MATRIX_SPARSE_COMPLEX].type =
    MATRIX_SPARSE_COMPLEX;
  subtract_method[MATRIX_SPARSE_COMPLEX][MATRIX_SPARSE_COMPLEX].op =
    (void *) msc_Subtract;

  subtract_method[MATRIX_SPARSE_COMPLEX][MATRIX_DENSE_COMPLEX].type =
    MATRIX_DENSE_COMPLEX;
  subtract_method[MATRIX_SPARSE_COMPLEX][MATRIX_DENSE_COMPLEX].op =
    (void *) msc_mdc_Subtract;

  subtract_method[MATRIX_DENSE_COMPLEX][MATRIX_SPARSE_COMPLEX].type =
    MATRIX_DENSE_COMPLEX;
  subtract_method[MATRIX_DENSE_COMPLEX][MATRIX_SPARSE_COMPLEX].op =
    (void *) mdc_msc_Subtract;

  /*
   * Object 
   */
  // addition
  addto_method[MATRIX_DENSE_REAL][MATRIX_DENSE_REAL].op = mdr_AddTo;
  addto_method[MATRIX_DENSE_REAL][MATRIX_DENSE_REAL].type = MATRIX_DENSE_REAL;

  addto_method[MATRIX_DENSE_REAL][MATRIX_DENSE_COMPLEX].op = md_mdc_AddTo;
  addto_method[MATRIX_DENSE_REAL][MATRIX_DENSE_COMPLEX].type = MATRIX_DENSE_COMPLEX;

  addto_method[MATRIX_DENSE_COMPLEX][MATRIX_DENSE_REAL].op = md_mdc_AddTo;
  addto_method[MATRIX_DENSE_COMPLEX][MATRIX_DENSE_REAL].type = MATRIX_DENSE_COMPLEX;

  addto_method[MATRIX_DENSE_COMPLEX][MATRIX_DENSE_COMPLEX].op = md_mdc_AddTo;
  addto_method[MATRIX_DENSE_COMPLEX][MATRIX_DENSE_COMPLEX].type = MATRIX_DENSE_COMPLEX;

  // subtraction
  subtractfrom_method[MATRIX_DENSE_REAL][MATRIX_DENSE_REAL].op = mdr_SubtractFrom;
  subtractfrom_method[MATRIX_DENSE_REAL][MATRIX_DENSE_REAL].type = MATRIX_DENSE_REAL;

  subtractfrom_method[MATRIX_DENSE_REAL][MATRIX_DENSE_COMPLEX].op = md_mdc_SubtractFrom;
  subtractfrom_method[MATRIX_DENSE_REAL][MATRIX_DENSE_COMPLEX].type = MATRIX_DENSE_COMPLEX;

  subtractfrom_method[MATRIX_DENSE_COMPLEX][MATRIX_DENSE_REAL].op = md_mdc_SubtractFrom;
  subtractfrom_method[MATRIX_DENSE_COMPLEX][MATRIX_DENSE_REAL].type = MATRIX_DENSE_COMPLEX;

  subtractfrom_method[MATRIX_DENSE_COMPLEX][MATRIX_DENSE_COMPLEX].op = md_mdc_SubtractFrom;
  subtractfrom_method[MATRIX_DENSE_COMPLEX][MATRIX_DENSE_COMPLEX].type = MATRIX_DENSE_COMPLEX;

  // element multiply
  elmultiplyby_method[MATRIX_DENSE_REAL][MATRIX_DENSE_REAL].type = MATRIX_DENSE_REAL;
  elmultiplyby_method[MATRIX_DENSE_REAL][MATRIX_DENSE_REAL].op = (void *) mdr_ElMultiplyBy;

  elmultiplyby_method[MATRIX_DENSE_REAL][MATRIX_DENSE_COMPLEX].op = md_mdc_ElMultiplyBy;
  elmultiplyby_method[MATRIX_DENSE_REAL][MATRIX_DENSE_COMPLEX].type = MATRIX_DENSE_COMPLEX;

  elmultiplyby_method[MATRIX_DENSE_COMPLEX][MATRIX_DENSE_REAL].op = md_mdc_ElMultiplyBy;
  elmultiplyby_method[MATRIX_DENSE_COMPLEX][MATRIX_DENSE_REAL].type = MATRIX_DENSE_COMPLEX;

  elmultiplyby_method[MATRIX_DENSE_COMPLEX][MATRIX_DENSE_COMPLEX].op = md_mdc_ElMultiplyBy;
  elmultiplyby_method[MATRIX_DENSE_COMPLEX][MATRIX_DENSE_COMPLEX].type = MATRIX_DENSE_COMPLEX;

  // element right divide
  elrdivideby_method[MATRIX_DENSE_REAL][MATRIX_DENSE_REAL].type = MATRIX_DENSE_REAL;
  elrdivideby_method[MATRIX_DENSE_REAL][MATRIX_DENSE_REAL].op = (void *) mdr_ElRdivideBy;

  elrdivideby_method[MATRIX_DENSE_REAL][MATRIX_DENSE_COMPLEX].op = md_mdc_ElRdivideBy;
  elrdivideby_method[MATRIX_DENSE_REAL][MATRIX_DENSE_COMPLEX].type = MATRIX_DENSE_COMPLEX;

  elrdivideby_method[MATRIX_DENSE_COMPLEX][MATRIX_DENSE_REAL].op = md_mdc_ElRdivideBy;
  elrdivideby_method[MATRIX_DENSE_COMPLEX][MATRIX_DENSE_REAL].type = MATRIX_DENSE_COMPLEX;

  elrdivideby_method[MATRIX_DENSE_COMPLEX][MATRIX_DENSE_COMPLEX].op = md_mdc_ElRdivideBy;
  elrdivideby_method[MATRIX_DENSE_COMPLEX][MATRIX_DENSE_COMPLEX].type = MATRIX_DENSE_COMPLEX;

  /*
   * Object addition.
   */
  add_method[MATRIX_DENSE_REAL][MATRIX_DENSE_REAL].type = MATRIX_DENSE_REAL;
  add_method[MATRIX_DENSE_REAL][MATRIX_DENSE_REAL].op = (void *) mdr_Add;

  add_method_weight[MATRIX_DENSE_REAL][MATRIX_DENSE_REAL].type = MATRIX_DENSE_REAL;
  add_method_weight[MATRIX_DENSE_REAL][MATRIX_DENSE_REAL].op = (void *) mdr_Add_weight;


  add_method[MATRIX_DENSE_REAL][MATRIX_DENSE_COMPLEX].type =
    MATRIX_DENSE_COMPLEX;
  add_method[MATRIX_DENSE_REAL][MATRIX_DENSE_COMPLEX].op = (void *) mdr_mdc_Add;

  add_method[MATRIX_DENSE_COMPLEX][MATRIX_DENSE_REAL].type =
    MATRIX_DENSE_COMPLEX;
  add_method[MATRIX_DENSE_COMPLEX][MATRIX_DENSE_REAL].op = (void *) mdc_mdr_Add;

  add_method[MATRIX_DENSE_COMPLEX][MATRIX_DENSE_COMPLEX].type =
    MATRIX_DENSE_COMPLEX;
  add_method[MATRIX_DENSE_COMPLEX][MATRIX_DENSE_COMPLEX].op = (void *) mdc_Add;

  add_method[MATRIX_DENSE_STRING][MATRIX_DENSE_STRING].type =
    MATRIX_DENSE_STRING;
  add_method[MATRIX_DENSE_STRING][MATRIX_DENSE_STRING].op = (void *) mds_Add;

  add_method[MATRIX_SPARSE_REAL][MATRIX_SPARSE_REAL].type = MATRIX_SPARSE_REAL;
  add_method[MATRIX_SPARSE_REAL][MATRIX_SPARSE_REAL].op = (void *) msr_Add;

  add_method[MATRIX_DENSE_REAL][MATRIX_SPARSE_REAL].type = MATRIX_DENSE_REAL;
  add_method[MATRIX_DENSE_REAL][MATRIX_SPARSE_REAL].op = (void *) mdr_msr_Add;

  add_method[MATRIX_SPARSE_REAL][MATRIX_DENSE_REAL].type = MATRIX_DENSE_REAL;
  add_method[MATRIX_SPARSE_REAL][MATRIX_DENSE_REAL].op = (void *) msr_mdr_Add;

  add_method[MATRIX_SPARSE_COMPLEX][MATRIX_SPARSE_COMPLEX].type =
    MATRIX_SPARSE_COMPLEX;
  add_method[MATRIX_SPARSE_COMPLEX][MATRIX_SPARSE_COMPLEX].op =
    (void *) msc_Add;

  add_method[MATRIX_DENSE_COMPLEX][MATRIX_SPARSE_COMPLEX].type =
    MATRIX_DENSE_COMPLEX;
  add_method[MATRIX_DENSE_COMPLEX][MATRIX_SPARSE_COMPLEX].op =
    (void *) mdc_msc_Add;

  add_method[MATRIX_SPARSE_COMPLEX][MATRIX_DENSE_COMPLEX].type =
    MATRIX_DENSE_COMPLEX;
  add_method[MATRIX_SPARSE_COMPLEX][MATRIX_DENSE_COMPLEX].op =
    (void *) msc_mdc_Add;

  add_method[MATRIX_SPARSE_REAL][MATRIX_SPARSE_COMPLEX].type =
    MATRIX_SPARSE_COMPLEX;
  add_method[MATRIX_SPARSE_REAL][MATRIX_SPARSE_COMPLEX].op =
    (void *) msr_msc_Add;

  add_method[MATRIX_SPARSE_COMPLEX][MATRIX_SPARSE_REAL].type =
    MATRIX_SPARSE_COMPLEX;
  add_method[MATRIX_SPARSE_COMPLEX][MATRIX_SPARSE_REAL].op =
    (void *) msc_msr_Add;

  /*
   * Logical scalar value.
   */

  logical_scalar_method[MATRIX_DENSE_REAL].type = 0;
  logical_scalar_method[MATRIX_DENSE_REAL].op = (void *) mdr_LogicalScalar;

  logical_scalar_method[MATRIX_DENSE_COMPLEX].type = 0;
  logical_scalar_method[MATRIX_DENSE_COMPLEX].op = (void *) mdc_LogicalScalar;

  /*
   * sizeof ()
   * All return size_t
   */

  sizeof_method[DOUBLE].type = 0;
  sizeof_method[DOUBLE].op = (void *) double_Sizeof;

  sizeof_method[BLTIN].type = 0;
  sizeof_method[BLTIN].op = (void *) bltin_Sizeof;

  sizeof_method[U_FUNCTION].type = 0;
  sizeof_method[U_FUNCTION].op = (void *) function_Sizeof;

  sizeof_method[BTREE].type = 0;
  sizeof_method[BTREE].op = (void *) btree_Sizeof;

  sizeof_method[MATRIX_DENSE_REAL].type = 0;
  sizeof_method[MATRIX_DENSE_REAL].op = (void *) mdr_Sizeof;

  sizeof_method[MATRIX_DENSE_COMPLEX].type = 0;
  sizeof_method[MATRIX_DENSE_COMPLEX].op = (void *) mdc_Sizeof;

  sizeof_method[MATRIX_DENSE_STRING].type = 0;
  sizeof_method[MATRIX_DENSE_STRING].op = (void *) mds_Sizeof;

  sizeof_method[MATRIX_SPARSE_REAL].type = 0;
  sizeof_method[MATRIX_SPARSE_REAL].op = (void *) msr_Sizeof;

  sizeof_method[MATRIX_SPARSE_COMPLEX].type = 0;
  sizeof_method[MATRIX_SPARSE_COMPLEX].op = (void *) msc_Sizeof;

  /*
   * Character pointer.
   */

  char_pointer_method[MATRIX_DENSE_REAL].type = 0;
  char_pointer_method[MATRIX_DENSE_REAL].op = (void *) mdr_CharPointer;

  char_pointer_method[MATRIX_DENSE_STRING].type = 0;
  char_pointer_method[MATRIX_DENSE_STRING].op = (void *) mds_CharPointer;

 /*
   * Double value.
   */

  double_method[MATRIX_DENSE_REAL].type = 0;
  double_method[MATRIX_DENSE_REAL].op = (void *) mdr_Double;

  double_method[MATRIX_DENSE_COMPLEX].type = 0;
  double_method[MATRIX_DENSE_COMPLEX].op = (void *) mdc_Double;

  /*
   * Int value.
   */

  int_method[MATRIX_DENSE_REAL].type = 0;
  int_method[MATRIX_DENSE_REAL].op = (void *) mdr_Integer;

  /*
   * Matrix-Real value.
   */

  matrix_real_method[MATRIX_DENSE_REAL].type = 0;
  matrix_real_method[MATRIX_DENSE_REAL].op = (void *) mdr_MatrixReal;

  matrix_real_method[MATRIX_DENSE_COMPLEX].type = 0;
  matrix_real_method[MATRIX_DENSE_COMPLEX].op = (void *) mdc_MatrixReal;

  /*
   * Matrix-String value.
   */

  matrix_string_method[MATRIX_DENSE_STRING].type = 0;
  matrix_string_method[MATRIX_DENSE_STRING].op = (void *) mds_MatrixString;

  /*
   * Append two objects.
   * [o1, o2]
   */

  append_method[MATRIX_DENSE_REAL][MATRIX_DENSE_REAL].type = MATRIX_DENSE_REAL;
  append_method[MATRIX_DENSE_REAL][MATRIX_DENSE_REAL].op = (void *) mdr_Append;

  append_method[MATRIX_DENSE_REAL][MATRIX_DENSE_COMPLEX].type =
    MATRIX_DENSE_COMPLEX;
  append_method[MATRIX_DENSE_REAL][MATRIX_DENSE_COMPLEX].op =
    (void *) mdr_mdc_Append;

  append_method[MATRIX_DENSE_COMPLEX][MATRIX_DENSE_COMPLEX].type =
    MATRIX_DENSE_COMPLEX;
  append_method[MATRIX_DENSE_COMPLEX][MATRIX_DENSE_COMPLEX].op =
    (void *) mdc_Append;

  append_method[MATRIX_DENSE_COMPLEX][MATRIX_DENSE_REAL].type =
    MATRIX_DENSE_COMPLEX;
  append_method[MATRIX_DENSE_COMPLEX][MATRIX_DENSE_REAL].op =
    (void *) mdc_mdr_Append;

  append_method[MATRIX_DENSE_STRING][MATRIX_DENSE_STRING].type =
    MATRIX_DENSE_STRING;
  append_method[MATRIX_DENSE_STRING][MATRIX_DENSE_STRING].op =
    (void *) mds_Append;

  append_method[MATRIX_DENSE_STRING][MATRIX_DENSE_REAL].type =
    MATRIX_DENSE_STRING;
  append_method[MATRIX_DENSE_STRING][MATRIX_DENSE_REAL].op =
    (void *) mds_mdr_Append;

  append_method[MATRIX_DENSE_REAL][MATRIX_DENSE_STRING].type =
    MATRIX_DENSE_STRING;
  append_method[MATRIX_DENSE_REAL][MATRIX_DENSE_STRING].op =
    (void *) mdr_mds_Append;

  append_method[MATRIX_SPARSE_REAL][MATRIX_SPARSE_REAL].type =
    MATRIX_SPARSE_REAL;
  append_method[MATRIX_SPARSE_REAL][MATRIX_SPARSE_REAL].op =
    (void *) msr_Append;

  append_method[MATRIX_DENSE_REAL][MATRIX_SPARSE_REAL].type = MATRIX_DENSE_REAL;
  append_method[MATRIX_DENSE_REAL][MATRIX_SPARSE_REAL].op =
    (void *) mdr_msr_Append;

  append_method[MATRIX_SPARSE_REAL][MATRIX_DENSE_REAL].type = MATRIX_DENSE_REAL;
  append_method[MATRIX_SPARSE_REAL][MATRIX_DENSE_REAL].op =
    (void *) msr_mdr_Append;

  append_method[MATRIX_SPARSE_COMPLEX][MATRIX_SPARSE_COMPLEX].type =
    MATRIX_SPARSE_COMPLEX;
  append_method[MATRIX_SPARSE_COMPLEX][MATRIX_SPARSE_COMPLEX].op =
    (void *) msc_Append;

  append_method[MATRIX_DENSE_COMPLEX][MATRIX_SPARSE_COMPLEX].type =
    MATRIX_DENSE_COMPLEX;
  append_method[MATRIX_DENSE_COMPLEX][MATRIX_SPARSE_COMPLEX].op =
    (void *) mdc_msc_Append;

  append_method[MATRIX_SPARSE_COMPLEX][MATRIX_DENSE_COMPLEX].type =
    MATRIX_DENSE_COMPLEX;
  append_method[MATRIX_SPARSE_COMPLEX][MATRIX_DENSE_COMPLEX].op =
    (void *) msc_mdc_Append;

  /*
   * Stack two objects.
   */

  stack_method[MATRIX_DENSE_REAL][MATRIX_DENSE_REAL].type = MATRIX_DENSE_REAL;
  stack_method[MATRIX_DENSE_REAL][MATRIX_DENSE_REAL].op = (void *) mdr_Stack;

  stack_method[MATRIX_DENSE_REAL][MATRIX_DENSE_COMPLEX].type =
    MATRIX_DENSE_COMPLEX;
  stack_method[MATRIX_DENSE_REAL][MATRIX_DENSE_COMPLEX].op =
    (void *) mdr_mdc_Stack;

  stack_method[MATRIX_DENSE_COMPLEX][MATRIX_DENSE_COMPLEX].type =
    MATRIX_DENSE_COMPLEX;
  stack_method[MATRIX_DENSE_COMPLEX][MATRIX_DENSE_COMPLEX].op =
    (void *) mdc_Stack;

  stack_method[MATRIX_DENSE_COMPLEX][MATRIX_DENSE_REAL].type =
    MATRIX_DENSE_COMPLEX;
  stack_method[MATRIX_DENSE_COMPLEX][MATRIX_DENSE_REAL].op =
    (void *) mdc_mdr_Stack;

  stack_method[MATRIX_DENSE_STRING][MATRIX_DENSE_STRING].type =
    MATRIX_DENSE_STRING;
  stack_method[MATRIX_DENSE_STRING][MATRIX_DENSE_STRING].op =
    (void *) mds_Stack;

  stack_method[MATRIX_DENSE_REAL][MATRIX_DENSE_STRING].type =
    MATRIX_DENSE_STRING;
  stack_method[MATRIX_DENSE_REAL][MATRIX_DENSE_STRING].op =
    (void *) mdr_mds_Stack;

  stack_method[MATRIX_DENSE_STRING][MATRIX_DENSE_REAL].type =
    MATRIX_DENSE_STRING;
  stack_method[MATRIX_DENSE_STRING][MATRIX_DENSE_REAL].op =
    (void *) mds_mdr_Stack;

  stack_method[MATRIX_SPARSE_REAL][MATRIX_SPARSE_REAL].type =
    MATRIX_SPARSE_REAL;
  stack_method[MATRIX_SPARSE_REAL][MATRIX_SPARSE_REAL].op = (void *) msr_Stack;

  stack_method[MATRIX_DENSE_REAL][MATRIX_SPARSE_REAL].type = MATRIX_DENSE_REAL;
  stack_method[MATRIX_DENSE_REAL][MATRIX_SPARSE_REAL].op =
    (void *) mdr_msr_Stack;

  stack_method[MATRIX_SPARSE_REAL][MATRIX_DENSE_REAL].type = MATRIX_DENSE_REAL;
  stack_method[MATRIX_SPARSE_REAL][MATRIX_DENSE_REAL].op =
    (void *) msr_mdr_Stack;

  stack_method[MATRIX_SPARSE_COMPLEX][MATRIX_SPARSE_COMPLEX].type =
    MATRIX_SPARSE_COMPLEX;
  stack_method[MATRIX_SPARSE_COMPLEX][MATRIX_SPARSE_COMPLEX].op =
    (void *) msc_Stack;

  stack_method[MATRIX_DENSE_COMPLEX][MATRIX_SPARSE_COMPLEX].type =
    MATRIX_DENSE_COMPLEX;
  stack_method[MATRIX_DENSE_COMPLEX][MATRIX_SPARSE_COMPLEX].op =
    (void *) mdc_msc_Stack;

  stack_method[MATRIX_SPARSE_COMPLEX][MATRIX_DENSE_COMPLEX].type =
    MATRIX_DENSE_COMPLEX;
  stack_method[MATRIX_SPARSE_COMPLEX][MATRIX_DENSE_COMPLEX].op =
    (void *) msc_mdc_Stack;

  /*
   * Sub-Matrix Expression Evaluation.
   */

  matrix_sub_method[MATRIX_DENSE_REAL].type = MATRIX_DENSE_REAL;
  matrix_sub_method[MATRIX_DENSE_REAL].op = (void *) mdr_MatrixSub;

  matrix_sub_method[MATRIX_DENSE_COMPLEX].type = MATRIX_DENSE_COMPLEX;
  matrix_sub_method[MATRIX_DENSE_COMPLEX].op = (void *) mdc_MatrixSub;

  matrix_sub_method[MATRIX_DENSE_STRING].type = MATRIX_DENSE_STRING;
  matrix_sub_method[MATRIX_DENSE_STRING].op = (void *) mds_MatrixSub;

  matrix_sub_method[MATRIX_SPARSE_REAL].type = MATRIX_SPARSE_REAL;
  matrix_sub_method[MATRIX_SPARSE_REAL].op = (void *) msr_MatrixSub;

  matrix_sub_method[MATRIX_SPARSE_COMPLEX].type = MATRIX_SPARSE_COMPLEX;
  matrix_sub_method[MATRIX_SPARSE_COMPLEX].op = (void *) msc_MatrixSub;

  /*
   * Sub-Matrix Expression Evaluation.
   * Here only row indices are used.
   */

  matrix_sub_r_method[MATRIX_DENSE_REAL].type = MATRIX_DENSE_REAL;
  matrix_sub_r_method[MATRIX_DENSE_REAL].op = (void *) mdr_MatrixSubR;

  matrix_sub_r_method[MATRIX_DENSE_COMPLEX].type = MATRIX_DENSE_COMPLEX;
  matrix_sub_r_method[MATRIX_DENSE_COMPLEX].op = (void *) mdc_MatrixSubR;

  matrix_sub_r_method[MATRIX_DENSE_STRING].type = MATRIX_DENSE_STRING;
  matrix_sub_r_method[MATRIX_DENSE_STRING].op = (void *) mds_MatrixSubR;

  matrix_sub_r_method[MATRIX_SPARSE_REAL].type = MATRIX_SPARSE_REAL;
  matrix_sub_r_method[MATRIX_SPARSE_REAL].op = (void *) msr_MatrixSubR;

  matrix_sub_r_method[MATRIX_SPARSE_COMPLEX].type = MATRIX_SPARSE_COMPLEX;
  matrix_sub_r_method[MATRIX_SPARSE_COMPLEX].op = (void *) msc_MatrixSubR;

  /*
   * Sub-Matrix Expression Evaluation.
   * Here only column indices are used.
   */

  matrix_sub_c_method[MATRIX_DENSE_REAL].type = MATRIX_DENSE_REAL;
  matrix_sub_c_method[MATRIX_DENSE_REAL].op = (void *) mdr_MatrixSubC;

  matrix_sub_c_method[MATRIX_DENSE_COMPLEX].type = MATRIX_DENSE_COMPLEX;
  matrix_sub_c_method[MATRIX_DENSE_COMPLEX].op = (void *) mdc_MatrixSubC;

  matrix_sub_c_method[MATRIX_DENSE_STRING].type = MATRIX_DENSE_STRING;
  matrix_sub_c_method[MATRIX_DENSE_STRING].op = (void *) mds_MatrixSubC;

  matrix_sub_c_method[MATRIX_SPARSE_REAL].type = MATRIX_SPARSE_REAL;
  matrix_sub_c_method[MATRIX_SPARSE_REAL].op = (void *) msr_MatrixSubC;

  matrix_sub_c_method[MATRIX_SPARSE_COMPLEX].type = MATRIX_SPARSE_COMPLEX;
  matrix_sub_c_method[MATRIX_SPARSE_COMPLEX].op = (void *) msc_MatrixSubC;

  /*
   * Matrix assignement.
   * Uses both row and column indices.
   * matrix_assign_method[LHS][RHS]
   */

  matrix_assign_method[MATRIX_DENSE_REAL][MATRIX_DENSE_REAL].type = MATRIX_DENSE_REAL;
  matrix_assign_method[MATRIX_DENSE_REAL][MATRIX_DENSE_REAL].op = (void *) mdr_MatrixAssign;

  matrix_assign_method[UNDEF][MATRIX_DENSE_REAL].type = MATRIX_DENSE_REAL;
  matrix_assign_method[UNDEF][MATRIX_DENSE_REAL].op = (void *) mdr_MatrixAssign;

  matrix_assign_method[MATRIX_DENSE_REAL][MATRIX_DENSE_COMPLEX].type = MATRIX_DENSE_COMPLEX;
  matrix_assign_method[MATRIX_DENSE_REAL][MATRIX_DENSE_COMPLEX].op = (void *) mdr_mdc_MatrixAssign;

  matrix_assign_method[MATRIX_DENSE_COMPLEX][MATRIX_DENSE_COMPLEX].type = MATRIX_DENSE_COMPLEX;
  matrix_assign_method[MATRIX_DENSE_COMPLEX][MATRIX_DENSE_COMPLEX].op = (void *) mdc_MatrixAssign;

  matrix_assign_method[UNDEF][MATRIX_DENSE_COMPLEX].type = MATRIX_DENSE_COMPLEX;
  matrix_assign_method[UNDEF][MATRIX_DENSE_COMPLEX].op = (void *) mdc_MatrixAssign;

  matrix_assign_method[MATRIX_DENSE_COMPLEX][MATRIX_DENSE_REAL].type = MATRIX_DENSE_COMPLEX;
  matrix_assign_method[MATRIX_DENSE_COMPLEX][MATRIX_DENSE_REAL].op = (void *) mdc_mdr_MatrixAssign;

  matrix_assign_method[MATRIX_DENSE_STRING][MATRIX_DENSE_STRING].type = MATRIX_DENSE_STRING;
  matrix_assign_method[MATRIX_DENSE_STRING][MATRIX_DENSE_STRING].op = (void *) mds_MatrixAssign;

  matrix_assign_method[UNDEF][MATRIX_DENSE_STRING].type = MATRIX_DENSE_STRING;
  matrix_assign_method[UNDEF][MATRIX_DENSE_STRING].op =
    (void *) mds_MatrixAssign;

  matrix_assign_method[MATRIX_DENSE_REAL][MATRIX_DENSE_STRING].type = MATRIX_DENSE_STRING;
  matrix_assign_method[MATRIX_DENSE_REAL][MATRIX_DENSE_STRING].op = (void *) mdr_mds_MatrixAssign;

  matrix_assign_method[MATRIX_SPARSE_REAL][MATRIX_SPARSE_REAL].type = MATRIX_SPARSE_REAL;
  matrix_assign_method[MATRIX_SPARSE_REAL][MATRIX_SPARSE_REAL].op = (void *) msr_MatrixAssign;

  matrix_assign_method[MATRIX_SPARSE_REAL][MATRIX_DENSE_REAL].type = MATRIX_SPARSE_REAL;
  matrix_assign_method[MATRIX_SPARSE_REAL][MATRIX_DENSE_REAL].op = (void *) msr_mdr_MatrixAssign;

  matrix_assign_method[MATRIX_SPARSE_COMPLEX][MATRIX_SPARSE_COMPLEX].type = MATRIX_SPARSE_COMPLEX;
  matrix_assign_method[MATRIX_SPARSE_COMPLEX][MATRIX_SPARSE_COMPLEX].op = (void *) msc_MatrixAssign;

  matrix_assign_method[MATRIX_SPARSE_COMPLEX][MATRIX_DENSE_COMPLEX].type = MATRIX_SPARSE_COMPLEX;
  matrix_assign_method[MATRIX_SPARSE_COMPLEX][MATRIX_DENSE_COMPLEX].op = (void *) msc_mdc_MatrixAssign;

  matrix_assign_method[MATRIX_SPARSE_COMPLEX][MATRIX_SPARSE_REAL].type = MATRIX_SPARSE_COMPLEX;
  matrix_assign_method[MATRIX_SPARSE_COMPLEX][MATRIX_SPARSE_REAL].op = (void *) msc_msr_MatrixAssign;

  matrix_assign_method[MATRIX_SPARSE_COMPLEX][MATRIX_DENSE_REAL].type = MATRIX_SPARSE_COMPLEX;
  matrix_assign_method[MATRIX_SPARSE_COMPLEX][MATRIX_DENSE_REAL].op = (void *) msc_mdr_MatrixAssign;

  /*
   * Matrix assignement.
   * Uses row indices.
   * matrix_assign_r_method[LHS][RHS]
   */

  matrix_assign_r_method[MATRIX_DENSE_REAL][MATRIX_DENSE_REAL].type =
    MATRIX_DENSE_REAL;
  matrix_assign_r_method[MATRIX_DENSE_REAL][MATRIX_DENSE_REAL].op =
    (void *) mdr_MatrixAssignR;

  matrix_assign_r_method[UNDEF][MATRIX_DENSE_REAL].type = MATRIX_DENSE_REAL;
  matrix_assign_r_method[UNDEF][MATRIX_DENSE_REAL].op =
    (void *) mdr_MatrixAssignR;

  matrix_assign_r_method[MATRIX_DENSE_REAL][MATRIX_DENSE_COMPLEX].type =
    MATRIX_DENSE_COMPLEX;
  matrix_assign_r_method[MATRIX_DENSE_REAL][MATRIX_DENSE_COMPLEX].op =
    (void *) mdr_mdc_MatrixAssignR;

  matrix_assign_r_method[MATRIX_DENSE_COMPLEX][MATRIX_DENSE_COMPLEX].type =
    MATRIX_DENSE_COMPLEX;
  matrix_assign_r_method[MATRIX_DENSE_COMPLEX][MATRIX_DENSE_COMPLEX].op =
    (void *) mdc_MatrixAssignR;

  matrix_assign_r_method[UNDEF][MATRIX_DENSE_COMPLEX].type =
    MATRIX_DENSE_COMPLEX;
  matrix_assign_r_method[UNDEF][MATRIX_DENSE_COMPLEX].op =
    (void *) mdc_MatrixAssignR;

  matrix_assign_r_method[MATRIX_DENSE_COMPLEX][MATRIX_DENSE_REAL].type =
    MATRIX_DENSE_COMPLEX;
  matrix_assign_r_method[MATRIX_DENSE_COMPLEX][MATRIX_DENSE_REAL].op =
    (void *) mdc_mdr_MatrixAssignR;

  matrix_assign_r_method[MATRIX_DENSE_STRING][MATRIX_DENSE_STRING].type =
    MATRIX_DENSE_STRING;
  matrix_assign_r_method[MATRIX_DENSE_STRING][MATRIX_DENSE_STRING].op =
    (void *) mds_MatrixAssignR;

  matrix_assign_r_method[MATRIX_DENSE_REAL][MATRIX_DENSE_STRING].type =
    MATRIX_DENSE_STRING;
  matrix_assign_r_method[MATRIX_DENSE_REAL][MATRIX_DENSE_STRING].op =
    (void *) mdr_mds_MatrixAssignR;

  matrix_assign_r_method[UNDEF][MATRIX_DENSE_STRING].type = MATRIX_DENSE_STRING;
  matrix_assign_r_method[UNDEF][MATRIX_DENSE_STRING].op =
    (void *) mds_MatrixAssignR;

  /*
   * Matrix assignement.
   * Uses column indices.
   * matrix_assign_c_method[LHS][RHS]
   */

  matrix_assign_c_method[MATRIX_DENSE_REAL][MATRIX_DENSE_REAL].type =
    MATRIX_DENSE_REAL;
  matrix_assign_c_method[MATRIX_DENSE_REAL][MATRIX_DENSE_REAL].op =
    (void *) mdr_MatrixAssignC;

  matrix_assign_c_method[UNDEF][MATRIX_DENSE_REAL].type = MATRIX_DENSE_REAL;
  matrix_assign_c_method[UNDEF][MATRIX_DENSE_REAL].op =
    (void *) mdr_MatrixAssignC;

  matrix_assign_c_method[MATRIX_DENSE_REAL][MATRIX_DENSE_COMPLEX].type =
    MATRIX_DENSE_COMPLEX;
  matrix_assign_c_method[MATRIX_DENSE_REAL][MATRIX_DENSE_COMPLEX].op =
    (void *) mdr_mdc_MatrixAssignC;

  matrix_assign_c_method[MATRIX_DENSE_COMPLEX][MATRIX_DENSE_COMPLEX].type =
    MATRIX_DENSE_COMPLEX;
  matrix_assign_c_method[MATRIX_DENSE_COMPLEX][MATRIX_DENSE_COMPLEX].op =
    (void *) mdc_MatrixAssignC;

  matrix_assign_c_method[UNDEF][MATRIX_DENSE_COMPLEX].type =
    MATRIX_DENSE_COMPLEX;
  matrix_assign_c_method[UNDEF][MATRIX_DENSE_COMPLEX].op =
    (void *) mdc_MatrixAssignC;

  matrix_assign_c_method[MATRIX_DENSE_COMPLEX][MATRIX_DENSE_REAL].type =
    MATRIX_DENSE_COMPLEX;
  matrix_assign_c_method[MATRIX_DENSE_COMPLEX][MATRIX_DENSE_REAL].op =
    (void *) mdc_mdr_MatrixAssignC;

  matrix_assign_c_method[MATRIX_DENSE_STRING][MATRIX_DENSE_STRING].type =
    MATRIX_DENSE_STRING;
  matrix_assign_c_method[MATRIX_DENSE_STRING][MATRIX_DENSE_STRING].op =
    (void *) mds_MatrixAssignC;

  matrix_assign_c_method[MATRIX_DENSE_REAL][MATRIX_DENSE_STRING].type =
    MATRIX_DENSE_STRING;
  matrix_assign_c_method[MATRIX_DENSE_REAL][MATRIX_DENSE_STRING].op =
    (void *) mdr_mds_MatrixAssignC;

  matrix_assign_c_method[UNDEF][MATRIX_DENSE_STRING].type = MATRIX_DENSE_STRING;
  matrix_assign_c_method[UNDEF][MATRIX_DENSE_STRING].op =
    (void *) mds_MatrixAssignC;

  /*
   * Sub-Vector Expression Evaluation.
   */

  matrix_vector_sub_method[MATRIX_DENSE_REAL].type = MATRIX_DENSE_REAL;
  matrix_vector_sub_method[MATRIX_DENSE_REAL].op = (void *) mdr_VectorSub;

  matrix_vector_sub_method[MATRIX_DENSE_COMPLEX].type = MATRIX_DENSE_COMPLEX;
  matrix_vector_sub_method[MATRIX_DENSE_COMPLEX].op = (void *) mdc_VectorSub;

  matrix_vector_sub_method[MATRIX_DENSE_STRING].type = MATRIX_DENSE_STRING;
  matrix_vector_sub_method[MATRIX_DENSE_STRING].op = (void *) mds_VectorSub;

  matrix_vector_sub_method[MATRIX_DENSE_ENTITY].type = MATRIX_DENSE_ENTITY;
  matrix_vector_sub_method[MATRIX_DENSE_ENTITY].op = (void *) mde_VectorSub;

  matrix_vector_sub_method[MATRIX_SPARSE_REAL].type = MATRIX_SPARSE_REAL;
  matrix_vector_sub_method[MATRIX_SPARSE_REAL].op = (void *) msr_VectorSub;

  matrix_vector_sub_method[MATRIX_SPARSE_COMPLEX].type = MATRIX_SPARSE_COMPLEX;
  matrix_vector_sub_method[MATRIX_SPARSE_COMPLEX].op = (void *) msc_VectorSub;

  /*
   * Sub-Matrix assignment (vector style).
   * [LHS][RHS]
   */

  matrix_assign_vector_method[MATRIX_DENSE_REAL][MATRIX_DENSE_REAL].type =
    MATRIX_DENSE_REAL;
  matrix_assign_vector_method[MATRIX_DENSE_REAL][MATRIX_DENSE_REAL].op =
    (void *) mdr_VectorAssign;

  matrix_assign_vector_method[MATRIX_DENSE_REAL][MATRIX_DENSE_COMPLEX].type =
    MATRIX_DENSE_COMPLEX;
  matrix_assign_vector_method[MATRIX_DENSE_REAL][MATRIX_DENSE_COMPLEX].op =
    (void *) mdr_mdc_VectorAssign;

  matrix_assign_vector_method[UNDEF][MATRIX_DENSE_REAL].type =
    MATRIX_DENSE_REAL;
  matrix_assign_vector_method[UNDEF][MATRIX_DENSE_REAL].op =
    (void *) mdr_VectorAssign;

  matrix_assign_vector_method[MATRIX_DENSE_COMPLEX][MATRIX_DENSE_COMPLEX].type =
    MATRIX_DENSE_COMPLEX;
  matrix_assign_vector_method[MATRIX_DENSE_COMPLEX][MATRIX_DENSE_COMPLEX].op =
    (void *) mdc_VectorAssign;

  matrix_assign_vector_method[MATRIX_DENSE_COMPLEX][MATRIX_DENSE_REAL].type =
    MATRIX_DENSE_COMPLEX;
  matrix_assign_vector_method[MATRIX_DENSE_COMPLEX][MATRIX_DENSE_REAL].op =
    (void *) mdc_mdr_VectorAssign;

  matrix_assign_vector_method[UNDEF][MATRIX_DENSE_COMPLEX].type = MATRIX_DENSE_COMPLEX;
  matrix_assign_vector_method[UNDEF][MATRIX_DENSE_COMPLEX].op = (void *) mdc_VectorAssign;

  //
  // string matrices
  // 
  matrix_assign_vector_method[MATRIX_DENSE_STRING][MATRIX_DENSE_STRING].type = MATRIX_DENSE_STRING;
  matrix_assign_vector_method[MATRIX_DENSE_STRING][MATRIX_DENSE_STRING].op = (void *) mds_VectorAssign;

  matrix_assign_vector_method[UNDEF][MATRIX_DENSE_STRING].type = MATRIX_DENSE_STRING;
  matrix_assign_vector_method[UNDEF][MATRIX_DENSE_STRING].op = (void *) mds_VectorAssign;

  matrix_assign_vector_method[MATRIX_DENSE_REAL][MATRIX_DENSE_STRING].type = MATRIX_DENSE_STRING;
  matrix_assign_vector_method[MATRIX_DENSE_REAL][MATRIX_DENSE_STRING].op = (void *) mdr_mds_VectorAssign;

  //
  // sparse matrices
  //
  matrix_assign_vector_method[MATRIX_SPARSE_REAL][MATRIX_SPARSE_REAL].type =
    MATRIX_SPARSE_REAL;
  matrix_assign_vector_method[MATRIX_SPARSE_REAL][MATRIX_SPARSE_REAL].op =
    (void *) msr_VectorAssign;

  matrix_assign_vector_method[MATRIX_SPARSE_REAL][MATRIX_DENSE_REAL].type =
    MATRIX_SPARSE_REAL;
  matrix_assign_vector_method[MATRIX_SPARSE_REAL][MATRIX_DENSE_REAL].op =
    (void *) msr_mdr_VectorAssign;

  matrix_assign_vector_method[MATRIX_SPARSE_COMPLEX][MATRIX_SPARSE_COMPLEX].
    type = MATRIX_SPARSE_COMPLEX;
  matrix_assign_vector_method[MATRIX_SPARSE_COMPLEX][MATRIX_SPARSE_COMPLEX].op =
    (void *) msc_VectorAssign;

  matrix_assign_vector_method[MATRIX_SPARSE_COMPLEX][MATRIX_DENSE_COMPLEX].
    type = MATRIX_SPARSE_COMPLEX;
  matrix_assign_vector_method[MATRIX_SPARSE_COMPLEX][MATRIX_DENSE_COMPLEX].op =
    (void *) msc_mdc_VectorAssign;

  matrix_assign_vector_method[MATRIX_SPARSE_COMPLEX][MATRIX_SPARSE_REAL].type =
    MATRIX_SPARSE_COMPLEX;
  matrix_assign_vector_method[MATRIX_SPARSE_COMPLEX][MATRIX_SPARSE_REAL].op =
    (void *) msc_msr_VectorAssign;

  matrix_assign_vector_method[MATRIX_SPARSE_COMPLEX][MATRIX_DENSE_REAL].type =
    MATRIX_SPARSE_COMPLEX;
  matrix_assign_vector_method[MATRIX_SPARSE_COMPLEX][MATRIX_DENSE_REAL].op =
    (void *) msc_mdr_VectorAssign;

  /*
   * Index coercion into int.
   * ind_coerce_int_method[LHS][INDEX]
   */

  ind_coerce_int_method[MATRIX_DENSE_REAL][MATRIX_DENSE_REAL].type = -1;
  ind_coerce_int_method[MATRIX_DENSE_REAL][MATRIX_DENSE_REAL].op =
    (void *) mdr_IndCoerceInt;

  ind_coerce_int_method[UNDEF][MATRIX_DENSE_REAL].type = -1;
  ind_coerce_int_method[UNDEF][MATRIX_DENSE_REAL].op = (void *) mdr_IndCoerceInt;

  ind_coerce_int_method[MATRIX_DENSE_COMPLEX][MATRIX_DENSE_REAL].type = -1;
  ind_coerce_int_method[MATRIX_DENSE_COMPLEX][MATRIX_DENSE_REAL].op = (void *) mdc_IndCoerceInt;

  ind_coerce_int_method[MATRIX_DENSE_STRING][MATRIX_DENSE_REAL].type = -1;
  ind_coerce_int_method[MATRIX_DENSE_STRING][MATRIX_DENSE_REAL].op = (void *) mds_IndCoerceInt;

  ind_coerce_int_method[MATRIX_DENSE_ENTITY][MATRIX_DENSE_REAL].type = -1;
  ind_coerce_int_method[MATRIX_DENSE_ENTITY][MATRIX_DENSE_REAL].op = (void *) mde_IndCoerceInt;

  ind_coerce_int_method[MATRIX_SPARSE_REAL][MATRIX_DENSE_REAL].type = -1;
  ind_coerce_int_method[MATRIX_SPARSE_REAL][MATRIX_DENSE_REAL].op = (void *) msr_IndCoerceInt;

  ind_coerce_int_method[MATRIX_SPARSE_COMPLEX][MATRIX_DENSE_REAL].type = -1;
  ind_coerce_int_method[MATRIX_SPARSE_COMPLEX][MATRIX_DENSE_REAL].op = (void *) msc_IndCoerceInt;

  /*
   * The size of an object.
   */

  size_method[MATRIX_DENSE_REAL].type = 0;
  size_method[MATRIX_DENSE_REAL].op = (void *) mdr_Size;

  size_method[MATRIX_DENSE_COMPLEX].type = 0;
  size_method[MATRIX_DENSE_COMPLEX].op = (void *) mdc_Size;

  size_method[MATRIX_DENSE_STRING].type = 0;
  size_method[MATRIX_DENSE_STRING].op = (void *) mds_Size;

  size_method[BTREE].type = 0;
  size_method[BTREE].op = (void *) btree_CountNodes;

  size_method[MATRIX_DENSE_ENTITY].type = 0;
  size_method[MATRIX_DENSE_ENTITY].op = (void *) mde_Size;


  /*
   * The "thing" a forloop needs for the index variable.
   */

  forloop_value_method[MATRIX_DENSE_REAL].type = MATRIX_DENSE_REAL;
  forloop_value_method[MATRIX_DENSE_REAL].op = (void *) mdr_ForLoopValue;

  forloop_value_method[MATRIX_DENSE_COMPLEX].type = MATRIX_DENSE_COMPLEX;
  forloop_value_method[MATRIX_DENSE_COMPLEX].op = (void *) mdc_ForLoopValue;

  forloop_value_method[MATRIX_DENSE_STRING].type = MATRIX_DENSE_STRING;
  forloop_value_method[MATRIX_DENSE_STRING].op = (void *) mds_ForLoopValue;

  /*
   * Create an empty something.
   */

  empty_method[MATRIX_DENSE_REAL].type = MATRIX_DENSE_REAL;
  empty_method[MATRIX_DENSE_REAL].op = (void *) mdr_Empty;

  /*
   * Return the result of a member reference.
   * At present, I cannot concieve how these functions
   * could always return the same type. If that is the
   * case, then they will have to return entities, which
   * contain the proper type indicator.
   */

  member_ref_method[MATRIX_DENSE_ENTITY].type = -1;
  member_ref_method[MATRIX_DENSE_ENTITY].op = (void *) mde_MemberRef;

  member_ref_method[BTREE].type = VAR;
  member_ref_method[BTREE].op = (void *) btree_MemberRef;

  member_ref_method[MATRIX_DENSE_REAL].type = -1;
  member_ref_method[MATRIX_DENSE_REAL].op = (void *) mdr_MemberRef;

  member_ref_method[MATRIX_SPARSE_REAL].type = -1;
  member_ref_method[MATRIX_SPARSE_REAL].op = (void *) msr_MemberRef;

  member_ref_method[MATRIX_DENSE_COMPLEX].type = -1;
  member_ref_method[MATRIX_DENSE_COMPLEX].op = (void *) mdc_MemberRef;

  member_ref_method[MATRIX_SPARSE_COMPLEX].type = -1;
  member_ref_method[MATRIX_SPARSE_COMPLEX].op = (void *) msc_MemberRef;

  member_ref_method[MATRIX_DENSE_STRING].type = -1;
  member_ref_method[MATRIX_DENSE_STRING].op = (void *) mds_MemberRef;

  member_ref_method[U_FUNCTION].type = -1;
  member_ref_method[U_FUNCTION].op = (void *) function_MemberRef;

  member_ref_method[BLTIN].type = -1;
  member_ref_method[BLTIN].op = (void *) bltin_MemberRef;

  /*
   * Compare (==) two objects.
   */

  eq_method[MATRIX_DENSE_REAL][MATRIX_DENSE_REAL].type = MATRIX_DENSE_REAL;
  eq_method[MATRIX_DENSE_REAL][MATRIX_DENSE_REAL].op = (void *) mdr_Eq;

  eq_method[MATRIX_DENSE_COMPLEX][MATRIX_DENSE_COMPLEX].type =
    MATRIX_DENSE_REAL;
  eq_method[MATRIX_DENSE_COMPLEX][MATRIX_DENSE_COMPLEX].op = (void *) mdc_Eq;

  eq_method[MATRIX_DENSE_REAL][MATRIX_DENSE_COMPLEX].type = MATRIX_DENSE_REAL;
  eq_method[MATRIX_DENSE_REAL][MATRIX_DENSE_COMPLEX].op = (void *) mdr_mdc_Eq;

  eq_method[MATRIX_DENSE_COMPLEX][MATRIX_DENSE_REAL].type = MATRIX_DENSE_REAL;
  eq_method[MATRIX_DENSE_COMPLEX][MATRIX_DENSE_REAL].op = (void *) mdc_mdr_Eq;

  eq_method[MATRIX_DENSE_STRING][MATRIX_DENSE_STRING].type = MATRIX_DENSE_REAL;
  eq_method[MATRIX_DENSE_STRING][MATRIX_DENSE_STRING].op = (void *) mds_Eq;

  eq_method[MATRIX_DENSE_REAL][MATRIX_DENSE_STRING].type = MATRIX_DENSE_REAL;
  eq_method[MATRIX_DENSE_REAL][MATRIX_DENSE_STRING].op = (void *) mdr_mds_Eq;

  eq_method[MATRIX_DENSE_STRING][MATRIX_DENSE_REAL].type = MATRIX_DENSE_REAL;
  eq_method[MATRIX_DENSE_STRING][MATRIX_DENSE_REAL].op = (void *) mds_mdr_Eq;

  eq_method[MATRIX_DENSE_COMPLEX][MATRIX_DENSE_STRING].type = MATRIX_DENSE_REAL;
  eq_method[MATRIX_DENSE_COMPLEX][MATRIX_DENSE_STRING].op = (void *) mdr_mds_Eq;

  eq_method[MATRIX_DENSE_STRING][MATRIX_DENSE_COMPLEX].type = MATRIX_DENSE_REAL;
  eq_method[MATRIX_DENSE_STRING][MATRIX_DENSE_COMPLEX].op = (void *) mds_mdr_Eq;

  eq_method[MATRIX_SPARSE_REAL][MATRIX_SPARSE_REAL].type = MATRIX_SPARSE_REAL;
  eq_method[MATRIX_SPARSE_REAL][MATRIX_SPARSE_REAL].op = (void *) msr_Eq;

  eq_method[MATRIX_DENSE_REAL][MATRIX_SPARSE_REAL].type = MATRIX_SPARSE_REAL;
  eq_method[MATRIX_DENSE_REAL][MATRIX_SPARSE_REAL].op = (void *) mdr_msr_Eq;

  eq_method[MATRIX_SPARSE_REAL][MATRIX_DENSE_REAL].type = MATRIX_SPARSE_REAL;
  eq_method[MATRIX_SPARSE_REAL][MATRIX_DENSE_REAL].op = (void *) msr_mdr_Eq;

  eq_method[MATRIX_SPARSE_COMPLEX][MATRIX_SPARSE_COMPLEX].type =
    MATRIX_SPARSE_REAL;
  eq_method[MATRIX_SPARSE_COMPLEX][MATRIX_SPARSE_COMPLEX].op = (void *) msc_Eq;

  eq_method[MATRIX_DENSE_COMPLEX][MATRIX_SPARSE_COMPLEX].type =
    MATRIX_SPARSE_REAL;
  eq_method[MATRIX_DENSE_COMPLEX][MATRIX_SPARSE_COMPLEX].op =
    (void *) mdc_msc_Eq;

  eq_method[MATRIX_SPARSE_COMPLEX][MATRIX_DENSE_COMPLEX].type =
    MATRIX_SPARSE_REAL;
  eq_method[MATRIX_SPARSE_COMPLEX][MATRIX_DENSE_COMPLEX].op =
    (void *) msc_mdc_Eq;

  eq_method[MATRIX_SPARSE_COMPLEX][MATRIX_DENSE_REAL].type = MATRIX_SPARSE_REAL;
  eq_method[MATRIX_SPARSE_COMPLEX][MATRIX_DENSE_REAL].op = (void *) msc_mdr_Eq;

  eq_method[MATRIX_DENSE_REAL][MATRIX_SPARSE_COMPLEX].type = MATRIX_SPARSE_REAL;
  eq_method[MATRIX_DENSE_REAL][MATRIX_SPARSE_COMPLEX].op = (void *) mdr_msc_Eq;

  eq_method[MATRIX_SPARSE_COMPLEX][MATRIX_SPARSE_REAL].type =
    MATRIX_SPARSE_REAL;
  eq_method[MATRIX_SPARSE_COMPLEX][MATRIX_SPARSE_REAL].op = (void *) msc_msr_Eq;

  eq_method[MATRIX_SPARSE_REAL][MATRIX_SPARSE_COMPLEX].type =
    MATRIX_SPARSE_REAL;
  eq_method[MATRIX_SPARSE_REAL][MATRIX_SPARSE_COMPLEX].op = (void *) msr_msc_Eq;


  /*
   * Compare (!=) two objects.
   */

  ne_method[MATRIX_DENSE_REAL][MATRIX_DENSE_REAL].type = MATRIX_DENSE_REAL;
  ne_method[MATRIX_DENSE_REAL][MATRIX_DENSE_REAL].op = (void *) mdr_Ne;

  ne_method[MATRIX_DENSE_COMPLEX][MATRIX_DENSE_COMPLEX].type =
    MATRIX_DENSE_REAL;
  ne_method[MATRIX_DENSE_COMPLEX][MATRIX_DENSE_COMPLEX].op = (void *) mdc_Ne;

  ne_method[MATRIX_DENSE_REAL][MATRIX_DENSE_COMPLEX].type = MATRIX_DENSE_REAL;
  ne_method[MATRIX_DENSE_REAL][MATRIX_DENSE_COMPLEX].op = (void *) mdr_mdc_Ne;

  ne_method[MATRIX_DENSE_COMPLEX][MATRIX_DENSE_REAL].type = MATRIX_DENSE_REAL;
  ne_method[MATRIX_DENSE_COMPLEX][MATRIX_DENSE_REAL].op = (void *) mdc_mdr_Ne;

  ne_method[MATRIX_DENSE_STRING][MATRIX_DENSE_STRING].type = MATRIX_DENSE_REAL;
  ne_method[MATRIX_DENSE_STRING][MATRIX_DENSE_STRING].op = (void *) mds_Ne;

  ne_method[MATRIX_DENSE_REAL][MATRIX_DENSE_STRING].type = MATRIX_DENSE_REAL;
  ne_method[MATRIX_DENSE_REAL][MATRIX_DENSE_STRING].op = (void *) mdr_mds_Ne;

  ne_method[MATRIX_DENSE_STRING][MATRIX_DENSE_REAL].type = MATRIX_DENSE_REAL;
  ne_method[MATRIX_DENSE_STRING][MATRIX_DENSE_REAL].op = (void *) mds_mdr_Ne;

  ne_method[MATRIX_DENSE_COMPLEX][MATRIX_DENSE_STRING].type = MATRIX_DENSE_REAL;
  ne_method[MATRIX_DENSE_COMPLEX][MATRIX_DENSE_STRING].op = (void *) mdr_mds_Ne;

  ne_method[MATRIX_DENSE_STRING][MATRIX_DENSE_COMPLEX].type = MATRIX_DENSE_REAL;
  ne_method[MATRIX_DENSE_STRING][MATRIX_DENSE_COMPLEX].op = (void *) mds_mdr_Ne;

  ne_method[MATRIX_SPARSE_REAL][MATRIX_SPARSE_REAL].type = MATRIX_SPARSE_REAL;
  ne_method[MATRIX_SPARSE_REAL][MATRIX_SPARSE_REAL].op = (void *) msr_Ne;

  ne_method[MATRIX_DENSE_REAL][MATRIX_SPARSE_REAL].type = MATRIX_SPARSE_REAL;
  ne_method[MATRIX_DENSE_REAL][MATRIX_SPARSE_REAL].op = (void *) mdr_msr_Ne;

  ne_method[MATRIX_SPARSE_REAL][MATRIX_DENSE_REAL].type = MATRIX_SPARSE_REAL;
  ne_method[MATRIX_SPARSE_REAL][MATRIX_DENSE_REAL].op = (void *) msr_mdr_Ne;

  ne_method[MATRIX_SPARSE_COMPLEX][MATRIX_SPARSE_COMPLEX].type =
    MATRIX_SPARSE_REAL;
  ne_method[MATRIX_SPARSE_COMPLEX][MATRIX_SPARSE_COMPLEX].op = (void *) msc_Ne;

  ne_method[MATRIX_SPARSE_COMPLEX][MATRIX_DENSE_COMPLEX].type =
    MATRIX_SPARSE_REAL;
  ne_method[MATRIX_SPARSE_COMPLEX][MATRIX_DENSE_COMPLEX].op =
    (void *) msc_mdc_Ne;

  ne_method[MATRIX_DENSE_COMPLEX][MATRIX_SPARSE_COMPLEX].type =
    MATRIX_SPARSE_REAL;
  ne_method[MATRIX_DENSE_COMPLEX][MATRIX_SPARSE_COMPLEX].op =
    (void *) mdc_msc_Ne;

  /*
   * Compare (<) two objects.
   */

  lt_method[MATRIX_DENSE_REAL][MATRIX_DENSE_REAL].type = MATRIX_DENSE_REAL;
  lt_method[MATRIX_DENSE_REAL][MATRIX_DENSE_REAL].op = (void *) mdr_Lt;

  lt_method[MATRIX_DENSE_COMPLEX][MATRIX_DENSE_COMPLEX].type =
    MATRIX_DENSE_REAL;
  lt_method[MATRIX_DENSE_COMPLEX][MATRIX_DENSE_COMPLEX].op = (void *) mdc_Lt;

  lt_method[MATRIX_DENSE_REAL][MATRIX_DENSE_COMPLEX].type = MATRIX_DENSE_REAL;
  lt_method[MATRIX_DENSE_REAL][MATRIX_DENSE_COMPLEX].op = (void *) mdr_mdc_Lt;

  lt_method[MATRIX_DENSE_COMPLEX][MATRIX_DENSE_REAL].type = MATRIX_DENSE_REAL;
  lt_method[MATRIX_DENSE_COMPLEX][MATRIX_DENSE_REAL].op = (void *) mdc_mdr_Lt;

  lt_method[MATRIX_DENSE_STRING][MATRIX_DENSE_STRING].type = MATRIX_DENSE_REAL;
  lt_method[MATRIX_DENSE_STRING][MATRIX_DENSE_STRING].op = (void *) mds_Lt;

  /*
   * Compare (<=) two objects.
   */

  le_method[MATRIX_DENSE_REAL][MATRIX_DENSE_REAL].type = MATRIX_DENSE_REAL;
  le_method[MATRIX_DENSE_REAL][MATRIX_DENSE_REAL].op = (void *) mdr_Le;

  le_method[MATRIX_DENSE_COMPLEX][MATRIX_DENSE_COMPLEX].type =
    MATRIX_DENSE_REAL;
  le_method[MATRIX_DENSE_COMPLEX][MATRIX_DENSE_COMPLEX].op = (void *) mdc_Le;

  le_method[MATRIX_DENSE_REAL][MATRIX_DENSE_COMPLEX].type = MATRIX_DENSE_REAL;
  le_method[MATRIX_DENSE_REAL][MATRIX_DENSE_COMPLEX].op = (void *) mdr_mdc_Le;

  le_method[MATRIX_DENSE_COMPLEX][MATRIX_DENSE_REAL].type = MATRIX_DENSE_REAL;
  le_method[MATRIX_DENSE_COMPLEX][MATRIX_DENSE_REAL].op = (void *) mdc_mdr_Le;

  le_method[MATRIX_DENSE_STRING][MATRIX_DENSE_STRING].type = MATRIX_DENSE_REAL;
  le_method[MATRIX_DENSE_STRING][MATRIX_DENSE_STRING].op = (void *) mds_Le;

  /*
   * Compare (>) two objects.
   */

  gt_method[MATRIX_DENSE_REAL][MATRIX_DENSE_REAL].type = MATRIX_DENSE_REAL;
  gt_method[MATRIX_DENSE_REAL][MATRIX_DENSE_REAL].op = (void *) mdr_Gt;

  gt_method[MATRIX_DENSE_REAL][MATRIX_DENSE_COMPLEX].type = MATRIX_DENSE_REAL;
  gt_method[MATRIX_DENSE_REAL][MATRIX_DENSE_COMPLEX].op = (void *) mdr_mdc_Gt;

  gt_method[MATRIX_DENSE_COMPLEX][MATRIX_DENSE_COMPLEX].type =
    MATRIX_DENSE_REAL;
  gt_method[MATRIX_DENSE_COMPLEX][MATRIX_DENSE_COMPLEX].op = (void *) mdc_Gt;

  gt_method[MATRIX_DENSE_COMPLEX][MATRIX_DENSE_REAL].type = MATRIX_DENSE_REAL;
  gt_method[MATRIX_DENSE_COMPLEX][MATRIX_DENSE_REAL].op = (void *) mdc_mdr_Gt;

  gt_method[MATRIX_DENSE_STRING][MATRIX_DENSE_STRING].type = MATRIX_DENSE_REAL;
  gt_method[MATRIX_DENSE_STRING][MATRIX_DENSE_STRING].op = (void *) mds_Gt;

  /*
   * Compare (>=) two objects.
   */

  ge_method[MATRIX_DENSE_REAL][MATRIX_DENSE_REAL].type = MATRIX_DENSE_REAL;
  ge_method[MATRIX_DENSE_REAL][MATRIX_DENSE_REAL].op = (void *) mdr_Ge;

  ge_method[MATRIX_DENSE_COMPLEX][MATRIX_DENSE_COMPLEX].type =
    MATRIX_DENSE_REAL;
  ge_method[MATRIX_DENSE_COMPLEX][MATRIX_DENSE_COMPLEX].op = (void *) mdc_Ge;

  ge_method[MATRIX_DENSE_REAL][MATRIX_DENSE_COMPLEX].type = MATRIX_DENSE_REAL;
  ge_method[MATRIX_DENSE_REAL][MATRIX_DENSE_COMPLEX].op = (void *) mdr_mdc_Ge;

  ge_method[MATRIX_DENSE_COMPLEX][MATRIX_DENSE_REAL].type = MATRIX_DENSE_REAL;
  ge_method[MATRIX_DENSE_COMPLEX][MATRIX_DENSE_REAL].op = (void *) mdc_mdr_Ge;

  ge_method[MATRIX_DENSE_STRING][MATRIX_DENSE_STRING].type = MATRIX_DENSE_REAL;
  ge_method[MATRIX_DENSE_STRING][MATRIX_DENSE_STRING].op = (void *) mds_Ge;

  /*
   * Compare (&&) two objects.
   */

  and_method[MATRIX_DENSE_REAL][MATRIX_DENSE_REAL].type = MATRIX_DENSE_REAL;
  and_method[MATRIX_DENSE_REAL][MATRIX_DENSE_REAL].op = (void *) mdr_And;

  and_method[MATRIX_DENSE_COMPLEX][MATRIX_DENSE_COMPLEX].type =
    MATRIX_DENSE_REAL;
  and_method[MATRIX_DENSE_COMPLEX][MATRIX_DENSE_COMPLEX].op = (void *) mdc_And;

  and_method[MATRIX_DENSE_REAL][MATRIX_DENSE_COMPLEX].type = MATRIX_DENSE_REAL;
  and_method[MATRIX_DENSE_REAL][MATRIX_DENSE_COMPLEX].op = (void *) mdr_mdc_And;

  and_method[MATRIX_DENSE_COMPLEX][MATRIX_DENSE_REAL].type = MATRIX_DENSE_REAL;
  and_method[MATRIX_DENSE_COMPLEX][MATRIX_DENSE_REAL].op = (void *) mdc_mdr_And;

  and_method[MATRIX_DENSE_STRING][MATRIX_DENSE_STRING].type = MATRIX_DENSE_REAL;
  and_method[MATRIX_DENSE_STRING][MATRIX_DENSE_STRING].op = (void *) mds_And;

  /*
   * Compare (||) two objects.
   */

  or_method[MATRIX_DENSE_REAL][MATRIX_DENSE_REAL].type = MATRIX_DENSE_REAL;
  or_method[MATRIX_DENSE_REAL][MATRIX_DENSE_REAL].op = (void *) mdr_Or;

  or_method[MATRIX_DENSE_COMPLEX][MATRIX_DENSE_COMPLEX].type =
    MATRIX_DENSE_REAL;
  or_method[MATRIX_DENSE_COMPLEX][MATRIX_DENSE_COMPLEX].op = (void *) mdc_Or;

  or_method[MATRIX_DENSE_STRING][MATRIX_DENSE_STRING].type = MATRIX_DENSE_REAL;
  or_method[MATRIX_DENSE_STRING][MATRIX_DENSE_STRING].op = (void *) mds_Or;

  /*
   * ! and object
   */

  not_method[MATRIX_DENSE_REAL].type = MATRIX_DENSE_REAL;
  not_method[MATRIX_DENSE_REAL].op = (void *) mdr_Not;

  not_method[MATRIX_DENSE_COMPLEX].type = MATRIX_DENSE_REAL;
  not_method[MATRIX_DENSE_COMPLEX].op = (void *) mdc_Not;

  not_method[MATRIX_DENSE_STRING].type = MATRIX_DENSE_REAL;
  not_method[MATRIX_DENSE_STRING].op = (void *) mds_Not;

  /*
   * Transpose an object.
   */

  transpose_method[MATRIX_DENSE_REAL].type = MATRIX_DENSE_REAL;
  transpose_method[MATRIX_DENSE_REAL].op = (void *) mdr_Transpose;

  transpose_method[MATRIX_DENSE_COMPLEX].type = MATRIX_DENSE_COMPLEX;
  transpose_method[MATRIX_DENSE_COMPLEX].op = (void *) mdc_Transpose;

  transpose_method[MATRIX_DENSE_STRING].type = MATRIX_DENSE_STRING;
  transpose_method[MATRIX_DENSE_STRING].op = (void *) mds_Transpose;

  transpose_method[MATRIX_SPARSE_REAL].type = MATRIX_SPARSE_REAL;
  transpose_method[MATRIX_SPARSE_REAL].op = (void *) msr_Transpose;

  transpose_method[MATRIX_SPARSE_COMPLEX].type = MATRIX_SPARSE_COMPLEX;
  transpose_method[MATRIX_SPARSE_COMPLEX].op = (void *) msc_Transpose;

  nc_transpose_method[MATRIX_DENSE_REAL].type = MATRIX_DENSE_REAL;
  nc_transpose_method[MATRIX_DENSE_REAL].op = (void *) mdr_Transpose;

  nc_transpose_method[MATRIX_DENSE_COMPLEX].type = MATRIX_DENSE_COMPLEX;
  nc_transpose_method[MATRIX_DENSE_COMPLEX].op = (void *) mdc_NcTranspose;

  nc_transpose_method[MATRIX_DENSE_STRING].type = MATRIX_DENSE_STRING;
  nc_transpose_method[MATRIX_DENSE_STRING].op = (void *) mds_Transpose;

  nc_transpose_method[MATRIX_SPARSE_REAL].type = MATRIX_SPARSE_REAL;
  nc_transpose_method[MATRIX_SPARSE_REAL].op = (void *) msr_Transpose;

  nc_transpose_method[MATRIX_SPARSE_COMPLEX].type = MATRIX_SPARSE_COMPLEX;
  nc_transpose_method[MATRIX_SPARSE_COMPLEX].op = (void *) msc_NcTranspose;

  /*
   * Reshape a matrix into a column vector.
   */

  reshape_col_method[MATRIX_DENSE_REAL].type = MATRIX_DENSE_REAL;
  reshape_col_method[MATRIX_DENSE_REAL].op = (void *) mdr_ReshapeCol;

  reshape_col_method[MATRIX_DENSE_COMPLEX].type = MATRIX_DENSE_COMPLEX;
  reshape_col_method[MATRIX_DENSE_COMPLEX].op = (void *) mdc_ReshapeCol;

  reshape_col_method[MATRIX_DENSE_STRING].type = MATRIX_DENSE_STRING;
  reshape_col_method[MATRIX_DENSE_STRING].op = (void *) mds_ReshapeCol;

  reshape_col_method[MATRIX_SPARSE_REAL].type = MATRIX_SPARSE_REAL;
  reshape_col_method[MATRIX_SPARSE_REAL].op = (void *) msr_ReshapeCol;

  reshape_col_method[MATRIX_SPARSE_COMPLEX].type = MATRIX_SPARSE_COMPLEX;
  reshape_col_method[MATRIX_SPARSE_COMPLEX].op = (void *) msc_ReshapeCol;

  /*
   * Increment (++) an object.
   */

  increment_method[MATRIX_DENSE_REAL].type = MATRIX_DENSE_REAL;
  increment_method[MATRIX_DENSE_REAL].op = (void *) mdr_Increment;

  increment_method[MATRIX_DENSE_COMPLEX].type = MATRIX_DENSE_COMPLEX;
  increment_method[MATRIX_DENSE_COMPLEX].op = (void *) mdc_Increment;

  /*
   * Decrement (--) an object.
   */

  decrement_method[MATRIX_DENSE_REAL].type = MATRIX_DENSE_REAL;
  decrement_method[MATRIX_DENSE_REAL].op = (void *) mdr_Decrement;

  decrement_method[MATRIX_DENSE_COMPLEX].type = MATRIX_DENSE_COMPLEX;
  decrement_method[MATRIX_DENSE_COMPLEX].op = (void *) mdc_Decrement;

  /*
   * Get an object's sublist (if it exists).
   * These functions _must_ return a BTREE!
   */

  attribute_method[BTREE].type = -1;
  attribute_method[BTREE].op = (void *) btree_Sublist;

  attribute_method[U_FUNCTION].type = -1;
  attribute_method[U_FUNCTION].op = (void *) function_Sublist;

  attribute_method[MATRIX_DENSE_REAL].type = -1;
  attribute_method[MATRIX_DENSE_REAL].op = (void *) mdr_Sublist;

  attribute_method[MATRIX_SPARSE_REAL].type = -1;
  attribute_method[MATRIX_SPARSE_REAL].op = (void *) msr_Sublist;

  attribute_method[MATRIX_DENSE_COMPLEX].type = -1;
  attribute_method[MATRIX_DENSE_COMPLEX].op = (void *) mdc_Sublist;

  attribute_method[MATRIX_SPARSE_COMPLEX].type = -1;
  attribute_method[MATRIX_SPARSE_COMPLEX].op = (void *) msc_Sublist;

  attribute_method[MATRIX_DENSE_STRING].type = -1;
  attribute_method[MATRIX_DENSE_STRING].op = (void *) mds_Sublist;

}

/*
 * This array allows rlab to tell the user,
 * in an understandable way what an entity
 * type is.
 */

static char *ent_desc[] =
{
  "Undefined",
  "Double",
  "Imaginary_Double",
  "List",
  "Builtin_Function",
  "User_Function",
  "",
  "",
  "Matrix_Sparse_Complex",
  "Matrix_Sparse_Real",
  "Matrix_Dense_Real",
  "Matrix_Dense_Complex",
  "Matrix_Dense_String",
  "Matrix_Dense_Entity"
};

char * etd (Ent * e)
{
  if (e)
  {
    if (ent_type(e) < NCL)
    {
      return (ent_desc[ent_type (e)]);
    }
  }

  return (ent_desc[0]);
}

/* **************************************************************
 * Copy an arbitrary object.
 * ************************************************************** */
void ent_data_class_copy(Ent *old, Ent *new_content)
{
  int type=UNDEF, rtype;
  void *rval=0;
  void *(*vfptr) ();

  type = ent_type (new_content);
  rtype = copy_method[type].type;
  vfptr = (VFPTR) copy_method[type].op;

  if (vfptr == 0)
    return;

  rval = (*vfptr) (ent_data (new_content));
  ent_data (old) = rval;
  ent_type (old) = rtype;
}


Ent *
class_copy (Ent * e)
{
  Ent *new;
  int type=UNDEF, rtype;
  void *rval=0;
  void *(*vfptr) ();

  if (e)
    type = ent_type (e);

  rtype = copy_method[type].type;
  vfptr = (VFPTR) copy_method[type].op;

  if (vfptr == 0)
  {
    fprintf (stderr, "Entity type: %s\n", etd (e));
    rerror ("copy operation not supported");
  }

  rval = (*vfptr) (ent_data (e));

  new = ent_Create ();
  ent_data (new) = rval;
  ent_type (new) = rtype;

  return (new);
}

/* **************************************************************
 * Destroy an arbitrary object.
 * ************************************************************** */
void class_destroy (Ent * e)
{
  int type=UNDEF;
  void *(*vfptr) ();

  if (e)
    type = ent_type (e);
  vfptr = (VFPTR) destroy_method[type].op;

  if (vfptr == 0)
  {
    fprintf (stderr, "Entity type: %s\n", etd (e));
    rerror ("destroy operation not supported");
  }

  (*vfptr) (ent_data (e));
}

/* **************************************************************
 * Print an object belonging to a defined class.
 * ************************************************************** */
#undef  THIS_SOLVER
#define THIS_SOLVER "class_print"
void class_print (Ent * e)
{
  FILE *rlab_stderr = (!RLAB_STDERR_DS) ? stderr : RLAB_STDERR_DS;
  int type=UNDEF;
  void *(*vfptr) ();
  FILE *diary_stream;

  if (e)
    type = ent_type (e);

  if (type != MATRIX_DENSE_ENTITY)
  {
    vfptr = (VFPTR) print_method[type].op;
    if (vfptr == 0)
    {
      fprintf (stderr, THIS_SOLVER ": Entity type: %s\n", etd (e));
      rerror (THIS_SOLVER ": " RLAB_ERROR_PRINT_FAILED "\n");
    }
    signal (SIGINT, intcatch);
    (*vfptr) (ent_data (e), rlab_stderr);
    if (get_write_diary ())
    {
      diary_stream = get_diary_file_ptr ();
      (*vfptr) (ent_data (e), diary_stream);
    }
    signal (SIGINT, intcatch_wait);
  }
  else
  {
    // for MDE's if it is of size 1-by-1 look what is it and then print
    MDE * edummy = ent_data(e);
    fprintf(rlab_stderr, "%s %i-by-%i:\n", etd(e), MNR(edummy), MNC(edummy));
    int i, j;
    for (i=0; i<MNR(edummy); i++)
    {
      for (j=0; j<MNC(edummy); j++)
      {
        Ent *e1 = Mde0(edummy,i,j);
        fprintf(rlab_stderr, "[%i;%i] -> %s\n", i+1, j+1, etd(e1));
/*        if (e1)
          fprintf(stdout, "[%i;%i] -> %s, refc = %i\n", i+1, j+1, etd(e1), ent_Ref(e1));
        else
          fprintf(stdout, "[%i;%i] -> %s, refc = 0\n", i+1, j+1, etd(e1));*/
        if (ent_type(e1)!=UNDEF)
        {
          class_print(e1);
        }
      }
    }
  }

}

/* **************************************************************
 * Negate an object.
 * ************************************************************** */

Ent *
class_negate (Ent * e)
{
  Ent *new;
  int type=UNDEF, rtype;
  void *rval=0;
  void *(*vfptr) ();

  if (e)
    type = ent_type (e);

  rtype = negate_method[type].type;
  vfptr = (VFPTR) negate_method[type].op;

  if (vfptr == 0)
  {
    fprintf (stderr, "Entity type: %s\n", etd (e));
    rerror ("negation operation not supported");
  }

  rval = (*vfptr) (ent_data (e));

  new = ent_Create ();
  ent_data (new) = rval;
  ent_type (new) = rtype;

  return (new);
}

#undef THIS_SOLVER
#undef METHOD
#undef METHOD_WEIGHT
#define THIS_SOLVER "add"
#define METHOD add_method
#define METHOD_WEIGHT add_method_weight
Ent * class_add (Ent * e1, Ent * e2)
{
  Ent *new=0, *ex1=0, *ex2=0;
  MDR *w=0, *w1=0, *w2=0;
  Btree *bw=0;
  ListNode *node1=0, *node2=0;
  int type1=UNDEF, type2=UNDEF, rtype;
  void *rval=0;
  void *(*vfptr) ();
  void *(*vfptr_wgt) ();

  if (e1)
    type1 = ent_type (e1);
  if (e2)
    type2 = ent_type (e2);

  // process first argument
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
      if (w1)
        if (!(MNR(w1)*MNC(w1)))
          w1=0;
    }
  }
  else if (   type1 == MATRIX_DENSE_REAL
           || type1 == MATRIX_DENSE_COMPLEX
           || type1 == MATRIX_SPARSE_COMPLEX
           || type1 == MATRIX_SPARSE_REAL
           || type1 == MATRIX_DENSE_STRING  )
  {
    ex1 = e1;
  }

  // process second argument
  if (type2 == BTREE)
  {
    node1 = btree_FindNode (ent_data (e2), RLAB_NAME_STAT_VALUE);
    if (!node1)
      rerror (THIS_SOLVER ": invalid second argument <<" RLAB_NAME_STAT_VALUE ";" RLAB_NAME_STAT_WEIGHT ">>");
    ex2   = var_ent (node1);
    type2 = ent_type(ex2);

    // fetch its weights
    node2 = btree_FindNode (ent_data (e2), RLAB_NAME_STAT_WEIGHT);
    if (node2)
    {
      if (ent_type(var_ent (node2)) != MATRIX_DENSE_REAL)
        rerror (THIS_SOLVER ": Second argument list entry " RLAB_NAME_STAT_WEIGHT " must be real matrix !");
      w2 = class_matrix_real (var_ent (node2));
      if (w2)
        if (!(MNR(w2)*MNC(w2)))
          w2=0;
    }
  }
  else if (  type2 == MATRIX_DENSE_REAL
          || type2 == MATRIX_DENSE_COMPLEX
          || type2 == MATRIX_SPARSE_COMPLEX
          || type2 == MATRIX_SPARSE_REAL
          || type2 == MATRIX_DENSE_STRING  )
  {
    ex2 = e2;
  }

  rtype = METHOD[type1][type2].type;
  vfptr = (VFPTR) METHOD[type1][type2].op;
  if (vfptr)
  {
    rval = (*vfptr) (ent_data (ex1),ent_data (ex2));
  }
  else
  {
    fprintf (stderr, "Entity types: %s and %s\n", etd (ex1), etd (ex2));
    rerror (THIS_SOLVER ": operation not supported");
  }
  new = ent_Create ();

  if (w1 || w2)
  {
    vfptr_wgt = (VFPTR) METHOD_WEIGHT[type1][type2].op;
    if (vfptr_wgt)
      w = (*vfptr_wgt) (ent_data (ex1), w1, ent_data (ex2), w2);

    if (w)
    {
      // calculation of weights was successful
      bw = btree_Create();

      // copy result to list entry 'val'
      install  (bw, RLAB_NAME_STAT_VALUE, ent_Assign_Rlab_MDR(rval));

      // copy result to list entry 'val'
      install  (bw, RLAB_NAME_STAT_WEIGHT, ent_Assign_Rlab_MDR(w));

      ent_data (new) = bw;
      ent_type (new) = BTREE;
    }
    else
    {
      // calculations of weights have failed
      ent_data (new) = rval;
      ent_type (new) = rtype;
    }
  }
  else
  {
    ent_data (new) = rval;
    ent_type (new) = rtype;
  }

  return (new);
}

#undef  THIS_SOLVER
#undef  METHOD
#define THIS_SOLVER "addto"
#define METHOD addto_method
Ent * class_addto (Ent * e1, Ent * e2)
{
  Ent *new=0, *ex1=0, *ex2=0;
  int type1=UNDEF, type2=UNDEF, rtype;
  void *rval=0;
  void *(*vfptr) ();

  if (e1)
    type1 = ent_type (e1);
  if (e2)
    type2 = ent_type (e2);

  // process first argument
  if(type1 == MATRIX_DENSE_REAL || type1 == MATRIX_DENSE_COMPLEX || type1 == MATRIX_DENSE_STRING  )
  {
    ex1 = e1;
  }

  // process second argument
  if (type2 == MATRIX_DENSE_REAL || type2 == MATRIX_DENSE_COMPLEX || type2 == MATRIX_DENSE_STRING  )
  {
    ex2 = e2;
  }

  rtype = METHOD[type1][type2].type;
  vfptr = (VFPTR) METHOD[type1][type2].op;
  if (vfptr)
  {
    rval = (*vfptr) (ent_data (ex1),ent_data (ex2));
  }
  else
  {
    fprintf (stderr, "Entity types: %s and %s\n", etd (ex1), etd (ex2));
    rerror (THIS_SOLVER ": operation not supported");
  }

  new = ent_Create ();
  ent_data (new) = rval;
  ent_type (new) = rtype;
  return (new);
}

#undef  THIS_SOLVER
#undef  METHOD
#define THIS_SOLVER "subtractfrom"
#define METHOD subtractfrom_method
Ent * class_subtractfrom (Ent * e1, Ent * e2)
{
  Ent *new=0, *ex1=0, *ex2=0;
  int type1=UNDEF, type2=UNDEF, rtype;
  void *rval=0;
  void *(*vfptr) ();

  if (e1)
    type1 = ent_type (e1);
  if (e2)
    type2 = ent_type (e2);

  // process first argument
  if(type1 == MATRIX_DENSE_REAL || type1 == MATRIX_DENSE_COMPLEX || type1 == MATRIX_DENSE_STRING  )
  {
    ex1 = e1;
  }

  // process second argument
  if (type2 == MATRIX_DENSE_REAL || type2 == MATRIX_DENSE_COMPLEX || type2 == MATRIX_DENSE_STRING  )
  {
    ex2 = e2;
  }

  rtype = METHOD[type1][type2].type;
  vfptr = (VFPTR) METHOD[type1][type2].op;
  if (vfptr)
  {
    rval = (*vfptr) (ent_data (ex1),ent_data (ex2));
  }
  else
  {
    fprintf (stderr, "Entity types: %s and %s\n", etd (ex1), etd (ex2));
    rerror (THIS_SOLVER ": operation not supported");
  }

  new = ent_Create ();
  ent_data (new) = rval;
  ent_type (new) = rtype;
  return (new);
}


#undef  THIS_SOLVER
#undef  METHOD
#define THIS_SOLVER "elmultiplyby"
#define METHOD elmultiplyby_method
Ent * class_elmultiplyby (Ent * e1, Ent * e2)
{
  Ent *new=0, *ex1=0, *ex2=0;
  int type1=UNDEF, type2=UNDEF, rtype;
  void *rval=0;
  void *(*vfptr) ();

  if (e1)
    type1 = ent_type (e1);
  if (e2)
    type2 = ent_type (e2);

  // process first argument
  if(type1 == MATRIX_DENSE_REAL || type1 == MATRIX_DENSE_COMPLEX || type1 == MATRIX_DENSE_STRING  )
  {
    ex1 = e1;
  }

  // process second argument
  if (type2 == MATRIX_DENSE_REAL || type2 == MATRIX_DENSE_COMPLEX || type2 == MATRIX_DENSE_STRING  )
  {
    ex2 = e2;
  }

  rtype = METHOD[type1][type2].type;
  vfptr = (VFPTR) METHOD[type1][type2].op;
  if (vfptr)
  {
    rval = (*vfptr) (ent_data (ex1),ent_data (ex2));
  }
  else
  {
    fprintf (stderr, "Entity types: %s and %s\n", etd (ex1), etd (ex2));
    rerror (THIS_SOLVER ": operation not supported");
  }

  new = ent_Create ();
  ent_data (new) = rval;
  ent_type (new) = rtype;
  return (new);
}


#undef  THIS_SOLVER
#undef  METHOD
#define THIS_SOLVER "elrdivideby"
#define METHOD elrdivideby_method
Ent * class_elrdivideby (Ent * e1, Ent * e2)
{
  Ent *new=0, *ex1=0, *ex2=0;
  int type1=UNDEF, type2=UNDEF, rtype;
  void *rval=0;
  void *(*vfptr) ();

  if (e1)
    type1 = ent_type (e1);
  if (e2)
    type2 = ent_type (e2);

  // process first argument
  if(type1 == MATRIX_DENSE_REAL || type1 == MATRIX_DENSE_COMPLEX || type1 == MATRIX_DENSE_STRING  )
  {
    ex1 = e1;
  }

  // process second argument
  if (type2 == MATRIX_DENSE_REAL || type2 == MATRIX_DENSE_COMPLEX || type2 == MATRIX_DENSE_STRING  )
  {
    ex2 = e2;
  }

  rtype = METHOD[type1][type2].type;
  vfptr = (VFPTR) METHOD[type1][type2].op;
  if (vfptr)
  {
    rval = (*vfptr) (ent_data (ex1),ent_data (ex2));
  }
  else
  {
    fprintf (stderr, "Entity types: %s and %s\n", etd (ex1), etd (ex2));
    rerror (THIS_SOLVER ": operation not supported");
  }

  new = ent_Create ();
  ent_data (new) = rval;
  ent_type (new) = rtype;
  return (new);
}

// **************************************************************
// Subtract two objects.
// **************************************************************
#undef THIS_SOLVER
#undef METHOD
#undef METHOD_WEIGHT
#define THIS_SOLVER "operation: subtract"
#define METHOD subtract_method
#define METHOD_WEIGHT subtract_method_weight
Ent *
class_subtract (Ent * e1, Ent * e2)
{
  Ent *new=0, *ex1=0, *ex2=0;
  MDR *w=0, *w1=0, *w2=0;
  Btree *bw=0;
  ListNode *node1=0, *node2=0;
  int type1=UNDEF, type2=UNDEF, rtype;
  void *rval=0;
  void *(*vfptr) ();
  void *(*vfptr_wgt) ();

  if (e1)
    type1 = ent_type (e1);
  if (e2)
    type2 = ent_type (e2);

  // process first argument
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
      if (w1)
        if (!(MNR(w1)*MNC(w1)))
          w1=0;
    }
  }
  else if (   type1 == MATRIX_DENSE_REAL
           || type1 == MATRIX_DENSE_COMPLEX
           || type1 == MATRIX_SPARSE_COMPLEX
           || type1 == MATRIX_SPARSE_REAL
           || type1 == MATRIX_DENSE_STRING  )
  {
    ex1 = e1;
  }

  // process second argument
  if (type2 == BTREE)
  {
    node1 = btree_FindNode (ent_data (e2), RLAB_NAME_STAT_VALUE);
    if (!node1)
      rerror (THIS_SOLVER ": invalid second argument <<" RLAB_NAME_STAT_VALUE ";" RLAB_NAME_STAT_WEIGHT ">>");
    ex2   = var_ent (node1);
    type2 = ent_type(ex2);

    // fetch its weights
    node2 = btree_FindNode (ent_data (e2), RLAB_NAME_STAT_WEIGHT);
    if (node2)
    {
      if (ent_type(var_ent (node2)) != MATRIX_DENSE_REAL)
        rerror (THIS_SOLVER ": Second argument list entry " RLAB_NAME_STAT_WEIGHT " must be real matrix !");
      w2 = class_matrix_real (var_ent (node2));
      if (w2)
        if (!(MNR(w2)*MNC(w2)))
          w2=0;
    }
  }
  else if (  type2 == MATRIX_DENSE_REAL
          || type2 == MATRIX_DENSE_COMPLEX
          || type2 == MATRIX_SPARSE_COMPLEX
          || type2 == MATRIX_SPARSE_REAL
          || type2 == MATRIX_DENSE_STRING  )
  {
    ex2 = e2;
  }

  rtype = METHOD[type1][type2].type;
  vfptr = (VFPTR) METHOD[type1][type2].op;

  if (vfptr)
  {
    rval = (*vfptr) (ent_data (ex1),ent_data (ex2));
  }
  else
  {
    fprintf (stderr, "Entity types: %s and %s\n", etd (ex1), etd (ex2));
    rerror (THIS_SOLVER ": operation not supported");
  }



  new = ent_Create ();

  if (w1 || w2)
  {
    vfptr_wgt = (VFPTR) METHOD_WEIGHT[type1][type2].op;
    if (vfptr_wgt)
      w = (*vfptr_wgt) (ent_data (ex1), w1, ent_data (ex2), w2);

    if (w)
    {
      // calculation of weights was successful
      bw = btree_Create();

      // copy result to list entry 'val'
      install  (bw, RLAB_NAME_STAT_VALUE, ent_Assign_Rlab_MDR(rval));

      // copy result to list entry 'val'
      install  (bw, RLAB_NAME_STAT_WEIGHT, ent_Assign_Rlab_MDR(w));

      ent_data (new) = bw;
      ent_type (new) = BTREE;
    }
    else
    {
      // calculations of weights have failed
      ent_data (new) = rval;
      ent_type (new) = rtype;
    }
  }
  else
  {
    ent_data (new) = rval;
    ent_type (new) = rtype;
  }

  return (new);
}

#undef THIS_SOLVER
#undef METHOD
#undef METHOD_WEIGHT
#define THIS_SOLVER "operation: multiply"
#define METHOD multiply_method
#define METHOD_WEIGHT multiply_method_weight
Ent *
class_multiply (Ent * e1, Ent * e2)
{
  Ent *new=0, *ex1=0, *ex2=0;
  MDR *w=0, *w1=0, *w2=0;
  Btree *bw=0;
  ListNode *node1=0, *node2=0;
  int type1=UNDEF, type2=UNDEF, rtype;
  void *rval=0;
  void *(*vfptr) ();
  void *(*vfptr_wgt) ();

  if (e1)
    type1 = ent_type (e1);
  if (e2)
    type2 = ent_type (e2);

  // process first argument
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
      if (w1)
        if (!(MNR(w1)*MNC(w1)))
          w1=0;
    }
  }
  else if (   type1 == MATRIX_DENSE_REAL
           || type1 == MATRIX_DENSE_COMPLEX
           || type1 == MATRIX_SPARSE_COMPLEX
           || type1 == MATRIX_SPARSE_REAL)
  {
    ex1 = e1;
  }

  // process second argument
  if (type2 == BTREE)
  {
    node1 = btree_FindNode (ent_data (e2), RLAB_NAME_STAT_VALUE);
    if (!node1)
      rerror (THIS_SOLVER ": invalid second argument <<" RLAB_NAME_STAT_VALUE ";" RLAB_NAME_STAT_WEIGHT ">>");
    ex2   = var_ent (node1);
    type2 = ent_type(ex2);

    // fetch its weights
    node2 = btree_FindNode (ent_data (e2), RLAB_NAME_STAT_WEIGHT);
    if (node2)
    {
      if (ent_type(var_ent (node2)) != MATRIX_DENSE_REAL)
        rerror (THIS_SOLVER ": Second argument list entry " RLAB_NAME_STAT_WEIGHT " must be real matrix !");
      w2 = class_matrix_real (var_ent (node2));
      if (w2)
        if (!(MNR(w2)*MNC(w2)))
          w2=0;
    }
  }
  else if (  type2 == MATRIX_DENSE_REAL
          || type2 == MATRIX_DENSE_COMPLEX
          || type2 == MATRIX_SPARSE_COMPLEX
          || type2 == MATRIX_SPARSE_REAL)
  {
    ex2 = e2;
  }

  vfptr = (VFPTR) METHOD[type1][type2].op;
  if (vfptr)
  {
    rval = (*vfptr) (ent_data (ex1),ent_data (ex2), &rtype);
  }
  else
  {
    fprintf (stderr, "%s: Entity types: %s and %s\n", THIS_SOLVER, etd (ex1), etd (ex2));
    fprintf (stderr, "%s: %s\n", THIS_SOLVER, RLAB_ERROR_OPERATION_FAILED);
    rerror (THIS_SOLVER);
  }

  new = ent_Create ();

  if (w1 || w2)
  {
    vfptr_wgt = (VFPTR) METHOD_WEIGHT[type1][type2].op;
    if (vfptr_wgt)
      w = (*vfptr_wgt) (ent_data (ex1), w1, ent_data (ex2), w2);

    if (w)
    {
      // calculation of weights was successful
      bw = btree_Create();

      // copy result to list entry 'val'
      install  (bw, RLAB_NAME_STAT_VALUE, ent_Assign_Rlab_MDR(rval));

      // copy result to list entry 'val'
      install  (bw, RLAB_NAME_STAT_WEIGHT, ent_Assign_Rlab_MDR(w));

      ent_data (new) = bw;
      ent_type (new) = BTREE;
    }
    else
    {
      // calculations of weights have failed
      ent_data (new) = rval;
      ent_type (new) = rtype;
    }
  }
  else
  {
    ent_data (new) = rval;
    ent_type (new) = rtype;
  }

  return (new);
}

#undef THIS_SOLVER
#undef METHOD
#undef METHOD_WEIGHT
#define THIS_SOLVER "operation: element multiply"
#define METHOD el_multiply_method
#define METHOD_WEIGHT el_multiply_method_weight
Ent *
class_el_multiply (Ent * e1, Ent * e2)
{
  Ent *new=0, *ex1=0, *ex2=0;
  MDR *w=0, *w1=0, *w2=0;
  Btree *bw=0;
  ListNode *node1=0, *node2=0;
  int type1=UNDEF, type2=UNDEF, rtype;
  void *rval=0;
  void *(*vfptr) ();
  void *(*vfptr_wgt) ();

  if (e1)
    type1 = ent_type (e1);
  if (e2)
    type2 = ent_type (e2);

  // process first argument
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
      if (w1)
        if (!(MNR(w1)*MNC(w1)))
          w1=0;
    }
  }
  else if (   type1 == MATRIX_DENSE_REAL
           || type1 == MATRIX_DENSE_COMPLEX
           || type1 == MATRIX_SPARSE_COMPLEX
           || type1 == MATRIX_SPARSE_REAL)
  {
    ex1 = e1;
  }

  // process second argument
  if (type2 == BTREE)
  {
    node1 = btree_FindNode (ent_data (e2), RLAB_NAME_STAT_VALUE);
    if (!node1)
      rerror (THIS_SOLVER ": invalid second argument <<" RLAB_NAME_STAT_VALUE ";" RLAB_NAME_STAT_WEIGHT ">>");
    ex2   = var_ent (node1);
    type2 = ent_type(ex2);

    // fetch its weights
    node2 = btree_FindNode (ent_data (e2), RLAB_NAME_STAT_WEIGHT);
    if (node2)
    {
      if (ent_type(var_ent (node2)) != MATRIX_DENSE_REAL)
        rerror (THIS_SOLVER ": Second argument list entry " RLAB_NAME_STAT_WEIGHT " must be real matrix !");
      w2 = class_matrix_real (var_ent (node2));
      if (w2)
        if (!(MNR(w2)*MNC(w2)))
          w2=0;
    }
  }
  else if (  type2 == MATRIX_DENSE_REAL
          || type2 == MATRIX_DENSE_COMPLEX
          || type2 == MATRIX_SPARSE_COMPLEX
          || type2 == MATRIX_SPARSE_REAL)
  {
    ex2 = e2;
  }

  rtype = METHOD[type1][type2].type;
  vfptr = (VFPTR) METHOD[type1][type2].op;

  if (vfptr)
  {
    rval = (*vfptr) (ent_data (ex1),ent_data (ex2));
  }
  else
  {
    fprintf (stderr, "Entity types: %s and %s\n", etd (ex1), etd (ex2));
    rerror (THIS_SOLVER ": Operation not supported");
  }


  new = ent_Create ();

  if (w1 || w2)
  {
    vfptr_wgt = (VFPTR) METHOD_WEIGHT[type1][type2].op;
    if (vfptr_wgt)
      w = (*vfptr_wgt) (ent_data (ex1), w1, ent_data (ex2), w2);

    if (w)
    {
      // calculation of weights was successful
      bw = btree_Create();

      // copy result to list entry 'val'
      install  (bw, RLAB_NAME_STAT_VALUE, ent_Assign_Rlab_MDR(rval));

      // copy result to list entry 'val'
      install  (bw, RLAB_NAME_STAT_WEIGHT, ent_Assign_Rlab_MDR(w));

      ent_data (new) = bw;
      ent_type (new) = BTREE;
    }
    else
    {
      // calculations of weights have failed
      ent_data (new) = rval;
      ent_type (new) = rtype;
    }
  }
  else
  {
    ent_data (new) = rval;
    ent_type (new) = rtype;
  }

  return (new);
}

/* **************************************************************
 * Right-divide two objects.
 * ************************************************************** */

Ent *
class_rdivide (Ent * e1, Ent * e2)
{
  Ent *new;
  int type1=UNDEF, type2=UNDEF, rtype;
  void *rval=0;
  void *(*vfptr) ();

  if (e1)
    type1 = ent_type (e1);
  if (e2)
    type2 = ent_type (e2);

  rtype = rdivide_method[type1][type2].type;
  vfptr = (VFPTR) rdivide_method[type1][type2].op;

  if (vfptr == 0)
  {
    fprintf (stderr, "Entity types: %s and %s\n", etd (e1), etd (e2));
    rerror ("right-divide operation not supported");
  }

  rval = (*vfptr) (ent_data (e1), ent_data (e2));

  new = ent_Create ();
  ent_data (new) = rval;
  ent_type (new) = rtype;

  return (new);
}

// **************************************************************
// *
// * Element right-divide two objects.
// *
// **************************************************************
#undef THIS_SOLVER
#undef METHOD
#undef METHOD_WEIGHT
#define THIS_SOLVER "operation: element right divide"
#define METHOD el_rdivide_method
#define METHOD_WEIGHT el_rdivide_method_weight
Ent *
class_el_rdivide (Ent * e1, Ent * e2)
{
  Ent *new=0, *ex1=0, *ex2=0;
  MDR *w=0, *w1=0, *w2=0;
  Btree *bw=0;
  ListNode *node1=0, *node2=0;
  int type1=UNDEF, type2=UNDEF, rtype;
  void *rval=0;
  void *(*vfptr) ();
  void *(*vfptr_wgt) ();

  if (e1)
    type1 = ent_type (e1);
  if (e2)
    type2 = ent_type (e2);

  // process first argument
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
      if (w1)
        if (!(MNR(w1)*MNC(w1)))
          w1=0;
    }
  }
  else if (   type1 == MATRIX_DENSE_REAL
           || type1 == MATRIX_DENSE_COMPLEX
           || type1 == MATRIX_SPARSE_COMPLEX
           || type1 == MATRIX_SPARSE_REAL)
  {
    ex1 = e1;
  }

  // process second argument
  if (type2 == BTREE)
  {
    node1 = btree_FindNode (ent_data (e2), RLAB_NAME_STAT_VALUE);
    if (!node1)
      rerror (THIS_SOLVER ": invalid second argument <<" RLAB_NAME_STAT_VALUE ";" RLAB_NAME_STAT_WEIGHT ">>");
    ex2   = var_ent (node1);
    type2 = ent_type(ex2);

    // fetch its weights
    node2 = btree_FindNode (ent_data (e2), RLAB_NAME_STAT_WEIGHT);
    if (node2)
    {
      if (ent_type(var_ent (node2)) != MATRIX_DENSE_REAL)
        rerror (THIS_SOLVER ": Second argument list entry " RLAB_NAME_STAT_WEIGHT " must be real matrix !");
      w2 = class_matrix_real (var_ent (node2));
      if (w2)
        if (!(MNR(w2)*MNC(w2)))
          w2=0;
    }
  }
  else if (  type2 == MATRIX_DENSE_REAL
          || type2 == MATRIX_DENSE_COMPLEX
          || type2 == MATRIX_SPARSE_COMPLEX
          || type2 == MATRIX_SPARSE_REAL)
  {
    ex2 = e2;
  }

  vfptr = (VFPTR) METHOD[type1][type2].op;
  if (vfptr)
  {
    rval = (*vfptr) (ent_data (ex1),ent_data (ex2), &rtype);
  }
  else
  {
    fprintf (stderr, "Entity types: %s and %s\n", etd (ex1), etd (ex2));
    rerror (THIS_SOLVER ": Operation not supported");
  }

  new = ent_Create ();

  if (w1 || w2)
  {
    vfptr_wgt = (VFPTR) METHOD_WEIGHT[type1][type2].op;
    if (vfptr_wgt)
      w = (*vfptr_wgt) (ent_data (ex1), w1, ent_data (ex2), w2);

    if (w)
    {
      // calculation of weights was successful
      bw = btree_Create();

      // copy result to list entry 'val'
      install  (bw, RLAB_NAME_STAT_VALUE, ent_Assign_Rlab_MDR(rval));

      // copy result to list entry 'val'
      install  (bw, RLAB_NAME_STAT_WEIGHT, ent_Assign_Rlab_MDR(w));

      ent_data (new) = bw;
      ent_type (new) = BTREE;
    }
    else
    {
      // calculations of weights have failed
      ent_data (new) = rval;
      ent_type (new) = rtype;
    }
  }
  else
  {
    ent_data (new) = rval;
    ent_type (new) = rtype;
  }

  return (new);
}

/* **************************************************************
 * Left-divide two objects.
 * ************************************************************** */

Ent *
class_ldivide (Ent * e1, Ent * e2)
{
  Ent *new;
  int type1=UNDEF, type2=UNDEF, rtype;
  void *rval=0;
  void *(*vfptr) ();

  if (e1)
    type1 = ent_type (e1);
  if (e2)
    type2 = ent_type (e2);

  rtype = ldivide_method[type1][type2].type;
  vfptr = (VFPTR) ldivide_method[type1][type2].op;

  if (vfptr == 0)
  {
    fprintf (stderr, "Entity types: %s and %s\n", etd (e1), etd (e2));
    rerror ("left-divide operation not supported");
  }

  rval = (*vfptr) (ent_data (e1), ent_data (e2));

  new = ent_Create ();
  ent_data (new) = rval;
  ent_type (new) = rtype;

  return (new);
}

// **************************************************************
// * Element left-divide two objects.
// * ************************************************************
#undef THIS_SOLVER
#undef METHOD
#undef METHOD_WEIGHT
#define THIS_SOLVER "operation: element left divide"
#define METHOD el_ldivide_method
#define METHOD_WEIGHT el_rdivide_method_weight
Ent *
class_el_ldivide (Ent * e1, Ent * e2)
{
  Ent *new=0, *ex1=0, *ex2=0;
  MDR *w=0, *w1=0, *w2=0;
  Btree *bw=0;
  ListNode *node1=0, *node2=0;
  int type1=UNDEF, type2=UNDEF, rtype;
  void *rval=0;
  void *(*vfptr) ();
  void *(*vfptr_wgt) ();

  if (e1)
    type1 = ent_type (e1);
  if (e2)
    type2 = ent_type (e2);

  // process first argument
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
      if (w1)
        if (!(MNR(w1)*MNC(w1)))
          w1=0;
    }
  }
  else if (   type1 == MATRIX_DENSE_REAL
           || type1 == MATRIX_DENSE_COMPLEX
           || type1 == MATRIX_SPARSE_COMPLEX
           || type1 == MATRIX_SPARSE_REAL)
  {
    ex1 = e1;
  }

  // process second argument
  if (type2 == BTREE)
  {
    node1 = btree_FindNode (ent_data (e2), RLAB_NAME_STAT_VALUE);
    if (!node1)
      rerror (THIS_SOLVER ": invalid second argument <<" RLAB_NAME_STAT_VALUE ";" RLAB_NAME_STAT_WEIGHT ">>");
    ex2   = var_ent (node1);
    type2 = ent_type(ex2);

    // fetch its weights
    node2 = btree_FindNode (ent_data (e2), RLAB_NAME_STAT_WEIGHT);
    if (node2)
    {
      if (ent_type(var_ent (node2)) != MATRIX_DENSE_REAL)
        rerror (THIS_SOLVER ": Second argument list entry " RLAB_NAME_STAT_WEIGHT " must be real matrix !");
      w2 = class_matrix_real (var_ent (node2));
      if (w2)
        if (!(MNR(w2)*MNC(w2)))
          w2=0;
    }
  }
  else if (  type2 == MATRIX_DENSE_REAL
          || type2 == MATRIX_DENSE_COMPLEX
          || type2 == MATRIX_SPARSE_COMPLEX
          || type2 == MATRIX_SPARSE_REAL)
  {
    ex2 = e2;
  }

  vfptr = (VFPTR) METHOD[type1][type2].op;
  if (vfptr)
  {
    rval = (*vfptr) (ent_data (ex1),ent_data (ex2), &rtype);
  }
  else
  {
    fprintf (stderr, "Entity types: %s and %s\n", etd (ex1), etd (ex2));
    rerror (THIS_SOLVER ": operation not supported");
  }
  new = ent_Create ();

  if (w1 || w2)
  {
    vfptr_wgt = (VFPTR) METHOD_WEIGHT[type1][type2].op;
    if (vfptr_wgt)
      w = (*vfptr_wgt) (ent_data (ex2), w2, ent_data (ex1), w1);

    if (w)
    {
      // calculation of weights was successful
      bw = btree_Create();

      // copy result to list entry 'val'
      install  (bw, RLAB_NAME_STAT_VALUE, ent_Assign_Rlab_MDR(rval));

      // copy result to list entry 'val'
      install  (bw, RLAB_NAME_STAT_WEIGHT, ent_Assign_Rlab_MDR(w));

      ent_data (new) = bw;
      ent_type (new) = BTREE;
    }
    else
    {
      // calculations of weights have failed
      ent_data (new) = rval;
      ent_type (new) = rtype;
    }
  }
  else
  {
    ent_data (new) = rval;
    ent_type (new) = rtype;
  }

  return (new);
}


/* **************************************************************
 * Power operation e1 ^ e2
 * ************************************************************** */

Ent *
class_power (Ent * e1, Ent * e2)
{
  Ent *new;
  int type1=UNDEF, type2=UNDEF, rtype;
  void *rval=0;
  void *(*vfptr) ();

  if (e1)
    type1 = ent_type (e1);
  if (e2)
    type2 = ent_type (e2);

  vfptr = (VFPTR) power_method[type1][type2].op;

  if (vfptr == 0)
  {
    fprintf (stderr, "Entity types: %s and %s\n", etd (e1), etd (e2));
    rerror ("power (^) operation not supported");
  }

  rval = (*vfptr) (ent_data (e1), ent_data (e2), &rtype);

  new = ent_Create ();
  ent_data (new) = rval;
  ent_type (new) = rtype;

  return (new);
}

/* **************************************************************
 * Element Power operation e1 .^ e2
 * ************************************************************** */
#undef THIS_SOLVER
#undef METHOD
#undef METHOD_WEIGHT
#define THIS_SOLVER "operation: element power"
#define METHOD el_power_method
#define METHOD_WEIGHT el_power_method_weight
Ent *
class_el_power (Ent * e1, Ent * e2)
{
  Ent *new=0, *ex1=0, *ex2=0;
  MDR *w=0, *w1=0, *w2=0;
  Btree *bw=0;
  ListNode *node1=0, *node2=0;
  int type1=UNDEF, type2=UNDEF, rtype;
  void *rval=0;
  void *(*vfptr) ();
  void *(*vfptr_wgt) ();

  if (e1)
    type1 = ent_type (e1);
  if (e2)
    type2 = ent_type (e2);

  // process first argument
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
      if (w1)
        if (!(MNR(w1)*MNC(w1)))
          w1=0;
    }
  }
  else if (   type1 == MATRIX_DENSE_REAL
              || type1 == MATRIX_DENSE_COMPLEX
              || type1 == MATRIX_SPARSE_COMPLEX
              || type1 == MATRIX_SPARSE_REAL)
  {
    ex1 = e1;
  }

  // process second argument
  if (type2 == BTREE)
  {
    node1 = btree_FindNode (ent_data (e2), RLAB_NAME_STAT_VALUE);
    if (!node1)
      rerror (THIS_SOLVER ": invalid second argument <<" RLAB_NAME_STAT_VALUE ";" RLAB_NAME_STAT_WEIGHT ">>");
    ex2   = var_ent (node1);
    type2 = ent_type(ex2);

    // fetch its weights
    node2 = btree_FindNode (ent_data (e2), RLAB_NAME_STAT_WEIGHT);
    if (node2)
    {
      if (ent_type(var_ent (node2)) != MATRIX_DENSE_REAL)
        rerror (THIS_SOLVER ": Second argument list entry " RLAB_NAME_STAT_WEIGHT " must be real matrix !");
      w2 = class_matrix_real (var_ent (node2));
      if (w2)
        if (!(MNR(w2)*MNC(w2)))
          w2=0;
    }
  }
  else if (  type2 == MATRIX_DENSE_REAL
             || type2 == MATRIX_DENSE_COMPLEX
             || type2 == MATRIX_SPARSE_COMPLEX
             || type2 == MATRIX_SPARSE_REAL)
  {
    ex2 = e2;
  }

  vfptr = (VFPTR) METHOD[type1][type2].op;

  if (vfptr)
  {
    rval = (*vfptr) (ent_data (ex1),ent_data (ex2), &rtype);
  }
  else
  {
    fprintf (stderr, "Entity types: %s and %s\n",etd (e1), etd (e2));
    rerror (THIS_SOLVER ": Operation not supported");
  }


  new = ent_Create ();

  if (w1 || w2)
  {
    vfptr_wgt = (VFPTR) METHOD_WEIGHT[type1][type2].op;
    if (vfptr_wgt)
      w = (*vfptr_wgt) (ent_data (ex1), w1, ent_data (ex2), w2);

    if (w)
    {
      // calculation of weights was successful
      bw = btree_Create();

      // copy result to list entry 'val'
      install  (bw, RLAB_NAME_STAT_VALUE, ent_Assign_Rlab_MDR(rval));

      // copy result to list entry 'wgt'
      install  (bw, RLAB_NAME_STAT_WEIGHT, ent_Assign_Rlab_MDR(w));

      ent_data (new) = bw;
      ent_type (new) = BTREE;
    }
    else
    {
      // calculations of weights have failed
      ent_data (new) = rval;
      ent_type (new) = rtype;
    }
  }
  else
  {
    ent_data (new) = rval;
    ent_type (new) = rtype;
  }

  return (new);
}

// {
//   Ent *new=0, *ex1=0, *ex2=0;
//   MDR *w=0, *w1=0, *w2=0;
//   Btree *bw=0;
//   ListNode *node1=0, *node2=0;
//
//   int type1=UNDEF, type2=UNDEF, rtype;
//   void *rval=0;
//   void *(*vfptr) ();
//   void *(*vfptr_wgt) ();
//
//   if (e1)
//     type1 = ent_type (e1);
//   if (e2)
//     type2 = ent_type (e2);
//
//
//   vfptr = (VFPTR) el_power_method[type1][type2].op;
//
//   if (vfptr == 0)
//   {
//     fprintf (stderr, "Entity types: %s and %s\n", etd (e1), etd (e2));
//     rerror ("element-power (.^) operation not supported");
//   }
//
//   rval = (*vfptr) (etd(e1), etd(e2), &rtype);
//
//   new = ent_Create ();
//   ent_data (new) = rval;
//   ent_type (new) = rtype;
//   return (new);
// }

/* **************************************************************
 * Return a logical (TRUE or FALSE) scalar value, given an
 * arbitrary object.
 * ************************************************************** */

int
class_logical_scalar (Ent * e)
{
  int type=UNDEF;
  int rval;
  int (*vfptr) ();

  if (e)
    type = ent_type (e);
  vfptr = (IFPTR) logical_scalar_method[type].op;

  if (vfptr == 0)
  {
    fprintf (stderr, "Entity type: %s\n", etd (e));
    rerror ("coerce-to-logical-scalar operation not supported");
  }

  rval = (int) (*vfptr) (ent_data (e));

  return (rval);
}

/* **************************************************************
 * Given an object, get its size. Return 0 for UNDEF elements.
 * ************************************************************** */

size_t
class_sizeof (Ent * e)
{
  int type=UNDEF;
//   int rtype;
  size_t (*vfptr) ();
  size_t rval;

  if (e)
    type = ent_type (e);

  if (type)
  {
//     rtype = sizeof_method[type].type;
    vfptr = (SFPTR) sizeof_method[type].op;

    if (vfptr == 0)
    {
      fprintf (stderr, "Entity type: %s\n", etd (e));
      rerror ("sizeof() operation not supported");
    }

    rval = (*vfptr) (ent_data (e));
  }
  else
    rval = 0;

  return (rval);
}

/* **************************************************************
 * Given an object, return character pointer if possible.
 * ************************************************************** */

char *
class_char_pointer (Ent * e)
{
  int type=UNDEF;
  char *rval=0;
  char *(*vfptr) ();

  if (e)
    type = ent_type (e);

  vfptr = (CFPTR) char_pointer_method[type].op;

  if (vfptr == 0)
  {
    fprintf (stderr, "Entity type: %s\n", etd (e));
    rerror ("coerce-to-character operation not supported");
  }

  rval = (char *) (*vfptr) (ent_data (e));

  return (rval);
}

/* **************************************************************
 * Given an object, return a double value, if possible...
 * ************************************************************** */
double class_double (Ent * e)
{
  int type=UNDEF;
  double *dval, dtmp;
  double *(*vfptr) ();

  if (e)
    type = ent_type (e);

  vfptr = (void *) double_method[type].op;

  if (vfptr == 0)
  {
    fprintf (stderr, "Entity type: %s\n", etd (e));
    rerror ("coerce-to-double operation not supported");
  }

  dval = (double *) (*vfptr) (ent_data (e));
  dtmp = *dval;
  GC_FREE (dval);

  return (dtmp);
}

/* **************************************************************
 * Given an object, return an integer value, if possible...
 * ************************************************************** */

int
class_int (Ent * e)
{
  int type=UNDEF;
  int *ival, itmp;
  int *(*vfptr) ();

  if (e)
    type = ent_type (e);

  vfptr = (void *) int_method[type].op;

  if (vfptr == 0)
  {
    fprintf (stderr, "Entity type: %s\n", etd (e));
    rerror ("coerce-to-integer operation not supported");
  }

  ival = (int *) (*vfptr) (ent_data (e));
  itmp = *ival;
  GC_FREE (ival);

  return (itmp);
}

Complex class_complex (Ent * e)
{
  Complex cval;
  MDC *c;

  if (!e)
    rerror ("coerce-to-double operation not supported");

  if (ent_type(e) == MATRIX_DENSE_REAL)
  {
    RE(cval) = class_double(e);
    IM(cval) = 0.0;
  }
  else if (ent_type(e) == MATRIX_DENSE_COMPLEX)
  {
    c = ent_data(e);
    cval = MdcV0(c,0);
  }
  else
    rerror ("class_complex: non-numeric argument");

  return cval;
}


/* **************************************************************
 * Given an object, return a Real matrix, if possible. The matrix
 * MUST be treated as READ ONLY. The user should NOT destroy or
 * modify it in any way, shape or form.
 * ************************************************************** */

MDR *
class_matrix_real (Ent * e)
{
  int type=UNDEF;
  void *(*vfptr) ();
  MDR *mreal;

  if (e)
    type = ent_type (e);

  vfptr = (VFPTR) matrix_real_method[type].op;

  if (vfptr == 0)
  {
    fprintf (stderr, "Entity type: %s\n", etd (e));
    rerror ("coerce-to-matrix-real operation not supported");
  }

  mreal = (MDR *) (*vfptr) (ent_data (e));
  return (mreal);
}

/* **************************************************************
 * Given an object, return a STRING matrix, if possible. The matrix
 * MUST be treated as READ ONLY. The user should NOT destroy or
 * modify it in any way, shape or form.
 * ************************************************************** */

MDS *
class_matrix_string (Ent * e)
{
  int type=UNDEF;
  void *(*vfptr) ();
  MDS *mstring;

  if (e)
    type = ent_type (e);

  vfptr = (VFPTR) matrix_string_method[type].op;

  if (vfptr == 0)
  {
    fprintf (stderr, "Entity type: %s\n", etd (e));
    rerror ("coerce-to-matrix-string operation not supported");
  }

  mstring = (MDS *) (*vfptr) (ent_data (e));
  return (mstring);
}

/* **************************************************************
 * Get the size of an object.
 * ************************************************************** */

int
class_size (Ent * e)
{
  int size, type=UNDEF;
  int (*vfptr) ();

  if (e)
    type = ent_type (e);

  vfptr = (IFPTR) size_method[type].op;

  if (vfptr == 0)
  {
    fprintf (stderr, "Entity type: %s\n", etd (e));
    rerror ("size operation not supported");
  }

  size = (int) (*vfptr) (ent_data (e));

  return (size);
}

/* **************************************************************
 * Append two objects. Operation originates from the syntax:
 *                [ obj1 , obj2 ]
 * ************************************************************** */
#undef  THIS_SOLVER
#define THIS_SOLVER "class_append"
Ent * class_append (Ent * e1, Ent * e2)
{
  Ent *new=0;
  MDE *m=0;
  int type1=UNDEF, type2=UNDEF, rtype;
  void *rval=0;
  void *(*vfptr) ();

  if (e1)
    type1 = ent_type (e1);
  if (e2)
    type2 = ent_type (e2);

  if ((type1 != MATRIX_DENSE_ENTITY) && (type2 != MATRIX_DENSE_ENTITY))
  {
    rtype = append_method[type1][type2].type;
    vfptr = (VFPTR) append_method[type1][type2].op;

    if (vfptr == 0)
    {
      fprintf (stderr, THIS_SOLVER ": Entity types: %s and %s\n", etd (e1), etd (e2));
      rerror (THIS_SOLVER ": " RLAB_ERROR_APPEND_FAILED "\n");
    }

    rval = (*vfptr) (ent_data (e1), ent_data (e2));

    new = ent_Create ();
    ent_data (new) = rval;
    ent_type (new) = rtype;
  }
  else
  {
    int nr1, nc1, nr2, nc2, nr, nc, i, j;

    if (type1 != MATRIX_DENSE_ENTITY)
    {
      nr1 = 1;
      nc1 = 1;
    }
    else
    {
      nr1 = MNR(ent_data(e1));
      nc1 = MNC(ent_data(e1));
    }
    if (type2 != MATRIX_DENSE_ENTITY)
    {
      nr2 = 1;
      nc2 = 1;
    }
    else
    {
      nr2 = MNR(ent_data(e2));
      nc2 = MNC(ent_data(e2));
    }

    nc = nc1 + nc2;

    if (nr1 == 0 || nr2==0)
    {
      nr = nr1 + nr2;
    }
    else if (nr1==nr2)
    {
      nr = nr1;
    }
    else
    {
      fprintf (stderr, THIS_SOLVER ": Entity types: %s and %s\n", etd (e1), etd (e2));
      rerror (THIS_SOLVER ": " RLAB_ERROR_APPEND_REQUIRES "\n");
    }

    m = mde_CreateEmpty(nr,nc);
    if (type1 == MATRIX_DENSE_ENTITY)
    {
      MDE * edt = ent_data(e1);
      for (i=0; i<nr1; i++)
        for (j=0; j<nc1; j++)
          Mde0(m,i,j) = ent_Copy(Mde0(edt,i,j));
    }
    else
    {
      Mde0(m,0,0) = ent_Copy(e1);
    }

    if (type2 == MATRIX_DENSE_ENTITY)
    {
      MDE * edt = ent_data(e2);
      for (i=0; i<nr2; i++)
        for (j=0; j<nc2; j++)
          Mde0(m,i,j+nc1) = ent_Copy(Mde0(edt,i,j));
    }
    else
    {
      Mde0(m,0,nc1) = ent_Copy(e2);
    }

    new = ent_Create();
    ent_data(new) = m;
    ent_type(new) = MATRIX_DENSE_ENTITY;
  }

  return (new);
}

/* **************************************************************
 * Stack two objects. Operation originates from the syntax:
 *                [ obj1 ; obj2 ]
 * ************************************************************** */
#undef  THIS_SOLVER
#define THIS_SOLVER "class_stack"
Ent * class_stack (Ent * e1, Ent * e2)
{
  Ent *new=0;
  MDE *m=0;
  int type1=UNDEF, type2=UNDEF, rtype;
  void *rval=0;
  void *(*vfptr) ();

  if (e1)
    type1 = ent_type (e1);
  if (e2)
    type2 = ent_type (e2);

  if ((type1 != MATRIX_DENSE_ENTITY) && (type2 != MATRIX_DENSE_ENTITY))
  {
    rtype = stack_method[type1][type2].type;
    vfptr = (VFPTR) stack_method[type1][type2].op;

    if (vfptr == 0)
    {
      fprintf (stderr, THIS_SOLVER ": Entity types: %s and %s\n", etd (e1), etd (e2));
      rerror (THIS_SOLVER ": " RLAB_ERROR_STACK_FAILED "\n");
    }

    rval = (*vfptr) (ent_data (e1), ent_data (e2));

    new = ent_Create ();
    ent_data (new) = rval;
    ent_type (new) = rtype;
  }
  else
  {
    int nr1, nc1, nr2, nc2, nr, nc, i, j;

    if (type1 != MATRIX_DENSE_ENTITY)
    {
      nr1 = 1;
      nc1 = 1;
    }
    else
    {
      nr1 = MNR(ent_data(e1));
      nc1 = MNC(ent_data(e1));
    }
    if (type2 != MATRIX_DENSE_ENTITY)
    {
      nr2 = 1;
      nc2 = 1;
    }
    else
    {
      nr2 = MNR(ent_data(e2));
      nc2 = MNC(ent_data(e2));
    }

    nr = nr1 + nr2;

    if (nc1 == 0 || nc2==0)
    {
      nc = nc1 + nc2;
    }
    else if (nc1==nc2)
    {
      nc = nc1;
    }
    else
    {
      fprintf (stderr, THIS_SOLVER ": Entity types: %s and %s\n", etd (e1), etd (e2));
      rerror (THIS_SOLVER ": " RLAB_ERROR_STACK_REQUIRES "\n");
    }

    m = mde_CreateEmpty(nr,nc);

    if (type1 == MATRIX_DENSE_ENTITY)
    {
      MDE * edt = ent_data(e1);
      for (i=0; i<nr1; i++)
        for (j=0; j<nc1; j++)
          Mde0(m,i,j) = ent_Copy(Mde0(edt,i,j));
    }
    else
    {
      Mde0(m,0,0) = ent_Copy(e1);
    }

    if (type2 == MATRIX_DENSE_ENTITY)
    {
      MDE * edt = ent_data(e2);
      for (i=0; i<nr2; i++)
        for (j=0; j<nc2; j++)
          Mde0(m,i+nr1,j) = ent_Copy(Mde0(edt,i,j));
    }
    else
    {
      Mde0(m,nr1,0) = ent_Copy(e2);
    }

    new = ent_Create();
    ent_data(new) = m;
    ent_type(new) = MATRIX_DENSE_ENTITY;
  }


  return (new);
}

/* **************************************************************
 * Evaluate a sub-matrix expression:
 *            VAR  [ exprI ; exprJ ]
 *
 * exprI, and exprJ could be almost anything, so before we figure
 * out who is going to evaluate the expression, we must coerce
 * exprI, and exprJ into integer vectors. Why integer vectors?
 * I am assuming that all matrices/arrays will ultimately be indexed
 * with integers.
 * ************************************************************** */
#undef  THIS_SOLVER
#define THIS_SOLVER "class_matrix_sub_1"
Ent * class_matrix_sub_1 (Ent * var, Ent * ei, Ent * ej)
{
  Ent *new=0;
  MDE *x=0, *v=0;
  int type=UNDEF, itype, jtype, rtype;
  int *i, *j, m, n;
  void *rval=0;
  void *(*vfptr) (), *(*ifptr) (), *(*jfptr) ();

#ifdef DEBUG_MDE
  fprintf (stderr, THIS_SOLVER ": Entity type: %s\n", etd (var));
#endif

  itype = ent_type (ei);
  ifptr = (VFPTR) ind_coerce_int_method[type][itype].op;
  if (ifptr == 0)
  {
    fprintf (stderr, THIS_SOLVER ": Entity type: %s and %s\n", etd (var), etd (ei));
    rerror (THIS_SOLVER ": " RLAB_ERROR_SUBMAT_INDEX_COERCION_FAILED "\n");
  }
  i = (int *) (*ifptr) (ent_data (var), ent_data (ei));

  jtype = ent_type (ej);
  jfptr = (VFPTR) ind_coerce_int_method[type][jtype].op;
  if (jfptr == 0)
  {
    fprintf (stderr, THIS_SOLVER ": Entity type: %s and %s\n", etd (var), etd (ej));
    rerror (THIS_SOLVER ": " RLAB_ERROR_SUBMAT_INDEX_COERCION_FAILED "\n");
  }
  j = (int *) (*jfptr) (ent_data (var), ent_data (ej));

  type = ent_type (var);
  if (type != MATRIX_DENSE_ENTITY)
  {
    vfptr = (VFPTR) matrix_sub_method[type].op;
    if (vfptr == 0)
    {
      fprintf (stderr, THIS_SOLVER  ": Entity type: %s\n", etd (var));
      rerror (THIS_SOLVER ": " RLAB_ERROR_SUBMAT_ASSIGNMENT_FAILED "\n");
    }
    rval = (*vfptr) (ent_data (var), i, j, &rtype);

    new = ent_Create ();
    ent_data (new) = rval;
    ent_type (new) = rtype;
  }
  else
  {
    v = ent_data(var);
    if (i[0]==1 && j[0]==1)
    {
      if (i[1] > MNR (v))
      {
        fprintf (stderr, "\t" RLAB_ERROR_INDEX_OUT_OF_BOUNDS "\n");
        fprintf (stderr, "\tindex value: %i, matrix size: %i\n", i[1], MNR (v));
        rerror (RLAB_ERROR_SUBMAT_EVAL_FAILED);
      }
      if (j[1] > MNC (v))
      {
        fprintf (stderr, "\t" RLAB_ERROR_INDEX_OUT_OF_BOUNDS "\n");
        fprintf (stderr, "\tindex value: %i, matrix size: %i\n", j[1], MNC (v));
        rerror (RLAB_ERROR_SUBMAT_EVAL_FAILED);
      }

//       fprintf(stderr, "refc = %i\n", ent_Ref(Mde1 (v,i[1],j[1])));
      new = Mde1 (v,i[1],j[1]);
    }
    else
    {
      x = mde_Create (i[0], j[0]);
      for (m = 1; m <= i[0]; m++)
      {
        if (i[m] > MNR (v))
        {
          mde_Destroy (x);
          fprintf (stderr, "\t" RLAB_ERROR_INDEX_OUT_OF_BOUNDS "\n");
          fprintf (stderr, "\tindex value: %i, matrix size: %i\n", i[m], MNR (v));
          rerror (RLAB_ERROR_SUBMAT_EVAL_FAILED);
        }
        for (n = 1; n <= j[0]; n++)
        {
          if (j[n] > MNC (v))
          {
            mde_Destroy (x);
            fprintf (stderr, "\t" RLAB_ERROR_INDEX_OUT_OF_BOUNDS "\n");
            fprintf (stderr, "\tindex value: %i, matrix size: %i\n", j[n], MNC (v));
            rerror (RLAB_ERROR_SUBMAT_EVAL_FAILED);
          }
          Mde1 (x, m, n) = Mde1 (v, i[m], j[n]);
        }
      }
      new = ent_Create();
      ent_data(new) = x;
      ent_type(new) = MATRIX_DENSE_ENTITY;
    }
  }

  if (i)
    GC_FREE (i);
  if (j)
    GC_FREE (j);

  return (new);
}

Ent *
class_matrix_sub_2 (Ent * var, Ent * ei)
{
  Ent *new;
  int type=UNDEF, itype, rtype;
  int *i;
  void *rval=0;
  void *(*vfptr) (), *(*ifptr) ();

  type = ent_type (var);
  vfptr = (VFPTR) matrix_sub_r_method[type].op;

  if (vfptr == 0)
  {
    fprintf (stderr, "Entity type: %s\n", etd (var));
    rerror ("sub-matrix operation not supported");
  }

  /*
   * First we must try and coerce the indices.
   */

  itype = ent_type (ei);
  ifptr = (VFPTR) ind_coerce_int_method[type][itype].op;
  if (ifptr == 0)
  {
    fprintf (stderr, "Entity type: %s and %s\n", etd (var), etd (ei));
    rerror ("sub-matrix index coercion operation not supported");
  }

  i = (int *) (*ifptr) (ent_data (var), ent_data (ei));

  /*
   * Now we can get the sub-matrix.
   */

  rval = (*vfptr) (ent_data (var), i, &rtype);

  new = ent_Create ();
  ent_data (new) = rval;
  ent_type (new) = rtype;

  if (i)
    GC_FREE (i);

  return (new);
}

Ent *
class_matrix_sub_3 (Ent * var, Ent * ej)
{
  Ent *new;
  int type=UNDEF, jtype, rtype;
  int *j;
  void *rval=0;
  void *(*vfptr) (), *(*jfptr) ();

  type = ent_type (var);
  vfptr = (VFPTR) matrix_sub_c_method[type].op;

  if (vfptr == 0)
  {
    fprintf (stderr, "Entity type: %s\n", etd (var));
    rerror ("sub-matrix operation not supported");
  }

  /*
   * First we must try and coerce the indices.
   */

  jtype = ent_type (ej);
  jfptr = (VFPTR) ind_coerce_int_method[type][jtype].op;
  if (jfptr == 0)
  {
    fprintf (stderr, "Entity type: %s and %s\n", etd (var), etd (ej));
    rerror ("sub-matrix index coercion operation not supported");
  }

  j = (int *) (*jfptr) (ent_data (var), ent_data (ej));

  /*
   * Now we can get the sub-matrix.
   */

  rval = (*vfptr) (ent_data (var), j, &rtype);

  new = ent_Create ();
  ent_data (new) = rval;
  ent_type (new) = rtype;

  if (j)
    GC_FREE (j);

  return (new);
}

/* **************************************************************
 * Handle assignment of a matrix:
 * var [ ei ; eij ] = erhs
 * ************************************************************** */
#undef THIS_SOLVER
#define THIS_SOLVER "class_matrix_assign_1"
Ent * class_matrix_assign_1 (Ent * var, Ent * ei, Ent * ej, Ent * erhs)
{
  int lhs_type=UNDEF, itype=UNDEF, jtype=UNDEF, rtype=UNDEF, rhs_type=UNDEF;
  int *i, *j;
  void *rval=0;
  void *(*vfptr) (), *(*ifptr) (), *(*jfptr) ();
  FILE *rlab_stderr = (!RLAB_STDERR_DS) ? stderr : RLAB_STDERR_DS;

  // Get the type of the assignment
  if (var)
    lhs_type = ent_type (var);
  if (erhs)
    rhs_type = ent_type (erhs);

#ifdef DEBUG_MDE
  fprintf (rlab_stderr, THIS_SOLVER ": Entity types: %s and %s\n", etd (var), etd (erhs));
#endif

  // Now we must try and coerce the indices.
  if (ei)
    itype = ent_type (ei);
  ifptr = (VFPTR) ind_coerce_int_method[lhs_type][itype].op;
  if (ifptr == 0)
  {
    fprintf (rlab_stderr, THIS_SOLVER ": Entity type: %s and %s\n", etd (var), etd (ei));
    rerror (THIS_SOLVER ": " RLAB_ERROR_SUBMAT_INDEX_COERCION_FAILED "\n");
  }
  i = (int *) (*ifptr) (ent_data (var), ent_data (ei));

  // Now we must try and coerce the indices.
  if (ej)
    jtype = ent_type (ej);
  jfptr = (VFPTR) ind_coerce_int_method[lhs_type][jtype].op;
  if (jfptr == 0)
  {
    fprintf (rlab_stderr, "Entity type: %s and %s\n", etd (var), etd (ej));
    rerror (THIS_SOLVER ": " RLAB_ERROR_SUBMAT_INDEX_COERCION_FAILED "\n");
  }
  j = (int *) (*jfptr) (ent_data (var), ent_data (ej));

  if (lhs_type != MATRIX_DENSE_ENTITY)
  {
    rtype = matrix_assign_method[lhs_type][rhs_type].type;
    vfptr = (VFPTR) matrix_assign_method[lhs_type][rhs_type].op;

    if (vfptr == 0)
    {
      fprintf (rlab_stderr, "Entity types: %s and %s\n", etd (var), etd (erhs));
      rerror (THIS_SOLVER ": " RLAB_ERROR_SUBMAT_ASSIGNMENT_FAILED "\n");
    }

    // Now we can do the assignment.
    rval = (*vfptr) (ent_data (var), i, j, ent_data (erhs));

    ent_data (var) = rval;
    ent_type (var) = rtype;
  }
  else
  {
    //
    // LHS is matrix dense entity:
    //  there are two possibilities:
    //    LHS is single entry, and RHS is not MDE
    //  and
    //    LHS 
    //  x[i;j] = RHS
    MDE * edummy = ent_data(var);
    int i1, j1;

    if (i[0]==1 && j[0]==1)
    {
      ent_Destroy( Mde1(edummy,i[1],j[1]) );
      if (ent_type(erhs)!=UNDEF)
        Mde1(edummy,i[1],j[1]) = ent_Copy(erhs);
      else
        Mde1(edummy,i[1],j[1]) = ent_Create();
    }
    else
    {
      if (rhs_type == MATRIX_DENSE_ENTITY)
      {
        MDE * rdummy = ent_data(erhs);
        for (i1=1; i1<=i[0]; i1++)
        {
          for (j1=1; j1<=j[0]; j1++)
          {
            ent_Destroy(Mde1(edummy,i[i1], j[j1]));
            if (ent_type(Mde1(rdummy,i1,j1))!=UNDEF)
              Mde1(edummy,i[i1], j[j1]) = ent_Copy(Mde1(rdummy,i1,j1));
            else
              Mde1(edummy,i[i1], j[j1]) = ent_Create();
          }
        }
      }
      else
      {
        fprintf (rlab_stderr, "Entity types: %s and %s\n", etd (var), etd (erhs));
        rerror (THIS_SOLVER ": " RLAB_ERROR_SUBMAT_ASSIGNMENT_FAILED "\n");
      }
    }
  }

  if (i)
    GC_FREE (i);
  if (j)
    GC_FREE (j);

  return (var);
}

#undef THIS_SOLVER
#define THIS_SOLVER "class_matrix_assign_2"
Ent * class_matrix_assign_2 (Ent * var, Ent * ei, Ent * erhs)
{
  int type=UNDEF, itype=UNDEF, rtype, rhs_type=UNDEF;
  int *i;
  void *rval=0;
  void *(*vfptr) (), *(*ifptr) ();
  FILE *rlab_stderr = (!RLAB_STDERR_DS) ? stderr : RLAB_STDERR_DS;

  // Get the type of the assignment
  if (var)
    type = ent_type (var);
  if (erhs)
    rhs_type = ent_type (erhs);

#ifdef DEBUG_MDE
  fprintf (rlab_stderr, "Entity types: %s and %s\n", etd (var), etd (erhs));
#endif

  // Now we must try and coerce the indices.
  if (ei)
    itype = ent_type (ei);
  ifptr = (VFPTR) ind_coerce_int_method[type][itype].op;
  if (ifptr == 0)
  {
    fprintf (rlab_stderr, THIS_SOLVER ": Entity type: %s and %s\n", etd (var), etd (ei));
    rerror (THIS_SOLVER ": " RLAB_ERROR_SUBMAT_INDEX_COERCION_FAILED "\n");
  }
  i = (int *) (*ifptr) (ent_data (var), ent_data (ei));

  //
  // Get the type of the assignment.
  //
  if (type != MATRIX_DENSE_ENTITY)
  {

    rtype = matrix_assign_r_method[type][rhs_type].type;
    vfptr = (VFPTR) matrix_assign_r_method[type][rhs_type].op;

    if (vfptr == 0)
    {
      fprintf (stderr, THIS_SOLVER ": Entity types: %i and %i\n", type, rhs_type);
      rerror (THIS_SOLVER ": " RLAB_ERROR_SUBMAT_EVAL_FAILED "\n");
    }

    rval = (*vfptr) (ent_data (var), i, ent_data (erhs));

    ent_data (var) = rval;
    ent_type (var) = rtype;
  }
//   else
//   {
//     //
//     // LHS is matrix dense entity:
//     //  there are two possibilities:
//     //    LHS is single entry, and RHS is not MDE
//     //  and
//     //    LHS 
//     //  x[i;j] = RHS
//     MDE * edummy = ent_data(var);
//     int i1;
// 
//     if (i[0]==1)
//     {
//       ent_Destroy( MdeV1(edummy,i[1]) );
//       if (ent_type(erhs)!=UNDEF)
//         MdeV1(edummy,i[1]) = ent_Copy(erhs);
//       else
//         MdeV1(edummy,i[1]) = ent_Create();
//     }
//     else
//     {
//       if (rhs_type == MATRIX_DENSE_ENTITY)
//       {
//         MDE * rdummy = ent_data(erhs);
//         for (i1=1; i1<=i[0]; i1++)
//         {
//           ent_Destroy(MdeV1(edummy,i[i1]));
//           if (ent_type(MdeV1(rdummy,i1))!=UNDEF)
//             MdeV1(edummy,i[i1]) = ent_Copy(MdeV1(rdummy,i1));
//           else
//             MdeV1(edummy,i[i1]) = ent_Create();
//         }
//       }
//       else
//       {
//         fprintf (stderr, THIS_SOLVER ": Entity types: %i and %i\n", type, rhs_type);
//         rerror (THIS_SOLVER ": " RLAB_ERROR_SUBMAT_EVAL_FAILED "\n");
//       }
//     }
//   }

  if (i)
    GC_FREE (i);
  return (var);
}

#undef THIS_SOLVER
#define THIS_SOLVER "class_matrix_assign_3"
Ent * class_matrix_assign_3 (Ent * var, Ent * ej, Ent * erhs)
{
  int type=UNDEF, jtype=UNDEF, rtype=UNDEF, rhs_type=UNDEF;
  int *j;
  void *rval=0;
  void *(*vfptr) (), *(*jfptr) ();
  FILE *rlab_stderr = (!RLAB_STDERR_DS) ? stderr : RLAB_STDERR_DS;

  // Get the type of the assignment
  if (var)
    type = ent_type (var);
  if (erhs)
    rhs_type = ent_type (erhs);

#ifdef DEBUG_MDE
  fprintf (rlab_stderr, "Entity types: %s and %s\n", etd (var), etd (erhs));
#endif

  // Now we must try and coerce the indices.
  if (ej)
    jtype = ent_type (ej);
  jfptr = (VFPTR) ind_coerce_int_method[type][jtype].op;
  if (jfptr == 0)
  {
    fprintf (rlab_stderr, THIS_SOLVER ": Entity type: %s and %s\n", etd (var), etd (ej));
    rerror (THIS_SOLVER ": " RLAB_ERROR_SUBMAT_EVAL_FAILED "\n");
  }
  j = (int *) (*jfptr) (ent_data (var), ent_data (ej));

  if (type != MATRIX_DENSE_ENTITY)
  {
    rtype = matrix_assign_c_method[type][rhs_type].type;
    vfptr = (VFPTR) matrix_assign_c_method[type][rhs_type].op;

    if (vfptr == 0)
    {
      fprintf (rlab_stderr, THIS_SOLVER ": Entity type: %s and %s\n", etd (var), etd (ej));
      rerror (THIS_SOLVER ": " RLAB_ERROR_SUBMAT_EVAL_FAILED "\n");
    }

    rval = (*vfptr) (ent_data (var), j, ent_data (erhs));

    ent_data (var) = rval;
    ent_type (var) = rtype;
  }
//   else
//   {
//     //
//     // LHS is matrix dense entity:
//     //  there are two possibilities:
//     //    LHS is single entry, and RHS is not MDE
//     //  and
//     //    LHS 
//     //  x[i;j] = RHS
//     MDE * edummy = ent_data(var);
//     int j1;
// 
//     if (j[0]==1)
//     {
//       ent_Destroy( MdeV1(edummy,j[1]) );
//       if (ent_type(erhs)!=UNDEF)
//         MdeV1(edummy,j[1]) = ent_Copy(erhs);
//       else
//         MdeV1(edummy,j[1]) = ent_Create();
//     }
//     else
//     {
//       if (rhs_type == MATRIX_DENSE_ENTITY)
//       {
//         MDE * rdummy = ent_data(erhs);
//         for (j1=1; j1<=j[0]; j1++)
//         {
//           ent_Destroy(MdeV1(edummy,j[j1]));
//           MdeV1(edummy,j[j1]) = ent_Copy(erhs);
// 
//           if (ent_type(MdeV1(rdummy,i1))!=UNDEF)
//             MdeV1(edummy,i[i1]) = ent_Copy(MdeV1(rdummy,i1));
//           else
//             MdeV1(edummy,i[i1]) = ent_Create();
// 
//         }
//       }
//       else
//       {
//         fprintf (rlab_stderr, THIS_SOLVER ": Entity type: %s and %s\n", etd (var), etd (ej));
//         rerror (THIS_SOLVER ": " RLAB_ERROR_SUBMAT_EVAL_FAILED "\n");
//       }
//     }
//   }


  if (j)
    GC_FREE (j);

  return (var);
}

/* *******************************************************
 * Evaluate a sub-matrix (vector) expression.
 * ******************************************************* */
#undef THIS_SOLVER
#define THIS_SOLVER "class_vector_sub"
Ent *
class_vector_sub (Ent * var, Ent * eind)
{
  Ent *new=0;
  int type=UNDEF, itype, rtype;
  int *i;
  void *rval=0;
  void *(*vfptr) (), *(*ifptr) ();

  if (var)
    type = ent_type (var);

#ifdef DEBUG_MDE
  fprintf (stderr, THIS_SOLVER ": Entity type: %s\n", etd (var));
#endif

  // First we must try and coerce the indices.
  itype = ent_type (eind);
  ifptr = (VFPTR) ind_coerce_int_method[type][itype].op;
  if (ifptr == 0)
  {
    fprintf (stderr, THIS_SOLVER ": Entity type: %s and %s\n", etd (var), etd (eind));
    rerror (THIS_SOLVER ": " RLAB_ERROR_SUBVEC_INDEX_COERCION_FAILED "\n");
  }
  i = (int *) (*ifptr) (ent_data (var), ent_data (eind));

  if (type != MATRIX_DENSE_ENTITY)
  {
    vfptr = (VFPTR) matrix_vector_sub_method[type].op;
    if (vfptr == 0)
    {
      fprintf (stderr, THIS_SOLVER ": Entity type: %s\n", etd (var));
      rerror (THIS_SOLVER ": " RLAB_ERROR_SUBVEC_OPERATION_FAILED "\n");
    }
    rval = (*vfptr) (ent_data (var), i, &rtype);

    new = ent_Create ();
    ent_data (new) = rval;
    ent_type (new) = rtype;
  }
  else
  {
    //
    // LHS is matrix dense entity:
    //  there are two possibilities:
    //    LHS is single entry, and RHS is not MDE
    //  and
    //    LHS 
    //  x[i;j] = RHS
    int i1;
    MDE * edummy = ent_data(var);

    if (i[0]==1)
    {
      if (SIZE(edummy)>= i[1])
      {
        new = ent_Copy(MdeV1(edummy,i[1]));
      }
      else
      {
        rerror (THIS_SOLVER ": " RLAB_ERROR_INDEX_OUT_OF_BOUNDS ": " RLAB_ERROR_SUBVEC_OPERATION_FAILED "\n");
      }
    }
    else
    {
      new = ent_Create();
      ent_data(new) = mde_Create(1, i[0]);
      ent_type(new) = MATRIX_DENSE_ENTITY;
      MDE * enew  = ent_data(new);

      for (i1=1; i1<=i[0]; i1++)
      {
        if (SIZE(enew)>= i[i1])
        {
          MdeV1(enew,i[i1]) = ent_Copy(MdeV1(edummy,i1));
        }
        else
        {
          rerror (THIS_SOLVER ": " RLAB_ERROR_INDEX_OUT_OF_BOUNDS ": " RLAB_ERROR_SUBVEC_OPERATION_FAILED "\n");
        }
      }
    }
  }

  if (i)
    GC_FREE (i);

  return (new);
}

/* **************************************************************
 * Assign to a matrix, like it was a vector.
 * var[i] = rhs
 * ************************************************************** */
#undef  THIS_SOLVER
#define THIS_SOLVER "class_matrix_vector_assign"
Ent *
class_matrix_vector_assign (Ent * var, Ent * ei, Ent * erhs)
{
  int type=UNDEF, itype=UNDEF, rtype=UNDEF, rhs_type=UNDEF;
  int *i;
  void *rval=0;
  void *(*vfptr) (), *(*ifptr) ();

  if (var)
    type = ent_type (var);
  if (erhs)
    rhs_type = ent_type (erhs);

#ifdef DEBUG_MDE
  fprintf (stderr, THIS_SOLVER ": Entity type: %s and %s\n", etd (var), etd (erhs));
#endif

  if (ei)
    itype = ent_type (ei);
  ifptr = (VFPTR) ind_coerce_int_method[type][itype].op;
  if (ifptr == 0)
  {
    fprintf (stderr, THIS_SOLVER ": Entity type: %s and %s\n", etd (var), etd (ei));
    rerror (THIS_SOLVER ": " RLAB_ERROR_SUBMAT_INDEX_COERCION_FAILED "\n");
  }
  i = (int *) (*ifptr) (ent_data (var), ent_data (ei));

  if (type != MATRIX_DENSE_ENTITY)
  {
    rtype = matrix_assign_vector_method[type][rhs_type].type;
    vfptr = (VFPTR) matrix_assign_vector_method[type][rhs_type].op;

    if (vfptr == 0)
    {
      fprintf (stderr, THIS_SOLVER  ": Entity types: %s and %s\n", etd (var), etd (erhs));
      rerror (THIS_SOLVER ": " RLAB_ERROR_INDEX_OUT_OF_BOUNDS ": " RLAB_ERROR_SUBMATVEC_ASSIGNMENT_FAILED "\n");
    }

    rval = (*vfptr) (ent_data (var), i, ent_data (erhs));
    ent_data (var) = rval;
    ent_type (var) = rtype;
  }
  else
  {
    int i1;

    MDE * edummy = ent_data(var);
    if (i[0]==1)
    {
      if (i[1] <= SIZE(edummy))
      {
        ent_Destroy(MdeV1(edummy,i[1]));
        MdeV1(edummy,i[1]) = ent_Copy(erhs);
      }
      else
        rerror (THIS_SOLVER ": " RLAB_ERROR_INDEX_OUT_OF_BOUNDS ": " RLAB_ERROR_SUBMATVEC_ASSIGNMENT_FAILED "\n");
    }
    else
    {
      if (rhs_type == MATRIX_DENSE_ENTITY)
      {
        MDE * edummy_rhs = ent_data(erhs);
        for (i1=1; i1<=i[0]; i1++)
        {
          if (i[i1] <= SIZE(edummy_rhs))
          {
            ent_Destroy(MdeV1(edummy,i[i1]));
            MdeV1(edummy,i[i1]) = ent_Copy(MdeV1(edummy_rhs,i[i1]));
          }
        }
      }
      else
      {
        fprintf (stderr, THIS_SOLVER  ": Entity types: %s and %s\n", etd (var), etd (erhs));
        rerror (THIS_SOLVER ": " RLAB_ERROR_INDEX_OUT_OF_BOUNDS ": " RLAB_ERROR_SUBMATVEC_ASSIGNMENT_FAILED "\n");
      }
    }
  }

  if (i)
    GC_FREE (i);

  return (var);
}

// **************************************************************
// Get the i'th value of an object.
//  modification (mk): return 0 pointer if it fails so that the calling
//  routine in code.c calls the error.
// **************************************************************
#undef  THIS_SOLVER
#define THIS_SOLVER "class_forloop_value"
Ent * class_forloop_value (Ent * e, int i)
{
  Ent *new=0;
  int rtype=UNDEF, type=UNDEF;
  void *rval=0;
  void *(*vfptr) ();

  if (e)
    type = ent_type (e);

  if (type != MATRIX_DENSE_ENTITY)
  {
    vfptr = (VFPTR) forloop_value_method[type].op;
    if (vfptr)
    {
      rtype = forloop_value_method[type].type;
      rval = (*vfptr) (ent_data (e), i);

      new = ent_Create ();
      ent_data (new) = rval;
      ent_type (new) = rtype;
    }
    else
    {
      fprintf (stderr, THIS_SOLVER ": " RLAB_ERROR_FORLOOPVAL_FAILED " '%s'\n", etd (e));
    }
  }
  else
  {
    MDE * edummy = ent_data(e);
    if (SIZE(edummy)>=i)
    {
      new = ent_Copy(MdeV1(edummy,i));
    }
    else
      new = ent_Create();

    rtype = ent_type(new);
  }

  return (new);
}

/* **************************************************************
 * Return an empty something. For now, just return an empty
 * dense-real-matrix. Go to all this trouble (using classes) to
 * make it easier to fix it later.
 * ************************************************************** */

Ent *
class_empty (void)
{
  Ent *new;
  int rtype, type;
  void *rval=0;
  void *(*vfptr) ();

  type = MATRIX_DENSE_REAL;

  rtype = empty_method[type].type;
  vfptr = (VFPTR) empty_method[type].op;

  if (vfptr == 0)
  {
    rerror ("empty operation not supported");
  }

  rval = (*vfptr) ();

  new = ent_Create ();
  ent_data (new) = rval;
  ent_type (new) = rtype;

  return (new);
}

/* **************************************************************
 * Check the object to see if it has a read-only member/attribute
 * with name `name'. If it does, return it. Other wise, return NULL
 * ************************************************************** */

void *
class_member_ref (Ent * e, char *name, int *rtype)
{
  int type=UNDEF;

  void *(*vfptr) ();

  void *ne=0;

  if (e)
    type = ent_type (e);

  vfptr = (VFPTR) member_ref_method[type].op;

  if (vfptr == 0)
  {
    fprintf (stderr, "Entity type: %s\n", etd (e));
    rerror ("member-ref operation not supported");
  }

  ne = (void *) (*vfptr) (ent_data (e), name, rtype);

  return (ne);
}

/* **************************************************************
 * Compare (==) two objects.
 * ************************************************************** */

Ent *
class_eq (Ent * e1, Ent * e2)
{
  Ent *new;
  int type1=UNDEF, type2=UNDEF, rtype;
  void *rval=0;
  void *(*vfptr) ();

  if (e1)
    type1 = ent_type (e1);
  if (e2)
    type2 = ent_type (e2);

  rtype = eq_method[type1][type2].type;
  vfptr = (VFPTR) eq_method[type1][type2].op;

  if (vfptr == 0)
  {
    fprintf (stderr, "Entity types: %s and %s\n", etd (e1), etd (e2));
    rerror ("equality operation not supported");
  }

  rval = (*vfptr) (ent_data (e1), ent_data (e2));

  new = ent_Create ();
  ent_data (new) = rval;
  ent_type (new) = rtype;

  return (new);
}

Ent *
class_ne (Ent * e1, Ent * e2)
{
  Ent *new;
  int type1=UNDEF, type2=UNDEF, rtype;
  void *rval=0;
  void *(*vfptr) ();

  if (e1)
    type1 = ent_type (e1);
  if (e2)
    type2 = ent_type (e2);

  rtype = ne_method[type1][type2].type;
  vfptr = (VFPTR) ne_method[type1][type2].op;

  if (vfptr == 0)
  {
    fprintf (stderr, "Entity types: %s and %s\n", etd (e1), etd (e2));
    rerror ("in-equality operation not supported");
  }

  rval = (*vfptr) (ent_data (e1), ent_data (e2));

  new = ent_Create ();
  ent_data (new) = rval;
  ent_type (new) = rtype;

  return (new);
}

Ent *
class_lt (Ent * e1, Ent * e2)
{
  Ent *new;
  int type1=UNDEF, type2=UNDEF, rtype;
  void *rval=0;
  void *(*vfptr) ();

  if (e1)
    type1 = ent_type (e1);
  if (e2)
    type2 = ent_type (e2);

  rtype = lt_method[type1][type2].type;
  vfptr = (VFPTR) lt_method[type1][type2].op;

  if (vfptr == 0)
  {
    fprintf (stderr, "Entity types: %s and %s\n", etd (e1), etd (e2));
    rerror ("less-than operation not supported");
  }

  rval = (*vfptr) (ent_data (e1), ent_data (e2));

  new = ent_Create ();
  ent_data (new) = rval;
  ent_type (new) = rtype;

  return (new);
}

Ent *
class_le (Ent * e1, Ent * e2)
{
  Ent *new;
  int type1=UNDEF, type2=UNDEF, rtype;
  void *rval=0;
  void *(*vfptr) ();

  if (e1)
    type1 = ent_type (e1);
  if (e2)
    type2 = ent_type (e2);

  rtype = le_method[type1][type2].type;
  vfptr = (VFPTR) le_method[type1][type2].op;

  if (vfptr == 0)
  {
    fprintf (stderr, "Entity types: %s and %s\n", etd (e1), etd (e2));
    rerror ("less-than-or-equal operation not supported");
  }

  rval = (*vfptr) (ent_data (e1), ent_data (e2));

  new = ent_Create ();
  ent_data (new) = rval;
  ent_type (new) = rtype;

  return (new);
}

Ent *
class_gt (Ent * e1, Ent * e2)
{
  Ent *new;
  int type1=UNDEF, type2=UNDEF, rtype;
  void *rval=0;
  void *(*vfptr) ();

  if (e1)
    type1 = ent_type (e1);
  if (e2)
    type2 = ent_type (e2);

  rtype = gt_method[type1][type2].type;
  vfptr = (VFPTR) gt_method[type1][type2].op;

  if (vfptr == 0)
  {
    fprintf (stderr, "Entity types: %s\t and %s\n", etd (e1), etd (e2));
    rerror ("greater-than operation not supported");
  }

  rval = (*vfptr) (ent_data (e1), ent_data (e2));

  new = ent_Create ();
  ent_data (new) = rval;
  ent_type (new) = rtype;

  return (new);
}

Ent *
class_ge (Ent * e1, Ent * e2)
{
  Ent *new;
  int type1=UNDEF, type2=UNDEF, rtype;
  void *rval=0;
  void *(*vfptr) ();

  if (e1)
    type1 = ent_type (e1);
  if (e2)
    type2 = ent_type (e2);

  rtype = ge_method[type1][type2].type;
  vfptr = (VFPTR) ge_method[type1][type2].op;

  if (vfptr == 0)
  {
    fprintf (stderr, "Entity types: %s\t and %s\n", etd (e1), etd (e2));
    rerror ("greater-than-or-equal operation not supported");
  }

  rval = (*vfptr) (ent_data (e1), ent_data (e2));

  new = ent_Create ();
  ent_data (new) = rval;
  ent_type (new) = rtype;

  return (new);
}

Ent *
class_and (Ent * e1, Ent * e2)
{
  Ent *new;
  int type1=UNDEF, type2=UNDEF, rtype;
  void *rval=0;
  void *(*vfptr) ();

  if (e1)
    type1 = ent_type (e1);
  if (e2)
    type2 = ent_type (e2);

  rtype = and_method[type1][type2].type;
  vfptr = (VFPTR) and_method[type1][type2].op;

  if (vfptr == 0)
  {
    fprintf (stderr, "Entity types: %s and %s\n", etd (e1), etd (e2));
    rerror ("and operation not supported");
  }

  rval = (*vfptr) (ent_data (e1), ent_data (e2));

  new = ent_Create ();
  ent_data (new) = rval;
  ent_type (new) = rtype;

  return (new);
}

Ent *
class_or (Ent * e1, Ent * e2)
{
  Ent *new;
  int type1=UNDEF, type2=UNDEF, rtype;
  void *rval=0;
  void *(*vfptr) ();

  if (e1)
    type1 = ent_type (e1);
  if (e2)
    type2 = ent_type (e2);

  rtype = or_method[type1][type2].type;
  vfptr = (VFPTR) or_method[type1][type2].op;

  if (vfptr == 0)
  {
    fprintf (stderr, "Entity types: %s and %s\n", etd (e1), etd (e2));
    rerror ("or operation not supported");
  }

  rval = (*vfptr) (ent_data (e1), ent_data (e2));

  new = ent_Create ();
  ent_data (new) = rval;
  ent_type (new) = rtype;

  return (new);
}

Ent *
class_not (Ent * e)
{
  Ent *new;
  int type=UNDEF, rtype;
  void *rval=0;
  void *(*vfptr) ();

  type = ent_type (e);

  rtype = not_method[type].type;
  vfptr = (VFPTR) not_method[type].op;

  if (vfptr == 0)
  {
    fprintf (stderr, "Entity type: %s\n", etd (e));
    rerror ("not operation not supported");
  }

  rval = (*vfptr) (ent_data (e));

  new = ent_Create ();
  ent_data (new) = rval;
  ent_type (new) = rtype;

  return (new);
}

/* **************************************************************
 * Transpose an object.
 * ************************************************************** */

Ent *
class_transpose (Ent * e)
{
  Ent *new;
  int type=UNDEF, rtype;
  void *rval=0;
  void *(*vfptr) ();

  type = ent_type (e);
  rtype = transpose_method[type].type;
  vfptr = (VFPTR) transpose_method[type].op;

  if (vfptr == 0)
  {
    fprintf (stderr, "Entity type: %s\n", etd (e));
    rerror ("transpose operation not supported");
  }

  rval = (*vfptr) (ent_data (e));

  new = ent_Create ();
  ent_data (new) = rval;
  ent_type (new) = rtype;

  return (new);
}

Ent *
class_nc_transpose (Ent * e)
{
  Ent *new;
  int type=UNDEF, rtype;
  void *rval=0;
  void *(*vfptr) ();

  type = ent_type (e);
  rtype = nc_transpose_method[type].type;
  vfptr = (VFPTR) nc_transpose_method[type].op;

  if (vfptr == 0)
  {
    fprintf (stderr, "Entity type: %s\n", etd (e));
    rerror ("non-conjugate transpose operation not supported");
  }

  rval = (*vfptr) (ent_data (e));

  new = ent_Create ();
  ent_data (new) = rval;
  ent_type (new) = rtype;

  return (new);
}

/* **************************************************************
 * Reshape an object into a column vector.
 * ************************************************************** */

Ent *
class_reshape_col (Ent * e)
{
  Ent *new;
  int type=UNDEF, rtype;
  void *rval=0;
  void *(*vfptr) ();

  type = ent_type (e);
  rtype = reshape_col_method[type].type;
  vfptr = (VFPTR) reshape_col_method[type].op;

  if (vfptr == 0)
  {
    fprintf (stderr, "Entity type: %s\n", etd (e));
    rerror ("column-reshape operation not supported");
  }

  rval = (*vfptr) (ent_data (e));

  new = ent_Create ();
  ent_data (new) = rval;
  ent_type (new) = rtype;

  return (new);
}

/* **************************************************************
 * Increment an object.
 * ************************************************************** */

Ent *
class_increment (Ent * e)
{
  int type=UNDEF, rtype;
  void *rval=0;
  void *(*vfptr) ();

  type = ent_type (e);
  rtype = increment_method[type].type;
  vfptr = (VFPTR) increment_method[type].op;

  if (vfptr == 0)
  {
    fprintf (stderr, "Entity type: %s\n", etd (e));
    rerror ("increment (++) operation not supported");
  }

  rval = (*vfptr) (ent_data (e));

  ent_data (e) = rval;
  ent_type (e) = rtype;

  return (e);
}

Ent *
class_decrement (Ent * e)
{
  int type=UNDEF, rtype;
  void *rval=0;
  void *(*vfptr) ();

  type = ent_type (e);
  rtype = decrement_method[type].type;
  vfptr = (VFPTR) decrement_method[type].op;

  if (vfptr == 0)
  {
    fprintf (stderr, "Entity type: %s\n", etd (e));
    rerror ("decrement (--) operation not supported");
  }

  rval = (*vfptr) (ent_data (e));

  ent_data (e) = rval;
  ent_type (e) = rtype;

  return (e);
}

/* **************************************************************
 * Get an objects "sublist". This is used by list_assign() to
 * get/check the validity of an object sub-list.
 *
 * This function will soon be obsolte (3/30/97).
 * ************************************************************** */

int
class_attribute (Ent * e, char *name)
{
  int type=UNDEF;
  int rval;
  int (*vfptr) ();

  if (e)
    type = ent_type (e);

  vfptr = (IFPTR) attribute_method[type].op;
  if (vfptr == 0)
  {
    fprintf (stderr, "Entity type: %s\n", etd (e));
    rerror ("attribute operation not supported");
  }

  rval = (int) (*vfptr) (ent_data (e), name);
  return (rval);
}
