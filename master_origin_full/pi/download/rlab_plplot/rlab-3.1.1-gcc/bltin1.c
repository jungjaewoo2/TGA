/* bltin1.c */

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
#include "listnode.h"
#include "ent.h"
#include "mem.h"
#include "bltin.h"
#include "util.h"
#include "ent.h"
#include "class.h"
#include "rfileio.h"
#include "mathl.h"
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <sys/time.h>

#ifdef __riscos
#include "riscos_ut.h"
#endif

/*
 * Header files from the available classes...
 */

#include "btree.h"
#include "btreef1.h"
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
#include "msc.h"
#include "mscf1.h"
#include "mde.h"
#include "sort.h"

#include <sys/types.h>
#include <unistd.h>

#include "rlab_solver_parameters_names.h"
#include "rlab_macros.h"
#include "rlab_macros_code.h"

#define THIS_FILE "bltin1.c"

extern FILE *RLAB_STDERR_DS;

static OpDef shift_method[NCL];
static OpDef flip_method[NCL];
static OpDef any_method[NCL];
static OpDef all_method[NCL];
static OpDef class_method[NCL];
static OpDef members_method[NCL];
static OpDef size_method[NCL];
static OpDef length_method[NCL];
static OpDef type_method[NCL];
static OpDef reshape_method[NCL];
static OpDef isinf_method[NCL];
static OpDef isnan_method[NCL];
static OpDef sum_method[NCL];
static OpDef cumsum_method[NCL];
static OpDef find_method[NCL];
static OpDef findvec_method[NCL][NCL];
static OpDef sort_method[NCL];
static OpDef vectorset_method[NCL];
static OpDef vectorunion_method[NCL];
static OpDef vectorintersect_method[NCL];
static OpDef vectorcomplement_method[NCL];
static OpDef merge_method[NCL];
static OpDef compact_method[NCL];
static OpDef strtol_method[NCL];
static OpDef sign_method[NCL];
static OpDef prod_method[NCL];
static OpDef cumprod_method[NCL];
static OpDef frexp_method[NCL];
static OpDef ldexp_method[NCL][NCL];
#ifdef HAVE_LOGB
static OpDef logb_method[NCL];
#endif /* HAVE_LOGB */

void
class_bltin1_init (void)
{
 /*
  * flipuplr ()
  */

  flip_method[MATRIX_DENSE_REAL].type = MATRIX_DENSE_REAL;
  flip_method[MATRIX_DENSE_REAL].op = (void *) &md_numeric_Flip;

  flip_method[MATRIX_DENSE_COMPLEX].type = MATRIX_DENSE_COMPLEX;
  flip_method[MATRIX_DENSE_COMPLEX].op = (void *) &md_numeric_Flip;

  flip_method[MATRIX_DENSE_STRING].type = MATRIX_DENSE_STRING;
  flip_method[MATRIX_DENSE_STRING].op = (void *) &mds_Flip;

 /*
  * shiftuplr ()
  */

  shift_method[MATRIX_DENSE_REAL].type = MATRIX_DENSE_REAL;
  shift_method[MATRIX_DENSE_REAL].op = (void *) &md_numeric_Shift;

  shift_method[MATRIX_DENSE_COMPLEX].type = MATRIX_DENSE_COMPLEX;
  shift_method[MATRIX_DENSE_COMPLEX].op = (void *) &md_numeric_Shift;

  shift_method[MATRIX_DENSE_STRING].type = MATRIX_DENSE_STRING;
  shift_method[MATRIX_DENSE_STRING].op = (void *) &mds_Shift;

  /*
   * any ()
   */

  any_method[MATRIX_DENSE_REAL].type = MATRIX_DENSE_REAL;
  any_method[MATRIX_DENSE_REAL].op = (void *) mdr_Any;

  any_method[MATRIX_DENSE_COMPLEX].type = MATRIX_DENSE_REAL;
  any_method[MATRIX_DENSE_COMPLEX].op = (void *) mdc_Any;

  any_method[MATRIX_SPARSE_REAL].type = MATRIX_DENSE_REAL;
  any_method[MATRIX_SPARSE_REAL].op = (void *) msr_Any;

  any_method[MATRIX_SPARSE_COMPLEX].type = MATRIX_DENSE_REAL;
  any_method[MATRIX_SPARSE_COMPLEX].op = (void *) msc_Any;

  /*
   * all ()
   */

  all_method[MATRIX_DENSE_REAL].type = MATRIX_DENSE_REAL;
  all_method[MATRIX_DENSE_REAL].op = (void *) mdr_All;

  all_method[MATRIX_DENSE_COMPLEX].type = MATRIX_DENSE_REAL;
  all_method[MATRIX_DENSE_COMPLEX].op = (void *) mdc_All;

  all_method[MATRIX_SPARSE_REAL].type = MATRIX_DENSE_REAL;
  all_method[MATRIX_SPARSE_REAL].op = (void *) msr_All;

  all_method[MATRIX_SPARSE_COMPLEX].type = MATRIX_DENSE_REAL;
  all_method[MATRIX_SPARSE_COMPLEX].op = (void *) msc_All;

  /*
   * class ()
   */

  class_method[UNDEF].type = -1;
  class_method[UNDEF].op = (void *) undef_Class;

  class_method[DOUBLE].type = -1;
  class_method[DOUBLE].op = (void *) double_Class;

  class_method[BTREE].type = -1;
  class_method[BTREE].op = (void *) btree_Class;

  class_method[U_FUNCTION].type = -1;
  class_method[U_FUNCTION].op = (void *) function_Class;

  class_method[BLTIN].type = -1;
  class_method[BLTIN].op = (void *) bltin_Class;

  class_method[MATRIX_DENSE_REAL].type = -1;
  class_method[MATRIX_DENSE_REAL].op = (void *) mdr_Class;

  class_method[MATRIX_DENSE_COMPLEX].type = -1;
  class_method[MATRIX_DENSE_COMPLEX].op = (void *) mdc_Class;

  class_method[MATRIX_DENSE_STRING].type = -1;
  class_method[MATRIX_DENSE_STRING].op = (void *) mds_Class;

  class_method[MATRIX_SPARSE_REAL].type = -1;
  class_method[MATRIX_SPARSE_REAL].op = (void *) msr_Class;

  class_method[MATRIX_SPARSE_COMPLEX].type = -1;
  class_method[MATRIX_SPARSE_COMPLEX].op = (void *) msc_Class;

  class_method[MATRIX_DENSE_ENTITY].type = -1;
  class_method[MATRIX_DENSE_ENTITY].op = (void *) mde_Class;
  /*
   * members ()
   * By definition, members returns a dense-string vector.
   */

  members_method[BTREE].type = -1;
  members_method[BTREE].op = (void *) btree_Members;

  members_method[MATRIX_DENSE_ENTITY].type = -1;
  members_method[MATRIX_DENSE_ENTITY].op = (void *) mde_Members;

  members_method[MATRIX_DENSE_REAL].type = -1;
  members_method[MATRIX_DENSE_REAL].op = (void *) mdr_Members;

  members_method[MATRIX_SPARSE_REAL].type = -1;
  members_method[MATRIX_SPARSE_REAL].op = (void *) msr_Members;

  members_method[MATRIX_DENSE_COMPLEX].type = -1;
  members_method[MATRIX_DENSE_COMPLEX].op = (void *) mdc_Members;

  members_method[MATRIX_SPARSE_COMPLEX].type = -1;
  members_method[MATRIX_SPARSE_COMPLEX].op = (void *) msc_Members;

  members_method[MATRIX_DENSE_STRING].type = -1;
  members_method[MATRIX_DENSE_STRING].op = (void *) mds_Members;

  members_method[U_FUNCTION].type = -1;
  members_method[U_FUNCTION].op = (void *) function_Members;

  members_method[BLTIN].type = -1;
  members_method[BLTIN].op = (void *) bltin_Members;

  /*
   * size ()
   * By definition, method returns a double pointer.
   */

  size_method[BTREE].type = MATRIX_DENSE_REAL;
  size_method[BTREE].op = (void *) btree_Size_BF;

  size_method[MATRIX_DENSE_REAL].type = MATRIX_DENSE_REAL;
  size_method[MATRIX_DENSE_REAL].op = (void *) mdr_Size_BF;

  size_method[MATRIX_DENSE_COMPLEX].type = MATRIX_DENSE_REAL;
  size_method[MATRIX_DENSE_COMPLEX].op = (void *) mdc_Size_BF;

  size_method[MATRIX_DENSE_STRING].type = MATRIX_DENSE_REAL;
  size_method[MATRIX_DENSE_STRING].op = (void *) mds_Size_BF;

  size_method[MATRIX_SPARSE_REAL].type = MATRIX_DENSE_REAL;
  size_method[MATRIX_SPARSE_REAL].op = (void *) msr_Size_BF;

  size_method[MATRIX_SPARSE_COMPLEX].type = MATRIX_DENSE_REAL;
  size_method[MATRIX_SPARSE_COMPLEX].op = (void *) msc_Size_BF;

  size_method[MATRIX_DENSE_ENTITY].type = MATRIX_DENSE_REAL;
  size_method[MATRIX_DENSE_ENTITY].op = (void *) mde_Size_BF;

  /*
   * length ()
   * By definition, members returns a dense-real vector.
   */

  length_method[BTREE].type = MATRIX_DENSE_REAL;
  length_method[BTREE].op = (void *) btree_Size_BF;

  length_method[MATRIX_DENSE_REAL].type = MATRIX_DENSE_REAL;
  length_method[MATRIX_DENSE_REAL].op = (void *) mdr_Length_BF;

  length_method[MATRIX_DENSE_COMPLEX].type = MATRIX_DENSE_REAL;
  length_method[MATRIX_DENSE_COMPLEX].op = (void *) mdc_Length_BF;

  length_method[MATRIX_DENSE_STRING].type = MATRIX_DENSE_REAL;
  length_method[MATRIX_DENSE_STRING].op = (void *) mds_Length_BF;

  length_method[MATRIX_SPARSE_REAL].type = MATRIX_DENSE_REAL;
  length_method[MATRIX_SPARSE_REAL].op = (void *) msr_Length_BF;

  length_method[MATRIX_SPARSE_COMPLEX].type = MATRIX_DENSE_REAL;
  length_method[MATRIX_SPARSE_COMPLEX].op = (void *) msc_Length_BF;

  length_method[MATRIX_DENSE_ENTITY].type = MATRIX_DENSE_REAL;
  length_method[MATRIX_DENSE_ENTITY].op = (void *) mde_Length_BF;

  /*
   * type ()
   * By definition, members returns a dense-string vector.
   */

  type_method[UNDEF].type = MATRIX_DENSE_STRING;
  type_method[UNDEF].op = (void *) undef_Type_BF;

  type_method[DOUBLE].type = MATRIX_DENSE_STRING;
  type_method[DOUBLE].op = (void *) double_Type_BF;

  type_method[BTREE].type = MATRIX_DENSE_STRING;
  type_method[BTREE].op = (void *) btree_Type_BF;

  type_method[U_FUNCTION].type = MATRIX_DENSE_STRING;
  type_method[U_FUNCTION].op = (void *) function_Type_BF;

  type_method[BLTIN].type = MATRIX_DENSE_STRING;
  type_method[BLTIN].op = (void *) bltin_Type_BF;

  type_method[MATRIX_DENSE_REAL].type = MATRIX_DENSE_STRING;
  type_method[MATRIX_DENSE_REAL].op = (void *) mdr_Type_BF;

  type_method[MATRIX_DENSE_COMPLEX].type = MATRIX_DENSE_STRING;
  type_method[MATRIX_DENSE_COMPLEX].op = (void *) mdc_Type_BF;

  type_method[MATRIX_DENSE_STRING].type = MATRIX_DENSE_STRING;
  type_method[MATRIX_DENSE_STRING].op = (void *) mds_Type_BF;

  type_method[MATRIX_SPARSE_REAL].type = MATRIX_DENSE_STRING;
  type_method[MATRIX_SPARSE_REAL].op = (void *) msr_Type_BF;

  type_method[MATRIX_SPARSE_COMPLEX].type = MATRIX_DENSE_STRING;
  type_method[MATRIX_SPARSE_COMPLEX].op = (void *) msc_Type_BF;

  /*
   * reshape ()
   */

  reshape_method[MATRIX_DENSE_REAL].type = MATRIX_DENSE_REAL;
  reshape_method[MATRIX_DENSE_REAL].op = (void *) mdr_Reshape;

  reshape_method[MATRIX_DENSE_COMPLEX].type = MATRIX_DENSE_COMPLEX;
  reshape_method[MATRIX_DENSE_COMPLEX].op = (void *) mdc_Reshape;

  reshape_method[MATRIX_DENSE_STRING].type = MATRIX_DENSE_STRING;
  reshape_method[MATRIX_DENSE_STRING].op = (void *) mds_Reshape;

  /*
   * isinf ()
   */

  isinf_method[MATRIX_DENSE_REAL].type = MATRIX_DENSE_REAL;
  isinf_method[MATRIX_DENSE_REAL].op = (void *) mdr_IsInf;

  isinf_method[MATRIX_DENSE_COMPLEX].type = MATRIX_DENSE_REAL;
  isinf_method[MATRIX_DENSE_COMPLEX].op = (void *) mdc_IsInf;

  /*
   * isnan ()
   */

  isnan_method[MATRIX_DENSE_REAL].type = MATRIX_DENSE_REAL;
  isnan_method[MATRIX_DENSE_REAL].op = (void *) mdr_IsNan;

  isnan_method[MATRIX_DENSE_COMPLEX].type = MATRIX_DENSE_REAL;
  isnan_method[MATRIX_DENSE_COMPLEX].op = (void *) mdc_IsNan;

  /*
   * sum ()
   */

  sum_method[MATRIX_DENSE_REAL].type = MATRIX_DENSE_REAL;
  sum_method[MATRIX_DENSE_REAL].op = (void *) mdr_Sum_BF;

  sum_method[MATRIX_DENSE_COMPLEX].type = MATRIX_DENSE_COMPLEX;
  sum_method[MATRIX_DENSE_COMPLEX].op = (void *) mdc_Sum_BF;

  sum_method[MATRIX_DENSE_STRING].type = MATRIX_DENSE_STRING;
  sum_method[MATRIX_DENSE_STRING].op = (void *) mds_Sum;

  sum_method[MATRIX_SPARSE_REAL].type = MATRIX_DENSE_REAL;
  sum_method[MATRIX_SPARSE_REAL].op = (void *) msr_Sum_BF;

  sum_method[MATRIX_SPARSE_COMPLEX].type = MATRIX_DENSE_COMPLEX;
  sum_method[MATRIX_SPARSE_COMPLEX].op = (void *) msc_Sum_BF;

  /*
   * cumsum ()
   */

  cumsum_method[MATRIX_DENSE_REAL].type = MATRIX_DENSE_REAL;
  cumsum_method[MATRIX_DENSE_REAL].op = (void *) mdr_CumSum_BF;

  cumsum_method[MATRIX_DENSE_COMPLEX].type = MATRIX_DENSE_COMPLEX;
  cumsum_method[MATRIX_DENSE_COMPLEX].op = (void *) mdc_CumSum_BF;

  /*
   * find ()
   */

  find_method[MATRIX_DENSE_REAL].type = MATRIX_DENSE_REAL;
  find_method[MATRIX_DENSE_REAL].op = (void *) mdr_Find_BF;

  find_method[MATRIX_DENSE_COMPLEX].type = MATRIX_DENSE_REAL;
  find_method[MATRIX_DENSE_COMPLEX].op = (void *) mdc_Find_BF;

  find_method[MATRIX_DENSE_STRING].type = MATRIX_DENSE_REAL;
  find_method[MATRIX_DENSE_STRING].op = (void *) mds_Find_BF;

  find_method[MATRIX_SPARSE_REAL].type = MATRIX_DENSE_REAL;
  find_method[MATRIX_SPARSE_REAL].op = (void *) msr_Find_BF;

  find_method[MATRIX_SPARSE_COMPLEX].type = MATRIX_DENSE_REAL;
  find_method[MATRIX_SPARSE_COMPLEX].op = (void *) msc_Find_BF;

  /*
  * findvec ()
  */

  findvec_method[MATRIX_DENSE_REAL][MATRIX_DENSE_REAL].type = MATRIX_DENSE_REAL;
  findvec_method[MATRIX_DENSE_REAL][MATRIX_DENSE_REAL].op = (void *) mdr_FindVector;

  /*
   * sort ()
   */

  sort_method[MATRIX_DENSE_REAL].type = BTREE;
  sort_method[MATRIX_DENSE_REAL].op = (void *) mdr_Sort_BF;

  sort_method[MATRIX_DENSE_COMPLEX].type = BTREE;
  sort_method[MATRIX_DENSE_COMPLEX].op = (void *) mdc_Sort_BF;

  sort_method[MATRIX_DENSE_STRING].type = BTREE;
  sort_method[MATRIX_DENSE_STRING].op = (void *) mds_Sort_BF;

  /*
   * vectorset
   */
  vectorset_method[MATRIX_DENSE_REAL].type  = MATRIX_DENSE_REAL;
  vectorset_method[MATRIX_DENSE_REAL].op    = (void *) mdr_VectorSet;
  vectorset_method[MATRIX_DENSE_STRING].type  = MATRIX_DENSE_STRING;
  vectorset_method[MATRIX_DENSE_STRING].op    = (void *) mds_VectorSet;

  /*
   * vectorunion
   */
  vectorunion_method[MATRIX_DENSE_REAL].type  = MATRIX_DENSE_REAL;
  vectorunion_method[MATRIX_DENSE_REAL].op    = (void *) mdr_VectorUnion;
  vectorunion_method[MATRIX_DENSE_STRING].type  = MATRIX_DENSE_STRING;
  vectorunion_method[MATRIX_DENSE_STRING].op    = (void *) mds_VectorUnion;

  /*
   * vectorintersect
   */
  vectorintersect_method[MATRIX_DENSE_REAL].type  = MATRIX_DENSE_REAL;
  vectorintersect_method[MATRIX_DENSE_REAL].op    = (void *) mdr_VectorIntersect;
  vectorintersect_method[MATRIX_DENSE_STRING].type  = MATRIX_DENSE_STRING;
  vectorintersect_method[MATRIX_DENSE_STRING].op    = (void *) mds_VectorIntersect;

  /*
   * vectorcomplement
   */
  vectorcomplement_method[MATRIX_DENSE_REAL].type  = MATRIX_DENSE_REAL;
  vectorcomplement_method[MATRIX_DENSE_REAL].op    = (void *) mdr_VectorComplement;
  vectorcomplement_method[MATRIX_DENSE_STRING].type  = MATRIX_DENSE_STRING;
  vectorcomplement_method[MATRIX_DENSE_STRING].op    = (void *) mds_VectorComplement;

  /*
   * merge
   */
  merge_method[MATRIX_DENSE_REAL].type  = MATRIX_DENSE_REAL;
  merge_method[MATRIX_DENSE_REAL].op    = (void *) mdr_Merge;

  /*
   * compact
   */
  compact_method[MATRIX_DENSE_REAL].type  = MATRIX_DENSE_REAL;
  compact_method[MATRIX_DENSE_REAL].op    = (void *) mdr_Compact;

  /*
   * strtol ()
   */

  strtol_method[MATRIX_DENSE_STRING].type = MATRIX_DENSE_REAL;
  strtol_method[MATRIX_DENSE_STRING].op = (void *) mds_Strtol_BF;

  /*
   * sign ()
   */

  sign_method[MATRIX_DENSE_REAL].type = MATRIX_DENSE_REAL;
  sign_method[MATRIX_DENSE_REAL].op = (void *) mdr_Sign_BF;

  sign_method[MATRIX_DENSE_COMPLEX].type = MATRIX_DENSE_COMPLEX;
  sign_method[MATRIX_DENSE_COMPLEX].op = (void *) mdc_Sign_BF;

  /*
   * prod ()
   */

  prod_method[MATRIX_DENSE_REAL].type = MATRIX_DENSE_REAL;
  prod_method[MATRIX_DENSE_REAL].op = (void *) mdr_Prod_BF;

  prod_method[MATRIX_DENSE_COMPLEX].type = MATRIX_DENSE_COMPLEX;
  prod_method[MATRIX_DENSE_COMPLEX].op = (void *) mdc_Prod_BF;

  /*
   * cumprod ()
   */

  cumprod_method[MATRIX_DENSE_REAL].type = MATRIX_DENSE_REAL;
  cumprod_method[MATRIX_DENSE_REAL].op = (void *) mdr_CumProd_BF;

  cumprod_method[MATRIX_DENSE_COMPLEX].type = MATRIX_DENSE_COMPLEX;
  cumprod_method[MATRIX_DENSE_COMPLEX].op = (void *) mdc_CumProd_BF;

  /*
   * frexp ()
   */

  frexp_method[MATRIX_DENSE_REAL].type = BTREE;
  frexp_method[MATRIX_DENSE_REAL].op = (void *) mdr_Frexp_BF;

  /*
   * ldexp ()
   */

  ldexp_method[MATRIX_DENSE_REAL][MATRIX_DENSE_REAL].type = MATRIX_DENSE_REAL;
  ldexp_method[MATRIX_DENSE_REAL][MATRIX_DENSE_REAL].op = (void *) mdr_Ldexp_BF;

#ifdef HAVE_LOGB
  /*
   * logb ()
   */

  logb_method[MATRIX_DENSE_REAL].type = MATRIX_DENSE_REAL;
  logb_method[MATRIX_DENSE_REAL].op = (void *) mdr_Logb_BF;
#endif /* HAVE_LOGB */

}

//
// change the size of a dense matrix in-place, that is, without creating
// a new entity:
//  x = rand(n1,n2);
//  x = reshape(x,1,n1*n2);   // makes a copy of 'x' which is then resized
//  resize(x,1,n1*n2);        // directly resizes, without making a copy first
//
#undef  THIS_SOLVER
#define THIS_SOLVER "resize"
Ent * ent_reshape_resize (int nargs, Datum args[])
{
  // call parameters:
  //  e1   - MDR
  Ent *e1=0, *e2=0, *e3=0;
  ListNode *var;

  int nr=0, nc=0, n=0;

  if (nargs!=3 && nargs!=2 )
    rerror(THIS_SOLVER ": two or three arguments required");

  //
  // get the matrix if it is used by more then one variable
  // make a copy of it
  //
  RLABCODE_GETARG_VAR_OR_ENTITY(THIS_SOLVER,1,e1,var);

//   var = (ListNode *) (args[0].u.ptr);
//   if (args[0].type != VAR)
//     rerror("\nresize: first argument has to be variable and not expression!");
//   e1 = var_ent (var);
//   if (e1->refc > 1)
//   {
//     ent_DecRef (e1);
//     e1 = ent_Duplicate (e1);
//     listNode_AttachEnt (var, e1);
//   }

  // it also suffices to provide one of the two sizes
  e2 = bltin_get_ent (args[1]);
  if (ent_type(e2) == MATRIX_DENSE_REAL)
    nr = (int) class_double(e2);
  if (nargs == 3)
  {
    e3 = bltin_get_ent (args[2]);
    if (ent_type(e3) == MATRIX_DENSE_REAL)
      nc = (int) class_double(e3);
  }

  // making things tard-proof
  if (nr || nc)
  {
    if (ent_type (e1) == MATRIX_DENSE_REAL)
    {
      MDR *x1 = ent_data (e1);
      n = x1->nrow * x1->ncol;
      if (n)
      {
        if(nr && !(n%nr) && !nc)
          nc = n/nr;
        else if (nc && !(n%nc) && !nr)
          nr = n/nc;

        if (nr && nc && n==nr*nc)
        {
          x1->nrow = nr;
          x1->ncol = nc;
        }
      }
    }
    else if (ent_type (e1) == MATRIX_DENSE_COMPLEX)
    {
      MDC * x1 = ent_data (e1);
      n = x1->nrow * x1->ncol;
      if (n)
      {
        if(nr && !(n%nr) && !nc)
          nc = n/nr;
        else if (nc && !(n%nc) && !nr)
          nr = n/nc;

        if (nr && nc && n==nr*nc)
        {
          x1->nrow = nr;
          x1->ncol = nc;
        }
      }
    }
    else if (ent_type (e1) == MATRIX_DENSE_STRING)
    {
      MDS * x1 = ent_data (e1);
      n = x1->nrow * x1->ncol;
      if (n)
      {
        if(nr && !(n%nr) && !nc)
          nc = n/nr;
        else if (nc && !(n%nc) && !nr)
          nr = n/nc;

        if (nr && nc && n==nr*nc)
        {
          x1->nrow = nr;
          x1->ncol = nc;
        }
      }
    }
  }

  ent_Clean (e1);
  ent_Clean (e2);
  ent_Clean (e3);

  return ent_Create_Rlab_Success();
}


#undef  THIS_SOLVER
#define THIS_SOLVER "chomp"
Ent * ent_Chomp (int nargs, Datum args[])
{
  // call parameters:
  Ent *e1=0;
  ListNode *var=0;
  MDS *s=0;

  int i, rval=RLAB_STATUS_FAILURE;

  if (nargs!=1)
    rerror(THIS_SOLVER ": " RLAB_ERROR_ONE_ARG_REQUIRED  "\n");

  //
  // get the matrix if it is used by more then one variable
  // make a copy of it
  //
  RLABCODE_GETARG_VAR_OR_ENTITY(THIS_SOLVER,1,e1,var);

  if (ent_type(e1) != MATRIX_DENSE_STRING)
    goto _exit_chomp;

  s = ent_data (e1);
  if (SIZE(s)<1)
    goto _exit_chomp;

  rval = RLAB_STATUS_SUCCESS;

  for(i=0; i<SIZE(s); i++)
  {
    chomp2_string(MdsV0(s,i));
  }

_exit_chomp:

  ent_Clean (e1);
  return ent_Create_Rlab_Double(rval);
}

#undef  THIS_SOLVER
#define THIS_SOLVER "fliplr"
#undef  METHOD
#define METHOD flip_method
Ent * ent_FlipLR(int nargs, Datum args[])
{
  Ent *e1=0;
  ListNode *var=0;
  void *(*vfptr) ();

  int type1;

  if (nargs!=1)
    rerror(THIS_SOLVER ": " RLAB_ERROR_ONE_ARG_REQUIRED "!\n");

  //
  // get the matrix if it is used by more then one variable
  // make a copy of it
  //
  RLABCODE_GETARG_VAR_OR_ENTITY(THIS_SOLVER,1,e1,var);

  type1 = ent_type (e1);
  vfptr = (VFPTR) METHOD[type1].op;

  if (vfptr == 0)
  {
    fprintf (stderr, "Entity type: %s: Operation not supported\n", etd (e1));
    return ent_Create_Rlab_Failure();
  }
  (*vfptr) (ent_data (e1),0,1);

  ent_Clean (e1);

  return ent_Create_Rlab_Success();
}



#undef  THIS_SOLVER
#define THIS_SOLVER "flipud"
#undef  METHOD
#define METHOD flip_method
Ent * ent_FlipUD(int nargs, Datum args[])
{
  Ent *e1=0;
  ListNode *var=0;
  void *(*vfptr) ();

  int type1;

  if (nargs!=1)
    rerror(THIS_SOLVER ": " RLAB_ERROR_ONE_ARG_REQUIRED "!\n");

  //
  // get the matrix if it is used by more then one variable
  // make a copy of it
  //
  RLABCODE_GETARG_VAR_OR_ENTITY(THIS_SOLVER,1,e1,var);

  type1 = ent_type (e1);
  vfptr = (VFPTR) METHOD[type1].op;

  if (vfptr == 0)
  {
    fprintf (stderr, "Entity type: %s: Operation not supported\n", etd (e1));
    return ent_Create_Rlab_Failure();
  }
  (*vfptr) (ent_data (e1),1,0);

  ent_Clean (e1);
  return ent_Create_Rlab_Success();
}

#undef  THIS_SOLVER
#define THIS_SOLVER "shiftu"
#undef  METHOD
#define METHOD shift_method
Ent * ent_ShiftU(int nargs, Datum args[])
{
  // call parameters:
  //  e1   - MDR
  Ent *e1=0, *e2=0;
  ListNode *var=0;
  void *(*vfptr) ();

  int type1, u=1;

  if ((nargs!=1)&&(nargs!=2))
    rerror(THIS_SOLVER ": " RLAB_ERROR_ONE_OR_TWO_ARG_REQUIRED "!\n");

  //
  // get the matrix if it is used by more then one variable
  // make a copy of it
  //
  RLABCODE_GETARG_VAR_OR_ENTITY(THIS_SOLVER,1,e1,var);

  type1 = ent_type (e1);
  vfptr = (VFPTR) METHOD[type1].op;

  if (nargs > 1)
  {
    e2 = bltin_get_ent (args[1]);
    if (ent_type(e2)==MATRIX_DENSE_REAL)
      u = class_int(e2);
  }

  if (vfptr == 0)
  {
    fprintf (stderr, "Entity type: %s: Operation not supported\n", etd (e1));
    return ent_Create_Rlab_Failure();
  }
  (*vfptr) (ent_data (e1),u,0);

  ent_Clean (e1);
  ent_Clean (e2);

  return ent_Create_Rlab_Success();
}

#undef  THIS_SOLVER
#define THIS_SOLVER "shiftd"
#undef  METHOD
#define METHOD shift_method
Ent * ent_ShiftD(int nargs, Datum args[])
{
  // call parameters:
  //  e1   - MDR
  Ent *e1=0, *e2=0;
  ListNode *var=0;
  void *(*vfptr) ();

  int type1, u=-1;

  if ((nargs!=1)&&(nargs!=2))
    rerror(THIS_SOLVER ": " RLAB_ERROR_ONE_OR_TWO_ARG_REQUIRED "!\n");

  //
  // get the matrix if it is used by more then one variable
  // make a copy of it
  //
  RLABCODE_GETARG_VAR_OR_ENTITY(THIS_SOLVER,1,e1,var);

  type1 = ent_type (e1);
  vfptr = (VFPTR) METHOD[type1].op;

  if (nargs > 1)
  {
    e2 = bltin_get_ent (args[1]);
    if (ent_type(e2)==MATRIX_DENSE_REAL)
      u = -class_int(e2);
  }

  if (vfptr == 0)
  {
    fprintf (stderr, "Entity type: %s: Operation not supported\n", etd (e1));
    return ent_Create_Rlab_Failure();
  }
  (*vfptr) (ent_data (e1),u,0);
  ent_Clean (e1);
  ent_Clean (e2);
  return ent_Create_Rlab_Success();
}

#undef  THIS_SOLVER
#define THIS_SOLVER "shiftl"
#undef  METHOD
#define METHOD shift_method
Ent * ent_ShiftL(int nargs, Datum args[])
{
  // call parameters:
  //  e1   - MDR
  Ent *e1=0, *e2=0;
  ListNode *var=0;
  void *(*vfptr) ();

  int type1, u=1;

  if ((nargs!=1)&&(nargs!=2))
    rerror(THIS_SOLVER ": " RLAB_ERROR_ONE_OR_TWO_ARG_REQUIRED "!\n");

  //
  // get the matrix if it is used by more then one variable
  // make a copy of it
  //
  RLABCODE_GETARG_VAR_OR_ENTITY(THIS_SOLVER,1,e1,var);

  type1 = ent_type (e1);
  vfptr = (VFPTR) METHOD[type1].op;

  if (nargs > 1)
  {
    e2 = bltin_get_ent (args[1]);
    if (ent_type(e2)==MATRIX_DENSE_REAL)
      u = class_int(e2);
  }

  if (vfptr == 0)
  {
    fprintf (stderr, "Entity type: %s: Operation not supported\n", etd (e1));
    return ent_Create_Rlab_Failure();
  }
  (*vfptr) (ent_data (e1),0,u);
  ent_Clean (e1);
  ent_Clean (e2);
  return ent_Create_Rlab_Success();
}

#undef  THIS_SOLVER
#define THIS_SOLVER "shiftr"
#undef  METHOD
#define METHOD shift_method
Ent * ent_ShiftR(int nargs, Datum args[])
{
  // call parameters:
  //  e1   - MDR
  Ent *e1=0, *e2=0;
  ListNode *var=0;
  void *(*vfptr) ();

  int type1, u=-1;

  if ((nargs!=1)&&(nargs!=2))
    rerror(THIS_SOLVER ": " RLAB_ERROR_ONE_OR_TWO_ARG_REQUIRED "!\n");

  //
  // get the matrix if it is used by more then one variable
  // make a copy of it
  //
  RLABCODE_GETARG_VAR_OR_ENTITY(THIS_SOLVER,1,e1,var);

  type1 = ent_type (e1);
  vfptr = (VFPTR) METHOD[type1].op;

  if (nargs > 1)
  {
    e2 = bltin_get_ent (args[1]);
    if (ent_type(e2)==MATRIX_DENSE_REAL)
      u = -class_int(e2);
  }

  if (vfptr == 0)
  {
    fprintf (stderr, "Entity type: %s: Operation not supported\n", etd (e1));
    return ent_Create_Rlab_Failure();
  }
  (*vfptr) (ent_data (e1),0,u);
  ent_Clean (e1);
  ent_Clean (e2);
  return ent_Create_Rlab_Success();
}

#undef  THIS_SOLVER
#define THIS_SOLVER "rot90"
#undef  METHOD
#define METHOD shift_method
Ent * ent_Rot90(int nargs, Datum args[])
{
  // call parameters:
  //  e1   - MDR
  Ent *e1=0, *e2=0;
  MD *m=0;
  ListNode *var=0;
  int n=1;


  if ((nargs!=1)&&(nargs!=2))
    rerror(THIS_SOLVER ": " RLAB_ERROR_ONE_OR_TWO_ARG_REQUIRED "!\n");

  //
  // get the matrix if it is used by more then one variable
  // make a copy of it
  //
  RLABCODE_GETARG_VAR_OR_ENTITY(THIS_SOLVER,1,e1,var);
  if ((ent_type(e1) != MATRIX_DENSE_REAL) &&
       (ent_type(e1) != MATRIX_DENSE_COMPLEX) &&
       (ent_type(e1) != MATRIX_DENSE_STRING) )
  {
    goto exit;
  }
  m = ent_data(e1);

  //
  // second argument is how many 90deg rotations of matrix to perform
  //
  if (nargs > 1)
  {
    e2 = bltin_get_ent (args[1]);
    if (ent_type(e2)==MATRIX_DENSE_REAL)
      n = class_int(e2);
  }
  while (n>4)
  {
    n -= 4;
  }
  while (n<0)
  {
    n += 4;
  }

  switch (n)
  {
    case 1:
      // 90deg = (m^t)^ud
      md_transpose_insitu((unsigned char *) MDPTR(m), MNR(m), MNC(m), rsizeof(m->type));
      SWAP(MNR(m),MNC(m),int);
      md_numeric_Flip ((MD *) m, 1, 0);
      break;

    case 2:
      // 180deg = (m^lr)^ud
      md_numeric_Flip ((MD *) m, 0, 1);
      md_numeric_Flip ((MD *) m, 1, 0);
      break;

    case 3:
      // 270deg = -90deg = (m^t)^lr
      md_transpose_insitu((unsigned char *) MDPTR(m), MNR(m), MNC(m), rsizeof(m->type));
      SWAP(MNR(m),MNC(m),int);
      md_numeric_Flip ((MD *) m, 0, 1);
      break;

    default:
      // 4*90=360, nothing to do
      break;
  }

exit:
  ent_Clean (e1);
  ent_Clean (e2);
  return ent_Create_Rlab_Success();
}

Ent * Any (int nargs, Datum args[])
{
  int type, rtype;
  void *rval;
  void *(*vfptr) ();
  Ent *e=0;

  if (nargs != 1)
  {
    rerror ("any: 1 argument allowed");
  }

  e = bltin_get_ent (args[0]);
  type = ent_type (e);

  rtype = any_method[type].type;
  vfptr = (VFPTR) any_method[type].op;

  if (vfptr == 0)
  {
    fprintf (stderr, "Entity type: %s\n", etd (e));
    rerror ("any() operation not supported");
  }

  rval = (*vfptr) (ent_data (e));

  ent_Clean (e);

  return ent_Assign_Rlab_Rtype(rval, rtype);
}

Ent * All (int nargs, Datum args[])
{
  int type, rtype;
  void *rval;
  void *(*vfptr) ();
  Ent *e=0;

  if (nargs != 1)
  {
    rerror ("all: 1 argument allowed");
  }

  e = bltin_get_ent (args[0]);
  type = ent_type (e);

  rtype = all_method[type].type;
  vfptr = (VFPTR) all_method[type].op;

  if (vfptr == 0)
  {
    fprintf (stderr, "Entity type: %s\n", etd (e));
    rerror ("all() operation not supported");
  }

  rval = (*vfptr) (ent_data (e));

  ent_Clean (e);

  return ent_Assign_Rlab_Rtype(rval, rtype);
}

/* **************************************************************
 * Get the class of an object.
 * ************************************************************** */
#undef THIS_SOLVER
#define THIS_SOLVER "class"
Ent * Group (int nargs, Datum args[])
{
  int type;
  char *rval=0;
  void *(*vfptr) ();
  Ent *e=0;

  if (nargs != 1)
  {
    rerror (THIS_SOLVER ": " RLAB_ERROR_ARG1_ONLY);
  }

  e = bltin_get_ent (args[0]);
  type = ent_type (e);

  vfptr = (VFPTR) class_method[type].op;

  if (vfptr == 0)
  {
    fprintf (stderr, "Entity type: %s\n", etd (e));
    rerror (THIS_SOLVER "() " RLAB_ERROR_OPERATION_FAILED);
  }

  rval = (*vfptr) (ent_data (e));

  ent_Clean (e);

  return ent_Create_Rlab_String(rval);
}

/* **************************************************************
 * Get the members of a list or an object. Each method _must_
 * return an array of string pointers. Members will re-format
 * the thing into a Matrix_Dense_String.
 * ************************************************************** */

Ent *
Members (int nargs, Datum args[])
{
  int n, type;
//   int rtype;
  char **rval, **rval2;
  void *(*vfptr) ();
  Ent *e=0;
  ListNode *listent, *var;
  MDS *smat, *smat2, *smat3;

  if (nargs != 1)
  {
    rerror ("members: 1 argument allowed");
  }

  /* We must be sure to get the _variable_. If we don't
   * there is no way to tell whats in the variable's listent.
   */

  if (args[0].type != VAR)
    rerror ("members: argument must be a variable");

  var = (ListNode *) args[0].u.ptr;
  if (!var)
  {
    return ent_Assign_Rlab_MDS(NULL);
  }
  e = var_ent (var);
  if (!e)
  {
    return ent_Assign_Rlab_MDS(NULL);
  }
  type = ent_type (e);
  if (type == UNDEF)
  {
    return ent_Assign_Rlab_MDS(NULL);
  }

  /* Get the object's attributes (elements). */
  vfptr = (VFPTR) members_method[type].op;
  if (vfptr == 0)
  {
    ent_Clean (e);
    return ent_Assign_Rlab_MDS(NULL);
  }

  rval = (*vfptr) (ent_data (e), &n);
  smat = mds_MakeCharMatrix (rval, n);

  /* Now, go back, and get the variable's listent if it exists. */
  if ((listent = var_listent (var)))
  {
    rval2 = btree_Members (ent_data (listent), &n);
    smat2 = mds_MakeCharMatrix (rval2, n);
    smat3 = mds_Append (smat, smat2);
    mds_Destroy (smat);
    mds_Destroy (smat2);
    smat = smat3;
  }

  ent_Clean (e);

  return ent_Assign_Rlab_MDS(smat);
}

/* **************************************************************
 * Get the size of an object.
 * ************************************************************** */

Ent *
Size (int nargs, Datum args[])
{
  int type, rtype;
  void *rval;
  void *(*vfptr) ();
  Ent *e=0;

  if (nargs != 1)
  {
    rerror ("size: 1 argument allowed");
  }

  e = bltin_get_ent (args[0]);
  type = ent_type (e);

  rtype = size_method[type].type;
  vfptr = (VFPTR) size_method[type].op;

  if (vfptr == 0)
  {
    fprintf (stderr, "Entity type: %s\n", etd (e));
    rerror ("size() operation not supported");
  }

  rval = (*vfptr) (ent_data (e));

  ent_Clean (e);

  return ent_Assign_Rlab_Rtype(rval, rtype);
}

/* **************************************************************
 * Get the Length of an object.
 * ************************************************************** */

Ent *
Length (int nargs, Datum args[])
{
  int type, rtype=UNDEF;
  void *rval;
  void *(*vfptr) ();
  Ent *rent=0, *e=0;

  if (nargs == 1)
  {
    e = bltin_get_ent (args[0]);
    type = ent_type (e);

    rtype = length_method[type].type;
    vfptr = (VFPTR) length_method[type].op;

    if (vfptr)
    {
      rval = (*vfptr) (ent_data (e));
      rent = ent_Assign_Rlab_Rtype(rval, rtype);
    }
  }

  ent_Clean (e);

  if (rent)
    return (rent);

  return ent_Create_Rlab_Error();
}

/* **************************************************************
 * Get the type of an object.
 * ************************************************************** */

Ent *
Type (int nargs, Datum args[])
{
  int type, rtype;
  void *rval;
  void *(*vfptr) ();
  Ent *e=0;

  if (nargs != 1)
  {
    rerror ("type: 1 argument allowed");
  }

  e = bltin_get_ent (args[0]);
  type = ent_type (e);

  rtype = type_method[type].type;
  vfptr = (VFPTR) type_method[type].op;

  if (vfptr == 0)
  {
    fprintf (stderr, "Entity type: %s\n", etd (e));
    rerror ("type() operation not supported");
  }

  rval = (*vfptr) (ent_data (e));

  ent_Clean (e);

  return ent_Assign_Rlab_Rtype(rval, rtype);
}

/* **************************************************************
 * Load an rlab file:
 *    comment [mk,2017-10-12]: the complication below is a result of the
 *    following scenario:
 *    load(fn) where fn contains name of the script in which exist the same
 *    variable fn with different value. then rlab used to crash as when exiting
 *    the load command the underlying entity has been changed.
 * ********
****************************************************** */

#undef  THIS_SOLVER
#define THIS_SOLVER "load"
extern int run_program (char *input);

Ent * Load (int nargs, Datum args[])
{
  Ent *e1=0;
  int rval=RLAB_STATUS_FAILURE;
  FILE *rlab_stderr = (!RLAB_STDERR_DS) ? stderr : RLAB_STDERR_DS;
  char *fn=0, *stmp=0;

  /* Check n_args */
  if (nargs != 1)
  {
    fprintf (rlab_stderr, THIS_SOLVER ": " RLAB_ERROR_ARG1_MDS_SCALAR "\n");
    ent_Clean (e1);
    return ent_Create_Rlab_Double(rval);
  }

  /* get arg from list */
  e1 = bltin_get_ent (args[0]);
  RLABCODE_PROCESS_ARG_MD_IGNORE_ERROR(THIS_SOLVER,0,e1,MATRIX_DENSE_STRING,fn,class_char_pointer,isvalidstring,<,1,NULL);
  if (fn)
    stmp = cpstr(fn);

  ent_Clean(e1); // get rid of e1 before executing the script: in the script e1 might change
  
  if (stmp)
  {
    run_program (stmp);       // Run rlab on the file
    rval = RLAB_STATUS_SUCCESS;
    GC_FREE(stmp);
  }

  return ent_Create_Rlab_Double(rval);
}

/* **************************************************************
 * Debug an rlab file...
 * ************************************************************** */

extern int run_program_debug (char *input);
int Rlab_Debug = 0;

Ent *
Debug (int nargs, Datum args[])
{
  char *stmp;
  Ent *e=0;

  /* Check n_args */
  if (nargs != 1)
  {
    fprintf (stderr, "debug: 1 argument allowed");
    return ent_Create_Rlab_Failure();
  }

  /* get arg from list */
  e = bltin_get_ent (args[0]);
  stmp = class_char_pointer (e);

  /* Run rlab on the file */
  Rlab_Debug = 1;
  run_program_debug (stmp);
  close_file_ds (stmp);
  Rlab_Debug = 0;

  ent_Clean (e);

  return ent_Create_Rlab_Success();
}

/* **************************************************************
 * Reshape an object.
 * ************************************************************** */

Ent *
Reshape (int nargs, Datum args[])
{
  int type, rtype;
  void *rval;
  void *(*vfptr) ();
  Ent *rent, *e, *e2, *e3;
  int nrow, ncol;

  if (nargs != 3)
  {
    rerror ("reshape: 3 arguments required");
  }

  e = bltin_get_ent (args[0]);
  type = ent_type (e);

  rtype = reshape_method[type].type;
  vfptr = (VFPTR) reshape_method[type].op;

  if (vfptr == 0)
  {
    fprintf (stderr, "Entity type: %s\n", etd (e));
    rerror ("reshape() operation not supported");
  }

  e2 = bltin_get_ent (args[1]);
  nrow = (int) class_double (e2);
  e3 = bltin_get_ent (args[2]);
  ncol = (int) class_double (e3);

  rval = (*vfptr) (ent_data (e), nrow, ncol);

  rent = ent_Create ();
  ent_data (rent) = rval;
  ent_type (rent) = rtype;

  ent_Clean (e);
  return (rent);
}

/* **************************************************************
 * Create a dense matrix of zeros...
 * ************************************************************** */
#undef  THIS_SOLVER
#define THIS_SOLVER "zeros"
Ent *
Zeros (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0;
  MD *m=0;
  int nr=1, nc=1;

  if (nargs == 1)
  {
    e1 = bltin_get_ent (args[0]);
    if ( (ent_type (e1) == MATRIX_DENSE_REAL) ||
          (ent_type (e1) == MATRIX_DENSE_COMPLEX) ||
          (ent_type (e1) == MATRIX_DENSE_STRING) )
    {
      m = (MD *) ent_data (e1);
      nr = MNR (m);
      nc = MNC (m);
    }
    else
      rerror (THIS_SOLVER ": " RLAB_ERROR_ARG1_MD "\n");
  }
  else if (nargs == 2)
  {
    e1 = bltin_get_ent (args[0]);
    e2 = bltin_get_ent (args[1]);

    if (ent_type (e1) != MATRIX_DENSE_REAL || ent_type (e2) != MATRIX_DENSE_REAL)
      rerror ("zeros: arguments must be real matrices!");
    nr = (int) class_double (e1);
    nc = (int) class_double (e2);

    if(nr < 0 || nc < 0 )
    {
      rerror ("zeros: dimensions of the matrix must be positive real integers!");
    }
  }

  ent_Clean (e1);
  ent_Clean (e2);

  m = mdr_Create (nr, nc);
  mdr_Zero (m);

  return ent_Assign_Rlab_MDR(m);
}


#undef  THIS_SOLVER
#define THIS_SOLVER "ones"
Ent * Ones (int nargs, Datum args[])
{
  // original by ian searle,
  // modification by marijan kostrun
  Ent *e1=0, *e2=0;
  MD *m = 0;
  int  nr = 1, nc = 1, i;

  if (nargs == 1)
  {
    e1 = bltin_get_ent (args[0]);
    if ( (ent_type (e1) == MATRIX_DENSE_REAL) ||
          (ent_type (e1) == MATRIX_DENSE_COMPLEX) ||
          (ent_type (e1) == MATRIX_DENSE_STRING) )
    {
      m = (MD *) ent_data (e1);
      nr = MNR (m);
      nc = MNC (m);
    }
    else
      rerror (THIS_SOLVER ": " RLAB_ERROR_ARG1_MD "\n");
  }
  else if (nargs == 2)
  {
    e1 = bltin_get_ent (args[0]);
    e2 = bltin_get_ent (args[1]);
    if (ent_type (e1) != MATRIX_DENSE_REAL || ent_type (e2) != MATRIX_DENSE_REAL)
      rerror ("zeros: arguments must be real matrices!");

    nr = (int) class_double (e1);
    nc = (int) class_double (e2);
    if(nr < 1 || nc < 1 )
    {
      nr = 0;
      nc = 0;
    }
  }

  ent_Clean (e1);
  ent_Clean (e2);

  m = mdr_Create (nr, nc);
  for (i=0; i<SIZE(m); i++)
    MdrV0 (m, i) = 1.0;

  return ent_Assign_Rlab_MDR(m);
}


//
// Create matrix of inf's and return it
//
#undef  THIS_SOLVER
#define THIS_SOLVER "inf"
Ent * Inf (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0;
  MD *m=0;
  int nr=1, nc=1;

  if (nargs == 1)
  {
    e1 = bltin_get_ent (args[0]);
    if ( (ent_type (e1) == MATRIX_DENSE_REAL) ||
          (ent_type (e1) == MATRIX_DENSE_COMPLEX) ||
          (ent_type (e1) == MATRIX_DENSE_STRING) )
    {
      m = (MD *) ent_data (e1);
      nr = MNR (m);
      nc = MNC (m);
    }
    else
      rerror (THIS_SOLVER ": " RLAB_ERROR_ARG1_MD "\n");
  }
  else if (nargs == 2)
  {
    e1 = bltin_get_ent (args[0]);
    e2 = bltin_get_ent (args[1]);
    if (ent_type (e1) != MATRIX_DENSE_REAL || ent_type (e2) != MATRIX_DENSE_REAL)
      rerror ("inf: arguments must be real matrices!");
    nr = (int) class_double (e1);
    nc = (int) class_double (e2);
  }

  ent_Clean (e1);
  ent_Clean (e2);

  return ent_Assign_Rlab_MDR(mdr_CreateInf (nr, nc));
}

//
// create matrix of nan's
//
#undef  THIS_SOLVER
#define THIS_SOLVER "nan"
Ent * Nan (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0;
  MD *m=0;
  int nr = 1, nc = 1;

  if (nargs == 1)
  {
    e1 = bltin_get_ent (args[0]);
    if ( (ent_type (e1) == MATRIX_DENSE_REAL) ||
          (ent_type (e1) == MATRIX_DENSE_COMPLEX) ||
          (ent_type (e1) == MATRIX_DENSE_STRING) )
    {
      m = (MD *) ent_data (e1);
      nr = MNR (m);
      nc = MNC (m);
    }
    else
      rerror (THIS_SOLVER ": " RLAB_ERROR_ARG1_MD "\n");
  }
  if (nargs == 2)
  {
    e1 = bltin_get_ent (args[0]);
    if (ent_type (e1) == MATRIX_DENSE_REAL)
      nr = (int) class_double (e1);

    e2 = bltin_get_ent (args[1]);
    if (ent_type (e2) == MATRIX_DENSE_REAL)
      nc = (int) class_double (e2);
  }

  ent_Clean(e1);
  ent_Clean(e2);

  return ent_Assign_Rlab_MDR(mdr_CreateNan (nr, nc));
}

#undef  THIS_SOLVER
#define THIS_SOLVER "isinf"
Ent * IsInf (int nargs, Datum args[])
{
  int type, rtype;
  void *rval=0;
  void *(*vfptr) ();
  Ent *e=0;

  if (nargs != 1)
  {
    rerror ("isinf: 1 argument allowed");
  }

  e = bltin_get_ent (args[0]);
  type = ent_type (e);

  rtype = isinf_method[type].type;
  vfptr = (VFPTR) isinf_method[type].op;

  if (vfptr == 0)
  {
    fprintf (stderr, "Entity type: %s\n", etd (e));
    rerror ("isinf() operation not supported");
  }

  rval = (*vfptr) (ent_data (e));

  ent_Clean (e);

  return ent_Assign_Rlab_Rtype (rval, rtype);
}

Ent * IsNan (int nargs, Datum args[])
{
  int type, rtype;
  void *rval=0;
  void *(*vfptr) ();
  Ent *e=0;

  if (nargs != 1)
  {
    rerror ("isnan: 1 argument allowed");
  }

  e = bltin_get_ent (args[0]);
  type = ent_type (e);

  rtype = isnan_method[type].type;
  vfptr = (VFPTR) isnan_method[type].op;

  if (vfptr == 0)
  {
    fprintf (stderr, "Entity type: %s\n", etd (e));
    rerror ("isnan() operation not supported");
  }

  rval = (*vfptr) (ent_data (e));

  ent_Clean (e);
  return ent_Assign_Rlab_Rtype (rval, rtype);
}

#undef  THIS_SOLVER
#define THIS_SOLVER "finite"
Ent * Finite (int nargs, Datum args[])
{
  Ent *e1=0;
  MDR *rval=0;

  if (nargs != 1)
  {
    rerror (THIS_SOLVER ":" RLAB_ERROR_ARG_1 );
  }

  e1 = bltin_get_ent (args[0]);
  if (ent_type(e1) == MATRIX_DENSE_REAL)
  {
    rval = mdr_Finite(ent_data(e1));
  }
  else if (ent_type(e1) == MATRIX_DENSE_COMPLEX)
  {
    rval = mdc_Finite(ent_data(e1));
  }
  else if (ent_type(e1) == MATRIX_DENSE_STRING)
  {
    rval = mdr_Create_SameSize(ent_data(e1));
    int i;
    for (i=0; i<SIZE(rval);i++)
    {
      MdrV0(rval,i) = -1.0;
    }
  }
  else
  {
    rval = mdr_CreateScalar(-1.0);
  }

  ent_Clean (e1);

  return ent_Assign_Rlab_MDR (rval);
}

/* **************************************************************
 * Clear (delete) a variable.
 * ************************************************************** */
#undef THIS_SOLVER
#define THIS_SOLVER "clear"
Ent * Clear (int nargs, Datum args[])
{
  int i = 1;
  int rval = 1;
  Datum d;
  ListNode *var;
  Ent *e=0, *enew=0;

  for (i = 0; i < nargs; i++)
  {
    d = args[i];

#ifdef DEBUG_MDE
    fprintf(stderr, "enttype = %i\n", d.type);
#endif

    if (d.type == VAR)
    {
      var = (ListNode *) d.u.ptr;
      e = var_ent (var);

      // skip protected lists
      if (ent_type (e) == BTREE)
      {
        Btree *bt = ent_data (e);
        if (bt->isconst)
        {
          ent_Clean (e);
          continue;
        }
      }

      if (ent_type(e)!=UNDEF)
      {
        ent_Destroy (e);
        rval = 0;
        if (var_listent (var))
        {
          ent_Destroy (var_listent (var));
          var->listent = 0;
        }
        enew = ent_Create ();
        ent_SetType (enew, UNDEF);
        ent_IncRef (enew);
        listNode_AttachEnt (var, enew);
      }
    }
    else if (d.type == ENTITY)
    {
      ent_Destroy (d.u.ptr);
      d.u.ptr = 0;
    }
  }

  return ent_Create_Rlab_Double(rval);
}

/* **************************************************************
 * Sum function...
 * ************************************************************** */
#undef  THIS_SOLVER
#define THIS_SOLVER "sum"
#undef  METHOD
#define METHOD sum_method
Ent *
Sum (int nargs, Datum args[])
{
  int type=UNDEF, rtype;
  void *rval;
  void *(*vfptr) ();
  Ent *e1=0, *e2=0;
  void *cond=0;

  if (nargs!=1 && nargs!=2)
  {
    rerror (THIS_SOLVER ": 1 or 2 arguments allowed");
  }

  // first argument:
  e1 = bltin_get_ent (args[0]);
  if (e1)
    type = ent_type (e1);

  // second argument:
  //  for string matrices it is a column-separator during joining
  if (nargs==2)
  {
    e2 = bltin_get_ent (args[1]);
    if ((ent_type(e1) == MATRIX_DENSE_STRING) && (ent_type(e2) == MATRIX_DENSE_STRING))
    {
      cond = (void *) class_char_pointer(e2);
      if (isvalidstring(cond) < 0)
        cond = 0;
    }
    else if ((ent_type(e1) == MATRIX_DENSE_REAL || ent_type(e1) == MATRIX_DENSE_COMPLEX) && (ent_type(e2) == MATRIX_DENSE_REAL))
      cond = (void *) ent_data(e2);
  }

  rtype = METHOD[type].type;
  vfptr = (VFPTR) METHOD[type].op;

  if (vfptr == 0)
  {
    fprintf (stderr, "Entity type: %s\n", etd (e1));
    rerror (THIS_SOLVER "() operation not supported");
  }

  rval = (*vfptr) (ent_data (e1), cond);

  ent_Clean (e1);
  ent_Clean (e2);

  return ent_Assign_Rlab_Rtype(rval, rtype);
}

/* **************************************************************
 * CumSum function...
 * ************************************************************** */

Ent *
CumSum (int nargs, Datum args[])
{
  int type, rtype;
  void *rval;
  void *(*vfptr) ();
  Ent *e=0;

  if (nargs != 1)
  {
    rerror ("cumsum: 1 argument allowed");
  }

  e = bltin_get_ent (args[0]);
  type = ent_type (e);

  rtype = cumsum_method[type].type;
  vfptr = (VFPTR) cumsum_method[type].op;

  if (vfptr == 0)
  {
    fprintf (stderr, "Entity type: %s\n", etd (e));
    rerror ("cumsum() operation not supported");
  }

  rval = (*vfptr) (ent_data (e));

  ent_Clean (e);

  return ent_Assign_Rlab_Rtype(rval, rtype);
}

/* **************************************************************
 * Find function...
 * ************************************************************** */

Ent * Find (int nargs, Datum args[])
{
  int type, rtype;
  void *rval;
  void *(*vfptr) ();
  Ent *e=0;

  if (nargs != 1)
  {
    rerror ("find: 1 argument allowed");
  }

  e = bltin_get_ent (args[0]);
  type = ent_type (e);

  rtype = find_method[type].type;
  vfptr = (VFPTR) find_method[type].op;

  if (vfptr == 0)
  {
    fprintf (stderr, "Entity type: %s\n", etd (e));
    rerror ("find() operation not supported");
  }

  rval = (*vfptr) (ent_data (e));

  ent_Clean (e);

  return ent_Assign_Rlab_Rtype(rval, rtype);
}


#undef THIS_SOLVER
#define THIS_SOLVER "findvec"
#undef METHOD
#define METHOD findvec_method
Ent *
FindVector (int nargs, Datum args[])
{
  int type1=UNDEF, type2=UNDEF;
  MDR *w;
  void *(*vfptr) ();
  Ent *e1=0, *e2=0;

  if (nargs != 2)
  {
    rerror (THIS_SOLVER ": 2 argument allowed");
  }

  e1 = bltin_get_ent (args[0]);
  type1 = ent_type (e1);

  e2 = bltin_get_ent (args[1]);
  type2 = ent_type (e2);

  vfptr = (VFPTR) METHOD[type1][type2].op;

  if (vfptr == 0)
  {
    fprintf (stderr, THIS_SOLVER ": Entities type: %s and %s\n", etd (e1), etd(e2));
    rerror (THIS_SOLVER ": operation not supported");
  }

  w = (MDR*) (*vfptr) (ent_data (e1),ent_data (e2));

  ent_Clean (e1);
  ent_Clean (e2);

  return ent_Assign_Rlab_MDR(w);
}

/* **************************************************************
 * Sort function...
 * ************************************************************** */
#undef THIS_SOLVER
#define THIS_SOLVER "sort"
#undef METHOD
#define METHOD sort_method
Ent *
Sort (int nargs, Datum args[])
{
  int type=UNDEF, old_sort=rlab_sort_get_method(), old_nans=rlab_sort_get_nans();
  void *rval=0;
  void *(*vfptr) ();
  Ent *e1=0, *e2=0;

  if (nargs < 1 || nargs > 2)
  {
    fprintf(stderr, THIS_SOLVER ": " RLAB_ERROR_ONE_OR_TWO_ARG_REQUIRED "\n");
    goto _exit_sort;
  }

  e1 = bltin_get_ent (args[0]);
  if (e1)
    type = ent_type (e1);

  vfptr = (VFPTR) METHOD[type].op;
  if (!vfptr)
  {
    fprintf (stderr, THIS_SOLVER ": Argument type '%s' is not supported!\n", etd (e1));
    goto _exit_sort;
  }

  if (nargs>1)
  {
    e2 = bltin_get_ent (args[1]);

    if (ent_type(e2)==BTREE)
    {
      ListNode *node=0;
      char *s_method=0, *s_nan=0;

      RLABCODE_PROCESS_BTREE_ENTRY_S(e2,node,RLAB_NAME_GEN_METHOD,s_method,1,0);
      if (s_method)
      {
        if (strcmp (s_method, RLAB_SORT_NAME_QUICKSORT)==0)
        {
          rlab_sort_set_method(RLAB_SORT_QUICK);
        }
        else if (strcmp (s_method, RLAB_SORT_NAME_HEAPSORT)==0)
        {
          rlab_sort_set_method(RLAB_SORT_HEAP);
        }
      }

      RLABCODE_PROCESS_BTREE_ENTRY_S(e2,node,RLAB_SORT_NAN_METHOD,s_nan,1,0);
      if (s_nan)
      {
        if (strcmp (s_nan, RLAB_SORT_NAN_NAME_INPLACE)==0)
        {
          rlab_sort_set_nans(RLAB_SORT_NAN_INPLACE);
        }
        else if (strcmp (s_nan, RLAB_SORT_NAN_NAME_ONBOTTOM)==0)
        {
          rlab_sort_set_nans(RLAB_SORT_NAN_ONBOTTOM);
        }
        else if (strcmp (s_nan, RLAB_SORT_NAN_NAME_ONTOP)==0)
        {
          rlab_sort_set_nans(RLAB_SORT_NAN_ONTOP);
        }
      }

    }
  }

  rval = (*vfptr) (ent_data (e1));

_exit_sort:

  // restore earlier state of sort 
  rlab_sort_set_method(old_sort);
  rlab_sort_set_nans(old_nans);

  // clean, clean, clean!
  ent_Clean (e1);
  ent_Clean (e2);
  return ent_Assign_Rlab_BTREE(rval);
}

#undef  THIS_SOLVER
#undef  METHOD
#define THIS_SOLVER "unique"
#define METHOD vectorset_method
Ent * VectorSet (int nargs, Datum args[])
{
  int type=UNDEF, rtype=UNDEF;
  void *rval=0;
  void *(*vfptr) ();
  Ent *rent, *e1=0;

  if (nargs != 1)
    rerror (THIS_SOLVER ": 1 argument allowed");

  e1 = bltin_get_ent (args[0]);
  if (e1)
    type = ent_type (e1);

  rtype = METHOD[type].type;
  vfptr = (VFPTR) METHOD[type].op;

  if (vfptr == 0)
  {
    fprintf (stderr, "Entity type: %s\n", etd (e1));
    rerror (THIS_SOLVER "() operation not supported");
  }

  rval = (*vfptr) (ent_data (e1));

  ent_Clean (e1);

  rent = ent_Create ();
  ent_data (rent) = rval;
  ent_type (rent) = rtype;
  return (rent);
}


#undef  THIS_SOLVER
#undef  METHOD
#define THIS_SOLVER "union"
#define METHOD vectorunion_method
Ent *
VectorUnion (int nargs, Datum args[])
{
  int type=UNDEF, rtype=UNDEF;
  void *rval=0;
  void *(*vfptr) ();
  Ent *e1=0, *e2=0;

  if (nargs != 2)
    rerror (THIS_SOLVER ": 2 argument allowed");

  e1 = bltin_get_ent (args[0]);
  e2 = bltin_get_ent (args[1]);

  type = ent_type (e1);
  if (ent_type(e2) != type)
    rerror (THIS_SOLVER ": only MDR arguments allowed!");

  rtype = METHOD[type].type;
  vfptr = (VFPTR) METHOD[type].op;

  if (vfptr == 0)
  {
    fprintf (stderr, "Entity type: %s and %s\n", etd (e1), etd (e2));
    rerror (THIS_SOLVER "() operation not supported");
  }

  rval = (*vfptr) (ent_data (e1), ent_data (e2));

  ent_Clean (e1);
  ent_Clean (e2);

  return ent_Assign_Rlab_Rtype(rval, rtype);
}

#undef  THIS_SOLVER
#undef  METHOD
#define THIS_SOLVER "intersect"
#define METHOD vectorintersect_method
Ent *
VectorIntersect (int nargs, Datum args[])
{
  int type=UNDEF, rtype=UNDEF;
  void *rval=0;
  void *(*vfptr) ();
  Ent *e1=0, *e2=0;

  if (nargs != 2)
    rerror (THIS_SOLVER ": 2 argument allowed");

  e1 = bltin_get_ent (args[0]);
  e2 = bltin_get_ent (args[1]);

  type = ent_type (e1);
  if (ent_type(e2) != type)
    rerror (THIS_SOLVER ": only MDR arguments allowed!");

  rtype = METHOD[type].type;
  vfptr = (VFPTR) METHOD[type].op;

  if (vfptr == 0)
  {
    fprintf (stderr, "Entity type: %s and %s\n", etd (e1), etd (e2));
    rerror (THIS_SOLVER "() operation not supported");
  }

  rval = (*vfptr) (ent_data (e1), ent_data (e2));

  ent_Clean (e1);
  ent_Clean (e2);

  return ent_Assign_Rlab_Rtype(rval, rtype);
}

#undef  THIS_SOLVER
#undef  METHOD
#define THIS_SOLVER "complement"
#define METHOD vectorcomplement_method
Ent *
VectorComplement (int nargs, Datum args[])
{
  int type=UNDEF, rtype=UNDEF;
  void *rval=0;
  void *(*vfptr) ();
  Ent *e1=0, *e2=0;

  if (nargs != 2)
    rerror (THIS_SOLVER ": " RLAB_ERROR_TWO_ARG_REQUIRED);

  e1 = bltin_get_ent (args[0]);
  if (!e1)
    rerror (THIS_SOLVER ": " RLAB_ERROR_TWO_ARG_REQUIRED);
  type = ent_type (e1);

  e2 = bltin_get_ent (args[1]);
  if (!e2)
    rerror (THIS_SOLVER ": " RLAB_ERROR_TWO_ARG_REQUIRED);

  if (ent_type(e2) != type)
    rerror (THIS_SOLVER ": " RLAB_ERROR_ARG1_ARG2_SAME_TYPE);

  rtype = METHOD[type].type;
  vfptr = (VFPTR) METHOD[type].op;

  if (vfptr == 0)
  {
    fprintf (stderr, "Entity type: %s and %s\n", etd (e1), etd (e2));
    rerror (THIS_SOLVER "() operation not supported");
  }

  rval = (*vfptr) (ent_data (e1), ent_data (e2));

  ent_Clean (e1);
  ent_Clean (e2);

  return ent_Assign_Rlab_Rtype(rval, rtype);
}

#undef  THIS_SOLVER
#undef  METHOD
#define THIS_SOLVER "merge"
#define METHOD merge_method
Ent * Merge (int nargs, Datum args[])
{
  int type=UNDEF, rtype;
  void *rval=0;
  void *(*vfptr) ();
  Ent *e1=0, *e2=0;

  if (nargs != 2)
    rerror (THIS_SOLVER ": 2 argument allowed");

  e1 = bltin_get_ent (args[0]);
  e2 = bltin_get_ent (args[1]);

  type = ent_type (e1);
  if (ent_type(e2) != type)
    rerror (THIS_SOLVER ": only MDR arguments allowed!");

  rtype = METHOD[type].type;
  vfptr = (VFPTR) METHOD[type].op;

  if (vfptr == 0)
  {
    fprintf (stderr, "Entity type: %s and %s\n", etd (e1), etd (e2));
    rerror (THIS_SOLVER "() operation not supported");
  }

  rval = (*vfptr) (ent_data (e1), ent_data (e2));

  ent_Clean (e1);
  ent_Clean (e2);

  return ent_Assign_Rlab_Rtype(rval, rtype);
}

#undef  THIS_SOLVER
#undef  METHOD
#define THIS_SOLVER "compact"
#define METHOD compact_method
Ent * Compact (int nargs, Datum args[])
{
  int type=UNDEF;
  void *rval=0;
  void *(*vfptr) ();
  Ent *e1=0, *e2=0,*e3=0,*e4=0;

  MDR *x=0, *sx=0, *wx=0;

  int col_idx = 1;

  if (nargs < 1 || nargs>4)
  {
    fprintf(stderr, THIS_SOLVER ": " RLAB_ERROR_ONE_TO_FOUR_ARG_REQUIRED "\n");
    goto _exit_compact;
  }

  //
  // first argument: data matrix
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type(e1)!=MATRIX_DENSE_REAL)
  {
    fprintf(stderr, THIS_SOLVER ": " RLAB_ERROR_ARG1_MDR_MATRIX "\n");
    goto _exit_compact;
  }
  type = ent_type (e1);
  x = ent_data(e1);
  if (EQNULL(x))
  {
    fprintf(stderr, THIS_SOLVER ": " RLAB_ERROR_ARG1_MDR_MATRIX "\n");
    goto _exit_compact;
  }

  //
  // col_index
  //
  if (nargs>=2)
  {
    e2 = bltin_get_ent (args[1]);
    if (ent_type(e2)==MATRIX_DENSE_REAL)
    {
      col_idx = (int) class_double(e2);
    }
  }
  if (col_idx < 1 || col_idx > MNC(x))
  {
    fprintf(stderr, THIS_SOLVER ": column index out of range!\n");
    goto _exit_compact;
  }

  //
  // std_index
  //
  if (nargs>=3)
  {
    e3 = bltin_get_ent (args[2]);
    if (ent_type(e3)==MATRIX_DENSE_REAL)
    {
      sx = ent_data(e3);
      if (MNC(sx)!=2)
        sx=0;
      if (sx)
        if(SIZE(sx)!=MNC(x)-1)
        {
          fprintf(stderr, THIS_SOLVER ": improperly sized array of value/deviation pairs!\n");
          goto _exit_compact;
        }
    }
  }

  //
  // wgt_index
  //
  if (nargs>=4)
  {
    e4 = bltin_get_ent (args[3]);
    if (ent_type(e4)==MATRIX_DENSE_REAL)
    {
      wx = ent_data(e4);
      if (MNC(wx)!=2)
        wx=0;
      if (wx)
        if(SIZE(wx)!=MNC(x)-1)
        {
          fprintf(stderr, THIS_SOLVER ": improperly sized array of value/weight pairs!\n");
          goto _exit_compact;
        }
    }
  }

  vfptr = (VFPTR) METHOD[type].op;
  if (vfptr)
    rval = (*vfptr) (x, col_idx, wx, sx);
  else
    fprintf (stderr, THIS_SOLVER ": Argument type '%s' is not supported!\n", etd (e1));

_exit_compact:

  ent_Clean (e1);
  ent_Clean (e2);
  ent_Clean (e3);
  ent_Clean (e4);

  return ent_Assign_Rlab_MDR(rval);
}

/* **************************************************************
 * Strtol (string-to-long-int) function...
 * ************************************************************** */

Ent * Strtol (int nargs, Datum args[])
{
  int base, type, rtype;
  void *rval;
  void *(*vfptr) ();
  Ent *e=0, *e2=0;

  if (nargs > 2 || nargs < 1)
  {
    rerror ("strtol: 1 or 2 arguments required");
  }

  e = bltin_get_ent (args[0]);
  type = ent_type (e);

  if (nargs == 2)
  {
    e2 = bltin_get_ent (args[1]);
    base = (int) class_double (e2);

    if (base == 0 || (base >= 2 && base <= 32))
      ;       /* OK */
    else
      rerror ("strtol: base must be 0, or 2 <= base <= 32");
  }
  else
  {
    base = 10;      /* Base-10 conversion, default. */
  }

  rtype = strtol_method[type].type;
  vfptr = (VFPTR) strtol_method[type].op;

  if (vfptr == 0)
  {
    fprintf (stderr, "Entity type: %s\n", etd (e));
    rerror ("strtol() operation not supported");
  }

  rval = (*vfptr) (ent_data (e), base);

  ent_Clean (e);

  return ent_Assign_Rlab_Rtype(rval, rtype);
}


/* **************************************************************
 * Sign function.
 * ************************************************************** */

Ent * Sign (int nargs, Datum args[])
{
  int type, rtype;
  void *rval;
  void *(*vfptr) ();
  Ent *e=0;

  if (nargs != 1)
  {
    rerror ("sign: 1 argument allowed");
  }

  e = bltin_get_ent (args[0]);
  type = ent_type (e);

  rtype = sign_method[type].type;
  vfptr = (VFPTR) sign_method[type].op;

  if (vfptr == 0)
  {
    fprintf (stderr, "Entity type: %s\n", etd (e));
    rerror ("sign() operation not supported");
  }

  rval = (*vfptr) (ent_data (e));

  ent_Clean (e);

  return ent_Assign_Rlab_Rtype(rval, rtype);
}

/* **************************************************************
 * Sizeof function.
 * ************************************************************** */

Ent * Sizeof (int nargs, Datum args[])
{
  Ent *e=0;

  if (nargs != 1)
  {
    rerror ("sizeof: 1 argument allowed");
  }

  e = convert_datum_to_ent (args[0]);

  ent_Clean (e);

  return ent_Create_Rlab_Double(class_sizeof (e));
}

#if !defined (winnt)

Ent * System (int nargs, Datum args[])
{
  int i, retval=0;
  Ent *e1=0;
  MDS *s=0;

  /* Check nargs */
  if (nargs != 1)
    rerror ("system: 1 argument allowed");

  e1 = bltin_get_ent (args[0]);
  if (ent_type(e1) == MATRIX_DENSE_STRING)
  {
    retval=0;
    s = ent_data(e1);
    for (i=0; i < SIZE(s); i++)
    {
      if (isvalidstring(MdsV0(s,i)) > 0)
        retval |= (int) system ( MdsV0(s,i) );
    }
  }

  ent_Clean (e1);

  return ent_Create_Rlab_Double(retval);
}

Ent * Fork (int nargs, Datum args[])
{
  int i, retval=0;
  Ent *e1=0;
  MDS *s=0;
  pid_t pid;

  /* Check nargs */
  if (nargs != 1)
    rerror ("fork: 1 argument allowed");

  e1 = bltin_get_ent (args[0]);
  if (ent_type(e1) == MATRIX_DENSE_STRING)
  {
    retval=0;
    s = ent_data(e1);
    if (SIZE(s) > 0)
    {
      for (i=0; i < SIZE(s); i++)
      {
        if (isvalidstring(MdsV0(s,i)) > 0)
        {
          pid = vfork();
          if(pid == 0)
          {
            // this is the child
            system( MdsV0(s,i) ); // start process
            exit(0);              // child exits
          }
          // parents does nothing
        }
      }
    }
  }

  ent_Clean (e1);

  return ent_Create_Rlab_Double(retval);
}

#else

/*
 * Windows-NT (Cygnus Win32 Stuff)
 */

#include <stdio.h>
#include <unistd.h>
#include <errno.h>

char *shell = "/bin/sh";

#undef  THIS_SOLVER
#define THIS_SOLVER "system"
Ent * System (int nargs, Datum args[])
{
  int pid, child_exit_status;
  double retval;
  char *s;
  unsigned retval;
  Ent *e=0;

  /* Check nargs */
  if (nargs != 1)
  {
    rerror (THIS_SOLVER ": 1 argument allowed");
  }

  e = bltin_get_ent (args[0]);
  s = class_char_pointer (e);
  if (isvalidstring(s)<0)
    goto _exit_system;

  switch (pid = fork ())
  {
    case -1:      /* fork failed */

      fprintf (stderr, THIS_SOLVER ": Could not create a new process, errno=%i\n", errno);
      retval = 127;
      break;

    case 0:     /* The child process. */
      execlp (shell, shell, "-c", "ls", (char *) 0);
      /* if get here, execl() failed */
      fprintf (stderr, THIS_SOLVER ": Execute of %s failed, errno=%i\n", shell, errno);
      fflush (stderr);
      rerror (THIS_SOLVER ": Fail!");

    default:      /* The parent process, wait for the child */
      retval = wait (&child_exit_status);
      break;
  }

_exit_system:

  ent_Clean (e);
  return ent_Create_Rlab_Double(retval);
}
#endif /* !winnt */

/* **************************************************************
 * Prod function...
 * ************************************************************** */

Ent * Prod (int nargs, Datum args[])
{
  int type, rtype;
  void *rval;
  void *(*vfptr) ();
  Ent *e=0;

  if (nargs != 1)
  {
    rerror ("prod: 1 argument allowed");
  }

  e = bltin_get_ent (args[0]);
  type = ent_type (e);

  rtype = prod_method[type].type;
  vfptr = (VFPTR) prod_method[type].op;

  if (vfptr == 0)
  {
    fprintf (stderr, "Entity type: %s\n", etd (e));
    rerror ("prod() operation not supported");
  }

  rval = (*vfptr) (ent_data (e));

  ent_Clean (e);

  return ent_Assign_Rlab_Rtype(rval,rtype);
}

/* **************************************************************
 * CumProd function...
 * ************************************************************** */

Ent * CumProd (int nargs, Datum args[])
{
  int type, rtype;
  void *rval;
  void *(*vfptr) ();
  Ent *e=0;

  if (nargs != 1)
  {
    rerror ("cumprod: 1 argument allowed");
  }

  e = bltin_get_ent (args[0]);
  type = ent_type (e);

  rtype = cumprod_method[type].type;
  vfptr = (VFPTR) cumprod_method[type].op;

  if (vfptr == 0)
  {
    fprintf (stderr, "Entity type: %s\n", etd (e));
    rerror ("cumprod() operation not supported");
  }

  rval = (*vfptr) (ent_data (e));

  ent_Clean (e);

  return ent_Assign_Rlab_Rtype (rval, rtype);
}

/* **************************************************************
 * frexp (ANSI) function...
 * ************************************************************** */

Ent * Frexp (int nargs, Datum args[])
{
  int type, rtype;
  void *rval;
  void *(*vfptr) ();
  Ent *e=0;

  if (nargs != 1)
  {
    rerror ("frexp: 1 argument allowed");
  }

  e = bltin_get_ent (args[0]);
  type = ent_type (e);

  rtype = frexp_method[type].type;
  vfptr = (VFPTR) frexp_method[type].op;

  if (vfptr == 0)
  {
    fprintf (stderr, "Entity type: %s\n", etd (e));
    rerror ("frexp() operation not supported");
  }

  rval = (*vfptr) (ent_data (e));

  ent_Clean (e);

  return ent_Assign_Rlab_Rtype (rval, rtype);
}

/* **************************************************************
 * ldexp (ANSI) function...
 * ************************************************************** */

Ent *Ldexp (int nargs, Datum args[])
{
  int type1, type2, rtype;
  void *rval;
  void *(*vfptr) ();
  Ent *e1=0, *e2=0;

  if (nargs != 2)
  {
    rerror ("ldexp: 2 arguments required");
  }

  e1 = bltin_get_ent (args[0]);
  type1 = ent_type (e1);

  e2 = bltin_get_ent (args[1]);
  type2 = ent_type (e2);

  rtype = ldexp_method[type1][type2].type;
  vfptr = (VFPTR) ldexp_method[type1][type2].op;

  if (vfptr == 0)
  {
    fprintf (stderr, "Entity types: %s and %s\n", etd (e1), etd (e2));
    rerror ("ldexp() operation not supported");
  }

  rval = (*vfptr) (ent_data (e1), ent_data (e2));

  ent_Clean (e1);
  ent_Clean (e2);

  return ent_Assign_Rlab_Rtype (rval, rtype);
}

#ifdef HAVE_LOGB
Ent * Logb (int nargs, Datum args[])
{
  int type, rtype;
  void *rval;
  void *(*vfptr) ();
  Ent *e=0;

  /* Check n_args */
  if (nargs != 1)
  {
    rerror ("logb: 1 argument allowed");
  }

  e = bltin_get_ent (args[0]);
  type = ent_type (e);

  rtype = logb_method[type].type;
  vfptr = (VFPTR) logb_method[type].op;

  if (vfptr == 0)
  {
    fprintf (stderr, "Entity type: %s\n", etd (e));
    rerror ("logb() operation not supported");
  }

  rval = (*vfptr) (ent_data (e));

  ent_Clean (e);

  return ent_Assign_Rlab_Rtype (rval, rtype);
}

#endif /* HAVE_LOGB */

/* **************************************************************
 * Timing functions (tic() and toc()).
 * ************************************************************** */

#include <time.h>

#define RLAB_TICTIMER_NO  32
#define RLAB_TICTIMER_DEFAULT 0

static struct timeval tictime[ RLAB_TICTIMER_NO ];
static char *tictime_id[ RLAB_TICTIMER_NO ] = { 0 };
static int calltic[ RLAB_TICTIMER_NO ] = { 0 };
static int rlab_default_tictimer = RLAB_TICTIMER_DEFAULT;

/*
 * Start the timer.
 */
Ent * Tic (int nargs, Datum args[])
{
  Ent *e1=0;
  MDR *xtic=0;

  int tic_no = rlab_default_tictimer;
  int i;

  /* Check n_args */
  if (nargs == 0)
  {
    gettimeofday( &tictime[tic_no], NULL);
    calltic[ tic_no ] = 1;
    if (tictime_id[ tic_no ])
    {
      GC_free(tictime_id[ tic_no ]);
      tictime_id[ tic_no ] = 0;
    }
  }
  else if (nargs == 1)
  {
    e1 = bltin_get_ent (args[0]);
    if (ent_type(e1) == MATRIX_DENSE_REAL)
    {
      xtic = class_matrix_real(e1);
      for (i=0; i<(xtic->nrow * xtic->ncol); i++)
      {
        tic_no = (int)(MdrV0(xtic,i) - 1.0);
        if (tic_no >= 0 && tic_no < RLAB_TICTIMER_NO)
        {
          gettimeofday( &tictime[tic_no], NULL);
          calltic[ tic_no ] = 1;
          if (tictime_id[ tic_no ])
          {
            GC_free(tictime_id[ tic_no ]);
            tictime_id[ tic_no ] = 0;
          }
        }
        else
        {
          rerror ("tic: index out of range!");
        }
      }
    }
  }

  ent_Clean (e1);

  return ent_Create_Rlab_Success();
}

/*
 * Report the elapsed time, in seconds,
 * since last call to tic().
 */
Ent * Toc (int nargs, Datum args[])
{
  Ent *e1=0;
  MDR *xtic=0, *elapsed=0;
  struct timeval toctime;

  int tic_no = rlab_default_tictimer;
  int i;

  /* Check n_args */
  if (nargs == 0)
  {
    if (calltic[tic_no])
    {
      elapsed = mdr_Create(1,1);

      gettimeofday(&toctime, NULL);
      MdrV0(elapsed,0) =
          (double) (toctime.tv_sec - tictime[tic_no].tv_sec) +
            ((double) (toctime.tv_usec - tictime[tic_no].tv_usec)) / 1e6;
    }
    else
    {
      elapsed = mdr_CreateScalar(-1.0);
    }
  }
  else if (nargs == 1)
  {
    //
    e1 = bltin_get_ent (args[0]);
    if (ent_type(e1) == MATRIX_DENSE_REAL)
    {
      xtic = class_matrix_real(e1);
      elapsed = mdr_Create(xtic->nrow, xtic->ncol);

      for (i=0; i<(xtic->nrow * xtic->ncol); i++)
      {
        tic_no = (int)(MdrV0(xtic,i) - 1.0);
        MdrV0(elapsed,i) = RLAB_STATUS_ERROR;
        if (tic_no >= 0 && tic_no < RLAB_TICTIMER_NO)
        {
          if (calltic[tic_no])
          {
            gettimeofday(&toctime, NULL);
            MdrV0(elapsed,i) =
              (double) (toctime.tv_sec - tictime[tic_no].tv_sec) +
                ((double) (toctime.tv_usec - tictime[tic_no].tv_usec)) / 1e6;
          }
        }
      }
    }
  }

  ent_Clean (e1);
  return ent_Assign_Rlab_MDR(elapsed);
}



/* **************************************************************
 * Sleep the RLaB process for a specified amount of time (relative)
 * or until a certain point in time (absolute)
 * ************************************************************** */
Ent * Sleep (int nargs, Datum args[])
{
  Ent *e1=0;
  MDR *x;

  struct timeval t;
  struct tm *time2;
  time_t tt1, tt2;

  double sleepval=0;

  if (nargs != 1)
    goto _exit_sleep;

  e1 = bltin_get_ent (args[0]);
  if (ent_type(e1)!=MATRIX_DENSE_REAL)
    goto _exit_sleep;

  x = ent_data(e1);

  if (EQSCAL(x))
  {
    // relative time
    sleepval = mdrV0(x,0);
  }
  else if (SIZE(x)==3)
  {
    // time now
    tt1 = time (NULL);
    // wake-up time
    time2 = localtime(&tt1);
    time2->tm_sec   =  mdiV0(x,2);
    time2->tm_min   =  mdiV0(x,1);
    time2->tm_hour  =  mdiV0(x,0);
    tt2 = mktime(time2);
    if (tt2 > tt1)
      sleepval = tt2 - tt1;
  }
  else if  (SIZE(x)==6)
  {
    // time now
    tt1 = time (NULL);
    // wake-up time
    time2 = localtime(&tt1);
    time2->tm_sec   =  mdiV0(x,5);
    time2->tm_min   =  mdiV0(x,4);
    time2->tm_hour  =  mdiV0(x,3);
    time2->tm_mday  =  mdiV0(x,2);
    time2->tm_mon   =  mdiV0(x,1);
    time2->tm_year  =  mdiV0(x,0);
    tt2 = mktime(time2);
    if (tt2 > tt1)
      sleepval = tt2 - tt1;
  }
  else if  (SIZE(x)==8)
  {
    // time now
    tt1 = time (NULL);
    // wake-up time
    time2 = localtime(&tt1);
    time2->tm_isdst  =  mdiV0(x,7);
    time2->tm_gmtoff = 3600 * mdiV0(x,6);
    time2->tm_sec   =  mdiV0(x,5);
    time2->tm_min   =  mdiV0(x,4);
    time2->tm_hour  =  mdiV0(x,3);
    time2->tm_mday  =  mdiV0(x,2);
    time2->tm_mon   =  mdiV0(x,1);
    time2->tm_year  =  mdiV0(x,0);
    tt2 = mktime(time2);
    if (tt2 > tt1)
      sleepval = tt2 - tt1;
  }

  if (sleepval)
  {
    t.tv_sec  = (long int) sleepval;
    t.tv_usec = (long int) (1e6 * (sleepval - (double) t.tv_sec) );
    select(0, NULL, NULL, NULL, &t);
  }

_exit_sleep:

  ent_Clean (e1);
  return ent_Create_Rlab_Double(sleepval);
}

/* **************************************************************
 * Change the numerical output format.
 * format()    should return the current values of the formats
 * ************************************************************** */

/* Control output print formats */
static int FWIDTH_DEFAULT = 9;
static int FPREC_DEFAULT = 3;
static int fwidth = 9;
static int fprec = 3;

int
get_fwidth (void)
{
  return (fwidth);
}

int
get_fprec (void)
{
  return (fprec);
}

Ent *
Format (int nargs, Datum args[])
{
  int size;
  Ent *e1=0, *e2=0, *rent;
  MDR *m, *rf=0;

  if (nargs > 2)
  {
    rerror ("format: takes two args at most");
  }
  else
  {

    /* Now reset the format(s). */
    if (nargs == 1)
    {
      /* Take a matrix, like size(). */
      e1 = bltin_get_ent (args[0]);
      m = class_matrix_real (e1);
      size = (int) (MNR (m) * MNC (m));

      if (size == 1)
      {
	/* Precision only. */
	fprec = (int) MdrV1 (m, 1);
      }
      else if (size > 1)
      {
	fwidth = (int) MdrV1 (m, 1);
	fprec = (int) MdrV1 (m, 2);
      }
    }
    else if (nargs == 2)
    {
      /* Width and Precision. */
      e1 = bltin_get_ent (args[0]);
      e2 = bltin_get_ent (args[1]);
      fwidth = (int) class_double (e1);
      fprec  = (int) class_double (e2);
    }
    else
    {
      /* Reset formats to default. */
      fwidth = FWIDTH_DEFAULT;
      fprec = FPREC_DEFAULT;
    }
  }

  /* Create the return matrix. */
  rf = mdr_Create (1, 2);
  Mdr1 (rf, 1, 1) = (double) fwidth;
  Mdr1 (rf, 1, 2) = (double) fprec;

  rent = ent_Create ();
  ent_data (rent) = rf;
  ent_type (rent) = MATRIX_DENSE_REAL;

  if (nargs >= 1)
    ent_Clean (e1);
  if (nargs == 2)
    ent_Clean (e2);

  return (rent);
}

// Ent *
// Tmpnam (int nargs, Datum args[])
// {
//   Ent *rent;
//   char fn[12] = "rlab_XXXXXX";
//   fn[11] = '\0';
//   MDS *ms;
//   if (nargs != 0)
//     rerror ("tmpnam: no arguments allowed");
//   tmpnam (fn);
//   ms = mds_CreateScalar (cpstr (fn));
//   rent = ent_Create ();
//   ent_data (rent) = ms;
//   ent_type (rent) = MATRIX_DENSE_STRING;
//   return (rent);
// }

#ifdef HAVE_GETENV
/* **************************************************************
 * RLaB interface to system getenv() function.
 * ************************************************************** */

Ent *
Getenv (int nargs, Datum args[])
{
  Ent *e=0;
  char *s=0, *retval=0;

  if (nargs != 1)
  {
    rerror ("getenv: requires 1 argument");
  }

  e = bltin_get_ent (args[0]);
  s = class_char_pointer (e);

  if (isvalidstring(s)>0)
    retval = getenv (s);

  ent_Clean (e);

  if (retval)
    return ent_Create_Rlab_String(retval);
  else
    return ent_Create_Rlab_String("");
}
#endif /* HAVE_GETENV */




#ifdef HAVE_PUTENV
/* **************************************************************
 * RLaB interface to system putenv() function.
 * ************************************************************** */

static char *putenv_string = 0;

Ent *
Putenv (int nargs, Datum args[])
{
  int retval;
  char *str;
  Ent *e=0;

  if (nargs != 1)
  {
    rerror ("putenv: requires 1 argument");
  }

  e = bltin_get_ent (args[0]);
  str = class_char_pointer (e);

  if (putenv_string != 0)
    GC_FREE (putenv_string);

  putenv_string = cpstr (str);

  // from linux man pages about putenv:
  // The string pointed to by string becomes part of the environment,
  // so altering the string changes the environment.
  // do not free 'putenv_string' !
  retval = putenv (putenv_string);

  ent_Clean (e);

  return ent_Create_Rlab_Double(retval);
}
#endif /* HAVE_PUTENV */

/* **************************************************************
 * RLaB interface to system chdir() function.
 * ************************************************************** */
#ifdef WIN32
#include <direct.h>
#define chdir _chdir
#endif

Ent *
Cd (int nargs, Datum args[])
{
  Ent *e, *rent;
  static char *cd_string;

  if (nargs != 1)
  {
    rerror ("cd: requires 1 argument");
  }

  e = bltin_get_ent (args[0]);
  if (cd_string == 0)
    GC_FREE (cd_string);
  cd_string = cpstr (class_char_pointer (e));

  if (chdir (cd_string))
  {
    switch (errno)
    {
    case EACCES:
      fprintf (stderr, "Search permission to: %sdenied\n", cd_string);
      errno = 0;
      break;

    case ENOTDIR:
      fprintf (stderr, "Part of the path is not a directory: %s\n", cd_string);
      errno = 0;
      break;

    case ENOENT:
      fprintf (stderr, "Part of the path does not exist: %s\n", cd_string);
      errno = 0;
      break;

    default:
      rerror ("cd: error during call");
    }
    rent = ent_Create ();
    ent_data (rent) = mdr_CreateScalar (0.0);
    ent_type (rent) = MATRIX_DENSE_REAL;
    return (rent);
  }

  rent = ent_Create ();
  ent_data (rent) = mdr_CreateScalar (1.0);
  ent_type (rent) = MATRIX_DENSE_REAL;

  ent_Clean (e);
  return (rent);
}

#undef  THIS_SOLVER
#define THIS_SOLVER "copy"
Ent * ent_copyentity (int nargs, Datum args[])
{
  Ent *e1=0, *rent;

  if (nargs != 1)
  { rerror (THIS_SOLVER ": 1 argument allowed"); }

  e1 = bltin_get_ent (args[0]);
  if (ent_type(e1) == UNDEF)
    rerror (THIS_SOLVER ": entity you are trying to copy is not defined");

  rent = ent_Duplicate ( e1 );
  rent->refc = 1;

  ent_Clean (e1);

  return (rent);
}

//
//
// B I T W I S E   O P E R A T O R S
//
//

#undef  THIS_SOLVER
#define THIS_SOLVER "bitshiftl"
Ent *
    ent_bit_shift_left (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *e3=0, *rval;
  MDR *a1=0, *a2=0, *a3=0;

  /* Check n_args */
  if (nargs != 3)
    rerror (THIS_SOLVER ": " RLAB_ERROR_THREE_ARG_REQUIRED);

  e1 = bltin_get_ent (args[0]);
  if (ent_type(e1) != MATRIX_DENSE_REAL)
    rerror (THIS_SOLVER ": " RLAB_ERROR_ARG1_MDR_INTEGER_MATRIX);
  a1 = ent_data(e1);

  e2 = bltin_get_ent (args[1]);
  if (ent_type(e2) != MATRIX_DENSE_REAL)
    rerror (THIS_SOLVER ": " RLAB_ERROR_ARG2_MDR_INTEGER_MATRIX);
  a2 = ent_data(e2);

  e3 = bltin_get_ent (args[2]);
  if (ent_type(e3) != MATRIX_DENSE_REAL)
    rerror (THIS_SOLVER ": " RLAB_ERROR_ARG3_MDR_INTEGER_MATRIX);
  a3 = ent_data(e3);

  rval = ent_Assign_Rlab_MDR(mdi_bit_shift_left(a1,a2,a3));

  ent_Clean(e1);
  ent_Clean(e2);
  ent_Clean(e3);

  return rval;
}

#undef  THIS_SOLVER
#define THIS_SOLVER "bitshiftr"
Ent *
    ent_bit_shift_right (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *e3=0, *rval;
  MDR *a1=0, *a2=0, *a3=0;

  /* Check n_args */
  if (nargs != 3)
    rerror (THIS_SOLVER ": " RLAB_ERROR_THREE_ARG_REQUIRED);

  e1 = bltin_get_ent (args[0]);
  if (ent_type(e1) != MATRIX_DENSE_REAL)
    rerror (THIS_SOLVER ": " RLAB_ERROR_ARG1_MDR_INTEGER_MATRIX);
  a1 = ent_data(e1);

  e2 = bltin_get_ent (args[1]);
  if (ent_type(e2) != MATRIX_DENSE_REAL)
    rerror (THIS_SOLVER ": " RLAB_ERROR_ARG2_MDR_INTEGER_MATRIX);
  a2 = ent_data(e2);

  e3 = bltin_get_ent (args[2]);
  if (ent_type(e3) != MATRIX_DENSE_REAL)
    rerror (THIS_SOLVER ": " RLAB_ERROR_ARG3_MDR_INTEGER_MATRIX);
  a3 = ent_data(e3);

  rval = ent_Assign_Rlab_MDR(mdi_bit_shift_right(a1,a2,a3));

  ent_Clean(e1);
  ent_Clean(e2);
  ent_Clean(e3);

  return rval;
}

#undef  THIS_SOLVER
#define THIS_SOLVER "bytesplit"
Ent * ent_byte_split (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *rval;
  int a, use_lsb=1;

  /* Check n_args */
  if (nargs != 1 && nargs != 2)
    rerror (THIS_SOLVER ": " RLAB_ERROR_ONE_OR_TWO_ARG_REQUIRED);

  e1 = bltin_get_ent (args[0]);
  if (ent_type(e1) != MATRIX_DENSE_REAL)
    rerror (THIS_SOLVER ": " RLAB_ERROR_ARG1_MDR_INTEGER_MATRIX);
  a = class_int(e1);

  if (nargs==2)
  {
    e2 = bltin_get_ent (args[1]);
    if (ent_type(e2) != MATRIX_DENSE_REAL)
      rerror (THIS_SOLVER ": " RLAB_ERROR_ARG1_MDR_INTEGER_MATRIX);
    use_lsb = class_int(e2);
  }
  else
    use_lsb = 1;

  rval = ent_Assign_Rlab_MDR(mdi_byte_split(a,use_lsb));

  ent_Clean(e1);
  ent_Clean(e2);

  return rval;
}

#undef  THIS_SOLVER
#define THIS_SOLVER "bitsplit"
Ent * ent_bit_split (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *e3=0, *rval;
  int nbytes=1, use_lsb=1;
  MDR *a=0;

  /* Check n_args */
  if (nargs != 1 &&nargs != 2 && nargs != 3)
    rerror (THIS_SOLVER ": " RLAB_ERROR_TWO_TO_THREE_ARG_REQUIRED);

  e1 = bltin_get_ent (args[0]);
  if (ent_type(e1) != MATRIX_DENSE_REAL)
    rerror (THIS_SOLVER ": " RLAB_ERROR_ARG1_MDR_INTEGER_MATRIX);
  a = ent_data(e1);

  if (nargs>=2)
  {
    e2 = bltin_get_ent (args[1]);
    if (ent_type(e2) == MATRIX_DENSE_REAL)
      nbytes = class_int(e2);
  }

  if (nargs==3)
  {
    e3 = bltin_get_ent (args[2]);
    if (ent_type(e3) == MATRIX_DENSE_REAL)
      use_lsb = class_int(e3);
  }

  rval = ent_Assign_Rlab_MDR(mdi_bit_split(a,nbytes,use_lsb));

  ent_Clean(e1);
  ent_Clean(e2);
  ent_Clean(e3);

  return rval;
}

#undef  THIS_SOLVER
#define THIS_SOLVER "bytejoin"
Ent * ent_byte_join (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *e3=0, *rval;
  MDR *a1=0;

  int stride=2, use_lsb=0;

  /* Check n_args */
  if ((nargs!=1) && (nargs!=2) && (nargs!=3))
    rerror (THIS_SOLVER ": " RLAB_ERROR_ONE_TO_THREE_ARG_REQUIRED);

  e1 = bltin_get_ent (args[0]);
  if (ent_type(e1) != MATRIX_DENSE_REAL)
    rerror (THIS_SOLVER ": " RLAB_ERROR_ONE_TO_THREE_ARG_REQUIRED);
  a1 = ent_data(e1);
  if (SIZE(a1)<1)
    rerror (THIS_SOLVER ": " RLAB_ERROR_ONE_TO_THREE_ARG_REQUIRED);
  if (!EQVECT(a1))
    rerror (THIS_SOLVER ": " RLAB_ERROR_ONE_TO_THREE_ARG_REQUIRED);

  if (nargs > 1)
  {
    e2 = bltin_get_ent (args[1]);
    if (ent_type(e2) != MATRIX_DENSE_REAL)
      rerror (THIS_SOLVER ": " RLAB_ERROR_ARG2_MDR_INTEGER);
    if ((stride < 5)&&(stride>0))
    {
      stride = class_int(e2);
    }
    else
    {
      printf (THIS_SOLVER ": parameter 'stride' has to be between 1 and 4!\n");
      rerror (THIS_SOLVER ": " RLAB_ERROR_ARG2_MDR_INTEGER);
    }
  }

  if (nargs==3)
  {
    e3 = bltin_get_ent (args[2]);
    if (ent_type(e3) != MATRIX_DENSE_REAL)
      rerror (THIS_SOLVER ": " RLAB_ERROR_ARG3_MDR_INTEGER);
    use_lsb = (class_int(e3) != 0);
  }

  rval = ent_Assign_Rlab_MDR(mdi_byte_join(a1,stride,use_lsb));

  ent_Clean(e1);
  ent_Clean(e2);

  return rval;
}

#undef  THIS_SOLVER
#define THIS_SOLVER "conrec"
Ent * ent_Contour(int nargs, Datum args[])
{
  Ent *e1=0, *e2=0, *e3=0, *e4=0, *rval=0;
  MDR *x=0, *y=0, *z=0, *d=0;
  MDR **r=0;
  MDE *ret=0;
  int i;

  /* Check n_args */
  if (nargs<2 && nargs>4)
    rerror (THIS_SOLVER ": " RLAB_ERROR_TWO_TO_FOUR_ARG_REQUIRED "\n");

  e1 = bltin_get_ent (args[0]);
  if (ent_type(e1) != MATRIX_DENSE_REAL)
    rerror (THIS_SOLVER ": " RLAB_ERROR_ARG1_MDR_VECTOR "\n");
  z = ent_data(e1);
  if (!EQVECT(z))
    rerror (THIS_SOLVER ": " RLAB_ERROR_ARG1_MDR_VECTOR);

  e2 = bltin_get_ent (args[1]);
  if (ent_type(e2) != MATRIX_DENSE_REAL)
    rerror (THIS_SOLVER ": " RLAB_ERROR_ARG2_MDR_MATRIX "\n");
  d = ent_data(e2);

  if (nargs > 2)
  {
    e3 = bltin_get_ent (args[2]);
    if (ent_type(e3) == MATRIX_DENSE_REAL)
    {
      x = ent_data(e3);
      if (!EQVECT(x))
        x=0;
      else if (SIZE(x)!=MNC(d))
      {
        rerror(THIS_SOLVER ": Dimension mismatch SIZE(x) != MNC(d)!\n");
      }
    }
  }

  if (nargs > 3)
  {
    e4 = bltin_get_ent (args[3]);
    if (ent_type(e4) == MATRIX_DENSE_REAL)
    {
      y = ent_data(e4);
      if (!EQVECT(y))
        y=0;
      else if (SIZE(y)!=MNR(d))
      {
        rerror(THIS_SOLVER ": Dimension mismatch SIZE(y) != MNR(d)!\n");
      }
    }
  }

  r = GC_MALLOC(SIZE(z) * sizeof(MDR));
  for (i=0; i<SIZE(z); i++)
  {
    r[i] = mdr_Create(20,2);
  }
  mdr_Contour(d, x, y, z, r);

  rval = ent_Create();
  ret = mde_Create(1,SIZE(z));
  for (i=0; i<SIZE(z);i++)
  {
    MdeV0(ret,i) = ent_Assign_Rlab_MDR(r[i]);
    r[i] = 0;
  }
  ent_data(rval) = ret;
  ent_type(rval) = MATRIX_DENSE_ENTITY;

  GC_FREE(r);
  ent_Clean(e1);
  ent_Clean(e2);
  ent_Clean(e3);
  ent_Clean(e4);
  return rval;
}

// bracketing of zero and force result to range
extern double rlab_op_zero_abs; // |res| <= zabs => res = 0
extern double rlab_op_zero_rel; // |res| <= zrel*max(arg1,arg2) => res=0
extern double rlab_op_min;      // arg1 = max(arg1, op_min)
extern double rlab_op_max;      // arg1 = min(arg1, op_max)

#undef  THIS_SOLVER
#define THIS_SOLVER "sys.op.max"
Ent * ent_Sys_Op_Max (int nargs, Datum args[])
{
  if (nargs==1)
  {
    Ent *e1 = bltin_get_ent (args[0]);
    if (ent_type(e1) == MATRIX_DENSE_REAL)
    {
      double dummy = class_double(e1);
      rlab_op_max = dummy;
    }
    ent_Clean(e1);
  }
  return ent_Create_Rlab_Double(rlab_op_max);
}

#undef  THIS_SOLVER
#define THIS_SOLVER "sys.op.min"
Ent * ent_Sys_Op_Min (int nargs, Datum args[])
{
  if (nargs==1)
  {
    Ent *e1 = bltin_get_ent (args[0]);
    if (ent_type(e1) == MATRIX_DENSE_REAL)
    {
      double dummy = class_double(e1);
      rlab_op_min = dummy;
    }
    ent_Clean(e1);
  }
  return ent_Create_Rlab_Double(rlab_op_min);
}

#undef  THIS_SOLVER
#define THIS_SOLVER "sys.op.zero.abs"
Ent * ent_Sys_Op_ZeroAbs (int nargs, Datum args[])
{
  if (nargs==1)
  {
    Ent *e1 = bltin_get_ent (args[0]);
    if (ent_type(e1) == MATRIX_DENSE_REAL)
    {
      double dummy = class_double(e1);
      if (dummy > 0)
        rlab_op_zero_abs = dummy;
    }
    ent_Clean(e1);
  }
  return ent_Create_Rlab_Double(rlab_op_zero_abs);
}

#undef  THIS_SOLVER
#define THIS_SOLVER "sys.op.zero.rel"
Ent * ent_Sys_Op_ZeroRel (int nargs, Datum args[])
{
  if (nargs==1)
  {
    Ent *e1 = bltin_get_ent (args[0]);
    if (ent_type(e1) == MATRIX_DENSE_REAL)
    {
      double dummy = class_double(e1);
      if (dummy > 0)
        rlab_op_zero_rel = dummy;
    }
    ent_Clean(e1);
  }
  return ent_Create_Rlab_Double(rlab_op_zero_rel);
}

#undef  THIS_SOLVER
#define THIS_SOLVER "isvector"
Ent * ent_IsVec (int nargs, Datum args[])
{
  double rval=0.0;
  int nr=0, nc=0;

  if (nargs!=1)
  {
    rerror(THIS_SOLVER ": " RLAB_ERROR_ARG1_MD "\n");
  }

  Ent *e1 = bltin_get_ent (args[0]);

  if (  (ent_type(e1) == MATRIX_DENSE_REAL)
         || (ent_type(e1) == MATRIX_DENSE_COMPLEX)
         || (ent_type(e1) == MATRIX_DENSE_STRING) )
  {
    nr = MNR(ent_data(e1));
    nc = MNC(ent_data(e1));
  }

  if (((nr == 1) && (nc > 0)) || ((nr > 0) && (nc == 1)))
    rval = 1.0;

  ent_Clean(e1);
  return ent_Create_Rlab_Double(rval);
}

#undef  THIS_SOLVER
#define THIS_SOLVER "ismatrix"
Ent * ent_IsMat(int nargs, Datum args[])
{
  double rval=0.0;
  int nr=0, nc=0;

  if (nargs!=1)
  {
    rerror(THIS_SOLVER ": " RLAB_ERROR_ARG1_MD "\n");
  }

  Ent *e1 = bltin_get_ent (args[0]);

  if (  (ent_type(e1) == MATRIX_DENSE_REAL)
         || (ent_type(e1) == MATRIX_DENSE_COMPLEX)
         || (ent_type(e1) == MATRIX_DENSE_STRING) )
  {
    nr = MNR(ent_data(e1));
    nc = MNC(ent_data(e1));
  }

  if ((nr > 0) && (nc > 0))
    rval = 1.0;

  ent_Clean(e1);
  return ent_Create_Rlab_Double(rval);
}

#undef  THIS_SOLVER
#define THIS_SOLVER "isscalar"
Ent * ent_IsScalar(int nargs, Datum args[])
{
  double rval=0.0;
  int nr=0, nc=0;

  if (nargs!=1)
  {
    rerror(THIS_SOLVER ": " RLAB_ERROR_ARG1_MD "\n");
  }

  Ent *e1 = bltin_get_ent (args[0]);

  if (  (ent_type(e1) == MATRIX_DENSE_REAL)
         || (ent_type(e1) == MATRIX_DENSE_COMPLEX)
         || (ent_type(e1) == MATRIX_DENSE_STRING) )
  {
    nr = MNR(ent_data(e1));
    nc = MNC(ent_data(e1));
  }

  if ((nr == 1) && (nc==1))
    rval = 1.0;

  ent_Clean(e1);
  return ent_Create_Rlab_Double(rval);
}

#undef  THIS_SOLVER
#define THIS_SOLVER "isnumber"
Ent * ent_IsNumber(int nargs, Datum args[])
{
  double rval=0.0;
  int nr=0, nc=0;

  if (nargs!=1)
  {
    rerror(THIS_SOLVER ": " RLAB_ERROR_ARG1_MD "\n");
  }

  Ent *e1 = bltin_get_ent (args[0]);

  if (  (ent_type(e1) == MATRIX_DENSE_REAL)
         || (ent_type(e1) == MATRIX_DENSE_COMPLEX)
         || (ent_type(e1) == MATRIX_SPARSE_REAL)
         || (ent_type(e1) == MATRIX_SPARSE_COMPLEX) )
  {
    nr = MNR(ent_data(e1));
    nc = MNC(ent_data(e1));

    if ((nr > 0) && (nc>0))
      rval = 1.0;
  }
  else if (ent_type(e1)==MATRIX_DENSE_STRING)
  {
    char *c = class_char_pointer(e1);
    char *p = 0;
    strtod(c, &p);
    if (c != p)
      rval = 1.0;
  }


  ent_Clean(e1);
  return ent_Create_Rlab_Double(rval);
}

#undef  THIS_SOLVER
#define THIS_SOLVER "isstring"
Ent * ent_IsString (int nargs, Datum args[])
{
  double rval=0.0;
  int nr=0, nc=0;

  if (nargs!=1)
  {
    rerror(THIS_SOLVER ": " RLAB_ERROR_ARG1_MD "\n");
  }

  Ent *e1 = bltin_get_ent (args[0]);

  if ( (ent_type(e1) == MATRIX_DENSE_STRING) )
  {
    nr = MNR(ent_data(e1));
    nc = MNC(ent_data(e1));
    if (nr * nc > 0)
      rval = 1.0;
  }

  ent_Clean(e1);
  return ent_Create_Rlab_Double(rval);
}


extern char * curr_file_name; /* main.c */
#undef  THIS_SOLVER
#define THIS_SOLVER "curr_file"
Ent * ent_Filename (int nargs, Datum args[])
{
  char *fn = 0;
  if (isvalidstring(curr_file_name)>0)
    fn = curr_file_name;
  return ent_Assign_Rlab_String (fn);
}

#undef  THIS_SOLVER
#define THIS_SOLVER "isempty"
Ent * ent_IsEmpty(int nargs, Datum args[])
{
  double rval=1.0;

  if (nargs!=1)
  {
    rerror(THIS_SOLVER ": " RLAB_ERROR_ARG1_MD "\n");
  }

  Ent *e1 = bltin_get_ent (args[0]);

  if (ent_type(e1)==UNDEF)
    goto _empty_exit;

  if (ent_type(e1)==BTREE)
  {
    rval = btree_GetRealNumNodes(ent_data(e1)) < 1;
  }
  else
  {
    rval = SIZE(ent_data(e1)) < 1;
  }

_empty_exit:

  ent_Clean(e1);

  return ent_Create_Rlab_Double(rval);
}
