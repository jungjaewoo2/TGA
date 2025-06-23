// Copyright (C) 2003-2013 Marijan Kostrun
//   part of rlab+4linux project on rlabplus.sourceforge.net
//
// loess wrapper for rlabplus
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
//
// See the file ./COPYING
// **********************************************************************


// rlab headers, located in variable $RLAB_SDK
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

// standard libraries
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <time.h>

//
// these functions are in flibs/loess/src/misc.c
//
extern double pf (double, double, double);
extern double qt (double, double);
extern int loess_(double  *, double  *, int *, double  *, double  *, int *,
                  int *, int *,
                  int *, char *, char *, double  *, char *, int *,
                  double  *, double  *, double  *,double  *, double  *, double  *,
                  double  *, double  *, double  *, double  *, double  *,
                  int *, int *, double  *, double  *, double  *);

extern int loess_pred_(double  *, double  *, double  *, int *, double  *, double  *,
                       double  *, double  *, int *,
                       int *, int *, int *, char *, double  *, char *,
                       int *, int *, double  *, double  *, double  *, double  *,
                       int *, double  *, double  *);

Ent *
ent_loess_init (int nargs, Datum args[])
{
  //
  // create a loess structure with default parameters
  //

  Ent *e1=0, *e2=0;
  MDR *y=0, *x=0, *wdata=0;
  int i, n, p;

  Btree *bw, *in, *model, *control;

  //
  // Load and Check arguments.
  //
  if (nargs != 2 )
  {
    fprintf (stdout,
             "loess.setup: Organizes data in a loess  compatible list.  Format:\n");
    fprintf (stdout,
             "loess.setup:   r=loess.setup(Y,X),\n");
    fprintf (stdout,
             "loess.setup: where  r=<<model;input;control>>, each a sublist, as well.\n");
    fprintf (stdout,
             "loess.setup: See the manual for more details.\n");
    rerror
        ("loess.setup: requires two arguments");
  }

  //
  // Get observations: y, or <<y;we>>
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) == MATRIX_DENSE_REAL)
  {
    // y
    MDR * y0 = class_matrix_real (e1);
    y = mdr_Float_BF (y0);
  }
  else if (ent_type (e1) == BTREE)
  {
    // <<val;wgt>>
    ListNode *node;
    node = btree_FindNode (ent_data (e1), "val");
    if (node != 0)
      y = class_matrix_real (var_ent (node));
    node = btree_FindNode (ent_data (e1), "wgt");
    if (node != 0)
    {
      wdata = mdr_Copy( class_matrix_real (var_ent (node)) );
      if (MNR (wdata) * MNC (wdata) != MNR (y) * MNC (y))
        rerror ("loess.setup: 'val' and 'wgt' must have the same size!");
      for (i = 0; i < MNR (wdata) * MNC (wdata); i++)
        if (MdrV0 (wdata, i) < 0)
          rerror ("loess.setup: 'wgt' must be positive !");
    }
  }
  else
    rerror ("loess.setup: 'Y' has to be a real matrix or a list <<val;wgt>>!");
  if (MNC (y) != 1 && MNR(y)!=1)
    rerror ("loess.setup: dependent variable 'y' has to be a real vector!");
  if (!wdata)
  {
    wdata = mdr_Create(MNR(y), 1);
    for (i=0; i < MNR(y); i++)
      MdrV0(wdata,i) = 1.0;
  }

  //
  // Get predictor variable: x
  //
  e2 = bltin_get_ent (args[1]);
  if (ent_type (e2) != MATRIX_DENSE_REAL)
    rerror ("loess.setup: 'X' has to be a real matrix!");
  MDR * x0 = class_matrix_real (e2);
  x = mdr_Float_BF (x0);
  n = x->nrow;
  p = x->ncol;

  // check x and y
  if (n != (y->nrow) * (y->ncol))
    rerror ("loess.setup: 'Y' and 'X' must have the same number of rows!");
  if (n * p == 0 || (y->nrow) * (y->ncol) == 0)
    rerror ("loess.setup: Cannot continue. 'Y' or 'X' is of zero-size!");

  //
  // in: x,y,weights
  //
  in = btree_Create ();
  // y
  install (in, "y", ent_Assign_Rlab_MDR(y));
  // x
  install (in, "x", ent_Assign_Rlab_MDR(x));
  // wgt
  install (in, "weights", ent_Assign_Rlab_MDR(wdata));

  //
  // model: span, deg, norm, par, drop2, family
  //
  model = btree_Create ();
  // span
  install (model, "span", ent_Create_Rlab_Double(0.75));
  // deg
  install (model, "degree", ent_Create_Rlab_Double( 2.0 ));
  // normalize
  install (model, "normalize", ent_Create_Rlab_Double( 1.0 ));
  // parameter
  MDR *par = mdr_Create(8,1);
  for(i = 0; i < 8; i++) MdrV0(par,i) = 0;
  install (model, "parametric", ent_Assign_Rlab_MDR(par));
  // drop_square
  MDR *ds = mdr_Create(8,1);
  for(i=0; i < 8; i++) MdrV0(ds,i) = 0;
  install (model, "drop_squares",ent_Assign_Rlab_MDR(ds));
  // family
  install (model, "family", ent_Create_Rlab_String("gaussian"));


  control = btree_Create ();
  // surface
  install (control, "surface", ent_Create_Rlab_String("interpolate"));
  // statistics
  install (control, "statistics", ent_Create_Rlab_String("approximate"));
  // trace_hat
  install (control, "trace_hat", ent_Create_Rlab_String("wait.to.decide"));
  // cell
  install (control, "cell", ent_Create_Rlab_Double( 0.2 ));
  // iterations
  install (control, "iterations", ent_Create_Rlab_Double( 4.0 ));


  // clean-up
  ent_Clean (e1);
  ent_Clean (e2);

  //
  // write result to output list
  //
  bw = btree_Create ();
  // in
  install (bw, "input", ent_Assign_Rlab_BTREE(in));
  // model
  install (bw, "model", ent_Assign_Rlab_BTREE(model));
  // control
  install (bw, "control", ent_Assign_Rlab_BTREE(control));
  //
  return ent_Assign_Rlab_BTREE(bw);
}

Ent *
ent_loess (int nargs, Datum args[])
{
  //
  // calculates
  //

  Ent *e1=0;

  Btree *bw, *out, *kdt;

  ListNode * node1, * node2;

  //
  // Load and Check arguments.
  //
  if (nargs != 1 )
  {
    fprintf (stdout,
             "loess.main: Computes the elements of loess compatible list.  Format:\n");
    fprintf (stdout,
             "loess.main:   r=loess.main(lo),\n");
    fprintf (stdout,
             "loess.main: where  r=<<input;model;control;kd_tree;output>>, each a list itself.\n");
    fprintf (stdout,
             "loess.main: See the manual for more information.\n");
    rerror
        ("loess.main: requires a single argument");
  }

  //
  // Get observations: y, or <<y;we>>
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) != BTREE)
    rerror("loess.main: 'lo' is a list created by previous call to loess.setup!");

  //
  // input
  //
  node1 = btree_FindNode (ent_data (e1), "input");
  if (node1 == 0)
    rerror("loess.main: lo=<<input;model;control>>, Missing 'input' sub-list!");
  // x:
  node2 = btree_FindNode (ent_data(var_ent (node1)), "x");
  if (node2 == 0)
    rerror("loess.main: lo=<<input;model;control>>, Missing 'input.x' explanatory data!");
  MDR *x1 = class_matrix_real (var_ent (node2));
  MDR *x  = mdr_Transpose (x1);
  // y:
  node2 = btree_FindNode (ent_data(var_ent (node1)), "y");
  if (node2 == 0)
    rerror("loess.main: lo=<<input;model;control>>, Missing 'input.y' response data!");
  MDR *y0 = class_matrix_real (var_ent (node2));
  MDR *y = mdr_Float_BF (y0);
  // weights:
  node2 = btree_FindNode (ent_data(var_ent (node1)), "weights");
  if (node2 == 0)
    rerror("loess.main: lo=<<input;model;control>>, Missing 'input.weights' weights!");
  MDR * weights0 = class_matrix_real (var_ent (node2));
  MDR * weights  = mdr_Float_BF (weights0);
  // n,p
  int n = MNR (x1);
  int p = MNC (x1);

  //
  // model
  //
  node1 = btree_FindNode (ent_data (e1), "model");
  if (node1 == 0)
    rerror("loess.main: lo=<<input;model;control>>, Missing 'model' sub-list!");
  // span:
  node2 = btree_FindNode (ent_data(var_ent (node1)), "span");
  if (node2 == 0)
    rerror("loess.main: lo=<<input;model;control>>, Missing 'model.span' parameter!");
  double span = class_double (var_ent (node2));
  // degree:
  node2 = btree_FindNode (ent_data(var_ent (node1)), "degree");
  if (node2 == 0)
    rerror("loess.main: lo=<<input;model;control>>, Missing 'model.degree' parameter!");
  int degree = (int) class_double (var_ent (node2));
  // normalize
  node2 = btree_FindNode (ent_data(var_ent (node1)), "normalize");
  if (node2 == 0)
    rerror("loess.main: lo=<<input;model;control>>, Missing 'model.normalize' parameter!");
  int normalize = (int) class_double (var_ent (node2));
  // parametric
  node2 = btree_FindNode (ent_data(var_ent (node1)), "parametric");
  if (node2 == 0)
    rerror("loess.main: lo=<<input;model;control>>, Missing 'model.parametric' parameter!");
  MDR * par  = class_matrix_real (var_ent (node2));
  MDR * ipar = mdr_Int_BF (par);
  // drop_squares
  node2 = btree_FindNode (ent_data(var_ent (node1)), "drop_squares");
  if (node2 == 0)
    rerror("loess.main: lo=<<input;model;control>>, Missing 'model.parametric' parameter!");
  MDR * ds  = class_matrix_real (var_ent (node2));
  MDR * ids = mdr_Int_BF (ds);
  // family
  node2 = btree_FindNode (ent_data(var_ent (node1)), "family");
  if (node2 == 0)
    rerror("loess.main: lo=<<input;model;control>>, Missing 'model.family' parameter!");
  MDS * fam = ent_data (var_ent (node2));

  //
  // control
  //
  node1 = btree_FindNode (ent_data (e1), "control");
  // surface
  node2 = btree_FindNode (ent_data(var_ent (node1)), "surface");
  if (node2 == 0)
    rerror("loess.main: lo=<<input;model;control>>, Missing 'control.surface' parameter!");
  MDS * surf = ent_data (var_ent (node2));
  // statistics
  node2 = btree_FindNode (ent_data(var_ent (node1)), "statistics");
  if (node2 == 0)
    rerror("loess.main: lo=<<input;model;control>>, Missing 'control.statistics' parameter!");
  MDS * stat = ent_data (var_ent (node2));
  // cell
  node2 = btree_FindNode (ent_data(var_ent (node1)), "cell");
  if (node2 == 0)
    rerror("loess.main: lo=<<input;model;control>>, Missing 'control.cell' parameter!");
  double cell = class_double (var_ent (node2));
  // trace_hat
  node2 = btree_FindNode (ent_data(var_ent (node1)), "trace_hat");
  if (node2 == 0)
    rerror("loess.main: lo=<<input;model;control>>, Missing 'control.trace_hat' parameter!");
  MDS * th = ent_data (var_ent (node2));
  // iterations
  node2 = btree_FindNode (ent_data(var_ent (node1)), "iterations");
  if (node2 == 0)
    rerror("loess.main: lo=<<input;model;control>>, Missing 'control.iterations' parameter!");
  int iterations = (int) class_double (var_ent (node2));

  //
  // output
  //
  out = btree_Create ();

  // fitval
  MDR * out_fitval = mdr_Create(n,1);
  install (out, "fitted_values", ent_Assign_Rlab_MDR(out_fitval));
  // fitted_residuals
  MDR * out_fitres = mdr_Create(n,1);
  install (out, "fitted_residuals", ent_Assign_Rlab_MDR(out_fitres));
  // pseudovalues
  MDR * out_pvals = mdr_Create (n,1);
  install (out, "pseudovalues", ent_Assign_Rlab_MDR(out_pvals));
  // diagonal
  MDR * out_diag = mdr_Create (n,1);
  install (out, "diagonal", ent_Assign_Rlab_MDR(out_diag));
  // robust
  MDR * out_rob = mdr_Create (n,1);
  install (out, "robust", ent_Assign_Rlab_MDR(out_rob));
  // divisor
  MDR * out_div = mdr_Create (n,1);
  install (out, "divisor", ent_Assign_Rlab_MDR(out_div));

  int max_kd = n > 200 ? n : 200;

  //
  // kdtree
  //
  kdt = btree_Create ();
  // xi
  MDR * kdt_xi = mdr_Create (max_kd,1);
  install (kdt, "xi", ent_Assign_Rlab_MDR(kdt_xi));
  // vert
  MDR * kdt_vert = mdr_Create (2*p, 1);
  install (kdt, "vert", ent_Assign_Rlab_MDR(kdt_vert));
  // vval
  MDR * kdt_vval = mdr_Create ((p+1)*max_kd, 1);
  install (kdt, "vval", ent_Assign_Rlab_MDR(kdt_vval));

  //
  // calculate kdtree and output using loess_
  //
  int  size_info[2];

  MDR * kdt_ipar = mdi_Create(7,1);
  MDR * kdt_a = mdi_Create (max_kd,1);
  //
  size_info[0] = p;
  size_info[1] = n;
  //
  int it2 = (!strcmp(MdsV0(fam,0), "gaussian")) ? 0 : iterations;
  if(!strcmp(mdstring(th) , "wait.to.decide"))
  {
    if(!strcmp(mdstring(surf), "interpolate"))
    {
      if (n < 0)
        mdstring(th)= cpstr("exact");
      else
        mdstring(th) = cpstr("approximate");
    }
    else
    {
      mdstring(th) = cpstr("exact");
    }
  }
  //
  double enp=0, s=0, oned=0, twod=0, tht=0;

  loess_(MDRPTR(y), MDRPTR(x), size_info, MDRPTR(weights),
         &span, &degree, MDIPTR(ipar), MDIPTR(ids), &normalize,
          mdstring(stat), mdstring(surf), &cell, mdstring(th),
         &it2,
         MDRPTR(out_fitval), MDRPTR(out_fitres),
        &enp, &s, &oned, &twod, MDRPTR(out_pvals), &tht,
         MDRPTR(out_diag), MDRPTR(out_rob), MDRPTR(out_div),
         MDIPTR(kdt_ipar), MDIPTR(kdt_a), MDRPTR(kdt_xi), MDRPTR(kdt_vert), MDRPTR(kdt_vval));

  //
  // insert the remaining elements into the list output and kdt
  //
  // output:
  // s
  install (out, "s", ent_Create_Rlab_Double(s));
  // one_delta
  install (out, "one_delta", ent_Create_Rlab_Double(oned));
  // two_delta
  install (out, "two_delta", ent_Create_Rlab_Double(twod));
  // trace_hat
  install (out, "trace_hat", ent_Create_Rlab_Double(tht));
  // enp
  install (out, "enp", ent_Create_Rlab_Double(enp));

  // kd_tree:
  //
  // parameter
  install (kdt, "parameter", ent_Assign_Rlab_MDR(kdt_ipar));
  // a
  install (kdt, "a", ent_Assign_Rlab_MDR(kdt_a));

  //
  // create a new, now complete, loess list
  //
  bw = btree_Copy ( ent_data (e1) );
  // kd_tree
  install (bw, "kd_tree", ent_Assign_Rlab_BTREE(kdt));
  // output
  install (bw, "output", ent_Assign_Rlab_BTREE(out));

  // clean-up
  ent_Clean (e1);
  mdr_Destroy (ipar);
  mdr_Destroy (ids);
  mdr_Destroy (x); mdr_Destroy (y); mdr_Destroy (weights);

  //
  return ent_Assign_Rlab_BTREE (bw);
}

Ent *
ent_loess_predict (int nargs, Datum args[])
{
  //
  // calculates
  //

  Ent *ex=0, *e1=0, *e3=0;

  int m;
  double  coverage=0.0;

  Btree *bw;

  ListNode * node1, * node2;

  //
  // Load and Check arguments.
  //
  if (nargs != 2 && nargs != 3)
  {
    fprintf (stdout,
             "loess.predict: Form a prediction based on the previously calculated  loess  structure.\n");
    fprintf (stdout,
             "loess.predict: Format:\n");
    fprintf (stdout,
             "loess.predict:   r=loess.predict(x,lo/,coverage/),\n");
    fprintf (stdout,
             "loess.predict: where  'lo'  is a loess compatible list, 'x' is a data set for which the\n");
    fprintf (stdout,
             "loess.predict: prediction is needed,  while 0 < coverage < 1 is the confidence interval\n");
    fprintf (stdout,
             "loess.predict: for prediction, if requested.  The result 'r' is a loess compatible\n");
    fprintf (stdout,
             "loess.predict: 'prediction' list, r=<<fit;se_fit;df;residual_scale/;ci_lower;ci_upper/>>,\n");
    fprintf (stdout,
             "loess.predict: where entries 'ci_..' are calculated only if coverage>0. See the manual\n");
    fprintf (stdout,
             "loess.predict: for more information.\n");
    rerror
        ("loess.predict: two or three arguments required.");
  }

  //
  // x: data points for predictions
  //
  ex = bltin_get_ent (args[0]);
  if (ent_type (ex) != MATRIX_DENSE_REAL)
    rerror("loess.predict: 'x' is a real matrix!");
  MDR * x0 = class_matrix_real (ex);
  MDR * xx = mdr_Transpose (x0);

  m = MNR (x0);

  //
  // loess list
  //
  e1 = bltin_get_ent (args[1]);
  if (ent_type (e1) != BTREE)
    rerror("loess.predict: 'lo' is a list created by previous call to loess_init!");

  //
  // input
  //
  node1 = btree_FindNode (ent_data (e1), "input");
  if (node1 == 0)
    rerror("loess.predict: lo=<<input;model;control>>, Missing 'input' list!");
  // x:
  node2 = btree_FindNode (ent_data(var_ent (node1)), "x");
  if (node2 == 0)
    rerror("loess.predict: lo=<<input;model;control>>, Missing 'input.x' explanatory data!");
  MDR *x1 = class_matrix_real (var_ent (node2));
  MDR *x  = mdr_Transpose (x1);
  // y:
  node2 = btree_FindNode (ent_data(var_ent (node1)), "y");
  if (node2 == 0)
    rerror("loess.predict: lo=<<input;model;control>>, Missing 'input.y' response data!");
  MDR *y = class_matrix_real (var_ent (node2));
  // weights:
  node2 = btree_FindNode (ent_data(var_ent (node1)), "weights");
  if (node2 == 0)
    rerror("loess.predict: lo=<<input;model;control>>, Missing 'input.weights' weights!");
  MDR * weights = class_matrix_real (var_ent (node2));
  // n,p
  int n = MNR (x1);
  int p = MNC (x1);

  //
  // model
  //
  node1 = btree_FindNode (ent_data (e1), "model");
  if (node1 == 0)
    rerror("loess.predict: lo=<<input;model;control>>, Missing 'model' list!");
  // span:
  node2 = btree_FindNode (ent_data(var_ent (node1)), "span");
  if (node2 == 0)
    rerror("loess.predict: lo=<<input;model;control>>, Missing 'model.span' parameter!");
  double mod_span = class_double (var_ent (node2));
  // degree:
  node2 = btree_FindNode (ent_data(var_ent (node1)), "degree");
  if (node2 == 0)
    rerror("loess.predict: lo=<<input;model;control>>, Missing 'model.degree' parameter!");
  int mod_degree = (int) class_double (var_ent (node2));
  // normalize
  node2 = btree_FindNode (ent_data(var_ent (node1)), "normalize");
  if (node2 == 0)
    rerror("loess.predict: lo=<<input;model;control>>, Missing 'model.normalize' parameter!");
  int mod_normalize = (int) class_double (var_ent (node2));
  // parametric
  node2 = btree_FindNode (ent_data(var_ent (node1)), "parametric");
  if (node2 == 0)
    rerror("loess.predict: lo=<<input;model;control>>, Missing 'model.parametric' parameter!");
  MDR * mod_par  = ent_data(var_ent (node2));
  MDR * mod_ipar = mdr_Int_BF (mod_par);
  // drop_squares
  node2 = btree_FindNode (ent_data(var_ent (node1)), "drop_squares");
  if (node2 == 0)
    rerror("loess.predict: lo=<<input;model;control>>, Missing 'model.parametric' parameter!");
  MDR * mod_ds  = ent_data (var_ent (node2));
  MDR * mod_ids = mdr_Int_BF (mod_ds);
  // family
  node2 = btree_FindNode (ent_data(var_ent (node1)), "family");
  if (node2 == 0)
    rerror("loess.predict: lo=<<input;model;control>>, Missing 'model.family' parameter!");
  MDS * mod_fam = ent_data (var_ent (node2));


  //
  // control
  //
  node1 = btree_FindNode (ent_data (e1), "control");
  if (node1 == 0)
    rerror("loess.predict: lo=<<input;model;control;output;kd_tree>>, missing 'control' list!");
  // surface
  node2 = btree_FindNode (ent_data(var_ent (node1)), "surface");
  if (node2 == 0)
    rerror("loess.predict: lo=<<input;model;control>>, Missing 'control.surface' parameter!");
  MDS * con_surf = ent_data (var_ent (node2));
  // statistics
//   node2 = btree_FindNode (ent_data(var_ent (node1)), "statistics");
//   if (node2 == 0)
//     rerror("loess.predict: lo=<<input;model;control>>, Missing 'control.statistics' parameter!");
//   MDS * con_stat = ent_data (var_ent (node2));
  // cell
  node2 = btree_FindNode (ent_data(var_ent (node1)), "cell");
  if (node2 == 0)
    rerror("loess.predict: lo=<<input;model;control>>, Missing 'control.cell' parameter!");
  double con_cell = class_double (var_ent (node2));
  // trace_hat
//   node2 = btree_FindNode (ent_data(var_ent (node1)), "trace_hat");
//   if (node2 == 0)
//     rerror("loess.predict: lo=<<input;model;control>>, Missing 'control.trace_hat' parameter!");
//   MDS * con_th = ent_data (var_ent (node2));
  // iterations
//   node2 = btree_FindNode (ent_data(var_ent (node1)), "iterations");
//   if (node2 == 0)
//     rerror("loess.predict: lo=<<input;model;control>>, Missing 'control.iterations' parameter!");
//   int con_iterations = (int) class_double (var_ent (node2));

  //
  // output
  //
  node1 = btree_FindNode (ent_data (e1), "output");
  if (node1 == 0)
    rerror("loess.predict: lo=<<input;model;control;output;kd_tree>>, missing 'output' list!");
  // divisor
  node2 = btree_FindNode (ent_data(var_ent (node1)), "divisor");
  if (node2 == 0)
    rerror("loess.predict: lo=<<input;model;control;output;kd_tree>>, Missing 'output.divisor' entry!");
  MDR * out_div = class_matrix_real (var_ent (node2));
  // s
  node2 = btree_FindNode (ent_data(var_ent (node1)), "s");
  if (node2 == 0)
    rerror("loess.predict: lo=<<input;model;control;output;kd_tree>>, Missing 'output.s' entry!");
  double out_s = class_double (var_ent (node2));
  // one_delta
  node2 = btree_FindNode (ent_data(var_ent (node1)), "one_delta");
  if (node2 == 0)
    rerror("loess.predict: lo=<<input;model;control;output;kd_tree>>, Missing 'output.one_delta' entry!");
  double out_1d = class_double (var_ent (node2));
  // two_delta
  node2 = btree_FindNode (ent_data(var_ent (node1)), "two_delta");
  if (node2 == 0)
    rerror("loess.predict: lo=<<input;model;control;output;kd_tree>>, Missing 'output.two_delta' entry!");
  double out_2d = class_double (var_ent (node2));
  // robust
  node2 = btree_FindNode (ent_data(var_ent (node1)), "robust");
  if (node2 == 0)
    rerror("loess.predict: lo=<<input;model;control;output;kd_tree>>, Missing 'output.robust' entry!");
  MDR * out_rob = class_matrix_real (var_ent (node2));

  //
  // kdtree
  //
  node1 = btree_FindNode (ent_data (e1), "kd_tree");
  if (node1 == 0)
    rerror("loess.predict: lo=<<input;model;control;output;kd_tree>>, missing 'kd_tree' list!");
  // parameter
  node2 = btree_FindNode (ent_data(var_ent (node1)), "parameter");
  if (node2 == 0)
    rerror("loess.predict: lo=<<input;model;control;output;kd_tree>>, Missing 'kd_tree.parameter' entry!");
  MDR * kdt_par  = ent_data (var_ent (node2));
  MDR * kdt_ipar = mdr_Int_BF(kdt_par);
  // a
  node2 = btree_FindNode (ent_data(var_ent (node1)), "a");
  if (node2 == 0)
    rerror("loess.predict: lo=<<input;model;control;output;kd_tree>>, Missing 'kd_tree.a' entry!");
  MDR * kdt_ra = ent_data (var_ent (node2));
  MDR * kdt_a  = mdr_Int_BF (kdt_ra);
  // xi
  node2 = btree_FindNode (ent_data(var_ent (node1)), "xi");
  if (node2 == 0)
    rerror("loess.predict: lo=<<input;model;control;output;kd_tree>>, Missing 'kd_tree.xi' entry!");
  MDR * kdt_xi = class_matrix_real (var_ent (node2));
  // vert
  node2 = btree_FindNode (ent_data(var_ent (node1)), "vert");
  if (node2 == 0)
    rerror("loess.predict: lo=<<input;model;control;output;kd_tree>>, Missing 'kd_tree.vert' entry!");
  MDR * kdt_vert = class_matrix_real (var_ent (node2));
  // vval
  node2 = btree_FindNode (ent_data(var_ent (node1)), "vval");
  if (node2 == 0)
    rerror("loess.predict: lo=<<input;model;control;output;kd_tree>>, Missing 'kd_tree.vval' entry!");
  MDR * kdt_vval = class_matrix_real (var_ent (node2));

  //
  // coverage
  //
  if (nargs == 3)
  {
    e3 = bltin_get_ent (args[2]);
    if (ent_type (e3) == MATRIX_DENSE_REAL)
      coverage = class_double (e3);
    coverage = coverage > 0.0 ? coverage : 0.0;
    coverage = 1.0 > coverage ? coverage : 1.0;
    ent_Clean(e3);
  }

  //
  // calculate kdtree and output using loess_
  //
  int  size_info[3];
  //
  size_info[0] = p;
  size_info[1] = n;
  size_info[2] = m;

  MDR *pre_fit    = mdr_Create(m, 1);
  MDR *pre_se_fit = mdr_Create(m, 1);
  double pre_res_scale = out_s;
  double pre_df = out_1d * out_1d / out_2d;
  int se = 1; // calculate standard errors as well

  //
  // call loess_pred_
  //
  loess_pred_(MDRPTR(y), MDRPTR(x), MDRPTR(xx), size_info,
              &out_s, MDRPTR(weights), MDRPTR(out_rob),
              &mod_span, &mod_degree, &mod_normalize,
              MDIPTR(mod_ipar), MDIPTR(mod_ids), mdstring(con_surf),
             &con_cell, mdstring(mod_fam),
              MDIPTR(kdt_ipar), MDIPTR(kdt_a), MDRPTR(kdt_xi), MDRPTR(kdt_vert), MDRPTR(kdt_vval),
              MDRPTR(out_div), &se, MDRPTR(pre_fit), MDRPTR(pre_se_fit));

  MDR * out_upper=0, *out_lower=0;
  if (coverage > 0.0)
  {
    //
    // calculate confidence intervals, as well (pointwise)
    //
    double limit, qt();
    int i;

    out_upper = mdr_Create(m,1);
    out_lower = mdr_Create(m,1);
    double t_dist = qt(1.0 - (1.0 - coverage)/2.0, pre_df);
    for(i = 0; i < m; i++)
    {
      limit = MdrV0(pre_se_fit,i) * t_dist;
      MdrV0(out_upper,i) = MdrV0(pre_fit,i) + limit;
      MdrV0(out_lower,i) = MdrV0(pre_fit,i) - limit;
    }
  }

  //
  // create a partial predict loess list
  //
  bw = btree_Create ();

  // fit
  install  (bw, "fit", ent_Assign_Rlab_MDR(pre_fit));
  // standard errors of the fit
  install  (bw, "se_fit", ent_Assign_Rlab_MDR(pre_se_fit));
  // residual scale
  install  (bw, "residual_scale", ent_Create_Rlab_Double(pre_res_scale));
  // degrees of freedom
  install  (bw, "df", ent_Create_Rlab_Double(pre_df));

  if (coverage > 0.0)
  {
    // confidence interval: upper limit
    install  (bw, "ci_upper", ent_Assign_Rlab_MDR(out_upper));
    // confidence interval: lower limit
    install  (bw, "ci_lower", ent_Assign_Rlab_MDR(out_lower));
  }

  // clean-up
  ent_Clean (ex);
  ent_Clean (e1);
  mdr_Destroy (mod_ipar);
  mdr_Destroy (mod_ids);
  mdr_Destroy (x);
  mdr_Destroy (kdt_ipar);
  mdr_Destroy (kdt_a);
  mdr_Destroy (xx);

  //
  return ent_Assign_Rlab_BTREE(bw);
}

Ent *
ent_loess_anova (int nargs, Datum args[])
{
  //
  // calculates
  //

  Ent *e1=0, *e2=0;

  Btree *bw;

  ListNode * node1, * node2;

  //
  // Load and Check arguments.
  //
  if (nargs != 2 )
  {
    fprintf (stdout,
             "loess.anova: Compare two loess models using analysis of variance.\n");
    fprintf (stdout,
             "loess.anova: Format:\n");
    fprintf (stdout,
             "loess.anova:   a=loess.anova(lo_1, lo_2),\n");
    fprintf (stdout,
             "loess.anova: where  'lo_1,2'  are two loess compatible lists. The result 'a' is\n");
    fprintf (stdout,
             "loess.anova: a loess compatible 'anova',  a=<<dfn;dfd;F_value;Pr_F>>.\n");
    fprintf (stdout,
             "loess.anova: See the manual for more information.\n");
    rerror
        ("loess.anova: requires two arguments");
  }

  //
  // loess list #1
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) != BTREE)
    rerror("loess.anova: 'lo' is a list created by previous call to loess_init!");

  //
  // output
  //
  node1 = btree_FindNode (ent_data (e1), "output");
  if (node1 == 0)
    rerror("loess.anova: lo=<<input;model;control;output;kd_tree>>, missing 'output' list!");
  // s
  node2 = btree_FindNode (ent_data(var_ent (node1)), "s");
  if (node2 == 0)
    rerror("loess.anova: lo=<<input;model;control;output;kd_tree>>, Missing 'output.s' entry!");
  double out1_s = class_double (var_ent (node2));
  // one_delta
  node2 = btree_FindNode (ent_data(var_ent (node1)), "one_delta");
  if (node2 == 0)
    rerror("loess.anova: lo=<<input;model;control;output;kd_tree>>, Missing 'output.one_delta' entry!");
  double out1_d1 = class_double (var_ent (node2));
  // two_delta
  node2 = btree_FindNode (ent_data(var_ent (node1)), "two_delta");
  if (node2 == 0)
    rerror("loess.anova: lo=<<input;model;control;output;kd_tree>>, Missing 'output.two_delta' entry!");
  double out1_d2 = class_double (var_ent (node2));
  // enp
  node2 = btree_FindNode (ent_data(var_ent (node1)), "enp");
  if (node2 == 0)
    rerror("loess.anova: lo=<<input;model;control;output;kd_tree>>, Missing 'output.enp' entry!");
  double out1_enp = class_double (var_ent (node2));

  //
  // loess list #2
  //
  e2 = bltin_get_ent (args[1]);
  if (ent_type (e1) != BTREE)
    rerror("loess.anova: 'lo' is a list created by previous call to loess.setup!");

  //
  // output
  //
  node1 = btree_FindNode (ent_data (e1), "output");
  if (node1 == 0)
    rerror("loess.anova: lo=<<input;model;control;output;kd_tree>>, missing 'output' list!");
  // s
  node2 = btree_FindNode (ent_data(var_ent (node1)), "s");
  if (node2 == 0)
    rerror("loess.anova: lo=<<input;model;control;output;kd_tree>>, Missing 'output.s' entry!");
  double out2_s = class_double (var_ent (node2));
  // one_delta
  node2 = btree_FindNode (ent_data(var_ent (node1)), "one_delta");
  if (node2 == 0)
    rerror("loess.anova: lo=<<input;model;control;output;kd_tree>>, Missing 'output.one_delta' entry!");
  double out2_d1 = class_double (var_ent (node2));
  // two_delta
  node2 = btree_FindNode (ent_data(var_ent (node1)), "two_delta");
  if (node2 == 0)
    rerror("loess.anova: lo=<<input;model;control;output;kd_tree>>, Missing 'output.two_delta' entry!");
  double out2_d2 = class_double (var_ent (node2));
  // enp
  node2 = btree_FindNode (ent_data(var_ent (node1)), "enp");
  if (node2 == 0)
    rerror("loess.anova: lo=<<input;model;control;output;kd_tree>>, Missing 'output.enp' entry!");
  double out2_enp = class_double (var_ent (node2));

  //
  // calculate anova
  //
  double  rssdiff, d1diff, tmp, pf();
  int     max_enp;

  rssdiff = fabs(out1_s * out1_s * out1_d1 - out2_s * out2_s * out2_d1);
  d1diff = fabs(out1_d1 - out2_d1);

  double dfn  = d1diff * d1diff / fabs(out1_d2 - out2_d2);
  max_enp = (out1_enp > out2_enp);

  tmp     = max_enp ? out1_d1 : out2_d1;
  double dfd  = tmp * tmp / (max_enp ? out1_d2 : out2_d2);

  tmp     = max_enp ? out1_s : out2_s;
  double F_value = (rssdiff / d1diff) / (tmp * tmp);
  double Pr_F    = 1 - pf(F_value, dfn, dfd);

  // clean-up
  ent_Clean (e1);
  ent_Clean (e2);

  //
  // create a partial anova loess list
  //
  bw = btree_Create ();
  // dfn
  install (bw, "dfn", ent_Create_Rlab_Double(dfn));
  // dfd
  install (bw, "dfd", ent_Create_Rlab_Double(dfd));
  // F_value
  install (bw, "F_value", ent_Create_Rlab_Double(F_value));
  // Pr_F
  install (bw, "Pr_F", ent_Create_Rlab_Double(Pr_F));
  //
  return ent_Assign_Rlab_BTREE (bw);
}

