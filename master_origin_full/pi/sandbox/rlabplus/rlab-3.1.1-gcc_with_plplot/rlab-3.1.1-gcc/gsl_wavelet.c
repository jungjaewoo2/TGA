// rlabplus (C) 2003-2007 Marijan Kostrun
//
// GSL Science Library - wavelet transforms
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
#include "mdc.h"
#include "list.h"
#include "btree.h"
#include "bltin.h"
#include "util.h"
#include "mathl.h"
#include "function.h"
#include "lp.h"


// gsl headers
#include <gsl/gsl_mode.h>
#include <gsl/gsl_precision.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_machine.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_wavelet.h>
#include <gsl/gsl_wavelet2d.h>

// standard libraries
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <time.h>

// naming convention for the solver parameters
#include "rlab_solver_parameters_names.h"


struct _wavelets
{
  int  index;
  char name[16];
};

static struct _wavelets bltin_wavelets[] =
{
  {1, "daub"},
  {2, "daub_c"},
  {3, "haar"},
  {4, "haar_c"},
  {5, "bspline"},
  {6, "bspline_c"},
  {0, "\0"}
};


Ent *
ent_gsl_dwt (int nargs, Datum args[])
{
  Ent *e1 = 0, *e2 = 0;
  ListNode *node;
  MDR *x1, *xd;
  int n, nr, nc, i, iwav=0, idx;

  gsl_wavelet *w = 0;
  gsl_wavelet_workspace * work = 0;

  if (nargs != 2)
    rerror("dwt: two arguments required");

  //
  // x
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) != MATRIX_DENSE_REAL)
    rerror("dwt: first argument has to be real vector");
  x1 = ent_data (e1);

  nr = x1->nrow;
  nc = x1->ncol;
  if (nr != 1 && nc != 1)
    rerror("dwt: first argument has to be real vector");
  n  = nr * nc;
  if (n == 0)
    rerror("dwt: first argument has to be real vector");

  n = (int) pow(2, (int) (log(n)/log(2)));

  xd = mdr_Float_BF ( x1 );

  //
  // parameters
  //
  e2 = bltin_get_ent (args[1]);
  if (ent_type (e2) == BTREE)
  {
    // family name
    node = btree_FindNode (ent_data (e2), RLAB_NAME_GSL_WV_FAMILY);
    if (node != 0)
    {
      char * fam = class_char_pointer (var_ent (node));
      if (!fam)
        rerror("dwt: incorrect entry for 'family'");
      iwav = 0;
      for (i = 0; bltin_wavelets[i].index; i++)
      {
        if (!strcmp (fam, bltin_wavelets[i].name))
        {
          iwav = bltin_wavelets[i].index;
          break;
        }
      }
      if (!iwav)
        rerror("dwt: incorrect entry for 'family'");
    }
    else
      rerror("dwt: incorrect entry for 'family'");

    // index
    node = btree_FindNode (ent_data (e2), RLAB_NAME_GSL_WV_INDEX);
    if (node != 0)
    {
      idx = 0;
      switch (iwav)
      {
        case 1:
        case 2:
          // daubechies family accepts index
          idx = (int) class_double (var_ent (node));
          if (iwav == 1)
            w = gsl_wavelet_alloc(gsl_wavelet_daubechies, idx);
          else
            w = gsl_wavelet_alloc(gsl_wavelet_daubechies_centered, idx);
          break;

        case 3:
        case 4:
          // haar family accepts index 2, ignore user's entry
          if (iwav == 3)
            w = gsl_wavelet_alloc(gsl_wavelet_haar, 2);
          else
            w = gsl_wavelet_alloc(gsl_wavelet_haar_centered, 2);
          break;

        case 5:
        case 6:
          // bspline family accepts index [i,j] or 100*i+j
          if (ent_type(var_ent (node)) != MATRIX_DENSE_REAL)
            rerror("dwt: incorrect entry for 'index' of family 'bspline/bspline_c'");
          MDR *ix = ent_data(var_ent(node));
          if (ix->type == RLAB_TYPE_INT32)
          {
            if (MNR(ix)*MNC(ix) == 1)
              idx = MdiV0(ix,0);
            else
              idx = 100 * MdiV0(ix,0) + MdiV0(ix,1);
          }
          else
          {
            if (MNR(ix)*MNC(ix) == 1)
              idx = MdrV0(ix,0);
            else
              idx = 100 * MdrV0(ix,0) + MdrV0(ix,1);
          }
          if (iwav == 5)
            w = gsl_wavelet_alloc(gsl_wavelet_bspline, idx);
          else
            w = gsl_wavelet_alloc(gsl_wavelet_bspline_centered, idx);
          break;
      }
      if (!w)
        rerror("dwt: incorrect entry for 'index' of wavelet family");

    }
    else
      rerror("dwt: missing entry for 'index' of wavelet family");
  }
  else
    rerror("dwt: missing wavelet options");

  work = gsl_wavelet_workspace_alloc (n);
  if (!work)
    rerror("dwt: terrible internal error. cannot allocate memory");

  gsl_wavelet_transform_forward(w, MDRPTR(xd), 1, n, work);

  gsl_wavelet_free (w);
  gsl_wavelet_workspace_free (work);

  ent_Clean (e1);
  ent_Clean (e2);

  return ent_Assign_Rlab_MDR(xd);
}

Ent *
ent_gsl_idwt (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0;
  ListNode *node;
  MDR *x1=0, *xd=0;
  int n, nr, nc, i, iwav=0, idx;

  gsl_wavelet *w = 0;
  gsl_wavelet_workspace * work = 0;

  if (nargs != 2)
    rerror("idwt: two arguments required");

  //
  // x
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) != MATRIX_DENSE_REAL)
    rerror("idwt: first argument has to be real vector");
  x1 = ent_data (e1);

  nr = x1->nrow;
  nc = x1->ncol;
  if (nr != 1 && nc != 1)
    rerror("idwt: first argument has to be real vector");
  n  = nr * nc;
  if (n == 0)
    rerror("idwt: first argument has to be real vector");

  n = (int) pow(2, (int) (log(n)/log(2)));

  xd = mdr_Float_BF ( x1 );

  //
  // parameters
  //
  e2 = bltin_get_ent (args[1]);
  if (ent_type (e2) == BTREE)
  {
    // family name
    node = btree_FindNode (ent_data (e2), "family");
    if (node != 0)
    {
      char * fam = class_char_pointer (var_ent (node));
      if (!fam)
        rerror("idwt: incorrect entry for 'family'");
      iwav = 0;
      for (i = 0; bltin_wavelets[i].index; i++)
      {
        if (!strcmp (fam, bltin_wavelets[i].name))
        {
          iwav = bltin_wavelets[i].index;
          break;
        }
      }
      if (!iwav)
        rerror("idwt: incorrect entry for 'family'");
    }
    else
      rerror("idwt: incorrect entry for 'family'");

    // index
    node = btree_FindNode (ent_data (e2), "index");
    if (node != 0)
    {
      idx = 0;
      switch (iwav)
      {
        case 1:
        case 2:
          // daubechies family accepts index
          idx = (int) class_double (var_ent (node));
          if (iwav == 1)
            w = gsl_wavelet_alloc(gsl_wavelet_daubechies, idx);
          else
            w = gsl_wavelet_alloc(gsl_wavelet_daubechies_centered, idx);
          break;

        case 3:
        case 4:
          // haar family accepts index 2, ignore user's entry
          if (iwav == 3)
            w = gsl_wavelet_alloc(gsl_wavelet_haar, 2);
          else
            w = gsl_wavelet_alloc(gsl_wavelet_haar_centered, 2);
          break;

        case 5:
        case 6:
          // bspline family accepts index [i,j] or 100*i+j
          if (ent_type(var_ent (node)) != MATRIX_DENSE_REAL)
            rerror("idwt: incorrect entry for 'index' of family 'bspline/bspline_c'");
          MDR *ix = ent_data(var_ent(node));
          if (ix->type == RLAB_TYPE_INT32)
          {
            if (MNR(ix)*MNC(ix) == 1)
              idx = MdiV0(ix,0);
            else
              idx = 100 * MdiV0(ix,0) + MdiV0(ix,1);
          }
          else
          {
            if (MNR(ix)*MNC(ix) == 1)
              idx = MdrV0(ix,0);
            else
              idx = 100 * MdrV0(ix,0) + MdrV0(ix,1);
          }
          if (iwav == 5)
            w = gsl_wavelet_alloc(gsl_wavelet_bspline, idx);
          else
            w = gsl_wavelet_alloc(gsl_wavelet_bspline_centered, idx);
          break;
      }
      if (!w)
        rerror("idwt: incorrect entry for 'index' of wavelet family");

    }
    else
      rerror("idwt: missing entry for 'index' of wavelet family");
  }

  work = gsl_wavelet_workspace_alloc (n);
  if (!work)
    rerror("idwt: terrible internal error. cannot allocate memory");

  gsl_wavelet_transform_inverse(w, MDRPTR(xd), 1, n, work);

  gsl_wavelet_free (w);
  gsl_wavelet_workspace_free (work);

  ent_Clean (e1);
  ent_Clean (e2);

  return ent_Assign_Rlab_MDR(xd);
}

// ******************************************************************
//
//  2 - D wavelet transform
//
// ******************************************************************

Ent *
ent_gsl_dwt2 (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0;
  ListNode *node;
  MDR *x1, *xd;
  int n, nr, nc, i, iwav=0, idx, ist = 1;

  gsl_wavelet *w = 0;
  gsl_wavelet_workspace * work = 0;

  if (nargs != 2)
    rerror("dwt2: two arguments required");

  //
  // x
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) != MATRIX_DENSE_REAL)
    rerror("dwt2: first argument has to be square matrix");
  x1 = ent_data (e1);

  nr = x1->nrow;
  nc = x1->ncol;
  if (nr!=nc || nr*nc==0)
    rerror("dwt2: first argument has to be square matrix");

  n = (int) pow(2, (int) (log(nr)/log(2)));
  if (n != nr)
    rerror("dwt2: size of square matrix has to be 2^J, J>0");

  xd = mdr_Float_BF ( x1 );

  //
  // parameters
  //
  e2 = bltin_get_ent (args[1]);
  if (ent_type (e2) == BTREE)
  {
    // family name
    node = btree_FindNode (ent_data (e2), RLAB_NAME_GSL_WV_FAMILY);
    if (node != 0)
    {
      char * fam = class_char_pointer (var_ent (node));
      if (!fam)
        rerror("dwt2: incorrect entry for 'family'");
      iwav = 0;
      for (i = 0; bltin_wavelets[i].index; i++)
      {
        if (!strcmp (fam, bltin_wavelets[i].name))
        {
          iwav = bltin_wavelets[i].index;
          break;
        }
      }
      if (!iwav)
        rerror("dwt2: incorrect entry for 'family'");
    }
    else
      rerror("dwt2: incorrect entry for 'family'");

    // index
    node = btree_FindNode (ent_data (e2), RLAB_NAME_GSL_WV_INDEX);
    if (node != 0)
    {
      idx = 0;
      switch (iwav)
      {
        case 1:
        case 2:
          // daubechies family accepts index
          idx = (int) class_double (var_ent (node));
          if (iwav == 1)
            w = gsl_wavelet_alloc(gsl_wavelet_daubechies, idx);
          else
            w = gsl_wavelet_alloc(gsl_wavelet_daubechies_centered, idx);
          break;

        case 3:
        case 4:
          // haar family accepts index 2, ignore user's entry
          if (iwav == 3)
            w = gsl_wavelet_alloc(gsl_wavelet_haar, 2);
          else
            w = gsl_wavelet_alloc(gsl_wavelet_haar_centered, 2);
          break;

        case 5:
        case 6:
          // bspline family accepts index [i,j] or 100*i+j
          if (ent_type(var_ent (node)) != MATRIX_DENSE_REAL)
            rerror("dwt2: incorrect entry for 'index' of family 'bspline/bspline_c'");
          MDR *ix = ent_data(var_ent(node));
          if (ix->type == RLAB_TYPE_INT32)
          {
            if (MNR(ix)*MNC(ix) == 1)
              idx = MdiV0(ix,0);
            else
              idx = 100 * MdiV0(ix,0) + MdiV0(ix,1);
          }
          else
          {
            if (MNR(ix)*MNC(ix) == 1)
              idx = MdrV0(ix,0);
            else
              idx = 100 * MdrV0(ix,0) + MdrV0(ix,1);
          }
          if (iwav == 5)
            w = gsl_wavelet_alloc(gsl_wavelet_bspline, idx);
          else
            w = gsl_wavelet_alloc(gsl_wavelet_bspline_centered, idx);
          break;
      }
      if (!w)
        rerror("dwt2: incorrect entry for 'index' of wavelet family");

      // procedure: standard (default) or non-standard (ist = 0)
      node = btree_FindNode (ent_data (e2), RLAB_NAME_GSL_WV_PROC);
      if (node != 0)
      {
        ist = (int) class_double (var_ent (node));
        if (ist != 1)
          ist = 0;
      }
    }
    else
      rerror("dwt2: missing entry for 'index' of wavelet family");
  }

  work = gsl_wavelet_workspace_alloc (n);
  if (!work)
    rerror("dwt2: terrible internal error. cannot allocate memory");

  if (ist)
    gsl_wavelet2d_transform_forward(w, MDRPTR(xd), n, n, n, work);
  else
    gsl_wavelet2d_nstransform_forward(w, MDRPTR(xd), n, n, n, work);

  gsl_wavelet_free (w);
  gsl_wavelet_workspace_free (work);

  ent_Clean (e1);
  ent_Clean (e2);

  return ent_Assign_Rlab_MDR(xd);
}

Ent *
ent_gsl_idwt2 (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0;
  ListNode *node;
  MDR *x1, *xd;
  int n, nr, nc, i, iwav=0, idx, ist = 1;

  gsl_wavelet *w = 0;
  gsl_wavelet_workspace * work = 0;

  if (nargs != 2)
    rerror("idwt2: two arguments required");

  //
  // x
  //
  e1 = bltin_get_ent (args[0]);
  if (ent_type (e1) != MATRIX_DENSE_REAL)
    rerror("idwt2: first argument has to be real vector");
  x1 = ent_data (e1);

  nr = x1->nrow;
  nc = x1->ncol;
  if (nr!=nc || nr*nc==0)
    rerror("idwt2: first argument has to be square matrix");

  n = (int) pow(2, (int) (log(nr)/log(2)));
  if (n != nr)
    rerror("idwt2: size of square matrix has to be 2^J, J>0");

  xd = mdr_Float_BF ( x1 );

  //
  // parameters
  //
  e2 = bltin_get_ent (args[1]);
  if (ent_type (e2) == BTREE)
  {
    // family name
    node = btree_FindNode (ent_data (e2), RLAB_NAME_GSL_WV_FAMILY);
    if (node != 0)
    {
      char * fam = class_char_pointer (var_ent (node));
      if (!fam)
        rerror("idwt2: incorrect entry for 'family'");
      iwav = 0;
      for (i = 0; bltin_wavelets[i].index; i++)
      {
        if (!strcmp (fam, bltin_wavelets[i].name))
        {
          iwav = bltin_wavelets[i].index;
          break;
        }
      }
      if (!iwav)
        rerror("idwt2: incorrect entry for 'family'");
    }
    else
      rerror("idwt2: incorrect entry for 'family'");

    // index
    node = btree_FindNode (ent_data (e2), RLAB_NAME_GSL_WV_INDEX);
    if (node != 0)
    {
      idx = 0;
      switch (iwav)
      {
        case 1:
        case 2:
          // daubechies family accepts index
          idx = (int) class_double (var_ent (node));
          if (iwav == 1)
            w = gsl_wavelet_alloc(gsl_wavelet_daubechies, idx);
          else
            w = gsl_wavelet_alloc(gsl_wavelet_daubechies_centered, idx);
          break;

        case 3:
        case 4:
          // haar family accepts index 2, ignore user's entry
          if (iwav == 3)
            w = gsl_wavelet_alloc(gsl_wavelet_haar, 2);
          else
            w = gsl_wavelet_alloc(gsl_wavelet_haar_centered, 2);
          break;

        case 5:
        case 6:
          // bspline family accepts index [i,j] or 100*i+j
          if (ent_type(var_ent (node)) != MATRIX_DENSE_REAL)
            rerror("idwt2: incorrect entry for 'index' of family 'bspline/bspline_c'");
          MDR *ix = ent_data(var_ent(node));
          if (ix->type == RLAB_TYPE_INT32)
          {
            if (MNR(ix)*MNC(ix) == 1)
              idx = MdiV0(ix,0);
            else
              idx = 100 * MdiV0(ix,0) + MdiV0(ix,1);
          }
          else
          {
            if (MNR(ix)*MNC(ix) == 1)
              idx = MdrV0(ix,0);
            else
              idx = 100 * MdrV0(ix,0) + MdrV0(ix,1);
          }
          if (iwav == 5)
            w = gsl_wavelet_alloc(gsl_wavelet_bspline, idx);
          else
            w = gsl_wavelet_alloc(gsl_wavelet_bspline_centered, idx);
          break;
      }
      if (!w)
        rerror("idwt2: incorrect entry for 'index' of wavelet family");

      // procedure: standard (default) or non-standard (ist = 0)
      node = btree_FindNode (ent_data (e2), "proc_standard");
      if (node != 0)
      {
        ist = (int) class_double (var_ent (node));
        if (ist != 1)
          ist = 0;
      }
    }
    else
      rerror("idwt2: missing entry for 'index' of wavelet family");
  }

  work = gsl_wavelet_workspace_alloc (n);
  if (!work)
    rerror("idwt2: terrible internal error. cannot allocate memory");

  if (ist)
    gsl_wavelet2d_transform_inverse(w, MDRPTR(xd), n, n, n, work);
  else
    gsl_wavelet2d_nstransform_inverse(w, MDRPTR(xd), n, n, n, work);

  gsl_wavelet_free (w);
  gsl_wavelet_workspace_free (work);

  ent_Clean (e1);
  ent_Clean (e2);

  return ent_Assign_Rlab_MDR(xd);
}



