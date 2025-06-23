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
#include "mde.h"
#include "sort.h"
#include "rlab_macros.h"

#include <stdio.h>
#include <fcntl.h>

#include "sep.h"
#define THIS_LIBRARY "sep"
#undef  THIS_FILE 
#define THIS_FILE "r_sep.c"

//
#include "rlab_solver_parameters_names.h"
//
// set default wand
//
#undef  THIS_SOLVER
#define THIS_SOLVER "sep.xtr"
Ent * ent_sep_extract (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0;
  MDR *x=0, *n=0, *m=0, *c=0;
  ListNode *node=0;

  double noiseval=0, gain=1.0, maskthresh=0.0;
  int status, i,j;

  if (nargs != 2)
  {
    fprintf (stdout,
             THIS_SOLVER ": sep library wrapper for RLaB\n");
    fprintf (stdout,
             THIS_SOLVER ": Format:\n");
    fprintf (stdout,
             THIS_SOLVER ":   " THIS_SOLVER "(data, options),\n");
    fprintf (stdout,
             THIS_SOLVER ": where\n");
    fprintf (stdout,
             THIS_SOLVER ":   'data' is data matrix, and options = <<>>.\n");
    rerror ( THIS_SOLVER ": requires 2 arguments");
  }

  // data matrix from multiple sources:
  //  e.g., pixel matrix from image magick
  //  or just raw data from who knows whom...
  e1 = bltin_get_ent (args[0]);
  if (ent_type(e1) != MATRIX_DENSE_REAL)
    rerror(THIS_SOLVER ": " RLAB_ERROR_ARG1_MDR_MATRIX "\n");
  x = ent_data(e1);
  if (SIZE(x)<1)
    rerror(THIS_SOLVER ": " RLAB_ERROR_ARG1_MDR_MATRIX "\n");

  e2 = bltin_get_ent (args[1]);
  if (ent_type(e2) != BTREE)
    rerror(THIS_SOLVER ": " RLAB_ERROR_ARG2_BTREE "\n");

  // populate default SEP image structure
  sep_image im;
  im.data = MDRPTR(x);
  im.dtype = MD_TYPE_INT32(x) ? SEP_TINT : SEP_TDOUBLE;
  im.w = MNR(x);
  im.h = MNC(x);

  //
  // art of noise
  // 

  // no noise
  im.noise = NULL;
  im.ndtype = 0;
  im.noise = NULL;
  im.noiseval = 0;
  im.noise_type = SEP_NOISE_NONE;

  node = btree_FindNode (ent_data (e2), RLAB_NAME_SEP_BKG_NOISESTD);
  if (node != 0)
  {
    if (ent_type(var_ent (node))!=MATRIX_DENSE_REAL)
      rerror(THIS_SOLVER ": options entry '"RLAB_NAME_SEP_BKG_NOISESTD"' must be positive scalar\n");

    im.noise = NULL;
    noiseval = class_double (var_ent (node));
    if (noiseval > 0)
    {
      im.noiseval = noiseval;
      im.noise_type = SEP_NOISE_STDDEV;
    }
  }
  node = btree_FindNode (ent_data (e2), RLAB_NAME_SEP_BKG_NOISEVAR);
  if (node != 0)
  {
    if (ent_type(var_ent (node))!=MATRIX_DENSE_REAL)
      rerror(THIS_SOLVER ": options entry '"RLAB_NAME_SEP_BKG_NOISEVAR"' must be positive scalar\n");

    im.noise = NULL;
    noiseval = class_double (var_ent (node));
    if (noiseval > 0)
    {
      im.noiseval = noiseval;
      im.noise_type = SEP_NOISE_VAR;
    }
  }
  node = btree_FindNode (ent_data (e2), RLAB_NAME_SEP_BKG_NOISEMAT);
  if (node != 0)
  {
    if (ent_type(var_ent (node))!=MATRIX_DENSE_REAL)
      rerror(THIS_SOLVER ": options entry '"RLAB_NAME_SEP_BKG_NOISEMAT"' must be matrix of the same size as data\n");
    n = ent_data(var_ent (node));
    if ((MNR(n)!=MNR(x))||(MNC(n)!=MNC(x)))
      rerror(THIS_SOLVER ": options entry '"RLAB_NAME_SEP_BKG_NOISEMAT"' must be matrix of the same size as data\n");
    im.noise  = MDRPTR(n);
    im.ndtype = MD_TYPE_INT32(n) ? SEP_TINT : SEP_TDOUBLE;
  }

  // gain
  node = btree_FindNode (ent_data (e2), RLAB_NAME_SEP_BKG_GAIN);
  if (node != 0)
  {
    if (ent_type(var_ent (node))!=MATRIX_DENSE_REAL)
      rerror(THIS_SOLVER ": options entry '"RLAB_NAME_SEP_BKG_GAIN"' must be positive scalar\n");
    gain = class_double (var_ent (node));
    if (gain <= 0.0)
      gain = 1.0;
  }
  im.gain = gain;

  // mask
  im.mask  = NULL;
  im.mdtype = 0;
  node = btree_FindNode (ent_data (e2), RLAB_NAME_SEP_BKG_MASK);
  if (node != 0)
  {
    m = ent_data(var_ent (node));
    if (ent_type(var_ent (node))!=MATRIX_DENSE_REAL)
      rerror(THIS_SOLVER ": options entry '"RLAB_NAME_SEP_BKG_MASK"' must be matrix\n");
    im.mask = MDRPTR(m);
    im.mdtype = MD_TYPE_INT32(m) ? SEP_TINT : SEP_TDOUBLE;
  }

  // mask threshold
  node = btree_FindNode (ent_data (e2), RLAB_NAME_SEP_BKG_MASK_THR);
  if (node != 0)
  {
    if (ent_type(var_ent (node))!=MATRIX_DENSE_REAL)
      rerror(THIS_SOLVER ": options entry '"RLAB_NAME_SEP_BKG_MASK_THR"' must be positive scalar\n");
    maskthresh = class_double (var_ent (node));
    if (maskthresh <= 0.0)
      maskthresh = 0.0;
  }
  im.maskthresh = maskthresh;

  //
  // prepare for source extraction
  //
  double dummy;
  float thresh = 1.5;
  int thresh_type = SEP_THRESH_REL;
  // relative threshold
  node = btree_FindNode (ent_data (e2), RLAB_NAME_SEP_XTR_REL_THR);
  if (node != 0)
  {
    if (ent_type(var_ent (node))!=MATRIX_DENSE_REAL)
      rerror(THIS_SOLVER ": options entry '"RLAB_NAME_SEP_XTR_REL_THR"' must be positive scalar\n");
    dummy = class_double (var_ent (node));
    if (dummy > 0 )
    {
      thresh = dummy;
      thresh_type = SEP_THRESH_REL;
    }
  }
  // absolute threshold
  node = btree_FindNode (ent_data (e2), RLAB_NAME_SEP_XTR_ABS_THR);
  if (node != 0)
  {
    if (ent_type(var_ent (node))!=MATRIX_DENSE_REAL)
      rerror(THIS_SOLVER ": options entry '"RLAB_NAME_SEP_XTR_ABS_THR"' must be positive scalar\n");
    dummy = class_double (var_ent (node));
    if (dummy > 0 )
    {
      thresh = dummy;
      thresh_type = SEP_THRESH_ABS;
    }
  }

  // minimum area 
  int idummy;
  int minarea = 5;
  node = btree_FindNode (ent_data (e2), RLAB_NAME_SEP_XTR_MINAREA);
  if (node != 0)
  {
    if (ent_type(var_ent (node))!=MATRIX_DENSE_REAL)
      rerror(THIS_SOLVER ": options entry '"RLAB_NAME_SEP_XTR_MINAREA"' must be positive scalar\n");
    idummy = class_double (var_ent (node));
    if (idummy > 0 )
    {
      minarea = idummy;
    }
  }

  // maximum deblend area 
  node = btree_FindNode (ent_data (e2), RLAB_NAME_SEP_XTR_MAXDEBAREA);
  if (node != 0)
  {
    if (ent_type(var_ent (node))!=MATRIX_DENSE_REAL)
      rerror(THIS_SOLVER ": options entry '"RLAB_NAME_SEP_XTR_MAXDEBAREA"' must be positive scalar\n");
    idummy = class_double (var_ent (node));
    if (idummy > 0 )
    {
      sep_lutz_max_deb_area(&idummy);
    }
  }

  // deblend threshold type 
  node = btree_FindNode (ent_data (e2), RLAB_NAME_SEP_XTR_DEBTHRESHTYPE);
  if (node != 0)
  {
    if (ent_type(var_ent (node))!=MATRIX_DENSE_REAL)
      rerror(THIS_SOLVER ": options entry '"RLAB_NAME_SEP_XTR_DEBTHRESHTYPE"' must be 0 (power) or 1 (linear)\n");
    idummy = (class_double (var_ent (node)) != 0.0);
    sep_lutz_threshold_type(&idummy);
  }
  // convolution array
  float *conv=0;
  int convw=0, convh=0;
  int filter_type = SEP_FILTER_MATCHED;
  node = btree_FindNode (ent_data (e2), RLAB_NAME_SEP_XTR_CONV);
  if (node != 0)
  {
    if (ent_type(var_ent (node))!=MATRIX_DENSE_REAL)
      rerror(THIS_SOLVER ": options entry '"RLAB_NAME_SEP_XTR_CONV"' must be real matrix\n");
    c = ent_data(var_ent (node));
    if (SIZE(c)<1)
      rerror(THIS_SOLVER ": options entry '"RLAB_NAME_SEP_XTR_CONV"' must be real matrix\n");
    convw = MNR(c);
    convh = MNC(c);
    conv = GC_MALLOC(convw * convh * sizeof(float));
    for (i=0; i<SIZE(c); i++)
      conv[i] = MdrV0(c,i);
    filter_type = SEP_FILTER_CONV;
  }

  // deblend
  int deblend_nthresh = 32;
  node = btree_FindNode (ent_data (e2), RLAB_NAME_SEP_XTR_DEBLEND);
  if (node != 0)
  {
    if (ent_type(var_ent (node))!=MATRIX_DENSE_REAL)
      rerror(THIS_SOLVER ": options entry '"RLAB_NAME_SEP_XTR_DEBLEND"' must be positive integer\n");
    idummy = class_double (var_ent (node));
    if (idummy > 0 )
    {
      deblend_nthresh = idummy;
    }
  }

  // contrast
  double deblend_cont = 0.005;
  node = btree_FindNode (ent_data (e2), RLAB_NAME_SEP_XTR_CONTRAST);
  if (node != 0)
  {
    if (ent_type(var_ent (node))!=MATRIX_DENSE_REAL)
      rerror(THIS_SOLVER ": options entry '"RLAB_NAME_SEP_XTR_CONTRAST"' must be positive scalar\n");
    dummy = class_double (var_ent (node));
    if (dummy > 0 && dummy<=1)
    {
      deblend_cont = dummy;
    }
  }

  // clean
  int clean_flag = 1;
  double clean_param = 1.0;
  node = btree_FindNode (ent_data (e2), RLAB_NAME_SEP_XTR_CLEAN);
  if (node != 0)
  {
    if (ent_type(var_ent (node))!=MATRIX_DENSE_REAL)
      rerror(THIS_SOLVER ": options entry '"RLAB_NAME_SEP_XTR_CLEAN"' must be positive scalar\n");
    dummy = class_double (var_ent (node));
    if (dummy <= 0 )
    {
      clean_flag = 0;
    }
    else
    {
      clean_flag = 1;
      clean_param = dummy;
    }
  }

  // set pixel stack
  int px_stack = 0;
  node = btree_FindNode (ent_data (e2), RLAB_NAME_SEP_XTR_PIXEL_STACK);
  if (node != 0)
  {
    if (ent_type(var_ent (node))!=MATRIX_DENSE_REAL)
      rerror(THIS_SOLVER ": options entry '"RLAB_NAME_SEP_XTR_PIXEL_STACK"' must be positive scalar\n");
    idummy = class_double (var_ent (node));
    if (idummy > 0 )
    {
      px_stack = idummy;
    }
  }

  // sort by
  int sort_by = 1; // 0-unsorted, 1-by flux, 2-by peak
  node = btree_FindNode (ent_data (e2), RLAB_NAME_SEP_XTR_SORT_BY);
  if (node != 0)
  {
    if (ent_type(var_ent (node))!=MATRIX_DENSE_STRING)
      rerror(THIS_SOLVER ": options entry '"RLAB_NAME_SEP_XTR_SORT_BY"' must be string scalar\n");
    char *s = class_char_pointer (var_ent (node));
    if (!strcmp(s,"flux"))
      sort_by=1;
    else if (!strcmp(s,"peak"))
      sort_by=2;
    else
      sort_by=0;
  }

  // do it
  sep_catalog *catalog = NULL;
  status = sep_extract(&im, thresh, thresh_type, minarea,
                        conv, convw, convh, filter_type,
                        deblend_nthresh, deblend_cont,
                        clean_flag, clean_param, &catalog);
  MDR *xy_val=0, *m2_val=0, *el_val=0, *peak_val=0, *cpeak_val=0;
  MDR *peak_xy_val=0, *cpeak_xy_val=0, *flux_val=0, *cflux_val=0;
  MDE *pix_arrays=0;
  Btree *rval=0;
  if (!status)
  {
    flux_val = mdr_Create(catalog->nobj, 1);
    peak_val = mdr_Create(catalog->nobj, 1);

    MDR *r_idx = mdr_Create(1, catalog->nobj);
    for (i=0; i<catalog->nobj; i++)
      MdrV0(r_idx,i) = i;
    if (sort_by == 1)
    {
      // sort in ascending order
      for (i=0; i<catalog->nobj; i++)
        MdrV0(flux_val,i) = catalog->flux[i];
      r_sort ((double *) MDRPTR(flux_val), 0, catalog->nobj-1, (double *)MDRPTR(r_idx));
      // now flipit upside down to be descending order
      for (i=0; i<(catalog->nobj/2); i++)
      {
        SWAP ( MdrV0(flux_val,i),  MdrV0(flux_val,catalog->nobj-1-i), double);
        SWAP ( MdrV0(r_idx,i),  MdrV0(r_idx,catalog->nobj-1-i), double);
      }
      for (i=0; i<catalog->nobj; i++)
        MdrV0(peak_val,i) = catalog->peak[ (int) MdrV0(r_idx,i) ];
    }
    else if (sort_by == 2)
    {
      for (i=0; i<catalog->nobj; i++)
        MdrV0(peak_val,i) = catalog->peak[i];
      r_sort ((double *) MDRPTR(peak_val), 0, catalog->nobj-1, (double *)MDRPTR(r_idx));
      // now flipit upside down to be descending order
      for (i=0; i<(catalog->nobj/2); i++)
      {
        SWAP ( MdrV0(peak_val,i),  MdrV0(peak_val,catalog->nobj-1-i), double);
        SWAP ( MdrV0(r_idx,i),  MdrV0(r_idx,catalog->nobj-1-i), double);
      }
      for (i=0; i<catalog->nobj; i++)
        MdrV0(flux_val,i) = catalog->flux[ (int) MdrV0(r_idx,i) ];
    }
    xy_val = mdr_Create(catalog->nobj, 2);
    m2_val = mdr_Create(catalog->nobj, 3);
    el_val = mdr_Create(catalog->nobj, 3);
    peak_xy_val = mdr_Create(catalog->nobj, 2);
    pix_arrays = mde_Create(catalog->nobj, 1);
    if (conv)
    {
      cpeak_val = mdr_Create(catalog->nobj, 1);
      cpeak_xy_val = mdr_Create(catalog->nobj, 2);
      cflux_val = mdr_Create(catalog->nobj, 1);
    }
    for (i=0; i<catalog->nobj; i++)
    {
      // barycenters
      Mdr0(xy_val,i,0) = catalog->x[ (int) MdrV0(r_idx,i) ] + 1.0;
      Mdr0(xy_val,i,1) = catalog->y[ (int) MdrV0(r_idx,i) ] + 1.0;
      // second momenta
      Mdr0(m2_val,i,0) = catalog->x2[ (int) MdrV0(r_idx,i) ];
      Mdr0(m2_val,i,1) = catalog->y2[ (int) MdrV0(r_idx,i) ];
      Mdr0(m2_val,i,2) = catalog->xy[ (int) MdrV0(r_idx,i) ];
      // ellipse
      Mdr0(el_val,i,0) = catalog->cxx[ (int) MdrV0(r_idx,i) ];
      Mdr0(el_val,i,1) = catalog->cyy[ (int) MdrV0(r_idx,i) ];
      Mdr0(el_val,i,2) = catalog->cxy[ (int) MdrV0(r_idx,i) ];
      // peak positions
      Mdr0(peak_xy_val,i,0) = catalog->xpeak[ (int) MdrV0(r_idx,i) ] + 1.0;
      Mdr0(peak_xy_val,i,1) = catalog->ypeak[ (int) MdrV0(r_idx,i) ] + 1.0;
      if (conv)
      {
        MdrV0(cpeak_val,i) = catalog->cpeak[ (int) MdrV0(r_idx,i) ];
        Mdr0(cpeak_xy_val,i,0) = catalog->xcpeak[ (int) MdrV0(r_idx,i) ] + 1.0;
        Mdr0(cpeak_xy_val,i,1) = catalog->ycpeak[ (int) MdrV0(r_idx,i) ] + 1.0;
        MdrV0(cflux_val,i) = catalog->cflux[ (int) MdrV0(r_idx,i) ];
      }
      // pixel arrays
      MDR *pix = mdi_Create(catalog->npix[ (int) MdrV0(r_idx,i) ],2);
      for (j=0; j<catalog->npix[ (int) MdrV0(r_idx,i) ]; j++)
      {
        Mdi0(pix,j,0) = (catalog->pix[ (int) MdrV0(r_idx,i) ][j] / im.w) + 1;
        Mdi0(pix,j,1) = (catalog->pix[ (int) MdrV0(r_idx,i) ][j] % im.w) + 1;
      }
      MdeV0(pix_arrays,i) = ent_Assign_Rlab_MDR(pix);
    } // for (i=0; i<catalog->nobj; i++)

    mdr_Destroy (r_idx);

    // assign!
    rval = btree_Create();
    install(rval, RLAB_NAME_SEP_XTR_CENTER,   ent_Assign_Rlab_MDR(xy_val));
    install(rval, RLAB_NAME_SEP_XTR_SECMOMNT, ent_Assign_Rlab_MDR(m2_val));
    install(rval, RLAB_NAME_SEP_XTR_ELLIPSE,  ent_Assign_Rlab_MDR(el_val));
    install(rval, RLAB_NAME_SEP_XTR_PEAK,     ent_Assign_Rlab_MDR(peak_val));
    install(rval, RLAB_NAME_SEP_XTR_PEAK_POS, ent_Assign_Rlab_MDR(peak_xy_val));
    install(rval, RLAB_NAME_SEP_XTR_FLUX,     ent_Assign_Rlab_MDR(flux_val));
    if (conv)
    {
      install(rval, RLAB_NAME_SEP_XTR_CPEAK,     ent_Assign_Rlab_MDR(cpeak_val));
      install(rval, RLAB_NAME_SEP_XTR_CPEAK_POS, ent_Assign_Rlab_MDR(cpeak_xy_val));
      install(rval, RLAB_NAME_SEP_XTR_CFLUX,     ent_Assign_Rlab_MDR(cflux_val));
    }
    install(rval, RLAB_NAME_SEP_XTR_PIXELS,     ent_Assign_Rlab_MDE(pix_arrays));
  }
  else
  {
    char errtext[61]={'\0'};
    sep_get_errmsg(status, errtext);
    printf(THIS_FILE ": " THIS_SOLVER ": Operation failed with status=%i: %s\n",
           status, errtext);
  }

  // cleanup sep
  sep_catalog_free(catalog);

  // clean up rlab
  if (conv)
    GC_FREE(conv);
  ent_Clean(e1);
  ent_Clean(e2);

  return ent_Assign_Rlab_BTREE(rval);
}

#undef  THIS_SOLVER
#define THIS_SOLVER "sep.bkg"
Ent * ent_sep_background (int nargs, Datum args[])
{
  Ent *e1=0, *e2=0;
  MDR *rval=0, *x=0, *n=0, *m=0, *t=0, *f=0;
  ListNode *node=0;

  double noiseval=0, gain=1.0, maskthresh=0.0, fthresh=0.0;
  int bw=64, bh=64, fw=3, fh=3, status;

  if (nargs != 2)
  {
    fprintf (stdout,
             THIS_SOLVER ": sep library wrapper for RLaB\n");
    fprintf (stdout,
             THIS_SOLVER ": Format:\n");
    fprintf (stdout,
             THIS_SOLVER ":   " THIS_SOLVER "(data, options),\n");
    fprintf (stdout,
             THIS_SOLVER ": where\n");
    fprintf (stdout,
             THIS_SOLVER ":   'data' is data matrix, and options = <<>>.\n");
    rerror ( THIS_SOLVER ": requires 2 arguments");
  }

  // data matrix from multiple sources:
  //  e.g., pixel matrix from image magick
  //  or just raw data from who knows whom...
  e1 = bltin_get_ent (args[0]);
  if (ent_type(e1) != MATRIX_DENSE_REAL)
    rerror(THIS_SOLVER ": " RLAB_ERROR_ARG1_MDR_MATRIX "\n");
  x = ent_data(e1);
  if (SIZE(x)<1)
    rerror(THIS_SOLVER ": " RLAB_ERROR_ARG1_MDR_MATRIX "\n");

  e2 = bltin_get_ent (args[1]);
  if (ent_type(e2) != BTREE)
    rerror(THIS_SOLVER ": " RLAB_ERROR_ARG2_BTREE "\n");

  // populate default SEP image structure
  sep_image im;
  im.data = MDRPTR(x);
  im.dtype = MD_TYPE_INT32(x) ? SEP_TINT : SEP_TDOUBLE;
  im.w = MNR(x);
  im.h = MNC(x);

  //
  // art of noise
  //

  // no noise
  im.noise = NULL;
  im.ndtype = 0;
  im.noise = NULL;
  im.noiseval = 0;
  im.noise_type = SEP_NOISE_NONE;

  node = btree_FindNode (ent_data (e2), RLAB_NAME_SEP_BKG_NOISESTD);
  if (node != 0)
  {
    if (ent_type(var_ent (node))!=MATRIX_DENSE_REAL)
      rerror(THIS_SOLVER ": options entry 'noise_std' must be positive scalar\n");

    im.noise = NULL;
    noiseval = class_double (var_ent (node));
    if (noiseval > 0)
    {
      im.noiseval = noiseval;
      im.noise_type = SEP_NOISE_STDDEV;
    }
  }
  node = btree_FindNode (ent_data (e2), RLAB_NAME_SEP_BKG_NOISEVAR);
  if (node != 0)
  {
    if (ent_type(var_ent (node))!=MATRIX_DENSE_REAL)
      rerror(THIS_SOLVER ": options entry 'noise_var' must be positive scalar\n");

    im.noise = NULL;
    noiseval = class_double (var_ent (node));
    if (noiseval > 0)
    {
      im.noiseval = noiseval;
      im.noise_type = SEP_NOISE_VAR;
    }
  }
  node = btree_FindNode (ent_data (e2), RLAB_NAME_SEP_BKG_NOISEMAT);
  if (node != 0)
  {
    if (ent_type(var_ent (node))!=MATRIX_DENSE_REAL)
      rerror(THIS_SOLVER ": options entry 'noise' must be matrix of the same size as data\n");
    n = ent_data(var_ent (node));
    if ((MNR(n)!=MNR(x))||(MNC(n)!=MNC(x)))
      rerror(THIS_SOLVER ": options entry 'noise' must be matrix of the same size as data\n");
    im.noise  = MDRPTR(n);
    im.ndtype = MD_TYPE_INT32(n) ? SEP_TINT : SEP_TDOUBLE;
  }

  // gain
  node = btree_FindNode (ent_data (e2), RLAB_NAME_SEP_BKG_GAIN);
  if (node != 0)
  {
    if (ent_type(var_ent (node))!=MATRIX_DENSE_REAL)
      rerror(THIS_SOLVER ": options entry 'gain' must be positive scalar\n");
    gain = class_double (var_ent (node));
    if (gain <= 0.0)
      gain = 1.0;
  }
  im.gain = gain;

  // mask
  im.mask  = NULL;
  im.mdtype = 0;
  node = btree_FindNode (ent_data (e2), RLAB_NAME_SEP_BKG_MASK);
  if (node != 0)
  {
    m = ent_data(var_ent (node));
    if (ent_type(var_ent (node))!=MATRIX_DENSE_REAL)
      rerror(THIS_SOLVER ": options entry 'mask' must be matrix\n");
    im.mask = MDRPTR(m);
    im.mdtype = MD_TYPE_INT32(m) ? SEP_TINT : SEP_TDOUBLE;
  }

  // mask threshold
  node = btree_FindNode (ent_data (e2), RLAB_NAME_SEP_BKG_MASK_THR);
  if (node != 0)
  {
    if (ent_type(var_ent (node))!=MATRIX_DENSE_REAL)
      rerror(THIS_SOLVER ": options entry 'mask_thresh' must be positive scalar\n");
    maskthresh = class_double (var_ent (node));
    if (maskthresh <= 0.0)
      maskthresh = 0.0;
  }
  im.maskthresh = maskthresh;

  //
  // prepare for background estimation
  //
  // background tile: 
  node = btree_FindNode (ent_data (e2), RLAB_NAME_SEP_BKG_1TILE);
  if (node != 0)
  {
    if (ent_type(var_ent (node))!=MATRIX_DENSE_REAL)
      rerror(THIS_SOLVER ": options entry '"RLAB_NAME_SEP_BKG_1TILE"' must be positive integer pair\n");
    t = ent_data(var_ent (node));
    if (SIZE(t) != 2)
      rerror(THIS_SOLVER ": options entry '"RLAB_NAME_SEP_BKG_1TILE"' must be positive integer pair\n");
    bw = mdiV0(t,0);
    bh = mdiV0(t,1);
  }
  // background filter in tiles: 
  node = btree_FindNode (ent_data (e2), RLAB_NAME_SEP_BKG_FILTER);
  if (node != 0)
  {
    if (ent_type(var_ent (node))!=MATRIX_DENSE_REAL)
      rerror(THIS_SOLVER ": options entry '"RLAB_NAME_SEP_BKG_FILTER"' must be positive integer pair\n");
    f = ent_data(var_ent (node));
    if (SIZE(f) != 2)
      rerror(THIS_SOLVER ": options entry '"RLAB_NAME_SEP_BKG_FILTER"' must be positive integer pair\n");
    fw = mdiV0(f,0);
    fh = mdiV0(f,1);
  }
  // filter threshold
  node = btree_FindNode (ent_data (e2), RLAB_NAME_SEP_BKG_FILTER_THR);
  if (node != 0)
  {
    if (ent_type(var_ent (node))!=MATRIX_DENSE_REAL)
      rerror(THIS_SOLVER ": options entry '"RLAB_NAME_SEP_BKG_FILTER_THR"' must be positive scalar\n");
    fthresh = class_double (var_ent (node));
    if (fthresh <= 0.0)
      fthresh = 0.0;
  }

  // do it
  sep_bkg *bkg = NULL;
  status = sep_background(&im, bw, bh, fw, fh, fthresh, &bkg);

  if (!status)
  {
    // return the same type of background as was the data matrix:
    if (MD_TYPE_INT32(x))
    {
      rval = mdi_Create(MNR(x),MNC(x));
      status = sep_bkg_array(bkg, MDRPTR(rval), SEP_TINT);
    }
    else
    {
      rval = mdr_Create(MNR(x),MNC(x));
      status = sep_bkg_array(bkg, MDRPTR(rval), SEP_TDOUBLE);
    }
  }

  // cleanup sep
  sep_bkg_free(bkg);

  // clean up rlab
  ent_Clean(e1);
  ent_Clean(e2);

  return ent_Assign_Rlab_MDR(rval);
}


