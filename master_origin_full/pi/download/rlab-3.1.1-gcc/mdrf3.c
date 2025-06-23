/*
 * mdrf3.c
 * Matrix-Dense-Real with and without FORTRAN77
 */

/*  This file is a part of rlabplus ("Our"-LaB)
    RLaB Copyright (C) 1995  Ian R. Searle,
    rlabplus (C) 2005-2008 M. Kostrun

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
#include "mdc.h"
#include "util.h"
#include "mathl.h"

#include "fi.h"
#include "lp.h"
#include "blas.h"

#include <stdio.h>
#include <math.h>

/*
 * Create Infs, and NaNs
 */
#include "mathl.h"
#include "rlab_macros.h"

// #define DEBUG

static double abs_err=1e-3, rel_err=1e-3;

/**
  * rotate vector [x,y,z] around vector [k_x,k_y,k_z]
  * and store result in [x,y,z]
  * using rodriguez formula
  *    https://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula
  *
  * \vec v_{rot} = \vec v \, \cos\theta +
  *     (\vec k \cross \vec v) \, \sin\theta +
  *     \vec k \, (\vec k \cdot \vec v)\,(1 - \cos\theta)
  */
int rotate3d(double *x, double *y, double *z, double kx, double ky, double kz, double theta)
{
  double dnk = kx*kx + ky*ky + kz*kz, xx, yy, zz;

  // check that the norm is not too small
  if (dnk <= abs_err)
    return 1;

  double cos_th=cos(theta), sin_th=sin(theta), one_minus_cos_th=1-cos_th;

  xx = (*x) * cos_th;
  yy = (*y) * cos_th;
  zz = (*z) * cos_th;

  double dotp = kx * (*x) + ky * (*y) + kz * (*z);
  dotp *= one_minus_cos_th / dnk;
  xx += kx * dotp;
  yy += ky * dotp;
  zz += kz * dotp;

  dnk = sin_th / sqrt(dnk);
  xx += (ky * (*z) - kz * (*y)) * dnk;
  yy += (kz * (*x) - kx * (*z)) * dnk;
  zz += (kx * (*y) - ky * (*x)) * dnk;

  // update return
  *x = xx;
  *y = yy;
  *z = zz;
  return 0;
}

//
// rotate vector around axis:
//    x -> rot(x,v)
//
int mdr_rotate_vector_axis(MDR *x, double theta, int axis_idx)
{
  if (MNC(x)!=3)
    return 1;

  int i;
  double kx=0.0, ky=0.0, kz=0.0;
  if (axis_idx==1)
  {
    kx=1.0;
  }
  else if (axis_idx==2)
  {
    ky=1.0;
  }
  else if (axis_idx==3)
  {
    kz=1.0;
  }
  else return 2;

  for (i=0; i<MNR(x); i++)
  {
    if (rotate3d(&Mdr0(x,i,0), &Mdr0(x,i,1), &Mdr0(x,i,2), kx, ky, kz, theta))
    {
      return 3;
    }
  }

  return 0;
}

//
//
//
//
//
//
int mdr_bwlabel (MDR *x, int conn, int * pixel_count, MDR **npix, MDR **xy_bbox)
{
  int *blob_label=NULL;
  int i, j, k, l, nr=MNR(x), nc=MNC(x);
  int curr_label, notfound_first=1;
  int x1,x2,x3,x4,m,sx,sx1,sx2,sx3,sx4;

  // default initiation of the return values:
  *npix = NULL;
  *xy_bbox = NULL;
  *pixel_count = 0;

  if (conn!=4 && conn!=8)
    return -1;

  blob_label = GC_MALLOC(nr * (nc>>2) * sizeof(int));

  //
  // process rows concluding with one in which the first nonzero pixel was found
  //
  i=0;
  curr_label = 1;
  while(notfound_first)
  {
    for (j=0; j<nc; j++)
    {
      if (Mdr0(x,i,j)==0)
        continue;

      (*blob_label)++;

      if (notfound_first)
      {
        notfound_first=0;
        blob_label[curr_label] = curr_label;
        Mdr0(x,i,j) = sign(Mdr0(x,i,j)) * curr_label;
        curr_label++;
        continue;
      }

      if ( sign(Mdr0(x,i,j)) == sign(Mdr0(x,i,j-1)) )
      {
        Mdr0(x,i,j) = Mdr0(x,i,j-1);
      }
      else
      {
        blob_label[curr_label] = curr_label;
        Mdr0(x,i,j) = sign(Mdr0(x,i,j)) * curr_label;
        curr_label++;
      }
    }
    i++;
    if (i==nr)
      break;
  }

  //
  // if this is last row then there are no pixels
  //
  if (notfound_first)
  {
    if (blob_label)
      GC_FREE(blob_label);
    return 0;
  }

  if (conn == 4)
  {
    for (;i<nr;i++)
    {
      // as i>=1, in each  row, first check the first column
      // above (conn=4) and above-right (conn=8)
      if (Mdr0(x,i,0) != 0)
      {
        (*blob_label)++; // count hot pixels

        if ( sign(Mdr0(x,i-1,0)) == sign(Mdr0(x,i,0)) )
          Mdr0(x,i,0) = Mdr0(x,i-1,0);
        else
        {
          blob_label[curr_label] = curr_label;
          Mdr0(x,i,0) = sign(Mdr0(x,i,0)) * curr_label;
          curr_label++;
        }
      }

      for(j=1;j<nc;j++)
      {
        if (Mdr0(x,i,j)==0)
          continue;

        sx = sign(Mdr0(x,i,j));
        x2  = Mdr0(x,i-1,j  );
        sx2 = sign(x2);
        x4  = Mdr0(x,i  ,j-1);
        sx4 = sign(x4);

        (*blob_label)++; // count hot pixels

        if ( sx == sx2 )
        {
          Mdr0(x,i,j) = x2;

          if ( sx == sx4 )
          {
            if (x2 != x4 )
            {
#ifdef DEBUG
              printf("min,max = %i, %i\n", (int) MIN(Mdr0(x,i,j-1), Mdr0(x,i,j)) , (int) MAX(Mdr0(x,i,j-1), Mdr0(x,i,j)) );
#endif
              blob_label[ (int) MAX( ABS(x4), ABS(x2) ) ] = (int) MIN( ABS(x4), ABS(x2) );
            }
          }
          continue;
        }

        if ( sx == sx4 )
        {
          Mdr0(x,i,j) = x4;
        }
        else
        {
          blob_label[curr_label] = curr_label;
          Mdr0(x,i,j) = sx * curr_label;
          curr_label++;
        }
      }
    }
  }
  else if (conn == 8)
  {
    for (;i<nr;i++)
    {
      // as i>=1, in each  row, first check the first column
      // above (conn=4) and above-right (conn=8)
      if (Mdr0(x,i,0) != 0)
      {
        (*blob_label)++; // count hot pixels

        sx  = sign(Mdr0(x,i,0));

        x1  = Mdr0(x,i-1,1); // above right
        sx1 = sign(x1);
        x2  = Mdr0(x,i-1,0); // above
        sx2 = sign(x2);


        if ((sx1 == sx) && (sx2 == sx))
        {
          Mdr0(x,i,0) = sx * MIN(ABS(x1),ABS(x2));
          blob_label[ (int) MAX(ABS(x1), ABS(x2)) ] = ABS((int) Mdr0(x,i,0));
        }
        else if ( sx2 == sx )
        {
          Mdr0(x,i,0) = x2;
        }
        else if ( sx1 == sx )
        {
          Mdr0(x,i,0) = x1;
        }
        else
        {
          blob_label[curr_label] = curr_label;
          Mdr0(x,i,0) = sx * curr_label;
          curr_label++;
        }
      }
      for(j=1;j<nc;j++)
      {
        if (Mdr0(x,i,j) == 0)
          continue;

        (*blob_label)++; // count hot pixels

        sx = sign(Mdr0(x,i,j));

        m  = 0;
        x1 = 0;
        if (j<nc-1)
        {
          x1 = Mdr0(x,i-1,j+1);
        }
        sx1 = sign(x1);
        x2  = Mdr0(x,i-1,j  );
        sx2 = sign(x2);
        x3  = Mdr0(x,i-1,j-1);
        sx3 = sign(x3);
        x4  = Mdr0(x,i  ,j-1);
        sx4 = sign(x4);

        if ((sx1 != sx) && (sx2 != sx)  && (sx3 != sx)  && (sx4 != sx) )
        {
          blob_label[curr_label] = curr_label;
          Mdr0(x,i,j) = sx * curr_label;
          curr_label++;
          continue;
        }

        if (j<nc-1)
          if (sx1 == sx)
            m = x1;

        if (sx2 == sx)
        {
          if (m==0)
            m = x2;
          else if (m!=x2)
          {
            blob_label[ (int) MAX(ABS(m), ABS(x2)) ] = (int) MIN(ABS(m),ABS(x2));
            if (ABS(m)>ABS(x2))
              m = x2;
          }
        }

        if (sx3 == sx)
        {
          if (m==0)
            m = x3;
          else if (m!=x3)
          {
            blob_label[ (int) MAX(ABS(m), ABS(x3)) ] = (int) MIN(ABS(m),ABS(x3));
            if (ABS(m) > ABS(x3))
              m = x3;
          }
        }

        if (sx4 == sx)
        {
          if (m==0)
            m = x4;
          else if (m!=x4)
          {
            blob_label[ (int) MAX(ABS(m), ABS(x4)) ] = (int) MIN(ABS(m),ABS(x4));
            if (ABS(m)>ABS(x4))
              m = x4;
          }
        }

        Mdr0(x,i,j) = m;
        continue;

      }
    }
  }
  curr_label--;

  // tell caller how many pixels are there
  *pixel_count = *blob_label;

  int new_count_labels = curr_label;

#ifdef DEBUG
  printf("\nhot pixels = %i\n", *blob_label);
  printf("before: %i\n", new_count_labels);
  for (k=1; k<=curr_label; k++)
  {
    printf("%i -> %i\n", k, blob_label[k]);
  }
#endif

  //
  // first pass: reduce all labels, and find their final number
  //
  new_count_labels = 0;
  for (k=1; k<=curr_label; k++)
  {
    while (blob_label[k] !=  blob_label[blob_label[k]])
    {
      blob_label[k] =  blob_label[blob_label[k]];
    }
    if (k == blob_label[k])
      new_count_labels++;
  }

  //
  // do we rewrite labels?
  //
  if (new_count_labels < curr_label)
  {
    //
    // second pass: find labels below new_count_labels that are not used
    //
    for (k=1; k<=new_count_labels; k++)
    {
      int goto_next_k;
      if (k != blob_label[k])
      {
        // label 'k' is not used
#ifdef DEBUG
        printf("label '%i' is not used\n", k);
#endif

        // find label 'j' above new_count_labels
        goto_next_k = 0;
        for (j=1; j<=curr_label; j++)
        {
          if ((j == blob_label[j]) && j>new_count_labels)
          {
            for (l=1; l<=curr_label; l++)
            {
              if (blob_label[l] == j)
                blob_label[l] = k;
            }
            goto_next_k = 1;
          }
          if (goto_next_k)
            break;
        }
      }
    }
  }

#ifdef DEBUG
  printf("after: %i\n", new_count_labels);
  for (k=1; k<=curr_label; k++)
  {
    printf("%i -> %i\n", k, blob_label[k]);
  }
#endif

  //
  // now we do it all
  //
  // process blob statistics without referring to the image from which
  // the mask was derived
  *npix=mdr_Create(new_count_labels, 1);
  mdr_Zero(*npix);
  *xy_bbox=mdr_Create(new_count_labels, 4);
  mdr_Ones(*xy_bbox);

  int num_replacements = *blob_label;
  for (i=0; ((i<nr) && (num_replacements>0)); i++)
  {
    for (j=0; ((j<nc) && (num_replacements>0)); j++)
    {
      if (Mdr0(x,i,j))
      {
        num_replacements--;

        k = blob_label[ ABS((int) Mdr0(x,i,j)) ];

        // update label: forget about sign of Mdr0(x,i,j)
        Mdr0(x,i,j) = k;

        k--;  // from [1:new_count_labels] to  [0:new_count_labels-1]

        // npix: count pixels
        MdrV0(*npix, k) = MdrV0(*npix, k) + 1;

        // bbox: min i, max i, min j, max j:
        if (MdrV0(*npix, k) == 1)
        {
          Mdr0(*xy_bbox, k, 0) = i+1;
          Mdr0(*xy_bbox, k, 1) = i+1;
          Mdr0(*xy_bbox, k, 2) = j+1;
          Mdr0(*xy_bbox, k, 3) = j+1;
          continue;
        }

        // min i
        if (Mdr0(*xy_bbox, k, 0) > i+1)
          Mdr0(*xy_bbox, k, 0) = i+1;

        // max i :
        if (Mdr0(*xy_bbox, k, 1) < i+1)
          Mdr0(*xy_bbox, k, 1) = i+1;

        // min j:
        if (Mdr0(*xy_bbox, k, 2) > j+1)
          Mdr0(*xy_bbox, k, 2) = j+1;

        // max j :
        if (Mdr0(*xy_bbox, k, 3) < j+1)
          Mdr0(*xy_bbox, k, 3) = j+1;

      }
    }
  }

#ifdef DEBUG
  printf("labels = \n");
  mdr_Print(x, stdout);
#endif

  if (blob_label)
    GC_FREE(blob_label);

  return (new_count_labels);
}

#undef DEBUG

















