/* This file is released as a part of rlabplus project, see rlabplus.sourceforge.net
   Copyright (C) 2006 M. Kostrun

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

//
// gnu_sammon.c:
//    reduce the dimension of a data set using sammon mapping.
//
// References : sprannlib's sammon.c by Dick de Ridder, 18-08-1996
//
// Note:  (M. Kostrun, 2006)
//        1. Changed memory management to Boehm's garbage collector.
//        2. Simplified interface to set and read current map (a caller
//        has to provide initial guess for the map).
//

static double sammon_set_distance (MDR *dset, MDR **set_distance)
{
  int i, k, num_samples=MNR(dset);
  double rval=-1.0;

  // euclidean metrics:
  double (*distance)(int, MDR*, MDR*, MDR*, int, int) =  distance_metric('e');

  if (! *set_distance)
    *set_distance = (MDR *) mdr_Create(num_samples,num_samples);

  if (! *set_distance)
    return rval;

  rval = 0.0;
  for (i=0; i<num_samples; i++)
  {
    Mdr0(*set_distance, i,i) = 0.0;
    for (k = i+1; k < num_samples; k++)
    {
      Mdr0(*set_distance, i,k) = sqrt(distance (1, dset, dset, NULL, i, k));
      Mdr0(*set_distance, k,i) = Mdr0(*set_distance, i,k);
      rval += Mdr0(*set_distance, i,k);
    }
  }

  return rval;
}

/*
 * SAMMON_STRESS: Calculate the Sammon stress betweem
 * 'set' and 'map'. The sammon_init() function has to be
 * called first.
 */
static double sammon_stress (MDR *map, MDR *set_distance, double constant)
{
  int i, j, num_samples=MNR(map);
  double E=0.0, d1, d2;

  // euclidean metrics:
  double (*distance)(int, MDR*, MDR*, MDR*, int, int) =  distance_metric('e');

  for (i=0; i<num_samples-1; i++)
  {
    for (j=i+1; j<num_samples; j++)
    {
      d1  = Mdr0(set_distance,i,j);
      d2  = sqrt(distance (1, map, map, NULL, i, j));
      E  += (d1 - d2) * (d1 - d2) / d1;
    }
  }

  return (E / constant);
}

/*
 * SAMMON_DSTRESS: Calculate the first and the second derivatives of Sammon's
 * stress w.r.t. sample i. Returns two vectors.
 */
static int sammon_dstress (int i, MDR *dE, MDR *d2E, MDR *map, MDR *set_distance, double constant)
{
  int j, k, map_dim=MNC(map), num_samples=MNR(map);
  double d1, d2, d3;

// euclidean metrics:
  double (*distance)(int, MDR*, MDR*, MDR*, int, int) =  distance_metric('e');

 /*
  * Initialize the vectors.
  */

  mdr_Zero(dE);
  mdr_Zero(d2E);

 /*
  * And calculate the first and second derivatives.
  */
  for (j=0; j<num_samples; j++)
  {
    if (j != i)
    {
      if (j > i)
        d1 = Mdr0(set_distance,i,j);
      else
        d1 = Mdr0(set_distance,j,i);

      d2 = sqrt(distance (1, map, map, NULL, i, j));

      for (k=0; k<map_dim; k++)
      {
        d3 = 2.0 * (Mdr0(map,i,k) - Mdr0(map,j,k));
        MdrV0(dE,k)  += (d3) * ((d1 - d2) / (d1 * d2));
        MdrV0(d2E,k) += (1.0 / (d1 * d2)) * ((d1 - d2) - ((0.25 * d3 * d3) / d2)
            * ((1.0 + d1 - d2) / d2));
      }
    }
  }

  for (k = 0; k < map_dim; k++)
  {
    MdrV0(dE,k)  = - MdrV0(dE,k) / constant;
    MdrV0(d2E,k) = - 2.0 * MdrV0(d2E,k) / constant;

//     fprintf(stderr, "sammon_dstress: dE [%i] = %g\n", k, MdrV0(dE,k) );
//     fprintf(stderr, "sammon_dstress: d2E[%i] = %g\n", k, MdrV0(d2E,k) );
  }

  return (RLAB_STATUS_SUCCESS);
}


int sammon_solve (MDR **map, MDR *set_distance, double constant,
                  int method, double learn_rate, double momentum, int max_steps, double min_max_update)
{
  double stress, previous_stress, delta=0, max_update=0;
  int    i, k, map_dim=MNC(*map), num_samples=MNR(*map), ival=RLAB_STATUS_SUCCESS;
  int   loops=0;
  MDR *dE=0, *d2E=0, *previous_update=0, *backup=0;

  dE  = mdr_Create(1, map_dim);
  d2E = mdr_Create(1, map_dim);
  backup = mdr_Create(num_samples,map_dim);

  previous_update = mdr_Create(num_samples, map_dim);
  mdr_Zero(previous_update);

  stress = sammon_stress (*map, set_distance, constant);

  previous_stress = stress + 1.0;

  while ( (stress < previous_stress) && (loops < max_steps) )
  {
    mdr_Copy1IntoPtr2(*map, &backup);

    max_update = 0.0;
    for (i = 0; i < num_samples; i++)
    {
      sammon_dstress (i, dE, d2E, *map, set_distance, constant);

      for (k = 0; k < map_dim; k++)
      {
        switch (method)
        {
          case 0:
          {
            delta = -learn_rate * MdrV0(dE,k);
            break;
          }
          case 1:
          {
            delta = -learn_rate * (MdrV0(dE,k) / abs(MdrV0(d2E,k)));
            break;
          }
        }

        double dummy = delta + momentum * Mdr0(previous_update,i,k);
        Mdr0(*map,i,k) += dummy;
        Mdr0(previous_update,i,k) = dummy;

        max_update = abs(dummy) > max_update ? abs(dummy) : max_update;
      }
    }

    previous_stress = stress;
    stress = sammon_stress (*map, set_distance, constant);
    loops++;

    if (max_update < min_max_update)
      break;
  }

  if (stress > previous_stress)
  {
    //
    // Restore the previous map: we're at the end of the
    // down-hill ride in the gradient descent.
    //
    mdr_Copy1IntoPtr2(backup, map);
  }

  mdr_Destroy(backup);
  mdr_Destroy(dE);
  mdr_Destroy(d2E);
  mdr_Destroy(previous_update);

  return ival;
}
