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
//
//
// F E E D   F O R W A R D   N E U R A L   N E T W O R K S
//
//
//
static double activation_inprod (MDR *y, MDR *p, MDR *q, int k)
{
  register double activation=0.0;

  int j;

  for(j=0; j<SIZE(y); j++)
    activation += Mdr0(p,k,j) * MdrV0(y,j);

  activation += MdrV0(q,k);

  return (activation);
}

static double activation_norm (MDR *y, MDR *p, MDR *q, int k)
{
  register double activation=1.0;

  int j, dim=MNC(y);

  for(j=0; j<dim; j++)
    activation *= MdrV0(y,j) + Mdr0(p,k,j);

  activation += MdrV0(q,k);

  return (activation);
}

static double activation_eucl (MDR *y, MDR *p, MDR *q, int k)
{
  register double activation=1.0;
  int j, dim=MNC(y);
  double tmp;

  for(j=0; j<dim; j++)
  {
    tmp = MdrV0(y,j) - Mdr0(p,k,j);
    activation += tmp * tmp;
  }

  tmp = 1 / MdrV0(q,k);

  activation *= tmp * tmp;

  return sqrt(activation);
}

static double activation_eucl2 (MDR *y, MDR *p, MDR *q, int k)
{
  register double activation=1.0;
  int j, dim=MNC(y);
  double tmp;

  for(j=0; j<dim; j++)
  {
    tmp = MdrV0(y,j) - Mdr0(p,j,k);
    activation += tmp * tmp;
  }

  tmp = 1 / MdrV0(q,k);

  activation *= tmp * tmp;

  return (activation);
}


static double(*activation_function(char act))(MDR *y, MDR *p, MDR *q, int p_col_idx)
{
  switch(act)
  {
    case 'i':
      // inner product: q_k + sum_j y_j * p(N_l x N_{l-1})_j,k
      return &activation_inprod;

    case 'n':
      // norm: q_k + prod_j y_j + p(N_l x N_{l-1})_j,k
      return &activation_norm;

    case '1':
      // euclidean: sqrt( 1/q_k^2 * sum_j (y_j - p(N_l x N_{l-1})_j,k)^2 )
      return &activation_eucl;

    case '2':
      // euclidean: 1/q_k^2 * sum_j (y_j - p(N_l x N_{l-1})_j,k)^2
      return &activation_eucl2;

    default:
      return &activation_inprod;
  }

  return NULL;
}

//
// linear
//
static double transfer_linear (double x)
{
  return (x);
}

static double d_transfer_linear (double x)
{
  return (x);
}

//
// sigmoid
//
static double transfer_sigmoid (double x)
{
  return (1.0 / (1.0 + exp(-x)));
}
static double d_transfer_sigmoid (double x)
{
  double tmp = 1.0 / (1.0 + exp(-x));
  return ((1.0 - tmp) * tmp);
}
static double d2_transfer_sigmoid (double x)
{
  double tmp = transfer_sigmoid(x);
  return (tmp * (1.0 - tmp) * (1.0 - 2.0 * tmp));
}

//
// double sigmoid
//
static double transfer_dblsigm (double x)
{
  return (2.0 * transfer_sigmoid(x) - 1.0);
}
static double d_transfer_dblsigm (double x)
{
  double tmp = 2.0 * transfer_sigmoid(x);
  tmp *= (1.0 - transfer_sigmoid(x));
  return tmp;
}
static double d2_transfer_dblsigm (double x)
{
  double tmp = 2.0 * transfer_sigmoid(x);
  tmp *= (1.0 - transfer_sigmoid(x));
  tmp *= (1.0 - 2.0 * transfer_sigmoid(x));
  return (tmp);
}

//
// tanh
//
static double transfer_tanh (double x)
{
  return tanh(x);
}
static double d_transfer_tanh (double x)
{
  double tmp=tanh(x);
  return (1.0 - tmp*tmp);
}

//
// gauss
//
static double transfer_gauss (double x)
{
  return exp(-0.5 * x * x);
}
static double d_transfer_gauss (double x)
{
  return -1.0 * x * exp(-0.5 * x * x);
}

//
// exponential sigmoid
//
static double transfer_expsigm (double x)
{
  return exp(transfer_sigmoid(x));
}
static double d_transfer_expsigm (double x)
{
  double tmp=transfer_sigmoid(x);
  return (tmp * (1.0 - tmp) * exp(tmp));
}

//
// exponential linear
//
static double transfer_exp (double x)
{
  return exp(x);
}

//
// log sigmoid
//
static double transfer_logsigm (double x)
{
  return log(transfer_sigmoid(x));
}
static double d_transfer_logsigm (double x)
{
  double tmp=transfer_sigmoid(x);
  return (1.0 - tmp);
}

//
// log linear
//
static double transfer_log (double x)
{
  return log(x);
}
static double d_transfer_log (double x)
{
  return 1/x;
}


static double(*transfer_function(char transf))(double x)
{
  switch(transf)
  {
    case '1':
      // unity
      return &transfer_linear;

    case 's':
      // sigmoid
      return &transfer_sigmoid;

    case 'd':
      // double sigmoid
      return &transfer_dblsigm;

    case 't':
      // tanh
      return &transfer_tanh;

    case 'g':
      // gauss
      return &transfer_gauss;

    case 'x':
      // exp sigmoid
      return &transfer_expsigm;

    case 'e':
      // exp
      return &transfer_exp;

    case 'o':
      // log sigmoid
      return &transfer_logsigm;


    default:
      return &transfer_linear;
  }

  return NULL;
}

static double(*d_transfer_function(char transf))(double x)
{
  switch(transf)
  {
    case '1':
      // unity
      return &d_transfer_linear;

    case 's':
      // unitary
      return &d_transfer_sigmoid;

    case 'd':
      // unitary
      return &d_transfer_dblsigm;

    case 't':
      // tanh
      return &d_transfer_tanh;

    case 'g':
      // tanh
      return &d_transfer_gauss;

    case 'x':
      // exp sigmoid
      return &d_transfer_expsigm;

    case 'e':
      // exp
      return &transfer_exp;

    case 'o':
      // log sigmoid
      return &d_transfer_logsigm;

    default:
      return &d_transfer_linear;
  }

  return NULL;
}

static double(*d2_transfer_function(char transf))(double x)
{
  switch(transf)
  {
    case 's':
      // unitary
      return &d2_transfer_sigmoid;

    case 'd':
      // unitary
      return &d2_transfer_dblsigm;

    default:
      return &d2_transfer_sigmoid;
  }

  return NULL;
}

//
static int
    ffnn_randinit_wgt_theta(MDR * Layers, MDR ** NetPtr_wgt, MDR **NetPtr_ifwgt,
                            MDR **NetPtr_theta, MDR **NetPtr_iftheta)
{
  int i,j,k,L=SIZE(Layers);
  int ival=RLAB_STATUS_FAILURE;

  for (i=0; i<(L-1); i++)
  {
    if (! NetPtr_wgt[i])
      goto _exit_ffnn_randinit_wgt_theta;

    int nr=mdiV0(Layers,i);
    int nc=mdiV0(Layers,i+1);

    if (MNR(NetPtr_wgt[i])!=nr || MNC(NetPtr_wgt[i])!=nc)
      goto _exit_ffnn_randinit_wgt_theta;

    //
    // fill-in weights
    //
    if (NetPtr_ifwgt[i])
    {
      if (MNR(NetPtr_ifwgt[i])!=nr || MNC(NetPtr_ifwgt[i])!=nc)
        goto _exit_ffnn_randinit_wgt_theta;

      for (k=0;k<nr; k++) for (j=0; j<nc; j++)
      {
        if (mdi0(NetPtr_ifwgt[i],k,j))
          continue;
        Mdr0(NetPtr_wgt[i],k,j) = 2 * nr * (gslrnguf_() - 0.5);
      }
    }
    else
    {
      // fix_wft is absent: assume no fixed weights
      for (k=0;k<nr; k++)
        for (j=0; j<nc; j++)
      {
        Mdr0(NetPtr_wgt[i],k,j) = 2 * nr * (gslrnguf_() - 0.5);
      }
    }

    //
    // fill-in biases
    //
    if (NetPtr_iftheta[i])
    {
      if (SIZE(NetPtr_iftheta[i])!=nr)
        goto _exit_ffnn_randinit_wgt_theta;

      for (k=0;k<nr; k++)
      {
        if (mdiV0(NetPtr_iftheta[i],k))
          continue;
        MdrV0(NetPtr_theta[i],k) = 2 * nr * (gslrnguf_() - 0.5);
      }
    }
    else
    {
      // fix_theta is absent: assume no fixed biases
      for (k=0;k<nr; k++)
        MdrV0(NetPtr_theta[i],k) = 2 * nr * (gslrnguf_() - 0.5);
    }
  }
  ival=RLAB_STATUS_SUCCESS;

_exit_ffnn_randinit_wgt_theta:

  return ival;
}

 /*
  * Normalize dataset such that the maximum distance in the dataset is
  * equal to the maximum distance in the map.
  * Now, since a network has outputs in [0, 1], the maximum distance in
  * a map of dimension <dim> = sqrt (dim). So, we calculate the maximum
  * distance <max> between samples a and b in the dataset and the mean of
  * a and b. Then, we subtract the mean and multiply by (<dim>/<max>).
  *
  * SCALING_FACTOR is used to further scale the dataset down, to prevent
  * goals of exactly 0 and 1. It should be set to something like 0.5 - 0.9.
  *
  * On the fly we can calculate Sammon's constant (rcp_lambda, Eqn. [16]).
  */
#include "gsl_stat.h"
static double
    samann_scale_dataset (MDR * dset, MDR * dset_dist, MDR ** mu, double *rcp_lambda)
{
  int i, k, num_samples=MNC(dset), nattr=MNR(dset); // dset is transposed!!!!

  // euclidean metrics:
  double (*distance)(int, MDR*, MDR*, MDR*, int, int) =  distance_metric('e');

  // find max distance between two points
  double max_dist=0.0;
  if (dset_dist)
  {
    *rcp_lambda=0;
    for (i=0; i<num_samples; i++)
    {
      Mdr0(dset_dist,i,i) = 0.0;
      for (k = i+1; k < num_samples; k++)
      {
        double tmp_dist = sqrt(distance (0, dset, dset, NULL, i, k)); // dset is transposed: columns are data
        Mdr0(dset_dist,k,i) = Mdr0(dset_dist,i,k) = tmp_dist;
        *rcp_lambda += tmp_dist;
        if (tmp_dist > max_dist)
          max_dist = tmp_dist;
      }
    }
    *rcp_lambda = 1 / (*rcp_lambda);
  }

  // find mean of the dset
  double factor = sqrt ((double) nattr) / max_dist;
  mdr_Avg(0, dset, mu);

  // shift/rescale dataset. don't forge to rescale the distances, as well
  for (i=0; i<num_samples; i++)
  {
    for (k=0; k<nattr; k++)
    {
      Mdr0(dset,k,i) = (Mdr0(dset,k,i) - MdrV0(*mu,k)) * factor;
    }
    for (k=0; k<i; k++)
    {
      Mdr0(dset_dist,k,i) = Mdr0(dset_dist,i,k) =  Mdr0(dset_dist,i,k) * factor;
    }
  }

  return (factor);
}


/*
 * SAMANN_INIT : initialise and set up a network for training with the
 * backpropagation learning rule. Also normalize the dataset.
 */
static double
    samann_init (MDR *set, MDR *set_dist,
                 MDR **mu, MDR *NetVect, MDS **NetPtr_act, MDS **NetPtr_transf,
                 double *rcp_lambda)
{
  double scale_factor=-1.0;
  int i,k,L=SIZE(NetVect);

  //
  // check all activation and transfer functions
  //
  for (i=0; i<(L-1); i++)
  {
    if (! NetPtr_act[i])
      goto _exit_samann_init;
    if (! NetPtr_transf[i])
      goto _exit_samann_init;

    int nr=mdiV0(NetVect,i);

    if (SIZE((MDS *)NetPtr_act[i])!=nr)
      goto _exit_samann_init;

    if (SIZE((MDS *)NetPtr_transf[i])!=nr)
      goto _exit_samann_init;

    for (k=0;k<nr; k++)
    {
      //
      // Check if the network is built up with inner-product or euclidian
      // distance activations and proper transfers (linear-unity or sigmoid)
      //
      if ((MdsV0(NetPtr_act[i],k)[0] != 'i') && (MdsV0(NetPtr_act[i],k)[0] != '1'))
        goto _exit_samann_init;

      if ((MdsV0(NetPtr_transf[i],k)[0] != '1') && (MdsV0(NetPtr_transf[i],k)[0] != 's'))
        goto _exit_samann_init;
    }
  }

  //
  // normalize learning set, calculate lambda and scale
  //
  scale_factor = samann_scale_dataset (set,set_dist,mu,rcp_lambda);

_exit_samann_init:

  return (scale_factor);
}

//
// evaluate feed forward neural network
//
static int  eval_ffnn (MDR *  NetVect,
                       MDR ** NetPtr_wgt, MDR ** NetPtr_theta,
                       MDS ** NetPtr_act, MDS ** NetPtr_transf, MDR **y, double * y_input, double * y_output)
{
  int ival=RLAB_STATUS_FAILURE,l,k;
  void * y0_save=0;

  int L = SIZE(NetVect);
  if (L<2)
    goto _exit_eval_ffnn;

  if (y_input)
    MDPTR(y[L-1]) = (void *) y_input;

  // we keep track of these two:
  if (y_output)
  {
    y0_save = MDPTR(y[0]);
    MDPTR(y[0]) = (void *) y_output;
  }

  for (l=L-2; l>=0; l--)
  {
    for (k=0; k<mdiV0(NetVect,l); k++)
    {
      double (*afunc)(MDR *y, MDR *p, MDR *q, int p_col_idx) = activation_function(MdsV0(NetPtr_act[l],k)[0]);
      if (!afunc)
        goto _exit_eval_ffnn;

      double (*tfunc)(double x) = transfer_function(MdsV0(NetPtr_transf[l],k)[0]);
      if (!tfunc)
        goto _exit_eval_ffnn;

      double tmp = afunc(y[l+1],NetPtr_wgt[l],NetPtr_theta[l],k);

      MdrV0(y[l],k) = tfunc(tmp);
    }
  }
  ival=RLAB_STATUS_SUCCESS;

_exit_eval_ffnn:

  if (y0_save)
  {
    MDPTR(y[0]) = y0_save;
    y0_save = 0;
  }

  return ival;
}

static double
    eval_samman_stress(MDR * dset, MDR * set_distance, MDR *  NetVect,
                       MDR ** NetPtr_wgt, MDR ** NetPtr_theta,
                       MDS ** NetPtr_act, MDS ** NetPtr_transf, MDR ** y_i, MDR ** y_j,
                       double rcp_lambda, double wgt_rcpdist)
{
  double d1, d2, E=0;

  // euclidean metrics:
  double (*distance)(int, MDR*, MDR*, MDR*, int, int) =  distance_metric('e');

  int num_samples=MNC(dset),i,j;

  for (i=0; i<num_samples-1; i++)
  {
    for (j=i+1; j<num_samples; j++)
    {
      // distance between points in initial space
      d1  = Mdr0(set_distance,i,j);

      // find their maps
      eval_ffnn (NetVect, NetPtr_wgt, NetPtr_theta, NetPtr_act, NetPtr_transf, y_i, &Mdr0(dset,0,i), NULL);
      eval_ffnn (NetVect, NetPtr_wgt, NetPtr_theta, NetPtr_act, NetPtr_transf, y_j, &Mdr0(dset,0,j), NULL);

      // then their distances
      d2  = sqrt(distance (0, y_i[0], y_j[0], NULL, 0, 0));
      //E  += (d1 - d2) * (d1 - d2) / d1;
      // test:
      E  += (d1 - d2) * (d1 - d2) * (wgt_rcpdist * d1 + 1/d1);
    }
  }

  return (E * (rcp_lambda));
}


static int
    samann_adapt_net (MDR * NetVect,
                      MDR ** NetPtr_wgt, MDR ** NetPtr_ifwgt, MDR ** NetPtr_theta, MDR ** NetPtr_iftheta,
                      MDS ** NetPtr_act, MDS ** NetPtr_transf, MDR **y1, MDR **y2, MDR **del1, MDR **del2,
                      MDR ** dw,
                      double eta, double alpha, int options,
                      double rcp_lambda, double wgt_rcpdist, double set_dist, double map_dist)
{
  // keep in mind that 'y1' and 'y2' contain all the intermediate nodes values for two input vectors
  double tmp;

  // we iterate over entire network once and propose modification of weights and biases
  // so that the sammon stress is not-increased with respect to two data y1 and y2
  int ival=RLAB_STATUS_FAILURE,l,k,j;
  int L = SIZE(NetVect);
  if (L<2)
    goto _exit_samann_adapt_net;

  //
  // 1. find delta's through backpropagation
  //
  for(l=0; l<L-1; l++)
  {
    if (l==0)
    {
      // output layer deserves special treatment: we are fitting network to sammon stress function
//       tmp = 2.0 * rcp_lambda * ((map_dist - set_dist) / (set_dist * map_dist));
      // test:
      tmp = 2.0 * rcp_lambda * ((map_dist - set_dist) / map_dist) * (wgt_rcpdist * set_dist +  1/set_dist);
    }
    else
    {
      // delta_l = NetPtr_wgt[l-1]^T * delta_{l-1}
      mdr_rgemtv(NetPtr_wgt[l-1], 1.0, del1[l-1], 0, &del1[l]);
      mdr_rgemtv(NetPtr_wgt[l-1], 1.0, del2[l-1], 0, &del2[l]);
    }

    for (k=0; k<mdiV0(NetVect,l); k++)
    {
      double (*afunc)(MDR *y, MDR *p, MDR *q, int p_col_idx) = activation_function(MdsV0(NetPtr_act[l],k)[0]);
      if (!afunc)
        goto _exit_samann_adapt_net;
      double (*tfunc)(double x) = transfer_function(MdsV0(NetPtr_transf[l],k)[0]);
      if (!tfunc)
        goto _exit_samann_adapt_net;
      double (*dtfunc)(double x) = d_transfer_function(MdsV0(NetPtr_transf[l],k)[0]);
      if (!dtfunc)
        goto _exit_samann_adapt_net;

      double a1_k = afunc(y1[l+1],NetPtr_wgt[l],NetPtr_theta[l],k);
      double a2_k = afunc(y2[l+1],NetPtr_wgt[l],NetPtr_theta[l],k);

      if (l==0)
      {
        // output layer deserves special treatment: we are fitting network to sammon stress function
        tmp *= (MdrV0(y1[l],k) - MdrV0(y2[l],k));
        MdrV0(del1[l],k) = tmp * dtfunc(a1_k);
        MdrV0(del2[l],k) = tmp * dtfunc(a2_k);
      }
      else
      {
        MdrV0(del1[l],k) *= dtfunc(a1_k);
        MdrV0(del2[l],k) *= dtfunc(a2_k);
      }

    } // for (k=0; k<mdiV0(NetVect,l); k++)

  } // for(l=0; l<L-1; l++)

  //
  // 2. adapt weights
  //
  for(l=0; l<L-1; l++)
  {
    for (k=0; k<mdiV0(NetVect,l); k++)
    {
      int adapt_theta=1;
      if (NetPtr_iftheta)
        if (NetPtr_iftheta[l])
          if (mdiV0(NetPtr_iftheta[l],k)==1)
            adapt_theta=0;

      if (adapt_theta)
      {
        MdrV0(NetPtr_theta[l],k) -= eta * (MdrV0(del1[l],k) - MdrV0(del2[l],k));
      }

      for (j=0; j<mdiV0(NetVect,l+1); j++)
      {
        int adapt_wgt=1;
        if (NetPtr_ifwgt)
          if (NetPtr_ifwgt[l])
            if (mdiV0(NetPtr_ifwgt[l],k)==1)
              adapt_wgt=0;

        if (adapt_wgt)
        {
          Mdr0(dw[l],k,j) = eta * (MdrV0(del1[l],k) * MdrV0(y1[l+1],j) - MdrV0(del2[l],k) * MdrV0(y2[l+1],j))
              + alpha * Mdr0(dw[l],k,j);
          Mdr0(NetPtr_wgt[l],k,j) -= Mdr0(dw[l],k,j);
        }
      }
    }
  }

  ival=RLAB_STATUS_SUCCESS;

_exit_samann_adapt_net:

  return ival;
}


int
    samann_learn_sample (MDR * dset, MDR *d_dist, MDR * NetVect,
                         MDR ** NetPtr_wgt, MDR ** NetPtr_ifwgt, MDR ** NetPtr_theta, MDR ** NetPtr_iftheta,
                         MDS ** NetPtr_act, MDS ** NetPtr_transf, MDR **y1, MDR **y2, MDR **dy1, MDR **dy2,
                         MDR ** dw,
                         int i1, int i2, double eta, double alpha, int options, double rcp_lambda, double wgt_rcpdist)
{
  int ival=RLAB_STATUS_FAILURE;

  // Set the metric function as indicated by dist
  double (*distance)(int, MDR*, MDR*, MDR*, int, int) =  distance_metric('e');


  //
  // evaluate nn on first sample
  //
  if (eval_ffnn (NetVect, NetPtr_wgt, NetPtr_theta, NetPtr_act, NetPtr_transf, y1,
                &Mdr0(dset,0,i1), NULL))
  {
    goto _exit_samann_learn_sample;
  }

  // evaluate nn on the second sample
  if (eval_ffnn (NetVect, NetPtr_wgt, NetPtr_theta, NetPtr_act, NetPtr_transf, y2,
                &Mdr0(dset,0,i2), NULL) )
  {
    goto _exit_samann_learn_sample;
  }

  double set_dist = Mdr0(d_dist, i1, i2);
  double map_dist = sqrt( distance (0, y1[0], y2[0], NULL, 0, 0) );

  // now do the updates
  if (samann_adapt_net( NetVect,
                        NetPtr_wgt, NetPtr_ifwgt, NetPtr_theta, NetPtr_iftheta,
                        NetPtr_act, NetPtr_transf, y1, y2, dy1, dy2, dw,
                        eta, alpha, options, rcp_lambda, wgt_rcpdist, set_dist, map_dist ) )
  {
    goto _exit_samann_learn_sample;
  }

  ival=RLAB_STATUS_SUCCESS;

_exit_samann_learn_sample:

  return ival;
}


static int
    samann_learn_backprop (MDR * dset, MDR *dset_dist, MDR * NetVect,
                           MDR ** NetPtr_wgt, MDR ** NetPtr_ifwgt, MDR ** NetPtr_theta, MDR ** NetPtr_iftheta,
                           MDR ** NetPtr_act, MDR ** NetPtr_transf,
                           double eta, double alpha, int maxi,
                           int options, double rcp_lambda, double wgt_rcpdist,
                           double * stress, FILE *fptr, MDR * y_out,
                           int ncompstress, double delta_stress_threshold, double stress_threshold)
{
  int i,i1,i2, ival=RLAB_STATUS_FAILURE, i1_2, i2_2, num_attempts=10;
  int num_samples = MNC(dset);
  size_t L=SIZE(NetVect), total_cycles=0;
  double stress_old;
  MDR **y1=NULL, **y2=NULL, **dy1=NULL, **dy2=NULL, **dw=NULL;

  //
  // this is computational space for internal nodes of the network for pairs of values
  //
  if (L<1)
    return ival;

  y1  = GC_MALLOC(L*sizeof(MDR*));
  y2  = GC_MALLOC(L*sizeof(MDR*));
  dy1 = GC_MALLOC(L*sizeof(MDR*));
  dy2 = GC_MALLOC(L*sizeof(MDR*));
  dw  = GC_MALLOC(L*sizeof(MDR*));

  for (i=0; i<L-1; i++)
  {
    y1[i]  = mdr_Create(mdiV0(NetVect,i),1);
    y2[i]  = mdr_Create(mdiV0(NetVect,i),1);
    dy1[i] = mdr_Create(mdiV0(NetVect,i),1);
    dy2[i] = mdr_Create(mdiV0(NetVect,i),1);
    dw[i]  = mdr_Create(mdiV0(NetVect,i), mdiV0(NetVect,i+1)); // for momentum acceleration
    mdr_Zero(dw[i]);
  }
  y1[L-1] = mdr_CreateEmpty(mdiV0(NetVect,L-1),1);  // input layer is special: we will use pointer from dset'
  y2[L-1] = mdr_CreateEmpty(mdiV0(NetVect,L-1),1);  // input layer is special: we will use pointer from dset'
  dy1[L-1] = mdr_Create(mdiV0(NetVect,L-1),1);       //  and layers need not be zeros
  dy2[L-1] = mdr_Create(mdiV0(NetVect,L-1),1);       //  and layers need not be zeros

  if (fptr)
  {
    fprintf (fptr,
             "--------------------------------------------------------------------------\n"
                 "RLaB: Using backpropagation algorithm for samann mapping calculation.\n"
            );
  }

  //
  // let's iterate
  //
  for (total_cycles=0; total_cycles<maxi; total_cycles++)
  {
    for (i1=0; i1<num_samples-1; i1++)
    {
      for (i2=i1+1; i2<num_samples; i2++)
      {
        if (options & RLAB_SAMANN_RANDOM_PAIR)
        {
          // go over two randomly chosen pairs and try to improve
          i2_2 = i1_2 = ((int) ((double)gslrnguf_() * (double) num_samples)) % num_samples;
          while (i2_2 == i1_2)
            i2_2 = ((int) ((double)gslrnguf_() * (double) num_samples)) % num_samples;
        }
        else
        {
          i1_2 = i1;
          i2_2 = i2;
        }

        // do 'num_attempts' iterations per pair of data
        i = 0;
        for(i=0; i<num_attempts; i++)
        {
          if (samann_learn_sample (dset, dset_dist, NetVect,
                                   NetPtr_wgt, NetPtr_ifwgt, NetPtr_theta, NetPtr_iftheta,
                                   NetPtr_act, NetPtr_transf, y1, y2, dy1, dy2, dw,
                                   i1_2, i2_2, eta, alpha, options, rcp_lambda, wgt_rcpdist) )
          {
            goto _exit_samann_learn;
          }
        }
      }
    }

    if (!(total_cycles%ncompstress) && total_cycles>0)
    {
      stress_old = *stress;

      *stress = eval_samman_stress(dset, dset_dist,
                                   NetVect, NetPtr_wgt, NetPtr_theta, NetPtr_act, NetPtr_transf,
                                   y1, y2, rcp_lambda, wgt_rcpdist);

      double delta_stress = ABS(stress_old - *stress);

      if (fptr)
        fprintf(fptr, "samann/backprop: Iteration %6i: Sammon's stress function = %g, delta stress = %g\n",
                total_cycles+1, *stress,  stress_old - *stress);

      if (delta_stress_threshold)
      {
        if (delta_stress < delta_stress_threshold)
          break;
      }

      if (stress_threshold)
      {
        if (*stress <= stress_threshold)
          break;

      }
    }

  } // for (total_cycles=0; total_cycles<maxi; total_cycles++)

  for (i=0;i<num_samples; i++)
  {
    eval_ffnn ( NetVect, NetPtr_wgt, NetPtr_theta, NetPtr_act, NetPtr_transf, y1,
               &Mdr0(dset,0,i), &Mdr0(y_out,0,i) );
  }
  ival = RLAB_STATUS_SUCCESS;

_exit_samann_learn:

  // clean-up computational space
  MDPTR(y1[L-1])=0;
  MDPTR(y2[L-1])=0;
  for (i=0; i<L; i++)
  {
    mdr_Destroy(dy1[i]); dy1[i] = 0;
    mdr_Destroy(dy2[i]); dy2[i] = 0;
    mdr_Destroy(y1[i]); y1[i] = 0;
    mdr_Destroy(y2[i]); y2[i] = 0;
    if (i<L-1)
    {
      mdr_Destroy(dw[i]); dw[i] = 0;
    }
  }
  GC_FREE(y1); GC_FREE(y2); GC_FREE(dy1); GC_FREE(dy2); GC_FREE(dw);
  return ival;
}


//
// proof of concept: use GSL minimizers without derivatives
//
struct RLAB_FFNN_EVAL_SAMANN
{
  MDR * dset;
  MDR * set_distance;
  MDR *  NetVect;
  int * idx_wgt;
  int * idx_the;
  MDR ** NetPtr_wgt;
  MDR ** NetPtr_theta;
  MDS ** NetPtr_act;
  MDS ** NetPtr_transf;
  int n_pairs;
  int  * i1;
  MDR ** y1;
  int  * i2;
  MDR ** y2;
  double rcp_lambda;
  double wgt_rcpdist;
};


#include <gsl/gsl_mode.h>
#include <gsl/gsl_precision.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_machine.h>
#include <gsl/gsl_min.h>
#include <gsl/gsl_multimin.h>

// scalar function of a vector variable
static double mins_gslrlab_f (const gsl_vector * x, void *params)
{
  struct RLAB_FFNN_EVAL_SAMANN * ffnn_samann = (struct RLAB_FFNN_EVAL_SAMANN *) params;
  int * i1 = ffnn_samann->i1;
  int * i2 = ffnn_samann->i2;
  int n_pairs = ffnn_samann->n_pairs;
  int L  = SIZE(ffnn_samann->NetVect);
  int * idx_wgt = ffnn_samann->idx_wgt;
  int * idx_the = ffnn_samann->idx_the;
  double rval = 0;
  int i;
  double d1,d2;

  // Set the metric function as indicated by dist
  double (*distance)(int, MDR*, MDR*, MDR*, int, int) =  distance_metric('e');

  // attach pointers from our NET definition to provided GSL vector with current weights
  for (i=0; i<L-1; i++)
  {
    MDPTR(ffnn_samann->NetPtr_wgt[i]) = &x->data[idx_wgt[i]];
    MDPTR(ffnn_samann->NetPtr_theta[i]) = &x->data[idx_the[i]];
  }

  rval = 0.0;
  for (i=0; i<n_pairs; i++)
  {
    // evaluate NET at two data points:

    // y1
    eval_ffnn ( ffnn_samann->NetVect, ffnn_samann->NetPtr_wgt, ffnn_samann->NetPtr_theta, ffnn_samann->NetPtr_act,
                ffnn_samann->NetPtr_transf, ffnn_samann->y1,
                &Mdr0(ffnn_samann->dset,0,i1[i]), NULL );
    // y2:
    eval_ffnn ( ffnn_samann->NetVect, ffnn_samann->NetPtr_wgt, ffnn_samann->NetPtr_theta, ffnn_samann->NetPtr_act,
                ffnn_samann->NetPtr_transf, ffnn_samann->y2,
                &Mdr0(ffnn_samann->dset,0,i2[i]), NULL );

    // distance between points in initial space
    d1  = Mdr0(ffnn_samann->set_distance,i1[i],i2[i]);

    // then their distances in mapped space
    d2  = sqrt(distance (0, ffnn_samann->y1[0], ffnn_samann->y2[0], NULL, 0, 0));

    rval += (d1 - d2) * (d1 - d2) * (ffnn_samann->wgt_rcpdist * d1 + 1/d1);
  }

//   for (i=0; i<x->size; i++)
//   {
//     printf ("  %g", gsl_vector_get (x, i));
//   }
//   printf (":%g\n", rval);


  return rval;
}



static int
    samann_learn_simplex (MDR * dset, MDR *dset_dist, MDR * NetVect,
                          MDR ** NetPtr_wgt, MDR ** NetPtr_ifwgt, MDR ** NetPtr_theta, MDR ** NetPtr_iftheta,
                          MDR ** NetPtr_act, MDR ** NetPtr_transf,
                          double s0, double abserr, int num_attempts, int maxi,
                          int options, double rcp_lambda, double wgt_rcpdist,
                          double * stress, FILE *fptr, MDR * y_out,
                          int ncompstress, double delta_stress_threshold, double stress_threshold)
{
  int i,k,i1,i2,ival=RLAB_STATUS_FAILURE;
  int num_samples = MNC(dset);
  int L=SIZE(NetVect), total_cycles=0;
  double stress_old;
  int batch_size = 4;

  MDR **NetPtr_wgt2, **NetPtr_theta2, **y1, **y2;


  //
  // this is computational space for internal nodes of the network for pairs of values
  //
  NetPtr_wgt2   = GC_MALLOC(L*sizeof(MDR*));
  NetPtr_theta2 = GC_MALLOC(L*sizeof(MDR*));
  y1 = GC_MALLOC(L*sizeof(MDR*));
  y2 = GC_MALLOC(L*sizeof(MDR*));

  int num_weights = 0;
  int num_thetas  = 0;
  for (i=0; i<L-1; i++)
  {
    //
    num_thetas  += mdiV0(NetVect,i);
    num_weights += mdiV0(NetVect,i) * mdiV0(NetVect,i+1);
    NetPtr_wgt2[i] = mdr_CreateEmpty(mdiV0(NetVect,i),mdiV0(NetVect,i+1));
    NetPtr_theta2[i] = mdr_CreateEmpty(mdiV0(NetVect,i),1);
    //
    y1[i]  = mdr_Create(mdiV0(NetVect,i),1);
    y2[i]  = mdr_Create(mdiV0(NetVect,i),1);
  }
  y1[L-1] = mdr_CreateEmpty(mdiV0(NetVect,L-1),1);  // input layer is special: we will use pointer from dset'
  y2[L-1] = mdr_CreateEmpty(mdiV0(NetVect,L-1),1);  // input layer is special: we will use pointer from dset'

  // these are 1-D arrays that will be used for f/df minimizers
  int cmx_dimx = num_weights + num_thetas;
  gsl_vector *ss = gsl_vector_alloc (cmx_dimx);
  gsl_vector *x  = gsl_vector_alloc (cmx_dimx);
  int * idx_wgt = GC_MALLOC(L*sizeof(int));
  int * idx_the = GC_MALLOC(L*sizeof(int));
  int * i1_pairs = GC_MALLOC(batch_size*sizeof(int));
  int * i2_pairs = GC_MALLOC(batch_size*sizeof(int));

  //
  // now copy wgt and thetas to p
  //
  i1=0;
  for (i=0; i<L-1; i++)
  {
    idx_wgt[i] = i1;
    int size = SIZE(NetPtr_wgt[i]);
    memcpy(&x->data[i1],MDPTR(NetPtr_wgt[i]),size*sizeof(double));
    for (k=0; k<size; k++)
    {
      int change_it = 1;
      if (NetPtr_ifwgt)
        if (NetPtr_ifwgt[i])
          if (mdiV0(NetPtr_ifwgt[i],k)==1)
            change_it = 0;
      if (change_it)
        gsl_vector_set (ss, i1+k, s0);
      else
        gsl_vector_set (ss, i1+k, 0.0);
    }
    i1 += size;
  }
  for (i=0; i<L-1; i++)
  {
    idx_the[i] = i1;
    int size = SIZE(NetPtr_theta[i]);
    memcpy(&x->data[i1],MDPTR(NetPtr_theta[i]),size*sizeof(double));
    for (k=0; k<size; k++)
    {
      int change_it = 1;
      if (NetPtr_iftheta)
        if (NetPtr_iftheta[i])
          if (mdiV0(NetPtr_iftheta[i],k)==1)
            change_it = 0;
      if (change_it)
        gsl_vector_set (ss, i1+k, s0);
      else
        gsl_vector_set (ss, i1+k, 0.0);
    }
    i1 += size;
  }

  // now wrap all of that for passing it as a pointer to the gsl solver
  // for the first call select only one pair:
  i1_pairs[0] = 0;
  i2_pairs[0] = 1;
  struct RLAB_FFNN_EVAL_SAMANN ffnn_samann;
  ffnn_samann.dset = dset;
  ffnn_samann.set_distance = dset_dist;
  ffnn_samann.NetVect = NetVect;
  ffnn_samann.NetPtr_wgt = NetPtr_wgt2;
  ffnn_samann.NetPtr_theta = NetPtr_theta2;
  ffnn_samann.NetPtr_act = NetPtr_act;
  ffnn_samann.NetPtr_transf = NetPtr_transf;
  ffnn_samann.idx_wgt = idx_wgt;
  ffnn_samann.idx_the = idx_the;
  ffnn_samann.n_pairs = 1;
  ffnn_samann.i1 = i1_pairs;
  ffnn_samann.i2 = i2_pairs;
  ffnn_samann.y1 = y1;
  ffnn_samann.y2 = y2;
  ffnn_samann.rcp_lambda = rcp_lambda;
  ffnn_samann.wgt_rcpdist = wgt_rcpdist;
  //
  // using GSL SIMPLEX solver
  //
  //
  gsl_set_error_handler_off();  // turn-off error handling inside GSL

  const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex;

  int status, i1_2, i2_2;

  //
  gsl_multimin_function f =
  {
    &mins_gslrlab_f,
    cmx_dimx,
    (void *) &ffnn_samann
  };

  //
  // do not write any run-time messages
  //
  if (fptr)
  {
    fprintf (fptr,
             "--------------------------------------------------------------------------\n"
                 "RLaB: Using simplex minimizer for samann mapping calculation.\n"
            );
  }

  //
  // let's iterate
  //
  int pair_count = 0;
  for (total_cycles=0; total_cycles<maxi; total_cycles++)
  {
    for (i1=0; i1<num_samples-1; i1++)
    {
      for (i2=i1+1; i2<num_samples; i2++)
      {
        if (pair_count < batch_size)
        {
          if (options & RLAB_SAMANN_RANDOM_PAIR)
          {
           // go over two randomly chosen pairs and try to improve
            i2_2 = i1_2 = ((int) ((double)gslrnguf_() * (double) num_samples)) % num_samples;
            while (i2_2 == i1_2)
              i2_2 = ((int) ((double)gslrnguf_() * (double) num_samples)) % num_samples;
          }
          else
          {
          // go over all pairs
            i1_2 = i1;
            i2_2 = i2;
          }
          i1_pairs[pair_count] = i1_2;
          i2_pairs[pair_count] = i2_2;
          pair_count++;
        }
        else
        {
          //
          // for each set of points we have to do minimization all over again:
          //  cannot reuse previous structure
          //
          pair_count = 0;

          gsl_multimin_fminimizer *s = gsl_multimin_fminimizer_alloc (T, cmx_dimx);
          gsl_multimin_fminimizer_set (s, &f, x, ss);

          // do 'num_attempts' iterations per pair of data
          for (i=0, status=0; (i<num_attempts) && (status==0); i++)
          {
            status = gsl_multimin_fminimizer_iterate (s);
            if (status)
              break;
            status = gsl_multimin_test_size (gsl_multimin_fminimizer_size (s), abserr);
          }
          gsl_vector_memcpy (x, s->x);
          gsl_multimin_fminimizer_free (s);
        }

      } // for (i2=i1+1; i2<num_samples; i2++)
    } // for (i1=0; i1<num_samples-1; i1++)

    if (!(total_cycles%ncompstress))
    {
      stress_old = *stress;

      for (i=0; i<L-1; i++)
      {
        MDPTR(NetPtr_wgt2[i]) = &x->data[idx_wgt[i]];
        MDPTR(NetPtr_theta2[i]) = &x->data[idx_the[i]];
      }
      *stress = eval_samman_stress(dset, dset_dist,
                                   NetVect, NetPtr_wgt2, NetPtr_theta2, NetPtr_act, NetPtr_transf,
                                   y1, y2, rcp_lambda, wgt_rcpdist);

      if (fptr)
        fprintf(fptr, "samann/simplex: Iteration %6i: Sammon's stress function = %g, delta stress = %g\n",
                total_cycles+1, *stress,  stress_old- *stress);

      if (delta_stress_threshold)
        if (ABS(stress_old - *stress) < delta_stress_threshold)
          break;

      if (stress_threshold)
        if (*stress <= stress_threshold)
          break;

    }

  } // for (total_cycles=0; total_cycles<maxi; total_cycles++)

  //
  // copy data from computational array back to book-keeping array
  //
  int size;
  for (i=0; i<L-1; i++)
  {
    size = SIZE(NetPtr_wgt[i]);
    if (size>0)
      memcpy(MDPTR(NetPtr_wgt[i]), &x->data[idx_wgt[i]], size * sizeof(double));
    size = SIZE(NetPtr_theta[i]);
    if (size>0)
      memcpy(MDPTR(NetPtr_theta[i]), &x->data[idx_the[i]], size * sizeof(double));
  }

  for (i=0;i<num_samples; i++)
  {
    eval_ffnn ( NetVect, NetPtr_wgt, NetPtr_theta, NetPtr_act, NetPtr_transf, y1,
                &Mdr0(dset,0,i), &Mdr0(y_out,0,i) );
  }

  gsl_vector_free (ss);
  gsl_vector_free (x);


  // clean-up computational space
  MDPTR(y1[L-1])=0;
  MDPTR(y2[L-1])=0;
  for (i=0; i<L; i++)
  {
    mdr_Destroy(y1[i]); y1[i] = 0;
    mdr_Destroy(y2[i]); y2[i] = 0;
  }
  GC_FREE(y1); GC_FREE(y2);

  // cleanup fake arrays
  for (i=0; i<L-1; i++)
  {
    MDPTR(NetPtr_wgt2[i]) = 0;
    mdr_Destroy(NetPtr_wgt2[i]);
    NetPtr_wgt2[i]=0;
    MDPTR(NetPtr_theta2[i]) = 0;
    mdr_Destroy(NetPtr_theta2[i]);
    NetPtr_theta2[i]=0;
  }
  GC_FREE(NetPtr_wgt2); GC_FREE(NetPtr_theta2);
  GC_FREE(idx_wgt); GC_FREE(idx_the);
  GC_FREE(i1_pairs); GC_FREE(i2_pairs);

  return ival;
}

//
// this function is used for passing pointers to the fortran function
//
static struct RLAB_FFNN_EVAL_SAMANN * newuoa_getptr( struct RLAB_FFNN_EVAL_SAMANN * netparams)
{
  static struct RLAB_FFNN_EVAL_SAMANN * ptr = 0;

  if (netparams)
    ptr = netparams;

  return ptr;
}

// scalar function of a vector variable
static int mins_newuoa_f (int * idummy, double * xval, double * f)
{
  struct RLAB_FFNN_EVAL_SAMANN * ffnn_samann = newuoa_getptr(NULL);
  if (!ffnn_samann)
    return 0;

  int n_pairs = ffnn_samann->n_pairs;
  int * i1 = ffnn_samann->i1;
  int * i2 = ffnn_samann->i2;
  int L  = SIZE(ffnn_samann->NetVect);
  int * idx_wgt = ffnn_samann->idx_wgt;
  int * idx_the = ffnn_samann->idx_the;
  double d1, d2;
  int i;

  // Set the metric function as indicated by dist
  double (*distance)(int, MDR*, MDR*, MDR*, int, int) =  distance_metric('e');

  // attach pointers from our NET definition to provided xval vector with current weights
  for (i=0; i<L-1; i++)
  {
    MDPTR(ffnn_samann->NetPtr_wgt[i])   = &xval[idx_wgt[i]];
    MDPTR(ffnn_samann->NetPtr_theta[i]) = &xval[idx_the[i]];
  }

  *f = 0.0;
  for (i=0; i<n_pairs; i++)
  {
    // evaluate NET at two data points:
    // y1
    eval_ffnn ( ffnn_samann->NetVect, ffnn_samann->NetPtr_wgt, ffnn_samann->NetPtr_theta, ffnn_samann->NetPtr_act,
                ffnn_samann->NetPtr_transf, ffnn_samann->y1,
                &Mdr0(ffnn_samann->dset,0,i1[i]), NULL );
    // y2:
    eval_ffnn ( ffnn_samann->NetVect, ffnn_samann->NetPtr_wgt, ffnn_samann->NetPtr_theta, ffnn_samann->NetPtr_act,
                ffnn_samann->NetPtr_transf, ffnn_samann->y2,
                &Mdr0(ffnn_samann->dset,0,i2[i]), NULL );

    // distance between points in initial space
    d1  = Mdr0(ffnn_samann->set_distance,i1[i],i2[i]);

    // then their distances
    d2  = sqrt(distance (0, ffnn_samann->y1[0], ffnn_samann->y2[0], NULL, 0, 0));

    *f += (d1 - d2) * (d1 - d2) * (ffnn_samann->wgt_rcpdist * d1 + 1/d1);
  }

  return 1;
}

#include "libmjdpowell.h"
static int
    samann_learn_newuoa (MDR * dset, MDR *dset_dist, MDR * NetVect,
                         MDR ** NetPtr_wgt, MDR ** NetPtr_ifwgt, MDR ** NetPtr_theta, MDR ** NetPtr_iftheta,
                         MDR ** NetPtr_act, MDR ** NetPtr_transf,
                         double rhobeg, double rhoend, int num_attempts, int maxi,
                         int options, double rcp_lambda, double wgt_rcpdist,
                         double * stress, FILE *fptr, MDR * y_out,
                         int ncompstress, double delta_stress_threshold, double stress_threshold)
{
  int i,i1,i2,ival=RLAB_STATUS_FAILURE;
  int num_samples = MNC(dset);
  int L=SIZE(NetVect), total_cycles=0;
  double stress_old;
  int i1_2, i2_2, iprint=0;
  int batch_size = 2;

  //
  // this is computational space for internal nodes of the network for pairs of values
  //
  MDR ** NetPtr_wgt2   = GC_MALLOC(L*sizeof(MDR*));
  MDR ** NetPtr_theta2 = GC_MALLOC(L*sizeof(MDR*));
  MDR ** y1 = GC_MALLOC(L*sizeof(MDR*));
  MDR ** y2 = GC_MALLOC(L*sizeof(MDR*));
  int num_weights = 0;
  int num_thetas  = 0;
  for (i=0; i<L-1; i++)
  {
    //
    num_thetas  += mdiV0(NetVect,i);
    num_weights += mdiV0(NetVect,i) * mdiV0(NetVect,i+1);
    NetPtr_wgt2[i] = mdr_CreateEmpty(mdiV0(NetVect,i),mdiV0(NetVect,i+1));
    NetPtr_theta2[i] = mdr_CreateEmpty(mdiV0(NetVect,i),1);
    //
    y1[i]  = mdr_Create(mdiV0(NetVect,i),1);
    y2[i]  = mdr_Create(mdiV0(NetVect,i),1);
  }
  y1[L-1] = mdr_CreateEmpty(mdiV0(NetVect,L-1),1);  // input layer is special: we will use pointer from dset'
  y2[L-1] = mdr_CreateEmpty(mdiV0(NetVect,L-1),1);  // input layer is special: we will use pointer from dset'

  // these are 1-D arrays that will be used for f/df minimizers
  int cmx_dimx  = num_weights + num_thetas;
  int npt       = cmx_dimx+num_attempts; // min, to 0.5 * cmx_dimx * (cmx_dimx+1) as max
  double *x  = GC_MALLOC(cmx_dimx * sizeof(double));
  double *wk = GC_MALLOC(((npt+13)*(npt+cmx_dimx)+3*cmx_dimx*(cmx_dimx+3)/2+1)*sizeof(double));
  int * idx_wgt = GC_MALLOC(L*sizeof(int));
  int * idx_the = GC_MALLOC(L*sizeof(int));
  int * i1_pairs = GC_MALLOC(batch_size*sizeof(int));
  int * i2_pairs = GC_MALLOC(batch_size*sizeof(int));

  //
  // now copy wgt and thetas to p
  //
  int size;
  i1=0;
  for (i=0; i<L-1; i++)
  {
    idx_wgt[i] = i1;
    size = SIZE(NetPtr_wgt[i]);
    if (size>0)
      memcpy((void*) &x[i1], (const void *) MDPTR(NetPtr_wgt[i]), (size_t) size*sizeof(double));
    i1 += size;
  }
  
  for (i=0; i<L-1; i++)
  {
    idx_the[i] = i1;
    size = SIZE(NetPtr_theta[i]);
    if (size>0)
      memcpy(&x[i1],MDPTR(NetPtr_theta[i]),size*sizeof(double));
    i1 += size;
  }

  // now wrap all of that for passing it as a pointer to the gsl solver
  struct RLAB_FFNN_EVAL_SAMANN ffnn_samann;
  ffnn_samann.dset = dset;
  ffnn_samann.set_distance = dset_dist;
  ffnn_samann.NetVect = NetVect;
  ffnn_samann.NetPtr_wgt = NetPtr_wgt2;
  ffnn_samann.NetPtr_theta = NetPtr_theta2;
  ffnn_samann.NetPtr_act = NetPtr_act;
  ffnn_samann.NetPtr_transf = NetPtr_transf;
  ffnn_samann.idx_wgt = idx_wgt;
  ffnn_samann.idx_the = idx_the;
  ffnn_samann.n_pairs =  batch_size;
  ffnn_samann.i1 = i1_pairs;
  ffnn_samann.i2 = i2_pairs;
  ffnn_samann.y1 = y1;
  ffnn_samann.y2 = y2;
  ffnn_samann.rcp_lambda = rcp_lambda;
  ffnn_samann.wgt_rcpdist = wgt_rcpdist;

  //
  // using NEWUOA solver
  //
  //
  newuoa_getptr(&ffnn_samann);

  //
  // do not write any run-time messages
  //
  if (fptr)
  {
    fprintf (fptr,
             "--------------------------------------------------------------------------\n"
                 "RLaB: Using NEWUOA simplex minimizer for samann mapping calculation.\n"
            );
  }

  //
  // let's iterate
  //
  int pair_count = 0;
  for (total_cycles=0; total_cycles<maxi; total_cycles++)
  {
    for (i1=0; i1<num_samples-1; i1++)
    {
      for (i2=i1+1; i2<num_samples; i2++)
      {
        if (pair_count < batch_size)
        {
          if (options & RLAB_SAMANN_RANDOM_PAIR)
          {
           // go over two randomly chosen pairs and try to improve
            i2_2 = i1_2 = ((int) ((double)gslrnguf_() * (double) num_samples)) % num_samples;
            while (i2_2 == i1_2)
              i2_2 = ((int) ((double)gslrnguf_() * (double) num_samples)) % num_samples;
          }
          else
          {
          // go over all pairs
            i1_2 = i1;
            i2_2 = i2;
          }
          i1_pairs[pair_count] = i1_2;
          i2_pairs[pair_count] = i2_2;
          pair_count++;
        }
        else
        {
          pair_count = 0;

          // do 'num_attempts' iterations per pair of data
          iprint=0;
          NEWUOA (&cmx_dimx, &npt, x, &rhobeg, &rhoend, &iprint, &num_attempts, wk, mins_newuoa_f, NULL);
        }

      } // for (i2=i1+1; i2<num_samples; i2++)
    } // for (i1=0; i1<num_samples-1; i1++)

    if (!(total_cycles%ncompstress))
    {
      stress_old = *stress;

      for (i=0; i<L-1; i++)
      {
        MDPTR(NetPtr_wgt2[i]) = &x[idx_wgt[i]];
        MDPTR(NetPtr_theta2[i]) = &x[idx_the[i]];
      }

      *stress = eval_samman_stress(dset, dset_dist,
                                   NetVect, NetPtr_wgt2, NetPtr_theta2, NetPtr_act, NetPtr_transf,
                                   y1, y2, rcp_lambda, wgt_rcpdist);

      if (fptr)
        fprintf(fptr, "samann/newuoa: Iteration %6i: Sammon's stress function = %g, delta stress = %g\n",
                total_cycles+1, *stress,  stress_old- *stress);

      if (delta_stress_threshold)
        if (ABS(stress_old - *stress) < delta_stress_threshold)
          break;

      if (stress_threshold)
        if (*stress <= stress_threshold)
          break;

    }

  } // for (total_cycles=0; total_cycles<maxi; total_cycles++)

  //
  // copy data from computational array back to book-keeping array
  //
  for (i=0; i<L-1; i++)
  {
    size = SIZE(NetPtr_wgt[i]);
    if (size>0)
      memcpy(MDPTR(NetPtr_wgt[i]), &x[idx_wgt[i]], size * sizeof(double));

    size = SIZE(NetPtr_theta[i]);
    if (size>0)
      memcpy(MDPTR(NetPtr_theta[i]), &x[idx_the[i]], size * sizeof(double));
  }

  for (i=0;i<num_samples; i++)
  {
    eval_ffnn ( NetVect, NetPtr_wgt, NetPtr_theta, NetPtr_act, NetPtr_transf, y1,
                &Mdr0(dset,0,i), &Mdr0(y_out,0,i) );
  }


  // clean-up computational space
  GC_FREE (wk);
  MDPTR(y1[L-1])=0;
  MDPTR(y2[L-1])=0;
  for (i=0; i<L; i++)
  {
    mdr_Destroy(y1[i]); y1[i] = 0;
    mdr_Destroy(y2[i]); y2[i] = 0;
  }
  GC_FREE(y1); GC_FREE(y2);

  // cleanup fake arrays
  for (i=0; i<L-1; i++)
  {
    MDPTR(NetPtr_wgt2[i]) = 0;
    mdr_Destroy(NetPtr_wgt2[i]);
    NetPtr_wgt2[i]=0;
    MDPTR(NetPtr_theta2[i]) = 0;
    mdr_Destroy(NetPtr_theta2[i]);
    NetPtr_theta2[i]=0;
  }
  GC_FREE(NetPtr_wgt2); GC_FREE(NetPtr_theta2);
  GC_FREE(idx_wgt); GC_FREE(idx_the);
  GC_FREE(i1_pairs); GC_FREE(i2_pairs);
  return ival;
}


//
//
//
//
//
//
static double
    eval_ffnn_mse(MDR * dset, MDR * tset, MDR *  NetVect,
                  MDR ** NetPtr_wgt, MDR ** NetPtr_theta,
                  MDS ** NetPtr_act, MDS ** NetPtr_transf, MDR ** y_i)
{
  double d, E=0;

  // euclidean metrics:
  double (*distance)(int, MDR*, MDR*, MDR*, int, int) =  distance_metric('e');

  int num_samples=MNC(dset),i;

  for (i=0; i<num_samples-1; i++)
  {
    // find nn map of point
    eval_ffnn (NetVect, NetPtr_wgt, NetPtr_theta, NetPtr_act, NetPtr_transf, y_i, &Mdr0(dset,0,i), NULL);

    // then the distances between target and map
    d  = distance (0, y_i[0], tset, NULL, 0, i);
    E += d;
  }

  E = sqrt(E);
  return (E);
}


static int
    ffnn_adapt_net (MDR *tset, int i1,
                    MDR * NetVect, MDR ** NetPtr_wgt, MDR ** NetPtr_ifwgt, MDR ** NetPtr_theta, MDR ** NetPtr_iftheta,
                    MDS ** NetPtr_act, MDS ** NetPtr_transf,
                    MDR **y1, MDR **del1, MDR ** dw,
                    double eta, double alpha, int options )
{
  // keep in mind that 'y1' and 'y2' contain all the intermediate nodes values for two input vectors
  double tmp;

  // we iterate over entire network once and propose modification of weights and biases
  // so that the sammon stress is not-increased with respect to two data y1 and y2
  int ival=RLAB_STATUS_FAILURE,l,k,j;
  int L = SIZE(NetVect);
  if (L<2)
    goto _exit_ffnn_adapt_net;

  //
  // 1. find delta's through backpropagation
  //
  for(l=0; l<L-1; l++)
  {
    if (l==0)
    {
      // output layer deserves special treatment: we are fitting network to sammon stress function
      tmp = 2.0;
    }
    else
    {
      // delta_l = NetPtr_wgt[l-1]^T * delta_{l-1}
      mdr_rgemtv(NetPtr_wgt[l-1], 1.0, del1[l-1], 0, &del1[l]);
    }

    for (k=0; k<mdiV0(NetVect,l); k++)
    {
      double (*afunc)(MDR *y, MDR *p, MDR *q, int p_col_idx) = activation_function(MdsV0(NetPtr_act[l],k)[0]);
      if (!afunc)
        goto _exit_ffnn_adapt_net;
      double (*tfunc)(double x) = transfer_function(MdsV0(NetPtr_transf[l],k)[0]);
      if (!tfunc)
        goto _exit_ffnn_adapt_net;
      double (*dtfunc)(double x) = d_transfer_function(MdsV0(NetPtr_transf[l],k)[0]);
      if (!dtfunc)
        goto _exit_ffnn_adapt_net;

      double a1_k = afunc(y1[l+1],NetPtr_wgt[l],NetPtr_theta[l],k);

      if (l==0)
      {
        // in the output layer: we are minimizing MSE
        tmp *= (MdrV0(y1[l],k) - Mdr0(tset,k,i1));
        MdrV0(del1[l],k) = tmp * dtfunc(a1_k);
      }
      else
      {
        MdrV0(del1[l],k) *= dtfunc(a1_k);
      }

    } // for (k=0; k<mdiV0(NetVect,l); k++)

  } // for(l=0; l<L-1; l++)

  //
  // 2. adapt weights
  //
  for(l=0; l<L-1; l++)
  {
    for (k=0; k<mdiV0(NetVect,l); k++)
    {
      int adapt_theta=1;
      if (NetPtr_iftheta)
        if (NetPtr_iftheta[l])
          if (mdiV0(NetPtr_iftheta[l],k)==1)
            adapt_theta=0;

      if (adapt_theta)
      {
        MdrV0(NetPtr_theta[l],k) -= eta * MdrV0(del1[l],k);
      }

      for (j=0; j<mdiV0(NetVect,l+1); j++)
      {
        int adapt_wgt=1;
        if (NetPtr_ifwgt)
          if (NetPtr_ifwgt[l])
            if (mdiV0(NetPtr_ifwgt[l],k)==1)
              adapt_wgt=0;

        if (adapt_wgt)
        {
          Mdr0(dw[l],k,j) = eta * MdrV0(del1[l],k) * MdrV0(y1[l+1],j) + alpha * Mdr0(dw[l],k,j);
          Mdr0(NetPtr_wgt[l],k,j) -= Mdr0(dw[l],k,j);
        }
      }
    }
  }

  ival=RLAB_STATUS_SUCCESS;

_exit_ffnn_adapt_net:

  return ival;
}



int
    ffnn_learn_sample (MDR * dset, MDR *tset,
                       MDR * NetVect, MDR ** NetPtr_wgt, MDR ** NetPtr_ifwgt, MDR ** NetPtr_theta, MDR ** NetPtr_iftheta,
                       MDS ** NetPtr_act, MDS ** NetPtr_transf,
                       MDR **y1, MDR **dy1, MDR ** dw,
                       int i1, double eta, double alpha, int options)
{
  int ival=RLAB_STATUS_FAILURE;

  //
  // evaluate nn on first sample
  //
  if (eval_ffnn (NetVect, NetPtr_wgt, NetPtr_theta, NetPtr_act, NetPtr_transf, y1,
                &Mdr0(dset,0,i1), NULL))
  {
    goto _exit_ffnn_learn_sample;
  }

  // now do the updates
  if (ffnn_adapt_net( tset, i1,
                      NetVect, NetPtr_wgt, NetPtr_ifwgt, NetPtr_theta, NetPtr_iftheta,
                      NetPtr_act, NetPtr_transf, y1, dy1, dw,
                      eta, alpha, options ) )
  {
    goto _exit_ffnn_learn_sample;
  }

  ival=RLAB_STATUS_SUCCESS;

_exit_ffnn_learn_sample:

  return ival;
}


static int
    ffnn_learn_backprop (MDR * dset, MDR * tset, MDR * NetVect,
                         MDR ** NetPtr_wgt, MDR ** NetPtr_ifwgt, MDR ** NetPtr_theta, MDR ** NetPtr_iftheta,
                         MDR ** NetPtr_act, MDR ** NetPtr_transf,
                         double eta, double alpha, int maxi,
                         int options, double * mse, FILE *fptr, MDR * y_out,
                         int ncompstress, double delta_mse_threshold, double mse_threshold)
{
  int i,i1,ival=RLAB_STATUS_FAILURE, i1_2, num_attempts=10;
  int num_samples = MNC(dset);
  int L=SIZE(NetVect), total_cycles=0;
  double mse_old;

  //
  // this is computational space for internal nodes of the network for pairs of values
  //
  MDR **y1  = (MDR **) GC_MALLOC(L*sizeof(MDR*));
  MDR **dy1 = (MDR **) GC_MALLOC(L*sizeof(MDR*));
  MDR **dw  = (MDR **) GC_MALLOC(L*sizeof(MDR*));
  for (i=0; i<L-1; i++)
  {
    y1[i]  = mdr_Create(mdiV0(NetVect,i),1);
    dy1[i] = mdr_Create(mdiV0(NetVect,i),1);
    dw[i]  = mdr_Create(mdiV0(NetVect,i), mdiV0(NetVect,i+1)); // for momentum acceleration
    mdr_Zero(dw[i]);
  }
  y1[L-1] = mdr_CreateEmpty(mdiV0(NetVect,L-1),1);  // input layer is special: we will use pointer from dset'
  dy1[L-1] = mdr_Create(mdiV0(NetVect,L-1),1);       //  and layers need not be zeros

  if (fptr)
  {
    fprintf (fptr,
             "--------------------------------------------------------------------------\n"
                 "RLaB: Using backpropagation algorithm for feed-forward neural network training.\n"
            );
  }

  //
  // let's iterate
  //
  for (total_cycles=0; total_cycles<maxi; total_cycles++)
  {

    for (i1=0; i1<num_samples; i1++)
    {

      if (options & RLAB_SAMANN_RANDOM_PAIR)
      {
        // choose random target of improvement
        i1_2 = ((int) ((double)gslrnguf_() * (double) num_samples)) % num_samples;
      }
      else
      {
        i1_2 = i1;
      }

      // do 'num_attempts' iterations per pair of data
      i = 0;
      for(i=0; i<num_attempts; i++)
      {
        if (ffnn_learn_sample (dset, tset,
                               NetVect, NetPtr_wgt, NetPtr_ifwgt, NetPtr_theta, NetPtr_iftheta,
                               NetPtr_act, NetPtr_transf,
                               y1, dy1, dw,
                               i1_2, eta, alpha, options) )
        {
          goto _exit_ffnn_learn;
        }
      }

    } // for (i1=0; i1<num_samples; i1++)

    if (!(total_cycles%ncompstress))
    {
      mse_old = *mse;

      *mse = eval_ffnn_mse(dset, tset, NetVect, NetPtr_wgt, NetPtr_theta, NetPtr_act, NetPtr_transf, y1);

      if (fptr)
        fprintf(fptr, "ffnn/backprop: Iteration %6i: training mean-square-error = %g, delta stress = %g\n",
                total_cycles+1, *mse,  mse_old - *mse);

      if (delta_mse_threshold)
        if (mse_old - *mse < delta_mse_threshold)
          break;

      if (mse_threshold)
        if (*mse <= mse_threshold)
          break;
    }

  } // for (total_cycles=0; total_cycles<maxi; total_cycles++)

  for (i=0;i<num_samples; i++)
  {
    eval_ffnn ( NetVect, NetPtr_wgt, NetPtr_theta, NetPtr_act, NetPtr_transf, y1,
               &Mdr0(dset,0,i), &Mdr0(y_out,0,i) );
  }
  ival = RLAB_STATUS_SUCCESS;

_exit_ffnn_learn:

  // clean-up computational space
  MDPTR(y1[L-1])=0;
  for (i=0; i<L; i++)
  {
    mdr_Destroy(dy1[i]); dy1[i] = 0;
    mdr_Destroy(y1[i]); y1[i] = 0;
    if (i<L-1)
    {
      mdr_Destroy(dw[i]); dw[i] = 0;
    }
  }
  GC_FREE(y1); GC_FREE(dy1); GC_FREE(dw);
  return ival;
}

