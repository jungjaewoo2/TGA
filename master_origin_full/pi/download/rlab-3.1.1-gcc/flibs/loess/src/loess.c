#include "S.h"
#include "loess.h"
#include "gc.h"

static char surf_stat[64];

int
loess_(y, x_, size_info, weights, span, degree, parametric, drop_square,
       normalize, statistics, surface, cell, trace_hat_in, iterations,
       fitted_values, fitted_residuals, enp, s, one_delta, two_delta,
       pseudovalues, trace_hat_out, diagonal, robust, divisor,
       parameter, a, xi, vert, vval)
    double  *y, *x_, *weights, *span, *cell, *pseudovalues,
  *fitted_values, *fitted_residuals, *enp, *s, *one_delta, *two_delta,
  *trace_hat_out, *diagonal, *robust, *divisor, *xi, *vert, *vval;
  int  *size_info, *degree, *parametric, *drop_square, *normalize,
  *iterations, *parameter, *a;
  char  *statistics, **surface, **trace_hat_in;
{
  double  *x=0, *x_tmp=0, new_cell, trL, delta1=0, delta2=0, sum_squares = 0,
  *pseudo_resid=0, *temp=0, *xi_tmp=0, *vert_tmp=0, *vval_tmp=0,
  *diag_tmp=0, trL_tmp = 0, d1_tmp = 0, d2_tmp = 0, sum, mean;
  int  i, j, k, p, N, D, sum_drop_sqr = 0, sum_parametric = 0,
  setLf,  nonparametric = 0, *order_parametric=0,
  *order_drop_sqr=0, zero = 0, max_kd, *a_tmp=0, *param_tmp=0;
  int  cut, comp();
  char new_stat[64];
  void condition();

  int jj=0;

  D = size_info[0]; // p
  N = size_info[1]; // n
  max_kd = (N > 200 ? N : 200);
  *one_delta = *two_delta = *trace_hat_out = 0;

  x = (double *) GC_malloc(D * N * sizeof(double));
  x_tmp = (double *) GC_malloc(D * N * sizeof(double));
  temp = (double *) GC_malloc(N * sizeof(double));
  a_tmp = (int *) GC_malloc(max_kd * sizeof(long));
  xi_tmp = (double *) GC_malloc(max_kd * sizeof(double));
  vert_tmp = (double *) GC_malloc(D * 2 * sizeof(double));
  vval_tmp = (double *) GC_malloc((D + 1) * max_kd * sizeof(double));
  diag_tmp = (double *) GC_malloc(N * sizeof(double));
  param_tmp = (int *) GC_malloc(N * sizeof(long));
  order_parametric = (int *) GC_malloc(D * sizeof(long));
  order_drop_sqr = (int *) GC_malloc(D * sizeof(long));
  if((*iterations) > 0)
    pseudo_resid = (double *) GC_malloc(N * sizeof(double));

  new_cell = (*span) * (*cell);
  for(i = 0; i < N; i++)
    robust[i] = 1;
  for(i = 0; i < (N * D); i++)
    x_tmp[i] = x_[i];
  if((*normalize) && (D > 1)) {
    cut = ceil(0.100000000000000000001 * N);
    for(i = 0; i < D; i++) {
      k = i * N;
      for(j = 0; j < N; j++)
        temp[j] = x_[k + j];
      qsort(temp, N, sizeof(double), comp);
      sum = 0;
      for(j = cut; j <= (N - cut - 1); j++)
        sum = sum + temp[j];
      mean = sum / (N - 2 * cut);
      sum = 0;
      for(j = cut; j <= (N - cut - 1); j++) {
        temp[j] = temp[j] - mean;
        sum = sum + temp[j] * temp[j];
      }
      divisor[i] = sqrt(sum / (N - 2 * cut - 1));
      for(j = 0; j < N; j++) {
        p = k + j;
        x_tmp[p] = x_[p] / divisor[i];
      }
    }
  }
  else
    for(i = 0; i < D; i++) divisor[i] = 1;
  j = D - 1;
  for(i = 0; i < D; i++) {
    sum_drop_sqr = sum_drop_sqr + drop_square[i];
    sum_parametric = sum_parametric + parametric[i];
    if(parametric[i])
      order_parametric[j--] = i;
    else
      order_parametric[nonparametric++] = i;
  }

  for(i = 0; i < D; i++) {
    order_drop_sqr[i] = 2 - drop_square[order_parametric[i]];
    k = i * N;
    p = order_parametric[i] * N;
    for(j = 0; j < N; j++)
      x[k + j] = x_tmp[p + j];
  }

  if((*degree) == 1 && sum_drop_sqr) {
    fprintf(stderr, "Specified the square of a factor predictor to be dropped when degree = 1");
    return 1;
  }

  if(D == 1 && sum_drop_sqr) {
    fprintf(stderr, "Specified the square of a predictor to be dropped with only one numeric predictor");
    return 1;
  }

  if(sum_parametric == D) {
    fprintf(stderr, "Specified parametric for all predictors");
    return 1;
  }

  for(j = 0; j <= (*iterations); j++)
  {
    if (j)
      sprintf(new_stat,"none");
    else
      sprintf(new_stat,statistics);

    for(i = 0; i < N; i++)
      robust[i] = weights[i] * robust[i];
    condition(surface, new_stat, trace_hat_in);
    setLf = !strcmp(surf_stat, "interpolate/exact");

    loess_raw(y, x, weights, robust, &D, &N, span, degree,
              &nonparametric, order_drop_sqr, &sum_drop_sqr,
              &new_cell, &surf_stat, fitted_values, parameter, a,
              xi, vert, vval, diagonal, &trL, &delta1, &delta2,
              &setLf);

    if(j == 0)
    {
      *trace_hat_out = trL;
      *one_delta = delta1;
      *two_delta = delta2;
    }
    for(i = 0; i < N; i++)
      fitted_residuals[i] = y[i] - fitted_values[i];
    if(j < (*iterations))
    {
      F77_SUB(lowesw)(fitted_residuals, &N, robust, temp);
    }
  }

  if((*iterations) > 0)
  {
    F77_SUB(lowesp)(&N, y, fitted_values, weights, robust, temp, pseudovalues);
    loess_raw(pseudovalues, x, weights, weights, &D, &N, span,
              degree, &nonparametric, order_drop_sqr, &sum_drop_sqr,
              &new_cell, &surf_stat, temp, param_tmp, a_tmp, xi_tmp,
              vert_tmp, vval_tmp, diag_tmp, &trL_tmp, &d1_tmp, &d2_tmp, &zero);
    for(i = 0; i < N; i++)
      pseudo_resid[i] = pseudovalues[i] - temp[i];
  }

  if((*iterations) == 0)
    for(i = 0; i < N; i++)
      sum_squares = sum_squares + weights[i] *
          fitted_residuals[i] * fitted_residuals[i];
  else
    for(i = 0; i < N; i++)
      sum_squares = sum_squares + weights[i] *
          pseudo_resid[i] * pseudo_resid[i];
  *enp = (*one_delta) + 2 * (*trace_hat_out) - N;
  *s = sqrt(sum_squares / (*one_delta));

  GC_free(x);
  GC_free(x_tmp);
  GC_free(temp);
  GC_free(xi_tmp);
  GC_free(vert_tmp);
  GC_free(vval_tmp);
  GC_free(diag_tmp);
  GC_free(a_tmp);
  GC_free(param_tmp);
  GC_free(order_parametric);
  GC_free(order_drop_sqr);
  if((*iterations) > 0)
    GC_free(pseudo_resid);

  return 0;
}

void
loess_free_mem(lo)
struct	loess_struct	*lo;
{
  GC_free(lo->in.x);
  GC_free(lo->in.y);
  GC_free(lo->in.weights);
  GC_free(lo->out.fitted_values);
  GC_free(lo->out.fitted_residuals);
  GC_free(lo->out.pseudovalues);
  GC_free(lo->out.diagonal);
  GC_free(lo->out.robust);
  GC_free(lo->out.divisor);
  GC_free(lo->kd_tree.parameter);
  GC_free(lo->kd_tree.a);
  GC_free(lo->kd_tree.xi);
  GC_free(lo->kd_tree.vert);
  GC_free(lo->kd_tree.vval);
}

void
condition (char *surface, char * new_stat, char *trace_hat_in)
{
  if(!strcmp(surface, "interpolate"))
  {
    if(!strcmp(new_stat, "none"))
      sprintf(surf_stat,"interpolate/none");
    else if(!strcmp(new_stat, "exact"))
      sprintf(surf_stat,"interpolate/exact");
    else if(!strcmp(new_stat, "approximate"))
    {
      if(!strcmp(trace_hat_in, "approximate"))
        sprintf(surf_stat,"interpolate/2.approx");
      else if(!strcmp(trace_hat_in, "exact"))
        sprintf(surf_stat,"interpolate/1.approx");
    }
  }
  else if(!strcmp(surface, "direct"))
  {
    if(!strcmp(new_stat, "none"))
      sprintf(surf_stat,"direct/none");
    else if(!strcmp(new_stat, "exact"))
      sprintf(surf_stat,"direct/exact");
    else if(!strcmp(new_stat, "approximate"))
      sprintf(surf_stat,"direct/approximate");
  }
}

int comp(double *d1, double *d2)
{
  if(*d1 < *d2)
    return(-1);
  else if(*d1 == *d2)
    return(0);
  else
    return(1);
}
