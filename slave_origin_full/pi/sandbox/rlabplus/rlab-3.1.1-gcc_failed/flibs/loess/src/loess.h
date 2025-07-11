/* for the meaning of these fields, see struct.m */
/* longs are used here so that the codes can be called from S */

#define TRUE  1
#define FALSE 0

//extern
struct loess_struct {
    struct {
      int    n;
      int    p;
      double  *y;
      double  *x;
      double  *weights;
    } in;
    struct {
      double  span;
      int    degree;
      int    normalize;
      int    parametric[8];
      int    drop_square[8];
      char    *family;
    } model;
    struct {
      char    *surface;
      char    *statistics;
      double  cell;
      char    *trace_hat;
      int    iterations;
    } control;
    struct {
      int *parameter;
      int *a;
      double  *xi;
      double  *vert;
      double  *vval;
    } kd_tree;
    struct {
      double  *fitted_values;
      double  *fitted_residuals;
      double  enp;
      double  s;
      double  one_delta;
      double  two_delta;
      double  *pseudovalues;
      double  trace_hat;
      double  *diagonal;
      double  *robust;
      double  *divisor;
    } out;
};

//extern
struct pred_struct {
  double  *fit;
  double  *se_fit;
  double  residual_scale;
  double  df;
};

//extern
struct anova_struct {
  double  dfn;
  double  dfd;
  double  F_value;
  double  Pr_F;
};

//extern
struct ci_struct {
  double  *fit;
  double  *upper;
  double  *lower;
};

