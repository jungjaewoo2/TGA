#include "S.h"
#include "loess.h"

#define min(x,y)  ((x) < (y) ? (x) : (y))
#define max(x,y)  ((x) > (y) ? (x) : (y))

double pf(double q, double df1, double df2)
{
  double ibeta();
  return (ibeta(q*df1/(df2+q*df1), df1/2, df2/2));
}

double qt(double p, double df)
{
  double t, invibeta();
  t = invibeta(fabs(2*p-1), 0.5, df/2);
  return((p>0.5?1:-1) * sqrt(t*df/(1-t)));
}

/**********************************************************************/
 /*
 * Incomplete beta function.
 * Reference:  Abramowitz and Stegun, 26.5.8.
 * Assumptions: 0 <= x <= 1; a,b > 0.
 */
#define DOUBLE_EPS      2.2204460492503131E-16
#define IBETA_LARGE     1.0e30
#define IBETA_SMALL     1.0e-30

double ibeta(double x, double a, double b)
{
  int flipped = 0, i, k, count;
  double I, temp, pn[6], ak, bk, next, prev, factor, val;

  if (x <= 0)
    return(0);
  if (x >= 1)
    return(1);

  /* use ibeta(x,a,b) = 1-ibeta(1-x,b,a) */
  if ((a+b+1)*x > (a+1)) {
    flipped = 1;
    temp = a;
    a = b;
    b = temp;
    x = 1 - x;
  }

  pn[0] = 0.0;
  pn[2] = pn[3] = pn[1] = 1.0;
  count = 1;
  val = x/(1.0-x);
  bk = 1.0;
  next = 1.0;
  do {
    count++;
    k = count/2;
    prev = next;
    if (count%2 == 0)
      ak = -((a+k-1.0)*(b-k)*val)/
          ((a+2.0*k-2.0)*(a+2.0*k-1.0));
    else
      ak = ((a+b+k-1.0)*k*val)/
          ((a+2.0*k)*(a+2.0*k-1.0));
    pn[4] = bk*pn[2] + ak*pn[0];
    pn[5] = bk*pn[3] + ak*pn[1];
    next = pn[4] / pn[5];
    for (i=0; i<=3; i++)
      pn[i] = pn[i+2];
    if (fabs(pn[4]) >= IBETA_LARGE)
      for (i=0; i<=3; i++)
        pn[i] /= IBETA_LARGE;
    if (fabs(pn[4]) <= IBETA_SMALL)
      for (i=0; i<=3; i++)
        pn[i] /= IBETA_SMALL;
  } while (fabs(next-prev) > DOUBLE_EPS*prev);
  factor = a*log(x) + (b-1)*log(1-x);
  factor -= gamma(a+1) + gamma(b) - gamma(a+b);
  I = exp(factor) * next;
  return(flipped ? 1-I : I);
}

/*
 * Rational approximation to inverse Gaussian distribution.
 * Absolute error is bounded by 4.5e-4.
 * Reference: Abramowitz and Stegun, page 933.
 * Assumption: 0 < p < 1.
 */

static double num[] = {
  2.515517,
  0.802853,
  0.010328
};

static double den[] = {
  1.000000,
  1.432788,
  0.189269,
  0.001308
};

double invigauss_quick(double p)
{
  int lower;
  double t, n, d, q;

  if(p == 0.5)
    return(0);
  lower = p < 0.5;
  p = lower ? p : 1 - p;
  t = sqrt(-2 * log(p));
  n = (num[2]*t + num[1])*t + num[0];
  d = ((den[3]*t + den[2])*t + den[1])*t + den[0];
  q = lower ? n/d - t : t - n/d;
  return(q);
}

/*
 * Inverse incomplete beta function.
 * Assumption: 0 <= p <= 1, a,b > 0.
 */

double invibeta(double p, double a, double b)
{
  int i;
  double ql, qr, qm, qdiff;
  double pl, pr, pm, pdiff;
  double invibeta_quick(), ibeta();

  /*        MEANINGFUL(qm);*/
  qm = 0;
  if(p == 0 || p == 1)
    return(p);

  /* initialize [ql,qr] containing the root */
  ql = qr = invibeta_quick(p, a, b);
  pl = pr = ibeta(ql, a, b);
  if(pl == p)
    return(ql);
  if(pl < p)
    while(1) {
    qr += 0.05;
    if(qr >= 1) {
      pr = qr = 1;
      break;
    }
    pr = ibeta(qr, a, b);
    if(pr == p)
      return(pr);
    if(pr > p)
      break;
    }
    else
      while(1) {
      ql -= 0.05;
      if(ql <= 0) {
        pl = ql = 0;
        break;
      }
      pl = ibeta(ql, a, b);
      if(pl == p)
        return(pl);
      if(pl < p)
        break;
      }

      /* a few steps of bisection */
      for(i = 0; i < 5; i++) {
        qm = (ql + qr) / 2;
        pm = ibeta(qm, a, b);
        qdiff = qr - ql;
        pdiff = pm - p;
        if(fabs(qdiff) < DOUBLE_EPS*qm || fabs(pdiff) < DOUBLE_EPS)
          return(qm);
        if(pdiff < 0) {
          ql = qm;
          pl = pm;
        } else {
          qr = qm;
          pr = pm;
        }
      }

      /* a few steps of secant */
      for(i = 0; i < 40; i++) {
        qm = ql + (p-pl)*(qr-ql)/(pr-pl);
        pm = ibeta(qm, a, b);
        qdiff = qr - ql;
        pdiff = pm - p;
        if(fabs(qdiff) < 2*DOUBLE_EPS*qm || fabs(pdiff) < 2*DOUBLE_EPS)
          return(qm);
        if(pdiff < 0) {
          ql = qm;
          pl = pm;
        } else {
          qr = qm;
          pr = pm;
        }
      }

      /* no convergence */
      return(qm);
}

/*
 * Quick approximation to inverse incomplete beta function,
 * by matching first two moments with the Gaussian distribution.
 * Assumption: 0 < p < 1, a,b > 0.
 */

double invibeta_quick(double p, double a, double b)
{
  //double x, m, s, fmax(), fmin(), invigauss_quick();
  double x, m, s, rval, invigauss_quick();

  x = a + b;
  m = a / x;
  s = sqrt((a*b) / (x*x*(x+1)));
  rval = min(1.0, invigauss_quick(p)*s + m);
  rval = max(0.0, rval);

  //return(fmax(0.0, fmin(1.0, invigauss_quick(p)*s + m)));
  return rval;
}

// static double
//     fmin(a, b)
//     double a, b;
// {
//   return(a < b ? a : b);
// }

// static double
//     fmax(a, b)
//     double a, b;
// {
//   return(a > b ? a : b);
// }

// int
//     max(a, b)
//     int a, b;
// {
//   return(a > b ? a : b);
// }

typedef double doublereal;
typedef int integer;

void Recover(a, b)
    char    *a;
int     *b;
{
  printf(a);
  //exit(1);
}

void
    Warning(a, b)
    char    *a;
int     *b;
{
  printf(a);
}
