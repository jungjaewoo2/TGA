//----------------------------------------------------------------------
//
// dlyap
//
// syntax: x = dlyap(a,c)
//
// Discrete Lyapunov equation solution.
// X = dlyap(A,C) solves the discrete Lyapunov equation:
//
//	A*X*A' + C = X
//
// See also lyap.
//
//----------------------------------------------------------------------
require lyap

dlyap=function(a,c)
{
  m = a.nr;
  n = a.nc;
  a = (a+eye(m,m))\(a-eye(m,m));
  c = (eye(m,m)-a)*c*(eye(m,m)-a')/2;
  x = lyap(a,,c);
  return x;
};
