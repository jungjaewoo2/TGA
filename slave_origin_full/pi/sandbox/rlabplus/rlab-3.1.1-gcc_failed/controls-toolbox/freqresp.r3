//----------------------------------------------------------------------
// freqresp
//
// Low level frequency response function.
//
//	G=freqresp(A,B,C,D,IU,S)
//	G=freqresp(NUM,DEN,S)
//----------------------------------------------------------------------
require polyval ltifr

freqresp=function(a,b,c,d,iu,s)
{
  if (nargs==3) {
    // It is in transfer function form.  Do directly, using Horner's method
    // of polynomial evaluation at the frequency points, for each row in
    // the numerator.  Then divide by the denominator.
    s = c[:];
    for (i in 1:a.nr) {
      g[;i] = polyval(a[i;],s);
    }
    g = polyval(b,s)*ones(1,a.nr).\g;
  } else {
    // It is in state space form.  Reduce to Hessenberg form then directly
    // evaluate frequency response.
    ny = d.nr;
    nu = d.nc;
    nx = a.nr;
    na = a.nc;
    nw = max(size(s));

    // Balance A
    tmp = balance(a);
    t = tmp.t;
    a = tmp.ab;

    b = t \ b;
    c = c * t;

    // Reduce A to Hesenburg form
    tmp = hess(a);
    p = tmp.p;
    a = tmp.h;

    // Apply similarity transformations from Hessenberg
    // reduction to B and C:
    if (nx>0) {
      b = p' * b[;iu];
      c = c * p;
      d = d[;iu];
      g = ltifr(a,b,s[:]);
      g = (c * g + diag(d) * ones(ny,nw)).';
    } else {
      d = d[;iu];
      g = (diag(d) * ones(ny,nw)).';
    }
  }
  return g;
};


