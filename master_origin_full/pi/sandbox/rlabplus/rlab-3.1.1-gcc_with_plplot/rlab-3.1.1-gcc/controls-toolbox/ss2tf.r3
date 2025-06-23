//--------------------------------------------------------------------------
//
// ss2tf
//
// syntax: </den;num/> = ss2tf(a,b,c,d,iu)
//
// State-space to transfer function conversion.
// </DEN;NUM/> = ss2tf(A,B,C,D,iu)  calculates the transfer function:
//
//	        NUM(s)          -1
//	H(s) = -------- = C(sI-A) B + D
//	        DEN(s)
// of the system:
//	.
//	x = Ax + Bu
//	y = Cx + Du
//
// from the iu'th input.  Vector DEN contains the coefficients of the
// denominator in descending powers of s.  The numerator coefficients
// are returned in matrix NUM with as many rows as there are 
// outputs y.
//
// See also: tf2ss
//
//--------------------------------------------------------------------------
require abcdchk isempty nargchk poly

ss2tf = function(a,b,c,d,iu)
{
  local (a,b,c,d,iu)
  
  msg = nargchk(4,5,nargs);
  if (msg!="") { error(msg); }
  msg = abcdchk(a,b,c,d);
  if (msg!="") { error(msg); }

  mc = d.nr;
  nu = d.nc;
  if (nargs == 4) {
    if (nu <= 1) {
      iu = 1;
    } else {
      error("iu must be specified for systems with more than one input.");
    }
  }

  den = poly(a);
  if (!isempty(b)) { b = b[;iu]; }
  if (!isempty(d)) { d = d[;iu]; }

  // System is just a gain or it has only a denominator:
  if (isempty(b) && isempty(c))
  {
     num = d;
     if (isempty(d) && isempty(a)) { den = []; }
     return <<den=den;num=num>>;
  }

  nc = length(a);
  num = ones(mc, nc+1);
  for (i in 1:mc) {
	num[i;] = poly(a-b*c[i;]) + (d[i] - 1) * den;
  }
  return <<den=den;num=num>>;
};
