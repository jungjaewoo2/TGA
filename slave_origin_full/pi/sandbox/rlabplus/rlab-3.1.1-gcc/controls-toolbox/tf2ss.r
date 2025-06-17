//----------------------------------------------------------------------
//
// tf2ss
//
// syntax: </a;b;c;d/> = tf2ss(num, den)
//
// Transfer function to state-space conversion.
// </A;B;C;D/> = tf2ss(NUM,DEN)  calculates the state-space 
// representation:
//	.
//	x = Ax + Bu
//	y = Cx + Du
//
// of the system:
//	        NUM(s) 
//	H(s) = --------
//	        DEN(s)
//
// from a single input.  Vector DEN must contain the coefficients of
// the denominator in descending powers of s.  Matrix NUM must 
// contain the numerator coefficients with as many rows as there are
// outputs y.  The A,B,C,D matrices are returned in controller 
// canonical form.  This calculation also works for discrete systems.
// To avoid confusion when using this function with discrete systems,
// always use a numerator polynomial that has been padded with zeros
// to make it the same length as the denominator.
//
// See also: ss2tf
//
//----------------------------------------------------------------------
tf2ss = function(num, den, typ)
{
  mnum = num.nr;
  nnum = num.nc;
  mden = den.nr;
  n    = den.nc;
  // Check for null systems
  if (n==0 && nnum==0) { 
     return<< a=[]; b=[]; c=[]; d=[]>>; 
  }

  // Strip leading zeros from denominator
  inz  = find(den != 0);
  den  = den[inz[1]:n];
  mden = den.nr;
  n    = den.nc;
  // Check for proper numerator
  if (nnum > n) {
     // Try to strip leading zeros to make proper
     if (all(all(num[;1:(nnum-n)] == 0))) {
        num = num[;(nnum-n+1):nnum];
        mnum = num.nr;
        nnum = num.nc;
     else
        error("Denominator must be higher or equal order than numerator.");
     }   
  }

  // Pad numerator with leading zeros, to make it have the same number of
  // columns as the denominator, and normalize it to den(1)
  num = [zeros(mnum,n-nnum), num]./den[1];

  // Do the D-matrix first
  if (length(num)) {
     d = num[;1];
  else
     d = [];
  }

  // Handle special constant case:
  if (n == 1) {
     return <<a = []; b = []; c = []; d = d>>;
  }

  // Now do the rest, starting by normalizing den to den(1),
  den = den[2:n] ./ den[1];
  a = [-den; eye(n-2,n-1)];
  b = eye(n-1,1);
  if (mnum > 0) {
     c = num[;2:n] - num[;1] * den;
  else
     c = [];
  }
  return <<a=a;b=b;c=c;d=d>>;
};
