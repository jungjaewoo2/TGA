//--------------------------------------------------------------------------
//
// are
//
// Syntax: X = are(A,B,C)
//
// Algebraic Riccati Equation solution.
// X = are(A, B, C) returns the stablizing solution (if it
// exists) to the continuous-time Riccati equation:
//
//        A'*X + X*A - X*B*X + C = 0
//
// assuming B is symmetric and nonnegative definite and C is
// symmetric.
//
// see also: ric
//--------------------------------------------------------------------------
require schord

are = function(A,B,C)
{
  global(eps)

  // check for correct input problem
  nr = A.nr;
  nc = A.nc; 
  n  = nr;
  if (nr != nc) {
     error("Nonsquare A matrix"); 
  }
  nr = B.nr;
  nc = B.nc;
  if (nr!=n || nc!=n) {
     error("Incorrectly dimensioned B matrix"); 
  }
  nr = C.nr;
  nc = C.nc;
  if (nr!=n || nc!=n) {
     error("Incorrectly dimensioned C matrix"); 
  }
  tmp = schur([A, -B; -C, -A']*(1.0+eps*eps*sqrt(-1))); 
  q = tmp.z;
  t = tmp.t; 
  tol = 10.0*eps*max(abs(diag(t)));	// ad hoc tolerance
  ns = 0;
  //
  //  Prepare an array called index to send message to ordering routine 
  //  giving location of eigenvalues with respect to the imaginary axis.
  //  -1  denotes open left-half-plane
  //   1  denotes open right-half-plane
  //   0  denotes within tol of imaginary axis
  //  
  index = [];
  for (i in 1:2*n)
  {
    if (real(t[i;i]) < -tol) {
    	index = [index, -1] ;
	ns = ns + 1;
    } else { if (real(t[i;i]) > tol) {
	index = [ index, 1 ];
    } else {
	index = [ index, 0 ];
    }}
  }
  if (ns != n) {
     error("No solution: (A,B) may be uncontrollable or no solution exists"); 
  }
  </ q; t /> = schord(q,t,index);
  X = real(q[n+1:n+n;1:n]/q[1:n;1:n]);
  
  return X;
};


