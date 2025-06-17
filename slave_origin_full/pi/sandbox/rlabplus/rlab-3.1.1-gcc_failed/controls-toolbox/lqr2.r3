//----------------------------------------------------------------------
//
// lqr2
//
// syntax: </k;p/> = lqr2(a,b,q,r,s)
//
// Linear-quadratic regulator design for continuous-time systems.
// </K;P/> = lqr2(A,B,Q,R)  calculates the optimal feedback gain matrix
// K such that the feedback law  u = -Kx  minimizes the cost function
//
//      J = Integral {x'Qx + u'Ru} dt
//	                                      .
// subject to the constraint equation:   x = Ax + Bu 
//       
// Also returned is P, the steady-state solution to the associated 
// algebraic Riccati equation:
//                          -1
//        0 = PA + A'P - PBR  B'P + Q
//
// </K;P/> = lqr2(A,B,Q,R,N) includes the cross-term 2x'Nu that 
// relates u to x in the cost functional. 
//
// The controller can be formed with reg.
//
// lqr2 uses the schur algorithm of [1] and is more numerically
// reliable than lqr, which uses eigenvector decomposition.
//
// See also: are, lqr, and lqe2.
//
// References:
//  A.J. Laub, "A Schur Method for Solving Algebraic Riccati
//  Equations", IEEE Transactions on Automatic Control, vol. AC-24,
//  1979, pp. 913-921.
//
//----------------------------------------------------------------------
require abcdchk are

lqr2 = function(a,b,q,r,s)
{
  // Convert data for linear-quadratic regulator problem to data for
  // the algebraic Riccati equation. 
  //	F = A - B*inv(R)*S'
  //	G = B*inv(R)*B'
  //	H = Q - S*inv(R)*S'
  // R must be symmetric positive definite.

  if (nargs <4 || nargs >5) {
     error("Wrong number of input arguments.");
  }
  msg = abcdchk(a,b);
  if (msg != "") { error(msg); }

  n  = b.nr; 
  m  = b.nc;
  ri = inv(r);
  rb = ri*b';
  g  = b*rb;
  if (nargs > 4) {
     rs = ri*s';
     f = a - b*rs;
     h = q - s*rs;
  } else {
     f = a;
     h = q;
     rs = zeros(m,n);
  }

  // Solve ARE:
  p = are(f,g,h);

  // Find gains:
  k = rs + rb*p;
 
  return <<k=k;p=p>>;
};

