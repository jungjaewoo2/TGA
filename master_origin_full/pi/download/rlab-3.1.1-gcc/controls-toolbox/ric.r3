//--------------------------------------------------------------------------
//
// ric
// 
// syntax: </kerr;serr/> = ric(a,b,q,r,k,s,n)
//
//	Riccati residual calculation.
//
//	</Kerr;Serr/> = ric(A,B,Q,R,K,S) Computes the error in the solution
//	to the Riccati equation.  Kerr is the error in the gain matrix
//	and Serr is the residual error in the Riccati equation.
//		          -1                           -1
//		Kerr = K-R  B'S;  Serr = SA + A'S - SBR  B'S + Q
//
//	</Kerr;Serr/> = ric(A,B,Q,R,K,S,N) computes the error in the 
//	solution to the Riccati equation with cross weighting term.
//	          -1                                    -1
//	Kerr = K-R  (N'+B'S);  Serr = SA + A'S - (SB+N)R  (N'+B'S) + Q
//
//	See also: dric, are, lqe, lqr
//
//--------------------------------------------------------------------------
ric = function(a,b,q,r,k,s,n)
{
  if (nargs == 7) {
     serr = s*a + a'*s - (s*b+n)/r*(n'+b'*s) + q;
     kerr = k - r\(n'+b'*s);
  } else {
     serr = s*a+a'*s-s*b/r*b'*s+q;
     kerr = k - r\b'*s;
  }
  return <<kerr=kerr;serr=serr>>;
};

