//--------------------------------------------------------------------------
//
// lqe
//
// syntax: </e;l;p/> = lqe(a,g,c,q,r,t)
//
//	Linear quadratic estimator design. For the continuous-time system:
//		.
//		x = Ax + Bu + Gw            {State equation}
//		z = Cx + Du + v             {Measurements}
//	with process noise and measurement noise covariances:
//		E{w} = E{v} = 0,  E{ww'} = Q,  E{vv'} = R, E{wv'} = 0
//
//	</E;L;P/> = lqe(A,G,C,Q,R) returns the gain matrix L such that the 
//	stationary Kalman filter:
//	      .
//	      x = Ax + Bu + L(z - Cx - Du)
//
//	produces an LQG optimal estimate of x. The estimator can be formed
//	with estim().
//
//	</E;L;P/> = lqe(A,G,C,Q,R) returns the gain matrix L, the Riccati
//	equation solution P which is the estimate error covariance, and 
//	the closed loop eigenvalues of the estimator: E = eig(A-L*C).val[:]
//
//	</E;L;P/> = lqe(A,G,C,Q,R,N) solves the estimator problem when the
//	process and sensor noise is correlated: E{wv'} = N.
//
//	See also: lqew, lqe2, and estim
//
//--------------------------------------------------------------------------
require lqr nargchk

lqe = function(a,g,c,q,r,t)
{
  msg = nargchk(5,6,nargs);
  if (msg!="") { error(msg); }

  // Calculate estimator gains using lqr() and duality:
  if (nargs==5) {
    </e;k;s/> = lqr(a',c',g*q*g',r);
  else
    </e;k;s/> = lqr(a',c',g*q*g',r,g*t);
  } 
  l=k';
  p=s';
  return <<e=e;l=l;p=p>>;
};

