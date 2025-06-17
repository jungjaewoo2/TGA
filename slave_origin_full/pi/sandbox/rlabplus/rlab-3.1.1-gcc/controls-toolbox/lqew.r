//----------------------------------------------------------------------
//
// lqew
//
// syntax: </e;l;p/> = lqew(a,g,c,j,q,r)
//
// Linear quadratic estimator design for the continuous-time 
// system with process noise feedthrough
//	.
//	x = Ax + Bu + Gw	    {State equation}
//	z = Cx + Du + Jw + v	    {Measurements}
//
// and with process noise and measurement noise covariances:
//    E{w} = E{v} = 0,  E{ww'} = Q,  E{vv'} = R,  E{wv'} = 0
//
// </E;L;P/> = lqew(A,G,C,J,Q,R) returns the gain matrix L such that the 
// stationary Kalman filter:  .
//                            x = Ax + Bu + L(z - Cx - Du)
//
// produces an LQG optimal estimate of x. The estimator can be formed
// with estim.  It also returns the Riccati equation solution P which 
// is the estimate error covariance, and the closed loop eigenvalues 
// of the estimator: E = eig(A-L*C).val[:]
//
// See also: lqe, lqe2, and estim.
//
//----------------------------------------------------------------------
require nargchk lqe

lqew=function(a,g,c,j,q,r)
{
  msg = nargchk(6,6,nargs);
  if (msg != "") { error(msg); }

  rr = r + j*q*j';
  nn = q*j';
  return lqe(a,g,c,q,rr,nn);

};

