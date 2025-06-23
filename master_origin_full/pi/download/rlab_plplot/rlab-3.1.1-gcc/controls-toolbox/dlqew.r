//----------------------------------------------------------------------
//
// dlqew
//
// syntax: </e;l;m;p/> = dlqew(a,g,c,j,q,r)
//
// Discrete linear quadratic estimator design for the system:
//	x[n+1] = Ax[n] + Bu[n] + Gw[n]	          {State equation}
//	z[n]   = Cx[n] + Du[n] + Jw[n] + v[n]     {Measurements}
// with process noise and measurement noise covariances:
//	E{w} = E{v} = 0,  E{ww'} = Q,  E{vv'} = R,  E{wv'} = 0
//
// returns the gain matrix l such that the discrete, stationary 
// Kalman filter with time and observation update equations:
//   _         *               *      _                _
//   x[n+1] = Ax[n] + Bu[n]    x[n] = x[n] + L(z[n] - Cx[n] - Du[n])
// produces an lqg optimal estimate of x.  The estimator can be
// formed using destim. Also returns the Riccati equation solution m, 
// the estimate error covariance after the measurement update: 
//                     *    *
//            P = E{[x-x][x-x]'}
// and the closed-loop eigenvalues of the estimator, e=eig(A-A*L*C).val
//
// See also: dlqe, lqed and destim.
//
//----------------------------------------------------------------------

require dlqe

dlqew = function(a,g,c,j,q,r)
{
  if (nargs != 6) {
     error("Need 6 arguments.");
  }

  rr = r + j*q*j';
  nn = q*j';
  return dlqe(a,g,c,q,rr,nn);
};


