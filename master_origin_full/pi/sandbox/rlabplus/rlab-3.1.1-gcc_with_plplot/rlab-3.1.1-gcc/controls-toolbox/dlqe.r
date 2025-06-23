//--------------------------------------------------------------------------
//
// dlqe
//
// Syntax: R = dlqe(a,g,c,q,r,t)
//         where R.l = gain matrix
//               R.m = the Riccati equation solution
//               R.p = the estimate error covariance
//               R.e = the closed-loop eigenvalues of the estimator
//
//	Discrete linear quadratic estimator design for the system:
//		x[n+1] = Ax[n] + Bu[n] + Gw[n]	  {State equation}
//		z[n]   = Cx[n] + Du[n] +  v[n]	  {Measurements}
//	with process noise and measurement noise covariances:
//		E{w} = E{v} = 0,  E{ww'} = Q,  E{vv'} = R,  E{wv'} = 0
//	R.l is the gain matrix L such that the discrete, stationary Kalman 
//	filter with time and observation update equations:
//	   _         *               *      _                _
//	   x[n+1] = Ax[n] + Bu[n]    x[n] = x[n] + L(z[n] - Cx[n] - Du[n])
//	produces an LQG optimal estimate of x.  The estimator can be 
//	formed using DESTIM.
//
//	It returns the gain matrix L(=R.l), the Riccati
//	equation solution M(=R.m), the estimate error covariance after the 
//	measurement update:            *    *
//                            R.p = E{[x-x][x-x]'}
//	and the closed-loop eigenvalues of the estimator, R.e=eig(A-A*L*C).val.
//
//	R = dlqe(A,G,C,Q,R,N) solves the discrete est. problem
//	when the process and sensor noise is correlated: E{wv'} = N.
//
//	See also: dlqew, lqed, destim.
//--------------------------------------------------------------------------

require dlqr

dlqe = function(a,g,c,q,r,t)
{
  global(eps,pi)
  
  if (nargs < 5 || nargs > 6) {
     error("Wrong number of arguments.");
  }
  // Use DLQR to calculate estimator gains using duality
  if (nargs==5) {
    </e;k;s/> = dlqr(a',c',g*q*g',r);
    m = s';
    l = (m*c')/(c*m*c' + r);
  else
    </e;k;s/> = dlqr(a',c',g*q*g',r,g*t);
    m = s';
    // Because of the notation used for discrete kalman filters, the kalman
    // gain matrix shows up in the observation update equation.  When designing
    // with cross terms, this requires that the A matrix be inverted when 
    // computing the L matrix.
    if (rcond(a)<eps) {
       disp("Warning: The A matrix must be non-singular for Kalman gain calculation.");
    } 
    l = (m*c'+a\g*t)/(c*m*c' + r);
    a = a-g*t/r*c;
    q = q-t/r*t';
  }
  p = a\(m-g*q*g')/a';
  
  return <<l=l;m=m;p=p;e=e>>;
  
};

