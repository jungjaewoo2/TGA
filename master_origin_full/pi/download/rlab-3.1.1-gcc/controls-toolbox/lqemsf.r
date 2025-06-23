//----------------------------------------------------------------------------
//
// lqe
//
// Syntax: G=lqemsf(A,D,M,Q,R,CT,Reps,iterm)
//
// This routine computes the Linear Quadratic Estimator Design for continuous
// systems. The system is given by the following differential equations:
//      .
//      x = Ax + Bu + Dw
//      z = Mx + Gu + v
//
// where z is the vector of sensor measurements, w is the vector if process
// noise, and v is the vector of sensor noise. The system has the following
// process noise and measurement noise covariances:
//
//           E(w) = E(v) = 0
//           E(ww') = Q
//           E(vv') = R
//           E(wv') = 0
//
// where E is the expectation operator. The gain matrix L is designed such
// that the stationary Kalman filter:
//     .
//     x = Ax + Bu + L(z - Mx - Gu)
//
// produces an LQG optimal estimate of x. If CT exists then the process and
// sensor noise are correlated: E(wv') = CT.
//
// This routine takes advantage of the dual nature of the optimal estimator,
// by using the lqrmsf routine to compute the optimal estimator.
//
// The routine lqrmsf uses the matrix sign function to solve the Riccati
// Equation to compute the estimator gains. The input Reps, is the tolerance
// for convergence, and iterm is the maximum number of iterations allowed
// in lqrmsf.
//
// Note: Three matrices are returned in a list.
//
//       G.l = Optimal Estimator Gain matrix
//       G.p = Riccatti Eqn. Solution, which is the estimate error covariance
//
// Ref: (1) Skelton, R. "Dynamic Systems Control Linear System Analysis and
//          Synthesis," John Wiley and Sons, 1988.
//
// Copyright (C), by Jeffrey B. Layton, 1994
// Version JBL 940918
//----------------------------------------------------------------------------

require lqrmsf

lqemsf = function(a,d,m,q,r,ct,Reps,iterm)
{

   if ( nargs < 5) {
       error("LQEMSF: Wrong number of input arguments.");
   }
 
// Set defaults if variables are not input.
// ========================================

   if (!exist(iterm)) {
       itermax=1000;
   else
       itermax=iterm;
   }

   if (!exist(Reps)) {
       eps=1.0e-10;
   else
       eps=Reps
   }

// Don't do any error checking, let lqr take care of that.
   if ( exist(ct)) {
        Dum=lqrmsf(a',m',d*q*d',r,d*ct,eps,itermax);
   else
        Dum=lqrmsf(a',m',d*q*d',r,eps,itermax);
   }
    
   return << l=Dum.G'; p=Dum.P' >>;
};

