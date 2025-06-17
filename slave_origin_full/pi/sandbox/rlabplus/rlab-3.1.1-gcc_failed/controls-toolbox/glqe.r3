//--------------------------------------------------------------------------
//
// glqe
//
// Syntax: G=glqe(A,D,M,Q,R,CT)
//         where G.dl = Gradient of the Optimal Estimator Gain matrix
//               G.dp = Gradient of the Riccatti Eqn. SOlution.
//
// This routine comptues the gradient of the Linear Quadratic Estimator with
// respect to a scalar parameters 'p' for continuous systems. The system is
// given by the following differential equations:
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
// where E is the expectation operator. The routine returns the gradient of
// the estimator gain, dl, as well as the gradient of the Riccatti solution.
//
// This routine takes advantage of the dual nature of the optimal estimator,
// by using the lqr routine to compute the optimal estimator.
//
// Note: Two matrices are returned in a list.
//
//       G.dl = Gradient of the Optimal Estimator Gain matrix
//       G.dp = Gradient of the Riccatti Eqn. SOlution.
//
// Ref: (1) Skelton, R. "Dynamic Systems Control Linear System Analysis and
//          Synthesis," John Wiley and Sons, 1988.
//
// Copyright (C), by Jeffrey B. Layton, 1994
// Version JBL 940918
//--------------------------------------------------------------------------

require glqr

glqe = function(a,da,d,dd,m,dm,q,dq,r,dr,ct,dct)
{

   if ( (nargs != 10) || (nargs != 12) ) {
       error("DLQE: Wrong number of input arguments.");
   }

// Create gradient of input matrices for mxglqr
   dalocal=da';
   dblocal=dm';
   dqlocal=dd*q*d' + d*dq*d' + d*q*dd';
   drlocal=dr;

// Don't do any error checking, let lqr take care of that.
   if ( exist(ct)) {
        dctlocal=dd*ct + d*dct';
        Dum=glqr(a',dalocal,m',dblocal,d*q*d',dqlocal,r,drlocal,d*ct,dctlocal);
   } else {
        Dum=glqr(a',dalocal,m',dblocal,d*q*d',dqlocal,r,drlocal);
   }
    
   return << l=d.dk'; p=d.dp' >>
};

