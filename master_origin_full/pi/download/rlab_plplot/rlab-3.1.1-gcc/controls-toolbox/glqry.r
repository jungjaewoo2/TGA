//--------------------------------------------------------------------------
//
// glqry
//
// Syntax: G=glqry(A,DA,B,DA,C,DC,D,DD,Q,DQ,R,DR)
//         where G.dk = Optimal Feedback Gain matrix
//               G.ds = Steady-State Solution to the Algebraic Riccatti Eqn.
//
// This routine computes the gradient of the Linear Quadratic Design with
// respect to a scalar parameter 'p' for continuous systems with weightings
// on the outputs. It compute the gradient of the feedback matrix K for the
// feedback law,
//
//   u = -Kx
//
// with respect to 'p'. The system is described by,
//   .
//   x = Ax + Bu
//   y = Cx + Du
//
// Note: Two matrices are returned in a list.
//
//       G.dk = Optimal Feedback Gain matrix
//       G.ds = Steady-State Solution to the Algebraic Riccatti Eqn.
//
// Copyright (C), by Jeffrey B. Layton, 1994
// Version JBL 940918
//--------------------------------------------------------------------------

require glqr

glqry = function(a,da,b,db,c,dc,d,dd,q,dq,r,dr)
{
  if (nargs != 12) {
       error("DLQRY: Wrong number of input arguments.");
   }

// Begin: use glqr to solve problem
   dqlocal=dc'*q*c + c'*dq*c + c'*q*dc;
   drlocal=dr+dd'*q*d + d'*dq*d + d'*q*dd;
   dctlocal=dc'*q*d + c'*dq*d + c'*q*dd;

   r=glqr(a,da,b,db,c'*q*c,dqlocal,r+d'*q*d,drlocal,c'*q*d,dctlocal);

   return << k=r.dk; s=r.dp >>;
};

