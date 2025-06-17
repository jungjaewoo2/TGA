//----------------------------------------------------------------------------
//
// lqry
//
// Syntax: G = lqry(A,B,Q,R,CT)
//         where G.e = Closed Loop eigenvalues
//               G.k = Optimal Feedback Gain matrix
//               G.s = Steady-State Solution to the Algebraic Riccatti Eqn.
//
// This routine computes the Linear Quadratic Design for Continuous Systems
// with weighting on the outputs. It calculates the optimal feedback gain
// matrix K such that the feedbakc law u = -Kx minimizes the following
// cost function:
//
//           Inf
//     J = Integral ( y'Qy + u'Ru ) dt
//            0
//
// subject to the constraint
//
//   .
//   x = Ax + Bu
//   y = Cx + Du
//
// Note: Three matrices are returned in a list.
//
//       G.k = Optimal Feedback Gain matrix
//       G.s = Steady-State Solution to the Algebraic Riccatti Eqn.
//       G.e = Closed Loop eigenvalues
//
//
// Copyright (C), by Jeffrey B. Layton, 1994
// Version JBL 940919
//----------------------------------------------------------------------------

require lqr

lqry = function(a,b,c,d,q,r)
{

  if (nargs != 6 ) {
     error("lqry: Wrong number of input arguments.");
  }
  qq = c'*q*c;
  rr = r + d'*q*d;
  nn = c'*q*d;
  return lqr(a,b,qq,rr,nn);
};

