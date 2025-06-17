//--------------------------------------------------------------------------
//
// dlqry
//
// Syntax: A=dlqry(a,b,c,d,q,r)
//         where A.e = Closed Loop eigenvalues
//               A.k = Optimal Feedback Gain matrix
//               A.s = Steady-State Solution to the Algebraic Riccatti Eqn.
//
// This routine computes the Linear Quadratic Design for Discrete Systems
// with weighting on the outputs. It calculates the optimal feedback gain
// matrix K such that the feedbakc law u(n) = -Kx(n) minimizes the following
// cost function:
//
//     J = Sum ( y'Qy + u'Ru ) dt
//
// subject to the constraint
//
//   x[n+1] = Ax[n] + Bu[n] 
//     y[n] = Cx[n] + Du[n]
//
// Note: Three matrices are returned in a list.
//
//       A.k = Optimal Feedback Gain matrix
//       A.s = Steady-State Solution to the Algebraic Riccatti Eqn.
//       A.e = Closed Loop eigenvalues
//
// Copyright (C), by Jeffrey B. Layton, 1994
// Version JBL 940918
//--------------------------------------------------------------------------

require dlqr

dlqry = function(a,b,c,d,q,r)
{

   if (nargs < 6) {
       error("DLQRY: Wrong number of input arguments.");
   }

   t=dlqr(a,b, c'*q*c, r+d'*q*d, c'*q*d);

   return << k=t.k; s=t.s; e=t.e >>;
};

