//----------------------------------------------------------------------------
//
// lqrymsf
//
// Symsfntax: G=lqrymsf(A,B,C,D,Q,R,Reps,iterm)
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
// It uses the matrix sign function (msf) technique to solve the Riccati
// equation necessary to compute the gains. The input Reps is the tolerance
// for the iteration procedure and the variable iterm is the maximum
// number of iterations.
//
//
// Note: Three matrices are returned in a list.
//
//       G.k = Optimal Feedback Gain matrix
//       G.s = Steady-State Solution to the Algebraic Riccatti Eqn.
//
//  Ref.: (1) Chapter 5 of Junkins & Kim, Dyn. & Ctrl. of Structures.
//        (2) Bierman, G.J.,
//            "Computational Aspects of the Matrix Sign Function to the ARE"
//            Proc. 23rd CDC, 1984, pp.514-519.
//
// Copyright (C), by Jeffrey B. Layton, 1994
// Version JBL 940918
//----------------------------------------------------------------------------

require lqrmsf

lqrymsf = function(a,b,c,d,q,r,Reps,iterm)
{

   if (nargs < 6 ) {
       error("LQRY: Wrong number of input arguments.");
   }

// Set defaults if variables are not input.
// ========================================

   if (!exist(iterm)) {
       itermax=1000;
   } else {
       itermax=iterm;
   }

   if (!exist(Reps)) {
       eps=1.0e-10;
   } else {
       eps=Reps
   }

   t=lqrmsf(a,b, c'*q*c, r+d'*q*d, eps, itermax, c'*q*d);

   return << k=t.G; s=t.P >>;
};

