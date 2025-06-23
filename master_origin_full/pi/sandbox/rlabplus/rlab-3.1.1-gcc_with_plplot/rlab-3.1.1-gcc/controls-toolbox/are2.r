//----------------------------------------------------------------------------
//
// are2
//
// Syntax: X=are2(A,B,C)
//
// This routine solves the Algebraic Riccati Equation. It returns a stabilizing
// solution (if one can be computed) to the continuous time Riccati equation
// given by,
//
//      A'*X + X*A - X*B*X + C = 0
//
//                 -1
//   A'P + PA - PBR  B'P + Q = 0
//
//
// The routine assumes that B is symmetric and nonnegative definite and
// that C is symmetric.
//
// Calling the routine as X=are2(A,B,C) returns the solution to the ARE
// in X.
//
// This routine is taken from the lqr solver so it uses Potter's algorithm
// to the solve the ARE.
//
// Ref: (1) Kailath, "Linear Systems," Prentice-Hall, 1980
//      (2) Laub, A. J., "A Schur Method for Solving Algebraic Riccati
//          Equations," IEEE Trans. AC-24, 1979, pp. 913-921.
//      (3) Junkins, J., and Kim, Y., "Introduction to Dynamics and Control
//          of Flexible Structures," AIAA Inc, Washington D.C., 1993.
//
// Copyright (C), by Jeffrey B. Layton, 1994
// Version JBL 940918
//----------------------------------------------------------------------------

are2 = function(A,B,C)
{
   if (nargs != 3) {
       error("ARE2: Wrong number of input arguments.");
   }

   // Error Traps for input.
   if (A.nr != A.nc) {
       error("ARE2: Nonsquare A matrix");
   }
   if ( (B.nr != A.nr) || (B.nc != A.nr) ) {
       error("ARE2: Incorrectly dimensioned B matrix");
   }
   if ( (C.nr != A.nr) || (C.nc != A.nr) ) {
       error("ARE2: Incorrectly dimensioned C matrix");
   }

   // Start solution by finding the spectral decomposition and eigenvectors
   // of Hamiltonian.

   e=eig( [ A, B; C, -A' ] );
   v=e.vec;
   d=e.val;

   // Sort eigenvectors by sorting eigenvalues (and storing the index).
   index=sort( real( d ) ).idx;
   d=real( d[ index ] );

   if ( !( d[A.nc] < 0 && d[A.nc+1] > 0 ) ) {
       error("ARE2: Can't order eigenvalues");
   }

   // Form the Partitions of the PHI matrix and solve for P
   Phi12=v[ 1:A.nc; index[1:A.nc] ];
   Phi22=v[ (A.nc+1):(2*A.nc); index[1:A.nc] ];
   X=-real(solve(Phi12',Phi22')');

   return X;
};

