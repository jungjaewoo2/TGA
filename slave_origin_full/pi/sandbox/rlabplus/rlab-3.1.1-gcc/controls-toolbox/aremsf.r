//--------------------------------------------------------------------------------
//
// aremsf
//
// Syntax: D=aremsf(A,B,C,eps,itermax)
//
// This routine solves the Algebraic Riccati Equation using the matrix
// Sign function.
//
//   A'P + PA - PBP + C = 0
//
// The input eps is the tolerance for the iteration procedure and the
// variable itermax is the maximum number of iterations.
//
// The routine returns the steady-state continuous ARE solution.
//
//  Ref.: (1) Chapter 5 of Junkins & Kim, Dyn. & Ctrl. of Structures.
//        (2) Bierman, G.J.,
//            "Computational Aspects of the Matrix Sign Function to the ARE"
//            Proc. 23rd CDC, 1984, pp.514-519.
//
// Copyright (C), by Jeffrey B. Layton, 1994
// Version JBL 940918
//--------------------------------------------------------------------------------
require abcdchk

aremsf = function(A,B,C,Reps,iterm)
{

   // Error Traps
   // ===========

   // Check if A and C are consistent.
   if ( A.nr != C.nr || A.nc != C.nc ) {
       error("AREMSF: A and C must be consistent");
   }
   
   if ( A.nr != A.nc || B.nc != A.nr ) {
       error("AREMSF: A and B must be consistent");
   }
   
   // Check if A and B are empty
   if ( (!length(A)) || (!length(B)) ) {
       error("AREMSF: A and B matrices cannot be empty.");
   }

   // A has to be square
   if (A.nr != A.nc) {
       error("AREMSF: A has to be square.");
   }

   // See if A and B are compatible.
   msg="";
   msg=abcdchk(A,B);
   if (msg != "") {
       estr="AREMSF: "+msg;
       error(estr);
   }

   // Check if C is positive semi-definite and symmetric
   if (!issymm(C)) {
       printf("aremsf: Warning: C is not symmetric.\n");
   else
       if (any(eig(C).val < -epsilon()*norm(C,"1")) ) { 
           printf("aremsf: Warning: C is not positive semi-definite.\n");
       }
   }

   // Set defaults if variables are not input.
   // ========================================

   if (nargs < 5) {
       itermax=1000;
   else
       itermax=iterm;
   }

   if (nargs < 4) {
       eps=1.0e-10;
   else
       eps=Reps;
   }

   n=A.nr;

   //  Iteration begins for computing matrix sign function of Hamiltonian
   //  (Use transformed Hamiltonian to get accurate inverse of Hamiltonian)

   Jhat=[zeros(n,n),-eye(n,n);eye(n,n),zeros(n,n)];
   Hhat=[C, A';A, -B];

   for (i in 1:itermax) {
        Hhatnew=0.5*(Hhat+Jhat*inv(Hhat)*Jhat);
        if (norm(Hhatnew-Hhat,"f") < eps) {
            break;
        }
        Hhat=Hhatnew;
   }

   SignH=-Jhat*Hhatnew;

   //  Compute the Riccati Solution

   W=0.5*(SignH+eye(2*n,2*n));
   W12=W[1:n;n+1:2*n];
   W11=W[1:n;1:n];
   P=-inv(W12)*W11;

   return P;
};

