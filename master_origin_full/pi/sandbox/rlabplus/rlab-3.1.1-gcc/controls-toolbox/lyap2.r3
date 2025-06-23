//-----------------------------------------------------------------------
//
// lyap2
//
// SyntaX: x=lyap2(A,B,C)
//
// This routine solves the lyapunov equation using an eigenvalue
// decomposition (spectral decomposition). The solution to the
// equation is returned in a matrix X (X=lyap2(A,B,C) ).
//
// If the routine is called with, X=lyap2(A,B) then it solves the
// particular form of the Lyapunov equation:
//
//       A*X+X*A' = -B
//
// Calling the routine with X=lyap(A,B,C) solves the general form
// of the lyapunov equation (sylvester equation):
//
//       A*X+X*B = -C
//
// ==================================================================
// Note: This routine assumes that there are no repeated eigenvalues.
//       If repeated eigenvalues exist, then the solution may be
//       unreliable.
// ==================================================================
//
// Ref: Junkins, J., and Kim, Y., "Introduction to Dynamics and Control
//      of Flexible Structures," AIAA Inc, Washington D.C., 1993,
//      Chapter 2.
//
// Copyright (C), by Jeffrey B. Layton, 1994
// Version JBL 940919
//-----------------------------------------------------------------------

lyap2 = function(A,B,C)
{

   if (nargs != 3 && nargs !=2) {
       error("LYAP2: Wrong number of input arguments.");
   }

// Dimension Check
   if (A.nr != A.nc) { error("lyap2: Dimensions of A and B do not agree"); }
   if (B.nr != B.nc) { error("lyap2: Dimensions of A and B do not agree"); }
   if (A.nr == 0) {
       X = [];
       return X;
   }

// Compute eigenvalues and the left and right eigenvectors of A
   E=eig(A);
   LamA=diag(E.val);
   DA=LamA;
   PhiA=E.vec;
   PsiA=eig(A').vec;

// Switch on Lyapunov solution or Sylvester solution.
   if ( exist (C)) {
       if (C.nr != A.nr) { error("lyap2: Dimensions of C do not agree"); }
       if (C.nc != B.nr) { error("lyap2: Dimensions of C do not agree"); }

// Computes eigenvalues and left and right eigenvectors of B
       E=eig(B);
       LamB=diag(E.val);
       DB=LamB;
       LamB=ones(A.nr,B.nr)*LamB;
       PhiB=E.vec;
       PsiB=eig(B').vec;
       LamA=LamA*ones(A.nr,B.nr);

       T=-(PsiA.'*C*PhiB) ./ (LamA+LamB);
       X=PhiA*T*PhiB';
   } else {
       LamA=LamA*ones(A.nr,A.nr);
       T=-(PsiA.'*B*PsiA) ./ (LamA+LamA.');
       X=PhiA*T*PhiA';
   }
   X=real(X);

   return X;
};

