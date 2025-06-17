//--------------------------------------------------------------------------------
// lqrmsf
//
// Syntax: D=lqrmsf(A,B,Q,R,CT,Reps,iterm)
//
// This routine solves the Algebraic Riccati Equation using the matrix
// Sign function. It also computes the optimal feedback gain matrix such that
// the feedback law u = -Kx  minimizes the following cost function:
//
//           Inf
//     J = Integral (x'Qx + u'Ru) dt
//            0
//
// subject to the constraint
//   .
//   x = Ax + Bu
//
// The solution P solves the following equation:
//
//                 -1
//   A'P + PA - PBR  B'P + Q = 0
//
// When the term CT is included, a Cross-Term is added that relates u to x
// in the following cost function:
//
//           Inf
//     J = Integral (x'Qx + u'Ru + 2*x'CTu) dt
//            0
//
// The input Reps is the tolerance for the iteration procedure and the
// variable iterm is the maximum number of iterations.
//
// Note: Three matrices are returned in a list.
//
//       D.G = Optimal Feedback Gain matrix
//       D.P = Steady-State Solution to the Algebraic Riccatti Eqn.
//
//
//  Ref.: (1) Chapter 5 of Junkins & Kim, Dyn. & Ctrl. of Structures.
//        (2) Bierman, G.J.,
//            "Computational Aspects of the Matrix Sign Function to the ARE"
//            Proc. 23rd CDC, 1984, pp.514-519.
//
// Copyright (C), by Jeffrey B. Layton, 1994
// Version JBL 940618
//--------------------------------------------------------------------------------

require abcdchk

static (lqrmsf_compute)    // Hide this function

lqrmsf = function(a,b,q,r,ct,Reps,iterm)
{
   local(nargs,itermax,eps,msg,estr)

// Count number of input arguments
   nargs=0;
   if (exist(a)) {nargs=nargs+1;}
   if (exist(b)) {nargs=nargs+1;}
   if (exist(q)) {nargs=nargs+1;}
   if (exist(r)) {nargs=nargs+1;}

// Error Traps
// ===========

// Check if Q and R are consistent.
   if ( a.nr != q.nr || a.nc != q.nc ) {
       error("LQRMSF: A and Q must be the same size");
   }
   
   if ( r.nr != r.nc || b.nc != r.nr ) {
       error("LQRMSF: B and R must be consistent");
   }
   
// Check if A and B are empty
   if ( (!length(a)) || (!length(b)) ) {
       error("LQRMSF: A and B matrices cannot be empty.");
   }

// A has to be square
   if (a.nr != a.nc) {
       error("LQRMSF: A has to be square.");
   }

// See if A and B are compatible.
   msg="";
   msg=abcdchk(a,b);
   if (msg != "") {
       estr="LQRMSF: "+msg;
       error(estr);
   }

// Check if Q is positive semi-definite and symmetric
   if (!issymm(q)) {
       printf("LQRMSF: Warning: Q is not symmetric.\n");
   } else {
       if (any(eig(q).val < -epsilon()*norm(q,"1")) ) { 
           printf("LQRMSF: Warning: Q is not positive semi-definite.\n");
       }
   }

// Check if R is positive definite and symmetric
   if (!issymm(r)) {
       printf("LQRMSF: Warning: R is not symmetric.\n");
   } else {
       if (any(eig(r).val < -epsilon()*norm(r,"1")) ) {
           printf("LQRMSF: Warning: R is not positive semi-definite.\n");
       }
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

//
// Call lqrmsf_compute with revised arguments. This prevents us from needlessly
// making a copy if ct does not exist, and also prevents us from overwriting
// A and Q if CT exists.
//

// If Cross-Term exists, modify Q and A.
   if (exist(ct)) {
       if ( ct.nr != a.nr || ct.nc != r.nc) {
            error("LQRMSF: CT must be consistent with Q and R");
       }
       return lqrmsf_compute (a-solve(r',b')'*ct', b, ...
                              q-solve(r',ct')'*ct', r, ct, ...
                              eps, itermax);
   } else {
       return lqrmsf_compute (a, b, q, r, zeros (a.nr, b.nc), ...
                              eps, itermax);
   }
};

//----------------------------------------------------------------------------
//
// This is where the computation is performed. Note that lqrmsf_compute is a
// static variable and is never seen from the global workspace.
//

lqrmsf_compute = function (A, B, Q, R, CT, eps, itermax)
{
   local(n,Jhat,Hhat,Hhatnew,SignH,W,W12,W11,errstr)

   n=A.nr;

//  Iteration begins for computing matrix sign function of Hamiltonian
//  (Use transformed Hamiltonian to get accurate inverse of Hamiltonian)

   Jhat=[zeros(n,n),-eye(n,n);eye(n,n),zeros(n,n)];
   Hhat=[Q, A';A,-B*inv(R)*B'];

   for (i in 1:itermax) {
        Hhatnew=0.5*(Hhat+Jhat*inv(Hhat)*Jhat);
        if (norm(Hhatnew-Hhat,"f") < eps) {
            break;
        }
        Hhat=Hhatnew;
   }

// Check iteration counter
   if (i == itermax) {
       errstr="lqrmsf: failed to converge in "+int2str(itermax)+" iterations"
       printf("%s",errstr);
   }

   SignH=-Jhat*Hhatnew;

//  Compute the Riccati Solution

   W=0.5*(SignH+eye(2*n,2*n));
   W12=W[1:n;n+1:2*n];
   W11=W[1:n;1:n];
   P=-inv(W12)*W11;
// G=inv(R)*B'*P;
   G=solve( R, CT'+B'*P );

   return << G=G; P=P >>
};

