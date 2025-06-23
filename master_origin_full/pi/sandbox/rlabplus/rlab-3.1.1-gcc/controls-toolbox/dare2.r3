//-------------------------------------------------------------------------------
//
// dare2
//
// Syntax: A=dare2(a,b,q,r,ct)
//
// This routine solve the steady-state Discrete Algebraic Riccati Equation.
//
//                              -1
//       P - A'PA + A'PB(R+B'PB)  BP'A - Q = 0
//
// The routine returns the steady-state Discrete ARE solution.
//
// Ref: (1) Laub, A. J., "A Schur Method for Solving Algebraic Riccati
//          Equations," IEEE Trans. AC-24, 1979, pp. 913-921.
//      (2) Kailath, "Linear Systems," Prentice-Hall, 1980.
//      (3) Pappas, T., Laub, A. J., and Sandell, N. R., "On the Numerical
//          Solution of the Discrete-Time Algebraic Riccati Equation," IEEE
//          Trans. AC-25, 1980, pp. 631-641.
//
// Copyright (C), by Jeffrey B. Layton, 1994
// Version JBL 940918
//-------------------------------------------------------------------------------

require abcdchk

static (dare2_compute)    // Hide this function

dare2 = function(a,b,q,r,ct)
{

   if (nargs < 4 ) {
       error("DARE2: Wrong number of inputs");
   }

// Check compatibility of A and B matrices
   msg="";
   msg=abcdchk(a,b);
   if (msg != "") {
       estr="DARE2: "+msg;
       error(estr);
   }

// Check if A and B are empty
   if ( (!length(a)) || (!length(b)) ) {
       error("DARE2: A and B matrices cannot be empty.");
   }

// Check of Q is if proper dimensions
   if ( (a.nr != q.nr) || (a.nc != q.nc) ) {
	error("DARE2: A and Q must be the same size");
   }

// Check if R is of proper dimensions
   if ( (r.nr != r.nc) || (b.nc != r.nr) ) {
	error("DARE2: B and R must be consistent");
   }

// Check if Q is positive semi-definite and symmetric
   if (!issymm(q)) {
       printf("dare2: Warning: Q is not symmetric.\n");
   } else {
       if (any(eig(q).val < -epsilon()*norm(q,"1")) ) {
           printf("dare2: Warning: Q is not positive semi-definite.\n");
       }
   }

// Check if R is positive definite and symmetric
   if (!issymm(r)) {
       printf("dare2: Warning: R is not symmetric.\n");
   } else {
       if (any(eig(r).val < -epsilon()*norm(r,"1")) ) {
           printf("dare2: Warning: R is not positive semi-definite.\n");
       }
   }

//
// Call dare2_compute with revised arguments. This prevents us from needlessly
// making a copy if ct does not exist, and also prevents us from overwriting
// a and q if ct exists.
//

   if (exist(ct)) {
       if ( ct.nr != a.nr || ct.nc != r.nc) {
            error("dare2: CT must be consistent with Q and R");
       }
       return dare2_compute (a-solve(r',b')'*ct', b, ...
                             q-solve(r',ct')'*ct', r, ct, a.nc);
   } else {
       return dare2_compute (a, b, q, r, zeros (a.nr, b.nc), a.nc );
   }
};

//----------------------------------------------------------------------------
//
// This is where the computation is performed. Note that dare2_compute is a
// static variable and is never seen from the global workspace.
//

dare2_compute = function (a, b, q, r, ct, n)
{
   local(G,A,v,d,index,e,Phi12,Phi22,p,Adum)

// Form G matrix (See the Ref.)
   G=solve(r',b')'*b';

// Form the Hamiltonian and compute the eigenvectors
   A=eig([a+G/a'*q,-G/a';-a'\q, inv(a)']);
   v=A.vec;
   d=A.val;

// Sort Eigenvectors by sorting eigenvalues and storing the index.
   Adum=sort(abs(d));   // sort on magnitude of eigenvalues
   index=Adum.idx;
   e=Adum.val;
   if ( !((e[n] < 1) && (e[n+1] > 1) ) ) {
	error("dare2: Can't order eigenvalues, (A,B) may be uncontrollable.");
   } else {
        e=d[index[1:n]];
   }
// select vectors with eigenvalues inside unit circle
   Phi12=v[1:n,index[1:n]];
   Phi22=v[[n+1]:[2*n],index[1:n]];
   p=real(Phi22/Phi12);

   return p;
};

