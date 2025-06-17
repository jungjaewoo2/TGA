//-------------------------------------------------------------------------------
//
// dlqr2
//
// Syntax: A = dlqr2(a,b,q,r,ct)
//         where A.e = Closed-Loop Eigenvalues
//               A.k = Optimal Feedback Gain Matrix.
//               A.p = Steady-State Solution to the Algebraic Riccatti Eqn.
//
// This routine computes the Linear Quadratic Regulator Design for Discrete
// Systems. It computes the optimal feedback gain matrix K such that the
// feedback law, u[n] = -Kx[n]  minimizes the following cost function:
//
//       J = Sum {x'Qx + u'Ru}
//
// subject to the constraint equation:   
//
//       x[n+1] = Ax[n] + Bu[n] 
//                
// The solution P solves the steady-state discrete Riccatti equation:
//
//                              -1
//       P - A'PA + A'PB(R+B'PB)  BP'A - Q = 0
//
// The routine returns the gains (G.k), the steady-state solution to the
// discrete Riccatti equation (G.p), and the eigenvalues ofthe closed loop
// system eig(A-B*K) (G.e).
//
// When the term CT is included, a Cross-Term is added to the cost
// function which relates u to x in the following cost function:
//
//       J = Sum {x'Qx + u'Ru + 2*x'Nu}
//
// Note: Three matrices are returned in a list.
//
//      A.k Optimal Feedback Gain Matrix.
//      A.p Steady-State Solution to the Algebraic Riccatti Eqn.
//      A.e Closed-Loop Eigenvalues
//
// Ref: (1) Laub, A. J., "A Schur Method for Solving Algebraic Riccati
//          Equations," IEEE Trans. AC-24, 1979, pp. 913-921.
//
// Copyright (C), by Jeffrey B. Layton, 1994
// Version JBL 940918
//-------------------------------------------------------------------------------

require abcdchk

static (dlqr2_compute)    // Hide this function

dlqr2 = function(a,b,q,r,ct)
{

   if (nargs < 4) {
       error("DLQR2: Wrong number of arguments.");
   }

// Check compatibility of A and B matrices
   msg="";
   msg=abcdchk(a,b);
   if (msg != "") {
       estr="DLQR2: "+msg;
       error(estr);
   }

// Check if A and B are empty
   if ( (!length(a)) || (!length(b)) ) {
       error("DLQR2: A and B matrices cannot be empty.");
   }

// Check of Q is if proper dimensions
   if ( (a.nr != q.nr) || (a.nc != q.nc) ) {
	error("DLQR2: A and Q must be the same size");
   }

// Check if R is of proper dimensions
   if ( (r.nr != r.nc) || (b.nc != r.nr) ) {
	error("DLQR2: B and R must be consistent");
   }

// Check if Q is positive semi-definite and symmetric
   if (!issymm(q)) {
       printf("DLQR2: Warning: Q is not symmetric.\n");
   } else {
       if (any(eig(q).val < -epsilon()*norm(q,1)) ) {
           printf("DLQR2: Warning: Q is not positive semi-definite.\n");
       }
   }

// Check if R is positive definite and symmetric
   if (!issymm(r)) {
       printf("DLQR2: Warning: R is not symmetric.\n");
   } else {
       if (any(eig(r).val < -epsilon()*norm(r,1)) ) {
           printf("DLQR2: Warning: R is not positive semi-definite.\n");
       }
   }

//
// Call dlqr2_compute with revised arguments. This prevents us from needlessly
// making a copy if ct does not exist, and also prevents us from overwriting
// a and q if ct exists.
//
 
   if (exist(ct)) {
       if ( ct.nr != a.nr || ct.nc != r.nc) {
            error("DLQR2: CT must be consistent with Q and R");
       }
       return dlqr2_compute (a-solve(r',b')'*ct', b, ...
                           q-solve(r',ct')'*ct', r, ct);
   } else {
       return dlqr2_compute (a, b, q, r, zeros (a.nr, b.nc));
   }
};

//----------------------------------------------------------------------------
//
// This is where the computation is performed. Note that dlqr_compute is a
// static variable and is never seen from the global workspace.
//

dlqr2_compute = function (a, b, q, r, ct)
{
   local(G,A,v,d,index,e,Phi12,Phi22,p,k)

// Form G matrix (See the Ref.)
   G=solve(r',b')'*b';

// Form the Hamiltonian and compute the eigenvectors
   A=eig([a+G/a'*q,-G/a';-a'\q, inv(a)']);
   v=A.vec;
   d=A.val;

// Sort Eigenvectors by sorting eigenvalues and storing the index.
   d=diag(d);
   index=sort(abs(d));   // sort on magnitude of eigenvalues
   e=real(d[ index ] );
   if ( !((e[a.nc] < 1) && (e[a.nc+1] > 1) ) ) {
	error("DLQR2: Can't order eigenvalues, (A,B) may be uncontrollable.");
   } else {
        e=d[index[1:a.nc]];
   }
// select vectors with eigenvalues inside unit circle
   Phi12=v[1:a.nc,index[1:a.nc]];
   Phi22=v([a.nc+1]:[2*a.nc],index[1:a.nc]);
   p=real(Phi22/Phi12);
   k=(r+b'*s*b)\b'*s*a + r\ct';

   return << k=k; p=p; e=e >>;
};

