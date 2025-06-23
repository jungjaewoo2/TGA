//----------------------------------------------------------------------------
//
// glqr
//
// Syntax: G = glqr(A,DA,B,DB,Q,DQ,R,DR,CT,DCT)
//         where G.dk = Gradient of Optimal Feedback Gain matrix
//               G.dp = Gradient of Steady-State Solution to ARE.
//
// This routine computes the gradient of the LQR control design with respect
// to a scalar parameter 'p' for Continuous Systems. It returns the gradient
// of the optimal feedback gain matrix with respect to 'p'. It also returns
// the gradient of the Riccatti Equation:
//
//                 -1
//   A'P + PA - PBR  B'P + Q = 0
//
// Note: Two matrices are returned in a list.
//
//       G.dk = Gradient of Optimal Feedback Gain matrix
//       G.dp = Gradient of Steady-State Solution to ARE.
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

require abcdchk rica gric

static (glqr_compute)	// Hide this function

glqr = function(a,da,b,db,q,dq,r,dr,ct,dct)
{

   if ( (nargs != 8) || (nargs != 10) ) {
       error("DLQR: Wrong number of input arguments.");
   }

// Check Matrix Dimensions
// -----------------------

// Check if Q and R are consistent.
   if ( a.nr != q.nr || a.nc != q.nc ) {
       error("GLQR: A and Q must be the same size");
   }
   
   if ( r.nr != r.nc || b.nc != r.nr ) {
       error("GLQR: B and R must be consistent");
   }
   
// Check if A and B are empty
   if ( (!length(a)) || (!length(b)) ) {
       error("GLQR: A and B matrices cannot be empty.");
   }

// A has to be square
   if (a.nr != a.nc) {
       error("GLQR: A has to be square.");
   }

// Check gradient matrices
// -----------------------

// Check if DQ and DR are consistent.
   if ( da.nr != dq.nr || da.nc != dq.nc ) {
       error("GLQR: DA and DQ must be the same size");
   }
   
   if ( dr.nr != dr.nc || db.nc != dr.nr ) {
       error("GLQR: DB and DR must be consistent");
   }
   
// Check if DA and DB are empty
   if ( (!length(da)) || (!length(db)) ) {
       error("GLQR: DA and DB matrices cannot be empty.");
   }

// DA has to be square
   if (da.nr != da.nc) {
       error("GLQR: DA has to be square.");
   }

// See if A and B are compatible.
   msg="";
   msg=abcdchk(a,b);
   if (msg != "") {
       estr="GLQR: "+msg;
       error(estr);
   }

// Check if Q is positive semi-definite and symmetric
   if (!issymm(q)) {
       printf("GLQR: Warning: Q is not symmetric.\n");
   else
       if (any(eig(q).val < -epsilon()*norm(q,"1")) ) {	
           printf("GLQR: Warning: Q is not positive semi-definite.\n");
       }
   }

// Check if R is positive definite and symmetric
   if (!issymm(r)) {
       printf("GLQR: Warning: R is not symmetric.\n");
   else
       if (any(eig(r).val < -epsilon()*norm(r,"1")) ) {
           printf("GLQR: Warning: R is not positive semi-definite.\n");
       }
   }

//
// Call glqr_compute with revised arguments. This prevents us from needlessly
// making a copy if ct does not exist, and also prevents us from overwriting
// a and q if ct exists.
//

// Call static routine to compute results
   if (exist(ct)) {
       if ( ct.nr != a.nr || ct.nc != r.nc) {
            error("GLQR: CT must be consistent with Q and R");
       }
       if ( dct.nr != da.nr || dct.nc != dr.nc) {
            error("GLQR: DCT must be consistent with DQ and DR");
       }
       return glqr_compute(a,da,b,db,q,dq,r,dr,ct,dct);
//     return glqr_compute (a-solve(r',b')'*ct', b, ...
//                         q-solve(r',ct')'*ct', r, ct);
   else
       return glqr_compute(a,da,b,db,q,dq,r,dr,zeros(a.nr,b.nc),zeros(a.nr,b.nc));
//     return glqr_compute (a, b, q, r, zeros (a.nr, b.nc));
   }
};

//----------------------------------------------------------------------------
//
// This is where the computation is performed. Note that glqr_compute is a
// static variable and is never seen from the global workspace.
//

glqr_compute = function(a,da,b,db,q,dq,r,dr,ct,dct)
{
   local (alocal, qlocal, dalocal, dqlocal, k, p, dp, dk)

// Modify the A and Q matrices
   if (exist(ct)) {
       alocal=a-solve(r',b')'*ct';
       qlocal=q-solve(r',ct')'*ct';
   else
       alocal=a;
       qlocal=q;
   }

// Solve the Riccatti equation
   p=rica(alocal,b,qlocal,r);

// Create matrices for input to gric
   if (exist(ct)) {
       dqlocal=dq-(dct/r)*ct' + (ct/r)*(dr/r)*ct' - (ct/r)*ct';
       dalocal=da-(db/r)*ct' + (b/r)*(dr/r)*ct' - (b/r)*ct';
   else
       dalocal=da;
       dqlocal=dq;
   }

// Compute Gradient of Riccatti solution
   dp = gric(a,dqlocallocal,b,db,q,dq,r,dr);

// Compute gradient of the gains
   if (exist(ct)) {
       dk=-(r\dr)*(r\(ct'+b'*p)) + r\(dct'+db'*p+b'*dp);
   else
       dk=-(r\dr)*(r\(b'*p)) + r\(db'*p) + r\(b'*dp);
   }
   
   return << dk=dk; dp=dp >>;
};

