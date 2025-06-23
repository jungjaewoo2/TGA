//----------------------------------------------------------------------------
//
// gric
//
// Syntax: DP=gric(Q,DA,B,DB,Q,DQ,R,DR)
//
// This routine computes the gradient of the solution of the Algebraic
// Riccatti Equation (ARE) with respect to a scalar parameter 'p'
// for Continuous systems. The ARE is,
//
//                 -1
//   A'P + PA - PBR  B'P + Q = 0
//
// The gradient of the input matrices (A,B,Q, and R) with respect to 'p'
// are required.
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

require abcdchk rica

static (gric_compute)	// Hide this function

gric = function(a,da,b,db,q,dq,r,dr)
{

   if (nargs != 8) {
      error("GRIC: Wrong number of input arguments.");
   }

// Check Dimensions
// ----------------

// Check if Q and R are consistent.
   if ( a.nr != q.nr || a.nc != q.nc ) {
       error("GRIC: A and Q must be the same size");
   }
   
   if ( r.nr != r.nc || b.nc != r.nr ) {
       error("GRIC: B and R must be consistent");
   }
   
// Check if A and B are empty
   if ( (!length(a)) || (!length(b)) ) {
       error("GRIC: A and B matrices cannot be empty.");
   }

// A has to be square
   if (a.nr != a.nc) {
       error("GRIC: A has to be square.");
   }

// Check gradient matrices
// -----------------------

// Check if DQ and DR are consistent.
   if ( da.nr != dq.nr || da.nc != dq.nc ) {
       error("GRIC: DA and DQ must be the same size");
   }
   
   if ( dr.nr != dr.nc || db.nc != dr.nr ) {
       error("GRIC: DB and DR must be consistent");
   }
   
// Check if DA and DB are empty
   if ( (!length(da)) || (!length(db)) ) {
       error("GRIC: DA and DB matrices cannot be empty.");
   }

// DA has to be square
   if (da.nr != da.nc) {
       error("GRIC: DA has to be square.");
   }

// See if A and B are compatible.
   msg="";
   msg=abcdchk(a,b);
   if (msg != "") {
       estr="GRIC: "+msg;
       error(estr);
   }

// Check if Q is positive semi-definite and symmetgric
   if (!issymm(q)) {
       printf("GRIC: Warning: Q is not symmetric.\n");
   } else {
       if (any(eig(q).val < -epsilon()*norm(q,"1")) ) {	
           printf("GRIC: Warning: Q is not positive semi-definite.\n");
       }
   }

// Check if R is positive definite and symmetgric
   if (!issymm(r)) {
       printf("GRIC: Warning: R is not symmetric.\n");
   } else {
       if (any(eig(r).val < -epsilon()*norm(r,"1")) ) {
           printf("GRIC: Warning: R is not positive semi-definite.\n");
       }
   }

//
// Call gric_compute to solve Riccatti.
//
    return gric_compute (a,da,b,db,q,dq,r,dr);
};

//----------------------------------------------------------------------------
//
// This is where the computation is performed. Note that gric_compute is a
// static variable and is never seen from the global workspace.
//

gric_compute = function (a,da,b,db,q,dq,r,dr)
{
   local (p, abar, bbar, cbar, term1, term2, term3, dp)

// Solve Riccatti by calling rica
   p=rica(a,b,q,r);

// Create matrices for solving sylvester equation.

// A bar:
   abar = a' - p*(b/r)*b';

// B bar:
   bbar = a - (b/r)*b'*p;

// C bar:
   term1 = p*(b/r)*db'*p;
   term2 = p*(b/r)*(dr/r)*b'*p;
   term3 = p*(db/r)*b'*p;
   cbar = da'*p + p*da - term3 + term2 - term1 + dq;

// Solve Sylvester's Equation (using lyap.r)
   dp=lyap(abar,bbar,cbar);
   
   return dp;
};

