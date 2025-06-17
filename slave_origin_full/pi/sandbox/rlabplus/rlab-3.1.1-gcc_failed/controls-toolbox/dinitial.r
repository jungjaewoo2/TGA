//--------------------------------------------------------------------------------
//
// dinitial
//
// Syntax: G=dinitial(a,b,c,d,x0,n)
//
// This routine plots the time response of the discrete linear system descried by,
//
//     x[n+1] = Ax[n] + Bu[n]
//     y[n]   = Cx[n] + Du[n]
//
// to an initial condition on the states. The routine may be called in two ways:
//
// (1)  G=dinitial(A,B,C,D)
//      (This produces the sequence with random initial conditions for
//       a state-space model using the default 200 sample points).
//
// (2)  G=dinitial(A,B,C,D,X0)
//      (This produces the sequence with initial conditions for
//       a state-space model using the default 200 sample points).
//
// (3)  G=dinitial(A,B,C,D,X0,N)
//      (This produces the sequence with initial conditions for
//       a state-space model using the user specified number of sample
//       points N).
//
// All 3 cases return the state and output responses. There are as
// many columns in the output response (G.y) as there are outputs.
// The state response matrix (G.x) has as many columns as there are
// states.
//
// Note: Two matrices are returned in a list.
//
//       G.x = X values in the plot.
//       G.y = Y values in the plot.
//
// Note: G.x has as many columns as there are states and N rows.
//       G.y has as many columns as there are outputs and N rows.
//
// Note: The default number of sample points is 200.
//
// Copyright(C), by Jeffrey B. Layton, 1994
// Version JBL 940912
//--------------------------------------------------------------------------------

require abcdchk  dlsim

dinitial = function(a,b,c,d,x0,n)
{
 
// Check number of arguments
   if (nargs < 4) {
       error("DINITIAL: Wrong number of input arguments.");
   }

// Check A,B,C, and D for compatibility
   msg="";
   msg=abcdchk(a,b,c,d);
   if (msg != "") {
       estr="lsim: "+msg;
       error(estr);
   }

// Check is X0 is input
   if (nargs == 4) {
       X0=rand(a.nr,1);
   else
       X0=x0;
       X0=X0[:];
   }

// Check if n was input
   if ( !exist(n) ) {
       N=200;
   else
       N=n;
   }

// Perform Simulation
   Dum=dlsim(a,b,c,d,zeros(N,b.nc),X0[:]);
   y=Dum.y;
   x=Dum.x;

   return << x=x; y=y >>
};

