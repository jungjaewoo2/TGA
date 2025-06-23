//---------------------------------------------------------------------------
//
// initial
//
// Syntax: G=initial(A,B,C,D,X0,T)
//
// This routine plots the time response of the linear system descried by,
//              .
//              x = Ax + Bu
//              y = Cx + Du
//
// to an initial condition on the states. The routine may be called in two ways:
//
// (1)  G=initial(A,B,C,D)
//      (This produces the time history with random initial conditions
//       and with an arbitary time T=0.0:20.0:0.1)
//
// (2)  G=initial(A,B,C,D,X0)
//      (This produces the time history with the given initial conditions
//       and with an arbitrary time T=0.0:20.0:0.1)
//
// (3)  G=initial(A,B,C,D,X0,T)
//      (This produces the time history with the given initial conditions
//       with the given time in T).
//
// For the third case, the times in the time vector must be regularly
// spaced.
//
// Note: Two matrices are returned in a list.
//
//       G.x = X values in the plot.
//       G.y = Y values in the plot.
//
// Copyright(C), by Jeffrey B. Layton, 1994
// Version JBL 940912
//---------------------------------------------------------------------------

require abcdchk c2d ltitr lsim

initial = function(a,b,c,d,x0,t)
{

// Check number of inputs
   if (nargs < 4) {
       error("INITIAL: Error in number of input arguments");
   }

// Check if A,B,C, and D are compatibile
   msg="";
   msg=abcdchk(a,b,c,d);
   if (msg != "") {
       estr="lsim: "+msg;
       error(estr);
   }

// No X0 specified, try random initial conditions
   if (nargs == 4) {
       X0=rand(a.nr,1);
   } else {
       X0=x0[:];
   }

// If no T is specified, then create default T
   if (nargs == 5) {
       T=[0.0:20.0:0.1];
       T=T.';
   }

// Check Time vector (if input)
   if (nargs == 6) {
       if (t.nr == 1) {
           T=t.';
       } else {
           T=t;
       }
   }

// Perform the simulation
   n=length(T);
   Dum=lsim(a,b,c,d,zeros(n,b.nc),T,X0[:]);
   x=Dum.x;
   y=Dum.y;

   return << x=x; y=y >>;
};

