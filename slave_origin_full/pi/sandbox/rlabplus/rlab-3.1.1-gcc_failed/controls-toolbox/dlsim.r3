//-----------------------------------------------------------------------
//
// dlsim
//
// Syntax: </x;y/> = dlsim(A,B,C,D,U,X0)
//         G       = dlsim(A,B,C,D,U,X0) 
//
// This routine plots the time response of the linear system descried by,
//
//     x[n+1] = Ax[n] + Bu[n]
//     y[n]   = Cx[n] + Du[n]
//
// to the input sequence U. The routine may be called in several fashions:
//
// (1)  G=dlsim(A,B,C,D,U)
//      (This produces the sequence with Zero initial conditions for
//       a state-space model).
//
// (2)  G=dlsim(A,B,C,D,U,X0)
//      (This produces the sequence with initial conditions for
//       a state-space model).
//
// (3)   G=dlsim(NUM,DEN,U)
//       (This produces the sequence with Zero initial conditions for
//        a transfer function model, G(z)=NUM(z)/DEN(z) ).
//
// For all 3 cases, the matrix U must have as many columns as there
// are inputs, U. Each row of U corresponds to a new time point, and
// U must have length(T) rows. The input T is the time vector which
// must be regularly spaced.
//
// Note: Two matrices are returned in a list.
//
//       G.x = X values in the plot.
//       G.y = Y values in the plot.
//
// Note: G.x has as many columns as there are states and length(U) rows.
//       G.y has as many columns as there are outputs and length(U) rows.
//
// Copyright(C), by Jeffrey B. Layton, 1994
// Version JBL 940907
//-----------------------------------------------------------------------

require tfchk tf2ss abcdchk ltitr stairs

dlsim = function(a,b,c,d,u,x0)
{

// If nargs < 3, error
   if ( (nargs < 3) || (nargs == 4) ) {
       error("DLSIM:  Incorrect number of arguments.");
   }

// Check for T.F. case
   if (nargs == 3) {
// Check if T.F. is proper
       Dum=tfchk(a,b);
       num=Dum.numc;
       den=Dum.denc;
       u=c;
// Convert T.F. to S.S.
       Dum=tf2ss(num,den);
       A=Dum.a;
       B=Dum.b;
       C=Dum.c;
       D=Dum.d;
   } else {
       A=a;
       B=b;
       C=c;
       D=d;
       msg="";
       msg=abcdchk(A,B,C,D);
       if (msg != "") {
           estr="DLSIM: "+msg;
           error(estr);
       }
   }

// If X0 doesn't exist, then create a zero vector for it
   if ( (nargs == 5) || (nargs == 3) ) {
       X0 = zeros(1,A.nr);
   } else {
       X0=x0;
   }

// Check Dimensions
   if (u.nc != D.nc) {
       error("DLSIM: U must have same number of cols. as inputs.");
   }

// Perform the simulation by using ltitr
   if (isempty(A) == 1) {
       x=[];
       y=u*D.';
   } else {
       x=ltitr(A,B,u,X0);
       y=x*C.'+u*D.';
   }

// Make the Stairs plot
   t = [0:length(y)-1];
   xlabel("No. of Samples");
   ylabel("Amplitude");
   Dum=stairs(t,y);

   return << x=x; y=y >>;
};

