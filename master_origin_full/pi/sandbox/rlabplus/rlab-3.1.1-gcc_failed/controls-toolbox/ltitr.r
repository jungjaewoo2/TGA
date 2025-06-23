//---------------------------------------------------------------------------
//
// ltitr
//
// Syntax: X=ltitr(A,B,U,X0)
//
// This routine computes the time response of the following system:
//
//  x[n+1] = Ax[n] + Bu[n]
//
// to the input U. The result is returned in a matrix X which has as
// many columns as there are outputs y (and with length(U) rows).
//
// This routine can also be used with initial conditions as:
//
//  X=ltitr(A,B,U,X0)
//
// Note: The input U must have as many columns as there are inputs u
// in the system. Each row of U corresponds to a new time point.
//
// Ref: Mathworks, "Control System Toolbox User's Guide," 1992.
//
// Copyright (C), by Jeffrey B. Layton, 1994
// Version JBL 940914
//---------------------------------------------------------------------------

ltitr = function(a,b,u,x0)
{

   U=u.';

   // Get problem dimensions
   n=U.nc;

   // Initialize x to the correct size
   x[a.nr;n]=0;

   // Check to see if x0 was input. If it wasn't make x0, zero.
   if (nargs == 3) {
       X0=zeros(a.nr,1);
   else
       X0=x0;
   }

   // Start computing
   if (nargs == 4) {
       X0=X0[:];
   }
   X0=X0[:];
   for (i in 1:n) {
        x[;i]=X0;
        X0=a*X0+b*U[;i];
   }
   x=x.';

   return x;
};

