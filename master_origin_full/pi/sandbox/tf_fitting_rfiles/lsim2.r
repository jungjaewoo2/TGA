//---------------------------------------------------------------------------------
//
// lsim -> lsim2 (no plotting
//
// Syntax: G=lsim(A,B,C,D,U,T,X0)
//
// This routine plots the time response of the linear system descried by,
//  		.
//  		x = Ax + Bu
//  		y = Cx + Du
//
// to the input time history U. The routine may be called in several fashions:
//
// (1)  G=lsim(A,B,C,D,U,T)
//      (This produces the time history with Zero initial conditions for
//       a state-space model).
//
// (2)  G=lsim(A,B,C,D,U,T,X0)
//      (This produces the time history with initial conditions for
//       a state-space model).
//
// (3)   G=lsim(NUM,DEN,U,T)
//       (This produces the time history with Zero initia conditions for
//        a transfer function model, G(s)=NUM(s)/DEN(s) ).
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
// Copyright(C), by Jeffrey B. Layton, 1994
// Version JBL 940906
//---------------------------------------------------------------------------------

require abcdchk tfchk tf2ss c2d ltitr

lsim2 = function(a,b,c,d,u,t,x0)
{
   global(_rlab_config,_ctb2_window)

   if (nargs < 4) {
       error("LSIM: Incorrect number of input arguments.");
   }

   // Check for T.F. case
   if (nargs == 4) {
       // Check if T.F. is proper
       Dum=tfchk(a,b);
       num=Dum.numc;
       den=Dum.denc;
       u=c;
       t=d;
       // Convert T.F. to S.S.
       Dum=tf2ss(num,den);
       A=Dum.a;
       B=Dum.b;
       C=Dum.c;
       D=Dum.d;
   else
       A=a;
       B=b;
       C=c;
       D=d;
       msg="";
       msg=abcdchk(A,B,C,D);
       if (msg != "") {
           estr="lsim: "+msg;
           error(estr);
       }
   }

   // If X0 doesn't exist, then create a zero vector for it
   if ( (nargs == 6) || (nargs == 4) ) {
       X0 = zeros(A.nr,1);
   else
       X0=x0;
   }

   // Perform dimension checks
   if (u.nc != D.nc) {
       error("LSIM: U has to have same number of cols. as inputs");
   }
   if (u.nr != length(t)) {
       error("LSIM: U must have same number of rows as length of T");
   }

   // Convert continuous to discrete using zero order hold.
   Dum=c2d(A,B,t[2]-t[1]);
   P=Dum.phi;
   G=Dum.gamma;

   // Propagate the simulation over the number of elements in t
   x = ltitr(P,G,u,X0);
   y = x * C.' + u * D.';

   #// Plot the results
   #if (_rlab_config.plot_support=="pgplot"||_rlab_config.plot_support=="plplot")
   #{
     #if (exist(_ctb2_window))
     ##{
        #plwin(_ctb2_window);
     #else
        #_ctb2_window = plstart();
     #}
   #}
   #xlabel("Time (secs)");
   #ylabel("Amplitude");
   #plot([t,y]);
   #xlabel("");
   #ylabel("");
   
   return << x=x; y=y >>;
};

