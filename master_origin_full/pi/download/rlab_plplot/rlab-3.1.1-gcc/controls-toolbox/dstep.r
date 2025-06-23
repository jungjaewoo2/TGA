//-----------------------------------------------------------------------------
//
// dstep
//
// G       = dstep(A,B,C,D,IU,N)
// </x;y/> = dstep(A,B,C,D,IU,N)
//
// This routine plots the unit step response of discrete-time linear systems.
// It plots the time response of the following linear discrete-time
// system
//
//     x[n+1] = Ax[n] + Bu[n]
//     y[n]   = Cx[n] + Du[n]
//
// to a unit step input applied to the input IU. The routine may be called
// in several fashions:
//
// (1)  G=dstep(A,B,C,D)
//      (This produces the dstep response of the system looping over the
//       inputs and looping over the outputs with an arbitrary sample
//       of sample points N=200).
//
// (2)  G=dstep(A,B,C,D,IU)
//      (This produces the dstep response of the system to the IU input
//       and looped over the outputs with an arbitrary number of sample
//       points N=200)
//
// (3)  G=dstep(A,B,C,D,IU,N)
//      (This produces the dstep response of the system to the IU input
//       and looped over the outputs for the specified number of sample
//       points N).
//
// (4)  G=dstep(NUM,DEN)
//      (This produces the step response of the transfer function model
//       to the IU input with an arbitrary number of sample points N=200)).
//
// (5)  G=dstep(NUM,DEN,T)
//      (This produces the step response of the transfer function model
//       with the specified number of samples points N).
//
// Note: Two matrices are returned in a list.
//
//       G.x = X values in the plot.
//       G.y = Y values in the plot.
//
// Note: The matrix G.y has as many columns as there are outputs and
//       has N rows. The matrix G.x has as many columns as states
//       and has N rows).
//
// Copyright(C), by Jeffrey B. Layton, 1994
// Version JBL 940915
//-----------------------------------------------------------------------------
require tfchk tf2ss abcdchk dlsim

dstep = function(a,b,c,d,iu,n)
{

// System check
// ------ -----
   if (nargs == 2) {
// T.F. no N
       Dum=tfchk(a,b);
       numc=Dum.numc;
       denc=Dum.denc;
       IU=1;
       Dum=tf2ss(numc,denc);
       A=Dum.a;
       B=Dum.b;
       C=Dum.c;
       D=Dum.d;
// Set default N
       N=200;
   else if (nargs == 3) {
// T.F. N specified
       Dum=tfchk(a,b);
       numc=Dum.numc;
       denc=Dum.denc;
       N=c;
       IU=1;
       Dum=tf2ss(numc,denc);
       A=Dum.a;
       B=Dum.b;
       C=Dum.c;
       D=Dum.d;
   else if (nargs >= 4) {
// State-Space
       A=a;
       B=b;
       C=c;
       D=d;
       msg="";
       msg=abcdchk(A,B,C,D);
       if (msg != "") {
           estr="DSTEP: "+msg;
           error(estr);
       }
       if (exist(iu)) {
           IU=iu;
           if (exist(n)) {
               N=n;
           else
               N=200;
           }
       else
           IU=B.nc;
           N=200;
       }
   }}}


// Perform the Simulation
   for (i in 1:IU) {

// Define B and D for ith input
        BI=B[;i];
        DI=D[;i];
        Dum=dlsim(A,BI,C,DI,ones(N,1));
        x=Dum.x;
        y=Dum.y;
   }

   return << x=x; y=y >>;
};

