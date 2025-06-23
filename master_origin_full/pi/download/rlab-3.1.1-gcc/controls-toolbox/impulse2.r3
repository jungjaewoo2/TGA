//-----------------------------------------------------------------------------
//
// impulse
//
// Syntax: impulse(a,b,c,d,ntim,iu)
//
// This routine finds the impulse response of a continuous-time
// linear system given by,
//         .
//         x = Ax + Bu
//         y = Cx + Du
//
// Calling the routine as,
//
//     impulse(A,B,C,D,NTIM)
//
// plots the time response of the linear system to an impulse
// response applied to each input. The variable NTIM is the number
// of time steps to be taken in the response.
//
// Calling the routine as,
//
//     impulse(A,B,C,D,NTIM,IU)
//
// plots the time response of the linear system to an impulse apllied
// to the single input IU.
//
// The routine can also be used to find the impulse response of a polynomial
// transfer function as,
//
//     impulse(NUM,DEN,NTIM)
//
// Version JBL 940519
//-----------------------------------------------------------------------------

require tfchk tf2ss abcdchk c2d simdata timeplot ltitr

impulse2 = function(a,b,c,d,ntim,iu)
{

   if ( (nargs <= 2) || (nargs > 6) ) {
       error("impulse: Error in number of input arguments.");
   }

   // Convert to S-S if in TF form
   if (nargs < 5) {
       Adum=tfchk(a,b);
       num=Adum.numc;
       den=Adum.denc;
       Adum=tf2ss(num,den);
       A=Adum.a;
       B=Adum.b;
       C=Adum.c;
       D=Adum.d;
       nu=1;
       ny=1;
       iiu=1;
   } else {
       // See if a,b,c, and d are compatible.
       A=a;
       B=b;
       C=c;
       D=d;
       msg=abcdchk(A,B,C,D);
       if (msg != "") {
           estr="IMPULSE: "+msg;
           error(estr);
       }
       ny=c.nr;
       nu=b.nc;
       if (!exist(iu)) {
           iiu=-1;
       } else {
           iiu=iu;
       }
   }

   // Determine sample frequency
   f=50*max(abs(eig(a).val)) / (2*pi);
   // Convert continuous model to discrete with 1/f as the sample time
   Adum=c2d(A,B,1/f);
   aa=Adum.phi;
   bb=Adum.gamma;
   cc=C;
   dd=zeros(ny,nu);

   // Set-up the time axis vector
   t=[0:(ntim-1)/f:1/f]';


   // Plot the simulation output

   if (iiu < 0) {
       for (i in 1:ny) {
            for (j in 1:nu) {
                 // Simulate the response by using ltitr
                 x=ltitr(aa,bb[;j],zeros(ntim,1),b[;j]);
                 y=x*C[i;].';
                 timeplot(i,j,f,ntim,y,"Impulse Response");
                 pause();
            }
       }
   } else {
       for (i in 1:ny) {
            x=ltitr(aa,bb[;iiu],zeros(ntim,1),b[;iiu]);
            y=x*C[i;].';
            timeplot(i,iiu,f,ntim,y,"Impulse Response");
            pause();
       }
   }

};

