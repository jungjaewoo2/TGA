//-------------------------------------------------------------------------------
// c2dm
//
// Syntax: G=c2dm(A,B,C,D,T,method,om1)
//
// This routine converts a continuous linear time invariant (LTI)
// system to a discrete time system. This routine allows the user
// to choose from several conversions methods.
//
//  "zoh"     Convert to discrete time assuming a zero order hold
//            on the inputs.
//  "foh"     Convert to discrete time assuming a first order hold
//            on the inputs.
//  "tustin"  Convert to discrete time using a bilinear (Tustin)
//            approximation to the derivative.
//  "prewarp" Convert to discrete time using a bilinear (Tustin)
//            approximation, but with frequency pre-warping.
//            The desired frequency is om1.
//
// The routine may be called as,
//
//      G=c2dm(A,B,C,D,T,"method",om1)
//
// where the state-space model is given by (A,B,C,D), the sample time
// is given by T, the desired conversion is given by "method," and the
// desired frequency for prewarping is given by om1.
//
// The results are passed back in a list:
//
//   G.ad = discrete A matrix
//   G.bd = discrete B matrix
//   G.cd = discrete C matrix
//   G.dd = discrete D matrix
//
// The routine can also be used to convert systems given in the Transfer
// Function form:
//
//     G=c2dm(num,den,T,"method")
//
// which converts a continuous polynomial G(s)=num(s)/den(s) to a discrete
// time polynomial G(z)=num(z)/den(z).
//
// Ths results are passed back in a list:
//
//   G.ad = discrete numerator
//   G.bd = discrete denominator
//
// Ref:  Frankin, G. F., Powell, J. D., Workman, M. L., "Digital Control of
//       Dynamic Systems," Addison-Wesley, Reading Mass., 1990.
//
// Copyright (C), by Jeffrey B. Layton, 1994
// Version JBL 940108
//-------------------------------------------------------------------------------

require tfchk tf2ss abcdchk c2d ss2tf isstr

static (c2dm_compute)    // Hide this function

c2dm = function(A,B,C,D,T,method,om1)
{
   // switch on which type of system is being input
   if (nargs == 2) {
       // Transfer Function with no method specified
       P=tfchk(A,B);
       P=tf2ss(P.numc,P.denc,2);
       tf=1;
       return c2dm_compute (P.a, P.b, P.c, P.d, C, "zoh", om1, tf);
   else if (nargs == 4) {
       // Transfer Function with method specified
       P=tfchk(A,B);
       P=tf2ss(P.numc,P.denc,2);
       tf=1;
       return c2dm_compute (P.a, P.b, P.c, P.d, C, D, om1, tf);
   else if (nargs == 5) {
       if (isstr(D)) {
           // Transfer function with method specified and prewarp freq.
           P=tfchk(A,B);
           P=tf2ss(P.numc,P.denc,2);
           tf=1;
           return c2dm_compute (P.a, P.b, P.c, P.d, C, D, T, tf);
       else
           // State-Space model without method specified
           msg="";
           msg=abcdchk(A,B,C,D);
           if (msg != "") {
               estr="C2DM: "+msg;
               error(estr);
           }
           return c2dm_compute (A, B, C, D, T, "zoh", om1, tf);
       }
   else
       // State Space system with specified method
       msg="";
       msg=abcdchk(A,B,C,D);
       if (msg != "") {
           estr="C2DM: "+msg;
           error(estr);
       }
       return c2dm_compute (A, B, C, D, T, method, om1, tf);
   }}}

};

//-------------------------------------------------------------------------------
//
// This is where the computation is performed. Note that c2dm_compute is a
// static variable and is never seen from the global workspace.
//

c2dm_compute = function(A, B, C, D, T, method, om1, tf)
{

   // Switch on Method

   if (strsplt(method)[1] == "z") {
       // Zero Order Hold
       // ---------------
       // Just call c2d which already performs the zero order hold conversion
       P=c2d(A,B,T);
       ad=P.phi;
       bd=P.gamma;
       cd=C;
       dd=D;
   else if (strsplt(method)[1] == "f") {
       // First Order Hold Using Triangle Approximation
       // ---------------------------------------------
       // Form FT matrix
       FT = [      A,                B,        zeros(A.nr,B.nc);
             zeros(B.nc,A.nr),zeros(B.nc,B.nc),eye(B.nc,B.nc)/T;
             zeros(B.nc,A.nr),zeros(B.nc,B.nc),zeros(B.nc,B.nc)];
       // Take Matrix exponential
       E=expm(FT*T);
       // Find Phi, Gamma1, Gamma2 partitions
       ad=E[1:A.nr;1:A.nr];
       gam1=E[1:A.nr;(A.nr+1):(A.nr+B.nc)];
       gam2=E[1:A.nr;(A.nr+B.nc+1):(A.nr+(2*B.nc))];
       // Form B,C,D discrete matrices
       bd=gam1+ad*gam2-gam2;
       cd=C;
       dd=D+C*gam2;
   else if (strsplt(method)[1] == "t") {
       // Bilinear (Tustin)
       // -----------------
       I=eye(A.nr,A.nr);
       E=inv(I-A.*T/2);
       ad=(I+A.*T/2)*E;
       bd=E*B;
       cd=T*C*E;
       dd=cd*B/2+D;
   else if (strsplt(method)[1] == "p") {
       // Bilinear (Tustin) with frequency prewarping
       // -------------------------------------------
       if ( nargs != 5) {
           if ( nargs != 7) {
               Tw=2*tan(om1*T/2)/om1;
               I=eye(A.nr,A.nr);
               E=inv(I-A.*Tw/2);
               ad=(I+A.*Tw/2)*E;
               bd=E*B;
               cd=Tw*C*E;
               dd=cd*B/2+D;
           else
               estr="c2dm: The critical frequency must be specified ";
               estr=estr+"when using the prewarp method.";
               error(estr);
           }
       }
   else
       estr="c2dm: The critical frequency must be specified ";
       estr=estr+"when using the prewarp method.";
       error(estr);
   }}}}

   if (tf) {
       Adum=ss2tf(ad,bd,cd,dd,1);
       ad=Adum.num;
       bd=Adum.den;
       return << ad=ad; bd=bd >>;
   }

   return << ad=ad; bd=bd; cd=cd; dd=dd >>;
};

