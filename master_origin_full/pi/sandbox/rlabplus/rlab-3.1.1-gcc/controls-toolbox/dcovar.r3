//-----------------------------------------------------------------------
//
// dcovar
//
// Syntax: R = dcovar(A,B,C,D,W)  
//         R = dcovar(NUM,DEN,W)
//
// This routine computes the covariance response of a discrete
// system to white noise input. Calling the routine as,
//
//    R=dcovar(A,B,C,D,W)
//
// computes the covariance response of the discrete state-space
// system defined by A,B,C,D to Gaussian white noise input with
// intensity W given by the following,
//  
//              E[w(k)w(tau)'] = W delta(k-tau)
//  
// where delta(k,tau) is the kronecker delta and E is the expectation
// operator.
//
// The routine returns both the state covariance matrix X and the
// output covariance matrix Y given by the following,
//
//     X=E[xx']
//     Y=E[yy']
//
// where x is the state vector and y is the output vector.
//
// This routine can also be used to compute the covariance response
// of a transfer function type system by calling it as,
//
//     R=covar(NUM,DEN,W)
//
// where NUM is a vector containing the coefficients of the numerator
// transfer function polynomial, DEN contains the coefficients of the
// denominator transfer function polynomial, and W is the intensity
// of the white noise input.
//
// The results are returned in a list. For example,  R=covar(A,B,C,D,W)
// produces,
//
//  R.Y = output covariance matrix (Square of Output RMS)
//  R.X = state covariance matrix
//
// Ref: Skelton, R. "Dynamic Systems Control Linear System Analysis and
//      Synthesis," John Wiley and Sons, 1988.
//
// Copyright (C), by Jeffrey B. Layton, 1994
// Version JBL 940918
//-----------------------------------------------------------------------

require abcdchk tfchk tf2ss dlyap 

dcovar = function(a,b,c,d,w)
{

   global (eps,pi)
   
   // Error trap for wrong number of inputs.
   if ( nargs != 3 && nargs != 5) {
      error("DCOVAR: Wrong number of input arguments.");
   }

   // Error Traps for matrix compatibility
   #if (a.nr != b.nr) {
   #    error("DCOVAR: A and B have different number of rows.");
   #}

   #if (w.nr != w.nc) {
   #    error("DCOVAR: W needs to be square.");
   #}

   #if (w.nr != b.nr) {
   #    error("DCOVAR: rows(W) not equal to rows(B).");
   #}

   #if (w.nc != b.nr) {
   #    error("DCOVAR: cols(W) not equal to rows(B).");
   #}

   if (nargs == 3) {   
       // Transfer Function Representation
       Dum=tfchk(a,b);
       num=Dum.numc;
       den=Dum.denc;
       W=c;
       // Convert to state-space form (using MATLAB compatible conversion)
       Dum=tf2ss(num,den,2);
       A=Dum.a;
       B=Dum.b;
       C=Dum.c;
       D=Dum.d;
   } else { if (nargs == 5) {
       // State-Space Representation
       A=a;
       B=b;
       C=c;
       D=d;
       W=w;
       msg=abcdchk(A,B,C,D);
       if (msg != "") {
           estr="DCOVAR: "+msg;
           error(estr);
       }
   } else {
       error("DCOVAR: Wrong number of input arguments");
       }
   }


   // Find state covariance by solving discrete lyapunov
   // Then find the output covariance matrix Y.
   X=dlyap(A,B*W*B');
   Y=C*X*C';
   
   
   // Systems with non-zero D matrix have infinite output covarianc
   if (any(abs(D[:])>eps)) {
      printf("Warning: Systems with non-zero D matrix have infinite output covariance\n");
      pinf = inf()*ones(length(Y));
      if (!isempty(Y)) {
        Y[(abs(D*W*D')>eps)] = pinf[abs(D*W*D')>eps]; 
      }
   }

   // A valid covariance must be positive semi-definite.
   if (min(real(eig(X).val)) < -eps) {
      printf("Warning: Invalid covariance - not positive semi-definite. Returning infinity.\n");
      X=inf()*ones(X.nr,X.nc);
      Y=inf()*ones(Y.nr,Y.nc);
      return << X=X; Y=Y >>
   }

   return << X=X; Y=Y >>
};

