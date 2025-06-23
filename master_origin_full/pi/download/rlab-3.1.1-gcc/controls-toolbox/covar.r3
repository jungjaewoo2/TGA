//-----------------------------------------------------------------------
//
// covar
//
// Syntax: R=covar(A,B,C,D,W)  or  R=covar(NUM,DEN,W)
//
// This routine computes the covariance response of a continuous
// system to white noise input. Calling the routine as,
//
//    R=covar(A,B,C,D,W)
//
// computes the covariance response of the continuous state-space
// system defined by A,B,C,D to Gaussian white noise input with
// intensity W given by the following,
//  
//    E[w(t)w(tau)'] = W delta(t-tau) 
//  
// where delta(t) is the dirac delta function and E is the expectation
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
//  ==================================================================
//  Warning: Unstable systems or systems with a non-zero D matrix have
//  an infinite covariance response.
//  ==================================================================
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
// Version JBL 940917
//-----------------------------------------------------------------------

require tfchk abcdchk tf2ss lyap 

covar = function(a,b,c,d,w)
{

   if ( nargs != 3 && nargs != 5) {
      error("covar: Wrong number of input arguments.");
   }

   if (nargs == 3) {   
       // T.F. syntax
       Dum=tfchk(a,b);
       num=Dum.numc;
       den=Dum.denc;
       W=c;
       // Convert to state-space form
       Dum=tf2ss(num,den,2);
       A=Dum.a;
       B=Dum.b;
       C=Dum.c;
       D=Dum.d;
   } else {
       A=a;
       B=b;
       C=c;
       D=d;
       W=w;
   }

   if (nargs == 5) {
       msg="";
       msg=abcdchk(A,B,C,D);
       if (msg != "") {
           estr="COVAR: "+msg;
           error(estr);
       }
   }

   X=lyap(A,,B*W*B');
   Y=C*X*C';

   // Systems with non-zero D matrix have infinite output covariance.
   eps = epsilon();
   if ( any( any( abs(d) > eps) ) ) {
       printf("Warning: Systems with non-zero D matrix have infinite output covariance.\n");
       pinf=inf()*ones(Y.nr,Y.nc);
       if (!isempty(Y)) {
           Y=inf()*ones(Y.nr,Y.nc);
       }
   }

   // A valid covariance must be positive semi-definite.
   if (min(real(eig(X).val)) < -eps) {
       printf("Warning: Invalid covaraince - not positive semi-definite. Returning infinity.\n");
       X=inf()*ones(X.nr,X.nc);
       Y=inf()*ones(Y.nr,Y.nc);
       return << X=X; Y=Y >>;
   }

   // The system must be stable for a valid covariance
   if (any(real(eig(a).val) > 0) ) {
       printf("Warning: Unstable system. Returning Infinity.\n");
       X=inf()*ones(length(X));
       Y=inf()*ones(length(Y));
       return << X=X; Y=Y >>;
   }

   return << X=X; Y=Y >>;
};

