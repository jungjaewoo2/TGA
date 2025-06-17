//--------------------------------------------------------------------------------
//
// housh
//
// Syntax: A=housh(u,j,eps)
//
// This routine constructs a Householder Transformation H=I-s*UU'
// that "mirrors" a vector u to the Jth unit vector. If Norm(U) < eps
// then zero = 1 (True) else zero = 0 (False).
//
// Used in tzero and tzreduce.
//
// REFERENCE: Ported to RLaB from the FORTRAN Code contained in:
// Emami-Naeini, A., and Van Dooren, A., "Computation of Zeros of Linear
// Multivaraible Systems," Automatica, Vol. 18, No. 4, April 1982, pp. 415-430.
//
// The routine returns a list.
//
// A.u     - Transformation Vector
// A.s     - S Vector
// A.zero  - Zero Flag
//
// Copyright (C), by Jeffrey B. Layton, 1994
// Version JBL 931018
//--------------------------------------------------------------------------------

housh = function(u,j,eps)
{
   local(s,alfa,dum1,u,s,zero)

   s=sum(u.*u);
   alfa=sqrt(s);
   if (alfa <= eps) {
       zero=1;
       return << u=u[:]; s=s; zero=zero >>
   }

   zero=0;
   dum1=u[j];
   if (dum1 > 0) {
       alfa=-alfa;
   }
   u[j]=dum1-alfa;
   s=1 ./(s-alfa*dum1);

   return << u=u[:]; s=s; zero=zero >>;
};

