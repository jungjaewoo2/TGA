//---------------------------------------------------------------------------
//
// gram
//
// Syntax: g=gram(a,b)
//
// This routine computes the Controllability and observability
// gramians. Calling the routine as A=gram(a,b) returns the
// controllability gramian:
//
//		Gc = integral {exp(at)bb'exp(a't)} dt
//
// Calling the routine as A=gram(a',c') returns the 
// observability gramian:
//
//		Go = integral {exp(ta')c'cexp(ta)} dt
//
// Ref: Laub, A., "Computation of Balancing Transformations", Proc. JACC
//      Vol.1, paper FA8-E, 1980.
// Ref: Skelton, R. "Dynamic Systems Control Linear System Analysis and
//      Synthesis," John Wiley and Sons, 1988.
//
// Copyright (C), by Jeffrey B. Layton, 1994
// Version JBL 940918
//---------------------------------------------------------------------------

require lyap

gram = function(a,b)
{

   if (nargs != 2)
   {
     error("GRAM: Wrong number of input arguments");
   }

   n=a.nr;
   A=svd(b,"A"); //  need a full decomposition
   u=A.u;

   s=[A.sigma;zeros(n,1)][1:n];

   g=u*lyap(u'*a*u,,s*s')*u';

   return g;
};


