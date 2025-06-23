//-------------------------------------------------------------------------
//
// dgram
//
// Syntax: g=dgram(a,b)
//
// This routine computes the discrete controllability and
// observability gramians. Calling the routine as g=dgram(a,b)
// returns the discrete controllability gramian. Calling the
// routine as g=dgram(a',c') returns the discrete observability
// gramian.
//
// For detail regarding how the gramians are computed, see
// the routine gram (continuous version of dgram).
//
// Ref: Kailath, T. "Linear Systems", Prentice-Hall, 1980.
// Ref: Laub, A., "Computation of Balancing Transformations", Proc. JACC
//      Vol.1, paper FA8-E, 1980.
// Ref: Skelton, R. "Dynamic Systems Control Linear System Analysis and
//      Synthesis," John Wiley and Sons, 1988.
//
// Copyright (C), by Jeffrey B. Layton, 1994
// Version JBL 940918
//-------------------------------------------------------------------------

require lyap

dgram = function(a,b)
{

   if (nargs != 2) {
       error("DGRAM: Wrong number of input arguments.");
   }

// Begin
   A=svd(b);
   u=A.u;
   s=A.sigma;

   g=u*lyap(u'*a*u,,s*s')*u';

   return g;
};

