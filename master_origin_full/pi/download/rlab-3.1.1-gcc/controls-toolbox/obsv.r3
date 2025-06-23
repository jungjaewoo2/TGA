//--------------------------------------------------------------
//
// obsv
//
// Syntax: ob=obsv(a,c)
//
// This routines forms the observability matrix.
// ob=obsv(a,c,) returns the observability matrix.
// ob = [c; ca; ca^2; ...]
//
// Copyright (C), by Jeffrey B. Layton, 1994
// Version JBL 940919
//--------------------------------------------------------------

obsv = function(a,c)
{
   if (nargs != 2) {
       error("OBSV: Wrong number of input arguments.");
   }

   ob=c;
   D=eye(a.nr,a.nc);
   for (i in 1:a.nc-1) {
        D=D*a;
        ob=[ob;c*D];
   }

   return ob;
};

