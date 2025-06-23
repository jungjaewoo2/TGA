//-----------------------------------------------------------
//
// ctrb
//
// Syntax: co=ctrb(a,b)
//
// This routine forms the controllability matrix.
// co = ctrb(a,b) returns the controllability matrix
// co = [B AB A^2B ...]
//
// Copyright (C), by Jeffrey B. Layton, 1994
// Version JBL 940920
//-----------------------------------------------------------

ctrb = function(a,b)
{
// Start by setting co to b
   co = b;
   D=eye(a.nr,a.nc);
   for (i in 1:a.nc-1) {
        D=D*a;
        co=[co,D*b];
   }

   return co;
};

