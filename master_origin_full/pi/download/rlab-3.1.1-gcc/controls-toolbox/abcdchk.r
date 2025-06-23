//--------------------------------------------------------------------------
//
// abcdchk
//
// Syntax: msg=abcdchk(a,b,c,d)
//
//      This routine checks the matrices A,B,C,D for consistency of
//      dimensions. If everything is consistent, then the string msg
//      is returned as a null string. Otherwise, an error message
//      is returned in msg.
//
//      Valid systems with empty matrices are allowed.  
//
// Copyright (C), by Jeffrey B. Layton, 1994
// Version JBL 931016
//--------------------------------------------------------------------------

abcdchk = function(a,b,c,d)
{
   msg="";

   if (a.nr != a.nc) {
       msg="The A matrix must be square";
   }

   if (nargs > 1) {
       if (b.nc) {
           if (a.nr != b.nr) {
               msg="The A and B matrices must have the same number of rows.";
           }
       }
   }
   if (nargs > 2) {
       if (c.nr) {
           if (c.nc != a.nc) {
               msg="The A and C matrices must have the same number of columns.";
           }
       }
   }
   if (nargs > 3) {
   // Check if a,b,c matrices have zero dimensions. If so, just return.
       if ((a.nr+b.nr+c.nr) == 0) {
            return msg;
       }
   // Check C and D matrix compatibilities
       if (d.nc || b.nc) {
           if (d.nr != c.nr) {
               msg="The C and D matrices must have the same number of rows.";
           }
       }
   // Check B and D matrix compatibilities
       if (d.nr || c.nr) {
           if (d.nc != b.nc) {
               msg="The B and D matrices must have the same number of columns.";
           }
       }
    }

   return msg;
};

