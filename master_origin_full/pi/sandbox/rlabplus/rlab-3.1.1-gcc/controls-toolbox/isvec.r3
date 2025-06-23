//-----------------------------------------------------------------------
//
// isvec
//
// Syntax: i=isvec(a)
//
// This routine checks to see if the input a is a vector
// (column or matrix). If it is a vector, then a 1 is returned, otherwise a
// 0 is returned.
// (column or matrix).
//
// Copyright (C), by Jeffrey B. Layton, 1994
// Version JBL 940918
//-----------------------------------------------------------------------

isvec = function(a)
{

   if (nargs == 1) {
       return ( ((a.nr == 1) && (a.nc > 1)) || ((a.nc == 1) && (a.nr > 1)) );
   } else {
       error("ISVEC: Wrong number of input arguments.");
   }

};

