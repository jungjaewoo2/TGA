//---------------------------------------------------------------------------
//
// isstr
//
// Syntax: i=isstr(a)
//
// This routine checks if the input argument a is a string. It
// returns a "1" if it is a string and a 0 otherwise.
//
// Copyright (C), by Jeffrey B. Layton, 1994
// Version JBL 940918
//---------------------------------------------------------------------------

isstr = function(a)
{

   if (class(a) == "string") {
       i=1;
   else
       i=0;
   }

   return i;
};

