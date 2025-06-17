//-------------------------------------------------------------------------
//
// augstate
//
// Syntax: G = augstate(a,b,c,d)
//         </aa;ba;ca;da/> = augstate(a,b,c,d)
//
// This routine augments the outputs of a state-space system with the
// states of the system. Calling the routine as G=augstate(A,B,C,D)
// produces the following augmented system:
// 	         .
// 	         x = Ax + Bu
// 
// 	        |y| = |C| x + |D| u
// 	        |x|   |I|     |0|
// 
// The results are returned in a list:
//
//     G.aa = augmented A matrix (no change from input)
//     G.ba = augmented B matrix (no change from input)
//     G.ca = augmented C matrix
//     G.da = augmented D matrix
//
// Copyright (C), by Jeffrey B. Layton
// Version JBL 940922
//-------------------------------------------------------------------------

rfile abcdchk

augstate = function(a,b,c,d)
{
   if (nargs != 4) {
       error("AUGSTATE: Wrong number of input arguments.");
   }

   // Check input system
   msg=abcdchk(a,b,c,d);
   if (msg != "") {
       estr="AUGSTATE: "+msg;
       error(estr);
   }

   // Augment
   aa=a;
   ba=b;
   ca=[c;eye(a.nr,a.nr)];
   da=[d;zeros(a.nr,b.nc)];

   return << aa=aa; ba=ba; ca=ca; da=da >>;
};

