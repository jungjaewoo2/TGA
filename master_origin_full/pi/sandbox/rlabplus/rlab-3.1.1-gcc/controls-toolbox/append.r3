//-------------------------------------------------------------------
//
// append
//
// Syntax: G = append(a1,b1,c1,d1,a2,b2,c2,d2)
//         </aa;ba;ca;da/> = append(a1,b1,c1,d1,a2,b2,c2,d2)
//
// This routine appends the dynamics of two state-space systems
// together. It combines the state-space matrices of two systems
// defined by the follwing:
//          .
//         |x1| = |A1 0| |x1| + |B1 0| |u1|
//         |x2|   |0 A2| |x2| + |0 B2| |u2|
//
//         |y1| = |C1 0| |x1| + |D1 0| |u1|
//         |y2|   |0 C2| |x2| + |0 D2| |u2|
//
// where system "1" and system "2" are combined to formed the
// appended system.
//
// The results are returned in a list:
//
//    G.aa = appended A matrix
//    G.ba = appended B matrix
//    G.ca = appended C matrix
//    G.da = appended D matrix
//
// Copyright (C), by Jeffrey B. Layton, 1994
// Version JBL 940922
//-------------------------------------------------------------------

require abcdchk

append = function(a1,b1,c1,d1,a2,b2,c2,d2)
{
   
   if (nargs != 8) {
       error("APPEND: Error in number of input arguments.");
   }
   // Check input systems for compatibility
   msg=abcdchk(a1,b1,c1,d1);
   if (msg != "") {
       estr="System 1: "+msg;
       error(estr);
   }

   msg=abcdchk(a2,b2,c2,d2);
   if (msg != "") {
       estr="System 2: "+msg;
       error(estr);
   }

   // Applend Dynamics
   aa=[a1,zeros(a1.nr,a2.nc);zeros(a2.nr,a1.nc),a2];
   ba=[b1,zeros(a1.nr,d2.nc);zeros(a2.nr,d1.nc),b2];
   ca=[c1,zeros(d1.nr,a2.nc);zeros(d2.nr,a1.nc),c2];
   da=[d1,zeros(d1.nr,d2.nc);zeros(d2.nr,d1.nc),d2];

   return << aa=aa; ba=ba; ca=ca; da=da >>;
};

