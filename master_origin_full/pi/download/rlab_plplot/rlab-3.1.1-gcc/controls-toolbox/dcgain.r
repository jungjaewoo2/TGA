//--------------------------------------------------------------------------
//
// dcgain
//
// Syntax: k=dcgain(a,b,c,d)
//
//
// This routine finds the D.C. gain (or low frequency gain) of a continuous
// system. Calling the routine as,
//
//     k=dcgain(a,b,c,d)
//
// computes the steady-state gain of the continuous state-space system.
//
// Calling the routine as,
//
//    k=dcgain(num,den)
//
// computes the steady-state gain of the continuous polynomial transfer
// function system,
//
//    G(s)= num(s)/den(s)
//
// where num contains the numerator polynomial coefficients and den
// contains the denominator polynomial coefficients in descending order.
//
// Copyright (C), by Jeffrey B. Layton, 1994
// Version JBL 940917
//--------------------------------------------------------------------------

require abcdchk  tfchk

dcgain = function(a,b,c,d)
{

// Transfer Function Representation
   if (nargs == 2) {
       A=tfchk(a,b);
       num=A.numc;
       den=A.denc;
       if ( (length(den) == 0) || (length(num) == 0) ) {
           k=[];
       else
           k=num[;length(den)]/den[length(den)];
       }

// State-Space Representation
   else if (nargs == 4) { 
       msg="";
       msg=abcdchk(a,b,c,d);
       if (msg != "") {
           estr="DCGAIN: "+msg;
           error(estr);
       }
       k=-c/a*b+d;
   else
       error("dcgain: Wrong Number of Input Arguments.");
   }}

   return k;
};

