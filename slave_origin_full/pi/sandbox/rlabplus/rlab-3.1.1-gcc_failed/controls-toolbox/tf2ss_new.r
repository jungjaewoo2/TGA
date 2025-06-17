//----------------------------------------------------------------------------
//
// tf2ss
//
// Syntax: A=tf2ss(num,den,meth)
//
// This routine calculates the state-space representation:
// 	.
// 	x = Ax + Bu
// 	y = Cx + Du
//
// from the system:
//
// 	        NUM(s) 
// 	H(s) = --------
// 	        DEN(s)
//
// for a single input.  Vector den must contain the coefficients of
// the denominator in descending powers of s.  Matrix num must 
// contain the numerator coefficients with as many rows as there are
// outputs y.
//
// There are 3 choices to perform the conversion:
//
//     meth = 1: Returns A,B,C,D in Phase-Variable Form (default)
//     meth = 2: Returns A,B,C,D in Controller Canonical Form
//               (MATLAB uses this method)
//     meth = 3: Returns A,B,C,D in Rectangular Form
//
// This calculation also works for discrete systems. To avoid confusion
// when using this function with discrete systems, always use a numerator
// polynomial that has been padded with zeros to make it the same length
// as the denominator.
//
// The state-space matrices are returned in a list. For example:
//
//      A=tf2ss(num,den,meth);
//
//      A.a = a matrix
//      A.b = b matrix
//      A.c = c matrix
//      A.d = d matrix
//
// Ref: (1) Skelton, R. "Dynamic Systems Control Linear System Analysis and
//          Synthesis," John Wiley and Sons, 1988.
//      (2) Kailath, "Linear Systems," Prentice-Hall, 1980
//
// Copyright (C), by Jeffrey B. Layton, 1994
// Version JBL 940106
//----------------------------------------------------------------------------

tf2ss = function(num,den,meth)
{

   mnum=num.nr;
   nnum=num.nc;
   nden=den.nc;

   // Check for meth
   if (! exist(meth)) {meth = 1;}

   // Define local version of den and num
   denl=den;
   numl=num;

   // Check for null systems
   if ((nden == 0) && (nnum == 0)) {
        a=[];
        b=[];
        c=[];
        d=[];
        return << a=a; b=b; c=c; d=d >>;
   }

   // Strip leading zeros from denominator
   inz=find(den != 0);
   denl=denl[inz[1]:nden];
   nden=den.nc;

   // Check for proper numerator
   if (nnum > nden) {
       // Try to strip leading zeros to make proper
       if (all(all(numl[;1:(nnum-1)] == 0))) {
           numl=numl[;(nnum-nden+1):nnum];
           mnum=numl.nr;
           nnum=numl.nc;
       else
           error("TF2SS: Denominator must be higher or equal order than numerator.");
       }
   }

   // Pad numerator with leading zeros, to make it have the same number of
   // columns as the denominator, and normalize it to den[1]
   numl=[zeros(mnum,nden-nnum),numl]./denl[1];
   nnum=numl.nc;

   //
   // ===================
   // Phase-Variable Form
   // ===================
   //
   if (meth == 1) {
       // Now do the rest, starting by normalizing den to den[1],
       denl=denl[2:nden]./denl[1];

       // Define constant term for Phase-Variable Form
       con=numl[;1];

       // Switch Order of Coefficients in num and den
       denl=denl[nden-1:1:-1];
       numl=numl[numl.nc:1:-1];

       // Create the State-Space Matrices
       a=[zeros(nden-2,1), eye(nden-2,nden-2)];
       a=[a;-denl];
       b=[zeros(nden-2,1);1];
       c=numl[1:nnum-1]-con*denl;
       d=con;

       return << a=a; b=b; c=c; d=d >>;

   //
   // =============================================
   // Controller Canonical Form (MATLAB compatible)
   // =============================================
   //
   else if (meth == 2) {

       // Do the D-matrix first
       if (length(numl)) {
           d=numl[;1];
       else
           d=[];
       }

       // Handle special constant case:
       if (nden == 1) {
          a=[];
          b=[];
          c=[];
          return << a=a; b=b; c=c; d=d >>;
       }

       // Now do the rest, starting by normalizing den to den[1],
       denl=denl[2:nden]./denl[1];
       a=[-denl;eye(nden-2,nden-1)];
       b=eye(nden-1,1);

       if (mnum > 0) {
           c=numl[;2:nden]-numl[;1]*denl;
       else
           c=[];
       }

       return << a=a; b=b; c=c; d=d >>;

   //
   // ================
   // Rectangular Form
   // ================
   //
   else if (meth == 3) {

       // First define the constant con[1]
       con=numl[1];

       numl=numl[2:nnum]';
       denl=denl[2:nden]';

       // Compute a,b,c,d
       b=numl-con*denl;
       c=zeros(1,nnum-1);
       c[1]=1;
       d=con;
       if (nnum == 2) {
           a=[-numl];
       else
           a=[eye(nnum-2,nnum-2);
              zeros(1,nnum-2)];
           a=[-numl,a];
       }

       return << a=a; b=b; c=c; d=d >>;
    }}}
};

