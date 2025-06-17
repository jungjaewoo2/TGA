//---------------------------------------------------------------
//
// jpinv
//
// Syntax: Ap = jpinv(A)
//
// This routine computes the pseudo-inverse using Greville's
// Iterative Method.
//
// A is the matrix to be inverted and the results is passed
// back as Ap.
//
// Copyright (C), by Jeffrey B. Layton, 1994
// Version JBL 940918
//---------------------------------------------------------------

jpinv = function(A)
{
   if (nargs != 1) {
       error("JPINV: Wrong number of input arguments.");
   }

   Thresh = max(size(A))*epsilon()*norm(A,"2");
   con = A[;1]'*A[;1];
   conn = norm(con,"2");
   Ap = zeros(A.nr,A.nc);

   if (conn >= Thresh) {
      Ap=A[;1]'/[A[;1]'*A[;1]];
   else
      b=size(A[;1]');
      mm=b[1];
      nn=b[2];
      Ap=zeros(mm,nn);
   }

   if (A.nc != 1) {
      for (k in 2:A.nc) {
          dk=Ap*A[;k];
          ck=A[;k]-A[;1:(k-1)]*dk;
          con=ck'*ck;
          conn=norm(con,"2");
          if (conn >= Thresh) {
             bkT=ck'/(ck'*ck);
          else
             bkT=dk'*Ap/(1+dk'*dk);
          }
          Ap=[Ap-dk*bkT;bkT];
      }
   }

   return Ap;
};

