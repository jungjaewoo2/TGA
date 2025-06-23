//--------------------------------------------------------------------
//
// ltifr
//
// syntax: g=ltifr(A,B,S)
//
// This routine computes the Linear time-invariant frequency response
// kernel. Calling the routine as g=ltifr(A,B,S) computes the frequency
// response of the following system:
//
//    g(s) = (sI-A)\b
//
// for the complex frequencies contained in the vector S. The column
// vector B must have as many rows as the matrix A. The results g is
// returned with size(A) rows and length(S) columns.
//
// Copyright (C), by Jeffrey B. Layton, 1994
// Version JBL 931011
//--------------------------------------------------------------------
            
ltifr = function(a,b,s)
{
   local(ns,na,e,g)

   ns=length(s);
   na=length(a);
   e=eye(na,na);
   for (i in 1:ns) {
        g[;i] = (s[i]*e-a)\b;
   }

   return g;
};

