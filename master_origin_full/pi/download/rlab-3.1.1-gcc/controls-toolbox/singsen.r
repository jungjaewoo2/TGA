//------------------------------------------------------------------------------
//
// singsen
//
// Syntax: D=singsen(A,dAdp)
//
// This routine computes the partial derivative of the singular values
// of a matrix with respect to a scalar parameter. All that is needed
// is the matrix and the partial of the matrix with respect to a parameter
// p ( d(A)/dp ).
//
// The routine returns a column vector containing the partial derivatives
// of the singular values of the matrix A which has a gradient of d(A)/dp.
// This routine computes the partial derivative of the singular values
//
//  Ref. Chapter 2 of Junkins and Kim, Dyn and Control of Flex. Structures
//
// Copyright (C), by Jeffrey B. Layton, 1994
// Version JBL 931120
//------------------------------------------------------------------------------

singsen = function(A,dAdp)
{
   if (nargs != 2) {
       error("SINGSEN: Wrong number of input arguments.");
   }

// Compute SVD of A
   Adum=svd(A,"A");
   U=Adum.u;
   S=diag(Adum.sigma);
   V=Adum.vt';

// Error Trap in SVD of A
   Acheck=U*S*Adum.vt;
   Aerr=norm(A-Acheck);
   if (Aerr > 1.0e-12) {
       error("SINGSEN: Error in taking svd of A");
   }

   n=max(size(S));
   svsen=zeros(n,1);

   for (i in 1:n) {
       svsen[i]=real(U[;i]'*dAdp*V[;i]);
   }

   return svsen;
};

