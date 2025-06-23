//-------------------------------------------------------------------//

// Synopsis:	Pre-multiply by random orthogonal matrix.

// Syntax:	B = qmult ( A )

// Description:

//	B is Q*A where Q is a random real orthogonal matrix from the
//	Haar distribution, of dimension the number of rows in
//	A. Special case: if A is a scalar then QMULT(A) is the same as

//		qmult(eye(a))

//       Called by RANDSVD.

//       Reference:
//       G.W. Stewart, The efficient generation of random
//       orthogonal matrices with an application to condition estimators,
//       SIAM J. Numer. Anal., 17 (1980), 403-409.

//	This file is a translation of qmult.m from version 2.0 of
//	"The Test Matrix Toolbox for Matlab", described in Numerical
//	Analysis Report No. 237, December 1993, by N. J. Higham.

//-------------------------------------------------------------------//

require norm

qmult = function ( A )
{
  n = A.nr; m = A.nc;

  //  Handle scalar A.
  if (max(n,m) == 1)
  {
    n = A;
    A = eye(n,n);
  }

  d = zeros(n,n);

  for (k in n-1:1:-1)
  {
    // Generate random Householder transformation.
    rand("normal", 0, 1);
    x = rand(n-k+1,1);
    s = norm(x, "2");
    sgn = sign(x[1]) + (x[1]==0);	// Modification for sign(1)=1.
    s = sgn*s;
    d[k] = -sgn;
    x[1] = x[1] + s;
    beta = s*x[1];

    // Apply the transformation to A.
    y = x'*A[k:n;];
    A[k:n;] = A[k:n;] - x*(y/beta);
  }

  // Tidy up signs.
  for (i in 1:n-1)
  {
    A[i;] = d[i]*A[i;];
  }

  A[n;] = A[n;]*sign(rand(1,1));
  B = A;

  return B;
};
