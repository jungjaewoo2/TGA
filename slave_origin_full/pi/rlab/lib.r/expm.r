//-------------------------------------------------------------------//

//  Syntax:	expm ( A )

//  Description:

//  The expm function computes the matrix exponential, or exp (A),
//  where A is a square matrix.

//  Reference: "Matrix Computations" Gene  Golub, and Charles Van
//  		Loan. Section 11.3.

//-------------------------------------------------------------------//

expm = function ( a_ )
{
  local (X, a, c, cX, d, k, N, n, q, j);

  if (a_.nr != a_.nc) { error ("expm: A must be square"); }

  //
  // Scale A if neccesary, such that the infinity-norm of A divided
  // by some power of 2 is less than 1/2.
  //

  j = norm(a_, "i");
  if (j > 0) 
  { 
    j = max ([0, fix (log (j)/log (2))+2]);
  }
  a = a_/2^j;
  n = a.nr;

  //
  // k = 1
  //

  c = 1/2;
  X = a; 
  d = eye (n, n) - c*a;
  N = eye (n, n) + c*a;
  q = 6;

  for (k in 2:q)
  {
    c = c*(q - k + 1)/(k*(2*q - k + 1));
    X = a*X;
    cX = c*X;
    N = N + cX;
    d = d + (-1)^k*cX;
  }

  //
  // Solve d*F = N for F, overwrite N for efficiency
  //

  N = d\N;

  for (k in 1:j)
  {
    N = N*N;
  }

  return N;
};
