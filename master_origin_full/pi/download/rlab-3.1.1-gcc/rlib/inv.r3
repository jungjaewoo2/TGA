//-------------------------------------------------------------------//

// Synopsis:	Compute the inverse of a matrix.

// Syntax:	inv ( P )

// Description:

//	inv ( A ) computes the inverse of A. The inverse function,
//	because it uses the solve function, uses a symmetric/hermitian
//	solver if the coefficient matrix is symmetric/hermitian.

//-------------------------------------------------------------------//

inv = function ( A )
{
  return solve (A, eye (size (A)));
};
