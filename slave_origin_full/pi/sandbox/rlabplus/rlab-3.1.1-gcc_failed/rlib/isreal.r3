//-------------------------------------------------------------------//

//  Synopsis:   Test if a matrix is ALL real.

//  Syntax:     isreal ( A )

//  Description:

//      isreal tests its argument, A, and returns TRUE (1) if A's
//      elements are all real, or FALSE (0) if any of A's elements are
//      complex. 

// See Also: any

//-------------------------------------------------------------------//

isreal = function ( A )
{
  return all (all (imag (A) == 0));
};
