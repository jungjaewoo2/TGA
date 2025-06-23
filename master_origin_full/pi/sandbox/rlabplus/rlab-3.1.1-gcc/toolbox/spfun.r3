//-------------------------------------------------------------------//

//  Synopsis:   Apply an arbitrary scalar function to each non-zero
//              element of a sparse matrix.

//  Syntax:     spfun ( A , FUN )

//  Description:

//  Apply an arbitrary function to each non-zero element of a sparse
//  matrix. The function may be either a builtin, or a user function,
//  it makes no difference.

//  FUN must be able to operate on vectors (Nx1).

//-------------------------------------------------------------------//

spfun = function ( A , fun )
{
  if (A.storage != "sparse")
  {
    error ("spfun: A must be sparse");
  }

  tmp = spconvert (A);
  tmp = [tmp[;1,2], fun (tmp[;3])];

  return spconvert ([tmp; A.nr, A.nc, 0]);
};
