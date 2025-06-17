//-------------------------------------------------------------------//

//  Synopsis:	Compute the matrix or vector norm.

//  Syntax:	norm ( A )
//              norm ( A , TYPE )
//              norm ( V , P )

//  Description:

//  This function is a cover-function for the two builtin functions
//  mnorm(), which computes matrix-norms, and vpnorm(), which computes
//  vector P-norms.

//  If A is a matrix, valid TYPE arguments are:

//       "M" or "m"
//         returns max(abs( MATRIX ))

//       "1", "O" or "o"
//         return the 1-norm (default)

//       "2" returns the matrix 2-norm (largest singular value)

//       "I" or "i"
//         returns the infinity-norm

//       "F", "f", "E" or "e"
//         returns the Frobenius norm.

//       inf()
//         returns the infinity-norm.

//  If A is a vector (either row or column) valid TYPE arguments
//  are numeric values of P for computing the vector P-norm. A sample
//  Rlab function which computes the vector P-norm is provided below,
//  although the actual implementation is via a builtin function.

//    pnorm = function ( V , P )
//    {
//       return (sum ( V.^P )).^(1/P);
//    }

//  See Also: mnorm, vpnorm
//-------------------------------------------------------------------//

norm = function ( A , ntype )
{

  #
  # This check doesn't fit neatly anywhere else,
  # just get it out of the way...
  #

  if (exist (ntype))
  {
    if (class(ntype) == "string")
    {
      # This has to be done with the matrix norm function.
      return mnorm ( A , ntype );
    }
  }

  #
  # Now, Check for matrix input.
  #

  if (min (size (A)) != 1)
  {
    #
    # Two possibilties...
    #

    if (!exist (ntype)) {
      return mnorm (A);
    }

    if (class (ntype) == "string" || class (ntype) == "num") {
      return mnorm (A, ntype);
    }

  else

    #
    # Vector input...
    #

    if (!exist (ntype)) {
      error ("norm: second argument required");
    }

    if (class (ntype) != "num") {
      error ("norm: second argument must be numeric");
    }

    return vpnorm (A, ntype);
  }
};
