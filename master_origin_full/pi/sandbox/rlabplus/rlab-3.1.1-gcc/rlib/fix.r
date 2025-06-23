//-------------------------------------------------------------------//

//  Synopsis:   Round array elements towards zero.

//  Syntax:     fix ( X )

//  Description:

//  Fix rounds the elements of X to the nearest integer towards zero.

//  See Also ceil, floor, round, sign, abs

//-------------------------------------------------------------------//

fix = function ( X )
{
  if (X.type == "real")
  {
    return sign (X) .* floor (abs (X));
  else if (X.type == "complex") {
    return (sign (real (X)) .* floor (abs (real (X)))) + ...
           ((0+1j)*(sign (imag (X))) .* floor (abs (imag (X))));
  else
    error("fix() only takes real or complex (scalar and matrix) arguments");
  }}
};
