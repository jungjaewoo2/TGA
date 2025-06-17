//-------------------------------------------------------------------//

//  Synopsis:   Create a sparse identity matrix.

//  Syntax:	speye ( S1, S2 )
//              speye ( A )

//  Description:

//  Create a sparse identity matrix with number of rows: S1, and
//  number of columns: S1. If the input A is a matrix, the first
//  element specifies the number of rows, and the second element
//  specifies the number of columns.

//  See Also: eye, ones, zeros
//-------------------------------------------------------------------//

speye = function ( m, n )
{
  #
  # Emulate Matlab's interface...
  #

  if (nargs == 1)
  {
    if (length (m) == 1)
    {
      n = m;
    else if (length (m) == 2) {
      n = m[2];
      m = m[1];
    else
      rerror ("speye: incorrect usage");
  }}}

  #
  # Create the matrix triplet form of the identity matrix.
  #

  s = min (m, n);
  k = (1:round (s))';
  K = [k, k, ones (s, 1)];
  if (m != n)
  {
    K = [K; m, n, 0];
  }

  #
  # Now convert it into compressed row-wise storage.
  #

  return spconvert (K);
};
