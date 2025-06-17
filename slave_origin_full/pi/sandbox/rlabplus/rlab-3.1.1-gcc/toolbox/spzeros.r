//-------------------------------------------------------------------//

//  Synopsis:   Create a sparse matrix with all zero elements.

//  Syntax:	spzeros ( M , N )
//              spzeros ( [ M, N ] )

//  Description:

//  Generate a sparse matrix with all zero elements.

//-------------------------------------------------------------------//

spzeros = function ( m, n )
{
  if (nargs == 1)
  {
    if (m.n >= 2)
    {
      n = m[2];
      m = m[1];
    else
      error ("spzeros: argument must be 2 element matrix");
    }
    spz = [m, n, 0];

  else if (nargs == 2) {

    spz = [m, n, 0];
  } }

  return spconvert (spz);
};
