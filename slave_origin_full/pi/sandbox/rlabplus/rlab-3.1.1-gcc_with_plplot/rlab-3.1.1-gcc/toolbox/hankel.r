//-------------------------------------------------------------------//

//  Syntax:	hankel ( C )
//		hankel ( C , R )

//  Description:

//  Hankel(C) is a square Hankel matrix whose first column is C and
//  whose elements are zero below the first anti-diagonal.

//  Hankel(C,R) is a Hankel matrix whose first column is C and whose
//  last row is R.

//  Hankel matrices are symmetric, constant across the anti-diagonals,
//  and have elements H[i;j] = R[i+j-1].  

//  See also toeplitz.r.
//-------------------------------------------------------------------//

hankel = function ( c_ , r_ ) 
{
  local(c, r, nc, h, j, nr);

  if (class(c_) != "num") { error("hankel: Input must be numeric"); }

  c = c_[:];
  nc = length (c);

  if (!exist (r_)) 
  {
    h = zeros (nc, nc);
    for (j in 1:nc) 
    {
      h[1:nc-j+1;j] = c[j:nc];
    }
  else
    if (class (r_) != "num") { error("hankel: Inputs must be numeric"); }
    if (length(c_) != length(r_)) 
    {
      error("hankel: Inputs must have the same length");
    }
    r = r_[:];
    nr = length (r);
    h = zeros (nc, nr);
    if (c[nc] != r[1]) 
    {
      error("hankel: First element of row must match last element of column");
    }
    for (j in 1:nc-1) 
    {
      h[;j] = [c[j:nc-1]; r[1:j]];
    }
    h[;nc] = r;
  }
  return h;
};
