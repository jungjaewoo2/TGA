//-------------------------------------------------------------------//

// Synopsis:    Toeplitz matrix.

// Syntax:	toeplitz ( C )
//		toeplitz ( C , R )

//  Description:

//      The toeplitz function returns a non-symmetric Toeplitz matrix
//      having C as its first column and R as its first
//      row. `toeplitz(C)' is a symmetric (or Hermitian) Toeplitz 
//      matrix. 

//      If the 1st element of C and R do not agree, then C[1] "wins".

//  See also hankel.r
//-------------------------------------------------------------------//

toeplitz = function ( c , r ) 
{
  local(c, r, nc, nr, t, i);

  if (class(c) != "num") 
  {
    error ("toeplitz: Inputs must be numeric");
  }

  c = c[:];                 // Force column-vector.
  nr = length (c);          // No. of rows in t.

  if (!exist (r)) 
  {
    r = c';                 // Symmetric / Hermitian
  else
    if (class(r) != "num") 
    { 
      error ("toeplitz: Inputs must be numeric"); 
    }
  }
  
  r = r[:].';               // Force row-vector.
  nc = length (r);          // No. of columns in t.
  t = zeros (nr, nc);

  for (i in 1:min(nr,nc))
  {
    t[i;i:t.nc] = r[1:r.n-i+1];  // Fill upper triangle rows.
    t[i:t.nr;i] = c[1:c.n-i+1];  // Fill lower triangle columns.
  }
  return t;
};
