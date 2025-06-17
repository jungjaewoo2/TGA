//-------------------------------------------------------------------//

//  Synopsis:   Return the lower triangular part of A.

//  Syntax:	tril ( A )
//		tril ( A , K )

//  Description:

//  tril(x) returns the lower triangular part of A.

//  tril(A,K) returns the elements on and below the K-th diagonal of
//  A.

//  K = 0: main diagonal
//  K > 0: above the main diag.
//  K < 0: below the main diag.

//  See Also: triu

//-------------------------------------------------------------------//

tril = function(x, k) 
{
  if (!exist (k)) { k = 0; }
  nr = x.nr; nc = x.nc;
  if(k > 0) 
  { 
    if (k > (nc - 1)) { error ("tril: invalid value for k"); }
  } else {
    if (abs (k) > (nr - 1)) { error ("tril: invalid value for k"); }
  }

  if (x.storage == "dense")
  {
    y = zeros(nr, nc);
    } else {
    y = spconvert ([nr, nc, 0]);
  }

  for(i in max( [1,1-k] ):nr) 
  {
    j = 1:min( [nc, i+k] );
    y[i;j] = x[i;j];
  }

  return y;
};
