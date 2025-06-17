//-------------------------------------------------------------------//

//  Synopsis:   Return the upper triangular part of A.

//  Syntax:	triu ( A )
//		triu ( A , K )

//  Description:

//  triu(x) returns the upper triangular part of A.

//  tril(x; k) returns the elements on and above the k-th diagonal of
//  A. 

//  K = 0: main diagonal
//  K > 0: above the main diag.
//  K < 0: below the main diag.

//  See Also: tril

//-------------------------------------------------------------------//

triu = function(x, k) 
{
  if (!exist (k)) { k = 0; }
  nr = x.nr; nc = x.nc;

  if(k > 0) 
  { 
    if (k > (nc - 1)) { error ("triu: invalid value for k"); }
  } else {
    if (abs (k) > (nr - 1)) { error ("triu: invalid value for k"); }
  }

  if (x.storage == "dense")
  {
    y = zeros(nr, nc);
    } else {
    y = spconvert ([ nr, nc, 0 ]);
  }

  for(j in max( [1,1+k] ):nc) 
  {
    i = 1:min( [nr, j-k] );
    y[i;j] = x[i;j];
  }

  return y;
};
