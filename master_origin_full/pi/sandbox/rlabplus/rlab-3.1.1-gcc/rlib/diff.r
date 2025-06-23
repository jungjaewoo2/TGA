//-------------------------------------------------------------------//

//  Synopsis:   Compute the difference between vector elements.

//  Syntax:	diff ( X )
//		diff ( X , k )

//  Description:

//  The function diff computes the difference between adjacent
//  elements in a vector. The function diff merely does:

//	x = x[2:n] - x[1:n-1];

//  If X is a matrix, then the operation is performed on each column
//  of X.

//  The default value for `k' is 1. If k is specified as other than 1,
//  then the operation is performed k times.

//-------------------------------------------------------------------//

diff = function ( X, k )
{
  local ( X , k )

  if (!exist (k))
  { k = 1; }

  for (i in 1:k)
  {
    m = X.nr; n = X.nc;
    if (m == 1 || n==1)
    {
      _isfin = find(finite(X) || X==inf() || X==-inf());
      if (isempty(_isfin))
      { return []; }
      n = length(_isfin);
      X = X[_isfin];
      X = X[2:n] - X[1:n-1];
    else
      if(any(any(isnan(X))))
      { return []; }
      X = X[2:m;] - X[1:m-1;];
    }
  }
  return X;
};
