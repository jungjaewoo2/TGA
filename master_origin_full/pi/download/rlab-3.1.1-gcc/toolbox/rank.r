//-------------------------------------------------------------------//
//
//  Syntax:	rank ( A )
//		rank ( A , tol )

//  Description:

//  Compute the rank of the matrix A. Rank returns the number of
//  singular values that are larger than:
//        max( size(x) ) * norm(x,"2") * eps. 

//  If the user specifies tol, the the number of singular values
//  larger than tol is returned.

//-------------------------------------------------------------------//

rank = function(x, tol)
{
  global (eps)

  s = svd(x);
  if (!exist (tol))
  { 
    tol = max(size(x)) * norm(x,"2") * eps;
  }
  return sum(s.sigma > tol);
};
