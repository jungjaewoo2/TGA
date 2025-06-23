//-------------------------------------------------------------------//
//
//  Syntax:	pinv ( A )
//		pinv ( A , tol )

//  Description:

//  The pinv function computes the pseudo-inverse of it's argument A. 

//  See Also: rank

//-------------------------------------------------------------------//

pinv = function(A, tol)
{
  // Pseudo-inverse, ignore singular values <= tol.
  // Default tol = max(size(A)) * s(1) * eps.

  s = svd(A);
  if(!exist (tol)) { tol = max(size(A)) * norm(A,"2") * epsilon(); }
 
  r = sum ((s.sigma > tol));
  if(r == 0) 
  { 
    X = zeros (size (A')); 
  } else {
    S = diag (ones (r,1) ./ s.sigma[1:r]');
    X = s.vt'[;1:r] * S * s.u[;1:r]';
  }
  return X;
};
