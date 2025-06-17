//----------------------------------------------------------------------
// orth
//
// syntax: Q = orth(A)
//
// Orthogonalization.  Q = orth(A) is an orthonormal basis for
// the range of  A .  Q'*Q = eye(A.nr,A.nc), the columns of Q span the
// same space as the columns of  A  and the number of  columns
// of Q is the rank of A. 
//
//----------------------------------------------------------------------
orth = function(A)
{
  global(eps,pi)
  
  // QR decomposition
  tmp = qr(A,"p");
  Q = tmp.q;
  R = tmp.r;
  E = tmp.p;

  // Determine r = effective rank
  tol = eps*norm(A,"f");
  r = sum(abs(diag(R)) > tol);
  r = r[1]; // fix for case where R is vector.
  // Use first r columns of Q.
  if (r > 0) {
     Q = Q[;1:r];
     // Cosmetic sign adjustment
     Q = -Q;
     Q[;r] = -Q[;r];
  else
     Q = [];
  }
  return Q;
};


