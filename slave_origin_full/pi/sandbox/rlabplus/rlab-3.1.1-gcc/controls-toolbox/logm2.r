//----------------------------------------------------------------------
// logm2
//
// syntax: B  = logm2(A, pert)
//         B.Y, B.warning
//
// logm2(X) is the matrix natural logarithm of  X .    Complex
// results  are  produced  if  X  has nonpositive eigenvalues.
// logm2  can be thought of as being computed using eigenvalues
// D  and eigenvectors  V,  such that if  </D;V/> = eig(X)  then  
// LOGM2(X) = V*log(D)/V.  
//
//----------------------------------------------------------------------
require rsf2csf norm

logm2 = function(A, pert)
{
  global (eps)
  
  warning = 0;
  if (!length(A)) {
     Y = []; 
     return << Y=Y; warning=warning>>; 
  }
  tmp = schur(A);
  Q = tmp.z;
  T = tmp.t;
  tmp = rsf2csf(tmp.z,tmp.t);
  Q = tmp.u; 
  T = tmp.t;
  dT = diag(T);
  F = log(diag(T));

  // Set log of zero to large negative number
  find_vec = find(!finite(F));
  F[find_vec] = -ones(length(find_vec),1).*(1/eps);

  F = diag(F);
  tol = eps*norm(T,"1");
  n = max(size(A));
  for (p in 1:n-1) {
    for (i in 1:n-p) {
      j = i + p;
      s = T[i;j]*(F[j;j]-F[i;i]);
      if (p > 1) {
         k = i+1:j-1;
         s = s + T[i;k]*F[k;j] - F[i;k]*T[k;j];
      }
      d = T[j;j] - T[i;i];
      if (abs(d) <= tol) {
	 warning = 1;
         d = tol;
      }
      F[i;j] = s/d;
    }
  }
  Y = Q*F*Q';

  // If diagonal elements are the same then check accurarcy against expm.
  // If accuracy is poor then perturb matrix and call logm2 again.
  if (warning) {
     tol = tol + eps;
     if  (norm(expm(Y) - A, "I") > 1e5*tol) { 
	if (nargs ==1) { 
           pert = tol; 
	else
           pert = pert * 100;
	}
	tmp = $self(A + pert*rand("normal",n,n), pert);
	Y = tmp.Y;
	warning = tmp.warning;
     }
  } 
  return <<Y=Y; warning=warning>>; 
};

