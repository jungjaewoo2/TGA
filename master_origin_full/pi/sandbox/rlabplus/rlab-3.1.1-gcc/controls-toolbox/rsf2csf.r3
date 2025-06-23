//----------------------------------------------------------------------
// rsf2csf
//
// Syntax: </t;u/> = rsf2csf(U,T)
//
//      </T;U/> = rsf2csf(U,T) converts a real, upper quasi-triangular
//	Schur form to a complex, upper triangular Schur form.
//
//----------------------------------------------------------------------
rsf2csf = function(U,T)
{
  // Find complex unitary similarities to zero subdiagonal elements.
  n = max(size(T));
  m = n;
  while (m > 1)
  {
    s = abs(T[m-1;m-1]) + abs(T[m;m]);
    if (s + abs(T[m;m-1]) > s)
    {
      k = m-1:m;
      mu = eig(T[k;k]).val - T[m;m];
      r = norm([mu[1], T[m;m-1]],"2");
      c = mu[1]/r;  s = T[m;m-1]/r;
      G = [c', s; -s, c];
      j = m-1:n;  T[k;j] = G*T[k;j];
      i = 1:m;    T[i;k] = T[i;k]*G';
      i = 1:n;    U[i;k] = U[i;k]*G';
    }
    T[m;m-1] = 0;
    m = m-1;
  }
  return <<u=U;t=T>>;
};

