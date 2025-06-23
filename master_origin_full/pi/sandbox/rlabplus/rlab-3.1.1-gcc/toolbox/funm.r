//-------------------------------------------------------------------//

//  Syntax:	funm ( A , FUN )

//  Description:

//  The funm function computes the matrix function of A as described
//  in scalar context by FUN. For example:

//	funm ( A , exp )

//  calculates the matrix exponential (using expm() is better though).

//  The second argument, FUN is NOT a string, it is a variable, in
//  this particular case, it is a variable representing either a
//  builtin, or user-function.

//  Funm uses Parlett's method. See Golub and VanLoan (1983), p. 384.

//  See Also: expm, logm
//-------------------------------------------------------------------//

funm = function ( A , fun )
{
  global (eps)

  S = schur (A + zeros (size (A))*1i);	// get Complex Schur form
  F = diag (fun (diag (S.t)));
  tol = eps*norm (S.t,"1");
  n = max (size (A));
  for (p in 1:n-1)
  {
    for (i in 1:n-p)
    {
      j = i+p;
      s = S.t[i;j]*(F[j;j] - F[i;i]);
      if (p > 1)
      {
        k = i+1:j-1;
        s = s + S.t[i;k]*F[k;j] - F[i;k]*S.t[k;j];
      }
      d = S.t[j;j] - S.t[i;i];
      if (abs (d) <= tol)
      {
        fprintf ("stderr", "WARNING: Result from FUNM is probably inaccurate.");
        d = tol;
      }
      F[i;j] = s/d;
    }
  }
  return S.z*F*S.z';
};
