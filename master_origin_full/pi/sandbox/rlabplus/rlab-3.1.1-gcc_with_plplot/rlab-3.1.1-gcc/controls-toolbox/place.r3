//----------------------------------------------------------------------
// place
//
// syntax:  K = place(A,B,P)
//
//	K = place(A,B,P)  computes the state feedback matrix K such that
//	the eigenvalues of  A-B*K  are those specified in vector P.
//	The complex eigenvalues in the vector P must appear in consecu-
//	tive complex conjugate pairs. No eigenvalue may be placed with
//	multiplicity greater than the number of inputs.  
//
//	The  displayed "ndigits" is an  estimate of how well the
//	eigenvalues were placed.   The value seems to give an estimate
//	of how many decimal digits in the eigenvalues of A-B*K match
//	the specified numbers given in the array P.
//
//	A warning message is printed if the nonzero closed loop poles
//	are greater than 10% from the desired locations specified in P.
//
//	See also: lqr and rlocus
//
// Ref: Kautsky, Nichols, Van Dooren, "Robust Pole Assignment in Linear 
//      State Feedback," Intl. J. Control, 41(1985)5, pp 1129-1155
//
//----------------------------------------------------------------------
place = function(A, B, P)
{
  global(eps);

  NTRY = 5;         // - number of iterations for "optimization" -
  nx = A.nr;
  na = A.nc;

  P = P[:];
  if (length(P) != nx)
  {
     error("P must have the same number of states as A."); 
  }

  n = B.nr;
  m = B.nc;
  if (nx == 0 || n == 0)
  {
     error("A and B matrices cannot be empty.");
  }
  nx = 0; 
  i = 1; pr = []; pi = [];
  while (i <= n)
  {
    if (imag(P[i]) != 0.0)
    {
       pr = [pr,real(P[i])]; pi = [pi, imag(P[i])];
       cmplx = [cmplx,1]; i = i + 2;
    } else {
       pr = [pr,real(P[i])]; pi = [pi, 0.0];
       cmplx = [cmplx,0]; i = i + 1;
    }
    nx = nx + 1;
  }

  m = rank(B);
  // Make sure there are more inputs than repeated poles:
  ps = sort(P).val;
  for (i in 1:n-m)
  {
      imax = min(n,i+m);
      if (all(ps[i:imax] == ps[i]))
      {
         error("Can''t place poles with multiplicity greater than the number of inputs.");
      }
  }
  nmmp1 = n-m+1; mp1 = m+1; jj = sqrt(-1);
  // [Qb,Rb] = qr(B);
  QRb = qr(B);
  q0  = QRb.q[;1:m]; 
  q1  = QRb.q[;mp1:n]; 
  Rb  = QRb.r[1:m;];
  //
  // - special case: (#inputs)==(#states) - efficient, but not clean
  if (m==n)
  {
     A = A - diag(real(P));
     i = 0; 
     for (j in 1:nx)
     {
       i = i + 1;
       if (cmplx(j))
       {
         A[i;i+1] = A[i;i+1] + pi[j];
         A[i+1;i] = A[i+1;i] - pi[j];
         i = i+1;
       }
     }
     printf("place: ndigits= %g\n", fix(log10(1.0/eps)));
     K = Rb\q0'*A;
     return K;	
  }
  //
  // - compute bases for eigenvectors -
  I = eye(n,n);
  for (i in 1:nx)
  {
    //[Q,R] = qr(((pr[i]+jj*pi[i])*I-A)'*q1);
    QR = qr(((pr[i]+jj*pi[i])*I-A)'*q1);
    Bx = [ Bx, QR.q[;nmmp1:n] ];
  }
  //
  // - choose basis set -
  // at each iteration of i pick the eigenvector Xj, j~=i, 
  // which is "most orthogonal" to the current eigenvector Xi
  for (i in 1:nx) 
  { 
     X[;i] = Bx[;(i-1)*m+1]; 
  }
  if (m>1)
  {
     for (k in 1:NTRY)
     {
         for (i in 1:nx)
         {
           S  = [ X[;1:i-1], X[;i+1:nx] ]; 
           S  = [ S, conj(S) ];
           Us = svd(S).u;
           Pr = Bx[;(i-1)*m+1:i*m]; 
           Pr = Pr*Pr';
           X[;i] = Pr*Us[;n]; 
           X[;i] = X[;i]/norm(X[;i]);
         }
     }
  }
  for (i in 1:nx)
  {
    if (cmplx(i))
    {
        Xf = [ Xf, X[;i], conj(X[;i]) ];
    } else {
        Xf = [ Xf, X[;i] ];
    }
  }
  cnd = cond(Xf);
  if (cnd*eps >= 1.0)
  {
    printf("place: can''t place eigenvalues there\n");
    return K;
  }
  printf("place: ndigits= %g\n", fix(log10(cnd/eps)));
  //
  // - compute feedback -
  K = Rb\q0'*(A-real(Xf*diag(P,0)/Xf));

  // Check results. Start by removing 0.0 pole locations
  P = sort(P).val;
  i = find(P != 0);
  P = P[i];
  Pc = sort(eig(A-B*K));
  Pc = Pc[i];
  if (max(abs(P-Pc)./abs(P)) > .1)
  {
     printf("Warning: Pole locations are more than 10% in error.\n");
  }

  return K;
};

