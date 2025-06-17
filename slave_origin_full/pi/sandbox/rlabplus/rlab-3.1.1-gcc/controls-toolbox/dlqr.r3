//--------------------------------------------------------------------------
//
// dlqr
//
// syntax:  </e;k;s/> = dlqr(a,b,q,r,nn)
//
//	Linear quadratic regulator design for discrete-time systems.
//	</E;K;S/> = dlqr(A,B,Q,R)  calculates the optimal feedback gain
//	matrix K such that the feedback law  u[n] = -Kx[n]  minimizes the
//	cost function
//		J = Sum {x'Qx + u'Ru}
//	subject to the constraint equation:
//		x[n+1] = Ax[n] + Bu[n]
//
//	Also returned is S, the steady-state solution to the associated
//	discrete matrix Riccati equation and the closed loop eigenvalues
//	E,                            -1
//           0 = S - A'SA + A'SB(R+B'SB) BS'A - Q      E=eig(A-B*K).val[:]
//
//	</E;K;S/> = dlqr(A,B,Q,R,N) includes the cross-term N that relates
//	u to x in the cost function:
//
//		J = Sum {x'Qx + u'Ru + 2*x'Nu}
//
//	The controller can be formed with dreg().
//
//	See also: dlqry, lqrd, dreg
//
//--------------------------------------------------------------------------
require nargchk abcdchk

dlqr=function(a,b,q,r,nn)
{
  global(eps)

  msg = nargchk(4,5,nargs);
  if (msg !="") { error(msg); }
  msg = abcdchk(a,b);
  if (msg !="") { error(msg); }

  if (!length(a) || !length(b)) {
	error("A and B matrices cannot be empty.");
  }

  m  = a.nr;
  n  = a.nc;
  mb = b.nr;
  nb = b.nc;
  mq = q.nr;
  nq = q.nc;
  if ((m != mq) || (n != nq)) {
	error("A and Q must be the same size");
  }
  mr = r.nr;
  nr = r.nc;
  if ((mr != nr) || (nb != mr)) {
	error("B and R must be consistent");
  }

  if (nargs==5) {
    mn = nn.nr;
    nnn= nn.nc;
    if ((mn != m) || (nnn != nr)) {
       error("N must be consistent with Q and R");
    }
    // Add cross term
    q = q - nn/r*nn';
    a = a - b/r*nn';
  } else {
    nn = zeros(m,nb);
  }

  // Check if q is positive semi-definite and symmetric
  nq = norm(q,"1");
  if (any(eig(q).val < -eps*nq) || (norm(q'-q,"1")/nq > eps)) {
	disp("Warning: Q is not symmetric and positive semi-definite");
  }
  // Check if r is positive definite and symmetric
  nr = norm(r,1);
  if (any(eig(r).val <= -eps*nr) || (norm(r'-r,"1")/nr > eps)) {
	disp("Warning: R is not symmetric and positive definite");
  }

  // eigenvectors of Hamiltonian
  tmp = eig([a+b/r*b'/a'*q,  -b/r*b'/a';
               -a'\q,          inv(a)']);
  v = tmp.vec;
  d = tmp.val[:];
  tmp = sort(abs(d)); // sort on magnitude of eigenvalues
  e = tmp.val;
  index = tmp.idx;
  if (!((e[n] < 1) && (e[n+1]>1))) {
	error("Can''t order eigenvalues, (A,B) may be uncontrollable.");
  } else {
    e = d[index[1:n]];
  }
  // select vectors with eigenvalues inside unit circle
  chi = v[1:n;index[1:n]];
  lambda = v[(n+1):(2*n);index[1:n]];
  s = real(lambda/chi);
  k = (r+b'*s*b)\b'*s*a + r\nn';

  return <<e=e;k=k;s=s>>;
};
