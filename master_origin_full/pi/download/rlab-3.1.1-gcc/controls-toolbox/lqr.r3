//----------------------------------------------------------------------
//
// lqr
//
// syntax: </e;k;s/> = lqr(a,b,q,r,nn)
//
// Linear quadratic regulator design for continuous systems.
// </E;K;S/> = lqr(A,B,Q,R)  calculates the optimal feedback gain
// matrix K such that the feedback law  u = -Kx  minimizes the cost
// function:
//	J = Integral {x'Qx + u'Ru} dt
//
// subject to the constraint equation:
//	.
//	x = Ax + Bu
//
// Also returned is S, the steady-state solution to the associated
// algebraic Riccati equation and the closed loop eigenvalues E:
//			  -1
//	0 = SA + A'S - SBR  B'S + Q     E = eig(A-B*K).val[:]
//
// </E;K;S/> = lqr(A,B,Q,R,N) includes the cross-term N that relates
// u to x in the cost function.
//
//	J = Integral {x'Qx + u'Ru + 2*x'Nu}
//
// The controller can be formed with reg.
//
// See also: lqry, lqr2, and reg.
//
//----------------------------------------------------------------------
require abcdchk norm

lqr=function(a,b,q,r,nn)
{

  global(eps,pi)

  if (nargs < 4 || nargs > 5) { error("Wrong number og arguments."); }
  msg = abcdchk(a,b);
  if (msg!="") { error(msg); }
  if (!length(a) || !length(b)) {
     error("A and B matrices cannot be empty.");
  }

  m  = a.nr;
  n  = a.nc;
  mb = b.nr;
  nb = b.nc;
  mq = q.nr;
  nq = q.nc;
  if (m != mq || n != nq) {
     error("A and Q must be the same size");
  }
  mr = r.nr; nr = r.nc;
  if (mr != nr || nb != mr) {
     error("B and R must be consistent");
  }

  if (nargs == 5) {
     mn = nn.nr; nnn = nn.nc;
     if (mn != m || nnn != nr) {
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
  nr = norm(r,"1");
  if (any(eig(r).val <= -eps*nr) || (norm(r'-r,"1")/nr > eps)) {
     disp("Warning: R is not symmetric and positive definite");
  }

  // Start eigenvector decomposition by finding eigenvectors of Hamiltonian:
  tmp = eig([a, b/r*b'; q, -a']);
  v = tmp.vec;
  d = tmp.val[:];
  tmp = sort(real(d));  //sort on real part of eigenvalues
  e = tmp.val;
  index = tmp.idx;
  if (!( (e[n]<0) && (e[n+1]>0) )) {
    error("Can't order eigenvalues, (A,B) may be uncontrollable.");
  } else {
    e = d[index[1:n]];		 // Return closed-loop eigenvalues
  }
  chi = v[1:n;index[1:n]];	 // select vectors with negative eigenvalues
  lambda = v[(n+1):(2*n);index[1:n]];
  s = -real(lambda/chi);
  k = r\(nn'+b'*s);

  return <<k=k; s=s; e=e>>;
};
