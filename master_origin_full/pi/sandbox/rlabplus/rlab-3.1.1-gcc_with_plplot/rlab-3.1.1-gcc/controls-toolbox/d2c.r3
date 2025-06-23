//--------------------------------------------------------------------------
// d2c
//
// Syntax: </a; b/> = d2c(phi, gamma, t)
//
//	Conversion of state space models from discrete to continuous time.
//	</A; B/> = D2C(Phi, Gamma, T)  converts the discrete-time system:
//
//		x[n+1] = Phi * x[n] + Gamma * u[n]
//
//	to the continuous-time state-space system:
//		.
//		x = Ax + Bu
//
//	assuming a zero-order hold on the inputs and sample time T.
//
//	See also: d2cm and c2d.
//--------------------------------------------------------------------------

require logm2 abcdchk norm

d2c = function(phi, gamma, t)
{
  global (eps)
  
  if (nargs != 3) { error("Wrong number of input arguments."); }
  error(abcdchk(phi,gamma));

  m = phi.nr; n = phi.nc;
  m = gamma.nr; nb = gamma.nc;

  // phi = 1 case cannot be computed through matrix logarithm.  Handle
  // as a special case.
  if (m == 1) {
	if (phi == 1) {
		a = 0; b = gamma/t;
		return <<a=a; b=b>>;
	}
  }

  // Remove rows in gamma that correspond to all zeros
  b = zeros(m,nb);
  nz = 0;
  nonzero = [];
  for (i in 1:nb)
  {
	if (any(gamma[;i]!=0)) { 
		nonzero = [nonzero, i];
		nz = nz + 1;
	}
  }

  // Do rest of cases using matrix logarithm.
 
  s = logm2([[phi, gamma[;nonzero]]; zeros(nz,n), eye(nz,nz)])/t;
  if (norm(imag(s),"I") > sqrt(eps)) { 
	disp("Warning: Accuracy of d2c conversion may be poor.")
  }
  s = real(s);
  a = s[1:n;1:n];
  if (length(b))
  {
	b[;nonzero] = s[1:n;n+1:n+nz];
  }
  return <<a=a;b=b>>;
};

