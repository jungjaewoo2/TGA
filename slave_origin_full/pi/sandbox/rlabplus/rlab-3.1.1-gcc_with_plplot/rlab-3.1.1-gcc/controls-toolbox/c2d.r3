//----------------------------------------------------------------------
//
// c2d
//
// syntax: </ gamma; phi /> = c2d(a, b, t)
//
// Conversion of state space models from continuous to discrete time.
// </Gamma;Phi/> = c2d(A,B,T)  converts the continuous-time system:
//	.
//	x = Ax + Bu
//
// to the discrete-time state-space system:
//
// 	x[n+1] = Phi * x[n] + Gamma * u[n]
//
// assuming a zero-order hold on the inputs and sample time T.
//
// See also: c2dm, d2c
//
//----------------------------------------------------------------------
require nargchk abcdchk expm

c2d = function(a, b, t)
{
  msg = nargchk(3,3,nargs);
  if (msg!="") {error(msg);}
  msg = abcdchk(a,b);
  if (msg!="") {error(msg);}

  n = a.nc;
  nb= b.nc;
  s = expm([[a, b]*t; zeros(nb,n+nb)]);
  phi = s[1:n;1:n];
  gamma = s[1:n;n+1:n+nb];

  return <<gamma=gamma; phi=phi>>;
};
