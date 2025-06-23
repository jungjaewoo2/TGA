//--------------------------------------------------------------------------
//
// acker
//
// Syntax:  K = acker(A,B,P)
//
//	Pole placement gain selection using Ackermann's formula.
//	K = ACKER(A,B,P)  calculates the feedback gain matrix K such that
//	the single input system
//	        .
//	        x = Ax + Bu 
//
//	with a feedback law of  u = -Kx  has closed loop poles at the 
//	values specified in vector P, i.e.,  P = eig(A-B*K).
//
//	See also: place
//
//	Note: This algorithm uses Ackermann's formula.  This method
//	is NOT numerically reliable and starts to break down rapidly
//	for problems of order greater than 10, or for weakly controllable
//	systems.  A warning message is printed if the nonzero closed loop
//	poles are greater than 10% from the desired locations specified 
//	in P.
//
//	Algorithm is from page 201 of:
//	Kailath, T.  "Linear Systems", Prentice-Hall, 1980.
//--------------------------------------------------------------------------

acker = function (a, b, p)
{
  require abcdchk ctrb poly polyval
  
  if (nargs != 3) { error("Need 3 arguments"); }
  msg = abcdchk(a,b);
  if (msg != "") { error(msg); }

  m = b.nr;
  n = b.nc;
  if (n != 1) { error("System must be single input"); }
  p = p[:]; // Make sure roots are in a column vector
  mp = p.nr; 
  np = p.nc;
  m = a.nr;
  n = a.nc;
  if (m != mp) { error("Vector p must have size(a) elements"); }

  // Form gains using Ackerman's formula
  k = ctrb(a,b)\polyvalm(real(poly(p)),a);
  k = k[n;];

  // Check results. Start by removing 0.0 pole locations
  p = sort(p).val;
  i = find(p != 0.0);
  p = p[i];
  pc = sort(eig(a-b*k).val[:]).val;
  pc = pc[i];
  if (max(abs(p-pc)./abs(p)) > .1)
  {
     printf("Warning: Pole locations are more than 10% in error.\n");
  }
  return k;
};

