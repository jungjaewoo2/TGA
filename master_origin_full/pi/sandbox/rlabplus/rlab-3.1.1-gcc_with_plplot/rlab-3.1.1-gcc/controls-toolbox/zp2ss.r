//----------------------------------------------------------------------
//
// zp2ss
//   
// Syntax:  </A;B;C;D/> = zp2ss(z,p,k)
//
//      Zero-pole to state-space conversion.
//
//	The rfile calculates a state-space representation:
//		.
//		x = Ax + Bu
//		y = Cx + Du
//
//	for SIMO systems given a set of pole locations in column vector p,
//	a matrix z with the zero locations in as many columns as there are
//	outputs, and the gains for each numerator transfer function in
//	vector k.  The A,B,C,D matrices are returned in block diagonal
//	form.
//
//	See also: ss2zp
//
//----------------------------------------------------------------------
require cplxpair poly tf2ss zp2tf 

zp2ss = function(z,p,k)
{
  global(eps)
  
  m = z.nr;
  n = z.nc;
  if (n==0) {
     n=length(k);  // Fix to handle multi-output when z is empty
  }   
  if (length(k) != n && (!isempty(z))) {
     error("Z should be a column vector or K should be SIMO.");
  }
  if (n > 1)  {
     // If it's multi-output, we can't use the nice algorithm
     // that follows, so use the numerically unreliable method
     // of going through polynomial form, and then return.
     tmp = zp2tf(z,p,k);  // Suppress compile-time diagnostics
     return tf2ss(tmp.num,tmp.den);
  }

  // Strip infinities and throw away.
  if (p!=[]) { p = p[find(abs(p) != inf())];}
  if (z!=[]) { z = z[find(abs(z) != inf())];}

  // Group into complex pairs
  np = length(p);
  nz = length(z);
  if (z!=[]) {z = cplxpair(z,1e6*nz*norm(z,"2")*eps + eps);}
  if (p!=[]) {p = cplxpair(p,1e6*np*norm(p,"2")*eps + eps);}

  // Initialize state-space matrices for running series
  a=[]; b=[]; c=[]; d=1;

  // If odd number of poles AND zeros, convert the pole and zero
  // at the end into state-space.
  //	H(s) = (s-z1)/(s-p1) = (s + num(2)) / (s + den(2))
  if (mod(np,2) && mod(nz,2))
  {
	a = p[np];
	b = 1;
	c = p[np] - z[nz];
	d = 1;
	np = np - 1;
	nz = nz - 1;
  }

  // If odd number of poles only, convert the pole at the
  // end into state-space.
  //  H(s) = 1/(s-p1) = 1/(s + den(2)) 
  if (mod(np,2))
  {
	a = p[np];
	b = 1;
	c = 1;
	d = 0;
	np = np - 1;
  }	

  // If odd number of zeros only, convert the zero at the
  // end, along with a pole-pair into state-space.
  //   H(s) = (s+num(2))/(s^2+den(2)s+den(3)) 
  if (mod(nz,2))
  {
	num = real(poly(z[nz]));
	den = real(poly(p[np-1:np]));
	wn = sqrt(prod(abs(p[np-1:np])));
	if (wn == 0) { wn = 1; }
	t = diag([1, 1/wn]);	// Balancing transformation
	a = t\[-den[2], -den[3]; 1, 0]*t;
	b = t\[1; 0];
	c = [1, num[2]]*t;
	d = 0;
	nz = nz - 1;
	np = np - 2;
  }

  // Now we have an even number of poles and zeros, although not 
  // necessarily the same number - there may be more poles.
  //   H(s) = (s^2+num(2)s+num(3))/(s^2+den(2)s+den(3))
  // Loop thru rest of pairs, connecting in series to build the model.
  i = 1;
  while (i < nz)
  {
	index = i:i+1;
	num = real(poly(z[index]));
	den = real(poly(p[index]));
	wn = sqrt(prod(abs(p[index])));
	if (wn == 0) { wn = 1; }
	t = diag([1,1/wn]);    // Balancing transformation
	a1 = t\[-den[2], -den[3]; 1, 0]*t;
	b1 = t\[1; 0];
	c1 = [num[2]-den[2], num[3]-den[3]]*t;
	d1 = 1;
        // [a,b,c,d] = series(a,b,c,d,a1,b1,c1,d1); 
        // Next lines perform series connection
	a = [a,zeros(a.nr,a1.nc); b1*c, a1];
	b = [b; b1*d];
	c = [d1*c, c1];
	d = d1*d;
	i = i + 2;
  }

  // Take care of any left over unmatched pole pairs.
  //   H(s) = 1/(s^2+den(2)s+den(3))
  while (i < np)
  {
	den = real(poly(p[i:i+1]));
	wn = sqrt(prod(abs(p[i:i+1])));
	if (wn == 0) { wn = 1; }
	t = diag([1,1/wn]);	// Balancing transformation
	a1 = t\[-den[2],-den[3]; 1,0]*t;
	b1 = t\[1; 0];
	c1 = [0,1]*t;
	d1 = 0;
        // [a,b,c,d] = series(a,b,c,d,a1,b1,c1,d1);
        // Next lines perform series connection 
	a = [a,zeros(a.nr,a1.nc); b1*c,a1];
	b = [b; b1*d];
	c = [d1*c, c1];
	d = d1*d;

	i = i + 2;
  }

  // Apply gain k:
  c = c*k;
  d = d*k;
  
  return <<a=a;b=b;c=c;d=d>>;
  
};

