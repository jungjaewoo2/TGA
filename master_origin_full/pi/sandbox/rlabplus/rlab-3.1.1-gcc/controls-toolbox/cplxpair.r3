//--------------------------------------------------------------------------
//
// cplxpair
//
// Syntax: z = cplxpair(x,tol)
//
//	Sort numbers into complex conjugate pairs.
//	z = cplxpair(x) rearranges the elements of vector X so that
//	complex numbers are collected into matched pairs of complex
//	conjugates.  The pairs are ordered by increasing real part.
//	Any purely real elements are placed after all the complex pairs.
//	z = cplxpair(x,tol) uses a relative tolerance of tol for
//	comparison purposes.  The default is tol = 100*eps.
//
//--------------------------------------------------------------------------

cplxpair = function (x,tol)
{
  global(eps)

  if (nargs < 2) {
     tol = 100*eps; // need to give a tolerance more generous than eps
  }
  j = sqrt(-1);
  m = x.nr;
  n = x.nc;
  z  = zeros(m,n); // return number in same shape vector as input
  x = x[:];	   // make sure x is a column vector
  nz = max(size(z));

  // first segregate reals
  ir = find(abs(imag(x)) <= tol*abs(x));
  nr = max(size(ir));
  if (!isempty(ir)) {
	r = sort(real(x[ir])).val;   // return sorted reals first
	x[ir] = [];	             // remove these values from input array
  }
  clear(ir);

  nc = max(size(x));
  if (nc == 0) {
       // no complex numbers
       z = r;
       return z;
  }
  if (mod(nc,2)) {
     error("Complex numbers can't be paired");
  }

  // copy complex x to a work vector c
  c = x;
  np = nc/2;	// number of pairs supposedly

  // get values to upper half plane
  c = real(c) + j*abs(imag(c));

  // ordered by real part (since sort is very sensitive to magnitudes)
  tmp = sort(real(c));
  cc = tmp.val;
  cind = tmp.idx;
  c = cc +j*imag(c[cind]);
  // check to see if real parts are at least in pair
  if (any(abs(cc[1:nc:2]-cc[2:nc:2]) > tol*abs(c[1:nc:2]))) {
     error("Complex numbers can't be paired");
  }
  x = x[cind];	// reorder x the same way c has been
  clear (cc);

  // check real part pairs to see if imag parts are paired by conjugates
  // be careful with multiple roots!
  ip = 1;	// initialize pair counter
  while (ip <= np)
  {
    // find indices for same real parts - but be careful because a real part can be 0
    ii = find(abs(real(c[1:nc:2])-real(c[2*ip])) <= tol*abs(c[2*ip]));
    if (isempty(ii))
    {
       error("Complex numbers can't be paired");
    }
    if (max(size(ii)) > 1)
    {
       // multiple pairs with same real part - sort on imag(c(ii)) for all with real part
       // ij below are indices with same real part
       ij = find(abs(real(c[1:nc]) - real(c[2*ip])) <= tol*abs(c[2*ip]));
       nn = max(size(ij));
       xtemp = x[ij];
       xind  = sort(imag(xtemp)).idx; //sort on imag parts - should be paired
       xtemp = xtemp[xind];
       // check pairing
       if (any((abs(xtemp-conj(xtemp[nn:1:-1])) > tol*abs(xtemp))))
       {
          error("Complex numbers can't be paired")
       } else {
          x[ij[1]:ij[nn-1]:2] = xtemp[1:nn/2];
          x[ij[2]:ij[nn]:2] = conj(xtemp[1:nn/2]);
       }
    } else {
       // only one pair with that real part
       if (x[2*ip]-conj(x[2*ip-1]) > tol*abs(x[2*ip]))
       {
          error("Complex numbers can't be paired")
       } else {
          xtmp = real(x[2*ip])-sqrt(-1)*abs(imag(x[2*ip]));
          x[2*ip-1:2*ip] = [xtmp, conj(xtmp)];
       }
     }
     ip = ip + max(size(ii)); // increment pair counter
  }

  // copy complex pairs into return vector
  z[1:nc] = x[1:nc];
  // append reals to this
  if(nr>0) { z[nc+1:nc+nr] = r; }

  return z;
};

