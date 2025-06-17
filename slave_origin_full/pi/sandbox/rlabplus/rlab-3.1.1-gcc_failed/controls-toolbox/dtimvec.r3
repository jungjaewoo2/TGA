//----------------------------------------------------------------------
//
// dtimvec
//
// syntax: npts = dtimvec(a,b,c,x0,st,precision)
//
// returns a suggestion for the number of samples, npts that should be 
// used with discrete state space systems (a,b,c) that have known 
// initial conditions x0. npts can be used, for instance, for step 
// and impulse responses. 
//	
// dtimvec attempts to produce a number of points that ensures the 
// output responses have decayed to approximate st% of their initial 
// or peak values. 
//
//----------------------------------------------------------------------
require chop

dtimvec = function(a,b,c,x0,st)
{
  global(eps,pi)
  
  n = c.nc;
  // Work out scale factors by looking at d.c. gain, 
  // eigenvectors and eigenvalues of the system.
  tmp = eig(a);
  m = tmp.vec;
  r = tmp.val[:];
  // Equate the  responses to the sum of a set of 
  // exponentials(r) multiplied by a vector of magnitude terms(vec).
  r1 = c.nr;
  if (r1>1) { c1=max(abs(c)); } else { c1=c; }
  // Cater for the case when m is singular
  if (rcond(m)<eps) { 
    vec = (c1*m).'.*(1.0e4*ones(n,1));
  } else {
    vec = (c1*m).'.*(m\x0);
  }
  // Cater for the case when poles and zeros cancel each other out:
  vec = vec + 1e-20*(vec==0);
  // d.c. gain (not including d matrix)
  dcgain =c1*x0; // or dcgain=c\(eye(n)-a)*b;
  ind = find(imag(r)>0);
  // If the d.c gain is small then base the scaling
  // factor on the estimated maximum peak in the response.
  if (abs(dcgain)<1e-8) {
    // Estimation of the maximum peak in the response.
    peak = 2*imag(vec[ind]).*exp(-abs(0.5*real(r[ind]).*pi./imag(r[ind])));
    pk = max([abs(dcgain);abs(peak);eps]);
  } else {
    pk = dcgain;
  }
  // Work out the st// settling time for the responses.
  // (or st// peak value settling time for low d.c. gain systems).
  lt = st*pk./abs(vec);
  n1 = log(lt)./log(abs(r)+(abs(r)==1));
  n  = max(real(n1));
  // For unstable problems n will be negative
  if (n<=3) { n = max(abs(n1)) + 2; }
  n = min([n,1000]);
  // Round the maximum time to an appropriate value for plotting.
  nn = chop(n,1,5);
  if (abs((n-nn)/n)>0.2) { 
     nn = chop(n,1,1);
     if (abs((n-nn)/n)>0.2) {
        nn = chop(n,2,2);
     }
  }
  npts=nn+1;
  
  return npts;
};

