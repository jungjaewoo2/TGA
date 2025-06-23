//---------------------------------------------------------------------
// 
// unwrap
//
// Syntax: y = unwrap(x)
//
//	Attempt to "unwrap" phase angles by adding +-2*Pi to them so
//	that instead of jumping from -Pi to Pi, they will have a
//	smooth transition across the branch cuts.
//	Y = unwrap(X) corrects the phase angles in vector X.   When X
//	is a matrix, the phase angles are corrected down each column.
//	The phase MUST be in radians.
//
//---------------------------------------------------------------------
unwrap = function(praw)
{
  global(eps,pi)

  m = praw.nr;
  n = praw.nc;
  mo = m;
  if (m == 1) {
     praw = praw.';
     m = praw.nr;
     n = praw.nc;
  }
  pfix = praw;
  closetopi = pi*170/180;    // Here is a tolerance.
  twopi = 2*pi;
  for (j in 1:n) {
	p = praw[;j];
	pd = diff(p);
	ijump = find(abs(pd) > closetopi);
	for (i in 1:max(size(ijump))) {
		k = ijump[i];
		p[k+1:m] = p[k+1:m] - twopi * sign(pd[k]);
	}
	pfix[;j] = p;
  }
  if (mo == 1) {
	pfix = pfix.';
  }
  return pfix;
};


