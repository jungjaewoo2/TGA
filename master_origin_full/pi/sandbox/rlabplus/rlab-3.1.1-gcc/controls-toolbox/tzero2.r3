//---------------------------------------------------------------------
//
// tzero2
//
// syntax:  z = tzero(A,B,C,D)
//
//	tzero2(A,B,C,D) returns the transmission zeros of the state-space
//	system:
//		.
//		x = Ax + Bu
//		y = Cx + Du
//
//	Care must be exercised to disregard any large zeros that may
//	actually be at infinity.  
//
//	This rfile finds the transmission zeros using an algorithm based
//	on QZ techniques.  It is not as numerically reliable as the 
//	algorithm in tzero.
//
//	See also: tzero
//
// For more information on this algorithm, see:
//   A.J Laub and B.C. Moore, Calculation of Transmission Zeros Using
//   QZ Techniques, Automatica, 14, 1978, p557
//
// For a better algorithm, not implemented here, see:
//   A. Emami-Naeini and P. Van Dooren, Computation of Zeros of Linear 
//   Multivariable Systems, Automatica, Vol. 14, No. 4, pp. 415-430, 1982.
//
//---------------------------------------------------------------------
require abcdchk

tzero = function(a,b,c,d)
{
  if (nargs != 4) {
     error("Need 4 arguments");
  }
  error(abcdchk(a,b,c,d));

  n = b.nr; m = b.nc;
  r = c.nr; n = c.nc;
  aa = [a,b;c,d];

  if (m == r)
  {
     bb = zeros(aa);
     bb[1:n;1:n] = eye(n,n);
     z = eig(aa,bb).val;
     z = z[find(finite(z))]; // Punch out NaN's and Inf's
  } else {
     nrm = norm(aa,"1");
     if (m > r)
     {
        aa1 = [aa;nrm*(rand(m-r,n+m)-.5)];
        aa2 = [aa;nrm*(rand(m-r,n+m)-.5)];
     } else {
        aa1 = [aa, nrm*(rand(n+r,r-m)-.5)];
        aa2 = [aa, nrm*(rand(n+r,r-m)-.5)];
     }
     bb = zeros(aa1);
     bb[1:n;1:n] = eye(n,n);
     z1 = eig(aa1,bb).val;
     z2 = eig(aa2,bb).val;
     z1 = z1[find(finite(z1))];   // Punch out NaN's and Inf's
     z2 = z2[find(finite(z2))];
     nz = length(z1);
     for (i in 1:nz)
     {
        if (any(abs(z1[i]-z2) < nrm*sqrt(eps)))
        {
           z = [z;z1[i]];
        }
     }
  }
  return z[:];
};

