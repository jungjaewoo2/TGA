//--------------------------------------------------------------------------
//
// zp2tf
//
// syntax: </den;num/> = zp2tf(z,p,k)
//
//	Zero-pole to transfer function conversion.
//	</DEN;NUM/> = zp2tf(Z,P,K)  forms the transfer function:
//
//	                NUM(s)
//	        H(s) = -------- 
//	                DEN(s)
//
//	given a set of zero locations in vector Z, a set of pole locations
//	in vector P, and a gain in scalar K.  Vectors NUM and DEN are 
//	returned with numerator and denominator coefficients in descending
//	powers of s.  
//
//	see also: tf2zp
//
//--------------------------------------------------------------------------
require isempty poly

zp2tf = function (z,p,k)
{
  local (z,p,k)
  
  den = real(poly(p[:]));
  md  = den.nr;
  nd  = den.nc;
  k   = k[:];
  mk  = k.nr;
  nk  = k.nc;
  if (isempty(z)) { 
     num = [zeros(mk,nd-1),k];
     return <<den=den;num=num>>;
  }
  m = z.nr;
  n = z.nc;
  if (mk != n) {
     if (m == 1) {
        error("z and p must be column vectors.");
     }
     error("k must have as many elements as z has columns.");
  }
  for (j in 1:n) {
    zj = z[;j];
    pj = real(poly(zj)*k[j]);
    num[j;] = [zeros(1,nd-length(pj)), pj];
  }
  
  return <<den=den;num=num>>;
};
