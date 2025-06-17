// -------------------------------------------------------------------------
//
// chop
//
// Syntax: X=chop(Xin,n,unit)
//
//       chop rounds the elements of X to n significant figures.
//       
//       e.g. chop(3.141592,5) returns 3.141600000..
//
//      chop(X,n,unit) rounds the elements of X to n significant
// 	figures whose digits (mantissa) are exactly divisible by unit. 
//
//       e.g. chop(3.141592,3,5) returns 3.150000000..
//            chop(3.141592,3,3) returns 3.150000000..
//            chop(3.141592,3,2) returns 3.140000000..
//
//--------------------------------------------------------------------------	   

chop = function(Xin,n,unit)
{
  // Set last sig. fig. rounding to 1 if only two input arguments.
  if (nargs<3) { unit=1; }

  // Cater for -ve numbers  and numbers = 0.
  X = abs(Xin) + (Xin==0);
  nx = X.nr;
  mx = X.nc;
  exponent = unit.*((10*ones(nx,mx)).^(floor(log10(X))-n+1));
  X = round(X./exponent).*exponent;

  // Put back sign and zeros
  X = sign(Xin).*X.*(Xin!=0);
  return X;
};

