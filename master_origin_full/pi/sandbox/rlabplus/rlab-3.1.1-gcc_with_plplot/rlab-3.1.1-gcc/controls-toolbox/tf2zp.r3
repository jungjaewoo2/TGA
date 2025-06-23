//----------------------------------------------------------------------
//
// tf2zp
//
// syntax:  </K;P;Z/> = tf2zp(NUM,DEN)
//
//	Transfer function to zero-pole conversion.
//
//	</K;P;Z/> = tf2zp(NUM,DEN) finds the zeros, poles, and gains:
//
//		          (s-z1)(s-z2)...(s-zn)
//		H(s) =  K ---------------------
//		          (s-p1)(s-p2)...(s-pn)
//
//	from a SIMO transfer function in polynomial form:
//
//		        NUM(s)
//		H(s) = -------- 
//		        den(s)
//
//	Vector DEN specifies the coefficients of the denominator in 
//	descending powers of s.  Matrix NUM indicates the numerator 
//	coefficients with as many rows as there are outputs.  The zero
//	locations are returned in the columns of matrix Z, with as many 
//	columns as there are rows in NUM.  The pole locations are returned
//	in column vector P, and the gains for each numerator transfer 
//	function in vector K. 
//	
//	See also: zp2tf
//
//----------------------------------------------------------------------
require roots tfchk 

tf2zp = function (num, den)
{
  global (eps)
  
  tmp = tfchk(num,den);
  num = tmp.numc;
  den = tmp.denc;

  // Normalize transfer function
  if (length(den))
  {
     coef = den[1];
  } else {
     coef = 1;
  }
  if (abs(coef)<eps)
  {
    error("Denominator must have non-zero leading coefficient.");
  }
  den = den./coef;
  num = num./coef;

  // Remove leading columns of zeros from numerator
  if (length(num))
  {
    while(all(num[;1]==0))
    {
      num[;1] = [];
    }
  }
  ny = num.nr;
  np = num.nc;

  // Poles
  p  = roots(den);
  // Zeros and Gain
  z = inf()*ones(np-1,ny);
  for (i in 1:ny)
  {
    zz = roots(num[i;]);
    if (length(zz))
    {
       z[1:length(zz); i] = zz; 
    }
    ndx = find(num[i;] != 0);
    if (length(ndx))
    {
       k[i;1] = num[i;ndx[1]]; 
    }
  }
  return <<z=z;p=p;k=k>>;
};


