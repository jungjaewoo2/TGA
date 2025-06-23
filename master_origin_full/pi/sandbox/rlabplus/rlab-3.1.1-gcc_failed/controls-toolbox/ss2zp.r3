//----------------------------------------------------------------------
//
// ss2zp
//
// Syntax: </k;p;z/> = ss2zp(A,B,C,D,IU) 
//
// </K;P;Z/> = ss2zp(A,B,C,D,IU)  calculates the transfer function in
// factored form:
//
//		     -1           (s-z1)(s-z2)...(s-zn)
//	H(s) = C(sI-A) B + D =  k ---------------------
//		                  (s-p1)(s-p2)...(s-pn)
// of the system:
//	.
//	x = Ax + Bu
//	y = Cx + Du
//
// from the single input IU.  The vector P contains the pole 
// locations of the denominator of the transfer function.  The 
// numerator zeros are returned in the columns of matrix Z with as 
// many columns as there are outputs y.  The gains for each numerator
// transfer function are returned in column vector K.
//
// See also: zp2ss, pzmap, tzero.
//----------------------------------------------------------------------
require abcdchk isempty tzero

ss2zp = function (a,b,c,d,iu)
{
  
  if (nargs < 4)
  {
     error("ss2zp: need at least 4 arguments.");
  }
  msg = abcdchk(a,b,c,d);
  if (msg != "") { error(msg); }
  nx = a.nr;
  ns = a.nc;

  if (nargs==4)
  {
     if (nx>0)
     {
        nb = b.nr;
        nu = b.nc;        
     } else {
        ny = d.nr;
        nu = d.nc;
     }
     if (nu<=1)
     {
        iu = 1;
     } else {
        error("IU must be specified for systems with more than one input.");
     }
  }

  // Remove relevant input:
  if (!isempty(b)) { b = b[;iu]; }
  if (!isempty(d)) { d = d[;iu]; }

  // Trap gain-only models
  if (nx==0 && !isempty(d))
  { 
     return <<z = []; p = []; k = d>>; 
  }

  // Do poles first, they're easy:
  p = eig(a).val[:];

  k = []; 
  // Compute zeros and gains using transmission zero calculation
  // Check if Control Toolbox is on path, in which case use the 
  // more accurate tzero method. 
  ny = d.nr;
  nu = d.nc;
  z = [];
  k = [];
  for (i in 1:ny)
  {
      tmp = tzero(a,b,c[i;],d[i;]);
      zi = tmp.z[:];
      gi = tmp.gain;
      mz = z.nr;
      nz = z.nc;
      nzi = length(zi);
      z = [[z;inf()*ones(max(0,nzi-mz),nz)],[zi;inf()*ones(max(0,mz-nzi),1)]];
      k = [k;gi];
  }
  // Now finish up by finding gains using Markov parameters
  if (isempty(k))
  {
     k = d;  
     CAn = c; 
     iter = 0;
     while (any(k==0))  // do until all k's are finished
     {
       i = find(k==0);
       if (iter > ns)
       {
          // Note: iter count is needed because when B is all zeros
          // it does not otherwise converge.
          k[i] = zeros(max(size(i)),1);
          break;
       }
       iter = iter + 1;
       markov = CAn*b;
       k[i] = markov[i];
       CAn = CAn*a;
     }
  }
  return <<z=z; p=p; k=k>>;
};


