//--------------------------------------------------------------------------
//
// ddcgain
//
// Syntax: k = ddcgain(a,b,c,d)
//
//	D.C. gain of discrete system.
//
//	K = ddcgain(A,B,C,D) computes the steady state (D.C. or low 
//	frequency) gain of the discrete state-space system (A,B,C,D).
//
//	K = ddcgain(NUM,DEN) computes the steady state gain of the 
//	discrete polynomial transfer function system G(z) = NUM(z)/DEN(z)
//	where NUM and DEN contain the polynomial coefficients in 
//	descending powers of z.
//
//	See also: dcgain
//
//--------------------------------------------------------------------------

require tfchk abcdchk

ddcgain = function(a,b,c,d)
{
  if (nargs==2) {
    // Transfer function discription
    </den;num/> = tfchk(a,b);
    if (length(num) && length(den)) {
      k = sum(num')'/sum(den);
    else
      k = [];
    }

  else if (nargs==4) {
    // State space description
    msg = abcdchk(a,b,c,d);
    if (msg != "") { error(msg); }
    k = c/(eye(a.nr,a.nr)-a)*b + d;

  else
    error("Wrong number of input arguments.");
  }}
  return k;  
};

