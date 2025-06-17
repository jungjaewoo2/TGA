//----------------------------------------------------------------------
//
// modred
//
// syntax: </ab;bb;cb;db/> = modred(a,b,c,d,elim)
//
// Model state reduction.
// </Ab;Bb;Cb;Db/> = modred(A,B,C,D,ELIM) reduces the order of a model
// by eliminating the states specified in vector ELIM.  The state
// vector is partioned into X1, to be kept, and X2, to be eliminated,
//
//	A = |A11  A12|		B = |B1|	C = |C1 C2|
//	    |A21  A22|		    |B2|
//	.
//	x = Ax + Bu,   y = Cx + Du
//
// The derivative of X2 is set to zero, and the resulting equations
// solved for X1.  The resulting system has length(ELIM) fewer states
// and can be envisioned as having set the ELIM states to be 
// infinitely fast.
//
// See also balreal and dmodred.
//
//----------------------------------------------------------------------
require abcdchk

modred=function(a,b,c,d,elim)
{
  msg = abcdchk(a,b,c,d);
  if (msg != "") { error(msg); }
  
  // Form keep vector:
  ns = b.nr; nu = b.nc;
  keep = 1:ns;
  keep[elim] = [];

  // Partition into x1, to be kept, and x2, to be eliminated:
  a11 = a[keep;keep];
  a12 = a[keep;elim];
  a21 = a[elim;keep];
  a22 = a[elim;elim];
  b1  = b[keep;];
  b2  = b[elim;];
  c1  = c[;keep];
  c2  = c[;elim];

  // Form final reduced matrices
  ab  = a11 - a12/a22*a21;
  bb  = b1 - a12/a22*b2;
  cb  = c1 - c2/a22*a21;
  db  = d - c2/a22*b2;

  return <<ab=ab;bb=bb;cb=cb;db=db>>;
};

