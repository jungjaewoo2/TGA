//----------------------------------------------------------------------
//
// dmodred
//
// syntax:  </ab;bb;cb;db/> = dmodred(a,b,c,d,elim)
//
// Discrete-time model state reduction.
//
// </Ab;Bb;Cb;Db/> = dmodred(A,B,C,D,elim) reduces the order of a model
// by eliminating the states specified in vector elim.  The state
// vector is partioned into X1, to be kept, and X2, to be eliminated,
//
//	A = |A11  A12|		B = |B1|	C = |C1 C2|
//	    |A21  A22|		    |B2|
//		
//	x[n+1] = Ax[n] + Bu[n],   y[n] = Cx[n] + Du[n]
//
// X2[n+1] is set to X2[n], and the resulting equations solved for
// X1.  The resulting system has length(elim) fewer states and can be
// envisioned as having set the elim states to be infinitely fast.
//
// See also dbalreal, balreal and modred.
//
//----------------------------------------------------------------------
require abcdchk

dmodred = function(a,b,c,d,elim)
{
  msg = abcdchk(a,b,c,d);
  if (msg != "") { error(msg); }
  
  // Form keep vector:
  ns = b.nr; 
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
  mk  = a11.nr; 
  a22 = a22 - eye(ns-mk,ns-mk);
  ab  = a11 - a12/a22*a21;
  bb  = b1 - a12/a22*b2;
  cb  = c1 - c2/a22*a21;
  db  = d - c2/a22*b2;
  
  return <<ab=ab;bb=bb;cb=cb;db=db>>;
};

