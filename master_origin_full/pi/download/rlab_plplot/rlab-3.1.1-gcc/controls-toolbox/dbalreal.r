//--------------------------------------------------------------------------
//
// dbalreal
//
// Syntax: </ab;bb;cb;m;T/> = dbalreal(a,b,c)
//
//	Discrete balanced state-space realization and model reduction.
//
//	</ab;bb;cb;m;T/> = dbalreal(A,B,C) returns a balanced state-space 
//	realization of the system (A,B,C). It also returns a vector m 
//	containing the diagonal of the gramian of the balanced realization
//	and matrix T, the similarity transformation used to convert 
//	(A,B,C) to (ab,bb,cb).  If the system (A,B,C) is normalized 
//	properly, small elements in gramian M indicate states that can be
//	removed to reduce the model to lower order.
//
//	See also dmodred, balreal and modred.
//
//--------------------------------------------------------------------------

require dgram

dbalreal = function(a,b,c)
{

  Gc = dgram(a,b);
  Go = dgram(a',c');
  R = chol(Gc);
  tmp = svd(R*Go*R');
  u = tmp.u;
  D = diag(tmp.sigma);
  V = tmp.vt';
  D = D.*sign(u'*V);
  T = R'*V*diag(diag(D).^(-.25));
  ab = T\a*T;
  bb = T\b;
  cb = c*T;
  m = diag(dgram(ab,bb))';

  return <<ab=ab;bb=bb;cb=cb;m=m;T=T>>;
};

