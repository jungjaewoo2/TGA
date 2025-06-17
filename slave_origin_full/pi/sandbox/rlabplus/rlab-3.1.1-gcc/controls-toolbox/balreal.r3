//--------------------------------------------------------------------------
//
// balreal
//
// Syntax: </ab;bb;cb;g;t/> = balreal(a,b,c)
//
//
//	</Ab;Bb;Cb;G;T/> = balreal(A,B,C) also returns a vector G containing
//	the diagonal of the gramian of the balanced realization, and
//	matrix T, the similarity transformation used to convert (A,B,C)
//	to (Ab,Bb,Cb).  If the system (A,B,C) is normalized properly,
//	small elements in gramian G indicate states that can be removed to
//	reduce the model to lower order.
//
//	See:
//	 1) Moore, B., Principal Component Analysis in Linear Systems:
//	    Controllability, Observability, and Model Reduction, IEEE
//	    Transactions on Automatic Control, 26-1, Feb. 1981.,
//	 2) Laub, A., "Computation of Balancing Transformations", Proc. JACC
//	    Vol.1, paper FA8-E, 1980.
//--------------------------------------------------------------------------
require gram

balreal = function(a,b,c)
{
  gc = gram(a,b);
  go = gram(a',c');
  r = chol(gc);
  rgr = r*go*r';
  rgr = tril(rgr) + tril(rgr,-1)';  // Make rgr exactly symmetric.
  tmp = eig(rgr);
  t = r'*tmp.vec*diag(tmp.val.^(-.25));
  ab = t\a*t;
  bb = t\b;
  cb = c*t;
  g = diag(gram(ab,bb))';

  // Sort so g is in descending order
  i = sort(-g).idx;
  ab = ab[i;i];
  bb = bb[i;];
  cb = cb[;i];
  t = t[i;i];
  g = g[i];

  return <<ab=ab;bb=bb;cb=cb;g=g;t=t>>;
};


