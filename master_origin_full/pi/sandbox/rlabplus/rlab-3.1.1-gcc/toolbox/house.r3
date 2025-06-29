//-------------------------------------------------------------------//

// Synopsis:	Create a Householder matrix.

// Syntax:	house ( X )

//	If HOUSE(x), which returns a list containing elements `v' and
//	`beta',  then H = EYE - beta*v*v' is a Householder matrix such
//	that Hx = -sign(x(1))*norm(x)*e_1. 
//	NB: If x = 0 then v = 0, beta = 1 is returned.
//            x can be real or complex.
//            sign(x) := exp(i*arg(x)) ( = x./abs(x) when x ~= 0).

//	Theory: (textbook references Golub & Van Loan 1989, 38-43;
//	         Stewart 1973, 231-234, 262; Wilkinson 1965, 48-50).
//      Hx = y: (I - beta*v*v')x = -s*e_1.
//      Must have |s| = norm(x), v = x+s*e_1, and
//      x'y = x'Hx =(x'Hx)' real => arg(s) = arg(x(1)).
//      So take s = sign(x(1))*norm(x) (which avoids cancellation).
//      v'v = (x(1)+s)^2 + x(2)^2 + ... + x(n)^2
//          = 2*norm(x)*(norm(x) + |x(1)|).

//	References:
//        G.H. Golub and C.F. Van Loan, Matrix Computations, second edition,
//           Johns Hopkins University Press, Baltimore, Maryland, 1989.
//        G.W. Stewart, Introduction to Matrix Computations, Academic Press,
//           New York, 1973,
//        J.H. Wilkinson, The Algebraic Eigenvalue Problem, Oxford University
//           Press, 1965.

//	This file is a translation of house.m from version 2.0 of
//	"The Test Matrix Toolbox for Matlab", described in Numerical
//	Analysis Report No. 237, December 1993, by N. J. Higham.

//-------------------------------------------------------------------//

house = function ( x )
{
  m = x.nr; n = x.nc;
  if (n > 1) { error ("Argument must be a column vector."); }

  s = norm(x,"2") * (sign(x[1]) + (x[1]==0)); // Modification for sign(0)=1.
  v = x;
  if (s == 0)	// Quit if x is the zero vector.
  {
    beta = 1; 
    return << beta = beta ; v = v >>;
  }

  v[1] = v[1] + s;
  beta = 1/(s'*v[1]);		// NB the conjugated s.

  // beta = 1/(abs(s)*(abs(s)+abs(x(1)) would guarantee beta real.
  // But beta as above can be non-real (due to rounding) only 
  // when x is complex.

  return << beta = beta ; v = v >>;
};
