//-------------------------------------------------------------------//

//  Synopsis:   Find the characteristic polynomial

//  Syntax:     charpoly(A)

//  Description:

//  Find characteristic polynomial of A
//
//  If A is an n by n matrix, poly(A) is an n+1 element row vector
//  whose elements are the coefficients of the characteristic polynomial
//      det(A-zI) = 0.
//  The output coefficients are stored in decending powers as
//      [C1,C2,....,Cn+1]
//  It represents polynomial
//      C[1]*z^n + C[2]*z^(n-1) + .... + C[n]*z + C[n+1]  
//
//  If A is an n by 1 column vector, elements of A are roots of 
//  polynomial
//      (z-A[1])(z-A[2]) ... (z-A[n]) = 0
//  The output coefficients represents the expanded form of above
//  polynomial.

//  Example:
//  > X=[5,-6;1,0]
//   X =
//          5         -6
//          1          0
//  > p=charpoly(X)
//   p =
//          1         -5          6
//  > r=roots(p)
//   r =
//          2
//          3
//  > p1=charpoly(r)
//   p1 =
//          1         -5          6
//

//  See Also: roots
//
//  Tzong-Shuoh Yang (tsyang@ce.berkeley.edu)  5/7/94
//
//-------------------------------------------------------------------//

charpoly = function(A)
{ 
  if (A.nc == A.nr) 
  {
    // A is a matrix, find its eigenvalues
    e = eig(A).val';
    l = inf();
    // trim INFs
    e = e[find(e != l && e != -l)];
    else if (A.nr == 1 || A.nc == 1) {
      // elements of A are roots of the polynomial
      e = A[:];	// force a column vector
    else
      error("Argument of poly() must be a square matrix or a vector");
  } }
  // construct char. polynomial
  c = [1,zeros(1,e.nr)];
  for (i in 1:e.nr) 
  {
    j = i + 1;
    c[2:j] = c[2:j] - e[i]*c[1:i];
  }
  
  if (all(imag(c) == 0)) 
  {
    c = real(c);
  }
  return c;
};
