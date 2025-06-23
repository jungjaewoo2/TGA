//-------------------------------------------------------------------//

// Synopsis:	Band reduction by two-sided unitary transformations.

// Syntax:	bandred ( A , KL, KU )

// Description:

//	bandred(A, KL, KU) is a matrix unitarily equivalent to A with
//	lower bandwidth KL and upper bandwidth KU (i.e. B(i,j) = 0 if
//	i > j+KL or j > i+KU).  The reduction is performed using
//	Householder transformations. If KU is omitted it defaults to
//	KL. 

//	Called by randsvd.
//	This is a `standard' reduction.  Cf. reduction to bidiagonal
//	form prior to computing the SVD.  This code is a little
//	wasteful in that it computes certain elements which are
//	immediately set to zero! 
//
//      Reference:
//      G.H. Golub and C.F. Van Loan, Matrix Computations, second edition,
//      Johns Hopkins University Press, Baltimore, Maryland, 1989.
//      Section 5.4.3.

//	This file is a translation of bandred.m from version 2.0 of
//	"The Test Matrix Toolbox for Matlab", described in Numerical
//	Analysis Report No. 237, December 1993, by N. J. Higham.

// Dependencies
   require house

//-------------------------------------------------------------------//

bandred = function ( A , kl , ku )
{
  if (!exist (ku)) { ku = kl; } else { ku = ku; }

  if (kl == 0 && ku == 0) {
    error ("You''ve asked for a diagonal matrix.  In that case use the SVD!");
  }

  // Check for special case where order of left/right transformations matters.
  // Easiest approach is to work on the transpose, flipping back at the end.

  flip = 0;
  if (ku == 0)
  {
    A = A';
    temp = kl; 
    kl = ku; 
    ku = temp; 
    flip = 1;
  }

  m = A.nr; n = A.nc; 

  for (j in 1 : min( min(m, n), max(m-kl-1, n-ku-1) ))
  {
    if (j+kl+1 <= m)
    {
      </ beta; v /> = house(A[j+kl:m;j]);
      temp = A[j+kl:m;j:n];
      A[j+kl:m;j:n] = temp - beta*v*(v'*temp);
      A[j+kl+1:m;j] = zeros(m-j-kl,1);
    }
    
    if (j+ku+1 <= n)
    {
      </ beta; v /> = house(A[j;j+ku:n]');
      temp = A[j:m;j+ku:n];
      A[j:m;j+ku:n] = temp - beta*(temp*v)*v';
      A[j;j+ku+1:n] = zeros(1,n-j-ku);
    }
  }
  
  if (flip) {
    A = A';
  }

  return A;
};
