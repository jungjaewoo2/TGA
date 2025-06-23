//-------------------------------------------------------------------//

// Synopsis:	Randomly sets matrix elements to zero.

// Syntax:	S = sparsify ( A , P )

// Description:

//	S is A with elements randomly set to zero (S = S' if A is
//	square and A = A', i.e. symmetry is preserved). Each element
//	has probability P of being zeroed. Thus on average 100*P
//	percent of the elements of A will be zeroed. 

//	Default: P = 0.25. 

//	This file is a translation of sparsify.m from version 2.0 of
//	"The Test Matrix Toolbox for Matlab", described in Numerical
//	Analysis Report No. 237, December 1993, by N. J. Higham.

//-------------------------------------------------------------------//

sparsify = function ( A , p )
{
  if (!exist (p)) { p = 0.25; }
  if (p < 0 || p > 1) {
    error("Second parameter must be between 0 and 1 inclusive.");
  }

  // Is A square and symmetric?
  symm = 0;
  if (min(size(A)) == max(size(A))) {
    if (norm(A-A',"1") == 0) { symm = 1; }
  }

  if (!symm)
  {
    A = A .* (rand(A) > p);		// Unsymmetric case
  else
    A = triu(A,1) .* (rand(A) > p);	// Preserve symmetry
    A = A + A';
    A = A + diag( diag(A) .* (rand(diag(A)) > p) );
  }

  return A;
};
