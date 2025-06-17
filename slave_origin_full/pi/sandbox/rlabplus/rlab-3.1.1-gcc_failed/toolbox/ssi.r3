//-------------------------------------------------------------------//

// Synopsis:  Compute the N eigenvalues and eigenvectors for the real
//            generalized eigenvalue problem, A*v = B*v*lambda, via
//            sub-space iteration.

// Syntax:    ssi ( A, B, V, N, eps )

// Description:

//      This routine uses subspace iteration to compute the first #
//      few eigenvalues and vectors for a real, symmetric, generalized
//      eigenvalue problem A*v = B*v*lambda. Where v are the
//      eigenvectors, and lambda are the eigenvalues.

//      Both A and B must be real and symmetric.
//      The matrix V contains the starting vectors, or the initial
//      guesses at the desired N eigenvectors.
//      N specifies the number of desired eigenvalues and
//      eigenvectors. If N is not specified, then ssi will calculate N
//      (usually the number of columns in V).
//      The user may also specify the error tolerance, eps. eps
//      defaults to 1.e-6.

//      This method is useful for large sparse eigenvalue problems
//      when only a some of the eigenvalues and vectors are needed.

//      Reference: Batte's finite element book

//      Original ssi() written by Scott Hunziker

//      Ported to Rlab by Ian Searle

// Dependencies

require symm

//-------------------------------------------------------------------//

ssi = function( A, B, V, n, eps )
{
  #
  # Error check arguments, and setup defaults if necessary.
  #

  if (!issymm (A) || !issymm (B))
  {
    error ("ssi: A and B must be symmetric");
  }

  if (A.type != "real" || B.type != "real" || V.type != "real")
  {
    error ("ssi: Arguments must be REAL");
  }

  if (!exist (eps)) { eps = 1.0E-6; }

  if (!exist (n))
  {
    if ( V.nc > 16 )
    {
      n = V.nc - 8;
      } else {
      n = int ((V.nc+1) / 2);
    }
    } else {
    if (n <= 0 || n > V.nc)
    {
      error ("ssi: illegal value for n");
    }
  }

  lambda = (diag(V'*A*V) ./ diag( V'*B*V))';

  n = 1:n;
  while ( 1 )
  {
    Vb = solve( A, B*V );
    e = eig(symm(Vb'*A*Vb), symm(Vb'*B*Vb));
    d = e.val;
    tmp = abs( lambda - d );
    for (i in n) { tmp[i] = tmp[i] ./ d[i]; }
    r = tmp[ maxi( tmp[n] ) ];
    lambda = e.val;
    V = Vb*e.vec;
    if ( r < eps ) { break; }
  }

  return << val=e.val; vec=V >>;
};
