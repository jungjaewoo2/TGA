//-------------------------------------------------------------------//

//  Synopsis:   Create a random sparse matrix.

//  Syntax:	sprand ( M , N , DENSITY )

//  Description:

//  Create a randomly populated and valued matrix. M and N specify the
//  matrix row and column sizes. DENSITY ( 0 < DENSITY < 1 ) specifies
//  the desired ratio of non-zero to zero elements. The location of
//  the non-zero elements, and their values are uniformly distributed
//  random numbers.

//-------------------------------------------------------------------//

sprand = function ( m, n, p )
{
  if (!exist (p)) { p = 0.1; }
  ntotal = m * n;
  nnz = max (p * ntotal, 1);

  #
  # Generate the row indices
  #

  # Set up the random number generator
  rng(1, "uniform", 1, m);
  row = int (rand (nnz, 1));

  #
  # Generate the column indices
  #

  # Set up the random number generator
  rng(1,"uniform", 1, n);
  col = int (rand (nnz, 1));

  #
  # Generate the elements.
  #
  rng (1, "uniform", 1, m*n);
  el = rand(nnz, 1);

  tmp = [row, col, el; m, n, 0];

  return spconvert (tmp);
};
