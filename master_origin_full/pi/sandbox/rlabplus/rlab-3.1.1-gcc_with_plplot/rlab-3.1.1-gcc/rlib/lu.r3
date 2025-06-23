//-------------------------------------------------------------------//

//  Synopsis:   Compute the L and U factors of the matrix A.

//  Syntax:	lu ( A )

//  Description:

//  Lu computes the "LU" factors of the input MATRIX. The input MATRIX
//  must be square. lu() returns l, u, and pvt in a LIST. The
//  factorization has the form: 

//  A = p*l*u			where A is the input matrix

//  The function lu is a user-function that utilizes the built-in
//  function factor. Factor does the factorization using the LAPACK
//  functions DGETRF/ZGETRF, and DGECON/ZGECON.

//  Example:

//  > a = [1,2,3;4,5,6;7,8,0]
//   a =
//   matrix columns 1 thru 3
//             1           2           3
//             4           5           6
//             7           8           0
//  > f = lu(a)
//    l            pvt            u            
//  > f.p
//   matrix columns 1 thru 3
//             0           0           1
//             1           0           0
//             0           1           0
//  > f.l
//   matrix columns 1 thru 3
//             1           0           0
//         0.143           1           0
//         0.571         0.5           1
//  > f.u
//   matrix columns 1 thru 3
//             7           8           0
//             0       0.857           3
//             0           0         4.5
//  > f.p*f.l*f.u
//   matrix columns 1 thru 3
//             1           2           3
//             4           5           6
//             7           8           0

//  See Also: backsub, factor, inv, solve

//-------------------------------------------------------------------//

static (swap);
require tril triu norm eye

lu = function ( A )
{
  if (A.nr != A.nc) { error ("lu() requires square A"); }

  x = factor (A, "g");	// Do the factorization

  //
  // Now create l, u, and pvt from lu and pvt.
  //

  l = tril (x.lu, -1) + eye (size (x.lu));
  u = triu (x.lu);
  pvt = eye (size (x.lu));

  //
  // Now re-arange the columns of pvt
  //

  for (i in 1:max (size (x.lu)))
  {
    pvt = pvt[ ; swap (1:pvt.nc, i, x.pvt[i]) ];
  }
  return << l = l; u = u; pvt = pvt >>;
};

//
//  In vector V, swap elements I, J
//

swap = function ( V, I, J )
{
  local (v, tmp);
  v = V;
  tmp = v[I];
  v[I] = v[J];
  v[J] = tmp;
  return v;
};
