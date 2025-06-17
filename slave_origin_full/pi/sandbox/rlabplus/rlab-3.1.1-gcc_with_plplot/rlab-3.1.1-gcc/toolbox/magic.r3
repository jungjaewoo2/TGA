//-------------------------------------------------------------------//

// Synopsis:	Create magic square matrix

// Syntax:      magic ( N )

// Description:

//      Given the matrix order, N, return the magic square
//      matrix. This matrix is called magic because of its special
//      properties:

//         The sum of each column is the same.
//         The sum of each row is the same.
//         The sum of elements on the diagonal is the same.

//      The elements of each matrix have values 1 through N^2.

//      Found this algorithm in a Fortran routine on an old VAX at the
//      University?  Don't know who the author is though?

//-------------------------------------------------------------------//

static (swap)

magic = function ( N )
{
  a = zeros (N, N);
   
  if (int(mod(N,4)) != 0)
  {
    if (int(mod(N,2)) == 0) { m = N/2; }
    if (int(mod(N,2)) != 0) { m = N; }

    # Odd order or upper corner of even order
  
    for (j in 1:m)
    {
      for (i in 1:m)
      {
        a[i;j] = 0;
      }
    }
  
    i = 1;
    j = int((m+1)/2);
    mm = m*m;

    for (k in 1:mm)
    {
      a[i;j] = k;
      i1 = i - 1;
      j1 = j + 1;
      if (i1 < 1) { i1 = m; }
      if (j1 > m) { j1 = 1; }
      if(int(a[i1;j1]) != 0)
      {
        i1 = i + 1;
        j1 = j;
      }
      i = i1;
      j = j1;
    }
    if (int(mod(N, 2)) != 0) { return a; }
  
    # Rest of even order
  
    t = m*m;
    for (i in 1:m)
    {
      for (j in 1:m)
      {
        im = i + m;
        jm = j + m;
        a[i;jm] = a[i;j] + 2*t;
        a[im;j] = a[i;j] + 3*t;
        a[im;jm] = a[i;j] + t;
      }
    }
  
    m1 = int((m-1)/2);
    if (m1 == 0) { return a; }
    for (j in 1:m1)
    {
      a = swap (a, m, 1, j, m+1, j);
    }
    m1 = int((m + 1)/2);
    m2 = m1 + m;

    a = swap (a, 1, m1, 1, m2, 1);
    a = swap (a, 1, m1, m1, m2, m1);

    m1 = N + 1 - int((m - 3)/2);
    if (m1 > N) { return a; }

    for (j in m1:N)
    {
      a = swap (a, m, 1, j, m+1, j);
    }
    return a;

  } else {

    k = 1;
    for (i in 1:N)
    {
      for (j in 1:N)
      {
        a[i;j] = k;
        if (int(mod(i,4)/2) == int(mod(j,4)/2))
        {
          a[i;j] = N*N + 1 - k;
        }
        k = k + 1;
      }
    }
  }
  
  return a;
};

#
# Swap elements in a matrix in column-major fashion.
# At present their is no error checking.
#

# A:	the input matirx
# N:	number of elements to swap
# R1:	the starting row number
# C1:	the starting column number
# R2:	the starting row number
# C2:	the stating column number
#

swap = function (A, N, R1, C1, R2, C2)
{
  local (i, nr, nc, tmp, s1, s2)

  nr = A.nr; nc = A.nc;
  s1 = (C1-1)*nr + R1;
  s2 = (C2-1)*nr + R2;

  for (i in 0:N-1)
  {
    tmp = A[s1+i];
    A[s1+i] = A[s2+i];
    A[s2+i] = tmp;
  }
  return A;
};
