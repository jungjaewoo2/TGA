//
// Do some setup...
//

clearall();

tic();   // Start timing the tests...

pi = atan (1.0)*4;

// should be 3 (heuristic).
X = 3;

// --------------------------------------------------------------------
show = function ( A )
{
  if (!exist (A))
  {
    return ("\tUNDEFINED");
  }

  if (class (A) != "list")
  {
    for (i in members (A))
    {
      printf ("\t%s:\t%s\n", i, A.[i]);
    }
  else
    for (i in members (A))
    {
      printf ("\t%-20s:%s\t%s\n", i, class (A.[i]), type(A.[i]));
    }
  }
};

epsilon = function()
{
  eps = 1.0;
  while((eps + 1.0) != 1.0)
  {
    eps = eps/2.0;
  }
  return 2*eps;
};

colors("red");
eps = 1e-14;  // numerical estimate for round-off error
printf("Note: Error bound was changed from the original %g to %g, as the latter\n", epsilon(), eps);
printf("      is more realistic estimate of round-off error\n");
colors()

// --------------------------------------------------------------------

clearall = function ( )
{
  for (i in members ($$))
  {
    if (class ($$.[i]) != "function")
    {
      if (i != "pi" && i != "eps" && i != "_rlab_search_path")
      {
        clear ($$.[i]);
      }
    }
  }
};

// --------------------------------------------------------------------

norm = function ( A , ntype )
{

  if (exist (ntype))
  {
    if (class(ntype) == "string")
    {
      // This has to be done with the matrix norm function.
      return mnorm ( A , ntype );
    }
  }

  //
  // Check for matrix input
  //

  if (min (size (A)) != 1)
  {
    //
    // Two possibilties...
    //

    if (!exist (ntype)) {
      return mnorm (A);
    }

    if (class (ntype) == "string" || class (ntype) == "num") {
      return mnorm (A, ntype);
    }

  else

    //
    // Vector input...
    //

    if (!exist (ntype)) {
      error ("norm: second argument required");
    }

    if (class (ntype) != "num") {
      error ("norm: second argument must be numeric");
    }

    return vpnorm (A, ntype);
  }
};

// --------------------------------------------------------------------

trace = function(m)
{
  if(m.class != "num")
  {
    error("must provide NUMERICAL input to trace()");
  }

  tr = 0;
  for(i in 1:min( [m.nr, m.nc] ))
  {
    tr = tr + m[i;i];
  }

  return tr;
};

// --------------------------------------------------------------------

symm = function( A )
{
  return (A + A')./2;
};

// --------------------------------------------------------------------

triu = function(x, k)
{
  if (!exist (k)) { k = 0; }
  nr = x.nr; nc = x.nc;

  if(k > 0)
  {
    if (k > (nc - 1)) { error ("triu: invalid value for k"); }
  else
    if (abs (k) > (nr - 1)) { error ("triu: invalid value for k"); }
  }

  y = zeros(nr, nc);

  for(j in max( [1,1+k] ):nc)
  {
    i = 1:min( [nr, j-k] );
    y[i;j] = x[i;j];
  }

  return y;
};

// --------------------------------------------------------------------

//
// Magic square
// magic() returns a square matrix whose elements
// are 1 through N^2
// Magic squares have equal row and column sums.
//
// Found this algorithm in a Fortran routine
// on an old VAX at the University? Seems to
// work quite well. Don't know who the original author is though?
//

magic = function ( N )
{
  a = zeros (N, N);

  if (int(mod(N,4)) != 0)
  {
    if (int(mod(N,2)) == 0) { m = N/2; }
    if (int(mod(N,2)) != 0) { m = N; }

    // Odd order or upper corner of even order

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

    // Rest of even order

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
      swap (a, m, 1, j, m+1, j);
    }
    m1 = int((m + 1)/2);
    m2 = m1 + m;

    swap (a, 1, m1, 1, m2, 1);
    swap (a, 1, m1, m1, m2, m1);

    m1 = N + 1 - int((m - 3)/2);
    if (m1 > N) { return a; }

    for (j in m1:N)
    {
      swap (a, m, 1, j, m+1, j);
    }
    return a;

  else

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

// --------------------------------------------------------------------

//
// Swap elements in a matrix in column-major fashion.
// At present their is no error checking.
//

// A:	the input matirx
// N:	number of elements to swap
// R1:	the starting row number
// C1:	the starting column number
// R2:	the starting row number
// C2:	the stating column number
//

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
  return 1;
};

// --------------------------------------------------------------------

tril = function(x, k)
{
  if (!exist (k)) { k = 0; }
  nr = x.nr; nc = x.nc;
  if(k > 0)
  {
    if (k > (nc - 1)) { error ("tril: invalid value for k"); }
  else
    if (abs (k) > (nr - 1)) { error ("tril: invalid value for k"); }
  }

  y = zeros(nr, nc);

  for(i in max( [1,1-k] ):nr)
  {
    j = 1:min( [nc, i+k] );
    y[i;j] = x[i;j];
  }

  return y;
};

// --------------------------------------------------------------------

eye = function( m , n )
{
  local(m, n)

  if (!exist (n))
  {
    if(m.n != 2) { error("only 2-el MATRIX allowed as eye() arg"); }
    new = zeros (m[1], m[2]);
    N = min ([m[1], m[2]]);
  else
    if (class (m) == "string" || class (n) == "string") {
      error ("eye(), string arguments not allowed");
    }
    if (max (size (m)) == 1 && max (size (n)) == 1)
    {
      new = zeros (m[1], n[1]);
      N = min ([m[1], n[1]]);
    else
      error ("matrix arguments to eye() must be 1x1");
    }
  }
  for(i in 1:N)
  {
    new[i;i] = 1.0;
  }
  return new;
};

// --------------------------------------------------------------------

speye = function ( m, n )
{
  //
  // Emulate Matlab's interface...
  //

  if (nargs == 1)
  {
    if (length (m) == 1)
    {
      n = m;
    else if (length (m) == 2) {
      n = m[2];
      m = m[1];
    else
      rerror ("speye: incorrect usage");
  }}}

  //
  // Create the matrix triplet form of the identity matrix.
  //

  s = min (m, n);
  k = (1:round (s))';
  K = [k, k, ones (s, 1)];
  if (m != n)
  {
    K = [K; m, n, 0];
  }

  //
  // Now convert it into compressed row-wise storage.
  //

  return spconvert (K);
};

// --------------------------------------------------------------------

sprand = function ( m, n, p )
{
  if (!exist (p)) { p = 0.1; }
  ntotal = m * n;
  nnz = max (p * ntotal, 1);

  //
  // Generate the row indices
  //

  // Set up the random number generator
  row = int (1 + m*urandom (nnz, 1));

  //
  // Generate the column indices
  //

  // Set up the random number generator
  col = int (1+ n * urandom (nnz, 1));

  //
  // Generate the elements.
  //
  el = 1 + m * n * urandom(nnz, 1);

  tmp = [row, col, el; m, n, 0];

  return spconvert (tmp);
};

// --------------------------------------------------------------------

ohess = function ( x )
{
  global (pi)

  if (any (imag (x))) { error("Parameter must be real."); }

  n = max(size(x));

  if (n == 1)
  {
    //  Handle scalar x.
    n = x;
    x = urandom(n-1, 1)*2*pi;
    H = eye(n,n);
    H[n;n] = sign(gaussian());
  else
    H = eye(n,n);
    H[n;n] = sign(x[n]) + (x[n]==0);   // Second term ensures H[n;n] nonzero.
  }

  for (i in n:2:-1)
  {
    // Apply Givens rotation through angle x(i-1).
    theta = x[i-1];
    c = cos(theta);
    s = sin(theta);
    H[ [i-1, i] ;] = [ c*H[i-1;]+s*H[i;] ;
                       -s*H[i-1;]+c*H[i;] ];
  }

  return H;
};

// --------------------------------------------------------------------
qmult = function ( A )
{
  local (A)

  n = A.nr; m = A.nc;

  //  Handle scalar A.
  if (max(n,m) == 1)
  {
    n = A;
    A = eye(n,n);
  }

  d = zeros(n,n);

  for (k in n-1:1:-1)
  {
    // Generate random Householder transformation.
    x = gaussian(n-k+1,1);
    s = norm(x, "2");
    sgn = sign(x[1]) + (x[1]==0);	// Modification for sign(1)=1.
    s = sgn*s;
    d[k] = -sgn;
    x[1] = x[1] + s;
    beta = s*x[1];

    // Apply the transformation to A.
    y = x'*A[k:n;];
    A[k:n;] = A[k:n;] - x*(y/beta);
  }

  // Tidy up signs.
  for (i in 1:n-1)
  {
    A[i;] = d[i]*A[i;];
  }

  A[n;] = A[n;]*sign(gaussian());
  B = A;

  return B;
};

// --------------------------------------------------------------------

house = function ( x )
{
  local (x)

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

// --------------------------------------------------------------------

bandred = function ( A , kl , ku )
{
  local (A, kl, ku)

  if (!exist (ku)) { ku = kl; else ku = ku; }

  if (kl == 0 && ku == 0) {
    error ("You''ve asked for a diagonal matrix.  In that case use the SVD!");
  }

  // Check for special case where order of left/right transformations matters.
  // Easiest approach is to work on the transpose, flipping back at the end.

  flip = 0;
  if (ku == 0)
  {
    A = A';
    temp = kl; kl = ku; ku = temp; flip = 1;
  }

  m = A.nr; n = A.nc;

  for (j in 1 : min( min(m, n), max(m-kl-1, n-ku-1) ))
  {
    if (j+kl+1 <= m)
    {
       </beta; v/> = house(A[j+kl:m;j]);
       temp = A[j+kl:m;j:n];
       A[j+kl:m;j:n] = temp - beta*v*(v'*temp);
       A[j+kl+1:m;j] = zeros(m-j-kl,1);
    }

    if (j+ku+1 <= n)
    {
       </beta; v/> = house(A[j;j+ku:n]');
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

// --------------------------------------------------------------------

randsvd = function (n, kappa, mode, kl, ku)
{
  global (eps)

  if (!exist (kappa)) { kappa = sqrt(1/eps); }
  if (!exist (mode))  { mode = 3; }
  if (!exist (kl))    { kl = n-1; }	// Full matrix.
  if (!exist (ku))    { ku = kl; } 	// Same upper and lower bandwidths.

  if (abs(kappa) < 1) {
    error("Condition number must be at least 1!");
  }

  posdef = 0;
  if (kappa < 0)			// Special case.
  {
    posdef = 1;
    kappa = -kappa;
  }

  p = min(n);
  m = n[1];		// Parameter n specifies dimension: m-by-n.
  n = n[max(size(n))];

  if (p == 1)		// Handle case where A is a vector.
  {
    A = gaussian(m, n);
    A = A/norm(A, "2");
    return A;
  }

  j = abs(mode);

  // Set up vector sigma of singular values.
  if (j == 3)
  {
    factor = kappa^(-1/(p-1));
    sigma = factor.^[0:p-1];

  else if (j == 4) {
    sigma = ones(p,1) - (0:p-1)'/(p-1)*(1-1/kappa);

  else if (j == 5) {	// In this case cond(A) <= kappa.
    sigma = exp( -urandom(p,1)*log(kappa) );

  else if (j == 2) {
    sigma = ones(p,1);
    sigma[p] = 1/kappa;

  else if (j == 1) {
    sigma = ones(p,1)./kappa;
    sigma[1] = 1;

  }}}}}


  // Convert to diagonal matrix of singular values.
  if (mode < 0) {
    sigma = sigma[p:1:-1];
  }

  sigma = diag(sigma);

  if (posdef)		// Handle special case.
  {
    Q = qmult(p);
    A = Q'*sigma*Q;
    A = (A + A')/2;	// Ensure matrix is symmetric.
    return A;
  }

  if (m != n)
  {
    sigma[m; n] = 0;	// Expand to m-by-n diagonal matrix.
  }

  if (kl == 0 && ku == 0)	// Diagonal matrix requested - nothing more to do.
  {
    A = sigma;
    return A;
  }

  // A = U*sigma*V, where U, V are random orthogonal matrices from the
  // Haar distribution.
  A = qmult(sigma');
  A = qmult(A');

  if (kl < n-1 || ku < n-1)	// Bandwidth reduction.
  {
   A = bandred(A, kl, ku);
  }

  return A;
};

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
    A = A .* (urandom(A) > p);		// Unsymmetric case
  else
    A = triu(A,1) .* (urandom(A) > p);	// Preserve symmetry
    A = A + A';
    A = A + diag( diag(A) .* (urandom(diag(A)) > p) );
  }

  return A;
};

//-------------------------------------------------------------------//

// Synopsis:    Lehmer matrix - symmetric positive definite.

// Syntax:      A = lehmer ( N )

// Description:

//      A is the symmetric positive definite N-by-N matrix with
//                     A(i,j) = i/j for j >= i.
//      A is totally nonnegative.  INV(A) is tridiagonal, and explicit
//      formulas are known for its entries.

//      N <= COND(A) <= 4*N*N.

//      References:
//        M. Newman and J. Todd, The evaluation of matrix inversion
//           programs, J. Soc. Indust. Appl. Math., 6 (1958), pp. 466-476.
//        Solutions to problem E710 (proposed by D.H. Lehmer): The inverse
//           of a matrix, Amer. Math. Monthly, 53 (1946), pp. 534-535.
//        J. Todd, Basic Numerical Mathematics, Vol. 2: Numerical Algebra,
//           Birkhauser, Basel, and Academic Press, New York, 1977, p. 154.

//	This file is a translation of lehmer.m from version 2.0 of
//	"The Test Matrix Toolbox for Matlab", described in Numerical
//	Analysis Report No. 237, December 1993, by N. J. Higham.

//-------------------------------------------------------------------//

lehmer = function ( n )
{
  local (n)
  global (tril)

  A = ones(n,1)*(1:n);
  A = A./A';
  A = tril(A) + tril(A,-1)';

  return A;
};

//
// -------------------- Test control flow -----------------------------
//

for (i in 1:1000)
{
  if (i == 33)
  {
    break;
  }
}

if (i != 33) { error(); }

i = 0;
while (i < 100)
{
  if (i == 34)
  {
    break;
  }
  i = i + 1;
}

if (i != 34) { error(); }

//
// Nested BREAK...
//

for (i in 1:1000)
{
  for (j in 1:1000)
  {
    if (i == 3 && j == 100)
    {
      break;
    }
  }
  if (i == 3 && j == 100)
  {
    break;
  }
}

if (i != 3 || j != 100) { error (); }

printf("\tpassed BREAK test...\n");

//
// CONTINUE Test...
//

j = 1;
for (i in 1:10)
{
  if (i == 3 || i == 7)
  {
    continue;
  }
  j = j + 1;
}

if (j != 9) { error (); }

j = 1; i = 1;
while (i < 10)
{
  if (i == 3 || i == 5 || i == 8)
  {
    i = i + 1;
    continue;
  }
  i = i + 1;
  j = j + 1;
}

if (j != 7) { error (); }

printf("\tpassed CONTINUE test...\n");

//
// -------------------- Test Relational Expressions -------------------
//

//    SCALAR CONSTANTS (REAL)

if( !(1<2) ) { error(); }
if( !(1<=2) ) { error(); }
if( 1>2 ) { error(); }
if( 1>=2 ) { error(); }
if( 1==2 ) { error(); }
if( !(1!=2) ) { error(); }
if( !1 ) { error(); }
if( !!!1) { error(); }

if( !([1]<[2]) ) { error(); }
if( !([1]<=[2]) ) { error(); }
if( [1]>[2] ) { error(); }
if( [1]>=[2] ) { error(); }
if( [1]==[2] ) { error(); }
if( !([1]!=[2]) ) { error(); }
if( ![1] ) { error(); }
if( !!![1]) { error(); }

printf ("\tfinished real-scalar constants...\n");

//    SCALAR CONSTANTS (COMPLEX)

if( !(1+2i<2+3i) ) { error(); }
if( !(1+2i<=2+3i) ) { error(); }
if( 1+2i>2+3i ) { error(); }
if( 1+2i>=2+3i ) { error(); }
if( 1+2i==2+3i ) { error(); }
if( !(1+2i!=2+3i) ) { error(); }
if( !1+2i ) { error(); }
if( !!!1+2i) { error(); }

if( !([1+2i]<[2+3i]) ) { error(); }
if( !([1+2i]<=[2+3i]) ) { error(); }
if( [1+2i]>[2+3i] ) { error(); }
if( [1+2i]>=[2+3i] ) { error(); }
if( [1+2i]==[2+3i] ) { error(); }
if( !([1+2i]!=[2+3i]) ) { error(); }
if( ![1+2i] ) { error(); }
if( !!![1+2i]) { error(); }

printf ("\tfinished complex-scalar constants...\n");

//    SCALAR ENTITIES (REAL)
a=1;b=2;
if( !(a<b) ) { error(); }
if( !(a<=b) ) { error(); }
if( a>b ) { error(); }
if( a>=b ) { error(); }
if( a==b ) { error(); }
if( !(a!=b) ) { error(); }
if( !a ) { error(); }
if( !!!a) { error(); }

if( !([a]<[b]) ) { error(); }
if( !([a]<=[b]) ) { error(); }
if( [a]>[b] ) { error(); }
if( [a]>=[b] ) { error(); }
if( [a]==[b] ) { error(); }
if( !([a]!=[b]) ) { error(); }
if( ![a] ) { error(); }
if( !!![a]) { error(); }

"\tfinished real-scalar entities..."

//    SCALAR ENTITIES (COMPLEX)
a=1+2i;b=2+3i;
if( !(a<b) ) { error(); }
if( !(a<=b) ) { error(); }
if( a>b ) { error(); }
if( a>=b ) { error(); }
if( a==b ) { error(); }
if( !(a!=b) ) { error(); }
if( !a ) { error(); }
if( !!!a) { error(); }

if( !([a]<[b]) ) { error(); }
if( !([a]<=[b]) ) { error(); }
if( [a]>[b] ) { error(); }
if( [a]>=[b] ) { error(); }
if( [a]==[b] ) { error(); }
if( !([a]!=[b]) ) { error(); }
if( ![a] ) { error(); }
if( !!![a]) { error(); }

printf ("\tfinished complex-scalar entities...\n");

//
// ------------------------- Test Misc Stuff ----- ---------------------
//

if(any(any(gaussian(3,3).^0 != ones(3,3)))) { error(); }
if(any(any(urandom(3,3).^0 != ones(3,3)))) { error(); }
if(any(any(gaussian(1,3).^0 != ones(1,3)))) { error(); }
if(any(any(urandom(3,1).^0 != ones(3,1)))) { error(); }
if(any(any(1.^zeros(3,1) != ones(3,1)))) { error(); }
if(any(any(1.^zeros(1,3) != ones(1,3)))) { error(); }

//
// --------------------- Test Misc Matrix Stuff  -----------------------
//

if (any (any ((x[;2] = [1;2;3]) != [0,1;0,2;0,3])))
{ error (); }
if (any (any ((x[;2] = [1+2i;2+3i;3+4i]) != [0,1+2i;0,2+3i;0,3+4i])))
{ error (); }

clear(s);
if (any (any ((s[;2] = ["a";"b";"c"]) != ["","a";"","b";"","c"]))) { error (); }

//
// --------------------- Test Matrix Multiply --------------------------
//

a = [1,2,3;4,5,6;7,8,9];
b = [4,5,6;7,8,9;10,11,12];
c = [ 48,  54,  60 ;
     111, 126, 141 ;
     174, 198, 222 ];

if (any (any (c != a*b))) { error ("failed Real-Real Multiply"); }

az = a + b*1i;
bz = b + a*1i;

cz = [-18+141*1i , -27+162*1i , -36+183*1i ;
        9+240*1i ,   0+279*1i ,  -9+318*1i ;
       36+339*1i ,  27+396*1i ,  18+453*1i ];

czz = [ 48+30*1i ,  54+36*1i  ,  60+42*1i ;
       111+66*1i , 126+81*1i  , 141+96*1i ;
       174+102*1i, 198+126*1i , 222+150*1i ];

czzz = [ 48+111*1i ,  54+126*1i ,  60+141*1i ;
        111+174*1i , 126+198*1i , 141+222*1i ;
        174+237*1i , 198+270*1i , 222+303*1i ];

if (any (any (cz != az*bz)))  { error ("failed Complex-Complex Multiply"); }
if (any (any (czz != a*bz)))  { error ("failed Real-Complex Multiply"); }
if (any (any (czzz != az*b))) { error ("failed Complex-Real Multiply"); }

a = [a,a];
b = [b;b];
c = [  96 , 108 , 120 ;
      222 , 252 , 282 ;
      348 , 396 , 444 ];

if (any (any (c != a*b))) { error ("failed Real-Real Multiply"); }

az = [az,az];
bz = [bz;bz];

cz = [  -36+282*1i ,  -54+324*1i ,  -72+366*1i ;
         18+480*1i ,    0+558*1i ,  -18+636*1i ;
         72+678*1i ,   54+792*1i ,   36+906*1i ];

czz = [  96+60*1i  , 108+72*1i  , 120+84*1i  ;
        222+132*1i , 252+162*1i , 282+192*1i ;
        348+204*1i , 396+252*1i , 444+300*1i ];

czzz = [  96+222*1i , 108+252*1i , 120+282*1i ;
         222+348*1i , 252+396*1i , 282+444*1i ;
         348+474*1i , 396+540*1i , 444+606*1i ];

if (any (any (cz != az*bz)))  { error ("failed Complex-Complex Multiply"); }
if (any (any (czz != a*bz)))  { error ("failed Real-Complex Multiply"); }
if (any (any (czzz != az*b))) { error ("failed Complex-Real Multiply"); }

printf ("\tfinished matrix multiply tests...\n");


//
// --------------------- Test Empty Matrix Behavior  ----------------
//

a = magic (5);

aa = a;
aa[1000,3,400,3,3,1;] =[];
b = a[2,4,5;];
if (all (all (aa != b))) { error ("error in empty matrix behavior"); }

aa = a;
aa[;100,3,1,1,3] = [];
c = a[;2,4,5];
if (all (all (aa != c))) { error ("error in empty matrix behavior"); }

aa = a + a'*1j;
b = aa[2,4,5;];
aa[1000,3,400,3,3,1;] =[];
if (all (all (aa != b))) { error ("error in empty matrix behavior"); }

aa = a + a'*1j;
c = aa[;2,4,5];
aa[;100,3,1,1,3] = [];
if (all (all (aa != c))) { error ("error in empty matrix behavior"); }

v = 1:10;
c = v;
c[6,6,3] = [];
if (all (all (c != [1,2,4,5,7,8,9,10]))) {
  error ("error in empty matrix behavior");
}
printf("\tpassed empty matrix tests...\n");

//
// ------------------- Element-by-element multiply tests -------------------
//

a = [1,2,3];
b = [1,2,3;4,5,6;7,8,9];
c = [1,4,9;4,10,18;7,16,27];
a
if (!all (all (a .* b == c))) { error (); }
if (!all (all (b .* a == c))) { error (); }

za = a + urandom(a)*1j;
zb = b + urandom(b)*1j;

if (!all (all (za .* zb == [za;za;za] .* zb))) { error (); }
if (!all (all (zb .* za == zb .* [za;za;za]))) { error (); }
if (!all (all (a .* zb == [a;a;a] .* zb))) { error (); }
if (!all (all (zb .* a == zb .* [a;a;a]))) { error (); }
if (!all (all (za .* b == [za;za;za] .* b))) { error (); }
if (!all (all (b .* za == b .* [za;za;za]))) { error (); }

printf ("\tpassed matrix row-wise multiplication test...\n");

a = [1;2;3];
b = [1,2,3;4,5,6;7,8,9];

if (!all (all (a .* b == [a,a,a] .* b))) { error (); }
if (!all (all (b .* a == b .* [a,a,a]))) { error (); }

za = a + urandom (a)*1j;
zb = b + urandom (b)*1j;

if (!all (all (za .* zb == [za,za,za] .* zb))) { error (); }
if (!all (all (zb .* za == zb .* [za,za,za]))) { error (); }
if (!all (all (za .* b == [za,za,za] .* b))) { error (); }
if (!all (all (b .* za == b .* [za,za,za]))) { error (); }
if (!all (all (a .* zb == [a,a,a] .* zb))) { error (); }
if (!all (all (zb .* a == zb .* [a,a,a]))) { error (); }

printf("\tpassed matrix column-wise multiplication test...\n");

//
// --------------------- Test Divide -------------------------------
//

a = [1,2,3;4,5,6;7,8,9];
b = [4,5,6;7,8,9;10,11,12];
c = [ 48,  54,  60 ;
     111, 126, 141 ;
     174, 198, 222 ];

if (any (any (a / 1 != a))) { error (); }
if (1/2 != 0.5) { error (); }
if (any (any ((2*a)./a != 2))) { error (); }

printf ("\tpassed divide test...\n");

a = [1,2,3];
b = [1,2,3;4,6,6;7,8,9];
c = [1,1,1;4,3,2;7,4,3];

if (!all (all (b ./ a == c))) { error (); }
if (!all (all ([a;a;a] ./ b == a ./ b))) { error (); }
if (!all (all (b ./ [a;a;a] == b ./ a))) { error (); }

za = a + urandom (a)*1j;
zb = b + urandom (b)*1j;

if (!all (all ([za;za;za] ./ zb == za ./ zb))) { error (); }
if (!all (all (zb ./ [za;za;za] == zb ./ za))) { error (); }
if (!all (all ([a;a;a] ./ zb == a ./ zb))) { error (); }
if (!all (all (zb ./ [a;a;a] == zb ./ a))) { error (); }
if (!all (all ([za;za;za] ./ b == za ./ b))) { error (); }
if (!all (all (b ./ [za;za;za] == b ./ za))) { error (); }

printf("\tpassed matrix row-wise divide test...\n");

a = [1;2;3];
b = [1,2,3;4,5,6;7,8,9];

if (!all (all ([a,a,a] ./ b == a ./ b))) { error (); }
if (!all (all (b ./ [a,a,a] == b ./ a))) { error (); }

za = a + urandom (a)*1j;
zb = b + urandom (b)*1j;

if (!all (all ([za,za,za] ./ zb == za ./ zb))) { error (); }
if (!all (all (zb ./ [za,za,za] == zb ./ za))) { error (); }
if (!all (all ([za,za,za] ./ b == za ./ b))) { error (); }
if (!all (all (b ./ [za,za,za] == b ./ za))) { error (); }
if (!all (all ([a,a,a] ./ zb == a ./ zb))) { error (); }
if (!all (all (zb ./ [a,a,a] == zb ./ a))) { error (); }

printf("\tpassed matrix column-wise division test...\n");

//
// ---------------------------- Power Tests  ------------------------
//

a = (-1:1);

if (!all((1).^1 == 1)) { error ("power problems!"); }
if (!all((0).^1 == 0)) { error ("power problems!"); }
if (!all((-1).^1 == -1)) { error ("power problems!"); }

if (!all((1).^(0i) == 1)) { error ("power problems!"); }
if (!all((0).^(0i) == 1)) { error ("power problems!"); }
if (!all((-1).^(0i) == 1)) { error ("power problems!"); }

if (!all((1)^1 == 1)) { error ("power problems!"); }
if (!all((0)^1 == 0)) { error ("power problems!"); }
if (!all((-1)^1 == -1)) { error ("power problems!"); }

if (!all((1)^(0i) == 1)) { error ("power problems!"); }
if (!all((0)^(0i) == 1)) { error ("power problems!"); }
if (!all((-1)^(0i) == 1)) { error ("power problems!"); }

if (!all(a.^1 == a)) { error ("power problems!"); }
if (!all(a.^0 == [1,1,1])) { error ("power problems!"); }
if (!all(a.^(0i) == [1,1,1])) { error ("power problems!"); }

printf("\tpassed power operation test...\n");

//
// ------------------------- Sparse-Real Tests ----------------------
//

// Create a sparse (but dense storage matrix)

a1 = [ 0,  0,  0,  0;
       2,  0,  0, -3;
       0,  0, 10,  0;
       0,  0,  0, -5 ];

a2 = [ 0,  0,  3,  0;
      19,  0, 14,  0;
       2,  0,  0,  0;
       0,  0, -1,  0 ];

a3 = randsvd(4,10);

a4 = sparsify (randsvd(20,10), 0.8);
a5 = sparsify (randsvd(20,10), 0.8);

// Make sparse versions...

s1 = sparse (a1);
s2 = sparse (a2);

s3 = sparse (a3);

s4 = sparse (a4);
s5 = sparse (a5);

z1 = a4./2 - 2.5*a5*1j;
sz1 = sparse (z1);

z2 = a4 + a5*1j;
sz2 = sparse (z2);


//
// Test sparse-dense conversions...
//

if (!all (all (full (s1) == (a1))))
{ error ("sparse-dense conversion error"); }

if (!all (all (full (sz1) == (z1))))
{ error ("sparse-dense conversion error"); }

printf ("\tpassed sparse-real-dense conversion test...\n");

if (!all (all (full (sparse(a1)) == (a1))))
{ error ("dense-sparse conversion error"); }

if (!all (all (full (sparse(a4)) == (a4))))
{ error ("dense-sparse conversion error"); }

if (!all (all (full (sparse(z1)) == (z1))))
{ error ("dense-sparse conversion error"); }

printf ("\tpassed dense-real-sparse conversion test...\n");

//
// Test sparse append and stack...
//

if (!all (all (full ([s1;s2]) == ([a1;a2]))))
{ error ("sparse-append conversion error"); }

if (!all (all (full ([s4;s5]) == ([a4;a5]))))
{ error ("sparse-append conversion error"); }

if (!all (all (full ([sz1;sz1]) == ([z1;z1]))))
{ error ("sparse-append conversion error"); }

printf ("\tpassed sparse-append test...\n");

if (!all (all (full ([s1,s2]) == ([a1,a2]))))
{ error ("sparse-stack conversion error"); }

if (!all (all (full ([s4,s5]) == ([a4,a5]))))
{ error ("sparse-stack conversion error"); }

if (!all (all (full ([sz1,sz1]) == ([z1,z1]))))
{ error ("sparse-stack conversion error"); }

printf ("\tpassed sparse-stack test...\n");

//
// Test sparse assignment...
//

a4[3:5;3:5] = [1:3;1:3;1:3];
s4[3:5;3:5] = [1:3;1:3;1:3];

z1[1,3;1] = [100;200];
sz1[1,3;1] = [100;200];

if (!all (all (full (s4) == a4)))
{ error ("sparse-assignment error"); }

if (!all (all (s4 == a4)))
{ error ("sparse-assignment error"); }

if (!all (all (full (sz1) == z1)))
{ error ("sparse-assignment error"); }

if (!all (all (sz1 == z1)))
{ error ("sparse-assignment error"); }

a4[6;3] = 100;
s4[6;3] = 100;

if (!all (all (full (s4) == a4)))
{ error ("sparse-assignment error"); }

a4[7,5,6;3,6,1] = [100:102;200:202;300:302];
s4[7,5,6;3,6,1] = [100:102;200:202;300:302];

if (!all (all (full (s4) == a4)))
{ error ("sparse-assignment error"); }

a5[1,2,3;3,6,1] = [100:102;200:202;300:302];
s5[1,2,3;3,6,1] = [100:102;200:202;300:302];

if (!all (all (full (s4) == a4)))
{ error ("sparse-assignment error"); }

printf ("\tpassed sparse-assignment test...\n");

//
// Test sparse add...
//

if (!all (all (full (s1 + s2) == (a1 + a2))))
{
  error ("sparse-add error");
}
printf ("\tpassed sparse-real-add test...\n");

if (!all (all (full (sz1 + sz2) == (z1 + z2))))
{
  error ("sparse-add error");
}
printf ("\tpassed sparse-complex-add test...\n");

//
// Test sparse subtract...
//

if (!all (all (full (s1 - s2) == (a1 - a2))))
{
  error ("sparse-add error");
}
printf ("\tpassed sparse-real-subtract test...\n");

if (!all (all (full (sz1 - sz2) == (z1 - z2))))
{
  error ("sparse-add error");
}
printf ("\tpassed sparse-complex-subtract test...\n");

// Test sparse multiply...

if (!all (all (full (s1 * s2) == (a1 * a2))))
{
  error ("sparse-multiply error");
}
if (!all (all (abs( full(s3*s2*s1)-(a3*a2*a1) )< 1e-14)))
{
  error ("sparse-multiply error");
}
printf ("\tpassed sparse-real-multiply test...\n");

if (max (max (abs (sz1*sz2 - z1*z2))) > 1e-12)
{
  error ("sparse-complex-multiply error");
}

printf ("\tpassed sparse-complex-multiply test...\n");

// Test sparse transpose...

if (!all (all (full (s1'') == (a1))))
{
  error ("sparse-transpose error");
}
printf ("\tpassed sparse-real-transpose test...\n");

if (!all (all (full (sz1'') == (z1))))
{
  error ("sparse-transpose error");
}
printf ("\tpassed sparse-complex-transpose test...\n");

//
// Test spconvert ()
//

a1 = [3, 3, 13;
      1, 4, 15;
      1, 2, 10;
      3, 1, 14;
      3, 4, 12;
      2, 4, 16;
      2, 1, 11];

c1 = [ 0, 10,  0, 15;
      11,  0,  0, 16;
      14,  0, 13, 12];

z1 = [3, 3, 13, 2;
      1, 4, 15, 3;
      1, 2, 10, 10;
      3, 1, 14, 10;
      3, 4, 12, 2;
      2, 4, 16, 100;
      2, 1, 11, -11];

cz1 = [ 0,      10+10i,     0, 15+3i;
       11-11j,       0,     0, 16+100j;
       14+10i,       0, 13+2i, 12+2j];

s1 = spconvert (a1);
if (!all (all (full(s1) == c1))) { error ("spconvert error"); }

sz1 = spconvert (z1);
if (!all (all (full(sz1) == cz1))) { error ("spconvert error"); }
printf ("\tpassed spconvert test...\n");

//
// REAL...
//

a = urandom(50,50);
b = urandom(50,50);

sa = sparse (a);
sb = sparse (b);

a1 = sparsify (a, 0.3);
sa1 = sparse (a1);

a2 = sparsify (a, 0.5);
sa2 = sparse (a2);

b1 = sparsify (b, 0.4);
sb1 = sparse (b1);

b2 = sparsify (b, 0.6);
sb2 = sparse (b2);

// ==

if (!all (all (sa == sa))) { error ("sparse equality error"); }
if (!all (all (sa1 == sa1))) { error ("sparse equality error"); }
if (!all (all (sa2 == sa2))) { error ("sparse equality error"); }

if (!all (all (sb == sb))) { error ("sparse equality error"); }
if (!all (all (sb1 == sb1))) { error ("sparse equality error"); }
if (!all (all (sb2 == sb2))) { error ("sparse equality error"); }

if (!all (all (sa == a))) { error ("sparse equality error"); }
if (!all (all (sa1 == a1))) { error ("sparse equality error"); }
if (!all (all (sa2 == a2))) { error ("sparse equality error"); }

if (!all (all (sb == b))) { error ("sparse equality error"); }
if (!all (all (sb1 == b1))) { error ("sparse equality error"); }
if (!all (all (sb2 == b2))) { error ("sparse equality error"); }

// check un-usual cases...
if (!all (all (sparse([]) == []))) { error ("sparse equality error"); }
if (!all (all (sparse([2.5]) == [2.5,2.5])))
{ error ("sparse equality error"); }
if (!all (all ([2.5,2.5] == sparse([2.5]))))
{ error ("sparse equality error"); }
if (!all (all (sparse([2.5]) == [2.5,2.5;2.5,2.5])))
{ error ("sparse equality error"); }

printf ("\tpassed sparse-real equality tests...\n");

// !=

if (all (all (sa != sa))) { error ("sparse in-equality error"); }
if (all (all (sa1 != sa1))) { error ("sparse in-equality error"); }
if (all (all (sa2 != sa2))) { error ("sparse in-equality error"); }
if (all (all (sb2 != sb2))) { error ("sparse in-equality error"); }

if (all (all (sb != sb))) { error ("sparse in-equality error"); }
if (all (all (sb1 != sb1))) { error ("sparse in-equality error"); }
if (all (all (sb2 != sb2))) { error ("sparse in-equality error"); }

if (all (all (sa != a))) { error ("sparse in-equality error"); }
if (all (all (sa1 != a1))) { error ("sparse in-equality error"); }
if (all (all (sa2 != a2))) { error ("sparse in-equality error"); }

if (all (all (sb != b))) { error ("sparse in-equality error"); }
if (all (all (sb1 != b1))) { error ("sparse in-equality error"); }
if (all (all (sb2 != b2))) { error ("sparse in-equality error"); }

// check un-usual cases...
if (all (all (sparse([]) != []))) { error ("sparse equality error"); }
if (all (all (sparse([2.5]) != [2.5,2.5])))
{ error ("sparse equality error"); }
if (all (all ([2.5,2.5] != sparse([2.5]))))
{ error ("sparse equality error"); }
if (all (all (sparse([2.5]) != [2.5,2.5;2.5,2.5])))
{ error ("sparse equality error"); }

printf ("\tpassed sparse-real in-equality tests...\n");

//
// COMPLEX...
//

za = urandom(50,50) + gaussian(50,50)*1j;
zb = urandom(50,50) + gaussian(50,50)*1j;

sza = sparse (za);
szb = sparse (zb);

za1 = sparsify (za, 0.3);
sza1 = sparse (za1);

za2 = sparsify (za, 0.5);
sza2 = sparse (za2);

zb1 = sparsify (zb, 0.4);
szb1 = sparse (zb1);

zb2 = sparsify (zb, 0.6);
szb2 = sparse (zb2);

if (!all (all (sza == sza))) { error ("sparse equality error"); }
if (!all (all (sza1 == sza1))) { error ("sparse equality error"); }
if (!all (all (sza2 == sza2))) { error ("sparse equality error"); }

if (!all (all (szb == szb))) { error ("sparse equality error"); }
if (!all (all (szb1 == szb1))) { error ("sparse equality error"); }
if (!all (all (szb2 == szb2))) { error ("sparse equality error"); }

if (!all (all (sza == za))) { error ("sparse equality error"); }
if (!all (all (sza1 == za1))) { error ("sparse equality error"); }
if (!all (all (sza2 == za2))) { error ("sparse equality error"); }

if (!all (all (szb == zb))) { error ("sparse equality error"); }
if (!all (all (szb1 == zb1))) { error ("sparse equality error"); }
if (!all (all (szb2 == zb2))) { error ("sparse equality error"); }

// check un-usual cases...
if (!all (all (sparse([] + []*1j) == []))) { error ("sparse equality error"); }
if (!all (all (sparse([2.5]+2.5j) == [2.5,2.5]+[2.5j,2.5j])))
{ error ("sparse equality error"); }
if (!all (all ([2.5,2.5]+[2.5j,2.5j] == sparse([2.5]+2.5j))))
{ error ("sparse equality error"); }
if (!all (all (sparse([2.5]+2.5j) == [2.5,2.5;2.5,2.5]+[2.5j,2.5j;2.5j,2.5j])))
{ error ("sparse equality error"); }

printf ("\tpassed sparse-complex equality tests...\n");

if (all (all (sza != sza))) { error ("sparse in-equality error"); }
if (all (all (sza1 != sza1))) { error ("sparse in-equality error"); }
if (all (all (sza2 != sza2))) { error ("sparse in-equality error"); }
if (all (all (szb2 != szb2))) { error ("sparse in-equality error"); }

if (all (all (szb != szb))) { error ("sparse in-equality error"); }
if (all (all (szb1 != szb1))) { error ("sparse in-equality error"); }
if (all (all (szb2 != szb2))) { error ("sparse in-equality error"); }

if (all (all (sza != za))) { error ("sparse in-equality error"); }
if (all (all (sza1 != za1))) { error ("sparse in-equality error"); }
if (all (all (sza2 != za2))) { error ("sparse in-equality error"); }

if (all (all (szb != zb))) { error ("sparse in-equality error"); }
if (all (all (szb1 != zb1))) { error ("sparse in-equality error"); }
if (all (all (szb2 != zb2))) { error ("sparse in-equality error"); }

// check un-usual cases...
if (all (all (sparse([]) != []))) { error ("sparse equality error"); }
if (all (all (sparse([2.5*1j]) != [2.5,2.5]*1j)))
{ error ("sparse equality error"); }
if (all (all ([2.5,2.5]*1j != sparse([2.5*1j]))))
{ error ("sparse equality error"); }
if (all (all (sparse([2.5*1j]) != [2.5,2.5;2.5,2.5]*1j)))
{ error ("sparse equality error"); }

printf ("\tpassed sparse-complex in-equality tests...\n");

//
// ---------------------------- LIST Tests --------------------------
//
//  List creation

m = [1,2,3;4,5,6;7,8,9];
listest = << << 11; 12 >>; << 21; 22>> >>;
if( listest.[1].[2] != 12 ) { error(); }
if(any(<<a=10;b=1:4;c=[1,2,3;4,5,6;7,8,9]>>.b != [1,2,3,4])) { error(); }
mlist.[0] = m;
if(any(any(mlist.[0] != m))) { error(); }

// Test list functions...

listest.fun = function ( a ) { return 2*a; };
if (!(listest.fun(0.5) == 2*0.5)) { error() ; }

// Test list members, etc...

a = << b = 1:3 ; c = [1:3;4:6;7:9] ; d = 1+2j >> ;
x = a;
x.c = 13;
if (any (any (a.c != [1:3;4:6;7:9]))) { error (); }

for (i in members (a))
{
  if (entinfo (a).refc != 1)
  {
    error ("list member reference count wrong");
  }
}

printf ("\tpassed list test...\n");

//
// Test sub-lists...
//

// Real dense matrix...
m = [1,2,3;4,5,6;7,8,9];
m.rid = 1:3;
m.cid = m.rid;
m.list = << pi=3.14; 2*3.14; urandom(3,3) >>;
x = members(m);
if (!all(m.[x[9]] == (1:3))) { error (); }

// Complex dense matrix...
mc = [1,2,3;4,5,6;7,8,9] + urandom(3,3)*1j;
mc.rid = 1:3;
mc.cid = mc.rid;
mc.list = << pi=3.14; 2*3.14; urandom(3,3) >>;
mc.rid[2] = 13;
x = members(mc);

if (!all(mc.[x[7]] == (1:3))) { x? m.[x[6]] ?error (); }
if (!all(mc.[x[9]] == [1,13,3])) { x? m.[x[8]]? error (); }

// String matrix...
s = ["a", "b", "c"];
s.rid = 1:3;
s.cid = s.rid;
s.list = << pi=3.14; 2*3.14; urandom(3,3) >>;
s.rid[2] = 14;
x = members(s);

if (!all(s.[x[7]] == (1:3))) { x? error (); }
if (!all(s.[x[9]] == [1,14,3])) { x? error (); }

// Function ...
f = function (a) { return 2*a; };
f.rid = 1:3;
f.cid = f.rid;
f.list = << pi=3.14; 2*3.14; urandom(3,3) >>;
f.rid[2] = 15;
x = members(f);

if (!all(f.[x[4]] == (1:3))) { x? error (); }
if (!all(f.[x[6]] == [1,15,3])) { x? error (); }

printf ("\tpassed sub-list test...\n");

//
// -------------------------- Test printf () --------------------------
//

sprintf (tmp, "%*.*d %*.*d %s %*.*f f\n", 5,3,2, 8,7,3, "string", 3, 4, 1234e-2);
if (!(tmp == "  002  0000003 string 12.3400 f\n")) { error ("sprintf() error"); }

sprintf (tmp, "%*.*d %*.*d %s %*.*f f\n", [5],[3],[2], [8],[7],[3], ...
         ["string"], [3], [4], [1234e-2]);
if (!(tmp == "  002  0000003 string 12.3400 f\n")) { error ("sprintf() error"); }

sprintf (tmp, "%*.*d %*.*d %s %*.*f f\n", [5,1],[3,2],[2,2], [8],[7],[3], ...
         ["string"], [3], [4], [1234e-2,4]);
if (!(tmp == "  002  0000003 string 12.3400 f\n")) { error ("sprintf() error"); }

sprintf (tmp, "%*.*d %*.*d %s %*.*f f\n", [5+2i,1],[3,2+4i],[2,2], [8],[7],[3], ...
         ["string"], [3+2i], [4], [1234e-2+12j,4]);
if (!(tmp == "  002  0000003 string 12.3400 f\n")) { error ("sprintf() error"); }

printf("\tpassed sprintf test...\n");

//
// ------------------------- Test num2str()  ----------------------------
//

num2str = function ( N )
{
  local (i, j, s, stmp)

  if (class (N) == "num" && type (N) == "real")
  {
    s = [];
    for (i in 1:N.nr)
    {
      for (j in 1:N.nc)
      {
        sprintf (stmp, "%.4g", N[i;j]);
        s[i;j] = stmp;
      }
    }
    return s;
  else
    error ("num2str: argument must be numeric-real");
  }
};

if ((x = num2str(123)) != "123") { error (); }
printf("\tpassed num2str test...\n");

//
// ------------------------- Test strtod()  ----------------------------
//

if (123.456 != strtod ("123.456")) { error (); }
if (!all (all ([1,2;3,4] == strtod (["1","2";"3","4"]))))
{
  error ();
}

printf("\tpassed strtod test...\n");

//
// ------------------------- Test getline()  ---------------------------
//

close( "test.getline" );

x = getline( "test.getline" );
if ( x.[1] !=  123.456 ) { error(); }
if ( x.[2] != -123.456 ) { error(); }
if ( x.[3] !=  123.456 ) { error(); }

x = getline( "test.getline" );
if ( x.[1] !=  .123 ) { error(); }
if ( x.[2] != -.123 ) { error(); }
if ( x.[3] !=  .123 ) { error(); }

x = getline( "test.getline" );
if ( x.[1] !=  123 ) { error(); }
if ( x.[2] != -123 ) { error(); }
if ( x.[3] !=  123 ) { error(); }

x = getline( "test.getline" );
if ( x.[1] !=  1e6 ) { error(); }
if ( x.[2] != -1e6 ) { error(); }
if ( x.[3] !=  1e6 ) { error(); }

x = getline( "test.getline" );
if ( x.[1] !=  1e5 ) { error(); }
if ( x.[2] != -1e5 ) { error(); }
if ( x.[3] !=  1e5 ) { error(); }

x = getline( "test.getline" );
if ( x.[1] !=  123.456e3 ) { error(); }
if ( x.[2] != -123.456e3 ) { error(); }
if ( x.[3] !=  123.456e3 ) { error(); }

x = getline( "test.getline" );
if ( x.[1] !=  123.456e3 ) { error(); }
if ( x.[2] != -123.456e3 ) { error(); }
if ( x.[3] !=  123.456e3 ) { error(); }

x = getline( "test.getline" );
if ( x.[1] !=  123.456e-3 ) { error(); }
if ( x.[2] != -123.456e-3 ) { error(); }
if ( x.[3] !=  123.456e-3 ) { error(); }

x = getline( "test.getline" );
if ( x.[1] !=  .123e3 ) { error(); }
if ( x.[2] != -.123e3 ) { error(); }
if ( x.[3] !=  .123e3 ) { error(); }

x = getline( "test.getline" );
if ( x.[1] !=  123e3 ) { error(); }
if ( x.[2] != -123e3 ) { error(); }
if ( x.[3] !=  123e3 ) { error(); }

x = getline( "test.getline" );
if ( x.[1] != "123abc" ) { error(); }
if ( x.[2] != "abc123" ) { error(); }
if ( x.[3] != "123.abc" ) { error(); }

x = getline( "test.getline" );
if ( x.[1] != "quoted string" ) { error(); }
if ( x.[2] != "q string with escapes \n \t \" " ) { error(); }

x = getline( "test.getline" );
if ( x.[1] != "quoted string" ) { error(); }
if ( x.[2] !=  1.23e3 ) { error(); }
if ( x.[3] !=  100 ) { error(); }
if ( x.[4] != "q string with escapes \n \t \" " ) { error(); }
if ( x.[5] !=  200 ) { error(); }

printf("\tpassed getline() test...\n");

//
// ---------------------- Test strsplt()  --------------------------
//

x = getline ("test.getline", 0);

if (strlen (x) != 79) { error (); }
if ((strsplt(x, 13)).n != 6) { error (); }
if ((strsplt(x,[13,12,14,13,15,11] )).n != 6) { error (); }
if ((strsplt(x, ".")).n != 7) { error (); }

printf("\tpassed strsplt() test...\n");

//
// ----------------------- Test findstr() --------------------------
//

s = "How much wood would a woodchuck chuck?";
if (any(findstr(s, "a") != 21)) { error("error in findstr"); }
if (any(findstr(s, "wood") != [10,23])) { error ("error in findstr"); }
if (any(findstr("wood", s) != [10,23])) { error ("error in findstr"); }
if (any(findstr(s, "Wood") != [] )) { error ("error in findstr"); }
if (any(findstr(s, " ") != [4,9,14,20,22,32])) { error ("error in findstr"); }

printf("\tpassed findstr() test...\n");

//
// ---------------------- Test fread/fwrite -------------------------
//

// Create some data...
x = 1:128;

// Write out some data...
fwrite ("testme.io", "int", x);

// Read it in...
xcheck = fread ("testme.io", 33333, "int");

// Check it out...
if (!all (x == xcheck))
{
  fprintf ("Warning!... possible problem with fwrte/fread\n");
}

// Write out some data...
fwrite ("testme.io", "short int", x);

// Read it in...
xcheck = fread ("testme.io", 245, "short int");

// Check it out...
if (!all (x == xcheck)) {
  printf ("Warning!... possible problem with fwrte/fread\n");
}

// Write out some data...
fwrite ("testme.io", "unsigned char", x);

// Read it in...
xcheck = fread ("testme.io", 312, "unsigned char");

// Check it out...
if (!all (x == xcheck)) {
  printf ("Warning!... possible problem with fwrte/fread\n");
}

// Write out some data...
fwrite ("testme.io", "float", x);

// Read it in...
xcheck = fread ("testme.io", 1111, "float");

// Check it out...
if (!all (x == xcheck)) {
  fprintf ("Warning!... possible problem with fwrte/fread\n");
}

// Write out some data...
fwrite ("testme.io", "double", x);

// Read it in...
xcheck = fread ("testme.io", 234, "double");


// Check it out...
if (!all (x == xcheck)) {
  printf ("Warning!... possible problem with fwrte/fread\n");
}

printf("\tpassed fwrite/fread test...\n");

//
// ------------------------- Test size()  ---------------------------
//

aaa = [1,2;3,4;5,6];
zzz = [1+2i, 3+4i, 5+6i];
sss = ["s1", "ss2"; "sss3", "ssss4"];

if (any ((size (aaa) != [3,2]))) { error (); }
if (any ((size (zzz) != [1,3]))) { error (); }
if (any ((size (sss) != [2,2]))) { error (); }
if (any ((size (listest) != 3))) { error (); }

printf("\tpassed size() test...\n");

//
// ------------------------- Test class()  ---------------------------
//

if (class (1) != "num") { error (); }
if (class ([1]) != "num") { error (); }
if (class ([1+2j]) != "num") { error (); }
if (class ("string") != "string") { error (); }
if (class (["string"]) != "string") { error (); }
if (class (members) != "function") { error (); }
if (class (listest) != "list") { error (); }

printf("\tpassed class() test...\n");

//
// ------------------------- Test type()  ---------------------------
//

if (type ([1]) != "real") { error (); }
if (type (1+2j) != "complex") { error (); }
if (type ("str") != "string") { error (); }
if (type (<< 1 >>) != "") { error (); }
if (type (<< 1 ; type="my-list">>) != "my-list") { error (); }
if (type (members) != "builtin") { error (); }
if (type (listest.fun) != "user") { error (); }

printf("\tpassed type() test...\n");

//
// -------------------------- Test find () ------------------------------
//

if (find ([0,1]) != 2) { error ("find() output incorrect"); }
if (find ([1,0]) != 1) { error ("find() output incorrect"); }
if (find ([0,1+1i]) != 2) { error ("find() output incorrect"); }
if (find ([1+0i,0]) != 1) { error ("find() output incorrect"); }
if (find ([0+1i,0]) != 1) { error ("find() output incorrect"); }

printf("\tpassed find() test....\n");

//
// ------------------------- Test reshape()  ---------------------------
//

if (any (any (reshape(1:4, 2, 2) != [1,3;2,4]))) { error (); }
if (any (any (reshape([1+2i,3+4i,5+6i,7+8i], 2, 2) != [1+2i,5+6i;3+4i,7+8i])))
{ error (); }
if (any (any (reshape(["1","2","3","4"], 2, 2) != ["1","3";"2","4"])))
{ error (); }

printf("\tpassed reshape() test...\n");

//
// ------------------------- Test mod()  ---------------------------
//

if (mod (5,10) != 5) { error (); }
if (mod (1,1) != 0) { error ("mod() output incorrect"); }
if (mod (4,2) != 0) { error ("mod() output incorrect"); }
if (mod (3,2) != 1) { error ("mod() output incorrect"); }
if (mod (5,3) != 2) { error ("mod() output incorrect"); }

printf("\tpassed mod() test...\n");

//
// -------------------------- Test floor () ----------------------------
//

tsm = spconvert([1,1,3.6;2,2,4.1;3,1,5.2]);

if (floor (1.9999) != 1) { error ("floor() output incorrect"); }
if (!all (all (floor ([1.99,1.99;2.99,2.99]) == [1,1;2,2]))) {
  error ("floor output incorrect");
}
if (!all (all (floor (tsm) == spconvert([1,1,3;2,2,4;3,1,5])))) {
  error ("floor output incorrect");
}

printf("\tpassed floor() test...\n");

//
// ------------------------- Test ceil()  ---------------------------
//

if (ceil (1.9999) != 2) { error ("ceil() output incorrect"); }
if (!all (all (ceil ([1.99,1.99;2.99,2.99]) == [2,2;3,3]))) {
  error ("ceil output incorrect");
}
if (!all (all (ceil (tsm) == spconvert([1,1,4;2,2,5;3,1,6])))) {
  error ("floor output incorrect");
}

printf("\tpassed ceil() test...\n");

//
// -------------------------- Test round () ----------------------------
//

if (round (1.8) != 2) { error ("round() output incorrect"); }
if (round (1.4) != 1) { error ("round() output incorrect"); }
if (!all (all (round ([1.99,1.99;2.4,2.4]) == [2,2;2,2]))) {
  error ("round output incorrect");
}
if (!all (all (round (tsm) == spconvert([1,1,4;2,2,4;3,1,5])))) {
  error ("floor output incorrect");
}

printf("\tpassed round() test.....\n");

//
// ------------------------- Test int()  ---------------------------
//

if (int (1.9999) != 1) { error ("int() output incorrect"); }
if (!all (all (int ([1.99,1.99;2.99,2.99]) == [1,1;2,2]))) {
  error ("int() output incorrect");
}
if (!all (all (int (tsm) == spconvert([1,1,3;2,2,4;3,1,5])))) {
  error ("floor output incorrect");
}

printf("\tpassed int() test...\n");

//
// -------------------------- Test abs () ------------------------------
//

A = urandom(5,5);
T = ( A == abs (A));
if (!all (all (T))) { error ("abs() incorrect"); }
printf("\tpassed abs() test...\n");

//
// -------------------------- Test max () ------------------------------
//

A = [1,10,100;2,20,200;1,2,3];
B = A./2;
ZA = A + urandom (3,3)*A*1i;
ZB = B + urandom (3,3)*B*1i;
if (!all (max (A) == [2,20,200])) { error( "max() incorrect"); }
if (max (max(A)) != 200) { error ("max() incorrect"); }
if (any (any (max (A, B) != A))) { error (); }
if (any (any (max (B, A) != max (A, B)))) { error (); }
if (any (any (max (ZB, ZA) != max (ZA, ZB)))) { error (); }
if (any (any (max (B, ZA) != max (ZA, B)))) { error (); }
if (any (any (max (ZB, A) != max (A, ZB)))) { error (); }

v = urandom(1,6);
sv = sparse(v);

if (max(v) != max(sv))
{
  error("sparse max() error");
  else
  printf("\tpassed sparse-vector-real max test...\n");
}
v = urandom(1,3) + urandom(1,3)*1j;
sv = sparse(v);

if (max(v) != max(sv))
{
  error("sparse max() error");
  else
  printf("\tpassed sparse-complex-vector max test...\n");
}

v = urandom(8,6);
sv = sparse(v);

if (any(max(v) != max(sv)))
{
  error("sparse max() error");
  else
  printf("\tpassed sparse-real-matrix max test...\n");
}

v = urandom(8,2) + gaussian(8,2)*1j;
sv = sparse(v);

if (any(max(v) != max(sv)))
{
  error("sparse max() error");
  else
  printf("\tpassed sparse-complex-matrix max test...\n");
}
printf("\tpassed max() test...\n");
//srand(1);

//
// -------------------------- Test min () ------------------------------
//

if (!all (min (A) == [1,2,3])) { error( "min() incorrect"); }
if (min (min(A)) != 1) { error ("min() incorrect"); }
if (any (any (min (A, B) != B))) { error (); }
if (any (any (min (B, A) != min (A, B)))) { error (); }
if (any (any (min (ZB, ZA) != min (ZA, ZB)))) { error (); }
if (any (any (min (B, ZA) != min (ZA, B)))) { error (); }
if (any (any (min (ZB, A) != min (A, ZB)))) { error (); }

v = urandom(1,6);
sv = sparse(v);

if (min(v) != min(sv))
{
  error("sparse min() error");
  else
  printf("\tpassed sparse-vector-real min test...\n");
}
v = urandom(1,3) + urandom(1,3)*1j;
sv = sparse(v);

if (min(v) != min(sv))
{
  error("sparse min() error");
  else
  printf("\tpassed sparse-complex-vector min test...\n");
}

v = urandom(8,6);
sv = sparse(v);

if (any(min(v) != min(sv)))
{
  error("sparse min() error");
  else
  printf("\tpassed sparse-real-matrix min test...\n");
}

v = urandom(8,2) + urandom(8,2)*1j;
sv = sparse(v);

if (any(min(v) != min(sv)))
{
  error("sparse min() error");
  else
  printf("\tpassed sparse-complex-matrix min test...\n");
}

printf("\tpassed min() test....\n");

//
// -------------------------- Test maxi () -----------------------------
//

if (!all (maxi (A) == [2,2,2])) { error( "maxi() incorrect"); }
if (maxi (maxi(A)) != 1) { error ("maxi() incorrect"); }
printf("\tpassed maxi() test....\n");

//
// -------------------------- Test mini () -----------------------------
//

if (!all (mini (A) == [1,3,3])) { error( "mini() incorrect"); }
if (mini (mini(A)) != 1) { error ("mini() incorrect"); }
printf("\tpassed mini() test....\n");

//
// -------------------------- Test sum () ------------------------------
//

S = [1:4; 4:7; 8:11];
if (sum (S[1;]) != 10) { error ("sum() incorrect"); }
if (sum (S[3;]) != 38) { error ("sum() incorrect"); }
if (!all (all (sum (S) == [13,16,19,22]))) { error ("sum() incorrect"); }
printf("\tpasswd sum() test....\n");


//
// -------------------------- Test cumsum () ------------------------------
//

S = [1:4; 4:7; 8:11];
if (any(cumsum (S[1;]) != [1,3,6,10])) { error ("cumsum() incorrect"); }
if (any(cumsum (S[;3]) != [3;9;19])) { error ("cumsum() incorrect"); }
if (!all (all (cumsum (S) == [1,2,3,4;5,7,9,11;13,16,19,22])))
{
  error ("cumsum() incorrect");
}
printf("\tpassed cumsum() test....\n");

//
// -------------------------- Test prod () ------------------------------
//

S = [1:4; 4:7; 8:11];
if (prod (S[1;]) != 24) { error ("prod() incorrect"); }
if (prod (S[;3]) != 180) { error ("prod() incorrect"); }
if (prod ([1024,1]) != 1024) { error ("prod() incorrect"); }
if (prod ([1,1024]) != 1024) { error ("prod() incorrect"); }
if (prod ([1;1024]) != 1024) { error ("prod() incorrect"); }
if (prod ([1024;1]) != 1024) { error ("prod() incorrect"); }

if (!all (all (prod (S) == [32,90,180,308])))
{
  error ("prod() incorrect");
}
printf("\tpassed prod() test....\n");

//
// -------------------------- Test cumprod () ------------------------------
//

S = [1:4; 4:7; 8:11];
if (any(cumprod (S[1;]) != [1,2,6,24])) { error ("cumprod() incorrect"); }
if (any(cumprod (S[;3]) != [3;18;180])) { error ("cumprod() incorrect"); }
if (!all (all (cumprod (S) == [1,2,3,4;4,10,18,28;32,90,180,308])))
{
  error ("cumprod() incorrect");
}
printf("\tpassed cumprod() test....\n");

//
// -------------------------- Test sort () -----------------------------
//

a = [1,4,3,9];
x = sort(a);

if (all (x.val != [1,3,4,9])) { error ("sort, val problem"); }
if (all (x.ind != [1,3,2,4])) { error ("sort, ind problem"); }

z = [(7+8i), (3+4i), (1+2i), (5+6i)];
x = sort(z);

if (all (x.val != [(1+2i), (3+4i), (5+6i), (7+8i)])) { error (); }
if (all (x.ind != [3,2,4,1])) { error (); }

s = ["a", "z", "m"];
x = sort(s);

if (all (x.val != ["a", "m", "z"])) { error (); }
if (all (x.ind != [1,3,2])) { error (); }

printf("\tpassed sort() test....\n");

//
// -------------------------- Test real () -----------------------------
//

a = urandom (10,10);
b = urandom (10,10);

if (!all(all (real (a + b*1j) == a))) { error (); }
printf("\tpassed real() test....\n");

//
// -------------------------- Test imag () -----------------------------
//

if (!all(all (imag (a + b*1j) == b))) { error (); }
printf("\tpassed imag() test....\n");

//
// -------------------------- Test conj () -----------------------------
//

if (!all(all (conj (a + b*1j) == (a - b*1j)))) { error (); }
printf("\tpassed conj() test....\n");


//
// -------------------------- Test fft () -----------------------------
//

if (100 != fft(ones(100,1))[1]) { error ("error in fft()"); }
printf("\tpassed fft test...\n");

//
// -------------------------- Test norm () -----------------------------
//

if (norm (magic (4)) != sum(magic(4)[;1])) { error("error in norm"); }
if (norm (magic (4),inf()) != sum(magic(4)[;1])) { error("error in norm: inf"); }
if (norm (magic (4),"i") != sum(magic(4)[;1])) { error("error in norm: i"); }
if (norm (magic (7),"i") != sum(magic(7)[;1])) { error("error in norm: i"); }

printf("\tpassed norm test...\n");

smagic = sparse(magic (4));
if (norm (smagic) != sum(magic(4)[;1])) { error("error in norm"); }
if (norm (smagic,inf()) != sum(magic(4)[;1])) { error("error in norm: inf"); }
if (norm (smagic,"i") != sum(magic(4)[;1])) { error("error in norm: i"); }

cmagic = magic(4) + magic(4)*1i;
scmagic = sparse(cmagic);

if (norm (scmagic) != sum(abs(cmagic[;1]))) { error("error in norm"); }
if (norm (scmagic,inf()) != sum(abs(cmagic[;1]))) { error("error in norm: inf"); }
if (norm (scmagic,"i") != sum(abs(cmagic[;1]))) { error("error in norm: i"); }

printf("\tpassed sparse-norm test...\n");

//
// -------------------------- Test solve () -----------------------------
//

//
// Real - General case
//

a = urandom(10,10);
b = [1;1;1;1;1;1;1;1;1;1]; // ones(10,1);
x = a\b;
if (max (max (abs(a*x - b)))/(X*norm (a)*a.nr) > eps)
{
  printf ("\tThe condition // of a: %d\n", 1/rcond (a));
  printf ("\tA*X - B:\n");
  abs (a*x - b)
  error ("possible solve() problems\n");
}

printf("\tpassed solve() (general) test...\n");

//
// Real - Symmetric case
//

s = symm (urandom(10,10));
x = s\b;
if (max (max (abs(s*x - b)))/(X*norm (s)*s.nr) > eps)
{
  printf ("\tThe condition // of s: %d\n", 1/rcond (s));
  printf ("\tA*X - B:\n");
  abs (s*x - b)
  error ("possible solve() problems\n");
}

printf("\tpassed solve() (symmetric) test...\n");


//
// Real-Sparse - General case (all there is)
//

a = sprand(10,10,0.4)/100 + spconvert([(1:10)',(1:10)',(1:10)']);
b = [1;1;1;1;1;1;1;1;1;1];
x = a\b;
if (max (max (abs(a*x - b)))/(X*norm (a)*a.nr) > eps)
{
  printf ("\tThe condition // of a: %d\n", 1/rcond (a));
  printf ("\tA*X - B:\n");
  abs (a*x - b)
  error ("possible solve() problems\n");
}

printf("\tpassed solve() (sparse-general) test...\n");

//
// Complex - General  case
//

az = randsvd(10,10) + randsvd(10,10)*1j;
bz = urandom(10,1) + urandom(10,1)*1j;
xz = az \ bz;
if (max (max (abs (az*xz - bz)))/(X*norm (az)*az.nr) > eps)
{
  printf ("\tThe condition // of z: %d\n", 1/rcond (az));
  printf ("\tA*X - B:\n");
  abs (az*xz - bz)
  error ("possible factor/backsub problems\n");
}

printf("\tpassed solve() (Complex-General) test....\n");

//
// Complex - Symmetric  case
//

az = symm(randsvd(10,10) + randsvd(10,10)*1j);
bz = urandom(10,1) + urandom(10,1)*1j;
xz = az \ bz;
if (max (max (abs (az*xz - bz)))/(X*norm (az)*az.nr) > eps)
{
  printf ("\tThe condition // of sz: %d\n", 1/rcond (az));
  printf ("\tA*X - B:\n");
  abs (az*xz - bz)
  error ("possible factor/backsub problems\n");
}

printf("\tpassed solve() (Complex-Symmetric) test....\n");


//
// ----------------------- Test factor() / backsub() --------------------------
//

//
// Real - General case
//
a = urandom(10,10);
b = ones(10,1);
f = factor (a);
x = backsub (f,b);
if (max (max (abs(a*x - b)))/(X*norm (a)*a.nr) > eps)
{
  printf ("\tThe condition // of a: %d\n", 1/rcond (a));
  printf ("\tA*X - B:\n");
  abs (a*x - b)
  error ("possible factor/backsub problems\n");
}

//
// Real - Symmetric case
//

s = symm (urandom(10,10));
f = factor (s);
x = backsub (f,b);
if (max (max (abs(s*x - b)))/(X*norm (s)*s.nr) > eps)
{
  printf ("\tThe condition // of s: %d\n", 1/rcond (s));
  printf ("\tA*X - B:\n");
  abs (s*x - b)
  error ("possible factor/backsub problems\n");
}

s = symm (urandom(10,10));
f = factor (s, "s");
x = backsub (f,b);
if (max (max (abs(s*x - b)))/(X*norm (s)*s.nr) > eps)
{
  printf ("\tThe condition // of s: %d\n", 1/rcond (s));
  printf ("\tmax(A*X - B) = %g\n", max(abs (s*x - b)));
  printf ("possible factor/backsub problems\n");
}

//
// Complex - General  case
//

az = randsvd(10,10) + randsvd(10,10)*1j;
bz = urandom(10,1) + urandom(10,1)*1j;
fz = factor (az);
xz = backsub (fz,bz);
if (max (max (abs (az*xz - bz)))/(X*norm (az)*az.nr) > eps)
{
  printf ("\tThe condition // of z: %d\n", 1/rcond (az));
  printf ("\tA*X - B:\n");
  abs (az*xz - bz)
  error ("possible factor/backsub problems\n");
}

az = randsvd(10,10) + randsvd(10,10)*1j;
bz = urandom(10,1) + urandom(10,1)*1j;
fz = factor (az, "g");
xz = backsub (fz,bz);
if (max (max (abs (az*xz - bz)))/(X*norm (az)*az.nr) > eps)
{
  printf ("\tThe condition // of z: %d\n", 1/rcond (az));
  printf ("\tA*X - B:\n");
  abs (az*xz - bz)
  error ("possible factor/backsub problems\n");
}

printf("\tpassed factor/backsub (Complex-General) test....\n");

//
// Complex - Symmetric case
//
sz = symm (urandom(10,10) + urandom(10,10)*1j);
fz = factor(sz);
xz = backsub (fz,bz);
if (max (max (abs (sz*xz - bz)))/(X*norm (sz)*sz.nr) > eps)
{
  printf ("\tThe condition // of sz: %d\n", 1/rcond (sz));
  printf ("\tA*X - B:\n");
  abs (sz*xz - bz)
  error ("possible factor/backsub problems\n");
}

sz = symm (urandom(10,10) + urandom(10,10)*1j);
fz = factor(sz, "s");
xz = backsub (fz,bz);
if (max (max (abs (sz*xz - bz)))/(X*norm (sz)*sz.nr) > eps)
{
  printf ("\tThe condition // of sz: %d\n", 1/rcond (sz));
  printf ("\tA*X - B:\n");
  abs (sz*xz - bz)
  error ("possible factor/backsub problems\n");
}

printf("\tpassed factor/backsub (Complex-Symmetric) test....\n");

//
// -------------------------- Test eig () -----------------------------
//

//
// REAL Standard symmetric eigenvalue problem
//

a = urandom(10,10);
ta = trace (a);
sa = symm (a);
tsa = trace (sa);

tol = 1.e-6;

if (!(ta < sum(eig(a).val) + tol && ta > sum(eig(a).val) - tol))
{
  error ("error in eig");
}
printf("\tpassed eig() (Real standard non-symmetric eigenvalue problem) test...\n");

if (!(tsa < sum(eig(sa).val) + tol && tsa > sum(eig(sa).val) - tol))
{
  error ("error in eig");
}
printf("\tpassed eig() (Real standard symmetric eigenvalue problem) test...\n");

//
// COMPLEX Standard symmetric eigenvalue problem
//

z = urandom(10,10) + urandom(10,10)*1j;
tz = trace (z);
sz = symm (z);
tsz = trace (sz);

tol = 1.e-6;

if (!(tz < sum(eig(z).val) + tol && tz > sum(eig(z).val) - tol))
{
  error ("error in eig");
}
printf("\tpassed eig() (Complex standard non-symmetric eigenvalue problem) test...\n");

if (!(tsz < sum(eig(sz).val) + tol && tsz > sum(eig(sz).val) - tol))
{
  error ("error in eig");
}
printf("\tpassed eig() (Complex standard symmetric eigenvalue problem) test...\n");

//
// REAL Generalized eigenvalue problem
//

b = urandom(10,10) + eye(10,10)*3;
sb = symm (b);
tb = trace (b);
tsb = trace (sb);

</ lambda ; x /> = eig (a, b);

for (i in 1:10)
{
  q[;i] = (a*x[;i]) - (lambda[i]*b*x[;i]);
}

if (abs(max (max (q))) > X*100*eps)
{
  abs(max (max (q)))
  error ("error in eig (a,b)");
}

printf("\tpassed eig() (Real generalized non-symmetric eigenvalue problem) test...\n");


//
// REAL Generalized symmetric eigenvalue problem
//

</ lambda ; x /> = eig (sa, sb);

for (i in 1:10)
{
  q[;i] = (sa*x[;i]) - (lambda[i]*sb*x[;i]);
}

if (abs(max (max (q))) > X*100*eps)
{
  abs(max (max (q)))
  error ("error in eig (a,b)");
}

printf("\tpassed eig() (Real generalized symmetric eigenvalue problem) test...\n");

//
// COMPLEX Generalized non-symmetric eigenvalue problem
//

az = randsvd(10,10) + randsvd(10,10)*3*1j;
saz = symm (az);

tsaz = trace (saz);
taz = trace (az);

bz = randsvd(10,10) + randsvd(10,10)*1j + 3*eye(10,10);
sbz = symm (bz);

tsbz = trace (sbz);
tbz = trace (bz);

</ lambda ; x /> = eig (az, bz);

for (i in 1:10)
{
  q[;i] = (az*x[;i]) - (lambda[i]*bz*x[;i]);
}

if (abs(max (max (q))) > X*100*eps)
{
  abs(max (max (q)))
  error ("error in eig (a,b)");
}

printf("\tpassed eig() (Complex generalized non-symmetric eigenvalue problem) test...\n");

</ lambda ; x /> = eig (saz, sbz);

for (i in 1:10)
{
  q[;i] = (saz*x[;i]) - (lambda[i]*sbz*x[;i]);
}

if (abs(max (max (q))) > X*100*eps)
{
  abs(max (max (q)))
  error ("error in eig (a,b)");
}

printf("\tpassed eig() (Complex generalized symmetric eigenvalue problem) test...\n");

//
// ------------------------------ Test lu() ---------------------------------
//

static (swap);

lu = function ( A )
{
  local (i, l, u, pvt, x)

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

a = urandom(10,10);
lua = lu (a);
if (max (max (abs(a - lua.pvt*lua.l*lua.u)))/(X*norm (a)*a.nr) > eps)
{
  printf ("\tThe condition // of a: %d\n", 1/rcond (a));
  printf ("\tA - p*l*u:\n");
  abs (a - lua.pvt*lua.l*lua.u)
  error ("possible lu()/factor() problems\n");
}

//
// Real

az = urandom(10,10) + urandom(10,10)*1j;
luz = lu (az);
if (max (max (abs (az - luz.pvt*luz.l*luz.u)))/(X*norm (az)*az.nr) > eps)
{
  printf ("\tThe condition // of z: %d\n", 1/rcond (az));
  printf ("\tA - p*l*u:\n");
  abs (az - luz.ovt*luz.l*luz.u)
  error ("possible lu()/factor()() problems\n");
}

printf("\tpassed lu/factor test....\n");

//
// -------------------------- Test svd ()   -----------------------------
//
a = randsvd(10,10);
s = svd (a);
if (max (max (abs (s.u*diag(s.sigma)*s.vt - a)))/(X*norm (a)*a.nr) > eps)
{
  error ("possible svd() problems");
}

z = randsvd(10,10) + randsvd(10,10)*1j;
sz = svd (z);
if (max (max (abs (sz.u*diag(sz.sigma)*sz.vt - z)))/(X*norm (z)*z.nr) > eps)
{
  error ("possible svd() problems");
}

printf("\tpassed svd test....\n");

//
// -------------------------- Test hess ()  -----------------------------
//

a = randsvd(10,10);
h = hess (a);
if (max (max (abs (h.p*h.h*h.p' - a)))/(X*norm (a)*a.nr) > eps)
{
  error ("possible hess() problems");
}

z = randsvd(10,10) + randsvd(10,10)*1j;
hz = hess (z);
if (max (max (abs (hz.p*hz.h*hz.p' - z)))/(X*norm (z)*z.nr) > eps)
{
  error ("possible hess() problems");
}

printf("\tpassed hess test....\n");

//
// -------------------------- Test qr () ------------------------------
//

a = ohess(4);
qa = qr (a);
if (max (max (abs (qa.q*qa.r - a)))/(X*norm (a)*a.nr) > eps)
  { error ("possible qr() problems"); }

qa = qr (a, "p");
if (max (max (abs (qa.q*qa.r - a*qa.p)))/(X*norm (a)*a.nr) > eps)
  { error ("possible qr() problems"); }

z = ohess (4) + ohess(4)*1i;
qz = qr (z);
if (max (max (abs (qz.q*qz.r -  z)))/(X*norm (z)*z.nr) > eps)
  { error ("possible qr() problems"); }

qz = qr (z, "p");
if (max (max (abs (qz.q*qz.r -  z*qz.p)))/(X*norm (z)*z.nr) > eps)
  { error ("possible qr() problems"); }

printf("\tpassed qr test....\n");

//
// -------------------------- Test chol () -----------------------------
//

c = lehmer(10);
u = chol (c);
if (max (max (abs (u'*u - c)))/(X*norm (c)*c.nr) > eps)
{
  error ("possible chol() problems");
}

cz = lehmer(10) + lehmer(10)*1j;
cz = symm(cz);
uz = chol (cz);
if (max (max (abs (uz'*uz - cz)))/(X*norm (cz)*cz.nr) > eps)
{
  error ("possible chol() problems");
}

printf("\tpassed chol test...\n");

//
// -------------------------- Test schur () ----------------------------
//

a = randsvd (10, 10);
sa = schur (a);
if (max (max (abs (sa.z*sa.t*sa.z' - a)))/(X*norm (a)*a.nr) > eps)
  { error ("possible schur() problems"); }

z = urandom (4,4) + urandom(4,4)*1i;
sz = schur (z);
if (max (max (abs (sz.z*sz.t*sz.z' - z)))/(X*norm (z)*z.nr) > eps)
  { error ("possible schur() problems"); }

a = randsvd (10, -10);
sa = schur (a);
if (max (max (abs (sa.z*sa.t*sa.z' - a)))/(X*norm (a)*a.nr) > eps)
  { error ("possible schur() problems"); }

z = urandom (4,4) + urandom(4,4)*1i;
sz = schur (z);
if (max (max (abs (sz.z*sz.t*sz.z' - z)))/(X*norm (z)*z.nr) > eps)
  { error ("possible schur() problems"); }

printf("\tpassed schur test...\n");

//
// -------------------------- Test lyap () ------------------------------
//

lyap = function ( A, B, C )
{
  if (!exist (B))
  {
    B = A';	// Solve the special form: A*X + X*A' = -C
  }

  if ((A.nr != A.nc) || (B.nr != B.nc) || (C.nr != A.nr) || (C.nc != B.nr)) {
    error ("Dimensions do not agree.");
  }

  //
  // Schur decomposition on A and B
  //

  sa = schur (A);
  sb = schur (B);

  //
  // transform C
  //

  tc = sa.z' * C * sb.z;

  X = sylv (sa.t, sb.t, tc);

  //
  // Undo the transformation
  //

  X = sa.z * X * sb.z';

  return X;
};

//srand();
a = randsvd (20,20);
b = urandom (20,20);
c = urandom (20,20);

x = lyap (a, b, c);
if (max (max (abs (a*x + x*b + c)))/(X*norm(c)*norm(a)*norm(b)) > eps)
{
  error ("possible problems with lyap() or sylv()");
}

printf("\tpassed sylv/lyap test...\n");

//
// ------------------ Test More Advanced Functions --------------------
//

printf("\tStarting the lqr/ode test..." );
printf("\tthis will take awhile\n" );

lqr = function( a, b, q, r )
{

  m = size(a)[1]; n = size(a)[2];
  mb = size(b)[1]; nb = size(b)[2];
  mq = size(q)[1]; nq = size(q)[2];

  if ( m != mq || n != nq )
  {
    fprintf( "stderr", "A and Q must be the same size.\n" );
    quit
  }

  mr = size(r)[1]; nr = size(r)[2];
  if ( mr != nr || nb != mr )
  {
    fprintf( "stderr", "B and R must be consistent.\n" );
    quit
  }

  nn = zeros( m, nb );

  // Start eigenvector decomposition by finding eigenvectors of Hamiltonian:

  e = eig( [ a, solve(r',b')'*b'; q, -a' ] );
  v = e.vec; d = e.val;

  index = sort( real( d ) ).ind;
  d = real( d[ index ] );

  if ( !( d[n] < 0 && d[n+1] > 0 ) )
  {
    fprintf( "stderr", "Can't order eigenvalues.\n" );
    quit
  }

  chi = v[ 1:n; index[1:n] ];
  lambda = v[ (n+1):(2*n); index[1:n] ];
  s = -real(solve(chi',lambda')');
  k = solve( r, nn'+b'*s );

  return << k=k; s=s >>;

};

// Now run a little test problem.

k = 1; m = 1; c = .1;
a = [0     ,1    ,0    , 0;
    -k/m, -c/m,  k/m,  c/m;
     0,     0,    0,    1;
     k/m,  c/m, -k/m, -c/m ];
b = [ 0; 1/m; 0; 0 ];
qxx = diag( [0, 0, 100, 0] );
ruu = [1];
K = lqr( a, b, qxx, ruu ).k;

dot = function( t, x )
{
  global (a, b, K)
  return (a-b*K)*x + b*K*([1,0,1,0]');
};
odeopts=<<>>;
odeopts.erel = 1e-5;
odeopts.eabs = 1e-5;
x = odeiv (dot, , [0:15:1/64], [0,0,0,0], odeopts);

m = maxi( x[;2] );
if ( (abs( x[m;2] - 1.195 ) > 0.001)  || ...
     any (abs( x[x.nr;2,4] - 1 ) > 0.001) )
{
  printf( "\tfailed***\n" );
  else
  printf( "\tpassed the lqr/ode test....\n" );
}

//
// --------------------- Nasty Function Test ------------------------
//

printf("\tStarting Nasty Function Test...");
printf("\tthis will take awhile\n");
check = function( a, b, c, d, e, f, g, h )
{
  if ( a+b+c+d == e+f+g+h && ...
       a^2+b^2+c^2+d^2 == e^2+f^2+g^2+h^2 && ...
       a^3+b^3+c^3+d^3 == e^3+f^3+g^3+h^3 )
  {
    return 1;
  else
    return 0;
  }
};

cnt = 0;

for(a in 8:10) {
  for(b in 7:(a-1)) {
    for(c in 6:(b-1)) {
      for(d in 5:(c-1)) {
        for(e in 4:(d-1)) {
          for(f in 3:(e-1)) {
            for(g in 2:(f-1)) {
              for(h in 1:(g-1)) {
	          if(check( a, b, c, d,  e, f, g, h ) || ...
                     check( a, e, c, d,  b, f, g, h ) || ...
                     check( a, f, c, d,  e, b, g, h ) || ...
                     check( a, g, c, d,  e, f, b, h ) || ...
                     check( a, h, c, d,  e, f, g, b ) || ...
                     check( a, b, e, d,  c, f, g, h ) || ...
                     check( a, b, f, d,  e, c, g, h ) || ...
                     check( a, b, g, d,  e, f, c, h ) || ...
                     check( a, b, h, d,  e, f, g, c ) || ...
                     check( a, b, c, e,  d, f, g, h ) || ...
                     check( a, b, c, f,  e, d, g, h ) || ...
                     check( a, b, c, g,  e, f, d, h ) || ...
                     check( a, b, c, h,  e, f, g, d ) || ...
                     check( a, e, f, d,  b, c, g, h ) || ...
                     check( a, e, g, d,  b, f, c, h ) || ...
                     check( a, e, h, d,  b, f, g, c ) || ...
                     check( a, f, g, d,  e, b, c, h ) || ...
                     check( a, f, h, d,  e, b, g, c ) || ...
                     check( a, g, h, d,  e, f, b, c ) || ...
                     check( a, b, e, f,  c, d, g, h ) || ...
                     check( a, b, e, g,  c, f, d, h ) || ...
                     check( a, b, e, h,  c, f, g, d ) || ...
                     check( a, b, f, g,  e, c, d, h ) || ...
                     check( a, b, f, h,  e, c, g, d ) || ...
                     check( a, b, g, h,  e, f, c, d ) || ...
                     check( a, e, f, g,  e, f, g, h ) || ...
                     check( a, e, f, h,  e, f, g, h ) || ...
                     check( a, e, g, h,  e, f, g, h ) || ...
                     check( a, f, g, h,  e, f, g, h ) ) { cnt++; }
              }
            }
          }
        }
      }
    }
  }
}

if (!cnt)
{
  // figure out the value of cnt, and check!
  printf("\tpassed nasty function test...\n");
else
  error();
}

// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //

printf("\tDONE, Elapsed time = %10.3f seconds\n", toc() );
clearall();
