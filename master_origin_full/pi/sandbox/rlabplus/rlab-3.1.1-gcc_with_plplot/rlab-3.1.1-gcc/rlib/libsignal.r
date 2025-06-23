//
// basic signal processing toolbox
//

rmrows  = function( data, del_idx )
{
  // check for input
  if ( !exist(data) || !exist(del_idx) )
  {
    printf ("rmrows: Remove rows from a data matrix\n");
    printf ("rmrows: Format:\n");
    printf ("rmrows:	  rmrows( data, rmidx )\n");
    printf ("rmrows: where 'rmidx' contains the indices of rows to be removed.\n");
    error ("two arguments required!");
  }
  data_idx = 1:data.nr;
  data_idx = complement(data_idx, del_idx);
  if (isempty(data_idx))
  {
    return [];
  else
    return data[data_idx;];
  }
};

rmcols  = function( data, del_idx )
{
  // check for input
  if ( !exist(data) || !exist(del_idx) )
  {
    printf ("rmrows: Remove rows from a data matrix\n");
    printf ("rmrows: Format:\n");
    printf ("rmrows:    rmrows( data, rmidx )\n");
    printf ("rmrows: where 'rmidx' contains the indices of rows to be removed.\n");
    error ("two arguments required!");
  }
  data_idx = 1:data.nc;
  data_idx = complement(data_idx, del_idx);
  if (isempty(data_idx))
  {
    return [];
  else
    return data[;data_idx];
  }
};

phase = function ( A )
{
  //-------------------------------------------------------------------//
  //  Syntax:     phase ( A )
  //  Description:
  //  Calculate the phase, or the angle, of a complex vector. Normally
  //  we are tempted to do this with a simple atan2() function
  //  call. However, this method often results in discontinuities at the
  //  pi boundaries. The method used herein attempts to keep the phase
  //  continuous over the pi boundaries.
  //-------------------------------------------------------------------//

  global(const);

  // Force a column vector  A = A[:];
  phi = atan2 (imag(A), real(A));
  N = phi.n;

  // Compute the difference between each element of phi...
  dphi = phi[1:N-1] - phi[2:N];

  // Now, find the indices where dphi is greater than pi.
  I = find (abs (dphi) > 3.5);

  // Now step through the indices, compensating for the pi boundary...
  for (i in I)
  {
    if (i != 0)
    { phi = phi + 2*(const.pi) * sign(dphi[i]) .* [zeros(1,i), ones(1,N-i)]; }
  }

  return phi;
};

nextpow2 = function ( n )
{
  // Description:
  //      nextpow2 returns the first integer such tht 2^P >= N. If N is
  //      a vector, then nextpow2 computes P such that 2^P >= length(N).
  // See Also: frexp
  //-------------------------------------------------------------------//

  if (length (n) > 1)
  { n = length (n); }

  if (n<0)
  { n = abs(n); }

  return round(log(ceil(n, <<base=2>>))/log(2));
};

hilbert = function ( x )
{
  //-------------------------------------------------------------------
  // Synopsis:    Hilbert transform
  // Syntax:      hilbert ( X )
  // Description:
  //      Hilbert computes the Hilbert transform of a real
  //      sequence. Hilbert returns a complex valued analytic function
  //      Z, such that:
  //              Z = X + jXh
  //      Where X is the original series, Xh
  //      Properties: y = hilbert (x);
  //                  A = abs (y) is an envelope of x.
  //                  ph= angle (y) is the instaneous phase of x.
  //                  -N*eps < real (y)' * imag (y) < N*eps
  // Reference: J. S. Bendat, A. G. Piersol, Random Data, Analysis
  //            and Measurement Procedures. 2nd Edition, 1986,
  //            J. Wiley and Sons.
  // See Also: fft, ifft
  //-------------------------------------------------------------------

  if (x.nr==1 || x.nc==1)
  {
    // x is a single vector
    y = fft (real (x));
    n = length(y);
    if (n != 1)
    {
      b = [1; ...
          2*ones ((n-1)/2,1);  ...
              ones (1-mod(n,2),1); ...
                  zeros ((n-1)/2,1) ];
                  y = b.*y;
    }
    return ifft(y);
    else
    // x is a matrix; consider it stacked
    // column vectors
      rval = zeros(x) + 0i;
    for (j in 1:x.nc)
    {
      y = fft (real (x[;j]));
      n = length(y);
      if (n != 1)
      {
        b = [1; 2*ones ((n-1)/2,1); ones (1-mod(n,2),1); zeros ((n-1)/2,1) ];
        y = b.*y;
      }
      rval[;j] = ifft(y);
    }
    return rval;
  }
};

deconv = function(b,a)
{
  //---------------------------------------------------------------------------
  // deconv
  // Syntax: A = deconv(b,a)
  //         </q;r/> = deconv(b,a)
  // This routine performs "de=convolution." If it is called as A=deconv(b,a)
  // then it deconvolves vector a out of vector b such that the following is
  // true:
  // b=conv(q,a) + r
  // where q is the result of the deconvolution and b is the remainder.
  // If b and a are vectors of polynomial coefficients, then using deconv
  // is equivalent to polynomial division (q is the quotient and r is the
  // remainder).
  //      The results are returned in a list:
  //      A.q = results of a being deconvolved out of b.
  //      A.r = remainder
  //
  // Copyright (C), by Jeffrey B. Layton, 1994
  // Version JBL 940918
  //--------------------------------------------------------------------------

  if (nargs != 2)
  { error("deconv: requires two arguments"); }
  if (class(b)!="num")
  { error("deconv: first argument 'b' must be numeric vector"); }
  if (class(a)!="num")
  { error("deconv: second argument 'a' must be numeric vector"); }

  nb=length(b);
  na=length(a);

  if (na > nb)
  {
    q=0;
    r=b;
    return << q=q; r=r >>;
  }

  if (a[1] == 0)
  { error("deconv: in second entry first coefficient must be non-zero"); }

  // Deconvolution and polynomial division are the same operations
  // as a digital filter's impulse response B(z)/A(z):
  q=filter([1,zeros(1,nb-na)], b, a).y;

  // if (mb != 1) {
  //     q=q(:);
  // }
  r = conv(q,a);
  resize(r, b.nr, b.nc);
  r=b-r;

  return << q=q; r=r>>;
};


conv = function( A , B, opts )
{
  if (!exist(opts))
  {
    //---------------------------------------------------------------------
    //  Synopsis:   Convolve two vectors.
    //  Syntax:     conv ( A , B )
    //  Description:
    //  This routine convolves the vectors a and b. This function uses a
    //  Finite Impulse Response filter to do the convolution. This method
    //  can be faster than using the DFT when the sizes of A and B are
    //  different.
    // --------------------------------------------------------------------

    if (class(A)!="num")
    { error("conv: first argument 'a' must be numeric vector"); }
    if (class(B)!="num")
    { error("conv: second argument 'b' must be numeric vector"); }

    na = max (size (A));
    nb = max (size (B));

    if (na > nb)
    {
      if (nb > 1)
      { A[na+nb-1] = 0; }
      c = filter (A, B, 1);
    else
      if (na > 1)
      { B[na+nb-1] = 0; }
      c = filter (B, A, 1);
    }
    return c.y;

  else
    //---------------------------------------------------------------------
    //
    // conv
    //
    // Syntax: c=conv(a,b)
    //
    // This routine convolves the vectors a and b into c. It is used for
    // polynomial multiplication (algebraically they are the same
    // operation).
    // Note: Currently no error checking.
    // Written by: Ian Searle
    //---------------------------------------------------------------------

    n = length(A) + length(B) - 1;
    X = fft( A, n );
    Y = fft( B, n );

    tmp = ifft( X .* Y );
    if( type(A) == "real" && type(B) == "real")
    { tmp = real(tmp); }

    return tmp;
  }
};

poly = function( r )
{
  //
  // construct a polynomial given its zeros: uses convolution function
  //
  if(type(r)!="real" && type(r)!="complex")
  {
    printf("poly: function requires a vector as an argument!\n");
    return [];
  }
  if(r.nc!=1 && r.nr!=1)
  {
    printf("poly: first argument 'r' has to be a vector!\n");
    return [];
  }
  _a = [1,-r[1]];
  for(j in 2:length(r))
  {
    _a = conv( _a, [1, -r[j]] );
  }

  return _a;
};


faxis = function ( X, _T, axis_type )
{
  //---------------------------------------------------------------------------
  //  faxis.r
  //  Syntax: faxis ( X )
  //    faxis ( X , T )
  //    faxis ( X , T, axis_type )
  //  Description:
  //  Faxis generates a frequency axis for FFT plots.
  //
  //  X = FFT data
  //  T = sampling period (optional argument)
  //  axis_type = type of axis to create:
  //              1 = Digital Rad/s [0,2pi], 2 = Analog Radians/s
  //              3 = Analog Hertz           4 = Normalized frequency [0,2]
  //
  //  Defaults if not specified:   T = 1, axis_type = 1
  //
  //---------------------------------------------------------------------------
  global (pi)

  if (!exist (axis_type))
  { axis_type = 1; }

  if (!exist (_T))
  { T = 1; axis_type = 1; else T = _T; }

  N = length(X);
  a = (0:N-1)/N;
  a = reshape (a, X.nr, X.nc);

  if (axis_type == 3)
  {
    a = a/T;
  else if (axis_type == 2)
  {
    a = a * (2*pi/T);
  else if (axis_type == 4)
  {
    a = 2 * a;
  else
    a = a * (2*pi);
  }}}

  return a
};

rot90 = function(A,k)
{
  if (! A.nr * A.nc)
  { return A; }

  if (!exist(k))
  { k = 1;}

  k = mod(int(k),4);
  if (k)
  {
    if (k==1)
    {
      return flipud(A');
    else
      rval = flipud(A');
      return $self(rval, k-1);
    }
  else
    return A;
  }
};

conv2 = function(p,q,w)
{
  m = p.nc;
  n = p.nr;
  a = q.nc;
  b = q.nr;

  if (!exist(w))
  { w = "full"; }
  if (class(w)!="string")
  { w = "full"; }

  if (w=="same")
  {
    // q has to have odd dimensions, not necessarily equal
    if (mod(a,2)!=1 || mod(b,2)!=1)
    { w = "full"; }
  }

  if (! m*n)
  { return q; }

  if (! a*b)
  { return p; }

  if (w == "full")
  {
    l_range = 1:(n+b-1);
    k_range = 1:(m+a-1);
  else if (w == "same")
  {
    l_range = (1 + int(b/2)):(n+int(b/2));
    k_range = (1 + int(a/2)):(m+int(a/2));
  else
    return [];
  }}
  r = zeros(length(l_range),length(k_range));

  for (l in l_range)
  {
    j_range = max(1,l+1-b) : min(n,l);
    true_l = l - l_range[1] + 1;
    for (k in k_range)
    {
      i_range = max(1,k+1-a) : min(m,k);
      true_k = k - k_range[1] + 1;
      for (j in j_range)
      {
        for (i in i_range)
        {
          r[true_l;true_k] = r[true_l;true_k] + p[j;i] * q[l+1-j;k+1-i];
        }
      }
    }
  }

  return r;
};

rldeconv2 = function(observed, psf, niter, miv)
{
  // direct copy of algorithm from
  //    http://en.wikipedia.org/wiki/Richardson%E2%80%93Lucy_deconvolution
  if (!exist(miv))
  {
    // miv is not given - construct initial guess
    // initial estimate is arbitrary - uniform 50% grey works fine
    latent_est = 0.5 .* ones(observed);
  else
    if (all(size(miv) == [1,1]))
    {
      // miv is scalar - use it to construct initial guess
      latent_est = miv .* ones(observed);
    else if (all(size(miv) == size(observed)))
    {
      latent_est = miv;
    else
      // miv is not given - construct initial guess
      // initial estimate is arbitrary - uniform 50% grey works fine
      latent_est = 0.5 .* ones(observed);
    }}
  }

  // create an inverted psf
  psf_hat = rot90(psf,2);

  // iterate towards ML estimate for the latent image
  for (i in 1:niter)
  {
    est_conv      = conv2(latent_est, psf, "same");
    relative_blur = observed ./ est_conv;
    error_est     = conv2(relative_blur, psf_hat, "same");
    latent_est    = latent_est .* error_est;
  }

  return latent_est;
};


