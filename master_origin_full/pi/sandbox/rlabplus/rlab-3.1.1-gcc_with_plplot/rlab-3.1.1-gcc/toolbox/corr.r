//-------------------------------------------------------------------//
// Synopsis:    Compute correlations faster than DSP toolbox, 
//              and with bandpassing

// Usage:

//   y = corr(x)
//   Compute the autocorrelation of vector x.  Does NOT handle
//   matrices  the way xcorr (in the Signal Processing toolbox) does;
//   instead, a matrix is treated as a vector (with reshape). 

//   y = corr(x1,x2)
//   Return the cross-correlation of vectors x1 and x2, which need not
//   be  the same length.  This is similar to xcorr, but faster (the
//   latter efficiently computes 4 correlations and throws 3 away).
//   Also, if  length(x1) is different from length(x2), xcorr returns
//   a wrong-length  result.

//   y = corr(x1,x2,f) 
//   As above, but also filter by the bandpass filter f.  f is a
//   two-element vector, with values in [0,1], that specifies the
//   passband; in f, the value 1 represents the Nyquist frequency.

// 2/94 David K. Mellinger  dkm1@cornell.edu
// This file is a translation of corr.m from the Osprey toolbox.

require isreal
require nextpow2 

//-------------------------------------------------------------------//

corr = function (in1, in2, f)
{
  local (in1, in2, f)

  if (!exist (in2)) { in2 = in1; }

  x1 = in1[:];
  x2 = in2[:];
  n1 = length (x1);
  n2 = length (x2);

  // pad with 0's for fft & circular corr
  nfft = 2 * 2^nextpow2 (max (n1,n2));

  y = conj (fft ([x1;zeros(nfft-n1,1)])) .* fft([x2;zeros(nfft-n2,1)]);
  if (exist (f))
  {
    if (length(f) == 1)
    {
      f = [0, f];    // one number means lowpass filter
    }
    n = nfft/2;
    i0 = round(f[1] * n);
    i1 = round(f[2] * n);
    y[   1 : i0] = zeros (i0,1);
    y[i1+1 : n]  = zeros (n-i1,1);
    y[nfft+1-i0 : nfft]    = zeros (i0,1);
    y[nfft+1-n  : nfft-i1] = zeros (n-i1, 1);
  }

  y = ifft (y);
  y = [y[nfft+2-n1:nfft]; y[1:n2]];
  
  // Make output have same form as input.
  if (all (isreal (in1)) && all (isreal (in2)))
  {
    y = real(y);
  }

  if (in1.nr == 1)
  {
    y = y.'; 
  }

  return y;
};
