//-------------------------------------------------------------------//
// czt.r

// Syntax:	czt ( X , M , W , A )

// Description:

//	czt returns the M element Chirp z-transform, where:
//	M: The number of points in the returned spectrum.
//	W: The "spacing" between points along the z-pane spiral.
//	   A complex number along the unit circle.
//	A: the complex starting point of the evaluated transform.
//	   A complex number along the unit circle.

//	At present X must be a row or column matrix. The column-wise
//	matrix czt operation will be added later.

//      Ian Searle, 10/10/93

//	References:
//	Rabiner, L.R. and B. Gold, Theory and Application of
//	 Digital Signal Processing, Prentice-Hall, Englewood
//	 Cliffs, New Jersey, pp. 393-399, 1975.
//-------------------------------------------------------------------//

czt = function ( x, M, W, A )
{
  global (pi)

  if (min (size (x)) != 1) { error ("czt: only vector X allowed"); }
  nr = x.nr; nc = x.nc;

  if (nc == 1) 
  { 
    x = x'; 
    nr = x.nr; nc = x.nc;
    column = 1;
  else
    column = 0;
  }
  N = nc;

  if (!exist (M)) { M = length (x); }
  if (!exist (W)) { W = exp (-1j .* 2 .* pi ./ M); }
  if (!exist (A)) { A = 1; }

  //
  // Compute a power of 2 value for fft()
  //

  L = 1;
  while (L < (N + M - 1)) { L = 2 .* L; }

  //
  // Form the L-point sequence y(n)
  //

  n = 0:(N-1);
  y = [ A.^(-n) .* W.^((n.^2)./2) .* x , zeros (size (N+1:L)) ];

  //
  // Form the L-point sequence v(n)
  //

  n = 0:(M-1);
  v = W.^(-(n.^2)./2);
  v = [v, zeros (size (M+1:(L - N + 1)))];
  n = (L-N+1):(L-1);
  v[(L-N+2):L] = W.^(-((L-n).^2)./2);

  //
  // Convolution
  //

  Y = fft (y, L);
  V = fft (v, L);
  G = V.*Y;
  g = ifft (G, L);

  k = 0:(M-1);
  X = g[1:M] .* W.^((k.^2)./2);

  if (column)
  {
    return X';
  else
    return X;
  }
};
