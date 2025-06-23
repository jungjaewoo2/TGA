//-------------------------------------------------------------------//

// Synopsis:	Compute Cross-Spectra of the input sequence

// Syntax:	crossspec ( X, Y, dt, M, WIN, raw );

// Description:

//      Compute the cross-spectral density of the input vector
//      X. crossspec returns a two-column matrix. The 1st column is the
//      frequency scale and the second column is the values of the
//      one-sided spectra.

//	X:	Input sequence
//      Y:      Input sequence
//	dt:	Sample interval
//              (default = 1)
//	M:	Length of sequence to use for computations. The
//		sequence X is "cut-up" into M-point sequences
//              (default = length (X))
//      WIN:    Window type: "hann", "rect", "hamm"
//              (default = "rect")
//      raw:    If raw == 1, then output the "raw" auto-spectra.
//              Otherwise output the auto-spectra from 0-Nyq/2
//              with a time-scale.
//              (default = 0)

// Dependencies

   require window
   require detrend
   require isreal

//-------------------------------------------------------------------//

crosspec = function ( X , Y, dt, M, WIN, raw )
{
  X = X[:];		// Make sure it is a row
  Y = Y[:];

  n = length (X);
  if (!isreal (X) && !isreal (Y)) {
    error ("crosspec: X and Y must be real");
  }

  //
  // Set the defaults.
  //

  if (!exist (dt)) { dt = 1; }
  if (!exist (M)) { M = n; }
  if (!exist (WIN)) { WIN = "rect"; }
  if (!exist (raw)) { raw = 0; }

  // Some error checking

  if (M > n) { error ("crosspec: M > n not allowed"); }
  if (M < 0) { error ("crosspec: M < 0 not allowed"); }

  // Create the frequncy scale (in Hertz) up to the Nyquist

  fscale = (((0:M-1)/M)/dt)';

  k = fix (n/M);                // The number of windows
  index = 1:M;

  // Create the window

  w = window (M, WIN);          // Window specification
  U = k*sum(w.^2)/w.n;		// Normalizing scale factor

  Pxy = zeros (M, 1);

  for (i in 1:k)
  {
    xw = w .* detrend (X[index]);
    yw = w .* detrend (Y[index]);
    index = index + M;
    Xxy = (2/(M*dt)) * (conj(dt*fft (xw)) .* (dt*fft (yw)));
    Pxy = Pxy + Xxy;
  }

  if (raw)
  {
    //
    // Return the "raw" spectrum.
    //

    return Pxy/U;

  else

    //
    // Select 1st half and eliminate DC value.
    // Return a matrix with the frequency scale in the 1st column
    // and the PSD in the 2nd.
    //

    Pxy = Pxy[2:(M/2)+1]/U;	// Select 1st half and eliminate
				// DC value

    return [fscale[2:(M/2)+1], Pxy];
  }
};
