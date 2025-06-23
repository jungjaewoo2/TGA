//-------------------------------------------------------------------//

// Synopsis:	Compute Auto-Spectra of the input sequence

// Syntax:	autospec ( X, dt, M, WIN, raw );

// Description:

//      Compute the auto-spectra (psd) of the input vector X. autospec
//      returns a two-column matrix. The 1st column is the frequency
//      scale and the second column is the values of the one-sided
//      spectra.

//	X:	Input sequence
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

autospec = function ( X , dt, M, WIN, raw )
{
  X = X[:];		// Make sure it is a row

  n = length (X);
  if (!isreal (X)) { error ("autospec: X must be real"); }

  //
  // Set the defaults.
  //

  if (!exist (dt)) { dt = 1; }
  if (!exist (M)) { M = n; }
  if (!exist (WIN)) { WIN = "rect"; }
  if (!exist (raw)) { raw = 0; }

  // Some error checking

  if (M > n) { error ("autospec: M > n not allowed"); }
  if (M < 0) { error ("autospec: M < 0 not allowed"); }

  // Create the frequncy scale (in Hertz) up to the Nyquist

  fscale = (((0:M-1)/M)/dt)';

  k = fix (n/M);                // The number of windows
  index = 1:M;

  // Create the window

  w = window (M, WIN);          // Window specification
  U = k*sqrt(sum(w.^2)/w.n);	// Normalizing scale factor

  Pxx = zeros (M, 1);

  for (i in 1:k)
  {
    xw = w .* detrend (X[index]);
    index = index + M;
    Xx = (2/(M*dt)) * abs (dt * fft (xw)) .^ 2;
    Pxx = Pxx + Xx;
  }

  if (raw)
  {

    //
    // Return the "raw" spectrum.
    //

    return Pxx/U;
  else

    //
    // Select 1st half and eliminate DC value.
    // Return a matrix with the frequency scale in the 1st column
    // and the PSD in the 2nd.
    //

    Pxx = Pxx[2:(M/2)+1]/U;
    return [fscale[2:(M/2)+1], Pxx];
  }
};
