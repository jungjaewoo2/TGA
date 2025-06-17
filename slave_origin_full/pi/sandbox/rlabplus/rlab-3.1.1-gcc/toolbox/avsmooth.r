//-------------------------------------------------------------------//

// Synopsis: Average smoothing filter.

// Syntax:   Y = avsmooth ( X , L ) 

// Description:

//      avsmooth smooths the input vector X using a median filter with
//      a rectangular window of "L" samples.

//      See also: steng, stmag, stzcr, mdsmooth

//      Matlab original by:
//       LT Dennis W. Brown 7-11-93, DWB 8-17-93
//       Naval Postgraduate School, Monterey, CA
//       May be freely distributed.
//       Not for use in commercial products.

//-------------------------------------------------------------------//

require mean

avsmooth = function ( x , L )
{
  local ( x , L )

  // default output
  y = [];

  // check args
  if (nargs != 2)
  {
    error("avsmooth: Invalid number of input arguments...");
  }

  // figure out if we have a vector
  if (min(size(x)) != 1)
  {
    error("avsmooth: Input arg \"x\" must be a 1xN or Nx1 vector.");
  }

  // work with Nx1 vectors
  x = x[:];

  // number of samples
  Ns = length(x);

  // room for output
  y = zeros(Ns,1);

  // pad ends to compensate for filter length
  x = [zeros(L/2,1); x ; 0 ; zeros(L/2+1,1)];

  // average filter
  for (k in 1:Ns)
  {
    y[k;1] = mean(x[k:k+L-1;1]);
  }

  return y;
};
