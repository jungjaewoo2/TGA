//
//
// library of useful functions for digital signal processing:
//    mostly from matlab, or my own
//

// poles and gain of an order n Butterworth analog lowpass filter prototype
buttap = function( n )
{
  global (const);
  if (!exist(n))
  { n = 0; }

  rval = <<>>;
  rval.zeros = [];

  if (n)
  {
    rval.poles = exp(1i .* (const.pi .* [1:2*n-1:2] ./ (2*n) + 0.5 * const.pi));
    rval.gain  = real(prod(-p));
  } else {
    rval.poles = [];
    rval.gain  = [];
  }

  return rval;
};
