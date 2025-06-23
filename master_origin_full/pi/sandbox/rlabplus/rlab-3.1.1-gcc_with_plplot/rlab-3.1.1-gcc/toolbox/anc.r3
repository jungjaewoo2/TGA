//-------------------------------------------------------------------//

//  Syntax:     anc ( SN, REF, mup, nw, delay )

//  Description:

//  Anc performs adaptive noise cancellation as described in:

//  "Adaptive Noise Cancelling: Principles and Applications" by
//  B. Widrow, et. al. Published in Proceedings of the IEEE, Vol 63,
//  No. 12, December 1975.

//  Inputs:
//    SN:    A vector of the signal + noise.
//    REF:   The reference input (noise).
//    mup:   The adaption rate factor (percentage of max).
//    nw:    The number of weights in the FIR filter.
//    delay: The number of samples to delay the primary
//           input (SN) by.

//  Returns:
//    e:     A vector of the error signal (the "cancelled" signal).
//    y:     A vector of the filter output.
//    w:     A matrix of the filter weights.

//  Ian Searle (5/25/95)
//  With guidance from LT Dennis W. Brown's SPC Toolbox

//  Dependencies
    require std

//-------------------------------------------------------------------//

anc = function ( SN, REF, mup, nw, delay )
{
  # Error checking, etc...
  if (!exist (mup)) { mup = 100; }
  if (!exist (nw)) { error ("anc: must supply argument NW"); }
  if (!exist (delay)) { delay = 0; }

  SN = SN[:]; REF = REF[:];
  n = SN.n;

  # Compute maximum nu...
  runs = std (REF)^2;
  printf ("mumax = %f\n", mumax = 2/(nw * runs));
  mu = mup/100 * mumax

  # Set things up...
  y = zeros (size(SN));  # The filter output
  e = zeros (size(SN));  # The error signal
  wold = zeros (nw,n);
  w = zeros (nw,1);

  # Adjust the reference and the input signal for the filter...
  SN = [zeros (delay,1); SN];
  REF = [zeros (nw-1,1); REF];

  # Loop over each element of the signal vector
  for (i in 1:n)
  {
    # 1st compute Y(k)
    y[i] = w' * REF[i+nw-1:i:-1];

    # Compute error: e(k) = SN(k) - Y(k)
    e[i] = SN[i] - y[i];

    # Update the filter weights (W(k+1))
    w = w + mu * e[i] * REF[i+nw-1:i:-1];

    # Save the old weights for information...
    wold[;i] = w;
  }

  return << e=e; y=y; w=wold >>;
};
