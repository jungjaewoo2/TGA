//-------------------------------------------------------------------//

//  Syntax:	detrend ( A )

//  Description:

//  Detrend removes a linear trend from a vector.

//  y = detrend (x) 

//  removes the best straight-line fit from the data in vector x and
//  returns it in vector y. Detrend is often used prior to fft().

//-------------------------------------------------------------------//

detrend = function ( x )
{
  local (a, m, n, y);

  m = x.nr; n = x.nc;
  m = max ( [m, n] );
  y = reshape (x, prod (size (x)), 1);
  a = [ (1:m)'/m, ones(m,1) ];
  y = y - a*(a\y);
  if (n != 1)
  {
    y = conj (y');
  }
  return y;
};
