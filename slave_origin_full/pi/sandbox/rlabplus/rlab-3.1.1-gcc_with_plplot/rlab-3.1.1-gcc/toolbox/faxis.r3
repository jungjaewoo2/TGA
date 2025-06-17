//---------------------------------------------------------------------------
//  faxis.r

//  Syntax:	faxis ( X )
//		faxis ( X , T )
//		faxis ( X , T, axis_type )

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

faxis = function ( X, _T, axis_type )
{
  global (pi)

  if (!exist (axis_type)) { axis_type = 1; }
  if (!exist (_T)) { T = 1; axis_type = 1; } else { T = _T; } 

  N = length (X);
  a = (0:N-1)/N;
  a = reshape (a, X.nr, X.nc);

  if (axis_type == 3)
  {
    a = a/T;
  }
  else if (axis_type == 2)
  {
    a = a * (2*pi/T);
  }
  else if (axis_type == 4)
  {
    a = 2 * a;
  }
  else
  {
    a = a * (2*pi);
  }

  return a
};
