//-------------------------------------------------------------------//

//  Synopsis:   Calculate the phase angle of a complex value.

//  Syntax:     angle ( Z )

//  Description:

//  Caclulate the phase angle of a complex value.

//-------------------------------------------------------------------//

angle = function(a)
{
  return atan2(imag(a), real(a));
};
