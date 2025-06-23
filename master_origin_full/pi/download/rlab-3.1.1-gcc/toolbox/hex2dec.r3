//-------------------------------------------------------------------//

// Synopsis:	Hexadecimal to decimal number conversion.

// Syntax:	hex2dec ( HEX )

// Description:

//      Hex2dec takes a string argument, and converts it to a decimal
//      number. 

// See Also: strtod, strtol

//-------------------------------------------------------------------//

hex2dec = function ( HEX )
{
  # Error checking...
  if (class (HEX) != "string") {
    error ("dec2hex: argument must be of class string");
  }

  return strtol (HEX, 16);
};
