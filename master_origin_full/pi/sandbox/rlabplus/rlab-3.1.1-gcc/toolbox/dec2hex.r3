//-------------------------------------------------------------------//

// Synopsis:	Decimal to Hexadecimal number conversion.

// Syntax:	dec2hex ( HEX )

// Description:

//      Dec2hex takes a numeric argument, and converts it to a 
//      string containing the hexa-decimal representation of the
//      input argument.

// See Also: strtod, strtol, hex2dec

//-------------------------------------------------------------------//

dec2hex = function ( DEC )
{
  # Error checking...
  if (class (DEC) != "num") {
    error ("dec2hex: argument must be class num");
  }

  sprintf (tmp, "%x", DEC);

  return tmp;
};
