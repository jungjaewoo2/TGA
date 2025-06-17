//-------------------------------------------------------------------//

// Synopsis:	Decimal to Hexadecimal number conversion.

// Syntax:	dec2hex ( HEX )

// Description:

//      Dec2hex takes a numeric argument, and converts it to a 
//      string containing the hexa-decimal representation of the
//      input argument.

// See Also: strtod, strtol, hex2dec
// [azg] Thu Jun 19 16:32:49 PDT 2008

//-------------------------------------------------------------------//

dec2hex = function ( DEC )
{
local(i, tmp)
  # Error checking...
  if (class (DEC) != "num") {
    error ("dec2hex: argument must be class num");
  }

for(i in 1:size(DEC)[1])
          {
            for(j in 1:size(DEC)[2])
            {
              sprintf (tmp, "%x", DEC[i;j]);
              x[i;j] = tmp;
            }
          }


  return x;
};
