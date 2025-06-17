//----------------------------------------------------------------------
//
// nargchk
//
// Check number of input arguments. Return error message if
// not between low and high.  If it is, return empty matrix.
//
//----------------------------------------------------------------------
nargchk = function(low,high,number)
{
  msg="";
  if (number < low) {
     msg = "Not enough input arguments.";
  else if (number > high) {
     msg = "Too many input arguments.";
  }}
  return msg;
};

