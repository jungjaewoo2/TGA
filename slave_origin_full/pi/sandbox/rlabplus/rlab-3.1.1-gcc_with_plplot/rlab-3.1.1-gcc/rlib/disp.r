//-------------------------------------------------------------------//

//  Synopsis:   Display a matrix.

//  Syntax:	disp ( S )

//  Description:

//  The disp function does two things,

//  1) If the argument is a matrix, it prints the matrix without the
//  variable label to the standard output.

//  2) If the argument is a string, it prints the string to the
//  standard output.

//-------------------------------------------------------------------//

disp = function ( S )
{
  if (class (S) == "string")
  {
    fprintf ("stderr", "%s\n", S[1]);
  else if (class (S) == "num") {
    S + 0 ?
  }}

  return 1;
};
