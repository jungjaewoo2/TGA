//-------------------------------------------------------------------//
//  Syntax:  log2 ( X )
//  Description:
//  RLaB equivalent for log2
// [azg] Tue Jul  1 19:35:03 PDT 2008
//  STATUS: New, untested
//

//-------------------------------------------------------------------//



log2 = function(x) 
{
  // Checking the class will automatically result in an error if
  // the argument doesn't exist.
  if( class(x) != "num") {
    error("erfc: argument is non-numeric.");
  }

return log(x)/log(2.0);
}
