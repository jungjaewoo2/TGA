//-------------------------------------------------------------------//
//  Syntax:  exp2 ( X )
//  Description:
//  RLaB equivalent for 2^matrix
// [azg] Fri Jun 20 12:54:59 PDT 2008
//

//-------------------------------------------------------------------//



exp2 = function(x) 
{
  // Checking the class will automatically result in an error if
  // the argument doesn't exist.
  if( class(x) != "num") {
    error("exp2: argument is non-numeric.");
  }

return  exp(log(2) * x) ;  
}
