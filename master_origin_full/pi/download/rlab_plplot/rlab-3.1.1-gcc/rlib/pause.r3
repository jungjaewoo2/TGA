//-------------------------------------------------------------------//

//  Synopsis:   Pause execution.

//  Syntax:	pause ( )

//  Description:

//  The pause function stops execution of an RLaB program until the
//  return key is pressed.

//  See Also: fprintf, getline

//-------------------------------------------------------------------//

pause = function ( message )
{
  if (!exist (message))
  { 
    message = "\tHit return to continue"; 
  }

  fprintf ("stdout", "%s\n", message);
  getline ("stdin");
};
