//-------------------------------------------------------------------//

//  Syntax:	logm ( A )

//  Description:

//  The logm function computes the matrix natural logarithm. Logm uses
//  Parlett's method in the rfile funm.r.

//  Another way to compute logm (A) is to:

//  logm ( X ) = e.vec * log (e.val) / e.vec

//  See Also: expm, funm
//-------------------------------------------------------------------//

//  Dependencies

require funm

logm = function ( X )
{
  global (log)
  return funm (X, log)
};
