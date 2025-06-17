//-------------------------------------------------------------------//
// lyap.r:
//
// Syntax:	lyap ( A , B , C )
//		lyap ( A ,   , C )
//
// Description:
//
// lyap solves the general form of the Lyapunov (Sylvester) equation:
//
//	A*X + X*B = -C
//
// To solve the special form of the Lyapunov equation:
//
//	A*X + X*A' = -C
//
// Skip the second argument when using lyap (`lyap (A,,C)').

// See Also: schur, sylv
//-------------------------------------------------------------------//

lyap = function ( A, B, C )
{
  if (!exist (B)) 
  { 
    B = A';	/* Solve the special form: A*X + X*A' = -C */
  }

  if ((A.nr != A.nc) || (B.nr != B.nc) || (C.nr != A.nr) || (C.nc != B.nr)) 
  {
    error ("Dimensions do not agree.");
  }

  //
  // Schur decomposition on A and B
  //

  sa = schur (A);
  sb = schur (B);

  //
  // transform C
  //

  tc = sa.z' * C * sb.z;

  X = sylv (sa.t, sb.t, tc);

  //
  // Undo the transformation
  //

  X = sa.z * X * sb.z';

  return X;
};
