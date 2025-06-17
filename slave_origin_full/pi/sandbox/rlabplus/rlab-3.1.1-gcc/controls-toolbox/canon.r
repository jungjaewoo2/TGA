//--------------------------------------------------------------------------
//
// canon
//
// Syntax: </ab;bb;cb;db;T/> = canon(a,b,c,d,Type)
//
//	State-space to canonical form transformation.
//	canon transforms the continuous state-space system (A,B,C,D) 
//      into the canonical form specified by
//	type = "modal" transforms the state-space system into modal form 
//	                where the system eigenvalues appear on the 
//	                diagonal.  The system must be diagonalizable.
//
//	type = "companion" transforms the state-space system into 
//	                companion canonical form where the characteristic
//	                polynomial appears in the right column.
//
//	the transformation matrix T is returned where z = Tx.
//
//	The modal form is useful for determining the relative controll-
//	abilty of the system modes.  Note: the companion form is ill-
//	conditioned and should be avoided if possible.
//
//	See also: ss2ss, ctrb, and ctrbf.
//--------------------------------------------------------------------------
require abcdchk isstr

canon = function(a,b,c,d,Type)
{
  if (nargs != 4 && nargs != 5) { error("cannon: need 4 or 5 arguments"); }  
  msg = abcdchk(a,b,c,d);
  if (msg != "") { error(msg);}  
  if (nargs==4)
  { 
    // No type specified, assume modal 
    Type = "modal";
  }
  
  if (!isstr(Type))
  {
     error("Type must be a string."); 
  }

  // Determine 'type' 
  if (Type == "modal")  
  { //Modal form
    eigen  = eig(a);
    n = length(eigen.val);
    k = 1;
    // Transformation to modal form based on eigenvectors
    while (k<=n)
    {
      if (imag(eigen.val[k]) != 0.0)
      {
        T[;k]   = real(eigen.vec[;k]); 
        T[;k+1] = imag(eigen.vec[;k]);
        k = k + 2;
      else
        T[;k] = eigen.vec[;k];
        k = k+1;
      }
    }
    T = real(T); 
    LU = factor(T);
    ab = backsub(LU,a)*T; bb = backsub(LU,b); cb = c*T; db = d;
  else if (Type == "companion") 
  { // Companion form
    // Transformation to companion form based on controllability matrix
    if (length(b))
    {
       T = ctrb(a,b[;1]);
    else
       T = [];
    }
    LU = factor(T);
    ab = backsub(LU,a)*T; bb = backsub(LU,b); cb = c*T; db = d;
  else
    error("TYPE must be either 'modal' or 'companion'.");
  }}

  // Return inverse of T to be compatible with ss2ss
  if (T!=[]) {T = inv(T);}

  return <<a=ab;b=bb;c=cb;d=db;T=T>>;
};

