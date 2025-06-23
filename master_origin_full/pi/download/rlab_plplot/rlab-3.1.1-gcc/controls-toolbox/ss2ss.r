//----------------------------------------------------------------------
//
// ss2ss
//
// Syntax: </at;bt;ct;dt/> = ss2ss(a,b,c,d,t)
//
//	</At;Bt;Ct;Dt/> = ss2ss(A,B,C,D,T) performs the similarity 
//	transform z = Tx.  The resulting state space system is:
//
//		.       -1        
//		z = [TAT  ] z + [TB] u
//		       -1
//		y = [CT   ] z + Du
//
//	See also: canon, balreal and balance.
//
//----------------------------------------------------------------------
require abcdchk

ss2ss = function(a,b,c,d,t)
{
  if (nargs != 5) {
     error("Wrong number of arguments");
  }
  error(abcdchk(a,b,c,d));

  at = t*a/t; 
  bt = t*b; 
  ct = c/t; 
  dt = d;

  return <<at=at; bt=bt; ct=ct; dt=dt>>;
};

