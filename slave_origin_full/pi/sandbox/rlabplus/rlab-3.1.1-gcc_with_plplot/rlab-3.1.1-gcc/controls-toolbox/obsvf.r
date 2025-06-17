//----------------------------------------------------------------------
//
// obsvf
//
// syntax:  </abar;bbar;cbar;k;t/> = obsvf(a,b,c)
//          </abar;bbar;cbar;k;t/> = obsvf(a,b,c,tol)
//
//	Observability staircase form.
//	</ABAR;BBAR;CBAR;K;T/> = obsvf(A,B,C) returns a decomposition
//	into the observable/unobservable subspaces.
//	</ABAR;BBAR;CBAR;K;T/> = obsvf(A,B,C,TOL) uses tolerance TOL.
//
//	If Ob=obsv(A,C) has rank r <= n, then there is a similarity
//	transformation T such that
//
//      Abar = T * A * T' ,  Bbar = T * B  ,  Cbar = C * T' .
//
//	and the transformed system has the form
//
//	      | Ano   A12|           |Bno|
//	Abar =  ----------  ,  Bbar =  ---  ,  Cbar = [ 0 | Co].
//	      |  0    Ao |           |Bo |
//
//	                                           -1           -1
//	where (Ao,Bo) is controllable, and Co(sI-Ao) Bo = C(sI-A) B.
//
//	See also: obsv
//
//----------------------------------------------------------------------
require ctrbf

obsvf = function(a,b,c,tol)
{
  // Use ctrbf and duality:

  if (nargs == 4) {
     tmp = ctrbf(a',c',b',tol);
  else
     tmp = ctrbf(a',c',b');
  }
  return <<abar=tmp.abar'; bbar=tmp.bbar'; cbar=tmp.cbar';t=tmp.t;k=tmp.k>>;
};

