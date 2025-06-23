//----------------------------------------------------------------------
//
// givens
//
// syntax: g = givens(x,y)
//      Givens rotation matrix.
//	G = givens(x,y) returns the complex Givens rotation matrix
//
//	    | c       s |                  | x |     | r | 
//	G = |           |   such that  G * |   |  =  |   |
//          |-conj(s) c |                  | y |     | 0 |
//	                                
//	where c is real, s is complex, and c^2 + |s|^2 = 1. 
//
//----------------------------------------------------------------------
givens = function(x,y)
{
  absx = abs(x);
  if (absx == 0.0)
  {
	c = 0.0; s = 1.0;
  } else {
	nrm = norm([x,y],"2");
	c = absx/nrm;
	s = x/absx*(conj(y)/nrm);
  }
  return [c, s;-conj(s), c];
};
