//-------------------------------------------------------------------//

//  Synopsis:   Compute the trace of a matrix.

//  Syntax:	trace ( A )

//  Description:

//  Compute the trace (the sum of the diagonal elements) for the input
//  matrix. `trace(A)' is the same as `sum( diag( A ) )'.

//-------------------------------------------------------------------//

trace = function(m) 
{
  if(m.class != "num") 
  { 
    error("must provide NUMERICAL input to trace()");
  }

  tr = 0;
  for(i in 1:min( [m.nr, m.nc] )) 
  {
    tr = tr + m[i;i];
  }
  
  return tr;
};
