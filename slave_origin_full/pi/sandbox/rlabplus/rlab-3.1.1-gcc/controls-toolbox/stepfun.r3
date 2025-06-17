//----------------------------------------------------------------------
//
// stepfun
//
// Syntax: y = stepfun(t,to)
//
// Unit step function.
// stepfun(T,T0), where T is a monotonically increasing vector,
// returns a vector the same length as T with zeros where T < T0
// and ones where T >= T0.
//
//----------------------------------------------------------------------

stepfun = function(t,to)
{
  m = t.nr;
  n = t.nc;
  y = zeros(m,n);
  i = find(t>=to);
  if (isempty(i)) { return y; }
  i = i[1];
  if (m == 1)
  {
	y[i:n] = ones(1,n-i+1);
  } else {
	y[i:m] = ones(m-i+1,1);
  }
  return y;
};


