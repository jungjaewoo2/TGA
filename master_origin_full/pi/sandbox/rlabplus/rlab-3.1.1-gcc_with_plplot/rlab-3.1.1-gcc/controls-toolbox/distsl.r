//--------------------------------------------------------------------------
// distsl
//
// Syntax: </xp;yp/> = distsl(x1,y1,a,b)
//	
//      returns the distance of a point (x1,y1) 
//	to a straight line of the form  y=ax+b
//
//	When x1 and y1 are vectors then a vector of points is returned. 
//--------------------------------------------------------------------------
distsl = function(x1,y1,a,b)
{
  yn = a.*x1+b;
  xn = (y1-b)./a;
  yp = (a.*(-yn.*y1 + y1.^2 + xn.*x1 - x1.^2) + b.*(xn-x1))./(xn- x1 + a.*y1 - a.*yn);
  xp = (yp-b)./a;
  xp = xp[:];
  yp = yp[:];
  ind = find(yp==nan());
  yp[ind] = zeros(length(ind),1);
  ind = find(xp==nan());
  xp[ind] = zeros(length(ind),1);
  return <<xp=xp;yp=yp>>;
};

