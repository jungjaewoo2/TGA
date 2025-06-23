// ---------------------------------------------------------------------
// perpxy
//
// syntax: </xp;yp/> = perpxy(x1,y1,x2,y2,x3,y3)
//
// peryxyFinds a point [xp,yp] which is the nearest point from
// [x3,y3] on the straight line formed between [x1,x2] and [x2,y2].
//
// </xp;yp/> = perpxy(x1,y1,x2,y2,x3,y3)
//
// Returns points [x4,y4] which are perpendicular to the
// the straight line formed between the points
// points [x1,y1] and [x2,y2], starting from the points
// [x3,y3]:
//
// i.e. (x3-x4)*(x1-x2) + (y3-y4)*(y1-y2) = 0
//      (x1-x4)*(y4-y2) - (x4-x2)*(y1-y4) = 0
//
//----------------------------------------------------------------------

perpxy = function (x1,y1,x2,y2,x3,y3)
{

  global (eps)

  xeq= abs(x1-x2) < 1e-10;
  x12= (x1-x2) + eps*(xeq);
  y12= (y1-y2);
  x4 = (x12.*(-x3.*x12-y3.*y12) + y12.*(x1.*y2-x2.*y1))./(-x12.*x12 - y12.*y12);
  y4 = (x1.*y2 + x4.*y12 - x2.*y1)./(x12);
  y4 = y4.*(!xeq) + y3.*(xeq);
  return <<x4=x4; y4=y4>>;
};



