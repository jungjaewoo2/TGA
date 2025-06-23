//-------------------------------------------------------------------//
// surspl.r
 
// Syntax:      surspl ( xi, yi, xj, yj, wj )
 
// Description:
 
//      Surface Spline Interpolation:

//	xi: x coordinates of points to interpolate to (vector).
//	yi: y coordinates of points to interpolate to (vector).
//      xj: known x coordinates (vector).
//      yj: known y coordinates (vector).
//	wj: known z coordinates of w[xj;yj] = f(xj,yj)

//	Return interpolated w[xi; yi] 
 
//      Eric Haugse, 8-24-94
 
//      References:
//      Robert L. Harder and Robert N. Desmarais, "Interpolation 
//      using Surface Splines"

//-------------------------------------------------------------------//

static (expand)
static (trans)

surspl = function (xi, yi, xj, yj, wj)
{
  // Make sure they are all column vectors...
  xi = xi[:]; yi = yi[:];
  xj = xj[:]; yj = yj[:];

  // Simple error checks...
  if (length (xj) != wj.nr) {
    error ("surspl: length (xj) must equal wj.nr");
  }
  if (length (yj) != wj.nc) {
    error ("surspl: length (yj) must equal wj.nc");
  }  

  // Get on with it...

  crdj = expand(xj,yj);
  Ljj = trans(crdj.x,crdj.y,crdj.x,crdj.y);
  A = [[zeros(3,3),Ljj.crd']; Ljj.tr];
  B = [zeros(3,1); wj[:]];
  coeff = solve(A,B,"S");
  
  crdi = expand(xi,yi);
  Lij = trans(crdi.x,crdi.y,crdj.x,crdj.y);
  return reshape(Lij.tr*coeff,xi.nr,yi.nr);
};

expand = function(x,y)
{
  nx = x.nr; ny = y.nr;
  y1 = zeros (nx*ny, 1);

  for (j in 1:ny)
  {
    if (j == 1)
    {
      x1 = x;
    } else {
      x1 = [x1;x];
    }
    for (i in 1:nx)
    {
      y1[(j-1)*nx+i] = y[j];
    }
  }
  return <<x = x1; y = y1>>;   
};

trans = function(xa,ya,xb,yb)
{
  na = xa.nr; 
  nb = xb.nr;

  k = zeros(na,nb);
  nones = ones(na,1);
  
  for (i in 1:na)
  {
    for (j in 1:nb)
    { 
      rsq = (xa[i] - xb[j])^2 + (ya[i] - yb[j])^2;
      if(rsq == 0)
      {
	k[i;j] = 0.0;
      } else {
	k[i;j] = rsq * log(rsq);
      }
    }
  }
  
  coord = [nones,xa,ya];
  transf = [coord,k];
  
  return <<crd=coord; tr=transf>>;
};
