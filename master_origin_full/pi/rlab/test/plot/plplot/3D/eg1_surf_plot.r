//
// file: main6.r
//
if (isfile("plplot.r"))
{
  rfile plplot
}
if (length(plwins().available)<1)
{ plwins (1); }
if (!exist(NITER))
{ NITER = 1; }

xmin = 0;
xmax = 1;
x = [xmin:xmax:1/16]';

ymin = -1;
ymax = 1;
y = [ymin:ymax:1/16]';

z = zeros(length(x),length(y));
for (i in 1:length(x))
{
  for (j in 1:length(y))
  { z[i;j] = 2 + x[i].^2 - y[j].^2; }
}

plwin (1);
plxlabel("x");
plylabel("y");
plzlabel("F(x,y)");
plimits (xmin,xmax,ymin,ymax);
plxtics (1/2, 5);
plytics (1/2, 5);
plztics (1, 5);
plgrid3 (,,"bnstu");
// plztics  (1/2, 5);
plegend ( ["Surface x^2+2*x+y^2"] );
plformat("with lines lw 1 lt 1 lc rgb red draw xy mag");
// plformat ("with points pt 1 ps 1 pc rgb red");
data = <<>>;
  // first data set: from rlab
  data.x = y;
  data.y = x;
  data.z = z';
plplot3(x,y,z);

