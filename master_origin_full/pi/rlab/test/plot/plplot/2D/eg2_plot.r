//
//
//
if (isfile("plplot.r"))
{
  rfile plplot
}
if (length(plwins().available)<1)
{ plwins (1); }
if (!exist(NITER))
{ NITER = 1; }

xwid = 200 + 600;
ywid = 100 + 400;

xoff = 600 * int(10*urandom()) / 10;
yoff = 400 * int(10*urandom()) / 10;
plwins(1,"xwin",[xwid,ywid],[xoff,yoff]);

for (i in 1:NITER)
{
  spinner();

  s = 1 + int(10*urandom()) / 10;
  x = [1:10:1/16]';
  y = sin(s * pi *x);
  z = cos(s * pi *x).^2;
  q = x .* cos(s * pi *x);

  plwin(1);
  plegend (["Just another function y=y(x)", "and another"], 0.75, "rits", [0.09,0.15]);
  plxlabel  ("x-axis");
  plylabel  ("y-axis");
  plformat  (["with lines lt 1 lc rgb blue using 1:2 axes x1y1", "with lines lt 1 lc rgb red using 1:3"]);
  plplot    ([x,y,z,q]);

  if (NITER > 1)
  { sleep (1); }
}
_plprint ("eg2.eps", "psc");
