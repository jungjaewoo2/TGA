//
//
//
if (isfile("plplot.r"))
{
  rfile plplot
}
if (length(plwins().available)<1)
{ plwins (1); }
NITER = 1;

xwid = 200 + 600;// * int(10*urandom()) / 10;
ywid = 100 + 400;// * int(10*urandom()) / 10;

if(_rlab_config.plot_support == "plplot")
{
  xoff = 600 * int(10*urandom()) / 10;
  yoff = 400 * int(10*urandom()) / 10;
  plwins(1,"xwin",[xwid,ywid],[xoff,yoff]);
else
  plwins(1,"xwin",[xwid,ywid]);
}


for (i in 1:NITER)
{
  spinner();

  s = 1 + int(10*urandom()) / 10;
  x = [1:10:1/16]';
  y = sin(s * pi *x);

  plwin(1);
  plimits   (,,-1.1,1.5);
  plegend   ("Just another function y=y(x)", 0.75, "rtis",[0.1,0.16]);
  plxlabel  ("x-axis");
  plylabel  ("y-axis");
  plformat  (["with yerrorlines lt 1 lc 1"]);
  plplot    ([x,y,y-0.1,y+0.1]);

  if (NITER>1)
  { sleep (1); }
}
_plprint ("eg3.eps", "psc");
