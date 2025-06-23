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

x = [1:10:1/16]';
y = sin(pi *x);

for (i in 1:NITER)
{
  spinner();

  xwid = 200 + 600 * int(10*urandom()) / 10;
  ywid = 100 + 400 * int(10*urandom()) / 10;

  if(_rlab_config.plot_support == "plplot")
  {
    xoff = 600 * int(10*urandom()) / 10;
    yoff = 400 * int(10*urandom()) / 10;
    plwins(1,"xwin",[xwid,ywid],[xoff,yoff]);
  else
    plwins(1,"/xwin",[xwid,ywid]);
    if (i==2)
    { break; }
  }

//   plwin (1);
  plwid(1);
  plegend ("Just another function y=y(x)", ceil(2 * uniform(),<<bin=0.5>>), "itr");
  plxlabel  ("x-axis");
  plylabel  ("y-axis");
  plformat  ("with lines lt 1 lc 1 lw 4 using 1:2");
  plplot    ([x,y]);
  if (NITER > 1)
  {
    sleep (1);
    plclose ();
  }
}
_plprint ("eg1.eps", "psc");
