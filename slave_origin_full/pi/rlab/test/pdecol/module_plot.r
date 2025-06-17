//
// file: module_plot.r
//
// plotting rutine for examples. At the end of execution type
// rfile module_plot
// to see the progress of calculation. changing the parameter
// of sleep command determines the speed.

plwins(1);

leg=blank(0,0);
for (i in 1:npde)
{
  leg = [leg, "u("+text(i)+")"];
}

tstr = "Fig.1: " + PDESOLVER + " for Example No. " + text(exno);

if (!isempty(wn))
{
  plwin(1);
  plegend( leg );
  xlabel ( "x" );
  ylabel ( "u(x)" );
  plimits(xlo, xhi, wmin, wmax);
  pltitle( tstr );
  for (i in 1:(wn.nc-npde+1):npde)
  {
    plot( [x, wn[;i:(i+npde-1)]] );
    sleep (0.1);
  }
}

