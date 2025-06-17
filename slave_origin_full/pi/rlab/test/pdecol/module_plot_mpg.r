//
// file: module_plot.r
//
// plotting rutine for examples. At the end of execution type
// rfile module_plot
// to see the progress of calculation. changing the parameter
// of sleep command determines the speed.

plwins(1);

leg = ["Displacement", "Velocity"];

tstr = "Driving the waves on a rope";

if (!isempty(wn))
{
  plwin(1);
  plegend( leg );
  xlabel ( "x" );
  ylabel ( "u(x)" );
  plimits(xlo, xhi, wmin, wmax);
  pltitle( tstr );
  j = 0;
  for (i in 1:(wn.nc-npde+1):npde)
  {
    data = [x, wn[;i:(i+npde-1)]];
    // plot on screen
    plot( data );
    // make a ppm frame of it
    if (i==1)
    {
      "your version of pgplot does not support ppm-output"
    }
    if (0)
    {
      psfname="./frames/"+gsub("0", " ",text(j,"%4g")).string+".ppm";
      plcopy(1,psfname+"/ppm");
      plegend(["u(x,t)", "u_t(x,t)"] );
      plwid([5,5,5,5]);
      plot( data );
      plclose();
      j = j + 1;
    }
    sleep (0.1);
  }
}

