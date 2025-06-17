//--------------------------------------------------------------------------
//
// bode_plot
//
// syntax: bode_plot(R)
//
// Plot Bode magnitude and phase angle diagrams (vs. frequency)
//
// R.mag   = magnitude of gain
// R.phase = phase angle        (in degrees)
// R.w     = frequency          (in rad/sec)
//
// calling sequence:
//
// R = bode(a,b,c,d,iu,w);  (or R = bode(num,den);)
// bode_plot(R);  
//
//--------------------------------------------------------------------------

bode_plot = function(R)
{
  // Bode plot
  // need to plot 2 graphs (top & bot) in one window

  global(eps,pi,_rlab_config,_bode_window)
  
  if (R.mag.nc != R.phase.nc || R.mag.nr != R.phase.nr) {
     error("Mag/Phase have compatibility problem.");
  }

  // check if a plot window is available 
  if (_rlab_config.plot_support=="pgplot"||_rlab_config.plot_support=="plplot")
  {
     if (exist(_bode_window)) {
        plwin(_bode_window);
     } else {
        _bode_window = plstart(1,2);
     }
  }
  if (_rlab_config.plot_support=="gnuplot")
  {
     multiplot (2,1);     
  }
  // plot mag diagram
  leg = [];
  for (i in 1:R.mag.nc) { sprintf(s,"%d",i); leg=[leg;s]; }
  leg1 = [leg;"0dB"];
  plegend(leg1);
  plaxis("log");
  pltitle("Bode Magnitude");
  xlabel("Frequency (rad/sec)");
  ylabel("Gain dB");
  plot(<<[R.w,20*log10(R.mag)];[R.w,zeros(R.w.n,1)]>>);


  // 180 degree line
  if (min(min(R.phase)) < -180) { 
    phase180 = -180;
  } else { if (max(max(R.phase)) < 180) {
    phase180 = -180;
  } else {
    phase180 = 180;
  }}

  // plot phase diagram
  sprintf(s,"%d",phase180);
  leg2 = [leg;s];
  plegend(leg2);
  plaxis("log");    
  pltitle("Bode Phase");
  xlabel("Frequency (rad/sec)");
  ylabel("Phase deg");
  plot(<<[R.w,R.phase];[R.w,phase180*ones(R.w.n,1)]>>);
  plegend("default");
  plaxis();
  pltitle("");
  xlabel("");
  ylabel("");
  if (_rlab_config.plot_support=="gnuplot")
  {
     nomultiplot ();     
  }

};



