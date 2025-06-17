//------------------------------------------------------------------------
//
// dbode
//
// Bode plot for discrete-time linear systems
// 
// calling sequence:
//
// R = dbode(a,b,c,d,ts,iu,w);
// dbode_plot(R);
//  
// R.mag   = magnitude of gain
// R.phase = phase angle        (in degrees)
// R.w     = frequency          (in rad/sec)
//
//------------------------------------------------------------------------
dbode_plot = function(R)
{
  require   bode_plot
  
  if (R.mag.nc != R.phase.nc || R.mag.nr != R.phase.nr) {
     error("Mag/Phase have compatibility problem.");
  }
  bode_plot(R);

};
