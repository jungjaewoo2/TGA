#
# plhold_demo.r
# A simple demonstration of plhold().
#

require ode4

#
# Start a normal plot window.
#

  #  plstart(,,);

#
# Set the limits cause plhold may not have enough
# information to do so itself.
#

  plimits (-3,3,-3,3);

#
# Take care of miscelaneous items
#

  pltitle ("Plhold Demo");
  xlabel ("X");
  ylabel ("Xdot");
  plegend ("Phase-Plane Trajectory");

#
# Set up for ode() usage.
#

  vdpol = function ( t , x ) 
  {
    xdot[1;1] = x[1] * (1 - x[2]^2) - x[2];
    xdot[2;1] = x[1];
    return xdot;
  };


  t0 = 0; dt = 0.1; x0 = [0.3; 0.25];

#
# Now go into a loop plotting the phase plane as we integrate.
# This is not a very efficient way to integrate, but it is kind
# of fun.
#

  sout = [];
  for (i in 1:100)
  {
    plhold ( (out = ode4 (vdpol, t0, tf=t0+dt, x0))[;2,3] );
    t0 = tf;
    x0 = out[out.nr;2,3];
    sout = [sout; out];
  }
  plptex ("STOP", out[out.nr;2], out[out.nr;3]);

#
# Calling plhold_off is VERY important, plotting will
# not work afterwords otherwise.
#

  plhold_off ();
  plegend ("default");
  plimits ();
