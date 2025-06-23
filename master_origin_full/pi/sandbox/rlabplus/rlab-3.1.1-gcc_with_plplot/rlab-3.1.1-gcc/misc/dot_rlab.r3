#
# .rlab: Rlab startup script.
# Put whatever you like in here, and it will be run
# automatically when you start Rlab.
#

# Calculate pi

// if(exist(const.pi))
// {
//   pi = const.pi;
// else
//   pi = atan(1.0)*4.0;
// }

# Calculate machine epsilon

epsilon = function()
{
  eps = 1.0;
  while((eps + 1.0) != 1.0)
  {
    eps = eps/2.0;
  }
  return 2*eps;
};

pi = const.pi;
eps = epsilon();
deg = pi/180;
rad = 180/pi;

#
# Rlab2 doesn't have read/write, so make up for it:
#
read = readb;
write = writeb;
variance=var;
covariance=covar;
urandom=uniform;
text=num2str;
set=unique;
complement=setdiff;

format(9,9);

#
# spice:
#
// if (exist(spice))
// {
//   if (class(spice)=="list")
//   {
//     rfile ngspice
//   }
// }

#
# end of startup script
#

