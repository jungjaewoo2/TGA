//-----------------------------------------------------------------------
//
// dimpulse_plot
//
// syntax: dimpulse_plot(R)
//
// plot impulse response of discrete-time linear system
//
// see also: dimpulse
//-----------------------------------------------------------------------
require stairs

dimpulse_plot = function(R)
{
  pltitle("Impulse response of discrete-time linear system");
  xlabel("No. of Samples"); 
  ylabel("Amplitude");
  stairs([0:R.n-1],R.y);
  
};


