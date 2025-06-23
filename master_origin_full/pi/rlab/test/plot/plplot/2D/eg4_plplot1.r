//--------------------------------------------------------------------------
// converted from
//    http://plplot.sourceforge.net/examples.php?demo=04&lbind=C
// Illustration of logarithmic scale and two y-axes
//--------------------------------------------------------------------------
if (isfile("plplot.r"))
{
  rfile plplot
}
if (length(plwins().available)<1)
{ plwins (1); }

// Set up data for log plot
N = 100;
exp_fmax=5;

f0 = 1.0; // corner frequency

freql = -2.0 + exp_fmax .* [0:N]' ./ N;
freq  = 10 .^ freql;
ampl  = 20.0 .* log10( 1.0 ./ sqrt(1.0 + (freq ./ f0 ).^2) );
phase = -rad * atan2(freq,f0);


//
//
//
plwin (1);
plaspect  (0.8);
plimits  (1e-2,1000,-60,0);
plimits2 (,,-90,0);
plscale   ("log");
plscale2  ("log");
plegend   (["Amplitude", "Phase"], 0.75, "rit", [0.05,0.05]);
plxlabel  ("Frequency (Hz)");
plylabel  ("Amplitude (dB)");
ply2label ("Phase shift (deg)");
plgrid    (2,0);
plformat ([...
    "with lines lt 1 lw 3 lc rgb blue using 1:2 axes x1y1", ...
    "with lines lt 1 lw 3 lc rgb green using 1:3 axes x1y2", ...
    []]);
pltext("-20 dB/decade", [50, -43], 1.25, <<color=10;incl=[1000,-90]>>);
plplot([freq,ampl,phase]);
_plprint ("eg4.eps", "psc");
