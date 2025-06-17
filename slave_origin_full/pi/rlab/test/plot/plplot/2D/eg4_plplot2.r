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

_plot_data = <<>>;
_plot_data.a = [freq,ampl];
_plot_data.b = [freq,phase];


//
//
//
plwin (1);
plaspect  (0.75);
plimits  (1e-2,1000,-60,0);
plimits2 (,,-90,0);
plscale   ("log");
plscale2  ("log");
plxlabel  ("Frequency (Hz)");
plylabel  ("Amplitude (dB)");
ply2label ("Phase shift (deg)");
plgrid    (2,0);
plformat ([...
    "with lines lt 1 lw 3 lc rgb blue axes x1y1", ...
    "with lines lt 1 lw 3 lc rgb green axes x1y2", ...
    []]);
plplot(_plot_data);

