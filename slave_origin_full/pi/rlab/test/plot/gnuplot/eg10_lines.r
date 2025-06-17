//
//
//
gnuwins(1);

t = [0:2:1/64] ';

N = 10;

_plot_data = <<>>;

for (i in 1:N)
{
  _l = text(i,"%02.0f");
  _y = real(exp(2*pi*1i*t*gaussian()) + exp(2*pi*1i*t*uniform()));
  _plot_data.[_l] = [t, _y];
}

// color definitions from
// https://www2.uni-hamburg.de/Wiss/FB/15/Sustainability/schneider/gnuplot/colors.htm

gnuwin(1);
gnuxlabel ("x");
gnuylabel ("y");
gnuformat([...
  "with lines lt 1 lw 1 lc rgb 0xf0f8ff", ... 
  "with lines lt 1 lw 1 lc rgb 'blue'", ...
  "with lines lt 1 lw 1 lc rgb 'green'", ...
  "with lines lt 1 lw 1 lc rgb 'black'", ...
  "with lines lt 1 lw 1 lc rgb 'orange'", ...
  "with lines lt 1 lw 1 lc rgb 'pink'", ...
  "with lines lt 1 lw 1 lc rgb 'magenta'", ...
  "with lines lt 1 lw 1 lc rgb 0x556B2F", ...
  "with lines lt 1 lw 1 lc rgb 0x006400", ...
  "with lines lt 1 lw 1 lc rgb 'red'", ...
[]]);
gnuplot(_plot_data);

