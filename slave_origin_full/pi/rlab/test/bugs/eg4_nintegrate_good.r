// file: eg_ninitegrate2.r
// double integration: take two


f = function(x)
{
  rval = 1;
  return rval;
};

x = [0:1:1/64]';
y1 = zeros(x);
y2 = 0.5 * x.^2;
for (i in 2:len(x))
{
  // construct a simplex
  s = [   0,    0; ...
          0, x[i]; ...
       x[i], x[i]];
  y1[i] = nintsimplex(f, ,s);
}

gnuwins(1);
gnuxlabel ("x");
gnuylabel ("y");
gnulegend (["numeric", "analytic"]);
gnuformat (["with points", "with lines"]);
gnuplot   ( [x, y1, y2] );




