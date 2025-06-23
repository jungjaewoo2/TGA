
gnuwins(1);

basename="func";
myfunc = "sin";
_leg = blank(0,0);

for (I in 1:3)
{
  t_I = text(I, "%.0f");
  cmd = basename + t_I + " = function(x) {" + ...
    "return " + t_I  + " + " + ...
    myfunc + "("+t_I+"*x);" + ...
    "};";
  cmd
  eval(cmd);
  _leg = [_leg, basename + t_I];
}

t = [0:2:1/64]';

gnuwin(1);
gnulegend(_leg);
gnuplot( <<a=[t, func1(t .* const.pi)];b=[t, func2(t .* const.pi)];c=[t, func3(t .* const.pi)]>> );
