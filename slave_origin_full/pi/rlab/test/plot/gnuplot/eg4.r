//
// file: main1.r
//

rng (1, "normal", [0,0.05]);

t = [];
x = const.pi * [;0:2:1/32]';
y = sin(x);
z = cos(x);
for (i in 1:10)
{
  dy1 = rand(x);
  dx1 = rand(x);
  y1  = sin(x+dx1) + dy1;
  s1  = 0.05*ones(x);
  t = [t; [(x+dx1)/pi, y1]];
}
t = t[sort(t[;1]).idx;];
//t = rmdbls(t, 1, [2,3]);

// choose two random colors from X11 pallete for plotting
// the data sets
_c  = members(gnudefault().color);
_i1 = 1 + int(length(_c) * urandom());
_i2 = 1 + int(length(_c) * urandom());
_i3 = 1 + int(length(_c) * urandom());

gnuwins(1,,stderr());

gnuwin(1);

gnutext    ("Text written in rgb '"+_c[_i3]+"'", [0.5, 0.5], 1,,...
    "tc @"+_c[_i3]+"@");

gnuxlabel ("x / {/Symbol p}");
gnuylabel ("y");
gnulimits (0,2,-1.1,1.1);
gnuytics  (1/4, 5);
gnuxtics  (1/4, 5);
gnucmd    ("unset grid; set grid xtics ytics;");
gnulegend (["function y=sin(x) in rgb '" + _c[_i1] + "'", ...
    "y with some random noise in rgb '" + _c[_i2] + "'"]);
// gnuformat (["with lines lw 3 lc rgb \"green\"", "with points lw 0.5 lc rgb \"red\""]);
gnuformat (["with lines lw 3 lc @"+_c[_i1]+"@", ...
    "with points lw 0.5 lc @"+_c[_i2]+"@"]);
gnuplot   (<<a=[x/pi,y];b=t[;1,2]>>);


//gnuclose(1);


