//
// file: main3.r
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
  y1  = sin(x) + dy1;
  s1  = 0.05*ones(x);
  t = [t; [(x+dx1)/pi, y1]];
}
t = t[sort(t[;1]).idx;];

gnuwins(1, "eps", stderr());
gnucmd ("set output 'gnuplot3.eps';");
gnucmd ("set size 2,2");

// big plot, takes the whole window
gnucmd("set multiplot");
gnucmd("set origin 0,0");
gnuxlabel ("{/=44 x} / {/Symbol=44 p}");
gnuxtics  (0.5,5);
gnuytics  (0.5,5);
gnucmd("unset grid");
gnuylabel ("{/=44 y}");
gnulimits (0,2,-1.5,1.5);
gnulegend (["y=sin(x)","random noise"]);
gnuformat (["with lines","with points"]);
gnuplot   (<<a=[x/pi,y];b=t[;1,2]>>);

// small plot that is inserted on top of the big plot
gnucmd("set origin 1.08,1.08");
gnucmd("set size   0.76,0.76");
gnucmd("set grid");
gnuxtics  (0.5,5);
gnuytics  (0.5,5);
gnucmd("clear");
gnulegend (["y=cos(x)"]);
//gnuxlabel ("x2 / {/Symbol p}");
//gnuylabel ("y2");
gnuplot   (<<a=[x/pi,z]>>);
gnucmd("unset size; unset origin");

gnuclose(1);


