//
// file: main2.r
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

//gnudefault("debug", 1);

if (1)
{
  W = gnustart(, "set term postscript eps enh color dashed colortext lw 2 \"Times-Roman\" 22 size 6.5,8");
  gnucmd("set output 'gnuplot2_6.5,8.eps'");
else
  W = gnustart(, "set term postscript eps enh color dashed colortext lw 2 \"Times-Roman\" 22 size 6.5,4");
  gnucmd("set output 'gnuplot2_6.5,4.eps'");
}
gnucmd("set multiplot");

// plot no.1 in the same device
gnucmd("set origin 0,0.55");
gnucmd("set size   1,0.45");
gnuxlabel ("x / {/Symbol p}");
gnuylabel ("y");
gnulimits (0,2,-1.5,1.5);
gnuxtics  (1/4,5);
gnuytics  (1/4,5);
gnulegend ("y=sin(x)");
gnuformat ("with lines");
gnuplot   (<<a=[x/pi,y]>>);

// plot no. 2 in the same device
gnucmd("set origin 0,0");
gnucmd("set size   1,0.45");
gnuxtics  (1/2,5);
gnuytics  (1/2,5);
gnuxlabel ("x / {/Symbol p}");
gnuylabel ("y");
gnulegend ("random noise");
gnuformat ("with points");
gnuplot   (<<b=t[;1,2]>>);

gnuclose(W); // this is necessary to flush the buffers



