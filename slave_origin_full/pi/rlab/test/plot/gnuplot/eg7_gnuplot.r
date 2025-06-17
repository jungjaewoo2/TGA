//
// file: main7.r
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

I = gnustart("./testfile7.gnu");
gnuwin(I);
gnuxlabel ("x / {/Symbol p}");
gnuylabel ("y");
gnulimits (0,2,-1.1,1.1);
gnuytics  (1/4, 5);
gnuxtics  (1/4, 5);
gnucmd    ("unset grid; set grid xtics ytics;");
gnulegend (["y=sin(x)", "random noise added"]);
gnuformat (["with lines", "with points"]);
gnuplot(<<a=[x/pi,y];b=t[;1,2]>>);

gnuclose(I);



