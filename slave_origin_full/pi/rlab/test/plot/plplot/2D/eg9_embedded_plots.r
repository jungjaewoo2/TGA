//
// file: main2.r
//
if (isfile("plplot.r"))
{
  rfile plplot
}
if (length(plwins().available)<1)
{ plwins (1); }
plwins(1);

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
  y2  = cos(x+dx1) + dy1;
  s1  = 0.05*ones(x);
  t = [t; [(x+dx1)/pi, y1]];
}
t = t[sort(t[;1]).idx;];

s = ifelse(uniform()>=0.5,1,-1);

//
// plot no. 1 in the same device
//
plwin(1); // this one resets the plot
plxlabel ("x / pi", 1);
plylabel ("y", 1);
plimits (0,2,-1.5,2.0);
plxtics  (1/4,5,1);
plytics  (1/4,5);
plegend ("y=sin(x)", 0.75, "lib",[0.05,0.05]);
plformat ("with lines lt 1 lc rgb red");
// plmultiplot([0.15,0.9,0.1,0.9]);
plplot([x/pi,s*y]);

//
// plot no. 2 in the same device
//
// erase area on the large plot where the smaller plot will go:
if (s==-1)
{
  _plfill([0.7,0.2,1.9,1.9],0,0);  // solid pattern in background color
}
// now plot it
plxtics  (1/2,5,0.75);
plytics  (1/2,5,0.75);
plylabel ("");
plxlabel ("");
plegend ("random noise",0.5, "ti", [0.0,0.075]);
plformat ("with points pc 1 pt 1");
plmultiplot([0.55,0.85,0.55,0.85]); // this is normalized area of the smaller plot
plplot   ([x/pi,y2]);
if (s==-1)
{
  _plprint ("eg9-case2.eps", "psc");
else
  _plprint ("eg9-case1.eps", "psc");
}




