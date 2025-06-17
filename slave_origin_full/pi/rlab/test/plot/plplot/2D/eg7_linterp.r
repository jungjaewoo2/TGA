//
// linear interpolation
//
if (isfile("plplot.r"))
{
  rfile plplot
}
if (length(plwins().available)<1)
{ plwins (1); }

opts = <<>>;
opts.order = 0;

if (!exist(NITER))
{ NITER = 1e3; }


m = [ ...
  0, 0 ;...
  1, 3 ;...
  2, 1 ;...
  3, 4  ...
  ];

x = [-1:4:1/64]';

for (k in 1:NITER)
{
  y1a = linterp(x, m, 0);
  y1b = linterp(m[;1], m, 0);
  y2a = linterp(x, m);
  y2b = linterp(m[;1], m);

}

//
// plot this
//
plwin(1);
plimits (-1,4,-3,8);
plxlabel ("X");
plylabel ("Y");
plxtics (1,1);
plytics (1,1);
plegend ([ ...
  "table of [x,y]" ,...
  "piece-wise constant and continuous at end-points: overall" ,...
  "at end-points" ,...
  "piece-wise linear and continuous: overall" ,...
  "at end-points" ...
  ], [0.75,1.25], "rit", [0.02,0.03]); // text size 0.75, line spacing 1.25
plformat ([ ...
  "with points ps 1.5 pt 6 lc rgb red" ,...
  "with lines lc rgb blue" ,...
  "with points ps 1 pt 5 lc rgb blue" ,...
  "with lines lc rgb green" ,...
  "with points ps 1 pt 4 lc rgb green" ...
]);
plplot( <<a=m; b=[x,y1a]; c=[m[;1],y1b]; d=[x,y2a]; e=[m[;1],y2b] >>);
_plprint ("eg7.eps", "psc");


