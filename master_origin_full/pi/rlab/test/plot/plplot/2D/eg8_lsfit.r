//
// eg_fit7.r: least squares fit to model   y(x) = A*exp(-lambda*x) + b
//
if (isfile("plplot.r"))
{
  rfile plplot
}
if (length(plwins().available)<1)
{ plwins (1); }

n = 20;

rng(1,"normal",[0,0.01]);

ffun = function(x,p)
{
  // p:
  // p[1] -> a
  // p[2] -> lambda
  // p[3] -> b
  rval = p[1]*exp(-p[2]*x)+ p[3];
  return rval;
};

dfundp = function(x,p)
{
  // p:
  // p[1] -> A
  // p[2] -> lambda
  // p[3] -> b
  rval = zeros(1,3);
  rval[;1] = exp(-p[2]*x);
  rval[;2] = -x*p[1]*exp(-p[2]*x);
  rval[;3] = 1;
  [p, rval];
  return rval;
};


x = [1:n]';
p0 = [5, 0.1, 1];
y1 = ffun(x,p0) + rand(n,1);

options=<<>>;
options.stdout = rconsole();

p = [6, 0.2, 1];
y = lsfit(y1, x, p, ffun, dfundp, options);

y2 = ffun(x,y.coef);

plwin (1);
plimits  (0,20,0,6);
plxlabel ("x");
plylabel ("y");
plimits(1,n,,);
plformat (["with points pt 2 ps 2 pc rgb red using 1:2", "with lines lt 1 lc rgb black using 1:3"]);
plegend(["raw data", "least-squares fit"], ...
    0.75, "rit", [0.05,0.05]);
plplot( [x,y1[;1],y2] );
_plprint ("eg8.eps", "psc");

printf("Observe: lsfit incorrectly calculates the covariance matrix!\n");
[y.coef; sqrt(diag(y.cov)')/n]


