//
// eg_fisher.r: or
//
NITER = 1;

plwins (3);

//
// create and plot training set
//
x0 = [0,0]; // 2 inputs
y0 = [  0]; // 1 output
x1 = [1,0; 0,1; 1,1]; // 2 inputs
y1 = [  1;   1;   1]; // 1 output

x0train = [];
y0train = [];
x1train = [];
y1train = [];
for (i in 1:25)
{
  rng (1, "normal", [0, 0.2]);
  x0train = [x0train; x0 + rand(x0)];
  y0train = [y0train;           y0 ];
  rng (1, "normal", [0, 0.1]);
  x1train = [x1train; x1 + rand(x1)];
  y1train = [y1train;           y1 ];
}

trainset=<<>>;
trainset.data     = [x0train; x1train];
trainset.feature  = [y0train; y1train];

plwin   (1);
plimits (-0.5,1.5,-0.5,1.5);
pltitle ( "Fisher discriminant: Training set for OR-function" );
// plformat (["with points pt 1 pc rgb black ps 3", "with points pt 1 pc rgb red ps 3"]);
plformat (["with points", "with points"]);
plegend ( ["0", "1"] );
plxlabel  ( "x#d1#u");
plylabel  ( "x#d2#u");
plplot    ( <<a=x0train;b=x1train>> );


//
// create testing set
//
rng (1, "uniform", [-0.5, 1.5]);
xs = rand (1000, 2);

for (i in 1:NITER)
{
  f1 = classify.fisher  (xs,trainset);
  f2 = classify.fisherq (xs,trainset);
}

idx = find (f1 == 0);
x0test = xs[idx;];
idx = find (f1 == 1);
x1test = xs[idx;];

plwin   (2);
plimits (-0.5,1.5,-0.5,1.5);
pltitle ( "Linear Fisher discriminant: Testing set for OR-function" );
plegend ( ["0", "1"] );
plxlabel  ( "x#d1#u");
plylabel  ( "x#d2#u");
plformat (["with points pt 1 pc rgb black ps 3", "with points pt 1 pc rgb red ps 3"]);
plplot    ( <<a=x0test;b=x1test>> );

idx = find (f2 == 0);
x0test = xs[idx;];
idx = find (f2 == 1);
x1test = xs[idx;];

plwin   (3);
plimits (-0.5,1.5,-0.5,1.5);
pltitle ( "Quadratic Fisher discriminant: Testing set for OR-function" );
plegend ( ["0", "1"] );
plxlabel  ( "x#d1#u");
plylabel  ( "x#d2#u");
plformat (["with points pt 1 pc rgb red ps 3", "with points pt 1 pc rgb black ps 3"]);
plplot    ( <<a=x0test;b=x1test>> );

