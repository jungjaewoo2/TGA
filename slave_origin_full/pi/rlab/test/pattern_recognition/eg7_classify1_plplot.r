//
//
//
//
plwins (1);

use_strtod = 1;

NITER = 1;
N = 500;
C = [0,0; 1,1; 0,1; 1,0];
s = [0.2, 0.2, 0.2,0.2];

t = [];
f = [];
tp = <<>>;
for (i in 1:C.nr)
{
  si = text(i, "%02.0f");
  tp.[si] = C[i;] + s[i]*gaussian(N,2);
  t = [t; tp.[si]];
  f = [f; i*ones(N,1)];
}
tset = <<>>;
tset.data = t;
if (use_strtod)
{
  tset.feature = num2str(f);
}
else
{
  tset.feature = f;
}

x  = uniform(50,2);

// classifier 1:
if (use_strtod)
{
  fx1 = strtod(classify.parzen(x, tset));
}
else
{
  fx1 = classify.parzen(x, tset);
}
for (i in 5:8)
{
  si = text(i, "%02.0f");
  tp.[si] = [];
}
for (i in 1:x.nr)
{
  si = text(4 + fx1[i], "%02.0f");
  tp.[ si ] = [tp.[ si ]; x[i;] ];
}

k = 5;
// classifier 2:
if (use_strtod)
{
  fx2 = strtod(classify.knn(x, tset, k));
}
else
{
  fx2 = classify.knn(x, tset, k);
}
for (i in 9:12)
{
  si = text(i, "%02.0f");
  tp.[si] = [];
}
for (i in 1:x.nr)
{
  si = text(8 + fx2[i], "%02.0f");
  tp.[ si ] = [tp.[ si ]; x[i;] ];
}

plwin  (1);
pltitle  ("Classification Demo");
plimits  (-1,2,-1,2);
plegend (["Training classes: 1", "2", "3", "4", ...
    "Classification: Parzen: 1", "2", "3", "4", ...
    "KNN("+text(k,"%.0f")+"): 1", "2", "3", "4"]);
plformat ([...
    "with points lc rgb 'red' ps 1 pt 1", "with points lc rgb 'blue' ps 1 pt 1", ...
    "with points lc rgb 'green' ps 1 pt 1", "with points lc rgb 'brown' ps 1 pt 1", ...
    "with points lc rgb 'red' ps 4 pt 4", "with points lc rgb 'blue' ps 4 pt 4", ...
    "with points lc rgb 'green' ps 4 pt 4", "with points lc rgb 'brown' ps 4 pt 4", ...
    "with points lc rgb 'red' ps 4 pt 3", "with points lc rgb 'blue' ps 4 pt 3", ...
    "with points lc rgb 'green' ps 4 pt 3", "with points lc rgb 'brown' ps 4 pt 3", ...
[]]);
plplot (tp);


