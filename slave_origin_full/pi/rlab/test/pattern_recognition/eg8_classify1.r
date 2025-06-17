//
//
//
//
gnuwins (1);

NITER = 100;
N0  = 50;
C0  = [0,0; 1,1; 0,1; 1,0];
s0  = 0.1 * ones(1,4);
f0  = [  0,   1,   0,   0]; // OR-operation

t = [];
f = [];
tp = <<>>;
for (i in 1:C0.nr)
{
  si = text(f0[i], "%02.0f");

  if (!exist(tp.[si]))
  { tp.[si] = []; }

  d = C0[i;] + s0[i]*gaussian(N0,2);

  tp.[si] = [tp.[si]; d];

  t = [t; d];
  f = [f; f0[i]*ones(N0,1)];
}
tset = <<>>;
tset.data = t;
tset.feature = text(f);

for (i in 1:NITER)
{
  spinner();
  x = uniform(5000,2);
//   c = classify.fisherq(x,tset,["0","1"]);
  c = classify.fisher(x,tset,["0","1"]);
}
tp.["02"] = x[ find(c=="0"); ];
tp.["03"] = x[ find(c=="1"); ];

gnuwin  (2);
gnutitle  ("Classification: Fisher linear/quadratic discriminant");
gnulegend (["Training classes: 0", "1", ...
    "Classification: Fisher: 0", "1", ...
[]]);
gnuformat ([...
    "with points lc rgb 'red' ps 1 pt 1", "with points lc rgb 'blue' ps 1 pt 1", ...
    "with points lc rgb 'red' ps 4 pt 3", "with points lc rgb 'blue' ps 4 pt 3", ...
[]]);
gnuplot (tp);

