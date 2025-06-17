//
//
//
//
gnuwins (2);

NITER = 1;
N0  = 50;
C0  = [0,0; 1,1; 1,0; 0,1];
s0  = 0.01 * ones(1,4);
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

opts = <<>>;
opts.lrate = 1e-2;
opts.maxi  = 1e3;
opts.map   = mean(t')';
for (i in 1:NITER)
{
  spinner();
  c = sammon(tset,1,opts);
}

c_hist = hist(c, 100);

gnuwin  (1);
gnuformat ([...
    "with points lc rgb 'red' ps 1 pt 1", "with points lc rgb 'blue' ps 1 pt 1", ...
    []]);
gnuxlabel ("x-coordinate");
gnuylabel ("y-coordinate");
gnuplot (tp);

gnuwin (2);
gnuxlabel ("sammon 1-dim coordinate");
gnuylabel ("frequency");
gnuplot(c_hist);

