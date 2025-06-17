//
// file: eg_bpath2.r
//
// MATLAB version by D.J. Higham, see README.
// rlab version by M. Kostrun, I-2007.
//
// Brownian path simulation

gnuwins(1);

rng(1, "normal", [0, 1]); // gaussian distribution with 0 mean and 1 std
N  = 1024;
t0 = 0;
tf = 1;
dt = (tf - t0)/N;
ti= [t0:tf:dt]';

//
// calculate brownian path
//
dw = [0; sqrt(dt) * rand (N,1)];
w  = cumsum(dw);

//
// plot brownian path
//
gnulimits(t0,tf,-4,4);
gnuxtics (1/4, 1);
gnuytics (1, 2);
gnuxlabel("Time t");
gnuylabel("Random Walk w(t)");
gnucmd   ("unset grid; set grid xtics ytics;");
gnuformat("with lines");
gnuplot  ([ti, w]);
