//
// file: eg_bpath3.r
//
// MATLAB version by D.J. Higham, see README.
// rlab version by M. Kostrun, I-2007.
//
// evaluating an sde along brownian path:
//  u(w(t)) = exp(t + 1/2*w(t))

gnuwins(1);

rng(1, "normal", [0, 1]); // gaussian distribution with 0 mean and 1 std
N  = 1024;
t0 = 0;
tf = 1;
dt = (tf - t0)/N;
ti = [t0:tf:dt]';

M  = 1024;            // number of concurrent paths being calculated

//
// calculate brownian paths
//
dw = [zeros(1,M); sqrt(dt) * rand (N,M)];
w  = cumsum(dw);

// find value of the function u along the path
u = exp(ti * ones(1, M) + 0.5*w);

// find the mean
mu = mean(u')';

// find analytic prediction for the mean
au = exp(9/8*ti);

//
// plot function u along brownian path
//
gnulimits(t0,tf,0,5);
gnuxtics (1/4, 1);
gnuytics (1, 2);
gnuxlabel("Time t");
gnuylabel("u(w(t)) = exp(t + 1/2*w(t))");
gnucmd   ("unset grid; set grid xtics ytics;");
gnuformat("with lines");
gnulegend(["u(w) along path no. 1", ...
    "u(w) along path no. 2", ...
    "u(w) along path no. 3", ...
    "u(w) along path no. 4", ...
    "u(w) along path no. 5", ...
    "mean u - numerical", ...
    "mean u - analytical"]);
gnuplot  ([ti, u[;1:5], mu, au]);
