//
// file: eg_emstrong.r
//
// MATLAB version by D.J. Higham, see README.
// rlab version by M. Kostrun, I-2007.
//
// Solves linear SDE
//   dx = lambda*x*dt + mu*x*dw,
// using Euler-Maruyama method, with
//   x(0) = x0 = 1,
//   lambda    = 2,
//   mu        = 1.
//
// Note: Brownian path over [0,1] had dt = 2^(-9). We use
// 5 different time steps, 16*dt, 8*dt, .., dt, and examine
// strong convergence at tf=1 by comparing analytical solution
// to numerical one.

gnuwins(1);

rng(1, "normal", [0, 1]); // gaussian distribution with 0 mean and 1 std
N  = 512;
t0 = 0;
tf = 1;
dt = (tf - t0)/N;
ti = [dt:tf:dt]';

// problem parameters
lambda = 2;
mu     = 1;
x0     = 1;

M = 1024; // number of paths
xerr = zeros(M,5);
data = zeros(5,2);

for (s in 1:M)
{
  spinner();
  // brownian path
  dw = sqrt(dt) * rand(N,1);
  w  = cumsum(dw);
  // analytical solution
  xtru = x0 * exp( (lambda-0.5*mu^2)*tf + mu * last(w) );

  // numerical solution for different E-M ti steps
  for (p in 1:5)
  {
    // numeric solution
    R    = 2^(p-1);
    Dt   = R * dt;
    L    = N / R;
    data[p;1] = Dt;
    xtmp = x0;
    for (j in 1:L)
    {
      winc   = sum( dw[ (R*(j-1)+1) : (R*j) ] );
      xtmp   = (1 + Dt*lambda + mu*winc) * xtmp ;
    }
    xerr[s;p] = abs(xtmp - xtru);
  }
}
data[;2] = mean(xerr)';

//
// plot it
//
gnutitle ("Strong convergence of E-M method");
gnuxlabel("Time step");
gnuylabel("Error");
gnulimits(0.001,0.1, 0.1,1);
gnuytics ([0.1,0.2,0.5,1]);
gnuxtics ([0.001,0.002,0.005,0.01, 0.02, 0.05, 0.1]);
gnucmd   ("unset grid; set grid xtics ytics;");
gnuscale ("log", "log");
gnuformat(["with lines"]);
gnuplot  (data);


