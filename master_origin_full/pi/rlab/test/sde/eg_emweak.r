//
// file: eg_emweak.r
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
// weak convergence at tf=1 by comparing the analytical solution
// to the numerical one.

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

M = 65536; // number of paths
data = zeros(5,2);

// numerical solution for different E-M ti steps
for (p in 1:5)
{
  // numeric solution
  Dt = 2^(p-10);
  L  = (tf-t0) / Dt;
  data[p;1] = Dt;
  xtmp = x0*ones(M,1);
  for (j in 1:L)
  {
    smiley();
    winc   = sqrt(Dt) * rand(M,1);
    // winc = sqrt(Dt) * sign(rand(M,1)); // use for weak E-M
    xtmp   = (1 + Dt*lambda + mu*winc) .* xtmp ;
  }
  data[p;2] = mean( xtmp ); // store solution at t=tf
}
data[;2] = abs( data[;2] - x0*exp(lambda*(tf-t0)) );

//
// plot it
//
gnutitle ("Weak convergence of E-M method");
gnuxlabel("Time step");
gnuylabel("Absolute Error");
gnulimits(0.001,0.1, 0.01,1);
gnuytics ([0.01,0.02,0.05,0.1,0.2,0.5,1]);
gnuxtics ([0.001,0.002,0.005,0.01, 0.02, 0.05, 0.1]);
gnucmd   ("unset grid; set grid xtics ytics;");
gnuscale ("log", "log");
gnuformat(["with lines"]);
gnuplot  (data);


