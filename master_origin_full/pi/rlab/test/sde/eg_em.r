//
// file: eg_em.r
//
// MATLAB version by D.J. Higham, see README.
// rlab version by M. Kostrun, I-2007.
//
// Euler-Maruyama method for linear SDE,
//   dx = lambda*x*dt + mu*x*dw,
// with
//   x(0) = x0 = 1,
//   lambda    = 2,
//   mu        = 1.

gnuwins(1);

rng(1, "normal", [0, 1]); // gaussian distribution with 0 mean and 1 std
N  = 1024;
t0 = 0;
tf = 1;
dt = (tf - t0)/N;
ti = [dt:tf:dt]';

sol=<<>>;

// problem parameters
lambda = 2;
mu     = 1;
x0     = 1;

//
// calculate brownian path
//
dw = sqrt(dt) * rand(N,1);
w  = cumsum(dw);

// analytic solution
ax = x0 * exp( (lambda-0.5*mu^2)*ti + mu * w );
sol.[1] = [[0; ti], [x0;ax]];

// numeric solution
R    = 4;
Dt   = R * dt;
L    = N / R;
xem  = zeros(L,1);
xtmp = x0;
for (j in 1:L)
{
  winc   = sum( dw[ (R*(j-1)+1) : (R*j) ] );
  xtmp   = (1 + Dt*lambda + mu*winc) * xtmp ;
  xem[j] = xtmp;
}
sol.[2] = [ [0:tf:Dt]', [x0; xem] ];

//
// plot brownian path
//
gnutitle ("Numerical vs. analytical integration of SDE");
gnulimits(t0,tf);
gnuxtics (1/4, 1);
gnuytics (1, 2);
gnuxlabel("Time t");
gnuylabel("Solution");
gnucmd   ("unset grid; set grid xtics ytics;");
gnuformat(["with lines", "with lines"]);
gnulegend(["analytical", "numerical"]);
gnuplot  (sol);
