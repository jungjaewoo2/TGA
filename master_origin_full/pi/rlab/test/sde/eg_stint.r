//
// file: eg_stint.r
//
// MATLAB version by D.J. Higham, see README.
// rlab version by M. Kostrun, I-2007.
//
// approximate stochastic integrals:
// Ito and Stratonovich

rng(1, "normal", [0, 1]); // gaussian distribution with 0 mean and 1 std
N  = 1024;
t0 = 0;
tf = 1;
dt = (tf - t0)/N;
ti = [t0:tf:dt]';

//
// calculate brownian path
//
dw = sqrt(dt) * rand (N,1);
w  = cumsum(dw);

ito    = sum( [0; w[1:N-1]] .* dw );
itoerr = abs(ito - 0.5*(w[N]^2-tf));
colors ("red");
printf( "ito = ");
colors  ();
printf( "%g +- %g\n", ito, itoerr);

strat    = sum( (0.5*([0; w[1:N-1]] + w) + ...
    0.5 * sqrt(dt) * rand(N,1)) .* dw );
straterr = abs( strat - 0.5 * w[N]^2);
colors ("green");
printf( "stratanovich = ");
colors  ();
printf( "%g +- %g\n", strat, straterr);
