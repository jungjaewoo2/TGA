//
// sparse7.r: test of iterative solver
//

n = 1000;

printf("Comparative performance of real sparse solvers\n");

NITER = 100;

dh = 0.1;

a  = zeros(n,n);

for(i in 1:n)
{ a[i;i]   = 1;}

for(i in 1:n-1)
{
  a[i;i+1] = dh*(1-2*rand());
  a[i+1;i] = dh*(1-2*rand());
}

for(i in 1:n-2)
{
  a[i;i+2] = dh*(1-2*rand());
  a[i+2;i] = dh*(1-2*rand());
}

sa = sparse(a);
b  = rand(n,1);
x  = solve(a,b);

sparams.realsolv("spkit");
sparams.iterator(1);
sparams.precond(1);

tic();
for(i in 1:NITER)
{ sx1 = solve(sa,b); }
printf("Calculation with spkit lasted %g sec.\n", toc());


sparams.realsolv("umfpack");
tic();
for(i in 1:NITER)
{ sx2 = solve(sa,b); }
printf("Calculation with umf lasted %g sec.\n", toc());

sparams.realsolv("superlu");
tic();
for(i in 1:NITER)
{ sx3 = solve(sa,b); }
printf("Calculation slu lasted %g sec.\n", toc());

printf("Differences from LAPACK solution:\n");
[max(abs(x-sx1)),max(abs(x-sx2)),max(abs(x-sx3))]

