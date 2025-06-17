//
//
//

NITER = 10;

n = 500;

a = rand(n,n);
b = rand(n,1);

sa = sparse( a );


tic();
for(i in 1:NITER)
{
  spinner();
  x  = solve(a,b);
}
printf("%g dense solution of a %gx%g matrix with LAPACK found in %g sec\n", NITER, n, n, toc());

sparams.realsolv("slu");
tic();

spsolve(sa);
for(i in 1:NITER)
{
  spinner();
  sx = spsolve(b);
}
printf("Sparse solution with 'spsolve' and SUPERLU found in %g sec\n", toc());

printf("Difference between the two = %g\n", norm(x-sx,2));

