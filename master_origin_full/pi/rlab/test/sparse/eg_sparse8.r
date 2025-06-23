//
//
//

NITER = 50;

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
printf("Dense solution found in %g sec\n", toc());

tic();
for(i in 1:NITER)
{
  spinner();
  sx = solve(sa,b);
}
printf("Sparse solution found in %g sec\n", toc());

[x,sx]



