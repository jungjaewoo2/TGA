//
// matrix power to an integer exponent
//

NITER = 1000;

rng(1,"uniform",[1,1000]);
n = int(rand());
//n=1000;

rng(2,"normal",[0,0.25]);
x = rand(8,8);
//x = diag([1,1,1,1]);

printf("Calculating x^%g %g times, where dim(x)=%g:\n", n, NITER, x.nr);

// brute force
tic();
for(k in 1:NITER)
{
  spinner();
  y1 = x;
  for (i in 2:n)
  { y1 = y1*x; }
}
printf("Brute force calculation lasted %g sec.\n", toc());

// builtin function
tic();
for(k in 1:NITER)
{
  spinner();
  y2 = mpow(x+0*1i,n);
}
printf("Built-in function took %g sec.\n", toc());

printf("Difference between the two = %g\n", max(max(abs(y2-y1)))/max(max(abs(y1))));



