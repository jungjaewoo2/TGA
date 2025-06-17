NITER = 1000;

xdim = 4;

printf("First LAPACK:\n");
x = rand(xdim,xdim) + 1i*rand(xdim,xdim);
x = 0.5*(x+x');

z = eig(x);

printf("Values:\n");
z.val
printf("Vectors:\n");
z.vec

printf("Now arpack:\n");

for(i in 1:NITER)
{
  spinner();
  y = eigs(x,1,1.0);
}

printf("Values:\n");
y.val
printf("Vectors:\n");
y.vec


