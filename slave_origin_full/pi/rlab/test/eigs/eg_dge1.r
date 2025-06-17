NITER = 1;


x = rand(7,7);

z = eig(x);

z.val
z.vec

printf("Now arpack:\n");

for(i in 1:NITER)
{
  spinner();
  y = eigs(x,3,0);
}

printf("Values:\n");
y.val
printf("Vectors:\n");
y.vec

