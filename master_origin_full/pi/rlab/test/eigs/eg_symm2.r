NITER = 1000;


x = rand(7,7);

x = 0.5*(x + x');

z = eig(x);

printf("LAPACK:\n");
printf("Values:\n");
z.val
printf("Vectors:\n");
z.vec


sx1 = sparse(x);
sx2 = sx1;

s=<<>>;
s.sigma = 1;
s.which = "LM";

printf("arpack:\n");

for(i in 1:NITER)
{
  spinner();
  y = eigs(sx1,2, 2);
  if( min(min(sx1==sx2))!=1)
  {
    printf("memory leak!\n");
    break;
  }
}

printf("Values:\n");
y.val
printf("Vectors:\n");
y.vec