NITER = 1000;

dima = 10;


a = rand(dima,dima);
a = 0.5*(a+a');
b = diag(ones(1:dima));

z = eig(a,b);

printf("LAPACK:\n");
printf("Values:\n");
z.val
printf("Vectors:\n");
z.vec



s=<<>>;
s.sigma = 1;
s.which = "LM";

printf("arpack:\n");

for(i in 1:NITER)
{
  spinner();
  y = eigs(a, b, 2, 2);
}

printf("Values:\n");
y.val
printf("Vectors:\n");
y.vec


