NITER = 10;

dima = 14;


a = rand(dima,dima) + 1i*rand(dima,dima);
a = 0.5*(a+a');
b = diag(ones(1:dima));

z = eig(a,b);

printf("LAPACK:\n");
printf("Values:\n");
real(z.val)
printf("Vectors:\n");
z.vec



s=<<>>;
s.sigma = 0;
s.which = "SM";

printf("arpack:\n");

for(i in 1:NITER)
{
  spinner();
  y = eigs(a, b, 2, s);
}

printf("Values:\n");
y.val
printf("Vectors:\n");
y.vec

