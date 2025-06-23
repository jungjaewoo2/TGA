NITER = 100;


x = rand(7,7);

x = 0.5*(x + x');

z = eig(x);

printf("Dense solver (LAPACK)\n");
printf("Eigenvalues:\n");
z.val
printf("Eigenvectors:\n");
z.vec

s=<<>>;
s.sigma = 1;
s.which = "LM";

printf("Dense solver (ARPACK)\n");

for(i in 1:NITER)
{
  spinner();
  y = eigs(x,2,s);
}

printf("Two eigenvalues closest to 1:\n");
y.val
printf("Their respective  eigenvectors:\n");
y.vec