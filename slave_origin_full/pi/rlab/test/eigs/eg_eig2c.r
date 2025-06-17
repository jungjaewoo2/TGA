//
// eg_eig.r: test of eigenvalue solver for dense matrices
//

//
// real
//
if (!exist(NITER))
{ NITER = 100; }

x  = rand(100,100) + 1i*rand(100,100);
sx = 0.5*(x+x');

b  = diag(ones(1,x.nr));

printf("a is symmetric matrix, b is positive definite\n");

tic();
for(i in 1:NITER)
{
  spinner();
  // solve it as a general eigenvalue problem
  y1 = eig(sx,b);
}
printf("\teig(a,b) took %g sec.\n", toc());

tic();
for(i in 1:NITER)
{
  spinner();
  // solve it as a general eigenvalue problem
  y2 = eig(sx,b,"sp");
}
printf("\teig(a,b,\"sp\") took %g sec.\n", toc());

y0=eig(sx).val;
y1=sort(real(y1.val)).val;
y2=sort(real(y2.val)).val;
