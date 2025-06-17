//
// eg_eig.r: test of eigenvalue solver for dense matrices
//
if (!exist(NITER))
{ NITER = 200; }

x  = rand(100,100);
sx = 0.5*(x+x');

printf("a is symmetric matrix\n");

tic();
for(i in 1:NITER)
{
  //spinner();
  // solve it as a general eigenvalue problem
  y1 = eig(sx, "G");
}
printf("\teig(a,\"G\") took %g sec.\n", toc());

tic();
for(i in 1:NITER)
{
  //spinner();
  // solve it as a general eigenvalue problem
  y2 = eig(sx);
}
printf("\teig(a) took %g sec.\n", toc());

tic();
for(i in 1:NITER)
{
  //spinner();
  // solve it as a general eigenvalue problem
  y3 = eig(sx, "S");
}
printf("\teig(a,\"S\") took %g sec.\n", toc());

