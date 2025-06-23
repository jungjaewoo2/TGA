//
//
//

NITER = 100;

n = 1000;

rfile bandred

a  = rand(n,n);
a  = bandred(a,20);
b  = rand(n,1);
sa = sparse( a );

x  = solve(a,b);

for(rs in ["umf","slu"])
{
  sparams.realsolv( rs );

  tic();
  for(i in 1:NITER)
  {
    spinner();
    x = solve(sa,b);
  }
  printf("%s: solve  lasted  %g sec.\n", rs, toc());

  tic();
  spsolve(sa);
  for(i in 1:NITER)
  { spinner(); sx = spsolve(b); }
  printf("%s: spsolve  lasted  %g sec.\n", rs, toc());
  printf("%s: maxErr = %g\n", rs, max(abs(x-sx)));
}


