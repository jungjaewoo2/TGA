//
// Test the Rlab solve function for sparse matrices, based
// on UMFPACK library.
//

NITER = 1;
n = 200;                             // Size of problems.

//
// LHS and RHS for the linear system
//
a   = rand(n,n); sa = sparse(a);
b   = rand (n,1);
a2  = rand(n,n)+1i*rand(n,n);
sa2 = sparse(a2);
b2  = rand (n,1);
a3  = rand(n,n);
sa3 = sparse(a3);
b3  = rand (n,1) +1i*rand(n,1);
a4  = rand(n,n)+1i*rand(n,n);
sa4 = sparse(a4);
b4  = rand (n,1) +1i*rand(n,1);

for(spsolver in ["superlu","umfpack", "jcg"])
{

  printf("Solver: %s\n", toupper(spsolver) );
  sparams.realsolv( spsolver );

  //
  // Test the "solution" of solve...
  //
  tic();

  for(i in 1:NITER)
  {
    xs = solve(sa,b);
    xd = solve(a,b);
  }
  printf("RR: the error between the results is %g\n", max(abs(xs-xd)));

  for(i in 1:NITER){
    zd = solve(a2,b2);
    zs = solve(sa2,b2);
  }
  printf("CR: the error between the results is %g\n", max(abs(zs-zd)));

  for(i in 1:NITER){
    zd = solve(a3,b3);
    zs = solve(sa3,b3);
  }
  printf("RC: the error between the results is %g\n", max(abs(zs-zd)));


  for(i in 1:NITER)
  {
    zd = solve(a4,b4);
    zs = solve(sa4,b4);
  }
  printf("CC: the error between the results is %g\n", max(abs(zs-zd)));

  printf("Solver %s took %g sec to complete the test.\n", spsolver, toc() );

}
