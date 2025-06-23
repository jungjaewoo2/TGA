// file: eg_mm1.r
//
// test of the matrix market read of a 4x3 matrix
 

N = 1000;

//
// real read
//
tic();
printf("%g real reads: ", N);
fname = "data/demo__1.mtx";
x1 = zeros(4,3);
for (i in 1:N)
{
  spinner();
  readmm(fname, x1);
  x2 = readmm(fname);
}
printf(" Done in %g sec\n", toc());

//
// symmetic read
//
tic();
printf("%g symmetric real reads: ", N);
fname = "data/demo__3.mtx";
x1 = zeros(4,4);
for (i in 1:N)
{
  spinner();
  readmm(fname, x1);
  x2 = readmm(fname);
}
printf(" Done in %g sec\n", toc());

//
// skew read
//
tic();
printf("%g skew real reads: ", N);
fname = "data/demo__4.mtx";
x1 = zeros(4,4);
for (i in 1:N)
{
  spinner();
  readmm(fname, x1);
  x2 = readmm(fname);
}
printf(" Done in %g sec\n", toc());

//
// complex read
//
tic();
printf("%g complex reads: ", N);
fname = "data/demo__2.mtx";
z1 = zeros(4,3) + 1i*zeros(4,3);
for (i in 1:N)
{
  spinner();
  readmm(fname, z1);
  z2 = readmm(fname);
}
printf(" Done in %g sec\n", toc());

