
// The RLAB benchmarks.
//
// This demo file runs a set of 6 benchmarks:
//
//    1) N=200 Real matrix multiply
//    2) N=200 Real matrix inverse
//    3) N=120 Real eigenvalues
//    4) 65536-point complex FFT
//    5) N=200 solve
//    6) 1000  iteration FOR loop
//
// The benchmarks illuminate computer architectural issues that
// affect the speed of software like RLAB on different machines.
//
// The time for your machine is measured and displayed versus times
// we've already taken from other standard machines.
//
// Thank you for your contributions:
// Ian Searle, Paul Bergman, Karl Storck, Andrew Grzegorek, and
// Chun-Ching Li
//
clearall();
iter = 10;  // evaluate them 10 times
            // change this number to 1 for slow machines

printf("Start our benchmark measurements...\n");
ts = zeros(1,6);

 // [azg] undefined? //srand(13);
rng(1,"uniform",[0,1]);

// rng(1,"uniform",0,1);
a = rand(200,200);

//  200 by 200 real multiply

tic();
for(i in 1:iter)
{
  b = a*a;
}
ts[1]=toc()/iter;

// N=200  real inverse

tic();
for(i in 1:iter)
{
    b = inv(a);
}
ts[2]=toc()/iter;

// N=120 real eigenvalues

a = rand(120,120);
tic();
for(i in 1:iter)
{
    b = eig(a);
}
ts[3]=toc()/iter;

// 65536 point complex FFT

a = rand(1,65536) + 1i;
tic();
for(i in 1:iter)
{
   b = fft(a);
}
ts[4]=toc()/iter;

// N = 200 solve
neq = 200;
a=rand(neq,neq);
b=rand(neq,1);
tic();
for(i in 1:iter)
{
   c = solve(a,b);
}
ts[5] = toc()/iter;


// 1000 FOR loops

a[1000] = 0;
iter2 = iter*100;
tic();
for(j in 1:iter2)
{
     for (i in 1:1000) { a[i] = 1; }
}
ts[6] = toc()/iter2;

// some benchmark data
times = [6.800e-01, 5.510e-01, 1.023e+00, 4.550e-01, 1.900e-01, 9.370e-03;...
         5.740e-01, 7.140e-01, 1.225e+00, 3.230e-01, 2.330e-01, 1.366e-02;...
         1.091e+00, 1.142e+00, 2.263e+00, 1.004e+00, 3.610e-01, 1.983e-02;...
         1.676e+00, 2.344e+00, 3.706e+00, 1.235e+00, 7.190e-01, 4.102e-02;...
         1.663e+00, 2.330e+00, 3.719e+00, 1.334e+00, 7.150e-01, 3.335e-02;...
         6.264e+00, 7.342e+00, 1.465e+01, 2.623e+00, 2.400e+00, 2.029e-01;...
         2.624e+00, 2.567e+00, 4.437e+00, 1.814e+00, 8.940e-01, 2.924e-02;...
         1.072e+00, 1.175e+00, 1.709e+00, 5.920e-01, 4.060e-01, 1.048e-02;...
         4.730e-01, 5.950e-01, 1.012e+00, 2.850e-01, 1.960e-01, 1.499e-02;...
         1.584e+00, 1.757e+00, 3.905e+00, 1.226e+00, 6.100e-01, 2.824e-02;...
         4.298e+00, 5.670e+00, 6.467e+00, 1.668e+00, 1.690e+00, 2.712e-02;...
         3.170e-01, 4.010e-01, 6.880e-01, 3.400e-01, 1.260e-01, 2.503e-02;...
         3.851e+00, 3.974e+00, 5.621e+00, 2.153e+00, 1.396e+00, 6.185e-02;...
         1.406e+02, 1.739e+02, 2.876e+02, 4.792e+01, 5.368e+01, 9.217e-01;...
         1.260e-01, 2.190e-01, 3.860e-01, 2.780e-01, 6.600e-02, 1.309e-02;...
         1.909e+01, 1.895e+01, 3.620e+01, 8.577e+00, 6.986e+00, 1.058e-01;...
         5.280e-01, 9.550e-01, 2.683e+00, 1.175e+00, 3.310e-01, 3.916e-02;...
         2.080e-01, 3.900e-01, 1.108e+00, 8.180e-01, 1.350e-01, 1.901e-02;...
         1.160e-01, 2.320e-01, 8.650e-01, 2.780e-01, 8.900e-02, 1.782e-02;...
         5.430e-01, 9.560e-01, 2.619e+00, 1.026e+00, 3.250e-01, 3.636e-02;...
         2.981e+00, 4.852e+00, 1.070e+01, 2.166e+00, 1.382e+00, 1.802e-01;...
         3.800e-01, 9.560e-01, 1.561e+00, 5.620e-01, 2.710e-01, 1.044e-02;...
         1.360e-01, 4.040e-01, 7.190e-01, 3.240e-01, 1.230e-01, 1.007e-02;...
         1.380e+00, 1.319e+00, 1.640e+00, 5.318e-01, 4.446e-01, 1.552e-02;...
         3.400e-02, 9.300e-02, 2.740e-01, 2.030e-01, 3.000e-02, 1.440e-02;...
         2.100e-02, 6.100e-02, 1.780e-01, 1.480e-01, 1.900e-02, 8.490e-03;...
         2.200e-02, 1.040e-01, 2.290e-01, 1.820e-01, 4.400e-02, 1.380e-02;...
         ts];
desc = ["PentiumPro 200", "Dell XPS-Pro200n 32MB Ram, 256k L2 Cache";...
        "SUN Ultra1 167", "SUN Ultra-1 Creator 3D, 640MB Ram";...
        "SUN SS-20 125",  "SUN SPARCstation 20, 2 CPUs, 256MB Ram";...
        "SUN SS-10 50",   "SUN SPARCserver 10, 128MB Ram";...
        "SUN SS-1000 50", "SUN SPARCserver 1000, 8 CPUs, 512MB Ram";...
        "SUN SS-LX 50",   "SUN SPARCstation LX, 64MB Ram";...
        "Pentium 90",     "NEC Versa 4050C Notebook, 40MB Ram";...
        "AMD K6 200",     "Home assem, ASUS TX97, 32MB Ram, 512K L2 cache";...
        "SUN Ultra2 200", "Ultra-Sparc-II (200 MHz), 128MB Ram";...
        "Pentium 120",    "Toshiba Tecra 500CDT, 256k L2 cache, 48MB Ram";...
        "PMac 7200 120",  "PowerMac 7200 120MHz (all optimizations off)";...
        "IBM RS6000",     "IBM RS6000";...
        "HP-715/50 50",   "HP-715/50 50MHz, 96 MB Ram";...
        "MacIIsi 20",     "Macintosh IIsi, 68030+68882fpu 20MHz, 17MB Ram";...
        "HP C160",        "HP 9000/780, 384MB Ram";...
        "Intel 80486 25", "486 25MHz PC clone 16MB RAM";...
        "SUN SS-1000 50*","SUN SPARCserver 1000, 8 CPUs, 512MB Ram";...
        "SUN SS-20 125*", "SUN SPARCstation 20, 2 CPUs, 256MB Ram";...
        "SUN Ultra1 167*","SUN Ultra-1 Creator 3D, 640MB Ram";...
        "SUN SS-10 50*",  "SUN SPARCserver 10, 128MB Ram";...
        "SUN SS-LX 50*",  "SUN SPARCstation LX, 64MB Ram";...
        "AMD K6 200*",    "Home assem, ASUS TX97, 32MB Ram, 512K L2 cache";...
        "PentiumPro 200*","Dell XPS-Pro200n 32MB Ram, 256k L2 Cache";...
        "PPro 200, WinNT","PentiumPro 200MHz, 64MB, 256K L2 Cache, Window NT";...
        "HP C160*",       "HP 9000/780, 384MB Ram, Convex MLIB library";...
        "HP J2240*",      "HP J2240 (2 processors) with HP optimized libraries (MLIB)";...
        "HP J2240",       "HP J2240 (2 processors) no optimized libraries";...
        "-Your machine-", "?"];

m = times.nr;
n = times.nc;
printf("      ------- Table of benchmark times (in sec) -------\n\n");
printf(" --Machine-MHz--     *       inv      eig      fft      solve    for\n");
for (i in 1:m)
{
  printf(" %-15s ",desc[i;1]);
  for (j in 1:n)
  {
    printf("%8.4f ",times[i;j]);
  }
  printf("\n");
}
printf("\n Note: * - use performance library\n");
pause();
ratio = ones(m,1)*times[7;]./times;

printf("\n      ------- Table of speed ratios to Pentium 90 -------\n\n");
printf(" --Machine-MHz--    *       inv      eig      fft      solve    for\n");
for (i in 1:m)
{
  printf(" %-15s ", desc[i;1]);
  for (j in 1:n)
  {
    printf("%8.4f ", ratio[i;j]);
  }
  printf("\n");
}
printf("\n Note: * - use performance library\n");
pause();
merit = prod(ratio').^(1/n);
ind = sort(-merit).idx;
printf(" To combine these numbers into a single 'figure of merit' for\n");
printf(" each machine, we compute the geometric mean. Here are the results:\n\n");
printf(" Rank --Machine-MHz--   Merit ----Detailed Hardware Info----\n");
for (i in 1:m)
{
  printf(" %4d %-15s %7.3f %s\n", i, desc[ind[i];1], merit[ind[i]], desc[ind[i];2]);
}
printf("\n Note: * - use performance library\n");
pause();

printf("The SOLVE benchmark compares the performance of different computer\n");
printf("systems while solving 200th order dense systems of linear equations\n\n");

printf("For solving a system of 200 equations, approximately\n");
nflops = 2*neq^3/3 + 2*neq^2;
printf("    nflops = 2/3*200^3 + 2*200^2 = %.5g\n",nflops);
printf("operations are performed.  Using our times from above, we find the\n");
printf("MFlop/second throughput for the various machines:\n\n");

MFlops = nflops./times[;5]'/1000000;
ind = sort(-MFlops).idx;
printf(" Rank --Machine-MHz--   MFlops ----Detailed Hardware Info----\n");
for (i in 1:m)
{
  printf(" %4d %-15s %8.2f %s\n", i, desc[ind[i];1], MFlops[ind[i]], desc[ind[i];2]);
}
printf("\n Note: * - use performance library\n");
printf("\nJust for reference, the Cray X-MP achieves 33 MFlops!\n");
printf("Since rlab does some extra calculations, this number is not the raw speed.\n");
printf("For example, the Pentium Pro achieves 60 Mflops running the linpack benchmark.\n");
printf("All machines were running a single threaded binary of rlab using one CPU,\n");
printf("even some machines have more than one CPU.\n\n");
printf("Please email this line and a short description of your machine to\n");
printf("             rlab-list@eskimo.com or yang@isec.com\n");
printf("Thanks.\n");
for (i in 1:ts.n) {
    printf("%.3e ", ts[i]);
}
printf("\n");
