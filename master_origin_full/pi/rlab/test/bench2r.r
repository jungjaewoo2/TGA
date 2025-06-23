printf("\nThe following benchmark program will print the average timings\n");
printf("to calculate the functions by 3 runs.\n");
printf("\n\n");
printf("!!! RLABPLUS - Benchmarkprogram !!!\n");
printf("=================================\n\n");
format(32);
rng(1,"uniform",[0,1]);

## Misc. operation ##

ergeb = 0;
zaehl = 3;
for (i in 1:zaehl)
{
  tic();
  b = abs(rand(1000,1000)/10);
  a = b';
  b = reshape(a,500,2000);
  a = b';
  zeit = toc();
  ergeb = ergeb+zeit;
}
ergeb = ergeb/zaehl;
printf("Creation, trans. & reshaping of a 1000x1000 matrix (sec.)___: %f\n",ergeb);
ergeb = 0;
for (i in 1:zaehl)
{
  b = abs(1+rand(1000,1000)/100);
  tic();
  a = b.^1000;
  zeit = toc();
  ergeb = ergeb+zeit;
}
ergeb = ergeb/zaehl;
printf("1000x1000 normal distributed random matrix^1000 (sec.)______: %f\n",ergeb);
ergeb = 0;
for (i in 1:zaehl)
{
  a = rand(2000000,1);
  tic();
  b = sort(a).val;
  zeit = toc();
  ergeb = ergeb+zeit;
}
ergeb = ergeb/zaehl;
printf("2000000 values sorted ascending (sec.)______________________: %f\n",ergeb);

## Analysis ##

ergeb = 0;
for (i in 1:zaehl)
{
  a = rand(1048576,1);
  tic();
  b = fft(a);
  zeit = toc();
  ergeb = ergeb+zeit;
}
ergeb = ergeb/zaehl;
printf("FFT over 1048576 values (sec.)______________________________: %f\n",ergeb);

printf("Double integration (sec.)___________________________________ : n.a.\n");
printf("Triple integration (sec.)___________________________________ : n.a.\n");

## Algebra ##

ergeb = 0;
for (i in 1:zaehl)
{
  a = rand(1000,1000)/5;
  tic();
  b = det(a);
  zeit = toc();
  ergeb = ergeb+zeit;
}
ergeb = ergeb/zaehl;
printf("Determinant of a 1000x1000 random matrix (sec.)_____________: %f\n",ergeb);
ergeb = 0;
for (k in 1:zaehl)
{
	tic();
	x = eye(1400,1400);
	for (i in 1:1400)
	{
	    for (j in 1:1400)
	    {
		x[j;i]=abs(i-j)+1;
	    }
	}
	zeit = toc();
	ergeb = ergeb+zeit;
}
ergeb = ergeb/zaehl;
printf("1400x1400 Toeplitzmatrix (sec.)_____________________________: %f\n",ergeb);
ergeb = 0;
for (i in 1:zaehl)
{
  a = rand(1000,1000);
  tic();
  b = inv(a);
  zeit = toc();
  ergeb = ergeb+zeit;
}
ergeb = ergeb/zaehl;
printf("Inverse of a 1000x1000 uniform distr. random matrix (sec.)__: %f\n",ergeb);
ergeb = 0;
for (i in 1:zaehl)
{
  a = rand(600,600);
  tic();
  c = eig(a).val;
  zeit = toc();
  ergeb = ergeb+zeit;
}
ergeb = ergeb/zaehl;
printf("Eigenval. of a normal distr. 600x600 randommatrix (sec.)____: %f\n",ergeb);
ergeb = 0;
for (i in 1:zaehl)
{
  a = rand(1000,1000);
  a = a'*a;
  tic();
  b = chol(a);
  zeit = toc();
  ergeb = ergeb+zeit;
}
ergeb = ergeb/zaehl;
printf("Cholesky decomposition of a 1000x1000-matrix (sec.)_________: %f\n",ergeb);
ergeb = 0;
for (i in 1:zaehl)
{
  a = rand(1000,1000);
  tic();
  b = a'*a;
  zeit = toc();
  ergeb = ergeb+zeit;
}
ergeb = ergeb/zaehl;
printf("1000x1000 cross-product matrix (sec.)_______________________: %f\n",ergeb);

## Number theory ##

phi = 1.61803398874989;
ergeb = 0;
for (i in 1:zaehl)
{
  a = floor(1000*rand(500000,1));
  tic();
  b = ceil((phi.^a-(-phi).^(-a))/sqrt(5));
  zeit = toc();
  ergeb = ergeb+zeit;
}
ergeb = ergeb/zaehl;
printf("Calculation of 500000 fibonacci numbers (sec.)______________: %f\n",ergeb);

## Stochastic-statistic ##

ergeb = 0;
a = 100*rand(1000,1000);
for (i in 1:zaehl)
{
  tic();
  b = Gamma(a);
  zeit = toc();
  ergeb = ergeb+zeit;
}
ergeb = ergeb/zaehl;
printf("Gamma function over a 1000x1000 matrix (sec.)_______________: %f\n",ergeb);


ergeb = 0;
a = 100*rand(1000,1000);
for (i in 1:zaehl)
{
  tic();
  b = Erf(a);
  zeit = toc();
  ergeb = ergeb+zeit;
}
ergeb = ergeb/zaehl;
printf("Gaussian error function over a 1000x1000 matrix (sec.)______: %f\n",ergeb);
ergeb = 0;
for (i in 1:zaehl)
{
  b = zeros(1000,1)+(1:1000)';
  for (j in 1:4) {b[j;1]=j;}
  a = rand(1000,1000);
  b = 1:1000;
  tic();
  b = a\b';
  zeit = toc();
  ergeb = ergeb+zeit;
}
ergeb = ergeb/zaehl;
printf("Linear regression over a 1000x1000 matrix (sec.)____________: %f\n",ergeb);
