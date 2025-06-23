//
// file: eg_intpow.r
//

rng(1,"uniform",[1,20]);

n = int(rand());
k = int(rand());
n = 2;
k = 2;

len2 = int( log(n) / log(2) );

b = 1;
a = k+0;
j = 0;
for (i in 0:len2-1)
{
  if(int(2^i) && int(n))
  {
    b = b * a;
    j++;
  }
  a = a * a;
  j++;
}
b = b * a;
j++;


printf(" %i ^ %i = %g in %g multiplications!\n", k, n, b, j);
printf(" diff = %g\n", b - (k+0)^n);

