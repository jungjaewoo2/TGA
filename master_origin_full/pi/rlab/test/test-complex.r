#
# Test complex division accuracy...
#

k = 1000000;
n = 10000;

er1r = 0; er1i = 0;
for (i in 1:n)
{
  z = k * (rand() + rand()*1j);
  b = k * (rand() + rand()*1j);
  c = z / b;
  a = c * b;
  er1r = er1r + abs (real(z) - real(a));
  er1i = er1i + abs (imag(z) - imag(a));
}

printf("%8.3g %8.3gj\n", er1r, er1i);
