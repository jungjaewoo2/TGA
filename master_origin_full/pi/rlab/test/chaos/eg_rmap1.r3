//
// file: recursive map
//

// recursive map:
//  x_{n+1} = f(x_n, a)
logistic = function(x, a)
{
  // logistic map function
  rval = a*x*(1-x);
  return rval;
};

a = 3.99999;
x0 = urandom();
n0 = 900000;
n1 = 100000;

// built-in function
tic();
y1  = rmap(logistic, a, x0, n1, n0);
printf("built-in[rmap]: %g iterations took %g sec\n", n0+n1, toc());

// script
x1 = x0;
y2 = zeros(n1,1);
tic();
for (i in 1:n0)
{ x1 = logistic(x1,a); }
for (i in 1:n1)
{
  x1 = logistic(x1,a);
  y2[i] = x1;
}
printf("scripted func : %g iterations took %g sec\n", n0+n1, toc());

max(y2-y1)


