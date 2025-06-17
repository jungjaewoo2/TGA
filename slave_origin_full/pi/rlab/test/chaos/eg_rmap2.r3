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

range_a = [0:4:1/32];

x0 = urandom();
n0 = 10000;
n1 = 100;

data = zeros(length(range_a),n1+1);

tic();
for (i in 1:length(range_a))
{
  spinner();
  data[i;1] = range_a[i];
  data[i;2:n1+1] = rmap(logistic, range_a[i], x0, n1, n0)';
}
printf("built-in[rmap]: %g iterations took %g sec\n", length(range_a)*(n0+n1), toc());

// flatten the matrix for gnuplot (cannot print 100 columns of data on a single plot)
xy = zeros(n1*length(range_a), 2);
for (i in 1:length(range_a))
{
  for (j in 1:n1)
  {
    xy[i + (j-1)*length(range_a);1] = data[i;1];
    xy[i + (j-1)*length(range_a);2] = data[i;j+1];
  }
}

gnuwins(1);
gnuformat("with points");
gnuxlabel("Parameter a");
gnuylabel("Limit n->inf() of x_n");
gnuplot  (xy);


