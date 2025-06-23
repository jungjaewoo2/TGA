// file: eg_ninitegrate1.r
// double integration: take one


printf("This produces a Segmentation Fault:\n");
printf("  rlab is written in c, and the built-in solvers CANNOT be\n");
printf("  called recursively. In this particular example - numerical\n");
printf("  integration - a double integral is represented as an integration\n");
printf("  that calls itself. The code fails but does not crash because I put signaling in.\n");
printf("  Don't do it.\n");

f1 = function(x)
{
  rval = 1;
  return rval;
};
f2 = function(x)
{
  global(_iopts, f1);
  x[1]
  rval = nintegrate(f1, ,[0,x[1]],_iopts);
  //rval = x;
  return rval;
};

# stop()

x = [0:1:1/64]';

// integration options
_iopts = <<>>;
_iopts.ikey = 5;
_iopts.eabs = 1e-9;
_iopts.erel = 1e-12;
_iopts.maxi = 2000;

y = zeros(x);
for (i in 2:len(x))
{
  y[i] = nintegrate(f2, ,[0,x[i]],_iopts);
}




