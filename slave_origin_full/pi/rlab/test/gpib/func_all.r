// rising edge formula
vfunc = function(t,p)
{
  // p:
  //  p[1] -> barV
  //  p[2] -> t0
  //  p[3] -> tau
  // model
  //  V(t) = barV * (1 - exp(-(t-t0)/tau))
  //
  res = p[1].*(1 - exp(-(t-p[2])./p[3]));
  return res;
};
