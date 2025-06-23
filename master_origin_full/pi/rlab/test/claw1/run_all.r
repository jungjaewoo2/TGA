//
//
//
JMAX = (int(1000*uniform()) + 1);
fns = ls("./eg*.r");
for (_jj in 1:JMAX)
{
  spinner();
  fn = shuffle(fns,1);
  printf("%g/%g: executing script %s\n", _jj, JMAX, fn);
  NITER = int(10*uniform()) + 5;
  load(fn);
  printf("%g/%g: Done\n\n", _jj, JMAX);
}


