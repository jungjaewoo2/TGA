//
//
//
//JMAX = int(1000*uniform()) + 1;
NOPLOTS=1;
JMAX = 200;

fns = ls("./eg*pde*.r");
if (isempty(fns))
{ stop (); }

z = zeros(JMAX,2);

for (_jj in 1:JMAX)
{
  m1 = strtod(getmemoryusage());
  fn = shuffle(fns,1);
  printf("%g/%g: executing script %s\n", _jj, JMAX, fn);
  NITER = int(10*uniform()) + 5;
  load(fn);
  m2 = strtod(getmemoryusage());
  z[_jj;] = [_jj, m2];
  printf("%g/%g: Memory consumption following script %s is %g%% (%g%%).\n\n", _jj, JMAX, fn, m2,  m2-m1);
}
colors("red");
printf("Total Memory consumption increase is %g%%.\n\n", z[JMAX;2]-z[1;2]);
colors();

// plot memory consumption
gnuwins (1);
gnuxlabel ("Script Iteration");
gnuylabel ("Memory Consumption (%%)");
gnuplot (z);

