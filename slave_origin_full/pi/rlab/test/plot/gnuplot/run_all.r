//
//
//
_JMAX = int(1000*uniform()) + 1;
_fns = ls("./eg*_*.r");
if (isempty(_fns))
{ stop (); }

_zz = zeros(_JMAX,2);

for (_jj in 1:_JMAX)
{
  _m1 = strtod(getmemoryusage());
  _fn = shuffle(_fns,1);
  printf("%g/%g: executing script %s\n", _jj, _JMAX, _fn);
  NITER = int(10*uniform()) + 5;
  load(_fn);
  _m2 = strtod(getmemoryusage());
  _zz[_jj;] = [_jj, _m2];
  printf("%g/%g: Memory consumption following script %s is %g%% (%g%%).\n\n", _jj, _JMAX, _fn, _m2,  _m2-_m1);
}
colors("red");
printf("Total Memory consumption increase is %g%%.\n\n", _zz[_JMAX;2]-_zz[1;2]);
colors();

// plot memory consumption
gnuwins (1);
gnuxlabel ("Script Iteration");
gnuylabel ("Memory Consumption (%%)");
gnuplot (_zz);

