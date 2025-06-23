//
// eg_keithley.r: short the terminals for voltage measurement
//
rfile libgpib.so

// initialize gpib devices
if (!exist(ud1))
{
  ud1 = ibdev(0,0,0); // gpib main board, here gpib-usb-hs, set to 0
}
if (!exist(ud2))
{
  ud2 = ibdev(0,7,0); // set keithley's gpib-address to 7
}

// number of readings
N = 100;

// data list for plotting
mydata = <<>>;

// use for storing the numerical values
x = zeros(N,2);

// initialize keithley:
//  1:  SYST:PRES         'Continuous measurement mode (INIT:CONT ON)'
//  2:  FUNC 'VOLT:DC'    'Select DCV function'
//  3:  VOLT:DC:NPLC 0.1  'Set maximum read rate (FAST)'
mycode = [ "SYST:PRES", "FUNC 'VOLT:DC'", "VOLT:DC:NPLC 0.1"];
ibwrt(ud2,mycode);

//
// read N data points through internal loop
//
r = ibqrd(ud2,["SENS:DATA:FRES?"],N);


# To above query Keithley replies with
#   -1.10276960E-06VDC,+6606.482SECS,+66586RDNG#
# this is to be extracted into a matrix with columns
# [time, voltage]
datarange = [1:strindex(r[1],"VDC")-1];
timerange = [strindex(r[1],"VDC")+4:strindex(r[1],"SECS")-1];

for (i in 1:r.nr)
{
  x[i;1] = strtod( substr(r[i], timerange) );
  x[i;2] = 1e6 * strtod( substr(r[i], datarange) );
}
mydata.a = x;

//
// single query reading within rlab loop
//
for (i in 1:N)
{
  r = ibqrd(ud2,["SENS:DATA:FRES?"]);
  x[i;1] = strtod( substr(r, timerange) );
  x[i;2] = 1e6 * strtod( substr(r, datarange) );
}
mydata.b = x;

//
// read/write query within rlab loop
//
for (i in 1:N)
{
  ibwrt(ud2, "SENS:DATA:FRES?");
  r = ibrd(ud2);
  x[i;1] = strtod( substr(r, timerange) );
  x[i;2] = 1e6 * strtod( substr(r, datarange) );
}
mydata.c = x;


// plot
xlabel ( "Instrument Time (sec)" );
ylabel ( "Instrument Reading (\gmV)" );
plegend( ["N queries", "1 query", "ibwrt/ibrd"] );
plot   ( mydata );

// reading times difference internal/external
ta = last(mydata.a)[1] - mydata.a[1;1];
tb = last(mydata.b)[1] - mydata.b[1;1];
tc = last(mydata.c)[1] - mydata.c[1;1];
printf("Times for %g measuremets\n", N);
printf("\tN queries : %g sec\n", ta);
printf("\t1 query   : %g sec\n", tb);
printf("\tibwrt/ibrd: %g sec\n", tc);












