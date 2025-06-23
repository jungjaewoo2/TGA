//
// eg_keithley.r: short the terminals for voltage measurement
//

// initialize serial port for communication with
// keithley integra 2700:
//  19200 baud rate, 8N1 with software control, leave default raw communication,
//  and write all communication from the port as characters (debug)
serial = "/dev/ttyS0";
//open(serial, "8N1", 19200, "xon|xoff",,"char");
open(serial, "8N1", 19200, "xon|xoff");


// initialize keithley:
//  1:  SYST:PRES         'Continuous measurement mode (INIT:CONT ON)'
//  2:  FUNC 'VOLT:DC'    'Select DCV function'
//  3:  VOLT:DC:NPLC 1  'Set medium read rate (MED)'
//  WARNING: setting high read rate will clog the port and the rlab will not be able
//           to read from port again.
//           To regain the control over serial port this is what I have to do
//           on my laptop:
//           1. logout
//           2. login on tty console as root (ctrl+alt+f1)
//           3. init 3
//           4. setserial start (if this returns green, you are fine, red is trouble)
//           5. init 5
//           6. logout, and try to login at graphical (ctrl+alt+f7)
//           7. open konsole as a superuser and check whether minicom starts
//           properly. If so, you are lucky, reset the serial port a few times.
//           Succesful reset will be announced by the instrument by beeping.
//           Go back to your rlab application.
initcode = ["SYST:PRES\r", "FUNC 'VOLT:DC'\r", "VOLT:DC:NPLC 1\r"];
readcode = ["SENS:DATA:FRES?\r"];

// number of readings
N = 100;

// data list for plotting
mydata = <<>>;

// use for storing the numerical values
x = zeros(N,2);

writem(serial, initcode);
sleep(0.5);

# To a query 'readcode' Keithley replies with
#   -1.10276960E-06VDC,+6606.482SECS,+66586RDNG#
# this is to be extracted into a matrix with columns
# [time(sec), voltage (uV)]

//
// read N data points through internal loop
//
r = readm(serial, readcode, N);

if (any(strlen(r)==0))
{ stop("Empty response from the port\n"); }

// determine the format of the output. See note above
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
  r = readm(serial, readcode);
  x[i;1] = strtod( substr(r, timerange) );
  x[i;2] = 1e6 * strtod( substr(r, datarange) );
}
mydata.b = x;

//
// read/write query within rlab loop
//
for (i in 1:N)
{
  spinner();
  writem(serial, readcode);
  r = readm(serial);
  x[i;1] = strtod( substr(r, timerange) );
  x[i;2] = 1e6 * strtod( substr(r, datarange) );
}

mydata.c = x;








// plot
xlabel ( "Instrument Time (sec)" );
ylabel ( "Instrument Reading (\gmV)" );
plegend( ["N queries", "1 query", "write/read"] );
plot   ( mydata );

// reading times difference internal/external
ta = last(mydata.a)[1] - mydata.a[1;1];
tb = last(mydata.b)[1] - mydata.b[1;1];
tc = last(mydata.c)[1] - mydata.c[1;1];
printf("Times for %g measuremets\n", N);
printf("\tN queries : %g sec\n", ta);
printf("\t1 query   : %g sec\n", tb);
printf("\twrite/read: %g sec\n", tc);












