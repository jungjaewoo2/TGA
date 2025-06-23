//
// eg8_tds644.r: ch1 terminal on probe compensation
//
rfile libgpib.so

// initialize gpib devices
if (!exist(ud1))
{
  ud1 = ibdev(0,12,0); // init tds644, pid,sid = 0,12
}

//
// get the responses of the device
//
reply = [];
for (i in 0:7)
{
  r = ibask(ud1, 2 .^ i);
  if (isempty(r)){ continue; }
  reply = [reply; i, r];
}
printf("Configuration bits for the device tds210 (ud1):\n");
reply

//
// basic commnunication to set-up waveform parameters
//
ibclr( ud1 );   // clear the osciloscope
mycode = [ ...
    "DATA:SOURCE CH1"; ...
    "DATA:ENCDG RIBINARY;WIDTH 1"; ...
    "HORIZONTAL:RECORDLENGTH 1000"; ...
    "DATA:START 1"; ...
    "DATA:STOP 2500"; ...
    "HEADER OFF"; ...
    "ACQUIRE:STATE RUN"; ...
    "CURVE?" ];
ibwrt ( ud1, mycode );

//
//
//
c = ibrd( ud1, 1 );
printf ("reply: %s\n", ascii(c));

stop()
c = ibrd( ud1, 1 );   // read second byte, the length of a string containing the length
n = strtod( c );      // convert it to number
L = ibrd( ud1, n );   // read string containing the number of bytes to transfer
N = strtod ( L );     // convert the string to the number
tic();
s = ibrd( ud1, N );   // read in N bytes, each representing a data point
toc()
x = [1:N]';
y = ascii(s)' +0;
plot( [x,y] );






