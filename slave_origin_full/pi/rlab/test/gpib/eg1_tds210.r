//
// eg_tds210.r: ch1 terminal on probe compensation
//
rfile libgpib.so

// initialize gpib devices
if (!exist(ud1))
{
  ud1 = ibdev(0,0,0); // init gpib main board, here gpib-usb-hs, pid,sid = 0,0
}
if (!exist(ud3))
{
  ud3 = ibdev(0,12,0); // init tds210, pid,sid = 0,6
}

//
// get the responses of the device
//
reply = [];
for (i in 1:255)
{
  r = ibask(ud3, i);
  if (isempty(r)){ continue; }
  reply = [reply; i, r];
}
printf("Configuration bits for the device tds210 (ud3):\n");
reply

//
// basic commnunication to set-up waveform parameters
//
ibclr( ud3 );   // clear the osciloscope
mycode = [ ...
    "DATA:SOURCE CH1"; ...
    "DATA:ENCDG RIBINARY;WIDTH 1"; ...
    "HORIZONTAL:RECORDLENGTH 1000"; ...
    "DATA:START 1"; ...
    "DATA:STOP 2500"; ...
    "HEADER OFF"; ...
    "ACQUIRE:STATE RUN"; ...
    "CURVE?" ];
ibwrt ( ud3, mycode );

//
//
//
c = ibrd( ud3, 1 );   // read first byte, '#'
c = ibrd( ud3, 1 );   // read second byte, the length of a string containing the length
n = strtod( c );      // convert it to number
L = ibrd( ud3, n );   // read string containing the number of bytes to transfer
N = strtod ( L );     // convert the string to the number
tic();
s = ibrd( ud3, N );   // read in N bytes, each representing a data point
toc()
x = [1:N]';
y = ascii(s)' +0;
plot( [x,y] );






