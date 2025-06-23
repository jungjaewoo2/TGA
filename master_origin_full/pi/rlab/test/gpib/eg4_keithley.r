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

// initialize keithley:
//  1:  SYST:PRES         'Continuous measurement mode (INIT:CONT ON)'
//  2:  FUNC 'VOLT:DC'    'Select DCV function'
//  3:  VOLT:DC:NPLC 0.1  'Set maximum read rate (FAST)'
mycode = [ "SYST:PRES", "FUNC 'VOLT:DC'", "VOLT:DC:NPLC 0.1"];
ibwrt(ud2,mycode);


//
// 7-day measurement of data points at 1sec interval
//

N = 500000;
tvdata=zeros(N,2);

i = 0;
while (i < N)
{
  tic();
  i++;
  // read the data
  r = ibqrd(ud2,["SENS:DATA:FRES?"]);
  // To above query Keithley replies with
  //   -1.10276960E-06VDC,+6606.482SECS,+66586RDNG#
  // this is to be extracted into a matrix with columns
  // [time, voltage]
  datarange = [1:strindex(r,"VDC")-1];
  timerange = [strindex(r,"VDC")+4:strindex(r,"SECS")-1];
  //
  ctime = strtod( substr(r, timerange) );
  if (i==1)
  { zerotime = ctime; }
  
  tvdata[i;1] = ctime - zerotime;
  tvdata[i;2] = 1e6 * strtod( substr(r, datarange) );
  
  sleep(1-toc());
}

writem("test.txt", tvdata);
close ("test.txt");
quit
