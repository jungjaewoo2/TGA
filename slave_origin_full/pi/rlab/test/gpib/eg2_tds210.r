//
// eg_tds210.r: tektronix tds210:ch1 oscilloscope to tm5003:pg506 calibration generator
//
rfile libgpib.so
rfile func_all

gnuwins (1);
plwins (1);

NITER = 100000;
// NITER = 2;

// use least squares fit for edge of the pulse
odropt=<<>>;
odropt.imethod = 2;


// initialize gpib devices
if (!exist(ud1))
{
  ud1 = ibdev(0,0,0); // init gpib main board, here gpib-usb-hs, pid,sid = 0,0
}
if (!exist(ud3))
{
  ud3 = ibdev(0,6,0); // init tds210, pid,sid = 0,6
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

// oscilloscope screen:
TIMEAX  = 10;   // horizontal division
VOLTAX  = 10;   // vertical division (only 8 visible)
// oscilloscope waveform
RECLEN  = 2500; // no. of points
WIDTH   = 1;    // resolution (255 points per 10 vertical divisions)

//
// basic commnunication to set-up waveform parameters
//
ibclr( ud3 );   // clear the osciloscope
mycode = [ ...
    "DATA:SOURCE CH1"; ...
    // single byte data expected
    "DATA:ENCDG RPBINARY"; ...
    "WIDTH " + text(WIDTH,"%g"); ...
    "HORIZONTAL:RECORDLENGTH " + text(RECLEN,"%g"); ...
    "DATA:START 1"; ...
    "DATA:STOP " + text(RECLEN,"%g"); ...
    "HEADER OFF"; ...
    // trigger oscilloscope on rising EDGE when voltage. use the
    // trigger level set by user on the oscilloscope
    "TRIG:MAI:TYP EDGE" ...
    ];
ibwrt ( ud3, mycode );

// obtain information about the axes: read it twice because
// for some, to me unknown reason, the first time an empty string
// might be returned
timescale_ns = strtod(lastr(ibqrd(ud3, , "HOR:MAI:SCA?", 2))) / 1e-9; // (ns)

// set vertical scales for conversion of integers to voltages
voltscale_div  = strtod(lastr(ibqrd(ud3, , "CH1:VOL?", 2))) / 1e-3;     // (mV)
voltscale_res  = 2^(8*WIDTH)-1;

// measurement query:
//  set the trigger on, and return the waveform once oscilloscope is triggered
myquery =  ["ACQUIRE:STATE ON", "CURVE?" ];

// find horizontal axis: time in nanoseconds
x  = ([1:RECLEN]'- RECLEN/2) * TIMEAX * timescale_ns / RECLEN;

s0  = zeros(RECLEN, 1);
sd2 = zeros(RECLEN, 1);
for (i in 1:NITER)
{
  //
  N1 = 0;
  s  = [];
  while (N1 < RECLEN * WIDTH)
  {
    spinner();
    ibclr (ud3);
    ibwrt (ud3, myquery);
    c = ibrd( ud3, 1 );       // read the first byte, should be '#'
    c = ibrd( ud3, 1 );       // read the second byte, this is ascii text of how much bytes follows
    n = strtod(char(c));      // convert it to character, and then to number: this is the number of bytes that follow
    L = ibrd( ud3, n );       // read next 'n' bytes
    N = strtod(char(L));      // convert the string to the number
    if (WIDTH == 1)
    { s = [s, ibrd( ud3, N )]; }  // read in N bytes, each representing a data point
    N1 = N1 + N;
  }
  resize(s, s.nc, s.nr);      // transpose in place

  // s contains integers in range 0-255, covering 10 vertical divisions
  // of which 8 are visible on the screen
  y = (VOLTAX/voltscale_res * (s - voltscale_res/2)) *voltscale_div;

  //
  s0  = s0  + y;
  sd2 = sd2 + y.^2;

  // plot it while we wait
  plimits (floor(min(x)), ceil(max(x)));
//   gnuxtics  ((ceil(max(x))-floor(min(x)))/10,10);
  xlabel ("Time (ns)");
  ylabel ("Voltage (mV)");
  plot   ( [x,y] );
}

// statistics of the waveforms:
// mean
s0 = 1/NITER * s0;
// std
sd = (sd2/NITER- s0.^2).^0.5;

// fit 'vl' the rising edge formula:
// 1. find vmax
_vmax = lastr(s0);
_vmin = s0[1];
_dt = max(diff(x));
// 2. find range of points during which vl increases from 0.1 to 0.9 of vmax
_ta = max(findroot([x,s0 - _vmin], 0.1 * (_vmax-_vmin)));
_tb = min(findroot([x,s0 - _vmin], 0.9 * (_vmax-_vmin)));
i4a = max(find(abs(x - _ta) < _dt));
i4b = min(find(abs(x - _tb) < _dt));
tf = x [i4a:i4b];
vf = s0[i4a:i4b] - _vmin;
pf = zeros(1,3);
pf[1] = _vmax - _vmin;
pf[2] = _ta;
pf[3] = 20;
// construct an array of :
//    rise time as a function of CS
cf = odrfit (vf, tf, pf, vfunc, odropt);

//
// _txt_pos1 = [0.6,0.25];
// _txt_pos2 = [0.6,0.2];

_txt_pos1 = [0.475,0.25];
_txt_pos2 = [0.475,0.2];


//
gnuwin (1);
gnulimits (floor(min(x)), ceil(max(x)));
gnuxtics  ((ceil(max(x))-floor(min(x)))/10,10);
// gnulimits (,,-50,100);
// gnuytics  (10,10);
//gnutext   ("TDS210 (Osc.) / P6139A / PG506 (Calib. Gen.)", _txt_pos1, 1, "screen");
//gnutext   ("TDS210 (Osc.) + 50{/Symbol W} Term./ PG506 (Calib. Gen.)", _txt_pos1, 1, "screen");
gnutext   ("TDS210 (Osc.) / PG506 (Calib. Gen.)", _txt_pos1, 1, "screen");
gnutext   ("Rise time = "+text(_tb - _ta,"%3.0f")+" ns ({/Symbol D}t), " ...
    +text(cf.coef[3],"%3.0f ns (exp)"), _txt_pos2, 1, "screen");
gnuxlabel ("Time (ns)");
gnuylabel ("Voltage (mV)");
gnulegend (["Waveform: spread", "mean", "Fit to raising exponential formula"]);
gnuformat (["using 1:2:3 with yerrorbars lt 1 lc rgb \"red\" lw 0.5", ...
    "with lines lt 1 lc rgb \"green\" lw 0.5", "with lines lt 3 lc rgb \"blue\" lw 2"] );
gnuplot   ( <<a=[x,s0,sd];b=[x,s0];c=[tf, _vmin + vfunc(tf, cf.coef)]>>, ...
    "fig/tds210_pg506-FRO.eps" );





