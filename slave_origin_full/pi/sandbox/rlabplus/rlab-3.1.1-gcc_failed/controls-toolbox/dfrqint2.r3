//-----------------------------------------------------------------------
//
// dfrqint2
//
// Syntax: w = dfrqint2(a,b,c,d,Ts,npts)
//         w = dfrqint2(num,den,Ts,npts)
//
// Discrete auto-ranging algorithm for Nyquist and Nichols plots.
//
//-----------------------------------------------------------------------
require tf2ss d2c freqint2

dfrqint2=function(a,b,c,d,Ts,npts)
{
  global(pi,eps)
  if (nargs==4) {
    Ts = c; npts = d;
    </a;b;c;d/> = tf2ss(a,b);
  }
  </a1;b1/> = d2c(a,b,Ts);
  w=freqint2(a1,b1,c,d,npts);
  w=w(find(w<=pi/Ts));
  if (!isempty(w)) { 
    w=sort([w,linspace(min(w),pi/Ts,128)]);
  } else {
    w=linspace(pi/Ts/10,pi/Ts,128);
  }
  return w;
};

