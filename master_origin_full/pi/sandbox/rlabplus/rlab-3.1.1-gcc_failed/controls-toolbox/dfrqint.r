//--------------------------------------------------------------------------
//
// dfrqint
//
// Syntax: w = dfrqint(a,b,c,d,Ts,npts)
//         w = dfrqint(num,den,Ts,npts)
//
// Discrete auto-ranging algorithm for dbode plots.
//
//--------------------------------------------------------------------------
require tf2ss d2cm freqint 

dfrqint = function(a,b,c,d,Ts,npts)
{
  global(eps,pi)
  
  if (nargs==4) {
    Ts = c; 
    npts = d;
    </a;b;c;d/> = tf2ss(a,b);
  }
  </a1;b1;c1;d1/> = d2cm(a,b,c,d,Ts,"tustin");
  w = freqint(a1,b1,c,d,npts);
  w = w[find(w<=pi/Ts)];
  if (!isempty(w)) { 
    w = sort([w',linspace(min(w),pi/Ts,128)]).val;
  else
    w = linspace(pi/Ts/10,pi/Ts,128);
  }
  return w[:];
};

