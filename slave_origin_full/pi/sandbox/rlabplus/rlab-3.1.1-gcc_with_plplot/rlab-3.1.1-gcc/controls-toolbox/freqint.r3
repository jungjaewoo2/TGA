//----------------------------------------------------------------------
//
// freqint
//
// syntax: W = freqint(A,B,C,D,Npts)
//         W = freqint(NUM,DEN,Npts)
//
// Auto-ranging algorithm for Bode frequency response.
// Generate more points where graph is changing rapidly.
// Calculate points based on eigenvalues and transmission zeros.
//
//----------------------------------------------------------------------
require roots tf2ss tzero

freqint=function(a,b,c,d,npts)
{
  global (eps)

  na = a.nr;
  ma = a.nc;

  if ((nargs==3)&&(na==1)) {
    // Transfer function form.
    npts = c;
    ep = roots(b);
    tz = roots(a);
  } else {
    // State space form
    if (nargs==3) {
       npts=c;
       </a;b;c;d/> = tf2ss(a,b);
    }
    ep = eig(a).val[:];
    tz = tzero(a,b,c,d).z;
  }

  if (isempty(ep)) { ep=-1000; }

  // Note: this algorithm does not handle zeros greater than 1e5
  ez = [ ep[find(imag(ep)>=0)]; tz[find(abs(tz)<1e5&&imag(tz)>=0)] ];
  // Round first and last frequencies to nearest decade
  integ    = abs(ez) < 1e-10; // Cater for systems with pure integrators
  highfreq = round(log10(max(3*abs(real(ez)+integ)+1.5*imag(ez)))+0.5);
  lowfreq  = round(log10(0.1*min(abs(real(ez+integ))+2*imag(ez)))-0.5);

  // Define a base range of frequencies
  diffzp = length(ep)-length(tz);
  itmp = 10*(sum(abs(imag(tz)) < abs(real(tz)))>0);
  if (itmp==[]) { itmp = 0; }
  w = logspace(lowfreq, highfreq, npts+diffzp+itmp);
  ez = ez[find(imag(ez) > abs(real(ez)))];

  // Oscillatory poles and zeros
  if (!isempty(ez)) {
    f=w;
    npts2=2+8/ceil((diffzp+eps)/10);
    ind=sort(-abs(real(ez))).idx;
    z=[];
    for (i in ind) {
      r1=max([0.8*imag(ez[i])-3*abs(real(ez[i])),10^lowfreq]);
      r2=1.2*imag(ez[i])+4*abs(real(ez[i]));
      z=z[find(z>r2 || z<r1)];
      indr=find(w<=r2 && w>=r1);
      f=f[find(f>r2 || f<r1)];
      z=[z,logspace(log10(r1),log10(r2),sum(w<=r2 && w>=r1)+npts2)];
    }
    w=sort([f,z]).val;
  }

  return w[:];

};


