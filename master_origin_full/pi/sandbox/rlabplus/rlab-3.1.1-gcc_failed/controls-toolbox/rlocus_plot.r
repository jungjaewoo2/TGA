// ---------------------------------------------------------------------
// rlocus, rlocus_plot:  Evans root locus.
//
// Syntax: R = rlocus(num, den)          transfer function
//         R = rlocus(num, den, k)       transfer function
//         R = rlocus(A, B, C, D)        state space
//         R = rlocus(A, B, C, D, k)     state space
//
//         rlocus_plot(R, title)         plotting
//
// Description:
//
// rlocus(num,den) calculates the locus of the roots of:
//                              num(s)
//                      1 + k ----------  =  0
//                              den(s)
// for a set of gains k which are adaptively calculated to produce a 
// smooth plot using rlocus_plot. Alternatively the vector k can be 
// specified with an optional right-hand argument R=rlocus(num,den,k). 
// Vectors num and den must contain the numerator and denominator 
// coeffiecients in descending powers of s or z.  This function returns
// a list R, 
//    R.rout  - length(k) rows and (length(den)-1) columns containing 
//              the complex root locations.  Each row of the matrix 
//              corresponds to a gain from vector k.  
//    R.k     - the gains k are returned.
//    R.poles - a vector containing the poles
//    R.zeros - a vector containing the zeros
// 
// rlocus(A, B, C, D) finds the root-locus from the equivalent SISO 
// state space system (A,B,C,D):       
//                  .
//                  x = Ax + Bu    u = -k*y
//                  y = Cx + Du
//
// rlocus_plot(R, title) plots the root locus calculated by rlocus.
// The argument 'title' is an optional string for graph title.
//
// Example:
//
//   num = [1,3];
//   den = [2, 1, 1, 2, 0];
//   R = rlocus(num,den);
//   rlocus_plot(R);
//      - OR -
//   rlocus_plot( rlocus(num,den) );
//
// ---------------------------------------------------------------------

rlocus_plot = function ( R, title )
{
  // root locus plot from R
  global (eps,pi,_rlab_config,_rloc_window)
  
  cros = <<[-1,1;1,-1]; [1,1;-1,-1]>>;
  circ = [cos(0:2*pi:0.4)[:], sin(0:2*pi:0.4)[:]];
  mep = max([eps;abs([real(R.poles);real(R.zeros);imag(R.poles);imag(R.zeros)])]);
  if (R.zeros.n==R.poles.n)
  {
     // Round up axis to units to the nearest 5
     ax = 1.2*mep;        
  else
     // Round graph axis    
     exponent = 5*10^(floor(log10(mep))-1);
     ax = 2*round(mep/exponent)*exponent;    
  }
  cros.[1] = cros.[1]*ax/50;
  cros.[2] = cros.[2]*ax/50;
  circ     = circ*ax/50;
  
  xy = <<>>;
  j = 1;
  n = R.rout.nc;
  for (i in 1:n)
  {
     xy.[j] = [real(R.rout[;i]), imag(R.rout[;i])];
     j = j + 1;
  }
  for (i in 1:R.poles.n)
  {
      xy.[j]   = cros.[1] + [real(R.poles[i]),imag(R.poles[i]);real(R.poles[i]),imag(R.poles[i])];
      xy.[j+1] = cros.[2] + [real(R.poles[i]),imag(R.poles[i]);real(R.poles[i]),imag(R.poles[i])];
      j = j + 2;
  }
  for (i in 1:R.zeros.n)
  {
      xy.[j] = [circ[;1] + real(R.zeros[i]), circ[;2] + imag(R.zeros[i])];
      j = j + 1;
  }

  if (_rlab_config.plot_support=="pgplot"||_rlab_config.plot_support=="plplot")
  {
    if (exist(_rloc_window))
    {
       plwin(_rloc_window);
    else
       _rloc_window = plstart();
    }
  }
  if (_rlab_config.plot_support == "pgplot") 
  {
     plcolor(ones(1,14)*7); // optional, here's for pgplot
     plline(ones(1,20));    // optional, here's for pgplot
  }
  plegend();
  plimits(-ax,ax,-ax,ax);
  plwid(4);
  if (exist(title))
  {
    pltitle(title);
  else
    pltitle("Root Locus");
  }
  plstyle("line");
  ylabel("Imag");
  xlabel("Real");
  plot(xy);
  if (_rlab_config.plot_support == "pgplot") 
  {
     plcolor(); 
     plline(); 
  }
  plimits();
  pltitle();
};


