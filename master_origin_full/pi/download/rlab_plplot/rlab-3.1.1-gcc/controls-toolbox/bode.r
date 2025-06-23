//--------------------------------------------------------------------------
//
// bode
//
// Syntax: </mag;phase;w/> = bode(a,b,c,d)
//         </mag;phase;w/> = bode(a,b,c,d,iu)
//         </mag;phase;w/> = bode(a,b,c,d,iu,w)
//         </mag;phase;w/> = bode(num,den)
//         </mag;phase;w/> = bode(num,den,w)
//
//	Bode frequency response for continuous-time linear systems.
//	bode(A,B,C,D,IU) produces Bode plot data from the single input IU 
//      to all the outputs of the continuous state-space system (A,B,C,D).
//	IU is an index into the inputs of the system and specifies which
//	input to use for the Bode response.  The frequency range and
//	number of points are chosen automatically.
//
//	bode(NUM,DEN) produces the Bode plot data for the polynomial transfer
//	function G(s) = NUM(s)/DEN(s) where NUM and DEN contain the 
//	polynomial coefficients in descending powers of s. 
//
//	bode(A,B,C,D,IU,W) or bode(NUM,DEN,W) uses the user-supplied 
//	frequency vector W which must contain the frequencies, in 
//	radians/sec, at which the Bode response is to be evaluated.  See 
//	LOGSPACE to generate logarithmically spaced frequency vectors. 
//	When invoked with left hand argument,
//		R = bode(A,B,C,D,...)
//		R = bode(NUM,DEN,...) 
//	returns the frequency vector R.w and matrices R.mag and R.phase (in 
//	degrees) with as many columns as outputs and length(w) rows.  No
//	plot is drawn on the screen. 
//
//      To plot the Bode diagrams, use command
//              R = bode(NUM,DEN,...)
//              bode_plot(R)
//
//	See also: bode_plot, fbode, logspace, margin, nichols, and nyquist.
//
//--------------------------------------------------------------------------
require abcdchk freqint freqresp unwrap

bode = function(a,b,c,d,iu,w)
{
  global(eps,pi)

  if (nargs < 2 || nargs > 6) {
     error("Wrong number of input arguments");
  }

  // --- Determine which syntax is being used ---
  if (nargs==2) {
    // Transfer function form without frequency vector
    num = a; 
    den = b; 
    w   = freqint(num,den,20);
    ny  = num.nr; 
    nn  = num.nc; 
    nu  = 1;

  else if (nargs==3) {
    // Transfer function form with frequency vector
    num = a; 
    den = b;
    w   = c[:];
    ny  = num.nr; 
    nn  = num.nc; 
    nu  = 1;

  else if (nargs==4) {
    // State space system, w/o iu or frequency vector
    msg = abcdchk(a,b,c,d);
    if (msg != "") { error(msg); }
    w = freqint(a,b,c,d,30);    
    ny = d.nr; 
    nu = d.nc;   
    if (ny*nu > 1) {
       //MIMO system 
       iu    = 0;
       mag   = []; 
       phase = [];
       for (i in 1:nu) {
         tmp   = $self(a,b,c,d,i,w);
         mag   = [mag,tmp.magout]; 
         phase = [phase,tmp.phase];
       }
    else
       //SISO system
       iu = 1; 
       nargs = 5;
    }
    if (!iu) {
       return <<mag=mag;phase=phase;w=w>>;
    }

  else if (nargs==5) {
    // State space system, with iu but w/o freq. vector
    msg = abcdchk(a,b,c,d);
    if (msg != "") { error(msg); }
    w   = freqint(a,b,c,d,30);
    ny  = d.nr; 
    nu  = d.nc;

  else
    msg = abcdchk(a,b,c,d);
    if (msg != "") { error(msg); }
    ny  = d.nr; 
    nu  = d.nc;

  }}}} 

  if (nu*ny==0) { return <<mag=[];phase=[];w=[]>>; }
  
  // Evaluate the frequency response
  if (nargs==2||nargs==3) {
    g = freqresp(num,den,sqrt(-1)*w);
  else
    g = freqresp(a,b,c,d,iu,sqrt(-1)*w);
  }

  mag = abs(g);
  phase = (180/pi)*unwrap(atan2(imag(g),real(g)));

  // Uncomment out the following statement if you don't want the phase to  
  // be unwrapped.  Note that phase unwrapping will not always work; it is
  // only a "guess" as to whether +-360 should be added to the phase 
  // to make it more aestheticly pleasing.  (See unwrap.r)

  # phase = (180/pi)*atan2(imag(g),real(g));

  // Uncomment the following line for decibels, but be warned that the
  // MARGIN function will not work with decibel data.
  # mag = 20*log10(mag);
 
  return <<mag=mag;phase=phase;w=w>>;
};
