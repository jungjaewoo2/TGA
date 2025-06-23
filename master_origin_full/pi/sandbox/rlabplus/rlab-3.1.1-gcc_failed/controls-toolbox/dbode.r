//--------------------------------------------------------------------------
// dbode
//
// Syntax: </mag;phase;w/> = dbode(a,b,c,d,Ts,iu,w)
//
//	Bode frequency response for discrete-time linear systems.
//
//	dbode(A,B,C,D,Ts,IU) produces a Bode plot from the single input IU
//	to all the outputs of the discrete state-space system (A,B,C,D).
//	IU is an index into the inputs of the system and specifies which
//	input to use for the Bode response.  Ts is the sample period.  
//	The frequency range and number of points are chosen automatically.
//
//	dbode(NUM,DEN,Ts) produces the Bode plot for the polynomial 
//	transfer function G(z) = NUM(z)/DEN(z) where NUM and DEN contain
//	the polynomial coefficients in descending powers of z. 
//
//	dbode(A,B,C,D,Ts,IU,W) or dbode(NUM,DEN,Ts,W) uses the user-
//	supplied freq. vector W which must contain the frequencies, in 
//	radians/sec, at which the Bode response is to be evaluated.  
//	Aliasing will occur at frequencies greater than the Nyquist 
//	frequency (pi/Ts).  See logspace to generate log. spaced frequency
//	vectors.  
//	It returns the frequency vector w and matrices mag and phase (in 
//	degrees) with as many columns as outputs and length(w) rows.
//
//	See also: dbode_plot, logspace, margin, and dnyquist.
//--------------------------------------------------------------------------

require abcdchk dfrqint freqresp

dbode = function(a,b,c,d,Ts,iu,w)
{
  global(eps,pi)
  
  if(nargs < 3 || nargs > 7) { 
    error("Wrong number of input argumaents.");
  }
  // --- Determine which syntax is being used ---
  if (nargs==3) {
    num = a; 
    den = b;
    if (length(c)==1)
    {
      // Transfer function without frequency vector
      Ts = c;
      w  = dfrqint(num,den,Ts,30);
    else  		
      // Assume this is the old syntax DBODE(NUM,DEN,W)
      disp("Warning: You are using the old syntax.  Use DBODE(NUM,DEN,Ts,W) instead.")
      Ts = 1; w=c;
    }
    ny = num.nr;
    nn = num.nc;
    nu = 1;

  else if (nargs==4) { 	
    // Transfer function form with frequency vector
    num = a; 
    den = b;
    Ts  = c; 
    w   = d;
    ny  = num.nr;
    nn  = num.nc;
    nu  = 1;

  else if (nargs==5) {
    // State space system without iu or freq. vector
    msg = abcdchk(a,b,c,d);
    if (msg != "") { error(msg); }
    w = dfrqint(a,b,c,d,Ts,30);
    ny = d.nr;
    nu = d.nc;
    if (ny*nu > 1) {
       //MIMO system 
       iu    = 0;
       mag   = []; 
       phase = [];
       for (i in 1:nu) {
         tmp   = $self(a,b,c,d,Ts,i,w);
         mag   = [mag,tmp.mag]; 
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

  else if (nargs==6) {
    error(abcdchk(a,b,c,d));
    if (length(iu)==1) {
      // State space system, with iu but w/o freq. vector
      w = dfrqint(a,b,c,d,Ts,30);
    else			// Assume this is the old syntax dbode(A,B,C,D,IU,W)
      disp("Warning: You are using the old syntax.  Use dbode(A,B,C,D,Ts,IU,W) instead.");
      w = iu; 
      iu= Ts; 
      Ts=1;
    }
    ny = d.nr;
    nu = d.nc;

  else
    msg = abcdchk(a,b,c,d);
    if (msg != "") { error(msg); }
    ny = d.nr;
    nu = d.nc;

  }}}}

  if (nu*ny==0) { return <<mag=[]; phase=[];w=[]>>; }

  // Compute frequency response
  if (nargs==3||nargs==4) {
    g = freqresp(num,den,exp(sqrt(-1)*w*Ts));
  else
    g = freqresp(a,b,c,d,iu,exp(sqrt(-1)*w*Ts));
  }

  mag = abs(g);
  phase = (180/pi)*unwrap(atan2(imag(g),real(g)));

  // Uncomment the following statement if you don't want the phase to  
  // be unwrapped.  Note that phase unwrapping may not always work; it
  // is only a "guess" as to whether +-360 should be added to the phase
  // to make it more aestheticly pleasing.
  // phase = (180/pi)*imag(log(g));

  // Uncomment the following line for decibels, but be warned that the
  // MARGIN function will not work with decibel data.
  // mag = 20*log10(mag);

  return <<mag=mag;phase=phase;w=w>>;
};

