//----------------------------------------------------------------------
//
// nyquist
//
// syntax: </ im; re; w /> = nyquist(a,b,c,d,iu)
//         </ im; re; w /> = nyquist(a,b,c,d,iu,w)
//         </ im; re; w /> = nyquist(num,den)
//         </ im; re; w /> = nyquist(num,den,w) 
//
// Nyquist frequency response for continuous-time linear systems.
// nyquist(A,B,C,D,IU) produces a Nyquist plot from the single input
// IU to all the outputs of the system:             
//        .                                    -1
//        x = Ax + Bu             G(s) = C(sI-A) B + D  
//        y = Cx + Du      RE(w) = real(G(jw)), IM(w) = imag(G(jw))
//
// The frequency range and number of points are chosen automatically.
//
// nyquist(NUM,DEN) produces the Nyquist plot for the polynomial 
// transfer function G(s) = NUM(s)/DEN(s) where NUM and DEN contain
// the polynomial coefficients in descending powers of s. 
//
// nyquist(A,B,C,D,IU,W) or nyquist(NUM,DEN,W) uses the user-supplied
// freq. vector W which must contain the frequencies, in radians/sec,
// at which the Nyquist response is to be evaluated.  When invoked 
// with left hand arguments,
//	</IM;RE;W/> = nyquist(A,B,C,D,...)
//	</IM;RE;W/> = nyuist(NUM,DEN,...) 
// returns the frequency vector W and matrices RE and IM with as many
// columns as outputs and length(W) rows.  No plot is drawn on the 
// screen.
//
// See also: logspace, margin, bode, and nichols.
//
//----------------------------------------------------------------------
require abcdchk freqint2 freqresp tfchk

nyquist = function(a,b,c,d,iu,w)
{

  global(eps,pi,_ctb2_window,_rlab_config)
  
  // --- Determine which syntax is being used ---
  if (nargs==1) {
    error("Wrong number of input arguments.");
  
  else if (nargs==2) {
    // Transfer function form without frequency vector
    </den;num/> = tfchk(a,b);
    w   = freqint2(num,den,30);
    ny  = num.nr; 
    nn  = num.nc;
    nu  = 1;

  else if (nargs==3) {	// Transfer function form with frequency vector
    </den;num/> = tfchk(a,b);
    w   = c;
    ny  = num.nr;
    nn  = num.nc;
    nu  = 1;

  else if (nargs==4) {
    // State space system, w/o iu or frequency vector
    msg = abcdchk(a,b,c,d);
    if (msg!="") { error(msg); }
    w  = freqint2(a,b,c,d,30);
    ny = d.nr; 
    nu = d.nc;   
    if (ny*nu > 1) {
       //MIMO system 
       iu = 0;
       re = []; 
       im = [];
       for (i in 1:nu) {
         tmp = $self(a,b,c,d,i,w);
         re  = [re,tmp.re]; 
         im  = [im,tmp.im];
         pause();
       }
    else
       //SISO system
       iu = 1; 
       nargs = 5;
    }
    
    
    if (!iu) { 
        return <<im=im;re=re;w=w>>;
    }

  else if (nargs==5) {
    // State space system, with iu but w/o freq. vector
    msg = abcdchk(a,b,c,d);
    if (msg!="") { error(msg); }
    w  = freqint2(a,b,c,d,30);
    ny = d.nr;
    nu = d.nc;

  else
    msg = abcdchk(a,b,c,d);
    if (msg!="") { error(msg); }
    ny = d.nr;
    nu = d.nc;

  }}}}}

  if (nu*ny==0) { return <<re=[]; im=[]; w=[]>>; }

  // Compute frequency response
  if ((nargs==2)||(nargs==3)) {
    g = freqresp(num,den,sqrt(-1)*w);
  else
    g = freqresp(a,b,c,d,iu,sqrt(-1)*w);
  }
  re = real(g); 
  im = imag(g);

  // plot graph.
  if (_rlab_config.plot_support=="pgplot"||_rlab_config.plot_support=="plplot")
  {
    if (exist(_ctb2_window))
    {
       plwin(_ctb2_window);
    else
       _ctb2_window = plstart();
    }
  }
  // Make arrows
  sx = max(re) - min(re);   
  sy = max(abs(2*im));
  sample = maxi(abs(2*im));
  arrow = [-1;0;-1]*sx/50 + sqrt(-1)*[1;0;-1]*sy/50;
  sample = sample + (sample==1);
  reim = diag(g[sample;]);
  d = diag(g[sample+1;]-g[sample-1;]); 
  rot = sign(d);  // Arrow is always horizonatl i.e. -1 or +1
  // rot=d./abs(d)  Use this when arrow is not horizontal
  xy =ones(3,1)*reim' + ones(3,1)*rot'.*arrow;
  xy2=ones(3,1)*reim' - ones(3,1)*rot'.*arrow;
  m = g.nr;
  n = g.nc; 
  // mark -1 + j0
  cros = <<[-1,1;1,-1]; [1,1;-1,-1]>>;
  mep = max([eps;abs(re);abs(im)]);
  ax  = 1.2*mep;
  cros.[1] = cros.[1]*ax/50 + [-1,0;-1,0];
  cros.[2] = cros.[2]*ax/50 + [-1,0;-1,0];
  if (_rlab_config.plot_support == "pgplot") 
  {  
     plcolor(ones(1,14)*7);  // optional, here's for pgplot
     plline(ones(1,20));     // optional, here's for pgplot
  }
  pltitle("Nyquist Plot");
  xlabel("Real Axis"); 
  ylabel("Imag Axis");
  plegend();

  if (n==1) { 
     plot(<<[re,im];[re,-im]; [real(xy),-imag(xy)]; [real(xy2),imag(xy2)];cros.[1];cros.[2]>>);
  else
     plot(<<[re,im];[re,-im]; [real(xy),-imag(xy)];[real(xy2),imag(xy2)];cros.[1];cros.[2]>>);
  }

  // clean up
  xlabel("");
  ylabel("");
  pltitle("");
  if (_rlab_config.plot_support == "pgplot")
  {
     plcolor();  // optional, here's for pgplot
     plline(); 
  }
  plegend("default");
  
  return <<re=re; im=im; w=w>>; 
};

