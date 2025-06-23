//----------------------------------------------------------------------
//
// pzmap
//
// syntax: </p;z/> = pzmap(num,den) 
//         </p;z/> = pzmap(a,b,c,d)
//
//	Plot pole-zero map of continuous-time linear system.
//	pzmap(A,B,C,D) computes the eigenvalues and transmission zeros of
//	the continuous-time state-space system (A,B,C,D) and plots them in
//	the complex s-plane.  The poles are plotted as x's and the zeros 
//	are plotted as o's.  
//
//	PZMAP(NUM,DEN) computes the poles and zeros of the SISO polynomial
//	transfer function G(s) = NUM(s)/DEN(s) where NUM and DEN contain 
//	the polynomial coefficients in descending powers of s.  If the 
//	system has more than one output, then the transmission zeros are 
//	computed.
//
//	pzmap(P,Z) plots the poles, P, and the zeros, Z, in the complex 
//	plane.  P and Z must be column vectors.
//
//		</P;Z/> = pzmap(NUM,DEN)  or  </P;Z/> = pzmap(A,B,C,D)
//
//	returns the poles and transmission zeros of the system in the 
//	column vectors P and Z.  No plot is drawn on the screen.  
//
//	The function sgrid or zgrid can be used to plot lines of constant
//	damping ratio and natural frequency in the s or z plane.
//
//	See also: rlocus, sgrid, zgrid, eig, tzero, ss2zp, and tf2zp.
//
//----------------------------------------------------------------------
require abcdchk nargchk roots tf2ss tfchk tzero

pzmap=function(a,b,c,d)
{
  global(eps,pi,_rlab_config,_rloc_window)

  msg = nargchk(2,4,nargs);
  if (msg!="")  { error(msg); }
  if (nargs==3) { error("Wrong number of input arguments."); }

  // --- Determine which syntax is being used ---
  if (nargs==2) {
    nd = b.nr; 
    md = b.nc;
    if (md<=1) {
      // Assume Pole-Zero form
      p = a; 
      z = b;
    else		
      // Transfer function form
      tmp = tfchk(a,b);
      num = tmp.numc;
      den = tmp.denc;
      p = roots(den);
      nn = num.nr; 
      mn = num.nc;
      if (nn==1) {
        z = roots(num);
      else
        </a;b;c;d/> = tf2ss(num,den);
        z = tzero(a,b,c,d).z;
      }
    }

  else			
    // State space system 
    msg = abcdchk(a,b,c,d);
    if (msg != "") { error(msg); }
    p = eig(a).val[:];
    z = tzero(a,b,c,d).z;

  }

  // plot graph
  if (_rlab_config.plot_support=="pgplot"||_rlab_config.plot_support=="plplot")
  {
    if (exist(_rloc_window))
    {
       plwin(_rloc_window);
    else
       _rloc_window = plstart();
    }
  }
  cros = <<[-1,1;1,-1]; [1,1;-1,-1]>>;
  circ = [cos(0:2*pi:0.4)[:], sin(0:2*pi:0.4)[:]];
  mep = max([eps;abs([real(p);real(z);imag(p);imag(z)])]);
  if (z.n==p.n)
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
  for (i in 1:p.n)
  {
      xy.[j]   = cros.[1] + [real(p[i]),imag(p[i]);real(p[i]),imag(p[i])];
      xy.[j+1] = cros.[2] + [real(p[i]),imag(p[i]);real(p[i]),imag(p[i])];
      j = j + 2;
  }
  for (i in 1:z.n)
  {
      xy.[j] = [circ[;1] + real(z[i]), circ[;2] + imag(z[i])];
      j = j + 1;
  }
  
  if (_rlab_config.plot_support == "pgplot") 
  {
      plcolor(ones(1,14)*7);    // optional, here's for pgplot
      plline(ones(1,20));       // optional, here's for pgplot
  }
  plegend();
  plimits(-ax,ax,-ax,ax);
  plwid(4);
  pltitle("Poles and Zeros Map");
  plstyle("line");
  ylabel("Imag Axis");
  xlabel("Real Axis");
  plot(xy);
  if (_rlab_config.plot_support == "pgplot") 
  {
    plcolor();
    plline();
  }
  plimits();
  pltitle();
  xlabel("");
  ylabel("");
  
  return <<pout=p; z=z>>;
  
};

