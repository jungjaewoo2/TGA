//--------------------------------------------------------------------------
// dimpulse
//
// Syntax: R = dimpulse(a,b,c,d,iu,n)
//         R = dimpulse(num,den,n)
//         where R.y = y state time history
//               R.x = x state time history
//               R.n = number of points
//
//	 Impulse response of discrete-time linear systems.
//	 dimpulse(A,B,C,D,IU)  plots the response of the discrete system:
//
//		x[n+1] = Ax[n] + Bu[n]
//		y[n]   = Cx[n] + Du[n]
//
//	to an unit sample applied to the single input IU.  The number of
//	points is determined automatically.
//
//	dimpulse(NUM,DEN)  plots the impulse response of the polynomial
//	transfer function  G(z) = NUM(z)/DEN(z)  where NUM and DEN contain
//	the polynomial coefficients in descending powers of z.
//
//	dimpulse(A,B,C,D,IU,N) or dimpulse(NUM,DEN,N) uses the user-
//	supplied number of points, N.  
// 
//	It returns the output and state time history in the matrices R.y and 
//      R.x.  R.y has as many columns as there are outputs and R.x has as 
//      many columns as there are states.
//
//	See also: dstep,dinitial,dlsim and impulse.
//--------------------------------------------------------------------------
require abcdchk dlsim dtimvec tf2ss tfchk

dimpulse = function(a,b,c,d,iu,n)
{
  global(eps,pi)
  
  if (nargs <= 2 || nargs >= 6) {
     error("Wrong number of input arguments.");
  }
  if (nargs==2) {
    // Transfer function without number of points
    tmp = tfchk(a,b);
    num = tmp.numc;
    den = tmp.denc;
    iu = 1;
    </a;b;c;d/> = tf2ss(num,den);

  } else { if (nargs==3) {
    // Transfer function with number of points
    tmp = tfchk(a,b);
    num = tmp.numc;
    den = tmp.denc;
    n = c;
    iu = 1;
    </a;b;c;d/> = tf2ss(num,den);

  } else { if (nargs>=4) {
    msg = abcdchk(a,b,c,d);
    if (msg != "") { error(msg); }

  }}}

  ny = d.nr; nu = d.nc;
  if (nu*ny==0) { 
     return <<yout=[];x=[];n=[]>>;
  }
     
  // Work out number of points if not supplied.
  if (nargs==4 || nargs==5 || nargs==2) {
    if (isempty(a)) {
      n = 3;
    } else {
      // The next line controls the number of samples in the plot if N not specified
      st=0.001; // Set settling time bound  = 0.1//
      n = b.nr; m = b.nc;
      x0=(eye(n,n)-a)\(b*ones(m,1));

      // Cater for pure integrator case
      infind  = find(!finite(x0));
      x0[infind]= ones(size(infind));
      n=dtimvec(a,b,c,x0,st);
    }
    if (nargs==4) {
        // Multivariable systems
	# [iu,nargin,y]=mulresp('dimpulse',a,b,c,d,n,nargout,0);
        if (ny*nu > 1) {
          //MIMO system 
          iu = 0;
          y  = []; 
          x  = [];
          for (i in 1:nu) {
            tmp = $self(a,b,c,d,i,n);
            y   = [y,tmp.y]; 
            x   = [x,tmp.x];
          }
        } else {
          //SISO system
          iu = 1; 
          nargs = 5;
        }	
	if (!iu) { return << yout=y; x=x; n=n >>; }
    }
  }
  
  if (nargs <= 3) {
    // Transfer function description
    G = dlsim(num,den,[1;zeros(n-1,1)]); // More efficient: uses FILTER
  } else {
    if (!isempty(b)) { b=b[;iu]; } 
    d = d[;iu];
    G = dlsim(a,b,c,d,[1;zeros(n-1,1)]);
  }

  #if nargout==0,		// With no output arguments, plot graph
  #  status = inquire('hold');
  #  stairs([0:n-1],y)
  #  hold on
  #  plot([0,n-1],[0;0],':')
  #  xlabel('No. of Samples'), ylabel('Amplitude')
  #
  #  if ~status, hold off, end	// Return hold to previous status
  #  return // Suppress output
  #end
  
  return <<y=G.y; x=G.x; n=n>>;
  
};

