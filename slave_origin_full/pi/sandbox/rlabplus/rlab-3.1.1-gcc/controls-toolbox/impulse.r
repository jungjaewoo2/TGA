//----------------------------------------------------------------------
//
// impulse
//
// syntax: </t;x;yout/> = impulse(a,b,c,d,iu,t)
//         </t;x;yout/> = impluse(num,den,t)
//
// Impulse response of continuous-time linear systems.
// impulse(A,B,C,D,IU)  plots the time response of the linear system
//	.
//	x = Ax + Bu
//	y = Cx + Du
//
// to an impulse applied to the single input IU.  The time vector is
// automatically determined.  
//
// impluse(NUM,DEN) plots the impulse response of the polynomial 
// transfer function  G(s) = NUM(s)/DEN(s)  where NUM and DEN contain
// the polynomial coefficients in descending powers of s.
//
// impulse(A,B,C,D,IU,T) or impulse(NUM,DEN,T) uses the user-supplied
// time vector T which must be regularly spaced.  When invoked with
// left hand arguments,
//	</T;X;Y/> = impulse(A,B,C,D,...)
//	</T;X;Y/> = impulse(NUM,DEN,...)
// returns the output and state time history in the matrices Y and X.
// No plot is drawn on the screen.  Y has as many columns as there 
// are outputs and length(T) rows.  X has as many columns as there 
// are states.
//
// See also: step,intitial, lsim and dimpulse.
//
//----------------------------------------------------------------------
require c2d checklist ltitr tfchk abcdchk tf2ss timvec

impulse = function(a,b,c,d,iu,t)
{
  global(eps,pi,_rlab_config,_ctb2_window)
  

  if (nargs < 2 || nargs > 6) {
     error("Wrong number of arguments");
  }
  if (nargs < 4) {
     // Convert to state space
     T = tfchk(a,b);
     num = T.numc;
     den = T.denc;
     if (nargs==3) { t = c; }
     iu = 1;
     </a;b;c;d/> = tf2ss(num,den);
     nargs = nargs + 3;
  else 
     msg = abcdchk(a,b,c,d);
     if (msg != "") { error(msg); }
  }

  ny = d.nr;
  nu = d.nc;
  if (nu*ny==0||isempty(a)) {
     return <<yout=[]; x=[]; t=[]>>;
  }

  if (exist(iu)) {
     if (!isempty(b)) {b = b[;iu]; }
  }

  // Workout time vector if not supplied.
  if (!exist(t)) {
     // The next three constants control the precision of the plot
     // and the time interval of the plot.
     st=0.005;     // Set settling time bound  = 0.5%
     tint=1;       // Set time interval to approx.  1*st% set. time
     precision=30; // Show approx 30 points for simple graph

     m=min(size(b));
     n=b.nr;
     m=b.nc;
     if (m>1) { x0=max(abs(b.')).'; else x0=b; }
     t = timvec(a,b,c,x0,st,precision);
  }

  //  Multivariable systems
  if (nargs==4) {
     yout = [];
     x = [];
     for (i in 1:nu) {
         tmp  = $self(a,b,c,d,i,t);
         yout = [yout,tmp.yout];
         x    = [x,tmp.x];
         pause();
     }
     return <<yout=yout; x=x; t=t>>;
  }

  dt = t[2]-t[1];
  T = c2d(a,b,dt);
  aa = T.phi;
  bb = T.gamma;
  n = length(t);
  x = ltitr(aa,bb,zeros(n,1),b);
  y = x * c.';
  
  // Plot Graph
  if (_rlab_config.plot_support=="pgplot"||_rlab_config.plot_support=="plplot")
  {
    if (exist(_ctb2_window))
    {
       plwin(_ctb2_window);
    else
       _ctb2_window = plstart();
    }
  }
  P = checklist([t,y]);
  P.[y.nc+1] = [t[1],0;t[length(t)],0];
  xlabel("Time (secs)"); 
  ylabel("Amplitude");
  plot(P);
  xlabel("");
  ylabel("");
  
  return <<yout=y; x=x; t=t>>;
};
