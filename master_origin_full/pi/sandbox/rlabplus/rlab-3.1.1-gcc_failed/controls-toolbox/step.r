//----------------------------------------------------------------------
// step
//
// syntax: </t;x;y/> = step(a,b,c,d,iu)
//         </t;x;y/> = step(a,b,c,d,iu,t)
//         </t;x;y/> = step(num,den)
//         </t;x;y/> = step(num,den,t)
//
// Step response of continuous-time linear systems.
// step(A,B,C,D,IU)  plots the time response of the linear system:
//	.
//	x = Ax + Bu
//	y = Cx + Du
//
// to a step applied to the input IU.  The time vector is auto-
// matically determined.  step(A,B,C,D,IU,T) allows the specification
// of a regularly spaced time vector T.
//
// </T;X;Y/> = step(A,B,C,D,IU) returns the output and state time 
// response in the matrices Y and X  respectively.  The matrix Y has 
// as many columns as there are outputs, and LENGTH(T) rows.  The 
// matrix X has as many columns as there are states.  If the time 
// vector is not specified, then the automatically determined time 
// vector is returned in T.
//
// </T;X;Y/> = step(NUM,DEN) calculates the 
// step response from the transfer function description 
// G(s) = NUM(s)/DEN(s) where NUM and DEN contain the polynomial 
// coefficients in descending powers of s.
//
// See also: initial, impulse, lsim and dstep.
//
//----------------------------------------------------------------------
require abcdchk c2d checklist ltitr nargchk tf2ss tfchk timvec

step = function(a,b,c,d,iu,t)
{
  global(_ctb2_window, _rlab_config)
  local (a,b,c,d,iu,t)

  msg = nargchk(2,6,nargs);
  if (msg!="") { error(msg); }

  if (nargs < 4) {
     // Convert transfer function to state space
     </den;num/> = tfchk(a,b);
     if (nargs==3) { t = c; }
     </a;b;c;d/> = tf2ss(num,den);
     iu = 1;
     nargs = nargs + 3;
  else
     msg = abcdchk(a,b,c,d);
     if (msg!="") {error(msg);}
  }

  ny = d.nr;
  nu = d.nc;
  if (nu*ny==0) {
     return <<t=[];x=[];y=[]>>; 
  }

  if (nargs>4) {
     // only iu-th input related items needed 
     if (!isempty(b)) { b = b[;iu]; }
     d = d[;iu];
  }

  // Workout time vector if not supplied.
  if (nargs==5 || nargs==4) {
    if (isempty(a)) {
      t = (0:1:.1)';
    else
      // The next two constants control the precision of the plot
      // and the time interval of the plot.
      st = 0.005; // Set settling time bound  = 0.5%
      precision = 30; // Show approx 30 points for simple graph
      // Step response is effectively equal to placing initial conditions
      // on the plant as follows:
      n = b.nr;
      m = b.nc;
      x0 = -a\(b*ones(m,1));
      // Cater for pure integrator case
      for (i in 1:x0.nc) {
        infind = find(!finite(x0[;i]));
        x0[infind;i] = ones(length(infind),1);
      }
      t = timvec(a,b,c,x0,st,precision);
    }
  }


  %  Multivariable systems
  if (nargs==4) {
     error("Sorry STEP won't able to deal with multiple input variable at this time");
  }


  // Simulation
  dt = t[2] - t[1];
  </bb;aa/> = c2d(a,b,dt);
  n  = t.n;
  nb = b.nr;
  mb = b.nc;
  x = ltitr(aa,bb,ones(n,1),zeros(nb,mb));
  if (isempty(a)) {
    x = [];
    y = ones(n,1)*d.';
  else
    y = x*c.'+ ones(n,1)*d.';
  }

  // plot graph   
  if (_rlab_config.plot_support=="pgplot"||_rlab_config.plot_support=="plplot")
  {
    if (exist( _ctb2_window )) 
    {
       plwin( _ctb2_window );
    else
       _ctb2_window = plstart();
    }
  }
  // matrix division next line often gives ill-condition error
  // any remedy ?
  dcgain = -c/a*b + d;

  // filter out INF/NAN, such as resonance
  P = checklist([t,y]);
  P.[y.nc+1] = [t[1],dcgain';t[n],dcgain'];
  plegend();
  xlabel("Time (secs)");
  ylabel("Amplitude");
  plot(P);
  xlabel("");
  ylabel("");
  plegend("default");
  
  return <<t=t;x=x;y=y>>;
};

