//----------------------------------------------------------------------
//
// sigma
//
// Syntax: </svout;w/> = sigma(a,b,c,d,w,invflag)
//
//	Singular value frequency response of continuous linear systems.
//	sigma(A,B,C,D) produces a singular value plot of the complex 
//	matrix                  -1
//                G(jw) = C(jwI-A) B + D 
//	as a function of frequency.  The singular values are an extension
//	of Bode magnitude response for MIMO systems.  The frequency range
//	and number of points are chosen automatically.  For square 
//	systems, sigma(A,B,C,D,'inv') produces the singular values of the
//	inverse complex matrix
//	     -1               -1      -1
//           G (jw) = [ C(jwI-A) B + D ]
//
//	sigma(A,B,C,D,W) or sigma(A,B,C,D,W,"inv") uses the user-supplied
//	frequency vector W which must contain the frequencies, in 
//	radians/sec, at which the singular value response is to be 
//	evaluated. When invoked with left hand arguments,
//		</SV;W/> = sigma(A,B,C,D,...)
//	returns the frequency vector W and the matrix SV with as many 
//	columns	as min(NU,NY) and length(W) rows, where NU is the number
//	of inputs and NY is the number of outputs.  No plot is drawn on 
//	the screen.  The singular values are returned in descending order.
//
//----------------------------------------------------------------------
require abcdchk sigma2 freqint

sigma = function(a,b,c,d,w,invflag)
{
  local(a,b,c,d,w,invflag)
  
  if (nargs==6) {
    // Trap call to RCT function
    if (!isstring(invflag)) {
      svout = sigma2(a,b,c,d,w,invflag);
      return <<svout=svout;w=w>>;
    }
    if (!length(invflag)) {
	nargs = nargs - 1;
    }
  }
  if (nargs < 4 || nargs > 6) {
     error("Wrong number of input arguments.");
  }
  
  msg=abcdchk(a,b,c,d);
  if (msg!="") { error(msg); }
  
  // Detect null systems
  if (!(length(d) || (length(b) &&  length(c)))) {
       return <<svout=[]; w=w>>;
  }

  // Determine status of invflag
  if (nargs==4) { 
    invflag = [];
    w = [];
  } else { if (nargs==5) {
    if (isstring(w)) {
      invflag = w;
      w = [];
      ny = d.nr; nu = d.nc;
      if (ny!=nu) { 
         error("The state space system must be square when using 'inv'."); 
      }
    } else {
      invflag = [];
    }
  } else {
    ny = d.nr; nu = d.nc;
    if (ny!=nu) { 
        error("The state space system must be square when using 'inv'.");
    }
  }}

  // Generate frequency range if one is not specified.

  // If frequency vector supplied then use Auto-selection algorithm
  // Fifth argument determines precision of the plot.
  if (!length(w)) {
    w=freqint(a,b,c,d,30);
  }

  nx = a.nr; na = a.nc;
  no = c.nr; ns = c.nc;
  nw = max(size(w));

  // Balance A
  tmp = balance(a);
  t = tmp.t;
  a = tmp.ab;
  b = t \ b;
  c = c * t;

  // Reduce A to Hessenberg form:
  tmp = hess(a);
  p = tmp.p;
  a = tmp.h;

  // Apply similarity transformations from Hessenberg
  // reduction to B and C:
  b = p' * b;
  c = c * p;

  s = w * sqrt(-1);
  I=eye(length(a),length(a));
  if (nx > 0) {
    for (i in 1:length(s)) {
      if (isempty(invflag)) {
        sv[;i] = svd(c*((s[i]*I-a)\b) + d).sigma[:];
      } else {
        sv[;i] = svd(inv(c*((s[i]*I-a)\b) + d)).sigma[:];
      }
    }
  } else {
    for (i in 1:length(s)) {
      if (isempty(invflag)) {
        sv[;i] = svd(d).sigma[:];
      } else {
        sv[;i] = svd(inv(d)).sigma[:];
      }
    }
  }
  sv = sv'; 
  
  // plot graph.
  #  //semilogx
  #  xlabel("Frequency (rad/sec)");
  #  ylabel("Singular Values dB");
  #  plot(<<[w,20*log10(sv)];[w,zeros(w.n,1)]>>);
  
  return <<svout=sv; w=w[:]>>;

};



