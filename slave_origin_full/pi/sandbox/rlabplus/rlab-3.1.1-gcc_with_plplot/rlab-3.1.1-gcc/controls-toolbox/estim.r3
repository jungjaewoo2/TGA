//----------------------------------------------------------------------
//
// estim
//
// syntax: </ae;be;ce;de/> = estim(a,b,c,d,l,e,f)
//
// Form continuous Kalman estimator.
// </Ae;Be;Ce;De/> = estim(A,B,C,D,L) produces the Kalman estimator 
// based on the continuous system (A,B,C,D) with Kalman gain matrix
// L assuming all the outputs of the system are sensor outputs.  The
// resulting state-space estimator is
//	  .
//	 xHat  = [A-LC] xHat + Ly
//
//	|yHat| = |C| xHat
//	|xHat|   |I|
//
// and has estimated sensors yHat and estimated states xHat as 
// outputs, and sensors y as inputs.  
//
// </Ae;Be;Ce;De/> = estim(A,B,C,D,L,SENSORS,KNOWN) forms the Kalman 
// estimator using the sensors specified by SENSORS, and the 
// additional known inputs specified by KNOWN.  The resulting system
// has estimated sensors and states as outputs, and the known inputs
// and sensors as inputs.  The KNOWN inputs are non-stochastic inputs
// of the plant and are usually control inputs.
//
// See also: reg, dcontrol, destim, lqr, dlqr, lqe and dlqe.
//
//----------------------------------------------------------------------
require abcdchk cloop ssselect

estim = function(a,b,c,d,l,e,f)
{
  
  if ((nargin!=5)&&(nargs!=7)) {
     error("Wrong number of input arguments."); 
  }
  error(abcdchk(a,b,c,d));

  nx = a.nr; na = a.nc;
  ny = d.nr; nu = d.nc;

  if (nargs==5) {
    sensors = [1:ny]; known = [];
  }
  if (nargs==7) {
    sensors = e; known = f;
  }

  nsens = length(sensors); nknown = length(known);

  // Check size of L with number of states and sensors
  nl = l.nr; ml = l.nc;
  if (ml!=nsens) {
    error("Number of sensors and size of L matrix don't match."); 
  }
  if (nl!=nx) {
    error("The A and L matrices must have the same number of rows."); 
  }

  inputs = [1:nsens] + nu; states = [1:nx] + ny;

  // Form continuous Kalman estimator
  ae = a;
  be = [b,l];
  ce = [c;eye(nx,nx)];
  de = [d, zeros(ny,nsens);zeros(nx,nu), zeros(nx,nsens)];

  // close sensor feedback loop
  </ae;be;ce;de/> = cloop(ae,be,ce,de,sensors,-inputs);
  return ssselect(ae,be,ce,de,[known,inputs],[sensors,states]);
};

