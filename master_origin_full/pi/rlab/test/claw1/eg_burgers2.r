//
// test for clawpack 1-d solver:
//
// burgers equation:
//      q_t + (0.5*q^2)_x = 0
// with
//  BC: rigid wall on the right, free on the left
//

gnuclose(gnuwins().available);
plwins (1);
if (!exist(NITER))
{ NITER = 50; }

xlo =  0;
xhi =  1;

nx=512;
dx = (xhi-xlo)/nx;
x =(xlo:xhi:dx)';
q0=0.5*sin(2*pi*x);
q = [x,q0,q0];


//
// Riemann solver function: finds waves and speeds in a Hyperbolic problem
//
R = function(t,x,ql,qr)
{
  // DF
  waves = diag(ql-qr);
  speed = 0.5*(qr+ql)';

  // entropy fix for transonic rarefactions:
  if (qr[1] < 0 && ql[1] > 0)
  {
    amdq =  -0.5*qr.^2;
    apdq =   0.5*ql.^2;
    return [waves,amdq,apdq;[speed,0,0]];
  }

  return [ waves ; speed ];
};

//
// Source function
//
S = function(t,x,q)
{
  return -q;
};


x=1;
data=<<>>;

for( i in 1:NITER)
{
  qn = q;

  for(dt in [1/16:1:1/16])
  {
    spinner();
    //qn = claw1(,R,S,qn,[dt-1/64,dt]);
    qn = claw1(,R,,qn,[dt-1/16,dt]);
    tl = "CLAWPACK 1-D: Burgers equation";
    plwin (1);
    pltitle(tl);
    data.[1] = q;
    data.[2] = qn;
    plimits(,,-1.5,1.5);
    plplot( data );
  }

  i?
}


