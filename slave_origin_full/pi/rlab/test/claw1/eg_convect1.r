//
// test for clawpack 1-d solver: acoustic example
//
// convection equation:
//      q_t + v(x) q_x = 0
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
q = [x,q0];

velocity = function(x)
{
  global(xlo,xhi,pi);
  return 1+2*sin(2*pi*(x-xlo)/(xhi-xlo));
};

//
// Riemann solver function: finds waves and speeds in a Hyperbolic problem
//
R = function(t,x,ql,qr)
{
  // DF
  waves = ql-qr;
  speed = velocity(x);

  // entropy fix for transonic rarefactions:
  amdq =  min(speed,0)*waves;
  apdq =  max(0,speed)*waves;
  return [waves,amdq,apdq;[speed,0,0]];
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

for (i in 1:NITER)
{
  qn = q;

  for(dt in [1/64:1:1/64])
  {
    spinner();

    //qn = claw1(,R,S,qn,[dt-0.01,dt]);
    qn = claw1(,R,,qn,[dt-0.01,dt]);
    tl = "CLAWPACK 1-D: Convection equation";
    pltitle(tl);
    //sleep(0.1);
    data.[1] = q;
    data.[2] = qn;
    plwin(1);
    plimits(,,-1.5,1.5);
    plplot( data );
  }

  i?
}


