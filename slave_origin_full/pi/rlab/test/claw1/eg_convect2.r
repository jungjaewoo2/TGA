//
// test for clawpack 1-d solver:
//
// convection equation:
//      q_t + (v(x)q)_x = 0
// with
//  BC: absorbing boundary on both sides
//
// OBSERVE: Riemann solver is the same as for eg_convect1.r,
// the difference is in source term. Conservation law is observed
// easily here.
//

gnuclose(gnuwins().available);
plwins (1);
if (!exist(NITER))
{ NITER = 50; }

xlo =  0;
xhi =  1;

nx=512;
dx = (xhi-xlo)/nx;
x  = (xlo:xhi:dx)';
q0 = 0.5*sin(2*pi*x);
q  = [x,q0];

velocity = function(x)
{
  global(xlo,xhi,pi);
  return 2+2*sin(2*pi*(x-xlo)/(xhi-xlo));
};

dvdx = function(x)
{
  global(xlo,xhi,pi);
  return 2*pi/(xhi-xlo)*2*cos(2*pi*(x-xlo)/(xhi-xlo));
};

//
// Riemann solver function: finds waves and speeds in a Hyperbolic problem
//
R = function(t,x,ql,qr)
{
  global(dx);

  // DF
  waves = ql-qr;
  speed = velocity(x);
  return [waves;speed];
};

//
// Source function
//
S = function(t,x,q)
{
  return -q;
};


data=<<>>;

for (i in 1:NITER)
{
  qn = q;

  for(dt in [1/64:0.25:1/64])
  {
    spinner();
    //qn = claw1(,R,S,qn,[dt-0.01,dt]);
    qn = claw1(,R,S,qn,[dt-1/64,dt]);
    tl = "CLAWPACK 1-D: Convection equation";
    pltitle(tl);
    //sleep(0.1);
    data.[1] = q;
    data.[2] = qn;
    plwin (1);
    plimits(,,-1.5,1.5);
    plplot( data );
  }

  i?
}


