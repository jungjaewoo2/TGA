//
// test for clawpack 1-d solver:
//
// burgers equation:
//      q_t + (0.5*q^2)_x = 0
// with
//  BC: absorbing boundary on the left and on the right
//

gnuclose(gnuwins().available);
plwins (1);
if (!exist(NITER))
{ NITER = 50; }

alfad=0.5;

xlo =  0;
xhi =  1;

nx=512;
dx = (xhi-xlo)/nx;
x =(xlo:xhi:dx)';
q0=ifelse(x<0.1,linterp(x,[0,0;0.1,1]),0);
q = [x,q0];


//
// Riemann solver function: finds waves and speeds in a Hyperbolic problem
//
R = function(t,x,ql,qr)
{
  // DF
  waves = ql-qr;
  speed = 0.5*(qr+ql);

  // entropy fix for transonic rarefactions:
  if (qr < 0 && ql > 0)
  {
    amdq =  -0.5*qr^2;
    apdq =   0.5*ql^2;
    return [waves,amdq,apdq;[speed,0,0]];
  }
  return [ waves ; speed ];
};

//
// Source function
//
S = function(t,x,q)
{
  global(alfad);
  return alfad*q-q.^2;
};


x=1;
data=<<>>;

for( i in 1:NITER)
{
  qn = q;
  for(dt in [1/64:12:1/64])
  {
    spinner();

    //qn = claw1(,R,S,qn,[dt-0.01,dt]);
    qn = claw1(,R,S,qn,[dt-0.01,dt]);
    tl = "CLAWPACK 1-D: Burgers equation";
    plwin (1);
    pltitle(tl);
    //sleep(0.1);
    data.[1] = q;
    data.[2] = qn;
    plimits(,,-0.5,3);
    plot( data );
  }

  i?
}


