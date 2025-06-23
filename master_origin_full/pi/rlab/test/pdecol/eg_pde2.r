//
// test for pdecol 1-d solver:
//
//  u_t = u_xx + pi*pi*sin(pi*x)    (1)
//
// x in [0,1], with boundary conditions u(0,t) = 1 and u(1,t) = 1
// and initial condition u(x,0) = 1.

PDESOLVER = "pdecol/epdcol";

npde = 1;

//
// f
//
f = function(t,x,u,ux,uxx)
{
  global(pi);
  rval = uxx + pi*pi*sin(pi*x);
  return rval;
};

//
// Df
//
Df = function(t,x,u,ux,uxx)
{
  //t
  rval=<<>>;
  rval.const  =  1;
  rval.dfdu   =  0;
  rval.dfdux  =  0;
  rval.dfduxx =  1;
  return rval;
};

//
// DbL,R, the boundary conditions on x=xL and x=xR
//
DbL = function(t,u,ux)
{
  // Left boundary fixed: u = 0
  // B(u) = u
  // Z(t) = 0
  rval=<<>>;
  rval.const = 1;
  rval.dbdu  = 1;
  rval.dbdux = 0;
  rval.dzdt  = 0;
  return rval;
};

DbR = function(t,u,ux)
{
  // Right boundary fixed: u = 0
  // B(u) = u
  // Z(t) = 0
  global(pi);
  rval=<<>>;
  rval.const = 1;
  rval.dbdu  = 1;
  rval.dbdux = 0;
  rval.dzdt  = 0;
  return rval;
};

w0 = function(x)
{
  global(pi);
  rval = -x*(1-x);
  return rval;
};


//
// time interval for integration
//
dt = 1/256;
ti = 0;
tf = 1;
T = [ti:tf:dt];

options=<<>>;
options.stdout = rconsole();
options.dt     = 1e-6;
options.imethod= 2;
options.erel   = 1e-4;

//
// resolution of the functions
//
sizeS = [128];
NITER = 1;
xlo = 0;
xhi = 1;

tic();

nx = sizeS[ 1 ];
dx = (xhi-xlo)/nx;
x  = (xlo:xhi:dx)';

wn = pdecol(npde,f,Df,DbL,DbR,x,w0,T,options);

if (!exist(NOPLOTS))
{
  wmax = max(max(wn));
  wmin = min(min(wn));

  exno = 2;

  rfile module_plot
}
