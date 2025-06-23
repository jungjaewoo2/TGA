//
// test for bacol 1-d solver:
//
//  u_t = u_xx + pi*pi*sin(pi*x)    (1)
//
// x in [0,1], with boundary conditions u(0,t) = 1 and u(1,t) = 1
// and initial condition u(x,0) = 1.

npde = 1;

PDESOLVER = "bacol";

NITER = 1;

//
// f
//
f = function(t,x,u,ux,uxx)
{
  //printf("t = %g, x = %g, u = %g, ux = %g, uxx = %g\n",t,x,u,ux,uxx);
  global(pi);
  rval = uxx + pi*pi*sin(pi*x);
  return rval;
};

//
// Df
//
Df = function(t,x,u,ux,uxx)
{
  //printf("df: t = %g, x = %g, u = %g, ux = %g, uxx = %g\n",t,x,u,ux,uxx);
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
bL = function(t,u,ux)
{
  //printf("bL: t = %g, u = %g, ux = %g\n",t,u,ux);
  // Left boundary fixed: u = 0
  // B(u) = u
  rval = u;
  return rval;
};
DbL = function(t,u,ux)
{
  //printf("DbL: t = %g, u = %g, ux = %g\n",t,u,ux);
  // Left boundary fixed: u = 0
  // B(u) = u
  // Z(t) = 0
  rval=<<>>;
  rval.const = 1;
  rval.dbdu  = 1;
  rval.dbdux = 0;
  rval.dbdt  = 0;
  return rval;
};
bR = function(t,u,ux)
{
  //printf("bR: t = %g, u = %g, ux = %g\n",t,u,ux);
  //printf("t = %g\n",t);
  // Right boundary fixed: u = 0
  // B(u) = u
  rval = u;
  return rval;
};
DbR = function(t,u,ux)
{
  //printf("DbR: t = %g, u = %g, ux = %g\n",t,u,ux);
  // Right boundary fixed: u = 0
  // B(u) = u
  // Z(t) = 0
  global(pi);
  rval=<<>>;
  rval.dbdu  = 1;
  rval.dbdux = 0;
  rval.dbdt  = 0;
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
options.stdout = stderr();
options.erel   = 1e-1;
options.eabs   = 1e-1;

//
// resolution of the functions
//
sizeS = [32];
xlo = 0;
xhi = 1;

tic();

nx = sizeS[ 1 ];
dx = (xhi-xlo)/nx;
x  = (xlo:xhi:dx)';

wn = [];

for (i in 2:length(T))
{
  wn = [wn, bacol(npde,f,Df,bL,bR,DbL,DbR,x,w0,T[1,i],options)];
}
// for (i in 1:NITER)
// {
//   spinner();
//   wn = bacol(npde,f,Df,bL,bR,DbL,DbR,x,w0,T,options);
// }

wmax = max(max(wn));
wmin = min(min(wn));

exno = 2;

rfile module_plot

