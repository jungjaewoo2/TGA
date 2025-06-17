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
  global(COEFF);
  rval = COEFF*uxx - u * ux;
  return rval;
};

//
// Df
//
Df = function(t,x,u,ux,uxx)
{
  global(COEFF);
  //printf("df: t = %g, x = %g, u = %g, ux = %g, uxx = %g\n",t,x,u,ux,uxx);
  //t
  rval=<<>>;
  rval.dfdu   = -ux;
  rval.dfdux  = -u;
  rval.dfduxx =  COEFF;
  return rval;
};

//
// DbL,R, the boundary conditions on x=xL and x=xR
//
bL = function(t,u,ux)
{
  global(COEFF);
  //printf("bL: t = %g, u = %g, ux = %g\n",t,u,ux);
  rval = u - 0.5 + 0.5 * tanh((-0.5 * t - 0.25)/(4.0 *COEFF));
  return rval;
};
DbL = function(t,u,ux)
{
  //printf("DbL: t = %g, u = %g, ux = %g\n",t,u,ux);
  global(COEFF);
  rval=<<>>;
  rval.dbdu  = 1;
  rval.dbdux = 0;
  rval.dbdt  = -1.0 / COEFF / sinh(-(0.5 * t + 0.25)/(4*COEFF)).^2;
  return rval;
};
bR = function(t,u,ux)
{
  global(COEFF);
  //printf("bR: t = %g, u = %g, ux = %g\n",t,u,ux);
  //printf("t = %g\n",t);
  rval = u - 0.5 + 0.5 * tanh((1.0 - 0.5 * t -0.25)/(4.0 * COEFF));
  return rval;
};
DbR = function(t,u,ux)
{
  //printf("DbR: t = %g, u = %g, ux = %g\n",t,u,ux);
  global(COEFF);
  rval=<<>>;
  rval.dbdu  = 1;
  rval.dbdux = 0;
  rval.dbdt  = -1 / COEFF / sinh(-(0.5 * t - 0.75)/(4*COEFF)).^2;
  return rval;
};

w0 = function(x)
{
  global(COEFF);
  rval = -0.5 -0.5 * tanh((x - 0.25)/(4*COEFF));
  return rval;
};


COEFF = 1e-3;

//
// time interval for integration
//
dt = 1/32;
ti = 0;
tf = 1;
T = [0:tf:dt];

options=<<>>;
options.stdout = stderr();
options.atol = 1e-3;
options.rtol = 1e-3;
options.idir = 1; //

//
// resolution of the functions
//
xlo = 0;
xhi = 1;
x  =[xlo:xhi:1/64]';

tic();

wn = [];

for (i in 2:length(T))
{
  spinner();
  wn = [wn, bacol(npde,f,Df,bL,bR,DbL,DbR,x,w0,T[1,i],options)];
}

wmax = max(max(wn));
wmin = min(min(wn));

exno = 1;

rfile module_plot

