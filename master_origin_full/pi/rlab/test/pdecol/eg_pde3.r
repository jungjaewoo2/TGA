//
// example 3
//
// u_tt = u_xx
// with boundary conditions
//    u(0) = 0
//    u(1) = sin(pi*t)
// and initial condition
//    u(x) = 0
//
// use w = [u, u_t]
// then
// w_t = [ w[2]; -w[1]_xx ]

PDESOLVER = "pdecol/epdcol";


//
// f
//
f = function(t,x,u,ux,uxx)
{
  rval = [ u[2]; uxx[1] ];
  return rval;
};

//
// Df
//
Df = function(t,x,u,ux,uxx)
{
  rval=<<>>;
  rval.const  = 1;
  rval.dfdu   = [0,1; 0,0];
  rval.dfdux  = [0,0; 0,0];
  rval.dfduxx = [0,0; 1,0];
  return rval;
};

//
// DbL,R, the boundary conditions on x=xL and x=xR
//
DbL = function(t,u,ux)
{
  // Left boundary fixed: w = 0
  // B(u) = [u(1);u(2)]
  // Z(t) = 0
  rval=<<>>;
  rval.const  = 1;
  rval.dbdu  = [1,0;0,1];
  rval.dbdux = [0,0;0,0];
  rval.dzdt  = [0;0];
  return rval;
};

DbR = function(t,u,ux)
{
  // Right boundary oscillatory
  // B(w) = w
  // Z(t) = [sin(pi*t); pi*cos(pi*t)]
  global(pi);
  rval=<<>>;
  rval.dbdu  = [1,0;0,1];
  rval.dbdux = [0,0;0,0];
  rval.dzdt  = [pi*cos(pi*t);-pi^2*sin(pi*t)];
  return rval;
};

w0 = function(x)
{
  global(pi);
  rval = [0;0];
  return rval;
};


//
// time interval for integration
//
dt = 1e-2;
ti = 0;
tf = 4;
T = [ti:tf:dt];

options=<<>>;
options.stdout = rconsole();
options.dt     = 1e-7;
options.imethod= 2;
options.erel   = 1e-6;

//
// resolution of the functions
//
sizeS = [128];
NITER = 1;

xlo = 0;
xhi = 1;

tic();

npde = 2;
nx = sizeS[ 1 ];
dx = (xhi-xlo)/nx;
x  = (xlo:xhi:dx)';

wn = pdecol(npde,f,Df,DbL,DbR,x,w0,T,options);

if (!exist(NOPLOTS))
{
  wmax = max(max(wn));
  wmin = min(min(wn));

  exno = 3;

  rfile module_plot_mpg
}

