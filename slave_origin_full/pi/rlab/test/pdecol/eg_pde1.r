//
// example 1 from section 7 of sincovec and madsen
//

PDESOLVER = "pdecol/epdcol";

//
// f
//
f = function(t,x,u,ux,uxx)
{
  rval = [...
      u[2]*u[2]*uxx[1] - u[1]*u[2] - u[1]*u[1] + 10.0 + 2.0*u[2]*ux[2]*ux[1];...
      u[1]*u[1]*uxx[2] + u[1]*u[2] - u[2]*u[2] + uxx[1] + 2.0*u[1]*ux[1]*ux[2] ];
  return rval;
};

//
// Df
//
Df = function(t,x,u,ux,uxx)
{
  rval=<<>>;
  rval.dfdu  =  [-u[2] -2.0*u[1],2.0*u[2]*uxx[1] - u[1] + 2.0*ux[2]*ux[1];...
      2.0*u[1]*uxx[2] + u[2] + 2.0*ux[1]*ux[2],u[1] - 2.0*u[2]];
  rval.dfdux = [2.0*u[2]*ux[2], 2.0*u[2]*ux[1];2.0*u[1]*ux[2],2.0*u[1]*ux[1]];
  rval.dfduxx= [u[2]*u[2], 0; 1, u[1]*u[1]];
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
  rval.const = 1;
  rval.dbdu  = [1,0;0,1];
  rval.dbdux = [0,0;0,0];
  rval.dzdt  = [0;0];
  return rval;
};

DbR = function(t,u,ux)
{
  // Right boundary oscillatory
  // B(u) = sin(u[1]*u[2])
  // Z(t) = 0
  global(pi);
  rval=<<>>;
  rval.dbdu  = [u[2]*cos(u[1]*u[2]), u[1]*cos(u[1]*u[2]); u[2]*sin(u[1]*u[2]), u[1]*sin(u[1]*u[2])];
  rval.dbdux = [1,0;0,1];
  rval.dzdt  = [0;0];
  return rval;
};

w0 = function(x)
{
  global(pi);
  rval = [0.5*(x + 1.0); pi];
  return rval;
};


//
// time interval for integration
//
dt = 1e-3;
ti = 0;
tf = 0.5;
T = [ti:tf:dt];

options=<<>>;
options.stdout = rconsole();
options.dt     = 1e-7;
options.imethod= 1;
options.erel   = 1e-4;

//
// resolution of the functions
//
sizeS = [128, 256, 512, 128, 1024];
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

  exno = 1;

  rfile module_plot
}

