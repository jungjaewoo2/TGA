//
// test for clawpack 1-d solver:
//
//  acoustic example
//
//  p_t + u_x = 0
//  u_t + p_x = 0
//
// with and without source term on [-1,1]
// BC: Rigid wall on the right, free on the left

gnuclose(gnuwins().available);
plwins (1);
if (!exist(NITER))
{ NITER = 50; }

//
// Riemann solver: DF is a constant so these can be calculated
// once for all calculation
//
DF = [0,1;1,0];
eigDF = eig(DF);
speed = real(eigDF.val);
vedf  = real(eigDF.vec);
iedf  = inv(vedf);

//
// Riemann solver function: finds waves and speeds in a Hyperbolic problem
//
R = function(t,x,ql,qr)
{
  // DF
  //x
  global(iedf,vedf,speed);
  aDF   = iedf * (ql-qr);
  waves = vedf .* aDF';
  return [ waves ; speed ];
};
//
// Source function
//
S = function(t,x,q){
    //return zeros(q);
    return -q;
    };
//
// Source function
//
bfsone = function(t,x,q){
    return q;
    };
//
// Capacity function
//
K = function(x){
    return 2.-x.^2;
    };

data=<<>>;

clawparams.bc(2,1,3);       // 2 layers, left absorbing, right reflecting
clawparams.tsrc("rk2");     // use rk2 to integrate the problem
clawparams.b4step1(bfsone); // use function q<-bfsone(t,x,q) to modify q


//
// time interval for integration
//
dt = 1/64;
ti = 0;
tf = 1;
T = [ti:tf:dt];

//
// resolution of the functions
//
sizeS = [64, 256, 512, 128, 1024];
xlo = -1;
xhi =  1;

tic();
for(i in 1:NITER)
{
  is = int(sizeS.nc*rand()+1);
  while (is>sizeS.nc)
  {
    is = is - sizeS.nc;
  }
  nx = sizeS[ is ];
  dx = (xhi-xlo)/nx;
  x  =(xlo:xhi:dx)';
  q0 = exp(-(x-.05).^2/0.1);
  q  = [x,q0,zeros(q0)];

  qn = q;

  r = rand();

  for(j in 1:(T.nc-1))
  {
    spinner();
    if (r <= 0.33)
    {
      qn = claw1(K,R,S,qn,[T[j],T[j+1]]);
      tl = "CLAWPACK 1-D: Acoustic Equation - with S and K";
    else if (r <= 0.66)
    {
      qn = claw1(K,R, ,qn,[T[j],T[j+1]]);
      tl = "CLAWPACK 1-D: Acoustic Equation - no S, with K";
    else
      qn = claw1( ,R,S,qn,[T[j],T[j+1]]);
      tl = "CLAWPACK 1-D: Acoustic Equation - with S, no K";
    }}
    plwin (1);
    pltitle(tl);
    plegend(["p(x,0)", "u(x,0)", "p(x,t)", "u(x,t)" ] );
    plformat(["with lines using 1:2", "with lines using 1:3", ...
        "with lines using 1:2", "with lines using 1:3"]);
    data.[1] = q;
    data.[2] = q;
    data.[3] = qn;
    data.[4] = qn;
    plimits(xlo,xhi,-1.5,1.5);
    plplot( data );
  }
}

printf("Calculation lasted %g sec.\n", toc());


