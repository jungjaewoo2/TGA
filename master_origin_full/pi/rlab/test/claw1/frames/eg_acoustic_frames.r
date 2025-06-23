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

//
// Note: (this requires netpbm package installed)
//  1) in rlab2 run this script. It should generate 255 ppm images. Exit rlab2.
//  2) In console type:
//      > ppmtompeg mpegpar.txt
//  3) Use kaffeine to watch the movie acoustic.mpg
//
gnuclose(gnuwins().available);
plwins (1);
if (!exist(NITER))
{ NITER = 50; }

xlo = -1;
xhi =  1;

nx=512;
dx = (xhi-xlo)/nx;
x =(xlo:xhi:dx)';
q0=exp(-(x-.05).^2/0.1);
q = [x,q0,zeros(q0)];


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
R = function(t,x,ql,qr){
    // DF
    global(iedf,vedf,speed);
    aDF   = iedf * (ql-qr);
    waves = vedf .* aDF';
    return [ waves ; speed ];
    };
//
// Source function
//
S = function(t,x,q){
    return -q;
    };
//
// before step one function: see manual
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

//
// prepare frames directory
//
mkdir("./frames");


x=1;
data=<<>>;

clawparams.bc(2,1,3);       // 2 layers, left absorbing, right reflecting
clawparams.tsrc("adams");   // use adams' to integrate the problem
clawparams.b4step1(bfsone); // use function q<-bfsone(t,x,q) to modify q

qn = q;
i = 0;

for(dt in [1/64:4:1/64])
{
  spinner();
  //qn = claw1(,R,S,qn,[dt-1/64,dt]);
  qn = claw1(K,R,S,qn,[dt-1/64,dt]);
  tl = "CLAWPACK 1-D: Acoustic Equation with damping";
  //
  // first plot
  //
  pltitle(tl);
  data.[1] = q[;1,2];
  data.[2] = q[;1,3];
  data.[3] = qn[;1,2];
  data.[4] = qn[;1,3];
  plimits(,,-1.5,1.5);
  plwid([5,5,5,5]);
  plegend(["p(x,0)", "u(x,0)", "p(x,t)", "u(x,t)" ] );
  plot( data );
  psfname="./frames/"+gsub("0", " ",text(i,"%4g")).string+".ppm";
  plcopy(1,psfname+"/ppm");
  plegend(["p(x,0)", "u(x,0)", "p(x,t)", "u(x,t)" ] );
  plwid([5,5,5,5]);
  plot( data );
  plclose();
  i=i+1;
}



