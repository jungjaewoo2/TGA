//---------------------------------------------------------------------
//
// xc12
//
// Syntax: R=xc12(X,Ap,Bp,Dp,Mp,W,X0,nc)
//         where R.C1 = First matrix assignability condition
//               R.C2 = Second matrix assignability condition
//
// This routine computes the two matrix covariance assignability
// constraints for a given plant.
//
// The input is:
// X = covariance to be assigned
// [Ap,Bp,Dp,Mp] = plant dynamic matrices
// W = applied initial disturbances
// X0 = applied initial conditions
// nc = order of desired compensator
//
// The routine returns a list:
//
//    R.C1 = First matrix assignability condition
//    R.C2 = Second matrix assignability condition
//
// Originally written by Lee D. Peterson, modified by Jeff Layton and
// ported to RLaB by Jeff Layton.
// Version 931026
//
// Note: current version does not check compatibility of dimensions
//---------------------------------------------------------------------

rfile jpinv

xc12 = function(X,Ap,Bp,Dp,Mp,W,X0,nc)
{
   local(nx,nu,nw,nz,A,B,M,D,Xp,Xpc,C1,C2,...
         Xc,X0p,X0c,T,u,s,D0,Dp0,Qp,Xpb,Xpbi,Qpb,L);

// Get plant dimensions

   nx=Ap.nr;
   nu=Bp.nc;
   nw=Dp.nc;
   nz=Mp.nr;

// Form up output feedback matrices

   A=[Ap,zeros(nx,nc);zeros(nc,nx),zeros(nc,nc)];

   B=[Bp,zeros(nx,nc);zeros(nc,nu),eye(nc,nc)];

   M=[Mp,zeros(nz,nc);zeros(nc,nx),eye(nc,nc)];

   D=[Dp;zeros(nc,nw)];

// Extract partitions of X

   Xp = X[1:nx;1:nx];
   Xpc = X[1:nx;nx+1:nx+nc];
   Xc = X[nx+1:nx+nc;nx+1:nx+nc];

// Extract partitions of X0

   X0p = X0[1:nx;1:nx];
   X0c = X0[nx+1:nx+nc;nx+1:nx+nc];

// Form D0
   T=svd(D*W*D'+X0);
   u=T.u;
   s=diag(T.sigma);
   D0 = u*sqrt(s);
   clear(u,s);
   Dp0=D0[1:nx;];

// Form defined matrices

//Qp = Xp*Ap' + Ap*Xp + Dp0*Dp0';
   Qp=Xp*Ap'+Ap*Xp+Dp*W*Dp'+X0p;
   Xpb=Xp-Xpc*inv(Xc)*Xpc';
   Xpbi=inv(Xpb);
   Qpb=Xpbi*(Xpb*Ap'+Ap*Xpb+Dp0*Dp0'+...
       Xpc*inv(Xc)*X0c*inv(Xc)*Xpc')*inv(Xpb);
   clear(Xpbi,Dp0,Xpb);
   L=(eye(nx+nc,nx+nc)-jpinv(M)*M)*inv(X)*B*jpinv(B);

// Form the three assignability conditions and return

   C1=(eye(nx,nx)-Bp*jpinv(Bp))*Qp*(eye(nx,nx)-Bp*jpinv(Bp));

   C2=(eye(nx,nx)-jpinv(Mp)*Mp)*Qpb*(eye(nx,nx)-jpinv(Mp)*Mp);

   return << C1=C1; C2=C2>>;
};

