// function S = xsopt(X,Mp,g1,g2,g3,R,nx,nu,nz,nc)
//-------------------------------------------------------------------------
//
// xsopt
//
// Syntax: S=xsopt(X,Mp,g1,g2,g3,R,nx,nu,nz,nc)
//
// Computes the optimal S=-S^* matrix with respect to which the
// control cost function
//
//   V = tr(UR)
//
// is stationary, all other variables being constant, while 
// assigning the covariance X to the closed-loop system
//
// g1, g2, g3 are matrices defined so that the 
// controller
//
//  xc_dot = Ac * xc + F * z
//       u =  G * xc + H * z
//
// is given by
//
//  G = [ H G    = g1 + g2 * S * g3
//        F Ac ]
//
// nx, nu, nz, nc == dimensions of plant, actuator, sensor 
//                   and controller state vectors
//
// Originally written by L. D. Peterson
// Modified by Jeff Layton 10/14/90
// Ported to RLaB, 10/13/93
//
// Version JBL 931013
// Note: current version does not check compatibility of dimensions
//-------------------------------------------------------------------------

rfile jpinv

xsopt = function(X,Mp,g1,g2,g3,R,nx,nu,nz,nc)
{
   local(rowsg2,ns,m,Xp,Xpc,Xc,tol,i,j,k,prodi,G,H,MpXpMptr,...
         MpXpc,G0Xc,G0,H0,b1,b2,b4,b,P1,P2,P4,P,alpha,b3,P3,S)

// Get dimension of S

   rowsg2=g2.nr;
   ns=g2.nc;

// Degrees of freedom in the S matrix

   m=ns*(ns-1) / 2;

// Extract partitions of X

   Xp=X[1:nx,1:nx];
   Xpc=X[1:nx,nx+1:nx+nc];
   Xc=X[nx+1:nx+nc,nx+1:nx+nc];

// Form all the Gi and Hi matrices

   G0=xgpart(g1,nu,nz,nc);
   H0=xhpart(g1,nu,nz,nc);

   tol=max([size(g2),size(g3)])*epsilon()*norm(g2)*norm(g3);

   for (i in 1:m) {
        prodi=g2*xsi(ns,i)*g3;
        for (j in 1:(nu+nc)) {
             for (k in 1:(nz+nc)) {
                  if ( abs(prod[j,k]) <= tol ) {
                      prod[j,k]=0;
                  }
             }
        }
        G[i]=xgpart(prodi,nu,nz,nc);
        H[i]=xhpart(prodi,nu,nz,nc);
   }

// Form some constant matrices 

   MpXpMptr=Mp*Xp*Mp';
   MpXpc=Mp*Xpc;
   G0Xc=G0*Xc;

// Form the b vectors

   b1=zeros(m,1);
   b2=zeros(m,1);
   b4=zeros(m,1);
   for (i in 1:m) {
        b1[i] = trace(H0*MpXpMptr*(H[i]')*R);
        b2[i] = trace(H0*MpXpc*(G[i]')*R) + trace(H[i]*MpXpc*G0'*R);
        b4[i] = trace(G0Xc*(G[i]')*R);
   }

   b3=b2;
   b=b1+b2/2+b3/2+b4;

// Form the P matrices

   P1=zeros(m,m);
   P2=zeros(m,m);
   P4=zeros(m,m);
   for (i in 1:m) {
        for (j in 1:m) {
             P1[i;j]=trace(H[i]*MpXpMptr*(H[j]')*R);
             P2[i;j]=trace(H[i]*MpXpc*(G[j]')*R);
             P4[i;j]=trace(G[i]*Xc*(G[j]')*R);
        }
   }
   P3=P2';
   P=P1+P2+P3+P4;

// Solve alpha equation

   alpha=-jpinv(P)*b;

// Formulate the optimal S matrix

   S=zeros(ns,ns);
   for (i in 1:m) {
        S=S+alpha[i]*xsi(ns,i);
   }

   return S;
};

