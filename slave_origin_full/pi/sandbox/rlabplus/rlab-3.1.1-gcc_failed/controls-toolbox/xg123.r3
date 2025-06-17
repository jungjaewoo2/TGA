//-----------------------------------------------------------------------------
//
// xg123
//
// Syntax: R=xg123(X,A,B,D,M,W,X0)
//
// This routine computes the 3 matrices that define the feedback
// gains for the Covariance Control Formulation.
//
// The input is:
// X = covariance to be assigned
// [A,B,D,M] = augmented dynamic matrices
// W = applied initial disturbances
// X0 = applied initial conditions
//
// Originally written by Lee D. Peterson, modified by Jeff Layton and
// Ported to RLaB by Jeff Layton
// Version 931026
//
// Note: current version does not check compatibility of dimensions
//-----------------------------------------------------------------------------

rfile jpinv

xg123 = function(X,A,B,D,M,W,X0)
{
   local(nx,ny,L,Q,Term1,Lplus,Xi,P,Gamma,G1,G2,G3)

// Get plant dimensions
   nx=X.nr;
   ny=X.nc;

// Form defined matrices

   L=(eye(nx,nx)-jpinv(M)*M)*inv(X)*B*jpinv(B);
   Q=X*A'+A*X+D*W*D'+X0;
   Term1=(eye(nx,nx)-jpinv(M)*M);
   Lplus=jpinv(L);
   Xi=inv(X);
   P=2.0*Lplus*Term1*Xi*Q+(eye(nx,nx)-Lplus*Term1*Xi)*Q*Lplus*L;
   Gamma=jpinv(M)*M*X*(eye(nx,nx)-B*jpinv(B));

// Form the three gain coefficient matrices and return

   G1=-0.5*jpinv(B)*Q*(2*eye(nx,nx)-B*jpinv(B))*...
      inv(X)*jpinv(M)+...
      0.5*jpinv(B)*(P'-P)*B*jpinv(B)*inv(X)*jpinv(M);

   G2=jpinv(B)*X*jpinv(M)*M*(eye(nx,nx)-Gamma*jpinv(Gamma));

   G3=(eye(nx,nx)-Gamma*jpinv(Gamma))*jpinv(M);

   return << G1=G1; G2=G2; G3=G3 >>;
};

