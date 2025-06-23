//----------------------------------------------------------------------
//
// c2dt
//
// syntax: </ gamma; h; j; phi /> = c2dt(A,B,C,T,lambda)
//
// Conversion of continuous state space models to discrete 
// models with pure time delay in the inputs.
// </Bd;Cd;Dd;Ad/> = C2DT(A,B,C,T,lambda) converts the continuous time
// system
//      .
//      x(t) = Ax(t) + Bu(t-lambda)
//      y(t) = Cx(t) + Du(t)
//
// to the discrete system with sample time T,
//
//      x(k+1) = Ad x(k) + Bd u(k)
//      y(k) = Cd x(k) + Dd u(k) 
//
// See also: pade
//
// From results in Franklin and Powell, Digital Control; Addison
// and Wesley; 1980 chapter 6, pages 171-177, with extension for
// multivariable controls and outputs.
//----------------------------------------------------------------------
require expm

c2dt = function(A,B,C,T,lambda)
{
  ns = B.nr; nc = B.nc;
  no = C.nr; ns = C.nc;
  l = ceil(lambda/T);
  m = l*T - lambda;

  s1 = expm([A*m, B*m; zeros(nc,ns+nc)]);
  s2 = expm([A*(T-m), B*(T-m); zeros(nc,ns+nc)]);
  s3 = eye((l-1)*nc,(l-1)*nc);
  s4 = zeros(ns,(l-2)*nc);
  s5 = zeros((l-1)*nc,ns+nc);
  s6 = zeros(nc,ns+l*nc);

  PHI1 = s1[1:ns;1:ns]*s2[1:ns;1:ns];
  GAMMA1 = s1[1:ns;1:ns]*s2[1:ns;ns+1:ns+nc];
  GAMMA2 = s1[1:ns;ns+1:ns+nc];
  if (l == 0)
  {
    // This is the modified z-transform case
    PHI = PHI1;
    GAMMA = PHI1 * GAMMA2 + GAMMA1;
    H     =  C;
    J     =  C * GAMMA2;
    //
    //	delay less than one period 
    //
  } else { if (l == 1)
  {
    PHI = [PHI1, GAMMA1; s6];
    GAMMA  = [GAMMA2;eye(nc,nc)];
    H      =  [C,zeros(no,nc)];
    J      =  zeros(no,nc);
  } else {
    PHI   =  [PHI1,GAMMA1,GAMMA2,s4;s5,s3;s6];
    GAMMA =  [zeros(ns+(l-1)*nc,nc);eye(nc,nc)];
    H     =  [C,zeros(no,l*nc)];
    J      =  zeros(no,nc);
  }}

  return <<phi=PHI; gamma=GAMMA; h=H; j=J>>;
};


