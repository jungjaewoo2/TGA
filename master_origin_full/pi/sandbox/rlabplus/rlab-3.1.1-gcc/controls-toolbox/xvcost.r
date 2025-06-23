//-------------------------------------------------------------------------
//
// V = xvcost(X,Mp,GG,R,nu,nz,nc)
//
// Computes the control cost function
//
//   V = tr(UR)
//
// Which assigns the covariance X to the closed-loop system with
// the feedback matrix GG
//
//  xc_dot = Ac * xc + F * z
//       u =  G * xc + H * z
//
//  GG = [ H G  
//         F Ac ]
//
// nu, nz, nc == dimensions of actuator, sensor 
//               and controller state vectors
//
// L.D. Peterson
// version 900305
//
// Note: current version does not check compatibility of dimensions
//-------------------------------------------------------------------------

rfile xgpart
rfile xhpart

xvcost = function(X,Mp,GG,R,nu,nz,nc)
{
   local(G,H,V)

// Get the G and H parts of GG

   G=xgpart(GG,nu,nz,nc);
   H=xhpart(GG,nu,nz,nc);

// Calculate the cost and return

   V=trace([H*Mp,G]*X*[H*Mp,G]'*R);

   return V;
};

