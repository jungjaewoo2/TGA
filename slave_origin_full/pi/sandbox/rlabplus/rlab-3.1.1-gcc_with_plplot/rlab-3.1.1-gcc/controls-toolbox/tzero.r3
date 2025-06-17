//---------------------------------------------------------------------
// 
// tzero
//
// syntax:  </ gain; z /> = tzero(A,B,C,D)
//
//	tzero(A,B,C,D) returns the transmission zeros z of the state-
//	space system:   .
//                       x = Ax + Bu   or   x[n+1] = Ax[n] + Bu[n]
//                       y = Cx + Du          y[n] = Cx[n] + Du[n]
//
//	also returns the transfer function gain if the system is SISO.
//
//	See also: pzmap, ss2zp, and zp2ss.
//
// Extracts from the system matrix of a state-space system [ A B C D ] 
// a regular pencil [ Lambda*Bf - Af ] which has the NU Invariant Zeros 
// of the system as generalized eigenvalues.
//  
//  Reference: Adapted from "Computation of Zeros of Linear Multivariable
//             Systems", A. Emami-Naeini, and P. Van Dooren; Automatica 
//             Vol. 18, No. 4, pp. 415-430, 1982.
//
//---------------------------------------------------------------------
require housh tzreduce

tzero=function(a,b,c,d)
{
  global(eps,pi)
  
  if (nargs!=4) { error("Wrong number of arguments"); }
  msg = abcdchk(a,b,c,d);
  if (msg!="") { error(msg); }

  if (isempty(b) && isempty(c)) {
     return <<z=[]; gain=0>>;
  }

  Zeps = 2*eps*norm(a,"f");  % Epsilon used for transmission zero calculation.

  nn = a.nr; nx = a.nc;
  pp = c.nr; na = c.nc;
  na = b.nr; mm = b.nc;
  if (pp*mm==1) { tf=1; } else { tf=0; }

  // *** Construct the Compound Matrix [ B A ] of Dimension (N+P)x(M+N) 
  //                                   [ D C ]

  bf = [b,a;d,c];

  // *** Reduce this system to one with the same invariant zeros and with
  //     D(*) full rank MU (The Normal Rank of the original system) ***

  </bf;mu;nu/> = tzreduce(bf,mm,nn,pp,Zeps,pp,0);
  Rank=mu;

  if (nu==0) { 
    z = [];
  } else {
    //  *** Pretranspose the system *** 
    mnu = mm + nu;
    numu = nu + mu;
    af = zeros(mnu,numu);
    af[mnu:1:-1;numu:1:-1] = bf[1:numu;1:mnu]';

    if (mu!=mm) {
      pp=mm;
      nn=nu;
      mm=mu;
      // *** Reduce the system to one with the same invariant zeros and with 
      //     D(*) square invertable *** 
      </af;mu;nu/> = tzreduce(af,mm,nn,pp,Zeps,pp-mm,mm);
      mnu=mm+nu;
    }

    if (nu==0) {
      z = [];
    } else {
      // *** Perform a unitary transformation on the columns of [ sI-A B ] in 
      //                           [ sBf-Af X ]                 [   -C D ] 
      //     order to reduce it to [   0    Y ] with Y & Bf square invertable 
      bf[1:nu;1:mnu]=[zeros(nu,mm),eye(nu,nu)];
      if (Rank!=0) {
        nu1 = nu + 1;
        i1 = nu + mu;
        i0 = mm;
        for (i in 1:mm) {
          i0  = i0-1;
          cols = i0 + [1:nu1];
          dummy = af[i1;cols];
          </s;dummy;zero/> = housh(dummy,nu1,Zeps);
          af[1:i1;cols] = af[1:i1;cols]*(eye(nu1,nu1)-s*dummy*dummy');
          bf[1:nu;cols] = bf[1:nu;cols]*(eye(nu1,nu1)-s*dummy*dummy');
          i1 = i1-1;
        }
      }
      // Solve Generalized zeros of sBF - AF
      z = eig(af[1:nu;1:nu]/bf[1:nu;1:nu]).val[:];
    }
  }

  // Compute transfer function gain
  if (nu==nx) {
    gain=bf[nu+1;1];
  } else {
    gain=bf[nu+1;1]*prod(diag(bf[nu+2:nx+1;nu+2:nx+1]));
  }
  
  return <<z=z;gain=gain>>;
};

