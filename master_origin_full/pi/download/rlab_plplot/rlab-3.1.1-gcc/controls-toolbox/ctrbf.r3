//--------------------------------------------------------------------------
//
// ctrbf
//
// Syntax: R = ctrbf(a, b, c, tol)
//         R is a list with components: R.Abar, R.Bbar, R.Cbar, R.T, R.k
//
//	Controllability staircase form.
//	R = CTRBF(A,B,C) returns a decomposition
//	into the controllable/uncontrollable subspaces.
//	R = CTRBF(A,B,C,TOL) uses tolerance TOL.
//
//	If Co=CTRB(A,B) has rank r <= n, then there is a similarity
//	transformation T (=R.T) such that
//
//	R.Abar = T * A * T' ,  R.Bbar = T * B ,  R.Cbar = C * T'
//
//	and the transformed system has the form
//
//	         | Anc    0 |             | 0 |
//	R.Abar =  ----------  ,  R.Bbar =  ---  ,  R.Cbar = [Cnc| Cc].
//	         | A21   Ac |             |Bc |
//	                                           -1          -1
//	where (Ac,Bc) is controllable, and Cc(sI-Ac)Bc = C(sI-A)B.
//
//	See also: CTRB.
//
// This file implements the Staircase Algorithm of Rosenbrock, 1968.
//--------------------------------------------------------------------------

require rank

ctrbf = function(a, b, c, tol)
{
  global (eps)
  
  ra = a.nr; ca = a.nc;
  rb = b.nr; cb = b.nc;
  //
  // initial conditions :
  //
  ptjn1 = eye(ra,ra);
  ajn1 = a;
  bjn1 = b;
  rojn1 = cb;
  deltajn1 = 0;
  sigmajn1 = ra;
  k = zeros(1,ra);
  if (nargs == 3) { tol = ra*norm(a,"1")*eps; }
 
  //
  // Major Loop 
  //
  for (jj in 1:ra)
  {
    tmp = svd(bjn1);
    uj = tmp.u;
    sj = diag(tmp.sigma);
    vj = tmp.vt';
    //[uj,sj,vj] = svd(bjn1);
    rsj = sj.nr; csj = sj.nc;
    p = rot90(eye(rsj,rsj),1);
    uj = uj*p;
    bb = uj'*bjn1;
    roj = rank(bb,tol);
    rbb = bb.nr; cbb = bb.nc;
    sigmaj = rbb - roj;
    sigmajn1 = sigmaj;
    k[jj] = roj;
    if (roj == 0) { break; }
    if (sigmaj == 0) { break; }
    abxy = uj' * ajn1 * uj;
    aj   = abxy[1:sigmaj;1:sigmaj];
    bj   = abxy[1:sigmaj;sigmaj+1:sigmaj+roj];
    ajn1 = aj;
    bjn1 = bj;
    ruj = uj.nr; cuj = uj.nc;
    ptj = ptjn1 * ...
          [uj, zeros(ruj,deltajn1); ...
           zeros(deltajn1,cuj), eye(deltajn1,deltajn1)];
    ptjn1 = ptj;
    deltaj = deltajn1 + roj;
    deltajn1 = deltaj;
  }
  //
  // transformation
  //
  t = ptjn1';
  abar = t * a * t';
  bbar = t * b;
  cbar = c * t';

  return <<abar=abar;bbar=bbar;cbar=cbar;t=t;k=k>>;
};

