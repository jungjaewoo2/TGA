//----------------------------------------------------------------------
//
// minreal
//
// syntax: </am;bm;cm;dm/> = minreal(a,b,c,d,tol)
//         </pm;zm/> = minreal(z,p,tol)
//         </DENm;NUMm/> = ninreal(NUM,DEN,tol)
//
// Minimal realization and pole-zero cancellation.
// </Am;Bm;Cm;Dm/> = minreal(A,B,C,D) returns a minimal realization
// of the state-space system {A,B,C,D}.  A message is displayed 
// indicating the number of states removed.
// </Am;Bm;Cm;Dm/> = minreal(A,B,C,D,TOL) uses the tolerance TOL
// in deciding which states to eliminate.
//
// </Pm;Zm/> = minreal(Z,P), where Z and P are column vectors
// containing poles and zeros, cancels the common roots that
// are within TOL = 10*SQRT(EPS)*ABS(Z(i)) of each other.
// </Pm;Zm/> = minreal(Z,P,TOL) uses tolerance TOL.
//
// For transfer functions, </DENm;NUMm/> = minreal(NUM,DEN), where
// NUM and DEN are row vectors of polynomial coefficients, cancels
// the common roots in the polynomials.
// </DENm;NUMm/> = minreal(NUM,DEN,TOL) uses tolerance TOL.
//
//----------------------------------------------------------------------
require ctrbf obsvf tf2zp zp2tf

minreal=function(a,b,c,d,tol)
{
  global(eps,pi)
  
  if (nargs == 2 || nargs == 3) {
     z = a;
     p = b;
     mz = z.nr; nz = z.nc;
     mp = p.nr; np = p.nc;
     if (mz == 1 && mp == 1) {
	// If transfer function, convert to zero-pole:
	tmp = tf2zp(z,p);
	z = tmp.z;
	p = tmp.p;
	k = tmp.k;
     }

     // Strip infinities from zeros and throw away.
     z = z[find(finite(z))];

     mz = max(size(z));
     mp = max(size(p));
     iz = ones(mz,1);

     // Loop through zeros, looking for matching poles:
     for (i in 1:mz) {
       zi = z[i];
       if (nargs == 2) {
          tol = 10*abs(zi)*sqrt(eps);
       else
          tol=c;
       }
       kk = find(abs(p-zi) <= tol);
       if (all(size(kk))) {
          p[kk[1]] = [];
          iz[i] = 0;
       }
     }

     // Eliminate matches in zeros:
     z = z[iz];

     // If transfer function, convert back.
     if (nz > 1 || np > 1) {
        tmp = zp2tf(z,p,k);
        z = tmp.num;
        p = tmp.den;
     }
     disp(int2str(mz-sum(iz))+" pole-zeros cancelled");
     return  <<num=z; den=p>>;
  }

  // Do state-space case
  ns = b.nr; nu = b.nc;
  if (nargs == 4) {
     tol = 10*ns*norm(a,"1")*eps;
  }
  tmp = ctrbf(a,b,c,tol);
  am = tmp.abar;
  bm = tmp.bbar;
  cm = tmp.cbar;
  t  = tmp.t;
  k  = tmp.k;
  kk = sum(k);
  nu = ns - kk;
  nn = nu;
  am = am[nu+1:ns;nu+1:ns];
  bm = bm[nu+1:ns;];
  cm = cm[;nu+1:ns];
  ns = ns - nu;
  if (ns) {
     tmp = obsvf(am,bm,cm,tol);
     am = tmp.abar;
     bm = tmp.bbar;
     cm = tmp.cbar;
     t  = tmp.t;
     k  = tmp.k;
     kk = sum(k);
     nu = ns - kk;
     nn = nn + nu;
     am = am[nu+1:ns;nu+1:ns];
     bm = bm[nu+1:ns;];
     cm = cm[;nu+1:ns];
  }
  disp(int2str(nn)+" states removed");
  
  return <<am=am;bm=bm;cm=cm;dm=d>>;
};

