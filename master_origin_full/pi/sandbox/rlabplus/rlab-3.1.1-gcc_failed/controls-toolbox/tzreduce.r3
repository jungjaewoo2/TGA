//---------------------------------------------------------------------
//
// tzreduce
//
// Syntax: </abcd;mu;nu/> = tzreduce(abcd,m,n,p,Zeps,ro,sigma)
//
// Extract from the system (A,B,C,D) the reduced system 
// (a,b,c,d) but with d full row rank.  Used in tzero.
//
// Extracts from the (N+P)x(M+N) system [ B A ],  a (NU+MU)x(M+NU) 'reduced' 
//        [ B' A' ]                     [ D C ]
// system [ D' C' ] having the same transmission zeros but with D' Full Row 
// Rank.  The system [ A' B' C' D' ] overwrites the old system.
//
// Reference: "Computation of Zeros of Linear Multivariable
//            Systems", A. Emami-Naeini, and P. Van Dooren; Automatica
//            Vol. 18, No. 4, pp. 415-430, 1982.
//
//---------------------------------------------------------------------
require housh

tzreduce=function(abcd,m,n,p,Zeps,ro,sigma)
{
  local(abcd,ro,sigma)
  
  Sum = zeros(1,max(p,m));
  mu = p;
  nu = n;
  while (mu!=0)
  {
    ro1 = ro;
    mnu = m + nu;
    numu= nu + mu;

    if (m!=0) {
      ro1 = ro1 + 1;
      irow = nu;
      // *** Compress rows of D(*).  First exploit triangular shape ***
      for (icol in 1:sigma-1) {
        rows = irow + [1:ro1];
        dummy = abcd[rows;icol];
        </s;dummy;zero/> = housh(dummy,1,Zeps);
        abcd[rows;icol:mnu] = (eye(ro1,ro1)-s*dummy*dummy')*abcd[rows;icol:mnu];
        irow = irow + 1;
      }

      // *** Continue householder with pivoting ***

      if (sigma==0) {
        sigma = 1;
        ro1 = ro1 - 1;
      }

      if (sigma!=m) {
        if (ro1==1) {
          Sum[sigma:m] = abcd[irow+1;sigma:m].*abcd[irow+1;sigma:m];
        } else {
          Sum[sigma:m] = sum(abcd[irow+1:irow+ro1;sigma:m].*abcd[irow+1:irow+ro1;sigma:m]);
        }
      }

      for (icol in sigma:m) {
        // *** Pivot if necessary ***
        if (icol!=m) {
          Rows = 1:numu;
          dum  = max (Sum[[1,icol+1:m]]);
          ibar = maxi(Sum[[1,icol+1:m]]);
          ibar = ibar+icol-1;
          if (ibar!=icol) {
            Sum[ibar]=Sum[icol]; 
            Sum[icol]=dum;
            dum=abcd[Rows;icol];
            abcd[Rows;icol]=abcd[Rows;ibar];
            abcd[Rows;ibar]=dum;
          }
        }

        // *** Perform Householder transformation ***

        dummy=abcd[irow+1:irow+ro1;icol];
        </s;dummy;zero/> = housh(dummy,1,Zeps);
        if (zero) { break; }
        if (ro1==1) { return <<abcd=abcd; mu=mu; nu=nu>>; }
        abcd[irow+1:irow+ro1;icol:mnu] = (eye(ro1,ro1)-s*dummy*dummy')*abcd[irow+1:irow+ro1;icol:mnu];
        irow = irow+1;
        ro1 = ro1-1;
        Sum[icol:m] = Sum[icol:m] - abcd[irow;icol:m] .* abcd[irow;icol:m];
      }

    }
    tau = ro1;
    sigma = mu-tau;

    // *** Compress the columns of C(*) ***
    if (nu<=0) {
       mu = sigma; 
       nu = 0; 
       return <<abcd=abcd; mu=mu; nu=nu>>;
    }

    i1 = nu + sigma;
    mm1= m+1;
    n1 = nu;
    if (tau!=1) {
      if (mm1==mnu) {
        Sum[1:tau] = (abcd[i1+1:i1+tau;mm1].*abcd[i1+1:i1+tau;mm1])';
      } else {
        Sum[1:tau] = sum((abcd[i1+1:i1+tau;mm1:mnu].*abcd[i1+1:i1+tau;mm1:mnu])');
      }
    }

    for (ro1 in 1:tau) {
      ro = ro1 - 1;
      i  = tau - ro;
      i2 = i + i1;

      % *** Pivot if necessary ***

      if (i!=1) {
        dum  = max (Sum[1:i]);
        ibar = maxi(Sum[1:i]);
        if (ibar!=i) {
          Sum[ibar] = Sum[i]; Sum[i] = dum;
          dum = abcd[i2;mm1:mnu];
          abcd[i2;mm1:mnu] = abcd[ibar+i1;mm1:mnu];
          abcd[ibar+i1;mm1:mnu] = dum;
        }
      }

      // *** Perform Householder Transformation ***
    
      cols = m + [1:n1];
      dummy = abcd[i2;cols];
      </s;dummy;zero/> = housh(dummy,n1,Zeps);
      if (zero) { break; }
      if (n1!=1) {
        abcd[1:i2;cols] = abcd[1:i2;cols]*(eye(n1,n1)-s*dummy*dummy');
        mn1 = m + n1;
        abcd[1:n1;1:mn1] = (eye(n1,n1)-s*dummy*dummy')*abcd[1:n1;1:mn1];
        Sum[1:i] = Sum[1:i]-abcd[i1+1:i1+i;mn1]' .* abcd[i1+1:i1+i;mn1]';
        mnu = mnu - 1;
      }
      n1 = n1 - 1;

    }

    if (!zero) { ro=tau; }

    nu = nu-ro;
    mu = sigma + ro;
  
    if (ro==0) { return <<abcd=abcd; mu=mu; nu=nu>>; }
  }
  
  return <<abcd=abcd; mu=mu; nu=nu>>;
};

