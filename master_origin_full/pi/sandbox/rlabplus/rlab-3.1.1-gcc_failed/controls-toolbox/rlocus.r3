// ---------------------------------------------------------------------
// rlocus, rlocus_plot:  Evans root locus.
//
// Syntax: R = rlocus(num, den)          transfer function
//         R = rlocus(num, den, k)       transfer function
//         R = rlocus(A, B, C, D)        state space
//         R = rlocus(A, B, C, D, k)     state space
//
//         rlocus_plot(R, title)         plotting
//
// Description:
//
// rlocus(num,den) calculates the locus of the roots of:
//                              num(s)
//                      1 + k ----------  =  0
//                              den(s)
// for a set of gains k which are adaptively calculated to produce a 
// smooth plot using rlocus_plot. Alternatively the vector k can be 
// specified with an optional right-hand argument R=rlocus(num,den,k). 
// Vectors num and den must contain the numerator and denominator 
// coeffiecients in descending powers of s or z.  This function returns
// a list R, 
//    R.rout  - length(k) rows and (length(den)-1) columns containing 
//              the complex root locations.  Each row of the matrix 
//              corresponds to a gain from vector k.  
//    R.k     - the gains k are returned.
//    R.poles - a vector containing the poles
//    R.zeros - a vector containing the zeros
// 
// rlocus(A, B, C, D) finds the root-locus from the equivalent SISO 
// state space system (A,B,C,D):       
//                  .
//                  x = Ax + Bu    u = -k*y
//                  y = Cx + Du
//
// rlocus_plot(R, title) plots the root locus calculated by rlocus.
// The argument 'title' is an optional string for graph title.
//
// Example:
//
//   num = [1,3];
//   den = [2, 1, 1, 2, 0];
//   R = rlocus(num,den);
//   rlocus_plot(R);
//      - OR -
//   rlocus_plot( rlocus(num,den) );
//
// Note: Current plot settings are for pgplot. You need to modify
//       rlocus_plot for the plotting package used.
//
// ---------------------------------------------------------------------

static (k_known)

rlocus = function (a,b,c,d,e,f)
{

  require tfchk tf2ss abcdchk perpxy poly polyval roots tzero vsort
  global (eps)
    
  // The next variable controls how many points are plotted on the graph.
  precision = 1; //Set to a higher number (e.g. 2 for more points on the graph).
  k = [];

  if ( nargs==3 || nargs==2 ) 	
  {  
     // Convert to state space
     if (nargs==3) { k = c; }
     tmp = tfchk(a,b);
     num = tmp.numc;
     den = tmp.denc;
     </a;b;c;d/> = tf2ss(a,b,2);
  }

  ny = d.nr;
  nu = d.nc;
  if ((ny*nu==0)||isempty(a))
  {  
    k = [];
    rout = [];
    return << rout=rout; k=k >>; 
  }

  //  Multivariable systems (not ready yet)
  msg = abcdchk(a,b,c,d);
  if (msg != "") { error(msg); }
  if (nargs==5)
  {
    k=e; 
  } else { 
    if (nargs==6)
    {
      b=b[;e]; d=d[;e]; k = f;
    }
  }

  // Trap MIMO case
  ny = d.nr;
  nu = d.nc;
  if (ny*nu != 1)
  {
     error("Root locus must be SISO.");       
  }

  sp = sum([1:length(a)].^2);

  // Work out scales and ranges
  ep = eig(a).val';    // poles
  p  = length(ep);
  if (nargs>=4)
  {
     tz  = tzero(a,b,c,d).z;  // zeros (state space)
     z   = length(tz);
     sp2 = [1:z].^2;
     if (ep == [])
     {
        den = 1;
     } else {
        den = real(poly(ep));
     }
     if (tz == [])
     {
        num = 1;
     } else {
        num = real(poly(tz));
     }
  } else {
     while(num[1]==0) { num[1]=[]; }
     if (num.n == 1)
     {
        tz = [];
     } else {
        tz = roots(num)[:];     // zeros (transfer function)
     }
     z = length(tz);
  }

  if (isempty(k))
  {
     den2   = den + 1e-20*(den==0); 
     dcgain = num[length(num)]/den2[length(den)];
     // Normalize for better numerical properties
     if (abs(den[1]) > eps && den[1] != 1)
     {
        den = den/den[1];
        num = num/den[1];
     }
     mep = max([eps;abs([real(ep);real(tz);imag(ep);imag(tz)])]);
     if (z==p)
     {
        // Round up axis to units to the nearest 5
        ax = 1.2*mep;        
     } else {
        // Round graph axis    
        exponent = 5*10^(floor(log10(mep))-1);
        ax = 2*round(mep/exponent)*exponent;    
     }
     // Set up initial gain based on D.C.Gain of open loop and 
     // positions of zeros and poles.
     // Since closed loop den = num + k*den, sensitivity is related
     // to difference between num and denominator coefficients.
     diff  = abs(num)./abs(den2[1:length(num)]);
     kinit = 0.01/(abs(dcgain) + polyval(diff,4)) + 1e-12;
     bc    = b*c;
     k     = [0,1e-4*kinit,kinit]; 
     dist  = kinit;
     r[;1] = ep;
     r[;2] = vsort(ep, eig(a-bc*(k[2]./(1+k[2]*d))).val', sp).v2;
     i     = 2; 
     ij    = sqrt(-1); 
     perr  = 1;
     terminate=0; 	
     while (!terminate)
     {
       i  = i+1;
       ki = k[i]./(1+k[i]*d);
       r1 = eig(a-bc*ki).val';
       // Sort out eigenvalues so that plotting appears continuous. 
       tmp = vsort(r[;i-1],r1,sp);
       rtmp = tmp.v2;
       pind = tmp.pind;
       r[;i] = rtmp;
       // Adjust values of k based on linearity of the roots:
       // First two points in line
       y1 = imag(r[;i-2]); y2 = imag(r[;i-1]);
       x1 = real(r[;i-2]); x2 = real(r[;i-1]);
       // Current  points
       y = imag(r[;i]); x = real(r[;i]);
       // Nearest x-y co-ordinates of new point to line
       </newx;newy/> = perpxy(x1,y1,x2,y2,x,y);
       // Error estimation
       err   = sqrt((newy-y).^2+(newx-x).^2);
       distm = norm((y-y2)+ij*(x-x2),2);
       fferr = find(!finite(err));
       err[fferr] = zeros(length(fferr),1);
       perr = precision*max(err[find(finite(err))])/(distm+ax/(1e3*(distm+eps))); 

       // If percentage error greater than threshold then go back and 
       // re-evaluate more roots:
       pind = any((y==0)!=(y2==0));
       if (perr>0.2 || pind)
       {
          // Decrement distance between gains
          npts = 3 + 5*round(min([5,perr*3])) + 17*pind;
          kval = logspace(log10(k[i-1]),log10(k[i]),npts);
          dist = (k[i]-k[i-1])/(npts-7*pind);
          i    = i-1;
          for (kcnt in kval[2:npts])
          {
            i = i + 1;
            k[i] = kcnt;
            ki = kcnt./(1+kcnt*d);
            r1 = eig(a-bc*ki).val';
            r[;i] = vsort(r[;i-1],r1,sp).v2;
            y = imag(r[;i]); y2 = imag(r[;i-1]);
            x = real(r[;i]); x2 = real(r[;i-1]);
            ind = 0; 
          }
       } else {
          // Fix for plot when rlocus goes comes from -inf to inf or -inf to inf
          infind = find(abs(x)>ax && sign(x2)!=sign(x) );
          if (length(infind)>0)
          {
            x[infind] = sign(x2(infind))*ax;
          } 
          // Increase/decrease distance to next k based on linearity estimate:
          dist = dist*(0.3/precision+exp(1-15*perr));
       }
       // Next gain value
       k[i+1] = k[i] + dist;
       // Termination criterion
       terminate=1;
       if (z==p)
       {
          tz = vsort(r1, tz).v2;
          terminate=max(abs(tz-r1)) < ax/100;
       } else { 
          // Make sure all loci tending to infinity are out of graph 
          // before terminating
          mx = max(abs([real(r1).';imag(r1).']));
          // Terminate when all poles have approached their corresponding zeros
          // Note: this may not work well if zeros are close together
          if (z > 0)
          {
             terminate = max(min(abs(r1*ones(1,z)-ones(p,1)*tz.')))<ax/100;
          }
          terminate = ( sum(mx>1.2*ax) >= (p-z) && terminate );
       }
       terminate = terminate || (abs(k[i]) > 1e30);
     }
     return << poles=ep; zeros=tz; rout=(r.'); k=k>>;
  } else {
     // When K is given:
     return k_known (a,b,c,d,ep,tz,sp,k);
  }
};

k_known = function (a,b,c,d,ep,tz,sp,k)
{
     // When K is given:
     ns = b.nr;
     nu = b.nc;
     nk = length(k); 
     i  = 1;
     r  = sqrt(-1) * ones(ns,nk);
     bc = b*c; 
     // Find eigenvalues of:  A - B*inv(I+k*D)*k*C:
     k = k./(1+k.*d);
     r[;1] = eig(a-bc*k[1]).val'; 
     k[1] = [];
     for (kk in k)
     {
        i = i + 1;
        r1 = eig(a-bc*kk).val';
        // Sort eigenvalues
        r[;i] = vsort(r[;i-1],r1,sp).v2;
     }
     return <<poles=ep; zeros=tz; rout=r.'; k=k>>;
};




