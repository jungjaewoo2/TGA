//----------------------------------------------------------------------
//
// c2dm
//
// syntax: </ad;bd;cd;dd/> = c2dm(a,b,c,d,Ts,"method",w)
//         </den;num/>     = c2dm(NUM,DEN,Ts,"method")
//
// Conversion of continuous LTI systems to discrete-time.
// </Ad;Bd;Cd;Dd/> = c2dm(A,B,C,D,Ts,"method") converts the continuous-
// time state-space system (A,B,C,D) to discrete time using "method":
//  "zoh"       Convert to discrete time assuming a zero order
//              hold on the inputs.
//  "foh"       Convert to discrete time assuming a first order 
//              hold on the inputs.
//  "tustin"    Convert to discrete time using the bilinear 
//              (Tustin) approximation to the derivative.
//  "prewarp"   Convert to discrete time using the bilinear 
//              (Tustin) approximation with frequency prewarping.
//              Specify the critical frequency with an additional
//              argument, i.e. c2dm(A,B,C,D,Ts,"prewarp",Wc)
//  "matched"   Convert the SISO system to discrete time using the
//              matched pole-zero method.
//
// </DENd;NUMd/> = c2dm(NUM,DEN,Ts,"method") converts the continuous-
// time polynomial transfer function G(s) = NUM(s)/DEN(s) to discrete
// time, G(z) = NUMd(z)/DENd(z), using "method".
//
// See also: c2d, and d2cm.
//
//----------------------------------------------------------------------
require abcdchk c2d nargchk roots
require series ss2tf ss2zp tfchk tf2ss zp2ss

c2dm = function(a,b,c,d,Ts,method,w)
{
  local(a,b,c,d,Ts,method,w)
  
  msg = nargchk(3,7,nargs);
  if (msg != "") { error (msg); }
  
  tf = 0;
  // determine which syntax is being used
  if (nargs==3) {
    // Transfer function without method, assume "zoh"
    </den;num/> = tfchk(a,b);
    Ts = c;
    method = "zoh";
    </a;b;c;d/> = tf2ss(num,den);
    tf = 1;

  else if (nargs==4){
    // Transfer function with method.
    </den;num/> = tfchk(a,b);
    Ts = c;
    method = d;
    </a;b;c;d/> = tf2ss(num,den);
    tf = 1;

  else if (nargs==5) {
    if (isstring(d)) {
      // Transfer function with method and prewarp const.
      </den;num/> = tfchk(a,b);
      w = Ts;
      Ts = c;
      method = d;
      </a;b;c;d/> = tf2ss(num,den);
      tf = 1;
    else			
      // State space system without method, assume "zoh"
      msg = abcdchk(a,b,c,d);
      if (msg != "") { error(msg); }
      method = "zoh";
    }

  else			
    // State space system with method.
    msg = abcdchk(a,b,c,d);
    if (msg != "") { error(msg); }
  }}}
  
  nx = a.nr; na = a.nc;
  nb = b.nr; nu = b.nc;

  // Determine conversion method
  if (method=="zoh") {
    // Zero order hold approximation.
    </bd;ad/> = c2d(a,b,Ts);
    cd = c; dd = d;

  else if (method=="foh") {
    // First order hold approximation.
    // Add nu integrators in series
    </ac;bc;cc;dc/> = series(zeros(nu,nu),eye(nu,nu),eye(nu,nu),zeros(nu,nu),a,b,c,d);
    </bd;ad/> = c2d(ac,bc,Ts);
    // Add (z-1)/(z*Ts) in series. This is equivalent to differentiating u.
    </ad;bd;cd;dd/> = series(zeros(nu,nu),eye(nu,nu),-eye(nu,nu)./Ts,eye(nu,nu)./Ts,ad,bd,cc,dc);

  else if (method=="tustin") {
    // Tustin approximation.
    I = eye(nx,nx);
    P = inv(I - a.*Ts/2);
    ad = (I + a.*Ts/2)*P;
    bd = P*b;
    cd = Ts*c*P;
    dd = cd*b/2 + d;

  else if (method=="prewarp") {
    // Tustin approximation with frequency prewarping.
    if (!((nargs==5)||(nargs==7))) {
      error("The critical frequency must be specified when using 'prewarp'.");
    }
    T = 2*tan(w*Ts/2)/w; // Prewarp
    I = eye(nx,nx);
    P = inv(I - a.*T/2);
    ad = (I + a.*T/2)*P;
    bd = P*b;
    cd = T*c*P;
    dd = cd*b/2 + d;
  
  else if (method=="matched") {
    // Matched pole-zero approximation.
    ny = d.nr; nu = d.nc;
    if (ny>1||nu>1) {
      error("System must be SISO for matched pole-zero method.");
    }
    if (tf && ny && nu) {
      z = roots(num); 
      p = roots(den);
      kc = num[length(num)]/den[length(den)];
    else
      </tmp;p;z/> = ss2zp(a,b,c,d);
      kc = -c/a*b + d;
    }
    pd = exp(p*Ts);
    zd = zeros(length(z),1);
    if (!isempty(z)) 
    {
      zd[z!=inf()] = exp(z[z!=inf()]*Ts);
      zd[z==inf()] = -1*ones(length(z[z==inf()]),1);
    }
    ndx = find(z==inf());
    if (!isempty(ndx)) { zd[ndx[1]]=inf(); }// Put one infinite zero at infinity.
    </ad;bd;cd;dd/> = zp2ss(zd,pd,1);
    // Match D.C. gain
    kc = -c/a*b + d;
    kd = cd/(eye(nx,nx)-ad)*bd + dd;
    km = sqrt(abs(kc/kd));
    sm = sign(kc/kd);
    bd = bd.*km;
    cd = cd.*km.*sm;
    dd = dd.*km.*km.*sm;

  else
    error("Conversion method is unknown.");

  }}}}}

  if (tf) {
    // Convert to TF form for output
    return ss2tf(ad,bd,cd,dd,1);
  }

  return << ad=ad; bd=bd; cd=cd; dd=dd>>;
};
