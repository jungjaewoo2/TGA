//----------------------------------------------------------------------
// drmodel
//
// syntax:  </a;b;c;d/> = drmodel(n,p,m,option)
//
// Generates random stable discrete nth order test models.
//
// </den;num/> = drmodel(N,,,"tf") forms an Nth order SISO transfer 
//               function model.
//
// </den;num/> = drmodel(N,P,,"tf") generates an Nth order single input, 
//               P output transfer function model.
//
// </A;B;C;D/> = drmodel(N,,,"ss") generates an Nth order SISO state 
//               space model
//
// </A;B;C;D/> = drmodel(N,P,M,"ss") generates an Nth order, P output, M 
//               input, state space model.
//		
// See also: rmodel
//
//----------------------------------------------------------------------
require orth zp2tf

drmodel=function(n,p,m,option)
{
  global(eps,pi)
  
  if (!exist(option)) {
     error("Must specify output option: tf or ss.");
  else
     if (option != "ss" && option != "tf") {
         error("Output option is either tf or ss.");
     }
  }
  if (nargs < 1 || nargs > 4) {
     error("Wrong number of input arguments.");
  }

  rand("normal",0,1);
  if (!exist(n)) {
    n=max([1,round(abs(10*rand()))]);
    if (option=="ss") {
      p=max([1,round(4*rand())]);
      m=max([1,round(4*rand())]);
    else
      m=1;
      p=1;
    }
  }
  if (!exist(p)) { p=1; }
  if (!exist(m)) { m=1; }

  rand("uniform",0,1);
  // Prob of an integrator is 0.10 for the first and 0.01 for all others
  nint = (rand(1,1)<0.10) + sum(rand(n-1,1)<0.01);
  // Prob of repeated roots is 0.05
  nrepeated = floor(sum(rand(n-nint,1)<0.05)/2);
  // Prob of complex roots is 0.5
  ncomplex = floor(sum(rand(n-nint-2*nrepeated,1)<0.5)/2);
  nreal = n-nint-2*nrepeated-2*ncomplex;

  // Determine random poles
  rep = 2*rand(nrepeated,1)-1;
  mag = rand(ncomplex,1);
  ang = pi*rand(ncomplex,1);
  jay = sqrt(-1);
  complex = mag.*exp(jay*ang);
  re = real(complex); 
  im = imag(complex);

  if (option=="tf") {
     // Random transfer function
     poles = [re+jay*im;re-jay*im;ones(nint,1);rep;rep;2*rand(nreal,1)-1];
  
     // Determine random zeros
     zer = inf()*ones(n,p);
     jay = sqrt(-1);
     for (i in 1:p) {
       rand("uniform",0,1);
       // Prob of complex zeros is 0.35
       ncomplex = floor(sum(rand(n,1)<0.35)/2);
       // Prob of real zero is 0.35
       nreal = sum(rand(n-2*ncomplex,1)<0.35);
       nzeros = 2*ncomplex+nreal;
       if (nzeros>0) {
         rand("normal",0,1);
         re = rand(ncomplex,1); 
         im = rand(ncomplex,1);
         zer[1:nzeros;i] = [re+jay*im;re-jay*im;rand(nreal,1)];
       }
     }
     return zp2tf(zer,poles,rand(p,1));
  else 
     // Random state space model
     a = zeros(n,n);
     for (i in 1:ncomplex) {
       ndx = [2*i-1,2*i];
       a[ndx;ndx] = [re[i],im[i];-im[i],re[i]];
     }
     ndx = [2*ncomplex+1:n];
     if (!isempty(ndx)) {
       a[ndx;ndx] = diag([ones(nint,1);rep;rep;2*rand(nreal,1)-1]);
     }
     T = orth(rand(n,n));
     a = T\a*T;
     b = rand(n,m);
     c = rand(p,n);
     d = rand(p,m);
     rand("uniform",0,1);
     b = b .* (rand(n,m)<0.75);
     c = c .* (rand(p,n)<0.75);
     d = d .* (rand(p,m)<0.5);
     return <<a=a;b=b;c=c;d=d>>;
  }

};

