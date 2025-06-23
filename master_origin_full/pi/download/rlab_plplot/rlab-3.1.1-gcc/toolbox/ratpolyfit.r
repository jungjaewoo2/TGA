ratpolyfit = function(x,y,kn,kd)
{
// Ratpolyfit Rational Polynomial Fitting
//
//  This programs finds 2 polynomials N(x) and D(x),
//  of user given order kn and kd respectively,
//  such that N(xi)/D(xi) ~= y(xi) in a least squares sense.
//
//usage: [N, D] = ratpolyfit(x, y, kn, kd)
//
//note: kn and kd must be large enough to get a good fit.
//      usually, kn=kd gives good results
//
//note: If you "overfit" the data, then you will usually have pole-zero
//      cancellations and/or poles and zeros with a very large magnitude.
//      If that happens, then reduce the values of kn and/or kd
//
//note: Often, if you have a good fit, you will find that your polynomials
//      have roots where the real function has zeros and poles.
//
//note: Polynomial curve fitting becomes ill conditioned
//      as the range of x increases and as kn and kd increase
//
//note: If you think that your function goes to infinity at some x value
//      then make sure y(xi) is set equal to Inf at that point.
//      The program will compensate for all +/- Inf values
//
//
//tested under version 6.5 (R13)
//
//see also: polyfit, padefit, polyval, vander
//

//Paul Godfrey
//pgodfrey@conexant.com
//May 31, 2006
//
//         Translated to rlab by Andrew Grzegorek
//         See examples in ratpolyfit_ex1.r

  //default length if none given
  if (exist(kn))
  {
    kn=round(real(kn));
    if (kn<0)
    { kn=0;}
  else
    kn=5;
  }
  if (exist(kd))
  {
    kd=round(real(kd));
    if (kd<0)
    { kd=0; }
  else
    kd=5;
  }

  x=x[:];
  y=y[:];

  // we must remove +/- Inf values from y first
  // and insert those as separate pol.s at the end
  p=find(!finite(y)); // find Nan 0/0 and Inf a/0 values
  // or find abs(y) values > than some really big number

  dinf=[];
  while (length(p)>0)
  {
    y=y.*(x-x[p[1]]); //adjust remaining y values
    y=rmrows(y,p[1]); // remove bad y value, now a NaN
    dinf=[dinf; x[p[1]]]; // remember where pole was
    x=rmrows(x,p[1]); // now remove that x value too
    if (kd>0)
    { kd=kd-1; }
    p=find(!finite(y)); // have all Inf values been removed yet?
  }

  yy=length(y);
  an=ones(yy,kn+1);
  if (kn>0)
  {
    for (k in kn:1:-1)
    {
      an[;k]=x.*an[;k+1]; // form vandermonde matrix
    }
  }
  ad=ones(yy,kd+1);
  if (kd>0)
  {
    for (k in kd:1:-1)
    {
      ad[;k]=x.*ad[;k+1];
    }
  }
  for (k in 1:yy)
  {
    ad[k;]=y[k]*ad[k;];
  }

  // A is basically N-y*D
  A=[an,-ad]; // LS solution is in the null space of A

  </s;u;vt/>=svd(A,"A");  //null space is in the cols of V
  v=vt';
  ND=v[;v.nr]; // use the "most null" vector

  N=ND[1:kn+1].';
  D=ND[(kn+2):ND.nr].';

  if(length(D) <1 )
  { D[1] = 0;}

  D1=D[1]; // kludge

  if (D1==0)
  { D1=1; }

  // we have to make the D polynomial monic since thats
  // what poly makes, so we have to first adjust N
  N = N ./ D1;
  D = D ./ D1;

  // and then add the removed +/- Inf poles back in
  if (length(D)>1)
  {
    if (length(dinf)>0)
    { D=real(poly([dinf; polyroots(D).roots]))'; }
  else
    // D is zero degree polynomial - it is unity
    if (length(dinf)>0)
    { D=real(poly(dinf))'; }
  }

  Dmax=max(abs(D));
  if (Dmax==0)
  { Dmax=1; }

  // normalize max Den value to be +/- 1
  N=N/Dmax;
  D=D/Dmax;

  rval = <<>>;
  rval.num = N;
  rval.den = D;
  return rval;

};
