//
// file: eg_sf1.r
//
NITER =10000;

erling = function(k , r, t)
{
  a = r.*t;
  rval = a.^k .* r .* exp(-r.*t) / Gamma(k+1);
  return rval;
};

erlingk = function(k, r0, s, t)
{
  global(const);
  a = r0 ./ s;
  b = 1 + 2 .* s .* t;
  rval = 2 * b.^-(1.5+k) .* exp(-0.5 * a) .* s .* (b-1).^k .* Gamma(1.5+k) ...
      / (const.sqrtpi .* Gamma(1+k)) ...
      .* Hypergeometric1F1(1.5+k, 0.5, 0.5*a./b);
  return rval;
};

//
erling0 = function(r0, s, t)
{
  a = r0 ./ s;
  b = 1 + 2 .* s .* t;
  rval = b.^(-2.5) .* (a+b) .* s .* exp(-0.5 * a + a ./ (2 .* b));
  return rval;
};

//
erling1 = function(r0, s, t)
{
  a = r0 ./ s;
  b = 1 + 2 .* s .* t;
  rval = 0.5 * (b-1) .* b.^(-4.5) .* (a.^2 + 6 * a .* b + 3 * b.^2) ...
      .* s .* exp(-0.5 * a + a ./ (2 .* b));
  return rval;
};

//
erling2 = function(r0, s, t)
{
  a = r0 ./ s;
  b = 1 + 2 .* s .* t;
  rval = 0.125 * (b-1).^2 .* b.^(-6.5) ...
      .* (a.^3 + 15 * a.^2 .* b + 45 * a .* b.^2 + 15 * b.^3) ...
      .* s .* exp(-0.5 * a + a ./ (2 .* b));
  return rval;
};

//
erling3 = function(r0, s, t)
{
  a = r0 ./ s;
  b = 1 + 2 .* s .* t;
  rval = 1/48 * (b-1).^3 .* b.^(-8.5) ...
      .* (a.^4 + 28 * a.^3 .* b + 210 * a.^2 .* b.^2 + 420 * a .* b.^3 + 105 * b.^4) ...
      .* s .* exp(-0.5 * a + a ./ (2 .* b));
  return rval;
};

//
erling4 = function(r0, s, t)
{
  a = r0 ./ s;
  b = 1 + 2 .* s .* t;
  rval = 1/384 * (b-1).^4 .* b.^(-10.5) ...
      .* (a.^5 + 45 * a.^4 .* b + 630 * a.^3 .* b.^2 + 3150 * a.^2 .* b.^3 + 4725 * a .*b.^4 ...
         + 945 * b.^5) ...
      .* s .* exp(-0.5 * a + a ./ (2 .* b));
  return rval;
};

//
erling5 = function(r0, s, t)
{
  a = r0 ./ s;
  b = 1 + 2 .* s .* t;
  rval = 1/3840 * (b-1).^5 .* b.^(-12.5) ...
      .* (a.^6 + 66 * a.^5 .* b + 1485 * a.^4 .* b.^2 + 13860 * a.^3 .* b.^3 + 51975 * a.^2 .*b.^4 ...
      + 62370 * a .* b.^5 + 10395 * b.^6) ...
      .* s .* exp(-0.5 * a + a ./ (2 .* b));
  return rval;
};

//
erling6 = function(r0, s, t)
{
  a = r0 ./ s;
  b = 1 + 2 .* s .* t;
  rval = 1/46080 * (b-1).^6 .* b.^(-14.5) ...
      .* (a.^7 + 91 * a.^6 .* b + 3003 * a.^5 .* b.^2 + 45045 * a.^4 .* b.^3 ...
          + 315315 * a.^3 .*b.^4 + 945945 * a.^2 .* b.^5 + 945945 * a .* b.^6 ...
          + 135135 * b.^7) ...
      .* s .* exp(-0.5 * a + a ./ (2 .* b));
  return rval;
};

//
erling7 = function(r0, s, t)
{
  a = r0 ./ s;
  b = 1 + 2 .* s .* t;
  rval = 1/645120 * (b-1).^7 .* b.^(-16.5) ...
      .* (a.^8 + 120 * a.^7 .* b + 5640 * a.^6 .* b.^2 + 120120 * a.^5 .* b.^3 ...
          + 1353150 * a.^4 .*b.^4 + 7567560 * a.^3 .* b.^5 + 18918900 * a.^2 .* b.^6 ...
          + 16216200 * a .* b.^7 + 2027025 * b.^8) ...
      .* s .* exp(-0.5 * a + a ./ (2 .* b));
  return rval;
};

//
erling8 = function(r0, s, t)
{
  a = r0 ./ s;
  b = 1 + 2 .* s .* t;
  rval = 1/10321920 * (b-1).^8 .* b.^(-18.5) ...
      .* (a.^9 + 153 * a.^8 .* b + 9180 * a.^7 .* b.^2 + 278460 * a.^6 .* b.^3 ...
      + 4594590 * a.^5 .*b.^4 + 41351310 * a.^4 .* b.^5 + 192972780 * a.^3 .* b.^6 ...
      + 413513100 * a.^2 .* b.^7 + 310134825 * a .* b.^8 + 34459425 * b.^9) ...
      .* s .* exp(-0.5 * a + a ./ (2 .* b));
  return rval;
};

//
erling9 = function(r0, s, t)
{
  a = r0 ./ s;
  b = 1 + 2 .* s .* t;
  rval = 1/185794560 * (b-1).^9 .* b.^(-20.5) ...
      .* (a.^10 + 190 * a.^9 .* b + 14535 * a.^8 .* b.^2 + 581400 * a.^7 .* b.^3 ...
      + 13226850 * a.^6 .*b.^4 + 174594420 * a.^5 .* b.^5 + 1309458150 * a.^4 .* b.^6 ...
      + 5237832600 * a.^3 .* b.^7 + 9820936125 * a.^2 .* b.^8 + 6547290750 * a .* b.^9 ...
      + 654729075 * b.^10) ...
      .* s .* exp(-0.5 * a + a ./ (2 .* b));
  return rval;
};






r = 1;
s = 0.1;
k = 15;
t = [0:20:1/8]';

y1 = erling (k, r, t);

tic();
for (i in 1:NITER)
{ y2 = erling5(r, s, t); }
printf("evaluation of erling%g lasted %g sec\n", k, toc());

tic();
for (i in 1:NITER)
{ y3 = erlingk(k, r, s, t); }
printf("evaluation of erlingk for k=%g lasted %g sec\n", k, toc());

tic();
for (i in 1:NITER)
{ y4 = ErlingN(k, r, s, t); }
printf("evaluation of ErlingN for k=%g lasted %g sec\n", k, toc());

gnulegend (["Erling function for k="+text(k), ...
            "erling"+text(k)+" function with s="+text(s), ...
            "erlingk function for k="+text(k)+" with s="+text(s), ...
            "ErlingN function for k="+text(k)+" with s="+text(s) ]);
gnuformat( ["with lines", "with lines", "with points", "with points"] );
gnuplot( [t,y1,y2,y3,y4] );
//gnuplot( <<a=[t,y1]; b=[t,y2]>> );


