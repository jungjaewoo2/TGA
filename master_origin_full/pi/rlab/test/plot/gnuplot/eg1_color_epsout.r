//
//
//
gnuwins(3);

color_func = function(k, k0, s, off)
{
  global(const);

  xp = mod(k - k0 - 0.5, 1)-0.5;

 rval = off + cos(const.pi .* xp) .* exp( -0.5 * (xp .^ 2) ./ s .^ 2 );
//  rval = xp;

  return rval;
};

kr = 1/6;
kg = 1/2;
kb = 5/6;

s = 0.1;
off = 0.05;

k=[0:1:1/128]';
yred_k = color_func(k,kr,s,off);
ygreen_k = color_func(k,kg,s,off);
yblue_k = color_func(k,kb,s,off);

a_r = yred_k    ./ (yred_k + ygreen_k + yblue_k);
a_g = ygreen_k  ./ (yred_k + ygreen_k + yblue_k);
a_b = yblue_k   ./ (yred_k + ygreen_k + yblue_k);

gnuwin (1);
gnuformat([...
  "with lines lw 1 lc @red@", ...
  "with lines lw 1 lc @green@", ...
  "with lines lw 1 lc @blue@", ...
[]]);
gnuplot(<<a=[k,yred_k];b=[k,ygreen_k];c=[k,yblue_k]>>);

gnuwin (2);
gnuformat([...
  "with lines lw 1 lc @red@", ...
  "with lines lw 1 lc @green@", ...
  "with lines lw 1 lc @blue@", ...
[]]);
gnuplot(<<a=[k,a_r];b=[k,a_g];c=[k,a_b]>>);

color_code = text(int(a_r * 255),"%02x") + text(int(a_g * 255),"%02x") + text(int(a_b * 255),"%02x");
c =<<>>;
c.data = text([k,1],"%g"," ") + " 0x" + color_code ;

gnuwin (3);
gnuformat ( "using 1:2:3 with points pt 5 ps 2 lc rgb variable $3" );
gnuplot( c, "gnuplot1_color.eps" );
