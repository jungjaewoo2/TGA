// -------------------------------------------------------------------------
//
// clxbode
//
// Syntax: g = clxbode(a,b,c,d,iu,w);
//
//     produces a complex gain instead of an absolute gain.
//
// -------------------------------------------------------------------------
require ltifr

clxbode = function(a,b,c,d,iu,w)
{
  cs = max(size(w));
  b = b[;iu];
  d = d[;iu];
  s = sqrt(-1) * w;
  g = ltifr(a,b,s);
  g = c * g + diag(d) * ones(c.nr,cs);

  return g;
};


