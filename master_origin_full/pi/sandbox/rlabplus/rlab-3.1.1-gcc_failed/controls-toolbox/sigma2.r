//----------------------------------------------------------------------
//
// sigma2
//
// Syntax: sv = sigma2(a,b,c,d,Type,w)
//
// Singular Value Frequency Response.
//
// SV = sigma2(A, B, C, D, TYPE, W) produces the matrix SV
// containing the singular values of the square system : 
//                .
//                x = Ax + Bu
//                y = Cx + Du
//                                                   -1
//     with the frequency response G(jw) = C(jwI - A)  B + D
//
// SIGMA calculates the SVD of one of the following types:
//
//     Type = 1   ----   G(jw) 
//     Type = 2   ----   inv(G(jw))
//     Type = 3   ----   I + G(jw)
//     Type = 4   ----   I + inv(G(jw)) 
//
// Vector W contains the frequencies at which the frequency response
// is to be evaluated. The SV matrix has rows which correspond to the 
// singular values in descending order.
//
// -------------------------------------------------------------------

require freqrc

sigma2 = function(a,b,c,d,Type,w)
{
  mg  = freqrc(a,b,c,d,w);
  rmg = mg.nr; cmg = mg.nc;
  rb  = b.nr;  cb  = b.nc;
  rc  = c.nr;  cc  = c.nc;
  gg  = ones(rc,cb);
  for (is in 1 : cmg) {
    gg = mg[;is][:];
    if (Type == 1) {
      sv[;is] = svd(gg).sigma[:];
    }
    if (Type == 2) {
      sv[;is] = svd(inv(gg)).sigma[:];
    }
    if (Type == 3) {
      sv[;is] = svd(eye(cb,cb) + gg).sigma[:];
    }
    if (Type == 4) {
      sv[;is] = svd(eye(cb,cb) + inv(gg)).sigma[:];
    }
  }
  return sv;
};


