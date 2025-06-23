//----------------------------------------------------------------------
//
// freqrc
//
// syntax: mg = freqrc(a,b,c,d,w)
//
// Generalized complex frequency response calculation.
//
// G = freqrc(A,B,C,D,W) produces a matrix G containing
// the complex frequency response :
//                                  -1
//     G(s) = Y(s)/U(s) = C(sI - A)   B + D
//
// of the linear system described in  state space as : 
//             .
//             x = Ax + Bu
//             y = Cx + Du
//
// when evaluated at the frequencies in complex vector W.
// Returned in G is a matrix where each column corresponds to a 
// frequency point in W, and each row corresponds to a paticular
// U-Y pair. The first ny rows, where ny is the size of the Y vector,
// correspond to the responses from the first input. And so on up to 
// ny * nu where nu is the size of the U vector.
//
//----------------------------------------------------------------------
require clxbode

freqrc = function(a,b,c,d,w)
{
  //
  // Calculating G(jw) :
  //
  cb = b.nc;
  tmp = hess(a);
  p = tmp.p;
  a = tmp.h;
  b = p'*b;
  c = c*p;
  
  mg = clxbode(a,b,c,d,1,w);
  for (iu in 2:cb)
  {
    mg = [mg;clxbode(a,b,c,d,iu,w)];
  } 
  
  return mg;
};

