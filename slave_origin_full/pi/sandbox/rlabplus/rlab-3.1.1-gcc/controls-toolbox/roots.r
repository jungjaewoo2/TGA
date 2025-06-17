//-------------------------------------------------------------------//

//  Synopsis:   Find the roots of a polynomial.

//  Syntax:     roots(P)

//  Description:

//  Find roots of polynomial P of degree n.
//
//  Case 1:
//    P is a vector (row or column) with n+1 components representing
//    polynomial coefficients in descending powers.
// 
//  Case 2:
//    P is a linked-list of n+1 square matrices, the n+1 matrices in P 
//    representing polynomial coefficients in descending powers.
//
//  Return a column vector containing roots of P.

//  Example 1:
//  The following script finds the roots of 
//  x^3 + 3*x^2 + 5 = 0
// 
//  > p=[1,3,0,5]
//   p =
//          1          3          0          5
//  > r=roots(p)
//   r =
//             0.213 - 1.19i
//             0.213 + 1.19i
//                -3.43 + 0i
//  > q=<<1;3;0;5>>;
//  > r=roots(q)
//   r =
//             0.213 - 1.19i
//             0.213 + 1.19i
//                -3.43 + 0i
//  
//  Example 2:
//  Find the damped and undamped frequencies of the following oscillation 
//  system (e.g. mass-damper-spring system, LRC circuit, etc.):
//           2
//          d u     du
//        M --- + C -- + K u = 0
//            2     dt
//          dt  
//
//  where M =   5  0      C =  8  -4      K =  50  -25
//              0  5          -4   8          -25   50
//
//  The characteristic equation corresponding to the governing equation
//  is
//        M*s^2 + C*s + K = 0
//
//  > M=[5,0;0,5]; C=[8,-4;-4,8]; K=[50,-25;-25,50];
//  > // undamped frequencies
//  > uf = roots(<<M;zeros(2,2);K>>)
//         -2.82e-16 + 2.24i
//         -2.82e-16 - 2.24i
//          6.53e-17 + 3.87i
//          6.53e-17 - 3.87i
//  The undamped vibration frequencies are 2.24 and 3.87 rad/s.
//  > // damped frequencies
//  > df = roots(<<M;C;K>>)
//               -0.4 + 2.2i
//               -0.4 - 2.2i
//              -1.2 + 3.68i
//              -1.2 - 3.68i
//  So we know the damped vibration frequencies are 2.2 and 3.68 rad/s.
//
//  Example 3:
//  Same as Example 2, but  C = 40  -10
//                             -10   40
//                            
//  > M=[5,0;0,5]; C=[40,-10;-10,40]; K=[50,-25;-25,50];
//  > f = roots(<<M;C;K>>)
//      -8.16  
//         -5  
//      -1.84  
//         -1  
//  All roots are negative real ===> no oscillation but stable. 

//  See Also: poly

//  Tzong-Shuoh Yang 
// (tsyang@cal.berkeley.edu)  5/7/94  scalar polynomial roots
//                            3/2/95  matrix polynomial roots
//                            5/4/98  fix bug roots([5,0])
//-------------------------------------------------------------------//

require compan

roots = function(P)
{
   local(P);

   if (!exist(P)) { return []; }

   if (class(P) == "num") {
     // scalar polynomial
     if (P.n <= 1) { return []; }
     if (min(size(P)) != 1) {
        error("Argument of roots() must be a column or a row vector.");
     }
     p = P[:]';
     // trim leading zeros
     z = find(p != 0);
     if (z==[]) { return []; }
     if (z[1] > 1) {
        p = p[z[1]:p.n];
     }
     // scalar polynomial
     if (p.n <= 1) { return []; }
     // find trailing zeros
     z  = find(p != 0);
     nz = p.n - z[z.n];
     p  = p[z[1]:z[z.n]];
     if (p.n == 1)
     {
        return zeros(nz,1);
     }
     // find companion matrix
     c = compan(p);
   else if (class(P) == "list") {
     // matrix polynomial
     n = max(strtod(members(P)));
     if (n <= 1) { return []; }
     ns = P.[1].nr;
     if (ns != P.[1].nc) {
        error("Argument list components of roots() must be square matrices.");
     }
     for (i in 2:n) {
        if (P.[i].nr != ns || P.[i].nc != ns) {
        error("Argument list components of roots() must be the same size.");        
        }
     }
     // find leading and trailing zero matrices
     z = [];
     for ( i in 1:n) {
       if (any((P.[i] != 0)[:])) {
          z = [z,i];
       }
     }
     nz = (n - z[z.n])*ns;
     if (z.n == 1)
     {
        if (nz == 0)
        {
           error ("Argument is not a polynomial.");
        }
        tmp = inv(P.[z[1]]); // check singularity
        return zeros(nz*ns,1);
     }
     // find companion matrix c
     n0 = z[1];
     n1 = z[z.n];
     P.[n0] = factor(P.[n0],"g");
     for (i in (n0+1):n1) {
         P.[i] = -backsub(P.[n0], P.[i]);
     }
     n2 = (n1-n0)*ns;
     c = zeros(n2,n2);
     for (i in 1:(n2-ns)) {
        c[i;ns+i] = 1;
     }
     j = (n2-ns+1):n2;
     k = j;
     for(i in (n0+1):n1) {
        c[j;k] = P.[i];
        k = k - ns;
     }
   else {
     error("Wrong argument class in roots().");
   }}}   
   // find roots and stuff with trailing zero(s)
   r = [zeros(nz,1);eig(c).val'];   
   // trim complex part if all real
   if (all(imag(r) == 0)) {
	  r = real(r);
   }
   // if you don't want the roots sorted in ascending order, 
   // comment out the next line
   r = sort(r).val;
   return r;
};

