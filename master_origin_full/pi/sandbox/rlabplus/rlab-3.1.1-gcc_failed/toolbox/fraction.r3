//-------------------------------------------------------------------//

// Synopsis:    convert a real number to a nearest proper fraction.

// Syntax:      fraction ( num , maxerr )

// Description:

// Fractional approximation of a real number num.  Return a row vector
// [a, b, c] where a is the integral part of num, b is the numerator,
// and c is the denominator.  The maximum error can be specified in
// maxerr such that abs(num - a - b/c) < maxerr.  If maxerr is not
// specified, then 1/64 wil be used as maxerr.
//
// Examples:
//
// > fraction(0.33)
//         0          1          3
// > fraction(0.33, 0.002)
//         0         21         64
// > fraction(0.33, 0.001)
//         0         26         79
// > fraction(1.33)
//         1          1          3
// > fraction(-1.25)
//        -1         -1          4
//
// Reference: L. J. Plebani, "Common-fraction approximation of real
//               numbers," Dr. Dobb's J., October, 1995.
//
// Tzong-Shuoh Yang  9/8/95   (yang@isec.com)
// (I love SI much much more than English unit system)
//-------------------------------------------------------------------//
fraction = function(num, maxerr)
{
   local(num, maxerr)

   if(!exist(maxerr))
   {
       maxerr = 1/64;
   }
   num = real(num);
   maxerr = abs(real(maxerr));

   if (num > 0) { s = 1; } else { s = -1; num = abs(num);}
   i = int(num);
   num = num - i;
   a = 0; b = 1; c = 1; d = 1;
   while(1)
   {
       h = a + c;
       k = b + d;
       e = h/k - num;
       if (abs(e) < maxerr)
       { 
           break;
       } else { if (e < 0)
       { 
           a = h;
           b = k;
       } else {
           c = h;
           d = k;
       }}
   }
   return [s*i,s*h,k];
};
