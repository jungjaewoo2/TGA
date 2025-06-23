//--------------------------------------------------------------------------------
//
// ddamp
//
// Syntax:  A = ddamp(a,Ts)
//          </mag;wn;z/> = ddamp(a,Ts)
//
// This routine compures the natural frequency and damping for discrete
// systems. It prints a table of the eiganvalues of the matrix a, the
// z-plane magnitude, the equivalent s-plane natural frequency (rad/s),
// the equivalent s-plane frequency (Hz.), and the damping factors.
//
// The routine can is called as,
//
//     A=ddamp(a,Ts)
//
// where a is the input variable and Ts is the sample time. The variable a
// can be in several forms:
//
//       1) If a is square, it is assumed to be the state-space
//          "a" matrix.
//       2) If a is a row vector, it is assumed to be a vector of
//          the polynomial coefficients from a transfer function.
//       3) If a is a column vector, it is assumed to contain
//          root locations.
//
// It returns the s-plane equivalent natural frequency (rad/s), the damping
// factors, and the z-plane magnitude. If the routine is called with only the
// input variable a (the sample time Ts is not input), then ddamp will not
// print the table and will return values for only the magnitude.
//
// For a discrete system eigenvalue, lambda, the equivalent s-plane
// natural frequency and damping ratio are
//
//           Wn = abs(log(lambda))/Ts
//
//            Z = -cos(angle(log(lambda)))
//
//
// The routine returns the natrural frequencies, damping, and the z-plane
// magnitudes in a list:
//
//  A.wn  = Natural Frequencies (rad/s)
//  A.z   = Damping Factors
//  A.mag = z-plane magnitudes
//
// Copyright (C), by Jeffrey B. Layton, 1994
// Version JBL 940917
//--------------------------------------------------------------------------------

require dsort roots

ddamp = function(a,Ts)
{
   global(eps,pi)
   
   // Start by sorting
   if (a.nr == a.nc) {
       r=dsort(eig(a).val[:]).s;
   } else { if (a.nr == 1) {
       r=dsort(roots(a)).s;
   } else { if (a.nc == 1) {
       r=a;
   } else {
       error("Must be a vector or a Square matrix.");
   }}}

   mag=abs(r);
   
   if (nargs == 2) {   
       // If sample time is given solve for equivalent s-plane roots
       s=log(r)/Ts;
       wn=abs(s);
       z=-cos(atan2(imag(s),real(s)));
       wnhz=sqrt(s)/(2.0*pi);

       printf("\n");
   printf("     Eigenvalue              Mag       Damping     Freq. (rad/s)    Freq. (Hz). \n");
   printf("     ----------              ---       -------     ------------     ----------- \n");
       for (i in 1:wn.nr) {
            if (imag(s[i]) < 0.0) {
                printf("%-10.6f -j%-10.6f", real (s[i]), abs(imag (s[i])));
            } else {
                printf("%-10.6f +j%-10.6f", real (s[i]), imag (s[i]));
            }
            printf("   ");
            printf("%-10.6f",mag[i]);
            printf("   ");
            printf("%-10.6f",z[i]);
            printf("     ");
            printf("%-10.6f",wn[i]);
            printf("      ");
            printf("%-10.6f",wnhz[i]);
            printf(" \n");
       }
       printf(" \n");
   } else {
       s=[];
       wn=[];
       z=[];
   }
   
   return << wn=wn; z=z; mag=mag >>;
};

