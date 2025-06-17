//---------------------------------------------------------------------------
//
// damp
//
// Syntax: s=damp(a)
//
//
// This routine computes the natural frequency and damping for
// continuous systems. It prints a table of the eigenvalues of
// the matrix a, the associated damping factors, the associated natural
// frequency (rad/s and Hz.).
//
// Also, the routine returns the natrural frequencies and damping
// factors in a list L:
//
//  L.omegan = Natural Frequencies (rad/s)
//  L.zeta   = Damping Factors
//
// Copyright (C), by Jeffrey B. Layton, 1994
// Version JBL 940918
//---------------------------------------------------------------------------

require roots esort

damp = function(a)
{
   global(pi)
   
   // Compute the eigenvalues of the matrix and sort them
   if (a.nr == a.nc) {
       eigvlsort = esort(eig(a).val[:]).s;
   } else { if (a.nr == 1) {
       eigvlsort = esort(roots(a)).s;
   } else { if (a.nc == 1) {
       eigvlsort = a;
   } else {
       error("damp: Must be a vector or a sqaure matrix.");
   }}}

   // Compute the natural frequency on rad/s.
   omegan=abs(eigvlsort);

   // Compute the damping factor
   zeta=-cos(atan2(imag(eigvlsort),real(eigvlsort)));

   // Compute the natural frequency (Hz).
   omegahz=sqrt(eigvlsort)/(2.0*pi);

   printf("\n");
   printf("     Eigenvalue            Damping     Freq. (rad/s)    Freq. (Hz). \n");
   printf("     ----------            -------     -------------    -----------\n");
   for (i in 1:omegan.nr) {
        if (imag(eigvlsort[i]) < 0.0) {
            printf("%-10.6f -%10.6fj", real (eigvlsort[i]), abs(imag (eigvlsort[i])));
        } else {
            printf("%-10.6f +%10.6fj", real (eigvlsort[i]), abs(imag (eigvlsort[i])));
        }
        printf("   ");
        printf("%-10.6f",zeta[i]);
        printf("     ");
        printf("%-10.6f",omegan[i]);
        printf("      ");
        printf("%-10.6f",omegahz[i]);
        printf("\n");
   }
   printf("\n");
   
   return << omegan=omegan; zeta=zeta >>;
};

