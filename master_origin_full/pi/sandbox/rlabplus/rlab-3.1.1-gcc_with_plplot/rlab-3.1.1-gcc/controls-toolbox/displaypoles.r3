//-----------------------------------------------------------------------
// displaypoles
//
// Syntax:
//          poles=displaypoles(A,description)
//
// This routine computes the poles of A and their corresponding damping
// ratio and damped circular frequency to the screen. The input, description,
// can be used as a heading for the output.
//
// The routine returns the poles of A.
// function poles = displaypoles(A,description,batch)
// 
// Originally Written by L. D. Peterson
// Modified and ported to RLaB by J. B. Layton
// Version JBL 940917
//-----------------------------------------------------------------------

require banner  isstr

displaypoles = function(A,description)
{

   if (nargs == 1) {
       iopt=0;
       str1="";
   } else { if (nargs == 2) {
       if (isstr(description)) {
           iopt=0;
       } else {
           iopt=description;
           str1="";
       }
   }}

   Adum=eig(A);              // poles of the design system
   poles=Adum.val;
   poles=sort(poles).val';
   cmplxpoles=find(abs(imag(poles))>1.e-10);
   realpoles=find(abs(imag(poles))<=1.e-10);
   poles_zeta=zeros(size(poles));
   poles_freq=zeros(size(poles));

   if (max(size(cmplxpoles)) != 0) {
       poles_zeta[cmplxpoles] = ...
   sin(atan(-real(poles[cmplxpoles])./abs(imag(poles[cmplxpoles]))));
       poles_freq[cmplxpoles] = ...
   -real(poles[cmplxpoles])./poles_zeta[cmplxpoles]./(2*pi);
   }

   if (max(size(realpoles)) != 0) {
       poles_zeta[realpoles]=ones(max(size(realpoles)),1);
       poles_freq[realpoles]=-real(poles[realpoles])./(2*pi);
   }
   for (j in 1:max(size(poles))) {
       if ((poles_zeta[j] != 0.) && (poles_freq[j] != 0.)) {
            poles_tc[j]=((poles_zeta[j]*poles_freq[j]*2*pi)^(-1));
       } else {
            poles_tc[j] = 0.;
       }
   }

   if ( (nargs >= 2) && (length(str1) != 0) ) {
        printf("/n");
        banner(str1+" Poles");
        printf("/n");
   }

   outstr="     re(rad/sec)   im(rad.sec)    zeta(//)";
   printf("%s",outstr);
   outstr="     f(hz)    time constant (sec) \n";
   printf("%s",outstr);
   outstr="---------------------------------------------------";
   printf("%s",outstr);
   outstr="-----------------------\n";
   printf("%s",outstr);
   line_count = 0;
   for (j in 1:max(size(poles)) ) {
        line_count=line_count+1;
        sprintf(outstring1,"%3.0f  %12.4e %12.4e ",...
                j,real(poles[j]),imag(poles[j]));
        if (poles_tc[j] != 0) {
            sprintf(outstring2,"  %6.2f   %12.4e %12.4e",...
                    poles_zeta[j]*100,poles_freq[j],poles_tc[j]);
        } else {
            sprintf(outstring2,"  %6.2f   %12.4e      --inf--",...
                    poles_zeta[j]*100,poles_freq[j]);
        }
        outstr=outstring1+outstring2+" \n";
        printf("%s",outstr);
   }

   return poles;
};

