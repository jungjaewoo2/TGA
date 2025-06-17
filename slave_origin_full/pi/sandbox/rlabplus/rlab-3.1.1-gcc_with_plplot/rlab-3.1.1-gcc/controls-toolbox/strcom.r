//----------------------------------------------------------------------
//
// strcom
//
// Syntax: b=strcom(a)
//
// This routines takes a string that has been split in to columns using
// strsplt and recomines them into the string b which is returned.
//
// Author: Jeff Layton 6/20/93
// Version JBL 930620
//----------------------------------------------------------------------

strcom = function(a)
{
   local(nc,i,b,bdum)

   nc=a.nc;

   if (nc == 0) {
       error("strcom: Input string is not in row form.");
   }

   b="";
   bdum="";
   for (i in 1:nc) {
        sprintf(bdum,"%s",a[1;i]);
        b=b+bdum;
   }

   return b;
};

