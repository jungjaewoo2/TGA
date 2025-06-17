//--------------------------------------------------------------------
//
// poly2str
//
// Syntax: a = poly2str(plynom,tvar)
//
// This routine returns a string (as part of a list) that is the
// representation of the polynomial coefficients contained in the
// variable plynom multipliedulitplied by the powers of the transform
// variables in tvar (usually either "s" or "z").
//
// This routine returns the following list:
//
//       a.sout = string
//       a.len = string length
//
// Copyright (C), by Jeffrey B. Layton, 1994
// Version JBL 940513
//--------------------------------------------------------------------

poly2str = function(plynom,tvar)
{
   local(polylen,lenout,first_flag,old,term,sout,j)

   if (!exist(tvar)) { tvar = "s"; }
// Determine length of polynomial
   polylen=length(plynom);

// Initialize string
   sout="";

// Initialize junk variables
   lenout=0;            // Length of output string
   first_flag=1;     // First term in polynomial flag
   old=0;            // flag

// Start looping for number of terms in polynomial
   for (j in 1:polylen) {

// Extract current term in polynomial
        term=plynom[j];

// If first term and term != 0 then add to output string
        if ((first_flag == 1) && (term != 0)) {
             if ((term == 1) && (j != polylen)) {
                  sout="  ";
             else
                  sout="   "+num2str(term);
             }
             first_flag=0;

// Else, add sign and current term to output string
        else
             if ((term == 1) && (j != polylen)) {
                  sout=sout+" + ";
             else if (term == 0) {
// In current version it does nothing if coeffiecient of
// term is zero.
             else if (term >= 0) {
                  sout=sout+" + "+num2str(term);
             else
                  sout=sout+" - "+num2str(abs(term));
             }}}
        }
// If term != 0 then add transform variable and sign for coefficient.
        if (term != 0) {
            if ((polylen-j) > 1) {
                 sout=sout+" "+tvar+"^"+num2str(polylen-j);
            else if ((polylen-j) == 1){
                 sout=sout+" "+tvar;
            }}
        }
        if (((length(sout)-old)>63) && (j != polylen)) {
              lenout=max([lenout,length(sout)-old])
              sout=sout;
              old=length(sout)-2;
        }
   }

   if (isempty(sout)) {
       sout="   0";
   }
   lenout=max([lenout,length(sout)-old]);

   return << sout=sout; lenout=lenout >>;
};

