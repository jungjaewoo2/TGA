//-------------------------------------------------------------------------------
//
// dsort
//
// Syntax: </ndx;s/> = dsort(p)
//
// Sort complex discrete eigenvalues in descending order.
//
// dsort(P) sorts the complex eigenvalues in the vector P in
// descending order by magnitude.  The unstable eigenvalues will
// appear first.
//
// </ndx;s/> = dsort(P) also returns the vector ndx containing the
// indexes used in the sort.
//
//-------------------------------------------------------------------------------

dsort = function(p)
{

   if (!exist(p)) {error("p doesn't exist");}

   if (p.nr == 1) {
       p=p.';
   }

   a=sort(-abs(p));
   s=a.val;
   ndx=a.idx';

   for (i in 1:p.nc) {
        k=1;
        while( k < length(s)) {
              if (imag(s[k;i]) != 0) {
                  if (imag(s[k;i]) < 0) {
                      s[k:k+1;i]=conj(s[k:k+1;i]);
                      swap=ndx[k;i];
                      ndx[k;i]=ndx[k+1;i];
                      ndx[k+1;i]=swap;
                  }
                  k=k+2;
              } else {
                  k=k+1;
              }
        }
   }

   return << s=s; ndx=ndx >>;
};

