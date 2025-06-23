//-------------------------------------------------------
//
// [AH] = xhpart(A,nu,nz,nc)
//
// Returns the corresponding "H" partition from A
//
// L.D. Peterson
// version 900305
//
//-------------------------------------------------------

xhpart = function(A,nu,nz,nc)
{
   local(AH)

   AH = A[1:nu;1:nz];

   return AH;
};

