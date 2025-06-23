//---------------------------------------------------
//
// [AG] = xgpart(A,nu,nz,nc)
//
// Returns the corresponding "G" partition from A
//
// L.D. Peterson
// version 900305
//
//---------------------------------------------------

xgpart = function(A,nu,nz,nc)
{
   local(AG);

   AG = A[1:nu;nz+1:nz+nc];

   return AG;
};

