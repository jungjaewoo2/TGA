//-------------------------------------------------------------------
//
// xsi
//
// Syntax: Si=xsi(ns,i)
//
// Computes the skew-symmetric basis matrix corresponding to
// entry i for a ns x ns matrix.
//
// Si = [  0    1
//               jk
//         -1   0   ]
//           jk
// for j<>k
//
// The ith entry is counted across rows from the first upper
// diagonal location.
//
// This function is used in the optimal S=-S^* computation for
// covariance control. (xsoptj)
//
// Originally written by L. D. Peterson
// Modified by Jeff Layton and Ported to RLaB, 10/13/93
//
// Version JBL 931013
//-------------------------------------------------------------------

xsi = function(ns,i)
{
   local(m,j,k,Si,icount)

// Begin with a full matrix of zeros

   Si=zeros(ns,ns);

// Number of skew-symmetric free parameters

   m=ns*(ns-1) / 2;
   if ( i > m) {
       error("i > m in xsi");
   }
   if (i == 0) {
       error("i == 0 in xsi");
   }

// Determine j and k by counting across rows

   j=1;
   k=2;

   for (icount in 1:m) {
        if (icount == i)  {icount=2*m;}
        if (k != ns) {
//
// Case 1: not at end of row
//
            k=k+1;
        else
//
// Case 2: at end of row
//
            j=j+1;
            k=j+1;
        }
   }

   Si[j;k]=1;
   Si[k;j]=-1;

   return Si
};

