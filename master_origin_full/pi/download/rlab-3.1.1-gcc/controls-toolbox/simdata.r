//-----------------------------------------------------------------------
//
// simdata
//
// Syntax:  Y=simdata(A,B,C,D,n)
//
// This routine forms the ERA data matrix Y from successive
// simulations of length n samples to the system
// 
//  x[k+1] = Ax[k] + Bu[k]
//    y[k] = Cx[k] + Du[k]
// 
// Originally written by Lee D. Peterson
// Modified and ported to RLaB by Jeffrey B. Layton
// Version JBL 940919
//-----------------------------------------------------------------------

rfile abcdchk

simdata = function(A,B,C,D,n)
{
   local(l,m,Y,akm1,ptr,estr,msg,k)

// See if A,B,C, and D are compatible.
   msg="";
   msg=abcdchk(A,B,C,D);
   if (msg != "") {
       estr="SIMDATA: "+msg;
       error(estr);
   }

// Get problem dimensions
   l=C.nr;
   m=B.nc;
   Y=zeros(n*l,m);

//
// Form each Markov parameter and stuff it into Y
//
   Y[1:l;]=D;
   Y[l+1:l+l;]=C*B;
   akm1=B;
   for (k in 2:n-1) {
        akm1=A*akm1;
        ptr=k*l;
        Y[ptr+1:ptr+l;]=C*akm1;
   }
   return Y;
};

