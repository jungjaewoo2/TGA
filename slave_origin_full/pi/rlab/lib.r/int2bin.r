/////////////////////////////////////////////////////////////////
// convert a decimal vector to n-bit binary
// [azg] needs work to 
//       -needs better error checking
//            eg. int2bin(8,3)
//			2             0             0 
//        - convert to string matrix to preserve dim(?)
//
// NOTE: if not specified defaults to miniumum number of bits 
// STATUS: new 
//         [azg] Sat Jun 21 00:23:51 PDT 2008
//         [azg] modified Tue Jul  1 21:24:32 PDT 2008
//
/////////////////////////////////////////////////////////////////
require exp2 log2
int2bin = function(M,n){
  
M=M[:];
col=M[:].nr;
  if ((nargs<2)){  n=max(floor(log2(abs(M)))+1); }
  if ((nargs == 2) && (max(M) >= 2^n)){ error(" Too few bits to represent M");}
  z=zeros(col,n); 
for (i in 1:col){
  for (j in 1:n){  
// new rlabplus requires additional *1.0 multiplication
// or dies with error "cannot assign integer to real matrix and viceversa"
// (is this Goto BLAS related?)
    z[i;j]=1.0*floor(M[i]./exp2(n-j)); 
    M[i]=M[i]-z[i;j].*exp2(n-j); 
    }
}
  
  return z; 
  };
