subsample = function( a, M){
// subsample vector a M times after truncating length of a to be a multiple of M
//
// This is brute force removal of samples. No fitering of interpolation
//
// STATUS: New UNTESTED

a=a[1:length(a) - mod(length(a),M)];
b=reshape(a,M,length(a)/M)[1;];

return b[:] ;
};
