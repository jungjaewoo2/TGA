extract_chunk = function(Im, M,threshold){
// finds index range of postive or negative response to square wave input
//
// M= chunk number
//
// STATUS: NEW untested [azg] Thu Mar 29 16:00:24 PDT 2018
// BUGS: triggered by all glitches. Needs gentle LPF

###########################################################################
#  extract_chunk( abs(  d[;1]   - sum(d[;1])/length(d[;1]))  ,3)
# Fri Apr 20 14:03:16 PDT 2018
Im = abs( Im -sum(Im)/length(Im));
###########################################################################

// find mid level threshold
#  if(! exist(threshold)){ th = (max(Im)+min(Im))/2; }
if(! exist(threshold)){ th = 50 + (sum(Im))/length(Im); }
// hack
if( exist(threshold)){ th = threshold;}

skip=100; # skip transition ringing/instability
# skip=0; # skip transition ringing/instability


N[1]=1;
for(i in 2:M+1:2){

N[i] = N[i-1]+ skip + find(Im[N[i-1]+skip:length(Im)] > th)[1];
N[i+1] = N[i] + skip + find(Im[N[i]+skip:length(Im)] < th)[1] ;
}

// skip initial and final 10? samples

# N=N+10;
N[M]=N[M] + 10;
N[M+1]=N[M+1] - 10;


# return N; # debug only
return N[M:M+1];
};
