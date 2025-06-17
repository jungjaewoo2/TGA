extract_chunk_with_slope = function(sig,k,dir,delta){

//	dir= +1 or -1 // for rising or falling slope
//	k = segment number
//
// [azg] Sun Jul 22 00:37:43 PDT 2018
// STATUS: somewhat flaky
// [azg] Fri Aug 24 23:52:38 PDT 2018
// Fixed bug causing ingnoring first falling slope
/////////////////////////////////////

require debounce_transitions

# d=diff(sign(sig - mean(sig)));
# d=diff(sign(sig - 1700));
#d=abs(diff(sign(sig - (max(sig) -min(sig))/2)));
# Fri Jan 25 10:33:13 KST 2019
d=abs(diff(sign(sig - (max(sig) +min(sig))/2)));

ind = debounce_transitions(find(d > 0),300);

if(!exist(dir)) {dir=1; }
## if((dir == +1) && (sig[ind[k+1]] < sig[ind[k+1]+1000]))  {k=k; else k=k+1;}
## if(dir == -1){ if((dir == -1) && (sig[ind[k+1]] > sig[ind[k+1]+1000]))  {k=k; else k=k+1;}}
if(dir == +1){ if (sig[ind[k+1]] < sig[ind[k+1]+1000]) {k=k; else k=k+1;}}
if(dir == -1){ if (sig[ind[k+1]] > sig[ind[k+1]+1000]) {k=k; else k=k+1;}}

if(!exist(delta)) {delta = max(diff(ind)); }
start=ind[k];
end  =ind[k+2];

# debug only
printf("Start = %i, End = %i, Delta = %i\n", start, end,delta );
return sig[start:end];
};
