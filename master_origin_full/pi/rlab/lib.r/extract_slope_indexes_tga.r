extract_slope_indexes_tga = function(sig){

//
// [azg] Sun Jul 22 00:37:43 PDT 2018
/////////////////////////////////////

require debounce_transitions

# d=diff(sign(sig - mean(sig)));
# d=diff(sign(sig - 1700));
#d=abs(diff(sign(sig - (max(sig) -min(sig))/2)));
# Fri Jan 25 10:33:13 KST 2019
# d=abs(diff(sign(sig - (max(sig) +min(sig))/2)));
# Sun Apr 28 06:42:50 KST 2019 adjust threshold to not detect moves to/from shipping position
d=abs(diff(sign(sig - 0.8*(max(sig) +min(sig))/2)));

ind = debounce_transitions(find(d > 0),300);

# debug only
# ind
return ind;
};
