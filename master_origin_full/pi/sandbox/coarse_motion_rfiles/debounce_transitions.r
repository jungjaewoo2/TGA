debounce_transitions = function (ind, debounce_range){

// remove all indexes of transitions within debounce_range except the first one
if(!exist(debounce_range)) { debounce_range=100;}
# read("junk_ind.txt");
# who();
# size(d_ind)
           # 1            60  
ind=ind[:];
# ind=d_ind[:];
k=1; while(k < length(ind)) { debounce= find(ind < (ind[k] +debounce_range)); ind=rmrows(ind, debounce[k+1:length(debounce)]); k++;}

return ind;
};
