single_slope = function(ds){
//
// returns a single monotonic slope of input chunk
// (trims beginning and end of waveform)
//
// STATUS: New (seems top work)
// [azg] Sat Jul 21 17:23:41 PDT 2018
////////////////////


ind[1] = (mini(ds));
ind[3] = (maxi(ds));
ind[4] = length(ds) - mini(ds[length(ds):1:-1]);
ind[2] = length(ds) - maxi(ds[length(ds):1:-1]);
# printf("ind[1] = %i,ind[2] = %i,ind[3] = %i,ind[4] = %i\n",  ind[1],  ind[2],  ind[3],  ind[4]); 
# plot(ds[ind2:ind1]);
i1 = sort(ind).val[2];
i2 = sort(ind).val[3];

ds=ds[i1:i2];

# calculate the slope 5% to 95%
minds = min(ds);
maxds = max(ds);

# This is ambiguous; depends on slope sign
# for positive slope
indmax = [find(ds > 0.95 *(maxds - minds)  + minds)][1];
indmin = max([find(ds < 0.05 *(maxds - minds)  + minds)]);

# for negative slope
indmax_n = max([find(ds > 0.95 *(maxds - minds)  + minds)]);
indmin_n = [find(ds < 0.05 *(maxds - minds)  + minds)][1];

slope = (ds[indmax] - ds[indmin]) / (indmax - indmin);

# correct the range to cover 0-100%
if(slope > 0) { 
	i1 = indmin - int(0.04*(indmax - indmin)) ; 
	i2 = indmax + int(0.04*(indmax - indmin)) ; 
}

# modified on Mon May 27 16:34:05 KST 2019 (0.05 -> 0.04) because of occasional negative indexes

if(slope < 0) { 
	i1 = indmax_n + int(0.04*(indmax_n - indmin_n)) ; 
	i2 = indmin_n - int(0.04*(indmax_n - indmin_n)) ; 
}



return ds[i1:i2];
};
