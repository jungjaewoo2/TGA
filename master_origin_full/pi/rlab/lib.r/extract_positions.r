extract_positions = function (fn) {
//
// fn -> filename or matrix with coarse laser data (from analog laser output)
// returns coarse motion locations: STEP0(HOME), STEP1, STEP2, STEP3, STEP4, STEP5(HOME), STEP6(END)
//
// Usage:  extract_positions("40073226.datcoarse.dat")
//         extract_positions(d)
// 
// Example:
//         for(i in ls("MAT16*.dat")){pltitle(i); err=extract_positions(i); pause();}
//
// 
////////////////////////////////////////////

DEBUG = 0;

require subsample runavg
global(d)

# Q: Why is this broken? A:(needs global(d))
if(class(fn)== "string"){if(isfile(fn)){ xx=read_ascii(fn);}}
# # flush the first row
# d=d[2:d.nr;];

	
if(class(fn)== "num"){d=fn;}

dsmooth = (runavg(d,100));
dsub=subsample(dsmooth,100);
signflat=sign( (abs(diff(dsub))) -11/2);

# [azg] needed? # signflat=signflat[2:length(signflat)];


 #plot(<<abs(max(dsmooth)*diff(signflat)); dsmooth>>);
 #pause();
flatsidx=find(abs(diff(signflat))> 0.5);

# [azg] # idx2 = find( diff(flatsidx)> 25 && diff(flatsidx)< 33);
idx2 = find( diff(flatsidx)> 25 && diff(flatsidx)< 40);
jnk2=zeros(dsub); jnk2[flatsidx[idx2]]=ones(flatsidx[idx2]);
# verify
# plot(<<2000*jnk2; dsub>>);




# [azg] Fri Oct 19 19:33:20 PDT 2018#      plimits();
k1 = 0.00545053888 ;//  convert ADCcount to  [mm]
if(DEBUG == 1){
"idx2"
idx2
"flatsidx"
flatsidx
"diff(flatsidx)"
diff(flatsidx)
plimits();
plot(<< k1*dsub; 12*jnk2>>);
}


H=k1*mean(dsub[flatsidx[idx2[1]+1]-25:flatsidx[idx2[1]+1]]);
S1= k1*mean(dsub[flatsidx[idx2[2]+1]-25:flatsidx[idx2[2]+1]]);
S2 =  k1*mean(dsub[flatsidx[idx2[3]+1]-25:flatsidx[idx2[3]+1]]);
S3 = k1*mean(dsub[flatsidx[idx2[4]+1]-25:flatsidx[idx2[4]+1]]);
S4 =  k1*mean(dsub[flatsidx[idx2[5]+1]-25:flatsidx[idx2[5]+1]]);
H2 = k1*mean(dsub[flatsidx[idx2[6]+1]-25:flatsidx[idx2[6]+1]]);
# END  =  k1*mean(dsub[flatsidx[idx2[6]+1]+100:flatsidx[idx2[6]+1]+125]);
# END  =  k1*mean(dsub[flatsidx[idx2[6]+1]+100:flatsidx[idx2[6]+1]+122]);
END  =  k1*mean(dsub[flatsidx[idx2[6]+1]+100: length(dsub)]);


# error in misssteps
STEP1 = 672.5 - (S1 - H)/0.0127083;
STEP2 = 740 - (S1 - S2)/0.0127083;
STEP3 = 740 - (S3 - S2)/0.0127083;
STEP4 = 740 - (S3 - S4)/0.0127083;
ENDmissstep = 63 -(END - H2)/0.0127083;

if(DEBUG == 1){
printf("%s\n",fn);
printf("  STEP1err = %f,  STEP2err = %f,  STEP3err = %f,  STEP4err = %f, ENDerr = %f\n", STEP1,STEP2,STEP3,STEP4,ENDmissstep);
}













# return <<dsub=dsub; flatsidx=flatsidx; idx2=idx2>>;
return <<dsub=dsub; flatsidx=flatsidx; idx2=idx2; STEP1err = STEP1;  STEP2err = STEP2;  STEP3err = STEP3;  STEP4err = STEP4; ENDerr = ENDmissstep
>>;
};
