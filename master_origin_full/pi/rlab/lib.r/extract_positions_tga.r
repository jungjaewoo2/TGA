extract_positions_tga = function (fn, SN) {
// extract_positions_tga = function (fn) {
//
// fn -> filename or matrix with coarse laser data (from analog laser output)
// returns coarse motion locations: STEP0(HOME), STEP1, STEP2, STEP3, STEP4, STEP5(HOME), STEP6(END)
//
// Usage:  extract_positions("40073226.datcoarse.dat","SN")
//         extract_positions(d,SN)
// 
// Example:
//         for(i in ls("MAT16*.dat")){pltitle(i); err=extract_positions(i); pause();}
//
// 
////////////////////////////////////////////
// Last change: Fri Jun 17 11:38:50 KST 2022

# DEBUG = 1;
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
# th=max((abs(diff(dsub[2:length(dsub)])))) /2;
# th=max((abs(diff(dsub[2:400])))) /2;
skip=200;
skip=50; // Fri Jun 17 11:38:50 KST 2022
th=max((abs(diff(dsub[skip:length(dsub)])))) /2;
# [azg] Fri Nov  8 12:40:08 KST 2019
th=mean((abs(diff(dsub[skip:length(dsub)])))) ;
# signflat=sign( (abs(diff(dsub))) -11/2);
#signflat=sign( abs(diff(dsub[2:length(dsub)])) -th);
signflat=sign( abs(diff(dsub[skip:length(dsub)])) -th);

# [azg] needed? # signflat=signflat[2:length(signflat)];


# plot(<<abs(max(dsmooth)*diff(signflat)); dsmooth>>);
# pause();
# flatsidx=find(abs(diff(signflat))> 0.5);
flatsidx=find(abs(diff(signflat))> 0.5)+skip;

# [azg] # idx2 = find( diff(flatsidx)> 25 && diff(flatsidx)< 33);
idx2 = find( diff(flatsidx)> 25 && diff(flatsidx)< 40);
jnk2=zeros(dsub); jnk2[flatsidx[idx2]]=ones(flatsidx[idx2]);
# verify
 # plot(<<2000*jnk2; dsub>>);




# [azg] Fri Oct 19 19:33:20 PDT 2018#      plimits();
k1 = 0.00545053888 ;//  convert ADCcount to  [mm]
# NEW Wed May  4 18:03:43 KST 2022
k1 = 0.00538941866;  //  convert ADCcount to  [mm]
#####################################################
# [azg] Wed Jun 10 00:44:47 PDT 2020
t=40e-5*(1:length(d));
t=t[:];
um=[t,k1*d];

sleep(1);
box=reads("|uname -n ");
# RPi version
# printf("box = %s\n",box); // DEBUG
if((box != "bonnie.grzegorek.com") && (box != "muppet.grzegorek.com")) {
writem("/tmp/"+SN+".dat",um);
system("/usr/bin/7za  a -aou /var/www/html/"+SN+".zip /tmp/"+SN+".dat 2>&1 >/dev/null ");
system("/bin/rm /tmp/"+SN+".dat");
}

# test only
if((box == "bonnie.grzegorek.com") || (box == "muppet.grzegorek.com")) {
writem("/tmp/"+SN+".dat",um);
system("/usr/bin/7za  a -aou /tmp/"+SN+".zip /tmp/"+SN+".dat 2>&1 >/dev/null ");
system("/bin/rm /tmp/"+SN+".dat");
}
#####################################################
if(DEBUG == 1){
"idx2"
idx2
"flatsidx"
flatsidx
"diff(flatsidx)"
diff(flatsidx)
# plimits();
# plot(<< k1*dsub; 12*jnk2>>); // this is broken in plplot but works in pgplot
}


# H  = k1*mean(dsub[flatsidx[idx2[1]+1]-25:flatsidx[idx2[1]+1]]);
# S1 = k1*mean(dsub[flatsidx[idx2[2]+1]-25:flatsidx[idx2[2]+1]]);
# S2 = k1*mean(dsub[flatsidx[idx2[3]+1]-25:flatsidx[idx2[3]+1]]);
# S3 = k1*mean(dsub[flatsidx[idx2[4]+1]-25:flatsidx[idx2[4]+1]]);
# S4 = k1*mean(dsub[flatsidx[idx2[5]+1]-25:flatsidx[idx2[5]+1]]);
# H2 = k1*mean(dsub[flatsidx[idx2[6]+1]-25:flatsidx[idx2[6]+1]]);
# # END  =  k1*mean(dsub[flatsidx[idx2[6]+1]+100:flatsidx[idx2[6]+1]+125]);
# # END  =  k1*mean(dsub[flatsidx[idx2[6]+1]+100:flatsidx[idx2[6]+1]+122]);
# END  =  k1*mean(dsub[flatsidx[idx2[6]+1]+100: length(dsub)]);


# H   = k1*mean(dsub[flatsidx[idx2[1]]+7:flatsidx[idx2[1]]+29]);
# S1  = k1*mean(dsub[flatsidx[idx2[2]]+7:flatsidx[idx2[2]]+29]);
# S2  = k1*mean(dsub[flatsidx[idx2[3]]+7:flatsidx[idx2[3]]+29]);
# S3  = k1*mean(dsub[flatsidx[idx2[4]]+7:flatsidx[idx2[4]]+29]);
# S4  = k1*mean(dsub[flatsidx[idx2[5]]+7:flatsidx[idx2[5]]+29]);

# Sun Apr 28 13:37:01 KST 2019
endi=maxi(flatsidx); 
# Mon Jul 15 10:59:26 KST 2019
# fix calculation of H crashes for endi == 12;
# if( endi > 12 ){H=k1*mean(dsub[flatsidx[endi - 12]+5:flatsidx[endi -11 ]]);}
if( endi < 12 ){colors("red");stop("ERROR: Unable to determine coarse motion end positions");
colors();}
H=k1*mean(dsub[flatsidx[endi - 11]-23:flatsidx[endi -11 ]-2]);
S0=k1*mean(dsub[flatsidx[endi - 10]+2:flatsidx[endi -9 ]-2]);
S1=k1*mean(dsub[flatsidx[endi - 8]+2:flatsidx[endi -7 ]-2]);
S2=k1*mean(dsub[flatsidx[endi - 6]+2:flatsidx[endi -5 ]-2]);
S3=k1*mean(dsub[flatsidx[endi - 4]+2:flatsidx[endi -3 ]-2]);
S4=k1*mean(dsub[flatsidx[endi - 2]+2:flatsidx[endi -1 ]-2]);




# # error in misssteps
# STEP1 = 672.5 - (S1 - H)/0.0127083;
# STEP2 = 740 - (S1 - S2)/0.0127083;
# STEP3 = 740 - (S3 - S2)/0.0127083;
# STEP4 = 740 - (S3 - S4)/0.0127083;
# ENDmissstep = 63 -(END - H2)/0.0127083;

# error in misssteps
# STEP1 = 1870 - (S1 - S0)/0.0051;
# STEP2 = 1870 - (S1 - S2)/0.0051;
# STEP3 = 1870 - (S3 - S2)/0.0051;
# STEP4 = 1870 - (S3 - S4)/0.0051;
#ENDmissstep = 63 -(END - H2)/0.0051;

# error in misssteps
STEP1 = 1569 - (S1 - S0)/0.00513850107;
STEP2 = 1569 - (S1 - S2)/0.00513850107;
STEP3 = 1569 - (S3 - S2)/0.00513850107;
STEP4 = 1569 - (S3 - S4)/0.00513850107;
#ENDmissstep = 63 -(END - H2)/0.00513850107;

if(DEBUG == 1){
printf("%s\n",fn);
printf("  H = %f,  S0 = %f,  S1 = %f,  S2 = %f,  S3 = %f ,  S4 = %f \n", H,S0,S1,S2,S3,S4);
printf("  STEP1err = %f,  STEP2err = %f,  STEP3err = %f,  STEP4err = %f \n", STEP1,STEP2,STEP3,STEP4);
# printf("  STEP1err = %f,  STEP2err = %f,  STEP3err = %f,  STEP4err = %f, ENDerr = %f\n", STEP1,STEP2,STEP3,STEP4,ENDmissstep);
}



# return <<dsub=dsub; flatsidx=flatsidx; idx2=idx2>>;
# return <<dsub=dsub; flatsidx=flatsidx; idx2=idx2; STEP1err = STEP1;  STEP2err = STEP2;  STEP3err = STEP3;  STEP4err = STEP4; ENDerr = ENDmissstep
return <<dsub=dsub; flatsidx=flatsidx; idx2=idx2; STEP1err = STEP1;  STEP2err = STEP2;  STEP3err = STEP3;  STEP4err = STEP4 >>;
};
