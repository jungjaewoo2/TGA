#!/home/pi/bin/rlab2
#rfile plplot
for(m in 1){
require runavg verify_wiring temperature 
# plot=plplot;
///////////////// NEW ////////////////////////////////////////////////
// while(1) {

# DEBUG2=1; // Testing at home without lasers 
DEBUG2=0; 

 


////////////////////////// AUG 27 ////////////////////////////////////
// add writing to a file:w
DATE=time2dstr(seconds(),"%Y-%m-%d");
TIME=time2dstr(seconds(),"%H:%M:%S");


# if(!isfile(fn1)){fprintf(fn1,"%s\n","SN,OPERATOR,DATE,TIME,STEP1,STEP2,STEP3,STEP4,Result");}
#
#if(isfile(fn1)){open(fn1,"a");}
#
#fprintf(fn1,"%s,%s,%s,%s,%2.2f,%2.2f,%2.2f,%2.2f,%s\n", SN,op,DATE,TIME,STEP1,STEP2,STEP3,STEP4,CMT);

// fprintf(fn1,"%s,%s,%s,%s,,\n", SN,op,DATE,TIME);
// fprintf(fn1,"%s,%s,%s,%s,%2.2f,%2.2f,%2.2f,%2.2f,%2.2f,%2.2f,%2.2f,\n", SN,op,DATE,TIME);
// close(fn1);
//////////////////////////////////////////////////////////////////////
# }
######################################################################
# printf("Twang test\n");


sleep(4); // make sure data is collected before processing starts

################# TWANG STARTS HERE ##################

# plclose(); // PGPLOT dies if window not closed ?

#!/home/pi/bin/rlab2
require expm damp 
# require expm damp plot 
# require rem expm lsim lsim2 tfeval ResRise ap5 damp plot subsample dfdp extract_chunk 
//////////////////////////////////////////////////////////////////////
# "line 2"
tic(4);
# tic();xxx=read_ascii("no_timestamp_production_vol=0.2.dat");
# xxx=read_ascii("/tmp/prod.dat");
# if( DEBUG2 == 1) { system("cp Original_with_reversed_wires_prod.dat /tmp/prod.dat;"); sleep 1; }
# if( DEBUG2 == 1) { system("cp OK.prod.dat /tmp/prod.dat;"); sleep 1; }
# if( DEBUG2 == 1) { system("cp prod3_slave.dat /tmp/prod.dat;"); sleep 1; }
if( DEBUG2 == 1) { system("cp slaveprod.dat /tmp/slaveprod.dat;"); sleep 1; }
dir="/tmp/";
file1="slaveprod.dat";
# xxx=read_ascii("/tmp/prodJune18_1.dat");
# dir="/home/azg/sandbox/usbreader/June26_2018/";
# file1="Original_with_reversed_wires_prod.dat";
xxx=read_ascii(dir + file1);
#######################################
# pol=verify_wiring(d);
# if(pol < 0 ) { 
# colors("red");
	# printf("%s ","STOP! WRONG POLARITY. Please check Coil #2 wires. Press " );
# colors("green");
	# printf("%s ","GREEN " );
# colors("red");
	# printf("%s \n","button to restart test from beginning" );

# blink_until_pressed();
# break; }
#######################################
DEBUG=0;
if(DEBUG == 1){printf("Data: %s\n",dir+file1);}

# printf("File read in %f [sec]\n",toc());
# flush first row
d=d[2:d.nr;];

dnew= d[;1] -sum(d[;1])/length(d[;1]);
th=100;
th=200;
# new to fix reversed wires 
# [azg] Wed Jun 27 00:52:26 PDT 2018
th=300 ;
# [azg] TGA needs lower th      Wed Apr 10 17:40:21 KST 2019
th=150;
dnew1=dnew - th;
dnew2=dnew + th;
dsign1=sign(dnew1);
dsign2=sign(dnew2);
ddiff1=abs(diff(dsign1));
ddiff2=abs(diff(dsign2));

ddiff = ddiff1 + ddiff2;

# find non-zero elements
ind=find(ddiff);
# " find((diff(ind) >5000) && (diff(ind) <6000 ))"
ind2=find((diff(ind) >5000) && (diff(ind) <6000 ));
 # for(i in 2:length(ind2)) { ind[ind2[i]] - ind[ind2[i-1]]  }

# " steps begin at:"
Rdc_extract_start=ind[find((diff(ind) >5000) && (diff(ind) <6000 ))];
# NOTE: the above finds only start of first and last long pulse
# Lets reconstruct all locations:
k1=Rdc_extract_start[1];
k4=Rdc_extract_start[2];
kstep=int((k4-k1)/3);
ind4=k1+(0:3)*kstep;


 # plot([4e-5*[1:length(d)]', [d;1]]);
 ##plot([4e-5*[ind[ind2[1]]:ind[ind2[4]]+7000]', [d[ind[ind2[1]]:ind[ind2[4]]+7000 ;1:2]]]);

ind[ind2]*4e-5;
k=1;

# for (idx in Rdc_extract_start) {
for (idx in ind4) {
	# limit the length to
	#LL = 800;
	# select part of settled waveform
	B2 = d[ idx +4000: idx+ 5000 ;];



	Vm=B2[;1];
	Vm=Vm * (511+1333)/511;
	Vm=Vm ./4096*3.3;
	Vp=B2[;3];
	Vp=Vp * (511+1333)/511;
	Vp=Vp ./4096*3.3;
	Im=B2[;4];
	dist=B2[;5]; # laser distance
	# Im=(Im - 2048) * (101+13+101)/(200*0.5*13)/4096*3.3;
	# Exact formula:
	Im = (Im -2048) /4096 *3.3 / (1/(1/0.5+1/(100+13+100)) * 13/(100+13+100) * 200 );
	Ip=B2[;2];
	# Ip=(Ip - 2048) * (101+13+101)/(200*0.5*13)/4096*3.3;
	# Exact formula:
	Ip = (Ip -2048) /4096 *3.3 / (1/(1/0.5+1/(100+13+100)) * 13/(100+13+100) * 200 );

        N =200;
	# L
	L[1]=1;
	L[2]=length(Vm);

	# "Rx="
	Rx[k] = abs(sum((Vp[L[2]-N:L[2]] -Vm[L[2]-N:L[2]])./(Im[L[2]-N:L[2]] )))/(N+1) -1.0;
	Ldist[k] = sum(dist[L[2]-N:L[2]] )/(N+1) ;
	Current[k] = sum(Im[L[2]-N:L[2]] )/(N+1) ;
	k=k+1;

}

for(j in 2:length(Rdc_extract_start)){ Rdc[j-1] = (Rx[j]+Rx[j-1])/2; } 

# wire resistance is 0.27 ohm
# Rdc_final=sum(Rx)/length(Rx) - 0.27;
# new value for improved PCB
# Wed Aug 29 20:18:25 PDT 2018
# Rdc_final=sum(Rx)/length(Rx) - 0.37;
# new value using 1% 13.0 ohm resistor
# Tue Oct 16 21:02:11 PDT 2018
# Rdc_final=sum(Rx)/length(Rx) - 0.63;
# Sun Mar 10 08:18:49 KST 2019
#RCAL=0.63;
#RDCMAX=14.0;
#RDCMIN=12.0;
# Rdc_final=sum(Rx)/length(Rx) - RCAL;
# RCAL will be subtracted in master
Rdc_final=sum(Rx)/length(Rx);



# DEBUG
# Rdc_final =  Rdc_final +1;


if( (Rdc_final > RDCMAX ) || ( Rdc_final < RDCMIN )){
	 colors("red");
	 printf ("Coil #2 resistance = %2.2f [Ohm]\t\tFAILED. Check connections.\n", Rdc_final );
	 colors();
 }



if( (Rdc_final < RDCMAX ) && ( Rdc_final > RDCMIN )){
	 colors("green");
	 printf ("Coil #2 resistance = %2.2f [Ohm]\t\tPASSED.\n", Rdc_final);
	 colors();
 }



# Leakage test

if (abs(sum(Ip-Im)) > 0.1* abs(sum(Im))) {
# if ((sum(Ip-Im)) > 0.1* sum(Im)) {
# DEBUG failure
# if ((sum(Ip-Im)) < 0.1* sum(Im)) {
 
	 colors("red");
	 printf ("Leakage test \t\t\t\tFAILED. Coil #2 shorted to ground.\n"); colors();
Short_pass="FAIL";
 }


if (abs(sum(Ip-Im)) < 0.1* abs(sum(Im))) {
	 colors("green");
	 printf ("Coil #2 leakage to ground\t\t\tPASSED.\n");
	 colors();
Short_pass="PASSED";
 }


# 	 printf ("Color check\n");
######################################################################

# " BEMF begin at:"
######################################################################
# BEMF_extract_start = ind[find((diff(ind) >2500) && (diff(ind) <3000 ))];
# the above failed when observed didtance was 2444
# NEW CODE START Sun Mar 10 09:50:36 KST 2019
# look only for short positive pulses
ind3=find(ddiff2);
#BEMF_extract_start = ind3[(find((diff(ind3) >24) && (diff(ind3) <32 )))];
# Wed Apr 24 13:26:18 KST 2019
BEMF_extract_start = ind3[(find((diff(ind3) >24) && (diff(ind3) <45 )))];
# k=max(find((diff(ind3) >24) && (diff(ind3) <32 ))) -1; // one before last
k = max(size(BEMF_extract_start))-1 ;
k = max(size(BEMF_extract_start))-2 ;

######################################################################
# " Laser begin at the same instant:"
Lsr_extract_start = BEMF_extract_start ;

# "Sanity check:"
# for(i in 2:length(BEMF_extract_start)) { BEMF_extract_start[i] - BEMF_extract_start[i-1]  }

skip=30;
skip=120; # *** last used
skip=80;
# skip=195;
# plot ([4e-5*[BEMF_extract_start[1]+skip : BEMF_extract_start[1]+1000 ]', [d[BEMF_extract_start[1]+skip : BEMF_extract_start[1]+1000;2]]]);

	                 
# k=6;
# k=7;
# [azg] Thu Oct 18 19:49:00 PDT 2018
# k = max(size(BEMF_extract_start))-3 ; // commented out Mon Mar 11 14:08:04 KST 2019

if(DEBUG == 1){printf("BEMF k= %i\n",k); }

# plot ([4e-5*[BEMF_extract_start[k]+skip : BEMF_extract_start[k]+1000 ]', [d[BEMF_extract_start[k]+skip : BEMF_extract_start[k]+1000;2]]]);

dc = sum(d[BEMF_extract_start[k]+500 : BEMF_extract_start[k]+1000;2]) / length(d[BEMF_extract_start[k]+500 : BEMF_extract_start[k]+1000;2]);
x = 4e-5*[BEMF_extract_start[k]+skip : BEMF_extract_start[k]+1000 ]';
x=x - x[1]; # reset time to zero
# y = d[BEMF_extract_start[k]+skip : BEMF_extract_start[k]+1000;2] -dc;

############################################################################
# [azg] Wed Apr 25 00:39:39 PDT 2018
# ss=5;
# y=subsample(y,ss);
# x=subsample(x,ss);

############################################################################




############################################################################
######################################################################
# plclose();
# " Lsr begin at:"
# Lsr_extract_start = ind[find((diff(ind) >2500) && (diff(ind) <3000 ))];

# Experimental data
# Conversion ratio of um to ADC count:
#                                      1200/(max(d[;5])-min(d[;5]))
um = 5.17;

skip=120;
skip=120+28; // +28 because now detecting leading and not trailing pulse edge 
# skip=75;
	                 

# Gain
# travel=5.45*abs(diff(Ldist))[1]; # in [um]
# gain = 1e-3*travel/abs(Current[1]);
# skew_travel=5.45*abs(diff(Ldist))[3]; # in [um]
# skew_gain = 1e-3*skew_travel/abs(Current[3]);
skew_travel=5.45*abs(Ldist[3] -Ldist[4]); # in [um]  May 08
skew_gain = 1e-3*skew_travel/abs(Current[3]-Current[4]);
beam_radius = (23.7 - 13.65)*1000;
skew_angle=skew_travel/beam_radius/pi*180;

// TRAVELmin=1200;
SKEW_TRAVELmin=193; # NOT CLEAR IF REQ'D
if(  skew_travel < SKEW_TRAVELmin ){
colors("red");
printf("Skew Travel range = %2.1f [um]\r\t\t\t\t\tFAILED\n",skew_travel);
colors();
}

if(  skew_travel  > SKEW_TRAVELmin ){
colors("green");
printf("Skew Travel range = %2.1f [um]\r\t\t\t\t\tPASSED.\n",skew_travel); colors();
}

// GAINmin=3.5;
SKEW_GAINmin=0.92; # [um/mA]
if(  skew_gain < SKEW_GAINmin ){
colors("red");
printf("Skew Gain = %2.2f [um/mA]\r\t\t\t\t\tFAILED\n",skew_gain);
colors();
}

if(  skew_gain > SKEW_GAINmin ){
colors("green");
printf("Skew Gain = %2.2f [um/mA]\r\t\t\t\t\tPASSED.\n",skew_gain); colors();
}


printf("Elapsed time %3.2f s\n",toc(4));
printf("########################################################################\n\n\n");

# fn1="/tmp/junk.csv";
# fn1="/home/pi/results/results.csv";
////////////////////////// AUG 28 ////////////////////////////////////
	// add writing to a file:w
DATE=time2dstr(seconds(),"%Y-%m-%d");
TIME=time2dstr(seconds(),"%H:%M:%S");
#*# op=reads("/var/www/html/Operator.txt");
# HACK to be removed
# SN=reads("/var/www/html/SN.txt");


# fn1="/home/pi/results/results_coarse_motion.csv";
# fn1="/var/www/html/"+DATE+"_twang.csv";
fn1="/tmp/slavetwang.csv";
# DEBUG
############ save the latest for debugging
#*# system("cp  /tmp/prod.dat /var/www/html/prod.dat"); sleep 1; 
############

////////////////////////// AUG 28 ////////////////////////////////////
////////////////////////// JAN 20 ////////////////////////////////////
#*# TADDR1="tbd";
#*# Tamb=temperature(TADDR1);
////////////////////////// JAN 20 ////////////////////////////////////



# if(!exist(fn1)){open(fn1,"a");fprintf(fn1,"%s\n","SN,,93.5698,3.96085,1.32148,13.3,,BEMFfres,BEMFdecay,BEMFresrise,Rdc_final,Short_pass,,LASERfres,LASERdecay,LASERresrise,Gain,travel");}

# if(!exist(fn1)){open(fn1,"a");fprintf(fn1,"%s\n","SN,,,,,,,BEMFfres,BEMFdecay,BEMFresrise,Rdc_final,Short_pass,,LASERfres,LASERdecay,LASERresrise,Gain,travel");}
# Tue Aug 28 21:13:44 PDT 2018

#*# if(isfile(fn1)){open(fn1,"a");}

#*# if(!isfile(fn1)){fprintf(fn1,"%s\n","SN,OPERATOR,DATE,TIME,BEMFfres,BEMFdecay,BEMFresrise,Rdc_final,Short_pass,,LASERfres,LASERdecay,LASERresrise,Gain,travel,Tambient");}


# LASERfres=fr;
# LASERdecay=1e3*zls.coef[1];
# LASERresrise=frrise;
Gain=gain;

# fprintf(fn1,"%s,,93.5698,3.96085,1.32148,13.3,,%2.2f,%2.2f,%2.2f,%2.2f,%s,,%2.2f,%2.2f,%2.2f,%2.2f,%2.2f,\n", SN, BEMFfres,BEMFdecay,BEMFresrise,Rdc_final,Short_pass,LASERfres,LASERdecay,LASERresrise,Gain,travel);

# fprintf(fn1,"%s,%s,%2.2f,%2.2f,%2.2f,%2.2f,%s,,%2.2f,%2.2f,%2.2f,%2.2f,%2.2f,\n", SN,op, BEMFfres,BEMFdecay,BEMFresrise,Rdc_final,Short_pass,LASERfres,LASERdecay,LASERresrise,Gain,travel);
# fprintf(fn1,"%s,%s,%s,%s,%2.2f,%2.2f,%2.2f,%2.2f,%s,,%2.2f,%2.2f,%2.2f,%2.2f,%2.2f,\n", SN,op, DATE, TIME, BEMFfres,BEMFdecay,BEMFresrise,Rdc_final,Short_pass,LASERfres,LASERdecay,LASERresrise,Gain,travel);
# fprintf(fn1,"%s,%s,%s,%s,%2.2f,%2.2f,%2.2f,%2.2f,%s,,%2.2f,%2.2f,%2.2f,%2.2f,%2.2f,%s,\n", SN,op, DATE, TIME, BEMFfres,BEMFdecay,BEMFresrise,Rdc_final,Short_pass,LASERfres,LASERdecay,LASERresrise,Gain,travel,Tamb);
#*# fprintf(fn1,"%s,%s,%s,%s,%2.2f,%2.2f,%2.2f,%2.2f,%s,,%2.2f,%2.2f,%2.2f,%2.2f,%2.2f,%s\n", SN,op, DATE, TIME, BEMFfres,BEMFdecay,BEMFresrise,Rdc_final,Short_pass,LASERfres,LASERdecay,LASERresrise,Gain,travel,Tamb);





#if(isfile(fn1)){open(fn1,"a");}

#if(!isfile(fn1)){fprintf(fn1,"%s\n","DATE,TIME,Rdc_final#2,Short_pass,SkewGain,SkewTravel");}

#fprintf(fn1,"%s,%s,%2.2f,%s,%2.2f,%2.2f\n", DATE, TIME, Rdc_final,Short_pass,Gain,travel);
##################################################
fn2="/tmp/slaveresults.r";
#if(isfile(fn2)){open(fn2,"a");}
fprintf(fn2,"Coil2Rdc=%2.2f;\nCoil2short_pass=\"%s\";\nSkewGain=%2.2f;\nSkewTravel=%2.2f;\n", Rdc_final,Short_pass,skew_gain,skew_travel);
close(fn2);
##################################################

# DEBUG ONLY
# printf("%s,%s,%2.2f,%2.2f,%2.2f,%2.2f,%s,,%2.2f,%2.2f,%2.2f,%2.2f,%2.2f,\n", SN,op, BEMFfres,BEMFdecay,BEMFresrise,Rdc_final,Short_pass,LASERfres,LASERdecay,LASERresrise,Gain,travel);

#close(fn1);
# sleep(10);

# clear(SN);
# clearall();

// DEBUG2=1; // Testing at home without lasers 

///////////////// NEW ////////////////////////////////////////////////
// }
# }
# }
# }

//////////////////////////////////////////////////////////////////////
}
