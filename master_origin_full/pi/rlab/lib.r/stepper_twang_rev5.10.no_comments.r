
# Last changed on: Thu 14 Jul 2022 03:49:40 PM PDT
DEBUG=0;
SAVEBAD=1;
# rfile plplot; # [azg] Sat Mar  5 11:19:45 KST 2022
colors();
CR=[]; # [azg] Wed Jul 20 13:05:50 KST 2022
LOOP=10;
LOOP=3;
LOOP=5;
rfile specs;
if ( DEBUG2 == 1 ) {
	system("cp /home/pi/sandbox/tf_fitting_rfiles/complete/Oct_21_test/OK.prod.dat /tmp/prod.dat");
	system("touch /tmp/ADCDONE");
	T=19.9;
	#system("cp /home/pi/sandbox/datcoarse.dat /tmp"); 
	#system("cp /home/pi/sandbox/2019-11-10_datcoarse.dat /tmp/datcoarse.dat"); 
}

box=reads("|uname -n "); 

SKIPCOARSE=0;
rfile loops;
if(RET == 2) {RET=1; SKIPCOARSE=1;}

LOOP=RET;
for(m in (LOOP:1:-1) ){
require savepdf3_tga coarse_motion_linearity5_tga extract_positions_tga runavg subsample blink_until_pressed verify_wiring temperature specs Rdc_20degC
plot=plplot;
///////////////// NEW ////////////////////////////////////////////////
// while(1) {
	run=LOOP-m+1;
	if( (SAVEBAD == 1) && (LOOP > 1) && exist(LOOP) ) {printf("Starting run %i out of %i\n",run,LOOP);}

// DEBUG2=1; // Testing at home without lasers (shoud be set in specs.r
///////////////// NEW DEBUG CRASHES///////////////////////////////////
// Save the results of previous run
if((run == 1) && (LOOP == 200000 ) && ( SAVEBAD == 1 )){

//"Save crash if())"
DATE=time2dstr(seconds(),"%Y-%m-%d");
TIME=time2dstr(seconds(),"%H:%M:%S");
ORUN=1;
if( isfile("/tmp/oldrun.mat") ) {ORUN=readm("/tmp/oldrun.mat");}

if((box == "bonnie.grzegorek.com") || (box == "muppet.grzegorek.com")) {
	if(!exist("/tmp/tmp")) {mkdir("/tmp/tmp");}
fb1="/tmp/tmp/"+DATE+"_"+TIME+"_run="+num2str(ORUN)+"_datcoarse.dat";
fb2="/tmp/tmp/"+DATE+"_"+TIME+"_run="+num2str(ORUN)+"_prod.dat";
}

if((box != "bonnie.grzegorek.com") && (box != "muppet.grzegorek.com")) { 
fb1="/home/pi/tmp/"+DATE+"_"+TIME+"_run="+num2str(ORUN)+"_datcoarse.dat";
fb2="/home/pi/tmp/"+DATE+"_"+TIME+"_run="+num2str(ORUN)+"_prod.dat";
}

if( isfile("/tmp/datcoarse.dat") ){system("cp /tmp/datcoarse.dat "+fb1);}
if( isfile("/tmp/prod.dat") ){system("cp /tmp/prod.dat "+fb2);}

}
if(!isfile("/tmp/ADCDONE") && isfile("/tmp/datcoarse.dat") ){sleep(60);} 

writem("/tmp/oldrun.mat",run);

///////////////// END DEBUG CRASHES///////////////////////////////////
if((box != "bonnie.grzegorek.com") && (box != "muppet.grzegorek.com")) { 
	############
TADDR1=grep(ls("/sys/bus/w1/devices/"),"28-"); # Temp sensor detection
	############
Tamb=temperature(TADDR1);
}
if(Tamb[1;1] == "NaN"){colors("red");pause("ERROR: Temperature probe unplugged? Press RED button to power off.")}
colors();


// Fake Temp for testing
if( DEBUG2 == 1) {Tamb="19.9999"; }

# orig dies when Tamb="N/A" (missing or not ready sensor) # T=eval(Tamb);
if( Tamb == "N/A" ){colors("red");printf("\nERROR: Failed to read temperature.\n");colors();}
if( Tamb != "N/A" ){T=eval(Tamb);}

if(DEBUG2 != 1){
//////////////////////////////////////////////////////////////////////

// Open 3 ports (Laser controller + stepper controller + barcode reader)
tic(1);
tic(2);

// Laser controller //////////////////////////////////////////////////
attr=<<>>;
attr.devname = "ttyUSB";
attr.id_vendor_id = "0403";

if(!exist(r1)){
r1 = usbdevname(attr);
}

 
if(DEBUG2 != 1){
if(isempty(r1)){colors("red");pause("ERROR: Laser controller unplugged? Press RED button to power off.")}
colors();
}


opts=<<>>;
opts.data_parity_stop="8N1";
opts.eol="\r";
opts.speed=115200; 
# opts.speed=9600; # misconfigured controller
opts.flow_control="n";

serial="serial://" + r1;

if(!open(serial,opts)) {
open(serial, opts);
}


readdist = ["MS,02\r"];
datarange = [7:14];
commmode =   ["Q0\r"];   // Measurements stop during communication mode
generalmode = ["R0\r"];  // Measurements control commands are accepted
analogchan21  = ["SW,CG,02,01\r"]; // laser 02 to analog out 1
analogchan11  = ["SW,CG,01,01\r"]; // laser 01 to analog out 1
analogchan22  = ["SW,CG,02,02\r"]; // laser 02 to analog out 2
analogchan12  = ["SW,CG,01,02\r"]; // laser 01 to analog out 2
//////////////////////////////////////////////////////////////////////
// Storage commands
storage_on = ["SW,OK,02,1\r"]; // ONLY in Communication Mode
storage_off= ["SW,OK,02,0\r"]; // for laser2 ONLY in Communication Mode
sstart = ["AS\r"];      // ONLY in General Mode
sstop  = ["AP\r"];      // ONLY in General Mode
sout   = ["AO,02\r"]; // response AO,hhhhhhhh...,hhhhhhhh ONLY in General Mode

//////////////////////////////////////////////////////////////////////
// switch analog output 1 to laser 2 /////////////////////////////////
// Do we need switch to command mode first?
 writem(serial, commmode);
      r = readm(serial);
      // verify if out is laser2
 writem(serial, analogchan11);
      r = readm(serial);
      //switch  back to measuremnt mode
 writem(serial, generalmode);
      r = readm(serial);
      
      
//////////////////////////////////////////////////////////////////////

// stepper controller ////////////////////////////////////////////////

attr=<<>>;
attr.devname = "ttyUSB";
attr.id_vendor_id = "10c4";

if(!exist(r2)){
r2 = usbdevname(attr);
}
 
if(isempty(r2)){colors("red");pause("ERROR: Stepper controller unplugged? Press RED button to power off.")}
colors();

opts=<<>>;
opts.data_parity_stop="8N1";
opts.eol="\r";
opts.speed=9600;
opts.flow_control="n";

stepper="serial://" + r2;

if(!open(stepper,opts)) {
open(stepper, opts);
}





total=["/1v100Z30000z400M1000V600P318D100M1000P3138M1000D3138M1000P3138M1000D3138M1000P1570R"]; # + 60%
total=["/1v100Z30000z400M1000V600P418D200M1000P3138M1000D3138M1000P3138M1000D3138M1000P1570R"]; 
#[azg]  accelerate home search Wed Jun 15 12:46:59 KST 2022
total=["/1v800Z30000z400M1000V600P418D200M1000P3138M1000D3138M1000P3138M1000D3138M1000P1570R"]; # 
# total=["/1v400Z30000z400M1000V600P418D200M1000P3138M1000D3138M1000P3138M1000D3138M1000P1570R"]; # 
# total=["/1v200Z30000z400M1000V600P418D200M1000P3138M1000D3138M1000P3138M1000D3138M1000P1570R"]; # 

# TEST VERSION to UNDO over extended elevator
# total=["/1D2000v800Z30000z400M1000V600P418D200M1000P3138M1000D3138M1000P3138M1000D3138M1000P1570R"]; # 

///////////// Barcode reader /////////////////////////////////////////
 # attr.id_vendor_id = "1eab";
 
 
attr=<<>>;
attr.devname = "ttyACM";
attr.id_vendor_id = "0720"; # Keyence
# BELOW TEST ONLY (TO BE DELETED)
# attr.id_vendor_id = "0483"; # Symcode 0483:5740

 
 
if(!exist(r3)){
r3 = usbdevname(attr);
}
 
if(isempty(r3)){colors("red");pause("ERROR: Barcode reader unplugged? Press RED button to power off.")}
colors();


} 

if( (!exist(SN)) && ( DEBUG2 == 1)) { SN="Fake1"; } ###################################### # [azg] BEGIN Sat May 14 11:49:57 KST 2022 if( !exist(SN) && isfile("/tmp/SN.txt") && ( LOOP == 5 || LOOP == 200000 )){SN=reads("/tmp/SN.txt");} # if( !exist(SN) && isfile("/tmp/SN.txt") && ( m != 1){SN=reads("/tmp/SN.txt");} # [azg] END Sat May 14 11:49:57 KST 2022 ######################################
 
if(!exist(SN)){
colors("red");
printf("Waiting for S/N...");
 
r=getline(r3,16);
SN=sum (strsplt(r)[1:length(strsplt(r)) -1]);

tic(3);
if( (LOOP == 200000) || (LOOP == 5 )){ writem( "/tmp/SN.txt",SN); } 
} 
colors("green");
printf("\rSN:              \t%s\n",SN);
colors();

////////////
 gmt= date2jd(gmtime());
 gmt= gmt -floor(gmt);
 gmt=100000*gmt;
 sprintf(dd,"%i.dat",gmt);
//////////////////////////////////////////////////////////////////////
	######### start of SKIPCOARSE ########
# SKIPCOARSE=1;
if( SKIPCOARSE != 1) {
printf("COARSE MOTION TEST:\n");
if(DEBUG2 != 1){
	# BROKEN if LOOP == 5
if(!exist(LOOP) || (LOOP == 1) || (LOOP == 5) ) {if(run == 1){blink_until_pressed();} }
}


//////////////////////////////////////////////////////////////////////


if(DEBUG2 != 1){
# system("sudo nice -20 /home/pi/bin/coarse_60s_MCP3202 > /tmp/datcoarse.dat&");
system("sudo nice -20 /home/pi/bin/coarse_40s_MCP3202 > /tmp/datcoarse.dat&");
}

if(DEBUG2 == 1){
}

sleep(0.5);



if(DEBUG2 != 1){

 writem(stepper, total);
      r = readm(stepper);
      #printf("%s\n",r);


 r = readm(stepper);
      #printf("%s\n",r);







sleep(1);
while(!isfile("/tmp/ADCDONE")){sleep(1);} 
// printf("%s\n", "ready");
 
// switch back analog output to laser 1
// sleep(10);
// switch back laser 1 to analog output 1
 writem(serial, commmode);
      r = readm(serial);
      #printf("%s\n",r);
      // verify if out is laser2
 writem(serial, analogchan11);
      r = readm(serial);
      #printf("%s\n",r);
 writem(serial, analogchan22);
      r = readm(serial);
      #printf("%s\n",r);
      //switch  back to measuremnt mode
 writem(serial, generalmode);
      r = readm(serial);
      #printf("%s\n",r);
 # remove below } if not DEBUG2
}
      # sleep(2);
/////////////////////////// AUG 27 ///////////////////////////////////
	printf("reading measured laser data\n");
if( DEBUG2 == 1) { 
	system("cp /home/pi/sandbox/2019-11-10_datcoarse.dat /tmp/datcoarse.dat"); 
	sleep(1); 
        system("touch /tmp/ADCDONE");
}


fbad="/home/pi/for_debugging/datcoarse_run="+num2str(run)+".dat;";
if( SAVEBAD == 1) { system("cp /tmp/datcoarse.dat "+fbad); }

xxx=read_ascii("/tmp/datcoarse.dat");
d=d[2:d.nr;];


// E= extract_positions_tga(d);
E= extract_positions_tga(d, SN);
if(DEBUG2 ==1){
printf("Before correcting laser error LCORR\n");
printf("S1=%f  S2=%f  S3=%f  S4=%f\n",E.STEP1err, E.STEP2err, E.STEP3err, E.STEP4err);
endi=maxi(E.flatsidx); printf("endi=%i \n\n",endi);
}



P4=coarse_motion_linearity5_tga(E.flatsidx ,SN);


savepdf3_tga(P4, SN);
if(!exist(LCORR)){LCORR = -4.43180556;}

STEP1= E.STEP1err + LCORR;
STEP2= E.STEP2err + LCORR;
STEP3= E.STEP3err + LCORR;
STEP4= E.STEP4err + LCORR;


////////////////////////// FEB 15 2019  ////////////////////////////////////
// UNDER = 12;
// OVER = -4;
CR=[]; // coarse motion result
if( ( STEP1 <= UNDER ) && (STEP1 >=  OVER  )){ CR[1]=1; else CR[1]=0;}
if( ( STEP2 <= UNDER ) && (STEP2 >=  OVER  )){ CR[2]=1; else CR[2]=0;}
if( ( STEP3 <= UNDER ) && (STEP3 >=  OVER  )){ CR[3]=1; else CR[3]=0;}
if( ( STEP4 <= UNDER ) && (STEP4 >=  OVER  )){ CR[4]=1; else CR[4]=0;}

if( prod(CR) == 1 ){ 
CMT="PASS"; 
colors("green");
printf("Coarse motion test:\t\t\tPASSED.\n");
else 
CMT="FAIL";
colors("red");
printf("Coarse motion test:\t\t\tFAILED\n");
}

colors();


////////////////////////// AUG 27 ////////////////////////////////////
// add writing to a file:w
DATE=time2dstr(seconds(),"%Y-%m-%d");
TIME=time2dstr(seconds(),"%H:%M:%S");
op="Unknown";
box=reads("|uname -n "); 
if((box != "bonnie.grzegorek.com") && (box != "muppet.grzegorek.com")) { 
if( (prod(CR) == 1) && ( SAVEBAD == 1) ){ system("/bin/rm "+fbad); }
op=reads("/var/www/html/Operator.txt");
}


fn1="/var/www/html/"+DATE+"_coarse_motion.csv";
if((box == "bonnie.grzegorek.com") || (box == "muppet.grzegorek.com")) { fn1="/tmp/"+DATE+"_coarse_motion.csv";}



if(!isfile(fn1)){fprintf(fn1,"%s\n","SN,OPERATOR,DATE,TIME,STEP1,STEP2,STEP3,STEP4,Result");}

if(isfile(fn1)){open(fn1,"a");}

fprintf(fn1,"%s,%s,%s,%s,%2.2f,%2.2f,%2.2f,%2.2f,%s\n", SN,op,DATE,TIME,STEP1,STEP2,STEP3,STEP4,CMT);

close(fn1);
//////////////////////////////////////////////////////////////////////
clear(d);
}
####### end of SKIPCOARSE ########

if( SKIPCOARSE == 1 ) { blink_until_pressed(); }
printf("Continue with Twang test\n");

if(DEBUG2 != 1){
 writem(serial, commmode);
      r = readm(serial);
      #printf("%s\n",r);
      // verify if out is laser2
writem(serial, analogchan12);
      //switch  back to measuremnt mode
 writem(serial, generalmode);
      r = readm(serial);
}



if((box != "bonnie.grzegorek.com") && (box != "muppet.grzegorek.com")) { 
system("ssh slave 'sudo nice -20 /home/pi/bin/prodMCP3202 > /tmp/slaveprod.dat 2>/dev/null'&");
system("sudo nice -20 /home/pi/bin/prodMCP3202 > /tmp/prod.dat 2>/dev/null&");
system(" ssh gen 'play -q bin/5Hz_1ms_TGA.flac vol 0.34' &");

sleep(5); // make sure data is collected before processing starts
system("ssh slave '/home/pi/bin/rlab2 /home/pi/bin/slave_twang.r'&");
}



plclose(); // PGPLOT dies if window not closed ?

require rem expm lsim lsim2 tfeval ResRise ap5 damp plot subsample dfdp extract_chunk 
//////////////////////////////////////////////////////////////////////
tic(4);
if( DEBUG2 == 1) { system("cp /home/pi/sandbox/tf_fitting_rfiles/complete/Jan_21_2019/OK.prod.dat /tmp/prod.dat;"); sleep(1); }
dir="/tmp/";
file1="prod.dat";
xxx=read_ascii(dir + file1);
pol=verify_wiring(d);
if(pol < 0 ) { 
colors("red");
	printf("%s ","STOP! WRONG DIRECTION OF MOTION. Suspect magnets or coils. Press " );
colors("green");
	printf("%s ","GREEN " );
colors("red");
	printf("%s \n","button to restart test from beginning" );

if(LOOP == 1){clear(SN);}
if(LOOP == 1){ blink_until_pressed(); break; }
}

if(DEBUG == 1){printf("Data: %s\n",dir+file1);}

d=d[2:d.nr;];

dnew= d[;1] -sum(d[;1])/length(d[;1]);
th=100;
th=200;
th=300 ;
th=150;
dnew1=dnew - th;
dnew2=dnew + th;
dsign1=sign(dnew1);
dsign2=sign(dnew2);
ddiff1=abs(diff(dsign1));
ddiff2=abs(diff(dsign2));

ddiff = ddiff1 + ddiff2;

ind=find(ddiff);
ind2=find((diff(ind) >5000) && (diff(ind) <6000 ));
 # for(i in 2:length(ind2)) { ind[ind2[i]] - ind[ind2[i-1]]  }

Rdc_extract_start=ind[find((diff(ind) >5000) && (diff(ind) <6000 ))];



ind[ind2]*4e-5;
k=1;

for (idx in Rdc_extract_start) {
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
	Rx[k] = sum((Vp[L[2]-N:L[2]] -Vm[L[2]-N:L[2]])./(Im[L[2]-N:L[2]] ))/(N+1) -1.0;
	Ldist[k] = sum(dist[L[2]-N:L[2]] )/(N+1) ;
	Current[k] = sum(Im[L[2]-N:L[2]] )/(N+1) ;
	k=k+1;

}

for(j in 2:length(Rdc_extract_start)){ Rdc[j-1] = (Rx[j]+Rx[j-1])/2; } 

Rdc_final=sum(Rx)/length(Rx) - RCAL1;
if( isempty(T) ) { colors("red"); printf("\n %s \n", "ERROR: Failed to read temperature"); colors();}
Coil1Rdc_20deg = Rdc_20degC( Rdc_final+RCAL1, T) - RCAL1_20C;

if( abs(diff(Ldist))[1] < abs(diff(Ldist))[3]) { 
colors("red");
	printf("%s ","STOP! WRONG DIRECTION OF MOTION. Suspect magnets or coils. Press " );
colors("green");
	printf("%s ","GREEN " );
colors("red");
	printf("%s \n","button to restart test from beginning" );

clear(SN);
blink_until_pressed();
break; 
}

if(( Coil1Rdc_20deg > RDCMAX ) || ( Coil1Rdc_20deg  < RDCMIN )){
	 colors("red");
	 printf ("Coil #1 resistance @20°C = %2.2f [Ohm]\tFAILED. Check connections.\n",Coil1Rdc_20deg  );
	 CR[length(CR)+1]=0;
	 CR[length(CR)+1]=0;
	 colors();
 }



if(( Coil1Rdc_20deg < RDCMAX ) && ( Coil1Rdc_20deg > RDCMIN )){
	 colors("green");
	 printf ("Coil #1 resistance @20°C = %2.2f [Ohm]\tPASSED.\n", Coil1Rdc_20deg );
	 CR[length(CR)+1]=1;
	 colors();
 }




if (abs(sum(Ip-Im)) > 0.1* abs(sum(Im))) {
 
	 colors("red");
	 printf ("Leakage test \t\t\t\tFAILED. Coil #1 shorted to ground.\n"); colors();
	 CR[length(CR)+1]=0;
Short_pass="FAIL";
 }


if (abs(sum(Ip-Im)) < 0.1* abs(sum(Im))) {
	 colors("green");
	 printf ("Coil #1 leakage to ground\t\tPASSED.\n");
	 CR[length(CR)+1]=1;
	 colors();
Short_pass="PASSED";
 }


ind3=find(ddiff2);
BEMF_extract_start = ind3[(find((diff(ind3) >24) && (diff(ind3) <32 )))];
k = max(size(BEMF_extract_start))-1 ;
k = max(size(BEMF_extract_start))-2 ;

Lsr_extract_start = BEMF_extract_start ;


skip=80;

	                 

if(DEBUG == 1){printf("BEMF k= %i\n",k); }

plot ([4e-5*[BEMF_extract_start[k]+skip : BEMF_extract_start[k]+1000 ]', [d[BEMF_extract_start[k]+skip : BEMF_extract_start[k]+1000;2]]]);

      # "First plot displayed"
      # "size(d)"
      # size(d)
dc = sum(d[BEMF_extract_start[k]+500 : BEMF_extract_start[k]+1000;2]) / length(d[BEMF_extract_start[k]+500 : BEMF_extract_start[k]+1000;2]);
x = 4e-5*[BEMF_extract_start[k]+skip : BEMF_extract_start[k]+1000 ]';
x=x - x[1]; # reset time to zero
y = d[BEMF_extract_start[k]+skip : BEMF_extract_start[k]+1000;2] -dc;

if(!exist(func)){
func = function(x,p){
	global(u2,pi)

    return p[4] + p[3]*exp(-(x)/p[1]).*sin(2*pi*p[2].*(x-p[5]));
    // return p[4] + p[3]*exp(-(x)/p[1]).*sin(2*pi*p[2].*(x-p[5]));
    };

}




// nonlinear fit
// initial guess
p=[ 0.0111863796,    65.4338127,    90.1738737,   -2.95681298,  -0.00143286012 ]';
p=pBEMF;



xy = [x,y];

options=<<>>;
options.imethod = 2;  // least squares fit

 zls=lsfit(y,x,p,func,dfdp);
BEMFzls=zls;




plclose(); // NEW Sat Apr 13 08:17:20 KST 2019
pltitle("Fitting measured data ");
plegend(["Data","Lsfit"]);
if (_rlab_config.plot_support=="pgplot")
	    {
              plaxis("lin",); # pgplot syntax
              xlabel("Time [s]");
              ylabel("Amplitude");
     else
              plscale("lin",); # plplot syntax
              plxlabel("Time [s]");
              plylabel("Amplitude");
   }

zzls=(func(x,zls.coef));


plot(<<[x,y];[x,zzls]>> );
plt1debug=(<<[x,y];[x,zzls]>> );
sleep(3);



lsnum2=[1,2*pi*zls.coef[2]];
lsden2=[1,2/zls.coef[1],1/(zls.coef[1]*zls.coef[1])+2*pi*2*pi*(zls.coef[2]*zls.coef[2])]; 

lsnum2=[zls.coef[1],0];
lsden2=[1,2/zls.coef[1],1/(zls.coef[1]*zls.coef[1])+2*pi*2*pi*(zls.coef[2]*zls.coef[2])]; 


dnum2=[zls.coef[1]];
dden2= lsden2;

freq = [10:3000:5];
inum2 = [1, 0];
itf = tfeval(inum2, lsden2, freq);
idb=20*log10(abs(itf));



plclose(); // NEW Sat Apr 13 08:17:20 KST 2019
if (_rlab_config.plot_support=="pgplot")
	    {
               plaxis("log",); # pgplot syntax
              xlabel("Frequency [Hz]");
              ylabel("Amplitude [dB]");
     else
              plscale("log",); # plplot syntax
              plxlabel("Frequency [Hz]");
              plylabel("Amplitude [dB]");
   }

plegend();
pltitle("Frequency response");
plegend(["Back EMF", "Input data"]);


dtf = tfeval(dnum2, dden2, freq);
dtf0 = tfeval(dnum2, dden2, freq[1]);

dispdb = 20*log10(abs(dtf)) - 20*log10(abs(dtf0));

plclose(); // NEW Sat Apr 13 08:17:20 KST 2019
pltitle("Displacement Frequency Response");
if (_rlab_config.plot_support=="pgplot")
	    {
               plaxis("log",); # pgplot syntax
               xlabel("Frequency [Hz]");
               ylabel("Amplitude [dB]");
     else
              plscale("log",); # plplot syntax
              plxlabel("Frequency [Hz]");
              plylabel("Amplitude [dB]");
   }

   
   # plscale("log",); # plplot syntax
plegend(["Input data","Displacement"]);
plegend(["Displacement"]);

plot(<< [freq;dispdb][;1:30]'>>);
sleep(3);


zeta_extracted = 1./(zls.coef[1]*sqrt(1/(zls.coef[1]*zls.coef[1])+2*pi*2*pi*zls.coef[2]*zls.coef[2]));


wn_extracted = zls.coef[2].*(2*pi) ./ sqrt(1- zeta_extracted[1]*zeta_extracted[1]) ;


den_extracted = [1, 2*zeta_extracted*wn_extracted,  wn_extracted^2 ]; 

X_extracted = damp (den_extracted);
X=damp(lsden2);


colors("yellow");
printf("From BackEMF Measurements:\n");
colors();

fr=X_extracted.wn[1]/(2*pi) ;
printf("Frequency Response Peak = %2.2f [dB]\n",X_extracted.FRpeaking[1]);




     #printf("\n****************************************************\n");
fr=X_extracted.wn[1]/(2*pi) ;
if( (fr <FMAX) && (fr > FMIN )){
colors("green");
printf("Resonant frequency = %2.1f [Hz]\t\tPASSED.\n",fr);
	 CR[length(CR)+1]=1;
colors();
}


if( (fr >FMAX) || (fr < FMIN )){
colors("red");
printf("Resonant frequency = %2.1f [Hz]\t\tFAILED.\n",fr);
	 CR[length(CR)+1]=0;
colors();
}



frrise = ResRise (1e3*zls.coef[1], fr);


printf("Resonant rise = %2.2f [dB] \n",frrise,,);

// # printf("Elapsed time %3.2f s\n",toc(4));

BEMFfres=fr;
BEMFdecay=1e3*zls.coef[1];
BEMFresrise=frrise;


plclose();

um = 5.17;

skip=120;
skip=120+28; // +28 because now detecting leading and not trailing pulse edge 
skip=80; // Wed Jul 17 12:12:14 KST 2019
	                 
if(DEBUG == 1){printf("Laser k= %i\n",k); }

plclose(); // NEW Sat Apr 13 08:17:20 KST 2019
pltitle("Fitting measured laser data ");
if (_rlab_config.plot_support=="pgplot")
	    {
               plaxis("lin",); # pgplot syntax
              xlabel("Time [s]");
              ylabel("Amplitude [um]");
     else
              plscale("lin",); # plplot syntax
              plxlabel("Time [s]");
              plylabel("Amplitude [um]");
   }
plot ([4e-5*[Lsr_extract_start[k]+skip : Lsr_extract_start[k]+1000 ]', [um*d[Lsr_extract_start[k]+skip : Lsr_extract_start[k]+1000;5]]]);

dc = sum(d[Lsr_extract_start[k]+500 : Lsr_extract_start[k]+1000;5]) / length(d[Lsr_extract_start[k]+500 : Lsr_extract_start[k]+1000;5]);
x = 4e-5*[Lsr_extract_start[k]+skip : Lsr_extract_start[k]+1000 ]';
x=x - x[1]; # reset time to zero
y = d[Lsr_extract_start[k]+skip : Lsr_extract_start[k]+1000;5] -dc;



// nonlinear fit
// initial guess
p=[ 0.0105749365,      65.70088,    56.3116544,   -1.47410849,  -0.00294813441 ]';
p=pLASER;



xy = [x,y];

options=<<>>;
options.imethod = 2;  // least squares fit

 zls=lsfit(y,x,p,func,dfdp);



plclose(); // NEW Sat Apr 13 08:17:20 KST 2019
pltitle("Fitting laser measured data ");
plegend(["Data","Lsfit"]);

if (_rlab_config.plot_support=="pgplot")
	    {
               plaxis("lin",); # pgplot syntax
              xlabel("Time [s]");
              ylabel("Amplitude [um]");
     else
              plscale("lin",); # plplot syntax
              plxlabel("Time [s]");
              plylabel("Amplitude [um]");
   }

zzls=(func(x,zls.coef));


plot(<<[x,um*y];[x,um*zzls]>> );
sleep(3);



lsnum2=[1,2*pi*zls.coef[2]];
lsden2=[1,2/zls.coef[1],1/(zls.coef[1]*zls.coef[1])+2*pi*2*pi*(zls.coef[2]*zls.coef[2])]; 

lsnum2=[zls.coef[1],0];
lsden2=[1,2/zls.coef[1],1/(zls.coef[1]*zls.coef[1])+2*pi*2*pi*(zls.coef[2]*zls.coef[2])]; 


dnum2=[zls.coef[1]];
dden2= lsden2;

freq = [10:3000:5];
inum2 = [1, 0];
itf = tfeval(inum2, lsden2, freq);
idb=20*log10(abs(itf));



plclose(); // NEW Sat Apr 13 08:17:20 KST 2019
if (_rlab_config.plot_support=="pgplot")
	    {
               plaxis("log",); # pgplot syntax
              xlabel("Frequency [Hz]");
              ylabel("Amplitude [dB]");
     else
              plscale("log",); # plplot syntax
              plxlabel("Frequency [Hz]");
              plylabel("Amplitude [dB]");
   }

plegend();
pltitle("Frequency response");


dtf = tfeval(dnum2, dden2, freq);
dtf0 = tfeval(dnum2, dden2, freq[1]);

dispdb = 20*log10(abs(dtf)) - 20*log10(abs(dtf0));

plclose(); // NEW Sat Apr 13 08:17:20 KST 2019
pltitle("Displacement Frequency Response");

if (_rlab_config.plot_support=="pgplot")
	    {
               plaxis("log",); # pgplot syntax
              xlabel("Frequency [Hz]");
              ylabel("Amplitude [dB]");
     else
              plscale("log",); # plplot syntax
              plxlabel("Frequency [Hz]");
              plylabel("Amplitude [dB]");
   }

plegend(["Input data","Displacement"]);
plegend(["Displacement"]);

plot(<< [freq;dispdb][;1:30]'>>);
sleep(2);
plclose();


zeta_extracted = 1./(zls.coef[1]*sqrt(1/(zls.coef[1]*zls.coef[1])+2*pi*2*pi*zls.coef[2]*zls.coef[2]));


wn_extracted = zls.coef[2].*(2*pi) ./ sqrt(1- zeta_extracted[1]*zeta_extracted[1]) ;


den_extracted = [1, 2*zeta_extracted*wn_extracted,  wn_extracted^2 ]; 

X_extracted = damp (den_extracted);
X=damp(lsden2);


colors("yellow");
printf("From Laser Measurements:\n");
colors();

     #printf("\n****************************************************\n");
fr=X_extracted.wn[1]/(2*pi) ;
if( (fr <FMAX) && (fr > FMIN )){
colors("green");
printf("Resonant frequency = %2.1f [Hz]\t\tPASSED.\n",fr);
	 CR[length(CR)+1]=1;
colors();
}


if( (fr >FMAX) || (fr < FMIN )){
colors("red");
printf("Resonant frequency = %2.1f [Hz]\t\tFAILED.\n",fr);
	 CR[length(CR)+1]=0;
colors();
}


frrise = ResRise (1e3*zls.coef[1], fr);


printf("Resonant rise = %2.2f [dB] \n",frrise,,);


travel=5.45*abs(diff(Ldist))[1]; # in [um]
gain = 1e-3*travel/abs(Current[1]-Current[2]);

// TRAVELmin=1200;


colors();
printf("Travel range = %2.1f [um]\n",travel); 

// GAINmin=3.5;
if(  gain < GAINmin ){
colors("red");
printf("Gain = %2.2f [um/mA]\r\t\t\t\t\tFAILED\n",gain);
colors();
}

if(  gain > GAINmin ){
colors("green");
printf("Gain = %2.2f [um/mA]\r\t\t\t\t\tPASSED.\n",gain); colors();
}


// # printf("Elapsed time %3.2f s\n",toc(4));

if((box != "bonnie.grzegorek.com") && (box != "muppet.grzegorek.com")) {
system("scp -q slave:/tmp/slaveresults.r /tmp;");
sleep(1);
 }



rfile slaveresults;

Coil2Rdc = Coil2Rdc - RCAL2;
Coil2Rdc_20deg = Rdc_20degC( Coil2Rdc + RCAL2, T) - RCAL2_20C;

if ( Coil2short_pass=="FAIL" ) {
	 colors("red");
	 printf ("Coil #2 leakage to ground\t\tFAILED.\n");
	 CR[length(CR)+1]=0;
	 colors();
 }


if ( Coil2short_pass=="PASSED" ) {
	 colors("green");
	 printf ("Coil #2 leakage to ground\t\tPASSED.\n");
	 CR[length(CR)+1]=1;
	 colors();
 }


SKEW_ANGLEmin=1.1 ;
if(  SkewAngle < SKEW_ANGLEmin ){
colors("red");
# printf("Skew Angle range = %2.1f [deg]\r\t\t\t\t\tFAILED.\n",SkewAngle);
colors();
}

if(  SkewAngle > SKEW_ANGLEmin ){
colors("green");
# printf("Skew Angle range = %2.1f [deg]\r\t\t\t\t\tPASSED.\n",SkewAngle);
colors();
}

printf("Skew Angle range = %2.1f [deg]\n",SkewAngle);


SKEW_GAINmin=0.92; # [um/mA]
if(  SkewGain < SKEW_GAINmin ){
colors("red");
printf("Skew Gain = %2.2f [um/mA]\r\t\t\t\t\tFAILED\n",SkewGain);
colors();
}

if(  SkewGain > SKEW_GAINmin ){
colors("green");
printf("Skew Gain = %2.2f [um/mA]\r\t\t\t\t\tPASSED.\n",SkewGain); colors();
}


if( (Coil2Rdc_20deg > RDCMAX ) || ( Coil2Rdc_20deg < RDCMIN )){
	 colors("red");
	 printf ("Coil#2 resistance @20°C = %2.2f [Ohm]\tFAILED. Check connections.\n", Coil2Rdc_20deg );
	 CR[length(CR)+1]=0;
	 colors();
 }



if( (Coil2Rdc_20deg < RDCMAX ) && ( Coil2Rdc_20deg > RDCMIN )){
	 colors("green");
	 printf ("Coil#2 resistance @20°C = %2.2f [Ohm]\tPASSED.\n", Coil2Rdc_20deg );
	 CR[length(CR)+1]=1;
	 colors();
 }




////////////////////////// AUG 28 ////////////////////////////////////
DATE=time2dstr(seconds(),"%Y-%m-%d");
TIME=time2dstr(seconds(),"%H:%M:%S");
if( DEBUG2 != 1){op=reads("/var/www/html/Operator.txt");}



fn1="/var/www/html/"+DATE+"_twang.csv";
box=reads("|uname -n "); 
if((box == "bonnie.grzegorek.com") || (box == "muppet.grzegorek.com")) {fn1="/tmp/"+DATE+"_twang.csv";}

box=reads("|uname -n ");
if((box != "bonnie.grzegorek.com") && (box != "muppet.grzegorek.com")) {
# system("cp  /tmp/prod.dat /var/www/html/prod.dat"); sleep(1); 
}





if(isfile(fn1)){open(fn1,"a");}

if(!isfile(fn1)){fprintf(fn1,"%s\n","SN,OPERATOR,DATE,TIME,BEMFfres,BEMFdecay,BEMFresrise,Rdc_Coil#1,Rdc_Coil#1_@20C,Rdc_Coil#2,Rdc_Coil#2_@20C,Coil #1 Short,Coil #2 Short,,LASERfres,LASERdecay,LASERresrise,Gain,travel,Skew Gain,Skew Angle,Tambient");}


LASERfres=fr;
LASERdecay=1e3*zls.coef[1];
LASERresrise=frrise;
Gain=gain;



fprintf(fn1,"%s,%s,%s,%s,%2.2f,%2.2f,%2.2f,%2.2f,%2.2f,%2.2f,%2.2f,%s,%s,,%2.2f,%2.2f,%2.2f,%2.2f,%2.2f,%2.2f,%2.2f,%s\n", SN,op, DATE, TIME, BEMFfres,BEMFdecay,BEMFresrise,Rdc_final,Coil1Rdc_20deg,Coil2Rdc,Coil2Rdc_20deg,Short_pass,Coil2short_pass,LASERfres,LASERdecay,LASERresrise,Gain,travel,SkewGain,SkewAngle,Tamb);


printf("Elapsed time %3.2f s\n",toc(1));
printf("Elapsed time after reading S/N %3.2f s\n",toc(3));

if( prod(CR) == 1 ){
colors("green");system("figlet -k PASS");colors();
}

if( prod(CR) == 0 ){
colors("red");system("figlet -k FAIL");colors();
}

printf("########################################################################\n\n\n");
close(fn1);
if(( DEBUG2 != 1) && ( LOOP != 200000 )){sleep(10);}

if(m == 1){clear(SN);
if(isfile("/tmp/SN.txt")){system("/bin/rm /tmp/SN.txt");}
}
if( DEBUG2 != 1){ clear(d);}

}

