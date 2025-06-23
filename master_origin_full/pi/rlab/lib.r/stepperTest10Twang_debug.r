#!/home/pi/bin/rlab2

require savepdf coarse_motion_linearity extract_positions runavg subsample blink_until_pressed
///////////////// NEW ////////////////////////////////////////////////
// while(1) {

// DEBUG2=1; // Testing at home without lasers 
DEBUG2=0;    // Normal opreation
if(DEBUG2 != 1){
//////////////////////////////////////////////////////////////////////

// Open 3 ports (Laser controller + stepper controller + barcode reader)
# for (i in 1:5){
# plclose(); // PGPLOT dies if window not closed ?
tic(1);
tic(2);
tic(3);

// Laser controller //////////////////////////////////////////////////
attr=<<>>;
attr.devname = "ttyUSB";
attr.id_vendor_id = "0403";

if(!exist(r1)){
r1 = usbdevname(attr);
}

#r1
# toc(1)

opts=<<>>;
opts.data_parity_stop="8N1";
opts.eol="\r";
# opts.debug="char";
opts.speed=115200;
opts.flow_control="n";

serial="serial://" + r1;

if(!open(serial,opts)) {
open(serial, opts);
}

# toc(2)

readdist = ["MS,02\r"];
# readcode = ["MS,01,+000.260"];
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
#      printf("%s\n",r);
      // verify if out is laser2
 writem(serial, analogchan12);
      r = readm(serial);
#      printf("%s\n",r);
# writem(serial, analogchan22);
 #     r = readm(serial);
  #    printf("%s\n",r);
      //switch  back to measuremnt mode
 writem(serial, generalmode);
      r = readm(serial);
#      printf("%s\n",r);
      
      
//////////////////////////////////////////////////////////////////////

// stepper controller ////////////////////////////////////////////////

attr=<<>>;
attr.devname = "ttyUSB";
attr.id_vendor_id = "10c4";

if(!exist(r2)){
r2 = usbdevname(attr);
}

#r2
#toc(1)
opts=<<>>;
opts.data_parity_stop="8N1";
opts.eol="\r";
# opts.debug="char";
opts.speed=9600;
opts.flow_control="n";

stepper="serial://" + r2;

if(!open(stepper,opts)) {
open(stepper, opts);
}




// opto = ["/1Z3000000R\r"];
// # move up 200
// # NOTE:  
// newhome = ["/1V200P400R\r"];
// # move up 672.5
// #step1 = ["/1P672.5R\r"];
// step1 = ["/1V200P1345R\r"];
// # move down 740
// step2 =["/1V200D1480R\r"];
// # move up 740
// step3 =["/1V200P1480R\r"];
// # move down 740
// step4 =["/1V200D1480R\r"];
// gohome = ["/1V200A400R\r"];
// # move up 63
// end =["/1V200P126R\r"];
// vel =["/1V200R\r"];
// VEL=num2str(400);
# max ping 65000
total=["/1v400Z3000000z200M1000V400P1345M1000D1480M1000P1480M1000D1480M1000v400Z300000M1000V400P126R"];

///////////// Barcode reader /////////////////////////////////////////
attr=<<>>;
attr.devname = "ttyACM";
attr.id_vendor_id = "1eab";
# attr.id_vendor_id = "0483";
# attr.id_vendor_id = "0403";
### [azg]seems to work without S/N
### attr.id_serial_short = "20383959424B";
# attr.id_serial_short = "A600CJ5X";
 
 
 
if(!exist(r3)){
r3 = usbdevname(attr);
}
 
# remove below } if not DEBUG2
} 

if( DEBUG2 == 1) { SN="blablabla"; }



# system("cat /dev/ttyACM0 >/dev/null"); // to be fixed

 
if(!exist(SN)){
colors("red");
printf("Waiting for S/N...");
 
r=getline(r3,16);
# remove trailing \n
SN=sum (strsplt(r)[1:length(strsplt(r)) -1]);
} 
colors("green");
# printf("\rS/N:              \t%s",r);
printf("\rSN:              \t%s\n",SN);
# printf("\rS/N:              \t%s\n",r);
colors();

////////////
 gmt= date2jd(gmtime());
 gmt= gmt -floor(gmt);
 gmt=100000*gmt;
 sprintf(dd,"%i.dat",gmt);
//////////////////////////////////////////////////////////////////////
printf("COARSE MOTION TEST:\n");
if(DEBUG2 != 1){
blink_until_pressed();
}


# pause("Press enter when ready to start coarse motion test");
//////////////////////////////////////////////////////////////////////

#system("sudo /home/pi/sandbox/EZstepper17/Jul25_2018/coarseMCP3202.short > /tmp/datcoarse.dat&");
#system("sudo /home/pi/sandbox/EZstepper17/Jul25_2018/coarse_25s_MCP3202 > /tmp/datcoarse.dat&");

if(DEBUG2 != 1){
system("sudo /home/pi/sandbox/EZstepper17/Jul25_2018/coarse_30s_MCP3202 > /tmp/datcoarse.dat&");
}

if(DEBUG2 == 1){
system("sudo /home/pi/sandbox/MCP3202_with_pigpio/coarse_30s_MCP3202 > /tmp/datcoarse.dat&");
}

# [azg] Thu Oct 18 20:39:11 PDT 2018
# sometimes ADC does not capture the beginning. Added sleep() below
sleep(0.5);
	# DEBUG only
# system("sudo /home/pi/sandbox/EZstepper17/Jul25_2018/coarseMCP3202.short > /home/pi/sandbox/EZstepper17/Aug_22_2018/"+SN+"_"+dd+"coarse.dat&");

printf("\n\n%s\n","**********************Stepper test starts NOW.");


if(DEBUG2 != 1){

 writem(stepper, total);
      r = readm(stepper);
      #printf("%s\n",r);


 r = readm(stepper);
      #printf("%s\n",r);

######################################################################
////////////////////// RS232 plot ////////////////////////////////////
######################################################################

# toc(3)




// switch back analog output to laser 1
// sleep(10);
// switch back laser 1 to analog output 1
# sleep(4)
#sleep(45);
# sleep(32);
sleep(35);  # Thu Oct 18 21:23:52 PDT 2018
# sleep(10);
 
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
# Fri Oct 19 21:47:51 PDT 2018
# sleep(10);
######################################################################
/////////////////////////// AUG 27 ///////////////////////////////////
	printf("reading measured laser data\n");
# coarse_motion_linearity("MAT1617CS01_70778.datcoarse.dat",+1,1,"MAT1617CS01")
# [azg] Thu Oct 18 19:44:35 PDT 2018

if( DEBUG2 == 1) { system("cp datcoarse.dat /tmp/datcoarse.dat;"); sleep 1; }

xxx=read_ascii("/tmp/datcoarse.dat");
# flush the first row
d=d[2:d.nr;];



# P1=coarse_motion_linearity("/tmp/datcoarse.dat",+1,1,SN);
P1=coarse_motion_linearity(d,+1,1,SN);
# sleep(1);
# plclose();
# P2=coarse_motion_linearity("/tmp/datcoarse.dat",-1,1,SN);
P2=coarse_motion_linearity(d,-1,1,SN);
# sleep(1);
# plclose();
savepdf(P1.plt1, P1.plt2, SN, 1);
savepdf(P2.plt1, P2.plt2, SN, -1);

# E=extract_positions("/tmp/datcoarse.dat");
E=extract_positions(d);

STEP1= E.STEP1err;
STEP2= E.STEP2err;
STEP3= E.STEP3err;
STEP4= E.STEP4err;


////////////////////////// AUG 27 ////////////////////////////////////
// add writing to a file:w
DATE=time2dstr(seconds(),"%Y-%m-%d");
TIME=time2dstr(seconds(),"%H:%M:%S");
op=reads("/var/www/html/Operator.txt");


# fn1="/home/pi/results/results_coarse_motion.csv";
fn1="/var/www/html/"+DATE+"_coarse_motion.csv";
# DEBUG
# fn1="junk.csv";op="me";


# if(!isfile(fn1)){open(fn1,"a");fprintf(fn1,"%s\n","SN,OPERATOR,DATE,TIME,STEP1,STEP2,STEP3,STEP4");}
if(!isfile(fn1)){fprintf(fn1,"%s\n","SN,OPERATOR,DATE,TIME,STEP1,STEP2,STEP3,STEP4");}

if(isfile(fn1)){open(fn1,"a");}

fprintf(fn1,"%s,%s,%s,%s,%2.2f,%2.2f,%2.2f,%2.2f\n", SN,op,DATE,TIME,STEP1,STEP2,STEP3,STEP4);

// fprintf(fn1,"%s,%s,%s,%s,,\n", SN,op,DATE,TIME);
// fprintf(fn1,"%s,%s,%s,%s,%2.2f,%2.2f,%2.2f,%2.2f,%2.2f,%2.2f,%2.2f,\n", SN,op,DATE,TIME);
// close(fn1);
//////////////////////////////////////////////////////////////////////
# }
######################################################################
printf("Continue with Twang test\n");
if(DEBUG2 != 1){
blink_until_pressed();
}

######################################################################
system("sudo /home/pi/sandbox/MCP3202_with_pigpio/prodMCP3202 > /tmp/prod.dat");

sleep(4); // make sure data is collected before processing starts

################# TWANG STARTS HERE ##################


plclose(); // PGPLOT dies if window not closed ?

#!/home/pi/bin/rlab2
require rem expm lsim lsim2 tfeval ResRise ap5 damp plot subsample dfdp extract_chunk psprint 
//////////////////////////////////////////////////////////////////////
##  attr=<<>>;
##  attr.devname = "ttyACM";
##  attr.id_vendor_id = "1eab";
##  # attr.id_vendor_id = "0483";
##  # attr.id_vendor_id = "0403";
##  ### [azg]seems to work without S/N
##  ### attr.id_serial_short = "20383959424B";
##  # attr.id_serial_short = "A600CJ5X";
 
 
 
##  if(!exist(r3)){
	##  r3 = usbdevname(attr);
##  }
 
 
# system("cat /dev/ttyACM0 >/dev/null"); // to be fixed

 
## #  if(!exist(SN)){
## if(isempty(SN)){
	## colors("red");
	## printf("Waiting for S/N...");
	 ## 
	## r=getline(r3,16);
	## # remove trailing \n
	## SN=sum (strsplt(r)[1:length(strsplt(r)) -1]);
## } 
## colors("green");
## # printf("\rS/N:              \t%s",r);
## printf("\rSN:              \t%s\n",SN);
## # printf("\rS/N:              \t%s\n",r);
## colors();
//////////////////////////////////////////////////////////////////////

# blink_until_pressed();


# pause("Press enter when ready to start Twang test");

//////////////////////////////////////////////////////////////////////
# "line 2"
tic(4);
# tic();xxx=read_ascii("no_timestamp_production_vol=0.2.dat");
# xxx=read_ascii("/tmp/prod.dat");
if( DEBUG2 == 1) { system("cp OK.prod.dat /tmp/prod.dat;"); sleep 1; }
dir="/tmp/";
file1="prod.dat";
# xxx=read_ascii("/tmp/prodJune18_1.dat");
# dir="/home/azg/sandbox/usbreader/June26_2018/";
# file1="Original_with_reversed_wires_prod.dat";
xxx=read_ascii(dir + file1);

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


 # plot([4e-5*[1:length(d)]', [d;1]]);
 ##plot([4e-5*[ind[ind2[1]]:ind[ind2[4]]+7000]', [d[ind[ind2[1]]:ind[ind2[4]]+7000 ;1:2]]]);

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

# wire resistance is 0.27 ohm
# Rdc_final=sum(Rx)/length(Rx) - 0.27;
# new value for improved PCB
# Wed Aug 29 20:18:25 PDT 2018
# Rdc_final=sum(Rx)/length(Rx) - 0.37;
# new value using 1% 13.0 ohm resistor
# Tue Oct 16 21:02:11 PDT 2018
# Rdc_final=sum(Rx)/length(Rx) - 0.63;
# Mon Oct 22 20:12:01 PDT 2018
Rdc_final=sum(Rx)/length(Rx) - 0.69;

# DEBUG
# Rdc_final =  Rdc_final +1;


if( (Rdc_final > 14.0) || ( Rdc_final < 12.0)){
	 colors("red");
	 printf ("Coil resistance = %2.2f [Ohm]\t\tFAILED. Check connections.\n", Rdc_final );
	 colors();
 }



if( (Rdc_final < 14.0) && ( Rdc_final > 12.0)){
	 colors("green");
	 printf ("Coil resistance = %2.2f [Ohm]\t\tPASSED.\n", Rdc_final);
	 colors();
 }



# Leakage test

if (abs(sum(Ip-Im)) > 0.1* abs(sum(Im))) {
# if ((sum(Ip-Im)) > 0.1* sum(Im)) {
# DEBUG failure
# if ((sum(Ip-Im)) < 0.1* sum(Im)) {
 
	 colors("red");
	 printf ("Leakage test \t\t\t\tFAILED. Coil shorted to ground.\n"); colors();
Short_pass="FAIL";
 }


if ((sum(Ip-Im)) < 0.1* sum(Im)) {
	 colors("green");
	 printf ("Coil leakage to ground\t\t\tPASSED.\n");
	 colors();
Short_pass="PASSED";
 }


# 	 printf ("Color check\n");
######################################################################

# " BEMF begin at:"
BEMF_extract_start = ind[find((diff(ind) >2500) && (diff(ind) <3000 ))];
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
k = max(size(BEMF_extract_start))-3 ;

if(DEBUG == 1){printf("BEMF k= %i\n",k); }

plot ([4e-5*[BEMF_extract_start[k]+skip : BEMF_extract_start[k]+1000 ]', [d[BEMF_extract_start[k]+skip : BEMF_extract_start[k]+1000;2]]]);

dc = sum(d[BEMF_extract_start[k]+500 : BEMF_extract_start[k]+1000;2]) / length(d[BEMF_extract_start[k]+500 : BEMF_extract_start[k]+1000;2]);
x = 4e-5*[BEMF_extract_start[k]+skip : BEMF_extract_start[k]+1000 ]';
x=x - x[1]; # reset time to zero
y = d[BEMF_extract_start[k]+skip : BEMF_extract_start[k]+1000;2] -dc;

############################################################################
# [azg] Wed Apr 25 00:39:39 PDT 2018
# ss=5;
# y=subsample(y,ss);
# x=subsample(x,ss);

############################################################################
func = function(x,p){
	global(u2,pi)

    return p[4] + p[3]*exp(-(x)/p[1]).*sin(2*pi*p[2].*(x-p[5]));
    // return p[4] + p[3]*exp(-(x)/p[1]).*sin(2*pi*p[2].*(x-p[5]));
    };

############################################################################




############################################################################

// nonlinear fit
// initial guess
#   T   f   Amp
# impulse only
p=[10e-3,100,200,0,0.2e-3]';
p=[10e-3,100,200,0,0.2e-3]';
p=[3e-3,100,-200,0,0.2e-3]';
p=[3e-3,50,-800,0,0.2e-3]';
p=[3e-3,100,800,0,0.2e-3]';
p=[3e-3,100,-800,0,1.2e-3]';
p=[0.00273891262,    135.484423,    88.5418362,  0.0607532924,  -0.00165477464 ]';
p=[0.00498460066,    96.9925509,    56.0774607,   0.737717201,  -0.00214367662]';
p=[0.00383467736,    81.2222657,    489.100125,   0.888078179,  0.000664467381 ]';
p=[0.0065547235,     62.114712,   -74.7312899,    2.28336221,  0.00327974194 ]';
p=[0.00786476649,    101.730612,   -100.142629,   0.399497618,  0.00231118448]';
p=[0.00495459908,    112.671236,   -96.8946583,  -0.803299967,  0.00234807956 ]';
p=[0.00495459908,    112.671236,   -96.8946583,  -0.803299967,  0.00234807956  ]';
p=[ 0.00218345898,    137.092441,   -64.5156475,   0.343337086,  0.00280853316  ]';
p=[0.00383467736,    81.2222657,    489.100125,   0.888078179,  0.000664467381 ]';
p=[0.00459778447,    95.4361057,    163.658406,   0.241976302,  -0.00181022603  ]';
p=[0.00461433261,    86.2906227,    170.583871,   0.736565834,  -0.00240747643 ]';
p=[0.00450409989,    85.3922691,    168.112407,   0.709823261,  -0.00244084615  ]';
p=[0.00407598237,    85.1165973,     178.45966,   0.558672074,  -0.00256723025  ]';
p=[ 0.00332153576,    87.3182479,    289.401213,   0.450137308,  -0.00117697422  ]';



xy = [x,y];

options=<<>>;
options.imethod = 2;  // least squares fit
# options.stdout = rconsole();

 zls=lsfit(y,x,p,func,dfdp);
# "BEMF zls.coef =" zls.coef


     #printf("\n****************************************************\n");
#if(1e3*zls.coef[1] <= 4.4 && 1e3*zls.coef[1] > 2.76 ){
	#printf("T decay = %2.2f [ms] \n",1e3*zls.coef[1]);
#}

#if(1e3*zls.coef[1] > 4.4 || 1e3*zls.coef[1] < 2.76 ){
	#colors("red");
	#printf("T decay = %2.2f [ms] \n",1e3*zls.coef[1],,);
	#colors();
#}



#printf("F res   = %f  [Hz]\n",zls.coef[2]);

pltitle("Fitting measured data ");
xlabel("Time [s]");
ylabel("Amplitude");
plegend(["Data","Lsfit"]);
##################################################
# if (_rlab_config.plot_support=="pgplot"||_rlab_config.plot_support=="plplot")
if (_rlab_config.plot_support=="pgplot")
	    {
               plaxis("lin",); # pgplot syntax
     else
              plscale("lin",); # plplot syntax
   }

##################################################
      # plaxis("lin",); # pgplot syntax
      # plscale("lin",); # plplot syntax
      # plaxis("lin","lin");
# " Debug 116"
zzls=(func(x,zls.coef));


# if(DEBUG2 == 1){printf("Line 624\n"); }
plot(<<[x,y];[x,zzls]>> );
#pause();
sleep(3);
# if(DEBUG2 == 1){printf("Line 627\n"); }

# _plprint ("time_data.eps", "psc");

# " Debug 126"

# from lsfit
lsnum2=[1,2*pi*zls.coef[2]];
lsden2=[1,2/zls.coef[1],1/(zls.coef[1]*zls.coef[1])+2*pi*2*pi*(zls.coef[2]*zls.coef[2])]; 

# if(DEBUG2 == 1){printf("Line 637\n"); }
# Need the same form as lsim2 input ( add 1/R factor to coef[1] 
# This is transfer function of back EMF
lsnum2=[zls.coef[1],0];
lsden2=[1,2/zls.coef[1],1/(zls.coef[1]*zls.coef[1])+2*pi*2*pi*(zls.coef[2]*zls.coef[2])]; 

# if(DEBUG2 == 1){printf("Line 643\n"); }
# the displacement is proportional to current
# since V = dI/dt * L
# current is Back EMF(s) * 1/sL

dnum2=[zls.coef[1]];
dden2= lsden2;

freq = [10:3000:10];
# input tf
#inum2 = num2[1];
inum2 = [1, 0];
#itf = tfeval(inum2, den2, freq);
itf = tfeval(inum2, lsden2, freq);
idb=20*log10(abs(itf));



# if(DEBUG2 == 1){printf("Line 661\n"); }
xlabel("Frequency [Hz]");
if (_rlab_config.plot_support=="pgplot")
	    {
               plaxis("log",); # pgplot syntax
     else
              plscale("log",); # plplot syntax
   }

# if(DEBUG2 == 1){printf("Line 670\n"); }
ylabel("Amplitude [dB]");
# plscale("log",); # plplot syntax
# plaxis("log",);
plegend();
pltitle("Frequency response");
# plegend(["Odr Fit","Lsfit","Displacement"]);
plegend(["Back EMF","Displacement"]);
plegend(["Back EMF", "Input data","Displacement"]);
plegend(["Back EMF", "Input data"]);

#** tf = tfeval(num2, den2, freq);
#** db=20*log10(abs(tf));
#** lstf = tfeval(lsnum2, lsden2, freq);
# lsdb=20*log10(abs(lstf));
#** lsdb=20*log10(num2[1]/lsnum2[1]) + 20*log10(abs(lstf));

# ******** # plot(<<[freq;db]';[freq;lsdb]'>>);
# ******** # plot(<<[freq;lsdb]'>>);
# sleep(3);

# if(DEBUG2 == 1){printf("Line 691\n"); }
dtf = tfeval(dnum2, dden2, freq);
dtf0 = tfeval(dnum2, dden2, freq[1]);

# WHy this worked before? [azg] Sat Apr 21 17:32:09 PDT 2018
# dispdb=20*log10(inum2/dnum2) + 20*log10(abs(dtf));
# dispdb= 20*log10(abs(dtf));
dispdb = 20*log10(abs(dtf)) - 20*log10(abs(dtf0));

xlabel("Frequency [Hz]");
ylabel("Amplitude [dB]");
pltitle("Displacement Frequency Response");
# plaxis("log",);
if (_rlab_config.plot_support=="pgplot")
	    {
               plaxis("log",); # pgplot syntax
     else
              plscale("log",); # plplot syntax
   }

   
   # plscale("log",); # plplot syntax
plegend(["Input data","Displacement"]);
plegend(["Displacement"]);

# if(DEBUG2 == 1){printf("Line 716\n"); }
# " Debug 199"
plot(<< [freq;dispdb][;1:30]'>>);
sleep(3);
# plot(<<[freq;lsdb]';[freq;dispdb]'>>);
# plot(<<[freq;dispdb]'>>);
# sleep(2);


# _plprint ("frequencyResponse.eps", "psc");



########## Extracte natural frew #########
#  zls.coef[2]) = wn./(2*pi) .* sqrt(1- zeta[1]*zeta[1])ma
#  2*pi*zls.coef[2])./  sqrt(1- zeta[1]*zeta[1]) = wn
# get the Q ffrom 2 consecutive peaks since tau and f is known
# tbc...

######################
# formula sanity check
# 
# from displacement transfer function (at resonanse)
# dispdb[10] - dispdb[1]=  -0.0346820712 
#
# from ResRise.r Tau_decay.rj
# ResRise(3.180000,100) = -0.00846014643  
# Tau_decay(-0.0346820712,100) = 3.17041434  

# from simulation fit:
# T decay = 3.180000 [ms] 
# lsfit F res   = 86.686933  [Hz]
#
######################
# zeta_extracted = 1./(2*pi*fringing*tau) # <- this is wrong
# zeta_extracted = 1./(2*pi*fres*tau)     # this is consistent with other formulas
# zeta_extracted = 1./(2*pi*zls.coef[1]*zls.coef[2]);  # <- this is wrong

zeta_extracted = 1./(zls.coef[1]*sqrt(1/(zls.coef[1]*zls.coef[1])+2*pi*2*pi*zls.coef[2]*zls.coef[2]));

# printf("zeta_extracted = %f\n",zeta_extracted);

wn_extracted = zls.coef[2].*(2*pi) ./ sqrt(1- zeta_extracted[1]*zeta_extracted[1]) ;
# printf("wn_extracted = %f\n",wn_extracted);

# ap5([freq;db]');


den_extracted = [1, 2*zeta_extracted*wn_extracted,  wn_extracted^2 ]; 

X_extracted = damp (den_extracted);
X=damp(lsden2);

# NOTE: we  print only first number from every member (eg. Poles have many numbers)
# printf("%s\n", "Input data");
# for(k in members(X)){printf("%s\r\t\t = %f\n",k,X.[k]);}

# printf("%s\n", "Extracted data");
# for(k in members(X_extracted)){printf("%s\r\t\t = %f\n",k,X_extracted.[k]);}

#printf("\n****************************************************\n");
fr=X_extracted.wn[1]/(2*pi) ;
#if( (fr <120) && (fr > 80 )){
#printf("Resonant frequency      = %2.2f [Hz]\n",fr);
#}
#if( (fr >120) || (fr < 80 )){
#colors("red");
##printf("Resonant frequency      = %2.2f [Hz]\n",fr);
#colors();
#}
# printf("Resonant frequency      = %f [Hz]\n",X_extracted.wn/(2*pi));
# printf("Frequency response peak = %2.2f [Hz]\n",X_extracted.Fringing);
printf("Frequency Response Peak = %2.2f [dB]\n",X_extracted.FRpeaking);
# printf("Coil current             = %f [A]\n",Isteady);
# printf("chunk number             = %i \n",chunk);
#"zls.coef"
#zls.coef
#printf("\n****************************************************\n");
# decay time constant
# DT[m;]=[Isteady*1e3, 1e3*zls.coef[1]]; 
# reject questionable fits
##if( (fr <120) && (fr > 80 )){ DT[m;]=[Isteady*1e3, 1e3*zls.coef[1]]; }
##m=m+1;
# printf("\n\n\t\t%s\t\t%s\t\t%s\n","   Input data","   Extracted","Error %");
# for(k in members(X)){ printf("%s\r\t\t = %f",k,X.[k]); printf("\t\t = %f\t\t%f\n",X_extracted.[k],(X.[k]-X_extracted.[k])/X.[k]*100);}
# toc()
# }

# TAUmin=2.76;
# TAUmax=4.4;
TAUmin=2.3;
TAUmax=5.9;


     #printf("\n****************************************************\n");
if(1e3*zls.coef[1] <= TAUmax && 1e3*zls.coef[1] > TAUmin ){
colors("green");
printf("T decay = %2.2f [ms] \t\t\tPASSED.\n",1e3*zls.coef[1]);
colors();
}

if(1e3*zls.coef[1] > TAUmax || 1e3*zls.coef[1] < TAUmin ){
colors("red");
printf("T decay = %2.2f [ms] \t\t\tFAILED.\n",1e3*zls.coef[1],,);
colors();
}


FMAX=120;
FMIN=80;

RESRISEmin = -1.2;
RESRISEmax = 2.8;

     #printf("\n****************************************************\n");
fr=X_extracted.wn[1]/(2*pi) ;
if( (fr <FMAX) && (fr > FMIN )){
colors("green");
printf("Resonant frequency = %2.1f [Hz]\t\tPASSED.\n",fr);
colors();
}


if( (fr >FMAX) || (fr < FMIN )){
colors("red");
printf("Resonant frequency = %2.1f [Hz]\t\tFAILED.\n",fr);
colors();
}



# frrise = ResRise (Tau_decay, fres);
# frrise = ResRise (1e3*z.coef[1], z.coef[2]);
# frrise = ResRise (1e3*zls.coef[1], zls.coef[2]);
frrise = ResRise (1e3*zls.coef[1], fr);
##printf("\n****************************************************\n");
# printf("Calculated Fres rise using IBM formula  = %f [dB]\n",frrise);
# res rise 0.8 =/- 2 dB
# +2.8 -1.2 dB
if(frrise <= RESRISEmax && frrise > RESRISEmin ){

colors("yellow");
printf("From BackEMF Measuremets:\n");
colors("green");
printf("Resonant rise = %2.2f [dB] \t\tPASSED.\n",frrise);
colors();
}

if(frrise > RESRISEmax || frrise < RESRISEmin ){
colors("red");
printf("Resonant rise = %2.2f [dB] \t\tFAILED.\n",frrise,,);
colors();
}
#printf("\n****************************************************\n");

printf("Elapsed time %3.2f s\n",toc(4));

BEMFfres=fr;
BEMFdecay=1e3*zls.coef[1];
BEMFresrise=frrise;

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%% laser begin here

######################################################################
plclose();
# " Lsr begin at:"
# Lsr_extract_start = ind[find((diff(ind) >2500) && (diff(ind) <3000 ))];

# Experimental data
# Conversion ratio of um to ADC count:
#                                      1200/(max(d[;5])-min(d[;5]))
um = 5.17;

skip=120;
# skip=75;
	                 
# k=6;
# [azg] Fri Oct 19 21:11:17 PDT 2018
# use the same k as BEMF
# k=7;
if(DEBUG == 1){printf("Laser k= %i\n",k); }

pltitle("Fitting measured laser data ");
xlabel("Time [s]");
ylabel("Amplitude [um]");
#plaxis("lin",); # pgplot syntax
if (_rlab_config.plot_support=="pgplot")
	    {
               plaxis("lin",); # pgplot syntax
     else
              plscale("lin",); # plplot syntax
   }
# plscale("lin",); # plplot syntax
plot ([4e-5*[Lsr_extract_start[k]+skip : Lsr_extract_start[k]+1000 ]', [um*d[Lsr_extract_start[k]+skip : Lsr_extract_start[k]+1000;5]]]);

dc = sum(d[Lsr_extract_start[k]+500 : Lsr_extract_start[k]+1000;5]) / length(d[Lsr_extract_start[k]+500 : Lsr_extract_start[k]+1000;5]);
x = 4e-5*[Lsr_extract_start[k]+skip : Lsr_extract_start[k]+1000 ]';
x=x - x[1]; # reset time to zero
y = d[Lsr_extract_start[k]+skip : Lsr_extract_start[k]+1000;5] -dc;


############################################################################
# func = function(x,p){
	# global(u2,pi)
# 
    # return p[4] + p[3]*exp(-(x)/p[1]).*sin(2*pi*p[2].*(x-p[5]));
    # };

############################################################################

############################################################################

// nonlinear fit
// initial guess
#   T   f   Amp
# impulse only
p=[ 0.00327633722,    84.8807812,    133.730521,   0.522658955,  -0.00134066124 ]'; # for laser data 
p=[ 0.00296779553,    86.1768055,    132.631621,   0.502642791,  -0.00124655545 ]'; # for laser data 


xy = [x,y];

options=<<>>;
options.imethod = 2;  // least squares fit
# options.stdout = rconsole();

 zls=lsfit(y,x,p,func,dfdp);
# zls.coef



pltitle("Fitting laser measured data ");
xlabel("Time [s]");
ylabel("Amplitude [um]");
plegend(["Data","Lsfit"]);

if (_rlab_config.plot_support=="pgplot")
	    {
               plaxis("lin",); # pgplot syntax
     else
              plscale("lin",); # plplot syntax
   }

#plscale("lin",); # plplot syntax
      #plaxis("lin",); # pgplot syntax
      # plaxis("lin","lin");
# " Debug 479"
zzls=(func(x,zls.coef));


plot(<<[x,um*y];[x,um*zzls]>> );
#pause();
sleep(3);

# _plprint ("time_data.eps", "psc");

# " Debug 489"

# from lsfit
lsnum2=[1,2*pi*zls.coef[2]];
lsden2=[1,2/zls.coef[1],1/(zls.coef[1]*zls.coef[1])+2*pi*2*pi*(zls.coef[2]*zls.coef[2])]; 

# Need the same form as lsim2 input ( add 1/R factor to coef[1] 
# This is transfer function of back EMF
lsnum2=[zls.coef[1],0];
lsden2=[1,2/zls.coef[1],1/(zls.coef[1]*zls.coef[1])+2*pi*2*pi*(zls.coef[2]*zls.coef[2])]; 

# the displacement is proportional to current
# since V = dI/dt * L
# current is Back EMF(s) * 1/sL

dnum2=[zls.coef[1]];
dden2= lsden2;

freq = [10:3000:10];
# input tf
#inum2 = num2[1];
inum2 = [1, 0];
#itf = tfeval(inum2, den2, freq);
itf = tfeval(inum2, lsden2, freq);
idb=20*log10(abs(itf));



xlabel("Frequency [Hz]");
ylabel("Amplitude [dB]");
if (_rlab_config.plot_support=="pgplot")
	    {
               plaxis("log",); # pgplot syntax
     else
              plscale("log",); # plplot syntax
   }

# plscale("log",); # plplot syntax
# plaxis("log",);
plegend();
pltitle("Frequency response");


dtf = tfeval(dnum2, dden2, freq);
dtf0 = tfeval(dnum2, dden2, freq[1]);

dispdb = 20*log10(abs(dtf)) - 20*log10(abs(dtf0));

xlabel("Frequency [Hz]");
ylabel("Amplitude [dB]");
pltitle("Displacement Frequency Response");
# plaxis("log",);
# plscale("log",); # plplot syntax

if (_rlab_config.plot_support=="pgplot")
	    {
               plaxis("log",); # pgplot syntax
     else
              plscale("log",); # plplot syntax
   }

plegend(["Input data","Displacement"]);
plegend(["Displacement"]);

# " Debug 199"
plot(<< [freq;dispdb][;1:30]'>>);
sleep(2);
plclose();
# plot(<<[freq;lsdb]';[freq;dispdb]'>>);
# plot(<<[freq;dispdb]'>>);
# sleep(2);


# _plprint ("frequencyResponse.eps", "psc");


zeta_extracted = 1./(zls.coef[1]*sqrt(1/(zls.coef[1]*zls.coef[1])+2*pi*2*pi*zls.coef[2]*zls.coef[2]));

# printf("zeta_extracted = %f\n",zeta_extracted);

wn_extracted = zls.coef[2].*(2*pi) ./ sqrt(1- zeta_extracted[1]*zeta_extracted[1]) ;
# printf("wn_extracted = %f\n",wn_extracted);


den_extracted = [1, 2*zeta_extracted*wn_extracted,  wn_extracted^2 ]; 

X_extracted = damp (den_extracted);
X=damp(lsden2);

     #printf("\n****************************************************\n");
if(1e3*zls.coef[1] <= TAUmax && 1e3*zls.coef[1] > TAUmin ){
colors("green");
printf("T decay = %2.2f [ms] \t\t\tPASSED.\n",1e3*zls.coef[1]);
colors();
}

if(1e3*zls.coef[1] > TAUmax || 1e3*zls.coef[1] < TAUmin ){
colors("red");
printf("T decay = %2.2f [ms] \t\t\tFAILED.\n",1e3*zls.coef[1],,);
colors();
}




     #printf("\n****************************************************\n");
fr=X_extracted.wn[1]/(2*pi) ;
if( (fr <FMAX) && (fr > FMIN )){
colors("green");
printf("Resonant frequency = %2.1f [Hz]\t\tPASSED.\n",fr);
colors();
}


if( (fr >FMAX) || (fr < FMIN )){
colors("red");
printf("Resonant frequency = %2.1f [Hz]\t\tFAILED.\n",fr);
colors();
}


#frrise = ResRise (1e3*zls.coef[1], zls.coef[2]);
frrise = ResRise (1e3*zls.coef[1], fr);
##printf("\n****************************************************\n");
# printf("Calculated Fres rise using IBM formula  = %f [dB]\n",frrise);
# res rise 0.8 =/- 2 dB
# +2.8 -1.2 dB
if(frrise <= RESRISEmax && frrise > RESRISEmin ){
colors("yellow");
printf("From Laser Measurements:\n");
colors("green");
printf("Resonant rise = %2.2f [dB] \t\tPASSED.\n",frrise);
colors();
}

if(frrise > RESRISEmax || frrise < RESRISEmin ){
colors("red");
printf("Resonant rise = %2.2f [dB] \t\tFAILED.\n",frrise,,);
colors();
}
#printf("\n****************************************************\n");

# Gain
travel=5.45*abs(diff(Ldist))[1]; # in [um]
gain = 1e-3*travel/abs(Current[1]);

if(  travel < 1200 ){
colors("red");
printf("Travel range = %2.1f [um]\r\t\t\t\t\tFAILED\n",travel);
colors();
}

if(  travel  > 1200 ){
colors("green");
printf("Travel range = %2.1f [um]\r\t\t\t\t\tPASSED.\n",travel); colors();
}

if(  gain < 3.5 ){
colors("red");
printf("Gain = %2.2f [um/mA]\r\t\t\t\t\tFAILED\n",gain);
colors();
}

if(  gain > 3.5 ){
colors("green");
printf("Gain = %2.2f [um/mA]\r\t\t\t\t\tPASSED.\n",gain); colors();
}


printf("Elapsed time %3.2f s\n",toc(4));
printf("########################################################################\n\n\n");

# fn1="/tmp/junk.csv";
# fn1="/home/pi/results/results.csv";
////////////////////////// AUG 28 ////////////////////////////////////
	// add writing to a file:w
DATE=time2dstr(seconds(),"%Y-%m-%d");
TIME=time2dstr(seconds(),"%H:%M:%S");
op=reads("/var/www/html/Operator.txt");
# HACK to be removed
# SN=reads("/var/www/html/SN.txt");


# fn1="/home/pi/results/results_coarse_motion.csv";
fn1="/var/www/html/"+DATE+"_twang.csv";
# DEBUG

////////////////////////// AUG 28 ////////////////////////////////////



# if(!exist(fn1)){open(fn1,"a");fprintf(fn1,"%s\n","SN,,93.5698,3.96085,1.32148,13.3,,BEMFfres,BEMFdecay,BEMFresrise,Rdc_final,Short_pass,,LASERfres,LASERdecay,LASERresrise,Gain,travel");}

# if(!exist(fn1)){open(fn1,"a");fprintf(fn1,"%s\n","SN,,,,,,,BEMFfres,BEMFdecay,BEMFresrise,Rdc_final,Short_pass,,LASERfres,LASERdecay,LASERresrise,Gain,travel");}
# Tue Aug 28 21:13:44 PDT 2018

if(isfile(fn1)){open(fn1,"a");}

if(!isfile(fn1)){fprintf(fn1,"%s\n","SN,OPERATOR,DATE,TIME,BEMFfres,BEMFdecay,BEMFresrise,Rdc_final,Short_pass,,LASERfres,LASERdecay,LASERresrise,Gain,travel");}


LASERfres=fr;
LASERdecay=1e3*zls.coef[1];
LASERresrise=frrise;
Gain=gain;

# fprintf(fn1,"%s,,93.5698,3.96085,1.32148,13.3,,%2.2f,%2.2f,%2.2f,%2.2f,%s,,%2.2f,%2.2f,%2.2f,%2.2f,%2.2f,\n", SN, BEMFfres,BEMFdecay,BEMFresrise,Rdc_final,Short_pass,LASERfres,LASERdecay,LASERresrise,Gain,travel);

# fprintf(fn1,"%s,%s,%2.2f,%2.2f,%2.2f,%2.2f,%s,,%2.2f,%2.2f,%2.2f,%2.2f,%2.2f,\n", SN,op, BEMFfres,BEMFdecay,BEMFresrise,Rdc_final,Short_pass,LASERfres,LASERdecay,LASERresrise,Gain,travel);
fprintf(fn1,"%s,%s,%s,%s,%2.2f,%2.2f,%2.2f,%2.2f,%s,,%2.2f,%2.2f,%2.2f,%2.2f,%2.2f,\n", SN,op, DATE, TIME, BEMFfres,BEMFdecay,BEMFresrise,Rdc_final,Short_pass,LASERfres,LASERdecay,LASERresrise,Gain,travel);

# DEBUG ONLY
# printf("%s,%s,%2.2f,%2.2f,%2.2f,%2.2f,%s,,%2.2f,%2.2f,%2.2f,%2.2f,%2.2f,\n", SN,op, BEMFfres,BEMFdecay,BEMFresrise,Rdc_final,Short_pass,LASERfres,LASERdecay,LASERresrise,Gain,travel);

close(fn1);
sleep(10);

clear(SN);
# clearall();

// DEBUG2=1; // Testing at home without lasers 

///////////////// NEW ////////////////////////////////////////////////
// }
# }
# }
# }

//////////////////////////////////////////////////////////////////////
