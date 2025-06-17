require savepdf2 coarse_motion_linearity3 extract_positions runavg subsample blink_until_pressed verify_wiring temperature
zzz=ls("*prod*.dat");
# for(FN in zzz){
SKIPCOARSE=1;

DEBUG2=1; // Testing at home without lasers 
if( SKIPCOARSE != 1 ){
if(DEBUG2 != 1){

tic(1);
tic(2);
tic(3);

attr=<<>>;
attr.devname = "ttyUSB";
attr.id_vendor_id = "0403";

if(!exist(r1)){
r1 = usbdevname(attr);
}


opts=<<>>;
opts.data_parity_stop="8N1";
opts.eol="\r";
opts.speed=115200;
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
storage_on = ["SW,OK,02,1\r"]; // ONLY in Communication Mode
storage_off= ["SW,OK,02,0\r"]; // for laser2 ONLY in Communication Mode
sstart = ["AS\r"];      // ONLY in General Mode
sstop  = ["AP\r"];      // ONLY in General Mode
sout   = ["AO,02\r"]; // response AO,hhhhhhhh...,hhhhhhhh ONLY in General Mode

 writem(serial, commmode);
      r = readm(serial);
      // verify if out is laser2
 writem(serial, analogchan12);
      r = readm(serial);
 #     r = readm(serial);
  #    printf("%s\n",r);
      //switch  back to measuremnt mode
 writem(serial, generalmode);
      r = readm(serial);
      
      


attr=<<>>;
attr.devname = "ttyUSB";
attr.id_vendor_id = "10c4";

if(!exist(r2)){
r2 = usbdevname(attr);
}

opts=<<>>;
opts.data_parity_stop="8N1";
opts.eol="\r";
opts.speed=9600;
opts.flow_control="n";

stepper="serial://" + r2;

if(!open(stepper,opts)) {
open(stepper, opts);
}




total=["/1v400Z3000000z200M1000V400P1345M1000D1480M1000P1480M1000D1480M1000v400Z300000M1000V400P126R"];

attr=<<>>;
attr.devname = "ttyACM";
attr.id_vendor_id = "1eab";
 
 
 
if(!exist(r3)){
r3 = usbdevname(attr);
}
 
} 

if( DEBUG2 == 1) { SN="Fake1"; }

if( (DEBUG2 == 1) && (length(strsplt(FN)) > 0 )) { SN=FN; }



 
if(!exist(SN)){
colors("red");
printf("Waiting for S/N...");
 
r=getline(r3,16);
SN=sum (strsplt(r)[1:length(strsplt(r)) -1]);
} 
colors("green");
printf("\rSN:              \t%s\n",SN);
colors();

 gmt= date2jd(gmtime());
 gmt= gmt -floor(gmt);
 gmt=100000*gmt;
 sprintf(dd,"%i.dat",gmt);
printf("COARSE MOTION TEST:\n");
if(DEBUG2 != 1){
blink_until_pressed();
}




if(DEBUG2 != 1){
}

if(DEBUG2 == 1){
system("sudo /home/pi/sandbox/MCP3202_with_pigpio/coarse_30s_MCP3202 > /tmp/datcoarse.dat&");
sleep(30); // to avoid "unable to lock.." message since coarse collectio of data takes 30 s (?)
}

sleep(0.5);
	# DEBUG only

printf("\n\n%s\n","**********************Stepper test starts NOW.");


if(DEBUG2 != 1){

 writem(stepper, total);
      r = readm(stepper);
      #printf("%s\n",r);


 r = readm(stepper);
      #printf("%s\n",r);






sleep(35);  # Thu Oct 18 21:23:52 PDT 2018
 
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
      # sleep(2);
	printf("reading measured laser data\n");

if( DEBUG2 == 1) { system("cp datcoarse.dat /tmp/datcoarse.dat;"); sleep 1; }

xxx=read_ascii("/tmp/datcoarse.dat");
d=d[2:d.nr;];





P4=coarse_motion_linearity3(d, SN);


savepdf2(P4, SN);

E=extract_positions(d);

STEP1= E.STEP1err;
STEP2= E.STEP2err;
STEP3= E.STEP3err;
STEP4= E.STEP4err;


DATE=time2dstr(seconds(),"%Y-%m-%d");
TIME=time2dstr(seconds(),"%H:%M:%S");
op=reads("/var/www/html/Operator.txt");


fn1="/var/www/html/"+DATE+"_coarse_motion.csv";


if(!isfile(fn1)){fprintf(fn1,"%s\n","SN,OPERATOR,DATE,TIME,STEP1,STEP2,STEP3,STEP4");}

if(isfile(fn1)){open(fn1,"a");}

fprintf(fn1,"%s,%s,%s,%s,%2.2f,%2.2f,%2.2f,%2.2f\n", SN,op,DATE,TIME,STEP1,STEP2,STEP3,STEP4);

}
printf("Continue with Twang test\n");
if(DEBUG2 != 1){
blink_until_pressed();
}

system("sudo /home/pi/sandbox/MCP3202_with_pigpio/prodMCP3202 > /tmp/prod.dat");

sleep(4); // make sure data is collected before processing starts



plclose(); // PGPLOT dies if window not closed ?

require rem expm lsim lsim2 tfeval ResRise ap5 damp plot subsample dfdp extract_chunk 
 
 
 
	##  r3 = usbdevname(attr);
 
 

 
	## colors("red");
	## printf("Waiting for S/N...");
	 ## 
	## r=getline(r3,16);
	## # remove trailing \n
	## SN=sum (strsplt(r)[1:length(strsplt(r)) -1]);




tic(4);
if( DEBUG2 == 1) { system("cp OK.prod.dat /tmp/prod.dat;"); sleep 1; }
dir="/tmp/";
file1="prod.dat";
xxx=read_ascii(dir + file1);
pol=verify_wiring(d);
if(pol < 0 ) { 
colors("red");
	printf("%s ","STOP! WRONG POLARITY. Please check coil wires. Press " );
colors("green");
	printf("%s ","GREEN " );
colors("red");
	printf("%s \n","button to restart test from beginning" );

blink_until_pressed();
break; }
DEBUG=0;
if(DEBUG == 1){printf("Data: %s\n",dir+file1);}

d=d[2:d.nr;];

dnew= d[;1] -sum(d[;1])/length(d[;1]);
th=100;
th=200;
th=300 ;
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

RCAL=0.63;
RDCMAX=14.0;
RDCMIN=12.0;
Rdc_final=sum(Rx)/length(Rx) - RCAL;



if( (Rdc_final > RDCMAX) || ( Rdc_final < RDCMIN)){
	 colors("red");
	 printf ("Coil resistance = %2.2f [Ohm]\t\tFAILED. Check connections.\n", Rdc_final );
	 colors();
 }



if( (Rdc_final < RDCMAX) && ( Rdc_final > RDCMIN)){
	 colors("green");
	 printf ("Coil resistance = %2.2f [Ohm]\t\tPASSED.\n", Rdc_final);
	 colors();
 }




if (abs(sum(Ip-Im)) > 0.1* abs(sum(Im))) {
 
	 colors("red");
	 printf ("Leakage test \t\t\t\tFAILED. Coil shorted to ground.\n"); colors();
Short_pass="FAIL";
 }


if (abs(sum(Ip-Im)) < 0.1* abs(sum(Im))) {
	 colors("green");
	 printf ("Coil leakage to ground\t\t\tPASSED.\n");
	 colors();
Short_pass="PASSED";
 }



ind3=find(ddiff2);
BEMF_extract_start = ind3[(find((diff(ind3) >24) && (diff(ind3) <32 )))];
k = max(size(BEMF_extract_start))-1 ;

Lsr_extract_start = BEMF_extract_start ;


skip=30;
skip=120; # *** last used
skip=80;

	                 

if(DEBUG == 1){printf("BEMF k= %i\n",k); }

plot ([4e-5*[BEMF_extract_start[k]+skip : BEMF_extract_start[k]+1000 ]', [d[BEMF_extract_start[k]+skip : BEMF_extract_start[k]+1000;2]]]);

dc = sum(d[BEMF_extract_start[k]+500 : BEMF_extract_start[k]+1000;2]) / length(d[BEMF_extract_start[k]+500 : BEMF_extract_start[k]+1000;2]);
x = 4e-5*[BEMF_extract_start[k]+skip : BEMF_extract_start[k]+1000 ]';
x=x - x[1]; # reset time to zero
y = d[BEMF_extract_start[k]+skip : BEMF_extract_start[k]+1000;2] -dc;


func = function(x,p){
	global(u2,pi)

    return p[4] + p[3]*exp(-(x)/p[1]).*sin(2*pi*p[2].*(x-p[5]));
    // return p[4] + p[3]*exp(-(x)/p[1]).*sin(2*pi*p[2].*(x-p[5]));
    };






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

 zls=lsfit(y,x,p,func,dfdp);


     #printf("\n****************************************************\n");
	#printf("T decay = %2.2f [ms] \n",1e3*zls.coef[1]);

	#colors("red");
	#printf("T decay = %2.2f [ms] \n",1e3*zls.coef[1],,);
	#colors();




pltitle("Fitting measured data ");
xlabel("Time [s]");
ylabel("Amplitude");
plegend(["Data","Lsfit"]);
if (_rlab_config.plot_support=="pgplot")
	    {
               plaxis("lin",); # pgplot syntax
     else
              plscale("lin",); # plplot syntax
   }

      # plaxis("lin",); # pgplot syntax
      # plscale("lin",); # plplot syntax
      # plaxis("lin","lin");
zzls=(func(x,zls.coef));


plot(<<[x,y];[x,zzls]>> );
sleep(3);



lsnum2=[1,2*pi*zls.coef[2]];
lsden2=[1,2/zls.coef[1],1/(zls.coef[1]*zls.coef[1])+2*pi*2*pi*(zls.coef[2]*zls.coef[2])]; 

lsnum2=[zls.coef[1],0];
lsden2=[1,2/zls.coef[1],1/(zls.coef[1]*zls.coef[1])+2*pi*2*pi*(zls.coef[2]*zls.coef[2])]; 


dnum2=[zls.coef[1]];
dden2= lsden2;

freq = [10:3000:10];
inum2 = [1, 0];
itf = tfeval(inum2, lsden2, freq);
idb=20*log10(abs(itf));



xlabel("Frequency [Hz]");
if (_rlab_config.plot_support=="pgplot")
	    {
               plaxis("log",); # pgplot syntax
     else
              plscale("log",); # plplot syntax
   }

ylabel("Amplitude [dB]");
plegend();
pltitle("Frequency response");
plegend(["Back EMF","Displacement"]);
plegend(["Back EMF", "Input data","Displacement"]);
plegend(["Back EMF", "Input data"]);



dtf = tfeval(dnum2, dden2, freq);
dtf0 = tfeval(dnum2, dden2, freq[1]);

dispdb = 20*log10(abs(dtf)) - 20*log10(abs(dtf0));

xlabel("Frequency [Hz]");
ylabel("Amplitude [dB]");
pltitle("Displacement Frequency Response");
if (_rlab_config.plot_support=="pgplot")
	    {
               plaxis("log",); # pgplot syntax
     else
              plscale("log",); # plplot syntax
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



fr=X_extracted.wn[1]/(2*pi) ;
printf("Frequency Response Peak = %2.2f [dB]\n",X_extracted.FRpeaking);

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



frrise = ResRise (1e3*zls.coef[1], fr);
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

printf("Elapsed time %3.2f s\n",toc(4));

BEMFfres=fr;
BEMFdecay=1e3*zls.coef[1];
BEMFresrise=frrise;


plclose();

um = 5.17;

skip=120;
	                 
if(DEBUG == 1){printf("Laser k= %i\n",k); }

pltitle("Fitting measured laser data ");
xlabel("Time [s]");
ylabel("Amplitude [um]");
if (_rlab_config.plot_support=="pgplot")
	    {
               plaxis("lin",); # pgplot syntax
     else
              plscale("lin",); # plplot syntax
   }
plot ([4e-5*[Lsr_extract_start[k]+skip : Lsr_extract_start[k]+1000 ]', [um*d[Lsr_extract_start[k]+skip : Lsr_extract_start[k]+1000;5]]]);

dc = sum(d[Lsr_extract_start[k]+500 : Lsr_extract_start[k]+1000;5]) / length(d[Lsr_extract_start[k]+500 : Lsr_extract_start[k]+1000;5]);
x = 4e-5*[Lsr_extract_start[k]+skip : Lsr_extract_start[k]+1000 ]';
x=x - x[1]; # reset time to zero
y = d[Lsr_extract_start[k]+skip : Lsr_extract_start[k]+1000;5] -dc;





p=[ 0.00327633722,    84.8807812,    133.730521,   0.522658955,  -0.00134066124 ]'; # for laser data 
p=[ 0.00296779553,    86.1768055,    132.631621,   0.502642791,  -0.00124655545 ]'; # for laser data 


xy = [x,y];

options=<<>>;
options.imethod = 2;  // least squares fit

 zls=lsfit(y,x,p,func,dfdp);



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

      #plaxis("lin",); # pgplot syntax
      # plaxis("lin","lin");
zzls=(func(x,zls.coef));


plot(<<[x,um*y];[x,um*zzls]>> );
sleep(3);



lsnum2=[1,2*pi*zls.coef[2]];
lsden2=[1,2/zls.coef[1],1/(zls.coef[1]*zls.coef[1])+2*pi*2*pi*(zls.coef[2]*zls.coef[2])]; 

lsnum2=[zls.coef[1],0];
lsden2=[1,2/zls.coef[1],1/(zls.coef[1]*zls.coef[1])+2*pi*2*pi*(zls.coef[2]*zls.coef[2])]; 


dnum2=[zls.coef[1]];
dden2= lsden2;

freq = [10:3000:10];
inum2 = [1, 0];
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

plegend();
pltitle("Frequency response");


dtf = tfeval(dnum2, dden2, freq);
dtf0 = tfeval(dnum2, dden2, freq[1]);

dispdb = 20*log10(abs(dtf)) - 20*log10(abs(dtf0));

xlabel("Frequency [Hz]");
ylabel("Amplitude [dB]");
pltitle("Displacement Frequency Response");

if (_rlab_config.plot_support=="pgplot")
	    {
               plaxis("log",); # pgplot syntax
     else
              plscale("log",); # plplot syntax
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


frrise = ResRise (1e3*zls.coef[1], fr);
colors("yellow");
printf("From Laser Measurements:\n");
if(frrise <= RESRISEmax && frrise > RESRISEmin ){
colors("green");
printf("Resonant rise = %2.2f [dB] \t\tPASSED.\n",frrise);
colors();
}

if(frrise > RESRISEmax || frrise < RESRISEmin ){
colors("red");
printf("Resonant rise = %2.2f [dB] \t\tFAILED.\n",frrise,,);
colors();
}

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

	// add writing to a file:w
DATE=time2dstr(seconds(),"%Y-%m-%d");
TIME=time2dstr(seconds(),"%H:%M:%S");
op=reads("/var/www/html/Operator.txt");


fn1="/var/www/html/"+DATE+"_twang.csv";

TADDR1="tbd";
Tamb=temperature(TADDR1);





if(isfile(fn1)){open(fn1,"a");}

if(!isfile(fn1)){fprintf(fn1,"%s\n","SN,OPERATOR,DATE,TIME,BEMFfres,BEMFdecay,BEMFresrise,Rdc_final,Short_pass,,LASERfres,LASERdecay,LASERresrise,Gain,travel,Tambient");}


LASERfres=fr;
LASERdecay=1e3*zls.coef[1];
LASERresrise=frrise;
Gain=gain;


fprintf(fn1,"%s,%s,%s,%s,%2.2f,%2.2f,%2.2f,%2.2f,%s,,%2.2f,%2.2f,%2.2f,%2.2f,%2.2f,%s,\n", SN,op, DATE, TIME, BEMFfres,BEMFdecay,BEMFresrise,Rdc_final,Short_pass,LASERfres,LASERdecay,LASERresrise,Gain,travel,Tamb);


close(fn1);aaa
sleep(10);

clear(SN);



}

