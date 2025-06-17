#!/home/pi/bin/rlab2
#!/home/azg/bin/rlab2
require rem expm lsim lsim2 tfeval ResRise ap5 damp plot subsample dfdp extract_chunk psprint 
# azg_barcode_reader_usb





# "line 2"
######################################################################
# rfile azg_barcode_reader_usb

# DEBUG and TESTING only
if(!exist("/tmp/no_timestamp_production_vol=0.2.dat")) {error ("Missing ADC data");}
######################################################################
tic();xxx=read_ascii("/tmp/no_timestamp_production_vol=0.2.dat");
# tic();xxx=read_ascii("no_timestamp_production_vol=0.2.dat");
# printf("File read in %f [sec]\n",toc());
# flush first row
d=d[2:d.nr;];

dnew= d[;1] -sum(d[;1])/length(d[;1]);
th=100;
th=200;
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
	# Im=(Im - 2048) * (101+13+101)/(200*0.5*13)/4096*3.3;
	# Exact formula:
	Im = (Im -2048) /4096 *3.3 / (1/(1/0.5+1/(100+13+100)) * 13/(100+13+100) * 200 );
	Ip=B2[;2];
	# Ip=(Ip - 2048) * (101+13+101)/(200*0.5*13)/4096*3.3;
	# Exact formula:
	Ip = (Ip -2048) /4096 *3.3 / (1/(1/0.5+1/(100+13+100)) * 13/(100+13+100) * 200 );


	N=200;
	# L
	L[1]=1;
	L[2]=length(Vm);

	# "Rx="
	Rx[k] = sum((Vp[L[2]-N:L[2]] -Vm[L[2]-N:L[2]])./(Im[L[2]-N:L[2]] ))/(N+1) -1.0;
	k=k+1;

}

for(j in 2:length(Rdc_extract_start)){ Rdc[j-1] = (Rx[j]+Rx[j-1])/2; } 

Rdc_final=sum(Rx)/length(Rx);

# DEBUG
# Rdc_final =  Rdc_final +1;


if( (Rdc_final > 14.0) || ( Rdc_final < 12.0)){
	 colors("red");
	 printf ("Coil resistance \t\t\tFAILED. Check connections.\n");
	 colors();
 }



if( (Rdc_final < 14.0) && ( Rdc_final > 12.0)){
	 colors("green");
	 printf ("Coil resistance = %2.2f [Ohm]\t\tPASSED.\n", Rdc_final);
	 colors();
 }



# Leakage test

if ((sum(Ip-Im)) > 0.1* sum(Im)) {
# DEBUG failure
# if ((sum(Ip-Im)) < 0.1* sum(Im)) {
 
	 colors("red");
	 printf ("Leakage test \t\t\tFAILED. Coil shorted to ground.\n");
	 colors();
 }


if ((sum(Ip-Im)) < 0.1* sum(Im)) {
	 colors("green");
	 printf ("Coil leakage to ground\t\t\tPASSED.\n");
	 colors();
 }


# 	 printf ("Color check\n");
######################################################################

# " BEMF begin at:"
BEMF_extract_start = ind[find((diff(ind) >2500) && (diff(ind) <3000 ))];

# "Sanity check:"
# for(i in 2:length(BEMF_extract_start)) { BEMF_extract_start[i] - BEMF_extract_start[i-1]  }

skip=30;
# plot ([4e-5*[BEMF_extract_start[1]+skip : BEMF_extract_start[1]+1000 ]', [d[BEMF_extract_start[1]+skip : BEMF_extract_start[1]+1000;2]]]);

	                 
k=6;
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


xy = [x,y];

options=<<>>;
options.imethod = 2;  // least squares fit
# options.stdout = rconsole();

 zls=lsfit(y,x,p,func,dfdp);
# zls.coef


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
      plaxis("lin",); # pgplot syntax
      #plscale("lin",); # plplot syntax
      # plaxis("lin","lin");
# " Debug 116"
zzls=(func(x,zls.coef));


plot(<<[x,y];[x,zzls]>> );
sleep(2);

# _plprint ("time_data.eps", "psc");

# " Debug 126"

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
# plscale("log",); # plplot syntax
plaxis("log",);
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

dtf = tfeval(dnum2, dden2, freq);
dtf0 = tfeval(dnum2, dden2, freq[1]);

# WHy this worked before? [azg] Sat Apr 21 17:32:09 PDT 2018
# dispdb=20*log10(inum2/dnum2) + 20*log10(abs(dtf));
# dispdb= 20*log10(abs(dtf));
dispdb = 20*log10(abs(dtf)) - 20*log10(abs(dtf0));

xlabel("Frequency [Hz]");
ylabel("Amplitude [dB]");
pltitle("Displacement Frequency Response");
plaxis("log",);
# plscale("log",); # plplot syntax
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


# frrise = ResRise (Tau_decay, fres);
# frrise = ResRise (1e3*z.coef[1], z.coef[2]);
frrise = ResRise (1e3*zls.coef[1], zls.coef[2]);
##printf("\n****************************************************\n");
#printf("Calculated Fres rise using IBM formula  = %f [dB]\n",frrise);
#printf("\n****************************************************\n");
# res rise 0.8 =/- 2 dB
# +2.8 -1.2 dB
if(frrise <= 2.8 && frrise > -1.2 ){
colors("green");
printf("Resonant rise = %2.2f [dB] \t\tPASSED.\n",frrise);
colors();
}

if(frrise > 2.8 || frrise < -1.2 ){
colors("red");
printf("Resonant rise = %2.2f [dB] \t\tFAILED.\n",frrise,,);
colors();
}
#printf("\n****************************************************\n");


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
#printf("Frequency response peak = %2.2f [Hz]\n",X_extracted.Fringing);
#printf("Actual Resonant Rise    = %2.2f [dB]\n",X_extracted.FRpeaking);
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

     #printf("\n****************************************************\n");
if(1e3*zls.coef[1] <= 4.4 && 1e3*zls.coef[1] > 2.76 ){
colors("green");
printf("T decay = %2.2f [ms] \t\t\tPASSED.\n",1e3*zls.coef[1]);
colors();
}

if(1e3*zls.coef[1] > 4.4 || 1e3*zls.coef[1] < 2.76 ){
colors("red");
printf("T decay = %2.2f [ms] \t\t\tFAILED.\n",1e3*zls.coef[1],,);
colors();
}




     #printf("\n****************************************************\n");
fr=X_extracted.wn[1]/(2*pi) ;
if( (fr <120) && (fr > 80 )){
colors("green");
printf("Resonant frequency = %2.2f [Hz]\t\tPASSED.\n",fr);
colors();
}


if( (fr >120) || (fr < 80 )){
colors("red");
printf("Resonant frequency = %2.2f [Hz]\t\t\tFAILED.\n",fr);
colors();
}
printf("Elapsed time %3.2f s\n",toc());
sleep(10);
