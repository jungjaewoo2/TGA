require rem expm lsim lsim2 tfeval ResRise ap5 damp plot subsample dfdp extract_chunk psprint



# Specs:
# Rdc    Fres    Decay [ms]      Fres rise
# 13±1	100±20	3.7+2.2/-1.4	0.8±2 dB
#               3.5+0.9/-0.74

###############################################
func = function(x,p){
	global(u2,pi)

    return p[4] + p[3]*exp(-(x)/p[1]).*sin(2*pi*p[2].*(x-p[5]));
    };
######################################################################


##p = function(x,p){
	##global(pi)
	##// jacobian df/dp
      ##y=[];
     ##y[;1] = (p[3].*x.*exp(-x/p[1]).*sin(2*pi*p[2]*x))/(p[1].*p[1]);
     ##y[;2] = (2*pi*p[3].*x.*exp(-x/p[1]) .* cos(2*pi*p[2].*x));
     ##y[;3] = (exp(-x/p[1]).* sin(2*pi*p[2].*x));
     ##y[;4] = 1;
     ##y[;5] =-2*p[2]* p[3]*pi*exp((-x)./p[1]) .* cos(2*p[2]*pi*(x-p[5]));

##return y;
##};

##dfdx = function(x,p){
	##global(pi)
	##// jacobian df/dx
	##y=[];
	##y[;1] = (2*pi*p[2]*p[3]*exp(-x/p[1]) .* cos(2*pi*p[2].*x) - (p[3]*exp(-x/p[1])).* sin(2*pi*p[2].*x)/p[1]);
##
##return y;
##};

#################################################################################


# for( i in 10:20:1){
#  plclose();
        #_ctb2_window = plstart(1,2);



# _ctb2_window = plstart(1,2); ## ??

#// Plot the results
if (_rlab_config.plot_support=="pgplot"||_rlab_config.plot_support=="plplot")
{
if (exist(_ctb2_window))
{
        plwin(_ctb2_window);
     else
        _ctb2_window = plstart(1,2);
     }
   }







subplot(1);
plegend();


############################################################################
# impulse only
############################################################################
# NEW Sat Apr 21 17:43:05 PDT 2018
# xxx=readb("impulseBEMF3.mat");
xxx=read_ascii("no_time_stamp_noise_debug8.dat");
 # idx= extract_chunk(d[;2],  3); diff(idx)*4e-5
# "d" name is hardcoded in ascii file

chunk=3;
m=1;
############################################################################
for(j in 1:60){
# for(j in 40){

############################################################################
dc1=[];
chunk = j;
# chunk = 5;
idx= extract_chunk(d[;1],  chunk, 30 ); 
# idx= extract_chunk(d[;1],  chunk ); 
if (  (diff (idx)*4e-5 > 0.09)  && (diff (idx)*4e-5 < 0.11)  ) {
# if (  (diff (idx)*4e-5 < 0.09)  || (diff (idx)*4e-5 > 0.11)  ) {error ( "wrong size of chunk of data\n"); }

# BEMF = d[;2][idx[1]:idx[2]];
# limit the length to
LL = 800;
BEMF = d[;2][ idx[1] : idx[1] + LL ];
#dc=sum( [d;2][idx[1]+1000:idx[2]-100] )  / length ( [d;2][idx[1]+1000:idx[2]-100] );
dc = sum(BEMF[400:800])/length(BEMF[400:800]);
if(! exist(dc1)){dc1=dc }


############################################################################
#x=[30:length(BEMF)]' *4e-5;
skip=1;
skip=20;
skip=25;
skip=30;
skip=40;
skip=50;
skip=100;
# skip=200; # dies
skip=1;
# x=[10:length(BEMF)]' *4e-5;
x=[skip:length(BEMF)]' *4e-5;
x=x - x[1]; # reset time to zero

#y=BEMF[30:length(BEMF)];
#y=BEMF[10:length(BEMF)];
y=BEMF[skip:length(BEMF)] -dc;
############################################################################
# [azg] Wed Apr 25 00:39:39 PDT 2018
ss=5;
y=subsample(y,ss);
x=subsample(x,ss);

# to get impulse response from step response we neet to differentiate
# squarewave starts after ~ chunk=12 ?
# this is correct but we have to deal with very noisy signal.

if(chunk > 12){y=diff(y); x=x[1:length(x)-1]; }

############################################################################
Vm=d[;1];
Vm=Vm * (511+1333)/511;
Vm=Vm ./4096*3.3;
Vp=d[;3];
Vp=Vp * (511+1333)/511;
Vp=Vp ./4096*3.3;
Im=d[;4];
Ip=d[;2];

# current in A
Im = (Im -2048) /4096 *3.3 / (1/(1/0.5+1/(100+13+100)) * 13/(100+13+100) * 200 );
Ip = (Ip -2048) /4096 *3.3 / (1/(1/0.5+1/(100+13+100)) * 13/(100+13+100) * 200 );

Isteady = (dc -2048) /4096 *3.3 / (1/(1/0.5+1/(100+13+100)) * 13/(100+13+100) * 200 );






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


xy = [x,y];

options=<<>>;
options.imethod = 2;  // least squares fit
# options.stdout = rconsole();

 zls=lsfit(y,x,p,func,dfdp);
# zls.coef


     printf("\n****************************************************\n");
if(1e3*zls.coef[1] <= 4.4 && 1e3*zls.coef[1] > 2.76 ){
printf("T decay = %f [ms] \n",1e3*zls.coef[1]);
}

if(1e3*zls.coef[1] > 4.4 || 1e3*zls.coef[1] < 2.76 ){
colors("red");
printf("T decay = %f [ms] \n",1e3*zls.coef[1],,);
colors();
}



printf("F res   = %f  [Hz]\n",zls.coef[2]);

pltitle("Fitting measured data ");
xlabel("Time [s]");
ylabel("Amplitude");
plegend(["Data","Lsfit"]);
      #plaxis("lin",); # pgplot syntax
      plscale("lin",); # plplot syntax
      # plaxis("lin","lin");
# " Debug 116"
zzls=(func(x,zls.coef));


plot(<<[x,y];[x,zzls]>> );
# sleep(2);

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
      plscale("log",); # plplot syntax
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
      plscale("log",); # plplot syntax
      plegend(["Input data","Displacement"]);
      plegend(["Displacement"]);

# " Debug 199"
plot(<< [freq;dispdb][;1:30]'>>);
#sleep(2);
# plot(<<[freq;lsdb]';[freq;dispdb]'>>);
# plot(<<[freq;dispdb]'>>);
     # sleep(2);


# _plprint ("frequencyResponse.eps", "psc");


# frrise = ResRise (Tau_decay, fres);
# frrise = ResRise (1e3*z.coef[1], z.coef[2]);
frrise = ResRise (1e3*zls.coef[1], zls.coef[2]);
     #printf("\n****************************************************\n");
printf("Calculated Fres rise using IBM formula  = %f [dB]\n",frrise);
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
if( (fr <120) && (fr > 80 )){
printf("Resonant frequency      = %f [Hz]\n",fr);
}
if( (fr >120) || (fr < 80 )){
colors("red");
printf("Resonant frequency      = %f [Hz]\n",fr);
colors();
}
# printf("Resonant frequency      = %f [Hz]\n",X_extracted.wn/(2*pi));
printf("Frequency response peak = %f [Hz]\n",X_extracted.Fringing);
printf("Actual Resonant Rise    = %f [dB]\n",X_extracted.FRpeaking);
printf("Coil current             = %f [A]\n",Isteady);
printf("chunk number             = %i \n",chunk);
"zls.coef"
zls.coef
     printf("\n****************************************************\n");


# decay time constant
DT[m;]=[Isteady*1e3, 1e3*zls.coef[1]]; 
m=m+1;
# printf("\n\n\t\t%s\t\t%s\t\t%s\n","   Input data","   Extracted","Error %");
# for(k in members(X)){ printf("%s\r\t\t = %f",k,X.[k]); printf("\t\t = %f\t\t%f\n",X_extracted.[k],(X.[k]-X_extracted.[k])/X.[k]*100);}
# toc()
# }
}}



plstart(1,1);
plstyle("point");
xlabel("Coil Current [mA]");
ylabel("Decay Time [ms]");
plegend();
pltitle("Decay Time as a Function of Displacement");
plpoint([2,2,2,2,2,2,2,2]);
plimits(-30,+30,0,6);
plot(DT);
 _plprint ("Decay_vs_Displacement.ps2", "psc");
psprint ("Decay_vs_Displacement.pdf");
