require rem expm lsim lsim2 tfeval ResRise ap5 damp plot

# ToDo
# use min(zzls)/max(zzls)
# or better use measured data input (instead:  min(J2.y)/max(J2.y)   )


# to determine if ringing is possible


# Specs:
# Rdc    Fres    Decay [ms]      Fres rise
# 13±1	100±20	3.7+2.2/-1.4	0.8±2 dB

rng(1,"normal", [0,1]);
plclose();



// RLC
L=15.9e-3;
C=.159e-3;
R=10;
R=13;
R=4;
# R=18;
num2=[1/L,0];
den2=[1,R/L,1/(L*C)];
t2=1e-3*(0.0:200.0:0.2)';
t2=1e-3*(0.0:40.0:0.2)';
period=100e-3*ones(t2.nr,t2.nc);
period=50e-3*ones(t2.nr,t2.nc);
period=200e-3*ones(t2.nr,t2.nc);
u2=(rem(t2,period) >= period ./ 2);
u2=-(u2-1);
# noise added to input
# u2=-(u2-1) + 0.1*rand(u2.nr,1); 

_ctb2_window = plstart(2,3); ## ??
subplot(1);
plegend();
# plaxis("lin","lin"); # pgplot syntax
J=lsim([1/L],den2,u2,t2);
# sleep(3);
J=lsim(num2,den2,u2,t2);
# sleep(3);
J2=lsim2(num2,den2,u2,t2);
# use as if measured data
x=t2;
y=100*J2.y;
# noise added to output
y=y + 0.1*rand(y.nr,1); 

###############################################
func = function(x,p){
	global(u2,pi)

    return p[3]*exp(-x/p[1]).*sin(2*pi*p[2].*x);
    };
    # Note: should have done it in omega instead of f
######################################################################


dfdp = function(x,p){
	global(pi)
	// jacobian df/dp
      y=[];
     y[;1] = (p[3].*x.*exp(-x/p[1]).*sin(2*pi*p[2]*x))/(p[1].*p[1]);
     y[;2] = (2*pi*p[3].*x.*exp(-x/p[1]) .* cos(2*pi*p[2].*x));
     y[;3] = (exp(-x/p[1]).* sin(2*pi*p[2].*x));

return y;
};

dfdx = function(x,p){
	global(pi)
	// jacobian df/dx
	y=[];
	y[;1] = (2*pi*p[2]*p[3]*exp(-x/p[1]) .* cos(2*pi*p[2].*x) - (p[3]*exp(-x/p[1])).* sin(2*pi*p[2].*x)/p[1]);

return y;
};


// nonlinear fit
// initial guess
p=[];
p=[3e-3,2*pi*100]';
p=[3e-3,100,10]';
p=[6e-3,200,20]';
#   T   f   Amp
p=[1e-3,200,20]';
p=[1e-4,200,20]';


xy = [x,y];

options=<<>>;
options.imethod = 2;  // least squares fit
# options.stdout = rconsole();

# z=odrfit(y,x,p,func,dfdp,dfdx,options);
#   z=odrfit(y,x,p,func,dfdp,options);
# z.coef'
# p

# z.coef-p
 zls=lsfit(y,x,p,func,dfdp);
"zls.coef lsfit"
zls.coef


printf("lsfit T decay = %f [ms] \n",1e3*zls.coef[1]);
printf("lsfit F res   = %f  [Hz]\n",zls.coef[2]);


# printf("Amplitude = %f\n",z.coef[3]);
# z.cov

# plot([x,y]);
# sleep(1);
# plot([x,func(x,z.coef)]);

pltitle("Fitting measured data + noise");
xlabel("Time [s]");
ylabel("Amplitude");
plegend(["Data","Lsfit"]);
      plaxis("lin",); # pgplot syntax
      #plscale("lin",); # plplot syntax
      # plaxis("lin","lin");
" Debug 116"
# zz=(func(x,z.coef));
zzls=(func(x,zls.coef));


plot(<<[x,y];[x,zzls]>> );
# sleep(2);

# _plprint ("time_data.eps", "psc");

" Debug 126"
# plwins(2);
# plwin(2);
" Debug 129"
# num2=[1,z.coef[1]];
# num2=[2*pi*z.coef[2]];
# den2=[1,z.coef[1],z.coef[2]]; 
# den2=[1,2/z.coef[1],1/(z.coef[1]*z.coef[1])+2*pi*2*pi*(z.coef[2]*z.coef[2])]; 

# from lsfit
lsnum2=[1,2*pi*zls.coef[2]];
lsden2=[1,2/zls.coef[1],1/(zls.coef[1]*zls.coef[1])+2*pi*2*pi*(zls.coef[2]*zls.coef[2])]; 
# missing 2*pi for s^1 den2
# lsden2=[1,2*2*pi'a/zls.coef[1],1/(zls.coef[1]*zls.coef[1])+2*pi*2*pi*(zls.coef[2]*zls.coef[2])]; 
# ?? Decay time constant is time (no 2*pi needed)
lsden2=[1,2/zls.coef[1],1/(zls.coef[1]*zls.coef[1])+2*pi*2*pi*(zls.coef[2]*zls.coef[2])]; 

# Need the same form as lsim2 input ( add 1/R factor to coef[1] 
# This is transfer function of back EMF
lsnum2=[1/R*zls.coef[1],0];
lsden2=[1,2/zls.coef[1],1/(zls.coef[1]*zls.coef[1])+2*pi*2*pi*(zls.coef[2]*zls.coef[2])]; 

# the displacement is proportional to current
# since V = dI/dt * L
# current is Back EMF(s) * 1/sL

dnum2=[1/(R*L)*zls.coef[1]];
dden2= lsden2;

freq = [10:3000:10];
# input tf
inum2 = num2[1];
      itf = tfeval(inum2, den2, freq);
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

      tf = tfeval(num2, den2, freq);
      db=20*log10(abs(tf));

      lstf = tfeval(lsnum2, lsden2, freq);
      # lsdb=20*log10(abs(lstf));
      lsdb=20*log10(num2[1]/lsnum2[1]) + 20*log10(abs(lstf));

plot(<<[freq;db]';[freq;lsdb]'>>);
# sleep(3);

      dtf = tfeval(dnum2, dden2, freq);
      dispdb=20*log10(inum2/dnum2) + 20*log10(abs(dtf));

   xlabel("Frequency [Hz]");
   ylabel("Amplitude [dB]");
   pltitle("Displacement Frequency Response");
      plaxis("log",);
      plegend(["Input data","Displacement"]);

plot(<< [freq;idb]' ;[freq;dispdb]'>>);
# plot(<<[freq;lsdb]';[freq;dispdb]'>>);
# plot(<<[freq;dispdb]'>>);
     sleep(2);


# _plprint ("frequencyResponse.eps", "psc");


# frrise = ResRise (Tau_decay, fres);
# frrise = ResRise (1e3*z.coef[1], z.coef[2]);
frrise = ResRise (1e3*zls.coef[1], zls.coef[2]);
printf("Calculated Fres rise from Lsfit  = %f [dB]\n",frrise);

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
# from ResRise.r Tau_decay.r 
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

printf("zeta_extracted = %f\n",zeta_extracted);

wn_extracted = zls.coef[2].*(2*pi) ./ sqrt(1- zeta_extracted[1]*zeta_extracted[1]) ;
printf("wn_extracted = %f\n",wn_extracted);

# ap5([freq;db]');


den_extracted = [1, 2*zeta_extracted*wn_extracted,  wn_extracted^2 ]; 

 X_extracted = damp (den_extracted);
X=damp(den2);

# NOTE: we  print only first number from every member (eg. Poles have many numbers)
# printf("%s\n", "Input data");
# for(k in members(X)){printf("%s\r\t\t = %f\n",k,X.[k]);}
# printf("%s\n", "Extracted data");
# for(k in members(X_extracted)){printf("%s\r\t\t = %f\n",k,X_extracted.[k]);}


printf("\n\n\t\t%s\t\t%s\t\t%s\n","   Input data","   Extracted","Error %");
for(k in members(X)){ printf("%s\r\t\t = %f",k,X.[k]); printf("\t\t = %f\t\t%f\n",X_extracted.[k],(X.[k]-X_extracted.[k])/X.[k]*100);}

