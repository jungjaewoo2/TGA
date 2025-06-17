for(k in 19:30){
d=read("R=12.04_chunk"+num2str(k)+".mat");

d=d1; # as stored in R=12.04_chunk19.mat
t=d[;1].*4e-5;
Vm=d[;2];
Vm=Vm * (511+1333)/511;
Vm=Vm ./4096*3.3;
Vp=d[;4];
Vp=Vp * (511+1333)/511;
Vp=Vp ./4096*3.3;
Im=d[;5];
# Im=(Im - 2048) * (101+13+101)/(200*0.5*13)/4096*3.3;
# Exact formula:
Im = (Im -2048) /4096 *3.3 / (1/(1/0.5+1/(100+13+100)) * 13/(100+13+100) * 200 );
Ip=d[;3];
# Ip=(Ip - 2048) * (101+13+101)/(200*0.5*13)/4096*3.3;
# Exact formula:
Ip = (Ip -2048) /4096 *3.3 / (1/(1/0.5+1/(100+13+100)) * 13/(100+13+100) * 200 );

# Exact formula:
# I = (ADCcount -2048) /4096 *3.3 / (1/(1/0.5+1/(100+13+100)) * 13/(100+13+100) * 200 )



#M=int(0.11099/4e-5)
N=200;
#  Ileakm=Vm/(511+1330) ;
#  Ileakp=Vp/(511+1330) ;
#  L=extract_chunk(Vm,21);
# "L"
# L
L[1]=1;
L[2]=length(t);

# "Rx="
# Rx = sum((Vp[L[2]-N:L[2]] -Vm[L[2]-N:L[2]])./(Im[L[2]-N:L[2]] + Ileakm[L[2]-N:L[2]]))/(N+1)
# Rx = sum((Vp[M:M+N] -Vm[M:M+N])./(Im[M:M+N] + Ileakm[M:M+N]))/(N+1)
Rx[k] = sum((Vp[L[2]-N:L[2]] -Vm[L[2]-N:L[2]])./(Im[L[2]-N:L[2]] ))/(N+1) -1.0;

# BackEMF=(Vp-Vm)-(Im+Ileakm)*Rx;


# pltitle("BackEMF");
# # xaxis("Time [s]");
# xlabel("Time [s]");
# ylabel("BackeEMF [V]");
# plegend();
# ylabel("BackEMF [V]");
# plot([t,BackEMF][L[1]:L[2];])
}
for(j in 20:30){(Rx[j]+Rx[j-1])/2} 
