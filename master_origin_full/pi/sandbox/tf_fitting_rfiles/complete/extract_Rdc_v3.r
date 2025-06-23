


tic();xxx=read_ascii("no_timestamp_production_vol=0.2.dat");
printf("File read in %f [sec]\n",toc());

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
" find((diff(ind) >5000) && (diff(ind) <6000 ))"
ind2=find((diff(ind) >5000) && (diff(ind) <6000 ))

for(i in 2:length(ind2)) { ind[ind2[i]] - ind[ind2[i-1]]  }

" steps begin at:"
Rdc_extract_start=ind[find((diff(ind) >5000) && (diff(ind) <6000 ))]


 # plot([4e-5*[1:length(d)]', [d;1]]);
 plot([4e-5*[ind[ind2[1]]:ind[ind2[4]]+7000]', [d[ind[ind2[1]]:ind[ind2[4]]+7000 ;1:2]]]);

ind[ind2]*4e-5
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
	 printf ("Coil resistance \tFAILED. Check connections.\n");
	 colors();
 }



if( (Rdc_final < 14.0) && ( Rdc_final > 12.0)){
	 colors("green");
	 printf ("Coil resistance = %2.2f\t\tPASSED.\n", Rdc_final);
	 colors();
 }



# Leakage test

if ((sum(Ip-Im)) > 0.1* sum(Im)) {
# DEBUG failure
# if ((sum(Ip-Im)) < 0.1* sum(Im)) {
 
	 colors("red");
	 printf ("Leakage test \t\tFAILED. Coil shorted to ground.\n");
	 colors();
 }


if ((sum(Ip-Im)) < 0.1* sum(Im)) {
	 colors("green");
	 printf ("Coil leakage to ground\t\tPASSED.\n");
	 colors();
 }


# 	 printf ("Color check\n");
