require extract_chunk


#for(k in 19:30){
k=1;
#d=read("R=12.04_chunk"+num2str(k)+".mat");

# xxx=read_ascii("no_timestamp_production.dat");
xxx=read_ascii("no_timestamp_production_vol=0.2.dat");


####################################
for(j in 1:100){
# for(j in 1:40){

dc1=[];
chunk = j;
# chunk = 5;
# idx= extract_chunk(d[;1],  chunk, 30 ); 
# idx= extract_chunk(d[;1],  chunk, 400 ); 
idx= extract_chunk(d[;1],  chunk, 100 );
# idx= extract_chunk(d[;1],  chunk ); 
# if (  (diff (idx)*4e-5 > 0.09)  && (diff (idx)*4e-5 < 0.11)  ) {
if (  (diff (idx) > 4000 )  && (diff (idx) < 6000 )  ) {
	# if (  (diff (idx)*4e-5 < 0.09)  || (diff (idx)*4e-5 > 0.11)  ) {error ( "wrong size of chunk of data\n"); }

# Modified for data without timestamps
####################################
	# BEMF = d[;2][idx[1]:idx[2]];
	# limit the length to
	LL = 800;
	B2 = d[ idx[2] - LL: idx[2] ;];


# d=d1; # as stored in R=12.04_chunk19.mat
### t=d[;1].*4e-5; # no time stamps in input file
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




#M=int(0.11099/4e-5)
N=200;
# L
L[1]=1;
L[2]=length(Vm);

# "Rx="
Rx[k] = sum((Vp[L[2]-N:L[2]] -Vm[L[2]-N:L[2]])./(Im[L[2]-N:L[2]] ))/(N+1) -1.0;
k=k+1;
}
# BackEMF=(Vp-Vm)-(Im+Ileakm)*Rx;


}
for(j in 2:4){(Rx[j]+Rx[j-1])/2} 
