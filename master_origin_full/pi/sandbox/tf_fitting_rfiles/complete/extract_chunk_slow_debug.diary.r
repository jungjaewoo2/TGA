// RLaB diary file: extract_chunk_slow_debug.diary.r. Opened Mon May 14 01:12:25 2018


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
ind[find((diff(ind) >5000) && (diff(ind) <6000 ))]


 # plot([4e-5*[1:length(d)]', [d;1]]);
 plot([4e-5*[ind[ind2[1]]:ind[ind2[4]]+7000]', [d[ind[ind2[1]]:ind[ind2[4]]+7000 ;1]]]);

ind[ind2]*4e-5
