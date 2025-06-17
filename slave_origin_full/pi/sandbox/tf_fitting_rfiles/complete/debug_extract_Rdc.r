require extract_chunk
tic();
xxx=read_ascii("no_timestamp_production_vol=0.2.dat");
" file read"
toc()

k=1;
####################################
for(j in 1:200){
j
tic();
	dc1=[];
	chunk = j;
	idx= extract_chunk(d[;1],  chunk, 100 );
	if (  (diff (idx) > 5000 )  && (diff (idx) < 6000 )  ) {
		        
		out[k;]=[j,idx];
		size([j,idx])

		k=k+1;
	}
toc()
}
