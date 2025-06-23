coarse_motion_linearity = function(fn,dir,seg) {

// calculate and display nonlineary of motion
// fn -> filename with measured data
// dir -> +1/-1 forward or reverse motion
// seg -> select which segment to analyze
//
// [azg] Mon Jul 23 17:39:19 PDT 2018
/////////////////////////////////////////////////////
# xxx= read_ascii("V200_coarse.dat"); V=num2str(200);
# xxx= read_ascii("V200_w_weight_coarse.dat"); 
V=num2str(200);
# xxx= read_ascii("V400_w_weight_coarse.dat"); V=num2str(400);
# xxx= read_ascii("V600_w_weight_coarse.dat"); V=num2str(600);

require extract_chunk_with_slope single_slope psprint
d=fn;
plstart(1,2);
ylabel("Travel [mm]");
pltitle("Coarse motion vs. constant slope [mm]");
# pltitle("Coarse motion vs. constant slope [mm] (without weight)");
#label("Steps");
plegend(["Measured travel","ideal"]);
xlabel("Time [s]");
/////////////////// NEW ////////////////////////////////////////////// 

# dir=-1; // forward
if(!exist(dir)){ dir=+1; }

if(dir == +1) { DIR="reverse"; else DIR="forward"; }
if(exist(seg)){ K=seg; else K=8;}

# d2=extract_chunk_with_slope(d,K,1); 
# d2=extract_chunk_with_slope(d,K,-1); 
d2=extract_chunk_with_slope(d,K,dir); 
d2=single_slope(d2); 
d2=6.31418919*d2; # convert to [um]
# plot(d2);
#ddd=extract_chunk_with_slope(d,K,1);
#ddd=single_slope(ddd); plot(ddd);

//////////////// END OF NEW //////////////////////////////////////////
#d2=(d[112000:132000]);
t=[1:length(d2)]*40e-5;
t=t[:];
# plot([t,d2]);

ind2= length(d2);
# size(ind2)
ind1= 1;
t2=t[ind1:ind2];
#d3=d2[ind1:ind2];
y=(d2[ind2] -d2[ind1])/(t[ind2]-t[ind1])*t[ind1:ind2] +d2[1];
subplot(1,1);
plot(<< [t2,d2/1e3]; [t2,y/1e3]  >>);
# sleep(1);
# ap1=[t2,d2/1e3];
# write_ascii("V"+V+"_reverse",ap1);
// plot(<< [t2,d2]; [t2,y]  >>);

subplot(1,2);

ylabel("Error [um]");
plegend(["Travel Error [um]"]);
pltitle("Deviation from linear motion of FH-MAA (Velocity="+V+")");
# ind3=[ind1:ind2] -ind1 +1;
ind3=[ind1:ind2];
#size(d3)
#       19319             1  
# size(t2)
#      19319             1  
# size(ind3)
#          1         19319  
ind3=ind3[:];
#plot(<< [ind3,d3-y]  >>);
xlabel("Time [s]");
 plot(<< [t2,d2-y]  >>);
sleep(10); 
plclose();
// psprint("Coarse_motion_linearity_V="+V+"_"+DIR+".pdf");
#psprint("Coarse_motion_linearity_V="+V+"_"+DIR+"unloaded.pdf");

#diary();

};
