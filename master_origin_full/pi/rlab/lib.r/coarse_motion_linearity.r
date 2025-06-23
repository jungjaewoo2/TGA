coarse_motion_linearity = function(fn,dir,seg,SN) {

// calculate and display nonlineary of motion
// fn -> filename with measured data
// dir -> +1/-1 forward or reverse motion
// seg -> select which segment to analyze
// SN  -> serial number (to label plots)
//
// [azg] Mon Jul 23 17:39:19 PDT 2018
// [azg] Sat Aug 25 00:39:44 PDT 2018
/////////////////////////////////////////////////////

require extract_chunk_with_slope single_slope psprint
global(d)

# Q: Why is this broken? A:(needs global(d))
if(class(fn)== "string"){if(isfile(fn)){ xx=read_ascii(fn);}}

if(class(fn)== "num"){d=fn;}

plstart(1,2);
ylabel("Travel [mm]");
pltitle("Coarse motion vs. constant slope [mm]");
plegend(["Measured travel","ideal"]);
xlabel("Time [s]");
# dir=-1; // forward
if(!exist(dir)){ dir=+1; }

if(dir == +1) { DIR="reverse"; else DIR="forward"; }
if(exist(seg)){ K=seg; else K=8;}

d2=extract_chunk_with_slope(d,K,dir); 
d2=single_slope(d2); 
# d2=6.31418919*d2; # convert to [um]
d2=6.35415*d2; # convert to [um] (IBM number)

t=[1:length(d2)]*40e-5;
t=t[:];
# plot([t,d2]);

ind2= length(d2);
ind1= 1;
t2=t[ind1:ind2];
y=(d2[ind2] -d2[ind1])/(t[ind2]-t[ind1])*t[ind1:ind2] +d2[1];
# y=(1400*6.35415)/(1400/200)*t[ind1:ind2] +d2[1];
# plot(<< [t2,d2/1e3]; [t2,y/1e3]  >>);
# Q: does pgplot support multiple plots stacking?
# plt1=<< [t2,d2/1e3]>>;
plt1=<< [t2,d2/1e3]; [t2,y/1e3]  >>;
subplot(1);
plot(plt1);
# ap1=[t2,d2/1e3];
# write_ascii("V"+V+"_reverse",ap1);
// plot(<< [t2,d2]; [t2,y]  >>);



ylabel("Error [um]");
plegend(["Travel Error [um]"]);
if(!exist(SN)){SN=" ";}
pltitle("Deviation from linear motion of FH-MAA "+SN);
# pltitle("Deviation from linear motion of FH-MAA ");
ind3=[ind1:ind2];
ind3=ind3[:];
xlabel("Time [s]");
 # plot(<< [t2,d2-y]  >>);
 plt2=<< [t2,d2-y]  >>;
plot(plt2);
sleep(1);
plclose();
# psprint("Coarse_motion_linearity_"+DIR+".pdf");
# psprint("Coarse_motion_linearity_V="+V+"_"+DIR+".pdf");
return << plt1=plt1; plt2=plt2 >>;
};
