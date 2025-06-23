coarse_motion_linearity5_tga = function(flatsidx,SN) {

// calculate and display nonlineary of motion
// using (global) "d" and "flatsidx" calculated previously by "extract_postions_tga"
// SN  -> serial number (to label plots)
//
// STATUS: New untested
// [azg] Fri Nov  8 00:32:55 PST 2019
/////////////////////////////////////////////////////
require extract_positions_tga plplot

global(d, _rlab_config)
# plot=plplot;

if (_rlab_config.plot_support=="pgplot") {
	    {  plstart(1,4);
               subplot(1,1); # pgplot syntax 
	    }
     else 
     	{     # plplot syntax
	        plwins(1);
	        plwin(1); 
	 }

   }

###########################################################################
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$



# jnk4=zeros(d); jnk4[100*flatsidx]=ones(100*flatsidx);
jnk4=zeros(d); jnk4[100*flatsidx]=ones(100*flatsidx);
# NOTE: extract_postions function subsampled data d 100x
# plclose();plot(d);plot(jnk4);
# plot(d);plot(jnk4);
endi=maxi(flatsidx); 
# printf("endi=%i \n\n",endi);
# sleep(2);

d2=d[100*flatsidx[endi-9]:100*flatsidx[endi-8]];
d3=d[100*flatsidx[endi-7]:100*flatsidx[endi-6]];
d4=d[100*flatsidx[endi-5]:100*flatsidx[endi-4]];
d5=d[100*flatsidx[endi-3]:100*flatsidx[endi-2]];

#######################################
d2=5.46822157*d2; # convert to [um]

t2=[1:length(d2)]*40e-5;
t2=t2[:];

ind2= int(0.9*length(d2));
ind1= int(0.1*length(d2));
# t2=t2[ind1:ind2];
t2=t2[1:length(d2)];
y2=(d2[ind2] -d2[ind1])/(t2[ind2]-t2[ind1])*t2[1:length(d2)] +d2[1];
# plt2=<< [(d2-d2[1])/1000 ,d2-y2]  >>;
 trunc=200;
p2= [(d2-d2[1])/1000 ,d2-y2][1:length(d2)-trunc;] ;
plt2 = << p2 >>;

# plclose(); plot(plt2);
# plot(plt2);
# sleep(2);

######################################################################

d5=5.46822157*d5; # convert to [um]

t5=[1:length(d5)]*40e-5;
t5=t5[:];

ind5= int(0.9*length(d5));
ind1= int(0.1*length(d5));
# t5=t5[ind1:ind5];
t5=t5[1:length(d5)];
y5=(d5[ind5] -d5[ind1])/(t5[ind5]-t5[ind1])*t5[1:length(d5)] +d5[1];
######################################################################
# plt5=<< [(d5-d5[1])/1000 ,d5-y5]  >>;
 trunc=200;
p5= [(d5-d5[1])/1000 ,d5-y5][1:length(d5)-trunc;] ;
plt5 = << p5 >>;

# plclose();plot(plt5);
# plot(plt5);
# sleep(2);



######################################################################
d3=5.46822157*d3; # convert to [um]

t3=[1:length(d3)]*40e-5;
t3=t3[:];

######################################################################

d4=5.46822157*d4; # convert to [um]

t4=[1:length(d4)]*40e-5;
t4=t4[:];

######################################################################

# revised 10-90% slope
ind3= int(0.9*length(d3));
ind1= int(0.1*length(d3));
t3=t3[1:length(d3)];
y3=(d3[ind3] -d3[ind1])/(t3[ind3]-t3[ind1])*t3[1:length(d3)] +d3[1];
# plt3=<< [(d3-d3[1])/1000 ,d3-y3]  >>;
 trunc=200;
p3= [(d3-d3[1])/1000 ,d3-y3][1:length(d3)-trunc;] ;
plt3 = << p3 >>;

# plclose();plot(plt3);
# plot(plt3);
# sleep(2);
######################################################################

# revised 10-90% slope
ind4= int(0.9*length(d4));
ind1= int(0.1*length(d4));
t4=t4[1:length(d4)];
y4=(d4[ind4] -d4[ind1])/(t4[ind4]-t4[ind1])*t4[1:length(d4)] +d4[1];
# plt4=<< [(d4-d4[1])/1000 ,d4-y4]  >>;
 trunc=200;
p4= [(d4-d4[1])/1000 ,d4-y4][1:length(d4)-trunc;] ;
plt4 = << p4 >>;

# plclose();plot(plt4);
# plot(plt4);
# sleep(2);

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
###########################################################################


# ylabel("Error [um]");
plegend(["FWD Travel Error [um]"]);
if(!exist(SN)){SN=" ";}
pltitle("Deviation from linear motion. SN= "+SN);
if (_rlab_config.plot_support=="pgplot")
	    {
	    {
               xlabel("Travel [mm]"); # pgplot syntax 
               ylabel("Error [um]");
	    }
     else
     {
               plxlabel("Travel [mm]"); 
               plylabel("Error [um]");
               plmultiplot([0.0,0.99,0.80,0.95]);
     }

   }

plot(plt2);
# ylabel("Error [um]");
plegend(["REV Travel Error [um]"]);
if(!exist(SN)){SN=" ";}
pltitle("Deviation from linear motion. SN= "+SN);
if (_rlab_config.plot_support=="pgplot")
	    {
	    {
               xlabel("Travel [mm]"); # pgplot syntax 
               ylabel("Error [um]");
	    }
     else
     {
               plxlabel("Travel [mm]"); # plplot syntax
               plylabel("Error [um]"); 
	       plmultiplot([0.0,0.99,0.55,0.7]);
     }
   }

plot(plt3);
plegend(["FWD Travel Error [um]"]);
if(!exist(SN)){SN=" ";}
pltitle("Deviation from linear motion. SN= "+SN);
if (_rlab_config.plot_support=="pgplot")
	    {
	    {
               xlabel("Travel [mm]"); # pgplot syntax 
               ylabel("Error [um]");
	    }
     else
     {
               plxlabel("Travel [mm]"); # plplot syntax
               plylabel("Error [um]"); 
	       plmultiplot([0.0,0.99,0.30,0.45]);
     }
   }

plot(plt4);
# ylabel("Error [um]");
plegend(["REV Travel Error [um]"]);
if(!exist(SN)){SN=" ";}
pltitle("Deviation from linear motion. SN= "+SN);
if (_rlab_config.plot_support=="pgplot")
	    {
	    {
               xlabel("Travel [mm]"); # pgplot syntax 
               ylabel("Error [um]");
	    }
     else
     {
               plxlabel("Travel [mm]"); # plplot syntax
               plylabel("Error [um]");
	       plmultiplot([0.0,0.99,0.05,0.2]);
     }
   }

#plimits(-0.5,8,,,,);
plot(plt5);
sleep(5);
# pause();
plclose();
return << plt3=plt3; plt2=plt2; plt4=plt4; plt5=plt5 >>;


};
