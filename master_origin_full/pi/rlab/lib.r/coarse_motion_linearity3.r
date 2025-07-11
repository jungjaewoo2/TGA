coarse_motion_linearity3 = function(fn,SN) {

// calculate and display nonlineary of motion
// fn -> filename with measured data
// SN  -> serial number (to label plots)
//
// [azg] Mon Jul 23 17:39:19 PDT 2018
// [azg] Sat Aug 25 00:39:44 PDT 2018
// [azg] Thu Mar 21 16:25:27 KST 2019 (subplot fiasco)
/////////////////////////////////////////////////////

require extract_slope_indexes single_slope psprint
global(d)

if(class(fn)== "string"){if(isfile(fn)){ xx=read_ascii(fn);}}

if(class(fn)== "num"){d=fn;}

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


# reverse"
ind=extract_slope_indexes(d);
d2=single_slope(d[1:ind[2]]); 
# [azg] Fri Nov 30 18:58:13 KST 2018
# re-caculated convesion coefficient of laser ADC count to um:
# 0.00546822157 * 1000
d2=5.46822157*d2; # convert to [um]

t=[1:length(d2)]*40e-5;
t=t[:];

ind2= length(d2);
ind1= 1;
t2=t[ind1:ind2];
y=(d2[ind2] -d2[ind1])/(t[ind2]-t[ind1])*t[ind1:ind2] +d2[1];


# "forward"
######################################################################
d3=single_slope(d[ind[1]:ind[3]]); 
d3=5.46822157*d3; # convert to [um]

t3=[1:length(d3)]*40e-5;
t3=t3[:];

ind3= length(d3);
ind1= 1;
t3=t3[ind1:ind3];
y3=(d3[ind3] -d3[ind1])/(t3[ind3]-t3[ind1])*t3[ind1:ind3] +d3[1];
######################################################################

######################################################################
d4=single_slope(d[ind[2]:ind[4]]); 
d4=5.46822157*d4; # convert to [um]

t4=[1:length(d4)]*40e-5;
t4=t4[:];

ind4= length(d4);
ind1= 1;
t4=t4[ind1:ind4];
y4=(d4[ind4] -d4[ind1])/(t4[ind4]-t4[ind1])*t4[ind1:ind4] +d4[1];
######################################################################

######################################################################
d5=single_slope(d[ind[3]:length(d)]); 
d5=5.46822157*d5; # convert to [um]

t5=[1:length(d5)]*40e-5;
t5=t5[:];

ind5= length(d5);
ind1= 1;
t5=t5[ind1:ind5];
y5=(d5[ind5] -d5[ind1])/(t5[ind5]-t5[ind1])*t5[ind1:ind5] +d5[1];
######################################################################


plegend(["FWD Travel Error [um]"]);
if(!exist(SN)){SN=" ";}
pltitle("Deviation from linear motion of FH-MAA "+SN);
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

ind3=[ind1:ind2];
ind3=ind3[:];
plimits(-0.5,8,,,,);
plt2=<< [(d2-d2[1])/1000 ,d2-y]  >>;
plot(plt2);
plegend(["REV Travel Error [um]"]);
if(!exist(SN)){SN=" ";}
pltitle("Deviation from linear motion of FH-MAA "+SN);
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

plimits(-0.5,8,,,,);
plt3=<< [(d3-d2[1])/1000 ,d3-y3]  >>;
plot(plt3);
plegend(["FWD Travel Error [um]"]);
if(!exist(SN)){SN=" ";}
pltitle("Deviation from linear motion of FH-MAA "+SN);
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

plimits(-0.5,8,,,,);
plt4=<< [(d4-d2[1])/1000 ,d4-y4]  >>;
plot(plt4);
plegend(["REV Travel Error [um]"]);
if(!exist(SN)){SN=" ";}
pltitle("Deviation from linear motion of FH-MAA "+SN);
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

plimits(-0.5,8,,,,);
plt5=<< [(d5-d2[1])/1000 ,d5-y5]  >>;
plot(plt5);
sleep(5);
# pause();
plclose();
return << plt3=plt3; plt2=plt2; plt4=plt4; plt5=plt5 >>;


};
