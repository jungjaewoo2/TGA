// RLaB diary file: RS232_debug_with_controller. Opened Tue Jun  5 20:54:55 2018
tic(1);
tic(2);
tic(3);

attr=<<>>;
attr.devname = "ttyUSB";
attr.id_vendor_id = "0403";

if(!exist(r1)){
r1 = usbdevname(attr);
}

r1
toc(1)
opts=<<>>;
opts.data_parity_stop="8N1";
opts.eol="\r";
# opts.debug="char";
opts.speed=115200;
opts.flow_control="x";
// opts.flow_control="n";

serial="serial://" + r1;

if(!open(serial,opts)) {
open(serial, opts);
}

toc(2)
# readcode = ["SENS:DATA:FRES?\r"];
# readcode = ["MS,01\r"];
readcode = ["MS,01,+000.260"];

// number of readings
N = 10;
N = 10000;
N = 1000;
// N = 100;
// use for storing the numerical values
x = zeros(N,2);

# initcode = ["R0\r"];
############################
# pause("serial link ready: Press Enter to continue" );
sleep(4);

//
// read/write query within rlab loop
//
// determine the format of the output. See note above
# datarange = [1:strindex(r[1],"VDC")-1];
# timerange = [strindex(r[1],"VDC")+4:strindex(r[1],"SECS")-1];

# datarange = [1:strindex(r[1],"MS,01")-1: ];
datarange = [7:14];

tic(4)
 for (i in 1:N)
  {
    spinner();
    writem(serial, readcode);
      r = readm(serial);
    printf("%i\t%s\n",i,r); // DEBUG ONLY
      x[i;1] = i;
      x[i;2] = 1e3 * strtod( substr(r, datarange) ); // mm to microns

}
toc(3)
printf("Collected %i datapoints in %f sec\n", N, toc(4));
printf("Average sample time = %f ms\n", toc(4)/N *1000);
dt=toc(4)/N; // in sec
      x[;1] = x[;1] * dt;

plegend();
xlabel( "Time [s]" );
ylabel( "Displacement [um]" );
pltitle( "Laser distance (fake loopback data from RS232)" );

plot(x);
sleep(2);


# pltitle( "Laser distance (previously measured real data)" );
# # read previously measured data from laser 1 using RS232
# zz=readm("laser1_data2_N=1000");
# plot(zz[1:200;])

# sleep(2);
# printf("Peak-to-peak motion swing is: %i [um] \n", max(zz[1:200;2])- min(zz[1:200;2]));
