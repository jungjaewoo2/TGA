// RLaB diary file: restore_stuck_elevator.diary.r. Opened Fri Jul  1 12:37:45 2022

// stepper controller ////////////////////////////////////////////////

attr=<<>>;
attr.devname = "ttyUSB";
attr.id_vendor_id = "10c4";

if(!exist(r2)){
	r2 = usbdevname(attr);
}
 


opts=<<>>;
opts.data_parity_stop="8N1";
opts.eol="\r";
opts.speed=9600;
opts.flow_control="n";

stepper="serial://" + r2;

if(!open(stepper,opts)) {
open(stepper, opts);
}


fwd=["/1P2000R"];
rev=["/1D2000R"];
# writem(stepper, rev);
# writem(stepper, fwd);
# writem(stepper, rev);
# writem(stepper, rev);
# diary();


fwd2=["/1V800P20000R"];
rev2=["/1V800D2000R"];
writem(stepper, rev2);
#writem(stepper, fwd2);
