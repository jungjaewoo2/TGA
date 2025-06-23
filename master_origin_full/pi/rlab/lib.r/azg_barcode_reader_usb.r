//
// find which port the device with properties below is attached to
//
attr=<<>>;
attr.devname = "ttyACM";
# at home works (Symbol?) #attr.id_vendor_id = "0483";
# idVendor=1eab, idProduct=0c06
attr.id_vendor_id = "1eab";

# attr.id_vendor_id = "0403";
### [azg]seems to work without S/N
### attr.id_serial_short = "20383959424B";
# attr.id_serial_short = "A600CJ5X";


# for( i in 1:5){

if(!exist(r1)){
r1 = usbdevname(attr);
# }



colors("red");
printf("Waiting for S/N...");

r=getline(r1,16);

colors("green");
printf("\rS/N:              \t%s",r);
# printf("\rS/N:              \t%s\n",r);
colors();

}
