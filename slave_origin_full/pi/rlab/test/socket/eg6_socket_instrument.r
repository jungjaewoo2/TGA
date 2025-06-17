//
//
//

sn = "tcp://192.168.0.4:1394";
open (sn);
sleep (0.1);
x = writem (sn, "beep.enable = 1\r\n");
x

close (sn);

