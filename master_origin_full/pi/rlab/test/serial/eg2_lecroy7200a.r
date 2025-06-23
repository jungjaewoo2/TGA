//
// eg_lecroy7200a.r: short the terminals for voltage measurement
//

// initialize serial port for communication with
// lecroy oscilloscope:
//  19200 baud rate, 8N1 with hardware control, leave default raw communication,
//  and write all communication from the port as characters (debug)
serial = "/dev/ttyS0";
open(serial, "8N1", 19200, "x");

firstquery=[...
    "ASET;DATE?",   ...
    "COMM_FORMAT?"  ...
        ];

writem(serial, firstquery)
sleep(0.5);



#writem(serial,["T1:waveform? dat1"]);
#x = readm(serial);
                    
x = readm(serial, "COMM_FORMAT?");

x








