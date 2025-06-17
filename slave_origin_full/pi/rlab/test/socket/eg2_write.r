//
// file: eg1_write.r
//
// In this example we write to a socket on the local machine.
// First we initiate a socket using nc(6) by executing in a separate
// konsole window:
//  > while true ; do nc6 -v -l 127.0.0.1 -p 5555; done
// or this one:
//  > while true ; do netcat -v -l 127.0.0.1 -p 5555; done
//
// then we execute this script which attempts to write some textual data
// to the same port.
//
// we have to close the port after we finish writing because this is default
// behavior of netcat

sn = "tcp://127.0.0.1:5555";

open(sn);

writem(sn, [0:255]);

// close (sn);


