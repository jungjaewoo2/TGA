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

sn = "tcp://127.0.0.1:5555";

if (!exist(j))
{ j=0; }

open(sn);
msg = text(j+[1:10], "Question No. %2.0f: ") + blank(1,10,"Is anybody listening on the other end?\n");
writem(sn, msg);
close (sn);

j=j+10;

