//
// file: eg3_read.r
//
// reading from a socket
// For the following example in the same directory open a new konsole, and in
// it type:
// > while true ; do nc6 -v -l 127.0.0.1 -p 5555 < eg1_write.r; done
// This will open a port on which netcat will listen, and upon connection dump
// the content of the file 'eg1_write.r' to the socket.
// The following script is supposed to pick up the content of the file and print it on screen.

sn = "tcp://127.0.0.1:5555";

open(sn);
sleep (0.01); // pause a little, otherwise you will get an empty

nbytes = 10;
while (1)
{
  x = readm(sn, nbytes);
  if (length(x)<nbytes)
  { break; }
  s = char(x);
  printf("%s", s);
}

close(sn);
