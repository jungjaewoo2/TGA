//
//
//

sn = "tcp://127.0.0.1:4242";
open(sn);
sleep (0.01);

x = readm(sn);
if (isempty(x))
{ x = readm(sn); }

x

    