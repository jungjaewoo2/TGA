//
// eg5_dg535.r: basic control of the srs pulse/delay generator dg535
//
rfile libgpib.so

// initialize gpib devices
if (!exist(ud1))
{
  ud1 = ibdev(0,15,0); // set dg's gpib-address to 15
}

mycode = [ ...
  "CL", ...
  "TZ 1,0", ...
  "DT 2,1,1", ...
  "OM 1,3", ...
  "OA 1,4", ...
  "OO 1,-0.5", ...
  "TM 2", ...
  "SS" ...
  ];
ibwrt(ud1,mycode + "\n");

_i = 0;
while ((_i++)<10)
{
  spinner();
  sleep (1);
  ibwrt(ud1,"SS\n");
}


