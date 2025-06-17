//
// eg_lecroy7200a.r: short the terminals for voltage measurement
//
rfile libgpib.so

// initialize gpib devices
if (!exist(ud1))
{
  ud1 = ibdev(0,0,0); // gpib main board, here gpib-usb-hs, set to 0
}
if (!exist(ud2))
{
  ud2 = ibdev(0,8,0); // set lecroy7200a's gpib-address to 8
}

firstquery=[...
    "ASET;DATE?",   ...
    "COMM_FORMAT?"  ...
        ];

for (i in 1:length(firstquery))
{
  r = ibqrd(ud2,firstquery[i]);
  printf("lecroy: %s\n", r);
}
ibwrt(ud2, "comm_format ind0,word,hex");



ibwrt (ud2,["T1:waveform? dat1"]);
x = ibrd(ud2,1024);

x








