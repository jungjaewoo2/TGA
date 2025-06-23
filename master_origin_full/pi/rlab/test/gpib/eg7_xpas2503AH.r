//
// eg_tds210.r: ch1 terminal on probe compensation
//
GPIB_CIC_ID = 0;

rfile libgpib.so


// initialize gpib devices
if (!exist(ud1))
{ ud1 = ibdev(0,2,0); }

ibsre (GPIB_CIC_ID,1)


ibqrd(ud1,,"READ=PWR[CHA]\r\n")


