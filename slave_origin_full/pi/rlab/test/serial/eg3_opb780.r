//
// eg_opb780.r
//
//
//
// Optek OPB780 Evaluation Kit
//
// Press Enter key to start.
//
//
//
// Optek OPB780 Evaluation Kit HW: V1.1, SW: V1.2
//
// Valid Commands are:
//
// R (ed)
//
// G (reen)
//
// B (lue)
//
// C (lear)
//
// A (ll)
//
// I (ntegration time)
//
// L (ed Enable)
//
// Enter Command (R G B C A I L):
//
//
s = "serial:///dev/ttyUSB0";

red_value = function( s )
{
  local(x,y,y1,y2);

  writem(s, "R\r\r");

  y = "";
  x = "";
  while (isempty(findstr(x,":")))
  {
    sleep (0.1);
    x = readm(s);
    y=y+x;
  }

  y1 = findstr(y,"Red=") + 4;
  y2 = findstr(y,"Press") - 2;
  rval = strtod(substr(y, y1:y2));

  return rval;
};

// open(s, "8N1", 9600, "h", ,"hex");
open(s, "8N1", 9600, "h");

data_t = [];

for (i in 1:100)
{
  spinner();

  writem(s, "R\r\r");

  data_t = [data_t; red_value(s)];
}











