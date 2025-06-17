//
// eg_arduino:
//  requires MATLAB serial port server to run on arduino
//
// see arduino.cc, and search for matlab on it.
//

s = "serial:///dev/ttyUSB0";
open(s, "8N1", 9600, "xon|xoff");

while (1)
{
  writem(s, "2l0\r");

  sleep( 10 * urandom() );

  writem(s, "2l1\r");

  sleep( 20 * urandom() );
}

