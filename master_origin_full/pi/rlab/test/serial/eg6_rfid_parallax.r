//
// eg6_rfid_parallax:
//    https://www.parallax.com/product/28340

static(attr, URL_OPTS, RFID_READER, RFID_TIMER, RFID_TIMEOUT);

attr=<<>>;
attr.devname = "tty";
attr.id_serial_short = "A704FCQI"; // put your serial number here dummy
url = usbdevname(attr);
RFID_READER = url;
printf("Trying to talk to device at %s\n", RFID_READER);

URL_OPTS = <<>>;
URL_OPTS.data_format  = "8N1";
URL_OPTS.speed = 2400;
RFID_TIMER = 10;
RFID_TIMEOUT = 10;

if (!open(RFID_READER,URL_OPTS))
{ sleep (3); }

read_rfid_tag = function()
{
  writem(RFID_READER,<<dtr=0>>);
  sleep(0.1);
  writem(RFID_READER,<<dtr=1>>);
  sleep(0.1);
  writem(RFID_READER,<<dtr=0>>);
  sleep(0.1);
  writem(RFID_READER,<<dtr=1>>);
  sleep(0.1);
  writem(RFID_READER,<<dtr=0>>);
  sleep(0.1);
  writem(RFID_READER,<<dtr=1>>);
  sleep(0.1);
  tic(RFID_TIMER);
  for (i in [1,2])
  {
    x  = [];
    while ((toc(RFID_TIMER)<RFID_TIMEOUT) && (length(x)<12))
    {
      x=[x,readm(RFID_READER,12)];
    }
  }
  if (length(x)==12)
  {
    _i = find(x!=0x0a && x!=0x0d);
    c = int(strtod("0x" + char(x[_i[3:length(_i)]])));
    printf("obtained reading = %i\n", c);
  }
  sleep(0.3);
  writem(RFID_READER,<<dtr=0>>);
  return c;
};




while (1)
{
  r = read_rfid_tag();
  if (!exist(r))
  { continue; }

  q = prompt("Accept reading [" + text(r,"%i") + "]", "y",["n", "y"]);

  if (q=="y")
  { break; }
}
