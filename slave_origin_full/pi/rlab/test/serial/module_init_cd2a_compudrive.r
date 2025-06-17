
//
//
static(URL, URL_OPTS);

# CHIP_ADD = "/dev/ttyS4";
CHIP_ADD = "/dev/ttyUSB1";
if(!exist(CHIP_ADD))
{
  x = last(reads("| sudo tail -10000 /var/log/messages | grep ttyUSB | grep attached | grep ttyUSB"));
  ix = strindex(x,"tty");
  if (ix)
  {
    CHIP_ADD = "/dev/" + substr(x,ix:ix+9);
  }
}
"\nTrying to talk to device at " + CHIP_ADD

URL_OPTS = <<>>;
URL_OPTS.data_format  = "8N1";
URL_OPTS.speed    = 2400;
URL_OPTS.flow_control = "xon|xoff";
URL_OPTS.eos = "\n";
URL = "serial://" + CHIP_ADD;

if (!open(URL,URL_OPTS))
{ sleep(3); }

url = URL;


static(WAIT_TIME_1CH,EOS,STX,ETX,ACK,NAK,CAN,DEBUG);
WAIT_TIME_1CH = 0.05;
EOS = "\r\n";
STX = char(2);
ETX = char(3);
ACK = char(6);
NAK = char(21);
CAN = char(24);
DEBUG=0;

got_two = 0;
printf("Please press \"Remote\" Button on CD2A Compudrive Console!\n");
while (1)
{
  spinner();
  x = readm(url,1);
  if(class(x)!="num")
  {
    sleep (0.1);
    continue;
  }
  if (isempty(x))
  {
    sleep (0.1);
    continue;    
  }
  
  if (x==ascii(ACK))
  { got_two = 1; }
  if (x==ascii(CAN) && got_two==1)
  {
    got_two = 2;
    break;
  }
  sleep (0.1);
}
printf("Received <ACK><CAN> from Console. Assuming system on-line!\n");

static(_send_cmd, _send_par,_set_position);
_send_cmd = function(cmd)
{
  if (!exist(cmd))
  { return 1;}

  if (class(cmd)!="string")
  { return 2; }

  sval = CAN + cmd + ETX + EOS;

  x = ascii(sval);
  if (DEBUG)
  { x }

  for (i in range(x))
  {
    y = x[i];
    writem(URL, y);
    sleep(WAIT_TIME_1CH);
  }  

  return 0;
};

_send_par = function(par)
{
  if (!exist(par))
  { return 1;}
  
  if (class(par)!="string")
  { return 2; }
  
  sval = STX + par + ETX + EOS;

  x = ascii(sval);
  if (DEBUG)
  { x }
  
  for (i in range(x))
  {
    y = x[i];
    writem(URL, y);
    sleep(WAIT_TIME_1CH);
  }
  
  return 0;
};

_set_position = function ( x )
{
  par = "SE" + text(x,"%8.2f");
  _send_par (par);
  _send_cmd ("P");

  i = 0;
  s_new = "";
  s_old = s_new;
  while (1)
  {
    i = readm(URL,14);
    j = find (i>31);
    if(!isempty(j))
    {s_new = char(i[j]);}
    if (s_new != s_old)
    {
      s_old = s_new;
      if (DEBUG)
      {
        printf("%s\n",s_old);
      else
        spinner();
      }
    }
    if (any(i==42))
    { break; }
  }
  return 0;
};

cdr_debug = function (x)
{
  if (!exist(x))
  { return DEBUG;}

  DEBUG = (x>0);
  return DEBUG;
};


cdr = <<>>;
cdr.cmd = function(cmd)
{
  return _send_cmd(cmd);
};

cdr.par = function(par)
{
  return _send_par(par);
};

cdr.halt = function()
{
  return _send_cmd("H");
};

cdr.wlen = function ( x )
{
  return _set_position(x);
};

cdr.scan = <<>>;
cdr.scan.init = function (a,b,c)
{
  // set scan mode to burst  
  _send_par ("TYB");
  
  // set to single scan
  _send_par ("NS1");

  // set start scan:
  par = "ST" + text(a,"%8.2f");
  _send_par (par);

  // set end scan:
  par = "EN" + text(b,"%8.2f");
  _send_par (par);

  // set scan increment  
  par = "BI" + text(c,"%8.2f");
  _send_par (par);

  // enable trigger scan
  _send_cmd("T");

  // wait until spectrometer positions itself
  i = 0;
  while ( all(i!=ascii("S")) )
  {
    i = readm(URL,14); 
    j = find (i>31);
    if(!isempty(j))
    {
      rval = char(i[j]);
    else
      rval = "";
    }
  } 
  
  return strtod(substr(rval,3:strlen(rval)));
};

cdr.scan.next = function()
{
  s = 1;
  _send_cmd("E");

  i = 0;
  while (all(i!=ascii("B")) && all(i!=ascii("E")))
  {
    i = readm(URL,14);
  }

  if (any(i==ascii("E")))
  { s = -1; }
  
  j = find (i>31);
  if(!isempty(j))
  {
    rval = char(i[j]);
  else
    rval = "";
  }
  
  return s * strtod(substr(rval,3:strlen(rval)));
};
