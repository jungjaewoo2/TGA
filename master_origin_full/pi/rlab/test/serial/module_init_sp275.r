
//
//
static(URL, URL_OPTS, _SER_N_TOC, _SER_N_TIMEOUT);

_SER_N_TOC = 30;        // use toc() timer No. 29
_SER_N_TIMEOUT = 0.25;  // how much to wait in seconds between the attempts

if(!exist(CHIP_ADD))
{
  x = last(reads("| sudo tail -10000 /var/log/messages | grep ttyUSB | grep attached | grep ttyUSB"));
  ix = strindex(x,"tty");
  if (ix)
  {
    CHIP_ADD = "/dev/" + substr(x,ix:ix+9);
  }
}
"Trying to talk to device at " + CHIP_ADD

URL_OPTS = <<>>;
URL_OPTS.data_format  = "8N2";
URL_OPTS.speed    = 9600;
URL_OPTS.flow_control = "xon|xoff";
URL_OPTS.eos = "\n";
URL = "serial://" + CHIP_ADD;

if (!open(URL,URL_OPTS))
{ sleep(5); }

url = URL;


