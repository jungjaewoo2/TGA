//
// eg_xech.r
//
//
//

//url = "serial:///dev/ttyUSB0";
url = "serial:///dev/ttyS0";

if(!exist(url))
{ stop("'url' required to continue!"); }

url_opts = <<>>;
url_opts.data_format  = "8N1";
url_opts.speed    = 9600;
url_opts.flow_control = "none";
url_opts.debug = "hex";
url_opts.raw   = 1;
url_opts.eos = "";

close (url);
open  (url, url_opts);
sleep (0.5);


xwrite = function( arg1 )
{
  global (url);
  xs = strsplt(arg1);
  for (x in xs)
  {
    writem(url,x);
    sleep (0.05);
  }
  return 1;
};

x_getfreq = function( )
{
  global (url);
  xwrite( " " );
  x=readm(url,7);

  if (x.n == 7)
  { x = x[2:7]; }

  rval = 0.01*(256*x[3]+x[4]);
  return rval;
};

x_setfreq = function( arg1 )
{
  f = text(arg1,"%.2fF");
  xwrite(f);
  return 1;
};

x_setvolt = function( arg1 )
{
  // check the bounds
  if (arg1<85)
  { arg1 = 85; }
  if (arg1>250)
  { arg1 = 250; }
  // set it
  v = text(arg1,"%.3fV");
  xwrite(v);

  return 1;
};

x_getvolt = function( )
{
  global (url);
  xwrite( " " );
  x=readm(url,7);

  if (x.n == 7)
  { x = x[2:7]; }

  rval = 0.1*(256*x[5]+x[6]);
  return rval;
};

x_setocp = function( arg1 )
{
  // check the bounds
  if (arg1<0.1)
  { arg1 = 0.1; }
  if (arg1>1.5)
  { arg1 = 1.5; }

  o = text(arg1,"%.3fP");
  xwrite(o);
  return 1;
};

x_getocp = function( )
{
  global (url);
  xwrite( " " );
  x=readm(url,7);

  if (x.n == 7)
  { x = x[2:7]; }

  rval = 0.001*(256*x[1]+x[2]);
  return rval;
};


x_stby = function()
{
  xwrite("N");
  return 1;
};

x_oper = function()
{
  xwrite("O");
  return 1;
};

x_toggle = function()
{
  xwrite("Y");
  return 1;
};

x_setvolt(230);
x_setfreq(50);
x_setocp(0.1);


stop()

v_Range = [85:250];
f_Range = [45:65];

drng(1, v_Range, ones(v_Range));
drng(2, f_Range, ones(f_Range));

for (i in 1:20)
{
  drng(1);
  v = drand();

  drng(2);
  f = drand();

  x_setvolt(v);
  sleep (0.3);

  x_setfreq(f);
  sleep (0.3);

  writem(url," "); sleep (0.2); x=readm(url,7);
  if (x[1]==32)
  { x = x[2:7]; }
  printf ("%3g = %3i %3i -> %g\n", v, x[5], x[6], 0.1*(256*x[5]+x[6]));
}


stop()

fqs = [45:65]';
fqs = [45,46,49,50,51]' + 0.25;
fqs = shuffle(fqs);


for (i in 87:90)
{
  char(i)
  writem(url, char(i));
  readm (url,16)
  sleep(1);
}


stop()


for (f in fqs)
{
  x_setfreq(f);
  sleep (0.5);
  writem(url," "); x=readm(url,7);
  if (x.n == 5)
  {
    x = [x[1:2], 18, x[3:5]];
  }
  printf ("%3g = %3i %3i -> %g\n", f, x[3], x[4], 0.01*(256*x[3]+x[4]));

}






