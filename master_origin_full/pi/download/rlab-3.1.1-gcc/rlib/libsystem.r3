// file: rlib/libsystem.r
//
// This file is a part of RLaB + rlabplus
// Copyright (C) 2013-2014  Marijan Kostrun
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
//
// See the file ../COPYING
//
// some useful functions that depend on system- and reads- call
// to linux system
static(INIT,WHICH);

if (exist(INIT))
{ EOF }

if (!exist(WHICH))
{ WHICH=reads("|which which"); }
if (strlen(WHICH)<1)
{
  INIT = 1;
  EOF
}

static(EMAIL,TOUCH,VI,KONSOLE,TEST,EMACS,DIFF,XTERM,TOP,ENV,FIND,UDEVADM,GREP,ECHO);

// check for 'xterm'
XTERM = reads("|" + WHICH + " xterm 2>/dev/null");
if (isempty(XTERM))
{ clear(XTERM); }
if (strlen(XTERM)<strlen("xterm"))
{ clear(XTERM);}

// check for 'diff'
DIFF = reads("|" + WHICH + " diff 2>/dev/null");
if (isempty(DIFF))
{ clear(DIFF); }
if (strlen(DIFF)<strlen("diff"))
{ clear(DIFF);}

// check for 'mail'
EMAIL = reads("|" + WHICH + " mail 2>/dev/null");
if (isempty(EMAIL))
{ clear(EMAIL); }
if (strlen(EMAIL)<strlen("EMAIL"))
{ clear(EMAIL);}

// check for 'touch'
TOUCH = reads("|" + WHICH + " touch 2>/dev/null");
if (isempty(TOUCH))
{ clear(TOUCH); }
if (strlen(TOUCH)<strlen("TOUCH"))
{ clear(TOUCH);}

// check for 'vi'
VI = reads("|" + WHICH + " vi 2>/dev/null");
if (isempty(VI))
{ clear(VI); }
if (strlen(VI)<strlen("VI"))
{ clear(VI);}

// check for 'emacs'
EMACS = reads("|" + WHICH + " emacs 2>/dev/null");
if (isempty(EMACS))
{ clear(EMACS); }
if (strlen(EMACS)<strlen("EMACS"))
{ clear(EMACS);}

// check for 'konsole'
KONSOLE = reads("|" + WHICH + " konsole 2>/dev/null");
if (isempty(KONSOLE))
{ clear(KONSOLE); }
if (strlen(KONSOLE)<strlen("KONSOLE"))
{ clear(KONSOLE);}

// check for 'test'
TEST = reads("|" + WHICH + " test 2>/dev/null");
if (isempty(TEST))
{ clear(TEST); }
if (strlen(TEST)<strlen("TEST"))
{ clear(TEST);}

// check for 'top'
TOP = reads("|" + WHICH + " top 2>/dev/null");
if (isempty(TOP))
{ clear(TOP); }
if (strlen(TOP)<strlen("TOP"))
{ clear(TOP);}

// check for 'env'
ENV = reads("|" + WHICH + " env 2>/dev/null");
if (isempty(ENV))
{ clear(ENV); }
if (strlen(ENV)<strlen("ENV"))
{ clear(ENV);}

// check for 'find'
FIND = reads("|" + WHICH + " find 2>/dev/null");
if (isempty(FIND))
{ clear(FIND); }
if (strlen(FIND)<strlen("FIND"))
{ clear(FIND);}

// check for 'udevadm'
UDEVADM = reads("|" + WHICH + " udevadm 2>/dev/null");
if (isempty(UDEVADM))
{ UDEVADM = reads("|" + WHICH + " /sbin/udevadm 2>/dev/null"); }
if (isempty(UDEVADM))
{ clear(UDEVADM); }
if (strlen(UDEVADM)<strlen("UDEVADM"))
{ clear(UDEVADM);}

// check for 'grep'
GREP = reads("|" + WHICH + " grep 2>/dev/null");
if (isempty(GREP))
{ clear(GREP); }
if (strlen(GREP)<strlen("GREP"))
{ clear(GREP);}

static(pty_init);
pty_init = function(titlestring, tty)
{
  // figure out the title string
  if (!exist(titlestring))
  { titlestring="RLaB@"+hostname()+":stderr"; }

  // did the user specify type of terminal/konsole
  if (!exist(tty))
  {
    if (exist(XTERM))
    {
      tty = XTERM;
    } else { if (exist(KONSOLE))
    {
      tty = KONSOLE;
    } else {
      tty = "/dev/null"
    }}
  }

  nextterm = "";
  if (tty != "/dev/null")
  {
    // check with system what is the name of the next konsole/terminal
    nextterm=openpty();

    // initialize terminal/konsole
    if (strindex(tty,"xterm"))
    {
      cmdstring = "("+XTERM+" -C -n Rlab -mesg -geometry 200+0+800 -T \""+titlestring+"\")&";
    } else {
      cmdstring = "("+KONSOLE+" -caption \"" + titlestring + "\")&";
    }
    system(cmdstring);
  }

  // report the name of the character device
  return nextterm;
};

rconsole = function(titlestring, tty)
{
  global(_rlab_config);

  // does stderr exist in _rlab_config and is it a character special file?
  if (exist(_rlab_config.stderr))
  {
    if(!ischarspec(_rlab_config.stderr))
    {
      // user must have closed the character device: close it and remove
      // it from the list of open files
      release (_rlab_config);
      clear (_rlab_config.stderr);
      protect (_rlab_config);
    }
  }

  // do we need to initiale character device
  if( !exist(_rlab_config.stderr) )
  {
    // initialize the stderr xterm/console
    ncurr_ptys = ls("/dev/pts/").n;
    release (_rlab_config);
    _rlab_config.stderr = pty_init(titlestring, tty);
    protect (_rlab_config);
    // go to sleep while the system sets up the terminals
    if(!ischarspec(_rlab_config.stderr))
    { sleep(0.1); }
    sleep (1);  // sleep one more for fortran solvers waiting for it
    _j=0;
    while(!fprintf(_rlab_config.stderr,"\n"))
    {
      sleep(0.1);
      if ((_j++) > 50)
      {
        printf("rlabplus: STDERR terminal failed to open!");
        return blank(0,0);
      }
    }
    // we are using character device and not file.
    fprintf(_rlab_config.stderr, "rlabplus: STDERR terminal");
    fprintf(_rlab_config.stderr, " on %s succesfully initialized!\n", _rlab_config.stderr);
    fprintf(_rlab_config.stderr, "rlabplus: Part of the rlabplus project.");
    fprintf(_rlab_config.stderr, " (c) Marijan Kostrun, 2004-2014.\n");
    fprintf(_rlab_config.stderr, "rlabplus: To access STDERR use:\n");
    fprintf(_rlab_config.stderr, "rlabplus:   (i) '_rlab_config.stderr'");
    fprintf(_rlab_config.stderr, " system variable; or,\n");
    fprintf(_rlab_config.stderr, "rlabplus:   (ii) output of  stderr()  function.\n");
    fprintf(_rlab_config.stderr, "rlabplus: as the first argument in  fprintf()\n");
    fprintf(_rlab_config.stderr, "rlabplus: To close STDERR, bring into the focus then press");
    fprintf(_rlab_config.stderr, " Ctrl+D, or press Ctrl-C in RLaB console.\n");
  }
  else
  {
    // stderr already initialized
    fprintf(_rlab_config.stderr, "rlabplus: STDERR X-terminal");
    fprintf(_rlab_config.stderr, " on %s already initialized!\n", _rlab_config.stderr);
  }

  return _rlab_config.stderr;
};
stderr = rconsole;

mail = function(msgaddress, msgsubject, message)
{
  if (!exist(EMAIL))
  { error("mail: System command not available. Cannot continue !"); }

  emailfs = "|"+EMAIL+" -v -s \"" + msgsubject + "\" -v "+ msgaddress;
  if(!min(isfile(message)))
  {
    open(emailfs, "w");
    writem(emailfs, message);
    close(emailfs);
  } else {
    emailfs = emailfs + " < " + message;
    open(emailfs, "w");
    close(emailfs);
  }
  return 1;
};


fdiff = function( file1, file2, flags )
{
  if (!exist(DIFF))
  { error("fdiff: System command 'diff' not available. Cannot continue!"); }

  if (!exist(file1))
  {
    printf("fdiff: Missing first file. Cannot continue!");
    return 0;
  }
  if (!exist(file2))
  {
    printf("fdiff: Missing second file. Cannot continue!");
    return 0;
  }

  if(type(file1) != "string")
  {
    printf("fdiff: Filename 1 has to be a non-zero length string!");
    return 0;
  }
  if (!strlen(file1))
  {
    printf("fdiff: Filename 1 has to be a non-zero length string!");
    return 0;
  }

  if(type(file2) != "string")
  {
    printf("fdiff: Filename 2 has to be a non-zero length string!");
    return 0;
  }
  if (!strlen(file2))
  {
    printf("fdiff: Filename 2 has to be a non-zero length string!");
    return 0;
  }

  if (!exist(flags))
  {
    // force '--brief' option
    flags = "-q";
  }

  res = reads("| " + DIFF + " "+ flags + " " + file1[1] + " " + file2[1]+ " 2>/dev/null");

  if (isempty(res))
  {
    rval = 0;
  } else {
    rval = 1;
  }

  return rval;
};

touch = function( fns, params )
{
  if (!exist(TOUCH))
  { error("touch: System command not available. Cannot continue !"); }

  if (!exist(fns))
  {
    printf ("touch: Nothing to be done!");
    return 0;
  }
  if(type(fns)!="string")
  {
    printf ("touch: Need string argument for filename!");
    return 0;
  }
  if (!exist(params))
  {
    params="";
  }
  if(type(params)!="string")
  {
    printf ("touch: Need string argument for parameters!");
    return 0;
  }
  for (i in 1:fns.nr)
  {
    for( j in 1:fns.nc)
    {
      system(TOUCH + " " + fns[i;j] + " " + params);
    }
  }
  return 1;
};


edit = function(variable, editor)
{
  if (!exist(VI))
  { error("vi: System command not available. Cannot continue !"); }
  if (!exist(KONSOLE))
  { error("konsole: System command not available. Cannot continue !"); }
  if (!exist(RM))
  { error("rm: System command not available. Cannot continue !"); }


  printf("\nEditing the variable content. When finished, will return to rlab2 prompt!\n\n");
  filename="~/tmp/rlab-edit"+text(int(rand(1,1)*10000))+".tmp";

  // choose editor
  if (!exist(editor))
  {
    editor="vi";
  }
  if (editor=="vi")
  {
    editorpath=VI;
    commandstring=KONSOLE + " -caption 'RLaB2:Edit' -e "+editorpath+" "+filename+" >/dev/null";
  }
  if (editor=="emacs")
  {
    editorpath=EMACS;
    commandstring=editorpath+" "+filename+" 2>/dev/null";
  }
  // write variable to the file in /tmp
  vartype=type(variable);
  if ((vartype=="real")||(vartype=="string"))
  {
    open(filename,"w");
    writem(filename, variable);
    close(filename);
  }
  // edit the file with the editor;
  system(commandstring);
  // when this is done read the same file and set it as the content of variable
  if (vartype=="real")
  {
    open(filename,"r");
    variable=readm(filename);
    close(filename);
  }
  if (vartype=="string")
  {
    variable=reads(filename);
  }
  syscmd=RM + " -rf "+filename+" 2>/dev/null";

  return variable
};


view = function(variable, editor)
{
  if (!exist(VI))
  { error("vi: System command not available. Cannot continue !"); }
  if (!exist(KONSOLE))
  { error("konsole: System command not available. Cannot continue !"); }
  if (!exist(RM))
  { error("rm: System command not available. Cannot continue !"); }

  printf("\nShowing the variable content. Its content cannot be changed!\n\n");
  filename="/tmp/rlab-edit_"+timestamp()+".tmp.dem";
  // choose editor
  if (!exist(editor))
  {
    editor="vi";
  }
  if (editor=="vi")
  {
    editorpath=VI;
    commandstring=KONSOLE + " -caption 'RLaB2:Demo' -e "+editorpath+" "+filename+" >/dev/null &";
  }
  if (editor=="emacs")
  {
    editorpath=EMACS;
    commandstring=editorpath+" "+filename+" 2>/dev/null &";
  }
  // write variable to the file in /tmp and limit the access to -r--------
  vartype=type(variable);
  if (vartype=="int")
  {
    writem(filename, variable, <<format="%.0f">>);
  }
  if (vartype=="real")
  {
    writem(filename, variable, <<format="%g">>);
  }
  if (vartype=="string")
  {
    writem(filename, variable);
  }
  if (vartype=="complex")
  {
    open(filename,"w");
    writem(filename, text(variable));
    close(filename);
  }
  // view the file with the editor;
  system(commandstring);
  return 0;
};

static(_teststat);
_teststat = function(filename, w)
{
  if (!exist(filename))
  { return 0; }

  if(filename.n == 0)
  { return 0; }

  filestatus=zeros(filename);

  if (!exist(w))
  { return filestatus; }

  if(class(filename)=="string")
  {
    for (i in 1:filename.nr)
    {
      for (j in 1:filename.nc)
      {
        x = stat(filename[i;j]);
        if (exist(x.[w]))
        { filestatus[i;j]=1; }
      }
    }
  }

  return filestatus;
};

isfile = function( filename )
{
  return _teststat(filename, "file");
};

isdir = function( filename )
{
  return _teststat(filename, "dir");
};

islink = function( filename )
{
  return _teststat(filename, "link");
};

isblkspec = function( filename )
{
  return _teststat(filename, "block");
};

ischarspec = function( filename )
{
  return _teststat(filename, "char");
};

md5sum = function( s )
{
  if (!exist(s))
  { return blank(0,0); }

  if(type(s)!="string")
  { return blank(0,0); }

  return hash("md5", s + "\n");
};

static(top_entries);
if(!exist(top_entries))
{ top_entries = ["pid", "user", "pr", "ni", "virt", "res", "shr", "s", "%cpu", "%mem", "time", "cmd"]; }
getmemoryusage = function(arg)
{

  if (!exist(arg))
  { arg = "%mem"; }

  n = find (top_entries == tolower(arg));
  if (isempty(n))
  { return blank(0,0); }

  p = text (getpid(), "%.0f");
  u = getenv("USER");
  s = reads("| /usr/bin/top -p " + p + " -b -n 1");
  if (isempty(s))
  {
    printf("getmemoryusage: 'top' failed");
    return blank(0,0);
  }
  s = grep (s, p);
  if (isempty(s))
  {
    printf("getmemoryusage: 'top' failed");
    return blank(0,0);
  }

  if (length(s)!=1)
  {
    printf("getmemoryusage: 'top' failed");
    for (i in range(s))
    {
      printf("\ts[%g] = %s\n", i, s[i]);
    }
    return blank();
  }

  chomp(s);
  if (strlen(s)>0)
  {
    s = strsplt(s,"'BLANK");
    chomp(s);
    if (length(s)>=n)
    { return s[n]; }    
  }
  return blank(0,0);
};

hostname = function()
{
  evar = "HOSTNAME";
  host = getenv(evar);
  if (host == "")
  {
    if (exist(ENV))
    {
      env = reads("|" + ENV);
      _i  = find(strindex(env,evar+"=")==1);
      if (!isempty(_i))
      {
        _i = min(_i);
        host = substr(env[_i],(strlen(evar)+2):1000);
      }
    }
  }

  if (host == "")
  {
    if (isfile("/etc/HOSTNAME"))
    {
      host = readm("/etc/HOSTNAME");
    }
  }

  return host;
};


//
//
//
static(udevadm_propeties);
if (!exist(udevadm_propeties))
{
  udevadm_propeties = [ ...
      "BUSNUM", "DEVNAME", "DEVNUM", "DEVPATH", "DEVTYPE", "DRIVER", "ID_BUS", ...
      "ID_FOR_SEAT", "ID_MODEL", "ID_MODEL_ENC", "ID_MODEL_FROM_DATABASE", ...
      "ID_MODEL_ID", "ID_PATH", "ID_PATH_TAG", "ID_REVISION", "ID_SERIAL", ...
      "ID_SERIAL_SHORT", "ID_USB_INTERFACES", "ID_VENDOR", "ID_VENDOR_ENC", ...
      "ID_VENDOR_FROM_DATABASE", "ID_VENDOR_ID", "MAJOR", "MINOR", ...
      "PRODUCT", "SUBSYSTEM", "TAGS", "TYPE" ...
  ];
}

usbdevname = function(idlist)
{
  if (!exist(FIND))
  { error("usbdevname: System command 'find' not available. Cannot continue !"); }
  if (!exist(UDEVADM))
  { error("usbdevname: System command 'udevadm' not available. Cannot continue !"); }
  if (!exist(GREP))
  { error("usbdevname: System command 'grep' not available. Cannot continue !"); }

  if (class(idlist)!="list")
  {
    printf("usbdevname: Find /dev/... of an usb device based on its attributes\n");
    printf("usbdevname: Format:\n");
    printf("usbdevname:   f = usbdevname(attr)\n");
    printf("usbdevname: 'attr' is a list <<... name='val'; ...>>, where available names are (case insensitive)\n");
    printf("usbdevname: 'f' is a file-link to the device, e.g., /dev/ttyUSB0, or some such\n");
    return blank(0,0);
  }

  for (m in members(idlist))
  {
    i = find(strindex(toupper(m), udevadm_propeties));
    if(isempty(i))
    {
      printf("usbfname: Field '%m' is not recognized!\n", m);
      printf("usbfname: Available fields are (case is not important):\n");
      return blank(0,0);
    }
  }

  PATTERNMATCHING = "";
  m = members(idlist);
  i = find(strindex("DEVNAME", toupper(m)));
  if (!isempty(i))
  {
    for (i1 in i)
    { PATTERNMATCHING = PATTERNMATCHING + " | " + GREP + " " + idlist.[m[i1]]; }
  }

  devs = reads("|"+FIND+" /sys/bus/usb/devices/usb*/ -name dev 2>/dev/null" + PATTERNMATCHING );
  if (isempty(devs))
  { return blank(0,0); }

  devs = rstrip(devs, "/dev");
  rval = blank(0,0);
  for (d in devs)
  {
    if (strlen(d)<1)
    { continue; }

    p = reads("|"+UDEVADM+" info -q property --export -p " + d + " 2>/dev/null");
    chomp(p);
    i = [];
    for (m in members(idlist))
    {
      i1 = find(strindex(p,toupper(m)));
      if (isempty(i1))
      { continue; }
      for (i2 in i1)
      {
        if (strindex(p[i2], idlist.[m]))
        {
          i = [i,1];
          continue;
        }
        else
        {
          i = [];
          break;
        }
      }
      if (isempty(i))
      { break; }
    }
    // provided list of attributes has to match 
    if (length(i)!=size(idlist))
    { continue; }

    j = strindex(p,"DEVNAME=");
    if (isempty(j))
    { continue; }
    k = min(find(j));
    rval = [rval, p[k]];
  }

  if (!isempty(rval))
  {
    rval = lstrip(rval, "'\"DEVNAME='\"");
    rval = rstrip(rval, "'\"'\"");
  }
  else
  {
    return rval;
  }

  // if the user has provided 'devname' or part of it
  // then use it to filter the output with
  m = members(idlist);
  i = find(strindex("DEVNAME", toupper(m)));
  if (!isempty(i))
  {
    j = find(strindex(rval, idlist.[m[i]]));
    if (isempty(j))
    {
      rval = blank(0,0);
    }
    else
    {
      rval = rval[j];
    }
  }

  return rval;
};

INIT = 1;


