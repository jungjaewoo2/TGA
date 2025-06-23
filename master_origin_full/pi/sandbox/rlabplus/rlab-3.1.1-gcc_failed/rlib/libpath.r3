//
//
//
static(INIT);

if (exist(INIT))
{ EOF }

path = function ()
{
  return getenv("RLAB2_PATH");
};

path_add = function(s)
{
  if (class(s)=="string")
  {
    if (length(s)>1)
    { s = join(s,":"); }
    current_path = getenv("RLAB2_PATH");
    if (strlen(current_path)>0)
    {
      putenv("RLAB2_PATH="+current_path+":"+s);
    }
    else
    {
      putenv("RLAB2_PATH="+s);
    }
  }

  return getenv("RLAB2_PATH");
};


path_extname = function(fns)
{
  if (isempty(fns))
  { return blank(); }

  ext = blank(fns);
  for (i in range(fns))
  {
    fn = fns[i];
    idx_dot = findstr(fn, ".");
    if (isempty(idx_dot))
    {
      ext[i] = "";
      continue;
    }
    idx_dot = last(idx_dot);
    idx_slash = findstr(fn, "/");
    if (isempty(idx_slash))
    { idx_slash = 0;}
    idx_slash = last(idx_slash);
    if (idx_slash > idx_dot)
    {
      ext[i] = "";
      continue;
    }
    ext[i] = substr(fn,(idx_dot+1):strlen(fn));
  }
  return ext;
};

path_sans_extname = function(fns)
{
  if (isempty(fns))
  { return blank(); }

  p = blank(fns);
  for (i in range(fns))
  {
    fn = fns[i];
    idx_dot = findstr(fn, ".");
    if (isempty(idx_dot))
    {
      p[i] = fn;
      continue;
    }
    idx_dot = last(idx_dot);
    idx_slash = findstr(fn, "/");
    if (isempty(idx_slash))
    { idx_slash = 0;}
    idx_slash = last(idx_slash);
    if (idx_slash > idx_dot)
    {
      p[i] = fn;
      continue;
    }
    p[i] = substr(fn,1:(idx_dot-1));
  }
  return p;
};


//
// system related path functions
//
static (WHICH);
if (!exist(WHICH))
{ WHICH=reads("|which which"); }
if (strlen(WHICH)<1)
{
  INIT = 1;
  EOF
}

// check for 'pwd'
static(PWD);
PWD = reads("|" + WHICH + " pwd 2>/dev/null");
if (isempty(PWD))
{ clear(PWD); }
if (strlen(PWD)<strlen("pwd"))
{ clear(PWD);}

pwd = function()
{
  rval = getenv("PWD");

  if (strlen(rval)<1)
  {
    if (exist(PWD))
    {
      env = reads("|" + PWD);
      _i  = find(strindex(env,evar+"=")==1);
      if (!isempty(_i))
      {
        rval = substr(env[_i],(strlen(evar)+2):1000);
      }
    }
  }

  return rval;
};

// check for 'ls'
static(LS);
LS = reads("|" + WHICH + " ls 2>/dev/null");
if (isempty(LS))
{ clear(LS); }
if (strlen(LS)<strlen("LS"))
{ clear(LS);}

// check for 'mv'
static(MV);
MV = reads("|" + WHICH + " mv 2>/dev/null");
if (isempty(MV))
{ clear(MV); }
if (strlen(MV)<strlen("MV"))
{ clear(MV);}

// check for 'cp'
static(CP);
CP = reads("|" + WHICH + " cp 2>/dev/null");
if (isempty(CP))
{ clear(CP); }
if (strlen(CP)<strlen("CP"))
{ clear(CP);}

// check for 'mkdir'
static(MKDIR);
MKDIR = reads("|" + WHICH + " mkdir 2>/dev/null");
if (isempty(MKDIR))
{ clear(MKDIR); }
if (strlen(MKDIR)<strlen("MKDIR"))
{ clear(MKDIR);}

// check for 'rm'
static (RM);
RM = reads("|" + WHICH + " rm 2>/dev/null");
if (isempty(RM))
{ clear(RM); }
if (strlen(RM)<strlen("RM"))
{ clear(RM);}

ls = function( dir, flags, grep )
{
  if (!exist(LS))
  {
    printf("ls: System command not available. Cannot continue !\n");
    return blank(0,0);
  }

  if (!exist(dir))
  { dir="./"; }

  if (!exist(flags))
  { flags = "-1A"; }

  if(type(dir)!="string")
  { return blank(0,0); }

  // process grep command
  if (!exist(grep))
  {
    grep = "";
  }
  if (class(grep)!="string")
  {
    grep = "";
  }
  if (strlen(grep) > 0)
  {  grep = " | " + GREP + " " + grep; }

  x = reads("|" + LS + " " + flags + " " + dir + " 2>/dev/null " + grep);
  if (isempty(x))
  { return blank(0,0); }

  //
  // chomp-off empty strings at the end
  //
  il = x.n;
  while (strlen(x[il]) == 0)
  {
    il--;
    if (il==0) { break; }
  }
  if (il > 0)
  {
    return x[1:il;];
  } else {
    return blank(0,0);
  }
};

mv = function( file1, file2, flags )
{
  if (!exist(MV))
  { error("mv: System command not available. Cannot continue !"); }

  if (!exist(file1))
  {
    printf("mv: missing file!");
    return 0;
  }
  if (!exist(file2))
  {
    printf("mv: missing destination !");
    return 0;
  }
  if(type(file1) != "string")
  {
    printf("mv: filename has to be a string!");
    return 0;
  }
  if (!exist(flags))
  {
    // force non-interactive mode
    flags = "-f";
  }
  else
  {
    // cannot move files in interactive mode
    flags = gsub("-f","-i",flags).string
        flags = gsub("f","i",flags).string;
  }
  return system(MV + " "+ flags + " " + file1 + " " + file2+ " 2>&1");
};

rm = function( file1, flags )
{
  if (!exist(RM))
  { error("rm: System command not available. Cannot continue !"); }

  if (!exist(file1))
  {
    printf("mv: missing file!");
    return 0;
  }
  if(type(file1) != "string")
  {
    printf("mv: filename has to be a string!");
    return 0;
  }
  if (!exist(flags))
  {
    // force non-interactive mode
    flags = "-rf";
  } else {
    // cannot move files in interactive mode
    flags = gsub("-f","-i",flags).string
        flags = gsub("f","i",flags).string;
  }
  return system(RM + " "+ flags + " " + file1 + " 2>&1");
};

cp = function( file1, file2, flags )
{
  if (!exist(CP))
  { stop("cp: System command not available. Cannot continue !"); }

  if (!exist(file1))
  {
    printf("cp: missing file!");
    return 0;
  }
  if(type(file1) != "string")
  {
    printf("cp: filename has to be a string!");
    return 0;
  }
  // force archive copy mode
  if (!exist(flags))
  { flags = "-a"; }

  return system(CP + " "+ flags + " " + file1 + " " + file2+ " 2>&1");
};


mkdir = function( dir )
{
  if (!exist(MKDIR))
  { error("mkdir: System command not available. Cannot continue !"); }

  local(i,j,cdir);
  if (!exist(dir))
  {
    return 0;
  }
  if (type(dir)!="string")
  {
    return 0;
  }
  rval = [];
  for(i in 1:dir.nr)
  {
    for (j in 1:dir.nc)
    {
      cdir=dir[i;j];
      if (length(cdir)>0)
      {
        pipecmd = MKDIR + " -p " + cdir + " 2>&1";
        sval = system(pipecmd);
        if (length(sval) > 0)
        {
          rval = [rval; sval];
        }
      }
    }
  }
  if(rval.nr * rval.nc != 0)
  { return rval; }

  return 1;
};




INIT = 1;


