//
// working with ngspice from rlab for parameter analysis and such
// marijan kostrun, VI-2012
//
// loosely based on script spice_readfile by Werner Hoch <werner.ho@gmx.de>, 2005


static (NGSPICE_CMD, RM_CMD, NGSPICE_RAWFILE_DEF, NGSPICE_TEMP_CKTFILE);
static (NGSPICE_LOGFILE_DEF, NGSPICE_ERRFILE_DEF, NGSPICE_FLAGS);

// ngspice executable
// NGSPICE_CMD_PIPE = "|ngspice -n -i 1>/dev/null 2>/dev/null";
// NGSPICE_CMD_FORK = "ngspice -n -i 1>/dev/null 2>/dev/null";
NGSPICE_RAWFILE_DEF = "./rawspice.raw";
NGSPICE_TEMP_CKTFILE = "./.tmp.spc";

// ngspice executable
NGSPICE_CMD = "ngspice";
NGSPICE_FLAGS = " -i -n -b ";

// os commands
RM_CMD = "/bin/rm -rf 2>/dev/null ";


//
// this is where our commands go
//
if (!exist(spice))
{ spice = <<>>; }

//
// execute spice script and return the result
//
spice.exec = function (script, fn_scp, fn_log, fn_err)
{
  eraseme = 0;
  my_ngspice_cmd = NGSPICE_CMD;

  // did user provide a file name to which she wants the ngspice script
  // to be written?
  if (exist(fn_scp))
  {
    f1 = fn_scp;
  else
    f1 = NGSPICE_TEMP_CKTFILE;
    eraseme = 1;
  }

  if (exist(script))
  {
    if (class(script)!="string")
    { return 1; }
    open(f1,"w");
    for (i in 1:script.n)
    {
      fprintf(f1, "%s\n", script[i]);
    }
    close(f1);
  }
  if(!isfile(f1))
  { return 2; }

  // did user provide log file?
  if (exist(fn_log))
  {
    if (class(fn_log)=="string")
    {
      my_ngspice_cmd = my_ngspice_cmd + " 1>" + fn_log[1];
    else
      return 3;
    }
  else
    my_ngspice_cmd = my_ngspice_cmd + " 1>/dev/null";
  }

  // did user provide error file?
  if (exist(fn_err))
  {
    if (class(fn_err)=="string")
    {
      my_ngspice_cmd = my_ngspice_cmd + " 2>" + fn_err[1];
    else
      return 3;
    }
  else
    my_ngspice_cmd = my_ngspice_cmd + " 2>/dev/null";
  }

  // call ngspice and execute the script, but do not move to the background
  system(my_ngspice_cmd + NGSPICE_FLAGS + f1);

  // clean-up
  if (eraseme)
  { system(RM_CMD + f1); }

  return 0;
};

spice.fork = function (script, fn_scp, fn_raw, fn_log)
{
  eraseme = 0;
  my_ngspice_cmd = NGSPICE_CMD;

  // did user provide a file name to which she wants the ngspice script
  // to be written?
  if (exist(fn_scp))
  {
    f1 = fn_scp;
  else
    f1 = NGSPICE_TEMP_CKTFILE;
    eraseme = 1;
  }

  if (exist(script))
  {
    if (class(script)!="string")
    { return 1; }
    open(f1,"w");
    for (i in 1:script.n)
    {
      fprintf(f1, "%s\n", script[i]);
    }
    close(f1);
  }
  if(!isfile(f1))
  { return 2; }


  if (exist(fn_raw))
  {
    if (class(fn_raw)=="string")
    {
      my_ngspice_cmd = my_ngspice_cmd + " -r " + fn_raw[1];
      NGSPICE_RAWFILE = fn_raw[1];
    else
      return 2;
    }
  }

  if (exist(fn_log))
  {
    if (class(fn_log)=="string")
    {
      my_ngspice_cmd = my_ngspice_cmd + " -o " + fn_log[1];
    else
      return 3;
    }
  }

  // call ngspice and execute the script in the background
  system("( setsid " + my_ngspice_cmd + " < " + f1 +" & disown)");

  // clean-up
  if (eraseme)
  { system(RM_CMD + f1); }

  return 0;
};


spice.raw = function(fn)
{
  // error handling
  if (!exist(fn))
  {
    if (exist(NGSPICE_RAWFILE))
    {
      fn = NGSPICE_RAWFILE;
    else
      fn = NGSPICE_RAWFILE_DEF;
    }
  }

  if (!isfile(fn))
  { error("File '"+fn+"' does not exist!"); }

  if (open(fn,"r"))
  {
    error("spice.raw: The file is already opened. Close it first!");
  }

  // defaultvalues
  realflag   = 1;   // number type is real, not complex
  paddedflag = 1;   // with zeros padded data
  binaryflag = 1;   // binary data type

  // this is were we save data to
  s = <<>>;

  // read the file header of the file
  // these are colonseperated key/value pairs
  while (1)
  {
    line = fgets(fn);

    // read the fixed string fields
    for (ff in ["title", "date", "plotname", "command", "option"])
    {
      if (strindex(tolower(line), ff + ": "))
      {
        s.[ff] = lstrip(substr(line,strlen(ff + ": "):strlen(line)), " ");
        s.[ff] = rstrip(s.[ff], " ");
      }
    }

    // flags
    if (strindex(tolower(line), "flags"))
    {
      sf = lstrip(substr(line,strlen("flags: "):strlen(line)), " ");
      sf = rstrip(sf, " ");
      if (sf == "")
      { printf("spice.raw: No flags specified!"); }
      sf = strsplt(sf," ");
      n = length(sf);
      if (n>2)
      { printf("spice.raw: %g flags given while only one or two expected!", n); }

      for (i in 1:n)
      {
        if (strindex(sf[i],"real"))
        {
          realflag=1;
        else if (strindex(sf[i],"complex"))
        {
          realflag=0;
        else if (strindex(sf[i],"padded"))
        {
          paddedflag=1;
        else if (strindex(sf[i],"unpadded"))
        {
          paddedflag=0;
          printf("spice.raw: unpadded data not handled yet");
        else
          printf("spice.raw: unknown flag %s", sf[i]);
        }}}}
      }

      s.flags = sf;
    }

    // number parameters
    for (ff in ["variables", "points"])
    {
      if (strindex(tolower(line), "no. " + ff + ":"))
      {
        s.[ff] = lstrip(substr(line,strlen("no. " + ff + ": "):strlen(line)), " ");
        s.[ff] = rstrip(s.[ff], " ");
        s.[ff] = strtod( s.[ff] );
      }
    }

    // did we get to the list of the variables?
    if (strindex(tolower(line), "variables:") == 1)
    {
      if (!exist(s.variables) || !exist(s.points))
      { error("spice.raw: missing declaration of variables and their length"); }

      names = blank(1,s.variables);
      units = blank(1,s.variables);

      // go over the rows and figure out the names of the variables and their units
      // field separator can be '\t' or 'space'
      // format is
      //    {i-1} 'name' 'unit'
      for (i in 1:s.variables)
      {
        line = fgets(fn);
        line = gsub(" ", "\t", line).string;
        line = lstrip(line, " ");
        line = rstrip(line, " ");

        dummy = strsplt(line," ");
        if (length(dummy)!=3)
        { error("spice.raw: list of variables improperly formatted"); }

        names[i] = dummy[2];
        units[i] = dummy[3];
      }

      s.names = names;
      s.units = units;
    }

    // if we get to 'values' or 'binary' we can start reading the data points
    if (strindex(tolower(line), "values:") == 1)
    {
      binaryflag = 0;
      break;
    else if (strindex(tolower(line), "binary:") == 1)
    {
      binaryflag = 1;
      break
    }}
  } // while(1)


  if (!exist(s.variables) || !exist(s.points))
  { error("spice.raw: missing declaration of variables and their length"); }

  // read the data
  if (binaryflag==1 && realflag==1)
  {
    // binary and real
    data = fread(fn, s.variables * s.points, "double");
    resize(data, s.variables, s.points);
    s.data = data';
  else if (binaryflag==0 && realflag==1)
  {
    // ascii and real
    s.data = zeros(s.points,s.variables);
    for (i in 1:(s.points+1))
    {
      sd = "";
      for (j in 1:(s.variables+1))
      { sd = sd + rstrip(fgets(fn),"\n"); }
      s.data[i;] = strtod(sd)[2:(s.variables+1)];
    }
  else if (binaryflag==1 && realflag==0)
  {
    // binary and complex
    data = fread(fn, 2 * s.variables * s.points, "double");
    resize(data, 2 * s.variables, s.points);
    data = data[1:2 * s.variables:2;] + 1i * data[2:2 * s.variables:2;];
    s.data = data';
  else
    // ascii and complex
    data = zeros(s.points,2*s.variables);
    for (i in 1:(s.points+1))
    {
      sd = "";
      for (j in 1:(2*s.variables+1))
      { sd = sd + rstrip(fgets(fn),"\n"); }
      data[i;] = strtod(sd)[2:(2*s.variables+1)];
    }
    s.data = data[; 1:2 * s.variables:2] + 1i * data[; 2:2 * s.variables:2];
  }}}

  close(fn);

  //
  // at this point we have a list with entries
  //
  rval = <<>>;
  rval.variables = s.variables;
  rval.points    = s.points;

  rval.data = <<>>;
  dj = 1;
  if (!strindex(s.flags,"real"))
  { dj = 2; }

  if (s.data.nc == dj * length(s.names))
  {
    for (i in range(s.names))
    {
      si = rstrip(gsub("_", "(", s.names[i]).string,")");
      si = gsub("_", ".", si).string;
      rval.data.[ si ] = s.data[;(dj * (i - 1) + 1) : (dj * i)];
    }
  else
    printf("error occured during processing of the data read from the 'raw' file\n");
    printf("expected %g columns, but found only %g\n", dj * length(s.names), s.data.nc);
    for (i in range(s.names))
    {
      si = rstrip(gsub("_", "(", s.names[i]).string,")");
      si = gsub("_", ".", si).string;
      rval.data.[ si ] = [];
    }
    rval.raw_data = s.data;
    printf("retrieved data can be found in entry 'raw_data' in the result list\n");
  }
  return rval;
};


