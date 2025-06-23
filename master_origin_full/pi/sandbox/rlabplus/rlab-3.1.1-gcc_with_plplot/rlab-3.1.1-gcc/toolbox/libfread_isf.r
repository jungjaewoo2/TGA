//
// working with ngspice from rlab for parameter analysis and such
// marijan kostrun, VI-2012
//
// loosely based on script spice_readfile by Werner Hoch <werner.ho@gmx.de>, 2005
//
static(ISF_KEY_RMVS, ISF_KEY_PARAMS, ISF_KEY_WORDS);
ISF_KEY_RMVS = [":WFMP:"];
ISF_KEY_PARAMS=["ENC", "BN_F", "BYT_O", "COMP", "PT_F"];
ISF_KEY_WORDS=["NR_P", "BYT_N", "BIT_N",  "WFI", ...
    "PT_F", "XUN", "XIN", "XZE", "PT_O", "YUN", "YMU", "YOF", "YZE", ...
        "VSCALE", "HSCALE", "VPOS", "VOFFSET", "HDELAY", "FILTERF"];

freadisf = function(fn)
{
  THIS_SOLVER = "readisf";

  // error handling
  if (!isfile(fn))
  { error("File '"+fn+"' does not exist!"); }

  if (open(fn,"r"))
  { error(THIS_SOLVER + ": The file is already opened. Close it first!"); }

  // this is were we save data to
  s = <<>>;

  // read the file header of the file
  // these are colonseperated key/value pairs
  l = char(fread(fn,,"char",,":CURV"));

  // these words are removed from file as they do not carry information
  for (rmv in ISF_KEY_RMVS)
  {
    if (strindex(l, rmv))
    {
      l = gsub("", rmv, l).string;
    else
      printf(THIS_SOLVER + ": keyword %s not found. Is there a problem with %s?\n", rmv, fn);
    }
  }


  // extract commands
  for (param in ISF_KEY_PARAMS)
  {
    i1 = strindex(l, param);
    if (!i1)
    { continue; }
    di = strindex(substr(l,[i1:(i1+10000)]),";")-1;
    p = substr(l,[i1:(i1+di)]);
    l = gsub("", p, l).string;
    s.[tolower(param)] = substr(p, (strlen(param)+2):(strlen(p)-1));
  }

  // process 'l' using ISF KEY WORDS:
  // numbers:
  //  ISF_KEY_WORD<space>value;
  // text
  //  ISF_KEY_WORD<space>"value";
  for (keyword in ISF_KEY_WORDS)
  {
    i1 = strindex(l, keyword);
    if (!i1)
    { continue; }
    di = strindex(substr(l,[i1:(i1+10000)]),";")-1;
    p = substr(l,[i1:(i1+di)]);
    l = gsub("", p, l).string;
    s.[tolower(keyword)] = eval(substr(p, (strlen(keyword)+2):(strlen(p)-1)));
  }

  // read in data:
  // :CURV{E}#
  if (s.enc == "BIN")
  {
    l = char(fread(fn,,"unsigned char",,"#"));     // gobble up to '#'
    n = strtod(char(fread(fn,1,"unsigned char"))); // read how long is data-length string
    c = strtod(char(fread(fn,n,"unsigned char"))); // read data-length string
    if (c != s.nr_p)
    {
      printf(THIS_SOLVER + ": nr_p = %g, which is different from %g\n", s.nr_p, c);
    }
    if (exist(s.xin) && exist(s.ymu) && exist(s.yof))
    {
      s.data = zeros(s.nr_p,2);
      s.data[;1] = [1:s.nr_p]' * s.xin;
      if (s.byt_n == 1)
      {
        s.data[;2] = s.ymu * (fread(fn,c,"char")' - s.yof);
      }
    else
      printf(THIS_SOLVER + ": Cannot retrieve data from the file: Missing 'xin', 'ymu' or 'yof' fields!\n");
    }
  }

  close(fn);

  return s;
};


