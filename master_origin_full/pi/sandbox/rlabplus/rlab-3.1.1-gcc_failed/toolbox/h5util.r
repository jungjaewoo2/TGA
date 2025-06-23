//
// library for work with HDF5 files
//  uses h5tools installed with the library for retrieving
//  information about the structure of a h5 file

static(H5LS_OPTIONS, H5DUMP_OPTIONS,H5CHECK_OPTIONS,H5REPACK_OPTIONS);
H5LS_OPTIONS = "--width=1000 -r";
H5DUMP_OPTIONS = "-H";
H5CHECK_OPTIONS = "-v2";
H5REPACK_OPTIONS = "-v";

static(WHICH,H5DUMP, H5LS, H5CHECK, H5REPACK);
WHICH=reads("|which which");
if (strlen(WHICH)>0)
{
  // check for 'h5ls'
  H5LS = reads("|" + WHICH + " h5ls 2>/dev/null");
  if (isempty(H5LS))
  { clear(H5LS); }
  if (strlen(H5LS)<strlen("H5LS"))
  { clear(H5LS);}

  // check for 'h5dump'
  H5DUMP = reads("|" + WHICH + " h5dump 2>/dev/null");
  if (isempty(H5DUMP))
  { clear(H5DUMP); }
  if (strlen(H5DUMP)<strlen("H5DUMP"))
  { clear(H5DUMP);}

  // check for 'h5check'
  H5CHECK = reads("|" + WHICH + " h5check 2>/dev/null");
  if (isempty(H5CHECK))
  { clear(H5CHECK); }
  if (strlen(H5CHECK)<strlen("H5CHECK"))
  { clear(H5CHECK);}

  // repack for 'h5repack'
  H5REPACK = reads("|" + WHICH + " h5repack 2>/dev/null");
  if (isempty(H5REPACK))
  { clear(H5REPACK); }
  if (strlen(H5REPACK)<strlen("H5REPACK"))
  { clear(H5REPACK);}
}

//
// h5list: list the content of an h5 file without opening it within rlab
//
h5list = function(file, object, options)
{
  if (!exist(file))
  { error("h5ls: File does not exist!"); }

  if (!isempty(findstr(file, "h5://")))
  {
    name = substr(file, (strlen("h5://")+1):strlen(file));
    if (strlen(name)==0)
    { stop("h5ls: File name cannot be of zero length!"); }
  else if (!isempty(findstr(file, "hdf5://")))
  {
    name = substr(file, (strlen("hdf5://")+1):strlen(file));
    if (strlen(name)==0)
    { stop("h5ls: File name cannot be of zero length!"); }
  else
    name = file;
  }}

  if (!exist(options))
  { options = H5LS_OPTIONS; }

  if (!exist(object))
  {
    object = "/";

    cmd = "|" + H5LS + " " + options + " " + name + object + " 2>/dev/null";
    x = reads(cmd);

    if (all(strlen(x)==0))
    {
      stop("h5ls: File '" + file + "' does not exist or is not proper HDF5 file!");
    }

    j = 0;
    for(i in 1:x.nr)
    {
      if (substr(x[i],1)== "/")
      { j++; }
    }
    rval = blank(j,3);
    j = 0;
    for(i in 1:x.nr)
    {
      if (substr(x[i],1)== "/")
      {
        j++;
        k = findstr (x[i], " ");
        n = strlen  (x[i]);
        rval[j;1] = substr(x[i], 1:(min(k)-1));
        // links are most easily identified in the output of 'h5ls -r'
        s2 = substr(x[i], (max(k)+1):n );
        if (substr(s2,1)=="/")
        {
          rval[j;2] = "HARDLINK";
          rval[j;3] = s2;
          continue;
        else
          rval[j;3] = "";
        }
        // type of object follows its name
        k2 = min(k);
        while (substr(x[i],k2)==" ")
        { k2++; }
        k3 = k2+1;
        while (substr(x[i],k3)!=" " && k3<=n)
        { k3++; }
        rval[j;2] = toupper(substr(x[i], k2:k3));
        if (rval[j;2]=="TYPE")
        { rval[j;2] = "DATATYPE"; }
      }
    }
    return rval;
  }

  // find the list of all objects in the file:
  x = $self(file);

  // make path to the object absolute by prepending '/' to it
  if(substr(object,1)!= "/")
  { object = "/" + object; }

  // on how many places does the object appear on the first place
  y  = strindex(x[;1], object);
  ny = sum(y==1);
  if (ny == 0)
  { return blank(0,0); }

  rval = blank(ny, 3);
  j = 0;
  for (i in 1:x.nr)
  {
    if (y[i] == 1)
    {
      j++;
      rval[j;] = x[i;];
    }
  }

  return rval;
};

h5dump = function(file, options)
{
  if (!exist(file))
  { error("h5dump: File does not exist!"); }

  if (!exist(options))
  { options = H5DUMP_OPTIONS; }

  cmd = H5DUMP + " " + options + " " + file;
  system(cmd);

  return 0;
};

h5check = function(file)
{
  if (!exist(file))
  { error("h5check: File does not exist!"); }

  if (!exist(options))
  { options = H5CHECK_OPTIONS; }

  cmd = "|" + H5CHECK + " " + options + " " + file + " | grep \"No non-compliance errors found\"";
  x = reads(cmd);

  if (strindex(x, "No non-compliance errors found"))
  { return 0; }

  return 1;
};

h5repack = function(file1, file2, options)
{
  if (!exist(file1))
  { error("h5repack: Source file does not exist!"); }

  if (!exist(file2))
  { error("h5repack: Need new file name!"); }

  if (file1==file2)
  { error("h5repack: Need two different names!"); }

  if (!exist(options))
  { options = H5REPACK_OPTIONS; }

  cmd = "|" + H5REPACK + " " + options + " " + file1 + " "+ file2 + " 2>/dev/null";

  system(cmd);

  return 0;
};

