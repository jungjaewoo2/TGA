//
// libgnuplot.so.r
// loader for enhanced gnuplot access library
//

static(_gnukey, _gnulegend, GNU_DEFAULTS, _gnuwins);
GNU_DEFAULTS = <<>>;
GNU_DEFAULTS.debug    = 0;
GNU_DEFAULTS.resetcmd = "unset multiplot; clear;";    // gnuwins
GNU_DEFAULTS.reseterm = "set term x11; set output";   // after printing to a file, this is how
                                                      // the output is redirected to the window
GNU_DEFAULTS.epscolor = "postscript eps enh color blacktext solid lw 2 \"Times-Roman\" 22";
GNU_DEFAULTS.epsmono  = "postscript eps enh mono blacktext solid lw 2 \"Times-Roman\" 22";

//
// load functions from libgnuplot.so
//
_HOME_ = getenv("HOME");
_LIBD_ = "/rlab/lib.so/gnuplot/rlabplus_libgnuplot.so";
fileaddr = _HOME_ + _LIBD_ ;


//
// instantiate gnuplot device (file or pipe)
//  I = gnustart(filename)
//  I = gnustart(), I = gnustart(,/term/,/stderr/)
//
if (!exist(gnustart))
{ gnustart = dlopen(fileaddr, "ent_gnuplot_init"); }

//
// close existing gnuplot device
//  gnuclose(I)
//
if (!exist(gnuclose))
{ gnuclose = dlopen(fileaddr, "ent_gnuplot_close"); }

//
// obtain or set default gnuplot device
// I = gnuwin()
// gnuwin(I)
//
if (!exist(gnuwin))
{ gnuwin = dlopen(fileaddr, "ent_gnuplot_win"); }

//
// obtain information about gnuplot devices
//  <<available;default;device>> = gnuwins()
//
if (!exist(_gnuwins))
{ _gnuwins = dlopen(fileaddr, "ent_gnuplot_allwin"); }

//
// send command to default gnuplot device
//  gnucmd(["cmd1" ...])
//  
if (!exist(gnucmd))
{
  gnucmd = dlopen(fileaddr, "ent_gnuplot_cmd");
}

//
// print data arrays to default gnuplot device
//  gnuprint(x)
//  gnuprint(x, y, z)
//  
if (!exist(gnuprint))
{
  gnuprint = dlopen(fileaddr, "ent_gnuplot_print");
}
clear(fileaddr,_HOME_,_LIBD_);

//
//
// scripted functions for manipulations with gnuplot
//
//

        
gnudebug = function (i)
{
  if (exist(i))
  {
    if (i)
    {
      GNU_DEFAULTS.debug = 1;
    else
      GNU_DEFAULTS.debug = 0;
    }
  }
  
  return GNU_DEFAULTS.debug;
};

gnutitle = function ( s )
{
  I = gnuwin();
  if (isempty(I))
  { gnuwins(1); }

  cmd = "set title;";
  if (exist(s))
  {
    if (strlen(s[1]) > 0)
    { cmd = "set title '" + s[1] + "';"; }
  }

  if (GNU_DEFAULTS.debug)
  { printf("gnutitle: %s\n", cmd); }
  
  gnucmd( cmd );
  return 1;
};


gnuscale2 = function ( sx, sy )
{
  I = gnuwin();
  if (isempty(I))
  { gnuwins(1); }

  cmd = "unset logscale x2;";
  if (exist(sx))
  {
    if (strindex(tolower(sx),"log"))
    { cmd = "set logscale x2;"; }
  }
  if (GNU_DEFAULTS.debug)
  { printf("gnuscale: %s\n", cmd); }
  gnucmd(cmd);

  cmd = "unset logscale y2;";
  if (exist(sy))
  {
    if (strindex(tolower(sy),"log"))
    { cmd = "set logscale y2;"; }
  }
  if (GNU_DEFAULTS.debug)
  { printf("gnuscale: %s\n", cmd); }
  gnucmd(cmd);

  return 1;
};


gnuscale = function ( sx, sy, sz )
{
  I = gnuwin();
  if (isempty(I))
  { gnuwins(1); }

  cmd = "unset logscale x;";
  if (exist(sx))
  {
    if (strindex(tolower(sx),"log"))
    { cmd = "set logscale x;"; }
  }
  if (GNU_DEFAULTS.debug)
  { printf("gnuscale: %s\n", cmd); }
  gnucmd(cmd);

  cmd = "unset logscale y;";
  if (exist(sy))
  {
    if (strindex(tolower(sy),"log"))
    { cmd = "set logscale y;"; }
  }
  if (GNU_DEFAULTS.debug)
  { printf("gnuscale: %s\n", cmd); }
  gnucmd(cmd);

  if (exist(sz))
  {
    cmd = "unset logscale z;";
    if (exist(sz))
    {
      if (strindex(tolower(sz),"log"))
      { cmd = "set logscale z"; }
    }
    if (GNU_DEFAULTS.debug)
    { printf("gnuscale: %s\n", cmd); }
    gnucmd(cmd);
  }

  return 1;
};


gnuztics = function ( x1, x2, fmt )
{
  I = gnuwin();
  if (isempty(I))
  { gnuwins(1); }
  
  if (!exist(fmt))
  { fmt = "%%g"; }
  if (class(fmt)!="string")
  { fmt = "%%g"; }

  cmd = "set format z '$1'";

  _gnuztics = "";

  if (exist(x1))
  {
    _gnuztics = _gnuztics + "set format z '" + fmt[1] + "'; ";
    if (class(x1) == "num")
    {
      // x1 is a vector
      if (length(x1) > 1)
      {
        _gnuztics = _gnuztics + "set ztics (";
        for (i in 1:length(x1)-1)
        { _gnuztics = _gnuztics + text(x1[i]) + ","; }
        _gnuztics = _gnuztics + text(x1[i+1]) + ");";
      else
      // x1 is a single value
        _gnuztics = _gnuztics + "set ztics " + text(x1) + ";";
      }
    }
    if (class(x1) == "string")
    { _gnuztics = _gnuztics + x1 + ";"; }
  else
    _gnuztics = _gnuztics + "set ztics;";
  }

  if (exist(x2))
  {
    if (class(x2) == "num")
    {
      // x2 is a single value
      _gnuztics = _gnuztics + "set mztics " + text(x2) + ";";
    }
    if (class(x2) == "string")
    { _gnuztics = _gnuztics + x2 +";"; }
  else
    _gnuztics = _gnuztics + "unset mztics;";
  }

  gnucmd(_gnuztics);
  return 1;
};


gnuxtics = function ( x1, x2, fmt )
{
  I = gnuwin();
  if (isempty(I))
  { error("gnuplot: non-device accessed"); }

  if (!exist(fmt))
  { fmt = "%%g"; }
  if (class(fmt)!="string")
  { fmt = "%%g"; }

  _gnuxtics = "unset x2tics; ";

  if (exist(x1))
  {
    _gnuxtics = _gnuxtics + "set format x '" + fmt[1] + "'; ";
    if (class(x1) == "num")
    {
      // x1 is a vector
      if (length(x1) > 1)
      {
        _gnuxtics = _gnuxtics + "set xtics (";
        for (i in 1:length(x1)-1)
        { _gnuxtics = _gnuxtics + text(x1[i]) + ","; }
        _gnuxtics = _gnuxtics + text(x1[i+1]) + ");";
        else
      // x1 is a single value
          _gnuxtics = _gnuxtics + "set xtics " + text(x1) + ";";
      }
    }
    if (class(x1) == "string")
    { _gnuxtics = _gnuxtics + x1 + ";"; }
  else
    _gnuxtics = _gnuxtics + "set xtics;";
  }

  if (exist(x2))
  {
    if (class(x2) == "num")
    {
      // x2 is a single value
      _gnuxtics = _gnuxtics + "set mxtics " + text(x2) + ";";
    }
    if (class(x2) == "string")
    { _gnuxtics = _gnuxtics + x2 +";"; }
  else
    _gnuxtics = _gnuxtics + "unset mxtics;";
  }

  gnucmd(_gnuxtics);
  return 1;
};


gnux2tics = function ( x1, x2, fmt )
{
  I = gnuwin();
  if (isempty(I))
  { error("gnuplot: non-device accessed"); }

  if (!exist(fmt))
  { fmt = "%%g"; }
  if (class(fmt)!="string")
  { fmt = "%%g"; }

  _gnuxtics = "";

  if (exist(x1))
  {
    _gnuxtics = _gnuxtics + "set format x2 '" + fmt[1] + "'; ";
    if (class(x1) == "num")
    {
      // x1 is a vector
      if (x1.nr*x1.nc > 1)
      {
        _gnuxtics = _gnuxtics + "set x2tics (";
        for (i in 1:x1.nr*x1.nc-1)
        { _gnuxtics = _gnuxtics + text(x1[i]) + ","; }
        _gnuxtics = _gnuxtics + text(x1[i+1]) + ");";
        else
      // x1 is a single value
          _gnuxtics = _gnuxtics + "set x2tics " + text(x1) + ";";
      }
    }
    if (class(x1) == "string")
    { _gnuxtics = _gnuxtics + x1 + ";"; }
  else
    _gnuxtics = _gnuxtics + "unset x2tics; ";
  }

  if (exist(x2))
  {
    if (class(x2) == "num")
    {
      // x2 is a single value
      _gnuxtics = _gnuxtics + "set mx2tics " + text(x2) + ";";
    }
    if (class(x2) == "string")
    { _gnuxtics = _gnuxtics ";" + x2; }
  else
    _gnuxtics = _gnuxtics + ";" + "unset mx2tics;";
  }

  gnucmd(_gnuxtics);
  return 1;
};


gnuytics = function ( x1, x2, fmt )
{
  I = gnuwin();
  if (isempty(I))
  { error("gnuplot: non-device accessed"); }
  
  if (!exist(fmt))
  { fmt = "%%g"; }
  if (class(fmt)!="string")
  { fmt = "%%g"; }

  _gnuytics = "unset y2tics; ";

  if (exist(x1))
  {
    _gnuytics = _gnuytics + "set format y '" + fmt[1] + "'; ";
    if (class(x1) == "num")
    {
      // x1 is a vector
      if (x1.nr*x1.nc > 1)
      {
        _gnuytics = _gnuytics + "set ytics (";
        for (i in 1:x1.nr*x1.nc-1)
        { _gnuytics = _gnuytics + text(x1[i]) + ","; }
        _gnuytics = _gnuytics + text(x1[i+1]) + ");";
      else
      // x1 is a single value
        _gnuytics = _gnuytics + "set ytics " + text(x1) + ";";
      }
    }
    if (class(x1) == "string")
    { _gnuytics = _gnuytics + x1 + ";"; }
  else
    _gnuytics = _gnuytics + "set ytics; ";
  }

  if (exist(x2))
  {
    if (class(x2) == "num")
    {
      // x2 is a single value
        _gnuytics = _gnuytics + "set mytics " + text(x2) + ";";
    }
    if (class(x2) == "string")
    { _gnuytics = _gnuytics + x2; }
  else
    _gnuytics = _gnuytics + "unset mytics;";
  }

  gnucmd(_gnuytics);
  return 1;
};


//
// gnuy2tics ([y1,..yK], DY)
// gnuy2tics ("set y2tics (*)", "set my2tics *")
//
gnuy2tics = function ( x1, x2, fmt )
{
  I = gnuwin();
  if (isempty(I))
  { error("gnuplot: non-device accessed"); }

  if (!exist(fmt))
  { fmt = "%%g"; }
  if (class(fmt)!="string")
  { fmt = "%%g"; }

  if (exist(x1))
  {
    _gnuytics = "set format y2 '" + fmt[1] + "'; ";
    if (class(x1) == "num")
    {
      // x1 is a vector
      if (x1.nr*x1.nc > 1)
      {
        _gnuytics = _gnuytics + "set y2tics (";
        for (i in 1:x1.nr*x1.nc-1)
        { _gnuytics = _gnuytics + text(x1[i]) + ","; }
        _gnuytics = _gnuytics + text(x1[i+1]) + ");";
        else
      // x1 is a single value
          _gnuytics = _gnuytics + "set y2tics " + text(x1) + ";";
      }
    }
    if (class(x1) == "string")
    { _gnuytics = _gnuytics + x1 + ";"; }
    else
      _gnuytics = "unset ytics; ";
  }
  if (exist(x2))
  {
    if (class(x2) == "num")
    {
      // x2 is a single value
      _gnuytics = _gnuytics + "set my2tics " + text(x2) + ";";
    }
    if (class(x2) == "string")
    { _gnuytics = _gnuytics + x2; }
  else
    _gnuytics = _gnuytics + "unset my2tics;";
  }

  gnucmd(_gnuytics);
  return 1;
};



//
// gnulimits (xlo,xhi,ylo,yhi)
// gnulimits ( "set xrange [*:*]",,"set yrange [*:*]" )
// gnulimits (xlo,xhi,ylo,yhi,zlo,zhi) 
// gnulimits ( "set xrange [*:*]",, "set yrange [*:*]",, "set zrange [*:*]" )
//
gnulimits = function (xlo, xhi, ylo, yhi, zlo, zhi)
{
  I = gnuwin();
  if (isempty(I))
  { gnuwins(1); }

  s = "set xrange [$1:$2];";
  if (exist(xlo))
  {
    if (class(xlo) == "num")
    { s = gsub(text(xlo),"$1",s).string; }
    if (class(xlo) == "string")
    { s = xlo; }
  }
  s = gsub("*","$1",s).string;
  if (exist(xhi))
  {
    if (class(xhi) == "num")
    { s = gsub(text(xhi),"$2",s).string;}
  }
  s = gsub("*","$2",s).string;
  if (GNU_DEFAULTS.debug)
  { printf("gnulimits: %s\n", s); }
  gnucmd( s );
  
  s = "set yrange [$1:$2];";
  if (exist(ylo))
  {
    if (class(ylo) == "num")
    { s = gsub(text(ylo),"$1",s).string; }
    if (class(ylo) == "string")
    { s = ylo; }
  }
  s = gsub("*","$1",s).string;
  if (exist(yhi))
  {
    if (class(yhi) == "num")
    { s = gsub(text(yhi),"$2",s).string;}
  }
  s = gsub("*","$2",s).string;
  if (GNU_DEFAULTS.debug)
  { printf("gnulimits: %s\n", s); }
  gnucmd( s );

  if(exist(zlo)||exist(zhi))
  {
   s = "set zrange [$1:$2];";
   if (exist(zlo))
   {
     if (class(zlo) == "num")
     { s = gsub(text(zlo),"$1",s).string; }
     if (class(zlo) == "string")
     { s = zlo; }
   }
   s = gsub("*","$1",s).string;
   if (exist(zhi))
   {
     if (class(zhi) == "num")
     { s = gsub(text(zhi),"$2",s).string;}
   }
   s = gsub("*","$2",s).string;
   if (GNU_DEFAULTS.debug)
   { printf("gnulimits: %s\n", s); }
   gnucmd( s );
  }

  return 1;
};


gnulimits2 = function (xlo, xhi, ylo, yhi)
{
  I = gnuwin();
  if (isempty(I))
  { error("gnuplot: non-device accessed"); }

  if (exist(xlo) && exist(xhi))
  {
    if (class(xlo) == "num" && class(xhi) == "num")
    { _gnurange = "set x2range [" + text(xlo) + ":" + text(xhi) + "];"; }
    if (class(xlo) == "string" && class(xhi) == "string")
    {
      _gnurange = xlo[1] + xhi[1];
      return 1;
    }
    else
      _gnurange = "";
  }
  if (exist(ylo) && exist(yhi))
  {
    if (class(ylo) == "num" && class(yhi) == "num")
    {
      _gnurange = _gnurange + ...
          "set y2range [" + text(ylo) + ":" + text(yhi) + "];";
    }
  }

  gnucmd(_gnurange);

  return 1;
};


gnuxlabel = function ( x1, x2 )
{
  I = gnuwin();
  if (isempty(I))
  { error("gnuplot: non-device accessed"); }

  xlab = "set xlabel;";
  if (exist(x1))
  {
    if (class(x1)=="string")
    {
      if (strlen(x1[1])>0)
      { xlab = "set xlabel '" + x1[1] + "';" ; }
    }
  }

  x2lab = "unset x2label";
  if (exist(x2))
  {
    if (class(x2)=="string")
    {
      if (strlen(x2[1])>0)
      { x2lab = " set x2label '" + x2[1] + "'" ; }
    }
  }

  gnucmd(xlab + x2lab);
  return 1;
};



gnuylabel = function ( x1, x2 )
{
  I = gnuwin();
  if (isempty(I))
  { error("gnuplot: non-device accessed"); }

  ylab = "set ylabel;";
  if (exist(x1))
  {
    if (class(x1)=="string")
    {
      if (strlen(x1[1])>0)
      { ylab = "set ylabel '" + x1[1] + "';" ; }
    }
  }

  y2lab = "unset y2label";
  if (exist(x2))
  {
    if (class(x2)=="string")
    {
      if (strlen(x2[1])>0)
      { y2lab = " set y2label '" + x2[1] + "'" ; }
    }
  }

  gnucmd(ylab + y2lab);
  return 1;
};



gnuzlabel = function ( x1 )
{
  I = gnuwin();
  if (isempty(I))
  { error("gnuplot: non-device accessed"); }

  zlab = "set zlabel;";
  if (exist(x1))
  {
    if (class(x1)=="string")
    {
      if (strlen(x1[1])>0)
      { zlab = "set zlabel '" + x1[1] + "';" ; }
    }
  }

  gnucmd(zlab);
  return 1;
};



gnuformat = function ( keys )
{
  //global (_gnukey,GNU_DEFAULTS);

  if (!exist(keys))
  { return 0; }

  if (class(keys) != "string")
  { return 0; }
  
  if (!exist(_gnukey))
  { _gnukey = blank(0,0); }

  j = 1;
  for (i in keys)
  {
    _gnukey[;j] = i;
    j ++;
  }

  return 1;
};


gnulegend = function ( keys )
{
  //global (_gnulegend,GNU_DEFAULTS);

  if (!exist(keys))
  { return 0; }

  if (class(keys) != "string")
  { return 0; }
  
  if (!exist(_gnulegend))
  { _gnulegend = blank(0,0); }

  j = 1;
  for (i in keys)
  {
    _gnulegend[;j] = i;
    j ++;
  }
  return 1;
};


//
// gnuwins: manage the number of open windows
//
gnuwins = function (N, s, fn)
{
  // check N: not reading the manual makes user
  // VERY creative
  if (!exist(N))
  { return _gnuwins(); }

  if (class(N)!="real" && class(N)!="num")
  {
    printf("gnuwins: improper first argument");
    return 0;
  }
  if (N<0 || N>32)
  {
    printf("gnuwins: improper first argument");
    return 0;
  }

  if (!exist(fn))
  { fn = ""; }

  if (!exist(s))
  { s  = ""; }

  w  = gnuwins().available;
  nw = length (w);

  if (N > nw)
  {
    for (i in (nw+1):N)
    {
      if (strlen(fn)>0)
      {
        gnustart(,,fn);
      else
        gnustart();
      }

      if (strlen(s[min(i,length(s))]) > 0)
      { gnucmd(s[min(i,length(s))]); }
    }
    return 1;
  }

  w  = gnuwins().available;
  nw = length (w);

  // close extra windows from the end
  if (N < nw)
  {
    for (i in nw:N+1:-1)
    { gnuclose(i); }
  }

  // restart windows
  for (i in gnuwins().available)
  {
    gnuwin(i);
    if (strlen(s[min(i,length(s))]) > 0)
    {
      gnucmd(s[min(i,length(s))]);
    else
      gnucmd(GNU_DEFAULTS.resetcmd);
    }
  }
  return 1;
};


gnuplot = function (data, file, fmt)
{
  //global (_gnukey,_gnulegend);

  I = gnuwin();
  if (isempty(I))
  { gnuwins(1); }

  if (exist (file))
  {
    if (class(file) == "string")
    {
      if (!exist(fmt))
      {
        fmt = GNU_DEFAULTS.epscolor;
      else
        if (fmt=="color")
        {
          fmt = GNU_DEFAULTS.epscolor;
        }
        if (fmt=="bw" || fmt=="mono")
        {
          fmt = GNU_DEFAULTS.epsmono;
        }
      }

      // forward the output of the terminal to a file
      // only if talking to an open gnuplot device
      _i = gnuwin();
      _j = find(_i == gnuwins().available);
      if (isempty(_j))
      { error("gnuplot: terrible internal error. what the heck are you doing?"); }
      _d = gnuwins().device[_j];
      if (_d == "|gnuplot")
      {
        printf ("gnuplot: forwarding output to file '%s'", file);
        gnucmd ("set term " + fmt);
        gnucmd ("set output '" + file + "'");
      }
    }
  }

  //
  // legend and format for each dataset
  //
  plottitle = "";
  if (exist (_gnulegend))
  {
    if (!isempty(_gnulegend))
    {
      if (class(_gnulegend) == "string")
      { plottitle = _gnulegend; }
    }
  }

  plotopts = "";
  if (exist (_gnukey))
  {
    if (class(_gnukey)=="string")
    { plotopts = _gnukey; }
  }

  if (class(data) == "num")
  {
    // plot
    //   data = [x, y1, y2, ....]
    //
    // first dataset
    plotcmd = "plot '-'";
    if (strlen(plotopts[1])>0)
    { plotcmd = plotcmd + plotopts[1]; }
    if(strlen(plottitle[1])>0)
    {
      plotcmd = plotcmd + " title \"" + plottitle[1] + "\" ";
    else
      plotcmd = plotcmd + " notitle ";
    }

    // rest of data sets
    for (j in 3:data.nc)
    {
      plotcmd = plotcmd + ",'-'";

      // is there an options string for the data set
      if (length(plotopts)>=j-1)
      {
        if (strlen(plotopts[j-1])>0)
        {
          plotcmd = plotcmd + " " + plotopts[j-1] + " ";
        else
          plotcmd = plotcmd + " with lines";
        }
      else
        plotcmd = plotcmd + " with lines";
      }

      // is there a title for the data set
      if (length(plottitle)>=j-1)
      {
        if (strlen(plottitle[j-1])>0)
        {
          plotcmd = plotcmd + " title \"" + plottitle[j-1] + "\" ";
        else
          plotcmd = plotcmd + " notitle ";
        }
      else
        plotcmd = plotcmd + " notitle ";
      }

    }
    gnucmd (plotcmd);

    // send data
    if (data.nc >=2)
    {
      gnuprint( data[;1,2] );
      for (j in 3:data.nc)
      { gnuprint( data[;1,j] ); }
    }
  }

  if (class(data) == "list")
  {
    // plot
    //   data = << [x1, y11 ..] ; [x2, y21 ..]; .. >>
    //
    if (size(data) == 0)
    { return 0; }
    if (class(data.[members(data)[1]]) == "list")
    { error("gnuplot: don't know how to plot a sublist");}

    // first dataset
    if (class(data.[members(data)[1]]) == "num")
    { plotcmd = "plot '-'"; }
    if (class(data.[members(data)[1]]) == "string")
    {
      if (isfile(data.[members(data)[1]]))
      {
        plotcmd = "plot '" + data.[members(data)[1]][1] + "'";
      else
        plotcmd = "plot "  + data.[members(data)[1]][1] ;
      }
    }

    if (strlen(plotopts[1])>0)
    { plotcmd = plotcmd + " " + plotopts[1]; }

    if(strlen(plottitle[1])>0)
    {
      plotcmd = plotcmd + " title \"" + plottitle[1] + "\" ";
    else
      plotcmd = plotcmd + " notitle ";
    }

    // the other datasets
    for (j in 2:size(data))
    {

      if (class(data.[members(data)[j]]) == "num")
      { plotcmd = plotcmd + ",'-' "; }

      if (class(data.[members(data)[j]]) == "string")
      {
        if (isfile(data.[members(data)[j]]))
        {
          plotcmd = plotcmd + ",'" + data.[members(data)[j]][1] ...
              + "' ";
        else
          plotcmd = plotcmd + ", " + data.[members(data)[j]][1] ...
              + " ";
        }
      }
      // is there an options string for the data set
      if (length(plotopts)>=j)
      {
        if (strlen(plotopts[j])>0)
        {
          plotcmd = plotcmd + plotopts[j] + " ";
        else
          plotcmd = plotcmd + "with lines";
        }
      else
        plotcmd = plotcmd + "with lines";
      }
      // is there a title for the data set
      if (length(plottitle)>=j)
      {
        if (strlen(plottitle[j])>0)
        {
          plotcmd = plotcmd + " title \"" + plottitle[j] + "\" ";
        else
          plotcmd = plotcmd + " notitle ";
        }
      else
        plotcmd = plotcmd + " notitle ";
      }
    }

    gnucmd (plotcmd);
    // send data
    for (j in members(data))
    {
      if (class(data.[j])=="num")
      { gnuprint( data.[j] ); }
    }
  }

  if (class(data) == "string")
  {
    //
    // plot file or a result of a gnuplot command
    //
    if (length(data) == 0)
    { return 0; }

    // first entry
    if(isfile(data[1]))
    {
      plotcmd = "plot '" + data[1] + "' ";
    else
      plotcmd = "plot "  + data[1] + " " ;
    }
    if (strlen(plotopts[1])>0)
    { plotcmd = plotcmd + plotopts[1]; }
    if(strlen(plottitle[1])>0)
    {
      plotcmd = plotcmd + " title \"" + plottitle[1] + "\" ";
    else
      plotcmd = plotcmd + " notitle ";
    }

    // the rest of the entries
    for (j in 2:length(data))
    {
      if (isfile(data[j]))
      {
        plotcmd = plotcmd + ",'" + data[j] + "' " ;
      else
        plotcmd = plotcmd + ", " + data[j] + " "  ;
      }

      // is there an options string for the data set
      if (length(plotopts)>=j)
      {
        if (strlen(plotopts[j])>0)
        {
          plotcmd = plotcmd + " " + plotopts[j];
        else
          plotcmd = plotcmd + " " + "with lines";
        }
      else
        plotcmd = plotcmd + " " + "with lines";
      }
      // is there a title for the entry
      if (length(plottitle)>=j)
      {
        if (strlen(plottitle[j])>0)
        {
          plotcmd = plotcmd + " title \"" + plottitle[j] + "\" ";
        else
          plotcmd = plotcmd + " notitle ";
        }
      else
        plotcmd = plotcmd + " notitle ";
      }
    }

    gnucmd (plotcmd);
  }

  clear (_gnukey,_gnulegend);

  if (exist (file))
  {
    // reset the terminal
    _i = gnuwin();
    _j = find(_i == gnuwins().available);
    _d = gnuwins().device[_j];
    if (_d == "|gnuplot" && class(file) == "string")
    {
      printf (" .. is finished\n");
      gnucmd (GNU_DEFAULTS.reseterm); // linux
    }
  }

  return 1;
};


//
// gnusplot: surface plot
// gnusplot(data, file, fmt), for a list of surfaces, where
//    i in members(data), and
//    data.[i] = 'file name' or 'gnuplot compatible command' or
//    a list of surface points, data.[i].x (vector), data.[i].y (vector),
//    data.[i].z (matrix)
//
gnusplot = function (data, file, fmt)
{
  //global (_gnukey, _gnulegend);
  
  I = gnuwin();
  if (isempty(I))
  { gnuwins(1); }

  if (exist (file))
  {
    if (class(file) == "string")
    {
      if (!exist(fmt))
      {
        fmt = GNU_DEFAULTS.epscolor;
      else
        if (fmt=="color")
        {
          fmt = GNU_DEFAULTS.epscolor;
        }
        if (fmt=="bw" || fmt=="mono")
        {
          fmt = GNU_DEFAULTS.epsmono;
        }
      }
      // forward the output of the terminal to a file
      // only if talking to an open gnuplot device
      _i = gnuwin();
      _j = find(_i == gnuwins().available);
      if (isempty(_j))
      { error("gnuplot: terrible internal error. what the heck are you doing?"); }
      _d = gnuwins().device[_j];
      if (_d == "|gnuplot")
      {
        printf ("gnuplot: forwarding output to file '%s'", file);
        gnucmd ("set term " + fmt);
        gnucmd ("set output '"+ file + "'");
      }
    }
  }

  //
  // title and format for each dataset
  //
  plottitle = "";
  if (exist(_gnulegend))
  {
    if (class(_gnulegend) == "string")
    { plottitle = _gnulegend; }
  }
  plotopts = "";
  if (exist (_gnukey))
  {
    if (class(_gnukey)=="string")
    { plotopts = _gnukey; }
  }

  if (class(data) == "list")
  {
    // splot of
    //   <<x;y;z>>                                  // list
    //   data.[i].x, data.[i].y, data.[i].z         // list
    //   data.[i] = 'file name'                     // string
    //   data.[i] = 'gnuplot compatible command'    // string
    //
    if (size(data) == 0)
    { return 0; }

    // check for <<x;y;z>>
    if (sum(members(data)=="x") && sum(members(data)=="y") && sum(members(data)=="z"))
    { plotcmd = "splot '-' "; }

    if (class(data.[members(data)[1]]) == "list")
    { plotcmd = "splot '-' "; }

    if (class(data.[members(data)[1]]) == "string")
    {
      if (isfile(data.[members(data)[1]]))
      {
        plotcmd = "splot '" ...
            + data.[members(data)[1]][1] + "' ";
      else
        plotcmd = "splot " ...
            + data.[members(data)[1]][1] + " ";
      }
    }
    if (strlen(plotopts[1]) > 0)
    { plotcmd = plotcmd + " " + plotopts[1]; }
    if (strlen(plottitle[1]) > 0)
    {
      plotcmd = plotcmd + " title \"" + plottitle[1] + "\"";
    else
      plotcmd = plotcmd + " notitle";
    }

    for (j in 2:size(data))
    {
      // check for <<x;y;z>>, if so stop immediatelly
      if (sum(members(data)=="x") && sum(members(data)=="y") && sum(members(data)=="z"))
      { break; }
        
      if (class(data.[members(data)[j]]) == "list")
      { plotcmd = plotcmd + ",'-' "; }
      if (class(data.[members(data)[j]]) == "string")
      {
        if (isfile(data.[members(data)[j]]))
        {
          plotcmd = plotcmd + ",'" + data.[members(data)[j]][1] + "' ";
        else
          plotcmd = plotcmd + ", " + data.[members(data)[j]][1] + " ";
        }
      }
      
      if (length(plotopts) >= j)
      {
        if(strlen(plotopts[j])>0)
        { plotcmd = plotcmd + plotopts[ j ]; }
      }
      if (length(plottitle) >= j)
      {
        if(strlen(plottitle[j])>0)
        {
          plotcmd = plotcmd + " title \"" + plottitle[j] + "\"";
        else
          plotcmd = plotcmd + " notitle";
        }
      else
        plotcmd = plotcmd + " notitle";
      }
    }

    if (GNU_DEFAULTS.debug)
    { printf("gnusplot: %s\n", plotcmd); }

    gnucmd (plotcmd);
    // send data
    for (j in members(data))
    {
      // <<x;y;z>>
      if (sum(members(data)=="x") && sum(members(data)=="y") && sum(members(data)=="z"))
      {
        if (exist(data.x) && exist(data.y) && exist(data.z))
        {
          gnuprint( data.x, data.y, data.z );
        else
          printf("gnusplot: <<x;y;z>> list is required\n");
          error ("gnusplot: cannot print!")
        }
        break;
      }
      // else
      if (class(data.[j])=="list")
      {
        if (exist(data.[j].x) && exist(data.[j].y) && exist(data.[j].z))
        {
          gnuprint( data.[j].x, data.[j].y, data.[j].z );
        else
          printf("Improperly formatted list!\n");
          error ("gnusplot: cannot print!")
        }
      }
    }
  }

  if (class(data) == "string")
  {
    if (length(data) == 0)
    { return 0; }

    if(isfile(data[1]))
    {
      plotcmd = "splot '" + data[1] + "' ";
    else
      plotcmd = "splot "  + data[1] + " ";
    }
    if (strlen(plotopts[1]) > 0)
    { plotcmd = plotcmd + " " + plotopts[1]; }
    if (strlen(plottitle[1]) > 0)
    {
      plotcmd = plotcmd + " title \"" + plottitle[1] + "\"";
    else
      plotcmd = plotcmd + " notitle";
    }

    for (j in 2:length(data))
    {
      if (isfile(data[j]))
      {
        plotcmd = plotcmd + ",'" + data[j] + "' ";
      else
        plotcmd = plotcmd + ", " + data[j] + " " ;
      }

      if (length(plotopts) >= j)
      {
        if(strlen(plotopts[j])>0)
        { plotcmd = plotcmd + plotopts[ j ]; }
      }
      if (length(plottitle) >= j)
      {
        if(strlen(plottitle[j])>0)
        {
          plotcmd = plotcmd + " title \"" + plottitle[j] + "\"";
        else
          plotcmd = plotcmd + " notitle";
        }
      else
        plotcmd = plotcmd + " notitle";
      }
    }
    gnucmd (plotcmd);
  }

  clear (_gnukey,_gnulegend);

  if (exist (file))
  {
    // reset the terminal
    _i = gnuwin();
    _j = find(_i == gnuwins().available);
    _d = gnuwins().device[_j];
    if (_d == "|gnuplot" && class(file) == "string")
    {
      printf (" .. is finished\n");
      gnucmd (GNU_DEFAULTS.reseterm); // linux
    }
  }

  return 1;
};



