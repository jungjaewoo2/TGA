// Copyright (C) 2003-2016 Marijan Kostrun
//   part of rlabplus for linux project on rlabplus.sourceforge.net
//
// gnuplot through pipe for rlabplus, this complements
// the functions in gnuplot.c / gnuplot.h in source directory
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
//
// See the file ./COPYING
//
require libstdio libsystem

static(_gnukey, _gnugrid, _gnulegend, GNU_DEFAULTS, _gnucolor_replace, PALETTE);
PALETTE = "/usr/share/X11/rgb.txt";

GNU_DEFAULTS = <<>>;
GNU_DEFAULTS.debug    = 0;  // say what you are doing at every step
GNU_DEFAULTS.info     = 1;  // just the major messages, e.g., forwarding plot to a file
if (isfile(PALETTE))
{
  GNU_DEFAULTS.palette  = PALETTE;
}
else
{
  GNU_DEFAULTS.palette  = "default";
}


// gnuwins, reset terminal after each call:
GNU_DEFAULTS.reset_wins = "reset;";

// reset output, after forwarding output to a file:
GNU_DEFAULTS.reset_term = "reset; set term x11 enhanced; set output";

// default terminal for redirection 'set output'
GNU_DEFAULTS.output = "eps";

static(gnu_object_count, gnu_arrow_count);
gnu_object_count = 0;
gnu_arrow_count = 0;

// font acrobatics: png supports only pfa and ttf fonts, while eps supports
// all the fonts that ghostscript supports. thus:
GNU_DEFAULTS.fontsize = 11;
GNU_DEFAULTS.font = <<>>;
GNU_DEFAULTS.font.eps = "Times-Roman";
GNU_DEFAULTS.font.png = "times";
GNU_DEFAULTS.font_symbol.eps = "Symbol";
GNU_DEFAULTS.font_symbol.png = "symbol";

// size of the plot terminal
GNU_DEFAULTS.term_size = <<>>;
GNU_DEFAULTS.term_size.eps = "6.5,4";   // in inches, m.k. preference
// GNU_DEFAULTS.term_size.eps = "5,3.5";   // in inches, gnuplot default
GNU_DEFAULTS.term_size.ps  = "10,7";    // in inches
GNU_DEFAULTS.term_size.png = "800,600"; // in pixels - it is a bitmap plot
GNU_DEFAULTS.term_size.jpeg = "640,480"; // in pixels - it is a bitmap plot

// set terminal command for different terminals:
GNU_DEFAULTS.set_term = <<>>;

GNU_DEFAULTS.set_term.eps = "postscript eps enh color dashed colortext lw 2 " ...
    + "\"" + GNU_DEFAULTS.font.eps + "\"" + num2str(2*GNU_DEFAULTS.fontsize, " %g") ...
    + " size " + GNU_DEFAULTS.term_size.eps;
GNU_DEFAULTS.set_term.eps_mono  = "postscript eps enh mono  dashed blacktext lw 2 " ...
    + "\"" + GNU_DEFAULTS.font.eps + "\"" + num2str(2*GNU_DEFAULTS.fontsize, " %g") ...
    + " size " + GNU_DEFAULTS.term_size.eps;
GNU_DEFAULTS.set_term.ps = "postscript enh color dashed colortext lw 2 " ...
    + "\"" + GNU_DEFAULTS.font.eps + "\"" + num2str(2*GNU_DEFAULTS.fontsize, " %g") ...
    + " size " + GNU_DEFAULTS.term_size.ps;
GNU_DEFAULTS.set_term.ps_mono  = "postscript enh mono  dashed blacktext lw 2 " ...
    + "\"" + GNU_DEFAULTS.font.eps + "\"" + num2str(2*GNU_DEFAULTS.fontsize, " %g") ...
    + " size " + GNU_DEFAULTS.term_size.ps;
GNU_DEFAULTS.set_term.png_t  = "png transparent nointerlace truecolor butt enh" ...
    + " font \"" + GNU_DEFAULTS.font.png + "\" " + num2str(2*GNU_DEFAULTS.fontsize,"%g") ...
    + " size " + GNU_DEFAULTS.term_size.png ...
    + " nocrop x000000 x202020 x404040 x606060 x808080 xA0A0A0 xC0C0C0 xE0E0E0";
GNU_DEFAULTS.set_term.png = "png notransparent nointerlace truecolor butt enh" ...
    + " font \"" + GNU_DEFAULTS.font.png + "\" " + num2str(2*GNU_DEFAULTS.fontsize,"%g") ...
    + " size " + GNU_DEFAULTS.term_size.png ...
    + " nocrop xffffff x000000 x404040 xff0000 xffa500 x66cdaa xcdb5cd" ...
    + " xadd8e6 x0000ff xdda0dd x9500d3";
GNU_DEFAULTS.set_term.jpeg = "jpeg enhanced normal" ...
    + " size " + GNU_DEFAULTS.term_size.jpeg ...
    + "font 'arial,12' fontscale 1.0" ...
    + " nocrop xffffff x000000 x404040 xff0000 xffa500 x66cdaa xcdb5cd" ...
    + " xadd8e6 x0000ff xdda0dd x9500d3";

// environment variables gnuplot might need for png
GNU_DEFAULTS.env = <<>>;
GNU_DEFAULTS.env.GDFONTPATH="/usr/share/fonts/truetype:/usr/share/ghostscript/fonts";

GNU_DEFAULTS.style_arrow = [ ...
  "head filled size screen 0.025,30,45 ls 1", ...
  "head nofilled size screen 0.03,15 ls 2", ...
  "head filled size screen 0.03,15,45 ls 1", ...
  "head filled size screen 0.03,15 ls 2", ...
  "heads filled size screen 0.03,15,135 ls 1", ...
  "head empty size screen 0.03,15,135 ls 2", ...
  "nohead ls 1", ...
  "heads size screen 0.008,90 ls 2", ...
[]];

// other setting for gnuplot
GNU_DEFAULTS.histmode = 1;

//
// configure environment variables for gnuplot
//
for (__i in members(GNU_DEFAULTS.env))
{
  if(getenv(__i)=="" && strlen(GNU_DEFAULTS.env.[__i])>0)
  {
    __cmd = __i + "=" + GNU_DEFAULTS.env.[__i];
    putenv(__cmd);
  }
}

//
// check if user has symbol.ttf installed for using with png
//
if(GNU_DEFAULTS.font_symbol.png == "")
{
  __i  = [0, findstr(GNU_DEFAULTS.env.GDFONTPATH,":"), ...
      strlen(GNU_DEFAULTS.env.GDFONTPATH)+1];
  for (__j in 2:length(__i))
  {
    __i1  = __i[__j -1] + 1;
    __i2  = __i[__j   ] - 1;
    __dir = substr(GNU_DEFAULTS.env.GDFONTPATH, [__i1:__i2]);
    __x   = ls(__dir);
    __X   = toupper(__x);
    __if  = find(__X == "SYMBOL.TTF");
    if (__if.n)
    {
      __fontname = substr(__x[__if], 1:(findstr(__x[__if],".")-1) );
      GNU_DEFAULTS.font_symbol.png = __fontname;
      if (GNU_DEFAULTS.debug)
      {
        printf (THIS_LIB + ": font %s in directory %s can be used for 'png' output\n", ...
            __x[__if], __dir);
        printf (THIS_LIB + ": e.g., as follows:\n");
        printf (THIS_LIB + ": {/%s 'char'}, or {/%s='size' 'char'}\n", ...
            __fontname, __fontname);
      }
    }
    else if (GNU_DEFAULTS.debug)
    {
      printf (THIS_LIB + ": No symbol font available for png terminal. Please find\n");
      printf (THIS_LIB + ": the font 'symbol.ttf' on the web and install it together\n");
      printf (THIS_LIB + ": with the libgd.\n");
    }
  }
}

//
// default palette for gnuplot
//
clear(GNU_DEFAULTS.color);
if (GNU_DEFAULTS.palette != "default")
{
  if (isfile(GNU_DEFAULTS.palette))
  {
    // list, colorname = <<space;xyz;hex>>
    GNU_DEFAULTS.color  = <<>>;
    __s = reads(GNU_DEFAULTS.palette);
    for(__i in length(__s):1:-1)
    {
      if(strlen(__s[__i]))
      { break; }
    }
    if (__i > 1)
    {
      // format of X11's rgb.txt file is
      // R G B rgbColorName
      GNU_DEFAULTS.color = <<>>;
      for (__j in 1:__i)
      {
        __y = strsplt(__s[__j],"'BLANK"); // since 2015-5-20/2.4.2.2
        if(length(__y)>=4)
        {
          __z = strtod(__y[1:3]);
          __name = tolower(__y[4]);
          if (length(__y)>4)
          {
            for (__k in 5:length(__y))
            { __name = __name + tolower(__y[__k]); }
          }
          if (exist(GNU_DEFAULTS.color.[__name]))
          { continue; }

          GNU_DEFAULTS.color.[__name] = <<>>;
          GNU_DEFAULTS.color.[__name].xyz   = __z;
          __no = int(65536 * __z[1] + 256 * __z[2] + __z[3]);
          GNU_DEFAULTS.color.[__name].hex   = toupper( num2str(__no, "\"#%x\"") );
          GNU_DEFAULTS.color.[__name].space = "rgb";
        }
      }
    }
  }
  else
  {
    printf(THIS_LIB + ": library is set to use X11 color specifications from a\n");
    printf(THIS_LIB + ": file  %s\n", GNU_DEFAULTS.palette);
    printf(THIS_LIB + ": but the file is not present!\n");
    printf(THIS_LIB + ": Please set 'GNU_DEFAULTS.palette' to 'default' to avoid\n");
    printf(THIS_LIB + ": seeing this message every time RLaB2 starts.\n");
  }
}

// clear all local variables before defining the
// library functions
for(__i in grep(members($$)',"__"))
{ clear($$.[__i]); }
clear(__i);

//
// replace 'lc @colorname@' with the color specification
// from GNU_DEFAULTS.color.[colorname]
// replace 'tc @colorname@' with the color specification
// from GNU_DEFAULTS.color.[colorname]
//
_gnucolor_replace = function( i, cmd )
{
  if (cmd == "gnuformat")
  {
    x1 = findstr(i, "lc @");  // do replacement only if near 'lc' command
    x2 = findstr(i, "lc  @"); // user typed an extra space by mistake. silly user!
  }
  else if (cmd == "gnutext")
  {
    x1 = findstr(i, "tc @");  // do replacement only if near 'tc' command
    x2 = findstr(i, "tc  @"); // user typed an extra space by mistake. silly user!
  }
  x  = findstr(i, "@");
  while (!isempty(x) && (!isempty(x1) || !isempty(x2)))
  {
    i1 = x[1];
    i2 = x[2];
    s  = substr(i, i1:i2);
    name = tolower(substr(i, (i1+1):(i2-1)));
    if (exist(GNU_DEFAULTS.color.[name]))
    {
      r = GNU_DEFAULTS.color.[name].space + " " ...
          + GNU_DEFAULTS.color.[name].hex;
    }
    else
    {
        r = "rgb 'black'";
    }
    i = gsub(r, s, i).string;
    x = findstr(i, "@");
  }
  return i;
};


//
//
// scripted functions for manipulations with gnuplot
//
//
gnudefault = function( x, val)
{
  if (!exist(x))
  { return GNU_DEFAULTS; }

  if (exist(GNU_DEFAULTS.[x]))
  {
    if (exist(val))
    {
      GNU_DEFAULTS.[x] = val;
    }

    return GNU_DEFAULTS.[x];
  }

  return [];
};

gnutext = function (s, loc, vel, cs, opts )
{
  _this_function = "gnutext";
  append = "";

  I = gnuwin();
  if (isempty(I))
  { gnuwins(1); }

  if (!exist(s) || !exist(loc))
  {
    gnucmd ("set label",_this_function);
    return 0;
  }

  if (!exist(cs))
  { cs = "first"; }

  if ( length(loc) != 2 )
  { return 0; }

  // location
  cmd = "set label '" + s + "' at " + cs + " ";
  cmd = cmd + num2str(loc[1]) + "," + num2str(loc[2]);

  if (exist(vel))
  {
    cmd = cmd + " font \"" + GNU_DEFAULTS.font.[GNU_DEFAULTS.output] ...
        + "," + num2str(2*vel*GNU_DEFAULTS.fontsize, "%g\"");
  }

  // fontname,fontsize
  if (exist(opts))
  {
    if (class(opts)=="list")
    {
      if (exist(opts.fontname))
      { fontname = opts.fontname; }

      if (exist(opts.fontsize))
      { fontsize = opts.fontsize; }

      if (exist(opts.append))
      { append = opts.append; }
      cmd = cmd + " font \"" + GNU_DEFAULTS.font.[GNU_DEFAULTS.output] ...
            + "," + num2str(2*vel*GNU_DEFAULTS.fontsize, "%g\"");
    }

    // if opts is a string just add it at the end
    if (exist(opts) && class(opts)=="string")
    {
      append = opts;
      append = _gnucolor_replace(append, "gnutext");
      cmd = cmd + " " + append;
    }
  }

  gnucmd( cmd, _this_function);
  return 1;
};

gnucir = function (center, r, a, cs, opts)
{
  _this_function = "gnucir";
  I = gnuwin();
  if (isempty(I))
  { gnuwins(1); }

  if ( center.nc != 2)
  { return 1; }

  // if opts is a string just add it at the end
  if (exist(opts) && class(opts)=="string")
  {
    a = " " + opts;
  }
  else
  {
    a = "";
  }

  if (!exist(cs))
  { cs = ""; }

  // location
  nr = max(center.nr, r.nr);
  for (i in 1:nr)
  {
    i1 = min(i, center.nr);
    i2 = min(i, r.nr);
    ia = min(i, a.n);

    gnu_object_count++;

    cmd = "set object " + num2str(gnu_object_count,"%.0f") + " circle at " ...
        + num2str(center[i1;],"%g",",") + " size " + cs + num2str(r[i2]," %g") + " " + a[ia];

    gnucmd( cmd, _this_function);
  }

  return 0;
};

gnurect = function (loc1, loc2, opts )
{
  _this_function = "gnurect";
  I = gnuwin();
  if (isempty(I))
  { gnuwins(1); }

  if ( loc1.nc != 2  || loc2.nc != 2)
  { return 1; }

  // if opts is a string just add it at the end
  if (exist(opts) && class(opts)=="string")
  {
    a = " " + opts;
  }
  else
  {
    a = "";
  }

  // location
  nr = max(loc1.nr, loc2.nr);
  for (i in 1:nr)
  {
    i1 = min(i, loc1.nr);
    i2 = min(i, loc2.nr);
    ia = min(i, a.n);

    gnu_object_count++;

    cmd = "set object " + num2str(gnu_object_count,"%.0f") + " rect from " ...
      + num2str(loc1[i1;],"%g",",") + " to " + num2str(loc2[i2;],"%g",",") + " " + a[ia];

    gnucmd( cmd, _this_function);
  }

  return 0;
};


gnuarrow = function (loc1, loc2, opts )
{
  _this_function = "gnuarrow";
  I = gnuwin();
  if (isempty(I))
  { gnuwins(1); }

  if ( loc1.nc != 2  || loc2.nc != 2)
  { return 1; }

  // if opts is a string just add it at the end
  if (exist(opts) && class(opts)=="string")
  { a = " " + opts; }

  // location
  nr = max(loc1.nr, loc2.nr);
  for (i in 1:nr)
  {
    i1 = min(i, loc1.nr);
    i2 = min(i, loc2.nr);

    gnu_arrow_count++;

    if (exist(a))
    {
      ia = min(i, a.n);
      sa = a[ia];
    }
    else
    {
      ia = mod(gnu_arrow_count, GNU_DEFAULTS.style_arrow.n);
      if (!ia)
      { ia = GNU_DEFAULTS.style_arrow.n; }
      sa = GNU_DEFAULTS.style_arrow[ ia ];
    }

    if (strindex(sa, " as ")>0)
    {
      cmd = "set arrow from " ...
        + num2str(loc1[i1;],"%g",",") + " to " + num2str(loc2[i2;],"%g",",") + " " + sa ;
    }
    else
    {
      cmd = "set arrow " + sa + " from " ...
        + num2str(loc1[i1;],"%g",",") + " to " + num2str(loc2[i2;],"%g",",");
    }

    gnucmd( cmd, _this_function);
  }

  return 0;
};

gnutitle = function ( s )
{
  _this_function = "gnutitle";
  I = gnuwin();
  if (isempty(I))
  { gnuwins(1); }

  cmd = "set title;";
  if (exist(s))
  {
    if (strlen(s[1]) > 0)
    { cmd = "set title '" + s[1] + "';"; }
  }

  gnucmd( cmd, _this_function);
  return 1;
};

gnuscale2 = function ( sx, sy )
{
  _this_function = "gnuscale2";
  I = gnuwin();
  if (isempty(I))
  { gnuwins(1); }

  cmd = "unset logscale x2;";
  if (exist(sx))
  {
    if (strindex(tolower(sx),"log"))
    { cmd = "set logscale x2;"; }
  }
  gnucmd( cmd, _this_function);

  cmd = "unset logscale y2;";
  if (exist(sy))
  {
    if (strindex(tolower(sy),"log"))
    { cmd = "set logscale y2;"; }
  }
  gnucmd( cmd, _this_function);

  return 1;
};


gnuscale = function ( sx, sy, sz )
{
  _this_function = "gnuscale";
  I = gnuwin();
  if (isempty(I))
  { gnuwins(1); }

  cmd = "unset logscale x;";
  if (exist(sx))
  {
    if (strindex(tolower(sx),"log"))
    { cmd = "set logscale x;"; }
  }
  
  gnucmd( cmd, _this_function);

  cmd = "unset logscale y;";
  if (exist(sy))
  {
    if (strindex(tolower(sy),"log"))
    { cmd = "set logscale y;"; }
  }
  gnucmd( cmd, _this_function);

  if (exist(sz))
  {
    cmd = "unset logscale z;";
    if (exist(sz))
    {
      if (strindex(tolower(sz),"log"))
      { cmd = "set logscale z"; }
    }
    gnucmd( cmd, _this_function);
  }

  return 1;
};

gnuztics = function ( x1, x2, fmt )
{
  _this_function = "gnuztics";
  I = gnuwin();
  if (isempty(I))
  { gnuwins(1); }

  if (!exist(fmt))
  { fmt = "%g"; }
  if (class(fmt)!="string")
  { fmt = "%g"; }

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
        { _gnuztics = _gnuztics + num2str(x1[i]) + ","; }
        _gnuztics = _gnuztics + num2str(x1[i+1]) + ");";
      }
      else
      {
        // x1 is a single value
        _gnuztics = _gnuztics + "set ztics " + num2str(x1) + ";";
      }
    }
    if (class(x1) == "string")
    { _gnuztics = _gnuztics + x1 + ";"; }
  }
  else
  {
    _gnuztics = _gnuztics + "set ztics;";
  }

  if (exist(x2))
  {
    if (class(x2) == "num")
    {
      // x2 is a single value
      _gnuztics = _gnuztics + "set mztics " + num2str(x2) + ";";
    }
    if (class(x2) == "string")
    { _gnuztics = _gnuztics + x2 +";"; }
  }
  else
  {
    _gnuztics = _gnuztics + "unset mztics;";
  }

  gnucmd( cmd, _this_function);
  return 1;
};


gnuxtics = function ( a1, a2, fmt, opt )
{
  _this_function = "gnuxtics";
  I = gnuwin();
  if (isempty(I))
  { error(_this_function+": non-device accessed"); }

  // option to be added at the end of the 'xtics ....'
  if (!exist(opt))
  { opt = ""; }
  if (class(opt)!="string")
  { opt = "mirror"; }

  _gnuxtics = "unset xtics;";
  if (exist(a1))
  {
    if (!exist(fmt))
    { fmt = "%g"; }
    if (class(fmt)!="string")
    { fmt = "%g"; }

    _gnuxtics = _gnuxtics + "set format x '" + fmt[1] + "'; ";
    _gnuxtics = _gnuxtics + "set xtics " + opt[1] + " ";
    if (class(a1) == "num")
    {
      // a1 is a numeric vector
      if (length(a1) > 1)
      {
        // a1 is a vector
        _gnuxtics = _gnuxtics + "(" + num2str(a1,,",") + ");";
      }
      else
      {
        // a1 is a single value
        _gnuxtics = _gnuxtics + num2str(a1) + ";";
      }
    }
    else if (class(a1) == "string")
    {
      // a1 is a string
      printf("gnuplot: (gnuxtics) Function called with string argument. Use 'gnucmd' instead.\n");
      error("");
    }
  }

  _gnuxtics = _gnuxtics + "unset mxtics; ";
  if (exist(a2))
  {
    if (class(a2) == "num")
    {
      // x2 is a single value
      _gnuxtics = _gnuxtics + "set mxtics " + num2str(a2) + ";";
    }
    else if (class(a2) == "string")
    {
      // a2 is a string, allow user to pass it along
      _gnuxtics = _gnuxtics + a2 +";";
    }
  }

  gnucmd(_gnuxtics, _this_function);
  return 1;
};


gnux2tics = function ( x1, x2, fmt )
{
  _this_function = "gnux2tics";
  I = gnuwin();
  if (isempty(I))
  { error(_this_function+": non-device accessed"); }

  if (!exist(fmt))
  { fmt = "%g"; }
  if (class(fmt)!="string")
  { fmt = "%g"; }

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
        { _gnuxtics = _gnuxtics + num2str(x1[i]) + ","; }
        _gnuxtics = _gnuxtics + num2str(x1[i+1]) + ");";
      }
      else
      {
          // x1 is a single value
          _gnuxtics = _gnuxtics + "set x2tics " + num2str(x1) + ";";
      }
    }
    if (class(x1) == "string")
    { _gnuxtics = _gnuxtics + x1 + ";"; }
  }
  else
  {
    _gnuxtics = _gnuxtics + "unset x2tics; ";
  }

  if (exist(x2))
  {
    if (class(x2) == "num")
    {
      // x2 is a single value
      _gnuxtics = _gnuxtics + "set mx2tics " + num2str(x2) + ";";
    }
    if (class(x2) == "string")
    { _gnuxtics = _gnuxtics ";" + x2; }
  }
  else
  {
    _gnuxtics = _gnuxtics + ";" + "unset mx2tics;";
  }

  gnucmd(_gnuxtics, _this_function);
  return 1;
};


gnuytics = function ( y1, y2, fmt, opt )
{
  _this_function = "gnuytics";
  I = gnuwin();
  if (isempty(I))
  { error(_this_function+": non-device accessed"); }

  // option to be added at the end of the 'ytics ....'
  if (!exist(opt))
  { opt = ""; }
  if (class(opt)!="string")
  { opt = "mirror"; }

  _gnuytics = "unset y2tics; ";
  if (exist(y1))
  {
    // format of the tics
    if (!exist(fmt))
    { fmt = "%g"; }
    if (class(fmt)!="string")
    { fmt = "%g"; }

    _gnuytics = _gnuytics + "set format y '" + fmt[1] + "'; ";
    _gnuytics = _gnuytics + "set ytics " + opt[1] + " ";
    if (class(y1) == "num")
    {
      if (length(y1) > 1)
      {
        // y1 is a vector
        _gnuytics = _gnuytics + "(" + num2str(y1,,",") + ");";
      }
      else
      {
        // y1 is a single value
        _gnuytics = _gnuytics + num2str(y1) + ";" ;
      }
    }
    if (class(y1) == "string")
    {
      printf("gnuplot: (gnuytics) Function called with string argument. Use 'gnucmd' instead.\n");
      error("");
    }
  }

  _gnuytics = _gnuytics + "unset mytics; ";
  if (exist(y2))
  {
    if (class(y2) == "num")
    {
      // y2 is a single value
      _gnuytics = _gnuytics + "set mytics " + num2str(y2) + ";";
    }
    else if (class(y2) == "string")
    {
      _gnuytics = _gnuytics + y2;
    }
  }

  gnucmd(_gnuytics, _this_function);
  return 1;
};


//
// gnuy2tics ([y1,..yK], DY)
// gnuy2tics ("set y2tics (*)", "set my2tics *")
//
gnuy2tics = function ( x1, x2, fmt )
{
  _this_function = "gnuy2tics";
  I = gnuwin();
  if (isempty(I))
  { error(_this_function+": non-device accessed"); }

  if (!exist(fmt))
  { fmt = "%g"; }
  if (class(fmt)!="string")
  { fmt = "%g"; }

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
        { _gnuytics = _gnuytics + num2str(x1[i]) + ","; }
        _gnuytics = _gnuytics + num2str(x1[i+1]) + ");";
      }
      else
      {
        // x1 is a single value
        _gnuytics = _gnuytics + "set y2tics " + num2str(x1) + ";";
      }
    }
    if (class(x1) == "string")
    { _gnuytics = _gnuytics + x1 + ";"; }
  }
  else
  {
    _gnuytics = "unset ytics; ";
  }
  if (exist(x2))
  {
    if (class(x2) == "num")
    {
      // x2 is a single value
      _gnuytics = _gnuytics + "set my2tics " + num2str(x2) + ";";
    }
    if (class(x2) == "string")
    { _gnuytics = _gnuytics + x2; }
  }
  else
  {
    _gnuytics = _gnuytics + "unset my2tics;";
  }

  gnucmd(_gnuytics, _this_function);
  return 1;
};



//
// gnulimits (xlo,xhi,ylo,yhi)
// gnulimits (xlo,xhi,ylo,yhi,zlo,zhi)
//
gnulimits = function (xlo, xhi, ylo, yhi, zlo, zhi)
{
  _this_function = "gnulimits";
  I = gnuwin();
  if (isempty(I))
  { gnuwins(1); }

  // each axis is dealt with separately:

  // x-axis: reset it first
  s = "unset xrange; set xrange [$1:$2];";
  if (exist(xlo))
  {
    if (class(xlo) == "num")
    {
      s = gsub(text(xlo[1]),"$1",s).string;
    }
    else if (class(xlo[1]) == "string")
    {
      printf("gnuplot: (gnulimits) Function called with string argument. Use 'gnucmd' instead.\n");
      error("");
    }
  }
  if (exist(xhi))
  {
    if (class(xhi) == "num")
    {
      s = gsub(text(xhi[1]),"$2",s).string;
    }
    else
    {
      printf("gnuplot: (gnulimits) Function called with string argument. Use 'gnucmd' instead.\n");
      error("");
    }
  }
  s = gsub("*","$1",s).string;
  s = gsub("*","$2",s).string;
  gnucmd( s, _this_function);

  // y-axis: reset it first
  s = "unset yrange; set yrange [$1:$2];";
  if (exist(ylo))
  {
    if (class(ylo) == "num")
    {
      s = gsub(text(ylo[1]),"$1",s).string;
    }
    else
    {
      printf("gnuplot: (gnulimits) Function called with string argument. Use 'gnucmd' instead.\n");
      error("");
    }
  }
  if (exist(yhi))
  {
    if (class(yhi) == "num")
    { s = gsub(text(yhi[1]),"$2",s).string;}
  }
  s = gsub("*","$1",s).string;
  s = gsub("*","$2",s).string;
  gnucmd( s, _this_function);

  // if z-axis is present, reset it first
  if(exist(zlo)||exist(zhi))
  {
   s = "unset zrange; set zrange [$1:$2];";
   if (exist(zlo))
   {
     if (class(zlo) == "num")
     {
       s = gsub(text(zlo[1]),"$1",s).string;
     }
     else
     {
       printf("gnuplot: (gnulimits) Function called with string argument. Use 'gnucmd' instead.\n");
       error("");
     }
   }
   if (exist(zhi))
   {
     if (class(zhi) == "num")
     {
       s = gsub(text(zhi[1]),"$2",s).string;
     }
     else
     {
       printf("gnuplot: (gnulimits) Function called with string argument. Use 'gnucmd' instead.\n");
       error("");
     }
   }
   s = gsub("*","$1",s).string;
   s = gsub("*","$2",s).string;
   gnucmd( s, _this_function);
  }

  return 1;
};


gnulimits2 = function (xlo, xhi, ylo, yhi)
{
  _this_function = "gnulimits2";
  I = gnuwin();
  if (isempty(I))
  { error(_this_function+": non-device accessed"); }

  if (exist(xlo) && exist(xhi))
  {
    if (class(xlo) == "num" && class(xhi) == "num")
    { _gnurange = "set x2range [" + num2str(xlo) + ":" + num2str(xhi) + "];"; }
    if (class(xlo) == "string" && class(xhi) == "string")
    {
      _gnurange = xlo[1] + xhi[1];
      return 1;
    }
  }
  else
  {
    _gnurange = "";
  }
  if (exist(ylo) && exist(yhi))
  {
    if (class(ylo) == "num" && class(yhi) == "num")
    {
      _gnurange = _gnurange + ...
          "set y2range [" + num2str(ylo) + ":" + num2str(yhi) + "];";
    }
  }

  gnucmd(_gnurange, _this_function);
  return 1;
};


gnuxlabel = function ( x1, x2 )
{
  _this_function = "gnuxlabel";
  I = gnuwin();
  if (isempty(I))
  { error(_this_function+": non-device accessed"); }

  xlab = "set xlabel;";
  if (exist(x1))
  {
    if (class(x1)=="string")
    {
      x1[1] = gsub("'", "\"", x1[1]).string;
      if (strlen(x1[1])>0)
      { xlab = "set xlabel \"" + x1[1] + "\";" ; }
    }
  }

  x2lab = "unset x2label";
  if (exist(x2))
  {
    if (class(x2)=="string")
    {
      x2[1] = gsub("\"", "'", x2[1]).string;
      if (strlen(x2[1])>0)
      { x2lab = " set x2label \"" + x2[1] + "\"" ; }
    }
  }

  gnucmd(xlab + x2lab, _this_function);
  return 1;
};



gnuylabel = function ( x1, x2 )
{
  _this_function = "gnuylabel";
  I = gnuwin();
  if (isempty(I))
  { error(_this_function+": non-device accessed"); }

  ylab = "set ylabel;";
  if (exist(x1))
  {
    if (class(x1)=="string")
    {
      if (strlen(x1[1])>0)
      { ylab = "set ylabel \"" + x1[1] + "\";" ; }
    }
  }

  y2lab = "unset y2label";
  if (exist(x2))
  {
    if (class(x2)=="string")
    {
      if (strlen(x2[1])>0)
      { y2lab = " set y2label \"" + x2[1] + "\";" ; }
    }
  }

  gnucmd(ylab + y2lab, _this_function);
  return 1;
};



gnuzlabel = function ( x1 )
{
  _this_function = "gnuzlabel";
  I = gnuwin();
  if (isempty(I))
  { error(_this_function+": non-device accessed"); }

  zlab = "set zlabel;";
  if (exist(x1))
  {
    if (class(x1)=="string")
    {
      if (strlen(x1[1])>0)
      { zlab = "set zlabel '" + x1[1] + "';" ; }
    }
  }

  gnucmd(zlab, _this_function);
  return 1;
};


gnuformat = function ( keys )
{
  //global (_gnukey,GNU_DEFAULTS);
  _this_function = "gnuformat";

  if (!exist(keys))
  { return 0; }

  if (class(keys) != "string")
  { return 0; }

  _gnukey = blank(0,0);
  j = 1;
  for (i in keys)
  {
    // check if 'i' contains statement
    //  ... lc @colorname@ ....
    // and if does, then do the replacement using
    // ...  lc 'colorspace' '#XXXXXX' ....
    if (exist(GNU_DEFAULTS.color))
    { i = _gnucolor_replace(i, "gnuformat"); }
    _gnukey[;j] = i;
    j++;
  }

  return 1;
};


gnulegend = function ( keys )
{
  _this_function = "gnulegend";
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
    _gnulegend[;j] = gsub( "'", "\"", i).string;
    j ++;
  }
  return 1;
};


//
// gnuwins: manage the number of open windows:
//    N - number of windows
//
//
gnuwins = function (N, s, fn)
{
  _this_function = "gnuwins";
  // check N: not reading the manual makes user
  // VERY creative
  //   global(_gnuwins);

  gw = _gnuwins();

  if (!exist(N))
  { return gw; }

  if (class(N)!="real" && class(N)!="num")
  {
    printf("gnuplot: (gnuwins) improper first argument\n");
    return <<>>;
  }
  if (N<0 || N>32)
  {
    printf("gnuplot: (gnuwins) improper first argument\n");
    return <<>>;
  }

  if (!exist(fn))
  { fn = "/dev/null"; }

  if (!exist(s))
  { s  = ""; }

  // count the windows and open/close windows if necessary
  w  = gw.win;
  nw = length(w);

  // create more windows
  if (N > nw)
  {
    for (i in (nw+1):N)
    {
      if (strlen(s[min(i,length(s))]) > 0)
      {
        __cmd = s[min(i,length(s))];
        if(any(__cmd == members(GNU_DEFAULTS.set_term)))
        { __cmd = "set term " + GNU_DEFAULTS.set_term.[ __cmd ]; }
//         "gnuwins: gnustart(,__cmd,fn)\n"?
        gnustart(,__cmd,fn);
      }
      else
      {
//         "gnuwins: gnustart(,,fn)\n"?
        gnustart(,,fn);
      }
    }
    return _gnuwins();
  } // if (N > nw)

  // close extra windows from the end
  if (N < nw)
  {
    for (i in nw:N+1:-1)
    {
      gnuclose(i);
    }
  }

  // update the count of windows
  gw = _gnuwins();
  w  = gw.win;

  // restart windows
  for (i in w)
  {
    gnuwin(i);
    if (strlen(s[min(i,length(s))]) > 0)
    {
      __cmd = s[min(i,length(s))];
      if(any(__cmd == members(GNU_DEFAULTS.set_term)))
      { __cmd = "set term " + GNU_DEFAULTS.set_term.[ __cmd ]; }

      gnucmd( __cmd, _this_function);

    }
    else
    {
      gnucmd(GNU_DEFAULTS.reset_wins, _this_function);
      gnu_object_count = 0;
      gnu_arrow_count = 0;
    }
  }
  return _gnuwins();
};

mingnuplotwins = function (n)
{
  _this_function = "mingnuplotwins";
  s = gnuwins();

  if (exist(n))
  {
    if (class(n)=="num")
    {
      n = n[1];
      if (length(s.win)<n)
      {
        gnuwins( n );
        s = gnuwins();
      }
      if (length(s.win)>n)
      {
        for (i in (n+1):length(s.win))
        {
          gnuwin (s.win[i]);
          gnucmd ("reset", _this_function);
        }
        gnuwin (s.win[1]);
      }
    }
  }

  return s;
};


gnuplot = function (data, file, fmt)
{
  _this_function = "gnuplot";
  //global (_gnukey,_gnulegend);
  global(_gnuwins);

  I = gnuwin();
  if (isempty(I))
  { gnuwins(1); }

  gw = _gnuwins();

  if (exist (file))
  {
    if (class(file) == "string")
    {
      if (!exist(fmt))
      {
        // can we figure out the format from file name
        _i1 = last(findstr(file, "."));
        if (_i1)
        {
          fmt = substr(file, (_i1+1):strlen(file));
        }
        else
        {
          fmt = "";
        }
        // check that the user is not providing something we have defined
        if (any(fmt==members(GNU_DEFAULTS.set_term)))
        {
          fmt = GNU_DEFAULTS.set_term.[fmt];
        }
        else
        {
          fmt = GNU_DEFAULTS.set_term.eps;
        }
      }
      else
      {
        // check that the user is not providing something we have defined
        if (any(fmt==members(GNU_DEFAULTS.set_term)))
        { fmt = GNU_DEFAULTS.set_term.[fmt]; }
      }

      // forward the output of the terminal to a file
      // only if talking to an open gnuplot device
      _i = gnuwin();
      _j = find(_i == gw.win);
      if (isempty(_j))
      { stop("gnuplot: terrible internal error. what the heck are you doing?"); }
      _d = gw.dev[_j];
      if (_d == "|gnuplot")
      {
        if (GNU_DEFAULTS.info)
        {
          printf ("gnuwin(%g): (gnuplot) forwarding output to file '%s'\n", I+0.0, file);
        }
        gnucmd ("set term " + fmt, _this_function);
        gnucmd ("set output '" + file + "'", _this_function);
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
  else
  {
    plotopts = "with lines";
  }

  if (class(data) == "UNDEF")
  {
    // first dataset
    plotcmd = "plot NaN";
    if (strlen(plotopts[1])>0)
    { plotcmd = plotcmd + " " + plotopts[1]; }
    if(strlen(plottitle[1])>0)
    {
      plotcmd = plotcmd + " title \"" + plottitle[1] + "\" ";
    }
    else
    {
      plotcmd = plotcmd + " notitle ";
    }

    gnucmd (plotcmd, _this_function);
  }


  if (class(data) == "num" || ishist(data))
  {
    // plot
    //   data = [x, y1, y2, ....]
    // or
    //   data = <<bin; range>>

    // handle histogram data right away.
    if (ishist(data))
    { data = hist_line(data,GNU_DEFAULTS.histmode); }

    // first dataset
    plotcmd = "plot '-'";
    if (strlen(plotopts[1])>0)
    { plotcmd = plotcmd + " " + plotopts[1]; }
    if(strlen(plottitle[1])>0)
    {
      plotcmd = plotcmd + " title \"" + plottitle[1] + "\" ";
    }
    else
    {
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
        }
        else
        {
          plotcmd = plotcmd + " with lines";
        }
      }
      else
      {
        plotcmd = plotcmd + " with lines";
      }

      // is there a title for the data set
      if (length(plottitle)>=j-1)
      {
        if (strlen(plottitle[j-1])>0)
        {
          plotcmd = plotcmd + " title \"" + plottitle[j-1] + "\" ";
        }
        else
        {
          plotcmd = plotcmd + " notitle ";
        }
      }
      else
      {
        plotcmd = plotcmd + " notitle ";
      }

    }
    gnucmd (plotcmd, _this_function);

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
    if (class(data.[members(data)[1]]) == "list" && !ishist(data.[members(data)[1]]))
    { error("gnuplot: don't know how to plot a sublist");}

    // first dataset
    if (class(data.[members(data)[1]]) == "list" && ishist(data.[members(data)[1]]))
    { plotcmd = "plot '-'";}
    if (class(data.[members(data)[1]]) == "num")
    { plotcmd = "plot '-'"; }
    if (class(data.[members(data)[1]]) == "string" && members(data)[1]!="data" )
    {
      if (isfile(data.[members(data)[1]]))
      {
        plotcmd = "plot '" + data.[members(data)[1]][1] + "'";
      }
      else
      {
        plotcmd = "plot "  + data.[members(data)[1]][1] ;
      }
    }
    else
    {
      plotcmd = "plot '-'";
    }


    if (strlen(plotopts[1])>0)
    { plotcmd = plotcmd + " " + plotopts[1]; }

    if(strlen(plottitle[1])>0)
    {
      plotcmd = plotcmd + " title \"" + plottitle[1] + "\" ";
    }
    else
    {
      plotcmd = plotcmd + " notitle ";
    }

    // the other datasets
    for (j in 2:size(data))
    {
      if (class(data.[members(data)[j]]) == "list" && ishist(data.[members(data)[j]]))
      { plotcmd = plotcmd + ",'-' ";}

      if (class(data.[members(data)[j]]) == "num")
      { plotcmd = plotcmd + ",'-' "; }

      if (class(data.[members(data)[j]]) == "string" && members(data)[j]!="data" )
      {
        if (isfile(data.[members(data)[j]]))
        {
          plotcmd = plotcmd + ",'" + data.[members(data)[j]][1] ...
              + "' ";
        }
        else
        {
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
        }
        else
        {
          plotcmd = plotcmd + "with lines";
        }
      }
      else
      {
        plotcmd = plotcmd + "with lines";
      }
      // is there a title for the data set
      if (length(plottitle)>=j)
      {
        if (strlen(plottitle[j])>0)
        {
          plotcmd = plotcmd + " title \"" + plottitle[j] + "\" ";
        }
        else
        {
          plotcmd = plotcmd + " notitle ";
        }
      }
      else
      {
        plotcmd = plotcmd + " notitle ";
      }
    }
    gnucmd (plotcmd, _this_function);

    // send data
    for (j in members(data))
    {
      if (ishist(data.[j]))
      {
        gnuprint(hist_line(data.[j],GNU_DEFAULTS.histmode));
      }

      if (class(data.[j])=="num")
      {
        gnuprint( data.[j] );
      }

      if (class(data.[j])=="string" && j=="data")
      {
        gnuprint( data.[j] );
      }
    } //for (j in members(data))
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
    }
    else
    {
      plotcmd = "plot "  + data[1] + " " ;
    }
    if (strlen(plotopts[1])>0)
    { plotcmd = plotcmd + plotopts[1]; }
    if(strlen(plottitle[1])>0)
    {
      plotcmd = plotcmd + " title \"" + plottitle[1] + "\" ";
    }
    else
    {
      plotcmd = plotcmd + " notitle ";
    }

    // the rest of the entries
    for (j in 2:length(data))
    {
      if (isfile(data[j]))
      {
        plotcmd = plotcmd + ",'" + data[j] + "' " ;
      }
      else
      {
        plotcmd = plotcmd + ", " + data[j] + " "  ;
      }

      // is there an options string for the data set
      if (length(plotopts)>=j)
      {
        if (strlen(plotopts[j])>0)
        {
          plotcmd = plotcmd + " " + plotopts[j];
        }
        else
        {
          plotcmd = plotcmd + " " + "with lines";
        }
      }
      else
      {
        plotcmd = plotcmd + " " + "with lines";
      }
      // is there a title for the entry
      if (length(plottitle)>=j)
      {
        if (strlen(plottitle[j])>0)
        {
          plotcmd = plotcmd + " title \"" + plottitle[j] + "\" ";
        }
        else
        {
          plotcmd = plotcmd + " notitle ";
        }
      }
      else
      {
        plotcmd = plotcmd + " notitle ";
      }
    }
    gnucmd (plotcmd, _this_function);
  }

  clear (_gnukey,_gnulegend);

  if (exist (file))
  {
    // reset the terminal
    _i = gnuwin();
    _j = find(_i == gw.win);
    _d = gw.dev[_j];
    if (_d == "|gnuplot" && class(file) == "string")
    {
      gnucmd (GNU_DEFAULTS.reset_term, _this_function); // linux
    }
  }

  return 0;
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
  _this_function = "gnusplot";
  //global (_gnukey, _gnulegend);
  global(_gnuwins);

  I = gnuwin();
  if (isempty(I))
  { gnuwins(1); }

  gw = _gnuwins();

  if (exist (file))
  {
    if (class(file) == "string")
    {
      if (!exist(fmt))
      {
        // can we figure out the format from file name
        _i1 = last(findstr(file, "."));
        if (_i1)
        {
          fmt = substr(file, (_i1+1):strlen(file));
        }
        else
        {
          fmt = "";
        }
        // check that the user is not providing something we have defined
        if (any(fmt==members(GNU_DEFAULTS.set_term)))
        {
          fmt = GNU_DEFAULTS.set_term.[fmt];
        }
        else
        {
          fmt = GNU_DEFAULTS.set_term.eps;
        }
      }
      else
      {
        // check that the user is not providing something we have defined
        if (any(fmt==members(GNU_DEFAULTS.set_term)))
        { fmt = GNU_DEFAULTS.set_term.[fmt]; }
      }

      // forward the output of the terminal to a file
      // only if talking to an open gnuplot device
      _i = gnuwin();
      _j = find(_i == gw.win);
      if (isempty(_j))
      { error("gnuplot: (gnuplot) terrible internal error. what the heck are you doing?"); }
      _d = gw.dev[_j];
      if (_d == "|gnuplot")
      {
        if (GNU_DEFAULTS.debug || GNU_DEFAULTS.info)
        { printf ("gnuwin(%g): (gnusplot) forwarding output to file '%s'\n", I+0.0, file); }
        gnucmd ("set term " + fmt, _this_function);
        gnucmd ("set output '"+ file + "'", _this_function);
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
    //   <<bin;xrange;yrange>>                      // list - histogram2d
    //   data.[i].x, data.[i].y, data.[i].z         // list
    //   data.[i] = 'file name'                     // string
    //   data.[i] = 'gnuplot compatible command'    // string
    //
    if (size(data) == 0)
    { return 0; }

    // check for <<x;y;z>>
    if (sum(members(data)=="x") && sum(members(data)=="y") && sum(members(data)=="z"))
    { plotcmd = "splot '-' "; }

    if (class(data) == "list" && ishist2(data))
    { plotcmd = "splot '-' "; }

    if (class(data.[members(data)[1]]) == "list" && !ishist2(data))
    { plotcmd = "splot '-' "; }

    if (class(data.[members(data)[1]]) == "string")
    {
      if (isfile(data.[members(data)[1]]))
      {
        plotcmd = "splot '" ...
            + data.[members(data)[1]][1] + "' ";
      }
      else
      {
        plotcmd = "splot " ...
            + data.[members(data)[1]][1] + " ";
      }
    }
    if (strlen(plotopts[1]) > 0)
    { plotcmd = plotcmd + " " + plotopts[1]; }
    if (strlen(plottitle[1]) > 0)
    {
      plotcmd = plotcmd + " title \"" + plottitle[1] + "\"";
    }
    else
    {
      plotcmd = plotcmd + " notitle";
    }

    for (j in 2:(size(data)*(!ishist2(data))))
    {
      // check for <<x;y;z>>, if so stop immediatelly
      if (sum(members(data)=="x") && sum(members(data)=="y") && sum(members(data)=="z"))
      { break; }

      if (class(data.[members(data)[j]]) == "list" || ishist2(data.[members(data)[j]]))
      { plotcmd = plotcmd + ",'-' "; }
      if (class(data.[members(data)[j]]) == "string")
      {
        if (isfile(data.[members(data)[j]]))
        {
          plotcmd = plotcmd + ",'" + data.[members(data)[j]][1] + "' ";
        }
        else
        {
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
        }
        else
        {
          plotcmd = plotcmd + " notitle";
        }
      }
      else
      {
        plotcmd = plotcmd + " notitle";
      }
    }
    gnucmd (plotcmd, _this_function);

    // send data
    if (class(data)=="list" && ishist2(data))
    {
      // input argument is histogram2d
      dummy = hist2_surf(data);
      { gnuprint( dummy.x, dummy.y, dummy.z ); }

    }
    else
    {
      // input argument is a list
      // <<x;y;z>>
      if (exist(data.x) && exist(data.y) && exist(data.z))
      {
        gnuprint( data.x, data.y, data.z );
        if (GNU_DEFAULTS.debug)
        { printf("gnuwin(%g): (gnusplot) print %s x,y,z\n", I+0.0, j); }
      }
      else
      {
        if (exist(data.x) || exist(data.y) || exist(data.z))
        {
          printf("gnuwin(%g): (gnusplot) <<x;y;z>> list is required\n", I+0.0);
          error ("gnuwin("+text(I+0.0,"%.0f")+"): (gnusplot) Did I say I can print?\n");
          return 1;
        }

        // input argument is a list of lists
        for (j in members(data))
        {
          // this has been passed to gnuplot earlier
          if (class(data.[j])=="string")
          { continue; }

          if (class(data.[j])!="list")
          {
            printf("gnuwin(%g): (gnusplot) Improperly formatted list!\n", I+0.0);
            error ("gnuwin("+text(I+0.0,"%.0f")+"): (gnusplot) Did I say I can print?\n");
            return 1;
          }

          // is it 2-D histrogram?
          if (ishist2(data.[j]))
          {
            dummy = hist2_surf(data.[j]);
            gnuprint( dummy.x, dummy.y, dummy.z );
            if (GNU_DEFAULTS.debug)
            { printf("gnuwin(%g): (gnusplot) print list entry %s as 2-D histogram\n", I+0.0, j); }
            continue;
          }

          // is it list of <<x;y;z>> lists?
          if (class(data.[j])=="list" && !ishist2(data.[j]))
          {
            if (exist(data.[j].x) && exist(data.[j].y) && exist(data.[j].z))
            {
              gnuprint( data.[j].x, data.[j].y, data.[j].z );
              if (GNU_DEFAULTS.debug)
              { printf("gnuwin(%g): (gnusplot) print list entry [%s] x,y,z\n", I+0.0, j); }
            }
            else
            {
              printf("gnuwin(%g): (gnusplot) Improperly formatted list!\n", I+0.0);
              error ("gnuwin("+text(I+0.0,"%.0f")+"): (gnusplot) Did I say I can print?\n");
              return 1;
            }
          }
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
    }
    else
    {
      plotcmd = "splot "  + data[1] + " ";
    }
    if (strlen(plotopts[1]) > 0)
    { plotcmd = plotcmd + " " + plotopts[1]; }
    if (strlen(plottitle[1]) > 0)
    {
      plotcmd = plotcmd + " title \"" + plottitle[1] + "\"";
    }
    else
    {
      plotcmd = plotcmd + " notitle";
    }

    for (j in 2:length(data))
    {
      if (isfile(data[j]))
      {
        plotcmd = plotcmd + ",'" + data[j] + "' ";
      }
      else
      {
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
        }
        else
        {
          plotcmd = plotcmd + " notitle";
        }
      }
      else
      {
        plotcmd = plotcmd + " notitle";
      }
    }

    gnucmd (plotcmd, _this_function);
  }

  clear (_gnukey,_gnulegend);

  if (exist (file))
  {
    // reset the terminal
    _i = gnuwin();
    _j = find(_i == gw.win);
    _d = gnuwins().dev[_j];
    if (_d == "|gnuplot" && class(file) == "string")
    {
      if (GNU_DEFAULTS.debug || GNU_DEFAULTS.info)
      {
        printf ("gnuwin(%g): (gnusplot) forwarding output to '%s' completed\n", I+0.0, file);
      }
      gnucmd (GNU_DEFAULTS.reset_term, _this_function); // linux
    }
  }

  return 1;
};



