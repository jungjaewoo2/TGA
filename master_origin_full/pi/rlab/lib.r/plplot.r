//
// New plplot.r for use with PLPLOT shared library.
// The help files for these functions are in
// misc/plhelp
//

// plplot.r
require libgeneral

// This file is a part of RLaB ("Our"-LaB)
// Version 2, Copyright (C) 2014-2016  Marijan Kostrun
// Version 1, Copyright (C) 1994  Ian R. Searle

// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.

// See the file ../COPYING

//
// If your system does not deal with Infs and NaNs well, then
// uncomment the following lines.
//
require libgeneral libstdio libsystem

// static (isinf, isnan)
//
// isinf = function ( A ) { return (0); };
// isnan = function ( A ) { return (0); };

//
// Book-keeping variables
//
static (PLPLOT_PLOT_WINDOWS); // The static plot window structure
if (!exist(PLPLOT_PLOT_WINDOWS))
{ PLPLOT_PLOT_WINDOWS = <<>>; }
static (PLPLOT_ACTIVE_WIN);   // The active/current plot window
if (!exist(PLPLOT_ACTIVE_WIN))
{ PLPLOT_ACTIVE_WIN = 0; }

//
// PLPLOT CONSTANTS
//
static (PLPLOT_COLORS, PLPLOT_LINE_STYLE, PLPLOT_POINT_STYLE, PLPLOT_FILL_STYLE);
PLPLOT_COLORS = 0:15;
PLPLOT_LINE_STYLE = 1:8;
PLPLOT_POINT_STYLE = 1:8;
PLPLOT_FILL_STYLE = 1:8;
static (PLPLOT_RESET_AFTER_PLPLOT);
PLPLOT_RESET_AFTER_PLPLOT = 1;
static(PLPLOT_DATASET_PLOT_FORMAT);
PLPLOT_DATASET_PLOT_FORMAT = [0, 1, 1, 1, 0, 1, 0, 0, 0, 0];
// Maintain the transformations for 3-D plots.
static (PLPLOT_BASE_X, PLPLOT_BASE_Y, PLPLOT_BASE_Z);
PLPLOT_BASE_X = 2;
PLPLOT_BASE_Y = 2;
PLPLOT_BASE_Z = 4;
// Defaults
//  2D grids
static (PLPLOT_GRID_X_DEFAULT, PLPLOT_ALT_GRID_X_DEFAULT);
PLPLOT_GRID_X_DEFAULT = "bcnst";
PLPLOT_ALT_GRID_X_DEFAULT = "cst";
static(PLPLOT_GRID_X_DEFAULT_LOG,PLPLOT_ALT_GRID_X_DEFAULT_LOG);
PLPLOT_GRID_X_DEFAULT_LOG = "bclnst";
PLPLOT_ALT_GRID_X_DEFAULT_LOG = "lcst";
static (PLPLOT_GRID_Y_DEFAULT,PLPLOT_ALT_GRID_Y_DEFAULT);
PLPLOT_GRID_Y_DEFAULT = "bnstv";
PLPLOT_ALT_GRID_Y_DEFAULT = "cmstv";
static(PLPLOT_GRID_Y_DEFAULT_LOG,PLPLOT_ALT_GRID_Y_DEFAULT_LOG);
PLPLOT_GRID_Y_DEFAULT_LOG = "lbcnstv";
PLPLOT_ALT_GRID_Y_DEFAULT_LOG = "lcstv";
// 3D grids
static (PLPLOT_GRID3D_X_DEFAULT, PLPLOT_GRID3D_Y_DEFAULT, PLPLOT_GRID3D_Z_DEFAULT);
PLPLOT_GRID3D_X_DEFAULT = "bnstu";
PLPLOT_GRID3D_Y_DEFAULT = "bnstu";
PLPLOT_GRID3D_Z_DEFAULT = "bcdmnstuv";
// PLPLOT_GRID3D_Z_DEFAULT = "bnstu";
// Initial windows size/config...
static (PLPLOT_WIN_LEN_X_DEFAULT, PLPLOT_WIN_LEN_Y_DEFAULT);
PLPLOT_WIN_LEN_X_DEFAULT = 560;
PLPLOT_WIN_LEN_Y_DEFAULT = 400;
static (PLPLOT_WIN_OFFSET_X_DEFAULT, PLPLOT_WIN_OFFSET_Y_DEFAULT);
PLPLOT_WIN_OFFSET_X_DEFAULT = 100;
PLPLOT_WIN_OFFSET_Y_DEFAULT = 100;
static (PLPLOT_WIN_XP_DEFAULT, PLPLOT_WIN_YP_DEFAULT);
PLPLOT_WIN_XP_DEFAULT = 0;
PLPLOT_WIN_YP_DEFAULT = 0;

//
// default plplot color0 scheme
//
static(PLPLOT_COLOR_SCHEME, PLPLOT_COLOR_SCHEME_DEFAULT, PLPLOT_COLOR_RGB_DEFAULT, PLPLOT_COLOR_RGB_DEFAULT_NAMES);
static(PLPLOT_COLOR_COUNT_DEFAULT);

//
// default color scheme
//
PLPLOT_COLOR_COUNT_DEFAULT = 15;
PLPLOT_COLOR_SCHEME_DEFAULT = [0:PLPLOT_COLOR_COUNT_DEFAULT]';
PLPLOT_COLOR_RGB_DEFAULT = _plgcol0([0:PLPLOT_COLOR_COUNT_DEFAULT]');
if (!exist(PLPLOT_COLOR_RGB_DEFAULT_NAMES))
{
  PLPLOT_COLOR_RGB_DEFAULT_NAMES = [ ...
      "black", "red", "yellow", "green", "aquamarine", "pink", "wheat", "grey", ...
      "brown", "blue", "blueviolet", "cyan", "turquoise", "magenta", "salmon", "white"];
}
static(PLPLOT_COLOR_RGB_NAMES);
if (!exist(PLPLOT_COLOR_RGB_NAMES))
{ PLPLOT_COLOR_RGB_NAMES = PLPLOT_COLOR_RGB_DEFAULT_NAMES; }

//
// color scheme after xmgrace, and gnuplot
//
static(PLPLOT_COLOR_COUNT_XMG, PLPLOT_COLOR_RGB_XMG_NAMES, PLPLOT_COLOR_RGB_XMG);
PLPLOT_COLOR_COUNT_XMG = 16;
PLPLOT_COLOR_RGB_XMG_NAMES = [ ...
    "white", "black", "red", "green", "blue", "yellow", ...
        "brown", "grey", "violet", "cyan", "magenta", "orange", ...
            "indigo", "maroon", "turquoise", "green5" ];
PLPLOT_COLOR_RGB_XMG = [ ...
    [255, 255, 255]; [0, 0, 0]; [255, 0, 0]; [0, 255, 0]; [0, 0, 255]; [255, 255, 0];...
    [188, 143, 143]; [153, 153, 153]; [148, 0, 211]; [0, 255, 255]; [255, 0, 255]; [255, 165, 0]; ...
    [114, 33, 188]; [103, 7, 72]; [64, 224, 208]; [0, 59, 0] ];

//
// decide on default color scheme
//
static(PLPLOT_COLOR_SCHEME,PLPLOT_COLOR_RGB);
PLPLOT_COLOR_SCHEME = "wb"; //plplot's white on black 
// PLPLOT_COLOR_SCHEME = "bw"; //plplot's black on white
// PLPLOT_COLOR_SCHEME = "xmg"; // xmgrace/gnuplot's black on white

//
// plot management
//
static(PLPLOT_AXES_ID, PLPLOT_AXES_DESC, PLPLOT_AXES_REV);
PLPLOT_AXES_ID    = [     1,      2,      3,      4];
PLPLOT_AXES_DESC  = ["x1y1", "x1y2", "x2y1", "x2y2"];
PLPLOT_AXES_REV   = <<x1y1=1; x1y2=2; x2y1=3; x2y2=4>>;

static(PLPLOT_DEFAULT_LINE_STYLE);
PLPLOT_DEFAULT_LINE_STYLE = 1;

//
// currently defined for
//  bw = black et c. on white background
//
static(modify_color0_scheme);
modify_color0_scheme = function(color_scheme)
{
  if (color_scheme == "bw")
  {
    PLPLOT_COLOR_RGB_NAMES = blank(PLPLOT_COLOR_RGB_DEFAULT_NAMES);

    PLPLOT_COLOR_RGB[1;] = PLPLOT_COLOR_RGB_DEFAULT[16;]; // background 0 <- white
    PLPLOT_COLOR_RGB_NAMES[1] = PLPLOT_COLOR_RGB_DEFAULT_NAMES[16];

    // shift all by one: update names
    for (i in 1:15)
    {
      PLPLOT_COLOR_RGB[i+1;] = PLPLOT_COLOR_RGB_DEFAULT[i;];
      PLPLOT_COLOR_RGB_NAMES[i+1] = PLPLOT_COLOR_RGB_DEFAULT_NAMES[i];
    }

    _plscol0  (PLPLOT_COLOR_SCHEME_DEFAULT, PLPLOT_COLOR_RGB);
    _plscolbg (PLPLOT_COLOR_RGB[1;]); // 1 is background color
    return 0;
  }

  if (color_scheme == "wb")
  {
    PLPLOT_COLOR_RGB_NAMES = blank(PLPLOT_COLOR_RGB_DEFAULT_NAMES);

    PLPLOT_COLOR_RGB = PLPLOT_COLOR_RGB_DEFAULT; // background 0 <- white
    PLPLOT_COLOR_RGB_NAMES = PLPLOT_COLOR_RGB_DEFAULT_NAMES;

    _plscol0  (PLPLOT_COLOR_SCHEME_DEFAULT, PLPLOT_COLOR_RGB);
    _plscolbg (PLPLOT_COLOR_RGB[1;]); // 1 is background color
    return 0;
  }

  if (color_scheme == "xmg")
  {
    PLPLOT_COLOR_RGB_NAMES = PLPLOT_COLOR_RGB_XMG_NAMES;
    PLPLOT_COLOR_RGB = PLPLOT_COLOR_RGB_XMG;

    _plscol0  (PLPLOT_COLOR_SCHEME_DEFAULT, PLPLOT_COLOR_RGB);
    _plscolbg (PLPLOT_COLOR_RGB[1;]); // 1 is background color
    return 0;
  }

  return 1;
};

//
// load X11 palette from the configuration file
//
static(PLPLOT_X11_PALETTE_FILE, PLPLOT_X11_PALETTE);
PLPLOT_X11_PALETTE_FILE = "/usr/share/X11/rgb.txt";
PLPLOT_X11_PALETTE = <<>>;
if (isfile(PLPLOT_X11_PALETTE_FILE))
{
  __s = reads(PLPLOT_X11_PALETTE_FILE);
  for(__i in length(__s):1:-1)
  {
    if(strlen(__s[__i])>0)
    { break; }
  }
  if (__i > 1)
  {
    // format of X11's rgb.txt file is
    // R G B rgbColorName
    for (__j in 1:__i)
    {
      __y = strsplt(__s[__j],"'BLANK"); // since 2015-5-20/2.4.2.2
      __l = length(__y);
      if(__l>=4)
      {
        if (__l > 4)
        {
          __name = tolower(join(__y[4:__l],""));
        else
          __name = __y[4];
        }
        if (!exist(PLPLOT_X11_PALETTE.[__name]))
        { PLPLOT_X11_PALETTE.[__name] = strtod(__y[1:3]); }
      }
    }
  }
  clear(__s, __i, __j, __l);
}

static (check_plot_object);
static (find_scales);
static (xy_scales);
static (x_scales);
static (y_scales);
static (z_scales);
static (XYZ_scales);
static (list_scales);
static (list_sort);
static (plot_matrix);
static (plot_list);
static (check_3d_list);
static (get_style);
static (plhold_first);


static (subplot_f)
subplot_f = 0;

//
// Static (private) functions. For use from within
// this file only.
//
static(convert_gnuformat_to_plformat);
convert_gnuformat_to_plformat = function(s)
{
  THIS_SOLVER = "convert_gnuformat_to_plformat";

  if(class(s)!="string")
  { return PLPLOT_DATASET_PLOT_FORMAT; }

  tokens = tolower(strsplt(s, " "));

  //
  // lines or points: use gnuplot
  //
  _iw = find(tokens == "w");
  if (isempty(_iw))
  { _iw = find(tokens == "with"); }

  lpb = 0;
  err = 0;
  err_xy = 0;
  err_lp = 0;
  if (length(_iw)==1 && _iw < length(tokens))
  {
    if (strindex(tolower(tokens[_iw+1]), "p"))
    {
      lpb = 1;
    else if (strindex(tolower(tokens[_iw+1]), "b"))
    {
      lpb = 2;
    }}

    if (strindex(tolower(tokens[_iw+1]), "err"))
    {
      err = 1;
      err_xy = 3;
      err_lp = 0;
      if (strindex(tolower(tokens[_iw+1]), "xy"))
      {
        err_xy=2;
      else if (strindex(tolower(tokens[_iw+1]), "x"))
      {
        err_xy=1;
      }}
      if (strindex(tolower(tokens[_iw+1]), "bar"))
      {
        lpb = 1;
        err_lp = 1;
      }
    }
  }

  //
  // line has three configurable elements:
  //    linetype  1:8
  //    linewidth device dependent
  //    linecolor 0:15
  lt = 0;
  lc = 0;
  lw = 0;
  if ((lpb==0) || (lpb==2))
  {
    //
    _ilt = find(tokens == "lt");
    if ((length(_ilt)==1) && _ilt < length(tokens))
    {
      lt = strtod(tokens[_ilt+1]);
      if (lt < 1 || lt>8)
      { lt = -1; }
    else
      lt = -1;
    }
    //
    _ilc = find(tokens == "lc");
    if ((length(_ilc)==1) && _ilc < length(tokens))
    {
      if (tokens[_ilc+1] == "rgb")
      {
        c = gsub("'",tokens[_ilc+2]).string;
        _idx_c = find(PLPLOT_COLOR_RGB_NAMES==c);
        if (length(_idx_c)==1)
        {
          lc = _idx_c - 1;
        else
          printf("PLPLOT: Warning! Requested color '%s' not available in color0 palette!\n");
          lc = -1; // unknown color
        }
      else
        lc = strtod(tokens[_ilc+1]);
      }
      while(lc<1)
      { lc = lc + 15; }
      while (lc > 15)
      { lc = lc - 15; }
      lc - mod(lc,15);
    else
      lc = -1;
    }

    //
    lw = PLPLOT_DATASET_PLOT_FORMAT[4];
    _ilw = find(tokens == "lw");
    if ((length(_ilw)==1) && _ilw < length(tokens))
    { lw = strtod(tokens[_ilw+1]); }

    if ((lc==-1)&&(lt!=-1))
    { lc = PLPLOT_DATASET_PLOT_FORMAT[3]; }
    if ((lt==-1)&&(lc!=-1))
    { lt = PLPLOT_DATASET_PLOT_FORMAT[2]; }
  }

  //
  // point has four configurable elements:
  //    pointtype  1:8
  //    linewidth device dependent
  //    pointcolor 1:14
  //    pointsize  any
  pt = 0;
  pc = 0;
  ps = 0;
  if ((lpb==1) || (lpb==2))
  {
    //
    _ipt = find(tokens == "pt");
    if ((length(_ipt)==1) && _ipt < length(tokens))
    {
      pt = strtod(tokens[_ipt+1]);
    else
      pt = -1;
    }
    //
    _ipc = find(tokens == "pc");
    if (isempty(_ipc))
    {
      _ipc = find(tokens == "lc");
    }
    if ((length(_ipc)==1) && _ipc < length(tokens))
    {
      if (tokens[_ipc+1] == "rgb")
      {
        c = gsub("'",tokens[_ipc+2]).string;
        _idx_c = find(PLPLOT_COLOR_RGB_NAMES==c);
        if (length(_idx_c)==1)
        {
          pc = _idx_c - 1;
        else
          printf("PLPLOT: Warning! Requested color '%s' is not available in palette color0!\n", c);
          printf("PLPLOT: Warning! Available colors are %s. Please choose another!\n", ...
            join(PLPLOT_COLOR_RGB_NAMES,","));
          pc = 1; // unknown color
        }
      else
        pc = strtod(tokens[_ipc+1]);
      }
      if (pc < 1 || pc>14)
      { pc = -1; }
    else
      pc = -1;
    }
    //
    _ips = find(tokens == "ps");
    if ((length(_ips)==1) && _ips < length(tokens))
    {
      ps = strtod(tokens[_ips+1]);
    else
      ps = PLPLOT_DATASET_PLOT_FORMAT[4];
    }

    if ((pc==-1)&&(pt!=-1))
    { pc = PLPLOT_DATASET_PLOT_FORMAT[6]; }
    if ((pt==-1)&&(pc!=-1))
    { pt = PLPLOT_DATASET_PLOT_FORMAT[5]; }
  }

  //
  // is user using key word 'using'?
  //
  _iu = find(tokens == "u");
  if (isempty(_iu))
  { _iu = find(tokens == "using"); }
  if ((length(_iu)==1) && max(_iu) < length(tokens))
  {
    // nobody reads manual anyway
    if ((strindex(tokens[_iu+1],"(")>0) || (strindex(tokens[_iu+1],")")>0) ...
         || (strindex(tokens[_iu+1],"*")>0) || (strindex(tokens[_iu+1],"/")>0) ...
         || (strindex(tokens[_iu+1],"+")>0) || (strindex(tokens[_iu+1],"-")>0))
    {
      printf("PLPLOT: " + THIS_SOLVER + ": Cannot convert algabraic expression in %s that involves columns!\n", ...
          tokens[_iu+1]);
      printf("PLPLOT: " + THIS_SOLVER + ": Why don't you read the f*!king manual, while I rest for a while!\n", ...
          tokens[_iu+1]);
      rerror("PLPLOT: " + THIS_SOLVER + "\n");
    }
    use_cols = strtod(tokens[_iu+1],<<csp=":">>);
    if (err)
    {
      if ((err_xy==1)&&(length(use_cols)!=4))
      {
        printf("PLPLOT: " + THIS_SOLVER + ": 'with xerror{lines,bars}' requires 4 data columns (x,y,xmin,xmax).\n" );
        printf("PLPLOT: " + THIS_SOLVER + ": 'using %s' provided only %g columns. Cannot continue!\n", ...
            tokens[_iu+1],  length(use_cols));
        rerror("PLPLOT: " + THIS_SOLVER + "\n");
      }
      if ((err_xy==3)&&(length(use_cols)!=4))
      {
        printf("PLPLOT: " + THIS_SOLVER + ": 'with [y]error{lines,bars}' requires 4 data columns (x,y,ymin,ymax).\n" );
        printf("PLPLOT: " + THIS_SOLVER + ": 'using %s' provided only %g columns. Cannot continue!\n", ...
            tokens[_iu+1],  length(use_cols));
        rerror("PLPLOT: " + THIS_SOLVER + "\n");
      }
      if ((err_xy==2)&&(length(use_cols)!=6))
      {
        printf("PLPLOT: " + THIS_SOLVER + ": 'with xyerror{lines,bars}' requires 6 data columns (x,y,xmin,xmax,ymin,ymax).\n" );
        printf("PLPLOT: " + THIS_SOLVER + ": 'using %s' provided only %g columns. Cannot continue!\n", ...
            tokens[_iu+1],  length(use_cols));
        rerror("PLPLOT: " + THIS_SOLVER + "\n");
      }
    }
  else
    if (err)
    {
      if ((err_xy==1) || (err_xy==3))
      {
        use_cols = [1:4]; // default configuration
      else
        use_cols = [1:6]; // default configuration
      }
    else
      use_cols = [1,2]; // default configuration
    }
  }

  //
  // is user using key word 'axes'?
  //
  axes = PLPLOT_AXES_REV.x1y1;  // default. in case user does not use it!
  _ia = find(tokens == "ax");
  if (isempty(_ia))
  { _ia = find(tokens == "axes"); }
  if (isempty(_ia))
  {
    _ia = find(tokens == "axis");
    if (!isempty(_ia))
    {
      printf("PLPLOT: " + THIS_SOLVER + ": If you haven't know that, plural of 'axis' is 'axes'.\n");
      printf("PLPLOT: " + THIS_SOLVER + ": I would send you away to practice grammar, but my boss told me to behave!\n");
    }
  }
  if ((length(_ia)==1) && max(_ia) < length(tokens))
  {
    // axis enforces use_cols=[1,2] in the absence of 'using 1:2...' in plformat
    if (isempty(use_cols))
    { use_cols = [1,2]; }
    if (exist(PLPLOT_AXES_REV.[tokens[_ia+1]]))
    {
      axes = PLPLOT_AXES_REV.[tokens[_ia+1]];
    else
      printf("PLPLOT: " + THIS_SOLVER + ": directive 'axes %s' is not supported!\n",  tokens[_ia+1]);
      printf("PLPLOT: " + THIS_SOLVER + ": Why don't you go away and read the manual for a change!\n");
      error ("PLPLOT: " + THIS_SOLVER + ": Cannot continue!\n");
    }
  }

  //
  // for 3D plots we use err, err_xy and err_lp flags
  // key command is
  //    draw ...
  _id = find(tokens == "draw");
  if (isempty(_id))
  { _id = find(tokens == "dr"); }
  if ((length(_id)==1) && max(_id) < length(tokens))
  {
    err_xy = 0;
    _nd = length(tokens);
    // identify messages:
    //   line{x}{y}
    if ( any(strindex(tokens[(_id+1):_nd], "xy")>0) )
    {
      err_xy = err_xy - 4;
    else if ( any(strindex(tokens[(_id+1):_nd], "x")>0) )
    {
      err_xy = err_xy - 1;
    else if ( any(strindex(tokens[(_id+1):_nd], "y")>0) )
    {
      err_xy = err_xy - 2;
    }}}
    // identify messages:
    //   {mag}
    if ( any(strindex(tokens[(_id+1):_nd], "mag")>0) )
    { err_xy = err_xy - 8; }
    // identify messages:
    //   {cont}
    if ( any(strindex(tokens[(_id+1):_nd], "con")>0) || any(strindex(tokens[(_id+1):_nd], "bas")>0))
    { err_xy = err_xy - 16; }
    // identify messages:
    //   {sid}
    if ( any(strindex(tokens[(_id+1):_nd], "sid")>0) )
    { err_xy = err_xy - 64; }
  }

  // finish gnuplot-style processing
  rval = <<>>;
  rval.axes = axes;
  rval.use_cols = use_cols;
  rval.plformat = [lpb, lt, lc, lw, pt, pc, ps, err, err_xy, err_lp];
  return rval;
};

static (reset_active_plot_object);
reset_active_plot_object = function ( nx, ny, xleng, yleng, xoff, yoff )
{
  if (length(PLPLOT_PLOT_WINDOWS) == 0)
  { return 1;}
  if (PLPLOT_ACTIVE_WIN == 0)
  { return 1;}

  pobj = PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN];

  if (!exist (nx))
  { nx = 1; }
  if (!exist (ny))
  { ny = 1; }

  if (!exist (xleng))
  {
    if (exist(pobj.xleng))
    {
      xleng = pobj.xleng;
    else
      xleng = PLPLOT_WIN_LEN_X_DEFAULT;
    }
  }
  if (!exist (yleng))
  {
    if (exist(pobj.yleng))
    {
      yleng = pobj.yleng;
    else
      yleng = PLPLOT_WIN_LEN_Y_DEFAULT;
    }
  }
  if (!exist (xoff))
  {
    if (exist(pobj.xoff))
    {
      xoff = pobj.xoff;
    else
      xoff = PLPLOT_WIN_OFFSET_X_DEFAULT;
    }
  }
  if (!exist (yoff))
  {
    if (exist(pobj.yoff))
    {
      yoff = pobj.yoff;
    else
      yoff = PLPLOT_WIN_OFFSET_Y_DEFAULT;
    }
  }

  pobj.plformat = <<>>;
  pobj.pltext = <<>>;
  pobj.axes = <<>>;
  pobj.cols = <<>>;
  pobj.style = <<>>;
  pobj.xtick = <<>>;
  pobj.ytick = <<>>;
  pobj.alt_xlabel_scale = <<>>;
  pobj.alt_ylabel_scale = <<>>;
  pobj.alt_zlabel_scale = <<>>;

  pobj.window.xp    = PLPLOT_WIN_XP_DEFAULT;      // Number of X pixels
  pobj.window.yp    = PLPLOT_WIN_YP_DEFAULT;      // Number of Y pixels
  pobj.window.xleng = xleng;        // Page length, X
  pobj.window.yleng = yleng;        // Page length, Y
  pobj.window.xoff  = xoff;           // Page offset, X
  pobj.window.yoff  = yoff;           // Page offset, Y

  pobj.subplot  = 0;      // The current subplot no.
  pobj.nplot    = nx*ny;  // Total no. of plots on window
  pobj.nx       = nx;
  pobj.ny       = ny;
  pobj.dev   = "xwin";
  pobj.fontld = 0;        // Loaded extended fonts?

  for (i in 1:(nx*ny))
  {
    // rlabplus:
    pobj.plformat.[i] = [];
    pobj.axes.[i] = [];
    pobj.xtick.[i] = [0,0];
    pobj.alt_xtick.[i] = [0,0];
    pobj.ytick.[i] = [0,0];
    pobj.alt_ytick.[i] = [0,0];
    pobj.style.[i] = "line";         // The type/style of plot to draw
    pobj.width[i] = 1;  // The pen width for current plot
    pobj.font[i] = 1;   // The current font
    pobj.xlabel[i] = "";
    pobj.xlabel_scale.[i] = [1.0, 3.0, 0.5, 0.5];
    pobj.alt_xlabel[i] = "";
    pobj.alt_xlabel_scale.[i] = [1.0, 2.0, 0.5, 0.5];
    pobj.ylabel[i] = "";
    pobj.ylabel_scale.[i] = [1.0, 5.0, 0.5, 0.5];
    pobj.alt_ylabel[i] = "";
    pobj.alt_ylabel_scale.[i] = [1.0, 5.0, 0.5, 0.5];
    pobj.zlabel[i] = "";
    pobj.alt_zlabel[i] = "";
    pobj.title[i] = "";
    pobj.title_scale.[i] = _plchr();
    pobj.orientation[i] = "portrait";
    pobj.pltext.[i] = cell();         // texts that go on top of the graph
    pobj.desc.[i] = blank();        // The legend description
    pobj.desc_pos.[i] = "itr";      // The legend default position is inside/top/right
    pobj.desc_pos_xy.[i] = [0,0];   // The legend default position offset from the position above
    pobj.desc_scale.[i] = 1.0;      // The legend default size and scaling
    pobj.gridx[i] =  PLPLOT_GRID_X_DEFAULT;	// Plot axes style, 2D-X
    pobj.alt_gridx[i] =  PLPLOT_ALT_GRID_X_DEFAULT; // Plot axes style, 2D-X
    pobj.gridy[i] =  PLPLOT_GRID_Y_DEFAULT;	// Plot axes style, 2D-Y
    pobj.alt_gridy[i] =  PLPLOT_ALT_GRID_Y_DEFAULT; // Plot axes style, 2D-Y
    pobj.grid3x[i] = PLPLOT_GRID3D_X_DEFAULT;	// Plot axes style, 3D-X
    pobj.grid3y[i] = PLPLOT_GRID3D_Y_DEFAULT;	// Plot axes style, 3D-Y
    pobj.grid3z[i] = PLPLOT_GRID3D_Z_DEFAULT;	// Plot axes style, 3D-Z
    pobj.aspect[i] = 0;		        // Plot aspect style
    pobj.alt[i] = 60;
    pobj.az[i] = 45;
    pobj.plvpor[i;] = nan(1,4);

    pobj.xmin[i] = 1j;
    pobj.xmax[i] = 1j;
    pobj.ymin[i] = 1j;
    pobj.ymax[i] = 1j;
    pobj.zmin[i] = 1j;
    pobj.zmax[i] = 1j;

    pobj.alt_xmin[i] = 1j;
    pobj.alt_xmax[i] = 1j;
    pobj.alt_ymin[i] = 1j;
    pobj.alt_ymax[i] = 1j;
    pobj.alt_zmin[i] = 1j;
    pobj.alt_zmax[i] = 1j;

    // load default values
    pobj.page.xp = PLPLOT_WIN_XP_DEFAULT;
    pobj.page.yp = PLPLOT_WIN_YP_DEFAULT;
    pobj.xleng   = xleng/nx;
    pobj.yleng   = yleng/ny;
    pobj.xoff    = xoff;
    pobj.yoff    = yoff;

    pobj.color[i;]  = 1:14;              // 14 possible colors...
    pobj.lstyle[i;] = 1:8;               // 8 possible line styles...
    pobj.pstyle[i;] = 1:8;               // 8 possible point styles...
    pobj.fstyle[i;] = 1:8;               // 8 possible fill patterns...
  }

  //
  // Save the newly generated plot-object
  // in a list of plot-objects.
  PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN] = pobj;
};

static (create_plot_object);
create_plot_object = function ( N, nx, ny, xleng, yleng, xoff, yoff )
{
  if (!exist (N))
  { N = 1; }
  if (!exist (nx))
  { nx = 1; }
  if (!exist (ny))
  { ny = 1; }
  if (!exist (xleng))
  { xleng = PLPLOT_WIN_LEN_X_DEFAULT; }
  if (!exist (yleng))
  { yleng = PLPLOT_WIN_LEN_Y_DEFAULT; }
  if (!exist (xoff))
  { xoff = PLPLOT_WIN_OFFSET_X_DEFAULT; }
  if (!exist (yoff))
  { yoff = PLPLOT_WIN_OFFSET_Y_DEFAULT; }

  pobj = <<>>;
  pobj.plformat = <<>>;
  pobj.pltext = <<>>;
  pobj.axes = <<>>;
  pobj.cols = <<>>;
  pobj.style = <<>>;
  pobj.xtick = <<>>;
  pobj.xtick_scale = <<>>;
  pobj.alt_xtick = <<>>;
  pobj.alt_xtick_scale = <<>>;
  pobj.ytick = <<>>;
  pobj.ytick_scale = <<>>;
  pobj.alt_ytick = <<>>;
  pobj.alt_ytick_scale = <<>>;
  pobj.alt_xlabel_scale = <<>>;
  pobj.alt_ylabel_scale = <<>>;
  pobj.alt_zlabel_scale = <<>>;
  pobj.ztick = <<>>;
  pobj.ztick_scale = <<>>;

  pobj.window.xp    = PLPLOT_WIN_XP_DEFAULT;      // Number of X pixels
  pobj.window.yp    = PLPLOT_WIN_YP_DEFAULT;      // Number of Y pixels
  pobj.window.xleng = xleng;        // Page length, X
  pobj.window.yleng = yleng;        // Page length, Y
  pobj.window.xoff  = xoff;           // Page offset, X
  pobj.window.yoff  = yoff;           // Page offset, Y

  pobj.subplot  = 0;      // The current subplot no.
  pobj.nplot    = nx*ny;  // Total no. of plots on window
  pobj.nx       = nx;
  pobj.ny       = ny;
  pobj.dev   = "xwin";
  pobj.fontld = 0;        // Loaded extended fonts?

  for (i in 1:(nx*ny))
  {
    // rlabplus:
    pobj.plformat.[i] = [];
    pobj.pltext.[i] = cell();         // texts that go on top of the graph
    pobj.axes.[i] = [];
    pobj.xtick.[i] = [0,0];
    pobj.xtick_scale.[i] = 1;
    pobj.alt_xtick.[i] = [0,0];
    pobj.alt_xtick_scale.[i] = 1;
    pobj.ytick.[i] = [0,0];
    pobj.ytick_scale.[i] = 1;
    pobj.alt_ytick.[i] = [0,0];
    pobj.alt_ytick_scale.[i] = 1;
    pobj.ztick.[i] = [0,0];
    pobj.ztick_scale.[i] = 1;
    pobj.style.[i] = "line";         // The type/style of plot to draw
    pobj.width[i] = 1;  // The pen width for current plot
    pobj.font[i] = 1;   // The current font
    pobj.xlabel[i] = "";
    pobj.xlabel_scale.[i] = [1.0, 3.0, 0.5, 0.5]; // size, display position, string pos, justification
    pobj.alt_xlabel[i] = "";
    pobj.alt_xlabel_scale.[i] = [1.0,2.0, 0.5, 0.5];
    pobj.ylabel[i] = "";
    pobj.ylabel_scale.[i] = [1.0, 5.0, 0.5, 0.5];
    pobj.alt_ylabel[i] = "";
    pobj.alt_ylabel_scale.[i] = [1.0, 5.0, 0.5, 0.5];
    pobj.zlabel[i] = "";
    pobj.alt_zlabel[i] = "";
    pobj.title[i] = "";
    pobj.title_scale.[i] = _plchr();
    pobj.orientation[i] = "portrait";
    pobj.desc.[i] = blank();        // The legend description
    pobj.desc_pos.[i] = "itr";      // The legend default position is inside/top/right
    pobj.desc_pos_xy.[i] = [0,0];   // The legend default position offset from the position above
    pobj.desc_scale.[i] = 1.0;      // The legend default size and scaling
    pobj.gridx[i] =  PLPLOT_GRID_X_DEFAULT; // Plot axes style, 2D-X
    pobj.alt_gridx[i] =  PLPLOT_ALT_GRID_X_DEFAULT; // Plot axes style, 2D-X
    pobj.gridy[i] =  PLPLOT_GRID_Y_DEFAULT; // Plot axes style, 2D-Y
    pobj.alt_gridy[i] =  PLPLOT_ALT_GRID_Y_DEFAULT; // Plot axes style, 2D-Y
    pobj.grid3x[i] = PLPLOT_GRID3D_X_DEFAULT; // Plot axes style, 3D-X
    pobj.grid3y[i] = PLPLOT_GRID3D_Y_DEFAULT; // Plot axes style, 3D-Y
    pobj.grid3z[i] = PLPLOT_GRID3D_Z_DEFAULT; // Plot axes style, 3D-Z
    pobj.aspect[i] = 0;           // Plot aspect style
    pobj.alt[i] = 60;
    pobj.az[i] = 45;
    pobj.plvpor[i;] = nan(1,4);

    pobj.xmin[i] = 1j;
    pobj.xmax[i] = 1j;
    pobj.ymin[i] = 1j;
    pobj.ymax[i] = 1j;
    pobj.zmin[i] = 1j;
    pobj.zmax[i] = 1j;

    pobj.alt_xmin[i] = 1j;
    pobj.alt_xmax[i] = 1j;
    pobj.alt_ymin[i] = 1j;
    pobj.alt_ymax[i] = 1j;
    pobj.alt_zmin[i] = 1j;
    pobj.alt_zmax[i] = 1j;

    // load default values
    pobj.page.xp = PLPLOT_WIN_XP_DEFAULT;
    pobj.page.yp = PLPLOT_WIN_YP_DEFAULT;
    pobj.xleng   = xleng/nx;
    pobj.yleng   = yleng/ny;
    pobj.xoff    = xoff;
    pobj.yoff    = yoff;

    pobj.color[i;]  = 1:14;              // 14 possible colors...
    pobj.lstyle[i;] = 1:8;               // 8 possible line styles...
    pobj.pstyle[i;] = 1:8;               // 8 possible point styles...
    pobj.fstyle[i;] = 1:8;               // 8 possible fill patterns...
  }

  //
  // Save the newly generated plot-object
  // in a list of plot-objects.
  PLPLOT_PLOT_WINDOWS.[N] = pobj;
};


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Check to make sure a plot-object exists. If one
// does not exist, create it.
//

check_plot_object = function ()
{
  if (PLPLOT_ACTIVE_WIN== 0)
  {
    plwins(1);
    plwin (1);
  }

  return 1;
};


//
// Set the current plot window
// Default value = 1
plwin = function ( N )
{
  this_function = "plwin";
  if (!exist (N))
  { return PLPLOT_ACTIVE_WIN; }


  // Check to make sure N is valid
  if (exist(PLPLOT_PLOT_WINDOWS.[N]))
  {
    PLPLOT_ACTIVE_WIN = N;
    _plsstrm (PLPLOT_ACTIVE_WIN);
    _pladv (0);

    // clean-up the settings
    if(PLPLOT_RESET_AFTER_PLPLOT)
    { reset_active_plot_object(); }

    return PLPLOT_ACTIVE_WIN;
  }

  return 0;
};

//
//rlabplus extension
//
newplwin = function ()
{

  if (length(PLPLOT_PLOT_WINDOWS) == 0)
  { return 1;}

  ps = plwins().win;
  newp = 0;
  j = 1;
  while(newp==0)
  {
    if(sum(j==ps)>0)
    {
      j++;
      continue;
    }
    newp = j;
  }

  return newp;
};

//
// rlabplus extension
//
plwins = function( NWIN, dev, sz, offs )
{
  //
  // process NWIN: an integer, number of windows requested,
  // or a matrix [nx(i);ny(i)]_i
  //
  if (!exist(NWIN))
  {
    // report what is available
    if (length(PLPLOT_PLOT_WINDOWS) == 0)
    { return <<win=[];act=[];dev=blank(0,0)>>; }

    retl = <<>>;
    retl.act   = PLPLOT_ACTIVE_WIN;
    retl.win = [];
    retl.dev    = blank(0,0);
    for(i in members(PLPLOT_PLOT_WINDOWS))
    {
      retl.win = [retl.win, strtod(i)];
      retl.dev    = [retl.dev, PLPLOT_PLOT_WINDOWS.[i].dev];
    }
    idx = sort(retl.win).idx;
    retl.dev    = retl.dev[idx];
    retl.win = retl.win[idx];
    return retl;
  }

  if (class(NWIN)!="real" && class(NWIN)!="num")
  {
    printf("plwins: improper first argument");
    return <<>>;
  }

  if (NWIN.nr == 1 && NWIN.nc==1)
  {
    nwin = NWIN;
  else if (NWIN.nc == 2)
  {
    nwin = NWIN.nr;
  else
    printf("plwins: improper first argument");
    return <<>>;
  }}

  if (nwin==0)
  { return plend(); }

  if (nwin>32)
  {
    printf("plwins: improper first argument");
    return <<>>;
  }

  //
  // process devices
  //
  if (!exist(dev))
  { dev = "xwin"; }

  //
  // process sizes
  //
  if (!exist(sz))
  { sz = [PLPLOT_WIN_LEN_X_DEFAULT, PLPLOT_WIN_LEN_Y_DEFAULT]; }
  if (exist(sz) && (sz.nr*sz.nc!=2 || class(sz) != "num"))
  { sz = [PLPLOT_WIN_LEN_X_DEFAULT, PLPLOT_WIN_LEN_Y_DEFAULT]; }

  //
  // process offsets
  //
  if (!exist(offs))
  { offs = [PLPLOT_WIN_OFFSET_X_DEFAULT, PLPLOT_WIN_OFFSET_Y_DEFAULT]; }
  if (exist(sz) && (sz.nr*sz.nc!=2 || class(sz) != "num"))
  { offs = [PLPLOT_WIN_OFFSET_X_DEFAULT, PLPLOT_WIN_OFFSET_Y_DEFAULT]; }

  //
  // too many open windows
  //
  d = plwins();
  wavl = d.win;
  davl = d.dev;
  if (nwin < length(wavl))
  {
    // kill extra windows
    for (i in (nwin+1):length(wavl))
    {
      plwin(i);
      plclose();
      clear(PLPLOT_PLOT_WINDOWS.[i]);
    }
  }

  //
  // open/modify requested number of plot windows
  //
  d = plwins();
  wavl = d.win;
  davl = d.dev;
  for (i in 1:nwin)
  {
    ndev = dev[ min(i,length(dev)) ];

    nx = 1;
    ny = 1;
    if (NWIN.nc==2)
    {
      nx = NWIN[i;1];
      ny = NWIN[i;2];
    }

    if (!exist(PLPLOT_PLOT_WINDOWS.[i]))
    {
      plstart(nx,ny,ndev,sz,offs);
    else
      // window exists and is open, check if the specifications have changed
      // and if so close it and open a new one
      if ( (PLPLOT_PLOT_WINDOWS.[i].nx != nx) || (PLPLOT_PLOT_WINDOWS.[i].ny != ny) )
      {
        plclose();
        clear(PLPLOT_PLOT_WINDOWS.[i]);
        plstart(nx,ny,ndev,sz,offs);
      }
    } // if (!exist(PLPLOT_PLOT_WINDOWS.[i]))
  } // for (i in 1:nwin)
  plwin(i);

  return 0;
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Show the current plot-window, and the possibilities
//

showplwin = function (w)
{
  if (length(PLPLOT_PLOT_WINDOWS) == 0)
  {
    printf ("No plot objects\n");
    return 0;
  }

  printf ("Current plot-window is:\t\t%i\n", PLPLOT_ACTIVE_WIN);
  printf ("Available plot windows are:\t");
  for (i in members (PLPLOT_PLOT_WINDOWS))
  {
    printf ("%s   ", i);
  }
  printf ("\n");

  if (exist (w))
  {
    for (i in members (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN]))
    {
      if ( class(PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].[i])=="num" || ...
           class(PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].[i])=="string" )
      {
        printf("%s: [%s]\n", i, join(PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].[i],","));
        continue;
      }
      printf("%s: <<%s>>\n", i, join(members(PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].[i]),";"));
    }
  }
};

getplot = function ( win_no )
{
  local (win_no)

      if (length(PLPLOT_PLOT_WINDOWS) != 0)
  {
    if (!exist (win_no))
    { win_no = PLPLOT_ACTIVE_WIN; }

    if (exist (PLPLOT_PLOT_WINDOWS.[win_no]))
    {
      return (PLPLOT_PLOT_WINDOWS.[win_no]);
    else
      return 0;
    }
  }
  return <<>>;
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Do a multiplot window
//
plmultiplot = function(s)
{
  THIS_SOLVER = "plmultiplot";
  if (!PLPLOT_ACTIVE_WIN)
  {
    printf("%s: Cannot be called before plwins()!\n", THIS_SOLVER);
    return 1;
  }

  if (class(s)!="num")
  { return 2; }

  if ((s.nc)!=4)
  { return 3; }

  if ((s.nr) != PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].nplot)
  { return 4; }

  PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].plvpor = s;

  return 0;
};


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Set the plot-window orientation (portrait, landscape, rotated).
//

plsori = function ( rot )
{
  _plsori (rot);
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Set/start/select the plot device
//

plstart = function ( nx, ny, dev, sz, offs )
{
  if (!exist (nx))
  { nx = 1; }
  if (!exist (ny))
  { ny = 1; }
  if (!exist (dev))
  { dev = "xwin"; }
  if (!exist (sz))
  { sz = [PLPLOT_WIN_LEN_X_DEFAULT, PLPLOT_WIN_LEN_Y_DEFAULT]; }
  if (!exist (offs))
  { offs = [PLPLOT_WIN_OFFSET_X_DEFAULT, PLPLOT_WIN_OFFSET_Y_DEFAULT]; }

  // Create the plot-object
  PLPLOT_ACTIVE_WIN = _plmkstrm ();
  create_plot_object (PLPLOT_ACTIVE_WIN, nx, ny);

  // Default window size and location for X-windows driver
  _plspage (PLPLOT_WIN_XP_DEFAULT, PLPLOT_WIN_YP_DEFAULT, sz[1], sz[2], offs[1], offs[2]);

  // update color scheme if such has changed
  modify_color0_scheme(PLPLOT_COLOR_SCHEME);

  // Start up the plot-window
  _plstart (dev, nx, ny);

  // Set the pen-width
  _plwid (8);

  // Turn between plot pause off
  _plspause (0);
  _pltext ();

  return PLPLOT_ACTIVE_WIN;
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Close a plot device. We must destroy the current plot-object
// And switch the output stream back to the default.
//

plclose = function ()
{
  if (size (PLPLOT_PLOT_WINDOWS) > 1)
  {
    //
    // Clear PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN] and reset P to 1st plot-window
    //

    clear (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN]);
    _plend1 ();
    _plsstrm (strtod (members (PLPLOT_PLOT_WINDOWS)[1]));
    PLPLOT_ACTIVE_WIN = strtod (members (PLPLOT_PLOT_WINDOWSIN)[1]);
    return PLPLOT_ACTIVE_WIN;

  else if (size (PLPLOT_PLOT_WINDOWS) == 1)
  {

    if (exist (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN]))
    {
      clear (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN]);
      PLPLOT_ACTIVE_WIN=0;
    }
    _plend1 ();
    return [];

  else if (size (PLPLOT_PLOT_WINDOWS) == 0)
  {
    PLPLOT_ACTIVE_WIN=0;
    return 0;
  }}}
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Close ALL the plot-windows
//
plend = function ()
{
  _plend ();
  if (exist (PLPLOT_PLOT_WINDOWS))
  { clear (PLPLOT_PLOT_WINDOWS); }

  PLPLOT_ACTIVE_WIN=0;
  PLPLOT_PLOT_WINDOWS = <<>>;

  return PLPLOT_ACTIVE_WIN;
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Change plot aspect ratio
//

plaspect = function ( aspect )
{
  check_plot_object ();
  i = mod (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].subplot, PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].nplot) + 1;
  if (!exist (aspect))
  {
    PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].aspect[i] = 0;
  else
    PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].aspect[i] = aspect;
  }
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Change plot line style
//

plstyle = function ( style )
{
  check_plot_object ();
  i = mod (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].subplot, PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].nplot) + 1;

  if (exist (style))
  {
    if (class (style) == "string")
    {
      for (j in 1:style.n)
      {
	if (style[j] != "line" && ...
            style[j] != "point" && ...
            style[j] != "line-point")
	{
	  error ("plstyle: STYLE must be either " + ...
	         "\"point\", \"line\" or \"line-point\"");
	}
      }
      PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].style.[i] = style;
    }
    return 1;
  }
  PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].style.[i] = "line";
  return 1;
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Control of the plot line style.
// There are 8 line styles

plline = function ( line_style )
{
  THIS_SOLVER = "plline";
  check_plot_object ();
  i = mod (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].subplot, PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].nplot) + 1;

  if (exist (line_style))
  {
    if (class (line_style) == "num")
    {
      if (line_style.nc != 8)
      {
        error ("plpoint: LVEC must be 1x8 in size");
      }
      PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].lstyle[i;] = line_style;
    }
    return 1;
  }
  PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].lstyle[i;] = 1:8;
  return 1;
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Control of the plot point style.
// There are 8 line styles

plpoint = function ( point_style )
{
  check_plot_object ();
  i = mod (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].subplot, PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].nplot) + 1;

  if (exist (point_style))
  {
    if (class (point_style) == "num")
    {
      if (point_style.nc != 8)
      {
	error ("plpoint: PVEC must be 1x8 in size");
      }
      PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].pstyle[i;] = point_style;
    }
    return 1;
  }
  PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].pstyle[i;] = 1:8;
  return 1;
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Get the right value of line-style
//

get_style = function ( STY, K )
{
  local (sty);
  sty = mod(K, STY.n);
  if(sty == 0)
  {
    sty = STY.n;
  }
  return STY[sty];
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Change fonts
//

plfont = function ( font )
{
  check_plot_object ();
  i = mod (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].subplot, PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].nplot) + 1;

  if (!exist (font)) { font = 1; }

  if (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].fontld == 0)
  {
    _plfontld (1);
    PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].fontld = 1;
  }

  PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].font[i] = font;
  return PLPLOT_ACTIVE_WIN;
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Change pen width
//

plwid = function ( width )
{
  check_plot_object ();
  i = mod (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].subplot, PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].nplot) + 1;

  if (!exist (width))
  {
    PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].width = ones(1,32);
  else
    PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].width[i;] = width;
  }

  return PLPLOT_ACTIVE_WIN;
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Place some text on the plot
//
pltext = function (s, loc, vel, opts )
{
  check_plot_object ();
  i = mod (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].subplot, PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].nplot) + 1;

  if (strlen(s)<1)
  { return 1; }

  if (length(loc)!=2)
  { return 1; }
  x = loc[1];
  y = loc[2];

  if (!exist(vel))
  { vel=1; }

  incl = [1,0];
  just = 0.5;
  c    = 1;
  if (exist(opts))
  {
    if (exist(opts.incl))
    {
      if (length(opts.incl)==2)
      { incl = opts.incl; }
    }
    if (exist(opts.just))
    {
      if (class(opts.just)=="num")
      {
        just = opts.just;
      }
      if (class(opts.just)=="string")
      {
        if (strindex(opts.just,"c"))
        { just = 0.5;}
        if (strindex(opts.just,"l"))
        { just = 0.0;}
        if (strindex(opts.just,"r"))
        { just = 1.0;}
      }
    }
    if (exist(opts.color))
    {
      c = opts.color;
    }
  }

  dx = x + incl[1];
  dy = incl[2];

  rval = <<>>;
  rval.x = x;
  rval.y = y;
  rval.dx = dx;
  rval.dy = dy;
  rval.scale = vel;
  rval.just  = just;
  rval.color = c;
  rval.text  = s;

  PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].pltext.[i] = ...
      [PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].pltext.[i], rval];

  return 0;
};

plptex = function ( text, x , y , dx , dy , just )
{
  if (!check_plot_object ())
  {
    printf ("Must use plot() before plptex()\n");
    return 0;
  }

  if (!exist (x)) { x = 0; }
  if (!exist (y)) { y = 0; }
  if (!exist (dx)) { dx = abs(x)+1; }
  if (!exist (dy)) { dy = 0; }
  if (!exist (just)) { just = 0; }

  _plptex (x, y, dx, dy, just, text);
  _plflush ();
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Set up the viewing altitude for 3-D plots
//
plalt = function ( ALT )
{
  check_plot_object ();
  i = mod (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].subplot, PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].nplot) + 1;

  if (exist (ALT))
  {
    PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].alt[i] = ALT;
  else
    PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].alt[i] = 60;
  }
  return PLPLOT_ACTIVE_WIN;
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Set the viewing azimuth for 3-D plots
//
plaz = function ( AZ )
{
  check_plot_object ();
  i = mod (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].subplot, PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].nplot) + 1;

  if (exist (AZ))
  {
    PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].az[i] = AZ;
  else
    PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].az[i] = 45;
  }
  return PLPLOT_ACTIVE_WIN;
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Sort list element names/labels by numeric order, then string order.
//

list_sort = function ( L )
{
  tl = <<>>;
  j = k = 1;

  for (i in members (L))
  {
    if (!isnan (strtod (i)))
    {
      num[j] = i;
      j++;
      else
      char[k] = i;
      k++;
    }
  }

  // Sort the numeric labels

  if (exist (num))
  {
    num = sort (strtod (num)).val;
    tl.num = num;
  }

  if (exist (char))
  {
    tl.char = char;
  }

  return tl;
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Set the subplot, this overides the action in plot().
//

subplot = function ( sub )
{
  check_plot_object ();

  if (!exist (sub))
  {
    subplot_f = 0;
    _pladv (0);
  else
    if (sub > PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].nplot)
    {
      error ("Current window does not have this many subplots");
    }
    if (sub > 0)
    {
      PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].subplot = sub - 1;
      subplot_f = 1;
      _pladv (sub);
    else
      if (sub == 0)
      {
        // Do not advance, stay at current subplot
        PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].subplot = PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].subplot - 1;
        subplot_f = 1;
      }
    }
  }
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// format data for plotting
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
plformat = function ( fmts )
{
  THIS_SOLVER = "plformat";

  if (!exist(fmts))
  { return 1; }
  if (class(fmts)!="string")
  { return 1; }

  if (!PLPLOT_ACTIVE_WIN)
  {
    plwins(1);
    plwin (1);
  }
  // The current index
  p = mod (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].subplot, PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].nplot) + 1;

  _plformat = zeros(length(fmts),length(PLPLOT_DATASET_PLOT_FORMAT));
  _axes = zeros(1,length(fmts));
  _cols = cell (1,length(fmts));
  for (i in 1:length(fmts))
  {
    rval = convert_gnuformat_to_plformat( fmts[i] );
    if (exist(rval.plformat))
    { _plformat[i;] = rval.plformat; }
    if (exist(rval.axes))
    { _axes[i] = rval.axes; }
    if (exist(rval.use_cols))
    {
      _cols[i] = rval.use_cols;
    else
      _cols[i] = []; // default is to use column 1 and 2
    }
  }

  // find all lines that have [-1,-1] for lt and lc
  // and change that to increasing type and color
  _il = find((_plformat[;2]==-1) && (_plformat[;3]==-1));
  if (length(_il) > 0)
  {
    _lt = PLPLOT_DATASET_PLOT_FORMAT[2];
    _lc = PLPLOT_DATASET_PLOT_FORMAT[3];
    k=0;
    for (_i in _il)
    {
      _plformat[_i;2] = _lt;
      _plformat[_i;3] = _lc;
      _lc++;
      if (_lc > 15)
      {
        _lc = 1;
        _lt = _lt + 1;
        if (_lt > 8)
        { _lt = 1; }
      }
    }
  }
  // find all points that have [-1,-1] for pt and pc
  // and change that to increasing type and color
  _il = find((_plformat[;5]==-1) && (_plformat[;6]==-1));
  if (length(_il) > 0)
  {
    _pt = PLPLOT_DATASET_PLOT_FORMAT[5];
    _pc = PLPLOT_DATASET_PLOT_FORMAT[6];
    k=0;
    for (_i in _il)
    {
      _plformat[_i;5] = _pt;
      _plformat[_i;6] = _pc;
      _pc++;
      if (_pc > 15)
      {
        _pc = 1;
        _pt = _pt + 1;
        if (_pt > 16)
        { _pt = 1; }
      }
    }
  }

  PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].plformat.[p] = _plformat;
  PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].axes.[p] = _axes;
  PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].cols.[p] = _cols;
  return 0;
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Plot data matrix or list as a scatter or line plot, or histogram (barplot)
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
plplot = function ( data, fn )
{
  THIS_SOLVER = "plplot";

  if (!PLPLOT_ACTIVE_WIN)
  {
    plwins(1);
    plwin (1);
  }
  p = mod (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].subplot, PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].nplot) + 1; // The current index

  _plgra ();
  _plfont (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].font[p]);
  _plwid  (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].width[p]);
  _plscolbg(PLPLOT_COLOR_RGB[1;]);
  _pllsty (PLPLOT_DEFAULT_LINE_STYLE);
  _plcol (1);

  if (subplot_f>0)
  {
    // Advance 1 subplot
    _pladv (0);
  else
    subplot_f = 0;     // The user has set the subplot

    //
    // Set the aspect ratio: _pladv (0) -> pushed into plwin()
    //
    if (all(!isnan(PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].plvpor)))
    {
      _plvpor(PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].plvpor);
    else
      if (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].aspect[p] != 0)
      {
        _plvasp (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].aspect[p]);
      else
        _plvsta ();
      }
    }
  }

  //
  // Draw the graph
  //  Step through the matrix plotting
  //  each column versus the 1st
  //
  if (class (data) == "num" || ishist(data))
  {
    // plot a single matrix or a special list: histogram
    if (ishist(data))
    { data = hist_line(data); }

    //
    // Set up the plot basics
    //
    if (data.nc == 1)
    { data = [(1:data.nr)', data]; }

    //
    // is there any plformat info on how to plot it
    // 
    if (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].plformat.[p].nr==0)
    {
      //
      // no format:
      // 
      //  plot columns 2:data.nc against the first
      _cols = cell(1, (data.nc)-1);
      _xcol_idx = 1;
      _ycol_idx = [2:(data.nc)];
      for (_i in _ycol_idx)
      { _cols[_i-1] = [_xcol_idx, _i]; }
      // no format:
      //  assume x1y1 coordinate system
      _axes = ones(1, (data.nc)-1);
      // no format: 
      //  construct plformats for each column
      _plformat = [];
      _fmt = PLPLOT_DATASET_PLOT_FORMAT;
      _lt = PLPLOT_DATASET_PLOT_FORMAT[2];
      _lc = PLPLOT_DATASET_PLOT_FORMAT[3];
      for (_i in _ycol_idx)
      {
        _fmt[2] = _lt;
        _fmt[3] = _lc;
        _lc++;
        if (_lc > 14)
        {
          _lc = 1;
          _lt = _lt + 1;
          if (_lt > 8)
          { _lt = 1; }
        }
        _plformat = [_plformat; _fmt ];
      }
    else
      // user knows what she is doing:
      //  assume that length of plformat determines number of plots
      //  and go along with that assumption
      _plformat = PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].plformat.[p];
      _axes = PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].axes.[p];
      _cols = PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].cols.[p];
      if (all(_axes == PLPLOT_AXES_REV.x1y1))
      {
        // all data sets on single axis. but user forgot to specify what
        // is plotted against what: 
        //  assume columns 2:data.nc against the first
        _idx_x1 = 1:(_plformat.nr);
        _idx_y1 = 1:(_plformat.nr);
        _xcol_idx = [];
        _ycol_idx = [];
        for (_i in 1:(_plformat.nr))
        {
          if (isempty(_cols[_i]))
          {
            if(_plformat[_i;8]==0)
            {
              // no error plots
              _cols[_i] = [1,min(_i+1,data.nc)];
            else if ((_plformat[_i;8]==1)||(_plformat[_i;8]==3))
            {
              // errorx or errory
              _cols[_i] = [1,min(_i+1,data.nc),min(_i+2,data.nc),min(_i+3,data.nc)];
            else
              // errorxy
              _cols[_i] = [1,min(_i+1,data.nc),min(_i+2,data.nc),min(_i+3,data.nc),...
                  min(_i+4,data.nc),min(_i+5,data.nc)];
            }}
          }
          _xcol_idx = unique([_xcol_idx, _cols[_i][1]]);
          _ycol_idx = unique([_ycol_idx, _cols[_i][2]]);
        }
        _idx_x2 = _idx_x1;
        _idx_y2 = _idx_y1;
        _alt_xcol_idx = _xcol_idx;
        _alt_ycol_idx = _ycol_idx;

      else
        // user has multiple axis: then everything has to be specified for
        // each plotted dataset
        _idx_x1 = find((_axes == PLPLOT_AXES_REV.x1y1) || (_axes == PLPLOT_AXES_REV.x1y2));
        _idx_y1 = find((_axes == PLPLOT_AXES_REV.x1y1) || (_axes == PLPLOT_AXES_REV.x2y1));
        _idx_x2 = find((_axes == PLPLOT_AXES_REV.x2y1) || (_axes == PLPLOT_AXES_REV.x2y2));
        _idx_y2 = find((_axes == PLPLOT_AXES_REV.x1y2) || (_axes == PLPLOT_AXES_REV.x2y2));
        _xcol_idx = [];
        for (_i in _idx_x1)
        {
          if (isempty(_cols[_i]))
          { _cols[_i] = [1,min(_i+1,data.nc)]; }
          _xcol_idx = unique([_xcol_idx, _cols[_i][1]]);
        }
        _ycol_idx = [];
        for (_i in _idx_y1)
        {
          if (isempty(_cols[_i]))
          { _cols[_i] = [1,min(_i+1,data.nc)]; }
          _ycol_idx = unique([_ycol_idx, _cols[_i][2]]);
        }
        _alt_xcol_idx = [];
        for (_i in _idx_x1)
        {
          if (isempty(_cols[_i]))
          { _cols[_i] = [1,min(_i+1,data.nc)]; }
          _alt_xcol_idx = unique([_alt_xcol_idx, _cols[_i][1]]);
        }
        _alt_ycol_idx = [];
        for (_i in _idx_y2)
        {
          if (isempty(_cols[_i]))
          { _cols[_i] = [1,min(_i+1,data.nc)]; }
          _alt_ycol_idx = unique([_alt_ycol_idx, _cols[_i][2]]);
        }
      }
    }

    //
    // now plot for each set of axes
    //
    _idx_x1y1 = find(_axes == PLPLOT_AXES_REV.x1y1);
    if (length(_idx_x1y1) < length(_axes))
    {
      _idx_x1y2 = find(_axes == PLPLOT_AXES_REV.x1y2);
      if (length(_idx_x1y1)+length(_idx_x1y2) < length(_axes))
      {
        _idx_x2y1 = find(_axes == PLPLOT_AXES_REV.x2y1);
        if (length(_idx_x1y1)+length(_idx_x1y2)+length(_idx_x2y1) < length(_axes))
        {
          _idx_x2y2 = find(_axes == PLPLOT_AXES_REV.x2y1);
          if (length(_idx_x1y1)+length(_idx_x1y2)+length(_idx_x2y1)+length(_idx_x2y2) != length(_axes))
          {
            printf("What is going on? Help! I need somebody! Help! Just anybody! Help!\n");
            return 1;
          }
        }
      }
    }
    if (!isempty(_idx_x1y1))
    {
      // process x range
      _x = <<>>;
      if ( (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].xmin[p]==1j) ...
           || (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].xmax[p]==1j) )
      { _x = find_scales ("x", data, _xcol_idx, p); }
      if (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].xmin[p]!=1j)
      { _x.min = PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].xmin[p]; }
      if (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].xmax[p]!=1j)
      { _x.max = PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].xmax[p]; }
      // process y range
      _y = <<>>;
      if ( (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].ymin[p]==1j) ...
            || (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].ymax[p]==1j) )
      { _y = find_scales ("y", data, _ycol_idx, p); }
      if (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].ymin[p]!=1j)
      { _y.min = PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].ymin[p]; }
      if (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].ymax[p]!=1j)
      { _y.max = PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].ymax[p]; }
      // logs
      _logx = 0;
      if(strindex(PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].gridx[p],"l")>0)
      {
        _logx = 1;
        _x.min = log10(_x.min);
        _x.max = log10(_x.max);
      }
      _logy = 0;
      if(strindex(PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].gridy[p],"l")>0)
      {
        _logy = 1;
        _y.min = log10(_y.min);
        _y.max = log10(_y.max);
      }

      // do we use x2? if not, plot the opposite edge of frame now
      _idx_x2 = find((_axes == PLPLOT_AXES_REV.x2y1) || (_axes == PLPLOT_AXES_REV.x2y2));
      _fmt_x = "";
      if (length(_idx_x2)==0)
      { _fmt_x = "c"; }
      // do we use y2? if not, plot the opposite edge of frame now
      _fmt_y = "";
      _idx_y2 = find((_axes == PLPLOT_AXES_REV.x1y2) || (_axes == PLPLOT_AXES_REV.x2y2));
      if (length(_idx_y2)==0)
      { _fmt_y = "c"; }

      // do we resize tic marks?
      scale = 1;
      if (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].xtick_scale.[p][1]!=1)
      { scale = PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].xtick_scale.[p][1]; }
      if (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].ytick_scale.[p][1]!=1)
      { scale = PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].ytick_scale.[p][1]; }
      if (scale!=1)
      {
        def_scale_default = _plchr();
        _plchr([0, scale]);
      }

      _plwind (_x.min, _x.max, _y.min, _y.max);
      _plbox (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].gridx[p] + _fmt_x, ...
          PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].xtick.[p][1], ...
              PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].xtick.[p][2], ...
                  PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].gridy[p] + _fmt_y, ...
                      PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].ytick.[p][1], ...
                          PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].ytick.[p][2]);

      if (scale !=1)
      { _plchr(def_scale_default); } // revert to the earlier text scaling

      for (_i in range(_idx_x1y1))
      { _plprintf(data, _plformat[_idx_x1y1[_i];], _logx, _logy, _cols[_idx_x1y1[_i]]); }
      _pllsty (PLPLOT_DEFAULT_LINE_STYLE);
      _plcol (1);
      _plwid (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].width[p]);
    }
    if (!isempty(_idx_x1y2))
    {
      // process x range
      _x = <<>>;
      if ( (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].xmin[p]==1j) ...
            || (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].xmax[p]==1j) )
      { _x = find_scales ("x", data, _xcol_idx, p); }
      if (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].xmin[p]!=1j)
      { _x.min = PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].xmin[p]; }
      if (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].xmax[p]!=1j)
      { _x.max = PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].xmax[p]; }
      // process y range
      _y = <<>>;
      if ( (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].alt_ymin[p]==1j) ...
            || (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].alt_ymax[p]==1j) )
      { _y = find_scales ("y", data, _alt_ycol_idx, p); }
      if (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].alt_ymin[p]!=1j)
      { _y.min = PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].alt_ymin[p]; }
      if (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].alt_ymax[p]!=1j)
      { _y.max = PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].alt_ymax[p]; }

      _logx = 0;
      if(strindex(PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].gridx[p],"l")>0)
      {
        _logx = 1;
        _x.min = log10(_x.min);
        _x.max = log10(_x.max);
      }
      _logy = 0;
      if(strindex(PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].alt_gridy[p],"l")>0)
      {
        _logy = 1;
        _y.min = log10(_y.min);
        _y.max = log10(_y.max);
      }

      // do we use x1? if not, plot the opposite edge of frame now
      _idx_x2 = find((_axes == PLPLOT_AXES_REV.x2y1) || (_axes == PLPLOT_AXES_REV.x2y2));
      _fmt_x = "";
      if (length(_idx_x2)==0)
      { _fmt_x = "c"; }
      // do we use y1? if not, plot the opposite edge of frame now
      _fmt_y = "";
      _idx_y1 = find((_axes == PLPLOT_AXES_REV.x1y1) || (_axes == PLPLOT_AXES_REV.x2y1));
      if (length(_idx_y1)==0)
      { _fmt_y = "c"; }

      scale = 1;
      if (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].xtick_scale.[p][1]!=1)
      { scale = PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].xtick_scale.[p][1]; }
      if (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].alt_ytick_scale.[p][1]!=1)
      { scale = PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].alt_ytick_scale.[p][1]; }
      if (scale!=1)
      {
        def_scale_default = _plchr();
        _plchr([0, scale]);
      }

      _plwind (_x.min, _x.max, _y.min, _y.max);
      _plbox (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].gridx[p] + _fmt_x, ...
          PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].xtick.[p][1], ...
              PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].xtick.[p][2], ...
                  PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].alt_gridy[p] + _fmt_y, ...
                      PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].alt_ytick.[p][1], ...
                          PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].alt_ytick.[p][2]);

      if (scale !=1)
      { _plchr(def_scale_default); } // revert to the earlier text scaling

      for (_i in range(_idx_x1y2))
      { _plprintf(data, _plformat[_idx_x1y2[_i];], _logx, _logy, _cols[_idx_x1y2[_i]]); }
      _pllsty (PLPLOT_DEFAULT_LINE_STYLE);
      _plcol (1);
      _plwid (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].width[p]);
    }
    if (!isempty(_idx_x2y1))
    {
      // process x range
      _x = <<>>;
      if ( (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].alt_xmin[p]==1j) ...
            || (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].alt_xmax[p]==1j) )
      { _x = find_scales ("x", data, _alt_xcol_idx, p); }
      if (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].alt_xmin[p]!=1j)
      { _x.min = PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].alt_xmin[p]; }
      if (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].alt_xmax[p]!=1j)
      { _x.max = PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].alt_xmax[p]; }
      // process y range
      _y = <<>>;
      if ( (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].ymin[p]==1j) ...
            || (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].ymax[p]==1j) )
      { _y = find_scales ("y", data, _ycol_idx, p); }
      if (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].ymin[p]!=1j)
      { _y.min = PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].ymin[p]; }
      if (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].ymax[p]!=1j)
      { _y.max = PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].ymax[p]; }
      // logs
      _logx = 0;
      if(strindex(PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].alt_gridx[p],"l")>0)
      {
        _logx = 1;
        _x.min = log10(_x.min);
        _x.max = log10(_x.max);
      }
      _logy = 0;
      if(strindex(PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].gridy[p],"l")>0)
      {
        _logy = 1;
        _y.min = log10(_y.min);
        _y.max = log10(_y.max);
      }

      // do we use x1? if not, plot the opposite edge of frame now
      _idx_x1 = find((_axes == PLPLOT_AXES_REV.x1y1) || (_axes == PLPLOT_AXES_REV.x1y2));
      _fmt_x = "";
      if (length(_idx_x1)==0)
      { _fmt_x = "c"; }
      // do we use y2? if not, plot the opposite edge of frame now
      _fmt_y = "";
      _idx_y2 = find((_axes == PLPLOT_AXES_REV.x1y2) || (_axes == PLPLOT_AXES_REV.x2y2));
      if (length(_idx_y2)==0)
      { _fmt_y = "c"; }

      scale = 1;
      if (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].alt_xtick_scale.[p][1]!=1)
      { scale = PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].alt_xtick_scale.[p][1]; }
      if (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].ytick_scale.[p][1]!=1)
      { scale = PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].ytick_scale.[p][1]; }
      if (scale!=1)
      {
        def_scale_default = _plchr();
        _plchr([0, scale]);
      }

      _plwind (_x.min, _x.max, _y.min, _y.max);
      _plbox (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].alt_gridx[p] + _fmt_x, ...
          PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].alt_xtick.[p][1], ...
              PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].alt_xtick.[p][2], ...
                  PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].gridy[p] + _fmt_y, ...
                      PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].ytick.[p][1], ...
                          PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].ytick.[p][2]);

      if (scale !=1)
      { _plchr(def_scale_default); } // revert to the earlier text scaling

      for (_i in range(_idx_x2y1))
      { _plprintf(data, _plformat[_idx_x1y1[_i];], _logx, _logy, _cols[_idx_x2y1[_i]]); }
      _pllsty (PLPLOT_DEFAULT_LINE_STYLE);
      _plcol (1);
      _plwid (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].width[p]);
    }
    _idx_x2y2 = find(_axes == PLPLOT_AXES_REV.x2y2);
    if (!isempty(_idx_x2y2))
    {
      // process x range
      _x = <<>>;
      if ( (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].alt_xmin[p]==1j) ...
            || (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].alt_xmax[p]==1j) )
      { _x = find_scales ("x", data, _alt_xcol_idx, p); }
      if (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].alt_xmin[p]!=1j)
      { _x.min = PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].alt_xmin[p]; }
      if (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].alt_xmax[p]!=1j)
      { _x.max = PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].alt_xmax[p]; }
      // process y range
      _y = <<>>;
      if ( (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].alt_ymin[p]==1j) ...
            || (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].alt_ymax[p]==1j) )
      { _y = find_scales ("y", data, _alt_ycol_idx, p); }
      if (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].alt_ymin[p]!=1j)
      { _y.min = PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].alt_ymin[p]; }
      if (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].alt_ymax[p]!=1j)
      { _y.max = PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].alt_ymax[p]; }
      // logs
      _logx = 0;
      if(strindex(PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].alt_gridx[p],"l")>0)
      {
        _logx = 1;
        _x.min = log10(_x.min);
        _x.max = log10(_x.max);
      }
      _logy = 0;
      if(strindex(PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].alt_gridy[p],"l")>0)
      {
        _logy = 1;
        _y.min = log10(_y.min);
        _y.max = log10(_y.max);
      }

      // do we use x1? if not, plot the opposite edge of frame now
      _idx_x1 = find((_axes == PLPLOT_AXES_REV.x1y1) || (_axes == PLPLOT_AXES_REV.x1y2));
      _fmt_x = "";
      if (length(_idx_x1)==0)
      { _fmt_x = "c"; }
      // do we use y1? if not, plot the opposite edge of frame now
      _fmt_y = "";
      _idx_y1 = find((_axes == PLPLOT_AXES_REV.x1y1) || (_axes == PLPLOT_AXES_REV.x2y1));
      if (length(_idx_y1)==0)
      { _fmt_y = "c"; }

      scale = 1;
      if (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].alt_xtick_scale.[p][1]!=1)
      { scale = PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].alt_xtick_scale.[p][1]; }
      if (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].alt_ytick_scale.[p][1]!=1)
      { scale = PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].alt_ytick_scale.[p][1]; }
      if (scale!=1)
      {
        def_scale_default = _plchr();
        _plchr([0, scale]);
      }

      _plwind (_x.min, _x.max, _y.min, _y.max);
      _plbox (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].alt_gridx[p] + _fmt_x, ...
          PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].alt_xtick.[p][1], ...
              PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].alt_xtick.[p][2], ...
                  PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].alt_gridy[p] + _fmt_y, ...
                      PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].alt_ytick.[p][1], ...
                          PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].alt_ytick.[p][2]);

      if (scale !=1)
      { _plchr(def_scale_default); } // revert to the earlier text scaling

      for (_i in range(_idx_x2y2))
      { _plprintf(data, _plformat[_idx_x2y2[_i];], _logx, _logy, _cols[_idx_x2y2[_i]]); }
      _pllsty (PLPLOT_DEFAULT_LINE_STYLE);
      _plcol (1);
      _plwid (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].width[p]);
    }

  else if (class (data) == "list")
  {
    if (length(data)==0)
    { error (THIS_SOLVER + ": How would you like me to plot an empty list?\n"); }

    //
    // is there any plformat info on how to plot it
    // 
    if (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].plformat.[p].nr==0)
    {
      //
      // no format:
      // 
      //  plot first two columns in each dataset
      _cols = cell(1, length(data));
      _xcol_idx = 1;
      _ycol_idx = 2;
      for (_i in 1:length(data))
      { _cols[_i] = [_xcol_idx, _ycol_idx]; }
      // no format:
      //  assume x1y1 coordinate system
      _axes = ones(1, length(data));
      // no format: 
      //  construct plformats for each column
      _fmt = PLPLOT_DATASET_PLOT_FORMAT;
      _plformat = _fmt;
      _lt = PLPLOT_DATASET_PLOT_FORMAT[2];
      _lc = PLPLOT_DATASET_PLOT_FORMAT[3];
      for (i in 2:length(data))
      {
        _lc++;
        _fmt[2] = _lt;
        _fmt[3] = _lc;
        if (_lc == PLPLOT_COLOR_COUNT_DEFAULT)
        {
          _lc = 1;
          _lt = _lt + 1;
          if (_lt > 8)
          { _lt = 1; }
        }
        _plformat = [_plformat; _fmt ];
      }
    else
      // user knows what she is doing:
      //  assume that length of plformat determines number of plots
      //  and go along with that assumption
      _plformat = PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].plformat.[p];
      _axes = PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].axes.[p];
      _cols = PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].cols.[p];
    }

    _s = <<>>;
    _m = members(data);
    for (_im in _m)
    {
      if (ishist(data.[_im]))
      {
        data.[_im] = hist_line(data.[_im]);
      }
    }

    //
    // now plot for each set of axes
    //
    //
    // now plot for each set of axes
    //
    _idx_x1y1 = find(_axes == PLPLOT_AXES_REV.x1y1);
    if (length(_idx_x1y1) < length(_axes))
    {
      _idx_x1y2 = find(_axes == PLPLOT_AXES_REV.x1y2);
      if (length(_idx_x1y1)+length(_idx_x1y2) < length(_axes))
      {
        _idx_x2y1 = find(_axes == PLPLOT_AXES_REV.x2y1);
        if (length(_idx_x1y1)+length(_idx_x1y2)+length(_idx_x2y1) < length(_axes))
        {
          _idx_x2y2 = find(_axes == PLPLOT_AXES_REV.x2y1);
          if (length(_idx_x1y1)+length(_idx_x1y2)+length(_idx_x2y1)+length(_idx_x2y2) != length(_axes))
          {
            printf("What is going on? Help! I need somebody! Help! Just anybody! Help!\n");
            return 0;
          }
        }
      }
    }

    _s = find_list_xy_scales(data, _cols, _axes, p);
    if (!isempty(_idx_x1y1))
    {
      _xmin = _s.xmin;
      _xmax = _s.xmax;
      _ymin = _s.ymin;
      _ymax = _s.ymax;

      // logs
      _logx = 0;
      if(strindex(PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].gridx[p],"l")>0)
      {
        _logx = 1;
        _xmin = log10(_xmin);
        _xmax = log10(_xmax);
      }
      _logy = 0;
      if(strindex(PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].gridy[p],"l")>0)
      {
        _logy = 1;
        _ymin = log10(_ymin);
        _ymax = log10(_ymax);
      }

      // do we use x2? if not, plot the opposite edge of frame now
      _idx_x2 = find((_axes == PLPLOT_AXES_REV.x2y1) || (_axes == PLPLOT_AXES_REV.x2y2));
      _fmt_x = "";
      if (length(_idx_x2)==0)
      { _fmt_x = "c"; }
      // do we use y2? if not, plot the opposite edge of frame now
      _fmt_y = "";
      _idx_y2 = find((_axes == PLPLOT_AXES_REV.x1y2) || (_axes == PLPLOT_AXES_REV.x2y2));
      if (length(_idx_y2)==0)
      { _fmt_y = "c"; }

      scale = 1;
      if (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].xtick_scale.[p][1]!=1)
      { scale = PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].xtick_scale.[p][1]; }
      if (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].ytick_scale.[p][1]!=1)
      { scale = PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].ytick_scale.[p][1]; }
      if (scale!=1)
      {
        def_scale_default = _plchr();
        _plchr([0, scale]);
      }

      _plwind (_xmin, _xmax, _ymin, _ymax);
      _plbox (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].gridx[p] + _fmt_x, ...
          PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].xtick.[p][1], ...
              PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].xtick.[p][2], ...
                  PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].gridy[p] + _fmt_y, ...
                      PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].ytick.[p][1], ...
                          PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].ytick.[p][2]);

      if (scale !=1)
      { _plchr(def_scale_default); } // revert to the earlier text scaling

      for (_i in range(_idx_x1y1))
      {
        _m_i = _m[_idx_x1y1[_i]];
        if (ishist(data.[_m_i]))
        {
          _plprintf(hist_line(data.[_m_i]), _plformat[_idx_x1y1[_i];], _logx, _logy, _cols[_idx_x1y1[_i]]);
        else
          _plprintf(data.[_m_i], _plformat[_idx_x1y1[_i];], _logx, _logy, _cols[_idx_x1y1[_i]]);
        }
      }
      _pllsty (PLPLOT_DEFAULT_LINE_STYLE);
      _plcol (1);
      _plwid (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].width[p]);
    }

    if (!isempty(_idx_x1y2))
    {
      _xmin = _s.xmin;
      _xmax = _s.xmax;
      _ymin = _s.alt_ymin;
      _ymax = _s.alt_ymax;

      // logs
      _logx = 0;
      if(strindex(PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].gridx[p],"l")>0)
      {
        _logx = 1;
        _xmin = log10(_xmin);
        _xmax = log10(_xmax);
      }
      _logy = 0;
      if(strindex(PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].alt_gridy[p],"l")>0)
      {
        _logy = 1;
        _ymin = log10(_ymin);
        _ymax = log10(_ymax);
      }

      // do we use x1? if not, plot the opposite edge of frame now
      _idx_x2 = find((_axes == PLPLOT_AXES_REV.x2y1) || (_axes == PLPLOT_AXES_REV.x2y2));
      _fmt_x = "";
      if (length(_idx_x2)==0)
      { _fmt_x = "c"; }
      // do we use y1? if not, plot the opposite edge of frame now
      _fmt_y = "";
      _idx_y1 = find((_axes == PLPLOT_AXES_REV.x1y1) || (_axes == PLPLOT_AXES_REV.x2y1));
      if (length(_idx_y1)==0)
      { _fmt_y = "c"; }

      scale = 1;
      if (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].xtick_scale.[p][1]!=1)
      { scale = PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].xtick_scale.[p][1]; }
      if (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].alt_ytick_scale.[p][1]!=1)
      { scale = PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].alt_ytick_scale.[p][1]; }
      if (scale!=1)
      {
        def_scale_default = _plchr();
        _plchr([0, scale]);
      }

      _plwind (_xmin, _xmax, _ymin, _ymax);
      _plbox (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].gridx[p] + _fmt_x, ...
          PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].xtick.[p][1], ...
              PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].xtick.[p][2], ...
                  PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].alt_gridy[p] + _fmt_y, ...
                      PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].alt_ytick.[p][1], ...
                          PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].alt_ytick.[p][2]);

      if (scale !=1)
      { _plchr(def_scale_default); } // revert to the earlier text scaling

      for (_i in range(_idx_x1y2))
      {
        _m_i = _m[_idx_x1y2[_i]];
        if (ishist(data.[_m_i]))
        {
          _plprintf(hist_line(data.[_m_i]),_plformat[_idx_x1y2[_i];], _logx, _logy, _cols[_idx_x1y2[_i]]);
        else
          _plprintf(data.[_m_i], _plformat[_idx_x1y2[_i];], _logx, _logy, _cols[_idx_x1y2[_i]]);
        }
      }
      _pllsty (PLPLOT_DEFAULT_LINE_STYLE);
      _plcol (1);
      _plwid (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].width[p]);
    }

    if (!isempty(_idx_x2y1))
    {
      _xmin = s.alt_xmin;
      _xmax = _s.alt_xmax;
      _ymin = _s.ymin;
      _ymax = _s.ymax;

      // logs
      _logx = 0;
      if(strindex(PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].alt_gridx[p],"l")>0)
      {
        _logx = 1;
        _xmin = log10(_xmin);
        _xmax = log10(_xmax);
      }
      _logy = 0;
      if(strindex(PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].alt_gridy[p],"l")>0)
      {
        _logy = 1;
        _ymin = log10(_ymin);
        _ymax = log10(_ymax);
      }

      // do we use x1? if not, plot the opposite edge of frame now
      _idx_x1 = find((_axes == PLPLOT_AXES_REV.x1y1) || (_axes == PLPLOT_AXES_REV.x1y2));
      _fmt_x = "";
      if (length(_idx_x1)==0)
      { _fmt_x = "c"; }
      // do we use y1? if not, plot the opposite edge of frame now
      _fmt_y = "";
      _idx_y2 = find((_axes == PLPLOT_AXES_REV.x1y2) || (_axes == PLPLOT_AXES_REV.x2y2));
      if (length(_idx_y2)==0)
      { _fmt_y = "c"; }

      scale = 1;
      if (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].alt_xtick_scale.[p][1]!=1)
      { scale = PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].alt_xtick_scale.[p][1]; }
      if (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].ytick_scale.[p][1]!=1)
      { scale = PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].ytick_scale.[p][1]; }
      if (scale!=1)
      {
        def_scale_default = _plchr();
        _plchr([0, scale]);
      }

      _plwind (_xmin, _xmax, _ymin, _ymax);
      _plbox (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].alt_gridx[p] + _fmt_x, ...
          PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].alt_xtick.[p][1], ...
              PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].alt_xtick.[p][2], ...
                  PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].gridy[p] + _fmt_y, ...
                      PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].ytick.[p][1], ...
                          PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].ytick.[p][2]);

      if (scale !=1)
      { _plchr(def_scale_default); } // revert to the earlier text scaling

      for (_i in range(_idx_x2y1))
      {
        _m_i = _m[_idx_x2y1[_i]];
        if (ishist(data.[_m_i]))
        {
          _plprintf(hist_line(data.[_m_i]), _plformat[_idx_x2y1[_i];], 0, 0, _cols[_idx_x2y1[_i]]);
        else
          _plprintf(data.[_m_i], _plformat[_idx_x2y1[_i];], 0, 0, _cols[_idx_x2y1[_i]]);
        }
      }
    }

    if (!isempty(_idx_x2y2))
    {
      _xmin = s.alt_xmin;
      _xmax = _s.alt_xmax;
      _ymin = _s.alt_ymin;
      _ymax = _s.alt_ymax;

      // logs
      _logx = 0;
      if(strindex(PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].alt_gridx[p],"l")>0)
      {
        _logx = 1;
        _xmin = log10(_xmin);
        _xmax = log10(_xmax);
      }
      _logy = 0;
      if(strindex(PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].alt_gridy[p],"l")>0)
      {
        _logy = 1;
        _ymin = log10(_ymin);
        _ymax = log10(_ymax);
      }

      // do we use x1? if not, plot the opposite edge of frame now
      _idx_x1 = find((_axes == PLPLOT_AXES_REV.x1y1) || (_axes == PLPLOT_AXES_REV.x1y2));
      _fmt_x = "";
      if (length(_idx_x1)==0)
      { _fmt_x = "c"; }
      // do we use y1? if not, plot the opposite edge of frame now
      _fmt_y = "";
      _idx_y1 = find((_axes == PLPLOT_AXES_REV.x1y1) || (_axes == PLPLOT_AXES_REV.x2y1));
      if (length(_idx_y1)==0)
      { _fmt_y = "c"; }

      scale = 1;
      if (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].alt_xtick_scale.[p][1]!=1)
      { scale = PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].alt_xtick_scale.[p][1]; }
      if (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].alt_ytick_scale.[p][1]!=1)
      { scale = PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].alt_ytick_scale.[p][1]; }
      if (scale!=1)
      {
        def_scale_default = _plchr();
        _plchr([0, scale]);
      }

      _plwind (_xmin, _xmax, _ymin, _ymax);
      _plbox (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].alt_gridx[p] + _fmt_x, ...
          PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].alt_xtick.[p][1], ...
              PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].alt_xtick.[p][2], ...
                  PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].alt_gridy[p] + _fmt_y, ...
                      PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].alt_ytick.[p][1], ...
                          PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].alt_ytick.[p][2]);

      if (scale !=1)
      { _plchr(def_scale_default); } // revert to the earlier text scaling

      for (_i in range(_idx_x2y2))
      {
        _m_i = _m[_idx_x2y2[_i]];
        if (ishist(data.[_m_i]))
        {
          _plprintf(hist_line(data.[_m_i]), _plformat[_idx_x2y2[_i];], 0, 0, _cols[_idx_x2y2[_i]]);
        else
          _plprintf(data.[_m_i], _plformat[_idx_x2y2[_i];], 0, 0, _cols[_idx_x2y2[_i]]);
        }
      }
    }

  else
    error ("plot: un-acceptable argument");
  }}

  //
  // Now do the legend
  //
  if (length(PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].desc.[p])>0)
  {
    desc = PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].desc.[p];
    idx_n = find(strlen(desc)>0);
    if (length(idx_n)>0)
    {
      desc_pos = PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].desc_pos.[p];
      desc_pos_xy = PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].desc_pos_xy.[p];
      desc_scale = PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].desc_scale.[p];
      _pllegend(desc[idx_n], _plformat[idx_n;], desc_scale, desc_pos, desc_pos_xy);
    }
  }

  _pllsty (PLPLOT_DEFAULT_LINE_STYLE);
  _plcol (1);
  _plwid (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].width[p]);

  // Do the axes labels in the end:
  if (strlen(PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].xlabel[p])>0)
  {
    _d = PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].xlabel_scale.[p];
    _d1 = 1;
    if (length(_d)>0)
    { _d1 = _d[1];}
    _d2 = 3.0;
    if (length(_d)>1)
    {
      _d2 = _d[2];
    else
      _d2 = _d2 / _d1;
    }
    _d3 = 0.5;
    if (length(_d)>2)
    { _d3 = _d[3];}
    _d4 = 0.5;
    if (length(_d)>3)
    { _d4 = _d[4];}

    if (_d1 !=1)
    {
      def_scale_default = _plchr();
      _plchr([0,_d1]); // revert to the earlier text scaling
    }

    _plmtex("b", _d2, _d3, _d4, PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].xlabel[p]);

    if (_d1 !=1)
    { _plchr(def_scale_default); }  // revert to the earlier text scaling
  }
  if (strlen(PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].ylabel[p])>0)
  {
    _d = PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].ylabel_scale.[p];
    _d1 = 1;
    if (length(_d)>0)
    { _d1 = _d[1];}
    _d2 = 5.0;
    if (length(_d)>1)
    {
      _d2 = _d[2];
    else
      _d2 = _d2 / _d1;
    }
    _d3 = 0.5;
    if (length(_d)>2)
    { _d3 = _d[3];}
    _d4 = 0.5;
    if (length(_d)>3)
    { _d4 = _d[4];}

    if (_d1 !=1)
    {
      def_scale_default = _plchr();
      _plchr([0,_d1]); // revert to the earlier text scaling
    }

    _plmtex("l", _d2, _d3, _d4, PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].ylabel[p]);

    if (_d1 !=1)
    { _plchr(def_scale_default); }  // revert to the earlier text scaling
  }

  // secondary x-label and the title are the same?
  if (strlen(PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].alt_xlabel[p])>0)
  {
    _d = PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].alt_xlabel_scale.[p];
    _d1 = 1;
    if (length(_d)>0)
    { _d1 = _d[1];}
    _d2 = 2.0;
    if (length(_d)>1)
    {
      _d2 = _d[2];
    else
      _d2 = _d2 / _d1;
    }
    _d3 = 0.5;
    if (length(_d)>2)
    { _d3 = _d[3];}
    _d4 = 0.5;
    if (length(_d)>3)
    { _d4 = _d[4];}

    // do we change the text scaling?
    if (_d1 !=1)
    {
      def_scale_default = _plchr();
      _plchr([0,_d1]); // revert to the earlier text scaling
    }

    _plmtex( "t", _d2, _d3, _d4, PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].alt_xlabel[p] );

    if (_d1 !=1)
    { _plchr(def_scale_default); }  // revert to the earlier text scaling

    if (strlen(PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].title[p])>0)
    { printf(THIS_SOLVER + ": Cannot have alt_xlabel and title at this time!\n"); }

  else
    // we plot the title separately from lab environment, so that we
    // can scale the text at will without modifying the axes-labels
    if (strlen(PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].title[p])>0)
    {
      // do we change the text scaling?
      def_scale_default = _plchr();
      def_scale = PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].title_scale.[p];
      _plchr(def_scale);

      _plmtex ("t", 1.0, 0.5, 0.5, PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].title[p]);

      _plchr(def_scale_default);  // revert to the earlier text scaling
    }
  }

  // secondary y-label
  if (strlen(PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].alt_ylabel[p])>0)
  {
    _d = PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].alt_ylabel_scale.[p];
    _d1 = 1;
    if (length(_d)>0)
    { _d1 = _d[1];}
    _d2 = 5.0;
    if (length(_d)>1)
    {
      _d2 = _d[2];
    else
      _d2 = _d2 / _d1;
    }
    _d3 = 0.5;
    if (length(_d)>2)
    { _d3 = _d[3];}
    _d4 = 0.5;
    if (length(_d)>3)
    { _d4 = _d[4];}

    if (_d1 !=1)
    {
      def_scale_default = _plchr();
      _plchr([0,_d1]); // revert to the earlier text scaling
    }

    _plmtex( "r", _d2, _d3, _d4, PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].alt_ylabel[p] );

    if (_d1 !=1)
    { _plchr(def_scale_default); }  // revert to the earlier text scaling
  }

  // graffitti?
  if (length(PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].pltext.[p])>0)
  {
    for (_i in 1:length(PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].pltext.[p]))
    {
      _s = PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].pltext.[p][_i];
      if (_s.scale !=1)
      {
        // do we change the text scaling?
        def_scale_default = _plchr();
        def_scale = [0, _s.scale];
        _plchr(def_scale);
      }
      _plcol (_s.color);

      if(strindex(PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].gridx[p],"l")>0)
      {
        _x  = log10(_s.x);
        _dx = log10(_s.dx);
      else
        _x  = _s.x;
        _dx = _s.dx;
      }
      _logy = 0;
      if(strindex(PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].gridy[p],"l")>0)
      {
        _y  = log10(_s.y);
        _dy = log10(_s.dy);
      else
        _y  = _s.y;
        _dy = _s.dy;
      }

      _plptex (_x, _y, _dx, _dy, _s.just, _s.text);
      _plcol (1);
      _plchr(def_scale_default);  // revert to the earlier text scaling
    }
  }

  _plflush ();

  // Increment the plot no. so that next time
  PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].subplot = PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].subplot + 1;

  if (strlen(fn)>0)
  {
    // isolate extension from fn
    ext = last(strsplt(fn, "."));
    if (strindex(ext, "ps"))
    { fmt = "psc"; }

    if (exist(fmt))
    {
      printf(THIS_SOLVER + ": Printing to file %s\n", fn);
      _plprint(fn, fmt);
      printf(THIS_SOLVER + ": Done\n");
    }
  }

  return PLPLOT_ACTIVE_WIN;
};

//
// plhold:
// Plot some data, and "hold" on for more.
// Plot the data, setting up the plot as usual the first time.
// On subsequent plots do not do any setup, just plot some
// more. plhold_off must be called to finish up.
//

plhold_first = 1; // True (1) if plhold() has NOT been used.
                  // Or if plhold_off() has been used.
                  // False (0) if plhold is in use

static (hxmin, hxmax, hymin, hymax)
plhold = function ( data, key )
{
  check_plot_object ();

  p = mod (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].subplot, PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].nplot) + 1;	// The current index
  if (!exist (key)) { key = 1; }
  if (abs (key) > data.nc)
  {
    error ("plot: KEY argument > M.nc");
  }

  if (plhold_first)
  {
    if (class (data) == "num")
    {
      //
      // Do the setup ONCE
      //
      hxmin = hxmax = hymin = hymax = 0;
      _plgra ();
      _plcol (1);
      _pllsty (PLPLOT_DEFAULT_LINE_STYLE);
      _plfont (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].font[p]);
      _plwid (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].width[p]);

      if (!subplot_f)
      {
	_pladv (0);        // Advance 1 subplot
      else
	subplot_f = 0;     // The user has set the subplot
      }

      if (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].aspect[p] != 0)
      {
	_plvasp (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].aspect[p]);
      else
	_plvsta ();
      }

      k = find ((1:data.nc) != abs (key));
      if (key > 0)
      {
	</hxmax;hxmin/> = x_scales (real(data)[;key], p);
	if (data.nc != 1)
	{
	  </hymax;hymin/> = y_scales (real(data)[;k],   p);
        else
	  </hymax;hymin/> = y_scales ((1:data.nr)',   p);
	}
      else if (key < 0) {
	</hxmax;hxmin/> = x_scales (real(data)[;k],   p);
	</hymax;hymin/> = y_scales (real(data)[;abs(key)], p);
      else
	</hxmax;hxmin/> = x_scales ((1:data.nr)', p);
	</hymax;hymin/> = y_scales (real(data),   p);
      } }

      _plwind (hxmin, hxmax, hymin, hymax);
      _plbox (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].gridx[p], ...
        PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].xtick.[p][1], ...
        PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].xtick.[p][2], ...
        PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].gridy[p], ...
        PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].ytick.[p][1], ...
        PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].ytick.[p][2]);

      _pllab (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].xlabel[p], PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].ylabel[p], PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].title[p]);
    else
      error ("plot: un-acceptable argument");
    }
    plhold_first = 0;
  }

  if (plot_matrix ( data, key, p, 0, hxmin, hxmax, hymin, hymax, hymax-hymin ) < 0)
  {
    return -1;
  }

  _plcol (1);
  _plflush ();
  _pltext ();

  return PLPLOT_ACTIVE_WIN;
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Clean up the plotting environment and get ready
// for normal interactive usage.
//

plhold_off = function ( )
{
  p = mod (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].subplot, PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].nplot) + 1;	// The current index
  plhold_first = 1;
  _plcol (1);
  _plflush ();
  _pltext ();
  PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].subplot = PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].subplot + 1;
  return PLPLOT_ACTIVE_WIN;
};


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Plot a 3-Dimensional surfaces and lines:
// 
// to plot single data set:
//    plplot( x, y, z )
//    plplot( [x,y,z] )
//    plplot( z )
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
plplot3 = function ( arg1, arg2, arg3 )
{
  _this_solver = "plplot3";

  data = <<>>;
  if (nargs==3)
  {
    //
    if ((class(arg1)=="num") && (class(arg2)=="num") && (class(arg3)=="num"))
    {
      //    plplot( x, y, z )
      data.[1] = <<>>;
      data.[1].x = arg1;
      data.[1].y = arg2;
      data.[1].z = arg3;
    else
      printf (_this_solver + ": Just to let you know I don't know how to plot this!\n");
      return 1;
    }
  else if (nargs == 1)
  {
    if (class(arg1)=="num")
    {
      if (arg1.nc > 3)
      {
        //    plplot( z )
        data.[1] = <<>>;
        data.[1].x = 1:(arg1.nr);
        data.[1].y = 1:(arg1.nc);
        data.[1].z = arg1;
      else if (arg1.nc == 3)
      {
        //    plplot( [x,y,z] )
        data.[1] = arg1;
      else
        printf (_this_solver + ": Just to let you know I don't know how to plot this!\n");
        return 1;
      }}
    else if (class(arg1)=="list")
    {
      if (exist(arg1.x) && exist(arg1.y) && exist(arg1.z))
      {
        data = <<>>;
        data.[1] = arg1;
      }
    }}
  }}

  if (!PLPLOT_ACTIVE_WIN)
  {
    plwins(1);
    plwin (1);
  }

  if (size(data) == 0)
  { return 1; }
  if (class(data) != "list")
  { return 2; }

  check_plot_object ();
  p = mod (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].subplot, PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].nplot) + 1; // The current index

  //
  // set up defaults once we start plotting data
  //
  _plgra ();
  _plfont (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].font[p]);
  _plwid  (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].width[p]);
  _plscolbg(PLPLOT_COLOR_RGB[1;]);
  _pllsty (PLPLOT_DEFAULT_LINE_STYLE);
  _plcol (1);

  if (subplot_f>0)
  {
    // Advance 1 subplot
    _pladv (0);
  else
    subplot_f = 0;     // The user has set the subplot

    //
    // Set the aspect ratio: _pladv (0) -> pushed into plwin()
    //
    if (all(!isnan(PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].plvpor[])))
    {
      _plvpor(PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].plvpor);
    else
      if (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].aspect[p] != 0)
      {
        _plvasp (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].aspect[p]);
      else
        _plvsta ();
      }
    }
  }

  _m = members(data);
  _n = length(_m);

  // replace histogrammatic entries with surface plots
  for (_im in _m)
  {
    if (ishist2(data.[_im]))
    { data.[_im] = hist2_surf(data.[_im]); }
  }

  // is there any plformat info on how to plot it
  if (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].plformat.[p].nr == 0)
  {
    _plformat = PLPLOT_DATASET_PLOT_FORMAT;
  else
    _plformat = PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].plformat.[p];
  }
  if (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].plformat.[p].nr < _n)
  {
    _fmt = _plformat[_plformat.nr; ];
    _lt = _fmt[2];
    _lc = _fmt[3];
    for (_i in 1:(_n - _plformat.nr))
    {
      _lc++;
      if (_lc > 14)
      {
        _lc = 1;
        _lt = _lt + 1;
        if (_lt > 8)
        { _lt = 1; }
      }
      _fmt[2] = _lt;
      _fmt[3] = _lc;
      _plformat = [_plformat; _fmt ];
    }
  }

  // ranges:
  // process x range
  _x = <<>>;
  if ( (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].xmin[p]==1j) ...
        || (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].xmax[p]==1j) )
  { _x = find_xyz_scales ("x", data, p); }
  if (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].xmin[p]!=1j)
  { _x.min = real(PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].xmin[p]); }
  if (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].xmax[p]!=1j)
  { _x.max = real(PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].xmax[p]); }

  // process y range
  _y = <<>>;
  if ( (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].ymin[p]==1j) ...
        || (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].ymax[p]==1j) )
  { _y = find_xyz_scales ("y", data, p); }
  if (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].ymin[p]!=1j)
  { _y.min = real(PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].ymin[p]); }
  if (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].ymax[p]!=1j)
  { _y.max = real(PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].ymax[p]); }
  // process z range
  _z = <<>>;
  if ( (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].zmin[p]==1j) ...
        || (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].zmax[p]==1j) )
  { _z = find_xyz_scales ("z", data, p); }
  if (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].zmin[p]!=1j)
  { _z.min = real(PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].zmin[p]); }
  if (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].zmax[p]!=1j)
  { _z.max = real(PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].zmax[p]); }

  // logs
  _logx = 0;
  if(strindex(PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].grid3x[p],"l")>0)
  {
    _logx = 1;
    _x.min = log10(_x.min);
    _x.max = log10(_x.max);
  }
  _logy = 0;
  if(strindex(PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].grid3y[p],"l")>0)
  {
    _logy = 1;
    _y.min = log10(_y.min);
    _y.max = log10(_y.max);
  }
  _logz = 0;
  if(strindex(PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].grid3z[p],"l")>0)
  {
    _logz = 1;
    _z.min = log10(_z.min);
    _z.max = log10(_z.max);
  }

  // do we resize tic marks?
  scale = 1;
  if (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].xtick_scale.[p][1]!=1)
  { scale = PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].xtick_scale.[p][1]; }
  if (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].ytick_scale.[p][1]!=1)
  { scale = PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].ytick_scale.[p][1]; }
  if (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].ztick_scale.[p][1]!=1)
  { scale = PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].ztick_scale.[p][1]; }

  // I don't know what this does
  xmin2d = -2.0; xmax2d = 2.0;
  ymin2d = -2.0; ymax2d = 4.0;
  _plenv (xmin2d, xmax2d, ymin2d, ymax2d, 0, -2);

  _plw3d (PLPLOT_BASE_X, PLPLOT_BASE_Y, PLPLOT_BASE_Z, ...
      _x.min, _x.max, _y.min, _y.max, _z.min, _z.max, ...
          PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].alt[p], ...
              PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].az[p]);

  if (scale!=1)
  {
    def_scale_default = _plchr();
    _plchr([0, scale]);
  }
  _plbox3 (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].grid3x[p], ...
      PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].xlabel[p], ...
          PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].xtick.[p][1], ...
              PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].xtick.[p][2], ...
                  PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].grid3y[p], ...
                      PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].ylabel[p], ...
                          PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].ytick.[p][1], ...
                              PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].ytick.[p][2], ...
                                  PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].grid3z[p], ...
                                      PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].zlabel[p], ...
                                          PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].ztick.[p][1], ...
                                              PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].ztick.[p][2]);
  if (scale !=1)
  { _plchr(def_scale_default); } // revert to the earlier text scaling

  // replace histogrammatic entries with surface plots
  for (_im in _m)
  {
    if (class(data.[_im])=="num")
    {
      _plprintf3(data.[_im], _plformat, _logx, _logy, _logz);
    else if (class(data.[_im])=="list")
    {
      if (exist(data.[_im].x) && exist(data.[_im].y) && exist(data.[_im].z))
      {
        _plprintf3(data.[_im].x, data.[_im].y, data.[_im].z, _plformat, _logx, _logy, _logz);
      }
    }}

    _pllsty (PLPLOT_DEFAULT_LINE_STYLE);
    _plcol (1);
    _plwid  (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].width[p]);
  }

  _plflush ();
  _pltext ();

  //
  // Increment the plot no. so that next time
  // we use the correct settings.
  //
  PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].subplot = PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].subplot + 1;

  return PLPLOT_ACTIVE_WIN;
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Plot contours. The data is composed in a list, with
// elements `x', `y', and `z'. x and y are single-dimension arrays
// (row or column matrices), and z is a two-dimensional array. The
// array z, is a function of x and y: z = f(x,y). Thus, the values in
// the array x can be thought of a "row-labels", and the values of y
// can be thought of as "column-lables" for the 2-dimensioal array z.
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

plcont = function ( data, myclevel )
{
  check_plot_object ();

  p = mod (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].subplot, PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].nplot) + 1;	// The current index

  _colors = PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].color;
  _idx_col = 2;

  //
  // is it a list <<x;y;z>> or a list which entries are the lists with entries <<x;y;z>>
  //
  if (    exist(data.x) && exist(data.y) && exist(data.z) ...
          &&    type(data.x)=="real" && type(data.y)=="real" && type(data.z)=="real"  )
  {

    if (data.x.n != data.z.nr || data.y.n != data.z.nc)
    {
      return 1;
    }

    xmin = xmax = ymin = ymax = zmin = zmax = 0;
    if (exist (data))
    {
      </Xmax;Xmin;Ymax;Ymin;Zmax;Zmin/> = XYZ_scales (data.x, data.y, data.z, p);
      if (Xmin < xmin) { xmin = Xmin; } if (Xmax > xmax) { xmax = Xmax; }
      if (Ymin < ymin) { ymin = Ymin; } if (Ymax > ymax) { ymax = Ymax; }
      if (Zmin < zmin) { zmin = Zmin; } if (Zmax > zmax) { zmax = Zmax; }
    }

    if (exist (myclevel))
    {
      clevel = myclevel;
    else
      clevel = linspace(zmin, zmax, 11);
    }


    _plgra ();
    _plcol  (_colors[p;1]);
    _pllsty (PLPLOT_DEFAULT_LINE_STYLE);
    _plfont (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].font[p]);
    _plwid (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].width[p]);

    //
    // Set up the 1st viewport for drawing the plot.
    //
    if (!subplot_f)
    {
      _pladv (0);        // Advance 1 subplot
    else
      subplot_f = 0;     // The user has set the subplot
    }

    _plvpas (0.15, 0.75, 0.15, 0.85, PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].aspect[p]);
    // _plvpor (0.15, 0.75, 0.15, 0.85);
    _plwind (xmin, xmax, ymin, ymax);
    _plbox (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].gridx[p], ...
        PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].xtick.[p][1], ...
        PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].xtick.[p][2], ...
        PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].gridy[p], ...
        PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].ytick.[p][1], ...
        PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].ytick.[p][2]);

    // Convert the data to log data if necessary.
    if (strindex (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].gridx[p], "l"))
    { x = log10 (real (data.x)); else x = real (data.x); }
    if (strindex (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].gridy[p], "l"))
    { y = log10 (real (data.y)); else y = real (data.y); }
    z = real (data.z);

    //
    // Draw the contours
    //

    llevel = [];
    for (i in 1:clevel.n)
    {
      k = mod (i-1, _colors.nc) + 1;
      j = mod (i-1, 8) + 1;
      _pllsty(j);
      _plcol  (_colors[p;k]);
      if (_plcont (x, y, z, 1, x.n, 1, y.n, clevel[i]))
      {
        llevel = [llevel, clevel[i]];
      }
    }

    //
    // Reset color and draw the labels.
    //
    _plcol  (_colors[p;1]);
    _pllab (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].xlabel[p], PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].ylabel[p], PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].title[p]);

    //
    // Draw the contour legend. Use a new viewport to the right
    // of the contour plot.
    //
    //_plvpas (0.75, 1.0, 0.15, 0.85, PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].aspect[p]);
    _plvpor (0.75, 1.0, 0.15, 0.85);
    _plwind (0, 1, 0, 1);

    v = 1 - 1/(2*llevel.n);
    for (i in 1:llevel.n)
    {
      xl = [0.1, 0.2, 0.3]';
      yl = [v, v, v]';
      v = v - 1/llevel.n;

      k = mod (i-1, 14) + 1;
      j = mod (i-1, 8) + 1;

      _plcol  (_colors[p;k]);
      _pllsty (j);

      _plline (3, xl, yl);
      sprintf (stmp, "%.2g", llevel[i]);
      plptex (stmp, xl[3]+.1, yl[3], , , 0);
    }

  else

    //
    // contour plot of multiple surfaces - uh oh
    //
    _ll = members(data);

    xmin = xmax = ymin = ymax = zmin = zmax = [];

    _nos_data = 0;
    for (idx_ll in range(_ll))
    {
      _ll1 = _ll[idx_ll];

      // silly user: ignore entries that are not properly formatted
      if (    !exist(data.[_ll1].x) || !exist(data.[_ll1].y) || !exist(data.[_ll1].z) ...
               ||  type(data.[_ll1].x)!="real" || type(data.[_ll1].y)!="real" || type(data.[_ll1].z)!="real"  )
      { continue; }
      if (data.[_ll1].x.n != data.[_ll1].z.nr || data.[_ll1].y.n != data.[_ll1].z.nc)
      { continue; }

      // figure out the maximum,minimum
      </Xmax;Xmin;Ymax;Ymin;Zmax;Zmin/> = XYZ_scales (data.[_ll1].x, data.[_ll1].y, data.[_ll1].z, p);
      xmin = min([xmin,Xmin]);
      xmax = max([xmax,Xmax]);
      ymin = min([ymin,Ymin]);
      ymax = max([ymax,Ymax]);
      zmin = min([zmin,Zmin]);
      zmax = max([zmax,Zmax]);

      _nos_data++;
    }

    if (!_nos_data)
    { return 1; }

    if (exist (myclevel))
    {
      clevel = myclevel;
    else
      clevel = linspace(zmin, zmax, 11);
    }

    _plgra ();
    _plcol  (_colors[p;1]);
    _pllsty (PLPLOT_DEFAULT_LINE_STYLE);
    _plfont (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].font[p]);
    _plwid (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].width[p]);

    //
    // Set up the 1st viewport for drawing the plot.
    //
    if (!subplot_f)
    {
      _pladv (0);        // Advance 1 subplot
    else
      subplot_f = 0;     // The user has set the subplot
    }

    _plvpas (0.15, 0.75, 0.15, 0.85, PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].aspect[p]);
    // _plvpor (0.15, 0.75, 0.15, 0.85);
    _plwind (xmin, xmax, ymin, ymax);
    _plbox (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].gridx[p], ...
        PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].xtick.[p][1], ...
        PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].xtick.[p][2], ...
        PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].gridy[p], ...
        PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].ytick.[p][1], ...
        PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].ytick.[p][2]);

    // now do the plots
    for (idx_ll in range(_ll))
    {
      _ll1 = _ll[idx_ll];

      // silly user: ignore entries that are not properly formatted
      if (    !exist(data.[_ll1].x) || !exist(data.[_ll1].y) || !exist(data.[_ll1].z) ...
          ||  type(data.[_ll1].x)!="real" || type(data.[_ll1].y)!="real" || type(data.[_ll1].z)!="real"  )
      { continue; }

      if (data.[_ll1].x.n != data.[_ll1].z.nr || data.[_ll1].y.n != data.[_ll1].z.nc)
      { continue; }

      // Convert the data to log data if necessary.
      if (strindex (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].gridx[p], "l"))
      { x = log10 (real (data.[_ll1].x)); else x = real (data.[_ll1].x); }
      if (strindex (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].gridy[p], "l"))
      { y = log10 (real (data.[_ll1].y)); else y = real (data.[_ll1].y); }
      z = real (data.[_ll1].z);

      //
      // Draw the contours
      //
      llevel=[];
      for (i in range(clevel))
      {
        k = mod (i-1, _colors.nc) + 1;
        j = mod (i-1, 8) + 1;
        _pllsty(j);
        _plcol  (_colors[p;k]);
        if (_plcont (x, y, z, 1, x.n, 1, y.n, clevel[i]))
        {
          llevel = [llevel, clevel[i]];
        }
      }

      //
      // Reset color and draw the labels only at the end
      //
      if (idx_ll == _nos_data)
      {
        _plcol  (_colors[p;1]);
        _pllab (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].xlabel[p], PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].ylabel[p], PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].title[p]);

        //
        // Draw the contour legend. Use a new viewport to the right
        // of the contour plot.
        //
        //_plvpas (0.75, 1.0, 0.15, 0.85, PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].aspect[p]);
        _plvpor (0.75, 1.0, 0.15, 0.85);
        _plwind (0, 1, 0, 1);

        v = 1 - 1/(2*llevel.n);
        for (i in 1:llevel.n)
        {
          xl = [0.1, 0.2, 0.3]';
          yl = [v, v, v]';
          v = v - 1/llevel.n;

          k = mod (i-1, _colors.nc) + 1;
          j = mod (i-1, 8) + 1;

          _plcol  (_colors[p;k]);
          _pllsty (j);

          _plline (3, xl, yl);
          plptex (text(llevel[i],"%.2g"), xl[3]+.1, yl[3], , , 0);
        }
      }
    }
  }

  // Flush  and go back to text mode.
  _plflush ();
  _pltext ();

  //
  // Increment the plot no. so that next time
  // we use the correct settings.
  //

  PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].subplot = PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].subplot + 1;

  return PLPLOT_ACTIVE_WIN;
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Plot 3-D lines, etc...
//

pl3d = function ( X, Y, Z, BR )
{
  local (X, Y, Z, BR)
  check_plot_object ();

  p = mod (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].subplot, PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].nplot) + 1;	// The current index

  //
  // Some basic checks
  //

  if ((N = X.n) != Y.n) { error ("pl3d: X and Y must have same length"); }
  if (N != Z.n) { error ("pl3d: X and Z must have same length"); }

  if (!exist (BR)) { BR = N; }
  if (mod (N, BR) != 0) { error ("pl3d: X.n must be divisible by BR"); }
  iBR = int (N / BR);
  if (iBR == 1) { k = N; else k = BR; }

  //
  // Figure out the scale limits.
  // Needs improvement!
  //

  xmin = xmax = ymin = ymax = zmin = zmax = 0;
  </xmax;xmin;ymax;ymin;zmax;zmin/> = XYZ_scales (X, Y, Z, p);

  _plgra ();
  _plcol (1);
  _pllsty (PLPLOT_DEFAULT_LINE_STYLE);
  _plfont (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].font[p]);
  _plwid (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].width[p]);

  if (!subplot_f)
  {
    _pladv (0);        // Advance 1 subplot
    else
    subplot_f = 0;     // The user has set the subplot
  }

  if (strindex (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].grid3x[p], "l"))
  { X = log10 (real (X)); else X = real (X); }
  if (strindex (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].grid3y[p], "l"))
  { Y = log10 (real (Y)); else Y = real (Y); }
  if (strindex (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].grid3z[p], "l"))
  { Z = log10 (real (Z)); else Z = real (Z); }

  // PLPLOT_BASE_X = 2; PLPLOT_BASE_Y = 2; PLPLOT_BASE_Z = 4;
  xmin2d = -2.0; xmax2d = 2.0;
  ymin2d = -3.0; ymax2d = 5.0;

  if (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].aspect[p] != 0)
  {
    _plvasp (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].aspect[p]);
  else
    _plvsta ();
  }

  _plwind (xmin2d, xmax2d, ymin2d, ymax2d);
  _plw3d (PLPLOT_BASE_X, PLPLOT_BASE_Y, PLPLOT_BASE_Z, xmin, xmax, ymin, ymax, ...
          zmin, zmax, PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].alt[p], PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].az[p]);
  _plbox3 (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].grid3x[p], PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].xlabel[p], 0, 0, ...
           PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].grid3y[p], PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].ylabel[p], 0, 0, ...
           PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].grid3z[p], PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].zlabel[p], 0, 0);
  _plmtex ("t", 1.0, 0.5, 0.5, PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].title[p]);

  _plcol (2);
  for (i in 1:iBR)
  {
    j = [(i-1)*k+1:i*k];
    _plline3 (k, X[j], Y[j], Z[j]);
  }
  _plflush ();
  _pltext ();

  //
  // Increment the plot no. so that next time
  // we use the correct settings.
  //

  PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].subplot = PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].subplot + 1;

  return PLPLOT_ACTIVE_WIN;
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Plot 3-D Impulses...
//

pl3dimp = function ( X, Y, Z )
{
  check_plot_object ();

  p = mod (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].subplot, PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].nplot) + 1;	// The current index

  //
  // Some basic checks
  //

  if ((N = X.n) != Y.n) { error ("pl3d: X and Y must have same length"); }
  if (N != Z.n) { error ("pl3d: X and Z must have same length"); }

  //
  // Figure out the scale limits.
  // Needs improvement!
  //

  xmin = xmax = ymin = ymax = zmin = zmax = 0;
  </xmax;xmin;ymax;ymin;zmax;zmin/> = XYZ_scales (X, Y, Z, p);

  _plgra ();
  _plcol (1);
  _pllsty (PLPLOT_DEFAULT_LINE_STYLE);
  _plfont (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].font[p]);
  _plwid (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].width[p]);

  if (!subplot_f)
  {
    _pladv (0);        // Advance 1 subplot
    else
    subplot_f = 0;     // The user has set the subplot
  }

  if (strindex (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].grid3x[p], "l"))
  { X = log10 (real (X)); else X = real (X); }
  if (strindex (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].grid3y[p], "l"))
  { Y = log10 (real (Y)); else Y = real (Y); }
  if (strindex (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].grid3z[p], "l"))
  { Z = log10 (real (Z)); else Z = real (Z); }

  // PLPLOT_BASE_X = 2; PLPLOT_BASE_Y = 2; PLPLOT_BASE_Z = 4;
  xmin2d = -2.0; xmax2d = 2.0;
  ymin2d = -3.0; ymax2d = 5.0;

  if (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].aspect[p] != 0)
  {
    _plvasp (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].aspect[p]);
  else
    _plvsta ();
  }

  _plenv (xmin2d, xmax2d, ymin2d, ymax2d, 0, -2);
  _plw3d (PLPLOT_BASE_X, PLPLOT_BASE_Y, PLPLOT_BASE_Z, xmin, xmax, ymin, ymax, ...
          zmin, zmax, PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].alt[p], PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].az[p]);
  _plbox3 (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].grid3x[p], PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].xlabel[p], 0, 0, ...
           PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].grid3y[p], PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].ylabel[p], 0, 0, ...
           PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].grid3z[p], PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].zlabel[p], 0, 0);
  _plmtex ("t", 1.0, 0.5, 0.5, PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].title[p]);

  _plcol (2);
  for (i in 1:X.nr)
  {
    _plline3 (2, [X[i],X[i]], [Y[i],Y[i]], [0,Z[i]]);
  }
  _plflush ();
  _pltext ();

  //
  // Increment the plot no. so that next time
  // we use the correct settings.
  //

  PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].subplot = PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].subplot + 1;

  return PLPLOT_ACTIVE_WIN;
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Various support functions for the WIN list
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//
// Replot
//

replot = function ( )
{
  check_plot_object ();
  _replot ();
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Set the X-axis label
//
plxlabel = function ( xstr, scale )
{
  check_plot_object ();
  i = mod (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].subplot, PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].nplot) + 1;

  if (!exist (xstr))
  { xstr = ""; }
  PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].xlabel[i] = xstr;

  if (exist (scale))
  { PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].xlabel_scale.[i] = scale; }

  return PLPLOT_ACTIVE_WIN;
};

plx2label = function ( xstr, scale )
{
  check_plot_object ();
  i = mod (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].subplot, PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].nplot) + 1;

  if (!exist (xstr))
  { xstr = ""; }
  PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].alt_xlabel[i] = xstr;

  if (exist (scale))
  { PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].alt_xlabel_scale.[i] = scale; }

  return PLPLOT_ACTIVE_WIN;
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Set the Y-axis label
//
plylabel = function ( xstr, scale )
{
  check_plot_object ();
  i = mod (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].subplot, PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].nplot) + 1;

  if (!exist (xstr))
  { xstr = ""; }

  PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].ylabel[i] = xstr;

  if (exist (scale))
  { PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].ylabel_scale.[i] = scale; }

  return PLPLOT_ACTIVE_WIN;
};

ply2label = function ( xstr, scale )
{
  check_plot_object ();
  i = mod (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].subplot, PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].nplot) + 1;

  if (!exist (xstr))
  { xstr = ""; }
  PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].alt_ylabel[i] = xstr;

  if (exist (scale))
  { PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].alt_ylabel_scale.[i] = scale; }

  return PLPLOT_ACTIVE_WIN;
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Set the Z-axis label
//

plzlabel = function ( xstr, alt_xstr )
{
  check_plot_object ();
  i = mod (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].subplot, PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].nplot) + 1;

  if (!exist (xstr))
  { xstr = ""; }
  PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].zlabel[i] = xstr;

  if (!exist (alt_xstr))
  { alt_xstr = ""; }
  PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].alt_zlabel[i] = alt_xstr;

  return PLPLOT_ACTIVE_WIN;
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Set the plot-title
//

pltitle = function ( xstr, scale )
{
  check_plot_object ();
  if (!exist (xstr))
  { xstr = ""; }

  if (!exist (scale))
  {
    def_scale = _plchr();
  else
  {
    if (length(scale)==1)
    {
      def_scale = [0, scale];
    else if (length(scale)==2)
    {
      def_scale = scale;
    }}
  }}

  p = mod (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].subplot, PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].nplot) + 1; // The current index

  if (strlen(xstr) >= 0)
  {
    PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].title[p] = xstr;
    PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].title_scale.[p] = def_scale;
  }

  return 0;
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Set the scale limits.
//

plimits = function ( xmin, xmax, ymin, ymax, zmin, zmax )
{
  check_plot_object ();
  i = mod (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].subplot, PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].nplot) + 1;

  if (exist (xmin))
  {
    PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].xmin[i] = xmin;
  else
    PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].xmin[i] = 1j;
  }
  if (exist (xmax))
  {
    PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].xmax[i] = xmax;
  else
    PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].xmax[i] = 1j;
  }

  if (exist (ymin))
  {
    PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].ymin[i] = ymin;
  else
    PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].ymin[i] = 1j;
  }
  if (exist (ymax))
  {
    PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].ymax[i] = ymax;
  else
    PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].ymax[i] = 1j;
  }

  if (exist (zmin))
  {
    PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].zmin[i] = zmin;
  else
    PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].zmin[i] = 1j;
  }
  if (exist (zmax))
  {
    PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].zmax[i] = zmax;
  else
    PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].zmax[i] = 1j;
  }
};

plimits2 = function ( xmin, xmax, ymin, ymax )
{
  check_plot_object ();
  i = mod (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].subplot, PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].nplot) + 1;

  if (exist (xmin))
  {
    PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].alt_xmin[i] = xmin;
  else
    PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].alt_xmin[i] = 1j;
  }
  if (exist (xmax))
  {
    PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].alt_xmax[i] = xmax;
  else
    PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].alt_xmax[i] = 1j;
  }

  if (exist (ymin))
  {
    PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].alt_ymin[i] = ymin;
  else
    PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].alt_ymin[i] = 1j;
  }
  if (exist (ymax))
  {
    PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].alt_ymax[i] = ymax;
  else
    PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].alt_ymax[i] = 1j;
  }
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Set 2-D grid styles. A not-so-friendly interface.
//
plgrid = function ( sty_x, sty_y )
{
  check_plot_object ();
  i = mod (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].subplot, PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].nplot) + 1;

  if (exist (sty_x))
  {
    if (class (sty_x) == "num")
    {
      if (sty_x>=1)
      {
        PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].gridx[i] = ...
            PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].gridx[i] + "g";
      }
      if (sty_x==2)
      {
        PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].gridx[i] = ...
            PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].gridx[i] + "h";
      }
    }
    if (class (sty_x) == "string")
    {
      PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].gridx[i] = sty_x[1];
    }
  else
    PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].gridx[i] = PLPLOT_GRID_X_DEFAULT;
  }

  if (exist (sty_y))
  {
    if (class (sty_y) == "num")
    {
      if (sty_y>=1)
      {
        PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].gridy[i] = ...
            PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].gridy[i] + "g";
      }
      if (sty_y==2)
      {
        PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].gridy[i] = ...
            PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].gridy[i] + "h";
      }
    }
    if (class (sty_y) == "string")
    {
      PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].gridy[i] = sty_y[1];
    }
  else
    PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].gridx[i] = PLPLOT_GRID_X_DEFAULT;
  }
};

plgrid2 = function ( sty_x, sty_y )
{
  check_plot_object ();
  i = mod (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].subplot, PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].nplot) + 1;

  if (exist (sty_x))
  {
    if (class (sty_x) == "num")
    {
      if (sty_x>=1)
      {
        PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].alt_gridx[i] = ...
            PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].alt_gridx[i] + "g";
      }
      if (sty_x==2)
      {
        PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].alt_gridx[i] = ...
            PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].alt_gridx[i] + "h";
      }
    }
    if (class (sty_x) == "string")
    {
      PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].alt_gridx[i] = sty_x[1];
    }
  else
    PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].alt_gridx[i] = PLPLOT_ALT_GRID_X_DEFAULT;
  }

  if (exist (sty_y))
  {
    if (class (sty_y) == "num")
    {
      if (sty_y>=1)
      {
        PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].alt_gridy[i] = ...
            PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].alt_gridy[i] + "g";
      }
      if (sty_y==2)
      {
        PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].alt_gridy[i] = ...
            PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].alt_gridy[i] + "h";
      }
    }
    if (class (sty_y) == "string")
    {
      PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].alt_gridy[i] = sty_y[1];
    }
  else
    PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].alt_gridx[i] = PLPLOT_ALT_GRID_X_DEFAULT;
  }
};


plscale = function ( scale_x, scale_y )
{
  check_plot_object ();
  i = mod (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].subplot, PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].nplot) + 1;

  if (exist(scale_x))
  {
    if (class(scale_x)=="string")
    {
      if (strindex(tolower(scale_x),"log")>0)
      {
        PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].gridx[i] = PLPLOT_GRID_X_DEFAULT_LOG;
      else if (strindex(tolower(scale_x),"lin")>0)
      {
        PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].gridx[i] = PLPLOT_GRID_X_DEFAULT;
      else
        PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].gridx[i] = scale_x;
      }}
    }
  else
    PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].gridx[i] = PLPLOT_GRID_X_DEFAULT;
  }

  if (exist(scale_y))
  {
    if (class(scale_y)=="string")
    {
      if (strindex(tolower(scale_y),"log")>0)
      {
        PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].gridy[i] = PLPLOT_GRID_Y_DEFAULT_LOG;
      else if (strindex(tolower(scale_y),"lin")>0)
      {
        PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].gridy[i] = PLPLOT_GRID_Y_DEFAULT;
      else
        PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].gridy[i] = scale_y;
      }}
    }
  else
    PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].gridy[i] = PLPLOT_GRID_Y_DEFAULT;
  }

  return 0;
};

plscale2 = function ( scale_x, scale_y )
{
  check_plot_object ();
  i = mod (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].subplot, PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].nplot) + 1;

  PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].alt_gridx[i] = PLPLOT_ALT_GRID_X_DEFAULT;
  if (exist(scale_x))
  {
    if (strindex(tolower(scale_x),"log"))
    { PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].alt_gridx[i] = PLPLOT_ALT_GRID_X_DEFAULT_LOG; }
  }

  PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].alt_gridy[i] = PLPLOT_ALT_GRID_Y_DEFAULT;
  if (exist(scale_y))
  {
    if (strindex(tolower(scale_y),"log"))
    { PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].alt_gridy[i] = PLPLOT_ALT_GRID_Y_DEFAULT_LOG; }
  }

  return 0;
};


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Set 3-D grid (axis) styles
//
plgrid3 = function ( sty_x, sty_y, sty_z )
{
  check_plot_object ();
  i = mod (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].subplot, PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].nplot) + 1;
  if (exist (sty_x))
  {
    if (class (sty_x) == "string")
    {
      PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].grid3x[i] = sty_x;
    else
      error ("plgrid3: requires string argument GRID_STY_X");
    }
  else
    PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].grid3x[i] = PLPLOT_GRID3D_X_DEFAULT;
  }
  if (exist (sty_y))
  {
    if (class (sty_y) == "string")
    {
      PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].grid3y[i] = sty_y;
    else
      error ("plgrid3: requires string argument GRID_STY_Y");
    }
  else
    PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].grid3y[i] = PLPLOT_GRID3D_Y_DEFAULT;
  }
  if (exist (sty_z))
  {
    if (class (sty_z) == "string")
    {
      PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].grid3z[i] = sty_z;
    else
      error ("plgrid3: requires string argument GRID_STY_Z");
    }
  else
    PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].grid3z[i] = PLPLOT_GRID3D_Z_DEFAULT;
  }
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Various internal support functions. Eventually these will be static.
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

find_scales = function (ax, M, use_cols, p)
{
  //
  // 1st check for un-plottable data
  //
  if (M==[])
  { error ("plot: Can't plot empty matrix"); }
  if (any (any (isinf (M))))
  { error ("plot: cannot plot Infs"); }
  if (any (any (isnan (M))))
  { error ("plot: cannot plot NaNs"); }

  if (!exist(ax))
  { ax = "x"; }

  //
  // Check computed scale limits against user's
  //
  if (ax == "x")
  {
    // do the x-axis
    if (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].xmin[p] != 1j)
    {
      xmin = real (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].xmin[p]);
    else
      xmin = min(min(M[;use_cols]));
    }
    if (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].xmax[p] != 1j)
    {
      xmax = PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].xmax[p];
    else
      xmax = max(max(M[;use_cols]));
    }
  else
    // do the y-axis
    if (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].ymin[p] != 1j)
    {
      xmin = real (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].ymin[p]);
    else
      xmin = min(min(M[;use_cols]));
    }
    if (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].ymax[p] != 1j)
    {
      xmax = PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].ymax[p];
    else
      xmax = max(max(M[;use_cols]));
    }
  }

  //
  // Check for potential errors
  //
  if (xmin == xmax)
  {
    // As good a guess as any
    xmin = xmin - 1;
    xmax = xmax + 1;
  }

  return <<min=xmin;max=xmax>>;
};


//
// from data list, column cell-array and axes array figure out
//  plot limits for first and second set of coordinates
//
find_list_xy_scales = function (data, cols, axes, p)
{
  if (!length(data))
  { error ("plot: Can't plot empty list"); }

  m = members(data);

  xmin = [];
  xmax = [];
  idx_x1 = find((axes == PLPLOT_AXES_REV.x1y1)||(axes == PLPLOT_AXES_REV.x1y2));
  if (!isempty(idx_x1))
  {
    // do the x-axis
    if (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].xmin[p] != 1j)
    {
      xmin = real (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].xmin[p]);
    else
      for (i in range(idx_x1))
      {
        m_i = m[idx_x1[i]];
        xmin = min(xmin, min(data.[m_i][;cols[idx_x1[i]][1]]));
      }
    }

    // do the x-axis
    if (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].xmax[p] != 1j)
    {
      xmax = real (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].xmax[p]);
    else
      for (i in range(idx_x1))
      {
        m_i = m[idx_x1[i]];
        xmax = max(xmax, max(data.[m_i][;cols[idx_x1[i]][1]]));
      }
    }
  }

  ymin = [];
  ymax = [];
  idx_y1 = find((axes == PLPLOT_AXES_REV.x1y1)||(axes == PLPLOT_AXES_REV.x2y1));
  if (!isempty(idx_y1))
  {
    // do the y-axis
    if (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].ymin[p] != 1j)
    {
      ymin = real (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].ymin[p]);
    else
      for (i in range(idx_y1))
      {
        m_i = m[idx_y1[i]];
        ymin = min(ymin, min(data.[m_i][;cols[idx_y1[i]][2]]));
      }
    }

    // do the y-axis
    if (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].ymax[p] != 1j)
    {
      ymax = real (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].ymax[p]);
    else
      for (i in range(idx_y1))
      {
        m_i = m[idx_y1[i]];
        ymax = max(ymax, max(data.[m_i][;cols[idx_y1[i]][2]]));
      }
    }
  }

  alt_ymin = [];
  alt_ymax = [];
  idx_y2 = find((axes == PLPLOT_AXES_REV.x1y2)||(axes == PLPLOT_AXES_REV.x2y2));
  if (!isempty(idx_y2))
  {
    // do the y-axis
    if (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].alt_ymin[p] != 1j)
    {
      alt_ymin = real (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].alt_ymin[p]);
    else
      for (i in range(idx_y2))
      {
        m_i = m[idx_y2[i]];
        alt_ymin = min(alt_ymin, min(data.[m_i][;cols[idx_y2[i]][2]]));
      }
    }

    // do the y-axis
    if (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].alt_ymax[p] != 1j)
    {
      alt_ymax = real (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].alt_ymax[p]);
    else
      for (i in range(idx_y2))
      {
        m_i = m[idx_y2[i]];
        alt_ymax = max(alt_ymax, max(data.[m_i][;cols[idx_y2[i]][2]]));
      }
    }
  }

  alt_xmin = [];
  alt_xmax = [];
  idx_x2 = find((axes == PLPLOT_AXES_REV.x2y1)||(axes == PLPLOT_AXES_REV.x2y2));
  if (!isempty(idx_x2))
  {
    // do the x-axis
    if (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].alt_xmin[p] != 1j)
    {
      alt_xmin = real (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].alt_xmin[p]);
    else
      for (i in range(idx_x2))
      {
        m_i = m[idx_x2[i]];
        alt_xmin = min(alt_xmin, min(data.[m_i][;cols[idx_x2[i]][1]]));
      }
    }

    // do the x-axis
    if (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].alt_xmax[p] != 1j)
    {
      alt_xmax = real (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].alt_xmax[p]);
    else
      for (i in range(idx_x2))
      {
        m_i = m[idx_x2[i]];
        alt_xmax = max(alt_xmax, max(data.[m_i][;cols[idx_x2[i]][1]]));
      }
    }
  }

  return <<xmin=xmin;xmax=xmax;alt_xmin=alt_xmin;alt_xmax=alt_xmax; ...
      ymin=ymin;ymax=ymax;alt_ymin=alt_ymin;alt_ymax=alt_ymax>>;
};

find_xyz_scales = function (ax, data, p)
{
  //
  // 1st check for un-plottable data
  //
  if (class(data) != "list")
  { error ("plot: Horrible internal error! List expected"); }
  if (length(data)==0)
  { error ("plot: Can't plot empty list"); }

  xmax = [];
  xmin = [];

  _m = members(data);
  for (_im in _m)
  {
    if (class(data.[_im])=="num")
    {
      if (ax == "x")
      {
        xmax = max(xmax, max(data.[_im][;1]));
        xmin = min(xmin, min(data.[_im][;1]));
        continue;
      }
      if (ax == "y")
      {
        xmax = max(xmax, max(data.[_im][;2]));
        xmin = min(xmin, min(data.[_im][;2]));
        continue;
      }
      if (ax == "z")
      {
        xmax = max(xmax, max(data.[_im][;3]));
        xmin = min(xmin, min(data.[_im][;3]));
        continue;
      }
    else if (class(data.[_im])=="list")
    {
      xmax = max([xmax, max(data.[_im].[ax])]);
      xmin = min([xmin, min(data.[_im].[ax])]);
    else
      error ("plot: Can't find min/max values\n");
    }}
  }


  //
  // Check computed scale limits against user's
  //
  if (ax == "x")
  {
    // do the x-axis
    if (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].xmin[p] != 1j)
    { xmin = real (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].xmin[p]); }
    if (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].xmax[p] != 1j)
    { xmax = PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].xmax[p]; }
  }
  if (ax == "y")
  {
    // do the y-axis
    if (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].ymin[p] != 1j)
    { xmin = real (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].ymin[p]); }
    if (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].ymax[p] != 1j)
    { xmax = PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].ymax[p]; }
  }
  if (ax == "z")
  {
    // do the Z-axis
    if (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].zmin[p] != 1j)
    { xmin = real (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].zmin[p]); }
    if (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].zmax[p] != 1j)
    { xmax = PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].zmax[p]; }
  }

  //
  // Check for potential errors
  //
  if (xmin == xmax)
  {
    // As good a guess as any
    xmin = xmin - 1;
    xmax = xmax + 1;
  }

  return <<min=xmin;max=xmax>>;
};

//
// Find the X or Y scale limits .
// M can be a multi-column matrix, all columns
// will be used.
//
x_scales = function ( M, p)
{
  //
  // 1st check for un-plottable data
  //
  if (M==[])
  { error ("plot: Can't plot empty matrix"); }
  if (any (any (isinf (M))))
  { error ("plot: cannot plot Infs"); }
  if (any (any (isnan (M))))
  { error ("plot: cannot plot NaNs"); }

  //
  // Check computed scale limits against user's
  //
  if (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].xmin[p] != 1j)
  {
    xmin = real (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].xmin[p]);
  else
    xmin = min(min(M));
  }
  if (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].xmax[p] != 1j)
  {
    xmax = PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].xmax[p];
  else
    xmax = max(max(M));
  }

  //
  // Check for potential errors
  //
  if (xmin == xmax)
  {
    // As good a guess as any
    xmin = xmin - 1;
    xmax = xmax + 1;
  }

  //
  // Finally, adjust if log-scales
  //
  if (strindex (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].gridx[p], "l")>0)
  {
    if (xmin <= 0 || xmax <= 0)
    { error ("cannot plot log(x<=0)"); }
    xmin = log10 (xmin);
    xmax = log10 (xmax);
  }

  return <<xmin=xmin;xmax=xmax>>;
};

y_scales = function ( M, p )
{
  //
  // 1st check for un-plottable data
  //
  if (M==[])
  { error ("plot: Can't plot empty matrix"); }
  if (any (any (isinf (M))))
  { error ("plot: cannot plot Infs"); }
  if (any (any (isnan (M))))
  { error ("plot: cannot plot NaNs"); }

  //
  // Check computed scale limits against user's
  //
  if (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].xmin[p] != 1j)
  {
    xmin = real (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].xmin[p]);
  else
    xmin = min(min(M));
  }
  if (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].xmax[p] != 1j)
  {
    xmax = PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].xmax[p];
  else
    xmax = max(max(M));
  }

  //
  // Check for potential errors
  //
  if (xmin == xmax)
  {
    // As good a guess as any
    xmin = xmin - 1;
    xmax = xmax + 1;
  }

  //
  // Finally, adjust if log-scales
  //
  if (strindex (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].gridy[p], "l"))
  {
    if (xmin <= 0 || xmax <= 0) { error ("cannot plot log(y<=0)"); }
    xmin = log10 (xmin);
    xmax = log10 (xmax);
  }

  return << xmin = xmin ; xmax = xmax >>;
};

z_scales = function ( M, p, xmin, xmax )
{

  //
  // 1st check for un-plottable data
  //

  if (any (any (isinf (M))))
  { error ("plot: cannot plot Infs"); }
  if (any (any (isnan (M))))
  { error ("plot: cannot plot NaNs"); }

  xmin = min (min (M));
  xmax = max (max (M));

  //
  // Check computed scale limits against user's
  //

  if (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].zmin[p] != 1j) { xmin = real (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].zmin[p]); }
  if (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].zmax[p] != 1j) { xmax = real (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].zmax[p]); }

  //
  // Check for potential errors
  //

  if (xmin == xmax)
  {
    // As good a guess as any
    xmin = xmin - 1;
    xmax = xmax + 1;
  }

  //
  // Finally, adjust if log-scales
  //

  if (strindex (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].gridz[p], "l"))
  {
    if (xmin <= 0 || xmax <= 0) { error ("plot: cannot plot log(z<=0)"); }
    xmin = log10 (xmin);
    xmax = log10 (xmax);
  }

  return << xmin = xmin ; xmax = xmax >>;
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Find the X and Y scales for a single matrix. (OLD)
//

xy_scales = function ( M, p, xmin, xmax, ymin, ymax )
{
  //
  // 1st check for un-plottable data
  //

  if (any (any (isinf (M))))
  { error ("plot: cannot plot infs"); }
  if (any (any (isnan (M))))
  { error ("plot: cannot plot NaNs"); }

  if (M.nc == 1)
  {
    assign(xmin,1);
    assign(xmax, M.nr);
    assign(ymin, min (M));
    assign(ymax, max (M));
  else if (M.nc==3)
  {
    assign(xmin, min (M[;1]));
    assign(xmax, max (M[;1]));
    assign(ymin, min (min (M[;2]-abs(M[;3])) ));
    assign(ymax, max (max (M[;2]+abs(M[;3])) ));
  else if (M.nc==4)
  {
    assign(xmin, min (M[;1]));
    assign(xmax, max (M[;1]));
    assign(ymin, min (min (M[;3]) ));
    assign(ymax, max (max (M[;4]) ));
  else
    assign(xmin, min (M[;1]));
    assign(xmax, max (M[;1]));
    assign(ymin, min (min (M[;2:M.nc]) ));
    assign(ymax, max (max (M[;2:M.nc]) ));
  }}}

  //
  // Check computed scale limits against user's
  //

  if (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].xmin[p] != 1j)
  { assign(xmin, real (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].xmin[p])); }
  if (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].xmax[p] != 1j)
  { assign(xmax, real (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].xmax[p])); }
  if (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].ymin[p] != 1j)
  { assign(ymin, real (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].ymin[p])); }
  if (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].ymax[p] != 1j) { assign(ymax, real (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].ymax[p])); }

  //
  // Check for potential errors
  //

  if (xmin == xmax)
  {
    assign(xmin, xmin - 1);
    assign(xmax, xmax + 1);
  }

  if (ymin == ymax)
  {
    assign(ymin, ymin - 1);
    assign(ymax, ymax + 1);
  }

  //
  // Finally, adjust if log-scales
  //

  if (strindex (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].gridx[p], "l"))
  {
    if (xmin <= 0 || xmax <= 0) { error ("plot: cannot plot log <= 0"); }
    assign(xmin, log10 (xmin));
    assign(xmax, log10 (xmax));
  }
  if (strindex (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].gridy[p], "l"))
  {
    if (ymin <= 0 || ymax <= 0) { error ("plot: cannot plot log <= 0"); }
    assign(ymin, log10 (ymin));
    assign(ymax, log10 (ymax));
  }

  return 1;
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Find the X, Y and Z scales for a single matrix.
//

XYZ_scales = function ( X, Y, Z, p )
{
  // X - scale
  if (any (any (isinf (X))))
  { error ("cannot plot infs"); }
  if (any (any (isnan (X))))
  { error ("cannot plot NaNs"); }

  xmin = min (real (X));
  xmax = max (real (X));

  // Y - scale
  if (any (any (isinf (Y))))
  { error ("cannot plot infs"); }
  if (any (any (isnan (Y))))
  { error ("cannot plot NaNs"); }

  ymin = min (real (Y));
  ymax = max (real (Y));

  // Z - scale
  if (any (any (isinf (Y))))
  { error ("cannot plot infs"); }
  if (any (any (isnan (Y))))
  { error ("cannot plot NaNs"); }

  zmin = min (min (real (Z)));
  zmax = max (max (real (Z)));

  //
  // Check computed scale limits against user's
  //

  if (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].xmin[p] != 1j) { xmin = real (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].xmin[p]); }
  if (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].xmax[p] != 1j) { xmax = real (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].xmax[p]); }
  if (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].ymin[p] != 1j) { ymin = real (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].ymin[p]); }
  if (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].ymax[p] != 1j) { ymax = real (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].ymax[p]); }
  if (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].zmin[p] != 1j) { zmin = real (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].zmin[p]); }
  if (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].zmax[p] != 1j) { zmax = real (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].zmax[p]); }

  //
  // Check for potential errors
  //

  if (xmin == xmax)
  {
    // As good a guess as any
    xmin = xmin - 1;
    xmax = xmax + 1;
  }

  if (ymin == ymax)
  {
    // As good a guess as any
    ymin = ymin - 1;
    ymax = ymax + 1;
  }

  if (zmin == zmax)
  {
    // As good a guess as any
    zmin = zmin - 1;
    zmax = zmax + 1;
  }

  //
  // Finally, adjust if log-scales
  //

  if (strindex (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].grid3x[p], "l"))
  {
    if (xmin <= 0 || xmax <= 0) { error ("plot: cannot plot log(x<=0)"); }
    xmin = log10 (xmin);
    xmax = log10 (xmax);
  }

  if (strindex (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].grid3y[p], "l"))
  {
    if (ymin <= 0 || ymax <= 0) { error ("plot: cannot plot log(y<=0)"); }
    ymin = log10 (ymin);
    ymax = log10 (ymax);
  }

  if (strindex (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].grid3z[p], "l"))
  {
    if (zmin <= 0 || zmax <= 0) { error ("plot: cannot plot log(z<=0)"); }
    zmin = log10 (zmin);
    zmax = log10 (zmax);
  }

  return <<xmin=xmin; xmax=xmax; ymin=ymin; ymax=ymax; zmin=zmin; zmax=zmax>>;
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Find the X and Y scales for a list of matrices
//

list_scales = function ( data, key, p )
{
  once = 1;

  for (i in members (data))
  {
    if(class(data.[i]) == "list")
    {
      if(ishist(data.[i]))
      { M = hist_line(data.[i]); }
    else if (class(data.[i]) == "num")
    {
      M = real (data.[i]);
    else
      continue;
    }}

    if (M.n == 0) { continue; }

    //
    // 1st check for un-plottable data
    //

    if (any (any (isinf (M))))
    { error ("plot: cannot plot infs"); }
    if (any (any (isnan (M))))
    { error ("plot: cannot plot NaNs"); }

    k = find ((1:M.nc) != abs (key));
    if (key > 0)
    {
      if (M.nc != 1)
      {
	      </xmax;xmin/> = x_scales (real(M)[;key], p);
	      </ymax;ymin/> = y_scales (real(M)[;k],   p);
      else
	      </xmax;xmin/> = x_scales ((1:M.nr)', p);
	      </ymax;ymin/> = y_scales (real(M),   p);
      }
    else if (key < 0)
    {
      </xmax;xmin/> = x_scales ( real(M)[;k],        p);
      </ymax;ymin/> = y_scales ( real(M)[;abs(key)], p);
    else
      </xmax;xmin/> = x_scales ( (1:M.nr)', p);
      </ymax;ymin/> = y_scales ( real(M),   p);
    }}

    if (once)
    {
      Xmin = xmin; Xmax = xmax; Ymin = ymin; Ymax = ymax;
      once = 0;
    }
    if (xmin < Xmin) { Xmin = xmin; }
    if (xmax > Xmax) { Xmax = xmax; }
    if (ymin < Ymin) { Ymin = ymin; }
    if (ymax > Ymax) { Ymax = ymax; }
  }

  return <<xmin=Xmin; xmax=Xmax; ymin=Ymin; ymax=Ymax >>;
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Plot the columns of a matrix (core function)
//
// Notes: This is the core function for plotting a matrix. If the
// matrix is a single column, then the matrix elements are plotted
// versus the row numbers. If it is a multi-column matrix, then
// columns 2:N are plotted versus column 1.
//
// p, K, k and l are indices for plot features.
// p: the current plot index (the plot //)
// K: usually 0. This index is used to start of the line style and
// color index (k = color index, l = line-style index). This is mostly
// used by plot_list, which may call plot_matrix repeatedly.
// k: the line color index. This value determines the line color used
// for each column of data. If K = 0, then k goes like 2:14, then
// flops back to 1:14.
// l: the line style inex. This value determines the line style used
// for each column of data - not the line-type (points, or lines). If
// K = 0, then l goes like 2:8, then flops back to 1:8.
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

plot_matrix = function ( M, key, p, K, xmin, xmax, ymin, ymax, v )
{
  np = M.nr;

  if (M.nc == 1)
  {
    x = 1:M.nr;
    y = real (M);
    k = mod (1+K, 14) + 1;
    l = mod (1+K, 8);

    if (strindex (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].gridx[p], "l"))
    { x = log10 (x); }
    if (strindex (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].gridy[p], "l"))
    { y = log10 (y); }

    _plcol (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].color[p;k]);
    _pllsty (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].lstyle[p;l]);

    if (get_style (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].style.[p], k-1) == "line")
    {
      _plline (M.nr, x, y);
    else if (get_style (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].style.[p], k-1) == "point") {
      _plpoin (M.nr, x, y, PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].pstyle[p;l]);
    else if (get_style (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].style.[p], k-1) == "line-point") {
      _plline (M.nr, x, y);
      _plpoin (M.nr, x, y, PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].pstyle[p;l]);
    else {
      _plline (M.nr, x, y);
    }}}}

    //
    // Now do the legend
    //
    if (!any (any (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].desc.[p] == 1j)))
    {
      // Use the default if necessary
      if (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].desc.[p][1] == "default")
      {
        desc = "c1";
      else if (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].desc.[p].n >= k-1)
      {
        desc = PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].desc.[p][k-1];
      else
        // Not sure what to do, user has messed up.
        desc = "";
      }}

      v = v - (ymax-ymin)/11;
      xl = (xmax-xmin)*[10.5/12, 11/12, 11.5/12]' + xmin;
      yl = [v, v, v]' + ymin;

      if (get_style (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].style.[p], k-1) == "line")
      {
        _plline (3, xl, yl);
      else if (get_style (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].style.[p], k-1) == "point")
      {
        _plpoin (3, xl, yl, PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].pstyle[p;l]);
      else if (get_style (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].style.[p], k-1) == "line-point")
      {
        _plline (3, xl, yl);
        _plpoin (3, xl, yl, PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].pstyle[p;l]);
      }}}

      plptex(desc, xl[1]-(xmax-xmin)/25, yl[3], , , 1);

    }

  else

    //
    // Check for large column dimension
    //

    if (M.nc > 3*M.nr)
    {
      printf (" Plot %i columns and %i rows, are you sure (y/n) ? ", M.nc, M.nr);
      ans = getline ("stdin");
      if (ans.[1] != "y")
      { return -1; }
    }

    ki = find ((1:M.nc) != abs (key));
    for (i in ki)
    {
      if (key > 0)
      {
        x = real (M[;key]);
        y = real (M[;i]);
      else if (key < 0)
      {
        x = real (M[;i]);
        y = real (M[;abs(key)]);
      else
        x = (1:M.nr)';
        y = real (M[;i]);
      }}

      // Check for log scales, adjust if necessary
      if (strindex (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].gridx[p], "l"))
      { x = log10 (x); }
      if (strindex (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].gridy[p], "l"))
      { y = log10 (y); }

      k = mod (i-1 + K, 14) + 1;
      l = mod (8 + i-2 + K, 8) + 1;

      _plcol (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].color[p;k]);
      _pllsty (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].lstyle[p;l]);

      if (get_style (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].style.[p], k-1) == "line")
      {
        _plline (np, x, y);
      else if (get_style (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].style.[p], k-1) == "point")
      {
        _plpoin (np, x, y, PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].pstyle[p;l]);
      else if (get_style (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].style.[p], k-1) == "line-point")
      {
        _plline (np, x, y);
        _plpoin (np, x, y, PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].pstyle[p;l]);
      else
        _plline (np, x, y);
      }}}

      //
      // Now do the legend
      //

      if (!any (any (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].desc.[p] == 1j)))
      {
        // Use the default if necessary
        if (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].desc.[p][1] == "default")
        {
          desc = "c" + num2str (i);
        else if (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].desc.[p].n >= k-1)
        {
          desc = PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].desc.[p][k-1];
        else
          // Not sure what to do, user has messed up.
          desc = "";
        }}

        v = v - (ymax-ymin)/11;
        xl = (xmax-xmin)*[10.5/12, 11/12, 11.5/12]' + xmin;
        yl = [v, v, v]' + ymin;

        if (get_style (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].style.[p], k-1) == "line")
        {
          _plline (3, xl, yl);
        else if (get_style (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].style.[p], k-1) == "point")
        {
          _plpoin (3, xl, yl, PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].pstyle[p;l]);
        else if (get_style (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].style.[p], k-1) == "line-point")
        {
          _plline (3, xl, yl);
          _plpoin (3, xl, yl, PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].pstyle[p;l]);
        }}}

        plptex(desc, xl[1]-(xmax-xmin)/25, yl[3], , , 1);

      }
    }
  }

  return k-1;
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Plot all of the matrices in a list on the same plot
//

plot_list = function ( L, key, p, xmin, xmax, ymin, ymax )
{
  k = 0;
  v = ymax - ymin;

  //
  // Sort out the list members
  //
  sl = list_sort (L);

  // Plot the list members with numeric labels 1st.
  if (exist (sl.num))
  {
    for (i in sl.num)
    {
      M = L.[i];
      if (class(M)=="list")
      {
        if(ishist(M))
        { M = hist_line(M); }
      }
      if (class (M) != "num")
      { continue; }
      if ((k = plot_matrix (M, key, p, k, xmin, xmax, ymin, ymax, v)) < 0)
      { return k; }
      v = v - (ymax-ymin)/11;
    }
  }

  // Now plot the list members with string labels.
  if (exist (sl.char))
  {
    for (i in sl.char)
    {
      M = L.[i];
      if (class(M)=="list")
      {
        if(ishist(M))
        { M = hist_line(M); }
      }
      if (class (M) != "num")
      { continue; }
      if ((k = plot_matrix (M, key, p, k, xmin, xmax, ymin, ymax, v)) < 0)
      { return k; }
      v = v - (ymax-ymin)/11;
    }
  }
  return 1;
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Check the elements of LIST.
// LIST must contain elements `x', `y',
// and `z'
//

check_3d_list = function ( LIST )
{
  //
  // Check existence and types
  //

  if (class (LIST) != "list") {
    error ("plot3: argument must be a list");
  }
  if (!exist (LIST.x)) {
    error ("plot3: arg must contain `x' member");
  else if (class (LIST.x) != "num") {
    error ("plot3: x must be numeric");
  } }
  if (!exist (LIST.y)) {
    error ("plot3: arg must contain `y' member");
  else if (class (LIST.y) != "num") {
    error ("plot3: y must be numeric");
  } }
  if (!exist (LIST.z)) {
    error ("plot3: arg must contain `z' member");
  else if (class (LIST.z) != "num") {
    error ("plot3: z must be numeric");
  } }

  //
  // Check sizes
  //

  if (LIST.x.n != LIST.z.nr)
  {
    error ("plot3: x.n != z.nr");
  }

  if (LIST.y.n != LIST.z.nc)
  {
    error ("plot3: y.n != z.nc");
  }

};



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Create a legend in the current plot window
//
// if pobj.desc.[p] = inf()		no legend
// if pobj.desc.[p] = "default"		default ("c1", "c2", ...)
// if pobj.desc.[p] = "string"		use "string" as description
//

//
// Set the current plot tics
//
plxtics = function ( x, nx, scale )
{
  check_plot_object ();
  p = mod (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].subplot, PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].nplot) + 1; // The current index

  if (!exist (x))
  { x = 0; }
  x = x[1];
  if (!exist (nx))
  { nx = 0; }
  nx = nx[1];
  if (!exist (scale))
  { scale = 1; }

  PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].xtick.[p][1] = x[1];
  PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].xtick.[p][2] = nx;
  PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].xtick_scale.[p] = scale;
  return PLPLOT_ACTIVE_WIN;
};

plx2tics = function ( x, nx, scale )
{
  check_plot_object ();
  p = mod (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].subplot, PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].nplot) + 1; // The current index

  if (!exist (x))
  { x = 0; }
  x = x[1];
  if (!exist (nx))
  { nx = 0; }
  nx = nx[1];
  if (!exist (scale))
  { scale = 1; }

  PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].alt_xtick.[p][1] = x[1];
  PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].alt_xtick.[p][2] = nx;
  PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].alt_xtick_scale.[p] = scale;
  return PLPLOT_ACTIVE_WIN;
};

plytics = function ( x, nx, scale )
{
  check_plot_object ();
  p = mod (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].subplot, PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].nplot) + 1; // The current index

  if (!exist (x))
  { x = 0; }
  x = x[1];
  if (!exist (nx))
  { nx = 0; }
  nx = nx[1];
  if (!exist (scale))
  { scale = 1; }

  PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].ytick.[p][1] = x[1];
  PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].ytick.[p][2] = nx;
  PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].ytick_scale.[p] = scale;
  return PLPLOT_ACTIVE_WIN;
};

ply2tics = function ( x, nx, scale )
{
  check_plot_object ();
  p = mod (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].subplot, PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].nplot) + 1; // The current index

  if (!exist (x))
  { x = 0; }
  x = x[1];
  if (!exist (nx))
  { nx = 0; }
  nx = nx[1];
  if (!exist (scale))
  { scale = 1; }

  PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].alt_ytick.[p][1] = x[1];
  PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].alt_ytick.[p][2] = nx;
  PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].alt_ytick_scale.[p] = scale;
  return PLPLOT_ACTIVE_WIN;
};

plztics = function ( x, nx, scale )
{
  check_plot_object ();
  p = mod (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].subplot, PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].nplot) + 1; // The current index

  if (!exist (x))
  { x = 0; }
  x = x[1];
  if (!exist (nx))
  { nx = 0; }
  nx = nx[1];
  if (!exist (scale))
  { scale = 1; }

  PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].ztick.[p][1] = x[1];
  PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].ztick.[p][2] = nx;
  PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].ztick_scale.[p] = scale;
  return PLPLOT_ACTIVE_WIN;
};

//
// Set the current plot legend string
//
plegend = function ( LEGEND, scale, pos, pos_xy )
{
  check_plot_object ();

  p = mod (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].subplot, PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].nplot) + 1;	// The current index

  if (!exist (LEGEND))
  {
    PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].desc.[p] = 1j;
    return PLPLOT_ACTIVE_WIN;
  }

  if (!exist (scale))
  { scale = 1.0;}
  if (class(scale)!="num")
  { scale = 1.0; }

  if (!exist (pos))
  { pos = "itr"; }
  if (class(pos)!="string")
  { pos = "itr"; }

  if (!exist (pos_xy))
  { pos_xy = [0,0]; }
  if (class(pos_xy)!="num")
  { pos_xy = [0,0]; }
  if (length(pos_xy)!=2)
  { pos_xy = [0,0]; }

  if (class (LEGEND) == "string")
  {
    PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].desc.[p] = LEGEND;
    PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].desc_pos.[p] = pos;
    PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].desc_pos_xy.[p] = pos_xy;
    PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].desc_scale.[p] = scale;
  }

  return PLPLOT_ACTIVE_WIN;
};

set3d = function (bx, by, h)
{
  if (!exist (bx)) { PLPLOT_BASE_X = 2; else PLPLOT_BASE_X = bx; }
  if (!exist (by)) { PLPLOT_BASE_Y = 2; else PLPLOT_BASE_Y = by; }
  if (!exist (h)) { PLPLOT_BASE_Z = 4; else PLPLOT_BASE_Z = h; }
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Set the symbol size and scaling.
//

plssym = function ( DEF, SCALE )
{
  if (!exist (DEF)) { DEF = 0.0; }
  if (!exist (SCALE)) { SCALE = 1.0; }

  _plssym (DEF, SCALE);
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// fill polygon
//
plfill = function ( data )
{
  check_plot_object ();

  p = mod (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].subplot, PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].nplot) + 1;	// The current index
  key = 1;
  textf = 1;

  //
  // Draw the polygon
  // Step through the matrix plotting
  // each column versus the 1st
  //

  if (class (data) == "num")
  {
    //
    // Set up the plot basics
    //

    _plgra ();                       // Switch to graphics mode
    _plcol (1);                      // Set color
    _pllsty (PLPLOT_DEFAULT_LINE_STYLE);                     // Set line style
    _plpsty (1);                     // Set fill style
    _plfont (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].font[p]);       // Set font
    _plwid (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].width[p]);

    //
    // Either advance the plot (subplot), or let the user specify
    //

    if (!subplot_f)
    {
      _pladv (0);        // Advance 1 subplot
    else
      subplot_f = 0;     // The user has set the subplot
    }

    //
    // Set the aspect ratio
    //

    if (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].aspect[p] != 0)
    {
      _plvasp (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].aspect[p]);
    else
      _plvsta ();
    }

    //
    // Compute scale limits
    //

    k = find ((1:data.nc) != abs (key));
    if (key > 0)
    {
      if (data.nc != 1)
      {
	</xmax;xmin/> = x_scales (real(data)[;key], p);
	</ymax;ymin/> = y_scales (real(data)[;k],   p);
      else
	</xmax;xmin/> = x_scales ((1:data.nr)', p);
	</ymax;ymin/> = y_scales (real(data),   p);
      }
    else if (key < 0) {
      </xmax;xmin/> = x_scales (real(data)[;k],   p);
      </ymax;ymin/> = y_scales (real(data)[;abs(key)], p);
    else
      </xmax;xmin/> = x_scales ((1:data.nr)', p);
      </ymax;ymin/> = y_scales (real(data),   p);
    } }

    _plwind (xmin, xmax, ymin, ymax);
    _plbox (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].gridx[p], ...
        PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].xtick.[p][1], ...
        PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].xtick.[p][2], ...
        PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].gridy[p], ...
        PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].ytick.[p][1], ...
        PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].ytick.[p][2]);

    for (i in 2:data.nc)
    {
       _plfill(data.nr, data[;1], data[;i]);
       _plline(data.nr+1, [data[;1];data[1;1]], [data[;i];data[1;i]]);
    }

    else if (class (data) == "list") {

      _plgra ();
      _plcol (1);
      _pllsty (PLPLOT_DEFAULT_LINE_STYLE);
      _plpsty (1);
      _plfont (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].font[p]);
      _plwid (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].width[p]);

      </xmax;xmin;ymax;ymin/> = list_scales ( data, key, p );

      if (!subplot_f)
      {
	_pladv (0);        // Advance 1 subplot
      else
	subplot_f = 0;     // The user has set the subplot
      }

      if (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].aspect[p] != 0)
      {
	_plvasp (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].aspect[p]);
      else
        _plvsta ();
      }

      _plwind (xmin, xmax, ymin, ymax);
      _plbox (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].gridx[p], ...
        PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].xtick.[p][1], ...
        PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].xtick.[p][2], ...
        PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].gridy[p], ...
        PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].ytick.[p][1], ...
        PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].ytick.[p][2]);

      k = 0;
      for (j in members(data))
      {
          k = k + 1;
          for (i in 2:data.[j].nc)
          {
              _plcol (mod(k,14)+2);
              _pllsty (mod(k,8)+1);
              _plpsty (mod(k,8)+1);
              _plfill(data.[j].nr, data.[j][;1], data.[j][;i]);
              _pllsty(PLPLOT_DEFAULT_LINE_STYLE);
              _plline(data.[j].nr+1, [data.[j][;1];data.[j][1;1]], ...
                                     [data.[j][;i];data.[j][1;i]]);
          }
      }
    else
      error ("plot: un-acceptable argument");
  } }

  _pllsty (PLPLOT_DEFAULT_LINE_STYLE);
  _plcol (1);
  _pllab (PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].xlabel[p], ...
      PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].ylabel[p], ...
          PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].title[p]);
  _plflush ();
  if (textf)
  {
    _pltext ();
  }

  //
  // Increment the plot no. so that next time
  // we use the correct settings.
  //

  PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].subplot = PLPLOT_PLOT_WINDOWS.[PLPLOT_ACTIVE_WIN].subplot + 1;
  return PLPLOT_ACTIVE_WIN;
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// plot bar chart: plbar(y) or plbar([x,y])
//
plbar = function ( data )
{
  D = <<>>;
  if (min(data.nr,data.nc) == 1)
  {
     // only one row or column
     n = max(data.nr,data.nc);
     if (n == 1) { error("plbar: no data to plot"); }
     for (i in 1:n)
     {
         D.[i] = [i-0.45, 0; i+0.45, 0; i+0.45, data[i]; i-0.45, data[i]];
     }
  else
     // two columns
     if (data.nc > 2) { error("plbar: more thab two columns"); }
     data = data[sort(data[;1]).idx;];
     nr = data.nr;
     for (i in 1:nr-1)
     {
        dx = data[i+1;1] - data[i;1];
        D.[i] = [data[i;1]-0.45*dx, 0;
                 data[i;1]+0.45*dx, 0;
                 data[i;1]+0.45*dx,data[i;2];
                 data[i;1]-0.45*dx,data[i;2]];
     }
     D.[nr]   = [data[nr;1]-0.45*dx, 0;
                 data[nr;1]+0.45*dx, 0;
                 data[nr;1]+0.45*dx,data[nr;2];
                 data[nr;1]-0.45*dx,data[nr;2]];
  }
  P = plfill(D);
  return PLPLOT_ACTIVE_WIN;
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Get any user-specific startup information.
//
// The user must have a function called:
//
//     plot_startup_stuff
//
// in their RLAB2_RC0 (.rlab) file. This next statement will
// execute that function, performing any specific operations
// the user has included.
//

if (exist (plot_startup_stuff))
{
  plot_startup_stuff ();
}

