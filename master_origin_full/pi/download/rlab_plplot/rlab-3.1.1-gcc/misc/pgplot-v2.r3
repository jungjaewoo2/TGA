//
// New plot.r for use with PGPLOT library.
// The help files for these functions are in
// misc/plhelp
//

// pgplot-v2.r

// This file is a part of RLaB ("Our"-LaB) + rlabplus extension
// RLaB (C) 1997  Ian R. Searle
// rlabplus (C) 2003-2017, Marijan Kostrun

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

require libsystem

static(THIS_LIBRARY);
if (!exist(THIS_LIBRARY))
{ THIS_LIBRARY = "pgplot"; }

static(THIS_LIBRARY_DEFAULT_DEVICE);
if (!exist(THIS_LIBRARY_DEFAULT_DEVICE))
{ THIS_LIBRARY_DEFAULT_DEVICE = "/XWIN"; }

static (PLOT_WINDOWS)   // The static plot window structure
if (!exist (PLOT_WINDOWS))
{ PLOT_WINDOWS = <<>>; }

static (PLOT_ACTIVE_WIN);   // The active/current plot window
if (!exist(PLOT_ACTIVE_WIN))
{ PLOT_ACTIVE_WIN = 0; }

static(THIS_LIBRARY_DEFAULT_CHAR_SIZE);
if (!exist(THIS_LIBRARY_DEFAULT_CHAR_SIZE))
{ THIS_LIBRARY_DEFAULT_CHAR_SIZE = 1.5; }

//
// PLOT CONSTANTS
//
static (PLOT_COLORS, PLOT_LINE_STYLE, PLOT_POINT_STYLE, PLOT_FILL_STYLE);
PLOT_COLORS = 0:15;
PLOT_LINE_STYLE = 1:8;
PLOT_POINT_STYLE = 1:8;
PLOT_FILL_STYLE = 1:8;
static (PLOT_RESET_AFTER_PLOT);
PLOT_RESET_AFTER_PLOT = 1;
static(PLOT_DATASET_PLOT_FORMAT);
PLOT_DATASET_PLOT_FORMAT = [0, 1, 1, 1, 0, 1, 0, 0, 0, 0];
// Maintain the transformations for 3-D plots.
static (PLOT_BASE_X, PLOT_BASE_Y, PLOT_BASE_Z);
PLOT_BASE_X = 2;
PLOT_BASE_Y = 2;
PLOT_BASE_Z = 4;
// Defaults
//  2D grids
static (PLOT_GRID_X_DEFAULT, PLOT_ALT_GRID_X_DEFAULT);
PLOT_GRID_X_DEFAULT = "bcnst";
PLOT_ALT_GRID_X_DEFAULT = "cmst";
static(PLOT_GRID_X_DEFAULT_LOG,PLOT_ALT_GRID_X_DEFAULT_LOG);
PLOT_GRID_X_DEFAULT_LOG = "bclnst";
PLOT_ALT_GRID_X_DEFAULT_LOG = "lcmst";
static (PLOT_GRID_Y_DEFAULT,PLOT_ALT_GRID_Y_DEFAULT);
PLOT_GRID_Y_DEFAULT = "bnstv"; // or "bcnstv"
PLOT_ALT_GRID_Y_DEFAULT = "cmstv";
static(PLOT_GRID_Y_DEFAULT_LOG,PLOT_ALT_GRID_Y_DEFAULT_LOG);
PLOT_GRID_Y_DEFAULT_LOG = "lbcnstv";
PLOT_ALT_GRID_Y_DEFAULT_LOG = "lcmstv";
// 3D grids
static (PLOT_GRID3D_X_DEFAULT, PLOT_GRID3D_Y_DEFAULT, PLOT_GRID3D_Z_DEFAULT);
PLOT_GRID3D_X_DEFAULT = "bnstu";
PLOT_GRID3D_Y_DEFAULT = "bnstu";
PLOT_GRID3D_Z_DEFAULT = "bcdmnstuv";
// Initial windows size/config...
static (PLOT_WIN_LEN_X_DEFAULT, PLOT_WIN_LEN_Y_DEFAULT);
PLOT_WIN_LEN_X_DEFAULT = 7;
PLOT_WIN_LEN_Y_DEFAULT = 5;
static (PLOT_WIN_OFFSET_X_DEFAULT, PLOT_WIN_OFFSET_Y_DEFAULT);
PLOT_WIN_OFFSET_X_DEFAULT = 0;
PLOT_WIN_OFFSET_Y_DEFAULT = 0;
static (PLOT_WIN_XP_DEFAULT, PLOT_WIN_YP_DEFAULT);
PLOT_WIN_XP_DEFAULT = 0;
PLOT_WIN_YP_DEFAULT = 0;

static(PLOT_COLORMAP);
if (!exist(PLOT_COLORMAP))
{
  PLOT_COLORMAP = <<>>;

  //
  // default PGPLOT color scheme
  //
  PLOT_COLORMAP.PG = <<>>;
  PLOT_COLORMAP.PG.COLOR_BKG = "black";
  PLOT_COLORMAP.PG.COLOR = [ ...
      "white", "red", "green", "blue", "cyan", ...
      "magenta", "yellow", "orange", "greenyellow", "greencyan", ...
          "bluecyan", "bluemagenta", "redmagenta", "darkgray", "lightgray"];
  PLOT_COLORMAP.PG.RGB_BKG = [0,0,0];
  PLOT_COLORMAP.PG.RGB = [ ...
      [255, 255, 255];  [255, 0, 0];    [0, 255, 0];    [0, 0, 255];    [0, 255, 255];...
      [255, 0, 255];    [255, 255, 0];  [255, 128, 0];  [128, 255, 0];  [0, 255, 128]; ...
      [0, 128, 255];    [128, 0, 255];  [255, 0, 128];  [83, 83, 83];   [166, 166, 166] ];

  //
  // color scheme after xmgrace, and gnuplot
  //
  PLOT_COLORMAP.XMG.COLOR_BKG = "white";
  PLOT_COLORMAP.XMG.COLOR = [ ...
      "black", "red", "green", "blue", "yellow", ...
      "brown", "grey", "violet", "cyan", "magenta", ...
          "orange", "indigo", "maroon", "turquoise", "green5" ];
  PLOT_COLORMAP.XMG.RGB_BKG = [255, 255, 255];
  PLOT_COLORMAP.XMG.RGB = [ ...
      [0, 0, 0];        [255, 0, 0];      [0, 255, 0];    [0, 0, 255];    [255, 255, 0];...
      [188, 143, 143];  [153, 153, 153];  [148, 0, 211];  [0, 255, 255];  [255, 0, 255]; ...
      [255, 165, 0];    [114, 33, 188];   [103, 7, 72];   [64, 224, 208]; [0, 59, 0] ];

}

//
// load X11 palette from the configuration file
//
static(PLOT_X11_PALETTE_FILE, PLOT_X11_PALETTE);
if (!exist(PLOT_X11_PALETTE_FILE))
{
  PLOT_X11_PALETTE_FILE = "/usr/share/X11/rgb.txt";
  if (!isfile(PLOT_X11_PALETTE_FILE))
  {
    printf("File '%s' is provided as default X11 palette but it doesn't exist. Please check manually!\n");
  }

  if (!exist(PLOT_X11_PALETTE))
  {
    PLOT_X11_PALETTE = <<>>;
    if (isfile(PLOT_X11_PALETTE_FILE))
    {
      __s = reads(PLOT_X11_PALETTE_FILE);
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
              __name = tolower(sum(__y[4:__l],""));
            }
            else
            {
              __name = __y[4];
            }
            if (!exist(PLOT_X11_PALETTE.[__name]))
            { PLOT_X11_PALETTE.[__name] = strtod(__y[1:3]); }
          }
        }
      }
      clear(__s, __i, __j, __l);
    }
  }
}

//
// decide on default color scheme
//
static(PLOT_COLORMAP_SCHEME);
// PLOT_COLOR_SCHEME = "PG"; //pgplot's white on black
PLOT_COLORMAP_SCHEME = "XMG"; // xmgrace/gnuplot's black on white

static(PLOT_COLOR_RGB_NAME,PLOT_COLOR_RGB_NAME_BKG,PLOT_COLOR_RGB,PLOT_COLOR_RGB_BKG);
PLOT_COLOR_RGB_NAME = PLOT_COLORMAP.[PLOT_COLORMAP_SCHEME].COLOR;
PLOT_COLOR_RGB_NAME_BKG = PLOT_COLORMAP.[PLOT_COLORMAP_SCHEME].COLOR_BKG;
PLOT_COLOR_RGB = PLOT_COLORMAP.[PLOT_COLORMAP_SCHEME].RGB;
PLOT_COLOR_RGB_BKG = PLOT_COLORMAP.[PLOT_COLORMAP_SCHEME].RGB_BKG;

//
// plot management
//
static(PLOT_AXES_ID, PLOT_AXES_DESC, PLOT_AXES_REV);
PLOT_AXES_ID    = [     1,      2,      3,      4];
PLOT_AXES_DESC  = ["x1y1", "x1y2", "x2y1", "x2y2"];
PLOT_AXES_REV   = <<x1y1=1; x1y2=2; x2y1=3; x2y2=4>>;

static(PLOT_DEFAULT_LINE_STYLE);
PLOT_DEFAULT_LINE_STYLE = 1;

//
// Static (private) functions. For use from within
// this file only.
//
static (newplwin);           // kmk 2005, instead of P=P+1; use P=newplwin();
static (create_plot_object);
static (check_plot_object);
static (x_scales);
static (y_scales);
static (z_scales);
static (XYZ_scales);
static (list_scales);
//static (ticks)
static (list_sort);
static (hist_scales);
static (plot_matrix);
static (plot_list);
static (find_char);
static (get_style);
static (make_legend);

static (subplot_spec)          // Pay attention to the specifications ?
static (subplot_horizontal)    // The horizontal specification.
static (subplot_vertical)      // The vertical specification.

subplot_spec = 0;

static (axis3d)
axis3d = [ 1, 0, 0;
           0, 1, 0;
           0, 0, 1];

static(is3dlist);
is3dlist = function ( list )
{
  if (class(list)!="list")
  { return 0; }

  for (i in ["x", "y", "z"])
  {
    if (!exist(list.[i]))
    { return 0; }

    if (class(list.[i])!="num")
    { return 0; }
  }

  if (length(list.x) != size(list.z)[1])
  { return 0; }

  if (length(list.y) != size(list.z)[2])
  { return 0; }

  return 1;
};

//
// Static (private) functions. For use from within
// this file only.
//
static(convert_gnuformat_to_plot_format);
convert_gnuformat_to_plot_format = function(s)
{
  this_function = "convert_gnuformat_to_plot_format";

  if(class(s)!="string")
  { return PLOT_DATASET_PLOT_FORMAT; }

  tokens = tolower(strsplt(s, <<csp="'BLANK";lstrip="'BLANK">>));
  lentok = length(tokens);

  //
  // contour map
  //
  _ip = find(strindex(tokens,"cont")>0);
  if (length(_ip)==1)
  {
    // identify statement:
    //    labels 'x1','x2','x3',...
    //      for sequential labeling contour levels
    //    ts nn
    //      for text size to be used with printing the labels
    //    spa n1,n2
    //      for spacing between the labels, and minimum spacing to have a label on the graph
    _ilb = find(strindex(tokens,"lab")>0);
    if (length(_ilb)==1)
    {
      labels = "%g";  // default format for printing the labels
      ts = THIS_LIBRARY_DEFAULT_CHAR_SIZE;  // default text size for printing labels
      sp = [150,80];  // default values for intval,minint
      if (_ilb < lentok)
      {
        _c = strsplt(tokens[_ilb+1],",");
        _c = gsub(" ","_", _c).string;
        _d = gsub("'",_c);
        _i = find (_d.count==2);
        if (length(_i)>0)
        {
          labels = _d.string[_i];
        }
      }
      _its = find(strindex(tokens,"ts")>0);
      if (length(_its)==1)
      {
        if (_its<lentok)
        {
          ts = strtod(tokens[_its+1]);
          if (isnan(ts))
          { ts = THIS_LIBRARY_DEFAULT_CHAR_SIZE * ts; }
        }
      }
      _isp = find(strindex(tokens,"spa")>0);
      if (length(_isp)==1)
      {
        if (_isp<lentok)
        {
          sp = strtod(tokens[_isp+1],<<csp=",">>);
          if (length(sp)==2)
          {
            if (isnan(sp[1]))
            { sp[1] = 150; }
            if (isnan(sp[2]))
            { sp[2] = 80; }
          }
        }
      }
    }

    // contour
    lpb = 4;
    lt = 1;
    lw = 1;
    lc = 1;
    zmin = nan();
    zmax = nan();

    // identify cbrange [zmin:zmax]
    _icb = find(strindex(tokens,"cbra")>0);
    if ((length(_icb)==1) && (max(_icb) < lentok))
    {
      if (strindex(tokens[_icb+1],"[")==1)
      {
        if (strindex(tokens[_icb+1],"]")==strlen(tokens[_icb+1]))
        {
          if (strindex(tokens[_icb+1],":") > 2)
          {
            _c = gsub(["[","]"],tokens[_icb+1]).string;
            _c = strtod(_c,<<csp=":">>);
            if (all(!isnan(_c)))
            {
              zmin = min(_c);
              zmax = max(_c);
            }
          }
        }
      }
    }

    // figure out the levels at which the contours will be plotted
    //    levels NN
    //    levels incr x_i,dx,x_f
    zlevels = [];
    nzlev = 0;
    _izl = find(strindex(tokens,"lev")>0);
    if ((length(_izl)==1) && (max(_izl) < lentok))
    {
      if (tokens[_izl+1] == "incr")
      {
        if (_izl+2 <= lentok)
        {
          lv = strtod(tokens[_izl+2]);
          if (length(lv)==3)
          {
            if (all(isnan(lv))==0)
            {
              zlevels = [lv[1]:lv[3]:lv[2]];
              nzlev = length(zlevels);
            }
          }
        }
      }
      else if (tokens[_izl+1] == "disc")
      {
        if (_izl+2 <= lentok)
        {
          lv = strtod(tokens[_izl+2]);
          if (all(isnan(lv))==0)
          {
            zlevels = lv;
            nzlev = length(zlevels);
          }
        }
      }
      else
      {
        lv = strtod(tokens[_izl+1]);
        if (length(lv)==1)
        {
          if (isnumber(lv))
          {
            nzlev = lv;
          }
        }
      }
    }

    // figure out the line properties to be used for plotting the contour
    // levels
    //
    _ilt = find(tokens == "lt");
    if ((length(_ilt)==1) && max(_ilt) < length(tokens))
    {
      lt = strtod(tokens[_ilt+1]);
      while (lt < 1)
      { lt = lt + 8; }
      while (lt > 8)
      { lt = lt - 8; }
    }
    //
    lw = PLOT_DATASET_PLOT_FORMAT[4];
    _ilw = find(tokens == "lw");
    if ((length(_ilw)==1) && max(_ilw) < length(tokens))
    { lw = strtod(tokens[_ilw+1]); }
    if (lw < 1)
    { lw = 1; }
    if (lw > 201)
    { lw = 201; }

    // [CONT_TYPE, LT, LW, LC1, LC2, LC3]
    pm = [nan(),lt,lw,nan(),nan(),nan()];

    _ilc = find(tokens == "lc");
    if ((length(_ilc)==1) && max(_ilc) < lentok)
    {
      // single color contour lines;
      //  lc N
      //  lc rgb name
      if (tokens[_ilc+1] == "map")
      {
        // rgb map line type: CONT_TYPE=2, gnuplot default
        // format:
        //    lc map r,g,b
        pm[1] = 3;
        pm[4] = 7;
        pm[5] = 5;
        pm[6] = 16;
        if (_ilc + 2 <= lentok)
        {
          lc = strtod(tokens[_ilc+2],<<csp=",">>);
          if (length(lc)==3)
          { pm[4:6] = lc; }
        }
      }
      else if ((tokens[_ilc+1] == "grey")||(tokens[_ilc+1] == "gray"))
      {
        // grey line type: CONT_TYPE=3, gnuplot default
        // format:
        //    lc grey idx0,idx1
        pm[1] = 2;
        pm[4] = 0;  // bg - white
        pm[5] = 1;  // fg - black
        if (_ilc + 2 <= lentok)
        {
          lc = strtod(tokens[_ilc+2],<<csp=",">>);
          if (length(lc)==2)
          { pm[4:5] = lc; }
        }
      }
      else
      {
        // single line type: CONT_TYPE=1
        pm[1] = 1;
        lc = 1;
        if (tokens[_ilc+1] == "rgb")
        {
          c = gsub("'",tokens[_ilc+2]).string;
          _idx_c = find(PLOT_COLOR_RGB_NAME==c);
          if (length(_idx_c)==1)
          {
            lc = _idx_c;
          }
          else
          {
            printf(THIS_LIBRARY + ": Warning! Requested color '%s' not available in the palette!\n");
            lc = 1; // unknown color
          }
        }
        else
        {
          lc = strtod(tokens[_ilc+1]);
          while(lc<1)
          { lc = lc + length(PLOT_COLOR_RGB_NAME); }
          while (lc > length(PLOT_COLOR_RGB_NAME))
          { lc = lc - length(PLOT_COLOR_RGB_NAME); }
        }

        // ser color
        pm[4] = lc;
      }
    }

    //
    // is user using key word 'axes'?
    //
    axes = PLOT_AXES_REV.x1y1;  // default. in case user does not use it!
    _ia = find(tokens == "ax" || tokens == "axes");
    if ((length(_ia)==1) && max(_ia) < length(tokens))
    {
      // axis enforces use_cols=[1,2] in the absence of 'using 1:2...' in plformat
      if (exist(PLOT_AXES_REV.[tokens[_ia+1]]))
      {
        axes = PLOT_AXES_REV.[tokens[_ia+1]];
      }
      else
      {
        printf(THIS_LIBRARY + ": " + this_function + ": directive 'axes %s' is not supported!\n",  tokens[_ia+1]);
        printf(THIS_LIBRARY + ": " + this_function + ": Why don't you go away and read the manual for a change!\n");
        error (THIS_LIBRARY + ": " + this_function + ": Cannot continue!\n");
      }
    }

    // finish gnuplot-style processing
    if ((nzlev==0) && isempty(zlevels))
    { nzlev = 10;}

    rval = <<>>;
    rval.axes = axes;
    rval.use_cols = [];
    rval.plformat = [lpb, zmin, zmax, nzlev, pm];
    if (!isempty(zlevels))
    { rval.contour_zlevels = zlevels; }
    if (exist(labels))
    {
      rval.contour_labels = labels;
      if (exist(ts))
      { rval.contour_ts = ts; }
      if (exist(sp))
      { rval.contour_sp = sp; }
    }
    return rval;
  }

  //
  // palette mapping in 3d: grey,map,grad
  //
  _ip = find(tokens == "pm3d");
  if ((length(_ip)==1) && (max(_ip) < length(tokens)))
  {
    //
    pm = [0,0,1,0,0,0,0]; // default is the grey map background (0) to foreground (1) color

    //
    lpb = 3; // for pm3d
    zmin = nan();
    zmax = nan();

    // identify cbrange [zmin:zmax]
    _icb = find(strindex(tokens,"cbra")>0);
    if ((length(_icb)==1) && (max(_icb) < lentok))
    {
      if (strindex(tokens[_icb+1],"[")==1)
      {
        if (strindex(tokens[_icb+1],"]")==strlen(tokens[_icb+1]))
        {
          if (strindex(tokens[_icb+1],":") > 2)
          {
            _c = gsub(["[","]"],tokens[_icb+1]).string;
            _c = strtod(_c,<<csp=":">>);
            if (all(!isnan(_c)))
            {
              zmin = min(_c);
              zmax = max(_c);
            }
          }
        }
      }
    }

    // grey:
    if (strindex(tokens[_ip+1],"grey") || strindex(tokens[_ip+1],"grey"))
    {
      pm[1] = 0;  // grey
      pm[2] = 0;  // white
      pm[3] = 1;  // black
      if (length(tokens)>=_ip+2)
      {
        if (strindex(tokens[_ip+2],","))
        {
          _d = strtod(tokens[_ip+2]);
          if(length(_d)>0)
          { pm[2] = _d[1]; }
          if(length(_d)>1)
          { pm[3] = _d[2]; }
        }
      }
    }
    // map
    if (strindex(tokens[_ip+1],"map"))
    {
      pm[1] = 1;  // palette map
      pm[2] = 7;  // R gnuplot function
      pm[3] = 5;  // G function
      pm[4] = 15;  // B function
      if (length(tokens)>=_ip+2)
      {
        if (strindex(tokens[_ip+2],","))
        {
          _d = strtod(tokens[_ip+2]);
          if(length(_d)>0)
          { pm[2] = _d[1]; }
          if(length(_d)>1)
          { pm[3] = _d[2]; }
          if(length(_d)>2)
          { pm[4] = _d[3]; }
        }
      }
    }
    // gradient
    if (strindex(tokens[_ip+1],"grad"))
    {
      pm[1] = 2;  // gradient
      pm[2] = 2;  // number of colors in the table
      pm[3] = 1;  // index of first color
      if (length(tokens)>=_ip+2)
      {
        if (strindex(tokens[_ip+2],","))
        {
          _d = strtod(tokens[_ip+2]);
          if(length(_d)>0)
          { pm[2] = _d[1]; }
          if(length(_d)>1)
          { pm[3] = _d[2]; }
        }
      }
    }

    //
    // is user using key word 'axes'?
    //
    axes = PLOT_AXES_REV.x1y1;  // default. in case user does not use it!
    _ia = find(tokens == "ax" || tokens == "axes");
    if ((length(_ia)==1) && max(_ia) < length(tokens))
    {
      // axis enforces use_cols=[1,2] in the absence of 'using 1:2...' in plformat
      if (exist(PLOT_AXES_REV.[tokens[_ia+1]]))
      {
        axes = PLOT_AXES_REV.[tokens[_ia+1]];
      }
      else
      {
        printf(THIS_LIBRARY + ": " + this_function + ": directive 'axes %s' is not supported!\n",  tokens[_ia+1]);
        printf(THIS_LIBRARY + ": " + this_function + ": Why don't you go away and read the manual for a change!\n");
        error (THIS_LIBRARY + ": " + this_function + ": Cannot continue!\n");
      }
    }

    // finish gnuplot-style processing
    rval = <<>>;
    rval.axes = axes;
    rval.use_cols = [];
    rval.plformat = [lpb, zmin, zmax, pm];
    return rval;
  }

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
    }
    else if (strindex(tolower(tokens[_iw+1]), "b"))
    {
      lpb = 2;
    }

    if (strindex(tolower(tokens[_iw+1]), "err"))
    {
      err = 1;
      err_xy = 3;
      err_lp = 0;
      if (strindex(tolower(tokens[_iw+1]), "xy"))
      {
        err_xy=2;
      }
      else if (strindex(tolower(tokens[_iw+1]), "x"))
      {
        err_xy=1;
      }

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
  //    linecolor 1:15
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
    }
    else
    {
      lt = -1;
    }
    //
    _ilc = find(tokens == "lc");
    if ((length(_ilc)==1) && _ilc < length(tokens))
    {
      if (tokens[_ilc+1] == "rgb")
      {
        c = gsub("'",tokens[_ilc+2]).string;
        _idx_c = find(PLOT_COLOR_RGB_NAME==c);
        if (length(_idx_c)==1)
        {
          lc = _idx_c;
        }
        else
        {
          printf(THIS_LIBRARY + ": Warning! Requested color '%s' not available in the palette!\n");
          lc = -1; // unknown color
        }
      }
      else
      {
        lc = strtod(tokens[_ilc+1]);
        while(lc<1)
        { lc = lc + length(PLOT_COLOR_RGB_NAME); }
        while (lc > length(PLOT_COLOR_RGB_NAME))
        { lc = lc - length(PLOT_COLOR_RGB_NAME); }
      }
    }
    else
    {
      lc = -1;
    }

    //
    lw = PLOT_DATASET_PLOT_FORMAT[4];
    _ilw = find(tokens == "lw");
    if ((length(_ilw)==1) && _ilw < length(tokens))
    { lw = strtod(tokens[_ilw+1]); }

    if ((lc==-1)&&(lt!=-1))
    { lc = PLOT_DATASET_PLOT_FORMAT[3]; }
    if ((lt==-1)&&(lc!=-1))
    { lt = PLOT_DATASET_PLOT_FORMAT[2]; }
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
    }
    else
    {
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
        _idx_c = find(PLOT_COLOR_RGB_NAME==c);
        if (length(_idx_c)==1)
        {
          pc = _idx_c;
        }
        else
        {
          printf(THIS_LIBRARY + ": Warning! Requested color '%s' is not available in palette !\n", c);
          printf(THIS_LIBRARY + ": Warning! Available colors are %s. Please choose another!\n", ...
              sum(PLOT_COLOR_RGB_NAME,","));
          pc = -1; // unknown color
        }
      }
      else
      {
        pc = strtod(tokens[_ipc+1]);
        while(pc<1)
        { pc = pc + length(PLOT_COLOR_RGB_NAME); }
        while (pc > length(PLOT_COLOR_RGB_NAME))
        { pc = pc - length(PLOT_COLOR_RGB_NAME); }
      }
    }
    else
    {
      pc = -1;
    }
    //
    _ips = find(tokens == "ps");
    if ((length(_ips)==1) && _ips < length(tokens))
    {
      ps = strtod(tokens[_ips+1]);
    }
    else
    {
      ps = PLOT_DATASET_PLOT_FORMAT[4];
    }

    if ((pc==-1)&&(pt!=-1))
    { pc = PLOT_DATASET_PLOT_FORMAT[6]; }
    if ((pt==-1)&&(pc!=-1))
    { pt = PLOT_DATASET_PLOT_FORMAT[5]; }
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
      printf(THIS_LIBRARY + ": " + this_function + ": Cannot convert algabraic expression in %s that involves columns!\n", ...
          tokens[_iu+1]);
      printf(THIS_LIBRARY + ": " + this_function + ": Why don't you read the f*!king manual, while I rest for a while!\n", ...
          tokens[_iu+1]);
      rerror(THIS_LIBRARY + ": " + this_function + "\n");
    }
    use_cols = strtod(tokens[_iu+1],<<csp=":">>);
    if (err)
    {
      if ((err_xy==1)&&(length(use_cols)!=4))
      {
        printf(THIS_LIBRARY + ": " + this_function + ": 'with xerror{lines,bars}' requires 4 data columns (x,y,xmin,xmax).\n" );
        printf(THIS_LIBRARY + ": " + this_function + ": 'using %s' provided only %g columns. Cannot continue!\n", ...
            tokens[_iu+1],  length(use_cols));
        rerror(THIS_LIBRARY + ": " + this_function + "\n");
      }
      if ((err_xy==3)&&(length(use_cols)!=4))
      {
        printf(THIS_LIBRARY + ": " + this_function + ": 'with [y]error{lines,bars}' requires 4 data columns (x,y,ymin,ymax).\n" );
        printf(THIS_LIBRARY + ": " + this_function + ": 'using %s' provided only %g columns. Cannot continue!\n", ...
            tokens[_iu+1],  length(use_cols));
        rerror(THIS_LIBRARY + ": " + this_function + "\n");
      }
      if ((err_xy==2)&&(length(use_cols)!=6))
      {
        printf(THIS_LIBRARY + ": " + this_function + ": 'with xyerror{lines,bars}' requires 6 data columns (x,y,xmin,xmax,ymin,ymax).\n" );
        printf(THIS_LIBRARY + ": " + this_function + ": 'using %s' provided only %g columns. Cannot continue!\n", ...
            tokens[_iu+1],  length(use_cols));
        rerror(THIS_LIBRARY + ": " + this_function + "\n");
      }
    }
  }
  else
  {
    if (err)
    {
      if ((err_xy==1) || (err_xy==3))
      {
        use_cols = [1:4]; // default configuration
      }
      else
      {
        use_cols = [1:6]; // default configuration
      }
    }
    else
    {
      use_cols = [1,2]; // default configuration
    }
  }

  //
  // is user using key word 'axes'?
  //
  axes = PLOT_AXES_REV.x1y1;  // default. in case user does not use it!
  _ia = find(tokens == "ax");
  if (isempty(_ia))
  { _ia = find(tokens == "axes"); }
  if (isempty(_ia))
  {
    _ia = find(tokens == "axis");
    if (!isempty(_ia))
    {
      printf(THIS_LIBRARY + ": " + this_function + ": If you haven't know that, plural of 'axis' is 'axes'.\n");
      printf(THIS_LIBRARY + ": " + this_function + ": I would send you away to practice grammar, but my boss told me to behave!\n");
    }
  }
  if ((length(_ia)==1) && max(_ia) < length(tokens))
  {
    // axis enforces use_cols=[1,2] in the absence of 'using 1:2...' in plformat
    if (isempty(use_cols))
    { use_cols = [1,2]; }
    if (exist(PLOT_AXES_REV.[tokens[_ia+1]]))
    {
      axes = PLOT_AXES_REV.[tokens[_ia+1]];
    }
    else
    {
      printf(THIS_LIBRARY + ": " + this_function + ": directive 'axes %s' is not supported!\n",  tokens[_ia+1]);
      printf(THIS_LIBRARY + ": " + this_function + ": Why don't you go away and read the manual for a change!\n");
      error (THIS_LIBRARY + ": " + this_function + ": Cannot continue!\n");
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
    }
    else if ( any(strindex(tokens[(_id+1):_nd], "x")>0) )
    {
        err_xy = err_xy - 1;
    }
    else if ( any(strindex(tokens[(_id+1):_nd], "y")>0) )
    {
      err_xy = err_xy - 2;
    }

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

//
static (reset_active_plot_object);
reset_active_plot_object = function ( nx, ny, xleng, yleng, xoff, yoff )
{
  if (length(PLOT_WINDOWS) == 0)
  { return 1;}
  if (PLOT_ACTIVE_WIN == 0)
  { return 1;}

  if (exist(PLOT_WINDOWS.[PLOT_ACTIVE_WIN]))
  {
    clear(PLOT_WINDOWS.[PLOT_ACTIVE_WIN]);
  }
  pobj = <<>>;

  if (!exist (nx))
  { nx = 1; }
  if (!exist (ny))
  { ny = 1; }

  if (!exist (xleng))
  {
    if (exist(pobj.xleng))
    {
      xleng = pobj.xleng;
    }
    else
    {
      xleng = PLOT_WIN_LEN_X_DEFAULT;
    }
  }
  if (!exist (yleng))
  {
    if (exist(pobj.yleng))
    {
      yleng = pobj.yleng;
    }
    else
    {
      yleng = PLOT_WIN_LEN_Y_DEFAULT;
    }
  }
  if (!exist (xoff))
  {
    if (exist(pobj.xoff))
    {
      xoff = pobj.xoff;
    }
    else
    {
      xoff = PLOT_WIN_OFFSET_X_DEFAULT;
    }
  }
  if (!exist (yoff))
  {
    if (exist(pobj.yoff))
    {
      yoff = pobj.yoff;
    }
    else
    {
      yoff = PLOT_WIN_OFFSET_Y_DEFAULT;
    }
  }

  pobj.plformat = <<>>;
  pobj.pltext = <<>>;
  pobj.contour_zlevels = <<>>;
  pobj.contour_labels = <<>>;
  pobj.contour_ts = <<>>;
  pobj.contour_sp = <<>>;
  pobj.axes = <<>>;
  pobj.cols = <<>>;
  pobj.style = <<>>;
  pobj.xtick = <<>>;
  pobj.xtick_scale = <<>>;
  pobj.ytick = <<>>;
  pobj.ytick_scale = <<>>;
  pobj.ztick = <<>>;
  pobj.ztick_scale = <<>>;
  pobj.alt_xtick = <<>>;
  pobj.alt_xtick_scale = <<>>;
  pobj.alt_xlabel_scale = <<>>;
  pobj.alt_ytick = <<>>;
  pobj.alt_ytick_scale = <<>>;
  pobj.alt_ylabel_scale = <<>>;
  pobj.alt_zlabel_scale = <<>>;

  pobj.window.xp    = PLOT_WIN_XP_DEFAULT;      // Number of X pixels
  pobj.window.yp    = PLOT_WIN_YP_DEFAULT;      // Number of Y pixels
  pobj.window.xleng = xleng;        // Page length, X
  pobj.window.yleng = yleng;        // Page length, Y
  pobj.window.xoff  = xoff;           // Page offset, X
  pobj.window.yoff  = yoff;           // Page offset, Y

  pobj.subplot  = 0;      // The current subplot no.
  pobj.nplot    = nx*ny;  // Total no. of plots on window
  pobj.nx       = nx;
  pobj.ny       = ny;
  pobj.dev   = THIS_LIBRARY_DEFAULT_DEVICE;
  pobj.fontld = 0;        // Loaded extended fonts?
  pobj.char_height = 1.5;       // Character height.

  for (i in 1:(nx*ny))
  {
    // rlabplus:
    pobj.plformat.[i] = [];
    pobj.pltext.[i] = <<>>;         // texts that go on top of the graph
    pobj.axes.[i] = [];
    pobj.xtick.[i] = [0,-1,-1,-1,0,-1,-1,-1];
    pobj.xtick_scale.[i] = 1;
    pobj.alt_xtick.[i] = [0,-1,-1,-1,0,-1,-1,-1];
    pobj.alt_xtick_scale.[i] = 1;
    pobj.ytick.[i] = [0,-1,-1,-1,0,-1,-1,-1];
    pobj.ytick_scale.[i] = 1;
    pobj.alt_ytick.[i] = [0,-1,-1,-1,0,-1,-1,-1];
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
    pobj.ylabel_scale.[i] = [1.0, 3.0, 0.5, 0.5];
    pobj.alt_ylabel[i] = "";
    pobj.alt_ylabel_scale.[i] = [1.0, 3.0, 0.5, 0.5];
    pobj.zlabel[i] = "";
    pobj.alt_zlabel[i] = "";
    pobj.title[i] = "";
    pobj.title_scale.[i] = [1.5,1];
    pobj.orientation[i] = "portrait";
    pobj.desc.[i] = blank();        // The legend description
    pobj.desc_pos.[i] = "itr";      // The legend default position is inside/top/right
    pobj.desc_pos_xy.[i] = [0,0];   // The legend default position offset from the position above
    pobj.desc_scale.[i] = 1.0;      // The legend default size and scaling
    pobj.gridx[i] =  PLOT_GRID_X_DEFAULT; // Plot axes style, 2D-X
    pobj.alt_gridx[i] =  PLOT_ALT_GRID_X_DEFAULT; // Plot axes style, 2D-X
    pobj.gridy[i] =  PLOT_GRID_Y_DEFAULT; // Plot axes style, 2D-Y
    pobj.alt_gridy[i] =  PLOT_ALT_GRID_Y_DEFAULT; // Plot axes style, 2D-Y
    pobj.grid3x[i] = PLOT_GRID3D_X_DEFAULT; // Plot axes style, 3D-X
    pobj.grid3y[i] = PLOT_GRID3D_Y_DEFAULT; // Plot axes style, 3D-Y
    pobj.grid3z[i] = PLOT_GRID3D_Z_DEFAULT; // Plot axes style, 3D-Z
    pobj.aspect[i] = 0;           // Plot aspect style
    pobj.alt[i] = 60;
    pobj.az[i] = 45;
    pobj.viewport[i;] = nan(1,4);

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
    pobj.page.xp = PLOT_WIN_XP_DEFAULT;
    pobj.page.yp = PLOT_WIN_YP_DEFAULT;
    pobj.xleng   = xleng/nx;
    pobj.yleng   = yleng/ny;
    pobj.xoff    = xoff;
    pobj.yoff    = yoff;

    pobj.color[i;]  = 1:length(PLOT_COLOR_RGB_NAME); // 14 possible colors...
    pobj.lstyle[i;] = 1:5;               // 5 possible line styles...
    pobj.pstyle[i;] = 1:8;               // 8 possible point styles...
    pobj.fstyle[i;] = 1:8;               // 8 possible fill patterns...
  }

  //
  // Save the newly generated plot-object
  // in a list of plot-objects.
  PLOT_WINDOWS.[PLOT_ACTIVE_WIN] = pobj;

  return 0;
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
  { xleng = PLOT_WIN_LEN_X_DEFAULT; }
  if (!exist (yleng))
  { yleng = PLOT_WIN_LEN_Y_DEFAULT; }
  if (!exist (xoff))
  { xoff = PLOT_WIN_OFFSET_X_DEFAULT; }
  if (!exist (yoff))
  { yoff = PLOT_WIN_OFFSET_Y_DEFAULT; }

  pobj = <<>>;
  pobj.plformat = <<>>;
  pobj.contour_zlevels = <<>>;
  pobj.contour_labels = <<>>;
  pobj.contour_ts = <<>>;
  pobj.contour_sp = <<>>;
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

  pobj.window.xp    = PLOT_WIN_XP_DEFAULT;      // Number of X pixels
  pobj.window.yp    = PLOT_WIN_YP_DEFAULT;      // Number of Y pixels
  pobj.window.xleng = xleng;        // Page length, X
  pobj.window.yleng = yleng;        // Page length, Y
  pobj.window.xoff  = xoff;           // Page offset, X
  pobj.window.yoff  = yoff;           // Page offset, Y

  pobj.subplot  = 0;      // The current subplot no.
  pobj.nplot    = nx*ny;  // Total no. of plots on window
  pobj.nx       = nx;
  pobj.ny       = ny;
  pobj.dev   = THIS_LIBRARY_DEFAULT_DEVICE;
  pobj.fontld = 0;        // Loaded extended fonts?
  pobj.char_height = 1.5;       // Character height.

  for (i in 1:(nx*ny))
  {
    // rlabplus:
    pobj.plformat.[i] = [];
    pobj.pltext.[i] = <<>>;         // texts that go on top of the graph
    pobj.axes.[i] = [];
    pobj.xtick.[i] = [0,-1,-1,-1,0,-1,-1,-1];
    pobj.xtick_scale.[i] = 1;
    pobj.alt_xtick.[i] = [0,-1,-1,-1,0,-1,-1,-1];
    pobj.alt_xtick_scale.[i] = 1;
    pobj.ytick.[i] = [0,-1,-1,-1,0,-1,-1,-1];
    pobj.ytick_scale.[i] = 1;
    pobj.alt_ytick.[i] = [0,-1,-1,-1,0,-1,-1,-1];
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
    pobj.ylabel_scale.[i] = [1.0, 3.0, 0.5, 0.5];
    pobj.alt_ylabel[i] = "";
    pobj.alt_ylabel_scale.[i] = [1.0, 3.0, 0.5, 0.5];
    pobj.zlabel[i] = "";
    pobj.alt_zlabel[i] = "";
    pobj.title[i] = "";
    pobj.title_scale.[i] = [1.5,1];
    pobj.orientation[i] = "portrait";
    pobj.desc.[i] = blank();        // The legend description
    pobj.desc_pos.[i] = "itr";      // The legend default position is inside/top/right
    pobj.desc_pos_xy.[i] = [0,0];   // The legend default position offset from the position above
    pobj.desc_scale.[i] = 1.0;      // The legend default size and scaling
    pobj.gridx[i] =  PLOT_GRID_X_DEFAULT; // Plot axes style, 2D-X
    pobj.alt_gridx[i] =  PLOT_ALT_GRID_X_DEFAULT; // Plot axes style, 2D-X
    pobj.gridy[i] =  PLOT_GRID_Y_DEFAULT; // Plot axes style, 2D-Y
    pobj.alt_gridy[i] =  PLOT_ALT_GRID_Y_DEFAULT; // Plot axes style, 2D-Y
    pobj.grid3x[i] = PLOT_GRID3D_X_DEFAULT; // Plot axes style, 3D-X
    pobj.grid3y[i] = PLOT_GRID3D_Y_DEFAULT; // Plot axes style, 3D-Y
    pobj.grid3z[i] = PLOT_GRID3D_Z_DEFAULT; // Plot axes style, 3D-Z
    pobj.aspect[i] = 0;           // Plot aspect style
    pobj.alt[i] = 60;
    pobj.az[i] = 45;
    pobj.viewport[i;] = nan(1,4);

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
    pobj.page.xp = PLOT_WIN_XP_DEFAULT;
    pobj.page.yp = PLOT_WIN_YP_DEFAULT;
    pobj.xleng   = xleng/nx;
    pobj.yleng   = yleng/ny;
    pobj.xoff    = xoff;
    pobj.yoff    = yoff;

    pobj.color[i;]  = 1:length(PLOT_COLOR_RGB_NAME); // 14 possible colors...
    pobj.lstyle[i;] = 1:5;               // 5 possible line styles...
    pobj.pstyle[i;] = 1:8;               // 8 possible point styles...
    pobj.fstyle[i;] = 1:8;               // 8 possible fill patterns...
  }

  //
  // Save the newly generated plot-object
  // in a list of plot-objects.
  PLOT_WINDOWS.[N] = pobj;
};


static(find_scales);
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
    if (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].xmin[p] != 1j)
    {
      xmin = real (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].xmin[p]);
    }
    else
    {
      xmin = min(min(M[;use_cols]));
    }
    if (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].xmax[p] != 1j)
    {
      xmax = PLOT_WINDOWS.[PLOT_ACTIVE_WIN].xmax[p];
    }
    else
    {
      xmax = max(max(M[;use_cols]));
    }
  }
  else
  {
    // do the y-axis
    if (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].ymin[p] != 1j)
    {
      xmin = real (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].ymin[p]);
    }
    else
    {
      xmin = min(min(M[;use_cols]));
    }
    if (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].ymax[p] != 1j)
    {
      xmax = PLOT_WINDOWS.[PLOT_ACTIVE_WIN].ymax[p];
    }
    else
    {
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




////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Check to make sure a plot-object exists. If one
// does not exist, create it.
//
check_plot_object = function ()
{
  if (PLOT_ACTIVE_WIN == 0)
  {
    plwins(1);
    plwin (1);
  }

  return PLOT_ACTIVE_WIN;
};


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Convert x1,x2,y1,y2 to box coordinates
//
rect2line = function(leg_box_xl,leg_box_xr,leg_box_yb,leg_box_yt,d)
{
  if (!exist(d))
  {
    d = 1;
  }
  box_corners = nan(5,2);
  if (nargs==4)
  {
    if (leg_box_xl == leg_box_xr)
    {
      leg_box_xl = leg_box_xl - 0.5 * d;
      leg_box_xr = leg_box_xr + 0.5 * d;
    }
    if (leg_box_yb == leg_box_yt)
    {
      leg_box_yb = leg_box_yb - 0.5 * d;
      leg_box_yt = leg_box_yt + 0.5 * d;
    }
    box_corners = [...
        leg_box_xl,leg_box_yb; ...
        leg_box_xl,leg_box_yt; ...
            leg_box_xr,leg_box_yt; ...
                leg_box_xr,leg_box_yb; ...
                    leg_box_xl,leg_box_yb];
  }
  else if (nargs == 1)
  {
    if (leg_box_xl.n == 4)
    {
      if (leg_box_xl[1] == leg_box_xl[2])
      {
        leg_box_xl[1] = leg_box_xl[1] - 0.5 * d;
        leg_box_xl[2] = leg_box_xl[2] + 0.5 * d;
      }
      if (leg_box_xl[3] == leg_box_xl[4])
      {
        leg_box_xl[3] = leg_box_xl[3] - 0.5 * d;
        leg_box_xl[4] = leg_box_xl[4] + 0.5 * d;
      }
      box_corners = [...
          leg_box_xl[1],leg_box_xl[3]; ...
          leg_box_xl[1],leg_box_xl[4]; ...
              leg_box_xl[2],leg_box_xl[4]; ...
                  leg_box_xl[2],leg_box_xl[3]; ...
                      leg_box_xl[1],leg_box_xl[3]];
    }
  }
  return box_corners;
};


//
// Set the properties current plot window
// Default value = 1
plwin = function ( N, arg2, arg3)
{
  this_function = "plwin";
  iclear=1;
  vp=[0.15,0.95,0.15,0.95];
 
  if (!exist (N))
  { return PLOT_ACTIVE_WIN; }

  if (length(arg2)==1)
  {
    iclear=arg2;
  }
  if (length(arg3)==1)
  {
    iclear=arg3;
  }

  if (length(arg2)==4)
  {
    vp = arg2;
  }
  if (length(arg3)==4)
  {
    vp = arg3;
  }

  // Check to make sure N is valid
  if (exist(PLOT_WINDOWS.[N]))
  {
    pgslct(N);
    pgsvp(vp);

    if (iclear)
    {
      pgpage();
    }

    PLOT_ACTIVE_WIN = N;
    return PLOT_ACTIVE_WIN;
  }

  return 0;
};



newplwin = function ()
{
  if (length (PLOT_WINDOWS) == 0)
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
  this_function = "plwins";

  //
  // process NWIN: an integer, number of windows requested,
  // or a matrix [nx(i);ny(i)]_i
  //
  if (!exist(NWIN))
  {
    // report what is available
    if (length(PLOT_WINDOWS) == 0)
    {
      return <<win=[];act=0;dev=blank(0,0)>>;
    }

    retl = <<>>;
    retl.act = PLOT_ACTIVE_WIN;
    retl.win = strtod(members(PLOT_WINDOWS));
    retl.dev = blank(0,0);
    for(i in members(PLOT_WINDOWS))
    {
      retl.dev    = [retl.dev, PLOT_WINDOWS.[i].dev];
    }
    idx = sort(retl.win).idx;
    retl.dev = retl.dev[idx];
    retl.win = retl.win[idx];
    return retl;
  }

  if (NWIN.nr == 1 && NWIN.nc==1)
  {
    nwin = NWIN;
  }
  else if (NWIN.nc == 2)
  {
    nwin = NWIN.nr;
  }
  else
  {
    printf(this_function + ": improper first argument");
    return <<>>;
  }

  if (nwin==0)
  {
    PLOT_ACTIVE_WIN = 0;
    PLOT_WINDOWS = <<>>;
    plend();
    return <<act=0;dev=blank();win=[]>>;
  }
  else if (nwin<0 || nwin>8)
  {
    printf(this_function + ": improper first argument");
    return <<>>;
  }

  //
  // process sizes
  //
  if (!exist(sz))
  { sz = [PLOT_WIN_LEN_X_DEFAULT, PLOT_WIN_LEN_Y_DEFAULT]; }
  if (exist(sz) && (sz.nr*sz.nc!=2 || class(sz) != "num"))
  { sz = [PLOT_WIN_LEN_X_DEFAULT, PLOT_WIN_LEN_Y_DEFAULT]; }

  //
  // process offsets
  //
  if (!exist(offs))
  { offs = [PLOT_WIN_OFFSET_X_DEFAULT, PLOT_WIN_OFFSET_Y_DEFAULT]; }
  if (exist(sz) && (sz.nr*sz.nc!=2 || class(sz) != "num"))
  { offs = [PLOT_WIN_OFFSET_X_DEFAULT, PLOT_WIN_OFFSET_Y_DEFAULT]; }

  //
  // kill extra windows
  //
  wavl = strtod(members(PLOT_WINDOWS));
  if (nwin < length(wavl))
  {
    // kill extra windows
    for (i in nwin+1:length(wavl))
    {
      plwin(wavl[i]);
      clear(PLOT_WINDOWS.[wavl[i]]);
      plclose();
    }

    // rearrange remaining windows so 1:nwin are open
    for (i in 1:nwin)
    {
      if (i == wavl[i])
      { continue; }
      if (!exist(PLOT_WINDOWS.[i]))
      {
        // copy wavl[i] -> i
        PLOT_WINDOWS.[i] = PLOT_WINDOWS.[wavl[i]];
        plcopy( wavl[i], PLOT_WINDOWS.[i].dev );  // this will get copied to device 'i'
        // kill wavl[i]
        plwin ( wavl[i] );
        plclose ();
        if (wavl[i] == PLOT_ACTIVE_WIN)
        { PLOT_ACTIVE_WIN = 0; }
      }
    }
  }

  //
  // open/modify requested number of plot windows
  //
  wavl = strtod(members(PLOT_WINDOWS));
  if (!exist(dev))
  { dev = THIS_LIBRARY_DEFAULT_DEVICE; }
  for(i in members(PLOT_WINDOWS))
  { dev = [dev, PLOT_WINDOWS.[i].dev]; }
  for (i in 1:nwin)
  {
    ndev = dev[ min(i,length(dev)) ];
    nwid = sz[ min(i,sz.nr); 1 ];
    nhig = sz[ min(i,sz.nr); 2 ];

    nx = 1;
    ny = 1;
    if (NWIN.nc==2)
    {
      nx = NWIN[i;1];
      ny = NWIN[i;2];
    }

    if (!exist(PLOT_WINDOWS.[i]))
    {
      plstart(nx,ny,ndev,nwid,nhig,offs);
      continue;
    }
    else
    {
      // window exists and is open, check if the specifications have changed
      // and if so close it and open a new one
      plwin(i);
      if (PLOT_WINDOWS.[i].nx != nx || PLOT_WINDOWS.[i].ny != ny)
      {
        plclose();
        plstart(nx,ny,ndev,nwid,nhig,offs);
      }
    } // if (!exist(PLOT_WINDOWS.[i]))

  } // for (i in 1:nwin)


  //
  // return list
  //
  retl = <<>>;
  retl.win = strtod(members(PLOT_WINDOWS));
  if (PLOT_ACTIVE_WIN==0)
  { PLOT_ACTIVE_WIN=retl.win[1]; }
  retl.dev = blank();
  for(i in members(PLOT_WINDOWS))
  { retl.dev = [retl.dev, PLOT_WINDOWS.[i].dev]; }
  retl.act = PLOT_ACTIVE_WIN;
  return retl;
};

getplot = function ( win_no )
{
  local (win_no)

  if (length (PLOT_WINDOWS) != 0)
  {
    if (!exist (win_no)) { win_no = P; }

    if (exist (PLOT_WINDOWS.[win_no]))
    {
      return (PLOT_WINDOWS.[win_no]);
    }
    else
    {
      return 0;
    }
  }
  return <<>>;
};

////////////
//
// Set/start/select the plot device
//
// plstart ( NX, NY, DEV, WIDTH, ASPECT, OFFSET )
//
// Plstart opens a plot-device (window) and prepares it for
// plotting. Plstart returns the window ID (a scalar). NX and NY
// specify the number of plots in the horizontal and vertical
// directions respectively. DEV specifies the plot-window device. Valid
// examples are "/XWIN", "/PS", "/CPS". WIDTH specifies the width of
// the plot window in inches, while ASPECT specifies the aspect ratio
// of the new plot window.  Defaults for each argument are:
//
//   NX:     1
//   NY:     1
//   DEV:    "/XWIN"
//   WIDTH:  (user's .Xdefaults or environment)
//   ASPECT: 1 (only if WIDTH is specified).
//
plstart = function ( nx, ny, dev, width, height, offs )
{
    if (!exist (nx))
    { nx = 1; }

    if (!exist (ny))
    { ny = 1; }

    if (!exist (dev))
    { dev = THIS_LIBRARY_DEFAULT_DEVICE; }

    // Create the plot-object
    // First, figure out the index
    if (!PLOT_ACTIVE_WIN)
    {
      PLOT_ACTIVE_WIN = 1;
    }
    else
    {
      //P = P + 1;
      PLOT_ACTIVE_WIN = newplwin();
    }

    create_plot_object (PLOT_ACTIVE_WIN, nx, ny);

    // Open the plot device.
    winid = pgopen (dev);
    if (winid != PLOT_ACTIVE_WIN)
    {
      printf("plstart: winid != PLOT_ACTIVE_WIN\n");
    }

    // Setup the window width/aspect ratio.
    if (exist (width) && exist(height))
    {
      if (all(offs!=0))
      {
        //pgvsiz(offs[1],offs[1]+width,offs[2],offs[2]+height);
        pgsvp(0,1,0,1);
      }
      else
      {
        pgpap (width, height/width);
      }
    }

    // Turn off page prompting.
    pgask (0);

    // Set up the number of plots per window.
    pgsubp (nx, ny);

    // let the WIN variable knows what kind of window is open
    PLOT_WINDOWS.[winid].dev = dev;

    // set the colors
    _pgcolor(0, PLOT_COLOR_RGB_BKG);
    for (i in 1:PLOT_COLOR_RGB.nr)
    { _pgcolor(i, PLOT_COLOR_RGB[i;]); }

    // Update the screen.
    pgpage ();

    return PLOT_ACTIVE_WIN;
};

////////////
//
// Close a plot device. We must destroy the current plot-object
// And switch the output stream back to the default.
//
plclose = function ()
{
  if (length(PLOT_WINDOWS) > 1)
  {
    //
    // Clear PLOT_WINDOWS.[PLOT_ACTIVE_WIN] and reset P to 1st plot-window
    //
    clear (PLOT_WINDOWS.[PLOT_ACTIVE_WIN]);
    pgclos ();
    PLOT_ACTIVE_WIN = strtod(members (PLOT_WINDOWS)[1]);
    pgslct (PLOT_ACTIVE_WIN);
    return PLOT_ACTIVE_WIN;
  }
  else if (length(PLOT_WINDOWS) == 1)
  {
    if (exist (PLOT_WINDOWS.[PLOT_ACTIVE_WIN]))
    {
      clear (PLOT_WINDOWS.[PLOT_ACTIVE_WIN]);
      PLOT_WINDOWS = <<>>;
      PLOT_ACTIVE_WIN = 0;
    }
    pgclos ();
    return 1;
  }

  return 0;
};

////////////
//
// Close ALL the plot-windows
//

plend = function ()
{
  pgend ();
  if (exist (PLOT_WINDOWS))
  { clear (PLOT_WINDOWS); }

  PLOT_ACTIVE_WIN = 0;
  PLOT_WINDOWS = <<>>;
};

////////////
//
// Change plot aspect ratio
//

plaspect = function ( aspect )
{
  check_plot_object ();
  i = mod (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].subplot, PLOT_WINDOWS.[PLOT_ACTIVE_WIN].nplot) + 1;
  if (!exist (aspect))
  {
    PLOT_WINDOWS.[PLOT_ACTIVE_WIN].aspect[i] = 0;
  }
  else
  {
    PLOT_WINDOWS.[PLOT_ACTIVE_WIN].aspect[i] = aspect;
  }
};

////////////
//
// Set the current suplot to draw to.
//
// subplot ( HORIZONTAL, VERTICAL )
//
// HORIZONTAL specifies the subpanel in the horizontal direction
// (columns) to use for the next plot. VERTICAL specifies the subpanel
// in the vertical direction (rows) to use for the next plot.
//
// Both HORIZONTAL and VERTICAL must be specified for subplot to have
// any effect at all. Otherwise, the next sub-panel is used for the
// next plot.
//
// subplot must be used to specifiy the placement of each and every
// plot whose placement differs from the normal sub-panel advancement
// (horizontal, then vertical).
//

subplot = function ( HORIZONTAL, VERTICAL )
{
  if (!exist (HORIZONTAL))
  {
    subplot_spec = 0;
    return 0;
  }

  if (!exist (VERTICAL))
  {
    subplot_spec = 0;
    return 0;
  }

  // Turn on the subplot indicator.
  subplot_spec = 1;

  // Set the subplot to use next.
  subplot_horizontal = HORIZONTAL;
  subplot_vertical = VERTICAL;

  return [subplot_horizontal, subplot_vertical];
};

////////////
//
// Change plot line style
//

plstyle = function ( style )
{
  check_plot_object ();
  i = mod (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].subplot, PLOT_WINDOWS.[PLOT_ACTIVE_WIN].nplot) + 1;

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
      PLOT_WINDOWS.[PLOT_ACTIVE_WIN].style.[i] = style;
    }
    return 1;
  }
  PLOT_WINDOWS.[PLOT_ACTIVE_WIN].style.[i] = "line";
  return 1;
};

////////////
//
// Control of the plot line style.
// There are 5 line styles

plline = function ( line_style )
{
  check_plot_object ();
  i = mod (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].subplot, PLOT_WINDOWS.[PLOT_ACTIVE_WIN].nplot) + 1;

  if (exist (line_style))
  {
    if (class (line_style) == "num")
    {
      if (line_style.n != 20)
      {
	error ("plline: LVEC must be 1x20 in size");
      }
      PLOT_WINDOWS.[PLOT_ACTIVE_WIN].lstyle[i;] = line_style;
    }
    return 1;
  }
  PLOT_WINDOWS.[PLOT_ACTIVE_WIN].lstyle[i;] = [1:5,1:5,1:5,1:5];
  return 1;
};

////////////
//
// Control of the plot line color
//

plcolor = function ( COLOR )
{
  check_plot_object ();
  i = mod (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].subplot, PLOT_WINDOWS.[PLOT_ACTIVE_WIN].nplot) + 1;

  if (exist (COLOR))
  {
    if (class (COLOR) == "num")
    {
      PLOT_WINDOWS.[PLOT_ACTIVE_WIN].color[i;] = COLOR;
    }
    return 1;
  }
  PLOT_WINDOWS.[PLOT_ACTIVE_WIN].COLOR[i;] = 1:14;
  return 1;
};

////////////
//
// Control of the plot point style.
// There are 8 line styles

plpoint = function ( point_style )
{
  check_plot_object ();
  i = mod (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].subplot, PLOT_WINDOWS.[PLOT_ACTIVE_WIN].nplot) + 1;

  if (exist (point_style))
  {
    if (class (point_style) == "num")
    {
      if (point_style.nc != 8)
      {
	error ("plpoint: PVEC must be 1x8 in size");
      }
      PLOT_WINDOWS.[PLOT_ACTIVE_WIN].pstyle[i;] = point_style;
    }
    return 1;
  }
  PLOT_WINDOWS.[PLOT_ACTIVE_WIN].pstyle[i;] = 1:8;
  return 1;
};

////////////
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

////////////
//
// Change fonts
//

plfont = function ( font )
{
  check_plot_object ();
  i = mod (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].subplot, PLOT_WINDOWS.[PLOT_ACTIVE_WIN].nplot) + 1;

  if (!exist (font)) { font = 1; }

  if (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].fontld == 0)
  {
    pgscf (1);
    PLOT_WINDOWS.[PLOT_ACTIVE_WIN].fontld = 1;
  }

  PLOT_WINDOWS.[PLOT_ACTIVE_WIN].font[i] = font;
  return PLOT_ACTIVE_WIN;
};

////////////
//
// Change character height
//

plch = function ( H )
{
  check_plot_object ();
  i = mod (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].subplot, PLOT_WINDOWS.[PLOT_ACTIVE_WIN].nplot) + 1;

  if (!exist (H)) { H = 1; }

  PLOT_WINDOWS.[PLOT_ACTIVE_WIN].char_height = H;

  return PLOT_ACTIVE_WIN;
};

////////////
//
// Change pen width
//

plwid = function ( width )
{
  check_plot_object ();
  i = mod (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].subplot, PLOT_WINDOWS.[PLOT_ACTIVE_WIN].nplot) + 1;

  if (!exist (width))
  {
    PLOT_WINDOWS.[PLOT_ACTIVE_WIN].width = ones(1,32);
  } else {
    PLOT_WINDOWS.[PLOT_ACTIVE_WIN].width[i;] = width;
  }
  return PLOT_ACTIVE_WIN;
};

plfillstyle = function ( FS )
{
  check_plot_object ();
  i = mod (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].subplot, PLOT_WINDOWS.[PLOT_ACTIVE_WIN].nplot) + 1;

  if (!exist (FS))
  {
    FS = 1;       // Default to solid.
  }

  // Error checking/massaging.
  if (class (FS) != "num") {
    error ("plfillstyle: FS must be scalar, numeric");
  }

  if (FS.n > 1) {
    FS = FS[1];
  }

  if (FS > 4 || FS < 1) {
    FS = 1;
  }

  PLOT_WINDOWS.[PLOT_ACTIVE_WIN].fillstyle[i] = FS;

  return PLOT_ACTIVE_WIN;
};

////////////
//
// Place some text on the plot
//

plptex = function ( text, x , y , angle , just )
{
  p = mod (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].subplot, PLOT_WINDOWS.[PLOT_ACTIVE_WIN].nplot) + 1;	// The current index

  if (!check_plot_object ())
  {
    printf ("Must use plot() before plptex()\n");
    return 0;
  }

  if (!exist (x)) { x = 0; }
  if (!exist (y)) { y = 0; }
  if (!exist (angle)) { angle = 0; }
  if (!exist (just)) { just = 0; }

  // Set line color.
  pgsci (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].color[p;1]);

  // Write the text.
  pgptxt (x, y, angle, just, text);

  // Update the screen.
  pgupdt ();
};

////////////
//
// Find a character in a string
//

find_char = function ( str , char )
{
  tmp = strsplt (str);
  for (i in 1:tmp.n)
  {
    if (tmp[i] == char)
    {
      return i;
    }
  }
  return 0;
};

////////////
//
// Sort list element names/labels by numeric order, then string order.
//

list_sort = function ( L ){
    tl = <<>>;
    k = 1;
    for (i in members (L)){
        char[k] = i;
        k++;
        }

    if (exist (char)){
        tl.char = char;
        }

  return tl;
};

////////////
//
// Find the maximum number of elements in a bin for a single
// column matrix.
//

hist_scales = function ( data, nbin )
{
  // allow data to have multi column
  dmin = min( min (real (data)));
  dmax = max( max (real (data)));

  dbin = linspace (dmin, dmax, nbin+1);
  binval = zeros (nbin, 1);

  binval = zeros(nbin,1);
  for (j in data.nc)
  {
    for (i in 1:nbin)
    {
      binval[i] = binval[i] + ...
          length (find (data[;j] >= dbin[i] && data[;j] < dbin[i+1]));
    }
  }

  return max (binval);
};

////////////
//
// Copy a plot-object
//
// ID, the plot window ID to copy.
// DEV, the type of output device to open.
//
plcopy = function ( ID, DEV )
{
  //
  // Make sure ID exists in WIN
  //
  if ( !exist (ID) || !exist (DEV))
  {
    printf ("plcopy: No windows specified. Copy what?\n");
    return 1;
  }

  if (!exist (PLOT_WINDOWS.[ID]))
  {
    printf ("plcopy: %g is not valid plot-object. Cannot continue!\n", ID);
    return 1;
  }

  //
  // Check that there are no other devices of the same,non-/XWIN type, open.
  // This is because non-/XWIN devices print to a file the name of which depends
  // on the type of device (e.g., /PS,/CPS -> pgplot.ps ).  kmk 2005.
  //
  if(DEV!="/XWIN")
  {
    ddev = strsplt(DEV,"/");
    rDEV = "/" + ddev[ ddev.nr * ddev.nc ];
    if(sum(strindex(plwins().dev,rDEV))>0)
    {
      printf("plcopy: Cannot perform copy. Device \"%s\" already opened as a file\n", rDEV);
      printf("plcopy: \"%s\" !  Ignoring request for another copy.  Please close\n", DEV);
      printf("plcopy: the device first with  plclose() .\n");
      return 1;   // unsuccessful
    }
  }

  //
  // Get the new object ID
  //
  PLOT_ACTIVE_WIN = pgopen (DEV);

  //
  // Perform the copy.
  //
  PLOT_WINDOWS.[PLOT_ACTIVE_WIN] = PLOT_WINDOWS.[ID];

  // Note the device.  kmk 2005.
  PLOT_WINDOWS.[PLOT_ACTIVE_WIN].dev = DEV;

  // Turn off page prompting.
  pgask (0);

  // Set up the number of plots per window.
  pgsubp (PLOT_WINDOWS.[ID].nx, PLOT_WINDOWS.[ID].ny);

  return 0;   // successful
};

////////////
//
// Set the X-axis label
//

plxlabel = function ( xstr, scale )
{
  if (strlen(xstr)<1)
  { return 1; }

  check_plot_object ();

  i = mod (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].subplot, PLOT_WINDOWS.[PLOT_ACTIVE_WIN].nplot) + 1;

  PLOT_WINDOWS.[PLOT_ACTIVE_WIN].xlabel[i] = xstr;

  if (exist (scale))
  {
    if (length(scale)>0)
    {
      PLOT_WINDOWS.[PLOT_ACTIVE_WIN].xlabel_scale.[i][1] = scale[1];
    }
    if (length(scale)>1)
    {
      PLOT_WINDOWS.[PLOT_ACTIVE_WIN].xlabel_scale.[i][2] = scale[2];
    }
    if (length(scale)>2)
    {
      PLOT_WINDOWS.[PLOT_ACTIVE_WIN].xlabel_scale.[i][3] = scale[3];
    }
    if (length(scale)>3)
    {
      PLOT_WINDOWS.[PLOT_ACTIVE_WIN].xlabel_scale.[i][4] = scale[4];
    }
  }

  return PLOT_ACTIVE_WIN;
};

plx2label = function ( xstr, scale )
{
  if (strlen(xstr)<1)
  { return 1; }

  check_plot_object ();

  i = mod (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].subplot, PLOT_WINDOWS.[PLOT_ACTIVE_WIN].nplot) + 1;

  PLOT_WINDOWS.[PLOT_ACTIVE_WIN].alt_xlabel[i] = xstr;

  if (exist (scale))
  {
    if (length(scale)>0)
    {
      PLOT_WINDOWS.[PLOT_ACTIVE_WIN].alt_xlabel_scale.[i][1] = scale[1];
    }
    if (length(scale)>1)
    {
      PLOT_WINDOWS.[PLOT_ACTIVE_WIN].alt_xlabel_scale.[i][2] = scale[2];
    }
    if (length(scale)>2)
    {
      PLOT_WINDOWS.[PLOT_ACTIVE_WIN].alt_xlabel_scale.[i][3] = scale[3];
    }
    if (length(scale)>3)
    {
      PLOT_WINDOWS.[PLOT_ACTIVE_WIN].alt_xlabel_scale.[i][4] = scale[4];
    }
  }

  return PLOT_ACTIVE_WIN;
};

////////////
//
// Set the Y-axis label
//

plylabel = function ( xstr, scale )
{
  if (strlen(xstr)<1)
  { return 1; }

  check_plot_object ();

  i = mod (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].subplot, PLOT_WINDOWS.[PLOT_ACTIVE_WIN].nplot) + 1;

  PLOT_WINDOWS.[PLOT_ACTIVE_WIN].ylabel[i] = xstr;

  if (exist (scale))
  {
    if (length(scale)>0)
    {
      PLOT_WINDOWS.[PLOT_ACTIVE_WIN].ylabel_scale.[i][1] = scale[1];
    }
    if (length(scale)>1)
    {
      PLOT_WINDOWS.[PLOT_ACTIVE_WIN].ylabel_scale.[i][2] = scale[2];
    }
    if (length(scale)>2)
    {
      PLOT_WINDOWS.[PLOT_ACTIVE_WIN].ylabel_scale.[i][3] = scale[3];
    }
    if (length(scale)>3)
    {
      PLOT_WINDOWS.[PLOT_ACTIVE_WIN].ylabel_scale.[i][4] = scale[4];
    }
  }

  return PLOT_ACTIVE_WIN;
};


ply2label = function ( xstr, scale )
{
  if (strlen(xstr)<1)
  { return 1; }

  check_plot_object ();

  i = mod (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].subplot, PLOT_WINDOWS.[PLOT_ACTIVE_WIN].nplot) + 1;

  PLOT_WINDOWS.[PLOT_ACTIVE_WIN].alt_ylabel[i] = xstr;

  if (exist (scale))
  {
    if (length(scale)>0)
    {
      PLOT_WINDOWS.[PLOT_ACTIVE_WIN].alt_ylabel_scale.[i][1] = scale[1];
    }
    if (length(scale)>1)
    {
      PLOT_WINDOWS.[PLOT_ACTIVE_WIN].alt_ylabel_scale.[i][2] = scale[2];
    }
    if (length(scale)>2)
    {
      PLOT_WINDOWS.[PLOT_ACTIVE_WIN].alt_ylabel_scale.[i][3] = scale[3];
    }
    if (length(scale)>3)
    {
      PLOT_WINDOWS.[PLOT_ACTIVE_WIN].alt_ylabel_scale.[i][4] = scale[4];
    }
  }

  return PLOT_ACTIVE_WIN;
};




////////////
//
// Set the Z-axis label
//

plzlabel = function ( xstr )
{
  if (strlen(xstr)<1)
  { return 1; }

  check_plot_object ();

  i = mod (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].subplot, PLOT_WINDOWS.[PLOT_ACTIVE_WIN].nplot) + 1;

  PLOT_WINDOWS.[PLOT_ACTIVE_WIN].zlabel[i] = xstr;

  return 0;
};

////////////
//
// Set the plot-title: title, [scale, offset]
//
pltitle = function ( xstr, scale )
{
  if (strlen(xstr)<1)
  { return 1; }

  check_plot_object ();

  i = mod (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].subplot, ...
      PLOT_WINDOWS.[PLOT_ACTIVE_WIN].nplot) + 1;	// The current index

  PLOT_WINDOWS.[PLOT_ACTIVE_WIN].title[i] = xstr;

  if (exist (scale))
  {
    if (length(scale)>0)
    {
      PLOT_WINDOWS.[PLOT_ACTIVE_WIN].title_scale.[i][1] = scale[1];
    }
    if (length(scale)>1)
    {
      PLOT_WINDOWS.[PLOT_ACTIVE_WIN].title_scale.[i][2] = scale[2];
    }
  }

  return 0;
};

////////////
//
// Set the scale limits.
//

plimits = function ( xmin, xmax, ymin, ymax, zmin, zmax )
{
  check_plot_object ();
  i = mod (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].subplot, PLOT_WINDOWS.[PLOT_ACTIVE_WIN].nplot) + 1;

  if (exist (xmin))
  {
    PLOT_WINDOWS.[PLOT_ACTIVE_WIN].xmin[i] = xmin;
  }
  else
  {
    PLOT_WINDOWS.[PLOT_ACTIVE_WIN].xmin[i] = 1j;
  }

  if (exist (xmax))
  {
    PLOT_WINDOWS.[PLOT_ACTIVE_WIN].xmax[i] = xmax;
  }
  else
  {
    PLOT_WINDOWS.[PLOT_ACTIVE_WIN].xmax[i] = 1j;
  }

  if (exist (ymin))
  {
    PLOT_WINDOWS.[PLOT_ACTIVE_WIN].ymin[i] = ymin;
  }
  else
  {
    PLOT_WINDOWS.[PLOT_ACTIVE_WIN].ymin[i] = 1j;
  }

  if (exist (ymax))
  {
    PLOT_WINDOWS.[PLOT_ACTIVE_WIN].ymax[i] = ymax;
  }
  else
  {
    PLOT_WINDOWS.[PLOT_ACTIVE_WIN].ymax[i] = 1j;
  }

  if (exist (zmin))
  {
    PLOT_WINDOWS.[PLOT_ACTIVE_WIN].zmin[i] = zmin;
  }
  else
  {
    PLOT_WINDOWS.[PLOT_ACTIVE_WIN].zmin[i] = 1j;
  }

  if (exist (zmax))
  {
    PLOT_WINDOWS.[PLOT_ACTIVE_WIN].zmax[i] = zmax;
  }
  else
  {
    PLOT_WINDOWS.[PLOT_ACTIVE_WIN].zmax[i] = 1j;
  }

  return 0;
};

plimits2 = function ( xmin, xmax, ymin, ymax )
{
  check_plot_object ();
  i = mod (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].subplot, PLOT_WINDOWS.[PLOT_ACTIVE_WIN].nplot) + 1;

  if (exist (xmin))
  {
    PLOT_WINDOWS.[PLOT_ACTIVE_WIN].alt_xmin[i] = xmin;
  }
  else
  {
    PLOT_WINDOWS.[PLOT_ACTIVE_WIN].alt_xmin[i] = 1j;
  }
  if (exist (xmax))
  {
    PLOT_WINDOWS.[PLOT_ACTIVE_WIN].alt_xmax[i] = xmax;
  }
  else
  {
    PLOT_WINDOWS.[PLOT_ACTIVE_WIN].alt_xmax[i] = 1j;
  }

  if (exist (ymin))
  {
    PLOT_WINDOWS.[PLOT_ACTIVE_WIN].alt_ymin[i] = ymin;
  }
  else
  {
    PLOT_WINDOWS.[PLOT_ACTIVE_WIN].alt_ymin[i] = 1j;
  }
  if (exist (ymax))
  {
    PLOT_WINDOWS.[PLOT_ACTIVE_WIN].alt_ymax[i] = ymax;
  }
  else
  {
    PLOT_WINDOWS.[PLOT_ACTIVE_WIN].alt_ymax[i] = 1j;
  }
};



////////////
//
// Set 3-D grid (axis) styles
//

plgrid3 = function ( sty_x, sty_y, sty_z )
{
  check_plot_object ();
  i = mod (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].subplot, PLOT_WINDOWS.[PLOT_ACTIVE_WIN].nplot) + 1;
  if (exist (sty_x))
  {
    if (class (sty_x) == "string")
    {
      PLOT_WINDOWS.[PLOT_ACTIVE_WIN].grid3x[i] = sty_x;
    } else {
      error ("plgrid3: requires string argument GRID_STY_X");
    }
  } else {
    PLOT_WINDOWS.[PLOT_ACTIVE_WIN].grid3x[i] = PLOT_GRID3D_X_DEFAULT;
  }
  if (exist (sty_y))
  {
    if (class (sty_y) == "string")
    {
      PLOT_WINDOWS.[PLOT_ACTIVE_WIN].grid3y[i] = sty_y;
    } else {
      error ("plgrid3: requires string argument GRID_STY_Y");
    }
  } else {
    PLOT_WINDOWS.[PLOT_ACTIVE_WIN].grid3y[i] = PLOT_GRID3D_Y_DEFAULT;
  }
  if (exist (sty_z))
  {
    if (class (sty_z) == "string")
    {
      PLOT_WINDOWS.[PLOT_ACTIVE_WIN].grid3z[i] = sty_z;
    } else {
      error ("plgrid3: requires string argument GRID_STY_Z");
    }
  } else {
    PLOT_WINDOWS.[PLOT_ACTIVE_WIN].grid3z[i] = PLOT_GRID3D_Z_DEFAULT;
  }
};

////////////
//
// A friendlier interface to changing 2-D grid/axis
// styles.
//

plaxis = function ( X_STR, Y_STR )
{
  check_plot_object ();
  i = mod (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].subplot, PLOT_WINDOWS.[PLOT_ACTIVE_WIN].nplot) + 1;

  if (exist (X_STR))
  {
    if (X_STR == "log") { PLOT_WINDOWS.[PLOT_ACTIVE_WIN].gridx[i] = "bcnstl"; }
  } else {
    PLOT_WINDOWS.[PLOT_ACTIVE_WIN].gridx[i] = PLOT_GRID_X_DEFAULT;
  }

  if (exist (Y_STR))
  {
    if (Y_STR == "log") { PLOT_WINDOWS.[PLOT_ACTIVE_WIN].gridy[i] = "bcnstlv"; }
  } else {
    PLOT_WINDOWS.[PLOT_ACTIVE_WIN].gridy[i] = PLOT_GRID_Y_DEFAULT;
  }
  return PLOT_ACTIVE_WIN;
};

////////////
//
// Various internal support functions. Eventually these will be static.
//
////////////

//
// Find the X or Y scale limits .
// M can be a multi-column matrix, all columns
// will be used.
//

x_scales = function ( M, p )
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

  if (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].xmin[p] != 1j) { xmin = real (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].xmin[p]); }
  if (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].xmax[p] != 1j) { xmax = real (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].xmax[p]); }

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

  if (find_char (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].gridx[p], "l"))
  {
    if (xmin <= 0 || xmax <= 0)
    {
      error ("cannot plot log(x<=0)");
    }
    xmin = log10 (xmin);
    xmax = log10 (xmax);
  }

  return << xmin = xmin ; xmax = xmax >>;
};

y_scales = function ( M, p )
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

  if (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].ymin[p] != 1j) { xmin = real (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].ymin[p]); }
  if (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].ymax[p] != 1j) { xmax = real (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].ymax[p]); }

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

  if (find_char (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].gridy[p], "l"))
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

  if (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].zmin[p] != 1j) { xmin = real (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].zmin[p]); }
  if (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].zmax[p] != 1j) { xmax = real (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].zmax[p]); }

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

  if (find_char (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].gridz[p], "l"))
  {
    if (xmin <= 0 || xmax <= 0) { error ("plot: cannot plot log(z<=0)"); }
    xmin = log10 (xmin);
    xmax = log10 (xmax);
  }

  return << xmin = xmin ; xmax = xmax >>;
};

////////////
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

  if (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].xmin[p] != 1j) { xmin = real (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].xmin[p]); }
  if (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].xmax[p] != 1j) { xmax = real (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].xmax[p]); }
  if (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].ymin[p] != 1j) { ymin = real (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].ymin[p]); }
  if (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].ymax[p] != 1j) { ymax = real (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].ymax[p]); }
  if (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].zmin[p] != 1j) { zmin = real (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].zmin[p]); }
  if (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].zmax[p] != 1j) { zmax = real (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].zmax[p]); }

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

  if (find_char (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].grid3x[p], "l"))
  {
    if (xmin <= 0 || xmax <= 0) { error ("plot: cannot plot log(x<=0)"); }
    xmin = log10 (xmin);
    xmax = log10 (xmax);
  }

  if (find_char (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].grid3y[p], "l"))
  {
    if (ymin <= 0 || ymax <= 0) { error ("plot: cannot plot log(y<=0)"); }
    ymin = log10 (ymin);
    ymax = log10 (ymax);
  }

  if (find_char (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].grid3z[p], "l"))
  {
    if (zmin <= 0 || zmax <= 0) { error ("plot: cannot plot log(z<=0)"); }
    zmin = log10 (zmin);
    zmax = log10 (zmax);
  }

  return <<xmin=xmin; xmax=xmax; ymin=ymin; ymax=ymax; zmin=zmin; zmax=zmax>>;
};

////////////
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
    } else {
      M = real (data.[i]);
    }
    if (class (M) != "num") { continue; }
    if (M.n == 0) { continue; }

    if (abs (key) > M.nc)
    { error ("plot: KEY argument > M.nc"); }

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
      	 </ymax;ymin/> = y_scales ( real(M)[;k],   p);
      } else {
	       </xmax;xmin/> = x_scales ( (1:M.nr)', p);
	       </ymax;ymin/> = y_scales ( real(M),   p);
      }
    } else { if (key < 0)
    {
      </xmax;xmin/> = x_scales ( real(M)[;k],        p);
      </ymax;ymin/> = y_scales ( real(M)[;abs(key)], p);
    } else {
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

////////////
//
// Find the maximum number of elements in a bin for a single
// column matrix.
//

// hist_scales = function ( data, nbin )
// {
//   dmin = min (real (data));
//   dmax = max (real (data));
//   dbin = linspace (dmin, dmax, nbin+1);
//   binval = zeros (nbin, 1);
//
//   for (i in 1:nbin)
//   {
//     binval[i] = length (find (data >= dbin[i] && data < dbin[i+1]));
//   }
//
//   return max (binval);
// };

////////////
//
// Find appropriate tick marks, and sub-intervals for a given data set.
//
// Input is the min and max values of the data range.
// Returns the tick intervals, and number of subdivisions for sub-ticks.
// Also returns the new min and max values (these may be reset by this function).
//

ticks = function ( dmin, dmax )
{
  //
  // Compute the range (total interval) of the data.
  //

  range = abs (dmax - dmin);

  //
  // Now, normalize the range into 1-10 interval.
  //

  sf = log10 (range);
  if (sf <= 0)
  {
    sf = 10.0;
  } else { if (sf <= 1) {
    sf =  1.0;
  } else { if (sf > 1) {
    sf =  1 / (10^(int (sf)+0));
  } } }

  //
  // Heuristics to determine number of ticks in interval.
  //

  range = range * sf;

  if (range <= 1)
  {
    tick = 0.2;
  } else {
    if (range <= 3)
    {
      tick = 0.5;
    } else {
      if (range <= 5)
      {
        tick = 1;
      } else {
        tick = 2;
      }
    }
  }

  //
  // Rescale tick to match original (unscaled data).
  //

  tick  = tick / sf;
  range = range / sf;

  //
  // Re-adjust max scale value.
  //

  i = 1;
  dmax_save = dmax;
  while (tick*i < range)
  {
    dmax = dmin + (tick*(i+1));
    i++;
  }

  //
  // Re-adjust min scale value.
  //

  i = 1;
  while (tick*i < range)
  {
    dmin = dmax_save - (tick*(i+1));
    i++;
  }

  // Default behavior
  subtick = 2;

  return << dmin=dmin; dmax=dmax; subtick=subtick; tick=tick >>;
};

////////////
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
////////////

plot_matrix = function ( M, key, p, K, xmin, xmax, ymin, ymax, v )
{
  np = M.nr;

  if (M.nc == 1)
  {
    x = 1:M.nr;
    y = real (M);
    k = mod (1+K, 14) + 1;
    l = mod (1+K, 8);

    if (find_char (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].gridx[p], "l"))
    { x = log10 (x); }
    if (find_char (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].gridy[p], "l"))
    { y = log10 (y); }

    pgsci (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].color[p;k]);     // Set line color.
    pgsls (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].lstyle[p;l]);    // Set line style.
    pgslw (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].width[p;l]);     // Set line width.

    if (get_style (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].style.[p], k-1) == "line")
    {
      pgline (M.nr, x, y);
    } else { if (get_style (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].style.[p], k-1) == "point") {
      pgpt (M.nr, x, y, PLOT_WINDOWS.[PLOT_ACTIVE_WIN].pstyle[p;l]);
    } else { if (get_style (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].style.[p], k-1) == "line-point") {
      pgline (M.nr, x, y);
      pgpt (M.nr, x, y, PLOT_WINDOWS.[PLOT_ACTIVE_WIN].pstyle[p;l]);
    } else { {
      pgline (M.nr, x, y);
    }}}}

    //
    // Now do the legend
    //
    if (!any (any (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].desc.[p] == 1j)))
    {
      // Use the default if necessary
      if (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].desc.[p][1] == "default")
      {
	       desc = "c1";
      } else { if (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].desc.[p].n >= k-1)
      { desc = PLOT_WINDOWS.[PLOT_ACTIVE_WIN].desc.[p][k-1];
      } else {
	       // Not sure what to do, user has messed up.
	       desc = "";
      }}

      v = v - (ymax-ymin)/11;
      xl = (xmax-xmin)*[10.5/12, 11/12, 11.5/12]' + xmin;
      yl = [v, v, v]' + ymin;

      if (get_style (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].style.[p], k-1) == "line")
      {
	       pgline (3, xl, yl);
      } else { if (get_style (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].style.[p], k-1) == "point")
      { pgpt (3, xl, yl, PLOT_WINDOWS.[PLOT_ACTIVE_WIN].pstyle[p;l]);
      } else { if (get_style (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].style.[p], k-1) == "line-point")
      {
        pgline (3, xl, yl);
	      pgpt (3, xl, yl, PLOT_WINDOWS.[PLOT_ACTIVE_WIN].pstyle[p;l]);
      }}}

      pgptxt (xl[1]-(xmax-xmin)/25, yl[3], 0, 1, desc);

    }

  } else {

    //
    // Check for large column dimension
    //
    if (M.nc > 3*M.nr && 0)
    {
      printf (" Plot %i columns and %i rows, are you sure (y/n) ? "...
               , M.nc, M.nr);
      ans = getline ("stdin");
      if (ans.[1] != "y") { return -1; }
    }

    ki = find ((1:M.nc) != abs (key));
    for (i in ki)
    {
      if (key > 0)
      {
	       x = real (M[;key]);
	       y = real (M[;i]);
      } else { if (key < 0) {
	       x = real (M[;i]);
	       y = real (M[;abs(key)]);
      } else {
	       x = (1:M.nr)';
	       y = real (M[;i]);
      }}

      // Check for log scales, adjust if necessary
      if (find_char (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].gridx[p], "l"))
      { x = log10 (x); }
      if (find_char (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].gridy[p], "l"))
      { y = log10 (y); }

      k = mod (i-1 + K, 14) + 1;
      l = mod (8 + i-2 + K, 8) + 1;

      pgsci (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].color[p;k]);     // Set line color.
      pgsls (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].lstyle[p;l]);    // Set line style.
      pgslw (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].width[p;l]);     // Set line width.

      if (get_style (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].style.[p], k-1) == "line")
      {
        pgline (np, x, y);
      }
      else if (get_style (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].style.[p], k-1) == "point")
      {
        pgpt (np, x, y, PLOT_WINDOWS.[PLOT_ACTIVE_WIN].pstyle[p;l]);
      }
      else if (get_style (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].style.[p], k-1) == "line-point")
      {
        pgline (np, x, y);
        pgpt (np, x, y, PLOT_WINDOWS.[PLOT_ACTIVE_WIN].pstyle[p;l]);
      }
      else
      {
        pgline (np, x, y);
      }

      //
      // Now do the legend
      //
      if (!any (any (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].desc.[p] == 1j)))
      {
	       // Use the default if necessary
	    if (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].desc.[p][1] == "default")
	    {
	      desc = "c" + num2str (i);
      } else { if (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].desc.[p].n >= k-1) {
        if (k > 1)
        {
          desc = PLOT_WINDOWS.[PLOT_ACTIVE_WIN].desc.[p][k-1];
        } else {
          desc = "";
        }
      } else {
	      // Not sure what to do, user has messed up.
	      desc = "";
      }}

      v = v - (ymax-ymin)/11;
      xl = (xmax-xmin)*[10.5/12, 11/12, 11.5/12]' + xmin;
      yl = [v, v, v]' + ymin;

      if (get_style (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].style.[p], k-1) == "line")
      {
        pgline (3, xl, yl);
      } else { if (get_style (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].style.[p], k-1) == "point") {
        pgpt (3, xl, yl, PLOT_WINDOWS.[PLOT_ACTIVE_WIN].pstyle[p;l]);
      } else { if (get_style (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].style.[p], k-1) == "line-point") {
        pgline (3, xl, yl);
        pgpt (3, xl, yl, PLOT_WINDOWS.[PLOT_ACTIVE_WIN].pstyle[p;l]);
      }}}

      pgslw  (5);                                         // set 5pt for legend text label
      pgptxt (xl[1]-(xmax-xmin)/25, yl[3], 0, 1, desc);
      pgslw (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].width[p;l]);                         // Set previous line width.

      }
    }
  }

  return k-1;
};

////////////
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
      if (class(L.[i])=="list")
      {
        if(ishist(L.[i]))
        { M = hist_line(L.[i]); }
      }
      if (class (M) != "num") { continue; }
      if ((k = plot_matrix (M, key, p, k, xmin, xmax, ymin, ymax, v)) < 0)
      {	return k; }
      v = v - (ymax-ymin)/11;
    }
  }

  // Now plot the list members with string labels.
  if (exist (sl.char))
  {
    for (i in sl.char)
    {
      M = L.[i];
      if (class(L.[i])=="list")
      {
        if(ishist(L.[i]))
        { M = hist_line(L.[i]); }
      }
      if (class (M) != "num")
      { continue; }
      if ((k = plot_matrix (M, key, p, k, xmin, xmax, ymin, ymax, v)) < 0)
      { return k;}
      v = v - (ymax-ymin)/11;
    }
  }
  return 1;
};

////////////
//
// Create a legend in the x1y1 worldview plot window
//
//  _pglegend = function(desc[idx_n], _plformat[idx_n;], desc_scale, desc_pos, desc_pos_xy);
// fmt = [lpb, lt, lc, lw, pt, pc, ps, err, err_xy, err_lp];
_pglegend = function(desc, fmt, sc, pos, pos_xy, zlevels)
{
  this_function = "_pglegend";

  // parameters for the legend in units of device viewport
  leg_color_idx = 1;  // text color of legend
  leg_char_spacing = 1.0; // spacing between the bbox and lines of text in legend, and between text lines
  leg_plot_wid_vp = 0.35; // in, width of the lines or points
  leg_box_top_margin_vp = 0.15; // in, how far from top edge
  leg_box_right_margin_vp = 0.15; // in, how far from right edge
  leg_box_bottom_margin_vp = 0.15; // in, how far from bottom edge
  leg_box_left_margin_vp = 0.15; // in, how far from left edge

  scale = sc * THIS_LIBRARY_DEFAULT_CHAR_SIZE;

  inside = 0;
  if (strindex(pos,"i"))
  { inside = 1; }

  // top <-> bottom
  top = 0;
  if (strindex(pos,"t"))
  { top = 1; }
  if (strindex(pos,"m"))
  { top = 0.5; }

  // left <-> right
  just = 0;
  if (strindex(pos,"r"))
  { just = 1; }
  if (strindex(pos,"c"))
  { just = 0.5; }

  bbox = 0;
  if (strindex(pos,"b"))
  { bbox = 1; }

  // get size of character in viewport units
  pgsch  (1); // ps
  s = max([pglen(1,"a"),pglen(1,"]")]);

  // device viewport
  desc_vp = pgqvp(1); // xl,xr,yb,yt

  // number of legend entries (do not count pm3d)
  nleg = length(desc) - find(fmt[;1]==3);
  if (nleg>0)
  {
    // world:
    desc_win = pgqwin();
    xmin = desc_win[1];
    xmax = desc_win[2];
    ymin = desc_win[3];
    ymax = desc_win[4];

    // scale world -> viewport
    vp_win_x = (desc_vp[2] - desc_vp[1]) / (xmax-xmin);
    vp_win_y = (desc_vp[4] - desc_vp[3]) / (ymax-ymin);

    // size of the legend text in characters in viewport
    leg_txt_wid_vp = (max(strlen(desc)) + 0.5) * scale * s; // this is the text
    leg_box_wid_vp = (max(strlen(desc)) + 2) * scale * s + leg_plot_wid_vp;
    leg_box_hei_vp = ((1 + leg_char_spacing) * nleg + leg_char_spacing) * scale * s;

    leg_box_wid = leg_box_wid_vp / vp_win_x;
    leg_box_hei = leg_box_hei_vp / vp_win_y;

      // on the right assuming short legend
    leg_box_xl = xmin + (leg_box_left_margin_vp ...
        + just * (desc_vp[2] - desc_vp[1] - leg_box_left_margin_vp ...
        - leg_box_right_margin_vp - leg_box_wid_vp)) / vp_win_x;
    leg_box_xr = leg_box_xl +  leg_box_wid;

    leg_box_yb = ymin + (leg_box_bottom_margin_vp ...
        + top * (desc_vp[4] - desc_vp[3] - leg_box_bottom_margin_vp - leg_box_top_margin_vp - leg_box_hei_vp)) / vp_win_y;
    leg_box_yt = leg_box_yb +  leg_box_hei;

    //
    // right aligned legend: text line/points
    //

    // draw legend in color 1
    // clear the rectangle where the legend will go
    pgsfs (1);
    pgsci (0);
    pgrect (leg_box_xl,leg_box_xr,leg_box_yb,leg_box_yt);
    if (bbox)
    {
      pgsci  (1);
      box_corners = [...
          leg_box_xl,leg_box_yb; ...
          leg_box_xl,leg_box_yt; ...
              leg_box_xr,leg_box_yt; ...
                  leg_box_xr,leg_box_yb; ...
                      leg_box_xl,leg_box_yb];
      pgline(5,box_corners[;1], box_corners[;2]);
    }

    y0 = leg_box_yt;
    dy   = (1 + leg_char_spacing) * scale * s / vp_win_y;
    dy_l = 0.5 * scale * s / vp_win_y;

    x0 = leg_txt_wid_vp/vp_win_x + leg_box_xl;
    x1 = (leg_txt_wid_vp + 0.5 * scale * s) / vp_win_x + leg_box_xl;
    x2 = (leg_txt_wid_vp + 0.5 * scale * s + leg_plot_wid_vp) / vp_win_x + leg_box_xl;
    dx = x2 - x1;

    //
    // draw normal legend
    //
    for (i in 1:length(desc))
    {
      // contour
      if (fmt[i;1]==4)
      {
        y0 = y0 - dy;

        // draw single color line
        if (fmt[i;5]==1)
        {
          // draw line
          pgsls  (fmt[i;6]);  // lt
          pgslw  (fmt[i;7]);  // lw
          pgsci  (fmt[i;8]);  // lc
          pgline ([x1 + 0.125*dx,y0+dy_l;x2 - 0.125 * dx,y0+dy_l]);
        }
        else if ((fmt[i;5]==2)||(fmt[i;5]==3))
        {
          // grey contour requires zmin,zmax
          _fmt = fmt[i;];
          if (isnan(_fmt[2])&&isnan(_fmt[3]))
          {
            p = mod (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].subplot, PLOT_WINDOWS.[PLOT_ACTIVE_WIN].nplot) + 1; // The current index
            if (!imag(PLOT_WINDOWS.[PLOT_ACTIVE_WIN].zmin[p]))
            { _fmt[2] = real(PLOT_WINDOWS.[PLOT_ACTIVE_WIN].zmin[p]); }
            if (!imag(PLOT_WINDOWS.[PLOT_ACTIVE_WIN].zmax[p]))
            { _fmt[3] = real(PLOT_WINDOWS.[PLOT_ACTIVE_WIN].zmax[p]); }
          }

          // gray lines in this area of the legend box:
          _wv = zeros(1,4);
          _wv[1] = x1 + 0.125 * dx;
          _wv[2] = x2 - 0.125 * dx;
          _wv[3] = y0 - 0.1 * dy;
          _wv[4] = y0 + 0.7 * dy;
          _pgdrawgrad(_fmt, desc[i], 0.9 * scale, _wv, zlevels[i]);
        }

        // write text
        pgstbg (0);
        pgsci  (leg_color_idx); // default line plot color
        pgslw  (1); // default line width
        pgsch  (scale);
        pgptxt (x0, y0, 0, 1, desc[i]);
      }
      else if (fmt[i;1] <= 2)
      {
        // draw single color line
        y0 = y0 - dy;

        if ((fmt[i;1]== 0)||(fmt[i;1]== 2))
        {
          pgsci  (fmt[i;3]);  // lc
          pgslw  (fmt[i;4]);  // lw
          pgline ([x1 + 0.125*dx,y0+dy_l;x2 - 0.125 * dx,y0+dy_l]);
        }
        if ((fmt[i;1]== 1)||(fmt[i;1]== 2))
        {
          pgsci  (fmt[i;6]);        // pc (=lc in absence of pc in format)
          pgsch  (0.75 * fmt[i;7]); // ps
          pgpt   ([x1+0.5*dx,y0+dy_l], fmt[i;5]); // pt
        }
        pgstbg (0);
        pgsci  (leg_color_idx); // default line plot color
        pgslw  (1); // default line width
        pgsch  (scale);
        pgptxt (x0, y0, 0, 1, desc[i]);
      }
    }
    pgsch  (1); // default character size

  } // if (nleg>0): end of normal 2-D dataset legends

  //
  // draw gnuplot style-legend for pm3d on the right
  //
  for (i in 1:length(desc))
  {
    // ignore pm3d's
    if (fmt[i;1]!=3)
    { continue; }

    // grey contour requires zmin,zmax
    _fmt = fmt[i;];
    if (isnan(_fmt[2])&&isnan(_fmt[3]))
    {
      p = mod (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].subplot, PLOT_WINDOWS.[PLOT_ACTIVE_WIN].nplot) + 1; // The current index
      PLOT_WINDOWS.[PLOT_ACTIVE_WIN].zmin[p]?
      if (!imag(PLOT_WINDOWS.[PLOT_ACTIVE_WIN].zmin[p]))
      { _fmt[2] = real(PLOT_WINDOWS.[PLOT_ACTIVE_WIN].zmin[p]); }
      if (!imag(PLOT_WINDOWS.[PLOT_ACTIVE_WIN].zmax[p]))
      { _fmt[3] = real(PLOT_WINDOWS.[PLOT_ACTIVE_WIN].zmax[p]); }
    }

    //
    // viewport for drawing scale of the pm3d on the right
    //  of the main plot window
    //
    _cvp0 = pgqvp(0);
    xmax = ceil(_cvp0[2],<<bin=0.5>>);
    _cvp1 = _cvp0;
    _cvp1[1] = _cvp0[2] + 0.225 * (xmax - _cvp0[2]);
    _cvp1[2] = _cvp0[2] + 0.450 * (xmax - _cvp0[2]);
    _cvp1[3] = _cvp0[3];
    _cvp1[4] = _cvp0[4];
    pgsvp(_cvp1); // set new worldview
    _pgdrawgrad(_fmt, desc[i], 0.9 * scale);
    pgsvp(_cvp0); // restor old worlview
  }
  pgsch  (1); // default character size

  return 0;
};


//
// Set the current plot legend string
//
plegend = function ( LEGEND, scale, pos, pos_xy )
{
  check_plot_object ();
  p = mod (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].subplot, PLOT_WINDOWS.[PLOT_ACTIVE_WIN].nplot) + 1; // The current index

  if (!exist (LEGEND))
  {
    PLOT_WINDOWS.[PLOT_ACTIVE_WIN].desc.[p] = 1j;
    return PLOT_ACTIVE_WIN;
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
    PLOT_WINDOWS.[PLOT_ACTIVE_WIN].desc.[p] = LEGEND;
    PLOT_WINDOWS.[PLOT_ACTIVE_WIN].desc_pos.[p] = pos;
    PLOT_WINDOWS.[PLOT_ACTIVE_WIN].desc_pos_xy.[p] = pos_xy;
    PLOT_WINDOWS.[PLOT_ACTIVE_WIN].desc_scale.[p] = scale;
  }

  return PLOT_ACTIVE_WIN;
};

set3d = function (bx, by, h)
{
  if (!exist (bx))
  { PLOT_BASE_X = 2; }
  else { PLOT_BASE_X = bx; }

  if (!exist (by))
  { PLOT_BASE_Y = 2; }
  else
  { PLOT_BASE_Y = by; }

  if (!exist (h))
  { PLOT_BASE_Z = 4; }
  else
  { PLOT_BASE_Z = h; }
};

//---------------------------------------------------------------------
//
// Plot a surface...
//

plsurf = function (x, y, z, theta, phi)
{
  global (pi)

  check_plot_object ();
  p = mod (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].subplot, PLOT_WINDOWS.[PLOT_ACTIVE_WIN].nplot) + 1;	// The current index

  // Force to column vectors.
  x = x[:];
  y = y[:];

  if (!exist (theta)) { theta = 60; }
  if (!exist (phi))   { phi = 35; }

  dtor = 2*pi / 360;

  //
  // Setup the 3 matrices of coordinates.
  //

  xm = ones (x.n, y.n) .* x;
  ym = ones (x.n, y.n) .* y';
  zm = z;

  //
  // First the Z rotation.
  // (Ignore the Z coordinate)
  //

  for (i in 1:x.n)
  {
    for (j in 1:y.n)
    {
      // Do the first axis...
      tmp = rotz ( [xm[i;j], ym[i;j], zm[i;j]], -theta*dtor );

      // Now do the second rotation...
      tmp = rotx ( tmp, -phi*dtor );

      xm[i;j] = tmp[1];
      ym[i;j] = tmp[2];
      zm[i;j] = tmp[3];
    }
  }

  //
  // Calculate the scale limits.
  //

  xmin = min (min (xm));
  xmax = max (max (xm));

  ymin = min (min (ym));
  ymax = max (max (ym));

  if (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].xmin[p] != 1j) { xmin = real (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].xmin[p]); }
  if (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].xmax[p] != 1j) { xmax = real (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].xmax[p]); }

  if (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].ymin[p] != 1j) { ymin = real (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].ymin[p]); }
  if (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].ymax[p] != 1j) { ymax = real (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].ymax[p]); }

  //
  // Set up the window/viewport.
  //

  pgpage ();
  pgvstd ();

  if (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].aspect[p] != 0)
  {
    pgwnad (xmin, xmax, ymin, ymax);
  } else {
    pgswin (xmin, xmax, ymin, ymax);
  }

  //
  // Now, pick the tickmark properties, maybe reset
  // scale min/max values.
  //

  </ xmax; xmin; nxsub; xtick /> = ticks (xmin, xmax);
  </ ymax; ymin; nysub; ytick /> = ticks (ymin, ymax);

  pgsci (10);                       // Set back to normal grid, label color
  pgslw (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].width[p]);         // Line width
  pgscf (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].font[p]);          // Set character font
  pgsch (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].char_height);      // Set character height

  pgbox (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].gridx[p], xtick, nxsub, ...
         PLOT_WINDOWS.[PLOT_ACTIVE_WIN].gridy[p], ytick, nysub);

  pgsci (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].color[p;1]);       // Set line color.

  //
  // Now, plot up the projected surface.
  //

  for (i in 1:xm.nc)
  {
    pgline (xm.nr, xm[;i], ym[;i]);
  }

  for (i in 1:xm.nr)
  {
    pgline (xm.nc, xm[i;], ym[i;]);
  }

  pgsls (1);
  pgsci (10);

  pglab(PLOT_WINDOWS.[PLOT_ACTIVE_WIN].xlabel[p], PLOT_WINDOWS.[PLOT_ACTIVE_WIN].ylabel[p], PLOT_WINDOWS.[PLOT_ACTIVE_WIN].title[p]);
  pgupdt ();
};

//---------------------------------------------------------------------
//
// Plot a surface with polygons...
//
plmesh = function (arg1, arg2, arg3, arg4, arg5, arg6)
{
  this_function = "plmesh";

  global (pi)

  check_plot_object ();
  p = mod (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].subplot, PLOT_WINDOWS.[PLOT_ACTIVE_WIN].nplot) + 1; // The current index

  if (!exist(arg1))
  {
    error(this_function + ": missing first argument 'x' or list <<x;y;z>>");
  }

  if (class(arg1)=="list")
  {
    // arg1 = <<x;y;z>>
    if (exist(arg1.x))
    {
      x = arg1.x[:];
    } else {
      error(this_function + ": entry 'x' missing in list that was provided as first argument");
    }

    if (exist(arg1.y))
    {
      y = arg1.y[:];
    } else {
      error(this_function + ": entry 'y' missing in list that was provided as first argument");
    }

    if (exist(arg1.z))
    {
      z = arg1.z;
    } else {
      error(this_function + ": entry 'z' missing in list that was provided as first argument");
    }

    // arg2 = theta
    if (exist(arg2))
    {
      theta = arg2;
    } else {
      theta = 60;
    }

    // arg3 = phi
    if (exist(arg3))
    {
      phi = arg3;
    } else {
      phi = 35;
    }

    // arg3 = phi
    if (exist(arg4))
    {
      axis = arg4;
    } else {
      axis = 0;
    }

  } else {

    if (!exist(arg1) || !exist(arg2) || !exist(arg3))
    { error(this_function + ": some of arguments 'x', 'y' and 'z' are missing"); }

    x = arg1[:];
    y = arg2[:];

    if (exist (arg4))
    {
      theta = arg4;
    } else {
      theta = 60;
    }

    if (exist(arg5))
    {
      phi = arg5;
    } else {
      phi = 35;
    }

    if (exist(arg6))
    {
      axis = arg6;
    } else {
      axis = 0;
    }

  }

  dtor = 2*pi / 360;

  //
  // Setup the 3 matrices of coordinates.
  //
  xm = ones (x.n, y.n) .* x;
  ym = ones (x.n, y.n) .* y';
  zm = z;

  //
  // First the Z rotation, then the X rotation.
  // (Ignore the Z coordinate)
  //

  for (i in 1:x.n)
  {
    for (j in 1:y.n)
    {
      // Do the first axis...
      tmp = rotz( [xm[i;j], ym[i;j], zm[i;j]], -theta*dtor );

      // Now do the second rotation...
      tmp = rotx ( tmp, -phi*dtor );

      xm[i;j] = tmp[1];
      ym[i;j] = tmp[2];
      zm[i;j] = tmp[3];
    }
  }

  //
  // Calculate the scale limits.
  //

  xmin = min (min (xm));
  xmax = max (max (xm));

  ymin = min (min (ym));
  ymax = max (max (ym));

  if (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].xmin[p] != 1j) { xmin = real (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].xmin[p]); }
  if (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].xmax[p] != 1j) { xmax = real (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].xmax[p]); }

  if (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].ymin[p] != 1j) { ymin = real (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].ymin[p]); }
  if (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].ymax[p] != 1j) { ymax = real (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].ymax[p]); }

  if (axis)
  {
    // Rotate the axes...
    tmp = rotz( xmax * axis3d[1;], -theta*dtor );
    ax1 = rotx( tmp, -phi*dtor );

    tmp = rotz( ymax * axis3d[2;], -theta*dtor );
    ax2 = rotx( tmp, -phi*dtor );

    tmp = rotz( ymax * axis3d[3;], -theta*dtor );
    ax3 = rotx( tmp, -phi*dtor );
  }

  //
  // Set up the window/viewport.
  //

  pgsls (1);
  pgsci (1);

  pgpage ();
  pgvstd ();

  if (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].aspect[p] != 0)
  {
    pgwnad (xmin, xmax, ymin, ymax);
  } else {
    pgswin (xmin, xmax, ymin, ymax);
  }

  //
  // Now, pick the tickmark properties, maybe reset
  // scale min/max values.
  //

  </ xmax; xmin; nxsub; xtick /> = ticks (xmin, xmax);
  </ ymax; ymin; nysub; ytick /> = ticks (ymin, ymax);

  pgsci (10);                        // Set back to normal grid, label color
  pgsch (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].char_height);       // Set character height
  pgbox (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].gridx[p], xtick, nxsub, ...
         PLOT_WINDOWS.[PLOT_ACTIVE_WIN].gridy[p], ytick, nysub);

  //
  // Plot the axes
  //

  if (axis)
  {
    pgline (2, [0, ax1[1]], [0, ax1[2]]);
    pgpt (1, ax1[1], ax1[2], 88);
    pgline (2, [0, ax2[1]], [0, ax2[2]]);
    pgpt (1, ax2[1], ax2[2], 89);
    pgline (2, [0, ax3[1]], [0, ax3[2]]);
    pgpt (1, ax3[1], ax3[2], 90);
  }

  pgsci (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].color[p;1]);       // Set line color.
  pgslw (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].width[p]);         // Line width.
  pgscf (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].font[p]);          // Set character font.
  pgsch (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].char_height);      // Set character height.
  pgsfs (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].fillstyle[p]);     // Set the fillstyle.

  //
  // Now, plot up the projected surface.
  //

  for (i in 1:xm.nr-1)
  {
    for (j in 1:xm.nc-1)
    {
      xpts = [ xm[i;j], xm[i;j+1], xm[i+1;j+1], xm[i+1;j] ];
      ypts = [ ym[i;j], ym[i;j+1], ym[i+1;j+1], ym[i+1;j] ];

      // Perform simple Backface Culling (hidden polygon removal).
      c = (xpts[1] - xpts[2]) * (ypts[3] - ypts[2]) - ...
      (ypts[1] - ypts[2]) * (xpts[3] - xpts[2]);
      if (c > 0)
      {
        pgpoly (4, xpts, ypts);
      }

    }
  }

  pgsls (1);
  pgsci (10);

  pglab(PLOT_WINDOWS.[PLOT_ACTIVE_WIN].xlabel[p], PLOT_WINDOWS.[PLOT_ACTIVE_WIN].ylabel[p], PLOT_WINDOWS.[PLOT_ACTIVE_WIN].title[p]);
  pgupdt ();
};

////////////
//
// Plot line-segments...
//

pllseg = function ( V )
{
  global (pi)

  check_plot_object ();
  p = mod (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].subplot, PLOT_WINDOWS.[PLOT_ACTIVE_WIN].nplot) + 1;	// The current index

  //
  // Error checking...
  //

  //
  // Plot scale limits...
  //

  xmin = min (V[;1]);
  xmax = max (V[;1]);

  ymin = min (V[;2]);
  ymax = max (V[;2]);

  if (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].xmin[p] != 1j) { xmin = real (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].xmin[p]); }
  if (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].xmax[p] != 1j) { xmax = real (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].xmax[p]); }

  if (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].ymin[p] != 1j) { ymin = real (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].ymin[p]); }
  if (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].ymax[p] != 1j) { ymax = real (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].ymax[p]); }

  pgslw (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].width[p]);         // Line width
  pgscf (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].font[p]);          // Set character font
  pgsch (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].char_height);      // Set character height

  pgsls (1);
  pgsci (1);

  pgenv (xmin, xmax, ymin, ymax, 0, -2);

  pgsci (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].color[p;1]);     // Set line color.
  pgsls (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].lstyle[p;1]);    // Set line style.
  pgslw (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].width[p;1]);     // Set line width.

  pgline (V.nr, V[;1], V[;2]);

  pgsls (1);
  pgsci (1);

  pglab(PLOT_WINDOWS.[PLOT_ACTIVE_WIN].xlabel[p], PLOT_WINDOWS.[PLOT_ACTIVE_WIN].ylabel[p], PLOT_WINDOWS.[PLOT_ACTIVE_WIN].title[p]);
  pgupdt ();
};

//---------------------------------------------------------------------
//
// Draw points in 3-D. pnts matrix consists of:
//
//    [ X , Y , Z ]
//
// A distinct point is drawn for each [X,Y,Z] triplet.
//

plpnt3 = function ( PNTS, theta, phi, axis, new )
{
  global (pi)

  check_plot_object ();
  p = mod (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].subplot, PLOT_WINDOWS.[PLOT_ACTIVE_WIN].nplot) + 1;	// The current index

  if (!exist (theta)) { theta = 60; }
  if (!exist (phi))   { phi = 35; }
  if (!exist (axis)) { axis = 0; }
  if (!exist (new)) { new = 1; }

  dtor = 2*pi / 360;

  //
  // First the Z rotation.
  // (Ignore the Z coordinate)
  //

  pnts = zeros (size (PNTS));

  for (i in 1:PNTS.nr)
  {
    // Do the first axis...
    tmp = rotz( PNTS[i;], -theta*dtor );

    // Now do the second rotation...
    tmp = rotx( tmp, -phi*dtor );

    // Load into the new pnts matrix...
    pnts[i;] = tmp';
  }

  //
  // Calculate the scale limits.
  //

  xmin = min (pnts[;1]);
  xmax = max (pnts[;1]);

  ymin = min (pnts[;2]);
  ymax = max (pnts[;2]);

  if (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].xmin[p] != 1j) { xmin = real (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].xmin[p]); }
  if (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].xmax[p] != 1j) { xmax = real (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].xmax[p]); }

  if (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].ymin[p] != 1j) { ymin = real (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].ymin[p]); }
  if (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].ymax[p] != 1j) { ymax = real (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].ymax[p]); }

  if (axis)
  {
    // Rotate the axes...
    tmp = rotz( xmax * axis3d[1;], -theta*dtor );
    ax1 = rotx( tmp, -phi*dtor );

    tmp = rotz( ymax * axis3d[2;], -theta*dtor );
    ax2 = rotx( tmp, -phi*dtor );

    tmp = rotz( ymax * axis3d[3;], -theta*dtor );
    ax3 = rotx( tmp, -phi*dtor );
  }

  //
  // Set up the window/viewport.
  //

  pgsls (1);
  pgsci (1);

  if (new)
  {
    pgpage ();
    pgvstd ();

    if (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].aspect[p] != 0)
    {
      pgwnad (xmin, xmax, ymin, ymax);
      } else {
      pgswin (xmin, xmax, ymin, ymax);
    }

    pgbox (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].gridx[p], 0, 0, PLOT_WINDOWS.[PLOT_ACTIVE_WIN].gridy[p], 0, 0);
  }

  //
  // Plot the axes
  //

  if (axis)
  {
    pgline (2, [0, ax1[1]], [0, ax1[2]]);
    pgpt (1, ax1[1], ax1[2], 88);
    pgline (2, [0, ax2[1]], [0, ax2[2]]);
    pgpt (1, ax2[1], ax2[2], 89);
    pgline (2, [0, ax3[1]], [0, ax3[2]]);
    pgpt (1, ax3[1], ax3[2], 90);
  }

  pgsci (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].color[p;1]);       // Set line color.
  pgslw (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].width[p]);         // Line width
  pgscf (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].font[p]);          // Set character font
  pgsch (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].char_height);      // Set character height

  //
  // Now, plot up the projected surface.
  //

  for (i in 1:pnts.nr)
  {
    pgpt (pnts.nr, pnts[;1], pnts[;2], PLOT_WINDOWS.[PLOT_ACTIVE_WIN].pstyle[p;1]);
    // pgpt (pnts.nr, pnts[;1], pnts[;2], 55);
  }

  pgsls (1);
  pgsci (1);

  pglab(PLOT_WINDOWS.[PLOT_ACTIVE_WIN].xlabel[p], PLOT_WINDOWS.[PLOT_ACTIVE_WIN].ylabel[p], PLOT_WINDOWS.[PLOT_ACTIVE_WIN].title[p]);
  pgupdt ();
};

//---------------------------------------------------------------------
//
// Draw points in 3-D. pnts matrix consists of:
//
//   [ X , Y , Z ]
//
// A line-segment is drawn for each pair of [X,Y,Z] triplets.
//
// PNTS:  Matrix of [X,Y,Z] triplets.
// theta: Defaults to 60 degrees. The angle (in degrees) of rotation
//        about the Z-axis.
// phi:   Defaults to 35 degrees. The angle (in degrees) of rotation
//        about the X-axis.
// axis:  If AXIS is TRUE (1), then the 3 orthogonal axes will
//        be drawn. Axis-1: [1, 0, 0], Axis-2: [0, 1, 0],
//        Axis-3: [0, 0, 1].  The axes are scaled to the plot window.
// new:
//

plln3 = function ( PNTS, theta, phi, axis, new )
{
  global (pi)

  check_plot_object ();
  p = mod (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].subplot, PLOT_WINDOWS.[PLOT_ACTIVE_WIN].nplot) + 1;	// The current index

  if (!exist (theta)) { theta = 60; }
  if (!exist (phi))   { phi = 35; }
  if (!exist (axis)) { axis = 0; }
  if (!exist (new)) { new = 1; }

  dtor = 2*pi / 360;

  //
  // First the Z rotation.
  // (Ignore the Z coordinate)
  //

  pnts = zeros (size (PNTS));

  for (i in 1:PNTS.nr)
  {
    // Do the first axis...
    tmp = rotz( PNTS[i;], -theta*dtor );

    // Now do the second rotation...
    tmp = rotx( tmp, -phi*dtor );

    // Load into the new pnts matrix...
    pnts[i;] = tmp';
  }

  //
  // Calculate the scale limits.
  //

  xmin = min (pnts[;1]);
  xmax = max (pnts[;1]);

  ymin = min (pnts[;2]);
  ymax = max (pnts[;2]);

  if (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].xmin[p] != 1j) { xmin = real (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].xmin[p]); }
  if (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].xmax[p] != 1j) { xmax = real (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].xmax[p]); }

  if (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].ymin[p] != 1j) { ymin = real (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].ymin[p]); }
  if (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].ymax[p] != 1j) { ymax = real (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].ymax[p]); }

  if (axis)
  {
    // Rotate the axes...
    tmp = rotz( xmax * axis3d[1;], -theta*dtor );
    ax1 = rotx( tmp, -phi*dtor );

    tmp = rotz( ymax * axis3d[2;], -theta*dtor );
    ax2 = rotx( tmp, -phi*dtor );

    tmp = rotz( ymax * axis3d[3;], -theta*dtor );
    ax3 = rotx( tmp, -phi*dtor );
  }

  //
  // Set up the window/viewport.
  //

  pgsls (1);
  pgsci (1);

  if (new)
  {
    // Clear the plot window, and start over.
    pgpage ();
    pgvstd ();

    if (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].aspect[p] != 0)
    {
      pgwnad (xmin, xmax, ymin, ymax);
      } else {
      pgswin (xmin, xmax, ymin, ymax);
    }

    pgbox (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].gridx[p], 0, 0, PLOT_WINDOWS.[PLOT_ACTIVE_WIN].gridy[p], 0, 0);
  }

  //
  // Plot the axes
  //

  if (axis)
  {
    pgline (2, [0, ax1[1]], [0, ax1[2]]);
    pgpt (1, ax1[1], ax1[2], 88);
    pgline (2, [0, ax2[1]], [0, ax2[2]]);
    pgpt (1, ax2[1], ax2[2], 89);
    pgline (2, [0, ax3[1]], [0, ax3[2]]);
    pgpt (1, ax3[1], ax3[2], 90);
  }

  pgsci (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].color[p;1]);       // Set line color.
  pgslw (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].width[p]);         // Line width
  pgscf (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].font[p]);          // Set character font
  pgsch (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].char_height);      // Set character height

  //
  // Now, plot up the projected line segments.
  //

  pgline (pnts.nr, pnts[;1], pnts[;2]);

  pgsls (1);
  pgsci (1);

  pglab(PLOT_WINDOWS.[PLOT_ACTIVE_WIN].xlabel[p], PLOT_WINDOWS.[PLOT_ACTIVE_WIN].ylabel[p], PLOT_WINDOWS.[PLOT_ACTIVE_WIN].title[p]);
  pgupdt ();
};

//---------------------------------------------------------------------
//
// Project a vector of 3D coordinates ( [ X, Y, Z ] ) into 2 Dimensions.
// This is accomplished via two axis rotations. One about Z (the first
// rotation), the second about X.
//
// THETA:  Z-axis rotation (degrees)
// PHI:    X-axis rotation (degrees)
//

project = function ( V, THETA, PHI )
{
  global (pi)
  dtor = 2*pi/360;

  //
  // Do Z-axis rotation...
  //

  for (i in 1:V.nr)
  {
    V[i;] = rotz( V[i;], -THETA*dtor )';
  }

  //
  // Do X-axis rotation...
  //

  for (i in 1:V.nr)
  {
    V[i;] = rotx( V[i;], -PHI*dtor )';
  }

  return V;
};

////////////
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

##############################################################################
#
# Set up the viewing altitude for 3-D plots
#

plalt = function ( ALT )
{
  check_plot_object ();
  i = mod (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].subplot, PLOT_WINDOWS.[PLOT_ACTIVE_WIN].nplot) + 1;

  if (exist (ALT))
  {
    PLOT_WINDOWS.[PLOT_ACTIVE_WIN].alt[i] = ALT;
    } else {
      PLOT_WINDOWS.[PLOT_ACTIVE_WIN].alt[i] = 60;
  }
  return PLOT_ACTIVE_WIN;
};

##############################################################################
#
# Set the viewing azimuth for 3-D plots
#

plaz = function ( AZ )
{
  check_plot_object ();
  i = mod (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].subplot, PLOT_WINDOWS.[PLOT_ACTIVE_WIN].nplot) + 1;

  if (exist (AZ))
  {
    PLOT_WINDOWS.[PLOT_ACTIVE_WIN].az[i] = AZ;
  }
  else
  {
    PLOT_WINDOWS.[PLOT_ACTIVE_WIN].az[i] = 45;
  }
  return PLOT_ACTIVE_WIN;
};


plxtics = function ( x, nx, scale )
{
  if (!exist (x))
  { x = 0; }
  if (!exist (nx))
  { nx = 0; }
  if (!exist (scale))
  { scale = 1; }

  check_plot_object ();

  p = mod (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].subplot, PLOT_WINDOWS.[PLOT_ACTIVE_WIN].nplot) + 1; // The current index

  PLOT_WINDOWS.[PLOT_ACTIVE_WIN].xtick.[p][1] = x[1];
  if (length(nx)>0)
  { PLOT_WINDOWS.[PLOT_ACTIVE_WIN].xtick.[p][5] = nx[1]; }

  PLOT_WINDOWS.[PLOT_ACTIVE_WIN].xtick_scale.[p] = scale;

  return PLOT_ACTIVE_WIN;
};

////////////
//
// Set 2-D grid styles. A not-so-friendly interface.
//
//
plxgrid = function ( sty_maj, sty_min )
{
  if (strlen(sty_maj)<1)
  {
    PLOT_WINDOWS.[PLOT_ACTIVE_WIN].gridx[i] = PLOT_GRID_X_DEFAULT;
    return 1;
  }

  check_plot_object ();

  i = mod (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].subplot, PLOT_WINDOWS.[PLOT_ACTIVE_WIN].nplot) + 1;

  // plot xgrid?
  if (!strindex(PLOT_WINDOWS.[PLOT_ACTIVE_WIN].gridx[i], "g"))
  {
    PLOT_WINDOWS.[PLOT_ACTIVE_WIN].gridx[i] = ...
        PLOT_WINDOWS.[PLOT_ACTIVE_WIN].gridx[i] + "g";
  }

  // [lpb, lt, lc, lw, pt, pc, ps, err, err_xy, err_lp];
  // default grid at major tics
  if (strindex(sty_maj,"with")<1)
  { sty_maj= "with lines " + sty_maj; }
  if (strindex(sty_maj,"lt")<1)
  { sty_maj = sty_maj + " lt 2"; }
  if (strindex(sty_maj,"lc")<1)
  { sty_maj = sty_maj + " lc rgb grey"; }
  _data = convert_gnuformat_to_plot_format( sty_maj ).plformat;
  PLOT_WINDOWS.[PLOT_ACTIVE_WIN].xtick.[i][2] = _data[3]; // color
  PLOT_WINDOWS.[PLOT_ACTIVE_WIN].xtick.[i][3] = _data[4]; // linewidth
  PLOT_WINDOWS.[PLOT_ACTIVE_WIN].xtick.[i][4] = _data[2]; // linetype

  if (strlen(sty_min) > 1)
  {
    // plot xgrid at minor tics?
    if (strindex(PLOT_WINDOWS.[PLOT_ACTIVE_WIN].gridx[i], "h")==0)
    {
      PLOT_WINDOWS.[PLOT_ACTIVE_WIN].gridx[i] = ...
          PLOT_WINDOWS.[PLOT_ACTIVE_WIN].gridx[i] + "h";
    }

    // default grid at minor tics
    if (strindex(sty_min,"with")<1)
    { sty_min= "with lines " + sty_min; }
    if (strindex(sty_min,"lt")<1)
    { sty_min= sty_min + " lt 3"; }
    if (strindex(sty_min,"lc")<1)
    { sty_min= sty_min + " lc rgb grey"; }
    _data = convert_gnuformat_to_plot_format( sty_min ).plformat;
    PLOT_WINDOWS.[PLOT_ACTIVE_WIN].xtick.[i][6] = _data[3]; // color
    PLOT_WINDOWS.[PLOT_ACTIVE_WIN].xtick.[i][7] = _data[4]; // linewidth
    PLOT_WINDOWS.[PLOT_ACTIVE_WIN].xtick.[i][8] = _data[2]; // linetype
  }

  return 0;
};

plx2grid = function ( sty_maj, sty_min )
{
  if (strlen(sty_maj)<1)
  {
    PLOT_WINDOWS.[PLOT_ACTIVE_WIN].alt_gridx[i] = PLOT_ALT_GRID_X_DEFAULT;
    return 1;
  }

  check_plot_object ();

  i = mod (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].subplot, PLOT_WINDOWS.[PLOT_ACTIVE_WIN].nplot) + 1;

  // plot xgrid?
  if (!strindex(PLOT_WINDOWS.[PLOT_ACTIVE_WIN].alt_gridx[i], "g"))
  {
    PLOT_WINDOWS.[PLOT_ACTIVE_WIN].alt_gridx[i] = ...
        PLOT_WINDOWS.[PLOT_ACTIVE_WIN].alt_gridx[i] + "g";
  }

  // [lpb, lt, lc, lw, pt, pc, ps, err, err_xy, err_lp];
  // default grid at major tics
  if (strindex(sty_maj,"with")<1)
  { sty_maj= "with lines " + sty_maj; }
  if (strindex(sty_maj,"lt")<1)
  { sty_maj = sty_maj + " lt 2"; }
  if (strindex(sty_maj,"lc")<1)
  { sty_maj = sty_maj + " lc rgb grey"; }
  _data = convert_gnuformat_to_plot_format( sty_maj ).plformat;
  PLOT_WINDOWS.[PLOT_ACTIVE_WIN].alt_xtick.[i][2] = _data[3]; // color
  PLOT_WINDOWS.[PLOT_ACTIVE_WIN].alt_xtick.[i][3] = _data[4]; // linewidth
  PLOT_WINDOWS.[PLOT_ACTIVE_WIN].alt_xtick.[i][4] = _data[2]; // linetype

  if (strlen(sty_min) > 1)
  {
    // plot xgrid at minor tics?
    if (strindex(PLOT_WINDOWS.[PLOT_ACTIVE_WIN].alt_gridx[i], "h")==0)
    {
      PLOT_WINDOWS.[PLOT_ACTIVE_WIN].alt_gridx[i] = ...
          PLOT_WINDOWS.[PLOT_ACTIVE_WIN].alt_gridx[i] + "h";
    }

    // default grid at minor tics
    if (strindex(sty_min,"with")<1)
    { sty_min= "with lines " + sty_min; }
    if (strindex(sty_min,"lt")<1)
    { sty_min= sty_min + " lt 3"; }
    if (strindex(sty_min,"lc")<1)
    { sty_min= sty_min + " lc rgb grey"; }
    _data = convert_gnuformat_to_plot_format( sty_min ).plformat;
    PLOT_WINDOWS.[PLOT_ACTIVE_WIN].alt_xtick.[i][6] = _data[3]; // color
    PLOT_WINDOWS.[PLOT_ACTIVE_WIN].alt_xtick.[i][7] = _data[4]; // linewidth
    PLOT_WINDOWS.[PLOT_ACTIVE_WIN].alt_xtick.[i][8] = _data[2]; // linetype
  }

  return 0;
};


plygrid = function ( sty_maj, sty_min )
{
  if (strlen(sty_maj)<1)
  {
    PLOT_WINDOWS.[PLOT_ACTIVE_WIN].gridy[i] = PLOT_GRID_Y_DEFAULT;
    return 1;
  }

  check_plot_object ();

  i = mod (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].subplot, PLOT_WINDOWS.[PLOT_ACTIVE_WIN].nplot) + 1;

  // plot xgrid at major tics?
  if (strindex(PLOT_WINDOWS.[PLOT_ACTIVE_WIN].gridy[i], "g")==0)
  {
    PLOT_WINDOWS.[PLOT_ACTIVE_WIN].gridy[i] = ...
        PLOT_WINDOWS.[PLOT_ACTIVE_WIN].gridy[i] + "g";
  }

  // [lpb, lt, lc, lw, pt, pc, ps, err, err_xy, err_lp];
  // default grid at major tics
  if (strindex(sty_maj,"with")<1)
  { sty_maj= "with lines " + sty_maj; }
  if (strindex(sty_maj,"lt")<1)
  { sty_maj = sty_maj + " lt 2"; }
  if (strindex(sty_maj,"lc")<1)
  { sty_maj = sty_maj + " lc rgb grey"; }

  _data = convert_gnuformat_to_plot_format( sty_maj ).plformat;
  PLOT_WINDOWS.[PLOT_ACTIVE_WIN].ytick.[i][2] = _data[3]; // color
  PLOT_WINDOWS.[PLOT_ACTIVE_WIN].ytick.[i][3] = _data[4]; // linewidth
  PLOT_WINDOWS.[PLOT_ACTIVE_WIN].ytick.[i][4] = _data[2]; // linetype

  if (strlen(sty_min) > 1)
  {
    // plot xgrid at minor tics?
    if (strindex(PLOT_WINDOWS.[PLOT_ACTIVE_WIN].gridy[i], "h")==0)
    {
      PLOT_WINDOWS.[PLOT_ACTIVE_WIN].gridy[i] = ...
          PLOT_WINDOWS.[PLOT_ACTIVE_WIN].gridy[i] + "h";
    }
    // default grid at minor tics
    if (strindex(sty_min,"with")<1)
    { sty_min= "with lines " + sty_min; }
    if (strindex(sty_min,"lt")<1)
    { sty_min= sty_min + " lt 3"; }
    if (strindex(sty_min,"lc")<1)
    { sty_min= sty_min + " lc rgb grey"; }

    // user wants something else?
    _data = convert_gnuformat_to_plot_format( sty_min ).plformat;
    PLOT_WINDOWS.[PLOT_ACTIVE_WIN].ytick.[i][6] = _data[3]; // color
    PLOT_WINDOWS.[PLOT_ACTIVE_WIN].ytick.[i][7] = _data[4]; // linewidth
    PLOT_WINDOWS.[PLOT_ACTIVE_WIN].ytick.[i][8] = _data[2]; // linetype
  }

  return 0;
};

ply2grid = function ( sty_maj, sty_min )
{
  if (strlen(sty_maj)<1)
  {
    PLOT_WINDOWS.[PLOT_ACTIVE_WIN].alt_gridy[i] = PLOT_ALT_GRID_Y_DEFAULT;
    return 1;
  }

  check_plot_object ();

  i = mod (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].subplot, PLOT_WINDOWS.[PLOT_ACTIVE_WIN].nplot) + 1;

  // plot xgrid?
  if (!strindex(PLOT_WINDOWS.[PLOT_ACTIVE_WIN].alt_gridy[i], "g"))
  {
    PLOT_WINDOWS.[PLOT_ACTIVE_WIN].alt_gridy[i] = ...
        PLOT_WINDOWS.[PLOT_ACTIVE_WIN].alt_gridy[i] + "g";
  }

  // [lpb, lt, lc, lw, pt, pc, ps, err, err_xy, err_lp];
  // default grid at major tics
  if (strindex(sty_maj,"with")<1)
  { sty_maj= "with lines " + sty_maj; }
  if (strindex(sty_maj,"lt")<1)
  { sty_maj = sty_maj + " lt 2"; }
  if (strindex(sty_maj,"lc")<1)
  { sty_maj = sty_maj + " lc rgb grey"; }
  _data = convert_gnuformat_to_plot_format( sty_maj ).plformat;
  PLOT_WINDOWS.[PLOT_ACTIVE_WIN].alt_ytick.[i][2] = _data[3]; // color
  PLOT_WINDOWS.[PLOT_ACTIVE_WIN].alt_ytick.[i][3] = _data[4]; // linewidth
  PLOT_WINDOWS.[PLOT_ACTIVE_WIN].alt_ytick.[i][4] = _data[2]; // linetype

  if (strlen(sty_min) > 1)
  {
    // plot xgrid at minor tics?
    if (strindex(PLOT_WINDOWS.[PLOT_ACTIVE_WIN].gridy[i], "h")==0)
    {
      PLOT_WINDOWS.[PLOT_ACTIVE_WIN].alt_gridy[i] = ...
          PLOT_WINDOWS.[PLOT_ACTIVE_WIN].alt_gridy[i] + "h";
    }
    // default grid at minor tics
    if (strindex(sty_min,"with")<1)
    { sty_min= "with lines " + sty_min; }
    if (strindex(sty_min,"lt")<1)
    { sty_min= sty_min + " lt 3"; }
    if (strindex(sty_min,"lc")<1)
    { sty_min= sty_min + " lc rgb grey"; }

    _data = convert_gnuformat_to_plot_format( sty_min ).plformat;
    PLOT_WINDOWS.[PLOT_ACTIVE_WIN].alt_ytick.[i][6] = _data[3]; // color
    PLOT_WINDOWS.[PLOT_ACTIVE_WIN].alt_ytick.[i][7] = _data[4]; // linewidth
    PLOT_WINDOWS.[PLOT_ACTIVE_WIN].alt_ytick.[i][8] = _data[2]; // linetype
  }

  return 0;
};

plx2tics = function ( x, nx, scale )
{
  check_plot_object ();
  p = mod (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].subplot, PLOT_WINDOWS.[PLOT_ACTIVE_WIN].nplot) + 1; // The current index

  if (!exist (x))
  { x = 0; }
  x = x[1];
  if (!exist (nx))
  { nx = 0; }
  nx = nx[1];
  if (!exist (scale))
  { scale = 1; }

  PLOT_WINDOWS.[PLOT_ACTIVE_WIN].alt_xtick.[p][1] = x[1];

  if (length(nx)>0)
  { PLOT_WINDOWS.[PLOT_ACTIVE_WIN].alt_xtick.[p][5] = nx[1]; }

  PLOT_WINDOWS.[PLOT_ACTIVE_WIN].alt_xtick_scale.[p] = scale;
  return PLOT_ACTIVE_WIN;
};

plytics = function ( x, nx, scale )
{
  check_plot_object ();
  p = mod (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].subplot, PLOT_WINDOWS.[PLOT_ACTIVE_WIN].nplot) + 1; // The current index

  if (!exist (x))
  { x = 0; }
  x = x[1];
  if (!exist (nx))
  { nx = 0; }
  nx = nx[1];
  if (!exist (scale))
  { scale = 1; }

  PLOT_WINDOWS.[PLOT_ACTIVE_WIN].ytick.[p][1] = x[1];

  if (length(nx)>0)
  { PLOT_WINDOWS.[PLOT_ACTIVE_WIN].ytick.[p][5] = nx[1]; }

  PLOT_WINDOWS.[PLOT_ACTIVE_WIN].ytick_scale.[p] = scale;
  return PLOT_ACTIVE_WIN;
};

ply2tics = function ( x, nx, scale )
{
  check_plot_object ();
  p = mod (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].subplot, PLOT_WINDOWS.[PLOT_ACTIVE_WIN].nplot) + 1; // The current index

  if (!exist (x))
  { x = 0; }
  x = x[1];
  if (!exist (nx))
  { nx = 0; }
  nx = nx[1];
  if (!exist (scale))
  { scale = 1; }

  PLOT_WINDOWS.[PLOT_ACTIVE_WIN].alt_ytick.[p][1] = x[1];

  if (length(nx)>0)
  { PLOT_WINDOWS.[PLOT_ACTIVE_WIN].alt_ytick.[p][5] = nx[1]; }

  PLOT_WINDOWS.[PLOT_ACTIVE_WIN].alt_ytick_scale.[p] = scale;
  return PLOT_ACTIVE_WIN;
};

plztics = function ( x, nx, scale )
{
  if (!exist (x))
  { x = 0; }
  if (!exist (nx))
  { nx = 0; }
  if (!exist (scale))
  { scale = 1; }

  check_plot_object ();

  p = mod (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].subplot, PLOT_WINDOWS.[PLOT_ACTIVE_WIN].nplot) + 1; // The current index

  PLOT_WINDOWS.[PLOT_ACTIVE_WIN].ztick.[p][1] = x[1];
  if (length(nx)>0)
  { PLOT_WINDOWS.[PLOT_ACTIVE_WIN].ztick.[p][2] = nx[1]; }

  PLOT_WINDOWS.[PLOT_ACTIVE_WIN].ztick_scale.[p] = scale;

  return PLOT_ACTIVE_WIN;
};


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// format data for plotting
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
plformat = function ( fmts )
{
  this_function = "plformat";

  if (!exist(fmts))
  { return 1; }
  if (class(fmts)!="string")
  { return 1; }

  if (!PLOT_ACTIVE_WIN)
  {
    plwins(1);
    plwin (1);
  }
  // The current index
  p = mod (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].subplot, PLOT_WINDOWS.[PLOT_ACTIVE_WIN].nplot) + 1;

  // line/symbol format for data
  _plformat = zeros(length(fmts),length(PLOT_DATASET_PLOT_FORMAT));
  _axes = zeros(1,length(fmts));
  _cols = cell (1,length(fmts));

  // contour parameters
  _zlevels = cell (1,length(fmts));
  _labels = cell (1,length(fmts));
  _ts = zeros(1,length(fmts));
  _sp = zeros(length(fmts),2);

  for (i in 1:length(fmts))
  {
    rval = convert_gnuformat_to_plot_format( fmts[i] );
    if (exist(rval.plformat))
    { _plformat[i;] = rval.plformat; }

    if (_plformat[i;1]==3)
    {
      if (imag(PLOT_WINDOWS.[PLOT_ACTIVE_WIN].zmin[p])==0)
      {
        _plformat[i;2] = real(PLOT_WINDOWS.[PLOT_ACTIVE_WIN].zmin[p]);
      }
      if (imag(PLOT_WINDOWS.[PLOT_ACTIVE_WIN].zmax[p]) == 0)
      {
        _plformat[i;3] = real(PLOT_WINDOWS.[PLOT_ACTIVE_WIN].zmax[p]);
      }
    }
    if (exist(rval.axes))
    { _axes[i] = rval.axes; }
    if (exist(rval.use_cols))
    {
      _cols[i] = rval.use_cols;
    }
    else
    {
      _cols[i] = []; // default is to use column 1 and 2
    }

    // for contour plots zlevels may be provided
    if (exist(rval.contour_zlevels))
    { _zlevels[i] = rval.contour_zlevels; }
    else
    { _zlevels[i] = []; }
    if (exist(rval.contour_labels))
    {
      _labels[i] = rval.contour_labels;
      if (exist(rval.contour_ts))
      { _ts[i] = rval.contour_ts; }
      if (exist(rval.contour_sp))
      { _sp[i;] = rval.contour_sp; }
    }
    else
    { _labels[i] = []; }
  }

  // find all lines that have [-1,-1] for lt and lc
  // and change that to increasing type and color
  _il = find((_plformat[;1] < 3) && (_plformat[;2]==-1) && (_plformat[;3]==-1));
  if (length(_il) > 0)
  {
    _lt = PLOT_DATASET_PLOT_FORMAT[2];
    _lc = PLOT_DATASET_PLOT_FORMAT[3];
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
  _il = find((_plformat[;1] < 3) && (_plformat[;5]==-1) && (_plformat[;6]==-1));
  if (length(_il) > 0)
  {
    _pt = PLOT_DATASET_PLOT_FORMAT[5];
    _pc = PLOT_DATASET_PLOT_FORMAT[6];
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

  PLOT_WINDOWS.[PLOT_ACTIVE_WIN].plformat.[p] = _plformat;
  PLOT_WINDOWS.[PLOT_ACTIVE_WIN].axes.[p] = _axes;
  PLOT_WINDOWS.[PLOT_ACTIVE_WIN].cols.[p] = _cols;
  PLOT_WINDOWS.[PLOT_ACTIVE_WIN].contour_zlevels.[p] = _zlevels;
  PLOT_WINDOWS.[PLOT_ACTIVE_WIN].contour_labels.[p] = _labels;
  PLOT_WINDOWS.[PLOT_ACTIVE_WIN].contour_ts.[p] = _ts;
  PLOT_WINDOWS.[PLOT_ACTIVE_WIN].contour_sp.[p] = _sp;

  return 0;
};

plscale = function ( scale_x, scale_y )
{
  check_plot_object ();
  i = mod (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].subplot, PLOT_WINDOWS.[PLOT_ACTIVE_WIN].nplot) + 1;

  if (exist(scale_x))
  {
    if (class(scale_x)=="string")
    {
      if (strindex(tolower(scale_x),"log")>0)
      {
        PLOT_WINDOWS.[PLOT_ACTIVE_WIN].gridx[i] = PLOT_GRID_X_DEFAULT_LOG;
      }
      else if (strindex(tolower(scale_x),"lin")>0)
      {
        PLOT_WINDOWS.[PLOT_ACTIVE_WIN].gridx[i] = PLOT_GRID_X_DEFAULT;
      }
      else
      {
        PLOT_WINDOWS.[PLOT_ACTIVE_WIN].gridx[i] = scale_x;
      }
    }
  }
  else
  {
    PLOT_WINDOWS.[PLOT_ACTIVE_WIN].gridx[i] = PLOT_GRID_X_DEFAULT;
  }

  if (exist(scale_y))
  {
    if (class(scale_y)=="string")
    {
      if (strindex(tolower(scale_y),"log")>0)
      {
        PLOT_WINDOWS.[PLOT_ACTIVE_WIN].gridy[i] = PLOT_GRID_Y_DEFAULT_LOG;
      }
      else if (strindex(tolower(scale_y),"lin")>0)
      {
        PLOT_WINDOWS.[PLOT_ACTIVE_WIN].gridy[i] = PLOT_GRID_Y_DEFAULT;
      }
      else
      {
        PLOT_WINDOWS.[PLOT_ACTIVE_WIN].gridy[i] = scale_y;
      }
    }
  }
  else
  {
    PLOT_WINDOWS.[PLOT_ACTIVE_WIN].gridy[i] = PLOT_GRID_Y_DEFAULT;
  }

  return 0;
};

plscale2 = function ( scale_x, scale_y )
{
  check_plot_object ();
  i = mod (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].subplot, PLOT_WINDOWS.[PLOT_ACTIVE_WIN].nplot) + 1;

  PLOT_WINDOWS.[PLOT_ACTIVE_WIN].alt_gridx[i] = PLOT_ALT_GRID_X_DEFAULT;
  if (exist(scale_x))
  {
    if (strindex(tolower(scale_x),"log"))
    { PLOT_WINDOWS.[PLOT_ACTIVE_WIN].alt_gridx[i] = PLOT_ALT_GRID_X_DEFAULT_LOG; }
  }

  PLOT_WINDOWS.[PLOT_ACTIVE_WIN].alt_gridy[i] = PLOT_ALT_GRID_Y_DEFAULT;
  if (exist(scale_y))
  {
    if (strindex(tolower(scale_y),"log"))
    { PLOT_WINDOWS.[PLOT_ACTIVE_WIN].alt_gridy[i] = PLOT_ALT_GRID_Y_DEFAULT_LOG; }
  }

  return 0;
};


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Place some text on the plot
//
pltext = function (s, loc, vel, opts )
{
  check_plot_object ();
  i = mod (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].subplot, PLOT_WINDOWS.[PLOT_ACTIVE_WIN].nplot) + 1;

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
  lc   = 1;
  cs = "w";

  if (exist(opts))
  {
    if (exist(opts.angle))
    {
      incl[2] = incl[1] * tan(opts.angle * 3.1415/180);
    }
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

    if (class(opts.color)=="num")
    {
      lc = opts.color;
    }
    else if (class(opts.color)=="string")
    {
      if (strindex(opts.color, "rgb")==1)
      {
        _c = strsplt(opts.color, "'BLANK");
        if (length(_c) > 1)
        {
          _c = gsub("'",_c[2]).string;
          _idx_c = find(PLOT_COLOR_RGB_NAME==_c);
          if (length(_idx_c)==1)
          {
            lc = _idx_c;
          }
        }
      }
    }

    // coordinate system in which x,y are given
    if (class(opts.cs)=="string")
    {
      if (strindex(opts.cs,"v") == 1)
      {
        cs = "v"; // viewport
      }
      else if (strindex(opts.cs,"n") == 1)
      {
        cs = "n"; // normalized device coordinates
      }
    }
  }

  dx = incl[1];
  dy = incl[2];

  rval = <<>>;
  rval.x = x;
  rval.y = y;
  rval.dx = dx;
  rval.dy = dy;
  rval.scale = vel;
  rval.just  = just;
  rval.color = lc;
  rval.text  = s;
  rval.cs    = cs;

  // putting stuff into the mde
  n_idx = length(PLOT_WINDOWS.[PLOT_ACTIVE_WIN].pltext.[i]);
  new_idx = num2str(n_idx+1, "%05.0f");
  PLOT_WINDOWS.[PLOT_ACTIVE_WIN].pltext.[i].[new_idx] = rval;
  return 0;
};



//
// from data list, column cell-array and axes array figure out
//  plot limits for first and second set of coordinates
//
static(find_list_xy_scales);
find_list_xy_scales = function (data, cols, axes, p)
{
  if (!length(data))
  { error ("plot: Can't plot empty list"); }

  m = members(data);

  xmin = [];
  xmax = [];
  idx_x1 = find((axes == PLOT_AXES_REV.x1y1)||(axes == PLOT_AXES_REV.x1y2));
  if (!isempty(idx_x1))
  {
    // do the x-axis
    if (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].xmin[p] != 1j)
    {
      xmin = real (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].xmin[p]);
    } else {
      for (i in range(idx_x1))
      {
        m_i = m[idx_x1[i]];
        xmin = min(xmin, min(data.[m_i][;cols[idx_x1[i]][1]]));
      }
    }

    // do the x-axis
    if (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].xmax[p] != 1j)
    {
      xmax = real (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].xmax[p]);
    } else {
      for (i in range(idx_x1))
      {
        m_i = m[idx_x1[i]];
        xmax = max(xmax, max(data.[m_i][;cols[idx_x1[i]][1]]));
      }
    }
  }

  ymin = [];
  ymax = [];
  idx_y1 = find((axes == PLOT_AXES_REV.x1y1)||(axes == PLOT_AXES_REV.x2y1));
  if (!isempty(idx_y1))
  {
    // do the y-axis
    if (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].ymin[p] != 1j)
    {
      ymin = real (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].ymin[p]);
    } else {
      for (i in range(idx_y1))
      {
        m_i = m[idx_y1[i]];
        ymin = min(ymin, min(data.[m_i][;cols[idx_y1[i]][2]]));
      }
    }

    // do the y-axis
    if (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].ymax[p] != 1j)
    {
      ymax = real (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].ymax[p]);
    } else {
      for (i in range(idx_y1))
      {
        m_i = m[idx_y1[i]];
        ymax = max(ymax, max(data.[m_i][;cols[idx_y1[i]][2]]));
      }
    }
  }

  alt_ymin = [];
  alt_ymax = [];
  idx_y2 = find((axes == PLOT_AXES_REV.x1y2)||(axes == PLOT_AXES_REV.x2y2));
  if (!isempty(idx_y2))
  {
    // do the y-axis
    if (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].alt_ymin[p] != 1j)
    {
      alt_ymin = real (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].alt_ymin[p]);
    } else {
      for (i in range(idx_y2))
      {
        m_i = m[idx_y2[i]];
        alt_ymin = min(alt_ymin, min(data.[m_i][;cols[idx_y2[i]][2]]));
      }
    }

    // do the y-axis
    if (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].alt_ymax[p] != 1j)
    {
      alt_ymax = real (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].alt_ymax[p]);
    } else {
      for (i in range(idx_y2))
      {
        m_i = m[idx_y2[i]];
        alt_ymax = max(alt_ymax, max(data.[m_i][;cols[idx_y2[i]][2]]));
      }
    }
  }

  alt_xmin = [];
  alt_xmax = [];
  idx_x2 = find((axes == PLOT_AXES_REV.x2y1)||(axes == PLOT_AXES_REV.x2y2));
  if (!isempty(idx_x2))
  {
    // do the x-axis
    if (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].alt_xmin[p] != 1j)
    {
      alt_xmin = real (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].alt_xmin[p]);
    } else {
      for (i in range(idx_x2))
      {
        m_i = m[idx_x2[i]];
        alt_xmin = min(alt_xmin, min(data.[m_i][;cols[idx_x2[i]][1]]));
      }
    }

    // do the x-axis
    if (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].alt_xmax[p] != 1j)
    {
      alt_xmax = real (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].alt_xmax[p]);
    } else {
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



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Plot data matrix or list as a scatter or line plot, or histogram (barplot)
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
plot = function ( data, export_fn )
{
  this_function = "plot";

  _zlevels = []; // for contour plotting

  if (!PLOT_ACTIVE_WIN)
  {
    plwins(1);
    plwin (1);
  }

  // The current index
  p = mod (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].subplot, PLOT_WINDOWS.[PLOT_ACTIVE_WIN].nplot) + 1;

  //
  // Draw single data set 3d special graph:
  //  pm3d, contour
  //
  if (is3dlist(data))
  {
    if (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].plformat.[p].nr==0)
    {
      printf(this_function + ": Don't know how to plot the dataset. Cannot continue !\n");
      return PLOT_ACTIVE_WIN;
    }

    // process x range
    _x = <<>>;
    if ( (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].xmin[p]==1j) ...
            || (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].xmax[p]==1j) )
    {
      _d = minmax(data.x);
      _x = <<min=_d[1];max=_d[2]>>;
    }
    if (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].xmin[p]!=1j)
    { _x.min = PLOT_WINDOWS.[PLOT_ACTIVE_WIN].xmin[p]; }
    if (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].xmax[p]!=1j)
    { _x.max = PLOT_WINDOWS.[PLOT_ACTIVE_WIN].xmax[p]; }

    // process y range
    _y = <<>>;
    if ( (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].ymin[p]==1j) ...
            || (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].ymax[p]==1j) )
    {
      _d = minmax(data.y);
      _y = <<min=_d[1];max=_d[2]>>;
    }
    if (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].ymin[p]!=1j)
    { _y.min = PLOT_WINDOWS.[PLOT_ACTIVE_WIN].ymin[p]; }
    if (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].ymax[p]!=1j)
    { _y.max = PLOT_WINDOWS.[PLOT_ACTIVE_WIN].ymax[p]; }

    // process z range for pm3d or contour plots
    _z = <<>>;
    if ( (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].zmin[p]==1j) ...
          || (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].zmax[p]==1j) )
    {
      _d = minmax(data.z);
      _z = <<min=min(_d[1;]);max=max(_d[2;])>>;
    }
    if (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].zmin[p]!=1j)
    { _z.min = real(PLOT_WINDOWS.[PLOT_ACTIVE_WIN].zmin[p]); }
    if (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].zmax[p]!=1j)
    { _z.max = real(PLOT_WINDOWS.[PLOT_ACTIVE_WIN].zmax[p]); }

    // logs
    _logx = 0;
    if(strindex(PLOT_WINDOWS.[PLOT_ACTIVE_WIN].gridx[p],"l")>0)
    {
      _logx = 1;
      _x.min = log10(_x.min);
      _x.max = log10(_x.max);
    }
    _logy = 0;
    if(strindex(PLOT_WINDOWS.[PLOT_ACTIVE_WIN].gridy[p],"l")>0)
    {
      _logy = 1;
      _y.min = log10(_y.min);
      _y.max = log10(_y.max);
    }

    // do we resize tic marks?
    scale = 1;
    if (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].xtick_scale.[p][1]!=1)
    { scale = PLOT_WINDOWS.[PLOT_ACTIVE_WIN].xtick_scale.[p][1]; }
    if (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].ytick_scale.[p][1]!=1)
    { scale = PLOT_WINDOWS.[PLOT_ACTIVE_WIN].ytick_scale.[p][1]; }

    pgsch (THIS_LIBRARY_DEFAULT_CHAR_SIZE * scale);
    pgswin(_x.min, _x.max, _y.min, _y.max);

    _plformat = PLOT_WINDOWS.[PLOT_ACTIVE_WIN].plformat.[p];
    _plformat[2] = _z.min;
    _plformat[3] = _z.max;
    _pgprintf(data, _plformat, _logx, _logy, ...
        PLOT_WINDOWS.[PLOT_ACTIVE_WIN].contour_zlevels.[p][1], ...
            PLOT_WINDOWS.[PLOT_ACTIVE_WIN].contour_labels.[p][1], ...
                PLOT_WINDOWS.[PLOT_ACTIVE_WIN].contour_ts.[p][1] );

    pgbox (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].gridx[p] + "c", ...
        PLOT_WINDOWS.[PLOT_ACTIVE_WIN].xtick.[p][1:4], ...
            PLOT_WINDOWS.[PLOT_ACTIVE_WIN].xtick.[p][5:8], ...
                PLOT_WINDOWS.[PLOT_ACTIVE_WIN].gridy[p] + "c", ...
                    PLOT_WINDOWS.[PLOT_ACTIVE_WIN].ytick.[p][1:4], ...
                        PLOT_WINDOWS.[PLOT_ACTIVE_WIN].ytick.[p][5:8]);

    pgsls (1);
    pgsci (1);
    pgslw (1);
  }
  else if ((class(data) == "num") || ishist(data))
  {
    //
    // Draw the graph
    //  Step through the matrix plotting
    //  each column versus the 1st
    //
    // plot a single matrix or a special list: histogram
    if (ishist(data))
    { data = hist_line(data); }

    if (data.nc == 1)
    { data = [(1:data.nr)', data]; }

    //
    // is there any plformat info on how to plot it
    //
    if (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].plformat.[p].nr==0)
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
      _fmt = PLOT_DATASET_PLOT_FORMAT;
      _lt = PLOT_DATASET_PLOT_FORMAT[2];
      _lc = PLOT_DATASET_PLOT_FORMAT[3];
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
    }
    else
    {
      // user knows what she is doing:
      //  assume that length of plformat determines number of plots
      //  and go along with that assumption
      _plformat = PLOT_WINDOWS.[PLOT_ACTIVE_WIN].plformat.[p];
      _axes = PLOT_WINDOWS.[PLOT_ACTIVE_WIN].axes.[p];
      _cols = PLOT_WINDOWS.[PLOT_ACTIVE_WIN].cols.[p];
      if (all(_axes == PLOT_AXES_REV.x1y1))
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
            }
            else if ((_plformat[_i;8]==1)||(_plformat[_i;8]==3))
            {
              // errorx or errory
              _cols[_i] = [1,min(_i+1,data.nc),min(_i+2,data.nc),min(_i+3,data.nc)];
            }
            else
            {
              // errorxy
                _cols[_i] = [1,min(_i+1,data.nc),min(_i+2,data.nc),min(_i+3,data.nc),...
                    min(_i+4,data.nc),min(_i+5,data.nc)];
            }
          }
          _xcol_idx = unique([_xcol_idx, _cols[_i][1]]);
          _ycol_idx = unique([_ycol_idx, _cols[_i][2]]);
        }
        _idx_x2 = _idx_x1;
        _idx_y2 = _idx_y1;
        _alt_xcol_idx = _xcol_idx;
        _alt_ycol_idx = _ycol_idx;

      }
      else
      {
        // user has multiple axis: then everything has to be specified for
        // each plotted dataset
        _idx_x1 = find((_axes == PLOT_AXES_REV.x1y1) || (_axes == PLOT_AXES_REV.x1y2));
        _idx_y1 = find((_axes == PLOT_AXES_REV.x1y1) || (_axes == PLOT_AXES_REV.x2y1));
        _idx_x2 = find((_axes == PLOT_AXES_REV.x2y1) || (_axes == PLOT_AXES_REV.x2y2));
        _idx_y2 = find((_axes == PLOT_AXES_REV.x1y2) || (_axes == PLOT_AXES_REV.x2y2));
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
    _idx_x1y2 = [];
    _idx_x2y1 = [];
    _idx_x2y2 = [];
    _idx_x1y1 = find(_axes == PLOT_AXES_REV.x1y1);
    if (length(_idx_x1y1) < length(_axes))
    {
      _idx_x1y2 = find(_axes == PLOT_AXES_REV.x1y2);
      if (length(_idx_x1y1)+length(_idx_x1y2) < length(_axes))
      {
        _idx_x2y1 = find(_axes == PLOT_AXES_REV.x2y1);
        if (length(_idx_x1y1)+length(_idx_x1y2)+length(_idx_x2y1) < length(_axes))
        {
          _idx_x2y2 = find(_axes == PLOT_AXES_REV.x2y1);
          if (length(_idx_x1y1)+length(_idx_x1y2)+length(_idx_x2y1)+length(_idx_x2y2) != length(_axes))
          {
            printf("What is going on? Help! I need somebody! Help! Just anybody! Help!\n");
            return 1;
          }
        }
      }
    }

    if (!isempty(_idx_x2y2))
    {
      // process x range
      _x = <<>>;
      if ( (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].alt_xmin[p]==1j) ...
            || (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].alt_xmax[p]==1j) )
      { _x = find_scales ("x", data, _alt_xcol_idx, p); }
      if (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].alt_xmin[p]!=1j)
      { _x.min = PLOT_WINDOWS.[PLOT_ACTIVE_WIN].alt_xmin[p]; }
      if (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].alt_xmax[p]!=1j)
      { _x.max = PLOT_WINDOWS.[PLOT_ACTIVE_WIN].alt_xmax[p]; }
      // process y range
      _y = <<>>;
      if ( (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].alt_ymin[p]==1j) ...
            || (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].alt_ymax[p]==1j) )
      { _y = find_scales ("y", data, _alt_ycol_idx, p); }
      if (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].alt_ymin[p]!=1j)
      { _y.min = PLOT_WINDOWS.[PLOT_ACTIVE_WIN].alt_ymin[p]; }
      if (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].alt_ymax[p]!=1j)
      { _y.max = PLOT_WINDOWS.[PLOT_ACTIVE_WIN].alt_ymax[p]; }
      // logs
      _logx = 0;
      if(strindex(PLOT_WINDOWS.[PLOT_ACTIVE_WIN].alt_gridx[p],"l")>0)
      {
        _logx = 1;
        _x.min = log10(_x.min);
        _x.max = log10(_x.max);
      }
      _logy = 0;
      if(strindex(PLOT_WINDOWS.[PLOT_ACTIVE_WIN].alt_gridy[p],"l")>0)
      {
        _logy = 1;
        _y.min = log10(_y.min);
        _y.max = log10(_y.max);
      }

      // do we use x1? if not, plot the opposite edge of frame now
      _idx_x1 = find((_axes == PLOT_AXES_REV.x1y1) || (_axes == PLOT_AXES_REV.x1y2));
      _fmt_x = "";
      if (length(_idx_x1)==0)
      { _fmt_x = "c"; }
      // do we use y1? if not, plot the opposite edge of frame now
      _fmt_y = "";
      _idx_y1 = find((_axes == PLOT_AXES_REV.x1y1) || (_axes == PLOT_AXES_REV.x2y1));
      if (length(_idx_y1)==0)
      { _fmt_y = "c"; }

      scale = 1;
      if (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].alt_xtick_scale.[p][1]!=1)
      { scale = PLOT_WINDOWS.[PLOT_ACTIVE_WIN].alt_xtick_scale.[p][1]; }
      if (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].alt_ytick_scale.[p][1]!=1)
      { scale = PLOT_WINDOWS.[PLOT_ACTIVE_WIN].alt_ytick_scale.[p][1]; }

      pgsch (THIS_LIBRARY_DEFAULT_CHAR_SIZE * scale);
      pgswin (_x.min, _x.max, _y.min, _y.max);
      pgbox (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].alt_gridx[p] + _fmt_x, ...
          PLOT_WINDOWS.[PLOT_ACTIVE_WIN].alt_xtick.[p][1], ...
              PLOT_WINDOWS.[PLOT_ACTIVE_WIN].alt_xtick.[p][2:5], ...
                  PLOT_WINDOWS.[PLOT_ACTIVE_WIN].alt_gridy[p] + _fmt_y, ...
                      PLOT_WINDOWS.[PLOT_ACTIVE_WIN].alt_ytick.[p][1], ...
                          PLOT_WINDOWS.[PLOT_ACTIVE_WIN].alt_ytick.[p][2;5]);

      for (_i in range(_idx_x2y2))
      { _pgprintf(data, _plformat[_idx_x2y2[_i];], _logx, _logy, _cols[_idx_x2y2[_i]]); }
      pgsls (1);
      pgsci (1);
      pgslw (1);
    }

    if (!isempty(_idx_x2y1))
    {
      // process x range
      _x = <<>>;
      if ( (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].alt_xmin[p]==1j) ...
            || (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].alt_xmax[p]==1j) )
      { _x = find_scales ("x", data, _alt_xcol_idx, p); }
      if (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].alt_xmin[p]!=1j)
      { _x.min = PLOT_WINDOWS.[PLOT_ACTIVE_WIN].alt_xmin[p]; }
      if (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].alt_xmax[p]!=1j)
      { _x.max = PLOT_WINDOWS.[PLOT_ACTIVE_WIN].alt_xmax[p]; }
      // process y range
      _y = <<>>;
      if ( (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].ymin[p]==1j) ...
            || (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].ymax[p]==1j) )
      { _y = find_scales ("y", data, _ycol_idx, p); }
      if (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].ymin[p]!=1j)
      { _y.min = PLOT_WINDOWS.[PLOT_ACTIVE_WIN].ymin[p]; }
      if (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].ymax[p]!=1j)
      { _y.max = PLOT_WINDOWS.[PLOT_ACTIVE_WIN].ymax[p]; }
      // logs
      _logx = 0;
      if(strindex(PLOT_WINDOWS.[PLOT_ACTIVE_WIN].alt_gridx[p],"l")>0)
      {
        _logx = 1;
        _x.min = log10(_x.min);
        _x.max = log10(_x.max);
      }
      _logy = 0;
      if(strindex(PLOT_WINDOWS.[PLOT_ACTIVE_WIN].gridy[p],"l")>0)
      {
        _logy = 1;
        _y.min = log10(_y.min);
        _y.max = log10(_y.max);
      }

      // do we use x1? if not, plot the opposite edge of frame now
      _idx_x1 = find((_axes == PLOT_AXES_REV.x1y1) || (_axes == PLOT_AXES_REV.x1y2));
      _fmt_x = "";
      if (length(_idx_x1)==0)
      { _fmt_x = "c"; }
      // do we use y2? if not, plot the opposite edge of frame now
      _fmt_y = "";
      _idx_y2 = find((_axes == PLOT_AXES_REV.x1y2) || (_axes == PLOT_AXES_REV.x2y2));
      if (length(_idx_y2)==0)
      { _fmt_y = "c"; }

      scale = 1;
      if (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].alt_xtick_scale.[p][1]!=1)
      { scale = PLOT_WINDOWS.[PLOT_ACTIVE_WIN].alt_xtick_scale.[p][1]; }
      if (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].ytick_scale.[p][1]!=1)
      { scale = PLOT_WINDOWS.[PLOT_ACTIVE_WIN].ytick_scale.[p][1]; }

      pgsch (THIS_LIBRARY_DEFAULT_CHAR_SIZE * scale);
      pgswin (_x.min, _x.max, _y.min, _y.max);
      pgbox (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].alt_gridx[p] + _fmt_x, ...
          PLOT_WINDOWS.[PLOT_ACTIVE_WIN].alt_xtick.[p][1], ...
              PLOT_WINDOWS.[PLOT_ACTIVE_WIN].alt_xtick.[p][2:5], ...
                  PLOT_WINDOWS.[PLOT_ACTIVE_WIN].gridy[p] + _fmt_y, ...
                      PLOT_WINDOWS.[PLOT_ACTIVE_WIN].ytick.[p][1], ...
                          PLOT_WINDOWS.[PLOT_ACTIVE_WIN].ytick.[p][2:5]);

      for (_i in range(_idx_x2y1))
      { _pgprintf(data, _plformat[_idx_x1y1[_i];], _logx, _logy, _cols[_idx_x2y1[_i]]); }
      pgsls (1);
      pgsci (1);
      pgslw (1);
    }

    if (!isempty(_idx_x1y2))
    {
      // process x range
      _x = <<>>;
      if ( (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].xmin[p]==1j) ...
            || (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].xmax[p]==1j) )
      { _x = find_scales ("x", data, _xcol_idx, p); }
      if (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].xmin[p]!=1j)
      { _x.min = PLOT_WINDOWS.[PLOT_ACTIVE_WIN].xmin[p]; }
      if (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].xmax[p]!=1j)
      { _x.max = PLOT_WINDOWS.[PLOT_ACTIVE_WIN].xmax[p]; }
      // process y range
      _y = <<>>;
      if ( (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].alt_ymin[p]==1j) ...
            || (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].alt_ymax[p]==1j) )
      { _y = find_scales ("y", data, _alt_ycol_idx, p); }
      if (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].alt_ymin[p]!=1j)
      { _y.min = PLOT_WINDOWS.[PLOT_ACTIVE_WIN].alt_ymin[p]; }
      if (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].alt_ymax[p]!=1j)
      { _y.max = PLOT_WINDOWS.[PLOT_ACTIVE_WIN].alt_ymax[p]; }

      _logx = 0;
      if(strindex(PLOT_WINDOWS.[PLOT_ACTIVE_WIN].gridx[p],"l")>0)
      {
        _logx = 1;
        _x.min = log10(_x.min);
        _x.max = log10(_x.max);
      }
      _logy = 0;
      if(strindex(PLOT_WINDOWS.[PLOT_ACTIVE_WIN].alt_gridy[p],"l")>0)
      {
        _logy = 1;
        _y.min = log10(_y.min);
        _y.max = log10(_y.max);
      }

      // do we use x1? if not, plot the opposite edge of frame now
      _idx_x2 = find((_axes == PLOT_AXES_REV.x2y1) || (_axes == PLOT_AXES_REV.x2y2));
      _fmt_x = "";
      if (length(_idx_x2)==0)
      { _fmt_x = "c"; }
      // do we use y1? if not, plot the opposite edge of frame now
      _fmt_y = "";
      _idx_y1 = find((_axes == PLOT_AXES_REV.x1y1) || (_axes == PLOT_AXES_REV.x2y1));
      if (length(_idx_y1)==0)
      { _fmt_y = "c"; }

      scale = 1;
      if (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].xtick_scale.[p][1]!=1)
      { scale = PLOT_WINDOWS.[PLOT_ACTIVE_WIN].xtick_scale.[p][1]; }
      if (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].alt_ytick_scale.[p][1]!=1)
      { scale = PLOT_WINDOWS.[PLOT_ACTIVE_WIN].alt_ytick_scale.[p][1]; }

      pgsch (THIS_LIBRARY_DEFAULT_CHAR_SIZE * scale);
      pgswin (_x.min, _x.max, _y.min, _y.max);
      pgbox (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].gridx[p] + _fmt_x, ...
          PLOT_WINDOWS.[PLOT_ACTIVE_WIN].xtick.[p][1:4], ...
              PLOT_WINDOWS.[PLOT_ACTIVE_WIN].xtick.[p][5:8], ...
                  PLOT_WINDOWS.[PLOT_ACTIVE_WIN].alt_gridy[p] + _fmt_y, ...
                      PLOT_WINDOWS.[PLOT_ACTIVE_WIN].alt_ytick.[p][1:4], ...
                          PLOT_WINDOWS.[PLOT_ACTIVE_WIN].alt_ytick.[p][5:8]);

      for (_i in range(_idx_x1y2))
      { _pgprintf(data, _plformat[_idx_x1y2[_i];], _logx, _logy, _cols[_idx_x1y2[_i]]); }
      pgsls (1);
      pgsci (1);
      pgslw (1);
    }

    if (!isempty(_idx_x1y1))
    {
      // process x range
      _x = <<>>;
      if ( (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].xmin[p]==1j) ...
            || (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].xmax[p]==1j) )
      { _x = find_scales ("x", data, _xcol_idx, p); }
      if (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].xmin[p]!=1j)
      { _x.min = PLOT_WINDOWS.[PLOT_ACTIVE_WIN].xmin[p]; }
      if (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].xmax[p]!=1j)
      { _x.max = PLOT_WINDOWS.[PLOT_ACTIVE_WIN].xmax[p]; }
      // process y range

      _y = <<>>;
      if ( (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].ymin[p]==1j) ...
            || (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].ymax[p]==1j) )
      { _y = find_scales ("y", data, _ycol_idx, p); }
      if (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].ymin[p]!=1j)
      { _y.min = PLOT_WINDOWS.[PLOT_ACTIVE_WIN].ymin[p]; }
      if (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].ymax[p]!=1j)
      { _y.max = PLOT_WINDOWS.[PLOT_ACTIVE_WIN].ymax[p]; }
      // logs
      _logx = 0;
      if(strindex(PLOT_WINDOWS.[PLOT_ACTIVE_WIN].gridx[p],"l")>0)
      {
        _logx = 1;
        _x.min = log10(_x.min);
        _x.max = log10(_x.max);
      }
      _logy = 0;
      if(strindex(PLOT_WINDOWS.[PLOT_ACTIVE_WIN].gridy[p],"l")>0)
      {
        _logy = 1;
        _y.min = log10(_y.min);
        _y.max = log10(_y.max);
      }

      // do we use x2? if not, plot the opposite edge of frame now
      _idx_x2 = find((_axes == PLOT_AXES_REV.x2y1) || (_axes == PLOT_AXES_REV.x2y2));
      _fmt_x = "";
      if (length(_idx_x2)==0)
      { _fmt_x = "c"; }
      // do we use y2? if not, plot the opposite edge of frame now
      _fmt_y = "";
      _idx_y2 = find((_axes == PLOT_AXES_REV.x1y2) || (_axes == PLOT_AXES_REV.x2y2));
      if (length(_idx_y2)==0)
      { _fmt_y = "c"; }

      // do we resize tic marks?
      scale = 1;
      if (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].xtick_scale.[p][1]!=1)
      { scale = PLOT_WINDOWS.[PLOT_ACTIVE_WIN].xtick_scale.[p][1]; }
      if (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].ytick_scale.[p][1]!=1)
      { scale = PLOT_WINDOWS.[PLOT_ACTIVE_WIN].ytick_scale.[p][1]; }

      pgsch (THIS_LIBRARY_DEFAULT_CHAR_SIZE * scale);
      pgswin(_x.min, _x.max, _y.min, _y.max);
      pgbox (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].gridx[p] + _fmt_x, ...
          PLOT_WINDOWS.[PLOT_ACTIVE_WIN].xtick.[p][1:4], ...
              PLOT_WINDOWS.[PLOT_ACTIVE_WIN].xtick.[p][5:8], ...
                  PLOT_WINDOWS.[PLOT_ACTIVE_WIN].gridy[p] + _fmt_y, ...
                      PLOT_WINDOWS.[PLOT_ACTIVE_WIN].ytick.[p][1:4], ...
                          PLOT_WINDOWS.[PLOT_ACTIVE_WIN].ytick.[p][5:8]);

      for (_i in range(_idx_x1y1))
      {
        _pgprintf(data, _plformat[_idx_x1y1[_i];], _logx, _logy, _cols[_idx_x1y1[_i]]);
      }
      pgsls (1);
      pgsci (1);
      pgslw (1);
    }

  }
  else if (class (data) == "list")
  {
    if (length(data)==0)
    { error (this_function + ": How would you like me to plot an empty list?\n"); }

    //
    // is there any plformat info on how to plot it
    //
    if (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].plformat.[p].nr==0)
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
      _fmt = PLOT_DATASET_PLOT_FORMAT;
      _plformat = _fmt;
      _lt = PLOT_DATASET_PLOT_FORMAT[2];
      _lc = PLOT_DATASET_PLOT_FORMAT[3];
      for (i in 2:length(data))
      {
        _lc++;
        _fmt[2] = _lt;
        _fmt[3] = _lc;
        if (_lc == PLOT_COLOR_COUNT_DEFAULT)
        {
          _lc = 1;
          _lt = _lt + 1;
          if (_lt > 8)
          { _lt = 1; }
        }
        _plformat = [_plformat; _fmt ];
      }
    }
    else
    {
      // user knows what she is doing:
      //  assume that length of plformat determines number of plots
      //  and go along with that assumption
      _plformat = PLOT_WINDOWS.[PLOT_ACTIVE_WIN].plformat.[p];
      _zlevels = PLOT_WINDOWS.[PLOT_ACTIVE_WIN].contour_zlevels.[p];
      _labels = PLOT_WINDOWS.[PLOT_ACTIVE_WIN].contour_labels.[p];
      _ts = PLOT_WINDOWS.[PLOT_ACTIVE_WIN].contour_ts.[p];
      _axes = PLOT_WINDOWS.[PLOT_ACTIVE_WIN].axes.[p];
      _cols = PLOT_WINDOWS.[PLOT_ACTIVE_WIN].cols.[p];
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
    _idx_x2y2 = [];
    _idx_x2y1 = [];
    _idx_x1y2 = [];
    _idx_x1y1 = find(_axes == PLOT_AXES_REV.x1y1);
    if (length(_idx_x1y1) < length(_axes))
    {
      _idx_x1y2 = find(_axes == PLOT_AXES_REV.x1y2);
      if (length(_idx_x1y1)+length(_idx_x1y2) < length(_axes))
      {
        _idx_x2y1 = find(_axes == PLOT_AXES_REV.x2y1);
        if (length(_idx_x1y1)+length(_idx_x1y2)+length(_idx_x2y1) < length(_axes))
        {
          _idx_x2y2 = find(_axes == PLOT_AXES_REV.x2y2);
          if (length(_idx_x1y1)+length(_idx_x1y2)+length(_idx_x2y1)+length(_idx_x2y2) != length(_axes))
          {
            printf("What is going on? Help! I need somebody! Help! Just anybody! Help!\n");
            return 0;
          }
        }
      }
    }
    _idx_x1 = union(_idx_x1y2, _idx_x1y1);
    _idx_x2 = union(_idx_x2y2, _idx_x2y1);
    _idx_y1 = union(_idx_x1y1, _idx_x2y1);
    _idx_y2 = union(_idx_x1y2, _idx_x2y2);


    _s = find_list_xy_scales(data, _cols, _axes, p);

    if (!isempty(_idx_x2y1))
    {
      _xmin = _s.alt_xmin;
      _xmax = _s.alt_xmax;
      _ymin = _s.ymin;
      _ymax = _s.ymax;

      // logs
      _logx = 0;
      if(strindex(PLOT_WINDOWS.[PLOT_ACTIVE_WIN].alt_gridx[p],"l")>0)
      {
        _logx = 1;
        _xmin = log10(_xmin);
        _xmax = log10(_xmax);
      }
      _logy = 0;
      if(strindex(PLOT_WINDOWS.[PLOT_ACTIVE_WIN].gridy[p],"l")>0)
      {
        _logy = 1;
        _ymin = log10(_ymin);
        _ymax = log10(_ymax);
      }

      // do we use x1? if not, plot the opposite edge of frame now
      _fmt_x = "";
      if (length(_idx_x1)==0)
      { _fmt_x = "c"; }
      // do we use y2? if not, plot the opposite edge of frame now
      _fmt_y = "";
      if (length(_idx_y2)==0)
      { _fmt_y = "c"; }

      scale = 1;
      if (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].alt_xtick_scale.[p][1]!=1)
      { scale = PLOT_WINDOWS.[PLOT_ACTIVE_WIN].alt_xtick_scale.[p][1]; }
      if (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].ytick_scale.[p][1]!=1)
      { scale = PLOT_WINDOWS.[PLOT_ACTIVE_WIN].ytick_scale.[p][1]; }

      pgsch (THIS_LIBRARY_DEFAULT_CHAR_SIZE * scale);
      pgswin (_xmin, _xmax, _ymin, _ymax);
      pgbox (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].alt_gridx[p] + _fmt_x, ...
          PLOT_WINDOWS.[PLOT_ACTIVE_WIN].alt_xtick.[p][1:4], ...
              PLOT_WINDOWS.[PLOT_ACTIVE_WIN].alt_xtick.[p][5:8], ...
                  PLOT_WINDOWS.[PLOT_ACTIVE_WIN].gridy[p] + _fmt_y, ...
                      PLOT_WINDOWS.[PLOT_ACTIVE_WIN].ytick.[p][1:4], ...
                          PLOT_WINDOWS.[PLOT_ACTIVE_WIN].ytick.[p][5:8]);

      for (_i in range(_idx_x2y1))
      {
        _m_i = _m[_idx_x2y1[_i]];
        if (ishist(data.[_m_i]))
        {
          _pgprintf(hist_line(data.[_m_i]), _plformat[_idx_x2y1[_i];], 0, 0, _cols[_idx_x2y1[_i]]);
        }
        else if (is3dlist(data.[_m_i]))
        {
          _pgprintf(data.[_m_i], _plformat[_idx_x2y1[_i];], _logx, _logy, ...
              PLOT_WINDOWS.[PLOT_ACTIVE_WIN].contour_zlevels.[p][_idx_x2y1[_i]], ...
                  PLOT_WINDOWS.[PLOT_ACTIVE_WIN].contour_labels.[p][_idx_x2y1[_i]], ...
                      PLOT_WINDOWS.[PLOT_ACTIVE_WIN].contour_ts.[p][_idx_x2y1[_i]] );
        }
        else
        {
          _pgprintf(data.[_m_i], _plformat[_idx_x2y1[_i];], 0, 0, _cols[_idx_x2y1[_i]]);
        }
      }
    }

    if (!isempty(_idx_x2y2))
    {
      _xmin = _s.alt_xmin;
      _xmax = _s.alt_xmax;
      _ymin = _s.alt_ymin;
      _ymax = _s.alt_ymax;

      // logs
      _logx = 0;
      if(strindex(PLOT_WINDOWS.[PLOT_ACTIVE_WIN].alt_gridx[p],"l")>0)
      {
        _logx = 1;
        _xmin = log10(_xmin);
        _xmax = log10(_xmax);
      }
      _logy = 0;
      if(strindex(PLOT_WINDOWS.[PLOT_ACTIVE_WIN].alt_gridy[p],"l")>0)
      {
        _logy = 1;
        _ymin = log10(_ymin);
        _ymax = log10(_ymax);
      }

      // do we use x1? if not, plot the opposite edge of frame now
      _fmt_x = "";
      if (length(_idx_x1)==0)
      { _fmt_x = "c"; }
      // do we use y1? if not, plot the opposite edge of frame now
      _fmt_y = "";
      if (length(_idx_y1)==0)
      { _fmt_y = "c"; }

      scale = 1;
      if (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].alt_xtick_scale.[p][1]!=1)
      { scale = PLOT_WINDOWS.[PLOT_ACTIVE_WIN].alt_xtick_scale.[p][1]; }
      if (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].alt_ytick_scale.[p][1]!=1)
      { scale = PLOT_WINDOWS.[PLOT_ACTIVE_WIN].alt_ytick_scale.[p][1]; }

      pgsch (THIS_LIBRARY_DEFAULT_CHAR_SIZE * scale);
      pgswin (_xmin, _xmax, _ymin, _ymax);
      pgbox (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].alt_gridx[p] + _fmt_x, ...
          PLOT_WINDOWS.[PLOT_ACTIVE_WIN].alt_xtick.[p][1:4], ...
              PLOT_WINDOWS.[PLOT_ACTIVE_WIN].alt_xtick.[p][5:8], ...
                  PLOT_WINDOWS.[PLOT_ACTIVE_WIN].alt_gridy[p] + _fmt_y, ...
                      PLOT_WINDOWS.[PLOT_ACTIVE_WIN].alt_ytick.[p][1:4], ...
                          PLOT_WINDOWS.[PLOT_ACTIVE_WIN].alt_ytick.[p][5:8]);

      for (_i in range(_idx_x2y2))
      {
        _m_i = _m[_idx_x2y2[_i]];
        if (ishist(data.[_m_i]))
        {
          _pgprintf(hist_line(data.[_m_i]), _plformat[_idx_x2y2[_i];], 0, 0, _cols[_idx_x2y2[_i]]);
        }
        else if (is3dlist(data.[_m_i]))
        {
          _pgprintf(data.[_m_i], _plformat[_idx_x2y2[_i];], _logx, _logy, ...
              PLOT_WINDOWS.[PLOT_ACTIVE_WIN].contour_zlevels.[p][_idx_x2y2[_i]], ...
                  PLOT_WINDOWS.[PLOT_ACTIVE_WIN].contour_labels.[p][_idx_x2y2[_i]], ...
                      PLOT_WINDOWS.[PLOT_ACTIVE_WIN].contour_ts.[p][_idx_x2y2[_i]] );
        }
        else
        {
          _pgprintf(data.[_m_i], _plformat[_idx_x2y2[_i];], 0, 0, _cols[_idx_x2y2[_i]]);
        }
      }
    }

    if (!isempty(_idx_x1y2))
    {
      _xmin = _s.xmin;
      _xmax = _s.xmax;
      _ymin = _s.alt_ymin;
      _ymax = _s.alt_ymax;

      // logs
      _logx = 0;
      if(strindex(PLOT_WINDOWS.[PLOT_ACTIVE_WIN].gridx[p],"l")>0)
      {
        _logx = 1;
        _xmin = log10(_xmin);
        _xmax = log10(_xmax);
      }
      _logy = 0;
      if(strindex(PLOT_WINDOWS.[PLOT_ACTIVE_WIN].alt_gridy[p],"l")>0)
      {
        _logy = 1;
        _ymin = log10(_ymin);
        _ymax = log10(_ymax);
      }

      // do we use x1? if not, plot the opposite edge of frame now
      _fmt_x = "";
      if (length(_idx_x2)==0)
      { _fmt_x = "c"; }
      // do we use y1? if not, plot the opposite edge of frame now
      _fmt_y = "";
      if (length(_idx_y1)==0)
      { _fmt_y = "c"; }

      scale = 1;
      if (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].xtick_scale.[p][1]!=1)
      { scale = PLOT_WINDOWS.[PLOT_ACTIVE_WIN].xtick_scale.[p][1]; }
      if (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].alt_ytick_scale.[p][1]!=1)
      { scale = PLOT_WINDOWS.[PLOT_ACTIVE_WIN].alt_ytick_scale.[p][1]; }

      pgsch (THIS_LIBRARY_DEFAULT_CHAR_SIZE * scale);
      pgswin (_xmin, _xmax, _ymin, _ymax);

      pgbox (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].gridx[p] + _fmt_x, ...
          PLOT_WINDOWS.[PLOT_ACTIVE_WIN].xtick.[p][1:4], ...
              PLOT_WINDOWS.[PLOT_ACTIVE_WIN].xtick.[p][5:8], ...
                  PLOT_WINDOWS.[PLOT_ACTIVE_WIN].alt_gridy[p] + _fmt_y, ...
                      PLOT_WINDOWS.[PLOT_ACTIVE_WIN].alt_ytick.[p][1:4], ...
                          PLOT_WINDOWS.[PLOT_ACTIVE_WIN].alt_ytick.[p][5:8]);

      for (_i in range(_idx_x1y2))
      {
        _m_i = _m[_idx_x1y2[_i]];
        if (ishist(data.[_m_i]))
        {
          _pgprintf(hist_line(data.[_m_i]),_plformat[_idx_x1y2[_i];], _logx, _logy, _cols[_idx_x1y2[_i]]);
        }
        else if (is3dlist(data.[_m_i]))
        {
          _pgprintf(data.[_m_i], _plformat[_idx_x1y2[_i];], _logx, _logy, ...
              PLOT_WINDOWS.[PLOT_ACTIVE_WIN].contour_zlevels.[p][_idx_x1y2[_i]], ...
                  PLOT_WINDOWS.[PLOT_ACTIVE_WIN].contour_labels.[p][_idx_x1y2[_i]], ...
                      PLOT_WINDOWS.[PLOT_ACTIVE_WIN].contour_ts.[p][_idx_x1y2[_i]] );
        }
        else
        {
          _pgprintf(data.[_m_i], _plformat[_idx_x1y2[_i];], _logx, _logy, _cols[_idx_x1y2[_i]]);
        }
      }
      pgsls (PLOT_DEFAULT_LINE_STYLE);
      pgsci (1);
      pgslw (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].width[p]);
    }

    if (!isempty(_idx_x1y1))
    {
      _xmin = _s.xmin;
      _xmax = _s.xmax;
      _ymin = _s.ymin;
      _ymax = _s.ymax;

      // logs
      _logx = 0;
      if(strindex(PLOT_WINDOWS.[PLOT_ACTIVE_WIN].gridx[p],"l")>0)
      {
        _logx = 1;
        _xmin = log10(_xmin);
        _xmax = log10(_xmax);
      }
      _logy = 0;
      if(strindex(PLOT_WINDOWS.[PLOT_ACTIVE_WIN].gridy[p],"l")>0)
      {
        _logy = 1;
        _ymin = log10(_ymin);
        _ymax = log10(_ymax);
      }

      // do we use x2? if not, plot the opposite edge of frame now
      _fmt_x = "";
      if (length(_idx_x2)==0)
      { _fmt_x = "c"; }
      // do we use y2? if not, plot the opposite edge of frame now
      _fmt_y = "";
      if (length(_idx_y2)==0)
      { _fmt_y = "c"; }

      scale = 1;
      if (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].xtick_scale.[p][1]!=1)
      { scale = PLOT_WINDOWS.[PLOT_ACTIVE_WIN].xtick_scale.[p][1]; }
      if (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].ytick_scale.[p][1]!=1)
      { scale = PLOT_WINDOWS.[PLOT_ACTIVE_WIN].ytick_scale.[p][1]; }

      pgsch (THIS_LIBRARY_DEFAULT_CHAR_SIZE * scale);
      pgswin (_xmin, _xmax, _ymin, _ymax);

      // first plot pm3d's
      for (_i in range(_idx_x1y1))
      {
        _m_i = _m[_idx_x1y1[_i]];
        if (is3dlist(data.[_m_i]) && (_plformat[_idx_x1y1[_i];1]==3))
        {
          _pgprintf(data.[_m_i], _plformat[ _idx_x1y1[_i]; ], _logx, _logy, ...
              _zlevels[ _idx_x1y1[_i] ], _labels[ _idx_x1y1[_i] ] );
        }
      }

      pgbox (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].gridx[p] + _fmt_x, ...
          PLOT_WINDOWS.[PLOT_ACTIVE_WIN].xtick.[p][1:4], ...
              PLOT_WINDOWS.[PLOT_ACTIVE_WIN].xtick.[p][5:8], ...
                  PLOT_WINDOWS.[PLOT_ACTIVE_WIN].gridy[p] + _fmt_y, ...
                      PLOT_WINDOWS.[PLOT_ACTIVE_WIN].ytick.[p][1:4], ...
                          PLOT_WINDOWS.[PLOT_ACTIVE_WIN].ytick.[p][5:8]);

      for (_i in range(_idx_x1y1))
      {
        _m_i = _m[_idx_x1y1[_i]];
        if (ishist(data.[_m_i]))
        {
          _pgprintf(hist_line(data.[_m_i]), _plformat[_idx_x1y1[_i];], _logx, _logy, _cols[_idx_x1y1[_i]]);
        }
        else if (is3dlist(data.[_m_i]) && (_plformat[_idx_x1y1[_i];1]!=3) )
        {
          _pgprintf(data.[_m_i], _plformat[_idx_x1y1[_i];], _logx, _logy, ...
              _zlevels[_idx_x1y1[_i]], _labels[_idx_x1y1[_i]], _ts[_idx_x1y1[_i]] );
        }
        else
        {
          _pgprintf(data.[_m_i], _plformat[_idx_x1y1[_i];], _logx, _logy, _cols[_idx_x1y1[_i]]);
        }
      }

      pgsls (PLOT_DEFAULT_LINE_STYLE);
      pgsci (1);
      pgslw (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].width[p]);
    }



  }
  else
  {
    error ("plot: un-acceptable argument");
  }

  //
  // Now do the legend
  //
  if (length(PLOT_WINDOWS.[PLOT_ACTIVE_WIN].desc.[p])>0)
  {
    desc = PLOT_WINDOWS.[PLOT_ACTIVE_WIN].desc.[p];
    idx_n = find(strlen(desc)>0);
    if (length(idx_n)>0)
    {
      desc_pos = PLOT_WINDOWS.[PLOT_ACTIVE_WIN].desc_pos.[p];
      desc_pos_xy = PLOT_WINDOWS.[PLOT_ACTIVE_WIN].desc_pos_xy.[p];
      desc_scale = PLOT_WINDOWS.[PLOT_ACTIVE_WIN].desc_scale.[p];
      _pglegend(desc[idx_n], _plformat[idx_n;], desc_scale, desc_pos, desc_pos_xy,_zlevels[idx_n]);
    }
  }

  pgsls (PLOT_DEFAULT_LINE_STYLE);
  pgsci (1);
  pgslw (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].width[p]);

  // Do the axes labels in the end:
  if (strlen(PLOT_WINDOWS.[PLOT_ACTIVE_WIN].xlabel[p])>0)
  {
    _d = PLOT_WINDOWS.[PLOT_ACTIVE_WIN].xlabel_scale.[p];
    pgsch (THIS_LIBRARY_DEFAULT_CHAR_SIZE * _d[1]);
    pgmtxt("B", _d[2], _d[3], _d[4], PLOT_WINDOWS.[PLOT_ACTIVE_WIN].xlabel[p]);
  }
  if (strlen(PLOT_WINDOWS.[PLOT_ACTIVE_WIN].ylabel[p])>0)
  {
    _d = PLOT_WINDOWS.[PLOT_ACTIVE_WIN].ylabel_scale.[p];
    pgsch (THIS_LIBRARY_DEFAULT_CHAR_SIZE * _d[1]);
    pgmtxt("L", _d[2], _d[3], _d[4], PLOT_WINDOWS.[PLOT_ACTIVE_WIN].ylabel[p]);
  }

  // secondary x-label and the title are the same?
  if (strlen(PLOT_WINDOWS.[PLOT_ACTIVE_WIN].alt_xlabel[p])>0)
  {
    _d = PLOT_WINDOWS.[PLOT_ACTIVE_WIN].alt_xlabel_scale.[p];
    pgsch (THIS_LIBRARY_DEFAULT_CHAR_SIZE * _d[1]);
    pgmtxt( "T", _d[2], _d[3], _d[4], PLOT_WINDOWS.[PLOT_ACTIVE_WIN].alt_xlabel[p] );
  }
  else if (strlen(PLOT_WINDOWS.[PLOT_ACTIVE_WIN].title[p])>0)
  {
    _d = PLOT_WINDOWS.[PLOT_ACTIVE_WIN].title_scale.[p];
    pgsch (THIS_LIBRARY_DEFAULT_CHAR_SIZE * _d[1]);
    pgmtxt("T", _d[2], 0.5, 0.5, PLOT_WINDOWS.[PLOT_ACTIVE_WIN].title[p]);
  }

  // secondary y-label
  if (strlen(PLOT_WINDOWS.[PLOT_ACTIVE_WIN].alt_ylabel[p])>0)
  {
    _d = PLOT_WINDOWS.[PLOT_ACTIVE_WIN].alt_ylabel_scale.[p];
    pgsch (THIS_LIBRARY_DEFAULT_CHAR_SIZE * _d[1]);
    pgmtxt( "R", _d[2], _d[3], _d[4], PLOT_WINDOWS.[PLOT_ACTIVE_WIN].alt_ylabel[p] );
  }

  // graffitti?
  if (length(PLOT_WINDOWS.[PLOT_ACTIVE_WIN].pltext.[p])>0)
  {
    for (_i in members(PLOT_WINDOWS.[PLOT_ACTIVE_WIN].pltext.[p]))
    {
      _s = PLOT_WINDOWS.[PLOT_ACTIVE_WIN].pltext.[p].[_i];
      pgsch (THIS_LIBRARY_DEFAULT_CHAR_SIZE * _s.scale);
      pgsci (_s.color);

      if (_s.cs == "w")
      {
        // world coordinate system - respect log and such
        if(strindex(PLOT_WINDOWS.[PLOT_ACTIVE_WIN].gridx[p],"l")>0)
        {
          _x  = log10(_s.x);
        }
        else
        {
          _x  = _s.x;
        }
        _logy = 0;
        if(strindex(PLOT_WINDOWS.[PLOT_ACTIVE_WIN].gridy[p],"l")>0)
        {
          _y  = log10(_s.y);
        }
        else
        {
          _y  = _s.y;
        }
      }
      else if (_s.cs == "n")
      {
        // normalized device coordinates
        _vp = pgqvp(0);   // [u1,v1,u2.v2], viewport
        _wo = pgqwin();   // [x1,x2,y1.y2], world
        // x:
        if(strindex(PLOT_WINDOWS.[PLOT_ACTIVE_WIN].gridx[p],"l")>0)
        {
          _x = _wo[1] + (_wo[2] - _wo[1])./log10(_vp[2]/_vp[1]).*log10(_s.x/_vp[1]);
        }
        else
        {
          _x = _wo[1] + (_wo[2] - _wo[1])./(_vp[2] - _vp[1]) .*  (_s.x - _vp[1]);
        }
        // y:
        if(strindex(PLOT_WINDOWS.[PLOT_ACTIVE_WIN].gridy[p],"l")>0)
        {
          _y = _wo[3] + (_wo[4] - _wo[3])./log10(_vp[4]/_vp[3]).*log10(_s.y/_vp[3]);
        }
        else
        {
          _y = _wo[3] + (_wo[4] - _wo[3])./(_vp[4] - _vp[3]) .*  (_s.y - _vp[3]);
        }
      }

      pgptxt(_x, _y, atan2(_s.dy,_s.dx) * 180/3.1415, _s.just, _s.text);
      pgsci (1);
      pgsch (1);
    }
  }

  // Increment the plot no. so that next time
  if (PLOT_WINDOWS.[PLOT_ACTIVE_WIN].nplot > 1)
  {
    PLOT_WINDOWS.[PLOT_ACTIVE_WIN].subplot = PLOT_WINDOWS.[PLOT_ACTIVE_WIN].subplot + 1;
  }

  if (strlen(export_fn)>1)
  {
    // make a printout
    act_win_save = PLOT_ACTIVE_WIN;
    // 1. get viewport in normalized device coordinates of the device to be copied
    desc_vp = pgqvp(0); // xl,xr,yb,yt
    // 2. open new file device
    plcopy(act_win_save, export_fn);
    // 3. set viewport of the new file device to be the same as the original
    plwin (plwin(),desc_vp);
    // 4. call thyself once more
    plot( data );
    // 55. close thyself
    plclose();
    // become active plot again
    plwin(act_win_save,desc_vp,0);
  }

  reset_active_plot_object();

  return PLOT_ACTIVE_WIN;
};









