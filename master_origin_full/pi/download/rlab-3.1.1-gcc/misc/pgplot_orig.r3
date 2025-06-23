//
// New plot.r for use with PGPLOT library.
// The help files for these functions are in
// misc/plhelp
//

// pgplot.r

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

static (PGPLOT_PLOT_WINDOWS)		// The static plot window structure
if (!exist (PGPLOT_PLOT_WINDOWS))
{ PGPLOT_PLOT_WINDOWS = <<>>; }

static (PGPLOT_ACTIVE_WIN);   // The active/current plot window
if (!exist(PGPLOT_ACTIVE_WIN))
{ PGPLOT_ACTIVE_WIN = 0; }

//
// Maintain the transformations for 3-D plots.
//

static (basex, basey, height)
basex = 2;
basey = 2;
height = 4;

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
static (check_3d_list);
static (find_char);
static (get_style);
static (make_legend);

//
// Defaults
//

static (grid_x_default, grid_y_default)
static (grid_3x_default, grid_3y_default, grid_3z_default)

grid_x_default = "bcnst";
grid_y_default = "bcnstv";
grid_3x_default = "bnstu";
grid_3y_default = "bnstu";
grid_3z_default = "bcdmnstuv";

static (subplot_spec)          // Pay attention to the specifications ?
static (subplot_horizontal)    // The horizontal specification.
static (subplot_vertical)      // The vertical specification.

subplot_spec = 0;

static (axis3d)

axis3d = [ 1, 0, 0;
           0, 1, 0;
           0, 0, 1];
//
// Initial windows size/config...
//

static (win_xp_df, win_yp_df);
static (win_xleng_df, win_yleng_df, win_wid, WIN_XLENG_DF, WIN_YLENG_DF);
static (win_xoff_df, win_yoff_df);
WIN_XLENG_DF = 560;
WIN_YLENG_DF = 400;

win_xp_df = 0;
win_yp_df = 0;
win_xleng_df = WIN_XLENG_DF;   // default width  of the plot window
win_yleng_df = WIN_YLENG_DF;   // default height of the plot window
win_xoff_df = 100;
win_yoff_df = 100;
win_wid = 80;         // default width of the konsole in X-windows
                      // used to scale all the widths

//
// Create the default plot-object.
// Initialize to all the default values
//


create_plot_object = function ( N, nx, ny )
{
  if (!exist (N))
  { N = 1; }

  pobj = <<>>;

  pobj.window.xp    = win_xp_df;      // Number of X pixels
  pobj.window.yp    = win_yp_df;      // Number of Y pixels
  pobj.window.xleng = win_xleng_df;   // Page length, X
  pobj.window.yleng = win_yleng_df;   // Page length, Y
  pobj.window.xoff  = win_xoff_df;    // Page offset, X
  pobj.window.yoff  = win_yoff_df;    // Page offset, Y
  pobj.subplot      = 0;        // The current subplot no.
  pobj.nplot        = nx*ny;    // Total no. of plots on window
  pobj.nx           = nx;
  pobj.ny           = ny;

  pobj.fontld = 0;              // Loaded extended fonts?

  pobj.char_height = 1.5;       // Character height.

  pobj.device = "/XWIN";

  for (i in 1:(nx*ny))
  {
    pobj.style.[i] = "line";    // The type/style of plot to draw
    pobj.nbin[i]   = 1j;          // The number of bins for a histogram
    pobj.font[i]   = 1;		        // The current font
    pobj.height[i] = 1;		      // The character height
    pobj.xlabel[i] = "";
    pobj.ylabel[i] = "";
    pobj.zlabel[i] = "";
    pobj.title[i]  = "";
    pobj.orientation[i] = "portrait";
    pobj.desc.[i]  = "default";		      // The legend description
    pobj.gridx[i]  = grid_x_default;	  // Plot axes style, 2D-X
    pobj.gridy[i]  = grid_y_default;	  // Plot axes style, 2D-Y
    pobj.grid3x[i] = grid_3x_default;	  // Plot axes style, 3D-X
    pobj.grid3y[i] = grid_3y_default;	  // Plot axes style, 3D-Y
    pobj.grid3z[i] = grid_3z_default;	  // Plot axes style, 3D-Z
    pobj.aspect[i] = 0;		              // Plot aspect style
    pobj.alt[i]    = 60;
    pobj.az[i]     = 45;

    pobj.xmin[i] = 1j;
    pobj.xmax[i] = 1j;
    pobj.ymin[i] = 1j;
    pobj.ymax[i] = 1j;
    pobj.zmin[i] = 1j;
    pobj.zmax[i] = 1j;

    // load default values
    pobj.page.xp = win_xp_df;
    pobj.page.yp = win_yp_df;
    pobj.xleng   = WIN_XLENG_DF;
    pobj.yleng   = WIN_YLENG_DF;
    pobj.xoff    = win_xoff_df;
    pobj.yoff    = win_yoff_df;

    pobj.color[i;]    = 1:14;               // 14 possible colors...
    pobj.lstyle[i;]   = [1:5,1:5,1:5,1:5];  // 5 possible line styles...
    pobj.pstyle[i;]   = 1:8;                // 8 possible point styles...
    pobj.width[i;]    = ones(1,32);         // 32 possible point line widths.
    pobj.fillstyle[i] = 1;                  // Default to solid.
  }

  //
  // Save the newly generated plot-object
  // in a list of plot-objects.
  //
  PGPLOT_PLOT_WINDOWS.[N] = pobj;
};

////////////////////////////////////////////////////////////////////////
//
// (Re) Set the initial plot window specs...

plwin_init = function (xp, yp, xleng, yleng, xoff, yoff)
{
  if (!exist (xp))
  {
    win_xp_df = 0;
  } else {
    win_xp_df = xp;
  }

  if (!exist (yp))
  {
    win_yp_df = 0;
  } else {
    win_yp_df = yp;
  }

  if (!exist (xleng))
  {
    win_xleng_df = WIN_XLENG_DF;
  } else {
    win_xleng_df = xleng;
  }

  if (!exist (yleng))
  {
    win_yleng_df = WIN_YLENG_DF;
  } else {
    win_yleng_df = yleng;
  }

  if (!exist (xoff))
  {
    win_xoff_df = 0;
  } else {
    win_xoff_df = xoff;
  }

  if (!exist (yoff))
  {
    win_yoff_df = 0;
  } else {
    win_yoff_df = yoff;
  }
};

////////////
#
# Check to make sure a plot-object exists. If one
# does not exist, create it.
#

check_plot_object = function ()
{
  if (length (PGPLOT_PLOT_WINDOWS) == 0)
  {
    plstart();
    return 0;
  }

  return 1;
};

////////////
#
# Set the current plot window
# Default value = 0
#

plwin = function ( N )
{
  check_plot_object ();
  if (!exist (N)) { N = 1; }

  # Check to make sure N is valid

  for (i in members (PGPLOT_PLOT_WINDOWS))
  {
    if (N == strtod (i))
    {
      pgslct ( N );
      PGPLOT_ACTIVE_WIN = N;
      return PGPLOT_ACTIVE_WIN;
    }
  }

  printf ("plwin: invalid argument, N = %i\n", N);
  printf ("      valid values are:\n");

  PGPLOT_PLOT_WINDOWS?
};

#
# rlabplus extension
#
newplwin = function ()
{
  if (length (PGPLOT_PLOT_WINDOWS) == 0)
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
plwins = function( NWIN, dev, sz )
{
  //
  // process NWIN: an integer, number of windows requested,
  // or a matrix [nx(i);ny(i)]_i
  //
  if (!exist(NWIN))
  {
    // report what is available
    if (length(PGPLOT_PLOT_WINDOWS) == 0)
    {
      return <<win=[];act=0;dev=blank(0,0)>>;
    }

    retl = <<>>;
    retl.act = PGPLOT_ACTIVE_WIN;
    retl.win = [];
    retl.dev = blank(0,0);
    for(i in members(PGPLOT_PLOT_WINDOWS))
    {
      retl.win = [retl.win, strtod(i)];
      retl.dev    = [retl.dev, PGPLOT_PLOT_WINDOWS.[i].dev];
    }
    idx = sort(retl.win).idx;
    retl.dev = retl.dev[idx];
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
  }
  else if (NWIN.nc == 2)
  {
    nwin = NWIN.nr;
  }
  else
  {
    printf("plwins: improper first argument");
    return <<>>;
  }

  if (nwin<0 || nwin>8)
  {
    printf("plwins: improper first argument");
    return <<>>;
  }

  //
  // process devices
  //
  if (!exist(dev))
  { dev = "/XWIN"; }

  //
  // process sizes
  //
  if (!exist(sz))
  { sz = [win_xleng_df, win_yleng_df]; }
  if (exist(sz) && (sz.nr*sz.nc!=2 || class(sz) != "num"))
  { sz = [win_xleng_df, win_yleng_df]; }

  //
  // too many open windows
  //
  wavl = plwins().win;
  davl = plwins().dev;
  if (nwin <= length(wavl))
  {
    // kill extra windows
    for (i in nwin+1:length(wavl))
    {
      plwin(i);
      plclose();
    }

    // rearrange remaining windows so 1:nwin are open
    for (i in 1:nwin)
    {
      if (i == wavl[i])
      { continue; }
      if (!exist(PGPLOT_PLOT_WINDOWS.[i]))
      {
        plcopy( wavl[i], davl[i] ); // this will get copied to device 'i'
        plwin ( wavl[i] );          // now choose old window as default
        plclose ();                 // kill it
      }
    }
  }

  //
  // open/modify requested number of plot windows
  //
  wavl = plwins().win;
  davl = plwins().dev;
  for (i in 1:nwin)
  {
    ndev = dev[ min(i,length(dev)) ];
    nwid = sz[ min(i,sz.nr); 1 ] / win_wid;
    nasp = sz[ min(i,sz.nr); 2 ] / sz[ min(i,sz.nr); 1 ];

    nx = 1;
    ny = 1;
    if (NWIN.nc==2)
    {
      nx = NWIN[i;1];
      ny = NWIN[i;2];
    }

    if (!exist(PGPLOT_PLOT_WINDOWS.[i]))
    {
      plstart(nx,ny,ndev,nwid,nasp);
      pgpage ();
      continue;
    }
    else
    {
      // window exists and is open, check if the specifications have changed
      // and if so close it and open a new one
      plwin(i);
      if (PGPLOT_PLOT_WINDOWS.[i].nx != nx || PGPLOT_PLOT_WINDOWS.[i].ny != ny)
      {
        plclose();
        plstart(nx,ny,ndev,nwid,nasp);
        pgpage ();
      }
    } // if (!exist(PGPLOT_PLOT_WINDOWS.[i]))
  } // for (i in 1:nwin)

  return 0;
};

getplot = function ( win_no )
{
  local (win_no)

  if (length (PGPLOT_PLOT_WINDOWS) != 0)
  {
    if (!exist (win_no)) { win_no = P; }

    if (exist (PGPLOT_PLOT_WINDOWS.[win_no]))
    {
      return (PGPLOT_PLOT_WINDOWS.[win_no]);
    } else {
      return 0;
    }
  }
  return <<>>;
};

////////////
//
// Set/start/select the plot device
//
// plstart ( NX, NY, DEV, WIDTH, ASPECT )
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

plstart = function ( nx, ny, dev, width, aspect )
{
    if (!exist (nx))
    { nx = 1; }

    if (!exist (ny))
    { ny = 1; }

    if (!exist (dev))
    { dev = "/XWIN"; }

    // Create the plot-object
    // First, figure out the index
    if (!PGPLOT_ACTIVE_WIN)
    {
      PGPLOT_ACTIVE_WIN = 1;
    }
    else
    {
      //P = P + 1;
      PGPLOT_ACTIVE_WIN = newplwin();
    }

    create_plot_object (PGPLOT_ACTIVE_WIN, nx, ny);

    // Open the plot device.
    winid = pgopen (dev);

    // Setup the window width/aspect ratio.
    if (exist (width))
    {
      if (!exist (aspect))
      { aspect = 1; }

      pgpap (width, aspect);
    }


    // Turn off page prompting.
    pgask (0);

    // Set up the number of plots per window.
    pgsubp (nx, ny);

    // let the WIN variable knows what kind of window is open
    PGPLOT_PLOT_WINDOWS.[winid].dev = dev;

    // Update the screen.
    pgpage ();

    return winid;
};

////////////
//
// Close a plot device. We must destroy the current plot-object
// And switch the output stream back to the default.
//
plclose = function ()
{
  if (length(PGPLOT_PLOT_WINDOWS) > 1)
  {
    //
    // Clear PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN] and reset P to 1st plot-window
    //
    clear (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN]);
    pgclos ();
    PGPLOT_ACTIVE_WIN = strtod(members (PGPLOT_PLOT_WINDOWS)[1]);
    pgslct (PGPLOT_ACTIVE_WIN);
    return PGPLOT_ACTIVE_WIN;
  }
  else if (length(PGPLOT_PLOT_WINDOWS) == 1)
  {
    if (exist (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN]))
    {
      clear (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN]);
      PGPLOT_PLOT_WINDOWS = <<>>;
      PGPLOT_ACTIVE_WIN = 0;
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
  if (exist (PGPLOT_PLOT_WINDOWS))
  { clear (PGPLOT_PLOT_WINDOWS); }

  PGPLOT_ACTIVE_WIN = 0;
  PGPLOT_PLOT_WINDOWS = <<>>;
};

////////////
//
// Change plot aspect ratio
//

plaspect = function ( aspect )
{
  check_plot_object ();
  i = mod (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].subplot, PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].nplot) + 1;
  if (!exist (aspect))
  {
    PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].aspect[i] = 0;
  }
  else
  {
    PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].aspect[i] = aspect;
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
  i = mod (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].subplot, PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].nplot) + 1;

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
      PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].style.[i] = style;
    }
    return 1;
  }
  PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].style.[i] = "line";
  return 1;
};

////////////
//
// Control of the plot line style.
// There are 5 line styles

plline = function ( line_style )
{
  check_plot_object ();
  i = mod (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].subplot, PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].nplot) + 1;

  if (exist (line_style))
  {
    if (class (line_style) == "num")
    {
      if (line_style.n != 20)
      {
	error ("plline: LVEC must be 1x20 in size");
      }
      PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].lstyle[i;] = line_style;
    }
    return 1;
  }
  PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].lstyle[i;] = [1:5,1:5,1:5,1:5];
  return 1;
};

////////////
//
// Control of the plot line color
//

plcolor = function ( COLOR )
{
  check_plot_object ();
  i = mod (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].subplot, PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].nplot) + 1;

  if (exist (COLOR))
  {
    if (class (COLOR) == "num")
    {
      PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].color[i;] = COLOR;
    }
    return 1;
  }
  PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].COLOR[i;] = 1:14;
  return 1;
};

////////////
//
// Control of the plot point style.
// There are 8 line styles

plpoint = function ( point_style )
{
  check_plot_object ();
  i = mod (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].subplot, PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].nplot) + 1;

  if (exist (point_style))
  {
    if (class (point_style) == "num")
    {
      if (point_style.nc != 8)
      {
	error ("plpoint: PVEC must be 1x8 in size");
      }
      PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].pstyle[i;] = point_style;
    }
    return 1;
  }
  PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].pstyle[i;] = 1:8;
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
  i = mod (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].subplot, PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].nplot) + 1;

  if (!exist (font)) { font = 1; }

  if (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].fontld == 0)
  {
    pgscf (1);
    PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].fontld = 1;
  }

  PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].font[i] = font;
  return PGPLOT_ACTIVE_WIN;
};

////////////
//
// Change character height
//

plch = function ( H )
{
  check_plot_object ();
  i = mod (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].subplot, PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].nplot) + 1;

  if (!exist (H)) { H = 1; }

  PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].char_height = H;

  return PGPLOT_ACTIVE_WIN;
};

////////////
//
// Change pen width
//

plwid = function ( width )
{
  check_plot_object ();
  i = mod (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].subplot, PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].nplot) + 1;

  if (!exist (width))
  {
    PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].width = ones(1,32);
  } else {
    PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].width[i;] = width;
  }
  return PGPLOT_ACTIVE_WIN;
};

plfillstyle = function ( FS )
{
  check_plot_object ();
  i = mod (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].subplot, PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].nplot) + 1;

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

  PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].fillstyle[i] = FS;

  return PGPLOT_ACTIVE_WIN;
};

////////////
//
// Place some text on the plot
//

plptex = function ( text, x , y , angle , just )
{
  p = mod (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].subplot, PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].nplot) + 1;	// The current index

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
  pgsci (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].color[p;1]);

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

//////////////////////////////////////////////////////////////////////////////
//
// Plot the columns of a matrix (X-Y plot).
//
//////////////////////////////////////////////////////////////////////////////

plot = function ( data )
{
  check_plot_object ();

  p = mod (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].subplot, PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].nplot) + 1;	// The current index
  key = 1;

  //
  // Draw the graph
  // Step through the matrix plotting
  // each column versus the 1st
  //

  if (subplot_spec)
  {
    subplot_horizontal
    pgpanl (subplot_horizontal, subplot_vertical);
    subplot_spec = 0;
  } else {
    pgpage ();
  }

  if (class (data) == "num" || ishist(data))
  {
    // plot a single matrix or a special list: histogram
    if (ishist(data))
    { data = hist_line(data); }

    //
    // Set up the plot basics
    //
    if (abs (key) > data.nc)
    { error ("plot: KEY argument > M.nc"); }

    //
    // Add row-values if necessary.
    //
    if (data.nc == 1)
    { data = [(1:data.nr)', data]; }

    //
    // Compute scale limits
    //
    </xmax;xmin/> = x_scales (real(data)[;1], p);
    </ymax;ymin/> = y_scales (real(data)[;2:data.nc], p);

    //
    // Now, pick the tickmark properties, maybe reset
    // scale min/max values.
    //
    </ xmax; xmin; nxsub; xtick /> = ticks (xmin, xmax);
    </ ymax; ymin; nysub; ytick /> = ticks (ymin, ymax);

    //
    // Set up the window/viewport.
    //
    pgsch (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].char_height);       // Set character height

    pgvstd ();

    if (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].aspect[p] != 0)
    {
      pgwnad (xmin, xmax, ymin, ymax);
    } else {
      pgswin (xmin, xmax, ymin, ymax);
    }

    pgsls (1);                         // Set line style.
    pgsci (14);                        // Set color light for background grids.
    pgbox ("g", xtick, nxsub,...
           "g", ytick, nysub);

    pgsci (10);                        // Set back to normal grid, label color
    pgbox (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].gridx[p], xtick, nxsub, ...
           PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].gridy[p], ytick, nysub);

    if (plot_matrix ( data, key, p, 0, xmin, xmax, ymin, ymax, ymax-ymin ) < 0)
    { return -1; }

  } else { if (class (data) == "list")
  {
    </xmax;xmin;ymax;ymin/> = list_scales ( data, key, p );

    //
    // Now, pick the tickmark properties, maybe reset
    // scale min/max values.
    //
    </ xmax; xmin; nxsub; xtick /> = ticks (xmin, xmax);
    </ ymax; ymin; nysub; ytick /> = ticks (ymin, ymax);

    //
    // Set up the window/viewport.
    //
    pgsch (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].char_height);       // Set character height
    pgvstd ();

    if (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].aspect[p] != 0)
    {
      pgwnad (xmin, xmax, ymin, ymax);
    } else {
      pgswin (xmin, xmax, ymin, ymax);
    }

    pgsls (1);                         // Set line style.
    pgsci (14);
    pgbox ("g", xtick, nxsub,...
           "g", ytick, nysub);

    pgsci (10);
    pgbox (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].gridx[p], xtick, nxsub,...
           PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].gridy[p], ytick, nysub);

    if (plot_list ( data, key, p, xmin, xmax, ymin, ymax ) < 0)
    {	return -1; }

  } else {
    error ("plot: un-acceptable argument");
  }}

  pgsci (10);                        // Set color back.
  pgsls (1);                         // Set line style.
  pgsch (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].char_height);       // Set character height

  // pglab(PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].xlabel[p], PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].ylabel[p], PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].title[p]);
  pglab(PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].xlabel[p], " ", PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].title[p]);
  pgmtxt ("L", 3, 0.5, 0.5, PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].ylabel[p]);

  pgupdt ();

  //
  // Increment the plot no. so that next time
  // we use the correct settings.
  //
  PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].subplot = PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].subplot + 1;
  return PGPLOT_ACTIVE_WIN;
};

//
// plot a data set to an eps file using the formatting information
// of the current plot window: modification of plot() that accounts
// for change in colors (green-on-black on terminal to black-on-white
// on paper) linewidths for labels and legends.
// Marijan Kostrun, IX-2006
//
epsplot = function ( data, fname )
{

  check_plot_object ();
  p = mod (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].subplot, PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].nplot) + 1; // The current index

  if (!exist(data))
  { stop ("epsplot: No data to plot!"); }

  if (!exist(fname))
  { stop ("epsplot: Missing filename to plot to!"); }

  plcopy (P, fname + "/CPS" );

  key = 1;

  //
  // Draw the graph
  // Step through the matrix plotting
  // each column versus the 1st
  //

  if (subplot_spec)
  {
    subplot_horizontal
    pgpanl (subplot_horizontal, subplot_vertical);
    subplot_spec = 0;

  } else {
    pgpage ();
  }

  if (class (data) == "num")
  {
    //
    // Set up the plot basics
    //
    if (abs (key) > data.nc)
    { error ("plot: KEY argument > M.nc"); }

    //
    // Add row-values if necessary.
    //
    if (data.nc == 1)
    { data = [(1:data.nr)', data]; }

    //
    // Compute scale limits
    //
    </xmax;xmin/> = x_scales (real(data)[;1], p);
    </ymax;ymin/> = y_scales (real(data)[;2:data.nc], p);

    //
        // Now, pick the tickmark properties, maybe reset
        // scale min/max values.
    //
    </ xmax; xmin; nxsub; xtick /> = ticks (xmin, xmax);
    </ ymax; ymin; nysub; ytick /> = ticks (ymin, ymax);

    //
    // Set up the window/viewport.
    //
    pgsch (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].char_height);       // Set character height

    pgvstd ();

    if (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].aspect[p] != 0)
    {
      pgwnad (xmin, xmax, ymin, ymax);

      } else {
        pgswin (xmin, xmax, ymin, ymax);
    }

    // black frame with gray grid
    // solid gray line for the tick lines
    pgsls (1);        // solid
    pgsci (15);       // light gray
    pgslw (3);        // 3pt thick
    pgbox ("g", xtick, nxsub,...
        "g", ytick, nysub);
    pgsci (1);        // black
    pgslw (5);        // 5pt thick
    pgbox (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].gridx[p], xtick, nxsub,...
        PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].gridy[p], ytick, nysub);

    if (plot_matrix ( data, key, p, 0, xmin, xmax, ymin, ymax, ymax-ymin ) < 0)
    { return -1; }

  } else { if (class (data) == "list")
  {

    </xmax;xmin;ymax;ymin/> = list_scales ( data, key, p );

    //
    // Now, pick the tickmark properties, maybe reset
    // scale min/max values.
    //
    </ xmax; xmin; nxsub; xtick /> = ticks (xmin, xmax);
    </ ymax; ymin; nysub; ytick /> = ticks (ymin, ymax);

    //
    // Set up the window/viewport.
    //
    pgsch (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].char_height);       // Set character height

    pgvstd ();

    if (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].aspect[p] != 0)
    {
      pgwnad (xmin, xmax, ymin, ymax);

    } else {
      pgswin (xmin, xmax, ymin, ymax);
    }

    // black frame with gray grid
    // solid gray line for the tick lines
    // solid gray line for the tick lines
    pgsls (1);        // solid
    pgsci (15);       // light gray
    pgslw (3);        // 3pt thick
    pgbox ("g", xtick, nxsub,...
        "g", ytick, nysub);
    pgsci (1);        // black
    pgslw (5);        // 5pt thick
    pgbox (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].gridx[p], xtick, nxsub,...
        PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].gridy[p], ytick, nysub);

    if (plot_list ( data, key, p, xmin, xmax, ymin, ymax ) < 0)
    { return -1; }

  } else {
    error ("plot: un-acceptable argument");
  }}

  pgsci (1);                      // Set color back.
  pgsls (1);                      // Set line style.
  pgslw (5);                      // 5pt thick
  pgsch (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].char_height);    // Set character height

  //pglab(PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].xlabel[p], PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].ylabel[p], PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].title[p]);
  pglab(PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].xlabel[p], " ", PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].title[p]);
  pgmtxt ("L", 3, 0.5, 0.5, PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].ylabel[p]);

  pgupdt ();

  //
  // Increment the plot no. so that next time
  // we use the correct settings.
  //
  PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].subplot = PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].subplot + 1;

  plclose ();

  return PGPLOT_ACTIVE_WIN;
};

////////////
//
// Plot the contours of a 2D data set.
//
// The data is composed in a list, with
// elements `x', `y', and `z'. x and y are single-dimension arrays
// (row or column matrices), and z is a two-dimensional array. The
// array z, is a function of x and y: z = f(x,y). Thus, the values in
// the array x can be thought of a "row-labels", and the values of y
// can be thought of as "column-lables" for the 2-dimensioal array z.
//
////////////

plcont = function ( CL, intval, minint )
{
  check_plot_object ();

  //
  // 1st check list contents
  //

  if (exist (CL)) { check_3d_list (CL); }
  if (!exist (intval)) { intval = 20; }
  if (!exist (minint)) { minint = 20; }

  p = mod (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].subplot, PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].nplot) + 1;	// The current index

  //
  // Figure out the scale limits.
  // Needs improvement!
  //

  xmin = xmax = ymin = ymax = zmin = zmax = 0;
  if (exist (CL))
  {
    </Xmax;Xmin;Ymax;Ymin;Zmax;Zmin/> = XYZ_scales (CL.x, CL.y, CL.z, p);
    if (Xmin < xmin) { xmin = Xmin; } if (Xmax > xmax) { xmax = Xmax; }
    if (Ymin < ymin) { ymin = Ymin; } if (Ymax > ymax) { ymax = Ymax; }
    if (Zmin < zmin) { zmin = Zmin; } if (Zmax > zmax) { zmax = Zmax; }
  }

  if (subplot_spec)
  {
    subplot_horizontal
    pgpanl (subplot_horizontal, subplot_vertical);
    subplot_spec = 0;
  } else {
    pgpage ();
  }

  //
  // Set up the 1st viewport for drawing the plot.
  //

  pgvstd ();

  if (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].aspect[p] != 0)
  {
    pgwnad (1, CL.x.n, 1, CL.y.n);
  } else {
    pgswin (1, CL.x.n, 1, CL.y.n);
  }

  //
  // Now, pick the tickmark properties, maybe reset
  // scale min/max values.
  //

  </ xmax; xmin; nxsub; xtick /> = ticks (xmin, xmax);
  </ ymax; ymin; nysub; ytick /> = ticks (ymin, ymax);

  pgsci (10);                        // Set back to normal grid, label color
  pgsch (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].char_height);       // Set character height
  pgbox (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].gridx[p], xtick, nxsub, ...
         PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].gridy[p], ytick, nysub);

  // Convert the data to log data if necessary.
  if (find_char (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].gridx[p], "l"))
  { x = log10 (real (CL.x)); } else { x = real (CL.x); }
  if (find_char (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].gridy[p], "l"))
  { y = log10 (real (CL.y)); } else { y = real (CL.y); }
  z = real (CL.z);

  if (exist (CL.clevel))
  {
    clevel = CL.clevel;
  } else {
    clevel = linspace(zmin, zmax, 10);
  }

  //
  // Draw the contours
  //

  pgsci (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].color[p;2]);          // Set line color.
  pgslw (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].width[p]);            // Set line width.

  pgcont (CL.z, CL.x.n, CL.y.n, 1, CL.x.n, 1, CL.y.n, ...
          clevel, clevel.n, [0,1,0,0,0,1]);

  for (i in 1:clevel.n)
  {
    cont_label = num2str (clevel[i]);
    pgconl (CL.z, CL.x.n, CL.y.n, 1, CL.x.n, 1, CL.y.n, ...
            clevel[i], [0,1,0,0,0,1], cont_label, intval, minint );
  }

  //
  // Reset color and draw the labels.
  //

  pgsci (10);
  pglab (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].xlabel[p], PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].ylabel[p], PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].title[p]);

  //
  // Increment the plot no. so that next time
  // we use the correct settings.
  //

  PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].subplot = PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].subplot + 1;

  return PGPLOT_ACTIVE_WIN;
};

////////////
// Plot a 2D image.
////////////

plimag = function ( CL )
{
  check_plot_object ();

  //
  // 1st check list contents
  //

  if (exist (CL)) { check_3d_list (CL); }

  p = mod (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].subplot, PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].nplot) + 1;	// The current index

  //
  // Figure out the scale limits.
  // Needs improvement!
  //

  xmin = xmax = ymin = ymax = zmin = zmax = 0;
  if (exist (CL))
  {
    </Xmax;Xmin;Ymax;Ymin;Zmax;Zmin/> = XYZ_scales (CL.x, CL.y, CL.z, p);
    if (Xmin < xmin) { xmin = Xmin; } if (Xmax > xmax) { xmax = Xmax; }
    if (Ymin < ymin) { ymin = Ymin; } if (Ymax > ymax) { ymax = Ymax; }
    if (Zmin < zmin) { zmin = Zmin; } if (Zmax > zmax) { zmax = Zmax; }
  }

  if (subplot_spec)
  {
    subplot_horizontal
    pgpanl (subplot_horizontal, subplot_vertical);
    subplot_spec = 0;
  } else {
    pgpage ();
  }

  //
  // Set up the 1st viewport for drawing the plot.
  //

  pgvstd ();

  if (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].aspect[p] != 0)
  {
    pgwnad (1, CL.x.n, 1, CL.y.n);
    } else {
    pgswin (1, CL.x.n, 1, CL.y.n);
  }

  //
  // Now, pick the tickmark properties, maybe reset
  // scale min/max values.
  //

  </ xmax; xmin; nxsub; xtick /> = ticks (xmin, xmax);
  </ ymax; ymin; nysub; ytick /> = ticks (ymin, ymax);

  pgsci (10);                        // Set back to normal grid, label color
  pgsch (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].char_height);       // Set character height
  pgbox (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].gridx[p], xtick, nxsub, ...
         PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].gridy[p], ytick, nysub);

  // Convert the data to log data if necessary.
  if (find_char (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].gridx[p], "l"))
  { x = log10 (real (CL.x)); } else { x = real (CL.x); }
  if (find_char (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].gridy[p], "l"))
  { y = log10 (real (CL.y)); } else { y = real (CL.y); }
  z = real (CL.z);

  //
  // Draw the image
  //

  pgimag (CL.z, CL.x.n, CL.y.n, 1, CL.x.n, 1, CL.y.n, ...
          zmin, zmax, [0,1,0,0,0,1]);

  //
  // Reset color and draw the labels.
  //

  pgsci (10);                        // Set color back.
  pgsls (1);                         // Set line style.
  pgsch (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].char_height);       // Set character height

  pglab (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].xlabel[p], PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].ylabel[p], PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].title[p]);

  //
  // Increment the plt no. so that next time
  // we use the correct settings.
  //

  PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].subplot = PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].subplot + 1;

  return PGPLOT_ACTIVE_WIN;
};

////////////
//
// Plot histogram (bar-chart) of binned data. That is, the data has already
// been organized into bins by the user, and is supplied in the form of:
//
//  [ X-axis values (center), Bin-height-1, Bin-height-2, ... ]
//
////////////

plbar = function ( data )
{
  check_plot_object ();

  p = mod (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].subplot, PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].nplot) + 1;	// The current index

  //
  // Draw the graph
  // Step through the matrix plotting
  // each column versus the 1st
  //

  if (subplot_spec)
  {
    subplot_horizontal
    pgpanl (subplot_horizontal, subplot_vertical);
    subplot_spec = 0;
  } else {
    pgpage ();
  }

  if (data.nc == 1)
  {
    //
    // Fix up the array so we only have to handle one case.
    //

    data = [ (1:data.nr)', data ];
  }

  //
  // Compute scale limits
  //

  // Check for log scales...
  if (find_char (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].gridx[p], "l"))
  {
    data[;1] = log10 (data[;1]);
  }

  // Compute our own xmin/xmax...
  delx = min (diff (data[;1]))/2;
  delx = delx * 0.8;
  xmin = min (data[;1]) - delx;
  xmax = max (data[;1]) + delx;

  if (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].xmin[p] != 1j) { xmin = real (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].xmin[p]); }
  if (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].xmax[p] != 1j) { xmax = real (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].xmax[p]); }

  // Use standard method for y-scales.
  </ymax; ymin/> = y_scales (real (data[;2:data.nc]), p);

  // Adjust Y for log-scale if necessary...
  if (find_char (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].gridy[p], "l"))
  {
    data[;2:data.nc] = log10 (data[;2:data.nc]);
  }

  //
  // Set up the window/viewport.
  //

  pgvstd ();

  if (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].aspect[p] != 0)
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
  pgsch (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].char_height);       // Set character height
  pgbox (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].gridx[p], xtick, nxsub, ...
         PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].gridy[p], ytick, nysub);

  pgsch (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].char_height);      // Set character height

  //
  // Plot/draw the rectangles/bins...
  //

  pgsfs (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].fillstyle[p]);     // Set the fillstyle.

  for (i in 1:data.nr)          // Loop over each row (X-axis).
  {
    x1 = data[i;1] - delx;
    x2 = data[i;1] + delx;

    for (j in 2:data.nc)        // Loop over each column (N-bins per X-value).
    {
      // Scootch each bin (after the first) over a little
      // so that we can see each bin.
      x1 = data[i;1] - (delx / (j/2));

      // Compute the Y values for each bin / column
      y1 = 0;
      y2 = data[i;j];

      pgsci (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].color[p;j]);     // Set line color.
      pgrect (x1, x2, y1, y2);
    }
  }

  pgsls (1);                         // Set line style.
  pgsci (10);                        // Set color.
  pgsch (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].char_height);       // Set character height.

  pglab(PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].xlabel[p], PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].ylabel[p], PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].title[p]);
  pgupdt ();

  //
  // Increment the plot no. so that next time
  // we use the correct settings.
  //

  PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].subplot = PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].subplot + 1;
  return PGPLOT_ACTIVE_WIN;
};

////////////
//
// Plot histogram of binned data.
//
////////////

plhist = function ( data, nbin )
{
  if (!exist (data) || isempty(data))
  {
    //error("plhist: Trying to plot emtpy or non-existing data set!");
    stop("plhist: Trying to plot emtpy or non-existing data set!");
  }
  if (!exist (nbin)) { nbin = 10; }

  check_plot_object ();

  p = mod (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].subplot, PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].nplot) + 1;	// The current index

  //
  // Draw the graph
  // Step through the matrix plotting
  // each column versus the 1st
  //

  if (subplot_spec)
  {
    subplot_horizontal
    pgpanl (subplot_horizontal, subplot_vertical);
    subplot_spec = 0;
  } else {
    pgpage ();
  }

  pgsch (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].char_height);      // Set character height
  pgsci (1);                        // Set color
  pgsls (1);                        // Set line style
  pgslw (1);
  pgscf (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].font[p]);          // Set font.

  //
  // Compute scale limits
  //

  xmin = min( min (data) );
  xmax = max( max (data) );
  ymin = 0;
  ymax = hist_scales (data, nbin);

  //
  // Check computed scale limits against user's
  //

  if (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].ymin[p] != 1j) { ymin = PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].ymin[p]; }
  if (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].ymax[p] != 1j) { ymax = PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].ymax[p]; }

  //
  // Set up the window/viewport.
  //

  pgvstd ();

  if (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].aspect[p] != 0)
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
  pgsch (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].char_height);       // Set character height
  pgbox (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].gridx[p], xtick, nxsub, ...
         PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].gridy[p], ytick, nysub);

  pgsci (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].color[p;2]);        // Set line color.

  pghist (data.n, data, min (data), max (data), nbin, 1);

  pgsls (1);
  pgsci (10);

  pglab(PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].xlabel[p], PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].ylabel[p], PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].title[p]);
  pgupdt ();

  //
  // Increment the plot no. so that next time
  // we use the correct settings.
  //

  PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].subplot = PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].subplot + 1;
  return PGPLOT_ACTIVE_WIN;
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
  if (!exist (ID))
  { ID = 1; }

  if (!exist (DEV))
  { error ("plcopy: must specify a plot device"); }

  if (!exist (PGPLOT_PLOT_WINDOWS.[ID]))
  { error ("plcopy: must specify valid plot-object"); }

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
      return 0;   // unsuccessful
    }
  }

  //
  // Get the new object ID
  //
  if (!exist (P))
  {
    PGPLOT_ACTIVE_WIN = 1;
  }
  else
  {
    //
    // figure out new P
    //
    //P = P + 1;
    PGPLOT_ACTIVE_WIN = newplwin();
  }

  //
  // Perform the copy.
  //
  PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN] = PGPLOT_PLOT_WINDOWS.[ID];

  // Open the plot device.
  winid = pgopen (DEV);

  // Note the device.  kmk 2005.
  PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].device = DEV;

  // Turn off page prompting.
  pgask (0);

  // Set up the number of plots per window.
  pgsubp (PGPLOT_PLOT_WINDOWS.[ID].nx, PGPLOT_PLOT_WINDOWS.[ID].ny);

  return 0;   // successful
};

////////////
//
// Set the X-axis label
//

xlabel = function ( xstr )
{
  check_plot_object ();
  if (!exist (xstr)) { xstr = ""; }
  i = mod (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].subplot, PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].nplot) + 1;
  PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].xlabel[i] = xstr;
};

////////////
//
// Set the Y-axis label
//

ylabel = function ( xstr )
{
  check_plot_object ();
  if (!exist (xstr)) { xstr = ""; }
  i = mod (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].subplot, PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].nplot) + 1;
  PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].ylabel[i] = xstr;
};

////////////
//
// Set the Z-axis label
//

zlabel = function ( xstr )
{
  check_plot_object ();
  if (!exist (xstr)) { xstr = ""; }
  i = mod (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].subplot, PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].nplot) + 1;
  PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].zlabel[i] = xstr;
};

////////////
//
// Set the plot-title
//

pltitle = function ( xstr )
{
  check_plot_object ();

  if (!exist (xstr))
  { xstr = ""; }

  p = mod (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].subplot, ...
      PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].nplot) + 1;	// The current index

  PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].title[p] = xstr;

  return 0;
};

////////////
//
// Set the scale limits.
//

plimits = function ( xmin, xmax, ymin, ymax, zmin, zmax )
{
  check_plot_object ();
  i = mod (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].subplot, PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].nplot) + 1;

  if (exist (xmin))
  {
    PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].xmin[i] = xmin;
  }
  else
  {
    PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].xmin[i] = 1j;
  }

  if (exist (xmax))
  {
    PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].xmax[i] = xmax;
  }
  else
  {
    PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].xmax[i] = 1j;
  }

  if (exist (ymin))
  {
    PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].ymin[i] = ymin;
  }
  else
  {
    PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].ymin[i] = 1j;
  }

  if (exist (ymax))
  {
    PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].ymax[i] = ymax;
  }
  else
  {
    PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].ymax[i] = 1j;
  }

  if (exist (zmin))
  {
    PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].zmin[i] = zmin;
  }
  else
  {
    PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].zmin[i] = 1j;
  }

  if (exist (zmax))
  {
    PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].zmax[i] = zmax;
  }
  else
  {
    PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].zmax[i] = 1j;
  }

  return 0;
};

////////////
//
// Set 2-D grid styles. A not-so-friendly interface.
//

plgrid = function ( sty_x, sty_y )
{
  check_plot_object ();
  i = mod (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].subplot, PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].nplot) + 1;

  if (exist (sty_x))
  {
    if (class (sty_x) == "string")
    {
      PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].gridx[i] = sty_x;
    } else {
      error ("plgrid: requires string argument GRID_STY_X");
    }
  } else {
    PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].gridx[i] = grid_x_default;
  }
  if (exist (sty_y))
  {
    if (class (sty_y) == "string")
    {
      PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].gridy[i] = sty_y;
    } else {
      error ("plgrid: requires string argument GRID_STY_Y");
    }
  } else {
    PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].gridy[i] = grid_y_default;
  }
};

////////////
//
// Set 3-D grid (axis) styles
//

plgrid3 = function ( sty_x, sty_y, sty_z )
{
  check_plot_object ();
  i = mod (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].subplot, PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].nplot) + 1;
  if (exist (sty_x))
  {
    if (class (sty_x) == "string")
    {
      PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].grid3x[i] = sty_x;
    } else {
      error ("plgrid3: requires string argument GRID_STY_X");
    }
  } else {
    PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].grid3x[i] = grid_3x_default;
  }
  if (exist (sty_y))
  {
    if (class (sty_y) == "string")
    {
      PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].grid3y[i] = sty_y;
    } else {
      error ("plgrid3: requires string argument GRID_STY_Y");
    }
  } else {
    PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].grid3y[i] = grid_3y_default;
  }
  if (exist (sty_z))
  {
    if (class (sty_z) == "string")
    {
      PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].grid3z[i] = sty_z;
    } else {
      error ("plgrid3: requires string argument GRID_STY_Z");
    }
  } else {
    PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].grid3z[i] = grid_3z_default;
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
  i = mod (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].subplot, PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].nplot) + 1;

  if (exist (X_STR))
  {
    if (X_STR == "log") { PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].gridx[i] = "bcnstl"; }
  } else {
    PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].gridx[i] = grid_x_default;
  }

  if (exist (Y_STR))
  {
    if (Y_STR == "log") { PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].gridy[i] = "bcnstlv"; }
  } else {
    PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].gridy[i] = grid_y_default;
  }
  return PGPLOT_ACTIVE_WIN;
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

  if (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].xmin[p] != 1j) { xmin = real (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].xmin[p]); }
  if (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].xmax[p] != 1j) { xmax = real (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].xmax[p]); }

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

  if (find_char (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].gridx[p], "l"))
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

  if (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].ymin[p] != 1j) { xmin = real (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].ymin[p]); }
  if (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].ymax[p] != 1j) { xmax = real (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].ymax[p]); }

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

  if (find_char (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].gridy[p], "l"))
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

  if (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].zmin[p] != 1j) { xmin = real (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].zmin[p]); }
  if (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].zmax[p] != 1j) { xmax = real (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].zmax[p]); }

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

  if (find_char (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].gridz[p], "l"))
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

  if (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].xmin[p] != 1j) { xmin = real (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].xmin[p]); }
  if (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].xmax[p] != 1j) { xmax = real (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].xmax[p]); }
  if (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].ymin[p] != 1j) { ymin = real (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].ymin[p]); }
  if (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].ymax[p] != 1j) { ymax = real (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].ymax[p]); }
  if (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].zmin[p] != 1j) { zmin = real (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].zmin[p]); }
  if (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].zmax[p] != 1j) { zmax = real (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].zmax[p]); }

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

  if (find_char (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].grid3x[p], "l"))
  {
    if (xmin <= 0 || xmax <= 0) { error ("plot: cannot plot log(x<=0)"); }
    xmin = log10 (xmin);
    xmax = log10 (xmax);
  }

  if (find_char (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].grid3y[p], "l"))
  {
    if (ymin <= 0 || ymax <= 0) { error ("plot: cannot plot log(y<=0)"); }
    ymin = log10 (ymin);
    ymax = log10 (ymax);
  }

  if (find_char (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].grid3z[p], "l"))
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

    if (find_char (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].gridx[p], "l"))
    { x = log10 (x); }
    if (find_char (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].gridy[p], "l"))
    { y = log10 (y); }

    pgsci (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].color[p;k]);     // Set line color.
    pgsls (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].lstyle[p;l]);    // Set line style.
    pgslw (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].width[p;l]);     // Set line width.

    if (get_style (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].style.[p], k-1) == "line")
    {
      pgline (M.nr, x, y);
    } else { if (get_style (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].style.[p], k-1) == "point") {
      pgpt (M.nr, x, y, PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].pstyle[p;l]);
    } else { if (get_style (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].style.[p], k-1) == "line-point") {
      pgline (M.nr, x, y);
      pgpt (M.nr, x, y, PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].pstyle[p;l]);
    } else { {
      pgline (M.nr, x, y);
    }}}}

    //
    // Now do the legend
    //
    if (!any (any (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].desc.[p] == 1j)))
    {
      // Use the default if necessary
      if (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].desc.[p][1] == "default")
      {
	       desc = "c1";
      } else { if (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].desc.[p].n >= k-1)
      { desc = PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].desc.[p][k-1];
      } else {
	       // Not sure what to do, user has messed up.
	       desc = "";
      }}

      v = v - (ymax-ymin)/11;
      xl = (xmax-xmin)*[10.5/12, 11/12, 11.5/12]' + xmin;
      yl = [v, v, v]' + ymin;

      if (get_style (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].style.[p], k-1) == "line")
      {
	       pgline (3, xl, yl);
      } else { if (get_style (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].style.[p], k-1) == "point")
      { pgpt (3, xl, yl, PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].pstyle[p;l]);
      } else { if (get_style (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].style.[p], k-1) == "line-point")
      {
        pgline (3, xl, yl);
	      pgpt (3, xl, yl, PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].pstyle[p;l]);
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
      if (find_char (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].gridx[p], "l"))
      { x = log10 (x); }
      if (find_char (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].gridy[p], "l"))
      { y = log10 (y); }

      k = mod (i-1 + K, 14) + 1;
      l = mod (8 + i-2 + K, 8) + 1;

      pgsci (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].color[p;k]);     // Set line color.
      pgsls (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].lstyle[p;l]);    // Set line style.
      pgslw (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].width[p;l]);     // Set line width.

      if (get_style (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].style.[p], k-1) == "line")
      {
        pgline (np, x, y);
      } else { if (get_style (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].style.[p], k-1) == "point") {
        pgpt (np, x, y, PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].pstyle[p;l]);
      } else { if (get_style (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].style.[p], k-1) == "line-point") {
        pgline (np, x, y);
        pgpt (np, x, y, PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].pstyle[p;l]);
      } else { {
        pgline (np, x, y);
      }}}}

      //
      // Now do the legend
      //
      if (!any (any (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].desc.[p] == 1j)))
      {
	       // Use the default if necessary
	    if (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].desc.[p][1] == "default")
	    {
	      desc = "c" + num2str (i);
      } else { if (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].desc.[p].n >= k-1) {
        if (k > 1)
        {
          desc = PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].desc.[p][k-1];
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

      if (get_style (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].style.[p], k-1) == "line")
      {
        pgline (3, xl, yl);
      } else { if (get_style (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].style.[p], k-1) == "point") {
        pgpt (3, xl, yl, PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].pstyle[p;l]);
      } else { if (get_style (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].style.[p], k-1) == "line-point") {
        pgline (3, xl, yl);
        pgpt (3, xl, yl, PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].pstyle[p;l]);
      }}}

      pgslw  (5);                                         // set 5pt for legend text label
      pgptxt (xl[1]-(xmax-xmin)/25, yl[3], 0, 1, desc);
      pgslw (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].width[p;l]);                         // Set previous line width.

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
  } else { if (class (LIST.x) != "num") {
    error ("plot3: x must be numeric");
  } }
  if (!exist (LIST.y)) {
    error ("plot3: arg must contain `y' member");
  } else { if (class (LIST.y) != "num") {
    error ("plot3: y must be numeric");
  } }
  if (!exist (LIST.z)) {
    error ("plot3: arg must contain `z' member");
  } else { if (class (LIST.z) != "num") {
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

////////////
//
// Create a legend in the current plot window
//
// if pobj.desc.[p] = inf()		no legend
// if pobj.desc.[p] = "default"		default ("c1", "c2", ...)
// if pobj.desc.[p] = "string"		use "string" as description
//

//
// Set the current plot legend string
//

plegend = function ( LEGEND )
{
  check_plot_object ();

  p = mod (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].subplot, PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].nplot) + 1;	// The current index

  if (!exist (LEGEND))
  {
    PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].desc.[p] = 1j;
    return PGPLOT_ACTIVE_WIN;
  }

  if (class (LEGEND) == "string")
  {
    PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].desc.[p] = LEGEND;
  }

  return PGPLOT_ACTIVE_WIN;
};

set3d = function (bx, by, h)
{
  if (!exist (bx)) { basex = 2; } else { basex = bx; }
  if (!exist (by)) { basey = 2; } else { basey = by; }
  if (!exist (h)) { height = 4; } else { height = h; }
};

//---------------------------------------------------------------------
//
// Plot a surface...
//

plsurf = function (x, y, z, theta, phi)
{
  global (pi)

  check_plot_object ();
  p = mod (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].subplot, PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].nplot) + 1;	// The current index

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
      tmp = rotate3 ( [xm[i;j], ym[i;j], zm[i;j]], theta*dtor );

      // Now do the second rotation...
      tmp = rotate1 ( tmp, phi*dtor );

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

  if (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].xmin[p] != 1j) { xmin = real (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].xmin[p]); }
  if (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].xmax[p] != 1j) { xmax = real (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].xmax[p]); }

  if (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].ymin[p] != 1j) { ymin = real (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].ymin[p]); }
  if (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].ymax[p] != 1j) { ymax = real (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].ymax[p]); }

  //
  // Set up the window/viewport.
  //

  pgpage ();
  pgvstd ();

  if (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].aspect[p] != 0)
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
  pgslw (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].width[p]);         // Line width
  pgscf (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].font[p]);          // Set character font
  pgsch (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].char_height);      // Set character height

  pgbox (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].gridx[p], xtick, nxsub, ...
         PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].gridy[p], ytick, nysub);

  pgsci (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].color[p;1]);       // Set line color.

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

  pglab(PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].xlabel[p], PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].ylabel[p], PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].title[p]);
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
  p = mod (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].subplot, PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].nplot) + 1; // The current index

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
      tmp = rotate3 ( [xm[i;j], ym[i;j], zm[i;j]], theta*dtor );

      // Now do the second rotation...
      tmp = rotate1 ( tmp, phi*dtor );

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

  if (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].xmin[p] != 1j) { xmin = real (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].xmin[p]); }
  if (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].xmax[p] != 1j) { xmax = real (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].xmax[p]); }

  if (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].ymin[p] != 1j) { ymin = real (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].ymin[p]); }
  if (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].ymax[p] != 1j) { ymax = real (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].ymax[p]); }

  if (axis)
  {
    // Rotate the axes...
    tmp = rotate3 ( xmax * axis3d[1;], theta*dtor );
    ax1 = rotate1 ( tmp, phi*dtor );

    tmp = rotate3 ( ymax * axis3d[2;], theta*dtor );
    ax2 = rotate1 ( tmp, phi*dtor );

    tmp = rotate3 ( ymax * axis3d[3;], theta*dtor );
    ax3 = rotate1 ( tmp, phi*dtor );
  }

  //
  // Set up the window/viewport.
  //

  pgsls (1);
  pgsci (1);

  pgpage ();
  pgvstd ();

  if (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].aspect[p] != 0)
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
  pgsch (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].char_height);       // Set character height
  pgbox (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].gridx[p], xtick, nxsub, ...
         PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].gridy[p], ytick, nysub);

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

  pgsci (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].color[p;1]);       // Set line color.
  pgslw (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].width[p]);         // Line width.
  pgscf (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].font[p]);          // Set character font.
  pgsch (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].char_height);      // Set character height.
  pgsfs (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].fillstyle[p]);     // Set the fillstyle.

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

  pglab(PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].xlabel[p], PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].ylabel[p], PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].title[p]);
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
  p = mod (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].subplot, PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].nplot) + 1;	// The current index

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

  if (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].xmin[p] != 1j) { xmin = real (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].xmin[p]); }
  if (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].xmax[p] != 1j) { xmax = real (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].xmax[p]); }

  if (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].ymin[p] != 1j) { ymin = real (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].ymin[p]); }
  if (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].ymax[p] != 1j) { ymax = real (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].ymax[p]); }

  pgslw (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].width[p]);         // Line width
  pgscf (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].font[p]);          // Set character font
  pgsch (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].char_height);      // Set character height

  pgsls (1);
  pgsci (1);

  pgenv (xmin, xmax, ymin, ymax, 0, -2);

  pgsci (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].color[p;1]);     // Set line color.
  pgsls (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].lstyle[p;1]);    // Set line style.
  pgslw (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].width[p;1]);     // Set line width.

  pgline (V.nr, V[;1], V[;2]);

  pgsls (1);
  pgsci (1);

  pglab(PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].xlabel[p], PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].ylabel[p], PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].title[p]);
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
  p = mod (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].subplot, PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].nplot) + 1;	// The current index

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
    tmp = rotate3 ( PNTS[i;], theta*dtor );

    // Now do the second rotation...
    tmp = rotate1 ( tmp, phi*dtor );

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

  if (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].xmin[p] != 1j) { xmin = real (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].xmin[p]); }
  if (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].xmax[p] != 1j) { xmax = real (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].xmax[p]); }

  if (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].ymin[p] != 1j) { ymin = real (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].ymin[p]); }
  if (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].ymax[p] != 1j) { ymax = real (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].ymax[p]); }

  if (axis)
  {
    // Rotate the axes...
    tmp = rotate3 ( xmax * axis3d[1;], theta*dtor );
    ax1 = rotate1 ( tmp, phi*dtor );

    tmp = rotate3 ( ymax * axis3d[2;], theta*dtor );
    ax2 = rotate1 ( tmp, phi*dtor );

    tmp = rotate3 ( ymax * axis3d[3;], theta*dtor );
    ax3 = rotate1 ( tmp, phi*dtor );
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

    if (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].aspect[p] != 0)
    {
      pgwnad (xmin, xmax, ymin, ymax);
      } else {
      pgswin (xmin, xmax, ymin, ymax);
    }

    pgbox (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].gridx[p], 0, 0, PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].gridy[p], 0, 0);
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

  pgsci (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].color[p;1]);       // Set line color.
  pgslw (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].width[p]);         // Line width
  pgscf (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].font[p]);          // Set character font
  pgsch (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].char_height);      // Set character height

  //
  // Now, plot up the projected surface.
  //

  for (i in 1:pnts.nr)
  {
    pgpt (pnts.nr, pnts[;1], pnts[;2], PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].pstyle[p;1]);
    // pgpt (pnts.nr, pnts[;1], pnts[;2], 55);
  }

  pgsls (1);
  pgsci (1);

  pglab(PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].xlabel[p], PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].ylabel[p], PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].title[p]);
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
  p = mod (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].subplot, PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].nplot) + 1;	// The current index

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
    tmp = rotate3 ( PNTS[i;], theta*dtor );

    // Now do the second rotation...
    tmp = rotate1 ( tmp, phi*dtor );

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

  if (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].xmin[p] != 1j) { xmin = real (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].xmin[p]); }
  if (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].xmax[p] != 1j) { xmax = real (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].xmax[p]); }

  if (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].ymin[p] != 1j) { ymin = real (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].ymin[p]); }
  if (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].ymax[p] != 1j) { ymax = real (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].ymax[p]); }

  if (axis)
  {
    // Rotate the axes...
    tmp = rotate3 ( xmax * axis3d[1;], theta*dtor );
    ax1 = rotate1 ( tmp, phi*dtor );

    tmp = rotate3 ( ymax * axis3d[2;], theta*dtor );
    ax2 = rotate1 ( tmp, phi*dtor );

    tmp = rotate3 ( ymax * axis3d[3;], theta*dtor );
    ax3 = rotate1 ( tmp, phi*dtor );
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

    if (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].aspect[p] != 0)
    {
      pgwnad (xmin, xmax, ymin, ymax);
      } else {
      pgswin (xmin, xmax, ymin, ymax);
    }

    pgbox (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].gridx[p], 0, 0, PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].gridy[p], 0, 0);
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

  pgsci (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].color[p;1]);       // Set line color.
  pgslw (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].width[p]);         // Line width
  pgscf (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].font[p]);          // Set character font
  pgsch (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].char_height);      // Set character height

  //
  // Now, plot up the projected line segments.
  //

  pgline (pnts.nr, pnts[;1], pnts[;2]);

  pgsls (1);
  pgsci (1);

  pglab(PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].xlabel[p], PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].ylabel[p], PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].title[p]);
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
    V[i;] = rotate3 ( V[i;], THETA*dtor )';
  }

  //
  // Do X-axis rotation...
  //

  for (i in 1:V.nr)
  {
    V[i;] = rotate1 ( V[i;], PHI*dtor )';
  }

  return V;
};

//---------------------------------------------------------------------
// Rotate about the 1st axis (X)
//

rotate1 = function ( V, TH )
{
  V = V[:];
  a = [ 1,      0,        0;
        0,   cos (TH), sin (TH);
        0,  -sin (TH), cos (TH)];

  return (a * V);
};

//---------------------------------------------------------------------
// Rotate about the 2st axis (Y)
//

rotate2 = function ( V, TH )
{
  V = V[:];
  a = [ cos (TH), 0,  sin (TH);
            0   , 1,     0;
       -sin (TH), 0,  cos (TH)];

  return (a * V);
};

//---------------------------------------------------------------------
// Rotate about the 2st axis (Z)
//

rotate3 = function ( V, TH )
{
  V = V[:];
  a = [ cos (TH), sin (TH), 0;
       -sin (TH), cos (TH), 0;
           0    ,     0,    1];

  return (a * V);
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
  i = mod (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].subplot, PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].nplot) + 1;

  if (exist (ALT))
  {
    PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].alt[i] = ALT;
    } else {
      PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].alt[i] = 60;
  }
  return PGPLOT_ACTIVE_WIN;
};

##############################################################################
#
# Set the viewing azimuth for 3-D plots
#

plaz = function ( AZ )
{
  check_plot_object ();
  i = mod (PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].subplot, PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].nplot) + 1;

  if (exist (AZ))
  {
    PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].az[i] = AZ;
  }
  else
  {
    PGPLOT_PLOT_WINDOWS.[PGPLOT_ACTIVE_WIN].az[i] = 45;
  }
  return PGPLOT_ACTIVE_WIN;
};





