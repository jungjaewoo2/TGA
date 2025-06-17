//  xmgrace interface for the RLaB2
//  rlabplus (c) 2000-2006 marijan kostrun
//  please report all problems to mkostrun@gmail.com
//
//  This file is a part of rlabplus
//
//  This program is free software; you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation; either version 2 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with this program; if not, write to the Free Software
//  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
//
//  See the file ./COPYING
//**********************************************************************


//
// Objects associated with the page
//
static( _xmgobjects, _xmgobjectypes );
_xmgobjectypes = [ ...
  "line", ...
  "string" ...
];
_xmgobjects=<<>>;

xmgobjectsclear = function()
{
  for (i in members( _xmgobjects ))
  {
    clear( _xmgobjects.[ i ] );
    _xmgobjects.[ i ] = <<>>;
    _xmgobjects.[ i ].cnt = 0;
  }
};

//
// write a graffiti across the xmg page or a graph
//
xmgtext = function( a1, a2, a3, a4, a5 )
{
  // a1 -> text
  // a2 -> location [x,y]
  // a3 -> size
  // a4 -> parameters: color, rotation, font, justification
  // a5 -> graph (if given then location is in graph coordinates, otherwise in world's)
  if ( !exist(a1) ) { return 0; }
  if (type(a1) != "string") { return 0; }
  if (length(a1) == 0) { return 0; }

  if ( !exist(_xmgobjects.["string"]) )
  {
    _xmgobjects.["string"] = <<>>;
    _xmgobjects.["string"].cnt = 0;
  }

  id = _xmgobjects.["string"].cnt;
  _xmgobjects.["string"].cnt = id + 1;
  _xmgobjects.["string"].[ id ] = <<>>;
  //
  // world or view from graph info
  //
  if (exist (a5))
  {
    if(type(a5)=="real" && a5.nr * a5.nc == 1)
    {
      _xmgobjects.["string"].[ id ].loctype = "world";
      _xmgobjects.["string"].[ id ].graph   = "g" + text(a5 + 0);
    }
  } else {
    _xmgobjects.["string"].[ id ].loctype = "view";
  }
  //
  // location
  //
  if (exist(a2))
  {
    if( a2.nr * a2.nc == 2)
    {
      _xmgobjects.["string"].[ id ].location = [ a2[1], a2[2] ];
    }
  } else {
    _xmgobjects.["string"].[ id ].location = [ 0.5, 0.5 ];
  }

  //
  // size
  //
  if (exist(a3))
  {
    if( a3.nr * a3.nc >= 1)
    {
      _xmgobjects.["string"].[ id ].size = a3[1];
    }
  } else {
    _xmgobjects.["string"].[ id ].size = 1.5;
  }

  //
  // parameters: color, rotation, font, justification
  //
  if (exist(a4))
  {
    if( a4.nr * a4.nc >= 1)
    {
      _xmgobjects.["string"].[ id ].color = a4[1];
    } else {
      _xmgobjects.["string"].[ id ].color = 1;
    }
    if( a4.nr * a4.nc >= 2)
    {
      _xmgobjects.["string"].[ id ].rot = a4[2];
    } else {
      _xmgobjects.["string"].[ id ].rot = 0;
    }
    if( a4.nr * a4.nc >= 3)
    {
      _xmgobjects.["string"].[ id ].font = a4[3];
    } else {
      _xmgobjects.["string"].[ id ].font = 0;
    }
    if( a4.nr * a4.nc >= 4)
    {
      _xmgobjects.["string"].[ id ].just = a4[4];
    } else {
      _xmgobjects.["string"].[ id ].just = 0;
    }
  } else {
    _xmgobjects.["string"].[ id ].color = 1;
    _xmgobjects.["string"].[ id ].rot   = 0;
    _xmgobjects.["string"].[ id ].font  = 0;
    _xmgobjects.["string"].[ id ].just  = 0;
  }

  _xmgobjects.["string"].[ id ].string = "\"" + a1 + "\"";
  return 1;
};

//
// draw a line across xmg page
//
xmgline = function( a1, a2, a3, a4, a5 )
{
  // a1 -> location [x1,y1, x2, y2]
  // a2 -> linewidth
  // a3 -> parameters: color, linestyle
  // a4 -> arrowmodes: arrow, atype, alength, alayout
  // a5 -> graph (if given then location is in graph coordinates, otherwise in world's)
  if ( !exist(a1) ) { return 0; }
  if ( type(a1) != "real" ) { return 0; }
  if ( a1.nr * a1.nc != 4 ) { return 0; }

  if ( !exist(_xmgobjects.["line"]) )
  {
    _xmgobjects.["line"] = <<>>;
    _xmgobjects.["line"].cnt = 0;
  }

  id = _xmgobjects.["line"].cnt;
  _xmgobjects.["line"].cnt    = id + 1;
  _xmgobjects.["line"].[ id ] = <<>>;
  //
  // world or view from graph info
  //
  if (exist (a5))
  {
    if(type(a5)=="real" && a5.nr * a5.nc == 1)
    {
      _xmgobjects.["line"].[ id ].loctype = "world";
      _xmgobjects.["line"].[ id ].graph   = "g" + text(a5 + 0);
    }
  } else {
    _xmgobjects.["line"].[ id ].loctype = "view";
  }

  //
  // location
  //
  _xmgobjects.["line"].[ id ].location = a1;

  //
  // linewidth
  //
  if (exist(a2))
  {
    if( a2.nr * a2.nc >= 1)
    {
      _xmgobjects.["line"].[ id ].width = a2[1];
    }
  } else {
    _xmgobjects.["line"].[ id ].width = 1.5;
  }

  //
  // parameters: color, linestyle
  //
  if (exist(a3))
  {
    if( a3.nr * a3.nc >= 1)
    {
      _xmgobjects.["line"].[ id ].color = a3[1];
    } else {
      _xmgobjects.["line"].[ id ].color = 1;
    }
    if( a3.nr * a3.nc >= 2)
    {
      _xmgobjects.["line"].[ id ].style = a3[2];
    } else {
      _xmgobjects.["line"].[ id ].style = 1;
    }
  } else {
    _xmgobjects.["line"].[ id ].color = 1;
    _xmgobjects.["line"].[ id ].style = 1;
  }

  //
  // arrow: arrow, atype, alength, alayout
  //
  if (exist(a4))
  {
    if( a4.nr * a4.nc >= 1)
    {
      _xmgobjects.["line"].[ id ].arrow = a4[1];
    } else {
      _xmgobjects.["line"].[ id ].arrow = 0;
    }
    if( a4.nr * a4.nc >= 2)
    {
      _xmgobjects.["line"].[ id ].atype = a4[2];
    } else {
      _xmgobjects.["line"].[ id ].atype = 0;
    }
    if( a4.nr * a4.nc >= 3)
    {
      _xmgobjects.["line"].[ id ].alen = a4[3];
    } else {
      _xmgobjects.["line"].[ id ].alen = 1;
    }
    if( a4.nr * a4.nc >= 5)
    {
      _xmgobjects.["line"].[ id ].alayout = a4[4:5];
    } else {
      _xmgobjects.["line"].[ id ].alayout = [1,1];
    }
    return 1;
  }

  //
  // default no-arrow
  //
  _xmgobjects.["line"].[ id ].arrow   = 0;
  _xmgobjects.["line"].[ id ].atype   = 0;
  _xmgobjects.["line"].[ id ].alen    = 1;
  _xmgobjects.["line"].[ id ].alayout = [1,1];

  return 1;
};


//
// initialize some default and control values
//
static( _xmgplottype, _xmgncol, _xmgcolormap);
_xmgplottype = [ ...
    "xy";...
    "xydx";...
    "xydy";...
    "xydxdx";...
    "xydydy";...
    "xydxdy";...
    "xydxdxdydy";...
    "bar";...
    "bardy";...
    "bardydy";...
    "xyhilo"...
    ];
_xmgncol = [ ...
    2;...
    3;...
    3;...
    4;...
    4;...
    4;...
    6;...
    1;...
    2;...
    3;...
    3 ...
    ];

//
// Default color map. Use function xmgcolormap to modify
//
_xmgcolormap=<<>>;

for(i in 0:15)
{ _xmgcolormap.[i] = <<>>; }

_xmgcolormap.[0].color = "white";
_xmgcolormap.[0].rgb   = [255, 255, 255];
_xmgcolormap.[1].color = "black";
_xmgcolormap.[1].rgb   = [0, 0, 0];
_xmgcolormap.[2].color = "red";
_xmgcolormap.[2].rgb   = [255, 0, 0];
_xmgcolormap.[3].color = "green";
_xmgcolormap.[3].rgb   = [0, 255, 0];
_xmgcolormap.[4].color = "blue";
_xmgcolormap.[4].rgb   = [0, 0, 255];
_xmgcolormap.[5].color = "yellow";
_xmgcolormap.[5].rgb   = [255, 255, 0];
_xmgcolormap.[6].color = "brown";
_xmgcolormap.[6].rgb   = [188, 143, 143];
_xmgcolormap.[7].color = "grey";
_xmgcolormap.[7].rgb   = [153, 153, 153];
_xmgcolormap.[8].color = "violet";
_xmgcolormap.[8].rgb   = [148, 0, 211];
_xmgcolormap.[9].color = "cyan";
_xmgcolormap.[9].rgb   = [0, 255, 255];
_xmgcolormap.[10].color = "magenta";
_xmgcolormap.[10].rgb   = [255, 0, 255];
_xmgcolormap.[11].color = "orange";
_xmgcolormap.[11].rgb   = [255, 165, 0];
_xmgcolormap.[12].color = "indigo";
_xmgcolormap.[12].rgb   = [114, 33, 188];
_xmgcolormap.[13].color = "maroon";
_xmgcolormap.[13].rgb   = [103, 7, 72];
_xmgcolormap.[14].color = "turquoise";
_xmgcolormap.[14].rgb   = [64, 224, 208];
_xmgcolormap.[15].color = "green5";       // replaces green4 which on the color laser printout
_xmgcolormap.[15].rgb   = [0, 59, 0];     // could not be distinguished from standard green


//
// xmgr service functions
//
xmgcolormap = function ( i, color, rgb )
{
  if(!exist(i))
  { error("xmgcolormap: need 'i' an integer, index of the color"); }

  if(exist(i)&&!exist(color)&&!exist(rgb))
  {
    if(exist(_xmgcolormap.[i]))
    { return _xmgcolormap.[i];}

    return <<>>;
  }

  if(type(color)!="string")
  { error("xmgcolormap: 'color' is not a string"); }

  if(type(rgb)!="real")
  { error("xmgcolormap: 'rgb' is a row-vector of the form [r,g,b] !"); }

  if(rgb.nr==1 && rgb.nc==3)
  {
    _xmgcolormap.[i] = <<>>;
    _xmgcolormap.[i].color = color;
    _xmgcolormap.[i].rgb   = rgb;
    return 1;
  }

  return 0;

};  // xmgcolormap

//
// initialize default xmg dataset, graph and paper, and their counters
//
static(xmg_paper);     xmg_paper   = <<>>;
static(xmg_graph);     xmg_graph   = <<>>;
static(xmg_dataset);   xmg_dataset = <<>>;
static(DefXMGdataset); DefXMGdataset = 0;
static(DefXMGgraph);   DefXMGgraph   = 0;
static(DefXMGpaper);   DefXMGpaper   = 0;

//
// initialize functions operating on xmg objects dataset, graph, paper:
//  defxmgpaper( P )
//  xmgpage_list( P )
//  defxmggraph( N )
//  xmggraph_list( N )
//  defxmgdataset( I )
//  xmgdataset_list( I )

defxmgpaper = function( P, _reset_ )
{
  if(exist(_reset_) && exist(xmg_paper.[P]))
  { clear(xmg_paper.[P]); }

  // set a default xmg_paper object
  if (exist(xmg_paper.[P]))
  {
    DefXMGpaper = P;
    return 1;

  } else {

    xmg_paper.[P]=<<>>;

    // default configuration
    xmg_paper.[P].papersize   = "letter";
    xmg_paper.[P].orientation = "landscape";
    xmg_paper.[P].graphs      = [ 0 ];
    xmg_paper.[P].filename    = "./tmpgraph.gr" ;
    DefXMGpaper = P;

    return 1;
  }

};  // defxmgpaper

defxmggraph = function( N, _reset_ )
{
  if(exist(_reset_) && exist(xmg_graph.[N]))
  { clear(xmg_graph.[N]); }

  // set a default xmg_graph object
  if (exist(xmg_graph.[N]))
  {
    DefXMGgraph = N;

  } else {

    xmg_graph.[N]=<<>>;

    // default configuration
    xmg_graph.[N].title        = "";
    xmg_graph.[N].titlesize    = 1.3;
    xmg_graph.[N].subtitle     = "";
    xmg_graph.[N].subtitlesize = 1.15;
    xmg_graph.[N].xlabel       = "X-axis";
    xmg_graph.[N].xlabelsize   = 1.0;
    xmg_graph.[N].xscale       = "Normal";
    xmg_graph.[N].xticks       = zeros(1,2);
    xmg_graph.[N].xticksize    = 1.0;
    xmg_graph.[N].xtickspec    = <<>>;
    xmg_graph.[N].ylabel       = "Y-axis";
    xmg_graph.[N].ylabelsize   = 1.0;
    xmg_graph.[N].yscale       = "Normal";
    xmg_graph.[N].yticks       = zeros(1,2);
    xmg_graph.[N].yticksize    = 1.0;
    xmg_graph.[N].ytickspec    = <<>>;
    xmg_graph.[N].xmax         = 1;
    xmg_graph.[N].xmin         = 0;
    xmg_graph.[N].ymax         = 1;
    xmg_graph.[N].ymin         = 0;
    xmg_graph.[N].legend       = <<>>;
    //xmg_graph.[N].legend.x     = 0.5;
    //xmg_graph.[N].legend.y     = 0.5;
    //xmg_graph.[N].legend.size  = 1.0;
    clear(xmg_graph.[N].xgrid);
    clear(xmg_graph.[N].ygrid);
    xmg_graph.[N].datasets     = [ 0 ];
    xmg_graph.[N].size         = [0,1,0,1];
    DefXMGgraph = N;
  }
};  // defxmggraph

defxmgdataset = function( I, _reset_ )
{
  if( exist(_reset_) && exist(xmg_dataset.[I]) )
  {
    clear(xmg_dataset.[I]);
  }
  // set a default xmg_dataset object
  if (exist(xmg_dataset.[I]))
  {
    DefXMGdataset = I;
  } else {
    xmg_dataset.[I] = <<>>;
    // default configuration:
    // line
    xmg_dataset.[I].dataset    = [];
    xmg_dataset.[I].dataname   = "";
    xmg_dataset.[I].plottype   = "";
    xmg_dataset.[I].linetype   = 0;
    xmg_dataset.[I].lineshape  = 1;
    xmg_dataset.[I].linecolor  = 1;
    xmg_dataset.[I].linesize   = 3;
    // symbol
    xmg_dataset.[I].symbshape  = 1;
    xmg_dataset.[I].symblwid   = 1;
    xmg_dataset.[I].symbcolor  = 1;
    xmg_dataset.[I].symbsize   = 0.75;
    xmg_dataset.[I].xmax       = 1;
    xmg_dataset.[I].xmin       = 0;
    xmg_dataset.[I].ymax       = 1;
    xmg_dataset.[I].ymin       = 0;
    DefXMGdataset = I;
  }
};

xmgpage_list = function ( P ){
    if(exist(P)){
        // return the xmg_paper object P
        if( exist(xmg_paper.[P]) ){
            return xmg_paper.[P];
        } else {{
            return <<>>;
            }}
        }
    return who(xmg_paper);
    };

xmggraph_list = function ( N ){
    if(exist(N)){
        // return the xmg_graph object N
        if( exist(xmg_graph.[N]) ){
            return xmg_graph.[N];
        } else {{
            return <<>>;
            }}
        }
    return who(xmg_graph);
    };

xmgdataset_list = function ( I )
{
  if(exist(I))
  {
    // Set the default xmg-object
    if( exist(xmg_dataset.[I]) )
    {
      return xmg_dataset.[I];
    } else {
    {
      return <<>>;
    }}
  }
  return who(xmg_dataset);
};

//
// functions operating on xmg_paper list
//

xmgpaper = function( papersize, orientation, P ){
    if (!exist(P)){P=DefXMGpaper;}
    defxmgpaper(P);
    if(!exist(xmg_paper.[P].papersize))  { xmg_paper.[P].papersize   = "letter"; }
    if(!exist(xmg_paper.[P].orientation)){ xmg_paper.[P].orientation = "landscape"; }
    if(exist(papersize)){
        if(papersize == "A4")    { xmg_paper.[P].papersize = "A4"; }
        if(papersize == "letter"){ xmg_paper.[P].papersize = "letter"; }
        }
    if(exist(orientation)){
        if(orientation == "landscape"){ xmg_paper.[P].orientation = "landscape"; }
        if(orientation == "portrait") { xmg_paper.[P].orientation = "portrait"; }
        }
    return 1;
    };

xmgfilename = function ( filename, P ){
    if (!exist(P)){P=DefXMGpaper;}
    defxmgpaper(P);
    if(!exist(xmg_paper.[P].filename)){ xmg_paper.[P].filename="./tmpgraph.gr"; }
    if (type(filename)=="string"){
        xmg_paper.[P].filename=filename;
        return 1;
        }
    if(!exist(filename)){
        xmg_paper.[P].filename="./tmpgraph.gr";
        return 1;
        }
    return 0;
    };


//
// functions operating on xmg_graph list
//

xmgtitle = function ( titletext, titlesize, N ){
    if (!exist(N)){N=DefXMGgraph;}
    defxmggraph(N);
    if(!exist(xmg_graph.[N].title))    { xmg_graph.[N].title     = ""; }
    if(!exist(xmg_graph.[N].titlesize)){ xmg_graph.[N].titlesize = 1.3; }
    if (type(titletext)=="string"){ xmg_graph.[N].title     = titletext; }
    if (type(titlesize)=="real")  { xmg_graph.[N].titlesize = titlesize; }
    return 1;
    };

xmgsubtitle = function ( subtitletext, subtitlesize, N ){
    if (!exist(N)){N=DefXMGgraph;}
    defxmggraph(N);
    if(!exist(xmg_graph.[N].subtitle))    { xmg_graph.[N].subtitle     = ""; }
    if(!exist(xmg_graph.[N].subtitlesize)){ xmg_graph.[N].subtitlesize = 1.15; }
    if (type(subtitletext)=="string"){ xmg_graph.[N].subtitle     = subtitletext; }
    if (type(subtitlesize)=="real")  { xmg_graph.[N].subtitlesize = subtitlesize; }
    return 1;
    };

xmgxylabels = function ( labels, labelsize, N )
{
  if( type(labels)!="string" )
  {
    return 0;
  }
  if (!exist(N))
  {
    N=DefXMGgraph;
  }
  defxmggraph(N);
  // set the default values, and then get the new labels
  if( !exist(xmg_graph.[N].xlabel) )
  {
    xmg_graph.[N].xlabel = "";
  }
  if( !exist(xmg_graph.[N].ylabel) )
  {
    xmg_graph.[N].ylabel = "";
  }
  if (labels.nr * labels.nc >= 1)
  {
    xmg_graph.[N].xlabel = labels[1];
  }
  if (labels.nr * labels.nc >= 2)
  {
    xmg_graph.[N].ylabel = labels[ 2 ];
  }
  // set the defaults, and then get the new label sizes
  if( !exist(xmg_graph.[N].xlabelsize) )
  {
    xmg_graph.[N].xlabelsize = 1.0;
  }
  if( !exist(xmg_graph.[N].ylabelsize) )
  {
    xmg_graph.[N].ylabelsize = 1.0;
  }
  if (labelsize.nr * labelsize.nc == 1)
  {
    xmg_graph.[N].xlabelsize = labelsize[1];
    xmg_graph.[N].ylabelsize = labelsize[1];
  }
  if (labelsize.nr * labelsize.nc >= 2)
  {
    xmg_graph.[N].xlabelsize = labelsize[1];
    xmg_graph.[N].ylabelsize = labelsize[2];
  }
  return 1;
};

xmgxyscale = function ( scale, N ){
    if (!exist(N)){N=DefXMGgraph;}
    defxmggraph(N);
    if(!exist(xmg_graph.[N].xscale)){ xmg_graph.[N].xscale = "normal"; }
    if(!exist(xmg_graph.[N].yscale)){ xmg_graph.[N].yscale = "normal"; }
    if (type(scale)!="string"){ return 0; }
    if (scale.nr*scale.nc!=2 && scale.nr*scale.nc!=1){ return 0; }
    if (scale[1]=="Log"    || scale[1]=="log")   { scale[1]="Logarithmic"; }
    if (scale[1]=="normal" || scale[1]=="Normal"){ scale[1]="normal"; }
    if (scale.nr*scale.nc==2){
        if (scale[2]=="Log"    || scale[2]=="log")   { scale[2]="Logarithmic"; }
        if (scale[2]=="normal" || scale[2]=="Normal"){ scale[2]="normal"; }
        }
    xmg_graph.[N].xscale = scale[1];
    xmg_graph.[N].yscale = scale[ min(2,max(size(scale))) ];
    return 1;
    };

xmgxticks = function ( ticks, ticksize, N )
{
  if( type(ticks) != "real" )
  { return 0; }
  if ( !exist(N) )
  { N=DefXMGgraph; }
  defxmggraph(N);

  xmg_graph.[N].xticks = zeros(1,2);
  if ( exist(ticks) )
  {
    if ( ticks.nr*ticks.nc >= 1)
    { xmg_graph.[N].xticks[1] = ticks[1]; }
    if ( ticks.nr*ticks.nc >= 2)
    { xmg_graph.[N].xticks[2] = ticks[2]; }
  }

  xmg_graph.[N].xticksize = 1.0;
  if ( exist(ticksize) )
  {
    if ( type(ticksize)=="real" )
    { xmg_graph.[N].xticksize = ticksize[1]; }
  }

  return 1;
};

xmgyticks = function ( ticks, ticksize, N )
{
  if( type(ticks) != "real" )
  { return 0; }
  if ( !exist(N) )
  { N=DefXMGgraph; }
  defxmggraph(N);

  xmg_graph.[N].yticks = zeros(1,2);
  if ( exist(ticks) )
  {
    if ( ticks.nr*ticks.nc >= 1)
    { xmg_graph.[N].yticks[1] = ticks[1]; }
    if ( ticks.nr*ticks.nc >= 2)
    { xmg_graph.[N].yticks[2] = ticks[2]; }
  }

  xmg_graph.[N].yticksize = 1.0;
  if ( exist(ticksize) )
  {
    if ( type(ticksize) == "real" )
    { xmg_graph.[N].yticksize = ticksize[1]; }
  }
  return 1;
};

xmgxtickspec = function ( a1, a2, a3 )
{
  if (!exist(a1) || !exist(a2))
  { return 0; }
  if( type(a1) != "real" || type(a2) == "complex")
  { return 0; }

  tick  = a1;
  label = blank(0,0);
  if (exist(a3))
  {
    label = a2;
    N     = a3[1];
    if (length(label)!= length(tick))
    {
      fprintf("xmgxtickspec: mismatch between 'tick' and 'label' sizes");
      return 0;
    }

  } else {
    N = a2[1];
  }

  if ( !exist(N) )
  { N=DefXMGgraph; }
  defxmggraph(N);

  xmg_graph.[N].xtickspec.tick = tick;
  if(length(label) > 0)
  {
    xmg_graph.[N].xtickspec.kind = "both";
    xmg_graph.[N].xtickspec.label = label;
  } else {
    xmg_graph.[N].xtickspec.kind = "ticks";
  }

  return 1;
};

xmgytickspec = function ( a1, a2, a3 )
{
  if (!exist(a1) || !exist(a2))
  { return 0; }
  if( type(a1) != "real" || type(a2) == "complex")
  { return 0; }

  tick  = a1;
  label = blank(0,0);
  if (exist(a3))
  {
    label = a2;
    N     = a3[1];
    if (length(label)!= length(tick))
    {
      fprintf("xmgytickspec: mismatch between 'tick' and 'label' sizes");
      return 0;
    }

    } else {
      N = a2[1];
  }

  if ( !exist(N) )
  { N=DefXMGgraph; }
  defxmggraph(N);

  xmg_graph.[N].ytickspec.tick = tick;
  if(length(label) > 0)
  {
    xmg_graph.[N].ytickspec.kind = "both";
    xmg_graph.[N].ytickspec.label = label;
    } else {
      xmg_graph.[N].ytickspec.kind = "ticks";
  }

  return 1;
};


xmgrange = function( range, N )
{
  if (!exist(N)){N=DefXMGgraph;}
  defxmggraph(N);
  if(!exist(xmg_graph.[N].xmax))
  { xmg_graph.[N].xmax = 1;}
  if(!exist(xmg_graph.[N].xmin))
  { xmg_graph.[N].xmin = 0;}
  if(!exist(xmg_graph.[N].ymax))
  { xmg_graph.[N].ymax = 1;}
  if(!exist(xmg_graph.[N].ymin))
  { xmg_graph.[N].ymin = 0;}
  if(type(range)!="real")
  { return 0;}
  if(range.nr*range.nc!=4)
  { return 0;}
  xmg_graph.[N].xmin=range[1];
  xmg_graph.[N].xmax=range[2];
  xmg_graph.[N].ymin=range[3];
  xmg_graph.[N].ymax=range[4];
  return 1;
};

xmglegend = function ( legendloc, legendsize, lopts, N )
{
  // default legend parameters that user might want to change
  _color  = 1;
  _length = 4;
  _vgap   = 2;
  _hgap   = 1;

  // user forgot to explicitely omit lopts
  if (class(lopts)=="num" && !exist(N))
  { N = lopts; }

  // does graph N exists?
  if (!exist(N))
  {N=DefXMGgraph;}
  defxmggraph(N);

  if(!exist(legendloc)||!exist(legendsize))
  {
    xmg_graph.[N].legend = <<>>;
    return 1;
  }
  if(type(legendloc)!="real")
  { return 0; }
  if(legendloc.nr*legendloc.nc != 2)
  { return 0; }

  //
  if (class(lopts) == "list")
  {
    if (exist(lopts.color))
    { _color = lopts.color; }
    if (exist(lopts.length))
    { _length = lopts.length; }
    if (exist(lopts.vgap))
    { _vgap = lopts.vgap; }
    if (exist(lopts.hgap))
    { _hgap = lopts.hgap; }
  }
  xmg_graph.[N].legend   = <<>>;
  xmg_graph.[N].legend.x    = legendloc[1];
  xmg_graph.[N].legend.y    = legendloc[2];
  xmg_graph.[N].legend.size = legendsize;
  xmg_graph.[N].legend.color  = _color;
  xmg_graph.[N].legend.length = _length;
  xmg_graph.[N].legend.vgap   = _vgap;
  xmg_graph.[N].legend.hgap   = _hgap;
  return 1;
};

xmgxgrid = function( gridmaj,gridmin, N ){
    if (!exist(N)){N=DefXMGgraph;}
    defxmggraph(N);
    if(!exist(gridmaj)){
        clear(xmg_graph.[N].xgrid);
        return 1;
        }
    if (type(gridmaj)!="real")  { return 0; }
    if (gridmaj.nr*gridmaj.nc!=3){ return 0; }
    xmg_graph.[N].xgrid=<<>>;
    xmg_graph.[N].xgrid.major=<<>>;
    xmg_graph.[N].xgrid.major.color = gridmaj[1];
    xmg_graph.[N].xgrid.major.width = gridmaj[2];
    xmg_graph.[N].xgrid.major.style = gridmaj[3];
    if (gridmin.nr*gridmin.nc==3 && type(gridmin)=="real"){
        xmg_graph.[N].xgrid.minor=<<>>;
        xmg_graph.[N].xgrid.minor.color = gridmin[1];
        xmg_graph.[N].xgrid.minor.width = gridmin[2];
        xmg_graph.[N].xgrid.minor.style = gridmin[3];
        }
    return 1;
    };

xmgygrid = function( gridmaj, gridmin, N ){
    if (!exist(N)){N=DefXMGgraph;}
    defxmggraph(N);
    if(!exist(gridmaj)){
        clear(xmg_graph.[N].ygrid);
        return 1;
        }
    if (type(gridmaj)!="real")  { return 0; }
    if (gridmaj.nr*gridmaj.nc!=3){ return 0; }
    xmg_graph.[N].ygrid=<<>>;
    xmg_graph.[N].ygrid.major=<<>>;
    xmg_graph.[N].ygrid.major.color = gridmaj[1];
    xmg_graph.[N].ygrid.major.width = gridmaj[2];
    xmg_graph.[N].ygrid.major.style = gridmaj[3];
    if (gridmin.nr*gridmin.nc==3 && type(gridmin)=="real"){
        xmg_graph.[N].ygrid.minor=<<>>;
        xmg_graph.[N].ygrid.minor.color = gridmin[1];
        xmg_graph.[N].ygrid.minor.width = gridmin[2];
        xmg_graph.[N].ygrid.minor.style = gridmin[3];
        }
    return 1;
    };

xmggraphsize = function( graphsize, N ){
    if (!exist(N)){N=DefXMGgraph;}
    defxmggraph(N);
    if(!exist(graphsize)){
        xmg_graph.[N].size = [0,1,0,1];
        return 1;
        }
    if (type(graphsize)!="real"){ return 0; }
    if (graphsize.nr*graphsize.nc!=4){ return 0; }
    xmg_graph.[N].size = graphsize;
    return 1;
    };


//
// functions operating on xmg_dataset list
//

xmgdataset = function ( dataset, plottype, dataname, I )
{
  if( isempty(dataset) )
  {
    printf("xmgdataset: dataset is empty!\n");
    return 0;
  }
  if(type(plottype)!="string")
  {
    return 0;
  }
  if(sum(plottype == _xmgplottype)!=1)
  {
    return 0;
  }
  if (!exist(I))
  {
    I=DefXMGdataset;
  }
  defxmgdataset(I);
  if(type(dataset)!="real")
  {
    printf("Invalid dataset: removing dataset %i from all graphs!\n");
    // remove dataset from the list
    clear(xmg_dataset.[I]);
    // remove dataset from all xmg_graph lists
    allgraphs = who(xmg_graph);
    for(graph in allgraphs)
    {
      if(sum(xmg_graph.[graph].datasets == I)>0)
      {
        newds = [];
        for(i in xmg_graph.[graph].datasets)
        {
          if(i != I)
          {
            newds = [newds, i];
          }
        }
        xmg_graph.[graph].datasets = newds;
      }
    }
    return 1;
  }
  xmg_dataset.[I].dataset = dataset;
  i = find( plottype == _xmgplottype );

  if(_xmgncol[ i[1] ]!= dataset.nc)
  {
    error("Invalid size of dataset: number of columns does not match its type!\n");
    return 0;
  }

  if(!exist(xmg_dataset.[I].linetype))
  { xmg_dataset.[I].linetype  = 0; }
  if(!exist(xmg_dataset.[I].lineshape))
  { xmg_dataset.[I].lineshape = 1; }
  if(!exist(xmg_dataset.[I].linecolor))
  { xmg_dataset.[I].linecolor = 1; }
  if(!exist(xmg_dataset.[I].linesize))
  { xmg_dataset.[I].linesize  = 3; }

  if(plottype != "bar" && plottype != "bardy" && plottype != "bardydy")
  {
    if(!exist(xmg_dataset.[I].xmax))
    { xmg_dataset.[I].xmax = max( xmg_dataset.[I].dataset[;1] ); }

    if(!exist(xmg_dataset.[I].xmin))
    { xmg_dataset.[I].xmin = min( xmg_dataset.[I].dataset[;1] ); }

    if(!exist(xmg_dataset.[I].ymax))
    { xmg_dataset.[I].ymax = max( xmg_dataset.[I].dataset[;2] ); }

    if(!exist(xmg_dataset.[I].ymin))
    { xmg_dataset.[I].ymin = min( xmg_dataset.[I].dataset[;2] ); }
  }
  xmg_dataset.[I].plottype = plottype;

  // name
  if (strlen(dataname) > 0)
  {
    xmg_dataset.[I].dataname = "\"" + dataname + "\"";
  } else {
    xmg_dataset.[I].dataname = "";
  }

  return 1;
};

xmgstyle = function ( linesymbol, setstyle, I )
{

  if( type(linesymbol) != "string" )
  { return 0; }
  if( isempty(setstyle) )
  { setstyle = []; }
  if( type(setstyle) != "real" )
  { return 0; }

  if (!exist(I))
  { I=DefXMGdataset; }
  defxmgdataset(I);
  xmg_dataset.[I].linetype   = 0;
  if(linesymbol == "line")
  { xmg_dataset.[I].linetype   = 0;}
  if(linesymbol == "symbol")
  { xmg_dataset.[I].linetype   = 1; }
  if(linesymbol == "both")
  { xmg_dataset.[I].linetype   = 2; }

  if (xmg_dataset.[I].linetype == 0 || xmg_dataset.[I].linetype == 2)
  {
    if ( setstyle.nc >= 1 )
    { xmg_dataset.[I].lineshape = setstyle[1;1]; }

    if ( setstyle.nc >= 2 )
    { xmg_dataset.[I].linecolor = setstyle[1;2]; }

    if ( setstyle.nc >= 3 )
    { xmg_dataset.[I].linesize  = setstyle[1;3]; }
  }

  if (xmg_dataset.[I].linetype == 1)
  {
    if ( setstyle.nc >= 1 )
    { xmg_dataset.[I].symbshape = setstyle[1;1]; }

    if ( setstyle.nc >= 2 )
    {
      xmg_dataset.[I].symbcolor = setstyle[1;2];
      xmg_dataset.[I].linecolor = setstyle[1;2];
    }

    if ( setstyle.nc >= 3 )
    { xmg_dataset.[I].symbsize  = setstyle[1;3]; }

    if ( setstyle.nc >= 4 )
    { xmg_dataset.[I].symblwid  = setstyle[1;4]; }
  }

  if (xmg_dataset.[I].linetype == 2)
  {
    if ( setstyle.nc >= 1 && setstyle.nr >=2 )
    {
      xmg_dataset.[I].symbshape = setstyle[2;1];
    }
    if ( setstyle.nc >= 2 && setstyle.nr >=2 )
    {
      xmg_dataset.[I].symbcolor = setstyle[2;2];
    }
    if ( setstyle.nc >= 3 && setstyle.nr >=2 )
    {
      xmg_dataset.[I].symbsize  = setstyle[2;3];
    }
    if ( setstyle.nc >= 4 && setstyle.nr >=2 )
    {
      xmg_dataset.[I].symblwid  = setstyle[2;4];
    }
  }
  return 1;
};


//
// functions associating datasets with graphs
//

xmggraph = function( datasets, N){
    if (!exist(N)){N=DefXMGgraph;}
    defxmggraph(N);
    if(type(datasets)!="real")   { return 0; }
    if(datasets.nr*datasets.nc<1){ return 0; }
    xmg_graph.[N].datasets = datasets;
    if(!exist(xmg_graph.[N].xmin)){
        xmg_graph.[N].xmin = 10^150;
        for(I in datasets){
            if(!exist(xmg_dataset.[I])){ continue; }
            xmg_graph.[N].xmin = min( [ xmg_graph.[N].xmin, xmg_dataset.[I].xmin] );
            }
        }
    if(!exist(xmg_graph.[N].xmax)){
        xmg_graph.[N].xmax = -10^150;
        for(I in datasets){
            if(!exist(xmg_dataset.[I])){ continue; }
            xmg_graph.[N].xmax = max( [ xmg_graph.[N].xmax, xmg_dataset.[I].xmax] );
            }
        }
    if(!exist(xmg_graph.[N].ymin)){
        xmg_graph.[N].ymin = 10^150;
        for(I in datasets){
            if(!exist(xmg_dataset.[I])){ continue; }
            xmg_graph.[N].ymin = min( [ xmg_graph.[N].ymin, xmg_dataset.[I].ymin] );
            }
        }
    if(!exist(xmg_graph.[N].ymax)){
        xmg_graph.[N].ymax = -10^150;
        for(I in datasets){
            if(!exist(xmg_dataset.[I])){ continue; }
            xmg_graph.[N].ymax = max( [ xmg_graph.[N].ymax, xmg_dataset.[I].ymax] );
            }
        }
    return 1;
    };


//
// functions associating graphs with papers (pages)
//

xmgpage = function( graphs, P ){
    if (!exist(P)){P=DefXMGpaper;}
    defxmgpaper(P);
    if(type(graphs)!="real") { return 0; }
    if(graphs.nr*graphs.nc<1){ return 0; }
    if(!exist(xmg_paper.[P].papersize))  { xmg_paper.[P].papersize   = "letter"; }
    if(!exist(xmg_paper.[P].orientation)){ xmg_paper.[P].orientation = "landscape"; }
    if(!exist(xmg_paper.[P].filename))   {
        xmg_paper.[P].filename    = "./tmpgraph.gr";
        }
    newg = [];
    for(N in graphs){
        if(exist(xmg_graph.[N])){ newg = [newg, N]; }
        }
    xmg_paper.[P].graphs = newg;
    return 1;
    };


//
// print a page with graphs to a file
//

xmgprint = function ( P ){
    if(class(P)=="list"){
        // P is a list. Assume a list of [x,y] datasets.
        j = 0;
        for(i in members(P)){
            xmgdataset( P.[i], "xy", "", j);
            j++;
            }
        xmggraph([0:j-1],0);
        P = 0;
        }
    if(class(P)=="num"){
        if(P.nr > 1){
            // Assume P is a matrix of datasets in format [x,y1,y2, ..]
            for(i in [2:P.nc]){ xmgdataset( P[;1,i], "xy", "", i-2); }
            xmggraph([0:P.nc-2],0);
            P = 0;
            }
        }
    // none of the above. P is an integer representing a page to which the
    // graphs have been previously defined using xmgdataset and xmggraph commands.
    if(!exist(P)){P=DefXMGpaper;}
    if(P.nr*P.nc!=1){ return 0; }
    if(!exist(xmg_paper.[P])){ return 0; }
    xmgfile=xmg_paper.[P].filename;
    open(xmgfile,"w");
    //
    // printing the header
    //
    fprintf(xmgfile,"%s\n","# Grace project file by RLaB2 using libxmgrace toolkit.");
    fprintf(xmgfile,"%s\n","# part of the rlabplus project, available from sourceforge.net");
    fprintf(xmgfile,"%s\n","# all bugs should be reported to kostrun@krampus.phys.uconn.edu");
    fprintf(xmgfile,"%s\n","@version 50102");
    if(xmg_paper.[P].papersize == "letter")
    {
      if(xmg_paper.[P].orientation == "landscape")
      {
        fprintf(xmgfile,"%s\n","@page size 792, 612");
      }
      if(xmg_paper.[P].orientation == "portrait")
      {
        fprintf(xmgfile,"%s\n","@page size 612, 792");
      }
    }
    if((xmg_paper.[P].papersize == "A4"))
    {
      if(xmg_paper.[P].orientation == "landscape")
      {
        fprintf(xmgfile,"%s\n","@page size 842, 595");
      }
      if(xmg_paper.[P].orientation == "portrait")
      {
        fprintf(xmgfile,"%s\n","@page size 595, 842");
      }
    }
    fprintf(xmgfile,"%s\n","@page scroll 5%");
    fprintf(xmgfile,"%s\n","@page inout 5%");
    fprintf(xmgfile,"%s\n","@link page off");
    fprintf(xmgfile,"%s\n","@map font 0 to \"Times-Roman\", \"Times-Roman\" ");
    fprintf(xmgfile,"%s\n","@map font 1 to \"Times-Italic\", \"Times-Italic\" ");
    fprintf(xmgfile,"%s\n","@map font 2 to \"Times-Bold\", \"Times-Bold\" ");
    fprintf(xmgfile,"%s\n","@map font 3 to \"Times-BoldItalic\", \"Times-BoldItalic\" ");
    fprintf(xmgfile,"%s\n","@map font 4 to \"Helvetica\", \"Helvetica\" ");
    fprintf(xmgfile,"%s\n","@map font 5 to \"Helvetica-Oblique\", \"Helvetica-Oblique\" ");
    fprintf(xmgfile,"%s\n","@map font 6 to \"Helvetica-Bold\", \"Helvetica-Bold\" ");
    fprintf(xmgfile,"%s\n","@map font 7 to \"Helvetica-BoldOblique\", \"Helvetica-BoldOblique\" ");
    fprintf(xmgfile,"%s\n","@map font 8 to \"Courier\", \"Courier\" ");
    fprintf(xmgfile,"%s\n","@map font 9 to \"Courier-Oblique\", \"Courier-Oblique\" ");
    fprintf(xmgfile,"%s\n","@map font 10 to \"Courier-Bold\", \"Courier-Bold\" ");
    fprintf(xmgfile,"%s\n","@map font 11 to \"Courier-BoldOblique\", \"Courier-BoldOblique\" ");
    fprintf(xmgfile,"%s\n","@map font 12 to \"Symbol\", \"Symbol\" ");
    fprintf(xmgfile,"%s\n","@map font 13 to \"ZapfDingbats\", \"ZapfDingbats\" ");
    for(i in 0:128)
    {
      if(!exist(_xmgcolormap.[i])){ continue; }
      fprintf(xmgfile,"@map color %i to (", i);
      fprintf(xmgfile,"%i,",  _xmgcolormap.[i].rgb[1] );
      fprintf(xmgfile,"%i,",  _xmgcolormap.[i].rgb[2] );
      fprintf(xmgfile,"%i), ",_xmgcolormap.[i].rgb[3] );
      fprintf(xmgfile,"\"%s\"\n", _xmgcolormap.[i].color );
    }
    fprintf(xmgfile,"%s\n","@reference date 0 ");
    fprintf(xmgfile,"%s\n","@date wrap off ");
    fprintf(xmgfile,"%s\n","@date wrap year 1950 ");
    fprintf(xmgfile,"%s\n","@default linewidth 1.5 ");
    fprintf(xmgfile,"%s\n","@default linestyle 1 ");
    fprintf(xmgfile,"%s\n","@default color 1 ");
    fprintf(xmgfile,"%s\n","@default pattern 1 ");
    fprintf(xmgfile,"%s\n","@default font 0 ");
    fprintf(xmgfile,"%s\n","@default char size 1.25000000 ");
    fprintf(xmgfile,"%s\n","@default symbol size 1.000000 ");
    fprintf(xmgfile,"%s\n","@default sformat \"%16.8g\" ");
    fprintf(xmgfile,"%s\n","@background color 0 ");
    fprintf(xmgfile,"%s\n","@page background fill on ");
    fprintf(xmgfile,"%s\n","@timestamp off ");
    fprintf(xmgfile,"%s\n","@timestamp 0.03, 0.03 ");
    fprintf(xmgfile,"%s\n","@timestamp color 1 ");
    fprintf(xmgfile,"%s\n","@timestamp rot 0 ");
    fprintf(xmgfile,"%s\n","@timestamp font 0 ");
    fprintf(xmgfile,"%s\n","@timestamp char size 1.000000 ");
    fprintf(xmgfile,"@timestamp def \"%s\"\n", time2dstr(seconds()));
    if(exist(_xmgobjects))
    {
      for( i in _xmgobjectypes )
      {
        if ( exist(_xmgobjects.[i]) )
        {
          if (i == "string")
          {
            for(j in 1:_xmgobjects.[i].cnt)
            {
              fprintf(xmgfile,"@with string\n");
              fprintf(xmgfile,"@    string on\n");
              if(_xmgobjects.[i].[j-1].loctype == "world")
              {
                fprintf(xmgfile,"@    string loctype world\n");
                fprintf(xmgfile,"@    string %s\n", _xmgobjects.[i].[j-1].graph);
              } else {
                fprintf(xmgfile,"@    string loctype view\n");
              }
              fprintf(xmgfile,"@    string %g, ", _xmgobjects.[i].[j-1].location[1]);
              fprintf(xmgfile,"%g\n", _xmgobjects.[i].[j-1].location[2]);
              fprintf(xmgfile,"@    string color %g\n", _xmgobjects.[i].[j-1].color);
              fprintf(xmgfile,"@    string rot %g\n", _xmgobjects.[i].[j-1].rot);
              fprintf(xmgfile,"@    string font %g\n", _xmgobjects.[i].[j-1].font);
              fprintf(xmgfile,"@    string just %g\n", _xmgobjects.[i].[j-1].just);
              fprintf(xmgfile,"@    string char size %g\n", _xmgobjects.[i].[j-1].size);
              fprintf(xmgfile,"@    string def %s\n",  _xmgobjects.[i].[j-1].string);
            }
          }
          if (i == "line")
          {
            for(j in 1:_xmgobjects.[i].cnt)
            {
              fprintf(xmgfile,"@with line\n");
              fprintf(xmgfile,"@    line on\n");
              if(_xmgobjects.[i].[j-1].loctype == "world")
              {
                fprintf(xmgfile,"@    line loctype world\n");
                fprintf(xmgfile,"@    line %s\n", _xmgobjects.[i].[j-1].graph);
              } else {
                fprintf(xmgfile,"@    line loctype view\n");
              }
              fprintf(xmgfile,"@    line ");
              fprintf(xmgfile,"%g, ", _xmgobjects.[i].[j-1].location[1]);
              fprintf(xmgfile,"%g, ", _xmgobjects.[i].[j-1].location[2]);
              fprintf(xmgfile,"%g, ", _xmgobjects.[i].[j-1].location[3]);
              fprintf(xmgfile,"%g\n", _xmgobjects.[i].[j-1].location[4]);
              fprintf(xmgfile,"@    line linewidth %g\n", _xmgobjects.[i].[j-1].width);
              fprintf(xmgfile,"@    line linestyle %g\n", _xmgobjects.[i].[j-1].style);
              fprintf(xmgfile,"@    line color %g\n", _xmgobjects.[i].[j-1].color);
              fprintf(xmgfile,"@    line arrow %g\n", _xmgobjects.[i].[j-1].arrow);
              fprintf(xmgfile,"@    line arrow length %g\n", _xmgobjects.[i].[j-1].alen);
              fprintf(xmgfile,"@    line arrow layout ");
              fprintf(xmgfile,"%g, ", _xmgobjects.[i].[j-1].alayout[1]);
              fprintf(xmgfile,"%g\n", _xmgobjects.[i].[j-1].alayout[2]);
              fprintf(xmgfile,"@line def\n");
            }
          }
        }
      }
      xmgobjectsclear();
    }
    noffset = min(xmg_paper.[P].graphs);
    for(n in xmg_paper.[P].graphs)
    {
      fprintf(xmgfile,"@g%1i on\n", n - noffset);
      fprintf(xmgfile,"@g%1i hidden false\n", n - noffset);
      fprintf(xmgfile,"@g%1i type XY\n", n - noffset);
      fprintf(xmgfile,"@g%1i stacked false\n", n - noffset);
      fprintf(xmgfile,"@g%1i bar hgap 0.000000\n", n - noffset);
      fprintf(xmgfile,"@g%1i fixedpoint off\n", n - noffset);
      fprintf(xmgfile,"@g%1i fixedpoint type 0\n", n - noffset);
      fprintf(xmgfile,"@g%1i fixedpoint xy 0.000000, 0.000000\n", n - noffset);
      fprintf(xmgfile,"@g%1i fixedpoint format general general\n", n - noffset);
      fprintf(xmgfile,"@g%1i fixedpoint prec 6, 6\n", n - noffset);
      fprintf(xmgfile,"@with g%1i\n", n - noffset);
      fprintf(xmgfile,"@    world xmin %g\n",xmg_graph.[n].xmin);
      fprintf(xmgfile,"@    world xmax %g\n",xmg_graph.[n].xmax);
      fprintf(xmgfile,"@    world ymin %g\n",xmg_graph.[n].ymin);
      fprintf(xmgfile,"@    world ymax %g\n",xmg_graph.[n].ymax);
      fprintf(xmgfile,"@    stack world 0, 0, 0, 0\n");
      fprintf(xmgfile,"@    znorm 1\n");
      if( xmg_paper.[P].orientation == "landscape")
      {
        fprintf(xmgfile,"@    view xmin %g\n", 0.150000 + xmg_graph.[n].size[1] );
        fprintf(xmgfile,"@    view xmax %g\n", 0.150000 + xmg_graph.[n].size[2] );
        fprintf(xmgfile,"@    view ymin %g\n", 0.150000 + 0.7*xmg_graph.[n].size[3] );
        fprintf(xmgfile,"@    view ymax %g\n", 0.150000 + 0.7*xmg_graph.[n].size[4] );
      } else {
        fprintf(xmgfile,"@    view xmin %g\n", 0.150000 + 0.7*xmg_graph.[n].size[1] );
        fprintf(xmgfile,"@    view xmax %g\n", 0.150000 + 0.7*xmg_graph.[n].size[2] );
        fprintf(xmgfile,"@    view ymin %g\n", 0.150000 + xmg_graph.[n].size[3] );
        fprintf(xmgfile,"@    view ymax %g\n", 0.150000 + xmg_graph.[n].size[4] );
      }
      fprintf(xmgfile,"@    title \"%s\"\n", xmg_graph.[n].title);
      fprintf(xmgfile,"@    title font 0\n");
      fprintf(xmgfile,"@    title size %g\n",xmg_graph.[n].titlesize);
      fprintf(xmgfile,"@    title color 1\n");
      fprintf(xmgfile,"@    subtitle \"%s\"\n", xmg_graph.[n].subtitle);
      fprintf(xmgfile,"@    subtitle font 0\n");
      fprintf(xmgfile,"@    subtitle size %g\n", xmg_graph.[n].subtitlesize);
      fprintf(xmgfile,"@    subtitle color 1\n");
      fprintf(xmgfile,"@    xaxes scale %s\n", xmg_graph.[n].xscale);
      fprintf(xmgfile,"@    yaxes scale %s\n", xmg_graph.[n].yscale);
      fprintf(xmgfile,"@    xaxes invert off\n");
      fprintf(xmgfile,"@    yaxes invert off\n");
      fprintf(xmgfile,"@    xaxis  on\n");
      fprintf(xmgfile,"@    xaxis  type zero false\n");
      fprintf(xmgfile,"@    xaxis  offset 0.000000 , 0.000000\n");
      fprintf(xmgfile,"@    xaxis  bar on\n");
      fprintf(xmgfile,"@    xaxis  bar color 1\n");
      fprintf(xmgfile,"@    xaxis  bar linestyle 1\n");
      fprintf(xmgfile,"@    xaxis  bar linewidth 1.0\n");
      fprintf(xmgfile,"@    xaxis  label \"%s\"\n", xmg_graph.[n].xlabel);
      fprintf(xmgfile,"@    xaxis  label layout para\n");
      fprintf(xmgfile,"@    xaxis  label place auto\n");
      fprintf(xmgfile,"@    xaxis  label char size %g\n", xmg_graph.[n].xlabelsize);
      fprintf(xmgfile,"@    xaxis  label font 0\n");
      fprintf(xmgfile,"@    xaxis  label color 1\n");
      fprintf(xmgfile,"@    xaxis  label place normal\n");

      if (sum(xmg_graph.[n].xticks)==0)
      {
        fprintf(xmgfile,"@    xaxis  tick off\n");
        fprintf(xmgfile,"@    xaxis  ticklabel off\n");
      } else {
        fprintf(xmgfile,"@    xaxis  tick on\n");
        fprintf(xmgfile,"@    xaxis  tick major %f\n", xmg_graph.[n].xticks[1]);
        fprintf(xmgfile,"@    xaxis  tick minor ticks %f\n", xmg_graph.[n].xticks[2]);
        fprintf(xmgfile,"@    xaxis  tick default 6\n");
        fprintf(xmgfile,"@    xaxis  tick place rounded true\n");
        fprintf(xmgfile,"@    xaxis  tick in\n");
        fprintf(xmgfile,"@    xaxis  tick major size 1.000000\n");
        if(exist(xmg_graph.[n].xgrid))
        {
          if(exist(xmg_graph.[n].xgrid.major))
          {
            fprintf(xmgfile,"@    xaxis  tick major color %i\n", ...
                xmg_graph.[n].xgrid.major.color);
            fprintf(xmgfile,"@    xaxis  tick major linewidth %g\n", ...
                xmg_graph.[n].xgrid.major.width );
            fprintf(xmgfile,"@    xaxis  tick major linestyle %i\n", ...
                xmg_graph.[n].xgrid.major.style );
            fprintf(xmgfile,"@    xaxis  tick major grid on\n");
          }
          if(exist(xmg_graph.[n].xgrid.minor))
          {
            fprintf(xmgfile,"@    xaxis  tick minor color %i\n", ...
                xmg_graph.[n].xgrid.minor.color);
            fprintf(xmgfile,"@    xaxis  tick minor linewidth %g\n", ...
                xmg_graph.[n].xgrid.minor.width);
            fprintf(xmgfile,"@    xaxis  tick minor linestyle %i\n", ...
                xmg_graph.[n].xgrid.minor.style);
            fprintf(xmgfile,"@    xaxis  tick minor grid on\n");
          } else {
            fprintf(xmgfile,"@    xaxis  tick minor color 1\n");
            fprintf(xmgfile,"@    xaxis  tick minor linewidth 1.0\n");
            fprintf(xmgfile,"@    xaxis  tick minor linestyle 1\n");
            fprintf(xmgfile,"@    xaxis  tick minor grid off\n");
          }
        } else {
          fprintf(xmgfile,"@    xaxis  tick major color %i\n", 1 );
          fprintf(xmgfile,"@    xaxis  tick major linewidth %g\n", 1.0 );
          fprintf(xmgfile,"@    xaxis  tick major linestyle %i\n", 1 );
          fprintf(xmgfile,"@    xaxis  tick major grid off\n");
          fprintf(xmgfile,"@    xaxis  tick minor color 1\n");
          fprintf(xmgfile,"@    xaxis  tick minor linewidth 1.0\n");
          fprintf(xmgfile,"@    xaxis  tick minor linestyle 1\n");
          fprintf(xmgfile,"@    xaxis  tick minor grid off\n");
        } // if(exist(xmg_graph.[n].xgrid))

        fprintf(xmgfile,"@    xaxis  tick minor size 0.500000\n");
        fprintf(xmgfile,"@    xaxis  ticklabel on\n");
        fprintf(xmgfile,"@    xaxis  ticklabel format general\n");
        fprintf(xmgfile,"@    xaxis  ticklabel prec 5\n");
        fprintf(xmgfile,"@    xaxis  ticklabel formula \"\"\n");
        fprintf(xmgfile,"@    xaxis  ticklabel append \"\"\n");
        fprintf(xmgfile,"@    xaxis  ticklabel prepend \"\"\n");
        fprintf(xmgfile,"@    xaxis  ticklabel angle 0\n");
        fprintf(xmgfile,"@    xaxis  ticklabel skip 0\n");
        fprintf(xmgfile,"@    xaxis  ticklabel stagger 0\n");
        fprintf(xmgfile,"@    xaxis  ticklabel place normal\n");
        fprintf(xmgfile,"@    xaxis  ticklabel offset auto\n");
        fprintf(xmgfile,"@    xaxis  ticklabel offset 0.000000 , 0.010000\n");
        fprintf(xmgfile,"@    xaxis  ticklabel start type auto\n");
        fprintf(xmgfile,"@    xaxis  ticklabel start 0.000000\n");
        fprintf(xmgfile,"@    xaxis  ticklabel stop type auto\n");
        fprintf(xmgfile,"@    xaxis  ticklabel stop 0.000000\n");
        fprintf(xmgfile,"@    xaxis  ticklabel char size %g\n", xmg_graph.[n].xticksize);
        fprintf(xmgfile,"@    xaxis  ticklabel font 0\n");
        fprintf(xmgfile,"@    xaxis  ticklabel color 1\n");
        fprintf(xmgfile,"@    xaxis  tick place both\n");
        if(isempty(xmg_graph.[n].xtickspec))
        {
          fprintf(xmgfile,"@    xaxis  tick spec type none\n");
        } else {
          if (xmg_graph.[n].xtickspec.kind == "ticks")
          {
            fprintf(xmgfile,"@    xaxis  tick spec type ticks\n");
            fprintf(xmgfile,"@    xaxis  tick spec %g\n", ...
                length(xmg_graph.[n].xtickspec.tick) );
            for (_i in 1:length(xmg_graph.[n].xtickspec.tick))
            {
              fprintf(xmgfile,"@    xaxis  tick major %g, %g\n", ...
                  _i-1, xmg_graph.[n].xtickspec.tick[_i] );
            }
          }
          if (xmg_graph.[n].xtickspec.kind == "both")
          {
            fprintf(xmgfile,"@    xaxis  tick spec type both\n");
            fprintf(xmgfile,"@    xaxis  tick spec %g\n", ...
                length(xmg_graph.[n].xtickspec.tick) );
            for (_i in 1:length(xmg_graph.[n].xtickspec.tick))
            {
              fprintf(xmgfile,"@    xaxis  tick major %g, %g\n", ...
                  _i-1, xmg_graph.[n].xtickspec.tick[_i] );
              fprintf(xmgfile,"@    xaxis  ticklabel %g, ", _i-1 );
              if (type(xmg_graph.[n].xtickspec.label) == "string")
              {
                fprintf(xmgfile,"\"%s\"\n", xmg_graph.[n].xtickspec.label[_i] );
              } else {
                fprintf(xmgfile,"\"%g\"\n", xmg_graph.[n].xtickspec.label[_i] );
              }
            }
          }
        } // if(isempty(xmg_graph.[n].xtickspec)
      }
      fprintf(xmgfile,"@    yaxis  on\n");
      fprintf(xmgfile,"@    yaxis  type zero false\n");
      fprintf(xmgfile,"@    yaxis  offset 0.000000 , 0.000000\n");
      fprintf(xmgfile,"@    yaxis  bar on\n");
      fprintf(xmgfile,"@    yaxis  bar color 1\n");
      fprintf(xmgfile,"@    yaxis  bar linestyle 1\n");
      fprintf(xmgfile,"@    yaxis  bar linewidth 1.0\n");
      fprintf(xmgfile,"@    yaxis  label \"%s\"\n", xmg_graph.[n].ylabel);
      fprintf(xmgfile,"@    yaxis  label layout para\n");
      fprintf(xmgfile,"@    yaxis  label place spec\n");
      fprintf(xmgfile,"@    yaxis  label place 0.000000, %g\n", ...
          0.08 + 0.02 * (xmg_graph.[n].ylabelsize - 1) );
      fprintf(xmgfile,"@    yaxis  label char size %g\n", xmg_graph.[n].ylabelsize);
      fprintf(xmgfile,"@    yaxis  label font 0\n");
      fprintf(xmgfile,"@    yaxis  label color 1\n");
      fprintf(xmgfile,"@    yaxis  label place normal\n");
      if (sum(xmg_graph.[n].yticks)==0)
      {
        fprintf(xmgfile,"@    yaxis  tick off\n");
        fprintf(xmgfile,"@    yaxis  ticklabel off\n");
      } else {
        fprintf(xmgfile,"@    yaxis  tick on\n");
        fprintf(xmgfile,"@    yaxis  tick major %f\n", xmg_graph.[n].yticks[1]);
        fprintf(xmgfile,"@    yaxis  tick minor ticks %f\n", xmg_graph.[n].yticks[2]);
        fprintf(xmgfile,"@    yaxis  tick default 6\n");
        fprintf(xmgfile,"@    yaxis  tick place rounded true\n");
        fprintf(xmgfile,"@    yaxis  tick in\n");
        fprintf(xmgfile,"@    yaxis  tick major size 1.000000\n");
        if(exist(xmg_graph.[n].ygrid))
        {
          if(exist(xmg_graph.[n].ygrid.major))
          {
            fprintf(xmgfile,"@    yaxis  tick major color %i\n", ...
                xmg_graph.[n].ygrid.major.color);
            fprintf(xmgfile,"@    yaxis  tick major linewidth %g\n", ...
                xmg_graph.[n].ygrid.major.width );
            fprintf(xmgfile,"@    yaxis  tick major linestyle %i\n", ...
                xmg_graph.[n].ygrid.major.style );
            fprintf(xmgfile,"@    yaxis  tick major grid on\n");
          }
          if(exist(xmg_graph.[n].ygrid.minor))
          {
            fprintf(xmgfile,"@    yaxis  tick minor color %i\n", ...
                xmg_graph.[n].ygrid.minor.color);
            fprintf(xmgfile,"@    yaxis  tick minor linewidth %g\n", ...
                xmg_graph.[n].ygrid.minor.width);
            fprintf(xmgfile,"@    yaxis  tick minor linestyle %i\n", ...
                xmg_graph.[n].ygrid.minor.style);
            fprintf(xmgfile,"@    yaxis  tick minor grid on\n");
          } else {
            fprintf(xmgfile,"@    yaxis  tick minor color 1\n");
            fprintf(xmgfile,"@    yaxis  tick minor linewidth 1.0\n");
            fprintf(xmgfile,"@    yaxis  tick minor linestyle 1\n");
            fprintf(xmgfile,"@    yaxis  tick minor grid off\n");
          }
        } else {
          fprintf(xmgfile,"@    yaxis  tick major color %i\n", 1 );
          fprintf(xmgfile,"@    yaxis  tick major linewidth %g\n", 1.0 );
          fprintf(xmgfile,"@    yaxis  tick major linestyle %i\n", 1 );
          fprintf(xmgfile,"@    yaxis  tick major grid off\n");
          fprintf(xmgfile,"@    yaxis  tick minor color 1\n");
          fprintf(xmgfile,"@    yaxis  tick minor linewidth 1.0\n");
          fprintf(xmgfile,"@    yaxis  tick minor linestyle 1\n");
          fprintf(xmgfile,"@    yaxis  tick minor grid off\n");
        }
        fprintf(xmgfile,"@    yaxis  tick minor size 0.500000\n");
        fprintf(xmgfile,"@    yaxis  ticklabel on\n");
        fprintf(xmgfile,"@    yaxis  ticklabel format general\n");
        fprintf(xmgfile,"@    yaxis  ticklabel prec 5\n");
        fprintf(xmgfile,"@    yaxis  ticklabel formula \"\"\n");
        fprintf(xmgfile,"@    yaxis  ticklabel append \"\"\n");
        fprintf(xmgfile,"@    yaxis  ticklabel prepend \"\"\n");
        fprintf(xmgfile,"@    yaxis  ticklabel angle 0\n");
        fprintf(xmgfile,"@    yaxis  ticklabel skip 0\n");
        fprintf(xmgfile,"@    yaxis  ticklabel stagger 0\n");
        fprintf(xmgfile,"@    yaxis  ticklabel place normal\n");
        fprintf(xmgfile,"@    yaxis  ticklabel offset auto\n");
        fprintf(xmgfile,"@    yaxis  ticklabel offset 0.000000 , 0.010000\n");
        fprintf(xmgfile,"@    yaxis  ticklabel start type auto\n");
        fprintf(xmgfile,"@    yaxis  ticklabel start 0.000000\n");
        fprintf(xmgfile,"@    yaxis  ticklabel stop type auto\n");
        fprintf(xmgfile,"@    yaxis  ticklabel stop 0.000000\n");
        fprintf(xmgfile,"@    yaxis  ticklabel char size %g\n", xmg_graph.[n].yticksize);
        fprintf(xmgfile,"@    yaxis  ticklabel font 0\n");
        fprintf(xmgfile,"@    yaxis  ticklabel color 1\n");
        fprintf(xmgfile,"@    yaxis  tick place both\n");
        fprintf(xmgfile,"@    yaxis  tick spec type none\n");
      }
      fprintf(xmgfile,"@    altxaxis  off\n");
      fprintf(xmgfile,"@    altyaxis  off\n");
      if (size(xmg_graph.[n].legend)==0)
      {
        fprintf(xmgfile,"@    legend off\n");
      } else {
        fprintf(xmgfile,"@    legend on\n");
        fprintf(xmgfile,"@    legend loctype view\n");
        fprintf(xmgfile,"@    legend ");
        fprintf(xmgfile,"%f,%f\n", xmg_graph.[n].legend.x, xmg_graph.[n].legend.y);
        fprintf(xmgfile,"@    legend box color 1\n");
        fprintf(xmgfile,"@    legend box pattern 1\n");
        fprintf(xmgfile,"@    legend box linewidth 1.0\n");
        fprintf(xmgfile,"@    legend box linestyle 1\n");
        fprintf(xmgfile,"@    legend box fill color 0\n");
        fprintf(xmgfile,"@    legend box fill pattern 1\n");
        fprintf(xmgfile,"@    legend font 0\n");
        fprintf(xmgfile,"@    legend char size %g\n", xmg_graph.[n].legend.size);
        fprintf(xmgfile,"@    legend color %g\n", xmg_graph.[n].legend.color);
        fprintf(xmgfile,"@    legend length %g\n", xmg_graph.[n].legend.length);
        fprintf(xmgfile,"@    legend vgap %g\n", xmg_graph.[n].legend.vgap);
        fprintf(xmgfile,"@    legend hgap %g\n", xmg_graph.[n].legend.hgap);
        fprintf(xmgfile,"@    legend invert false\n");
      }
      fprintf(xmgfile,"@    frame type 0\n");
      fprintf(xmgfile,"@    frame linestyle 1\n");
      fprintf(xmgfile,"@    frame linewidth 2.0\n");
      fprintf(xmgfile,"@    frame color 1\n");
      fprintf(xmgfile,"@    frame pattern 1\n");
      fprintf(xmgfile,"@    frame background color 0\n");
      fprintf(xmgfile,"@    frame background pattern 0\n");
      for (i in xmg_graph.[n].datasets)
      {
        if( !exist(xmg_dataset.[i]) )
        { continue; }
        if(isempty(xmg_dataset.[i]))
        {
          printf("xmgprint: Dataset %i is empty! Not trying to print it.\n", i);
          continue;
        }
        fprintf(xmgfile,"@    s%1i hidden false\n", i );
        fprintf(xmgfile,"@    s%1i type %s\n",i, xmg_dataset.[i].plottype);
        if( xmg_dataset.[i].linetype == 1 || xmg_dataset.[i].linetype == 2)
        {
          fprintf(xmgfile,"@    s%1i symbol %i\n",i, xmg_dataset.[i].symbshape);
          fprintf(xmgfile,"@    s%1i symbol size %g\n" , i, xmg_dataset.[i].symbsize );
          fprintf(xmgfile,"@    s%1i symbol color %i\n", i, xmg_dataset.[i].symbcolor);
          fprintf(xmgfile,"@    s%1i symbol pattern 1\n", i );
          fprintf(xmgfile,"@    s%1i symbol fill color 1\n",i);
          fprintf(xmgfile,"@    s%1i symbol fill pattern 0\n",i);
          fprintf(xmgfile,"@    s%1i symbol linewidth %g\n",i,xmg_dataset.[i].symblwid);
          fprintf(xmgfile,"@    s%1i symbol linestyle 1\n",i);
          fprintf(xmgfile,"@    s%1i symbol char 65\n",i);
          fprintf(xmgfile,"@    s%1i symbol char font 0\n",i);
          fprintf(xmgfile,"@    s%1i symbol skip 0\n",i);
        } else {
          fprintf(xmgfile,"@    s%1i symbol 0\n", i);
        }
        if( xmg_dataset.[i].linetype == 0 || xmg_dataset.[i].linetype == 2)
        {
          fprintf(xmgfile,"@    s%1i line type 1\n",i );
          fprintf(xmgfile,"@    s%1i line linestyle %1i\n",i,xmg_dataset.[i].lineshape );
          fprintf(xmgfile,"@    s%1i line linewidth %g\n", i,xmg_dataset.[i].linesize );
          fprintf(xmgfile,"@    s%1i line color %1i\n", i, xmg_dataset.[i].linecolor );
          fprintf(xmgfile,"@    s%1i line pattern 1\n",i);
        } else {
          fprintf(xmgfile,"@    s%1i line type 0\n",i);
        }
        fprintf(xmgfile,"@    s%1i baseline type 0\n",i);
        fprintf(xmgfile,"@    s%1i baseline off\n",i);
        fprintf(xmgfile,"@    s%1i dropline off\n",i);
        fprintf(xmgfile,"@    s%1i fill type 0\n",i);
        fprintf(xmgfile,"@    s%1i fill rule 0\n",i);
        fprintf(xmgfile,"@    s%1i fill color 1\n",i);
        fprintf(xmgfile,"@    s%1i fill pattern 1\n",i);
        fprintf(xmgfile,"@    s%1i avalue off\n",i);
        fprintf(xmgfile,"@    s%1i avalue type 2\n",i);
        fprintf(xmgfile,"@    s%1i avalue char size 1.000000\n",i);
        fprintf(xmgfile,"@    s%1i avalue font 0\n",i);
        fprintf(xmgfile,"@    s%1i avalue color 1\n",i);
        fprintf(xmgfile,"@    s%1i avalue rot 0\n",i);
        fprintf(xmgfile,"@    s%1i avalue format general\n",i);
        fprintf(xmgfile,"@    s%1i avalue prec 3\n",i);
        fprintf(xmgfile,"@    s%1i avalue prepend \"\"\n",i);
        fprintf(xmgfile,"@    s%1i avalue append \"\"\n",i);
        fprintf(xmgfile,"@    s%1i avalue offset 0.000000 , 0.000000\n",i);
        fprintf(xmgfile,"@    s%1i errorbar on\n",i);
        fprintf(xmgfile,"@    s%1i errorbar place both\n",i);
        fprintf(xmgfile,"@    s%1i errorbar color %1i\n",i, xmg_dataset.[i].linecolor);
        fprintf(xmgfile,"@    s%1i errorbar pattern 1\n",i);
        fprintf(xmgfile,"@    s%1i errorbar size 1.000000\n",i);
        fprintf(xmgfile,"@    s%1i errorbar linewidth 1.0\n",i);
        fprintf(xmgfile,"@    s%1i errorbar linestyle 1\n",i);
        fprintf(xmgfile,"@    s%1i errorbar riser linewidth 1.0\n",i);
        fprintf(xmgfile,"@    s%1i errorbar riser linestyle 1\n",i);
        fprintf(xmgfile,"@    s%1i errorbar riser clip off\n",i);
        fprintf(xmgfile,"@    s%1i errorbar riser clip length 0.100000\n",i);
        fprintf(xmgfile,"@    s%1i comment \"RLaB xmgr output\"\n",i);
        fprintf(xmgfile,"@    s%1i legend  ",i);
        if (strlen(xmg_dataset.[i].dataname) > 0)
        {
          fprintf(xmgfile,"%s\n", xmg_dataset.[i].dataname);
        } else {
          fprintf(xmgfile,"\"\"\n");
        }
      }
    }
    for(n in xmg_paper.[P].graphs)
    {
      for (i in xmg_graph.[n].datasets)
      {
        if(!exist(xmg_dataset.[i]))
        { continue; }
        if(xmg_dataset.[i].dataset.nr*xmg_dataset.[i].dataset.nc == 0)
        {
          printf("xmgprint: dataset %i is empty. Skipping it!\n", i);
          continue;
        }
        fprintf(xmgfile,"@target G%1i.S%1i\n", n - noffset, i);
        fprintf(xmgfile,"@type %s\n",xmg_dataset.[i].plottype);
        writem(xmgfile,xmg_dataset.[i].dataset);
      }
      fprintf(xmgfile,"&\n");
    }
    close(xmgfile);
    printf("The xmgrace file %s sucessfuly created.\n", xmg_paper.[P].filename);
    return 1;
};

