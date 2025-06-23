//  xfig interface for the RLaB2/RLaB3
//    rlabplus (c) 2000-2017 marijan kostrun
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
static(_INIT_);
if (exist(_INIT_))
{
  EOF
}

//
//
// L I B X F I G :
//    global variables and functions for the library
//
//

//
// clip world data set to world bbox
//
static(__intersection_line_bbox, __intersection_line_line, __clip2world);

__intersection_line_line = function (ls1, ls2)
{
  // ls1,ls2 = [x1, y1; x2, y2]
  // check arguments
  if (class(ls1)!="num")
  { return []; }
  if (ls1.nr!=2 || ls1.nc!=2)
  { return []; }
  if (class(ls2)!="num")
  { return []; }
  if (ls2.nr!=2 || ls2.nc!=2)
  { return []; }

  m = [(ls1[1;] - ls1[2;])',(ls2[2;] - ls2[1;])'];
  if (det(m)!=0)
  {
    b = (ls2[2;] - ls1[2;])';
    a = solve(m,b);
    if ((min(a)<0)||(max(a)>1))
    { return []; }

    rval = a[1] .* ls1[1;] + (1-a[1]).* ls1[2;];
  }

  // lines are parallel or overlapping: assume no intersection
  return [];
};

__intersection_line_bbox = function (ls1, bbox)
{
  // go over 4 lines comprising the bounding box
  // 1:
  ls2 = bbox;
  ls2[2;2] = ls2[1;2];
  rval = __intersection_line_line(ls1, ls2);
  if (length(rval)==2)
  { return rval; }

  // 2:
  shiftu(ls2);
  ls2[2;] = bbox[2;];
  rval = __intersection_line_line(ls1, ls2);
  if (length(rval)==2)
  { return rval; }

  // 3:
  shiftu(ls2);
  ls2[2;] = [bbox[1;1], bbox[2;2]];
  rval = __intersection_line_line(ls1, ls2);
  if (length(rval)==2)
  { return rval; }

  shiftu(ls2);
  ls2[2;] = bbox[1;];
  rval = __intersection_line_line(ls1, ls2);
  return rval;
};

__clip2world = function(xy, bbox)
{
  if (class(xy)!="num")
  { return []; }
  if (xy.nc != 2)
  { return []; }

  // rval will feature value > 0 for each point that is plotted
  // and 0 for the points that are not plotted
  // each set of points with same value, say, 1
  // will be plotted as one line segment
  // this is interrupted by 0 for points that are not plotted, 
  // then for next value, a new segment is created, and so forth
  rval = ones(1,xy.nr);

  if (class(bbox)!="num")
  { return rval; }
  if (bbox.nr!=2 || bbox.nc!=2)
  { return rval; }

  segment_idx = 1;
  point_idx = 1;
  while (point_idx <= xy.nr)
  {
    if ( (xy[point_idx;1]>=bbox[1;1]) && (xy[point_idx;1]<=bbox[2;1]))
    {
      if ( (xy[point_idx;2]>=bbox[1;2]) && (xy[point_idx;2]<=bbox[2;2]))
      {
        rval[point_idx] = segment_idx;
        point_idx++;
        continue;
      }
    }

    rval[point_idx] = 0; // out of bounds
    if (point_idx > 1)
    {
      if (rval[point_idx-1] != 0)
      { segment_idx++; }
    }

    point_idx++;
    continue;
  }

  return rval;
};

//
// XFIG built-in linestyles
//
static(xfig_line_style);

if (!exist(xfig_line_style))
{
  xfig_line_style = <<>>;
  xfig_line_style.desc_text = ["default", "solid", "dash", "dot", "dash dot", ...
      "dash dot dot","dash dot dot dot"];
  xfig_line_style.val = [-1,0,1,2,3,4,5];
  xfig_line_style.desc_symbol = ["", "_", "-", ".", "-.", "-..", "-..."];
}

//
// XFIG built-in color table
//
static(xfig_color_table_default,xfig_color_table,__find_or_manage_colors_in_color_table);
if (!exist(xfig_color_table_default))
{
  xfig_color_table_default = <<>>;
  xfig_color_table_default.name = [ ...
      "default", ...
      "black", "blue", "green", "cyan", ...
          "red", "magenta", "yellow", "white", ...
              "blue4", "blue3", "blue2", "lightblue", ...
                  "green4", "green3", "green2", "cyan4", ...
                      "cyan3", "cyan2", "red4", "red3", ...
                          "red2", "magenta4", "magenta3", "magenta2", ...
                              "brown4", "brown3", "brown2", "pink4", ...
                                  "pink3", "pink2", "pink", "gold" , ...
      []];
      xfig_color_table_default.id = [-1:31];
}

if (!exist(xfig_color_table))
{
  xfig_color_table= <<>>;
  xfig_color_table.name = blank();
  xfig_color_table.id = [];
  xfig_color_table.rgb = [];
}

_xfig_colors = function()
{
  rval = <<>>;
  rval.name = [xfig_color_table_default.name,xfig_color_table.name];
  rval.id = [xfig_color_table_default.id, xfig_color_table.id];
  rval.rgb  = [zeros(length(xfig_color_table_default.id),3); xfig_color_table.rgb];
  return rval;
};

__find_or_manage_colors_in_color_table = function( n, rgb )
{
  // is 'n' number?
  if (isnumber(n))
  {
    if (n>=-1 && n <= 31)
    { return n; }

    // three possibility for i_n:
    i_n = find(xfig_color_table.id==n);

    if (length(i_n) > 1)
    {
      printf("horrible internal error: help! help!\n");
      printf("refusing to make a choice: reverting to default color!\n");
      return -1;
    }

    if (length(i_n) == 1)
    {
      // assert or replace the color
      if (exist(rgb))
      {
        if (rgb.nc == 3)
        {
          if (any(rgb!=xfig_color_table.rgb[i_n;]))
          { xfig_color_table.rgb[i_n;] = rgb; }
        }
      }
      return xfig_color_table.id[i_n];
    }

    // add color with the number user provided
    if (!exist(rgb))
    { return -1; }

    rval = n;
    xfig_color_table.rgb = [xfig_color_table.rgb; rgb];
    xfig_color_table.name = [xfig_color_table.name; num2str(rval,"user%.0f") ];
    xfig_color_table.id = [xfig_color_table.id, rval];
    return rval;
  }

  // is 'n' string?
  if (strlen(n)<1)
  { return -1; }

  // check if it is among built-in colors
  i_n = find(xfig_color_table_default.name == tolower(n));

  // three possibility for i_n:
  if (length(i_n)==1)
  {
    return xfig_color_table_default.id[i_n];
  }
  if (length(i_n)>1)
  {
    printf("color name '%s' can be one of these colors %s.\n", n, ...
        join(xfig_color_table_default.id[i_n],","));
    printf("refusing to make a choice: reverting to default color!\n");
    return -1;
  }

  // check if it is among additional colors
  i_a = find(xfig_color_table.name == tolower(n));

  // three possibility for i_a:
  if (length(i_a)==1)
  {
    if (exist(rgb))
    {
      if (rgb.nc == 3)
      {
        if (any(rgb!=xfig_color_table.rgb[i_a;]))
        { xfig_color_table.rgb[i_a;] = rgb; }
      }
    }
    return xfig_color_table.id[i_a];
  }
  if (length(i_a)>1)
  {
    printf("color name '%s' can be one of these colors %s.\n", n, ...
        join(xfig_color_table.id[i_a],","));
    printf("refusing to make a choice: reverting to default color!\n");
    return -1;
  }

  // i_a is empty - no such color
  if (!exist(rgb))
  { return -1; }

  rval = 32 + xfig_color_table.name.n;
  xfig_color_table.rgb = [xfig_color_table.rgb; rgb];
  xfig_color_table.name = [xfig_color_table.name; tolower(n)];
  xfig_color_table.id = [xfig_color_table.id, rval];
  return rval;
};


//
// XFIG arrow styles
//
static(xfig_arrow_type);
if (!exist(xfig_arrow_type))
{
  xfig_arrow_type = <<>>;
  xfig_arrow_type.val = [0,1,2,3];
  xfig_arrow_type.desc_text = ["stick", "triangle", "indented", "diamond"];
  xfig_arrow_type.desc_symbol = [">", "|>", ">>", "<>"];
}
static(xfig_arrow_style);
if (!exist(xfig_arrow_style))
{
  xfig_arrow_style = <<>>;
  xfig_arrow_style.val = [0,1];
  xfig_arrow_style.desc_text = ["hollow", "filled"];
}

_xfig_arrow = function (typ,sty,thick,wid,hei)
{
  y = 0;
  s = 0;
  t = 1;
  w = 4;
  h = 8;

  //
  // type of the arrow
  // 
  if (isnumber(typ))
  {
    if (typ>=0 && typ<=3)
    { y = floor(typ); }
  }
  if (strlen(typ)>0)
  {
    i = find(strindex(xfig_arrow_type.desc_text,typ));
    if (length(i)==1)
    { y = xfig_arrow_type.val[i]; }
    if (length(i)==0)
    {
      i = find(strindex(xfig_arrow_type.desc_symbol,typ));
      if (length(i)==1)
      { y = xfig_arrow_type.val[i]; }
    }
  }

  //
  // style of arrow
  //
  if (isnumber(sty))
  {
    if (sty==0 || sty==1)
    { s = sty; }
  }
  if (strlen(sty)>0)
  {
    i = find(strindex(xfig_arrow_style.desc_text,sty));
    if (length(i)==1)
    { s = xfig_arrow_style.val[i]; }
  }

  if (isnumber(thick))
  {
    if (thick >= 0)
    { t = thick; }
  }

  if (isnumber(wid))
  {
    if (wid >= 0)
    { w = wid; }
  }

  if (isnumber(hei))
  {
    if (hei >= 0)
    { h = hei; }
  }

  return [y,s,t,w,h];
};

static (__label, __papersize, __ellipse, __line, __box, __line_style, __picture, __arc);

__picture = function (fn, xy_0, xy_d, orientation, param)
{
  if (strlen(fn)<1)
  { return []; }

  if (!exist(xy_0))
  { return []; }
  if (xy_0.nc != 2)
  { return []; }

  if (!exist(xy_d))
  { return []; }
  if (xy_d.nc != 2)
  { return []; }

  if (!exist(orientation))
  { orientation = 0; }

  sub_type = 5; // imported picture bounding box

  // default parameters
  line_style = -1;
  thickness = 0;
  pen_color = 0;
  fill_color = 0;
  depth = 50;
  pen_style = 0;
  area_fill = -1;
  style_val = 1;
  join_style = 0;
  cap_style = 0;
  radius = 0; // arc_box only
  forward_arrow_flag = 0;
  forward_arrow = nan(1,5);
  backward_arrow_flag = 0;
  backward_arrow = nan(1,5);

  //
  npoints = 5;
  dx = xy_d[1];
  dy = xy_d[2];
  flipped = 0;

  // return value is a list of line-type objects
  rval = <<>>;
  rval.file = fn;

  // compute then flatten bounding box
  dx = abs(dx);
  dy = abs(dy);
  switch(round(orientation,<<bin=90>>))
  {
    case 90:
      bb = xy_0 + [dy,-dx; dy,0; 0,0; 0,-dx; dy,-dx] + [0,dx];
      break;

    case 180:
      bb = xy_0 + [dx,0; 0,0; 0,-dy; dx,-dy; dx,0] + [-dx,dy];
      break;

    case 270:
      bb = xy_0 + [0,0; 0,-dx; dy,-dx; dy,0; 0,0] + [-dy,0];
      break;

    default:
      bb = xy_0 + [0,-dy; dx,-dy; dx,0; 0,0; 0,-dy];
      break;
  }
  data = zeros(1,10);
  for (i in 1:5)
  {
    data[i*2-1] = bb[i;1];
    data[i*2  ] = bb[i;2];
  }

  //
  // process parameters
  //
  if (class(param)!="list")
  {
    rval.desc = [sub_type, line_style, thickness, pen_color, fill_color, depth, ...
        pen_style, area_fill, style_val, join_style, ...
            cap_style, radius, 0, 0, npoints, flipped, data];
    return rval;
  }

  if (exist(param.line_style))
  { line_style = __line_style(param.line_style); }

  if (class(param.thickness)=="num")
  {
    if (param.thickness >= 0)
    { thickness = floor(param.thickness); }
  }

  if (exist(param.pen_color))
  {
    pen_color = __find_or_manage_colors_in_color_table(param.pen_color);
  }

  if (exist(param.fill_color))
  {
    fill_color = __find_or_manage_colors_in_color_table(param.fill_color);
    area_fill = 20;
  }

  if (class(param.depth)=="num")
  {
    if (param.depth >= 1)
    { depth = floor(param.depth); }
  }

  if (class(param.pen_style)=="num")
  {
    if (param.pen_style >= -1)
    { pen_style = floor(param.pen_style); }
  }

  if (class(param.area_fill)=="num")
  {
    if (param.area_fill >= -1)
    { area_fill = floor(param.area_fill); }
  }

  if (class(param.style_val)=="num")
  {
    if (param.style_val > 0)
    { style_val = floor(param.style_val); }
  }

  if (class(param.join_style)=="num")
  {
    if (param.join_style >= 0)
    { join_style = floor(param.join_style); }
  }

  if (class(param.cap_style)=="num")
  {
    if (param.cap_style >= 0)
    { cap_style = floor(param.cap_style); }
  }

  if (class(param.radius)=="num")
  {
    if (param.radius > 0)
    { radius = floor(param.radius); }
  }

  if (class(param.flipped)=="num")
  {
    if (param.flipped > 0)
    { flipped = 1; }
  }

  // let's get out of here
  rval.desc = [sub_type, line_style, thickness, pen_color, fill_color, depth, ...
      pen_style, area_fill, style_val, join_style, ...
          cap_style, radius, 0, 0, npoints, flipped, data];

  return rval;
};



__line_style = function ( val )
{
  if (!exist(val))
  { return -1; }

  if (class(val)=="num")
  {
    if (val>=-1 && val <=5)
    { return floor(val); }
  }

  if (class(val)=="string")
  {
    if (strindex(val, "sol") || strindex(val, "_"))
    {
      return 0;
    }
    idot = length(findstr(val, "dot"));
    if (!idot)
    {
      idot = length(findstr(val, "."));
    }
    idas = length(findstr(val, "das"));
    if (!idas)
    {
      idas = length(findstr(val, "-"));
    }

    if (idot && !idas)
    {
      return 1;  
    }
    if (!idot && idas)
    {
      return 2; 
    }
    if (idot==1 && idas==1)
    {
      return 3;
    }
    if (idot==2 && idas==1)
    {
      return 4;
    }
    if (idot==3 && idas==1)
    {
      return 5;
    }
  }

  return -1;
};


__arc = function(center_xy, xy, sub_type, param)
{
  // sub_type: 1-open, 2-pie wedge

  //
  if (!exist(center_xy))
  { return []; }
  if (center_xy.n != 2)
  { return []; }

  //
  if (!exist(xy))
  { return []; }
  if (xy.nc != 2)
  { return []; }
  if (xy.nr != 3)
  { return []; }

  // default parameters
  line_style = -1;
  thickness = 1;
  pen_color = 0;
  fill_color = 0;
  depth = 50;
  pen_style = 0;
  area_fill = -1;
  style_val = 1;
  cap_style = 0;
  dir = 0;
  forward_arrow_flag = 0;
  backward_arrow_flag = 0;
  xy_data = [xy[1;], xy[2;], xy[3;]];
  forward_arrow = nan(1,5);
  backward_arrow = nan(1,5);

  // return value is a vector

  //
  // process parameters
  //
  if (class(param)!="list")
  {
    rval = [sub_type, line_style, thickness, pen_color, fill_color, depth, ...
        pen_style, area_fill, style_val, ...
            cap_style, dir, forward_arrow_flag, backward_arrow_flag, ...
                center_xy, xy_data, forward_arrow, backward_arrow];
    return rval;
  }

  if (exist(param.line_style))
  { line_style = __line_style(param.line_style); }

  if (class(param.thickness)=="num")
  {
    if (param.thickness >= 0)
    { thickness = floor(param.thickness); }
  }

  if (exist(param.pen_color))
  {
    pen_color = __find_or_manage_colors_in_color_table(param.pen_color);
  }

  if (exist(param.fill_color))
  {
    fill_color = __find_or_manage_colors_in_color_table(param.fill_color);
    area_fill = 20;
  }

  if (class(param.depth)=="num")
  {
    if (param.depth >= 1)
    { depth = floor(param.depth); }
  }

  if (class(param.pen_style)=="num")
  {
    if (param.pen_style >= -1)
    { pen_style = floor(param.pen_style); }
  }

  if (class(param.area_fill)=="num")
  {
    if (param.area_fill >= -1)
    { area_fill = floor(param.area_fill); }
  }

  if (class(param.style_val)=="num")
  {
    if (param.style_val > 0)
    { style_val = floor(param.style_val); }
  }

  if (class(param.cap_style)=="num")
  {
    if (param.cap_style >= 0)
    { cap_style = floor(param.cap_style); }
  }

  if (strlen(param.dir)>0)
  {
    dir = 0;
    if (strindex(param.dir,"count")>0)
    { dir = 1; }
  }

  if (class(param.forward_arrow)=="num")
  {
    if (param.forward_arrow.n == 5)
    {
      forward_arrow = param.forward_arrow;
    }
  }
  if (all(!isnan(forward_arrow)))
  {
    forward_arrow_flag = 1;
  }

  if (class(param.backward_arrow)=="num")
  {
    if (param.backward_arrow.n == 5)
    {
      backward_arrow = param.backward_arrow;
    }
  }
  if (all(!isnan(backward_arrow)))
  {
    backward_arrow_flag = 1;
  }

  // let's get out of here
  rval = [sub_type, line_style, thickness, pen_color, fill_color, depth, ...
      pen_style, area_fill, style_val, ...
          cap_style, dir, forward_arrow_flag, backward_arrow_flag, ...
              center_xy, xy_data, forward_arrow, backward_arrow];
  return rval;
};


__line = function(xy, param)
{
  //
  if (!exist(xy))
  { return []; }
  if (xy.nc != 2)
  { return []; }

  sub_type = 1; // polyline: immutable

  if (all(xy[1;] == xy[xy.nr;]))
  {
    sub_type = 3; // polygon
  }

  // default parameters
  line_style = -1;
  thickness = 1;
  pen_color = 0;
  fill_color = 0;
  depth = 50;
  pen_style = 0;
  area_fill = -1;
  style_val = 1;
  join_style = 0;
  cap_style = 0;
  radius = 0; // arc_box only
  forward_arrow_flag = 0;
  forward_arrow = nan(1,5);
  backward_arrow_flag = 0;
  backward_arrow = nan(1,5);

  //
  npoints = xy.nr;

  // return value is a list of line-type objects
  rval = <<>>;
  rval.data = xy;

  //
  // process parameters
  //
  if (class(param)!="list")
  {
    rval.desc = [sub_type, line_style, thickness, pen_color, fill_color, depth, ...
        pen_style, area_fill, style_val, join_style, ...
            cap_style, radius, forward_arrow_flag, backward_arrow_flag, npoints, ...
                forward_arrow, backward_arrow];
                return rval;
  }

  if (exist(param.line_style))
  { line_style = __line_style(param.line_style); }

  if (class(param.thickness)=="num")
  {
    if (param.thickness >= 0)
    { thickness = floor(param.thickness); }
  }

  if (exist(param.pen_color))
  {
    pen_color = __find_or_manage_colors_in_color_table(param.pen_color);
  }

  if (exist(param.fill_color))
  {
    fill_color = __find_or_manage_colors_in_color_table(param.fill_color);
    area_fill = 20;
  }

  if (class(param.depth)=="num")
  {
    if (param.depth >= 1)
    { depth = floor(param.depth); }
  }

  if (class(param.pen_style)=="num")
  {
    if (param.pen_style >= -1)
    { pen_style = floor(param.pen_style); }
  }

  if (class(param.area_fill)=="num")
  {
    if (param.area_fill >= -1)
    { area_fill = floor(param.area_fill); }
  }

  if (class(param.style_val)=="num")
  {
    if (param.style_val > 0)
    { style_val = floor(param.style_val); }
  }

  if (class(param.join_style)=="num")
  {
    if (param.join_style >= 0)
    { join_style = floor(param.join_style); }
  }

  if (class(param.cap_style)=="num")
  {
    if (param.cap_style >= 0)
    { cap_style = floor(param.cap_style); }
  }

  if (class(param.radius)=="num")
  {
    if (param.radius > 0)
    { radius = floor(param.radius); }
  }

  if (class(param.forward_arrow)=="num")
  {
    if (param.forward_arrow.n == 5)
    {
      forward_arrow = param.forward_arrow;
    }
  }
  if (all(!isnan(forward_arrow)))
  {
    forward_arrow_flag = 1;
  }

  if (class(param.backward_arrow)=="num")
  {
    if (param.backward_arrow.n == 5)
    {
      backward_arrow = param.backward_arrow;
    }
  }
  if (all(!isnan(backward_arrow)))
  {
    backward_arrow_flag = 1;
  }

  // let's get out of here
  rval.desc = [sub_type, line_style, thickness, pen_color, fill_color, depth, ...
      pen_style, area_fill, style_val, join_style, ...
          cap_style, radius, forward_arrow_flag, backward_arrow_flag, npoints, ...
              forward_arrow, backward_arrow];
              return rval;
};



__box = function(xy, param)
{
  //
  if (!exist(xy))
  { return []; }
  if (xy.nc != 2)
  { return []; }
  if (xy.nr != 5)
  { return []; }

  sub_type = 2; // box: immutable

  // default parameters
  line_style = -1;
  thickness = 1;
  pen_color = 0;
  fill_color = 0;
  depth = 50;
  pen_style = 0;
  area_fill = -1;
  style_val = 1;
  join_style = 0;
  cap_style = 0;
  radius = 0; // arc_box only
  forward_arrow_flag = 0;
  forward_arrow = nan(1,5);
  backward_arrow_flag = 0;
  backward_arrow = nan(1,5);

  //
  npoints = 5;

  // return value is a list of line-type objects
  rval = <<>>;
  rval.data = xy;

  //
  // process parameters
  //
  if (class(param)!="list")
  {
    rval.desc = [sub_type, line_style, thickness, pen_color, fill_color, depth, ...
        pen_style, area_fill, style_val, join_style, ...
            cap_style, radius, forward_arrow_flag, backward_arrow_flag, npoints, ...
                forward_arrow, backward_arrow];
                return rval;
  }

  if (exist(param.line_style))
  { line_style = __line_style(param.line_style); }

  if (class(param.thickness)=="num")
  {
    if (param.thickness >= 0)
    { thickness = floor(param.thickness); }
  }

  if (exist(param.pen_color))
  {
    pen_color = __find_or_manage_colors_in_color_table(param.pen_color);
  }

  if (exist(param.fill_color))
  {
    fill_color = __find_or_manage_colors_in_color_table(param.fill_color);
    area_fill = 20;
  }

  if (class(param.depth)=="num")
  {
    if (param.depth >= 1)
    { depth = floor(param.depth); }
  }

  if (class(param.pen_style)=="num")
  {
    if (param.pen_style >= -1)
    { pen_style = floor(param.pen_style); }
  }

  if (class(param.area_fill)=="num")
  {
    if (param.area_fill >= -1)
    { area_fill = floor(param.area_fill); }
  }

  if (class(param.style_val)=="num")
  {
    if (param.style_val > 0)
    { style_val = floor(param.style_val); }
  }

  if (class(param.join_style)=="num")
  {
    if (param.join_style >= 0)
    { join_style = floor(param.join_style); }
  }

  if (class(param.cap_style)=="num")
  {
    if (param.cap_style >= 0)
    { cap_style = floor(param.cap_style); }
  }

  if (class(param.radius)=="num")
  {
    if (param.radius > 0)
    { radius = floor(param.radius); }
  }

  if (class(param.forward_arrow)=="num")
  {
    if (param.forward_arrow.n == 5)
    {
      forward_arrow = param.forward_arrow;
    }
  }
  if (all(!isnan(forward_arrow)))
  {
    forward_arrow_flag = 1;
  }

  if (class(param.backward_arrow)=="num")
  {
    if (param.backward_arrow.n == 5)
    {
      backward_arrow = param.backward_arrow;
    }
  }
  if (all(!isnan(backward_arrow)))
  {
    backward_arrow_flag = 1;
  }

  // let's get out of here
  rval.desc = [sub_type, line_style, thickness, pen_color, fill_color, depth, ...
      pen_style, area_fill, style_val, join_style, ...
          cap_style, radius, forward_arrow_flag, backward_arrow_flag, npoints, ...
              forward_arrow, backward_arrow];
              return rval;
};



//
// create data set for ellipse/circle
//
__ellipse = function(xy, radii, arg1, arg2)
{
  //
  if (!exist(xy))
  { return []; }
  if (xy.nc != 2)
  { return []; }
  center_xy = xy;

  if (!exist(radii))
  { return []; }
  if ((radii.n != 2) && (radii.n != 1))
  { return []; }

  // default parameters
  line_style = -1;
  thickness = 1;
  pen_color = 0;
  fill_color = 0;
  depth = 50;
  pen_style = 0;
  area_fill = -1;
  style_val = 1;
  direction = 1;
  rr = [];
  start_xy = [0,0];
  end_xy = [0,0];

  // circle or ellipse
  sub_type = 1;
  rr = abs(radii);
  if (radii.n == 1)
  {
    rr = [rr, rr];
    sub_type = 3;
    if (exist(arg1)&&!exist(arg2))
    {
      angle = 0;
      param = arg1;
    }
    else if(exist(arg1)&&exist(arg2))
    {
      printf("ellipse: ignoring input parameter 'angle' for drawing circle !\n");
      angle = 0;
      param = arg2;
    }
  }
  else
  {
    if (exist(arg1))
    { angle = arg1; }
    if (exist(arg2))
    { param = arg2; }
  }

  // angle
  if (!exist(angle))
  {
    angle = 0;
  }

  //
  // process parameters
  //
  if (class(param)!="list")
  {
    // minimal data set:
    return [sub_type, line_style, thickness, pen_color, fill_color, depth, ...
        pen_style, area_fill, style_val, direction, angle, ...
            center_xy, rr, start_xy, end_xy];
  }

  if (exist(param.line_style))
  { line_style = __line_style(param.line_style); }

  if (class(param.thickness)=="num")
  {
    if (param.thickness >= 0)
    { thickness = floor(param.thickness); }
  }

  if (exist(param.pen_color))
  {
    pen_color = __find_or_manage_colors_in_color_table(param.pen_color);
  }

  if (exist(param.fill_color))
  {
    fill_color = __find_or_manage_colors_in_color_table(param.fill_color);
    area_fill = 20;
  }

  if (class(param.depth)=="num")
  {
    if (param.depth >= 1)
    { depth = floor(param.depth); }
  }

  if (class(param.pen_style)=="num")
  {
    if (param.pen_style >= -1)
    { pen_style = floor(param.pen_style); }
  }

  if (class(param.area_fill)=="num")
  {
    if (param.area_fill >= -1)
    { area_fill = floor(param.area_fill); }
  }

  if (class(param.style_val)=="num")
  {
    if (param.style_val > 0)
    { style_val = floor(param.style_val); }
  }

  if (class(param.dir)=="num")
  {
    if (param.dir> 0)
    { direction = floor(param.dir); }
  }

  if (class(param.start)=="num")
  {
    if (param.start.n == 2)
    { start_xy = floor(param.start); }
  }

  if (class(param.end)=="num")
  {
    if (param.end.n == 2)
    { end_xy = floor(param.end); }
  }

  return [sub_type, line_style, thickness, pen_color, fill_color, depth, ...
      pen_style, area_fill, style_val, direction, angle, ...
          center_xy, rr, start_xy, end_xy];
};


//
//
//
//
//
__label = function(s, xy, param)
{
  // default values:
  just = 0;
  colour = -1;
  depth = 50;
  pen_style = 0; // unused
  font_flags = 4; // post script
  font = -1;
  font_size = 10;
  angle = 0;
  length = -1;
  height = -1;

    // what did user provide
  if (strlen(s)<1)
  {
    return [];
  }
  if (!exist(xy))
  {
    return [];
  }
  if (xy.nc!=2)
  {
    return [];
  }

  if (class(param)!="list")
  {
    // add piece of text:
    return <<data=[just, colour, depth, pen_style, font, font_size, angle, font_flags, height, length, xy];text=s>>;
  }

  if (exist(param.just))
  {
    if (param.just == 0 || param.just == 1 || param.just == 2)
    { just = param.just; }
  }

  if (exist(param.font_color))
  {
    colour = __find_or_manage_colors_in_color_table(param.font_color);
  }

  if (exist(param.depth))
  {
    if ((param.depth >= 1)&& (param.depth <= 999))
    { depth = param.depth; }
  }

  if (exist(param.font_flags))
  {
    font_flags = 1;
    if (strindex(tolower(param.font_flags),"lat"))
    {
      font_flags = 2;
    }
    else if (strindex(tolower(param.font_flags),"pos"))
    {
      font_flags = 4;
    }
    if (strindex(tolower(param.font_flags),"hid"))
    {
      font_flags = font_flags + 8;
    }
  }

  if (exist(param.font))
  {
    if (font_flags == 2 || font_flags==10)
    {
      // latex fonts
      if (strindex(tolower(param.font),"rom"))
      {
        font = 1;
      }
      else if (strindex(tolower(param.font),"bo"))
      {
        font = 2;
      }
      else if (strindex(tolower(param.font),"it"))
      {
        font = 3;
      }
      else if (strindex(tolower(param.font),"sa"))
      {
        font = 4;
      }
      else if (strindex(tolower(param.font),"ty"))
      {
        font = 5;
      }
      else
      {
        font = 0;
      }
    }
    else if (font_flags == 4 || font_flags==12)
    {
      // postscript fonts
      if (strindex(tolower(param.font),"tim"))
      {
        // times
        font = -2;
        font_base = 0;
      }
      else if (strindex(tolower(param.font),"ava"))
      {
        // AvantGarde
        font = -2;
        font_base = 4;
      }
      else if (strindex(tolower(param.font),"boo"))
      {
        // Bookman
        font = -2;
        font_base = 8;
      }
      else if (strindex(tolower(param.font),"cou"))
      {
        // Courier
        font = -2;
        font_base = 12;
      }
      else if (strindex(tolower(param.font),"helv")&&strindex(tolower(param.font),"nar"))
      {
        // Helvetica Narrow
        font = -2;
        font_base = 20;
      }
      else if (strindex(tolower(param.font),"helv"))
      {
        // Helvetica
        font = -2;
        font_base = 16;
      }
      else if (strindex(tolower(param.font),"new"))
      {
        // New Century
        font = -2;
        font_base = 24;
      }
      else if (strindex(tolower(param.font),"pal"))
      {
        // Palatino
        font = -2;
        font_base = 28;
      }
      else if (strindex(tolower(param.font),"sym"))
      {
        // Helvetica
        font_base = -1;
        font = 32;
      }
      else if (strindex(tolower(param.font),"zap"))
      {
        // Helvetica
        font = -2;
        font_base = 33;
      }
    }

    if (font == -2)
    {
      font = font_base;
      if (strindex(tolower(param.font),"it") && strindex(tolower(param.font),"bo"))
      {
        font = font + 3;
      }
      else if (strindex(tolower(param.font),"it") || strindex(tolower(param.font),"ob") || strindex(tolower(param.font),"di"))
      {
        font = font + 1;
      }
      else if (strindex(tolower(param.font),"bo") || strindex(tolower(param.font),"de"))
      {
        font = font + 2;
      }
    }
  }

  if (exist(param.font_size))
  {
    if (param.font_size >= 1)
    { font_size = floor(param.font_size); }
  }

  if (exist(param.angle))
  {
    angle = param.angle;
  }

  if (exist(param.height))
  {
    if (param.height>= -1)
    { height = param.height; }
  }

  if (exist(param.length))
  {
    if (param.length>= -1)
    { length = param.length; }
  }

  return <<desc=[just, colour, depth, pen_style, font, font_size, 3.1415/180*angle, font_flags, ...
      height, length, xy];data=s>>;
};


__papersize = function(val)
{
  rval = blank();
  if (strlen(val)>0)
  {
    if (strindex(tolower(val), "let"))
    {
      rval = "Letter";
    }
    else if (strindex(tolower(val), "leg"))
    {
      rval = "Legal";
    }
    else if (strindex(tolower(val), "led"))
    {
      rval = "Ledger";
    }
    else if (strindex(tolower(val), "tab"))
    {
      rval = "Tabloid";
    }
    else if (tolower(val)=="a")
    {
      rval = "A";
    }
    else if (tolower(val)=="b")
    {
      rval = "B";
    }
    else if (tolower(val)=="c")
    {
      rval = "C";
    }
    else if (tolower(val)=="d")
    {
      rval = "D";
    }
    else if (tolower(val)=="e")
    {
      rval = "E";
    }
    else if (tolower(val)=="a4")
    {
      rval = "A4";
    }
    else if (tolower(val)=="a3")
    {
      rval = "A3";
    }
    else if (tolower(val)=="a2")
    {
      rval = "A2";
    }
    else if (tolower(val)=="a1")
    {
      rval = "A1";
    }
    else if (tolower(val)=="a0")
    {
      rval = "A0";
    }
    else if (tolower(val)=="b5")
    {
      rval = "B5";
    }
  }

  return rval;
};

// ******************************************************
//
//
// X F I G   C L A S S:
//
//
// ******************************************************
xfig = classdef( )
{
  static(_orientation,_justification,_units,_papersize,_magnification,_multiple_pages,_transparent_color);
  static(_resolution,_coord_system);
  static(_filename,text_table_data,text_table_text);
  static(ellipse_data,polyline_data,polyline_desc,picture_desc,picture_file);
  static(world_origin,world_bbox,world_xfig_origin,world_xfig_dir,world_xfig_scale,arc_data,_comment);

  //
  _comment = "";
  _orientation = "Landscape";
  _justification = "Center";
  _units = "Inches";
  _papersize = "Letter";
  _magnification = 100;
  _multiple_pages = "Single";
  _transparent_color = -2;
  _resolution = 1200;
  _coord_system = 2;
  text_table_desc = [];
  text_table_data = blank();
  ellipse_data = [];
  //
  arc_data = [];
  //
  polyline_data = cell();
  polyline_desc = [];
  //
  picture_desc = [];
  picture_file = blank();

  // default is that paper is the world
  world_origin = [0.0, 0.0];
  world_bbox = [];
  world_xfig_origin = [0.0, 0.0];
  world_xfig_scale = 1.0;
  world_xfig_dir = [1, 1];

  //
  public(orientation,justification,units,papersize,magnification,transparent_color);
  public(resolution,color,export,label,ellipse,line,picture,draw_scale,world,box,arc,comment);

  world = function (world_xy_0, xfig_xy_0, xfig_xy_dir, xy_scale, bbox)
  {
    // what is the world origin 
    if (class(world_xy_0)=="num" && length(world_xy_0)==2)
    { world_origin = world_xy_0; }

    // where is world origin is on xfig plane:
    if (class(xfig_xy_0)=="num" && length(xfig_xy_0)==2)
    { world_xfig_origin = xfig_xy_0; }

    // direction
    if (class(xfig_xy_dir)=="num" && length(xfig_xy_dir)==2)
    { world_xfig_dir = xfig_xy_dir; }

    // what about scale
    if (isnumber(xy_scale))
    { world_xfig_scale = xy_scale; }

    // what about bounding box
    if (isnumber(bbox) && length(bbox)==2)
    { world_bbox = [world_origin; world_origin + bbox]; }

    return <<origin=world_origin;xfig_origin=world_xfig_origin;xfig_direction=world_xfig_dir;...
        xfig_scale=world_xfig_scale;bbox=world_bbox>>;
  };

  //
  // functions that modify class-member properties:
  //
  orientation = function(val)
  {
    if (strlen(val)>0)
    {
      if (strindex(tolower(val), "land"))
      {
        _orientation = "Landscape";
      }
      if (strindex(tolower(val), "por"))
      {
        _orientation = "Portrait";
      }
    }
    return _orientation;
  };

  comment = function(val)
  {
    if (strlen(val)>0)
    {
      _comment = val;
    }
    return _comment;
  };

  justification = function(val)
  {
    if (strlen(val)>0)
    {
      if (strindex(tolower(val), "cen"))
      {
        _justification = "Center";
      }
      else if (strindex(tolower(val), "lef"))
      {
        _justification = "Flush Left";
      }
    }
    return _justification;
  };

  units = function(val)
  {
    if (strlen(val)>0)
    {
      if (strindex(tolower(val), "m"))
      {
        _units = "Metric";
      }
      else if (strindex(tolower(val), "in"))
      {
        _units = "Inches";
      }
    }
    return _units;
  };

  papersize = function ( val )
  {
    if (strlen(val)>0)
    {
      rval = __papersize(val);
      if (strlen(val)>0)
      {
        _papersize = rval;
      }
    }

    return _papersize;
  };

  magnification = function(val)
  {
    if (type(val)=="real" || type(val)=="int")
    {
      if (val>0)
      {
        _magnification = val;
      }
    }
    return _magnification;
  };

  transparent_color = function(val)
  {
    if (type(val)=="real" || type(val)=="int")
    {
      val = floor(val);
      if ((val>=-2)&&(val<=512))
      {
        _transparent_color = val;
      }
    }
    return _transparent_color;
  };

  resolution = function(val)
  {
    if (type(val)=="real" || type(val)=="int")
    {
      val = floor(val);
      if (val>=1)
      {
        _resolution = val;
      }
    }
    return _resolution;
  };

  color = __find_or_manage_colors_in_color_table;

  label = function(s, xy, param)
  {
    use_world = 0;
    if (type(param) == "list")
    {
      if (class(param.world_coord)=="num")
      {
        use_world = (param.use_world > 0);
      }
    }
    if (use_world)
    {
      belong = __clip2world(xy, world_bbox);
      i_min = 1;
      i_max = max(belong);
      for (i in i_min:i_max)
      {
        idx_i = find(belong == i);
        xy_i = world_xfig_origin + (xy[idx_i;] - world_origin) .* world_xfig_dir .* world_xfig_scale;
        rval = __label(s, xy_i, param);
        if (type(rval)=="list")
        {
          text_table_desc = [text_table_desc; rval.desc];
          text_table_data = [text_table_data; blank(rval.desc.nr,1,rval.data)];
        }
      }
    }
    else
    {
      rval = __label(s, xy, param);
      if (type(rval)=="list")
      {
        text_table_desc = [text_table_desc; rval.desc];
        text_table_data = [text_table_data; blank(rval.desc.nr,1,rval.data)];
      }
    }

    return <<data=text_table_data;desc=text_table_desc>>;
  };

  ellipse = function(xy, radii, angle, params)
  {
    use_world = 0;
    if (type(param) == "list")
    {
      if (class(param.world_coord)=="num")
      {
        use_world = (param.world_coord > 0);
      }
    }
    if (use_world)
    {
      xy = world_xfig_origin + (xy - world_origin) .* world_xfig_dir .* world_xfig_scale;
      radii = radii * world_scale;
    }
    e = __ellipse(xy, radii, angle, params);
    if (!isempty(e))
    {
      ellipse_data = [ellipse_data; e];
    }
    return ellipse_data;
  };

  line = function(xy, param)
  {
    // if xy is in world coordinates, convert that to
    // xfig coordinates
    use_world = 0;
    if (type(param) == "list")
    {
      if (class(param.world_coord)=="num")
      {
        use_world = (param.world_coord > 0);
      }
    }
    if (use_world)
    {
      // do clipping: belong contains sequences of increasing integers interrupted by 0's
      belong = __clip2world(xy, world_bbox);
      i_min = 1;
      i_max = max(belong);
      for (i in i_min:i_max)
      {
        idx_i = find(belong == i);
        xy_0_i = [];
        xy_f_i = [];
        if (min(idx_i)>1)
        {
          // segment does not start with the first point in dataset
          idcs = min(idx_i) + [-1, 0];
          ls1 = xy[idcs;];
          xy_0_i = __intersection_line_bbox(ls1, world_bbox);
        }
        if (max(idx_i)<xy.nr)
        {
          // segment does not end with the last point in dataset
          idcs = max(idx_i) + [0,1];
          ls1 = xy[idcs;];
          xy_f_i = __intersection_line_bbox(ls1, world_bbox);
        }

        // do not plot single points
        if ((length(idx_i) < 2) && (length(xy_0_i)==0) && (length(xy_f_i)==0))
        { continue; }

        // convert to xfig coordinates where the lines start from the bounding box edges
        xy_i = world_xfig_origin + ([xy_0_i;xy[idx_i;];xy_f_i] - world_origin) .* world_xfig_dir .* world_xfig_scale;
        l = __line(xy_i,param);
        if (!isempty(l))
        {
          polyline_data = [polyline_data; l.data];
          polyline_desc = [polyline_desc; l.desc];
        }
      }
    }
    else
    {
      l = __line(xy,param);
      if (!isempty(l))
      {
        polyline_data = [polyline_data; l.data];
        polyline_desc = [polyline_desc; l.desc];
      }
    }
    return <<data=polyline_data;desc=polyline_desc>>;
  };

  box = function(xy_box, param)
  {
    // xy_box = [x_min, y_min; x_max,y_max]
    xy = zeros(5,2);
    xy[1;] = xy_box[1;];
    xy[2;] = [xy_box[2;1], xy_box[1;2]];
    xy[3;] = xy_box[2;];
    xy[4;] = [xy_box[1;1], xy_box[2;2]];
    xy[5;] = xy_box[1;];

    // if xy is in world coordinates, convert that to
    // xfig coordinates
    use_world = 0;
    if (type(param) == "list")
    {
      if (class(param.world_coord)=="num")
      {
        use_world = (param.world_coord > 0);
      }
    }
    if (use_world)
    {
      // do clipping: belong contains sequences of increasing integers interrupted by 0's
      belong = __clip2world(xy, world_bbox);
      i_min = 1;
      i_max = max(belong);
      for (i in i_min:i_max)
      {
        idx_i = find(belong == i);
        xy_0_i = [];
        xy_f_i = [];
        if (min(idx_i)>1)
        {
          // segment does not start with the first point in dataset
          idcs = min(idx_i) + [-1, 0];
          ls1 = xy[idcs;];
          xy_0_i = __intersection_line_bbox(ls1, world_bbox);
        }
        if (max(idx_i)<xy.nr)
        {
          // segment does not end with the last point in dataset
          idcs = max(idx_i) + [0,1];
          ls1 = xy[idcs;];
          xy_f_i = __intersection_line_bbox(ls1, world_bbox);
        }

        // do not plot single points
        if ((length(idx_i) < 2) && (length(xy_0_i)==0) && (length(xy_f_i)==0))
        { continue; }

        // convert to xfig coordinates where the lines start from the bounding box edges
        xy_i = world_xfig_origin + ([xy_0_i;xy[idx_i;];xy_f_i] - world_origin) .* world_xfig_dir .* world_xfig_scale;
        if (i_min != i_max)
        { l = __line(xy_i,param); }
        else
        { l = __box(xy_i,param); }
        if (!isempty(l))
        {
          polyline_data = [polyline_data; l.data];
          polyline_desc = [polyline_desc; l.desc];
        }
      }
    }
    else
    {
      l = __box(xy,param);
      if (!isempty(l))
      {
        polyline_data = [polyline_data; l.data];
        polyline_desc = [polyline_desc; l.desc];
      }
    }
    return <<data=polyline_data;desc=polyline_desc>>;
  };

  arc = function(center_xy, xy, param)
  {
    // if xy is in world coordinates, convert that to xfig coordinates - no clipping possible
    use_world = 0;
    sub_type = 1;
    if (type(param) == "list")
    {
      if (class(param.world_coord)=="num")
      {
        use_world = (param.world_coord > 0);
      }
      if (strlen(param.arc_type)>0)
      {
        if (strindex(param.arc_type,"we")>0 || strindex(param.arc_type,"cl")>0)
        {
          sub_type = 2;
        }
      }
    }
    if (use_world)
    {
      // convert to xfig coordinates 
      center_xy = world_xfig_origin + (center_xy - world_origin) .* world_xfig_dir .* world_xfig_scale;
      xy = world_xfig_origin + (xy - world_origin) .* world_xfig_dir .* world_xfig_scale;
    }
    l = __arc(center_xy,xy,sub_type,param);
    if (!isempty(l))
    {
      arc_data = [arc_data; l];
    }

    return arc_data;
  };

  picture = function(fn, xy_0, xy_d, orient, param)
  {
    p = __picture (fn, xy_0, xy_d, orient, param);
    if (!isempty(p))
    {
      picture_desc = [picture_desc; p.desc];
      picture_file = [picture_file; p.file];
    }
    return <<file=picture_file;desc=picture_desc>>;
  };

  draw_scale = function(pos, len_xfig, len_world, param)
  {
    unit_name = "";
    frac_up = 0.25;
    //
    if (!all(isnumber(pos)))
    { return 1; }
    if (!isnumber(len_xfig))
    { return 1; }
    if (!isnumber(len_world))
    { return 1; }

    if (class(param)=="list")
    {
      if (exist(param.frac_up))
      { frac_up = param.frac_up; }
      if (strlen(param.unit_name)>0)
      { unit_name = param.unit_name; }
    }

    // create scale
    pts = zeros(4,2);
    pts[1;] = pos - [0, frac_up * len_xfig];
    pts[2;] = pos;
    pts[3;] = pos + [len_xfig, 0];
    pts[4;] = pts[3;] - [0, frac_up * len_xfig];
    // call line with this
    l = __line(pts,param);
    if (!isempty(l))
    {
      polyline_data = [polyline_data; l.data];
      polyline_desc = [polyline_desc; l.desc];
    }

    // put text at two ends
    t = __label("0", pts[1;] - [0, frac_up * len_xfig], param);
    if (type(t)=="list")
    {
      text_table_desc = [text_table_desc; t.desc];
      text_table_data = [text_table_data; blank(t.desc.nr,1,t.data)];
    }
    t = __label( num2str(len_world,"%g" + unit_name), pts[4;] - [0, frac_up * len_xfig], param);
    if (type(t)=="list")
    {
      text_table_desc = [text_table_desc; t.desc];
      text_table_data = [text_table_data; blank(t.desc.nr,1,t.data)];
    }

    return 0;
  };

  export = function ( fn, fmt )
  {
    use_fig2dev=0;
    if (strlen(fn) > 0)
    {
      _filename = fn;
      if (strlen(fmt)>0)
      { use_fig2dev=1; }
    }

    if (strlen(_filename)<1)
    {
      _filename = "stdout";
    }

    if (_filename != "stdout")
    { open(_filename,"w"); }

    //
    // standard preamble
    //
    fprintf(_filename,"#FIG 3.2 %s\n", _comment);
    fprintf(_filename,"%s\n", _orientation);
    fprintf(_filename,"%s\n", _justification);
    fprintf(_filename,"%s\n", _units);
    fprintf(_filename,"%s\n", _papersize);
    fprintf(_filename,"%.0f\n", _magnification);
    fprintf(_filename,"%s\n", _multiple_pages);
    fprintf(_filename,"%.0f\n", _transparent_color);
    fprintf(_filename,"%.0f %.0f\n\n", _resolution, _coord_system);

    // colors: 0 32 #aabbcc
    if (xfig_color_table.id.n>0)
    {
      for (i in 1:(xfig_color_table.id.n))
      {
        fprintf(_filename,"0 %.0f #", xfig_color_table.id[i]);
        fprintf(_filename,"%02x", int(xfig_color_table.rgb[i;1]));
        fprintf(_filename,"%02x", int(xfig_color_table.rgb[i;2]));
        fprintf(_filename,"%02x\n", int(xfig_color_table.rgb[i;3]));
      }
      fprintf(_filename,"\n");
    }

    // print all circles
    for (i in 1:(ellipse_data.nr))
    {
      fprintf(_filename,"1 ");
      for (k in 1:8)
      {
        fprintf(_filename," %.0f", ellipse_data[i;k]); // 
      }
      fprintf(_filename," %g", ellipse_data[i;9]);  // style_val
      fprintf(_filename," %.0f", ellipse_data[i;10]);
      fprintf(_filename," %g", 3.1415/180*ellipse_data[i;11]); // angle
      fprintf(_filename," %.0f", _resolution * ellipse_data[i;12]); // center
      fprintf(_filename," %.0f", _resolution * ellipse_data[i;13]);
      fprintf(_filename," %.0f", _resolution * ellipse_data[i;14]); // radius
      fprintf(_filename," %.0f", _resolution * ellipse_data[i;15]);
      fprintf(_filename," %.0f", _resolution * ellipse_data[i;16]); // start
      fprintf(_filename," %.0f", _resolution * ellipse_data[i;17]);
      fprintf(_filename," %.0f", _resolution * ellipse_data[i;18]); // end
      fprintf(_filename," %.0f\n\n", _resolution * ellipse_data[i;19]);
    }

    // print all lines
    for (i in 1:(polyline_desc.nr))
    {
      fprintf(_filename,"2 ");
      for (k in 1:8)
      {
        fprintf(_filename," %.0f", polyline_desc[i;k]); //
      }
      fprintf(_filename," %g", polyline_desc[i;9]);  // style_val
      for (k in 10:15)
      {
        fprintf(_filename," %.0f", polyline_desc[i;k]); //
      }
      // forward arrow
      if (polyline_desc[i;13])
      {
        fprintf(_filename,"\n\t");
        for (k in 16:20)
        {
          fprintf(_filename," %.0f", polyline_desc[i;k]); //
        }
      }
      // backward arrow
      if (polyline_desc[i;14])
      {
        fprintf(_filename,"\n\t");
        for (k in 21:25)
        {
          fprintf(_filename," %.0f", polyline_desc[i;k]); //
        }
      }
      fprintf(_filename,"\n"); // new line
      for (k in 1:(polyline_data[i].nr))
      {
        fprintf(_filename,"    %.0f", _resolution * polyline_data[i][k;1]); // x_n
        fprintf(_filename,   " %.0f\n", _resolution * polyline_data[i][k;2]); // y_n
      }
      fprintf(_filename,"\n"); // new line
    }

    // print all pictures
    for (i in 1:(picture_desc.nr))
    {
      fprintf(_filename,"2 ");
      for (k in 1:8)
      {
        fprintf(_filename," %.0f", picture_desc[i;k]); //
      }
      fprintf(_filename," %g", picture_desc[i;9]);  // style_val
      for (k in 10:15)
      {
        fprintf(_filename," %.0f", picture_desc[i;k]); //
      }
      fprintf(_filename,"\n %.0f", picture_desc[i;16]); //
      fprintf(_filename," %s\n", picture_file[i]); // new line
      for (k in 1:5)
      {
        fprintf(_filename,"    %.0f",   _resolution * picture_desc[i;16+2*k-1]); // x_n
        fprintf(_filename,   " %.0f\n", _resolution * picture_desc[i;16+2*k  ]); // y_n
      }
      fprintf(_filename,"\n");
    }

    // print all texts
    for (i in 1:(text_table_desc.nr))
    {
      fprintf(_filename,"4 ");
      for (k in 1:5)
      {
        fprintf(_filename," %.0f", text_table_desc[i;k]); // just, colour, depth, pen_style, font
      }
      fprintf(_filename," %g", text_table_desc[i;6]);   // font_size
      fprintf(_filename," %g", text_table_desc[i;7]);   // angle
      fprintf(_filename," %.0f", text_table_desc[i;8]); // font_flags
      fprintf(_filename," %g", text_table_desc[i;9]);   // height
      fprintf(_filename," %g", text_table_desc[i;10]);  //length
      fprintf(_filename," %.0f", _resolution * text_table_desc[i;11] - 0.6 * text_table_desc[i;6]);
      fprintf(_filename," %.0f", _resolution * text_table_desc[i;12] + 0.5 * text_table_desc[i;6]);
      fprintf(_filename," %s\001\n\n", text_table_data[i]);
    }

    // print all arcs
    for (i in 1:(arc_data.nr))
    {
      fprintf(_filename,"5 ");
      for (k in 1:8)
      {
        fprintf(_filename," %.0f", arc_data[i;k]); //
      }
      fprintf(_filename," %g", arc_data[i;9]);  // style_val
      for (k in 10:13)
      {
        fprintf(_filename," %.0f", arc_data[i;k]); //
      }
      // center
      fprintf(_filename," %g", _resolution * arc_data[i;14]);
      fprintf(_filename," %g", _resolution * arc_data[i;15]);
      // data pts
      fprintf(_filename," %.0f", _resolution * arc_data[i;16]); // x_1
      fprintf(_filename," %.0f", _resolution * arc_data[i;17]); // y_1
      fprintf(_filename," %.0f", _resolution * arc_data[i;18]); // x_2
      fprintf(_filename," %.0f", _resolution * arc_data[i;19]); // y_2
      fprintf(_filename," %.0f", _resolution * arc_data[i;20]); // x_3
      fprintf(_filename," %.0f", _resolution * arc_data[i;21]); // y_3
      // forward arrow
      if (arc_data[i;12])
      {
        fprintf(_filename,"\n\t");
        for (k in 22:26)
        {
          fprintf(_filename," %.0f", arc_data[i;k]); //
        }
      }
      // backward arrow
      if (arc_data[i;13])
      {
        fprintf(_filename,"\n\t");
        for (k in 27:31)
        {
          fprintf(_filename," %.0f", arc_data[i;k]); //
        }
      }
      fprintf(_filename,"\n\n"); // new line
    }

    if (_filename != "stdout")
    { close(_filename); }

    if (use_fig2dev)
    {
      sleep(1);
      _fig2dev_fn = gsub("."+fmt, ".fig", _filename).string;
      system("fig2dev -L " + fmt + " " + _filename + " " + _fig2dev_fn + " 1>/dev/null 2>/dev/null");
    }
  };

};









_INIT_ = 1;

