// libooo.r - rlabplus library to create OpenOffice.org documents
// based on ooolib project by J. Colton
// # ooolib - This perl library is built to create OpenOffice.org documents.
// # Copyright (C) 2003  Joseph Colton
//
// # This library is free software; you can redistribute it and/or
// # modify it under the terms of the GNU Lesser General Public
// # License as published by the Free Software Foundation; either
// # version 2.1 of the License, or (at your option) any later version.
//
// # This library is distributed in the hope that it will be useful,
// # but WITHOUT ANY WARRANTY; without even the implied warranty of
// # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// # Lesser General Public License for more details.
//
// # You should have received a copy of the GNU Lesser General Public
// # License along with this library; if not, write to the Free Software
// # Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA
//
// # You can contact me by email at josephcolton@gmail.com
//
// rlabplus version of the library by Marijan Kostrun, XI-2008
// mkostrun@gmail.com
//
// user functions provided by the library:
//  ooo.init         -> initialize the library for a single ooo document
//  ooo.special      -> specifies special structures in the document, lists and such
//  ooo.set          -> set values of different parameters of the document
//  ooo.data         -> set content of the document
//  ooo.generate     -> write the complete document to a file system


//
// local functions: ignorance is a bliss
//
static(_oooWriteMeta, _oooTimeStamp, _oooWriteManifest, _oooWriteStyles);
static(_oooWriteSettings, _oooWriteContent, _oooWriteMimetype);
static(_ooosxcCellLocConv, _oooStyleName, _ooolib_sxc_cell_StyleName);
static(_oooStyleName, _ooolib_sxc_cell_StyleName);

//
// global variables for the library
//
static(_builddir, _ooodebug, _ooolibrary, _oooversion, _ooolanguage, _BUILDDIR_, _OOOEXEC_, _ooolib_FloatFormat);
_BUILDDIR_ = getenv("PWD") + "/.ooo_tmp";
_builddir = _BUILDDIR_;
_ooodebug = 0;
_ooolibrary = "rlabplus/libooo";
_oooversion = "1.0.2";
_ooolanguage = "en-US";
_ooolib_FloatFormat = "%g";

if (!exist(_OOOEXEC_))
{
  WHICH=reads("|which which");
  if (strlen(WHICH)>0)
  {
    for (soff in ["soffice", "libreoffice4.2", "libreoffice4.4", "libreoffice"] )
    {
      // check for office executable
      _OOOEXEC_ = reads("|" + WHICH + " " + soff + " 2>/dev/null");
      if (isempty(_OOOEXEC_))
      {
        clear(_OOOEXEC_);
      }
    then:
      if (isempty(_OOOEXEC_))
      {
        printf("HIE: Cannot find the open office executable necessary for this library!\n");
        stop("HIE: Please locate its name, then modify the entry in the library!\n");
      }
    }
  }
}


// -----------------------------------------------------------------------
//
// document related variables
//
// -----------------------------------------------------------------------

// lists / hashes
static(_ooo_i, _cellhash, _options);
_ooo_i = <<>>;
_cellhash = <<>>;   // %cellhash
_options  = <<>>;   // %options

// text arrays
static(_fontdecls, _keywords, _documenttext, _autostyles);
_fontdecls    = blank(0,0);   //  @fontdecls
_keywords     = blank(0,0);   //  @keywords
_documenttext = blank(0,0);   //  @documenttext
_autostyles   = blank(0,0);   //  @autostyles

// other internal variables
static(_minx, _miny, _maxx, _maxy);
_minx = 1;
_miny = 1;
_maxx = 32000;
_maxy = 32000;

static(_ooolib_sxc_cell_styles);
_ooolib_sxc_cell_styles = blank(0,0);   // @ooolib_sxc_cell_styles

// main function list
ooo = <<>>;

ooo.init = function (mykind, opt)
{
  if (mykind != "sxw" && mykind!="sxc" && mykind!="xls")
  { error("oooInit: the library can generate 'sxc' or 'sxw' type files"); }

  if (exist(opt))
  {
    if (opt == "debug")
    { _ooodebug = 1;}
  }

  // Set default values
  _cellhash.x = 1;
  _cellhash.y = 1;
  _cellhash.xmax = 1;
  _cellhash.ymax = 1;
  _fontdecls    = blank(0,0);
  _keywords     = blank(0,0);
  _documenttext = blank(0,0);
  _autostyles   = blank(0,0);
  _minx = 1;
  _miny = 1;
  _maxx = 32000;
  _maxy = 32000;
  _ooolib_sxc_cell_styles = blank(0,0);

  // Formatting
  _options.["nextstyle"] = 1;
  _options.["nextstyle table-column"] = 2;
  _options.["nextstyle table-cell"] = 1;
  _options.nextlist = 1;
  _options.justify = "left";
  _options.bold = "off";
  _options.italic = "off";
  _options.underline = "off";
  _options.textcolor = "000000";
  _options.textbgcolor = "FFFFFF";
  _options.textsize = "12";
  _options.textsize_default = "12";
  _options.textsize_min = "6";
  _options.textsize_max = "96";

  // Default Font
  _fontdecls = "<style:font-decl style:name=\"Times\" fo:" ...
      +"font-family=\"Times\" style:font-family-generic=\"roman\" style:" ...
      +"font-pitch=\"variable\"/>";
  _options.fontname = "Times";
  _options.["font-decl Times"] = _options.fontname;

  // Used for auto-incrementing cells after giving values
  _cellhash.xauto = 0;
  _cellhash.yauto = 0;

  // User defined meta variables
  _options.["info1 name"] = "Info 1";
  _options.["info2 name"] = "Info 2";
  _options.["info3 name"] = "Info 3";
  _options.["info4 name"] = "Info 4";

  // Set all knowns
  _options.kind = mykind;
  _options.time = time2dstr(seconds(),"%c");

  return 0;
};

_ooolib_sxc_cell_StyleName = function()
{
  _prop = blank(0,0);

  // Default style
  if (_options.bold!="on" && _options.italic!="on" && _options.textcolor=="000000")
  { return "ro1"; }

  if (_options.bold == "on")
  { _prop = [_prop, "fo:font-weight=\"bold\""]; }
  if (_options.italic == "on")
  { _prop = [_prop, "fo:font-style=\"italic\""]; }
  if (_options.textcolor != "000000")
  { _prop = [_prop, "fo:color=\"#"+_options.textcolor+"\""]; }

  if (length(_prop)>0)
  {
    _propline = "";
    for (i in _prop)
    { _propline = _propline + " " + i; }
  } else {
    _propline = " ";
  }

  if (exist(_options.["cell-style"+_propline]))
  { return _options.["cell-style"+_propline]; }

  _num = _options.["nextstyle table-cell"];
  _options.["nextstyle table-cell"] = _num + 1;
  _style = "ce"+text(_num);

  _ooolib_sxc_cell_styles = [_ooolib_sxc_cell_styles, ...
      "<style:style style:name=\""+_style+"\" style:family=\"table-cell\"" ...
      +" style:parent-style-name=\"Default\">"];
  _ooolib_sxc_cell_styles = [_ooolib_sxc_cell_styles, "<style:properties"+_propline+"/>"];
  _ooolib_sxc_cell_styles = [_ooolib_sxc_cell_styles, "</style:style>"];

  return _style;
};


ooo.data = function (style, origtext, numval, numfmt)
{
  if(!exist(style))
  { error("ooo.data: requires at least two arguments"); }

  // origtext can be a string array
  if (!exist(origtext))
  {
    _text = blank(0,0);
  } else {
    _text = blank(origtext);
  }

  if (exist(origtext))
  {
    for (i in range(origtext))
    { _text[i] = _oooCleanText( origtext[i] ); }
  }
  if (exist(numval)&&!exist(numfmt))
  { numfmt = _ooolib_FloatFormat; }

  if (style == "h")
  {
    _stylename = _oooStyleName("Heading");
    _documenttext = [_documenttext, ...
        "<text:p text:style-name=\""+_stylename+"\">"+_text+"</text:p>"];
  }
  if (substr(style,1)=="h" && ascii(substr(style,2))>=49 && ...
      ascii(substr(style,2))<=57)
  {
    _stylename = _oooStyleName("Heading "+ substr(style,2));
    _documenttext = [_documenttext, "<text:h text:style-name=\"" ...
        +_stylename+"\" text:level=\""+substr(style,2)+"\">"+_text[1]+"</text:h>"];
  }
  if (style == "default")
  {
    _stylename = _oooStyleName("Standard");
    if (exist(_text))
    {
      // Default Text
      _documenttext = [_documenttext, "<text:p text:style-name=\""...
          +_stylename+"\">"+_text+"</text:p>"];
    } else {
      // Blank line
      _documenttext = [_documenttext,  "<text:p text:style-name=\""...
          +_stylename+"\"/>"];
    }
  }
  if (style == "textbody")
  {
    // textbody does some special formatting stuff.
    _stylename = _oooStyleName("Text body");
    if (any(_text!=""))
    {
      // A paragraph with text
      for (_paragraph in _text)
      {
        // Send Begin Paragraph tag
        if(strlen(_paragraph)==0)
        { continue; }
        // replace spaces in text with their xml
        _documenttext = [_documenttext, "<text:p text:style-name=\""+_stylename+"\">"];
        _documenttext = [_documenttext, gsub(". <text:s/>", ". ", _paragraph).string];
        _documenttext = [_documenttext, "</text:p>"];
      }
    } else {
      // This is a blank paragraph tag.
      _documenttext = [_documenttext, "<text:p text:style-name=\""+_stylename+"\"/>"];
    }
  }

  if (style == "cell-float")
  {
    x = _cellhash.x;
    y = _cellhash.y;

    _cellhash.[text(x)+" "+text(y)+" type"]   = "float";
    _cellhash.[text(x)+" "+text(y)+" value"]  = _text[1];
    _cellhash.[text(x)+" "+text(y)+" format"] = _ooolib_FloatFormat;
    _cellhash.[text(x)+" "+text(y)+" style"]  = _ooolib_sxc_cell_StyleName();

    _oooCellUpdate();
  }
  if (style == "cell-text")
  {
    x = _cellhash.x;
    y = _cellhash.y;

    _cellhash.[text(x)+" "+text(y)+" type"]   = "text";
    _cellhash.[text(x)+" "+text(y)+" value"]  = _text[1];
    _cellhash.[text(x)+" "+text(y)+" format"] = "";
    _cellhash.[text(x)+" "+text(y)+" style"]  = _ooolib_sxc_cell_StyleName();

    _oooCellUpdate();
  }
  if (style == "cell-formula")
  {
    x = _cellhash.x;
    y = _cellhash.y;

    _cellhash.[text(x)+" "+text(y)+" type"]     = "formula";
    _cellhash.[text(x)+" "+text(y)+" formula"]  = _text[1];
    if (exist(numval) && exist(numfmt))
    {
      _cellhash.[text(x)+" "+text(y)+" value"]  = text(numval,"%e");
      _cellhash.[text(x)+" "+text(y)+" value-print"]  = text(numval,numfmt);
    }
    _cellhash.[text(x)+" "+text(y)+" format"] = _ooolib_FloatFormat;
    _cellhash.[text(x)+" "+text(y)+" style"]  = _ooolib_sxc_cell_StyleName();

    _oooCellUpdate();
  }
  if (style == "cell-skip")
  {
    // Does not do anything to the cell, just skips it
    _oooCellUpdate();
  }

  return "ok";
};


_oooStyleName = function (style)
{
    // Uses these numbers to keep track of which styles have been created
    _j = _options.justify;
    _b = _options.bold;
    _i = _options.italic;
    _u = _options.underline;
    _tc = _options.textcolor;
    _tb = _options.textbgcolor;
    _ts = _options.textsize;
    _fn = _options.fontname;

    // does style starts with h/H?: if ($style =~ /^h/)
    if (tolower(substr(style,1))=="h" && _j=="left")
    {
      return "";
    }

    if (_j=="left" && _b=="off" && _i=="off" && _u=="off" && ...
        _tc=="000000" && _tb=="FFFFFF" && _ts==_options.textsize_default && ...
        _fn=="Times")
    {
      // Nothing special needs to be done
      return "";
    }

    // A style needs to be looked up or created
    _idx = "style "+_j+" "+_b+" "+_i+" "+_u+" "+_tc+" "+_tb+" "+_ts+" "+_fn;
    // does style starts with h/H?: if ($style =~ /^h/)
    if (tolower(substr(style,1))=="h")
    {_idx = "style "+ _j;}

    if (exist(_options.[_idx]))
    { return _options.[_idx]; }

    _stylenum = "P"+text(_options.nextstyle);
    _options.nextstyle = _options.nextstyle + 1;
    _options.[_idx] = _stylenum;

    // justification:
    if (_j == "left")
    {_jt = "fo:text-align=\"start\" style:justify-single-word=\"false\" ";}
    if (_j == "right")
    {_jt = "fo:text-align=\"end\" style:justify-single-word=\"false\" ";}
    if (_j == "center")
    {_jt = "fo:text-align=\"center\" style:justify-single-word=\"false\" ";}
    if (_j == "block")
    {_jt = "fo:text-align=\"justify\" style:justify-single-word=\"false\" ";}

    // bold:
    _bt = "";
    if (_options.bold == "on")
    {
      _bt = "fo:font-weight=\"bold\" style:font-weight-asian=\"bold\" style:" ...
          +"font-weight-complex=\"bold\" ";
    }

    _it = "";
    if (_options.italic == "on")
    {
      _it = "fo:font-style=\"italic\" style:font-style-asian" ...
            +"=\"italic\" style:font-style-complex=\"italic\" ";
    }

    _ut = "";
    if (_options.underline == "on")
    {
      _ut = "style:text-underline=\"single\" style:text-underline-color=\"font-color\" ";
    }

    _tct = "";
    if (_tc != "000000")
    { _tct = "fo:color=\"#"+_tc+"\" "; }

    _tbt = "";
    if (_tb != "FFFFFF")
    { _tbt = "style:text-background-color=\"#"+_tb+"\" "; }

    _tst = "";
    if (_ts != _options.textsize_default)
    { _tst = "fo:font-size=\""+_ts+"pt\" "; }

    _fnt = "";
    if (exist(_fn))
    { _fnt = "style:font-name=\""+_fn+"\" "; }

    // does style starts with h/H?: if ($style =~ /^h/)
    if (tolower(substr(style,1))=="h")
    {
      _autostyles=[_autostyles, "<style:style style:name=\"" ...
          +_stylenum+"\" style:" ...
          +"family=\"paragraph\" style:parent-style-name=\""+style+"\">"];
      _autostyles=[_autostyles, "<style:properties "+_jt+"/>"];
      _autostyles=[_autostyles, "</style:style>"];
    } else {
      _autostyles=[_autostyles, "<style:style style:name=\""...
          +_stylenum+"\" style:family=\"paragraph\" style:parent-style-name=\"" ...
          +style+"\">"];
      _autostyles=[_autostyles, "<style:properties "+_tct+_tbt+_tst+_bt+_it+_ut+_jt+_fnt+"/>"];
      _autostyles=[_autostyles, "</style:style>"];
    }

    return _stylenum;
};

ooo.set = function (name, val, val2)
{
  // Set variables in the options hash
  if (!exist(name))
  {
    printf("ooo.set: requires two arguments\n");
    return "error";
  }

  // we use global variables to communicate some of the options
  if (name == "format" && strlen(val)>0)
  { _ooolib_FloatFormat = val; }

  if ((name=="title" || name=="author" || name=="subject" || ...
       name=="comments") && class(val)=="string")
  {
    _options.[name] = _oooCleanText( val );
  }
  if (name=="builddir" && class(val)=="string")
  { _builddir = val; }

  if (name=="keyword" && class(val)=="string")
  { _keywords = [_keywords, val[:]']; }

  if (strindex(name,"meta"))
  {
    id = substr(name,5);
    if (strindex(name, "name"))
    { _options.["info"+id+" name"] = _oooCleanText( val ); }
    if (strindex(name, "value"))
    { _options.["info"+id+" value"] = _oooCleanText( val ); }
  }

  if (name == "cell-loc")
  {
    x = nan();
    y = nan();
    if (class(val)=="string")
    {
      // val provides the writing cell
      xy = _ooosxcCellLocConv( val );
      x = xy[1];
      y = xy[2];
    }
    if (class(val)=="num")
    {
      // is val=[x,y] ?
      if(length(val)==2)
      {
        x=val[1];
        y=val[2];
      }
      // if val=x, then check for val2
      if (length(val)==1)
      {
        x=val;
        if (exist(val2))
        {
          if (class(val2)=="num")
          { y = val2; }
        } else {
          error("ooo.set: 'cell-loc' requires two coordinates");
        }
      }
    }

    if(isnan(x)||isnan(y))
    { error("ooo.set: 'cell-loc' requires two coordinates"); }

    _cellhash.x = x;
    _cellhash.y = y;
    _oooCellCheck();
  }

  if (name == "column-width")
  {
    if (!exist(val) || !exist(val2))
    { error("ooo.set: column-width requires 'column' and 'width' arguments\n"); }
    if (class(val)!="num" || class(val2)!="num")
    { error("ooo.set: column-width requires 'column' and 'width' arguments\n"); }

    // Set width of column
    col = val;
    wid = val2;
    _cellhash.["column "+text(col)+" width"] = wid;
  }
  if (name == "cell-left")
  {
    // Subtract 1 from x
    _cellhash.x = _cellhash.x - 1;
    _oooCellCheck();
  }
  if (name == "cell-right")
  {
    // Add 1 to x
    _cellhash.x = _cellhash.x + 1;
    _oooCellCheck();
  }
  if (name == "cell-up")
  {
    // Subtract 1 from y
    _cellhash.y = _cellhash.y - 1;
    _oooCellCheck();
  }
  if (name == "cell-down")
  {
    // Add 1 to y
    _cellhash.y = _cellhash.y + 1;
    _oooCellCheck();
  }
  if (name == "cell-auto")
  {
    if (!exist(val) || !exist(val2))
    { error("ooo.set: 'cell-auto' requires 'xauto' and 'yauto' arguments\n"); }
    if (class(val)!="num" || class(val2)!="num")
    { error("ooo.set: 'cell-auto' requires 'xauto' and 'yauto' arguments\n"); }

    // Set the auto-increment values.
    _cellhash.xauto = val;
    _cellhash.yauto = val2;
  }

  if (name == "justify")
  {
    if (!exist(val))
    { error("ooo.set: 'justify' requires string argument\n"); }
    if (class(val)!="string")
    { error("ooo.set: 'justify' requires string argument\n"); }

    // Set the justification of the text
    if (val=="right" || val=="left" || val=="center" || val=="block")
    { _options.justify = val; }
  }
  if (name == "bold")
  {
    if (!exist(val))
    { error("ooo.set: 'bold' requires string argument\n"); }
    if (class(val)!="string")
    { error("ooo.set: 'bold' requires string argument\n"); }

    // Adjust the bold properties
    if (val == "on" || val == "off")
    { _options.bold = val; }
  }
  if (name == "italic")
  {
    if (!exist(val))
    { error("ooo.set: 'italic' requires string argument\n"); }
    if (class(val)!="string")
    { error("ooo.set: 'italic' requires string argument\n"); }

    // Adjust the italic properties
    if (val == "on" || val == "off")
    { _options.italic = val; }
  }
  if (name == "underline")
  {
    if (!exist(val))
    { error("ooo.set: 'underline' requires string argument\n"); }
    if (class(val)!="string")
    { error("ooo.set: 'underline' requires string argument\n"); }

    // Adjust the underline properties
    if (val == "on" || val == "off")
    { _options.underline = val; }
  }

  if (name == "text-color")
  {
    if (!exist(val))
    { error("ooo.set: 'text-color' requires string argument\n"); }
    if (class(val)!="string")
    { error("ooo.set: 'text-color' requires string argument\n"); }

    // allow user shortcuts to predefined colors
    if (val == "default")
    {
      _options.textcolor = "000000";
      return "ok";
    }
    if (val == "red")
    {
      // red = FF0000
      _options.textcolor = "FF0000";
      return "ok";
    }
    if (val == "green")
    {
      // Green = 00FF00
      _options.textcolor = "00FF00";
      return "ok";
    }
    if (val == "blue")
    {
      // Blue == 0000FF
      _options.textcolor = "0000FF";
      return "ok";
    }
    if (val == "black")
    {
      // Defaults back to black
      _options.textcolor = "000000";
      return "ok";
    }
    if (val == "white")
    {
      // White = FFFFFF
      _options.textcolor = "FFFFFF";
      return "ok";
    }

    // is val a 6 digit valid hexadecimal string?
    if (length(val)==6)
    {
      mi = [48:57,65:70]; // 0,1..9,A,..F
      ia = ascii(toupper(val));
      s = 0;
      for (i in 1:6)
      { s = s + sum(ia[i]==mi); }
      if (s!=6)
      {
        printf("ooo.set: 'text-color' argument %s is not a valid color\n", val);
        return "error";
      }
    }
    _options.textcolor = toupper(val);
  }

  if (name == "text-bgcolor")
  {
    if (!exist(val))
    { error("ooo.set: 'text-bgcolor' requires string argument\n"); }
    if (class(val)!="string")
    { error("ooo.set: 'text-bgcolor' requires string argument\n"); }

    // allow user shortcuts to predefined colors
    if (val == "default")
    {
      _options.textbgcolor = "000000";
      return "ok";
    }
    if (val == "red")
    {
      // red = FF0000
      _options.textbgcolor = "FF0000";
      return "ok";
    }
    if (val == "green")
    {
      // Green = 00FF00
      _options.textbgcolor = "00FF00";
      return "ok";
    }
    if (val == "blue")
    {
      // Blue == 0000FF
      _options.textbgcolor = "0000FF";
      return "ok";
    }
    if (val == "black")
    {
      // Defaults back to black
      _options.textbgcolor = "000000";
      return "ok";
    }
    if (val == "white")
    {
      // White = FFFFFF
      _options.textbgcolor = "FFFFFF";
      return "ok";
    }

    // is val a 6 digit valid hexadecimal string?
    if (length(val)==6)
    {
      mi = [48:57,65:70]; // 0,1..9,A,..F
      ia = ascii(toupper(val));
      s = 0;
      for (i in 1:6)
      { s = s + sum(ia[i]==mi); }
      if (s!=6)
      {
        printf("ooo.set: 'text-bgcolor' argument %s is not a valid color\n", val);
        return "error";
      }
    }
    _options.textbgcolor = toupper(val);
    return "ok";
  }

  if (name == "text-size")
  {
    if (!exist(val))
    { error("ooo.set: 'text-size' requires numeric argument or \"default\"\n"); }
    if (class(val)!="num" || val=="default")
    { error("ooo.set: 'text-size' requires numeric argument or \"default\"\n"); }

    if (val == "default")
    {
      _options.textsize = _options.textsize_default;
      return "ok";
    }

    if (class(val)=="string")
    {val = text(val); }
    if (val >= strtod(_options.textsize_min) && val <= strtod(_options.textsize_max))
    {

      _options.textsize = text(val);
      return "ok";
    }

    return "error";
  }

  if (name == "text-font")
  {
    if (!exist(val))
    { error("ooo.set: 'text-font' requires string argument\n"); }
    if (class(val)!="string")
    { error("ooo.set: 'text-font' requires string argument\n"); }

    if (exist(_options.["font-decl "+val]))
    {
      _options.fontname = _options.["font-decl "+val];

    } else {

      if (val == "Arial")
      {
        _fontdecls = [_fontdecls, "<style:font-decl style:name=\"Arial\" fo:" ...
            +"font-family=\"Arial\" style:font-family-generic" ...
                +"=\"swiss\" style:font-pitch=\"variable\"/>"];
        _options.fontname = "Arial";
        _options.["font-decl "+val] = _options.fontname;
      }
      if (val == "Bitstream Vera Sans")
      {
        _fontdecls = [_fontdecls, "<style:font-decl style:name" ...
            +"=\"Bitstream Vera Sans\" fo:font-family" ...
                +"=\"\&apos;Bitstream Vera Sans\&apos;\" style:font-pitch" ...
                +"=\"variable\"/>"];
        _options.fontname = "Bitstream Vera Sans";
        _options.["font-decl "+val] = _options.fontname;
      }
      if (val == "Bitstream Vera Serif") {
        _fontdecls = [_fontdecls, "<style:font-decl style:name" ...
            +"=\"Bitstream Vera Serif\" fo:font-family" ...
                +"=\"\&apos;Bitstream Vera Serif\&apos;\" style:" ...
                +"font-family-generic=\"roman\" style:font-pitch=\"variable\"/>"];
        _options.fontname = "Bitstream Vera Serif";
        _options.["font-decl "+val] = _options.fontname;
      }
      if (val == "Bookman") {
        _fontdecls = [_fontdecls, "<style:font-decl style:name=\"Bookman\" fo:font-family" ...
            +"=\"Bookman\" style:font-pitch=\"variable\"/>"];
        _options.fontname = "Bookman";
        _options.["font-decl "+val] = _options.fontname;
      }
      if (val == "Courier") {
        _fontdecls = [_fontdecls, "<style:font-decl style:name=\"Courier\" fo:font-family" ...
            +"=\"Courier\" style:font-family-generic=\"modern\" style:" ...
                +"font-pitch=\"fixed\"/>"];
        _options.fontname = "Courier";
        _options.["font-decl "+val] = _options.fontname;
      }
      if (val == "Courier 10 Pitch") {
        _fontdecls = [_fontdecls, "<style:font-decl style:name=\"Courier 10 Pitch\" fo:" ...
            +"font-family=\"\&apos;Courier 10 Pitch\&apos;\" style:font-pitch" ...
                +"=\"fixed\"/>"];
        _options.fontname = "Courier 10 Pitch";
        _options.["font-decl "+val] = _options.fontname;
      }
      if (val == "Helvetica") {
        _fontdecls = [_fontdecls, "<style:font-decl style:name=\"Helvetica\" fo:font-family" ...
            +"=\"Helvetica\" style:font-family-generic=\"swiss\" style:font-pitch" ...
                +"=\"variable\"/>"];
        _options.fontname = "Helvetica";
        _options.["font-decl "+val] = _options.fontname;
      }
      if (val == "Lucidabright") {
        _fontdecls = [_fontdecls, "<style:font-decl style:name=\"Lucidabright\" fo:" ...
            +"font-family=\"Lucidabright\" style:font-pitch=\"variable\"/>"];
        _options.fontname = "Lucidabright";
        _options.["font-decl "+val] = _options.fontname;
      }
      if (val == "Lucidasans") {
        _fontdecls = [_fontdecls, "<style:font-decl style:name=\"Lucidasans\" fo:font-family" ...
            +"=\"Lucidasans\" style:font-pitch=\"variable\"/>"];
        _options.fontname = "Lucidasans";
        _options.["font-decl "+val] = _options.fontname;
      }
      if (val == "Lucida Sans Unicode") {
        _fontdecls = [_fontdecls, "<style:font-decl style:name=\"Lucida Sans Unicode\" fo:" ...
            +"font-family=\"\&apos;Lucida Sans Unicode\&apos;\" style:font-pitch=\"variable\"/>"];
        _options.fontname = "Lucida Sans Unicode";
        _options.["font-decl "+val] = _options.fontname;
      }
      if (val == "Lucidatypewriter") {
        _fontdecls = [_fontdecls, "<style:font-decl style:name=\"Lucidatypewriter\" fo:" ...
            +"font-family=\"Lucidatypewriter\" style:font-pitch=\"fixed\"/>"];
        _options.fontname = "Lucidatypewriter";
        _options.["font-decl "+val] = _options.fontname;
      }
      if (val == "Luxi Mono") {
        _fontdecls = [_fontdecls, "<style:font-decl style:name=\"Luxi Mono\" fo:font-family" ...
            +"=\"\&apos;Luxi Mono\&apos;\" style:font-pitch" ...
                +"=\"fixed\" style:font-charset=\"x-symbol\"/>"];
        _options.fontname = "Luxi Mono";
        _options.["font-decl "+val] = _options.fontname;
      }
      if (val == "Luxi Sans") {
        _fontdecls = [_fontdecls, "<style:font-decl style:name=\"Luxi Sans\" fo:" ...
            +"font-family=\"\&apos;Luxi Sans\&apos;\" style:font-pitch" ...
                +"=\"variable\" style:font-charset=\"x-symbol\"/>"];
        _options.fontname = "Luxi Sans";
        _options.["font-decl "+val] = _options.fontname;
      }
      if (val == "Luxi Serif") {
        _fontdecls = [_fontdecls,  "<style:font-decl style:name=\"Luxi Serif\" fo:" ...
            +"font-family=\"\&apos;Luxi Serif\&apos;\" style:font-pitch" ...
                +"=\"variable\" style:font-charset=\"x-symbol\"/>"];
        _options.fontname = "Luxi Serif";
        _options.["font-decl "+val] = _options.fontname;
      }
      if (val == "Symbol") {
        _fontdecls = [_fontdecls, "<style:font-decl style:name=\"Symbol\" fo:font-family" ...
            +"=\"Symbol\" style:font-pitch=\"variable\" style:font-charset" ...
                +"=\"x-symbol\"/>"];
        _options.fontname = "Symbol";
        _options.["font-decl "+val] = _options.fontname;
      }
      if (val == "Tahoma") {
        _fontdecls = [_fontdecls, "<style:font-decl style:name=\"Tahoma\" fo:" ...
            +"font-family=\"Tahoma\" style:font-pitch=\"variable\"/>"];
        _options.fontname = "Tahoma";
        _options.["font-decl "+val] = _options.fontname;
      }
      if (val == "Times") {
        _fontdecls = [_fontdecls, "<style:font-decl style:name=\"Times\" fo:" ...
            +"font-family=\"Times\" style:font-family-generic" ...
                +"=\"roman\" style:font-pitch=\"variable\"/>"];
        _options.fontname = "Times";
        _options.["font-decl "+val] = _options.fontname;
      }
      if (val == "Times New Roman") {
        _fontdecls = [_fontdecls, "<style:font-decl style:name=\"Times New Roman\" fo:" ...
            +"font-family=\"\&apos;Times New Roman\&apos;\" style:" ...
                +"font-family-generic=\"roman\" style:font-pitch=\"variable\"/>"];
        _options.fontname = "Times New Roman";
        _options.["font-decl "+val] = _options.fontname;
      }
      if (val == "Utopia") {
        _fontdecls = [_fontdecls, "<style:font-decl style:name=\"Utopia\" fo:" ...
            +"font-family=\"Utopia\" style:font-family-generic" ...
                +"=\"roman\" style:font-pitch=\"variable\"/>"];
        _options.fontname = "Utopia";
        _options.["font-decl "+val] = _options.fontname;
      }
      if (val == "Zapf Chancery") {
        _fontdecls = [_fontdecls, "<style:font-decl style:name=\"Zapf Chancery\" fo:" ...
            +"font-family=\"\&apos;Zapf Chancery\&apos;\" style:font-pitch=\"variable\"/>"];
        _options.fontname = "Zapf Chancery";
        _options.["font-decl "+val] = _options.fontname;
      }
      if (val == "Zapf Dingbats") {
        _fontdecls = [_fontdecls, "<style:font-decl style:name=\"Zapf Dingbats\" fo:" ...
            +"font-family=\"\&apos;Zapf Dingbats\&apos;\" style:font-pitch" ...
                +"=\"variable\" style:font-charset=\"x-symbol\"/>"];
        _options.fontname = "Zapf Dingbats";
        _options.["font-decl "+val] = _options.fontname;
      }
    }
  }
  return "ok";
};


_oooCellUpdate = function()
{
  // Updates max values and current values for x and y
  if (_cellhash.x > _cellhash.xmax)
  { _cellhash.xmax = _cellhash.x; }
  if (_cellhash.y > _cellhash.ymax)
  { _cellhash.ymax = _cellhash.y; }

  _cellhash.x = _cellhash.x + _cellhash.xauto;
  _cellhash.y = _cellhash.y + _cellhash.yauto;

  _oooCellCheck();
};

_oooCellCheck = function ()
{
  // Tests to make sure the x and y values are legal
  if (_cellhash.x < _minx)
  { _cellhash.x = _minx; }
  if (_cellhash.y < _miny)
  { _cellhash.y = _miny; }
  if (_cellhash.x > _maxx)
  { _cellhash.x = _maxx; }
  if (_cellhash.y > _maxy)
  { _cellhash.y = _maxy; }
};

_oooCleanText = function (val)
{
  if (class(val)=="num")
  {
    if (isnan(val))
    {
      return "";
    } else {
      return text(val,_ooolib_FloatFormat);
    }
  }

  // & -> &amp;
  val = gsub("&amp;","&",val).string;
  // ' -> &apos;
  val = gsub("&apos;","'",val).string;
  // \" -> &quot;
  val = gsub("&quot;","\"",val).string;
  // < -> &lt;
  val = gsub("&lt;","<",val).string;
  // > -> &gt;
  val = gsub("&gt;",">",val).string;
  // \t -> <text:tab-stop/>
  val = gsub("<text:tab-stop/>","\t",val).string;

  // i gotta figure out what to do with this?
  // from  http://perl-xml.sourceforge.net/faq/#encoding_conversion
  // $text =~ s/([^\x20-\x7F])/'&#' . ord($1) . ';'/gse;

  return val;
};

ooo.generate = function(filename_final)
{
  // Create the document and return a filename
  mykind   = _options.kind;

  // allow user to specify the builddir
  if (_builddir == "")
  { _builddir = _BUILDDIR_; }

  if (!exist(mykind))
  {
    printf("_oooGenerate: Document type not selected\n");
    return "error";
  }

  currdir = pwd();
  tempfname = timestamp();

  tempfnext=".sxw";
  if (mykind=="sxc" || mykind=="xls")
  { tempfnext  = ".sxc"; }

  // Create the files
  if (_ooodebug)
  { printf("Writing XML files in \"%s\"\n", _builddir); }
  _oooWriteMimetype (mykind);
  _oooWriteContent  (mykind);
  _oooWriteStyles   (mykind);
  _oooWriteMeta     (mykind);
  _oooWriteSettings (mykind);
  _oooWriteManifest (mykind);

  // change local sheel directory
  cd(_builddir);

  if (_ooodebug)
  { printf("Building $type file from XML\n"); }

  system("zip -r '" + tempfname + tempfnext + "' mimetype > /dev/null 2> /dev/null");
  system("zip -r '" + tempfname + tempfnext + "' content.xml > /dev/null 2> /dev/null");
  system("zip -r '" + tempfname + tempfnext + "' styles.xml > /dev/null 2> /dev/null");
  system("zip -r '" + tempfname + tempfnext + "' meta.xml > /dev/null 2> /dev/null");
  system("zip -r '" + tempfname + tempfnext + "' settings.xml > /dev/null 2> /dev/null");
  system("zip -r '" + tempfname + tempfnext + "' META-INF/manifest.xml > /dev/null 2> /dev/null");
  if (_ooodebug)
  {
    system("cp " + tempfname + tempfnext + " ../");
  }
  else
  {
    system("mv " + tempfname + tempfnext + " ../");
  }
  cd("../");

  // rm files that were created in the process
  if (_ooodebug)
  {
    printf("Skipping Clean-up\n");
  }
  else
  {
    system("rm -rf " + _builddir);
  }

  if (strindex(mykind,"xls"))
  {
    _cmd = _OOOEXEC_ + " --headless --convert-to xls:\"MS Excel 97\" " + tempfname + tempfnext;
    _cmd?
    if (_ooodebug)
    { printf(_cmd + "\n"); }
    system(_cmd); // convert it to excel
    _j=0;
    while(!isfile(tempfname + ".xls") && _j<20)
    {
      sleep(0.1);
      printf("*");
      _j++;
    }
    if (_ooodebug)
    {
      system("cp " + tempfname + ".xls"    + " " + filename_final + ".xls");
    }
    else
    {
      system("mv " + tempfname + ".xls"    + " " + filename_final + ".xls");
    }
  }

  if (_ooodebug)
  {
    system("cp " + tempfname + tempfnext + " " + filename_final);
  }
  else
  {
    system("mv " + tempfname + tempfnext + " " + filename_final);
  }


  // Return the filename
  return 0;
};


_ooosxcCellLocConv = function(cellid)
{
  // Convert A1 to 1,1
  cellid = toupper(cellid);
  lc = strlen(cellid);
  for (i in 1:lc)
  {
    // pick i-th character from cellid
    ic = ascii(substr(cellid, i));
    // is it a character?
    if (ic >= 65 && ic <= 90)
    { continue; }
    // is it a digit ?
    if (ic >= 48 && ic <= 57)
    { break; }
    // if neither, stop
    printf("_ooosxcCellLocConv: %s is not proper cell id\n", cellid);
    error("_ooosxcCellLocConv: improper 'cellid'. Check your code");
  }
  // x:
  xid = substr(cellid,1:(i-1));
  pw  = ones(xid);
  for (j in (i-1:1:-1))
  { pw[i-j] = 26^(j-1); }
  x   = sum((ascii(xid)-64).*pw);

  // y:
  yid = substr(cellid,i:lc);
  y   = strtod(yid);

  // report
  return [x, y];
};


_ooosxcCellLeft = function (cellx)
{
  // Subtract 1 from cellx
  x = cellx - 1;
  if (x <= 0)
  { x = 1; }
  return x;
};

_ooosxcCellRight = function(cellx)
{
  // Add 1 to cellx
  return (cellx + 1);
};

_ooosxcCellUp = function (celly)
{
  // Subtract 1 from celly
  x = celly - 1;
  if (x <= 0)
  { x = 1; }
  return x;
};

_ooosxcCellDown = function (celly)
{
  // Add 1 to celly
  return (celly + 1);
};

//
_oooTimeStamp = function()
{
  return time2dstr(seconds(),"%c");
};


_oooWriteMimetype = function( mykind )
{
  // This writes the mimetype file
  if (!isdir(_builddir))
  { mkdir(_builddir); }

  mimefile = _builddir + "/mimetype";
  if (_ooodebug)
  { printf ("_oooWriteMimetype: using file %s\n", mimefile); }

  // Open file for writing
  open(mimefile,"w");

  if (mykind == "sxw")
  {
    fprintf(mimefile, "application/vnd.sun.xml.writer");
  } else { if (mykind == "sxc")
  {
    fprintf(mimefile, "application/vnd.sun.xml.calc");
  }}
  close(mimefile);
};

//
// local function to write the metadata for the ooo document
// in a subdirectory of a working directory
//
_oooWriteMeta = function(mykind)
{
  if (!isdir(_builddir))
  { mkdir(_builddir); }

  if (!exist(mykind))
  { error("_oooWriteMeta: requires file type"); }

  metafile = _builddir + "/meta.xml";
  if (_ooodebug)
  { printf ("oooWriteMeta: using metafile %s\n", metafile); }

  metadate = _oooTimeStamp();

  open(metafile,"w");
  // Encoding Information
  fprintf(metafile, "<?xml version=\"1.0\" encoding=\"UTF-8\"?>");
  // DOCTYPE line
  fprintf(metafile, "<!DOCTYPE office:document-meta PUBLIC" ...
      + " \"-//OpenOffice.org//DTD OfficeDocument 1.0//EN\" \"office.dtd\">");
  // Information about decoding the document and XML
  fprintf(metafile, "<office:document-meta xmlns:office=\"http" ...
      + "://openoffice.org/2000/office\" xmlns:xlink=\"" ...
      +"http://www.w3.org/1999/xlink\" xmlns:" ...
      + "dc=\"http://purl.org/dc/elements/1.1/\" xmlns:meta=\"" ...
      + "http://openoffice.org/2000/meta\" office:version=\"1.0\">");
  // Beginning of the meta information
  fprintf(metafile, "<office:meta>");
  fprintf(metafile, "<meta:generator> %s v.%s Generator</meta:generator>", ...
         _ooolibrary, _oooversion);

  // Load information from _options, but prepare for missing data
  mytitle   = "rlabplus";
  mycomment = "OOO file generated by libooo/rlabplus";
  mysubject = "No particular subject";
  myauthor  = "Anonymous";

  // read header info from _options list
  if (exist(_options.title))
  { mytitle = _options.title; }
  if (exist(_options.comment))
  { mycomment = _options.comment; }
  if (exist(_options.subject))
  { mysubject = _options.subject; }
  if (exist(_options.author))
  { myauthor = _options.author; }

  // write header information
  fprintf(metafile, "<dc:title>%s</dc:title>", mytitle);
  fprintf(metafile, "<dc:description>%s</dc:description>", mycomment);
  fprintf(metafile, "<dc:subject>%s</dc:subject>", mysubject);
  fprintf(metafile, "<meta:initial-creator>%s</meta:initial-creator>", myauthor);
  fprintf(metafile, "<meta:creation-date>%s</meta:creation-date>", metadate);
  fprintf(metafile, "<dc:creator>%s</dc:creator>", myauthor);
  fprintf(metafile, "<dc:date>%s</dc:date>", metadate);
  fprintf(metafile, "<meta:keywords>");
  if (length(_keywords) > 0)
  {
    for (k in _keywords)
    { fprintf(metafile,  "<meta:keyword>%s</meta:keyword>", k); }
  }
  fprintf(metafile, "</meta:keywords>");

  // Document Language
  fprintf(metafile, "<dc:language>%s</dc:language>", _ooolanguage);
  // This is the document version number - pull an uniform random integer
  fprintf(metafile, "<meta:editing-cycles>%u</meta:editing-cycles>", irand());
  // Time the document has been worked on - pull an uniform random real
  fprintf(metafile,  "<meta:editing-duration>PT0S</meta:editing-duration>");

  // This is user defined variables.
  for (i in [1:4])
  {
    if (exist(_options.["info"+text(i)+" value"]))
    {
      name = _options.["info"+text(i)+" name"];
      val  = _options.["info"+text(i)+" value"];
      fprintf(metafile, "<meta:user-defined" ...
          + " meta:name=\"%s\">%s</meta:user-defined>", name, val);
    } else {
      name = _options.["info"+text(i)+" name"];
      fprintf(metafile, "<meta:user-defined meta:name=\"%s\"/>", name);
    }
  } // next i

  // Statistical Data
  if (mykind == "sxw")
  {
    // oowriter -  text
    fprintf(metafile,  "<meta:document-statistic meta:table-count=\"0\"" ...
        + " meta:image-count=\"0\" meta:object-count=\"0\" meta:page-count=\"1\"" ...
        + " meta:paragraph-count=\"1\" meta:word-count=\"0\"" ...
        + " meta:character-count=\"0\"/>");
  } else { if (mykind == "sxc")
  {
    // oocalc - spreadsheet
    cellcount = _cellhash.ymax * _cellhash.xmax;
    fprintf(metafile,  "<meta:document-statistic meta:table-count=\"1\" " ...
        + "meta:cell-count=\"%g\"/>", cellcount);
  }}

  fprintf(metafile, "</office:meta>");
  fprintf(metafile, "</office:document-meta>");

  close(metafile);
};

//
// Writes the META-INF/manifest.xml file
//
_oooWriteManifest = function(mykind)
{
  if(!exist(mykind))
  { error("_oooWriteManifest: missing file type\n"); }

  if (!isdir(_builddir))
  { mkdir(_builddir); }
  if (!isdir(_builddir + "/META-INF"))
  { mkdir(_builddir + "/META-INF"); }
  manfile = _builddir + "/META-INF/manifest.xml";

  if (_ooodebug)
  { printf ("_oooWriteManifest: using manifest file %s\n", manfile); }

  open(manfile,"w");

  // encoding Information
  fprintf(manfile, "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
  // DOCTYPE line
  fprintf(manfile, "<!DOCTYPE manifest:manifest PUBLIC \"-//OpenOffice.org" ...
      + "//DTD Manifest 1.0//EN\" \"Manifest.dtd\">\n");
  // reading information
  fprintf(manfile, "<manifest:manifest xmlns:manifest=\"" ...
      + "http://openoffice.org/2001/manifest\">\n");
  // Mime type
  if (mykind == "sxw")
  {
    // oowriter
    fprintf(manfile, "<manifest:file-entry manifest:media-type=\"" ...
        + "application/vnd.sun.xml.writer\" manifest:full-path=\"/\"/>\n");
  } else { if (mykind == "sxc")
  {
    // oocalc
    fprintf(manfile, " <manifest:file-entry manifest:media-type=\"" ...
        + "application/vnd.sun.xml.calc\" manifest:full-path=\"/\"/>\n");
  }}

  // For pictures
  fprintf(manfile, " <manifest:file-entry manifest:media-type=\"\"" ...
      + " manifest:full-path=\"Pictures/\"/>\n");
  // Contents File
  fprintf(manfile, " <manifest:file-entry manifest:media-type=\"text/xml\"" ...
      + " manifest:full-path=\"content.xml\"/>\n");

  // Styles File
  fprintf(manfile, " <manifest:file-entry manifest:media-type=\"text/xml\"" ...
      + " manifest:full-path=\"styles.xml\"/>\n");

  // Meta File
  fprintf(manfile, " <manifest:file-entry manifest:media-type=\"text/xml\"" ...
      + " manifest:full-path=\"meta.xml\"/>\n");

  // Settings File
  fprintf(manfile, " <manifest:file-entry manifest:media-type=\"text/xml\"" ...
      + " manifest:full-path=\"settings.xml\"/>\n");

  // End of file
  fprintf(manfile, "</manifest:manifest>\n");
  close(manfile);
};


_oooWriteStyles = function(mykind)
{
  if(!exist(mykind))
  { error("_oooWriteStyles: missing file type\n"); }

  if (!isdir(_builddir))
  { mkdir(_builddir); }

  stylefile = _builddir + "/styles.xml";
  if (_ooodebug)
  { printf ("_oooWriteStyles: using style file %s\n", stylefile); }

  open(stylefile,"w");

  if (mykind == "sxw")
  {
    // oo writer
    fprintf(stylefile,  "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
    fprintf(stylefile,  "<!DOCTYPE office:document-styles PUBLIC \"" ...
        + "-//OpenOffice.org//DTD OfficeDocument 1.0//EN\" \"office.dtd\">");
    fprintf(stylefile,  "<office:document-styles xmlns:office=\"" ...
        + "http://openoffice.org/2000/office\" xmlns:style=\"" ...
            + "http://openoffice.org/2000/style\" xmlns:text=\"" ...
            + "http://openoffice.org/2000/text\" xmlns:table=\"" ...
            + "http://openoffice.org/2000/table\" xmlns:draw=\"" ...
            + "http://openoffice.org/2000/drawing\" xmlns:fo=\"" ...
            + "http://www.w3.org/1999/XSL/Format\" xmlns:xlink=\"" ...
            + "http://www.w3.org/1999/xlink\" xmlns:number=\"" ...
            + "http://openoffice.org/2000/datastyle\" xmlns:svg=\"" ...
            + "http://www.w3.org/2000/svg\" xmlns:chart=\"" ...
            + "http://openoffice.org/2000/chart\" xmlns:dr3d=\"" ...
            + "http://openoffice.org/2000/dr3d\" xmlns:math=\"" ...
            + "http://www.w3.org/1998/Math/MathML\" xmlns:form=\"" ...
            + "http://openoffice.org/2000/form\" xmlns:script=\"" ...
            + "http://openoffice.org/2000/script\" office:version=\"1.0\">");

    if (length(_fontdecls)>0)
    {
      fprintf(stylefile, "<office:font-decls>");
      for (i in _fontdecls)
      { fprintf(stylefile, "%s", i); }
      fprintf(stylefile, "</office:font-decls>");
    }
    fprintf(stylefile,  "<office:styles>");
    fprintf(stylefile,  "<style:default-style style:family=\"graphics\">");
    fprintf(stylefile,  "<style:properties draw:" ...
        + "start-line-spacing-horizontal=\"0.283cm\" draw:" ...
            + "start-line-spacing-vertical=\"0.283cm\" draw:" ...
            + "end-line-spacing-horizontal=\"0.283cm\" draw:" ...
            + "end-line-spacing-vertical=\"0.283cm\" style:" ...
            + "use-window-font-color=\"true\" style:" ...
            + "font-name=\"Times New Roman\" fo:" ...
            + "font-size=\"12pt\" fo:language=\"en\" fo:" ...
            + "country=\"US\" style:font-name-asian=\"Times\" style:" ...
            + "font-size-asian=\"12pt\" style:" ...
            + "language-asian=\"ja\" style:country-asian=\"JP\" style:" ...
            + "font-name-complex=\"Times\" style:font-size-complex=\"12pt\" style:" ...
            + "language-complex=\"none\" style:country-complex=\"none\" style:" ...
            + "text-autospace=\"ideograph-alpha\" style:line-break=\"strict\" style:" ...
            + "writing-mode=\"lr-tb\">");
    fprintf(stylefile,  "<style:tab-stops/>");
    fprintf(stylefile,  "</style:properties>");
    fprintf(stylefile,  "</style:default-style>");
    fprintf(stylefile,  "<style:default-style style:family=\"paragraph\">");
    fprintf(stylefile,  "<style:properties style:use-window-font-color=" ...
        + "\"true\" style:font-name=\"Times New Roman\" fo:" ...
            + "font-size=\"12pt\" fo:language=\"en\" fo:" ...
            + "country=\"US\" style:font-name-asian=\"Times\" style:" ...
            + "font-size-asian=\"12pt\" style:" ...
            + "language-asian=\"ja\" style:country-asian=\"JP\" style:" ...
            + "font-name-complex=\"Times\" style:font-size-complex=\"12pt\" style:" ...
            + "language-complex=\"none\" style:country-complex=\"none\" fo:" ...
            + "hyphenate=\"false\" fo:hyphenation-remain-char-count=\"2\" fo:" ...
            + "hyphenation-push-char-count=\"2\" fo:" ...
            + "hyphenation-ladder-count=\"no-limit\" style:" ...
            + "text-autospace=\"ideograph-alpha\" style:" ...
            + "punctuation-wrap=\"hanging\" style:" ...
            + "line-break=\"strict\" style:" ...
            + "tab-stop-distance=\"1.251cm\" style:writing-mode=\"page\"/>");
    fprintf(stylefile,  "</style:default-style>");
    fprintf(stylefile,  "<style:style style:name=\"Standard\" style:" ...
        + "family=\"paragraph\" style:class=\"text\"/>");
    fprintf(stylefile,  "<style:style style:name=\"Text body\" style:" ...
        + "family=\"paragraph\" style:parent-style-name=\"Standard\" style:" ...
            + "class=\"text\">");
    fprintf(stylefile,  "<style:properties fo:margin-top=\"0cm\" fo:" ...
        + "margin-bottom=\"0.212cm\"/>");
    fprintf(stylefile,  "</style:style>");
    fprintf(stylefile,  "<style:style style:name=\"List\" style:" ...
        + "family=\"paragraph\" style:parent-style-name=\"Text body\" style:" ...
            + "class=\"list\">");
    fprintf(stylefile,  "<style:properties style:font-name-asian=\"Times\"/>");
    fprintf(stylefile,  "</style:style>");
    fprintf(stylefile,  "<style:style style:name=\"Caption\" style:" ...
        + "family=\"paragraph\" style:parent-style-name=\"Standard\" style:" ...
            + "class=\"extra\">");
    fprintf(stylefile,  "<style:properties fo:margin-top=\"0.212cm\" fo:" ...
        + "margin-bottom=\"0.212cm\" fo:font-size=\"10pt\" fo:" ...
            + "font-style=\"italic\" style:" ...
            + "font-name-asian=\"Times\" style:" ...
            + "font-size-asian=\"10pt\" style:font-style-asian=\"italic\" style:" ...
            + "font-size-complex=\"10pt\" style:font-style-complex=\"italic\" text:" ...
            + "number-lines=\"false\" text:line-number=\"0\"/>");
    fprintf(stylefile,  "</style:style>");
    fprintf(stylefile,  "<style:style style:name=\"Index\" style:" ...
        + "family=\"paragraph\" style:parent-style-name=\"Standard\" style:" ...
            + "class=\"index\">");
    fprintf(stylefile,  "<style:properties style:" ...
        + "font-name-asian=\"Times\" text:number-lines=\"false\" text:" ...
            + "line-number=\"0\"/>");
    fprintf(stylefile,  "</style:style>");
    fprintf(stylefile,  "<text:outline-style>");
    fprintf(stylefile,  "<text:outline-level-style text:level=\"1\" style:num-format=\"\"/>");
    fprintf(stylefile,  "<text:outline-level-style text:level=\"2\" style:num-format=\"\"/>");
    fprintf(stylefile,  "<text:outline-level-style text:level=\"3\" style:num-format=\"\"/>");
    fprintf(stylefile,  "<text:outline-level-style text:level=\"4\" style:num-format=\"\"/>");
    fprintf(stylefile,  "<text:outline-level-style text:level=\"5\" style:num-format=\"\"/>");
    fprintf(stylefile,  "<text:outline-level-style text:level=\"6\" style:num-format=\"\"/>");
    fprintf(stylefile,  "<text:outline-level-style text:level=\"7\" style:num-format=\"\"/>");
    fprintf(stylefile,  "<text:outline-level-style text:level=\"8\" style:num-format=\"\"/>");
    fprintf(stylefile,  "<text:outline-level-style text:level=\"9\" style:num-format=\"\"/>");
    fprintf(stylefile,  "<text:outline-level-style text:" ...
        + "level=\"10\" style:num-format=\"\"/>");
    fprintf(stylefile,  "</text:outline-style>");
    fprintf(stylefile,  "<text:footnotes-configuration style:" ...
        + "num-format=\"1\" text:start-value=\"0\" text:" ...
            + "footnotes-position=\"page\" text:start-numbering-at=\"document\"/>");
    fprintf(stylefile,  "<text:endnotes-configuration style:" ...
        + "num-format=\"i\" text:start-value=\"0\"/>");
    fprintf(stylefile,  "<text:linenumbering-configuration text:" ...
        + "number-lines=\"false\" text:offset=\"0.499cm\" style:" ...
            + "num-format=\"1\" text:number-position=\"left\" text:" ...
            + "increment=\"5\"/>");
    fprintf(stylefile,  "</office:styles>");
    fprintf(stylefile,  "<office:automatic-styles>");
    fprintf(stylefile,  "<style:page-master style:name=\"pm1\">");
    fprintf(stylefile,  "<style:properties fo:page-width=\"20.999cm\" fo:" ...
        + "page-height=\"29.699cm\" style:num-format=\"1\" style:" ...
            + "print-orientation=\"portrait\" fo:margin-top=\"2cm\" fo:" ...
            + "margin-bottom=\"2cm\" fo:margin-left=\"2cm\" fo:" ...
            + "margin-right=\"2cm\" style:writing-mode=\"lr-tb\" style:" ...
            + "footnote-max-height=\"0cm\">");
    fprintf(stylefile,  "<style:footnote-sep style:width=\"0.018cm\" style:" ...
        + "distance-before-sep=\"0.101cm\" style:distance-after-sep=\"0.101cm\" style:" ...
            + "adjustment=\"left\" style:rel-width=\"25%%\" style:color=\"#000000\"/>");
    fprintf(stylefile,  "</style:properties>");
    fprintf(stylefile,  "<style:header-style/>");
    fprintf(stylefile,  "<style:footer-style/>");
    fprintf(stylefile,  "</style:page-master>");
    fprintf(stylefile,  "</office:automatic-styles>");
    fprintf(stylefile,  "<office:master-styles>");
    fprintf(stylefile,  "<style:master-page style:" ...
        + "name=\"Standard\" style:page-master-name=\"pm1\"/>");
    fprintf(stylefile,  "</office:master-styles>");
    fprintf(stylefile,  "</office:document-styles>");

  } else { if (mykind == "sxc")
  {
    // oocalc
    fprintf(stylefile,  "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
    fprintf(stylefile,  "<!DOCTYPE office:" ...
        + "document-styles PUBLIC \"-//OpenOffice." ...
            + "org//DTD OfficeDocument 1.0//EN\" \"office.dtd\">");
    fprintf(stylefile,  "<office:document-styles xmlns:" ...
        + "office=\"http://openoffice.org/2000/office\" xmlns:" ...
            + "style=\"http://openoffice.org/2000/style\" xmlns:" ...
            + "text=\"http://openoffice.org/2000/text\" xmlns:" ...
            + "table=\"http://openoffice.org/2000/table\" xmlns:" ...
            + "draw=\"http://openoffice.org/2000/drawing\" xmlns:" ...
            + "fo=\"http://www.w3.org/1999/XSL/Format\" xmlns:" ...
            + "xlink=\"http://www.w3.org/1999/xlink\" xmlns:" ...
            + "number=\"http://openoffice.org/2000/datastyle\" xmlns:" ...
            + "svg=\"http://www.w3.org/2000/svg\" xmlns:" ...
            + "chart=\"http://openoffice.org/2000/chart\" xmlns:" ...
            + "dr3d=\"http://openoffice.org/2000/dr3d\" xmlns:" ...
            + "math=\"http://www.w3.org/1998/Math/MathML\" xmlns:" ...
            + "form=\"http://openoffice.org/2000/form\" xmlns:" ...
            + "script=\"http://openoffice.org/2000/script\" office:version=\"1.0\">");

    if (length(_fontdecls)>0)
    {
      fprintf(stylefile, "<office:font-decls>");
      for (i in _fontdecls)
      {  fprintf(stylefile, "%s", i); }
      fprintf(stylefile, "</office:font-decls>");
    }

    fprintf(stylefile,  "<office:styles>");
    fprintf(stylefile,  "<style:default-style style:" ...
        + "family=\"table-cell\">");
    fprintf(stylefile,  "<style:properties style:" ...
        + "decimal-places=\"2\" style:font-name=\"Times\" fo:language=\"en\" fo:" ...
            + "country=\"US\" style:font-name-asian=\"Lucida Sans Unicode\" style:" ...
            + "language-asian=\"ja\" style:country-asian=\"JP\" style:" ...
            + "font-name-complex=\"Times\" style:language-complex=\"none\" style:" ...
            + "country-complex=\"none\" style:tab-stop-distance=\"0.5inch\"/>");
    fprintf(stylefile,  "</style:default-style>");
    fprintf(stylefile,  "<number:number-style style:name=\"N0\" style:" ...
        + "family=\"data-style\">");
    fprintf(stylefile,  "<number:number number:min-integer-digits=\"1\"/>");
    fprintf(stylefile,  "</number:number-style>");
    fprintf(stylefile,  "<number:currency-style style:" ...
        + "name=\"N104P0\" style:family=\"data-style\" style:volatile=\"true\">");
    fprintf(stylefile,  "<number:currency-symbol number:" ...
        + "language=\"en\" number:country=\"US\">\$</number:currency-symbol>");
    fprintf(stylefile,  "<number:number number:" ...
        + "decimal-places=\"2\" number:min-integer-digits=\"1\" number:" ...
            + "grouping=\"true\"/>");
    fprintf(stylefile,  "</number:currency-style>");
    fprintf(stylefile,  "<number:currency-style style:name=\"N104\" style:" ...
        + "family=\"data-style\">");
    fprintf(stylefile,  "<style:properties fo:color=\"#ff0000\"/>");
    fprintf(stylefile,  "<number:text>-</number:text>");
    fprintf(stylefile,  "<number:currency-symbol number:" ...
        + "language=\"en\" number:country=\"US\">\$</number:currency-symbol>");
    fprintf(stylefile,  "<number:number number:decimal-places=\"2\" number:" ...
        + "min-integer-digits=\"1\" number:grouping=\"true\"/>");
    fprintf(stylefile,  "<style:map style:" ...
        + "condition=\"value()&gt;=0\" style:apply-style-name=\"N104P0\"/>");
    fprintf(stylefile,  "</number:currency-style>");
    fprintf(stylefile,  "<style:style style:name=\"Default\" style:family=\"table-cell\">");
    fprintf(stylefile,  "<style:properties style:font-name-asian=\"Times\"/>");
    fprintf(stylefile,  "</style:style>");
    fprintf(stylefile,  "<style:style style:name=\"Result\" style:" ...
        + "family=\"table-cell\" style:parent-style-name=\"Default\">");
    fprintf(stylefile,  "<style:properties fo:font-style=\"italic\" style:" ...
        + "text-underline=\"single\" style:text-underline-color=\"font-color\" fo:" ...
            + "font-weight=\"bold\"/>");
    fprintf(stylefile,  "</style:style>");
    fprintf(stylefile,  "<style:style style:name=\"Result2\" style:" ...
        + "family=\"table-cell\" style:parent-style-name=\"Result\" style:" ...
            + "data-style-name=\"N104\"/>");
    fprintf(stylefile,  "<style:style style:name=\"Heading\" style:" ...
        + "family=\"table-cell\" style:parent-style-name=\"Default\">");
    fprintf(stylefile,  "<style:properties fo:text-align=\"center\" style:" ...
        + "text-align-source=\"fix\" fo:font-size=\"16pt\" fo:font-style=\"italic\" fo:" ...
            + "font-weight=\"bold\"/>");
    fprintf(stylefile,  "</style:style>");
    fprintf(stylefile,  "<style:style style:name=\"Heading1\" style:" ...
        + "family=\"table-cell\" style:parent-style-name=\"Heading\">");
    fprintf(stylefile,  "<style:properties fo:direction=\"ltr\" style:" ...
        + "rotation-angle=\"90\"/>");
    fprintf(stylefile,  "</style:style>");
    fprintf(stylefile,  "</office:styles>");
    fprintf(stylefile,  "<office:automatic-styles>");
    fprintf(stylefile,  "<style:page-master style:name=\"pm1\">");
    fprintf(stylefile,  "<style:properties style:writing-mode=\"lr-tb\"/>");
    fprintf(stylefile,  "<style:header-style>");
    fprintf(stylefile,  "<style:properties fo:" ...
        + "min-height=\"0.2957inch\" fo:margin-left=\"0inch\" fo:" ...
            + "margin-right=\"0inch\" fo:margin-bottom=\"0.0984inch\"/>");
    fprintf(stylefile,  "</style:header-style>");
    fprintf(stylefile,  "<style:footer-style>");
    fprintf(stylefile,  "<style:properties fo:" ...
        + "min-height=\"0.2957inch\" fo:margin-left=\"0inch\" fo:" ...
            + "margin-right=\"0inch\" fo:margin-top=\"0.0984inch\"/>");
    fprintf(stylefile,  "</style:footer-style>");
    fprintf(stylefile,  "</style:page-master>");
    fprintf(stylefile,  "<style:page-master style:name=\"pm2\">");
    fprintf(stylefile,  "<style:properties style:writing-mode=\"lr-tb\"/>");
    fprintf(stylefile,  "<style:header-style>");
    fprintf(stylefile,  "<style:properties fo:min-height=\"0.2957inch\" fo:" ...
        + "margin-left=\"0inch\" fo:margin-right=\"0inch\" fo:" ...
            + "margin-bottom=\"0.0984inch\" fo:border=\"0.0346inch solid #000000\" fo:" ...
            + "padding=\"0.0071inch\" fo:background-color=\"#c0c0c0\">");
    fprintf(stylefile,  "<style:background-image/>");
    fprintf(stylefile,  "</style:properties>");
    fprintf(stylefile,  "</style:header-style>");
    fprintf(stylefile,  "<style:footer-style>");
    fprintf(stylefile,  "<style:properties fo:min-height=\"0.2957inch\" fo:" ...
        + "margin-left=\"0inch\" fo:margin-right=\"0inch\" fo:" ...
            + "margin-top=\"0.0984inch\" fo:border=\"0.0346inch solid #000000\" fo:" ...
            + "padding=\"0.0071inch\" fo:background-color=\"#c0c0c0\">");
    fprintf(stylefile,  "<style:background-image/>");
    fprintf(stylefile,  "</style:properties>");
    fprintf(stylefile,  "</style:footer-style>");
    fprintf(stylefile,  "</style:page-master>");
    fprintf(stylefile,  "</office:automatic-styles>");
    fprintf(stylefile,  "<office:master-styles>");
    fprintf(stylefile,  "<style:master-page style:" ...
        + "name=\"Default\" style:page-master-name=\"pm1\">");
    fprintf(stylefile,  "<style:header>");
    fprintf(stylefile,  "<text:p>");
    fprintf(stylefile,  "<text:sheet-name>???</text:sheet-name>");
    fprintf(stylefile,  "</text:p>");
    fprintf(stylefile,  "</style:header>");
    fprintf(stylefile,  "<style:header-left style:display=\"false\"/>");
    fprintf(stylefile,  "<style:footer>");
    fprintf(stylefile,  "<text:p>Page <text:page-number>1</text:page-number>");
    fprintf(stylefile,  "</text:p>");
    fprintf(stylefile,  "</style:footer>");
    fprintf(stylefile,  "<style:footer-left style:display=\"false\"/>");
    fprintf(stylefile,  "</style:master-page>");
    fprintf(stylefile,  "<style:master-page style:" ...
        + "name=\"Report\" style:page-master-name=\"pm2\">");
    fprintf(stylefile,  "<style:header>");
    fprintf(stylefile,  "<style:region-left>");
    fprintf(stylefile,  "<text:p>");
    fprintf(stylefile,  "<text:sheet-name>???</text:sheet-name> (<text:" ...
        + "title>???</text:title>)</text:p>");
    fprintf(stylefile,  "</style:region-left>");
    fprintf(stylefile,  "<style:region-right>");
    fprintf(stylefile,  "<text:p>");

    // Date
    fprintf(stylefile,  "<text:date style:data-style-name=\"N2\" text:" ...
        + "date-value=\"%s-%s-%s\">%s/%s/%s</text:date>, <text:" ...
        + "time>%s:%s:%s</text:time>", ...
            time2dstr(seconds(),"%Y"),time2dstr(seconds(),"%m"), ...
            time2dstr(seconds(),"%d"),...
            time2dstr(seconds(),"%Y"),time2dstr(seconds(),"%m"), ...
            time2dstr(seconds(),"%d"),...
            time2dstr(seconds(),"%H"),time2dstr(seconds(),"%M"), ...
            time2dstr(seconds(),"%S"));
    fprintf(stylefile,  "</text:p>");
    fprintf(stylefile,  "</style:region-right>");
    fprintf(stylefile,  "</style:header>");
    fprintf(stylefile,  "<style:header-left style:display=\"false\"/>");
    fprintf(stylefile,  "<style:footer>");
    fprintf(stylefile,  "<text:p>Page <text:page-number>1</text:" ...
        + "page-number> / <text:" ...
        + "page-count>99</text:page-count>");
    fprintf(stylefile,  "</text:p>");
    fprintf(stylefile,  "</style:footer>");
    fprintf(stylefile,  "<style:footer-left style:display=\"false\"/>");
    fprintf(stylefile,  "</style:master-page>");
    fprintf(stylefile,  "</office:master-styles>");
    fprintf(stylefile,  "</office:document-styles>");
  }}

  close(stylefile);
};


//
// write the settings file
//
_oooWriteSettings = function(mykind)
{
  if(!exist(mykind))
  { error("_oooWriteSettings: missing file type\n"); }

  if (!isdir(_builddir))
  { mkdir(_builddir); }

  setfile = _builddir + "/settings.xml";
  if (_ooodebug)
  { printf ("_oooWriteSettings: using settings file %s\n", setfile); }

  open(setfile,"w");

  if (mykind == "sxw")
  {
    fprintf(setfile, "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
    fprintf(setfile, "<!DOCTYPE office:document-settings PUBLIC \"-//OpenOffice." ...
        + "org//DTD OfficeDocument 1.0//EN\" \"office.dtd\">");
    fprintf(setfile, "<office:document-settings xmlns:office=\"http:" ...
        + "//openoffice.org/2000/office\" xmlns:xlink=\"http:" ...
            + "//www.w3.org/1999/xlink\" xmlns:config=\"http:" ...
            + "//openoffice.org/2001/config\" office:version=\"1.0\">");
    fprintf(setfile, "<office:settings>");
    fprintf(setfile, "<config:config-item-set config:name=\"view-settings\">");
    fprintf(setfile, "<config:config-item config:" ...
        + "name=\"ViewAreaTop\" config:type=\"int\">0</config:config-item>");
    fprintf(setfile, "<config:config-item config:name=\"ViewAreaLeft\" config:" ...
        + "type=\"int\">0</config:config-item>");
    fprintf(setfile, "<config:config-item config:name=\"ViewAreaWidth\" config:" ...
        + "type=\"int\">19846</config:config-item>");
    fprintf(setfile, "<config:config-item config:name=\"ViewAreaHeight\" config:" ...
        + "type=\"int\">10691</config:config-item>");
    fprintf(setfile, "<config:config-item config:name=\"ShowRedlineChanges\" config:" ...
        + "type=\"boolean\">true</config:config-item>");
    fprintf(setfile, "<config:config-item config:name=\"ShowHeaderWhileBrowsing\" config:" ...
        + "type=\"boolean\">false</config:config-item>");
    fprintf(setfile, "<config:config-item config:name=\"ShowFooterWhileBrowsing\" config:" ...
        + "type=\"boolean\">false</config:config-item>");
    fprintf(setfile, "<config:config-item config:name=\"InBrowseMode\" config:" ...
        + "type=\"boolean\">false</config:config-item>");
    fprintf(setfile, "<config:config-item-map-indexed config:name=\"Views\">");
    fprintf(setfile, "<config:config-item-map-entry>");
    fprintf(setfile, "<config:config-item config:name=\"ViewId\" config:" ...
        + "type=\"string\">view2</config:config-item>");
    fprintf(setfile, "<config:config-item config:name=\"ViewLeft\" config:" ...
        + "type=\"int\">3002</config:config-item>");
    fprintf(setfile, "<config:config-item config:name=\"ViewTop\" config:" ...
        + "type=\"int\">3002</config:config-item>");
    fprintf(setfile, "<config:config-item config:name=\"VisibleLeft\" config:" ...
        + "type=\"int\">0</config:config-item>");
    fprintf(setfile, "<config:config-item config:name=\"VisibleTop\" config:" ...
        + "type=\"int\">0</config:config-item>");
    fprintf(setfile, "<config:config-item config:name=\"VisibleRight\" config:" ...
        + "type=\"int\">19844</config:config-item>");
    fprintf(setfile, "<config:config-item config:name=\"VisibleBottom\" config:" ...
        + "type=\"int\">10689</config:config-item>");
    fprintf(setfile, "<config:config-item config:name=\"ZoomType\" config:" ...
        + "type=\"short\">0</config:config-item>");
    fprintf(setfile, "<config:config-item config:name=\"ZoomFactor\" config:" ...
        + "type=\"short\">100</config:config-item>");
    fprintf(setfile, "<config:config-item config:name=\"IsSelectedFrame\" config:" ...
        + "type=\"boolean\">false</config:config-item>");
    fprintf(setfile, "</config:config-item-map-entry>");
    fprintf(setfile, "</config:config-item-map-indexed>");
    fprintf(setfile, "</config:config-item-set>");
    fprintf(setfile, "<config:config-item-set config:name=\"configuration-settings\">");
    fprintf(setfile, "<config:config-item config:name=\"AddParaTableSpacing\" config:" ...
        + "type=\"boolean\">false</config:config-item>");
    fprintf(setfile, "<config:config-item config:name=\"PrintReversed\" config:" ...
        + "type=\"boolean\">false</config:config-item>");
    fprintf(setfile, "<config:config-item config:name=\"LinkUpdateMode\" config:" ...
        + "type=\"short\">1</config:config-item>");
    fprintf(setfile, "<config:config-item config:name=\"CharacterCompressionType\" config:" ...
        + "type=\"short\">0</config:config-item>");
    fprintf(setfile, "<config:config-item config:name=\"PrintSingleJobs\" config:" ...
        + "type=\"boolean\">false</config:config-item>");
    fprintf(setfile, "<config:config-item config:name=\"UpdateFromTemplate\" config:" ...
        + "type=\"boolean\">false</config:config-item>");
    fprintf(setfile, "<config:config-item config:name=\"PrintPaperFromSetup\" config:" ...
        + "type=\"boolean\">false</config:config-item>");
    fprintf(setfile, "<config:config-item config:name=\"PrintLeftPages\" config:" ...
        + "type=\"boolean\">true</config:config-item>");
    fprintf(setfile, "<config:config-item config:name=\"PrintTables\" config:" ...
        + "type=\"boolean\">true</config:config-item>");
    fprintf(setfile, "<config:config-item config:name=\"ChartAutoUpdate\" config:" ...
        + "type=\"boolean\">true</config:config-item>");
    fprintf(setfile, "<config:config-item config:name=\"PrintControls\" config:" ...
        + "type=\"boolean\">true</config:config-item>");
    fprintf(setfile, "<config:config-item config:name=\"PrinterSetup\" config:type=\"base64Binary\">");
    fprintf(setfile, "fAL+/0VQU09OIFBNLTc3MEMAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAARVBTT04gUE0tNzcwQwAAAAAAAAAAAAAAAAAAAAAAAAAWAAEAwgEAAAAAAAABAAhSAAAEdAAAM1ROVwEACABFUFNPTiBQTS03NzBDAAAAAAAAAAAAAAAAAAAAAAAAAAAEJAKUACYBD4uABwEACQCaCzQIAAABAAcAaAECAAAAaAEAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAQAAAAIAAAABAAAAIQEAAAAAAAAAAAAARExMTmFtZTE2PUVQSUdVSjNRLkRMTAAAAAAAAAAAAABETExOYW1lMzI9RVBJREEyMzAuRExMAAAAAAAAAAAAAEVQU09OIFBNLTc3MEMAAAAAAAAAAAAAAAAAAAAAAAAAAAMCAGgBaAEBAAAAAAAAAAABAAAJAEwLgQ9MC4EPZABoAWgBoAtxECoAKgAqAMYAoAtxECoAKgAqAMYAAAAAABQAAAAAAAAAAAAyAAAA/wAAAAAAAAAAAAAAAAABAAAAAAACAAAAAgAAAAEAAQAGAAYAAgAAAAAAAAACAAAAAAAAAAAAAAAAAAUAAQAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA0CJoLAQAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA</config:config-item>");
    fprintf(setfile, "<config:config-item config:name=\"PrintAnnotationMode\" config:" ...
        + "type=\"short\">0</config:config-item>");
    fprintf(setfile, "<config:config-item config:name=\"ApplyUserData\" config:" ...
        + "type=\"boolean\">true</config:config-item>");
    fprintf(setfile, "<config:config-item config:name=\"FieldAutoUpdate\" config:" ...
        + "type=\"boolean\">true</config:config-item>");
    fprintf(setfile, "<config:config-item config:name=\"SaveVersionOnClose\" config:" ...
        + "type=\"boolean\">false</config:config-item>");
    fprintf(setfile, "<config:config-item config:name=\"SaveGlobalDocumentLinks\" config:" ...
        + "type=\"boolean\">false</config:config-item>");
    fprintf(setfile, "<config:config-item config:name=\"IsKernAsianPunctuation\" config:" ...
        + "type=\"boolean\">false</config:config-item>");
    fprintf(setfile, "<config:config-item config:name=\"AlignTabStopPosition\" config:" ...
        + "type=\"boolean\">false</config:config-item>");
    fprintf(setfile, "<config:config-item config:name=\"CurrentDatabaseDataSource\" config:" ...
        + "type=\"string\"/>");
    fprintf(setfile, "<config:config-item config:name=\"PrinterName\" config:" ...
        + "type=\"string\">EPSON PM-770C</config:config-item>");
    fprintf(setfile, "<config:config-item config:name=\"PrintFaxName\" config:" ...
        + "type=\"string\"/>");
    fprintf(setfile, "<config:config-item config:name=\"PrintRightPages\" config:" ...
        + "type=\"boolean\">true</config:config-item>");
    fprintf(setfile, "<config:config-item config:name=\"AddParaTableSpacingAtStart\" config:" ...
        + "type=\"boolean\">false</config:config-item>");
    fprintf(setfile, "<config:config-item config:name=\"PrintProspect\" config:" ...
        + "type=\"boolean\">false</config:config-item>");
    fprintf(setfile, "<config:config-item config:name=\"PrintGraphics\" config:" ...
        + "type=\"boolean\">true</config:config-item>");
    fprintf(setfile, "<config:config-item config:name=\"CurrentDatabaseCommandType\" config:" ...
        + "type=\"int\">0</config:config-item>");
    fprintf(setfile, "<config:config-item config:name=\"PrintPageBackground\" config:" ...
        + "type=\"boolean\">true</config:config-item>");
    fprintf(setfile, "<config:config-item config:name=\"CurrentDatabaseCommand\" config:" ...
        + "type=\"string\"/>");
    fprintf(setfile, "<config:config-item config:name=\"PrintDrawings\" config:" ...
        + "type=\"boolean\">true</config:config-item>");
    fprintf(setfile, "<config:config-item config:name=\"PrintBlackFonts\" config:" ...
        + "type=\"boolean\">false</config:config-item>");
    fprintf(setfile, "</config:config-item-set>");
    fprintf(setfile, "</office:settings>");
    fprintf(setfile, "</office:document-settings>");

  } else { if (mykind == "sxc")
  {
    fprintf(setfile, "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");

    fprintf(setfile, "<!DOCTYPE office:document-settings PUBLIC \"-//OpenOffice." ...
        + "org//DTD OfficeDocument 1.0//EN\" \"office.dtd\">");

    fprintf(setfile, "<office:document-settings xmlns:office=\"http:" ...
        + "//openoffice.org/2000/office\" xmlns:xlink=\"http:" ...
            + "//www.w3.org/1999/xlink\" xmlns:config=\"http:" ...
            + "//openoffice.org/2001/config\" office:version=\"1.0\">");
    fprintf(setfile, "<office:settings>");
    fprintf(setfile, "<config:config-item-set config:name=\"view-settings\">");
    fprintf(setfile, "<config:config-item config:name=\"VisibleAreaTop\" config:" ...
        + "type=\"int\">0</config:config-item>");
    fprintf(setfile, "<config:config-item config:name=\"VisibleAreaLeft\" config:" ...
        + "type=\"int\">0</config:config-item>");
    fprintf(setfile, "<config:config-item config:name=\"VisibleAreaWidth\" config:" ...
        + "type=\"int\">2258</config:config-item>");
    fprintf(setfile, "<config:config-item config:name=\"VisibleAreaHeight\" config:" ...
        + "type=\"int\">451</config:config-item>");
    fprintf(setfile, "<config:config-item-map-indexed config:name=\"Views\">");
    fprintf(setfile, "<config:config-item-map-entry>");
    fprintf(setfile, "<config:config-item config:name=\"ViewId\" config:" ...
        + "type=\"string\">View1</config:config-item>");
    fprintf(setfile, "<config:config-item-map-named config:name=\"Tables\">");
    fprintf(setfile, "<config:config-item-map-entry config:name=\"Sheet1\">");
    fprintf(setfile, "<config:config-item config:name=\"CursorPositionX\" config:" ...
        + "type=\"int\">0</config:config-item>");
    fprintf(setfile, "<config:config-item config:name=\"CursorPositionY\" config:" ...
        + "type=\"int\">0</config:config-item>");
    fprintf(setfile, "<config:config-item config:name=\"HorizontalSplitMode\" config:" ...
        + "type=\"short\">0</config:config-item>");
    fprintf(setfile, "<config:config-item config:name=\"VerticalSplitMode\" config:" ...
        + "type=\"short\">0</config:config-item>");
    fprintf(setfile, "<config:config-item config:name=\"HorizontalSplitPosition\" config:" ...
        + "type=\"int\">0</config:config-item>");
    fprintf(setfile, "<config:config-item config:name=\"VerticalSplitPosition\" config:" ...
        + "type=\"int\">0</config:config-item>");
    fprintf(setfile, "<config:config-item config:name=\"ActiveSplitRange\" config:" ...
        + "type=\"short\">2</config:config-item>");
    fprintf(setfile, "<config:config-item config:name=\"PositionLeft\" config:" ...
        + "type=\"int\">0</config:config-item>");
    fprintf(setfile, "<config:config-item config:name=\"PositionRight\" config:" ...
        + "type=\"int\">0</config:config-item>");
    fprintf(setfile, "<config:config-item config:name=\"PositionTop\" config:" ...
        + "type=\"int\">0</config:config-item>");
    fprintf(setfile, "<config:config-item config:name=\"PositionBottom\" config:" ...
        + "type=\"int\">0</config:config-item>");
    fprintf(setfile, "</config:config-item-map-entry>");
    fprintf(setfile, "</config:config-item-map-named>");
    fprintf(setfile, "<config:config-item config:name=\"ActiveTable\" config:type=\"string\">");
    fprintf(setfile, "Sheet1</config:config-item>");
    fprintf(setfile, "<config:config-item config:name=\"HorizontalScrollbarWidth\" config:" ...
        + "type=\"int\">270</config:config-item>");
    fprintf(setfile, "<config:config-item config:name=\"ZoomType\" config:" ...
        + "type=\"short\">0</config:config-item>");
    fprintf(setfile, "<config:config-item config:name=\"ZoomValue\" config:" ...
        + "type=\"int\">100</config:config-item>");
    fprintf(setfile, "<config:config-item config:name=\"PageViewZoomValue\" config:" ...
        + "type=\"int\">60</config:config-item>");
    fprintf(setfile, "<config:config-item config:name=\"ShowPageBreakPreview\" config:" ...
        + "type=\"boolean\">false</config:config-item>");
    fprintf(setfile, "<config:config-item config:name=\"ShowZeroValues\" config:" ...
        + "type=\"boolean\">true</config:config-item>");
    fprintf(setfile, "<config:config-item config:name=\"ShowNotes\" config:" ...
        + "type=\"boolean\">true</config:config-item>");
    fprintf(setfile, "<config:config-item config:name=\"ShowGrid\" config:" ...
        + "type=\"boolean\">true</config:config-item>");
    fprintf(setfile, "<config:config-item config:name=\"GridColor\" config:" ...
        + "type=\"long\">12632256</config:config-item>");
    fprintf(setfile, "<config:config-item config:name=\"ShowPageBreaks\" config:" ...
        + "type=\"boolean\">true</config:config-item>");
    fprintf(setfile, "<config:config-item config:name=\"HasColumnRowHeaders\" config:" ...
        + "type=\"boolean\">true</config:config-item>");
    fprintf(setfile, "<config:config-item config:name=\"HasSheetTabs\" config:" ...
        + "type=\"boolean\">true</config:config-item>");
    fprintf(setfile, "<config:config-item config:name=\"IsOutlineSymbolsSet\" config:" ...
        + "type=\"boolean\">true</config:config-item>");
    fprintf(setfile, "<config:config-item config:name=\"IsSnapToRaster\" config:" ...
        + "type=\"boolean\">false</config:config-item>");
    fprintf(setfile, "<config:config-item config:name=\"RasterIsVisible\" config:" ...
        + "type=\"boolean\">false</config:config-item>");
    fprintf(setfile, "<config:config-item config:name=\"RasterResolutionX\" config:" ...
        + "type=\"int\">1000</config:config-item>");
    fprintf(setfile, "<config:config-item config:name=\"RasterResolutionY\" config:" ...
        + "type=\"int\">1000</config:config-item>");
    fprintf(setfile, "<config:config-item config:name=\"RasterSubdivisionX\" config:" ...
        + "type=\"int\">1</config:config-item>");
    fprintf(setfile, "<config:config-item config:name=\"RasterSubdivisionY\" config:" ...
        + "type=\"int\">1</config:config-item>");
    fprintf(setfile, "<config:config-item config:name=\"IsRasterAxisSynchronized\" config:" ...
        + "type=\"boolean\">true</config:config-item>");
    fprintf(setfile, "</config:config-item-map-entry>");
    fprintf(setfile, "</config:config-item-map-indexed>");
    fprintf(setfile, "</config:config-item-set>");
    fprintf(setfile, "<config:config-item-set config:name=\"configuration-settings\">");
    fprintf(setfile, "<config:config-item config:name=\"ShowZeroValues\" config:" ...
        + "type=\"boolean\">true</config:config-item>");
    fprintf(setfile, "<config:config-item config:name=\"ShowNotes\" config:" ...
        + "type=\"boolean\">true</config:config-item>");
    fprintf(setfile, "<config:config-item config:name=\"ShowGrid\" config:" ...
        + "type=\"boolean\">true</config:config-item>");
    fprintf(setfile, "<config:config-item config:name=\"GridColor\" config:" ...
        + "type=\"long\">12632256</config:config-item>");
    fprintf(setfile, "<config:config-item config:name=\"ShowPageBreaks\" config:" ...
        + "type=\"boolean\">true</config:config-item>");
    fprintf(setfile, "<config:config-item config:name=\"LinkUpdateMode\" config:" ...
        + "type=\"short\">3</config:config-item>");
    fprintf(setfile, "<config:config-item config:name=\"HasColumnRowHeaders\" config:" ...
        + "type=\"boolean\">true</config:config-item>");
    fprintf(setfile, "<config:config-item config:name=\"HasSheetTabs\" config:" ...
        + "type=\"boolean\">true</config:config-item>");
    fprintf(setfile, "<config:config-item config:name=\"IsOutlineSymbolsSet\" config:" ...
        + "type=\"boolean\">true</config:config-item>");
    fprintf(setfile, "<config:config-item config:name=\"IsSnapToRaster\" config:" ...
        + "type=\"boolean\">false</config:config-item>");
    fprintf(setfile, "<config:config-item config:name=\"RasterIsVisible\" config:" ...
        + "type=\"boolean\">false</config:config-item>");
    fprintf(setfile, "<config:config-item config:name=\"RasterResolutionX\" config:" ...
        + "type=\"int\">1000</config:config-item>");
    fprintf(setfile, "<config:config-item config:name=\"RasterResolutionY\" config:" ...
        + "type=\"int\">1000</config:config-item>");
    fprintf(setfile, "<config:config-item config:name=\"RasterSubdivisionX\" config:" ...
        + "type=\"int\">1</config:config-item>");
    fprintf(setfile, "<config:config-item config:name=\"RasterSubdivisionY\" config:" ...
        + "type=\"int\">1</config:config-item>");
    fprintf(setfile, "<config:config-item config:name=\"IsRasterAxisSynchronized\" config:" ...
        + "type=\"boolean\">true</config:config-item>");
    fprintf(setfile, "<config:config-item config:name=\"AutoCalculate\" config:" ...
        + "type=\"boolean\">true</config:config-item>");
    fprintf(setfile, "<config:config-item config:name=\"PrinterName\" config:" ...
        + "type=\"string\">EPSON PM-770C</config:config-item>");
    fprintf(setfile, "<config:config-item config:name=\"PrinterSetup\" config:type=\"base64Binary\">");
    fprintf(setfile, "fAL+/0VQU09OIFBNLTc3MEMAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAARVBTT04gUE0tNzcwQwAAAAAAAAAAAAAAAAAAAAAAAAAWAAEAwgEAAAAAAAABAAhSAAAEdAAAM1ROVwEACABFUFNPTiBQTS03NzBDAAAAAAAAAAAAAAAAAAAAAAAAAAAEJAKUACYBD4uABwEACQCaCzQIAAABAAcAaAECAAAAaAEAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAQAAAAIAAAABAAAAIQEAAAAAAAAAAAAARExMTmFtZTE2PUVQSUdVSjNRLkRMTAAAAAAAAAAAAABETExOYW1lMzI9RVBJREEyMzAuRExMAAAAAAAAAAAAAEVQU09OIFBNLTc3MEMAAAAAAAAAAAAAAAAAAAAAAAAAAAMCAGgBaAEBAAAAAAAAAAABAAAJAEwLgQ9MC4EPZABoAWgBoAtxECoAKgAqAMYAoAtxECoAKgAqAMYAAAAAABQAAAAAAAAAAAAyAAAA/wAAAAAAAAAAAAAAAAABAAAAAAACAAAAAgAAAAEAAQAGAAYAAgAAAAAAAAACAAAAAAAAAAAAAAAAAAUAAQAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA0CJoLAQAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA</config:config-item>");
    fprintf(setfile, "<config:config-item config:name=\"ApplyUserData\" config:" ...
        + "type=\"boolean\">true</config:config-item>");
    fprintf(setfile, "<config:config-item config:name=\"CharacterCompressionType\" config:" ...
        + "type=\"short\">0</config:config-item>");
    fprintf(setfile, "<config:config-item config:name=\"IsKernAsianPunctuation\" config:" ...
        + "type=\"boolean\">false</config:config-item>");
    fprintf(setfile, "<config:config-item config:name=\"SaveVersionOnClose\" config:" ...
        + "type=\"boolean\">false</config:config-item>");
    fprintf(setfile, "<config:config-item config:name=\"UpdateFromTemplate\" config:" ...
        + "type=\"boolean\">false</config:config-item>");
    fprintf(setfile, "</config:config-item-set>");
    fprintf(setfile, "</office:settings>");
    fprintf(setfile, "</office:document-settings>");
  }}

  close(setfile);
};

ooo.special = function(style, listnum, mytext)
{
  // Adds special data to the document: brakes, tables?
  if(!exist(style))
  { error("ooo.special: requires at least one argument"); }

  // Styles
  if (style == "pagebreak")
  {
    if (exist(_options.pagebreak))
    {
      stylenum = _options.pagebreak;
      _documenttext = [_documenttext, "<text:p text:style-name=\"" + stylenum + "\"/>" ];
    } else {
    {
      stylenum = "P" + text(_options.nextstyle);
      _options.nextstyle = _options.nextstyle + 1;
      _options.pagebreak = stylenum;
      _autostyles = [_autostyles, ...
        "<style:style style:name=\""+ stylenum + ...
        "\" style:family=\"paragraph\" style:parent-style-name=\"Standard\">"];
      _autostyles = [_autostyles, "<style:properties fo:break-before=\"page\"/>"];
      _autostyles = [_autostyles, "</style:style>"];
      _documenttext = [_documenttext, "<text:p text:style-name=\""+stylenum+"\"/>"];
    }}
    return "ok";
  }

  if (style == "list-ordered")
  {
    mytext = _oooCleanText(mytext);

    // listnum="L[0-9]"
    if (exist(listnum))
    {
      // Good listnum
      stylenum = _options.[listnum];

      if(!exist(stylenum))
      {
        _errormessage = "Invalid list number \""+listnum+"\".";
        return "error";
      }
      _documenttext = [_documenttext, ...
          "<text:ordered-list text:style-name=\""+listnum+"\" text:continue-numbering=\"true\">"];
      _documenttext=[_documenttext, "<text:list-item>"];
      _documenttext=[_documenttext, ...
          "<text:p text:style-name=\""+stylenum+"\">$text</text:p>"];
      _documenttext=[_documenttext, "</text:list-item>"];
      _documenttext=[_documenttext, "</text:ordered-list>"];
      return listnum;
    } else {
    {
      // Needs a listnum
      stylenum = "P" + text(_options.nextstyle);
      _options.nextstyle = _options.nextstyle + 1;
      listnum = "L" + text(_options.nextlist);
      _options.nextlist=_options.nextlist+1;
      _options.[listnum] = stylenum;
      // The style
      _autostyles=[_autostyles, ...
          "<style:style style:name=\""+stylenum+"\" style:family=\"paragraph\" style:" ...
          + "parent-style-name=\"Standard\" style:list-style-name=\""+listnum+"\"/>"];
      // The list information
      _autolists=[_autolists, "<text:list-style style:name=\""+listnum+"\">"];
      _autolists=[_autolists, "<text:list-level-style-number text:" ...
          +"level=\"1\" text:style-name=\"Numbering Symbols\" style:" ...
          +"num-suffix=\".\" style:num-format=\"1\">"];
      _autolists=[_autolists, "<style:properties text:min-label-width=\"0.499cm\"/>"];
      _autolists=[_autolists, "</text:list-level-style-number>"];
      _autolists=[_autolists, "<text:list-level-style-number text:level=\"2\" text:" ...
          +"style-name=\"Numbering Symbols\" style:num-suffix=\".\" style:num-format=\"1\">"];
      _autolists=[_autolists, "<style:properties text:space-before=\"0.501cm\" text:" ...
          +"min-label-width=\"0.499cm\"/>"];
      _autolists=[_autolists, "</text:list-level-style-number>"];
      _autolists=[_autolists, "<text:list-level-style-number text:level=\"3\" text:" ...
          +"style-name=\"Numbering Symbols\" style:num-suffix=\".\" style:num-format=\"1\">"];
          _autolists=[_autolists, "<style:properties text:space-before=\"1cm\" text:" ...
              +"min-label-width=\"0.499cm\"/>"];
      _autolists=[_autolists, "</text:list-level-style-number>"];
      _autolists=[_autolists, "<text:list-level-style-number text:level=\"4\" text:" ...
          +"style-name=\"Numbering Symbols\" style:num-suffix=\".\" style:num-format=\"1\">"];
      _autolists=[_autolists, "<style:properties text:space-before=\"1.501cm\" text:" ...
          +"min-label-width=\"0.499cm\"/>"];
      _autolists=[_autolists, "</text:list-level-style-number>"];
      _autolists=[_autolists, "<text:list-level-style-number text:level=\"5\" text:" ...
          +"style-name=\"Numbering Symbols\" style:num-suffix=\".\" style:num-format=\"1\">"];
      _autolists=[_autolists, "<style:properties text:space-before=\"2cm\" text:" ...
          +"min-label-width=\"0.499cm\"/>"];
      _autolists=[_autolists, "</text:list-level-style-number>"];
      _autolists=[_autolists, "<text:list-level-style-number text:level=\"6\" text:" ...
          +"style-name=\"Numbering Symbols\" style:num-suffix=\".\" style:num-format=\"1\">"];
      _autolists=[_autolists, "<style:properties text:space-before=\"2.501cm\" text:" ...
          +"min-label-width=\"0.499cm\"/>"];
      _autolists=[_autolists, "</text:list-level-style-number>"];
      _autolists=[_autolists, "<text:list-level-style-number text:level=\"7\" text:" ...
          +"style-name=\"Numbering Symbols\" style:num-suffix=\".\" style:num-format=\"1\">"];
      _autolists=[_autolists, "<style:properties text:space-before=\"3.001cm\" text:" ...
          +"min-label-width=\"0.499cm\"/>"];
      _autolists=[_autolists, "</text:list-level-style-number>"];
      _autolists=[_autolists, "<text:list-level-style-number text:level=\"8\" text:" ...
          +"style-name=\"Numbering Symbols\" style:num-suffix=\".\" style:num-format=\"1\">"];
      _autolists=[_autolists, "<style:properties text:space-before=\"3.502cm\" text:" ...
          +"min-label-width=\"0.499cm\"/>"];
      _autolists=[_autolists, "</text:list-level-style-number>"];
      _autolists=[_autolists, "<text:list-level-style-number text:level=\"9\" text:" ...
          +"style-name=\"Numbering Symbols\" style:num-suffix=\".\" style:num-format=\"1\">"];
      _autolists=[_autolists, "<style:properties text:space-before=\"4.001cm\" text:" ...
          +"min-label-width=\"0.499cm\"/>"];
      _autolists=[_autolists, "</text:list-level-style-number>"];
      _autolists=[_autolists, "<text:list-level-style-number text:level=\"10\" text:" ...
          +"style-name=\"Numbering Symbols\" style:num-suffix=\".\" style:num-format=\"1\">"];
      _autolists=[_autolists, "<style:properties text:space-before=\"4.502cm\" text:" ...
          +"min-label-width=\"0.499cm\"/>"];
      _autolists=[_autolists, "</text:list-level-style-number>"];
      _autolists=[_autolists, "</text:list-style>"];

      // Text for the document
      _documenttext=[_documenttext, "<text:ordered-list text:style-name=\""+listnum+"\">"];
      _documenttext=[_documenttext, "<text:list-item>"];
      _documenttext=[_documenttext, "<text:p text:style-name=\""+stylenum+"\">"];
      _documenttext=[_documenttext, mytext];
      _documenttext=[_documenttext, "</text:p>"];
      _documenttext=[_documenttext, "</text:list-item>"];
      _documenttext=[_documenttext, "</text:ordered-list>"];
      return listnum;
    }}}

  _errormessage = "I did not understand the oooSpecial style";
  return "error";
};


_oooWriteContent = function( mykind )
{
  // writes the content.xml file
  if(!exist(mykind))
  { error("_oooWriteContent: missing type of the file\n"); }

  if (!isdir(_builddir))
  { mkdir(_builddir); }

  contentfile = _builddir + "/content.xml";
  if (_ooodebug)
  { printf ("_oooWriteContent: using settings file %s\n", contentfile); }

  open(contentfile,"w");

  if (mykind == "sxw")
  {
    // writer file
    // Encoding information
    fprintf(contentfile, "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
    // DOCTYPE information
    fprintf(contentfile, "<!DOCTYPE office:document-content PUBLIC \"-//OpenOffice." ...
        + "org//DTD OfficeDocument 1.0//EN\" \"office.dtd\">");
    // Information about reading this file
    fprintf(contentfile, "<office:document-content xmlns:office=\"http:" ...
        + "//openoffice.org/2000/office\" xmlns:style=\"http://openoffice.org/2000/style\" xmlns:" ...
            + "text=\"http://openoffice.org/2000/text\" xmlns:table=\"http:" ...
            + "//openoffice.org/2000/table\" xmlns:draw=\"http://openoffice." ...
            + "org/2000/drawing\" xmlns:fo=\"http://www.w3.org/1999/XSL/Format\" xmlns:" ...
            + "xlink=\"http://www.w3.org/1999/xlink\" xmlns:" ...
            + "number=\"http://openoffice.org/2000/datastyle\" xmlns:svg=\"http://www.w3." ...
            + "org/2000/svg\" xmlns:chart=\"http://openoffice.org/2000/chart\" xmlns:" ...
            + "dr3d=\"http://openoffice.org/2000/dr3d\" xmlns:" ...
            + "math=\"http://www.w3.org/1998/Math/MathML\" xmlns:" ...
            + "form=\"http://openoffice.org/2000/form\" xmlns:" ...
            + "script=\"http://openoffice.org/2000/script\" office:class=\"text\" office:" ...
            + "version=\"1.0\">");
    // No scripts
    fprintf(contentfile, "<office:script/>");
    // Fonts
    if (length(_fontdecls) > 0)
    {
      fprintf(contentfile, "<office:font-decls>");
      for (line in _fontdecls)
      {
        if (_ooodebug)
        { printf("_oooWriteContent: %s\n", line); }
        fprintf(contentfile, "%s", line);

      }
      fprintf(contentfile, "</office:font-decls>");
    }

    // Styles line bold, center, etc.
    if (length(_autostyles) > 0)
    {
      fprintf(contentfile, "<office:automatic-styles>");
      for (_autostyle in _autostyles)
      {
        if (_ooodebug)
        { printf("_oooWriteContent: %s\n", _autostyle); }
        fprintf(contentfile, "%s", _autostyle);
      }
      fprintf(contentfile, "</office:automatic-styles>");
    } else {
      // No automatic styles
      fprintf(contentfile, "<office:automatic-styles/>");
    }

    // Body of context
    fprintf(contentfile, "<office:body>");
    fprintf(contentfile, "<text:sequence-decls>");
    fprintf(contentfile, "<text:sequence-decl text:" ...
        + "display-outline-level=\"0\" text:name=\"Illustration\"/>");
    fprintf(contentfile, "<text:sequence-decl text:display-outline-level" ...
        + "=\"0\" text:name=\"Table\"/>");
    fprintf(contentfile, "<text:sequence-decl text:display-outline-level" ...
        + "=\"0\" text:name=\"Text\"/>");
    fprintf(contentfile, "<text:sequence-decl text:display-outline-level" ...
        + "=\"0\" text:name=\"Drawing\"/>");
    fprintf(contentfile, "</text:sequence-decls>");

    // Text
    if (length(_documenttext) > 0)
    {
      for (i in _documenttext)
      {
          if (_ooodebug)
          { printf("_oooWriteContent: %s\n", i); }
          fprintf(contentfile, "%s", gsub("&#10;", "\n", i).string);
      }
    } else {
      fprintf(contentfile, "<text:p text:style-name=\"Standard\"/>");
    }
    fprintf(contentfile, "</office:body>");
    fprintf(contentfile, "</office:document-content>");

  } else { if (mykind == "sxc")
  {
    // Encoding information
    fprintf(contentfile, "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
    // DOCTYPE information
    fprintf(contentfile, "<!DOCTYPE office:document-content PUBLIC \"-//OpenOffice." ...
        + "org//DTD OfficeDocument 1.0//EN\" \"office.dtd\">");
    // Information about decoding XML
    fprintf(contentfile, "<office:document-content xmlns:office=\"http:" ...
        + "//openoffice.org/2000/office\" xmlns:style=\"http:" ...
            + "//openoffice.org/2000/style\" xmlns:text=\"http://openoffice." ...
            + "org/2000/text\" xmlns:table=\"http://openoffice.org/2000/table\" xmlns:" ...
            + "draw=\"http://openoffice.org/2000/drawing\" xmlns:" ...
            + "fo=\"http://www.w3.org/1999/XSL/Format\" xmlns:xlink=\"http://www.w3." ...
            + "org/1999/xlink\" xmlns:number=\"http://openoffice.org/2000/datastyle\" xmlns:" ...
            + "svg=\"http://www.w3.org/2000/svg\" xmlns:chart=\"http://openoffice." ...
            + "org/2000/chart\" xmlns:dr3d=\"http://openoffice.org/2000/dr3d\" xmlns:" ...
            + "math=\"http://www.w3.org/1998/Math/MathML\" xmlns:" ...
            + "form=\"http://openoffice.org/2000/form\" xmlns:" ...
            + "script=\"http://openoffice.org/2000/script\" office:class=\"spreadsheet\" office:" ...
            + "version=\"1.0\">");
    // Script information tag?
    fprintf(contentfile, "<office:script/>");

    // Fonts
    if (length(_fontdecls)>0)
    {
      fprintf(contentfile, "<office:font-decls>");
      for (line in _fontdecls)
      {
        if (_ooodebug)
        { printf("_oooWriteContent: %s\n", line); }
        fprintf(contentfile, "%s", line);
      }
      fprintf(contentfile, "</office:font-decls>");
    }

    // Styles
    fprintf(contentfile, "<office:automatic-styles>");

    // table-column styles
    fprintf(contentfile, "<style:style style:name=\"co1\" style:family=\"table-column\">");
    fprintf(contentfile, "<style:properties fo:break-before=\"auto\" style:column-width=\"2.267cm\"/>");
    fprintf(contentfile, "</style:style>");

    // Default values
    columnwidth = "2.267";
    for (x in 1:_cellhash.xmax)
    {
      width = _cellhash.["column " + text(x) + " width"];
      if(!exist(width))
      { width = columnwidth; }
      if (width != columnwidth)
      {
        // Create a new style
        mynum = _options.["nextstyle table-column"];
        _options.["nextstyle table-column"] = _options.["nextstyle table-column"] + 1;
        fprintf(contentfile, "<style:style style:name=\"co"+text(mynum) ...
            + "\" style:family=\"table-column\">");
        fprintf(contentfile, "<style:properties fo:break-before=\"auto\" style:" ...
            + "column-width=\""+text(width)+"cm\"/>");
        fprintf(contentfile, "</style:style>");
        _cellhash.["column " + text(x) + " style"] = "co"+text(mynum);
      } else {
      {
        _cellhash.["column "+text(x)+" style"] = "co1";
      }}
    }

    // table-row styles
    fprintf(contentfile, "<style:style style:name=\"ro1\" style:family=\"table-row\">");
    fprintf(contentfile, "<style:properties fo:break-before=\"auto\"/>");
    fprintf(contentfile, "</style:style>");

    // table styles
    fprintf(contentfile, "<style:style style:name=\"ta1\" style:family=\"table\" style:" ...
        + "master-page-name=\"Default\">");
    fprintf(contentfile, "<style:properties table:display=\"true\"/>");
    fprintf(contentfile, "</style:style>");

    // table-cell styles
    if (length(_ooolib_sxc_cell_styles) > 0)
    {
      for (line in _ooolib_sxc_cell_styles)
      { fprintf(contentfile, "%s", line); }
    }
    fprintf(contentfile, "</office:automatic-styles>");

    // Beginning of document content
    fprintf(contentfile, "<office:body>");

    fprintf(contentfile, "<table:table table:name=\"Sheet1\" table:style-name=\"ta1\">");

    // First tell the document what each of the column style codes are.
    for(x in 1:_cellhash.xmax)
    {
      // Each Column
      style = _cellhash.["column "+text(x)+" style"];
      fprintf(contentfile, "<table:table-column table:style-name=\"" + style ...
          + "\" table:number-columns-repeated=\"1\" table:default-cell-style-" ...
          + "name=\"Default\"/>");
    }

    // Row by row print the cells
    for (y in 1:_cellhash.ymax)
    {
      // One column cell at a time down the row
      fprintf(contentfile, "<table:table-row table:style-name=\"ro1\">");
      for (x in 1:_cellhash.xmax)
      {
        // Load information about the cell
        celltype = _cellhash.[text(x)+" "+text(y)+" type"];
        if (exist(_cellhash.[text(x)+" "+text(y)+" value"]))
        {
          val      = _cellhash.[text(x)+" "+text(y)+" value"];
        } else {
          val = "";
        }
        style    = _cellhash.[text(x)+" "+text(y)+" style"];

        if (!exist(style))
        { style = "ro1"; }

        if (!exist(celltype))
        { celltype = "";}

        if (celltype == "text")
        {
          fprintf(contentfile, "<table:table-cell table:style-name=\""+style+"\">");
          fprintf(contentfile, "<text:p>"+val+"</text:p>");
          fprintf(contentfile, "</table:table-cell>");
        } else { if (celltype == "float")
        {
          fprintf(contentfile, "<table:table-cell table:value-type=\"float\" table:" ...
              + "value=\"%s\" table:style-name=\"%s\">", val, style);
          fprintf(contentfile, "<text:p>%s</text:p>", val);
          fprintf(contentfile, "</table:table-cell>");

        } else { if (celltype == "formula")
        {
          formula = _cellhash.[text(x)+" "+text(y)+" formula"];

          // cell-formula code from Vladimir Vuksan
          fprintf(contentfile, "<table:table-cell table:value-type=\"float\" table:" ...
              + "formula=\"" + formula + "\" table:value=\""+val+"\" table:style-name=\"" ...
              + style + "\">");
          fprintf(contentfile, "<text:p>");
          if (exist(_cellhash.[text(x)+" "+text(y)+" value-print"]))
          { fprintf(contentfile, "%s", _cellhash.[text(x)+" "+text(y)+" value-print"]); }
          fprintf(contentfile, "</text:p>");
          fprintf(contentfile, "</table:table-cell>");

        } else {
          fprintf(contentfile, "<table:table-cell/>");
        }}}
      }
      fprintf(contentfile, "</table:table-row>");
    }

    // Closing document
    fprintf(contentfile, "</table:table>");
    fprintf(contentfile, "</office:body>");
    fprintf(contentfile, "</office:document-content>");
  }}

  close(contentfile);
};


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// utility functions
//
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++

// write a matrix
ooo.writem = function(data, loc_11, row_label, col_label)
{
  global (ooo);

  _x1 = 1;
  _y1 = 1;

  _this_func = "ooo.writem";
  if (!exist(data))
  {
    printf(_this_func + ": " + "missing data matrix\n");
    return (1);
  }

  if (!exist(loc_11))
  {
    printf(_this_func + ": " + "missing location of the top-left corner where matrix goes\n");
    return (2);
  }

  if (class(loc_11) != "num" || loc_11.n!=2)
  {
    printf(_this_func + ": " + "top-left corner location is not correct\n");
    return (3);
  }

  _x1 = loc_11[2];
  _y1 = loc_11[1];

  // write data matrix
  if (data.n >0)
  {
    _xn = _x1 + data.nc - 1;
    _yn = _y1 + data.nr - 1;
    for (_x in _x1:_xn)
    {
      for (_y in _y1:_yn)
      {
        _d = data[_y - _y1 + 1; _x - _x1 + 1];
        if (class(data) == "string")
        {
          ooo.set ("cell-loc",  _x, _y);
          ooo.data("cell-text", _d);
        } else { if (class(data) == "num")
        {
          if (!isnan(_d))
          {
            ooo.set ("cell-loc",  _x, _y);
            ooo.data("cell-float", _d);
          }
        }}
        ooo.set("justify", "center");
      }
    }
  }

  // write row_labels
  if (exist(row_label))
  {
    if (class(row_label)!="list")
    {
      // default formatting for text/number row labels
      ooo.set ("bold", "on");
      if(_x1 > 1)
      {
        for (_j in range(row_label))
        {
          ooo.set ("cell-loc",  _x1 - 1, _y1 + _j - 1);
          if (class(row_label) == "string")
          {
            ooo.data("cell-text", row_label[_j]);
          } else { if (class(row_label) == "num")
          {
            ooo.data("cell-float", row_label[_j]);
          }}
        }
      }
      ooo.set("bold","off");
    }
  }

  // write column labels
  if (exist(col_label))
  {
    ooo.set ("bold", "on");

    if(_y1 > 1)
    {
      for (_i in range(col_label))
      {
        ooo.set ("cell-loc",  _x1 + _i - 1, _y1 - 1);
        if (class(col_label) == "string")
        {
          ooo.data("cell-text", col_label[_i]);
        } else { if (class(col_label) == "num")
        {
          ooo.data("cell-float", col_label[_i]);
        }}
      }
    }

    ooo.set("bold","off");
  }

};














