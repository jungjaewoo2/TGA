// file: rlabplus/lib.r/libstdio.r
//
// This file is a part of RLaB + rlabplus
// Copyright (C) 2013  Marijan Kostrun
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
//
// See the file ../COPYING
//
// some useful functions that depend on system- and reads- call
// to linux system
static(INIT);
if (exist(INIT))
{ EOF }

static(_char_table);

_char_table = [32:122];

inputs = function(timeout, echo, lim)
{
  DT = 0.1; // minimum time interval for tty device

  if (!exist(lim))
  { lim = 0; }

  if (!exist(echo))
  { echo = 1; }

  if (!exist(timeout))
  { timeout = 5; }

  // we use this system timer
  ig = 32;

  x = 0;
  r = int([]);

  tic(ig);
  while(x != 10)
  {
    dt = max(timeout - toc(ig), DT);
    if (dt == DT)
    { break; }

    x = showkey(dt);  // zero implies nothing was pressed during time 'dt'

    if (x!=10)
    {
      if (length(lim) == 1)
      {
        // accept it only if greater than 'lim'
        if (x > lim)
        {
          r = [r, x];
        }
        else
        {
          continue;
        }
      }
      else
      {
        // accept only if input exists in the 'lim' array
        _i = find(x == lim);
        if(isempty(_i))
        { continue; }

        r = [r, x];
      }

      switch(echo)
      {
        case 1:
          if (x >= 32 && x<=122)
          {
            printf(char(x));
          }
          else
          {
            printf("\\" + text(x,"{%i}"));
          }
          break;

        case 2:
          printf("*");
          break;

        case 3:
          j = int(4 * uniform());
          for ( j1 in 1:j)
          { printf("*"); }
          break;

        case 4:
          printf( char(_char_table[length(_char_table) * uniform() + 1]) );
          break;

        case 5:
          j = int(4*uniform());
          for ( j1 in 1:j)
          { printf( char(_char_table[length(_char_table) * uniform() + 1]) ); }
          break;
      }
    }
  }

  if (isempty(r))
  {
    rval = "";
  }
  else
  {
    rval = char(r);
  }

  return rval;
};

//
// prompt user for input and wait for response
//  p -> text of the prompt
//  d -> default choice
//  c -> available choices
//  opt -> options
//
prompt = function ( p, d, c, opt)
{
  _this_function = "prompt: ";
  _pco = 5;   // print that many choice options in ()
  _qm = "?";  // question mark after prompt
  _cs = 0;    // don't clear screen before printing the prompt
  _nc = 0;    // number of screen clearings
  _to = 0;    // do not timeout the user:
  _ec = 1;    //    default echo mode for timed input
  _lm = 31;   //    accept only printable characters, ignore the characters below
  _pc = "";   // print choices (alternatively, user can disable them!)

  // figure the class of the input:
  //  if default value is provided then this overrides 'c'
  if (!exist(c)||!exist(d))
  { error(_this_function + ": 'default' and 'choice' must be provided even if empty!\n"); }

  if (!isempty(c) && isempty(d))
  {
    _dc = type(c);
    if (type(c)=="real")
    { d = []; }
    if (type(c)=="string")
    { d = blank(0,0); }
  }
  else if (isempty(c) && !isempty(d))
  {
    _dc = type(d);
    if (type(d)=="real")
    { c = []; }
    if (type(d)=="string")
    { c = blank(0,0); }
  }
  else
  {
    if (class(c)!=class(d))
    { error("prompt: 'default' and 'choice' cannot have different classes!\n"); }

    _dc = type(c);
  }

  // special case: d is the same as c
  // then just return
  if (all(d==c) && length(d)==1)
  { return c; }

  // default options
  if (class(opt)=="list")
  {
    if (exist(opt.timeout))
    {
      if (opt.timeout >= 0)
      { _to = opt.timeout; }
    }

    // this one supersedes all the other values
    // expected for the input
    if (exist(opt.input_class))
    {
      if (opt.input_class == "int")
      {
        _dc = "int";
      }
      else if (opt.input_class == "string")
      {
        _dc = "string";
      }
    }

    if (class(opt.question_mark)=="string")
    { _qm = opt.question_mark; }

    if (exist(opt.clear_screen))
    {
      if (opt.clear_screen==0 || opt.clear_screen==1)
      {
        _cs = opt.clear_screen;
        _nc = 1;
      }
    }

    if (exist(opt.echo))
    {
      if (opt.echo>=1 && opt.echo<=5)
      { _ec = floor(opt.echo); }
    }

    if (exist(opt.limit))
    {
      if (all(opt.limit>=1) && all(opt.echo<=255))
      { _lm = floor(opt.limit); }
    }

    if (class(opt.print_choices) == "string")
    { _pc = opt.print_choices; }

    if (class(opt.use_eval) == "num")
    { _ue = opt.use_eval; }
  }

  // print default value if provided, or (n/a) if none such
  _sd = "";
  if (exist(d))
  {
    if (isempty(d))
    {
      _sd = " (n/a)";
    }
    else if (class(d)=="string")
    {
      _sd = "[" + d[1] + "]";
    }
    else
    {
      _fmt = "%g";
      if (type(d) == "int")
      { _fmt = "%i"; }
      _sd = "[" + num2str(d[1],_fmt) + "]";
    }
  }
  else
  {
    _sd = " (n/a)";
  }

  // did user provide prompt text?
  if (!exist(p))
  { p = "Your input"; }
  if (class(p) != "string")
  { error(_this_function + "string scalar expected!\n"); }

  // wait indefinitely for valid keyboard input
  while (1)
  {
    // do we clear screen before asking question
    if ((_nc))
    {
      clrscr();
      _nc = _nc - _cs;
    }
    // what prompt do we print when asking user for
    // an input
    printf (p);
    if (!strlen(_pc) && !isempty(c))
    {
      if (_dc == "string")
      {
        if (c.n < 8)
        {
          printf(" (%s)", join(c[:]',","));
        }
        else
        {
          printf(" (%s .. %s)", c[1], last(c));
        }
      }
      else if (c.n < 8)
      {
        printf(" (%s)", num2str(c,"%g",","));
      }
      else
      {
        printf(" (%g .. %g)", min(c), max(c));
      }
    }
    printf (_sd);
    printf (" %s ", _qm);

    // read input
    if (_to)
    {
      x = inputs(_to, _ec, _lm);
    }
    else
    {
      x = reads();
    }

    // convert input
    if (x=="" || x=="\n")
    {
      // user pressed enter
      if (_dc == "string")
      { x = blank(0,0); }
      else
      { x = []; }
    }
    else
    {
      // there is some input
      if (_dc == "string")
      {
        x = strsplt(x,",");
      }
      else
      {
        x = "[" + x + "]";
        x = eval(x);
        //x = strtod(x,<<csp=",";lstrip="'BLANK">>);
      }
    }

    if (isempty(x) && isempty(d))
    { continue; }

    // empty 'c' means that everything but empty string goes
    if (isempty(x) && !isempty(d))
    {
      xval = d[1];
      break;
    }

    if (isempty(c))
    {
      if (class(x)=="num")
      {
        if (_dc == "int")
        {
          xval = int(x);
          break;
        }
        else if (_dc == "real")
        {
          xval = x;
          break;
        }
        continue;
      }
    }

    // 'c' is provided, so compare xval against 'c'
    if (class(c)!=class(x))
    { continue; }

    _f = 1;
    for (_j in range(x))
    {
      if (isempty(find(c == x[_j])))
      {
        _f = 0;
        break;
      }
    }
    if (_f)
    {
      xval = x;
      break;
    }
  }

  return xval;
};



//
// fgets: read a line from input file as a text
//
fgets = function(fn)
{
  x = fread(fn,,"char",,"\n");
  if(!isempty(x))
  {
    rval = char(x);
  }
  else
  {
    rval = blank(0,0);
  }
  return rval;
};

beep = function()
{
  printf("%s", char(7));
};


INIT = 1;



