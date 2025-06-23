//
// libgeneral.r:  unspecific functions
//
static(_INIT_);

if (exist(_INIT_))
{ EOF }


progress = function(idx, n, opts)
{
  static(prevlen);
  if (!exist(prevlen))
  {
    prevlen = 0;
  }

  this_solver = "counter";

  if (class(idx)!="num")
  {
    error(this_solver + ": first argument has to be positive integer scalar!\n");
  }
  if (class(n)!="num")
  {
    error(this_solver + ": second argument has to be positive integer scalar!\n");
  }

  x_0 = "";
  fmt = "%.0f of %.0f";
  completed_msg = "";
  fin = 0;
  if (class(opts)=="list")
  {
    if (class(opts.fmt)=="string")
    {
      fmt = opts.fmt;
    }
    if (class(opts.preamble)=="string")
    {
      x_0 = opts.preamble;
    }
    if (class(opts.completed)=="string")
    {
      completed_msg = opts.completed;
    }
    if (class(opts.terminate)=="num")
    {
      fin = (opts.terminate>0);
    }
  }

  for (i in 1:prevlen)
  {
    printf("\b");
  }
  prevlen = sprintf(x,fmt, idx, n);

  if (idx == 1)
  {
    printf(x_0);
  }

  printf(x);

  if ( ((idx == n) || (fin)) && (strlen(completed_msg)>0) )
  {
    for (i in 1:prevlen)
    {
      printf("\b");
    }
    clrpos();
    printf(completed_msg);
  }

  return 0;
};

gcd = function(i1, i2)
{
  if(!exist(i1) || !exist(i2))
  { return []; }

  if (any(i1 !=int(i1)))
  { return []; }

  if (any(i2 !=int(i2)))
  { return []; }

  rmax = max(i1.nr,i2.nr);
  cmax = max(i1.nc,i2.nc);
  rval = zeros(rmax,cmax);
  for (r in 1:rmax)
  {
    for (c in 1:cmax)
    {
      r1 = min(r,i1.nr);
      c1 = min(c,i1.nc);
      r2 = min(r,i2.nr);
      c2 = min(c,i2.nc);

      a = i1[r1;c1];
      b = i2[r2;c2];
      while (b != 0)
      {
        t = b;
        b = mod(a, t);
        a = t;
      }
      rval[r;c] = a;
    }
  }

  return rval;
};

if (!exist( create_charactertable ))
{
  //
  // a column vector of characters of the latin alphabet
  //
  create_charactertable = function()
  {
    return char([65:90]');
  }; // end of create_charactertable = function()
}

//
// access the data matrix from the end
//
last = function (x, n)
{
  if (!exist(x))
  {
    printf("last: Retrieve rows from the end of a data matrix 'x'.\n");
    printf("last: Format:\n");
    printf("last:   y = last(x, n)\n");
    printf("last: which is equivalent to   y = x[x.nr-1+n;]  .\n");
    return [];
  }
  if (isempty(x))
  { return []; }

  if (!exist(n))
  { n = 1; }
  if (x.nc != 1 && x.nr != 1)
  {
      // full matrix
    if (max(n) > x.nr)
    { error ("last: index out of bounds!"); }
    if (min(n) < 1)
    { error ("last: index out of bounds!"); }
    return x[ x.nr+1-n ;];
  }

  // vector
  if (max(n) > length(x))
  { error ("last: index out of bounds!"); }
  if (min(n) < 1)
  { error ("last: index out of bounds!"); }
  return x[ length(x)+1-n ];
}; // end of last = function (x, n)

//
// access the last row of a data matrix
//
lastr = function (x, n)
{
  if (!exist(x))
  {
    printf("lastr: Retrieve rows from the end of a data matrix 'x'.\n");
    printf("lastr: Format:\n");
    printf("lastr:   y = lastr(x, n)\n");
    return [];
  }
  if (isempty(x))
  { return []; }

  if (x.storage != "dense")
  {
    printf ("lastr: Improper argument.");
    error ();
  }
  if (!exist(n))
  { n = 1; }
  if (!isreal(n))
  { n = real(n); }
  if (n.n != 1)
  { n = n[1]; }
  return x[max(1,x.nr-int(n)+1):x.nr;];
}; // end of last = function (x, n)

//
// access the last column of a data matrix
//
lastc = function (x, n)
{
  if (!exist(x))
  {
    printf("lastc: Retrieve columns from the end of a data matrix 'x'.\n");
    printf("lastc: Format:\n");
    printf("lastc:   y = lastc(x, n)\n");
    return [];
  }
  if (isempty(x))
  { return []; }

  if (x.storage != "dense")
  {
    printf ("lastr: Improper argument.");
    error ();
  }
  if (!exist(n))
  { n = 1; }
  if (!exist(n))
  { n = 1; }
  if (!isreal(n))
  { n = real(n); }
  if (n.n != 1)
  { n = n[1]; }
  return x[;max(1,x.nc-n+1):x.nc];
}; // end of last = function (x, n)

//
// create range of inidices necessary to access a vector or a matrix
//
range = function(x,s)
{
  //
  // s = <<dir;dir2;range;range2>>
  //
  _d  = 1;
  _r  = [];
  if (exist(s))
  {
    if (exist(s.dir))
    {
      if (s.dir == -1)
      { _d = -1; }
    }
    if (exist(s.range))
    {
      _r = s.range;
    }
  }

  if (isscalar(x))
  { return 1; }

  if (isvector(x))
  {
    n = length(x);

    if (n==0)
    { return []; }

    if (_d == 1)
    {
      rval = [1:n];
    }

    if (_d == -1)
    {
      rval = [n:1:-1];
    }

    if (isempty(_r))
    {
      return rval;
    }

    return rval[_r];
  }

  _d2 = 1;
  _r2 = [];
  if (exist(s))
  {
    if (exist(s.dir2))
    {
      if (s.dir2 == -1)
      { _d2 = -1; }
    }
    if (exist(s.range2))
    {
      _r2 = s.range2;
    }
  }

  if (ismatrix(x))
  {
    nr = x.nr;
    nc = x.nc;

    if ((nr==0)||(nc==0))
    { return []; }

    if (_d == 1)
    {
      rval1 = [1:nr];
    }
    if (_d == -1)
    {
      rval1 = [nr:1:-1];
    }

    if (!isempty(_r))
    {
      rval1 = rval1[_r];
    }

    rval1?

    if (_d2 == 1)
    {
      rval2 = [1:nc];
    }
    if (_d2 == -1)
    {
      rval2 = [nc:1:-1];
    }

    if (!isempty(_r2))
    {
      rval2 = rval2[_r2];
    }

    rval = zeros(rval1.n * rval2.n, 2);
    k = 0;
    for (i in 1:rval1.n)
    {
      for (j in 1:rval2.n)
      {
        k++;
        rval[k;] = [rval1[i], rval2[j]];
      }
    }

    return rval;
  }

  return [];
};

names = function( s )
{
  if (!exist(s))
  { return blank(0,0); }

  if (class(s)!="string")
  { return blank(0,0); }

  rval = grep(members($$)', s[1]);

  return rval;
};

join = function (arg1, js, fmt)
{
  rval = "[]";
  if (isempty(arg1))
  { return rval; }

  if(class(arg1)=="string")
  {
    if (!exist(js))
    { js = ","; }
    rval = sum(arg1, js);
  }

  if (class(arg1)=="num")
  {
    if (!exist(f))
    { fmt = "%g"; }
    if (!exist(js))
    { js  = ",";}
    rval = num2str(arg1, fmt, js);
  }

  return rval;
};

isrange = function(x)
{
  if (!isvector(x))
  { return 0; }

  if (x.n <= 2)
  { return 1; }

  d2x = diff(diff(x));
  if (all(d2x==0))
  { return 1; }

  return 0;
};

nrange = function(x)
{
  if (!isrange(x))
  { return []; }

  if (x.n <= 2)
  { return x; }

  x1 = min(x);
  x2 = max(x);
  dx = diff(x)[1];

  if (dx == 1)
  { return [x1,x2]; }

  return [x1,x2,dx];
};




_INIT_ = 1;


