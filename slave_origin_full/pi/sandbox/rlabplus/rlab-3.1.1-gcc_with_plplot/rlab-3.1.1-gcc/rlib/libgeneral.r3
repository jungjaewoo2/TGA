//
// libgeneral.r:  unspecific functions
//
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
  } else {
      // vector
    if (max(n) > length(x))
    { error ("last: index out of bounds!"); }
    if (min(n) < 1)
    { error ("last: index out of bounds!"); }
    return x[ length(x)+1-n ];
  }
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
range = function(x,bof,eof,bof2,eof2)
{
  if (!exist(x))
  {
    error("range: Cannot find range of undefined argument\n");
  }
  if (!exist(bof))
  { bof=0; }
  if (!exist(eof))
  { eof=0; }

  if (x.nr==1 || x.nc==1)
  {
    rval = [(1+bof):(x.nr*x.nc+eof)];
    return rval;
  }

  if (!exist(bof2))
  { bof2=0; }
  if (!exist(eof2))
  { eof2=0; }

  rval  = [];
  rvaly = [(1+bof):(x.nc+eof)]';
  for (i in (1+bof2):(x.nr+eof2))
  {
    rvalx = i .* ones(rvaly);
    rval = [rval; [rvalx, rvaly]];
  }
  return rval;
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
  if (!isvec(x))
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


