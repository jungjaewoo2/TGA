//
// functions that help in statistical analysis
//
isstat = function (x)
{
  if (!exist(x))
  { return 0; }

  if (class(x)!="list")
  { return 0;}

  if (!exist(x.val))
  { return 0; }

  if (!exist(x.wgt))
  { return 0; }

  if (class(x.wgt)!="num")
  { return 0; }

  if (type(x.wgt)!="real")
  { return 0; }

  // weights cannot be empty
  if (isempty(x.wgt))
  { return 0; }

  if (isempty(x.val))
  { return 0; }

  // weight has to be positive
  if (any(isnan(x.wgt)))
  { return 0; }
  if (any(isnan(x.val)))
  { return 0; }
  if (any(x.wgt<=0))
  { return 0; }

  s_val = size(x.val);
  s_wgt = size(x.wgt);

  // cannot have weights for non-existing values
  if (any(s_wgt > s_val))
  { return 0; }

  // return dimension of val:
  if (any(s_val==1))
  { return 1;}

  // return no. of columns in x.val
  return s_val[2];
};

rattle = function (x, e, d, d_avg)
{
  global(uniform);

  if (class(x)!="num")
  { return x; }

  if (class(d)!="function")
  { d = uniform; }

  if (!exist(d_avg))
  { d_avg = 0.5; }

  if (!exist(e))
  { e = 0.1; }

  rval = (1 + e *(d(x) - d_avg)) .* x;
  return rval;
};

countunique = function (x)
{
  global(range);

  if (!exist(x))
  { return []; }

  s_x = unique(x);

  rval = <<>>;
  rval.unique = s_x;
  rval.count  = zeros(s_x);

  for (i in range(s_x))
  {
    rval.count[i] = sum(x == s_x[i]);
  }

  
  return rval;
};



