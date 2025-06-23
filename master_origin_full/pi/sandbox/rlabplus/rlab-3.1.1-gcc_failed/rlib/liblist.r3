//
//
//

listrm = function (x, nodes)
{
  global(rval);

  if (!exist(nodes))
  { return x; }

  if (class(nodes)!="string")
  { return x; }

  if (!exist(x))
  { return x; }

  if (class(x)!="list")
  { return x; }

  rval = x;
  for (n in nodes)
  {
    if (!strindex(n,"[") && !strindex(n,"]"))
    { n =  "[\"" + n + "\"]"; }
    eval("clear(rval."+n+")");
  }

  return rval;
};

listmv = function (x, from_node, to_node)
{
  global(rval,rval_from);

  if (!exist(from_node))
  { return x; }

  if (class(from_node)!="string")
  { return x; }

  if (!exist(to_node))
  { return x; }

  if (class(to_node)!="string")
  { return x; }

  if (!exist(x))
  { return x; }

  if (class(x)!="list")
  { return x; }

  rval = x;

  // does from node exist?
  f = strsplt(from_node, ".");
  if (length(f) > 1)
  {
    fn = "";
    for (i in 1:length(f))
    {
      if (strindex(f[i],"[") && strindex(f[i],"]"))
      {
        fn = fn + "." + f[i];
      } else {
        fn = fn + ".[\"" + f[i] + "\"]";
      }
      if (eval("exist(rval"+fn+")"))
      {
        continue;
      } else {
        // from_node does not exist: exit gracefully, nothing to do
        return rval;
      }
    }
    from_node = fn;
  } else {
    from_node = "." + f;
  }
  eval("rval_from = rval"+from_node);

  if (exist(rval_from))
  {
    // remove source node content of which is being moved to new
    // location in the tree
    eval("clear(rval"+from_node+")");

    // check that to_node exists, and if it does not create all necessary
    // supernodes
    t = strsplt(to_node, ".");
    if (!isempty(t))
    {
      if (length(t) > 1)
      {
        tn = "";
        for (i in 1:(length(t)-1))
        {
          if (strindex(t[i],"[") && strindex(t[i],"]"))
          {
            tn = tn + "." + t[i];
          } else {
            tn = tn + ".[\"" + t[i] + "\"]";
          }
          if (!eval("exist(rval"+tn+")"))
          { eval("rval" + tn + "=<<>>"); }
        }
      }
    }
    if (!eval("exist(rval."+to_node+")"))
    {
      eval("rval."+to_node+"=rval_from");
    } else {
      m = eval("members(rval_from)");
      for (mi in m)
      {
        if (!strindex(mi,"[") && !strindex(mi,"]"))
        { mi =  "[\"" + mi + "\"]"; }

        eval("rval." +to_node+"."+mi + "=rval_from."+mi);
      }
      eval("clear(rval_from)");
    }
  }

  return rval;
};

listcp = function (x, from_node, to_node)
{
  global(rval,rval_from);

  if (!exist(from_node))
  { return x; }

  if (class(from_node)!="string")
  { return x; }

  if (!exist(to_node))
  { return x; }

  if (class(to_node)!="string")
  { return x; }

  if (!exist(x))
  { return x; }

  if (class(x)!="list")
  { return x; }

  rval = x;

  // does from node exist?
  f = strsplt(from_node, ".");
  if (length(f) > 1)
  {
    fn = "";
    for (i in 1:length(f))
    {
      if (strindex(f[i],"[") && strindex(f[i],"]"))
      {
        fn = fn + "." + f[i];
      } else {
        fn = fn + ".[\"" + f[i] + "\"]";
      }
      if (eval("exist(rval"+fn+")"))
      {
        continue;
      } else {
        // from_node does not exist: exit gracefully, nothing to do
        return rval;
      }
    }
    from_node = fn;
  } else {
    from_node = "." + f;
  }
  eval("rval_from = rval"+from_node);

  if (exist(rval_from))
  {
    // copy:
    //  do not remove source node content of which is being moved to new
    //  location in the tree
//     eval("clear(rval"+from_node+")");

    // check that to_node exists, and if it does not create all necessary
    // supernodes
    t = strsplt(to_node, ".");
    if (!isempty(t))
    {
      if (length(t) > 1)
      {
        tn = "";
        for (i in 1:(length(t)-1))
        {
          if (strindex(t[i],"[") && strindex(t[i],"]"))
          {
            tn = tn + "." + t[i];
          } else {
            tn = tn + ".[\"" + t[i] + "\"]";
          }
          if (!eval("exist(rval"+tn+")"))
          { eval("rval" + tn + "=<<>>"); }
        }
      }
    }
    if (!eval("exist(rval."+to_node+")"))
    {
      eval("rval."+to_node+"=rval_from");
    } else {
      m = eval("members(rval_from)");
      for (mi in m)
      {
        if (!strindex(mi,"[") && !strindex(mi,"]"))
        { mi =  "[\"" + mi + "\"]"; }

        eval("rval." +to_node+"."+mi + "=rval_from."+mi);
      }
      eval("clear(rval_from)");
    }
  }

  return rval;
};

