//
// recursive printout of the members of a list
//

listinfo = function(x)
{
  if (class(x) == "list")
  {
    rval = blank(0,0);
    for (i in members(x))
    {
      if (!isempty(strtod(i)))
      {
        is = "[" + i + "]";
      else
        is = i;
      }
      if (class(x.[i])=="list")
      {
        rval = [rval; is + "." + listinfo(x.[i])];

      else
        rval = [rval; is];
      }
    }
    return rval;
  }
};

x = <<>>;
x.a = 1;
x.b = rand(3,3);
x.c = <<>>;
x.c.[1] = 1;
x.c.[2] = 2;
x.d = <<>>;
x.d.[1] = 1;
x.d.[2] = <<>>;
x.d.[2].a = 1;
x.d.[2].b = <<>>;
x.d.[2].b.[1] = rand(3,3);
x.d.[2].b.[2] = 1;

printf("Outline of the list variable 'x' is:\n");
"x." + listinfo(x)



